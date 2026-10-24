! Module for 1d evp dynamics
! Mimics the 2d B grid solver
! functions in this module includes conversion from 1d to 2d and vice versa.
! cpp flag _OPENMP_TARGET is for gpu. Otherwize optimized for cpu
! FIXME: For now it allocates all water point, which in most cases could be avoided.
!===============================================================================
! Created by Till Rasmussen (DMI), Mads Hvid Ribergaard (DMI), and Jacob W. Poulsen, Intel

module ice_dyn_evp1d

  !- modules -------------------------------------------------------------------
  use ice_kinds_mod
  use ice_blocks, only: nx_block, ny_block, nghost
  use ice_constants
  use ice_communicate, only: my_task, master_task
  use ice_domain_size, only: max_blocks, nx_global, ny_global
  use ice_fileunits, only: nu_diag
  use ice_exit, only: abort_ice
  use ice_dyn_shared, only: rheo_area_min, free_drift_approx, evp1d_reorder

  !- directives ----------------------------------------------------------------
  implicit none
  private

  !- public routines -----------------------------------------------------------
  public :: dyn_evp1d_init, dyn_evp1d_run, dyn_evp1d_finalize

  !- private routines ----------------------------------------------------------

  !- private vars --------------------------------------------------------------
  ! nx and ny are module variables for arrays after gather (G_*) Dimension according to CICE is
  ! nx_global+2*nghost, ny_global+2*nghost
  ! nactive are number of active points (both t and u). navel is number of active
  ! nActive and navel are computed on the master alone, yet are read on every
  ! task as the trip count handed to gather_slot_r8.  Zero-initialised so a
  ! non-master task passes a defined value rather than whatever was on the
  ! stack; the helpers return before using it.
  integer(kind=int_kind), save :: nx, ny, nallocated
  integer(kind=int_kind), save :: nActive = 0, navel = 0
  ! The slots the kernels actually visit, aList(1:nIce), ascending.
  !
  ! With evp1d_reorder off, nIce == nActive and aList is the identity, so the
  ! loops are the contiguous ones they always were.  With it on, aList holds
  ! only the slots where iceTmask .or. iceUmask, so the 240 subcycles never
  ! stream the ice-free ocean -- and because ice is geographically clustered
  ! and slots run in cell order, the entries come in long contiguous runs
  ! rather than scattered singletons.
  !
  ! ⚠ It costs no extra pass.  set_skipMe already walks 1..nActive every
  ! timestep testing exactly this union -- it was counting the result into a
  ! local and throwing it away.  So whatever fraction of cells changes state
  ! from one step to the next, the list is rebuilt inside work that was
  ! already being done.
  integer(kind=int_kind), save :: nIce = 0
  logical(kind=log_kind), save :: first_reorder = .true.
  integer(kind=int_kind), allocatable, dimension(:) :: aList

  ! indexes
  integer(kind=int_kind), allocatable, dimension(:,:) :: iwidx
  logical(kind=log_kind), allocatable, dimension(:)   :: skipTcell,skipUcell
  integer(kind=int_kind), allocatable, dimension(:)   :: nbrW,nbrSW,nbrS,nbrE,nbrNE,nbrN ! arrays for neighbour points
  integer(kind=int_kind), allocatable, dimension(:)   :: indxti, indxtj, indxTij
  ! Scratch map over the gathered grid, linear index i+(j-1)*nx -> slot.
  ! -1 marks a cell the solver needs before slots are handed out, 0 a cell it
  ! does not. Built by calc_navel, consumed and released by convert_2d_1d_init,
  ! so it exists only during init and only on the master task.
  integer(kind=int_kind), allocatable, dimension(:)   :: ijslot

  ! ONE global block scratch, shared by every gather that lands straight in a
  ! 1d array -- see gather_slot_r8.  It replaces the 23 (nx,ny) automatics
  ! dyn_evp1d_run used to put on the stack every timestep.  Allocated on every
  ! task, not just the master: it is the actual argument of a collective, so it
  ! has to be a valid array everywhere even though only the master reads it.
  real   (kind=dbl_kind), allocatable, dimension(:,:) :: G_work

  ! 1D arrays to allocate

  ! Grid
  real   (kind=dbl_kind), allocatable, dimension(:)   :: &
     HTE_1d,HTN_1d, HTEm1_1d,HTNm1_1d, dxT_1d, dyT_1d, uarear_1d, &
     DminTarea_1d, dxhy_1d, dyhx_1d, fcor_1d

  ! time varying
  real(kind=dbl_kind)   , allocatable, dimension(:)   ::                  &
    cdn_ocn,aiu,uocn,vocn,waterxU,wateryU,forcexU,forceyU,umassdti,fmU,   &
    rheofactU,                                                            &
    umass,                                                                &
    strintxU_1d,strintyU_1d,uvel_init,vvel_init, strength, uvel, vvel,          &
    stressp_1, stressp_2, stressp_3, stressp_4, stressm_1, stressm_2,     &
    stressm_3, stressm_4, stress12_1, stress12_2, stress12_3, stress12_4, &
    str1, str2, str3, str4, str5, str6, str7, str8, Tbu, Cb

  ! Every boundary condition reduces to dst = s1 + w*(s1 - s2):
  !   cyclic         s1 is the matching cell on the opposite edge, w = 0
  !   zero_gradient  s1 is the cell just inside the edge,          w = 0
  !   linear_extrap  s1, s2 the two cells just inside,   w = the ghost depth
  ! One ordered list for all of them, built once in halo_sweep.
  integer(kind=int_kind), allocatable, dimension(:)   ::                       &
    halo_bc_dst, halo_bc_s1, halo_bc_s2
  real   (kind=dbl_kind), allocatable, dimension(:)   :: halo_bc_w
  integer(kind=int_kind)                              :: n_halo_bc

  ! The 1d stresses are subcycle state and persist between timesteps, so
  ! gather_dyn does not fetch them.  They must still be seeded once from the
  ! block arrays, because dyn_evp1d_init runs BEFORE the restart read: on a
  ! restarted run the block copies carry the file values while the 1d arrays
  ! are still zero.  One gather on the first call covers that and a cold start.
  logical(kind=log_kind), save :: stresses_seeded = .false.

!=============================================================================
  contains
!=============================================================================
! module public subroutines
! In addition all water points are assumed to be active and allocated thereafter.
!=============================================================================

  subroutine dyn_evp1d_init

    use ice_grid, only: G_HTE, G_HTN

    implicit none

    ! local variables

    real(kind=dbl_kind)   , allocatable, dimension(:,:) :: G_dyT, G_dxT, G_uarear, &
                                                           G_DminTarea, G_fcor
    logical(kind=log_kind), allocatable, dimension(:,:) :: G_tmask

    integer(kind=int_kind) :: ierr

    character(len=*), parameter :: subname = '(dyn_evp1d_init)'

    nx=nx_global+2*nghost
    ny=ny_global+2*nghost

    allocate(G_dyT(nx,ny),G_dxT(nx,ny),G_uarear(nx,ny),G_tmask(nx,ny), &
             G_DminTarea(nx,ny), G_fcor(nx,ny),G_work(nx,ny),stat=ierr)
    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

    ! gather from blks to global
    call gather_static(G_uarear, G_dxT, G_dyT, G_DminTarea, G_fcor, G_tmask)

    ! calculate number of water points (T and U). Only needed for the static version
    ! tmask in ocean/ice
    if (my_task == master_task) then
      call calc_nActiveTU(G_tmask,nActive)
      call evp1d_alloc_static_na(nActive)
      call calc_2d_indices_init(G_tmask)
      call calc_navel(nActive, navel)
      call evp1d_alloc_static_navel(navel)
      call numainit(1,nActive,navel+1)
      call convert_2d_1d_init(nActive,G_HTE, G_HTN, G_uarear, G_dxT, G_dyT, &
                              G_DminTarea, G_fcor)
      call evp1d_alloc_static_halo()
      ! The halo geometry -- indxTij, na0, navel and the boundary types -- does
      ! not change, so the list is built once here.
      call halo_sweep(.false.)
      ! ijslot has served its purpose: slots are assigned and both halo lists
      ! are built.
      deallocate(ijslot, stat=ierr)
      if (ierr/=0) then
         call abort_ice(subname//' ERROR: deallocating', file=__FILE__, line=__LINE__)
      endif

    endif

    deallocate(G_dyT,G_dxT,G_uarear,G_tmask,G_DminTarea,G_fcor,stat=ierr)
    if (ierr/=0) then
       call abort_ice(subname//' ERROR: deallocating', file=__FILE__, line=__LINE__)
    endif

  end subroutine dyn_evp1d_init

!=============================================================================

  subroutine dyn_evp1d_run(dt          ,                                           &
                           L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
                           L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
                           L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4, &
                           L_umass     , L_uvel      , &
                           L_vvel      , L_icetmask  , L_iceUmask)

    use ice_calendar, only : write_history, write_restart, write_ic
    use ice_dyn_shared, only : ndte, strintxU_slot, strintyU_slot, &
           stressp_1_slot, stressp_2_slot, stressp_3_slot, stressp_4_slot, &
           stressm_1_slot, stressm_2_slot, stressm_3_slot, stressm_4_slot, &
           stress12_1_slot, stress12_2_slot, stress12_3_slot, stress12_4_slot
    use ice_grid      , only : dynMap, dynOff, dynCnt, unpack_from_slots, &
           pack_to_slots
    use ice_domain    , only : nblocks
    ! renamed: this module's own strintxU is the 1-D GATHERED array
    use ice_flux      , only : strintxU, strintyU
    use ice_dyn_core1d, only : stress_1d, stepu_1d, calc_diag_1d
    ! strintxU/strintyU and taubxU/taubyU are diagnostics: history is their
    ! only consumer.  When no stream wants them this step, calc_diag_1d is
    ! skipped and four of scatter_dyn's six collectives go with it.
    use ice_history_shared, only: hist_field_needed, &
        n_strintx, n_strinty, n_siforceintstrx, n_siforceintstry, &
        n_taubx  , n_tauby  , n_siforceintstrx, n_siforceintstry
    use ice_timers    , only : ice_timer_start, ice_timer_stop, timer_evp1dcore

    use icepack_intfc , only : icepack_query_parameters, icepack_warnings_flush, &
      icepack_warnings_aborted

    implicit none

    ! nx_block, ny_block, max_blocks
    ! uvel/vvel have their i and j dimensions flattened; gather_global and
    ! scatter_global resolve to their _flat specifics
    real(kind=dbl_kind)   , dimension(:,:)  , intent(inout) :: &
      L_uvel      , L_vvel

    ! these have their i and j dimensions flattened, i.e. element (i,j,iblk)
    ! lives at ((j-1)*nx_block + i, iblk)
    real(kind=dbl_kind)   , dimension(:,:)  , intent(inout) ::  &
      L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 ,  &
      L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 ,  &
      L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4
    real(kind=dbl_kind), intent(in) :: dt   ! dynamics timestep

    ! flattened as above.  umass is the last block field this solver still
    ! takes: strength, cdn_ocnU, aiU, uocnU, vocnU, forcexU and forceyU used
    ! to come in too, and were never read -- dyn_prep2 fills the slots as
    ! well as the blocks, and the gather has taken the slot copies since
    ! 2677687.  Only the declaration outlived the use.
    real(kind=dbl_kind)   , dimension(:,:)  , intent(in) ::    &
      L_umass
    ! L_iceTmask and L_iceUmask have their i and j dimensions flattened,
    ! i.e. element (i,j,iblk) lives at ((j-1)*nx_block + i, iblk)
    logical(kind=log_kind), dimension(:,:), intent(in) ::      &
      L_iceTmask  , L_iceUmask

    ! local variables

    ! nx, ny
    ! Only the globals that still have a reason to exist.  uvel and vvel
    ! round-trip: gathered in, written at this solver's cells, scattered back,
    ! so every cell evp1d does not own has to keep the value that arrived.
    ! The four diagnostics do the same, on the steps history asks for them.
    ! Everything else now goes straight to the 1d arrays through G_work --
    ! 23 (nx,ny) automatics off the stack, which on a 1440x1080 grid is about
    ! 345 MB per call.
    real(kind=dbl_kind),    dimension(nx,ny) ::                &
      G_strintxU  , G_strintyU  , G_taubxU    , G_taubyU    ,  &
      G_uvel      , G_vvel
    ! G_taubxU and G_taubyU are post processed from Cb
    logical(kind=log_kind), dimension (nx,ny) ::                 &
      G_iceUmask  , G_iceTmask

    ! does any history stream want strint/taub this step?
    logical(kind=log_kind) :: lstrint, ltaub, lsiforce
    integer(kind=int_kind) :: iblk   ! block index for the siforce expansion

    character(len=*), parameter :: subname = '(dyn_evp1d_run)'

    integer(kind=int_kind) :: ksub, iw
    logical(kind=log_kind) :: do_seed

    real   (kind=dbl_kind) :: rhow

    ! From 3d to 2d on master task
    ! ⚠ These MUST be set before gather_dyn, not after it.  They decide
    ! whether four of its collectives are entered at all, and they were being
    ! read here while still undefined -- an uninitialised stack logical
    ! steering a collective.  On one task that merely risked gathering, or
    ! not gathering, the four diagnostics at random; on more than one it is a
    ! hang, because the ranks need not agree on an undefined value.  Nothing
    ! in hist_field_needed depends on the gather, so the answer is simply to
    ! ask first.
    lstrint = hist_field_needed(n_strintx) .or. hist_field_needed(n_strinty) .or. &
              hist_field_needed(n_siforceintstrx) .or. &
              hist_field_needed(n_siforceintstry)
    ltaub   = hist_field_needed(n_taubx  ) .or. hist_field_needed(n_tauby  )
    lsiforce = hist_field_needed(n_siforceintstrx) .or. &
               hist_field_needed(n_siforceintstry)

    call gather_dyn(&
                    L_umass     , L_uvel      , L_vvel      ,               &
                    L_iceTmask  , L_iceUmask  ,                             &
                    G_strintxU  , G_strintyU  , G_taubxU    , G_taubyU    , &
                    lstrint     , ltaub       ,                             &
                    G_uvel      , G_vvel      ,                             &
                    G_iceTmask,  G_iceUmask)

    ! One-time seed; see stresses_seeded above.
    !
    ! do_seed is taken on EVERY task and the module flag is cleared on every
    ! task too, immediately after the collective.  Setting the flag inside the
    ! master-only block below leaves it .false. on every other rank, so they
    ! keep calling a collective that the master has stopped entering -- which
    ! is a hang on any run with more than one task, and invisible on one.
    do_seed = .not. stresses_seeded
    if (do_seed) then
       ! Writes the 1d stresses itself now, so the twelve globals and the
       ! loop that copied out of them are both gone.  Still called on every
       ! task -- it is twelve collectives -- and still master-only in effect,
       ! because gather_slot_r8 returns before the copy on everyone else.
       call gather_stresses()
    endif
    stresses_seeded = .true.

    if (my_task == master_task) then
       call set_skipMe(G_iceTmask, G_iceUmask,nActive)

       ! Map from 2d to 1d
       call convert_2d_1d_dyn(nActive, dt, G_uvel, G_vvel)

       ! map from cpu to gpu (to) and back.
       ! This could be optimized considering which variables change from time step to time step
       ! and which are constant.
       ! in addition initialization of Cb and str1, str2, str3, str4, str5, str6, str7, str8
       call icepack_query_parameters(rhow_out=rhow)
       call icepack_warnings_flush(nu_diag)
       if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

       call ice_timer_start(timer_evp1dcore)
#ifdef _OPENMP_TARGET
       !$omp target data map(to: nbrW, nbrSW, nbrS, nbrE, nbrNE, nbrN, skipUcell, skipTcell,&
       !$omp                 strength, dxT_1d, dyT_1d, DminTarea_1d,           &
       !$omp                 dxhy_1d, dyhx_1d,                                 &
       !$omp                 HTE_1d,HTN_1d,HTEm1_1d,                           &
       !$omp                 HTNm1_1d,forcexU, forceyU, umassdti, fmU,         &
       !$omp                 uarear_1d,uvel_init, vvel_init, Tbu, Cb,          &
       !$omp                 str1, str2, str3, str4, str5, str6, str7, str8,   &
       !$omp                 cdn_ocn, aiu, uocn, vocn, waterxU, wateryU, rhow, &
       !$omp                 rheofactU                                            &
       !$omp             map(tofrom: uvel,vvel,                                &
       !$omp                 stressp_1, stressp_2, stressp_3, stressp_4,       &
       !$omp                 stressm_1, stressm_2, stressm_3, stressm_4,       &
       !$omp                 stress12_1,stress12_2,stress12_3,stress12_4)
       !$omp target update to(arlx1i,denom1,capping,deltaminEVP,e_factor,epp2i,brlx)
#endif
       ! initialization of str? in order to avoid influence from old time steps
       str1(1:navel+1)=c0
       str2(1:navel+1)=c0
       str3(1:navel+1)=c0
       str4(1:navel+1)=c0
       str5(1:navel+1)=c0
       str6(1:navel+1)=c0
       str7(1:navel+1)=c0
       str8(1:navel+1)=c0

       do ksub = 1,ndte        ! subcycling
          call stress_1d (nbrW, nbrSW, nbrS, 1, nIce, aList,                                    &
                          uvel, vvel, dxT_1d, dyT_1d, DminTarea_1d,                  &
                          dxhy_1d, dyhx_1d, skipTcell, strength,                     &
                          HTE_1d, HTN_1d, HTEm1_1d, HTNm1_1d,                        &
                          stressp_1,  stressp_2,  stressp_3,  stressp_4,             &
                          stressm_1,  stressm_2,  stressm_3,  stressm_4,             &
                          stress12_1, stress12_2, stress12_3, stress12_4,            &
                          str1, str2, str3, str4, str5, str6, str7, str8)

          call stepu_1d  (1, nIce, aList, cdn_ocn, aiu, uocn, vocn,                         &
                          waterxU, wateryU, forcexU, forceyU, umassdti, fmU, uarear_1d, &
                          uvel_init, vvel_init, uvel, vvel,                             &
                          str1, str2, str3, str4, str5, str6, str7, str8,               &
                          nbrE, nbrNE, nbrN, skipUcell, Tbu, Cb, rhow, rheofactU)
          call evp1d_halo_update()
       enddo
       ! Skipped outright when no history stream wants strintx/strinty --
       ! which is what the comment here used to say was possible.
       if (lstrint) &
       call calc_diag_1d(1        , nIce     , aList, &
                         uarear_1d, skipUcell, &
                         str1     , str2     , &
                         str3     , str4     , &
                         str5     , str6     , &
                         str7     , str8     , &
                         nbrE       , nbrNE       , &
                         nbrN      , rheofactU, &
                         strintxU_1d, strintyU_1d)

       call ice_timer_stop(timer_evp1dcore)

#ifdef _OPENMP_TARGET
       !$omp end target data
#endif
       ! Map results back to 2d
       call convert_1d_2d_dyn(nActive, navel,                                         &
                              G_strintxU  , G_strintyU  , &
                              G_uvel      , G_vvel      , G_taubxU    , &
                              G_taubyU)

    endif ! master_task

    call scatter_dyn(&
                     L_uvel      , L_vvel      ,                             &
                     G_strintxU  , G_strintyU  , G_uvel      , G_vvel      , &
                     G_taubxU    , G_taubyU    , lstrint     , ltaub)

    ! ice_history reads strintxU/strintyU and taubxU/taubyU from the slots the
    ! scatter above just filled, so the block forms are not needed -- except
    ! by siforceintstrx/y, the one field that still forms a block-layout
    ! product with aice.  Expand for that case, and only for it.
    if (lsiforce) then
       !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
       do iblk = 1, nblocks
          call unpack_from_slots(strintxU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintxU(:,iblk))
          call unpack_from_slots(strintyU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintyU(:,iblk))
       enddo
       !$OMP END PARALLEL DO
    endif

    ! The twelve stresses go back only when a file is due.  Nothing outside
    ! the dynamics reads them in between -- ice_history touches them only
    ! inside its own write_history gate -- and the 1d arrays carry the state
    ! onward, so scattering every timestep was writing values no one read.
    !
    ! strintxU, strintyU, taubxU and taubyU cannot be deferred to a write
    ! step the way the stresses are -- with hist_avg the accumulator sums
    ! them every timestep -- but they CAN be skipped entirely when no stream
    ! asks for them, which is what lstrint/ltaub above do.  hist_field_needed
    ! also covers the snapshot case, where only the write step's value
    ! survives.
    if (any(write_history) .or. write_restart == 1 .or. write_ic) then
       call scatter_stresses(L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
                             L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
                             L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4)

       ! ⚠ The restart and the history read the SLOT copies now, so the
       ! scatter above is only half the job on this path: the block arrays
       ! it just wrote have to reach the slots too.  Same I/O gate as the
       ! scatter, so this costs nothing on an ordinary timestep.
       !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
       do iblk = 1, nblocks
          call pack_to_slots(L_stressp_1(:,iblk), dynMap(:,iblk), &
                             stressp_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
          call pack_to_slots(L_stressp_2(:,iblk), dynMap(:,iblk), &
                             stressp_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
          call pack_to_slots(L_stressp_3(:,iblk), dynMap(:,iblk), &
                             stressp_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
          call pack_to_slots(L_stressp_4(:,iblk), dynMap(:,iblk), &
                             stressp_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
          call pack_to_slots(L_stressm_1(:,iblk), dynMap(:,iblk), &
                             stressm_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
          call pack_to_slots(L_stressm_2(:,iblk), dynMap(:,iblk), &
                             stressm_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
          call pack_to_slots(L_stressm_3(:,iblk), dynMap(:,iblk), &
                             stressm_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
          call pack_to_slots(L_stressm_4(:,iblk), dynMap(:,iblk), &
                             stressm_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
          call pack_to_slots(L_stress12_1(:,iblk), dynMap(:,iblk), &
                             stress12_1_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
          call pack_to_slots(L_stress12_2(:,iblk), dynMap(:,iblk), &
                             stress12_2_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
          call pack_to_slots(L_stress12_3(:,iblk), dynMap(:,iblk), &
                             stress12_3_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
          call pack_to_slots(L_stress12_4(:,iblk), dynMap(:,iblk), &
                             stress12_4_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
       enddo
       !$OMP END PARALLEL DO
    endif
    ! calculate number of active points. allocate if initial or if array size should increase
    ! call calc_nActiveTU(iceTmask_log,nActive, iceUmask)
    ! if (nActiveold ==0) then ! first
    !     call evp_1d_alloc(nActive, nActive,nx,ny)
    !     nactiveold=nActive+buf1d ! allocate
    !     call init_unionTU(nx, ny, iceTmask_log,iceUmask)
    ! else if (nactiveold < nActive) then
    !     write(nu_diag,*) 'Warning nActive is bigger than old allocation. Need to re allocate'
    !     call evp_1d_dealloc() ! only deallocate if not first time step
    !     call evp_1d_alloc(nActive, nActive,nx,ny)
    !     nactiveold=nActive+buf1d ! allocate
    !     call init_unionTU(nx, ny, iceTmask_log,iceUmask)
    ! endif
    ! call cp_2dto1d(nActive)
    ! FIXME THIS IS THE LOGIC FOR RE ALLOCATION IF NEEDED
    ! call add_1d(nx, ny, natmp, iceTmask_log, iceUmask, ts)

  end subroutine dyn_evp1d_run

!=============================================================================

  subroutine dyn_evp1d_finalize()
    implicit none

    character(len=*), parameter :: subname = '(dyn_evp1d_finalize)'

    if (my_task == master_task) then
       write(nu_diag,*) 'Close evp 1d log'
    endif

    ! G_work outlives dyn_evp1d_init, unlike the other gathered globals there.
    if (allocated(G_work)) deallocate(G_work)

  end subroutine dyn_evp1d_finalize

!=============================================================================

  subroutine evp1d_alloc_static_na(na0)
    implicit none

    integer(kind=int_kind), intent(in) :: na0
    integer(kind=int_kind) :: ierr
    character(len=*), parameter :: subname = '(evp1d_alloc_static_na)'

    allocate(skipTcell(1:na0),  &
             skipUcell(1:na0),  &
             iwidx(1:nx,1:ny),  &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

    allocate(indxTi(1:na0), &
             indxTj(1:na0), &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

    allocate(nbrW(1:na0) , &
             nbrSW(1:na0) , &
             nbrS(1:na0) , &
             nbrE(1:na0) , &
             nbrNE(1:na0) , &
             nbrN(1:na0), &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

    allocate( HTE_1d    (1:na0), &
              HTN_1d    (1:na0), &
              HTEm1_1d  (1:na0), &
              HTNm1_1d  (1:na0), &
              dxT_1d    (1:na0), &
              dyT_1d    (1:na0), &
              DminTarea_1d(1:na0), &
              dxhy_1d   (1:na0), &
              fcor_1d   (1:na0), &
              dyhx_1d   (1:na0), &
              strength  (1:na0), &
              stressp_1 (1:na0), &
              stressp_2 (1:na0), &
              stressp_3 (1:na0), &
              stressp_4 (1:na0), &
              stressm_1 (1:na0), &
              stressm_2 (1:na0), &
              stressm_3 (1:na0), &
              stressm_4 (1:na0), &
              stress12_1(1:na0), &
              stress12_2(1:na0), &
              stress12_3(1:na0), &
              stress12_4(1:na0), &
              stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

    allocate(cdn_ocn  (1:na0), aiu      (1:na0), rheofactU(1:na0), &
             uocn     (1:na0), vocn     (1:na0), &
             waterxU  (1:na0), wateryU  (1:na0), &
             forcexU  (1:na0), forceyU  (1:na0), &
             umassdti (1:na0), fmU      (1:na0), &
             umass    (1:na0),                   &
             uarear_1d(1:na0),                   &
             strintxU_1d (1:na0), strintyU_1d (1:na0), &
             Tbu      (1:na0), Cb       (1:na0), &
             uvel_init(1:na0), vvel_init(1:na0), &
             aList    (1:na0),                   &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

    ! The identity, which is what the kernels walk when evp1d_reorder is off.
    ! Written once here rather than rebuilt every timestep: with the flag off
    ! the list never changes.
    do ierr = 1, na0
       aList(ierr) = ierr
    enddo
    nIce = na0

  end subroutine evp1d_alloc_static_na

!=============================================================================

  subroutine evp1d_alloc_static_navel(navel0)
    implicit none

    integer(kind=int_kind), intent(in) :: navel0
    integer(kind=int_kind) :: ierr
    character(len=*), parameter :: subname = '(evp1d_alloc_static_na)'

    ! navel0+1: the last slot is the dead slot.  Directions that leave the
    ! grid point at it, so the kernels can read a neighbour unconditionally
    ! and get zero -- the same convention the compressed 2-D layout uses.
    allocate(str1(1:navel0+1)   , str2(1:navel0+1), str3(1:navel0+1), &
             str4(1:navel0+1)   , str5(1:navel0+1), str6(1:navel0+1), &
             str7(1:navel0+1)   , str8(1:navel0+1),                   &
             indxTij(1:navel0+1), uvel(1:navel0+1), vvel(1:navel0+1), &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

  end subroutine evp1d_alloc_static_navel

!=============================================================================

  subroutine evp1d_alloc_static_halo()

    implicit none
    integer(kind=int_kind) :: ierr
    character(len=*), parameter :: subname = '(evp1d_alloc_static_halo)'

    ! One entry per ghost cell per edge, for the cyclic pass and again for
    ! the zero_gradient/linear_extrap pass.
    allocate(halo_bc_dst(4*nghost*(nx+ny)), halo_bc_s1(4*nghost*(nx+ny)), &
             halo_bc_s2(4*nghost*(nx+ny)), halo_bc_w (4*nghost*(nx+ny)), &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

  end subroutine evp1d_alloc_static_halo

!=============================================================================

  subroutine calc_nActiveTU(Tmask,na0, Umask)

    ! Calculate number of active points with a given mask.

    implicit none
    logical(kind=log_kind), intent(in) :: Tmask(:,:)
    logical(kind=log_kind), optional, intent(in) :: Umask(:,:)
    integer(kind=int_kind), intent(out)  :: na0
    integer(kind=int_kind)              :: i,j
    character(len=*), parameter :: subname = '(calc_nActivceTU)'

    ! Same range as calc_2d_indices_init, which fills the index arrays this
    ! counts: ilo..ihi+1, one cell into the ghosts, where ihi is nx-nghost.
    ! The two must agree or nActive exceeds what gets filled and the tail of
    ! indxti/indxtj is left at zero.
    na0=0
    if (present(Umask)) then
       do i=1+nghost,nx-nghost+1
       do j=1+nghost,ny-nghost+1
          if ((Tmask(i,j)) .or. (Umask(i,j))) then
             na0=na0+1
          endif
       enddo
       enddo
    else
       do i=1+nghost,nx-nghost+1
       do j=1+nghost,ny-nghost+1
          if (Tmask(i,j)) then
             na0=na0+1
          endif
       enddo
       enddo
    endif

  end subroutine calc_nActiveTU

!=============================================================================

  subroutine set_skipMe(iceTmask, iceUmask,na0)

    implicit none

    logical(kind=log_kind), intent(in) :: iceTmask(:,:), iceUmask(:,:)
    integer(kind=int_kind), intent(in) :: na0
    integer(kind=int_kind)              :: iw, i, j, niw
    character(len=*), parameter :: subname = '(set_skipMe)'

    skipUcell=.false.
    skipTcell=.false.
    niw=0
    ! first count
    do iw=1, na0
      i = indxti(iw)
      j = indxtj(iw)
      if ( iceTmask(i,j) .or. iceUmask(i,j)) then
         niw=niw+1
         ! ⚠ The UNION, not iceTmask alone.  iceTmask comes from aice and
         ! iceUmask from aice_init -- two different snapshots -- so iceUmask
         ! is true at a few cells where iceTmask is false.  Dropping those
         ! would leave stepu never updating their velocity, silently.  When
         ! that artifact is eventually removed the union collapses to one
         ! mask, but removing it is not bit-for-bit, so not yet.
         if (evp1d_reorder) aList(niw) = iw
      endif
      if (.not. (iceTmask(i,j))) skipTcell(iw)=.true.
      if (.not. (iceUmask(i,j))) skipUcell(iw)=.true.
      ! No U point outside ilo..ihi, jlo..jhi, which is where the 2-D solver
      ! steps U.  ihi is nx-nghost, so every ghost column past it is skipped,
      ! not just the outermost one -- the same one-ghost assumption as above.
      if (i > nx - nghost) skipUcell(iw)=.true.
      if (j > ny - nghost) skipUcell(iw)=.true.
    enddo
    if (evp1d_reorder) then
       nIce = niw
       ! Once, on the first step: the compaction ratio is the whole point of
       ! the flag and it is entirely domain-dependent.  A basin that is mostly
       ! ice gains nothing and pays one extra integer stream; a global domain,
       ! where the ice-free tropics are most of the ocean, is where it pays.
       if (first_reorder) then
          first_reorder = .false.
          write(nu_diag,*) subname//' evp1d_reorder: ', niw, ' of ', na0, &
                           ' slots carry ice'
       endif
    else
       nIce = na0
    endif
    !    write(nu_diag,*) 'number of points and Active points', na0, niw

  end subroutine set_skipMe

!=============================================================================

  subroutine calc_2d_indices_init(Tmask)
    ! All points are active. Need to find neighbors.
    ! This should include de selection of u points.

    implicit none

    ! nx, ny
    logical(kind=log_kind), dimension(:,:), intent(in) :: Tmask

    ! local variables

    integer(kind=int_kind) :: i, j, Nmaskt
    character(len=*), parameter :: subname = '(calc_2d_indices_init)'

    indxti(:) = 0
    indxtj(:) = 0
    Nmaskt = 0
    ! The T cells the 2-D solver computes stress at: ilo..ihi+1, jlo..jhi+1,
    ! one cell into the ghosts.  ihi is nx-nghost, so the upper bound is
    ! nx-nghost+1 -- which is nx only when nghost is 1.  Running to nx took an
    ! extra ring of T cells for nghost > 1 and the two solvers diverged.
    do j = 1 + nghost, ny - nghost + 1
    do i = 1 + nghost, nx - nghost + 1
       if (Tmask(i,j)) then
          Nmaskt = Nmaskt + 1
          indxti(Nmaskt) = i
          indxtj(Nmaskt) = j
       end if
    end do
    end do

  end subroutine calc_2d_indices_init

  subroutine gather_static(G_uarear, G_dxT, G_dyT, G_DminTarea, G_fcor, G_Tmask)

     ! In standalone  distrb_info is an integer. Not needed anyway
     use ice_communicate, only : master_task
     use ice_gather_scatter, only : gather_global, gather_globalR8_compressed
     use ice_domain, only : distrb_info
     use ice_grid, only: dyT, dxT, uarear, tmask, &
             dynMap, dynOff, uarear_slot, dxT_slot, dyT_slot
     use ice_dyn_shared, only: DminTarea, fcor_blk, DminTarea_slot
     implicit none

     ! nx, ny
     real(kind=dbl_kind)   , dimension(:,:), intent(out) :: G_uarear, G_dxT, G_dyT, &
                                                            G_DminTarea, G_fcor
     logical(kind=log_kind), dimension(:,:), intent(out) :: G_Tmask

     character(len=*), parameter :: subname = '(gather_static)'

     G_uarear = c0
     G_dyT = c0
     G_dxT = c0
     G_DminTarea = c0
     G_fcor = c0
     G_tmask = .false.

     ! copy from distributed I_* to G_*
     ! ⚠ These four are METRICS -- nonzero over land -- so a compressed
     ! gather, which leaves c0 wherever a cell has no slot, is safe only
     ! because the 1-D kernels read them at the OWN cell and nowhere else:
     ! dxT(iw), dyT(iw), uarear(iw), DminTarea(iw).  The neighbour reads in
     ! those kernels are all on str1..str8.  Every cell whose metric is read
     ! is active, and active implies iceTmask implies tmask implies a slot.
     ! If a metric is ever read at a neighbour this must be revisited --
     ! that is the uarea bug which broke every real-grid gridc test.
     call gather_globalR8_compressed(G_uarear, uarear_slot, dynMap, dynOff, &
                 master_task, distrb_info, c0, grid_ext=.true.)
     call gather_globalR8_compressed(G_dxT, dxT_slot, dynMap, dynOff, &
                 master_task, distrb_info, c0, grid_ext=.true.)
     call gather_globalR8_compressed(G_dyT, dyT_slot, dynMap, dynOff, &
                 master_task, distrb_info, c0, grid_ext=.true.)
     ! DminTarea is deltamin*tarea, set once in init_dyn_shared.  It is gathered
     ! rather than rebuilt from dxT*dyT because tarea is only that product on the
     ! init_grid2 path; on a MOM supergrid it is a sum of four sub-areas.
     call gather_globalR8_compressed(G_DminTarea, DminTarea_slot, dynMap, dynOff, &
                 master_task, distrb_info, c0, grid_ext=.true.)
     ! Coriolis is constant, so it is a static like the metrics above.  With it
     ! here, fmU and umassdti need not be gathered every timestep: dyn_prep2
     ! forms them as fcor*umass and umass/dt, so umass alone yields both.
     ! ⚠ fcor and Tmask keep the block gather: fcor_blk has no slot form,
     ! and Tmask is logical where gather_globalR8_compressed is real.
     call gather_global(G_fcor     , fcor_blk , master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_Tmask , Tmask , master_task, distrb_info, grid_ext=.true.)

  end subroutine gather_static


!=============================================================================
! Gather only the twelve stresses.  Used once per run, to seed the 1d arrays
! from the block layout -- see stresses_seeded.  Separate from gather_dyn so
! the steady-state path carries no stress gathers at all.

  subroutine gather_stresses()

     use ice_dyn_shared, only : &
          stressp_1_slot, stressp_2_slot, stressp_3_slot, stressp_4_slot, &
          stressm_1_slot, stressm_2_slot, stressm_3_slot, stressm_4_slot, &
          stress12_1_slot, stress12_2_slot, stress12_3_slot, stress12_4_slot
     implicit none

     character(len=*), parameter :: subname = '(gather_stresses)'

     ! ⚠ Gathered from the SLOT copies, which have been the live stress
     ! state since 118d511 -- the block arrays are I/O staging now and are
     ! only current just after a restart read.  The seed runs once, on the
     ! first call, which is exactly when stresses_to_slots has filled them.
     ! Zero at a slotless cell is right here: a cell with no slot carries
     ! no stress in either layout.
     !
     ! Straight into the 1d arrays, so the twelve (nx,ny) automatics and the
     ! seed loop that copied out of them are both gone.
     call gather_slot_r8(stressp_1 , stressp_1_slot , nActive)
     call gather_slot_r8(stressp_2 , stressp_2_slot , nActive)
     call gather_slot_r8(stressp_3 , stressp_3_slot , nActive)
     call gather_slot_r8(stressp_4 , stressp_4_slot , nActive)
     call gather_slot_r8(stressm_1 , stressm_1_slot , nActive)
     call gather_slot_r8(stressm_2 , stressm_2_slot , nActive)
     call gather_slot_r8(stressm_3 , stressm_3_slot , nActive)
     call gather_slot_r8(stressm_4 , stressm_4_slot , nActive)
     call gather_slot_r8(stress12_1, stress12_1_slot, nActive)
     call gather_slot_r8(stress12_2, stress12_2_slot, nActive)
     call gather_slot_r8(stress12_3, stress12_3_slot, nActive)
     call gather_slot_r8(stress12_4, stress12_4_slot, nActive)

  end subroutine gather_stresses

!=============================================================================

  ! The twelve stresses are NOT gathered.  They are subcycle state that lives
  ! in the module-level 1d arrays and persists between timesteps, so gathering
  ! them every step only re-fetched what was already there -- 12 of the 46
  ! collectives this routine and scatter_dyn cost per timestep.  The one thing
  ! that did arrive by that route is dyn_prep2's zeroing of cells with no ice,
  ! and convert_2d_1d_dyn now does that directly on the 1d arrays instead.
  ! waterxU and wateryU are not gathered either: dyn_prep2 forms them as pure
  ! functions of uocn, vocn and fmU, all three of which arrive here anyway, so
  ! convert_2d_1d_dyn recomputes them from the same operands instead.
  ! Gather a compressed field and land it straight in its 1d array, through
  ! the single shared G_work instead of a global of its own.
  !
  ! This is the whole point of the change: the value written is bit-for-bit
  ! what convert_2d_1d_dyn used to copy out of a per-field G_ array, because
  ! indxTij(iw) decodes to exactly the (indxti(iw), indxtj(iw)) that copy
  ! indexed with.  What goes away is the storage -- one (nx,ny) array in
  ! place of 23 -- and on a production grid that is the difference between
  ! ~360 MB of stack per timestep and ~12 MB.  That stack footprint, not the
  ! collective count, is what has kept evp1d off real domains.
  !
  ! ⚠ n must be the 1d array's own extent.  Everything routed through here is
  ! allocated (1:nActive); uvel and vvel are navel+1 and are NOT, because they
  ! round-trip through a global of their own for the scatter back.
  subroutine gather_slot_r8(A_1d, A_slot, n)

     use ice_communicate, only : my_task, master_task
     use ice_gather_scatter, only : gather_globalR8_compressed
     use ice_grid, only : dynMap, dynOff
     use ice_domain, only : distrb_info
     implicit none

     ! ⚠ allocatable, and that is not decoration.  evp1d_alloc_static_na runs
     ! inside the master-only block of dyn_evp1d_init, so on every other task
     ! this array is UNALLOCATED -- and this routine is a collective, called
     ! on all of them.  Passing an unallocated allocatable to a plain
     ! assumed-shape dummy is illegal; to an allocatable dummy it is legal and
     ! well defined, which is what lets the early return below be the whole
     ! story.  NO_MPI cannot show this: with one task my_task IS master_task,
     ! so the array is always allocated and the return never fires.
     real(kind=dbl_kind), allocatable, dimension(:), intent(inout) :: A_1d
     real(kind=dbl_kind), dimension(:), intent(in)    :: A_slot
     integer(kind=int_kind), intent(in) :: n

     integer(kind=int_kind) :: iw, i, j

     character(len=*), parameter :: subname = '(gather_slot_r8)'

     call gather_globalR8_compressed(G_work, A_slot, dynMap, dynOff, &
            master_task, distrb_info, c0, grid_ext=.true.)

     ! indxTij, nActive and the 1d arrays exist on the master alone -- see
     ! dyn_evp1d_init, where everything after gather_static is master-only.
     if (my_task /= master_task) return

     do iw = 1, n
        j = int((indxTij(iw) - 1) / nx) + 1
        i =      indxTij(iw) - (j - 1) * nx
        A_1d(iw) = G_work(i,j)
     end do

  end subroutine gather_slot_r8

!=============================================================================

  ! The same, for a field that still has only a block form.  umass is the one
  ! left: fmU is fcor*umass and umassdti is umass/dt, so gathering it covers
  ! both.
  subroutine gather_blk_slot_r8(A_1d, L_blk, n)

     use ice_communicate, only : my_task, master_task
     use ice_gather_scatter, only : gather_global
     use ice_domain, only : distrb_info
     implicit none

     ! allocatable for the same reason as in gather_slot_r8
     real(kind=dbl_kind), allocatable, dimension(:), intent(inout) :: A_1d
     ! i and j flattened: element (i,j,iblk) at ((j-1)*nx_block + i, iblk)
     real(kind=dbl_kind), dimension(:,:), intent(in)    :: L_blk
     integer(kind=int_kind), intent(in) :: n

     integer(kind=int_kind) :: iw, i, j

     character(len=*), parameter :: subname = '(gather_blk_slot_r8)'

     call gather_global(G_work, L_blk, master_task, distrb_info, grid_ext=.true.)

     if (my_task /= master_task) return

     do iw = 1, n
        j = int((indxTij(iw) - 1) / nx) + 1
        i =      indxTij(iw) - (j - 1) * nx
        A_1d(iw) = G_work(i,j)
     end do

  end subroutine gather_blk_slot_r8

!=============================================================================

  ! The mirror of gather_slot_r8, for a field that goes out and is not read
  ! back in: zero the shared global, write this solver's cells into it, scatter.
  !
  ! Zeroing is right here and wrong for the round-trip fields.  A stress at a
  ! cell evp1d does not own is zero in both layouts -- that is what
  ! convert_1d_2d_dyn did to the twelve stress globals before scattering them.
  ! uvel, vvel and the four diagnostics are the opposite case: they must keep
  ! the value the gather brought in at those cells, which is why they still
  ! have globals of their own and do NOT come through here.
  subroutine scatter_slot_r8(L_blk, A_1d, n)

     use ice_communicate, only : my_task, master_task
     use ice_gather_scatter, only : scatter_global
     use ice_domain, only : distrb_info
     implicit none

     ! i and j flattened: element (i,j,iblk) at ((j-1)*nx_block + i, iblk)
     real(kind=dbl_kind), dimension(:,:), intent(out) :: L_blk
     ! allocatable for the same reason as in gather_slot_r8
     real(kind=dbl_kind), allocatable, dimension(:), intent(in) :: A_1d
     integer(kind=int_kind), intent(in) :: n

     integer(kind=int_kind) :: iw, i, j

     character(len=*), parameter :: subname = '(scatter_slot_r8)'

     if (my_task == master_task) then
        G_work = c0
        do iw = 1, n
           j = int((indxTij(iw) - 1) / nx) + 1
           i =      indxTij(iw) - (j - 1) * nx
           G_work(i,j) = A_1d(iw)
        end do
     endif

     call scatter_global(L_blk, G_work, master_task, distrb_info, grid_ext=.true.)

  end subroutine scatter_slot_r8

!=============================================================================

  subroutine gather_dyn(L_umass     , L_uvel      , L_vvel      ,               &
                        L_icetmask  , L_iceUmask  ,                             &
                        G_strintxU  , G_strintyU  , G_taubxU    , G_taubyU    , &
                        lstrint     , ltaub       ,                             &
                        G_uvel      , G_vvel      ,                             &
                        G_iceTmask,  G_iceUmask)
     use ice_communicate, only : master_task
     use ice_gather_scatter, only : gather_global, gather_globalR8_compressed
     use ice_grid, only : dynMap, dynOff
     use ice_dyn_shared, only : strintxU_slot, strintyU_slot, &
          taubxU_slot, taubyU_slot, TbU_slot, &
           strength_slot, cdn_ocnU_slot, aiU_slot, uocnU_slot, vocnU_slot, &
           forcexU_slot, forceyU_slot, uvel_slot, vvel_slot
     use ice_domain, only : distrb_info
     implicit none
     ! uvel/vvel flattened, as L_umass below
     real(kind=dbl_kind)   , dimension(:,:)  , intent(in)  :: &
        L_uvel      , L_vvel
     ! these have their i and j dimensions flattened, i.e. element (i,j,iblk)
     ! lives at ((j-1)*nx_block + i, iblk)
     real(kind=dbl_kind)   , dimension(:,:)  , intent(in)  :: &
        L_umass
     ! L_iceTmask and L_iceUmask have their i and j dimensions flattened,
     ! i.e. element (i,j,iblk) lives at ((j-1)*nx_block + i, iblk)
     logical(kind=log_kind), dimension(:,:), intent(in)  :: &
        L_iceTmask  , L_iceUmask
     ! the four diagnostics, gathered so the round trip leaves alone any cell
     ! this solver does not own -- see the note in convert_1d_2d_dyn
     real(kind=dbl_kind)   , dimension(:,:), intent(inout) :: &
        G_strintxU  , G_strintyU  , G_taubxU    , G_taubyU
     logical(kind=log_kind), intent(in) :: lstrint, ltaub
     ! nx, ny.  uvel and vvel keep a global of their own: unlike everything
     ! else here they go back out again, and scatter_dyn needs the gathered
     ! value still sitting at every cell this solver does not own.
     real(kind=dbl_kind)   , dimension(:,:), intent(out) ::     &
        G_uvel      , G_vvel
     logical(kind=log_kind), dimension(:,:), intent(out) :: &
        G_iceUmask  , G_iceTmask
     character(len=*), parameter :: subname = '(gather_dyn)'

     ! Straight into the 1d arrays.  These eleven used to be gathered into
     ! eleven (nx,ny) automatics and then copied out one cell at a time by
     ! convert_2d_1d_dyn; the copy is now the tail of gather_slot_r8 and the
     ! globals are gone.
     call gather_slot_r8(strength, strength_slot , nActive)
     call gather_slot_r8(cdn_ocn , cdn_ocnU_slot , nActive)
     call gather_slot_r8(aiu     , aiU_slot      , nActive)
     call gather_slot_r8(uocn    , uocnU_slot    , nActive)
     call gather_slot_r8(vocn    , vocnU_slot    , nActive)
     call gather_slot_r8(forcexU , forcexU_slot  , nActive)
     call gather_slot_r8(forceyU , forceyU_slot  , nActive)
     ! umass alone: fmU is fcor*umass and umassdti is umass/dt, and fcor is a
     ! static gathered once, so one field here replaces the two that were.
     call gather_blk_slot_r8(umass, L_umass      , nActive)
     ! TbU has no block form any more -- seabed_stress_factor_* write the
     ! slots directly -- so it is gathered the way strintxU already is.
     call gather_slot_r8(Tbu     , TbU_slot      , nActive)
     call gather_globalR8_compressed(G_uvel, uvel_slot, dynMap, dynOff, &
            master_task, distrb_info, c0, grid_ext=.true.)
     call gather_globalR8_compressed(G_vvel, vvel_slot, dynMap, dynOff, &
            master_task, distrb_info, c0, grid_ext=.true.)
     ! Only when history wants them, matching the scatter in scatter_dyn.
     if (lstrint) then
        call gather_globalR8_compressed(G_strintxU  , strintxU_slot, dynMap, dynOff, &
             master_task, distrb_info, c0, grid_ext=.true.)
        call gather_globalR8_compressed(G_strintyU  , strintyU_slot, dynMap, dynOff, &
             master_task, distrb_info, c0, grid_ext=.true.)
     endif
     if (ltaub) then
        call gather_globalR8_compressed(G_taubxU    , taubxU_slot, dynMap, dynOff, &
             master_task, distrb_info, c0, grid_ext=.true.)
        call gather_globalR8_compressed(G_taubyU    , taubyU_slot, dynMap, dynOff, &
             master_task, distrb_info, c0, grid_ext=.true.)
     endif
     call gather_global(G_iceTmask  ,     L_iceTmask  ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_iceUmask  ,     L_iceUmask  ,     master_task, distrb_info, grid_ext=.true.)
  end subroutine gather_dyn


!=============================================================================
! Scatter the twelve stresses back to the block layout.  Called only when a
! history or restart file is due: nothing outside the dynamics reads them in
! between -- ice_history touches them only inside its write_history gate, and
! the 1d arrays carry the state onward regardless.

  subroutine scatter_stresses(L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
                              L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
                              L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4)

     use ice_communicate, only : master_task
     use ice_gather_scatter, only : scatter_global
     use ice_domain, only : distrb_info
     implicit none

     real(kind=dbl_kind), dimension(:,:), intent(out) ::        &
        L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
        L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
        L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4

     character(len=*), parameter :: subname = '(scatter_stresses)'

     ! Straight out of the 1d arrays, each through the shared G_work, so the
     ! twelve (nx,ny) globals this used to need are gone.  convert_1d_2d_dyn
     ! zeroed them and filled them at the active cells every timestep even
     ! though only an I/O step ever scattered them; that work moves here, to
     ! the gate that actually wants it.
     call scatter_slot_r8(L_stressp_1 , stressp_1 , nActive)
     call scatter_slot_r8(L_stressp_2 , stressp_2 , nActive)
     call scatter_slot_r8(L_stressp_3 , stressp_3 , nActive)
     call scatter_slot_r8(L_stressp_4 , stressp_4 , nActive)
     call scatter_slot_r8(L_stressm_1 , stressm_1 , nActive)
     call scatter_slot_r8(L_stressm_2 , stressm_2 , nActive)
     call scatter_slot_r8(L_stressm_3 , stressm_3 , nActive)
     call scatter_slot_r8(L_stressm_4 , stressm_4 , nActive)
     call scatter_slot_r8(L_stress12_1, stress12_1, nActive)
     call scatter_slot_r8(L_stress12_2, stress12_2, nActive)
     call scatter_slot_r8(L_stress12_3, stress12_3, nActive)
     call scatter_slot_r8(L_stress12_4, stress12_4, nActive)

  end subroutine scatter_stresses

!=============================================================================

  ! The twelve stresses are NOT scattered here -- scatter_stresses does that,
  ! and only when a file is due.  strintxU, strintyU, taubxU and taubyU ARE,
  ! every timestep: ice_history accumulates them into a running average
  ! (accum_hist_field, before the write_history gate), so unlike the stresses
  ! they are consumed every step rather than only at write time.
  subroutine scatter_dyn(                         L_uvel      , L_vvel      ,                          &
                         G_strintxU  , G_strintyU  , G_uvel      , G_vvel      ,                          &
                         G_taubxU    , G_taubyU    , lstrint     , ltaub )

     use ice_communicate, only : master_task
     use ice_gather_scatter, only : scatter_global, scatter_globalR8_compressed
     use ice_grid, only : dynMap, dynOff, dynCnt, pack_to_slots
     use ice_dyn_shared, only : strintxU_slot, strintyU_slot, &
          taubxU_slot, taubyU_slot, uvel_slot, vvel_slot
     use ice_domain, only : distrb_info, nblocks
     implicit none

     ! nx_block, ny_block, max_blocks
     ! uvel/vvel flattened
     real(kind=dbl_kind), dimension(:,:)  , intent(out) :: &
        L_uvel      , L_vvel

     ! these have their i and j dimensions flattened, i.e. element (i,j,iblk)
     ! lives at ((j-1)*nx_block + i, iblk)

     ! nx, ny
     real(kind=dbl_kind), dimension(:,:), intent(in) ::         &
        G_strintxU  , G_strintyU  , G_uvel      , G_vvel      , &
        G_taubxU    , G_taubyU

     ! whether history wants the four diagnostics this step
     logical(kind=log_kind), intent(in) :: lstrint, ltaub

     integer(kind=int_kind) :: iblk

     character(len=*), parameter :: subname = '(scatter_dyn)'




     ! uvel and vvel always; the four diagnostics only when a history
     ! stream wants them, which is four collectives per timestep saved in
     ! the common case where it does not.
     if (lstrint) then
        call scatter_globalR8_compressed(strintxU_slot, G_strintxU  , dynMap, dynOff, &
             master_task, distrb_info, grid_ext=.true.)
        call scatter_globalR8_compressed(strintyU_slot, G_strintyU  , dynMap, dynOff, &
             master_task, distrb_info, grid_ext=.true.)
     endif
     call scatter_global(L_uvel      , G_uvel      , master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_vvel      , G_vvel      , master_task, distrb_info, grid_ext=.true.)
     if (ltaub) then
        call scatter_globalR8_compressed(taubxU_slot, G_taubxU    , dynMap, dynOff, &
             master_task, distrb_info, grid_ext=.true.)
        call scatter_globalR8_compressed(taubyU_slot, G_taubyU    , dynMap, dynOff, &
             master_task, distrb_info, grid_ext=.true.)
     endif

     ! The compressed velocities, which deformations reads after the solver
     ! returns.  They belong here rather than in the caller: this routine is
     ! what produces the velocities, so it is what must leave BOTH forms
     ! current.  Having the caller remember is how they came to hold zeros
     ! for a whole run.
     !
     ! ⚠ Not a compressed scatter.  Doing that would need the block forms
     ! unpacked back out for restart, history, transport and the couplers,
     ! and the cells with no slot would then be decided here rather than
     ! round-tripped through the gather unchanged -- which is exactly the
     ! ownership question 1e93e13 settled.  It would also move bytes, not
     ! collectives: the count is the same either way.  This is local traffic
     ! over two arrays and costs nothing on the wire.
     !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
     do iblk = 1, nblocks
        call pack_to_slots(L_uvel(:,iblk), dynMap(:,iblk), uvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
        call pack_to_slots(L_vvel(:,iblk), dynMap(:,iblk), vvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
     enddo
     !$OMP END PARALLEL DO

  end subroutine scatter_dyn

!=============================================================================

  subroutine convert_2d_1d_init(na0, G_HTE, G_HTN, G_uarear,  G_dxT, G_dyT, &
                                G_DminTarea, G_fcor)

     implicit none

     integer(kind=int_kind), intent(in) ::  na0
     real (kind=dbl_kind), dimension(:, :), intent(in) :: G_HTE, G_HTN, G_uarear, G_dxT, G_dyT, &
                                                          G_DminTarea, G_fcor

     ! local variables

     integer(kind=int_kind) :: iw, lo, up, j, i, ij, islot, ierr

     character(len=*), parameter :: subname = '(convert_2d_1d_init)'

     ! Hand out the slots.  The layout is the one the kernels rely on: the
     ! active T and U points first, 1..na0, in the order calc_2d_indices_init
     ! found them, then the cells that are only ever read.  stress_1d and
     ! stepu_1d both run over 1..na0, so nothing downstream has to know where
     ! the halo tail begins beyond na0 itself.
     do iw = 1, na0
        ij = indxti(iw) + (indxtj(iw) - 1) * nx
        ijslot(ij)  = iw
        indxTij(iw) = ij
     end do

     ! The read-only cells take the slots after them, in ascending linear
     ! index -- the same order the sorted union used to produce.
     islot = na0
     do ij = 1, nx * ny
        if (ijslot(ij) == -1) then
           islot = islot + 1
           ijslot(ij)     = islot
           indxTij(islot) = ij
        endif
     end do
     if (islot /= navel) then
        call abort_ice(subname//' ERROR: slot count does not match navel', &
                       file=__FILE__, line=__LINE__)
     endif

     ! Neighbour slots.  A direction that leaves the grid takes the dead slot
     ! at navel+1, which holds zero and is never written.  That only happens
     ! for nbrE/nbrNE/nbrN at i == nx or j == ny, where the cell has no U point to
     ! step -- skipUcell is already .true. there, matching the 2-D solver,
     ! which computes stress over jlo..jhi+1 but steps U only over jlo..jhi.
     ! nbrW/nbrSW/nbrS step towards i-1 and j-1 and active T cells start at
     ! 1+nghost, so for nghost >= 1 they always land inside the grid.
     do iw = 1, na0
        i = indxti(iw)
        j = indxtj(iw)
        nbrW (iw) = slot_of(i - 1, j    , navel + 1)
        nbrSW(iw) = slot_of(i - 1, j - 1, navel + 1)
        nbrS (iw) = slot_of(i    , j - 1, navel + 1)
        nbrE (iw) = slot_of(i + 1, j    , navel + 1)
        nbrNE(iw) = slot_of(i + 1, j + 1, navel + 1)
        nbrN (iw) = slot_of(i    , j + 1, navel + 1)
     end do

     !tar      i$OMP PARALLEL PRIVATE(iw, lo, up, j, i)
     ! write 1D arrays from 2D arrays (target points)
     !tar      call domp_get_domain(1, na0, lo, up)
     lo=1
     up=na0
     do iw = 1, na0
        ! get 2D indices
        i = indxti(iw)
        j = indxtj(iw)
        ! map
        uarear_1d(iw)  = G_uarear(i, j)
        dxT_1d(iw)     = G_dxT(i, j)
        dyT_1d(iw)     = G_dyT(i, j)
        DminTarea_1d(iw) = G_DminTarea(i, j)
        fcor_1d(iw)      = G_fcor(i, j)
        HTE_1d(iw)     = G_HTE(i, j)
        HTN_1d(iw)     = G_HTN(i, j)
        HTEm1_1d(iw)   = G_HTE(i - 1, j)
        HTNm1_1d(iw)   = G_HTN(i, j - 1)
        ! stress_1d used to form these inline each subcycle.  Building them
        ! here from the same gathered grid keeps every value identical while
        ! letting the kernel take them as arguments, the way the 2d one does.
        !
        ! NOT gathered from ice_dyn_shared's dxhy/dyhx.  Both the allocation
        ! (ice_dyn_shared.F90:465) and the computation (:914) are guarded by
        ! evp_algorithm == "standard_2d", so on this path those arrays are
        ! never ALLOCATED -- gathering them is not "reads zeros", it is a
        ! reference to an unallocated allocatable.  Where they do exist they
        ! also differ from what is built here: computed over ilo:ihi and then
        ! halo-updated with fillValue=c1 (:930), so at a true domain edge they
        ! hold 1.0 rather than anything derived from HTE/HTN.
        !
        ! Both of those are upstream behaviour, unchanged at 12181ec.  Do not
        ! "fix" them to make one canonical copy serve both solvers: the
        ! fillValue would change 2d EVP results at domain-edge cells and move
        ! the branch off bit-for-bit against the baseline.
        dxhy_1d(iw)    = p5 * (G_HTE(i, j) - G_HTE(i - 1, j))
        dyhx_1d(iw)    = p5 * (G_HTN(i, j) - G_HTN(i, j - 1))
     end do

  end subroutine convert_2d_1d_init

!=============================================================================

  ! What is left here after the gathers started landing straight in the 1d
  ! arrays: the per-cell work that is NOT a copy.  The zeroing that depends on
  ! skipTcell/skipUcell -- which gather_dyn runs too early to know, since
  ! set_skipMe comes after it -- and the four quantities derived from fields
  ! gathered above.  uvel and vvel are still copied, because they alone keep a
  ! global of their own for the scatter back.
  subroutine convert_2d_1d_dyn(na0         , dt          , G_uvel      , G_vvel      )
     use ice_dyn_shared, only: cosw, sinw
     implicit none
     integer(kind=int_kind), intent(in) ::  na0
     real(kind=dbl_kind), intent(in) :: dt   ! dynamics timestep
     ! nx, ny
     real(kind=dbl_kind), dimension(:, :), intent(in) :: &
        G_uvel      , G_vvel
     integer(kind=int_kind) ::  lo, up, iw, i, j
     character(len=*), parameter :: subname = '(convert_2d_1d_dyn)'
     lo=1
     up=na0
     do iw = 1, na0
        ! get 2D indices
        i = indxti(iw)
        j = indxtj(iw)
        ! The stresses are not copied in: they persist in the 1d arrays, which
        ! is why gather_dyn no longer gathers them.  What did arrive that way is
        ! dyn_prep2's zeroing of cells carrying no ice, so do it here instead --
        ! skipTcell is exactly dyn_prep2's .not. iceTmask condition.  Without
        ! this a cell that loses its ice would keep its stress and resume from
        ! it on re-icing.
        if (skipTcell(iw)) then
           stressp_1 (iw) = c0
           stressp_2 (iw) = c0
           stressp_3 (iw) = c0
           stressp_4 (iw) = c0
           stressm_1 (iw) = c0
           stressm_2 (iw) = c0
           stressm_3 (iw) = c0
           stressm_4 (iw) = c0
           stress12_1(iw) = c0
           stress12_2(iw) = c0
           stress12_3(iw) = c0
           stress12_4(iw) = c0
        endif
        ! strength, cdn_ocn, aiu, uocn, vocn, forcexU, forceyU, umass and Tbu
        ! were copied out of their globals here.  gather_dyn now writes them
        ! directly, with the same indxTij -> (i,j) decode, so the values are
        ! unchanged and only the storage differs.
        !
        ! The rule dyn_prep2 uses: 1, or 0 where the U cell carries less ice
        ! than rheo_area_min.  #1133 added this factor to the 2-D stepu but
        ! not to this kernel, so the two solvers disagreed on any grid with
        ! marginal ice -- gx3 shows it, a solid-ice box grid does not.
        if (free_drift_approx .and. aiu(iw) <= rheo_area_min) then
           rheofactU(iw) = c0
        else
           rheofactU(iw) = c1
        endif
        strintxU_1d(iw)   = C0
        strintyU_1d(iw)   = C0
        Cb(iw)         = c0
        uvel(iw)       = G_uvel(i,j)
        vvel(iw)       = G_vvel(i,j)
        uvel_init(iw)  = G_uvel(i,j)
        vvel_init(iw)  = G_vvel(i,j)

        ! Derived from the fields copied in just above, so this MUST follow
        ! them: fmU and umassdti from umass and the static fcor, then waterxU
        ! and wateryU from uocn, vocn and that fmU.  Same expressions and
        ! operands as dyn_prep2, so bit-for-bit:
        !   fm       = fcor*umass          umassdti = umass/dt
        !   waterx   = uocn*cosw - vocn*sinw*sign(c1,fm)
        !   watery   = vocn*cosw + uocn*sinw*sign(c1,fm)
        ! dyn_prep2 fills all four only at its active U cells, so the inactive
        ! ones read zero here.  It does not zero fm itself, so a block-layout
        ! inactive cell keeps a stale value there -- but nothing reads fm at a
        ! skipped cell, since stepu cycles on skipUcell, so zeroing is both
        ! safe and clearer than reproducing a staleness no one observes.
        if (skipUcell(iw)) then
           fmU     (iw) = c0
           umassdti(iw) = c0
           waterxU (iw) = c0
           wateryU (iw) = c0
        else
           fmU     (iw) = fcor_1d(iw) * umass(iw)
           umassdti(iw) = umass(iw) / dt
           waterxU (iw) = uocn(iw)*cosw - vocn(iw)*sinw*sign(c1,fmU(iw))
           wateryU (iw) = vocn(iw)*cosw + uocn(iw)*sinw*sign(c1,fmU(iw))
        endif
     end do
     ! Halos can potentially have values of u and v
     do iw=na0+1,navel
        j = int((indxTij(iw) - 1) / (nx)) + 1
        i =      indxTij(iw) - (j - 1) * nx
        uvel(iw)=G_uvel(i,j)
        vvel(iw)=G_vvel(i,j)
     end do
  end subroutine convert_2d_1d_dyn

!=============================================================================

  subroutine convert_1d_2d_dyn(na0         , navel0      ,                             &
                               G_strintxU  , G_strintyU  , &
                               G_uvel      , G_vvel      , G_taubxU    , &
                               G_taubyU)

     implicit none

     integer(kind=int_kind), intent(in) ::  na0, navel0
     ! nx, ny.  Only the round-trip fields are left: the ones that have to
     ! keep, at every cell this solver does not own, the value the gather
     ! brought in.  The twelve stresses went out through scatter_slot_r8
     ! instead, and ten more globals -- G_strength, G_cdn_ocn, G_aiu, G_uocn,
     ! G_vocn, G_waterxU, G_wateryU, G_forcexU, G_forceyU, G_Tbu -- turned out
     ! to be zeroed here and read by nobody at all.
     real(kind=dbl_kind), dimension(:, :), intent(inout) ::     &
        G_strintxU  , G_strintyU  , &
        G_uvel      , G_vvel      , G_taubxU    , &
        G_taubyU

     integer(kind=int_kind) ::  lo, up, iw, i, j
     character(len=*), parameter :: subname = '(convert_1d_2d_dyn)'

     ! G_strintxU, G_strintyU, G_taubxU and G_taubyU are deliberately NOT
     ! zeroed, for the same reason as G_uvel and G_vvel below: gather_dyn has
     ! filled them with the block field, and scatter_dyn writes them straight
     ! back.  The 2-D solver's unpack_from_slots writes only the cells it has
     ! a slot for and leaves the rest, so zeroing here made evp1d clear
     ! diagnostics at cells it does not own where the 2-D solver does not.
     ! G_uvel and G_vvel are deliberately NOT zeroed.  Unlike the others they
     ! arrive holding the gathered block field -- gather_dyn fills them with
     ! grid_ext, ghosts included -- and they are scattered straight back.  A
     ! cell with no slot is one this solver does not own, typically a ghost
     ! whose interior neighbour is land; zeroing it here made the scatter
     ! overwrite the block value with zero, so evp1d silently cleared
     ! velocities that the 2-D solver leaves alone.  Leaving them means the
     ! scatter writes back what it read at every cell the solver does not
     ! touch.  Both loops below write every slot, so nothing stale survives
     ! where the solver does own the cell.

     lo=1
     up=na0
     do iw = lo, up
        ! get 2D indices
        i = indxti(iw)
        j = indxtj(iw)
        ! map to 2d
        G_strintxU(i,j)   = strintxU_1d(iw)
        G_strintyU(i,j)   = strintyU_1d (iw)
        G_taubxU(i,j)     = -uvel(iw)*Cb(iw)
        G_taubyU(i,j)     = -vvel(iw)*Cb(iw)
        G_uvel(i,j)       = uvel(iw)
        G_vvel(i,j)       = vvel(iw)
     end do

     do iw=na0+1,navel0
        j = int((indxTij(iw) - 1) / (nx)) + 1
        i = indxTij(iw) - (j - 1) * nx
        G_uvel(i,j)       = uvel(iw)
        G_vvel(i,j)       = vvel(iw)
     end do

  end subroutine convert_1d_2d_dyn

  subroutine calc_navel(na0, navel0)
     ! Calculate number of active points, including halo points.
     !
     ! Marks every cell the 1-D solver touches -- each active T cell, plus the
     ! six neighbours its kernels read -- and counts them.  The marks stay in
     ! ijslot for convert_2d_1d_init, which turns them into slot numbers.
     !
     ! This replaces a sorted-merge union of seven index vectors built from
     ! i + (j-1)*nx arithmetic.  That arithmetic has no way to say "there is no
     ! such neighbour": at i == nx the (+1,.) directions wrapped round into the
     ! next row, and at j == ny they ran off the end of the grid entirely, so
     ! the union carried indices naming cells that do not exist.  A scan with
     ! explicit bounds guards cannot generate them.  It is how set_neighbours
     ! builds dynNbr for the compressed 2-D layout, and the point of doing it
     ! the same way here is that the two solvers then agree by construction
     ! rather than by two separate pieces of index arithmetic agreeing.

     implicit none

     integer(kind=int_kind), intent(in) :: na0
     integer(kind=int_kind), intent(out) :: navel0

     ! local variables

     integer(kind=int_kind) :: iw, i, j, ij, ierr

     character(len=*), parameter :: subname = '(calc_navel)'

     allocate(ijslot(1:nx*ny), stat=ierr)
     if (ierr/=0) then
        call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
     endif

     ijslot(:) = 0

     do iw = 1, na0
        i = indxti(iw)
        j = indxtj(iw)
        ! No trailing backslashes in these comments: .F90 is preprocessed and
        ! cpp splices a line ending in one into the comment above it, which
        ! silently deletes the call.
        call mark_cell(i    , j    )   ! ( 0,  0) this cell
        ! the U points stress_1d reads
        call mark_cell(i - 1, j    )   ! (-1,  0)
        call mark_cell(i - 1, j - 1)   ! (-1, -1)
        call mark_cell(i    , j - 1)   ! ( 0, -1)
        ! the T cells stepu_1d reads
        call mark_cell(i + 1, j    )   ! (+1,  0)
        call mark_cell(i + 1, j + 1)   ! (+1, +1)
        call mark_cell(i    , j + 1)   ! ( 0, +1)
     end do

     ! the halo update reads and writes cells of its own, beyond anything the
     ! kernels reach -- see halo_sweep
     call halo_sweep(.true.)

     navel0 = 0
     do ij = 1, nx * ny
        if (ijslot(ij) /= 0) navel0 = navel0 + 1
     end do

  end subroutine calc_navel

!=======================================================================

  subroutine mark_cell(i, j)
     ! Mark (i,j) as a cell the solver needs, if it is inside the gathered
     ! grid.  Outside means the neighbour does not exist, which is precisely
     ! what the index arithmetic this replaces could not express.

     implicit none

     integer(kind=int_kind), intent(in) :: i, j

     character(len=*), parameter :: subname = '(mark_cell)'

     if (i < 1 .or. i > nx) return
     if (j < 1 .or. j > ny) return
     ijslot(i + (j - 1) * nx) = -1

  end subroutine mark_cell

!=======================================================================

  integer(kind=int_kind) function slot_of(i, j, dead)
     ! Slot holding cell (i,j), or the dead slot if it is off the grid.

     implicit none

     integer(kind=int_kind), intent(in) :: i, j, dead

     character(len=*), parameter :: subname = '(slot_of)'

     if (i < 1 .or. i > nx .or. j < 1 .or. j > ny) then
        slot_of = dead
     else
        slot_of = ijslot(i + (j - 1) * nx)
     endif

  end function slot_of

!=======================================================================

  subroutine numainit(lo,up,uu)

     implicit none
     integer(kind=int_kind),intent(in) :: lo,up,uu
     integer(kind=int_kind) :: iw
     character(len=*), parameter :: subname = '(numainit)'

     !$omp parallel do schedule(runtime) private(iw)
     do iw = lo,up
        skipTcell(iw)=.false.
        skipUcell(iw)=.false.
        nbrW(iw)=0
        nbrSW(iw)=0
        nbrS(iw)=0
        nbrE(iw)=0
        nbrNE(iw)=0
        nbrN(iw)=0
        aiu(iw)=c0
        rheofactU(iw)=c1
        Cb(iw)=c0
        cdn_ocn(iw)=c0
        dxT_1d(iw)=c0
        dyT_1d(iw)=c0
        DminTarea_1d(iw)=c0
        dxhy_1d(iw)=c0
        dyhx_1d(iw)=c0
        fmU(iw)=c0
        forcexU(iw)=c0
        forceyU(iw)=c0
        HTE_1d(iw)=c0
        HTEm1_1d(iw)=c0
        HTN_1d(iw)=c0
        HTNm1_1d(iw)=c0
        strength(iw)= c0
        stress12_1(iw)=c0
        stress12_2(iw)=c0
        stress12_3(iw)=c0
        stress12_4(iw)=c0
        stressm_1(iw)=c0
        stressm_2(iw)=c0
        stressm_3(iw)=c0
        stressm_4(iw)=c0
        stressp_1(iw)=c0
        stressp_2(iw)=c0
        stressp_3(iw)=c0
        stressp_4(iw)=c0
        strintxU_1d(iw)= c0
        strintyU_1d(iw)= c0
        Tbu(iw)=c0
        uarear_1d(iw)=c0
        umassdti(iw)=c0
        uocn(iw)=c0
        uvel_init(iw)=c0
        uvel(iw)=c0
        vocn(iw)=c0
        vvel_init(iw)=c0
        vvel(iw)=c0
        waterxU(iw)=c0
        wateryU(iw)=c0
     enddo
     !$omp end parallel do
     !$omp parallel do schedule(runtime) private(iw)
     do iw = lo,uu
        uvel(iw)=c0
        vvel(iw)=c0
        str1(iw)=c0
        str2(iw)=c0
        str3(iw)=c0
        str4(iw)=c0
        str5(iw)=c0
        str6(iw)=c0
        str7(iw)=c0
        str8(iw)=c0
     enddo
     !$omp end parallel do

  end subroutine numainit

!=======================================================================

  subroutine halo_sweep(domark)
     ! Walk every ghost point the boundary conditions fill, in the order the
     ! 2-D halo fills them.
     !
     ! One sweep, two uses.  With domark it marks the cells the halo update
     ! reads and writes, so calc_navel gives them slots; without it, it builds
     ! the list.  Sharing the loop is the point: a cell the build pass needs
     ! but the mark pass missed would drop out of the list with no sign of it.
     !
     ! Order follows ice_boundary.F90.  The cyclic exchange runs first in both
     ! directions, then zero_gradient/linear_extrap east/west over the full
     ! column and north/south over the full row.  Running each pass the full
     ! length of the edge is what fills the corners: the north/south pass
     ! reads a ghost the east/west pass has just written, so a bi-cyclic
     ! corner ends up with the diagonal cell, which is what the 2-D exchange
     ! gets from its diagonal neighbour.  The list is therefore ordered, and
     ! the loop applying it must stay sequential.

     use ice_blocks, only: ew_boundary_type, ns_boundary_type

     implicit none

     logical(kind=log_kind), intent(in) :: domark

     ! local variables

     integer(kind=int_kind) :: i, j, g, ilo, ihi, jlo, jhi
     real   (kind=dbl_kind) :: w

     character(len=*), parameter :: subname = '(halo_sweep)'

     if (.not. domark) n_halo_bc = 0

     ! the gathered field is one block, so its physical domain is the whole
     ! grid less the ghost rim
     ilo = 1  + nghost
     ihi = nx - nghost
     jlo = 1  + nghost
     jhi = ny - nghost

     ! the cyclic exchange, both directions
     if (trim(ew_boundary_type) == 'cyclic') then
        do j = 1, ny
        do g = 1, nghost
           call halo_point(g      , j, ihi-nghost+g, j, ihi-nghost+g, j, c0, domark)
           call halo_point(ihi + g, j, ilo+g-1     , j, ilo+g-1     , j, c0, domark)
        end do
        end do
     endif

     if (trim(ns_boundary_type) == 'cyclic') then
        do i = 1, nx
        do g = 1, nghost
           call halo_point(i, g      , i, jhi-nghost+g, i, jhi-nghost+g, c0, domark)
           call halo_point(i, jhi + g, i, jlo+g-1     , i, jlo+g-1     , c0, domark)
        end do
        end do
     endif

     ! then the extrapolating types, east/west before north/south
     if (trim(ew_boundary_type) == 'zero_gradient' .or. &
         trim(ew_boundary_type) == 'linear_extrap') then
        do j = 1, ny
        do g = 1, nghost
           w = c0
           if (trim(ew_boundary_type) == 'linear_extrap') w = real(nghost-g+1, dbl_kind)
           call halo_point(g      , j, ilo, j, ilo+1, j, w, domark)   ! west
           w = c0
           if (trim(ew_boundary_type) == 'linear_extrap') w = real(g, dbl_kind)
           call halo_point(ihi + g, j, ihi, j, ihi-1, j, w, domark)   ! east
        end do
        end do
     endif

     if (trim(ns_boundary_type) == 'zero_gradient' .or. &
         trim(ns_boundary_type) == 'linear_extrap') then
        do i = 1, nx
        do g = 1, nghost
           w = c0
           if (trim(ns_boundary_type) == 'linear_extrap') w = real(nghost-g+1, dbl_kind)
           call halo_point(i, g      , i, jlo, i, jlo+1, w, domark)   ! south
           w = c0
           if (trim(ns_boundary_type) == 'linear_extrap') w = real(g, dbl_kind)
           call halo_point(i, jhi + g, i, jhi, i, jhi-1, w, domark)   ! north
        end do
        end do
     endif

  end subroutine halo_sweep

!=======================================================================

  subroutine halo_point(id, jd, i1, j1, i2, j2, w, domark)
     ! Mark, or append, one halo point.
     !
     ! Marking is what lets the append succeed: these cells are the halo
     ! update's own inputs and outputs, and a ghost whose interior neighbour
     ! is land is reached by no kernel, so without marking it has no slot and
     ! the point would be dropped.  The 2-D halo fills it regardless, which is
     ! where evp1d used to differ.

     implicit none

     integer(kind=int_kind), intent(in) :: id, jd, i1, j1, i2, j2
     real   (kind=dbl_kind), intent(in) :: w
     logical(kind=log_kind), intent(in) :: domark

     ! local variables

     integer(kind=int_kind) :: sd, s1, s2

     character(len=*), parameter :: subname = '(halo_point)'

     if (domark) then
        call mark_cell(id, jd)
        call mark_cell(i1, j1)
        if (w /= c0) call mark_cell(i2, j2)
        return
     endif

     sd = ijslot(id + (jd - 1) * nx)
     s1 = ijslot(i1 + (j1 - 1) * nx)
     if (sd == 0 .or. s1 == 0) return

     s2 = ijslot(i2 + (j2 - 1) * nx)
     if (s2 == 0) then
        ! w == 0 does not read s2, so a missing one is no obstacle
        if (w /= c0) return
        s2 = s1
     endif

     n_halo_bc = n_halo_bc + 1
     halo_bc_dst(n_halo_bc) = sd
     halo_bc_s1 (n_halo_bc) = s1
     halo_bc_s2 (n_halo_bc) = s2
     halo_bc_w  (n_halo_bc) = w

  end subroutine halo_point

!=======================================================================

  subroutine evp1d_halo_update()

     implicit none
     integer(kind=int_kind) :: iw

     character(len=*), parameter :: subname = '(evp1d_halo_update)'

! One sequential pass over the list built in halo_sweep.  Not parallelised
! and not to be: north/south destinations are east/west sources at the
! corners, so the entries are order-dependent.  See halo_sweep.
     do iw = 1, n_halo_bc
        uvel(halo_bc_dst(iw)) = uvel(halo_bc_s1(iw)) + halo_bc_w(iw) *      &
                               (uvel(halo_bc_s1(iw)) - uvel(halo_bc_s2(iw)))
        vvel(halo_bc_dst(iw)) = vvel(halo_bc_s1(iw)) + halo_bc_w(iw) *      &
                               (vvel(halo_bc_s1(iw)) - vvel(halo_bc_s2(iw)))
     end do

  end subroutine evp1d_halo_update

end module ice_dyn_evp1d

