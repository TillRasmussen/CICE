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
  use ice_blocks, only: nx_block, ny_block, nghost, ns_boundary_type
  use ice_constants
  use ice_communicate, only: my_task, master_task
  use ice_domain_size, only: max_blocks, nx_global, ny_global
  use ice_fileunits, only: nu_diag
  use ice_exit, only: abort_ice
  use ice_dyn_shared, only: rheo_area_min, free_drift_approx

  !- directives ----------------------------------------------------------------
  implicit none
  private

  !- public routines -----------------------------------------------------------
  public :: dyn_evp1d_init, dyn_evp1d_run, dyn_evp1d_finalize

  !- private routines ----------------------------------------------------------

  !- private vars --------------------------------------------------------------
  ! nx and ny are module variables for arrays after gather (G_*) Dimension according to CICE is
  ! nx_global+2*nghost, ny_global+2*nghost
  ! nactive are number of active points (both t and u). nStored_TU_glob is number of active
  integer(kind=int_kind), save :: nx, ny, nActive_TU_glob, nStored_TU_glob, nallocated

  ! indexes
  integer(kind=int_kind), allocatable, dimension(:,:) :: iwidx
  logical(kind=log_kind), allocatable, dimension(:)   :: skipTcell,skipUcell
  integer(kind=int_kind), allocatable, dimension(:)   :: ee,ne,se,nw,sw,sse ! arrays for neighbour points
  integer(kind=int_kind), allocatable, dimension(:)   :: indxti, indxtj, indxTij
  ! Scratch map over the gathered grid, linear index i+(j-1)*nx -> slot.
  ! -1 marks a cell the solver needs before slots are handed out, 0 a cell it
  ! does not. Built by calc_nStoredTU, consumed and released by convert_2d_1d_init,
  ! so it exists only during init and only on the master task.
  integer(kind=int_kind), allocatable, dimension(:)   :: ijslot

  ! 1D arrays to allocate

  ! Grid
  real   (kind=dbl_kind), allocatable, dimension(:)   :: &
     HTE_1d,HTN_1d, HTEm1_1d,HTNm1_1d, dxT_1d, dyT_1d, uarear_1d, &
     DminTarea_1d

  ! time varying
  real(kind=dbl_kind)   , allocatable, dimension(:)   ::                  &
    cdn_ocn,aiu,uocn,vocn,waterxU,wateryU,forcexU,forceyU,umassdti,fmU,   &
    rheofactU,                                                            &
    strintxU,strintyU,uvel_init,vvel_init, strength, uvel, vvel,          &
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
  real   (kind=dbl_kind), allocatable, dimension(:)   :: halo_bc_sgn
  integer(kind=int_kind), save :: n_tri_pair = 0, n_tri_neg = 0
  integer(kind=int_kind), allocatable, dimension(:) :: tri_a, tri_b, tri_neg
  integer(kind=int_kind), allocatable, dimension(:) :: tst_d, tst_s
  integer(kind=int_kind), save :: n_tri_str = 0

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
                                                           G_DminTarea
    logical(kind=log_kind), allocatable, dimension(:,:) :: G_tmask

    integer(kind=int_kind) :: ierr

    character(len=*), parameter :: subname = '(dyn_evp1d_init)'

    nx=nx_global+2*nghost
    ny=ny_global+2*nghost

    allocate(G_dyT(nx,ny),G_dxT(nx,ny),G_uarear(nx,ny),G_DminTarea(nx,ny), &
             G_tmask(nx,ny),stat=ierr)
    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

    ! gather from blks to global
    call gather_static(G_uarear, G_dxT, G_dyT, G_DminTarea, G_tmask)

    ! calculate number of water points (T and U). Only needed for the static version
    ! tmask in ocean/ice
    if (my_task == master_task) then
      call calc_nActiveTU(G_tmask)
      call evp1d_alloc_static_nActive
      call calc_2d_indices_init(G_tmask)
      call calc_nStoredTU
      call evp1d_alloc_static_nStored
      call numainit(1,nActive_TU_glob,nStored_TU_glob)
      call convert_2d_1d_init(G_HTE, G_HTN, G_uarear, G_dxT, G_dyT, G_DminTarea)
      call evp1d_alloc_static_halo()
      ! The halo geometry -- indxTij, nActive_TU_glob, nStored_TU_glob and the boundary types -- does
      ! not change, so the list is built once here.
      call halo_sweep(.false.)
      ! ijslot has served its purpose: slots are assigned and both halo lists
      ! are built.
      deallocate(ijslot, stat=ierr)
      if (ierr/=0) then
         call abort_ice(subname//' ERROR: deallocating', file=__FILE__, line=__LINE__)
      endif
    endif

    deallocate(G_dyT,G_dxT,G_uarear,G_DminTarea,G_tmask,stat=ierr)
    if (ierr/=0) then
       call abort_ice(subname//' ERROR: deallocating', file=__FILE__, line=__LINE__)
    endif

  end subroutine dyn_evp1d_init

!=============================================================================

  subroutine dyn_evp1d_run(L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
                           L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
                           L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4, &
                           L_strength,                                             &
                           L_cdn_ocn   , L_aiu       , L_uocn      , L_vocn      , &
                           L_waterxU   , L_wateryU   , L_forcexU   , L_forceyU   , &
                           L_umassdti  , L_fmU       , L_strintxU  , L_strintyU  , &
                           L_Tbu       , L_taubxU    , L_taubyU    , L_uvel      , &
                           L_vvel      , L_icetmask  , L_iceUmask)

    use ice_dyn_shared, only : ndte
    use ice_dyn_core1d, only : stress_1d, stepu_1d, calc_diag_1d
    use ice_timers    , only : ice_timer_start, ice_timer_stop, timer_evp1dcore

    use icepack_intfc , only : icepack_query_parameters, icepack_warnings_flush, &
      icepack_warnings_aborted

    implicit none

    ! nx_block, ny_block, max_blocks
    real(kind=dbl_kind)   , dimension(:,:,:), intent(inout) :: &
      L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 ,  &
      L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 ,  &
      L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4,  &
      L_strintxU  , L_strintyU  , L_uvel      , L_vvel      ,  &
      L_taubxU    , L_taubyU
    real(kind=dbl_kind)   , dimension(:,:,:), intent(in) ::    &
      L_strength  ,                                            &
      L_cdn_ocn   , L_aiu       , L_uocn     , L_vocn   ,      &
      L_waterxU   , L_wateryU   , L_forcexU  , L_forceyU,      &
      L_umassdti  , L_fmU       , L_Tbu
    logical(kind=log_kind), dimension(:,:,:), intent(in) ::    &
      L_iceUmask  , L_iceTmask

    ! local variables

    ! nx, ny
    real(kind=dbl_kind),    dimension(nx,ny) ::                &
      G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 ,  &
      G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 ,  &
      G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4,  &
      G_strength,                                              &
      G_cdn_ocn   , G_aiu       , G_uocn      , G_vocn      ,  &
      G_waterxU   , G_wateryU   , G_forcexU   , G_forceyU   ,  &
      G_umassdti  , G_fmU       , G_strintxU  , G_strintyU  ,  &
      G_Tbu       , G_uvel     , G_vvel       , G_taubxU    ,  &
      G_taubyU                  ! G_taubxU and G_taubyU are post processed from Cb
    logical(kind=log_kind), dimension (nx,ny) ::                 &
      G_iceUmask  , G_iceTmask

    logical(kind=log_kind) :: ltripole
    character(len=*), parameter :: subname = '(dyn_evp1d_run)'

    integer(kind=int_kind) :: ksub

    real   (kind=dbl_kind) :: rhow

    ! From 3d to 2d on master task
    call gather_dyn(L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
                    L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
                    L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4, &
                    L_strength,                                             &
                    L_cdn_ocn   , L_aiu       , L_uocn      , L_vocn      , &
                    L_waterxU   , L_wateryU   , L_forcexU   , L_forceyU   , &
                    L_umassdti  , L_fmU       ,                             &
                    L_Tbu       , L_uvel      , L_vvel      ,               &
                    L_icetmask  , L_iceUmask  ,                             &
                    L_strintxU  , L_strintyU  , L_taubxU    , L_taubyU    , &
                    G_strintxU  , G_strintyU  , G_taubxU    , G_taubyU    , &
                    G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
                    G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
                    G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
                    G_strength  ,                                           &
                    G_cdn_ocn   , G_aiu       , G_uocn      , G_vocn      , &
                    G_waterxU   , G_wateryU   , G_forcexU   , G_forceyU   , &
                    G_umassdti  , G_fmU       ,                             &
                    G_Tbu       , G_uvel      , G_vvel      ,               &
                    G_iceTmask,  G_iceUmask)

    if (my_task == master_task) then
       call set_skipMe(G_iceTmask, G_iceUmask)
       ! Map from 2d to 1d
       call convert_2d_1d_dyn(                                                                   &
                              G_stressp_1 , G_stressp_2 , G_stressp_3 ,  G_stressp_4,     &
                              G_stressm_1 , G_stressm_2 , G_stressm_3 ,  G_stressm_4,     &
                              G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4,     &
                              G_strength,                                                 &
                              G_cdn_ocn   , G_aiu       , G_uocn     ,  G_vocn     ,      &
                              G_waterxU   , G_wateryU   , G_forcexU  , G_forceyU   ,      &
                              G_umassdti  , G_fmU       ,                                 &
                              G_Tbu       , G_uvel     , G_vvel)


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
       !$omp target data map(to: ee, ne, se, nw, sw, sse, skipUcell, skipTcell,&
       !$omp                 strength, dxT_1d, dyT_1d, DminTarea_1d,           &
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
       str1(1:nStored_TU_glob)=c0
       str2(1:nStored_TU_glob)=c0
       str3(1:nStored_TU_glob)=c0
       str4(1:nStored_TU_glob)=c0
       str5(1:nStored_TU_glob)=c0
       str6(1:nStored_TU_glob)=c0
       str7(1:nStored_TU_glob)=c0
       str8(1:nStored_TU_glob)=c0

    ltripole = (trim(ns_boundary_type) == 'tripole')
       do ksub = 1,ndte        ! subcycling
          call stress_1d (ee, ne, se, 1, nActive_TU_glob,                                &
                          uvel, vvel, dxT_1d, dyT_1d, skipTcell, strength,           &
                          DminTarea_1d, HTE_1d, HTN_1d, HTEm1_1d, HTNm1_1d,          &
                          stressp_1,  stressp_2,  stressp_3,  stressp_4,             &
                          stressm_1,  stressm_2,  stressm_3,  stressm_4,             &
                          stress12_1, stress12_2, stress12_3, stress12_4,            &
                          str1, str2, str3, str4, str5, str6, str7, str8)

          call stepu_1d  (1, nActive_TU_glob, cdn_ocn, aiu, uocn, vocn,                     &
                          waterxU, wateryU, forcexU, forceyU, umassdti, fmU, uarear_1d, &
                          uvel_init, vvel_init, uvel, vvel,                             &
                          str1, str2, str3, str4, str5, str6, str7, str8,               &
                          nw, sw, sse, skipUcell, Tbu, Cb, rhow, rheofactU)
          if (ltripole) call evp1d_tripole_seam()
          call evp1d_halo_update()
       enddo
       if (ltripole) call evp1d_tripole_stress()
       ! This can be skipped if diagnostics of strintx and strinty is not needed
       ! They will either both be calculated or not.
       call calc_diag_1d(1   , nActive_TU_glob, &
                         uarear_1d, skipUcell, &
                         str1     , str2     , &
                         str3     , str4     , &
                         str5     , str6     , &
                         str7     , str8     , &
                         nw       , sw       , &
                         sse      , rheofactU, &
                         strintxU, strintyU)

       call ice_timer_stop(timer_evp1dcore)

#ifdef _OPENMP_TARGET
       !$omp end target data
#endif
       ! Map results back to 2d
       call convert_1d_2d_dyn(                                                                          &
                              G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
                              G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
                              G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
                              G_strength  ,                                           &
                              G_cdn_ocn   , G_aiu       , G_uocn      , G_vocn      , &
                              G_waterxU   , G_wateryU   , G_forcexU   , G_forceyU   , &
                              G_umassdti  , G_fmU       , G_strintxU  , G_strintyU  , &
                              G_Tbu       , G_uvel      , G_vvel      , G_taubxU    , &
                              G_taubyU)

    endif ! master_task

    call scatter_dyn(L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
                     L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
                     L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4, &
                     L_strintxU  , L_strintyU  ,  L_uvel     , L_vvel      , &
                     L_taubxU    , L_taubyU    ,                             &
                     G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
                     G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
                     G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
                     G_strintxU  , G_strintyU  , G_uvel      , G_vvel      , &
                     G_taubxU    , G_taubyU)
    ! calculate number of active points. allocate if initial or if array size should increase
    ! call calc_nActiveTU(iceTmask_log,nActive_TU_glob, iceUmask)
    ! if (nActiveold ==0) then ! first
    !     call evp_1d_alloc(nActive_TU_glob, nActive_TU_glob,nx,ny)
    !     nactiveold=nActive_TU_glob+buf1d ! allocate
    !     call init_unionTU(nx, ny, iceTmask_log,iceUmask)
    ! else if (nactiveold < nActive_TU_glob) then
    !     write(nu_diag,*) 'Warning nActive_TU_glob is bigger than old allocation. Need to re allocate'
    !     call evp_1d_dealloc() ! only deallocate if not first time step
    !     call evp_1d_alloc(nActive_TU_glob, nActive_TU_glob,nx,ny)
    !     nactiveold=nActive_TU_glob+buf1d ! allocate
    !     call init_unionTU(nx, ny, iceTmask_log,iceUmask)
    ! endif
    ! call cp_2dto1d(nActive_TU_glob)
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

  end subroutine dyn_evp1d_finalize

!=============================================================================

  subroutine evp1d_alloc_static_nActive
    implicit none

    integer(kind=int_kind) :: ierr
    character(len=*), parameter :: subname = '(evp1d_alloc_static_nActive)'

    allocate(skipTcell(1:nActive_TU_glob),  &
             skipUcell(1:nActive_TU_glob),  &
             iwidx(1:nx,1:ny),  &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

    allocate(indxTi(1:nActive_TU_glob), &
             indxTj(1:nActive_TU_glob), &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

    allocate(ee(1:nActive_TU_glob) , &
             ne(1:nActive_TU_glob) , &
             se(1:nActive_TU_glob) , &
             nw(1:nActive_TU_glob) , &
             sw(1:nActive_TU_glob) , &
             sse(1:nActive_TU_glob), &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

    allocate( HTE_1d    (1:nActive_TU_glob), &
              HTN_1d    (1:nActive_TU_glob), &
              HTEm1_1d  (1:nActive_TU_glob), &
              HTNm1_1d  (1:nActive_TU_glob), &
              dxT_1d    (1:nActive_TU_glob), &
              dyT_1d    (1:nActive_TU_glob), &
              DminTarea_1d(1:nActive_TU_glob), &
              strength  (1:nActive_TU_glob), &
              stressp_1 (1:nActive_TU_glob), &
              stressp_2 (1:nActive_TU_glob), &
              stressp_3 (1:nActive_TU_glob), &
              stressp_4 (1:nActive_TU_glob), &
              stressm_1 (1:nActive_TU_glob), &
              stressm_2 (1:nActive_TU_glob), &
              stressm_3 (1:nActive_TU_glob), &
              stressm_4 (1:nActive_TU_glob), &
              stress12_1(1:nActive_TU_glob), &
              stress12_2(1:nActive_TU_glob), &
              stress12_3(1:nActive_TU_glob), &
              stress12_4(1:nActive_TU_glob), &
              stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

    allocate(cdn_ocn  (1:nActive_TU_glob), aiu      (1:nActive_TU_glob), rheofactU(1:nActive_TU_glob), &
             uocn     (1:nActive_TU_glob), vocn     (1:nActive_TU_glob), &
             waterxU  (1:nActive_TU_glob), wateryU  (1:nActive_TU_glob), &
             forcexU  (1:nActive_TU_glob), forceyU  (1:nActive_TU_glob), &
             umassdti (1:nActive_TU_glob), fmU      (1:nActive_TU_glob), &
             uarear_1d(1:nActive_TU_glob),                   &
             strintxU (1:nActive_TU_glob), strintyU (1:nActive_TU_glob), &
             Tbu      (1:nActive_TU_glob), Cb       (1:nActive_TU_glob), &
             uvel_init(1:nActive_TU_glob), vvel_init(1:nActive_TU_glob), &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

  end subroutine evp1d_alloc_static_nActive

!=============================================================================

  subroutine evp1d_alloc_static_nStored
    implicit none

    integer(kind=int_kind) :: ierr
    character(len=*), parameter :: subname = '(evp1d_alloc_static_nStored)'

    ! nStored_TU_glob: the last slot is the dead slot.  Directions that leave the
    ! grid point at it, so the kernels can read a neighbour unconditionally
    ! and get zero -- the same convention the compressed 2-D layout uses.
    allocate(str1(1:nStored_TU_glob)   , str2(1:nStored_TU_glob), str3(1:nStored_TU_glob), &
             str4(1:nStored_TU_glob)   , str5(1:nStored_TU_glob), str6(1:nStored_TU_glob), &
             str7(1:nStored_TU_glob)   , str8(1:nStored_TU_glob),                   &
             indxTij(1:nStored_TU_glob), uvel(1:nStored_TU_glob), vvel(1:nStored_TU_glob), &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

  end subroutine evp1d_alloc_static_nStored

!=============================================================================

  subroutine evp1d_alloc_static_halo()

    implicit none
    integer(kind=int_kind) :: ierr
    character(len=*), parameter :: subname = '(evp1d_alloc_static_halo)'

    ! One entry per ghost cell per edge, for the cyclic pass and again for
    ! the zero_gradient/linear_extrap pass.
    allocate(halo_bc_sgn(4*nghost*(nx+ny)), &
             halo_bc_dst(4*nghost*(nx+ny)), halo_bc_s1(4*nghost*(nx+ny)), &
             halo_bc_s2(4*nghost*(nx+ny)), halo_bc_w (4*nghost*(nx+ny)), &
             stat=ierr)

    if (ierr/=0) then
       call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
    endif

  end subroutine evp1d_alloc_static_halo

!=============================================================================

  subroutine calc_nActiveTU(Tmask, Umask)

    ! Calculate number of active points with a given mask.

    implicit none
    logical(kind=log_kind), intent(in) :: Tmask(:,:)
    logical(kind=log_kind), optional, intent(in) :: Umask(:,:)
    integer(kind=int_kind)              :: i,j
    character(len=*), parameter :: subname = '(calc_nActivceTU)'

    ! Same range as calc_2d_indices_init, which fills the index arrays this
    ! counts: ilo..ihi+1, one cell into the ghosts, where ihi is nx-nghost.
    ! The two must agree or nActive_TU_glob exceeds what gets filled and the tail of
    ! indxti/indxtj is left at zero.
    nActive_TU_glob=0
    if (present(Umask)) then
       do i=1+nghost,nx-nghost+1
       do j=1+nghost,ny-nghost+1
          if ((Tmask(i,j)) .or. (Umask(i,j))) then
             nActive_TU_glob=nActive_TU_glob+1
          endif
       enddo
       enddo
    else
       do i=1+nghost,nx-nghost+1
       do j=1+nghost,ny-nghost+1
          if (Tmask(i,j)) then
             nActive_TU_glob=nActive_TU_glob+1
          endif
       enddo
       enddo
    endif

  end subroutine calc_nActiveTU

!=============================================================================

  subroutine set_skipMe(iceTmask, iceUmask)

    implicit none

    logical(kind=log_kind), intent(in) :: iceTmask(:,:), iceUmask(:,:)
    integer(kind=int_kind)              :: is, i, j
    character(len=*), parameter :: subname = '(set_skipMe)'

    skipUcell=.false.
    skipTcell=.false.
    do is=1, nActive_TU_glob
      i = indxti(is)
      j = indxtj(is)
      if (.not. (iceTmask(i,j))) skipTcell(is)=.true.
      if (.not. (iceUmask(i,j))) skipUcell(is)=.true.
      ! No U point outside ilo..ihi, jlo..jhi, which is where the 2-D solver
      ! steps U.  ihi is nx-nghost, so every ghost column past it is skipped,
      ! not just the outermost one -- the same one-ghost assumption as above.
      if (i > nx - nghost) skipUcell(is)=.true.
      if (j > ny - nghost) skipUcell(is)=.true.
    enddo

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

  subroutine gather_static(G_uarear, G_dxT, G_dyT, G_DminTarea, G_Tmask)

     ! In standalone  distrb_info is an integer. Not needed anyway
     use ice_communicate, only : master_task
     use ice_gather_scatter, only : gather_global
     use ice_domain, only : distrb_info
     use ice_grid, only: dyT, dxT, uarear, tmask
     use ice_dyn_shared, only: DminTarea
     implicit none

     ! nx, ny
     real(kind=dbl_kind)   , dimension(:,:), intent(out) :: G_uarear, G_dxT, G_dyT, &
                                                            G_DminTarea
     logical(kind=log_kind), dimension(:,:), intent(out) :: G_Tmask

     character(len=*), parameter :: subname = '(gather_static)'

     G_uarear = c0
     G_dyT = c0
     G_dxT = c0
     G_DminTarea = c0
     G_tmask = .false.

     ! copy from distributed I_* to G_*
     call gather_global(G_uarear, uarear, master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_dxT   , dxT   , master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_dyT   , dyT   , master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_DminTarea, DminTarea, master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_Tmask , Tmask , master_task, distrb_info, grid_ext=.true.)

  end subroutine gather_static

!=============================================================================

  subroutine gather_dyn(L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
                        L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
                        L_stress12_1, L_stress12_2, L_stress12_3,L_stress12_4 , &
                        L_strength  ,                                           &
                        L_cdn_ocn   , L_aiu       , L_uocn      , L_vocn      , &
                        L_waterxU   , L_wateryU   , L_forcexU   , L_forceyU   , &
                        L_umassdti  , L_fmU       ,                             &
                        L_Tbu       , L_uvel      , L_vvel      ,               &
                        L_icetmask  , L_iceUmask  ,                             &
                        L_strintxU  , L_strintyU  , L_taubxU    , L_taubyU    , &
                        G_strintxU  , G_strintyU  , G_taubxU    , G_taubyU    , &
                        G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
                        G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
                        G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
                        G_strength,                                             &
                        G_cdn_ocn   , G_aiu       , G_uocn      , G_vocn      , &
                        G_waterxU   , G_wateryU   , G_forcexU   , G_forceyU   , &
                        G_umassdti  , G_fmU       ,                             &
                        G_Tbu       , G_uvel      , G_vvel      ,               &
                        G_iceTmask,  G_iceUmask)

     use ice_communicate, only : master_task
     use ice_gather_scatter, only : gather_global
     use ice_domain, only : distrb_info
     implicit none

     ! nx_block, ny_block, max_blocks
     real(kind=dbl_kind)   , dimension(:,:,:), intent(in)  ::   &
        L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
        L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
        L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4, &
        L_strength  ,                                           &
        L_cdn_ocn   , L_aiu       , L_uocn      , L_vocn      , &
        L_waterxU   , L_wateryU   , L_forcexU   , L_forceyU   , &
        L_umassdti  , L_fmU       ,                             &
        L_Tbu       , L_uvel      , L_vvel
     logical(kind=log_kind), dimension(:,:,:), intent(in)  ::   &
        L_iceUmask  , L_iceTmask
     ! the four diagnostics, gathered so the round trip leaves alone any cell
     ! this solver does not own -- see the note in convert_1d_2d_dyn
     real(kind=dbl_kind)   , dimension(:,:,:), intent(in)  ::   &
        L_strintxU  , L_strintyU  , L_taubxU    , L_taubyU
     real(kind=dbl_kind)   , dimension(:,:), intent(inout) ::   &
        G_strintxU  , G_strintyU  , G_taubxU    , G_taubyU

     ! nx, ny
     real(kind=dbl_kind)   , dimension(:,:), intent(out) ::     &
        G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
        G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
        G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
        G_strength,                                             &
        G_cdn_ocn   , G_aiu       , G_uocn      , G_vocn      , &
        G_waterxU   , G_wateryU   , G_forcexU   , G_forceyU   , &
        G_umassdti  , G_fmU       ,                             &
        G_Tbu       , G_uvel      , G_vvel
     logical(kind=log_kind), dimension(:,:), intent(out) ::     &
        G_iceUmask  , G_iceTmask

     character(len=*), parameter :: subname = '(gather_dyn)'

     ! copy from distributed I_* to G_*
     call gather_global(G_stressp_1 ,     L_stressp_1,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stressp_2 ,     L_stressp_2,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stressp_3 ,     L_stressp_3,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stressp_4 ,     L_stressp_4,     master_task, distrb_info,c0, grid_ext=.true.)

     call gather_global(G_stressm_1 ,     L_stressm_1,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stressm_2 ,     L_stressm_2,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stressm_3 ,     L_stressm_3,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stressm_4 ,     L_stressm_4,     master_task, distrb_info,c0, grid_ext=.true.)

     call gather_global(G_stress12_1,     L_stress12_1,    master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stress12_2,     L_stress12_2,    master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stress12_3,     L_stress12_3,    master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_stress12_4,     L_stress12_4,    master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_strength  ,     L_strength  ,    master_task, distrb_info,c0, grid_ext=.true.)

     call gather_global(G_cdn_ocn   ,     L_cdn_ocn   ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_aiu       ,     L_aiu       ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_uocn      ,     L_uocn      ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_vocn      ,     L_vocn      ,     master_task, distrb_info, grid_ext=.true.)

     call gather_global(G_waterxU   ,     L_waterxU   ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_wateryU   ,     L_wateryU   ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_forcexU   ,     L_forcexU   ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_forceyU   ,     L_forceyU   ,     master_task, distrb_info, grid_ext=.true.)

     call gather_global(G_umassdti  ,     L_umassdti  ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_fmU       ,     L_fmU       ,     master_task, distrb_info, grid_ext=.true.)

     call gather_global(G_Tbu       ,     L_Tbu       ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_uvel      ,     L_uvel      ,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_vvel      ,     L_vvel      ,     master_task, distrb_info,c0, grid_ext=.true.)
     ! Symmetric with the scatter in scatter_dyn, so the round trip is an
     ! identity at every cell this solver does not write.
     call gather_global(G_strintxU  ,     L_strintxU  ,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_strintyU  ,     L_strintyU  ,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_taubxU    ,     L_taubxU    ,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_taubyU    ,     L_taubyU    ,     master_task, distrb_info,c0, grid_ext=.true.)
     call gather_global(G_iceTmask  ,     L_iceTmask  ,     master_task, distrb_info, grid_ext=.true.)
     call gather_global(G_iceUmask  ,     L_iceUmask  ,     master_task, distrb_info, grid_ext=.true.)

  end subroutine gather_dyn

!=============================================================================

  subroutine scatter_dyn(L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
                         L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
                         L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4, &
                         L_strintxU  , L_strintyU  , L_uvel      , L_vvel      , &
                         L_taubxU    , L_taubyU    ,                             &
                         G_stressp_1 , G_stressp_2 , G_stressp_3 ,  G_stressp_4, &
                         G_stressm_1 , G_stressm_2 , G_stressm_3 ,  G_stressm_4, &
                         G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
                         G_strintxU  , G_strintyU  , G_uvel      , G_vvel      , &
                         G_taubxU    , G_taubyU )

     use ice_communicate, only : master_task
     use ice_gather_scatter, only : scatter_global
     use ice_domain, only : distrb_info
     implicit none

     ! nx_block, ny_block, max_blocks
     real(kind=dbl_kind), dimension(:,:,:), intent(out) :: &
        L_stressp_1 , L_stressp_2 , L_stressp_3 , L_stressp_4 , &
        L_stressm_1 , L_stressm_2 , L_stressm_3 , L_stressm_4 , &
        L_stress12_1, L_stress12_2, L_stress12_3, L_stress12_4, &
        L_strintxU  , L_strintyU  , L_uvel      , L_vvel      , &
        L_taubxU    , L_taubyU

     ! nx, ny
     real(kind=dbl_kind), dimension(:,:), intent(in) ::         &
        G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
        G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
        G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
        G_strintxU  , G_strintyU  , G_uvel      , G_vvel      , &
        G_taubxU    , G_taubyU

     character(len=*), parameter :: subname = '(scatter_dyn)'

     call scatter_global(L_stressp_1,  G_stressp_1,  master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stressp_2,  G_stressp_2,  master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stressp_3,  G_stressp_3,  master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stressp_4,  G_stressp_4,  master_task, distrb_info, grid_ext=.true.)

     call scatter_global(L_stressm_1,  G_stressm_1,  master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stressm_2,  G_stressm_2,  master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stressm_3,  G_stressm_3,  master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stressm_4,  G_stressm_4,  master_task, distrb_info, grid_ext=.true.)

     call scatter_global(L_stress12_1, G_stress12_1, master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stress12_2, G_stress12_2, master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stress12_3, G_stress12_3, master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_stress12_4, G_stress12_4, master_task, distrb_info, grid_ext=.true.)

     call scatter_global(L_strintxU  , G_strintxU  , master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_strintyU  , G_strintyU  , master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_uvel      , G_uvel      , master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_vvel      , G_vvel      , master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_taubxU    , G_taubxU    , master_task, distrb_info, grid_ext=.true.)
     call scatter_global(L_taubyU    , G_taubyU    , master_task, distrb_info, grid_ext=.true.)

  end subroutine scatter_dyn

!=============================================================================

  subroutine convert_2d_1d_init(G_HTE, G_HTN, G_uarear,  G_dxT, G_dyT, G_DminTarea)

     implicit none

     real (kind=dbl_kind), dimension(:, :), intent(in) :: G_HTE, G_HTN, G_uarear, G_dxT, G_dyT, &
                                                          G_DminTarea

     ! local variables

     integer(kind=int_kind) :: is, lo, up, j, i, ij, islot

     character(len=*), parameter :: subname = '(convert_2d_1d_init)'

     ! Hand out the slots.  The layout is the one the kernels rely on: the
     ! active T and U points first, 1..nActive_TU_glob, in the order calc_2d_indices_init
     ! found them, then the cells that are only ever read.  stress_1d and
     ! stepu_1d both run over 1..nActive_TU_glob, so nothing downstream has to know where
     ! the halo tail begins beyond nActive_TU_glob itself.
     do is = 1, nActive_TU_glob
        ij = indxti(is) + (indxtj(is) - 1) * nx
        ijslot(ij)  = is
        indxTij(is) = ij
     end do

     ! The read-only cells take the slots after them, in ascending linear
     ! index -- the same order the sorted union used to produce.
     islot = nActive_TU_glob
     do ij = 1, nx * ny
        if (ijslot(ij) == -1) then
           islot = islot + 1
           ijslot(ij)     = islot
           indxTij(islot) = ij
        endif
     end do
     if (islot /= nStored_TU_glob - 1) then
        call abort_ice(subname//' ERROR: slot count does not match nStored_TU_glob', &
                       file=__FILE__, line=__LINE__)
     endif

     ! Neighbour slots.  A direction that leaves the grid takes the dead slot
     ! at nStored_TU_glob, which holds zero and is never written.  That only happens
     ! for nw/sw/sse at i == nx or j == ny, where the cell has no U point to
     ! step -- skipUcell is already .true. there, matching the 2-D solver,
     ! which computes stress over jlo..jhi+1 but steps U only over jlo..jhi.
     ! ee/ne/se step towards i-1 and j-1 and active T cells start at
     ! 1+nghost, so for nghost >= 1 they always land inside the grid.
     do is = 1, nActive_TU_glob
        i = indxti(is)
        j = indxtj(is)
        ee (is) = slot_of(i - 1, j    , nStored_TU_glob)
        ne (is) = slot_of(i - 1, j - 1, nStored_TU_glob)
        se (is) = slot_of(i    , j - 1, nStored_TU_glob)
        nw (is) = slot_of(i + 1, j    , nStored_TU_glob)
        sw (is) = slot_of(i + 1, j + 1, nStored_TU_glob)
        sse(is) = slot_of(i    , j + 1, nStored_TU_glob)
     end do

     !tar      i$OMP PARALLEL PRIVATE(is, lo, up, j, i)
     ! write 1D arrays from 2D arrays (target points)
     !tar      call domp_get_domain(1, nActive_TU_glob, lo, up)
     lo=1
     up=nActive_TU_glob
     do is = 1, nActive_TU_glob
        ! get 2D indices
        i = indxti(is)
        j = indxtj(is)
        ! map
        uarear_1d(is)  = G_uarear(i, j)
        dxT_1d(is)     = G_dxT(i, j)
        dyT_1d(is)     = G_dyT(i, j)
        DminTarea_1d(is) = G_DminTarea(i, j)
        HTE_1d(is)     = G_HTE(i, j)
        HTN_1d(is)     = G_HTN(i, j)
        HTEm1_1d(is)   = G_HTE(i - 1, j)
        HTNm1_1d(is)   = G_HTN(i, j - 1)
     end do

  end subroutine convert_2d_1d_init

!=============================================================================

  subroutine convert_2d_1d_dyn(                                                                    &
                               G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
                               G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
                               G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
                               G_strength  , G_cdn_ocn   , G_aiu       , G_uocn      , &
                               G_vocn      , G_waterxU   , G_wateryU   , G_forcexU   , &
                               G_forceyU   , G_umassdti  , G_fmU       ,  G_Tbu      , &
                               G_uvel      , G_vvel      )

     implicit none


     ! nx, ny
     real(kind=dbl_kind), dimension(:, :), intent(in) :: &
        G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4, &
        G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4, &
        G_stress12_1, G_stress12_2, G_stress12_3,G_stress12_4, &
        G_strength  , G_cdn_ocn   , G_aiu       , G_uocn     , &
        G_vocn      , G_waterxU   , G_wateryU   , G_forcexU  , &
        G_forceyU   , G_umassdti  , G_fmU       , G_Tbu      , &
        G_uvel      , G_vvel

     integer(kind=int_kind) ::  lo, up, is, i, j
     character(len=*), parameter :: subname = '(convert_2d_1d_dyn)'

     lo=1
     up=nActive_TU_glob
     do is = 1, nActive_TU_glob
        ! get 2D indices
        i = indxti(is)
        j = indxtj(is)
        ! map
        stressp_1(is)  = G_stressp_1(i, j)
        stressp_2(is)  = G_stressp_2(i, j)
        stressp_3(is)  = G_stressp_3(i, j)
        stressp_4(is)  = G_stressp_4(i, j)
        stressm_1(is)  = G_stressm_1(i, j)
        stressm_2(is)  = G_stressm_2(i, j)
        stressm_3(is)  = G_stressm_3(i, j)
        stressm_4(is)  = G_stressm_4(i, j)
        stress12_1(is) = G_stress12_1(i, j)
        stress12_2(is) = G_stress12_2(i, j)
        stress12_3(is) = G_stress12_3(i, j)
        stress12_4(is) = G_stress12_4(i, j)
        strength(is)   = G_strength(i,j)
        cdn_ocn(is)    = G_cdn_ocn(i, j)
        aiu(is)        = G_aiu(i, j)
        ! The rule dyn_prep2 uses: 1, or 0 where the U cell carries less ice
        ! than rheo_area_min.  #1133 added this factor to the 2-D stepu but
        ! not to this kernel, so the two solvers disagreed on any grid with
        ! marginal ice -- gx3 shows it, a solid-ice box grid does not.
        if (free_drift_approx .and. G_aiu(i, j) <= rheo_area_min) then
           rheofactU(is) = c0
        else
           rheofactU(is) = c1
        endif
        uocn(is)       = G_uocn(i, j)
        vocn(is)       = G_vocn(i, j)
        waterxU(is)    = G_waterxU(i, j)
        wateryU(is)    = G_wateryU(i, j)
        forcexU(is)    = G_forcexU(i, j)
        forceyU(is)    = G_forceyU(i, j)
        umassdti(is)   = G_umassdti(i, j)
        fmU(is)        = G_fmU(i, j)
        strintxU(is)   = C0
        strintyU(is)   = C0
        Tbu(is)        = G_Tbu(i, j)
        Cb(is)         = c0
        uvel(is)       = G_uvel(i,j)
        vvel(is)       = G_vvel(i,j)
        uvel_init(is)  = G_uvel(i,j)
        vvel_init(is)  = G_vvel(i,j)
     end do

     ! Halos can potentially have values of u and v
     do is=nActive_TU_glob+1,nStored_TU_glob-1
        j = int((indxTij(is) - 1) / (nx)) + 1
        i =      indxTij(is) - (j - 1) * nx
        uvel(is)=G_uvel(i,j)
        vvel(is)=G_vvel(i,j)
     end do

  end subroutine convert_2d_1d_dyn

!=============================================================================

  subroutine convert_1d_2d_dyn(                                                                                &
                               G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
                               G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
                               G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
                               G_strength,                                             &
                               G_cdn_ocn   , G_aiu       , G_uocn      , G_vocn      , &
                               G_waterxU   , G_wateryU   , G_forcexU   , G_forceyU   , &
                               G_umassdti  , G_fmU       , G_strintxU  , G_strintyU  , &
                               G_Tbu       , G_uvel      , G_vvel      , G_taubxU    , &
                               G_taubyU)

     implicit none

     ! nx, ny
     real(kind=dbl_kind), dimension(:, :), intent(inout) ::     &
        G_stressp_1 , G_stressp_2 , G_stressp_3 , G_stressp_4 , &
        G_stressm_1 , G_stressm_2 , G_stressm_3 , G_stressm_4 , &
        G_stress12_1, G_stress12_2, G_stress12_3, G_stress12_4, &
        G_strength,                                             &
        G_cdn_ocn   , G_aiu       , G_uocn      , G_vocn      , &
        G_waterxU   , G_wateryU   , G_forcexU   , G_forceyU   , &
        G_umassdti  , G_fmU       , G_strintxU  , G_strintyU  , &
        G_Tbu       , G_uvel      , G_vvel      , G_taubxU    , &
        G_taubyU

     integer(kind=int_kind) ::  lo, up, is, i, j
     character(len=*), parameter :: subname = '(convert_1d_2d_dyn)'

     G_stressp_1  = c0
     G_stressp_2  = c0
     G_stressp_3  = c0
     G_stressp_4  = c0
     G_stressm_1  = c0
     G_stressm_2  = c0
     G_stressm_3  = c0
     G_stressm_4  = c0
     G_stress12_1 = c0
     G_stress12_2 = c0
     G_stress12_3 = c0
     G_stress12_4 = c0
     G_strength   = c0
     G_cdn_ocn    = c0
     G_aiu        = c0
     G_uocn       = c0
     G_vocn       = c0
     G_waterxU    = c0
     G_wateryU    = c0
     G_forcexU    = c0
     G_forceyU    = c0
     G_umassdti   = c0
     G_fmU        = c0
     ! G_strintxU, G_strintyU, G_taubxU and G_taubyU are deliberately NOT
     ! zeroed, for the same reason as G_uvel and G_vvel below: gather_dyn has
     ! filled them with the block field, and scatter_dyn writes them straight
     ! back.  The 2-D solver's unpack_from_slots writes only the cells it has
     ! a slot for and leaves the rest, so zeroing here made evp1d clear
     ! diagnostics at cells it does not own where the 2-D solver does not.
     G_Tbu        = c0
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
     up=nActive_TU_glob
     do is = lo, up
        ! get 2D indices
        i = indxti(is)
        j = indxtj(is)
        ! map to 2d
        G_stressp_1 (i,j) = stressp_1(is)
        G_stressp_2 (i,j) = stressp_2(is)
        G_stressp_3 (i,j) = stressp_3(is)
        G_stressp_4 (i,j) = stressp_4(is)
        G_stressm_1 (i,j) = stressm_1(is)
        G_stressm_2 (i,j) = stressm_2(is)
        G_stressm_3 (i,j) = stressm_3(is)
        G_stressm_4 (i,j) = stressm_4(is)
        G_stress12_1(i,j) = stress12_1(is)
        G_stress12_2(i,j) = stress12_2(is)
        G_stress12_3(i,j) = stress12_3(is)
        G_stress12_4(i,j) = stress12_4(is)
        G_strintxU(i,j)   = strintxU(is)
        G_strintyU(i,j)   = strintyU (is)
        G_taubxU(i,j)     = -uvel(is)*Cb(is)
        G_taubyU(i,j)     = -vvel(is)*Cb(is)
        G_uvel(i,j)       = uvel(is)
        G_vvel(i,j)       = vvel(is)
     end do

     do is=nActive_TU_glob+1,nStored_TU_glob-1
        j = int((indxTij(is) - 1) / (nx)) + 1
        i = indxTij(is) - (j - 1) * nx
        G_uvel(i,j)       = uvel(is)
        G_vvel(i,j)       = vvel(is)
     end do

  end subroutine convert_1d_2d_dyn

  subroutine calc_nStoredTU
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

 
     ! local variables

     integer(kind=int_kind) :: is, i, j, ij, ierr

     character(len=*), parameter :: subname = '(calc_nStoredTU)'

     allocate(ijslot(1:nx*ny), stat=ierr)
     if (ierr/=0) then
        call abort_ice(subname//' ERROR: allocating', file=__FILE__, line=__LINE__)
     endif

     ijslot(:) = 0

     do is = 1, nActive_TU_glob
        i = indxti(is)
        j = indxtj(is)
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

     nStored_TU_glob = 0
     do ij = 1, nx * ny
        if (ijslot(ij) /= 0) nStored_TU_glob = nStored_TU_glob + 1
     end do

     ! ⚠ Plus the dead slot.  nStored_TU_glob counts what is STORED, which
     ! includes the one slot at the end that owns no cell and absorbs every
     ! direction leaving the domain -- the same convention nStored_TU uses
     ! for the 2-D layout.  Folding it in here removes the +1 that was
     ! otherwise written out at every allocation and every slot_of fallback.
     nStored_TU_glob = nStored_TU_glob + 1

  end subroutine calc_nStoredTU

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
     integer(kind=int_kind) :: is
     character(len=*), parameter :: subname = '(numainit)'

     !$omp parallel do schedule(runtime) private(is)
     do is = lo,up
        skipTcell(is)=.false.
        skipUcell(is)=.false.
        ee(is)=0
        ne(is)=0
        se(is)=0
        nw(is)=0
        sw(is)=0
        sse(is)=0
        aiu(is)=c0
        rheofactU(is)=c1
        Cb(is)=c0
        cdn_ocn(is)=c0
        dxT_1d(is)=c0
        dyT_1d(is)=c0
        DminTarea_1d(is)=c0
        fmU(is)=c0
        forcexU(is)=c0
        forceyU(is)=c0
        HTE_1d(is)=c0
        HTEm1_1d(is)=c0
        HTN_1d(is)=c0
        HTNm1_1d(is)=c0
        strength(is)= c0
        stress12_1(is)=c0
        stress12_2(is)=c0
        stress12_3(is)=c0
        stress12_4(is)=c0
        stressm_1(is)=c0
        stressm_2(is)=c0
        stressm_3(is)=c0
        stressm_4(is)=c0
        stressp_1(is)=c0
        stressp_2(is)=c0
        stressp_3(is)=c0
        stressp_4(is)=c0
        strintxU(is)= c0
        strintyU(is)= c0
        Tbu(is)=c0
        uarear_1d(is)=c0
        umassdti(is)=c0
        uocn(is)=c0
        uvel_init(is)=c0
        uvel(is)=c0
        vocn(is)=c0
        vvel_init(is)=c0
        vvel(is)=c0
        waterxU(is)=c0
        wateryU(is)=c0
     enddo
     !$omp end parallel do
     !$omp parallel do schedule(runtime) private(is)
     do is = lo,uu
        uvel(is)=c0
        vvel(is)=c0
        str1(is)=c0
        str2(is)=c0
        str3(is)=c0
        str4(is)=c0
        str5(is)=c0
        str6(is)=c0
        str7(is)=c0
        str8(is)=c0
     enddo
     !$omp end parallel do

  end subroutine numainit

!=======================================================================

  subroutine halo_sweep(domark)
     ! Walk every ghost point the boundary conditions fill, in the order the
     ! 2-D halo fills them.
     !
     ! One sweep, two uses.  With domark it marks the cells the halo update
     ! reads and writes, so calc_nStoredTU gives them slots; without it, it builds
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

     integer(kind=int_kind) :: i, j, g, ilo, ihi, jlo, jhi, w, isrc

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
           call halo_point(g      , j, ihi-nghost+g, j, ihi-nghost+g, j, 0, domark)
           call halo_point(ihi + g, j, ilo+g-1     , j, ilo+g-1     , j, 0, domark)
        end do
        end do
     endif

     ! ⚠ AFTER the east/west cyclic pass, not before.  This is the
     ! north/south pass and the ordering rule above applies to it: the
     ! north/south pass must read ghosts the east/west pass has already
     ! written, or the corners come out wrong.  Running the fold first
     ! let the cyclic pass overwrite it at the two ghost columns with an
     ! unfolded source and sgn +1 -- measured on gbox80 tripoleT as the
     ! four corner cells (ghost columns x seam row and north ghost row).
     if (trim(ns_boundary_type) == 'tripole') then
        call tripole_sweep(domark)

        ! u(i, jhi+g) = -u(nx-i, jhi-g).  The source row is jhi-g, NOT the seam,
        ! so these are independent of the seam symmetry and may run either side
        ! of it.
        !
        ! ⚠ i runs to nx-1, not nx-2.  The earlier bound assumed the cyclic
        ! pass would fill what the fold skipped, but cyclic writes only the
        ! GHOST columns 1 and nx -- column nx-1 is PHYSICAL, so nothing wrote
        ! (nx-1, jhi+g) at all, and (1, jhi+g) then inherited that stale value
        ! cyclically.  Measured as the entire residual: with ndte=1 the only
        ! cells differing from the 2-D solver were uvel and vvel at exactly
        ! those two.
        !
        ! At i = nx-1 the mirror column nx-i is 1, a ghost.  Under a cyclic
        ! east-west boundary that ghost IS column nx-1, so the source is taken
        ! there directly rather than relying on the cyclic pass having run.
        ! Without cyclic there is no such image and the column is left to the
        ! old behaviour.
        do g = 1, nghost
        ! ⚠ to nx, not nx-1.  The east ghost column is a fold destination
        ! like any other; stopping short left it to the cyclic pass, which
        ! has no fold and no sign.
        do i = 1, nx
           isrc = nx - i
           if (isrc < 1 + nghost) then
              if (trim(ew_boundary_type) /= 'cyclic') cycle
              isrc = isrc + (nx - 2*nghost)
           endif
           call halo_point(i, jhi + g, isrc, jhi - g, isrc, jhi - g, &
                           0, domark, -c1)
        end do
        end do
     endif

     if (trim(ns_boundary_type) == 'cyclic') then
        do i = 1, nx
        do g = 1, nghost
           call halo_point(i, g      , i, jhi-nghost+g, i, jhi-nghost+g, 0, domark)
           call halo_point(i, jhi + g, i, jlo+g-1     , i, jlo+g-1     , 0, domark)
        end do
        end do
     endif

     ! then the extrapolating types, east/west before north/south
     if (trim(ew_boundary_type) == 'zero_gradient' .or. &
         trim(ew_boundary_type) == 'linear_extrap') then
        do j = 1, ny
        do g = 1, nghost
           w = 0
           if (trim(ew_boundary_type) == 'linear_extrap') w = nghost-g+1
           call halo_point(g      , j, ilo, j, ilo+1, j, w, domark)   ! west
           w = 0
           if (trim(ew_boundary_type) == 'linear_extrap') w = g
           call halo_point(ihi + g, j, ihi, j, ihi-1, j, w, domark)   ! east
        end do
        end do
     endif

     if (trim(ns_boundary_type) == 'zero_gradient' .or. &
         trim(ns_boundary_type) == 'linear_extrap') then
        do i = 1, nx
        do g = 1, nghost
           w = 0
           if (trim(ns_boundary_type) == 'linear_extrap') w = nghost-g+1
           call halo_point(i, g      , i, jlo, i, jlo+1, w, domark)   ! south
           w = 0
           if (trim(ns_boundary_type) == 'linear_extrap') w = g
           call halo_point(i, jhi + g, i, jhi, i, jhi-1, w, domark)   ! north
        end do
        end do
     endif

  end subroutine halo_sweep

!=======================================================================

  subroutine halo_point(id, jd, i1, j1, i2, j2, w, domark, sgn)
     ! Mark, or append, one halo point.
     !
     ! Marking is what lets the append succeed: these cells are the halo
     ! update's own inputs and outputs, and a ghost whose interior neighbour
     ! is land is reached by no kernel, so without marking it has no slot and
     ! the point would be dropped.  The 2-D halo fills it regardless, which is
     ! where evp1d used to differ.

     implicit none

     integer(kind=int_kind), intent(in) :: id, jd, i1, j1, i2, j2, w
     logical(kind=log_kind), intent(in) :: domark

     ! local variables

     integer(kind=int_kind) :: sd, s1, s2

     real   (kind=dbl_kind), intent(in), optional :: sgn
     character(len=*), parameter :: subname = '(halo_point)'

     if (domark) then
        call mark_cell(id, jd)
        call mark_cell(i1, j1)
        if (w /= 0) call mark_cell(i2, j2)
        return
     endif

     sd = ijslot(id + (jd - 1) * nx)
     s1 = ijslot(i1 + (j1 - 1) * nx)
     if (sd == 0 .or. s1 == 0) return

     s2 = ijslot(i2 + (j2 - 1) * nx)
     if (s2 == 0) then
        ! w == 0 does not read s2, so a missing one is no obstacle
        if (w /= 0) return
        s2 = s1
     endif

     n_halo_bc = n_halo_bc + 1
     halo_bc_dst(n_halo_bc) = sd
     halo_bc_s1 (n_halo_bc) = s1
     halo_bc_s2 (n_halo_bc) = s2
     halo_bc_w  (n_halo_bc) = real(w, kind=dbl_kind)
     ! +1 everywhere unless the caller says otherwise.  A tripole fold
     ! negates both velocity components across the seam and the entry form
     ! had no way to say so.  x*1.0 is exact, so this changes nothing on any
     ! configuration that never passes sgn.
     halo_bc_sgn(n_halo_bc) = c1
     if (present(sgn)) halo_bc_sgn(n_halo_bc) = sgn

  end subroutine halo_point

!=======================================================================

  subroutine evp1d_halo_update()

     implicit none
     integer(kind=int_kind) :: is

     character(len=*), parameter :: subname = '(evp1d_halo_update)'

! One sequential pass over the list built in halo_sweep.  Not parallelised
! and not to be: north/south destinations are east/west sources at the
! corners, so the entries are order-dependent.  See halo_sweep.
     do is = 1, n_halo_bc
        uvel(halo_bc_dst(is)) = halo_bc_sgn(is) *                           &
                               (uvel(halo_bc_s1(is)) + halo_bc_w(is) *      &
                               (uvel(halo_bc_s1(is)) - uvel(halo_bc_s2(is))))
        vvel(halo_bc_dst(is)) = halo_bc_sgn(is) *                           &
                               (vvel(halo_bc_s1(is)) + halo_bc_w(is) *      &
                               (vvel(halo_bc_s1(is)) - vvel(halo_bc_s2(is))))
     end do

  end subroutine evp1d_halo_update

  subroutine evp1d_tripole_seam()

     ! Symmetry ALONG the tripole seam row, applied before the halo list.
     !
     ! ⚠ This cannot be a halo_bc entry, for two independent reasons:
     !
     !   the arithmetic.  An entry computes sgn*(s1 + w*(s1-s2)); the seam
     !   needs p5*(x1-x2), and solving the coefficients for that gives
     !   sgn = 0.  Two independent coefficients would be required.
     !
     !   the ordering.  Both ends of a pair are written FROM THE ORIGINAL
     !   PAIR.  The halo list is applied sequentially and deliberately so --
     !   its north/south entries read ghosts the east/west ones just wrote --
     !   so the second write of a pair would read the first's result.
     !
     ! Visiting each pair once and writing both ends from saved locals is
     ! what makes it correct.
     !
     ! Mapping measured against the 2-D halo on a one-block tripole grid,
     ! not transcribed from ice_boundary: for a pair (i, nx-i) on row jhi,
     !     xavg = p5*(u(i) - u(nx-i)),  u(i) = xavg,  u(nx-i) = -xavg

     implicit none

     integer(kind=int_kind) :: k, ka, kb
     real   (kind=dbl_kind) :: ua, ub, va, vb

     character(len=*), parameter :: subname = '(evp1d_tripole_seam)'

     do k = 1, n_tri_pair
        ka = tri_a(k)
        kb = tri_b(k)
        ua = uvel(ka)
        ub = uvel(kb)
        va = vvel(ka)
        vb = vvel(kb)
        uvel(ka) =  p5*(ua - ub)
        uvel(kb) = -p5*(ua - ub)
        vvel(ka) =  p5*(va - vb)
        vvel(kb) = -p5*(va - vb)
     end do

     ! The two columns with no distinct partner -- i = nx/2 and i = nx-1 --
     ! are negated outright.  Measured, not assumed: both come back as
     ! -(original) from the 2-D halo.
     do k = 1, n_tri_neg
        uvel(tri_neg(k)) = -uvel(tri_neg(k))
        vvel(tri_neg(k)) = -vvel(tri_neg(k))
     end do

  end subroutine evp1d_tripole_seam

  subroutine evp1d_tripole_stress()

     ! The tripole fold for the twelve stresses, applied ONCE PER TIMESTEP
     ! after the subcycle -- where evp's epilogue applies
     ! ice_HaloUpdate_stress to the block/slot copies.
     !
     ! ⚠ Why evp1d needs its own.  That epilogue updates stressp_1_slot and
     ! friends, but evp1d subcycles on its OWN stressp_1, seeded once by
     ! gather_stresses and never refreshed.  Without this the seam never
     ! reaches the working stresses and the error enters once per timestep --
     ! which is exactly what the measurement showed: bit-for-bit within a
     ! timestep at any ndte, diverging only across them.
     !
     ! ⚠ Each component is filled from its DIAGONAL PARTNER, 1<->3 and 2<->4,
     ! matching the pairing in evp's twelve calls.  Sign is +1: these go in as
     ! field_type_scalar, unlike the velocities.
     !
     ! Re-gathering from the slot copies would be the other way to do it, and
     ! would cost twelve collectives per timestep -- the ones 2677687 removed.

     implicit none
     integer(kind=int_kind) :: k, d, sr

     character(len=*), parameter :: subname = '(evp1d_tripole_stress)'

     ! ⚠ A one-way copy into the GHOST, not a swap: the handler fills
     ! array1's ghost cells from array2's interior and leaves the interior
     ! alone.  Both partners are read before either is written, because the
     ! two directed updates of a pair would otherwise see each other.
     do k = 1, n_tri_str
        d = tst_d(k)
        sr = tst_s(k)
        call swap_pair(stressp_1 , stressp_3 , d, sr)
        call swap_pair(stressp_2 , stressp_4 , d, sr)
        call swap_pair(stressm_1 , stressm_3 , d, sr)
        call swap_pair(stressm_2 , stressm_4 , d, sr)
        call swap_pair(stress12_1, stress12_3, d, sr)
        call swap_pair(stress12_2, stress12_4, d, sr)
     end do

  contains

     subroutine swap_pair(f1, f3, d, sr)
        real   (kind=dbl_kind), dimension(:), intent(inout) :: f1, f3
        integer(kind=int_kind), intent(in) :: d, sr
        real   (kind=dbl_kind) :: v1, v3
        v1 = f1(sr)
        v3 = f3(sr)
        f1(d) = v3
        f3(d) = v1
     end subroutine swap_pair

  end subroutine evp1d_tripole_stress

  subroutine tripole_sweep(domark)

     ! Mark, or build, the seam-row symmetry lists.  Same two-pass shape as
     ! halo_sweep: the marking pass is what gives these cells slots at all.

     implicit none

     logical(kind=log_kind), intent(in) :: domark

     integer(kind=int_kind) :: i, ip, g, jhi, sa, sb, ierr

     character(len=*), parameter :: subname = '(tripole_sweep)'

     jhi = ny - nghost

     if (domark) then
        do i = 2, nx/2 - 1
           call mark_cell(i     , jhi)
           call mark_cell(nx - i, jhi)
        end do
        call mark_cell(nx/2  , jhi)
        call mark_cell(nx - 1, jhi)
        ! and the ghost rows the stress fold writes, plus their sources.
        ! Columns only within calc_2d_indices_init's range: marking outside
        ! it just grows the read-only tail, which the stresses cannot index.
        do g = 1, nghost
        do i = 1 + nghost, nx - nghost + 1
           call mark_cell(i, jhi + g)
           call mark_cell(nx + 1 - i, jhi + 1 - g)
        end do
        end do
        return
     endif

     if (.not. allocated(tri_a)) then
        allocate(tri_a(nx), tri_b(nx), tri_neg(2), &
                 tst_d(nx*nghost), tst_s(nx*nghost), stat=ierr)
        if (ierr /= 0) call abort_ice(subname//' ERROR: allocating', &
                                      file=__FILE__, line=__LINE__)
     endif

     ! Ghost-row fold for the stresses.  field_loc_center and
     ! field_type_scalar, unlike the velocities' NEcorner/vector, so the
     ! offsets are 0,0 instead of 1,1 and the sign is +1 instead of -1:
     !
     !     stressX_1(i, jhi+g) = stressX_3(nx+1-i, jhi+1-g)
     !
     ! A centre-located field has no degenerate seam row -- that symmetry
     ! exists only for NEcorner and Nface -- so there is no stage 1 here.
     n_tri_str = 0
     do g = 1, nghost
     do i = 1, nx
        ! Measured against ice_HaloUpdate_stress on a one-block tripole grid:
        !     dst (i, jhi+g) <- partner (nx+1-i, jhi+1-g)
        ! verified i=8 -> 75, i=15 -> 68, i=78 -> 5, source row jhi.
        !
        ! ⚠ At the ends nx+1-i leaves the physical range and the handler uses
        ! the CYCLIC IMAGE, not the ghost: measured i=1 -> 2, where nx+1-i is
        ! 82, the east ghost, whose image is column 2.  Reading the ghost
        ! instead was the whole of the residual the derived mapping left.
        ip = nx + 1 - i
        if (ip > nx - nghost) ip = ip - (nx - 2*nghost)
        if (ip < 1 + nghost)  ip = ip + (nx - 2*nghost)
        if (ip < 1 .or. ip > nx) cycle
        sa = ijslot(i  + (jhi + g - 1) * nx)
        sb = ijslot(ip + (jhi + 1 - g - 1) * nx)
        ! ⚠ Both ends must be ACTIVE slots.  The stress arrays are
        ! dimension(nActive_TU_glob); slots above that are the read-only neighbour
        ! tail and have no stress storage at all.  Cells outside
        ! calc_2d_indices_init's range -- the ghost COLUMNS i=1 and i=nx --
        ! land there, and indexing them overran f1 by exactly nx*ny-nActive_TU_glob.
        if (sa == 0 .or. sb == 0) cycle
        if (sa > nActive_TU_glob .or. sb > nActive_TU_glob) cycle
        n_tri_str = n_tri_str + 1
        tst_d(n_tri_str) = sa
        tst_s(n_tri_str) = sb
     end do
     end do

     n_tri_pair = 0
     do i = 2, nx/2 - 1
        ip = nx - i
        sa = ijslot(i  + (jhi - 1) * nx)
        sb = ijslot(ip + (jhi - 1) * nx)
        if (sa == 0 .or. sb == 0) cycle
        n_tri_pair = n_tri_pair + 1
        tri_a(n_tri_pair) = sa
        tri_b(n_tri_pair) = sb
     end do

     ! The two columns with no distinct partner: nx/2 and nx-1.
     n_tri_neg = 0
     do ip = 1, 2
        i = nx/2
        if (ip == 2) i = nx - 1
        sa = ijslot(i + (jhi - 1) * nx)
        if (sa == 0) cycle
        n_tri_neg = n_tri_neg + 1
        tri_neg(n_tri_neg) = sa
     end do

  end subroutine tripole_sweep

     subroutine swap_pair(f1, f3, d, sr)
        real   (kind=dbl_kind), dimension(:), intent(inout) :: f1, f3
        integer(kind=int_kind), intent(in) :: d, sr
        real   (kind=dbl_kind) :: v1, v3
        v1 = f1(sr)
        v3 = f3(sr)
        f1(d) = v3
        f3(d) = v1
     end subroutine swap_pair

end module ice_dyn_evp1d

