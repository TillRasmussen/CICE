!=======================================================================
!
! Elastic-viscous-plastic sea ice dynamics model
! Computes ice velocity and deformation
!
! See:
!
! Hunke, E. C., and J. K. Dukowicz (1997). An elastic-viscous-plastic model
! for sea ice dynamics. J. Phys. Oceanogr., 27, 1849-1867.
!
! Hunke, E. C. (2001).  Viscous-Plastic Sea Ice Dynamics with the EVP Model:
! Linearization Issues. J. Comput. Phys., 170, 18-38.
!
! Hunke, E. C., and J. K. Dukowicz (2002).  The Elastic-Viscous-Plastic
! Sea Ice Dynamics Model in General Orthogonal Curvilinear Coordinates
! on a Sphere - Incorporation of Metric Terms. Mon. Weather Rev.,
! 130, 1848-1865.
!
! Hunke, E. C., and J. K. Dukowicz (2003).  The sea ice momentum
! equation in the free drift regime.  Los Alamos Tech. Rep. LA-UR-03-2219.
!
! Hibler, W. D. (1979). A dynamic thermodynamic sea ice model. J. Phys.
! Oceanogr., 9, 817-846.
!
! Bouillon, S., T. Fichefet, V. Legat and G. Madec (2013).  The
! elastic-viscous-plastic method revisited.  Ocean Model., 71, 2-12.
!
! author: Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb (LANL)
! 2004: Block structure added by William Lipscomb
! 2005: Removed boundary calls for stress arrays (WHL)
! 2006: Streamlined for efficiency by Elizabeth Hunke
!       Converted to free source form (F90)

      module ice_dyn_evp

      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: field_loc_center, field_loc_NEcorner, &
          field_loc_Nface, field_loc_Eface, &
          field_type_scalar, field_type_vector
      use ice_constants, only: c0, p027, p055, p111, p166, &
          p222, p25, p333, p5, c1
      use ice_dyn_shared, only: stepu, stepuv_CD, stepu_C, stepv_C, &
          calc_seabed_stress_slot, calc_strint_slot, &
          str1_slot, str2_slot, str3_slot, str4_slot, &
          str5_slot, str6_slot, str7_slot, str8_slot, &
          dyn_prep1, dyn_prep2, dyn_finish, &
          ndte, yield_curve, ecci, denom1, arlx1i, &
          uvel_init, vvel_init, uvelE_init, vvelE_init, uvelN_init, vvelN_init, &
          dynOff, dynCnt, dynMap, Nbr_Local, &
        NbrW_Task , NbrSW_Task, NbrS_Task , NbrE_Task , &
        NbrE_Task, NbrNE_Task, NbrN_Task, &
        NbrNE_Task, NbrN_Task , NbrNW_Task, NbrSE_Task, &
          uvel_slot, vvel_slot, aiU_slot, fmU_slot, TbU_slot, &
          rheofactU_slot, &
          umassdti_slot, waterxU_slot, wateryU_slot, &
          forcexU_slot, forceyU_slot, &
          uocnU_slot, vocnU_slot, cdn_ocnU_slot, uarear_slot, tarear_slot, dynIJ, &
          pack_to_slots, unpack_from_slots, &
          nbr_W, nbr_S, nbr_SW, nbr_E, nbr_N, nbr_NE, &
          iceTmask_slot, iceUmask_slot, iceNmask_slot, iceEmask_slot, &
          stressp_1_slot , stressp_2_slot , stressp_3_slot , stressp_4_slot , &
          stressm_1_slot , stressm_2_slot , stressm_3_slot , stressm_4_slot , &
          stress12_1_slot, stress12_2_slot, stress12_3_slot, stress12_4_slot, &
          dxhy_slot, dyhx_slot, cxp_slot, cyp_slot, cxm_slot, cym_slot, &
          DminTarea_slot, strength_slot, dxT_slot, dyT_slot, &
          seabed_stress_factor_LKD, seabed_stress_factor_prob, seabed_stress_method, &
          seabed_stress, Ktens, revp
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_ice_strength, icepack_query_parameters

      implicit none
      private

      ! Grid ratios for the BCs, built once in slots.  The OWN metric comes
      ! from the slot's cell, but the NEIGHBOUR metric is read from the block
      ! array: dxN at the E neighbour is often ordinary land, which carries a
      ! real metric yet has no slot, so Nbr*_Task would resolve it to the
      ! shared dead slot and read c1 instead.
      real (kind=dbl_kind), dimension (:), allocatable :: &
         ratiodxN_slot , ratiodyE_slot , &
         ratiodxNr_slot, ratiodyEr_slot

      ! rheofact compressed.  dyn_prep2 fills these in the same loop that
      ! computes the block forms, so there is no second traversal.
      real (kind=dbl_kind), dimension (:), allocatable :: &
         rheofactE_slot, rheofactN_slot

      real (kind=dbl_kind), allocatable :: &
         divergU  (:)   , & ! div array on U points, differentiate from divu
         tensionU (:)   , & ! tension array on U points
         shearU   (:,:) , & ! shear array on U points
         deltaU   (:)      ! delta array on U points

      ! B grid
      real (kind=dbl_kind), allocatable :: &
         waterxU  (:,:)   , & ! for ocean stress calculation, x (m/s)
         wateryU  (:,:)   , & ! for ocean stress calculation, y (m/s)
         forcexU  (:,:)   , & ! work array: combined atm stress and ocn tilt, x
         forceyU  (:,:)       ! work array: combined atm stress and ocn tilt, y

      ! the i and j dimensions are flattened: element (i,j,iblk) of the block
      ! layout lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), allocatable, target :: &
         umass    (:,:)   , & ! total mass of ice and snow (u grid)
         umassdti (:,:)       ! mass of U-cell/dte (kg/m^2 s)

      ! the i and j dimensions are flattened: element (i,j,iblk) of the block
      ! layout lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), allocatable, target :: &
         tmass    (:,:)       ! total mass of ice and snow (kg/m^2)

      public :: evp, init_evp

!=======================================================================

      contains

!=======================================================================
! Elastic-viscous-plastic dynamics driver
!
      subroutine init_evp
      use ice_blocks, only: get_block, nx_block, ny_block, block
      use ice_domain_size, only: max_blocks
      use ice_domain, only: nblocks, blocks_ice
      use ice_grid, only: grid_ice, dxN, dyE
      use ice_grid, only: dynIJ, dynBlk
      use ice_calendar, only: dt_dyn
      use ice_dyn_shared, only: init_dyn_shared, evp_algorithm, &
          dynTot, dynOff, dynCnt, dynMap, dxN_slot, dyE_slot
      use ice_grid, only: pack_to_slots
      use ice_dyn_evp1d, only: dyn_evp1d_init

!allocate c and cd grid var. Follow structucre of eap
      integer (int_kind) :: ierr

      character(len=*), parameter :: subname = '(init_evp)'

     type (block) :: &
         this_block   ! block information for current block

      integer (kind=int_kind) :: &
         i, j, iblk            , & ! block index
         ijb                   , & ! 1D index of (i,j) in the flattened block layout
         is                    , & ! slot index
         ilo,ihi,jlo,jhi           ! beginning and end of physical domain


      call init_dyn_shared(dt_dyn)

      ! Guarded on grid_ice as well as evp_algorithm, to match the run-time
      ! call below.  dyn_evp1d_run sits inside "if (grid_ice == B)", but this
      ! init did not, so a C or CD run with shared_mem_1d used to build the
      ! whole 1d machinery -- the static gathers, convert_2d_1d_init, the halo
      ! pairing and some forty persistent arrays of length navel on the master
      ! task -- for a solver that is then never called.  The 1d solver has
      ! never been implemented for C/CD; until it is, do not pay for it.
      if (evp_algorithm == "shared_mem_1d" .and. grid_ice == "B") then
         call dyn_evp1d_init
      endif

      allocate( tmass    (nx_block*ny_block,max_blocks), & ! total mass of ice and snow (kg/m^2)
                waterxU  (nx_block*ny_block,max_blocks), & ! for ocean stress calculation, x (m/s)
                wateryU  (nx_block*ny_block,max_blocks), & ! for ocean stress calculation, y (m/s)
                forcexU  (nx_block*ny_block,max_blocks), & ! work array: combined atm stress and ocn tilt, x
                forceyU  (nx_block*ny_block,max_blocks), & ! work array: combined atm stress and ocn tilt, y
                umass    (nx_block*ny_block,max_blocks), & ! total mass of ice and snow (u grid)
                umassdti (nx_block*ny_block,max_blocks), & ! mass of U-cell/dte (kg/m^2 s)
                stat=ierr)
      if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory B-Grid evp')


      if (grid_ice == 'CD' .or. grid_ice == 'C') then

         allocate( divergU  (dynTot), &
                   tensionU (dynTot), &
                   shearU   (nx_block*ny_block,max_blocks), &
                   deltaU   (dynTot), &
                   stat=ierr)
         if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory U evp')


         allocate( ratiodxN_slot (dynTot), ratiodyE_slot (dynTot), &
                   ratiodxNr_slot(dynTot), ratiodyEr_slot(dynTot), &
                   rheofactE_slot(dynTot), rheofactN_slot(dynTot), &
                   stat=ierr)
         if (ierr/=0) call abort_ice(subname//' ERROR: Out of memory ratio slot')
         ratiodxN_slot  = c0
         ratiodyE_slot  = c0
         ratiodxNr_slot = c0
         ratiodyEr_slot = c0
         ! c1, not c0: since #1133 dyn_prep2 fills these only when
         ! free_drift_approx is set, so a c0 default would zero the internal
         ! stress everywhere with the flag off.  This initialisation is now
         ! the ONLY thing that makes the flag-off case correct -- the block
         ! copies that used to carry the same c1 are gone.
         rheofactE_slot = c1
         rheofactN_slot = c1

         ! Constant for the run, so built once here, in a slot loop.
         !
         ! The neighbour metric now comes from the slots.  It could not
         ! before: dxN at the E neighbour is often ordinary LAND, which
         ! carries a real metric, and until the neighbour ring joined
         ! mask_any that cell had no slot -- Nbr*_Task resolved it to the
         ! shared dead slot and read c1 instead of the cell's value.
         !
         ! The interior test stays.  It is what keeps the loop off the
         ! padding, where the metric is genuinely zero and the division
         ! would trap; it is not about neighbours.
         do is = 1, dynTot
            ijb = dynIJ(is)
            if (ijb == 0) cycle                  ! dead slot
            iblk = dynBlk(is)
            this_block = get_block(blocks_ice(iblk),iblk)
            j = (ijb-1)/nx_block + 1
            i = ijb - (j-1)*nx_block
            if (i < this_block%ilo .or. i > this_block%ihi .or. &
                j < this_block%jlo .or. j > this_block%jhi) cycle
            ratiodxN_slot (is) = - dxN_slot(NbrE_Task(is)) / dxN_slot(is)
            ratiodyE_slot (is) = - dyE_slot(NbrN_Task(is)) / dyE_slot(is)
            ratiodxNr_slot(is) =   c1 / ratiodxN_slot(is)
            ratiodyEr_slot(is) =   c1 / ratiodyE_slot(is)
         enddo

      endif

      end subroutine init_evp

!=======================================================================
#ifdef CICE_IN_NEMO
! Wind stress is set during this routine from the values supplied
! via NEMO (unless calc_strair is true).  These values are supplied
! rotated on u grid and multiplied by aice.  strairxT = 0 in this
! case so operations in dyn_prep1 are pointless but carried out to
! minimise code changes.
#endif
!
! author: Elizabeth C. Hunke, LANL

      subroutine evp (dt)

      use ice_arrays_column, only: Cdn_ocn
      use ice_boundary, only: ice_halo, ice_HaloMask, ice_HaloUpdate, &
          ice_HaloUpdate2DR8_compressed, ice_HaloUpdate3DR8_compressed, &
          ice_HaloDestroy, ice_HaloUpdate_stress, &
          ice_HaloUpdate_stress_compressed
      use ice_blocks, only: block, get_block, nx_block, ny_block, ns_boundary_type
      use ice_domain, only: nblocks, blocks_ice, halo_info, maskhalo_dyn
      use ice_domain_size, only: max_blocks
      use ice_flux, only: rdg_conv, rdg_shear, strairxT, strairyT, &
          strairxU, strairyU, uocn, vocn, ss_tltx, ss_tlty, fmU, &
          strtltxU, strtltyU, strocnxU, strocnyU, strintxU, strintyU, taubxU, taubyU, &
          strax, stray, &
          hwater, hwater_slot, &
          strairxN, strairyN, fmN, &
          strtltxN, strtltyN, strocnxN, strocnyN, strintxN, strintyN, taubxN, taubyN, &
          strairxE, strairyE, fmE, &
          strtltxE, strtltyE, strocnxE, strocnyE, strintxE, strintyE, taubxE, taubyE, &
          stressp_1, stressp_2, stressp_3, stressp_4, &
          stressm_1, stressm_2, stressm_3, stressm_4, &
          stress12_1, stress12_2, stress12_3, stress12_4, &
          stresspT, stressmT, stress12T, &
          stresspU, stressmU, stress12U
      use ice_grid, only: tmask, umask, umaskCD, nmask, emask, uvm, epm, npm, &
          dxE, dxN, dxT, dxU, dyE, dyN, dyT, dyU, &
          tarear, uarear, earear, narear, grid_average_X2Y, uarea, &
          grid_ice, grid_atm_dynu, grid_atm_dynv, grid_ocn_dynu, grid_ocn_dynv
      use ice_restoring, only: ice_restoring_halo, restore_ice
      use ice_state, only: aice, aiU, vice, vsno, aice_slot, vice_slot, vsno_slot, uvel, vvel, uvelN, vvelN, &
          aicen_slot, vicen_slot, &
          uvelE, vvelE, divu, shear, vort, &
          aice_init, aice0, aicen, vicen, strength
      use ice_timers, only: timer_dynamics, timer_bound, &
          ice_timer_start, ice_timer_stop, timer_evp
      use ice_dyn_shared, only: evp_algorithm, stack_fields, unstack_fields, &
          DminTarea, visc_method, deformations, deformationsC_T, deformationsCD_T, &
          strain_rates_U, dxhy, dyhx, cxp, cyp, cxm, cym, &
          iceTmask, iceUmask, iceEmask, iceNmask, &
          dyn_haloUpdate, fld2, fld3, fld4, fld2c
      ! the compressed layout, and the C/CD slot arrays the subcycle runs on
      use ice_dyn_shared, only: dynMap, dynOff, dynCnt, dynTot, Nbr_Local, &
                               NbrW_Task , NbrSW_Task, NbrS_Task , NbrE_Task , &
                               NbrNE_Task, NbrN_Task , NbrNW_Task, NbrSE_Task
      ! strint* and taub* are diagnostics: history is their only consumer,
      ! nothing in the dynamics reads them.  So when no stream wants them this
      ! step, the work that builds them is skipped outright.  Gated coarsely
      ! -- one flag for every strint field, one for every taub field --
      ! because a run either asks for these or it does not, and a coarse test
      ! cannot miss a consumer the way a per-face one can.
      use ice_history_shared, only: hist_field_needed, &
          n_strintx , n_strinty , n_strintxN, n_strintyN, &
          n_strintxE, n_strintyE, n_siforceintstrx, n_siforceintstry, &
          n_taubx   , n_tauby   , n_taubxN  , n_taubyN  , &
          n_taubxE  , n_taubyE
      use ice_dyn_shared, only: uvelE_slot, vvelE_slot, uvelN_slot, vvelN_slot, &
          shearU_slot, stresspT_slot, stressmT_slot, stress12T_slot, &
          stress12U_slot, stresspU_slot, stressmU_slot, &
          uvel_slot, vvel_slot, &
          strintxU_slot, strintyU_slot, taubxU_slot, taubyU_slot, &
          strintxE_slot, strintyE_slot, strintxN_slot, strintyN_slot, &
          taubxE_slot, taubyE_slot, taubxN_slot, taubyN_slot, &
          aiE_slot, cdn_ocnE_slot, uocnE_slot, vocnE_slot, &
          ss_tltxU_slot, ss_tltyU_slot, &
          umass_slot, emass_slot, nmass_slot, &
          fcor_slot, fcorE_slot, fcorN_slot, &
          cdn_ocnU_slot, cdn_ocnn_slot, cdn_ocne_slot, &
          taubxU_slot, taubyU_slot, taubxE_slot, taubyE_slot, taubxN_slot, taubyN_slot, &
          ss_tltxE_slot, ss_tltyE_slot, &
          ss_tltxN_slot, ss_tltyN_slot, &
          waterxE_slot, wateryE_slot, forcexE_slot, forceyE_slot, &
          emassdti_slot, fmE_slot, TbE_slot, &
          aiN_slot, cdn_ocnN_slot, uocnN_slot, vocnN_slot, &
          waterxN_slot, wateryN_slot, forcexN_slot, forceyN_slot, &
          nmassdti_slot, fmN_slot, TbN_slot
      use ice_dyn_shared, only: uarea_slot, dxN_slot, dyE_slot, &
          earear_slot, narear_slot, &
          dxU_slot, dyT_slot, dxT_slot, dyU_slot, &
          dxE_slot, dyN_slot, epm_slot, npm_slot, uvm_slot, &
          tarea_slot, hm_slot, earea_slot, narea_slot, &
          uareaavg_slot, wsumN_slot, wsumE_slot, &
          wsumNW_slot, wsumSE_slot, wsumNE_slot, &
          zetax2T_slot, etax2T_slot, zetax2U_slot, etax2U_slot, &
          strength_slot, strengthU_slot
      use ice_dyn_shared, only: grid_average_X2YS_compressed, &
          grid_average_X2YA_compressed, fld3c, fld4c
      use ice_grid, only: pack_to_slots, unpack_from_slots
      use ice_boundary, only: ice_HaloUpdate2DR8_compressed, &
          ice_HaloUpdate3DR8_compressed
      use ice_domain, only: maskhalo_dyn
      use ice_dyn_evp1d, only: dyn_evp1d_run

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         ms               , & ! slot of this cell
         is               , & ! slot within this block
         ksub           , & ! subcycle step
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         i, j, ij           ! local indices



      logical (kind=log_kind) :: &
         calc_strair  ! calculate air/ice stress

      logical (kind=log_kind) :: &
         lstrint , & ! a history stream wants strintx/strinty this step
         ltaub       ! a history stream wants taubx/tauby this step

      integer (kind=int_kind), dimension (nx_block,ny_block,max_blocks) :: &
         halomask     ! generic halo mask

      type (ice_halo) :: &
         halo_info_mask !  ghost cell update info for masked halo

      type (block) :: &
         this_block   ! block information for current block

      character(len=*), parameter :: subname = '(evp)'

      call ice_timer_start(timer_dynamics) ! dynamics

      lstrint = hist_field_needed(n_strintx ) .or. hist_field_needed(n_strinty ) .or. &
                hist_field_needed(n_strintxN) .or. hist_field_needed(n_strintyN) .or. &
                hist_field_needed(n_strintxE) .or. hist_field_needed(n_strintyE) .or. &
                hist_field_needed(n_siforceintstrx) .or. &
                hist_field_needed(n_siforceintstry)
      ltaub   = hist_field_needed(n_taubx   ) .or. hist_field_needed(n_tauby   ) .or. &
                hist_field_needed(n_taubxN  ) .or. hist_field_needed(n_taubyN  ) .or. &
                hist_field_needed(n_taubxE  ) .or. hist_field_needed(n_taubyE  )


      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      if (grid_ice == 'CD' .or. grid_ice == 'C') then

         divergU  (:) = c0
         tensionU (:) = c0
         shearU   (:,:) = c0
         deltaU   (:) = c0

      endif

       ! This call is needed only if dt changes during runtime.
!      call set_evp_parameters (dt)

      !-----------------------------------------------------------------
      ! boundary updates
      ! commented out because the ghost cells are freshly
      ! updated after cleanup_itd
      !-----------------------------------------------------------------

!      call ice_timer_start(timer_bound)
!      call ice_HaloUpdate (aice,              halo_info &
!                           field_loc_center,  field_type_scalar)
!      call ice_HaloUpdate (vice,              halo_info &
!                           field_loc_center,  field_type_scalar)
!      call ice_HaloUpdate (vsno,              halo_info &
!                           field_loc_center,  field_type_scalar)
!      call ice_timer_stop(timer_bound)

      ! Compressed: one pass over dynTot, not nblocks x nx_block x ny_block.
      rdg_conv  = c0
      rdg_shear = c0
      divu      = c0
      shear     = c0
      vort      = c0

      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ilo,ihi,jlo,jhi,this_block) SCHEDULE(runtime)
      do iblk = 1, nblocks

         !-----------------------------------------------------------------
         ! preparation for dynamics
         !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call dyn_prep1 (nx_block,           ny_block,           &
                         ilo, ihi,           jlo, jhi,           &
                           aice_slot,          vice_slot,          &
                           vsno_slot,                              &
                           dynMap  (:,iblk),   tmask   (:,  iblk), &
                         dynOff(iblk),       dynCnt(iblk),       &
                         tmass   (:,iblk),   iceTmask(:,iblk))

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (iceTmask,          halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_timer_stop(timer_bound)

      !-----------------------------------------------------------------
      ! convert fields from T to U grid
      !-----------------------------------------------------------------

      call stack_fields(tmass, aice_init, cdn_ocn, fld3)
      call ice_HaloUpdate (fld3,             halo_info, &
                           field_loc_center, field_type_scalar)
      call stack_fields(uocn, vocn, ss_tltx, ss_tlty, fld4)
      call ice_HaloUpdate (fld4,             halo_info, &
                           field_loc_center, field_type_vector)
      call unstack_fields(fld3, tmass, aice_init, cdn_ocn)
      call unstack_fields(fld4, uocn, vocn, ss_tltx, ss_tlty)

      call grid_average_X2Y('S', tmass    , 'T'          , umass   , 'U')
      call grid_average_X2Y('S', tmass    , 'T'          , umass_slot, 'U')
      call grid_average_X2Y('S', cdn_ocn  , 'T'          , cdn_ocnU_slot, 'U')
      call grid_average_X2Y('S', aice_init, 'T'          , aiU     , 'U')
      call grid_average_X2Y('S', ss_tltx  , grid_ocn_dynu, ss_tltxU_slot, 'U')
      call grid_average_X2Y('S', ss_tlty  , grid_ocn_dynv, ss_tltyU_slot, 'U')

      ! ⚠ The compressed forms are produced HERE, by the average itself, not
      ! packed from the block result afterwards and not inside the
      ! grid_ice=="B" branch below.  dyn_finish reads them unconditionally --
      ! the C and CD paths call it for U as well as for N and E -- so filling
      ! them under "B" left them unwritten on the C path while the block
      ! originals were filled for every grid.
      !
      ! Resolution is by the RANK of the destination: a rank-1 actual selects
      ! grid_average_X2Y_base_toslot, a rank-2/3 one the block specifics.  A
      ! mistake here is a compile error, not a silent one.
      call grid_average_X2Y('S', aice_init, 'T'          , aiU_slot  , 'U')
      call grid_average_X2Y('S', uocn     , grid_ocn_dynu, uocnU_slot, 'U')
      call grid_average_X2Y('S', vocn     , grid_ocn_dynv, vocnU_slot, 'U')

      if (grid_ice == 'CD' .or. grid_ice == 'C') then
         call grid_average_X2Y('S', tmass    , 'T'          , emass_slot, 'E')
         call grid_average_X2Y('S', cdn_ocn  , 'T'          , cdn_ocne_slot, 'E')
         call grid_average_X2Y('S', ss_tltx  , grid_ocn_dynu, ss_tltxE_slot, 'E')
         call grid_average_X2Y('S', ss_tlty  , grid_ocn_dynv, ss_tltyE_slot, 'E')
         call grid_average_X2Y('S', tmass    , 'T'          , nmass_slot, 'N')
         call grid_average_X2Y('S', cdn_ocn  , 'T'          , cdn_ocnn_slot, 'N')
         call grid_average_X2Y('S', ss_tltx  , grid_ocn_dynu, ss_tltxN_slot, 'N')
         call grid_average_X2Y('S', ss_tlty  , grid_ocn_dynv, ss_tltyN_slot, 'N')

         ! Same here: produced compressed rather than packed from the block
         ! result further down.  The block forms of these are gone: nothing
         ! reads them any more.
         call grid_average_X2Y('S', aice_init, 'T'          , aiE_slot  , 'E')
         call grid_average_X2Y('S', uocn     , grid_ocn_dynu, uocnE_slot, 'E')
         call grid_average_X2Y('S', vocn     , grid_ocn_dynv, vocnE_slot, 'E')
         call grid_average_X2Y('S', aice_init, 'T'          , aiN_slot  , 'N')
         call grid_average_X2Y('S', uocn     , grid_ocn_dynu, uocnN_slot, 'N')
         call grid_average_X2Y('S', vocn     , grid_ocn_dynv, vocnN_slot, 'N')
      endif
      !----------------------------------------------------------------
      ! Set wind stress to values supplied via NEMO or other forcing
      ! Map T to U, N, E as needed
      ! This wind stress is rotated on u grid and multiplied by aice in coupler
      !----------------------------------------------------------------
      call icepack_query_parameters(calc_strair_out=calc_strair)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (.not. calc_strair) then
         call grid_average_X2Y('F', strax, grid_atm_dynu, strairxU, 'U')
         call grid_average_X2Y('F', stray, grid_atm_dynv, strairyU, 'U')
      else
         call ice_HaloUpdate (strairxT,         halo_info, &
                              field_loc_center, field_type_vector)
         call ice_HaloUpdate (strairyT,         halo_info, &
                              field_loc_center, field_type_vector)
         call grid_average_X2Y('F', strairxT, 'T', strairxU, 'U')
         call grid_average_X2Y('F', strairyT, 'T', strairyU, 'U')
      endif

      if (grid_ice == 'CD' .or. grid_ice == 'C') then
         if (.not. calc_strair) then
            call grid_average_X2Y('F', strax   , grid_atm_dynu, strairxN, 'N')
            call grid_average_X2Y('F', stray   , grid_atm_dynv, strairyN, 'N')
            call grid_average_X2Y('F', strax   , grid_atm_dynu, strairxE, 'E')
            call grid_average_X2Y('F', stray   , grid_atm_dynv, strairyE, 'E')
         else
            call grid_average_X2Y('F', strairxT, 'T'          , strairxN, 'N')
            call grid_average_X2Y('F', strairyT, 'T'          , strairyN, 'N')
            call grid_average_X2Y('F', strairxT, 'T'          , strairxE, 'E')
            call grid_average_X2Y('F', strairyT, 'T'          , strairyE, 'E')
         endif
      endif

      if (trim(grid_ice) == 'B') then
         !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,ij,i,j,is,ms) SCHEDULE(runtime)
         do iblk = 1, nblocks

            !-----------------------------------------------------------------
            ! more preparation for dynamics
            !-----------------------------------------------------------------

            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            call dyn_prep2 (nx_block,             ny_block,             &
                            ilo, ihi,             jlo, jhi,             &
                            aiU_slot,   umass_slot,   &
                            umassdti_slot, &
                            fcor_slot, &
                            umask     (:,iblk)  ,                       &
                            uocnU_slot, vocnU_slot, &
                            strairxU, &
                            strairyU, &
                            ss_tltxU_slot, &
                            ss_tltyU_slot, &
                            iceTmask  (:,iblk),   iceUmask  (:,iblk),   &
                            fmU_slot, dt                  , &
                            strtltxU, &
                            strtltyU, &
                            strintxU  (:,iblk), strintyU  (:,iblk), &
                            taubxU_slot  , taubyU_slot  , &
                            waterxU_slot, &
                            wateryU_slot, &
                            forcexU_slot, &
                            forceyU_slot, &
                            stressp_1 (:,iblk), stressp_2 (:,iblk), &
                            stressp_3 (:,iblk), stressp_4 (:,iblk), &
                            stressm_1 (:,iblk), stressm_2 (:,iblk), &
                            stressm_3 (:,iblk), stressm_4 (:,iblk), &
                            stress12_1(:,iblk), stress12_2(:,iblk), &
                            stress12_3(:,iblk), stress12_4(:,iblk), &
                            uvel_init,   vvel_init,   &
                            uvel(:,iblk), vvel(:,iblk), &
                           dynMap(:,iblk), dynOff(iblk), dynCnt(iblk), iceUmask_slot, &
                           iceTmask_slot, &
                            TbU_slot, &
                            rheofactU_slot)

            !-----------------------------------------------------------------
            ! ice strength
            !-----------------------------------------------------------------

            strength(:,:,iblk) = c0  ! initialize
            ! ⚠ iceTmask_slot, not skipPhys.  icellT is built over jlo..jhi+1 x
            ! ilo..ihi+1 -- the physical domain PLUS one north and east ghost
            ! row, deliberately, because stresses are needed there -- and it
            ! selects on iceTmask.  dyn_prep1 records exactly that set in
            ! iceTmask_slot, the same mask the stress kernel walks, so the slot loop
            ! reproduces icellT cell for cell with no dynMap lookup.
            do is = 1, dynCnt(iblk)
               if (.not. iceTmask_slot(dynOff(iblk)+is)) cycle
               ms = dynOff(iblk) + is
               ij = dynIJ(ms)
               j  = (ij-1)/nx_block + 1
               i  = ij - (j-1)*nx_block
               call icepack_ice_strength(aice     = aice    (i,j,  iblk), &
                                         vice     = vice    (i,j,  iblk), &
                                         aice0    = aice0   (i,j,  iblk), &
                                         aicen    = aicen_slot   (ms,:), &
                                         vicen    = vicen_slot   (ms,:), &
                                         strength = strength(i,j,  iblk))
            enddo  ! is

         enddo  ! iblk
         !$OMP END PARALLEL DO
      elseif (trim(grid_ice) == 'CD' .or. grid_ice == 'C') then
         !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,ij,i,j,is,ms) SCHEDULE(runtime)
         do iblk = 1, nblocks

            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            call dyn_prep2 (nx_block,             ny_block,             &
                            ilo, ihi,             jlo, jhi,             &
                            aiU_slot,   umass_slot,   &
                            umassdti_slot, &
                            fcor_slot, &
                            umaskCD   (:,iblk)  ,                       &
                            uocnU_slot, vocnU_slot, &
                            strairxU, &
                            strairyU, &
                            ss_tltxU_slot, &
                            ss_tltyU_slot, &
                            iceTmask  (:,iblk),   iceUmask  (:,iblk),   &
                            fmU_slot, dt,                   &
                            strtltxU, &
                            strtltyU, &
                            strintxU  (:,iblk), strintyU  (:,iblk), &
                            taubxU_slot  , taubyU_slot  , &
                            waterxU_slot, &
                            wateryU_slot, &
                            forcexU_slot, &
                            forceyU_slot, &
                            stressp_1 (:,iblk), stressp_2 (:,iblk), &
                            stressp_3 (:,iblk), stressp_4 (:,iblk), &
                            stressm_1 (:,iblk), stressm_2 (:,iblk), &
                            stressm_3 (:,iblk), stressm_4 (:,iblk), &
                            stress12_1(:,iblk), stress12_2(:,iblk), &
                            stress12_3(:,iblk), stress12_4(:,iblk), &
                            uvel_init,   vvel_init,   &
                            uvel(:,iblk), vvel(:,iblk), &
                           dynMap(:,iblk), dynOff(iblk), dynCnt(iblk), iceUmask_slot, &
                           iceTmask_slot, &
                            TbU_slot)

            !-----------------------------------------------------------------
            ! ice strength
            !-----------------------------------------------------------------

            strength(:,:,iblk) = c0  ! initialize
            ! ⚠ iceTmask_slot, not skipPhys.  icellT is built over jlo..jhi+1 x
            ! ilo..ihi+1 -- the physical domain PLUS one north and east ghost
            ! row, deliberately, because stresses are needed there -- and it
            ! selects on iceTmask.  dyn_prep1 records exactly that set in
            ! iceTmask_slot, the same mask the stress kernel walks, so the slot loop
            ! reproduces icellT cell for cell with no dynMap lookup.
            do is = 1, dynCnt(iblk)
               if (.not. iceTmask_slot(dynOff(iblk)+is)) cycle
               ms = dynOff(iblk) + is
               ij = dynIJ(ms)
               j  = (ij-1)/nx_block + 1
               i  = ij - (j-1)*nx_block
               call icepack_ice_strength(aice     = aice    (i,j,  iblk), &
                                         vice     = vice    (i,j,  iblk), &
                                         aice0    = aice0   (i,j,  iblk), &
                                         aicen    = aicen_slot   (ms,:), &
                                         vicen    = vicen_slot   (ms,:), &
                                         strength = strength(i,j,  iblk) )
            enddo  ! is


            !-----------------------------------------------------------------
            ! more preparation for dynamics on N grid
            !-----------------------------------------------------------------

            call dyn_prep2 (nx_block,             ny_block,             &
                            ilo, ihi,             jlo, jhi,             &
                            aiN_slot,   nmass_slot,   &
                            nmassdti_slot, &
                            fcorN_slot, &
                            nmask     (:,iblk)  ,                       &
                            uocnN_slot, vocnN_slot, &
                            strairxN, &
                            strairyN, &
                            ss_tltxN_slot, &
                            ss_tltyN_slot, &
                            iceTmask  (:,iblk),   iceNmask  (:,iblk),   &
                            fmN_slot, dt                  , &
                            strtltxN, &
                            strtltyN, &
                            strintxN  (:,:,iblk), strintyN  (:,:,iblk), &
                            taubxN_slot  , taubyN_slot  , &
                            waterxN_slot, &
                            wateryN_slot, &
                            forcexN_slot, &
                            forceyN_slot, &
                            stressp_1 (:,iblk), stressp_2 (:,iblk), &
                            stressp_3 (:,iblk), stressp_4 (:,iblk), &
                            stressm_1 (:,iblk), stressm_2 (:,iblk), &
                            stressm_3 (:,iblk), stressm_4 (:,iblk), &
                            stress12_1(:,iblk), stress12_2(:,iblk), &
                            stress12_3(:,iblk), stress12_4(:,iblk), &
                            uvelN_init,   vvelN_init,   &
                            uvelN     (:,:,iblk), vvelN     (:,:,iblk), &
                           dynMap(:,iblk), dynOff(iblk), dynCnt(iblk), iceNmask_slot, &
                           iceTmask_slot, &
                            TbN_slot, &
                            rheofactN_slot)

            !-----------------------------------------------------------------
            ! more preparation for dynamics on E grid
            !-----------------------------------------------------------------

            call dyn_prep2 (nx_block,             ny_block,             &
                            ilo, ihi,             jlo, jhi,             &
                            aiE_slot,   emass_slot,   &
                            emassdti_slot, &
                            fcorE_slot, &
                            emask     (:,iblk)  ,                       &
                            uocnE_slot, vocnE_slot, &
                            strairxE, &
                            strairyE, &
                            ss_tltxE_slot, &
                            ss_tltyE_slot, &
                            iceTmask  (:,iblk),   iceEmask  (:,iblk),   &
                            fmE_slot, dt                  , &
                            strtltxE, &
                            strtltyE, &
                            strintxE  (:,:,iblk), strintyE  (:,:,iblk), &
                            taubxE_slot  , taubyE_slot  , &
                            waterxE_slot, &
                            wateryE_slot, &
                            forcexE_slot, &
                            forceyE_slot, &
                            stressp_1 (:,iblk), stressp_2 (:,iblk), &
                            stressp_3 (:,iblk), stressp_4 (:,iblk), &
                            stressm_1 (:,iblk), stressm_2 (:,iblk), &
                            stressm_3 (:,iblk), stressm_4 (:,iblk), &
                            stress12_1(:,iblk), stress12_2(:,iblk), &
                            stress12_3(:,iblk), stress12_4(:,iblk), &
                            uvelE_init,   vvelE_init,   &
                            uvelE     (:,:,iblk), vvelE     (:,:,iblk), &
                           dynMap(:,iblk), dynOff(iblk), dynCnt(iblk), iceEmask_slot, &
                           iceTmask_slot, &
                            TbE_slot, &
                            rheofactE_slot)

            do i=1,nx_block
            do j=1,ny_block
               if (.not.iceUmask((j-1)*nx_block+i,iblk)) then
                  stresspU ((j-1)*nx_block+i,iblk) = c0
                  stressmU ((j-1)*nx_block+i,iblk) = c0
                  stress12U((j-1)*nx_block+i,iblk) = c0
               endif
               if (.not.iceTmask((j-1)*nx_block+i,iblk)) then
                  stresspT ((j-1)*nx_block+i,iblk) = c0
                  stressmT ((j-1)*nx_block+i,iblk) = c0
                  stress12T((j-1)*nx_block+i,iblk) = c0
               endif
            enddo
            enddo
         enddo  ! iblk
         !$OMP END PARALLEL DO

      endif ! grid_ice

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (grid_ice == 'CD' .or. grid_ice == 'C') then

         call ice_timer_start(timer_bound)
         call ice_HaloUpdate (uvelE,           halo_info, &
                              field_loc_Eface, field_type_vector)
         call ice_HaloUpdate (vvelN,           halo_info, &
                              field_loc_Nface, field_type_vector)
         call ice_timer_stop(timer_bound)

         if (grid_ice == 'C') then
            call grid_average_X2Y('A', uvelE, 'E', uvelN, 'N')
            call grid_average_X2Y('A', vvelN, 'N', vvelE, 'E')
            uvelN(:,:,:) = uvelN(:,:,:)*reshape(npm, shape(uvelN))
            vvelE(:,:,:) = vvelE(:,:,:)*reshape(epm, shape(vvelE))
         endif

         call ice_timer_start(timer_bound)
         call ice_HaloUpdate (uvelN,           halo_info, &
                              field_loc_Nface, field_type_vector)
         call ice_HaloUpdate (vvelE,           halo_info, &
                              field_loc_Eface, field_type_vector)
         call ice_timer_stop(timer_bound)

         call grid_average_X2Y('A', uvelE, 'E', uvel, 'U')
         call grid_average_X2Y('A', vvelN, 'N', vvel, 'U')
         uvel(:,:) = uvel(:,:)*reshape(uvm, shape(uvel))
         vvel(:,:) = vvel(:,:)*reshape(uvm, shape(vvel))
      endif

      ! strength_slot is filled HERE, ahead of the halo, because the halo is
      ! now compressed.  The packs that used to run after it -- one per solver
      ! branch -- would have re-packed block ghosts the compressed halo no
      ! longer updates, straight over the haloed values.
      !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
      do iblk = 1, nblocks
         call pack_to_slots(strength(:,:,iblk), dynMap(:,iblk), strength_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
      enddo
      !$OMP END PARALLEL DO

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate2DR8_compressed (strength_slot,    halo_info,         &
                                          field_loc_center, field_type_scalar, &
                                          dynMap,           dynOff)

      ! velocities may have changed in dyn_prep2
      call stack_fields(uvel, vvel, fld2)
      call ice_HaloUpdate (fld2,               halo_info, &
                           field_loc_NEcorner, field_type_vector)
      call unstack_fields(fld2, uvel, vvel)
      call ice_timer_stop(timer_bound)

      if (maskhalo_dyn) then
         halomask = 0
         if (grid_ice == 'B') then
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  if (iceUmask((j-1)*nx_block+i,iblk)) halomask(i,j,iblk) = 1
               enddo
               enddo
            enddo
         elseif (grid_ice == 'C' .or. grid_ice == 'CD') then
            !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,i,j,ij) SCHEDULE(runtime)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j = jlo,jhi
               do i = ilo,ihi
                  ij = (j-1)*nx_block + i
                  if (iceTmask(ij         ,iblk) .or. &
                      iceTmask(ij-1       ,iblk) .or. &
                      iceTmask(ij+1       ,iblk) .or. &
                      iceTmask(ij-nx_block,iblk) .or. &
                      iceTmask(ij+nx_block,iblk)) then
                     halomask(i,j,iblk) = 1
                  endif
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
         endif
         call ice_timer_start(timer_bound)
         call ice_HaloUpdate (halomask,         halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_timer_stop(timer_bound)
         call ice_HaloMask(halo_info_mask, halo_info, halomask)
      endif

      !-----------------------------------------------------------------
      ! seabed stress factor TbU (TbU is part of Cb coefficient)
      !-----------------------------------------------------------------

      if (seabed_stress) then

         if (grid_ice == "B") then

            if ( seabed_stress_method == 'LKD' ) then
               call seabed_stress_factor_LKD (1               , dynTot, iceUmask_slot        , &
                                              NbrE_Task       , NbrNE_Task    , &
                                              NbrN_Task       ,                 &
                                              vice_slot       , aice_slot     , &
                                              hwater_slot     , TbU_slot)

            elseif ( seabed_stress_method == 'probabilistic' ) then
               call seabed_stress_factor_prob (1               , dynTot        , &
                                               iceTmask_slot           , iceUmask_slot         , &
                                               NbrE_Task       , NbrNE_Task    , &
                                               NbrN_Task       ,                 &
                                               aicen_slot      , vicen_slot    , &
                                               hwater_slot     , TbU_slot)
            endif

         elseif (grid_ice == "C" .or. grid_ice == "CD") then

            if ( seabed_stress_method == 'LKD' ) then
               call seabed_stress_factor_LKD (1               , dynTot, iceEmask_slot        , &
                                              NbrE_Task       , NbrNE_Task    , &
                                              NbrN_Task       ,                 &
                                              vice_slot       , aice_slot     , &
                                              hwater_slot     , TbE_slot, &
                                              'E')
               call seabed_stress_factor_LKD (1               , dynTot, iceNmask_slot        , &
                                              NbrE_Task       , NbrNE_Task    , &
                                              NbrN_Task       ,                 &
                                              vice_slot       , aice_slot     , &
                                              hwater_slot     , TbN_slot, &
                                              'N')

            elseif ( seabed_stress_method == 'probabilistic' ) then
               call seabed_stress_factor_prob (1               , dynTot        , &
                                               iceTmask_slot           , iceUmask_slot         , &
                                               NbrE_Task       , NbrNE_Task    , &
                                               NbrN_Task       ,                 &
                                               aicen_slot      , vicen_slot    , &
                                               hwater_slot     , TbU_slot      , &
                                               TbE_slot        , TbN_slot      , &
                                               iceEmask_slot           , iceNmask_slot)
            endif

         endif

      endif

      call ice_timer_start(timer_evp)

      if (grid_ice == "B") then

         ! ⚠ HOISTED above the evp_algorithm branch so BOTH solvers consume the
         ! same compressed inputs.  It used to sit inside the standard_2d arm,
         ! which left these slots stale on the evp1d path -- harmless while
         ! evp1d read only the block arrays, and a silent wrong answer the
         ! moment it reads a slot instead.  Homogenising means one producer.
         ! Pack the velocities into the compressed layout once, subcycle on
         ! them there, and unpack once at the end.  uvel/vvel themselves
         ! stay in the block layout for restart, history, transport and the
         ! couplers; only the ndte subcycles run compressed.
         !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
         do iblk = 1, nblocks
            call pack_to_slots(uvel(:,iblk), dynMap(:,iblk), uvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vvel(:,iblk), dynMap(:,iblk), vvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            ! constant through the subcycling, so packed once and never
            ! unpacked; the block-layout originals stay for the consumers
            ! outside the dynamics
            ! umassdti, waterxU, wateryU, forcexU and forceyU are not packed
            ! here: dyn_prep2 already filled their compressed copies in the
            ! same loop that computed them, which costs one store instead of
            ! a second traversal of the block.
            ! strength is constant through the subcycling; the twelve
            ! stress components are the subcycle's own state, read and
            ! written every pass and unpacked once at the end for the
            ! restart.
            ! ⚠ The stresses are NOT packed from the block arrays any more.  The
            ! slot copies are the live state: they persist between timesteps,
            ! dyn_prep2 applies the no-ice reset to them directly, and the
            ! restart read seeds them through stresses_to_slots.  Packing here
            ! would overwrite all of that with block values nothing maintains.
         enddo
         !$OMP END PARALLEL DO

         if (evp_algorithm == "shared_mem_1d" ) then

            call dyn_evp1d_run(dt        ,                                     &
                               stressp_1 , stressp_2, stressp_3  , stressp_4 , &
                               stressm_1 , stressm_2 , stressm_3 , stressm_4 , &
                               stress12_1, stress12_2, stress12_3, stress12_4, &
                               uvel      , vvel      ,              &
                               iceTmask  , iceUmask)

            ! scatter_dyn leaves uvel_slot and vvel_slot current along with
            ! the block forms, which is what the shared deformations call
            ! below reads.

            ! ⚠ The compressed fold writes only cells that OWN A SLOT.  The
            ! ghost columns and the land cells on the seam row have none, so
            ! they keep what the scatter left -- +0.0 -- where the 2-D path's
            ! per-subcycle halo stamps -0.0 into them.  Numerically identical,
            ! bitwise not, and md5/complog both see it: on gx3 it was exactly
            ! ten words of uvel/vvel at five slotless cells.
            !
            ! Running the real handler once makes the two paths agree BY
            ! CONSTRUCTION rather than by reproducing its sign-of-zero
            ! convention in the compressed fold.  One halo update per
            ! timestep, against ndte=240 subcycles in the branch above.
            !
            ! Slotted ghost cells are unaffected: the compressed halo has
            ! already given them the folded value, so the handler writes back
            ! what is there and uvel_slot/vvel_slot cannot go stale.  Only the
            ! slotless cells change, and those have no slot to keep current.
            if (trim(ns_boundary_type) == 'tripole' .or. &
                trim(ns_boundary_type) == 'tripoleT') then
               call ice_timer_start(timer_bound)
               call ice_HaloUpdate (uvel, halo_info, &
                                    field_loc_NEcorner, field_type_vector)
               call ice_HaloUpdate (vvel, halo_info, &
                                    field_loc_NEcorner, field_type_vector)
               call ice_timer_stop(timer_bound)
            endif

         else ! evp_algorithm == standard_2d (Standard CICE)


            ! Dead slots too: stepu reads them for absent neighbours.
            ! Once per timestep, not per subcycle -- stress writes all eight
            ! unconditionally at every non-iceTmask_slot slot, and iceTmask_slot is fixed for
            ! the whole timestep, so passes 2..ndte only rewrote zeros over
            ! zeros.  It cannot move out of the timestep: iceTmask_slot changes as the
            ! ice edge moves, and a slot that was active last step would hand
            ! stepu a stale value.
            str1_slot(:) = c0
            str2_slot(:) = c0
            str3_slot(:) = c0
            str4_slot(:) = c0
            str5_slot(:) = c0
            str6_slot(:) = c0
            str7_slot(:) = c0
            str8_slot(:) = c0

            do ksub = 1,ndte        ! subcycling

               !-----------------------------------------------------------------
               ! stress tensor equation, total surface stress
               !
               ! One task-wide call: every argument is already a whole array, so
               ! the block loop only supplied (lo,hi).  Collapsing it is exact --
               ! the slots the per-block ranges omit are the blocks' dead slots,
               ! and those are permanently iceTmask_slot.  Threading moved inside.
               !-----------------------------------------------------------------
                  call stress (1, dynTot, iceTmask_slot, &
                               uvel_slot, vvel_slot, &
                               dxT_slot, dyT_slot, &
                               dxhy_slot, dyhx_slot, &
                               cxp_slot, cyp_slot, &
                               cxm_slot, cym_slot, &
                               DminTarea_slot, &
                               strength_slot, &
                               stressp_1_slot, stressp_2_slot, &
                               stressp_3_slot, stressp_4_slot, &
                               stressm_1_slot, stressm_2_slot, &
                               stressm_3_slot, stressm_4_slot, &
                               stress12_1_slot, stress12_2_slot, &
                               stress12_3_slot, stress12_4_slot, &
                               NbrW_Task, NbrSW_Task, NbrS_Task, &
                               str1_slot, str2_slot, &
                               str3_slot, str4_slot, &
                               str5_slot, str6_slot, &
                               str7_slot, str8_slot)

               !-----------------------------------------------------------------
               ! momentum equation.  Task-wide, as stress above; this was the
               ! last kernel in the block loop, so the loop is gone.
               !-----------------------------------------------------------------
                  call stepu (&
                              1, dynTot,                     &
                              cdn_ocnU_slot,  &
                              iceUmask_slot, &
                              aiU_slot , &
                              str1_slot, str2_slot, &
                              str3_slot, str4_slot, &
                              str5_slot, str6_slot, &
                              str7_slot, str8_slot, &
                              uocnU_slot, vocnU_slot, &
                              waterxU_slot, wateryU_slot, &
                              forcexU_slot, forceyU_slot, &
                              umassdti_slot, fmU_slot, &
                              uarear_slot,    &
                              taubxU_slot  , &
                              uvel_init,   vvel_init,  &
                              uvel_slot, vvel_slot, &
                              NbrE_Task, NbrNE_Task, NbrN_Task, &
                              TbU_slot , &
                              rheofactU_slot)

               ! U fields at NE corner.  Exchanged in the compressed layout,
               ! so the velocities stay compressed across all ndte subcycles,
               ! and bundled, so the pair costs one message rather than two --
               ! the compressed counterpart of what stack_fields does for the
               ! block layout.
               call ice_timer_start(timer_bound)
               call stack_fields(uvel_slot, vvel_slot, fld2c)
               if (maskhalo_dyn) then
                  call ice_HaloUpdate3DR8_compressed(fld2c, halo_info_mask, &
                       field_loc_NEcorner, field_type_vector, dynMap, dynOff)
               else
                  call ice_HaloUpdate3DR8_compressed(fld2c, halo_info, &
                       field_loc_NEcorner, field_type_vector, dynMap, dynOff)
               endif
               call unstack_fields(fld2c, uvel_slot, vvel_slot)
               call ice_timer_stop(timer_bound)

               call ice_restoring_halo(setfld='velocity', &
                                       uvel_s=uvel_slot, vvel_s=vvel_slot)

            enddo  ! sub cycling

            !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
            do iblk = 1, nblocks
               call unpack_from_slots(uvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), uvel(:,iblk))
               call unpack_from_slots(vvel_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), vvel(:,iblk))
               ! Both diagnostics are built here, once, rather than on
               ! every one of the ndte passes: str1..str8 still hold what
               ! the final pass used, and taubxU_slot holds its Cb.
               ! strintxU/strintyU and taubxU/taubyU are outputs only; stepu
               ! writes strint to slots every pass and this is where the block
               ! arrays history reads get their one and only write -- and only
               ! when a stream is going to read them.
               if (lstrint) then
                  call calc_strint_slot(dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk), iceUmask_slot, &
                       uarear_slot, &
                       rheofactU_slot, &
                       str1_slot, str2_slot, &
                       str3_slot, str4_slot, &
                       str5_slot, str6_slot, &
                       str7_slot, str8_slot, &
                       NbrE_Task, NbrNE_Task, NbrN_Task, &
                       strintxU_slot, &
                       strintyU_slot)
                  call unpack_from_slots(strintxU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintxU(:,iblk))
                  call unpack_from_slots(strintyU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintyU(:,iblk))
               endif
               if (ltaub) then
                  call calc_seabed_stress_slot(dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk), iceUmask_slot, &
                       uvel_slot, &
                       taubxU_slot, &
                       vvel_slot, &
                       taubyU_slot)
                  call unpack_from_slots(taubxU_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubxU  (:,iblk))
                  call unpack_from_slots(taubyU_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubyU  (:,iblk))
               endif
            enddo
            !$OMP END PARALLEL DO
         endif ! evp algorithm

         !-----------------------------------------------------------------
         ! save quantities for mechanical redistribution
         !-----------------------------------------------------------------

            call deformations (                                                                                                              &
                               1, dynTot, iceTmask_slot,                                                &
                               NbrW_Task, NbrSW_Task, NbrS_Task, dynIJ,          &
                               uvel_slot, vvel_slot, &
                               dxT_slot, dyT_slot,   &
                               dxU_slot, dyU_slot,   &
                               cxp_slot, cyp_slot,   &
                               cxm_slot, cym_slot,   &
                               tarear_slot, vort,                                      &
                               shear   , divu     ,                                                                      &
                               rdg_conv, rdg_shear )

      elseif (grid_ice == "C") then

         ! ------------------------------------------------------------------
         ! Pack once, subcycle compressed, unpack once -- the shape the EAP
         ! conversion uses.  The block arrays stay authoritative outside the
         ! subcycle, which is what restart, history and the couplers read.
         !
         ! shearU is not packed: strain_rates_U is the only thing
         ! that writes it and it writes straight into the slots.
         ! ------------------------------------------------------------------
         !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
         do iblk = 1, nblocks
            call pack_to_slots(uvelE(:,:,iblk), dynMap(:,iblk), uvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vvelE(:,:,iblk), dynMap(:,iblk), vvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(uvelN(:,:,iblk), dynMap(:,iblk), uvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vvelN(:,:,iblk), dynMap(:,iblk), vvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            ! ⚠ uvel and vvel too.  strain_rates_U reads them on the FIRST
            ! subcycle pass -- as uvelU/vvelU -- before the grid_average
            ! shifts write them later in the same pass.  Without this the
            ! first pass used whatever uvel_slot held from the last timestep,
            ! which is the previous values in a continuous run and zeros
            ! after a restart: the restart test could not reproduce.
            call pack_to_slots(uvel(:,iblk), dynMap(:,iblk), uvel_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vvel(:,iblk), dynMap(:,iblk), vvel_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            ! constant through the subcycling, so packed once
            ! The step kernels read these every pass but nothing changes
            ! them during the subcycle, so one pack each is enough.
         enddo
         !$OMP END PARALLEL DO

         do ksub = 1,ndte        ! subcycling

            !-----------------------------------------------------------------
            ! strain rates at U point
            ! NOTE these are actually strain rates * area  (m^2/s)
            !-----------------------------------------------------------------
            call strain_rates_U (1             , dynTot        , &
                                 iceUmask_slot         ,                 &
                                 NbrE_Task, NbrN_Task      ,                 &
                                 uvelE_slot    , vvelE_slot    , &
                                 uvelN_slot    , vvelN_slot    , &
                                 uvel_slot     , vvel_slot     , &
                                 dxE_slot      , dyN_slot      , &
                                 dxU_slot      , dyU_slot      , &
                                 ratiodxN_slot , ratiodxNr_slot, &
                                 ratiodyE_slot , ratiodyEr_slot, &
                                 epm_slot      , npm_slot      , &
                                 divergU       , tensionU      , &
                                 shearU_slot   , deltaU          )

            call ice_timer_start(timer_bound)
            if (maskhalo_dyn) then
               call ice_HaloUpdate2DR8_compressed(shearU_slot, halo_info_mask, &
                    field_loc_NEcorner, field_type_scalar, dynMap, dynOff)
            else
               call ice_HaloUpdate2DR8_compressed(shearU_slot, halo_info, &
                    field_loc_NEcorner, field_type_scalar, dynMap, dynOff)
            endif
            call ice_timer_stop(timer_bound)

               call stressC_T (1, dynTot       , iceTmask_slot, &
                               NbrW_Task, NbrSW_Task, NbrS_Task, &
                               uvelE_slot, &
                               vvelN_slot, &
                               dxN_slot  , &
                               dyE_slot  , &
                               dxT_slot, dyT_slot, &
                               uarea_slot   , &
                               uareaavg_slot, &
                               DminTarea_slot, &
                               strength_slot, &
                               shearU_slot, &
                               zetax2T_slot, &
                               etax2T_slot , &
                               stresspT_slot , &
                               stressmT_slot , &
                               stress12T_slot)


            ! All four go compressed in ONE exchange.  zetax2T and etax2T
            ! have no block form any more -- nothing outside the dynamics
            ! reads them -- so the block halo that used to carry them is gone.
            call ice_timer_start(timer_bound)
            call stack_fields(zetax2T_slot, etax2T_slot, &
                              stresspT_slot, stressmT_slot, fld4c)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld4c, halo_info_mask, &
                    field_loc_center, field_type_scalar, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld4c, halo_info, &
                    field_loc_center, field_type_scalar, dynMap, dynOff)
            endif
            call unstack_fields(fld4c, zetax2T_slot, etax2T_slot, &
                                stresspT_slot, stressmT_slot)
            call ice_timer_stop(timer_bound)

            ! The T-to-U shift on the compressed layout.  'NE' is the compass
            ! direction of the T2US stencil; hm is a real 0/1 mask there, not
            ! a logical, and tarea is its weight.
            if (visc_method == 'avg_strength') then
               call grid_average_X2YS_compressed('NE', strength_slot, &
                    tarea_slot, hm_slot, wsumNE_slot, strengthU_slot)
            elseif (visc_method == 'avg_zeta') then
               call grid_average_X2YS_compressed('NE', etax2T_slot, &
                    tarea_slot, hm_slot, wsumNE_slot, etax2U_slot)
            endif

            call stressC_U (1             , dynTot        , &
                            iceUmask_slot         ,                 &
                            uarea_slot    , etax2U_slot   , &
                            deltaU        , strengthU_slot, &
                            shearU_slot   , stress12U_slot)

            call ice_timer_start(timer_bound)
            if (maskhalo_dyn) then
               call ice_HaloUpdate2DR8_compressed(stress12U_slot, halo_info_mask, &
                    field_loc_NEcorner, field_type_scalar, dynMap, dynOff)
            else
               call ice_HaloUpdate2DR8_compressed(stress12U_slot, halo_info, &
                    field_loc_NEcorner, field_type_scalar, dynMap, dynOff)
            endif
            call ice_timer_stop(timer_bound)


               call div_stress_Ex (1, dynTot, iceEmask_slot, &
                                   NbrS_Task, NbrE_Task, &
                                   dxE_slot, dyE_slot, &
                                   dxU_slot  , &
                                   dyT_slot  , &
                                   earear_slot, rheofactE_slot, &
                                   stresspT_slot , &
                                   stressmT_slot , &
                                   stress12U_slot, &
                                   strintxE_slot )

               call div_stress_Ny (1, dynTot, iceNmask_slot, &
                                   NbrW_Task, NbrN_Task, &
                                   dxN_slot, dyN_slot, &
                                   dxT_slot  , &
                                   dyU_slot  , &
                                   narear_slot, rheofactN_slot, &
                                   stresspT_slot , &
                                   stressmT_slot , &
                                   stress12U_slot, &
                                   strintyN_slot )


                call stepu_C (1            , dynTot       , & ! u, E point
                              iceEmask_slot        , cdn_ocnE_slot, &
                              aiE_slot     , uocnE_slot   , &
                              vocnE_slot   , waterxE_slot , &
                              forcexE_slot , emassdti_slot, &
                              fmE_slot     , strintxE_slot, &
                              taubxE_slot  , uvelE_init   , &
                              uvelE_slot   , vvelE_slot   , &
                              TbE_slot)

                call stepv_C (1            , dynTot       , & ! v, N point
                              iceNmask_slot        , cdn_ocnN_slot, &
                              aiN_slot     , uocnN_slot   , &
                              vocnN_slot   , wateryN_slot , &
                              forceyN_slot , nmassdti_slot, &
                              fmN_slot     , strintyN_slot, &
                              taubyN_slot  , vvelN_init   , &
                              uvelN_slot   , vvelN_slot   , &
                              TbN_slot)

            ! Shifts, masks and halo exchanges all on the compressed layout,
            ! so the velocities never leave slot form inside the subcycle.
            ! Directions are each stencil's compass sense: E2NA is NW, N2EA is
            ! SE, E2UA is N, N2UA is E; the weight is the source cell's area.
            call ice_timer_start(timer_bound)
            if (maskhalo_dyn) then
               call ice_HaloUpdate2DR8_compressed(uvelE_slot, halo_info_mask, &
                    field_loc_Eface, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate2DR8_compressed(uvelE_slot, halo_info, &
                    field_loc_Eface, field_type_vector, dynMap, dynOff)
            endif
            if (maskhalo_dyn) then
               call ice_HaloUpdate2DR8_compressed(vvelN_slot, halo_info_mask, &
                    field_loc_Nface, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate2DR8_compressed(vvelN_slot, halo_info, &
                    field_loc_Nface, field_type_vector, dynMap, dynOff)
            endif
            call ice_timer_stop(timer_bound)

            call grid_average_X2YA_compressed('NW', uvelE_slot, earea_slot, wsumNW_slot, uvelN_slot)
            call grid_average_X2YA_compressed('SE', vvelN_slot, narea_slot, wsumSE_slot, vvelE_slot)
            uvelN_slot(:) = uvelN_slot(:)*npm_slot(:)
            vvelE_slot(:) = vvelE_slot(:)*epm_slot(:)

            call ice_timer_start(timer_bound)
            if (maskhalo_dyn) then
               call ice_HaloUpdate2DR8_compressed(uvelN_slot, halo_info_mask, &
                    field_loc_Nface, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate2DR8_compressed(uvelN_slot, halo_info, &
                    field_loc_Nface, field_type_vector, dynMap, dynOff)
            endif
            if (maskhalo_dyn) then
               call ice_HaloUpdate2DR8_compressed(vvelE_slot, halo_info_mask, &
                    field_loc_Eface, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate2DR8_compressed(vvelE_slot, halo_info, &
                    field_loc_Eface, field_type_vector, dynMap, dynOff)
            endif
            call ice_timer_stop(timer_bound)

            call grid_average_X2YA_compressed('N', uvelE_slot, earea_slot, wsumN_slot, uvel_slot)
            call grid_average_X2YA_compressed('E', vvelN_slot, narea_slot, wsumE_slot, vvel_slot)

            uvel_slot(:) = uvel_slot(:)*uvm_slot(:)
            vvel_slot(:) = vvel_slot(:)*uvm_slot(:)

            ! U fields at NE corner, bundled into one compressed exchange
            call ice_timer_start(timer_bound)
            call stack_fields(uvel_slot, vvel_slot, fld2c)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info_mask, &
                    field_loc_NEcorner, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info, &
                    field_loc_NEcorner, field_type_vector, dynMap, dynOff)
            endif
            call unstack_fields(fld2c, uvel_slot, vvel_slot)
            call ice_timer_stop(timer_bound)

            call ice_restoring_halo(setfld='velocity', &
                                    uvel_s=uvel_slot, vvel_s=vvel_slot)

            ! ⚠ stepu_C/stepv_C and the four grid_average_X2Y shifts above

         enddo                     ! subcycling

         ! Unpack what the subcycle carried compressed.  deformationsC_T below
         ! reads the block arrays, so this has to come first.
         !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
         do iblk = 1, nblocks
            call unpack_from_slots(shearU_slot   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), shearU   (:,iblk))
            call unpack_from_slots(uvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), uvelE(:,:,iblk))
            call unpack_from_slots(vvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), vvelE(:,:,iblk))
            call unpack_from_slots(uvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), uvelN(:,:,iblk))
            call unpack_from_slots(vvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), vvelN(:,:,iblk))
            call unpack_from_slots(uvel_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), uvel(:,iblk))
            call unpack_from_slots(vvel_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), vvel(:,iblk))
            ! strint* and taub* are outputs only: the kernels write slots
            ! every pass and this is their one and only block write.
            ! Cb was carried in taubxE_slot / taubyN_slot through the
            ! subcycling; the stresses are formed from it here, once.
            if (lstrint) then
               call unpack_from_slots(strintxE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintxE(:,:,iblk))
               call unpack_from_slots(strintyN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintyN(:,:,iblk))
            endif
            if (ltaub) then
               call calc_seabed_stress_slot(dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk), iceEmask_slot, &
                    uvelE_slot, &
                    taubxE_slot)
               call calc_seabed_stress_slot(dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk), iceNmask_slot, &
                    vvelN_slot, &
                    taubyN_slot)
               call unpack_from_slots(taubxE_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubxE  (:,iblk))
               call unpack_from_slots(taubyN_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubyN  (:,iblk))
            endif
         enddo
         !$OMP END PARALLEL DO

         !-----------------------------------------------------------------
         ! save quantities for mechanical redistribution
         !-----------------------------------------------------------------

            call deformationsC_T (1, dynTot, iceTmask_slot,               &
                                  NbrW_Task, NbrSW_Task, NbrS_Task, dynIJ,           &
                                  uvelE_slot, vvelE_slot,  &
                                  uvelN_slot, vvelN_slot,  &
                                  dxN_slot, dyE_slot,      &
                                  dxT_slot, dyT_slot,      &
                                  tarear_slot, uarea_slot, &
                                  uareaavg_slot,                                                      &
                                  shearU_slot, vort,    &
                                  shear   , divu     , &
                                  rdg_conv, rdg_shear)

      elseif (grid_ice == "CD") then

         ! Pack once, subcycle compressed, unpack once -- see the C branch.
         !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
         do iblk = 1, nblocks
            call pack_to_slots(uvelE(:,:,iblk), dynMap(:,iblk), uvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vvelE(:,:,iblk), dynMap(:,iblk), vvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(uvelN(:,:,iblk), dynMap(:,iblk), uvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vvelN(:,:,iblk), dynMap(:,iblk), vvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            ! ⚠ uvel and vvel too.  strain_rates_U reads them on the FIRST
            ! subcycle pass -- as uvelU/vvelU -- before the grid_average
            ! shifts write them later in the same pass.  Without this the
            ! first pass used whatever uvel_slot held from the last timestep,
            ! which is the previous values in a continuous run and zeros
            ! after a restart: the restart test could not reproduce.
            call pack_to_slots(uvel(:,iblk), dynMap(:,iblk), uvel_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vvel(:,iblk), dynMap(:,iblk), vvel_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            ! constant through the subcycling, so packed once
            ! The step kernels read these every pass but nothing changes
            ! them during the subcycle, so one pack each is enough.
         enddo
         !$OMP END PARALLEL DO

         do ksub = 1,ndte        ! subcycling

               call stressCD_T (1, dynTot       , iceTmask_slot, &
                                NbrW_Task, NbrS_Task, &
                                uvelE_slot, &
                                vvelE_slot, &
                                uvelN_slot, &
                                vvelN_slot, &
                                dxN_slot  , &
                                dyE_slot  , &
                                dxT_slot, dyT_slot, &
                                                      DminTarea_slot, &
                                strength_slot, &
                                zetax2T_slot , &
                                etax2T_slot  , &
                                stresspT_slot , &
                                stressmT_slot , &
                                stress12T_slot )


            ! calls ice_haloUpdate, controls bundles and masks
            call ice_timer_start(timer_bound)
            call stack_fields(zetax2T_slot, etax2T_slot, fld2c)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info_mask, &
                    field_loc_center, field_type_scalar, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info, &
                    field_loc_center, field_type_scalar, dynMap, dynOff)
            endif
            call unstack_fields(fld2c, zetax2T_slot, etax2T_slot)
            call ice_timer_stop(timer_bound)

            ! T-to-U on the compressed layout; 'NE' is the T2US stencil's
            ! compass direction, tarea its weight and hm its (real 0/1) mask.
            if (visc_method == 'avg_strength') then
               call grid_average_X2YS_compressed('NE', strength_slot, &
                    tarea_slot, hm_slot, wsumNE_slot, strengthU_slot)
            elseif (visc_method == 'avg_zeta') then
               call grid_average_X2YS_compressed('NE', zetax2T_slot, &
                    tarea_slot, hm_slot, wsumNE_slot, zetax2U_slot)
               call grid_average_X2YS_compressed('NE', etax2T_slot, &
                    tarea_slot, hm_slot, wsumNE_slot, etax2U_slot)
            endif

            !-----------------------------------------------------------------
            ! strain rates at U point
            ! NOTE these are actually strain rates * area  (m^2/s)
            !-----------------------------------------------------------------
            call strain_rates_U (1             , dynTot        , &
                                 iceUmask_slot         ,                 &
                                 NbrE_Task, NbrN_Task      ,                 &
                                 uvelE_slot    , vvelE_slot    , &
                                 uvelN_slot    , vvelN_slot    , &
                                 uvel_slot     , vvel_slot     , &
                                 dxE_slot      , dyN_slot      , &
                                 dxU_slot      , dyU_slot      , &
                                 ratiodxN_slot , ratiodxNr_slot, &
                                 ratiodyE_slot , ratiodyEr_slot, &
                                 epm_slot      , npm_slot      , &
                                 divergU       , tensionU      , &
                                 shearU_slot   , DeltaU          )

            call stressCD_U     (1             , dynTot        , &
                                 iceUmask_slot         ,                 &
                                 uarea_slot    ,                 &
                                 zetax2U_slot  , etax2U_slot   , &
                                 strengthU_slot,                 &
                                 divergU       , tensionU      , &
                                 shearU_slot   , DeltaU        , &
                                 stresspU_slot , stressmU_slot , &
                                 stress12U_slot)

            ! calls ice_haloUpdate, controls bundles and masks
            ! Six stresses, two compressed exchanges -- one per field
            ! location, which is why they cannot be merged into one.
            call ice_timer_start(timer_bound)
            call stack_fields(stresspT_slot, stressmT_slot, &
                              stress12T_slot, fld3c)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld3c, halo_info_mask, &
                    field_loc_center, field_type_scalar, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld3c, halo_info, &
                    field_loc_center, field_type_scalar, dynMap, dynOff)
            endif
            call unstack_fields(fld3c, stresspT_slot, stressmT_slot, &
                                stress12T_slot)

            call stack_fields(stresspU_slot, stressmU_slot, &
                              stress12U_slot, fld3c)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld3c, halo_info_mask, &
                    field_loc_NEcorner, field_type_scalar, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld3c, halo_info, &
                    field_loc_NEcorner, field_type_scalar, dynMap, dynOff)
            endif
            call unstack_fields(fld3c, stresspU_slot, stressmU_slot, &
                                stress12U_slot)
            call ice_timer_stop(timer_bound)


               call div_stress_Ex (1, dynTot, iceEmask_slot, &
                                   NbrS_Task, NbrE_Task, &
                                   dxE_slot, dyE_slot, &
                                   dxU_slot  , &
                                   dyT_slot  , &
                                   earear_slot, rheofactE_slot, &
                                   stresspT_slot , &
                                   stressmT_slot , &
                                   stress12U_slot, &
                                   strintxE_slot )

               call div_stress_Ey (1, dynTot, iceEmask_slot, &
                                   NbrS_Task, NbrE_Task, &
                                   dxE_slot, dyE_slot, &
                                   dxU_slot  , &
                                   dyT_slot  , &
                                   earear_slot, rheofactE_slot, &
                                   stresspU_slot , &
                                   stressmU_slot , &
                                   stress12T_slot, &
                                   strintyE_slot )

               call div_stress_Nx (1, dynTot, iceNmask_slot, &
                                   NbrW_Task, NbrN_Task, &
                                   dxN_slot, dyN_slot, &
                                   dxT_slot  , &
                                   dyU_slot  , &
                                   narear_slot, rheofactN_slot, &
                                   stresspU_slot , &
                                   stressmU_slot , &
                                   stress12T_slot, &
                                   strintxN_slot )

               call div_stress_Ny (1, dynTot, iceNmask_slot, &
                                   NbrW_Task, NbrN_Task, &
                                   dxN_slot, dyN_slot, &
                                   dxT_slot  , &
                                   dyU_slot  , &
                                   narear_slot, rheofactN_slot, &
                                   stresspT_slot , &
                                   stressmT_slot , &
                                   stress12U_slot, &
                                   strintyN_slot )


               call stepuv_CD (1            , dynTot       , & ! E point
                               iceEmask_slot        , cdn_ocnE_slot, &
                               aiE_slot     , uocnE_slot   , &
                               vocnE_slot   , waterxE_slot , &
                               wateryE_slot , forcexE_slot , &
                               forceyE_slot , emassdti_slot, &
                               fmE_slot     , strintxE_slot, &
                               strintyE_slot, taubxE_slot  , &
                               uvelE_init   , vvelE_init   , &
                               uvelE_slot   , vvelE_slot   , &
                               TbE_slot)

               call stepuv_CD (1            , dynTot       , & ! N point
                               iceNmask_slot        , cdn_ocnN_slot, &
                               aiN_slot     , uocnN_slot   , &
                               vocnN_slot   , waterxN_slot , &
                               wateryN_slot , forcexN_slot , &
                               forceyN_slot , nmassdti_slot, &
                               fmN_slot     , strintxN_slot, &
                               strintyN_slot, taubxN_slot  , &
                               uvelN_init   , vvelN_init   , &
                               uvelN_slot   , vvelN_slot   , &
                               TbN_slot)

            ! CD evolves all four velocities independently, so unlike the C
            ! branch there is no E2NA/N2EA pair here -- only the two shifts
            ! that build the U-point velocities.  Each face pair goes in one
            ! compressed bundle.
            call ice_timer_start(timer_bound)
            call stack_fields(uvelE_slot, vvelE_slot, fld2c)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info_mask, &
                    field_loc_Eface, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info, &
                    field_loc_Eface, field_type_vector, dynMap, dynOff)
            endif
            call unstack_fields(fld2c, uvelE_slot, vvelE_slot)
            call stack_fields(uvelN_slot, vvelN_slot, fld2c)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info_mask, &
                    field_loc_Nface, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info, &
                    field_loc_Nface, field_type_vector, dynMap, dynOff)
            endif
            call unstack_fields(fld2c, uvelN_slot, vvelN_slot)
            call ice_timer_stop(timer_bound)

            call grid_average_X2YA_compressed('N', uvelE_slot, earea_slot, wsumN_slot, uvel_slot)
            call grid_average_X2YA_compressed('E', vvelN_slot, narea_slot, wsumE_slot, vvel_slot)

            uvel_slot(:) = uvel_slot(:)*uvm_slot(:)
            vvel_slot(:) = vvel_slot(:)*uvm_slot(:)

            ! U fields at NE corner, bundled into one compressed exchange
            call ice_timer_start(timer_bound)
            call stack_fields(uvel_slot, vvel_slot, fld2c)
            if (maskhalo_dyn) then
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info_mask, &
                    field_loc_NEcorner, field_type_vector, dynMap, dynOff)
            else
               call ice_HaloUpdate3DR8_compressed(fld2c, halo_info, &
                    field_loc_NEcorner, field_type_vector, dynMap, dynOff)
            endif
            call unstack_fields(fld2c, uvel_slot, vvel_slot)
            call ice_timer_stop(timer_bound)

            call ice_restoring_halo(setfld='velocity', &
                                    uvel_s=uvel_slot, vvel_s=vvel_slot)


         enddo                     ! subcycling

         ! Unpack what the subcycle carried compressed; deformationsCD_T below
         ! reads the block arrays.
         !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
         do iblk = 1, nblocks
            call unpack_from_slots(shearU_slot   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), shearU   (:,iblk))
            call unpack_from_slots(uvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), uvelE(:,:,iblk))
            call unpack_from_slots(vvelE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), vvelE(:,:,iblk))
            call unpack_from_slots(uvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), uvelN(:,:,iblk))
            call unpack_from_slots(vvelN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), vvelN(:,:,iblk))
            call unpack_from_slots(uvel_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), uvel(:,iblk))
            call unpack_from_slots(vvel_slot (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), vvel(:,iblk))
            ! strint* and taub* are outputs only: the kernels write slots
            ! every pass and this is their one and only block write.
            ! Cb was carried in taubxE_slot / taubxN_slot through the
            ! subcycling; both components are formed from it here, once.
            if (lstrint) then
               call unpack_from_slots(strintxE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintxE(:,:,iblk))
               call unpack_from_slots(strintyE_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintyE(:,:,iblk))
               call unpack_from_slots(strintxN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintxN(:,:,iblk))
               call unpack_from_slots(strintyN_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), strintyN(:,:,iblk))
            endif
            if (ltaub) then
               call calc_seabed_stress_slot(dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk), iceEmask_slot, &
                    uvelE_slot, &
                    taubxE_slot, &
                    vvelE_slot, &
                    taubyE_slot)
               call calc_seabed_stress_slot(dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk), iceNmask_slot, &
                    uvelN_slot, &
                    taubxN_slot, &
                    vvelN_slot, &
                    taubyN_slot)
               call unpack_from_slots(taubxE_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubxE  (:,iblk))
               call unpack_from_slots(taubyE_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubyE  (:,iblk))
               call unpack_from_slots(taubxN_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubxN  (:,iblk))
               call unpack_from_slots(taubyN_slot  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), taubyN  (:,iblk))
            endif
         enddo
         !$OMP END PARALLEL DO

         !-----------------------------------------------------------------
         ! save quantities for mechanical redistribution
         !-----------------------------------------------------------------

            call deformationsCD_T (1, dynTot, iceTmask_slot,              &
                                   NbrW_Task, NbrS_Task, dynIJ,          &
                                   uvelE_slot, vvelE_slot, &
                                   uvelN_slot, vvelN_slot, &
                                   dxN_slot, dyE_slot,     &
                                   dxT_slot, dyT_slot,     &
                                   tarear_slot, vort,   &
                                   shear   , divu     , &
                                   rdg_conv, rdg_shear)
      endif   ! grid_ice

      call ice_timer_stop(timer_evp)

      if (maskhalo_dyn) then
         call ice_HaloDestroy(halo_info_mask)
      endif

      ! Force symmetry across the tripole seam
      if (trim(ns_boundary_type) == 'tripole' .or. &
          trim(ns_boundary_type) == 'tripoleT') then
         ! TODO: C/CD-grid
      if (maskhalo_dyn) then
         !-------------------------------------------------------
         ! set halomask to zero because ice_HaloMask always keeps
         ! local copies AND tripole zipper communication
         !-------------------------------------------------------
         halomask = 0
         call ice_HaloMask(halo_info_mask, halo_info, halomask)

         call ice_HaloUpdate_stress_compressed(stressp_1_slot, stressp_3_slot, halo_info_mask, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stressp_3_slot, stressp_1_slot, halo_info_mask, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stressp_2_slot, stressp_4_slot, halo_info_mask, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stressp_4_slot, stressp_2_slot, halo_info_mask, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)

         call ice_HaloUpdate_stress_compressed(stressm_1_slot, stressm_3_slot, halo_info_mask, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stressm_3_slot, stressm_1_slot, halo_info_mask, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stressm_2_slot, stressm_4_slot, halo_info_mask, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stressm_4_slot, stressm_2_slot, halo_info_mask, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)

         call ice_HaloUpdate_stress_compressed(stress12_1_slot, stress12_3_slot, halo_info_mask, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stress12_3_slot, stress12_1_slot, halo_info_mask, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stress12_2_slot, stress12_4_slot, halo_info_mask, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stress12_4_slot, stress12_2_slot, halo_info_mask, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)

         call ice_HaloDestroy(halo_info_mask)
      else
         call ice_HaloUpdate_stress_compressed(stressp_1_slot, stressp_3_slot, halo_info, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stressp_3_slot, stressp_1_slot, halo_info, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stressp_2_slot, stressp_4_slot, halo_info, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stressp_4_slot, stressp_2_slot, halo_info, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)

         call ice_HaloUpdate_stress_compressed(stressm_1_slot, stressm_3_slot, halo_info, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stressm_3_slot, stressm_1_slot, halo_info, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stressm_2_slot, stressm_4_slot, halo_info, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stressm_4_slot, stressm_2_slot, halo_info, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)

         call ice_HaloUpdate_stress_compressed(stress12_1_slot, stress12_3_slot, halo_info, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stress12_3_slot, stress12_1_slot, halo_info, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stress12_2_slot, stress12_4_slot, halo_info, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
         call ice_HaloUpdate_stress_compressed(stress12_4_slot, stress12_2_slot, halo_info, &
                                     field_loc_center,  field_type_scalar, &
                                     dynMap, dynOff)
      endif   ! maskhalo
      endif   ! tripole

      ! ⚠ No unpack here either.  The block stress arrays exist only for
      ! the restart writer, the history and evp1d's gather, and each of
      ! those now takes the slot copies or converts at its own boundary.

      !-----------------------------------------------------------------
      ! ice-ocean stress
      !-----------------------------------------------------------------

      call dyn_finish (1             , dynTot        , &
                       iceUmask_slot         , cdn_ocnU_slot , &
                       uvel_slot     , vvel_slot     , &
                       uocnU_slot    , vocnU_slot    , &
                       aiU_slot      , fmU_slot      , &
                       strocnxU      , strocnyU)

      if (grid_ice == 'CD' .or. grid_ice == 'C') then

         call dyn_finish (1             , dynTot        , &
                          iceNmask_slot         , cdn_ocnn_slot , &
                          uvelN_slot    , vvelN_slot    , &
                          uocnN_slot    , vocnN_slot    , &
                          aiN_slot      , fmN_slot      , &
                          strocnxN      , strocnyN)

         call dyn_finish (1             , dynTot        , &
                          iceEmask_slot         , cdn_ocne_slot , &
                          uvelE_slot    , vvelE_slot    , &
                          uocnE_slot    , vocnE_slot    , &
                          aiE_slot      , fmE_slot      , &
                          strocnxE      , strocnyE)

      endif

      ! strintxU and strintyU are diagnostics and nothing else: two halo
      ! exchanges and two averages whose only product is a history field, so
      ! they are skipped when no stream asks for one.
      if ((grid_ice == 'CD' .or. grid_ice == 'C') .and. lstrint) then
         call ice_timer_start(timer_bound)
         call ice_HaloUpdate (strintxE,        halo_info, &
                              field_loc_Eface, field_type_vector)
         call ice_HaloUpdate (strintyN,        halo_info, &
                              field_loc_Nface, field_type_vector)
         call ice_timer_stop(timer_bound)
         call grid_average_X2Y('S', strintxE, 'E', strintxU, 'U')    ! diagnostic
         call grid_average_X2Y('S', strintyN, 'N', strintyU, 'U')    ! diagnostic
      endif

      call ice_timer_stop(timer_dynamics)    ! dynamics

      end subroutine evp

!=======================================================================
! Computes the rates of strain and internal stress components for
! each of the four corners on each T-grid cell.
! Computes stress terms for the momentum equation
!
! author: Elizabeth C. Hunke, LANL

      subroutine stress (lo, hi,     icemask_s,     &
                         uvel,       vvel,       &
                         dxT,        dyT,        &
                         dxhy,       dyhx,       &
                         cxp,        cyp,        &
                         cxm,        cym,        &
                                     DminTarea,  &
                         strength,               &
                         stressp_1,  stressp_2,  &
                         stressp_3,  stressp_4,  &
                         stressm_1,  stressm_2,  &
                         stressm_3,  stressm_4,  &
                         stress12_1, stress12_2, &
                         stress12_3, stress12_4, &
                         nbrW, nbrSW, nbrS,                    &
                         str1, str2, str3, str4, &
                         str5, str6, str7, str8 )

      use ice_dyn_shared, only: strain_rates, visc_replpress

      integer (kind=int_kind), intent(in) :: &
         lo, hi                ! slot range to process

      ! uvel and vvel are COMPRESSED: this block's slice, 1-based, so cells
      ! are 1..cnt and the last element is the block's dead slot
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uvel     , & ! x-component of velocity (m/s)
         vvel         ! y-component of velocity (m/s)

      ! nbr(is,d) is the block-local slot one step in direction d from is
      integer (kind=int_kind), dimension (:), intent(in), contiguous :: &
         nbrW     , & ! neighbour slot one step W
         nbrSW    , & ! neighbour slot one step SW
         nbrS        ! neighbour slot one step S

      ! Every one of these is COMPRESSED and indexed by this block's slot, so
      ! nothing here needs nx_block or ny_block -- the kernel iterates over
      ! slots, as ice_dyn_core1d's stress_1d does.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         strength , & ! ice strength (N/m)
         dxT      , & ! width of T-cell through the middle (m)
         dyT      , & ! height of T-cell through the middle (m)
         dxhy     , & ! 0.5*(HTE - HTW)
         dyhx     , & ! 0.5*(HTN - HTS)
         cyp      , & ! 1.5*HTE - 0.5*HTW
         cxp      , & ! 1.5*HTN - 0.5*HTS
         cym      , & ! 0.5*HTE - 1.5*HTW
         cxm          ! 0.5*HTN - 1.5*HTS

      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         stressp_1, stressp_2, stressp_3, stressp_4 , & ! sigma11+sigma22
         stressm_1, stressm_2, stressm_3, stressm_4 , & ! sigma11-sigma22
         stress12_1,stress12_2,stress12_3,stress12_4    ! sigma12

      ! str1..str8 are COMPRESSED, indexed by BLOCK-LOCAL slot 1..cnt; the LAST
      ! element is the "no such neighbour" dead slot stepu falls back on
      real (kind=dbl_kind), dimension(:), intent(out), contiguous :: &
         str1, str2, str3, str4, str5, str6, str7, str8

      ! icemask_s(is) is .true. where the slot carries no ice this step; it
      ! replaces the icellT / indxTi / indxTj chain, and dyn_prep2 fills it
      ! over exactly that chain's range.
      logical (kind=log_kind), dimension(:), intent(in), contiguous :: &
         icemask_s

      ! local variables

      integer (kind=int_kind) :: &
         is             , & ! compressed slot
         isW, isS, isSW     ! its W, S and SW neighbour slots

      real (kind=dbl_kind) :: &
        divune, divunw, divuse, divusw            , & ! divergence
        tensionne, tensionnw, tensionse, tensionsw, & ! tension
        shearne, shearnw, shearse, shearsw        , & ! shearing
        Deltane, Deltanw, Deltase, Deltasw        , & ! Delt
        zetax2ne, zetax2nw, zetax2se, zetax2sw    , & ! 2 x zeta (bulk visc)
        etax2ne, etax2nw, etax2se, etax2sw        , & ! 2 x eta (shear visc)
        rep_prsne, rep_prsnw, rep_prsse, rep_prssw, & ! replacement pressure
!        puny                                      , & ! puny
        ssigpn, ssigps, ssigpe, ssigpw            , &
        ssigmn, ssigms, ssigme, ssigmw            , &
        ssig12n, ssig12s, ssig12e, ssig12w        , &
        ssigp1, ssigp2, ssigm1, ssigm2, ssig121, ssig122, &
        csigpne, csigpnw, csigpse, csigpsw        , &
        csigmne, csigmnw, csigmse, csigmsw        , &
        csig12ne, csig12nw, csig12se, csig12sw    , &
        str12ew, str12we, str12ns, str12sn        , &
        strp_tmp, strm_tmp

      ! DminTarea is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         DminTarea    ! deltaminEVP*tarea

      character(len=*), parameter :: subname = '(stress)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------


      ! schedule(static), not runtime: this walks uniform-cost slots, and CICE
      ! ships ICE_OMPSCHED="static,1", which round-robins individual cells and
      ! leaves every cache line holding str1..str8 write-shared.  Same reason
      ! ice_dyn_core1d pins its kernels.
      !$omp parallel do schedule(static)                           &
      !$omp default(none)                                          &
      !$omp private(is, isW, isS, isSW                           , &
      !$omp         divune, divunw, divuse, divusw               , &
      !$omp         tensionne, tensionnw, tensionse, tensionsw   , &
      !$omp         shearne, shearnw, shearse, shearsw           , &
      !$omp         Deltane, Deltanw, Deltase, Deltasw           , &
      !$omp         zetax2ne, zetax2nw, zetax2se, zetax2sw       , &
      !$omp         etax2ne, etax2nw, etax2se, etax2sw           , &
      !$omp         rep_prsne, rep_prsnw, rep_prsse, rep_prssw   , &
      !$omp         ssigpn, ssigps, ssigpe, ssigpw               , &
      !$omp         ssigmn, ssigms, ssigme, ssigmw               , &
      !$omp         ssig12n, ssig12s, ssig12e, ssig12w, ssigp1   , &
      !$omp         ssigp2, ssigm1, ssigm2, ssig121, ssig122     , &
      !$omp         csigpne, csigpnw, csigpse, csigpsw           , &
      !$omp         csigmne, csigmnw, csigmse, csigmsw           , &
      !$omp         csig12ne, csig12nw, csig12se, csig12sw       , &
      !$omp         str12ew, str12we, str12ns, str12sn           , &
      !$omp         strp_tmp, strm_tmp)                            &
      !$omp  shared(lo, hi, icemask_s, uvel, vvel, dxT, dyT         , &
      !$omp         dxhy, dyhx, cxp, cyp, cxm, cym, DminTarea    , &
      !$omp         strength, nbrW, nbrSW, nbrS                                , &
      !$omp         stressp_1, stressp_2, stressp_3, stressp_4   , &
      !$omp         stressm_1, stressm_2, stressm_3, stressm_4   , &
      !$omp         stress12_1, stress12_2, stress12_3, stress12_4, &
      !$omp         str1, str2, str3, str4, str5, str6, str7, str8, &
      !$omp         arlx1i, denom1, revp)
      do is = lo, hi
         if (.not. icemask_s(is)) cycle

          ! each neighbour slot was fetched twice, once for uvel and once for
          ! vvel; read them once instead
          isW  = nbrW(is)
          isS  = nbrS(is)
          isSW = nbrSW(is)

         !-----------------------------------------------------------------
         ! strain rates
         ! NOTE these are actually strain rates * area  (m^2/s)
         !-----------------------------------------------------------------

         call strain_rates (                            uvel(is),   uvel(isW),  &
                            uvel(isS),  uvel(isSW), &
                            vvel(is),   vvel(isW),  &
                            vvel(isS),  vvel(isSW), &
                            dxT(is),    dyT(is),    &
                            cxp(is),    cyp(is),    &
                            cxm(is),    cym(is),    &
                            divune,     divunw,     &
                            divuse,     divusw,     &
                            tensionne,  tensionnw,  &
                            tensionse,  tensionsw,  &
                            shearne,    shearnw,    &
                            shearse,    shearsw,    &
                            Deltane,    Deltanw,    &
                            Deltase,    Deltasw     )

         !-----------------------------------------------------------------
         ! viscosities and replacement pressure
         !-----------------------------------------------------------------

         call visc_replpress (strength(is), DminTarea(is), Deltane, &
                              zetax2ne, etax2ne, rep_prsne)

         call visc_replpress (strength(is), DminTarea(is), Deltanw, &
                              zetax2nw, etax2nw, rep_prsnw)

         call visc_replpress (strength(is), DminTarea(is), Deltasw, &
                              zetax2sw, etax2sw, rep_prssw)

         call visc_replpress (strength(is), DminTarea(is), Deltase, &
                              zetax2se, etax2se, rep_prsse)

         !-----------------------------------------------------------------
         ! the stresses                            ! kg/s^2
         ! (1) northeast, (2) northwest, (3) southwest, (4) southeast
         !-----------------------------------------------------------------

         ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code
         stressp_1 (is) = (stressp_1 (is)*(c1-arlx1i*revp) &
                           + arlx1i*(zetax2ne*divune - rep_prsne)) * denom1
         stressp_2 (is) = (stressp_2 (is)*(c1-arlx1i*revp) &
                           + arlx1i*(zetax2nw*divunw - rep_prsnw)) * denom1
         stressp_3 (is) = (stressp_3 (is)*(c1-arlx1i*revp) &
                           + arlx1i*(zetax2sw*divusw - rep_prssw)) * denom1
         stressp_4 (is) = (stressp_4 (is)*(c1-arlx1i*revp) &
                           + arlx1i*(zetax2se*divuse - rep_prsse)) * denom1

         stressm_1 (is) = (stressm_1 (is)*(c1-arlx1i*revp) &
                           + arlx1i*etax2ne*tensionne) * denom1
         stressm_2 (is) = (stressm_2 (is)*(c1-arlx1i*revp) &
                           + arlx1i*etax2nw*tensionnw) * denom1
         stressm_3 (is) = (stressm_3 (is)*(c1-arlx1i*revp) &
                           + arlx1i*etax2sw*tensionsw) * denom1
         stressm_4 (is) = (stressm_4 (is)*(c1-arlx1i*revp) &
                           + arlx1i*etax2se*tensionse) * denom1

         stress12_1(is) = (stress12_1(is)*(c1-arlx1i*revp) &
                           + arlx1i*p5*etax2ne*shearne) * denom1
         stress12_2(is) = (stress12_2(is)*(c1-arlx1i*revp) &
                           + arlx1i*p5*etax2nw*shearnw) * denom1
         stress12_3(is) = (stress12_3(is)*(c1-arlx1i*revp) &
                           + arlx1i*p5*etax2sw*shearsw) * denom1
         stress12_4(is) = (stress12_4(is)*(c1-arlx1i*revp) &
                           + arlx1i*p5*etax2se*shearse) * denom1

         !-----------------------------------------------------------------
         ! Eliminate underflows.
         ! The following code is commented out because it is relatively
         ! expensive and most compilers include a flag that accomplishes
         ! the same thing more efficiently.  This code is cheaper than
         ! handling underflows if the compiler lacks a flag; uncomment
         ! it in that case.  The compiler flag is often described with the
         ! phrase "flush to zero".
         !-----------------------------------------------------------------

!         call icepack_query_parameters(puny_out=puny)
!         call icepack_warnings_flush(nu_diag)
!         if (icepack_warnings_aborted()) call abort_ice(error_message=subname &
!            file=__FILE__, line=__LINE__)

!         stressp_1(is) = sign(max(abs(stressp_1(is)),puny),stressp_1(is))
!         stressp_2(is) = sign(max(abs(stressp_2(is)),puny),stressp_2(is))
!         stressp_3(is) = sign(max(abs(stressp_3(is)),puny),stressp_3(is))
!         stressp_4(is) = sign(max(abs(stressp_4(is)),puny),stressp_4(is))

!         stressm_1(is) = sign(max(abs(stressm_1(is)),puny),stressm_1(is))
!         stressm_2(is) = sign(max(abs(stressm_2(is)),puny),stressm_2(is))
!         stressm_3(is) = sign(max(abs(stressm_3(is)),puny),stressm_3(is))
!         stressm_4(is) = sign(max(abs(stressm_4(is)),puny),stressm_4(is))

!         stress12_1(is) = sign(max(abs(stress12_1(is)),puny),stress12_1(is))
!         stress12_2(is) = sign(max(abs(stress12_2(is)),puny),stress12_2(is))
!         stress12_3(is) = sign(max(abs(stress12_3(is)),puny),stress12_3(is))
!         stress12_4(is) = sign(max(abs(stress12_4(is)),puny),stress12_4(is))

         !-----------------------------------------------------------------
         ! combinations of the stresses for the momentum equation ! kg/s^2
         !-----------------------------------------------------------------

         ssigpn  = stressp_1(is) + stressp_2(is)
         ssigps  = stressp_3(is) + stressp_4(is)
         ssigpe  = stressp_1(is) + stressp_4(is)
         ssigpw  = stressp_2(is) + stressp_3(is)
         ssigp1  =(stressp_1(is) + stressp_3(is))*p055
         ssigp2  =(stressp_2(is) + stressp_4(is))*p055

         ssigmn  = stressm_1(is) + stressm_2(is)
         ssigms  = stressm_3(is) + stressm_4(is)
         ssigme  = stressm_1(is) + stressm_4(is)
         ssigmw  = stressm_2(is) + stressm_3(is)
         ssigm1  =(stressm_1(is) + stressm_3(is))*p055
         ssigm2  =(stressm_2(is) + stressm_4(is))*p055

         ssig12n = stress12_1(is) + stress12_2(is)
         ssig12s = stress12_3(is) + stress12_4(is)
         ssig12e = stress12_1(is) + stress12_4(is)
         ssig12w = stress12_2(is) + stress12_3(is)
         ssig121 =(stress12_1(is) + stress12_3(is))*p111
         ssig122 =(stress12_2(is) + stress12_4(is))*p111

         csigpne = p111*stressp_1(is) + ssigp2 + p027*stressp_3(is)
         csigpnw = p111*stressp_2(is) + ssigp1 + p027*stressp_4(is)
         csigpsw = p111*stressp_3(is) + ssigp2 + p027*stressp_1(is)
         csigpse = p111*stressp_4(is) + ssigp1 + p027*stressp_2(is)

         csigmne = p111*stressm_1(is) + ssigm2 + p027*stressm_3(is)
         csigmnw = p111*stressm_2(is) + ssigm1 + p027*stressm_4(is)
         csigmsw = p111*stressm_3(is) + ssigm2 + p027*stressm_1(is)
         csigmse = p111*stressm_4(is) + ssigm1 + p027*stressm_2(is)

         csig12ne = p222*stress12_1(is) + ssig122 &
                  + p055*stress12_3(is)
         csig12nw = p222*stress12_2(is) + ssig121 &
                  + p055*stress12_4(is)
         csig12sw = p222*stress12_3(is) + ssig122 &
                  + p055*stress12_1(is)
         csig12se = p222*stress12_4(is) + ssig121 &
                  + p055*stress12_2(is)

         str12ew = p5*dxT(is)*(p333*ssig12e + p166*ssig12w)
         str12we = p5*dxT(is)*(p333*ssig12w + p166*ssig12e)
         str12ns = p5*dyT(is)*(p333*ssig12n + p166*ssig12s)
         str12sn = p5*dyT(is)*(p333*ssig12s + p166*ssig12n)

         !-----------------------------------------------------------------
         ! for dF/dx (u momentum)
         !-----------------------------------------------------------------
         strp_tmp  = p25*dyT(is)*(p333*ssigpn  + p166*ssigps)
         strm_tmp  = p25*dyT(is)*(p333*ssigmn  + p166*ssigms)

         ! northeast (i,j)
         str1(is) = -strp_tmp - strm_tmp - str12ew &
              + dxhy(is)*(-csigpne + csigmne) + dyhx(is)*csig12ne

         ! northwest (i+1,j)
         str2(is) = strp_tmp + strm_tmp - str12we &
              + dxhy(is)*(-csigpnw + csigmnw) + dyhx(is)*csig12nw

         strp_tmp  = p25*dyT(is)*(p333*ssigps  + p166*ssigpn)
         strm_tmp  = p25*dyT(is)*(p333*ssigms  + p166*ssigmn)

         ! southeast (i,j+1)
         str3(is) = -strp_tmp - strm_tmp + str12ew &
              + dxhy(is)*(-csigpse + csigmse) + dyhx(is)*csig12se

         ! southwest (i+1,j+1)
         str4(is) = strp_tmp + strm_tmp + str12we &
              + dxhy(is)*(-csigpsw + csigmsw) + dyhx(is)*csig12sw

         !-----------------------------------------------------------------
         ! for dF/dy (v momentum)
         !-----------------------------------------------------------------
         strp_tmp  = p25*dxT(is)*(p333*ssigpe  + p166*ssigpw)
         strm_tmp  = p25*dxT(is)*(p333*ssigme  + p166*ssigmw)

         ! northeast (i,j)
         str5(is) = -strp_tmp + strm_tmp - str12ns &
              - dyhx(is)*(csigpne + csigmne) + dxhy(is)*csig12ne

         ! southeast (i,j+1)
         str6(is) = strp_tmp - strm_tmp - str12sn &
              - dyhx(is)*(csigpse + csigmse) + dxhy(is)*csig12se

         strp_tmp  = p25*dxT(is)*(p333*ssigpw  + p166*ssigpe)
         strm_tmp  = p25*dxT(is)*(p333*ssigmw  + p166*ssigme)

         ! northwest (i+1,j)
         str7(is) = -strp_tmp + strm_tmp + str12ns &
              - dyhx(is)*(csigpnw + csigmnw) + dxhy(is)*csig12nw

         ! southwest (i+1,j+1)
         str8(is) = strp_tmp - strm_tmp + str12sn &
              - dyhx(is)*(csigpsw + csigmsw) + dxhy(is)*csig12sw

      enddo                     ! is
      !$omp end parallel do

      end subroutine stress

!=======================================================================
! Computes the strain rates and internal stress components for C grid
!
! author: JF Lemieux, ECCC
! updated: D. Bailey, NCAR
! Nov 2021
!
! Bouillon, S., T. Fichefet, V. Legat and G. Madec (2013). The
! elastic-viscous-plastic method revisited. Ocean Model., 71, 2-12.
!
! Kimmritz, M., S. Danilov and M. Losch (2016). The adaptive EVP method
! for solving the sea ice momentum equation. Ocean Model., 101, 59-67.

      subroutine stressC_T  (lo, hi  ,   icemask_s   , &
                             nbrW, nbrSW, nbrS,              &
                             uvelE_s    , vvelN_s   , &
                             dxN_s      , dyE_s     , &
                             dxT_s        , dyT_s       , &
                             uarea_s    , uareaavg_s, &
                             DminTarea_s  ,             &
                             strength_s , shearU_s  , &
                             zetax2T_s  , etax2T_s  , &
                             stresspT_s , stressmT_s, &
                             stress12T_s)

      use ice_dyn_shared, only: visc_replpress, e_factor, &
                                nbr_W, nbr_S, nbr_SW

      integer (kind=int_kind), intent(in) :: &
         lo, hi                ! slot range to process

      ! .true. where the slot carries no ice this step
      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         icemask_s

      ! nbr(is,d) is the block-local slot one step in direction d from slot is
      integer (kind=int_kind), dimension (:), intent(in), contiguous :: &
         nbrW     , & ! neighbour slot one step W
         nbrSW    , & ! neighbour slot one step SW
         nbrS        ! neighbour slot one step S

      ! Compressed, block-local slots 1..cnt with the dead slot last.  shearU
      ! and uarea are read at W, S and SW here; the velocities and their
      ! metrics at W and S in the fused strain rates.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uvelE_s  , & ! x-component of velocity (m/s) at the E point
         vvelN_s  , & ! y-component of velocity (m/s) at the N point
         dxN_s    , & ! width of N-cell through the middle (m)
         dyE_s    , & ! height of E-cell through the middle (m)
         uarea_s  , & ! area of U-cell (m^2)
         uareaavg_s,& ! precomputed sum of the 4 U areas round this T cell
         shearU_s     ! shearU local for this routine

      ! Pointwise, so block-indexed.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         dxT_s      , & ! width of T-cell through the middle (m)
         dyT_s          ! height of T-cell through the middle (m)

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         strength_s   ! ice strength (N/m), compressed

      ! Written here, read one cell away by div_stress_*, so compressed.
      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         stresspT_s , & ! sigma11+sigma22
         stressmT_s , & ! sigma11-sigma22
         stress12T_s    ! sigma12

      ! Compressed because grid_average_X2YS_compressed shifts etax2T from T
      ! to U.  Nothing outside the dynamics reads zetax2T or etax2T -- checked
      ! across the tree -- so they have no block form at all.
      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         zetax2T_s , & ! zetax2 = 2*zeta (bulk viscosity)
         etax2T_s      ! etax2  = 2*eta  (shear viscosity)

      ! local variables

      integer (kind=int_kind) :: &
         is, isW, isS, isSW ! own, W, S and SW slots

      real (kind=dbl_kind) :: &
        divT      , & ! divergence at T point
        tensionT      ! tension at T point

      real (kind=dbl_kind) :: &
        shearTsqr , & ! strain rates squared at T point
        shearT    , & ! strain rate at T point
        DeltaT    , & ! delt at T point
        uareaavgr , & ! 1 / uarea avg
        rep_prsT      ! replacement pressure at T point

      ! DminTarea_s is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         DminTarea_s    ! deltaminEVP*tarea

      character(len=*), parameter :: subname = '(stressC_T)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      ! schedule(static), not runtime -- see the note in stress.
      !$omp parallel do schedule(static) default(none)              &
      !$omp private(is, isW, isS, isSW, divT, tensionT, shearT, shearTsqr, DeltaT, uareaavgr, rep_prsT) &
      !$omp  shared(lo, hi, icemask_s, nbrW, nbrSW, nbrS, uvelE_s, vvelN_s, dxN_s, &
      !$omp         dyE_s, dxT_s, dyT_s, uarea_s, uareaavg_s, &
      !$omp         DminTarea_s, strength_s, shearU_s, zetax2T_s, &
      !$omp         etax2T_s, stresspT_s, stressmT_s, stress12T_s, &
      !$omp         e_factor, arlx1i, revp, denom1)
      do is = lo, hi
         if (.not. icemask_s(is)) cycle

         !-----------------------------------------------------------------
         ! Square of shear strain rate at T obtained from interpolation of
         ! U point values (Bouillon et al., 2013, Kimmritz et al., 2016
         !-----------------------------------------------------------------
         isS  = nbrS(is)
         isSW = nbrSW(is)
         isW  = nbrW(is)

         ! divergence and tension, formerly a separate pass over every slot
         ! writing divT/tensionT.  Both are read only at this slot, so they
         ! are scalars here: two array writes and two reads per slot gone.
         divT     = dyE_s(is)*uvelE_s(is) - dyE_s(isW)*uvelE_s(isW) &
                  + dxN_s(is)*vvelN_s(is) - dxN_s(isS)*vvelN_s(isS)

         tensionT = (dyT_s(is)**2)*(uvelE_s(is)/dyE_s(is) - uvelE_s(isW)/dyE_s(isW)) &
                  - (dxT_s(is)**2)*(vvelN_s(is)/dxN_s(is) - vvelN_s(isS)/dxN_s(isS))

         ! The denominator is precomputed in BLOCK space -- see ice_grid.
         ! Summing it here would contribute 0 for an unmapped neighbour, where
         ! the block path uses that cell's real area; uarea is geometry, not
         ! state, and is nonzero even over land.  Numerator order is self, S,
         ! SW, W and must stay that way.
         uareaavgr = c1/uareaavg_s(is)

         shearTsqr = (shearU_s(is  )**2 * uarea_s(is  )  &
                    + shearU_s(isS )**2 * uarea_s(isS )  &
                    + shearU_s(isSW)**2 * uarea_s(isSW)  &
                    + shearU_s(isW )**2 * uarea_s(isW )) &
                    * uareaavgr

         shearT    = (shearU_s(is  )    * uarea_s(is  )  &
                    + shearU_s(isS )    * uarea_s(isS )  &
                    + shearU_s(isSW)    * uarea_s(isSW)  &
                    + shearU_s(isW )    * uarea_s(isW )) &
                    * uareaavgr

         DeltaT = sqrt(divT**2 + e_factor*(tensionT**2 + shearTsqr))

         !-----------------------------------------------------------------
         ! viscosities and replacement pressure at T point
         !-----------------------------------------------------------------

         call visc_replpress (strength_s(is), DminTarea_s(is), DeltaT, &
                              zetax2T_s (is), etax2T_s(is), rep_prsT)

         !-----------------------------------------------------------------
         ! the stresses                            ! kg/s^2
         !-----------------------------------------------------------------

         ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code

         stresspT_s(is)  = (stresspT_s (is)*(c1-arlx1i*revp) &
                          + arlx1i*(zetax2T_s(is)*divT - rep_prsT)) * denom1

         stressmT_s(is)  = (stressmT_s (is)*(c1-arlx1i*revp) &
                          + arlx1i*etax2T_s(is)*tensionT) * denom1

         stress12T_s(is) = (stress12T_s(is)*(c1-arlx1i*revp) &
                          + arlx1i*p5*etax2T_s(is)*shearT    ) * denom1

      enddo                     ! ij
      !$omp end parallel do

      end subroutine stressC_T

!=======================================================================
!
! Computes the strain rates and internal stress components for U points
!
! author: JF Lemieux, ECCC
! Nov 2021
!
! Bouillon, S., T. Fichefet, V. Legat and G. Madec (2013). The
! elastic-viscous-plastic method revisited. Ocean Model., 71, 2-12.
!
! Kimmritz, M., S. Danilov and M. Losch (2016). The adaptive EVP method
! for solving the sea ice momentum equation. Ocean Model., 101, 59-67.

      subroutine stressC_U  (lo       , hi        ,&
                             icemask_s    ,           &
                             uarea_s  ,            &
                             etax2U_s , deltaU    ,&
                             strengthU_s, shearU_s,&
                             stress12U_s)

      use ice_dyn_shared, only: visc_replpress, &
                                visc_method, deltaminEVP

      integer (kind=int_kind), intent(in) :: &
         lo, hi                ! slot range to process

      ! .true. where the slot carries no ice this step
      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         icemask_s

      ! Everything here is POINTWISE.  The compressed arrays appear only
      ! because other kernels read them at a neighbour -- shearU by stressC_T
      ! at W/S/SW, stress12U by div_stress_* -- so they live in slots and this
      ! routine has to address them that way.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         deltaU       ! deltaU array

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uarea_s    , & ! area of U point
         shearU_s   , & ! shearU array
         etax2U_s   , & ! 2*eta at the U point
         strengthU_s    ! ice strength at the U point

      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         stress12U_s  ! sigma12

      ! local variables

      integer (kind=int_kind) :: &
         is                 ! compressed slot of (i,j)

      real (kind=dbl_kind) :: &
         lzetax2U , & ! bulk viscosity at U point
         letax2U  , & ! shear viscosity at U point
         lrep_prsU, & ! replacement pressure at U point
         DminUarea    ! Dmin on U

      character(len=*), parameter :: subname = '(stressC_U)'

      !-----------------------------------------------------------------
      ! viscosities and replacement pressure at U point
      ! avg_zeta: Bouillon et al. 2013, C1 method of Kimmritz et al. 2016
      ! avg_strength: C2 method of Kimmritz et al. 2016
      ! if outside do and stress12U equation repeated in each loop for performance
      !-----------------------------------------------------------------

      if (visc_method == 'avg_zeta') then
         !$omp parallel do schedule(static) default(none) &
         !$omp private(is)                                &
         !$omp  shared(lo, hi, icemask_s, etax2U_s, shearU_s, &
         !$omp         stress12U_s, arlx1i, revp, denom1)
         do is = lo, hi
         if (.not. icemask_s(is)) cycle
            stress12U_s(is) = (stress12U_s(is)*(c1-arlx1i*revp) &
                             + arlx1i*p5*etax2U_s(is)*shearU_s(is)) * denom1
         enddo

      elseif (visc_method == 'avg_strength') then
         !$omp parallel do schedule(static) default(none)              &
         !$omp private(is, DminUarea, lzetax2U, letax2U, lrep_prsU)    &
         !$omp  shared(lo, hi, icemask_s, uarea_s, strengthU_s, deltaU,   &
         !$omp         shearU_s, stress12U_s, arlx1i, revp, denom1,    &
         !$omp         deltaminEVP)
         do is = lo, hi
         if (.not. icemask_s(is)) cycle
            DminUarea = deltaminEVP*uarea_s(is)
            ! only need etax2U here, but other terms are calculated with etax2U
            ! minimal extra calculations here even though it seems like there is
            call visc_replpress (strengthU_s(is), DminUarea, deltaU(is), &
                                 lzetax2U      , letax2U  , lrep_prsU)
            stress12U_s(is) = (stress12U_s(is)*(c1-arlx1i*revp) &
                             + arlx1i*p5*letax2U*shearU_s(is)) * denom1
         enddo

      endif

      end subroutine stressC_U

!=======================================================================
! Computes the strain rates and internal stress components for T points
!
! author: JF Lemieux, ECCC
! Nov 2021

      subroutine stressCD_T (lo, hi  ,   icemask_s   , &
                             nbrW, nbrS,              &
                             uvelE_s ,   vvelE_s  , &
                             uvelN_s ,   vvelN_s  , &
                             dxN_s   ,   dyE_s    , &
                             dxT_s     ,   dyT_s      , &
                                         DminTarea_s, &
                             strength_s,            &
                             zetax2T_s,  etax2T_s , &
                             stresspT_s, stressmT_s, &
                             stress12T_s)

      use ice_dyn_shared, only: visc_replpress, e_factor, &
                                nbr_W, nbr_S

      integer (kind=int_kind), intent(in) :: &
         lo, hi                ! slot range to process

      ! .true. where the slot carries no ice this step
      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         icemask_s

      ! nbr(is,d) is the block-local slot one step in direction d from slot is
      integer (kind=int_kind), dimension (:), intent(in), contiguous :: &
         nbrW     , & ! neighbour slot one step W
         nbrS        ! neighbour slot one step S

      ! Compressed: the fused strain rates reach W and S on these.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uvelE_s  , & ! x-component of velocity (m/s) at the E point
         vvelE_s  , & ! y-component of velocity (m/s) at the E point
         uvelN_s  , & ! x-component of velocity (m/s) at the N point
         vvelN_s  , & ! y-component of velocity (m/s) at the N point
         dxN_s    , & ! width of N-cell through the middle (m)
         dyE_s        ! height of E-cell through the middle (m)

      ! Pointwise, so block-indexed.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         dxT_s      , & ! width of T-cell through the middle (m)
         dyT_s          ! height of T-cell through the middle (m)

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         strength_s   ! ice strength (N/m), compressed

      ! Compressed because grid_average_X2YS_compressed shifts etax2T from T
      ! to U.  Nothing outside the dynamics reads zetax2T or etax2T -- checked
      ! across the tree -- so they have no block form at all.
      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         zetax2T_s , & ! zetax2 = 2*zeta (bulk viscosity)
         etax2T_s      ! etax2  = 2*eta  (shear viscosity)

      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         stresspT_s , & ! sigma11+sigma22
         stressmT_s , & ! sigma11-sigma22
         stress12T_s    ! sigma12

      ! local variables

      integer (kind=int_kind) :: &
         is                , & ! compressed slot of (i,j)
         isW, isS              ! its W and S neighbour slots

      real (kind=dbl_kind) :: &
        divT      , & ! divergence at T point
        tensionT  , & ! tension at T point
        shearT    , & ! shear at T point
        DeltaT        ! delt at T point

      real (kind=dbl_kind) :: &
        rep_prsT      ! replacement pressure at T point

      ! DminTarea_s is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         DminTarea_s    ! deltaminEVP*tarea

      character(len=*), parameter :: subname = '(stressCD_T)'

      !-----------------------------------------------------------------
      ! strain rates at T point
      ! NOTE these are actually strain rates * area  (m^2/s)
      !-----------------------------------------------------------------

      ! schedule(static), not runtime -- see the note in stress.
      !$omp parallel do schedule(static) default(none)              &
      !$omp private(is, isW, isS, divT, tensionT, shearT, DeltaT, rep_prsT) &
      !$omp  shared(lo, hi, icemask_s, nbrW, nbrS, uvelE_s, vvelE_s, uvelN_s, &
      !$omp         vvelN_s, dxN_s, dyE_s, dxT_s, dyT_s, DminTarea_s, &
      !$omp         strength_s, zetax2T_s, etax2T_s, stresspT_s, &
      !$omp         stressmT_s, stress12T_s, e_factor, arlx1i, revp, &
      !$omp         denom1)
      do is = lo, hi
         if (.not. icemask_s(is)) cycle

         ! strain rates at T, formerly two passes over every slot writing
         ! divT/tensionT/shearT/DeltaT.  All four are read only at this slot,
         ! so they are scalars: four array writes and six reads per slot gone,
         ! along with two whole traversals.
         isW = nbrW(is)
         isS = nbrS(is)

         divT     = dyE_s(is)*uvelE_s(is) - dyE_s(isW)*uvelE_s(isW) &
                  + dxN_s(is)*vvelN_s(is) - dxN_s(isS)*vvelN_s(isS)

         tensionT = (dyT_s(is)**2)*(uvelE_s(is)/dyE_s(is) - uvelE_s(isW)/dyE_s(isW)) &
                  - (dxT_s(is)**2)*(vvelN_s(is)/dxN_s(is) - vvelN_s(isS)/dxN_s(isS))

         shearT   = (dxT_s(is)**2)*(uvelN_s(is)/dxN_s(is) - uvelN_s(isS)/dxN_s(isS)) &
                  + (dyT_s(is)**2)*(vvelE_s(is)/dyE_s(is) - vvelE_s(isW)/dyE_s(isW))

         DeltaT   = sqrt(divT**2 + e_factor*(tensionT**2 + shearT**2))

         !-----------------------------------------------------------------
         ! viscosities and replacement pressure at T point
         !-----------------------------------------------------------------

         call visc_replpress (strength_s(is), DminTarea_s(is), DeltaT, &
                              zetax2T_s (is), etax2T_s(is), rep_prsT)

         !-----------------------------------------------------------------
         ! the stresses                            ! kg/s^2
         !-----------------------------------------------------------------

         ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code

         stresspT_s(is)  = (stresspT_s (is)*(c1-arlx1i*revp) &
                          + arlx1i*(zetax2T_s(is)*divT - rep_prsT)) * denom1

         stressmT_s(is)  = (stressmT_s (is)*(c1-arlx1i*revp) &
                          + arlx1i*etax2T_s(is)*tensionT) * denom1

         stress12T_s(is) = (stress12T_s(is)*(c1-arlx1i*revp) &
                          + arlx1i*p5*etax2T_s(is)*shearT) * denom1

      enddo                     ! ij
      !$omp end parallel do

      end subroutine stressCD_T

!=======================================================================
! Computes the strain rates and internal stress components for U points
!
! author: JF Lemieux, ECCC
! Nov 2021

      subroutine stressCD_U (lo       , hi        ,&
                             icemask_s    ,           &
                             uarea_s   ,           &
                             zetax2U_s , etax2U_s, &
                             strengthU_s,          &
                             divergU   , tensionU, &
                             shearU_s  , deltaU  , &
                             stresspU_s, stressmU_s, &
                             stress12U_s)

      use ice_dyn_shared, only: visc_replpress, &
                                visc_method, deltaminEVP

      integer (kind=int_kind), intent(in) :: &
         lo, hi                ! slot range to process

      ! .true. where the slot carries no ice this step
      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         icemask_s

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         zetax2U_s  , & ! 2*zeta at U point, compressed
         etax2U_s   , & ! 2*eta at U point, compressed
         strengthU_s    ! ice strength at U point, compressed

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         divergU  , & ! div at U point
         tensionU , & ! tension at U point
         deltaU       ! delt at U point

      ! Pointwise here; compressed because div_stress_* read them a cell away.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         uarea_s  , & ! area of U-cell (m^2)
         shearU_s     ! shear at U point

      real (kind=dbl_kind), dimension (:), intent(inout), contiguous :: &
         stresspU_s , & ! sigma11+sigma22
         stressmU_s , & ! sigma11-sigma22
         stress12U_s    ! sigma12

      ! local variables

      integer (kind=int_kind) :: &
         is                 ! compressed slot of (i,j)

      real (kind=dbl_kind) :: &
        lzetax2U  , & ! bulk viscosity at U point
        letax2U   , & ! shear viscosity at U point
        lrep_prsU , & ! replacement pressure at U point
        DminUarea     ! Dmin on U

      character(len=*), parameter :: subname = '(stressCD_U)'

      !$omp parallel do schedule(static) default(none)                 &
      !$omp private(is, DminUarea, lzetax2U, letax2U, lrep_prsU)       &
      !$omp  shared(lo, hi, icemask_s, uarea_s, zetax2U_s, etax2U_s,      &
      !$omp         strengthU_s, divergU, tensionU, shearU_s, deltaU,  &
      !$omp         stresspU_s, stressmU_s, stress12U_s,               &
      !$omp         arlx1i, revp, denom1, deltaminEVP, visc_method, ktens)
      do is = lo, hi
         if (.not. icemask_s(is)) cycle

         !-----------------------------------------------------------------
         ! viscosities and replacement pressure at U point
         ! avg_zeta: Bouillon et al. 2013, C1 method of Kimmritz et al. 2016
         ! avg_strength: C2 method of Kimmritz et al. 2016
         !-----------------------------------------------------------------

         if (visc_method == 'avg_zeta') then
            lzetax2U  = zetax2U_s(is)
            letax2U   = etax2U_s(is)
            lrep_prsU = (c1-Ktens)/(c1+Ktens)*lzetax2U*deltaU(is)

         elseif (visc_method == 'avg_strength') then
            DminUarea = deltaminEVP*uarea_s(is)
            ! only need etax2U here, but other terms are calculated with etax2U
            ! minimal extra calculations here even though it seems like there is
            call visc_replpress (strengthU_s(is), DminUarea, deltaU(is), &
                                 lzetax2U      , letax2U  , lrep_prsU )
         endif

         !-----------------------------------------------------------------
         ! the stresses                            ! kg/s^2
         !-----------------------------------------------------------------

         ! NOTE: for comp. efficiency 2 x zeta and 2 x eta are used in the code

         stresspU_s(is)  = (stresspU_s (is)*(c1-arlx1i*revp) &
                          + arlx1i*(lzetax2U*divergU(is) - lrep_prsU)) * denom1

         stressmU_s(is)  = (stressmU_s (is)*(c1-arlx1i*revp) &
                          + arlx1i*letax2U*tensionU(is)) * denom1

         stress12U_s(is) = (stress12U_s(is)*(c1-arlx1i*revp) &
                          + arlx1i*p5*letax2U*shearU_s(is)) * denom1

      enddo                     ! ij

      end subroutine stressCD_U

!=======================================================================
! Computes divergence of stress tensor at the E or N point for the mom equation
!
! author: JF Lemieux, ECCC
! Nov 2021
!
! Hunke, E. C., and J. K. Dukowicz (2002).  The Elastic-Viscous-Plastic
! Sea Ice Dynamics Model in General Orthogonal Curvilinear Coordinates
! on a Sphere - Incorporation of Metric Terms. Mon. Weather Rev.,
! 130, 1848-1865.
!
! Bouillon, S., M. Morales Maqueda, V. Legat and T. Fichefet (2009). An
! elastic-viscous-plastic sea ice model formulated on Arakawa B and C grids.
! Ocean Model., 27, 174-184.

      subroutine div_stress_Ex(lo, hi  , icemask_s   , &
                               nbrS, nbrE,            &
                               dxE_s     , dyE_s      , &
                               dxU_s   , dyT_s    , &
                               arear_s   , rheofactE_s, &
                               stressp_s, stressm_s, &
                               stress12_s,          &
                               strintx_s )

      use ice_dyn_shared, only: nbr_E, nbr_S

      integer (kind=int_kind), intent(in) :: &
         lo, hi                ! slot range to process

      ! .true. where the slot carries no ice this step
      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         icemask_s

      ! nbr(is,d) is the block-local slot one step in direction d from slot is
      integer (kind=int_kind), dimension (:), intent(in), contiguous :: &
         nbrS     , & ! neighbour slot one step S
         nbrE        ! neighbour slot one step E

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         dxE_s     , & ! width of E-cell through the middle (m)
         dyE_s     , & ! height of E-cell through the middle (m)
         rheofactE_s   ! mult. factor = 1, set to 0 if aiE <= rheo_area_min

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         arear_s       ! earear or narear

      ! Read at a neighbour -- dyT to the E, dxU to the S -- so compressed.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         dxU_s   , & ! width of U-cell through the middle (m), compressed
         dyT_s       ! height of T-cell through the middle (m), compressed

      real (kind=dbl_kind), optional, dimension (:), intent(in), contiguous :: &
         stressp_s , & ! stressp  (U or T) used for strintx calculation
         stressm_s , & ! stressm  (U or T) used for strintx calculation
         stress12_s    ! stress12 (U or T) used for strintx calculation

      ! Output only: the step kernel reads it in the same pass and nothing
      ! else does, so it stays compressed.  The block array that history
      ! reads is written once, after the subcycle.
      real (kind=dbl_kind), optional, dimension (:), intent(out), contiguous :: &
         strintx_s   ! div of stress tensor for u component

      ! local variables

      integer (kind=int_kind) :: &
         is, isA, isB

      character(len=*), parameter :: subname = '(div_stress_Ex)'

      ! schedule(static), not runtime -- see the note in stress.
      !$omp parallel do schedule(static) default(none)              &
      !$omp private(is, isA, isB)                                   &
      !$omp  shared(lo, hi, icemask_s, nbrS, nbrE, dxE_s, dyE_s, dxU_s, dyT_s, &
      !$omp         arear_s, rheofactE_s, stressp_s, stressm_s, &
      !$omp         stress12_s, strintx_s)
      do is = lo, hi
         if (.not. icemask_s(is)) cycle
         isA = nbrE(is)   ! (i+1,j  )
         isB = nbrS(is)   ! (i  ,j-1)
         strintx_s(is) = rheofactE_s(is) * arear_s(is) * &
              ( p5 * dyE_s(is)  * ( stressp_s(isA)  - stressp_s (is ) ) &
              + (p5/ dyE_s(is)) * ( (dyT_s(isA)**2) * stressm_s (isA)   &
                                  -(dyT_s(is )**2) * stressm_s (is ) ) &
              + (c1/ dxE_s(is)) * ( (dxU_s(is )**2) * stress12_s(is )   &
                                  -(dxU_s(isB)**2) * stress12_s(isB) ) )
      enddo
      !$omp end parallel do

      end subroutine div_stress_Ex

!=======================================================================

      subroutine div_stress_Ey(lo, hi  , icemask_s   , &
                               nbrS, nbrE,            &
                               dxE_s     , dyE_s      , &
                               dxU_s   , dyT_s    , &
                               arear_s   , rheofactE_s, &
                               stressp_s, stressm_s, &
                               stress12_s,          &
                               strinty_s )

      use ice_dyn_shared, only: nbr_S, nbr_E

      integer (kind=int_kind), intent(in) :: &
         lo, hi                ! slot range to process

      ! .true. where the slot carries no ice this step
      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         icemask_s

      ! nbr(is,d) is the block-local slot one step in direction d from slot is
      integer (kind=int_kind), dimension (:), intent(in), contiguous :: &
         nbrS     , & ! neighbour slot one step S
         nbrE        ! neighbour slot one step E

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         dxE_s     , & ! width of E-cell through the middle (m)
         dyE_s     , & ! height of E-cell through the middle (m)
         rheofactE_s   ! mult. factor = 1, set to 0 if aiE <= rheo_area_min

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         arear_s       ! earear or narear

      ! Read at a neighbour -- dxU to the S, dyT to the E -- so compressed.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         dxU_s   , & ! width of U-cell through the middle (m), compressed
         dyT_s       ! height of T-cell through the middle (m), compressed

      real (kind=dbl_kind), optional, dimension (:), intent(in), contiguous :: &
         stressp_s , & ! stressp  (U or T) used for strinty calculation
         stressm_s , & ! stressm  (U or T) used for strinty calculation
         stress12_s    ! stress12 (U or T) used for strinty calculation

      ! Output only: the step kernel reads it in the same pass and nothing
      ! else does, so it stays compressed.  The block array that history
      ! reads is written once, after the subcycle.
      real (kind=dbl_kind), optional, dimension (:), intent(out), contiguous :: &
         strinty_s   ! div of stress tensor for v component

      ! local variables

      integer (kind=int_kind) :: &
         is, isA, isB

      character(len=*), parameter :: subname = '(div_stress_Ey)'

      ! schedule(static), not runtime -- see the note in stress.
      !$omp parallel do schedule(static) default(none)              &
      !$omp private(is, isA, isB)                                   &
      !$omp  shared(lo, hi, icemask_s, nbrS, nbrE, dxE_s, dyE_s, dxU_s, dyT_s, &
      !$omp         arear_s, rheofactE_s, stressp_s, stressm_s, &
      !$omp         stress12_s, strinty_s)
      do is = lo, hi
         if (.not. icemask_s(is)) cycle
         isA = nbrS(is)   ! (i  ,j-1)
         isB = nbrE(is)   ! (i+1,j  )
         strinty_s(is) = rheofactE_s(is) * arear_s(is) * &
              ( p5 * dxE_s(is)  * ( stressp_s(is )  - stressp_s (isA) ) &
              - (p5/ dxE_s(is)) * ( (dxU_s(is )**2) * stressm_s (is )   &
                                  -(dxU_s(isA)**2) * stressm_s (isA) ) &
              + (c1/ dyE_s(is)) * ( (dyT_s(isB)**2) * stress12_s(isB)   &
                                  -(dyT_s(is )**2) * stress12_s(is ) ) )
      enddo
      !$omp end parallel do

      end subroutine div_stress_Ey

!=======================================================================

      subroutine div_stress_Nx(lo, hi  , icemask_s   , &
                               nbrW, nbrN,            &
                               dxN_s     , dyN_s      , &
                               dxT_s   , dyU_s    , &
                               arear_s   , rheofactN_s, &
                               stressp_s, stressm_s, &
                               stress12_s,          &
                               strintx_s )

      use ice_dyn_shared, only: nbr_W, nbr_N

      integer (kind=int_kind), intent(in) :: &
         lo, hi                ! slot range to process

      ! .true. where the slot carries no ice this step
      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         icemask_s

      ! nbr(is,d) is the block-local slot one step in direction d from slot is
      integer (kind=int_kind), dimension (:), intent(in), contiguous :: &
         nbrW     , & ! neighbour slot one step W
         nbrN        ! neighbour slot one step N

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         dxN_s     , & ! width of N-cell through the middle (m)
         dyN_s     , & ! height of N-cell through the middle (m)
         rheofactN_s   ! mult. factor = 1, set to 0 if aiN <= rheo_area_min

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         arear_s       ! earear or narear

      ! Read at a neighbour -- dyU to the W, dxT to the N -- so compressed.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         dxT_s   , & ! width of T-cell through the middle (m), compressed
         dyU_s       ! height of U-cell through the middle (m), compressed

      real (kind=dbl_kind), optional, dimension (:), intent(in), contiguous :: &
         stressp_s , & ! stressp  (U or T) used for strintx calculation
         stressm_s , & ! stressm  (U or T) used for strintx calculation
         stress12_s    ! stress12 (U or T) used for strintx calculation

      ! Output only: the step kernel reads it in the same pass and nothing
      ! else does, so it stays compressed.  The block array that history
      ! reads is written once, after the subcycle.
      real (kind=dbl_kind), optional, dimension (:), intent(out), contiguous :: &
         strintx_s   ! div of stress tensor for u component

      ! local variables

      integer (kind=int_kind) :: &
         is, isA, isB

      character(len=*), parameter :: subname = '(div_stress_Nx)'

      ! schedule(static), not runtime -- see the note in stress.
      !$omp parallel do schedule(static) default(none)              &
      !$omp private(is, isA, isB)                                   &
      !$omp  shared(lo, hi, icemask_s, nbrW, nbrN, dxN_s, dyN_s, dxT_s, dyU_s, &
      !$omp         arear_s, rheofactN_s, stressp_s, stressm_s, &
      !$omp         stress12_s, strintx_s)
      do is = lo, hi
         if (.not. icemask_s(is)) cycle
         isA = nbrW(is)   ! (i-1,j  )
         isB = nbrN(is)   ! (i  ,j+1)
         strintx_s(is) = rheofactN_s(is) * arear_s(is) * &
              ( p5 * dyN_s(is)  * ( stressp_s(is )  - stressp_s (isA) ) &
              + (p5/ dyN_s(is)) * ( (dyU_s(is )**2) * stressm_s (is )   &
                                  -(dyU_s(isA)**2) * stressm_s (isA) ) &
              + (c1/ dxN_s(is)) * ( (dxT_s(isB)**2) * stress12_s(isB)   &
                                  -(dxT_s(is )**2) * stress12_s(is ) ) )
      enddo
      !$omp end parallel do

      end subroutine div_stress_Nx

!=======================================================================

      subroutine div_stress_Ny(lo, hi  , icemask_s   , &
                               nbrW, nbrN,            &
                               dxN_s     , dyN_s      , &
                               dxT_s   , dyU_s    , &
                               arear_s   , rheofactN_s, &
                               stressp_s, stressm_s, &
                               stress12_s,          &
                               strinty_s )

      use ice_dyn_shared, only: nbr_N, nbr_W

      integer (kind=int_kind), intent(in) :: &
         lo, hi                ! slot range to process

      ! .true. where the slot carries no ice this step
      logical (kind=log_kind), dimension (:), intent(in), contiguous :: &
         icemask_s

      ! nbr(is,d) is the block-local slot one step in direction d from slot is
      integer (kind=int_kind), dimension (:), intent(in), contiguous :: &
         nbrW     , & ! neighbour slot one step W
         nbrN        ! neighbour slot one step N

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         dxN_s     , & ! width of N-cell through the middle (m)
         dyN_s     , & ! height of N-cell through the middle (m)
         rheofactN_s   ! mult. factor = 1, set to 0 if aiN <= rheo_area_min

      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         arear_s       ! earear or narear

      ! Read at a neighbour -- dxT to the N, dyU to the W -- so compressed.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous :: &
         dxT_s   , & ! width of T-cell through the middle (m), compressed
         dyU_s       ! height of U-cell through the middle (m), compressed

      real (kind=dbl_kind), optional, dimension (:), intent(in), contiguous :: &
         stressp_s , & ! stressp  (U or T) used for strinty calculation
         stressm_s , & ! stressm  (U or T) used for strinty calculation
         stress12_s    ! stress12 (U or T) used for strinty calculation

      ! Output only: the step kernel reads it in the same pass and nothing
      ! else does, so it stays compressed.  The block array that history
      ! reads is written once, after the subcycle.
      real (kind=dbl_kind), optional, dimension (:), intent(out), contiguous :: &
         strinty_s   ! div of stress tensor for v component

      ! local variables

      integer (kind=int_kind) :: &
         is, isA, isB

      character(len=*), parameter :: subname = '(div_stress_Ny)'

      ! schedule(static), not runtime -- see the note in stress.
      !$omp parallel do schedule(static) default(none)              &
      !$omp private(is, isA, isB)                                   &
      !$omp  shared(lo, hi, icemask_s, nbrW, nbrN, dxN_s, dyN_s, dxT_s, dyU_s, &
      !$omp         arear_s, rheofactN_s, stressp_s, stressm_s, &
      !$omp         stress12_s, strinty_s)
      do is = lo, hi
         if (.not. icemask_s(is)) cycle
         isA = nbrN(is)   ! (i  ,j+1)
         isB = nbrW(is)   ! (i-1,j  )
         strinty_s(is) = rheofactN_s(is) * arear_s(is) * &
              ( p5 * dxN_s(is)  * ( stressp_s(isA)  - stressp_s (is ) ) &
              - (p5/ dxN_s(is)) * ( (dxT_s(isA)**2) * stressm_s (isA)   &
                                  -(dxT_s(is )**2) * stressm_s (is ) ) &
              + (c1/ dyN_s(is)) * ( (dyU_s(is )**2) * stress12_s(is )   &
                                  -(dyU_s(isB)**2) * stress12_s(isB) ) )
      enddo
      !$omp end parallel do

      end subroutine div_stress_Ny

!=======================================================================

      end module ice_dyn_evp

!=======================================================================
