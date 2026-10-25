!=======================================================================
!
! Drivers for remapping and upwind ice transport
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! 2004: Revised by William Lipscomb from ice_transport_mpdata.
!       Stripped out mpdata, retained upwind, and added block structure.
! 2006: Incorporated remap transport driver and renamed from
!       ice_transport_upwind.
! 2011: ECH moved edgearea arrays into ice_transport_remap.F90

      module ice_transport_driver

      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, c1, p5, &
          field_loc_center, &
          field_type_scalar, field_type_vector, &
          field_loc_NEcorner, &
          field_loc_Nface, field_loc_Eface
      use ice_diagnostics, only: diagnostic_abort
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_compute_tracers
      use icepack_intfc, only: icepack_query_tracer_flags, &
          icepack_query_tracer_sizes, icepack_query_tracer_indices, &
          icepack_query_parameters

      implicit none
      private
      public :: init_transport, transport_remap, transport_upwind

      character (len=char_len), public ::     &
         advection   ! type of advection scheme used
                     ! 'upwind' => 1st order donor cell scheme
                     ! 'remap' => remapping scheme

! NOTE: For remapping, hice and hsno are considered tracers.
!       ntrace is not equal to ntrcr!

      integer (kind=int_kind) ::                      &
         ntrace              ! number of tracers in use

      integer (kind=int_kind), dimension(:), allocatable, public ::      &
         tracer_type     , & ! = 1, 2, or 3 (depends on 0, 1 or 2 other tracers)
         depend              ! tracer dependencies (see below)

      logical (kind=log_kind), dimension (:), allocatable, public ::     &
         has_dependents      ! true if a tracer has dependent tracers

      logical (kind=log_kind), public ::     &
         conserv_check       ! if true, check conservation

      integer (kind=int_kind), parameter ::                      &
         integral_order = 3  ! polynomial order of quadrature integrals
                             ! linear=1, quadratic=2, cubic=3

      logical (kind=log_kind), parameter ::     &
         l_dp_midpt = .true. ! if true, find departure points using
                             ! corrected midpoint velocity

!=======================================================================

      contains

!=======================================================================
!
! This subroutine is a wrapper for init_remap, which initializes the
! remapping transport scheme.  If the model is run with upwind
! transport, no initializations are necessary.
!
! authors William H. Lipscomb, LANL

      subroutine init_transport

      use ice_state, only: trcr_depend
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_advect
      use ice_transport_remap, only: init_remap

      integer (kind=int_kind) ::       &
         k, nt, nt1     ! tracer indices

      integer (kind=int_kind) :: &
         ntrcr    , nt_Tsfc  , nt_qice   , nt_qsno , &
         nt_sice  , nt_fbri  , nt_iage   , nt_FY   , &
         nt_alvl  , nt_vlvl  ,                       &
         nt_apnd  , nt_hpnd  , nt_ipnd   , nt_fsd  , &
         nt_smice , nt_smliq , nt_rhos   , nt_rsnw , &
         nt_isosno, nt_isoice, nt_bgc_Nit

      character(len=*), parameter :: subname = '(init_transport)'

      call ice_timer_start(timer_advect)  ! advection

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_qice_out=nt_qice, &
         nt_qsno_out=nt_qsno, nt_sice_out=nt_sice, nt_fbri_out=nt_fbri, &
         nt_iage_out=nt_iage, nt_FY_out=nt_FY, nt_fsd_out=nt_fsd, &
         nt_alvl_out=nt_alvl, nt_vlvl_out=nt_vlvl, &
         nt_apnd_out=nt_apnd, nt_hpnd_out=nt_hpnd, nt_ipnd_out=nt_ipnd, &
         nt_smice_out=nt_smice, nt_smliq_out=nt_smliq, nt_rhos_out=nt_rhos, &
         nt_rsnw_out=nt_rsnw, nt_bgc_Nit_out=nt_bgc_Nit, &
         nt_isosno_out=nt_isosno, nt_isoice_out=nt_isoice)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ntrace = 2 + ntrcr ! hice,hsno,trcr

      if (allocated(tracer_type))    deallocate(tracer_type)
      if (allocated(depend))         deallocate(depend)
      if (allocated(has_dependents)) deallocate(has_dependents)

      allocate (tracer_type   (ntrace), &
                depend        (ntrace), &
                has_dependents(ntrace))

      ! define tracer dependency arrays
      ! see comments in remapping routine

      depend(1:2)         = 0 ! hice, hsno
      tracer_type(1:2)    = 1 ! no dependency

      k = 2

      do nt = 1, ntrcr
         depend(k+nt) = trcr_depend(nt) ! 0 for ice area tracers
                                        ! 1 for ice volume tracers
                                        ! 2 for snow volume tracers
         tracer_type(k+nt) = 2          ! depends on 1 other tracer
         if (trcr_depend(nt) == 0) then
            tracer_type(k+nt) = 1       ! depends on no other tracers
         elseif (trcr_depend(nt) > 2) then
            if (trcr_depend(trcr_depend(nt)-2) > 0) then
               tracer_type(k+nt) = 3    ! depends on 2 other tracers
            endif
         endif
      enddo

      has_dependents = .false.
      do nt = 1, ntrace
         if (depend(nt) > 0) then
            nt1 = depend(nt)
            has_dependents(nt1) = .true.
            if (nt1 > nt) then
               write(nu_diag,*)     &
                  'Tracer nt2 =',nt,' depends on tracer nt1 =',nt1
               call abort_ice(subname//       &
                  'ERROR: remap transport: Must have nt2 > nt1')
            endif
         endif
      enddo                 ! ntrace

      ! diagnostic output
      if (my_task == master_task) then
         write (nu_diag, *) 'tracer          index  depend  type has_dependents'
            nt = 1
               write(nu_diag,1000) 'hi          ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            nt = 2
               write(nu_diag,1000) 'hs          ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
         k=2
         do nt = k+1, k+ntrcr
            if (nt-k==nt_Tsfc) &
               write(nu_diag,1000) 'nt_Tsfc     ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_qice) &
               write(nu_diag,1000) 'nt_qice     ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_qsno) &
               write(nu_diag,1000) 'nt_qsno     ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_sice) &
               write(nu_diag,1000) 'nt_sice     ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_fbri) &
               write(nu_diag,1000) 'nt_fbri     ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_iage) &
               write(nu_diag,1000) 'nt_iage     ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_FY) &
               write(nu_diag,1000) 'nt_FY       ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_alvl) &
               write(nu_diag,1000) 'nt_alvl     ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_vlvl) &
               write(nu_diag,1000) 'nt_vlvl     ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_apnd) &
               write(nu_diag,1000) 'nt_apnd     ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_hpnd) &
               write(nu_diag,1000) 'nt_hpnd     ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_ipnd) &
               write(nu_diag,1000) 'nt_ipnd     ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_smice) &
               write(nu_diag,1000) 'nt_smice    ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_smliq) &
               write(nu_diag,1000) 'nt_smliq    ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_rhos) &
               write(nu_diag,1000) 'nt_rhos     ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_rsnw) &
               write(nu_diag,1000) 'nt_rsnw     ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_fsd) &
               write(nu_diag,1000) 'nt_fsd      ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_isosno) &
               write(nu_diag,1000) 'nt_isosno   ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_isoice) &
               write(nu_diag,1000) 'nt_isoice   ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
            if (nt-k==nt_bgc_Nit) &
               write(nu_diag,1000) 'nt_bgc_Nit  ',nt,depend(nt),tracer_type(nt),&
                                                  has_dependents(nt)
         enddo
         write(nu_diag,*) ' '
      endif ! master_task
 1000 format (1x,a,2x,i6,2x,i6,2x,i4,4x,l4)

      if (trim(advection)=='remap') call init_remap ! grid quantities

      call ice_timer_stop(timer_advect)  ! advection

      end subroutine init_transport

!=======================================================================
!
! This subroutine solves the transport equations for one timestep
! using the conservative remapping scheme developed by John Dukowicz
! and John Baumgardner (DB) and modified for sea ice by William
! Lipscomb and Elizabeth Hunke.
!
! This scheme preserves monotonicity of ice area and tracers.  That is,
! it does not produce new extrema.  It is second-order accurate in space,
! except where gradients are limited to preserve monotonicity.
!
! authors William H. Lipscomb, LANL

      subroutine transport_remap (dt)

      use ice_blocks, only: nx_block, ny_block
      use ice_boundary, only: ice_HaloUpdate
      use ice_global_reductions, only: global_sum, global_sum_prod
      use ice_domain, only: nblocks, distrb_info, blocks_ice, halo_info
      use ice_domain_size, only: ncat, max_blocks
      use ice_blocks, only: nx_block, ny_block, block, get_block, nghost
      use ice_state, only: aice0, aice0_slot, aicen, vicen, vsnon, trcrn, &
          uvel, vvel, uvelE, vvelN
      use ice_bound_state, only: bound_state
      use ice_grid, only: tarea, grid_ice, &
          dynMap, dynOff, dynCnt, dynTot, unpack_from_slots
      use ice_calendar, only: istep1
      use ice_timers, only: ice_timer_start, ice_timer_stop, &
          timer_advect, timer_bound
      use ice_transport_remap, only: horizontal_remap, make_masks

      real (kind=dbl_kind), intent(in) ::     &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n              , & ! ice category index
         nt, nt1, nt2       ! tracer indices

      ! Compressed: (slot, category) and (slot, tracer, category) over all
      ! local blocks.  trm was the largest array in transport --
      ! nx_block*ny_block*ntrace*ncat*max_blocks -- and is now dynTot deep
      ! instead of the full block rectangle.
      real (kind=dbl_kind), dimension (dynTot,0:ncat) :: &
         aim          , & ! mean ice category areas in each grid cell
         aimask           ! = 1. if ice is present, = 0. otherwise

      real (kind=dbl_kind), dimension (dynTot,ntrace,ncat) :: &
         trm          , & ! mean tracer values in each grid cell
         trmask           ! = 1. if tracer is present, = 0. otherwise

      logical (kind=log_kind) :: &
         ckflag           ! if true, abort the model

      integer (kind=int_kind) :: &
         istop, jstop     ! indices of grid cell where model aborts

      integer (kind=int_kind) :: &
         ntrcr            ! number of tracers

      type (block) :: &
         this_block       ! block information for current block

      ! variables related to optional bug checks

      logical (kind=log_kind), parameter :: &
         l_monotonicity_check = .false.   ! if true, check monotonicity

      real (kind=dbl_kind), dimension(0:ncat) :: &
         asum_init    , & ! initial global ice area
         asum_final       ! final global ice area

      real (kind=dbl_kind), dimension(ntrace,ncat) :: &
         atsum_init   , & ! initial global ice area*tracer
         atsum_final      ! final global ice area*tracer

      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable :: &
         tmin         , & ! local min tracer
         tmax         , & ! local max tracer
         trm_b        , & ! trm expanded for the monotonicity check
         trmask_b         ! trmask expanded for the monotonicity check

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable :: &
         aim_b        , & ! aim expanded for the monotonicity check
         aimask_b         ! aimask expanded for the monotonicity check

      integer (kind=int_kind) :: &
         alloc_error

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1        , & ! expanded tracer-like field for the global sums
         work2            ! expanded area weight for the global sums

      character(len=char_len_long) :: &
         fieldid

      character(len=*), parameter :: subname = '(transport_remap)'

      call ice_timer_start(timer_advect)  ! advection
      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-------------------------------------------------------------------
      ! Prepare for remapping.
      ! Initialize, update ghost cells, fill tracer arrays.
      !-------------------------------------------------------------------

      ckflag = .false.
      istop = 0
      jstop = 0

      !-------------------------------------------------------------------
      ! Compute open water area in each grid cell.
      ! Note: An aggregate_area call is needed only if the open
      !       water area has changed since the previous call.
      !       Here we assume that aice0 is up to date.
      !-------------------------------------------------------------------

!      !$OMP PARALLEL DO PRIVATE(i,j,iblk) SCHEDULE(runtime)
!      do iblk = 1, nblocks
!      do j = 1, ny_block
!      do i = 1, nx_block
!         call aggregate_area (ncat,
!                              aicen(i,j,:,iblk),     &
!                              aice (i,j,  iblk),     &
!                              aice0(i,j,  iblk))
!      enddo
!      enddo
!      enddo
!      !$OMP END PARALLEL DO

      !-------------------------------------------------------------------
      ! Ghost cell updates for state variables.
      ! Commented out because ghost cells are updated after cleanup_itd.
      !-------------------------------------------------------------------
!      call ice_timer_start(timer_bound)

!      call ice_HaloUpdate (aice0,            halo_info,     &
!                           field_loc_center, field_type_scalar)

!      call bound_state (aicen,        &
!                        vicen, vsnon, &
!                        ntrcr, trcrn)

!      call ice_timer_stop(timer_bound)

      !-------------------------------------------------------------------
      ! Ghost cell updates for ice velocity.
      ! Commented out because ghost cell velocities are computed
      !  in ice_dyn_evp.
      !-------------------------------------------------------------------

!      call ice_timer_start(timer_bound)
!      call ice_HaloUpdate (uvel,               halo_info,     &
!                           field_loc_NEcorner, field_type_vector)
!      call ice_HaloUpdate (vvel,               halo_info,     &
!                           field_loc_NEcorner, field_type_vector)
!      call ice_timer_stop(timer_bound)


      !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
      do iblk = 1, nblocks

      !-------------------------------------------------------------------
      ! Fill arrays with fields to be remapped.
      !-------------------------------------------------------------------

         call state_to_tracers(nx_block,            ny_block,            &
                               dynMap(:,iblk),                           &
                               ntrcr,               ntrace,              &
                               aice0_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), aicen(:,:,:,  iblk), &
                               trcrn(:,:,:,:,iblk),                      &
                               vicen(:,:,:,  iblk), vsnon(:,:,:,  iblk), &
                               aim(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),  &
                               trm(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:,:))

      enddo
      !$OMP END PARALLEL DO

      !-------------------------------------------------------------------
      ! Optional conservation and monotonicity checks.
      !-------------------------------------------------------------------

      if (conserv_check) then

         !-------------------------------------------------------------------
         ! Compute initial values of globally conserved quantities.
         !-------------------------------------------------------------------

         ! ⚠ global_sum takes the block layout, so the compressed fields are
         ! expanded here.  Cells with no slot come back as c0 where the block
         ! form carried their old value -- land, where aicen is zero anyway
         ! and aice0 is untouched by transport.  Whatever they contribute is
         ! identical in the initial and final sums, so it cancels out of the
         ! comparison global_conservation makes, which is the only thing
         ! these numbers are used for.
         do n = 0, ncat
            call expand_slot_field(aim(:,n), work1)
            asum_init(n) = global_sum(work1,            distrb_info,       &
                                      field_loc_center, tarea)
         enddo

         do n = 1, ncat
            call expand_slot_field(aim(:,n), work2)
            do nt = 1, ntrace
               if (tracer_type(nt)==1) then ! does not depend on another tracer
                  call expand_slot_field(trm(:,nt,n), work1)
                  atsum_init(nt,n) =      &
                      global_sum_prod(work1(:,:,:), work2(:,:,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               elseif (tracer_type(nt)==2) then ! depends on another tracer
                  nt1 = depend(nt)
                  call expand_slot_field(trm(:,nt,n)*trm(:,nt1,n), work1)
                  atsum_init(nt,n) =     &
                      global_sum_prod(work1(:,:,:), work2(:,:,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               elseif (tracer_type(nt)==3) then ! depends on two tracers
                  nt1 = depend(nt)
                  nt2 = depend(nt1)
                  call expand_slot_field(trm(:,nt,n)*trm(:,nt1,n)          &
                                                    *trm(:,nt2,n), work1)
                  atsum_init(nt,n) =     &
                      global_sum_prod(work1(:,:,:), work2(:,:,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               endif            ! tracer_type
            enddo               ! nt
         enddo                  ! n

      endif                     ! conserv_check

      if (l_monotonicity_check) then

         ! This path is a debug tool -- l_monotonicity_check is a parameter,
         ! fixed .false., so the branch is eliminated at compile time.  It
         ! keeps the block layout: tmin/tmax go through the block
         ! ice_HaloUpdate, and converting them would buy nothing that ever
         ! runs.  aim and trm are expanded for it instead.
         allocate(tmin(nx_block,ny_block,ntrace,ncat,max_blocks),     &
                  tmax(nx_block,ny_block,ntrace,ncat,max_blocks),     &
                  aim_b(nx_block,ny_block,0:ncat,max_blocks),         &
                  trm_b(nx_block,ny_block,ntrace,ncat,max_blocks),    &
                  aimask_b(nx_block,ny_block,0:ncat,max_blocks),      &
                  trmask_b(nx_block,ny_block,ntrace,ncat,max_blocks), &
                  STAT=alloc_error)

         if (alloc_error /= 0)      &
              call abort_ice (subname//'ERROR: allocation error')

         tmin(:,:,:,:,:) = c0
         tmax(:,:,:,:,:) = c0

         !$OMP PARALLEL DO PRIVATE(iblk,n) SCHEDULE(runtime)
         do iblk = 1, nblocks

            !-------------------------------------------------------------------
            ! Compute masks.
            ! Masks are used to prevent tracer values in cells without ice
            !  from being used in the monotonicity check.
            !-------------------------------------------------------------------

            call make_masks (dynCnt(iblk),                                       &
                             ntrace,                                             &
                             has_dependents,                                     &
                             aim   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),   &
                             aimask(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),   &
                             trm   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:,:), &
                             trmask(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:,:))
         enddo
         !$OMP END PARALLEL DO

         do n = 0, ncat
            call expand_slot_field(aim   (:,n), aim_b   (:,:,n,:))
            call expand_slot_field(aimask(:,n), aimask_b(:,:,n,:))
         enddo
         do n = 1, ncat
            do nt = 1, ntrace
               call expand_slot_field(trm   (:,nt,n), trm_b   (:,:,nt,n,:))
               call expand_slot_field(trmask(:,nt,n), trmask_b(:,:,nt,n,:))
            enddo
         enddo

         !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,n) SCHEDULE(runtime)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            !-------------------------------------------------------------------
            ! Compute local max and min of tracer fields.
            !-------------------------------------------------------------------

            do n = 1, ncat
               call local_max_min                                        &
                            (nx_block,             ny_block,             &
                             ilo, ihi,             jlo, jhi,             &
                             trm_b (:,:,:,n,iblk),                       &
                             tmin  (:,:,:,n,iblk), tmax    (:,:,:,n,iblk), &
                             aimask_b(:,:,n,iblk), trmask_b(:,:,:,n,iblk))
            enddo
         enddo
         !$OMP END PARALLEL DO

         call ice_timer_start(timer_bound)
         call ice_HaloUpdate (tmin,             halo_info,     &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (tmax,             halo_info,     &
                              field_loc_center, field_type_scalar)
         call ice_timer_stop(timer_bound)

         !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,n) SCHEDULE(runtime)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do n = 1, ncat
               call quasilocal_max_min (nx_block, ny_block,  &
                                        ilo, ihi, jlo, jhi,  &
                                        tmin(:,:,:,n,iblk),  &
                                        tmax(:,:,:,n,iblk))
            enddo
         enddo
         !$OMP END PARALLEL DO

      endif                     ! l_monotonicity_check

      !-------------------------------------------------------------------
      ! Main remapping routine: Step ice area and tracers forward in time.
      !-------------------------------------------------------------------

      if (grid_ice == 'CD' .or. grid_ice == 'C') then
          call horizontal_remap (dt,             ntrace,         &
                                 uvel(:,:), vvel(:,:), &
                                 aim,            trm,            &
                                 tracer_type,    depend,         &
                                 has_dependents, integral_order, &
                                 l_dp_midpt,                     &
                                 uvelE  (:,:,:), vvelN  (:,:,:))
      else
          call horizontal_remap (dt,             ntrace,         &
                                 uvel(:,:), vvel(:,:), &
                                 aim,            trm,            &
                                 tracer_type,    depend,         &
                                 has_dependents, integral_order, &
                                 l_dp_midpt)
      endif

      !-------------------------------------------------------------------
      ! Given new fields, recompute state variables.
      !-------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
      do iblk = 1, nblocks

         call tracers_to_state (nx_block,            ny_block,            &
                                dynMap(:,iblk),                           &
                                ntrcr,               ntrace,              &
                                aim(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:),  &
                                trm(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:,:),&
                               aice0(:,:,    iblk), aice0_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), aicen(:,:,:,  iblk), &
                                trcrn(:,:,:,:,iblk),                      &
                                vicen(:,:,:,  iblk), vsnon(:,:,  :,iblk))

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      !-------------------------------------------------------------------
      ! Ghost cell updates for state variables.
      !-------------------------------------------------------------------

      call ice_timer_start(timer_bound)

      call bound_state (aicen,        &
                        vicen, vsnon, &
                        ntrcr, trcrn)

      call ice_timer_stop(timer_bound)

      !-------------------------------------------------------------------
      ! Optional conservation and monotonicity checks
      !-------------------------------------------------------------------

      !-------------------------------------------------------------------
      ! Compute final values of globally conserved quantities.
      ! Check global conservation of area and area*tracers.  (Optional)
      !-------------------------------------------------------------------

      if (conserv_check) then

         do n = 0, ncat
            call expand_slot_field(aim(:,n), work1)
            asum_final(n) = global_sum(work1,            distrb_info,      &
                                       field_loc_center, tarea)
         enddo

         do n = 1, ncat
            call expand_slot_field(aim(:,n), work2)
            do nt = 1, ntrace
               if (tracer_type(nt)==1) then ! does not depend on another tracer
                  call expand_slot_field(trm(:,nt,n), work1)
                  atsum_final(nt,n) =      &
                      global_sum_prod(work1(:,:,:), work2(:,:,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               elseif (tracer_type(nt)==2) then ! depends on another tracer
                  nt1 = depend(nt)
                  call expand_slot_field(trm(:,nt,n)*trm(:,nt1,n), work1)
                  atsum_final(nt,n) =     &
                      global_sum_prod(work1(:,:,:), work2(:,:,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               elseif (tracer_type(nt)==3) then ! depends on two tracers
                  nt1 = depend(nt)
                  nt2 = depend(nt1)
                  call expand_slot_field(trm(:,nt,n)*trm(:,nt1,n)          &
                                                    *trm(:,nt2,n), work1)
                  atsum_final(nt,n) =     &
                      global_sum_prod(work1(:,:,:), work2(:,:,:),          &
                                      distrb_info,  field_loc_center,      &
                                      tarea)
               endif            ! tracer_type
            enddo               ! nt
         enddo                  ! n

         if (my_task == master_task) then
            fieldid = subname//':000'
            call global_conservation (ckflag, fieldid,     &
                                      asum_init(0), asum_final(0))

            if (ckflag) then
               write (nu_diag,*) 'istep1, my_task =',     &
                                  istep1, my_task
               write (nu_diag,*) 'transport: conservation error, cat 0'
               call abort_ice(subname//'ERROR: conservation error1')
            endif

            do n = 1, ncat
               write(fieldid,'(a,i3.3)') subname,n
               call global_conservation                                 &
                                     (ckflag, fieldid,                  &
                                      asum_init(n),    asum_final(n),   &
                                      atsum_init(:,n), atsum_final(:,n))

               if (ckflag) then
                  write (nu_diag,*) 'istep1, my_task, cat =',     &
                                     istep1, my_task, n
                  write (nu_diag,*) 'transport: conservation error, cat ',n
                  call abort_ice(subname//'ERROR: conservation error2')
               endif
            enddo               ! n

         endif                  ! my_task = master_task

      endif                     ! conserv_check

      !-------------------------------------------------------------------
      ! Check tracer monotonicity.  (Optional)
      !-------------------------------------------------------------------

      if (l_monotonicity_check) then

         ! Re-expanded: aim and trm have been stepped forward by the remap
         ! since the pre-transport expansion above.
         do n = 0, ncat
            call expand_slot_field(aim(:,n), aim_b(:,:,n,:))
         enddo
         do n = 1, ncat
            do nt = 1, ntrace
               call expand_slot_field(trm(:,nt,n), trm_b(:,:,nt,n,:))
            enddo
         enddo

         !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,n,ckflag,istop,jstop) SCHEDULE(runtime)
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            ckflag = .false.
            istop = 0
            jstop = 0

            do n = 1, ncat
               call check_monotonicity (nx_block,           ny_block,           &
                                        ilo, ihi,           jlo, jhi,           &
                                        tmin(:,:,:,n,iblk), tmax(:,:,:,n,iblk), &
                                        aim_b(:,:, n,iblk), trm_b(:,:,:,n,iblk),&
                                        ckflag,                                 &
                                        istop,              jstop)

               if (ckflag) then
                  write (nu_diag,*) 'istep1, my_task, iblk, cat =',     &
                                     istep1, my_task, iblk, n
                  call diagnostic_abort(istop,jstop,iblk,' monotonicity error')
               endif
            enddo               ! n

         enddo                  ! iblk
         !$OMP END PARALLEL DO

         deallocate(tmin, tmax, aim_b, trm_b, aimask_b, trmask_b, STAT=alloc_error)
         if (alloc_error /= 0) call abort_ice (subname//'ERROR: deallocation error')

      endif                     ! l_monotonicity_check

      call ice_timer_stop(timer_advect)  ! advection

      end subroutine transport_remap

!=======================================================================
! Scatter a compressed field back to the block layout.
!
! Only the diagnostics need this: global_sum and the monotonicity check
! take the (i,j,block) shape.  Cells with no slot come back as c0, so the
! caller gets a field that is zero wherever the compressed layout does not
! reach -- see the note at the conservation sums for why that is harmless
! there.

      subroutine expand_slot_field (comp, work)

      use ice_domain, only: nblocks
      use ice_grid, only: dynMap, dynOff, dynCnt, unpack_from_slots

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         comp         ! compressed field, all local blocks

      real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
         work         ! block layout

      integer (kind=int_kind) :: iblk

      character(len=*), parameter :: subname = '(expand_slot_field)'

      work(:,:,:) = c0
      do iblk = 1, nblocks
         call unpack_from_slots(comp(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                                dynMap(:,iblk), work(:,:,iblk))
      enddo

      end subroutine expand_slot_field

!=======================================================================
!
! Computes the transport equations for one timestep using upwind. Sets
! several fields into a work array and passes it to upwind routine.

      subroutine transport_upwind (dt)

      use ice_boundary, only: ice_HaloUpdate, ice_HaloExtrapolate, &
          ice_HaloExtrapolate2DR8_compressed, ice_HaloUpdate2DR8_compressed
      use ice_blocks, only: nx_block, ny_block, block, get_block, nx_block, ny_block, &
          ew_boundary_type, ns_boundary_type
      use ice_domain, only: blocks_ice, halo_info, nblocks, distrb_info
      use ice_domain_size, only: ncat, max_blocks
      use ice_state, only: aice0, aice0_slot, aicen, vicen, vsnon, trcrn, &
          uvel, vvel, trcr_depend, trcr_base, &
          n_trcr_strata, nt_strata, uvelE, vvelN
      use ice_bound_state, only: bound_state
      use ice_flux, only: Tf
      use ice_grid, only: tmask, grid_ice, &
          dynMap, dynOff, dynCnt, dynNbr, dynTot, &
          skipEdgeE, skipEdgeN, skipPhys, &
          HTE_slot, HTN_slot, tarea_slot, pack_to_slots
      use ice_timers, only: ice_timer_start, ice_timer_stop, &
          timer_bound, timer_advect

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         ntrcr          , & !
         narr               ! max number of state variable arrays

      integer (kind=int_kind) :: &
         i, j, iblk     , & ! horizontal indices
         ms             , & ! slot of a cell
         nc             , & ! number of slots this block holds
         ilo,ihi,jlo,jhi    ! beginning and end of physical domain

      ! Compressed.  uee and vnn are the only per-timestep users of
      ! ice_HaloExtrapolate; every other caller is one-time grid setup in
      ! ice_grid on block-shaped metrics, which must stay block-shaped.
      real (kind=dbl_kind), dimension (dynTot) :: &
         uee, vnn           ! cell edge velocities

      ! Compressed: (slot, array).  narr is 1 + ncat*(3+ntrcr), so the block
      ! form was nx_block*ny_block*narr*max_blocks.
      real (kind=dbl_kind), dimension (:,:), allocatable :: &
         works              ! work array

      type (block) ::     &
         this_block         ! block information for current block

      character(len=*), parameter :: subname = '(transport_upwind)'

      call ice_timer_start(timer_advect)  ! advection

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      narr = 1 + ncat*(3+ntrcr) ! max number of state variable arrays

      allocate (works(dynTot,narr))

      !-------------------------------------------------------------------
      ! Get ghost cell values of state variables.
      ! (Assume velocities are already known for ghost cells, also.)
      !-------------------------------------------------------------------
!      call bound_state (aicen,        &
!                        vicen, vsnon, &
!                        ntrcr, trcrn)

!      call ice_timer_start(timer_bound)
!      call ice_HaloUpdate (uvel,               halo_info,     &
!                           field_loc_NEcorner, field_type_vector)
!      call ice_HaloUpdate (vvel,               halo_info,     &
!                           field_loc_NEcorner, field_type_vector)
!      call ice_timer_stop(timer_bound)

      !-------------------------------------------------------------------
      ! Average corner velocities to edges.
      !-------------------------------------------------------------------
      if (grid_ice == 'CD' .or. grid_ice == 'C') then
         !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
         do iblk = 1, nblocks
            call pack_to_slots(uvelE(:,:,iblk), dynMap(:,iblk), &
                               uee(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
            call pack_to_slots(vvelN(:,:,iblk), dynMap(:,iblk), &
                               vnn(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1))
         enddo
         !$OMP END PARALLEL DO
      else
         ! On C and CD the branch above is natural: uvelE and vvelN ARE
         ! E- and N-point velocities, defined over the whole array.  On B
         ! there is no such velocity -- it has to be INTERPOLATED from the
         ! corner values -- and that interpolation only ran over
         ! ilo..ihi, jlo..jhi, while ice_HaloExtrapolate below reads
         ! ilo-nghost..ihi+nghost.  The ghost columns it read had never been
         ! written.
         !
         ! ice_HaloUpdate overwrites most ghosts straight afterwards, so the
         ! only values that survived were at a true domain boundary, where
         ! there is no neighbour -- and those were extrapolated from
         ! whatever was in memory.  Reachable with advection='upwind',
         ! grid_ice='B' and a non-cyclic boundary, which is what
         ! set_nml.alt04 runs on gx3; invisible without -finit-real=snan.
         !
         ! So interpolate into the ghosts too.  The inputs are there: the
         ! commented-out halo calls above record that velocities are already
         ! known for ghost cells.  Only the outermost line each way cannot
         ! be done -- uee needs uvel one row south, vnn needs vvel one
         ! column west -- and those keep c0.
         uee(:) = c0
         vnn(:) = c0

         !$OMP PARALLEL DO PRIVATE(iblk,i,j,ms) SCHEDULE(runtime)
         do iblk = 1, nblocks

            do j = 2, ny_block          ! uee needs uvel one row south
            do i = 1, nx_block
               ms = dynMap((j-1)*nx_block + i, iblk)
               if (ms < 1) cycle
               uee(dynOff(iblk)+ms) = p5*(uvel((j-1)*nx_block+i,iblk) + uvel((j-1-1)*nx_block+i,iblk))
            enddo
            enddo

            do j = 1, ny_block
            do i = 2, nx_block          ! vnn needs vvel one column west
               ms = dynMap((j-1)*nx_block + i, iblk)
               if (ms < 1) cycle
               vnn(dynOff(iblk)+ms) = p5*(vvel((j-1)*nx_block+i,iblk) + vvel((j-1)*nx_block+i-1,iblk))
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO

         call ice_timer_start(timer_bound)
         call ice_HaloExtrapolate2DR8_compressed(uee, distrb_info, &
                                  ew_boundary_type, ns_boundary_type, &
                                  dynMap, dynOff)
         call ice_HaloExtrapolate2DR8_compressed(vnn, distrb_info, &
                                  ew_boundary_type, ns_boundary_type, &
                                  dynMap, dynOff)
         call ice_HaloUpdate2DR8_compressed (uee,             halo_info,     &
                              field_loc_Eface, field_type_vector, dynMap, dynOff)
         call ice_HaloUpdate2DR8_compressed (vnn,             halo_info,     &
                              field_loc_Nface, field_type_vector, dynMap, dynOff)
         call ice_timer_stop(timer_bound)
      endif

      ! ⚠ nc is per-iteration state and MUST be PRIVATE.  It is declared in
      ! this routine, not inside the region, so leaving it off the clause
      ! makes it shared -- which is exactly how 8da46d2 reached boreas.  No
      ! compiler warns.
      !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,nc) &
      !$OMP SCHEDULE(runtime)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         !-----------------------------------------------------------------
         ! fill work arrays with fields to be advected
         !-----------------------------------------------------------------

         nc = dynCnt(iblk)

         call state_to_work (nx_block,             ny_block,             &
                             dynMap(:,iblk),                             &
                             ntrcr,                                      &
                             narr,                 trcr_depend,          &
                             aicen (:,:,  :,iblk), trcrn (:,:,:,:,iblk), &
                             vicen (:,:,  :,iblk), vsnon (:,:,  :,iblk), &
                             aice0_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                             works (dynOff(iblk)+1:dynOff(iblk)+nc+1,:))

         !-----------------------------------------------------------------
         ! advect
         !-----------------------------------------------------------------

         call upwind_field (nc,                                            &
                            dynNbr  (dynOff(iblk)+1:dynOff(iblk)+nc,:),    &
                            skipEdgeE(dynOff(iblk)+1:dynOff(iblk)+nc),     &
                            skipEdgeN(dynOff(iblk)+1:dynOff(iblk)+nc),     &
                            skipPhys (dynOff(iblk)+1:dynOff(iblk)+nc),     &
                            dt,                                            &
                            narr,                                          &
                            works(dynOff(iblk)+1:dynOff(iblk)+nc+1,:),     &
                            uee(dynOff(iblk)+1:dynOff(iblk)+nc+1),         &
                            vnn(dynOff(iblk)+1:dynOff(iblk)+nc+1),         &
                            HTE_slot(dynOff(iblk)+1:dynOff(iblk)+nc+1),    &
                            HTN_slot(dynOff(iblk)+1:dynOff(iblk)+nc+1),    &
                            tarea_slot(dynOff(iblk)+1:dynOff(iblk)+nc+1))

         !-----------------------------------------------------------------
         ! convert work arrays back to state variables
         !-----------------------------------------------------------------

         call work_to_state (nx_block,            ny_block,             &
                             dynMap(:,iblk),                            &
                             ntrcr,               narr,                 &
                             trcr_depend(:),      trcr_base(:,:),       &
                             n_trcr_strata(:),    nt_strata(:,:),       &
                             tmask(:,      iblk), Tf    (:,:,iblk),     &
                             aicen(:,:,  :,iblk), trcrn (:,:,:,:,iblk), &
                             vicen(:,:,  :,iblk), vsnon (:,:,  :,iblk), &
                             aice0(:,:,    iblk), aice0_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                             works(dynOff(iblk)+1:dynOff(iblk)+nc+1,:))

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      deallocate (works)

      !-------------------------------------------------------------------
      ! Ghost cell updates for state variables.
      !-------------------------------------------------------------------

      call ice_timer_start(timer_bound)

      call bound_state (aicen,        &
                        vicen, vsnon, &
                        ntrcr, trcrn)

      call ice_timer_stop(timer_bound)

      call ice_timer_stop(timer_advect)  ! advection

      end subroutine transport_upwind

!=======================================================================
! The next few subroutines (through check_monotonicity) are called
! by transport_remap.
!=======================================================================
!
! Fill ice area and tracer arrays.
! Assume that the advected tracers are hicen, hsnon, trcrn,
!  qicen(1:nilyr), and qsnon(1:nslyr).
! This subroutine must be modified if a different set of tracers
!   is to be transported.  The rule for ordering tracers
!   is that a dependent tracer (such as qice) must have a larger
!   tracer index than the tracer it depends on (i.e., hice).
!
! author William H. Lipscomb, LANL

      subroutine state_to_tracers (nx_block, ny_block,   &
                                   map,                  &
                                   ntrcr,    ntrace,     &
                                   aice0,    aicen,      &
                                   trcrn,                &
                                   vicen,    vsnon,      &
                                   aim,      trm)

      use ice_domain_size, only: ncat, nslyr

      ! The state variables stay block-shaped: they are model state, read
      ! and written by I/O, and the compressed layout stops at this
      ! boundary.  aim and trm are the work arrays, so they are slots.
      !
      ! Every loop here ran over the WHOLE block, so no range mask is
      ! needed -- a cell with no slot is simply dropped, and nothing
      ! downstream reads one.

      integer (kind=int_kind), intent(in) ::     &
         nx_block, ny_block, & ! block dimensions
         ntrcr             , & ! number of tracers in use
         ntrace                ! number of tracers in use incl. hi, hs

      integer (kind=int_kind), dimension (:), intent(in) :: &
         map         ! this block's dynMap column

      ! ⚠ Compressed.  update_state writes aice0_slot in the same loop that
      ! calls icepack_aggregate, so it is current here, and this routine no
      ! longer walks the block to find it.  The slice is cnt+1 long: real
      ! cells then the dead slot, which matches aim's first dimension.
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aice0       ! fractional open water area, this block's slots

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), intent(in) :: &
         aicen   , & ! fractional ice area
         vicen   , & ! volume per unit area of ice          (m)
         vsnon       ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), intent(in) :: &
         trcrn       ! ice area tracers

      real (kind=dbl_kind), dimension (:,0:), intent(out) :: &
         aim         ! mean ice area in each grid cell

      real (kind=dbl_kind), dimension (:,:,:), intent(out) :: &
         trm         ! mean tracer values in each grid cell

      ! local variables

      integer (kind=int_kind) ::     &
         nt_qsno , & !
         i, j, n , & ! standard indices
         ms      , & ! slot of this cell
         it, kt      ! tracer indices

      real (kind=dbl_kind) ::     &
         puny    , & !
         rhos    , & ! snow density (km/m^3)
         Lfresh  , & ! latent heat of melting fresh ice (J/kg)
         w1          ! work variable

      character(len=*), parameter :: subname = '(state_to_tracers)'

      call icepack_query_parameters(puny_out=puny, rhos_out=rhos, &
           Lfresh_out=Lfresh)
      call icepack_query_tracer_indices(nt_qsno_out=nt_qsno)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! The dead slot is written too: it is the last element of the slice,
      ! and a stencil that lands on it must read c0.
      ! A straight copy now: slot order on both sides, dead slot included,
      ! so the block walk and its map lookup are gone.
      aim(:,0) = c0
      aim(1:size(aice0),0) = aice0(:)

      do n = 1, ncat

         trm(:,:,n) = c0
         aim(:,n)   = c0

         !-------------------------------------------------------------------
         ! Find grid cells where ice is present and fill area array.
         !-------------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            ms = map((j-1)*nx_block + i)
            if (ms < 1) cycle
            aim(ms,n) = aicen(i,j,n)
         enddo
         enddo

         !-------------------------------------------------------------------
         ! Fill tracer array
         ! Note: If aice > 0, then hice > 0, but we can have hsno = 0.
         ! Alse note: We transport qice*nilyr rather than qice, so as to
         !  avoid extra operations here and in tracers_to_state.
         !
         ! The icells list is gone: it selected aim > puny over the whole
         ! block, which is the test now made per cell as it is visited.
         !-------------------------------------------------------------------

         kt = 2

         do j = 1, ny_block
         do i = 1, nx_block
            ms = map((j-1)*nx_block + i)
            if (ms < 1) cycle
            if (aim(ms,n) <= puny) cycle

            w1 = c1 / aim(ms,n)
            trm(ms,1,n) = vicen(i,j,n) * w1 ! hice
            trm(ms,2,n) = vsnon(i,j,n) * w1 ! hsno

            do it = 1, ntrcr
               if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
                  trm(ms,kt+it,n) = trcrn(i,j,it,n) + rhos*Lfresh ! snow enthalpy
               else
                  trm(ms,kt+it,n) = trcrn(i,j,it,n) ! other tracers
               endif
            enddo
         enddo
         enddo
      enddo                     ! ncat

      end subroutine state_to_tracers

!=======================================================================
!
! Convert area and tracer arrays back to state variables.
!
! author William H. Lipscomb, LANL

      subroutine tracers_to_state (nx_block, ny_block,   &
                                   map,                  &
                                   ntrcr,    ntrace,     &
                                   aim,      trm,        &
                                   aice0,    aice0_s,  aicen,      &
                                   trcrn,                &
                                   vicen,    vsnon)

      use ice_domain_size, only: ncat, nslyr

      integer (kind=int_kind), intent(in) ::     &
         nx_block, ny_block, & ! block dimensions
         ntrcr             , & ! number of tracers in use
         ntrace                ! number of tracers in use incl. hi, hs

      integer (kind=int_kind), dimension (:), intent(in) :: &
         map       ! this block's dynMap column

      real (kind=dbl_kind), dimension (:,0:), intent(in) :: &
         aim       ! fractional ice area

      real (kind=dbl_kind), dimension (:,:,:), intent(in) :: &
         trm       ! mean tracer values in each grid cell

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         aice0     ! fractional ice area

      ! ⚠ BOTH layouts are written below.  The slot copy is what the
      ! transport boundary now reads on the way in; the block copy stays
      ! because history, restart and the thermo loops have not moved yet.
      ! Written in the SAME loop -- one extra store per cell, not a
      ! second traversal.
      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         aice0_s   ! the same, this block's slots

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), intent(inout) :: &
         aicen , & ! fractional ice area
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), intent(inout) :: &
         trcrn     ! tracers

      ! local variables

      integer (kind=int_kind) ::     &
         nt_qsno    , & !
         i, j, n    , & ! standard indices
         ms         , & ! slot of this cell
         it, kt         ! tracer indices

      real (kind=dbl_kind) :: &
         rhos       , & !
         Lfresh         !

      character(len=*), parameter :: subname = '(tracers_to_state)'

      call icepack_query_parameters(rhos_out=rhos, Lfresh_out=Lfresh)
      call icepack_query_tracer_indices(nt_qsno_out=nt_qsno)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      ! ⚠ The block form wrote aice0 over the whole block.  A cell with no
      ! slot keeps its old value instead -- which is the same value: it was
      ! copied into aim by state_to_tracers and update_fields cannot have
      ! changed it, because that cell is outside the physical domain or its
      ! fluxes resolved to the dead slot and carried c0.
      do j = 1, ny_block
      do i = 1, nx_block
         ms = map((j-1)*nx_block + i)
         if (ms < 1) cycle
         aice0(i,j) = aim(ms,0)
            aice0_s(ms)  = aim(ms,0)
      enddo
      enddo

      do n = 1, ncat

         !-------------------------------------------------------------------
         ! Compute state variables.
         !
         ! The icells list is gone: it selected aim > c0 over the whole
         ! block, which is the test now made per cell as it is visited.
         !-------------------------------------------------------------------

         kt = 2

         do j = 1, ny_block
         do i = 1, nx_block
            ms = map((j-1)*nx_block + i)
            if (ms < 1) cycle
            if (aim(ms,n) <= c0) cycle

            aicen(i,j,n) = aim(ms,n)
            vicen(i,j,n) = aim(ms,n)*trm(ms,1,n) ! aice*hice
            vsnon(i,j,n) = aim(ms,n)*trm(ms,2,n) ! aice*hsno

            do it = 1, ntrcr
               if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
                  trcrn(i,j,it,n) = trm(ms,kt+it,n) - rhos*Lfresh ! snow enthalpy
               else
                  trcrn(i,j,it,n) = trm(ms,kt+it,n)  ! other tracers
               endif
            enddo
         enddo
         enddo
      enddo                     ! ncat

      end subroutine tracers_to_state

!=======================================================================
!
! Check whether values of conserved quantities have changed.
! An error probably means that ghost cells are treated incorrectly.
!
! author William H. Lipscomb, LANL

      subroutine global_conservation (ckflag, fieldid,            &
                                      asum_init,  asum_final,     &
                                      atsum_init, atsum_final)

      character(len=*), intent(in) ::     &
         fieldid       ! field information string

      real (kind=dbl_kind), intent(in) ::     &
         asum_init , & ! initial global ice area
         asum_final    ! final global ice area

      real (kind=dbl_kind), dimension(ntrace), intent(in), optional :: &
         atsum_init, & ! initial global ice area*tracer
         atsum_final   ! final global ice area*tracer

      logical (kind=log_kind), intent(inout) ::     &
         ckflag        ! if true, abort on return

      ! local variables

      integer (kind=int_kind) ::     &
         nt            ! tracer index

      real (kind=dbl_kind) ::     &
         puny      , & !
         diff          ! difference between initial and final values

      character(len=*), parameter :: subname = '(global_conservation)'

      call icepack_query_parameters(puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (asum_init > puny) then
         diff = asum_final - asum_init
         if (abs(diff/asum_init) > puny) then
            ckflag = .true.
            write (nu_diag,*)
            write (nu_diag,*) subname,'Ice area conserv error ', trim(fieldid)
            write (nu_diag,*) subname,'  Initial global area  =', asum_init
            write (nu_diag,*) subname,'  Final global area    =', asum_final
            write (nu_diag,*) subname,'  Fractional error     =', abs(diff)/asum_init
            write (nu_diag,*) subname,'  asum_final-asum_init =', diff
         endif
      endif

      if (present(atsum_init)) then
         do nt = 1, ntrace
            if (abs(atsum_init(nt)) > puny) then
               diff = atsum_final(nt) - atsum_init(nt)
               if (abs(diff/atsum_init(nt)) > puny) then
                  ckflag = .true.
                  write (nu_diag,*)
                  write (nu_diag,*) subname,'Ice area*tracer conserv error ', trim(fieldid),nt
                  write (nu_diag,*) subname,'  Tracer index               =', nt
                  write (nu_diag,*) subname,'  Initial global area*tracer =', atsum_init(nt)
                  write (nu_diag,*) subname,'  Final global area*tracer   =', atsum_final(nt)
                  write (nu_diag,*) subname,'  Fractional error           =', abs(diff)/atsum_init(nt)
                  write (nu_diag,*) subname,'  atsum_final-atsum_init     =', diff
               endif
            endif
         enddo
      endif                     ! present(atsum_init)

      end subroutine global_conservation

!=======================================================================
!
! At each grid point, compute the local max and min of a scalar
! field phi: i.e., the max and min values in the nine-cell region
! consisting of the home cell and its eight neighbors.
!
! To extend to the neighbors of the neighbors (25 cells in all),
! follow this call with a call to quasilocal_max_min.
!
! author William H. Lipscomb, LANL

      subroutine local_max_min (nx_block, ny_block,     &
                                ilo, ihi, jlo, jhi,     &
                                trm,                    &
                                tmin,     tmax,         &
                                aimask,   trmask)

      integer (kind=int_kind), intent(in) ::     &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), intent(in), dimension(nx_block,ny_block) :: &
         aimask         ! ice area mask

      real (kind=dbl_kind), intent(in), dimension (nx_block,ny_block,ntrace) :: &
         trm        , & ! tracer fields
         trmask         ! tracer mask

      real (kind=dbl_kind), intent(out), dimension (nx_block,ny_block,ntrace) :: &
         tmin       , & ! local min tracer
         tmax           ! local max tracer

      ! local variables

      integer (kind=int_kind) ::     &
         i, j       , & ! horizontal indices
         nt, nt1        ! tracer indices

      real (kind=dbl_kind), dimension(nx_block,ny_block) ::     &
         phimask        ! aimask or trmask, as appropriate

      real (kind=dbl_kind) ::     &
         phi_nw, phi_n, phi_ne , & ! field values in 8 neighbor cells
         phi_w , phi_e         , &
         phi_sw, phi_s, phi_se

      character(len=*), parameter :: subname = '(local_max_min)'

      do nt = 1, ntrace

         if (tracer_type(nt)==1) then  ! does not depend on another tracer

            do j = 1, ny_block
            do i = 1, nx_block
               phimask(i,j) = aimask(i,j)
            enddo
            enddo

         else   ! depends on another tracer

            nt1 = depend(nt)
            do j = 1, ny_block
            do i = 1, nx_block
               phimask(i,j) = trmask(i,j,nt1)
            enddo
            enddo

         endif

         !-----------------------------------------------------------------------
         !  Store values of trm in the 8 neighbor cells.
         !  If aimask = 1, use the true value; otherwise use the home cell value
         !  so that non-physical values of phi do not contribute to the gradient.
         !-----------------------------------------------------------------------

         do j = jlo, jhi
         do i = ilo, ihi

            phi_nw = phimask(i-1,j+1) * trm(i-1,j+1,nt)     &
               + (c1-phimask(i-1,j+1))* trm(i,  j,  nt)
            phi_n  = phimask(i,  j+1) * trm(i,  j+1,nt)     &
               + (c1-phimask(i,  j+1))* trm(i,  j,  nt)
            phi_ne = phimask(i+1,j+1) * trm(i+1,j+1,nt)     &
               + (c1-phimask(i+1,j+1))* trm(i,  j,  nt)
            phi_w  = phimask(i-1,j)   * trm(i-1,j,  nt)     &
               + (c1-phimask(i-1,j))  * trm(i,  j,  nt)
            phi_e  = phimask(i+1,j)   * trm(i+1,j,  nt)     &
               + (c1-phimask(i+1,j))  * trm(i,  j,  nt)
            phi_sw = phimask(i-1,j-1) * trm(i-1,j-1,nt)     &
               + (c1-phimask(i-1,j-1))* trm(i,  j,  nt)
            phi_s  = phimask(i,  j-1) * trm(i,  j-1,nt)     &
               + (c1-phimask(i,  j-1))* trm(i,  j,  nt)
            phi_se = phimask(i+1,j-1) * trm(i+1,j-1,nt)     &
               + (c1-phimask(i+1,j-1))* trm(i,  j,  nt)

            !-----------------------------------------------------------------------
            !     Compute the minimum and maximum among the nine local cells.
            !-----------------------------------------------------------------------

            tmax(i,j,nt) = max (phi_nw, phi_n,  phi_ne, phi_w,     &
                   trm(i,j,nt), phi_e,  phi_sw, phi_s,  phi_se)

            tmin(i,j,nt) = min (phi_nw, phi_n,  phi_ne, phi_w,     &
                   trm(i,j,nt), phi_e,  phi_sw, phi_s,  phi_se)

         enddo         ! i
         enddo         ! j

      enddo            ! nt

      end subroutine local_max_min

!=======================================================================
!
! Extend the local max and min by one grid cell in each direction.
! Incremental remapping is monotone for the "quasilocal" max and min,
! but in rare cases may violate monotonicity for the local max and min.
!
! author William H. Lipscomb, LANL

      subroutine quasilocal_max_min (nx_block, ny_block,     &
                                     ilo, ihi, jlo, jhi,     &
                                     tmin,     tmax)

      integer (kind=int_kind), intent(in) ::     &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), intent(inout), dimension (nx_block,ny_block,ntrace) :: &
         tmin       , & ! local min tracer
         tmax           ! local max tracer

      ! local variables

      integer (kind=int_kind) ::     &
         i, j        , & ! horizontal indices
         nt              ! tracer index

      character(len=*), parameter :: subname = '(quasilocal_max_min)'

      do nt = 1, ntrace

         do j = jlo, jhi
         do i = ilo, ihi

            tmax(i,j,nt) =     &
              max (tmax(i-1,j+1,nt), tmax(i,j+1,nt), tmax(i+1,j+1,nt),     &
                   tmax(i-1,j,  nt), tmax(i,j,  nt), tmax(i+1,j,  nt),     &
                   tmax(i-1,j-1,nt), tmax(i,j-1,nt), tmax(i+1,j-1,nt))

            tmin(i,j,nt) =     &
              min (tmin(i-1,j+1,nt), tmin(i,j+1,nt), tmin(i+1,j+1,nt),     &
                   tmin(i-1,j,  nt), tmin(i,j,  nt), tmin(i+1,j,  nt),     &
                   tmin(i-1,j-1,nt), tmin(i,j-1,nt), tmin(i+1,j-1,nt))

         enddo                  ! i
         enddo                  ! j

      enddo

      end subroutine quasilocal_max_min

!======================================================================
!
! At each grid point, make sure that the new tracer values
! fall between the local max and min values before transport.
!
! author William H. Lipscomb, LANL

      subroutine check_monotonicity (nx_block, ny_block,     &
                                     ilo, ihi, jlo, jhi,     &
                                     tmin,     tmax,         &
                                     aim,      trm,          &
                                     ckflag,                 &
                                     istop,    jstop)

      integer (kind=int_kind), intent(in) ::     &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), intent(in), dimension (nx_block,ny_block) ::     &
         aim            ! new ice area

      real (kind=dbl_kind), intent(in), dimension (nx_block,ny_block,ntrace) ::     &
         trm            ! new tracers

      real (kind=dbl_kind), intent(in), dimension (nx_block,ny_block,ntrace) ::     &
         tmin       , & ! local min tracer
         tmax           ! local max tracer

      logical (kind=log_kind), intent(inout) ::     &
         ckflag         ! if true, abort on return

      integer (kind=int_kind), intent(inout) ::     &
         istop, jstop   ! indices of grid cell where model aborts

      ! local variables

      integer (kind=int_kind) ::     &
         i, j       , & ! horizontal indices
         nt, nt1, nt2   ! tracer indices

      real (kind=dbl_kind) ::     &
         puny       , & !
         w1, w2         ! work variables

      logical (kind=log_kind), dimension (nx_block, ny_block) ::   &
         l_check        ! if true, check monotonicity

      character(len=*), parameter :: subname = '(check_monotonicity)'

      call icepack_query_parameters(puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do nt = 1, ntrace

         !-------------------------------------------------------------------
         ! Load logical array to identify tracers that need checking.
         !-------------------------------------------------------------------

         if (tracer_type(nt)==1) then ! does not depend on another tracer

            do j = jlo, jhi
            do i = ilo, ihi
               if (aim(i,j) > puny) then
                  l_check(i,j) = .true.
               else
                  l_check(i,j) = .false.
               endif
            enddo
            enddo

         elseif (tracer_type(nt)==2) then ! depends on another tracer

            nt1 = depend(nt)
            do j = jlo, jhi
            do i = ilo, ihi
               if (abs(trm(i,j,nt1)) > puny) then
                  l_check(i,j) = .true.
               else
                  l_check(i,j) = .false.
               endif
            enddo
            enddo

         elseif (tracer_type(nt)==3) then ! depends on two tracers

            nt1 = depend(nt)
            nt2 = depend(nt1)
            do j = jlo, jhi
            do i = ilo, ihi
               if (abs(trm(i,j,nt1)) > puny .and.     &
                   abs(trm(i,j,nt2)) > puny) then
                  l_check(i,j) = .true.
               else
                  l_check(i,j) = .false.
               endif
            enddo
            enddo
         endif

         !-------------------------------------------------------------------
         ! Make sure new values lie between tmin and tmax
         !-------------------------------------------------------------------

         do j = jlo, jhi
         do i = ilo, ihi

            if (l_check(i,j)) then
               ! w1 and w2 allow for roundoff error when abs(trm) is big
               w1 = max(c1, abs(tmin(i,j,nt)))
               w2 = max(c1, abs(tmax(i,j,nt)))
               if (trm(i,j,nt) < tmin(i,j,nt)-w1*puny) then
                  ckflag = .true.
                  istop = i
                  jstop = j
                  write (nu_diag,*) ' '
                  write (nu_diag,*) 'new tracer < tmin'
                  write (nu_diag,*) 'i, j, nt =', i, j, nt
                  write (nu_diag,*) 'new tracer =', trm (i,j,nt)
                  write (nu_diag,*) 'tmin ='      , tmin(i,j,nt)
                  write (nu_diag,*) 'ice area ='  , aim(i,j)
               elseif (trm(i,j,nt) > tmax(i,j,nt)+w2*puny) then
                  ckflag = .true.
                  istop = i
                  jstop = j
                  write (nu_diag,*) ' '
                  write (nu_diag,*) 'new tracer > tmax'
                  write (nu_diag,*) 'i, j, nt =', i, j, nt
                  write (nu_diag,*) 'new tracer =', trm (i,j,nt)
                  write (nu_diag,*) 'tmax ='      , tmax(i,j,nt)
                  write (nu_diag,*) 'ice area ='  , aim(i,j)
               endif
            endif

         enddo                  ! i
         enddo                  ! j

      enddo                     ! nt

      end subroutine check_monotonicity

!=======================================================================
! The remaining subroutines are called by transport_upwind.
!=======================================================================
!
! Fill work array with state variables in preparation for upwind transport

      subroutine state_to_work (nx_block, ny_block,        &
                                map,                       &
                                ntrcr,                     &
                                narr,     trcr_depend,     &
                                aicen,    trcrn,           &
                                vicen,    vsnon,           &
                                aice0,    works)

      use ice_domain_size, only: ncat

      ! The state variables stay block-shaped -- they are model state, read
      ! and written by I/O.  works is the transport work array, so it is
      ! slots.  Every loop here runs over the whole block, so no range mask
      ! is needed: a cell with no slot is simply dropped.

      integer (kind=int_kind), intent(in) ::     &
         nx_block, ny_block, & ! block dimensions
         ntrcr             , & ! number of tracers in use
         narr                  ! number of 2D state variable arrays in works array

      integer (kind=int_kind), dimension (:), intent(in) :: &
         map         ! this block's dynMap column

      integer (kind=int_kind), dimension (ntrcr), intent(in) ::     &
         trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), intent(in) ::     &
         aicen   , & ! concentration of ice
         vicen   , & ! volume per unit area of ice          (m)
         vsnon       ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), intent(in) ::     &
         trcrn       ! ice tracers

      ! compressed, as in state_to_tracers
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aice0       ! concentration of open water, this block's slots

      real (kind=dbl_kind), dimension(:,:), intent (out) ::      &
         works       ! work array, (slot, array)

      ! local variables

      integer (kind=int_kind) :: &
         nt_alvl, nt_apnd, nt_fbri

      logical (kind=log_kind) :: &
         tr_pond_lvl, tr_pond_topo, tr_pond_sealvl

      integer (kind=int_kind) ::      &
         i, j, n, it, & ! counting indices
         ms         , & ! slot of this cell
         narrays        ! counter for number of state variable arrays

      character(len=*), parameter :: subname = '(state_to_work)'

      call icepack_query_tracer_flags(tr_pond_lvl_out=tr_pond_lvl, &
           tr_pond_topo_out=tr_pond_topo, tr_pond_sealvl_out=tr_pond_sealvl)
      call icepack_query_tracer_indices(nt_alvl_out=nt_alvl, nt_apnd_out=nt_apnd, &
           nt_fbri_out=nt_fbri)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-----------------------------------------------------------------
      ! This array is used for performance (balance memory/cache vs
      ! number of bound calls);  a different number of arrays may perform
      ! better depending on the machine used, number of processors, etc.
      ! --tested on SGI R2000, using 4 pes for the ice model under MPI
      !-----------------------------------------------------------------

      ! The dead slot is never reached by the traversal below, and
      ! upwind_field reads it wherever a neighbour is missing, so it has to
      ! start at c0.
      works(:,:) = c0

      ! straight copy: slot order on both sides
      works(1:size(aice0),1) = aice0(:)
      narrays = 1

      do n=1, ncat

         do j = 1, ny_block
         do i = 1, nx_block
            ms = map((j-1)*nx_block + i)
            if (ms < 1) cycle
            works(ms,narrays+1) = aicen(i,j,n)
            works(ms,narrays+2) = vicen(i,j,n)
            works(ms,narrays+3) = vsnon(i,j,n)
         enddo                  ! i
         enddo                  ! j
         narrays = narrays + 3

         do it = 1, ntrcr
            if (trcr_depend(it) == 0) then
               do j = 1, ny_block
               do i = 1, nx_block
                  ms = map((j-1)*nx_block + i)
                  if (ms < 1) cycle
                  works(ms,narrays+it) = aicen(i,j,n)*trcrn(i,j,it,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 1) then
               do j = 1, ny_block
               do i = 1, nx_block
                  ms = map((j-1)*nx_block + i)
                  if (ms < 1) cycle
                  works(ms,narrays+it) = vicen(i,j,n)*trcrn(i,j,it,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 2) then
               do j = 1, ny_block
               do i = 1, nx_block
                  ms = map((j-1)*nx_block + i)
                  if (ms < 1) cycle
                  works(ms,narrays+it) = vsnon(i,j,n)*trcrn(i,j,it,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 2+nt_alvl) then
               do j = 1, ny_block
               do i = 1, nx_block
                  ms = map((j-1)*nx_block + i)
                  if (ms < 1) cycle
                  works(ms,narrays+it) = aicen(i,j        ,n) &
                                        * trcrn(i,j,nt_alvl,n) &
                                        * trcrn(i,j,it     ,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 2+nt_apnd .and. &
                    (tr_pond_topo .or. tr_pond_sealvl)) then
               do j = 1, ny_block
               do i = 1, nx_block
                  ms = map((j-1)*nx_block + i)
                  if (ms < 1) cycle
                  works(ms,narrays+it) = aicen(i,j        ,n) &
                                        * trcrn(i,j,nt_apnd,n) &
                                        * trcrn(i,j,it     ,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 2+nt_apnd .and. &
                    tr_pond_lvl) then
               do j = 1, ny_block
               do i = 1, nx_block
                  ms = map((j-1)*nx_block + i)
                  if (ms < 1) cycle
                  works(ms,narrays+it) = aicen(i,j        ,n) &
                                        * trcrn(i,j,nt_alvl,n) &
                                        * trcrn(i,j,nt_apnd,n) &
                                        * trcrn(i,j,it     ,n)
               enddo
               enddo
            elseif (trcr_depend(it) == 2+nt_fbri) then
               do j = 1, ny_block
               do i = 1, nx_block
                  ms = map((j-1)*nx_block + i)
                  if (ms < 1) cycle
                  works(ms,narrays+it) = vicen(i,j        ,n) &
                                        * trcrn(i,j,nt_fbri,n) &
                                        * trcrn(i,j,it     ,n)
               enddo
               enddo
            endif
         enddo
         narrays = narrays + ntrcr

      enddo                     ! n

      if (narr /= narrays) write(nu_diag,*)      &
           "Wrong number of arrays in transport bound call"

      end subroutine state_to_work

!=======================================================================
!
! Convert work array back to state variables

      subroutine work_to_state (nx_block,      ny_block, &
                                map,                     &
                                ntrcr,         narr,     &
                                trcr_depend,             &
                                trcr_base,               &
                                n_trcr_strata,           &
                                nt_strata,               &
                                tmask,         Tf,       &
                                aicen,         trcrn,    &
                                vicen,         vsnon,    &
                                aice0,         aice0_s,  works)

      use ice_domain_size, only: ncat

      integer (kind=int_kind), intent (in) :: &
         nx_block, ny_block, & ! block dimensions
         ntrcr             , & ! number of tracers in use
         narr                  ! number of 2D state variable arrays in works array

      integer (kind=int_kind), dimension (:), intent(in) :: &
         map                   ! this block's dynMap column

      integer (kind=int_kind), dimension (ntrcr), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (ntrcr,3), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (ntrcr,2), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      ! tmask is assumed shape, 1D over the block in Fortran storage order,
      ! i.e. element (i,j) is stored at (j-1)*nx_block+i
      logical (kind=log_kind), dimension (:), intent (in) :: &
         tmask

      real (kind=dbl_kind), intent (in) :: &
         Tf    (nx_block,ny_block)

      ! Compressed: (slot, array).  A slot row is contiguous over the array
      ! index, which is what icepack_compute_tracers wants for atrcrn.
      real (kind=dbl_kind), dimension (:,:), intent (in) :: &
         works

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), intent(out) :: &
         aicen      , & ! concentration of ice
         vicen      , & ! volume per unit area of ice          (m)
         vsnon          ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat),intent(out) :: &
         trcrn          ! ice tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         aice0          ! concentration of open water

      ! ⚠ BOTH layouts are written below.  The slot copy is what the
      ! transport boundary now reads on the way in; the block copy stays
      ! because history, restart and the thermo loops have not moved yet.
      ! Written in the SAME loop -- one extra store per cell, not a
      ! second traversal.
      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         aice0_s   ! the same, this block's slots

      ! local variables

      ! The nx_block*ny_block x narr copy of works that used to live here is
      ! gone.  It existed only to hand icepack_compute_tracers a contiguous
      ! (cell, tracer) row; a slot row of works already is one, so the array
      ! and the identity indxi/indxj traversal that filled it both go.
      integer (kind=int_kind) ::  &
         i, j, n    , & ! counting indices
         ms         , & ! slot of this cell
         narrays    , & ! counter for number of state variable arrays
         nt_Tsfc        ! Tsfc tracer number

      character(len=*), parameter :: subname = '(work_to_state)'

      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      do j=1,ny_block
      do i=1,nx_block
         ms = map((j-1)*nx_block + i)
         if (ms < 1) cycle
         aice0(i,j) = works(ms,1)
            aice0_s(ms)  = works(ms,1)
      enddo
      enddo
      narrays = 1               ! aice0 is first array

      do n=1,ncat

         do j=1,ny_block
         do i=1,nx_block
            ms = map((j-1)*nx_block + i)
            if (ms < 1) cycle

            aicen(i,j,n) = works(ms,narrays+1)
            vicen(i,j,n) = works(ms,narrays+2)
            vsnon(i,j,n) = works(ms,narrays+3)

            call icepack_compute_tracers(trcr_depend   = trcr_depend(:),   &
                                         atrcrn        = works(ms,narrays+4:narrays+3+ntrcr), &
                                         aicen         = aicen(i,j,n),     &
                                         vicen         = vicen(i,j,n),     &
                                         vsnon         = vsnon(i,j,n),     &
                                         trcr_base     = trcr_base(:,:),   &
                                         n_trcr_strata = n_trcr_strata(:), &
                                         nt_strata     = nt_strata(:,:),   &
                                         trcrn         = trcrn(i,j,:,n),   &
                                         Tf            = Tf(i,j))

            ! tcraig, don't let land points get non-zero Tsfc
            if (.not.tmask((j-1)*nx_block+i)) then
               trcrn(i,j,nt_Tsfc,n) = c0
            endif

         enddo
         enddo

         narrays = narrays + 3 + ntrcr

      enddo                     ! ncat

      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      end subroutine work_to_state

!=======================================================================
!
! upwind transport algorithm

      subroutine upwind_field (ncells,   nbr,        &
                               skipE,    skipN,      &
                               skipme,               &
                               dt,                   &
                               narrays,  phi,        &
                               uee,      vnn,        &
                               HTE,      HTN,        &
                               tarea)

      use ice_grid, only: nbr_W, nbr_E, nbr_N, nbr_S

      integer (kind=int_kind), intent (in) ::     &
         ncells            , & ! number of slots this block holds
         narrays               ! number of 2D arrays to be transported

      integer (kind=int_kind), dimension (:,:), intent(in) ::   &
         nbr                   ! neighbour slot in each compass direction

      ! ⚠ Three ranges, three masks.  The block form computed worka and
      ! workb together over 1..ihi x 1..jhi, which is a SUPERSET of what it
      ! then read: worka is read at (i,j) and (i-1,j) over the physical
      ! domain, so it is needed on ilo-1..ihi x jlo..jhi, and workb at
      ! (i,j) and (i,j-1), so on ilo..ihi x jlo-1..jhi.  Those two ranges
      ! are exactly skipEdgeE and skipEdgeN, which is what they were built
      ! for, so the single loop splits in two and each computes only what is
      ! read.  The update itself is the physical domain, skipPhys.
      logical (kind=log_kind), dimension (:), intent(in) ::   &
         skipE             , & ! .true. outside the EAST-edge range
         skipN             , & ! .true. outside the NORTH-edge range
         skipme                ! .true. outside the physical domain

      real (kind=dbl_kind), intent(in) ::         &
         dt                    ! time step

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         phi                   ! scalar field, (slot, array)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         uee, vnn              ! cell edge velocities

      ! ⚠ all read at this cell only, never summed across neighbours
      real (kind=dbl_kind), dimension(:), intent(in) :: &
         HTE               , & ! length of east cell edge
         HTN               , & ! length of north cell edge
         tarea                 ! grid cell area

      ! local variables

      integer (kind=int_kind) :: &
         is, n             , & ! standard indices
         isE, isN, isW, isS    ! neighbour slots

      real (kind=dbl_kind), dimension (size(phi,1)) :: &
         worka, workb

      character(len=*), parameter :: subname = '(upwind_field)'

      !-------------------------------------------------------------------
      ! upwind transport
      !-------------------------------------------------------------------

      do n = 1, narrays

         worka(:) = c0
         workb(:) = c0

         do is = 1, ncells
            if (skipE(is)) cycle
            isE = nbr(is,nbr_E)
            worka(is) =     &
               upwind(phi(is,n),phi(isE,n),uee(is),HTE(is),dt)
         enddo

         do is = 1, ncells
            if (skipN(is)) cycle
            isN = nbr(is,nbr_N)
            workb(is) =     &
               upwind(phi(is,n),phi(isN,n),vnn(is),HTN(is),dt)
         enddo

         do is = 1, ncells
            if (skipme(is)) cycle
            isW = nbr(is,nbr_W)
            isS = nbr(is,nbr_S)
            phi(is,n) = phi(is,n) - ( worka(is)-worka(isW)    &
                                    + workb(is)-workb(isS) )  &
                                    / tarea(is)
         enddo

      enddo                     ! narrays

      end subroutine upwind_field

!=======================================================================
!
! Define upwind function
!

      real(kind=dbl_kind) function upwind(y1,y2,a,h,dt)

      real(kind=dbl_kind), intent(in) :: y1,y2,a,h,dt

      upwind = p5*dt*h*((a+abs(a))*y1+(a-abs(a))*y2)

      end function upwind

!=======================================================================

      end module ice_transport_driver

!=======================================================================
