!=======================================================================
!
! Primary state variables in various configurations
! Note: other state variables are at the end of this...
! The primary state variable names are:
!-------------------------------------------------------------------
! for each category   aggregated over     units
!                       categories
!-------------------------------------------------------------------
! aicen(i,j,n)         aice(i,j)           ---
! vicen(i,j,n)         vice(i,j)           m
! vsnon(i,j,n)         vsno(i,j)           m
! trcrn(i,j,it,n)      trcr(i,j,it)
!
! Area is dimensionless because aice is the fractional area
! (normalized so that the sum over all categories, including open
! water, is 1.0).  That is why vice/vsno have units of m instead of m^3.
!
! Variable names follow these rules:
!
! (1) For 3D variables (indices i,j,n), write 'ice' or 'sno' or
!     'sfc' and put an 'n' at the end.
! (2) For 2D variables (indices i,j) aggregated over all categories,
!     write 'ice' or 'sno' or 'sfc' without the 'n'.
! (3) For 2D variables (indices i,j) associated with an individual
!     category, write 'i' or 's' instead of 'ice' or 'sno' and put an 'n'
!     at the end: e.g. hin, hsn.  These are not declared here
!     but in individual modules (e.g., icepack_therm_vertical).
!
! authors C. M. Bitz, UW
!         Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free form source (F90) by Elizabeth Hunke

      module ice_state

      use ice_kinds_mod
      use ice_constants, only: field_loc_center, field_type_scalar, c0
      use ice_domain_size, only: max_blocks, ncat
      use ice_blocks, only: nx_block, ny_block
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_diag
      use icepack_intfc, only: icepack_query_tracer_sizes
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

      implicit none
      private
      public :: alloc_state, alloc_state_slots

      !-----------------------------------------------------------------
      ! state of the ice aggregated over all categories
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(:,:,:), allocatable, public :: &
         aice  , & ! concentration of ice on T grid
         vice  , & ! volume per unit area of ice          (m)
         vsno      ! volume per unit area of snow         (m)

      ! the i and j dimensions are flattened: element (i,j,iblk) of the block
      ! layout lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), dimension(:,:), allocatable, target, public :: &
         aiU       ! concentration of ice on U grid

      real (kind=dbl_kind), dimension(:,:,:,:), allocatable, public :: &
         trcr      ! ice tracers
                   ! 1: surface temperature of ice/snow (C)

      !-----------------------------------------------------------------
      ! state of the ice for each category
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (:,:,:), allocatable, public:: &
         aice0     ! concentration of open water

      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable, public :: &
         trcrn     ! tracers
                   ! 1: surface temperature of ice/snow (C)

      !-----------------------------------------------------------------
      ! The same state in the compressed layout.
      !
      ! ⚠ SLOT FIRST, then category and tracer.  Three separate pieces of
      ! infrastructure require it, which is a stronger reason than the one
      ! this comment used to give:
      !
      !   - ice_HaloUpdate3DR8_compressed takes array(:,:) with the SLOT in
      !     dimension 1 and a BUNDLE index in dimension 2.  bound_state
      !     halo-updates aicen/vicen/vsnon/trcrn every timestep, so this is
      !     a per-step requirement, not a corner case -- and slot-first lets
      !     all ncat levels go in ONE bundled call.
      !   - aim(dynTot,0:ncat) and trm(dynTot,ntrace,ncat), the transport
      !     driver's work arrays, are slot-first, so state_to_tracers copies
      !     run contiguous on both sides.
      !   - the aggregate slot arrays above are already slot-indexed.
      !
      ! ⚠ The cost this pays is icepack's, and it is cache locality rather
      ! than copying.  icepack takes a whole cell's categories as one
      ! assumed-shape dummy, so with slot first those ncat (or ntrcr*ncat)
      ! values lie dynTot apart -- a separate cache line each.  It is NOT a
      ! copy: icepack's dummies are dimension(:) with no CONTIGUOUS
      ! attribute anywhere in it (282 assumed-shape against 25
      ! explicit-shape, and the one explicit dimension(0:ncat) is hin_max,
      ! not state), so a strided actual is passed by descriptor.  An earlier
      ! version of this comment claimed a copy-in/copy-out per call; that
      ! was wrong, and -check noarg_temp_created in the ifx flags is
      ! unrelated to it.
      !
      ! Category-first was tried on paper and rejected: it would make the
      ! halo need a transpose into scratch every step, so it trades one
      ! per-cell locality cost for two per-step whole-state costs.
      !
      ! ⚠ These are NOT yet the live copy.  They are the storage the
      ! loop-by-loop conversion needs before it can start: a loop that calls
      ! icepack per cell passes eight or nine of these at once, so the unit
      ! of conversion is a LOOP, not a variable, and no single variable can
      ! move on its own.  Until a given loop is converted, the block form
      ! above stays authoritative for it.
      !
      ! ⚠ Allocated by alloc_state_slots, NOT alloc_state: dynTot does not
      ! exist yet when alloc_state runs.  init_grid2 builds the ragged layout
      ! as its last act and calls alloc_state_slots there.
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(:), allocatable, public :: &
         aice_slot  , & ! concentration of ice on T grid, compressed
         vice_slot  , & ! volume per unit area of ice  (m), compressed
         vsno_slot  , & ! volume per unit area of snow (m), compressed
         aice0_slot     ! concentration of open water, compressed

      real (kind=dbl_kind), dimension(:,:), allocatable, public :: &
         trcr_slot  , & ! ice tracers              (slot,ntrcr)
         aicen_slot , & ! concentration of ice     (slot,ncat)
         vicen_slot , & ! volume per unit area of ice  (m) (slot,ncat)
         vsnon_slot     ! volume per unit area of snow (m) (slot,ncat)

      ! ⚠ target: ice_HaloUpdate3DR8_compressed takes rank-2 (slot,bundle),
      ! so trcrn_slot is halo-updated through trcrn_flat, a contiguous
      ! rank-remapped view of it -- one bundled call instead of ncat.
      real (kind=dbl_kind), dimension(:,:,:), allocatable, target, public :: &
         trcrn_slot     ! tracers                  (slot,ntrcr,ncat)

      real (kind=dbl_kind), dimension(:,:), pointer, public :: &
         trcrn_flat => null()   ! trcrn_slot viewed as (slot, ntrcr*ncat)


      !-----------------------------------------------------------------
      ! tracers infrastructure arrays
      !-----------------------------------------------------------------

      integer (kind=int_kind), dimension (:), allocatable, public :: &
         trcr_depend   ! = 0 for ice area tracers
                       ! = 1 for ice volume tracers
                       ! = 2 for snow volume tracers

      integer (kind=int_kind), dimension (:), allocatable, public :: &
         n_trcr_strata ! number of underlying tracer layers

      integer (kind=int_kind), dimension (:,:), allocatable, public :: &
         nt_strata     ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), allocatable, public :: &
         trcr_base     ! = 0 or 1 depending on tracer dependency
                       ! argument 2:  (1) aice, (2) vice, (3) vsno

      !-----------------------------------------------------------------
      ! dynamic variables closely related to the state of the ice
      !-----------------------------------------------------------------

      ! These have their i and j dimensions flattened: element (i,j,iblk)
      ! lives at ((j-1)*nx_block + i, iblk), as aiU and aice_init below
      ! already do.  divu, shear and vort joined them so the deformation
      ! kernels can take a rank-1 (:,iblk) slice instead of relying on
      ! sequence association from a rank-2 one.
      real (kind=dbl_kind), dimension(:,:), allocatable, public :: &
         uvel     , & ! x-component of velocity on U grid (m/s)
         vvel     , & ! y-component of velocity on U grid (m/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         shear    , & ! strain rate II component (1/s)
         vort         ! vorticity (1/s)

      real (kind=dbl_kind), dimension(:,:,:), allocatable, public :: &
         uvelE    , & ! x-component of velocity on E grid (m/s)
         vvelE    , & ! y-component of velocity on E grid (m/s)
         uvelN    , & ! x-component of velocity on N grid (m/s)
         vvelN    , & ! y-component of velocity on N grid (m/s)
         strength     ! ice strength (N/m)

      !-----------------------------------------------------------------
      ! ice state at start of time step, saved for later in the step
      !-----------------------------------------------------------------

      ! the i and j dimensions are flattened: element (i,j,iblk) of the block
      ! layout lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), dimension(:,:), allocatable, target, public :: &
         aice_init       ! initial concentration of ice, for diagnostics

      ! rank-3 alias of aice_init, for routines that take (i,j,iblk) arrays
      real (kind=dbl_kind), dimension(:,:,:), pointer, public :: &
         aice_init3d => null()   ! initial concentration of ice, for diagnostics

      real (kind=dbl_kind), dimension(:,:,:,:), allocatable, public :: &
         aicen_init  , & ! initial ice concentration, for linear ITD
         vicen_init  , & ! initial ice volume (m), for linear ITD
         vsnon_init  , & ! initial snow volume (m), for aerosol
         Tsfcn_init      ! initial ice surface temperature (degC)

!=======================================================================

      contains

!=======================================================================
!
! Allocate the compressed copies of the state.
!
! Separate from alloc_state because of ordering: CICE_InitMod calls
! alloc_state before init_grid2, and dynTot is not known until
! init_ragged_layout runs at the end of init_grid2.  Rather than add a call
! to all fifteen drivers' CICE_InitMod, init_ragged_layout calls this itself
! once the layout exists.
!
      subroutine alloc_state_slots(nslot)

      integer (int_kind), intent(in) :: nslot   ! dynTot

      integer (int_kind) :: ierr, ntrcr
      character(len=*),parameter :: subname='(alloc_state_slots)'

      ! ⚠ Only the aggregates, deliberately.  Holding both layouts is a
      ! transitional cost and it is not small: the full set is
      ! 4 + ntrcr + 3*ncat + ntrcr*ncat words per slot, which at ntrcr=19,
      ! ncat=5 is 139 -- about 1.2 GB of duplicate state on a single-task
      ! quarter-degree run, and trcrn_slot alone is 95 of those 139.
      !
      ! So each array is allocated by the commit that converts the loop
      ! which first reads it, not up front.  The declarations above are
      ! complete; this list grows as the conversion proceeds and the block
      ! forms it replaces are retired.
      ! ntrcr is a runtime size, so it is queried here exactly as
      ! alloc_state does rather than threaded through the caller.
      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      allocate ( &
         aice_slot (nslot) , & ! concentration of ice on T grid
         vice_slot (nslot) , & ! volume per unit area of ice  (m)
         vsno_slot (nslot) , & ! volume per unit area of snow (m)
         aice0_slot(nslot) , & ! concentration of open water
         aicen_slot(nslot,ncat)       , & ! concentration of ice, by category
         vicen_slot(nslot,ncat)       , & ! ice volume  (m), by category
         vsnon_slot(nslot,ncat)       , & ! snow volume (m), by category
         trcrn_slot(nslot,ntrcr,ncat) , & ! tracers
         trcr_slot (nslot,ntrcr)      , & ! tracers, aggregated
         stat=ierr)
      if (ierr/=0) call abort_ice(error_message=trim(subname)//' ERROR: Out of memory', &
         file=__FILE__, line=__LINE__)

      aice_slot  = c0
      vice_slot  = c0
      vsno_slot  = c0
      aice0_slot = c0
      aicen_slot = c0
      vicen_slot = c0
      vsnon_slot = c0
      trcrn_slot = c0
      trcr_slot  = c0

      ! A contiguous rank-2 view for the bundled compressed halo.  Rank
      ! remapping needs the target contiguous, which an allocatable is.
      trcrn_flat(1:nslot,1:ntrcr*ncat) => trcrn_slot

      end subroutine alloc_state_slots


!=======================================================================
!
! Allocate space for all state variables
!
      subroutine alloc_state
      integer (int_kind) :: ntrcr, ierr
      character(len=*),parameter :: subname='(alloc_state)'

      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
          file=__FILE__, line=__LINE__)

      allocate ( &
         aice      (nx_block,ny_block,max_blocks) , & ! concentration of ice T grid
         aiU       (nx_block*ny_block,max_blocks) , & ! concentration of ice U grid
         vice      (nx_block,ny_block,max_blocks) , & ! volume per unit area of ice (m)
         vsno      (nx_block,ny_block,max_blocks) , & ! volume per unit area of snow (m)
         aice0     (nx_block,ny_block,max_blocks) , & ! concentration of open water
         uvel      (nx_block*ny_block,max_blocks) , & ! x-component of velocity on U grid (m/s)
         vvel      (nx_block*ny_block,max_blocks) , & ! y-component of velocity on U grid (m/s)
         uvelE     (nx_block,ny_block,max_blocks) , & ! x-component of velocity on E grid (m/s)
         vvelE     (nx_block,ny_block,max_blocks) , & ! y-component of velocity on E grid (m/s)
         uvelN     (nx_block,ny_block,max_blocks) , & ! x-component of velocity on N grid (m/s)
         vvelN     (nx_block,ny_block,max_blocks) , & ! y-component of velocity on N grid (m/s)
         divu      (nx_block*ny_block,max_blocks) , & ! strain rate I component, velocity divergence (1/s)
         shear     (nx_block*ny_block,max_blocks) , & ! strain rate II component (1/s)
         vort      (nx_block*ny_block,max_blocks) , & ! vorticity (1/s)
         strength  (nx_block,ny_block,max_blocks) , & ! ice strength (N/m)
         aice_init (nx_block*ny_block,max_blocks) , & ! initial concentration of ice, for diagnostics
         aicen     (nx_block,ny_block,ncat,max_blocks) , & ! concentration of ice
         vicen     (nx_block,ny_block,ncat,max_blocks) , & ! volume per unit area of ice (m)
         vsnon     (nx_block,ny_block,ncat,max_blocks) , & ! volume per unit area of snow (m)
         aicen_init(nx_block,ny_block,ncat,max_blocks) , & ! initial ice concentration, for linear ITD
         vicen_init(nx_block,ny_block,ncat,max_blocks) , & ! initial ice volume (m), for linear ITD
         vsnon_init(nx_block,ny_block,ncat,max_blocks) , & ! initial snow volume (m), for aerosol
         Tsfcn_init(nx_block,ny_block,ncat,max_blocks) , & ! initial snow/ice surface temperature(degC)
         trcr      (nx_block,ny_block,ntrcr,max_blocks) , & ! ice tracers: 1: surface temperature of ice/snow (C)
         trcrn     (nx_block,ny_block,ntrcr,ncat,max_blocks) , & ! tracers: 1: surface temperature of ice/snow (C)
         stat=ierr)
      if (ierr/=0) call abort_ice(error_message=trim(subname)//' ERROR: Out of memory1', &
         file=__FILE__, line=__LINE__)

      ! rank-3 view of the flattened array; no copy is made
      aice_init3d(1:nx_block,1:ny_block,1:max_blocks) => aice_init

      aice  = c0
      aiU   = c0
      vice  = c0
      vsno  = c0
      aice0 = c0
      uvel  = c0
      vvel  = c0
      uvelE = c0
      vvelE = c0
      uvelN = c0
      vvelN = c0
      divu  = c0
      shear = c0
      vort  = c0
      strength   = c0
      aice_init  = c0
      aicen = c0
      vicen = c0
      vsnon = c0
      aicen_init = c0
      vicen_init = c0
      vsnon_init = c0
      Tsfcn_init = c0
      trcr  = c0
      trcrn = c0

      allocate ( &
         trcr_depend(ntrcr)   , & !
         n_trcr_strata(ntrcr) , & ! number of underlying tracer layers
         nt_strata(ntrcr,2)   , & ! indices of underlying tracer layers
         trcr_base(ntrcr,3)   , & ! = 0 or 1 depending on tracer dependency, (1) aice, (2) vice, (3) vsno
         stat=ierr)
      if (ierr/=0) call abort_ice(error_message=trim(subname)//' ERROR: Out of memory2', &
         file=__FILE__, line=__LINE__)

      trcr_depend = 0
      n_trcr_strata = 0
      nt_strata = 0
      trcr_base = c0

      end subroutine alloc_state

!=======================================================================

      end module ice_state

!=======================================================================
