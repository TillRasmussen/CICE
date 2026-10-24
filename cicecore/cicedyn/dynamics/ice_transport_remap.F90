!=======================================================================
! Transports quantities using the second-order conservative remapping
! scheme developed by John Dukowicz and John Baumgardner (DB) and modified
! for sea ice by William Lipscomb and Elizabeth Hunke.
!
! References:
!
! Dukowicz, J. K., and J. R. Baumgardner, 2000: Incremental
!  remapping as a transport/advection algorithm, J. Comput. Phys.,
!  160, 318-335.
!
! Lipscomb, W. H., and E. C. Hunke, 2004: Modeling sea ice
!  transport using incremental remapping, Mon. Wea. Rev., 132,
!  1341-1354.
!
! Lemieux, J.-F., Lipscomb, W. H., Craig, A., Bailey, D. A.,
! Hunke, E. C., Blain, P., Rasmussen, T. A. S., Bentsen, M.,
! Dupont, F., Hebert, D., and Allard, R., 2024: CICE on a
! C-grid: new momentum, stress, and transport schemes for CICEv6.5,
! Geosci. Model Dev., 17, 6703-6724.
!
! authors William H. Lipscomb, LANL
!         John Baumgardner, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004-05: Block structure added (WHL)
! 2006: Moved remap driver to ice_transport_driver
!       Geometry changes:
!       (1) Reconstruct fields in stretched logically rectangular coordinates
!       (2) Modify geometry so that the area flux across each edge
!           can be specified (following an idea of Mats Bentsen)
! 2010: ECH removed unnecessary grid arrays and optional arguments from
!       horizontal_remap
! 2023: TAR, DMI Remove commented code and unnecessary arrays
! 2024: JFL and WHL improved robustness of Bentsen's approach for area flux.

      module ice_transport_remap

      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_calendar, only: istep1
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c0, c1, c2, c12, p333, p4, p5, p6, &
          eps13, eps16, &
          field_loc_center, field_type_scalar, &
          field_loc_NEcorner, field_type_vector
      use ice_diagnostics, only: diagnostic_abort
      use ice_domain_size, only: max_blocks, ncat
      use ice_fileunits, only: nu_diag
      use ice_exit, only: abort_ice
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters
      use ice_grid, only : grid_ice

      implicit none
      private
      public :: init_remap, horizontal_remap, make_masks

      integer (kind=int_kind), parameter ::     &
         ngroups  = 6      ,&! number of groups of triangles that
                             ! contribute transports across each edge
         nvert = 3           ! number of vertices in a triangle

      ! for triangle integral formulas
      real (kind=dbl_kind), parameter ::  &
         p5625m = -9._dbl_kind/16._dbl_kind    ,&
         p52083 = 25._dbl_kind/48._dbl_kind

      logical :: &
         l_edge_flux_adj ! if true, prescribe area flux across each edge
                         ! if false, area flux is determined internally
                         ! and is passed out

      logical (kind=log_kind), parameter :: bugcheck = .false.

!=======================================================================
! Here is some information about how the incremental remapping scheme
! works in CICE and how it can be adapted for use in other models.
!
! The remapping routine is designed to transport a generic mass-like
! field (in CICE, the ice fractional area) along with an arbitrary number
! of tracers in two dimensions.  The velocity components are assumed
! to lie at grid cell corners and the transported scalars at cell centers.
! Incremental remapping has the following desirable properties:
!
! (1) Tracer monotonicity is preserved.  That is, no new local
!     extrema are produced in fields like ice thickness or internal
!     energy.
! (2) The reconstucted mass and tracer fields vary linearly in x and y.
!     This means that remapping is 2nd-order accurate in space,
!     except where horizontal gradients are limited to preserve
!     monotonicity.
! (3) There are economies of scale.  Transporting a single field
!     is rather expensive, but additional fields have a relatively
!     low marginal cost.
!
! The following generic conservation equations may be solved:
!
!            dm/dt = del*(u*m)             (0)
!       d(m*T1)/dt = del*(u*m*T1)          (1)
!    d(m*T1*T2)/dt = del*(u*m*T1*T2)       (2)
! d(m*T1*T2*T3)/dt = del*(u*m*T1*T2*T3)    (3)
!
! where d is a partial derivative, del is the 2D divergence operator,
! u is the horizontal velocity, m is the mass density field, and
! T1, T2, and T3 are tracers.
!
! In CICE, these equations have the form
!
!               da/dt = del*(u*a)          (4)
! dv/dt =   d(a*h)/dt = del*(u*a*h)        (5)
! de/dt = d(a*h*q)/dt = del*(u*a*h*q)      (6)
!            d(aT)/dt = del*(u*a*t)        (7)
!
! where a = fractional ice area, v = ice/snow volume, h = v/a = thickness,
! e = ice/snow internal energy (J/m^2), q = e/v = internal energy per
! unit volume (J/m^3), and T is a tracer.  These equations express
! conservation of ice area, volume, internal energy, and area-weighted
! tracer, respectively.
!
! (Note: In CICE, a, v and e are prognostic quantities from which
!  h and q are diagnosed.  The remapping routine works with tracers,
!  which means that h and q must be derived from a, v, and e before
!  calling the remapping routine.)
!
! Earlier versions of CICE assumed fixed ice and snow density.
! Beginning with CICE 4.0, the ice and snow density can be variable.
! In this case, equations (5) and (6) are replaced by
!
! dv/dt =        d(a*h)/dt = del*(u*a*h)          (8)
! dm/dt =    d(a*h*rho)/dt = del*(u*a*h*rho)      (9)
! de/dt = d(a*h*rho*qm)/dt = del*(u*a*h*rho*qm)   (10)
!
! where rho = density and qm = internal energy per unit mass (J/kg).
! Eq. (9) expresses mass conservation, which in the variable-density
! case is no longer equivalent to volume conservation (8).
!
! Tracers satisfying equations of the form (1) are called "type 1."
! In CICE the paradigmatic type 1 tracers are hi and hs.
!
! Tracers satisfying equations of the form (2) are called "type 2".
! The paradigmatic type 2 tracers are qi and qs (or rhoi and rhos
!  in the variable-density case).
!
! Tracers satisfying equations of the form (3) are called "type 3."
! The paradigmatic type 3 tracers are qmi and qms in the variable-density
! case.  There are no such tracers in the constant-density case.
!
! The fields a, T1, and T2 are reconstructed in each grid cell with
! 2nd-order accuracy.  T3 is reconstructed with 1st-order accuracy
! (i.e., it is transported in upwind fashion) in order to avoid
! additional mathematical complexity.
!
! The mass-like field lives in the array "mm" (shorthand for mean
! mass) and the tracers fields in the array "tm" (mean tracers).
! In order to transport tracers correctly, the remapping routine
! needs to know the tracers types and relationships.  This is done
! as follows:
!
! Each field in the "tm" array is assigned an index, 1:ntrace.
! (Note: ntrace is not the same as ntrcr, the number of tracers
! in the trcrn state variable array.  For remapping purposes we
! have additional tracers hi and hs.)
!
! The tracer types (1,2,3) are contained in the "tracer_type" array.
! For standard CICE:
!
!     tracer_type = (1 1 1 2 2 2 2 2)
!
! Type 2 and type 3 tracers are said to depend on type 1 tracers.
! For instance, qi depends on hi, which is to say that
! there is a conservation equation of the form (2) or (6).
! Thus we define a "depend" array.  For standard CICE:
!
!          depend = (0 0 0 1 1 1 1 2)
!
! which implies that elements 1-3 (hi, hs, Ts) are type 1,
! elements 4-7 (qi) depend on element 1 (hi), and element 8 (qs)
! depends on element 2 (hs).
!
! We also define a logical array "has_dependents".  In standard CICE:
!
!  has_dependents = (T T F F F F F F),
!
! which means that only elements 1 and 2 (hi and hs) have dependent
! tracers.
!
! For the variable-density case, things are a bit more complicated.
! Suppose we have 4 variable-density ice layers and one variable-
! density snow layer.  Then the indexing is as follows:
! 1    = hi
! 2    = hs
! 3    = Ts
! 4-7  = rhoi
! 8    = rhos
! 9-12 = qmi
! 13   = qms
!
! The key arrays are:
!
!    tracer_type = (1 1 1 2 2 2 2 2 3 3 3 3 3)
!
!         depend = (0 0 0 1 1 1 1 2 4 5 6 7 8)
!
! has_dependents = (T T F T T T T T F F F F F)
!
! which imply that hi and hs are type 1 with dependents rhoi and rhos,
! while rhoi and rhos are type 2 with dependents qmi and qms.
!
! Tracers added to the ntrcr array are handled automatically
! by the remapping with little extra coding.  It is necessary
! only to provide the correct type and dependency information.
!
! When using this routine in other models, most of the tracer dependency
! apparatus may be irrelevant.  In a layered ocean model, for example,
! the transported fields are the layer thickness h (the mass density
! field) and two or more tracers (T, S, and various trace species).
! Suppose there are just two tracers, T and S.  Then the tracer arrays
! have the values:
!
!    tracer_type = (1 1)
!         depend = (0 0)
! has_dependents = (F F)
!
! which is to say that all tracer transport equations are of the form (1).
!
! The tracer dependency arrays are optional input arguments for the
! main remapping subroutine.  If these arrays are not passed in, they
! take on the default values tracer_type(:) = 1, depend(:) = 0, and
! has_dependents(:) = F, which are appropriate for most purposes.
!
! Another optional argument is integral_order.  If integral_order = 1,
! then the triangle integrals are exact for linear functions of x and y.
! If integral_order = 2, these integrals are exact for both linear and
! quadratic functions.  If integral_order = 3, integrals are exact for
! cubic functions as well.  If all tracers are of type 1, then the
! integrals of mass*tracer are quadratic, and integral_order = 2 is
! sufficient.  In CICE, where there are type 2 tracers, we integrate
! functions of the form mass*tracer1*tracer2.  Thus integral_order = 3
! is required for exactness, though integral_order = 2 may be good enough
! in practice.
!
! Finally, a few words about the edgearea fields:
!
! In earlier versions of this scheme, the divergence of the velocity
! field implied by the remapping was, in general, different from the
! value of del*u computed in the dynamics.  For energetic consistency
! (in CICE as well as in layered ocean models such as HYPOP),
! these two values should agree.  This can be ensured by using the
! edge flux adjustment (EFA) method by setting l_edge_flux_adj = T.
! The EFA method specifies the area transported across each grid cell
! edge in the arrays edgearea_e and edgearea_n.  The departure regions
! are then tweaked, following an idea by Mats Bentsen, such that they
! have the desired area.  If l_edge_flux_adj = F, these regions are
! not tweaked, and the edgearea arrays are output variables.
!
!=======================================================================

      contains

!=======================================================================
!
! Grid quantities used by the remapping transport scheme
!
! Note: Arrays needed for nonuniform grids has been deleted.
!       They can be found in version 6.5 and earlier
!
! author William H. Lipscomb, LANL

      subroutine init_remap

      use ice_blocks, only: nghost, ew_boundary_type, ns_boundary_type

      character(len=*), parameter :: subname = '(init_remap)'

      ! check that if open boundaries are used, nghost > 1
      ! this was moved from initialization to support unit testing
      ! with open bcs and nghost=1 and support for upwind advection and nghost=1
      if ((ew_boundary_type == 'zero_gradient' .or. &
           ew_boundary_type == 'linear_extrap' .or. &
           ns_boundary_type == 'zero_gradient' .or. &
           ns_boundary_type == 'linear_extrap') .and. &
          (nghost == 1)) then
         if (my_task == master_task) then
            write(nu_diag,*) subname//' ERROR: regional bcs require nghost>=2 in remap advection'
            write(nu_diag,*) subname//' boundary types = '//trim(ew_boundary_type)//' '//trim(ns_boundary_type)
            write(nu_diag,*) subname//' with nghost = ', nghost
         endif
         call abort_ice(subname//' ERROR: regional bcs require nghost >= 2 in remap advection', file=__FILE__, line=__LINE__)
      endif

      !-------------------------------------------------------------------
      ! Set logical l_edge_flux_adj depending of the grid type.
      !
      ! If l_edge_flux_adj is true, the area of each departure region is
      ! computed in advance (e.g., by taking the divergence of the
      ! velocity field and passed to locate_triangles.  The departure
      ! regions are adjusted to obtain the desired area.
      ! If false, edgearea is computed in locate_triangles and passed out.
      !
      ! l_edge_flux_adj = .false. has been the default approach in CICE.
      ! It is used like this for the B-grid. However, idealized tests with
      ! the C-grid have shown that l_edge_flux_adj = .false. leads to a
      ! checkerboard pattern in prognostic fields such as aice Using
      ! l_edge_flux_adj = .true. eliminates the checkerboard pattern in
      ! C-grid simulations (Lemieux et al. 2024).
      !-------------------------------------------------------------------

      if (grid_ice == 'CD' .or. grid_ice == 'C') then
         l_edge_flux_adj = .true.
      else
         l_edge_flux_adj = .false.
      endif

      end subroutine init_remap

!=======================================================================
!
! Solve the transport equations for one timestep using the incremental
! remapping scheme developed by John Dukowicz and John Baumgardner (DB)
! and modified for sea ice by William Lipscomb and Elizabeth Hunke.
!
! This scheme preserves monotonicity of ice area and tracers.  That is,
! it does not produce new extrema.  It is second-order accurate in space,
! except where gradients are limited to preserve monotonicity.
!
! This version of the remapping allows the user to specify the areal
! flux across each edge, based on an idea developed by Mats Bentsen.
!
! author William H. Lipscomb, LANL
! 2006: Moved driver (subroutine transport_remap) into separate module.
!       Geometry changes (logically rectangular coordinates, fixed
!        area fluxes)

      subroutine horizontal_remap (dt,             ntrace,   &
                                   uvel,           vvel,     &
                                   mm,             tm,       &
                                   tracer_type,    depend,   &
                                   has_dependents,           &
                                   integral_order,           &
                                   l_dp_midpt,               &
                                   uvelE,          vvelN)

      use ice_boundary, only: ice_halo, ice_HaloMask, ice_HaloUpdate, &
          ice_HaloDestroy
      use ice_domain, only: nblocks, blocks_ice, halo_info, maskhalo_remap
      use ice_blocks, only: block, get_block, nghost
      use ice_boundary, only: ice_HaloUpdate2DR8_compressed, &
                              ice_HaloUpdate3DR8_compressed
      use ice_grid, only: HTE, HTN, dxu, dyu,       &
                          earea, narea, hm,   &
                          dynMap, dynOff, dynCnt, dynNbr, dynIJ, dynTot, skipRmp, &
                          skipPhys, tarear_slot,                          &
                          dxU_slot, dyU_slot, HTN_slot, HTE_slot,         &
                          HTNe_slot, HTEn_slot, pack_to_slots
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound

      real (kind=dbl_kind), intent(in) ::     &
         dt       ! time step

      integer (kind=int_kind), intent(in) :: &
         ntrace   ! number of tracers in use

      ! uvel/vvel have their i and j dimensions flattened: element (i,j,iblk)
      ! lives at ((j-1)*nx_block + i, iblk)
      real (kind=dbl_kind), intent(in), dimension(nx_block*ny_block,max_blocks) :: &
         uvel,  & ! x-component of velocity (m/s) ugrid
         vvel     ! y-component of velocity (m/s) ugrid

      real (kind=dbl_kind), intent(in), optional, dimension(nx_block,ny_block,max_blocks) :: &
         uvelE, & ! x-component of velocity (m/s) egrid
         vvelN    ! y-component of velocity (m/s) ngrid

      ! Compressed: (slot, category) and (slot, tracer, category), over all
      ! local blocks.  These are the driver's aim/trm.
      real (kind=dbl_kind), intent(inout), dimension (dynTot,0:ncat) :: &
         mm       ! mean mass values in each grid cell

      real (kind=dbl_kind), intent(inout), dimension (dynTot,ntrace,ncat) :: &
         tm       ! mean tracer values in each grid cell

      integer (kind=int_kind), dimension (ntrace), intent(in) :: &
         tracer_type    , & ! = 1, 2, or 3 (see comments above)
         depend             ! tracer dependencies (see above)

      logical (kind=log_kind), dimension (ntrace), intent(in) :: &
         has_dependents     ! true if a tracer has dependent tracers

      integer (kind=int_kind), intent(in) :: &
         integral_order     ! polynomial order for triangle integrals

      logical (kind=log_kind), intent(in) :: &
         l_dp_midpt         ! if true, find departure points using
                            ! corrected midpoint velocity
      ! local variables

      integer (kind=int_kind) ::     &
         i, j           , & ! horizontal indices
         iblk           , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         nc,              & ! slots held by the current block
         ms,              & ! slot of a cell, for the edgearea fill
         is, ijb,         & ! slot index and the cell it owns
         n, m               ! ice category, tracer indices

      ! ⚠ Rank-2 alias of tc/tx/ty for the compressed halo, which takes
      ! (slot, bundle).  Pointer remapping rather than a copy, so the three
      ! tracer halos stay one collective each.
      real (kind=dbl_kind), pointer :: p2(:,:)

      ! ⚠ These are PRIVATE to the block loop below, so each thread gets its
      ! own copy -- and on the STACK they are what made that loop fail.  At
      ! a production block size with ntrace ~ 30 and ncat 5, mtflxe and
      ! mtflxn alone are tens of MB per thread; xp and yp add nx*ny*4*6.
      ! Allocatable means each thread's copy comes from the heap instead,
      ! the same remedy 8bb517d applied to VP's Krylov bases.
      real (kind=dbl_kind), dimension(:), allocatable ::  &
         edgearea_e     , & ! area of departure regions for east edges
         edgearea_n         ! area of departure regions for north edges, slots

      ! The velocities, packed per block for departure_points.
      !
      ! ⚠ nx_block*ny_block + 1, not nx_block*ny_block.  A block whose every
      ! cell is mapped has dynCnt = nx_block*ny_block -- boxislands with the
      ! full ghost rim does exactly that -- and the slice passed is
      ! 1:dynCnt+1, because pack_to_slots zeroes the LAST element it is given
      ! as the dead slot.  Sized to the cell count this overran by one.
      real (kind=dbl_kind), dimension(nx_block*ny_block+1) ::  &
         uvel_s, vvel_s

      ! Compressed, and they stay that way across the halo exchange: the
      ! whole path from departure_points to locate_triangles is slots now.
      ! (dynTot) rather than (nx_block,ny_block,max_blocks).
      real (kind=dbl_kind), dimension (dynTot) ::     &
         dpx            , & ! x coordinates of departure points at cell corners
         dpy                ! y coordinates of departure points at cell corners

      ! Compressed.  construct_fields already computed these as slots and
      ! then unpacked them; transport_integrals reads them at the DONOR
      ! cell, which is now a slot, so the block form has no reader left.
      real (kind=dbl_kind), dimension(dynTot,0:ncat) :: &
         mc             , & ! mass at geometric center of cell
         mx, my             ! limited derivative of mass wrt x and y

      ! Compressed, one block at a time, so allocated inside the loop like
      ! the rest of the per-block working set -- at dynTot they would be a
      ! per-thread stack array.
      real (kind=dbl_kind), dimension(:,:), allocatable :: &
         mmask              ! = 1. if mass is present, = 0. otherwise

      ! ⚠ TARGET so the halo call below can alias (slot,ntrace,ncat) to the
      ! rank 2 (slot, bundle) the compressed halo takes -- the same pointer
      ! remapping ice_HaloUpdate3DR8_compressed uses internally.  One call
      ! per field, not one per category, so the collective count is
      ! unchanged.
      real (kind=dbl_kind), dimension (dynTot,ntrace,ncat), target :: &
         tc             , & ! tracer values at geometric center of cell
         tx, ty             ! limited derivative of tracer wrt x and y

      real (kind=dbl_kind), dimension (:,:,:), allocatable ::     &
         tmask              ! = 1. if tracer is present, = 0. otherwise

      ! Compressed: (slot, category) and (slot, tracer, category).
      real (kind=dbl_kind), dimension (:,:), allocatable ::     &
         mflxe, mflxn       ! mass transports across E and N cell edges

      real (kind=dbl_kind), dimension (:,:,:), allocatable :: &
         mtflxe, mtflxn     ! mass*tracer transports across E and N cell edges

      ! Slot-indexed now: (slot, group) and (slot, vertex, group).
      real (kind=dbl_kind), dimension (:,:), allocatable ::     &
         triarea            ! area of east-edge departure triangle

      real (kind=dbl_kind), dimension (:,:,:), allocatable ::  &
         xp, yp             ! x and y coordinates of special triangle points
                            ! (need 4 points for triangle integrals)
      integer (kind=int_kind), dimension (:,:), allocatable ::     &
         iflux          , & ! i index of cell contributing transport
         jflux          , & ! j index of cell contributing transport
         isflux             ! and that donor cell's slot

      integer (kind=int_kind), dimension(ngroups,max_blocks) ::     &
         icellsng           ! number of cells with ice

      integer (kind=int_kind), dimension(:,:), allocatable ::     &
         indxsng              ! the slot each of those cells owns

      integer (kind=int_kind), dimension(nx_block,ny_block,max_blocks) :: &
         halomask           ! temporary mask for fast halo updates

      logical (kind=log_kind) ::     &
         l_stop             ! if true, abort the model

      integer (kind=int_kind) ::     &
         istop, jstop       ! indices of grid cell where model aborts

      character (len=char_len) ::   &
         edge               ! 'north' or 'east'

      type (ice_halo) :: &
         halo_info_tracer   ! masked halo

      type (block) ::     &
         this_block         ! block information for current block

      character(len=*), parameter :: subname = '(horizontal_remap)'

      !-------------------------------------------------------------------
      ! Remap the ice area and associated tracers.
      ! Remap the open water area (without tracers).
      !-------------------------------------------------------------------

      ! ⚠ nc, uvel_s and vvel_s are per-iteration state and MUST be PRIVATE.
      ! They are declared in this routine, not inside the region, so leaving
      ! them off the clause makes them shared and the packs race -- which is
      ! exactly how 8da46d2 reached boreas.  No compiler warns.
      !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block,n,   &
      !$OMP          mmask,tmask,istop,jstop,l_stop, &
      !$OMP          nc,uvel_s,vvel_s) &
      !$OMP SCHEDULE(runtime)
      do iblk = 1, nblocks

         ! PRIVATE, so each thread allocates its own -- see the second block
         ! loop below for why these are not stack arrays.
         allocate(mmask(dynCnt(iblk)+1,0:ncat),         &
                  tmask(dynCnt(iblk)+1,ntrace,ncat))

         l_stop = .false.
         istop = 0
         jstop = 0

         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         !-------------------------------------------------------------------
         ! Compute masks and count ice cells.
         ! Masks are used to prevent tracer values in cells without ice from
         !  being used to compute tracer gradients.
         !-------------------------------------------------------------------

         call make_masks (dynCnt(iblk),                                      &
                          ntrace,                                           &
                          has_dependents,                                   &
                          mm(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:), &
                          mmask,                                            &
                          tm(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:,:),&
                          tmask)

         !-------------------------------------------------------------------
         ! Construct linear fields, limiting gradients to preserve monotonicity.
         ! Note: Pass in unit arrays instead of true distances HTE, HTN, etc.
         !       The resulting gradients are in scaled coordinates.
         !-------------------------------------------------------------------

         ! open water
         call construct_fields(iblk,              ntrace,            &
                               tracer_type,       depend,            &
                               has_dependents,                       &
                               hm     (:, iblk),                    &
                               mm(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,0),                     &
                               mc(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,0),                     &
                               mx(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,0),                     &
                               my(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,0),                     &
                               mmask(:,0) )

         ! ice categories

         do n = 1, ncat

            call construct_fields(iblk,                ntrace,              &
                                  tracer_type,         depend,              &
                                  has_dependents,                          &
                                  hm       (:, iblk),                      &
                                  mm(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,n),                     &
                                  mc(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,n),                     &
                                  mx(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,n),                     &
                                  my(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,n),                     &
                                  mmask  (:,n),                            &
                                  tm(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:,n),                   &
                                  tc(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:,n),                   &
                                  tx(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:,n),                   &
                                  ty(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:,n),                   &
                                  tmask(:,:,n) )

         enddo                  ! n

         !-------------------------------------------------------------------
         ! Given velocity field at cell corners, compute departure points
         ! of trajectories.
         !-------------------------------------------------------------------

         ! uvel and vvel are packed here rather than taken from the
         ! dynamics' uvel_slot/vvel_slot: reaching across to those assumes a
         ! producer kept them current, which is the assumption that cost
         ! 797772a.  The metrics already have slot forms in ice_grid.
         nc = dynCnt(iblk)
         call pack_to_slots(uvel(:,iblk), dynMap(:,iblk), uvel_s(1:nc+1))
         call pack_to_slots(vvel(:,iblk), dynMap(:,iblk), vvel_s(1:nc+1))

         call departure_points(nc,                                            &
                               skipRmp(dynOff(iblk)+1:dynOff(iblk)+nc),       &
                               dynNbr (dynOff(iblk)+1:dynOff(iblk)+nc,:),     &
                               dynIJ  (dynOff(iblk)+1:dynOff(iblk)+nc),       &
                               nx_block,                                      &
                               dt,                                            &
                               uvel_s(1:nc+1),   vvel_s(1:nc+1),              &
                               dxU_slot(dynOff(iblk)+1:dynOff(iblk)+nc+1),    &
                               dyU_slot(dynOff(iblk)+1:dynOff(iblk)+nc+1),    &
                               HTN_slot(dynOff(iblk)+1:dynOff(iblk)+nc+1),    &
                               HTE_slot(dynOff(iblk)+1:dynOff(iblk)+nc+1),    &
                               HTNe_slot(dynOff(iblk)+1:dynOff(iblk)+nc+1),   &
                               HTEn_slot(dynOff(iblk)+1:dynOff(iblk)+nc+1),   &
                               dpx(dynOff(iblk)+1:dynOff(iblk)+nc+1),         &
                               dpy(dynOff(iblk)+1:dynOff(iblk)+nc+1),         &
                               l_dp_midpt,       l_stop,                      &
                               istop,            jstop)

         if (l_stop) then
            call diagnostic_abort(istop,jstop,iblk,'bad departure points')
         endif

         deallocate(mmask, tmask)

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      !-------------------------------------------------------------------
      ! Ghost cell updates
      ! If nghost >= 2, these calls are not needed
      !-------------------------------------------------------------------

      if (nghost==1) then

         call ice_timer_start(timer_bound)

         ! departure points
         call ice_HaloUpdate2DR8_compressed(dpx, halo_info, &
                              field_loc_NEcorner, field_type_vector, dynMap, dynOff)
         call ice_HaloUpdate2DR8_compressed(dpy, halo_info, &
                              field_loc_NEcorner, field_type_vector, dynMap, dynOff)

         ! mass field -- (slot, category) is already the rank 2 the
         ! compressed halo takes, so these need no aliasing
         call ice_HaloUpdate3DR8_compressed(mc, halo_info, &
                              field_loc_center, field_type_scalar, dynMap, dynOff)
         call ice_HaloUpdate3DR8_compressed(mx, halo_info, &
                              field_loc_center, field_type_vector, dynMap, dynOff)
         call ice_HaloUpdate3DR8_compressed(my, halo_info, &
                              field_loc_center, field_type_vector, dynMap, dynOff)

         ! tracer fields
         if (maskhalo_remap) then
            halomask(:,:,:) = 0
            !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,n,m,j,i) SCHEDULE(runtime)
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi

               ! ⚠ The slot loop must carry the block form's bounds: it ran
               ! ilo..ihi, jlo..jhi, and a slot loop would widen that to the
               ! whole block including the ghost rim.  skipRmp is exactly
               ! that range at nghost 1 and wider at nghost 2, so the cell
               ! is tested through dynIJ rather than through the mask.
               do n = 1, ncat
               do m = 1, ntrace
               do is = dynOff(iblk)+1, dynOff(iblk)+dynCnt(iblk)
                  ijb = dynIJ(is)
                  i = mod(ijb-1, nx_block) + 1
                  j = (ijb-1) / nx_block + 1
                  if (i < ilo .or. i > ihi .or. j < jlo .or. j > jhi) cycle
                  if (tc(is,m,n) /= c0) halomask(i,j,iblk) = 1
                  if (tx(is,m,n) /= c0) halomask(i,j,iblk) = 1
                  if (ty(is,m,n) /= c0) halomask(i,j,iblk) = 1
               enddo
               enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
            call ice_HaloUpdate (halomask,         halo_info, &
                                 field_loc_center, field_type_scalar)
            call ice_HaloMask(halo_info_tracer, halo_info, halomask)

            p2(1:dynTot,1:ntrace*ncat) => tc
            call ice_HaloUpdate3DR8_compressed(p2, halo_info_tracer, &
                                 field_loc_center, field_type_scalar, dynMap, dynOff)
            p2(1:dynTot,1:ntrace*ncat) => tx
            call ice_HaloUpdate3DR8_compressed(p2, halo_info_tracer, &
                                 field_loc_center, field_type_vector, dynMap, dynOff)
            p2(1:dynTot,1:ntrace*ncat) => ty
            call ice_HaloUpdate3DR8_compressed(p2, halo_info_tracer, &
                                 field_loc_center, field_type_vector, dynMap, dynOff)
            call ice_HaloDestroy(halo_info_tracer)
         else
            p2(1:dynTot,1:ntrace*ncat) => tc
            call ice_HaloUpdate3DR8_compressed(p2, halo_info, &
                                 field_loc_center, field_type_scalar, dynMap, dynOff)
            p2(1:dynTot,1:ntrace*ncat) => tx
            call ice_HaloUpdate3DR8_compressed(p2, halo_info, &
                                 field_loc_center, field_type_vector, dynMap, dynOff)
            p2(1:dynTot,1:ntrace*ncat) => ty
            call ice_HaloUpdate3DR8_compressed(p2, halo_info, &
                                 field_loc_center, field_type_vector, dynMap, dynOff)
         endif
         call ice_timer_stop(timer_bound)

      endif  ! nghost

      ! tcraig, this OMP loop sometimes fails with cce/14.0.3, compiler bug??
      ! TILL I can trigger the same with ifort (IFORT) 18.0.0 20170811
      !$OMP PARALLEL DO PRIVATE(iblk,i,j,ms,ilo,ihi,jlo,jhi,this_block,n,  &
      !$OMP                     edgearea_e,edgearea_n,edge,iflux,jflux,isflux, &
      !$OMP                     xp,yp,indxsng,mflxe,mflxn, &
      !$OMP                     mtflxe,mtflxn,triarea,istop,jstop,l_stop) &
      !$OMP SCHEDULE(runtime)
      do iblk = 1, nblocks

         ! Each thread allocates its own: these are PRIVATE, and on the
         ! stack they overflowed it.  One allocation per block per timestep
         ! against the triangle geometry costs nothing.
         allocate(edgearea_e(dynCnt(iblk)+1),                     &
                  edgearea_n(dynCnt(iblk)+1),                     &
                  mflxe (dynCnt(iblk)+1,0:ncat),                  &
                  mflxn (dynCnt(iblk)+1,0:ncat),                  &
                  mtflxe(dynCnt(iblk)+1,ntrace,ncat),             &
                  mtflxn(dynCnt(iblk)+1,ntrace,ncat),             &
                  triarea(dynCnt(iblk)+1,ngroups),                &
                  xp(dynCnt(iblk)+1,0:nvert,ngroups),             &
                  yp(dynCnt(iblk)+1,0:nvert,ngroups),             &
                  iflux(dynCnt(iblk)+1,ngroups),                  &
                  jflux(dynCnt(iblk)+1,ngroups),                  &
                  isflux(dynCnt(iblk)+1,ngroups),                 &
                  indxsng(nx_block*ny_block,ngroups))

         l_stop = .false.
         istop = 0
         jstop = 0

         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         !-------------------------------------------------------------------
         ! If l_edge_flux_adj is true, compute edgearea by taking the
         ! divergence of the velocity field. Otherwise, initialize edgearea.
         !-------------------------------------------------------------------

         edgearea_e(:) = c0
         edgearea_n(:) = c0

         if (l_edge_flux_adj) then
            if (grid_ice == 'CD' .or. grid_ice == 'C') then ! velocities are already on the center
               if (.not.present(uvelE).or..not.present(vvelN)) then
                  call abort_ice (subname//'ERROR: uvelE,vvelN required with C|CD and l_edge_flux_adj')
               endif

               do j = jlo, jhi
               do i = ilo-1, ihi
                  ms = dynMap((j-1)*nx_block + i, iblk)
                  if (ms < 1) cycle
                  edgearea_e(ms) = uvelE(i,j,iblk) * HTE((j-1)*nx_block+i,iblk) * dt
               enddo
               enddo

               do j = jlo-1, jhi
               do i = ilo, ihi
                  ms = dynMap((j-1)*nx_block + i, iblk)
                  if (ms < 1) cycle
                  edgearea_n(ms) = vvelN(i,j,iblk)*HTN((j-1)*nx_block+i,iblk) * dt
               enddo
               enddo

            else
               do j = jlo, jhi
               do i = ilo-1, ihi
                  ms = dynMap((j-1)*nx_block + i, iblk)
                  if (ms < 1) cycle
                  edgearea_e(ms) = (uvel((j-1)*nx_block+i,iblk) + uvel((j-1-1)*nx_block+i,iblk)) &
                                        * p5 * HTE((j-1)*nx_block+i,iblk) * dt
               enddo
               enddo

               do j = jlo-1, jhi
               do i = ilo, ihi
                  ms = dynMap((j-1)*nx_block + i, iblk)
                  if (ms < 1) cycle
                  edgearea_n(ms) = (vvel((j-1)*nx_block+i,iblk) + vvel((j-1)*nx_block+i-1,iblk)) &
                                        * p5 * HTN((j-1)*nx_block+i,iblk) * dt
               enddo
               enddo
            endif
         endif

         !-------------------------------------------------------------------
         ! Transports for east cell edges.
         !-------------------------------------------------------------------

         !-------------------------------------------------------------------
         ! Compute areas and vertices of departure triangles.
         !-------------------------------------------------------------------

         edge = 'east'
         call locate_triangles(nx_block,          ny_block,           &
                               ilo, ihi,          jlo, jhi,           &
                               nghost,            edge,               &
                               icellsng(:,iblk),                      &
                               indxsng(:,:),                          &
                               dpx(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),            &
                               dpy(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), &
                               dynIJ(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),  &
                               dxU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               dyU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               earea (:,  iblk),  narea (:,  iblk),   &
                               xp,                yp,                 &
                               iflux,             jflux,              &
                               isflux,                                &
                               triarea,           edgearea_e(:))

         !-------------------------------------------------------------------
         ! Given triangle vertices, compute coordinates of triangle points
         !  needed for transport integrals.
         !-------------------------------------------------------------------

         call triangle_coordinates (integral_order,  icellsng(:,iblk), &
                                    indxsng(:,:),                      &
                                    xp,              yp)

         !-------------------------------------------------------------------
         ! Compute the transport across east cell edges by summing contributions
         ! from each triangle.
         !-------------------------------------------------------------------

         ! open water
         call transport_integrals(                                     &
                                  ntrace,          icellsng (:,iblk),  &
                                  indxsng(:,:),                        &
                                  tracer_type,     depend,             &
                                  integral_order,  triarea,            &
                                  isflux,                                &
                                  xp,              yp,                 &
                                  mc(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,0),          mx(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,0),            &
                                  my(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,0),  mflxe(:,0))

         ! ice categories
         do n = 1, ncat
            call transport_integrals                                     &
                               (ntrace,            icellsng (:,iblk),    &
                                indxsng(:,:),                        &
                                tracer_type,       depend,               &
                                integral_order,    triarea,              &
                                isflux,                                    &
                                xp,                yp,                   &
                                mc(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,n),           mx(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,n),              &
                                my(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,n),  mflxe (:,  n),      &
                                tc(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:,n),         tx(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:,n),            &
                                ty(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:,n),  mtflxe(:,:,n))

         enddo

         !-------------------------------------------------------------------
         ! Repeat for north edges
         !-------------------------------------------------------------------

         edge = 'north'
         call locate_triangles(nx_block,          ny_block,           &
                               ilo, ihi,          jlo, jhi,           &
                               nghost,            edge,               &
                               icellsng(:,iblk),                      &
                               indxsng(:,:),                          &
                               dpx(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),            &
                               dpy(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), dynMap(:,iblk), &
                               dynIJ(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1),  &
                               dxU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               dyU_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               earea (:,  iblk),  narea (:,  iblk),   &
                               xp,                yp,                 &
                               iflux,             jflux,              &
                               isflux,                                &
                               triarea,           edgearea_n(:))

         call triangle_coordinates (integral_order,  icellsng(:,iblk), &
                                    indxsng(:,:),                      &
                                    xp,              yp)

         ! open water
         call transport_integrals(                                    &
                                  ntrace,          icellsng (:,iblk), &
                                  indxsng(:,:),                        &
                                  tracer_type,     depend,            &
                                  integral_order,  triarea,           &
                                  isflux,                               &
                                  xp,              yp,                &
                                  mc(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,0),          mx(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,0),            &
                                  my(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,0),  mflxn(:,0))

         ! ice categories
         do n = 1, ncat
            call transport_integrals                                     &
                               (ntrace,            icellsng (:,iblk),    &
                                indxsng(:,:),                        &
                                tracer_type,       depend,               &
                                integral_order,    triarea,              &
                                isflux,                                    &
                                xp,                yp,                   &
                                mc(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,n),           mx(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,n),              &
                                my(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,n),  mflxn (:,  n),      &
                                tc(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:,n),         tx(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:,n),            &
                                ty(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:,n),  mtflxn(:,:,n))

         enddo                  ! n

         !-------------------------------------------------------------------
         ! Update the ice area and tracers.
         !-------------------------------------------------------------------

         ! open water
         call update_fields (dynCnt(iblk),                                   &
                             skipPhys(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                             dynNbr  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk),:), &
                             dynIJ   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                             nx_block,                              &
                             ntrace,                                &
                             tracer_type,        depend,            &
                             tarear_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                             l_stop,            &
                             istop,              jstop,             &
                             mflxe(:,0),       mflxn(:,0),      &
                             mm(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,0))

         if (l_stop) then
            call diagnostic_abort(istop,jstop,iblk,'negative area (open water)')
         endif

         ! ice categories
         do n = 1, ncat

            call update_fields(dynCnt(iblk),                                   &
                               skipPhys(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                               dynNbr  (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk),:), &
                               dynIJ   (dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)), &
                               nx_block,                               &
                               ntrace,                                 &
                               tracer_type,          depend,           &
                               tarear_slot(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1), &
                               l_stop,                                 &
                               istop,                jstop,            &
                               mflxe (:,  n),      mflxn (:,  n),  &
                               mm(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,n),      &
                               mtflxe(:,:,n),      mtflxn(:,:,n),  &
                               tm(dynOff(iblk)+1:dynOff(iblk)+dynCnt(iblk)+1,:,n))

            if (l_stop) then
               write (nu_diag,*) 'istep1, my_task, iblk, cat =',     &
                                  istep1, my_task, iblk, n
               call diagnostic_abort(istop,jstop,iblk,'negative area (ice)')
            endif
         enddo                  ! n

         deallocate(edgearea_e, edgearea_n, mflxe, mflxn, mtflxe, mtflxn, &
                    triarea, xp, yp, iflux, jflux, isflux, indxsng)

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      end subroutine horizontal_remap

!=======================================================================
!
! Make area and tracer masks.
!
! If an area is masked out (mm < puny), then the values of tracers
!  in that grid cell are assumed to have no physical meaning.
!
! Similarly, if a tracer with dependents is masked out
!  (abs(tm) < puny), then the values of its dependent tracers in that
!  grid cell are assumed to have no physical meaning.
! For example, the enthalpy value has no meaning if the thickness
!  is zero.
!
! author William H. Lipscomb, LANL

      subroutine make_masks (ncells,                   &
                             ntrace,                   &
                             has_dependents,           &
                             mm,             mmask,    &
                             tm,             tmask)

      ! The index lists this used to return -- icells, indxi, indxj -- are
      ! gone.  construct_fields was their only consumer and it now tests the
      ! mmask per slot instead, so they were being built, redefined over a
      ! second range, and never read.
      !
      ! Every loop here ran over the WHOLE block, 1..nx_block x 1..ny_block,
      ! so the slot form needs no range mask: a slot loop covers exactly the
      ! cells the block loop did, minus the ones with no slot, and nothing
      ! reads a mask outside mask_any.

      integer (kind=int_kind), intent(in) ::     &
         ncells            , & ! number of slots this block holds
         ntrace                ! number of tracers in use

      logical (kind=log_kind), dimension (ntrace), intent(in) ::     &
         has_dependents        ! true if a tracer has dependent tracers

      ! This block's slice of the compressed layout, category 0 upwards.
      real (kind=dbl_kind), dimension (:,0:), intent(in) :: &
         mm                    ! mean ice area in each grid cell

      real (kind=dbl_kind), dimension (:,0:), intent(out) :: &
         mmask                 ! = 1. if ice is present, else = 0.

      real (kind=dbl_kind), dimension (:,:,:), intent(in), optional :: &
         tm                    ! mean tracer values in each grid cell

      real (kind=dbl_kind), dimension (:,:,:), intent(out), optional :: &
         tmask                 ! = 1. if tracer is present, else = 0.

      ! local variables

      integer (kind=int_kind) ::     &
         ms                , & ! slot index
         n                 , & ! ice category index
         nt                    ! tracer index

      real (kind=dbl_kind) :: &
         puny                  !

      character(len=*), parameter :: subname = '(make_masks)'

      call icepack_query_parameters(puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      !-------------------------------------------------------------------
      ! open water mask
      !-------------------------------------------------------------------

      do ms = 1, ncells
         if (mm(ms,0) > puny) then
            mmask(ms,0) = c1
         else
            mmask(ms,0) = c0
         endif
      enddo

      do n = 1, ncat

         !-------------------------------------------------------------------
         ! ice area mask: grid cells where ice is present
         !-------------------------------------------------------------------

         mmask(:,n) = c0
         do ms = 1, ncells
            if (mm(ms,n) > puny) mmask(ms,n) = c1
         enddo

         !-------------------------------------------------------------------
         ! tracer masks
         !-------------------------------------------------------------------

         if (present(tm)) then

            tmask(:,:,n) = c0

            ! Same cells the ice-area mask just marked, which is what the
            ! icells list carried: mm > puny over the WHOLE block, before it
            ! was redefined over the narrower range below.
            do nt = 1, ntrace
               if (has_dependents(nt)) then
                  do ms = 1, ncells
                     if (mmask(ms,n) > puny .and. abs(tm(ms,nt,n)) > puny) then
                        tmask(ms,nt,n) = c1
                     endif
                  enddo
               endif
            enddo

         endif                     ! present(tm)

      enddo ! n

      end subroutine make_masks

!=======================================================================
!
! Construct fields of ice area and tracers.
!
! authors William H. Lipscomb, LANL
!         John R. Baumgardner, LANL

      subroutine construct_fields (iblk,           ntrace,     &
                                   tracer_type,    depend,     &
                                   has_dependents,             &
                                   hm,                         &
                                   mm,             mc,         &
                                   mx,             my,         &
                                   mmask,                      &
                                   tm,             tc,         &
                                   tx,             ty,         &
                                   tmask)

      use ice_grid, only: dynMap, dynOff, dynCnt, dynNbr, skipRmp, &
                          pack_to_slots

      integer (kind=int_kind), intent(in) ::   &
         iblk              , & ! block index, to slice the compressed layout
         ntrace                ! number of tracers in use

      integer (kind=int_kind), dimension (ntrace), intent(in) ::     &
         tracer_type       , & ! = 1, 2, or 3 (see comments above)
         depend                ! tracer dependencies (see above)

      logical (kind=log_kind), dimension (ntrace), intent(in) ::     &
         has_dependents        ! true if a tracer has dependent tracers

      ! Assumed shape, so this routine does not name the block dimensions.
      ! ⚠ hm is rank 1: the caller passes hm(:,iblk), which is already
      ! flattened.  The others arrive as genuine rank-2/rank-3 sections.
      ! contiguous because they are handed on to pack_to_slots, whose dummy
      ! is explicit-shape by design -- see ice_grid.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous ::   &
         hm                    ! land/boundary mask, thickness (T-cell)

      ! Compressed: this block's slice, one category.  They arrive already
      ! in the slot layout, so the packs this routine used to do are gone.
      real (kind=dbl_kind), dimension (:), intent(in), contiguous ::   &
         mm                , & ! mean value of mass field
         mmask                 ! = 1. if ice is present, = 0. otherwise

      real (kind=dbl_kind), dimension (:,:), intent(in), optional, contiguous ::   &
         tm                , & ! mean tracer            (slot, tracer)
         tmask                 ! = 1. if tracer is present, = 0. otherwise

      ! Compressed, written straight out.
      real (kind=dbl_kind), dimension (:), intent(out), contiguous ::   &
         mc                , & ! mass value at geometric center of cell
         mx, my                ! limited derivative of mass wrt x and y

      real (kind=dbl_kind), dimension (:,:), intent(out), optional, contiguous ::   &
         tc                , & ! tracer at geometric center of cell
         tx, ty                ! limited derivative of tracer wrt x and y

      ! local variables

      integer (kind=int_kind) ::   &
         nt, nt1               ! tracer indices


      real (kind=dbl_kind), parameter :: xxav=c1/c12 ! mean T-cell values of xx
      real (kind=dbl_kind), parameter :: yyav=c1/c12 ! mean T-cell values of yy

      ! (slot, tracer), this block only -- so the same size the block form
      ! carried, not dynTot x ntrace.  Read at a DIFFERENT tracer than the
      ! one being computed (nt1 = depend(nt)), which is why they are rank 2
      ! while the rest of the working set is rank 1.
      real (kind=dbl_kind), dimension (dynCnt(iblk)+1,ntrace) ::  &
         mtxav             , & ! x coordinate of center of mass*tracer
         mtyav                 ! y coordinate of center of mass*tracer

      real (kind=dbl_kind) ::   &
         puny, &
         w2, w3, w7        ! work variables

      integer (kind=int_kind) :: ncells, is

      ! The compressed working set.  mm, mmask, tm and tmask are no longer
      ! here: they arrive compressed, so only the genuinely local vectors
      ! remain.
      real (kind=dbl_kind), dimension (dynCnt(iblk)+1) :: &
         hm_s, zero_s, &
         mxav_s, myav_s, &
         tc_s, tx_s, ty_s

      character(len=*), parameter :: subname = '(construct_fields)'

      !-------------------------------------------------------------------
      ! Compute field values at the geometric center of each grid cell,
      ! and compute limited gradients in the x and y directions.
      !
      ! For second order accuracy, each state variable is approximated as
      ! a field varying linearly over x and y within each cell.  For each
      ! category, the integrated value of m(x,y) over the cell must
      ! equal mm(i,j,n)*tarea(i,j), where tarea(i,j) is the cell area.
      ! Similarly, the integrated value of m(x,y)*t(x,y) must equal
      ! the total mass*tracer, mm(i,j,n)*tm(i,j,n)*tarea(i,j).
      !
      ! These integral conditions are satisfied for linear fields if we
      ! stipulate the following:
      ! (1) The mean mass, mm, is equal to the mass at the cell centroid.
      ! (2) The mean value tm1 of type 1 tracers is equal to the value
      !     at the center of mass.
      ! (3) The mean value tm2 of type 2 tracers is equal to the value
      !     at the center of mass*tm1, where tm2 depends on tm1.
      !     (See comments at the top of the module.)
      !
      ! We want to find the value of each state variable at a standard
      ! reference point, which we choose to be the geometric center of
      ! the cell.  The geometric center is located at the intersection
      ! of the line joining the midpoints of the north and south edges
      ! with the line joining the midpoints of the east and west edges.
      ! To find the value at the geometric center, we must know the
      ! location of the cell centroid/center of mass, along with the
      ! mean value and the gradients with respect to x and y.
      !
      ! The cell gradients are first computed from the difference between
      ! values in the neighboring cells, then limited by requiring that
      ! no new extrema are created within the cell.
      !
      ! For rectangular coordinates the centroid and the geometric
      ! center coincide, which means that some of the equations in this
      ! subroutine could be simplified.  However, the full equations
      ! are retained for generality.
      !-------------------------------------------------------------------

      !-------------------------------------------------------------------
      ! Initialize
      !-------------------------------------------------------------------

      call icepack_query_parameters(puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      mc(:) = c0
      mx(:) = c0
      my(:) = c0

      if (present(tm)) then
         tc(:,:) = c0
         tx(:,:) = c0
         ty(:,:) = c0
      endif

      ! limited gradient of mass field in each cell (except masked cells)
      ! Note: The gradient is computed in scaled coordinates with
      !       dxT = dyT = hte = htn = 1.

      ncells = dynCnt(iblk)

      ! xav and yav are never assigned -- they are c0 for the whole routine,
      ! which is what the "echmod: grid arrays = 0" note below refers to --
      ! so they are passed as a zero vector.
      call pack_to_slots(hm   , dynMap(:,iblk), hm_s)
      zero_s(:) = c0

      call limited_gradient (ncells,                                            &
                             skipRmp(dynOff(iblk)+1:dynOff(iblk)+ncells),       &
                             dynNbr (dynOff(iblk)+1:dynOff(iblk)+ncells,:),     &
                             mm,       hm_s,                                    &
                             zero_s,   zero_s,                                  &
                             mx,     my)

      ! ⚠ mmask is the indicator of the icells list: make_masks sets it to c1
      ! at exactly the cells it counted, mm > puny.  So testing it per slot
      ! selects the same cells the indxi/indxj loop did, and every such cell
      ! has a slot -- ice implies ocean implies tmask implies mask_any.
      mc(:) = c0
      do is = 1, ncells
         if (mmask(is) <= puny) cycle

         ! mass field at geometric center
         mc(is) = mm(is)

      enddo                     ! is



      ! tracers

      if (present(tm)) then

         mxav_s(:) = c0
         myav_s(:) = c0
         do is = 1, ncells       ! cells with mass
            if (mmask(is) <= puny) cycle

            ! center of mass (mxav,myav) for each cell

            mxav_s(is) = mx(is)*xxav / mm(is)
            myav_s(is) = my(is)*yyav / mm(is)

         enddo

         do nt = 1, ntrace

            ! The block form relied on tc/tx/ty having been zeroed for every
            ! tracer up front; a slot vector is reused across nt, so it has
            ! to be cleared here or a tracer_type the chain does not handle
            ! would inherit the previous tracer's values.
            tc_s(:) = c0
            tx_s(:) = c0
            ty_s(:) = c0

            if (tracer_type(nt)==1) then ! independent of other tracers

               call limited_gradient(ncells,                                      &
                                     skipRmp(dynOff(iblk)+1:dynOff(iblk)+ncells), &
                                     dynNbr (dynOff(iblk)+1:dynOff(iblk)+ncells,:),&
                                     tm(:,nt), mmask,                             &
                                     mxav_s,   myav_s,                            &
                                     tx_s,     ty_s)

               if (has_dependents(nt)) then   ! need center of area*tracer

                  mtxav(:,nt) = c0
                  mtyav(:,nt) = c0

                  tc_s(:) = c0
                  do is = 1, ncells  ! Note: no tx or ty in ghost cells
                                     ! (bound calls are later)
                     if (mmask(is) <= puny) cycle

                     ! tracer value at geometric center
                     tc_s(is) = tm(is,nt) - tx_s(is)*mxav_s(is)   &
                                          - ty_s(is)*myav_s(is)

                     if (tmask(is,nt) > puny) then

                        ! center of area*tracer
                        w2 = mc(is)*tx_s(is)   &
                           + mx(is)*tc_s(is)
                        w3 = mc(is)*ty_s(is)   &
                           + my(is)*tc_s(is)
                        w7 = c1 / (mm(is)*tm(is,nt))
                        ! echmod: grid arrays = 0
                        mtxav(is,nt) = w2*xxav *w7
                        mtyav(is,nt) = w3*yyav * w7
                     endif         ! tmask

                  enddo            ! is

               else                ! no dependents

                  tc_s(:) = c0
                  do is = 1, ncells      ! mass is present
                     if (mmask(is) <= puny) cycle

                     ! tracer value at geometric center
                     tc_s(is) = tm(is,nt) - tx_s(is)*mxav_s(is)   &
                                          - ty_s(is)*myav_s(is)
                  enddo            ! is

               endif               ! has_dependents

            elseif (tracer_type(nt)==2) then   ! tracer nt depends on nt1
               nt1 = depend(nt)

               call limited_gradient(ncells,                                      &
                                     skipRmp(dynOff(iblk)+1:dynOff(iblk)+ncells), &
                                     dynNbr (dynOff(iblk)+1:dynOff(iblk)+ncells,:),&
                                     tm(:,nt), tmask(:,nt1),                      &
                                     mtxav(:,nt1), mtyav(:,nt1),                  &
                                     tx_s,     ty_s)

               tc_s(:) = c0
               do is = 1, ncells     ! ice is present
                  if (mmask(is) <= puny) cycle
                  tc_s(is) = tm(is,nt)                       &
                           - tx_s(is) * mtxav(is,nt1)        &
                           - ty_s(is) * mtyav(is,nt1)
               enddo               ! is

            elseif (tracer_type(nt)==3) then  ! upwind approx; gradient = 0

               tc_s(:) = c0
               do is = 1, ncells
                  if (mmask(is) <= puny) cycle

                  tc_s(is) = tm(is,nt)
               enddo               ! is

            endif                  ! tracer_type

            tc(:,nt) = tc_s(:)
            tx(:,nt) = tx_s(:)
            ty(:,nt) = ty_s(:)

         enddo                    ! ntrace

      endif                     ! present (tm)

      end subroutine construct_fields

!=======================================================================
! Compute a limited gradient of the scalar field phi in scaled coordinates.
! "Limited" means that we do not create new extrema in phi.  For
! instance, field values at the cell corners can neither exceed the
! maximum of phi(i,j) in the cell and its eight neighbors, nor fall
! below the minimum.
!
! authors William H. Lipscomb, LANL
!         John R. Baumgardner, LANL

      subroutine limited_gradient (ncells,   skipme,   &
                                   nbr,                &
                                   phi,      phimask,  &
                                   cnx,      cny,      &
                                   gx,       gy)

      use ice_grid, only: nbr_W, nbr_E, nbr_N, nbr_S, &
                          nbr_NW, nbr_NE, nbr_SW, nbr_SE

      integer (kind=int_kind), intent(in) :: &
         ncells                ! slots held by this block

      ! .true. outside the remap loop range.  STATIC block geometry --
      ! skipRmp in ice_grid -- not a statement about where the ice is.  The
      ! block form ran ilo-nghost+1..ihi+nghost-1, so at nghost = 1 it
      ! computed the physical cells only and a halo exchange filled the
      ! ghosts; iterating every slot would compute where that left c0.
      logical (kind=log_kind), dimension(:), intent(in), contiguous :: &
         skipme

      ! (slot, direction) -> neighbour slot; a missing neighbour points at
      ! the block's dead slot, which holds zero.
      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nbr

      ! compressed: this block's slice, 1-based, cells 1..ncells, dead slot last
      real (kind=dbl_kind), dimension (:), intent (in), contiguous :: &
         phi      , & ! input tracer field (mean values in each grid cell)
         phimask  , & ! phimask = 1 if phi has physical meaning, = 0 otherwise.
                      ! For instance, aice has no physical meaning in land cells,
                      ! and hice no physical meaning where aice = 0.
         cnx      , & ! x-coordinate of phi relative to geometric center of cell
         cny          ! y-coordinate of phi relative to geometric center of cell

      real (kind=dbl_kind), dimension (:), intent(out), contiguous :: &
         gx       , & ! limited x-direction gradient
         gy           ! limited y-direction gradient

      ! local variables

      integer (kind=int_kind) :: &
         is                    , & ! slot counter
         isW, isE, isN, isS    , & ! the eight neighbour slots
         isNW, isNE, isSW, isSE

      real (kind=dbl_kind) ::   &
         phi_nw, phi_n, phi_ne , & ! values of phi in 8 neighbor cells
         phi_w,         phi_e  , &
         phi_sw, phi_s, phi_se , &
         qmn, qmx              , & ! min and max value of phi within grid cell
         pmn, pmx              , & ! min and max value of phi among neighbor cells
         w1, w2, w3, w4            ! work variables

      real (kind=dbl_kind) ::   &
         puny        , & !
         gxtmp, gytmp    ! temporary term for x- and y- limited gradient

      character(len=*), parameter :: subname = '(limited_gradient)'

      call icepack_query_parameters(puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      gx(:) = c0
      gy(:) = c0

      do is = 1, ncells
         if (skipme(is)) cycle
         if (phimask(is) <= puny) cycle

         isW  = nbr(is,nbr_W )
         isE  = nbr(is,nbr_E )
         isN  = nbr(is,nbr_N )
         isS  = nbr(is,nbr_S )
         isNW = nbr(is,nbr_NW)
         isNE = nbr(is,nbr_NE)
         isSW = nbr(is,nbr_SW)
         isSE = nbr(is,nbr_SE)

         ! Store values of phi in the 8 neighbor cells.
         ! Note: phimask = 1. or 0.  If phimask = 1., use the true value;
         !  if phimask = 0., use the home cell value so that non-physical
         !  values of phi do not contribute to the gradient.
         phi_nw = phimask(isNW) * phi(isNW)  &
            + (c1-phimask(isNW))* phi(is  )
         phi_n  = phimask(isN ) * phi(isN )  &
            + (c1-phimask(isN ))* phi(is  )
         phi_ne = phimask(isNE) * phi(isNE)  &
            + (c1-phimask(isNE))* phi(is  )
         phi_w  = phimask(isW ) * phi(isW )  &
            + (c1-phimask(isW ))* phi(is  )
         phi_e  = phimask(isE ) * phi(isE )  &
            + (c1-phimask(isE ))* phi(is  )
         phi_sw = phimask(isSW) * phi(isSW)  &
            + (c1-phimask(isSW))* phi(is  )
         phi_s  = phimask(isS ) * phi(isS )  &
            + (c1-phimask(isS ))* phi(is  )
         phi_se = phimask(isSE) * phi(isSE)  &
            + (c1-phimask(isSE))* phi(is  )

         ! unlimited gradient components
         ! (factors of two cancel out)

         gxtmp = (phi_e - phi_w) * p5
         gytmp = (phi_n - phi_s) * p5

         ! minimum and maximum among the nine local cells
         pmn = min (phi_nw, phi_n,  phi_ne, phi_w, phi(is),   &
                    phi_e,  phi_sw, phi_s,  phi_se)
         pmx = max (phi_nw, phi_n,  phi_ne, phi_w, phi(is),   &
                    phi_e,  phi_sw, phi_s,  phi_se)

         pmn = pmn - phi(is)
         pmx = pmx - phi(is)

         ! minimum and maximum deviation of phi within the cell
         w1  =  (p5 - cnx(is)) * gxtmp   &
              + (p5 - cny(is)) * gytmp
         w2  =  (p5 - cnx(is)) * gxtmp   &
              - (p5 + cny(is)) * gytmp
         w3  = -(p5 + cnx(is)) * gxtmp   &
              - (p5 + cny(is)) * gytmp
         w4  =  (p5 - cny(is)) * gytmp   &
              - (p5 + cnx(is)) * gxtmp

         qmn = min (w1, w2, w3, w4)
         qmx = max (w1, w2, w3, w4)

         ! the limiting coefficient
         if ( abs(qmn) > abs(pmn) ) then ! 'abs(qmn) > puny' not sufficient
            w1 = max(c0, pmn/qmn)
         else
            w1 = c1
         endif

         if ( abs(qmx) > abs(pmx) ) then
            w2 = max(c0, pmx/qmx)
         else
            w2 = c1
         endif

         w1 = min(w1, w2)

         ! Limit the gradient components
         gx(is) = w1 * gxtmp
         gy(is) = w1 * gytmp

      enddo                     ! is

      end subroutine limited_gradient

!=======================================================================
!
! Given velocity fields on cell corners, compute departure points
! of back trajectories in nondimensional coordinates.
!
! author William H. Lipscomb, LANL

      subroutine departure_points (ncells,     skipme,   &
                                   nbr,        ijslot,   &
                                   nx_block,             &
                                   dt,                   &
                                   uvel,       vvel,     &
                                   dxu,        dyu,      &
                                   HTN,        HTE,      &
                                   HTNe,       HTEn,     &
                                   dpx,        dpy,      &
                                   l_dp_midpt, l_stop,   &
                                   istop,      jstop)

      use ice_grid, only: nbr_W, nbr_E, nbr_N, nbr_S, &
                          nbr_NW, nbr_NE, nbr_SW, nbr_SE

      integer (kind=int_kind), intent(in) ::   &
         ncells            , & ! slots held by this block
         nx_block              ! only to turn ijslot into i and j for the abort

      ! .true. outside the remap loop range -- skipRmp, static geometry.
      logical (kind=log_kind), dimension(:), intent(in), contiguous ::   &
         skipme

      ! (slot, direction) -> neighbour slot; a missing neighbour points at
      ! the block's dead slot, which holds zero.
      integer (kind=int_kind), dimension (:,:), intent(in) ::   &
         nbr

      ! (j-1)*nx_block+i of the cell owning each slot, for the abort message
      integer (kind=int_kind), dimension (:), intent(in) ::   &
         ijslot

      real (kind=dbl_kind), intent(in) ::   &
         dt               ! time step (s)

      ! compressed: this block's slice, 1-based, cells 1..ncells, dead slot last
      real (kind=dbl_kind), dimension (:), intent(in), contiguous ::   &
         uvel         , & ! x-component of velocity (m/s)
         vvel         , & ! y-component of velocity (m/s)
         dxu          , & ! E-W dimensions of U-cell (m)
         dyu          , & ! N-S dimensions of U-cell (m)
         HTN          , & ! length of north face of T-cell (m)
         HTE          , & ! length of east face of T-cell (m)
         HTNe         , & ! HTN one cell EAST, precomputed in block space
         HTEn             ! HTE one cell NORTH, precomputed in block space

      ! Compressed, and written straight out: nothing between here and
      ! locate_triangles needs a block form any more, the halo exchange
      ! included.
      real (kind=dbl_kind), dimension (:), intent(out), contiguous ::   &
         dpx          , & ! coordinates of departure points (m)
         dpy              ! coordinates of departure points (m)

      logical (kind=log_kind), intent(in) ::   &
         l_dp_midpt       ! if true, find departure points using
                          ! corrected midpoint velocity

      logical (kind=log_kind), intent(inout) ::   &
         l_stop           ! if true, abort on return

      integer (kind=int_kind), intent(inout) ::   &
         istop, jstop     ! indices of grid cell where model aborts

      ! local variables

      integer (kind=int_kind) ::   &
         is           , & ! slot counter
         isstop       , & ! slot where the departure point went out of bounds
         ijb          , & ! that slot's cell, to report i and j
         isa, isb     , & ! the four slots the bilinear interpolation reads,
         isc, isd         ! in the order (i2-1,j2-1) (i2,j2-1) (i2,j2) (i2-1,j2)

      real (kind=dbl_kind) ::                  &
         mpx,  mpy    , & ! coordinates of midpoint of back trajectory,
                          ! relative to cell corner
         mpxt, mpyt   , & ! midpoint coordinates relative to cell center
         ump,  vmp        ! corrected velocity at midpoint

      character(len=*), parameter :: subname = '(departure_points)'

      !-------------------------------------------------------------------
      ! Estimate departure points.
      ! This estimate is 1st-order accurate in time; improve accuracy by
      !  using midpoint approximation (to add later).
      ! The loop range is carried by skipme: for nghost = 1 that is the
      ! physical cells, with the ghost values filled by a halo exchange
      ! afterwards; for nghost = 2 it reaches a layer into the ghosts.
      !-------------------------------------------------------------------

      dpx(:) = c0
      dpy(:) = c0
      isstop = 0

      do is = 1, ncells
         if (skipme(is)) cycle

         dpx(is) = -dt*uvel(is)
         dpy(is) = -dt*vvel(is)

         ! Check for values out of bounds (more than one grid cell away).
         !
         ! ⚠ HTNe and HTEn are the shifted lengths, precomputed in block
         ! space.  Reading HTN at nbr(is,nbr_E) would give c0 wherever that
         ! neighbour has no slot, and every positive dpx would then look out
         ! of bounds -- HTN is an edge length, nonzero over land.
         if ((dpx(is) > c0 .and. (dpx(is) < -HTN(is) .or. dpx(is) > HTNe(is))) .or.   &
             (dpy(is) > c0 .and. (dpy(is) < -HTE(is) .or. dpy(is) > HTEn(is)))) then
            l_stop = .true.
            isstop = is
         endif

      enddo

      if (l_stop) then
         ijb   = ijslot(isstop)
         istop = mod(ijb-1, nx_block) + 1
         jstop = (ijb-1) / nx_block + 1
         write (nu_diag,*) ' '
         write (nu_diag,*)   &
                    'Warning: Departure points out of bounds in remap'
         write (nu_diag,*) 'my_task, i, j =', my_task, istop, jstop
         write (nu_diag,*) 'dpx, dpy =', dpx(isstop), dpy(isstop)
         write (nu_diag,*) 'HTN(i,j), HTN(i+1,j) =', HTN(isstop), HTNe(isstop)
         write (nu_diag,*) 'HTE(i,j), HTE(i,j+1) =', HTE(isstop), HTEn(isstop)
         return
      endif

      if (l_dp_midpt) then ! find dep pts using corrected midpt velocity

         do is = 1, ncells
            if (skipme(is)) cycle
            if (uvel(is)/=c0 .or. vvel(is)/=c0) then

               !-------------------------------------------------------------------
               ! Scale departure points to coordinate system in which grid cells
               ! have sides of unit length.
               !-------------------------------------------------------------------

               dpx(is) = dpx(is) / dxu(is)
               dpy(is) = dpy(is) / dyu(is)

               !-------------------------------------------------------------------
               ! Estimate midpoint of backward trajectory relative to corner (i,j).
               !-------------------------------------------------------------------

               mpx = p5 * dpx(is)
               mpy = p5 * dpy(is)

               !-------------------------------------------------------------------
               ! Determine which cell the trajectory lies in, and the midpoint
               !  coordinates relative to that cell's center in a stretch
               !  coordinate system with vertices at (1/2, 1/2), (1/2, -1/2), etc.
               !
               ! The block form picked a cell (i2,j2) and then read the four
               ! cells (i2-1,j2-1) (i2,j2-1) (i2,j2) (i2-1,j2).  Every one of
               ! the sixteen combinations is self or one of the eight
               ! neighbours, so the selection is a choice of slots -- the
               ! arithmetic that follows is unchanged.
               !-------------------------------------------------------------------

               if (mpx >= c0 .and. mpy >= c0) then    ! cell (i+1,j+1)
                  isa  = is
                  isb = nbr(is,nbr_E )
                  isc  = nbr(is,nbr_NE); isd = nbr(is,nbr_N )
                  mpxt = mpx - p5
                  mpyt = mpy - p5
               elseif (mpx < c0 .and. mpy < c0) then  ! cell (i,j)
                  isa  = nbr(is,nbr_SW); isb = nbr(is,nbr_S )
                  isc  = is
                  isd = nbr(is,nbr_W )
                  mpxt = mpx + p5
                  mpyt = mpy + p5
               elseif (mpx >= c0 .and. mpy < c0) then ! cell (i+1,j)
                  isa  = nbr(is,nbr_S ); isb = nbr(is,nbr_SE)
                  isc  = nbr(is,nbr_E ); isd = is
                  mpxt = mpx - p5
                  mpyt = mpy + p5
               elseif (mpx < c0 .and. mpy >= c0) then ! cell (i,j+1)
                  isa  = nbr(is,nbr_W ); isb = is
                  isc  = nbr(is,nbr_N ); isd = nbr(is,nbr_NW)
                  mpxt = mpx + p5
                  mpyt = mpy - p5
               endif

               !-------------------------------------------------------------------
               ! Using a bilinear approximation, estimate the velocity at the
               ! trajectory midpoint in the (i2,j2) reference frame.
               !
               ! A neighbour with no slot reads the dead slot, zero, which is
               ! what the block layout carries for a velocity at an unmapped
               ! cell.  Velocities are state, not geometry, so unlike HTN
               ! above these need no shifted copy.
               !-------------------------------------------------------------------

               ump = uvel(isa)*(mpxt-p5)*(mpyt-p5)     &
                   - uvel(isb)*(mpxt+p5)*(mpyt-p5)     &
                   + uvel(isc)*(mpxt+p5)*(mpyt+p5)     &
                   - uvel(isd)*(mpxt-p5)*(mpyt+p5)

               vmp = vvel(isa)*(mpxt-p5)*(mpyt-p5)     &
                   - vvel(isb)*(mpxt+p5)*(mpyt-p5)     &
                   + vvel(isc)*(mpxt+p5)*(mpyt+p5)     &
                   - vvel(isd)*(mpxt-p5)*(mpyt+p5)

               !-------------------------------------------------------------------
               ! Use the midpoint velocity to estimate the coordinates of the
               !  departure point relative to corner (i,j).
               !-------------------------------------------------------------------

               dpx(is) = -dt * ump
               dpy(is) = -dt * vmp

            endif               ! nonzero velocity

         enddo                 ! is

      endif                  ! l_dp_midpt

      end subroutine departure_points

!=======================================================================
!
! Compute areas and vertices of transport triangles for north or
!  east cell edges.
!
! authors William H. Lipscomb, LANL
!         John R. Baumgardner, LANL

      subroutine locate_triangles (nx_block,     ny_block,   &
                                   ilo, ihi,     jlo, jhi,   &
                                   nghost,       edge,       &
                                   icells,                   &
                                   indxs,                    &
                                   dpx,          dpy,  map,  &
                                   ijslot,                   &
                                   dxu,          dyu,        &
                                   earea,        narea,      &
                                   xp,           yp,         &
                                   iflux,        jflux,      &
                                   isflux,                   &
                                   triarea,      edgearea)

      integer (kind=int_kind), intent(in) ::   &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi   , & ! beginning and end of physical domain
         nghost                ! number of ghost cells

      character (len=char_len), intent(in) ::   &
         edge             ! 'north' or 'east'

      ! dpx and dpy arrive COMPRESSED -- this block's slice, 1-based -- and
      ! map turns a cell into its slot.  The loops below keep their (j,i)
      ! order deliberately: the icellsd list they build is walked again by
      ! transport_integrals, and slot order is not (j,i) order, so iterating
      ! slots here would reorder floating-point sums downstream.
      real (kind=dbl_kind), dimension(:), intent(in), contiguous ::  &
         dpx          , & ! x coordinates of departure points at cell corners
         dpy              ! y coordinates of departure points at cell corners

      integer (kind=int_kind), dimension(:), intent(in) ::  &
         map              ! this block's dynMap column

      ! ⚠ The reverse of map: the flattened (i,j) each slot sits at.  Only
      ! the vertex shift below needs it -- ishift/jshift are the donor's
      ! (i,j) minus this cell's -- and deriving i and j from it costs two
      ! integer operations against the two array loads the indxi/indxj
      ! lists used to cost.
      integer (kind=int_kind), dimension (:), intent(in) ::   &
         ijslot           ! flattened (i,j) of each slot

      ! dxu and dyu are read only to scale the departure points, which is
      ! now a slot loop, so they arrive compressed.
      real (kind=dbl_kind), dimension(:), intent(in), contiguous ::  &
         dxu          , & ! E-W dimension of U-cell (m)
         dyu              ! N-S dimension of U-cell (m)

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) ::  &
         earea        , & ! area of E-cell
         narea            ! area of N-cell

      ! Slot-indexed: (slot, vertex, group).  The triangle geometry is
      ! per cell and these are only ever read at the cell that owns them,
      ! so the block layout bought nothing.
      real (kind=dbl_kind), dimension (:,0:,:), intent(out) :: &
         xp, yp           ! coordinates of triangle vertices

      real (kind=dbl_kind), dimension (:,:), intent(out) :: &
         triarea          ! area of departure triangle

      ! ⚠ Slot-indexed storage, but they still HOLD (i,j) of the donor cell.
      ! mc, mx, my, tc, tx and ty are read there and are still block arrays;
      ! turning these into donor SLOTS is a separate change that converts
      ! those six and nine halo calls with them.
      integer (kind=int_kind), dimension (:,:), intent(out) :: &
         iflux        , & ! i index of cell contributing transport
         jflux        , & ! j index of cell contributing transport
         isflux           ! and that donor cell's SLOT

      integer (kind=int_kind), dimension (ngroups), intent(out) ::   &
         icells           ! number of cells where triarea > puny

      integer (kind=int_kind), dimension (:,:), intent(out) :: &
         indxs            ! the slot each of those cells owns

      ! Compressed.  Read at this cell only, so a cell with no slot reads the
      ! dead slot -- c0, which is what the block form held there: edgearea is
      ! zeroed wholesale and then filled from uvelE/vvelN, which are masked
      ! to zero wherever a cell is outside mask_any.
      real (kind=dbl_kind), dimension(:), intent(inout) ::   &
         edgearea         ! area of departure region for each edge
                          ! edgearea > 0 for eastward/northward flow

      ! local variables

      integer (kind=int_kind) ::   &
         i, j, ij, ic         , & ! horizontal indices
         ijb                  , & ! flattened (i,j) of a slot
         is                   , & ! slot counter
         ms, mw               , & ! slot of this cell and of its W or S neighbour
         mdead                , & ! this block's dead slot, where map gives 0
         ib, jb               , & ! limits for loops for bugcheck
         ng, nv               , & ! triangle indices
         ishift   , jshift    , & ! differences between neighbor cells
         ishift_tl, jshift_tl , & ! i,j indices of TL cell relative to edge
         ishift_bl, jshift_bl , & ! i,j indices of BL cell relative to edge
         ishift_tr, jshift_tr , & ! i,j indices of TR cell relative to edge
         ishift_br, jshift_br , & ! i,j indices of BR cell relative to edge
         ishift_tc, jshift_tc , & ! i,j indices of TC cell relative to edge
         ishift_bc, jshift_bc , & ! i,j indices of BC cell relative to edge
         is_l, js_l           , & ! i,j shifts for TL1, BL2 for area consistency
         is_r, js_r           , & ! i,j shifts for TR1, BR2 for area consistency
         ise_tl, jse_tl       , & ! i,j of TL other edge relative to edge
         ise_bl, jse_bl       , & ! i,j of BL other edge relative to edge
         ise_tr, jse_tr       , & ! i,j of TR other edge relative to edge
         ise_br, jse_br           ! i,j of BR other edge relative to edge

      integer (kind=int_kind) ::   &
         icellsd          ! number of cells where departure area > 0.

      integer (kind=int_kind), dimension (nx_block*ny_block) ::  &
         indxid       , & ! compressed index in i-direction
         indxjd           ! compressed index in j-direction

      ! The scaled departure points, compressed like dpx/dpy themselves.
      real (kind=dbl_kind), dimension(size(dpx)) ::   &
         dx, dy           ! scaled departure points

      ! ⚠ These two STAY BLOCK-SHAPED, deliberately, and the reason is
      ! measured rather than assumed.
      !
      ! They are earea/narea -- AREAS, so nonzero over land -- and they are
      ! used as multipliers.  Most reads are at the cell itself and would
      ! convert trivially, but twelve are at NAMED shifts:
      !
      !     areafact(ms,ng) = -areafac_ce(i+ise_tl, j+jse_tl)
      !
      ! Every shift is in {-1,0,+1} in both i and j, so each is a compass
      ! neighbour and dynNbr would cover them mechanically.  The catch is
      ! that dynNbr resolves a missing neighbour to the DEAD SLOT, which
      ! holds c0, while the block form reads the real area that land cell
      ! carries.
      !
      ! Instrumented on the rmap configuration (kmt_type='boxislands', so
      ! interior land surrounded by ocean) counting shifted reads whose cell
      ! has no slot AND whose areafac is nonzero -- exactly the divergent
      ! case:
      !
      !     tested 3,369,024    misses 33,696    (~1%)
      !
      ! So the hazard is REACHABLE, not theoretical.  Shifted slot copies
      ! would avoid it, but there are six distinct shifts across the two
      ! arrays per edge, and rebuilding six vectors per block per edge per
      ! step costs more bandwidth than the two block arrays it saves --
      ! these are only nx_block*ny_block each, trivial next to mtflxe/mtflxn.
      !
      ! Converting these buys almost nothing and risks silently wrong
      ! transported mass next to land.  Leave them alone.
      real (kind=dbl_kind), dimension(nx_block,ny_block) ::   &
         areafac_c    , & ! earea or narea
         areafac_ce       ! areafac_c on other edge (narea or earea)

      real (kind=dbl_kind) ::   &
         xcl, ycl     , & ! coordinates of left corner point
                          ! (relative to midpoint of edge)
         xdl, ydl     , & ! left departure point
         xil, yil     , & ! left intersection point
         xcr, ycr     , & ! right corner point
         xdr, ydr     , & ! right departure point
         xir, yir     , & ! right intersection point
         xic, yic     , & ! x-axis intersection point
         xicl, yicl   , & ! left-hand x-axis intersection point
         xicr, yicr   , & ! right-hand x-axis intersection point
         xdm, ydm     , & ! midpoint of segment connecting DL and DR;
                          ! shifted if l_edge_flux_adj = T
         md           , & ! slope of line connecting DL and DR
         mdl          , & ! slope of line connecting DL and DM
         mdr          , & ! slope of line connecting DR and DM
         area1, area2 , & ! temporary triangle areas
         area3, area4 , & !
         area_c       , & ! center polygon area
         puny         , & !
         w1, w2           ! work variables

      ! Slot-indexed.  Both are per-cell temporaries of the icellsd
      ! traversal, and every icellsd cell has a real slot -- dpx is nonzero
      ! only where departure_points ran -- so neither can alias the dead
      ! slot across iterations.
      real (kind=dbl_kind), dimension (size(dpx),ngroups) ::   &
         areafact         ! = 1 for positive flux, -1 for negative

      real (kind=dbl_kind), dimension(size(dpx)) ::   &
         areasum          ! sum of triangle areas for a given edge

      character(len=*), parameter :: subname = '(locate_triangles)'

      !-------------------------------------------------------------------
      ! Triangle notation:
      ! For each edge, there are 20 triangles that can contribute,
      ! but many of these are mutually exclusive.  It turns out that
      ! at most 5 triangles can contribute to transport integrals at once.
      !
      ! See Figure 3 in DB for pictures of these triangles.
      ! See Table 1 in DB for logical conditions.
      !
      ! For the north edge, DB refer to these triangles as:
      ! (1) NW, NW1, W, W2
      ! (2) NE, NE1, E, E2
      ! (3) NW2, W1, NE2, E1
      ! (4) H1a, H1b, N1a, N1b
      ! (5) H2a, H2b, N2a, N2b
      !
      ! For the east edge, DB refer to these triangles as:
      ! (1) NE, NE1, N, N2
      ! (2) SE, SE1, S, S2
      ! (3) NE2, N1, SE2, S1
      ! (4) H1a, H1b, E1a, E2b
      ! (5) H2a, H2b, E2a, E2b
      !
      ! The code below works for either north or east edges.
      ! The respective triangle labels are:
      ! (1) TL,  TL1, BL,  BL2
      ! (2) TR,  TR1, BR,  BR2
      ! (3) TL2, BL1, TR2, BR1
      ! (4) BC1a, BC1b, TC1a, TC2b
      ! (5) BC2a, BC2b, TC2a, TC2b
      !
      ! where the cell labels are:
      !
      !          |        |
      !     TL   |   TC   |   TR     (top left, center, right)
      !          |        |
      !   ------------------------
      !          |        |
      !     BL   |   BC   |   BR     (bottom left, center, right)
      !          |        |
      !
      ! and the transport is across the edge between cells TC and BC.
      !
      ! Departure points are scaled to a local coordinate system
      !  whose origin is at the midpoint of the edge.
      ! In this coordinate system, the lefthand corner CL = (-0.5,0)
      !  and the righthand corner CR = (0.5, 0).
      !-------------------------------------------------------------------

      !-------------------------------------------------------------------
      ! Initialize
      !-------------------------------------------------------------------

      call icepack_query_parameters(puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      areafac_c(:,:)  = c0
      areafac_ce(:,:) = c0

      triarea(:,:)   = c0
      xp     (:,:,:) = c0
      yp     (:,:,:) = c0

      ! iflux and jflux default to the cell itself, so they are seeded per
      ! cell rather than wholesale.  areafact is a slot vector now, so its
      ! zeroing is one assignment rather than a traversal.
      areafact(:,:) = c0

      do ng = 1, ngroups
         do j = 1, ny_block
         do i = 1, nx_block
            ms = map((j-1)*nx_block + i)
            if (ms > 0) then
               iflux(ms,ng) = i
               jflux(ms,ng) = j
            endif
         enddo
         enddo
      enddo

      if (trim(edge) == 'north') then

         ! index shifts for neighbor cells

         ishift_tl = -1
         jshift_tl =  1
         ishift_bl = -1
         jshift_bl =  0
         ishift_tr =  1
         jshift_tr =  1
         ishift_br =  1
         jshift_br =  0
         ishift_tc =  0
         jshift_tc =  1
         ishift_bc =  0
         jshift_bc =  0

         ! index shifts for TL1, BL2, TR1 and BR2 for area consistency

         is_l = -1
         js_l =  0
         is_r =  1
         js_r =  0

         ! index shifts for neighbor east edges

         ise_tl = -1
         jse_tl =  1
         ise_bl = -1
         jse_bl =  0
         ise_tr =  0
         jse_tr =  1
         ise_br =  0
         jse_br =  0

         ! area scale factor
         ! earea, narea valid on halo

         do j = 1, ny_block
         do i = 1, nx_block
            areafac_c(i,j) = narea(i,j)
         enddo
         enddo

         ! area scale factor for other edge (east)

         do j = 1, ny_block
         do i = 1, nx_block
            areafac_ce(i,j) = earea(i,j)
         enddo
         enddo

      else                      ! east edge

         ! index shifts for neighbor cells

         ishift_tl =  1
         jshift_tl =  1
         ishift_bl =  0
         jshift_bl =  1
         ishift_tr =  1
         jshift_tr = -1
         ishift_br =  0
         jshift_br = -1
         ishift_tc =  1
         jshift_tc =  0
         ishift_bc =  0
         jshift_bc =  0

         ! index shifts for TL1, BL2, TR1 and BR2 for area consistency

         is_l =  0
         js_l =  1
         is_r =  0
         js_r = -1

         ! index shifts for neighbor north edges

         ise_tl =  1
         jse_tl =  0
         ise_bl =  0
         jse_bl =  0
         ise_tr =  1
         jse_tr = -1
         ise_br =  0
         jse_br = -1

         ! area scale factors
         ! earea, narea valid on halo

         do j = 1, ny_block
         do i = 1, nx_block
            areafac_c(i,j) = earea(i,j)
         enddo
         enddo

         ! area scale factor for other edge (north)

         do j = 1, ny_block
         do i = 1, nx_block
            areafac_ce(i,j) = narea(i,j)
         enddo
         enddo

      endif

      !-------------------------------------------------------------------
      ! Compute mask for edges with nonzero departure areas and for
      ! one grid-cell wide channels
      !-------------------------------------------------------------------

      ! ⚠ dynMap returns 0 for a cell with no slot, where dynNbr would have
      ! returned the dead slot.  Reading a slot array at 0 is out of bounds,
      ! so resolve it here: the dead slot is the last element of the slice,
      ! and it holds c0 -- which is what dpx carries at an unmapped cell.
      mdead = size(dpx)

      icellsd = 0
      if (trim(edge) == 'north') then
         do j = jlo-1, jhi
         do i = ilo, ihi
            ms = map((j-1)*nx_block + i)
            if (ms < 1) ms = mdead
            mw = map((j-1)*nx_block + i-1)
            if (mw < 1) mw = mdead
            if (dpx(mw)/=c0 .or. dpy(mw)/=c0   &
                            .or.               &
                dpx(ms)/=c0 .or. dpy(ms)/=c0) then
               icellsd = icellsd + 1
               indxid(icellsd) = i
               indxjd(icellsd) = j
            else
               if ( abs(edgearea(ms)) > c0 ) then ! 1 grid-cell wide channel: dpx,y = 0, edgearea /= 0
                  icellsd = icellsd + 1
                  indxid(icellsd) = i
                  indxjd(icellsd) = j
               endif
            endif
         enddo
         enddo
      else       ! east edge
         do j = jlo, jhi
         do i = ilo-1, ihi
            ms = map((j-1)*nx_block + i)
            if (ms < 1) ms = mdead
            mw = map((j-1-1)*nx_block + i)
            if (mw < 1) mw = mdead
            if (dpx(mw)/=c0 .or. dpy(mw)/=c0   &
                            .or.               &
                dpx(ms)/=c0 .or. dpy(ms)/=c0) then
               icellsd = icellsd + 1
               indxid(icellsd) = i
               indxjd(icellsd) = j
            else
               if ( abs(edgearea(ms)) > c0 ) then ! 1 grid-cell wide channel: dpx,y = 0, edgearea /= 0
                  icellsd = icellsd + 1
                  indxid(icellsd) = i
                  indxjd(icellsd) = j
               endif
            endif
         enddo
         enddo
      endif       ! edge = north/east

      !-------------------------------------------------------------------
      ! Scale the departure points
      !-------------------------------------------------------------------

      ! ⚠ The BOUNDS matter, not just the dead slot.  The block form ran
      ! 1..ihi x 1..jhi -- not the whole block -- and dxu is zero in the
      ! outer ghost ring on a real grid, so iterating every slot divides
      ! 0/0 there and -fpe0 traps it.  A uniform box grid has dxu = 16e5
      ! everywhere including ghosts, which is why no local configuration
      ! shows this and decomp_suite does.
      !
      ! Cells outside these bounds keep c0: the block form left dx
      ! undefined there and never read it, so any value is safe.
      dx(:) = c0
      dy(:) = c0
      do j = 1, jhi
      do i = 1, ihi
         ms = map((j-1)*nx_block + i)
         if (ms > 0) then
            dx(ms) = dpx(ms) / dxu(ms)
            dy(ms) = dpy(ms) / dyu(ms)
         endif
      enddo
      enddo

      !-------------------------------------------------------------------
      ! Compute departure regions, divide into triangles, and locate
      !  vertices of each triangle.
      ! Work in a nondimensional coordinate system in which lengths are
      !  scaled by the local metric coefficients (dxu and dyu).
      ! Note: The do loop includes north faces of the j = 1 ghost cells
      !       when edge = 'north'.  The loop includes east faces of i = 1
      !       ghost cells when edge = 'east'.
      !-------------------------------------------------------------------

      do ij = 1, icellsd
         i = indxid(ij)
         j = indxjd(ij)

         xcl = -p5
         ycl =  c0

         xcr =  p5
         ycr =  c0

         ! Departure points

         ms = map((j-1)*nx_block + i)
         if (ms < 1) ms = mdead
         if (trim(edge) == 'north') then ! north edge
            mw = map((j-1)*nx_block + i-1)
            if (mw < 1) mw = mdead
            xdl = xcl + dx(mw)
            ydl = ycl + dy(mw)
            xdr = xcr + dx(ms)
            ydr = ycr + dy(ms)
         else                   ! east edge; rotate trajectory by pi/2
            mw = map((j-1-1)*nx_block + i)
            if (mw < 1) mw = mdead
            xdl = xcl - dy(ms)
            ydl = ycl + dx(ms)
            xdr = xcr - dy(mw)
            ydr = ycr + dx(mw)
         endif

         xdm = p5 * (xdr + xdl)
         ydm = p5 * (ydr + ydl)

         ! Intersection points

         xil = xcl
         yil = (xcl*(ydm-ydl) + xdm*ydl - xdl*ydm) / (xdm - xdl)

         xir = xcr
         yir = (xcr*(ydr-ydm) - xdm*ydr + xdr*ydm) / (xdr - xdm)

         md = (ydr - ydl) / (xdr - xdl)

         if (abs(md) > puny) then
            xic = xdl - ydl/md
         else
            xic = c0
         endif
         yic = c0

         xicl = xic
         yicl = yic
         xicr = xic
         yicr = yic

         !-------------------------------------------------------------------
         ! Locate triangles in TL cell (NW for north edge, NE for east edge)
         ! and BL cell (W for north edge, N for east edge).
         !
         ! areafact_c or areafac_ce (areafact_c for the other edge) are used
         ! (with shifted indices) to make sure that a flux area on one edge
         ! is consistent with the analogous area on the other edge and to
         ! ensure that areas add up when using l_edge_flux_adj = T.
         ! See PR #849 for details.
         !-------------------------------------------------------------------

         if (yil > c0 .and. xdl < xcl .and. ydl >= c0) then

            ! TL (group 1)

            ng = 1
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xil
            yp    (ms,2,ng) = yil
            xp    (ms,3,ng) = xdl
            yp    (ms,3,ng) = ydl
            iflux   (ms,ng) = i + ishift_tl
            jflux   (ms,ng) = j + jshift_tl
            areafact(ms,ng) = -areafac_ce(i+ise_tl,j+jse_tl)

         elseif (yil < c0 .and. xdl < xcl .and. ydl < c0) then

            ! BL (group 1)

            ng = 1
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xdl
            yp    (ms,2,ng) = ydl
            xp    (ms,3,ng) = xil
            yp    (ms,3,ng) = yil
            iflux   (ms,ng) = i + ishift_bl
            jflux   (ms,ng) = j + jshift_bl
            areafact(ms,ng) = areafac_ce(i+ise_bl,j+jse_bl)

         elseif (yil < c0 .and. xdl < xcl .and. ydl >= c0) then

            ! TL1 (group 1)

            ng = 1
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xdl
            yp    (ms,2,ng) = ydl
            xp    (ms,3,ng) = xic
            yp    (ms,3,ng) = yic
            iflux   (ms,ng) = i + ishift_tl
            jflux   (ms,ng) = j + jshift_tl
            areafact(ms,ng) = areafac_c(i+is_l,j+js_l)

            ! BL1 (group 3)

            ng = 3
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xic
            yp    (ms,2,ng) = yic
            xp    (ms,3,ng) = xil
            yp    (ms,3,ng) = yil
            iflux   (ms,ng) = i + ishift_bl
            jflux   (ms,ng) = j + jshift_bl
            areafact(ms,ng) = areafac_ce(i+ise_bl,j+jse_bl)

         elseif (yil > c0 .and. xdl < xcl .and. ydl < c0) then

            ! TL2 (group 3)

            ng = 3
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xil
            yp    (ms,2,ng) = yil
            xp    (ms,3,ng) = xic
            yp    (ms,3,ng) = yic
            iflux   (ms,ng) = i + ishift_tl
            jflux   (ms,ng) = j + jshift_tl
            areafact(ms,ng) = -areafac_ce(i+ise_tl,j+jse_tl)

            ! BL2 (group 1)

            ng = 1
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xic
            yp    (ms,2,ng) = yic
            xp    (ms,3,ng) = xdl
            yp    (ms,3,ng) = ydl
            iflux   (ms,ng) = i + ishift_bl
            jflux   (ms,ng) = j + jshift_bl
            areafact(ms,ng) = -areafac_c(i+is_l,j+js_l)

         endif                  ! TL and BL triangles

         !-------------------------------------------------------------------
         ! Locate triangles in TR cell (NE for north edge, SE for east edge)
         ! and in BR cell (E for north edge, S for east edge).
         !-------------------------------------------------------------------

         if (yir > c0 .and. xdr >= xcr .and. ydr >= c0) then

            ! TR (group 2)

            ng = 2
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xdr
            yp    (ms,2,ng) = ydr
            xp    (ms,3,ng) = xir
            yp    (ms,3,ng) = yir
            iflux   (ms,ng) = i + ishift_tr
            jflux   (ms,ng) = j + jshift_tr
            areafact(ms,ng) = -areafac_ce(i+ise_tr,j+jse_tr)

         elseif (yir < c0 .and. xdr >= xcr .and. ydr < c0) then

            ! BR (group 2)

            ng = 2
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xir
            yp    (ms,2,ng) = yir
            xp    (ms,3,ng) = xdr
            yp    (ms,3,ng) = ydr
            iflux   (ms,ng) = i + ishift_br
            jflux   (ms,ng) = j + jshift_br
            areafact(ms,ng) = areafac_ce(i+ise_br,j+jse_br)

         elseif (yir < c0 .and. xdr >= xcr  .and. ydr >= c0) then

            ! TR1 (group 2)

            ng = 2
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xic
            yp    (ms,2,ng) = yic
            xp    (ms,3,ng) = xdr
            yp    (ms,3,ng) = ydr
            iflux   (ms,ng) = i + ishift_tr
            jflux   (ms,ng) = j + jshift_tr
            areafact(ms,ng) = areafac_c(i+is_r,j+js_r)

            ! BR1 (group 3)

            ng = 3
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xir
            yp    (ms,2,ng) = yir
            xp    (ms,3,ng) = xic
            yp    (ms,3,ng) = yic
            iflux   (ms,ng) = i + ishift_br
            jflux   (ms,ng) = j + jshift_br
            areafact(ms,ng) = areafac_ce(i+ise_br,j+jse_br)

         elseif (yir > c0 .and. xdr >= xcr .and. ydr < c0) then

            ! TR2 (group 3)

            ng = 3
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xic
            yp    (ms,2,ng) = yic
            xp    (ms,3,ng) = xir
            yp    (ms,3,ng) = yir
            iflux   (ms,ng) = i + ishift_tr
            jflux   (ms,ng) = j + jshift_tr
            areafact(ms,ng) = -areafac_ce(i+ise_tr,j+jse_tr)

            ! BR2 (group 2)

            ng = 2
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xdr
            yp    (ms,2,ng) = ydr
            xp    (ms,3,ng) = xic
            yp    (ms,3,ng) = yic
            iflux   (ms,ng) = i + ishift_br
            jflux   (ms,ng) = j + jshift_br
            areafact(ms,ng) = -areafac_c(i+is_r,j+js_r)

         endif                  ! TR and BR triangles

         !-------------------------------------------------------------------
         ! Redefine departure points if not located in central cells (TC or BC)
         !-------------------------------------------------------------------

         if (xdl < xcl) then
            xdl = xil
            ydl = yil
         endif

         if (xdr > xcr) then
            xdr = xir
            ydr = yir
         endif

         !-------------------------------------------------------------------
         ! For l_edge_flux_adj = T, shift the midpoint so that the departure
         ! region has the prescribed area
         !-------------------------------------------------------------------

         if (l_edge_flux_adj) then

            ! Sum the areas of the left and right triangles.
            ! Note that yp(ms,1,ng) = 0 for all triangles, so we can
            !  drop those terms from the area formula.

            ng = 1
            area1 = p5 * ( (xp(ms,2,ng)-xp(ms,1,ng)) *   &
                            yp(ms,3,ng)                   &
                         -  yp(ms,2,ng) *                 &
                           (xp(ms,3,ng)-xp(ms,1,ng)) )   &
                         * areafact(ms,ng)

            ng = 2
            area2 = p5 * ( (xp(ms,2,ng)-xp(ms,1,ng)) *   &
                            yp(ms,3,ng)                   &
                         -  yp(ms,2,ng) *                 &
                           (xp(ms,3,ng)-xp(ms,1,ng)) )   &
                         * areafact(ms,ng)

            ng = 3
            area3 = p5 * ( (xp(ms,2,ng)-xp(ms,1,ng)) *   &
                            yp(ms,3,ng)                   &
                         -  yp(ms,2,ng) *                 &
                           (xp(ms,3,ng)-xp(ms,1,ng)) )   &
                         * areafact(ms,ng)

            !-----------------------------------------------------------
            ! Check whether the central triangles lie in one grid cell or two.
            ! If all are in one grid cell, then adjust the area of the central
            !  region so that the sum of all triangle areas is equal to the
            !  prescribed value.
            ! If two triangles are in one grid cell and one is in the other,
            !  then compute the area of the lone triangle. Then adjust
            !  the area of the remaining two-triangle region so that the sum of
            !  all triangle areas has the prescribed value.
            !-----------------------------------------------------------

            if (ydl*ydr >= c0) then   ! Both DPs lie on same side of x-axis

               ! compute required area of central departure region
               area_c  = edgearea(ms) - area1 - area2 - area3

               ! shift midpoint so that the area of remaining triangles = area_c
               w1 = c2*area_c/areafac_c(i,j)    &
                    + (xdr-xcl)*ydl + (xcr-xdl)*ydr
               w2 = (xdr-xdl)**2 + (ydr-ydl)**2
               w1 = w1/w2
               xdm = xdm + (ydr - ydl) * w1
               ydm = ydm - (xdr - xdl) * w1

               ! compute left and right intersection points
               mdl = (ydm - ydl) / (xdm - xdl)
               mdr = (ydr - ydm) / (xdr - xdm)

               if (abs(mdl) > puny) then
                  xicl = xdl - ydl/mdl
               else
                  xicl = c0
               endif
               yicl = c0

               if (abs(mdr) > puny) then
                  xicr = xdr - ydr/mdr
               else
                  xicr = c0
               endif
               yicr = c0

            elseif (xic < c0 .and. xic > xcl) then  ! fix ICL = IC

               xicl = xic
               yicl = yic

               ! compute midpoint between ICL and DR
               xdm = p5 * (xdr + xicl)
               ydm = p5 *  ydr

               ! compute area of (lone) triangle adjacent to left corner
               area4 = p5 * (xcl - xic) * ydl * areafac_c(i,j)
               area_c  = edgearea(ms) - area1 - area2 - area3 - area4

               ! shift midpoint so that area of remaining triangles = area_c
               w1 = c2*area_c/areafac_c(i,j) + (xcr-xic)*ydr
               w2 = (xdr-xic)**2 + ydr**2
               w1 = w1/w2
               xdm = xdm + ydr*w1
               ydm = ydm - (xdr - xic) * w1

               ! compute ICR
               mdr = (ydr - ydm) / (xdr - xdm)
               if (abs(mdr) > puny) then
                  xicr = xdr - ydr/mdr
               else
                  xicr = c0
               endif
               yicr = c0

            elseif (xic >= c0 .and. xic < xcr) then  ! fix ICR = IR

               xicr = xic
               yicr = yic

               ! compute midpoint between ICR and DL
               xdm = p5 * (xicr + xdl)
               ydm = p5 *  ydl

               ! compute area of (lone) triangle adjacent to right corner
               area4 = p5 * (xic - xcr) * ydr * areafac_c(i,j)
               area_c  = edgearea(ms) - area1 - area2 - area3 - area4

               ! shift midpoint so that area of remaining triangles = area_c
               w1 = c2*area_c/areafac_c(i,j) + (xic-xcl)*ydl
               w2 = (xic-xdl)**2 + ydl**2
               w1 = w1/w2
               xdm = xdm - ydl*w1
               ydm = ydm - (xic - xdl) * w1

               ! compute ICL

               mdl = (ydm - ydl) / (xdm - xdl)
               if (abs(mdl) > puny) then
                  xicl = xdl - ydl/mdl
               else
                  xicl = c0
               endif
               yicl = c0

            endif   ! ydl*ydr >= c0

         endif  ! l_edge_flux_adj

         !-------------------------------------------------------------------
         ! Locate triangles in BC cell (H for both north and east edges)
         ! and TC cell (N for north edge and E for east edge).
         !-------------------------------------------------------------------

         ! Start with cases where both DPs lie in the same grid cell

         if (ydl >= c0 .and. ydr >= c0 .and. ydm >= c0) then

            ! TC1a (group 4)

            ng = 4
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xcr
            yp    (ms,2,ng) = ycr
            xp    (ms,3,ng) = xdl
            yp    (ms,3,ng) = ydl
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

            ! TC2a (group 5)

            ng = 5
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xdr
            yp    (ms,2,ng) = ydr
            xp    (ms,3,ng) = xdl
            yp    (ms,3,ng) = ydl
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

            ! TC3a (group 6)

            ng = 6
            xp    (ms,1,ng) = xdl
            yp    (ms,1,ng) = ydl
            xp    (ms,2,ng) = xdr
            yp    (ms,2,ng) = ydr
            xp    (ms,3,ng) = xdm
            yp    (ms,3,ng) = ydm
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

         elseif (ydl >= c0 .and. ydr >= c0 .and. ydm < c0) then  ! rare

            ! TC1b (group 4)

            ng = 4
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xicl
            yp    (ms,2,ng) = yicl
            xp    (ms,3,ng) = xdl
            yp    (ms,3,ng) = ydl
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

            ! TC2b (group 5)

            ng = 5
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xdr
            yp    (ms,2,ng) = ydr
            xp    (ms,3,ng) = xicr
            yp    (ms,3,ng) = yicr
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

            ! BC3b (group 6)

            ng = 6
            xp    (ms,1,ng) = xicr
            yp    (ms,1,ng) = yicr
            xp    (ms,2,ng) = xicl
            yp    (ms,2,ng) = yicl
            xp    (ms,3,ng) = xdm
            yp    (ms,3,ng) = ydm
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

         elseif (ydl <= c0 .and. ydr <= c0 .and. ydm <= c0) then

            ! BC1a (group 4)

            ng = 4
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xdl
            yp    (ms,2,ng) = ydl
            xp    (ms,3,ng) = xcr
            yp    (ms,3,ng) = ycr
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

            ! BC2a (group 5)

            ng = 5
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xdl
            yp    (ms,2,ng) = ydl
            xp    (ms,3,ng) = xdr
            yp    (ms,3,ng) = ydr
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

            ! BC3a (group 6)

            ng = 6
            xp    (ms,1,ng) = xdl
            yp    (ms,1,ng) = ydl
            xp    (ms,2,ng) = xdm
            yp    (ms,2,ng) = ydm
            xp    (ms,3,ng) = xdr
            yp    (ms,3,ng) = ydr
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

         elseif (ydl <= c0 .and. ydr <= c0 .and. ydm > c0) then  ! rare

            ! BC1b (group 4)

            ng = 4
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xdl
            yp    (ms,2,ng) = ydl
            xp    (ms,3,ng) = xicl
            yp    (ms,3,ng) = yicl
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

            ! BC2b (group 5)

            ng = 5
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xicr
            yp    (ms,2,ng) = yicr
            xp    (ms,3,ng) = xdr
            yp    (ms,3,ng) = ydr
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

            ! TC3b (group 6)

            ng = 6
            xp    (ms,1,ng) = xicl
            yp    (ms,1,ng) = yicl
            xp    (ms,2,ng) = xicr
            yp    (ms,2,ng) = yicr
            xp    (ms,3,ng) = xdm
            yp    (ms,3,ng) = ydm
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

         ! Now consider cases where the two DPs lie in different grid cells

         elseif (ydl > c0 .and. ydr < c0 .and. xic >= c0  &
                                         .and. ydm >= c0) then

            ! TC1b (group 4)

            ng = 4
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xicr
            yp    (ms,2,ng) = yicr
            xp    (ms,3,ng) = xdl
            yp    (ms,3,ng) = ydl
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

            ! BC2b (group 5) lone triangle

            ng = 5
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xicr
            yp    (ms,2,ng) = yicr
            xp    (ms,3,ng) = xdr
            yp    (ms,3,ng) = ydr
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

            ! TC3b (group 6)

            ng = 6
            xp    (ms,1,ng) = xdl
            yp    (ms,1,ng) = ydl
            xp    (ms,2,ng) = xicr
            yp    (ms,2,ng) = yicr
            xp    (ms,3,ng) = xdm
            yp    (ms,3,ng) = ydm
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

         elseif (ydl > c0 .and. ydr < c0 .and. xic >= c0  &
                                         .and. ydm < c0 ) then  ! less common

            ! TC1b (group 4)

            ng = 4
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xicl
            yp    (ms,2,ng) = yicl
            xp    (ms,3,ng) = xdl
            yp    (ms,3,ng) = ydl
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

            ! BC2b (group 5) lone triangle

            ng = 5
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xicr
            yp    (ms,2,ng) = yicr
            xp    (ms,3,ng) = xdr
            yp    (ms,3,ng) = ydr
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

            ! BC3b (group 6)

            ng = 6
            xp    (ms,1,ng) = xicr
            yp    (ms,1,ng) = yicr
            xp    (ms,2,ng) = xicl
            yp    (ms,2,ng) = yicl
            xp    (ms,3,ng) = xdm
            yp    (ms,3,ng) = ydm
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

         elseif (ydl > c0 .and. ydr < c0 .and. xic < c0   &
                                         .and. ydm < c0) then

            ! TC1b (group 4) lone triangle

            ng = 4
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xicl
            yp    (ms,2,ng) = yicl
            xp    (ms,3,ng) = xdl
            yp    (ms,3,ng) = ydl
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

            ! BC2b (group 5)

            ng = 5
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xicl
            yp    (ms,2,ng) = yicl
            xp    (ms,3,ng) = xdr
            yp    (ms,3,ng) = ydr
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

            ! BC3b (group 6)

            ng = 6
            xp    (ms,1,ng) = xdr
            yp    (ms,1,ng) = ydr
            xp    (ms,2,ng) = xicl
            yp    (ms,2,ng) = yicl
            xp    (ms,3,ng) = xdm
            yp    (ms,3,ng) = ydm
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

         elseif (ydl > c0 .and. ydr < c0 .and. xic <  c0  &
                                         .and. ydm >= c0) then  ! less common

            ! TC1b (group 4) lone triangle

            ng = 4
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xicl
            yp    (ms,2,ng) = yicl
            xp    (ms,3,ng) = xdl
            yp    (ms,3,ng) = ydl
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

            ! BC2b (group 5)

            ng = 5
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xicr
            yp    (ms,2,ng) = yicr
            xp    (ms,3,ng) = xdr
            yp    (ms,3,ng) = ydr
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

            ! TC3b (group 6)

            ng = 6
            xp    (ms,1,ng) = xicl
            yp    (ms,1,ng) = yicl
            xp    (ms,2,ng) = xicr
            yp    (ms,2,ng) = yicr
            xp    (ms,3,ng) = xdm
            yp    (ms,3,ng) = ydm
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

         elseif (ydl < c0 .and. ydr > c0 .and. xic <  c0  &
                                         .and. ydm >= c0) then

            ! BC1b (group 4) lone triangle

            ng = 4
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xdl
            yp    (ms,2,ng) = ydl
            xp    (ms,3,ng) = xicl
            yp    (ms,3,ng) = yicl
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

            ! TC2b (group 5)

            ng = 5
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xdr
            yp    (ms,2,ng) = ydr
            xp    (ms,3,ng) = xicl
            yp    (ms,3,ng) = yicl
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

            ! TC3b (group 6)

            ng = 6
            xp    (ms,1,ng) = xicl
            yp    (ms,1,ng) = yicl
            xp    (ms,2,ng) = xdr
            yp    (ms,2,ng) = ydr
            xp    (ms,3,ng) = xdm
            yp    (ms,3,ng) = ydm
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

         elseif (ydl < c0 .and. ydr > c0 .and. xic < c0  &
                                         .and. ydm < c0) then ! less common

            ! BC1b (group 4) lone triangle

            ng = 4
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xdl
            yp    (ms,2,ng) = ydl
            xp    (ms,3,ng) = xicl
            yp    (ms,3,ng) = yicl
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

            ! TC2b (group 5)

            ng = 5
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xdr
            yp    (ms,2,ng) = ydr
            xp    (ms,3,ng) = xicr
            yp    (ms,3,ng) = yicr
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

            ! BC3b (group 6)

            ng = 6
            xp    (ms,1,ng) = xicr
            yp    (ms,1,ng) = yicr
            xp    (ms,2,ng) = xicl
            yp    (ms,2,ng) = yicl
            xp    (ms,3,ng) = xdm
            yp    (ms,3,ng) = ydm
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

         elseif (ydl < c0 .and. ydr > c0 .and. xic >= c0  &
                                         .and. ydm <  c0) then

            ! BC1b (group 4)

            ng = 4
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xdl
            yp    (ms,2,ng) = ydl
            xp    (ms,3,ng) = xicr
            yp    (ms,3,ng) = yicr
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

            ! TC2b (group 5) lone triangle

            ng = 5
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xdr
            yp    (ms,2,ng) = ydr
            xp    (ms,3,ng) = xicr
            yp    (ms,3,ng) = yicr
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

            ! BC3b (group 6)

            ng = 6
            xp    (ms,1,ng) = xicr
            yp    (ms,1,ng) = yicr
            xp    (ms,2,ng) = xdl
            yp    (ms,2,ng) = ydl
            xp    (ms,3,ng) = xdm
            yp    (ms,3,ng) = ydm
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

         elseif (ydl < c0 .and. ydr > c0 .and. xic >= c0   &
                                         .and. ydm >= c0) then  ! less common

            ! BC1b (group 4)

            ng = 4
            xp    (ms,1,ng) = xcl
            yp    (ms,1,ng) = ycl
            xp    (ms,2,ng) = xdl
            yp    (ms,2,ng) = ydl
            xp    (ms,3,ng) = xicl
            yp    (ms,3,ng) = yicl
            iflux   (ms,ng) = i + ishift_bc
            jflux   (ms,ng) = j + jshift_bc
            areafact(ms,ng) = areafac_c(i,j)

            ! TC2b (group 5) lone triangle

            ng = 5
            xp    (ms,1,ng) = xcr
            yp    (ms,1,ng) = ycr
            xp    (ms,2,ng) = xdr
            yp    (ms,2,ng) = ydr
            xp    (ms,3,ng) = xicr
            yp    (ms,3,ng) = yicr
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

            ! TC3b (group 6)

            ng = 6
            xp    (ms,1,ng) = xicl
            yp    (ms,1,ng) = yicl
            xp    (ms,2,ng) = xicr
            yp    (ms,2,ng) = yicr
            xp    (ms,3,ng) = xdm
            yp    (ms,3,ng) = ydm
            iflux   (ms,ng) = i + ishift_tc
            jflux   (ms,ng) = j + jshift_tc
            areafact(ms,ng) = -areafac_c(i,j)

         endif                  ! TC and BC triangles

      enddo                     ! ij

      !-------------------------------------------------------------------
      ! Compute triangle areas with appropriate sign.
      ! These are found by computing the area in scaled coordinates and
      !  multiplying by a scale factor (areafact).
      ! Note that the scale factor is positive for fluxes out of the cell
      !  and negative for fluxes into the cell.
      !
      ! Note: The triangle area formula below gives A >=0 iff the triangle
      !        points x1, x2, and x3 are taken in counterclockwise order.
      !       These points are defined above in such a way that the
      !        order is nearly always CCW.
      !       In rare cases, we may compute A < 0.  In this case,
      !        the quadrilateral departure area is equal to the
      !        difference of two triangle areas instead of the sum.
      !        The fluxes work out correctly in the end.
      !
      ! Also compute the cumulative area transported across each edge.
      ! If l_edge_flux_adj = T, this area is compared to edgearea as a bug check.
      ! If l_edge_flux_adj = F, this area is passed as an output array.
      !-------------------------------------------------------------------

      areasum(:) = c0

      do ng = 1, ngroups
         icells(ng) = 0

         do ij = 1, icellsd
            i = indxid(ij)
            j = indxjd(ij)
            ms = map((j-1)*nx_block + i)
            ! icellsd cells always have a slot: dpx is nonzero only where
            ! departure_points ran, which is the slotted range.
            if (ms < 1) cycle

            triarea(ms,ng) = p5 * ( (xp(ms,2,ng)-xp(ms,1,ng)) *   &
                                     (yp(ms,3,ng)-yp(ms,1,ng))   &
                                   - (yp(ms,2,ng)-yp(ms,1,ng)) *   &
                                     (xp(ms,3,ng)-xp(ms,1,ng)) )   &
                                   * areafact(ms,ng)

            if (abs(triarea(ms,ng)) < eps16*areafac_c(i,j)) then
               triarea(ms,ng) = c0
            else
               icells(ng) = icells(ng) + 1
               ic = icells(ng)
               indxs(ic,ng) = ms
            endif

            areasum(ms) = areasum(ms) + triarea(ms,ng)

         enddo                  ! ij
      enddo                     ! ng

      if (l_edge_flux_adj) then
         if (bugcheck) then   ! set bugcheck = F to speed up code
            do ij = 1, icellsd
               i = indxid(ij)
               j = indxjd(ij)
               ms = map((j-1)*nx_block + i)
               ! icellsd cells always have a slot: dpx is nonzero only where
               ! departure_points ran, which is the slotted range.
               if (ms < 1) cycle
               if ( abs(areasum(ms) - edgearea(ms)) > eps13*areafac_c(i,j) .and. abs(edgearea(ms)) > c0 ) then
                  write(nu_diag,*) ''
                  write(nu_diag,*) 'Areas do not add up: m, i, j, edge =',   &
                           my_task, i, j, trim(edge)
                  write(nu_diag,*) 'edgearea =', edgearea(ms)
                  write(nu_diag,*) 'areasum =', areasum(ms)
                  write(nu_diag,*) 'areafac_c =', areafac_c(i,j)
                  write(nu_diag,*) ''
                  write(nu_diag,*) 'Triangle areas:'
                  do ng = 1, ngroups   ! not vector friendly
                     if (abs(triarea(ms,ng)) > eps16*abs(areafact(ms,ng))) then
                        write(nu_diag,*) ng, triarea(ms,ng)
                     endif
                 enddo
               endif
            enddo
         endif          ! bugcheck

      else            ! l_edge_flux_adj = F
         do ij = 1, icellsd
            i = indxid(ij)
            j = indxjd(ij)
            edgearea(ms) = areasum(ms)
         enddo
      endif     ! l_edge_flux_adj

      ! The donor cell as a SLOT, derived once from the (i,j) form rather
      ! than at each of the 51 places iflux/jflux are assigned.  iflux and
      ! jflux stay because the vertex transform below needs the integer
      ! shift, which a slot cannot give.
      !
      ! A donor with no slot leaves 0 here; transport_integrals resolves
      ! that to the dead slot, where mc and tc are zero -- which is what the
      ! block layout carries at a land cell.
      do ng = 1, ngroups
         do ij = 1, icells(ng)
            ms = indxs(ij,ng)
            isflux(ms,ng) = map((jflux(ms,ng)-1)*nx_block + iflux(ms,ng))
         enddo
      enddo

      !-------------------------------------------------------------------
      ! Transform triangle vertices to a scaled coordinate system centered
      !  in the cell containing the triangle.
      !-------------------------------------------------------------------

      if (trim(edge) == 'north') then
         do ng = 1, ngroups
            do nv = 1, nvert
               do ij = 1, icells(ng)
                  ms = indxs(ij,ng)
                  ijb = ijslot(ms)
                  i = mod(ijb-1, nx_block) + 1
                  j = (ijb-1)/nx_block + 1
                  ishift = iflux(ms,ng) - i
                  jshift = jflux(ms,ng) - j
                  xp(ms,nv,ng) = xp(ms,nv,ng) - c1*ishift
                  yp(ms,nv,ng) = yp(ms,nv,ng) + p5 - c1*jshift
               enddo            ! ij
            enddo               ! nv
         enddo                  ! ng
      else                      ! east edge
         do ng = 1, ngroups
            do nv = 1, nvert
               do ij = 1, icells(ng)
                  ms = indxs(ij,ng)
                  ijb = ijslot(ms)
                  i = mod(ijb-1, nx_block) + 1
                  j = (ijb-1)/nx_block + 1
                  ishift = iflux(ms,ng) - i
                  jshift = jflux(ms,ng) - j
                  ! Note rotation of pi/2 here
                  w1 = xp(ms,nv,ng)
                  xp(ms,nv,ng) =  yp(ms,nv,ng) + p5 - c1*ishift
                  yp(ms,nv,ng) = -w1 - c1*jshift
               enddo            ! ij
            enddo               ! nv
         enddo                  ! ng
      endif

      if (bugcheck) then
         if (trim(edge) == 'north') then
            ib = ilo
            jb = jlo-1
         else ! east edge
            ib = ilo-1
            jb = jlo
         endif
         do ng = 1, ngroups
         do nv = 1, nvert
            do j = jb, jhi
            do i = ib, ihi
               if (abs(triarea(ms,ng)) > puny) then
                  if (abs(xp(ms,nv,ng)) > p5+puny) then
                     write(nu_diag,*) ''
                     write(nu_diag,*) 'WARNING: xp =', xp(ms,nv,ng)
                     write(nu_diag,*) 'm, i, j, ng, nv =', my_task, i, j, ng, nv
                  endif
                  if (abs(yp(ms,nv,ng)) > p5+puny) then
                     write(nu_diag,*) ''
                     write(nu_diag,*) 'WARNING: yp =', yp(ms,nv,ng)
                     write(nu_diag,*) 'm, i, j, ng, nv =', my_task, i, j, ng, nv
                  endif
               endif   ! triarea
            enddo
            enddo
         enddo
         enddo
      endif  ! bugcheck

      end subroutine locate_triangles

!=======================================================================
!
! For each triangle, find the coordinates of the quadrature points needed
!  to compute integrals of linear, quadratic, or cubic polynomials,
!  using formulas from A.H. Stroud, Approximate Calculation of Multiple
!  Integrals, Prentice-Hall, 1971.  (Section 8.8, formula 3.1.)
! Linear functions can be integrated exactly by evaluating the function
!  at just one point (the midpoint).  Quadratic functions require
!  3 points, and cubics require 4 points.
! The default is cubic, but the code can be sped up slightly using
!  linear or quadratic integrals, usually with little loss of accuracy.
!
! The formulas are as follows:
!
! I1 = integral of f(x,y)*dA
!    = A * f(x0,y0)
! where A is the traingle area and (x0,y0) is the midpoint.
!
! I2 = A * (f(x1,y1) + f(x2,y2) + f(x3,y3))
! where these three points are located halfway between the midpoint
! and the three vertics of the triangle.
!
! I3 = A * [ -9/16 *  f(x0,y0)
!           + 25/48 * (f(x1,y1) + f(x2,y2) + f(x3,y3))]
! where (x0,y0) is the midpoint, and the other three points are
! located 2/5 of the way from the midpoint to the three vertices.
!
! author William H. Lipscomb, LANL

      subroutine triangle_coordinates (integral_order, icells,    &
                                       indxs,                     &
                                       xp,             yp)

      ! Entirely pointwise on xp and yp, so the slot list is all it needs --
      ! no i, no j, and no block dimensions.

      integer (kind=int_kind), intent(in) ::   &
         integral_order      ! polynomial order for quadrature integrals

      integer (kind=int_kind), dimension (ngroups), intent(in) ::     &
         icells              ! number of cells where triarea > puny

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         indxs   ! slot of each cell in the list

      real (kind=dbl_kind), intent(inout), dimension (:,0:,:) :: &
         xp, yp          ! coordinates of triangle points

      ! local variables

      integer (kind=int_kind) ::   &
         ij                , & ! list index
         ms                , & ! slot that list entry owns
         ng                  ! triangle index

      character(len=*), parameter :: subname = '(triangle_coordinates)'

      if (integral_order == 1) then ! linear (1-point formula)

         do ng = 1, ngroups
         do ij = 1, icells(ng)
            ms = indxs(ij,ng)

            ! coordinates of midpoint
            xp(ms,0,ng) = p333   &
                        * (xp(ms,1,ng) + xp(ms,2,ng) + xp(ms,3,ng))
            yp(ms,0,ng) = p333   &
                        * (yp(ms,1,ng) + yp(ms,2,ng) + yp(ms,3,ng))

         enddo                  ! ij
         enddo                  ! ng

      elseif (integral_order == 2) then ! quadratic (3-point formula)

         do ng = 1, ngroups
         do ij = 1, icells(ng)
            ms = indxs(ij,ng)

            ! coordinates of midpoint
            xp(ms,0,ng) = p333   &
                        * (xp(ms,1,ng) + xp(ms,2,ng) + xp(ms,3,ng))
            yp(ms,0,ng) = p333   &
                        * (yp(ms,1,ng) + yp(ms,2,ng) + yp(ms,3,ng))

            ! coordinates of the 3 points needed for integrals

            xp(ms,1,ng) = p5*xp(ms,1,ng) + p5*xp(ms,0,ng)
            yp(ms,1,ng) = p5*yp(ms,1,ng) + p5*yp(ms,0,ng)

            xp(ms,2,ng) = p5*xp(ms,2,ng) + p5*xp(ms,0,ng)
            yp(ms,2,ng) = p5*yp(ms,2,ng) + p5*yp(ms,0,ng)

            xp(ms,3,ng) = p5*xp(ms,3,ng) + p5*xp(ms,0,ng)
            yp(ms,3,ng) = p5*yp(ms,3,ng) + p5*yp(ms,0,ng)

         enddo                  ! ij
         enddo                  ! ng

      else                      ! cubic (4-point formula)

         do ng = 1, ngroups
         do ij = 1, icells(ng)
            ms = indxs(ij,ng)

            ! coordinates of midpoint
            xp(ms,0,ng) = p333   &
                        * (xp(ms,1,ng) + xp(ms,2,ng) + xp(ms,3,ng))
            yp(ms,0,ng) = p333   &
                        * (yp(ms,1,ng) + yp(ms,2,ng) + yp(ms,3,ng))

            ! coordinates of the other 3 points needed for integrals

            xp(ms,1,ng) = p4*xp(ms,1,ng) + p6*xp(ms,0,ng)
            yp(ms,1,ng) = p4*yp(ms,1,ng) + p6*yp(ms,0,ng)

            xp(ms,2,ng) = p4*xp(ms,2,ng) + p6*xp(ms,0,ng)
            yp(ms,2,ng) = p4*yp(ms,2,ng) + p6*yp(ms,0,ng)

            xp(ms,3,ng) = p4*xp(ms,3,ng) + p6*xp(ms,0,ng)
            yp(ms,3,ng) = p4*yp(ms,3,ng) + p6*yp(ms,0,ng)

         enddo                  ! ij
         enddo                  ! ng

      endif

      end subroutine triangle_coordinates

!=======================================================================
!
! Compute the transports across each edge by integrating the mass
! and tracers over each departure triangle.
! Input variables have the same meanings as in the main subroutine.
! Repeated use of certain sums makes the calculation more efficient.
! Integral formulas are described in triangle_coordinates subroutine.
!
! author William H. Lipscomb, LANL

      subroutine transport_integrals (ntrace,         icells,      &
                                      indxs,                       &
                                      tracer_type,    depend,      &
                                      integral_order, triarea,     &
                                      isflux,                      &
                                      xp,             yp,          &
                                      mc,             mx,          &
                                      my,             mflx,        &
                                      tc,             tx,          &
                                      ty,             mtflx)

      integer (kind=int_kind), intent(in) ::   &
         ntrace              , & ! number of tracers in use
         integral_order   ! polynomial order for quadrature integrals

      integer (kind=int_kind), dimension (ntrace), intent(in) ::     &
         tracer_type       , & ! = 1, 2, or 3 (see comments above)
         depend              ! tracer dependencies (see above)

      integer (kind=int_kind), dimension (ngroups), intent(in) ::     &
         icells           ! number of cells where triarea > puny

      ! One list now.  indxi/indxj are gone: they existed to index the sum
      ! accumulators and the block-shaped flux outputs, and both of those
      ! are slot-indexed, so nothing in here names a cell as (i,j) any more.
      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         indxs   ! slot of each cell in the list

      real (kind=dbl_kind), intent(in), dimension (:,0:,:) :: &
         xp, yp           ! coordinates of triangle points

      real (kind=dbl_kind), intent(in), dimension (:,:) :: &
         triarea          ! triangle area

      ! The donor cell as a slot.  0 where that cell has none, which the
      ! loops below resolve to the dead slot.
      integer (kind=int_kind), intent(in), dimension (:,:) :: &
         isflux

      ! Read at the DONOR slot, so compressed.  mflx below stays block:
      ! update_fields reads it at neighbours and has not been converted.
      real (kind=dbl_kind), intent(in), dimension (:), contiguous :: &
         mc, mx, my

      ! Compressed: update_fields reads these at the W and S neighbours,
      ! so they have to be slot arrays rather than reached through dynIJ.
      real (kind=dbl_kind), intent(out), dimension (:), contiguous :: &
         mflx

      real (kind=dbl_kind), intent(in), dimension (:,:), optional, contiguous :: &
         tc, tx, ty

      real (kind=dbl_kind), intent(out), dimension (:,:), optional, contiguous :: &
         mtflx

      ! local variables

      integer (kind=int_kind) ::   &
         ij, ms     , & ! list index and slot of edge
         is2        , & ! slot of the cell contributing transport
         ng            , & ! triangle index
         nt, nt1         ! tracer indices

      real (kind=dbl_kind) ::   &
         m0, m1, m2, m3         , & ! mass field at internal points
         w0, w1, w2, w3           ! work variables

      ! Slot-indexed: each is written and read at the same edge cell, which
      ! already carries its slot, so they never needed the block rectangle.
      real (kind=dbl_kind), dimension (size(mc)) ::   &
         msum, mxsum, mysum     , & ! sum of mass, mass*x, and mass*y
         mxxsum, mxysum, myysum   ! sum of mass*x*x, mass*x*y, mass*y*y

      real (kind=dbl_kind), dimension (size(mc), ntrace) ::   &
         mtsum            , & ! sum of mass*tracer
         mtxsum           , & ! sum of mass*tracer*x
         mtysum             ! sum of mass*tracer*y

      character(len=*), parameter :: subname = '(transport_integrals)'

      !-------------------------------------------------------------------
      ! Initialize
      !-------------------------------------------------------------------

      mflx(:) = c0
      if (present(mtflx)) then
         do nt = 1, ntrace
            mtflx(:,nt) = c0
         enddo
      endif

      !-------------------------------------------------------------------
      ! Main loop
      !-------------------------------------------------------------------

      do ng = 1, ngroups

         if (integral_order == 1) then  ! linear (1-point formula)

            do ij = 1, icells(ng)
               ms = indxs(ij,ng)
               ms = indxs(ij,ng)

               is2 = isflux(ms,ng)

               if (is2 < 1) is2 = size(mc)

               ! mass transports

               m0 = mc(is2) + xp(ms,0,ng)*mx(is2)   &
                              + yp(ms,0,ng)*my(is2)
               msum(ms) = m0

               mflx(ms) = mflx(ms) + triarea(ms,ng)*msum(ms)

               ! quantities needed for tracer transports
               mxsum(ms)  =         m0*xp(ms,0,ng)
               mxxsum(ms) = mxsum(ms)*xp(ms,0,ng)
               mxysum(ms) = mxsum(ms)*yp(ms,0,ng)
               mysum(ms)  =         m0*yp(ms,0,ng)
               myysum(ms) = mysum(ms)*yp(ms,0,ng)
            enddo               ! ij

         elseif (integral_order == 2) then  ! quadratic (3-point formula)

            do ij = 1, icells(ng)
               ms = indxs(ij,ng)
               ms = indxs(ij,ng)

               is2 = isflux(ms,ng)

               if (is2 < 1) is2 = size(mc)

               ! mass transports
               ! Weighting factor of 1/3 is incorporated into the ice
               ! area terms m1, m2, and m3.
               m1 = p333 * (mc(is2) + xp(ms,1,ng)*mx(is2)   &
                                      + yp(ms,1,ng)*my(is2))
               m2 = p333 * (mc(is2) + xp(ms,2,ng)*mx(is2)   &
                                      + yp(ms,2,ng)*my(is2))
               m3 = p333 * (mc(is2) + xp(ms,3,ng)*mx(is2)   &
                                      + yp(ms,3,ng)*my(is2))
               msum(ms) = m1 + m2 + m3
               mflx(ms) = mflx(ms) + triarea(ms,ng)*msum(ms)

               ! quantities needed for mass_tracer transports
               w1 = m1 * xp(ms,1,ng)
               w2 = m2 * xp(ms,2,ng)
               w3 = m3 * xp(ms,3,ng)

               mxsum(ms) = w1 + w2 + w3

               mxxsum(ms) = w1*xp(ms,1,ng) + w2*xp(ms,2,ng)   &
                           + w3*xp(ms,3,ng)

               mxysum(ms) = w1*yp(ms,1,ng) + w2*yp(ms,2,ng)   &
                           + w3*yp(ms,3,ng)

               w1 = m1 * yp(ms,1,ng)
               w2 = m2 * yp(ms,2,ng)
               w3 = m3 * yp(ms,3,ng)

               mysum(ms) = w1 + w2 + w3

               myysum(ms) = w1*yp(ms,1,ng) + w2*yp(ms,2,ng)   &
                           + w3*yp(ms,3,ng)
            enddo               ! ij

         else                   ! cubic (4-point formula)

            do ij = 1, icells(ng)
               ms = indxs(ij,ng)
               ms = indxs(ij,ng)

               is2 = isflux(ms,ng)

               if (is2 < 1) is2 = size(mc)

               ! mass transports

               ! Weighting factors are incorporated into the
               ! terms m0, m1, m2, and m3.
               m0 = p5625m * (mc(is2) + xp(ms,0,ng)*mx(is2)   &
                                        + yp(ms,0,ng)*my(is2))
               m1 = p52083 * (mc(is2) + xp(ms,1,ng)*mx(is2)   &
                                        + yp(ms,1,ng)*my(is2))
               m2 = p52083 * (mc(is2) + xp(ms,2,ng)*mx(is2)   &
                                        + yp(ms,2,ng)*my(is2))
               m3 = p52083 * (mc(is2) + xp(ms,3,ng)*mx(is2)   &
                                        + yp(ms,3,ng)*my(is2))
               msum(ms) = m0 + m1 + m2 + m3
               mflx(ms) = mflx(ms) + triarea(ms,ng)*msum(ms)

               ! quantities needed for tracer transports
               w0 = m0 * xp(ms,0,ng)
               w1 = m1 * xp(ms,1,ng)
               w2 = m2 * xp(ms,2,ng)
               w3 = m3 * xp(ms,3,ng)

               mxsum(ms) = w0 + w1 + w2 + w3

               mxxsum(ms) = w0*xp(ms,0,ng) + w1*xp(ms,1,ng)   &
                           + w2*xp(ms,2,ng) + w3*xp(ms,3,ng)

               mxysum(ms) = w0*yp(ms,0,ng) + w1*yp(ms,1,ng)   &
                           + w2*yp(ms,2,ng) + w3*yp(ms,3,ng)

               w0 = m0 * yp(ms,0,ng)
               w1 = m1 * yp(ms,1,ng)
               w2 = m2 * yp(ms,2,ng)
               w3 = m3 * yp(ms,3,ng)

               mysum(ms) = w0 + w1 + w2 + w3

               myysum(ms) = w0*yp(ms,0,ng) + w1*yp(ms,1,ng)   &
                           + w2*yp(ms,2,ng) + w3*yp(ms,3,ng)

            enddo               ! ij

         endif                  ! integral_order

         ! mass * tracer transports

         if (present(mtflx)) then

            do nt = 1, ntrace
               if (tracer_type(nt)==1) then ! does not depend on another tracer

                  do ij = 1, icells(ng)
                     ms = indxs(ij,ng)
                     ms = indxs(ij,ng)

                     is2 = isflux(ms,ng)

                     if (is2 < 1) is2 = size(mc)

                     mtsum(ms,nt) =  msum(ms) * tc(is2,nt)   &
                                   + mxsum(ms) * tx(is2,nt)   &
                                   + mysum(ms) * ty(is2,nt)

                     mtflx(ms,nt) = mtflx(ms,nt)   &
                                 + triarea(ms,ng) * mtsum(ms,nt)

                     ! quantities needed for dependent tracers

                     mtxsum(ms,nt) =  mxsum(ms) * tc(is2,nt)   &
                                    + mxxsum(ms) * tx(is2,nt)   &
                                    + mxysum(ms) * ty(is2,nt)

                     mtysum(ms,nt) =  mysum(ms) * tc(is2,nt)   &
                                    + mxysum(ms) * tx(is2,nt)   &
                                    + myysum(ms) * ty(is2,nt)
                  enddo         ! ij

               elseif (tracer_type(nt)==2) then ! depends on another tracer
                  nt1 = depend(nt)

                  do ij = 1, icells(ng)
                     ms = indxs(ij,ng)
                     ms = indxs(ij,ng)

                     is2 = isflux(ms,ng)

                     if (is2 < 1) is2 = size(mc)

                     mtsum(ms,nt) =  mtsum(ms,nt1) * tc(is2,nt)   &
                                   + mtxsum(ms,nt1) * tx(is2,nt)   &
                                   + mtysum(ms,nt1) * ty(is2,nt)

                     mtflx(ms,nt) = mtflx(ms,nt)   &
                                   + triarea(ms,ng) * mtsum(ms,nt)
                  enddo         ! ij


               elseif (tracer_type(nt)==3) then ! depends on two tracers
                  nt1 = depend(nt)

                  do ij = 1, icells(ng)
                     ms = indxs(ij,ng)
                     ms = indxs(ij,ng)

                     is2 = isflux(ms,ng)

                     if (is2 < 1) is2 = size(mc)

                     ! upwind approx (tx=ty=0) for type 3 tracers
                     mtsum(ms,nt) =  mtsum(ms,nt1) * tc(is2,nt)

                     mtflx(ms,nt) = mtflx(ms,nt)   &
                                   + triarea(ms,ng) * mtsum(ms,nt)
                  enddo         ! ij

               endif            ! tracer type
            enddo               ! ntrace
         endif                  ! present(mtflx)
      enddo                     ! ng

      end subroutine transport_integrals

!=======================================================================
!
! Given transports through cell edges, compute new area and tracers.
!
! author William H. Lipscomb, LANL

      subroutine update_fields (ncells,      skipme,     &
                                nbr,         ijslot,     &
                                nx_block,                &
                                ntrace,                  &
                                tracer_type, depend,     &
                                tarear,                  &
                                l_stop,                  &
                                istop,       jstop,      &
                                mflxe,       mflxn,      &
                                mm,                      &
                                mtflxe,      mtflxn,     &
                                tm)

      use ice_grid, only: nbr_W, nbr_S

      integer (kind=int_kind), intent(in) ::   &
         ncells            , & ! number of slots this block holds
         nx_block          , & ! block width, to report (i,j) on abort
         ntrace                ! number of tracers in use

      ! ⚠ The block form ran over jlo..jhi x ilo..ihi, the physical domain,
      ! NOT the whole block.  A slot loop has no bounds of its own, so the
      ! range has to be carried explicitly or it silently widens to the ghost
      ! rim -- the fault class 6681fca cost.  skipPhys is exactly that range.
      logical (kind=log_kind), dimension (:), intent(in) ::   &
         skipme                ! .true. outside the physical domain

      ! ⚠ dynNbr resolves a missing neighbour to the DEAD SLOT, which holds
      ! c0, so it already does what the "if (mw < 1) mw = mdead" guard did
      ! around the old map lookup.
      integer (kind=int_kind), dimension (:,:), intent(in) ::   &
         nbr                   ! neighbour slot in each compass direction

      integer (kind=int_kind), dimension (:), intent(in) ::   &
         ijslot                ! flattened (i,j) of each slot, to report an abort

      integer (kind=int_kind), dimension (ntrace), intent(in) ::     &
         tracer_type       , & ! = 1, 2, or 3 (see comments above)
         depend              ! tracer dependencies (see above)

      real (kind=dbl_kind), dimension (:), intent(in), contiguous ::   &
         mflxe, mflxn     ! mass transport across east and north cell edges

      ! ⚠ read at this cell only, never at a neighbour -- tarea is nonzero
      ! over land, so summing it across slots would be the uarea bug again.
      real (kind=dbl_kind), dimension (:), intent(in) ::   &
         tarear           ! 1/tarea

      real (kind=dbl_kind), dimension (:), intent(inout) ::   &
         mm               ! mass field (mean)

      real (kind=dbl_kind), dimension (:,:), intent(in), optional, contiguous ::   &
         mtflxe, mtflxn   ! mass*tracer transport across E and N cell edges

      real (kind=dbl_kind), dimension (:,:), intent(inout), optional ::   &
         tm               ! tracer fields

      logical (kind=log_kind), intent(inout) ::   &
         l_stop           ! if true, abort on return

      integer (kind=int_kind), intent(inout) ::   &
         istop, jstop     ! indices of grid cell where model aborts

      ! local variables

      integer (kind=int_kind) ::   &
         is, mw, msth   , & ! slot of this cell and of its W and S neighbours
         isstop         , & ! slot where the mass went negative
         ijb            , & ! flattened (i,j) of that slot
         nt, nt1, nt2     ! tracer indices

      real (kind=dbl_kind), dimension(ncells,ntrace) ::   &
         mtold            ! old mass*tracer

      real (kind=dbl_kind) ::   &
         puny, &          !
         w1               ! work variable

      character(len=*), parameter :: subname = '(update_fields)'

      !-------------------------------------------------------------------
      ! Save starting values of mass*tracer
      !-------------------------------------------------------------------

      call icepack_query_parameters(puny_out=puny)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      if (present(tm)) then
         do nt = 1, ntrace
            if (tracer_type(nt)==1) then ! does not depend on other tracers
               do is = 1, ncells
                  if (skipme(is)) cycle
                  mtold(is,nt) = mm(is) * tm(is,nt)
               enddo            ! is
            elseif (tracer_type(nt)==2) then  ! depends on another tracer
               nt1 = depend(nt)
               do is = 1, ncells
                  if (skipme(is)) cycle
                  mtold(is,nt) = mm(is) * tm(is,nt1) * tm(is,nt)
               enddo            ! is
            elseif (tracer_type(nt)==3) then  ! depends on two tracers
               nt1 = depend(nt)
               nt2 = depend(nt1)
               do is = 1, ncells
                  if (skipme(is)) cycle
                  mtold(is,nt) = mm(is)    &
                            * tm(is,nt2) * tm(is,nt1) * tm(is,nt)
               enddo            ! is
            endif               ! depend(nt) = 0
         enddo                  ! nt
      endif                     ! present(tm)

      !-------------------------------------------------------------------
      ! Update mass field
      !-------------------------------------------------------------------

      isstop = 0

      do is = 1, ncells
         if (skipme(is)) cycle

         mw   = nbr(is,nbr_W)
         msth = nbr(is,nbr_S)
         w1 = mflxe(is) - mflxe(mw)   &
            + mflxn(is) - mflxn(msth)
         mm(is) = mm(is) - w1*tarear(is)

         if (mm(is) < -puny) then    ! abort with negative value
            l_stop = .true.
            isstop = is
         elseif (mm(is) < c0) then   ! set to zero
            mm(is) = c0
         endif

      enddo

      if (l_stop) then
         is = isstop
         ijb = ijslot(is)
         istop = mod(ijb-1, nx_block) + 1
         jstop = (ijb-1)/nx_block + 1
         mw   = nbr(is,nbr_W)
         msth = nbr(is,nbr_S)
         w1 = mflxe(is) - mflxe(mw)   &
            + mflxn(is) - mflxn(msth)
         write (nu_diag,*) ' '
         write (nu_diag,*) 'New mass < 0, i, j =', istop, jstop
         write (nu_diag,*) 'Old mass =', mm(is) + w1*tarear(is)
         write (nu_diag,*) 'New mass =', mm(is)
         write (nu_diag,*) 'Net transport =', -w1*tarear(is)
         return
      endif

      !-------------------------------------------------------------------
      ! Update tracers
      !-------------------------------------------------------------------

      if (present(tm)) then

         ! The icells list is gone: it selected mm > c0 over exactly the
         ! range skipme already carries, so the test is inlined per slot.

         do nt = 1, ntrace

            do is = 1, ncells
               if (skipme(is)) cycle
               tm(is,nt) = c0
            enddo

            if (tracer_type(nt)==1) then ! does not depend on other tracers

               do is = 1, ncells
                  if (skipme(is)) cycle
                  if (mm(is) <= c0) cycle

                  mw   = nbr(is,nbr_W)
                  msth = nbr(is,nbr_S)
                  w1  = mtflxe(is,nt) - mtflxe(mw,nt)   &
                      + mtflxn(is,nt) - mtflxn(msth,nt)
                  tm(is,nt) = (mtold(is,nt) - w1*tarear(is))   &
                                / mm(is)
               enddo            ! is

            elseif (tracer_type(nt)==2) then ! depends on another tracer
               nt1 = depend(nt)

               do is = 1, ncells
                  if (skipme(is)) cycle
                  if (mm(is) <= c0) cycle

                  if (abs(tm(is,nt1)) > c0) then
                     mw   = nbr(is,nbr_W)
                     msth = nbr(is,nbr_S)
                     w1  = mtflxe(is,nt) - mtflxe(mw,nt)   &
                         + mtflxn(is,nt) - mtflxn(msth,nt)
                     tm(is,nt) = (mtold(is,nt) - w1*tarear(is))   &
                                 / (mm(is) * tm(is,nt1))
                  endif

               enddo            ! is

            elseif (tracer_type(nt)==3) then ! depends on two tracers
               nt1 = depend(nt)
               nt2 = depend(nt1)

               do is = 1, ncells
                  if (skipme(is)) cycle
                  if (mm(is) <= c0) cycle

                  if (abs(tm(is,nt1)) > c0 .and.   &
                      abs(tm(is,nt2)) > c0) then
                     mw   = nbr(is,nbr_W)
                     msth = nbr(is,nbr_S)
                     w1  = mtflxe(is,nt) - mtflxe(mw,nt)   &
                         + mtflxn(is,nt) - mtflxn(msth,nt)
                     tm(is,nt) = (mtold(is,nt) - w1*tarear(is))   &
                              / (mm(is) * tm(is,nt2) * tm(is,nt1))
                  endif
               enddo            ! is

            endif               ! tracer_type
         enddo                  ! nt
      endif                     ! present(tm)

      end subroutine update_fields

!=======================================================================

      end module ice_transport_remap

!=======================================================================
