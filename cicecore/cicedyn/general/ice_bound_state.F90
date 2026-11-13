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

      module ice_bound_state

      use ice_kinds_mod
      use ice_constants, only: field_loc_center, field_type_scalar, c0
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_diag
!      use icepack_intfc, only: icepack_query_tracer_sizes
!      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

      implicit none
      private
      public :: bound_state, pack_state_slots

!=======================================================================

      contains

!=======================================================================
!
! Get ghost cell values for ice state variables in each thickness category.
! NOTE: This subroutine cannot be called from inside a block loop!
!
! author: William H. Lipscomb, LANL

      subroutine bound_state (aicen,        &
                              vicen, vsnon, &
                              ntrcr, trcrn, &
                              restore, slots_current)

      use ice_boundary, only: ice_halo, ice_HaloMask, ice_HaloUpdate, &
          ice_HaloDestroy, ice_HaloUpdate3DR8_compressed
      use ice_domain_size, only: max_blocks, ncat
      use ice_blocks, only: nx_block, ny_block
      use ice_domain, only: halo_info, maskhalo_bound, nblocks
      use ice_grid, only: dynMap, dynOff,       &
          pack_cat_to_slots, pack_trc_to_slots, &
          unpack_cat_from_slots, unpack_trc_from_slots
      use ice_state, only: aicen_slot, vicen_slot, vsnon_slot, trcrn_slot, &
          trcrn_flat
      use ice_restoring, only: ice_restoring_halo

      integer (kind=int_kind), intent(in) :: &
         ntrcr     ! number of tracers in use

      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat,max_blocks), intent(inout) :: &
         aicen , & ! fractional ice area
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), intent(inout), dimension(:,:,:,:,:) :: &  ! (nx_block,ny_block,ntrcr,ncat,max_blocks)
         trcrn     ! ice tracers

      logical (kind=log_kind), intent(in), optional :: &
         restore      , & !
         ! ⚠ PHASE, not a tuning switch.  The compressed branch packs the
         ! block state into the slots before the halo exchange, which is
         ! needed at INIT -- init_state calls this before anything has filled
         ! the slots.  Every per-step caller has just unpacked the slots it
         ! computed on, so there the pack copies a value onto itself.  Those
         ! callers say so, and the pack is skipped.
         slots_current    ! .true. if the slots already hold the state

      ! local variables

      integer (kind=int_kind) :: i, j, n, iblk

      integer (kind=int_kind), &
         dimension(nx_block,ny_block,max_blocks) :: halomask

      type (ice_halo) :: halo_info_aicemask  ! halo mask

      logical (kind=log_kind) :: lrestore  ! local restart flag
      logical (kind=log_kind) :: lpack     ! pack block -> slots first

      character(len=*), parameter :: subname = '(bound_state)'

      if (present(restore)) then
         lrestore = restore
      else
         lrestore = .true.
      endif

      ! Default is to pack: a caller that has not thought about it gets the
      ! safe behaviour, and only one that knows the slots are current opts out.
      if (present(slots_current)) then
         lpack = .not. slots_current
      else
         lpack = .true.
      endif

      ! ⚠ The non-masked branch halos the COMPRESSED forms.  That is only
      ! bit-for-bit because init_ragged_layout forces a 2*nghost frame:
      ! the compressed halo writes c0 when a copy's SOURCE has no slot,
      ! and trcrn carries a non-zero default at land, so without the
      ! frame 5500 ghost values differed on evpB.  See the comment there
      ! before narrowing that frame.
      !
      ! The block forms are still authoritative for every unconverted
      ! consumer, so the slots are unpacked back afterwards.  That makes
      ! this a pack + halo + unpack where it used to be a halo + pack --
      ! the same number of passes.  The saving arrives when the exit
      ! unpacks in ice_step_mod can be gated, which this unblocks.
      if (maskhalo_bound) then
         ! ⚠ This branch halos the BLOCK arrays and never touches the slots,
         ! so it needs the block form to be current on entry.  The step_*
         ! exit unpacks used to provide that; they are gone, because on the
         ! compressed path below bound_state re-packs anyway.  Refresh here
         ! instead, so the two branches do not depend on each other.
         if (.not. lpack) then
            !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
            do iblk = 1, nblocks
               call unpack_cat_from_slots(aicen_slot, dynMap(:,iblk), aicen(:,:,:,iblk))
               call unpack_cat_from_slots(vicen_slot, dynMap(:,iblk), vicen(:,:,:,iblk))
               call unpack_cat_from_slots(vsnon_slot, dynMap(:,iblk), vsnon(:,:,:,iblk))
               call unpack_trc_from_slots(trcrn_slot, dynMap(:,iblk), trcrn(:,:,:,:,iblk))
            enddo
            !$OMP END PARALLEL DO
         endif

         call ice_HaloUpdate (aicen,            halo_info, &
                              field_loc_center, field_type_scalar)
         halomask(:,:,:) = 0

         !$OMP PARALLEL DO PRIVATE(iblk,n,i,j)
         do iblk = 1, nblocks
         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            if (aicen(i,j,n,iblk) > c0) halomask(i,j,iblk) = 1
         enddo
         enddo
         enddo
         enddo
         !$OMP END PARALLEL DO

         call ice_HaloMask(halo_info_aicemask, halo_info, halomask)

         call ice_HaloUpdate (trcrn(:,:,:,:,:), halo_info_aicemask, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vicen,            halo_info_aicemask, &
                              field_loc_center, field_type_scalar)
         call ice_HaloUpdate (vsnon,            halo_info_aicemask, &
                              field_loc_center, field_type_scalar)
         call ice_HaloDestroy(halo_info_aicemask)

      else
         if (lpack) call pack_state_slots

         call ice_HaloUpdate3DR8_compressed(aicen_slot, halo_info,      &
                  field_loc_center, field_type_scalar, dynMap, dynOff)
         call ice_HaloUpdate3DR8_compressed(vicen_slot, halo_info,      &
                  field_loc_center, field_type_scalar, dynMap, dynOff)
         call ice_HaloUpdate3DR8_compressed(vsnon_slot, halo_info,      &
                  field_loc_center, field_type_scalar, dynMap, dynOff)
         ! trcrn_slot is (slot,ntrcr,ncat); trcrn_flat is its contiguous
         ! (slot, ntrcr*ncat) view, so all levels go in ONE call.
         call ice_HaloUpdate3DR8_compressed(trcrn_flat, halo_info,      &
                  field_loc_center, field_type_scalar, dynMap, dynOff)

         !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
         do iblk = 1, nblocks
            call unpack_cat_from_slots(aicen_slot, &
                 dynMap(:,iblk), aicen(:,:,:,iblk))
            call unpack_cat_from_slots(vicen_slot, &
                 dynMap(:,iblk), vicen(:,:,:,iblk))
            call unpack_cat_from_slots(vsnon_slot, &
                 dynMap(:,iblk), vsnon(:,:,:,iblk))
            call unpack_trc_from_slots(trcrn_slot, &
                 dynMap(:,iblk), trcrn(:,:,:,:,iblk))
         enddo
         !$OMP END PARALLEL DO
      endif

      if (lrestore) then
         call ice_restoring_halo(setfld='state')
      endif

      ! ⚠ The halo above touched the BLOCK forms only.  The compressed
      ! copies carry the same state for every consumer that has converted
      ! -- transport reads aicen_slot and its remap is a +/-1 stencil, so a
      ! stale ghost slot is a wrong answer, not a rounding difference -- so
      ! they are refreshed here rather than at each call site.
      !
      ! Transitional.  The end state is to halo the slot forms directly:
      ! ice_HaloUpdate3DR8_compressed takes (slot,bundle), so ncat levels
      ! go in ONE call and trcrn through its trcrn_flat view.  That needs
      ! every reader of the block forms converted first, and the
      ! maskhalo_bound branch needs its mask built from the slots.
      call pack_state_slots

      end subroutine bound_state

!=======================================================================
!
! Refresh the compressed copies of the state from the block forms.
!
! bound_state calls this as its last act, so the slots track the block
! arrays through the timestep.  ⚠ It is also needed OUTSIDE the step
! loop: init_restart populates trcrn for fsd, bgc, aerosols, isotopes
! and brine AFTER init_state has already run bound_state, so without a
! refresh there the slot copies miss those tracers entirely and the
! initial-condition history writes zeros for every field computed from
! them.  That is a silent wrong answer, not a crash.

      subroutine pack_state_slots

      use ice_domain, only: nblocks
      use ice_grid, only: dynMap, dynOff,                      &
          pack_cat_to_slots, pack_trc_to_slots, pack_to_slots, &
          unpack_cat_from_slots, unpack_trc_from_slots
      use ice_state, only: aicen, vicen, vsnon, trcrn,    &
          aicen_slot, vicen_slot, vsnon_slot, trcrn_slot, &
          aice, vice, vsno, aice0, trcr,                  &
          aice_slot, vice_slot, vsno_slot, aice0_slot, trcr_slot
      use ice_state, only: uvel, vvel, uvelE, vvelE, uvelN, vvelN
      use ice_dyn_shared, only: uvel_slot, vvel_slot, &
          uvelE_slot, vvelE_slot, uvelN_slot, vvelN_slot

      integer (kind=int_kind) :: iblk

      character(len=*), parameter :: subname = '(pack_state_slots)'

      !$OMP PARALLEL DO PRIVATE(iblk) SCHEDULE(runtime)
      do iblk = 1, nblocks
         call pack_cat_to_slots(aicen(:,:,:,iblk), dynMap(:,iblk), &
              aicen_slot)
         call pack_cat_to_slots(vicen(:,:,:,iblk), dynMap(:,iblk), &
              vicen_slot)
         call pack_cat_to_slots(vsnon(:,:,:,iblk), dynMap(:,iblk), &
              vsnon_slot)
         call pack_trc_to_slots(trcrn(:,:,:,:,iblk), dynMap(:,iblk), &
              trcrn_slot)
         ! ⚠ The aggregates too.  update_state is their only other
         ! producer and CICE_RunMod is its only caller, so without this
         ! they are c0 at the initial-condition write and every history
         ! field guarded by aice silently reads zero.
         call pack_to_slots(aice (:,:,iblk), dynMap(:,iblk), &
              aice_slot)
         call pack_to_slots(vice (:,:,iblk), dynMap(:,iblk), &
              vice_slot)
         call pack_to_slots(vsno (:,:,iblk), dynMap(:,iblk), &
              vsno_slot)
         call pack_to_slots(aice0(:,:,iblk), dynMap(:,iblk), &
              aice0_slot)
         call pack_cat_to_slots(trcr(:,:,:,iblk), dynMap(:,iblk), &
              trcr_slot)

         ! ⚠ The velocities too.  Every producer of uvel_slot is inside a
         ! solver, so on a continue run nothing has packed them when the
         ! initial-condition history is written and accum_hist -- which
         ! takes n_uvel and n_vvel from the slots -- records zero.
         ! allocated() covers kdyn <= 0, where init_dyn_shared never ran.
         if (allocated(uvel_slot)) then
            call pack_to_slots(uvel(:,iblk), dynMap(:,iblk), &
                 uvel_slot)
            call pack_to_slots(vvel(:,iblk), dynMap(:,iblk), &
                 vvel_slot)
         endif

         ! Allocated for grid_ice C and CD only.
         if (allocated(uvelE_slot)) then
            call pack_to_slots(uvelE(:,:,iblk), dynMap(:,iblk), &
                 uvelE_slot)
            call pack_to_slots(vvelE(:,:,iblk), dynMap(:,iblk), &
                 vvelE_slot)
            call pack_to_slots(uvelN(:,:,iblk), dynMap(:,iblk), &
                 uvelN_slot)
            call pack_to_slots(vvelN(:,:,iblk), dynMap(:,iblk), &
                 vvelN_slot)
         endif
      enddo
      !$OMP END PARALLEL DO

      end subroutine pack_state_slots

!=======================================================================

      end module ice_bound_state

!=======================================================================
