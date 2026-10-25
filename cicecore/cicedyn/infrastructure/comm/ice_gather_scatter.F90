!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ice_gather_scatter

!  This module contains routines for gathering data to a single
!  processor from a distributed array, and scattering data from a
!  single processor to a distributed array.
!
!  NOTE: The arrays gathered and scattered are assumed to have
!        horizontal dimensions (nx_block, ny_block).
!
! author: Phil Jones, LANL
! Oct. 2004: Adapted from POP version by William H. Lipscomb, LANL
! Jan. 2008: Elizabeth Hunke replaced old routines with new POP
!              infrastructure, added specialized routine scatter_global_stress
! Jul, 2026: T. Craig significant cleanup

#ifndef NO_MPI
   use mpi   ! MPI Fortran module
#endif
   use ice_kinds_mod
   use ice_communicate, only: my_task
#ifndef NO_MPI
   use ice_communicate, only: mpiR8, mpiR4, mpitag_gs, MPI_COMM_ICE, &
       ice_barrier, add_mpi_barriers
#endif
   use ice_constants, only: spval_dbl, spval, spval_int, c0, c1, &
       field_loc_center, field_loc_NEcorner, field_loc_Nface, field_loc_Eface, &
       field_loc_noupdate, &
       field_type_scalar, field_type_vector, field_type_angle, &
       field_type_noupdate
   use ice_blocks, only: block, nx_block, ny_block, nblocks_tot, get_block, &
       nblocks_x, nblocks_y, nghost, i_global, j_global, &
       ew_boundary_type, ns_boundary_type
   use ice_distribution, only: distrb
   use ice_domain_size, only: nx_global, ny_global
   use ice_exit, only: abort_ice
   use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted

   implicit none
   private

   public :: gather_global,      &
             scatter_global,     &
             gather_globalR8_compressed,  &
             scatter_globalR8_compressed, &
             slot_gather_setup,  &
             slot_gather_filter, &
             gather_slots_r8

!-----------------------------------------------------------------------
!
!  Descriptors for the slot-compacted gather (slot_gather_setup /
!  gather_slots_r8).  gather_globalR8_compressed expands a slot array back
!  to full block shape and sends nx_block*ny_block per block; these move
!  only the slots that actually feed the destination 1-D array.
!
!  ⚠ Which block owns which global cell -- interior plus the grid_ext rims
!  -- is NOT re-derived here.  Duplicating that rule is how the two paths
!  silently drift.  slot_gather_setup MEASURES it by pushing a field of
!  unique (block,slot) ids through gather_globalR8_compressed itself and
!  reading back where each one landed, so the compacted path is correct by
!  construction for whatever that routine does, rims and corners included.
!
!-----------------------------------------------------------------------

   integer (int_kind), allocatable, target, dimension(:,:) :: &
      sg_sendIdx       ! (slot, local blk) index into the caller's slot array
   integer (int_kind), allocatable, target, dimension(:)   :: &
      sg_nSend         ! (local blk) how many this task sends for that block
   integer (int_kind), allocatable, target, dimension(:,:) :: &
      sg_recvIW        ! (slot, global blk) destination 1-D index, dst_task only
   integer (int_kind), allocatable, target, dimension(:)   :: &
      sg_nRecv         ! (global blk) how many arrive for that block, dst only
   integer (int_kind), allocatable, dimension(:,:) :: &
      sg_sendIJ        ! (slot, local blk) block-local cell of each send slot,
                       ! so a task can test its own mask without the master
   integer (int_kind), allocatable, target, dimension(:,:) :: &
      sg_fSend         ! (slot, local blk) the surviving subset of sg_sendIdx
   integer (int_kind), allocatable, target, dimension(:)   :: &
      sg_nfSend        ! (local blk) how many survive
   integer (int_kind), allocatable, target, dimension(:,:) :: &
      sg_fRecvIW       ! (slot, global blk) matching destinations, dst only
   integer (int_kind), allocatable, target, dimension(:)   :: &
      sg_nfRecv        ! (global blk) how many survive, dst only
   integer (int_kind) :: &
      sg_maxSend = 0   ! longest message, both sides
   logical (log_kind) :: &
      sg_ready = .false., &
      sg_firstFilter = .true., &
      sg_filter = .false.  ! are the filtered lists the live ones?

!-----------------------------------------------------------------------
!
!  overload module functions
!
!-----------------------------------------------------------------------

   interface gather_global
     module procedure gather_globalR8, &
                      gather_globalR4, &
                      gather_globalI4, &
                      gather_globalL1, &
                      gather_globalR8_flat, &
                      gather_globalL1_flat
   end interface

   interface scatter_global
     module procedure scatter_globalR8, &
                      scatter_globalR4, &
                      scatter_globalI4, &
                      scatter_globalL1, &
                      scatter_globalR8_flat
   end interface

!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

!***********************************************************************

 contains

!***********************************************************************
! autogenerate start
!***********************************************************************

 subroutine gather_globalR8(ARRAY_G, ARRAY, dst_task, src_dist, fillValue, grid_ext)

!  Generated by ./generate_GatherScatter.sh on 2026-08-27
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!  If grid_ext is true, then ARRAY_G is extended grid
!  If grid_ext is false, then ARRAY_G is 1:nx_global,1:ny_global

   integer (int_kind), intent(in) :: &
      dst_task   ! task to which array should be gathered

   type (distrb), intent(in) :: &
      src_dist   ! distribution of blocks in the source array

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
      ARRAY      ! array containing horizontal slab of distributed field

   real (dbl_kind), intent(in), optional :: &
      fillValue  ! special value

   logical (log_kind), intent(in), optional :: &
      grid_ext   ! extended global grid flag

   real (dbl_kind), dimension(:,:), intent(inout) :: &
      ARRAY_G    ! array containing global horizontal field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,n,          &! dummy loop counters
      iblock,jblock,  &! block index
      ilo,ihi,jlo,jhi  ! block sizes

#ifndef NO_MPI
   integer (int_kind) :: &
      nsends,         &! number of actual sends
      src_block,      &! block locator for send
      ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
      status

   integer (int_kind), dimension(:), allocatable :: &
      snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
      snd_status

   real (dbl_kind), dimension(:,:), allocatable :: &
      msg_buffer
#endif

   real (dbl_kind) :: &
      lfillValue  ! local special value

   logical (log_kind) :: &
      lext        ! local grid_ext global grid flag

   type (block) :: &
      this_block  ! block info for current block

   real (dbl_kind), dimension(:,:), allocatable :: &
      LARRAY_G    ! array containing global horizontal field on dst_task
                  ! but with indexing that matches nx_global, ny_global
                  ! to simplify overlapping ext and non-ext coding.
                  ! for non-ext, ARRAY_G and LARRAY_G are (nx_global,ny_global).
                  ! for ext incoming ARRAY_G is (1:nx_global+2*nghost,1:ny_global+2*nghost)
                  ! LARRAY_G is same size, but indexing is
                  !            (1-nghost:nx_global+nghost,1-nghost:ny_global+nghost)
                  ! Could maybe use a pointer but a separate array and copy is more robust

   character(len=*), parameter :: subname = '(gather_globalR8)'

   if (present(fillValue)) then
      lfillValue = fillValue
   else
      lfillValue = spval_dbl
   endif

   if (present(grid_ext)) then
      lext = grid_ext
   else
      lext = .false.
   endif

!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (my_task == dst_task) then

      if (lext) then
         allocate(LARRAY_G(1-nghost:nx_global+nghost,1-nghost:ny_global+nghost))
         LARRAY_G(1-nghost:nx_global+nghost,1-nghost:ny_global+nghost) = ARRAY_G(1:nx_global+2*nghost,1:ny_global+2*nghost)
      else
         allocate(LARRAY_G(1:nx_global,1:ny_global))
         LARRAY_G(:,:) = ARRAY_G(:,:)
      endif

      do n=1,nblocks_tot

         this_block = get_block(n,n)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         iblock = this_block%iblock
         jblock = this_block%jblock

         !*** copy local blocks

         if (src_dist%blockLocation(n) == my_task+1) then

            do j = jlo, jhi
            do i = ilo, ihi
               LARRAY_G(this_block%i_glob(i),this_block%j_glob(j)) = ARRAY(i,j,src_dist%blockLocalID(n))
            end do
            end do

            !*** fill outside halo for extended grid gather

            if (lext) then
               ! south
               if (jblock == 1) then
                  do j = 1, nghost
                  do i = ilo, ihi
                     LARRAY_G(this_block%i_glob(i),j-nghost) = ARRAY(i,j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! north
               if (jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = ilo, ihi
                     LARRAY_G(this_block%i_glob(i),ny_global+j) = ARRAY(i,jhi+j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! west
               if (iblock == 1) then
                  do j = jlo, jhi
                  do i = 1, nghost
                     LARRAY_G(i-nghost,this_block%j_glob(j)) = ARRAY(i,j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! east
               if (iblock == nblocks_x) then
                  do j = jlo, jhi
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,this_block%j_glob(j)) = ARRAY(ihi+i,j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! southwest
               if (iblock == 1 .and. jblock == 1) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(i-nghost,j-nghost) = ARRAY(i,j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! northwest
               if (iblock == 1 .and. jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(i-nghost,ny_global+j) = ARRAY(i,jhi+j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! southeast
               if (iblock == nblocks_x .and. jblock == 1) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,j-nghost) = ARRAY(ihi+i,j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! northeast
               if (iblock == nblocks_x .and. jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,ny_global+j) = ARRAY(ihi+i,jhi+j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif
            endif

         !*** fill land blocks with special values

         else if (src_dist%blockLocation(n) == 0) then

            do j = jlo, jhi
            do i = ilo, ihi
               LARRAY_G(this_block%i_glob(i),this_block%j_glob(j)) = lfillValue
            end do
            end do

         endif

      end do

#ifndef NO_MPI
      !*** copy blocks from other tasks

      allocate (msg_buffer(nx_block,ny_block))

      do n=1,nblocks_tot
         if (src_dist%blockLocation(n) > 0 .and. &
             src_dist%blockLocation(n) /= my_task+1) then

            this_block = get_block(n,n)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            iblock = this_block%iblock
            jblock = this_block%jblock

            call MPI_RECV(msg_buffer, size(msg_buffer), &
                 mpiR8, src_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                 MPI_COMM_ICE, status, ierr)

            do j = jlo, jhi
            do i = ilo, ihi
               LARRAY_G(this_block%i_glob(i),this_block%j_glob(j)) = msg_buffer(i,j)
            end do
            end do

            !*** fill outside halo for extended grid gather

            if (lext) then
               ! south
               if (jblock == 1) then
                  do j = 1, nghost
                  do i = ilo, ihi
                     LARRAY_G(this_block%i_glob(i),j-nghost) = msg_buffer(i,j)
                  enddo
                  enddo
               endif

               ! north
               if (jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = ilo, ihi
                     LARRAY_G(this_block%i_glob(i),ny_global+j) = msg_buffer(i,jhi+j)
                  enddo
                  enddo
               endif

               ! west
               if (iblock == 1) then
                  do j = jlo, jhi
                  do i = 1, nghost
                     LARRAY_G(i-nghost,this_block%j_glob(j)) = msg_buffer(i,j)
                  enddo
                  enddo
               endif

               ! east
               if (iblock == nblocks_x) then
                  do j = jlo, jhi
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,this_block%j_glob(j)) = msg_buffer(ihi+i,j)
                  enddo
                  enddo
               endif

               ! southwest
               if (iblock == 1 .and. jblock == 1) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(i-nghost,j-nghost) = msg_buffer(i,j)
                  enddo
                  enddo
               endif

               ! northwest
               if (iblock == 1 .and. jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(i-nghost,ny_global+j) = msg_buffer(i,jhi+j)
                  enddo
                  enddo
               endif

               ! southeast
               if (iblock == nblocks_x .and. jblock == 1) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,j-nghost) = msg_buffer(ihi+i,j)
                  enddo
                  enddo
               endif

               ! northeast
               if (iblock == nblocks_x .and. jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,ny_global+j) = msg_buffer(ihi+i,jhi+j)
                  enddo
                  enddo
               endif
            endif ! lext
         endif
      end do

      deallocate(msg_buffer)
#endif

      ! copy LARRAY_G to ARRAY_G, offset extended grid index

      if (lext) then
         ARRAY_G(1:nx_global+2*nghost,1:ny_global+2*nghost) = LARRAY_G(1-nghost:nx_global+nghost,1-nghost:ny_global+nghost)
      else
         ARRAY_G(1:nx_global,1:ny_global) = LARRAY_G(1:nx_global,1:ny_global)
      endif
      deallocate(LARRAY_G)

#ifndef NO_MPI
!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else

      allocate(snd_request(nblocks_tot), &
               snd_status (MPI_STATUS_SIZE, nblocks_tot))

      nsends = 0
      do n=1,nblocks_tot
         if (src_dist%blockLocation(n) == my_task+1) then
            nsends = nsends + 1
            src_block = src_dist%blockLocalID(n)
            call MPI_ISEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                 mpiR8, dst_task, 3*mpitag_gs+n, &
                 MPI_COMM_ICE, snd_request(nsends), ierr)
         endif
      end do

      if (nsends > 0) call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
      deallocate(snd_request, snd_status)
#endif

   endif

#ifndef NO_MPI
   if (add_mpi_barriers) then
      call ice_barrier()
   endif
#endif

!-----------------------------------------------------------------------

 end subroutine gather_globalR8

!***********************************************************************

 subroutine scatter_globalR8(ARRAY, ARRAY_G, src_task, dst_dist, &
                             field_loc, field_type, fillValue, &
                             grid_ext)

!  Generated by ./generate_GatherScatter.sh on 2026-08-27
!  This subroutine scatters a global sized array from src_task to all
!  processors.
!  If grid_ext is true, then ARRAY_G is extended grid and updates outer halo
!  If grid_ext is false, then ARRAY_G is 1:nx_global,1:ny_global and outer
!    halo is updated for cyclic and tripole bcs.
!  If grid_ext is false and interior is true, then only the interior data will 
!    be scattered, no outer halo bcs will be applied

   integer (int_kind), intent(in) :: &
      src_task        ! task from which array should be scattered

   type (distrb), intent(in) :: &
      dst_dist        ! distribution of resulting blocks

   real (dbl_kind), dimension(:,:), intent(in) :: &
      ARRAY_G         ! array containing global field on src_task

   integer (int_kind), intent(in), optional :: &
      field_type,    &! id for type of field (scalar, vector, angle)
      field_loc       ! id for location on horizontal grid
                      !  (center, NEcorner, Nface, Eface)

   real (dbl_kind), intent(in), optional :: &
      fillValue       ! special value

   logical (log_kind), intent(in), optional :: &
      grid_ext        ! extended grid flag

   real (dbl_kind), dimension(:,:,:), intent(inout) :: &
      ARRAY           ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,n,         &! dummy loop indices
      ilo,ihi,jlo,jhi,&! block indices
      isrc, jsrc,    &! source addresses
      ioffset, joffset,   &! offsets for tripole boundary conditions
      isign,         &! sign factor for tripole boundary conditions
      dst_block       ! location of block in dst array

   real (dbl_kind) :: &
      efac            ! linear extrapolation factor

   type (block) :: &
      this_block      ! block info for current block

   real (dbl_kind) :: &
      lfillValue      ! local special value

   logical (log_kind) :: &
      lext,          &! local extended global grid flag
      noupdate        ! flag if noupdate is on

   real (dbl_kind), dimension(:,:), allocatable :: &
      LARRAY_G    ! array containing global horizontal field on dst_task
                  ! but with indexing that matches nx_global, ny_global
                  ! to simplify overlapping ext and non-ext coding.
                  ! for non-ext, ARRAY_G and LARRAY_G are (nx_global,ny_global).
                  ! for ext incoming ARRAY_G is (1:nx_global+2*nghost,1:ny_global+2*nghost)
                  ! LARRAY_G is same size, but indexing is
                  !            (1-nghost:nx_global+nghost,1-nghost:ny_global+nghost)
                  ! Could maybe use a pointer but a separate array and copy is more robust

#ifndef NO_MPI
   integer (int_kind) :: &
      nrecvs,        &! actual number of messages received
      ierr            ! MPI error flag

   integer (int_kind), dimension(:), allocatable :: &
      rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
      rcv_status      ! status array for receives

   real (dbl_kind), dimension(:,:), allocatable :: &
      msg_buffer      ! buffer for sending blocks
#endif

   character(len=*), parameter :: subname = '(scatter_globalR8)'

   if (present(fillValue)) then
      lfillValue = fillValue
   else
      lfillValue = c0
   endif

   if (present(grid_ext)) then
      lext = grid_ext
   else
      lext = .false.
   endif

   noupdate = .false.
   if (present(field_loc)) then
      if (field_loc == field_loc_noupdate) then
         noupdate = .true.
      endif
   endif
   if (present(field_type)) then
      if (field_type == field_type_noupdate) then
         noupdate = .true.
      endif
   endif

   ! check size of global grid
   if (my_task == src_task) then
      if (lext) then
         if (size(ARRAY_G,dim=1)/= nx_global+2*nghost .or. size(ARRAY_G,dim=2) /= ny_global+2*nghost) then
            call abort_ice(subname//'ERROR: grid_ext = true, wrong ARRAY_G size')
         endif
      else
         if (size(ARRAY_G,dim=1)/= nx_global .or. size(ARRAY_G,dim=2) /= ny_global) then
            call abort_ice(subname//'ERROR: grid_ext = false, wrong ARRAY_G size')
         endif
      endif
   endif

   ARRAY = lfillValue

   isign = -99
   ioffset = -999
   joffset = -999
   if (.not.lext .and. .not.noupdate .and. &
      (ns_boundary_type == 'tripole' .or. ns_boundary_type == 'tripoleT')) then

      select case (field_type)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case (field_type_noupdate) ! ghost cells never used - use cell center
         isign =  1
      case default
         call abort_ice(subname//'ERROR: Unknown field kind')
      end select

      if (ns_boundary_type == 'tripole') then
         ioffset = 0
         joffset = 0
      else ! tripoleT fold
         ioffset = -1
         joffset = 1
      endif

      select case (field_loc)
        case (field_loc_center)     ! cell center location
           !ioffset = ioffset
           !joffset = joffset
        case (field_loc_NEcorner)   ! cell corner (velocity) location
           ioffset = ioffset + 1
           joffset = joffset + 1
        case (field_loc_Eface)      ! cell east face
           ioffset = ioffset + 1
           !joffset = joffset
        case (field_loc_Nface)      ! cell north face
           !ioffset = ioffset
           joffset = joffset + 1
        case (field_loc_noupdate)   ! no update
           !ioffset = ioffset
           !joffset = joffset
        case default
           call abort_ice(subname//'ERROR: Unknown field location')
      end select

   endif ! not lext

!-----------------------------------------------------------------------
!
!  if this task is the src_task, compute updated LARRAY_G which contains
!  outer halo extended grid from ARRAY_G (extended or not), fill outer
!  halo in global array then scatter
!
!  i_global, j_global are indexed 1-nghost:n*_global+nghost
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

      allocate(LARRAY_G(1-nghost:nx_global+nghost,1-nghost:ny_global+nghost))

      LARRAY_G = lfillValue
      if (lext) then
         do j = 1-nghost, ny_global+nghost
         do i = 1-nghost, nx_global+nghost
            LARRAY_G(i,j) = ARRAY_G(i+nghost,j+nghost)
         enddo
         enddo
      else

         ! fill cells that can be directly copied

         do j = 1-nghost, ny_global+nghost
         do i = 1-nghost, nx_global+nghost
            isrc = i_global(i)
            jsrc = j_global(j)
            if (isrc >= 1 .and. isrc <= nx_global .and. &
                jsrc >= 1 .and. jsrc <= ny_global) then
               ! copy valid cell directly
               LARRAY_G(i,j) = ARRAY_G(isrc,jsrc)
            endif
         enddo
         enddo

         ! fill outer halo, only if noupdate is false
         ! do north/south boundary first (partly due to tripole feature) (1:nx_global)
         ! then do east/west including corners (1-nghost:ny_global+nghost)

         if (.not.noupdate) then
            if (ns_boundary_type == 'tripole' .or. ns_boundary_type == 'tripoleT') then
               do j = ny_global, ny_global+nghost
               ! north of active cells
               do i = 1, nx_global
                  isrc = nx_global - i + 1 - ioffset   ! ioffset = 0 for u center, ioffset = 1 for T center
                  jsrc = ny_global - (j-ny_global) - joffset + 1  ! joffset = 0 for u center, joffset = 1 for T center
                  if (isrc < 1        ) isrc = isrc + nx_global
                  if (isrc > nx_global) isrc = isrc - nx_global

                  !*** for center and Eface on u-fold, and NE corner and Nface
                  !*** on T-fold, do not need to replace
                  !*** top row of physical domain, so jsrc should be greater than j

                  if (jsrc > j) then
                     ! do nothing
                  elseif (jsrc == j) then
                     ! don't average in the scatter even though halo update does
                     ! LARRAY_G(i,j) = 0.5_dbl_kind * (ARRAY_G(i,jsrc) + isign*ARRAY_G(isrc,jsrc))
                  else
                     ! copy
                     LARRAY_G(i,j) = isign * ARRAY_G(isrc,jsrc)
                  endif
               enddo
               enddo
            else
               !--- bottom edge
               do j = 1-nghost, 0
               do i = 1, nx_global
                  if (ns_boundary_type == 'cyclic') then
                     LARRAY_G(i,j) = ARRAY_G(i,j+ny_global)
                  elseif (ns_boundary_type == 'zero_gradient') then
                     LARRAY_G(i,j) = ARRAY_G(i,1)
                  elseif (ns_boundary_type == 'linear_extrap') then
                     efac = real(2-j,dbl_kind)
                     LARRAY_G(i,j) = efac*ARRAY_G(i,1) - (efac-c1)*ARRAY_G(i,2)
                  endif
               enddo
               enddo
               !--- top edge
               do j = ny_global+1, ny_global+nghost
               do i = 1, nx_global
                  if (ns_boundary_type == 'cyclic') then
                     LARRAY_G(i,j) = ARRAY_G(i,j-ny_global)
                  elseif (ns_boundary_type == 'zero_gradient') then
                     LARRAY_G(i,j) = ARRAY_G(i,ny_global)
                  elseif (ns_boundary_type == 'linear_extrap') then
                     efac = real(j-ny_global+1,dbl_kind)
                     LARRAY_G(i,j) = efac*ARRAY_G(i,ny_global) - (efac-c1)*ARRAY_G(i,ny_global-1)
                  endif
               enddo
               enddo
            endif

            do j = 1-nghost, ny_global+nghost
            !--- left edge
            do i = 1-nghost, 0
               if (ew_boundary_type == 'cyclic') then
                  LARRAY_G(i,j) = LARRAY_G(i+nx_global,j)
               elseif (ew_boundary_type == 'zero_gradient') then
                  LARRAY_G(i,j) = LARRAY_G(1,j)
               elseif (ew_boundary_type == 'linear_extrap') then
                  efac = real(2-i,dbl_kind)
                  LARRAY_G(i,j) = efac*LARRAY_G(1,j) - (efac-c1)*LARRAY_G(2,j)
               endif
            enddo
            !--- right edge
            do i = nx_global+1, nx_global+nghost
               if (ew_boundary_type == 'cyclic') then
                  LARRAY_G(i,j) = LARRAY_G(i-nx_global,j)
               elseif (ew_boundary_type == 'zero_gradient') then
                  LARRAY_G(i,j) = LARRAY_G(nx_global,j)
               elseif (ew_boundary_type == 'linear_extrap') then
                  efac = real(i-nx_global+1,dbl_kind)
                  LARRAY_G(i,j) = efac*LARRAY_G(nx_global,j) - (efac-c1)*LARRAY_G(nx_global-1,j)
               endif
            enddo
            enddo

         endif   ! noupdate
      endif   ! lext

      !*** send non-local blocks away

#ifndef NO_MPI
      allocate (msg_buffer(nx_block,ny_block))

      do n=1,nblocks_tot
         if (dst_dist%blockLocation(n) > 0 .and. &
             dst_dist%blockLocation(n)-1 /= my_task) then

            msg_buffer = lfillValue
            this_block = get_block(n,n)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j=jlo-nghost,jhi+nghost
            do i=ilo-nghost,ihi+nghost
               ! fill entire block
               isrc = this_block%i_glob(ilo) + i - ilo
               jsrc = this_block%j_glob(jlo) + j - jlo
               msg_buffer(i,j) = LARRAY_G(isrc,jsrc)
            enddo
            enddo

            call MPI_SEND(msg_buffer, nx_block*ny_block, &
                          mpiR8, dst_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                          MPI_COMM_ICE, ierr)

         endif
      end do

      deallocate(msg_buffer)
#endif
      !*** copy any local blocks

      do n=1,nblocks_tot
         if (dst_dist%blockLocation(n) == my_task+1) then
            dst_block = dst_dist%blockLocalID(n)
            this_block = get_block(n,n)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j=jlo-nghost,jhi+nghost
            do i=ilo-nghost,ihi+nghost
               ! fill entire block
               isrc = this_block%i_glob(ilo) + i - ilo
               jsrc = this_block%j_glob(jlo) + j - jlo
               ARRAY(i,j,dst_block) = LARRAY_G(isrc,jsrc)
            enddo
            enddo
         endif
      end do

#ifndef NO_MPI
!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else   ! my_task /= src_task

      allocate (rcv_request(nblocks_tot), &
                rcv_status(MPI_STATUS_SIZE, nblocks_tot))

      rcv_request = 0
      rcv_status  = 0

      nrecvs = 0
      do n=1,nblocks_tot
         if (dst_dist%blockLocation(n) == my_task+1) then
            nrecvs = nrecvs + 1
            dst_block = dst_dist%blockLocalID(n)
            call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                          mpiR8, src_task, 3*mpitag_gs+n, &
                          MPI_COMM_ICE, rcv_request(nrecvs), ierr)
         endif
      end do

      if (nrecvs > 0) &
         call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

      deallocate(rcv_request, rcv_status)
#endif
   endif   ! my_task == src_task

#ifndef NO_MPI
   if (add_mpi_barriers) then
     call ice_barrier()
   endif
#endif

!-----------------------------------------------------------------------

 end subroutine scatter_globalR8

!***********************************************************************
!
!  NOT generated by ./generate_GatherScatter.sh -- hand-written, must be
!  preserved if that script is re-run.  Scatters into a distributed field
!  whose i and j dimensions are flattened, i.e.
!  dimension(nx_block*ny_block,max_blocks), where element (i,j,iblk) lives
!  at ((j-1)*nx_block + i, iblk).  ARRAY is aliased to rank 3 and handed to
!  scatter_globalR8, so no data is copied.
!
!***********************************************************************

 subroutine scatter_globalR8_flat(ARRAY, ARRAY_G, src_task, dst_dist, &
                                  field_loc, field_type, fillValue, &
                                  grid_ext)

   integer (int_kind), intent(in) :: &
      src_task        ! task from which array should be scattered

   type (distrb), intent(in) :: &
      dst_dist        ! distribution of resulting blocks

   real (dbl_kind), dimension(:,:), intent(in) :: &
      ARRAY_G         ! array containing global field on src_task

   integer (int_kind), intent(in), optional :: &
      field_type,    &! id for type of field (scalar, vector, angle)
      field_loc       ! id for location on horizontal grid

   real (dbl_kind), intent(in), optional :: &
      fillValue       ! special value

   logical (log_kind), intent(in), optional :: &
      grid_ext        ! extended grid flag

   ! TARGET and CONTIGUOUS are required to alias this to rank 3 below
   real (dbl_kind), dimension(:,:), intent(inout), target, contiguous :: &
      ARRAY           ! array containing distributed field

   real (dbl_kind), pointer :: &
      ARRAY_3d(:,:,:) ! rank-3 alias of ARRAY

   character(len=*), parameter :: subname = '(scatter_globalR8_flat)'

   if (size(ARRAY,dim=1) /= nx_block*ny_block) then
      call abort_ice(subname//' ERROR: first dimension is not nx_block*ny_block')
   endif

   ARRAY_3d(1:nx_block,1:ny_block,1:size(ARRAY,dim=2)) => ARRAY

   call scatter_globalR8(ARRAY_3d, ARRAY_G, src_task, dst_dist, &
                         field_loc, field_type, fillValue, grid_ext)

 end subroutine scatter_globalR8_flat

!***********************************************************************

 subroutine scatter_globalR8_compressed(ARRAY, ARRAY_G, map, off,      &
                                        src_task, dst_dist,            &
                                        field_loc, field_type,         &
                                        fillValue, grid_ext)

!  Scatter a global field directly into the compressed (slot) layout, so a
!  caller whose only consumer reads slots never needs the block array.
!
!  This cannot alias the way scatter_globalR8_flat does: the flat form is a
!  reshape of the block array, but a slot vector is a mapping through `map`,
!  so the data really has to be packed.  One block-shaped temporary is used
!  for that and discarded.

   integer (int_kind), intent(in) :: &
      src_task        ! task from which array should be scattered

   type (distrb), intent(in) :: &
      dst_dist        ! distribution of resulting blocks

   real (dbl_kind), dimension(:,:), intent(in) :: &
      ARRAY_G         ! array containing global field on src_task

   integer (int_kind), dimension(:,:), intent(in) :: &
      map             ! block-local slot of each cell, 0 if not stored

   integer (int_kind), dimension(:), intent(in) :: &
      off             ! first slot of each block, minus one

   integer (int_kind), intent(in), optional :: &
      field_type,    &! id for type of field (scalar, vector, angle)
      field_loc       ! id for location on horizontal grid

   real (dbl_kind), intent(in), optional :: &
      fillValue       ! special value

   logical (log_kind), intent(in), optional :: &
      grid_ext        ! extended grid flag

   real (dbl_kind), dimension(:), intent(inout) :: &
      ARRAY           ! compressed field, dimension(dynTot)

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (dbl_kind), dimension(:,:,:), allocatable :: &
      ARRAY_b         ! block-shaped scratch, scattered into then packed

   integer (int_kind) :: i, j, ij, iblk, ierr

   character(len=*), parameter :: subname = '(scatter_globalR8_compressed)'

!-----------------------------------------------------------------------

   allocate(ARRAY_b(nx_block,ny_block,size(map,dim=2)), stat=ierr)
   if (ierr /= 0) call abort_ice(subname//' ERROR: out of memory')
   ARRAY_b = c0

   call scatter_globalR8(ARRAY_b, ARRAY_G, src_task, dst_dist, &
                         field_loc, field_type, fillValue, grid_ext)

   !  Pack into the slots.  Cells with no slot are dropped -- they carry
   !  nothing the compressed consumer reads.
   do iblk = 1, size(map,dim=2)
      do j = 1, ny_block
      do i = 1, nx_block
         ij = (j-1)*nx_block + i
         if (map(ij,iblk) > 0) then
            ARRAY(off(iblk) + map(ij,iblk)) = ARRAY_b(i,j,iblk)
         endif
      enddo
      enddo
   enddo

   deallocate(ARRAY_b)

 end subroutine scatter_globalR8_compressed

!***********************************************************************

 subroutine gather_globalR4(ARRAY_G, ARRAY, dst_task, src_dist, fillValue, grid_ext)

!  Generated by ./generate_GatherScatter.sh on 2026-08-27
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!  If grid_ext is true, then ARRAY_G is extended grid
!  If grid_ext is false, then ARRAY_G is 1:nx_global,1:ny_global

   integer (int_kind), intent(in) :: &
      dst_task   ! task to which array should be gathered

   type (distrb), intent(in) :: &
      src_dist   ! distribution of blocks in the source array

   real (real_kind), dimension(:,:,:), intent(in) :: &
      ARRAY      ! array containing horizontal slab of distributed field

   real (real_kind), intent(in), optional :: &
      fillValue  ! special value

   logical (log_kind), intent(in), optional :: &
      grid_ext   ! extended global grid flag

   real (real_kind), dimension(:,:), intent(inout) :: &
      ARRAY_G    ! array containing global horizontal field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,n,          &! dummy loop counters
      iblock,jblock,  &! block index
      ilo,ihi,jlo,jhi  ! block sizes

#ifndef NO_MPI
   integer (int_kind) :: &
      nsends,         &! number of actual sends
      src_block,      &! block locator for send
      ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
      status

   integer (int_kind), dimension(:), allocatable :: &
      snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
      snd_status

   real (real_kind), dimension(:,:), allocatable :: &
      msg_buffer
#endif

   real (real_kind) :: &
      lfillValue  ! local special value

   logical (log_kind) :: &
      lext        ! local grid_ext global grid flag

   type (block) :: &
      this_block  ! block info for current block

   real (real_kind), dimension(:,:), allocatable :: &
      LARRAY_G    ! array containing global horizontal field on dst_task
                  ! but with indexing that matches nx_global, ny_global
                  ! to simplify overlapping ext and non-ext coding.
                  ! for non-ext, ARRAY_G and LARRAY_G are (nx_global,ny_global).
                  ! for ext incoming ARRAY_G is (1:nx_global+2*nghost,1:ny_global+2*nghost)
                  ! LARRAY_G is same size, but indexing is
                  !            (1-nghost:nx_global+nghost,1-nghost:ny_global+nghost)
                  ! Could maybe use a pointer but a separate array and copy is more robust

   character(len=*), parameter :: subname = '(gather_globalR4)'

   if (present(fillValue)) then
      lfillValue = fillValue
   else
      lfillValue = spval
   endif

   if (present(grid_ext)) then
      lext = grid_ext
   else
      lext = .false.
   endif

!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (my_task == dst_task) then

      if (lext) then
         allocate(LARRAY_G(1-nghost:nx_global+nghost,1-nghost:ny_global+nghost))
         LARRAY_G(1-nghost:nx_global+nghost,1-nghost:ny_global+nghost) = ARRAY_G(1:nx_global+2*nghost,1:ny_global+2*nghost)
      else
         allocate(LARRAY_G(1:nx_global,1:ny_global))
         LARRAY_G(:,:) = ARRAY_G(:,:)
      endif

      do n=1,nblocks_tot

         this_block = get_block(n,n)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         iblock = this_block%iblock
         jblock = this_block%jblock

         !*** copy local blocks

         if (src_dist%blockLocation(n) == my_task+1) then

            do j = jlo, jhi
            do i = ilo, ihi
               LARRAY_G(this_block%i_glob(i),this_block%j_glob(j)) = ARRAY(i,j,src_dist%blockLocalID(n))
            end do
            end do

            !*** fill outside halo for extended grid gather

            if (lext) then
               ! south
               if (jblock == 1) then
                  do j = 1, nghost
                  do i = ilo, ihi
                     LARRAY_G(this_block%i_glob(i),j-nghost) = ARRAY(i,j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! north
               if (jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = ilo, ihi
                     LARRAY_G(this_block%i_glob(i),ny_global+j) = ARRAY(i,jhi+j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! west
               if (iblock == 1) then
                  do j = jlo, jhi
                  do i = 1, nghost
                     LARRAY_G(i-nghost,this_block%j_glob(j)) = ARRAY(i,j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! east
               if (iblock == nblocks_x) then
                  do j = jlo, jhi
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,this_block%j_glob(j)) = ARRAY(ihi+i,j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! southwest
               if (iblock == 1 .and. jblock == 1) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(i-nghost,j-nghost) = ARRAY(i,j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! northwest
               if (iblock == 1 .and. jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(i-nghost,ny_global+j) = ARRAY(i,jhi+j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! southeast
               if (iblock == nblocks_x .and. jblock == 1) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,j-nghost) = ARRAY(ihi+i,j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! northeast
               if (iblock == nblocks_x .and. jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,ny_global+j) = ARRAY(ihi+i,jhi+j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif
            endif

         !*** fill land blocks with special values

         else if (src_dist%blockLocation(n) == 0) then

            do j = jlo, jhi
            do i = ilo, ihi
               LARRAY_G(this_block%i_glob(i),this_block%j_glob(j)) = lfillValue
            end do
            end do

         endif

      end do

#ifndef NO_MPI
      !*** copy blocks from other tasks

      allocate (msg_buffer(nx_block,ny_block))

      do n=1,nblocks_tot
         if (src_dist%blockLocation(n) > 0 .and. &
             src_dist%blockLocation(n) /= my_task+1) then

            this_block = get_block(n,n)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            iblock = this_block%iblock
            jblock = this_block%jblock

            call MPI_RECV(msg_buffer, size(msg_buffer), &
                 mpiR4, src_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                 MPI_COMM_ICE, status, ierr)

            do j = jlo, jhi
            do i = ilo, ihi
               LARRAY_G(this_block%i_glob(i),this_block%j_glob(j)) = msg_buffer(i,j)
            end do
            end do

            !*** fill outside halo for extended grid gather

            if (lext) then
               ! south
               if (jblock == 1) then
                  do j = 1, nghost
                  do i = ilo, ihi
                     LARRAY_G(this_block%i_glob(i),j-nghost) = msg_buffer(i,j)
                  enddo
                  enddo
               endif

               ! north
               if (jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = ilo, ihi
                     LARRAY_G(this_block%i_glob(i),ny_global+j) = msg_buffer(i,jhi+j)
                  enddo
                  enddo
               endif

               ! west
               if (iblock == 1) then
                  do j = jlo, jhi
                  do i = 1, nghost
                     LARRAY_G(i-nghost,this_block%j_glob(j)) = msg_buffer(i,j)
                  enddo
                  enddo
               endif

               ! east
               if (iblock == nblocks_x) then
                  do j = jlo, jhi
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,this_block%j_glob(j)) = msg_buffer(ihi+i,j)
                  enddo
                  enddo
               endif

               ! southwest
               if (iblock == 1 .and. jblock == 1) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(i-nghost,j-nghost) = msg_buffer(i,j)
                  enddo
                  enddo
               endif

               ! northwest
               if (iblock == 1 .and. jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(i-nghost,ny_global+j) = msg_buffer(i,jhi+j)
                  enddo
                  enddo
               endif

               ! southeast
               if (iblock == nblocks_x .and. jblock == 1) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,j-nghost) = msg_buffer(ihi+i,j)
                  enddo
                  enddo
               endif

               ! northeast
               if (iblock == nblocks_x .and. jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,ny_global+j) = msg_buffer(ihi+i,jhi+j)
                  enddo
                  enddo
               endif
            endif ! lext
         endif
      end do

      deallocate(msg_buffer)
#endif

      ! copy LARRAY_G to ARRAY_G, offset extended grid index

      if (lext) then
         ARRAY_G(1:nx_global+2*nghost,1:ny_global+2*nghost) = LARRAY_G(1-nghost:nx_global+nghost,1-nghost:ny_global+nghost)
      else
         ARRAY_G(1:nx_global,1:ny_global) = LARRAY_G(1:nx_global,1:ny_global)
      endif
      deallocate(LARRAY_G)

#ifndef NO_MPI
!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else

      allocate(snd_request(nblocks_tot), &
               snd_status (MPI_STATUS_SIZE, nblocks_tot))

      nsends = 0
      do n=1,nblocks_tot
         if (src_dist%blockLocation(n) == my_task+1) then
            nsends = nsends + 1
            src_block = src_dist%blockLocalID(n)
            call MPI_ISEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                 mpiR4, dst_task, 3*mpitag_gs+n, &
                 MPI_COMM_ICE, snd_request(nsends), ierr)
         endif
      end do

      if (nsends > 0) call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
      deallocate(snd_request, snd_status)
#endif

   endif

#ifndef NO_MPI
   if (add_mpi_barriers) then
      call ice_barrier()
   endif
#endif

!-----------------------------------------------------------------------

 end subroutine gather_globalR4

!***********************************************************************

 subroutine scatter_globalR4(ARRAY, ARRAY_G, src_task, dst_dist, &
                             field_loc, field_type, fillValue, &
                             grid_ext)

!  Generated by ./generate_GatherScatter.sh on 2026-08-27
!  This subroutine scatters a global sized array from src_task to all
!  processors.
!  If grid_ext is true, then ARRAY_G is extended grid and updates outer halo
!  If grid_ext is false, then ARRAY_G is 1:nx_global,1:ny_global and outer
!    halo is updated for cyclic and tripole bcs.
!  If grid_ext is false and interior is true, then only the interior data will 
!    be scattered, no outer halo bcs will be applied

   integer (int_kind), intent(in) :: &
      src_task        ! task from which array should be scattered

   type (distrb), intent(in) :: &
      dst_dist        ! distribution of resulting blocks

   real (real_kind), dimension(:,:), intent(in) :: &
      ARRAY_G         ! array containing global field on src_task

   integer (int_kind), intent(in), optional :: &
      field_type,    &! id for type of field (scalar, vector, angle)
      field_loc       ! id for location on horizontal grid
                      !  (center, NEcorner, Nface, Eface)

   real (real_kind), intent(in), optional :: &
      fillValue       ! special value

   logical (log_kind), intent(in), optional :: &
      grid_ext        ! extended grid flag

   real (real_kind), dimension(:,:,:), intent(inout) :: &
      ARRAY           ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,n,         &! dummy loop indices
      ilo,ihi,jlo,jhi,&! block indices
      isrc, jsrc,    &! source addresses
      ioffset, joffset,   &! offsets for tripole boundary conditions
      isign,         &! sign factor for tripole boundary conditions
      dst_block       ! location of block in dst array

   real (dbl_kind) :: &
      efac            ! linear extrapolation factor

   type (block) :: &
      this_block      ! block info for current block

   real (real_kind) :: &
      lfillValue      ! local special value

   logical (log_kind) :: &
      lext,          &! local extended global grid flag
      noupdate        ! flag if noupdate is on

   real (real_kind), dimension(:,:), allocatable :: &
      LARRAY_G    ! array containing global horizontal field on dst_task
                  ! but with indexing that matches nx_global, ny_global
                  ! to simplify overlapping ext and non-ext coding.
                  ! for non-ext, ARRAY_G and LARRAY_G are (nx_global,ny_global).
                  ! for ext incoming ARRAY_G is (1:nx_global+2*nghost,1:ny_global+2*nghost)
                  ! LARRAY_G is same size, but indexing is
                  !            (1-nghost:nx_global+nghost,1-nghost:ny_global+nghost)
                  ! Could maybe use a pointer but a separate array and copy is more robust

#ifndef NO_MPI
   integer (int_kind) :: &
      nrecvs,        &! actual number of messages received
      ierr            ! MPI error flag

   integer (int_kind), dimension(:), allocatable :: &
      rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
      rcv_status      ! status array for receives

   real (real_kind), dimension(:,:), allocatable :: &
      msg_buffer      ! buffer for sending blocks
#endif

   character(len=*), parameter :: subname = '(scatter_globalR4)'

   if (present(fillValue)) then
      lfillValue = fillValue
   else
      lfillValue = 0._real_kind
   endif

   if (present(grid_ext)) then
      lext = grid_ext
   else
      lext = .false.
   endif

   noupdate = .false.
   if (present(field_loc)) then
      if (field_loc == field_loc_noupdate) then
         noupdate = .true.
      endif
   endif
   if (present(field_type)) then
      if (field_type == field_type_noupdate) then
         noupdate = .true.
      endif
   endif

   ! check size of global grid
   if (my_task == src_task) then
      if (lext) then
         if (size(ARRAY_G,dim=1)/= nx_global+2*nghost .or. size(ARRAY_G,dim=2) /= ny_global+2*nghost) then
            call abort_ice(subname//'ERROR: grid_ext = true, wrong ARRAY_G size')
         endif
      else
         if (size(ARRAY_G,dim=1)/= nx_global .or. size(ARRAY_G,dim=2) /= ny_global) then
            call abort_ice(subname//'ERROR: grid_ext = false, wrong ARRAY_G size')
         endif
      endif
   endif

   ARRAY = lfillValue

   isign = -99
   ioffset = -999
   joffset = -999
   if (.not.lext .and. .not.noupdate .and. &
      (ns_boundary_type == 'tripole' .or. ns_boundary_type == 'tripoleT')) then

      select case (field_type)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case (field_type_noupdate) ! ghost cells never used - use cell center
         isign =  1
      case default
         call abort_ice(subname//'ERROR: Unknown field kind')
      end select

      if (ns_boundary_type == 'tripole') then
         ioffset = 0
         joffset = 0
      else ! tripoleT fold
         ioffset = -1
         joffset = 1
      endif

      select case (field_loc)
        case (field_loc_center)     ! cell center location
           !ioffset = ioffset
           !joffset = joffset
        case (field_loc_NEcorner)   ! cell corner (velocity) location
           ioffset = ioffset + 1
           joffset = joffset + 1
        case (field_loc_Eface)      ! cell east face
           ioffset = ioffset + 1
           !joffset = joffset
        case (field_loc_Nface)      ! cell north face
           !ioffset = ioffset
           joffset = joffset + 1
        case (field_loc_noupdate)   ! no update
           !ioffset = ioffset
           !joffset = joffset
        case default
           call abort_ice(subname//'ERROR: Unknown field location')
      end select

   endif ! not lext

!-----------------------------------------------------------------------
!
!  if this task is the src_task, compute updated LARRAY_G which contains
!  outer halo extended grid from ARRAY_G (extended or not), fill outer
!  halo in global array then scatter
!
!  i_global, j_global are indexed 1-nghost:n*_global+nghost
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

      allocate(LARRAY_G(1-nghost:nx_global+nghost,1-nghost:ny_global+nghost))

      LARRAY_G = lfillValue
      if (lext) then
         do j = 1-nghost, ny_global+nghost
         do i = 1-nghost, nx_global+nghost
            LARRAY_G(i,j) = ARRAY_G(i+nghost,j+nghost)
         enddo
         enddo
      else

         ! fill cells that can be directly copied

         do j = 1-nghost, ny_global+nghost
         do i = 1-nghost, nx_global+nghost
            isrc = i_global(i)
            jsrc = j_global(j)
            if (isrc >= 1 .and. isrc <= nx_global .and. &
                jsrc >= 1 .and. jsrc <= ny_global) then
               ! copy valid cell directly
               LARRAY_G(i,j) = ARRAY_G(isrc,jsrc)
            endif
         enddo
         enddo

         ! fill outer halo, only if noupdate is false
         ! do north/south boundary first (partly due to tripole feature) (1:nx_global)
         ! then do east/west including corners (1-nghost:ny_global+nghost)

         if (.not.noupdate) then
            if (ns_boundary_type == 'tripole' .or. ns_boundary_type == 'tripoleT') then
               do j = ny_global, ny_global+nghost
               ! north of active cells
               do i = 1, nx_global
                  isrc = nx_global - i + 1 - ioffset   ! ioffset = 0 for u center, ioffset = 1 for T center
                  jsrc = ny_global - (j-ny_global) - joffset + 1  ! joffset = 0 for u center, joffset = 1 for T center
                  if (isrc < 1        ) isrc = isrc + nx_global
                  if (isrc > nx_global) isrc = isrc - nx_global

                  !*** for center and Eface on u-fold, and NE corner and Nface
                  !*** on T-fold, do not need to replace
                  !*** top row of physical domain, so jsrc should be greater than j

                  if (jsrc > j) then
                     ! do nothing
                  elseif (jsrc == j) then
                     ! don't average in the scatter even though halo update does
                     ! LARRAY_G(i,j) = 0.5_dbl_kind * (ARRAY_G(i,jsrc) + isign*ARRAY_G(isrc,jsrc))
                  else
                     ! copy
                     LARRAY_G(i,j) = isign * ARRAY_G(isrc,jsrc)
                  endif
               enddo
               enddo
            else
               !--- bottom edge
               do j = 1-nghost, 0
               do i = 1, nx_global
                  if (ns_boundary_type == 'cyclic') then
                     LARRAY_G(i,j) = ARRAY_G(i,j+ny_global)
                  elseif (ns_boundary_type == 'zero_gradient') then
                     LARRAY_G(i,j) = ARRAY_G(i,1)
                  elseif (ns_boundary_type == 'linear_extrap') then
                     efac = real(2-j,dbl_kind)
                     LARRAY_G(i,j) = efac*ARRAY_G(i,1) - (efac-c1)*ARRAY_G(i,2)
                  endif
               enddo
               enddo
               !--- top edge
               do j = ny_global+1, ny_global+nghost
               do i = 1, nx_global
                  if (ns_boundary_type == 'cyclic') then
                     LARRAY_G(i,j) = ARRAY_G(i,j-ny_global)
                  elseif (ns_boundary_type == 'zero_gradient') then
                     LARRAY_G(i,j) = ARRAY_G(i,ny_global)
                  elseif (ns_boundary_type == 'linear_extrap') then
                     efac = real(j-ny_global+1,dbl_kind)
                     LARRAY_G(i,j) = efac*ARRAY_G(i,ny_global) - (efac-c1)*ARRAY_G(i,ny_global-1)
                  endif
               enddo
               enddo
            endif

            do j = 1-nghost, ny_global+nghost
            !--- left edge
            do i = 1-nghost, 0
               if (ew_boundary_type == 'cyclic') then
                  LARRAY_G(i,j) = LARRAY_G(i+nx_global,j)
               elseif (ew_boundary_type == 'zero_gradient') then
                  LARRAY_G(i,j) = LARRAY_G(1,j)
               elseif (ew_boundary_type == 'linear_extrap') then
                  efac = real(2-i,dbl_kind)
                  LARRAY_G(i,j) = efac*LARRAY_G(1,j) - (efac-c1)*LARRAY_G(2,j)
               endif
            enddo
            !--- right edge
            do i = nx_global+1, nx_global+nghost
               if (ew_boundary_type == 'cyclic') then
                  LARRAY_G(i,j) = LARRAY_G(i-nx_global,j)
               elseif (ew_boundary_type == 'zero_gradient') then
                  LARRAY_G(i,j) = LARRAY_G(nx_global,j)
               elseif (ew_boundary_type == 'linear_extrap') then
                  efac = real(i-nx_global+1,dbl_kind)
                  LARRAY_G(i,j) = efac*LARRAY_G(nx_global,j) - (efac-c1)*LARRAY_G(nx_global-1,j)
               endif
            enddo
            enddo

         endif   ! noupdate
      endif   ! lext

      !*** send non-local blocks away

#ifndef NO_MPI
      allocate (msg_buffer(nx_block,ny_block))

      do n=1,nblocks_tot
         if (dst_dist%blockLocation(n) > 0 .and. &
             dst_dist%blockLocation(n)-1 /= my_task) then

            msg_buffer = lfillValue
            this_block = get_block(n,n)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j=jlo-nghost,jhi+nghost
            do i=ilo-nghost,ihi+nghost
               ! fill entire block
               isrc = this_block%i_glob(ilo) + i - ilo
               jsrc = this_block%j_glob(jlo) + j - jlo
               msg_buffer(i,j) = LARRAY_G(isrc,jsrc)
            enddo
            enddo

            call MPI_SEND(msg_buffer, nx_block*ny_block, &
                          mpiR4, dst_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                          MPI_COMM_ICE, ierr)

         endif
      end do

      deallocate(msg_buffer)
#endif
      !*** copy any local blocks

      do n=1,nblocks_tot
         if (dst_dist%blockLocation(n) == my_task+1) then
            dst_block = dst_dist%blockLocalID(n)
            this_block = get_block(n,n)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j=jlo-nghost,jhi+nghost
            do i=ilo-nghost,ihi+nghost
               ! fill entire block
               isrc = this_block%i_glob(ilo) + i - ilo
               jsrc = this_block%j_glob(jlo) + j - jlo
               ARRAY(i,j,dst_block) = LARRAY_G(isrc,jsrc)
            enddo
            enddo
         endif
      end do

#ifndef NO_MPI
!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else   ! my_task /= src_task

      allocate (rcv_request(nblocks_tot), &
                rcv_status(MPI_STATUS_SIZE, nblocks_tot))

      rcv_request = 0
      rcv_status  = 0

      nrecvs = 0
      do n=1,nblocks_tot
         if (dst_dist%blockLocation(n) == my_task+1) then
            nrecvs = nrecvs + 1
            dst_block = dst_dist%blockLocalID(n)
            call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                          mpiR4, src_task, 3*mpitag_gs+n, &
                          MPI_COMM_ICE, rcv_request(nrecvs), ierr)
         endif
      end do

      if (nrecvs > 0) &
         call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

      deallocate(rcv_request, rcv_status)
#endif
   endif   ! my_task == src_task

#ifndef NO_MPI
   if (add_mpi_barriers) then
     call ice_barrier()
   endif
#endif

!-----------------------------------------------------------------------

 end subroutine scatter_globalR4

!***********************************************************************

 subroutine gather_globalI4(ARRAY_G, ARRAY, dst_task, src_dist, fillValue, grid_ext)

!  Generated by ./generate_GatherScatter.sh on 2026-08-27
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!  If grid_ext is true, then ARRAY_G is extended grid
!  If grid_ext is false, then ARRAY_G is 1:nx_global,1:ny_global

   integer (int_kind), intent(in) :: &
      dst_task   ! task to which array should be gathered

   type (distrb), intent(in) :: &
      src_dist   ! distribution of blocks in the source array

   integer (int_kind), dimension(:,:,:), intent(in) :: &
      ARRAY      ! array containing horizontal slab of distributed field

   integer (int_kind), intent(in), optional :: &
      fillValue  ! special value

   logical (log_kind), intent(in), optional :: &
      grid_ext   ! extended global grid flag

   integer (int_kind), dimension(:,:), intent(inout) :: &
      ARRAY_G    ! array containing global horizontal field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,n,          &! dummy loop counters
      iblock,jblock,  &! block index
      ilo,ihi,jlo,jhi  ! block sizes

#ifndef NO_MPI
   integer (int_kind) :: &
      nsends,         &! number of actual sends
      src_block,      &! block locator for send
      ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
      status

   integer (int_kind), dimension(:), allocatable :: &
      snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
      snd_status

   integer (int_kind), dimension(:,:), allocatable :: &
      msg_buffer
#endif

   integer (int_kind) :: &
      lfillValue  ! local special value

   logical (log_kind) :: &
      lext        ! local grid_ext global grid flag

   type (block) :: &
      this_block  ! block info for current block

   integer (int_kind), dimension(:,:), allocatable :: &
      LARRAY_G    ! array containing global horizontal field on dst_task
                  ! but with indexing that matches nx_global, ny_global
                  ! to simplify overlapping ext and non-ext coding.
                  ! for non-ext, ARRAY_G and LARRAY_G are (nx_global,ny_global).
                  ! for ext incoming ARRAY_G is (1:nx_global+2*nghost,1:ny_global+2*nghost)
                  ! LARRAY_G is same size, but indexing is
                  !            (1-nghost:nx_global+nghost,1-nghost:ny_global+nghost)
                  ! Could maybe use a pointer but a separate array and copy is more robust

   character(len=*), parameter :: subname = '(gather_globalI4)'

   if (present(fillValue)) then
      lfillValue = fillValue
   else
      lfillValue = spval_int
   endif

   if (present(grid_ext)) then
      lext = grid_ext
   else
      lext = .false.
   endif

!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (my_task == dst_task) then

      if (lext) then
         allocate(LARRAY_G(1-nghost:nx_global+nghost,1-nghost:ny_global+nghost))
         LARRAY_G(1-nghost:nx_global+nghost,1-nghost:ny_global+nghost) = ARRAY_G(1:nx_global+2*nghost,1:ny_global+2*nghost)
      else
         allocate(LARRAY_G(1:nx_global,1:ny_global))
         LARRAY_G(:,:) = ARRAY_G(:,:)
      endif

      do n=1,nblocks_tot

         this_block = get_block(n,n)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         iblock = this_block%iblock
         jblock = this_block%jblock

         !*** copy local blocks

         if (src_dist%blockLocation(n) == my_task+1) then

            do j = jlo, jhi
            do i = ilo, ihi
               LARRAY_G(this_block%i_glob(i),this_block%j_glob(j)) = ARRAY(i,j,src_dist%blockLocalID(n))
            end do
            end do

            !*** fill outside halo for extended grid gather

            if (lext) then
               ! south
               if (jblock == 1) then
                  do j = 1, nghost
                  do i = ilo, ihi
                     LARRAY_G(this_block%i_glob(i),j-nghost) = ARRAY(i,j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! north
               if (jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = ilo, ihi
                     LARRAY_G(this_block%i_glob(i),ny_global+j) = ARRAY(i,jhi+j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! west
               if (iblock == 1) then
                  do j = jlo, jhi
                  do i = 1, nghost
                     LARRAY_G(i-nghost,this_block%j_glob(j)) = ARRAY(i,j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! east
               if (iblock == nblocks_x) then
                  do j = jlo, jhi
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,this_block%j_glob(j)) = ARRAY(ihi+i,j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! southwest
               if (iblock == 1 .and. jblock == 1) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(i-nghost,j-nghost) = ARRAY(i,j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! northwest
               if (iblock == 1 .and. jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(i-nghost,ny_global+j) = ARRAY(i,jhi+j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! southeast
               if (iblock == nblocks_x .and. jblock == 1) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,j-nghost) = ARRAY(ihi+i,j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif

               ! northeast
               if (iblock == nblocks_x .and. jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,ny_global+j) = ARRAY(ihi+i,jhi+j,src_dist%blockLocalID(n))
                  enddo
                  enddo
               endif
            endif

         !*** fill land blocks with special values

         else if (src_dist%blockLocation(n) == 0) then

            do j = jlo, jhi
            do i = ilo, ihi
               LARRAY_G(this_block%i_glob(i),this_block%j_glob(j)) = lfillValue
            end do
            end do

         endif

      end do

#ifndef NO_MPI
      !*** copy blocks from other tasks

      allocate (msg_buffer(nx_block,ny_block))

      do n=1,nblocks_tot
         if (src_dist%blockLocation(n) > 0 .and. &
             src_dist%blockLocation(n) /= my_task+1) then

            this_block = get_block(n,n)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            iblock = this_block%iblock
            jblock = this_block%jblock

            call MPI_RECV(msg_buffer, size(msg_buffer), &
                 MPI_INTEGER, src_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                 MPI_COMM_ICE, status, ierr)

            do j = jlo, jhi
            do i = ilo, ihi
               LARRAY_G(this_block%i_glob(i),this_block%j_glob(j)) = msg_buffer(i,j)
            end do
            end do

            !*** fill outside halo for extended grid gather

            if (lext) then
               ! south
               if (jblock == 1) then
                  do j = 1, nghost
                  do i = ilo, ihi
                     LARRAY_G(this_block%i_glob(i),j-nghost) = msg_buffer(i,j)
                  enddo
                  enddo
               endif

               ! north
               if (jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = ilo, ihi
                     LARRAY_G(this_block%i_glob(i),ny_global+j) = msg_buffer(i,jhi+j)
                  enddo
                  enddo
               endif

               ! west
               if (iblock == 1) then
                  do j = jlo, jhi
                  do i = 1, nghost
                     LARRAY_G(i-nghost,this_block%j_glob(j)) = msg_buffer(i,j)
                  enddo
                  enddo
               endif

               ! east
               if (iblock == nblocks_x) then
                  do j = jlo, jhi
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,this_block%j_glob(j)) = msg_buffer(ihi+i,j)
                  enddo
                  enddo
               endif

               ! southwest
               if (iblock == 1 .and. jblock == 1) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(i-nghost,j-nghost) = msg_buffer(i,j)
                  enddo
                  enddo
               endif

               ! northwest
               if (iblock == 1 .and. jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(i-nghost,ny_global+j) = msg_buffer(i,jhi+j)
                  enddo
                  enddo
               endif

               ! southeast
               if (iblock == nblocks_x .and. jblock == 1) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,j-nghost) = msg_buffer(ihi+i,j)
                  enddo
                  enddo
               endif

               ! northeast
               if (iblock == nblocks_x .and. jblock == nblocks_y) then
                  do j = 1, nghost
                  do i = 1, nghost
                     LARRAY_G(nx_global+i,ny_global+j) = msg_buffer(ihi+i,jhi+j)
                  enddo
                  enddo
               endif
            endif ! lext
         endif
      end do

      deallocate(msg_buffer)
#endif

      ! copy LARRAY_G to ARRAY_G, offset extended grid index

      if (lext) then
         ARRAY_G(1:nx_global+2*nghost,1:ny_global+2*nghost) = LARRAY_G(1-nghost:nx_global+nghost,1-nghost:ny_global+nghost)
      else
         ARRAY_G(1:nx_global,1:ny_global) = LARRAY_G(1:nx_global,1:ny_global)
      endif
      deallocate(LARRAY_G)

#ifndef NO_MPI
!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else

      allocate(snd_request(nblocks_tot), &
               snd_status (MPI_STATUS_SIZE, nblocks_tot))

      nsends = 0
      do n=1,nblocks_tot
         if (src_dist%blockLocation(n) == my_task+1) then
            nsends = nsends + 1
            src_block = src_dist%blockLocalID(n)
            call MPI_ISEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                 MPI_INTEGER, dst_task, 3*mpitag_gs+n, &
                 MPI_COMM_ICE, snd_request(nsends), ierr)
         endif
      end do

      if (nsends > 0) call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
      deallocate(snd_request, snd_status)
#endif

   endif

#ifndef NO_MPI
   if (add_mpi_barriers) then
      call ice_barrier()
   endif
#endif

!-----------------------------------------------------------------------

 end subroutine gather_globalI4

!***********************************************************************

 subroutine scatter_globalI4(ARRAY, ARRAY_G, src_task, dst_dist, &
                             field_loc, field_type, fillValue, &
                             grid_ext)

!  Generated by ./generate_GatherScatter.sh on 2026-08-27
!  This subroutine scatters a global sized array from src_task to all
!  processors.
!  If grid_ext is true, then ARRAY_G is extended grid and updates outer halo
!  If grid_ext is false, then ARRAY_G is 1:nx_global,1:ny_global and outer
!    halo is updated for cyclic and tripole bcs.
!  If grid_ext is false and interior is true, then only the interior data will 
!    be scattered, no outer halo bcs will be applied

   integer (int_kind), intent(in) :: &
      src_task        ! task from which array should be scattered

   type (distrb), intent(in) :: &
      dst_dist        ! distribution of resulting blocks

   integer (int_kind), dimension(:,:), intent(in) :: &
      ARRAY_G         ! array containing global field on src_task

   integer (int_kind), intent(in), optional :: &
      field_type,    &! id for type of field (scalar, vector, angle)
      field_loc       ! id for location on horizontal grid
                      !  (center, NEcorner, Nface, Eface)

   integer (int_kind), intent(in), optional :: &
      fillValue       ! special value

   logical (log_kind), intent(in), optional :: &
      grid_ext        ! extended grid flag

   integer (int_kind), dimension(:,:,:), intent(inout) :: &
      ARRAY           ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,n,         &! dummy loop indices
      ilo,ihi,jlo,jhi,&! block indices
      isrc, jsrc,    &! source addresses
      ioffset, joffset,   &! offsets for tripole boundary conditions
      isign,         &! sign factor for tripole boundary conditions
      dst_block       ! location of block in dst array

   real (dbl_kind) :: &
      efac            ! linear extrapolation factor

   type (block) :: &
      this_block      ! block info for current block

   integer (int_kind) :: &
      lfillValue      ! local special value

   logical (log_kind) :: &
      lext,          &! local extended global grid flag
      noupdate        ! flag if noupdate is on

   integer (int_kind), dimension(:,:), allocatable :: &
      LARRAY_G    ! array containing global horizontal field on dst_task
                  ! but with indexing that matches nx_global, ny_global
                  ! to simplify overlapping ext and non-ext coding.
                  ! for non-ext, ARRAY_G and LARRAY_G are (nx_global,ny_global).
                  ! for ext incoming ARRAY_G is (1:nx_global+2*nghost,1:ny_global+2*nghost)
                  ! LARRAY_G is same size, but indexing is
                  !            (1-nghost:nx_global+nghost,1-nghost:ny_global+nghost)
                  ! Could maybe use a pointer but a separate array and copy is more robust

#ifndef NO_MPI
   integer (int_kind) :: &
      nrecvs,        &! actual number of messages received
      ierr            ! MPI error flag

   integer (int_kind), dimension(:), allocatable :: &
      rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
      rcv_status      ! status array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
      msg_buffer      ! buffer for sending blocks
#endif

   character(len=*), parameter :: subname = '(scatter_globalI4)'

   if (present(fillValue)) then
      lfillValue = fillValue
   else
      lfillValue = 0
   endif

   if (present(grid_ext)) then
      lext = grid_ext
   else
      lext = .false.
   endif

   noupdate = .false.
   if (present(field_loc)) then
      if (field_loc == field_loc_noupdate) then
         noupdate = .true.
      endif
   endif
   if (present(field_type)) then
      if (field_type == field_type_noupdate) then
         noupdate = .true.
      endif
   endif

   ! check size of global grid
   if (my_task == src_task) then
      if (lext) then
         if (size(ARRAY_G,dim=1)/= nx_global+2*nghost .or. size(ARRAY_G,dim=2) /= ny_global+2*nghost) then
            call abort_ice(subname//'ERROR: grid_ext = true, wrong ARRAY_G size')
         endif
      else
         if (size(ARRAY_G,dim=1)/= nx_global .or. size(ARRAY_G,dim=2) /= ny_global) then
            call abort_ice(subname//'ERROR: grid_ext = false, wrong ARRAY_G size')
         endif
      endif
   endif

   ARRAY = lfillValue

   isign = -99
   ioffset = -999
   joffset = -999
   if (.not.lext .and. .not.noupdate .and. &
      (ns_boundary_type == 'tripole' .or. ns_boundary_type == 'tripoleT')) then

      select case (field_type)
      case (field_type_scalar)
         isign =  1
      case (field_type_vector)
         isign = -1
      case (field_type_angle)
         isign = -1
      case (field_type_noupdate) ! ghost cells never used - use cell center
         isign =  1
      case default
         call abort_ice(subname//'ERROR: Unknown field kind')
      end select

      if (ns_boundary_type == 'tripole') then
         ioffset = 0
         joffset = 0
      else ! tripoleT fold
         ioffset = -1
         joffset = 1
      endif

      select case (field_loc)
        case (field_loc_center)     ! cell center location
           !ioffset = ioffset
           !joffset = joffset
        case (field_loc_NEcorner)   ! cell corner (velocity) location
           ioffset = ioffset + 1
           joffset = joffset + 1
        case (field_loc_Eface)      ! cell east face
           ioffset = ioffset + 1
           !joffset = joffset
        case (field_loc_Nface)      ! cell north face
           !ioffset = ioffset
           joffset = joffset + 1
        case (field_loc_noupdate)   ! no update
           !ioffset = ioffset
           !joffset = joffset
        case default
           call abort_ice(subname//'ERROR: Unknown field location')
      end select

   endif ! not lext

!-----------------------------------------------------------------------
!
!  if this task is the src_task, compute updated LARRAY_G which contains
!  outer halo extended grid from ARRAY_G (extended or not), fill outer
!  halo in global array then scatter
!
!  i_global, j_global are indexed 1-nghost:n*_global+nghost
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

      allocate(LARRAY_G(1-nghost:nx_global+nghost,1-nghost:ny_global+nghost))

      LARRAY_G = lfillValue
      if (lext) then
         do j = 1-nghost, ny_global+nghost
         do i = 1-nghost, nx_global+nghost
            LARRAY_G(i,j) = ARRAY_G(i+nghost,j+nghost)
         enddo
         enddo
      else

         ! fill cells that can be directly copied

         do j = 1-nghost, ny_global+nghost
         do i = 1-nghost, nx_global+nghost
            isrc = i_global(i)
            jsrc = j_global(j)
            if (isrc >= 1 .and. isrc <= nx_global .and. &
                jsrc >= 1 .and. jsrc <= ny_global) then
               ! copy valid cell directly
               LARRAY_G(i,j) = ARRAY_G(isrc,jsrc)
            endif
         enddo
         enddo

         ! fill outer halo, only if noupdate is false
         ! do north/south boundary first (partly due to tripole feature) (1:nx_global)
         ! then do east/west including corners (1-nghost:ny_global+nghost)

         if (.not.noupdate) then
            if (ns_boundary_type == 'tripole' .or. ns_boundary_type == 'tripoleT') then
               do j = ny_global, ny_global+nghost
               ! north of active cells
               do i = 1, nx_global
                  isrc = nx_global - i + 1 - ioffset   ! ioffset = 0 for u center, ioffset = 1 for T center
                  jsrc = ny_global - (j-ny_global) - joffset + 1  ! joffset = 0 for u center, joffset = 1 for T center
                  if (isrc < 1        ) isrc = isrc + nx_global
                  if (isrc > nx_global) isrc = isrc - nx_global

                  !*** for center and Eface on u-fold, and NE corner and Nface
                  !*** on T-fold, do not need to replace
                  !*** top row of physical domain, so jsrc should be greater than j

                  if (jsrc > j) then
                     ! do nothing
                  elseif (jsrc == j) then
                     ! don't average in the scatter even though halo update does
                     ! LARRAY_G(i,j) = 0.5_dbl_kind * (ARRAY_G(i,jsrc) + isign*ARRAY_G(isrc,jsrc))
                  else
                     ! copy
                     LARRAY_G(i,j) = isign * ARRAY_G(isrc,jsrc)
                  endif
               enddo
               enddo
            else
               !--- bottom edge
               do j = 1-nghost, 0
               do i = 1, nx_global
                  if (ns_boundary_type == 'cyclic') then
                     LARRAY_G(i,j) = ARRAY_G(i,j+ny_global)
                  elseif (ns_boundary_type == 'zero_gradient') then
                     LARRAY_G(i,j) = ARRAY_G(i,1)
                  elseif (ns_boundary_type == 'linear_extrap') then
                     efac = real(2-j,dbl_kind)
                     LARRAY_G(i,j) = efac*ARRAY_G(i,1) - (efac-c1)*ARRAY_G(i,2)
                  endif
               enddo
               enddo
               !--- top edge
               do j = ny_global+1, ny_global+nghost
               do i = 1, nx_global
                  if (ns_boundary_type == 'cyclic') then
                     LARRAY_G(i,j) = ARRAY_G(i,j-ny_global)
                  elseif (ns_boundary_type == 'zero_gradient') then
                     LARRAY_G(i,j) = ARRAY_G(i,ny_global)
                  elseif (ns_boundary_type == 'linear_extrap') then
                     efac = real(j-ny_global+1,dbl_kind)
                     LARRAY_G(i,j) = efac*ARRAY_G(i,ny_global) - (efac-c1)*ARRAY_G(i,ny_global-1)
                  endif
               enddo
               enddo
            endif

            do j = 1-nghost, ny_global+nghost
            !--- left edge
            do i = 1-nghost, 0
               if (ew_boundary_type == 'cyclic') then
                  LARRAY_G(i,j) = LARRAY_G(i+nx_global,j)
               elseif (ew_boundary_type == 'zero_gradient') then
                  LARRAY_G(i,j) = LARRAY_G(1,j)
               elseif (ew_boundary_type == 'linear_extrap') then
                  efac = real(2-i,dbl_kind)
                  LARRAY_G(i,j) = efac*LARRAY_G(1,j) - (efac-c1)*LARRAY_G(2,j)
               endif
            enddo
            !--- right edge
            do i = nx_global+1, nx_global+nghost
               if (ew_boundary_type == 'cyclic') then
                  LARRAY_G(i,j) = LARRAY_G(i-nx_global,j)
               elseif (ew_boundary_type == 'zero_gradient') then
                  LARRAY_G(i,j) = LARRAY_G(nx_global,j)
               elseif (ew_boundary_type == 'linear_extrap') then
                  efac = real(i-nx_global+1,dbl_kind)
                  LARRAY_G(i,j) = efac*LARRAY_G(nx_global,j) - (efac-c1)*LARRAY_G(nx_global-1,j)
               endif
            enddo
            enddo

         endif   ! noupdate
      endif   ! lext

      !*** send non-local blocks away

#ifndef NO_MPI
      allocate (msg_buffer(nx_block,ny_block))

      do n=1,nblocks_tot
         if (dst_dist%blockLocation(n) > 0 .and. &
             dst_dist%blockLocation(n)-1 /= my_task) then

            msg_buffer = lfillValue
            this_block = get_block(n,n)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j=jlo-nghost,jhi+nghost
            do i=ilo-nghost,ihi+nghost
               ! fill entire block
               isrc = this_block%i_glob(ilo) + i - ilo
               jsrc = this_block%j_glob(jlo) + j - jlo
               msg_buffer(i,j) = LARRAY_G(isrc,jsrc)
            enddo
            enddo

            call MPI_SEND(msg_buffer, nx_block*ny_block, &
                          MPI_INTEGER, dst_dist%blockLocation(n)-1, 3*mpitag_gs+n, &
                          MPI_COMM_ICE, ierr)

         endif
      end do

      deallocate(msg_buffer)
#endif
      !*** copy any local blocks

      do n=1,nblocks_tot
         if (dst_dist%blockLocation(n) == my_task+1) then
            dst_block = dst_dist%blockLocalID(n)
            this_block = get_block(n,n)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j=jlo-nghost,jhi+nghost
            do i=ilo-nghost,ihi+nghost
               ! fill entire block
               isrc = this_block%i_glob(ilo) + i - ilo
               jsrc = this_block%j_glob(jlo) + j - jlo
               ARRAY(i,j,dst_block) = LARRAY_G(isrc,jsrc)
            enddo
            enddo
         endif
      end do

#ifndef NO_MPI
!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else   ! my_task /= src_task

      allocate (rcv_request(nblocks_tot), &
                rcv_status(MPI_STATUS_SIZE, nblocks_tot))

      rcv_request = 0
      rcv_status  = 0

      nrecvs = 0
      do n=1,nblocks_tot
         if (dst_dist%blockLocation(n) == my_task+1) then
            nrecvs = nrecvs + 1
            dst_block = dst_dist%blockLocalID(n)
            call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                          MPI_INTEGER, src_task, 3*mpitag_gs+n, &
                          MPI_COMM_ICE, rcv_request(nrecvs), ierr)
         endif
      end do

      if (nrecvs > 0) &
         call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

      deallocate(rcv_request, rcv_status)
#endif
   endif   ! my_task == src_task

#ifndef NO_MPI
   if (add_mpi_barriers) then
     call ice_barrier()
   endif
#endif

!-----------------------------------------------------------------------

 end subroutine scatter_globalI4

!***********************************************************************
! autogenerate end
!***********************************************************************

 subroutine gather_globalL1(ARRAY_G, ARRAY, dst_task, src_dist, fillValue, grid_ext)

!  Convert to integer and call gather_globalI4
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!  If grid_ext is true, then ARRAY_G is extended grid
!  If grid_ext is false, then ARRAY_G is 1:nx_global,1:ny_global

   integer (int_kind), intent(in) :: &
      dst_task   ! task to which array should be gathered

   type (distrb), intent(in) :: &
      src_dist   ! distribution of blocks in the source array

   logical (log_kind), dimension(:,:,:), intent(in) :: &
      ARRAY      ! array containing horizontal slab of distributed field

   logical (log_kind), intent(in), optional :: &
      fillValue  ! special value

   logical (log_kind), intent(in), optional :: &
      grid_ext   ! extended global grid flag

   logical (log_kind), dimension(:,:), intent(inout) :: &
      ARRAY_G    ! array containing global horizontal field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:,:), allocatable :: &
      IARRAY     ! integer array containing horizontal slab of distributed field

   integer (int_kind), dimension(:,:), allocatable :: &
      IARRAY_G   ! integer array containing global horizontal field on dst_task

   integer (int_kind) :: &
      ifillValue ! integer fill value

   character(len=*), parameter :: subname = '(gather_globalL1)'

   allocate(IARRAY(size(ARRAY,dim=1),size(ARRAY,dim=2),size(ARRAY,dim=3)))
   allocate(IARRAY_G(size(ARRAY_G,dim=1),size(ARRAY_G,dim=2)))

   IARRAY = 0   ! 0 = false
   where (ARRAY) IARRAY = 1   ! 1 = true
   IARRAY_G = 0
   where (ARRAY_G) IARRAY_G = 1

   if (present(fillValue)) then
      ifillValue = 0
      if (fillValue) ifillValue = 1
      call gather_global(IARRAY_G, IARRAY, dst_task, src_dist, ifillValue, grid_ext)
   else
      call gather_global(IARRAY_G, IARRAY, dst_task, src_dist, grid_ext = grid_ext)
   endif

   ARRAY_G = .false.  ! 0 = false
   where (mod(abs(IARRAY_G),2) == 1) ARRAY_G = .true.  ! .ne.0 = true

!-----------------------------------------------------------------------

 end subroutine gather_globalL1

!***********************************************************************
!
!  The two routines below are NOT generated by ./generate_GatherScatter.sh --
!  they are hand-written and must be preserved if that script is re-run.
!
!  They gather a distributed field whose i and j dimensions are flattened,
!  i.e. dimension(nx_block*ny_block,max_blocks), where element (i,j,iblk)
!  lives at ((j-1)*nx_block + i, iblk).  ARRAY is aliased to rank 3 and
!  handed to the corresponding generated routine, so no data is copied.
!
!***********************************************************************

 subroutine gather_globalR8_flat(ARRAY_G, ARRAY, dst_task, src_dist, fillValue, grid_ext)

   integer (int_kind), intent(in) :: &
      dst_task   ! task to which array should be gathered

   type (distrb), intent(in) :: &
      src_dist   ! distribution of blocks in the source array

   ! TARGET and CONTIGUOUS are required to alias this to rank 3 below
   real (dbl_kind), dimension(:,:), intent(in), target, contiguous :: &
      ARRAY      ! array containing horizontal slab of distributed field

   real (dbl_kind), intent(in), optional :: &
      fillValue  ! special value

   logical (log_kind), intent(in), optional :: &
      grid_ext   ! extended global grid flag

   real (dbl_kind), dimension(:,:), intent(inout) :: &
      ARRAY_G    ! array containing global horizontal field on dst_task

   real (dbl_kind), pointer :: &
      ARRAY_3d(:,:,:)   ! rank-3 alias of ARRAY

   character(len=*), parameter :: subname = '(gather_globalR8_flat)'

   if (size(ARRAY,dim=1) /= nx_block*ny_block) then
      call abort_ice(subname//' ERROR: first dimension is not nx_block*ny_block')
   endif

   ARRAY_3d(1:nx_block,1:ny_block,1:size(ARRAY,dim=2)) => ARRAY

   call gather_globalR8(ARRAY_G, ARRAY_3d, dst_task, src_dist, fillValue, grid_ext)

 end subroutine gather_globalR8_flat

!***********************************************************************

 subroutine gather_globalR8_compressed(ARRAY_G, ARRAY, map, off,       &
                                       dst_task, src_dist, spc_val,    &
                                       grid_ext)

!  Gather a compressed (slot) field into a global block-shaped array.  The
!  mirror of scatter_globalR8_compressed, and needed for the same round trip:
!  evp1d gathers the diagnostics in so that cells it does not own are written
!  back unchanged rather than zeroed.
!
!  Cells with no slot are left at whatever spc_val/the gather puts there --
!  they hold nothing the compressed producer wrote.

   integer (int_kind), intent(in) :: &
      dst_task        ! task to which array should be gathered

   type (distrb), intent(in) :: &
      src_dist        ! distribution of blocks in the source array

   real (dbl_kind), dimension(:), intent(in) :: &
      ARRAY           ! compressed field, dimension(dynTot)

   integer (int_kind), dimension(:,:), intent(in) :: &
      map             ! block-local slot of each cell, 0 if not stored

   integer (int_kind), dimension(:), intent(in) :: &
      off             ! first slot of each block, minus one

   real (dbl_kind), intent(in), optional :: &
      spc_val

   logical (log_kind), intent(in), optional :: &
      grid_ext

   real (dbl_kind), dimension(:,:), intent(inout) :: &
      ARRAY_G         ! array containing global field on dst_task

!-----------------------------------------------------------------------

   real (dbl_kind), dimension(:,:,:), allocatable :: &
      ARRAY_b         ! block-shaped scratch, unpacked into then gathered

   integer (int_kind) :: i, j, ij, iblk, ierr

   character(len=*), parameter :: subname = '(gather_globalR8_compressed)'

!-----------------------------------------------------------------------

   allocate(ARRAY_b(nx_block,ny_block,size(map,dim=2)), stat=ierr)
   if (ierr /= 0) call abort_ice(subname//' ERROR: out of memory')
   ARRAY_b = c0

   do iblk = 1, size(map,dim=2)
      do j = 1, ny_block
      do i = 1, nx_block
         ij = (j-1)*nx_block + i
         if (map(ij,iblk) > 0) then
            ARRAY_b(i,j,iblk) = ARRAY(off(iblk) + map(ij,iblk))
         endif
      enddo
      enddo
   enddo

   call gather_globalR8(ARRAY_G, ARRAY_b, dst_task, src_dist, spc_val, grid_ext)

   deallocate(ARRAY_b)

 end subroutine gather_globalR8_compressed

!***********************************************************************

 subroutine slot_gather_setup(map, off, cnt, indxTij, nActive, &
                              nxext, nyext, dst_task, src_dist)

    use ice_fileunits, only: nu_diag

!  Build the descriptors gather_slots_r8 uses.  Call once, after the ragged
!  layout and the destination 1-D cell list both exist.
!
!  ⚠ The block-ownership rule -- which block supplies which global cell,
!  interior plus the grid_ext rims and their corners -- is MEASURED here,
!  not restated.  Every slot on every task is stamped with a unique
!  (block,slot) id, the id field is pushed through
!  gather_globalR8_compressed, and dst_task reads back which id landed at
!  each of its 1-D cells.  Whatever that routine does, this reproduces,
!  because it is that routine that answered.  Restating the rule instead is
!  how the two paths drift apart silently.
!
!  Both sides then agree to order a block's contribution by ASCENDING SLOT,
!  so no index travels with the data afterwards.

    integer (int_kind), dimension(:,:), intent(in) :: &
       map          ! dynMap: block-local slot of each cell, 0 if none
    integer (int_kind), dimension(:), intent(in) :: &
       off,       & ! dynOff: first slot of each block, minus one
       cnt          ! dynCnt: slots held by each block
    integer (int_kind), dimension(:), allocatable, intent(in) :: &
       indxTij      ! flattened extended-global cell of each 1-D slot;
                    ! allocated on dst_task alone, hence allocatable
    integer (int_kind), intent(in) :: &
       nActive,   & ! 1-D slots to be filled
       nxext,     & ! extended global dimensions, matching indxTij's encoding
       nyext,     &
       dst_task
    type (distrb), intent(in) :: &
       src_dist

    real (dbl_kind), dimension(:), allocatable :: probe
    real (dbl_kind), dimension(:,:), allocatable :: G_probe
    integer (int_kind), dimension(:,:), allocatable :: tmpIW, sLst
    integer (int_kind), dimension(:), allocatable :: mySlot
    integer (int_kind) :: n, s, k, i, j, iw, iblk, nb, maxSlots, ierr, big
    real (dbl_kind) :: v
#ifndef NO_MPI
    integer (int_kind), dimension(:), allocatable :: snd_request
    integer (int_kind), dimension(:,:), allocatable :: snd_status
    integer (int_kind) :: nsends, status(MPI_STATUS_SIZE)
#endif
    character(len=*), parameter :: subname = '(slot_gather_setup)'

    nb       = size(map, dim=2)
    maxSlots = nx_block*ny_block
    big      = maxSlots + 1

    !-----------------------------------------------------------------
    ! stamp every slot with (global block, slot) and push it through the
    ! existing path
    !-----------------------------------------------------------------
    allocate(probe(maxval(off(1:nb)+cnt(1:nb))), stat=ierr)
    if (ierr /= 0) call abort_ice(subname//' ERROR: out of memory (probe)')
    probe = c0
    do n = 1, nblocks_tot
       if (src_dist%blockLocation(n) == my_task+1) then
          iblk = src_dist%blockLocalID(n)
          do s = 1, cnt(iblk)
             probe(off(iblk)+s) = real(n,dbl_kind)*real(big,dbl_kind) &
                                + real(s,dbl_kind)
          enddo
       endif
    enddo

    if (my_task == dst_task) then
       allocate(G_probe(nxext,nyext), stat=ierr)
    else
       allocate(G_probe(1,1), stat=ierr)
    endif
    if (ierr /= 0) call abort_ice(subname//' ERROR: out of memory (G_probe)')
    G_probe = c0

    call gather_globalR8_compressed(G_probe, probe, map, off, &
                                    dst_task, src_dist, c0, grid_ext=.true.)

    !-----------------------------------------------------------------
    ! dst_task: decode where each 1-D cell was fed from
    !-----------------------------------------------------------------
    allocate(sLst(maxSlots,nblocks_tot), stat=ierr)
    if (ierr /= 0) call abort_ice(subname//' ERROR: out of memory (sLst)')
    sLst = 0

    if (my_task == dst_task) then
       allocate(tmpIW(maxSlots,nblocks_tot), &
                sg_recvIW(maxSlots,nblocks_tot), &
                sg_fRecvIW(maxSlots,nblocks_tot), &
                sg_nfRecv(nblocks_tot), &
                sg_nRecv(nblocks_tot), stat=ierr)
       if (ierr /= 0) call abort_ice(subname//' ERROR: out of memory (recv)')
       tmpIW = 0
       sg_nRecv = 0

       do iw = 1, nActive
          j = (indxTij(iw) - 1)/nxext + 1
          i =  indxTij(iw) - (j-1)*nxext
          v = G_probe(i,j)
          if (v > c0) then
             n = int(v/real(big,dbl_kind))
             s = nint(v - real(n,dbl_kind)*real(big,dbl_kind))
             if (n < 1 .or. n > nblocks_tot .or. s < 1 .or. s > maxSlots) &
                call abort_ice(subname//' ERROR: probe decode out of range')
             tmpIW(s,n) = iw
          endif
       enddo

       ! compact to ascending slot order -- the order the sender will pack
       do n = 1, nblocks_tot
          k = 0
          do s = 1, maxSlots
             if (tmpIW(s,n) > 0) then
                k = k + 1
                sLst(k,n)      = s
                sg_recvIW(k,n) = tmpIW(s,n)
             endif
          enddo
          sg_nRecv(n) = k
       enddo
       deallocate(tmpIW)
    endif

    !-----------------------------------------------------------------
    ! hand each task the slot list for its own blocks.  One-off, so the
    ! message is the padded full-length list rather than a count exchange.
    !-----------------------------------------------------------------
    allocate(mySlot(maxSlots), sg_sendIdx(maxSlots,nb), sg_nSend(nb), &
             sg_sendIJ(maxSlots,nb), sg_fSend(maxSlots,nb), sg_nfSend(nb), &
             stat=ierr)
    if (ierr /= 0) call abort_ice(subname//' ERROR: out of memory (send)')
    sg_sendIdx = 0
    sg_sendIJ  = 0
    sg_nSend   = 0

#ifndef NO_MPI
    if (my_task == dst_task) then
       allocate(snd_request(nblocks_tot), &
                snd_status(MPI_STATUS_SIZE,nblocks_tot), stat=ierr)
       if (ierr /= 0) call abort_ice(subname//' ERROR: out of memory (req)')
       nsends = 0
    endif
#endif

    do n = 1, nblocks_tot
       if (my_task == dst_task) then
          if (src_dist%blockLocation(n) == my_task+1) then
             iblk = src_dist%blockLocalID(n)
             call store_send_list(sLst(:,n), iblk)
#ifndef NO_MPI
          else if (src_dist%blockLocation(n) > 0) then
             nsends = nsends + 1
             call MPI_ISEND(sLst(1,n), maxSlots, MPI_INTEGER, &
                            src_dist%blockLocation(n)-1, 4*mpitag_gs+n, &
                            MPI_COMM_ICE, snd_request(nsends), ierr)
#endif
          endif
#ifndef NO_MPI
       else if (src_dist%blockLocation(n) == my_task+1) then
          iblk = src_dist%blockLocalID(n)
          call MPI_RECV(mySlot, maxSlots, MPI_INTEGER, dst_task, &
                        4*mpitag_gs+n, MPI_COMM_ICE, status, ierr)
          call store_send_list(mySlot, iblk)
#endif
       endif
    enddo

#ifndef NO_MPI
    if (my_task == dst_task .and. nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
    if (allocated(snd_request)) deallocate(snd_request, snd_status)
#endif

    sg_maxSend = 0
    if (nb > 0) sg_maxSend = maxval(sg_nSend(1:nb))
    sg_ready   = .true.

    ! Once: what the compaction is actually worth here.  It is entirely
    ! domain-dependent -- it is the ratio of slots that feed the 1-D array
    ! to the full block volume the old path sent -- so a number from one
    ! grid says nothing about another.
    if (my_task == dst_task) then
       write(nu_diag,*) subname//' moves ', sum(sg_nRecv), ' of ', &
                        nblocks_tot*maxSlots, ' values per field'
    endif

    deallocate(probe, G_probe, sLst, mySlot)

  contains

    subroutine store_send_list(lst, lblk)
       integer (int_kind), dimension(:), intent(in) :: lst
       integer (int_kind), intent(in) :: lblk
       integer (int_kind) :: kk, ijc
       integer (int_kind), dimension(maxSlots) :: inv
       ! slot -> block-local cell, so slot_gather_filter can test this task's
       ! own mask at the right place without asking the master for anything
       inv = 0
       do ijc = 1, maxSlots
          if (map(ijc,lblk) > 0) inv(map(ijc,lblk)) = ijc
       enddo
       sg_nSend(lblk) = 0
       do kk = 1, maxSlots
          if (lst(kk) == 0) exit
          sg_nSend(lblk) = kk
          sg_sendIdx(kk,lblk) = off(lblk) + lst(kk)
          sg_sendIJ (kk,lblk) = inv(lst(kk))
       enddo
    end subroutine store_send_list

 end subroutine slot_gather_setup

!***********************************************************************

 subroutine slot_gather_filter(active_only, maskT, maskU, G_maskT, G_maskU, &
                               indxTij, nxext, dst_task, src_dist)

    use ice_fileunits, only: nu_diag

!  Narrow what gather_slots_r8 moves to the slots that carry ice.
!
!  Call once per step, AFTER the masks have been gathered and BEFORE the
!  data fields.  That ordering is the whole trick: both sides then decide
!  membership from the same mask -- the task from its own block copy, the
!  master from the gathered global one -- so the two lists agree entry for
!  entry and NO INDEX HAS TO TRAVEL with the data.
!
!  Both sides walk the existing ascending-slot lists in order and keep the
!  survivors in that order, so the k-th value sent is the k-th expected.
!
!  ⚠ Safe only because gather_slots_r8 pre-zeroes 1:nActive and the 1-D
!  kernels read these fields at their OWN cell exclusively -- every use in
!  ice_dyn_core1d is indexed (iw), none through nbrW/nbrS/... -- while
!  cycling on skipme at exactly the cells dropped here.  A field read at a
!  neighbour could not be filtered this way: it would read c0 where the
!  block path supplies a real value, the de0b9c3 fault class.

    logical (log_kind), intent(in) :: &
       active_only  ! .false. restores the full slot list
    logical (log_kind), dimension(:,:), intent(in) :: &
       maskT, maskU     ! this task's blocks, i and j flattened
    logical (log_kind), dimension(:,:), intent(in) :: &
       G_maskT, G_maskU ! extended-global, meaningful on dst_task alone
    integer (int_kind), dimension(:), allocatable, intent(in) :: &
       indxTij          ! dst_task only, hence allocatable
    integer (int_kind), intent(in) :: nxext, dst_task
    type (distrb), intent(in) :: src_dist

    integer (int_kind) :: n, k, kk, i, j, iw, ij, iblk
    character(len=*), parameter :: subname = '(slot_gather_filter)'

    if (.not. sg_ready) call abort_ice(subname//' ERROR: setup not run')

    sg_filter = active_only
    if (.not. active_only) return

    do n = 1, nblocks_tot
       if (src_dist%blockLocation(n) == my_task+1) then
          iblk = src_dist%blockLocalID(n)
          kk = 0
          do k = 1, sg_nSend(iblk)
             ij = sg_sendIJ(k,iblk)
             if (maskT(ij,iblk) .or. maskU(ij,iblk)) then
                kk = kk + 1
                sg_fSend(kk,iblk) = sg_sendIdx(k,iblk)
             endif
          enddo
          sg_nfSend(iblk) = kk
       endif
    enddo

    if (my_task == dst_task) then
       do n = 1, nblocks_tot
          kk = 0
          do k = 1, sg_nRecv(n)
             iw = sg_recvIW(k,n)
             j  = (indxTij(iw) - 1)/nxext + 1
             i  =  indxTij(iw) - (j-1)*nxext
             if (G_maskT(i,j) .or. G_maskU(i,j)) then
                kk = kk + 1
                sg_fRecvIW(kk,n) = iw
             endif
          enddo
          sg_nfRecv(n) = kk
       enddo

       ! Once: what the ice filter is worth on top of the slot compaction.
       if (sg_firstFilter) then
          sg_firstFilter = .false.
          write(nu_diag,*) subname//' moves ', sum(sg_nfRecv), ' of ', &
                           sum(sg_nRecv), ' slots (', &
                           nblocks_tot*nx_block*ny_block, ' before compaction)'
       endif
    endif

 end subroutine slot_gather_filter

!***********************************************************************

 subroutine gather_slots_r8(A_1d, A_slot, nActive, dst_task, src_dist)

!  Gather a compressed field straight into the destination 1-D array,
!  moving only the slots that feed it.  gather_globalR8_compressed sends
!  nx_block*ny_block per block whatever the layout holds; this sends
!  sg_nSend(iblk).
!
!  ⚠ A_1d is pre-zeroed over 1:nActive.  The old path reached a cell with
!  no slot through ARRAY_b = c0, so it read back c0; here such a cell is
!  simply never sent, and the zero has to be supplied instead.  Without
!  this the stale value of the previous step survives there.

    real (dbl_kind), dimension(:), allocatable, intent(inout) :: &
       A_1d         ! destination, dst_task alone -- hence allocatable
    real (dbl_kind), dimension(:), intent(in) :: &
       A_slot       ! compressed source field, dimension(dynTot)
    integer (int_kind), intent(in) :: nActive, dst_task
    type (distrb), intent(in) :: src_dist

    integer (int_kind) :: n, k, iblk, ierr
    integer (int_kind), dimension(:,:), pointer :: sndIdx, rcvIW
    integer (int_kind), dimension(:),   pointer :: nSnd, nRcv
    real (dbl_kind), dimension(:), allocatable :: rbuf
#ifndef NO_MPI
    real (dbl_kind), dimension(:,:), allocatable :: sbuf
    integer (int_kind), dimension(:), allocatable :: snd_request
    integer (int_kind), dimension(:,:), allocatable :: snd_status
    integer (int_kind) :: nsends, status(MPI_STATUS_SIZE)
#endif
    character(len=*), parameter :: subname = '(gather_slots_r8)'

    if (.not. sg_ready) call abort_ice(subname//' ERROR: setup not run')

    ! which descriptor set is live this step
    if (sg_filter) then
       nSnd => sg_nfSend ; sndIdx => sg_fSend
       if (my_task == dst_task) then
          nRcv => sg_nfRecv ; rcvIW => sg_fRecvIW
       endif
    else
       nSnd => sg_nSend  ; sndIdx => sg_sendIdx
       if (my_task == dst_task) then
          nRcv => sg_nRecv  ; rcvIW => sg_recvIW
       endif
    endif

    if (my_task == dst_task) then
       A_1d(1:nActive) = c0
       allocate(rbuf(max(1,maxval(nRcv))), stat=ierr)
       if (ierr /= 0) call abort_ice(subname//' ERROR: out of memory (rbuf)')
    endif

#ifndef NO_MPI
    if (my_task /= dst_task) then
       allocate(sbuf(max(1,sg_maxSend),size(nSnd)), &
                snd_request(size(nSnd)), &
                snd_status(MPI_STATUS_SIZE,size(nSnd)), stat=ierr)
       if (ierr /= 0) call abort_ice(subname//' ERROR: out of memory (sbuf)')
       nsends = 0
    endif
#endif

    do n = 1, nblocks_tot
       if (my_task == dst_task) then
          if (src_dist%blockLocation(n) == my_task+1) then
             iblk = src_dist%blockLocalID(n)
             do k = 1, nRcv(n)
                A_1d(rcvIW(k,n)) = A_slot(sndIdx(k,iblk))
             enddo
#ifndef NO_MPI
          else if (src_dist%blockLocation(n) > 0 .and. nRcv(n) > 0) then
             call MPI_RECV(rbuf, nRcv(n), mpiR8, &
                           src_dist%blockLocation(n)-1, 5*mpitag_gs+n, &
                           MPI_COMM_ICE, status, ierr)
             do k = 1, nRcv(n)
                A_1d(rcvIW(k,n)) = rbuf(k)
             enddo
#endif
          endif
#ifndef NO_MPI
       else if (src_dist%blockLocation(n) == my_task+1) then
          iblk = src_dist%blockLocalID(n)
          if (nSnd(iblk) > 0) then
             do k = 1, nSnd(iblk)
                sbuf(k,iblk) = A_slot(sndIdx(k,iblk))
             enddo
             nsends = nsends + 1
             call MPI_ISEND(sbuf(1,iblk), nSnd(iblk), mpiR8, dst_task, &
                            5*mpitag_gs+n, MPI_COMM_ICE, &
                            snd_request(nsends), ierr)
          endif
#endif
       endif
    enddo

#ifndef NO_MPI
    if (my_task /= dst_task) then
       if (nsends > 0) call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
       deallocate(sbuf, snd_request, snd_status)
    endif
#endif
    if (allocated(rbuf)) deallocate(rbuf)

 end subroutine gather_slots_r8

!***********************************************************************

 subroutine gather_globalL1_flat(ARRAY_G, ARRAY, dst_task, src_dist, fillValue, grid_ext)

   integer (int_kind), intent(in) :: &
      dst_task   ! task to which array should be gathered

   type (distrb), intent(in) :: &
      src_dist   ! distribution of blocks in the source array

   ! TARGET and CONTIGUOUS are required to alias this to rank 3 below
   logical (log_kind), dimension(:,:), intent(in), target, contiguous :: &
      ARRAY      ! array containing horizontal slab of distributed field

   logical (log_kind), intent(in), optional :: &
      fillValue  ! special value

   logical (log_kind), intent(in), optional :: &
      grid_ext   ! extended global grid flag

   logical (log_kind), dimension(:,:), intent(inout) :: &
      ARRAY_G    ! array containing global horizontal field on dst_task

   logical (log_kind), pointer :: &
      ARRAY_3d(:,:,:)   ! rank-3 alias of ARRAY

   character(len=*), parameter :: subname = '(gather_globalL1_flat)'

   if (size(ARRAY,dim=1) /= nx_block*ny_block) then
      call abort_ice(subname//' ERROR: first dimension is not nx_block*ny_block')
   endif

   ARRAY_3d(1:nx_block,1:ny_block,1:size(ARRAY,dim=2)) => ARRAY

   call gather_globalL1(ARRAY_G, ARRAY_3d, dst_task, src_dist, fillValue, grid_ext)

 end subroutine gather_globalL1_flat

!***********************************************************************

 subroutine scatter_globalL1(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type, fillValue, &
                               grid_ext)

!  Convert to integer and call scatter_globalI4
!  This subroutine scatters a global sized array from src_task to all
!  processors.
!  If grid_ext is true, then ARRAY_G is extended grid and updates outer halo
!  If grid_ext is false, then ARRAY_G is 1:nx_global,1:ny_global and outer
!    halo is updated for cyclic and tripole bcs.
!  If grid_ext is false and interior is true, then only the interior data will
!    be scattered, no outer halo bcs will be applied

   integer (int_kind), intent(in) :: &
      src_task        ! task from which array should be scattered

   type (distrb), intent(in) :: &
      dst_dist        ! distribution of resulting blocks

   logical (log_kind), dimension(:,:), intent(in) :: &
      ARRAY_G         ! array containing global field on src_task

   integer (int_kind), intent(in), optional :: &
      field_type,    &! id for type of field (scalar, vector, angle)
      field_loc       ! id for location on horizontal grid
                      !  (center, NEcorner, Nface, Eface)

   logical (log_kind), intent(in), optional :: &
      fillValue       ! special value

   logical (log_kind), intent(in), optional :: &
      grid_ext        ! extended grid flag

   logical (log_kind), dimension(:,:,:), intent(inout) :: &
      ARRAY           ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:,:), allocatable :: &
      IARRAY     ! integer array containing horizontal slab of distributed field

   integer (int_kind), dimension(:,:), allocatable :: &
      IARRAY_G   ! integer array containing global horizontal field on dst_task

   integer (int_kind) :: &
      ifillValue ! integer fill value

   character(len=*), parameter :: subname = '(scatter_globalL1)'

   allocate(IARRAY(size(ARRAY,dim=1),size(ARRAY,dim=2),size(ARRAY,dim=3)))
   allocate(IARRAY_G(size(ARRAY_G,dim=1),size(ARRAY_G,dim=2)))

   IARRAY_G = 0   ! 0 = false
   where (ARRAY_G) IARRAY_G = 1   ! 1 = true

   if (present(fillValue)) then
      ifillValue = 0
      if (fillValue) ifillValue = 1
      call scatter_global(IARRAY, IARRAY_G, src_task, dst_dist, &
                          field_loc, field_type, ifillValue, grid_ext)
   else
      call scatter_global(IARRAY, IARRAY_G, src_task, dst_dist, &
                          field_loc, field_type, grid_ext=grid_ext)
   endif

   ARRAY = .false.  ! 0 = false
   where (mod(abs(IARRAY),2) == 1) ARRAY = .true.  ! .ne.0 = true

!-----------------------------------------------------------------------

 end subroutine scatter_globalL1

!***********************************************************************

 end module ice_gather_scatter

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
