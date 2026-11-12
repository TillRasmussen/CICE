      program transportchk

      ! This tests state_to_tracers and tracers_to_state in ice_transport_driver.
      ! CICE_InitMod (from the standalone model) is used to initialize the
      ! tracer configuration.  The routines are then called on one block of
      ! hardwired inputs and the results are checked against known values:
      !   state_to_tracers fills aim and trm, trm is zero where aicen <= puny
      !   round trip state -> tracers -> state recovers the state where aicen > puny
      !   round trip where 0 < aicen <= puny keeps aicen, zeroes volumes
      !   tracers_to_state where ice is removed (aim = 0) zeroes aicen, vicen, vsnon
      !     and leaves trcrn unchanged
      !   tracers_to_state where ice is added to an ice-free cell
      !   tracers_to_state where a cell is ice-free before and after
      !   round trip for a block with no ice
      !
      ! ⚠ Both routines work in the COMPRESSED layout: their state arguments
      ! are this block's slot slice, nc+1 long -- the stored cells followed by
      ! the dead slot.  The hardwired state is still built from (i,j), because
      ! that is what makes the pattern cycle through the cases, but it is built
      ! AT each slot via dynIJ rather than over the whole block.  Cells with no
      ! slot are not part of either routine's contract and are not tested.
      !
      ! The dead slot is tested, because a stencil that lands on it must read
      ! c0, and because the two routines treat it differently: state_to_tracers
      ! writes it, tracers_to_state leaves it alone.

      use CICE_InitMod
      use ice_kinds_mod, only: int_kind, dbl_kind
      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: c0, c1, p5
      use ice_communicate, only: my_task, master_task, MPI_COMM_ICE
      use ice_domain_size, only: ncat, nslyr
      use ice_exit, only: abort_ice, end_run
      use ice_fileunits, only: nu_diag
      use ice_global_reductions, only: global_maxval
      use ice_grid, only: dynMap, dynOff, dynCnt, dynIJ
      use ice_transport_driver, only: state_to_tracers, tracers_to_state
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_query_parameters
      use icepack_intfc, only: icepack_query_tracer_sizes, icepack_query_tracer_indices

      implicit none

      integer(int_kind) :: i, j, n, it, ntrcr, ntrace, nt_qsno
      integer(int_kind) :: nc, off, ms, ij, nd
      integer(int_kind) :: cnt(3)
      real(dbl_kind) :: puny, rhos, Lfresh, errmax

      ! this block's slots: 1..nc are stored cells, nc+1 is the dead slot
      real(dbl_kind), allocatable :: &
         aice0(:), aicen(:,:), vicen(:,:), vsnon(:,:), trcrn(:,:,:), &
         aice0r(:), aicenr(:,:), vicenr(:,:), vsnonr(:,:), trcrnr(:,:,:), &
         aim(:,:), trm(:,:,:)

      ! the one block-shaped argument left: tracers_to_state still writes
      ! aice0 in both layouts, and that write is the only remaining user of
      ! dynMap in either routine
      real(dbl_kind), allocatable :: aice0b(:,:)

      integer(int_kind), parameter :: ntests = 8
      integer(int_kind) :: errorflag(ntests), errorflag0, gflag
      character(len=64) :: stringflag(ntests)
      integer(int_kind), parameter ::  &
         passflag = 0, &
         failflag = 1

      real(dbl_kind), parameter :: fillval = -1.0e36_dbl_kind
      real(dbl_kind), parameter :: errtol = 1.0e-14_dbl_kind
      real(dbl_kind), parameter :: aimadd = 0.2_dbl_kind   ! area added to ice-free cells
      character(len=*), parameter :: subname='(transportchk)'

      !-----------------------------------------------------------------
      ! Initialize CICE
      !-----------------------------------------------------------------

      call CICE_Initialize

      call icepack_query_parameters(puny_out=puny, rhos_out=rhos, Lfresh_out=Lfresh)
      call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
      call icepack_query_tracer_indices(nt_qsno_out=nt_qsno)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)
      ntrace = 2 + ntrcr

      nc  = dynCnt(1)
      off = dynOff(1)
      nd  = nc + 1              ! the dead slot

      errorflag = passflag
      stringflag(1) = 'state_to_tracers aim and trm'
      stringflag(2) = 'round trip where aicen > puny'
      stringflag(3) = 'round trip where 0 < aicen <= puny'
      stringflag(4) = 'tracers_to_state ice removed'
      stringflag(5) = 'tracers_to_state ice added to ice-free cell'
      stringflag(6) = 'tracers_to_state ice-free before and after'
      stringflag(7) = 'round trip with no ice in block'
      stringflag(8) = 'dead slot and the block aice0'

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) '=========================================================='
         write(6,*) ' '
         write(6,*) 'RunningUnitTest TRANSPORTCHK'
         write(6,*) ' '
         write(6,*) ' nx_block, ny_block = ', nx_block, ny_block
         write(6,*) ' ncat, ntrcr, nslyr = ', ncat, ntrcr, nslyr
         write(6,*) ' nt_qsno, puny      = ', nt_qsno, puny
         write(6,*) ' slots in block 1   = ', nc
         write(6,*) ' '
      endif

      if (nc < 1) call abort_ice(error_message=subname//' block 1 has no slots', &
         file=__FILE__, line=__LINE__)

      allocate(aice0 (nd), aicen (nd,ncat), vicen (nd,ncat), vsnon (nd,ncat), &
               trcrn (nd,ntrcr,ncat), &
               aice0r(nd), aicenr(nd,ncat), vicenr(nd,ncat), vsnonr(nd,ncat), &
               trcrnr(nd,ntrcr,ncat), &
               aim(nd,0:ncat), trm(nd,ntrace,ncat), &
               aice0b(nx_block,ny_block))

      !-----------------------------------------------------------------
      ! Hardwired state: cycle through ice-free, aicen < puny, aicen = puny,
      ! and normal ice, with some snow-free cells
      !-----------------------------------------------------------------

      call init_state(aice0, aicen, vicen, vsnon, trcrn)

      !-----------------------------------------------------------------
      ! Test 1: state_to_tracers
      ! aim is the area, trm is zero where aicen <= puny, otherwise
      ! hice, hsno, snow enthalpy + rhos*Lfresh, other tracers
      !-----------------------------------------------------------------

      aim = fillval
      trm = fillval
      call state_to_tracers (off+1, off+nc, ntrcr, ntrace, &
                             aice0, aicen, trcrn, vicen, vsnon, aim, trm)

      if (any(aim(:,0) /= aice0)) errorflag(1) = failflag
      if (any(aim(1:nc,1:ncat) /= aicen(1:nc,:))) errorflag(1) = failflag
      errmax = c0
      do n = 1, ncat
      do ms = 1, nc
         if (aicen(ms,n) > puny) then
            errmax = max(errmax, relerr(trm(ms,1,n), vicen(ms,n)/aicen(ms,n)))
            errmax = max(errmax, relerr(trm(ms,2,n), vsnon(ms,n)/aicen(ms,n)))
            do it = 1, ntrcr
               if (trm(ms,2+it,n) /= trcrn(ms,it,n) + snwoff_of(it)) errorflag(1) = failflag
            enddo
         else
            if (any(trm(ms,:,n) /= c0)) errorflag(1) = failflag
         endif
      enddo
      enddo
      if (errmax > errtol) errorflag(1) = failflag
      if (my_task == master_task) write(6,*) ' test 1 max relative error hice, hsno = ', errmax

      ! state_to_tracers writes the dead slot: the areas and tracers it
      ! carries must be c0, whatever the state arrays held there
      if (any(aim(nd,1:ncat) /= c0)) errorflag(8) = failflag
      if (any(trm(nd,:,:) /= c0))    errorflag(8) = failflag

      !-----------------------------------------------------------------
      ! Tests 2 and 3: round trip state -> tracers -> state
      ! The state arrays start from fill values, so every stored cell must
      ! be set (trcrn is not set where aim = 0, so it is not checked there).
      ! tracers_to_state does NOT write the dead slot, which is checked below.
      !-----------------------------------------------------------------

      aice0r = fillval
      aicenr = fillval
      vicenr = fillval
      vsnonr = fillval
      trcrnr = fillval
      aice0b = fillval
      call tracers_to_state (off, nc, nx_block, ny_block, dynMap(:,1), ntrcr, ntrace, &
                             aim, trm, aice0b, aice0r, aicenr, trcrnr, vicenr, vsnonr)

      if (any(aice0r(1:nc) /= aice0(1:nc))) errorflag(2) = failflag
      if (any(aicenr(1:nc,:) /= aicen(1:nc,:))) errorflag(2) = failflag
      errmax = c0
      cnt = 0
      do n = 1, ncat
      do ms = 1, nc
         if (aicen(ms,n) > puny) then
            cnt(1) = cnt(1) + 1
            errmax = max(errmax, relerr(vicenr(ms,n), vicen(ms,n)))
            errmax = max(errmax, relerr(vsnonr(ms,n), vsnon(ms,n)))
            do it = 1, ntrcr
               errmax = max(errmax, relerr(trcrnr(ms,it,n), trcrn(ms,it,n)))
            enddo
         elseif (aicen(ms,n) > c0) then
            ! tiny ice keeps its area, volumes and tracers are zeroed
            cnt(2) = cnt(2) + 1
            if (vicenr(ms,n) /= c0 .or. vsnonr(ms,n) /= c0) errorflag(3) = failflag
            do it = 1, ntrcr
               if (trcrnr(ms,it,n) /= c0 - snwoff_of(it)) errorflag(3) = failflag
            enddo
         else
            cnt(3) = cnt(3) + 1
            if (vicenr(ms,n) /= c0 .or. vsnonr(ms,n) /= c0) errorflag(2) = failflag
         endif
      enddo
      enddo
      if (errmax > errtol) errorflag(2) = failflag
      if (cnt(1) == 0) errorflag(2) = failflag
      if (cnt(2) == 0) errorflag(3) = failflag
      if (my_task == master_task) then
         write(6,*) ' test 2 cells aicen > puny, max relative error = ', cnt(1), errmax
         write(6,*) ' test 3 cells 0 < aicen <= puny = ', cnt(2)
         write(6,*) ' test 2 cells aicen = 0 = ', cnt(3)
      endif

      ! the dead slot is outside the loop bound tracers_to_state uses, so it
      ! must still hold the fill value
      if (aice0r(nd) /= fillval)       errorflag(8) = failflag
      if (any(aicenr(nd,:) /= fillval)) errorflag(8) = failflag
      if (any(vicenr(nd,:) /= fillval)) errorflag(8) = failflag
      if (any(vsnonr(nd,:) /= fillval)) errorflag(8) = failflag

      ! the block copy of aice0 is written through dynMap, at every cell that
      ! owns a slot and nowhere else
      do ms = 1, nc
         ij = dynIJ(off+ms)
         j  = (ij-1)/nx_block + 1
         i  = ij - (j-1)*nx_block
         if (aice0b(i,j) /= aim(ms,0)) errorflag(8) = failflag
      enddo

      !-----------------------------------------------------------------
      ! Tests 4, 5 and 6: tracers_to_state after a change in aim similar
      ! to transport.  aim, trm are from state_to_tracers above, then
      !   ice removed (aim = 0) from some cells with aicen > puny
      !   ice added (aim = aimadd) to some cells with aicen = 0
      !   other cells with aicen = 0 stay ice-free
      ! The state arrays start from the original state.
      !-----------------------------------------------------------------

      do n = 1, ncat
      do ms = 1, nc
         ij = dynIJ(off+ms)
         j  = (ij-1)/nx_block + 1
         i  = ij - (j-1)*nx_block
         if (aicen(ms,n) > puny .and. mod(i + 2*j + n, 5) == 0) then
            aim(ms,n) = c0
         elseif (aicen(ms,n) == c0 .and. mod(i+j,2) == 0) then
            aim(ms,n) = aimadd
            trm(ms,1,n) = 1.1_dbl_kind
            trm(ms,2,n) = 0.05_dbl_kind
            do it = 1, ntrcr
               trm(ms,2+it,n) = -0.5_dbl_kind*real(it,dbl_kind) + snwoff_of(it)
            enddo
         endif
      enddo
      enddo

      aice0r = aice0
      aicenr = aicen
      vicenr = vicen
      vsnonr = vsnon
      trcrnr = trcrn
      call tracers_to_state (off, nc, nx_block, ny_block, dynMap(:,1), ntrcr, ntrace, &
                             aim, trm, aice0b, aice0r, aicenr, trcrnr, vicenr, vsnonr)

      cnt = 0
      do n = 1, ncat
      do ms = 1, nc
         if (aicen(ms,n) > c0 .and. aim(ms,n) == c0) then
            ! ice removed: area and volumes zeroed, tracers unchanged
            cnt(1) = cnt(1) + 1
            if (aicenr(ms,n) /= c0 .or. vicenr(ms,n) /= c0 .or. vsnonr(ms,n) /= c0) &
               errorflag(4) = failflag
            if (any(trcrnr(ms,:,n) /= trcrn(ms,:,n))) errorflag(4) = failflag
         elseif (aicen(ms,n) == c0 .and. aim(ms,n) > c0) then
            ! ice added to an ice-free cell
            cnt(2) = cnt(2) + 1
            if (aicenr(ms,n) /= aimadd) errorflag(5) = failflag
            if (vicenr(ms,n) /= aimadd*1.1_dbl_kind)  errorflag(5) = failflag
            if (vsnonr(ms,n) /= aimadd*0.05_dbl_kind) errorflag(5) = failflag
            do it = 1, ntrcr
               if (trcrnr(ms,it,n) /= trm(ms,2+it,n) - snwoff_of(it)) errorflag(5) = failflag
            enddo
         elseif (aicen(ms,n) == c0 .and. aim(ms,n) == c0) then
            ! ice-free before and after: area and volumes zero, tracers unchanged
            cnt(3) = cnt(3) + 1
            if (aicenr(ms,n) /= c0 .or. vicenr(ms,n) /= c0 .or. vsnonr(ms,n) /= c0) &
               errorflag(6) = failflag
            if (any(trcrnr(ms,:,n) /= trcrn(ms,:,n))) errorflag(6) = failflag
         endif
      enddo
      enddo
      if (any(aice0r(1:nc) /= aim(1:nc,0))) errorflag(4) = failflag
      if (cnt(1) == 0) errorflag(4) = failflag
      if (cnt(2) == 0) errorflag(5) = failflag
      if (cnt(3) == 0) errorflag(6) = failflag
      if (my_task == master_task) then
         write(6,*) ' test 4 cells ice removed = ', cnt(1)
         write(6,*) ' test 5 cells ice added = ', cnt(2)
         write(6,*) ' test 6 cells ice-free before and after = ', cnt(3)
      endif

      !-----------------------------------------------------------------
      ! Test 7: round trip for a block with no ice, tracers keep their
      ! (arbitrary) values, everything else is zero
      !-----------------------------------------------------------------

      aicen = c0
      vicen = c0
      vsnon = c0
      aice0 = c1
      aim = fillval
      trm = fillval
      call state_to_tracers (off+1, off+nc, ntrcr, ntrace, &
                             aice0, aicen, trcrn, vicen, vsnon, aim, trm)
      if (any(aim(:,0) /= c1) .or. any(aim(1:nc,1:ncat) /= c0)) errorflag(7) = failflag
      if (any(trm /= c0)) errorflag(7) = failflag

      aice0r = fillval
      aicenr = fillval
      vicenr = fillval
      vsnonr = fillval
      trcrnr = trcrn
      call tracers_to_state (off, nc, nx_block, ny_block, dynMap(:,1), ntrcr, ntrace, &
                             aim, trm, aice0b, aice0r, aicenr, trcrnr, vicenr, vsnonr)
      if (any(aice0r(1:nc) /= c1)) errorflag(7) = failflag
      if (any(aicenr(1:nc,:) /= c0) .or. any(vicenr(1:nc,:) /= c0) &
          .or. any(vsnonr(1:nc,:) /= c0)) errorflag(7) = failflag
      if (any(trcrnr /= trcrn)) errorflag(7) = failflag

      !-----------------------------------------------------------------
      ! Report
      !-----------------------------------------------------------------

      errorflag0 = maxval(errorflag)
      gflag = global_maxval(errorflag0, MPI_COMM_ICE)
      errorflag0 = gflag
      do n = 1, ntests
         gflag = global_maxval(errorflag(n), MPI_COMM_ICE)
         errorflag(n) = gflag
      enddo

      if (my_task == master_task) then
         write(6,*) ' '
         write(6,*) 'TRANSPORTCHK COMPLETED SUCCESSFULLY'
         do n = 1, ntests
            if (errorflag(n) == passflag) then
               write(6,*) 'PASS ',trim(stringflag(n))
            else
               write(6,*) 'FAIL ',trim(stringflag(n))
            endif
         enddo
         if (errorflag0 == passflag) then
            write(6,*) 'TRANSPORTCHK TEST COMPLETED SUCCESSFULLY'
         else
            write(6,*) 'TRANSPORTCHK TEST FAILED'
         endif
         write(6,*) ' '
         write(6,*) '=========================================================='
         write(6,*) ' '
      endif

      !-----------------------------------------------------------------
      ! Gracefully end
      !-----------------------------------------------------------------

      call end_run()

      contains

!=======================================================================
! Hardwired state, built at each slot of block 1.  The pattern is the same
! function of (i,j,n) as the block form used; dynIJ recovers (i,j) from the
! slot, so the cases still cycle cell by cell.  The dead slot carries c0.

      subroutine init_state(aice0, aicen, vicen, vsnon, trcrn)

      real(dbl_kind), intent(out) :: aice0(:), aicen(:,:), vicen(:,:), &
                                     vsnon(:,:), trcrn(:,:,:)

      integer(int_kind) :: i, j, n, it, ms, ij

      aicen = c0
      vicen = c0
      vsnon = c0
      trcrn = c0

      do n = 1, ncat
      do ms = 1, nc
         ij = dynIJ(off+ms)
         j  = (ij-1)/nx_block + 1
         i  = ij - (j-1)*nx_block
         select case (mod(i + 3*j + 5*n, 8))
         case (0)
            aicen(ms,n) = c0
         case (1)
            aicen(ms,n) = p5*puny
         case (2)
            aicen(ms,n) = puny
         case (3)
            aicen(ms,n) = c0
         case default
            aicen(ms,n) = 0.01_dbl_kind + 0.1_dbl_kind*real(mod(i*j+n,9),dbl_kind)
         end select
         vicen(ms,n) = aicen(ms,n)*(0.3_dbl_kind + 0.7_dbl_kind*real(n,dbl_kind) &
                                  + 0.01_dbl_kind*real(i,dbl_kind))
         if (mod(i+j,3) == 0) then
            vsnon(ms,n) = c0
         else
            vsnon(ms,n) = aicen(ms,n)*0.03_dbl_kind*real(mod(j,7)+1,dbl_kind)
         endif
         do it = 1, ntrcr
            if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
               trcrn(ms,it,n) = -rhos*Lfresh - 1.0e5_dbl_kind*real(it+n,dbl_kind)
            else
               trcrn(ms,it,n) = -1.7_dbl_kind + 0.37_dbl_kind*real(it,dbl_kind) &
                              + 0.013_dbl_kind*real(i-j,dbl_kind) + 0.11_dbl_kind*real(n,dbl_kind)
            endif
         enddo
      enddo
      enddo
      aice0(:) = c1 - sum(aicen,dim=2)

      end subroutine init_state

!=======================================================================
! Offset added to snow enthalpy tracers by state_to_tracers

      real(dbl_kind) function snwoff_of(it)

      integer(int_kind), intent(in) :: it

      if (it >= nt_qsno .and. it < nt_qsno+nslyr) then
         snwoff_of = rhos*Lfresh
      else
         snwoff_of = c0
      endif

      end function snwoff_of

!=======================================================================

      real(dbl_kind) function relerr(a, b)

      real(dbl_kind), intent(in) :: a, b

      if (b == c0) then
         relerr = abs(a)
      else
         relerr = abs(a-b)/abs(b)
      endif

      end function relerr

!=======================================================================

      end program transportchk

!=======================================================================
