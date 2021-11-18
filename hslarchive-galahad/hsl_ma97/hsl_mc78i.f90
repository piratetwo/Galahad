! (c) STFC 2010-2011
! Originating author: Jonathan Hogg
!
! Given a pivot order, this package performs common tasks
! required in the analyse phase of a symmetric sparse direct solver.
! Either the entire analyse may be performed or individual tasks.
! The matrix may be hled in assembled form or in elemental form.
!
! Version 1.4.0
! See ChangeLog for version history

! To convert to long:
! s/_integer/_long_integer
! Set pkg_type to long
module hsl_mc78_integer
   implicit none

   private
   public :: mc78_control
   public :: mc78_analyse, mc78_supervars, mc78_compress_by_svar, mc78_etree, &
      mc78_elt_equiv_etree, mc78_postorder, mc78_col_counts, mc78_supernodes, &
      mc78_stats, mc78_row_lists, mc78_optimize_locality

   integer, parameter :: dp = kind(0d0) ! not package type
   integer, parameter :: long = selected_int_kind(18)

   integer, parameter :: minsz_ms = 16 ! minimum size to use merge sort

   integer, parameter :: pkg_type = kind(0) ! package type - integer or long

   type mc78_control
      integer :: heuristic = 1 ! 1=ma77 2=cholmod
      integer :: nrelax(3) = (/ 4, 16, 48 /) ! CHOLMOD-like
      real(dp) :: zrelax(3) = (/ 0.8, 0.1, 0.05 /) ! CHOLMOD-like
      integer :: nemin = 16  ! Node amalgamation parameter

      integer :: unit_error = 6
      integer :: unit_warning = 6
      logical :: ssa_abort = .false. ! If .true., then return with an error if
         ! an assembled matrix is detected as symbolically singular (we do
         ! not garuntee to detect _all_ symbolically singular matrices).
         ! If .false., then a warning is raised instead.

      logical :: svar = .false. ! If .true. then supervariables are used in
         ! the assembled case, otherwise they are not. Supervaraibles are
         ! always used in the elemental case.
      logical :: sort = .false. ! If .true. then entries within each supernode's
         ! row lists are sorted. Otherwise they might not be.
      logical :: lopt = .false. ! If .true. then variable ordering is optimized
         ! for cache locality. Otherwise it is not.
   end type mc78_control

   integer, parameter :: MC78_ERROR_ALLOC = -1 ! allocate error
   integer, parameter :: MC78_ERROR_SSA   = -2 ! symbolically singular assembled
   integer, parameter :: MC78_ERROR_ROW_SMALL = -3 ! supplied row array to short
   integer, parameter :: MC78_ERROR_UNKNOWN = -99 ! internal/unknown error

   ! Warning flags are treated as bit masks, add together if multiple occour
   integer, parameter :: MC78_WARNING_SSA = 1 ! symbolically singular assembled
   integer, parameter :: MC78_WARNING_BLK_SVAR = 2 ! svar and blk pivs requested

   interface mc78_analyse
      module procedure mc78_analyse_assembled_integer
      module procedure mc78_analyse_elemental_integer
   end interface mc78_analyse

   interface mc78_supervars
      module procedure mc78_supervars_integer
   end interface mc78_supervars

   interface mc78_compress_by_svar
      module procedure mc78_compress_by_svar_integer
   end interface mc78_compress_by_svar

   interface mc78_etree
      module procedure mc78_etree_integer
   end interface mc78_etree

   interface mc78_elt_equiv_etree
      module procedure mc78_elt_equiv_etree_integer
   end interface

   interface mc78_postorder
      ! Note: cannot distinguish postorder_std between integer and long versions
      module procedure mc78_postorder_std
      module procedure mc78_postorder_detect
   end interface mc78_postorder

   interface mc78_col_counts
      module procedure mc78_col_counts_integer
   end interface mc78_col_counts

   ! Note: cannot distinguish mc78_supernodes between integer and long versions
   ! Note: cannot distinguish mc78_stats between integer and long versions

   interface mc78_row_lists
      module procedure mc78_row_lists_nosvar_integer
      module procedure mc78_row_lists_svar_integer
   end interface mc78_row_lists

   ! Note: cannot distinguish mc78_optimize_locality between integer and long
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            Main analysis routines   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! For assembled matrix input, this subroutine performs a full analysis.
! This is essentially a wrapper around the rest of the package.
!
! Performance might be improved by:
! * Improving the sort algorithm used in find_row_idx
!
subroutine mc78_analyse_assembled_integer(n, ptr, row, perm, nnodes, sptr, &
      sparent, rptr, rlist, control, info, stat, nfact, nflops, piv_size)
   integer, intent(in) :: n ! Dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! Column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! Row indices
   integer, dimension(n), intent(inout) :: perm
      ! perm(i) must hold position of i in the pivot sequence. 
      ! On exit, holds the pivot order to be used by factorization.
   integer, intent(out) :: nnodes ! number of supernodes found
   integer, dimension(:), allocatable, intent(out) :: sptr ! supernode pointers
   integer, dimension(:), allocatable, intent(out) :: sparent ! assembly tree
   integer(long), dimension(:), allocatable, intent(out) :: rptr
      ! pointers to rlist
   integer, dimension(:), allocatable, intent(out) :: rlist ! row lists
   ! For details of control, info : see derived type descriptions
   type(mc78_control), intent(in) :: control
   integer, intent(out) :: info
   integer, optional, intent(out) :: stat
   integer(long), optional, intent(out) :: nfact
   integer(long), optional, intent(out) :: nflops
   integer, dimension(n), optional, intent(inout) :: piv_size ! If
      ! present, then matches matrix order and specifies block pivots. 
      ! piv_size(i) is the number of entries pivots in the block pivot
      ! containing column i.

   integer(pkg_type), dimension(:), allocatable :: bptr ! copy of matrix with
      ! added entries for block pivots - column pointers
   integer, dimension(:), allocatable :: brow ! copy of matrix with added
      ! entries for block pivots - row indices
   integer :: flag ! return status flag for call to compress_by_svar
   integer(pkg_type) :: sz
   integer, dimension(:), allocatable :: svara ! array for supervariables
   integer :: i
   integer, dimension(:), allocatable :: invp ! inverse permutation of perm
   integer, dimension(:), allocatable :: sinvp
   integer :: j
   integer :: k
   integer :: nsvar ! number of supervariables
   integer, dimension(:), allocatable :: sperm
   integer(pkg_type), dimension(:), allocatable :: ptr2
   integer :: realn ! number of variables with an actual entry present
   integer, dimension(:), allocatable :: row2
   integer :: st ! stat argument in allocate calls
   logical :: svar_r

   integer :: svar_type ! 0=none, 1=col, 2=compressed form

   integer, dimension(:), allocatable :: scc

   ! initialise
   info = 0

   svar_r = control%svar

   ! Ensure allocatable output arguments are deallocated
   deallocate(sptr, stat=st)
   deallocate(sparent, stat=st)
   deallocate(rptr, stat=st)
   deallocate(rlist, stat=st)

   ! Initialise inverse permutation and check for duplicates
   allocate(invp(n), stat=st)
   if(st.ne.0) goto 490
   do i = 1, n
      j = perm(i)
      invp(j) = i
   end do

   if(present(piv_size)) then
      call convert_to_blk_piv(n, invp, piv_size)
      allocate(bptr(n+1), brow(ptr(n+1)-1+2*n), stat=st)
      if(st.ne.0) goto 490
      call mc78_block_prep(n, ptr, row, bptr, brow, perm, invp, piv_size, &
         st)
      if(st.ne.0) goto 490
      if(svar_r) then
         ! Supervariables don't interact well with block pivots, so don't do it
         svar_r = .false.
         info = info + MC78_WARNING_BLK_SVAR
      endif
   endif

   ! Determine supervariables (if required)
   if(svar_r) then
      allocate(svara(n), stat=st)
      if(st.ne.0) goto 490
      realn = n
      call mc78_supervars(realn, ptr, row, perm, invp, nsvar, svara, st)
      if(st.ne.0) goto 490
      if(n.ne.realn) then
         if(control%ssa_abort) then
            if(control%unit_error.gt.0) write(control%unit_error, "(2a)") &
               "HSL_MC78: Error, matrix is symbolically singular and ", &
               "control%ssa_abort=.true.."
            info = MC78_ERROR_SSA
            return
         else
            if(control%unit_warning.gt.0) write(control%unit_warning, "(a)") &
               "HSL_MC78: Warning, matrix is symbolically singular."
            info = info + MC78_WARNING_SSA
         endif
      endif
      if(3*nsvar.lt.2*n) then
         svar_type = 2 ! Use compressed form
      else
         svar_type = 0 ! Do not use supervariables
      endif

      select case(svar_type)
      case(0) ! do not use supervariables
         ! release resources
         deallocate(svara, stat=st)
      case(2) ! Compressed form
         ! It is worth using the compressed form
         ! Determine upper bound on size of data for compressed array
         sz = 0
         k = 1
         do i = 1, nsvar
            j = invp(k)
            sz = sz + ptr(j+1) - ptr(j)
            k = k + svara(i)
         end do
         allocate(ptr2(nsvar+1), row2(sz), sperm(nsvar), sinvp(nsvar), stat=st)
         if(st.ne.0) goto 490
         call mc78_compress_by_svar(n, ptr, row, invp, nsvar, svara, &
            ptr2, sz, row2, flag, st)
         select case(flag)
         case(0) ! Everything OK
            ! Do nothing
         case(-1) ! Allocate failure
            goto 490
         case default ! Should never happen
            info = MC78_ERROR_UNKNOWN
            return
         end select
         ! Compressed matrix is in pivot order already
         do i = 1, nsvar
            sperm(i) = i
            sinvp(i) = i
         end do
      end select
   else
      svar_type = 0
      realn = n ! Assume full rank
   endif

   select case(svar_type)
   case(0)
      if(present(piv_size)) then
         call mc78_inner_analyse(n, realn, bptr, brow, perm, invp, nnodes, &
            sptr, sparent, scc, rptr, rlist, control, info, st, &
            block_pivots=piv_size)
      else
         call mc78_inner_analyse(n, realn, ptr, row, perm, invp, nnodes, &
            sptr, sparent, scc, rptr, rlist, control, info, st)
      endif
   case(2)
      call mc78_inner_analyse(nsvar, i, ptr2, row2, sperm, sinvp, nnodes, &
         sptr, sparent, scc, rptr, rlist, control, info, st, wt=svara, &
         block_pivots=piv_size)
      if(st.ne.0) goto 490
      if(info.lt.0) return
      if(i.ne.nsvar) then
         ! Note: This code should NEVER execute
         if(control%unit_error.gt.0) &
            write(control%unit_error, "(a,2(a,i8))") "MC78_ANALYSE Internal ", &
               "Error: supervariable matrix is rank deficient: i = ", i, &
               "nsvar = ", nsvar
         info = MC78_ERROR_UNKNOWN
         return
      endif
      call svar_unmap(n, nsvar, svara, perm, invp, nnodes, sinvp, sptr, st)
   end select
   if(st.ne.0) goto 490
   if(info.lt.0) return

   ! Calculate info%num_factor and info%num_flops
   call mc78_stats(nnodes, sptr, scc, nfact=nfact, nflops=nflops)

   if(control%lopt) then
      ! Reorder elimination variables for better cache locality
      if(control%sort) then
         call mc78_optimize_locality(n, realn, perm, invp, nnodes, sptr, &
            sparent, rptr, rlist, st, sort=.true.)
      else
         call mc78_optimize_locality(n, realn, perm, invp, nnodes, sptr, &
            sparent, rptr, rlist, st)
      endif
      if(st.ne.0) goto 490
   else if(control%sort) then
      call dbl_tr_sort(n, nnodes, rptr, rlist, st)
      if(st.ne.0) goto 490
   endif

   if(present(piv_size)) &
      call convert_from_blk_piv(n, perm, piv_size)

   !print *, "n, realn = ", n, realn
   !print *, "ptr = ", ptr
   !do i = 1, n
   !   print *, "row(", i, ") = ", row(ptr(i):ptr(i+1)-1)
   !end do
   !print *, "nnodes = ", nnodes
   !print *, "sptr = ", sptr(1:nnodes+1)
   !print *, "sparent = ", sparent(1:nnodes)
   !print *, "rptr = ", rptr(1:nnodes+1)
   !print *, "rlist = ", rlist(1:rptr(nnodes+1)-1)
   !print *, "invp = ", invp
   !print *, "perm = ", perm
   !print *, "piv_size = ", piv_size

   return

   !!!!!!!!!!!!!!!!!!
   ! Error handlers !
   !!!!!!!!!!!!!!!!!!

   490 continue
   info = MC78_ERROR_ALLOC
   if(present(stat)) stat = st
   return
end subroutine mc78_analyse_assembled_integer

!
! Inner core for assembled analyse routine, used to make calls in compressed
! (supervariable) case and standard case uniform
!
subroutine mc78_inner_analyse(n, realn, ptr, row, perm, invp, nnodes, sptr, &
      sparent, scc, rptr, rlist, control, info, st, wt, block_pivots)
   integer, intent(in) :: n ! Dimension of system
   integer, intent(out) :: realn ! Symbolic dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! Column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! Row indices
   integer, dimension(n), intent(inout) :: perm
      ! perm(i) must hold position of i in the pivot sequence. 
      ! On exit, holds the pivot order to be used by factorization.
   integer, dimension(n), intent(inout) :: invp ! inverse of perm
   integer, intent(out) :: nnodes ! number of supernodes found
   integer, dimension(:), allocatable, intent(out) :: sptr ! supernode pointers
   integer, dimension(:), allocatable, intent(out) :: sparent ! assembly tree
   integer, dimension(:), allocatable, intent(out) :: scc ! supernodal col cnt
   integer(long), dimension(:), allocatable, intent(out) :: rptr
      ! pointers to rlist
   integer, dimension(:), allocatable, intent(out) :: rlist ! row lists
   ! For details of control, info : see derived type descriptions
   type(mc78_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(out) :: st
   integer, dimension(n), optional, intent(in) :: wt ! Weights of columns
      ! (i.e. size of each supervariable they represent)
   integer, dimension(n), optional, intent(inout) :: block_pivots ! If
      ! present, then matches pivot order and specifies block pivots. 
      ! block_pivots(i) may take one of the following values:
      ! 0 - pivot i is in the middle of a block pivot
      ! 1 - pivot i is the first pivot of a block pivot
      ! 2 - pivot i is the last pivot of a block pivot
      ! 3 - pivot i is a 1x1 pivot

   integer, dimension(:), allocatable :: cc ! number of entries in each column
   integer :: ntot ! total number of variables
   integer, dimension(:), allocatable :: parent ! parent of each node in etree
   integer, dimension(:), allocatable :: tperm ! temporary permutation vector

   ! Build elimination tree
   allocate(parent(n), stat=st)
   if(st.ne.0) return
   call mc78_etree(n, ptr, row, perm, invp, parent, st)
   if(st.ne.0) return

   ! Postorder tree (modifies perm!)
   call mc78_postorder(n, realn, ptr, perm, invp, parent, st, block_pivots)
   if(st.ne.0) return

   if(n.ne.realn) then
      if(control%ssa_abort) then
         if(control%unit_error.gt.0) write(control%unit_error, "(2a)") &
            "HSL_MC78: Error, matrix is symbolically singular and ", &
            "control%ssa_abort=.true.."
         info = MC78_ERROR_SSA
         return
      else
         if(control%unit_warning.gt.0) write(control%unit_warning, "(a)") &
            "HSL_MC78: Warning, matrix is symbolically singular."
         info = info + MC78_WARNING_SSA
      endif
   endif

   ! Determine column counts
   allocate(cc(n+1), stat=st)
   if(st.ne.0) return
   call mc78_col_counts(n, ptr, row, perm, invp, parent, cc, st, wt=wt)
   if(st.ne.0) return

   ! Identify supernodes
   allocate(tperm(n), sptr(n+1), sparent(n), scc(n), stat=st)
   if(st.ne.0) return
   call mc78_supernodes(n, realn, parent, cc, tperm, nnodes, sptr, sparent, &
      scc, invp, control, info, st, wt=wt, block_pivots=block_pivots)
   if(info.lt.0) return

   ! Apply permutation to obtain final elimination order
   call apply_perm(n, tperm, perm, invp, cc, block_pivots=block_pivots)

   ! Determine column patterns - keep%nodes(:)%index
   allocate(rptr(nnodes+1), rlist(sum_long(scc(1:nnodes))), stat=st)
   if(st.ne.0) return
   if(present(wt)) then
      ntot = sum(wt)
      call mc78_row_lists(n, wt, ntot, ptr, row, perm, invp, nnodes, sptr, &
         sparent, scc, rptr, rlist, control, info, st)
   else
      call mc78_row_lists(n, ptr, row, perm, invp, nnodes, sptr, &
         sparent, scc, rptr, rlist, control, info, st)
   endif
   if(st.ne.0) return
end subroutine mc78_inner_analyse

!
! This subroutine performs full analyse when A is in elemental form.
! This is essentially a wrapper around the rest of the package.
!
subroutine mc78_analyse_elemental_integer(n, nelt, starts, vars, perm, &
      eparent, nnodes, sptr, sparent, rptr, rlist, control, info, stat, &
      nfact, nflops, piv_size)
   integer, intent(in) :: n ! Maximum integer used to index an element
   integer, intent(in) :: nelt ! Number of elements
   integer(pkg_type), dimension(nelt+1), intent(in) :: starts ! Element pointers
   integer, dimension(starts(nelt+1)-1), intent(in) :: vars ! Variables
      !assoicated with each element. Element i has vars(starts(i):starts(i+1)-1)
   integer, dimension(n), intent(inout) :: perm
      ! perm(i) must hold position of i in the pivot sequence. 
      ! On exit, holds the pivot order to be used by factorization.
   integer, dimension(nelt), intent(out) :: eparent ! On exit, eparent(i) holds
      ! node of assembly that element i is a child of.
   integer, intent(out) :: nnodes ! number of supernodes found
   integer, dimension(:), allocatable, intent(out) :: sptr ! supernode pointers
   integer, dimension(:), allocatable, intent(out) :: sparent ! assembly tree
   integer(long), dimension(:), allocatable, intent(out) :: rptr
      ! pointers to rlist
   integer, dimension(:), allocatable, intent(out) :: rlist ! row lists
   ! For details of control, info : see derived type descriptions
   type(mc78_control), intent(in) :: control
   integer, intent(out) :: info
   integer, optional, intent(out) :: stat
   integer(long), optional, intent(out) :: nfact ! If present, then on exit
      ! contains the number of entries in L
   integer(long), optional, intent(out) :: nflops ! If present, then on exit
      ! contains the number of floating point operations in factorize.
   integer, dimension(n), optional, intent(inout) :: piv_size ! If
      ! present, then matches matrix order and specifies block pivots. 
      ! piv_size(i) is the number of entries pivots in the block pivot
      ! containing column i.

   integer, dimension(:), allocatable :: cc ! number of entries in each column
   integer :: i
   integer, dimension(:), allocatable :: invp ! inverse permutation of perm
   integer :: j
   integer :: nsvar ! number of supervariables
   integer, dimension(:), allocatable :: parent ! parent of each node in etree
   integer, dimension(:), allocatable :: perm2 ! temporary permutation vector
   integer(pkg_type), dimension(:), allocatable :: ptr ! column pointers for
      ! equivilent matrix
   integer :: realn ! Set to actual number of variables present
   integer, dimension(:), allocatable :: row ! row indices for equivilent matrix
   integer :: st ! stat argument in allocate calls
   integer, dimension(:), allocatable :: svara ! array for supervariables
   integer, dimension(:), allocatable :: sinvp ! inverse permutation of svars
   integer, dimension(:), allocatable :: sperm ! permutation vector of svars
   integer(long) :: sz ! temporary var for size of arrays at allocation

   integer, dimension(:), allocatable :: scc

   ! initialise
   info = 0

   ! Ensure allocatable output arguments are deallocated
   deallocate(sptr, stat=st)
   deallocate(sparent, stat=st)
   deallocate(rptr, stat=st)
   deallocate(rlist, stat=st)

   ! Initialise inverse permutation and check for duplicates
   allocate(invp(n), stat=st)
   if(st.ne.0) goto 490
   do i = 1, n
      j = perm(i)
      invp(j) = i
   end do

   if(present(piv_size)) &
      call convert_to_blk_piv(n, invp, piv_size)

   ! Determine supernodes, build equivilant lwr matrix and find elimination tree
   sz = starts(nelt+1)-1
   if(present(piv_size)) sz = sz + n
   allocate(ptr(n+2), row(sz), svara(n+1), parent(n), stat=st)
   if(st.ne.0) goto 490
   realn = n
   call mc78_elt_equiv_etree(realn, nelt, starts, vars, perm, invp, nsvar, &
      svara, ptr, row, eparent, parent, st, block_pivots=piv_size)
   if(st.ne.0) goto 490

   ! Set up permutations of supervariables (initially the idenity)
   allocate(sperm(nsvar), sinvp(nsvar), stat=st)
   if(st.ne.0) goto 490
   sperm(1:nsvar) = (/ (i, i=1,nsvar) /)
   sinvp(1:nsvar) = (/ (i, i=1,nsvar) /)

   ! Postorder tree (modifies perm!)
   call mc78_postorder(nsvar, sperm, sinvp, parent, st, &
      block_pivots=piv_size)
   if(st.ne.0) goto 490

   ! Determine column counts
   allocate(cc(nsvar+1), stat=st)
   if(st.ne.0) goto 490
   call mc78_col_counts(nsvar, ptr, row, sperm, sinvp, parent, cc, st, wt=svara)
   if(st.ne.0) goto 490

   ! Identify supernodes
   allocate(perm2(nsvar), sptr(nsvar+1), sparent(nsvar), scc(nsvar), stat=st)
   if(st.ne.0) goto 490
   call mc78_supernodes(nsvar, nsvar, parent, cc, perm2, nnodes, sptr, &
      sparent, scc, sinvp, control, info, st, wt=svara, block_pivots=piv_size)
   if(info.eq.MC78_ERROR_ALLOC) goto 490
   if(info.lt.0) return

   ! Apply permutation to obtain final elimination order
   call apply_perm(nsvar, perm2, sperm, sinvp, cc, block_pivots=piv_size)

   ! Determine column patterns - keep%nodes(:)%index
   allocate(rptr(nnodes+1), rlist(sum_long(scc(1:nnodes))), stat=st)
   if(st.ne.0) goto 490
   call mc78_row_lists(nsvar, svara, n, ptr, row, sperm, sinvp, &
      nnodes, sptr, sparent, scc, rptr, rlist, control, info, st)
   if(st.ne.0) goto 490

   ! Unmap from supervariables to real variables
   call svar_unmap(n, nsvar, svara, perm, invp, nnodes, sinvp, &
      sptr, st)
   if(st.ne.0) goto 490

   ! Calculate info%num_factor and info%num_flops
   call mc78_stats(nnodes, sptr, scc, nfact=nfact, nflops=nflops)

   if(control%lopt) then
      ! Reorder elimination variables for better cache locality
      if(control%sort) then
         call mc78_optimize_locality(n, realn, perm, invp, nnodes, sptr, &
            sparent, rptr, rlist, st, sort=.true.)
      else
         call mc78_optimize_locality(n, realn, perm, invp, nnodes, sptr, &
            sparent, rptr, rlist, st)
      endif
      if(st.ne.0) goto 490
   else if(control%sort) then
      call dbl_tr_sort(n, nnodes, rptr, rlist, st)
      if(st.ne.0) goto 490
   endif

   ! Adjust eparent with final permuation. Use svara to contain a mapping
   ! from original variables to supernodes
   do i = 1, nnodes
      do j = sptr(i), sptr(i+1)-1
         svara(invp(j)) = i
      end do
   end do
   svara(n+1) = n+1
   do i = 1, nelt
      eparent(i) = svara(eparent(i))
   end do

   if(present(piv_size)) &
      call convert_from_blk_piv(n, perm, piv_size)

   !print *, "n, realn, nelt = ", n, realn, nelt
   !print *, "starts = ", starts
   !print *, "vars = ", vars
   !print *, "nnodes = ", nnodes
   !print *, "sptr = ", sptr(1:nnodes+1)
   !print *, "sparent = ", sparent(1:nnodes)
   !print *, "eparent = ", eparent(:)
   !print *, "rptr = ", rptr(1:nnodes+1)
   !print *, "rlist = ", rlist(1:rptr(nnodes+1)-1)
   !print *, "invp = ", invp
   !print *, "perm = ", perm

   return

   !!!!!!!!!!!!!!!!!!
   ! Error handlers !
   !!!!!!!!!!!!!!!!!!

   490 continue
   info = MC78_ERROR_ALLOC
   if(present(stat)) stat = st
   return
end subroutine mc78_analyse_elemental_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Supervariable routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine find supervariables of A using the algorithm of [1].
!
! [1] I.S. Duff and J.K. Reid. "Exploiting zeros on the diagonal in the direct
!     solution of indefinite sparse symmetric linear systems".
!     ACM TOMS 22(2), pp227-257. 1996.
!
subroutine mc78_supervars_integer(n, ptr, row, perm, invp, nsvar, svar, st)
   integer, intent(inout) :: n ! Dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! Column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! Row indices
   integer, dimension(n), intent(inout) :: perm
      ! perm(i) must hold position of i in the pivot sequence. 
      ! On exit, holds the pivot order to be used by factorization.
   integer, dimension(n), intent(inout) :: invp ! inverse of perm
   integer, intent(out) :: nsvar ! number of supervariables
   integer, dimension(n), intent(out) :: svar ! number of vars in each svar
   integer, intent(out) :: st

   logical :: full_rank ! flags if supervariable 1 has ever become empty.
      ! If it has not, then the varaibles in s.v. 1 are those that never
      ! occur
   integer :: i
   integer(long) :: ii
   integer :: j
   integer :: idx ! current index
   integer :: next_sv ! head of free sv linked list
   integer :: nsv ! new supervariable to move j to
   integer :: piv ! current pivot
   integer :: col ! current column of A
   integer :: sv ! current supervariable
   integer :: svc ! temporary holding supervariable count
   integer, dimension(:), allocatable :: sv_new  ! Maps each supervariable to
      ! a new supervariable with which it is associated.
   integer, dimension(:), allocatable :: sv_seen ! Flags whether svariables have
      ! been seen in the current column. sv_seen(j) is set to col when svar j
      ! has been encountered.
   integer, dimension(:), allocatable :: sv_count ! number of variables in sv.

   allocate(sv_new(n+1), sv_seen(n+1), sv_count(n+1), stat=st)
   if(st.ne.0) return

   svar(:) = 1
   sv_count(1) = n
   sv_seen(1) = 0

   ! Setup linked list of free super variables
   next_sv = 2
   do i = 2, n
      sv_seen(i) = i+1
   end do
   sv_seen(n+1) = -1

   ! Determine supervariables using modified Duff and Reid algorithm
   full_rank = .false.
   do col = 1, n
      if(ptr(col+1).ne.ptr(col)) then
         ! If column is not empty, add implicit diagonal entry
         j = col
         sv = svar(j)
         if(sv_count(sv).eq.1) then ! Are we only (remaining) var in sv
            full_rank = full_rank .or. (sv.eq.1)
            ! MUST BE the first time that sv has been seen for this
            ! column, so just leave j in sv, and go to next variable.
            ! (Also there can be no other vars in this block pivot)
         else
            ! There is at least one other variable remaining in sv
            ! MUST BE first occurence of sv in the current row/column,
            ! so define a new supervariable and associate it with sv.
            sv_seen(sv) = col
            sv_new(sv) = next_sv
            nsv = next_sv
            next_sv = sv_seen(next_sv)
            sv_new(nsv) = nsv ! avoids problems with duplicates
            sv_seen(nsv) = col
            ! Now move j from sv to nsv
            nsv = sv_new(sv)
            svar(j) = nsv
            sv_count(sv) = sv_count(sv) - 1
            sv_count(nsv) = 1
            ! This sv cannot be empty as initial sv_count was > 1
         endif
      endif
      do ii = ptr(col), ptr(col+1)-1
         j = row(ii)
         sv = svar(j)
         if(sv_count(sv).eq.1) then ! Are we only (remaining) var in sv
            full_rank = full_rank .or. (sv.eq.1)
            ! If so, and this is first time that sv has been seen for this
            ! column, then we can just leave j in sv, and go to next variable.
            if(sv_seen(sv).lt.col) cycle
            ! Otherwise, we have already defined a new supervariable associated
            ! with sv. Move j to this variable, then retire (now empty) sv.
            nsv = sv_new(sv)
            if(sv.eq.nsv) cycle
            svar(j) = nsv
            sv_count(nsv) = sv_count(nsv) + 1
            ! Old sv is now empty, add it to top of free stack
            sv_seen(sv) = next_sv
            next_sv = sv
         else
            ! There is at least one other variable remaining in sv
            if(sv_seen(sv).lt.col) then
               ! this is the first occurence of sv in the current row/column,
               ! so define a new supervariable and associate it with sv.
               sv_seen(sv) = col
               sv_new(sv) = next_sv
               sv_new(next_sv) = next_sv ! avoids problems with duplicates
               next_sv = sv_seen(next_sv)
               sv_count(sv_new(sv)) = 0
               sv_seen(sv_new(sv)) = col
            endif
            ! Now move j from sv to nsv
            nsv = sv_new(sv)
            svar(j) = nsv
            sv_count(sv) = sv_count(sv) - 1
            sv_count(nsv) = sv_count(nsv) + 1
            ! This sv cannot be empty as sv_count was > 1
         endif
      end do
   end do

   ! Note: block pivots do not mix well with supervariables as any significant
   ! number (unless aligned to s.v.s) will demolish any gain from using them.
   ! Converting vlock pivots to s.v.s results in potentially large amount of
   ! unneeded fillin to left of block pivot.
   !! If block pivots are being used, we force all pivots of a block pivot
   !! to be in either the same supervariable, or in supervariables of size 1
   !if(present(block_pivots)) then
   !   piv = 1
   !   do while(piv.le.n)
   !      ! Check if we need to split pivots
   !      split = .false.
   !      sv = svar(piv)
   !      do i = piv+1, piv+block_pivots(piv)
   !         j = invp(i)
   !         if(svar(j).ne.sv) then
   !            split = .true.
   !            exit
   !         endif
   !      end do
   !      ! Do split if required
   !      if(split) then
   !         j = invp(i)
   !         do i = piv, piv+block_pivots(piv)
   !            sv = svar(j)
   !            if(sv_count(sv).eq.1) cycle ! Already a singleton
   !            ! Otherwise create a new sv and move j to it
   !            nsv = next_sv
   !            next_sv = sv_seen(next_sv)
   !            svar(j) = nsv
   !            sv_count(nsv) = 1
   !            sv_count(sv) = sv_count(sv) - 1
   !         end do
   !      endif
   !      piv = piv + block_pivots(piv) + 1
   !   end do
   !endif

   ! Now modify pivot order such that all variables in each supervariable are
   ! consecutive. Do so by iterating over pivots in elimination order. If a
   ! pivot has not already been listed, then order that pivot followed by
   ! any other pivots in that supervariable.

   ! We will build a new inverse permutation in invp, and then find perm
   ! afterwards. First copy invp to perm:
   perm(:) = invp(:)
   ! Next we iterate over the pivots that have not been ordered already
   ! Note: as we begin, all entries of sv_seen are less than or equal to n+1
   ! hence we can use <=n+1 or >n+1 as a flag to indicate if a variable has been
   ! ordered.
   idx = 1
   nsvar = 0
   do piv = 1, n
      if(sv_seen(piv).gt.n+1) cycle ! already ordered
      ! Record information for supervariable
      sv = svar(perm(piv))
      if(.not.full_rank .and. sv.eq.1) cycle ! Don't touch unused vars
      nsvar = nsvar + 1
      svc = sv_count(sv)
      sv_new(nsvar) = svc ! store # vars in s.v. to copy to svar later
      j = piv
      ! Find all variables that are members of sv and order them.
      do while(svc.gt.0)
         do j = j, n
            if(svar(perm(j)).eq.sv) exit
         end do
         sv_seen(j) = n+2 ! flag as ordered
         invp(idx) = perm(j)
         idx = idx + 1
         svc = svc - 1
         j = j + 1
      end do
   end do
   ! Push unused variables to end - these are those vars still in s.v. 1
   if(.not.full_rank) then
      svc = sv_count(1)
      ! Find all variables that are members of sv and order them.
      j = 1
      do while(svc.gt.0)
         do j = j, n
            if(svar(perm(j)).eq.1) exit
         end do
         invp(idx) = perm(j)
         idx = idx + 1
         svc = svc - 1
         j = j + 1
      end do
      n = n - sv_count(1)
   end if
   ! Recover perm as inverse of invp
   do piv = 1, n
      perm(invp(piv)) = piv
   end do
   ! sv_new has been used to store number of variables in each svar, copy into
   ! svar where it is returned.
   svar(1:nsvar) = sv_new(1:nsvar)
end subroutine mc78_supervars_integer

!
! This subroutine takes a set of supervariables and compresses the supplied
! matrix using them.
!
! As we would need a full scan of the matrix to calculate the correct size of
! row2, we instead allow the user to make a guess at a good size and return
! an error if this turns out to be incorrect. An upper bound on the required
! size may be obtained by summing the number of entries in the first column of
! each supervariable.
!
! Error returns:
!   MC78_ERROR_ALLOC      Failed to allocate memory
!   MC78_ERROR_ROW_SMALL  row2 too small
subroutine mc78_compress_by_svar_integer(n, ptr, row, invp, nsvar, svar, ptr2, &
      lrow2, row2, info, st)
   integer, intent(in) :: n ! Dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! Column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! Row indices
   integer, dimension(n), intent(in) :: invp ! inverse of perm
   integer, intent(in) :: nsvar
   integer, dimension(nsvar), intent(in) :: svar ! super variables of A
   integer(pkg_type), dimension(nsvar+1), intent(out) :: ptr2
   integer(pkg_type), intent(in) :: lrow2
   integer, dimension(lrow2), intent(out) :: row2
   integer, intent(out) :: info
   integer, intent(out) :: st

   integer :: piv, svc, sv, col
   integer(pkg_type) :: j, idx
   integer, dimension(:), allocatable :: flag, sv_map

   info = 0 ! by default completed succefully

   allocate(flag(nsvar), sv_map(n), stat=st)
   if(st.ne.0) then
      info = MC78_ERROR_ALLOC
      return
   endif
   flag(:) = 0

   ! Setup sv_map
   piv = 1
   do svc = 1, nsvar
      do piv = piv, piv + svar(svc) - 1
         sv_map( invp(piv) ) = svc
      end do
   end do

   piv = 1
   idx = 1
   do svc = 1, nsvar
      col = invp(piv)
      ptr2(svc) = idx
      do j = ptr(col), ptr(col+1)-1
         sv = sv_map(row(j))
         if(flag(sv).eq.piv) cycle ! Already dealt with this supervariable
         if(idx.gt.lrow2) then
            ! oops, row2 is too small
            info = MC78_ERROR_ROW_SMALL
            return
         endif
         ! Add row entry for this sv
         row2(idx) = sv
         flag(sv) = piv
         idx = idx + 1
      end do
      piv = piv + svar(svc)
   end do
   ptr2(svc) = idx
end subroutine mc78_compress_by_svar_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Elimination tree routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine determines the elimination tree of a PAP^T where A is a
! sparse symmetric matrix stored in compressed sparse column form with
! entries both above and below the diagonal present in the argument matrix.
! P is a permutation stored in order such that order(i) gives the pivot
! position of column i. i.e. order(3) = 5 means that the fifth pivot is
! A_33.
!
! The elimination tree is returned in the array parent. parent(i) gives the
! parent in the elimination tree of pivot i.
!
! The algorithm used is that of Liu [1].
!
! [1] Liu, J. W. 1986. A compact row storage scheme for Cholesky factors using
!     elimination trees. ACM TOMS 12, 2, 127--148.
!
subroutine mc78_etree_integer(n, ptr, row, perm, invp, parent, st)
   integer, intent(in) :: n ! dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! column pointers of A
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! row indices of A
   integer, dimension(n), intent(in) :: perm ! perm(i) is the pivot position
      ! of column i
   integer, dimension(n), intent(in) :: invp ! inverse of perm
   integer, dimension(n), intent(out) :: parent ! parent(i) is the
      ! parent of pivot i in the elimination tree
   integer, intent(out) :: st ! stat parmeter for allocate calls

   integer(pkg_type) :: i ! next index into row
   integer :: j ! current entry in row
   integer :: k ! current ancestor
   integer :: l ! next ancestor
   integer :: piv ! current pivot
   integer :: rowidx ! current column of A = invp(piv)
   integer, dimension(:), allocatable :: vforest ! virtual forest, used for
      ! path compression (shortcuts to top of each tree)

   ! Allocate virtual forest and initialise it
   allocate(vforest(n), stat=st)
   if(st.ne.0) return
   vforest(:) = n+1

   ! Loop over rows of A in pivot order
   piv = 1
   do while(piv.le.n)
      !print *, "row ", piv
      rowidx = invp(piv)
      ! Loop over entries in row in lower triangle of PAP^T
      do i = ptr(rowidx), ptr(rowidx+1)-1
         j = perm(row(i))
         if(j.ge.piv) cycle ! not in lower triangle
         !print *, "  entry ", j
         k = j
         do while(vforest(k).lt.piv)
            l = vforest(k)
            vforest(k) = piv
            k = l
         end do
         ! Check if we have already done this pivot
         if(vforest(k).eq.piv) cycle 
         parent(k) = piv
         vforest(k) = piv
      end do
      parent(piv) = n + 1 ! set to be a root if not overwritten
      piv = piv + 1 ! move on to next pivot
   end do
end subroutine mc78_etree_integer

!
! This subroutine identifies supervariables of A using a modified variant of
! the algorithm of Duff and Reid [1]. A lower triangular equivilant matrix
! is returned that is expressed in terms of these supervariables. The grouping
! of variables into supervaribles is returned through a modified pivot order
! and an array specifying the number of variables in each supervariable in
! elimination order. Finally the vector eparent is also returned. This contains
! the variable (in natural numbering) that corresponds to the least pivot in
! each supervariable.
!
! [1] I.S. Duff and J.K. Reid. "Exploiting zeros on the diagonal in the direct
!     solution of indefinite sparse symmetric linear systems".
!     ACM TOMS 22(2), pp227-257. 1996.
!
! Note: If block pivots are present they have priority over supervariables - 
! members of same block pivot must remain in same supervariable. This is
! enforced by moving them all at once.
subroutine mc78_elt_equiv_etree_integer(n, nelt, starts, vars, perm, invp, &
      nsvar, svar, ptr, row, eparent, parent, st, block_pivots)
   integer, intent(inout) :: n ! dimension of system
   integer, intent(in) :: nelt ! number of elements
   integer(pkg_type), dimension(nelt+1), intent(in) :: starts ! variable
      ! pointers of elements
   integer, dimension(starts(nelt+1)-1), intent(in) :: vars ! variables of
      ! elements
   integer, dimension(n), intent(inout) :: perm ! perm(i) is the pivot position
      ! of column i
   integer, dimension(n), intent(inout) :: invp ! inverse of perm
   integer, intent(out) :: nsvar ! number of supervariables found
   integer, dimension(n), intent(out) :: svar ! size of each supervariable
   integer(pkg_type), dimension(n+1), intent(out) :: ptr ! column pointers
      ! for equivilant lower triangular form
   integer, dimension(:), intent(out) :: row ! row indices
      ! for equivilant lower triangular form
   integer, dimension(nelt), intent(out) :: eparent ! parent nodes of each
      ! element - i.e. the least pivot in each element
   integer, dimension(n), intent(out) :: parent ! parent(i) is parent of node
      ! i in the elimination tree
   integer, intent(out) :: st ! stat parmeter for allocate calls
   integer, dimension(n), optional, intent(inout) :: block_pivots ! Used for
      ! block pivots, see description in analyse phase.

   integer :: csv ! column supervariable in loop
   integer :: elt ! current element
   logical :: full_rank ! flags if supervariable 1 has ever become empty.
      ! If it has not, then the varaibles in s.v. 1 are those that never
      ! occour
   integer :: i
   integer(pkg_type) :: ii
   integer(pkg_type) :: idx ! Next insert position
   integer :: j
   integer :: k
   integer :: minpiv ! minimum pivot of current element
   integer, dimension(:), allocatable :: mp_head, mp_next ! mp_head and mp_next
      ! store a linked list of elements for which a given variable is the
      ! minimum pivot.
   integer :: next_sv ! Top of stack of free supervariables (stored as a linked
      ! list in unused part of sv_seen)
   integer :: orign ! original system dimension
   integer :: nsv ! temporary variable storing new supervariable to move var to
   integer :: piv ! current pivot
   integer :: sv ! current supervariable
   integer :: svc ! temporary variable storing supervariable count remaining
   integer, dimension(:), allocatable :: sv_count ! sv_count(s) is the number
      ! of variables in supervariable s.
   integer, dimension(:), allocatable :: sv_map ! sv_map(v) is the current
      ! supervariable to which variable v belongs.
   integer, dimension(:), allocatable :: sv_new ! sv_map(s) is new
      ! supervariable for variables currently in supervariable s.
   integer, dimension(:), allocatable :: sv_seen ! sv_seen(s) is used to flag
      ! if supervariable s has been found in the current element before.
      ! In addition the part corresponding to unused supervariables is used
      ! to store a stack (as a linked list) of empty supervariables.
   integer(pkg_type), dimension(:), allocatable :: uprptr ! column pointers for
      ! upper triangular equivilant form.
   integer, dimension(:), allocatable :: uprrow ! row indices for upper
      ! triangular equivilant form.
   logical :: used ! flag if a variable has been used

   ! Initialise supervariable representation
   allocate(sv_new(max(nelt+1,n+1)), sv_seen(max(nelt,n)+1), &
      sv_map(max(nelt,n)+1), sv_count(max(nelt,n)+1), stat=st)
   if(st.ne.0) return
   sv_map(:) = 1 ! All vars are intially in supervariable 1
   sv_count(1) = n ! ... which thus has all variables
   sv_seen(1) = 0 ! Flag supervariable 1 as unseen on first iteration
   orign = n

   if(present(block_pivots)) then
      ! Do not mix supervariables and block pivots

      ! Still need to determine minimum pivots and rank
      sv_seen(:) = 0
      do elt = 1, nelt
         minpiv = n+1
         do ii = starts(elt), starts(elt+1)-1
            j = vars(ii)
            ! Mark variable as used
            sv_seen(j) = 1
            ! Determine minimum pivot
            minpiv = min(minpiv, perm(j))
         end do
         ! Store the minimum pivot as original variable (avoids messy remapping)
         if(minpiv.le.n) eparent(elt) = invp(minpiv)
      end do

      ! Build invp that pushes unsued vars to the end. Be careful of unused vars
      ! that are in fact part of block pivots, split them out.
      perm(:) = invp(:)
      piv = 1
      j = 1
      ! Handle variables that are actually used
      do while(piv.le.n)
         used = .true.
         do i = piv, n
            used = used .and. (sv_seen(perm(i)).eq.1)
            if(block_pivots(i).ge.2) exit ! end of block pivot
         end do

         if(used) then
            ! Block pivot is entirely composed of used variables
            do piv = piv, i
               invp(j) = perm(piv)
               sv_seen(perm(piv)) = 2
               j = j + 1
            end do
         else
            ! Block pivot has some unused variables in it
            k = 0
            do piv = piv, i
               if(sv_seen(perm(piv)).eq.1) then
                  invp(j) = perm(piv)
                  sv_seen(perm(piv)) = 2
                  j = j + 1
                  if(k.eq.0) then
                     ! This is the new start of the block pivot
                     select case(block_pivots(piv))
                     case(0) ! was in the middle. now a start
                        block_pivots(piv) = 1
                     case(2) ! was the end. now a 1x1
                        block_pivots(piv) = 3
                     end select
                  endif
                  k = piv
               endif
            end do
            if(k.ne.0) then
               ! The was at least one used variable in the block pivot
               select case(block_pivots(k))
               case(0) ! was the middle. now an end
                  block_pivots(k) = 2
               case(1) ! was the start. now a 1x1
                  block_pivots(k) = 3
               end select
            endif
         endif
         piv = i + 1
      end do
      ! Handle unused variables; build supervariables
      nsvar = 0
      do piv = 1, n
         i = perm(piv)
         if(sv_seen(i).eq.0) then
            invp(j) = i
            j = j + 1
            block_pivots(piv) = 3 ! Force to 1x1
         else
            nsvar = nsvar + 1
            svar(nsvar) = 1
            sv_new(i) = nsvar
         endif
      end do
      ! Map block_pivots in original variable order into sv_map
      do i = 1, n
         sv_map(perm(i)) = block_pivots(i)
      end do
      ! Map sv_map in new pivot order back into block_pivots
      do i = 1, n
         block_pivots(i) = sv_map(invp(i))
      end do
      ! Reestablish perm
      do i = 1, n
         perm(invp(i)) = i
      end do
   else
      ! Setup linked list of free supervariables - we utilise the unused part
      ! of sv_seen for this
      next_sv = 2
      do i = 2, n
         sv_seen(i) = i+1
      end do
      sv_seen(n+1) = -1
      
      ! Determine supervariables. At the same time find the least pivot
      ! associated with each element.
      nsvar = 1
      full_rank = .false.
      do elt = 1, nelt
         minpiv = n+1
         do ii = starts(elt), starts(elt+1)-1
            j = vars(ii)
            ! Determine minimum pivot
            minpiv = min(minpiv, perm(j))
            sv = sv_map(j)
            if(sv_count(sv).eq.1) then ! Are we only (remaining) var in sv
               full_rank = full_rank .or. (sv.eq.1)
               ! If so, and this is first time that sv has been seen for this
               ! element, then we can just leave j in sv, and go to next
               ! variable.
               if(sv_seen(sv).lt.elt) cycle
               ! Otherwise, we have already defined a new supervariable
               ! associated with sv. Move j to this variable, then retire (now
               ! empty) sv.
               ! Note: as only var in sv, cannot have fellows in block pivot
               nsv = sv_new(sv)
               if(sv.eq.nsv) cycle  ! don't delete a variable because of a
                                    ! duplicate
               sv_map(j) = nsv
               sv_count(nsv) = sv_count(nsv) + 1
               ! Old sv is now empty, add it to top of free stack
               sv_seen(sv) = next_sv
               next_sv = sv
               nsvar = nsvar - 1
            else
               ! There is at least one other variable remaining in sv
               if(sv_seen(sv).lt.elt) then
                  ! this is the first occurence of sv in the current element,
                  ! so define a new supervariable and associate it with sv.
                  sv_seen(sv) = elt
                  sv_new(sv) = next_sv
                  sv_new(next_sv) = next_sv  ! ensure we are tolerant of
                                             ! duplicates
                  next_sv = sv_seen(next_sv)
                  sv_count(sv_new(sv)) = 0
                  sv_seen(sv_new(sv)) = elt
                  nsvar = nsvar + 1
               endif
               ! Now move j from sv to nsv
               nsv = sv_new(sv)
               sv_map(j) = nsv
               sv_count(sv) = sv_count(sv) - 1
               sv_count(nsv) = sv_count(nsv) + 1
               ! We know sv can't be empty, so it doesn't need adding to free
               ! stack
            endif
         end do
         ! Store the minimum pivot as original variable (avoids messy remapping)
         if(minpiv.le.n) then
            eparent(elt) = invp(minpiv)
         else
            eparent(elt) = n+1
         endif
      end do

      ! Now modify pivot order such that all variables in each supervariable are
      ! consecutive. Do so by iterating over pivots in elimination order. If a
      ! pivot has not already been listed, then order that pivot followed by
      ! any other pivots in that supervariable.

      ! We will build a new inverse permutation in invp, and then find perm
      ! afterwards. First copy invp to perm:
      perm(:) = invp(:)
      ! Next we iterate over the pivots that have not been ordered already
      ! Note: as we begin, all entries of sv_seen are less than or equal to n+1
      ! hence we can use <=n+1 or >n+1 as a flag to indicate if a variable has
      ! been ordered.
      idx = 1
      nsvar = 0
      do piv = 1, n
         if(sv_seen(piv).gt.n+1) cycle ! already ordered
         ! Record information for supervariable
         sv = sv_map(perm(piv))
         if(.not.full_rank .and. sv.eq.1) cycle ! Don't touch unused vars
         svc = sv_count(sv)
         nsvar = nsvar + 1
         svar(nsvar) = svc
         ! Find all variables that are members of sv and order them.
         j = piv
         do while(svc.gt.0)
            do j = j, n
               if(sv_map(perm(j)).eq.sv) exit
            end do
            sv_seen(j) = n+2 ! flag as ordered
            sv_new(perm(j)) = nsvar ! new mapping to sv
            invp(idx) = perm(j)
            idx = idx + 1
            svc = svc - 1
            j = j + 1
         end do
      end do
      sv_new(n+1) = nsvar+1
      ! Push unused variables to end - these are those vars still in s.v. 1
      if(.not.full_rank) then
         svc = sv_count(1)
         ! Find all variables that are members of sv and order them.
         j = 1
         do while(svc.gt.0)
            do j = j, n
               if(sv_map(perm(j)).eq.1) exit
            end do
            invp(idx) = perm(j)
            idx = idx + 1
            svc = svc - 1
            j = j + 1
         end do
         n = n - sv_count(1)
      end if
      ! Recover perm as inverse of invp
      do piv = 1, n
         perm(invp(piv)) = piv
      end do
   endif

   ! build linked lists by supervariable
   allocate(mp_head(nsvar+1), mp_next(nelt), stat=st)
   if(st.ne.0) return
   mp_head(:) = -1
   do elt = 1, nelt
      if(eparent(elt).gt.orign) cycle
      minpiv = sv_new(eparent(elt))
      if(present(block_pivots) .and. minpiv.ne.1) then
         do while(block_pivots(minpiv-1).lt.2)
            minpiv = minpiv - 1
            if(minpiv.eq.1) exit
         end do
      endif
      ! Store element in linked list for minpiv
      mp_next(elt) = mp_head(minpiv)
      mp_head(minpiv) = elt
   end do

   ! Iterate over columns in pivot order, storing the lower triangular
   ! equivilant matrix as we go. At the same time, build the column counts for
   ! the upper triangle in uprptr, but offset by 2 (ie uprptr(i+2) for col i).
   ! Observe that all the pivots associated with the supervariable
   ! to which minpiv belongs _must_ appear in each element that minpiv does.
   ! Note: This only generates the lower triangular part of the matrix!
   allocate(uprptr(nsvar+2), stat=st)
   if(st.ne.0) return
   uprptr(:) = 0
   sv_seen(:) = 0
   idx = 1
   ptr(:) = -1
   do csv = 1, nsvar
      elt = mp_head(csv)
      ptr(csv) = idx
      sv_seen(csv) = csv ! Mark diagonal as seen, as it is implicit.
      do while(elt.ne.-1)
         do ii = starts(elt), starts(elt+1)-1
            sv = sv_new(vars(ii))
            ! Skip this sv if it is already included (or is implicit)
            if(sv_seen(sv).ge.csv) cycle
            sv_seen(sv) = csv ! Mark as seen
            ! If we can't skip it, then add entry (sv,csv) to lwr matrix
            row(idx) = sv
            idx = idx + 1
            ! Add count in upper triangle for (csv,sv)
            uprptr(sv+2) = uprptr(sv+2) + 1
         end do
         ! Move on to next element for which this is the minimum pivot
         elt = mp_next(elt)
      end do
      if(present(block_pivots) .and. csv.ne.nsvar) then
         ! Add entry (csv+1,csv) to ensure elimination tree correct
         if(block_pivots(csv).lt.2 .and. sv_seen(csv+1).ne.csv) then
            sv_seen(csv+1) = csv
            row(idx) = csv+1
            idx = idx + 1
            ! Add count in upper triangle for (csv, csv+1)
            uprptr(csv+1+2) = uprptr(csv+1+2) + 1
         endif
      endif
   end do
   ptr(nsvar+1) = idx

   ! Build upper form - work out column start for col i in uprptr(i+1)
   uprptr(1:2) = 1
   do i = 1, nsvar
      uprptr(i+2) = uprptr(i+1) + uprptr(i+2)
   end do

   ! Now iterate over lwr form, droppping entries into upr form
   allocate(uprrow(uprptr(nsvar+2)), stat=st)
   if(st.ne.0) return
   do csv = 1, nsvar
      do ii = ptr(csv), ptr(csv+1) - 1
         sv = row(ii)
         uprrow(uprptr(sv+1)) = csv
         uprptr(sv+1) = uprptr(sv+1) + 1
      end do
   end do

   ! Now determine supervariable elimination tree
   call etree_no_perm(nsvar, uprptr, uprrow, parent, st)
   if(st.ne.0) return
end subroutine mc78_elt_equiv_etree_integer

!
! Specialised version of mc78_etree that assumes elimination order is identity
!
subroutine etree_no_perm(n, ptr, row, parent, st)
   integer, intent(in) :: n ! dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! column pointers of A
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! row indices of A
   integer, dimension(n), intent(out) :: parent ! parent(i) is the
      ! parent of pivot i in the elimination tree
   integer, intent(out) :: st ! stat parmeter for allocate calls

   integer(pkg_type) :: ii ! next index into row
   integer :: k ! current ancestor
   integer :: l ! next ancestor
   integer :: piv ! current pivot
   integer, dimension(:), allocatable :: vforest ! virtual forest, used for
      ! path compression (shortcuts to top of each tree)

   ! Allocate virtual forest and initialise it
   allocate(vforest(n), stat=st)
   if(st.ne.0) return
   vforest(:) = n+1

   ! Loop over rows of A in pivot order
   do piv = 1, n
      ! Loop over entries in row in lower triangle of PAP^T
      do ii = ptr(piv), ptr(piv+1)-1
         k = row(ii)
         do while(vforest(k).lt.piv)
            l = vforest(k)
            vforest(k) = piv
            k = l
         end do
         if(vforest(k).eq.piv) cycle ! Already done from here, don't overwrite
         parent(k) = piv
         vforest(k) = piv
      end do
      parent(piv) = n + 1 ! set to be a root if not overwritten
   end do
end subroutine etree_no_perm

!
! This subroutine will postorder the elimination tree. That is to say it will
! reorder the nodes of the tree such that they are in depth-first search order.
!
! This is done by performing a depth-first search to identify mapping from the
! original pivot order to the new one. This map is then applied to order, invp
! and parent to enact the relabelling.
!
subroutine mc78_postorder_std(n, perm, invp, parent, st, block_pivots)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: perm ! perm(i) is the pivot
      ! position of column i
   integer, dimension(n), intent(inout) :: invp ! inverse of perm
   integer, dimension(n), intent(inout) :: parent ! parent(i) is the
      ! parent of pivot i in the elimination tree
   integer, intent(out) :: st ! stat parmeter for allocate calls
   integer, dimension(n), optional, intent(inout) :: block_pivots ! If
      ! present, then matches pivot order and specifies block pivots. 
      ! block_pivots(i) may take one of the following values:
      ! 0 - pivot i is in the middle of a block pivot
      ! 1 - pivot i is the first pivot of a block pivot
      ! 2 - pivot i is the last pivot of a block pivot
      ! 3 - pivot i is a 1x1 pivot

   integer, dimension(:), allocatable :: chead ! chead(i) is first child of i
   integer, dimension(:), allocatable :: cnext ! cnext(i) is next child of i
   integer :: i
   integer :: id
   integer :: j
   integer, dimension(:), allocatable :: map ! mapping from original pivot
      ! order to new one
   integer :: node
   integer :: shead ! pointer to top of stack
   integer, dimension(:), allocatable :: stack ! stack for depth first search

   !
   ! Build linked lists of children for each node
   !
   allocate(chead(n+1), cnext(n+1), stat=st)
   if(st.ne.0) return
   chead(:) = -1 ! no parent if necessary
   do i = n, 1, -1 ! do in reverse order so they come off in original order
      j = parent(i)
      cnext(i) = chead(j)
      chead(j) = i
   end do

   !
   ! Perform depth first search to build map
   !
   allocate(map(n+1), stack(n), stat=st)
   if(st.ne.0) return
   ! Place virtual root on top of stack
   shead = 1
   stack(shead) = n+1
   id = n + 1 ! next node id
   do while(shead.ne.0)
      ! Get node from top of stack
      node = stack(shead)
      shead = shead - 1

      ! Number it
      map(node) = id
      id = id - 1

      ! Place all its children on the stack such that the last child is
      ! at the top of the stack and first child closest to the bottom
      i = chead(node)
      do while(i.ne.-1)
         shead = shead + 1
         stack(shead) = i
         i = cnext(i)
      end do
   end do

   !
   ! Apply map to perm, invp and parent (and block_pivots if present)
   !

   ! invp is straight forward, use stack as a temporary
   stack(1:n) = invp(1:n)
   do i = 1, n
      j = map(i)
      invp(j) = stack(i)
   end do

   ! perm can be easily done as the inverse of invp
   do i = 1, n
      perm(invp(i)) = i
   end do

   ! parent is done in two stages. The first copies it to stack and permutes
   ! parent(i), but not the locations. i.e. if 1 is a parent of 3, and
   ! map(1)=2 and map(3)=4, then the first stage sets stack(1) = 4.
   ! The second stage then permutes the entries of map back into parent
   do i = 1, n
      stack(i) = map(parent(i))
   end do
   do i = 1, n
      parent(map(i)) = stack(i)
   end do

   ! permute block_pivots if required
   if(present(block_pivots)) then
      stack(1:n) = block_pivots(1:n)
      do i = 1, n
         block_pivots(map(i)) = stack(i)
      end do
   endif
end subroutine mc78_postorder_std

!
! This subroutine will postorder the elimination tree. That is to say it will
! reorder the nodes of the tree such that they are in depth-first search order.
!
! This is done by performing a depth-first search to identify mapping from the
! original pivot order to the new one. This map is then applied to order, invp
! and parent to enact the relabelling.
!
subroutine mc78_postorder_detect(n, realn, ptr, perm, invp, parent, st, &
      block_pivots)
   integer, intent(in) :: n
   integer, intent(out) :: realn
   integer(pkg_type), dimension(n+1), intent(in) :: ptr
   integer, dimension(n), intent(inout) :: perm ! perm(i) is the pivot
      ! position of column i
   integer, dimension(n), intent(inout) :: invp ! inverse of perm
   integer, dimension(n), intent(inout) :: parent ! parent(i) is the
      ! parent of pivot i in the elimination tree
   integer, intent(out) :: st ! stat parmeter for allocate calls
   integer, dimension(n), optional, intent(inout) :: block_pivots ! If
      ! present, then matches pivot order and specifies block pivots. 
      ! block_pivots(i) may take one of the following values:
      ! 0 - pivot i is in the middle of a block pivot
      ! 1 - pivot i is the first pivot of a block pivot
      ! 2 - pivot i is the last pivot of a block pivot
      ! 3 - pivot i is a 1x1 pivot

   integer, dimension(:), allocatable :: chead ! chead(i) is first child of i
   integer, dimension(:), allocatable :: cnext ! cnext(i) is next child of i
   integer :: i
   integer :: id
   integer :: j
   integer, dimension(:), allocatable :: map ! mapping from original pivot
      ! order to new one
   integer :: node
   integer :: shead ! pointer to top of stack
   integer, dimension(:), allocatable :: stack ! stack for depth first search

   realn = n

   !
   ! Build linked lists of children for each node
   !
   allocate(chead(n+1), cnext(n+1), stat=st)
   if(st.ne.0) return
   chead(:) = -1 ! no parent if necessary
   do i = n, 1, -1 ! do in reverse order so they come off in original order
      j = parent(i)
      cnext(i) = chead(j)
      chead(j) = i
   end do

   !
   ! Perform depth first search to build map
   !
   allocate(map(n+1), stack(n), stat=st)
   if(st.ne.0) return
   ! Place virtual root on top of stack
   shead = 1
   stack(shead) = n+1
   id = n + 1 ! next node id
   do while(shead.ne.0)
      ! Get node from top of stack
      node = stack(shead)
      shead = shead - 1

      ! Number it
      map(node) = id
      id = id - 1

      ! Place all its children on the stack such that the last child is
      ! at the top of the stack and first child closest to the bottom
      if(node.eq.n+1) then
         ! Virtual root node, detect children with no entries at same time
         ! placing those that are empty at the top of the stack
         ! First do those which are proper roots
         i = chead(node)
         do while(i.ne.-1)
            if(ptr(invp(i)+1)-ptr(invp(i)).eq.0) then
               i = cnext(i)
               cycle
            endif
            shead = shead + 1
            stack(shead) = i
            i = cnext(i)
         end do
         ! Second do those which are null roots
         i = chead(node)
         do while(i.ne.-1)
            if(ptr(invp(i)+1)-ptr(invp(i)).ne.0) then
               i = cnext(i)
               cycle
            endif
            realn = realn - 1
            shead = shead + 1
            stack(shead) = i
            i = cnext(i)
         end do
      else ! A normal node
         i = chead(node)
         do while(i.ne.-1)
            shead = shead + 1
            stack(shead) = i
            i = cnext(i)
         end do
      endif
   end do

   !
   ! Apply map to perm, invp and parent (and block_pivots if present)
   !

   ! invp is straight forward, use stack as a temporary
   stack(1:n) = invp(1:n)
   do i = 1, n
      j = map(i)
      invp(j) = stack(i)
   end do

   ! perm can be easily done as the inverse of invp
   do i = 1, n
      perm(invp(i)) = i
   end do

   ! parent is done in two stages. The first copies it to stack and permutes
   ! parent(i), but not the locations. i.e. if 1 is a parent of 3, and
   ! map(1)=2 and map(3)=4, then the first stage sets stack(1) = 4.
   ! The second stage then permutes the entries of map back into parent
   do i = 1, n
      stack(i) = map(parent(i))
   end do
   do i = 1, n
      parent(map(i)) = stack(i)
   end do

   ! permute block_pivots if required
   if(present(block_pivots)) then
      stack(1:n) = block_pivots(1:n)
      do i = 1, n
         block_pivots(map(i)) = stack(i)
      end do
   endif
end subroutine mc78_postorder_detect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Column count routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine determines column counts given the elimination tree and
! pattern of the matrix PAP^T.
!
! The algorithm is a specialisation of that given by Gilbert, Ng and Peyton [1],
! to only determine column counts. It is also described in Section 4.4 "Row
! counts" of [2].
!
! The essential technique is to determine the net number of entries introduced
! at a node (the "weight" in [1]). This is composed over the following terms:
!  wt[i] = [ - #children of node i
!            - #common indices between children
!            + #additional "new" row indices from column of A ]
!
! The clever part of this algorithm is how to determine the number of common
! indices between the children. This is accomplished by storing the last column
! at which an index was encountered, and a partial elimination tree. This
! partial elimination tree consists of all nodes processed so far, plus their
! parents. As we have a postorder on the tree, the current top of the tree
! containing node i is the least common ancestor of node i and the current node.
! We then observe that the first time an index will be double counted is at the
! least common ancestor of the current node and the last node where it was
! encountered.
!
! [1] Gilbert, Ng, Peyton, "An efficient algorithm to compute row and column
!     counts for sparse Cholesky factorization", SIMAX 15(4) 1994.
!
! [2] Tim Davis's book "Direct Methods for Sparse Linear Systems", SIAM 2006.
!
subroutine mc78_col_counts_integer(n, ptr, row, perm, invp, parent, cc, st, wt)
   integer, intent(in) :: n ! dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! column pointers of A
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! row indices of A
   integer, dimension(n), intent(in) :: perm ! perm(i) is the pivot
      ! position of column i
   integer, dimension(n), intent(in) :: invp ! inverse of perm
   integer, dimension(n), intent(in) :: parent ! parent(i) is the
      ! parent of pivot i in the elimination tree
   integer, dimension(n+1), intent(out) :: cc ! On exit, cc(i) is the
      ! number of entries in the lower triangular part of L (includes diagonal)
      ! for the column containing pivot i. For most of the routine however, it
      ! is used as a work space to track the net number of entries appearing
      ! for the first time at node i of the elimination tree (this may be
      ! negative).
   integer, intent(out) :: st ! stat parmeter for allocate calls
   integer, dimension(:), optional, intent(in) :: wt ! weights (eg number of
      ! variables in each supervariable)
   
   integer :: col ! column of matrix associated with piv
   integer, dimension(:), allocatable :: first ! first descendants
   integer :: i
   integer(pkg_type) :: ii
   integer :: totalwt
   integer, dimension(:), allocatable :: last_nbr ! previous neighbour
   integer, dimension(:), allocatable :: last_p ! previous p?
   integer :: par ! parent node of piv
   integer :: piv ! current pivot
   integer :: pp ! last pivot where u was encountered
   integer :: lca ! least common ancestor of piv and pp
   integer :: u ! current entry in column col
   integer :: uwt ! weight of u
   integer, dimension(:), allocatable :: vforest ! virtual forest

   !
   ! Determine first descendants, and set cc = 1 for leaves and cc = 0 for
   ! non-leaves.
   !
   allocate(first(n+1), stat=st)
   if(st.ne.0) return
   do i = 1, n+1
      first(i) = i
   end do
   if(present(wt)) then
      totalwt = 0 ! Find sum of weights so we can determine non-physical value
      do i = 1, n
         par = parent(i)
         first(par) = min(first(i), first(par)) ! first descendant
         if(first(i).eq.i) then ! is it a leaf or not?
            cc(i) = wt(invp(i))
         else
            cc(i) = 0
         endif
         totalwt = totalwt + wt(invp(i))
      end do
      cc(n+1) = totalwt + 1 ! Set to non-physical value
   else
      do i = 1, n
         par = parent(i)
         first(par) = min(first(i), first(par)) ! first descendant
         if(first(i).eq.i) then ! is it a leaf or not?
            cc(i) = 1
         else
            cc(i) = 0
         endif
      end do
      cc(n+1) = n + 1 ! Set to non-physical value
   endif

   !
   ! We store the partial elimination trees in a virtual forest. It is
   ! initialised such that each node is in its own tree to begin with.
   !
   allocate(vforest(n+1), stat=st)
   if(st.ne.0) return
   vforest(:) = 0

   !
   ! Initialise previous pivot and neightbour arrays to indicate no previous
   ! pivot or neightbour.
   !
   allocate(last_p(n+1), last_nbr(n+1), stat=st)
   if(st.ne.0) return
   last_p(:) = 0
   last_nbr(:) = 0

   !
   ! Determine cc(i), the number of net new entries to pass up tree from
   ! node i.
   !
   do piv = 1, n
      ! Loop over entries in column below the diagonal
      col = invp(piv)
      do ii = ptr(col), ptr(col+1)-1
         u = perm(row(ii))
         if(u.le.piv) cycle ! not in lower triangular part

         ! Check if entry has been seen by a descendant of this pivot, if
         ! so we skip the tests that would first add one to the current
         ! pivot's weight before then subtracting it again.
         if(first(piv).gt.last_nbr(u)) then
            ! Count new entry in current column
            uwt = 1
            if(present(wt)) uwt = wt(invp(u))
            cc(piv) = cc(piv) + uwt

            ! Determine least common ancestor of piv and the node at which
            ! u was last encountred
            pp = last_p(u)
            if(pp.ne.0) then
               ! u has been seen before, find top of partial elimination
               ! tree for node pp
               lca = FIND(vforest, pp)
               ! prevent double counting of u at node lca
               cc(lca) = cc(lca) - uwt
            endif

            ! Update last as u has now been seen at piv.
            last_p(u) = piv
         endif

         ! Record last neighbour of u so we can determine if it has been
         ! seen in this subtree before
         last_nbr(u) = piv
      end do
      ! Pass uneliminated variables up to parent
      par = parent(piv)
      if(present(wt)) then
         cc(par) = cc(par) + cc(piv) - wt(invp(piv))
      else
         cc(par) = cc(par) + cc(piv) - 1
      endif

      ! place the parent of piv into the same partial elimination tree as piv
      vforest(piv) = par ! operation "UNION" from [1]
   end do
end subroutine mc78_col_counts_integer

! Return top most element of tree containing u.
! Implements path compression to speed up subsequent searches.
integer function FIND(vforest, u)
   integer, dimension(:), intent(inout) :: vforest
   integer, intent(in) :: u

   integer :: current, prev

   prev = -1
   current = u
   do while(vforest(current).ne.0)
      prev = current
      current = vforest(current)
      if(vforest(current).ne.0) vforest(prev) = vforest(current)
   end do

   FIND = current
end function FIND

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Supernode amalgamation routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine identifies (relaxed) supernodes from the elimination tree
! and column counts.
!
! A node, u, and its parent, v, are merged if:
! (a) No new fill-in is introduced i.e. cc(v) = cc(u)-1
! (b) The number of columns in both u and v is less than nemin
!
! Note: assembly tree must be POSTORDERED on output
subroutine mc78_supernodes(n, realn, parent, cc, sperm, nnodes, sptr, sparent, &
      scc, invp, control, info, st, wt, block_pivots)
   integer, intent(in) :: n
   integer, intent(in) :: realn
   integer, dimension(n), intent(in) :: parent ! parent(i) is the
      ! parent of supernode i in the elimination/assembly tree. 
   integer, dimension(n), intent(in) :: cc ! cc(i) is the column count
      ! of supernode i, including elements eliminated at supernode i.
   integer, dimension(n), intent(out) :: sperm ! on exit contains a permutation
      ! from pivot order to a new pivot order with contigous supernodes
   integer, intent(out) :: nnodes ! number of supernodes
   integer, dimension(n+1), intent(out) :: sptr
   integer, dimension(n), intent(out) :: sparent
   integer, dimension(n), intent(out) :: scc
   integer, dimension(n), intent(in) :: invp
   type(mc78_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(out) :: st ! stat paremter from allocate calls
   integer, dimension(n), optional, intent(in) :: wt ! weights (number of vars
      ! in each initial node)
   integer, dimension(n), optional, intent(in) :: block_pivots ! If
      ! present, then matches pivot order and specifies block pivots. 
      ! block_pivots(i) may take one of the following values:
      ! 0 - pivot i is in the middle of a block pivot
      ! 1 - pivot i is the first pivot of a block pivot
      ! 2 - pivot i is the last pivot of a block pivot
      ! 3 - pivot i is a 1x1 pivot

   integer :: i, j, k
   integer :: flag
   integer, dimension(:), allocatable :: height ! used to track height of tree
   logical, dimension(:), allocatable :: mark ! flag array for nodes to finalise
   integer, dimension(:), allocatable :: map ! map vertex idx -> supernode idx
   integer, dimension(:), allocatable :: nelim ! number of eliminated variables
   integer, dimension(:), allocatable :: nvert ! number of elimd supervariables
   integer :: node
   integer, dimension(:), allocatable :: npar ! temporary array of snode pars
   integer :: par ! parent of current node
   integer :: shead ! current head of stack
   integer, dimension(:), allocatable :: stack ! used to navigate tree
   integer :: v
   integer, dimension(:), allocatable :: vhead ! heads of vertex linked lists
   integer, dimension(:), allocatable :: vnext ! next element in linked lists
   integer(long), dimension(:), allocatable :: ezero ! number of explicit zeros
   integer, dimension(:), allocatable :: chead ! chead(i) is first child of i
   integer, dimension(:), allocatable :: cnext ! cnext(i) is next child of i
   integer, dimension(:), allocatable :: child
   integer :: nchild
   integer :: start ! First pivot in block pivot
   integer :: totalwt ! sum of weights

   !
   ! Initialise supernode representation
   !
   allocate(nelim(n+1), nvert(n+1), vhead(n+1), vnext(n+1), stack(n), &
      height(n+1), mark(n), stat=st)
   if(st.ne.0) goto 490
   vnext(:) = -1
   vhead(:) = -1
   height(:) = 1

   ! Initialise number of variables in each node
   if(present(wt)) then
      totalwt = 0
      do i = 1, n
         nelim(i) = wt(invp(i))
         totalwt = totalwt + wt(invp(i))
      end do
   else ! All nodes initially contain a single variable
      nelim(:) = 1
      totalwt = n
   endif
   nvert(:) = 1

   allocate(map(n+1), npar(n+1), ezero(n+1), stat=st)
   if(st.ne.0) goto 490

   ezero(:) = 0 ! Initially no explicit zeros
   ezero(n+1) = huge(ezero) ! ensure root is not merged

   ! Ensure virtual root never gets amlgamated
   nelim(n+1) = totalwt+1 + control%nemin

   !
   ! Build child linked lists for nodes; merge block pivots if needed
   !
   allocate(chead(n+1), cnext(n+1), child(n), stat=st)
   if(st.ne.0) goto 490
   if(present(block_pivots)) then
      chead(:) = -1 ! no parent if necessary
      do i = realn, 1, -1 ! do in reverse order so come off in original order
         if(block_pivots(i).lt.2) cycle
         j = parent(i)
         if(j.ne.n+1) then
            do while(block_pivots(j).lt.2)
               j = parent(j)
            end do
         end if
         cnext(i) = chead(j)
         chead(j) = i
      end do
   else
      chead(:) = -1 ! no parent if necessary
      do i = realn, 1, -1 ! do in reverse order so come off in original order
         j = parent(i)
         cnext(i) = chead(j)
         chead(j) = i
      end do
   endif

   !
   ! Merge supernodes.
   !
   flag = 0
   v = 1
   nnodes = 0
   start=n+2
   do par = 1, n+1
      if(present(block_pivots) .and. par.lt.n+1) then
         if(block_pivots(par).lt.2) then
            if(start.ge.n+1) start = par
            cycle
         endif
         ! Merge pivots start to par, but don't add to vertex list (yet)

         do node = start, par-1
            ! Add together eliminated variables
            nelim(par) = nelim(par) + nelim(node)
            nvert(par) = nvert(par) + nvert(node)

            ! nodes have same height
            height(par) = max(height(par), height(node))
         end do
      endif

      nchild = 0
      node = chead(par)
      do while(node.ne.-1)
         nchild = nchild + 1
         child(nchild) = node
         node = cnext(node)
      end do
      call sort_by_val(nchild, child, cc, st)
      if(st.ne.0) goto 490

      do j = 1, nchild
         node = child(j)
         if(do_merge(node, par, nelim, cc, ezero, control, invp, flag, wt)) then
            ! Merge contents of node into par. Delete node.
            call merge_nodes(node, par, nelim, nvert, vhead, vnext, height, &
               ezero, cc)
            mark(node) = .false.
         else
            mark(node) = .true.
         endif
      end do

      if(present(block_pivots) .and. par.lt.n+1) then
         if(block_pivots(par).ge.2) then
            ! Add vertices start to par-1 into par
            do node = start, par-1
               vnext(node) = vhead(par)
               vhead(par) = node
               mark(node) = .false.
            end do
            start = n+2
         endif
      endif
   end do

   if(flag.ne.0) then
      if(control%unit_error.gt.0) write(control%unit_error, "(a)") &
         "MC78 Internal Error: Unrecognised amalgamation heuristic."
      info = MC78_ERROR_UNKNOWN
      return
   endif

   do node = 1, realn
      if(.not.mark(node)) cycle
      ! Node not merged, now a complete supernode

      ! Record start of supernode
      nnodes = nnodes + 1
      sptr(nnodes) = v
      npar(nnodes) = parent(node)
      if(present(wt)) then
         scc(nnodes) = cc(node) + nelim(node) - wt(invp(node))
      else
         scc(nnodes) = cc(node) + nelim(node) - 1
      endif

      ! Record height in tree of parent vertices
      height(parent(node)) = max(height(parent(node)), height(node) + 1)

      ! Determine last vertex of node so we can number backwards
      v = v + nvert(node)
      k = v

      ! Loop over member vertices of node and number them
      shead = 1
      stack(shead) = node
      do while(shead.gt.0)
         i = stack(shead)
         shead = shead - 1

         ! Order current vertex
         k = k - 1
         sperm(i) = k
         map(i) = nnodes

         ! Stack successor, if any
         if(vnext(i).ne.-1) then
            shead = shead + 1
            stack(shead) = vnext(i)
         endif

         ! Descend into tree rooted at i
         if(vhead(i).ne.-1) then
            shead = shead + 1
            stack(shead) = vhead(i)
         endif
      end do
   end do
   sptr(nnodes+1) = v ! Record end of final supernode
   map(n+1) = nnodes + 1 ! virtual root vertex maps to virtual root sn
   npar(nnodes+1) = n + 1

   ! Handle permutation of empty columns
   do i = realn+1, n
      sperm(i) = i
   end do

   ! Allocate arrays for return and copy data into them correctly
   do node = 1, nnodes
      par = npar(node) ! parent /vertex/ of supernode
      par = map(par)   ! parent /node/   of supernode
      sparent(node) = par ! store parent
   end do

   return

   490 continue
   info = MC78_ERROR_ALLOC
   return
end subroutine mc78_supernodes

!
! Sort n items labelled by idx into decreasing order of val(idx(i))
!
recursive subroutine sort_by_val(n, idx, val, st)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: idx
   integer, dimension(:), intent(in) :: val
   integer, intent(out) :: st

   integer :: ice_idx, ice_val, ik_idx, ik_val
   integer :: klo,kor,k,kdummy

   st = 0

   if(n.ge.minsz_ms) then
      call sort_by_val_ms(n, idx, val, st)
   else
      klo = 2
      kor = n
      do kdummy = klo,n
         ! items kor, kor+1, .... ,nchild are in order
         ice_idx = idx(kor-1)
         ice_val = val(ice_idx)
         do k = kor,n
            ik_idx = idx(k)
            ik_val = val(ik_idx)
            if (ice_val >= ik_val) exit
            idx(k-1) = ik_idx
         end do
         idx(k-1) = ice_idx
         kor = kor - 1
      end do
   endif
end subroutine sort_by_val

! Sort n items labelled by idx into decreasing order of val(idx(i))
!
! Merge sort version, dramatically improves performance for nodes with large
! numbers of children
! (Passes to simple sort for small numbers of entries)
recursive subroutine sort_by_val_ms(n, idx, val, st)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: idx
   integer, dimension(:), intent(in) :: val
   integer, intent(out) :: st

   integer :: i, j, jj, jj2, k, kk, kk2
   integer :: mid
   integer, dimension(:), allocatable :: work

   if(n.le.1) return
   if(n.lt.minsz_ms) then
      call sort_by_val(n, idx, val, st)
      return
   endif
   mid = (n-1)/2 + 1

   ! Recurse to order half lists
   call sort_by_val_ms(mid, idx(1:mid), val, st)
   if(st.ne.0) return
   call sort_by_val_ms(n - mid, idx(mid+1:n), val, st)
   if(st.ne.0) return

   ! Merge two half lists
   ! (Take a copy of the first half list so we don't overwrite it)
   allocate(work(mid), stat=st)
   if(st.ne.0) return
   work(:) = idx(1:mid)
   j = 1
   k = mid+1
   jj = work(j)
   jj2 = val(jj)
   kk = idx(k)
   kk2 = val(kk)
   do i = 1, n
      if(jj2.ge.kk2) then
         idx(i) = jj
         j = j + 1
         if(j.gt.mid) exit
         jj = work(j)
         jj2 = val(jj)
      else
         idx(i) = kk
         k = k + 1
         if(k.gt.n) exit
         kk = idx(k)
         kk2 = val(kk)
      endif
   end do
   if(j.le.mid) idx(i+1:n) = work(j:mid)
end subroutine sort_by_val_ms

!
! Return .true. if we should merge node and par, .false. if we should not
!
logical function do_merge(node, par, nelim, cc, ezero, control, invp, info, wt)
   integer, intent(in) :: node ! node to merge and delete
   integer, intent(in) :: par ! parent to merge into
   integer, dimension(:), intent(in) :: nelim
   integer, dimension(:), intent(in) :: cc
   integer(long), dimension(:), intent(in) :: ezero
   type(mc78_control), intent(in) :: control
   integer, dimension(:), intent(in) :: invp
   integer, intent(out) :: info
   integer, dimension(:), optional, intent(in) :: wt

   real(dp) :: z, ne

   info = 0

   if(ezero(par).eq.huge(ezero)) then
      do_merge = .false.
      return
   endif

   select case(control%heuristic)
   case(1)
      !
      ! HSL_MA77 style nemin
      !
      if(present(wt)) then
         do_merge = (cc(par).eq.cc(node)-wt(invp(node)) .and. &
            nelim(par).eq.wt(invp(par))) .or. &
            (nelim(par).lt.control%nemin .and. nelim(node).lt.control%nemin)
      else
         do_merge = (cc(par).eq.cc(node)-1 .and. nelim(par).eq.1) .or. &
            (nelim(par).lt.control%nemin .and. nelim(node).lt.control%nemin)
      endif
   case(2)
      !
      ! CHOLMOD style nrelax/zrelax
      !
      ! FIXME: currently assumes nodes are square, not trapezoidal

      ! calculate number of non-zeros in new node
      z = ezero(par) + ezero(node) + &
         (cc(par)-1+nelim(par) - cc(node)+1) * nelim(par)
      ! find this as a fraction of total non-zeros in new node
      ne = nelim(par) + nelim(node)
      z = z / ( (cc(par)-1+ne)*ne )

      do_merge = (ne .le. control%nrelax(1)) .or. &
         (cc(par).eq.cc(node)-1 .and. nelim(par).eq.1) .or. &
         (ne .le. control%nrelax(2) .and. z .lt. control%zrelax(1)) .or. &
         (ne .le. control%nrelax(3) .and. z .lt. control%zrelax(2)) .or. &
         (z .lt. control%zrelax(3))
   case default
      ! Note: This bit of code should NEVER execute
      do_merge = .false.
      info = MC78_ERROR_UNKNOWN
   end select
end function do_merge

!
! This subroutine merges node with its parent, deleting node in the process.
!
subroutine merge_nodes(node, par, nelim, nvert, vhead, vnext, height, ezero, cc)
   integer, intent(in) :: node ! node to merge and delete
   integer, intent(in) :: par ! parent to merge into
   integer, dimension(:), intent(inout) :: nelim
   integer, dimension(:), intent(inout) :: nvert
   integer, dimension(:), intent(inout) :: vhead
   integer, dimension(:), intent(inout) :: vnext
   integer, dimension(:), intent(inout) :: height
   integer(long), dimension(:), intent(inout) :: ezero
   integer, dimension(:), intent(in) :: cc

   ! Add node to list of children merged into par
   vnext(node) = vhead(par)
   vhead(par) = node

   ! Work out number of explicit zeros in new node
   ! FIXME: probably wrong now with weights and block pivots
   ezero(par) = ezero(par) + ezero(node) + &
      (cc(par)-1+nelim(par) - cc(node) + 1_long) * nelim(par)

   ! Add together eliminated variables
   nelim(par) = nelim(par) + nelim(node)
   nvert(par) = nvert(par) + nvert(node)

   ! nodes have same height
   height(par) = max(height(par), height(node))
end subroutine merge_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Statistics routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine merely calculates interesting statistics
!
subroutine mc78_stats(nnodes, sptr, scc, nfact, nflops)
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer, dimension(nnodes), intent(in) :: scc
   integer(long), optional, intent(out) :: nfact
   integer(long), optional, intent(out) :: nflops

   integer :: j
   integer :: m ! number of entries in retangular part of ndoe
   integer :: nelim ! width of node
   integer :: node ! current node of assembly tree
   integer(long) :: r_nfact, r_nflops

   if(.not.present(nfact) .and. .not.present(nflops)) return ! nothing to do

   r_nfact = 0
   r_nflops = 0
   do node = 1, nnodes
      nelim = sptr(node+1) - sptr(node)
      m = scc(node) - nelim

      ! number of entries
      r_nfact = r_nfact + (nelim * (nelim+1)) / 2 ! triangular block
      r_nfact = r_nfact + nelim * m ! below triangular block

      ! flops
      do j = 1, nelim
         r_nflops = r_nflops + (m+j)**2
      end do
   end do

   if(present(nfact)) nfact = r_nfact
   if(present(nflops)) nflops = r_nflops

   !print *, "n = ", n
   !print *, "nnodes = ", nnodes
   !print *, "nfact = ", nfact
   !print *, "sum cc=", sum_long(cc(1:n))
   !print *, "nflops = ", nflops
end subroutine mc78_stats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Row list routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine determines the row indices for each supernode
!
subroutine mc78_row_lists_nosvar_integer(n, ptr, row, perm, invp, nnodes, &
      sptr, sparent, scc, rptr, rlist, control, info, st)
   integer, intent(in) :: n
   integer(pkg_type), dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   integer, dimension(n), intent(in) :: perm
   integer, dimension(n), intent(in) :: invp
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer, dimension(nnodes), intent(in) :: sparent
   integer, dimension(nnodes), intent(in) :: scc
   integer(long), dimension(nnodes+1), intent(out) :: rptr
   integer, dimension(sum_long(scc(1:nnodes))), intent(out) :: rlist
   type(mc78_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(out) :: st

   integer :: child ! current child of node
   integer :: col ! current column of matrix corresponding to piv
   integer(long) :: i
   integer(long) :: idx ! current insert position into nodes(node)%index
   integer :: j
   integer :: node ! current node of assembly tree
   integer :: piv ! current pivot position
   integer, dimension(:), allocatable :: seen ! tag last time index was seen
   integer, dimension(:), allocatable :: chead ! head of child linked lists
   integer, dimension(:), allocatable :: cnext ! pointer to next child

   ! Allocate and initialise memory
   allocate(seen(n), chead(nnodes+1), cnext(nnodes+1), stat=st)
   if(st.ne.0) then
      info = MC78_ERROR_ALLOC
      return
   endif
   seen(:) = 0
   chead(:) = -1

   ! Build child linked lists (backwards so pop off in good order)
   do node = nnodes, 1, -1
      i = sparent(node)
      cnext(node) = chead(i)
      chead(i) = node
   end do

   ! Loop over nodes from bottom up building row lists.
   rptr(1) = 1
   do node = 1, nnodes

      ! Allocate space for row indices
      rptr(node+1) = rptr(node) + scc(node)
      idx = rptr(node) ! insert position

      ! Add entries eliminated at this node
      do piv = sptr(node), sptr(node+1)-1
         seen(piv) = node
         rlist(idx) = piv
         idx = idx + 1
      end do

      ! Find indices inherited from children
      child = chead(node)
      do while (child.ne.-1)
         do i = rptr(child), rptr(child+1)-1
            j = rlist(i)
            if(j.lt.sptr(node)) cycle ! eliminated
            if(seen(j).eq.node) cycle ! already seen
            seen(j) = node
            rlist(idx) = j
            idx = idx + 1
         end do
         child = cnext(child)
      end do

      ! Find new indices from A
      do piv = sptr(node), sptr(node+1)-1
         col = invp(piv)
         do i = ptr(col), ptr(col+1)-1
            j = perm(row(i))
            if(j.lt.piv) cycle ! in upper triangle
            if(seen(j).eq.node) cycle ! already seen in this snode
            ! Otherwise, this is a new entry
            seen(j) = node
            rlist(idx) = j
            idx = idx + 1
         end do
      end do

      ! Note: following error check won't work with block pivots
      !if(idx .ne. rptr(node+1)) then
      !   ! Note: This bit of code should NEVER execute
      !  if(control%unit_error.gt.0) write(control%unit_error, "(3(a,i8))") &
      !      "MC78 Internal Error: node ", node, ": found ", idx-rptr(node), &
      !      " entries, but expected to find ", rptr(node+1)-rptr(node)
      !   info = MC78_ERROR_UNKNOWN
      !   !print *, rlist(1:idx-1)
      !   return
      !endif
   end do
end subroutine mc78_row_lists_nosvar_integer

subroutine mc78_row_lists_svar_integer(nsvar, svar, n, ptr, row, perm, invp, &
      nnodes, sptr, sparent, scc, rptr, rlist, control, info, st)
   integer, intent(in) :: nsvar
   integer, dimension(nsvar), intent(in) :: svar
   integer, intent(in) :: n
   integer(pkg_type), dimension(nsvar+1), intent(in) :: ptr
   integer, dimension(ptr(nsvar+1)-1), intent(in) :: row
   integer, dimension(nsvar), intent(in) :: perm
   integer, dimension(nsvar), intent(in) :: invp
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer, dimension(nnodes), intent(in) :: sparent
   integer, dimension(nnodes), intent(in) :: scc
   integer(long), dimension(nnodes+1), intent(out) :: rptr
   integer, dimension(sum_long(scc(1:nnodes))), intent(out) :: rlist
   type(mc78_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(out) :: st

   integer :: child ! current child of node
   integer :: col ! current column of matrix corresponding to piv
   integer(long) :: i
   integer(long) :: idx ! current insert position into nodes(node)%index
   integer :: j
   integer :: k
   integer :: node ! current node of assembly tree
   integer :: piv ! current pivot position
   integer, dimension(:), allocatable :: seen ! tag last time index was seen
   integer, dimension(:), allocatable :: chead ! head of child linked lists
   integer, dimension(:), allocatable :: cnext ! pointer to next child
   integer, dimension(:), allocatable :: svptr ! pointers for row list starts

   ! Allocate and initialise memory
   allocate(seen(n), chead(nnodes+1), cnext(nnodes+1), svptr(nsvar+1), stat=st)
   if(st.ne.0) return
   seen(:) = 0
   chead(:) = -1

   ! Build svptr array
   svptr(1) = 1
   do i = 1, nsvar
      svptr(i+1) = svptr(i) + svar(invp(i))
   end do

   ! Build child linked lists (backwards so pop off in good order)
   do node = nnodes, 1, -1
      i = sparent(node)
      cnext(node) = chead(i)
      chead(i) = node
   end do

   ! Loop over nodes from bottom up building row lists.
   rptr(1) = 1
   do node = 1, nnodes

      ! Allocate space for row indices
      rptr(node+1) = rptr(node) + scc(node)
      idx = rptr(node) ! insert position

      ! Add entries eliminated at this node
      do i = sptr(node), sptr(node+1)-1
         do piv = svptr(i), svptr(i+1)-1
            seen(piv) = nnodes+1
            rlist(idx) = piv
            idx = idx + 1
         end do
      end do

      ! Find indices inherited from children
      child = chead(node)
      do while (child.ne.-1)
         do i = rptr(child), rptr(child+1)-1
            j = rlist(i)
            if(seen(j).ge.node) cycle ! already seen (or eliminated already)
            seen(j) = node
            rlist(idx) = j
            idx = idx + 1
         end do
         child = cnext(child)
      end do

      ! Find new indices from A
      do piv = sptr(node), sptr(node+1)-1
         col = invp(piv)
         do i = ptr(col), ptr(col+1)-1
            j = perm(row(i))
            if(seen(svptr(j)).ge.node) cycle ! already seen (or eliminated)
            ! Otherwise, this is a new entry
            ! Iterate over variables in supervariable
            do k = svptr(j), svptr(j+1)-1
               seen(k) = node
               rlist(idx) = k
               idx = idx + 1
            end do
         end do
      end do

      if(idx .ne. rptr(node+1)) then
         ! Note: This bit of code should NEVER execute
         if(control%unit_error.gt.0) write(control%unit_error, "(3(a,i8))") &
            "MC78 Internal Error: node ", node, ": found ", idx-rptr(node), &
            " entries, but expected to find ", rptr(node+1)-rptr(node)
         info = MC78_ERROR_UNKNOWN
         !print *, rlist(1:idx-1)
         return
      endif
   end do
end subroutine mc78_row_lists_svar_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Optimize cache locality routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! The following subroutine reorders the elimination order within each node in
! such a way that the order of the variables in
! the "primary" child has ordering agreement with the node. Consider
! the following tree:
!          A
!         / \
!        B   C
! 
! and assume B has more partially summed variables than C. Then
! B is the primary child and the ordering of the
! corresponding fully summed variables in the parent A matches the ordering
! of the partially summed variables in B (but not in C). Any additional
! partially summed variables present in C but not in B are then ordered in A
! such that they match C.
!
! This is done by two passes of the tree.
! The first builds a map from variables to the nodes at which
! they are eliminated, and orders the children of each node such
! that the first has the largest number of partially summed variables.
! The second pass uses a depth first search of the now ordered tree. It loops
! over non-fully summed variables and when it first encounters each it will
! place it as the next variable at its elimination node.
!
subroutine mc78_optimize_locality(n, realn, perm, invp, nnodes, sptr, sparent, &
      rptr, rlist, st, sort)

   integer, intent(in) :: n ! dimension of system
   integer, intent(in) :: realn ! symbolic dimension of system
   integer, dimension(n), intent(inout) :: perm ! on exit, will have been
      ! reordered for better cache locality
   integer, dimension(n), intent(inout) :: invp ! inverse of perm. on exit
      ! will have been changed to match new perm
   integer, intent(in) :: nnodes ! number of supernodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer, dimension(nnodes), intent(in) :: sparent
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, dimension(rptr(nnodes+1)-1), intent(inout) :: rlist
   integer, intent(out) :: st ! stat parameter
   logical, optional, intent(in) :: sort

   integer :: i
   integer(long) :: ii ! loop index
   integer :: id ! current insert position into list
   integer :: k
   integer :: j ! temporary variable
   integer, allocatable :: list(:) ! list of nodes in a weighted depth-first
      ! order such that children are visitied in decreasing number of
      ! partially summed variables (ie child with most p.s.v. visited first)
   integer, allocatable :: map(:) ! maps variables to nodes where
      ! they are eliminated
   integer :: node ! current node
   integer, allocatable :: ord(:) ! tracks number of variables ordered at
      ! each node
   integer, allocatable :: perm2(:) ! permutation to apply to perm
   integer :: pnode ! parent node
   integer :: shead ! current top of stack
   integer, allocatable :: stack(:) ! used for depth first walk of tree
   integer :: start ! first entry on stack of child from current node
   integer, allocatable :: chead(:) ! heads of child linked lists
   integer, allocatable :: cnext(:) ! tails of child linked lists

   ! Allocate arrays for depth first search of tree
   allocate (map(n), list(nnodes+1), stack(nnodes), chead(nnodes+1), &
      cnext(nnodes), stat=st)
   if (st /= 0) return

   !
   ! Build elimination map
   !
   do node = 1, nnodes
      do ii = sptr(node), sptr(node+1)-1
         map(ii) = node
      end do
   end do

   !
   ! Build child linked lists
   !
   chead(:) = -1 ! no child if necessary
   do i = nnodes, 1, -1 ! do in reverse order so they come off in original order
      j = sparent(i)
      cnext(i) = chead(j)
      chead(j) = i
   end do

   !
   ! Perform depth first search of tree, such that children of a node are
   ! visited in order of the number of partially summer variables, largest
   ! first.
   !
   shead = 1
   stack(shead) = nnodes + 1
   id = nnodes+1
   do while(shead.ne.0)
      ! Get node from top of stack
      node = stack(shead)
      shead = shead - 1

      ! Number it
      list(id) = node
      id = id - 1

      ! Place all its children on the stack
      start = shead + 1
      i = chead(node)
      do while(i.ne.-1)
         shead = shead + 1
         stack(shead) = i
         i = cnext(i)
      end do
      ! Order children just placed on stack such that child with least partially
      ! summed variables is at the top
      call order_children(shead-start+1, stack(start:shead), nnodes, sptr, &
         rptr, st)
      if(st.ne.0) return
   end do

   !
   ! Next loop over children reordering partially summed variables.
   !
   allocate(ord(nnodes),perm2(n),stat=st)
   if (st.ne.0) return

   do node = 1, nnodes
      ord(node) = sptr(node)
   end do

   do k = 1, nnodes
      node = list(k)

      ! Order variables first encountered at this node
      do ii = rptr(node), rptr(node+1)-1
         j = rlist(ii)
         pnode = map(j)
         if(pnode .ne. -1) then ! check if we have ordered j already
            ! order at parent
            perm2(j) = ord(pnode)
            ord(pnode) = ord(pnode) + 1
            map(j) = -1 ! mark as ordered
         endif
         rlist(ii) = perm2(j)
      end do
   end do

   do i = realn+1, n
      perm2(i) = i
   end do

   !
   ! Apply permutation to perm and invp
   !
   ! Use perm as a temporary variable to permute invp.
   perm(1:n) = invp(1:n)
   do i = 1, n
      j = perm2(i)
      invp(j) = perm(i)
   end do

   ! Recover invp as inverse of perm
   do i = 1, n
      perm(invp(i)) = i
   end do

   if(present(sort)) then
      if(sort) then
         call dbl_tr_sort(n, nnodes, rptr, rlist, st)
         if(st.ne.0) return
      endif
   endif
end subroutine mc78_optimize_locality

! Orders nodes stored in child(1:n) such that the number of partially summed
! variables at each node is decreasing (ie one with least is in posn n)
!
! Simple sort version, good for nodes with small numbers of children
! (Passes to mergesort for large numbers of entries)
recursive subroutine order_children(n, child, nnodes, sptr, rptr, st)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: child
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, intent(out) :: st

   integer :: nelim, m
   integer :: k, kdummy, klo, kor
   integer :: ice_idx, ice_psum, ik_idx, ik_psum

   st = 0

   if(n.ge.minsz_ms) then
      call order_children_ms(n, child, nnodes, sptr, rptr, st)
   else
      klo = 2
      kor = n
      do kdummy = klo, n
         ! items kor, kor+1, .... ,n are in order
         ice_idx = child(kor-1)
         nelim = sptr(ice_idx+1) - sptr(ice_idx)
         m = int(rptr(ice_idx+1) - rptr(ice_idx))
         ice_psum = m - nelim
         do k = kor, n
            ik_idx = child(k)
            nelim = sptr(ik_idx+1) - sptr(ik_idx)
            m = int(rptr(ik_idx+1) - rptr(ik_idx))
            ik_psum = m - nelim
            if (ice_psum .ge. ik_psum) exit
            child(k-1) = ik_idx
         end do
         child(k-1) = ice_idx
         kor = kor - 1
      end do
   endif
end subroutine order_children

! Orders nodes stored in child(1:n) such that the number of partially summed
! variables at each node is decreasing (ie one with least is in posn n)
!
! Merge sort version, dramatically improves performance for nodes with large
! numbers of children
! (Passes to simple sort for small numbers of entries)
recursive subroutine order_children_ms(n, child, nnodes, sptr, rptr, st)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: child
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, intent(out) :: st

   integer :: i, j, jj, jj2, k, kk, kk2, m, nelim
   integer :: mid
   integer, dimension(:), allocatable :: work

   if(n.le.1) return
   if(n.lt.minsz_ms) then
      call order_children(n, child, nnodes, sptr, rptr, st)
      return
   endif
   mid = (n-1)/2 + 1

   ! Recurse to order half lists
   call order_children_ms(mid, child(1:mid), nnodes, sptr, rptr, st)
   if(st.ne.0) return
   call order_children_ms(n - mid, child(mid+1:n), nnodes, sptr, rptr, st)
   if(st.ne.0) return

   ! Merge two half lists
   ! (Take a copy of the first half list so we don't overwrite it)
   allocate(work(mid), stat=st)
   if(st.ne.0) return
   work(:) = child(1:mid)
   j = 1
   k = mid+1
   jj = work(j)
   nelim = sptr(jj+1) - sptr(jj)
   m = int(rptr(jj+1) - rptr(jj))
   jj2 = m - nelim
   kk = child(k)
   nelim = sptr(kk+1) - sptr(kk)
   m = int(rptr(kk+1) - rptr(kk))
   kk2 = m - nelim
   do i = 1, n
      if(jj2.ge.kk2) then
         child(i) = jj
         j = j + 1
         if(j.gt.mid) exit
         jj = work(j)
         nelim = sptr(jj+1) - sptr(jj)
         m = int(rptr(jj+1) - rptr(jj))
         jj2 = m - nelim
      else
         child(i) = kk
         k = k + 1
         if(k.gt.n) exit
         kk = child(k)
         nelim = sptr(kk+1) - sptr(kk)
         m = int(rptr(kk+1) - rptr(kk))
         kk2 = m - nelim
      endif
   end do
   if(j.le.mid) child(i+1:n) = work(j:mid)
end subroutine order_children_ms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assorted auxilary routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Converts piv_size to block_pivots:
! piv_size(i) is size of block pivot containing column i of A
! block_pivots(j) is a flag for pivot j of L and has one of the following values
!  0 - pivot i is in the middle of a block pivot
!  1 - pivot i is the first pivot of a block pivot
!  2 - pivot i is the last pivot of a block pivot
!  3 - pivot i is a 1x1 pivot
!
subroutine convert_to_blk_piv(n, invp, block)
   integer, intent(in) :: n
   integer, dimension(n), intent(in) :: invp
   integer, dimension(n), intent(inout) :: block

   integer, dimension(:), allocatable :: blk2
   integer :: i, cnt

   allocate(blk2(n))

   ! Take a copy of block so we can change it
   blk2(1:n) = block(1:n)

   ! Iterate over pivots in elimination order, recording starts and ends of blks
   cnt = blk2(invp(1))-1 ! Initialise for first pivot
   block(1) = 1 ! First pivot is start of a block
   do i = 2, n
      block(i) = 0
      if(cnt.eq.0) then
         ! this is first pivot of a block, previous is last pivot of a block
         cnt = blk2(invp(i))
         block(i-1) = block(i-1) + 2
         block(i) = block(i) + 1
      endif
      cnt = cnt - 1
   end do
   block(n) = block(n) + 2 ! end of matrix must end a block pivot

end subroutine convert_to_blk_piv

!
! Converts block_pivots back to piv_size:
! block_pivots(j) is a flag for pivot j of L and has one of the following values
!  0 - pivot i is in the middle of a block pivot
!  1 - pivot i is the first pivot of a block pivot
!  2 - pivot i is the last pivot of a block pivot
!  3 - pivot i is a 1x1 pivot
! piv_size(i) is size of block pivot containing column i of A
!
subroutine convert_from_blk_piv(n, perm, block)
   integer, intent(in) :: n
   integer, dimension(n), intent(in) :: perm
   integer, dimension(n), intent(inout) :: block

   integer, dimension(:), allocatable :: blk2
   integer :: i, sa, cnt

   allocate(blk2(n))

   ! convert first/last notation to block size notation
   cnt = -1; sa = -1 ! these values should never actually be used
   do i = 1, n
      select case(block(i))
      case (0) ! middle pivot of a block
         cnt  = cnt + 1
      case (1) ! first pivot of a block
         sa = i
         cnt = 1
      case (2) ! end pivot of a block
         cnt = cnt + 1
         block(sa:i) = cnt
      case (3) ! only pivot of a block
         block(i) = 1
      end select
   end do

   ! Permute back to original matrix order
   blk2(1:n) = block(1:n)
   do i = 1, n
      block(i) = blk2(perm(i))
   end do
end subroutine convert_from_blk_piv

!
! This subroutine copies a matrix pattern and adds subdiagonal entries as
! needed to force block pivots to have a parent-child relation in the
! elimination tree.
!
subroutine mc78_block_prep(n, ptr, row, bptr, brow, perm, invp, block_pivots, &
      st)
   integer, intent(in) :: n ! Dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! Column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! Row indices
   integer(pkg_type), dimension(n+1), intent(out) :: bptr ! Column pointers
   integer, dimension(:), intent(out) :: brow ! Row indices
   integer, dimension(n), intent(inout) :: perm
      ! perm(i) must hold position of i in the pivot sequence. 
      ! On exit, holds the pivot order to be used by factorization.
   integer, dimension(n), intent(inout) :: invp ! inverse permutation of perm
   integer, dimension(n), intent(inout) :: block_pivots ! Matches pivot order
      ! and specifies block pivots.
      ! block_pivots(i) may take one of the following values:
      ! 0 - pivot i is in the middle of a block pivot
      ! 1 - pivot i is the first pivot of a block pivot
      ! 2 - pivot i is the last pivot of a block pivot
      ! 3 - pivot i is a 1x1 pivot
   integer, intent(out) :: st

   integer :: col
   integer :: i, j, k
   integer(pkg_type) :: ii
   integer :: idx
   integer :: piv
   integer, dimension(:), allocatable :: seen

   allocate(seen(n), stat=st)
   if(st.ne.0) return

   ! First pass through ptr and ensure that all block pivots contain no
   ! empty columns
   perm(:) = invp(:)
   piv = 1
   j = 1
   ! Handle variables that are actually used
   do while(piv.le.n)
      do i = piv, n
         if(block_pivots(i).ge.2) exit ! end of block pivot
      end do

      k = 0
      do piv = piv, i
         if(ptr(perm(piv)).eq.ptr(perm(piv)+1)) cycle
         invp(j) = perm(piv)
         j = j + 1
         if(k.eq.0) then
            ! This is the new start of the block pivot
            select case(block_pivots(piv))
            case(0) ! was in the middle. now a start
               block_pivots(piv) = 1
            case(2) ! was the end. now a 1x1
               block_pivots(piv) = 3
            end select
         endif
         k = piv
      end do
      if(k.ne.0) then
         ! The was at least one used variable in the block pivot
         select case(block_pivots(k))
         case(0) ! was the middle. now an end
            block_pivots(k) = 2
         case(1) ! was the start. now a 1x1
            block_pivots(k) = 3
         end select
      endif
      piv = i + 1
   end do
   ! Handle unused variables
   do piv = 1, n
      i = perm(piv)
      if(ptr(i).eq.ptr(i+1)) then
         invp(j) = i
         j = j + 1
         block_pivots(piv) = 3 ! Force to 1x1
      endif
   end do
   ! Map block_pivots in original variable order into sv_map
   do i = 1, n
      seen(perm(i)) = block_pivots(i)
   end do
   ! Map sv_map in new pivot order back into block_pivots
   do i = 1, n
      block_pivots(i) = seen(invp(i))
   end do
   ! Reestablish perm
   do i = 1, n
      perm(invp(i)) = i
   end do

   ! Now iterate over cleaned up block pivot sequence
   seen(:) = 0
   piv = 1
   idx = 1
   do col = 1, n
      piv = perm(col)
      bptr(col) = idx
      if(block_pivots(piv).eq.3) then
         ! 1x1 pivot, just copy the column
         idx = idx + int(ptr(col+1) - ptr(col))
         brow(bptr(col):idx-1) = row(ptr(col):ptr(col+1)-1)
      else
         ! copy the column, but add an entry on subdiagonal(s)
         do ii = ptr(col), ptr(col+1)-1
            j = row(ii)
            seen(j) = col
            brow(idx) = j
            idx = idx + 1
         end do
         if(block_pivots(piv).ne.1) then
            ! Not the first column, add an entry above the diagonal
            j = invp(piv-1)
            if(seen(j).lt.col) then
               brow(idx) = j
               idx = idx + 1
            endif
         endif
         if(block_pivots(piv).ne.2) then
            ! Not the last column, add an entry below the diagonal
            j = invp(piv+1)
            if(seen(j).lt.col) then
               brow(idx) = j
               idx = idx + 1
            endif
         endif
      endif
   end do
   bptr(n+1) = idx
end subroutine mc78_block_prep

!
! This subroutine will take information concerning a compressed matrix and a
! supervariable map, and will decompress the information so it relates to
! the original matrix
!
subroutine svar_unmap(n, nsvar, svar, perm, invp, nnodes, sinvp, &
      snptr, st)
   integer, intent(in) :: n
   integer, intent(in) :: nsvar
   integer, dimension(nsvar), intent(in) :: svar
   integer, dimension(n), intent(out) :: perm
   integer, dimension(n), intent(inout) :: invp
   integer, intent(in) :: nnodes
   integer, dimension(nsvar), intent(in) :: sinvp
   integer, dimension(nnodes+1), intent(inout) :: snptr
   integer, intent(out) :: st

   integer, dimension(:), allocatable :: svptr
   integer :: i, j, k
   integer :: j1, j2
   integer :: idx

   ! Set up svptr
   allocate(svptr(nsvar+1), stat=st)
   if(st.ne.0) return
   svptr(1) = 1
   do i = 1, nsvar
      svptr(i+1) = svptr(i) + svar(i)
   end do

   ! Take a copy of invp in perm to ease remapping
   perm(:) = invp(:)

   ! Remap invp
   idx = 1
   do i = 1, nsvar
      j = sinvp(i)
      do k = svptr(j), svptr(j+1)-1
         invp(idx) = perm(k)
         idx = idx + 1
      end do
   end do

   ! Expand supernode pointer
   j1 = snptr(1)
   do i = 1, nnodes
      j2 = snptr(i+1)
      snptr(i+1) = snptr(i)
      do j = j1, j2-1
         snptr(i+1) = snptr(i+1) + svar(sinvp(j))
      end do
      j1 = j2
   end do

   ! Finally, recover perm as inverse of invp
   do i = 1, n
      perm(invp(i)) = i
   end do
end subroutine svar_unmap


!
! This subroutine performs a double transpose sort on the row indices of sn
!
subroutine dbl_tr_sort(n, nnodes, rptr, rlist, st)
   integer, intent(in) :: n
   integer, intent(in) :: nnodes
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, dimension(rptr(nnodes+1)-1), intent(inout) :: rlist
   integer, intent(out) :: st

   integer :: node
   integer :: i, j
   integer(long) :: ii, jj
   integer(long), dimension(:), allocatable :: ptr
   integer(long), dimension(:), allocatable :: nptr
   integer, dimension(:), allocatable :: col

   allocate(ptr(n+2), stat=st)
   if(st.ne.0) return
   ptr(:) = 0

   ! Count number of entries in each row. ptr(i+2) = #entries in row i
   do node = 1, nnodes
      do ii = rptr(node), rptr(node+1)-1
         j = rlist(ii) ! row entry
         ptr(j+2) = ptr(j+2) + 1
      end do
   end do

   ! Determine row starts. ptr(i+1) = start of row i
   ptr(1:2) = 1
   do i = 1, n
      ptr(i+2) = ptr(i+1) + ptr(i+2)
   end do

   jj = ptr(n+2)-1 ! total number of entries
   allocate(col(jj), stat=st)
   if(st.ne.0) return

   ! Now fill in col array
   do node = 1, nnodes
      do ii = rptr(node), rptr(node+1)-1
         j = rlist(ii) ! row entry
         col( ptr(j+1) ) = node
         ptr(j+1) = ptr(j+1) + 1
      end do
   end do

   ! Finally transpose back into nodes
   allocate(nptr(nnodes))
   nptr(:) = rptr(1:nnodes)
   do i = 1, n
      do jj = ptr(i), ptr(i+1)-1
         node = col(jj)
         rlist(nptr(node)) = i
         nptr(node) = nptr(node) + 1
      end do
   end do
end subroutine dbl_tr_sort

!
! This subroutine applies the permutation perm to order, invp and cc
!
subroutine apply_perm(n, perm, order, invp, cc, block_pivots)
   integer, intent(in) :: n
   integer, dimension(n), intent(in) :: perm
   integer, dimension(n), intent(inout) :: order
   integer, dimension(n), intent(inout) :: invp
   integer, dimension(n), intent(inout) :: cc
   integer, dimension(n), optional, intent(inout) :: block_pivots

   integer :: i
   integer :: j

   ! Use order as a temporary variable to permute cc. Don't care about cc(n+1)
   order(1:n) = cc(1:n)
   do i = 1, n
      j = perm(i)
      cc(j) = order(i)
   end do

   ! Use order as a temporary variable to permute invp.
   order(1:n) = invp(1:n)
   do i = 1, n
      j = perm(i)
      invp(j) = order(i)
   end do

   ! Use order as a temporary variable to permute block_pivots if present
   if(present(block_pivots)) then
      order(1:n) = block_pivots(1:n)
      do i = 1, n
         j = perm(i)
         block_pivots(j) = order(i)
      end do
   endif

   ! Recover order as inverse of invp
   do i = 1, n
      order(invp(i)) = i
   end do
end subroutine apply_perm

! calculate 64-bit sum of 32-bit ints
! (because the F2003 sum intrinsic doens't have a kind argument)
integer(long) pure function sum_long(array)
   integer, dimension(:), intent(in) :: array

   integer :: i

   sum_long = 0
   do i = 1, size(array)
      sum_long = sum_long + array(i)
   end do
end function sum_long

end module hsl_mc78_integer
