! COPYRIGHT (c) 2012 Science and Technology Facilities Council
! Original date 6 June 2012 as version 1.0.0
!
! Version 1.1.0
! See ChangeLog for version history.
!
! Written by: Jonathan Hogg and Jennifer Scott

! Given a sparse symmetric  matrix A, this package 
! uses a matching algorithm to compute an elimination
! order that is suitable for use with a sparse direct solver. 
! It optionally computes scaling factors.

! Note: this version uses mc64 to compute the matching
!       and follows approach of Duff and Pralet (2005).
!       It does not compute the optimal matching in the
!       structurally singular case as we found the extra work
!       required to do this did not result in a better
!       pivot order.

! To convert from double to single:
! * Change wp
! * Replace _double by _single
! * Replace call to mc64wd to mc64w

module hsl_mc80_double
   use hsl_mc34_double
   use hsl_mc68_integer
   ! Also calls mc64
   implicit none

   private
   public :: mc80_order, mc80_order_full, mc80_control, mc80_info

   integer, parameter :: wp = kind(0d0)
   real(wp), parameter :: one = 1.0_wp
   real(wp), parameter :: zero = 0.0_wp
   real(wp), parameter :: rinf = huge(0.0_wp)

   ! Error flags
   integer, parameter :: MC80_SUCCESS               = 0
   integer, parameter :: MC80_ERROR_ALLOCATION      = -1
   integer, parameter :: MC80_ERROR_A_N_OOR         = -2
   integer, parameter :: MC80_ERROR_SINGULAR        = -3
   integer, parameter :: MC80_ERROR_MC68            = -4
   integer, parameter :: MC80_ERROR_ORD_OOR         = -5
   integer, parameter :: MC80_ERROR_NO_METIS        = -6

   ! warning flags
   integer, parameter :: MC80_WARNING_SINGULAR      = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Data type for information returned by code
   !
   type mc80_info
      integer :: compress_rank = 0 ! order of compressed matrix (ordering
         ! applied to this matrix)
      integer :: flag = 0 ! Takes one of the enumerated flag values:
         ! Possible  following values:
         !    0 : successful entry (for structurally nonsingular matrix).
         !   -1 : allocation error
         !   -2 : Metis not available when required
         !   -3 : either n or ord is out of range (immediate return)
         !   -4 : unexpected error returned by hsl_mc68
         !   +1 : successful entry (for structurally singular matrix).
      integer :: flag68 = 0 ! error flag from hsl_mc68
      integer :: max_cycle = 0 ! maximum cycle length
      integer :: struct_rank = 0 ! structural rank 
      integer :: stat = 0 ! stat parameter

   end type mc80_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   type mc80_control
      logical :: action = .true. ! if set to .true. and the matrix
         ! is found to be structurally singular, the code exits immediately
         ! with an error flag set.
      logical :: unmatched_scale_zero = .false. ! If true, then in the singular
         ! case, the scaling factors associated with the unmatched part of the
         ! matrix are set to zero. This will breakdown if the matched submatrix
         ! is numerically singular, but if no any factorization will likely have
         ! far fewer delayed pivots. If false, then Duff and Pralet are followed
         ! and the value is set such that the corresponding scaled entries are
         ! <= 1 in absolute value.
      logical :: unmatched_last = .false. ! If true, then in singular case, rows
         ! and columns associated with the unmatched part are ordered last in
         ! the elimination order. If false, they are placed in a fill-minimising
         ! position.
   end type mc80_control

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   interface mc80_order
      ! to be used if lower triangle of A available
      module procedure mc80_order_double
   end interface

   interface mc80_order_full
      ! to be used if lower and upper triangles of A available
      module procedure mc80_order_full_double
   end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This routine computes the matching, splits long cycles,
! compresses the matrix, applies AMD, Min Deg or Metis to the compressed
! matrix and then returns an ordering for the original
! matrix. This ordering flags 2x2 pivots using negative signs
! (as used on input to hsl_ma77).
! If scale is present, the mc64 scaling factors are returned.
! The matrix may be singular.
!
! Input (ptr, row , val) is the ** lower triangular part ** of the matrix.
! Diagonal entries need not be present.
! There is no matrix data checking.
!
subroutine mc80_order_double(ord, n, ptr, row, val, order, &
      control, info, scale)
   integer, intent(in) :: ord ! controls ordering on compressed matrix
         ! 1 An approximate minimum degree ordering is used.
         ! 2 A minimum degree ordering is used (as in MA27).
         ! 3 METIS ordering with default settings is used.
   integer, intent(in) :: n
   integer, dimension(:), intent(in) :: ptr
   integer, dimension(:), intent(in) :: row
   real(wp), dimension(:), intent(in) :: val
   integer, dimension(:), intent(out) :: order ! |order(i)|  holds the position
      ! of variable i in the elimination order (pivot sequence). If a
      ! 1x1 pivot i is obtained,  order(i)>0. If a
      ! 2x2 pivot involving  i and j is obtained, 
      ! order(i)<0, order(j)<0 and |order(j)|=|order(i)|+1. 
   type(mc80_control), intent(in) :: control
   type (mc80_info), intent(out) :: info ! used to hold information
   real(wp), dimension(n), intent(out), optional :: scale ! if present,
      ! returns the mc64 symmetric scaling

   integer, dimension(:), allocatable :: cperm ! used to hold matching
   integer, dimension(:), allocatable :: ptr2 ! column pointers for expanded
      ! matrix. 
   integer, dimension(:), allocatable :: row2 ! row indices for expanded matrix.
   real(wp), dimension(:), allocatable :: val2 ! entries of expanded matrix.
   real(wp), dimension(:), allocatable :: scale2 ! holds scaling factors
      ! if scale is not present.

   integer :: i, j, k, ne

   info%compress_rank = 0
   info%flag = 0
   info%flag68 = 0
   info%stat = 0
   info%max_cycle = 0 
   info%struct_rank = n

   ! check n has valid value
   if (n < 0) then
      info%flag = MC80_ERROR_A_N_OOR
      return
   end if

   ! check ord has valid value
   if (ord < 1 .or. ord > 3) then
      info%flag = MC80_ERROR_ORD_OOR
      return
   end if

   ! just return with no action if n = 0
   if (n.eq. 0) return

   !
   ! Take absolute values, expand out, removing any explicit zeroes
   !
   ne = ptr(n+1) - 1
   allocate(ptr2(n+1), row2(2*ne), val2(2*ne), cperm(n), stat=info%stat)
   if (info%stat.ne.0) then
      info%flag = MC80_ERROR_ALLOCATION
      return
   end if

   k = 1
   do i = 1, n
      ptr2(i) = k
      do j = ptr(i), ptr(i+1)-1
         if (val(j).eq.zero) cycle
         row2(k) = row(j)
         val2(k) = abs(val(j))
         k = k + 1
      end do
   end do
   ptr2(n+1) = k

   call mc34_expand(n, row2, ptr2, cperm, a=val2)
   
   ! Compute matching and scaling

   if (present(scale)) then 
      call mc80_scale(n,ptr2,row2,val2,scale,control,info, &
           perm=cperm)
   else
      allocate(scale2(n), stat=info%stat)
      if (info%stat.ne.0) then
         info%flag = MC80_ERROR_ALLOCATION
         return
      end if
      call mc80_scale(n,ptr2,row2,val2,scale2,control,info, &
           perm=cperm)
      deallocate(scale2, stat=info%stat)
   end if
   deallocate(val2, stat=info%stat)

   if (info%flag.lt.0) return

   ! Note: row j is matched with column cperm(j)
   !
   ! Split matching into 1- and 2-cycles only and then
   ! compress matrix and order.

   call mc80_split(ord,n,row2,ptr2,order,cperm,control,info)

   if(present(scale)) scale(1:n) = exp( scale(1:n) )

end subroutine mc80_order_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This routine is the same as mc80_order EXCEPT on
! input ptr, row , val hold the ** lower AND upper ** triangular 
! parts of the matrix.
! this reduces amount of copies of matrix required (so slightly
! more efficient on memory and does not need to expand supplied matrix)
!  
subroutine mc80_order_full_double(ord, n, ptr, row, val, order,  &
      control, info, scale)
   integer, intent(in) :: ord ! controls ordering on compressed matrix
         ! 1 An approximate minimum degree ordering is used.
         ! 2 A minimum degree ordering is used (as in MA27).
         ! 3 METIS ordering with default settings is used.
   integer, intent(in) :: n
   integer, dimension(:), intent(in) :: ptr
   integer, dimension(:), intent(in) :: row
   real(wp), dimension(:), intent(in) :: val
   integer, dimension(:), intent(out) :: order ! |order(i)|  holds the position
      ! of variable i in the elimination order (pivot sequence). If a
      ! 1x1 pivot i is obtained,  order(i)>0. If a
      ! 2x2 pivot involving  i and j is obtained, 
      ! order(i)<0, order(j)<0 and |order(j)|=|order(i)|+1. 
   type(mc80_control), intent(in) :: control
   type (mc80_info), intent(out) :: info ! used to hold information
   real(wp), dimension(n), intent(out), optional :: scale ! if present,
      ! returns the mc64 symmetric scaling

   integer, dimension(:), allocatable :: cperm ! used to hold matching
   integer, dimension(:), allocatable :: ptr2 ! column pointers for expanded 
      ! matrix. 
   integer, dimension(:), allocatable :: row2 ! row indices for expanded matrix 
   real(wp), dimension(:), allocatable :: val2 ! entries of expanded matrix.
   real(wp), dimension(:), allocatable :: scale2 ! holds scaling factors
      ! if scale is not present.

   integer :: i, j, k, ne

   info%compress_rank = 0
   info%flag = 0
   info%flag68 = 0
   info%stat = 0
   info%max_cycle = 0 
   info%struct_rank = n

   ! check n has valid value
   if (n < 0) then
     info%flag = MC80_ERROR_A_N_OOR
     return
   end if

   ! check ord has valid value
   if (ord < 1 .or. ord > 3) then
     info%flag = MC80_ERROR_ORD_OOR
     return
   end if

   ! just return with no action if n = 0
   if (n.eq.0) return

   !
   ! Remove any explicit zeroes and take absolute values
   !
   ne = ptr(n+1) - 1
   allocate(ptr2(n+1), row2(ne), val2(ne), cperm(n), stat=info%stat)
   if (info%stat.ne.0) then
      info%flag = MC80_ERROR_ALLOCATION
      return
   end if

   k = 1
   do i = 1, n
      ptr2(i) = k
      do j = ptr(i), ptr(i+1)-1
         if (val(j).eq.zero) cycle
         row2(k) = row(j)
         val2(k) = abs(val(j))
         k = k + 1
      end do
   end do
   ptr2(n+1) = k

   ! Compute matching and scaling

   if (present(scale)) then 
      call mc80_scale(n,ptr2,row2,val2,scale,control,info, &
         perm=cperm)
   else
      allocate(scale2(n), stat=info%stat)
      if (info%stat.ne.0) then
         info%flag = MC80_ERROR_ALLOCATION
         return
      end if
      call mc80_scale(n,ptr2,row2,val2,scale2,control,info, &
         perm=cperm)
      deallocate(scale2, stat=info%stat)
   end if
   deallocate(val2, stat=info%stat)

   if (info%flag.lt.0) return

   ! Note: row j is matched with column cperm(j)
   ! write (*,'(a,15i4)') 'cperm',cperm(1:min(15,n))
   !
   ! Split matching into 1- and 2-cycles only and then
   ! compress matrix and order.

   call mc80_split(ord,n,row2,ptr2,order,cperm,control,info)

   if(present(scale)) scale(1:n) = exp( scale(1:n) )

end subroutine mc80_order_full_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Split matching into 1- and 2-cycles only and then
! compress matrix and order using mc68.
!
! Input (ptr2, row2 , val2) holds the ** lower and upper triangles ** 
! of the matrix (with explicit zeros removed).
! Overwritten in the singular case
!
subroutine mc80_split(ord,n,row2,ptr2,order,cperm,control,info)
   integer, intent(in) :: ord ! controls choice of ordering algorithm 
      ! on compressed matrix
   integer, intent(in) :: n
   integer, dimension(:), intent(in) :: ptr2 
   integer, dimension(:), intent(in) :: row2 
   integer, dimension(n), intent(out) :: order ! used to hold ordering
   integer, dimension(n), intent(inout) :: cperm ! used to hold matching 
   type (mc80_control), intent(in) :: control
   type (mc80_info), intent(inout) :: info ! used to hold information

   type (mc68_control) :: control68
   type (mc68_info) :: info68

   integer, dimension(:), allocatable :: iwork ! work array
   integer, dimension(:), allocatable :: old_to_new, new_to_old
      ! holds mapping between original matrix indices and those in condensed
      ! matrix.
   integer, dimension(:), allocatable :: ptr3 ! column pointers for condensed 
      ! matrix.
   integer, dimension(:), allocatable :: row3 ! row indices for condensed 
      ! matrix.


   integer :: csz ! current cycle length
   integer :: i, j, j1, j2, jj, k, krow
   integer :: max_csz ! maximum cycle length
   integer :: ncomp ! order of compressed matrix
   integer :: ncomp_matched ! order of compressed matrix (matched entries only)
   integer :: ne ! number of non zeros

   ! Use iwork to track what has been matched:
   ! -2 unmatched
   ! -1 matched as singleton
   !  0 not yet seen
   ! >0 matched with specified node

   ne = ptr2(n+1) - 1
   allocate(ptr3(n+1), row3(ne), old_to_new(n), new_to_old(n), iwork(n), &
      stat=info%stat)
   if (info%stat.ne.0) return

   iwork(1:n) = 0
   max_csz = 0
   do i = 1, n
      if (iwork(i).ne.0) cycle
      j = i
      csz = 0
      do
         if (cperm(j).eq.-1) then
            ! unmatched by MC64
            iwork(j) = -2
            csz = csz + 1
            exit
         else if (cperm(j).eq.i) then
            ! match as singleton, unmatched or finished
            iwork(j) = -1
            csz = csz + 1
            exit
         end if
         ! match j and cperm(j)
         jj = cperm(j)
         iwork(j) = jj
         iwork(jj) = j
         csz = csz + 2
         ! move onto next start of pair
         j = cperm(jj)
         if (j.eq.i) exit
      end do
      max_csz = max(max_csz, csz)
   end do
   info%max_cycle = max_csz

   ! Overwrite cperm with new matching
   cperm(1:n) = iwork(1:n)

   !
   ! Build maps for new numbering schemes
   !
   k = 1
   do i = 1,n
      j = cperm(i)
      if (control%unmatched_last .and. j.eq.-2) cycle
      if (j<i .and. j.gt.0) cycle
      old_to_new(i) = k
      new_to_old(k) = i ! note: new_to_old only maps to first of a pair
      if (j.gt.0) old_to_new(j) = k   
      k = k + 1
   end do
   ncomp_matched = k-1
   
   ! Place unmatched columns to rear
   if(control%unmatched_last) then
      do i = 1, n
         if(cperm(i).eq.-2) then
            old_to_new(i) = k
            new_to_old(k) = i
            k = k + 1
         endif
      end do
   end if

   !
   ! Produce a condensed version of the matrix for ordering.
   ! Hold pattern using ptr3 and row3.
   !
   ptr3(1) = 1
   iwork(:) = 0 ! Use to indicate if entry is in a paired column
   ncomp = 1
   jj = 1
   do i = 1, n
      j = cperm(i)
      if (j<i .and. j.gt.0) cycle ! already seen
      if (control%unmatched_last .and. j.eq.-2) cycle ! column not participating
      do k = ptr2(i), ptr2(i+1)-1
         krow = old_to_new(row2(k))
         if (iwork(krow).eq.i) cycle ! already added to column
         if (krow>ncomp_matched) cycle ! unmatched row not participating
         row3(jj) = krow
         jj = jj + 1
         iwork(krow) = i
      end do
      if (j.gt.0) then
         ! Also check column cperm(i)
         do k = ptr2(j), ptr2(j+1)-1
            krow = old_to_new(row2(k))
            if (iwork(krow).eq.i) cycle ! already added to column
            if (krow>ncomp_matched) cycle ! unmatched row not participating
            row3(jj) = krow
            jj = jj + 1
            iwork(krow) = i
         end do
      end if
      ptr3(ncomp+1) = jj
      ncomp = ncomp + 1
   end do
   ncomp = ncomp - 1
   info%compress_rank = ncomp

   if(control%unmatched_last) ncomp = ncomp_matched

   ! store just lower triangular part for input to hsl_mc68
   ptr3(1) = 1
   jj = 1
   j1 = 1
   do i = 1, ncomp
      j2 = ptr3(i+1)
      do k = j1, j2-1
         krow = row3(k)
         if ( krow.lt.i ) cycle ! already added to column
         row3(jj) = krow
         jj = jj + 1
      end do
      ptr3(i+1) = jj
      j1 = j2
   end do

   ! reorder the compressed matrix using hsl_mc68.
   ! switch off hsl_mc68 printing
   control68%lp = -1
   control68%wp = -1
   control68%mp = -1
   control68%print_level = -1
   call mc68_order(ord,ncomp,ptr3,row3,order,control68,info68)

   if (info68%flag < 0) then
      info%flag68 = info68%flag
      select case(info68%flag)
      case(-1)
         info%flag = MC80_ERROR_ALLOCATION
         info%stat = info68%stat
      case(-5)
         info%flag = MC80_ERROR_NO_METIS
      case default
         info%flag = MC80_ERROR_MC68
      end select
      return
   end if

   do i = 1, ncomp
      j = order(i)
      iwork(j) = i
   end do

   !
   ! Translate inverse permutation in iwork back to 
   ! permutation for original variables.
   ! Set negative signs for 2x2 pivots (exploited by hsl_ma77).
   !
   k = 1
   do i = 1, ncomp
      j = new_to_old( iwork(i) )
      order(j) = k
      k = k + 1
      if (cperm(j).gt.0) then
         order(j) = -order(j)
         j = cperm(j)
         order(j) = -k
         k = k + 1
      end if
   end do

   ! Place unmatched columns last
   if(control%unmatched_last) then
      do i = ncomp+1, ncomp + (n-info%struct_rank)
         j = new_to_old( i )
         order(j) = k
         k = k + 1
      end do
   endif

end subroutine mc80_split

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Scale the matrix using MC64, accounting for singular matrices using the
! approach of Duff and Pralet
!
! Expects a full matrix as input
!
subroutine mc80_scale(n, ptr, row, val, scale, control, info, perm)
   integer, intent(in) :: n
   integer, dimension(:), intent(in) :: ptr
   integer, dimension(:), intent(in) :: row
   real(wp), dimension(:), intent(in) :: val
   real(wp), dimension(n), intent(out) :: scale ! returns the symmetric scaling
   type(mc80_control), intent(in) :: control
   type(mc80_info), intent(inout) :: info ! used to hold information
   integer, dimension(n), intent(out), optional :: perm ! if present, returns
      ! the matching

   integer, dimension(:), allocatable :: ptr2 ! column pointers after 
      ! zeros removed.
   integer, dimension(:), allocatable :: row2 ! row indices after zeros
      !  removed. 
   real(wp), dimension(:), allocatable :: val2 ! matrix of absolute values
      ! (zeros removed).
   real(wp), dimension(:), allocatable :: cscale ! temporary copy of scaling
      ! factors. Only needed if A rank deficient. allocated to have size n.

   integer :: i, j, k, ne

   info%struct_rank = n

   !
   ! Remove any explicit zeroes and take absolute values
   !
   ne = ptr(n+1) - 1
   allocate(ptr2(n+1), row2(ne), val2(ne),stat=info%stat)
   if (info%stat.ne.0) then
      info%flag = MC80_ERROR_ALLOCATION
      return
   end if

   k = 1
   do i = 1, n
      ptr2(i) = k
      do j = ptr(i), ptr(i+1)-1
         if (val(j).eq.zero) cycle
         row2(k) = row(j)
         val2(k) = val(j)
         k = k + 1
      end do
   end do
   ptr2(n+1) = k

   call mc80_match(n,row2,ptr2,val2,scale,control,info,perm=perm)
   if (info%flag.lt.0) return

   if (info%struct_rank.ne.n) then
      ! structurally singular case. At this point, scaling factors
      ! for rows in corresponding to rank deficient part are set to 
      ! zero. The following is to set them according to Duff and Pralet.
      deallocate(ptr2, stat=info%stat)
      deallocate(row2, stat=info%stat)
      deallocate(val2, stat=info%stat)
      if(.not.control%unmatched_scale_zero) then
         allocate(cscale(n),stat=info%stat)
         if (info%stat.ne.0) then
            info%flag = MC80_ERROR_ALLOCATION
            return
         end if
         cscale(1:n) = scale(1:n)
         do i = 1,n
            if (cscale(i).ne.-huge(scale)) cycle
            do j = ptr(i), ptr(i+1)-1
               k = row(j)
               if (cscale(k).eq.-huge(scale)) cycle
               scale(i) = max(scale(i), val(j)+scale(k))
            end do
            if(scale(i).eq.-huge(scale)) then
               scale(i) = zero
            else
               scale(i) = -scale(i)
            endif
         end do
      end if
   end if

end subroutine mc80_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Input (ptr2, row2 , val2) holds the ** lower and upper triangles ** 
! of the matrix (with explicit zeros removed).
! val2 holds absolute values of matrix entries.
! Overwritten in the singular case
!
subroutine mc80_match(n,row2,ptr2,val2,scale,control,info,perm)
   integer, intent(in) :: n
   integer, dimension(:), intent(inout) :: ptr2 ! In singular case, overwritten
      ! by column pointers for non singular part of matrix.
   integer, dimension(:), intent(inout) :: row2 ! In singular case, overwritten
      ! by row indices for non singular part of matrix.
   real(wp), dimension(:), intent(inout) :: val2 ! In singular case, overwritten
      ! by entries for non singular part of matrix.
   real(wp), dimension(n), intent(out) :: scale ! returns the symmetric scaling
   type(mc80_control), intent(in) :: control
   type (mc80_info), intent(inout) :: info ! used to hold information
   integer, dimension(n), intent(out), optional :: perm ! if present, returns
      ! the matching

   integer, dimension(:), allocatable :: iwork ! work array
   integer, dimension(:), allocatable :: cperm ! used to hold matching
   integer, dimension(:), allocatable :: old_to_new, new_to_old
      ! holds mapping between original matrix indices and those in reduced
      ! non singular matrix. 
   real(wp), dimension(:), allocatable :: cmax ! (log) column maximum
   real(wp), dimension(:), allocatable :: dw ! array used by mc64

   integer :: i, j, j1, j2, jj, k
   integer :: ne ! number of non zeros
   integer :: nn ! Holds number of rows/cols in non singular part of matrix
   integer :: nne ! Only used in singular case. Holds number of non zeros
     ! in non-singular part of matrix.
   integer :: rank ! returned by mc64
   real(wp) :: colmax ! max. entry in col. of expanded matrix

   allocate(iwork(5*n), cperm(n), dw(2*n), cmax(n), stat=info%stat)
   if (info%stat.ne.0) then
      info%flag = MC80_ERROR_ALLOCATION
      return
   end if
    
   ! Compute column maximums    
   do i = 1,n
      colmax = max(zero,maxval(val2(ptr2(i):ptr2(i+1)-1)))
      if (colmax.ne.zero) colmax = log(colmax)
      cmax(i) = colmax
   end do

   do i = 1,n
      val2(ptr2(i):ptr2(i+1)-1) = cmax(i) - log(val2(ptr2(i):ptr2(i+1)-1))
   end do

   ne = ptr2(n+1)-1
   call mc64wd(n,ne,ptr2,row2,val2,cperm,rank, &
      iwork(1),iwork(n+1),iwork(2*n+1),iwork(3*n+1),iwork(4*n+1), &
      dw(1),dw(n+1))

   if (rank.eq.n) then
      do i = 1,n
         scale(i) = (dw(i)+dw(n+i)-cmax(i))/2
      end do
      if (present(perm)) perm(1:n) = cperm(1:n)
      return
   end if

   !!!! we have to handle the singular case. Either immediate exit
   ! or set warning, squeeze out the unmatched entries and recall mc64wd.

   info%struct_rank = rank
   if (.not.control%action) then
      info%flag = MC80_ERROR_SINGULAR
      return
   end if

   info%flag = MC80_WARNING_SINGULAR

   allocate(old_to_new(n), new_to_old(n),stat=info%stat)
   if (info%stat.ne.0) then
      info%flag = MC80_ERROR_ALLOCATION
      return
   end if

   k = 0
   do i = 1,n
      if (cperm(i) < 0) then
         ! row i and col j are not part of the matching
         old_to_new(i) = -1
      else
         k = k + 1
         ! old_to_new(i) holds the new index for variable i after
         ! removal of singular part and new_to_old(k) is the
         ! original index for k
         old_to_new(i) = k
         new_to_old(k) = i
      end if
   end do

   ! Overwrite ptr2, row2 and val2
   nne = 0
   k = 0
   ptr2(1) = 1
   j2 = 1
   do i = 1,n
      j1 = j2
      j2 = ptr2(i+1)
      ! skip over unmatched entries
      if (cperm(i) < 0) cycle
      k = k + 1
      do j = j1,j2-1
         jj = row2(j)
         if (cperm(jj) < 0) cycle
         nne = nne + 1
         row2(nne) = old_to_new(jj)
         val2(nne) = val2(j)
      end do
      ptr2(k+1) = nne + 1
    end do
    ! nn is order of non-singular part.
    nn = k
    call mc64wd(nn,nne,ptr2,row2,val2,cperm,rank, &
       iwork(1),iwork(nn+1),iwork(2*nn+1),iwork(3*nn+1),iwork(4*nn+1), &
       dw(1),dw(nn+1))

    do i = 1,n
       j = old_to_new(i)
       if (j < 0) then
          scale(i) = -huge(scale)
       else
         ! Note: we need to subtract col max using old matrix numbering
         scale(i) = (dw(j)+dw(nn+j)-cmax(i))/2
      end if
   end do

   if (present(perm)) then
      perm(1:n) = -1
      do i = 1,nn
         j = cperm(i)
         perm(new_to_old(i)) = new_to_old(j)
      end do
   end if

end subroutine mc80_match

!**********************************************************************
end module hsl_mc80_double
