! COPYRIGHT (c) 2013-4 Science and Technology Facilities Council,
!                      Academy of Sciences of Czech Republic.
! Version 2.2.2
! Original date February 2013
!
! Written by:  Jennifer Scott    STFC Rutherford Appleton Laboratory
!              Miroslav Tuma     Academy of Sciences of Czech Republic

! To convert from double to single and visa versa:
! * Change wp
! * Change _double
! * Change BLAS calls: dcopy
! * Change HSL calls: mc77id, mc77ad, mc61id, mc61ad
! * Change default for control%small from 10**(-20) to 10**(-12)
!
module hsl_mi28_double

  use hsl_mc64_double
  use hsl_mc68_integer
  use hsl_mc69_double

  implicit none

   private
   public :: mi28_keep, mi28_control, mi28_info
   public :: mi28_factorize, mi28_finalise, mi28_precondition, mi28_solve

  integer, parameter  :: wp = kind(0.0d0)
  integer, parameter :: long = selected_int_kind(18)
  real(wp), parameter :: zero = 0.0_wp
  real(wp), parameter :: one = 1.0_wp
  real(wp), parameter :: sfact = 2.0_wp
  real(wp), parameter :: sfact2 = 4.0_wp
  real(wp), parameter :: alpham = 0.001_wp

   ! Error flags
   integer, parameter :: MI28_ERROR_ALLOCATION      = -1
   integer, parameter :: MI28_ERROR_ROW_TOO_SMALL   = -2
   integer, parameter :: MI28_ERROR_VAL_TOO_SMALL   = -3
   integer, parameter :: MI28_ERROR_N_OOR           = -4
   integer, parameter :: MI28_ERROR_PTR             = -5
   integer, parameter :: MI28_ERROR_MISS_DIAG       = -6
   integer, parameter :: MI28_ERROR_MC77            = -7
   integer, parameter :: MI28_ERROR_MC64            = -8
   integer, parameter :: MI28_ERROR_SINGULAR        = -9
   integer, parameter :: MI28_ERROR_SCALE           = -10
   integer, parameter :: MI28_ERROR_USER_PERM       = -11
   integer, parameter :: MI28_ERROR_MC61            = -12
   integer, parameter :: MI28_ERROR_MC68            = -13
   integer, parameter :: MI28_ERROR_DEALLOCATION    = -14

   ! warning flags
   integer, parameter :: MI28_WARNING_OOR_IDX        = 1
   integer, parameter :: MI28_WARNING_DUP_IDX        = 2
   integer, parameter :: MI28_WARNING_MC64           = 3
   integer, parameter :: MI28_WARNING_MC61           = 4
   integer, parameter :: MI28_WARNING_NEG_DIAG       = 5

  interface mi28_factorize
      module procedure mi28_factorize_double
  end interface

  interface mi28_precondition
      module procedure mi28_precondition_double
  end interface

  interface mi28_solve
      module procedure mi28_solve_double
  end interface

  interface mi28_finalise
      module procedure mi28_finalise_double
  end interface

!****************************************************************************
  type mi28_keep ! The scalar keep holds the preconditioner, scaling and pivot
                 ! sequence for passing from factorize to precondition
   integer(long), allocatable ::  fact_ptr(:) ! column pointers for incomplete 
      ! factor
   integer, allocatable ::  fact_row(:) 
      ! First fact_ptr(n+1)-1 entries hold row indices for incomplete factor
   real(wp), allocatable ::  fact_val(:) 
      ! First fact_ptr(n+1)-1 entries hold entries incomplete factor
   real(wp), allocatable :: scale(:) ! Allocated size n if
      ! scaling required (control%scale >0). Used to hold scaling factors of A.
      ! of A.
      ! The matrix is permuted and then scaled so that
      !  S Q^T A Q S   is factorized.
      ! To get S, we permute the entries of scale at end of factorize. 
   integer, allocatable :: invp(:) ! specifies the pivot order such that invp(j)
      ! holds the j-th pivot column.
   integer, allocatable :: perm(:) ! specifies the pivot order such that 
      ! perm(i) holds the position of i in the pivot sequence (perm and invp
      ! are inverses of each other) 
   real(wp), allocatable :: w(:) ! Allocated size n if control%iorder>0.

  end type mi28_keep

!****************************************************************************
  type mi28_control ! The scalar control of this type controls the action

    real(wp) :: alpha = zero ! initial shift
    logical :: check = .true. ! if set to true, user's data is checked.
         ! Otherwise, no checking and may fail in unexpected way if
         ! there are duplicates/out-of-range entries.
    integer :: iorder = 6 ! controls ordering of A. Options:
!         <=0  no ordering
!           1  RCM
!           2  AMD
!           3  user-supplied ordering
!           4  ascending degree
!           5  Metis
!          >= 6, Sloan (MC61)
     integer :: iscale = 1 ! controls whether scaling is used.
        ! iscale = 1 is Lin and More scaling (l2 scaling)
        ! iscale = 2 is mc77 scaling
        ! iscale = 3 is mc64 scaling
        ! iscale = 4 is diagonal scaling
        ! iscale = 5 user-supplied scaling
        ! iscale <= 0, no scaling
        ! iscale >= 6, Lin and More
     real(wp) :: lowalpha = alpham ! Shift after first breakdown is
        ! max(shift_factor*alpha,lowalpha)
     integer :: maxshift = 3 ! During search for shift, we decrease
        ! the lower bound max(alpha,lowalpha) on the shift by
        ! shift_factor2 at most maxshift times (so this limits the
        ! number of refactorizations that are performed ... idea is
        ! reducing alpha as much as possible will give better preconditioner
        ! but reducing too far will lead to breakdown and then a refactorization
        ! is required (expensive so limit number of reductions)
        ! Note: Lin and More set this to 3.
     logical :: rrt = .false. ! controls whether entries of RR^T that cause no
        ! additional fill are allowed. They are allowed if
        ! rrt = .true. and not otherwise.
     real(wp) :: shift_factor = sfact ! if the current shift is found
        ! to be too small, it is increased by at least a factor of shift_factor.
        ! Values <= 1.0 are treated as default.
     real(wp) :: shift_factor2 = sfact2 ! if factorization is successful
        ! with current (non zero) shift, the shift
        ! is reduced by a factor of shift_factor2.
        ! Values <= 1.0 are treated as default.
     real(wp) :: small = 10.0_wp**(-20)
     real(wp) :: tau1 = 0.001_wp  ! used to select "small" entries that
         ! are dropped from L (but may be included in R). 
     real(wp) :: tau2 = 0.0001_wp ! used to select "tiny" entries that are
         ! dropped from R.  Require
         ! tau2 < tau1 (otherwise, tau2 = 0.0 is used locally).
     integer :: unit_error = 6 ! unit number for error messages.
        ! Printing is suppressed if unit_error  <  0.
     integer :: unit_warning = 6 ! unit number for warning messages.
        ! Printing is suppressed if unit_warning  <  0.

  end type mi28_control

!****************************************************************************

  type mi28_info ! The scalar info of this type returns information to user.

      integer :: band_after = 0 ! semibandwidth after MC61
      integer :: band_before = 0 ! semibandwidth before MC61
      integer :: dup = 0 ! number of duplicated entries found in row.
      integer :: flag = 0 ! error flag
      integer :: flag61 = 0 ! error flag from mc61
      integer :: flag64 = 0 ! error flag from hsl_mc64
      integer :: flag68 = 0 ! error flag from hsl_mc68
      integer :: flag77 = 0 ! error flag from mc77
      integer :: nrestart = 0 ! number of restarts (after reducing the shift)
      integer :: nshift = 0 ! number of non-zero shifts used
      integer :: oor = 0 ! number of out-of-range entries found in row.
      real(wp) :: profile_before = 0 ! semibandwidth before MC61
      real(wp) :: profile_after = 0 ! semibandwidth after MC61
      integer(long) :: size_r = 0_long ! size of arrays jr and ar that are used 
         ! for r   
      integer :: stat = 0 ! Fortran stat parameter
      real(wp) :: alpha = zero ! on successful exit, holds shift used

  end type mi28_info

!****************************************************************************
contains

   subroutine mi28_factorize_double(n, ptr, row, val, lsize, rsize, keep,  &
      control, info, scale, perm)

   ! input matrix
   integer, intent(in) :: n  ! order of matrix
   integer, intent(inout) ::  ptr(n+1) ! column pointers ! 
   integer, intent(inout) ::  row(:) ! row indices of lower triangular part 
     ! with diagonal entry first in each column. Only
     ! changed  if control%check = .true. (any duplicates
     ! are summed and out-of-range are removed and diagonal first in each col.).
   real(wp), intent(inout) ::  val(:) ! entries of lower triangular part.
     ! on exit, reordered in the same way as row.
   integer, intent(in) :: lsize ! if lsize > 0, max. number of entries in
     ! column of factor is lsize + original entries.
   integer, intent(in) :: rsize ! if rsize > 0, max. number of entries in
     ! r is rsize. If rsize = 0, no r is computed.

   type(mi28_keep), intent(out) :: keep ! See derived-type declaration
   type(mi28_control), intent(in) :: control ! See derived-type declaration
   type(mi28_info), intent(out) :: info      ! See derived-type declaration
   real(wp), intent(in), optional :: scale(n) ! Must be present 
      ! and hold scaling factors for A if control%scale =5. 
   integer, intent(in), optional :: perm(n) ! Must be present and 
      ! hold pivot order if control%order = 3.  
      ! perm(i) must hold the position of i in the pivot sequence
   
   !!!!!!!!!!!!!!!!!!!!!
   real(wp) :: alpha_in
   character(50)  :: context      ! Procedure name (used when printing).
   real(wp) :: tau1 ! set to control%tau1
   real(wp) :: tau2 ! set to control%tau2 but set to zero if
      ! control%tau2 > control%tau1.
   integer :: i
   integer :: iorder ! controls ordering
   integer :: iscale ! controls scaling
   integer :: j
   integer :: jm  ! set to 0 if control%rrt = 0 and to 2 otherwise
   integer :: k
   real(wp) :: lowalpha
   integer :: liw
   integer :: nout ! stream for errors
   integer :: ne ! set to ptr(n+1) - 1
   integer :: nout1 ! stream for warnings
   integer :: st ! stat parameter
   real(wp) :: shift_factor
   real(wp) :: shift_factor2
   integer :: job ! job parameter for mc61 (controls ordering algorithm)
   integer :: icntl61(10) ! controls for MC61
   integer :: info61(10) ! information from MC61
   real(wp) :: cntl61(5)
   real(wp) :: rinfo61(15)

   type(mc68_control) :: control68
   type(mc68_info) :: info68

   integer, allocatable, dimension(:) :: iw1,iw ! work array
   integer, allocatable, dimension(:) :: row_perm,ptr_perm
   real(wp), allocatable, dimension(:) :: val_perm

   !!!!!!!!!!!!!!!!!!!!!

   context = 'mi28_factorize'

   ! Set stream numbers
   nout = control%unit_error
   nout1 = control%unit_warning

! initialise data held in info.
   info%alpha = max(zero,control%alpha)
   info%band_before = 0
   info%band_after = 0
   info%dup = 0
   info%flag = 0
   info%flag61 = 0
   info%flag64 = 0
   info%flag68 = 0
   info%flag77 = 0
   info%nrestart = 0
   info%nshift = 0
   info%oor = 0
   info%profile_before = zero
   info%profile_after = zero
   info%size_r = 0_long
   info%stat = 0

   ! Check matrix data
   if (n < 1) then
      info%flag = MI28_ERROR_n_OOR
      call mi28_print_flag(context,nout,info%flag)
      return
   end if

   ne = ptr(n+1) - 1
   if (ne < 1) then
      info%flag = MI28_ERROR_PTR
      call mi28_print_flag(context,nout,info%flag)
      return
   end if

   if (size(row).lt.ne) then
      info%flag = MI28_ERROR_ROW_TOO_SMALL
      call mi28_print_flag(context,nout,info%flag)
      return
   end if

   if (size(val).lt.ne) then
      info%flag = MI28_ERROR_VAL_TOO_SMALL
      call mi28_print_flag(context,nout,info%flag)
      return
   end if

   ! allocate keep%fact_ptr and work space
   deallocate (keep%fact_ptr,stat=st)
   allocate (keep%fact_ptr(n+1),iw1(n+1),stat=st)
   if (st .ne. 0) go to 100

   iorder = control%iorder
   if (iorder.gt.6) iorder = 6 ! Sloan
   if (iorder.lt.0) iorder = 0 ! no ordering

   if (iorder.gt.0) then
      deallocate (keep%invp,stat=st)
      deallocate (keep%perm,stat=st)
      deallocate (keep%w,stat=st)
      allocate (keep%perm(n), keep%invp(n), keep%w(n), stat=st)
      if (st.ne.0) go to 100
   end if

   if (iorder == 3) then
     ! check user-supplied perm and set keep%invp
     if (.not. present(perm)) then
        info%flag = MI28_ERROR_USER_PERM
        call mi28_print_flag(context,nout,info%flag)
        return
     end if
     keep%invp(:) = 0
     do i = 1, n
        j = perm(i)
        if (j < 1 .or. j > n) exit
        if (keep%invp(j) /= 0) exit ! Duplicate found
        keep%invp(j) = i
     end do
     if (i-1 /= n) then
        info%flag = MI28_ERROR_USER_PERM
        call mi28_print_flag(context,nout,info%flag)
        return
     end if
     keep%perm(1:n) = perm(1:n)
   end if

   if (control%check) then
      ! check and clean data (order by row indices within each column)
      call check_matrix(n,ne,ptr,row,val,info)

      if (info%flag.lt.0) then
         call mi28_print_flag(context,nout,info%flag)
         return
      else if (info%flag.gt.0) then
         call mi28_print_flag(context,nout1,info%flag)
      end if

      ! row,ptr,val hold lower triangle of matrix in CSC format,
      ! with diagonal first (and other entries ordered by increasing row index).
      ne = ptr(n+1) - 1
   end if

   if (iorder == 1 .or. iorder == 6) then

       ! Use mc61 for RCM  or Sloan algorithm
       liw = 8*n + 2
       allocate (row_perm(2*ne),val_perm(ne),ptr_perm(n+1), &
                 iw(liw),stat=st)
       if (st.ne.0) go to 100

       row_perm(1:ne) = row(1:ne)
       ptr_perm(1:n+1) = ptr(1:n+1)

       call mc61id(icntl61,cntl61)
       icntl61(1:2) = -1 ! suppress errors/warnings
       ! note: we do not handle duplicates or out-of-range ... they
       ! are treated as an error as we have already checked for these
       ! (and if we handle them here we won't be dealing with the reals)

       job = 1
       if (iorder == 1) job = 2 ! RCM
       call mc61ad(job,n,2*ne,row_perm,ptr_perm,keep%perm,liw,iw,keep%w, &
            icntl61,cntl61,info61,rinfo61)

      deallocate (iw,stat=st)

      if (info61(1) < 0) then 
         info%flag = MI28_ERROR_MC61
         info%flag61 = info61(1)
         call mi28_print_flag(context,nout,info%flag,flag61=info%flag61)
         return
      end if

      info%band_before = int(rinfo61(3))
      info%band_after = int(rinfo61(7))
      info%profile_before = rinfo61(1)
      info%profile_after = rinfo61(5)

      if (job.eq.1 .and. rinfo61(1).le.rinfo61(5)) then
        ! failed to reduce profile so identity permutation and work
        ! with original matrix.
         iorder = 0
      else if (job.eq.2 .and. rinfo61(3).le.rinfo61(7)) then
        ! failed to reduce semibandwidth so identity permutation and work
        ! with original matrix.
         iorder = 0
      end if
      if (iorder.eq.0) then
        do i = 1,n
          keep%invp(i) = i
          keep%perm(i) = i
        end do
        info%flag = MI28_WARNING_MC61
        call mi28_print_flag(context,nout1,info%flag)
      end if

   else if (iorder == 2 .or. iorder == 5) then

! AMD or Metis using hsl_mc68
      job = 1
      if (iorder == 5) job = 3 ! Metis

      ! suppress printing
      control68%lp = -1
      control68%mp = -1
      control68%wp = -1

      call mc68_order(job,n,ptr,row,keep%perm,control68,info68)

      ! check for errors (but ignore deallocation error)
      if (info68%flag < 0 .and. info68%flag.ne.-2) then
         info%flag68 = info68%flag
         info%flag = MI28_ERROR_MC68
         call mi28_print_flag(context,nout,info%flag,flag68=info%flag68)
         return
      end if

    else if (iorder == 4) then
       ! order by ascending degree (used by Edmond Chow).
       ! compute degrees (hold in iw1)
       iw1(1:n) = 0
       do i = 1,n
          do j = ptr(i),ptr(i+1)-1
             k = row(j)
             iw1(i) = iw1(i) + 1
             if (k.ne.i) iw1(k) = iw1(k) + 1
          end do
       end do

       call kb07ai(iw1,n,keep%invp)
       do i = 1,n
          j = keep%invp(i)
          keep%perm(j) = i
       end do

    end if

    if (iorder >= 1) then
       ! Permute the matrix (and put diagonal entry first in each column).
       ! Allocate space for permuted matrix (we do not overwrite user's
       ! data as user may require original matrix).
       if (.not. allocated(row_perm)) then
         allocate (row_perm(ne),val_perm(ne),ptr_perm(n+1),stat=st)
         if (st.ne.0) go to 100
       end if

! Permutation held in keep%perm. keep%perm(i) holds position of variable i in 
! list of permuted variables. We return the inverse permutation (so that
! invp(j) holds index of row of A that is permuted to row j ie
! if invp(j) = i, j-th pivot is i).
      if (iorder.ne.4 .and. iorder.ne.3) then
        do i = 1,n
          j = keep%perm(i)
          keep%invp(j) = i
        end do
      end if

      call permute_matrix(n,ne,ptr,row,val,keep%perm,ptr_perm,row_perm, &
            val_perm,iw1)

   else if (.not. control%check) then
         ! no checking of data but will check diag. entry is positive.
        do i = 1,n
          j = ptr(i)
          if (ptr(i+1).eq.ptr(i)) then
            info%flag = MI28_ERROR_MISS_DIAG
            call mi28_print_flag(context,nout1,info%flag)
            return
          else if (row(j).ne.i) then
            ! diagonal is not first entry in column
            info%flag = MI28_ERROR_MISS_DIAG
            call mi28_print_flag(context,nout1,info%flag)
            return
          end if
          if (val(j).le.zero) &
            info%flag = MI28_WARNING_NEG_DIAG
       end do
       if (info%flag.eq.MI28_WARNING_NEG_DIAG) &
          call mi28_print_flag(context,nout1,info%flag)
   end if

   deallocate(iw1,stat=st)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! compute scaling if requested.

   iscale = control%iscale
   if (iscale.gt.5) iscale = 1 ! use default
   if (iscale.lt.0) iscale = 0

   deallocate (keep%scale,stat=st)
   if (iscale.gt.0) then
      allocate (keep%scale(n),stat=st)
   else 
      ! scaling not used so allocate keep%scale to have size 1
      allocate(keep%scale(1),stat=st)
   end if
   if (st.ne.0) go to 100

   if (iscale.eq.5) then
     if (.not. present(scale)) then
       info%flag = MI28_ERROR_SCALE
       call MI28_print_flag(context,nout,info%flag)
       return
     end if
     ! user supplied scaling is for A. 
     ! Permute to give scaling for permuted matrix
     if (iorder.gt.0) then
        do i = 1,n
           keep%scale(i) = scale(keep%invp(i))
        end do
      else
         keep%scale(1:n) = scale(1:n)
      end if
   end if

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! check control%tau2
   tau1 = abs(control%tau1)
   tau2 = abs(control%tau2)
   if (control%tau2 .gt. tau1) tau2 = zero

   shift_factor = control%shift_factor
   if (shift_factor.le.one) shift_factor = sfact
   shift_factor2 =control%shift_factor2
   if (shift_factor2.le.one) shift_factor = sfact2

   alpha_in = max(zero,control%alpha)
   lowalpha = control%lowalpha
   if (lowalpha.le.zero) lowalpha = alpham
   if (alpha_in > zero) lowalpha = alpha_in

   ! set jm (to control what happens to entries of RR^T)
   jm = 2
   if (control%rrt) jm = 0

   if (iorder.eq.0) then
     ! we work with original matrix ptr,row,val
     call ichsd(n, ptr, row, val, keep%fact_ptr, keep%fact_row,           &
          keep%fact_val, tau1, tau2, jm,                                  &
          control%maxshift, shift_factor, shift_factor2,control%small,    &
          lsize, rsize, iscale, nout, info, alpha_in, lowalpha,           &
          keep%scale)
  else
     ! we work with permuted matrix ptr_perm,row_perm,val_perm
     call ichsd(n, ptr_perm, row_perm, val_perm, keep%fact_ptr, keep%fact_row, &
          keep%fact_val, tau1, tau2, jm,                                  &
          control%maxshift, shift_factor, shift_factor2,control%small,    &
          lsize, rsize, iscale, nout, info, alpha_in, lowalpha,           &
          keep%scale)
  end if

    if (st .ne. 0) go to 100

    if (info%flag .lt. 0) return

    if (iscale.eq.0) then
       deallocate(keep%scale,stat=st)
    else if (iorder.gt.0) then
       ! at this point, keep%scale has scaling factors for the permuted
       ! matrix. To allow the user to input the computed scaling on
       ! a subsequent call, we need to return scaling factors for A.
       do i = 1,n
          keep%w(i) = keep%scale(keep%perm(i))
       end do
       keep%scale(1:n) = keep%w(1:n)
    end if

100 if (st .ne. 0) then
      info%stat=st
      info%flag = MI28_ERROR_ALLOCATION
      call mi28_print_flag(context,nout,info%flag,st = info%stat)
      return
    end if

   end subroutine mi28_factorize_double

!****************************************************************************

   subroutine mi28_precondition_double(n, keep, z, y, info)

! Purpose: to perform the preconditioning operation y = (LL^T)^{-1}z
! ie solve LL^Ty = z, where L is incomplete factor held in CSC format.

      integer, intent(in) :: n
      type(mi28_keep), intent(inout) :: keep
      real(wp), intent(in) :: z(n)
      real(wp), intent(out) :: y(n)
      type(mi28_info), intent(inout) :: info

      integer :: i
      integer(long) :: j,jstrt,jstop
      integer :: k
      real(wp) :: temp

      info%flag = 0

      if (allocated(keep%scale) .and. allocated(keep%invp)) then
       ! recall that keep%scale has scaling factors for A so permute
       ! to get scaling factors for the permuted matrix that has been factorized
        do i = 1,n
          j = keep%invp(i)
          keep%w(i) = keep%scale(j) * z(j)
        end do
!  Forward substitution
        do i = 1,n
          jstrt = keep%fact_ptr(i)
          jstop = keep%fact_ptr(i+1) - 1
          temp = keep%w(i)*keep%fact_val(jstrt)
          do j = jstrt+1,jstop
            k = keep%fact_row(j)
            keep%w(k) = keep%w(k) - keep%fact_val(j)*temp
          end do
          keep%w(i) = temp
        end do

! Back substitution
        do i = n,1,-1
          jstrt = keep%fact_ptr(i)
          jstop = keep%fact_ptr(i+1) - 1
          temp = zero
          do j = jstrt+1,jstop
            k = keep%fact_row(j)
            temp = temp - keep%fact_val(j)*keep%w(k)
          end do
          keep%w(i) = (keep%w(i) + temp)*keep%fact_val(jstrt)
        end do

        do i = 1,n
          j = keep%invp(i)
          y(j) = keep%scale(j) * keep%w(i)
        end do
!!!!!!!!!!!!!!
      else if (allocated(keep%scale)) then
          do i = 1,n
            y(i) = keep%scale(i) * z(i)
          end do
!  Forward substitution
        do i = 1,n
          jstrt = keep%fact_ptr(i)
          jstop = keep%fact_ptr(i+1) - 1
          temp = y(i)*keep%fact_val(jstrt)
          do j = jstrt+1,jstop
            k = keep%fact_row(j)
            y(k) = y(k) - keep%fact_val(j)*temp
          end do
          y(i) = temp
        end do

! Back substitution
        do i = n,1,-1
          jstrt = keep%fact_ptr(i)
          jstop = keep%fact_ptr(i+1) - 1
          temp = zero
          do j = jstrt+1,jstop
            k = keep%fact_row(j)
            temp = temp - keep%fact_val(j)*y(k)
          end do
          y(i) = (y(i) + temp)*keep%fact_val(jstrt)
        end do

        do i = 1,n
          y(i) = keep%scale(i) * y(i)
        end do
!!!!!!!!!!!!!!
      else if (allocated(keep%invp)) then
        do i = 1,n
          j = keep%invp(i)
          keep%w(i) = z(j)
        end do

!  Forward substitution
        do i = 1,n
          jstrt = keep%fact_ptr(i)
          jstop = keep%fact_ptr(i+1) - 1
          temp = keep%w(i)*keep%fact_val(jstrt)
          do j = jstrt+1,jstop
            k = keep%fact_row(j)
            keep%w(k) = keep%w(k) - keep%fact_val(j)*temp
          end do
          keep%w(i) = temp
        end do

! Back substitution
        do i = n,1,-1
          jstrt = keep%fact_ptr(i)
          jstop = keep%fact_ptr(i+1) - 1
          temp = zero
          do j = jstrt+1,jstop
            k = keep%fact_row(j)
            temp = temp - keep%fact_val(j)*keep%w(k)
          end do
          keep%w(i) = (keep%w(i) + temp)*keep%fact_val(jstrt)
        end do

        do i = 1,n
          j = keep%invp(i)
          y(j) = keep%w(i)
        end do
!!!!!!!!!!!!!!
      else
         ! copy z into y
        call dcopy(n,z,1,y,1)

!  Forward substitution
        do i = 1,n
          jstrt = keep%fact_ptr(i)
          jstop = keep%fact_ptr(i+1) - 1
          temp = y(i)*keep%fact_val(jstrt)
          do j = jstrt+1,jstop
            k = keep%fact_row(j)
            y(k) = y(k) - keep%fact_val(j)*temp
          end do
          y(i) = temp
        end do

! Back substitution
        do i = n,1,-1
          jstrt = keep%fact_ptr(i)
          jstop = keep%fact_ptr(i+1) - 1
          temp = zero
          do j = jstrt+1,jstop
            k = keep%fact_row(j)
            temp = temp - keep%fact_val(j)*y(k)
          end do
          y(i) = (y(i) + temp)*keep%fact_val(jstrt)
        end do

      end if

   end subroutine mi28_precondition_double

!****************************************************************************

   subroutine mi28_solve_double(trans, n, keep, z, y, info)

! Purpose: to solve Ly = SQ^Tz or L^TS^{-1}Q^Ty = z
! where L is incomplete factor held in CSC format, 
! S is scaling and Q permutation used in computing L.

      logical, intent(in) :: trans !true if L^T system required, false otherwise
      integer, intent(in) :: n
      type(mi28_keep), intent(inout) :: keep
      real(wp), intent(in) :: z(n)
      real(wp), intent(out) :: y(n)
      type(mi28_info), intent(inout) :: info

      integer :: i
      integer(long) :: j,jstrt,jstop
      integer :: k
      real(wp) :: temp

      info%flag = 0

      if (.not.trans) then
!  Forward substitution

        if (allocated(keep%scale) .and. allocated(keep%invp)) then
          do i = 1,n
            j = keep%invp(i)
            y(i) = keep%scale(j) * z(j)
          end do
        else if (allocated(keep%scale)) then
          do i = 1,n
            y(i) = keep%scale(i) * z(i)
          end do
        else if (allocated(keep%invp)) then
          do i = 1,n
            j = keep%invp(i)
            y(i) = z(j)
          end do
        else
          call dcopy(n,z,1,y,1)
        end if

        do i = 1,n
          jstrt = keep%fact_ptr(i)
          jstop = keep%fact_ptr(i+1) - 1
          temp = y(i)*keep%fact_val(jstrt)
          do j = jstrt+1,jstop
            k = keep%fact_row(j)
            y(k) = y(k) - keep%fact_val(j)*temp
          end do
          y(i) = temp
        end do

      else

        if (allocated(keep%invp)) then
! Back substitution
          call dcopy(n,z,1,keep%w,1)
          do i = n,1,-1
            jstrt = keep%fact_ptr(i)
            jstop = keep%fact_ptr(i+1) - 1
            temp = zero
            do j = jstrt+1,jstop
              k = keep%fact_row(j)
              temp = temp - keep%fact_val(j)*keep%w(k)
            end do
            keep%w(i) = (keep%w(i) + temp)*keep%fact_val(jstrt)
          end do

          if (allocated(keep%scale)) then
            do i = 1,n
              j = keep%invp(i)
              y(j) = keep%scale(j) * keep%w(i)
            end do
          else
            do i = 1,n
              j = keep%invp(i)
              y(j) = keep%w(i)
            end do
          end if

        else
          call dcopy(n,z,1,y,1)
          do i = n,1,-1
            jstrt = keep%fact_ptr(i)
            jstop = keep%fact_ptr(i+1) - 1
            temp = zero
            do j = jstrt+1,jstop
              k = keep%fact_row(j)
              temp = temp - keep%fact_val(j)*y(k)
            end do
            y(i) = (y(i) + temp)*keep%fact_val(jstrt)
          end do

          if (allocated(keep%scale)) then
            do i = 1,n
              y(i) = keep%scale(i) * y(i)
            end do
          end if

        end if

      end if

   end subroutine mi28_solve_double

!****************************************************************************

   subroutine mi28_finalise_double(keep, info)

! Purpose: to deallocate components of keep.

      type(mi28_keep), intent(inout) :: keep
      type(mi28_info), intent(inout) :: info

      info%flag = 0

      if (allocated (keep%invp)) &
        deallocate (keep%invp,stat=info%stat)
      if (info%stat.ne.0) go to 10
      if (allocated (keep%perm)) &
        deallocate (keep%perm,stat=info%stat)
      if (info%stat.ne.0) go to 10
      if (allocated (keep%scale)) &
        deallocate (keep%scale,stat=info%stat)
      if (info%stat.ne.0) go to 10
      if (allocated (keep%w)) &
        deallocate (keep%w,stat=info%stat)
      if (info%stat.ne.0) go to 10
      if (allocated (keep%fact_row)) &
        deallocate (keep%fact_row,stat=info%stat)
      if (info%stat.ne.0) go to 10
      if (allocated (keep%fact_ptr)) &
        deallocate (keep%fact_ptr,stat=info%stat)
      if (info%stat.ne.0) go to 10
      if (allocated (keep%fact_val)) &
        deallocate (keep%fact_val,stat=info%stat)
      if (info%stat.ne.0) go to 10

  10  if (info%stat.ne.0) info%flag = MI28_ERROR_DEALLOCATION

   end subroutine mi28_finalise_double
!****************************************************************************
!
!   incomplete Cholesky decomposition wrapper

      subroutine ichsd(n,ia,ja,aa,il,jl,al,tau1,tau2,jm,             &
        maxshift,shift_factor,shift_factor2,small,                   &
        lsize,rsize,iscale,nout,info,alpha_in,lowalpha,scalesj)

      ! input matrix
      integer, intent(in) :: n  ! order of matrix
      integer, intent(in) ::  ia(n+1) ! column pointers for lower triangle.
      integer, intent(in) ::  ja(:) ! row indices of lower triangular part.
        ! diagonal entry must be first entry in each col. 
      real(wp), intent(in) ::  aa(:) ! entries of lower triangular part

      ! computed incomplete factorization
      integer(long), intent(out) ::  il(n+1) ! column pointers for incomplete 
         ! factor
      integer, intent(out), allocatable ::  jl(:) ! row indices for incomplete 
         ! factor
      real(wp), intent(out), allocatable ::  al(:) ! entries of incomplete
         ! factor

      ! arguments that control computation of incomplete factorization
      real(wp), intent(in) :: tau1 ! used to select "small" 
         ! entries to drop.
      real(wp), intent(in) :: tau2 ! used to select "tiny" entries.
      integer, intent(in) :: jm ! controls whether RR^T entries computed
      integer, intent(in) :: maxshift! max number of shifts allowed in 
         ! Lin and More
      real(wp), intent(in) :: shift_factor
      real(wp), intent(in) :: shift_factor2
      real(wp), intent(in) :: small ! = control%small
      integer, intent(in) :: lsize ! max. number of entries in
        ! column of factor is lsize + original entries.
      integer, intent(in) :: rsize ! if rsize > 0, max. number of entries in
        ! r is rsize. If rsize = 0, no r is computed.

      integer, intent(in) :: iscale ! controls scaling. 

      integer, intent(in) :: nout ! stream for errors
      type(mi28_info), intent(inout) :: info    ! See derived-type declaration
      real(wp), intent(in) :: alpha_in ! Holds initial shift
      real(wp), intent(inout) :: lowalpha ! Lower bound on shift (gets altered
         ! if shift altered)
      real(wp), intent(inout) :: scalesj(:) ! Holds scaling factors.

      ! local arrays
      integer, allocatable, dimension(:) :: list ! used by ictkl3. 
         ! allocated to have size n+1
      integer(long), allocatable, dimension(:) :: start ! used by ictkl3. 
         ! allocated to have size n+1
      integer, allocatable, dimension(:) :: wn02 ! used by ictkl3. 
         ! allocated to have size n+1
      real(wp), allocatable, dimension(:) :: wr01 ! used by ictkl3. 
         ! allocated to have size n+1

      ! the following are allocated if rsize .ne. 0
      integer, allocatable, dimension(:) :: listr ! used by ictkl3. 
         ! allocated to have size n+1
      integer(long), allocatable, dimension(:) :: startr ! used by ictkl3. 
         ! allocated to have size n+1
      integer(long), allocatable ::  ir(:) ! column pointers for r. allocated to
         ! have size n+1
      integer, allocatable ::  jr(:) ! row indices for r. allocated to
         ! have size max_r 
      real(wp), allocatable ::  ar(:) ! entries of r. allocated to
         ! have size max_r 

      real(wp), allocatable ::  d(:) ! holds the diagonal entries of matrix
         ! and update as factorization proceeds (detect breakdown asap)

      ! if we have a shift that gives a successful factorization but we hope
      ! to reduce the shift, we copy the successful incomplete factorization
      ! into the followinf arrays (note: if we can't successfully allocate them
      ! then we must recompute the factorization)
      integer(long), allocatable ::  il_copy(:) 
      integer, allocatable ::  jl_copy(:) 
      real(wp), allocatable ::  al_copy(:) 


      ! local scalars
      real(wp) :: alpha ! manteuffel shift
      real(wp) :: alpha_old ! initialise to 0 and set to an alpha that works
         ! if alpha is reduced
      character(50)  :: context      ! Procedure name (used when printing).
     ! real(wp) :: dsfl
      integer :: flag  ! error flag
      integer :: flag_previous  ! error flag using previous shift
      integer :: i
      integer(long) :: ii
      integer :: lsize_local ! set to max(0,lsize). allows us to handle lsize <0
      integer(long) :: max_p ! Holds size of jl and al.
         ! jl and al will be allocated to have size max_p.
      integer(long) :: max_r ! Holds size of jr and ar.
         ! jr and ar will be allocated to have size max_r
      integer :: idiag ! set to index of diagonal entry ie first entry in column
      integer :: j
      integer :: ne ! number of entries in A (lower triangle)
      integer(long) :: nec ! set to il_copy(n+1)-1
      integer :: st  ! stat parameter
      logical :: success ! set to .true. if we have a successful factorization
         ! held in il_copy, jl_copy, al_copy
      real(wp) :: temp

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      context = 'mi28_factorize'

      if (iscale.gt.0 .and. iscale.lt.5) then

         call compute_scaling(n,ia,ja,aa,iscale,scalesj,st,info)

         if (st.ne.0) go to 500
         if (info%flag.lt.0) then
            if (info%flag77.ne.0) &
              call mi28_print_flag(context,nout,info%flag,flag77=info%flag77)
            if (info%flag64.ne.0) &
              call mi28_print_flag(context,nout,info%flag,flag64=info%flag64)
            return
         end if

         if (info%flag.eq.MI28_WARNING_MC64) then
            call mi28_print_flag(context,nout,info%flag)
         end if

      end if

!  -- allocate space for the preconditioner (ensure not already allocated first)
!     Treat lsize < 0 as 0.
      lsize_local = max(0,lsize)
      lsize_local = min(lsize_local,n-1)
      ne = ia(n+1) - 1
      max_p = ne + lsize_local*int(n,long)

      deallocate (jl,stat=st)
      deallocate (al,stat=st)
      allocate(jl(max_p),al(max_p),stat=st)
      if (st.ne.0) go to 500
 
      allocate(start(n+1),list(n+1),wn02(n+1),wr01(n+1),d(n),stat=st)
      if (st.ne.0) go to 500

      ! allocate arrays for R (not needed if rsize = 0)
      if (rsize.ne.0) then
         allocate(startr(n+1),listr(n+1),stat=st)
         if (st.ne.0) go to 500
         max_r = max(0,rsize)*int(n,long)
         allocate(ir(n+1),jr(max_r),ar(max_r),stat=st)
         if (st.ne.0) go to 500
         info%size_r = max_r
      end if

       !  -- get the Manteuffel shift alpha

      ! compute the initial shift. Note: we have already checked
      ! that the digaonal (which is the first entry in the column) is present.
      if (alpha_in.ne.zero) then
         alpha = max(alpha_in,zero)
      else
         alpha = zero
         do i = 1,n
            idiag = ia(i)
            temp = aa(idiag)
            if (temp.eq.zero) then
               alpha = lowalpha
            else
               if (iscale.ne.0) temp = temp*(scalesj(i)**2)
               ! if there is a negative diag entry, temp < 0
               ! and so we set alpha > 0
               alpha = max(alpha,-temp)
            end if
         end do
         ! at this point, alpha = zero unless one or more
         ! negative diagonal entries was found. 
         if (alpha.gt.zero) alpha = max(alpha,lowalpha)
      end if

      flag = 0
      alpha_old  = zero
      success = .false.

shift1: do
        flag_previous = flag

        ! store the (scaled and shifted) diagonal entries
        ! Note: we have a separate copy of the diagonal entries as
        ! we keep them up-to-date as the factorization progresses. This
        ! allows us to abort the factorization as soon as a non positive
        ! entry is encountered.
        if (iscale.gt.0) then
          do j = 1,n
            ii = ia(j)
            d(j) = scalesj(j)*aa(ii)*scalesj(j) + alpha
          end do
        else
          do j = 1,n
            ii = ia(j)
            d(j) = aa(ii) + alpha
          end do
        end if

!    -- perform the decomposition
        if (rsize.eq.0) then
           call ictkl3(n,ia,ja,aa,iscale,scalesj,d,il,jl,al,lsize_local,rsize, &
             tau1,tau2,small,jm,wn02,wr01,start,list,flag)

        else
           call ictkl3(n,ia,ja,aa,iscale,scalesj,d,il,jl,al,lsize_local,rsize, &
             tau1,tau2,small,jm,wn02,wr01,start,list,flag, &
             ir,jr,ar,startr,listr)

        end if

!    -- check the output flag and alpha

           if (flag.eq.0 .and. alpha.eq.lowalpha) then

              ! preconditioner computed and the shift is at the lower bound.
              ! might get away with smaller shift so try decreasing shift

              if (info%nrestart .lt. maxshift) then
                 alpha_old = alpha
                 lowalpha = lowalpha/shift_factor2
                 if (info%nrestart.gt.0 .and. flag_previous.eq.0) &
                    lowalpha = lowalpha/shift_factor2
                 alpha = lowalpha
                 info%nrestart = info%nrestart + 1
                 info%nshift   = info%nshift + 1

                 ! If we have space, keep the factorization we already have
                 st = 0
                 if (.not.allocated(il_copy)) then
                    allocate (il_copy(n+1),stat=st)
                    if (st.ne.0) cycle
                 end if
                 il_copy(1:n+1) = il(1:n+1)
                 nec = il_copy(n+1) - 1
                 if (.not.allocated(jl_copy)) then
                    allocate (jl_copy(nec),stat=st)
                    if (st.ne.0) cycle
                 else if (size(jl_copy).lt.nec) then
                    deallocate (jl_copy,stat=st)
                    allocate (jl_copy(nec),stat=st)
                    if (st.ne.0) cycle
                 end if
                 jl_copy(1:nec) = jl(1:nec)
                 if (.not.allocated(al_copy)) then
                    allocate (al_copy(nec),stat=st)
                    if (st.ne.0) cycle
                 else if (size(al_copy).lt.nec) then
                    deallocate (al_copy,stat=st)
                    allocate (al_copy(nec),stat=st)
                    if (st.ne.0) cycle
                 end if
                 al_copy(1:nec) = al(1:nec)
                 success = .true.
                 cycle
              else
                 ! reached the limit on number of restarts
                 exit shift1
              end if

           else if (flag.lt.0) then

              ! preconditioner has not been computed successfully.
              ! if we already have a factorization, restore it

              if (success) then
                 nec = il_copy(n+1) - 1
                 il(1:n+1) = il_copy(1:n+1) 
                 jl(1:nec) = jl_copy(1:nec)
                 al(1:nec) = al_copy(1:nec)
                 flag = 0
                 alpha = alpha_old
                 exit shift1
              end if

              ! Increase the shift

              if (alpha_old.ne.zero) then
                 ! alpha_old worked ok so go back to it
                 alpha = alpha_old
              else

               alpha = max(shift_factor*alpha,lowalpha)

              ! if error is at same step as last time, increase
              ! shift more rapidly (perhaps also if flag is close
              ! to flag_previous so little progress has been made)
              if (flag_previous.ne.0) then
                if (abs(flag-flag_previous).le.n/10 .and. &
                    n+flag.ge.n/10) &
                    alpha = max(shift_factor*alpha,lowalpha)
              end if

              end if
              info%nshift = info%nshift + 1

           end if

         if (flag.ge.0) exit shift1

      end do shift1

      info%alpha = alpha

 500  info%stat=st
      if (info%stat .ne. 0) then
         info%flag = MI28_ERROR_ALLOCATION
         call mi28_print_flag(context,nout,info%flag,st = info%stat)
      end if

      end subroutine ichsd
!****************************************************************************
!
! routine to print errors and warnings
!
subroutine mi28_print_flag(context,nout,flag,st,flag61,flag64,flag68,flag77)
   integer, intent(in) :: flag, nout
   integer, intent(in), optional :: st
   integer, intent(in), optional :: flag61,flag64,flag68,flag77
   character (len = *), optional, intent(in) :: context

   if (nout < 0) return
   if (flag < 0) then
      write (nout,'(/3a,i3)') ' Error return from ',trim(context),&
         '. Error flag = ', flag
   else
      write (nout,'(/3a,i3)') ' Warning from ',trim(context),&
         '. Warning flag = ',flag
   end if

   select case(flag)
   case(MI28_ERROR_ALLOCATION)
      if (present(st)) then
         write (nout,'(a,i6)') ' Allocation error. stat parameter = ', st
      else
         write (nout,'(a)') ' Allocation error'
      end if

   case(MI28_ERROR_ROW_TOO_SMALL)
         write (nout,'(a)') ' Size of array row is less than ptr(n+1)-1'

   case(MI28_ERROR_VAL_TOO_SMALL)
         write (nout,'(a)') ' Size of array val is less than ptr(n+1)-1'

   case(MI28_ERROR_n_OOR)
         write (nout,'(a)') ' n is out of range'

   case(MI28_ERROR_PTR)
         write (nout,'(a)') ' ptr(n+1)-1 is less than 1 or ptr is not monotonic'

   case(MI28_ERROR_MC77)
      write (nout,'(a,i3)') &
      ' Unexpected error return from MC77 (scaling routine). flag77 = ',flag77

   case(MI28_ERROR_MC64)
      write (nout,'(a,i3)') &
      ' Unexpected error return from HSL_MC64 (scaling routine). flag64 =', &
        flag64

   case(MI28_ERROR_SINGULAR)
      write (nout,'(a)') &
         ' Matrix found by HSL_MC64 to be structurally singular'

   case(MI28_ERROR_SCALE)
         write (nout,'(a)') ' scale should be present'

   case(MI28_ERROR_USER_PERM)
         write (nout,'(a)') ' Either perm not present when expected'
         write (nout,'(a)') ' or perm does not hold a permutation'

   case(MI28_ERROR_MC61)
      write (nout,'(a,i3)') &
      ' Unexpected error return from  MC61 (RCM or Sloan ordering). flag61 =', &
        flag61

   case(MI28_ERROR_MC68)
      write (nout,'(a,i3)') &
       ' Unexpected error return from HSL_MC68 (ordering routine). flag68 =', &
        flag68

   case(MI28_ERROR_MISS_DIAG)
         write (nout,'(a)') &
      ' One or more diagonal entries is missing'

   case(MI28_WARNING_OOR_IDX)
         write (nout,'(a)') ' out-of-range indices in array row removed'

   case(MI28_WARNING_DUP_IDX)
         write (nout,'(a)') &
 ' Duplicated indices in array row removed; corresponding entries in val summed'

   case(MI28_WARNING_MC64)
         write (nout,'(a)') &
       ' Computed scaling factors are large and so scaling is not used'

   case(MI28_WARNING_MC61)
         write (nout,'(a)') &
       ' MC61 did not improve ordering. Original ordering retained'

   case(MI28_WARNING_NEG_DIAG)
         write (nout,'(a)') &
       ' One or more diagonal entries of matrix is negative.'

   end select

end subroutine mi28_print_flag

!************************************************
 
!  left-looking incomplete Cholesky (ic) preconditioner by value
!  with positive semidefinite stabilization.
!  sparse left-looking implementation.
!  uses two drop tolerances:
!  1/ tau1: used to select "small" entries that are dropped from L
!  2/ tau2: used to select "tiny" entries that are dropped from R
!  So all entries of L are at least tau1 in absolute value, and those in
!  R are all smaller than those in L and are at least tau2.

      subroutine ictkl3(n,ia,ja,aa,iscale,scalesj,d,il,jl,al,lsize,rsize, &
        tau1,tau2,small,jm,wn02,wr01,startl,listl,info,ir,jr,ar,startr,listr)

      integer, intent(in) :: n ! order of A

      integer, intent(in) ::  ia(n+1) ! column pointers for lower triangle.
      integer, intent(in) ::  ja(:) ! row indices of lower triangular part.
        ! diagonal entry must be first entry in each col. 
      real(wp), intent(in) ::  aa(:) ! entries of lower triangular part
      integer, intent(in) :: iscale ! controls scaling
      real(wp), intent(in) :: scalesj(:) ! scaling factors
      real(wp), intent(inout) :: d(n) ! scaled and shifted diagonal entries  
         ! (updated as factorization proceeds)
      integer(long), intent(out) :: il(n+1) ! column pointers for incomplete
         ! factor. 
      integer, intent(out) :: jl(:) ! row indices for incomplete factor. 
      real(wp), intent(out) :: al(:) ! entries of incomplete factor. 
      integer, intent(in) :: lsize ! max. number of entries in
        ! column of factor is lsize + original entries.
      integer, intent(in) :: rsize ! if rsize > 0, max. number of entries in
        ! R is rsize. If rsize = 0, no R is computed.

      real(wp), intent(in) :: tau1 ! larger of the two drop tolerances
         ! Entries of magnitude at least tau2 but less than tau1 are
         ! dropped from the factor but may be included in R.
         ! Entries of magnitude at least tau1 are kept in the factor
         ! (if there is room for them ... otherwise just largest are kept
         ! and those that are dropped may be used in R)
      real(wp), intent(in) :: tau2 ! smaller of the two drop tolerances
         ! (no entries with magnitude < tau2 are kept)
      real(wp), intent(in) :: small ! determines when a pivot is too small
      integer :: jm ! controls whether entries of RR^T are computed.

      integer, intent(out) :: listl(n+1) ! linked list of columns
      integer(long), intent(out) :: startl(n+1) ! starting pointers for 
         ! the fully formed columns.  used also as an auxiliary vector 
         ! for the fresh column
      integer, intent(out) :: wn02(n+1)
      real(wp), intent(out) :: wr01(n+1)

      integer, intent(out) :: info ! error flag 
 
      ! These are used in the case rsize > 0 (if rsize = 0
      ! they are never allocated). As they are not needed if rsize = 0,
      ! they have been made optional
      integer, optional, intent(inout) :: listr(n+1) ! linked list of 
         ! columns for r
      integer(long), optional, intent(inout) :: startr(n+1) ! 
      integer(long), optional, intent(inout) :: ir(n+1) ! column pointers for R.
      integer, optional, intent(inout) :: jr(:) ! row indices for R.  
      real(wp), optional, intent(inout) :: ar(:)! values of entries in R.

      integer :: a_jstrt,a_jstop ! start and end of column j of A
      integer :: asize ! number of entries in current column of A
      real(wp) :: diagj ! diag entry col j
      real(wp) :: diagjtmp ! set to diagj or sqrt(diagj)
      integer :: i
      integer(long) :: istrt,istop,ii
      integer :: ind ! number of entries in current column after dropping
      integer :: ind2  ! number of entries in current column of A (not
        ! including the diagonal) ... then updated to entries in current
        ! column of L
      integer :: indr
      integer :: j
      integer(long) :: jj ! position of diagonal of current column within il
      integer :: jsize ! number of entries to be included oin col. j of L
      integer(long) :: k1,kstrt
      integer :: k
      integer :: largep !
      integer :: larger !
      integer(long) :: lind !
      integer :: lnext,lnextrow 
      integer :: lcurr ! current column in factor
      integer :: rcurr ! current column in r
      integer :: radd
      integer(long) :: rind
      integer(long) :: rstrt,rstop
      integer :: rnext,rnextrow
      integer :: saved,savedr
      real(wp) :: temp

!     initialise
      info = 0

!  -- startl: starting pointers for the fully formed columns
!  -- used also as an auxiliary vector for the fresh column
!  -- listl: linked list 
      startl(1:n) = 0_long
      listl(1:n) = 0

!  -- initialize the loop of computed columns
      saved = 0 
      al(1) = aa(1) ! copy the diagonal entry into place

      a_jstrt = ia(1) 
      ! as columns of factor are computed, they will be put into the
      ! front of the data structure and so set il for column 1 of factor 
      il(1) = 1

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! start with the special case rsize = 0 (we do not need any arrays
      ! for holding R in this case and so it is simple and fast but
      ! quality of preconditioner not as good as rsize > 0)
 
      if (rsize == 0) then 
        ! -- perform the Cholesky decomposition 
  main1:do j = 1,n

          ! check size of diagonal entry.
          ! (note: diaganal entry is up-to-date)
          ! If necessary, exit for another shift
          diagj = d(j)
          if (diagj.le.small) then
            info = -j
            return
          end if

          ind2 = 0
          a_jstop = ia(j+1) - 1
          ! Loop over off-diagonal entries in column j of a.
          ! Scatter the (scaled) entries into full array wr01.
          ! Set startl(i) = 1 for each row i that has an entry
          ! in col. j and set wn01(ind2) to hold the row index of the
          ! ind2-th entry in column j
          if (iscale.gt.0) then 
            do k = a_jstrt+1,a_jstop
              i = ja(k)
              wr01(i) = scalesj(j)*aa(k)*scalesj(i)
              startl(i) = 1_long
              ind2 = ind2 + 1
              wn02(ind2) = i
            end do
          else
            do k = a_jstrt+1,a_jstop
              i = ja(k)
              wr01(i) = aa(k)
              startl(i) = 1_long
              ind2 = ind2 + 1
              wn02(ind2) = i
            end do
          end if

!      -- loop over previous columns of L that update current column j
 
          lcurr = listl(j)
          do 
            if (lcurr.eq.0) exit
            kstrt = startl(lcurr)
            temp = al(kstrt)
            istrt = kstrt + 1
            istop = il(lcurr+1) - 1

!        -- update startl and listl for the row lnextrow
            lnext = listl(lcurr)
            if (istrt.lt.istop) then
              startl(lcurr) = istrt
              lnextrow = jl(istrt)
              listl(lcurr) = listl(lnextrow)
              listl(lnextrow) = lcurr
            end if
 
!        -- perform the updates on col j of L from column lcurr
            ! loop over entries in col. lcurr
            do k1 = istrt,istop
              i = jl(k1)
              if (startl(i).ne.0_long) then 
!          -- update the entry with row index i (it is already non zero) 
                wr01(i) = wr01(i) - temp*al(k1)
              else
!            -- we have new fill in in row i of col j.
             !  Increment count ind2 of entries in col j
             !  Store index (i) of fill entry and compute its value
                ind2 = ind2 + 1
                startl(i) = int(ind2,long)
                wn02(ind2) = i
                wr01(i) = -temp*al(k1)
              end if
            end do

            lcurr = lnext
 
          end do
          ! ind2 is now the number of entries in col j of L after fill in.
          ! The indices of these entries are in wn02.
 
!      -- scale column j of L, optionally drop small off-diagonal
          ! entries and reset startl to 0

          diagjtmp = one/sqrt(diagj)

          jj = il(j)
          jl(jj) = j
          al(jj) = diagjtmp

          if (tau1.gt.zero) then
            ind = 0
            do i = 1,ind2
              k = wn02(i)
              temp = wr01(k)
              temp = temp*diagjtmp
              if (abs(temp).ge.tau1 .or. k.eq.j) then
                ! keep entries of size at least tau1 and diagonal entry
                ind = ind + 1
                wn02(ind) = k
                wr01(k) = temp
              end if
              startl(k) = 0_long
            end do
            ! ind2 is the number of entries in column j after dropping
            ind2 = ind
          else
            do i = 1,ind2
              k = wn02(i)
              wr01(k) = wr01(k)*diagjtmp
              startl(k) = 0_long
            end do
          end if
 
!      -- get the space for the entries to be stored.
          ! jsize is number of (off-diagonal) entries to be included in L.

          asize = a_jstop - a_jstrt
          jsize = min(asize+lsize+saved,ind2)
          saved = saved + asize + lsize - jsize

          !!!! if we set saved = 0, we are not passing on any spare space
          ! eg get spare space in col. 1 since never any fill in col. 1
          ! so that col 1 of L has same pattern as A (independent of lsize).
          ! using saved = 0 allows direct comparison with Lin and More
          ! also for some problems, seem to get worse preconditioner
          ! if saved > 0 allowed (eg DNVS/m_t1)
           saved = 0
 
          largep = ind2 - jsize + 1
          if (ind2.ge.1) then
 
!          -- select jsize large entries at the end of wn02.
!          On exit, entries in wr01 corresponding to indices wn02(largep:ind2) 
!          are all larger than or equal
!          to the entries corresponding to indices in wn02(1:largep).
!          (indirect addressing used ... wr01 is not altered).
 
            call split2(ind2,largep,wn02,wr01)

!          -- sort accepted entries
            call srtshai2(jsize,wn02(largep))
 
          end if

!        -- gather the accepted scattered entries into col. j of factor
!         And update diagonal entries and check for breakdown.

          lind = il(j) + 1
          do i = largep,ind2
             k = wn02(i)
             temp = wr01(k)
             d(k) = d(k) - temp**2
             if (d(k) .le. small .and. k > j) then
                ! diagonal entry too small so return for new shift
                info = -j
                return
             end if
             al(lind) = temp
             jl(lind) = k
             lind = lind + 1
          end do

!      -- finalize the loop. 

          if (j.lt.n) then 
          ! set a_jstrt, il for next column
            a_jstrt = a_jstop + 1
            il(j+1) = lind

            if (jsize.gt.0) then
              startl(j) = il(j) + 1
              lnext = jl(il(j) + 1)
              listl(j) = listl(lnext)
              listl(lnext) = j
            end if

          else
            il(n+1) = lind
          end if

        end do main1

        return   ! we have finished the special case rsize = 0

      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! we are now dealing with the case rsize > 0.

        savedr = 0
!  -- initialize for matrix R
        ir(1) = 1
        listr(1:n) = 0

!    -- loop over the columns 
 main2 :do j = 1,n

          diagj = d(j)
          ! test for breakdown
          if (diagj.le.small) then
             info = -j
             return
          end if

          ind2 = 0
          a_jstop = ia(j+1) - 1
          ! Loop over off-diagonal entries in column j of a.
          ! Scatter the (scaled) entries into full array wr01.
          ! Set startl(i) = 1 for each row i that has an entry
          ! in col. j and set wn01(ind2) to hold the row index of the
          ! ind2-th entry in column j
          if (iscale.gt.0) then 
            do k = a_jstrt+1,a_jstop
              i = ja(k)
              wr01(i) = scalesj(j)*aa(k)*scalesj(i)
              startl(i) = 1_long
              ind2 = ind2 + 1
              wn02(ind2) = i
            end do
          else
            do k = a_jstrt+1,a_jstop
              i = ja(k)
              wr01(i) = aa(k)
              startl(i) = 1_long
              ind2 = ind2 + 1
              wn02(ind2) = i
            end do
          end if


          lcurr = listl(j)
          do 
            if (lcurr.eq.0) exit
            lnext = listl(lcurr)
            kstrt = startl(lcurr)
            temp = al(kstrt)
            istrt = kstrt + 1
            istop = il(lcurr+1) - 1

!        -- update listl for the row lnextrow 
            if (istrt.le.istop) then
              startl(lcurr) = istrt
              lnextrow = jl(istrt)
              listl(lcurr) = listl(lnextrow)
              listl(lnextrow) = lcurr
            end if
 
!        -- gather updates to L from previous L columns 
            do ii = istrt,istop
              k = jl(ii)
              if (startl(k).ne.0_long) then
                wr01(k) = wr01(k) - al(ii)*temp
              else
                ! create fill in
                ind2 = ind2 + 1
                wn02(ind2) = k
                startl(k) = 1_long
                wr01(k) = -al(ii)*temp
              end if
            end do
 
!        -- gather updates to R from previous L columns
            rstrt = startr(lcurr)
            rstop = ir(lcurr+1) - 1 
            do ii = rstrt,rstop
              k = jr(ii)
              if (k.gt.j) then
                if (startl(k).ne.0_long) then
                  wr01(k) = wr01(k) - ar(ii)*temp
                else
                  ind2 = ind2 + 1
                  wn02(ind2) = k
                  startl(k) = 1_long
                  wr01(k) = -ar(ii)*temp
                end if
              end if
            end do
 
!        -- back to next column having a nonzero in a row i
            lcurr = lnext
 
          end do

          rcurr = listr(j)
          do 
            if (rcurr.eq.0) exit
            rnext = listr(rcurr)
            kstrt = startr(rcurr)
            temp = ar(kstrt)
            rstrt = kstrt + 1
            rstop = ir(rcurr+1) - 1
 
!        -- update listr for the row rnextrow 
            if (rstrt.le.rstop) then
              startr(rcurr) = rstrt
              rnextrow = jr(rstrt)
              listr(rcurr) = listr(rnextrow)
              listr(rnextrow) = rcurr
            end if
 
!        -- gather L updates from previous R columns
            istrt = startl(rcurr)
            istop = il(rcurr+1) - 1
            do ii = istrt,istop
              k = jl(ii)
              if (k.gt.j) then
                if (startl(k).ne.0_long) then
                  wr01(k) = wr01(k) - al(ii)*temp
                else
                  ind2 = ind2 + 1
                  wn02(ind2) = k
                  startl(k) = 1_long
                  wr01(k) = -al(ii)*temp
                end if
              end if
            end do
 
            if (jm == 0) then
!             -- gather R updates from previous R columns (RR^T)
              do ii = kstrt,rstop
                k = jr(ii)
                if (k.gt.j) then
                  if (startl(k).ne.0_long) &
                    ! entry all already non zero so allow update.
                    wr01(k) = wr01(k) - ar(ii)*temp
                else if (k.eq.j) then
                  diagj = diagj - ar(ii)*temp
                end if
              end do
            end if
!        -- back to next column of R having a nonzero in a row i
            rcurr = rnext
 
          end do
 
!      -- Store diag. entry of L, scale column j and reset startl to 0
          if (jm.eq.0 .and. diagj.le.small) then
             info = j
             return
          end if
          diagjtmp = one/sqrt(diagj) 

          jj = il(j) 
          jl(jj) = j
          al(jj) = diagjtmp

          do ii = 1,ind2
            k = wn02(ii)
            wr01(k) = wr01(k)*diagjtmp
            startl(k) = 0_long
          end do

          ! jsize is number of entries to be included in L.
          asize = a_jstop - a_jstrt
          jsize = min(asize+lsize+saved,ind2)
          saved = saved + asize + lsize - jsize

          largep = ind2 - jsize + 1

          lind = il(j) + 1
          rind = ir(j)

          if (ind2.ge.1) then
!            -- select jsize large entries at the end of wn02
              call split2(ind2,largep,wn02,wr01)

!            -- sort accepted entries
              call srtshai2(jsize,wn02(largep))

!            -- even accepted entries are filtered
              indr = largep
              ! we will test diagonal entry as soon as updated
              if (tau1.eq.zero) then 
                 ! tau1 = 0 so tau2 = 0 and no dropping of small entries
                 do i = largep,ind2
                    k = wn02(i)
                    temp = wr01(k)
                    d(k) = d(k) - temp**2
                    if (d(k) .le. small .and. k > j) then
                       ! diagonal entry too small so return for new shift
                       info = -j
                       return
                    end if
                    al(lind) = temp
                    jl(lind) = k
                    lind = lind + 1
                 end do

              else
                 do i = largep,ind2
                    k = wn02(i)
                    temp = wr01(k)
                    if (abs(temp).ge.tau1) then
                       ! entry is large enough to keep in L
                       d(k) = d(k) - temp**2
                       if (d(k) .le. small .and. k > j) then
                          ! diagonal entry too small so return for new shift
                          info = -j
                          return
                       end if
                       al(lind) = temp
                       jl(lind) = k
                       lind = lind + 1
                    else if (abs(temp).ne.zero .and. abs(temp).ge.tau2) then
                       ! entry is small but not too tiny for R 
                       wn02(indr) = k
                       indr = indr + 1
                    end if
                 end do
              end if

!           -- sort R entries by values
             radd = min(indr-1,rsize+savedr)
             larger = indr - 1 - radd + 1

             call split2(indr-1,larger,wn02,wr01)

             savedr = savedr + rsize - radd

             do ii = larger,indr-1
                k= wn02(ii)
                temp = wr01(k)
                wr01(k) = zero
                ar(rind) = temp
                jr(rind) = k
                rind = rind + 1
             end do

             ! throw away the smallest entries
             do ii = 1,larger-1
                 k = wn02(ii)
                 wr01(k) = zero
             end do
          end if

!        -- sort the remaining R entries
          call uxvsr4(int(rind-ir(j)),jr(ir(j):rind-1),ar(ir(j):rind-1))

          if (j.lt.n) then
            ! set a_jstrt, il and jl for next column
            a_jstrt = a_jstop + 1
            il(j+1) = lind
            ir(j+1) = rind

            ! define startl, startr and linked list entries (remember no 
            ! diagonal entry for R)
            startl(j) = il(j) + 1
            if (lind.gt.il(j)+1) then
              lnext = jl(il(j)+1)
              listl(j) = listl(lnext)
              listl(lnext) = j
            end if

            if (rind.gt.ir(j)) then
              startr(j) = ir(j)
              rnext = jr(ir(j))
              listr(j) = listr(rnext)
              listr(rnext) = j
            else
              startr(j) = ir(j) + 1
            end if
          else
            il(n+1) = lind
          end if

        end do main2

      end subroutine ictkl3

!****************************************************************************

!   ascending shellsort of an integer vector.
!   additional real vector is permuted in the same way.

      subroutine uxvsr4(k,arrayi,arrayr)

      integer, intent(in) :: k
      integer, intent(inout) :: arrayi(*) ! vector to be sorted
      real(wp), intent(inout) :: arrayr(*) ! vector to be permuted in same way

      integer :: is,la,lt,ls,lls,i,j,js,khalf,ii
      real(wp) :: lb

      if (k.le.1) return

!  -- shellsort initialization
 
      khalf = k/2
      call mylog2(2*k,ii)
      ls = 2**ii-1

!  -- shellsort loop
 
      do lt = 1,ii
        if (ls.gt.khalf) then
          ls = ls/2
        else
          lls = k - ls
          do i = 1,lls
            is = i + ls
            la = arrayi(is)
            lb = arrayr(is)
            j = i
            js = is
 100        continue
            if (la.ge.arrayi(j)) then
              arrayi(js) = la
              arrayr(js) = lb
              cycle
            else
              arrayi(js) = arrayi(j)
              arrayr(js) = arrayr(j)
              js = j
              j = j - ls
            end if
            if (j.ge.1) go to 100
            arrayi(js) = la
            arrayr(js) = lb
          end do
          ls = ls/2
        end if
      end do

      end subroutine uxvsr4

!****************************************************************************

!   split array x by magnitudes of its entries so that
!   entries with indices from ix(ncut:n) are all larger or equal
!   than the entries with indices in ix(1:ncut).
!   reals are indirectly addressed
!   ideas from Saad, Freund & Nachtigal, Lin & More.
!   this routine is especially useful if we have to get only a few
!   large entries as we typically need for preconditioning

      subroutine split2(n,ncut,ix,x)

      integer, intent(in) :: n,ncut
      integer, intent(inout) :: ix(n)
      real(wp), intent(in) :: x(*)

      real(wp) :: midval
      integer :: first,mid,last
      integer :: j,itemp

      if (n.le.1 .or. ncut.le.0 .or. ncut.gt.n) return

        first = 1
        last = n

          do  
            if (first.ge.last) exit
            mid = last
            midval = abs(x(ix(mid)))
            do j = last-1,first,-1
              if (abs(x(ix(j))).gt.midval) then
                mid = mid - 1
!            -- swap
                itemp = ix(mid)
                ix(mid) = ix(j)
                ix(j) = itemp
              end if
            end do
!        -- final swap
            itemp = ix(mid)
            ix(mid) = ix(last)
            ix(last) = itemp
!        -- final test
            if (mid.gt.ncut) then
              last = mid - 1
            else if (mid.lt.ncut) then
              first = mid + 1
            else
              exit
            end if
          end do

      end subroutine split2
!****************************************************************************
!   ascending shellsort of an integer array.

      subroutine srtshai2(n,ix)

      integer, intent(in) :: n ! size of array to be sorted
      integer, intent(inout) :: ix(n) ! array to be sorted

      integer :: is,la,lt,ls,lls,i,j,js,nhalf,ii
 
        if (n.eq.1) return

        call mylog2(2*n,ii)
        ls = 2**ii - 1
        nhalf = n/2
        do lt = 1,ii
          if (ls.gt.nhalf) then
            ls = ls/2
          else
            lls = n - ls
            do i = 1,lls
              is = i + ls
              la = ix(is)
              j = i
              js = is
 100          continue
              if (la.ge.ix(j)) then
                ix(js) = la
                cycle
              else
                ix(js) = ix(j)
                js = j
                j = j - ls
              end if
              if (j.ge.1) go to 100
              ix(js) = la
            end do
            ls = ls/2
          end if
        end do

      end subroutine srtshai2
!****************************************************************************

!   upper whole part of logarithm with base two.
!   the name is chosen in order to avoid conflicts
!   with some compilers and their intrinsic funcions.

      subroutine mylog2(n,k)

      integer,intent(in) :: n
      integer,intent(out) :: k
      integer :: j

      if (n.eq.0) then
        k = 0
        return
      end if

      k = 0
      j = n
      do 
        if (j.eq.1) exit
        k = k + 1
        j = j/2 + mod(j,2)
      end do

      end subroutine mylog2

!*******************************
!  Check matrix for out of range entries (remove) and duplicates (sum)
!  and order by entries in each column. Lower triangular part
!  of matrix held (anything in upper triangle is out of range).
!  The diagonal must be present.

      subroutine check_matrix(n,ne,ptr,row,val,info)

      integer, intent(in) :: n,ne
      integer, intent(inout) :: ptr(n+1)
      integer, intent(inout) :: row(ne)
      real(wp), intent(inout) :: val(ne)
      type(mi28_info), intent(inout) :: info

      integer :: flag69,i,j,mx_type

      mx_type = 4 ! do not want to insist positive diagonal 
      call mc69_cscl_clean(mx_type,n,n,ptr,row,flag69,val=val,noor=info%oor, &
           ndup=info%dup)


      ! check for warnings/errors and translate to mi28 error flags

      if (flag69.eq.-1) then
         info%flag = MI28_ERROR_ALLOCATION
      else if (flag69.eq.-5 .or. flag69.eq.-6) then
         info%flag = MI28_ERROR_PTR
      else if (flag69.eq.-10 .or. flag69.eq.4 .or. flag69.eq.5) then
         info%flag = MI28_ERROR_MISS_DIAG
      else if (info%dup > 0) then
         info%flag = MI28_WARNING_DUP_IDX
      else if (info%oor > 0) then
         info%flag = MI28_WARNING_OOR_IDX
      end if

      if (info%flag.lt.0) return

      ! check diagonal entries are positive
      do i = 1,n
        j = ptr(i)
        if (val(j).le.zero) then
          info%flag = MI28_WARNING_NEG_DIAG
          exit
        end if
     end do
     
    end subroutine check_matrix

!*******************************
! Permute symmetric matrix (only lower triangle held and want
! lower triangle of permuted matrix). put diagonal entry
! first in each column.

    subroutine permute_matrix(n,ne,ptr,row,val,perm,ptr_perm,row_perm, &
        val_perm,iw)

    integer, intent(in) :: n,ne
    ! original matrix (lower triangle)
    integer, intent(in) :: ptr(n+1)
    integer, intent(in) :: row(ne)
    real(wp), intent(in) :: val(ne)

    integer, intent(in) :: perm(n)

    ! permuted matrix (lower triangle)
    integer, intent(out) :: ptr_perm(n+1)
    integer, intent(out) :: row_perm(ne)
    real(wp), intent(out) :: val_perm(ne)

    integer, intent(out) :: iw(n) !work array

    integer :: i,i2,j,j2,k,k1

      ! compute column counts for permuted matrix (lower triangle only)
      iw(1:n) = 0
      do j = 1,n
        j2 = perm(j)
        do k = ptr(j),ptr(j+1)-1
           i = row(k)
           i2 = perm(i)
           if (i2 < j2) then
              iw(i2) = iw(i2) + 1
           else
              iw(j2) = iw(j2) + 1
           end if
        end do
      end do

      ! set column pointer for permuted matrix
      ptr_perm(1) = 1
      do j = 1,n
         ptr_perm(j+1) = ptr_perm(j) + iw(j)
         ! set iw(j) to point to one entry past first entry in col j
         iw(j) = ptr_perm(j) + 1
      end do

      ! drop the entries into permuted matrix, with diag first in col.
      do j = 1,n
        j2 = perm(j)
        do k = ptr(j),ptr(j+1)-1
           i = row(k)
           i2 = perm(i)
           if (i2 == j2) then
              k1 = ptr_perm(j2)
              row_perm(k1) = j2
           else if (i2 < j2) then
              k1 = iw(i2)
              iw(i2) = iw(i2) + 1
              row_perm(k1) = j2
           else
              k1 = iw(j2)
              iw(j2) = iw(j2) + 1
              row_perm(k1) = i2
           end if
           val_perm(k1) = val(k)
         end do
       end do

      ! we now have row_perm, ptr_perm and val_perm holding
      ! lower triangular part of permuted matrix.
    end subroutine permute_matrix

!*******************************

      subroutine compute_scaling(n,ia,ja,aa,iscale,scalesj,st,info)

      integer, intent(in) :: n  ! order of matrix
      integer, intent(in) ::  ia(n+1) ! column pointers for lower triangle.
      integer, intent(in) ::  ja(:) ! row indices of lower triangular part.
        ! diagonal entry must be first entry in each col. 
      real(wp), intent(in) ::  aa(:) ! entries of lower triangular part
      integer, intent(in) :: iscale ! controls choice of scaling
      real(wp), intent(out) :: scalesj(n)
      integer, intent(out) :: st
      type(mi28_info), intent(inout) :: info    ! See derived-type declaration

      integer :: i
      integer :: j, jstrt, jstop
      integer :: k
      real(wp) :: temp,temp2
      logical :: sing

      st = 0
      if (iscale.eq.1) then
        ! compute l2 scaling factors as in Lin and More.
        ! Remember: only lower triangular part of A is held.
        scalesj(1:n) = zero
        do i = 1,n
          jstrt = ia(i)
          jstop = ia(i+1) - 1
          do j = jstrt,jstop
            k = ja(j)
            temp = aa(j)
            temp2 = temp*temp
            scalesj(k) = scalesj(k) + temp2
            if (k.eq.i) cycle
            scalesj(i) = scalesj(i) + temp2
          end do
        end do

        do i = 1,n
          scalesj(i) = sqrt(scalesj(i))
        end do

        do i = 1,n
           if (scalesj(i).gt.zero) then
              scalesj(i) = one/sqrt(scalesj(i))
           else
              scalesj(i) = one
           end if
        end do

      else if (iscale.eq.2) then
         call mc77_scale(n, ia, ja, aa, scalesj, info%flag77, st)

         if (st.ne.0) return

         if (info%flag77 < 0) info%flag = MI28_ERROR_MC77

      else if (iscale.eq.3) then
         call mc64_scale(n, ia, ja, aa, scalesj, st, sing, info%flag64)

         if (st.ne.0) return

         if (sing) then
            info%flag = MI28_ERROR_SINGULAR
            return
         end if

         if (info%flag64 < 0) info%flag = MI28_ERROR_MC64

         if (info%flag64 == 2) info%flag = MI28_WARNING_MC64

      else if (iscale.eq.4) then
         ! diagonal scaling
         do i = 1,n
            jstrt = ia(i)
            scalesj(i) = aa(jstrt)
            if (scalesj(i).gt.zero) then
              scalesj(i) = one/sqrt(scalesj(i))
            else
              scalesj(i) = one
            end if
         end do

      end if

      end subroutine compute_scaling
!*******************************
!
! Following Ruiz and Ucar
! "A symmetry preserving algorithm for matrix scaling"
! we do one iteration of the inf norm and then 3 of the one norm.
!
subroutine mc77_scale(n, ptr, row, val, scalesj, info, st)

   integer, intent(in) :: n ! order of system
   integer, intent(in) :: ptr(n+1) ! column pointers of A
   integer, intent(in) :: row(*) ! row indices of A (lower triangle)
   real(wp), intent(in) :: val(*) ! entries of A (in same order as in row).
   real(wp), dimension(n), intent(out) :: scalesj
   integer, intent(out) :: info
   integer, intent(out) :: st

   integer :: ne, i, j, k

   integer, allocatable, dimension(:) :: iw
   real(wp), allocatable, dimension(:) :: dw
   real(wp), dimension(:), allocatable :: val2

   integer :: icntl77(10), info77(10)
   real(wp) :: cntl77(10), rinfo77(10)

   ! Take absolute value of matrix to avoid overheads
   ne = ptr(n+1)-1
   allocate(val2(ne), iw(2*n), dw(2*n), stat=st)
   if (st .ne. 0) return

   val2(1:ne) = abs(val(1:ne))

   ! Set controls
   call mc77id(icntl77, cntl77)
   icntl77(1) = -1 ! error messages
   icntl77(2) = -1 ! warning messages
   icntl77(3) = -1 ! diagnostic messages
   icntl77(4) = -1 ! disable checking
   icntl77(5) = 1 ! absolute value precomputed
   icntl77(6) = -1 ! symmetric matrix

   ! Single iteration of inf norm
   icntl77(7) = 1 ! max number of iterations
   call mc77ad(0, n, n, ne, ptr, row, val2, iw, size(iw), &
      dw, size(dw), icntl77, cntl77, info77, rinfo77)

   info = info77(1)
   if (info < 0) return

   do i = 1,n
      scalesj(i) = 1/dw(i)
   end do

   ! Apply scaling
   ! Note: we could modify mc77 to take an input vector of scaling to make
   ! this step unnecessary
   do i = 1,n
      do j = ptr(i), ptr(i+1)-1
         k = row(j)
         val2(j) = scalesj(i) * val2(j) * scalesj(k)
      end do
   end do

   ! Up to 3 iterations of one norm
   icntl77(7) = 3 ! max number of iterations
   call mc77ad(1, n, n, ptr(n+1)-1, ptr, row, val2, iw, size(iw), &
      dw, size(dw), icntl77, cntl77, info77, rinfo77)

   info = info77(1)
   if (info < 0) return

   do i = 1,n
      scalesj(i) = scalesj(i) / dw(i)
   end do
   
end subroutine mc77_scale

!**************************************************************
!
! Call mc64 to get a scaling
!
subroutine mc64_scale(n, ptr, row, val, scalesj, st, sing, flag)

   integer, intent(in) :: n ! order of system
   integer, intent(in) :: ptr(n+1) ! column pointers of A
   integer, intent(in) :: row(*) ! row indices of A (lower triangle)
   real(wp), intent(in) :: val(*) ! entries of A (in same order as in row).
   real(wp), dimension(n), intent(out) :: scalesj
   integer, intent(out) :: st ! stat parameter
   logical, intent(out) :: sing ! set to true if matrix found to be singular
   integer, intent(out) :: flag ! error flag from mc64

   type(mc64_control) :: control64
   type(mc64_info) :: info64

   integer :: i
   integer, dimension(:), allocatable :: perm64
   real(wp), dimension(:), allocatable :: cscale
   
   allocate(perm64(2*n), cscale(2*n), stat=st)
   if (st .ne. 0) return
   sing = .false.

   control64%checking = 1 ! checking disabled
   control64%ldiag = -1 ! printing disabled
   call mc64_matching(5, 3, n, n, ptr, row, val, &
      control64, info64, perm64, scale=cscale)

   flag = info64%flag
   select case(flag)
   case(0)
      ! success; do nothing
   case(1)
      ! structually singular matrix
      sing = .true.
      return
   case(2)
      ! large scaling factors; do nothing (don't scale)
      do i = 1, n
         scalesj(i) = one
      end do
      return
   case(-5)
      ! allocation error
      st = info64%stat
      return
   case default
      ! some other error; should never happen
      ! signal through stat=-99
      st = -99
      return
   end select

   do i = 1, n
      scalesj(i) = exp(cscale(i))
   end do

end subroutine mc64_scale

!****************************************************************************
end module hsl_mi28_double
     
