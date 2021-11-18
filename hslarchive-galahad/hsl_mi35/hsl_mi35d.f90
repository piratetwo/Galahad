! COPYRIGHT (c) 2015 Science and Technology Facilities Council,
!                      Academy of Sciences of Czech Republic.
! Version 1.2.0
! Original date April 2015

! Written by:  Jennifer Scott    STFC Rutherford Appleton Laboratory
!              Miroslav Tuma     Academy of Sciences of Czech Republic

! Incomplete factorization preconditioner for sparse least squares problems

! This is based on hsl_mi28 but is for least squares so
! that it computes an incomplete factorization of  C = A^T* W^2 * A
! where W is a diagonal matrix of weights.
! Two options:
! (a) The user forms and inputs C.
! (b) Forming and storing C explicitly is avoided
! (although, if ordering is required, the pattern of C does
! have to be computed and stored)
! Using the l2 scaling, we can also avoid forming C but do have
! to form each col of C in turn so expensive. Other scalings not offered
! for option (b).

! Note: for (b), as C is not held explicitly, if a shift is needed during the 
! factorization, the restart code is more expensive than in (a)
! as the columns of C will have to be recomputed. Possibly remedy: store
! as many cols of C as we can to avoid recomputing ... but we don't
! do this at present.

! Note also that lsize has a (slightly) different meaning from hsl_mi28: here
! is the max. number of entries in a column of the incomplete factor
! (in hsl_mi28 it was the maximum number of fill entries in a column but here
! a column of C could be close to dense and we want L to have sparse columns)
!

! To convert from double to single and visa versa:
! * Change wp
! * Change _double
! * Change BLAS calls: dcopy
! * Change HSL calls: mc77id, mc77ad, mc61id, mc61ad
! * Change default for control%small from 10**(-20) to 10**(-12)
!
module hsl_mi35_double

   use hsl_mc64_double
   use hsl_mc68_integer

   implicit none

   private
   public :: mi35_keep, mi35_control, mi35_info
   public :: mi35_factorize, mi35_finalise, mi35_precondition, mi35_solve
   public :: mi35_check_matrix, mi35_factorizeC, mi35_formC

  integer, parameter  :: wp = kind(0.0d0)
  integer, parameter  :: long = selected_int_kind(18)
  real(wp), parameter :: zero = 0.0_wp
  real(wp), parameter :: one = 1.0_wp
  real(wp), parameter :: sfact = 2.0_wp
  real(wp), parameter :: sfact2 = 4.0_wp
  real(wp), parameter :: alpham = 0.001_wp

   ! Error flags
   integer, parameter :: mi35_ERROR_ALLOCATION      = -1
   integer, parameter :: mi35_ERROR_ROW_TOO_SMALL   = -2
   integer, parameter :: mi35_ERROR_VAL_TOO_SMALL   = -3
   integer, parameter :: mi35_ERROR_N_OOR           = -4
   integer, parameter :: mi35_ERROR_PTR             = -5
   integer, parameter :: mi35_error_all_null_cols   = -6
   integer, parameter :: mi35_ERROR_MC77            = -7
   integer, parameter :: mi35_ERROR_MC64            = -8
   integer, parameter :: mi35_ERROR_SCALE           = -10
   integer, parameter :: mi35_ERROR_USER_PERM       = -11
   integer, parameter :: mi35_ERROR_MC61            = -12
   integer, parameter :: mi35_ERROR_MC68            = -13
   integer, parameter :: mi35_ERROR_DEALLOCATION    = -14
   integer, parameter :: mi35_error_colC            = -15
   integer, parameter :: mi35_error_nzc             = -16
   integer, parameter :: mi35_error_rowA            = -17

   ! warning flags
   integer, parameter :: mi35_warning_OOR_IDX        = 1
   integer, parameter :: mi35_warning_DUP_IDX        = 2
   integer, parameter :: mi35_warning_NZERO          = 3
   integer, parameter :: mi35_warning_null_rows      = 4
   integer, parameter :: mi35_warning_null_cols      = 5
   integer, parameter :: mi35_warning_null_weights   = 6
   integer, parameter :: mi35_warning_order          = 7
   integer, parameter :: mi35_WARNING_MC64           = 8
   integer, parameter :: mi35_WARNING_MC61           = 9



  interface mi35_check_matrix
      module procedure mi35_check_matrix_double
  end interface

  interface mi35_factorize
      module procedure mi35_factorize_double
  end interface

  interface mi35_formC
      module procedure mi35_formC_double
  end interface

  interface mi35_factorizeC
      module procedure mi35_factorizeC_double
  end interface

  interface mi35_precondition
      module procedure mi35_precondition_double
  end interface

  interface mi35_solve
      module procedure mi35_solve_double
  end interface

  interface mi35_finalise
      module procedure mi35_finalise_double
  end interface

!****************************************************************************
  type mi35_keep ! The scalar keep holds the preconditioner, scaling and pivot
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
      !  S Q^T C Q S   is factorized.
      ! To get S, we permute the entries of scale at end of factorize. 
   integer, allocatable :: invp(:) ! specifies the pivot order such that invp(j)
      ! holds the j-th pivot column.
   integer, allocatable :: perm(:) ! specifies the pivot order such that 
      ! perm(i) holds the position of i in the pivot sequence (perm and invp
      ! are inverses of each other) 
   real(wp), allocatable :: w(:) ! Allocated size n if control%iorder>0.

  end type mi35_keep

!****************************************************************************
  type mi35_control ! The scalar control of this type controls the action

    real(wp) :: alpha = zero ! initial shift
    integer :: iorder = 6 ! controls ordering of C. Options:
!         <=0  no ordering
!           1  RCM
!           2  AMD
!           3  user-supplied ordering
!           4  ascending degree
!           5  Metis
!          >= 6, Sloan (MC61)
     integer :: iscale = 4 ! controls whether scaling is used.
        ! iscale <= 0 no scaling
        ! iscale = 1 is Lin and More scaling (l2 scaling)
            ! If C is not formed explicitly, this is potentially expensive 
            ! since it involves computing each col. of C.
        ! iscale = 2 is mc77 scaling ONLY AVAILABLE if C is formed.
        ! iscale = 3 is mc64 scaling ONLY AVAILABLE if C is formed.
        ! iscale = 4 is diagonal scaling
        ! iscale = 5 user-supplied scaling
        ! For all other iscale, use default
     integer :: limit_rowA = -1
       ! If a row of A has more than limit_rowA entries, computation terminates
       ! If negative, no limit (or at least, limit is n).
     integer :: limit_colC = -1
       ! If a col of C has more than limit_colC entries, computation terminates
       ! If negative, no limit (or at least, limit is n).
     integer :: limit_C = -1
       ! If C (upper and lower triangular parts)
       ! has more than limit_C entries, computation terminates.
       ! If negative, no limit (or at least, limit is huge(1)). 
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

  end type mi35_control

!****************************************************************************

  type mi35_info ! The scalar info of this type returns information to user.

      real(wp) :: avlenC = zero ! average number of entries in a row/col of C
      integer :: band_after = 0 ! semibandwidth after MC61
      integer :: band_before = 0 ! semibandwidth before MC61
      integer :: dup = 0 ! number of duplicated entries found in array row.
      integer :: flag = 0 ! error flag
      integer :: flag61 = 0 ! error flag from mc61. Note: if error found, we
        ! just use original order (this includes if failure of allocation of
        ! memory for calling mc61)
      integer :: flag64 = 0 ! error flag from hsl_mc64
      integer :: flag68 = 0 ! error flag from hsl_mc68. Note: if error found, we
        ! just use original order.
      integer :: flag77 = 0 ! error flag from mc77
      integer :: maxlen = 0 ! largest number of entries in a row of A
      integer :: maxlenC = 0 ! largest number of entries in a row/col of C
      integer :: nrestart = 0 ! number of restarts (after reducing the shift)
      integer :: nshift = 0 ! number of non-zero shifts used
      integer :: nnz_C = 0 ! number of entries in lower triangular part of C
      integer :: nzero = 0 ! number of explicit zero entries found and removed.
      integer :: nzero_weight = 0 ! number of zero weights (which are removed).
      integer :: oor = 0 ! number of out-of-range entries found in array row.
      real(wp) :: profile_before = 0 ! semibandwidth before MC61
      real(wp) :: profile_after = 0 ! semibandwidth after MC61
      integer(long) :: size_r = 0_long ! size of arrays jr and ar that are used 
         ! for r   
      integer :: stat = 0 ! Fortran stat parameter
      real(wp) :: alpha = zero ! on successful exit, holds shift used

  end type mi35_info

!****************************************************************************
contains

!  A is m x n (m >= n) sparse matrix in CSC format.
!  weight holds (optional) weights.
!  Check matrix A for out of range entries (remove) and duplicates (sum).
!  Remove any entries with value zero.
!  Also remove any null rows/cols (this is essential for hsl_mi35).
!  Also remove rows corresponding to zero weights.
!  This routine is OPTIONAL but behaviour of factorization routine not
!  guaranteed if this routine has not been called.
!  
      subroutine mi35_check_matrix_double(m,n,ptr,row,val,control,info,weight,b)

      ! input matrix. 
      integer, intent(inout) :: m  ! number of rows in matrix A. On exit, 
         ! row dimension after removal of null rows
         ! and rows corresponding to zero weights
      integer, intent(inout) :: n  ! number of cols in matrix A.  On exit, 
         ! row dimension after removal of cols rows
      integer, intent(inout) ::  ptr(n+1) ! column pointers. May be changed
      integer, intent(inout) ::  row(:) ! row indices. May be changed
      real(wp), intent(inout) ::  val(:) ! entries of matrix. May be changed
         ! (zeros are removed)
      type(mi35_info), intent(out) :: info
      type(mi35_control), intent(in) :: control
      real(wp), intent(inout), optional :: weight(m) ! if present, must hold 
         ! weights. Rows of A corresponding to zero weights are removed and zero
         ! weights are removed from weight.
      real(wp), intent(inout), optional :: b(m) ! if present and null rows are
         ! removed from A, the corrresponding entries are removed from b.

      character(50) :: context ! Procedure name (used when printing).
      integer :: i,j,k,k1,k2,kne
      integer :: m_new ! row dimension after removal of null rows
      integer :: ne,nout,nout1
      integer :: n_new ! col dimension after removal of null cols
      integer :: noor,ndup,nzero,nzj
      integer :: st
      integer, allocatable :: lenrow(:), iw(:) ! both allocated to have size m

   info%avlenC = zero
   info%alpha = max(zero,control%alpha)
   info%band_before = 0
   info%band_after = 0
   info%dup = 0
   info%flag = 0
   info%flag61 = 0
   info%flag68 = 0
   info%maxlen = 0
   info%nrestart = 0
   info%nshift = 0
   info%nnz_C = 0
   info%nzero = 0
   info%nzero_weight = 0
   info%oor = 0
   info%profile_before = zero
   info%profile_after = zero
   info%size_r = 0_long
   info%stat = 0

      context = 'mi35_check_matrix'

      ! Set stream numbers
      nout = control%unit_error
      nout1 = control%unit_warning

      ! Simple check on matrix data
      if (n < 1) then
         info%flag = mi35_ERROR_n_OOR
         call mi35_print_flag(context,nout,info%flag)
         return
      end if

      ne = ptr(n+1) - 1
      if (ne < 1) then
         info%flag = mi35_ERROR_PTR
         call mi35_print_flag(context,nout,info%flag)
         return
      end if

      if (size(row).lt.ne) then
         info%flag = mi35_ERROR_ROW_TOO_SMALL
         call mi35_print_flag(context,nout,info%flag)
         return
      end if

      if (size(val).lt.ne) then
         info%flag = mi35_ERROR_VAL_TOO_SMALL
         call mi35_print_flag(context,nout,info%flag)
         return
      end if

      allocate(iw(m),lenrow(m),stat=st)
      if (st /= 0) then
         info%stat = st
         info%flag = mi35_ERROR_ALLOCATION
         return
      end if
      iw(:) = 0
      lenrow(:) = 0

      ! identify out of range entries and duplicates and count
      ! the number of non-null cols.
        ndup = 0; noor = 0; nzero = 0; kne = 0
        k1 = ptr(1)
        n_new = 0
        do j = 1, n
           k2 = ptr(j+1)
           if (k1.gt.k2) then
              ! ptr is not monotonic
              info%flag = mi35_ERROR_PTR
              return
           else if (k1.eq.k2) then
              ! null col.
              cycle
           end if
           nzj = kne
           do k = k1, k2-1
              i = row(k)
              if (val(k).eq.zero) then
                 ! exclude entry if it is zero
                 nzero = nzero + 1
              else if (i.gt.m .or. i.lt.1) then
                 ! Out-of-range entry
                 noor = noor + 1
              else if (iw(i).le.nzj) then
                 kne = kne + 1
                 row(kne) = i
                 val(kne) = val(k)
                 iw(i) = kne
                 lenrow(i) = lenrow(i) + 1
              else
                ! We have a duplicate in column j
                ndup = ndup + 1
                val(iw(i)) = val(iw(i)) + val(k)
              end if
           end do
           if (kne.eq.nzj) then
              ! no non zero entries in the column so exclude column
              k1 = k2
              cycle
           end if
           n_new = n_new + 1
           ptr(n_new+1) = kne + 1
           k1 = k2
        end do
        info%oor = noor
        info%dup = ndup
        info%nzero = nzero

        if (n_new.eq.0) then
           info%flag = mi35_error_all_null_cols
           call mi35_print_flag(context,nout,info%flag)
           return
        end if

        info%maxlen = maxval(lenrow)
        if (control%limit_rowA.gt.0) then
           if (info%maxlen.gt.control%limit_rowA) then
              info%flag = mi35_error_rowA
              call mi35_print_flag(context,nout,info%flag)
              return
           end if
        end if

        ! the number of cols is now n_new.

        ! if iw(i) = 0 then there was no entry in row i.
        ! check also for zero weights.
        if (present(weight)) then
           do i = 1,m
              if (weight(i).eq.zero) then
                 iw(i) = 0
                 info%nzero_weight = info%nzero_weight + 1
              end if
           end do
        end if

      ! Squeeze out null rows and rows corresponding to zero weights.
      m_new = 0
      if (present(weight)) then
        if (present(b)) then
           do i = 1,m
              k = iw(i)
              if (k.ne.0) then
                 m_new = m_new + 1
                 iw(i) = m_new
                 b(m_new) = b(i)
                 weight(m_new) = weight(i)
              end if
           end do
        else
           do i = 1,m
              k = iw(i)
              if (k.ne.0) then
                 m_new = m_new + 1
                 iw(i) = m_new
                 weight(m_new) = weight(i)
              end if
           end do
        end if
      else
        if (present(b)) then
           do i = 1,m
              k = iw(i)
              if (k.ne.0) then
                 m_new = m_new + 1
                 iw(i) = m_new
                 b(m_new) = b(i)
              end if
           end do
        else
           do i = 1,m
              k = iw(i)
              if (k.ne.0) then
                 m_new = m_new + 1
                 iw(i) = m_new
              end if
           end do
        end if
      end if

        ! update the array row to use new indices
        do k = 1,kne
          row(k) = iw(row(k))
        end do

      ! number of non null rows is m_new and iw(i) holds new index for row i

      ! Simple checks on new matrix data
      if (n_new < 1 .or. m_new < n_new) then
         info%flag = mi35_ERROR_n_OOR
         call mi35_print_flag(context,nout,info%flag)
         n = n_new
         m = m_new
         return
      end if

       ! set warning flags (higher numbers overwrite lower ones)
       if (info%oor.gt.0) then
         info%flag = mi35_warning_OOR_IDX
         call mi35_print_flag(context,nout1,info%flag)
       end if
       if (info%dup.gt.0) then
        info%flag = mi35_warning_DUP_IDX
          call mi35_print_flag(context,nout1,info%flag)
       end if
       if (info%nzero.gt.0) then
        info%flag = mi35_warning_NZERO
          call mi35_print_flag(context,nout1,info%flag)
       end if
       if (m_new.ne.m) then
          info%flag = mi35_warning_null_rows
          call mi35_print_flag(context,nout1,info%flag)
       end if
       if (n_new.ne.n) then
         info%flag = mi35_warning_null_cols
         call mi35_print_flag(context,nout1,info%flag)
       end if
       if (info%nzero_weight.gt.0) then
         info%flag = mi35_warning_null_weights
         call mi35_print_flag(context,nout1,info%flag)
       end if

       ! overwrite original row/col
       n = n_new
       m = m_new

    end subroutine mi35_check_matrix_double

!****************************************************************************

! Incomplete factorization of C = A^T *W^2 * A. User supplies matrix A in CSC
! format; the code generates A^T but no
! checks are made on A or W. C is only formed one col at a time (although
! pattern formed if reordering requested); this limits memory requirements
! (but could be more expensive in terms of time, esp. if restarts needed)

   subroutine mi35_factorize_double(m, n, ptr, row, val, lsize, rsize, keep,  &
      control, info, weight, scale, perm)

   ! input matrix
   integer, intent(in) :: m  ! rows in matrix A
   integer, intent(in) :: n  ! cols in matrix A
   integer, intent(in) ::  ptr(n+1) ! column pointers for A 
   integer, intent(in) ::  row(:) ! row indices for A.
   real(wp), intent(in) ::  val(:) ! entries of A in CSC format.
   integer, intent(in) :: lsize ! if lsize > 0, max. number of entries in a
     ! column of factor is lsize.
   integer, intent(in) :: rsize ! if rsize > 0, max. number of entries in
     ! r is rsize. If rsize = 0, no R is computed.

   type(mi35_keep), intent(out) :: keep ! See derived-type declaration
   type(mi35_control), intent(in) :: control ! See derived-type declaration
   type(mi35_info), intent(out) :: info      ! See derived-type declaration
   real(wp), intent(in), optional :: weight(m) ! Non-zero weights (entries of W)
   real(wp), intent(in), optional :: scale(n) ! Must be present
      ! and hold scaling factors if control%scale =5. 
   integer, intent(in), optional :: perm(n) ! Must be present and 
      ! hold pivot order if control%order = 3.  
      ! perm(i) must hold the position of i in the pivot sequence
   
   !!!!!!!!!!!!!!!!!!!!!
   real(wp) :: alpha_in
   character(50)  :: context      ! Procedure name (used when printing).
   integer :: flag
   integer :: i
   integer :: iorder ! controls ordering
   integer :: iscale ! controls scaling
   integer :: j, jj
   integer :: jm  ! set to 0 if control%rrt = 0 and to 2 otherwise
   integer :: k
   real(wp) :: lowalpha
   integer :: limit(2) ! used to hold control%limit_colC and limit_C
   integer :: liw
   integer :: nout ! stream for errors
   integer :: ne ! set to ptr(n+1) - 1
   integer :: nout1 ! stream for warnings
   integer :: nzc ! entries in (lower triangle of) C (pattern)
   integer :: st ! stat parameter
   real(wp) :: shift_factor
   real(wp) :: shift_factor2
   real(wp) :: tau1 ! set to control%tau1
   real(wp) :: tau2 ! set to control%tau2 but set to zero if
      ! control%tau2 > control%tau1.
   integer :: job ! job parameter for mc61 (controls ordering algorithm)
   integer :: icntl61(10) ! controls for MC61
   integer :: info61(10) ! information from MC61
   real(wp) :: cntl61(5)
   real(wp) :: rinfo61(15)

   type(mc68_control) :: control68
   type(mc68_info) :: info68

   integer, allocatable, dimension(:) :: iw,iw1,mask ! work arrays
   integer, allocatable, dimension(:) :: rowptr,col! holds A in CSR format
   integer, allocatable, dimension(:) :: row_perm,ptr_perm
   integer, allocatable, dimension(:) :: rowc,ptrc
   real(wp), allocatable, dimension(:) :: val_perm
   real(wp), allocatable, dimension(:) :: val_AT !values for W^2*A in CSR format
       ! (= A^T*W^2 in CSC format)

   !!!!!!!!!!!!!!!!!!!!!

   context = 'mi35_factorize'

   ! Set stream numbers
   nout = control%unit_error
   nout1 = control%unit_warning

! initialise data held in info.
   info%avlenC = zero
   info%alpha = max(zero,control%alpha)
   info%band_before = 0
   info%band_after = 0
   info%flag = 0
   info%flag61 = 0
   info%flag68 = 0
   info%maxlenC = 0
   info%nrestart = 0
   info%nshift = 0
   info%profile_before = zero
   info%profile_after = zero
   info%size_r = 0_long
   info%stat = 0

   limit(1) = min(n,control%limit_colC)
   if (limit(1).lt.0) limit(1) = n
   limit(2) = control%limit_C
   if (limit(2).lt.0) limit(2) = huge(1)

   ne = ptr(n+1) - 1

   ! allocate keep%fact_ptr
   deallocate (keep%fact_ptr,stat=st)
   allocate (keep%fact_ptr(n+1),stat=st)
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
        info%flag = mi35_ERROR_USER_PERM
        call mi35_print_flag(context,nout,info%flag)
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
        info%flag = mi35_ERROR_USER_PERM
        call mi35_print_flag(context,nout,info%flag)
        return
     end if
     keep%perm(1:n) = perm(1:n)
   end if

  ! Set rowptr, colA, val_out to hold W^2 * A in CSR format 
!   (= A^T * W^2 in CSC format)
   allocate (mask(m),rowptr(m+1),col(ne),val_AT(ne),stat=st)
   if (st.ne.0) go to 100
   mask(1:m) = 0

   do j = 1, n
      do i = ptr(j), ptr(j+1)-1
         k = row(i)
         mask(k) = mask(k) + 1
      end do
   end do

   rowptr(1) = 1
   do i = 1, m
      rowptr(i+1) = rowptr(i) + mask(i)
      mask(i) = rowptr(i)
   end do

   if (present(weight)) then
      do j = 1, n
         do i = ptr(j), ptr(j+1)-1
            k = row(i)
            jj = mask(k)
            col(jj) = j
            val_AT(jj) = weight(k) * weight(k) * val(i)
            mask(k) = jj + 1
         end do
      end do
   else
      do j = 1, n
         do i = ptr(j), ptr(j+1)-1
            k = row(i)
            jj = mask(k)
            col(jj) = j
            val_AT(jj) = val(i)
            mask(k) = jj + 1
         end do
      end do
   end if

   if (iorder.gt.0 .and. iorder.ne.3) then
     ! to compute ordering, we need the pattern of lower triangle of C = A^T*A

       allocate (ptrc(n+1),stat=st)
       if (st.ne.0) go to 100

      ! we use keep%perm just as a work array here
      if (iorder == 1 .or. iorder == 6) then
         ! allow space for whole of C
         call ata_pattern(2,m,n,ptr,row,rowptr,col,limit,nzc,ptrc,rowc,mask,&
                 flag,st,count=keep%perm)
      else
         call ata_pattern(1,m,n,ptr,row,rowptr,col,limit,nzc,ptrc,rowc,mask,&
                 flag,st,count=keep%perm)
      end if

      ! if failed because of allocation or deallocation error (problem
      ! with rowc, which could happen if C very large/dense), then we will
      ! continue computation without reordering ... issue a warning

      if (flag.eq.-1 .or. flag.eq.-2) then
        iorder = 0
        do i = 1,n
          keep%invp(i) = i
          keep%perm(i) = i
        end do
        info%flag = mi35_WARNING_order
        call mi35_print_flag(context,nout1,info%flag)
      else if (flag.eq.-3) then
        info%flag = mi35_error_colC
        call mi35_print_flag(context,nout,info%flag)
        return
      else if (flag.eq.-4) then
         info%flag = mi35_error_nzc
         call mi35_print_flag(context,nout,info%flag)
         return  
      end if

   end if

   if (iorder == 1 .or. iorder == 6) then

       ! Use mc60 for RCM  or Sloan algorithm
       liw = 8*n + 2
       allocate (iw(liw),stat=st)
       if (st.ne.0) go to 100

       call mc61id(icntl61,cntl61)
       icntl61(1:2) = -1 ! suppress errors/warnings
       ! note: we do not handle duplicates or out-of-range ... they
       ! are treated as an error as we have already checked for these
       ! (and if we handle them here we won't be dealing with the reals)

       job = 1
       if (iorder == 1) job = 2 ! RCM
       call mc61ad(job,n,size(rowc),rowc,ptrc,keep%perm,liw,iw,keep%w, &
            icntl61,cntl61,info61,rinfo61)

      deallocate (iw,stat=st)

      ! if there is a problem with ordering, we will retain original order
      if (info61(1) < 0) then 
         info%flag61 = info61(1)
         iorder = 0
 
      else

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

      end if

      if (iorder.eq.0) then
        do i = 1,n
          keep%invp(i) = i
          keep%perm(i) = i
        end do
        info%flag = mi35_WARNING_order
        call mi35_print_flag(context,nout1,info%flag)
      end if

   else if (iorder == 2 .or. iorder == 5) then

      ! AMD or Metis using hsl_mc68
      job = 1
      if (iorder == 5) job = 3 ! Metis

      ! suppress printing
      control68%lp = -1
      control68%mp = -1
      control68%wp = -1

      call mc68_order(job,n,ptrc,rowc,keep%perm,control68,info68)

      ! check for errors. Return if METIS not linked.
      ! Otherwise, if there is an error, we will just use original order

      if (info68%flag .eq. -5) then
         info%flag68 = info68%flag
         info%flag = mi35_ERROR_MC68
         call mi35_print_flag(context,nout,info%flag,flag68=info%flag68)
         return
      end if

      if (info68%flag < 0 .and. info68%flag.ne.-2) then
         info%flag68 = info68%flag
         iorder = 0
         do i = 1,n
            keep%invp(i) = i
            keep%perm(i) = i
         end do
         info%flag = mi35_WARNING_order
         call mi35_print_flag(context,nout1,info%flag)
      end if

    else if (iorder == 4) then
       ! order by ascending degree (used by Edmond Chow).
       ! compute degrees
       allocate (iw1(n),stat=st)
       if (st.ne.0) go to 100
 
       iw1(1:n) = 0
       do i = 1,n
          do j = ptrc(i),ptrc(i+1)-1
             k = rowc(j)
             iw1(i) = iw1(i) + 1
             if (k.ne.i) iw1(k) = iw1(k) + 1
          end do
       end do

       call kb07ai(iw1,n,keep%invp)
       do i = 1,n
          j = keep%invp(i)
          keep%perm(j) = i
       end do

       deallocate(iw1,stat=st)

    end if

    deallocate(rowc,stat=st)
    deallocate(ptrc,stat=st)

    if (iorder >= 1) then
       ! Permute A. We are going to factorize Q^T C Q = (AQ)^T W^2 (AQ) so
       ! we need to compute AQ, which is a column permutation of A.

       ! Allocate space for permuted matrix.
       allocate (row_perm(ne),val_perm(ne),ptr_perm(n+1),stat=st)
       if (st.ne.0) go to 100

! Permutation held in keep%perm. keep%perm(i) holds position of variable i in 
! list of permuted variables. We also want the inverse permutation
       if (iorder.ne.4 .and. iorder.ne.3) then
          do i = 1,n
             j = keep%perm(i)
             keep%invp(j) = i
           end do
       end if

       if (present(weight)) then
          call permute_matrix(m,n,ne,ptr,row,val,rowptr,col,val_AT, &
                  keep%invp,ptr_perm,row_perm,val_perm,mask,weight)
       else
          call permute_matrix(m,n,ne,ptr,row,val,rowptr,col,val_AT, &
                  keep%invp,ptr_perm,row_perm,val_perm,mask)
       end if

    end if

    deallocate (mask,stat=st)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! compute scaling if requested. Only offer Lin and More or diagonal
   ! scaling (as we do not want to form C explicitly and mc64 or mc77 would
   ! require this). Default is digaonal scaling (iscale = 4)

   iscale = control%iscale
   if (iscale.gt.5) iscale = 4 ! 
   if (iscale.eq.2) iscale = 4 ! use default
   if (iscale.eq.3) iscale = 4 ! use default
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
       info%flag = mi35_ERROR_SCALE
       call mi35_print_flag(context,nout,info%flag)
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
     if (present(weight)) then
       call ichsd_ls(m, n, ptr, row, val, rowptr, col, val_AT,            &
          keep%fact_ptr, keep%fact_row, keep%fact_val, tau1, tau2,jm,     &
          control%maxshift, shift_factor, shift_factor2,control%small,    &
          lsize, rsize, iscale, nout, info, limit, alpha_in, lowalpha,    &
          keep%scale,weight)
     else
       call ichsd_ls(m, n, ptr, row, val, rowptr, col, val_AT,            &
          keep%fact_ptr, keep%fact_row, keep%fact_val, tau1, tau2,jm,     &
          control%maxshift, shift_factor, shift_factor2,control%small,    &
          lsize, rsize, iscale, nout, info, limit, alpha_in, lowalpha,    &
          keep%scale)
     end if
   else
     ! we work with permuted matrix ptr_perm,row_perm,val_perm
     if (present(weight)) then
       call ichsd_ls(m, n, ptr_perm, row_perm, val_perm, rowptr, col, val_AT, &
          keep%fact_ptr, keep%fact_row, keep%fact_val, tau1, tau2,jm,     &
          control%maxshift, shift_factor, shift_factor2,control%small,    &
          lsize, rsize, iscale, nout, info, limit, alpha_in, lowalpha,    &
          keep%scale,weight)
     else
       call ichsd_ls(m, n, ptr_perm, row_perm, val_perm, rowptr, col, val_AT, &
          keep%fact_ptr, keep%fact_row, keep%fact_val, tau1, tau2,jm,     &
          control%maxshift, shift_factor, shift_factor2,control%small,    &
          lsize, rsize, iscale, nout, info, limit, alpha_in, lowalpha,    &
          keep%scale)
     end if
   end if

    if (st .ne. 0) go to 100

    if (info%flag .lt. 0) return

    if (iscale.eq.0) then
       deallocate(keep%scale,stat=st)
    else if (iorder.gt.0) then
       ! at this point, keep%scale has scaling factors for the permuted
       ! C. To allow the user to input the computed scaling on
       ! a subsequent call, we need to return scaling factors for C.
       do i = 1,n
          keep%w(i) = keep%scale(keep%perm(i))
       end do
       keep%scale(1:n) = keep%w(1:n)
    end if

    do i = 1,n
       k = int(keep%fact_ptr(i+1)-keep%fact_ptr(i))
       info%maxlenC = max(info%maxlenC,k)
       info%avlenC = info%avlenC + real(k)
    end do
    info%avlenC = info%avlenC/real(n)

100 if (st .ne. 0) then
      info%stat=st
      info%flag = mi35_ERROR_ALLOCATION
      call mi35_print_flag(context,nout,info%flag,st = info%stat)
      return
    end if

   end subroutine mi35_factorize_double

!****************************************************************************
! Compute the entries in the lower triangular part of = A^T *W^2 * A.

! A is m x n, held in CSC format.
! There is no checking of the user's data.

   subroutine mi35_formC_double(m,n,ptrA,rowA,valA,ptrC,rowC,valC,&
              control,info,weight)

   integer, intent(in) :: m  ! rows in A
   integer, intent(in) :: n  ! cols in A
   integer, intent(in) ::  ptrA(n+1) ! column pointers for A 
   integer, intent(in) ::  rowA(:) ! row indices of A
   real(wp), intent(in) ::  valA(:) ! entries of A

   integer, intent(out) ::  ptrC(n+1) ! column pointers for C 
   integer, intent(out), allocatable ::  rowC(:) ! row indices of lower 
     ! triangular part of C. 
   real(wp), intent(out), allocatable ::  valC(:) ! entries of lower 
     ! triangular part of C. 

   type(mi35_control), intent(in) :: control ! See derived-type declaration
   type(mi35_info), intent(inout) :: info      ! See derived-type declaration
   real(wp), optional, intent(in) ::  weight(m) ! weights (not checked)

   ! local arrays
   integer, allocatable  :: count(:) ! work array, length n
   integer, allocatable :: iw(:) ! work array, length m

   integer, allocatable :: rowptr(:) !
   integer, allocatable :: colA(:) !
   real(wp), allocatable :: val_out(:) !   

  ! on input (ptrA, rowA, valA) hold A in CSC format  = A^T in CSR
  ! We need (rowptr, rowA, val_out) to hold A in CSR format

   character(50)  :: context
   integer :: i,i1,i2,ii,j,j1,j2,jcol,jj,k,ne,st
   integer(long) :: lenc
   integer :: limit(2)
   integer :: nout
   integer :: nzc
   real(wp) :: atemp

   context = 'mi35_formC'
   info%flag = 0

   ! Set stream numbers
   nout = control%unit_error

   info%avlenC = zero
   info%maxlenC = 0

   info%stat = 0
   st = 0

   limit(1) = control%limit_colC
   if (limit(1).lt.0) limit(1) = huge(1)
   limit(2) = control%limit_C
   if (limit(2).lt.0) limit(2) = huge(1)

  ! Set (rowptr, colA, val_out) to hold W^2 * A in CSR format 

   ne = ptrA(n+1) - 1
   allocate(iw(m),count(n),rowptr(m+1),colA(ne),val_out(ne),stat=st)
   if (st.ne.0) go to 100

   iw(1:m) = 0
   do j = 1, n
      do i = ptrA(j), ptrA(j+1)-1
         k = rowA(i)
         iw(k) = iw(k) + 1
      end do
   end do

   rowptr(1) = 1
   do i = 1, m
      rowptr(i+1) = rowptr(i) + iw(i)
      iw(i) = rowptr(i)
   end do

   if (present(weight)) then
      do j = 1, n
         do i = ptrA(j), ptrA(j+1)-1
            k = rowA(i)
            jj = iw(k)
            colA(jj) = j
            val_out(jj) = weight(k) * weight(k) * valA(i)
            iw(k) = jj + 1
         end do
      end do
   else
      do j = 1, n
         do i = ptrA(j), ptrA(j+1)-1
            k = rowA(i)
            jj = iw(k)
            colA(jj) = j
            val_out(jj) = valA(i)
            iw(k) = jj + 1
         end do
      end do
   end if

   ! Now we can compute C = A^T *W^2* A (lower triangle only)

   ! First compute how many entries in lower triangular part of C.

   iw(1:m) = 0
   count(1:n) = 0

   nzc = 0
   ptrC(1) = 1
   do j = 1,n
      i1 = ptrA(j)
      i2 = ptrA(j+1) - 1
      do i = i1,i2
         ii = rowA(i)
         j1 = rowptr(ii)
         j2 = rowptr(ii+1) - 1
         do jj = j1,j2
            jcol = colA(jj)
            if (jcol.le.j .and. iw(jcol).lt.j) then
              iw(jcol) = j
              nzc = nzc + 1
              ! by symmetry, there is an entry in col jcol, row j
              if (jcol.ne.j) count(jcol) = count(jcol) + 1
            end if
          end do
       end do
       ptrC(j+1) = nzc + 1

       lenc = ptrC(j+1) - ptrC(j)
       lenc = lenc + count(j) ! this allows us to count the upper triangular 
                              ! part (by symmetry, equal to entries in row j)
       if (lenc.gt.limit(1)) then
         info%flag = mi35_error_colC
         call mi35_print_flag(context,nout,info%flag)
         return
       end if

       ! test number of entries in C (upper and lower parts)
       if (nzc.gt.(limit(2)-n)/2) then
          info%flag = mi35_error_nzc
          call mi35_print_flag(context,nout,info%flag)
          return
       end if

    end do

    info%nnz_C = ptrC(n+1) - 1

   ! now ready to allocate space and put in the entries
     deallocate(rowC,stat=st)
     allocate(rowC(nzc),stat=st)
     if (st.ne.0) go to 100

     deallocate(valC,stat=st)
     allocate(valC(nzc),stat=st)
     if (st.ne.0) go to 100

   ! Perform multiplication C = A^T *W^2* A; only lower triangle of C required
   ! Array count is length n so use here as work array.
   ! matmulC computes C = X*Y where X and Y are in CSR format
   ! X = A^T is held in CSR format using ptrA, rowA, valA
   ! Y = W^2*A is held in CSR format using rowptr, colA, val_out

   call matmulC(n, ptrA, rowA, valA, rowptr, colA, val_out, &
       ptrC, rowC, valC, count)

   ! ensure the diagonal entry is the first in each column of C
   do j = 1,n
      k = ptrC(j+1) - ptrC(j)
      info%maxlenC = max(info%maxlenC,k)
      info%avlenC = info%avlenC + real(k)
      i1 = ptrC(j)
      k = rowC(i1)
      ! cycle if first entry is the diagonal
      !!!!! I think as we have set things up, the diagonal will
      ! always be first entry????
      if (k.eq.j) cycle
      atemp = valC(i1)
      do i = ptrC(j)+1,ptrC(j+1) - 1
         ii = rowC(i)
         if (ii.eq.j) then
            rowC(i1) = ii
            rowC(i) = k
            valC(i1) = valC(i)
            valC(i) = atemp
            exit
         end if
      end do
   end do
   info%avlenC = info%avlenC/real(n)

   return

100 if (st .ne. 0) then
      info%stat=st
      info%flag = mi35_ERROR_ALLOCATION
      call mi35_print_flag(context,nout,info%flag,st = info%stat)
      return
    end if

   end subroutine mi35_formC_double

!****************************************************************************

  ! compute incomplete factorization of C= A^T *W^2* A. User must input lower
  ! triangle of C.

   subroutine mi35_factorizeC_double(n, ptr, row, val, lsize, rsize, keep,  &
      control, info, scale, perm)

   ! input matrix
   integer, intent(in) :: n  ! order of matrix
   integer, intent(in) ::  ptr(n+1) ! column pointers ! 
   integer, intent(in) ::  row(:) ! row indices of lower triangular part.
     ! If no reordering, diagonal entry must be first in col.
   real(wp), intent(in) ::  val(:) ! entries of lower triangular part.
   integer, intent(in) :: lsize ! max. number of entries in a
     ! column of factor is lsize.
   integer, intent(in) :: rsize ! if rsize > 0, max. number of entries in
     ! r is rsize. If rsize = 0, no r is computed.

   type(mi35_keep), intent(out) :: keep ! See derived-type declaration
   type(mi35_control), intent(in) :: control ! See derived-type declaration
   type(mi35_info), intent(inout) :: info      ! See derived-type declaration
   real(wp), intent(in), optional :: scale(n) ! Must be present 
      ! and hold scaling factors if control%scale =5. 
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
   integer :: iscale ! controls scaling. Note: more options available in this
      ! case where C has been formed.
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

   context = 'mi35_factorizeC'

   ! Set stream numbers
   nout = control%unit_error
   nout1 = control%unit_warning

! initialise data held in info.
   info%alpha = max(zero,control%alpha)
   info%band_before = 0
   info%band_after = 0
   info%flag = 0
   info%flag61 = 0
   info%flag64 = 0
   info%flag68 = 0
   info%flag77 = 0
   info%nrestart = 0
   info%nshift = 0
   info%profile_before = zero
   info%profile_after = zero
   info%size_r = 0_long
   info%stat = 0

   ne = ptr(n+1) - 1

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
        info%flag = mi35_ERROR_USER_PERM
        call mi35_print_flag(context,nout,info%flag)
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
        info%flag = mi35_ERROR_USER_PERM
        call mi35_print_flag(context,nout,info%flag)
        return
     end if
     keep%perm(1:n) = perm(1:n)
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
         info%flag = mi35_ERROR_MC61
         info%flag61 = info61(1)
         call mi35_print_flag(context,nout,info%flag,flag61=info%flag61)
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
        info%flag = mi35_WARNING_MC61
        info%flag61 = info61(1)
        call mi35_print_flag(context,nout1,info%flag,flag61=info%flag61)
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

      ! check for errors. Return if METIS not linked.
      ! Otherwise, if there is an error, we will just use original order

      if (info68%flag .eq. -5) then
         ! Metis no available
         info%flag68 = info68%flag
         info%flag = mi35_ERROR_MC68
         call mi35_print_flag(context,nout,info%flag,flag68=info%flag68)
         return
      end if

      if (info68%flag < 0 .and. info68%flag.ne.-2) then
         info%flag68 = info68%flag
         iorder = 0
         do i = 1,n
            keep%invp(i) = i
            keep%perm(i) = i
         end do
         info%flag = mi35_WARNING_order
         call mi35_print_flag(context,nout1,info%flag)
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

      call permute_matrixC(n,ne,ptr,row,val,keep%perm,ptr_perm,row_perm, &
            val_perm,iw1)

   end if

   deallocate(iw1,stat=st)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! compute scaling if requested.

   iscale = control%iscale
   if (iscale.gt.5) iscale = 4 ! use default
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
       info%flag = mi35_ERROR_SCALE
       call mi35_print_flag(context,nout,info%flag)
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
       ! C. To allow the user to input the computed scaling on
       ! a subsequent call, we need to return scaling factors for C.
       do i = 1,n
          keep%w(i) = keep%scale(keep%perm(i))
       end do
       keep%scale(1:n) = keep%w(1:n)
    end if

100 if (st .ne. 0) then
      info%stat=st
      info%flag = mi35_ERROR_ALLOCATION
      call mi35_print_flag(context,nout,info%flag,st = info%stat)
      return
    end if

   end subroutine mi35_factorizeC_double

!****************************************************************************

   subroutine mi35_precondition_double(n, keep, z, y, info)

! Purpose: to perform the preconditioning operation y = (LL^T)^{-1}z
! ie solve LL^Ty = z, where L is incomplete factor held in CSC format.

      integer, intent(in) :: n
      type(mi35_keep), intent(inout) :: keep
      real(wp), intent(in) :: z(n)
      real(wp), intent(out) :: y(n)
      type(mi35_info), intent(inout) :: info

      integer :: i
      integer(long) :: j,jstrt,jstop
      integer :: k
      real(wp) :: temp

      info%flag = 0

      if (allocated(keep%scale) .and. allocated(keep%invp)) then
       ! recall that keep%scale has scaling factors for C so permute
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

   end subroutine mi35_precondition_double

!****************************************************************************

   subroutine mi35_solve_double(trans, n, keep, z, y, info)

! Purpose: to solve Ly = SQ^Tz or L^TS^{-1}Q^Ty = z
! where L is incomplete factor held in CSC format, 
! S is scaling and Q permutation used in computing L.

      logical, intent(in) :: trans !true if L^T system required, false otherwise
      integer, intent(in) :: n
      type(mi35_keep), intent(inout) :: keep
      real(wp), intent(in) :: z(n)
      real(wp), intent(out) :: y(n)
      type(mi35_info), intent(inout) :: info

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

   end subroutine mi35_solve_double

!****************************************************************************

   subroutine mi35_finalise_double(keep, info)

! Purpose: to deallocate components of keep.

      type(mi35_keep), intent(inout) :: keep
      type(mi35_info), intent(inout) :: info

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

  10  if (info%stat.ne.0) info%flag = mi35_ERROR_DEALLOCATION

   end subroutine mi35_finalise_double

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
      integer, intent(in) :: lsize ! max. number of entries in a
        ! column of factor is lsize.
      integer, intent(in) :: rsize ! if rsize > 0, max. number of entries in
        ! r is rsize. If rsize = 0, no r is computed.

      integer, intent(in) :: iscale ! controls scaling. 

      integer, intent(in) :: nout ! stream for errors
      type(mi35_info), intent(inout) :: info    ! See derived-type declaration
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
!      integer :: ne ! number of entries in A (lower triangle)
      integer(long) :: nec ! set to il_copy(n+1)-1
      integer :: st  ! stat parameter
      logical :: success ! set to .true. if we have a successful factorization
         ! held in il_copy, jl_copy, al_copy
      real(wp) :: temp

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      context = 'mi35_factorizeC'

      if (iscale.gt.0 .and. iscale.lt.5) then

         call compute_scaling(n,ia,ja,aa,iscale,scalesj,st,info)

         if (st.ne.0) go to 500
         if (info%flag.lt.0) then
            if (info%flag77.ne.0) &
              call mi35_print_flag(context,nout,info%flag,flag77=info%flag77)
            if (info%flag64.ne.0) &
              call mi35_print_flag(context,nout,info%flag,flag64=info%flag64)
            return
         end if

         if (info%flag.eq.mi35_WARNING_MC64) then
            call mi35_print_flag(context,nout,info%flag)
         end if

      end if

!  -- allocate space for the preconditioner (ensure not already allocated first)
!     Treat lsize < 0 as 0.
      lsize_local = max(0,lsize)
      lsize_local = min(lsize_local,n-1)
      max_p = n + lsize_local*int(n,long)

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
      ! that the digaonal (which is the first entry in the column) is present
      ! and it should not be zero or negative because C = A^T*A and A does not
      ! have null rows/cols).
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
             tau1,tau2,small,jm,wn02,wr01,start,list,st,flag)

        else
           call ictkl3(n,ia,ja,aa,iscale,scalesj,d,il,jl,al,lsize_local,rsize, &
             tau1,tau2,small,jm,wn02,wr01,start,list,st,flag, &
             ir,jr,ar,startr,listr)

        end if

        if (st.ne.0) go to 500

       !  write (*,*) 'flag,alpha',flag,alpha

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
         info%flag = mi35_ERROR_ALLOCATION
         call mi35_print_flag(context,nout,info%flag,st = info%stat)
      end if

      end subroutine ichsd

!****************************************************************************
!
!   incomplete Cholesky decomposition wrapper for normal equations
!   where A and A^T held, but not C= A^T*W^2*A

      subroutine ichsd_ls(m,n,ia,ja,aa,rowptr,col,val_AT,            &
        il,jl,al,tau1,tau2,jm,                                       &
        maxshift,shift_factor,shift_factor2,small,                   &
        lsize,rsize,iscale,nout,info,limit,alpha_in,lowalpha,scalesj,weight)

      ! input matrix A
      integer, intent(in) :: m  ! rows in A
      integer, intent(in) :: n  ! cols in A
      integer, intent(in) ::  ia(n+1) ! column pointers for A.
      integer, intent(in) ::  ja(:) ! row indices of A.
      real(wp), intent(in) :: aa(:) ! entries of A in CSC format
      integer, intent(in) ::  rowptr(m+1) ! row pointers for A.
      integer, intent(in) ::  col(:) ! col indices of A.
      real(wp), intent(in) :: val_AT(:) ! entries of W^2*A in CSR format

      ! computed incomplete factorization of C = S*A^T*W^2*A*S
      integer(long), intent(out) ::  il(n+1) ! column pointers for incomplete 
         ! factor of C 
      integer, intent(out), allocatable ::  jl(:) ! row indices for incomplete 
         ! factor of C
      real(wp), intent(out), allocatable ::  al(:) ! entries of incomplete
         ! factor of C

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
      integer, intent(in) :: lsize ! if lsize > 0, max. number of entries in a
        ! column of factor is lsize.
      integer, intent(in) :: rsize ! if rsize > 0, max. number of entries in
        ! r is rsize. If rsize = 0, no r is computed.

      integer, intent(in) :: iscale ! controls scaling. 

      integer, intent(in) :: nout ! stream for errors
      type(mi35_info), intent(inout) :: info    ! See derived-type declaration
      integer, intent(in) :: limit(2) ! limit on number of entries in C
      real(wp), intent(in) :: alpha_in ! Holds initial shift
      real(wp), intent(inout) :: lowalpha ! Lower bound on shift (gets altered
         ! if shift altered)
      real(wp), intent(inout) :: scalesj(:) ! Holds scaling factors.
      real(wp), optional, intent(in) :: weight(m) ! Holds weights.

      ! local arrays
      integer, allocatable, dimension(:) :: list ! used by ictkl3_ls. 
         ! allocated to have size n+1
      integer(long), allocatable, dimension(:) :: start ! used by ictkl3_ls. 
         ! allocated to have size n+1
      integer, allocatable, dimension(:) :: wn02 ! used by ictkl3_ls. 
         ! allocated to have size n+1
      real(wp), allocatable, dimension(:) :: wr01 ! used by ictkl3_ls. 
         ! allocated to have size n+1

      ! the following are allocated if rsize .ne. 0
      integer, allocatable, dimension(:) :: listr ! used by ictkl3_ls. 
         ! allocated to have size n+1
      integer(long), allocatable, dimension(:) :: startr ! used by ictkl3_ls. 
         ! allocated to have size n+1
      integer(long), allocatable ::  ir(:) ! column pointers for r. allocated to
         ! have size n+1
      integer, allocatable ::  jr(:) ! row indices for r. allocated to
         ! have size max_r 
      real(wp), allocatable ::  ar(:) ! entries of r. allocated to
         ! have size max_r 

      real(wp), allocatable ::  d(:) ! holds the diagonal entries of C
         ! and update as factorization proceeds (detect breakdown asap)
      real(wp), allocatable ::  dcopy(:) ! holds the diagonal entries of C

      ! if we have a shift that gives a successful factorization but we hope
      ! to reduce the shift, we copy the successful incomplete factorization
      ! into the followinf arrays (note: if we can't successfully allocate them
      ! then we must recompute the factorization)
      integer(long), allocatable ::  il_copy(:) 
      integer, allocatable ::  jl_copy(:) 
      real(wp), allocatable ::  al_copy(:)

      integer, allocatable :: mask(:) ! allocate to have size n. 
        ! use in computing entries of C

      ! local scalars
      real(wp) :: alpha ! manteuffel shift
      real(wp) :: alpha_old ! initialise to 0 and set to an alpha that works
         ! if alpha is reduced
      character(50)  :: context      ! Procedure name (used when printing).
     ! real(wp) :: dsfl
      integer :: flag  ! error flag
      integer :: flag_previous  ! error flag using previous shift
      integer :: i
      integer :: lsize_local ! set to max(0,lsize). allows us to handle lsize <0
      integer(long) :: max_p ! Holds size of jl and al.
         ! jl and al will be allocated to have size max_p.
      integer(long) :: max_r ! Holds size of jr and ar.
         ! jr and ar will be allocated to have size max_r
      integer :: j
      integer(long) :: nec ! set to il_copy(n+1)-1
      integer :: st  ! stat parameter
      logical :: success ! set to .true. if we have a successful factorization
         ! held in il_copy, jl_copy, al_copy

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      context = 'mi35_factorize'

!  -- allocate space for the preconditioner (ensure not already allocated first)
!     Treat lsize < 0 as 0.
      lsize_local = max(0,lsize)
      lsize_local = min(lsize_local,n-1)
      max_p = n + lsize_local*int(n,long)

      deallocate (jl,stat=st)
      deallocate (al,stat=st)
      allocate(jl(max_p),al(max_p),stat=st)
      if (st.ne.0) go to 500
 
      allocate(start(n+1),list(n+1),wn02(n+1),wr01(n+1),d(n),dcopy(n), &
               mask(n),stat=st)
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

      ! compute scaling of C, if requested

      if (iscale.eq.1) then

         call scalel2(m,n,ia,ja,aa,rowptr,col,val_AT,scalesj, &
              wr01,wn02,mask,limit,flag)

         if (flag.lt.0) then
            info%flag = flag
            call mi35_print_flag(context,nout,info%flag)
            return
         end if

      end if

      ! compute (scaled) diagonal entries of C
      if (iscale.eq.1 .or. iscale.eq.5) then
         if (present(weight)) then
           call ata_diag(m,n,ia,ja,aa,d,scale=scalesj,weight=weight)
         else
           call ata_diag(m,n,ia,ja,aa,d,scale=scalesj)
         end if

      else
         if (present(weight)) then
           call ata_diag(m,n,ia,ja,aa,d,weight=weight)
         else
           call ata_diag(m,n,ia,ja,aa,d)
         end if
         if (iscale.eq.4) then
            ! diagonal scaling
            do i = 1,n
              scalesj(i) = one/sqrt(d(i))
              d(i) = one
            end do
        end if

      end if
      dcopy(1:n) = d(1:n)

       !  -- get the Manteuffel shift alpha

      ! compute the initial shift. 
      if (alpha_in.ne.zero) then
         alpha = max(alpha_in,zero)
      else
         alpha = zero
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
        ! entry is encountered
        do j = 1,n
          d(j) = dcopy(j) + alpha
        end do


!    -- perform the decomposition
        if (rsize.eq.0) then
           call ictkl3_ls(m,n,ia,ja,aa,rowptr,col,val_AT,info%nnz_C, &
             iscale,scalesj,d,il,jl,al,                              &
             lsize_local,rsize,tau1,tau2,small,jm,wn02,wr01,mask,start,list, &
             limit,flag)

        else
           call ictkl3_ls(m,n,ia,ja,aa,rowptr,col,val_AT,info%nnz_C, &
             iscale,scalesj,d,il,jl,al,                              &
             lsize_local,rsize,tau1,tau2,small,jm,wn02,wr01,mask,start,list,&
             limit,flag,ir,jr,ar,startr,listr)

        end if

!    -- check the output flag and alpha

        !  write (6,*) 'flag,alpha',flag,alpha

           if (flag.lt.0) then
              info%flag = flag
              call mi35_print_flag(context,nout,info%flag)
              return
           end if

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

           else if (flag.gt.0) then

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

         if (flag.eq.0) exit shift1

      end do shift1

      info%alpha = alpha

 500  info%stat=st
      if (info%stat .ne. 0) then
         info%flag = mi35_ERROR_ALLOCATION
         call mi35_print_flag(context,nout,info%flag,st = info%stat)
      end if


      end subroutine ichsd_ls
!****************************************************************************
!
! routine to print errors and warnings
!
subroutine mi35_print_flag(context,nout,flag,st,flag61,flag64,flag68,flag77)
   integer, intent(in) :: flag, nout
   integer, intent(in), optional :: st
   integer, intent(in), optional :: flag61,flag64,flag68,flag77
   character (len = *), intent(in) :: context

   if (nout < 0) return
   if (flag < 0) then
      write (nout,'(/3a,i3)') ' Error return from ',trim(context),&
         '. Error flag = ', flag
   else
      write (nout,'(/3a,i3)') ' Warning from ',trim(context),&
         '. Warning flag = ',flag
   end if

   select case(flag)
   case(mi35_ERROR_ALLOCATION)
      if (present(st)) then
         write (nout,'(a,i6)') ' Allocation error. stat parameter = ', st
      else
         write (nout,'(a)') ' Allocation error'
      end if

   case(mi35_ERROR_ROW_TOO_SMALL)
         write (nout,'(a)') ' Size of array row is less than ptr(n+1)-1'

   case(mi35_ERROR_VAL_TOO_SMALL)
         write (nout,'(a)') ' Size of array val is less than ptr(n+1)-1'

   case(mi35_ERROR_n_OOR)
         write (nout,'(a)') &
         ' After removal of null rows/cols, m > n >= 1 is not satisfied'

   case(mi35_error_all_null_cols)
         write (nout,'(a)') ' All columns of A are null columns'

   case(mi35_ERROR_PTR)
         write (nout,'(a)') ' ptr(n+1)-1 is less than 1 or ptr is not monotonic'

   case(mi35_ERROR_SCALE)
         write (nout,'(a)') ' scale should be present'

   case(mi35_ERROR_USER_PERM)
         write (nout,'(a)') ' Either perm not present when expected'
         write (nout,'(a)') ' or perm does not hold a permutation'

   case(mi35_ERROR_MC77)
      write (nout,'(a,i3)') &
      ' Unexpected error return from MC77 (scaling routine). flag77 = ',flag77

   case(mi35_ERROR_MC64)
      write (nout,'(a,i3)') &
      ' Unexpected error return from HSL_MC64 (scaling routine). flag64 =', &
        flag64

   case(mi35_ERROR_MC61)
      write (nout,'(a,i3)') &
      ' Unexpected error return from  MC61 (RCM or Sloan ordering). flag61 =', &
        flag61

   case(mi35_ERROR_MC68)
      write (nout,'(a,i3)') &
       ' METIS ordering requested but METIS not linked. flag68 =',flag68


   case(mi35_ERROR_rowA)
         write (nout,'(a)') &
        ' One or more rows of A has more than control%limit_rowA entries'

   case(mi35_ERROR_colC)
         write (nout,'(a)') &
        ' One or more columns of C has more than control%limit_colC entries'

   case(mi35_ERROR_nzc)
         write (nout,'(a)') &
        ' C has more than control%limit_C entries'



   case(mi35_WARNING_OOR_IDX)
         write (nout,'(a)') ' out- of-range indices in row removed'

   case(mi35_WARNING_DUP_IDX)
         write (nout,'(a)') &
      ' Duplicated indices in row removed; corresponding entries in val summed'

   case(mi35_WARNING_order)
         write (nout,'(a)') &
       ' Problem with computing an ordering. Original ordering retained'

  case(mi35_WARNING_NZERO)
         write (nout,'(a)') &
      ' Explicit zeros have been removed'

   case(mi35_WARNING_NULL_ROWS)
         write (nout,'(a)') &
      ' Null row(s) found and removed'

   case(mi35_WARNING_NULL_COLS)
         write (nout,'(a)') &
      ' Null col(s) found and removed'

   case(mi35_WARNING_NULL_WEIGHTS)
         write (nout,'(a)') &
      ' Null weight(s) found and removed. Corresponding rows removed.'

   case(mi35_WARNING_MC61)
      write (nout,'(a,i3)') &
      ' Warning error return from  MC61 (Original order retained). flag61 =', &
        flag61

   end select

end subroutine mi35_print_flag

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
        tau1,tau2,small,jm,wn02,wr01,startl,listl,st,info,ir,jr,ar,startr,listr)

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
      integer, intent(in) :: lsize ! max. number of entries in a
        ! column of factor is lsize.
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

      integer, intent(out) :: st ! stat parameter
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
      st = 0

!  -- startl: starting pointers for the fully formed columns
!  -- used also as an auxiliary vector for the fresh column
!  -- listl: linked list 
      startl(1:n) = 0_long
      listl(1:n) = 0

!  -- initialize the loop of computed columns
      saved = 0

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
          ! (note: diagonal entry is up-to-date)
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

!&&&&&&
!          asize = a_jstop - a_jstrt
!          jsize = min(asize+lsize+saved,ind2)
!          saved = saved + asize + lsize - jsize

          jsize = min(lsize+saved,ind2)
          saved = saved + lsize - jsize
!&&&&&&
          !!!! if we set saved = 0, we are not passing on any spare space
          ! eg get spare space in col. 1 since never any fill in col. 1
          ! so that col 1 of L has same pattern as A (independent of lsize).
          ! using saved = 0 allows direct comparison with Lin and More
          ! also for some problems, seem to get worse preconditioner
          ! if saved > 0 allowed (eg DNVS/m_t1)
      !     saved = 0
 
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
             info = -j
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
!&&&&&&
!          asize = a_jstop - a_jstrt
!          jsize = min(asize+lsize+saved,ind2)
!          saved = saved + asize + lsize - jsize

          jsize = min(lsize+saved,ind2)
          saved = saved + lsize - jsize
!&&&&&&

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
!************************************************ 
!  left-looking incomplete Cholesky (ic) preconditioner of C=A^T*W^2*A  by value
!  with positive semidefinite stabilization.
!  sparse left-looking implementation.
!  uses two drop tolerances:
!  1/ tau1: used to select "small" entries that are dropped from L
!  2/ tau2: used to select "tiny" entries that are dropped from R
!  So all entries of L are at least tau1 in absolute value, and those in
!  R are all smaller than those in L and are at least tau2.

      subroutine ictkl3_ls(m,n,ia,ja,aa,rowptr,col,val_AT,nzc,      &
         iscale,scalesj,d,il,jl,al,lsize,rsize,                     &
         tau1,tau2,small,jm,wn02,wr01,mask,startl,listl,limit,info, &
         ir,jr,ar,startr,listr)

      integer, intent(in) :: m  ! rows in A
      integer, intent(in) :: n  ! cols in A

      integer, intent(in) ::  ia(n+1) ! column pointers for A.
      integer, intent(in) ::  ja(:) ! row indices of A.
      real(wp), intent(in) ::  aa(:) ! entries of A in CSC format.
      integer, intent(in) ::  rowptr(m+1) ! col pointers for A.
      integer, intent(in) ::  col(:) ! col indices of A.
      real(wp), intent(in) :: val_AT(:) ! entries of W^2*A in CSR format
      integer, intent(out) :: nzc ! number of entries in C (lower triangle)
      integer, intent(in) :: iscale ! controls scaling
      real(wp), intent(in) :: scalesj(:) ! scaling factors
      real(wp), intent(inout) :: d(n) ! scaled and shifted diagonal entries  
         ! of C (updated as factorization proceeds)
      integer(long), intent(out) :: il(n+1) ! column pointers for incomplete
         ! factor of C. 
      integer, intent(out) :: jl(:) ! row indices for incomplete factor of C. 
      real(wp), intent(out) :: al(:) ! entries of incomplete factor of C. 
      integer, intent(in) :: lsize ! max. number of entries in a
        ! column of factor is lsize.
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
      integer, intent(out) :: mask(n)

      integer, intent(in) :: limit(2) ! limit on number of entries in C
      integer, intent(out) :: info ! error flag that is set to j
        ! if pivot candidate j (diagonal entry of col j) is too small.
        ! It is set neagtive if the limit on number of entries
        ! in a col of C or in total in C is exceeded.
 
      ! These are used in the case rsize > 0 (if rsize = 0
      ! they are never allocated). As they are not needed if rsize = 0,
      ! they have been made optional
      integer, optional, intent(inout) :: listr(n+1) ! linked list of 
         ! columns for r
      integer(long), optional, intent(inout) :: startr(n+1) ! 
      integer(long), optional, intent(inout) :: ir(n+1) ! column pointers for R.
      integer, optional, intent(inout) :: jr(:) ! row indices for R.  
      real(wp), optional, intent(inout) :: ar(:)! values of entries in R.

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
      nzc = n ! allow for diagonal

!  -- startl: starting pointers for the fully formed columns
!  -- used also as an auxiliary vector for the fresh column
!  -- listl: linked list 
      startl(1:n) = 0_long
      listl(1:n) = 0

      mask(1:n) = 0

!  -- initialize the loop of computed columns
      saved = 0

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
          ! (note: diagonal entry is up-to-date)
          ! If necessary, exit for another shift
          diagj = d(j)
          if (diagj.le.small) then
            info = j
            return
          end if

          ! compute col j of C.
          ! Scatter the (scaled) entries into full array wr01 with indices of
          ! entries held in wn02(1:ind2).
          if (iscale.gt.0) then
             call ata_j(j,m,n,ia,ja,aa,rowptr,col,val_AT, &
                  ind2,wr01,wn02,mask,scalesj)
          else
             call ata_j(j,m,n,ia,ja,aa,rowptr,col,val_AT, &
                  ind2,wr01,wn02,mask)
          end if

          ! check ind2 does not exceed limit(1) and that the sum of entries in C
          ! does not exceed limit(2)
          nzc = nzc + ind2
          if (ind2.gt.limit(1)) then
             info = mi35_error_colC
             return
          else if (nzc.gt.(limit(2)-n)/2) then
             info = mi35_error_nzc
             return
          end if

          ! Set startl(i) = 1 for each row i of C that has an entry in col. j
          do k = 1, ind2
            i = wn02(k)
            startl(i) = 1_long
          end do

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

          jsize = min(lsize+saved,ind2)
          saved = saved + lsize - jsize

          !!!! if we set saved = 0, we are not passing on any spare space
          ! eg get spare space in col. 1 since never any fill in col. 1
          ! so that col 1 of L has same pattern as A (independent of lsize).
          ! using saved = 0 allows direct comparison with Lin and More
          ! also for some problems, seem to get worse preconditioner
          ! if saved > 0 allowed (eg DNVS/m_t1)
      !     saved = 0
 
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
                info = j
                return
             end if
             al(lind) = temp
             jl(lind) = k
             lind = lind + 1
          end do

!      -- finalize the loop. 

          if (j.lt.n) then 
          ! set il for next column
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
             info = j
             return
          end if

          ! compute col j of C.
          ! Scatter the (scaled) entries into full array wr01 with indices of
          ! entries held in wn02(1:ind2).

          if (iscale.gt.0) then
             call ata_j(j,m,n,ia,ja,aa,rowptr,col,val_AT, &
                  ind2,wr01,wn02,mask,scalesj)
          else
             call ata_j(j,m,n,ia,ja,aa,rowptr,col,val_AT, &
                  ind2,wr01,wn02,mask)
          end if

          ! check ind2 does not exceed limit(1) and that the sum of entries in C
          ! does not exceed limit(2)
          nzc = nzc + ind2

          if (ind2.gt.limit(1)) then
             info = mi35_error_colC
             return
          else if (nzc.gt.(limit(2)-n)/2) then
             info = mi35_error_nzc
             return
          end if

          ! Set startl(i) = 1 for each row i of C that has an entry in col. j
          do k = 1,ind2
            i = wn02(k)
            startl(i) = 1_long
          end do

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
                    ! entry already non zero so allow update.
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

          jsize = min(lsize+saved,ind2)
          saved = saved + lsize - jsize

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
                       info = j
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
                          info = j
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
            ! set il and jl for next column
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

      end subroutine ictkl3_ls

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
! Column permutation of matrix A.
! A held by cols using row, ptr, val and held by rows using col,rowptr.
! rowptr is not changed and we overwrite col.
! We do not want to change user's data (row,ptr,val) so we write
! permute matrix to row_perm,ptr_perm,val_perm

    subroutine permute_matrix(m,n,ne,ptr,row,val,rowptr,col,val_AT, &
        perm,ptr_perm,row_perm,val_perm,iw,weight)

    integer, intent(in) :: m,n,ne ! original matrix
    integer, intent(in) :: ptr(n+1)
    integer, intent(in) :: row(ne)
    real(wp), intent(in) :: val(ne)

    integer, intent(in) :: rowptr(m+1)
    integer, intent(out) :: col(ne)
    real(wp), intent(out) :: val_AT(ne)

    integer, intent(in) :: perm(n) ! perm(j) holds jth pivot so that if 
      ! k = perm(j), col. j of permuted matrix is col.k of original matrix

    ! permuted matrix
    integer, intent(out) :: ptr_perm(n+1)
    integer, intent(out) :: row_perm(ne)
    real(wp), intent(out) :: val_perm(ne)

    integer, intent(out) :: iw(m) ! work array

    real(wp), optional, intent(in) :: weight(m)


    integer :: i,j,jj,k

      ! Set ptr_perm(j) to point to start of col j in permuted matrix
      ptr_perm(1) = 1
      do j = 1,n
         k = perm(j)
         ptr_perm(j+1) = ptr_perm(j) + (ptr(k+1) - ptr(k))
      end do

      ! drop the entries into permuted matrix.
      do j = 1,n
         k = perm(j)
         row_perm(ptr_perm(j):ptr_perm(j+1)-1) = &
             row(ptr(k):ptr(k+1)-1)
         val_perm(ptr_perm(j):ptr_perm(j+1)-1) = &
             val(ptr(k):ptr(k+1)-1)
       end do

      ! we now have row_perm, ptr_perm and val_perm holding permuted matrix
      ! in CSC format. 

    ! set rowptr, col and val_AT to hold permuted matrix in CSR
    ! format (=A^T*W^2 in CSC)
    iw(1:m) = rowptr(1:m)
    if (present(weight)) then
      do j = 1, n
         do i = ptr_perm(j), ptr_perm(j+1)-1
            k = row_perm(i)
            jj = iw(k)
            col(jj) = j
            val_AT(jj) = weight(k) * weight(k) * val_perm(i)
            iw(k) = jj + 1
         end do
      end do
   else
      do j = 1, n
         do i = ptr_perm(j), ptr_perm(j+1)-1
            k = row_perm(i)
            jj = iw(k)
            col(jj) = j
            val_AT(jj) = val_perm(i)
            iw(k) = jj + 1
         end do
      end do
   end if

    end subroutine permute_matrix

!*******************************
! Permute symmetric matrix (only lower triangle held and want
! lower triangle of permuted matrix). put diagonal entry
! first in each column.

    subroutine permute_matrixC(n,ne,ptr,row,val,perm,ptr_perm,row_perm, &
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

    end subroutine permute_matrixC


!*******************************

! Code to compute the pattern of lower triangle of C= A^T * A, without forming C
! C held in CSC format

! A is m x n, held in CSC format.
! There is no checking of the user's data.

   subroutine ata_pattern(job,m,n,ptr,row,rowptr,col,limit,nzc,ptrc,rowc,mask,&
              flag,stat,count)

   integer, intent(in) :: job ! 
      ! 1 to true if only lower triangular part of C is required
      ! 2 to true if only lower triangular part of C is required but rowc
      !   allocated to be large enough for upper and lower parts
   integer, intent(in) :: m  ! rows in A
   integer, intent(in) :: n  ! cols in A
   integer, intent(in) ::  ptr(n+1) ! column pointers for A 
   integer, intent(in) ::  row(:) ! row indices of A
   integer, intent(out) :: rowptr(m+1) ! row pointers for A
   integer, intent(out) :: col(:) ! col indices of A 
   integer, intent(in) :: limit(2) ! used to terminate the routine if C
     ! found to be too dense. 
     ! If a col of C has more than limit(1) entries, computation terminates
     ! If C has more than limit(2) entries, computation terminates.
     ! If negative, no limit
   integer, intent(out) :: nzc ! no. of entries in (lower triangle of) C
   integer, intent(out) ::  ptrc(n+1) ! column pointers for C 
   integer, intent(out), allocatable :: rowc(:)
   integer, intent(out) :: mask(m) ! work array
   integer, intent(out) :: flag ! error flag
     ! -1 : allocation error
     ! -2 : deallocation error
     ! -3 : a column of C has more than limit(1) entries (if lower = .true. this
     !      is the lower triangular part of C only)
     ! -4 : number of entries in C exceeds limit(2). 
   integer, intent(out) :: stat ! stat parameter
   integer, intent(out), optional :: count(n) ! only used if job = 1,2. In this
     ! case, accumulates counts in rows so that we can get a count for
     ! the number of entries in a col of C (upper and lower parts)

   integer :: i,i1,i2,ii,j,j1,j2,jcol,jj,k
   integer :: lenrowc,lenc

   flag = 0
   stat = 0

  ! Set rowptr, col to hold A in CSR format
   mask(1:m) = 0

   do j = 1, n
      do i = ptr(j), ptr(j+1)-1
         k = row(i)
         mask(k) = mask(k) + 1
      end do
   end do

   rowptr(1) = 1
   do i = 1, m
      rowptr(i+1) = rowptr(i) + mask(i)
      mask(i) = rowptr(i)
   end do

   do j = 1, n
      do i = ptr(j), ptr(j+1)-1
         k = row(i)
         jj = mask(k)
         col(jj) = j
         mask(k) = jj + 1
      end do
   end do

   ! First compute how many entries in (lower triangular part of) C.
   mask(1:m) = 0
   if (present(count)) count(1:n) = 0

   ! loop over cols of C
   nzc = 0
   ptrc(1) = 1
      do j = 1,n
         i1 = ptr(j)
         i2 = ptr(j+1) - 1
         do i = i1,i2
            ii = row(i)
            j1 = rowptr(ii)
            j2 = rowptr(ii+1) - 1
            do jj = j1,j2
               jcol = col(jj)
                if (jcol.le.j .and. mask(jcol).lt.j) then
                   mask(jcol) = j
                   nzc = nzc + 1
                   ! by symmetry, there is an entry in col jcol, row j
                   if (jcol.ne.j) count(jcol) = count(jcol) + 1
                end if
             end do
          end do
          ptrC(j+1) = nzc + 1
          lenc = ptrC(j+1) - ptrC(j)
          lenc = lenc + count(j) ! this allows us to include upper triangular 
                                 ! part (by symmetry, equal to entries in row j)
          if (lenc.gt.limit(1)) then
             flag = -3
             return
          end if

         ! check size of C (upper and lower triangular parts)
         if (nzc.gt.(limit(2)-n)/2) then
           flag = -4
           return
         end if

       end do

   ! reinitialise mask
   mask(1:m) = 0

   lenrowc = nzc
   if (job.eq.2) lenrowc = 2*nzc - n

   ! allocate space for the pattern of C
     allocate(rowc(lenrowc),stat=stat)
     if (stat.ne.0) then
       flag = -1
       return
     end if

   call matmulC_pattern(n, ptr, row, rowptr, col, ptrc, rowc, count)

   end subroutine ata_pattern

!****************************************************************************
   subroutine matmulC_pattern(n,ia,ja,ib,jb,ic,jc,mask)

    ! compute the product of two sparse matrices A=(ia,ja,a) and
    ! B = (ib,jb,b) in CSR format. Output C in CSR format.
    ! Here we compute only the upper triangle of C
    ! since we are interested in A = B^T so that C = B^T * B is symmetric.
    ! Of course, upper triangle of C in CSR format = lower triangle in CSC
    ! format
    !!! note: as we only want the case where C is square, we have 
    ! taken row size of A = col size of B

       ! n: row size of A = column size of B
    integer, intent (in) :: n
       ! ia: row pointer of A
       ! ja: column indices of A
       ! ib: row pointer of B
       ! jb: column indices of B
    integer, intent(in) :: ia(*),ib(*),ja(*),jb(*)
       ! ic: row pointer of C
       ! jc: column indices of C
    integer, intent(out) :: ic(*),jc(*)
       ! a: entry values of A
       ! b: entry values of B
       ! c: entry values of C
       ! working array
    integer, intent(out) :: mask(*) ! length n

    integer :: nz,i,j,k,icol,icol_add,neigh

    ! initialise the mask array which is an array that
    ! has non-zero value if the
    ! column index already exist, in which case
    ! the value is the index of that column

    ! code computes one row of C at a time
    mask(1:n) = 0
    nz = 0
    ic(1) = 1
    do i = 1,n
       do j = ia(i),ia(i+1)-1
          neigh = ja(j)
          do k = ib(neigh),ib(neigh+1)-1
             icol_add = jb(k)

             ! cycle if entry is in lower triangle
             if (icol_add.lt.i) cycle

             icol = mask(icol_add)
             if (icol.eq.0) then
                nz = nz + 1
                jc(nz) = icol_add
                mask(icol_add) = nz     ! add mask
             end if
          end do
       end do
       ! done with row i, so set mask to zero again
       do j = ic(i),nz
          mask(jc(j)) = 0
       end do
       ic(i+1) = nz + 1
    end do

  end subroutine matmulC_pattern 

!****************************************************************************
! Code to compute the entries in the strict upper triangular 
! part of a single row of C= S * A^T *W^2 * A * S, without forming C.
! As C is symmetric, strict upper == strict lower so
! we are equivalently computing entries in strict lower triangular
! part of a column of C.
! Here S is an optional diagonal scaling matrix, supplied by user

! A is m x n, input in CSC format (ptr,row,val) ie A^T in CSR format
! W^2*A input in CSR format (rowptr,col,val_AT).
! There is no checking of the user's data.

   subroutine ata_j(q,m,n,ptr,row,val,rowptr,col,val_AT, &
         nc,wr01,wn02,mask,scale)

   integer, intent(in) :: q ! index of row that is required
   integer, intent(in) :: m  ! rows in A
   integer, intent(in) :: n  ! cols in A
   integer, intent(in) ::  ptr(n+1) ! column pointers for A = row pointers A^T 
   integer, intent(in) ::  row(:) ! row indices of A
   real(wp), intent(in) ::  val(:) ! entries of A in CSC format (=A^T in CSR)
   integer, intent(in) ::  rowptr(m+1) ! row pointers for A 
   integer, intent(in) ::  col(:) ! col indices of A
   real(wp), intent(in) ::  val_AT(:) ! entries of W^2*A in CSR format
   integer, intent(out) :: nc ! no. of entries in strict upper triangular  
     ! part of the sought-after column of C
   real(wp), intent(out) ::  wr01(n) ! on exit, holds strict upper triangular  
     ! part of row of C
   integer, intent(out) ::  wn02(n) ! on exit, holds col indices of entries in 
     ! strict upper triangular part of row of C
   integer, intent(inout) :: mask(n) ! work array that must be set to zero 
     ! on entry and is reset to zero before exit
   real(wp), optional, intent(in) :: scale(n) ! scaling factors 

   integer :: i,icol,icol_add,j,k,kcol,neigh
   real(wp) :: aij

       i = q
       nc = 0
       do j = ptr(i),ptr(i+1)-1
          aij = val(j)
          neigh = row(j)
          do k = rowptr(neigh),rowptr(neigh+1)-1
             icol_add = col(k)

             ! cycle if entry is in lower triangle
             if (icol_add.le.i) cycle

             icol = mask(icol_add)
             if (icol.eq.0) then
                nc = nc + 1
                wr01(icol_add) =  aij*val_AT(k)
                wn02(nc) = icol_add
                mask(icol_add) = nc     ! add mask
             else
                kcol = wn02(icol)
                wr01(kcol) = wr01(kcol) + aij*val_AT(k)
             end if
          end do
       end do

       if (present(scale)) then
          do j = 1,nc
             k = wn02(j)
             wr01(k) = scale(i) * wr01(k) * scale(k)
             mask(k) = 0
          end do
       else
          do j = 1,nc
             k = wn02(j)
             mask(k) = 0
          end do
       end if

   end subroutine ata_j

!****************************************************************************
! Code to compute the diagonal entries of C = S * A^T * W^2* A * S, 
! without forming C.
! Here S is an optional diagonal scaling matrix, supplied by user

! A is m x n, held in CSC format.
! There is no checking of the user's data.

   subroutine ata_diag(m,n,ptr,row,val,d,weight,scale)

   integer, intent(in) :: m  ! rows in A
   integer, intent(in) :: n  ! cols in A
   integer, intent(in) ::  ptr(n+1) ! column pointers for A = row pointers A^T 
   integer, intent(in) ::  row(:) ! row indices for A = col indices A^T 
   real(wp), intent(in) ::  val(:) ! entries of A
   real(wp), intent(out) :: d(n) ! on exit, holds (scaled) diagonal entries 
   real(wp), optional, intent(in) :: weight(m) ! weights
   real(wp), optional, intent(in) :: scale(n) ! scaling factors 

   integer :: i,j,q
   real(wp) :: diag

    if (present(weight)) then
      do q = 1,n
        diag = zero
        do j = ptr(q),ptr(q+1)-1
          i = row(j)
          diag = diag + val(j)*val(j)*weight(i)*weight(i)
        end do
        if (present(scale)) diag = scale(q)* diag *scale(q)
        d(q) = diag
      end do
    else
      do q = 1,n
        diag = zero
        do j = ptr(q),ptr(q+1)-1
          i = row(j)
          diag = diag + val(j)*val(j)
        end do
        if (present(scale)) diag = scale(q)* diag *scale(q)
        d(q) = diag
      end do
    end if

   end subroutine ata_diag
!****************************************************************************
! Code to scale the LS matrix C= A^T * W^2 * A, without forming C
! (at least, we just form one col of C at a time).
! The scaling is l2 scaling (ie 2 norm of each col of C is computed).

! A is m x n, held in CSC format = A^T in CSR format
! W^2*A is  held in CSR  format = A^T*W^2 in CSC format
! There is no checking of the data

   subroutine scalel2(m,n,ptr,row,val,rowptr,col,val_AT,scale, &
       wr01,wn02,mask,limit,flag)

   integer, intent(in) :: m  ! rows in A
   integer, intent(in) :: n  ! cols in A
   integer, intent(in) :: ptr(n+1) ! column pointers for A 
   integer, intent(in) :: row(:) ! row indices of A
   real(wp), intent(in) :: val(:) ! entries of A in CSC format
   integer, intent(in) :: rowptr(m+1) ! column pointers for A 
   integer, intent(in) :: col(:) ! row indices of A
   real(wp), intent(in) :: val_AT (:) ! entries of W^2*A in CSR format
   real(wp), intent(out) :: scale(n) ! scaling factors
   real(wp), intent(out) :: wr01(n) ! work array
   integer, intent(out) :: wn02(n) ! work array
   integer, intent(out) :: mask(n) ! work array
   integer, intent(in) :: limit(2)
   integer, intent(out) :: flag ! error flag

   integer :: i,icol,icol_add,j,k,kcol,nc,neigh
   integer(long) :: nzc
   real(wp) :: aij
   real(wp) :: temp

   flag = 0
   mask(1:n) = 0

   ! loop over rows of C
   nzc = 0_long
   do i = 1,n
       nc = 0
       do j = ptr(i),ptr(i+1)-1
          aij = val(j)
          neigh = row(j)
          do k = rowptr(neigh),rowptr(neigh+1)-1
             icol_add = col(k)
             icol = mask(icol_add)
             if (icol.eq.0) then
                nc = nc + 1
                wr01(icol_add) =  aij*val_AT(k)
                wn02(nc) = icol_add
                mask(icol_add) = nc     ! add mask
             else
                kcol = wn02(icol)
                wr01(kcol) = wr01(kcol) + aij*val_AT(k)
             end if
          end do
       end do

       temp = zero
       do j = 1,nc
          k = wn02(j)
          temp = temp + wr01(k)* wr01(k)
          mask(k) = 0
       end do

       scale(i) = sqrt(temp)

       nzc = nzc + nc

       if (nc.gt.limit(1)) then
          flag = mi35_error_colC
          return
       else if (nzc.gt.(limit(2)-n)/2) then
          flag = mi35_error_nzc
          return
       end if

    end do

    ! note: we can't have scale(i) = zero if null cols have been removed!
    do i = 1,n
   !   if (scale(i).gt.zero) then
         scale(i) = one/sqrt(scale(i))
   !   else
   !      scale(i) = one
   !   end if
    end do

    end subroutine scalel2

!****************************************************************************
   subroutine matmulC(n,ia,ja,a,ib,jb,b,ic,jc,c,mask)

    ! compute the product C=A*B of two sparse matrices A=(ia,ja,a) and
    ! B = (ib,jb,b) in CSR format. Output C in CSR format.
    ! Here we compute only the upper triangle of C
    ! since we are interested in A = B^T so that C = B^T * B is symmetric.
    ! Of course, upper triangle C in CSR format = lower triangle in CSC format
    !!! note: as we only want the case where C is square, we have 
    ! taken row size of A = col size of B = n

       ! n: row size of A = column size of B
    integer, intent (in) :: n
       ! ia: row pointer of A
       ! ja: column indices of A
       ! ib: row pointer of B
       ! jb: column indices of B
    integer, intent(in) :: ia(n+1),ib(*),ja(*),jb(*)
       ! ic: row pointer of C
       ! jc: column indices of C
    integer, intent(out) :: ic(n+1),jc(*)
       ! a: entry values of A
       ! b: entry values of B
       ! c: entry values of C
    real(wp), intent(in) :: a(*),b(*)
    real(wp), intent(out) :: c(*)
    integer, intent(out) :: mask(n) ! work array

    integer :: i,j,k,kcol,icol,icol_add,neigh,nz
    real(wp) :: aij

    ! initialise the mask array which is an array that
    ! has non-zero value if the
    ! column index already exists, in which case
    ! the value it holds is the index of that column

    ! code computes one row of C at a time (== one col)
    mask(1:n) = 0
    nz = 0
    ic(1) = 1
    do i = 1,n
       do j = ia(i),ia(i+1)-1
          aij = a(j)
          neigh = ja(j)
          do k = ib(neigh),ib(neigh+1)-1
             icol_add = jb(k)

             ! cycle if entry is in lower triangle
             if (icol_add.lt.i) cycle

             icol = mask(icol_add)
             if (icol.eq.0) then
                nz = nz + 1
                jc(nz) = icol_add
                c(nz) =  aij*b(k)
                mask(icol_add) = nz     ! add mask
             else
                c(icol) = c(icol) + aij*b(k)
             end if
          end do
       end do
       ! done with row i, so set mask to zero again
       do j = ic(i),nz
          kcol = jc(j)
          mask(kcol) = 0
       end do
       ic(i+1) = nz + 1
    end do

  end subroutine matmulC 
!*******************************

      subroutine compute_scaling(n,ia,ja,aa,iscale,scalesj,st,info)

      integer, intent(in) :: n  ! order of matrix
      integer, intent(in) ::  ia(n+1) ! column pointers for lower triangle.
      integer, intent(in) ::  ja(:) ! row indices of lower triangular part.
        ! diagonal entry must be first entry in each col. 
      real(wp), intent(in) ::  aa(:) ! entries of lower triangular part 
        ! (diagonal entry is non zero)
      integer, intent(in) :: iscale ! controls choice of scaling
      real(wp), intent(out) :: scalesj(n)
      integer, intent(out) :: st
      type(mi35_info), intent(inout) :: info    ! See derived-type declaration

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

        ! note: we can't have scale(i) = zero if null cols have been removed!
        do i = 1,n
        !   if (scalesj(i).gt.zero) then
              scalesj(i) = one/sqrt(scalesj(i))
        !   else
        !      scalesj(i) = one
        !   end if
        end do

      else if (iscale.eq.2) then
         call mc77_scale(n, ia, ja, aa, scalesj, info%flag77, st)

         if (st.ne.0) return

         if (info%flag77 < 0) info%flag = mi35_ERROR_MC77

      else if (iscale.eq.3) then
         call mc64_scale(n, ia, ja, aa, scalesj, st, sing, info%flag64)

         if (st.ne.0) return

         if (sing) then
            info%flag = mi35_WARNING_MC64
            return
         end if

         if (info%flag64 < 0) info%flag = mi35_ERROR_MC64

         if (info%flag64 == 2) info%flag = mi35_WARNING_MC64

      else if (iscale.eq.4) then
         ! diagonal scaling
         do i = 1,n
            jstrt = ia(i)
            scalesj(i) = aa(jstrt)
            ! note: cannot have zero diagonal entry if matrix has been checked
         !   if (scalesj(i).gt.zero) then
              scalesj(i) = one/sqrt(scalesj(i))
         !   else
         !     scalesj(i) = one
         !   end if
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
end module hsl_mi35_double
     
