! COPYRIGHT (c) 1995, 2007 Science & Technology Facilities Council
! COPYRIGHT (c) 1995 Timothy A. Davis and Patrick Amestoy
!
! Version: 3.3.1
! For version history see ChangeLog
!
MODULE hsl_mc68_integer
      USE hsl_zb01_integer
      IMPLICIT NONE
      PRIVATE

! ---------------------------------------------------
! Precision
! ---------------------------------------------------
      INTEGER, PARAMETER :: myreal_mc68 = kind(1.0D0)
      INTEGER, PARAMETER :: myint = kind(1)
      INTEGER, PARAMETER :: long = selected_int_kind(18)

! ---------------------------------------------------
! Error flags
! ---------------------------------------------------
      INTEGER (myint), PARAMETER :: &
        mc68_err_memory_alloc = -1, & ! memory alloc error
        mc68_err_memory_dealloc = -2, & ! memory dealloc error
        mc68_err_n = -3, & ! n<1
        mc68_err_ord = -4, & ! ord not associated with an ordering
        mc68_err_metis = -5, & ! MeTiS ordering requested but not linked
        mc68_err_zb01 = -6 ! None (de)allocation error from call to zb01_expand1

! ---------------------------------------------------
! Warning flags
! ---------------------------------------------------
      INTEGER (myint), PARAMETER :: mc68_warn_diag = 1, & ! No diags and ord=4
        mc68_warn_rank = 2, & ! Matrix rank deficient when using ord=4
        mc68_warn_rank_diag = 3 ! Matrix rank deficient, no diags and ord=4

! ---------------------------------------------------
! Derived type definitions
! ---------------------------------------------------
      TYPE, PUBLIC :: mc68_control
        INTEGER :: lp = 6 ! stream number for error messages
        INTEGER :: wp = 6 ! stream number for warning messages
        INTEGER :: mp = 6 ! stream number for diagnostic messages
        INTEGER :: nemin = 1 ! stream number for diagnostic messages
        INTEGER :: print_level = 0 ! amount of informational output required
        INTEGER :: row_full_thresh = 100 ! percentage threshold for full row
        INTEGER :: row_search = 10 ! Num of rows searched for pivot with ord=6
      END TYPE mc68_control

      TYPE, PUBLIC :: mc68_info
        INTEGER :: flag = 0 ! error/warning flag
        INTEGER :: iostat = 0 ! holds Fortran iostat parameter
        INTEGER :: stat = 0 ! holds Fortran stat parameter
        INTEGER :: out_range = 0 ! holds number of out of range entries ignored
        INTEGER :: duplicate = 0 ! holds number of duplicate entries
        INTEGER :: n_compressions = 0 ! holds number of compressions in order
        INTEGER :: n_zero_eigs = -1 ! holds the number of zero eigs from ma47
        INTEGER(long) :: l_workspace = 0 ! holds length of workspace iw used in
          ! order
        INTEGER :: zb01_info = 0 ! holds flag from zb01_expand1 call
        INTEGER :: n_dense_rows = 0 ! holds number of dense rows from amdd
      END TYPE mc68_info

      INTERFACE mc68_order
        MODULE PROCEDURE mc68_order_integer
      END INTERFACE

      PUBLIC mc68_order

CONTAINS

! ---------------------------------------------------------------

! subroutine mc68_order_integer(ord,A,perm,control,info) constructs
! an elimination order PERM for a symmetric matrix A using a chosen
! ordering ORD
SUBROUTINE mc68_order_integer(ord,n,ptr,row,perm,control,info, min_l_workspace)
   INTEGER, INTENT (IN) :: ord ! Specifies which ordering is to be used.
      ! The choice is as follows:
      ! ord = 1 : Approximate minimum degree with provision for dense rows
      ! 2 : Minimum degree
      ! 3 : MeTiS
      ! 4 : MA47 ordering for indefinite matrices
   INTEGER, INTENT (IN) :: n ! Must hold the number of rows in A
   INTEGER, INTENT (IN) :: ptr(n+1) ! ptr(j) holds position in row of start of
      ! row indices for column j. ptr(n)+1 must equal the number of entries
      ! stored + 1. Only the lower triangular entries are stored with no
      ! duplicates or out-of-range entries
   INTEGER, INTENT (IN) :: row(:) ! size at least ptr(n+1)-1
   INTEGER (myint), INTENT (OUT) :: perm(n) ! Holds elimination order on output
   TYPE (mc68_control), INTENT (IN) :: control
   TYPE (mc68_info), INTENT (OUT) :: info ! information type:
      ! info%flag values:
      ! = 0 if successful
      ! = MC68_ERR_MEMORY_ALLOC if memory allocation failed
      ! = MC68_ERR_MEMORY_DEALLOC if memory deallocation failed
      ! = MC68_ERR_N if n<1
      ! = MC68_ERR_METIS if MeTiS ordering is requested but not linked
   INTEGER, INTENT (IN), OPTIONAL :: min_l_workspace ! If present, specifies
      ! the minimum amount of workspace to be use

   ! ---------------------------------------------
   ! Local variables
   ! ---------------------------------------------
   INTEGER, ALLOCATABLE :: ipe(:) ! copy of pointers which is later modified
   INTEGER, ALLOCATABLE :: iw(:) ! copy of row indices
   INTEGER, ALLOCATABLE :: work1(:) ! work array
   INTEGER, ALLOCATABLE :: work2(:) ! work array
   INTEGER, ALLOCATABLE :: work3(:) ! work array
   INTEGER, ALLOCATABLE :: work4(:) ! work array
   INTEGER, ALLOCATABLE :: work6(:) ! work array
   INTEGER, ALLOCATABLE :: work7(:) ! work array
   INTEGER, ALLOCATABLE :: work8(:) ! work array
   INTEGER, ALLOCATABLE :: work9(:) ! work array
   INTEGER :: lp ! stream number for error messages
   INTEGER :: wp ! stream number for warning messages
   INTEGER :: mp ! stream number for diagnostic messages
   INTEGER (long) :: iwlen ! length of iw
   INTEGER :: i, k1, iwfr, j, diag
   REAL (myreal_mc68) :: thresh
   LOGICAL :: printe ! errors to be printed?
   LOGICAL :: printw ! warnings to be printed?
   LOGICAL :: printi ! basic diagnostic to be printed?
   LOGICAL :: printd ! additional diagnostic to be printed?

   ! ---------------------------------------------
   ! Set stream numbers
   ! ---------------------------------------------
   lp = control%lp
   wp = control%wp
   mp = control%mp

   ! ---------------------------------------------
   ! Printing levels
   ! ---------------------------------------------
   printe = (control%print_level>=0 .AND. lp>=0)
   printw = (control%print_level>=0 .AND. wp>=0)
   printi = (control%print_level>=1 .AND. mp>=0)
   printd = (control%print_level>=2 .AND. mp>=0)
   IF (printi) THEN
      WRITE (mp,'(a)') ' '
      WRITE (mp,'(a)') 'MC68_order:'
   END IF

   ! Initialise info
   info%flag = 0
   info%stat = 0

   ! ---------------------------------------------
   ! Check that restrictions are adhered to
   ! ---------------------------------------------
   IF (n<1) THEN
      info%flag = mc68_err_n
      IF (printe) &
         CALL mc68_print_message(info%flag,lp,context='mc68_order')
      RETURN
   END IF

   SELECT CASE (ord)
   CASE (1)
      IF (printi) WRITE (mp,'(a60)') 'Approximate minimum degree ordering'
   CASE (2)
      IF (printi) &
         WRITE (mp,'(a60)') 'Minimum degree ordering using methodology of MA27'
   CASE (3)
      IF (printi) WRITE (mp,'(a60)') 'Nested bisection ordering using MeTiS'
   CASE (4)
      IF (printi) &
         WRITE (mp,'(a60)') &
            'Ordering for indefinite matrices using methodology of MA47'
   CASE DEFAULT
      info%flag = mc68_err_ord
      IF (printe) &
         CALL mc68_print_message(info%flag,lp,context='mc68_order')
      RETURN
   END SELECT

   IF (printi) WRITE (mp,'(a,i15)') 'n  =  ', n

   IF (n==1) THEN
      ! ---------------------------------------------
      ! Matrix is of order 1 so no need to call ordering subroutines
      ! ---------------------------------------------
      perm(1) = 1
      IF (printi) THEN
         WRITE (mp,'(a)') ' '
         WRITE (mp,'(a)') 'Matrix of order 1'
      END IF
      GO TO 20
   END IF

   IF (ord==3) THEN
      IF (.not.check_for_metis()) THEN
         info%flag = mc68_err_metis
         IF (printe) &
            CALL mc68_print_message(info%flag,lp,context='mc68_order')
         RETURN
      END IF
   END IF

   SELECT CASE (ord)
   CASE (1)
      ! ----------------------------------
      ! Approximate minimum degree ordering with provision for
      ! dense rows
      ! ----------------------------------
      iwlen = 2*ptr(n+1) + n
      IF (present(min_l_workspace)) iwlen = max(iwlen,min_l_workspace+0_long)
      info%l_workspace = iwlen

      ALLOCATE (work1(n),work8(10),ipe(n+1),iw(iwlen),STAT=info%stat)
      IF (info%stat/=0) THEN
         info%flag = mc68_err_memory_alloc
         IF (printe) &
            CALL mc68_print_message(info%flag,lp,context='mc68_order')
         RETURN
      END IF

      call half_to_full(n, ptr, row, ipe, iw)
      iwfr = ipe(n+1) ! Index of next entry to be filled in iw

      ! Form ordering
      work1(1:n) = ipe(2:n+1) - ipe(1:n)
      work8(1) = 6
      work8(2) = 6
      work8(3) = -1
      work8(4) = 1 ! enable dense row detection
      work8(5) = huge(0)
      work8(6:10) = 0
      CALL amdd(n,iwlen,ipe,iwfr,work1,iw,perm,work8,info)

      ! Deallocate arrays
      DEALLOCATE (ipe,iw,work1,work8,STAT=info%stat)
      IF (info%stat/=0) THEN
         info%flag = mc68_err_memory_dealloc
         IF (printe) &
            CALL mc68_print_message(info%flag,lp,context='mc68_order')
         RETURN
       END IF

   CASE (2)
      ! ----------------------------------
      ! Minimum degree ordering using methodology of MA27
      ! ----------------------------------

      ! Set length of iw
      iwlen = 2*ptr(n+1) - 1 + n
      IF (present(min_l_workspace)) iwlen = max(iwlen,min_l_workspace+0_long)

      ! Allocate required arrays
      ALLOCATE (work2(n),ipe(n+1),iw(iwlen),STAT=info%stat)
      IF (info%stat/=0) THEN
         info%flag = mc68_err_memory_alloc
         IF (printe) &
            CALL mc68_print_message(info%flag,lp,context='mc68_order')
         RETURN
      END IF

      ! IPE and IW hold matrix in VARIANT of CSC format
      ! iw(ipe(i)) must be length of column i.
      call half_to_full_variant(n, ptr, row, ipe, iw, .true.)
      ! Additionally:
      ! ipe(i) must be 0 if col i has no off-diagonal non-zeros.
      DO j = 1, n
         IF (iw(ipe(j))==0) ipe(j) = 0
      END DO

      thresh = float(control%row_full_thresh)/100.0
      iwfr = ipe(n+1)

      CALL mc68_min_deg_anal(n,ipe,iw,iwlen,iwfr,work2,huge(0), &
         info%n_compressions,thresh,info)

      ! set ipe correctly
      DO i = 1, n
         IF (work2(i).ne.0) CYCLE
         k1 = i
         DO WHILE(work2(k1).eq.0)
            k1 = -ipe(k1)
         END DO
         ipe(i) = -k1
      END DO
      info%n_zero_eigs = -1

      CALL mc68_min_deg_tree_search(n,ipe,work2,perm,control%nemin,info)

      IF (info%flag<0 .AND. printe) THEN
         CALL mc68_print_message(info%flag,lp,context='mc68_order')
         RETURN
      END IF

      ! Deallocate required arrays
      DEALLOCATE (work2,ipe,iw,STAT=info%stat)
      IF (info%stat/=0) THEN
         info%flag = mc68_err_memory_dealloc
         IF (printe) &
            CALL mc68_print_message(info%flag,lp,context='mc68_order')
         RETURN
      END IF

   CASE (3)
      ! ----------------------------------
      ! MeTiS ordering
      ! ----------------------------------

      ! Set length of iw
      iwlen = 2*ptr(n+1) - 2
      IF (present(min_l_workspace)) iwlen = max(iwlen,min_l_workspace+0_long)
      info%l_workspace = iwlen

      ! Allocate arrays
      ALLOCATE (work1(8),ipe(n+1),iw(iwlen),work3(n),STAT=info%stat)
      IF (info%stat/=0) THEN
         info%flag = mc68_err_memory_alloc
         IF (printe) &
            CALL mc68_print_message(info%flag,lp,context='mc68_order')
         RETURN
      END IF

      ! Fill iw and ipe removing any diagonal entries
      call half_to_full(n, ptr, row, ipe, iw)

      ! Carry out ordering
      work1(1) = 0 ! MeTiS defaults
      CALL metis_nodend(n,ipe,iw,1,work1,work3,perm)

      ! Compression information not returned from meTiS
      info%n_compressions = 0
      info%n_zero_eigs = -1
      info%n_dense_rows = -1

      ! Deallocate arrays
      DEALLOCATE (ipe,iw,work1,work3,STAT=info%stat)
      IF (info%stat/=0) THEN
         info%flag = mc68_err_memory_dealloc
         IF (printe) &
            CALL mc68_print_message(info%flag,lp,context='mc68_order')
         RETURN
      END IF

   CASE (4)
      ! ----------------------------------
      ! Ordering for indefinite matrices using methodology of MA47
      ! ----------------------------------
      ! Set length of iw
      iwlen = 4*( ptr(n+1) - 1 ) + n
      IF (present(min_l_workspace)) iwlen = max(iwlen,min_l_workspace+0_long)

      ! Allocate required arrays
      ALLOCATE (work2(n),ipe(n+1),iw(iwlen),work3(n),work4(n),work6(n), &
         work7(n),work8(n),work9(n),STAT=info%stat)
      IF (info%stat/=0) THEN
         info%flag = mc68_err_memory_alloc
         IF (printe) CALL mc68_print_message(info%flag,lp, &
           context='mc68_order')
         RETURN
      END IF

      ! IPE and IW hold matrix in VARIANT of CSC format
      ! iw(ipe(i)) must be length of column i.
      call half_to_full_variant(n, ptr, row, ipe, iw, .false., ndiag=diag, &
         dflag=work2)

      iwfr = ipe(n+1) ! Index of next entry to be filled in iw
      ! Warn if no diagonal entries present
      IF (diag==0) info%flag = mc68_warn_diag

      ! Carryout ordering
      CALL mc68_ma47_analyse(n,ipe,iw,iwlen,iwfr,control%row_search,work3, &
         work2,work4,info)
      info%l_workspace = iwlen
      info%n_dense_rows = -1
      IF (info%flag<0) THEN
         IF (printe) CALL mc68_print_message(info%flag,lp, &
           context='mc68_order')
         RETURN
      END IF

      CALL mc68_ma47_treesearch(n,work4,work3,work2,perm,control%nemin,info)

      IF (info%flag>0 .and. printw) &
         CALL mc68_print_message(info%flag,lp,context='mc68_order')

      IF (info%flag<0) THEN
         IF (printe) &
            CALL mc68_print_message(info%flag,lp,context='mc68_order')
         RETURN
      END IF

      ! Deallocate required arrays
      DEALLOCATE (work3,work4,work6,work7,work8,work9,work2,ipe,iw, &
         STAT=info%stat)
      IF (info%stat<0) THEN
         info%flag = mc68_err_memory_dealloc
         IF (printe) &
            CALL mc68_print_message(info%flag,lp,context='mc68_order')
         RETURN
      END IF
   END SELECT

   20 continue
   IF (printd) THEN
      ! ---------------------------------------------
      ! Print out perm
      ! ---------------------------------------------
      WRITE (mp,'(a7)') 'perm = '
      WRITE (mp,'(5i15)') (perm(i),i=1,n)
   ELSE IF (printi) THEN
      ! ---------------------------------------------
      ! Print out first few entries of perm
      ! ---------------------------------------------
      WRITE (mp,'(a21)') 'perm(1:min(5,n)) = '
      WRITE (mp,'(5i15)') (perm(i),i=1,min(5,n))
   END IF

   IF (printi) THEN
      CALL mc68_print_message(info%flag,mp,context='mc68_order')
   ELSE IF (printw .AND. (info%stat>0)) THEN
      CALL mc68_print_message(info%flag,wp,context='mc68_order')
   END IF
END SUBROUTINE mc68_order_integer

! Convert a matrix in half storage to one in full storage using
! VARIANT CSC format.
! Optionally drop any diagonal entries if drop=.true.
subroutine half_to_full_variant(n, ptr, row, ptr2, row2, drop, ndiag, dflag)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   integer, dimension(*), intent(out) :: ptr2
   integer, dimension(*), intent(out) :: row2
   logical, intent(in) :: drop
   integer, optional, intent(out) :: ndiag ! number of diagonal entries
      ! (only if drop=true)
   integer, optional, dimension(n), intent(out) :: dflag ! flag array for
      ! diagonal entries (only if drop=true)
      ! dflag(i) = 1 if present in col i, -1 if not

   integer :: i, j, k

   if(present(dflag)) dflag(:) = -1
   if(present(ndiag)) ndiag = 0

   ! Set ptr2(j) to hold no. nonzeros in column j
   ptr2(1:n+1) = 0
   DO j = 1, n
      DO k = ptr(j), ptr(j+1) - 1
         i = row(k)
         IF (j/=i) THEN
            ptr2(i) = ptr2(i) + 1
            ptr2(j) = ptr2(j) + 1
         ELSEIF(.not.drop) THEN
            ptr2(i) = ptr2(i) + 1
            if(present(ndiag)) ndiag = ndiag + 1
            if(present(dflag)) dflag(i) = 1
         END IF
      END DO
   END DO

   ! Set ptr2(j) to point to where row indices will end in row2
   row2(1) = ptr2(1)
   ptr2(1) = ptr2(1) + 1
   DO j = 2, n
      row2(ptr2(j-1)+1) = ptr2(j)
      ptr2(j) = ptr2(j-1) + ptr2(j) + 1
   END DO
   ptr2(n+1) = ptr2(n) + 1

   ! Fill ptr2 and row2
   DO j = 1, n
      DO k = ptr(j), ptr(j+1) - 1
         i = row(k)
         IF (j/=i) THEN
            row2(ptr2(i)) = j
            row2(ptr2(j)) = i
            ptr2(i) = ptr2(i) - 1
            ptr2(j) = ptr2(j) - 1
         ELSEIF(.not.drop) THEN
            row2(ptr2(i)) = j
            ptr2(i) = ptr2(i) - 1
         END IF
      END DO
   END DO

end subroutine half_to_full_variant

! Convert a matrix in half storage to one in full storage.
! Drops any diagonal entries.
subroutine half_to_full(n, ptr, row, ptr2, row2)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   integer, dimension(*), intent(out) :: ptr2
   integer, dimension(*), intent(out) :: row2

   integer :: i, j, k

   ! Set ptr2(j) to hold no. nonzeros in column j
   ptr2(1:n+1) = 0
   DO j = 1, n
      DO k = ptr(j), ptr(j+1) - 1
         i = row(k)
         IF (j/=i) THEN
            ptr2(i) = ptr2(i) + 1
            ptr2(j) = ptr2(j) + 1
         END IF
      END DO
   END DO

   ! Set ptr2(j) to point to where row indices will end in row2
   DO j = 2, n
      ptr2(j) = ptr2(j-1) + ptr2(j)
   END DO
   ptr2(n+1) = ptr2(n) + 1

   ! Fill ptr2 and row2
   DO j = 1, n
      DO k = ptr(j), ptr(j+1) - 1
         i = row(k)
         IF (j/=i) THEN
            row2(ptr2(i)) = j
            row2(ptr2(j)) = i
            ptr2(i) = ptr2(i) - 1
            ptr2(j) = ptr2(j) - 1
         END IF
      END DO
   END DO
   DO j = 1, n
      ptr2(j) = ptr2(j) + 1
   END DO

end subroutine half_to_full

! This subroutine returns true if a working MeTiS routine is linked, or
! false if the dummy metis.f is used.
logical function check_for_metis()
   integer :: ptr(3), row(2), icntl(8), invp(2), perm(2)

   ptr = (/ 1, 2, 3 /)
   row = (/ 2, 1 /)
   icntl(1) = 0 ! MeTiS defaults

   CALL metis_nodend(2,ptr,row,1,icntl,invp,perm)

   check_for_metis = (perm(1).ne.-1)
end function check_for_metis

SUBROUTINE mc68_print_message(flag,unit,context)
! Prints out errors and warnings according to value of flag

! flag: is an integer scaler of intent(in). It is the information flag
! whose corresponding error message is printed
        INTEGER (myint), INTENT (IN) :: flag

! unit: is an integer scaler of intent(in). It is the unit number the
! error message should be printed on
        INTEGER (myint), INTENT (IN) :: unit

! context: is an optional assumed size character array of intent(in).
! It describes the context under which the error occured
        CHARACTER (len=*), OPTIONAL, INTENT (IN) :: context

        INTEGER (myint) :: length

        IF (unit<=0) RETURN

        IF (flag>0) THEN
          WRITE (unit,advance='yes',fmt='('' WARNING: '')')
        ELSE IF (flag<0) THEN
          WRITE (unit,advance='yes',fmt='('' ERROR: '')')
        END IF

        IF (present(context)) THEN
          length = len_trim(context)
          WRITE (unit,advance='no',fmt='('' '', a,'': '')') context(1:length)
        END IF

        SELECT CASE (flag)
        CASE (0)
          WRITE (unit,'(A)') 'successful completion'

        CASE (mc68_err_memory_alloc)
          WRITE (unit,'(A)') 'memory allocation failure'

        CASE (mc68_err_memory_dealloc)
          WRITE (unit,'(A)') 'memory deallocation failure'

        CASE (mc68_err_n)
          WRITE (unit,'(A)') 'restriction n>=1 violated'

        CASE (mc68_err_ord)
          WRITE (unit,'(A)') 'ord is not associated with an ordering'

        CASE (mc68_err_metis)
          WRITE (unit,'(A)') 'MeTiS ordering requested but not linked'

        CASE (mc68_err_zb01)
          WRITE (unit,'(A)') 'temporary file failure'

        CASE (mc68_warn_diag)
          WRITE (unit,'(A)') 'no diagonal entries'

        CASE (mc68_warn_rank)
          WRITE (unit,'(A)') 'matrix rank deficient'

        CASE (mc68_warn_rank+mc68_warn_diag)
          WRITE (unit,'(A)') 'no diagonal entries and matrix rank deficient'

        END SELECT

END SUBROUTINE mc68_print_message




SUBROUTINE mc68_min_deg_anal(n,ipe,iw,lw,iwfr,nv,iovflo,ncmpa,fratio,info)
! F90 version of subroutine MA27HD.

! ANALYSIS SUBROUTINE

! GIVEN REPRESENTATION OF THE WHOLE MATRIX (EXCLUDING DIAGONAL)
! PERFORM MINIMUM DEGREE ORDERING, CONSTRUCTING TREE POINTERS.
! IT WORKS WITH SUPERVARIABLES WHICH ARE COLLECTIONS OF ONE OR MORE
! VARIABLES, STARTING WITH SUPERVARIABLE I CONTAINING VARIABLE I FOR
! I = 1,2,...,N. ALL VARIABLES IN A SUPERVARIABLE ARE ELIMINATED
! TOGETHER. EACH SUPERVARIABLE HAS AS NUMERICAL NAME THAT OF ONE
! OF ITS VARIABLES (ITS PRINCIPAL VARIABLE).

! N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
! IPE(I) MUST BE SET TO POINT TO THE POSITION IN IW OF THE
! START OF ROW I OR HAVE THE VALUE ZERO IF ROW I HAS NO OFF-
! DIAGONAL NON-ZEROS. DURING EXECUTION IT IS USED AS FOLLOWS. IF
! SUPERVARIABLE I IS ABSORBED INTO SUPERVARIABLE J THEN IPE(I)=-J.
! IF SUPERVARIABLE I IS ELIMINATED THEN IPE(I) EITHER POINTS TO THE
! LIST OF SUPERVARIABLES FOR CREATED ELEMENT I OR IS ZERO IF
! THE CREATED ELEMENT IS NULL. IF ELEMENT I
! IS ABSORBED INTO ELEMENT J THEN IPE(I)=-J.
! IW MUST BE SET ON ENTRY TO HOLD LISTS OF VARIABLES BY
! ROWS, EACH LIST BEING HEADED BY ITS LENGTH.
! DURING EXECUTION THESE LISTS ARE REVISED AND HOLD
! LISTS OF ELEMENTS AND SUPERVARIABLES. THE ELEMENTS
! ALWAYS HEAD THE LISTS. WHEN A SUPERVARIABLE
! IS ELIMINATED ITS LIST IS REPLACED BY A LIST OF SUPERVARIABLES
! IN THE NEW ELEMENT.
! LW MUST BE SET TO THE LENGTH OF IW. IT IS NOT ALTERED.
! IWFR MUST BE SET TO THE POSITION IN IW OF THE FIRST FREE VARIABLE.
! IT IS REVISED DURING EXECUTION AND CONTINUES TO HAVE THIS MEANING.
! NV(JS) NEED NOT BE SET. DURING EXECUTION IT IS ZERO IF
! JS IS NOT A PRINCIPAL VARIABLE AND IF IT IS IT HOLDS
! THE NUMBER OF VARIABLES IN SUPERVARIABLE JS. FOR ELIMINATED
! VARIABLES IT IS SET TO THE DEGREE AT THE TIME OF ELIMINATION.
! IOVFLO should be set to a high legitimate integer.
! It is used as a flag.
! NCMPA number of compresses.
! FRATIO is the density of rows regarded as dense.

        REAL (myreal_mc68), INTENT (IN) :: fratio
        INTEGER (long), INTENT (IN) :: lw
        INTEGER, INTENT (IN) :: n, iovflo
        INTEGER, INTENT (INOUT) :: iwfr, ipe(n), iw(lw)
        INTEGER, INTENT (OUT) :: nv(n), ncmpa
        TYPE (mc68_info), INTENT (INOUT) :: info

! Local variables

! NXT(JS) NEED NOT BE SET. DURING EXECUTION IT IS THE NEXT
! SUPERVARIABLE HAVING THE SAME DEGREE AS JS, OR ZERO
! IF IT IS LAST IN ITS LIST.
! LST(JS) NEED NOT BE SET. DURING EXECUTION IT IS THE
! LAST SUPERVARIABLE HAVING THE SAME DEGREE AS JS OR
! -(ITS DEGREE) IF IT IS FIRST IN ITS LIST.
! IPD(ID) NEED NOT BE SET. DURING EXECUTION IT
! IS THE FIRST SUPERVARIABLE WITH DEGREE ID OR ZERO
! IF THERE ARE NONE.
! FLAG IS USED AS WORKSPACE FOR ELEMENT AND SUPERVARIABLE FLAGS.
! WHILE THE CODE IS FINDING THE DEGREE OF SUPERVARIABLE IS
! FLAG HAS THE FOLLOWING VALUES.
! A) FOR THE CURRENT PIVOT/NEW ELEMENT ME
! FLAG(ME)=-1
! B) FOR VARIABLES JS
! FLAG(JS)=-1 IF JS IS NOT A PRINCIPAL VARIABLE
! FLAG(JS)=0 IF JS IS A SUPERVARIABLE IN THE NEW ELEMENT
! FLAG(JS)=NFLG IF JS IS A SUPERVARIABLE NOT IN THE NEW
! ELEMENT THAT HAS BEEN COUNTED IN THE DEGREE
! CALCULATION
! FLAG(JS).GT.NFLG IF JS IS A SUPERVARIABLE NOT IN THE NEW
! ELEMENT THAT HAS NOT BEEN COUNTED IN THE DEGREE
! CALCULATION
! C) FOR ELEMENTS IE
! FLAG(IE)=-1 IF ELEMENT IE HAS BEEN MERGED INTO ANOTHER
! FLAG(IE)=-NFLG IF ELEMENT IE HAS BEEN USED IN THE DEGREE
! CALCULATION FOR IS.
! FLAG(IE).LT.-NFLG IF ELEMENT IE HAS NOT BEEN USED IN THE
! DEGREE CALCULATION FOR IS
! ..
! .. Array Arguments ..
        INTEGER, ALLOCATABLE, DIMENSION (:) :: flag, ipd, lst, nxt

! ..
! .. Local Scalars ..
! LIMIT  Limit on number of variables for putting node in root.
! NVROOT Number of variables in the root node
! ROOT   Index of the root node (N+1 if none chosen yet).
        INTEGER :: i, id, idl, idn, ie, ip, is, jp, jp1, jp2, js, k, k1, k2, &
          ke, kp, kp0, kp1, kp2, ks, l, len, limit, ln, ls, lwfr, md, me, ml, &
          ms, nel, nflg, np, np0, ns, nvpiv, nvroot, root
! ..
! ..
! .. Intrinsic Functions ..
        INTRINSIC abs, min, nint
! ..
! If a column of the reduced matrix has relative density greater than
! CNTL(2), it is forced into the root. All such columns are taken to
! have sparsity pattern equal to their merged patterns, so the fill
! and operation counts may be overestimated.

! IS,JS,KS,LS,MS,NS ARE USED TO REFER TO SUPERVARIABLES.
! IE,JE,KE ARE USED TO REFER TO ELEMENTS.
! IP,JP,KP,K,NP ARE USED TO POINT TO LISTS OF ELEMENTS.
! OR SUPERVARIABLES.
! ID IS USED FOR THE DEGREE OF A SUPERVARIABLE.
! MD IS USED FOR THE CURRENT MINIMUM DEGREE.
! IDN IS USED FOR THE NO. OF VARIABLES IN A NEWLY CREATED ELEMENT
! NEL IS USED TO HOLD THE NO. OF VARIABLES THAT HAVE BEEN
! ELIMINATED.
! ME=MS IS THE NAME OF THE SUPERVARIABLE ELIMINATED AND
! OF THE ELEMENT CREATED IN THE MAIN LOOP.
! NFLG IS USED FOR THE CURRENT FLAG VALUE IN ARRAY FLAG. IT STARTS
! WITH THE VALUE IOVFLO AND IS REDUCED BY 1 EACH TIME IT IS USED
! UNTIL IT HAS THE VALUE 2 WHEN IT IS RESET TO THE VALUE IOVFLO.

! INITIALIZATIONS

        ALLOCATE (flag(n),ipd(n),lst(n),nxt(n),STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_alloc
          RETURN
        END IF
        info%n_dense_rows = 0
        ipd(:) = 0
        nv(:) = 1
        flag(:) = iovflo
        md = 1
        ncmpa = 0
        nflg = iovflo
        nel = 0
        root = n + 1
        nvroot = 0

! LINK TOGETHER VARIABLES HAVING SAME DEGREE
        DO is = 1, n
          k = ipe(is)
          IF (k<=0) THEN
! WE HAVE A VARIABLE THAT CAN BE ELIMINATED AT ONCE BECAUSE THERE
! IS
! NO OFF-DIAGONAL NON-ZERO IN ITS ROW.
            nel = nel + 1
            flag(is) = -1
            nxt(is) = 0
            lst(is) = 0
          ELSE
            id = iw(k) + 1
            ns = ipd(id)
            IF (ns>0) lst(ns) = is
            nxt(is) = ns
            ipd(id) = is
            lst(is) = -id
          END IF
        END DO

! START OF MAIN LOOP
        DO ml = 1, n
! LEAVE LOOP IF ALL VARIABLES HAVE BEEN ELIMINATED.
          IF (nel+nvroot+1>=n) THEN
            GO TO 140
          ELSE
! FIND NEXT SUPERVARIABLE FOR ELIMINATION.
            DO id = md, n
              ms = ipd(id)
              IF (ms>0) GO TO 10
            END DO
10          md = id
! NVPIV HOLDS THE NUMBER OF VARIABLES IN THE PIVOT.
            nvpiv = nv(ms)

! REMOVE CHOSEN VARIABLE FROM LINKED LIST
            ns = nxt(ms)
            nxt(ms) = 0
            lst(ms) = 0
            IF (ns>0) lst(ns) = -id
            ipd(id) = ns
            me = ms
            nel = nel + nvpiv
! IDN HOLDS THE DEGREE OF THE NEW ELEMENT.
            idn = 0

! RUN THROUGH THE LIST OF THE PIVOTAL SUPERVARIABLE, SETTING TREE
! POINTERS AND CONSTRUCTING NEW LIST OF SUPERVARIABLES.
! KP IS A POINTER TO THE CURRENT POSITION IN THE OLD LIST.
            kp = ipe(me)
            flag(ms) = -1
! IP POINTS TO THE START OF THE NEW LIST.
            ip = iwfr
! LEN HOLDS THE LENGTH OF THE LIST ASSOCIATED WITH THE PIVOT.
            len = iw(kp)
            DO kp1 = 1, len
              kp = kp + 1
              ke = iw(kp)
! JUMP IF KE IS AN ELEMENT THAT HAS NOT BEEN MERGED INTO
! ANOTHER.
              IF (flag(ke)<=-2) THEN
! SEARCH VARIABLE LIST OF ELEMENT KE, USING JP AS A POINTER TO
! IT.
                ie = ke
                jp = ipe(ie)
                ln = iw(jp)
              ELSE
! JUMP IF KE IS AN ELEMENT THAT HAS BEEN MERGED INTO ANOTHER
! OR
! IS
! A SUPERVARIABLE THAT HAS BEEN ELIMINATED.
                IF (flag(ke)<=0) THEN
                  IF (ipe(ke)/=-root) THEN
                    GO TO 30
                  ELSE
! KE has been merged into the root
                    ke = root
                    IF (flag(ke)<=0) GO TO 30
                  END IF
                END IF
! WE HAVE A SUPERVARIABLE. PREPARE TO SEARCH REST OF LIST.
                jp = kp - 1
                ln = len - kp1 + 1
                ie = ms
              END IF

! SEARCH FOR DIFFERENT SUPERVARIABLES AND ADD THEM TO THE NEW
! LIST,
! COMPRESSING WHEN NECESSARY. THIS LOOP IS EXECUTED ONCE FOR
! EACH ELEMENT IN THE LIST AND ONCE FOR ALL THE SUPERVARIABLES
! IN THE LIST.
              DO jp1 = 1, ln
                jp = jp + 1
                is = iw(jp)
! JUMP IF IS IS NOT A PRINCIPAL VARIABLE OR HAS ALREADY BEEN
! COUNTED.
                IF (flag(is)<=0) THEN
                  IF (ipe(is)==-root) THEN
! IS has been merged into the root
                    is = root
                    iw(jp) = root
                    IF (flag(is)<=0) GO TO 20
                  ELSE
                    GO TO 20
                  END IF
                END IF
                flag(is) = 0
                IF (iwfr>=lw) THEN
! PREPARE FOR COMPRESSING IW BY ADJUSTING POINTERS AND
! LENGTHS SO THAT THE LISTS BEING SEARCHED IN THE INNER AND
! OUTER
! LOOPS CONTAIN ONLY THE REMAINING ENTRIES.
                  ipe(ms) = kp
                  iw(kp) = len - kp1
                  ipe(ie) = jp
                  iw(jp) = ln - jp1
! COMPRESS IW
                  CALL mc68_compress(n,ipe,iw,ip-1,lwfr,ncmpa)
! COPY NEW LIST FORWARD
                  jp2 = iwfr - 1
                  iwfr = lwfr
                  IF (ip<=jp2) THEN
                    DO jp = ip, jp2
                      iw(iwfr) = iw(jp)
                      iwfr = iwfr + 1
                    END DO
                  END IF
! ADJUST POINTERS FOR THE NEW LIST AND THE LISTS BEING
! SEARCHED.
                  ip = lwfr
                  jp = ipe(ie)
                  kp = ipe(me)
                END IF
! STORE IS IN NEW LIST.
                iw(iwfr) = is
                idn = idn + nv(is)
                iwfr = iwfr + 1
! REMOVE IS FROM DEGREE LINKED LIST
                ls = lst(is)
                lst(is) = 0
                ns = nxt(is)
                nxt(is) = 0
                IF (ns>0) lst(ns) = ls
                IF (ls<0) THEN
                  ls = -ls
                  ipd(ls) = ns
                ELSE IF (ls>0) THEN
                  nxt(ls) = ns
                END IF
20              CONTINUE
              END DO
! JUMP IF WE HAVE JUST BEEN SEARCHING THE VARIABLES AT THE END
! OF
! THE LIST OF THE PIVOT.
              IF (ie==ms) THEN
                GO TO 40
              ELSE
! SET TREE POINTER AND FLAG TO INDICATE ELEMENT IE IS ABSORBED
! INTO
! NEW ELEMENT ME.
                ipe(ie) = -me
                flag(ie) = -1
              END IF
30            CONTINUE
            END DO

! STORE THE DEGREE OF THE PIVOT.
40          nv(ms) = idn + nvpiv
! JUMP IF NEW ELEMENT IS NULL.
            IF (iwfr==ip) THEN
              ipe(me) = 0
            ELSE
              k1 = ip
              k2 = iwfr - 1

! RUN THROUGH NEW LIST OF SUPERVARIABLES REVISING EACH
! ASSOCIATED
! LIST,
! RECALCULATING DEGREES AND REMOVING DUPLICATES.
              limit = nint(fratio*(n-nel))
              DO k = k1, k2
                is = iw(k)
                IF (is/=root) THEN
                  IF (nflg<=2) THEN
! RESET FLAG VALUES TO +/-IOVFLO.
                    DO i = 1, n
                      IF (flag(i)>0) flag(i) = iovflo
                      IF (flag(i)<=-2) flag(i) = -iovflo
                    END DO
                    nflg = iovflo
                  END IF
! REDUCE NFLG BY ONE TO CATER FOR THIS SUPERVARIABLE.
                  nflg = nflg - 1
! BEGIN WITH THE DEGREE OF THE NEW ELEMENT. ITS VARIABLES
! MUST
! ALWAYS
! BE COUNTED DURING THE DEGREE CALCULATION AND THEY ARE
! ALREADY
! FLAGGED WITH THE VALUE 0.
                  id = idn
! RUN THROUGH THE LIST ASSOCIATED WITH SUPERVARIABLE IS
                  kp1 = ipe(is) + 1
! NP POINTS TO THE NEXT ENTRY IN THE REVISED LIST.
                  np = kp1
                  kp2 = iw(kp1-1) + kp1 - 1
                  DO kp = kp1, kp2
                    ke = iw(kp)
! TEST WHETHER KE IS AN ELEMENT, A REDUNDANT ENTRY OR A
! SUPERVARIABLE.
                    IF (flag(ke)==-1) THEN
                      IF (ipe(ke)/=-root) THEN
                        GO TO 60
                      ELSE
! KE has been merged into the root
                        ke = root
                        iw(kp) = root
                        IF (flag(ke)==-1) GO TO 60
                      END IF
                    END IF
                    IF (flag(ke)>=0) THEN
                      GO TO 70
                    ELSE
! SEARCH LIST OF ELEMENT KE, REVISING THE DEGREE WHEN
! NEW
! VARIABLES
! FOUND.
                      jp1 = ipe(ke) + 1
                      jp2 = iw(jp1-1) + jp1 - 1
                      idl = id
                      DO jp = jp1, jp2
                        js = iw(jp)
! JUMP IF JS HAS ALREADY BEEN COUNTED.
                        IF (flag(js)>nflg) THEN
                          id = id + nv(js)
                          flag(js) = nflg
                        END IF
                      END DO
! JUMP IF ONE OR MORE NEW SUPERVARIABLES WERE FOUND.
                      IF (id<=idl) THEN
! CHECK WHETHER EVERY VARIABLE OF ELEMENT KE IS IN NEW
! ELEMENT ME.
                        DO jp = jp1, jp2
                          js = iw(jp)
                          IF (flag(js)/=0) GO TO 50
                        END DO
! SET TREE POINTER AND FLAG TO INDICATE THAT ELEMENT
! KE
! IS ABSORBED
! INTO NEW ELEMENT ME.
                        ipe(ke) = -me
                        flag(ke) = -1
                        GO TO 60
                      END IF
! STORE ELEMENT KE IN THE REVISED LIST FOR SUPERVARIABLE
! IS AND FLAG IT.
50                    iw(np) = ke
                      flag(ke) = -nflg
                      np = np + 1
                    END IF
60                  CONTINUE
                  END DO
                  np0 = np
                  GO TO 90
! TREAT THE REST OF THE LIST ASSOCIATED WITH SUPERVARIABLE
! IS.
! IT
! CONSISTS ENTIRELY OF SUPERVARIABLES.
70                kp0 = kp
                  np0 = np
                  DO kp = kp0, kp2
                    ks = iw(kp)
                    IF (flag(ks)<=nflg) THEN
                      IF (ipe(ks)==-root) THEN
                        ks = root
                        iw(kp) = root
                        IF (flag(ks)<=nflg) GO TO 80
                      ELSE
                        GO TO 80
                      END IF
                    END IF
! ADD TO DEGREE, FLAG SUPERVARIABLE KS AND ADD IT TO NEW
! LIST.
                    id = id + nv(ks)
                    flag(ks) = nflg
                    iw(np) = ks
                    np = np + 1
80                  CONTINUE
                  END DO
! MOVE FIRST SUPERVARIABLE TO END OF LIST, MOVE FIRST
! ELEMENT
! TO END
! OF ELEMENT PART OF LIST AND ADD NEW ELEMENT TO FRONT OF
! LIST.
90                IF (id>=limit) THEN
! Treat IS as full. Merge it into the root node.
                    info%n_dense_rows = info%n_dense_rows + 1
                    IF (nvroot==0) THEN
                      root = is
                      ipe(is) = 0
                    ELSE
                      iw(k) = root
                      ipe(is) = -root
                      nv(root) = nv(root) + nv(is)
                      nv(is) = 0
                      flag(is) = -1
                    END IF
                    nvroot = nv(root)
                  ELSE
                    iw(np) = iw(np0)
                    iw(np0) = iw(kp1)
                    iw(kp1) = me
! STORE THE NEW LENGTH OF THE LIST.
                    iw(kp1-1) = np - kp1 + 1

! CHECK WHETHER ROW IS IS IDENTICAL TO ANOTHER BY LOOKING
! IN
! LINKED
! LIST OF SUPERVARIABLES WITH DEGREE ID AT THOSE WHOSE
! LISTS
! HAVE
! FIRST ENTRY ME. NOTE THAT THOSE CONTAINING ME COME FIRST
! SO THE
! SEARCH CAN BE TERMINATED WHEN A LIST NOT STARTING WITH
! ME
! IS
! FOUND.
                    js = ipd(id)
                    DO l = 1, n
                      IF (js<=0) THEN
                        GO TO 120
                      ELSE
                        kp1 = ipe(js) + 1
                        IF (iw(kp1)/=me) THEN
                          GO TO 120
                        ELSE
! JS HAS SAME DEGREE AND IS ACTIVE. CHECK IF
! IDENTICAL
! TO IS.
                          kp2 = kp1 - 1 + iw(kp1-1)
                          DO kp = kp1, kp2
                            ie = iw(kp)
! JUMP IF IE IS A SUPERVARIABLE OR AN ELEMENT NOT
! IN
! THE LIST OF IS.
                            IF (abs(flag(ie)+0)>nflg) GO TO 100
                          END DO
                          GO TO 110

100                       js = nxt(js)
                        END IF
                      END IF
                    END DO
! SUPERVARIABLE AMALGAMATION. ROW IS IS IDENTICAL TO ROW
! JS.
! REGARD ALL VARIABLES IN THE TWO SUPERVARIABLES AS BEING
! IN
! IS. SET
! TREE POINTER, FLAG AND NV ENTRIES.
110                 ipe(js) = -is
                    nv(is) = nv(is) + nv(js)
                    nv(js) = 0
                    flag(js) = -1
! REPLACE JS BY IS IN LINKED LIST.
                    ns = nxt(js)
                    ls = lst(js)
                    IF (ns>0) lst(ns) = is
                    IF (ls>0) nxt(ls) = is
                    lst(is) = ls
                    nxt(is) = ns
                    lst(js) = 0
                    nxt(js) = 0
                    IF (ipd(id)==js) ipd(id) = is
                    GO TO 130
! INSERT IS INTO LINKED LIST OF SUPERVARIABLES OF SAME
! DEGREE.
120                 ns = ipd(id)
                    IF (ns>0) lst(ns) = is
                    nxt(is) = ns
                    ipd(id) = is
                    lst(is) = -id
                    md = min(md,id)
                  END IF
                END IF
130             CONTINUE
              END DO

! RESET FLAGS FOR SUPERVARIABLES IN NEWLY CREATED ELEMENT AND
! REMOVE THOSE ABSORBED INTO OTHERS.
              DO k = k1, k2
                is = iw(k)
                IF (nv(is)/=0) THEN
                  flag(is) = nflg
                  iw(ip) = is
                  ip = ip + 1
                END IF
              END DO
              iwfr = k1
              flag(me) = -nflg
! MOVE FIRST ENTRY TO END TO MAKE ROOM FOR LENGTH.
              iw(ip) = iw(k1)
              iw(k1) = ip - k1
! SET POINTER FOR NEW ELEMENT AND RESET IWFR.
              ipe(me) = k1
              iwfr = ip + 1
            END IF
          END IF
        END DO

! Absorb any remaining variables into the root
140     DO is = 1, n
          IF (nxt(is)/=0 .OR. lst(is)/=0) THEN
            IF (nvroot==0) THEN
              root = is
              ipe(is) = 0
            ELSE
              ipe(is) = -root
            END IF
            nvroot = nvroot + nv(is)
            nv(is) = 0
          END IF
        END DO
! Link any remaining elements to the root
        DO ie = 1, n
          IF (ipe(ie)>0) ipe(ie) = -root
        END DO
        IF (nvroot>0) nv(root) = nvroot


        DEALLOCATE (flag,ipd,lst,nxt,STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_dealloc
          RETURN
        END IF


END SUBROUTINE mc68_min_deg_anal


SUBROUTINE mc68_min_deg_tree_search(n,ipe,nv,perm,nemin,info)

! Tree search
! Given son to father tree pointers, reorder so that eldest son has
! smallest degree and perform depth-first
! search to find pivot order and number of eliminations
! and assemblies at each stage.
! N must be set to the matrix order. It is not altered.
! IPE(I) must be set equal to -(father of node I) or zero if
! node is a root, if NV(I) > 0. If NV(I) = 0, then I is
! subordinate variable of a supervariable and -IPE(I) points to
! principal variable.  It is altered to point to its next
! younger brother if it has one, but otherwise is not changed.
! NV(I) must be set to zero if variable is a subordinate variable
! of a supervariable and to the degree otherwise.
! PERM is set to the new permutation after dfs of tree.  PERM(I) is
! the position of variable I in the pivot order.

        INTEGER, INTENT (IN) :: n, nemin
        INTEGER, INTENT (INOUT) :: ipe(n), nv(n)
        TYPE (mc68_info), INTENT (INOUT) :: info
        INTEGER, INTENT (OUT) :: perm(n)

! Local variables
! IPS(I) is used temporarily to hold
! -(eldest son of node I) if it has one and 0 otherwise. It is
! finally set to hold the position of node I in the order.
! NE(IS) is set to the number of variables
! eliminated at stage IS of the elimination.
! NA(IS) is set to the number of elements
! assembled at stage IS of the elimination.
! NODE (I) is used during the code
! to hold the number of subordinate variables for variable I and
! on output it holds
! the node (in dfs ordering) at which variable I is eliminated.
! It is also defined for subordinate variables.
! ND(IS) is set to the degree at stage IS of
! the elimination.
! NSTEPS is set to the number of elimination steps.
! NEMIN is used to control the amalgamation process between
! a son and its father (if the number of fully summed
! variables of both nodes is smaller than NEMIN).
! SUBORD(I) holds the first subordinate variable
! for variable I if I
! is a principal variable and holds the next subordinate variable
! if otherwise.  It is zero at the end of the chain.

! .. Array Arguments ..
        INTEGER, ALLOCATABLE, DIMENSION (:) :: fils, frere, ips, na, nd, ne, &
          node, subord
! ..
! .. Scalars ..
        INTEGER :: i, ib, if, ifson, il, in, inb, inf, infs, inl, ino, inos, &
          ins, insw, int, iperm, is, ison, k, l, nr, nr1

! Initialisations

        ALLOCATE (fils(n),frere(n),ips(n),na(n),nd(n),ne(n),node(n),subord(n), &
          STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_alloc
          RETURN
        END IF

        ips(:) = 0
        ne(:) = 0
        node(:) = 0
        subord(:) = 0


! Set IPS(I) to -(eldest son of node I) and IPE(I) to next younger
! brother of node I if it has one.
        nr = n + 1
        DO i = 1, n
          if = -ipe(i)
          IF (nv(i)==0) THEN
! I is a subordinate node, principal variable is IF
            IF (subord(if)/=0) subord(i) = subord(if)
            subord(if) = i
            node(if) = node(if) + 1
! Node IF is the father of node I.
          ELSE IF (if/=0) THEN
! IS is younger brother of node I.
! IPS(IF) will eventually point to - eldest son of IF.
            is = -ips(if)
            IF (is>0) ipe(i) = is
            ips(if) = -i
          ELSE
! I is a root node
            nr = nr - 1
            ne(nr) = i
          END IF
        END DO

! We reorganize the tree so that the eldest son has maximum number of
! variables.  We combine nodes when the number of variables in a son
! is greater than or equal to the number of variables in the father.
! If the eldest son has the maximum number of variables,
! and if a combination is possible, it has to be possible with
! the eldest son.
! FILS is just used as workspace during this reorganization and is
! reset
! afterwards.

        DO i = 1, n
          fils(i) = ips(i)
        END DO
        nr1 = nr
        ins = 0
10      CONTINUE
! Jump if all roots processed.
        IF (nr1<=n) THEN
! Get next root
          ins = ne(nr1)
          nr1 = nr1 + 1
20        CONTINUE
! Depth first search through eldest sons.
          inl = fils(ins)
          IF (inl<0) THEN
            ins = -inl
            GO TO 20
          ELSE
            DO WHILE (ipe(ins)<0)
! INS is youngest son otherwise IPE value would be positive.
              ins = -ipe(ins)
! INS is now the father of the reorganized son so we can
! clear the pointer to the sons.
              fils(ins) = 0
! Continue backtracking until we encounter node with younger
! brother.
            END DO

            IF (ipe(ins)/=0) THEN
! INB is younger brother of INS.
              inb = ipe(ins)
              IF (nv(inb)>=nv(ins)) THEN
                ins = inb
! Do depth first search from younger brother
              ELSE
! Exchange INB and INS
! Find previous brother of INS (could be the father)
! then we do depth first search with INS = INB
                inf = inb
30              CONTINUE
                inf = ipe(inf)
                IF (inf>0) GO TO 30
! -INF IS THE FATHER
                inf = -inf
                infs = -fils(inf)
! INFS is eldest son of INF
                IF (infs==ins) THEN
! INS is eldest brother .. a role which INB now assumes
                  fils(inf) = -inb
                  ips(inf) = -inb
                  ipe(ins) = ipe(inb)
                  ipe(inb) = ins
                ELSE
                  insw = infs
40                CONTINUE
                  infs = ipe(insw)
                  IF (infs/=ins) THEN
                    insw = infs
                    GO TO 40
                  END IF
                  ipe(ins) = ipe(inb)
                  ipe(inb) = ins
                  ipe(insw) = inb
                END IF
                ins = inb
! Depth first search from moved younger brother
              END IF
              GO TO 20
            END IF
          END IF
! INS is a root, check for next one.
          ins = 0
          GO TO 10
        END IF
! Set FRERE and FILS
        DO i = 1, n
          frere(i) = ipe(i)
          fils(i) = ips(i)
        END DO

! Depth-first search.
! IL holds the current tree level. Roots are at level N, their sons
! are at level N-1, etc.
! IS holds the current elimination stage. We accumulate the number
! of eliminations at stage is directly in NE(IS). The number of
! assemblies is accumulated temporarily in NA(IL), for tree
! level IL, and is transferred to NA(IS) when we reach the
! appropriate stage IS.
        is = 1
! I is the current node.
        i = 0
! IPERM is used as pointer to setting permutation vector
        iperm = 1
        DO k = 1, n
          IF (i<=0) THEN
! Pick up next root.
! Stop if all roots used (needed because of subordinate variables)
            IF (nr>n) THEN
              GO TO 130
            ELSE
              i = ne(nr)
              ne(nr) = 0
              nr = nr + 1
              il = n
              na(n) = 0
            END IF
          END IF
! Go to son for as long as possible, clearing father-son pointers
! in IPS as each is used and setting NA(IL)=0 for all levels
! reached.
          DO l = 1, n
            IF (ips(i)>=0) THEN
              GO TO 50
            ELSE
              ison = -ips(i)
              ips(i) = 0
              i = ison
              il = il - 1
              na(il) = 0
            END IF
          END DO
! Record position of node I in the order.
50        ips(i) = k
! Add number of subordinate variables to variable I
          ne(is) = ne(is) + node(i) + 1
          IF (il<n) na(il+1) = na(il+1) + 1
          na(is) = na(il)
          nd(is) = nv(i)
          node(i) = is
          perm(i) = iperm
          iperm = iperm + 1
! Order subordinate variables to node I
          in = i
60        CONTINUE
          IF (subord(in)/=0) THEN
            in = subord(in)
            node(in) = is
            perm(in) = iperm
            iperm = iperm + 1
            GO TO 60
          END IF
! Check for static condensation
          IF (na(is)==1) THEN
            IF (nd(is-1)-ne(is-1)==nd(is)) GO TO 70
          END IF
! Check for small numbers of eliminations in both last two steps.
          IF (ne(is)<nemin) THEN
            IF (na(is)/=0) THEN
              IF (ne(is-1)<nemin) GO TO 70
            END IF
          END IF
          is = is + 1
          GO TO 120

! Combine the last two steps
70        na(is-1) = na(is-1) + na(is) - 1
          nd(is-1) = nd(is) + ne(is-1)
          ne(is-1) = ne(is) + ne(is-1)
          ne(is) = 0
          node(i) = is - 1
! Find eldest son (IFSON) of node I (IS)
! Note that node I must have a son (node IS-1 is youngest)
          ifson = -fils(i)
! Now find youngest son INO (he is node IS-1)
          in = ifson
80        CONTINUE
          ino = in
          in = frere(in)
          IF (in>0) GO TO 80
! Cannot be root node .. so points to father
! Merge node IS-1 (INO) into node IS (I)
          nv(ino) = 0
! IPE already set .. was father pointer now principal variable
! pointer
! Now make subsidiary nodes of INO into subsidiary nodes of I.
! Subordinate nodes of INO become subordinate nodes of I
          in = i
90        CONTINUE
          IF (subord(in)/=0) THEN
            in = subord(in)
            node(in) = is - 1
            GO TO 90
          END IF
          subord(in) = ino
          in = ino
          IF (subord(in)/=0) THEN
            in = subord(in)
            ipe(in) = -i
          END IF

! INOS is eldest son of INO
          inos = -fils(ino)

! Find elder brother of node INO
! First check to see if he is only son
          IF (ifson/=ino) THEN
            in = ifson
100         CONTINUE
            ins = in
            in = frere(in)
            IF (in/=ino) GO TO 100
! INS is older brother .. make him brother of first son of INO (ie
! INOS)
! and  make INOS point to I now as father.
! Jump if there is no son of INO
            IF (inos==0) THEN
! Elder brother of INO just points to (new) father.
              frere(ins) = -i
              GO TO 120
            ELSE
              frere(ins) = inos
            END IF
          END IF
! INT is youngest brother of INOS.  Make him point to (new) father.
          in = inos
          IF (in/=0) THEN
110         CONTINUE
            int = in
            in = frere(in)
            IF (in>0) GO TO 110
            frere(int) = -i
          END IF
120       ib = ipe(i)
          IF (ib>=0) THEN
! Node I has a younger brother or is a root
            IF (ib>0) na(il) = 0
            i = ib
          ELSE
! I has no brothers. Go to father of node I
            i = -ib
            il = il + 1
          END IF
        END DO
130     continue


        DEALLOCATE (fils,frere,ips,na,nd,ne,node,subord,STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_dealloc
          RETURN
        END IF

END SUBROUTINE mc68_min_deg_tree_search

SUBROUTINE mc68_compress(n,ipe,iw,lw,iwfr,ncmpa)
! Is identical to subroutine MA27UD.  Internal version for MA57.
! COMPRESS LISTS HELD BY MA27H/HD (MA57H/HD) IN IW AND ADJUST POINTERS
! IN IPE TO CORRESPOND.
! N IS THE MATRIX ORDER. IT IS NOT ALTERED.
! IPE(I) POINTS TO THE POSITION IN IW OF THE START OF LIST I OR IS
! ZERO IF THERE IS NO LIST I. ON EXIT IT POINTS TO THE NEW POSITION.
! IW HOLDS THE LISTS, EACH HEADED BY ITS LENGTH. ON OUTPUT THE SAME
! LISTS ARE HELD, BUT THEY ARE NOW COMPRESSED TOGETHER.
! LW HOLDS THE LENGTH OF IW. IT IS NOT ALTERED.
! IWFR NEED NOT BE SET ON ENTRY. ON EXIT IT POINTS TO THE FIRST FREE
! LOCATION IN IW.
! ON RETURN IT IS SET TO THE FIRST FREE LOCATION IN IW.
! NCMPA is number of compresses.

! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: lw, n
        INTEGER, INTENT (OUT) :: iwfr
        INTEGER, INTENT (INOUT) :: ncmpa
! ..
! .. Array Arguments ..
        INTEGER, INTENT (INOUT) :: ipe(n), iw(lw)
! ..
! .. Local Scalars ..
        INTEGER i, ir, k, k1, k2, lwfr
! ..
        ncmpa = ncmpa + 1
! PREPARE FOR COMPRESSING BY STORING THE LENGTHS OF THE
! LISTS IN IPE AND SETTING THE FIRST ENTRY OF EACH LIST TO
! -(LIST NUMBER).
        DO i = 1, n
          k1 = ipe(i)
          IF (k1>0) THEN
            ipe(i) = iw(k1)
            iw(k1) = -i
          END IF
        END DO

! COMPRESS
! IWFR POINTS JUST BEYOND THE END OF THE COMPRESSED FILE.
! LWFR POINTS JUST BEYOND THE END OF THE UNCOMPRESSED FILE.
        iwfr = 1
        lwfr = iwfr
        DO ir = 1, n
          IF (lwfr>lw) THEN
            GO TO 20
          ELSE
! SEARCH FOR THE NEXT NEGATIVE ENTRY.
            DO k = lwfr, lw
              IF (iw(k)<0) GO TO 10
            END DO
            GO TO 20
! PICK UP ENTRY NUMBER, STORE LENGTH IN NEW POSITION, SET NEW
! POINTER
! AND PREPARE TO COPY LIST.
10          i = -iw(k)
            iw(iwfr) = ipe(i)
            ipe(i) = iwfr
            k1 = k + 1
            k2 = k + iw(iwfr)
            iwfr = iwfr + 1
            IF (k1<=k2) THEN
! COPY LIST TO NEW POSITION.
              DO k = k1, k2
                iw(iwfr) = iw(k)
                iwfr = iwfr + 1
              END DO
            END IF
            lwfr = k2 + 1
          END IF
        END DO
20      RETURN

END SUBROUTINE mc68_compress
! -----------------------------------------------------------

SUBROUTINE mc68_ma47_analyse(n,ipe,iw,lw,iwfr,cntl,count,nv,next,info)
        INTEGER, PARAMETER :: wp = kind(0.0D0)

! Analysis subroutine

! Given representation of the whole matrix, perform Markowitz
! ordering,
! constructing tree pointers.  It works with supervariables which
! are collections of one or more variables whose rows have the same
! pattern, starting with supervariable I containing variable I for I
! = 1,2,...,N.  Each supervariable has as numerical name that of one
! of its variables (its principal variable). When a supervariable is
! eliminated, its name is used for the name of the element created.


! N must be set to the matrix order. It is not altered.
! Constraint: N <= HUGE/3.
! IPE(I) must be set to the position in IW of the list for row I,
! I=1,N.
! During execution, IPE(I) is used as follows:
! if I is a supervariable, IPE(I) is the position in IW of its list;
! if I is an element, IPE(I) is the position in IW of its list, or
! 0 if the list has zero length;
! if I is a variable that has been absorbed in a supervariable,
! IPE(I) = 0.
! IW must be set on entry to hold lists of entries by rows, each list
! being headed by its length. The lists include any diagonal entries
! and the entries of both the upper and lower triangular parts.
! During execution, it holds lists for supervariables and lists for
! elements. A list for a supervariable represents its row of the
! reduced matrix and has the form:
! length ( = total number of elements and supervariables)
! list of element parts in the form JS+N*PART, where JS is an
! element name and PART has the value
! 2 for the part with leading zeros;
! 0 for the full part;
! 1 for the part with trailing zeros;
! list of supervariables.
! A list for an element has the form:
! length ( = 2 + total number of supervariables)
! number of supervariables in the part with leading zeros
! number of supervariables in the part with trailing zeros
! list of supervariables.
! LW must be set to the length of IW. It may be altered.
! IWFR must be set to the position in IW of the first unused location.
! It is revised during execution and continues to have
! this meaning.
! COUNT need not be set on entry. For supervariables, it is used to
! hold the row counts in the reduced matrix. COUNT(I) is negated
! if supervariable I has been eliminated as the trailing
! part of a tile or oxo pivot.
! NV must be set on entry so that NV(I) = 1 for a nondefective
! variable
! (nonzero diagonal entry) and NV(I) = -1 otherwise, I = 1,..., N.
! During execution, NV(JS) is used as follows:
! if JS is a variable that has been absorbed into another,
! NV(JS) = 0;
! if JS is a supervariable, NV(JS) holds its number of
! variables, negated for a defective supervariable.
! NEXT need not be set. During execution, it is used as follows:
! if JS is a supervariable that has not been eliminated or
! absorbed, NEXT(JS) is the next such supervariable having the
! same row count, or zero if JS is last in its list;
! if supervariable IS has been absorbed and JS is the next
! variable in its supervariable, NEXT(IS)=-JS;
! if supervariable IS was eliminated in a block pivot with
! supervariable JS, NEXT(IS)=-JS;
! if IE is an element whose parent in the tree is JE, NEXT(IE)=-JE;
! if IE is an element with no parent in the tree, NEXT(IE)=0.
! NCMP is number of times garbage collection routine is called.
! ICNTL must be set by the user as follows and is not altered.
! ICNTL(1)  must be set to the stream number for error messages.
! A value less than 1 suppresses output.
! ICNTL(2) must be set to the stream number for diagnostic output.
! A value less than 1 suppresses output.
! ICNTL(3) must be set to control the amount of output:
! 0 None.
! 1 Error messages only.
! 2 Error and warning messages.
! 3 As 2, plus scalar parameters and a few entries of array
! parameters on entry and exit.
! 4  As 2, plus all parameters on entry and exit.
! CNTL must set to one if the user wants the pivot order chosen
! by Markowitz strategy,, and greater than 1 for Markowits strategy
! with
! each search for a structured pivot limited to this number of rows.
! INFO is left unaltered except
! 8 is added to INFO%flag and
! INFO%n_zero_eigs is set to the number of zero eigenvalues found.

! Local constants

! Local variables
! BOUND True if row count is to be bounded rather than recomputed.
! CMIN holds the minimum Markowitz cost of a potential structured
! pivot
! found so far.
! COST holds the Markowitz cost of a potential structured pivot.
! FLG Flag value.
! FLG1 Flag value.
! I Temporary variable.
! IE is an element.
! IEL is an element.
! IP is used to point to a list of elements or supervariables.
! IR is used for row counts.
! IRL last value of IR
! IS is a supervariable.
! ITHR Index for do loop in case threshold needs to be raised.
! JP is used to point to a list of elements or supervariables.
! JP1 start of range for JP.
! JP2 end of range for JP.
! JP3 start of range for JP.
! JP4 end of range for JP.
! JS is a supervariable.
! K is used to point to a list of elements or supervariables.
! KE is an element.
! KIND has the value
! 1 for a full pivot,
! 2 for a tile pivot, or
! 3 for an oxo pivot
! KP is used to point to a list of elements or supervariables.
! KP1 start of range for KP.
! KP2 end of range for KP.
! KS is a supervariable.
! K1 start of range for K.
! K2 end of range for K.
! LIST is used to search a list.
! LOOP Row of the pivot when constructing the index list of a new
! element.
! LS is a supervariable.
! ME is the element created by the current pivot step.
! MINF is lower bound for the row count of a nondefective variable.
! MINR is used for the current minimum row count.
! ML is used for the main loop index.
! MP is stream for warning and diagnostic messages.
! MROWS holds the max. number of rows to be searched for a structured
! pivot.
! MS is a supervariable.
! NEL is used to hold the number of variables that have been
! eliminated.
! NFLG is used for the current flag value in array FLAG. It starts
! with the value N*3 and is reduced by 1 each time it is used
! until it has the value 4 when it is reset to the value N*3.
! NP is used to point to a list of elements or supervariables.
! NP0 is used to point to a list of elements or supervariables.
! NR is used for the row counts in the generated element:
! NR(1) for a row of the middle block
! NR(2) for a row of the last block
! NR(3) for a row of the first block
! NROWS holds the number of rows searched for a structured pivot.
! NS is a supervariable.
! NSC is used for supervariable counts:
! NSC(1) first zero block
! NSC(2) full block
! NSC(3) second zero block
! NSVARS Supervariable counts.
! NVC Number of variables in both pivot rows, excluding the pivot.
! NVPIV holds the number of variables in the pivot (each half in the
! case of a tile or oxo).
! PART Part of element: 0 for full part, 1 for trailing part,
! 2 for leading part.
! PIV holds the name(s) of the pivot supervariable(s), ordered in the
! case of a tile pivot so that the first diagonal entry is nonzero.
! PIVOT As PIV, except that this holds the original names when one
! has been split.
! PIVT In the case of a structured pivot with unbalanced
! supervariables,
! PIVT is equal to whichever of PIVOT(1) and PIVOT(2) has more
! variables. Zero otherwise.
! RANK Upper bound for the rank.
! SHIFT Movement of current list after compress.
! THRESH Row counts greater than THRESH are lower bounds (fill-ins may
! have been ignored.
! LAST need not be set. During execution, it is used as follows:
! if JS is a supervariable that has not been eliminated or
! absorbed, LAST(JS) is the previous such supervariable having
! the same row count, or zero if JS is first in its list;
! otherwise, if NEXT(JS) < 0, LAST(JS) is one of the values
! NEXT(JS), NEXT(-NEXT(JS)), NEXT(-NEXT(-NEXT(JS))), ...
! IPR need not be set. During execution IPR(IR) is the first
! supervariable with row count IR or zero if there are none.
! FLAG is used as workspace for element and supervariable flags.
! if IS is a pivotal supervariable, FLAG(IS) = -1;
! if IS is a non-pivotal supervariable, FLAG(IS) >= 0;
! if IS is a supervariable involved in the current pivot row (other
! than in the pivot):
! for a full pivot:           FLAG(IS) = 1
! for a tile or oxo pivot:
! if in both rows,       FLAG(IS) = 0
! if in first row only,  FLAG(IS) = 1
! if in second row only, FLAG(IS) = 2
! if IE is an element, FLAG(IE) < -1;
! if I is neither a supervariable nor an element, FLAG(I) = -1.
! LEAF need not be set. During execution, if IS is a supervariable,
! its
! final variable is LEAF(IS). LEAF(IS) is linked to IS through NEXT
! in a chain that covers the other variables of the supervariable.
! For an element IE to which the current supervariable IS belongs,
! LEAF(IE) holds the part to which the supervariable belongs.
! SVARS is used as workspace for the supervariables of the new
! element.
! .. Parameters ..
        REAL (wp) :: zero
        PARAMETER (zero=0E0_wp)
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: n, cntl
        INTEGER, INTENT (INOUT) :: iwfr
        INTEGER (long), INTENT (INOUT) :: lw

! ..
! .. Array Arguments ..
        INTEGER, INTENT (INOUT) :: ipe(n), nv(n)
        INTEGER, DIMENSION (:), ALLOCATABLE, INTENT (INOUT) :: iw
        INTEGER, INTENT (OUT) :: count(n), next(n)

        TYPE (mc68_info), INTENT (INOUT) :: info
! ..
! .. Local Scalars ..

        REAL (wp) :: cmin, cost

        INTEGER :: flg, flg1, ie, iel, ip, ir, irl, is, ithr, jp, jp1, jp2, &
          jp3, jp4, js, k, k1, k2, ke, kind1, kp, kp1, kp2, ks, list, loop, &
          ls, me, minf, minr, ml, mrows, ms, nel, nflg, np, np0, nrows, &
          ns, nsvars, nvc, nvpiv, part, pivt, rank, shift, thresh
        INTEGER (long) :: szw
        LOGICAL :: bound
! ..
! .. Local Arrays ..
        INTEGER :: nr(3), nsc(3), piv(2), pivot(2), icntl4
        INTEGER, ALLOCATABLE, DIMENSION (:) :: flag, ipr, last, leaf, svars

        TYPE (zb01_info) :: info_zb01
! ..
! .. Intrinsic Functions ..
        INTRINSIC abs, max, min, real, sign, sqrt

        IF (cntl<=1) THEN
          icntl4 = 0
        ELSE
          icntl4 = cntl
        END IF
        mrows = n
        IF (icntl4>1) mrows = icntl4

        ALLOCATE (flag(n),ipr(n),last(n),leaf(n),svars(n),STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_alloc
          RETURN
        END IF

! Create supervariables from variables with identical rows
        CALL mc68_ma47_merge(n,ipe,iw,lw,nv,next,last,leaf,flag,count,svars)

! Initializations

        thresh = int( max(sqrt(real(n,kind=wp)),3._wp) )
        rank = n
        minf = 1
        minr = 1
        info%n_compressions = 0
        nflg = n*3
        nel = 0
        DO is = 1, n
          ipr(is) = 0
          flag(is) = nflg
        END DO

! Link together supervariables having same row count
        DO is = 1, n
          k = ipe(is)
          IF (k==0) THEN
! Row IS is identical to another row and has been absorbed.
            flag(is) = -1
          ELSE
            ir = iw(k)
            IF (ir==0) THEN
              rank = rank + nv(is)
! Treat zero block row as if it had nonzero diagonal entries.
              nv(is) = -nv(is)
              ir = 1
            END IF
! Store the row count and add the row to its linked list.
            count(is) = ir
            ns = ipr(ir)
            IF (ns>0) last(ns) = is
            next(is) = ns
            ipr(ir) = is
            last(is) = 0
          END IF
        END DO

! Start of main loop
        DO ml = 1, n
! Leave loop if all variables have been eliminated.
          IF (nel>=n) THEN
            GO TO 280
          ELSE
! Find minimum row count
            ir = minr
            DO minr = ir, n
              IF (ipr(minr)/=0) GO TO 10
            END DO
! Find next pivot.
10          DO ithr = 1, n
! Outer loop on ITHR needed in case of the threshold needing to
! be
! increased.
              nrows = 0
              cmin = real(n,kind=wp)**2
              DO ir = minr, n
                IF (minf<=ir) THEN
! Look for a full pivot
                  ms = ipr(ir)
                  DO list = 1, n
                    IF (ms==0) THEN
                      GO TO 20
! If this is a full pivot, accept it or raise the
! threshold and try
! again.
                    ELSE IF (nv(ms)>0) THEN
                      GO TO 60
                    ELSE
                      ms = next(ms)
                    END IF
                  END DO
                  GO TO 30
20                minf = ir + 1
                END IF
! Look for a tile or oxo pivot with MS as its defective
! variable.
30              ms = ipr(ir)
                DO list = 1, n
                  IF (ms==0) THEN
                    GO TO 50
                  ELSE
! We need not look for a tile with MS as its non-defective
! variable
! because if the row count of its mate is IS, its M-cost
! is
! (IR+IS-3)*(IS-1)
! >= (IS-1)**2 since IR >= 2
! >= (IR-1)**2 since IS >= IR
! Reduce NFLG by one to cater for this supervariable.
                    IF (nflg<=4) CALL mc68_ma47_reset(n,flag,nflg)
                    nrows = nrows + 1
                    nflg = nflg - 1
                    kp = ipe(ms)
                    kp1 = kp + 1
                    kp2 = kp + iw(kp)
                    DO kp = kp1, kp2
                      part = (iw(kp)-1)/n
                      ke = iw(kp) - part*n
                      IF (flag(ke)/=-1) THEN
                        IF (flag(ke)<=-2) THEN
! KE is an element touched by MS.
                          jp = ipe(ke)
                          IF (part==2) THEN
! Supervariable in leading part of element
                            jp1 = jp + 3 + iw(jp+1)
                            jp2 = jp + iw(jp)
                          ELSE
! Supervariable in trailing part of element
                            jp1 = jp + 3
                            jp2 = jp + iw(jp) - iw(jp+2)
                          END IF
                        ELSE
! Supervariable MS is defective so cannot lie in
! full
! part of element
! We have reached the list of variables
                          jp1 = kp
                          jp2 = kp2
                        END IF
! Search the variable list.
                        DO jp = jp1, jp2
                          is = iw(jp)
                          IF (flag(is)>nflg) THEN
! This is the first time IS has been encountered
! in
! this search.
                            flag(is) = nflg
                            IF (nv(is)<0) THEN
! Potential oxo pivot
                              cost = real(count(is)-1,kind=wp)* &
                                real(ir-1,kind=wp)
                            ELSE
! Potential tile pivot
                              cost = real(ir+count(is)-3,kind=wp)* &
                                real(ir-1,kind=wp)
                            END IF

                            IF (cost<cmin) THEN
! This is the best potential pivot so far found.
                              cmin = cost
                              pivot(1) = is
                              pivot(2) = ms
! Exit loop no better pivot is possible.
                              IF (cmin<=real(ir-1,kind=wp)**2) GO TO 70
                            END IF
                          END IF
                        END DO

! Exit loop if we have searched the list of variables.
                        IF (jp2==kp2) GO TO 40
                      END IF
                    END DO
! Exit loop if enough rows have been searched.
40                  IF (nrows>=mrows) THEN
                      GO TO 70
                    ELSE
                      ms = next(ms)
                    END IF
                  END IF
                END DO
! Exit loop no better pivot is possible.
50              IF (cmin<=real(ir,kind=wp)**2) GO TO 70
              END DO
              GO TO 70
60            IF (ir<=thresh) THEN
                GO TO 90
              ELSE
                GO TO 80
              END IF
70            ir = max(count(pivot(1)),count(pivot(2)))
              IF (ir<=thresh) GO TO 100
! Revise the threshold and repeat the choice of pivot.
80            CALL mc68_ma47_thresh(thresh,ir+n/10,n,ipe,iw,lw,count,nv,next, &
                last,ipr,flag,nflg)
            END DO

! Accept full pivot
90          kind1 = 1
            me = ms
            pivot(1) = ms
            pivot(2) = ms
            pivt = 0
            nvpiv = nv(ms)
            cmin = real(ir-1,kind=wp)**2
            flag(ms) = -1
            nel = nel + nvpiv
            GO TO 110

! Accept tile or oxo pivot
! PIVOT(2) cannot have greater row count than PIVOT(1).
100         kind1 = 2
            IF (nv(pivot(1))<0) kind1 = 3
            flag(pivot(1)) = -1
            flag(pivot(2)) = -1
            piv(1) = pivot(1)
            piv(2) = pivot(2)
            nvpiv = abs(nv(pivot(1)))
            IF (nvpiv==abs(nv(pivot(2)))) THEN
              pivt = 0
            ELSE
              IF (nvpiv>abs(nv(pivot(2)))) THEN
                pivt = pivot(1)
                nvpiv = abs(nv(pivot(2)))
              ELSE
                pivt = pivot(2)
              END IF
! Split PIVT
              ms = leaf(pivt)
              DO k = 2, nvpiv
                ms = -next(ms)
              END DO
              leaf(ms) = leaf(pivt)
              leaf(pivt) = -next(ms)
! PIVT remains as an active supervariable, but with less
! variables.
              flag(pivt) = nflg
              next(ms) = 0
              nv(pivt) = sign(abs(nv(pivt))-nvpiv,nv(pivt))
              nv(ms) = sign(nvpiv,nv(pivt))
              count(ms) = count(pivt)
              IF (pivt==pivot(1)) THEN
                piv(1) = ms
              ELSE
                piv(2) = ms
              END IF
            END IF

            me = piv(1)
            nsc(2) = 0
            nvc = 0
            nel = nel + nvpiv*2

! Find variables of new element
110         nsvars = 0
            ir = 0
            DO loop = min(kind1,2), 1, -1
              ms = pivot(loop)
              nsc(1) = nsvars
              nr(1) = ir
! Remove chosen variable from linked list (unless it has been
! split).
              IF (ms/=pivt) THEN
                ns = next(ms)
                ls = last(ms)
                next(ms) = 0
                IF (ns>0) last(ns) = ls
                IF (ls>0) THEN
                  next(ls) = ns
                ELSE
                  ipr(count(ms)) = ns
                END IF
              END IF

! Run through the list of the pivotal supervariable, setting
! tree
! pointers and constructing new list of supervariables.
! KP is a pointer to the current position in the old list.
              kp = ipe(ms)
              kp1 = kp + 1
              kp2 = kp + iw(kp)
              DO kp = kp1, kp2
                part = (iw(kp)-1)/n
                ke = iw(kp) - part*n
                IF (flag(ke)/=-1) THEN
                  IF (flag(ke)<=-2) THEN
! KE is an element.
! Link KE or its ancestor to tree if not already linked.
                    ie = ke
                    DO list = 1, n
                      IF (next(ie)==0) THEN
                        GO TO 120
                      ELSE
                        iel = ie
                        ie = -last(ie)
                        IF (ie==me) THEN
                          GO TO 130
                        ELSE
                          last(iel) = -me
                        END IF
                      END IF
                    END DO
120                 next(ie) = -me
                    last(ie) = -me
! Find the relevant part of the list of element KE.
130                 jp = ipe(ke)
                    jp1 = jp + 3
                    jp2 = jp + iw(jp)
                    IF (part/=0) THEN
! Supervariable in full part of element
                      IF (part==2) THEN
! Supervariable in leading part of element
                        jp1 = jp1 + iw(jp+1)
                      ELSE
! Supervariable in trailing part of element
                        jp2 = jp2 - iw(jp+2)
                      END IF
                    END IF
                  ELSE
! We have reached the list of variables
                    jp1 = kp
                    jp2 = kp2
                  END IF

! Search for different supervariables and add them to the
! new
! list.
! This loop is executed once for each element in the list
! and
! once
! for all the supervariables in the list.
                  DO jp = jp1, jp2
                    is = iw(jp)
                    IF (flag(is)>loop) THEN
! IS is not a supervariable or has already been counted.
                      IF (flag(is)==2) THEN
! Supervariable in both rows of the block pivot
                        nvc = nvc + abs(nv(is))
                        nsc(2) = nsc(2) + 1
                        flag(is) = 0
                        count(is) = count(is) - nvpiv
                      ELSE
! New supervariable
                        ir = ir + abs(nv(is))
                        flag(is) = loop
! Store IS in new list.
                        nsvars = nsvars + 1
                        svars(nsvars) = is
! Remove IS from row count linked list
                        ls = last(is)
                        last(is) = 0
                        ns = next(is)
                        next(is) = 0
                        IF (ns>0) last(ns) = ls
                        IF (ls>0) THEN
                          next(ls) = ns
                        ELSE
                          ipr(count(is)) = ns
                        END IF
                        count(is) = count(is) - nvpiv
                      END IF
                    END IF
                  END DO

                  IF (jp2/=kp2 .AND. loop==1) THEN
! See if the pivot is covered by element KE. If so,
! element
! KE will
! be covered by the new element and is no longer needed.
                    IF (kind1/=1) THEN
                      DO jp = jp1, jp2
                        IF (iw(jp)==pivot(2)) GO TO 140
                      END DO
                      GO TO 150
140                   flag(ke) = -1
                      GO TO 160
                    ELSE IF (part==0) THEN
                      flag(ke) = -1
                    END IF
150                 CONTINUE
                  END IF
                END IF
160             CONTINUE
              END DO
            END DO

            IF (kind1==1) THEN
! Complete calculation of counts for a full pivot.
              nsc(2) = nsvars
              nr(1) = ir
              nsc(3) = 0
              nr(3) = 0
            ELSE
! Link PIV(2) into tree.
              next(piv(2)) = -me
              last(piv(2)) = -me
              count(piv(2)) = -count(piv(2))
              IF (kind1==2) THEN
! Complete calculation of counts for a tile pivot.
                nsc(3) = nsvars - nsc(1)
                nr(3) = ir - nr(1)
                nsc(2) = nsc(1)
                nr(2) = nr(1)
                nsc(1) = 0
                nr(1) = ir
              ELSE
! Complete calculation of counts for an oxo pivot.
                nsc(3) = nsvars - nsc(1)
                nr(3) = ir - nr(1) + nvc
                nsc(1) = nsc(1) - nsc(2)
                nr(2) = nr(1)
                nr(1) = ir
! Bring the variables of the leading part to the front.
                k1 = 1
                DO k = 1, nsc(1) + nsc(2)
                  IF (flag(svars(k))==2) THEN
                    ks = svars(k)
                    svars(k) = svars(k1)
                    svars(k1) = ks
                    k1 = k1 + 1
                  END IF
                END DO
              END IF
            END IF

! Deal with the case where the list of the new element has zero
! length.
            IF (nsvars==0) THEN
              ipe(me) = 0
            ELSE
! Run through new list of supervariables looking for elements
! that
! may be absorbed into the new element
! Reduce NFLG by one to cater for these tests.
              IF (nflg<=4) CALL mc68_ma47_reset(n,flag,nflg)
              nflg = nflg - 1
              DO k = 1, nsvars
                is = svars(k)
! Run through the list associated with supervariable IS
                kp = ipe(is)
                kp1 = kp + 1
                kp2 = iw(kp) + kp
                DO kp = kp1, kp2
                  part = (iw(kp)-1)/n
                  ke = iw(kp) - n*part
! Exit if we have reached the list of supervariables
                  IF (flag(ke)>=0) THEN
                    GO TO 190
! Cycle if this is a dummy or has already been checked
                  ELSE IF (flag(ke)/=-1) THEN
                    IF (flag(ke)/=-nflg) THEN
                      flag(ke) = -nflg
                      jp = ipe(ke)
                      jp1 = jp + 3
                      jp2 = jp1 + iw(jp+1)
                      jp4 = jp + iw(jp)
                      jp3 = jp4 - iw(jp+2)
                      IF (kind1==1) THEN
! Check for inclusion in a new full element
                        DO jp = jp1, jp4
                          IF (flag(iw(jp))>2) GO TO 180
                        END DO

                      ELSE IF (kind1==2) THEN
! Check for inclusion in a new tile element.
! Check the full part
                        DO jp = jp2, jp3
                          IF (flag(iw(jp))>2) THEN
                            GO TO 180
                          ELSE IF (flag(iw(jp))==1) THEN
                            GO TO 180
                          END IF
                        END DO
                        flg1 = 0
! Check the trailing part
                        DO jp = jp3 + 1, jp4
                          flg = flag(iw(jp))
                          IF (flg>2) THEN
                            GO TO 180
                          ELSE IF (flg>0) THEN
                            flg1 = flg
                          END IF
                        END DO
! Check the leading part
                        DO jp = jp1, jp2 - 1
                          flg = flag(iw(jp))
                          IF (flg>2) THEN
                            GO TO 180
                          ELSE IF (flg>0) THEN
                            IF (flg1/=0) GO TO 180
                          END IF
                        END DO
                      ELSE

! Check for inclusion in a new oxo element.
! Check the full part
                        DO jp = jp2, jp3
                          IF (flag(iw(jp))>0) GO TO 180
                        END DO
                        flg1 = 0
! Check the trailing part
                        DO jp = jp3 + 1, jp4
                          flg = flag(iw(jp))
                          IF (flg>2) THEN
                            GO TO 180
                          ELSE IF (flg>0) THEN
                            IF (flg1==0) flg1 = flg
                            IF (flg/=flg1) GO TO 180
                          END IF
                        END DO
                        flg1 = 3 - flg1
! Check the leading part
                        DO jp = jp1, jp2 - 1
                          flg = flag(iw(jp))
                          IF (flg>2) THEN
                            GO TO 180
                          ELSE IF (flg>0) THEN
                            IF (flg1==3) flg1 = flg
                            IF (flg/=flg1) GO TO 180
                          END IF
                        END DO
                      END IF
! Element KE is absorbed into new element ME.
                      flag(ke) = -1
! Link KE or its ancestor to tree if not already linked.
                      ie = ke
                      DO list = 1, n
                        IF (next(ie)==0) THEN
                          GO TO 170
                        ELSE
                          iel = ie
                          ie = -last(ie)
                          IF (ie==me) THEN
                            GO TO 180
                          ELSE
                            last(iel) = -me
                          END IF
                        END IF
                      END DO
170                   next(ie) = -me
                      last(ie) = -me
                    END IF
                  END IF
180               CONTINUE
                END DO
190             CONTINUE
              END DO

! Run through new list of supervariables revising each
! associated
! list,
! recalculating row counts and removing duplicates.
              DO loop = 1, kind1
! Find the range to be treated and set flags appropriately.
                IF (loop==1) THEN
! Treat middle block
                  k1 = 1 + nsc(1)
                  k2 = k1 + nsc(2) - 1

                ELSE IF (loop==2) THEN
! Treat last block
                  k1 = k2 + 1
                  k2 = nsvars
                  DO k = k1, k2
                    is = svars(k)
                    IF (nv(is)/=0) flag(is) = nflg
                  END DO
                ELSE

! Treat first block
                  DO k = k1, k2
                    is = svars(k)
                    IF (nv(is)/=0) flag(is) = 1
                  END DO
                  k1 = 1
                  k2 = nsc(1)
                  DO k = k1, k2
                    is = svars(k)
                    IF (nv(is)/=0) flag(is) = nflg
                  END DO
                END IF
! Run through the list of supervariables.
                DO k = k1, k2
                  is = svars(k)
                  bound = count(is) > thresh
                  IF (loop==1) nv(is) = abs(nv(is))
! Reduce NFLG by one to cater for this supervariable.
                  IF (nflg<=4) CALL mc68_ma47_reset(n,flag,nflg)
                  nflg = nflg - 1
! Begin with the row count of the new element. Its variables
! must always be counted during the row count calculation
! and they are already flagged with the value 0, 1, or 2.
                  ir = nr(loop)
! Run through the list associated with supervariable IS
                  kp = ipe(is)
! NP points to the next entry in the revised list.
                  np = kp + 1
                  kp1 = kp + 1
                  kp2 = iw(kp) + kp
                  DO kp = kp1, kp2
                    part = (iw(kp)-1)/n
                    ke = iw(kp) - n*part
                    IF (flag(ke)<=-2) THEN
! KE is an element.  Flag it.
                      flag(ke) = -nflg
                      IF ( .NOT. bound) THEN
! Store the part involved in preparation for the test
! for identical
! variables.
                        leaf(ke) = part
! Find the bounds for the search of the element list
                        jp = ipe(ke)
                        jp1 = jp + 3
                        jp2 = jp + iw(jp)
                        IF (part/=0) THEN
! Supervariable in full part of element
                          IF (part==2) THEN
! Supervariable in leading part of element
                            jp1 = jp1 + iw(jp+1)
                          ELSE

! Supervariable in trailing part of element
                            jp2 = jp2 - iw(jp+2)
                          END IF
                        END IF
! Search list of element KE, revising the row count
! when
! new
! variables found.
                        irl = ir
                        DO jp = jp1, jp2
                          js = iw(jp)
! Cycle if JS has been eliminated or already
! counted.
                          IF (flag(js)>nflg) THEN
                            ir = ir + abs(nv(js))
                            flag(js) = nflg
                          END IF
                        END DO
! KP need not be stored in the new list if it is
! linked
! to ME and KE
! made no contribution to the count.
                        IF (ir==irl .AND. last(ke)==-me) GO TO 200
                      END IF
! Store KP in the new list
                      iw(np) = iw(kp)
                      np = np + 1

                    ELSE IF (flag(ke)>=0) THEN
                      GO TO 210
                    END IF
200                 CONTINUE
                  END DO

                  np0 = np
                  GO TO 220
! KE is a supervariable
210               CONTINUE
! Treat the rest of the list associated with supervariable
! IS. It consists entirely of supervariables.
                  np0 = np
                  kp1 = kp
                  DO kp = kp1, kp2
                    ks = iw(kp)
                    IF (flag(ks)>nflg) THEN
! Add to row count, flag supervariable KS and add it to
! new list.
                      ir = ir + abs(nv(ks))
                      flag(ks) = nflg
                      iw(np) = ks
                      np = np + 1
                    END IF
                  END DO
220               IF (bound) ir = count(is)
                  IF (np>kp2) THEN
! List is longer than it was. Copy it to free space.
                    kp = ipe(is)
                    IF (np+iwfr-kp>=lw) THEN
! Compress IW
                      CALL mc68_ma47_compress(n,ipe,flag,iw,iwfr-1,iwfr, &
                        info%n_compressions)
                      IF (np+iwfr-kp>=lw) THEN
! Expand the length of iw
                        szw = lw
                        lw = np + iwfr - kp + 1
                        CALL zb01_resize1(iw,szw,lw,info_zb01)
                        IF (info_zb01%flag<0) THEN
                          IF (info_zb01%flag==-11) THEN
                            info%flag = mc68_err_memory_alloc
                            info%stat = info_zb01%stat
                          ELSE IF (info_zb01%flag==-12) THEN
                            info%flag = mc68_err_memory_dealloc
                            info%stat = info_zb01%stat
                          ELSE IF (info_zb01%flag<-12) THEN
                            info%flag = mc68_err_zb01
                            info%iostat = info_zb01%iostat
                            info%zb01_info = info_zb01%flag
                          END IF
                          RETURN
                        END IF
                      END IF

                      shift = ipe(is) - kp
                      kp = kp + shift
                      kp2 = kp2 + shift
                      np = np + shift
                      np0 = np0 + shift
                    END IF

                    np = np + iwfr - kp
                    np0 = np0 + iwfr - kp
                    ipe(is) = iwfr
                    kp1 = kp
                    DO kp = kp1, kp2
                      iw(iwfr) = iw(kp)
                      iwfr = iwfr + 1
                    END DO
                    iw(iwfr) = 0
                    iwfr = iwfr + 1
                  END IF
! Move first supervariable to end of list, move first
! element
! to end
! of element part of list and add new element to front of
! list.
                  iw(np) = iw(np0)
                  kp = ipe(is)
                  iw(np0) = iw(kp+1)
                  iw(kp+1) = me + (loop-1)*n
! Store the new length of the list.
                  iw(kp) = np - kp
                  IF (ir==0) THEN
! Treat zero row as if it had a nonzero diagonal block
                    ir = -nv(is)
                    nv(is) = ir
                    rank = rank - ir
                    ir = 1
                  END IF

! Unless the Markowitz cost is zero (which means that rows
! not
! previously identical cannot be identical now) or the row
! count is
! above the threshold, check whether row IS is identical to
! another by
! looking in the linked list of supervariables with row
! count
! IR at
! those whose lists have first entry ME and occur in the
! same
! part
! of it. Note that those containing ME come first so the
! search can
! be terminated when a list not starting with ME is found.
                  IF (cmin/=zero) THEN
                    IF (ir<=thresh) THEN
                      js = ipr(ir)
                      DO list = 1, n
                        IF (js<=0) THEN
                          GO TO 260
                        ELSE
                          kp = ipe(js)
                          IF (iw(kp+1)/=me+(loop-1)*n) THEN
                            GO TO 260
                          ELSE
! JS has same row count and is active. Check if
! identical to IS.
                            IF (sign(1,nv(js))==sign(1,nv(is))) THEN
                              kp1 = kp
                              DO kp = kp1 + 2, kp1 + iw(kp1)
                                part = (iw(kp)-1)/n
                                ie = iw(kp) - part*n
! Jump if IE (which may be a
! supervariable or an element)
! is not in the list of IS.
                                IF (abs(flag(ie))>nflg) THEN
                                  GO TO 230
! Jump if JS belongs to another part
! of IE
                                ELSE IF (flag(ie)==-nflg) THEN
                                  IF (part/=leaf(ie)) GO TO 230
                                END IF
                              END DO
                              GO TO 240
                            END IF

230                         js = next(js)
                          END IF
                        END IF
                      END DO
                      GO TO 250

240                   CONTINUE
! Supervariable amalgamation. Row IS is identical to row
! JS.  Regard all variables in the two supervariables as
! being in IS. Set tree pointer, FLAG and NV entries.
250                   ipe(js) = 0
                      nv(is) = nv(is) + nv(js)
                      nv(js) = 0
                      flag(js) = -1
! Replace JS by IS in linked list.
                      ns = next(js)
                      ls = last(js)
                      IF (ns>0) last(ns) = is
                      IF (ls>0) next(ls) = is
                      last(is) = ls
                      next(is) = ns
                      IF (ipr(ir)==js) ipr(ir) = is
                      count(is) = ir
                      next(js) = -leaf(is)
                      leaf(is) = leaf(js)
                      last(js) = -is
                      GO TO 270
                    END IF
                  END IF
! Insert IS into linked list of supervariables of same row
! count.
260               ns = ipr(ir)
                  IF (ns>0) last(ns) = is
                  next(is) = ns
                  ipr(ir) = is
                  last(is) = 0
                  minr = min(minr,ir)
                  IF (nv(is)>0) minf = min(minf,ir)
                  count(is) = ir
270               CONTINUE
                END DO
              END DO

! Reset flags for supervariables in newly created element,
! remove those absorbed into others, and store the new list.
              IF (iwfr+nsvars+3>=lw) THEN
! Compress IW
                CALL mc68_ma47_compress(n,ipe,flag,iw,iwfr-1,iwfr, &
                  info%n_compressions)
                IF (iwfr+nsvars+3>=lw) THEN
                  szw = lw
                  lw = iwfr + nsvars + 4
                  CALL zb01_resize1(iw,szw,lw,info_zb01)
                  IF (info_zb01%flag<0) THEN
                    IF (info_zb01%flag==-11) THEN
                      info%flag = mc68_err_memory_alloc
                      info%stat = info_zb01%stat
                    ELSE IF (info_zb01%flag==-12) THEN
                      info%flag = mc68_err_memory_dealloc
                      info%stat = info_zb01%stat
                    ELSE IF (info_zb01%flag<-12) THEN
                      info%flag = mc68_err_zb01
                      info%iostat = info_zb01%iostat
                      info%zb01_info = info_zb01%flag
                    END IF
                    RETURN
                  END IF
! GO TO 280
                ELSE
                END IF
              END IF

              ip = iwfr
              iwfr = iwfr + 3
              k2 = 0
              DO loop = 1, 3
                k1 = k2 + 1
                k2 = k1 + nsc(loop) - 1
                DO k = k1, k2
                  is = svars(k)
                  IF (nv(is)==0) THEN
                    nsc(loop) = nsc(loop) - 1
                  ELSE
                    flag(is) = nflg
                    iw(iwfr) = is
                    iwfr = iwfr + 1
                  END IF
                END DO
              END DO
              iw(ip) = iwfr - ip - 1
              iw(ip+1) = nsc(1)
              iw(ip+2) = nsc(3)
              flag(me) = -nflg
! Set pointer for new element.
              ipe(me) = ip
            END IF
          END IF
        END DO

280     IF (rank<n) info%flag = info%flag + mc68_warn_rank
        info%n_zero_eigs = n - rank


        DEALLOCATE (flag,ipr,last,leaf,svars,STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_dealloc
          RETURN
        END IF

END SUBROUTINE mc68_ma47_analyse

SUBROUTINE mc68_ma47_treesearch(n,father,count,ne,perm,nemin,info)
! Given the tree provided by MC68_MA47_ANALYSE/KD, which has a node
! for every
! variable, perform depth-first search to find pivot order and
! construct a new tree whose nodes correspond to block eliminations,
! and whose depth-first-search post ordering is 1,2,...,NODES.

! N must be set to the matrix order. It is not altered.
! FATHER(I) must be set equal to -(father of node I) in the given tree
! or zero if node I is a root. It is altered to point to its next
! younger brother if it has one. Finally, it is set to the father
! of node I in the new tree, or NODES+1 if node I is a root.
! COUNT(I) must be set to the row count at the time of elimination if
! I
! is a node at which an elimination takes place, negated for the
! second node of a tile or oxo pivot. It is revised for a node at
! which an amalgamaton occurs, but is otherwise unchanged.
! NE(I) must be set to the number of variables eliminated at node I,
! negated for a defective variable. It is revised for a node at
! which an amalgamaton occurs, but is otherwise unchanged.
! PERM need not be set.

! Local variables
! SON(I) need not be set. It is used to hold
! (eldest son of node I) if it has one and 0 otherwise.
! NODE(I) need not be set. It is used temporarily to hold the row
! count of the secondary variable of a tile or oxo , -1 for a
! node merged into its father and 0 otherwise. It is finally
! set to hold the elimination stage of variable I, negated for a
! defective variable.
! NA need not be set. NA(I) is used temporarily to hold the next
! younger
! brother of node I if it has one or -(father) if it does not.
! NA(STAGE) and NA(LEVEL) are used to hold the numbers of elements
! assembled at stage STAGE and level LEVEL of the elimination.
! MARK need not be set. The indices of the root nodes are stored in
! MARK(I), I=NR,N. MARK(STAGE) eventually is set to the Markowitz
! cost at stage STAGE of the elimination.
! NODES need not be set. It is set to the number of elimination steps.
! ICNTL must be set by the user as follows and is not altered. The
! only entry used is
! ICNTL(6) Two nodes of the assembly tree are merged only if both
! involve less than ICNTL(6) eliminations.
! COUNT2 Row count of the second row of a tile or oxo.
! I   Index of current tree node.
! IBRTHR  Brother of node I.
! IFATHR  Father of node I.
! ISON Son of node I.
! J Temporary variable.
! K Position in elimination sequence
! L Temporary variable.
! LDIAG Control for amount of information output.
! LEVEL holds the current tree level. Roots are at level N, their sons
! are at level N-1, etc.
! MP Stream number for warnings and diagnostics.
! NDE Do index when looping through new tree.
! NEMIN is used to control the amalgamation process between
! a son and its father (if the number of fully summed
! variables of both nodes is smaller than NEMIN).
! NST Number of structured pivots in the active chain of the
! depth-first search.
! NR Position of next root.
! NRL The indices of the root nodes are stored in MARK(I), I=NRL,N.
! NVPIV holds the number of variables in the pivot (each half in the
! case of an oxo).
! PERM1(LEVEL) is used temporarily to hold the node at level LEVEL of
! the
! active chain. PERM1(K) is eventually set to the variable that is in
! position K of the elimination.
! STAGE holds the current elimination stage.
! Set local print variables
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: n, nemin
! ..
! .. Array Arguments ..
        INTEGER, INTENT (INOUT) :: count(n), father(n), ne(n)
        INTEGER, INTENT (OUT) :: perm(n)
        TYPE (mc68_info), INTENT (INOUT) :: info
! ..
! .. Local Scalars ..
        INTEGER :: count2, i, ibrthr, ifathr, ison, j, k, l, ldiag, level, mp, &
          nr, nrl, nst, nvpiv, stage
        INTEGER, ALLOCATABLE, DIMENSION (:) :: mark, na, node, son, perm1
! ..
        mp = -1
        ldiag = 0
        IF (mp<=0) ldiag = 0


        ALLOCATE (mark(n),na(n),node(n),son(n),perm1(n),STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_alloc
          RETURN
        END IF


! Set SON(I) to (eldest son of node I) and NA(I) to next younger
! brother of node I if it has one or -(father) if it does not.
! Find the root nodes.
! Load all the eliminations of any oxo and tile pivot on its first
! node and store the second count in NODE.
        DO i = 1, n
          son(i) = 0
          node(i) = 0
        END DO
        nrl = n + 1
        DO i = 1, n
          ifathr = -father(i)
          na(i) = -ifathr
          IF (ifathr==0) THEN
! We have a root
            nrl = nrl - 1
            mark(nrl) = i
          ELSE
            ibrthr = son(ifathr)
! If node I has a next younger brother, set pointer to him.
            IF (ibrthr>0) na(i) = ibrthr
! Set father-to-eldest-son pointer
            son(ifathr) = i
            IF (count(i)<0) THEN
! Second node of tile or oxo. Load the eliminations on the first
! node
! and store the second count.
              ne(ifathr) = ne(ifathr)*2
              ne(i) = 0
              node(ifathr) = -count(i)
            END IF
          END IF
        END DO

! Depth-first search looking for father-son pairs that can be merged.
! Adjust NE, COUNT, and NODE to correspond.
! PERM1(LEVEL) is used temporarily to hold the node at level
! LEVEL of the active chain.
        nr = nrl
        i = 0
        DO k = 1, n
          IF (i<=0) THEN
! Pick up next root.
            i = mark(nr)
            nr = nr + 1
            level = n
            nst = 0
            perm1(level) = i
          END IF
! Go to son for as long as possible, clearing father-son pointers
! in SON as each is used.
          DO l = 1, n
            IF (son(i)<=0) THEN
              GO TO 10
            ELSE
              ison = son(i)
              son(i) = 0
              i = ison
              IF (node(i)/=0) nst = nst + 1
              level = level - 1
              perm1(level) = i
            END IF
          END DO
! Jump if node has no eliminations
10        nvpiv = ne(i)
          IF (nvpiv/=0) THEN
! Jump if at a root.
            IF (level/=n) THEN
              ifathr = perm1(level+1)
              IF (node(i)/=0 .OR. node(ifathr)/=0) THEN
                IF (node(i)==0 .OR. node(ifathr)==0) THEN
                  GO TO 20
                ELSE
                  IF (nvpiv>0 .AND. ne(ifathr)>0) THEN
! Both nodes are tiles.
! Do not merge if this would cause fill.
                    IF (node(i)-nvpiv/2/=node(ifathr)) THEN
                      GO TO 20
                    ELSE IF (count(i)-nvpiv/=count(ifathr)) THEN
                      GO TO 20
                    END IF
                  ELSE IF (nvpiv<0 .AND. ne(ifathr)<0) THEN
! Both nodes are oxos.
                    nvpiv = -nvpiv/2
! Do not merge unless sure that this would not cause fill.
                    IF (node(i)-nvpiv/=node(ifathr)) THEN
                      GO TO 20
                    ELSE IF (count(i)-nvpiv/=count(ifathr)) THEN
                      GO TO 20
                    ELSE IF (count(i)==node(i)) THEN
                      GO TO 20
                    END IF
                  ELSE
                    GO TO 20
                  END IF
! Merge two tiles or two oxos.
                  node(ifathr) = node(i)
                  nst = nst - 1
                END IF
! Both nodes are full
! Merge the nodes if this would cause no more fill
              ELSE IF (count(i)-nvpiv/=count(ifathr)) THEN
! Do not merge the nodes if there is a structured pivot on
! path
! to root
                IF (nst>0) THEN
                  GO TO 20
! Do not merge the nodes if either node has more than NEMIN
! eliminations
                ELSE IF (nvpiv>=nemin) THEN
                  GO TO 20
                ELSE IF (ne(ifathr)>=nemin) THEN
                  GO TO 20
                END IF
              END IF
! Merge the two nodes
              ne(ifathr) = ne(ifathr) + ne(i)
              ne(i) = 0
              IF (ldiag>4) WRITE (mp,fmt='(A,2I5)') ' Merging nodes', i, &
                ifathr
              count(ifathr) = count(ifathr) + nvpiv
              node(i) = -1
            END IF
          END IF
! Go to next brother, or failing this to the father.
20        ibrthr = na(i)
          IF (node(i)>0) nst = nst - 1
          IF (ibrthr>0) THEN
! Node I has a younger brother. Go to him.
            perm1(level) = ibrthr
            i = ibrthr
            IF (node(i)>0) nst = nst + 1
          ELSE
! Go to father of node I
            level = level + 1
            i = -ibrthr
          END IF
        END DO
        DO i = 1, n
          son(i) = 0
        END DO

! Set SON(I) to (eldest son of node I) and NA(I) to next younger
! brother of node I if it has one or -(father) if it does not.
        DO i = 1, n
          ifathr = -father(i)
          na(i) = -ifathr
          IF (ifathr/=0) THEN
            ibrthr = son(ifathr)
! If node I has a next younger brother, set pointer to him.
            IF (ibrthr>0) na(i) = ibrthr
! Set father-to-son pointer
            son(ifathr) = i
          END IF
        END DO

! Depth-first search in which FATHER is revised to accord with the
! merges. PERM1(LEVEL) is used temporarily to hold the node at level
! LEVEL of the active chain.
        i = 0
        nr = nrl
        DO k = 1, n
          IF (i<=0) THEN
! Pick up next root.
            i = mark(nr)
            nr = nr + 1
            level = n
            perm1(n) = i
          END IF
! Go to son for as long as possible, clearing father-son pointers
! in SON as each is used.
          DO l = 1, n
            IF (son(i)<=0) THEN
              GO TO 30
            ELSE
              ison = son(i)
              son(i) = 0
              i = ison
              father(i) = -perm1(level)
              IF (node(i)>=0) THEN
                level = level - 1
                perm1(level) = i
              END IF
            END IF
          END DO

! Go to next brother, or failing this to the father.
30        ibrthr = na(i)
          IF (ibrthr>0) THEN
! Node I has a younger brother. Go to him.
            IF (node(i)<0) level = level - 1
            i = ibrthr
            perm1(level) = i
            father(i) = -perm1(level+1)
            IF (node(i)<0) level = level + 1
          ELSE
! Go to father of node I
            IF (node(i)>=0) level = level + 1
            i = -ibrthr
          END IF
        END DO

! Set SON(I) to (eldest son of node I) and FATHER(I) to next younger
! brother of node I if it has one.
        DO i = 1, n
          son(i) = 0
        END DO
! First pass is for nodes without eliminations.
        DO i = 1, n
          IF (ne(i)==0 .AND. count(i)>=0) THEN
            ifathr = -father(i)
            ibrthr = son(ifathr)
! If node I has a next younger brother, set pointer to him.
            IF (ibrthr>0) father(i) = ibrthr
! Set father-to-eldest-son pointer
            son(ifathr) = i
          END IF
        END DO

! Second pass is for second nodes of tile and oxo pivots.
        DO i = 1, n
          IF (count(i)<0) THEN
            ifathr = -father(i)
            ibrthr = son(ifathr)
! If node I has a next younger brother, set pointer to him.
            IF (ibrthr>0) father(i) = ibrthr
! Set father-to-eldest-son pointer
            son(ifathr) = i
          END IF
        END DO

! Third pass is for elimination nodes.
        DO i = 1, n
          IF (ne(i)/=0) THEN
            ifathr = -father(i)
            IF (ifathr/=0) THEN
              ibrthr = son(ifathr)
! If node I has a next younger brother, set pointer to him.
              IF (ibrthr>0) father(i) = ibrthr
! Set father-to-son pointer
              son(ifathr) = i
            END IF
          END IF
        END DO

! Depth-first search.
! STAGE holds the current elimination stage. The number
! of assemblies for nodes of the active chain is accumulated
! temporarily in NA(LEVEL), for tree level N, N-1,..., and is
! transfered to NA(STAGE) when we reach the appropriate stage STAGE.
        stage = 1
! I is the current node.
        i = 0
        nr = nrl
        DO k = 1, n
          IF (i<=0) THEN
! Pick up next root.
            i = mark(nr)
            nr = nr + 1
            level = n
            na(n) = 0
          END IF
! Go to son for as long as possible, clearing father-son pointers
! in SON as each is used and setting NA(LEVEL)=0 for all levels
! reached.
          l = level
          DO level = l, 1, -1
            IF (son(i)<=0) THEN
              GO TO 40
            ELSE
              ison = son(i)
              son(i) = 0
              i = ison
              na(level-1) = 0
            END IF
          END DO
! Record variable in position I in the order.
40        perm1(k) = i
          count2 = node(i)
          node(i) = stage
! Jump if node has no eliminations
          IF (ne(i)/=0) THEN
            IF (level<n) na(level+1) = na(level+1) + 1
            na(stage) = na(level)
            IF (count2==0) THEN
! Full pivot
              mark(stage) = (count(i)-1)**2
            ELSE IF (ne(i)>0) THEN
! Tile pivot
              DO j = k - ne(i) + 1, k - ne(i)/2
                node(perm1(j)) = -node(perm1(j))
              END DO
              mark(stage) = (count(i)+count2-3)*(count2-1)
            ELSE
! Oxo pivot
              DO j = k + ne(i) + 1, k
                node(perm1(j)) = -node(perm1(j))
              END DO
              mark(stage) = (count(i)-1)*(count2-1)
            END IF
            stage = stage + 1
          END IF
! Go to next brother, or failing this to the father.
          ibrthr = father(i)
          IF (ibrthr>0) THEN
! Node I has a younger brother. Go to him.
            na(level) = 0
            i = ibrthr
          ELSE
! Go to father of node I
            level = level + 1
            ison = i
            i = -ibrthr
          END IF
        END DO

! Search for 2x2 pivots
        i = 1
        DO WHILE (i<n)
          IF (node(perm1(i))>0) THEN
            i = i + 1
          ELSE
            perm1(i) = -perm1(i)
            perm1(i+1) = -perm1(i+1)
            i = i + 2
          END IF
        END DO
        DO i = 1, n
          IF (perm1(i)>0) THEN
            perm(perm1(i)) = i
          ELSE
            perm(-perm1(i)) = -i
          END IF
        END DO

        DEALLOCATE (mark,na,node,son,perm1,STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_dealloc
          RETURN
        END IF


END SUBROUTINE mc68_ma47_treesearch

SUBROUTINE mc68_ma47_merge(n,ipe,iw,lw,nv,next,last,leaf,flag,var,svar)

! Merge variables whose rows are identical into supervariables. This
! is done by updating the supervariable structure for the submatrix
! of the first J-1 columns to that for the first J columns,
! J=1,2,...,N. If a nontrivial supervariable IS is involved in
! column J, all of its variables that are involved in column J are
! removed to make a new supervariable.

! N must be set to the matrix order. It is not altered.
! IPE(I) must be set to the position in IW of the list for row I,
! I=1,N. On return, it is unchanged except that if I is a variable
! that has been absorbed into a supervariable, IPE(I) = 0.
! IW must be set on entry to hold lists of entries by rows, each list
! being headed by its length. There must be no duplicate indices
! in a list. IW is not altered.
! LW must be set to the length of IW. It is not altered.
! NV must be set on entry so that NV(I) = 1 for a nondefective
! variable
! (nonzero diagonal entry) and NV(I) = -1 otherwise, I = 1,..., N.
! On return:
! if I is a variable that has been absorbed, NV(I) = 0;
! if I is a supervariable, NV(I) holds its number of variables,
! negated for a defective supervariable.
! NEXT need not be set. During execution, NEXT(I) is the next
! variable in a circular list of variables in a supervariable.
! On return, if variable I has been absorbed and J is the next
! variable in its supervariable, NEXT(I)=-J;
! LAST need not be set. During execution, LAST(I) is the previous
! variable in a circular list of variables in a supervariable.
! On return, if variable I has been absorbed and J is the next
! variable in its supervariable, LAST(I)=-J;
! LEAF need not be set. On return, if IS is a supervariable, LEAF(IS)
! holds its first variable.
! FLAG is used as workspace for supervariable flags. If IS has already
! been encountered in column J, FLAG(IS) = J.
! VAR is used as workspace. If IS is a supervariable involved in
! column
! J, VAR(IS) is the first variable to be removed from supervariable
! IS. If IS is not a supervariable, VAR(IS) is the next free
! supervariable index in a chain of such indices.
! SVAR is used as workspace. SVAR(I) is the supervariable to which I
! belongs.

! Local variables

! FREE The first free supervariable index (head of chain)
! I    Row index.
! IS   Supervariable.
! JS   Supervariable.
! J    Column index.
! K    DO index.
! KK   Temporary variable.
! LS   Last (previous) variable in circular list.
! NS   Next variable in circular list.

! Begin by setting SVAR, LAST, and NEXT to represent all variables
! belonging to supervariable 1. Also initialize FLAG and set
! FREE and VAR as a chain of indices that are free for use as
! supervariable names.
! .. Scalar Arguments ..
        INTEGER (long), INTENT (IN) :: lw
        INTEGER, INTENT (IN) :: n
! ..
! .. Array Arguments ..
        INTEGER, INTENT (IN) :: iw(lw)
        INTEGER, INTENT (INOUT) :: ipe(n), nv(n)
        INTEGER, INTENT (OUT) :: flag(n), last(n), leaf(n), next(n), svar(n), &
          var(n)
! ..
! .. Local Scalars ..
        INTEGER :: free, i, is, j, js, k, kk, ls, ns
! ..
! .. Intrinsic Functions ..
        INTRINSIC sign
! ..
        DO i = 1, n
          svar(i) = 1
          last(i) = i - 1
          next(i) = i + 1
          flag(i) = 0
          var(i) = i + 1
        END DO
        last(1) = n
        next(n) = 1
        free = 2

! Scan the columns in turn, splitting the supervariables that are
! involved.
        DO j = 1, n
          kk = ipe(j)
          DO k = kk + 1, kk + iw(kk)
            i = iw(k)
            is = svar(i)
            IF (flag(is)/=j) THEN
! First occurrence of supervariable IS for column J
              flag(is) = j
              IF (next(i)/=i) THEN
! No action needed since IS has I as its only variable.
! Establish new supervariable
                js = free
                free = var(js)
                var(is) = i
                svar(i) = js
! Remove I from old circular list
                ns = next(i)
                ls = last(i)
                next(ls) = ns
                last(ns) = ls
! Make new circular list
                next(i) = i
                last(i) = i
              END IF
            ELSE
! Subsequent occurrence of IS for column J
! Remove I from old circular list
              IF (next(i)==i) THEN
! Supervariable now empty
                ns = var(is)
                var(is) = free
                free = is
              ELSE
                ns = next(i)
                ls = last(i)
                next(ls) = ns
                last(ns) = ls
                ns = var(is)
              END IF
! Add I to new list
              ls = last(ns)
              next(ls) = i
              next(i) = ns
              last(ns) = i
              last(i) = ls
              svar(i) = svar(ns)
            END IF
          END DO
        END DO

! Set the data in final format.
        DO is = 1, n
          leaf(is) = is
! No action for a trivial supervariable
          IF (last(is)/=is) THEN
! No action if already treated
            IF (last(is)>=0) THEN
! IS is a nontrivial supervariable
              ls = last(is)
              leaf(is) = ls
              DO k = 1, n
                i = ls
                IF (i==is) THEN
                  GO TO 10
                ELSE
                  ipe(i) = 0
                  ls = last(i)
                  next(i) = -ls
                  last(i) = -ls
                  nv(i) = 0
                END IF
              END DO
10            nv(is) = sign(k,nv(is))
            END IF
          END IF
        END DO
END SUBROUTINE mc68_ma47_merge

SUBROUTINE mc68_ma47_reset(n,flag,nflg)
! Reset flag values to +/- N*3.
! N Matrix order. Unchanged.
! FLAG Unchaged if FLAG(I) = -1, 0, 1, or 2.
! Otherwise, reset to +/- N*3.

! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: n
        INTEGER, INTENT (OUT) :: nflg
! ..
! .. Array Arguments ..
        INTEGER, INTENT (INOUT) :: flag(n)
! ..
! .. Local Scalars ..
        INTEGER :: i, n3
! ..
        n3 = n*3
        DO i = 1, n
          IF (flag(i)>2) flag(i) = n3
          IF (flag(i)<=-2) flag(i) = -n3
        END DO
        nflg = n3
END SUBROUTINE mc68_ma47_reset

SUBROUTINE mc68_ma47_thresh(thresh,newthr,n,ipe,iw,lw,count,nv,next, &
          last,ipr,flag,nflg)
! Alter the threshold for calculation of row counts.

! THRESH Old threshold. Changed to NEWTHR.
! NEWTHR New threshold. It is not altered.
! N must be set to the matrix order. It is not altered.
! IPE(I) must be set to the position in IW of the list for row I,
! I=1,N.
! It is not altered.
! IW must be set on entry to hold lists of entries by rows, each list
! being headed by its length. It is not altered.
! LW must be set to the length of IW. It is not altered.
! COUNT is used to hold the row counts in the reduced matrix.
! COUNT(IS) is recalculated if THRESH < COUNT(IS) <= NEWTHR.
! NV: if JS is a supervariable, NV(JS) holds its number of variables,
! negated for a defective supervariable. It is not altered.
! NEXT: if JS is a supervariable that has not been eliminated or
! absorbed, NEXT(JS) is the next such supervariable having the
! same row count, or zero if JS is last in its list. Revised
! if a recalculated row count makes it necessary.
! LAST: if JS is a supervariable that has not been eliminated or
! absorbed, LAST(JS) is the previous such supervariable having
! the same row count, or zero if JS is last in its list.
! Revised if a recalculated row count makes it necessary.
! IPR: IPR(IR) is the first supervariable with row count IR or zero
! if there are none. Revised if a recalculated row count makes it
! necessary.
! FLAG:  if IS is a supervariable, FLAG(IS) >= 0;
! if IE is an element, FLAG(IE) < -1;
! if I is neither a supervariable nor an element, FLAG(I) = -1.
! The values of the supervariable flags may change while remaining
! positive. Otherwise, unchanged.
! NFLG is used for the current flag value in array FLAG, as in
! MC68_MA47_ANALYSE.
! Local variables

! .. Scalar Arguments ..
        INTEGER (long), INTENT (IN) :: lw
        INTEGER, INTENT (IN) :: n, newthr
        INTEGER, INTENT (INOUT) :: nflg, thresh
! ..
! .. Array Arguments ..
        INTEGER, INTENT (IN) :: ipe(n), iw(lw), nv(n)
        INTEGER, INTENT (INOUT) :: count(n), flag(n), ipr(n), last(n), next(n)
! ..
! .. Local Scalars ..
        INTEGER :: ir, is, jp, jp1, jp2, k, ke, kp, kp2, ls, ms, ns, part
! ..
! ..
! .. Intrinsic Functions ..
        INTRINSIC abs
! ..
        DO ms = 1, n
          IF (flag(ms)>=0) THEN
            IF (count(ms)>thresh) THEN
              IF (count(ms)<=newthr) THEN
! Find the row count afresh
                ir = 0
                IF (nflg<=4) CALL mc68_ma47_reset(n,flag,nflg)
! Reduce NFLG by one to cater for this supervariable.
                nflg = nflg - 1
                k = ipe(ms)
                kp2 = k + iw(k)
                DO kp = k + 1, kp2
                  part = (iw(kp)-1)/n
                  ke = iw(kp) - part*n
                  IF (flag(ke)/=-1) THEN
                    IF (flag(ke)<=-2) THEN
! KE is an element.
                      jp = ipe(ke)
                      jp1 = jp + 3
                      jp2 = jp + iw(jp)
                      IF (part/=0) THEN
! Supervariable in full part of element
                        IF (part==2) THEN
! Supervariable in leading part of element
                          jp1 = jp1 + iw(jp+1)
                        ELSE
! Supervariable in trailing part of element
                          jp2 = jp2 - iw(jp+2)
                        END IF
                      END IF
                    ELSE
! We have reached the list of variables
                      jp1 = kp
                      jp2 = kp2
                    END IF
! Search the variable list.
                    DO jp = jp1, jp2
                      is = iw(jp)
                      IF (flag(is)>nflg) THEN
                        flag(is) = nflg
                        ir = ir + abs(nv(is))
                      END IF
                    END DO
                    IF (jp2==kp2 .OR. ir>newthr) GO TO 10
                  END IF
                END DO
10              IF (ir/=count(ms)) THEN
! Remove MS from linked list
                  ns = next(ms)
                  ls = last(ms)
                  next(ms) = 0
                  IF (ns>0) last(ns) = ls
                  IF (ls>0) THEN
                    next(ls) = ns
                  ELSE
                    ipr(count(ms)) = ns
                  END IF
! Insert MS into linked list of supervariables of same row
! count.
                  ns = ipr(ir)
                  IF (ns>0) last(ns) = ms
                  next(ms) = ns
                  ipr(ir) = ms
                  last(ms) = 0
                  count(ms) = ir
                END IF
              END IF
            END IF
          END IF
        END DO
        thresh = newthr
END SUBROUTINE mc68_ma47_thresh

SUBROUTINE mc68_ma47_compress(n,ipe,flag,iw,lw,iwfr,ncmpa)
! Compress lists held by MC68_MA47_ANALYSE and MA47KD in IW, removing
! inactive
! variables from element lists. Adjust pointers in IPE to correspond.
! N is the matrix order. It is not altered.
! IPE(I) points to the position in IW of the start of list I or is
! zero if there is no list I. On exit it points to the new position.
! FLAG holds element and supervariable flags:
! if IS is an active supervariable, FLAG(IS) >= 0;
! if IS is inactive, FLAG(IS) = -1;
! if IE is an element, FLAG(IE) < -1;
! It is not altered.
! IW holds the lists, each headed by its length. Every entry must be
! nonnegative. On output, the same lists are held, but they are
! now compressed together.
! LW holds the length of IW. It is not altered.
! IWFR need not be set on entry. On exit it points to the first free
! location in IW.

! Local variables
! I    Row index.
! IR   DO index.
! K    Temporary variable.
! L    Temporary variable.
! LEN1 Length of the leading zero part of element list
! LEN2 Length of the middle part of element list
! LEN3 Length of the trailing zero part of element list
! LWFR points just beyond the end of the uncompressed file.

! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: lw, n
        INTEGER, INTENT (INOUT) :: ncmpa
        INTEGER, INTENT (OUT) :: iwfr
! ..
! .. Array Arguments ..
        INTEGER, INTENT (IN) :: flag(n)
        INTEGER, INTENT (INOUT) :: ipe(n), iw(lw)
! ..
! .. Local Scalars ..
        INTEGER :: i, ir, k, l, len1, len2, len3, lwfr
! ..
        ncmpa = ncmpa + 1
! Prepare for compressing by storing the lengths of the lists in IPE
! and
! setting the first entry of each list to -(list number).
        DO i = 1, n
          l = ipe(i)
          IF (l>0 .AND. flag(i)/=-1) THEN
            ipe(i) = iw(l)
            iw(l) = -i
          END IF
        END DO

! Compress
! IWFR points just beyond the end of the compressed file.
        iwfr = 1
        lwfr = 1
        DO ir = 1, n
! Search for the next negative entry.
          DO k = lwfr, lw
            IF (iw(k)<0) GO TO 10
          END DO
          RETURN
! Pick up entry number, store length in new position, set new
! pointer
! and prepare to copy list.
10        i = -iw(k)
          iw(iwfr) = ipe(i)
          ipe(i) = iwfr
          iwfr = iwfr + 1
          IF (flag(i)<=-2) THEN
! We have an element list. Remove dummy entries.
            l = iwfr - 1
            len1 = iw(k+1)
            iw(l+1) = len1
            len3 = iw(k+2)
            iw(l+2) = len3
            len2 = iw(l) - len1 - len3 - 2
            iwfr = l + 3
            DO lwfr = k + 3, k + 2 + len1
              IF (flag(iw(lwfr))<0) THEN
                iw(l+1) = iw(l+1) - 1
                iw(l) = iw(l) - 1
              ELSE
                iw(iwfr) = iw(lwfr)
                iwfr = iwfr + 1
              END IF
            END DO

            k = lwfr
            DO lwfr = k, k - 1 + len2
              IF (flag(iw(lwfr))<0) THEN
                iw(l) = iw(l) - 1
              ELSE
                iw(iwfr) = iw(lwfr)
                iwfr = iwfr + 1
              END IF
            END DO

            k = lwfr
            DO lwfr = k, k - 1 + len3
              IF (flag(iw(lwfr))<0) THEN
                iw(l+2) = iw(l+2) - 1
                iw(l) = iw(l) - 1
              ELSE
                iw(iwfr) = iw(lwfr)
                iwfr = iwfr + 1
              END IF
            END DO
          ELSE
! Copy list to new position.
            DO lwfr = k + 1, k + iw(iwfr-1)
              iw(iwfr) = iw(lwfr)
              iwfr = iwfr + 1
            END DO
          END IF
        END DO

END SUBROUTINE mc68_ma47_compress

SUBROUTINE amdd(n,iwlen,pe,pfree,len,iw,elen,icntl,info)

! -------------------------------------------------------------------
! AMDD is a modified version of MC47 version 1. Dense rows are detected
! and removed. AMDD is then applied to the remaining matrix. The dense rows are
! then appended to the end of this ordering.

! We use the term Le to denote the set of all supervariables in element
! E.
! A row is declared as dense if removing it from the current matrix
! will result in a significant reduction in the mean degree of the remaining
! rows
! A row is sparse if it is not dense.
! -------------------------------------------------------------------

! .. Parameters ..
        INTEGER, PARAMETER :: wp = kind(0.0D0)
        REAL (wp) :: zero
        PARAMETER (zero=0E0_wp)
! ..
! .. Scalar Arguments ..
        INTEGER (long), INTENT (IN) :: iwlen
        INTEGER, INTENT (IN) :: n
        INTEGER, INTENT (INOUT) :: pfree
! ..
! .. Array Arguments ..
        INTEGER, INTENT (IN) :: icntl(10)
        INTEGER, INTENT (INOUT) :: iw(iwlen), len(n), pe(n)
        INTEGER, INTENT (OUT) :: elen(n)

        TYPE (mc68_info), INTENT (INOUT) :: info

! N must be set to the matrix order.
! Restriction:  N .ge. 1

! IWLEN must be set to the length of IW. It is not altered. On input,
! the matrix is stored in IW (1..PFREE-1).
! *** We do not recommend running this algorithm with ***
! ***      IWLEN .LT. PFREE + N.                      ***
! *** Better performance will be obtained if          ***
! ***      IWLEN .GE. PFREE + N                       ***
! *** or better yet                                   ***
! ***      IWLEN .GT. 1.2 * PFREE                     ***
! Restriction: IWLEN .GE. PFREE-1

! PE(i) must be set to the the index in IW of the start of row I, or be
! zero if row I has no off-diagonal entries. During execution,
! it is used for both supervariables and elements:
! * Principal supervariable I:  index into IW of the
! list of supervariable I.  A supervariable
! represents one or more rows of the matrix
! with identical pattern.
! * Non-principal supervariable I:  if I has been absorbed
! into another supervariable J, then PE(I) = -J.
! That is, J has the same pattern as I.
! Note that J might later be absorbed into another
! supervariable J2, in which case PE(I) is still -J,
! and PE(J) = -J2.
! * Unabsorbed element E:  the index into IW of the list
! of element E.  Element E is created when
! the supervariable of the same name is selected as
! the pivot.
! * Absorbed element E:  if element E is absorbed into element
! E2, then PE(E) = -E2.  This occurs when one of its
! variables is eliminated and when the pattern of
! E (that is, Le) is found to be a subset of the pattern
! of E2 (that is, Le2).  If element E is "null" (it has
! no entries outside its pivot block), then PE(E) = 0.

! On output, PE holds the assembly tree/forest, which implicitly
! represents a pivot order with identical fill-in as the actual
! order (via a depth-first search of the tree). If NV(I) .GT. 0,
! then I represents a node in the assembly tree, and the parent of
! I is -PE(I), or zero if I is a root. If NV(I)=0, then (I,-PE(I))
! represents an edge in a subtree, the root of which is a node in
! the assembly tree.

! PFREE must be set to the position in IW of the first free variable.
! During execution, additional data is placed in IW, and PFREE is
! modified so that components  of IW from PFREE are free.
! On output, PFREE is set equal to the size of IW that would have
! caused no compressions to occur.  If NCMPA is zero, then
! PFREE (on output) is less than or equal to IWLEN, and the space
! IW(PFREE+1 ... IWLEN) was not used. Otherwise, PFREE (on output)
! is greater than IWLEN, and all the memory in IW was used.

! LEN(I) must be set to hold the number of entries in row I of the
! matrix, excluding the diagonal.  The contents of LEN(1..N) are
! undefined on output.

! IW(1..PFREE-1) must be set to  hold the patterns of the rows of
! the matrix.  The matrix must be symmetric, and both upper and
! lower triangular parts must be present.  The diagonal must not be
! present.  Row I is held as follows:
! IW(PE(I)...PE(I) + LEN(I) - 1) must hold the list of
! column indices for entries in row I (simple
! supervariables), excluding the diagonal.  All
! supervariables start with one row/column each
! (supervariable I is just row I). If LEN(I) is zero on
! input, then PE(I) is ignored on input. Note that the
! rows need not be in any particular order, and there may
! be empty space between the rows.
! During execution, the supervariable I experiences fill-in. This
! is represented by constructing a list of the elements that cause
! fill-in in supervariable I:
! IE(PE(i)...PE(I) + ELEN(I) - 1) is the list of elements
! that contain I. This list is kept short by removing
! absorbed elements. IW(PE(I)+ELEN(I)...PE(I)+LEN(I)-1)
! is the list of supervariables in I. This list is kept
! short by removing nonprincipal variables, and any entry
! J that is also contained in at least one of the
! elements in the list for I.
! When supervariable I is selected as pivot, we create an element E
! of the same name (E=I):
! IE(PE(E)..PE(E)+LEN(E)-1) is the list of supervariables
! in element E.
! An element represents the fill-in that occurs when supervariable
! I is selected as pivot.
! CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
! The contents of IW are undefined on output.

! ELEN(I) need not be set. See the description of IW above. At the
! start of execution, ELEN(I) is set to zero. For a supervariable,
! ELEN(I) is the number of elements in the list for supervariable
! I. For an element, ELEN(E) is the negation of the position in the
! pivot sequence of the supervariable that generated it. ELEN(I)=0
! if I is nonprincipal.
! On output ELEN(1..N) holds the inverse permutation (the same
! as the 'INVP' argument in Sparspak). That is, if K = ELEN(I),
! then row I is the Kth pivot row.  Row I of A appears as the
! (ELEN(I))-th row in the permuted matrix, PAP^T.

! ICNTL is an INTEGER array of length 10 that contains control
! parameters and must be set by the user. 
!     ICNTL(1) - ICNTL(3) are not used.

!     ICNTL(4) controls the choice of AMD algorithm
!   <= 0 No checking for dense rows performed
!   > 0 Corresponds to automatic setting of the minimum density
!     requirement. mc68_order sets this to 1.

!     ICNTL(5) defines the largest positive
!     integer that your computer can represent (-iovflo should also
!     be representable). HUGE(1) in Fortran 95. 

! Local arrays:
! ---------------

! NV(I) During execution, ABS(NV(I)) is equal to the
! number of rows represented by the principal supervariable I. If I
! is a nonprincipal variable, then NV(I) = 0. Initially, NV(I) = 1
! for all I.  NV(I) .LT. 0 signifies that I is a principal variable
! in the pattern Lme of the current pivot element ME. On termination,
! NV(E) holds the true degree of element E at the time it was
! created (including the diagonal part).

! LAST(I) In a degree list, LAST(I) is the
! supervariable preceding I, or zero if I is the head of the list.
! In a hash bucket, LAST(I) is the hash key for I. LAST(HEAD(HASH))
! is also used as the head of a hash bucket if HEAD(HASH) contains
! a degree list (see HEAD, below).
! On output, LAST(1..N) holds the permutation (the same as the
! 'PERM' argument in Sparspak). That is, if I = LAST(K), then row I
! is the Kth pivot row.  Row LAST(K) of A is the K-th row in the
! permuted matrix, PAP^T.


! DEGREE If I is a supervariable and sparse,
! then DEGREE(I) holds the current approximation of the external
! degree of row I (an upper bound). The external degree is the
! number of entries in row I, minus ABS(NV(I)) (the diagonal
! part). The bound is equal to the external degree if ELEN(I) is
! less than or equal to two. We also use the term "external degree"
! for elements E to refer to |Le \ Lme|. If I is full in the reduced
! matrix, then DEGREE(I)=N+1. If I is dense in the reduced matrix,
! then DEGREE(I)=N+1+last_approximate_external_deg of I.
! All dense rows are stored in the list pointed by HEAD(N).
! Quasi dense rows are stored first, and are followed by full rows
! in the reduced matrix. LASTD holds the last row in
! this list of dense rows or is zero if the list is empty.

! HEAD(DEG) is used for degree lists.
! HEAD(DEG) is the first supervariable in a degree list (all
! supervariables I in a degree list DEG have the same approximate
! degree, namely, DEG = DEGREE(I)). If the list DEG is empty then
! HEAD(DEG) = 0.
! During supervariable detection HEAD(HASH) also serves as a
! pointer to a hash bucket.
! If HEAD(HASH) .GT. 0, there is a degree list of degree HASH. The
! hash bucket head pointer is LAST(HEAD(HASH)).
! If HEAD(HASH) = 0, then the degree list and hash bucket are
! both empty.
! If HEAD(HASH) .LT. 0, then the degree list is empty, and
! -HEAD(HASH) is the head of the hash bucket.
! After supervariable detection is complete, all hash buckets are
! empty, and the (LAST(HEAD(HASH)) = 0) condition is restored for
! the non-empty degree lists.

! DENXT(I)  For supervariable I, DENXT(I) is
! the supervariable following I in a link list, or zero if I is
! the last in the list. Used for two kinds of lists: degree lists
! and hash buckets (a supervariable can be in only one kind of
! list at a time). For element E, DENXT(E) is the number of
! variables with dense or full rows in the element E.

! W(I) The flag array W determines the status
! of elements and variables, and the external degree of elements.
! For elements:
! if W(E) = 0, then the element E is absorbed.
! if W(E) .GE. WFLG, then W(E)-WFLG is the size of the set
! |Le \ Lme|, in terms of nonzeros (the sum of ABS(NV(I))
! for each principal variable I that is both in the
! pattern of element E and NOT in the pattern of the
! current pivot element, ME).
! if WFLG .GT. WE(E) .GT. 0, then E is not absorbed and has
! not yet been seen in the scan of the element lists in
! the computation of |Le\Lme| in loop 150 below.
! ***SD: change comment to remove reference to label***
! For variables:
! during supervariable detection, if W(J) .NE. WFLG then J is
! not in the pattern of variable I.
! The W array is initialized by setting W(I) = 1 for all I, and by
! setting WFLG = 2. It is reinitialized if WFLG becomes too large
! (to ensure that WFLG+N does not cause integer overflow).


! Local variables:
! ---------------

! DEG:        the degree of a variable or element
! DEGME:      size (no. of variables), |Lme|, of the current element,
! ME (= DEGREE(ME))
! DEXT:       external degree, |Le \ Lme|, of some element E
! DMAX:       largest |Le| seen so far
! E:          an element
! ELENME:     the length, ELEN(ME), of element list of pivotal var.
! ELN:        the length, ELEN(...), of an element list
! EMP1:       stores number of rows found to be empty
! HASH:       the computed value of the hash function
! HMOD:       the hash function is computed modulo HMOD = MAX(1,N-1)
! I:          a supervariable
! IDUMMY:     loop counter
! ILAST:      the entry in a link list preceding I
! INEXT:      the entry in a link list following I
! IOVFLO:     local copy of ICNTL(5)
! J:          a supervariable
! JDUMMY:     loop counter
! JLAST:      the entry in a link list preceding J
! JNEXT:      the entry in a link list, or path, following J
! K:          the pivot order of an element or variable
! KNT1:       loop counter used during element construction
! KNT2:       loop counter used during element construction
! KNT3:       loop counter used during element construction
! LASTD:      index of the last row in the list of dense rows
! LENJ:       LEN(J)
! LN:         length of a supervariable list
! MAXMEM:     amount of memory needed for no compressions
! ME:         current supervariable being eliminated, and the
! current element created by eliminating that
! supervariable
! MEM:        memory in use assuming no compressions have occurred
! MINDEG:     current approximate minimum degree
! MU:         current mean of external degrees during dense detection
! NBD:        total number of dense rows selected
! NCMPA:      counter for the number of times IW was compressed
! NEL:        number of pivots selected so far
! NEWMEM:     amount of new memory needed for current pivot element
! NLEFT:      N-NEL, the number of nonpivotal rows/columns remaining
! NRLADU:     counter for the forecast number of reals in matrix factor
! NVI:        the number of variables in a supervariable I (= NV(I))
! NVJ:        the number of variables in a supervariable J (= NV(J))
! NVPIV:      number of pivots in current element
! P:          pointer into lots of things
! P1:         pe (i) for some variable i (start of element list)
! P2:         pe (i) + elen (i) -  1 for some var. i (end of el. list)
! P3:         index of first supervariable in clean list
! PJ:         pointer into an element or variable
! PDST:       destination pointer, for compression
! PEND:       end of memory to compress
! PME:        pointer into the current element (PME1...PME2)
! PME1:       the current element, ME, is stored in IW(PME1...PME2)
! PME2:       the end of the current element
! PN:         pointer into a "clean" variable, also used to compress
! PSRC:       source pointer, for compression
! SDEN:       used to remember whether dense rows occur
! SLENME:     number of variables in variable list of pivotal variable
! THRESH:     local copy of ICNTL(4)
! THRESM :    local integer holding the threshold used to detect quasi
! dense rows. When quasi dense rows are reintegrated in the
! graph to be processed then THRESM is modified.
! WE:         W(E)
! WFLG:       used for flagging the W array.  See description of W.
! WNVI:       WFLG-NV(I)
! X:          either a supervariable or an element

! OPS:        counter for forecast number of flops
! RELDEN :    holds average density to set THRESM automatically

! IDENSE is true if supervariable I is dense

! -------------------------------------------------------------------
! FUNCTIONS CALLED:
! -------------------------------------------------------------------

! ====================================================================
! INITIALIZATIONS
! ====================================================================

! ..
! .. Local Arrays ..
        INTEGER, ALLOCATABLE, DIMENSION (:) :: nv, last, degree, head, denxt, &
          w

! ..
! .. Local Scalars ..
        REAL (myreal_mc68) mu, relden
        INTEGER deg, degme, dext, dmax, e, elenme, eln, emp1, hash, hmod, i, &
          ilast, inext, iovflo, j, jdummy, jlast, jnext, k, knt1, &
          knt2, knt3, lastd, lenj, ln, maxmem, me, mem, mindeg, nbd, ncmpa, &
          nel, newmem, nleft, nvi, nvj, nvpiv, p, p1, p2, p3, pdst, &
          pend, pj, pme, pme1, pme2, pn, psrc, sden, slenme, thresh, thresm, &
          we, wflg, wnvi, x
        INTEGER(long) :: idummy
! ..
! .. Intrinsic Functions ..
        INTRINSIC abs, int, log, max, min, mod, real, sqrt
! ..

        dmax = 0
        hmod = max(1,n-1)
        iovflo = icntl(5)
        lastd = 0
        mem = pfree - 1
        maxmem = mem
        mindeg = 1
        nbd = 0
        ncmpa = 0
        nel = 0
        thresh = icntl(4)
        wflg = 2
        emp1 = 0
        sden = 0

        ALLOCATE (nv(n),last(n),degree(n),head(n),denxt(n),w(n), &
          STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_alloc
          RETURN
        END IF

        IF (thresh>0) THEN
          thresm = int(sqrt(real(n)))

! ---------------------------------------------------------------------
! initialize head(n) and next(n)
! ---------------------------------------------------------------------
          relden = 0
          head(1:n) = 0
          denxt(1:n) = 0

! ---------------------------------------------------------------------
! create the degree hash buckets and linked lists
! for the dense nodes
! ---------------------------------------------------------------------
          relden = 0
          degree(1:n) = len(1:n)
          DO i = 1, n
            deg = len(i)
            IF (deg>0) THEN
              sden = 1
              relden = relden + deg
! insert node in degree list
              denxt(i) = head(deg)
              head(deg) = i
            ELSE
              emp1 = emp1 + 1
            END IF
          END DO

          IF (n==emp1) THEN
            mu = 0
          ELSE
            mu = relden/(n-emp1)
          END IF

! ---------------------------------------------------------------------

! 1) Recalculate the degree length of all nodes adjacent to
! the dense nodes in the degree list.  (Note:  Many of the
! dense nodes in the degree list will no longer be dense after
! this section.)

! 2) Constuct the ordering for the nodes not sent to AMD by
! selecting the most dense node in the degree list and
! then reduce the lengths of all adjacent nodes. Repeat this
! until no nodes are left with length higher than dense.
! The dense nodes are placed in the last(n) array.
! NOTE:  1) nodes are placed after the final value
! of lastnode in the last(n) array
! 2) the AMD routine will not effect anything after
! last node in the last(n) array.
! 3) nodes are saved in degree order and in their
! original state, i.e., no reverse mapping is
! needed on these.
! ---------------------------------------------------------------------

          IF (sden==1) THEN
            sden = 0
            k = n
10          CONTINUE

! ** get node from bucket
            me = head(k)

! ** main loop control
            IF (mu/=zero .AND. (n-1/=sden+emp1) .AND. (n-2/=sden+emp1)) THEN
              IF (me==0) THEN
                k = k - 1
                IF ((mu-((relden-2*k)/real(n-sden-1-emp1))>=(40.0_wp* &
                  log(real(n-sden-1-emp1)))/real(n-sden-1-emp1)+mu/real(n- &
                  sden-2-emp1)) .AND. k/=0) GO TO 10
              ELSE

! ** remove node from bucket
                head(k) = denxt(me)

! ** get degree of current node
                deg = degree(me)

! ** skip this node if degree was changed to less than dense
                IF (deg>=1) THEN

! ** check if degree was changed
                  IF (deg<k .AND. deg>0) THEN

! ** insert back into linked list at the lower degree
                    denxt(me) = head(deg)
                    head(deg) = me
                  ELSE

! ** update degree lengths of adjacent nodes
                    p1 = pe(me) + len(me) - 1
                    DO i = pe(me), p1
                      j = iw(i)
                      IF (degree(j)==1) THEN
                        degree(j) = 0
                        emp1 = emp1 + 1
                      ELSE IF (degree(j)<2*n+1) THEN
                        degree(j) = degree(j) - 1
                      END IF
                    END DO
                    sden = sden + 1
                    relden = relden - 2*deg
                    IF (n==sden+emp1) THEN
                      mu = 0
                    ELSE
                      mu = relden/(n-sden-emp1)
                    END IF
! Mark as dense row
                    degree(me) = deg + 2*n + 1
                    pe(me) = 0
                  END IF
                END IF
                GO TO 10
              END IF
            END IF
          END IF

          DO i = 1, n
            IF (degree(i)<2*n+1) THEN
              IF (degree(i)>0) THEN
                p1 = pe(i)
                p2 = pe(i) + len(i) - 1
                pj = p1
                DO j = p1, p2
! Search I for dense variables
                  IF (degree(iw(j))>2*n+1) THEN
                    len(i) = len(i) - 1
                  ELSE
                    iw(pj) = iw(j)
                    pj = pj + 1
                  END IF
                END DO
                degree(i) = len(i)
                IF (len(i)==0) pe(i) = 0
                len(i) = pj - p1
              ELSE
                len(i) = 0
                pe(i) = 0
              END IF
            END IF
          END DO
        END IF

        IF (thresh>0) THEN
! ----------------------------------------------------------
! initialize arrays and eliminate rows with no off-diag. nz.
! ----------------------------------------------------------
          last(1:n) = 0
          head(1:n) = 0
          denxt(1:n) = 0
          nv(1:n) = 1
          DO i = 1, n
            IF (len(i)==0 .AND. degree(i)==0) THEN
              nel = nel + 1
              elen(i) = -nel
              pe(i) = 0
              w(i) = 0
            ELSE
              w(i) = 1
              elen(i) = 0
            END IF
          END DO

          IF (n==nel) THEN
            GO TO 100
          ELSE
            thresm = n + 1
          END IF
        ELSE

          thresm = thresh
          last(1:n) = 0
          head(1:n) = 0
          nv(1:n) = 1
          degree(1:n) = len(1:n)
          DO i = 1, n
            IF (degree(i)==0) THEN
              nel = nel + 1
              elen(i) = -nel
              pe(i) = 0
              w(i) = 0
            ELSE
              w(i) = 1
              elen(i) = 0
            END IF
          END DO
        END IF

! ----------------------------------------------------------------
! initialize degree lists
! ----------------------------------------------------------------
        DO i = 1, n
          deg = degree(i)
          IF (deg>0) THEN
! ----------------------------------------------------------
! place i in the degree list corresponding to its degree
! or in the dense row list if i is dense
! ----------------------------------------------------------
! test for row density
            IF ((thresm>=0) .AND. (deg>=2*n+1)) THEN
! I is dense and will be inserted in the degree
! list of N
              deg = n
              inext = head(deg)
              IF (inext/=0) last(inext) = i
              denxt(i) = inext
              head(deg) = i
              last(i) = 0
              IF (lastd==0) lastd = i
            ELSE
! place i in the degree list corresponding to its degree
              inext = head(deg)
              IF (inext/=0) last(inext) = i
              denxt(i) = inext
              head(deg) = i
            END IF
          END IF
        END DO
        nbd = sden

! We suppress dense row selection if none of them was found in A
! in the 1st pass
        IF (nbd==0 .AND. thresh>0) thresm = -1

        DO WHILE (nel<n)

! ==================================================================
! GET PIVOT OF MINIMUM APPROXIMATE DEGREE
! ==================================================================
! -------------------------------------------------------------
! find next supervariable for elimination
! -------------------------------------------------------------
          DO deg = mindeg, n
            me = head(deg)
            IF (me>0) GO TO 20
          END DO
20        mindeg = deg

          IF (deg<n) THEN
! -------------------------------------------------------------
! remove chosen variable from linked list
! -------------------------------------------------------------
            inext = denxt(me)
            IF (inext/=0) last(inext) = 0
            head(deg) = inext
! -------------------------------------------------------------
! me represents the elimination of pivots nel+1 to nel+nv(me).
! place me itself as the first in this set.  It will be moved
! to the nel+nv(me) position when the permutation vectors are
! computed.
! -------------------------------------------------------------
            elenme = elen(me)
            elen(me) = -(nel+1)
            nvpiv = nv(me)
            nel = nel + nvpiv
            denxt(me) = 0

! ====================================================================
! CONSTRUCT NEW ELEMENT
! ====================================================================

! -------------------------------------------------------------
! At this point, me is the pivotal supervariable.  It will be
! converted into the current element.  Scan list of the
! pivotal supervariable, me, setting tree pointers and
! constructing new list of supervariables for the new element,
! me.  p is a pointer to the current position in the old list.
! -------------------------------------------------------------

! flag the variable "me" as being in the front by negating nv(me)
            nv(me) = -nvpiv
            degme = 0
            IF (elenme==0) THEN
! ----------------------------------------------------------
! There are no elements involved.
! Construct the new element in place.
! ----------------------------------------------------------
              pme1 = pe(me)
              pme2 = pme1 - 1
              DO p = pme1, pme1 + len(me) - 1
                i = iw(p)
                nvi = nv(i)
                IF (nvi>0) THEN
! ----------------------------------------------------
! i is a principal variable not yet placed in the
! generated element. Store i in new list
! ----------------------------------------------------
                  degme = degme + nvi
! flag i as being in Lme by negating nv (i)
                  nv(i) = -nvi
                  pme2 = pme2 + 1
                  iw(pme2) = i

! ----------------------------------------------------
! remove variable i from degree list.
! ----------------------------------------------------
                  ilast = last(i)
                  inext = denxt(i)
                  IF (inext/=0) last(inext) = ilast
                  IF (ilast/=0) THEN
                    denxt(ilast) = inext
                  ELSE
! i is at the head of the degree list
                    head(degree(i)) = inext
                  END IF
                END IF
              END DO
! this element takes no new memory in iw:
              newmem = 0
            ELSE
! ----------------------------------------------------------
! construct the new element in empty space, iw (pfree ...)
! ----------------------------------------------------------
              p = pe(me)
              pme1 = pfree
              slenme = len(me) - elenme
              DO knt1 = 1, elenme
! search the elements in me.
                e = iw(p)
                p = p + 1
                pj = pe(e)
                ln = len(e)
! -------------------------------------------------------
! search for different supervariables and add them to the
! new list, compressing when necessary.
! -------------------------------------------------------
                DO knt2 = 1, ln
                  i = iw(pj)
                  pj = pj + 1
                  nvi = nv(i)
                  IF (nvi>0) THEN
! -------------------------------------------------
! compress iw, if necessary
! -------------------------------------------------
                    IF (pfree>iwlen) THEN
! prepare for compressing iw by adjusting
! pointers and lengths so that the lists being
! searched in the inner and outer loops contain
! only the remaining entries.
! ***** SD: Seperate compression subroutine tried
! but found to be inefficient in comparison ****
                      pe(me) = p
                      len(me) = len(me) - knt1
! Check if anything left in supervariable ME
                      IF (len(me)==0) pe(me) = 0
                      pe(e) = pj
                      len(e) = ln - knt2
! Check if anything left in element E
                      IF (len(e)==0) pe(e) = 0
                      ncmpa = ncmpa + 1
! store first item in pe
! set first entry to -item
                      DO j = 1, n
                        pn = pe(j)
                        IF (pn>0) THEN
                          pe(j) = iw(pn)
                          iw(pn) = -j
                        END IF
                      END DO

! psrc/pdst point to source/destination
                      pdst = 1
                      psrc = 1
                      pend = pme1 - 1

! while loop:
                      DO idummy = 1, iwlen
                        IF (psrc>pend) THEN
                          GO TO 30
                        ELSE
! search for next negative entry
                          j = -iw(psrc)
                          psrc = psrc + 1
                          IF (j>0) THEN
                            iw(pdst) = pe(j)
                            pe(j) = pdst
                            pdst = pdst + 1
! copy from source to destination
                            lenj = len(j)
                            DO knt3 = 0, lenj - 2
                              iw(pdst+knt3) = iw(psrc+knt3)
                            END DO
                            pdst = pdst + lenj - 1
                            psrc = psrc + lenj - 1
                          END IF
                        END IF
                      END DO

! move the new partially-constructed element
30                    p1 = pdst
                      DO psrc = pme1, pfree - 1
                        iw(pdst) = iw(psrc)
                        pdst = pdst + 1
                      END DO
                      pme1 = p1
                      pfree = pdst
                      pj = pe(e)
                      p = pe(me)
                    END IF

! -------------------------------------------------
! i is a principal variable not yet placed in Lme
! store i in new list
! -------------------------------------------------
                    degme = degme + nvi
! flag i as being in Lme by negating nv (i)
                    nv(i) = -nvi
                    iw(pfree) = i
                    pfree = pfree + 1

! -------------------------------------------------
! remove variable i from degree link list
! -------------------------------------------------
                    ilast = last(i)
                    inext = denxt(i)
                    IF (inext/=0) last(inext) = ilast
                    IF (ilast/=0) THEN
                      denxt(ilast) = inext
                    ELSE
! i is at the head of the degree list
                      head(degree(i)) = inext
                    END IF
                  END IF
                END DO

! set tree pointer and flag to indicate element e is
! absorbed into new element me (the parent of e is me)
                pe(e) = -me
                w(e) = 0
              END DO

! search the supervariables in me.
              knt1 = elenme + 1
              e = me
              pj = p
              ln = slenme

! -------------------------------------------------------
! search for different supervariables and add them to the
! new list, compressing when necessary.
! -------------------------------------------------------
              DO knt2 = 1, ln
                i = iw(pj)
                pj = pj + 1
                nvi = nv(i)
                IF (nvi>0) THEN
! -------------------------------------------------
! compress iw, if necessary
! -------------------------------------------------
                  IF (pfree>iwlen) THEN
! prepare for compressing iw by adjusting
! pointers and lengths so that the lists being
! searched in the inner and outer loops contain
! only the remaining entries.
                    pe(me) = p
                    len(me) = len(me) - knt1
! Check if anything left in supervariable ME
                    IF (len(me)==0) pe(me) = 0
                    pe(e) = pj
                    len(e) = ln - knt2
! Check if anything left in element E
                    IF (len(e)==0) pe(e) = 0
                    ncmpa = ncmpa + 1
! store first item in pe
! set first entry to -item
                    DO j = 1, n
                      pn = pe(j)
                      IF (pn>0) THEN
                        pe(j) = iw(pn)
                        iw(pn) = -j
                      END IF
                    END DO

! psrc/pdst point to source/destination
                    pdst = 1
                    psrc = 1
                    pend = pme1 - 1

! while loop:
! 122              CONTINUE
                    DO idummy = 1, iwlen
                      IF (psrc>pend) THEN
                        GO TO 40
                      ELSE
! search for next negative entry
                        j = -iw(psrc)
                        psrc = psrc + 1
                        IF (j>0) THEN
                          iw(pdst) = pe(j)
                          pe(j) = pdst
                          pdst = pdst + 1
! copy from source to destination
                          lenj = len(j)
                          DO knt3 = 0, lenj - 2
                            iw(pdst+knt3) = iw(psrc+knt3)
                          END DO
                          pdst = pdst + lenj - 1
                          psrc = psrc + lenj - 1
                        END IF
                      END IF
                    END DO

! move the new partially-constructed element
40                  p1 = pdst
                    DO psrc = pme1, pfree - 1
                      iw(pdst) = iw(psrc)
                      pdst = pdst + 1
                    END DO
                    pme1 = p1
                    pfree = pdst
                    pj = pe(e)
                    p = pe(me)
                  END IF

! -------------------------------------------------
! i is a principal variable not yet placed in Lme
! store i in new list
! -------------------------------------------------
                  degme = degme + nvi
! flag i as being in Lme by negating nv (i)
                  nv(i) = -nvi
                  iw(pfree) = i
                  pfree = pfree + 1

! -------------------------------------------------
! remove variable i from degree link list
! -------------------------------------------------
                  ilast = last(i)
                  inext = denxt(i)
                  IF (inext/=0) last(inext) = ilast
                  IF (ilast/=0) THEN
                    denxt(ilast) = inext
                  ELSE
! i is at the head of the degree list
                    head(degree(i)) = inext
                  END IF
                END IF
              END DO

              pme2 = pfree - 1
! this element takes newmem new memory in iw (possibly zero)
              newmem = pfree - pme1
              mem = mem + newmem
              maxmem = max(maxmem,mem)
            END IF

! -------------------------------------------------------------
! me has now been converted into an element in iw (pme1..pme2)
! -------------------------------------------------------------
! degme holds the external degree of new element
            degree(me) = degme
            pe(me) = pme1
            len(me) = pme2 - pme1 + 1

! -------------------------------------------------------------
! make sure that wflg is not too large.  With the current
! value of wflg, wflg+n must not cause integer overflow
! -------------------------------------------------------------
            IF (wflg>iovflo-n) THEN
              DO x = 1, n
                IF (w(x)/=0) w(x) = 1
              END DO
              wflg = 2
            END IF

! ====================================================================
! COMPUTE (w(e) - wflg) = |Le(G')\Lme(G')| FOR ALL ELEMENTS
! where G' is the subgraph of G containing just the sparse rows)
! ====================================================================
! -------------------------------------------------------------
! Scan 1:  compute the external degrees of elements touched
! with respect to the current element.  That is:
! (w (e) - wflg) = |Le \ Lme|
! for each element e involving a supervariable in Lme.
! The notation Le refers to the pattern (list of
! supervariables) of a previous element e, where e is not yet
! absorbed, stored in iw (pe (e) + 1 ... pe (e) + iw (pe (e))).
! The notation Lme refers to the pattern of the current element
! (stored in iw (pme1..pme2)).
! -------------------------------------------------------------
            DO pme = pme1, pme2
              i = iw(pme)
              eln = elen(i)
              IF (eln>0) THEN
! note that nv (i) has been negated to denote i in Lme:
                nvi = -nv(i)
                wnvi = wflg - nvi
                DO p = pe(i), pe(i) + eln - 1
                  e = iw(p)
                  we = w(e)
                  IF (we>=wflg) THEN
! unabsorbed element e has been seen in this loop
                    we = we - nvi
                  ELSE IF (we/=0) THEN
! e is an unabsorbed element - this is
! the first we have seen e in all of Scan 1
                    we = degree(e) + wnvi
                  END IF
                  w(e) = we
                END DO
              END IF
            END DO

! ====================================================================
! DEGREE UPDATE AND ELEMENT ABSORPTION
! ====================================================================

! -------------------------------------------------------------
! Scan 2:  for each sparse i in Lme, sum up the external degrees
! of each Le for the elements e appearing within i, plus the
! supervariables in i.  Place i in hash list.
! -------------------------------------------------------------

            DO pme = pme1, pme2
              i = iw(pme)
! remove absorbed elements from the list for i
              p1 = pe(i)
              p2 = p1 + elen(i) - 1
              pn = p1
              hash = 0
              deg = 0

! -------------------------------------------------------
! scan the element list associated with supervariable i
! -------------------------------------------------------
              DO p = p1, p2
                e = iw(p)
! dext = | Le | - | (Le \cap Lme)\D | - DENXT(e)
                dext = w(e) - wflg
                IF (dext>0) THEN
                  deg = deg + dext
                  iw(pn) = e
                  pn = pn + 1
                  hash = hash + e
                ELSE IF (dext==0) THEN
! aggressive absorption: e is not adjacent to me, but
! |Le(G') \ Lme(G')| is 0, so absorb it into me
                  pe(e) = -me
                  w(e) = 0
                END IF
              END DO

! count the number of elements in i (including me):
              elen(i) = pn - p1 + 1

! ----------------------------------------------------------
! scan the supervariables in the list associated with i
! ----------------------------------------------------------
              p3 = pn
              DO p = p2 + 1, p1 + len(i) - 1
                j = iw(p)
                nvj = nv(j)
                IF (nvj>0) THEN
! j is unabsorbed, and not in Lme.
! add to degree and add to new list
                  deg = deg + nvj
                  iw(pn) = j
                  pn = pn + 1
                  hash = hash + j
                END IF
              END DO

! ----------------------------------------------------------
! update the degree and check for mass elimination
! ----------------------------------------------------------
              IF (deg==0) THEN
! -------------------------------------------------------
! mass elimination - supervariable i can be eliminated
! -------------------------------------------------------
                pe(i) = -me
                nvi = -nv(i)
                degme = degme - nvi
                nvpiv = nvpiv + nvi
                nel = nel + nvi
                nv(i) = 0
                elen(i) = 0
              ELSE
! -------------------------------------------------------
! update the upper-bound degree of i
! A bound for the new external degree is the old bound plus
! the size of the generated element
! -------------------------------------------------------

! the following degree does not yet include the size
! of the current element, which is added later:
                degree(i) = min(deg,degree(i))

! -------------------------------------------------------
! add me to the list for i
! -------------------------------------------------------
! move first supervariable to end of list
                iw(pn) = iw(p3)
! move first element to end of element part of list
                iw(p3) = iw(p1)
! add new element to front of list.
                iw(p1) = me
! store the new length of the list in len (i)
                len(i) = pn - p1 + 1

! -------------------------------------------------------
! place in hash bucket.  Save hash key of i in last (i).
! -------------------------------------------------------
                hash = abs(mod(hash,hmod)) + 1
                j = head(hash)
                IF (j<=0) THEN
! the degree list is empty, hash head is -j
                  denxt(i) = -j
                  head(hash) = -i
                ELSE
! degree list is not empty - has j as its head
! last is hash head
                  denxt(i) = last(j)
                  last(j) = i
                END IF
                last(i) = hash
              END IF
            END DO
            degree(me) = degme

! -------------------------------------------------------------
! Clear the counter array, w (...), by incrementing wflg.
! -------------------------------------------------------------
            dmax = max(dmax,degme)
            wflg = wflg + dmax

! make sure that wflg+n does not cause integer overflow
            IF (wflg>=iovflo-n) THEN
              DO x = 1, n
                IF (w(x)/=0) w(x) = 1
              END DO
              wflg = 2
            END IF
! at this point, w (1..n) .lt. wflg holds

! ====================================================================
! SUPERVARIABLE DETECTION
! ====================================================================
            DO pme = pme1, pme2
              i = iw(pme)
              IF ((nv(i)<0) .AND. (degree(i)<=n)) THEN
! only done for sparse rows
! replace i by head of its hash bucket, and set the hash
! bucket header to zero

! -------------------------------------------------------
! examine all hash buckets with 2 or more variables.  We
! do this by examing all unique hash keys for super-
! variables in the pattern Lme of the current element, me
! -------------------------------------------------------
                hash = last(i)
! let i = head of hash bucket, and empty the hash bucket
                j = head(hash)
                IF (j/=0) THEN
                  IF (j<0) THEN
! degree list is empty
                    i = -j
                    head(hash) = 0
                  ELSE
! degree list is not empty, restore last () of head
                    i = last(j)
                    last(j) = 0
                  END IF
                  IF (i/=0) THEN

! while loop:
                    DO jdummy = 1, n
                      IF (denxt(i)==0) THEN
                        GO TO 80
                      ELSE
! ----------------------------------------------------
! this bucket has one or more variables following i.
! scan all of them to see if i can absorb any entries
! that follow i in hash bucket.  Scatter i into w.
! ----------------------------------------------------
                        ln = len(i)
                        eln = elen(i)
! do not flag the first element in the list (me)
                        DO p = pe(i) + 1, pe(i) + ln - 1
                          w(iw(p)) = wflg
                        END DO

! ----------------------------------------------------
! scan every other entry j following i in bucket
! ----------------------------------------------------
                        jlast = i
                        j = denxt(i)

! while loop:
                        DO idummy = 1, n
                          IF (j==0) THEN
                            GO TO 70
                          ELSE

! -------------------------------------------------
! check if j and i have identical nonzero pattern
! -------------------------------------------------
! jump if i and j do not have same size data structure
                            IF (len(j)==ln) THEN
! jump if i and j do not have same number adj elts
                              IF (elen(j)==eln) THEN
! do not flag the first element in the list (me)

                                DO p = pe(j) + 1, pe(j) + ln - 1
! jump if an entry (iw(p)) is in j but not in i
                                  IF (w(iw(p))/=wflg) GO TO 50
                                END DO

! -------------------------------------------------
! found it!  j can be absorbed into i
! -------------------------------------------------
                                pe(j) = -i
! both nv (i) and nv (j) are negated since they
! are in Lme, and the absolute values of each
! are the number of variables in i and j:
                                nv(i) = nv(i) + nv(j)
                                nv(j) = 0
                                elen(j) = 0
! delete j from hash bucket
                                j = denxt(j)
                                denxt(jlast) = j
                                GO TO 60
                              END IF
                            END IF

! -------------------------------------------------
50                          CONTINUE
! j cannot be absorbed into i
! -------------------------------------------------
                            jlast = j
                            j = denxt(j)
                          END IF
60                        CONTINUE
                        END DO

! ----------------------------------------------------
! no more variables can be absorbed into i
! go to next i in bucket and clear flag array
! ----------------------------------------------------
70                      wflg = wflg + 1
                        i = denxt(i)
                        IF (i==0) GO TO 80
                      END IF
                    END DO
                  END IF
                END IF
              END IF
80            CONTINUE
            END DO

! ====================================================================
! RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM ELEMENT
! Squeeze out absorbed variables
! ====================================================================
            p = pme1
            nleft = n - nel
            DO pme = pme1, pme2
              i = iw(pme)
              nvi = -nv(i)
              IF (nvi>0) THEN
! i is a principal variable in Lme
! restore nv (i) to signify that i is principal
                nv(i) = nvi
                IF (degree(i)<=n) THEN
! -------------------------------------------------------
! compute the external degree (add size of current elem)
! -------------------------------------------------------
                  deg = min(degree(i)+degme-nvi,nleft-nvi)
                  degree(i) = deg

! -------------------------------------------------------
! place the supervariable at the head of the degree list
! -------------------------------------------------------
                  inext = head(deg)
                  IF (inext/=0) last(inext) = i
                  denxt(i) = inext
                  last(i) = 0
                  head(deg) = i
! -------------------------------------------------------
! save the new degree, and find the minimum degree
! -------------------------------------------------------
                  mindeg = min(mindeg,deg)
                END IF
! -------------------------------------------------------
! place the supervariable in the element pattern
! -------------------------------------------------------
                iw(p) = i
                p = p + 1
              END IF
            END DO

! =====================================================================
! FINALIZE THE NEW ELEMENT
! =====================================================================
            nv(me) = nvpiv + degme
! nv (me) is now the degree of pivot (including diagonal part)
! save the length of the list for the new element me
            len(me) = p - pme1
            IF (len(me)==0) THEN
! there is nothing left of the current pivot element
              pe(me) = 0
              w(me) = 0
            END IF
            IF (newmem/=0) THEN
! element was not constructed in place: deallocate part
! of it (final size is less than or equal to newmem,
! since newly nonprincipal variables have been removed).
              pfree = p
              mem = mem - newmem + len(me)
            END IF

! =====================================================================
! END WHILE (selecting pivots)
          ELSE
! DEGREE(ME).GT.N+1 so ME is dense
! RESTARTING STRATEGY
! FOR EACH  dense row d
! 1/ insert d in the degree list according to the
! value degree(d)-(N+1) (updating MINDEG)
! 2/ ME is assumed to have no adjacent variables because just
! sorting according to order removed from matrix at initialisation
! 3/ get back to min degree process

! While loop: ME is the current dense row
! make sure that WFLG is not too large
            IF (wflg>iovflo-nbd-1) THEN
              DO x = 1, n
                IF (w(x)/=0) w(x) = 1
              END DO
              wflg = 2
            END IF
            wflg = wflg + 1
            DO idummy = 1, n
              pe(me) = 0
              len(me) = 0

! ---------------------------------------------------------
! remove chosen variable from link list
! ---------------------------------------------------------
              inext = denxt(me)
              IF (inext/=0) THEN
                last(inext) = 0
              ELSE
                lastd = 0
              END IF
! ----------------------------------------------------------
! build adjacency list of ME in quotient graph
! and calculate its external degree in ndense(me)
! ----------------------------------------------------------
              denxt(me) = 0
! Flag ME as having been considered in this calculation
              w(me) = wflg
              p1 = pe(me)
              p2 = p1 + len(me) - 1
! LN-1 holds the pointer in IW to last elt/var in adj list
! of ME.  LEN(ME) will then be set to LN-P1
! ELN-1 hold the pointer in IW to  last elt in in adj list
! of ME.  ELEN(ME) will then be set to ELN-P1
! element adjacent to ME
              ln = p1
              eln = p1

! ----------------------------------------------
! DEGREE(ME)-(2*N+1) holds last external degree computed
! when ME was detected as dense
! DENXT(ME) is the exact external degree of ME
! ----------------------------------------------
              wflg = wflg + 1
              len(me) = ln - p1
              elen(me) = eln - p1
              IF (denxt(me)==0) denxt(me) = 1
! ---------------------------------------------------------
! place ME in the degree list of DENXT(ME), update DEGREE
! ---------------------------------------------------------
              deg = degree(me) - 2*n - 1
              degree(me) = denxt(me)
              mindeg = min(deg,mindeg)
              jnext = head(deg)
              IF (jnext/=0) last(jnext) = me
              denxt(me) = jnext
              head(deg) = me
!  nel = nel+1

! ------------------------------
! process dense row
! ------------------------------
              me = inext
              IF (me==0) THEN
                GO TO 90
              ELSE IF (degree(me)<=(n+1)) THEN
                GO TO 90
              END IF
            END DO
90          head(n) = me
! get back to min degree elimination loop
          END IF
        END DO
! =====================================================================

100     CONTINUE
! ===================================================================
! COMPUTE THE PERMUTATION VECTORS
! ===================================================================

! ----------------------------------------------------------------
! The time taken by the following code is O(n).  At this
! point, elen (e) = -k has been done for all elements e,
! and elen (i) = 0 has been done for all nonprincipal
! variables i.  At this point, there are no principal
! supervariables left, and all elements are absorbed.
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! compute the ordering of unordered nonprincipal variables
! ----------------------------------------------------------------

        DO i = 1, n
          IF (elen(i)==0) THEN
! ----------------------------------------------------------
! i is an un-ordered row.  Traverse the tree from i until
! reaching an element, e.  The element, e, was the
! principal supervariable of i and all nodes in the path
! from i to when e was selected as pivot.
! ----------------------------------------------------------
            j = -pe(i)
! while (j is a variable) do:
            DO jdummy = 1, n
              IF (elen(j)<0) THEN
                GO TO 110
              ELSE
                j = -pe(j)
              END IF
            END DO
110         e = j
! ----------------------------------------------------------
! get the current pivot ordering of e
! ----------------------------------------------------------
            k = -elen(e)

! ----------------------------------------------------------
! traverse the path again from i to e, and compress the
! path (all nodes point to e).  Path compression allows
! this code to compute in O(n) time.  Order the unordered
! nodes in the path, and place the element e at the end.
! ----------------------------------------------------------
            j = i
! while (j is a variable) do:
            DO idummy = 1, n
              IF (elen(j)<0) THEN
                GO TO 120
              ELSE
                jnext = -pe(j)
                pe(j) = -e
                IF (elen(j)==0) THEN
! j is an unordered row
                  elen(j) = k
                  k = k + 1
                END IF
                j = jnext
              END IF
            END DO
! leave elen (e) negative, so we know it is an element
120         elen(e) = -k
          END IF
        END DO

! ----------------------------------------------------------------
! reset the inverse permutation (elen (1..n)) to be positive,
! and compute the permutation (last (1..n)).
! ----------------------------------------------------------------
        DO i = 1, n
          k = abs(elen(i))
          last(k) = i
          elen(i) = k
        END DO

! ====================================================================
! RETURN THE MEMORY USAGE IN IW AND SET INFORMATION ARRAYS
! ====================================================================
! If maxmem is less than or equal to iwlen, then no compressions
! occurred, and iw (maxmem+1 ... iwlen) was unused.  Otherwise
! compressions did occur, and iwlen would have had to have been
! greater than or equal to maxmem for no compressions to occur.
! Return the value of maxmem in the pfree argument.

        DEALLOCATE (nv,last,degree,head,denxt,w,STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_dealloc
          RETURN
        END IF

        info%n_compressions = ncmpa
        info%n_zero_eigs = -1
        info%n_dense_rows = nbd
        pfree = maxmem

END SUBROUTINE amdd

! -------------------------------------------------------------------


END MODULE hsl_mc68_integer

! Additional modules to provide backwards compatibility.
! These modules are deprecated and may be removed at a later date.
module hsl_mc68_double
   use hsl_mc68_integer
end module hsl_mc68_double
module hsl_mc68_single
   use hsl_mc68_integer
end module hsl_mc68_single
