! COPYRIGHT (c) 2007 University of Manchester and
!           Council for the Central Laboratory of the Research Councils

! Version 2.0.0
! See ChangeLog for version history

!  Purpose: AMG preconditioner. Given an sparse matrix A and an vector z,
!  HSL_MI20 computes the vector x = Mz, where M is an
!  algebraic multigrid (AMG) v-cycle preconditioner for
!  A. A classical AMG method is used.
!  A must have positive diagonal entries and (most of) the
!  off-diagonals entries must be negative (the diagonal should be large
!  compared to the sum of the off-diagonals).
!  During the multigrid
!  coarsening process, positive off-diagonal entries are ignored and, when
!  calculating the interpolation weights, positive off-diagonal entries
!  are added to the diagonal.

!  To generate the single precision version from the double:
!  (1) change _double to _single
!  (2) change kind(1.0d0)  to kind(1.0)
!  (3) replace dgetrf by sgetrf and dgetrs by sgetrs (lapack calls)

module hsl_mi20_double ! 11 sept 2007

 use hsl_zd11_double
 use hsl_MC65_double
 use hsl_ma48_double

 implicit none

  integer, parameter, private :: myreal = kind(1.0d0)
  real (myreal), parameter, private :: one = 1.0_myreal
  real (myreal), parameter, private :: zero = 0.0_myreal

 ! errors
 ! == == ==
 !
 !  matrix errors: -1 to -10
 ! ---------------
 ! -1 out of range entries in col
 ! -2 missing diagonal entry
 ! -3 diagonal entry <= 0
 ! -4 ptr wrong size
 ! -5 col wrong size
 ! -6 val wrong size
 ! -7 out of range values in ptr
 ! -8 duplicate entries in col
 ! -9 m < 1
 !
 ! errors from mi20 subroutines: -10 to -99
 ! ----------------------------
 ! -10 allocation error in first level of coarsening
 ! -11 deallocation error
 ! -12 coarsening failure due to positive rows
 ! -14 final vector norm > err_tol * initial vector norm
 ! -15 preconditioner called without successful setup
 ! -16 size of vectors smaller than size of matrix in preconditioner
 ! -17 failure in _getrf
 ! -18 failure in hsl_ma48
 ! -19 user requested hsl_ma48 solver without supplying hsl_ma48 control type
 !
 ! out of range mi20_control errors: -100 to -200
 ! --------------------------------
 ! -100 testing
 ! -101 st_parameter
 ! -102 err_tol
 ! -103 max_points
 ! -104 st_method
 ! -105 aggressive
 ! -106 c_fail
 ! -107 v_iterations
 ! -108 smoother
 ! -109 pre_smoothing
 ! -110 post_smoothing
 ! -111 no smoothing
 ! -112 coarse_solver
 ! -113 coarse_solver_its
 ! -114 print_level
 ! -115 damping
 ! -116 max_levels
 ! -117 ma48
 ! -118 trunc_parameter
 ! -119 reduction
 !
 ! warnings
 ! --------
 ! 1  st_method changed
 ! 10 coarsening stopped because of allocation error
 ! 11 coarsening stopped because of deallocation error
 ! 12 coarsening stopped because of positive rows
 ! 13 coarsening stagnation


 ! == == == == == == == == == == =
 ! ==  DERIVED TYPES  ==
 ! == == == == == == == == == == =


 ! --------------------
 ! --  mi20_control  --
 ! --------------------

 ! Derived type to hold control parameters for hsl_mi20
 type mi20_control

  ! number of coarsening steps per coarse level
  integer :: aggressive = 1

  ! determines the conditions for coarsening failure
  integer :: c_fail = 1

  ! size of mi20_data object
  integer :: max_levels = 100

  ! maximum number of points allowed on a coarse level
  integer :: max_points = 1

  ! controls when the coarsening is considered to have stagnated
  real (kind=myreal) :: reduction = 0.8

  ! defines the method used to find strong transpose connections
  integer :: st_method = 2

  ! parameter defines which connections are considered 'strong'
  real(kind=myreal) :: st_parameter = 0.25

  ! determines whether the matrix input to mi20_setup is tested for validity
  integer :: testing = 1

  ! interpolation truncation parameter
  real (kind=myreal) :: trunc_parameter = zero

  ! determines which coarse solver is used
  integer :: coarse_solver = 3

  ! determines the number of coarse solver iterations (for iterative methods)
  integer :: coarse_solver_its = 10

  ! damping factor for Jacobi smoother
  real(kind=myreal) :: damping = 0.8

  ! error tolerance for preconditioner
  real(kind=myreal) :: err_tol = 1.0e10

  ! determines the number of coarse levels used in the preconditioning
  integer :: levels = -1

  ! controls the use of ma48 setup by the preconditioner
  integer :: ma48 = 0

  ! number of pre-smoothing iterations
  integer :: pre_smoothing = 2

  ! controls the smoother type used
  integer :: smoother = 2

  ! number of post-smoothing iterations
  integer :: post_smoothing = 2

  ! number of AMG iterations when the preconditoner is called
  integer :: v_iterations = 1

  ! determines the levels of messages output
  integer :: print_level = 1

  ! determines the output stream for general messages
  integer :: print = 6

  ! determines the output_stream for error messages
  integer :: error = 6

  ! determines whether one pass coarsening is used
  logical :: one_pass_coarsen = .false.

!!  ! tolerance - NOT USED IN PRECONDITIONER
!!  real(kind=myreal) ::  tol = -one

!!  ! .true. if tol is relative to starting value - NOT USED IN PRECONDITIONER
!!  logical ::  tol_relative = .true.



 end type mi20_control


 ! -----------------
 ! --  mi20_info  --
 ! -----------------

 ! Derived type to communucate errors and information to the calling program.
 type mi20_info

  ! error/warning information
  integer :: flag = 0

  ! number of levels actually generated
  integer :: clevels = 0

  ! number of points on the coarsest level
  integer :: cpoints = 0

  ! number of non-zeros in coarsest matrix
  integer :: cnnz = 0

  ! holds the fortran stat parameter from calls to allocation and deallocation
  integer :: stat

  ! returns getrf info
  integer :: getrf_info

  ! number of iterations
  integer :: iterations
  
  ! norm of residual
  real(kind=myreal) :: residual

  ! ma48 derived types used to return information
  type(ma48_ainfo):: ma48_ainfo
  type(ma48_finfo):: ma48_finfo
  type(ma48_sinfo):: ma48_sinfo

 end type mi20_info


 ! -----------------
 ! --  mi20_keep  --
 ! -----------------

 ! Derived type to hold internal arrays and parameters for hsl_mi20 package.
 type mi20_keep

  ! set to .true. when have new preconditioner
  logical :: new_preconditioner = .true.

  ! number of coarse levels generated in the setup phase
  integer :: clevels = 0

  ! array to store factors required by LAPACK direct solver
  real(kind=myreal), dimension(:,:), allocatable :: lapack_factors

  ! array to store pivots required by LAPACK direct solver
  integer, dimension(:), allocatable :: lapack_pivots

  ! local copy of parameter with same name in mi20_control
  integer :: st_method

  ! true if data has been generated for LAPACK direct solver
  logical :: lapack_data = .false.

  ! stores how far generation of ma48 data has proceeded:
  ! set to 1 if ma48_initialise has been performed
  ! set to 2 if ma48_analyse has been performed
  ! set to 3 if ma48_factorize has been performed
  integer :: ma48_data = 0

  ! holds ma48 factors
  type(ma48_factors) :: ma48_factors

  ! ma48_control type
  type(ma48_control) :: ma48_cntrl

  ! holds matrix to be factorized by ma48
  type(zd11_type) :: ma48_matrix

  ! set to true if ma48_matrix exists
  logical :: ma48_matrix_exists = .false.

  ! if positive stores the level at which the last direct solve was performed
  integer :: dsolve_level = -1
  
  ! keep the max_its for the solve phase
  integer :: max_its = 0 

  ! set a flag for if conversion routine has been used
  logical :: zd11_internal_conversion = .false.

  ! the matrix A, for if mi20 does the conversion internally
  type( zd11_type ) :: A

 end type mi20_keep


 ! -----------------
 ! --  mi20_data  --
 ! -----------------

 ! Derived type holds coarse level data generated in mi20_setup

 type mi20_data

  ! sparse coarse level coefficient matrices (held as zd11 type)
  type(zd11_type) :: A_mat

   ! sparse interpolation/restriction matrices (held as zd11 type)
  type(zd11_type) :: I_mat

 end type mi20_data


 ! ------------------------
 ! -- mi20_solve_control --
 ! ------------------------
 
 ! Derived type holds controls for the solve phase
 
 type mi20_solve_control
    
    ! absolute convergence tolerance
    real(myreal):: abs_tol = 0
    
    ! breakdown tolerance (if using CG or biCGStab)
    real(myreal):: breakdown_tol = epsilon(1.0d0)
    
    ! restart (if using GMRES)
    integer :: gmres_restart = 100
    
    ! use initial guess?
    logical :: init_guess = .false.
    
    ! Which iterative method to use?
    ! set to 0 for pure AMG
    ! set to 1 for CG
    ! set to 2 for GMRES
    ! set to 3 for BiCGStab
    integer :: krylov_solver = 2
    
    ! max number of iterations (<0 means 2*matrix%m)
    integer :: max_its = -1
    
    ! which side preconditioning (GMRES) -- left or right?
    integer :: preconditioner_side = 1
    
    ! relative convergence tolerance
    real(myreal):: rel_tol = 1e-6
    
 end type mi20_solve_control

 ! -----------------
 ! --  list_node  --
 ! -----------------

 ! Derived type used to generate lists (of points with the same weight) and
 ! used in the first coarsening pass

 type list_node

  ! point number
  integer :: point = 0

  ! weight for this point
  integer :: weight = 0

  ! pointer to next point in the list
  type(list_node), pointer :: next => null()

  ! pointer to previous point in the list
  type(list_node), pointer :: last => null()

 end type list_node


 ! --------------------------
 ! --  precon_vec_storage  --
 ! --------------------------

 ! Derived type used to store vectors used in the preconditioner

 type precon_vec_storage

  ! vector to hold solution for level k, i.e. error at level k-1
  real(kind=myreal), dimension(:), allocatable :: vec_sol

  ! stores error or residual at level k
  real(kind=myreal), dimension(:), allocatable :: vec_er

  ! holds the right hand side vector at level k, i.e. residual at level k-1
  real(kind=myreal), dimension(:), allocatable :: vec_rhs

 end type precon_vec_storage

 ! == == == == == ==
 ! == INTERFACES  ==
 ! == == == == == ==

 interface mi20_precondition
    module procedure mi20_precondition_withA
    module procedure mi20_precondition_noA
 end interface

 interface mi20_solve
    module procedure mi20_solve_withA
    module procedure mi20_solve_noA
 end interface


contains

 ! == == == == == == == == == =
 ! ==  SUBROUTINES  ==
 ! == == == == == == == == == =

 ! ----------------
 ! -- mi20_setup --
 ! ----------------

 ! A subroutine to generate coarse level data required
 ! by the preconditioner.
 ! This also resets mi20_keep to default values

 subroutine mi20_setup(matrix, coarse_data, keep, control, info)

  ! matrix is the input matrix used to generate coarse level data.
  type(zd11_type), intent(inout) :: matrix

  ! derived type to hold coarse level data

  type(mi20_data), dimension(:), allocatable, intent(inout) :: coarse_data
! CHANGE BACK TO INTENT OUT (had to set to intent(inout) because of bug
! in ifort compiler.

  ! derived type containing mi20 internal values and arrays
  type(mi20_keep), intent(inout) :: keep

  ! derived type containing control parameters
  type(mi20_control), intent(in) :: control

  ! derived type to communicate information to the calling program
  type(mi20_info), intent(out) :: info

  ! Terminology
  ! -----------
  ! In multigrid it is usual to consider each matrix row as
  ! corresponding to a point, and to consider the connections between these
  ! points, i.e. point/row i is connected to point/row j if the matrix element
  ! at row i, column j is non-zero. Such connectivity is directional.
  ! During each coarsening step, we start with some matrix, and use this to
  ! produce a matrix at the coarse level.
  ! Method
  ! ------
  ! The steps to generate coarse level are:
  ! 1. Divide the points on the current level into C points (those which exist
  ! on the coarse level), F points (those which don't exist on the coarse
  ! level and so need to interpolate from the C points) and unconnected
  ! points. This requires three coarsening passes (after first removing
  ! unconnected points):
  ! Pass 1: starts by giving each point a weight, which is equal to the
  ! number of strong transpose connections. Then the following process is
  ! repeated until the maximum weight is zero: set an undecided point with
  ! maximum weight to be a C point, set all undecided points with a strong
  ! transpose connection to the new C to be F points, and increase the weight
  ! for undeciced points for each strong connection to a new F point.
  ! Pass 2: add C points to ensure that all F to F strongly connected pairs
  ! share a common C point.
  ! Pass 3: add C points to allow any remaining undecided points to become F.
  ! 2. Generate the interpolation/restriction matrix.
  ! 3. Generate the coarse level matrix

  ! storage for strong transpose connections
  type(zd11_type) :: str_t

  ! storage for temporary matrices during aggressive coarsening
  type(zd11_type) :: temp_matrix, temp_imatrix

  ! this array of type(list_node) is used to create linked lists
  type(list_node), allocatable, dimension(:), target :: lists

  ! this array of type(list_node) is used to point to the ends of each list
  type(list_node), allocatable, dimension(:), target :: lists_ends

  integer :: level            ! the current level being coarsened
  integer :: coarse_step      ! the current coarsening step
  integer :: size_matrix      ! size of matrix being coarsened
  integer :: size_cmatrix     ! size of new coarse matrix
  integer :: n_str            ! number of strong connections
  integer :: i, j             ! general purpose integers
  integer :: size_lists_ends  ! length of the array lists_ends
  integer :: max_weight       ! largest weight
  integer :: MC65_status      ! to get error status from MC65 subroutines
  integer :: old_matrix       ! position in coarse_data of last matrix
                              ! to be generated

  ! amount by which to increment the weight of undecided points for
  ! each strong connection to a F point
  integer :: weight_inc

  integer :: alloc_stat    ! used to find (de)allocation status
  integer :: copy_flag     ! copy of error flag

  ! array containing thresholds for each row
  real(kind=myreal), dimension(:), allocatable :: thresholds

  logical :: coarse_level_formed
  logical :: print_out    ! .true. if outputting diagnostics
  logical :: error_out    ! .true. if outputting errors

  ! point_type
  ! ----------
  ! array to store point type information as follows:
  !
  ! subroutine test_connections:
  ! sets point_type to -2 for any point with no connections and stores the
  ! number of strong transpose connections in point_type for all other points
  !
  ! subroutine generate_lists:
  ! sets point_type to 0 for any point with connections
  !
  ! subroutine coarsen_first_pass and beyond:
  ! point_type == -1 denotes a fine (F) point
  ! point_type == 0 denotes an undecided point
  ! point_type == -2 denotes a point with no connections
  ! point_type > 0 denotes a coarse (C) point
  !
  ! and before calling subroutine gen_interpolation the positive values in
  ! point_type denote the index of the C point on the next coarse level
  integer, allocatable, dimension(:) :: point_type


  ! reset keep%new_preconditioner
  ! -----------------------------
  keep%new_preconditioner = .true.

  ! reset other values before starting coarsening
  ! ---------------------------------------------
  coarse_level_formed = .false.
  info%flag = 0
  info%clevels = 0
  keep%st_method = control%st_method
  keep%clevels = -1  ! when set > 0 we know setup called successfully

  ! check control parameters
  ! ------------------------
  ! if they are not in range return to calling program
  call test_control_setup(control, info)
  if (info%flag < 0) return

  ! switch on/off output and error messages
  ! ---------------------------------------
  error_out = .false.
  if ( (control%print_level .gt. 0) .and. (control%error .ge. 0) ) &
    error_out = .true.

  print_out = .false.
  if ( (control%print_level .gt. 1) .and. (control%print .ge. 0) ) &
    print_out = .true.

  ! set m and type for matrix
  ! -------------------------
  matrix%n = matrix%m
  call zd11_put(matrix%type, "general", alloc_stat)
  if (alloc_stat .ne. 0) then
   info%stat = alloc_stat; info%flag = -10
   if (error_out) then
    write(control%error,'(a,i5)')" Error return from mi20_setup. info%flag = " &
                                   ,info%flag
    write(control%error,'(a)') &
     " failure to allocate type in zd11 matrix"
   end if
   return
  end if

  ! test matrix if required
  ! -----------------------
  ! return to calling program if error occurs
  if (control%testing == 1) then
   call test_matrix(matrix, control, info)
   if (info%flag < 0) return
  end if

  ! set initial values
  ! ------------------
  level = 0
  size_cmatrix = 0
  coarse_step = 0

  if (print_out) then
   write(control%print,'(/a)') " == == == == == == == == == == == == == == == ="
   write(control%print,'(a)') " mi20_setup called"
  end if

  ! allocate storage for coarse level data
  ! --------------------------------------
  if ( allocated(coarse_data) ) then
   info%flag = -10
   if (error_out) then
    write(control%error,'(a,i5)') &
   " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
     " user supplied mi20_data array has already been allocated"
   end if
   return
  end if
  allocate(coarse_data(control%max_levels), stat=alloc_stat)
  if (alloc_stat .ne. 0) then
   info%stat = alloc_stat; info%flag = -10
   if (error_out) then
    write(control%error,'(a,i5)') &
   " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
     " failure to allocate mi20_data array"
   end if
   return
  end if

  ! deal with the case when matrix%m == 1
  ! -------------------------------------
  if (matrix%m .eq. 1) then
   if (print_out) &
     write(control%print,'(a)') " mi20_setup completed. matrix%m = 1"
   keep%clevels = 0
   return
  end if

  ! one_coarsening_level
  ! --------------------
  ! attempts to generate a new coarse level (which is stored in level_data)
  one_coarsening_level: do

   if (print_out) then
    write(control%print,'(/a)') "-------------------------------"
    write(control%print,'(a,i5)') " coarsening level", level
    write(control%print,'(a)') "-------------------------------"
   end if

   ! one_coarsening_step
   ! -------------------
   ! performs one coarsening iteration
   coarse_step = 1
   one_coarsening_step: do

    ! find position of matrix to be coarsened
    if (coarse_step>1) then
     old_matrix = level+1
    else
     old_matrix = level
    end if

    ! output some info
    if ( print_out .and. (control%aggressive>1) ) &
     write(control%print,'(/a,i5)') " coarsening step", coarse_step

    ! find size of matrix to be coarsened
    if ( (level == 0) .and. (coarse_step == 1) ) then
     size_matrix = matrix%m
    else
     size_matrix = coarse_data(old_matrix)%A_mat%m
    end if

    ! Identify connection type between points (i.e. strong or weak)
    ! -------------------------------------------------------------

    ! Reminder of terminology
    ! -----------------------
    ! Point i is strongly coupled/connected to point j means
    ! the value at point i is strongly influenced by point j
    ! i.e. matrix_ij is large

    ! allocate some storage arrays
    ! ----------------------------

    ! array to store point type information
    allocate(point_type(size_matrix), stat=alloc_stat)
    if (alloc_stat .ne. 0) then
     info%stat = alloc_stat
     if (coarse_level_formed) then
      info%flag = 10
     else
      info%flag = -10
     end if
     exit one_coarsening_level
    end if
    point_type = 0

    ! array to hold threshold for each row
     allocate(thresholds(size_matrix), stat=alloc_stat)
     if (alloc_stat .ne. 0) then
      info%stat = alloc_stat
      if (coarse_level_formed) then
       info%flag = 10
      else
       info%flag = -10
      end if
      exit one_coarsening_level
    end if

    ! test conections are suitable for coarsening
    ! -------------------------------------------
    if ( (level == 0) .and. (coarse_step == 1) ) then
     call test_connections(matrix, control, &
       info, keep, point_type, thresholds, n_str)
     if (info%flag == 1) then
      if (error_out) &
        write(control%error,'(a,i5)')  &
       " Warning from mi20_setup. info%flag = " ,info%flag
     end if
    else
     call test_connections(coarse_data(old_matrix)%A_mat, control, &
       info, keep, point_type, thresholds, n_str)
    end if

    ! deal with failure
    ! -----------------
    if (info%flag < 0) then
     if (coarse_level_formed) info%flag = -info%flag
     exit one_coarsening_level
    end if

    ! find strong transpose connections (if required)
    ! -----------------------------------------------
    if (keep%st_method == 2) then
     ! find strong transpose connections
     if ( (level == 0) .and. (coarse_step == 1) ) then
      call find_strong_t(matrix, info, thresholds, point_type, &
                                n_str, str_t)
     else
      call find_strong_t(coarse_data(old_matrix)%A_mat, info, &
        thresholds, point_type, n_str, str_t)
     end if
    end if

    ! deal with errors
    if (info%flag < 0) then
     if (coarse_level_formed) info%flag = -info%flag
     exit one_coarsening_level
    end if

    ! Generate lists of points with equal weights
    ! -------------------------------------------
    ! these are linked lists, created using an array of type(list_nodes),
    ! each linked list corresponds to one value of weight
    !
    ! in addition an array of type(list_nodes) called lists_ends is used
    ! to point to the start and end of each list

    ! allocate lists
    allocate(lists(size_matrix), stat=alloc_stat)

    ! deal with failure
     if (alloc_stat .ne. 0) then
      info%stat = alloc_stat
      if (coarse_level_formed) then
       info%flag = 10
      else
       info%flag = -10
      end if
      exit one_coarsening_level
    end if

    ! allocate lists_ends - initial size is determined by average number
    ! of strong connections per matrix row
    size_lists_ends = 2*int( real(n_str)/size_matrix + 1 )
    allocate(lists_ends(size_lists_ends), stat=alloc_stat)

    ! deal with error
    if (alloc_stat .ne. 0) then
     info%stat = alloc_stat
     if (coarse_level_formed) then
       info%flag = 10
      else
      info%flag = -10
     end if
      exit one_coarsening_level
    end if

    ! set up lists
    call generate_lists(point_type, info, lists, lists_ends, max_weight)

    ! deal with errors
    if (info%flag < 0) then
     if (coarse_level_formed) info%flag = -info%flag
     exit one_coarsening_level
    end if

     ! Find C & F points - first pass
    ! ------------------------------
    weight_inc = 2
    if ( (level == 0) .and. (coarse_step == 1) ) then
     call coarsen_first_pass(matrix, str_t, &
       thresholds, weight_inc, lists, lists_ends, point_type, &
       max_weight, size_cmatrix, info, keep)
    else
      call coarsen_first_pass(coarse_data(old_matrix)%A_mat, str_t, &
       thresholds, weight_inc, lists, lists_ends, point_type, &
       max_weight, size_cmatrix, info, keep)
    end if

    ! deal with errors
    if (info%flag < 0) then
     if (coarse_level_formed) info%flag = -info%flag
     exit one_coarsening_level
    end if

    ! Deallocate lists_ends
    deallocate(lists_ends, stat=alloc_stat)
    if (alloc_stat .ne. 0) then
     info%stat = alloc_stat
     if (coarse_level_formed) then
       info%flag = 11
     else
       info%flag = -11
     end if
     exit one_coarsening_level
    end if


    ! Bypass 2nd and 3rd pass of coarsening if requested
    ! --------------------------------------------------
    if(.not.(control%one_pass_coarsen)) then

     ! Find C & F points - second pass
     ! -------------------------------
     if ( (level == 0) .and. (coarse_step == 1) ) then
       call coarsen_second_pass(matrix, thresholds, &
          str_t, keep, size_cmatrix, point_type)
      else
       call coarsen_second_pass(coarse_data(old_matrix)%A_mat, thresholds, &
          str_t, keep, size_cmatrix, point_type)
      end if

     ! Find C & F points - third pass
     ! -------------------------------
     ! can bypass if no points are undecided
     if ( any(point_type(1:size_matrix) == 0) ) then
      if ( (level == 0) .and. (coarse_step == 1) ) then
       call coarsen_third_pass(matrix, &
          str_t, thresholds, keep, size_cmatrix, point_type)
      else
       call coarsen_third_pass(coarse_data(old_matrix)%A_mat, &
          str_t, thresholds, keep, size_cmatrix, point_type)
      end if
     end if

    ! end of 2nd and 3rd pass coarsening
    ! ----------------------------------
    end if

    ! If the 2nd and 3rd pass of coarsening has not been bypassed
    ! at this point any remaining undecided points must be those which
    ! have no strong transpose connections, and only strong connections
    ! to other undecided points.

    ! Test for coarsening stagnation
    ! ------------------------------
    if (real(size_cmatrix) >= control%reduction*real(size_matrix)) then
     if (coarse_level_formed) then
       info%flag = 13
       exit one_coarsening_level
     end if
    end if

    ! find position of the C points on next level
    ! -------------------------------------------
    i = 0;
    do j = 1, size_matrix
     if ( point_type(j) .gt. 0 ) then
      i = i + 1
      point_type(j) = i
     end if
    end do

    ! coarsening has completed successfully
    ! -------------------------------------
    ! output some info
    if (print_out) then
     write(control%print,'(a)') &
       " Coarsening has completed successfully"
     write(control%print,'(a,i8)') &
       " Number of C points", count(point_type(1:size_matrix) .gt. 0)
     write(control%print,'(a,i8)') &
       " Number of F points", count(point_type(1:size_matrix) == -1)
     write(control%print,'(a,i8)') &
       " Number of unconnected points", &
       count(point_type(1:size_matrix) == -2)
     write(control%print,'(a,i8)') &
       " Number of undecided points  ", &
       count(point_type(1:size_matrix) == 0)
    end if

    ! Find interpolation matrix
    ! -------------------------
    ! call subroutine to generate interpolation matrix
    if (coarse_step>1)then
     ! after coarse step 1 put new interpolation matrix in temp_imatrix
     call gen_interpolation(coarse_data(level+1)%A_mat, &
         point_type, thresholds, temp_imatrix, control, info)
    else
     if (level == 0) then
      call gen_interpolation(matrix, point_type, &
         thresholds, coarse_data(level+1)%I_mat, control, info)
     else
      call gen_interpolation(coarse_data(level)%A_mat, &
         point_type, thresholds, coarse_data(level+1)%I_mat, control, info)
     end if
    end if

    ! deal with failure
    if (info%flag < 0) then
     if (coarse_level_formed) info%flag = -info%flag
     exit one_coarsening_level
    end if

     ! Find coarse level matrix
     ! ------------------------
     if (coarse_step>1) then
      ! after coarse step 1 put new coarse matrix in temp_matrix
      call gen_coarse_matrix(coarse_data(level+1)%A_mat, &
        temp_imatrix, temp_matrix, info)
     else
     if (level == 0) then
      call gen_coarse_matrix(matrix, &
        coarse_data(level+1)%I_mat, coarse_data(level+1)%A_mat, info)
     else
      call gen_coarse_matrix(coarse_data(level)%A_mat, &
        coarse_data(level+1)%I_mat, coarse_data(level+1)%A_mat, info)
     end if
    end if

    ! deal with failure
    if (info%flag < 0) then
     if (coarse_level_formed) info%flag = -info%flag
     exit one_coarsening_level
    end if

    ! Deal with aggressive coarsening
    ! -------------------------------
    if (coarse_step>1) then
     ! delete (level+1)
     call MC65_matrix_destruct(coarse_data(level+1)%A_mat, MC65_status, &
          stat=info%stat)
     ! deal with any MC65 error
     if (MC65_status == MC65_ERR_MEMORY_DEALLOC) then
      ! deal with the case of failure on level 1, i.e. where we may have
      ! wiped the matrix formed on the previous coarsening step, so:
      ! if level == 1 we are not able to form any coarse levels
      ! and reset coarse step to 1
      if (level == 1) then
        coarse_level_formed = .false.
      end if
      coarse_step = 1
      if (coarse_level_formed) then
        info%flag = 11
      else
        info%flag = -11
      end if
      exit one_coarsening_level
     end if

     call MC65_matrix_destruct(coarse_data(level+1)%I_mat, MC65_status, &
          stat=info%stat)
     ! deal with any MC65 error
     if (MC65_status == MC65_ERR_MEMORY_DEALLOC) then
      ! deal with the case of failure on level 1, i.e. where we may have
      ! wiped the matrix formed on the previous coarsening step:
      ! if we are on level 1, there are no coarse matrices
      if (level == 1) coarse_level_formed = .false.
      coarse_step = 1
      if (coarse_level_formed) then
        info%flag = 11
      else
        info%flag = -11
      end if
      exit one_coarsening_level
     end if

     ! copy temp matrices to (level+1)
     temp_matrix%ne  = -1 ! Need to define %ne to keep nagfor -C=undefined happy
     temp_imatrix%ne = -1 ! Need to define %ne to keep nagfor -C=undefined happy
     coarse_data(level+1)%A_mat = temp_matrix
     coarse_data(level+1)%I_mat = temp_imatrix

     ! deallocate temporary matrices
     call MC65_matrix_destruct(temp_matrix, MC65_status, stat=info%stat)
     ! deal with any MC65 error
     if (MC65_status == MC65_ERR_MEMORY_DEALLOC) then
       if (coarse_level_formed) then
         info%flag = 11
       else
         info%flag = -11
       end if
       exit one_coarsening_level
     end if

     call MC65_matrix_destruct(temp_imatrix, MC65_status, stat=info%stat)
     ! deal with any MC65 error
     if (MC65_status == MC65_ERR_MEMORY_DEALLOC) then
       if (coarse_level_formed) then
         info%flag = 11
       else
         info%flag = -11
       end if
       exit one_coarsening_level
     end if

    end if ! end of dealing with aggressive coarsening

    ! generation of coarse level matrices has completed successfully
    ! --------------------------------------------------------------
    ! output some info
    if (print_out) then
     write(control%print,'(a)') &
       " Coarse level matrices generated successfully"
    end if
    coarse_level_formed = .true.

    ! deallocation
    call setup_deallocation(info, point_type, thresholds, &
                         lists, lists_ends, str_t)

    ! deal with failure
    if (info%flag < 0) then
     if (coarse_level_formed) info%flag = -info%flag
     exit one_coarsening_level
    end if

    ! stop coarsening at this level if there is only one point
    if (size_cmatrix == 1) exit one_coarsening_step

    ! stop coarsening if we have performed the right number of steps
    if (coarse_step == control%aggressive) exit one_coarsening_step

    coarse_step = coarse_step + 1
   end do one_coarsening_step

   ! coarsening success for this level
   ! ---------------------------------
   level = level + 1

   ! output info on coarsening
   if (print_out) then
    write(control%print,'(a,i4,a,i8,a)') " Level", level, " generated with", &
    coarse_data(level)%I_mat%n, " points"
   end if

   ! set values in info
   info%clevels = level
   keep%clevels = level
   info%cnnz = size(coarse_data(level)%A_mat%val)
   info%cpoints = coarse_data(level)%A_mat%m

   ! test for completion of coarsening
   if ( (size_cmatrix <= control%max_points) &
        .or. (level == control%max_levels) ) exit

  end do one_coarsening_level

  ! deal with coarsening problems
  ! -----------------------------
  if ( info%flag .ne. 0 ) then

   ! output message
   if (error_out) &
    call error_message_setup(control, info, level, coarse_step)

   ! deal with coarsening failure
   ! ----------------------------
   if (info%flag .ne. 1) then

    ! We do not want to overwrite the error flag that has been
    ! set so take a copy (we are not interested in any deallocation
    ! errors that occur at this point as they will obscure the error
    ! that has already happened)

    copy_flag = info%flag

    ! attempt to deallocate memory
    call setup_deallocation(info, point_type, thresholds, &
         lists, lists_ends, str_t)
    info%flag = copy_flag

    ! If this is coarsening step 1 (ie non aggressive coarsening) we can't
    ! form the next coarse level. So deallocate data stored in
    ! coarse_data(level+1).
    ! If this is level 0 we can deallocate coarse_data too.

    if (coarse_step == 1) then
     call MC65_matrix_destruct(coarse_data(level+1)%I_mat, MC65_status)

     call MC65_matrix_destruct(coarse_data(level+1)%A_mat, MC65_status)

     if (level == 0) then
      deallocate(coarse_data)
      ! return to the user
      return
     end if
    else
     ! for aggressive coarsening if coarse_step > 1 we can use the coarse
     ! level already created. So set necessary fields in info
!!$     if (coarse_step>1) then
!!$      ! after coarse step 1 put new coarse matrix in temp_matrix
!!$      call gen_coarse_matrix(coarse_data(level+1)%A_mat, &
!!$        temp_imatrix, temp_matrix, info)
!!$     else
!!$     if (level == 0) then
!!$      call gen_coarse_matrix(matrix, &
!!$        coarse_data(level+1)%I_mat, coarse_data(level+1)%A_mat, info)
!!$     else
!!$      call gen_coarse_matrix(coarse_data(level)%A_mat, &
!!$        coarse_data(level+1)%I_mat, coarse_data(level+1)%A_mat, info)
!!$     end if
!!$
     level = level + 1
     info%clevels = level
     keep%clevels = level
     info%cnnz = size(coarse_data(level)%A_mat%val)
     info%cpoints = coarse_data(level)%A_mat%m
    end if
   end if ! end of dealing with coarsening errors
  end if ! end of dealing with coarsening problems

  ! set diagonals first in last A_mat generated so Gauss-Seidel and Jacobi
  ! direct solvers work on this level
  call MC65_matrix_diagonal_first(coarse_data(level)%A_mat, MC65_status)

  if (print_out) then
   write(control%print,'(/a)') " coarsening finished"
   write(control%print,'(a)') &
     " == == == == == == == == == == == == == == == == ="
  end if

 end subroutine mi20_setup

 ! ---------------------------
 ! --     mi20_setup_csr    --
 ! -- (T. Rees, April 2015) --
 ! ---------------------------

 subroutine mi20_setup_csr(ptr, col, val, ne, n, & 
      coarse_data, keep, control, info)

  use hsl_mc69_double

     ! matrix in
   integer, intent(in) :: col(:), ptr(:)
   integer, intent(in) :: ne, n
   real(myreal), intent(in) :: val(:)

  ! derived type to hold coarse level data
  type(mi20_data), dimension(:), allocatable, intent(inout) :: coarse_data
! CHANGE BACK TO INTENT OUT (had to set to intent(inout) because of bug
! in ifort compiler.

  ! derived type containing mi20 internal values and arrays
  type(mi20_keep), intent(inout) :: keep

  ! derived type containing control parameters
  type(mi20_control), intent(in) :: control

  ! derived type to communicate information to the calling program
  type(mi20_info), intent(out) :: info

  integer :: matrix_type = 2 ! real, unsymmetric
  integer :: flag 

  logical :: error_out
  
  error_out = .false.
  if ( (control%print_level .gt. 0) .and. (control%error .ge. 0) ) &
       error_out = .true.
  
  ! convert to csc format using hsl_mc69
  allocate(keep%A%ptr(n+1))
  allocate(keep%A%val(ne),keep%A%col(ne))
  call mc69_cscl_convert(matrix_type, &
       n, n, ptr, col, &
       keep%A%ptr, keep%A%col, &
       flag, & 
       val, keep%A%val)
  if (flag < 0) goto 100
  if (flag > 0) info%flag = 0
  
  keep%A%m = n
  keep%A%n = n
  keep%A%ne = ne
  
  call mi20_setup(keep%A,coarse_data,keep,control,info)
  
  keep%zd11_internal_conversion = .true.

  RETURN
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! E r r o r   M e s s a g e s !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
100 continue
  select case (flag)
  case (-1)
     if (error_out) then
        write(control%error,'(a)') 'Allocation error in mi20_setup_csr'
     end if
     info%flag = -301
  case (-3)
    if (error_out) then 
       write(control%error,'(a)') 'm supplied in mi20_setup_csr < 1'
    end if
    info%flag = -302
 case (-10)
    if (error_out) then
       write(control%error,'(a)') 'all entries supplied are out-of-range'
    end if
    info%flag = -303
 case (-5)
    if (error_out) then 
       write(control%error,'(a)') 'ptr(1) < 1 in array passed to mi20_setup_csr'
    end if
    info%flag = -304
 case (-6)
    if (error_out) then
       write(control%error,'(a)') 'entries of ptr array passed to ', &
            'mi20_setup_csr not monotonic increasing'
    end if
    info%flag = -305
 end select

 ! deallocate if needed...
 if (allocated(keep%A%col)) deallocate(keep%A%col)
 if (allocated(keep%A%val)) deallocate(keep%A%val)
 if (allocated(keep%A%ptr)) deallocate(keep%A%ptr)

 return
    
 end subroutine mi20_setup_csr

 ! ---------------------------
 ! --    mi20_setup_csc     --
 ! -- (T. Rees, April 2015) --
 ! ---------------------------

 subroutine mi20_setup_csc(ptr, row, val, ne, n, & 
      coarse_data, keep, control, info)
   
   use hsl_mc69_double

   ! matrix in
   integer, intent(in) :: row(:), ptr(:)
   integer, intent(in) :: ne, n
   real(myreal), intent(in) :: val(:)

  ! derived type to hold coarse level data
  type(mi20_data), dimension(:), allocatable, intent(inout) :: coarse_data
! CHANGE BACK TO INTENT OUT (had to set to intent(inout) because of bug
! in ifort compiler.

  ! derived type containing mi20 internal values and arrays
  type(mi20_keep), intent(inout) :: keep

  ! derived type containing control parameters
  type(mi20_control), intent(in) :: control

  ! derived type to communicate information to the calling program
  type(mi20_info), intent(out) :: info

  integer :: matrix_type = 2 ! real, unsymmetric
  integer :: flag 

  logical :: error_out



      error_out = .false.
    if ( (control%print_level .gt. 0) .and. (control%error .ge. 0) ) &
         error_out = .true.
    
    ! convert to csc format using hsl_mc69
    allocate(keep%A%ptr(n+1))
    allocate(keep%A%val(ne),keep%A%col(ne))
    call mc69_csrl_convert(matrix_type, &
         n, n, ptr, row, &
         keep%A%ptr, keep%A%col, &
         flag, & 
         val, keep%A%val)
    if (flag < 0) goto 100
    if (flag > 0) info%flag = 0
    
    keep%A%m = n
    keep%A%n = n
    keep%A%ne = ne
    
    call mi20_setup(keep%A,coarse_data,keep,control,info)

    keep%zd11_internal_conversion = .true.

    RETURN
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! E r r o r   M e s s a g e s !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
100 continue
    select case (flag)
    case (-1)
       if (error_out) then
          write(control%error,'(a)') 'Allocation error in mi20_setup_csc'
       end if
       info%flag = -301
    case (-3)
       if (error_out) then 
          write(control%error,'(a)') 'm supplied in mi20_setup_csc < 1'
       end if
       info%flag = -302
    case (-10)
       if (error_out) then
          write(control%error,'(a)') 'all entries supplied are out-of-range'
       end if
       info%flag = -303
    case (-5)
       if (error_out) then 
          write(control%error,'(a)') 'ptr(1) < 1 in array passed', & 
                                     'to mi20_setup_csc'
       end if
       info%flag = -304
    case (-6)
       if (error_out) then
          write(control%error,'(a)') 'entries of ptr array passed to ', &
               'mi20_setup_csc not monotonic increasing'
       end if
       info%flag = -305
    end select

    ! deallocate if needed...
    if (allocated(keep%A%col)) deallocate(keep%A%col)
    if (allocated(keep%A%val)) deallocate(keep%A%val)
    if (allocated(keep%A%ptr)) deallocate(keep%A%ptr)
     
    return



 end subroutine mi20_setup_csc

 ! ---------------------------
 ! --    mi20_setup_coord   -- 
 ! -- (T. Rees, April 2015) --
 ! ---------------------------

 subroutine mi20_setup_coord(row, col, val, ne, n, & 
      coarse_data, keep, control, info)

   use hsl_mc69_double

   ! matrix in
   integer, intent(in) :: row(:), col(:)
   integer, intent(in) :: ne, n
   real(myreal), intent(in) :: val(:)

  ! derived type to hold coarse level data
  type(mi20_data), dimension(:), allocatable, intent(inout) :: coarse_data
! CHANGE BACK TO INTENT OUT (had to set to intent(inout) because of bug
! in ifort compiler.

  ! derived type containing mi20 internal values and arrays
  type(mi20_keep), intent(inout) :: keep

  ! derived type containing control parameters
  type(mi20_control), intent(in) :: control

  ! derived type to communicate information to the calling program
  type(mi20_info), intent(out) :: info
  
  integer :: matrix_type = 2 ! real, unsymmetric
  integer :: flag
  logical :: error_out

   error_out = .false.
    if ( (control%print_level .gt. 0) .and. (control%error .ge. 0) ) &
         error_out = .true.
    
    ! convert to csc format using hsl_mc69
    allocate(keep%A%ptr(n+1))
    allocate(keep%A%val(ne),keep%A%col(ne))
    call mc69_coord_convert(matrix_type, &
         n, n, ne, & 
         col, row, & ! Swap row and col, as mc69 gives CSC
         keep%A%ptr, keep%A%col, &
         flag, & 
         val, keep%A%val)
    if (flag < 0) goto 100
    if (flag > 0) info%flag = 0

    keep%A%m = n
    keep%A%n = n
    keep%A%ne = ne
            
    call mi20_setup(keep%A, coarse_data,keep,control,info)

    keep%zd11_internal_conversion = .true.

    RETURN
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! E r r o r   M e s s a g e s !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
100 continue
    select case (flag)
    case (-1)
       if (error_out) then 
          write(control%error,'(a)') 'Allocation error in mi20_setup_coord'
       end if
       info%flag = -301
    case (-3)
       if (error_out) then 
          write(control%error,'(a)') 'm supplied in mi20_setup_coord < 1'
       end if
       info%flag = -302
    case (-10)
       if (error_out) then 
          write(control%error,'(a)') 'all entries supplied are out-of-range'
       end if
       info%flag = -303
    end select

    ! deallocate if needed...
    if (allocated(keep%A%col)) deallocate(keep%A%col)
    if (allocated(keep%A%val)) deallocate(keep%A%val)
    if (allocated(keep%A%ptr)) deallocate(keep%A%ptr)

    return
   
 end subroutine mi20_setup_coord
 

 ! ---------------------------
 ! --   test_control_setup  --
 ! ---------------------------

 ! A subroutine to test that setup control paramteres are in range.
 ! If they are not info%flag is set to a negative value and an
 ! appropriate message to sent to the output stream defined by
 ! control%error

 subroutine test_control_setup(control, info)

  type(mi20_control), intent(in) :: control
  type(mi20_info), intent(inout) :: info

  logical :: error_out ! .true. if error messages are to be issued

  ! switch on/off error messages
  ! ----------------------------
  if ( (control%print_level .gt. 0) .and. (control%error .ge. 0) )then
   error_out = .true.
  else
   error_out = .false.
  end if

  ! perform testing
  ! ---------------
  if ( (control%trunc_parameter .lt. zero) &
   .or. (control%trunc_parameter .ge. one) ) then
   info%flag = -118
   if (error_out) then
    write(control%error,'(a,i5)')" Error return from mi20_setup. info%flag = "&
                                   ,info%flag
    write(control%error,'(a)') &
      " mi20_control error: trunc_parameter must lie in the range 0.0 to 1.0"
   end if
  end if

  if ( (control%print_level<0) .or. (control%print_level>2) ) then
   info%flag = -114
   if (control%error .ge. 0) then
    write(control%error,'(a,i5)')&
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
      " mi20_control error: print_level must lie in the range 0 to 2"
   end if
  end if

  if (control%aggressive .lt. 1) then
   info%flag = -105
   if (error_out) then
    write(control%error,'(a,i5)')&
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
      " mi20_control error: aggressive must be greater than 0"
   end if
  end if

  if ( (control%c_fail<1) .or. (control%c_fail>2) ) then
   info%flag = -106
   if (error_out) then
    write(control%error,'(a,i5)')&
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
      " mi20_control error: c_fail must lie in the range 1 to 2"
   end if

  end if

  if (control%max_levels<1) then
   info%flag = -116
   if (error_out) then
    write(control%error,'(a,i5)')&
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
      " mi20_control error: max_levels must be greater than 0"
   end if

  end if

  if (control%reduction < 0.5_myreal .or. control%reduction > 1.0_myreal) then
   info%flag = -119
   if (error_out) then
    write(control%error,'(a,i5)')&
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
      " mi20_control error: reduction must lie between 0.5 and 1.0"
   end if

  end if

  if (control%max_points<1) then
   info%flag = -103
   if (error_out) then
    write(control%error,'(a,i5)')&
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
     " mi20_control error: max_points must be greater than 0"
   end if
  end if

  if ( (control%st_method<1) .or. (control%st_method>2) ) then
   info%flag = -104
   if (error_out) then
    write(control%error,'(a,i5)')&
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
      " mi20_control error: st_method must lie in the range 1 to 2"
   end if
  end if

  if ( (control%st_parameter<zero) .or. (control%st_parameter>one) ) then
   info%flag = -101
   if (error_out) then
    write(control%error,'(a,i5)')&
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
      " mi20_control error: st_parameter must lie between 0.0 & 1.0"
   end if
  end if

  if ( (control%testing<0) .or. (control%testing>1) ) then
   info%flag = -100
   if (error_out) then
    write(control%error,'(a,i5)')&
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
     " mi20_control error: testing must lie in the range 0 to 1"
   end if

  end if

 end subroutine test_control_setup


 ! -------------------------------
 ! -- test_control_precondition --
 ! -------------------------------

 ! A subroutine to test that precondition control parameters are in range.
 ! If they are not info%flag is set to a negative value and an
 ! appropriate message to sent to the output stream defined by
 ! control%error

subroutine test_control_precondition(control, info)

  type(mi20_control), intent(in) :: control
  type(mi20_info), intent(inout) :: info

  logical :: error_out ! .true. if error messages are to be issued

  ! switch on/off error messages
  ! ----------------------------
  if ( (control%print_level .gt. 0) .and. (control%error .ge. 0) )then
   error_out = .true.
  else
   error_out = .false.
  end if

  ! perform testing
  ! ---------------

  if ( (control%print_level<0) .or. (control%print_level>2) ) then
   info%flag = -114
   if (control%error .ge. 0) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_precondition. info%flag = " ,info%flag
    write(control%error,'(a)') &
      " mi20_control error: print_level must lie in the range 0 to 2"
   end if
  end if

  if ( (control%coarse_solver<1) .or. (control%coarse_solver>4) ) then
   info%flag = -112
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_precondition. info%flag = " ,info%flag
    write(control%error,'(a)') &
      " mi20_control error: coarse_solver must lie in the range 1 to 4"
   end if
  end if

  if ( (control%coarse_solver<3) .and. (control%coarse_solver_its<1) ) then
   info%flag = -113
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_precondition. info%flag = " ,info%flag
    write(control%error,'(a)') &
    " mi20_control error: coarse_solver_its must be greater than 0"
   end if
  end if

  if ( (control%smoother == 1) .and. &
   ( (control%damping < zero) .or. (control%damping > one) ) ) then
   info%flag = -115
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_precondition. info%flag = " ,info%flag
    write(control%error,'(a)') &
     " mi20_control error: damping must lie between 0.0 and 1.0"
   end if
  end if

  if ( control%err_tol .le. zero ) then
   info%flag = -102
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_precondition. info%flag = " ,info%flag
    write(control%error,'(a)') &
     " mi20_control error: err_tol must be greater than 0.0"
   end if
  end if

  if ( (control%coarse_solver .eq. 3) .and. &
    (control%ma48 .gt. 4) .or. (control%ma48 .lt. 0) ) then
   info%flag = -117
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_precondition. info%flag = " ,info%flag
    write(control%error,'(a)') &
     " mi20_control error: ma48 must lie in the range 0 to 4"
   end if
  end if

  if ( (control%smoother<1) .or. (control%smoother>2) ) then
   info%flag = -108
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_precondition. info%flag = " ,info%flag
    write(control%error,'(a)') &
     " mi20_control error: smoother must lie in the range 1 to 2"
   end if
  end if

  if (control%pre_smoothing<0) then
   info%flag = -109
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_precondition. info%flag = " ,info%flag
    write(control%error,'(a)') &
     " mi20_control error: pre_smoothing can not be negative"
   end if
  end if

  if ( control%post_smoothing<0 ) then
   info%flag = -110
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_precondition. info%flag = " ,info%flag
    write(control%error,'(a)') &
     " mi20_control error: post_smoothing can not be negative"
   end if
  end if

  if (control%pre_smoothing+control%post_smoothing == 0) then
   info%flag = -111
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_precondition. info%flag = " ,info%flag
    write(control%error,'(a)') &
     " mi20_control error: no smoothing set"
   end if
  end if


  if (control%v_iterations<1) then
   info%flag = -107
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_precondition. info%flag = " ,info%flag
    write(control%error,'(a)') &
     " mi20_control error: v_iterations should be greater than 0"
   end if
  end if

 end subroutine test_control_precondition


 ! -----------------
 ! -- test_matrix --
 ! -----------------

 ! A subroutine to test the matrix conforms to certain properties
 ! 1. storage arrays are allocated
 ! 2. all values in matrix%col .lt. matrix%m .and. .gt. 0
 ! 3. matrix%ptr, matrix%val and matrix%col are long enough
 ! 4. no duplicate column entries in a row
 ! 5. the diagonal exists and is positive
 ! 6. matrix%ptr(1) == 1 .and. matrix%ptr(i+1)>ptr(i)
 ! 7. matrix%ptr(m+1)-1 .le. size(col) .and. size(val)

subroutine test_matrix(matrix, control, info)

  type(zd11_type), intent(in) :: matrix  ! zd11 object holding matrix to test
  type(mi20_control), intent(in) :: control
  type(mi20_info), intent(inout) :: info

  integer :: row, row_start, row_end, col
  integer :: i
  integer :: nnz ! number of non zero entries in the matrix
  integer :: st

  integer, allocatable :: iw(:)

  logical :: diag_missing
  logical :: error_out ! .true. if error messages are to be issued

  ! switch on/off error messages
  ! ----------------------------
  error_out = .false.
  if ( (control%print_level .gt. 0) .and. (control%error .ge. 0) ) &
   error_out = .true.

  ! check matrix%m
  ! --------------
  if (matrix%m .lt. 1) then
   info%flag = -9
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a,i12)') &
    " zd11_type input argument has component m set as ", matrix%m
   end if
   return
  end if

  ! check size of ptr
  ! -----------------
  if (.not. allocated(matrix%ptr)) then
   info%flag = -4
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
    " zd11_type input argument does not have an allocated ptr component"
   end if
   return
  end if

  if (size(matrix%ptr) .lt. matrix%m+1) then
   info%flag = -4
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
    " zd11_type input argument has component ptr which is too small"
   end if
   return
  end if

  ! set number of non zeros
  nnz = matrix%ptr(1+matrix%m)-1

  ! check size of val
  ! -----------------
  if (.not. allocated(matrix%val)) then
   info%flag = -6
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
    " zd11_type input argument does not have an allocated val component"
   end if
   return
  end if

  if (size(matrix%val) .lt. nnz) then
   info%flag = -6
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
    " zd11_type input argument has component val which is too small"
   end if
   return
  end if

  ! check matrix%col
  ! ----------------
  if (.not. allocated(matrix%col)) then
   info%flag = -5
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
    " zd11_type input argument does not have an allocated col component"
   end if
   return
  end if

  if (size(matrix%col) .lt. nnz) then
   info%flag = -5
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
    " zd11_type input argument has component col which is too small"
   end if
   return
  end if

  if ( any(matrix%col(1:nnz) .lt. 1)) then
   info%flag = -1
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') &
    " zd11_type input argument has &
               &component col which contains values less than 1"
   end if
   return
  end if

  if ( any(matrix%col(1:nnz) .gt. matrix%m)) then
   info%flag = -1
   if (error_out) then
    write(control%error,'(a,i5)') &
    " Error return from mi20_setup. info%flag = ",info%flag
    write(control%error,'(a)') " zd11_type input argument has &
         &component col containing values greater", &
         " than number of columns in the matrix"
   end if
   return
  end if

  ! run through matrix and do remainder of testing
  ! ----------------------------------------------
  ! first check first entry in matrix%ptr
  if ( matrix%ptr(1) .ne. 1 ) then
   info%flag = -7
  ! otherwise run through matrix rows
  else

! array iw will be used to test for duplicates
   deallocate (iw,stat=st)
   allocate (iw(1:matrix%m),stat=st)
   if (st /= 0) then
     info%stat = st; info%flag = -10
     if (error_out) write(control%error,'(a,i5)') &
    " Error return from mi20_setup. info%flag = ",info%flag
     return
   end if
   iw(1:matrix%m) = 0
   row_loop: do row = 1, matrix%m

    row_start = matrix%ptr(row)
    row_end = matrix%ptr(row+1)-1

    if (row_start .gt. row_end) then
     info%flag = -7
     exit
    end if

    ! run through the row
    diag_missing = .true.

    do i = row_start, row_end
     col = matrix%col(i)

     ! test diagonal
     if (col .eq. row) then
      diag_missing = .false.

      ! test for positive entry
      if (matrix%val(i) .le. zero ) then
       info%flag = -3
       exit row_loop
      end if
     end if

     ! check for duplicate entries in row
     if ( iw(col) .eq. row ) then
      info%flag = -8
      exit row_loop
     end if

     iw(col) = row
    end do

    if (diag_missing) then
     info%flag = -2
     exit row_loop
    end if

   end do row_loop
  end if
  deallocate (iw,stat=st)

  ! output errors from running through matrix
  if (error_out) then
   if (info%flag < 0)  write(control%error,'(a,i5)') &
      " Error return from mi20_setup. info%flag = ", info%flag
   if (info%flag .eq. -2) then
     write(control%error,'(a)') " zd11_type input argument has &
            &missing diagonal entries"
   else if (info%flag .eq. -7) then
     write(control%error,'(a)') " zd11_type input argument has &
            &component ptr with out of range entries"
   else if (info%flag .eq. -3) then
     write(control%error,'(a)') " zd11_type input argument has &
            &non positive diagonal entries"
   else if (info%flag .eq. -8) then
     write(control%error,'(a)') " zd11_type input argument has &
            &component col with duplicate columns ", &
            "in the same row"
   end if
  end if

 end subroutine test_matrix

 ! ----------------------
 ! -- test_connections --
 ! ----------------------

 ! A subroutine to test connections, this:
 !
 ! 1. Tests for inability to coarsen this level. Two
 ! coarsening problems may occur:
 ! 1) Rows occur which have only positive off diagonals. Such
 !    a row has connections but no strong connections since
 !    only negative connections can be strong. These are
 !    treated as unconnected and take no further part in
 !    coarsening.
 ! 2) Rows occur which have only strong connections
 !    to rows with no connections, i.e. the connection is to
 !    a row with no negative off diagonals. These are
 !    treated as unconnected and take no further part
 !    in the coarsening
 !
 ! Two failure criteria for coarsening are defined:
 ! a) if (control%c_fail == 1) then failure occurs if ANY row has
 !    only positive off diagonals.
 ! b) if (control%c_fail == 2) then failure occurs if ALL rows have
 !    only positive off diagonals.
 ! In addition failure may occur if all connections are to unconnected
 ! points (as a result of the connection being directional).
 ! If conditions for failure are met info%flag is set to -12
 !
 ! 2. if (keep%st_method == 1) also tests if matrix pattern is symmetric
 ! if it is not then keep%st_method is changed to 2 and
 ! info%flag is set to 1
 !
 ! This subroutine also:
 ! 1. moves diagonals to the first position within the row
 ! 2. for a point with no off-diagonals sets the point_type to -2
 ! 3. stores the number of strong transpose connections in
 !    for all other points
 ! 4. finds the threshold for each row

 subroutine test_connections(matrix, control, info, keep, point_type, &
                           thresholds, n_str)

  ! matrix to be tested
  type(zd11_type), intent(inout) :: matrix
  type(mi20_control), intent(in) :: control
  type(mi20_info), intent(inout) :: info
  type(mi20_keep), intent(inout) :: keep
  integer, intent(inout), dimension(:) :: point_type
  real(kind=myreal), intent(inout), dimension(:) :: thresholds
  integer, intent(out) :: n_str  ! number of strong connections

  integer :: MC65_status         ! error flag for MC65 subroutines
  integer :: row_start, row_end  ! start and end of a row
  integer :: row                 ! current row
  integer :: i                   ! general purpose integer

  logical :: pattern       ! used in call to MC65_pattern_is_symmetric
  logical :: symmetric     ! used in call to MC65_pattern_is_symmetric

  pattern = .true.

  ! set diagonals so they are first in each row
  call MC65_matrix_diagonal_first(matrix, MC65_status)

  n_str = 0

  ! run through matrix rows
  row_loop: do row = 1, matrix%m

   ! get location of columns for this row
   row_start = matrix%ptr(row)
   row_end = matrix%ptr(row+1)-1

   ! check first entry is the diagonal
   if (row .ne. matrix%col(row_start) ) then
    info%flag = -2
    return
   end if

   ! find rows with only one entry and set as unconnected
   if ( row_start == row_end ) then
    ! subtract from n_str any strong connections
    ! already included from this point
    n_str = n_str - point_type(row)
    !  set point type to -2 (i.e. an unconnected point)
    point_type(row) = -2
    thresholds(row) = 0
   else ! this row has off diagonals
    ! find threshold for this row and deal with the case of
    ! rows with only positive off diagonals
    thresholds(row) = &
      control%st_parameter * minval( matrix%val(row_start:row_end) )
    if (thresholds(row) .ge. 0) then  ! this row has no negative connections
     if (control%c_fail == 1) then ! coarsening has failed
      info%flag = -12
      return
     else
      ! subtract from n_str any strong connections
      ! already counted for this point
      n_str = n_str - point_type(row)
      ! make this an unconnected point
      point_type(row) = -2
      thresholds(row) = 0
      cycle row_loop       ! leave this row
     end if
    end if

    ! this row has off-diagonals and strong connections
    do i = row_start+1, row_end  ! run through off-diagonal columns
     ! test if connection is strong
     if ( matrix%val(i) .le. thresholds(row) ) then  ! a strong connection
      ! store number of strong transpose connections in point_type if
      ! the connection isn't to an unconnected point
      if ( point_type(matrix%col(i)) .ne. -2) then
       point_type(matrix%col(i)) = point_type(matrix%col(i)) + 1
       ! count the number of strong connections
       n_str = n_str + 1
      end if
    end if
   end do
   end if
  end do row_loop

  ! test for the case when all rows have only positive values
  if (n_str == 0) then
   info%flag = -12
   return
  end if

  if (keep%st_method == 1) then
   call MC65_matrix_is_symmetric(matrix, symmetric, MC65_status, pattern, &
        stat=info%stat)
   if (MC65_status .ne. 0) then
    if (MC65_status == MC65_ERR_MEMORY_ALLOC) then
     info%flag = -10
    else if (MC65_status == MC65_ERR_MEMORY_DEALLOC) then
     info%flag = -11
    end if
    return
   end if

   if (.not. symmetric) then
    keep%st_method=2
    info%flag = 1
   end if
  end if

 end subroutine test_connections

 ! --------------------
 ! -- find_strong_t --
 ! --------------------

 ! A subroutine to find the strong transpose connections
 ! and store these in a zd11 object

 subroutine find_strong_t(matrix, info, thresholds, point_type, n_str, str_t)

  type(zd11_type), intent(in) :: matrix
  type(mi20_info), intent(inout) :: info
  integer, intent(in) :: n_str
  integer, intent(in), dimension(:) :: point_type
  real(kind=myreal), intent(in), dimension(:) :: thresholds
  type(zd11_type), intent(out) :: str_t

  ! use to point at next position to be filled in str_t for each row
  integer, allocatable, dimension(:) :: ptr

  integer :: row_start, row_end  ! position of start and end of row in matrix
  integer :: row                 ! current row
  integer :: i                   ! general purpose integer
  integer :: strong
  integer :: next_position
  integer :: alloc_stat      ! used to find de(allocation) status

  ! allocate storage for ptr in str_t
  str_t%m = matrix%m
  allocate(str_t%ptr(1+matrix%m), stat=alloc_stat)
  if (alloc_stat .ne. 0) go to 20
  str_t%ptr(:) = 0

  ! first loop finds the number of str_t connections per row and stores the
  ! value in str_t%ptr, so that str_t%ptr(i+1) stores the value of the
  ! number of str_t connections in row i
  do row = 1, matrix%m
   if (point_type(row) == -2) cycle    ! bypass unconnected rows

   ! get number of strong transpose connections from point_types
   str_t%ptr(row+1) = point_type(row)
  end do

  ! this loop calculate the actual entries for str_t%ptr
  str_t%ptr(1) = 1
  do row = 2, matrix%m
   str_t%ptr(row) = str_t%ptr(row) + str_t%ptr(row-1)
  end do
  str_t%ptr(1+matrix%m) = n_str+1

  ! allocate str_t%col
  allocate(str_t%col(n_str), stat=alloc_stat)
  if (alloc_stat .ne. 0) go to 20

  ! allocate ptr (keeps track of next position to be filled in each row)
  allocate(ptr(matrix%m), stat = alloc_stat)
  if (alloc_stat .ne. 0) go to 20
  ptr = 0

  ! loop to put values in str_t%col
  do row = 1, matrix%m
   if (point_type(row) .ne. -2) then ! bypass unconnected rows
    ! find position of columns for this row
    row_start = matrix%ptr(row)
    row_end = matrix%ptr(row+1)-1
    ! run through columns
    do i = row_start+1, row_end  ! bypass diagonal
     ! bypass connections to unconnected points
     if ( point_type(matrix%col(i)) .ne. -2) then
      if ( matrix%val(i) .le. thresholds(row)) then  ! strong connection
       ! store 'row' in str_t for point 'strong'
       strong = matrix%col(i)
       next_position = str_t%ptr(strong) + ptr(strong)
       str_t%col(next_position) = row
       ptr(strong) = ptr(strong) + 1
      end if
     end if
    end do
   end if
  end do

  ! deallocation
  deallocate(ptr, stat=alloc_stat)
  if (alloc_stat .ne. 0) then
   info%stat = alloc_stat; info%flag = -11
   return
  end if

 20 continue
  ! allocation error
  if (alloc_stat .ne. 0) then
   info%stat = alloc_stat; info%flag = -10
  end if

 end subroutine find_strong_t


 ! --------------------
 ! -- generate_lists --
 ! --------------------

 ! A subroutine to generate linked lists of points (i.e. rows)
 ! each list contains points with the same weight

 subroutine generate_lists(point_type, info, lists, lists_ends, &
                            max_weight)

  ! array containing point type data
  integer, dimension(:), intent(inout) :: point_type

  type(mi20_info), intent(inout) :: info

  ! lists is an array of list_nodes, with one list_node object for
  ! each row of the matrix (i.e. for each point). The aim is to
  ! use this array to set up a series of linked lists, where each
  ! list containts points with the same weight
  type(list_node), intent(inout), dimension(:), target :: lists

  ! The list_node at the w th position in the lists_ends array points
  ! to the first and last list_node in the array lists with weight w
  type(list_node), intent(inout), allocatable, dimension(:), &
     target :: lists_ends

  integer, intent(out) :: max_weight ! returns the maximum value of weight

  integer :: n_points ! the number of points
  integer :: point    ! the current point

  ! set starting value of weight and max_weight
  n_points =size(point_type)
  lists(1:n_points)%weight = 0
  max_weight = 0

  ! set up the lists...

  do point = 1,n_points

   ! bypass points with no connections
   if (point_type(point) .ne. -2) then

    ! set the starting weight (=number of strong transpose connections)
    lists(point)%weight = point_type(point)

    ! set point_type to zero
    point_type(point) = 0
    lists(point)%point = point

   end if

   ! bypass points with no strong transpose connections
   if (lists(point)%weight>0) then

    ! add point to list
    call list_add_point(point, lists, lists_ends, info)

    ! stop if an error occurs
    if (info%flag < 0) return

    ! update max_weight
    if (lists(point)%weight > max_weight) then
     max_weight = lists(point)%weight
    end if

   end if

  end do

 end subroutine generate_lists


 ! --------------------
 ! -- list_add_point --
 ! --------------------

 ! A subroutine to add a point to its associated list

 subroutine list_add_point(point, lists, lists_ends, info)

  integer, intent(in) :: point

  type(list_node), intent(inout), dimension(:), target :: lists
  type(list_node), intent(inout), allocatable, dimension(:), &
     target :: lists_ends
  type(mi20_info), intent(inout) :: info

  integer :: weight            ! weight of point to be added to a list
  integer :: size_lists_ends   ! size of lists_ends
  integer :: alloc_stat        ! used to capture (de)allocation status

  ! temporary array used when extending lists_ends
  type(list_node), allocatable, dimension(:) :: lists_temp

  size_lists_ends = size(lists_ends)
  weight = lists(point)%weight

  ! increase size of list_ends if necessary
  if (weight .gt. size_lists_ends) then

   ! allocate lists_temps
   allocate(lists_temp(size_lists_ends), stat=alloc_stat)
   if (alloc_stat .ne. 0) then
    info%stat = alloc_stat; info%flag = -10
    return
   end if

   ! copy lists_ends to lists_temp
   lists_temp(1:size_lists_ends) = lists_ends(1:size_lists_ends)

   ! deallocate lists_ends
   deallocate(lists_ends, stat=alloc_stat)
   if (alloc_stat .ne. 0) then
    info%stat = alloc_stat; info%flag = -11
    deallocate(lists_temp)
    return
   end if

   ! allocate new lists_ends
   allocate(lists_ends(2*weight), stat=alloc_stat)
   if (alloc_stat .ne. 0) then
    info%stat = alloc_stat;  info%flag = -10
    deallocate(lists_temp)
    return
   end if

   ! copy from lists_temp to lists_ends
   lists_ends(1:size_lists_ends) = lists_temp(1:size_lists_ends)

   ! deallocate lists_temp
   deallocate(lists_temp, stat=alloc_stat)
   if (alloc_stat .ne. 0) then
    info%stat = alloc_stat;  info%flag = -11
    return
   end if
   size_lists_ends = weight*2
  end if

  if ( associated(lists_ends(weight)%next) ) then
   ! i.e. a list for this weight already exists
   ! point 'point' to end of list
   lists(point)%last => lists_ends(weight)%last
   ! nullify next pointer from 'point' since it is the last point
   lists(point)%next => null()
   ! point last point in the list to 'point'
   lists_ends(weight)%last%next => lists(point)
   ! point list_end to 'point'
   lists_ends(weight)%last => lists(point)
  else
   ! a list for this weight has not been started
   lists_ends(weight)%next => lists(point)
   lists_ends(weight)%last => lists(point)
   lists(point)%last => null()
   lists(point)%next => null()
  end if

 end subroutine list_add_point


 ! ------------------------
 ! -- list_promote_point --
 ! ------------------------

 ! A subroutine to increase the weight of point by the amount weight_inc
 ! and then promot point to the correct list

 subroutine list_promote_point(point, weight_inc, lists, lists_ends, info)

  integer, intent(in) :: point
  integer, intent(in) :: weight_inc
  type(list_node), intent(inout), dimension(:), target :: lists
  type(list_node), intent(inout), allocatable, dimension(:), target &
                                 :: lists_ends
  type(mi20_info), intent(inout) :: info

  integer :: weight

  weight = lists(point)%weight

  ! remove point from list
  call list_remove_point(point, lists, lists_ends)

  lists(point)%weight = weight + weight_inc

  ! add point to next list
  call list_add_point(point, lists, lists_ends, info)

 end subroutine list_promote_point


 ! -----------------------
 ! -- list_remove_point --
 ! -----------------------

 ! A subroutine to remove a point from a list

 subroutine list_remove_point(point, lists, lists_ends)

  integer, intent(in) :: point
  type(list_node), intent(inout), dimension(:), target :: lists
  type(list_node), intent(inout), dimension(:), target :: lists_ends

  integer :: weight

  weight = lists(point)%weight

  ! if weights == 0 point is not in a list so do nothing
  if (weight .ne. 0) then

   ! case with only one point in the list
   if ( (.not. associated(lists(point)%last)) .and. &
         (.not. associated(lists(point)%next)) ) then
     lists_ends(weight)%next => null()
     lists_ends(weight)%last => null()

    ! case where it is the first point in the list
    else if (.not. associated(lists(point)%last)) then
     ! point list start at next 'point' in list
     lists_ends(weight)%next => lists(point)%next
     ! make last of 'first point' in list point to null()
     lists(point)%next%last => null()

    ! case where it is the last point in list
    else if (.not. associated( lists(point)%next)) then
     ! point list end at previous 'point' in list
     lists_ends(weight)%last => lists(point)%last
     ! make last of 'first point' in list point to null()
     lists(point)%last%next => null()
    else
     lists(point)%next%last => lists(point)%last
     lists(point)%last%next => lists(point)%next
    end if
  end if

 end subroutine list_remove_point


!! ------------------
! ! -- output_lists --
! ! ------------------

! ! A subroutine to output list to screen for testing

!! subroutine output_lists(lists, lists_ends, control)

!!  type(list_node), intent(in), dimension(:), target :: lists
!!  type(list_node), intent(in), dimension(:), target :: lists_ends
!!  type(mi20_control), intent(in) :: control

!!  integer :: int1, int2, size_lists_ends

!!  size_lists_ends = size(lists_ends)
!!  write(control%print,*)
!!  write(control%print,*) "size of lists_ends is", size_lists_ends

!!  ! output weights
!!   do int1 = 1,size_lists_ends
!!    if ( associated(lists_ends(int1)%next) ) then
!!     int2 = lists_ends(int1)%next%point
!!     write(control%print,*)
!!     write(control%print,*) "weight:", int1
!!     write(control%print,*) "-----------"
!!     write(control%print,*) "point:", int2
!!     do
!!      if ( .not. associated(lists(int2)%next) ) exit
!!      int2 =  lists(int2)%next%point
!!      write(control%print,*) "point:", int2
!!     end do
!!    end if
!!   end do

!! end subroutine output_lists


 ! ------------------------
 ! -- coarsen_first_pass --
 ! ------------------------

 ! A subroutine to perform the first coarsening pass:
 ! An undecided point with maximum weight becomes a new C point
 ! All undecided points with strong transpose connection to the new C become
 ! new F points
 ! The weights of the undecided points are increased for each strong
 ! connection to a new F point.

 subroutine coarsen_first_pass(matrix, str_t, thresholds, weight_inc, &
   lists, lists_ends, point_type, max_weight, size_cmatrix, info, keep)

  type(zd11_type), intent(in) :: matrix
  type(zd11_type), intent(in) :: str_t
  real(kind=myreal), intent(in), dimension(:) :: thresholds
  integer, intent(in) :: weight_inc
  type(list_node), intent(inout), dimension(:), target :: lists
  type(list_node), intent(inout), allocatable, dimension(:), &
      target :: lists_ends
  integer, intent(inout), dimension(:) :: point_type
  integer, intent(inout) :: max_weight    ! current maximum weight
  integer, intent(inout) :: size_cmatrix  ! size of coarse level matrix
  type(mi20_info), intent(inout) :: info
  type(mi20_keep), intent(in) :: keep

  integer :: row_i_start       ! start of row associated with point i
  integer :: row_i_end         ! end of row associated with point i
  integer :: new_C             ! point number of new C
  integer :: index_i, index_k  ! loop indicies
  integer :: point_i, point_k  ! point numbers
  integer :: pos               ! current position in the row

  size_cmatrix = 0

  do
   ! check there are points with the current max_weight, if not reduce weight
   ! and test again
   do
    if ( associated( lists_ends(max_weight)%next ) ) exit
     max_weight = max_weight - 1
    if (max_weight == 0) exit
   end do

   ! first pass ends when the maximum weight is zero
   if (max_weight == 0) exit

   ! find the new C point
   new_C = lists_ends(max_weight)%next%point
   size_cmatrix = size_cmatrix + 1
   point_type(new_C) = 1

   ! remove new_C from lists
   call list_remove_point(new_C, lists, lists_ends )
   lists(new_C)%next => null()
   lists(new_C)%last => null()
   lists(new_C)%point = -lists(new_C)%point
   lists(new_C)%weight = -1

   ! Find points with a strong connection to new_C, i.e. points which new_C
   ! has a strong transpose connection to. Make these point F. Then for each
   ! strong connection from the new F points to an undecided point, increment
   ! the weight of the undecided point.

   ! for keep%st_method=2, the strong transpose connections have already been
   ! found and are stored in str_t
   if (keep%st_method == 2) then
    ! Run through points 'point_i' with a strong connection to new_C,
    ! ie str_t connections for new_C
    do index_i = str_t%ptr(new_C), str_t%ptr(new_C+1)-1
     point_i = str_t%col(index_i)
     ! continue only for undecided points
     if (point_type(point_i) .ne. 0) cycle
     ! set point type for this new F point
     point_type(point_i) = -1

     ! remove point_i from lists
     call list_remove_point(point_i, lists, lists_ends)
     lists(point_i)%next => null()
     lists(point_i)%last => null()
     lists(point_i)%point = -lists(point_i)%point
     lists(point_i)%weight = -2

     ! increment weights of all undecided points 'point_k' which point_i
     ! has a strong connection to
     row_i_start = matrix%ptr(point_i)
      row_i_end = matrix%ptr(point_i+1)-1
     do index_k = row_i_start+1, row_i_end ! bypass diagonal
      point_k = matrix%col(index_k)
      ! test that the point is undecided and the connection is striong
      if ( (point_type(point_k) == 0) .and. &
           ( matrix%val(index_k) .le. thresholds(point_i)) ) then
       ! increment the weight
       call list_promote_point(point_k, weight_inc, &
                          lists, lists_ends, info)
       ! stop if an error occurs
       if (info%flag<0) return

       if ( lists(point_k)%weight > max_weight ) then
        max_weight = lists(point_k)%weight
       end if
      end if
     end do
    end do

   ! this method finds the strong transpose connections on the fly, and
   ! requires the matrix storage to be symmetric
   else if (keep%st_method == 1) then
    ! for each column i corresponding to an entry in row new_C,
    ! run through the corresponding row i and find column new_C,
    ! then test the value at this entry to see if the connection from
    ! row i to row new_C is strong, if it is then we have found a
    ! strong transpose connection to new_C
    do index_i = matrix%ptr(new_C)+1, matrix%ptr(new_C+1)-1
     point_i = matrix%col(index_i)
     ! continue only for undecided points
     if (point_type(point_i) .ne. 0) cycle

     ! run through row_i and find column new_C
     row_i_start = matrix%ptr(point_i)
     row_i_end = matrix%ptr(point_i+1)-1
     do pos = row_i_start+1, row_i_end
      if ( matrix%col(pos) == new_C ) then  ! we have found new_C
       ! test if connection is strong
       if ( matrix%val(pos) .le. thresholds(point_i) ) then
        ! set point_i to be F point
        point_type(point_i) = -1

         ! remove point_i from lists
        call list_remove_point(point_i, lists, lists_ends)
        lists(point_i)%next => null()
        lists(point_i)%last => null()
        lists(point_i)%point = -lists(point_i)%point
        lists(point_i)%weight = -2

        ! next increment weights of all undecided points 'point_k'
        ! which point_i has a strong connection to
        do index_k = row_i_start+1, row_i_end
         point_k = matrix%col(index_k)
         if ( (point_type(point_k) == 0) .and. &
            ( matrix%val(index_k) .le. thresholds(point_i)) ) then
          call list_promote_point(point_k, weight_inc, &
                          lists, lists_ends, info)
          ! stop if an error occurs
          if (info%flag<0) return

          if ( lists(point_k)%weight > max_weight ) then
            max_weight = lists(point_k)%weight
           end if
          end if
        end do
       end if
       exit
      end if
     end do
    end do
   end if
  end do

 end subroutine coarsen_first_pass

 ! -------------------------
 ! -- coarsen_second_pass --
 ! -------------------------

 ! A subroutine to ensure that all pairs of strongly connected F points
 ! share a common C point

 subroutine coarsen_second_pass(matrix, thresholds, &
                str_t, keep, size_cmatrix, point_type)

  type(zd11_type), intent(in) :: matrix
  real(kind=myreal), intent(in), dimension(:) :: thresholds
  type(zd11_type), intent(in) :: str_t
  type(mi20_keep), intent(in) :: keep
  integer, intent(inout) :: size_cmatrix  ! size of coarse level matrix
  integer, intent(inout), dimension(:) :: point_type

  integer :: point_i         ! 1st F point in pair being tested
  integer :: point_j         ! 2nd F point in pair being tested
  ! i.e. we perform the test when point_j has a strong connection to point_i

  integer :: row_i_start, row_i_end  ! start and end of row i

  integer :: row_j_start, row_j_end  ! start and end of row j
  integer :: index_j         ! index associated with point_j


  integer :: point_ii     ! point connected to point_i when testing
                          ! for common C point
  integer :: point_jj     ! point connected to point_j when testing
                          ! for common C point
  ! i.e. if point_ii and point_jj are C points, with strong connections
  ! to point_i and point_j respectively, then the test is successful

  integer :: index_ii     ! index associated with point_ii
  integer :: index_jj     ! index associated with point_ii
  integer :: size_matrix  ! size of current matrix
  integer :: new_C        ! value of new C point

  logical :: test         ! .true. when test is passed


  ! number of times the current pair of F points has failed
  integer :: no_fails

  size_matrix = matrix%m

  ! Run through all points 'point_i'
  do point_i = 1, size_matrix
   ! bypass if point_i isn't an F point
   if ( point_type(point_i) == -1 ) then
    ! set number of fails for point_i to be 0
    no_fails = 0

    ! Run through each F point 'point_j' with a strong connection to point_i
    ! and test for common C
    row_i_start = matrix%ptr(point_i)
    row_i_end = matrix%ptr(point_i+1)-1
     do index_j = row_i_start+1, row_i_end ! bypass diagonal
     point_j = matrix%col(index_j)
     ! continue for only F points and strong connections to point_i
      if ( (point_type(point_j) == -1) .and. &
            (matrix%val(index_j) .le. thresholds(point_i)) ) then

      ! perform test for this pair of F points
      test = .false.
      ! run through all points 'point_ii' connected to point_i
      test_loop: do index_ii = row_i_start+1, row_i_end
       point_ii = matrix%col(index_ii)
       ! continue if point_ii is a C point
       ! with strong connection to point_i
       if ( (point_type(point_ii) .gt. 0) .and. &
               (matrix%val(index_ii) .le. thresholds(point_i)) ) then
        ! run through all points 'point_jj' connected to point_j
        row_j_start = matrix%ptr(point_j)
        row_j_end = matrix%ptr(point_j+1)-1
        do index_jj = row_j_start+1, row_j_end
         point_jj = matrix%col(index_jj)
         ! continue if point_jj is a C point with a strong connection
         ! to point_j and look for a match to point_ii
         if ( (point_type(point_jj) .gt. 0) .and. &
               (matrix%val(index_jj) .le. thresholds(point_j)) ) then
           ! test for common C point
          test = (point_ii == point_jj)
          if ( test ) exit test_loop
         end if
        end do
       end if
      end do test_loop

      ! if test fails we need a new C point
      if (.not. test) then
       no_fails = no_fails + 1
       ! if this is the first failure for point_i then
       ! tentatively set point_j to be the new_C point
       if (no_fails == 1) then
        new_C = point_j
       end if
      end if
     end if
     ! No further testing is needed if 2 fails have already occured
     ! for point_i, otherwise continue testing point_i against any
     ! remaining strong connections to F points
     if (no_fails == 2) exit
    end do
    ! end of testing for point_i

    if (no_fails > 0) then  ! add a new C point
     ! if 2 fails have occurred, change the new_C point to be point_i
      if (no_fails == 2) then
       new_C = point_i
      end if

     call add_new_C(new_C, matrix, str_t, keep, &
                thresholds, point_type, size_cmatrix)

    end if
   end if
  end do

 end subroutine coarsen_second_pass


 ! ------------------
 ! -- add_new_C --
 ! ------------------

 ! A subroutine which allows subroutine coarsen_second_pass to convert a
 ! point to a new C point, and also and convert any undecided points to
 ! be F points if they have a strong connection to the new C point, whilst
 ! also ensuring that strongly connected F-F pairs without strong
 ! connections to a common C point are not generated.
 ! By the second pass an undecided point has no strong transpose
 ! connections (else it would have become a C point by now) and has no strong
 ! connections to a C point (else it would have become an F point earlier).
 ! Therefore if an undecided point is converted to an F point, and if it has
 ! a strong connection to another F point, it is not possible this pair of F
 ! points share a common C point.
 ! Also any new F point can't become a C point since it has no strong
 ! transpose connections.
 ! Therefore if there is a connection from a new F point to an existing
 ! F point, we must convert the existing F point to a C point. After this
 ! it is necessary to convert any undecided point with a strong transpose
 ! connections to the new C point to an F points, and so on, hence the
 ! need for a recursive subroutine.

 recursive subroutine add_new_C(new_C, matrix, str_t, keep, &
                thresholds, point_type, size_cmatrix)

  integer, intent(in) :: new_C          ! new C point to be added
  type(zd11_type), intent(in) :: matrix ! matrix being coarsened
  type(zd11_type), intent(in) :: str_t  ! strong transpose connections
                                        ! if st_method 2 is used
  type(mi20_keep), intent(in) :: keep   ! keep
  real(kind=myreal), intent(in), dimension(:) :: thresholds
  integer, intent(inout), dimension(:) :: point_type
  integer, intent(inout) :: size_cmatrix   ! number of C points

  integer :: next_new_C    ! next new C point
  integer :: new_F         ! new F point

  ! generate purpose integers used when looping through rows
  integer :: index_k, index_l
  integer :: point_k, point_l
  integer :: row_k_start, row_k_end, pos

  point_type(new_C) = 1
  size_cmatrix = size_cmatrix+1

  ! for st_method 2 we have strong transpose connections stored in str_t
  if (keep%st_method == 2) then
   do index_k = str_t%ptr(new_C), str_t%ptr(new_C+1)-1
    point_k = str_t%col(index_k)

    if (point_type(point_k) .ne. 0) cycle
     ! set undecided points to be F
     new_F = point_k
     point_type(new_F) = -1
     ! now find any F points with a strong connection to new_F
     ! and set them to be C points
     do index_l = matrix%ptr(new_F)+1, matrix%ptr(new_F+1)-1
      point_l = matrix%col(index_l)
      ! continue for only F points and strong connections ko new_F
       if ( (point_type(point_l) == -1) .and. &
          (matrix%val(index_l) .le. thresholds(new_F)) ) then
       ! set point_l to be new C
       next_new_C = point_l
       call add_new_C(next_new_C, matrix, str_t, keep, thresholds, &
           point_type, size_cmatrix)
      end if
     end do

   end do
  ! for st_method 1 we find strong transpose connections on the fly
  else if (keep%st_method == 1) then
   ! for each column k corresponding to an entry in row new_C,
   ! run through the corresponding row k and find column new_C,
   ! then test the value at this entry to see if the connection from
   ! row i to row new_C is strong, if it is then we have found a
   ! strong transpose connection to new_C

   do index_k = matrix%ptr(new_C)+1, matrix%ptr(new_C+1)-1
    point_k = matrix%col(index_k)

    ! continue only for undecided points
    if (point_type(point_k) .ne. 0) cycle

    ! run through row_k and find column new_C
    row_k_start = matrix%ptr(point_k)
    row_k_end = matrix%ptr(point_k+1)-1
    do pos = row_k_start+1, row_k_end
     if ( matrix%col(pos) == new_C ) then  ! we have found new_C
      ! test if connection is strong
      if ( matrix%val(pos) .le. thresholds(point_k) ) then
       ! set undecided points to be F
       new_F = point_k
       point_type(new_F) = -1
       ! now find any F points with a strong connection to new_F
       ! and set them to be C points
       do index_l = matrix%ptr(new_F)+1, matrix%ptr(new_F+1)-1
        point_l = matrix%col(index_l)
        ! continue for only F points and strong connections to new_F
         if ( (point_type(point_l) == -1) .and. &
              (matrix%val(index_l) .le. thresholds(new_F)) ) then
         ! set point_l to be C
         next_new_C = point_l
         call add_new_C(next_new_C, matrix, str_t, keep, &
           thresholds, point_type, size_cmatrix)
        end if
       end do
      end if
     end if
    end do
   end do
  end if

 end subroutine add_new_C


 ! ------------------------
 ! -- coarsen_third_pass --
 ! ------------------------

! A subroutine to deal with any remaining undecided points.
! Any undecided point must have a strong connection, otherwise it would
! have been flagged as unconnected in subroutine test_connections.
! However any undecided point has no strong transpose connections,
! else it would have become a C point.
! And any undecided point has no strong connection to a C point,
! else the undecided point would have become an F point by now.
! So either:
!   the undecided point must have a strong connection to an F point, in
!   which case set the F point to be a new C point and the undecided
!   point to be F. Note: Stuben recommends indirect interpolation
!   for such points
! or:
!   the undecided point has strong connections only to points which
!   have been set as unconnected earlier, in which case set the
!   undecided point also to be unconnected

 subroutine coarsen_third_pass(matrix, str_t, thresholds, keep, &
                           size_cmatrix, point_type)

  type(zd11_type), intent(in) :: matrix
  type(zd11_type), intent(in) :: str_t
  real(kind=myreal), intent(inout), dimension(:) :: thresholds
  type(mi20_keep), intent(in) :: keep
  integer, intent(inout) :: size_cmatrix
  integer, intent(inout), dimension(:) :: point_type

  integer :: size_matrix  ! size of matrix being coarsened
  integer :: weight       ! weight used to decide which F point becomes C
  integer :: max_weight   ! current maximum weight
  integer :: new_C        ! value of new_C point
  integer :: point_i      ! undecided point currently under incestigation
  integer :: row_i_start  ! starting position of row i
  integer :: row_i_end    ! end if row i
  integer :: point_j      ! point with a connection to point_i
  integer :: index_j      ! index associated with point_j
  integer :: row_j_start  ! start of row associated with point_j
  integer :: row_j_end    ! end of row associated with point_j
  integer :: point_k      ! point with a transpose connection to point_j
  integer :: index_k      ! index associated with point_k
  integer :: pos          ! current position in row k
  integer :: row_k_start  ! start of row k
  integer :: row_k_end    ! end of row k

  size_matrix = matrix%m
  do point_i = 1,size_matrix
   ! continue only for undecided points
   if ( point_type(point_i) == 0 ) then
    ! run through all points 'point_j' with a connection to point_i
    ! for strong F connections find the F point which influences the
    ! most undecided points
    max_weight = 0
    row_i_start = matrix%ptr(point_i)
    row_i_end = matrix%ptr(point_i+1)-1

    do index_j = row_i_start+1, row_i_end
     point_j = matrix%col(index_j)
     ! continue only for strong F connections to point_i
     if ( (point_type(point_j) == -1) .and. &
            ( matrix%val(index_j) .le. thresholds(point_i)) ) then
      ! for each point_j, find the number of undecided points which have
      ! a strong connection to point_j, set the point_j with the most
      ! to be a C point
      weight = 0

      if (keep%st_method == 1) then
       ! for each undecided 'point_k' with a strong transpose
       ! connection to 'point_j' increase the weight
       row_j_start = matrix%ptr(point_j)
       row_j_end = matrix%ptr(point_j+1)-1
       do index_k = row_j_start+1, row_j_end   ! run through this row
        ! find associated point number
        point_k = matrix%col(index_k)
        ! continue only for undecided points
        if (point_type(point_k) .ne. 0) cycle
        ! run through row associated with point_k
        row_k_start = matrix%ptr(point_k)
        row_k_end = matrix%ptr(point_k+1)-1
        do pos = row_k_start+1, row_k_end
         ! find connection to point_j
         if ( matrix%col(pos) == point_j ) then
           ! test whether point_k has a strong connection to point j
           if ( matrix%val(pos) .le. thresholds(point_k) ) then
           ! if it does increment the weight for point_j
           weight = weight + 1
          end if
          exit
         end if
        end do
       end do

      else if (keep%st_method == 2) then
       ! for each undecided 'point_k' with a strong transpose
       ! connection to 'point_j' increase the weight
       row_j_start = str_t%ptr(point_j)
       row_j_end = str_t%ptr(point_j+1)-1
       do index_k = row_j_start, row_j_end   ! run through this row
        ! find point number
        point_k = str_t%col(index_k)
        ! if this point is undecided increment the weight
        if (point_type(point_k) == 0) then
         weight = weight + 1
        end if
       end do
      end if

      ! if this point_j has the most strong transpose connections
      ! to undecided points, make it the new C point
      if (weight > max_weight) then
       max_weight = weight
       new_C = point_j
      end if
     end if
    end do

    ! if no F points were found then then set point_i as unconnected
    if (max_weight .eq. 0) then
     point_type(point_i) = -2
     thresholds(point_i) = 0
    else
     ! at this point we have a new C point for point_i
     point_type(new_C) = 1
     size_cmatrix = size_cmatrix + 1

     ! change any undecided points with a strong connection to new_C to be
     ! F points
     if (keep%st_method == 2) then
      do index_k = str_t%ptr(new_C), str_t%ptr(new_C+1)-1
       point_k = str_t%col(index_k)
       ! continue only for undecided points
       if (point_type(point_k) .eq. 0) then
        point_type(point_k) = -1
       end if
      end do

     else if (keep%st_method == 1) then
      ! find rows which contain a connection to new_C
      row_j_start = matrix%ptr(new_C)
      row_j_end = matrix%ptr(new_C+1)-1
      ! run through these rows
      do index_k = row_j_start+1, row_j_end
       ! find associated point number
       point_k = matrix%col(index_k)
       ! continue only for undecided points
       if (point_type(point_k) .ne. 0) cycle
       row_k_start = matrix%ptr(point_k)
       row_k_end = matrix%ptr(point_k+1)-1
       ! run through points and find the new_C point
       do pos = row_k_start+1, row_k_end
        if ( matrix%col(pos) == new_C ) then
         ! test whether point_k has a strong connection to new_C
         if ( matrix%val(pos) .le. thresholds(point_k) ) then
          ! set undecided points to be F
          point_type(point_k) = -1
         end if
         exit
        end if
       end do
      end do
     end if
    end if
   end if
  end do

 end subroutine coarsen_third_pass


 ! -------------------------
 ! -- error_message_setup --
 ! -------------------------

 ! A subroutine to output warning and error messages relating to mi20_setup,
 ! excluding messages relating to mi20_control and the input zd11
 ! matrix, since subroutines test_control_setup and test_matrix are better
 ! suited for this.

subroutine error_message_setup(control, info, level, step)

  type(mi20_info), intent(in) :: info
  type(mi20_control), intent(in) :: control
  integer :: level             ! current coarsening level
  integer, optional :: step    ! current coarsening step

  if (info%flag .ne. 0) then

   if(control%aggressive > 1) then
    if (info%flag < 0) then
     write(control%error,'(a/3(a,i4))') " Error return from mi20_setup", &
     " level ", level,", step ", step, ". info%flag = " ,info%flag
    else
     write(control%error,'(a/3(a,i4))') &
         " Warning. Coarsening failure in mi20_setup", &
         " level ", level, ", step ", step, ". info%flag = " ,info%flag
    end if
   else
    if (info%flag < 0) then
     write(control%error,'(a/3(a,i4))')" Error return from mi20_setup", &
     " level ", level,", info%flag = " ,info%flag
    else
     write(control%error,'(a/3(a,i4))') &
         " Warning. Coarsening failure in mi20_setup", &
         " level ", level, ", info%flag = " ,info%flag
    end if
   end if

   ! error messages
   if (abs(info%flag) == 10) then                 ! allocation error
     write(control%error,'(a)') " allocation failure"
   else if (abs(info%flag) == 11) then                 ! deallocation error
     write(control%error,'(a)') " deallocation failure"
   else if (abs(info%flag) == 2) then                  ! missing diagonal
     write(control%error,'(a)') " zd11_type input argument has &
            &missing diagonal entries"
   ! coarsening errors from positive rows
   else if (abs(info%flag) == 12) then
     if (control%c_fail == 1) then
     write(control%error,'(a)') &
       " row found with either no negative off-diagonals or ", &
       " all connections are to rows with no negative off-diagonals"
    else
     write(control%error,'(a)') &
       " all rows have either no negative off-diagonals or ", &
       " all connections are to rows with no negative off-diagonals"
    end if
   else if (info%flag == 13) then                   ! coarsening stagnation
     write(control%error,'(a)') " coarsening stagnation"
   end if

  end if

 end subroutine error_message_setup


 ! -------------------------------
 ! -- error_message_precondition--
 ! -------------------------------

 ! A subroutine to output warning and error messages relating to
 ! mi20_precondition, exclusing messages relating to mi20_control as
 ! subroutines test_control_precondition is better suited for this.

subroutine error_message_precondition(control, info)

  type(mi20_info), intent(in) :: info
  type(mi20_control), intent(in) :: control

  if (info%flag .ne. 0) then

   if (info%flag < 0) then
    write(control%error,'(a,i5)') &
        " Error return from mi20_precondition, info%flag = " ,info%flag
   else
    write(control%error,'(a,i5)') &
        " Warning from mi20_precondition, info%flag = " ,info%flag
   end if

   ! error messages
   if (info%flag .eq. -10) then                   ! allocation error
     write(control%error,'(a)') " allocation failure"
   else if (info%flag .eq. -11) then              ! deallocation error
     write(control%error,'(a)') " deallocation failure"
    else if (info%flag .eq. -14) then             ! vector norm too large
     write(control%error,'(a)') " vector norm exceeded error criteria"
    else if (info%flag .eq. -15) then             ! setup required
     write(control%error,'(a)') &
    " mi20_setup must be called successfully before mi20_precondition"
    else if (info%flag .eq. -16) then             ! vectors too small
     write(control%error,'(a)') " size of input vector is too small"
    else if (info%flag .eq. -17) then             ! lapack error
     write(control%error,'(a)') &
      " error from lapack dense solver _getrf on coarsest level"
    else if (info%flag .eq. -18) then             ! hsl_ma48 error
     write(control%error,'(a)') &
      " error from dense solver hsl_ma48 on coarsest level"
    else if (info%flag .eq. -19) then             ! ma48_control required
     write(control%error,'(a)') " user must supply ma48_control type"
    else if (info%flag .eq. 20) then   !  user has requested too many levels
     write(control%error,'(a)') " number of levels requested exceeds &
                           &number of levels available"
   end if

  end if

 end subroutine error_message_precondition


 ! ------------------------
 ! -- setup_deallocation --
 ! ------------------------

 ! A subroutine to deallocate any memory allocated in mi20_setup

 subroutine setup_deallocation(info, point_type, thresholds, &
                          lists, lists_ends, str_t)

  type(mi20_info), intent(inout) :: info
  integer, allocatable, dimension(:), intent(inout) :: point_type
  real(kind=myreal), allocatable, dimension(:), intent(inout) :: thresholds
  type(list_node), allocatable, dimension(:), intent(inout), target :: lists
  type(list_node), allocatable, dimension(:), intent(inout), &
     target :: lists_ends
  type(zd11_type), intent(inout) :: str_t

  integer :: alloc_stat ! used to find deallocation status

  if ( allocated(point_type) ) then
   deallocate(point_type, stat=alloc_stat)
   if (alloc_stat .ne. 0) go to 10
  end if
  if ( allocated(lists) ) then
   deallocate(lists, stat=alloc_stat)
   if (alloc_stat .ne. 0) go to 10
  end if
  if ( allocated(lists_ends) ) then
   deallocate(lists_ends, stat=alloc_stat)
   if (alloc_stat .ne. 0) go to 10
  end if
  if ( allocated(thresholds) ) then
   deallocate(thresholds, stat=alloc_stat)
   if (alloc_stat .ne. 0) go to 10
  end if
  if ( allocated(str_t%ptr) ) then
   deallocate(str_t%ptr, stat=alloc_stat)
   if (alloc_stat .ne. 0) go to 10
  end if
  if ( allocated(str_t%col) ) then
   deallocate(str_t%col, stat=alloc_stat)
   if (alloc_stat .ne. 0) go to 10
  end if
  return

 10 info%stat = alloc_stat
    info%flag=-11

 end subroutine setup_deallocation


 ! -------------------
 ! -- mi20_finalize --
 ! -------------------

 ! A subroutine to deallocate memory associated with objects of type
 ! mi20_data and mi20_keep

 subroutine mi20_finalize(coarse_data, keep, control, info, ma48_cntrl)

  type(mi20_data), allocatable, dimension(:), intent(inout) :: coarse_data
  type(mi20_keep), intent(inout) :: keep
  type(mi20_control), intent(in) ::control
  type(mi20_info), intent(inout) :: info
  type(ma48_control), optional, intent(inout) :: ma48_cntrl

  integer :: level          ! current level being deallocated
  integer :: n_levels       ! length of array
  integer :: alloc_stat     ! used to find deallocation status
  logical :: deallocate_ma48_matrix = .true.
  logical :: error_out      ! .true. if error messages are to be issued

  info%flag = 0
  ! switch on/off error messages
  ! ----------------------------
  error_out = .false.
  if ( (control%print_level .gt. 0) .and. (control%error .ge. 0) ) &
   error_out = .true.

  ! deallocate direct solver data in keep
  if ( present(ma48_cntrl) .and. (control%ma48 .gt. 0) ) then
   call dsolver_deallocation(keep, control, info, deallocate_ma48_matrix, &
                                  ma48_cntrl)
  else
   call dsolver_deallocation(keep, control, info, deallocate_ma48_matrix)
  end if

  ! set keep%ma48_data to 1 to ensure ma48_analyze is used next time the
  ! preconditioner is used
  keep%ma48_data = 1

  ! set the internal conversion routine flag to zero
  if (keep%zd11_internal_conversion) then
     if (allocated(keep%A%col)) deallocate(keep%A%col)
     if (allocated(keep%A%val)) deallocate(keep%A%val)
     if (allocated(keep%A%ptr)) deallocate(keep%A%ptr)
     keep%zd11_internal_conversion = .false.
  end if
  
  if (allocated(coarse_data)) then
   n_levels = size(coarse_data)
   info%flag = 0

   ! run through array
    do level = 1,n_levels
     ! deallocate A_mat
     if ( allocated(coarse_data(level)%A_mat%val) ) then
      deallocate(coarse_data(level)%A_mat%val, stat = alloc_stat)
      if (alloc_stat .ne. 0) go to 10
     end if
     if ( allocated(coarse_data(level)%A_mat%ptr) ) then
      deallocate(coarse_data(level)%A_mat%ptr, stat = alloc_stat)
      if (alloc_stat .ne. 0) go to 10
     end if
     if ( allocated(coarse_data(level)%A_mat%col) ) then
      deallocate(coarse_data(level)%A_mat%col, stat = alloc_stat)
      if (alloc_stat .ne. 0) go to 10
     end if
     if ( allocated(coarse_data(level)%A_mat%type) ) then
      deallocate(coarse_data(level)%A_mat%type, stat = alloc_stat)
      if (alloc_stat .ne. 0) go to 10
     end if
     ! deallocate I_mat
     if ( allocated(coarse_data(level)%I_mat%val) ) then
      deallocate(coarse_data(level)%I_mat%val, stat = alloc_stat)
      if (alloc_stat .ne. 0) go to 10
     end if
     if ( allocated(coarse_data(level)%I_mat%ptr) ) then
      deallocate(coarse_data(level)%I_mat%ptr, stat = alloc_stat)
      if (alloc_stat .ne. 0) go to 10
     end if
     if ( allocated(coarse_data(level)%I_mat%col) ) then
      deallocate(coarse_data(level)%I_mat%col, stat = alloc_stat)
      if (alloc_stat .ne. 0) go to 10
     end if
     if ( allocated(coarse_data(level)%I_mat%type) ) then
      deallocate(coarse_data(level)%I_mat%type, stat = alloc_stat)
      if (alloc_stat .ne. 0) go to 10
     end if
    end do

   ! deallocate mi20_data array
   deallocate(coarse_data, stat=alloc_stat)

   if (alloc_stat .ne. 0) go to 10
  end if

  return

 10 info%stat = alloc_stat
    info%flag = -11
    if ((error_out) ) then
       write(control%error,'(a,i5)') &
      " Error return from mi20_finalize, info%flag = " ,info%flag
     write(control%error,'(a)') " deallocation error"
    end if

 end subroutine mi20_finalize


 ! -----------------------
 ! -- gen_interpolation --
 ! -----------------------

 ! A subroutine to generate the next interpolation matrix
 subroutine gen_interpolation(matrix, point_type, thresholds, i_matrix, &
                               control, info)

  ! matrix being coarsened
  type(zd11_type), intent(in) :: matrix

  ! array containing the point type data produced during coarsening - for
  ! C points this stores the index on the next coarse level
  integer, intent(in), dimension(:) :: point_type

  ! array of thresholds defining strong connections
  real(kind=myreal), intent(in), dimension(:) :: thresholds

  ! interpolation matrix to be generated
  type(zd11_type), intent(out) :: i_matrix

  ! control parameters
  type(mi20_control), intent(in) :: control

  ! derived type containing info
  type(mi20_info), intent(inout) :: info

  ! arrays to temporarily hold col and val when using interpolation truncation
  integer, allocatable, dimension(:) :: col
  real(kind=myreal), allocatable, dimension(:) :: val

  integer :: size_matrix      ! size of matrix being coarsened
  integer :: max_nnz_i_matrix ! maximum possible number non zeros in i_matrix
  integer :: ptr       ! pointers to next position to be filled in the matrix
  integer :: point_i   ! general purpose integers for looping
  integer :: point_j   ! through matrix rows, etc
  integer :: index_j
  integer :: row_start
  integer :: row_end
  integer :: alloc_stat

  ! parameters used to calculate interpolation weights
  real(kind=myreal) :: alpha     ! stores an intermediate result
  real(kind=myreal) :: i_weight  ! current interpolation weight
  real(kind=myreal) :: sum_C     ! sum of values associated with the C points
  real(kind=myreal) :: sum_non_C ! sum of negative values associated with non C
                                 ! points
  real(kind=myreal) :: sum_pos   ! sum of positive values
  real(kind=myreal) :: value     ! current value from matrix

  ! parameters used in interpolation truncation
  real(kind=myreal) :: sum_weights       ! sum of the weights in current row
  real(kind=myreal) :: max_weight        ! maximum weight in current row
  real(kind=myreal) :: sum_trunc_weights ! sum of the truncated set of weights
  real(kind=myreal) :: trunc_threshold   ! defines which weights are removed
  real(kind=myreal) :: scaling           ! scaling applied to weights

  ! find the size of the matrix
  size_matrix = matrix%m

  ! set i_matrix parameters

  ! set number of rows
  i_matrix%m = size_matrix

  ! set number of columns, which is the number of C points
  i_matrix%n = count(point_type(1:size_matrix) .gt. 0)

  ! set type
  call zd11_put(i_matrix%type, "general", alloc_stat)
  if (alloc_stat .ne. 0) go to 10

  ! find maximum possible number of non zeros in i_matrix

  ! first find non-zeros resulting from C points, i.e. one non-zero per C
  max_nnz_i_matrix = i_matrix%n

  ! then run through matrix and find the non-zeros associated with F points
  ! i.e. for each F point there is one non-zero for
  ! each strong connection to a C point
  do point_i = 1,size_matrix
   ! find rows associated with F points
   if (point_type(point_i) == -1) then
    row_start = matrix%ptr(point_i)
    row_end = matrix%ptr(point_i+1)-1
    ! run through the row and add to nnz_i_matrix for
    ! every strong connection
    do index_j = row_start+1, row_end
     point_j = matrix%col(index_j)
     if ( (point_type(point_j) .gt. 0) .and. &
            (matrix%val(index_j) .le. thresholds(point_i)) ) then
      max_nnz_i_matrix = max_nnz_i_matrix + 1
     end if
    end do
   end if
  end do

  ! allocate col, val and ptr to store the interpolation matrix - if
  ! interpolation is to be truncated use temporary arrays for val and col,
  ! otherwise we can store values in i_matrix directly
  allocate(i_matrix%ptr(size_matrix+1), stat = alloc_stat)
  if (alloc_stat .ne. 0) go to 10

  if (control%trunc_parameter .gt. zero) then
   allocate(val(max_nnz_i_matrix), stat = alloc_stat)
   if (alloc_stat .ne. 0) go to 10

   allocate(col(max_nnz_i_matrix), stat = alloc_stat)
   if (alloc_stat .ne. 0) go to 10
  else
   allocate(i_matrix%val(max_nnz_i_matrix), stat = alloc_stat)
   if (alloc_stat .ne. 0) go to 10

   allocate(i_matrix%col(max_nnz_i_matrix), stat = alloc_stat)
   if (alloc_stat .ne. 0) go to 10
  end if

  ! run through the points and fill in values for
  ! i_matrix%ptr, i_matrix%col and i_matrix%val
  ptr = 1
  do point_i = 1,size_matrix
  ! set value of ptr for this row
  i_matrix%ptr(point_i) = ptr

  ! for C points put a 1 on the diagonal
  if ( point_type(point_i) .gt. 0 ) then
   if (control%trunc_parameter .gt. zero) then
    val(ptr) = one
    col(ptr) = point_type(point_i)
   else
    i_matrix%val(ptr) = one
    i_matrix%col(ptr) = point_type(point_i)
   end if

   ! update pointer
   ptr = ptr + 1

  ! for F points - generate the interpolation weights
  else if ( point_type(point_i) == -1 ) then

   ! for interpolation truncation must store sum of interpolation weights and
   ! maximum interpolation weight
   if (control%trunc_parameter .gt. zero) then
    sum_weights = zero
    max_weight = zero
   end if

   ! run through corresponding row in matrix and find the sums necessary to
   ! calculate the interpolation weights
   sum_C = zero
   sum_non_C = zero
   sum_pos = zero
   row_start = matrix%ptr(point_i)
   row_end = matrix%ptr(point_i+1)-1
   do index_j = row_start, row_end
    point_j = matrix%col(index_j)
    value = matrix%val(index_j)
    ! C points
    if ( (point_type(point_j) .gt. 0) .and. (value .le. thresholds(point_i)) ) &
     then
     sum_C = sum_C + value
     ! negative values
    else if ( value < 0 ) then
     sum_non_C = sum_non_C + value
     ! positive values
    else if ( value > 0 ) then
     sum_pos = sum_pos + value
    end if
   end do

   ! find alpha
   alpha = -( sum_C + sum_non_C ) / ( sum_C * sum_pos )

   ! run through C points, calculate interpolation weights and store
   do index_j = row_start+1, row_end
    point_j = matrix%col(index_j)
    value = matrix%val(index_j)
    ! bypass weak connections
    if ( (point_type(point_j) .gt. 0 ) .and. &
              ( value .le. thresholds(point_i)) ) then
     ! find interpolation weight
     i_weight = value * alpha

     ! put values in val and col
     if (control%trunc_parameter .gt. zero) then
      val(ptr) = i_weight
      col(ptr) = point_type(point_j)
      sum_weights = sum_weights + i_weight
      if (i_weight .gt. max_weight) then
       max_weight = i_weight
      end if
     else
      i_matrix%val(ptr) = i_weight
      i_matrix%col(ptr) = point_type(point_j)
     end if

     ! update pointer
     ptr = ptr + 1
    end if
   end do

   ! for interpolation truncation run through the values and truncate
   if (control%trunc_parameter .gt. zero) then
    ! find truncation threshold
    trunc_threshold = max_weight * control%trunc_parameter

    ! keep track of the sum of the truncated row sum
    sum_trunc_weights = zero

    ! find start and end of this row
    row_start = i_matrix%ptr(point_i)
    row_end = ptr-1

    ! set ptr back to start of the row
    ptr = row_start

    ! run through the values and remove small weights
    do index_j = row_start, row_end
     ! get the weight
     i_weight = val(index_j)

     ! if the weight is above threshold, move it to the new position in
     ! temp_matrix and update nnz_i_matrix
     if (i_weight .gt. trunc_threshold) then
      sum_trunc_weights = sum_trunc_weights + i_weight
      val(ptr) = val(index_j)
      ptr = ptr + 1
     end if
    end do

    ! find scaling
    scaling = sum_weights/sum_trunc_weights

    ! run through the remaining values and scale weights
    row_end = ptr-1
    do index_j = row_start, row_end
    val(index_j) = val(index_j) * scaling
    end do
   end if
  end if
 end do

 ! For interpolation truncation, generate col and val in i_matrix, copy across
 ! values into i_matrix and delete the temporary arrays
 if (control%trunc_parameter .gt. zero) then

  allocate(i_matrix%val(ptr-1), stat = alloc_stat)
  if (alloc_stat .ne. 0) go to 10

  i_matrix%val = val(1:ptr-1)

  deallocate(val, stat=alloc_stat)
  if (alloc_stat .ne. 0) then
   info%stat = alloc_stat; info%flag=-11
  end if

  allocate(i_matrix%col(ptr-1), stat = alloc_stat)
  if (alloc_stat .ne. 0) go to 10

  i_matrix%col = col(1:ptr-1)

  deallocate(col, stat=alloc_stat)
  if (alloc_stat .ne. 0) then
   info%stat = alloc_stat; info%flag=-11
  end if

 end if

 ! set the last value of i_matrix%ptr
 i_matrix%ptr(size_matrix+1) = ptr

 return

  10 info%stat = alloc_stat
     info%flag = -10

 end subroutine gen_interpolation


 ! -----------------------
 ! -- gen_coarse_matrix --
 ! -----------------------

 ! Subroutine to generate the coarse level coefficient matrix

 subroutine gen_coarse_matrix(matrix, i_matrix, c_matrix, info)

  type(zd11_type), intent(in) :: matrix    ! matrix being coarsened
  type(zd11_type), intent(in) :: i_matrix  ! interpolation matrix
  type(zd11_type), intent(out) :: c_matrix ! new coarse level matrix
  type(mi20_info), intent(inout) :: info

  integer :: MC65_status        ! error flag for MC65 subroutines
  type(zd11_type) :: temp_zd11  ! hold an intermediate result
  type(zd11_type) :: i_matrix_t ! transpose of interpolation matrix

  ! Find temp_zd11 = matrix x i_matrix
  call MC65_matrix_multiply(matrix, i_matrix, temp_zd11, MC65_status, &
        stat=info%stat)
  if (MC65_status == MC65_ERR_MEMORY_ALLOC) then
   info%flag = -10
   return
  else if (MC65_status == MC65_ERR_MEMORY_DEALLOC) then
   info%flag = -11
   return
  end if

  ! Get i_matrix transpose
  call MC65_matrix_transpose(i_matrix, i_matrix_t, MC65_status, &
        stat=info%stat)
  if (MC65_status == MC65_ERR_MEMORY_ALLOC) then
   info%flag = -10
   return
  else if (MC65_status == MC65_ERR_MEMORY_DEALLOC) then
   info%flag = -11
   return
  end if

  ! Find matrix_c = i_matrix_t x temp_zd11
  call MC65_matrix_multiply(i_matrix_t, temp_zd11, c_matrix, MC65_status, &
        stat=info%stat)
  if (MC65_status == MC65_ERR_MEMORY_ALLOC) then
   info%flag = -10
   return
  else if (MC65_status == MC65_ERR_MEMORY_DEALLOC) then
   info%flag = -11
   return
  end if

  ! Deallocate i_matrix_t and temp_mat
  call MC65_matrix_destruct(temp_zd11, MC65_status, stat=info%stat)
  if (MC65_status == MC65_ERR_MEMORY_DEALLOC) then
   info%flag = -11
   return
  end if

  call MC65_matrix_destruct(i_matrix_t, MC65_status, stat=info%stat)
  if (MC65_status == MC65_ERR_MEMORY_DEALLOC) then
   info%flag = -11
  end if

 end subroutine gen_coarse_matrix

! ---------------------------
! -- mi20_precondition_noA --
! -- (T. Rees, April 2015) --
! ---------------------------

 subroutine mi20_precondition_noA(coarse_data, rhs, solution, keep, &
      control, info, ma48_cntrl)

   ! coarse level data generated using mi20_setup
   type(mi20_data), dimension(:), intent(in) :: coarse_data

   ! the vector on which to apply the preconditioner
   real(kind=myreal), dimension(:), intent(in) :: rhs

   ! the vector after preconditioning
   real(kind=myreal), dimension(:), intent(out) :: solution

   ! derived type containing mi20 internal values and arrays
   type(mi20_keep), intent(inout) :: keep

   ! control parameters
   type(mi20_control), intent(in) :: control

   ! derived type to pass information to the calling program
   type(mi20_info), intent(inout) :: info

   ! optional user supplied derived type to control ma48
   type(ma48_control), optional, intent(inout) :: ma48_cntrl


   !------------------
   ! run the code...

   if ( keep%zd11_internal_conversion ) then
         ! pull A from keep, and send back 
      call mi20_precondition_withA(keep%A,coarse_data, rhs, solution, & 
           keep, control, info, ma48_cntrl)
   else
      if ( (control%print_level .gt. 0) .and. (control%error .ge. 0) ) then
         write(control%error,'(a)') 'mi20_precondition must be called with', & 
              ' ''A'' if mi20_setup was called'
      end if
      info%flag = -20
      ! flag an error
   end if

 end subroutine mi20_precondition_noA

 ! -----------------------------
 ! -- mi20_precondition_withA --
 ! -----------------------------

 ! A subroutine to perform the preconditioning operation

 subroutine mi20_precondition_withA(matrix, coarse_data, rhs, solution, keep, &
      control, info, ma48_cntrl)

   ! the matrix used to generate the preconditioner
   type(zd11_type), intent(in) :: matrix

   ! coarse level data generated using mi20_setup
   type(mi20_data), dimension(:), intent(in) :: coarse_data

   ! the vector on which to apply the preconditioner
   real(kind=myreal), dimension(:), intent(in) :: rhs

   ! the vector after preconditioning
   real(kind=myreal), dimension(:), intent(out) :: solution

   ! derived type containing mi20 internal values and arrays
   type(mi20_keep), intent(inout) :: keep

   ! control parameters
   type(mi20_control), intent(in) :: control

   ! derived type to pass information to the calling program
   type(mi20_info), intent(inout) :: info

   ! optional user supplied derived type to control ma48
   type(ma48_control), optional, intent(inout) :: ma48_cntrl

   ! vector storage required for the preconditioner
   type(precon_vec_storage), dimension(:), allocatable :: vec_storage

   real(kind=myreal) :: rn    ! norm of residual
   real(kind=myreal) :: rn0   ! initial residual norm
   real(kind=myreal) :: val   ! a value in matrix

   integer :: level           ! current level
   integer :: vec_length      ! vector length on current level
   integer :: cnt             ! current iteration
   integer :: MC65_status     ! error flag for MC65 subroutines
   integer :: max_level       ! max number of levels in preconditioner
   integer :: i               ! general purpose integer
   integer :: size_clevel     ! size of coarse level matrix
   integer :: row
   integer :: col
   integer :: row_start
   integer :: row_end
   integer :: length
   integer :: alloc_stat

   logical :: deallocate_ma48_matrix ! .true. if ma48_matrix in keep is to be
   ! deallocated
   logical :: error_out        ! .true. if outputting errors
   logical :: print_out        ! .true. if outputting diagnostics
   logical :: delete_csolver   ! .true. if need to delete coarse solver data

   ! reset info%flag
   info%flag = 0

   ! switch on/off output and error messages
   ! ---------------------------------------
   error_out = .false.
   if ( (control%print_level .gt. 0) .and. (control%error .ge. 0) ) &
        error_out = .true.

   print_out = .false.
   if ( (control%print_level .gt. 1) .and. (control%print .ge. 0) ) &
        print_out = .true.

   ! check setup was sucessful
   if (keep%clevels == -1) then
      info%flag = -15
      if (error_out) &
           call error_message_precondition(control, info)
      return
   end if

   ! check control parameters are in range
   ! -------------------------------------
   call test_control_precondition(control, info)
   if (info%flag < 0) return

   ! deal with the case when matrix%m == 1
   ! -------------------------------------
   if (matrix%m .eq. 1) then
      solution(1) = rhs(1) / matrix%val(1)
      return
   end if


   ! check vectors are long enough
   ! -----------------------------
   if ( (size(rhs)<matrix%n) .or. (size(solution)<matrix%n) ) then
      info%flag = -16
      if (error_out) call error_message_precondition(control, info)
      return
   end if

   ! perform other checks
   ! --------------------

   ! check if ma48_cntrl present (if required for ma48 solver)
   if ( (control%coarse_solver .eq. 3) .and. (control%ma48 .gt. 0) ) then
      if (.not. present(ma48_cntrl) ) then
         info%flag = -19
         if (error_out) call error_message_precondition(control, info)
         return
      end if
   end if

   ! find the value of max_level (ie the level of the coarsest matrix used)
   ! ---------------------------
   if ( (control%levels .lt. 0) .or. (keep%clevels .lt. control%levels) ) then
      max_level = keep%clevels
      if (keep%clevels < control%levels) then
         info%flag = 20
         if (error_out) call error_message_precondition(control, info)
      end if
   else
      max_level = control%levels
   end if

   ! check if direct solver data can be deleted
   ! ------------------------------------------
   delete_csolver = .false.
   deallocate_ma48_matrix = .false.
   ! check 1 - using an iterative coarse solver
   if (control%coarse_solver .le. 2) then
      delete_csolver = .true.
      deallocate_ma48_matrix = .true.
      ! check 2 - change of max_level
   else if (max_level .ne. keep%dsolve_level) then
      ! force analyze in hsl_ma48 if levels have changed
      if (keep%ma48_data .gt. 1) then
         keep%ma48_data = 1
      end if
      delete_csolver = .true.
      deallocate_ma48_matrix = .true.
      ! check 3 - using lapack solver and hsl_ma48 data exists or
   else if ((control%coarse_solver .eq. 4) .and. ( keep%ma48_data .eq. 3)) then
      delete_csolver = .true.
      deallocate_ma48_matrix = .true.
      ! check 4 - using hsl_ma48 solver and lapack data exists
   else if ( (control%coarse_solver .eq. 3) .and. (keep%lapack_data) )then
      delete_csolver = .true.
      ! check 5 - using hsl_ma48 solver and the calling program requests
      ! initialize, analyse or factorize
   else if ( (control%coarse_solver .eq. 3) .and. &
        (control%ma48 .gt. 0 ) .and. (control%ma48 .lt. 4) )then
      delete_csolver = .true.
      ! check 6 - new preconditoner
   else if (keep%new_preconditioner) then
      delete_csolver = .true.
      deallocate_ma48_matrix = .true.
   end if

   ! reset necessary values in keep if this is a new preconditioner
   ! --------------------------------------------------------------
   if (keep%new_preconditioner) then
      keep%lapack_data = .false.
      keep%ma48_data = 0
      call ma48_initialize(keep%ma48_factors,keep%ma48_cntrl)
      keep%ma48_cntrl%lp = control%error
      keep%ma48_cntrl%wp = control%error
      keep%ma48_cntrl%mp = control%print
      select case (control%print_level)
      case (0)
         keep%ma48_cntrl%ldiag = 0
      case (1, 2)
         keep%ma48_cntrl%ldiag = 2
      end select
      keep%ma48_matrix_exists = .false.
      keep%dsolve_level = -1
      keep%new_preconditioner = .false.
   end if

   ! delete direct solver data if necessary
   ! --------------------------------------
   if (delete_csolver) then
      if (control%ma48 .eq. 0) then
         call dsolver_deallocation(keep, control, info, deallocate_ma48_matrix)
      else
         call dsolver_deallocation(keep, control, info, &
              deallocate_ma48_matrix, ma48_cntrl)
      end if
      ! deal with errors
      if ( info%flag < 0 ) then
         if (error_out) call error_message_precondition(control, info)
         return
      end if
   end if

   ! output for solver only
   ! ----------------------
   !!  if ( print_out .and. control%tol > 0 ) &
   !!   write(control%print,*) "mi20_solver called"


   ! Generate direct solver data (if required)
   ! -----------------------------------------
   ! find size of coarse level matrix
   if (max_level == 0) then
      size_clevel = matrix%m
   else
      size_clevel = coarse_data(max_level)%A_mat%m
   end if

   if (size_clevel .gt. 1) then

      ! LAPACK direct solver setup
      ! --------------------------
      if (control%coarse_solver .eq. 4) then

         ! check if new LAPACK data is required
         if (.not. keep%lapack_data) then

            ! allocate new storage
            allocate(keep%lapack_factors(size_clevel,size_clevel), &
               stat=alloc_stat)
            if (alloc_stat .ne. 0) then
               info%stat = alloc_stat; info%flag = -10
               if (error_out) call error_message_precondition(control, info)
               return
            end if

            allocate(keep%lapack_pivots(size_clevel), stat=alloc_stat)
            if (alloc_stat .ne. 0) then
               info%flag = -10; info%stat = alloc_stat
               if (error_out) call error_message_precondition(control, info)
               return
            end if

            ! copy the max_level matix into keep%lapack_factors
            keep%lapack_factors = 0
            if (max_level == 0) then
               do row = 1, size_clevel
                  row_start = matrix%ptr(row)
                  row_end = matrix%ptr(row+1)-1
                  do i = row_start, row_end
                     col = matrix%col(i)
                     val = matrix%val(i)
                     keep%lapack_factors(row, col) = val
                  end do
               end do
            else
               do row = 1, size_clevel
                  row_start = coarse_data(max_level)%A_mat%ptr(row)
                  row_end = coarse_data(max_level)%A_mat%ptr(row+1)-1
                  do i = row_start, row_end
                     col = coarse_data(max_level)%A_mat%col(i)
                     val = coarse_data(max_level)%A_mat%val(i)
                     keep%lapack_factors(row, col) = val
                  end do
               end do
            end if

            ! call lapack factorization subroutine
            call dgetrf(size_clevel, size_clevel, keep%lapack_factors, &
                 size_clevel, keep%lapack_pivots, info%getrf_info)

            ! deal with failure
            if (info%getrf_info .ne. 0) then
               info%flag = -17
               if (error_out) call error_message_precondition(control, info)
               return
            end if

            keep%lapack_data = .true.
            keep%dsolve_level = max_level
         end if
      end if

      ! hsl_ma48 direct solver setup
      ! ----------------------------
      if (control%coarse_solver .eq. 3) then


         ! check for the following:
         ! - factors don't exist
         ! or
         ! - the calling program requires initialize, analyse or factorize
         if ( (keep%ma48_data .ne. 3) .or. &
              ((control%ma48 .gt. 0 ).and.(control%ma48 .lt. 4 )) ) then

            ! first - check if need call ma48_initialize (and update
            ! keep%ma48_data)
            if ((keep%ma48_data .eq. 0) .or. (control%ma48 .eq. 1)) then

               ! if user supplied ma48_cntrl requires initializing then do so
               if ( (present(ma48_cntrl)) .and. (keep%ma48_data .eq. 0) &
                    .or. (control%ma48 .eq. 1) ) then
                  call ma48_initialize(keep%ma48_factors, ma48_cntrl)
                  ma48_cntrl%lp = control%error
                  ma48_cntrl%wp = control%error
                  ma48_cntrl%mp = control%print
                  select case (control%print_level)
                  case (0)
                     ma48_cntrl%ldiag = 0
                  case (1, 2)
                     ma48_cntrl%ldiag = 2
                  end select
               end if

               ! update keep%ma48_data
               if (keep%ma48_data .eq. 0) keep%ma48_data = 1

               ! if requested then return
               if (control%ma48 .eq. 1) return
            end if

            ! second - check if we need to analyse, i.e. analyse not been
            ! performed or an analyse is requested by control%ma48
            if ((keep%ma48_data .eq. 1) .or. (control%ma48 .eq. 2))  then

               ! if ma48_matrix does not exist then allocate storage and copy in
               ! values from the max_level matrix (ma48 requires a different
               ! sparse storage format)
               if (.not. keep%ma48_matrix_exists) then
                  if (max_level == 0) then
                     length = matrix%ptr(matrix%n+1)-1
                     keep%ma48_matrix%ne = length
                     keep%ma48_matrix%m = matrix%m
                     keep%ma48_matrix%n = matrix%n
                  else
                     length = coarse_data(max_level)%A_mat% &
                          ptr(coarse_data(max_level)%A_mat%n+1)-1
                     keep%ma48_matrix%ne = length
                     keep%ma48_matrix%m = coarse_data(max_level)%A_mat%m
                     keep%ma48_matrix%n = coarse_data(max_level)%A_mat%n
                  end if

                  ! allocate arrays
                  allocate(keep%ma48_matrix%col(length), stat=alloc_stat)
                  if (alloc_stat .ne. 0) go to 10

                  allocate(keep%ma48_matrix%row(length), stat=alloc_stat)
                  if (alloc_stat .ne. 0) go to 10

                  allocate(keep%ma48_matrix%val(length), stat=alloc_stat)
                  if (alloc_stat .ne. 0) go to 10

                  ! fill arrays
                  if (max_level == 0) then
                     ! The two loops below have been changed from array slice
                     ! notation to do loops due to a reported issue with the
                     ! intel compiler.
                     do i = 1, length
                        keep%ma48_matrix%val(i) = matrix%val(i)
                     end do
                     do i = 1, length
                        keep%ma48_matrix%col(i) = matrix%col(i)
                     end do

                     do row = 1, matrix%m
                        row_start = matrix%ptr(row)
                        row_end = matrix%ptr(row+1)-1
                        keep%ma48_matrix%row(row_start:row_end) = row
                     end do

                  else
                     ! The two loops below have been changed from array slice
                     ! notation to do loops due to a reported issue with the
                     ! intel compiler.
                     do i = 1, length
                        keep%ma48_matrix%val(i) = &
                             coarse_data(max_level)%A_mat%val(i)
                     end do
                     do i = 1, length
                        keep%ma48_matrix%col(i) = &
                             coarse_data(max_level)%A_mat%col(i)
                     end do

                     do row = 1, coarse_data(max_level)%A_mat%m
                        row_start = coarse_data(max_level)%A_mat%ptr(row)
                        row_end = coarse_data(max_level)%A_mat%ptr(row+1)-1
                        keep%ma48_matrix%row(row_start:row_end) = row
                     end do
                  end if

               end if

               ! perform analyse (and optional factorize)
               if (control%ma48 .eq. 0) then
                  call ma48_analyse(keep%ma48_matrix, keep%ma48_factors, &
                       keep%ma48_cntrl, info%ma48_ainfo, info%ma48_finfo)
               else if (control%ma48 .ge. 3) then
                  call ma48_analyse(keep%ma48_matrix, keep%ma48_factors, &
                       ma48_cntrl, info%ma48_ainfo, info%ma48_finfo)
               else if (control%ma48 .eq. 2) then
                  call ma48_analyse(keep%ma48_matrix, keep%ma48_factors, &
                       ma48_cntrl, info%ma48_ainfo)
               end if

               ! deal with error
               if (info%ma48_ainfo%flag .lt. 0) then
                  info%flag = -18
                  if (error_out) call error_message_precondition(control, info)
                  return
               end if

               ! update keep%ma48_data
               if ( (control%ma48 .eq. 0) .or. (control%ma48 .ge. 3) ) then
                  keep%ma48_data = 3
                  keep%dsolve_level = max_level
               else if (control%ma48 .eq. 2) then
                  keep%ma48_data = 2
               end if

               ! check if need to return to calling program
               if ( (control%ma48 .eq. 2) .or. (control%ma48 .eq. 3) ) return

            end if

            ! third - factorize if required and deal with factorize errors
            if (keep%ma48_data .eq. 2) then
               call ma48_factorize(keep%ma48_matrix, keep%ma48_factors, &
                    ma48_cntrl, info%ma48_finfo)

               ! deal with error
               if (info%ma48_finfo%flag .lt. 0) then
                  info%flag = -18
                  if (error_out) call error_message_precondition(control, info)
                  return
               end if

               keep%ma48_data = 3
               keep%dsolve_level = max_level

               ! return to calling program if required
               if (control%ma48 .eq. 3) return

            end if

         end if
      end if

   end if ! end of setting up direct solvers


   ! Generate vector storage
   ! -----------------------
   allocate(vec_storage(0:max_level), stat=alloc_stat)
   if (alloc_stat .ne. 0) go to 10

   do level = 0,max_level
      if (level == 0) then
         vec_length = matrix%m
      else
         vec_length = coarse_data(level)%A_mat%m
      end if

      allocate(vec_storage(level)%vec_sol (vec_length), stat=alloc_stat)
      if (alloc_stat .ne. 0) go to 10
      allocate(vec_storage(level)%vec_er (vec_length), stat=alloc_stat)
      if (alloc_stat .ne. 0) go to 10
      allocate(vec_storage(level)%vec_rhs (vec_length), stat=alloc_stat)
      if (alloc_stat .ne. 0) go to 10
   end do

10 if (alloc_stat .ne. 0) then
      info%stat = alloc_stat; info%flag = -10
      if (error_out) call error_message_precondition(control, info)
      return
   end if

   ! set initial values
   ! ------------------
   vec_length = matrix%n

   vec_storage(0)%vec_rhs(1:vec_length) = rhs(1:vec_length)

   ! find starting residual norm rn
   rn = 0
   do i = 1,vec_length
      rn = rn + rhs(i)**2
   end do
   rn = sqrt(rn)
   rn0 = rn

   cnt = 1

   ! iteration loop
   ! --------------
   do

      ! check for error from excessive divergenve
      if (rn > rn0 * control%err_tol) then
         info%flag = -14
         if (error_out) call error_message_precondition(control, info)
         call precon_vec_deallocation(vec_storage, max_level, info)
         return
      end if

      !!   ! check for termination against convergence
      !!   ! (not needed for preconditioner)
      !!   if (control%tol>0) then
      !!    if (control%tol_relative) then
      !!     if (rn < control%tol * rn0) exit
      !!    else
      !!     if (rn < control%tol) exit
      !!    end if
      !!   end if

      ! check for convergence against number of iterations
      if (cnt > control%v_iterations) exit

      ! perform 1 v-cycle
      if (control%coarse_solver .eq. 3) then
         call amg_v_cycle(matrix, coarse_data, max_level, vec_storage, &
              keep, control, info, ma48_cntrl)
      else
         call amg_v_cycle(matrix, coarse_data, max_level, vec_storage, &
              keep, control, info)
      end if

      ! deal with errors
      if (info%flag < 0) then
         if (error_out) call error_message_precondition(control, info)
         call precon_vec_deallocation(vec_storage, max_level, info)
         return
      end if

      ! update solution
      if  (cnt == 1) then
         solution(1:vec_length) = vec_storage(0)%vec_sol(1:vec_length)
      else
         solution(1:vec_length) = solution(1:vec_length) &
       + vec_storage(0)%vec_sol(1:vec_length)
      end if

      ! find new residual and store in vec_storage
      call MC65_matrix_multiply_vector(matrix, solution, &
           vec_storage(0)%vec_rhs, MC65_status)

      vec_storage(0)%vec_rhs(1:vec_length) = &
         rhs(1:vec_length) - vec_storage(0)%vec_rhs(1:vec_length)

      ! find new residual norm
      rn = 0
      do i = 1,vec_length
         rn = rn + vec_storage(0)%vec_rhs(i)**2
      end do
      rn = sqrt(rn)
      ! for solver only
      !!   if ((print_out).and.(control%tol>0)) &
      !!    write(control%print,*) "after v-cycle: vector norm=", rn

      cnt = cnt + 1
   end do

   ! deallocation
   ! ------------
   call precon_vec_deallocation(vec_storage, max_level, info)
   if ( (info%flag < 0) .and. (error_out) ) then
      call error_message_precondition(control, info)
   end if

 end subroutine mi20_precondition_withA


 ! --------------------------
 ! -- dsolver_deallocation --
 ! --------------------------

 ! Subroutine to free up memory allocated for the direct solvers

 subroutine dsolver_deallocation(keep, control, info, &
      deallocate_ma48_matrix_in, ma48_cntrl)

  type(mi20_keep), intent(inout) :: keep
  type(mi20_control), intent(in) :: control
  type(mi20_info), intent(inout) :: info
  logical, intent(in), optional :: deallocate_ma48_matrix_in
  type(ma48_control), optional, intent(inout) :: ma48_cntrl

  logical :: deallocate_ma48_matrix ! .true. if ma48_matrix
                                    ! is to be deallocated
  integer :: alloc_stat

  if (present(deallocate_ma48_matrix_in)) then
   deallocate_ma48_matrix = deallocate_ma48_matrix_in
  else
   deallocate_ma48_matrix = .true.
  end if

  ! deallocate memory allocated for lapack solver
  ! ---------------------------------------------
  if (allocated(keep%lapack_factors)) then
   deallocate(keep%lapack_factors, stat=alloc_stat)
   if (alloc_stat .ne. 0) then
    info%stat = alloc_stat; info%flag=-11
   end if
  end if

  if (allocated(keep%lapack_pivots)) then
   deallocate(keep%lapack_pivots, stat=alloc_stat)
   if (alloc_stat .ne. 0) then
    info%stat = alloc_stat; info%flag=-11
   end if
  end if
  keep%lapack_data = .false.

  ! deallocate date associated with ma48
  ! ------------------------------------
  if ( (control%ma48 .gt. 0) .and. present(ma48_cntrl) ) then
   call ma48_finalize(keep%ma48_factors, ma48_cntrl, alloc_stat)
  else
   call ma48_finalize(keep%ma48_factors, keep%ma48_cntrl, alloc_stat)
  end if
  if (alloc_stat .ne. 0) then
   info%stat = alloc_stat; info%flag = -11
  end if
  if (keep%ma48_data == 3) then
   keep%ma48_data = 2
  end if

  ! deallocate ma48_matrix if required
  if (deallocate_ma48_matrix) then
   call ma48_matrix_deallocation(keep%ma48_matrix, info)
   keep%ma48_matrix_exists = .false.
  end if

 end subroutine dsolver_deallocation


 ! ------------------------------
 ! -- ma48_matrix_deallocation --
 ! ------------------------------

 ! Subroutine to free up memory allocated for the ma48_matrix_type required by
 ! ma48 direct solver

 subroutine ma48_matrix_deallocation(matrix, info)

  type(zd11_type), intent(inout) :: matrix
  type(mi20_info), intent(inout) :: info

  integer :: alloc_stat

  if (allocated(matrix%row)) then
   deallocate(matrix%row, stat=alloc_stat)
   if (alloc_stat .ne. 0) go to 10
  end if
  if (allocated(matrix%col)) then
   deallocate(matrix%col, stat=alloc_stat)
   if (alloc_stat .ne. 0) go to 10
  end if
  if (allocated(matrix%val)) then
   deallocate(matrix%val, stat=alloc_stat)
   if (alloc_stat .ne. 0) go to 10
  end if
  return

10 info%stat = alloc_stat
   info%flag = -11

 end subroutine ma48_matrix_deallocation



 ! -----------------------------
 ! -- precon_vec_deallocation --
 ! -----------------------------

 ! Subroutine to free up memory allocated in a precon_vec_storage
 ! array
 ! vec_storage contains the allocatable arrays to be deallocated and
 ! max_level gives the number of levels allocated

 subroutine precon_vec_deallocation(vec_storage, max_level, info)

  type(precon_vec_storage), dimension(:), allocatable, intent(inout) :: &
       vec_storage
  integer, intent(in) :: max_level
  type(mi20_info), intent(inout) :: info

  integer :: level, alloc_stat

  ! deallocate memory
  if (allocated(vec_storage)) then
   do level = 0,max_level
    if (allocated(vec_storage(level)%vec_sol)) then
     deallocate(vec_storage(level)%vec_sol, stat=alloc_stat)
     if (alloc_stat .ne. 0) go to 10
    end if
    if (allocated(vec_storage(level)%vec_rhs)) then
     deallocate(vec_storage(level)%vec_rhs, stat=alloc_stat)
     if (alloc_stat .ne. 0) go to 10
    end if
    if (allocated(vec_storage(level)%vec_er)) then
     deallocate(vec_storage(level)%vec_er, stat=alloc_stat)
     if (alloc_stat .ne. 0) go to 10
    end if
   end do
   deallocate(vec_storage, stat=alloc_stat)
   if (alloc_stat .ne. 0) go to 10
  end if
  return

 10 info%stat = alloc_stat
    info%flag=-11

 end subroutine precon_vec_deallocation

 ! -----------------
 ! -- amg_v_cycle --
 ! -----------------

 ! Subroutine to perform one amg v-cycle

 ! Given a system Ax=b to be solved, level k of the v-cycle attempts to find
 ! the error for an approximate solution y by solving Ae=r, where r=b-Ay. The
 ! steps at each level are as follows:
 ! 1. presmooth y using iterative solver, e.g. Gauss-Seidel
 ! 2. calculate the residual
 ! 4. restrict residual to next coarse level
 ! 5. perform coarse level solve on Ae=r to give approximation to e, using
 !    either:
 !    -direct solver
 !    -simple iterative method, e.g. Gauss-Seidel or
 !    -next multigrid level
 ! 6. prolong error e back to previous level
 ! 7. update solution y=y+e
 ! 8. postsmooth solution vector

 subroutine amg_v_cycle(matrix, coarse_data, max_level, vec_storage, &
      keep, control, info, ma48_cntrl)

   ! matrix A used to generate preconditioner
   type(zd11_type), intent(in) :: matrix

   ! coarse level data generated in mi20_setup
   type(mi20_data), dimension(:), intent(in) :: coarse_data

   ! the last coarse level used in the preconditioner
   integer, intent(in) :: max_level

   ! for storage of v-cycle data
   type(precon_vec_storage), dimension(0:), intent(inout) :: vec_storage
   ! use of vec_storage
   ! ------------------
   ! vec_sol: this holds the solution for level k (and solution at level k is
   ! equal to the error at level k-1)
   ! vec_er: stores error or residual at level k
   ! vec_rhs: the rhs at level k (this is equal to the residual at level k-1)

   ! derived type to store internal data for the preconditioner
   type(mi20_keep), intent(inout) :: keep

   ! derived type storing control parameters
   type(mi20_control), intent(in) :: control

   ! derived type to pass information back to the calling program
   type(mi20_info), intent(inout) :: info

   ! optional derived type containing hsl_ma48 control parameters
   type(ma48_control), optional, intent(in) :: ma48_cntrl


   integer :: level         ! current level
   integer :: MC65_status   ! error flag for MC65 subroutines
   integer :: pre_sweeps    ! number of pre smoothing sweeps
   integer :: post_sweeps   ! number of post smoothing sweeps
   integer :: coarse_sweeps ! number of coarse level iterations (for iterative
   ! solvers only)
   integer :: vec_length    ! vector length on current level

   ! directions of pre and post smoothers: 1 for forward and -1 for reverse
   integer :: pre_dir      ! direction of pre smoothing
   integer :: post_dir     ! direction of post smoothing

   logical :: sol_is_zero  ! .true. if the solution is currently zero
   logical :: error_out    ! .true. if errors should be printed

   integer :: getrs_info   ! used in call to getrs to return solver information
   integer :: i            ! temporary variable
   
   error_out = .false.
   if ( (control%print_level .gt. 0) .and. (control%error .ge. 0) ) &
        error_out = .true.

   ! check that ma48_cntrl exists if required

   ! start at level 0
   level = 0

   ! set the direction of the smoother sweeps
   pre_dir = -1
   post_dir = 1

   ! find the number of sweeps
   pre_sweeps = control%pre_smoothing
   post_sweeps = control%post_smoothing
   coarse_sweeps = control%coarse_solver_its

   ! storage in vec_storage
   ! ----------------------
   ! vec_sol: solution for level k, i.e. error at level k-1
   ! vec_er: error/residual at level k
   ! vec_rhs: rhs at level k, i.e. residual at level k-1

   ! loop performing pre-smoothing and restriction of residual
   ! ---------------------------------------------------------
   do
      ! find length of vectors
      if (level == 0) then
         vec_length = matrix%n
      else
         vec_length = coarse_data(level)%A_mat%n
      end if

      ! check for termination of this loop
      if ( (vec_length == 1) .or. (level == max_level) ) exit

      ! set solution vector to zero
      vec_storage(level)%vec_sol = 0
      sol_is_zero = .true.

      ! perform pre-smoothing

      ! Gauss Seidel smoother
      if (control%smoother == 2) then
         if (level == 0) then
            call gs_smoother(matrix, vec_storage(level)%vec_rhs, &
                 pre_sweeps, pre_dir, vec_storage(level)%vec_sol, sol_is_zero)
         else
            call gs_smoother(coarse_data(level)%A_mat, &
                 vec_storage(level)%vec_rhs, pre_sweeps, pre_dir, &
                 vec_storage(level)%vec_sol, sol_is_zero)
         end if

         ! damped Jacobi smoother - use vec_er for temporary vector
      else if (control%smoother == 1) then
         if (level == 0) then
            call dj_smoother(matrix, vec_storage(level)%vec_rhs,       &
                 pre_sweeps, vec_storage(level)%vec_sol, sol_is_zero, &
                 vec_storage(level)%vec_er, control%damping)
         else
            call dj_smoother(coarse_data(level)%A_mat, &
                 vec_storage(level)%vec_rhs, pre_sweeps, &
                 vec_storage(level)%vec_sol, sol_is_zero, &
                 vec_storage(level)%vec_er, control%damping)
         end if
      end if

      ! find residual
      if (level == 0) then
         call MC65_matrix_multiply_vector( matrix, vec_storage(level)%vec_sol, &
              vec_storage(level)%vec_er, MC65_status)
      else
         call MC65_matrix_multiply_vector( coarse_data(level)%A_mat, &
              vec_storage(level)%vec_sol, vec_storage(level)%vec_er, &
              MC65_status)
      end if

      vec_storage(level)%vec_er = &
           vec_storage(level)%vec_rhs - vec_storage(level)%vec_er

      ! restrict residual
      call MC65_matrix_multiply_vector( coarse_data(level+1)%I_mat,  &
           vec_storage(level)%vec_er, vec_storage(level+1)%vec_rhs, &
           MC65_status, trans=.true.)

      ! go to next level
      level = level + 1
   end do

   ! coarse level solver
   ! -------------------
   ! for a vector length of 1 always perform a direct solve
   if (vec_length == 1) then
      vec_storage(level)%vec_sol = &
           vec_storage(level)%vec_rhs / coarse_data(level)%A_mat%val(1)
   else
      ! Gauss Seidel solver
      ! -------------------
      if (control%coarse_solver == 2) then

         ! set solution vector to zero
         vec_storage(level)%vec_sol = 0
         sol_is_zero = .true.

         if (level == 0) then
            call gs_smoother(matrix, vec_storage(level)%vec_rhs, &
                 coarse_sweeps, pre_dir, vec_storage(level)%vec_sol, &
                 sol_is_zero)
            call gs_smoother(matrix, vec_storage(level)%vec_rhs, &
                 coarse_sweeps, post_dir, vec_storage(level)%vec_sol, &
                 sol_is_zero)
         else
            call gs_smoother(coarse_data(level)%A_mat, &
                 vec_storage(level)%vec_rhs, coarse_sweeps, pre_dir, &
                 vec_storage(level)%vec_sol, sol_is_zero)
            call gs_smoother(coarse_data(level)%A_mat, &
                 vec_storage(level)%vec_rhs, coarse_sweeps, post_dir, &
                 vec_storage(level)%vec_sol, sol_is_zero)
         end if

         ! damped Jacobi solver  - use vec_er for temporary vector
         ! --------------------
      else if (control%coarse_solver == 1) then

         ! set solution vector to zero
         vec_storage(level)%vec_sol = 0
         sol_is_zero = .true.

         if (level == 0) then
            call dj_smoother(matrix, vec_storage(level)%vec_rhs, &
                 coarse_sweeps, vec_storage(level)%vec_sol, sol_is_zero, &
                 vec_storage(level)%vec_er, control%damping)
         else
            call dj_smoother(coarse_data(level)%A_mat, &
                 vec_storage(level)%vec_rhs, coarse_sweeps, &
                 vec_storage(level)%vec_sol, sol_is_zero, &
                 vec_storage(level)%vec_er, control%damping)
         end if

         ! LAPACK direct solver
         ! --------------------
      else if (control%coarse_solver .eq. 4) then

         vec_storage(level)%vec_sol = vec_storage(level)%vec_rhs
         call dgetrs('n', vec_length, 1, keep%lapack_factors, vec_length,  &
              keep%lapack_pivots, vec_storage(level)%vec_sol, vec_length, &
              getrs_info)

         ! ma48 direct solver
         ! ------------------
      else if (control%coarse_solver .eq. 3) then

         if ( (control%ma48 .gt. 0) .and. present(ma48_cntrl) ) then
            call ma48_solve(keep%ma48_matrix, keep%ma48_factors,          &
                 vec_storage(level)%vec_rhs, vec_storage(level)%vec_sol, &
                 ma48_cntrl, info%ma48_sinfo)
         else
            call ma48_solve(keep%ma48_matrix, keep%ma48_factors,          &
                 vec_storage(level)%vec_rhs, vec_storage(level)%vec_sol, &
                 keep%ma48_cntrl, info%ma48_sinfo)
         end if

         ! deal with error
         if (info%ma48_sinfo%flag .lt. 0) then
            info%flag = -18
            if (error_out) call error_message_precondition(control, info)
            return
         end if

      end if
   end if

   ! loop to perform prolongation of error and post smooth
   ! -----------------------------------------------------
   do
      if ( level == 0 ) exit

      level = level - 1

      ! prolong error
      call MC65_matrix_multiply_vector( coarse_data(level+1)%I_mat,  &
           vec_storage(level+1)%vec_sol, vec_storage(level)%vec_er, &
           MC65_status, trans=.false.)

      ! update solution
      ! JDH 10/04/2013: Limited length of vector to what is actually generated
      !                 in preceding call to avoid nagfor undef variable error
      i = coarse_data(level+1)%I_mat%m
      vec_storage(level)%vec_sol(1:i) = &
           vec_storage(level)%vec_sol(1:i) + vec_storage(level)%vec_er(1:i)

      ! post smoothing

      ! Gauss Seidel smoother
      if (control%smoother == 2) then
         if (level == 0) then
            call gs_smoother(matrix, vec_storage(level)%vec_rhs, &
                 post_sweeps, post_dir, vec_storage(level)%vec_sol, sol_is_zero)
         else
            call gs_smoother(coarse_data(level)%A_mat, &
                 vec_storage(level)%vec_rhs, post_sweeps, post_dir, &
                 vec_storage(level)%vec_sol, sol_is_zero)
         end if

         ! damped Jacobi smoother - use vec_er for temporary vector
      else if (control%smoother == 1) then
         if (level == 0) then
            call dj_smoother(matrix, vec_storage(level)%vec_rhs,        &
                 post_sweeps, vec_storage(level)%vec_sol, sol_is_zero, &
                 vec_storage(level)%vec_er, control%damping)
         else
            call dj_smoother(coarse_data(level)%A_mat, &
                 vec_storage(level)%vec_rhs, post_sweeps, &
                 vec_storage(level)%vec_sol, sol_is_zero, &
                 vec_storage(level)%vec_er, control%damping)
         end if
      end if
   end do

 end subroutine amg_v_cycle


 ! -----------------
 ! -- gs_smoother --
 ! -----------------

 ! Subroutine to perform Gauss_Seidel smoothing

 subroutine gs_smoother(matrix, rhs, sweeps, dir, vector, vector_is_zero)

  ! matrix for Gauss Seidel smoother
  type(zd11_type), intent(in) :: matrix

  ! right hand side vector for Gauss Seidel smoother
  real(kind=myreal), intent(in), dimension(:) :: rhs

  ! number of sweeps
  integer, intent(in) :: sweeps

  ! direction of sweep: forward if +1 & backward if -1
  integer, intent(in) :: dir

  ! vector to be smoothed
  real(kind=myreal), intent(inout), dimension(:) :: vector

  ! if vector_is_zero is .true. assume vector is zero on entry,
  ! and change vector_is_zero to .false. when vector becomes non-zero
  logical, intent(inout) :: vector_is_zero


  integer :: i, j         ! loop indices
  integer :: row          ! current row
  integer :: vec_length   ! length of vectors
  integer :: col          ! current column
  real(kind=myreal) :: sm ! some sum of values
  integer :: start, end   ! start and end for do loop through matrix rows
  integer :: step         ! step in do loop through matrix rows
  integer :: s, e         ! start and end location of current row

  vec_length = matrix%n

  if (dir > 0) then
   start = 1
   end = vec_length
   step = 1
  else
    start = vec_length
    end = 1
    step = -1
  end if

  do j=1,sweeps
   ! if vector is zero we can bypass some multiplications
   if (vector_is_zero) then
    ! run though rows
    do row = start, end, step

     ! find start and end of current row
     s = matrix%ptr(row)
     e = matrix%ptr(row+1)-1

     sm = rhs(row)

     ! run through the row
     do i = s+1, e
      ! get current column
      col = matrix%col(i)
      if (((dir>0) .and. (col<row)) .or. ((dir<0) .and. (col>row))) then
       sm = sm - vector(col) * matrix%val(i)
      end if
     end do
     ! update value in vector
     vector(row) = sm / matrix%val(s)
    end do
    vector_is_zero = .false.
   else
    ! run through the rows
    do row = start, end, step
     ! get start and ends of current row
     s = matrix%ptr(row)
     e = matrix%ptr(row+1)-1

     sm = rhs(row)

     ! run through the row
     do i = s+1, e
      col = matrix%col(i)
      sm = sm - vector(col) * matrix%val(i)
     end do

     ! update value in vector
     vector(row) = sm / matrix%val(s)
    end do
   end if
  end do

 end subroutine gs_smoother


 ! -----------------
 ! -- dj_smoother --
 ! -----------------

 ! Subroutine to perform damped Jacobi smoothing

 subroutine dj_smoother(matrix, rhs, sweeps, vector, vector_is_zero, &
                         new_vector, damping)

  ! matrix for damped Jacobi smoother
  type(zd11_type), intent(in) :: matrix

  ! right hand side vector for damped Jacobi smoother
  real(kind=myreal), intent(in), dimension(:) :: rhs

  ! number of sweeps to perform
  integer, intent(in) :: sweeps

  ! vector to be smoothed
  real(kind=myreal), intent(inout), dimension(:) :: vector

  ! if vector_is_zero is .true. assume vector is zero on entry,
  ! and change vector_is_zero to .false. when vector becomes non-zero
  logical, intent(inout) :: vector_is_zero

  ! vector to hold the new values
  real(kind=myreal), intent(inout), dimension(:) :: new_vector

  ! Jacobi damping factor
  real(kind=myreal), intent(in) :: damping

  integer :: i, j            ! loop indices
  integer :: row             ! current row
  integer :: vec_length      ! length of vector
  integer :: col             ! current column
  real(kind=myreal) :: sm    ! some sum of values
  integer :: s, e            ! start and end f current row
  real(kind=myreal) :: omd   ! 1-damping factor

  omd = one - damping
  vec_length = matrix%n

  do j=1,sweeps

   ! if solution is zero we can bypass some multiplications
   if (vector_is_zero) then
    ! run though rows
    do row = 1, vec_length

    ! get pointer to start of row - start of row holds
    ! the diagonal entry
    s = matrix%ptr(row)

    if (damping .lt. one) then
     new_vector(row) = damping * rhs(row) / matrix%val(s)
    else
     new_vector(row) = rhs(row) / matrix%val(s)
    end if
   end do
   vector_is_zero = .false.

   ! case when vector isn't zero
   else
    ! run through the rows
    do row = 1, vec_length

    ! find start and end of row
    s = matrix%ptr(row)
    e = matrix%ptr(row+1)-1

    sm = rhs(row)

    ! run through the row except for diagonal entry
    do i = s+1, e
     col = matrix%col(i)
     sm = sm - vector(col) * matrix%val(i)
    end do

    if (damping .lt. one) then
     col = matrix%col(s)
     new_vector(row) = damping * sm / matrix%val(s) + omd * vector(col)
    else
     new_vector(row) = sm / matrix%val(s)
    end if

    end do
   end if

   ! set values in vector
   vector = new_vector

   end do

 end subroutine dj_smoother

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    S o l v e    S u b r o u t i n e s        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! T Rees, Dec 2013 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 ! ----------------------
 ! --- Solve no A -------
 ! ----------------------
 
 subroutine mi20_solve_noA(coarse_data,rhs,sol, & 
                           keep,control, solve_control, info)
   type( mi20_data ), intent(in) :: coarse_data(:)
   real(myreal), intent(in) :: rhs(:)
   real(myreal), intent(inout) :: sol(:)
   type( mi20_keep ), intent(inout) :: keep
   type( mi20_control), intent(in) :: control
   type( mi20_solve_control), intent(in) :: solve_control
   type( mi20_info ), intent(out) :: info

   !------------------
   ! run the code...
   if ( keep%zd11_internal_conversion ) then
      ! pull A from keep, and send back 
      call mi20_solve_withA(keep%A,coarse_data,rhs,sol, & 
           keep,control, solve_control, info)
   else
      if ( (control%print_level .gt. 0) .and. (control%error .ge. 0) ) then
         write(control%error,'(a)') 'mi20_solve must be called with ''A'' ', & 
              'if mi20_setup was called'
      end if
      info%flag = -203
      ! flag an error
   end if

 end subroutine mi20_solve_noA


 ! -----------------------
 ! ---- Solve with A -----
 ! -----------------------
 subroutine mi20_solve_withA(matrix,coarse_data,rhs,sol, & 
                             keep,control, solve_control, info)
   type( zd11_type ), intent(in) :: matrix
   type( mi20_data ), intent(in) :: coarse_data(:)
   real(myreal), intent(in) :: rhs(:)
   real(myreal), intent(inout) :: sol(:)
   type( mi20_keep ), intent(inout) :: keep
   type( mi20_control), intent(in) :: control
   type( mi20_solve_control), intent(in) :: solve_control
   type( mi20_info ), intent(out) :: info

   if (.not. solve_control%init_guess) then 
      sol = 0.0
   end if

   keep%max_its = 2*matrix%m
   if (solve_control%max_its > 0) then
      keep%max_its = solve_control%max_its
   end if
! todo :: print out krylov solver being used
   select case (solve_control%krylov_solver)
   case (0) ! simple iteration
      call simple_iteration(matrix,coarse_data,rhs,sol, & 
                            keep,control,solve_control,info)
   case (1) ! CG (MI21)
      call solve_cg(matrix,coarse_data,rhs,sol,keep,control,solve_control,info)
   case (2) ! GMRES (MI24)
      call solve_gmres(matrix,coarse_data,rhs,sol, & 
                       keep,control,solve_control,info)
   case (3) ! BiCGStab (MI26)
      call solve_bicgstab(matrix,coarse_data,rhs,sol, &
                          keep,control,solve_control,info)
   case (4) ! MINRES (HSL_MI32)
      call solve_minres(matrix,coarse_data,rhs,sol, & 
                        keep,control,solve_control,info)
   case default
      write(control%error,'(a)') 'Out of range solve_control%krylov_method'
      info%flag = -120
      return
   end select
 end subroutine mi20_solve_withA
 
 !!!!!!!!!!!!!!!!!!!!!!!!
 !!! Simple iteration !!!
 !!!!!!!!!!!!!!!!!!!!!!!!
 subroutine simple_iteration(matrix,coarse_data,rhs,sol, & 
                             keep,control,solve_control,info)
    type( zd11_type ), intent(in) :: matrix
    type( mi20_data ), intent(in) :: coarse_data(:)
    real(myreal), intent(in) :: rhs(:)
    real(myreal), intent(inout) :: sol(:)
    type( mi20_keep ), intent(inout) :: keep
    type( mi20_control ), intent(in) :: control
    type( mi20_solve_control), intent(in) :: solve_control
    type( mi20_info ), intent(out) :: info
    
    real(myreal) :: res(matrix%m)
    real(myreal) :: norm_res, norm_res_0
    real(myreal) :: dnrm2     ! BLAS function
    
    real(myreal) :: temp(matrix%m)
    integer :: i !,j
    logical :: error_out, print_out
    
    error_out = .false.
    if ( (control%print_level .gt. 0) .and. (control%error .ge. 0) ) &
         error_out = .true.
    
    print_out = .false.
    if ( (control%print_level .gt. 1) .and. (control%print .ge. 0) ) &
         print_out = .true.

    
    call A_prod(matrix%m,matrix%col,matrix%ptr,matrix%val, & 
                sol(1:matrix%m),temp(1:matrix%m))

    res(1:matrix%m) =  rhs(1:matrix%m) - temp(1:matrix%m)
    !    call daxpy(matrix%m,-1.0,temp,1,res,1)
    norm_res_0 = dnrm2(matrix%m,res,1)
    if (norm_res_0 < solve_control%abs_tol) then
       info%residual = norm_res_0
       info%iterations = 0
       if (print_out) then
          write(control%print,'(a)') 'Convergence on entry to solve routine'
          write(control%print,'(a,es12.4)') '2-norm of residual = ', norm_res_0
       end if
       return
    end if
    do i = 1,keep%max_its
       call mi20_precondition(matrix,coarse_data, &
            res,temp,keep,control,info)
       sol(1:matrix%m) = sol(1:matrix%m) + temp(1:matrix%m)
       call A_prod(matrix%m,matrix%col,matrix%ptr,matrix%val,sol,temp)
       res(1:matrix%m) = rhs(1:matrix%m) - temp(1:matrix%m)
!       call daxpy(matrix%m,-1.0,temp,1,res,1)
       norm_res = dnrm2(matrix%m,res,1)
       if (norm_res < & 
            (solve_control%abs_tol + solve_control%rel_tol*norm_res_0)) then
          info%residual = norm_res
          info%iterations = i
          if (print_out) then 
             write(control%print,'(a,i0,a)') 'Convergence in ', i,' iterations'
             write(control%print,'(a,es12.4)') '2-norm of residual = ', norm_res
          end if
          return
       end if
    end do
    info%residual = norm_res
    info%iterations = i
    if (print_out) then
       write(control%print,'(a,es12.4)') '2-norm of residual = ', norm_res
    end if
    if (error_out) then 
       write(control%error,'(a,a,i0)') 'Number of iterations required ', &
        'exceeds the maximum of ', keep%max_its
     end if
    info%flag = -200
    
  end subroutine simple_iteration

  subroutine solve_cg(matrix,coarse_data,rhs,sol, & 
                      keep,control,solve_control,info)

    type( zd11_type ), intent(in) :: matrix
    type( mi20_data ), intent(in) :: coarse_data(:)
    real(myreal), intent(in) :: rhs(:)
    real(myreal), intent(inout) :: sol(:)
    type( mi20_keep ), intent(inout) :: keep
    type( mi20_control ), intent(in) :: control
    type( mi20_solve_control), intent(in) :: solve_control
    type( mi20_info ), intent(out) :: info
    
    ! arrays needed by mi21
    real(myreal) :: cntl(5), rsave(6)
    integer :: icntl(8), isave(10), info21(4)
    real(myreal) :: w(matrix%m,4)
    real(myreal) :: resid
    integer :: locy, locz, iact
    logical :: print_out, error_out

    external mi21id, mi21ad
    
    call mi21id(icntl, cntl, isave, rsave)

    call initialize_controls(w,rhs,sol,cntl,icntl,solve_control,& 
         keep,control%error,matrix%m)
    ! use preconditioning
    icntl(3) = 1
    ! breakdown tolerance
    icntl(7) = 1
    cntl(3) = solve_control%breakdown_tol
    
    error_out = .false.
    if ( (control%print_level .gt. 0) .and. (control%error .ge. 0) ) &
         error_out = .true.

    print_out = .false.
    if ( (control%print_level .gt. 0) .and. (control%print .ge. 0) ) &
         print_out = .true.


    iact = 0
    do
       call mi21ad(iact,matrix%m,w,matrix%m,locy,locz,resid,icntl,cntl,& 
            info21,isave, rsave)
       if (iact == -1) then
          if (error_out) then 
             write(control%error,'(a)') 'Error in solver loop'
          end if
          select case (info21(1))
          case (-3)
             ! the algorithm has broken down
             info%flag = -201
             if (error_out) then
                write(control%error,'(a)') 'Conjugate gradients has broken down'
             end if
          case (-4)
             ! maximum number of iterations reached
             info%flag = -200
          case default
             ! -1 and -2 should never occur because of pre-processing
             if (error_out) then 
                write(control%error,'(a,i0)') 'flag = ', info21(1)
             end if
             info%flag = -222
          end select
          exit
       else if (iact == 1) then 
          info%iterations = info21(2)
          info%residual = resid
          if (print_out) then
             write(control%print,'(a,i0,a)') 'Convergence in ', info21(2),&
                                             ' iterations'
             write(control%print,'(a,es12.4)') '2-norm of residual = ', resid
          end if
          sol(1:matrix%m) = w(1:matrix%m,2)
          exit
       else if (iact == 2) then
          call A_prod(matrix%m,matrix%col,matrix%ptr,matrix%val, & 
                      w(:,locz),w(:,locy))
       else if (iact == 3) then 
          call mi20_precondition(matrix,coarse_data, &
            w(1:matrix%m,locz),w(1:matrix%m,locy),keep,control,info)
          if (info%flag < 0) then 
             if (print_out) then
                write(control%print,'(a)') 'Error return from mi20_precondition'
             end if
          end if
       end if
    end do
    
  end subroutine solve_cg
  
  subroutine solve_gmres(matrix,coarse_data,rhs,sol, & 
                         keep,control,solve_control,info)
    type( zd11_type ), intent(in) :: matrix
    type( mi20_data ), intent(in) :: coarse_data(:)
    real(myreal), intent(in) :: rhs(:)
    real(myreal), intent(inout) :: sol(:)
    type( mi20_keep ), intent(inout) :: keep
    type( mi20_control ), intent(in) :: control
    type( mi20_solve_control), intent(in) :: solve_control
    type( mi20_info ), intent(out) :: info
    
!    integer :: preconditioner_side

    ! arrays needed by mi24
    real(myreal) ::  cntl(5), rsave(9)
    integer :: icntl(8), isave(17), info24(4)
    real(myreal) :: w(matrix%m,solve_control%gmres_restart + 7)
    real(myreal) :: h(solve_control%gmres_restart + 1, &
                      solve_control%gmres_restart + 2)
    real(myreal) :: resid
    integer :: locy, locz, iact, ldh
!    integer :: restart
    logical :: lsave(4)
    logical :: print_out, error_out
    
    integer :: precond, no_precond!, i
    
    external mi24id, mi24ad
    
    call mi24id(icntl, cntl, isave, rsave, lsave)

    call initialize_controls(w,rhs,sol,cntl,icntl,solve_control,& 
         keep,control%error,matrix%m)
    ! use preconditioning
    icntl(3) = 3

    ldh = solve_control%gmres_restart + 1
   
    error_out = .false.
    if ( (control%print_level .gt. 0) .and. (control%error .ge. 0) ) &
         error_out = .true.
    
    print_out = .false.
    if ( (control%print_level .gt. 1) .and. (control%print .ge. 0) ) &
         print_out = .true. 
    
    if (solve_control%preconditioner_side < 0) then
       ! left precondition
       precond = 3
       no_precond = 4
    else
       ! right precondition
       precond = 4
       no_precond = 3
    end if
    
    iact = 0
    do
       call mi24ad(iact,matrix%m,solve_control%gmres_restart, & 
                   w,matrix%m,locy,locz,h,ldh,resid,&
            icntl,cntl,info24,isave, rsave, lsave)
       if (iact == -1) then
          if (error_out) then 
             write(control%error,'(a)') 'Error in solver loop'
          end if
          select case (info24(1))
          case (-5)
             ! maximum number of iterations reached
             info%flag = -200
          case default
             ! -1 to -4 will never occur because of pre-processing
             if (error_out) then 
                write(control%error,'(a,i0)') 'flag = ', info24(1)
             end if
             info%flag = -222
          end select
          exit
       else if (iact == 1) then 
          info%iterations = info24(2)
          info%residual = resid
          if (print_out) then
             write(control%print,'(a,i0,a)') 'Convergence in ', info24(2), & 
                                             ' iterations'
             write(control%print,'(a,es12.4)') '2-norm of residual = ', resid
          end if
          sol(1:matrix%m) = w(1:matrix%m,2)
          exit
       else if (iact == 2) then
          call A_prod(matrix%m,matrix%col,matrix%ptr,matrix%val, & 
                      w(1:matrix%m,locz),w(1:matrix%m,locy))
       else if (iact == precond) then 
          call mi20_precondition(matrix,coarse_data, &
            w(1:matrix%m,locz),w(1:matrix%m,locy),keep,control,info)
          if (info%flag < 0) then 
             if (error_out) then 
                write(control%error,'(a)') 'Error return from mi20_precondition'
             end if
          end if
       else if (iact == no_precond) then
          w(1:matrix%m,locy) = w(1:matrix%m,locz)
       end if
    end do

  end subroutine solve_gmres
  
  subroutine solve_bicgstab(matrix,coarse_data,rhs,sol, & 
                            keep,control,solve_control,info)
    type( zd11_type ), intent(in) :: matrix
    type( mi20_data ), intent(in) :: coarse_data(:)
    real(myreal), intent(in) :: rhs(:)
    real(myreal), intent(inout) :: sol(:)
    type( mi20_keep ), intent(inout) :: keep
    type( mi20_control ), intent(in) :: control
    type( mi20_solve_control), intent(in) :: solve_control
    type( mi20_info ), intent(out) :: info

    ! arrays needed by mi26
    real(myreal) ::  cntl(5), rsave(9)
    integer :: icntl(8), isave(14), info26(4)
    real(myreal) :: w(matrix%m,8)
    real(myreal) :: resid
    integer :: locy, locz, iact
    logical :: print_out, error_out

!    integer :: restart
    
    external mi26id, mi26ad
    
    call mi26id(icntl, cntl, isave, rsave)

    call initialize_controls(w,rhs,sol,cntl,icntl,solve_control, &
         keep,control%error,matrix%m)
    ! use preconditioning
    icntl(3) = 1
    ! breakdown tolerance
    cntl(3) = solve_control%breakdown_tol
    
    error_out = .false.
    if ( (control%print_level .gt. 0) .and. (control%error .ge. 0) ) &
         error_out = .true.

    print_out = .false.
    if ( (control%print_level .gt. 1) .and. (control%print .ge. 0) ) &
         print_out = .true.

    iact = 0
    do
       call mi26ad(iact,matrix%m,w,matrix%m,locy,locz,resid,&
            icntl,cntl,info26,isave, rsave)
       if (iact == -1) then
          if (error_out) then 
             write(control%error,'(a)') 'Error in solver loop'
          end if
          select case (info26(1))
          case (-3)
             ! the algorithm has broken down
             info%flag = -201
          case (-4)
             ! maximum number of iterations reached
             info%flag = -200
          case default
             ! -1 and -2 will never occur because of pre-processing
             if (error_out) then 
                write(control%error,'(a,i0)') 'flag = ', info26(1)
             end if
             info%flag = -222
          end select
          exit
       else if (iact == 1) then 
          info%iterations = info26(2)
          info%residual = resid

          if (print_out) then  
             write(control%print,'(a,i0,a)') 'Convergence in ', info26(2),& 
                                             ' iterations'
             write(control%print,'(a,es12.4)') '2-norm of residual = ', resid
          end if
          sol(1:matrix%m) = w(1:matrix%m,2)
          exit
       else if (iact == 2) then
          call A_prod(matrix%m,matrix%col,matrix%ptr,matrix%val, & 
                      w(1:matrix%m,locz),w(1:matrix%m,locy))
       else if (iact == 3) then 
          call mi20_precondition(matrix,coarse_data, &
            w(1:matrix%m,locz),w(1:matrix%m,locy),keep,control,info)
          if (info%flag < 0) then 
             if (error_out) then
                write(control%error,'(a)') 'Error return from mi20_precondition'
             end if
          end if
       end if
    end do
    
  end subroutine solve_bicgstab

  subroutine solve_minres(matrix,coarse_data,rhs,sol,& 
                          keep,control,solve_control,info)
    
    use hsl_mi32_double
    
    type( zd11_type ), intent(in) :: matrix
    type( mi20_data ), intent(in) :: coarse_data(:)
    real(myreal), intent(in) :: rhs(:)
    real(myreal), intent(inout) :: sol(:)
    type( mi20_keep ), intent(inout) :: keep
    type( mi20_control ), intent(in) :: control
    type( mi20_solve_control), intent(in) :: solve_control
    type( mi20_info ), intent(out) :: info

    ! arrays needed by hsl_mi32
    real(myreal) :: V_in( matrix%n )
    real(myreal), pointer :: V_out( : )
    type ( mi32_keep ) :: mi32keep
    type ( mi32_control ) :: mi32control
    type ( mi32_info ) :: mi32info
    logical :: print_out, error_out
    integer :: action 

!    integer :: restart
    
!    CALL MI32_INITIALIZE( mi32keep, mi32control )

    error_out = .false.
    if ( (control%print_level .gt. 0) .and. (control%error .ge. 0) ) &
         error_out = .true.

    print_out = .false.
    if ( (control%print_level .gt. 1) .and. (control%print .ge. 0) ) &
         print_out = .true.

    ! calculate the residual
    ! V_in <- A*x
    call A_prod(matrix%m,matrix%col,matrix%ptr,matrix%val, &
                sol(1:matrix%m),V_in(1:matrix%m))
    ! V_in <- V_in  - b
    V_in = V_in - rhs(1:matrix%m)

    action = 1

    mi32control%stop_relative = solve_control%rel_tol
!    mi32control%stop_absolute = epsilon(1.0d0) ! put this in to guard against 
    mi32control%itmax = solve_control%max_its

    do
       call mi32_minres(action,matrix%m, sol(1:matrix%m), V_in, V_out, & 
            mi32keep, mi32control, mi32info)
       select case (action)
          case (-3)
             ! the algorithm has broken down
             info%flag = -201
             exit
          case (-2)
             ! maximum number of iterations reached
             info%flag = -200
             exit
          case ( -1, :-4)
             ! other errors should never occur because of preprocessing
             if (error_out) then 
                write(control%error,'(a,i0)') 'flag = ', action
             end if
             info%flag = -222
             exit
          case ( 2 ) 
             ! perform preconditioning
             call mi20_precondition(matrix,coarse_data, &
                                    V_out,V_in,keep,control,info)
          case ( 3 ) 
             ! perform matrix-vector product
             call A_prod(matrix%m,matrix%col,matrix%ptr,matrix%val, &
                         V_out,V_in)
          case ( 0 ) 
             ! converged!
             info%iterations = mi32info%iter
             info%residual = mi32info%rnorm
             if (print_out) then 
                write(control%print,'(a,i0,a)') 'Convergence in ', &
                     mi32info%iter,' iterations'
                write(control%print,'(a,es12.4)') & 
                     'Preconditioned norm of residual = ', mi32info%rnorm
             end if
             exit
          end select

    end do
    
  end subroutine solve_minres


  subroutine A_prod(m,col,ptr,val,x,y)
    integer, intent(in) :: m
    integer, intent(in) :: col(:)
    real(myreal), intent(in) :: val(:)
    integer, intent(in) :: ptr(:)
    real(myreal), intent(in) :: x(:)
    real(myreal), intent(out) :: y(:)
    
    integer :: i, j

    y = 0.0
    do i = 1, m
       do j = ptr(i),ptr(i+1)-1
          y(i) = y(i) + val(j)*x(col(j))
       end do
    end do
    
  end subroutine A_prod

  subroutine initialize_controls(w,rhs,sol,cntl,icntl, & 
                                 solve_control,keep,error,m)
    real(myreal), intent(out) :: w(:,:)
    real(myreal), intent(in) :: rhs(:)
    real(myreal), intent(in) :: sol(:)
    real(myreal), intent(out) :: cntl(:)
    integer, intent(out) :: icntl(:)
    type( mi20_solve_control ), intent(in) :: solve_control
    type( mi20_keep ), intent(in) :: keep
    integer :: error, m
    
    w(1:m,1) = rhs(1:m)
    ! set the stream number of error messages...
    icntl(1) = error
    ! set the stream number for warnings...
    icntl(2) = error
    ! use the default convergence test
    icntl(4) = 0
    cntl(1) = solve_control%rel_tol
    cntl(2) = solve_control%abs_tol
    ! set the starting vector...
    icntl(5) = 1
    w(1:m,2) = sol(1:m)
    ! set the max its
    icntl(6) = keep%max_its
    
  end subroutine initialize_controls

end module hsl_mi20_double
