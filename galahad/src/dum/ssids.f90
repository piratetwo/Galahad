! THIS VERSION: 05/11/2013 AT 13:35:00 GMT.

!-*-*-*-*-*-  G A L A H A D  -  D U M M Y   S S I D S   M O D U L E  -*-*-*-*-*-

MODULE SPRAL_SSIDS
    
   USE GALAHAD_SYMBOLS
!  use iso_c_binding

! Parameters (all private)
  integer, parameter, private  :: short = kind(0)
  integer(short), parameter, private  :: wp = kind(0.0d0)
  integer(short), parameter, private  :: long = selected_int_kind(18)
  real(wp), parameter, private :: one = 1.0_wp
  real(wp), parameter, private :: zero = 0.0_wp
  integer(short), parameter, private :: nemin_default = 8

  interface SSIDS_analyse
      module procedure analyse_double
  end interface

  interface SSIDS_analyse_coord
      module procedure SSIDS_analyse_coord_double
  end interface

  interface SSIDS_factor
      module procedure SSIDS_factor_double
  end interface

!  interface SSIDS_factor_solve
!     module procedure SSIDS_factor_solve_one_double
!     module procedure SSIDS_factor_solve_mult_double
! end interface

  interface SSIDS_solve
      module procedure SSIDS_solve_one_double
      module procedure SSIDS_solve_mult_double
  end interface

  interface SSIDS_solve_fredholm
      module procedure SSIDS_solve_fredholm_double
   end interface SSIDS_solve_fredholm

  interface SSIDS_enquire_posdef
    module procedure SSIDS_enquire_posdef_double
  end interface

  interface SSIDS_enquire_indef
    module procedure SSIDS_enquire_indef_double
  end interface

  interface SSIDS_alter
    module procedure SSIDS_alter_double
  end interface

  interface SSIDS_sparse_fwd_solve
      module procedure SSIDS_sparse_fwd_solve_double
   end interface SSIDS_sparse_fwd_solve

   interface SSIDS_free
      module procedure free_akeep_double
      module procedure free_fkeep_double
   end interface SSIDS_free

  interface SSIDS_finalise
      module procedure finalise_both_double
  end interface

   type auction_options
      integer :: max_iterations = 30000
      integer :: max_unchanged(3) = (/ 10,   100, 100 /)
      real :: min_proportion(3) = (/ 0.90, 0.0, 0.0 /)
      real :: eps = 0.01
   end type auction_options

   type auction_inform
      integer :: matched = 0 ! #matched rows/cols
      integer :: iterations = 0 ! #iterations
   end type auction_inform

  type SSIDS_options ! The scalar control of this type controls the action
    logical :: action = .true. ! pos_def = .false. only.
    integer(short) :: nemin = nemin_default
    real(wp) :: multiplier = 1.1 
    integer :: ordering = 1 ! controls choice of ordering
    integer(short) :: print_level = 0 ! Controls diagnostic printing
    integer :: scaling = 0 ! controls use of scaling. 
    real(wp) :: small = tiny(one) ! Minimum pivot size
    real(wp) :: u = 0.01 ! Initial relative pivot threshold
    integer(short) :: unit_diagnostics = 6 ! unit for diagnostic printing.
    integer(short) :: unit_error = 6 ! unit number for error messages
    integer(short) :: unit_warning = 6 ! unit number for warning messages
    type(auction_options) :: auction
    integer(long) :: min_subtree_work = 1e5 ! Minimum amount of work
    integer :: min_ldsrk_work = 1e4 ! Minimum amount of work to aim
    logical :: use_gpu_solve = .true.
    integer :: nstream = 1
    integer :: presolve = 0
!   type(C_PTR) :: cublas
  end type SSIDS_options

  type SSIDS_inform ! The scalar info of this type returns information to user.
    integer(short) :: flag = 0     ! Value zero after successful entry.
    integer(short) :: flag77 = 0 ! error flag from mc77
    integer(short) :: matrix_dup = 0 ! Number of duplicated entries.
    integer(short) :: matrix_rank = 0 ! Rank of factorized matrix.
    integer(short) :: matrix_outrange = 0 ! Number of out-of-range entries.
    integer(short) :: matrix_missing_diag = 0 ! Number of missing diag. entries
    integer(short) :: maxdepth = 0 ! Maximum depth of the tree.
    integer(short) :: maxfront = 0 ! Maximum front size.
    integer(long)  :: num_factor = 0_long ! Number of entries in the factor.
    integer(long)  :: num_flops = 0_long ! Number of flops to calculate L.
    integer(short) :: num_delay = 0 ! Number of delayed eliminations.
    integer(short) :: num_neg = 0 ! Number of negative eigenvalues.
    integer(short) :: num_sup = 0 ! Number of supervariables. 
    integer(short) :: num_two = 0 ! Number of 2x2 pivots.
    integer :: ordering = 0 ! ordering actually used
    integer(short) :: stat = 0 ! STAT value (when available).
    integer :: cuda_error = 0 ! cuda error value
    type(auction_inform) :: auction
   end type SSIDS_inform

   type ssids_akeep
      logical :: check ! copy of check as input to analyse phase
      integer :: flag ! copy of error flag.
      integer :: maxmn ! maximum value of blkm or blkn
      integer :: n ! Dimension of matrix
      integer :: ne ! Set to number of entries input by user.
      integer(long) :: nfactor 
      integer :: nnodes = -1 ! Number of nodes in assembly tree
      integer :: num_two ! # 2x2 pivots, set to 0 as we ignore any -ve signs
      integer, dimension(:), allocatable :: child_ptr 
      integer, dimension(:), allocatable :: child_list
      integer, dimension(:), allocatable :: invp ! inverse of pivot order
      integer, dimension(:), allocatable :: level ! level(i) of assembly tree
      integer, dimension(:,:), allocatable :: nlist ! map from A to factors
      integer, dimension(:), allocatable :: nptr ! Entries into nlist for nodes
      integer, dimension(:), allocatable :: rlist ! rlist(rptr(i):rptr(i+1)-1)
      integer(long), dimension(:), allocatable :: rptr ! Pointers into rlist
      integer, dimension(:), allocatable :: sparent ! sparent(i) is parent
      integer, dimension(:), allocatable :: sptr ! (super)node pointers.
      integer(long), dimension(:), allocatable :: subtree_work ! For each node,
      integer, dimension(:), allocatable :: rlist_direct ! List of rows
      integer, allocatable :: ptr(:) ! column pointers
      integer, allocatable :: row(:)! row indices
      integer :: lmap ! used by hsl_mc69
      integer, allocatable :: map(:) ! used by hsl_mc69
      integer :: matrix_dup
      integer :: matrix_outrange
      integer :: matrix_missing_diag
      integer :: maxdepth
      integer(long) :: num_flops ! determine if parallelism should be used
      integer :: num_sup
      integer :: ordering
      real(wp), dimension(:), allocatable :: scaling
!     type(C_PTR) :: gpu_nlist = C_NULL_PTR
!     type(C_PTR) :: gpu_rlist = C_NULL_PTR
!     type(C_PTR) :: gpu_rlist_direct = C_NULL_PTR
   end type ssids_akeep
   
   type smalloc_type
      real(wp), dimension(:), allocatable :: rmem ! real memory
      integer(long) :: rmem_size ! needed as size(rmem,kind=long) is f2003
      integer(long) :: rhead = 0 ! last location containing useful info in mem
      integer, dimension(:), allocatable :: imem ! integer memory
      integer(long) :: imem_size ! needed as size(imem,kind=long) is f2003
      integer(long) :: ihead = 0 ! last location containing useful infor in mem
      type(smalloc_type), pointer :: next_alloc => null()
      type(smalloc_type), pointer :: top_real => null() ! Last page r alloc ok
      type(smalloc_type), pointer :: top_int => null() ! Last page i alloc ok
!$    integer(omp_lock_kind) :: lock
   end type smalloc_type

   type node_type ! Data type for storing each node of the factors
      integer :: nelim
      integer :: ndelay
      integer(long) :: rdptr ! entry into (rebuilt) rlist_direct
      integer :: ncpdb ! #contrib. to parent's diag. block
!     type(C_PTR) :: gpu_lcol
      real(wp), dimension(:), pointer :: lcol ! values in factors
      integer, dimension(:), pointer :: perm ! permutation of columns at node
      type(smalloc_type), pointer :: rsmptr, ismptr
      integer(long) :: rsmsa, ismsa
   end type node_type

   type eltree_level
!     type(C_PTR) :: ptr_levL ! device pointer
      integer :: lx_size
      integer :: lc_size
      integer :: max_nch
      integer :: total_nch
      integer :: off_col_ind
      integer :: off_row_ind
      integer :: ncp_pre
      integer :: ncb_asm_pre
      integer :: ncp_post
      integer :: ncb_asm_post
      integer :: ncb_slv_n
      integer :: ncb_slv_t
!     type(C_PTR) :: gpu_cpdata_pre
!     type(C_PTR) :: gpu_blkdata_pre
!     type(C_PTR) :: gpu_cpdata_post
!     type(C_PTR) :: gpu_blkdata_post
!     type(C_PTR) :: gpu_solve_n_data
!     type(C_PTR) :: gpu_solve_t_data
   end type eltree_level

   type gpu_type
      integer :: n = 0
      integer :: nnodes = 0
      integer :: num_levels ! number of levels
      integer, dimension(:), allocatable :: lvlptr ! pointers into lvllist
      integer, dimension(:), allocatable :: lvllist ! list of nodes at level
      integer, dimension(:), allocatable :: off_L ! offsets for each node
      integer, dimension(:), allocatable :: off_lx ! node offsets for fwd solve
      integer, dimension(:), allocatable :: off_lc ! offsets for node contrib.
      integer, dimension(:), allocatable :: off_ln ! node offsets for bwd solve
      integer, dimension(:), allocatable :: rlist_direct
      integer(long) :: rd_size = 0
      integer :: max_lx_size = 0
      integer :: max_lc_size = 0
      type(eltree_level), dimension(:), allocatable :: values_L(:) ! data
!     type(C_PTR) :: gpu_rlist_direct
!     type(C_PTR) :: gpu_col_ind
!     type(C_PTR) :: gpu_row_ind
!     type(C_PTR) :: gpu_diag
!     type(C_PTR) :: gpu_sync
   end type gpu_type

!  type gpu_slv_type
!     integer(C_SIZE_T), dimension(:), allocatable :: rptr
!     type(C_PTR) :: rlist
!  end type gpu_slv_type

   type ssids_fkeep
      integer :: flag ! copy of error flag.
      real(wp), dimension(:), allocatable :: scaling ! scaling for enties
      type(node_type), dimension(:), allocatable :: nodes ! Stores node pointers
      type(smalloc_type), pointer :: alloc=>null() ! Linked list of memory pages
      logical :: pos_def ! set to true if user indicates matrix pos. definite
      integer :: matrix_rank
      integer :: maxfront
      integer :: num_delay
      integer(long) :: num_factor
      integer(long) :: num_flops
      integer :: num_neg
      integer :: num_two
!     type(C_PTR), dimension(:), allocatable :: stream_handle
      type(gpu_type), dimension(:), allocatable :: stream_data
      type(gpu_type) :: top_data
!     type(gpu_slv_type) :: gpu_slv
      logical :: host_factors = .false.
   end type ssids_fkeep

contains

  subroutine analyse_double(check, n, ptr, row, akeep,                         &
                            options, inform, order, val)
   logical, intent(in) :: check
   integer, intent(in) :: n
   integer, intent(in) :: row(:), ptr(:)
   type (SSIDS_akeep), intent (out) :: akeep
   type (SSIDS_options), intent (in) :: options
   type (SSIDS_inform), intent (out) :: inform
   integer(short), OPTIONAL, intent (inout) :: order(:)
   real(wp), OPTIONAL, intent(in) :: val(:) 

   IF ( options%unit_error >= 0 .AND. options%print_level > 0 )                &
     WRITE( options%unit_error,                                                &
         "( ' We regret that the solution options that you have ', /,          &
  &         ' chosen are not all freely available with GALAHAD.', /,           &
  &         ' If you have SPRAL, this option may be enabled by', /,            &
  &         ' replacing the dummy subroutine SSIDS_analyse', /,                &
  &         ' with its SPRAL namesake and dependencies. See ', /,              &
  &         '   $GALAHAD/src/makedefs/packages for details.' )" )
   inform%flag = GALAHAD_error_unknown_solver

  end subroutine analyse_double

  subroutine SSIDS_analyse_coord_double( n, ne, row, col, akeep,               &
                                         options, inform, order, val)
   integer, intent(in) :: n, ne
   integer, intent(in) :: row(:), col(:)
   type (SSIDS_akeep), intent (out) :: akeep
   type (SSIDS_options), intent(in) :: options
   type (SSIDS_inform), intent(out) :: inform
   integer(short), OPTIONAL, intent (inout) :: order(:)
   real(wp), OPTIONAL, intent(in) :: val(:) 

   IF ( options%unit_error >= 0 .AND. options%print_level > 0 )                &
     WRITE( options%unit_error,                                                &
         "( ' We regret that the solution options that you have ', /,          &
  &         ' chosen are not all freely available with GALAHAD.', /,           &
  &         ' If you have SPRAL, this option may be enabled by', /,            &
  &         ' replacing the dummy subroutine SSIDS_analyse_coord', /,          &
  &         ' with its SPRAL namesake and dependencies. See ', /,              &
  &         '   $GALAHAD/src/makedefs/packages for details.' )" )
   inform%flag = GALAHAD_error_unknown_solver

  end subroutine SSIDS_analyse_coord_double

  subroutine SSIDS_factor_double(posdef,val,akeep,fkeep,options,inform,        &
                                 scale,ptr,row)
   logical, intent(in) :: posdef
   real(wp), intent(in), target :: val(*)
   type (SSIDS_akeep), intent (in) :: akeep
   type (SSIDS_fkeep), intent (out) :: fkeep
   type (SSIDS_options), intent (in) :: options
   type (SSIDS_inform), intent (inout) :: inform
   real(wp), intent(inout), optional :: scale(:)
   integer(short), intent(in), optional :: ptr(akeep%n+1)
   integer(short), intent(in), optional :: row(*)

   IF ( options%unit_error >= 0 .AND. options%print_level > 0 )                &
     WRITE( options%unit_error,                                                &
         "( ' We regret that the solution options that you have ', /,          &
  &         ' chosen are not all freely available with GALAHAD.', /,           &
  &         ' If you have SPRAL, this option may be enabled by', /,            &
  &         ' replacing the dummy subroutine SSIDS_factor', /,                 &
  &         ' with its SPRAL namesake and dependencies. See ', /,              &
  &         '   $GALAHAD/src/makedefs/packages for details.' )" )
   inform%flag = GALAHAD_error_unknown_solver

  end subroutine SSIDS_factor_double

! subroutine SSIDS_factor_solve_double(matrix_type,val,nrhs,x,lx,akeep,fkeep,  &
!                                      options,inform,scale,ptr,row)
!  real(wp), intent(in) :: val(*)
!  integer(short) :: lx, nrhs
!  real(wp), intent(inout) :: x(lx,nrhs)
!  type (SSIDS_akeep), intent (in) :: akeep
!  type (SSIDS_fkeep), intent (out) :: fkeep
!  type (SSIDS_options), intent (in) :: options
!  type (SSIDS_inform), intent (inout) :: inform
!  real(wp), intent(inout), optional :: scale(:) 
!  integer(short), intent(in), optional :: ptr(akeep%n+1)
!  integer(short), intent(in), optional :: row(*)

!  IF ( options%unit_error >= 0 .AND. options%print_level > 0 )                &
!    WRITE( options%unit_error,                                                &
!        "( ' We regret that the solution options that you have ', /,          &
! &         ' chosen are not all freely available with GALAHAD.', /,           &
! &         ' If you have SPRAL, this option may be enabled by', /,            &
! &         ' replacing the dummy subroutine SSIDS_factor_solve', /,           &
! &         ' with its SPRAL namesake and dependencies. See ', /,              &
! &         '   $GALAHAD/src/makedefs/packages for details.' )" )
!  inform%flag = GALAHAD_error_unknown_solver

! end subroutine SSIDS_factor_solve_double

! subroutine SSIDS_factor_solve_one_double(matrix_type,val,x1,akeep,fkeep,     &
!                                          options,inform,scale,ptr,row)
!  integer(short), intent(in) :: matrix_type 
!  real(wp), intent(in) :: val(*)
!  real(wp), intent(inout) :: x1(:) 
!  type (SSIDS_akeep), intent (in) :: akeep
!  type (SSIDS_fkeep), intent(out) :: fkeep
!  type (SSIDS_options), intent (in) :: options
!  type (SSIDS_inform), intent (inout) :: inform
!  real(wp), intent(inout), optional :: scale(:)
!  integer(short), intent(in), optional :: ptr(akeep%n+1)
!  integer(short), intent(in), optional :: row(*)

!  IF ( options%unit_error >= 0 .AND. options%print_level > 0 )                &
!    WRITE( options%unit_error,                                                &
!        "( ' We regret that the solution options that you have ', /,          &
! &         ' chosen are not all freely available with GALAHAD.', /,           &
! &         ' If you have SPRAL, this option may be enabled by', /,            &
! &         ' replacing the dummy subroutine SSIDS_factor_solve_one', /,       &
! &         ' and dependencies. See ', /,                                      &
! &         '   $GALAHAD/src/makedefs/packages for details.' )" )
!  inform%flag = GALAHAD_error_unknown_solver

! end subroutine SSIDS_factor_solve_one_double

  subroutine SSIDS_solve_mult_double(nrhs,x,ldx,akeep,fkeep,options,inform,job)
   integer(short), intent (in) :: nrhs, ldx
   real(wp), intent (inout) :: x(ldx,nrhs)
   type (SSIDS_akeep), intent (in) :: akeep
   type (SSIDS_fkeep), intent (in) :: fkeep
   type (SSIDS_options), intent (in) :: options
   type (SSIDS_inform), intent (inout) :: inform
   integer(short), optional, intent (in) :: job

   IF ( options%unit_error >= 0 .AND. options%print_level > 0 )                &
     WRITE( options%unit_error,                                                &
         "( ' We regret that the solution options that you have ', /,          &
  &         ' chosen are not all freely available with GALAHAD.', /,           &
  &         ' If you have SPRAL, this option may be enabled by', /,            &
  &         ' replacing the dummy subroutine SSIDS_solve', /,                  &
  &         ' with its SPRAL namesake and dependencies. See ', /,              &
  &         '   $GALAHAD/src/makedefs/packages for details.' )" )
   inform%flag = GALAHAD_error_unknown_solver

  end subroutine SSIDS_solve_mult_double

  subroutine SSIDS_solve_one_double(x,akeep,fkeep,options,inform,job)
   real(wp), intent (inout) :: x(:)
   type (SSIDS_akeep), intent (in) :: akeep
   type (SSIDS_fkeep), intent (in) :: fkeep
   type (SSIDS_options), intent (in) :: options
   type (SSIDS_inform), intent (inout) :: inform
   integer(short), optional, intent (in) :: job

   IF ( options%unit_error >= 0 .AND. options%print_level > 0 )                &
     WRITE( options%unit_error,                                                &
         "( ' We regret that the solution options that you have ', /,          &
  &         ' chosen are not all freely available with GALAHAD.', /,           &
  &         ' If you have SPRAL, this option may be enabled by', /,            &
  &         ' replacing the dummy subroutine SSIDS_solve_one ', /,             &
  &         ' with its SPRAL namesake and dependencies. See ', /,              &
  &         '   $GALAHAD/src/makedefs/packages for details.' )" )
   inform%flag = GALAHAD_error_unknown_solver

  end subroutine SSIDS_solve_one_double

  subroutine SSIDS_solve_fredholm_double( nrhs, flag_out, x, ldx,              &
                                          akeep, fkeep, options, inform )
   integer, intent(in) :: nrhs
   logical, intent(out) :: flag_out(nrhs)
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,2*nrhs), intent(inout) :: x
   type(SSIDS_akeep), intent(in) :: akeep
   type(SSIDS_fkeep), intent(in) :: fkeep
   type(SSIDS_options), intent(in) :: options
   type(SSIDS_inform), intent(out) :: inform

   IF ( options%unit_error >= 0 .AND. options%print_level > 0 )                &
     WRITE( options%unit_error,                                                &
         "( ' We regret that the solution options that you have ', /,          &
  &         ' chosen are not all freely available with GALAHAD.', /,           &
  &         ' If you have SPRAL, this option may be enabled by', /,            &
  &         ' replacing the dummy subroutine SSIDS_solve_fredholm', /,         &
  &         ' with its SPRAL namesake and dependencies. See ', /,              &
  &         '   $GALAHAD/src/makedefs/packages for details.' )" )
   inform%flag = GALAHAD_error_unknown_solver
  end subroutine SSIDS_solve_fredholm_double

  subroutine SSIDS_enquire_posdef_double(akeep,fkeep,options,inform,d)
    real(wp), dimension( : ), intent(out) :: d
    type (SSIDS_akeep), intent (in) :: akeep
    type (SSIDS_fkeep), intent(in) :: fkeep
    type (SSIDS_options), intent (inout) :: options
    type (SSIDS_inform), intent (inout) :: inform

   IF ( options%unit_error >= 0 .AND. options%print_level > 0 )                &
     WRITE( options%unit_error,                                                &
         "( ' We regret that the solution options that you have ', /,          &
  &         ' chosen are not all freely available with GALAHAD.', /,           &
  &         ' If you have SPRAL, this option may be enabled by', /,            &
  &         ' replacing the dummy subroutine SSIDS_enquire_posdef', /,         &
  &         ' with its SPRAL namesake and dependencies. See ', /,              &
  &         '   $GALAHAD/src/makedefs/packages for details.' )" )
   inform%flag = GALAHAD_error_unknown_solver

  end subroutine SSIDS_enquire_posdef_double

  subroutine SSIDS_enquire_indef_double(akeep,fkeep,options,inform,piv_order,d)
    integer(short), optional, intent(out) :: piv_order(:)
    real(wp), optional, intent(out) :: d(:,:)
    type (SSIDS_akeep), intent (in) :: akeep
    type (SSIDS_fkeep), intent (in) :: fkeep
    type (SSIDS_options), intent (inout) :: options
    type (SSIDS_inform), intent (inout) :: inform

   IF ( options%unit_error >= 0 .AND. options%print_level > 0 )                &
     WRITE( options%unit_error,                                                &
         "( ' We regret that the solution options that you have ', /,          &
  &         ' chosen are not all freely available with GALAHAD.', /,           &
  &         ' If you have SPRAL, this option may be enabled by', /,            &
  &         ' replacing the dummy subroutine SSIDS_enquire_indef', /,          &
  &         ' with its SPRAL namesake and dependencies. See ', /,              &
  &         '   $GALAHAD/src/makedefs/packages for details.' )" )
   inform%flag = GALAHAD_error_unknown_solver

  end subroutine SSIDS_enquire_indef_double

  subroutine SSIDS_alter_double(d,akeep,fkeep,options,inform)
    real(wp), intent (in) :: d(:,:)
    type (SSIDS_akeep), intent (in) :: akeep
    type (SSIDS_fkeep), intent (in) :: fkeep
    type (SSIDS_options), intent (inout) :: options
    type (SSIDS_inform), intent (inout) :: inform

   IF ( options%unit_error >= 0 .AND. options%print_level > 0 )                &
     WRITE( options%unit_error,                                                &
         "( ' We regret that the solution options that you have ', /,          &
  &         ' chosen are not all freely available with GALAHAD.', /,           &
  &         ' If you have SPRAL, this option may be enabled by', /,            &
  &         ' replacing the dummy subroutine SSIDS_alter', /,                  &
  &         ' with its SPRAL namesake and dependencies. See ', /,              &
  &         '   $GALAHAD/src/makedefs/packages for details.' )" )
   inform%flag = GALAHAD_error_unknown_solver

  end subroutine SSIDS_alter_double

  subroutine SSIDS_sparse_fwd_solve_double(nbi, bindex, b, order, lflag,       &
      nxi, xindex, x, akeep, fkeep, options, inform)
   integer, intent(in) :: nbi
   integer, intent(in) :: bindex(:)
   real(wp), intent(in) :: b(:)
   integer, intent(in) :: order(:)
   logical, intent(inout), dimension(:) :: lflag
   integer, intent(out) :: nxi
   integer, intent(out) :: xindex(:)
   real(wp), intent(inout) :: x(:)
   type(SSIDS_akeep), intent(in) :: akeep
   type(SSIDS_fkeep), intent(in) :: fkeep
   type(SSIDS_options), intent(in) :: options
   type(SSIDS_inform), intent(out) :: inform
   IF ( options%unit_error >= 0 .AND. options%print_level > 0 )                &
     WRITE( options%unit_error,                                                &
         "( ' We regret that the solution options that you have ', /,          &
  &         ' chosen are not all freely available with GALAHAD.', /,           &
  &         ' If you have SPRAL, this option may be enabled by', /,            &
  &         ' replacing the dummy subroutine SSIDS_sparse_fwd_solve ', /,      &
  &         ' with its SPRAL namesake and dependencies. See ', /,              &
  &         '   $GALAHAD/src/makedefs/packages for details.' )" )
   inform%flag = GALAHAD_error_unknown_solver

  end subroutine SSIDS_sparse_fwd_solve_double

  subroutine free_akeep_double(akeep)
     type(ssids_akeep), intent(inout) :: akeep
  end subroutine free_akeep_double

  subroutine free_fkeep_double(fkeep)
   type(ssids_fkeep), intent(inout) :: fkeep
  end subroutine free_fkeep_double

  subroutine finalise_both_double(akeep,fkeep)
    type (SSIDS_akeep), intent (inout) :: akeep
    type (SSIDS_fkeep), intent (inout) :: fkeep
  end subroutine finalise_both_double

END MODULE SPRAL_SSIDS
