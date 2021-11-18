! COPYRIGHT (c) 2006 Council for the Central Laboratory
!               of the Research Councils
! Original date 21 August 2006. Version 1.0.0.
!
! Version 1.4.1
! See ChangeLog for history

module hsl_ma54_double

  implicit none
  EXTERNAL dgemm, dtrsm, dsyrk, dcopy, dtpsv, dgemv

  integer,parameter,private :: wp = kind(1.0d0) ! Precision parameter.
  integer,parameter,private :: long = selected_int_kind(18) ! Long integer.
  real(wp),parameter,private :: unity = 1.0_wp
  integer(long),parameter,private :: one = 1_long

  interface ma54_to_block
      module procedure ma54_to_block_double
  end interface

  interface ma54_from_block
      module procedure ma54_from_block_double
  end interface

  interface ma54_factor
      module procedure ma54_factor_double
  end interface

  interface ma54_solve1
      module procedure ma54_solve1_double
  end interface

  interface ma54_solve2
      module procedure ma54_solve2_double
  end interface

  interface ma54_forward1
      module procedure ma54_forward1_double
  end interface

  interface ma54_forward2
      module procedure ma54_forward2_double
  end interface

  interface ma54_back1
      module procedure ma54_back1_double
  end interface

  interface ma54_back2
      module procedure ma54_back2_double
  end interface

  interface ma54_diag
      module procedure ma54_diag_double
  end interface

contains

subroutine ma54_to_block_double(n,p,nb,ap,buf,info,partial)
! Rearrange a packed triangular matrix to blocked hybrid format
   integer, intent (in) :: n ! Specifies the matrix order
   integer, intent (in) :: p ! Column p is be at the end of a block
   integer, intent (in) :: nb ! Block size for the blocked hybrid format.
   real (wp), intent (inout) :: ap((n*(n+one))/2)
!                Holds the matrix, which is rearranged to the new format.
   real (wp) :: buf(nb*int(n,long)) ! Work array.
   integer, intent (out) :: info ! set to unity of these values:
!                               0 Successful solution.
!                             < 0 Failure
   logical, optional, intent(in) :: partial ! Ignored.

!  Local variables
   integer :: jb ! Length of block
   integer :: j1 ! First column of block
   integer :: i, j ! Row, column index within the block column
   integer(long) :: ijap, ijbuf ! Positions of entry in ap and buf
   integer(long)  :: ac ! Position in ap of start of block column
   integer :: l ! Length of row or column
   integer :: lb ! Length of block column

   intrinsic min

   info = 0
   if (p > n)  info = -3
   if (n < 0)  info = -1
   if (p < 0)  info = -2
   if (nb < 1) info = -5
   if (info/=0 .or. n==0) return

   ijap = 1
! Main loop over block columns
   do j1 = 1,p,nb
      jb = min(nb,p-j1+1)
      lb = n - j1 + 1
      ac = ijap
! Copy to the buffer
      ijbuf = 1
      do j = 1, jb
         l = lb - j + 1
         call dcopy(l,ap(ijap),1,buf(ijbuf),1)
         ijbuf = ijbuf + lb + 1
         ijap = ijap + l
      end do

! Copy array back from buffer
      ijap = ac
      do i = 1, lb
         l = min(i,jb)
         call dcopy(l,buf(i),lb,ap(ijap),1)
        ijap = ijap + l
      end do
   end do

end subroutine ma54_to_block_double

subroutine ma54_factor_double(n,p,nb,ap,buf,info,n_threads)
! Cholesky part factorization of a symmetric positive-definite matrix held
! in packed blocked hybrid format.

!$     use omp_lib

  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: p ! Factorize to column p
  integer, intent (in) :: nb ! Block size for the blocked hybrid format.
  real (wp), intent (inout) :: ap(0:(n*(n+one))/2-1) ! Holds the matrix
!              A in packed blocked hybrid format. The ma54_factor of A is
!              contained in packed blocked hybrid format on exit.
  real (wp) :: buf(nb*int(n,long)) ! Work array.
  integer, intent (out) :: info ! set to unity of these values:
!                       0 Successful solution.
!                    /= 0 Failure:
!                     > 0. The leading minor of order info is not positive
!                       definite, and the factorization could not be completed.
   integer, intent (inout), optional :: n_threads ! Controls the number of
!              threads used in OpenMP parallel regions. Its value may be changed
!              by another thread during the execution of this subroutine.
!              0  Execute on omp_get_num_threads() threads.
!              1  Execute on a single thread without reference to OpenMP.
!             >1  Execute on n_threads threads.
!             Absence is treated as equivalent to presence with the value 0.

!  Local variables
    integer i ! Temporary variable
    integer ib ! Length of block row
    integer(long) ii ! Position in ap of diagonal entry
    integer j ! First column of current block column
    integer jb ! Number of columns in current block column
    integer(long) jbt ! Size of packed triangle of order jb
    integer(long) jj ! Position in ap of start of current diagonal block
    integer(long) jk ! Position in ap of start of current block
    integer kb ! Number of columns in earlier block column
    integer k ! First column of current block
    integer l ! Length of column
    integer lb ! Length of block column
    integer(long) kk ! Position in buf of diagonal entry
    integer nt  ! Number of threads

    intrinsic min

    info = 0
    if (p > n)  info = -3
    if (n < 0)  info = -1
    if (p < 0)  info = -2
    if (nb <= 0) info = -5
    if (info/=0 .or. n==0) return

    jj = 0
! Loop over leading block columns
    do j = 1, p, nb
      jb = min(nb,p-j+1)
      jbt = (jb*(jb+one))/2
      ii = jj
      kk = 1
! Move diagonal block to buf in full upper format
      do i = 1, jb
! Transfer i-th row of ap to i-th col of buf, if n > nb
        call dcopy(i,ap(ii),1,buf(kk),1)
        ii = ii + i
        kk = kk + jb
      end do

! Choose the number of threads for the next parallel loop
      nt = 1
!$    nt = omp_get_max_threads()
      if( nt>1 ) then
         if ( present(n_threads) ) then
            if( n_threads>0 ) nt = n_threads
         end if
      end if
!$    if (nt>0) call omp_set_num_threads(nt)

! Apply the previous elimination operations to the block column
!$omp parallel default(none) &
!$omp private(i, k, jk, kb, ib) &
!$omp shared(ap, j, jb, n, nb, jj, jbt, p, buf, info, lb)
! Note that the implied barrier at end of single avoids write-write
! data races in the syrk updates, allowing a NOWAIT on the gemm loop
      jk = 0
      do k = 1, min(j-1,p), nb ! Loop over block columns
        kb = min(nb,p-k+1)
        jk = jk + (j-k)*int(kb,long) - (kb*(kb-one))/2
! Update the diagonal block, held in buf
!$omp single
        call dsyrk('U','T',jb,kb,-unity,ap(jk),kb,unity,buf,jb)
!$omp end single
!$omp do
        do i = j+jb, n, nb
           ib = min(n+1-i, nb)
           call dgemm('T','N',jb,ib,kb,-unity,ap(jk),kb, &
              ap(jk+kb*int(i-j,long)),kb,unity, &
              ap(jj+jbt+jb*int(i-j-jb)),jb)
        end do
!$omp end do nowait
        jk = jk + (n+one-j)*kb
      end do
! Cholesky factorization of diagonal block
!$omp single
        call dkcf(jb,buf,jb,i)
!       call dpotrf('u',jb,buf,jb,i)
        if (i>0) info = i + j - 1
!$omp end single
! Important: implied barrier at end of single flushes info
        if(info.eq.0) then
! Calculate the sub-diagonal blocks of L in the current block column
!$omp do
          do i = j+jb, n, nb
            ib = min(n+1-i, nb)
            call dtrsm('L','U','T','N',jb,ib,unity,buf,jb, &
              ap(jj+jbt+jb*int(i-j-jb)),jb)
          end do
!$omp end do
        end if
!$omp end parallel
      if(info.ne.0) return

      ii = jj
      kk = 1
! Move  diagonal block back from buf to ap
      do i = 1, jb
        call dcopy(i,buf(kk),1,ap(ii),1)
        ii = ii + i
        kk = kk + jb
      end do

      jj = jj + jb*(n-j+one) - (jb*(jb-one))/2
    end do

! Loop over trailing block columns
    do j = p+1, n, nb
      jb = min(nb,n-j+1)
      jbt = (jb*(jb+one))/2
      ii = jj
      kk = 1
! Move block column to buf in rectangular format
      lb = n - j + 1
      do i = 1, jb
        l = lb - i + 1
        call dcopy(l,ap(ii),1,buf(kk),1)
        ii = ii + l
        kk = kk + lb + 1
      end do

! Choose the number of threads for the next parallel loop
      nt = 1
!$    nt = omp_get_max_threads()
      if( nt>1 ) then
         if ( present(n_threads) ) then
            if( n_threads>0 ) nt = n_threads
         end if
      end if
!$    if (nt>0) call omp_set_num_threads(nt)

! Apply the previous elimination operations to the block column
!$omp parallel default(none) &
!$omp private(i, k, jk, kb, ib) &
!$omp shared(ap, j, jb, n, nb, jj, jbt, p, buf, info, lb)
! Note that the implied barrier at end of single avoids write-write
! data races in the syrk updates, allowing a NOWAIT on the gemm loop
      jk = 0
      do k = 1, min(j-1,p), nb ! Loop over block columns
        kb = min(nb,p-k+1)
        jk = jk + (j-k)*int(kb,long) - (kb*(kb-one))/2
! Update the diagonal block, held in buf
!$omp single
        call dsyrk('L','T',jb,kb,-unity,ap(jk),kb,unity,buf,lb)
!$omp end single
! Update the off-diagonal part of the block column
!$omp do
        do i = j+jb, n, nb
          ib = min(n+1-i, nb)
          call dgemm('T','N',ib,jb,kb,-unity, &
            ap(jk+kb*int(i-j,long)),kb,ap(jk),kb, unity,buf(i-j+1),lb)
        end do
!$omp end do nowait

        jk = jk + (n+one-j)*kb
      end do
!$omp end parallel
      if(info.ne.0) return

      ii = jj
      kk = 1
! Move block column back from buf to ap
      do i = 1, jb
         l = lb - i + 1
         call dcopy(l,buf(kk),1,ap(ii),1)
         ii = ii + l
         kk = kk + lb + 1
      end do

      jj = jj + jb*(n-j+one) - (jb*(jb-one))/2
    end do

end subroutine ma54_factor_double


subroutine ma54_from_block_double(n,p,nb,ap,buf,info,partial)
! Rearrange blocked dhybrid packed matrix to packed triangular format
   integer, intent (in) :: n ! Specifies the matrix order
   integer, intent (in) :: p ! Column p is at the end of a block
   integer, intent (in) :: nb ! Block size  for the blocked hybrid format.
   real (wp), intent (inout) :: ap((n*(n+one))/2) ! Holds the matrix, which
!                       is rearranged to the packed triangular format.
   real (wp) :: buf(nb*int(n,long)) ! Work array.
   integer, intent (out) :: info ! set to unity of these values:
!                           0 Successful solution.
!                         < 0 Failure:
   logical, optional, intent(in) :: partial ! Ignored

   integer :: jb ! Column length of block
   integer :: j1 ! First column of block
   integer :: i, j ! Row, column index within the block
   integer(long) :: ijap, ijbuf ! Positions of ap(i1+i-1,j1+j-1), buf(i,j)
   integer(long) :: ac ! Position in ap of start of block
   integer :: l ! Length of row or column
   integer :: lb ! Length of block column
   intrinsic min

   info = 0
   if (p > n)  info = -3
   if (n < 0)  info = -1
   if (p < 0)  info = -2
   if (nb < 1) info = -5
   if (info/=0 .or. n==0) return

   ijap = 1
! Main loop over block columns
   do j1 = 1, p, nb
      jb = min(nb,p-j1+1)
      lb = n - j1 + 1

! Copy to buffer
      ac = ijap
      do i = 1, lb
         l = min(i,jb)
         call dcopy(l,ap(ijap),1,buf(i),lb)
         ijap = ijap + l
      end do

! Copy from the buffer
      ijap = ac
      ijbuf = 1
      do j = 1, jb
         l = lb - j + 1
         call dcopy(l,buf(ijbuf),1,ap(ijap),1)
         ijbuf = ijbuf + lb + 1
         ijap = ijap + l
      end do
   end do

end subroutine ma54_from_block_double


subroutine ma54_forward2_double(n,p,nbo,nrhs,ap,b,ldb,mbo,buf,info)
! Partial forward substitution for unity or more sets of equations, given a
!                partial Cholesky factorization in blocked hybrid format.
  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: p ! Number of columns
  integer, intent (in) :: nbo ! Block size for the blocked hybrid format.
  integer, intent (in) :: nrhs ! The number of right-hand sides.
  real (wp), intent (in) :: ap(0:(p*(n*2-p+one))/2-1) ! Holds the matrix.
  integer, intent (in) :: ldb ! The first dimension of the array b.
  real (wp), intent (inout) :: b(ldb,*) ! Holds the right-hand sides on
!                              entry and is overwritten by the solution.
  integer, intent (in) :: mbo ! Block size for the right-hand sides.
  real (wp) :: buf(int(p,long)*nbo + int(n,long)*mbo)
                              ! Work array (needed if nrhs >=4)
  integer, intent (out) :: info ! set to unity of these values:
!                         0 Successful solution.
!                       < 0 Failure

! Local variables
  integer(long) :: bufb ! Position in buf that holds a buffer for b
  integer :: i ! Row index
  integer :: ib ! Number of columns in current block column of b
  integer(long) :: ibp ! Current position in bufb
  integer(long) :: ii ! Position in ap of diagonal entry
  integer(long) :: ij ! Position in ap
  integer :: j ! Column index
  integer :: jb ! Number of columns in current block column of ap
  integer(long) :: jbp ! Current position in bufb
  integer(long) :: jbd ! nb*int(kb,long), increment for jbp
  integer(long) :: jd ! size of current trapezoid / increment for jj
  integer(long) :: jj ! Position in ap of start of current diagonal block
  integer :: k ! First column of current block
  integer :: kb ! ! Number of columns in current block column of b
  integer(long) :: kk ! Position in buf of diagonal entry
  integer(long) :: kks ! Current position in buf
  integer :: mb ! block size used for the blocked rhs:min(m,mbo)
  integer :: nb ! block size used for the blocked hybrid format: min(n,nbo)
  integer(long) :: nb2 ! nb*int(nb,long)
  integer(long) :: nbt ! Size of packed triangle of order nb
  integer(long) :: nb1 ! Size of packed triangle of order nb-1

  info = 0
  if (p>n) info = -3
  if (n<0) info = -1
  if (p<0) info = -2
  if (nrhs<0) info = -4
  if (nbo<1) info = -5
  if (ldb<max(1,n)) info = -6
  if (mbo<1) info = -7
  nb = min(n,nbo)
  mb = min(mbo,nrhs)
  if (info/=0 .or. n==0) return

  if (nrhs<4) then
    do j = 1, nrhs
      call ma54_forward1_double(n,p,nbo,ap,b(1,j),info)
    end do
    return
  end if

  bufb = p*int(nb,long)

  nb2 = nb*int(nb,long)
  nbt = (nb2+nb)/2
  nb1 = nbt - nb

!  Load the diagonal blocks (triangles) of ap into buf.
  kks = 1
  jd = n*int(nb,long) - nb1
  jj = 0
  do j = 1, p, nb
    jb = min(p-j+1,nb)
    ii = jj
    kk = kks
    do i = 1, jb
      call dcopy(i,ap(ii),1,buf(kk),1)
      ii = ii + i
      kk = kk + jb
    end do
    jj = jj + jd ! -> next triangle in ap
    jd = jd - nb2 ! next delta for jj
    kks = kks + nb2 ! -> next triangle of ap in buf
  end do ! do j

! Main loop over the block columns of b
  do k = 1, nrhs, mb
    kks = 1
    jd = n*int(nb,long) - nb1
    jj = 0
    kb = min(nrhs-k+1,mb)

! Copy kb RHS from  b(1:n,k:k+kb-1) to buffer bufb
    jbp = 1
    do i = 1, p, nb
      ib = min(p-i+1,nb)
      do j = k, k + kb - 1
        call dcopy(ib,b(i,j),1,buf(bufb+jbp),1)
        jbp = jbp + ib
      end do
    end do
    do i = p+1, n, nb
      ib = min(n-i+1,nb)
      do j = k, k + kb - 1
        call dcopy(ib,b(i,j),1,buf(bufb+jbp),1)
        jbp = jbp + ib
      end do
    end do

! Solve U'*Y=B
    jbd = nb*int(kb,long) ! delta for jbp
    jbp = 1
    do j = 1, p, nb
      jb = min(p-j+1,nb)
      call dtrsm('L','U','T','N',jb,kb,unity,buf(kks),jb,buf(bufb+jbp),jb)
      ibp = jbp + jb*int(kb,long) ! -> starting block of B
      ij = jj + (jb*(jb+one))/2 ! -> starting block of AP
      do i = j + nb, p, nb
        ib = min(p-i+1,nb)
        call dgemm('T','N',ib,kb,jb,-unity,ap(ij),jb,buf(bufb+jbp),jb,unity, &
          buf(bufb+ibp),ib)
        ibp = ibp + ib*int(kb,long) ! -> next block of B
        ij = ij + ib*int(jb,long) ! -> next block of AP
      end do ! do i
      do i = p+1, n, nb
        ib = min(n-i+1,nb)
        call dgemm('T','N',ib,kb,jb,-unity,ap(ij),jb,buf(bufb+jbp),jb,unity, &
          buf(bufb+ibp),ib)
        ibp = ibp + ib*int(kb,long) ! -> next block of B
        ij = ij + ib*int(jb,long) ! -> next block of AP
      end do ! do i
      jj = jj + jd
      jd = jd - nb2
      jbp = jbp + jbd ! -> next block of B
      kks = kks + nb2 ! -> next triangle of ap
    end do ! do j

!       Copy solution to b(1:n,k:k+kb+1)
    jbp = 1 ! -> to buf
    do i = 1, p, nb
      ib = min(p-i+1,nb)
      do j = k, k + kb - 1
        call dcopy(ib,buf(bufb+jbp),1,b(i,j),1)
        jbp = jbp + ib
      end do
    end do
    do i = p+1, n, nb
      ib = min(n-i+1,nb)
      do j = k, k + kb - 1
        call dcopy(ib,buf(bufb+jbp),1,b(i,j),1)
        jbp = jbp + ib
      end do
    end do
  end do ! do k

end subroutine ma54_forward2_double

subroutine ma54_forward1_double(n,p,nb,ap,b,info)
! Partial forward substitution for unity set of equations, given a
!          partial Cholesky factorization in blocked hybrid format.
  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: p ! Number of columns
  integer, intent (in) :: nb ! Block size for the blocked hybrid format.
  real (wp), intent (in) :: ap(0:(p*(n*2-p+one))/2-1) ! Holds the matrix.
  real (wp), intent (inout) :: b(n) ! b(1:n) holds the right-hand
!                   sides on entry and is overwritten by the solution.
  integer, intent (out) :: info ! set to unity of these values:
!          0 Successful solution.
!        < 0 The value of argument -info is not valid.

! Local variables
  integer j ! First column of current block column
  integer jb ! Number of columns in current block column
  integer(long) jbt ! Size of packed triangle of order jb
  integer(long) jd ! size of current trapezoid / increment for jj
  integer(long) jj ! Position in ap of start of current diagonal block
  integer(long) nb2 ! nb*int(nb,long)
  integer(long) nbt ! Size of packed triangle of order nb

! .. Executable Statements ..
  info = 0
  if (p>n) info = -3
  if (n<0) info = -1
  if (p<0) info = -2
  if (nb<1) info = -5

  if (info/=0 .or. n==0)  return
  jj = 0
  nb2 = nb*int(nb,long)
  nbt = (nb2+nb)/2
  jd = n*int(nb,long) - nbt + nb

!     Solve U'*Y=B
  do j = 1, p - nb, nb
    call dtpsv('U','T','N',nb,ap(jj),b(j),1)
    call dgemv('T',nb,n+1-j-nb,-unity,ap(jj+nbt),nb,b(j),1,unity,b(j+nb),1)
    jj = jj + jd
    jd = jd - nb2
  end do
  jb = p - j + 1
  jbt = (jb*(jb+one))/2
  call dtpsv('U','T','N',jb,ap(jj),b(j),1)
  if (n>=j+jb) then
    call dgemv('T',jb,n+1-j-jb,-unity,ap(jj+jbt),jb,b(j),1,unity,b(j+jb),1)
  end if
end subroutine ma54_forward1_double


subroutine ma54_back2_double(n,p,nb,nrhs,ap,b,ldb,mb,buf,info)
! Partial back substitution for unity or more sets of equations, given a
!              partial Cholesky factorization in blocked hybrid format.
  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: p ! Number of columns
  integer, intent (in) :: nb ! Block size for the blocked hybrid format.
  integer, intent (in) :: nrhs ! The number of right-hand sides.
  real (wp), intent (in) :: ap(0:(p*(n*2-p+one))/2-1) ! Holds the matrix.
  integer, intent (in) :: ldb ! The first dimension of the array b.
  real (wp), intent (inout) :: b(ldb,*) ! Holds the right-hand sides on
!                             entry and is overwritten by the solution.
  integer, intent (in) :: mb ! Block size for the right-hand sides.
  real (wp) :: buf(p*int(nb,long) + n*int(mb,long)) ! Work array
  integer, intent (out) :: info ! set to unity of these values:
!                 0 Successful solution.
!               < 0 Failure

! Local variables
  integer(long) :: bufb
  integer :: i ! Row index
  integer :: ib ! Number of columns in current block column of b
  integer(long) :: ibp ! Current position in bufb
  integer(long) :: ii ! Position in ap of diagonal entry
  integer(long) :: ij ! Position in ap
  integer :: j ! Column index
  integer :: jb ! Number of columns in current block column of ap
  integer(long) :: jbp ! Current position in bufb
  integer(long) :: jbd ! nb*int(kb,long), increment for jbp
  integer(long) :: jd ! size of current trapezoid / increment for jj
  integer(long) :: jj ! Position in ap of start of current diagonal block
  integer :: k ! First column of current block
  integer :: kb ! ! Number of columns in current block column of b
  integer(long) :: kk ! Position in buf of diagonal entry
  integer(long) :: kks ! Current position in buf
  integer(long) :: nb2 ! nb*int(nb,long)
  integer(long) :: nbt ! Size of packed triangle of order nb
  integer(long) :: nb1 ! Size of packed triangle of order nb-1

  info = 0
  if (p>n) info = -3
  if (n<0) info = -1
  if (p<0) info = -2
  if (nrhs<0) info = -4
  if (nb<1) info = -5
  if (ldb<max(1,n)) info = -6
  if (mb<1) info = -7
  if (info/=0 .or. n==0) return

  if (nrhs<4) then
    do j = 1, nrhs
      call ma54_back1_double(n,p,nb,ap,b(1,j),info)
    end do
    return
  end if

  bufb = p*int(nb,long)

  nb2 = nb*int(nb,long)
  nbt = (nb2+nb)/2
  nb1 = nbt - nb

!  Load the diagonal blocks (triangles) of ap into buf.
  kks = 1
  jd = n*int(nb,long) - nb1
  jj = 0
  do j = 1, p, nb
    jb = min(p-j+1,nb)
    ii = jj
    kk = kks
    do i = 1, jb
      call dcopy(i,ap(ii),1,buf(kk),1)
      ii = ii + i
      kk = kk + jb
    end do
    jj = jj + jd ! -> next triangle in ap
    jd = jd - nb2 ! next delta for jj
    kks = kks + nb2 ! -> next triangle of ap in buf
  end do ! do j

!     Main loop over the block columns of b
  do k = 1, nrhs, mb
    kb = min(nrhs-k+1,mb)

! Copy kb RHS from  b(1:n,k:k+kb-1) to buffer bufb
    jbp = 1
    do i = 1, p, nb
      ib = min(p-i+1,nb)
      do j = k, k + kb - 1
        call dcopy(ib,b(i,j),1,buf(bufb+jbp),1)
        jbp = jbp + ib
      end do
    end do
    do i = p+1, n, nb
      ib = min(n-i+1,nb)
      do j = k, k + kb - 1
        call dcopy(ib,b(i,j),1,buf(bufb+jbp),1)
        jbp = jbp + ib
      end do
    end do

! Solve U*X=Y
    j = ((p-1)/nb)*int(nb,long)
    jj = (j*(n+n-j+one))/2
    jbp = 1+j*int(kb,long)
    jbd = nb*int(kb,long) ! delta for jbp
    jd = (nb*(2*(n-j)+nb+one))/2
    kks = 1+j*int(nb,long)
    do j = j+1, 1, -nb
      jb = min(p-j+1,nb)
      ij = jj + (jb*(jb+one))/2 ! -> starting block of AP
      ibp = jbp + jb*int(kb,long) ! -> starting update block of B
      do i = j + nb, p, nb
        ib = min(p-i+1,nb)
        call dgemm('N','N',jb,kb,ib,-unity,ap(ij),jb,buf(bufb+ibp),ib,unity, &
          buf(bufb+jbp),jb)
        ibp = ibp + ib*int(kb,long) ! -> next block of B
        ij = ij + ib*int(jb,long) ! -> next block of AP
      end do ! do i
      do i = p+1, n, nb
        ib = min(n-i+1,nb)
        call dgemm('N','N',jb,kb,ib,-unity,ap(ij),jb,buf(bufb+ibp),ib,unity, &
          buf(bufb+jbp),jb)
        ibp = ibp + ib*int(kb,long) ! -> next block of B
        ij = ij + ib*int(jb,long)  ! -> next block of AP
      end do ! do i
      call dtrsm('L','U','N','N',jb,kb,unity,buf(kks),jb,buf(bufb+jbp),jb)
      jj = jj - jd ! -> next triangle of AP
      jd = jd + nb2 ! next delta for jj
      jbp = jbp - jbd ! -> current block of B
      kks = kks - nb2 ! -> current triangle of ap
    end do ! do j

! Copy solution to b(1:n,k:k+kb+1)
    jbp = 1 ! -> to buf
    do i = 1, p, nb
      ib = min(p-i+1,nb)
      do j = k, k + kb - 1
        call dcopy(ib,buf(bufb+jbp),1,b(i,j),1)
        jbp = jbp + ib
      end do
    end do
    do i = p+1, n, nb
      ib = min(n-i+1,nb)
      do j = k, k + kb - 1
        call dcopy(ib,buf(bufb+jbp),1,b(i,j),1)
        jbp = jbp + ib
      end do
    end do
  end do ! do k

end subroutine ma54_back2_double

subroutine ma54_back1_double(n,p,nb,ap,b,info)
! Partial back substitution for unity set of equations, given a
!             partial Cholesky factorization in blocked hybrid format.
  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: p ! Number of columns
  integer, intent (in) :: nb ! Block size for the blocked hybrid format.
  real (wp), intent (in) :: ap(0:(p*(n*2-p+one))/2-1) ! Holds the matrix.
  real (wp), intent (inout) :: b(n) ! b(1:n) holds the right-hand
! sides on entry and is overwritten by the solution.
  integer, intent (out) :: info ! set to unity of these values:
!                   0 Successful solution.
!                 < 0 The value of argument -info is not valid.

  integer i ! Temporary variable
  integer j ! First column of current block column
  integer jb ! Number of columns in current block column
  integer(long) jbt ! Size of packed triangle of order jb
  integer(long) jd ! size of current trapezoid / increment for jj
  integer(long) jj ! Position in ap of start of current diagonal block
  integer(long) nb2 ! nb*int(nb,long)
  integer(long) nbt ! Size of packed triangle of order nb

  info = 0
  if (p>n) info = -3
  if (n<0) info = -1
  if (p<0) info = -2
  if (nb<1) info = -5
  if (info/=0 .or. n==0) return

  jj = 0
  nb2 = nb*int(nb,long)
  nbt = (nb2+nb)/2
  jd = n*int(nb,long) - nbt + nb
  do j = 1, p - nb, nb
    jj = jj + jd
    jd = jd - nb2
  end do
  jb = p - j + 1
  jbt = (jb*(jb+one))/2
  if (n>=j+jb) then
    call dgemv('N',jb,n+1-j-jb,-unity,ap(jj+jbt),jb,b(j+jb),1,unity,b(j),1)
  end if
  call dtpsv('U','N','N',jb,ap(jj),b(j),1)
  i = j - nb
  do j = i, 1, -nb
    jd = jd + nb2
    jj = jj - jd
    call dgemv('N',nb,n+1-j-nb,-unity,ap(jj+nbt),nb,b(j+nb),1,unity,b(j),1)
    call dtpsv('U','N','N',nb,ap(jj),b(j),1)
  end do

end subroutine ma54_back1_double

subroutine ma54_solve1_double(n,nb,ap,b,info)
! Solve unity or more sets of equations, given the Cholesky factorization
!                  of its matrix in blocked hybrid format.
! For the character argument, the case is insignificant.
  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: nb ! Block size for the blocked hybrid format.
  real (wp), intent (in) :: ap(0:(n*(n+one))/2-1) ! Holds the matrix.
  real (wp), intent (inout) :: b(n) !  holds the right-hand side on entry
!            and is overwritten by the solution.
  integer, intent (out) :: info ! set to unity of these values:
!                0 Successful solution.
!              < 0 Failure

    integer i ! Temporary variable
    integer j ! First column of current block column
    integer jb ! Number of columns in current block column
    integer(long) jd ! size of current trapezoid / increment for jj
    integer(long) jj ! Position in ap of start of current diagonal block
    integer(long) nb2 ! nb*int(nb,long)
    integer(long) nbt ! Size of packed triangle of order nb

  info = 0
  if (n<0) info = -1
  if (nb<1) info = -5

  if (info/=0 .or. n==0) return

    jj = 0
    nb2 = nb*int(nb,long)
    nbt = (nb2+nb)/2
    jd = n*int(nb,long) - nbt + nb

!     Solve U'*Y=B
    do j = 1, n - nb, nb
      call dtpsv('U','T','N',nb,ap(jj),b(j),1)
      call dgemv('T',nb,n+1-j-nb,-unity,ap(jj+nbt),nb,b(j),1,unity,b(j+nb),1)
      jj = jj + jd
      jd = jd - nb2
    end do
    jb = n - j + 1
    call dtpsv('U','T','N',jb,ap(jj),b(j),1)

!     Solve U*X=Y
    call dtpsv('U','N','N',jb,ap(jj),b(j),1)
    i = j - nb
    do j = i, 1, -nb
      jd = jd + nb2
      jj = jj - jd
      call dgemv('N',nb,n+1-j-nb,-unity,ap(jj+nbt),nb,b(j+nb),1,unity,b(j),1)
      call dtpsv('U','N','N',nb,ap(jj),b(j),1)
    end do
 end subroutine ma54_solve1_double

subroutine ma54_solve2_double(n,nb,nrhs,ap,b,ldb,mb,buf,info)
! Solve unity or more sets of equations, given the Cholesky factorization
!              of its matrix in blocked hybrid format.
  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: nb ! Block size for the blocked hybrid format.
  real (wp), intent (in) :: ap(0:(n*(n+one))/2-1) ! Holds the matrix.
  integer, intent (in) :: ldb ! The first dimension of the array b.
  integer, intent (in) :: nrhs ! The number of right-hand sides.
  real (wp), intent (inout) :: b(ldb,*) ! b(1:n,1:nrhs) holds the
!             right-hand sides on entry and is overwritten by the solution.
  integer, intent (in) :: mb ! Block size for the right-hand sides.
  real (wp) :: buf(*) ! work array
  integer, intent (out) :: info ! set to unity of these values:
!                                0 Successful solution.
!                              < 0 Failure

    integer(long) :: bufb
    integer :: i ! Row index
    integer :: ib ! Number of columns in current block column of b
    integer(long) :: ibp ! Current position in bufb
    integer(long) :: ii ! Position in ap of diagonal entry
    integer(long) :: ij ! Position in ap
    integer :: j ! Column index
    integer :: jb ! Number of columns in current block column of ap
    integer(long) :: jbp ! Current position in bufb
    integer(long) :: jbd ! nb*int(kb,long), increment for jbp
    integer(long) :: jd ! size of current trapezoid / increment for jj
    integer(long) :: jj ! Position in ap of start of current diagonal block
    integer :: k ! First column of current block
    integer :: kb ! ! Number of columns in current block column of b
    integer(long) :: kk ! Position in buf of diagonal entry
    integer(long) :: kks ! Current position in buf
    integer(long) :: nb2 ! nb*int(nb,long)
    integer(long) :: nbt ! Size of packed triangle of order nb
    integer(long) :: nb1 ! Size of packed triangle of order nb-1
    intrinsic min

  info = 0
  if (n<0) info = -1
  if (nrhs<0) info = -4
  if (nb<1) info = -5
  if (ldb<max(1,n)) info = -6
  if (mb<1) info = -7

  if (info/=0 .or. n==0 .or. nrhs==0) return

! If there are few rhs, solve each separately
    if (nrhs<4) then
      do j = 1, nrhs
        call ma54_solve1_double(n,nb,ap,b(1,j),info)
      end do
      return
    end if

    bufb = n*int(nb,long)

    nb2 = nb*int(nb,long)
    nbt = (nb2+nb)/2
    nb1 = nbt - nb

!     Load the diagonal blocks (triangles) of ap into buf.
    kks = 1 ! note that kks is not restored to 1 here.
    jd = n*int(nb,long) - nb1 ! jd is not restored to this value here.
    jj = 0 ! is not restored to zero here.
    do j = 1, n, nb
      jb = min(n-j+1,nb)
      ii = jj
      kk = kks
      do i = 1, jb
        call dcopy(i,ap(ii),1,buf(kk),1)
        ii = ii + i
        kk = kk + jb
      end do
      jj = jj + jd ! -> next triangle in ap
      jd = jd - nb2 ! next delta for jj
      kks = kks + nb2 ! -> next triangle of ap in buf
    end do ! do j

!     Main loop over the block columns of b
    kks = 1 ! is changed but restored to 1
    jd = n*int(nb,long) - nb1 ! is changed but restored to this value
    jj = 0 ! is changed but restored to 0
    do k = 1, nrhs, mb
      kb = min(nrhs-k+1,mb)

!       Copy kb RHS from  b(1:n,k:k+kb-1) to buffer bufb
      jbp = 1 ! -> to bufb
      do i = 1, n, nb
        ib = min(n-i+1,nb)
        do j = k, k + kb - 1
          call dcopy(ib,b(i,j),1,buf(bufb+jbp),1)
          jbp = jbp + ib
        end do
      end do

!      Solve U'*Y=B
      jbd = nb*int(kb,long) ! delta for jbp
      jbp = 1 ! is changed but restored to 1
      do j = 1, n - nb, nb
        call dtrsm('L','U','T','N',nb,kb,unity,buf(kks),nb,buf(bufb+jbp),nb)
        ibp = jbp + jbd ! -> starting block of B
        ij = jj + nbt ! -> starting block of AP
        do i = j + nb, n, nb
          jb = min(n-i+1,nb)
          call dgemm('T','N',jb,kb,nb,-unity,ap(ij),nb,buf(bufb+jbp),nb,unity, &
            buf(bufb+ibp),jb)
          ibp = ibp + jbd ! -> next block of B
          ij = ij + nb2 ! -> next block of AP
        end do ! do i
        jj = jj + jd
        jd = jd - nb2
        jbp = jbp + jbd ! -> next block of B
        kks = kks + nb2 ! -> next triangle of ap
      end do ! do j
      jb = n - j + 1
      call dtrsm('L','U','T','N',jb,kb,unity,buf(kks),jb,buf(bufb+jbp),jb)

!       Solve U*X=Y
      call dtrsm('L','U','N','N',jb,kb,unity,buf(kks),jb,buf(bufb+jbp),jb)
      do j = j - nb, 1, -nb
        jd = jd + nb2
        jj = jj - jd ! -> next triangle of AP
        ij = jj + nbt ! -> starting block of AP
        ibp = jbp ! -> starting update block of B
        jbp = jbp - jbd ! -> current block of B
        kks = kks - nb2 ! -> current triangle of ap
        do i = j + nb, n, nb
          jb = min(n-i+1,nb)
          call dgemm('N','N',nb,kb,jb,-unity,ap(ij),nb,buf(bufb+ibp),jb,unity, &
            buf(bufb+jbp),nb)
          ibp = ibp + jbd ! -> next block of B
          ij = ij + nb2 ! -> next block of AP
        end do ! do i
        call dtrsm('L','U','N','N',nb,kb,unity,buf(kks),nb,buf(bufb+jbp),nb)
      end do ! do j

!       Copy solution to b(1:n,k:k+kb+1)
      jbp = 1 ! -> to buf
      do i = 1, n, nb
        ib = min(n-i+1,nb)
!         copy buf(is:is+ib+k-1) to B(i:i+ib-1,0:k-1); K <= NB
        do j = k, k + kb - 1
          call dcopy(ib,buf(bufb+jbp),1,b(i,j),1)
          jbp = jbp + ib
        end do
      end do
    end do ! do k

  end subroutine ma54_solve2_double

 subroutine ma54_diag_double(n,p,nb,ap,diag)
! Extract the diagonal of a partial Cholesky factorization in blocked
! hybrid format.
  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: p ! Number of columns
  integer, intent (in) :: nb ! Block size for the blocked hybrid format.
  real (wp), intent (in) :: ap((n*(n+one))/2) ! Holds the matrix.
  real (wp), intent (out) :: diag(p) ! Diagonal of factorization.

  integer js ! First column of current block column
  integer j  ! Column index within the block
  integer(long) k  ! Index within ap
  integer l  ! Length of rectangular part of the block column

  k = 0
  l = n-nb
  do js = 1, p, nb
     do j = 1, min(nb,p-js+1)
       k = k + j
       diag(js+j-1) = ap(k)
     end do
     k = k + l*int(nb,long)
     l = l - nb
  end do
end subroutine ma54_diag_double


subroutine dkcf(n,a,lda,info)
! Kernel subroutine for Cholesky factorization
  integer, intent (in) :: n ! Specifies the matrix order.
  integer, intent (in) :: lda ! Leading extent of array a.
! The inequality lda >= n must hold.
! For efficiency, the value n for lda is preferable.
  real (wp), intent (inout) :: a(0:lda-1,0:n-1) ! The upper-triangular part,
! a(i,j), i<=j<=n, must be set to hold the upper-triangular part
! of the matrix and is overwritten by the Cholesky ma54_factor.
! For efficient execution, a(:,1:n) should fit into
! level-1 cache.
  integer, intent (out) :: info ! set to unity of these values:
! 0 Successful solution.
! i>0. Pivot i is not positive.
! .. Locals ..
  real (wp), parameter :: zero = 0.0_wp, unity = 1.0_wp
  integer :: i ! Temporary variable
  integer :: j ! Temporary variable
  integer :: k ! Temporary variable
  integer :: nr ! Temporary variable
  real (wp) :: rd0 ! scalar to hold reciprocal of diagonal A(K,K) for 0<=K<N
  real (wp) :: t00, t01, t02, t03 ! regs to hold a(j,i:i+3)
  real (wp) :: ai0, ai1, ai2, ai3 ! regs to hold a(k,i:i+3)
!     or A(K:K+3,J) or A(K:K+3,J+1)

!     this code uses 8 FP regs



!     This is a 1 x 4 blocking kernel
!     It is ONLY optimized for the case MOD(N,4) = 0

  info = 0
  nr = mod(n,4)

!     Main loop: ma54_factor A(0:n-nr-1,0:n-nr-1) = U'*U

  do j = 0, n - nr - 4, 4

!       Process row J = A(J,J:N-1)

!       Process A(J  ,J  )

    t00 = zero
    t01 = zero
    t02 = zero
    t03 = zero

    do k = 0, j - 1, 4
      ai0 = a(k,j)
      ai1 = a(k+1,j)
      ai2 = a(k+2,j)
      ai3 = a(k+3,j)
      ai0 = ai0*ai0
      ai1 = ai1*ai1
      ai2 = ai2*ai2
      ai3 = ai3*ai3
      t00 = t00 - ai0
      t01 = t01 - ai1
      t02 = t02 - ai2
      t03 = t03 - ai3
    end do
    ai0 = a(j,j)
    t00 = t00 + t01
    t02 = t02 + t03
    t00 = t00 + t02
    ai0 = ai0 + t00

    if (ai0<=zero) go to 2
    ai0 = sqrt(ai0)
    a(j,j) = ai0
    rd0 = unity/ai0


!       process 1 by 3 rectangle A(J,J+1:J+3)

    t00 = a(j,j+1)
    t01 = a(j,j+2)
    t02 = a(j,j+3)
    do k = 0, j - 1
      ai0 = a(k,j+1)
      ai1 = a(k,j+2)
      ai2 = a(k,j+3)
      ai0 = ai0*a(k,j)
      ai1 = ai1*a(k,j)
      ai2 = ai2*a(k,j)
      t00 = t00 - ai0
      t01 = t01 - ai1
      t02 = t02 - ai2
    end do

!       scale and store A(J,J+1:J+3)

    t00 = t00*rd0
    t01 = t01*rd0
    t02 = t02*rd0
    a(j,j+1) = t00
    a(j,j+2) = t01
    a(j,j+3) = t02

!       Process the remainder of rows J = A(J,J+4:N-1)
!       using 1 by 4 rectangles

    do i = j + 4, n - 4, 4
      t00 = a(j,i)
      t01 = a(j,i+1)
      t02 = a(j,i+2)
      t03 = a(j,i+3)
      do k = 0, j - 1
        ai0 = a(k,i)
        ai1 = a(k,i+1)
        ai2 = a(k,i+2)
        ai3 = a(k,i+3)
        ai0 = ai0*a(k,j)
        ai1 = ai1*a(k,j)
        ai2 = ai2*a(k,j)
        ai3 = ai3*a(k,j)
        t00 = t00 - ai0
        t01 = t01 - ai1
        t02 = t02 - ai2
        t03 = t03 - ai3
      end do

!         scale and store A(J,I+0:I+3)

      t00 = t00*rd0
      t01 = t01*rd0
      t02 = t02*rd0
      t03 = t03*rd0
      a(j,i) = t00
      a(j,i+1) = t01
      a(j,i+2) = t02
      a(j,i+3) = t03
    end do

!       Process row J+1 = A(J+1,J+1:N-1)

!       Process A(J+1,J+1)

    t00 = zero
    t01 = zero
    t02 = zero
    t03 = zero

    do k = 0, j - 1, 4
      ai0 = a(k,j+1)
      ai1 = a(k+1,j+1)
      ai2 = a(k+2,j+1)
      ai3 = a(k+3,j+1)
      ai0 = ai0*ai0
      ai1 = ai1*ai1
      ai2 = ai2*ai2
      ai3 = ai3*ai3
      t00 = t00 - ai0
      t01 = t01 - ai1
      t02 = t02 - ai2
      t03 = t03 - ai3
    end do
    ai0 = a(j,j+1)
    ai0 = ai0*ai0
    t00 = t00 - ai0

    ai0 = a(j+1,j+1)
    t00 = t00 + t01
    t02 = t02 + t03
    t00 = t00 + t02
    ai0 = ai0 + t00

    if (ai0<=zero) go to 3
    ai0 = sqrt(ai0)
    a(j+1,j+1) = ai0
    rd0 = unity/ai0


!       process 1 by 2 rectangle A(J+1,J+2:J+3)

    t00 = a(j+1,j+2)
    t01 = a(j+1,j+3)
    do k = 0, j
      ai0 = a(k,j+2)
      ai1 = a(k,j+3)
      ai0 = ai0*a(k,j+1)
      ai1 = ai1*a(k,j+1)
      t00 = t00 - ai0
      t01 = t01 - ai1
    end do

!       scale and store A(J+1,J+2:J+3)

    t00 = t00*rd0
    t01 = t01*rd0
    a(j+1,j+2) = t00
    a(j+1,j+3) = t01

!       Process the remainder of row J+1 = A(J+1,J+4:N-1)
!       using 1 by 4 rectangles

    do i = j + 4, n - 4, 4
      t00 = a(j+1,i)
      t01 = a(j+1,i+1)
      t02 = a(j+1,i+2)
      t03 = a(j+1,i+3)
      do k = 0, j
        ai0 = a(k,i)
        ai1 = a(k,i+1)
        ai2 = a(k,i+2)
        ai3 = a(k,i+3)
        ai0 = ai0*a(k,j+1)
        ai1 = ai1*a(k,j+1)
        ai2 = ai2*a(k,j+1)
        ai3 = ai3*a(k,j+1)
        t00 = t00 - ai0
        t01 = t01 - ai1
        t02 = t02 - ai2
        t03 = t03 - ai3
      end do

!         scale and store A(J+1,I+0:I+3)

      t00 = t00*rd0
      t01 = t01*rd0
      t02 = t02*rd0
      t03 = t03*rd0
      a(j+1,i) = t00
      a(j+1,i+1) = t01
      a(j+1,i+2) = t02
      a(j+1,i+3) = t03
    end do

!       Process row J+2 = A(J+2,J+2:N-1)

!       Process A(J+2,J+2)

    t00 = zero
    t01 = zero
    t02 = zero
    t03 = zero

    do k = 0, j - 1, 4
      ai0 = a(k,j+2)
      ai1 = a(k+1,j+2)
      ai2 = a(k+2,j+2)
      ai3 = a(k+3,j+2)
      ai0 = ai0*ai0
      ai1 = ai1*ai1
      ai2 = ai2*ai2
      ai3 = ai3*ai3
      t00 = t00 - ai0
      t01 = t01 - ai1
      t02 = t02 - ai2
      t03 = t03 - ai3
    end do
    ai0 = a(j,j+2)
    ai1 = a(j+1,j+2)
    ai0 = ai0*ai0
    ai1 = ai1*ai1
    t00 = t00 - ai0
    t01 = t01 - ai1

    ai0 = a(j+2,j+2)
    t00 = t00 + t01
    t02 = t02 + t03
    t00 = t00 + t02
    ai0 = ai0 + t00

    if (ai0<=zero) go to 4
    ai0 = sqrt(ai0)
    a(j+2,j+2) = ai0
    rd0 = unity/ai0


!       process 1 by 1 rectangle A(J+2,J+3)

    t00 = a(j+2,j+3)
    do k = 0, j + 1
      ai0 = a(k,j+3)
      ai0 = ai0*a(k,j+2)
      t00 = t00 - ai0
    end do

!       scale and store A(J+2,J+3)

    t00 = t00*rd0
    a(j+2,j+3) = t00

!       Process the remainder of row J+2 = A(J+2,J+4:N-1)
!       using 1 by 4 rectangles

    do i = j + 4, n - 4, 4
      t00 = a(j+2,i)
      t01 = a(j+2,i+1)
      t02 = a(j+2,i+2)
      t03 = a(j+2,i+3)
      do k = 0, j + 1
        ai0 = a(k,i)
        ai1 = a(k,i+1)
        ai2 = a(k,i+2)
        ai3 = a(k,i+3)
        ai0 = ai0*a(k,j+2)
        ai1 = ai1*a(k,j+2)
        ai2 = ai2*a(k,j+2)
        ai3 = ai3*a(k,j+2)
        t00 = t00 - ai0
        t01 = t01 - ai1
        t02 = t02 - ai2
        t03 = t03 - ai3
      end do

!         scale and store A(J+2,I+0:I+3)

      t00 = t00*rd0
      t01 = t01*rd0
      t02 = t02*rd0
      t03 = t03*rd0
      a(j+2,i) = t00
      a(j+2,i+1) = t01
      a(j+2,i+2) = t02
      a(j+2,i+3) = t03
    end do

!       Process row J+3 = A(J+3,J+3:N-1)

!       Process A(J+3,J+3)

    t00 = zero
    t01 = zero
    t02 = zero
    t03 = zero

    do k = 0, j - 1, 4
      ai0 = a(k,j+3)
      ai1 = a(k+1,j+3)
      ai2 = a(k+2,j+3)
      ai3 = a(k+3,j+3)
      ai0 = ai0*ai0
      ai1 = ai1*ai1
      ai2 = ai2*ai2
      ai3 = ai3*ai3
      t00 = t00 - ai0
      t01 = t01 - ai1
      t02 = t02 - ai2
      t03 = t03 - ai3
    end do
    ai0 = a(j,j+3)
    ai1 = a(j+1,j+3)
    ai2 = a(j+2,j+3)
    ai0 = ai0*ai0
    ai1 = ai1*ai1
    ai2 = ai2*ai2
    t00 = t00 - ai0
    t01 = t01 - ai1
    t02 = t02 - ai2

    ai0 = a(j+3,j+3)
    t00 = t00 + t01
    t02 = t02 + t03
    t00 = t00 + t02
    ai0 = ai0 + t00

    if (ai0<=zero) go to 5
    ai0 = sqrt(ai0)
    a(j+3,j+3) = ai0
    rd0 = unity/ai0


!       Process the remainder of row J+3 = A(J+3,J+4:N-1)
!       using 1 by 4 rectangles

    do i = j + 4, n - 4, 4
      t00 = a(j+3,i)
      t01 = a(j+3,i+1)
      t02 = a(j+3,i+2)
      t03 = a(j+3,i+3)
      do k = 0, j + 2
        ai0 = a(k,i)
        ai1 = a(k,i+1)
        ai2 = a(k,i+2)
        ai3 = a(k,i+3)
        ai0 = ai0*a(k,j+3)
        ai1 = ai1*a(k,j+3)
        ai2 = ai2*a(k,j+3)
        ai3 = ai3*a(k,j+3)
        t00 = t00 - ai0
        t01 = t01 - ai1
        t02 = t02 - ai2
        t03 = t03 - ai3
      end do

!         scale and store A(J+3,I+0:I+3)

      t00 = t00*rd0
      t01 = t01*rd0
      t02 = t02*rd0
      t03 = t03*rd0
      a(j+3,i) = t00
      a(j+3,i+1) = t01
      a(j+3,i+2) = t02
      a(j+3,i+3) = t03
    end do
  end do
  if (nr==0) return

!     here 0 < NR < 4
!     scale last NR cols of A

  call dtrsm('L','U','T','N',n-nr,nr,unity,a,lda,a(0,n-nr),lda)

!     update triangle A(N-NR:N-1,N-NR,N-1) with A(0:N-NR-1,N-NR:N-1)

  call dsyrk('U','T',nr,n-nr,-unity,a(0,n-nr),lda,unity,a(n-nr,n-nr),lda)

!     Cholesky ma54_factor triangle A(N-NR:N-1,N-NR,N-1)

  j = n - nr
  if (nr==3) then ! ma54_factor row J and update A(J+1:J+2,J+1:J+2)
    t00 = a(j,j)
    if (t00<=zero) go to 2
    t00 = sqrt(t00)
    a(j,j) = t00
    rd0 = unity/t00
    ai0 = a(j,j+1)
    ai1 = a(j,j+2)
    ai0 = ai0*rd0
    ai1 = ai1*rd0
    a(j,j+1) = ai0
    a(j,j+2) = ai1
    t00 = a(j+1,j+1)
    t01 = a(j+1,j+2)
    t02 = a(j+2,j+2)
    ai2 = ai0*ai0
    ai3 = ai0*ai1
    ai1 = ai1*ai1
    t00 = t00 - ai2
    t01 = t01 - ai3
    t02 = t02 - ai1
    a(j+1,j+1) = t00
    a(j+1,j+2) = t01
    a(j+2,j+2) = t02
    j = j + 1
  end if
  if (nr>=2) then ! ma54_factor row J and update A(J+1,J+1)
    t00 = a(j,j)
    if (t00<=zero) go to 2
    t00 = sqrt(t00)
    a(j,j) = t00
    rd0 = unity/t00
    ai0 = a(j,j+1)
    ai0 = ai0*rd0
    a(j,j+1) = ai0
    t00 = a(j+1,j+1)
    ai1 = ai0*ai0
    t00 = t00 - ai1
    a(j+1,j+1) = t00
    j = j + 1
  end if
  if (nr>=1) then ! A(J,J) = SQRT( A(J,J) )
    t00 = a(j,j)
    if (t00<=zero) go to 2
    t00 = sqrt(t00)
    a(j,j) = t00
  end if

  return

2 info = j + 1
  return
3 info = j + 2
  return
4 info = j + 3
  return
5 info = j + 4
  return
end subroutine dkcf

end module hsl_ma54_double
