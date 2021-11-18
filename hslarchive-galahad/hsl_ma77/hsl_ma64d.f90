! COPYRIGHT (c) 2007 Council for the Central Laboratory
!               of the Research Councils
!
! Version 6.3.0
! For version history see ChangeLog
!
! To change precision:
!    Change dgemm, dcopy, dswap, dgemv, dtpsv, daxpy
!    Change _double, kind(1.0d0)

module hsl_ma64_double

  implicit none
  EXTERNAL dgemm, dcopy, dswap, dgemv, dtpsv, daxpy
  private
  public ma64_factor, ma64_solveL1, ma64_solveD1, &
         ma64_solveDLT1, ma64_solveLT1, &
         ma64_solveL2, ma64_solveD2, ma64_solveDLT2, &
         ma64_solveLT2, ma64_control, ma64_info
  integer, parameter :: wp = kind(1.0d0) ! Precision parameter.
  integer,parameter:: long = selected_int_kind(18) ! Long integer.
  real (wp), parameter :: unity = 1.0_wp

  interface ma64_factor
      module procedure ma64_factor_double
  end interface

  interface ma64_solveL1
     module procedure ma64_solveL1_double
  end interface

  interface ma64_solveD1
     module procedure ma64_solveD1_double
  end interface

  interface ma64_solveDLT1
     module procedure ma64_solveDLT1_double
  end interface

  interface ma64_solveLT1
     module procedure ma64_solveLT1_double
  end interface

  interface ma64_solveL2
     module procedure ma64_solveL2_double
  end interface

  interface ma64_solveD2
     module procedure ma64_solveD2_double
  end interface

  interface ma64_solveDLT2
     module procedure ma64_solveDLT2_double
  end interface

  interface ma64_solveLT2
     module procedure ma64_solveLT2_double
  end interface

  type ma64_control
     integer :: p_thresh=32 ! If p<=p_thresh, execute on a single thread.
     real (wp) :: small=1e-20_wp ! Diagonal entries of D of modulus less
!                      than this are replaced by zero.
     real (wp) :: static=0.0_wp ! If static>0 and p stable pivots are not found,
!                      the 1x1 pivot that is nearest to satisfying the test is
!                      chosen and info%num_nothresh is incremented by one.
!                      If its absolute value is less than cntl%static, it is
!                      replaced by the nearer of -static and static and
!                      info%num_perturbed is incremented by one.
     logical :: twos=.false. ! If true, the signs of perm indicate recommended
!                      2x2 pivots.
     real (wp) :: u=0.1_wp ! Initial relative pivot tolerance. Values greater
!                      than 0.5 are treated as 0.5 and values less than zero are
!                      treated as zero.
     real (wp) :: umin=1.0_wp ! Minimum relative pivot tolerance. Values greater
!                      than u are treated as u and values less than zero are
!                      treated as zero. If p stable pivots have not been not
!                      found and the candidate pivot with greatest relative
!                      pivot tolerance has tolerance v>=umin, this is accepted
!                      as a pivot and u is reset to v. If p=n> and both u and
!                      umin are greater than 0.5, umin is treated as having the
!                      value 0.5.
end type ma64_control

  type ma64_info
     real (wp) :: detlog=0 ! log of the absolute value of the determinant
!                          of D or zero if the determinant is zero.
     integer :: detsign=0 !  the sign of the determinant of D or zero
!                         if the determinant of D is zero.
!                         Not used in the complex symmetric case.

     integer :: flag=0 ! Zero after successful entry, negative after error.
     integer :: num_neg=0 ! Number of negative eigenvalues of D.
!                           Not used in the complex symmetric case.
     integer :: num_nothresh=0 !  the number diagonal entries of D that
!                      were chosen as 1x1 pivots without satisfying the
!                      relative pivot threshold.
     integer :: num_perturbed=0 ! Number diagonal entries of D that
!                      were perturbed to cntl%static or -cntl%static.
     integer :: num_zero=0 ! Number of zero eigenvalues of D.
     integer :: num_2x2=0 ! Number of 2x2 blocks in D.
     real (wp) :: usmall ! Set to -1 if num_perturbed > 0. Otherwise,
!              if q = p, it holds the smallest relative pivot value of the
!                  chosen pivots, or
!              if q < p and a positive value of cntl%umin would have led to a
!                  greater value of q, it holds the largest such value;
!                  otherwise, it holds zero.
     real (wp) :: u=0 ! Set to the final value of the relative pivot
!                 threshold u.
  end type ma64_info


contains

subroutine ma64_factor_double &
             (n,p,nb,nbi,a,la,cntl,q,ll,perm,d,buf,info,s,n_threads)
! Partial factorization of a symmetric indefinite matrix held in lower
! packed format to the form
!                 P A P' = ( L   ) ( D   ) (L' M')
!                          ( M  I) (   S ) (   I )
! where P is a permutation matrix that involves only the first p rows and
! columns, L is unit lower triangular of size q and D is block diagonal of
! order q with blocks of size 1 and 2.

!$     use omp_lib

  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: p ! Pivots to be chosen in rows and cols 1:p.
  integer, intent (in) :: nb ! Block size
  integer, intent (in) :: nbi ! Inner block size
  integer(long), intent (in) :: la ! Size of a. la>= min (n*n,(n*(n+nb+1))/2).
  real (wp), intent (inout) :: a(la) ! On entry, a(la+1-ln:la) holds the
!              matrix in lower packed format, where ln=(n*(n+1))/2. On exit,
!              a(1:ll) holds ( L ) as a sequence of block columns, each of which
!                            ( M )
!              consists of a triangular matrix packed by columns followed by a
!              rectangular matrix packed by columns. There are nb columns in
!              each block column except the last, which may have fewer. Also,
!              a(la+1-ls:la) holds the matrix S in lower packed format,
!              where ls=((n-q)*(n-q+1))/2.
  type(ma64_control), intent (in) :: cntl ! Controls the actions.
  integer, intent (out) :: q ! Order of D.
  integer(long), intent (out) :: ll ! Size of (L M)' in block-column form,
!              that is, (q*(n+n-q+1))/2.
  integer, intent (inout) :: perm(p) ! The input value is ignored if
!              cntl%twos is false; otherwise, each sequence
!              perm(i)<0, perm(i+1)<0, ... perm(i+k)<0 must be of even length
!              and flags recommended 2x2 pivots.
!              On return, for i = 1, 2, ..., p, perm(i) is set to the index
!              of the row of A that was permuted to row i.
  real (wp), intent (out) :: d(2*p) ! d(1:2*q) is set to hold the inverse of
!              D, except that zero diagonal blocks (pivots) are not inverted.
!              Diagonal entries are in d(1:2*q-1:2) and entries below
!              the diagonal are in d(2:2*q-2:2). d(2*q) is set to zero.
  real (wp) :: buf(nb*n+n) ! Work array for the current part of -LD.
  type(ma64_info), intent (inout) :: info
  integer, intent (in), optional :: s ! Columns 1:s should be searched last for
!              pivots.
  integer, intent (inout), optional, target :: n_threads ! Controls the number
!              of threads used in OpenMP parallel regions. Its value may be
!              changed by another thread during the execution of this
!              subroutine.
!              0  Execute on omp_get_num_threads() threads.
!             >1  Execute on n_threads threads.
!             Absence is treated as equivalent to presence with the value 0.

!     .. Local Scalars ..
    real (wp) amax  ! Entry of largest absolute value to left of diagonal in
!                     row m.
    real (wp) amax2 ! Second largest absolute value to left of diagonal in row m
    real (wp) amaxb ! Largest absolute value below diagonal in column m
    real (wp) amaxt ! Largest absolute value in column t.
    integer   deti   ! Determinant held as detr*radix**deti
    real (wp) detpiv ! Determinant of candidate pivot is detpiv/detscale.
    real (wp) detpiv0 ! First term in calculation of detpiv.
    real (wp) detpiv1 ! Second term in calculation of detpiv.
    real (wp) detpiv2 ! Value of detpiv for best candidate 2x2 pivot.
    real (wp) detr   ! Scaled determinant
    real (wp) detscale ! Inverse of the largest entry in the candidate pivot.
    real (wp) detscale2 ! Value of detscale for best candidate 2x2 pivot.
    integer i ! Row index
    integer i1 ! First row of block
    integer j ! Column index
    integer j1 ! First column of block
    integer j2 ! Last column of block
    integer js ! Last column available for swapping
    integer(long) k  ! Position in a
    integer(long) kdiag  ! Position in a of diagonal np+1
    integer(long) kkj ! Position of diagonal of column j
    integer(long) kkt ! Position of diagonal of column t
    integer(long) kkr1 ! Position of diagonal of column r+1
    integer(long) kkq !  Start of the pivotal block column
    integer(long) kkq1 ! Position of diagonal of column q+1
    integer(long) kkm ! Position of diagonal of column m
    integer(long) kkn ! Position of diagonal of column np
    integer(long) kkmb ! Position of diagonal of column mbest
    integer(long) kkm1 ! Position of diagonal of column m-1
    integer(long) kq1 ! Position in buf of diagonal of column q+1
    integer lj ! Height of block containing column j
    integer lm ! Height of block containing column m
    integer lq ! Height of the pivot block
    integer m  ! Column searched most recently
    integer m1  ! First column of the block in which m appears
    integer mbest ! Column with best relative pivot value for a 1x1 pivot
    integer mbest1, mbest2 ! Pair with best relative pivot value for a 2x2 pivot
    integer mdummy  ! Loop execution count
    integer mlast ! Last column of the block in which m appears
    integer np ! Columns 1 to np are rearranged to block column form.
!                np = p or np = n.
    integer nt  ! Number of threads
    integer pivsiz ! Size of the chosen pivot, 0 if none chosen, or -1
!                    if column is essentially zero
    integer q1 ! First column of the block containing column q+1
    integer q2 ! Last column of the inner block containing column q+1
    integer r ! Number of pivot operations applied to columns q+1:p
    integer rm ! Number of pivot operations applied to the (outer) block
               ! containing column m
    real (wp) rmax ! Largest entry in column m
    real (wp) rmax1 ! Largest entry in column m-1 outwith rows m,m-1 if (m,m-1)
                    ! is a recommended 2x2 pivot. Otherwise, -1.
    real (wp) rmax2 ! Largest entry in column m outwith rows m,m-1.
    integer t ! Candidate 2x2 pivot is in columns t and m
    real (wp) ubest1 ! Relative pivot value of best candidate 1x1 pivot
    real (wp) ubest2 ! Relative pivot value of best candidate 2x2 pivot
    real (wp) urel ! Relative pivot value
    real (wp) urel1 ! Relative pivot value in column m-1 if (m,m-1) is a
!                    recommended 2x2 pivot.
    real (wp) u ! Relative pivot threshold
    real (wp) umin ! Minimum relative pivot threshold

    info%flag = 0
    if (n < 0) then
       info%flag = -1
    else if (p < 0) then
       info%flag = -2
    else if (p > n) then
       info%flag = -3
    else if (nbi <= 1) then
       info%flag = -4
    else if ( la < min(n*int(n,long),(n*(n+nb+1_long))/2) ) then
       info%flag = -7
    else if (cntl%static < cntl%small .and. cntl%static/=0.0_wp) then
       info%flag = -10
    else if (mod(nb,nbi) /= 0) then
       info%flag = -12
    end if
    info%detlog = 0.0_wp
    info%num_neg = 0
    info%num_nothresh = 0
    info%num_perturbed = 0
    info%num_zero = 0
    info%num_2x2 = 0
    u = min(max(cntl%u,0.0_wp),0.5_wp)
    info%u = u
    umin = min(max(cntl%umin,0.0_wp),u)
    if (p==n) umin = min(umin,0.5_wp)
    info%usmall = 1.0_wp
    q = 0
    if (info%flag/=0 .or. p==0) return
    deti = 0
    detr = 1.0_wp
    nt = 1
!$  if (p>cntl%p_thresh) nt = omp_get_max_threads()
    np = n
    if( nt==1 .and. p<=nb ) np = p
    lq = n
    kkq1 = 1
    kkq = 1
    m = p
! m is updated at the start of the main loop so initializing it to p causes
! it to have the value 1 during the first execution of that loop.
    js = p
    q1 = 1
    q2 = min(nbi,p)
    rmax1 = -1.0_wp

! For j = 1,p: set perm(j) = -j for the first part of a 2x2 pivot; otherwise,
! set perm(j) = j
    if (cntl%twos) then
      js = 0
      do j = 1, p
        if (perm(j) >= 0) then
          perm(j) = j
        else
! First part of a 2x2 pivot
          perm(j) = -j
          if (j==p) then
             info%flag = -11
             return
          else if (perm(j+1) >= 0) then
             info%flag = -11
             return
          else
             perm(j+1) = j+1
          end if
        end if
      end do
    else
      do j = 1, p
        perm(j) = j
      end do
    end if

! Rearrange first np columns to block column form
    call to_block

  if (present(s)) then
     if (s>0 .and. s<p) then
        js = js - s
        i = min(s,p-s)
! Ignore recommended 2x2 pivots that will be split
        if (perm(i)<0) perm(i) = i
        if (perm(p-i)<0) perm(p-i) = p-i
! Make first i columns be last
        do j = 1,i
           call swap_cols(j,p+1-j)
        end do
! Ensure that the first index of a 2x2 pair is labelled
        do j = 2,i
           if (perm(j)<0) then
               perm(j) = -perm(j)
               perm(j-1) = -perm(j-1)
           end if
        end do
        do j = p+2-i,p
           if (perm(j)<0) then
               perm(j) = -perm(j)
               perm(j-1) = -perm(j-1)
           end if
        end do
     end if
  end if

  pivot: do
! Perform a pivotal operation
    pivsiz = 0
    ubest1 = 0.0_wp
    ubest2 = 0.0_wp
    sweep: do mdummy = 1, p-q ! Look for a pivotal column or pair of columns

! Update m and the scalars associated with column m
      m = m+1
      if (m>p) then
! Go back to column q+1
         m = q+1
         m1 = q1
         kkm = kkq1
         lm = lq
         r = q
         rm = q
      else if (m<m1+nb) then
! Within the current block column
         kkm = kkm + lm + 1
      else
! At the the start of another block column
         lm = lm - nb
         m1 = m1 + nb
         rm = r
         kkm = kkm + lm + 1
      end if

! If there are more than nbi active columns and there are enough columns
! not yet considered, update all candidate columns and perform swaps.
      if (m-q>nbi) then
        if (js-m > m-q) then
           mlast = min(p,m1+nb-1)
!     Inner update: of columns (m:mlast) for pivots rm+1:q (BLAS3)
           if (m<=mlast .and. rm<q) call update (m,mlast,rm+1,q)
!     Outer update: of columns (m1+nb:p) for pivots r+1:q (BLAS3)
           if (m1+nb<=p .and. r<q) call update (m1+nb,p,r+1,q)
           r = q
           rm = r
           do j = q+1,m
              call swap_cols(j,js)
              js = js -1
           end do
        end if
      end if

! Update column m
      kkr1 = kkq1 + (rm-q)*(lq+1_long)
      k = kkr1+m-rm-1
      if(q>rm) then
         call dgemv('NoTrans',n-m+1,q-rm,unity,a(k),lq, &
                    buf(n*(rm+1_long-q1)+m),n,unity,a(kkm),1)
      end if


! Find largest and second largest entry to left of diagonal in row m.
       j = q + 1
       lj = lq
       k = kkq1 + m - j - lj
       amax = 0.0_wp
       amax2 = 0.0_wp
       t = 0
       if (j<m) then
          t = j
          kkt = kkq1
       end if
       do j1 = q1,m-1,nb
          do j = j, min(j1+nb-1,m-1)
             k = k + lj
             if (abs(a(k))>abs(amax)) then
                t = j
                amax2 = abs(amax)
                amax = a(k)
                kkt = k - (m-j)
             else
                amax2 = max(abs(a(k)),amax2)
             end if
          end do
          lj = lj - nb
       end do

! Now calculate largest entry below the diagonal of column m.
       amaxb = 0.0_wp
       do i = m+1,n
          amaxb = max(abs(a(kkm+i-m)),amaxb)
       end do

! Now calculate largest entry in the whole of column m and make sure that is
! is neither small nor infinity.
       rmax = max(abs(amax),abs(a(kkm)),amaxb)
       if (rmax<=cntl%small) then
! All entries of the column are small
           a(kkm) = 0.0_wp
           info%num_zero = info%num_zero + 1
           pivsiz = -1
           perm(m) = abs(perm(m))
           rmax1 = -1.0_wp
           exit sweep
       else if (rmax > huge(1.0_wp)) then
! There is an infinity in the column
           info%flag = -13
           return
       end if

! Calculate the relative pivot value and see if it is the best so far
       if (abs(a(kkm))>cntl%small) then
          urel = abs(a(kkm))/rmax
       else
          urel = 0.0_wp
       end if
       if (urel >= ubest1) then
          ubest1 = urel
          mbest = m
          kkmb = kkm
       end if

! If (m,m+1) recommended as 2x2 pivot store needed data and cycle so that
! column m+1 is updated
       if (perm(m)<0) then
          rmax1 = abs(amax)
          do i = m+2,n
             rmax1 = max(abs(a(kkm+i-m)),rmax1)
          end do
          urel1 = urel
          kkm1 = kkm
          perm(m) = abs(perm(m))
          cycle sweep
       end if

! If (m-1,m) recommended as 2x2 pivot, look for pivot in columns m-1,m
       lookm1m: if (rmax1>=0.0_wp) then

! Find largest entry in column m outwith rows m-1 and m
          rmax2 = amaxb
          j = q + 1
          lj = lq
          k = kkq1 + m - j - lj
          do j1 = q1,m-2,nb
             do j = j, min(j1+nb-1,m-2)
                k = k + lj
                rmax2 = max(abs(a(k)),rmax2)
             end do
             lj = lj - nb
          end do
          k = kkm1 + 1

! If both diagonal entries >= cntl%small or the off-diagonal entry
! >= cntl%small, consider 2x2 pivot
          test2x2: if ( min(abs(a(kkm)),abs(a(kkm1))).ge.cntl%small .or. &
              abs(a(k)).ge.cntl%small ) then
             detscale = 1.0_wp/max( abs(a(kkm)), abs(a(kkm1)), abs(a(k)) )
             detpiv1 =  (abs(a(k))*detscale)*abs(a(k))
             detpiv0 = (a(kkm)*detscale)*a(kkm1)
             detpiv = detpiv0 - detpiv1
! Make sure that the 2x2 pivot is not singular and that there is little
! cancellation in calculating detpiv. Bearing in mind the way detscale
! is calculated, if the largest entry of the matrix is chosen as pivot,
! the one entry of the reduced matrix has absolute value abs(detpiv).
             if (abs(detpiv)> &
                 max(cntl%small,abs(detpiv0)/2,abs(detpiv1)/2)) then
! OK as 2x2 pivot if all entries in the rest of the columns are small
                if(max(rmax2,rmax1)<=cntl%small) then
                   pivsiz = 2
                   rmax1 = -1.0_wp
                   t = m-1
                   exit sweep
                end if

! Calculate the relative pivot value (as 2x2 pivot)
                urel = abs(detpiv)/max(                         &
                    abs(a(kkm)*detscale)*rmax1+abs(a(k)*detscale)*rmax2,   &
                    abs(a(kkm1)*detscale)*rmax2+abs(a(k)*detscale)*(rmax1) )
                rmax1 = -1.0_wp

! OK as 2x2 pivot if relative pivot value is big enough
                if (urel>u) then
                   pivsiz = 2
                   t = m-1
                   info%usmall = min(urel,info%usmall)
                   exit sweep
                end if

! If this has the best relative pivot value so far, record this
                if(urel>ubest2)then
                   ubest2 = urel
                   detpiv2 = detpiv
                   detscale2 = detscale
                   mbest1 = m
                   mbest2 = m-1
                end if

             end if

          end if test2x2

! If m-1 OK as 1x1 pivot, accept this
          rmax1 = -1.0_wp
          if (urel1>u) then
             if (abs(a(kkm1))>cntl%small) then
                pivsiz = 1
                info%usmall = min(urel1,info%usmall)
                call swap_cols(m-1,m)
                exit sweep
             end if
          end if

       end if lookm1m

! If there is a candidate 2x2 pivot, try it.
! Look for a 1x1 pivot only if the 2x2 pivot is unacceptable.
       tgt0: if (t>0) then
       if ( min(abs(a(kkm)),abs(a(kkt))).ge.cntl%small .or. &
              abs(amax).ge.cntl%small ) then
! Store value of the largest entry in whole of column m outwith rows m and t
          rmax2 = max(amax2,amaxb)

! Find largest entry in column t outwith rows m and t
          amaxt = 0.0_wp
          j = q + 1
          lj = lq
          k = kkq1 + t - j - lj
          do j1 = q1,t-1,nb
            do j = j, min(j1+nb-1,t-1)
               k = k + lj
               amaxt = max(abs(a(k)),amaxt)
            end do
            lj = lj - nb
          end do
          if (j1/=t) lj = lj + nb
          k = k + lj
          do i = t+1,m-1
            amaxt = max(abs(a(k+i-t)),amaxt)
          end do
          do i = m+1,n
            amaxt = max(abs(a(k+i-t)),amaxt)
          end do

           detscale = 1.0_wp/max( abs(a(kkm)), abs(a(kkt)), abs(amax) )
          detpiv1 =  (amax*detscale)*amax
          detpiv0 = a(kkm)*detscale*a(kkt)
          detpiv = detpiv0 - detpiv1
! Make sure that the 2x2 pivot is not singular and that there is little
! cancellation in calculating detpiv. Bearing in mind the way detscale
! is calculated, if the largest entry of the matrix is chosen as pivot,
! the one entry of the reduced matrix has absolute value abs(detpiv).
          left2x2:if (abs(detpiv)> &
              max(cntl%small,abs(detpiv0)/2,abs(detpiv1)/2)) then

! OK as 2x2 pivot if all entries in the rest of the columns are small
            if(max(rmax2,amaxt)<=cntl%small) then
                pivsiz = 2
                exit sweep
            end if

! Calculate the relative pivot value (as 2x2 pivot)
            urel = abs(detpiv)/max( &
               abs(a(kkm)*detscale)*amaxt+(abs(amax)*detscale)*rmax2, &
               abs(a(kkt)*detscale)*rmax2+(abs(amax)*detscale)*amaxt )

! OK as 2x2 pivot if relative pivot value is big enough
            if (urel>u) then
              pivsiz = 2
               info%usmall = min(urel,info%usmall)
               exit sweep
            end if

! If this has the best relative pivot value so far, record this
            if(urel>ubest2)then
               ubest2 = urel
               detpiv2 = detpiv
               detscale2 = detscale
               mbest1 = m
               mbest2 = t
            end if
         end if left2x2
      end if
      end if tgt0

! If 2x2 pivot rejected or only one column left, take best 1x1 pivot if is OK.
      if (t>0 .or. m==p) then
         if (ubest1>u) then
            pivsiz = 1
            info%usmall = min(ubest1,info%usmall)
            if (mbest/=m) call swap_cols(mbest,m)
            exit sweep
         end if
      end if

    end do sweep

! No pivot found in search of all available columns
    pivsiz0: if (pivsiz==0) then
! Since all the columns have been updated, revise m to q+1
        m = q+1
        kkm = kkq1
        m1 = q1
        lm = lq
        rm = q
        r = q
! Perform relaxed pivoting if the best pivot is good enough
        if (max(ubest1,ubest2)>=umin) then
           if (ubest1>=ubest2) then
! Accept 1x1 pivot
              u = min(u,ubest1)
              info%usmall = min(ubest1,info%usmall)
              pivsiz = 1
              if (mbest/=m) call swap_cols(m,mbest)
           else
! Accept 2x2 pivot
              u = min(u,ubest2)
              info%usmall = min(ubest2,info%usmall)
              pivsiz = 2
! Revise m to q+1
              m = q+2
              kkm = kkq1+lq+1
              if (m==m1+nb) then
                 m1 = m1+nb
                 lm = lm - nb
              end if
              detpiv = detpiv2
              detscale = detscale2
              t = min(mbest1,mbest2)
              if (t/=q+1) call swap_cols(q+1,t)
              t = max(mbest1,mbest2)
              if (t/=q+2) call swap_cols(q+2,t)
              t = q+1
           end if

! Perform static pivoting if this has been requested
        else if (cntl%static>0) then
           if (mbest/=m) call swap_cols(m,mbest)
           pivsiz = 1
           info%num_nothresh = info%num_nothresh + 1
           if (abs(a(kkm)) < cntl%static) then
               a(kkm) = sign(cntl%static,real(a(kkm),wp))
               info%num_perturbed = info%num_perturbed + 1
               info%usmall = -1
           else
               info%usmall = min(ubest1,info%usmall)
           end if
        end if
    end if pivsiz0

    pivsizes: if (pivsiz==1) then
! Swap columns if m not q+1
      if (q+1/=m) call swap_cols(q+1,m)
! We store D**-1.
      d(2*q+1) = 1.0_wp/a(kkq1)
      d(2*q+2) = 0.0_wp
! Update info
!      info%detlog = info%detlog + log(abs(a(kkq1)))
      deti = deti + exponent(detr) + exponent(abs(a(kkq1)))
      detr = fraction(detr)*fraction(abs(a(kkq1)))
      if (a(kkq1)<0.0_wp) info%num_neg = info%num_neg + 1
! Store L in A and -LD or -conjg(LD) in buf
      kq1 = q+1+(q+1_long-q1)*n
      buf(kq1) = -a(kkq1)
      a(kkq1) = 1.0_wp
      do i = 1, n-q-1
         buf(kq1+i) = -a(kkq1+i)
         a(kkq1+i) = d(2*q+1)*a(kkq1+i)
      end do
      j1 = q1
      lj = lq
! Update columns q+2 to m
      j = q+2
      kkj = kkq1
      do j1 = j1,m,nb
         do j = j, min(j1+nb-1,m)
            kkj = kkj + lj + 1
            call daxpy(n-j+1,buf(kq1+j-q-1),a(kkq1+j-q-1),1,a(kkj),1)
         end do
         lj = lj - nb
      end do
! Update q and kkq1
      kkq1 = kkq1 + lq + 1
      if(q+2==q1+nb)kkq1 = kkq1 - nb
      q = q+1

    else if (pivsiz==2) then pivsizes
! Swap columns unless t==q+1 and m==q+2
      if (q+2/=m) then
         if (q+1/=t) call swap_cols(q+1,t)
         call swap_cols(q+2,m)
      end if
! We store D**-1
      k = kkq1 + lq + 1
      if(q+2==q1+nb)k = k - nb
      d(2*q+1) = (a(k)*detscale)/detpiv
      d(2*q+3) = (a(kkq1)*detscale)/detpiv
      d(2*q+2) = -(a(kkq1+1)*detscale)/detpiv
      d(2*q+4) = 0.0_wp

! Update info
      info%num_2x2 = info%num_2x2 + 1
!     info%detlog = info%detlog + log(abs(detpiv)) - log(detscale)
      deti = deti + exponent(detr) - exponent(detscale)
      detr = fraction(detr)*abs(detpiv)/fraction(detscale)
      if (detpiv<0.0_wp) then
         info%num_neg = info%num_neg + 1
      else if (a(kkq1)+a(k)<0.0_wp) then
         info%num_neg = info%num_neg + 2
      end if
! Store L in A and -LD or -conjg(LD) in buf
      kq1 = q+1+(q+1_long-q1)*n
      buf(kq1) = -a(kkq1)
      buf(kq1+1) = -a(kkq1+1)
      buf(kq1+n) = -a(kkq1+1)
      buf(kq1+n+1) = -a(k)
      a(kkq1) = 1.0_wp
      a(k) = 1.0_wp
      a(kkq1+1) = 0.0_wp
      do i = 2, n-q-1
         buf(kq1+i  ) = -a(kkq1+i)
         buf(kq1+i+n) = -a(k+i-1)
         a(kkq1+i) = d(2*q+1)*a(kkq1+i) + d(2*q+2)*a(k+i-1)
         a(k+i-1) = d(2*q+3)*a(k+i-1) - d(2*q+2)*buf(kq1+i)
      end do
! Update columns q+3 to m
      j1 = q1
      lj = lq
      if(q+3>=q1+nb) then
! Pivot is at end of block
         lj = lq - nb
         j1 = q1 + nb
      end if
      kkj = k + lj + 1
      do j1 = j1,m,nb
         j = max(q+3,j1)
         j2 = min(j1+nb-1,m)
         if(j2>=j) then
            call dgemm('n','t',n-j+1,j2-j+1,2,unity, &
                   a(kkq1+j-q-1),n-((q+1)/nb)*nb,buf(kq1+j-q-1),n, &
                   unity,a(kkj),lj)
         end if
         kkj = kkj + (lj+1_long)*(j2-j+1) - nb
         lj = lj - nb
      end do
! Update q and kkq1
      kkq1 = k + lq + 1
      if (q+3>=q1+nb) kkq1 = kkq1 - nb ! Pivot is at end of block
      q = q + 2

    else if (pivsiz==-1) then pivsizes
! Handle a row that is zero
      if (q+1/=m) call swap_cols(q+1,m)
      d(2*q+1) = 0.0_wp
      d(2*q+2) = 0.0_wp
! Store L in A and -LD in buf
      kq1 = q+1+(q+1_long-q1)*n
      buf(kq1) = 0.0_wp
      a(kkq1) = 1.0_wp
      do i = 1, n-q-1
         buf(kq1+i) = 0.0_wp
         a(kkq1+i) = 0.0_wp
      end do
      kkq1 = kkq1 + lq + 1
      if(q+2==q1+nb)kkq1 = kkq1 - nb
      q = q+1
    end if pivsizes

    j = m1+nb
    if (q>=q2) then
!     Inner block is complete
!     See if the inner update can be done with the outer update
      if( q>=q1+nb-1 .and. r==rm ) then
         j = m+1
      else
!     Update columns (m+1:min(p,m1+nb-1)) for pivots rm+1:q2 (BLAS3)
        mlast = min(p,m1+nb-1)
        if (m<mlast) call update (m+1,mlast,rm+1,q2)
      end if
      rm = q2
      q2 = min (p,q2+nbi)
    end if

    if (q>=q1+nb-1) then
!     Outer block is complete
!     Update columns (j:p) for pivots r+1:q1+nb-1 (BLAS3)
      if (j<=p) call update (j,p,r+1,q1+nb-1)
      r = q1+nb-1
      rm = r

      if (p<n) then
!     Update columns (p+1:n) for pivots q1:q1+nb-1
         if(np<n) then
            call update2 (q1+nb-1)
         else
            call update (p+1,n,q1,q1+nb-1)
         end if
      end if

      q1 = q1 + nb
      if (q1==q) then
!     If a 2x2 pivot spans the blocks, move final column of buf forward
         call dcopy(n-q+1,buf(q+n*int(nb,long)),1,buf(q),1)
      end if
      kkq = kkq + int(nb,long)*lq
      lq = lq - nb

    end if

    if (q == p .or. pivsiz==0) exit pivot
  end do pivot

  if (q>=q1 .and. p<n) then
!  Update columns (p+1:n) for pivots q1:q (BLAS3)
     if(np<n) then
        call update2 (q)
     else
        call update (p+1,n,q1,q)
     end if
  end if

! Update info
  if(q<p .and. info%num_nothresh==0) &
            info%usmall = min(info%usmall, max(ubest1,ubest2))
  if(info%num_zero/=0) then
     info%detsign = 0
     info%detlog = 0.0_wp
  else
     info%detsign = -1
     if (mod(info%num_neg,2)==0)info%detsign = 1
     info%detlog = log(detr)+deti*log(1.0_wp*radix(1.0_wp))
  end if
  info%u = u

! Rearrage columns 1 to np
  call from_block




contains

    subroutine update(jl,jr,pl,pr)
!     Update columns jl:jr for pivots pl:pr when jr<=np.
       integer jl,jr,pl,pr
! Local variables
       integer ij ! Index for looping over blocks
       integer i2 ! Last row of block
       integer j  ! Column index
       integer nc  ! Number of block columns to be updated
       integer nr  ! Number of block rows to be updated
       integer nn  ! Number of blocks to be updated
       integer j0  ! Last column of block ahead of jl
       integer p0  ! First column of pivot block

       j0 = nb*((jl-1)/nb)
       p0 = nb*((pl-1)/nb)+1
       nr = 1 + (n-j0-1)/nb
       nc = 1 + (jr-j0-1)/nb
       nn = (nc*(2*nr-nc+1))/2

! Decide on the number of threads to use
       nt = 1
       if (nn>1 .and. p>cntl%p_thresh .and. np==n) then
!$        nt = omp_get_max_threads()
          if( nt>1 ) then
             if ( present(n_threads) ) then
                if( n_threads>0 ) nt = n_threads
             end if
          end if
       end if

       if(nt==1) then  ! No parallelization
          do j1 = j0+1, jr, nb
             j = max(jl,j1)
             j2 = min(j1+nb-1,jr)
             if(n.ge.j .and. j2.ge.j) then
                kkj = j1-1
                kkj = 1 + ((kkj)*(n+n-kkj+nb))/2
                if(j.ne.j1) kkj = kkj + (j-j1)*(n-j0+1)
                call dgemm &
                    ('n','t',n-j+1,j2-j+1,pr-pl+1,unity, &
                    a(kkq+j-q1+lq*(pl-p0)),lq,buf(j+n*(pl-p0)), &
                    n,unity,a(kkj),n-j1+1)
            endif
          end do

       else ! Parallelize by blocks
!$        if (nt>0) call omp_set_num_threads(nt)
!$OMP     PARALLEL DO &
!$OMP     DEFAULT(SHARED) PRIVATE(I, J, I1, J1, I2, J2, KKJ) &
!$OMP     SCHEDULE(STATIC)
          do ij = 1, nn
! Find the block indices
             j = int( nr + 1.4999d0 - sqrt( (nr+0.5d0)**2 -ij*2 ) )
             i = ij - ((j-1)*(nr+nr-j))/2
! Find the leading and trailing row and column
             i1 = j0+1+(i-1)*nb
             j1 = j0+1+(j-1)*nb
             i2 = min(i1+nb-1,n)
             j2 = min(j1+nb-1,jr)
! Find the first row and column to be updated
             j = max(jl,j1)
             i = max(jl,i1)
! Find position of diagonal of column j
             kkj = j1-1
             kkj = 1 + ((kkj)*(n+n-kkj+nb))/2
             if(j/=j1) kkj = kkj + (j-j1)*(n-j0+1)
             call dgemm &
                 ('n','t',i2-i+1,j2-j+1,pr-pl+1,unity,a(kkq+i-q1+lq*(pl-p0)) &
                        ,lq,buf(j+n*(pl-p0)),n,unity,a(kkj+i-j),n-j1+1)
          end do
!OMP      END PARALLEL DO
       end if
    end subroutine update

    subroutine update2(pr)
!     Update columns p+1:n for pivots q1:pr within the current block column
!     when np<n. Columns p+1:n are held in packed format.
      integer pr
! Local variables
       integer j1  ! First column of block column
       integer j2  ! Last column of block column
       integer(long) k ! Column start in packed form
       integer(long) kfree ! Start of free space in a
       integer(long) l ! Column start in unpacked form
       integer p0  ! First column of pivot block

       p0 = nb*((q1-1)/nb)+1
       kkj = kdiag
       kfree = kkn + n+1-np
       do j1 = p+1, n, nb/2
          j2 = min(j1+nb/2-1,n)
! Copy block of columns to rectangular format
          k = kkj
          l = kfree
          do j = j1, j2
             call dcopy(n+1-j,a(k),1,a(l),1)
             k = k + n+1-j
             l = l + n+2-j1
          end do
! Perform the update
          call dgemm &
              ('n','t',n+1-j1,j-j1,pr-q1+1,unity,a(kkq+j1-q1+lq*(q1-p0)), &
                     lq,buf(j1+n*(q1-p0)),n,unity,a(kfree),n+1-j1)
! Copy block of columns back
          k = kkj
          l = kfree
          do j = j1, j2
             call dcopy(n+1-j,a(l),1,a(k),1)
             k = k + n+1-j
             l = l + n+2-j1
          end do
          kkj = k
       end do
    end subroutine update2

    subroutine to_block
! Rearrange the first np columns of a packed triangular matrix to block packed
! triangular format
      integer j ! Column index
      integer j1 ! First column of block
      integer(long)  k ! Position in packed triangular format
      integer(long)  l ! Position in block packed triangular format
      integer lj ! Block length
      lj = n
      k = la + 1 - (n*(n+1_long))/2
      l = -lj
      do j1 =  1, np, nb
         do j = j1, min(j1+nb-1,np)
! Move column
            l = l + lj + 1
            a(l-j+j1:l-1) = 0.0_wp
            call dcopy(n+1-j,a(k),1,a(l),1)
            k = k + n+1-j
         end do
         lj = lj - nb
      end do
      kkn = l
      kdiag = k
      if (np<n) then
! To avoid traps on undefined values in calls of gemm from update2, set zeros
! in the part of a that is used as the first buffer and corresponds to
! unwanted values that are nevertheless passed to gemm.
        k = kkn + n+1-np
        do j = 1,min(n-p,nb/2)
          a(k:k+j-2) = 0.0_wp
          k = k + n-p
        end do
      end if
   end subroutine to_block

   subroutine from_block
! Rearrange columns 1 to q of a block packed triangular matrix so that each
! block column consists of the diagonal block packed by columns followed by
! the off-diagonal part as a rectangular matrix held by columns, and rearrange
! columns q+1 to np to packed triangular format.
      integer j ! Column index
      integer jb ! Number of columns in the block
      integer j1 ! First column of block
      integer(long) k ! Position in packed triangular format
      integer(long) k1 ! Position of the start of the block column
      integer(long) l ! Position in block packed triangular format
      integer(long) l1 ! New position in block packed triangular format
      integer lj ! Block length

      lj = n
      k1 =  1
      l1 =  1
      do j1 =  1, q, nb
         jb = min(nb,q-j1+1)
! Copy diagonal block to buf
         k =  1
         l = k1
         do j = 1, jb
            call dcopy(jb+1-j,a(l),1,buf(k),1)
            k = k + jb+1-j
            l = l + lj + 1
         end do
! Copy off-diagonal block to buf
         l = k1 + jb
         do j = 1, jb
            call dcopy(lj-jb,a(l),1,buf(k),1)
            k = k + lj-jb
            l = l + lj
         end do
         k1 = l - jb
! Copy the whole block column back
         k = 1
         do j = 1, jb
            call dcopy(lj+1-j,buf(k),1,a(l1),1)
            k = k+lj+1-j
            l1 = l1+lj+1-j
         end do
         lj = lj - nb
      end do
      ll = l1 - 1

! Rearrange columns q+1 to np to packed triangular format
      m = nb*(np/nb)
      q1 = 1 + nb*((q-1)/nb)
      lj = n - m
      l = kkn
      k = kdiag
      do j1 = m+1, q1, -nb
         do j = min(j1+nb-1,np), max(j1,q+1), -1
! Move column
            k = k - (n+1-j)
            call dcopy(n+1-j,a(l),1,a(k),1)
            l = l - lj - 1
         end do
         lj = lj + nb
       end do
   end subroutine from_block

   subroutine swap_cols(j1,j2)
   integer j1,j2
! Swap columns. It may be assumed that 1<=j1<j2<=p.

   integer(long) d1,d2 ! Positions of the diagonal entries for columns j1,j2
   integer j ! Index of the leading column of the block
   integer jb ! Size of block
   integer(long) k1 ! Start of current part of column j1
   integer(long) k2 ! Start of current part of column j2
   integer l !
   integer lb ! No. rows in block column
   real(wp) temp

   j = perm(j1)
   perm(j1) = perm(j2)
   perm(j2) = j

! Swap columns of buffer
   call dswap(q-q1+1,buf(j1),n,buf(j2),n)

   k1 = j1
   k2 = j2
   lb = n

! Swap rows
   do j = 1, j1, nb
      jb = min(nb,n-j+1)
      l = min(nb,j1-j)
      if (l>0) call dswap(l,a(k1),lb,a(k2),lb)
      k1 = k1 + l*int(lb,long) - jb
      k2 = k2 + l*int(lb,long) - jb
      lb = lb - nb
   end do

   j = j - nb
   lb = lb + nb
   d1 = k1 + jb
   k1 = d1 + 1
   k2 = k2 + lb + jb

! Swap columns with rows
   jb = min(nb,j2-j1-1)
   do j = j, j2, nb
      l = min(jb,j+nb-1-j1,j2-j)
      if (l>0) call dswap(l,a(k1),1,a(k2),lb)
      k1 = k1 + l
      k2 = k2 + l*int(lb,long) - nb
      lb = lb - nb
   end do

   j = j - nb
   lb = lb + nb
   d2 = k2 + nb
   k2 = d2 + 1
   k1 = k1 + 1

! Swap the diagonals
   temp = a(d1)
   a(d1) = a(d2)
   a(d2) = temp

! Swap columns
   if (n>j2) call dswap(n-j2,a(k1),1,a(k2),1)

   end subroutine swap_cols
end subroutine ma64_factor_double



  subroutine ma64_solveL1_double(n,q,nb,b,flag,a,ll)
! Solves                ( L   ) x = b
!                       ( M  I)
! where  L is unit lower triangular of size q.

    integer, intent (in) :: n ! Matrix order
    integer, intent (in) :: q ! Order of L (and D)
    integer, intent (in) :: nb ! The block size used for the blocked hybrid
!                                format.
    real (wp), intent (inout) :: b(n) ! b(n) holds the right-hand
!                           sides on entry and is overwritten by the solution.
    integer, intent (out) :: flag
    integer(long), intent (in) :: ll ! Size of a
    real (wp), intent (in) :: a( ll )
!                               Holds ( L ) packed by columns.
!                                     ( M )

! .. Locals ..
    integer j ! Column index
    integer jb ! Width of block
    integer(long) k ! Position in array
    integer lj ! Length of block

    flag = 0
    if (n < 0) then
       flag = -1
    else if (nb <= 1) then
       flag = -4
    else if (q < 0) then
       flag = -8
    else if (q > n) then
       flag = -9
    end if
    if (flag/=0 .or. q==0) return
    k = 1
    lj = n
    do j = 1, q, nb
      jb = min(nb,q-j+1)
      call dtpsv('L','N','U',jb,a(k),b(j),1)
      k = k + (jb*(jb+1_long))/2
      lj = lj - jb
      if (lj>0) call dgemv('N',lj,jb,-unity,a(k),lj,b(j),1,unity,b(j+jb),1)
      k = k + jb*lj
    end do

  end subroutine ma64_solveL1_double

  subroutine ma64_solveD1_double(n,q,b,flag,d)
! Solves                ( D   ) x = b
!                       (   I )
! where D is block diagonal of size q with blocks of size 1 and 2.

   integer, intent (in) :: n ! Matrix order
    integer, intent (in) :: q ! Order of L (and D)
    real (wp), intent (inout) :: b(n) ! b(n) holds the right-hand
!                           sides on entry and is overwritten by the solution.
    integer, intent (out) :: flag
    real (wp), intent (in) :: d(2*q) !  Holds the inverse of D,
!              except that zero diagonal blocks (pivots) are not inverted.
!              Diagonal entries are in d(1:2*q-1:2) and entries below
!              the diagonal are in d(2:2*q-2:2).

! .. Locals ..
    integer j ! Column index
    real (wp) temp

    flag = 0
    if (n < 0) then
       flag = -1
    else if (q < 0) then
       flag = -8
    else if (q > n) then
       flag = -9
    end if
    if (flag/=0 .or. q==0) return

! Apply operations associated with D
    j = 1
    do
      if (j>q) exit
      if (d(2*j)==0.0_wp) then
! 1x1 pivot
        b(j) = b(j)*d(2*j-1)
        j = j + 1
      else
! 2x2 pivot
        temp = b(j)*d(2*j-1) + b(j+1)*d(2*j)
        b(j+1) = b(j)*d(2*j) + b(j+1)*d(2*j+1)
        b(j) = temp
        j = j + 2
      end if
    end do

  end subroutine ma64_solveD1_double

  subroutine ma64_solveLT1_double(n,q,nb,b,flag,a,ll)
! Solves                ( L' M') x = b
!                       (    I )
! where  L is unit lower triangular of size q.

   integer, intent (in) :: n ! Matrix order
    integer, intent (in) :: q ! Order of L (and D)
    integer, intent (in) :: nb ! The block size used for the blocked hybrid
!                                format.
    real (wp), intent (inout) :: b(n) ! b(n) holds the right-hand
!                           sides on entry and is overwritten by the solution.
    integer, intent (out) :: flag
    integer(long), intent (in) :: ll ! Size of a
    real (wp), intent (in) :: a( ll )
!                               Holds ( L ) packed by columns.
!                                     ( M )

! .. Locals ..
    integer j ! Column index
    integer jb ! Width of block
    integer j1 ! First column of block
    integer(long) k ! Position in array
    integer lj ! Length of block
    integer m !

    flag = 0
! Back substitution
    flag = 0
    if (n < 0) then
       flag = -1
    else if (nb <= 1) then
       flag = -4
    else if (q < 0) then
       flag = -8
    else if (q > n) then
       flag = -9
    end if
    if (flag/=0 .or. q==0) return
    m = nb*((q-1)/nb)
    k = 1 + ll
    j1 = 1 + m
    lj = n - q
    do j = j1, 1, -nb
      jb = min(nb,q-j+1)
      k = k - lj*jb
      if (lj>0) call dgemv('T',lj,jb,-unity,a(k),lj,b(j+jb),1,unity,b(j),1)
      k = k - (jb*(jb+1_long))/2
      call dtpsv('L','T','U',jb,a(k),b(j),1)
      lj = lj + jb
    end do

  end subroutine ma64_solveLT1_double

  subroutine ma64_solveDLT1_double(n,q,nb,b,flag,a,ll,d)
! Solves       ( D   ) ( L' M') x = b
!              (   I ) (    I )
! where D is block diagonal of size q with blocks of size 1 and 2
! and L is unit lower triangular of size q.

    integer, intent (in) :: n ! Matrix order
    integer, intent (in) :: q ! Order of L (and D)
    integer, intent (in) :: nb ! The block size used for the blocked hybrid
!                                format.
    real (wp), intent (inout) :: b(n) ! b(n) holds the right-hand
!                           sides on entry and is overwritten by the solution.
    integer, intent (out) :: flag
    integer(long), intent (in) :: ll ! Size of a
    real (wp), intent (in) :: a( ll )
!                               Holds ( L ) packed by columns.
!                                     ( M )
    real (wp), intent (in) :: d(2*q) !  Holds the inverse of D,
!              except that zero diagonal blocks (pivots) are not inverted.
!              Diagonal entries are in d(1:2*q-1:2) and entries below
!              the diagonal are in d(2:2*q-2:2).

    if (nb <= 1) then
       flag = -4
       return
    end if
    call ma64_solveD1_double(n,q,b,flag,d)
    call ma64_solveLT1_double(n,q,nb,b,flag,a,ll)

  end subroutine ma64_solveDLT1_double

  subroutine ma64_solveL2_double(n,q,nb,nrhs,b,ldb,flag,a,ll)
! Solves                ( L   ) x = b
!                       ( M  I)
! where  L is unit lower triangular of size q.

    integer, intent (in) :: n ! Matrix order
    integer, intent (in) :: q ! Order of L (and D)
    integer, intent (in) :: nb ! The block size used for the blocked hybrid
!                                format.
    integer, intent (in) :: nrhs ! Number of right-hand sides
    integer, intent (in) :: ldb ! Leading extent of b
    real (wp), intent (inout) :: b(ldb,nrhs) ! Holds the right-hand
!                           sides on entry and is overwritten by the solution.
    integer, intent (out) :: flag
    integer(long), intent (in) :: ll ! Size of a
    real (wp), intent (in) :: a( ll )
!                               Holds ( L ) packed by columns.
!                                     ( M )

! .. Locals ..
    integer j ! Column index
    integer jb ! Width of block
    integer(long) k ! Position in array
    integer lj ! Length of block
    integer r ! Index for rhs

    flag = 0
    if (n < 0) then
       flag = -1
    else if (nb <= 1) then
       flag = -4
    else if (nrhs < 0) then
       flag = -5
    else if (ldb < n) then
       flag = -6
    else if (q < 0) then
       flag = -8
    else if (q > n) then
       flag = -9
    end if
    if (flag/=0 .or. q==0) return

     k = 1
     lj = n
     do j = 1, q, nb
        jb = min(nb,q-j+1)
        do r = 1,nrhs
           call dtpsv('L','N','U',jb,a(k),b(j,r),1)
        end do
        k = k + (jb*(jb+1_long))/2
        lj = lj - jb
        if (lj>0) call dgemm &
           ('N','N',lj,nrhs,jb,-unity,a(k),lj,b(j,1),ldb,unity,b(j+jb,1),ldb)
        k = k + jb*int(lj,long)
     end do

  end subroutine ma64_solveL2_double

  subroutine ma64_solveD2_double(n,q,nrhs,b,ldb,flag,d)
! Solves                ( D   ) x = b
!                       (   I )
! where D is block diagonal of size q with blocks of size 1 and 2.

   integer, intent (in) :: n ! Matrix order
    integer, intent (in) :: q ! Order of L (and D)
    integer, intent (in) :: nrhs ! Number of right-hand sides
    integer, intent (in) :: ldb ! Leading extent of b
    real (wp), intent (inout) :: b(ldb,nrhs) ! Holds the right-hand
!                           sides on entry and is overwritten by the solution.
    integer, intent (out) :: flag
    real (wp), intent (in) :: d(2*q) !  Holds the inverse of D,
!              except that zero diagonal blocks (pivots) are not inverted.
!              Diagonal entries are in d(1:2*q-1:2) and entries below
!              the diagonal are in d(2:2*q-2:2).

! .. Locals ..
    integer j ! Column index

    flag = 0
    if (n < 0) then
       flag = -1
    else if (nrhs < 0) then
       flag = -5
    else if (ldb < n) then
       flag = -6
    else if (q < 0) then
       flag = -8
    else if (q > n) then
       flag = -9
    end if
    if (flag/=0 .or. q==0) return

    do j = 1, nrhs
       call ma64_solveD1_double(n,q,b(1,j),flag,d)
    end do

  end subroutine ma64_solveD2_double

  subroutine ma64_solveLT2_double(n,q,nb,nrhs,b,ldb,flag,a,ll)
! Solves                ( L' M') x = b
!                       (    I )
! where  L is unit lower triangular of size q.

    integer, intent (in) :: n ! Matrix order
    integer, intent (in) :: q ! Order of L (and D)
    integer, intent (in) :: nb ! The block size used for the blocked hybrid
!                                format.
    integer, intent (in) :: nrhs ! Number of right-hand sides
    integer, intent (in) :: ldb ! Leading extent of b
    real (wp), intent (inout) :: b(ldb,nrhs) ! Holds the right-hand
!                           sides on entry and is overwritten by the solution.
    integer, intent (out) :: flag
    integer(long), intent (in) :: ll ! Size of a
    real (wp), intent (in) :: a( ll )
!                               Holds ( L ) packed by columns.
!                                     ( M )

! .. Locals ..
    integer j ! Column index
    integer jb ! Width of block
    integer j1 ! First column of block
    integer(long) k ! Position in array
    integer lj ! Length of block
    integer m !
    integer r ! Index for rhs

    flag = 0
    if (n < 0) then
       flag = -1
    else if (nb <= 1) then
       flag = -4
    else if (nrhs < 0) then
       flag = -5
    else if (ldb < n) then
       flag = -6
    else if (q < 0) then
       flag = -8
    else if (q > n) then
       flag = -9
    end if
    if (flag/=0 .or. q<1) return

! Back substitution
     m = nb*((q-1)/nb)
     k = 1 + ll
     j1 = 1 + m
     lj = n - q
     do j = j1, 1, -nb
        jb = min(nb,q-j+1)
        k = k - jb*int(lj,long)
        if (lj>0) call dgemm &
            ('T','N',jb,nrhs,lj,-unity,a(k),lj,b(j+jb,1),ldb,unity,b(j,1),ldb)
        k = k - (jb*(jb+1_long))/2
        do r = 1, nrhs
           call dtpsv('L','T','U',jb,a(k),b(j,r),1)
        end do
        lj = lj + jb
     end do

  end subroutine ma64_solveLT2_double

  subroutine ma64_solveDLT2_double(n,q,nb,nrhs,b,ldb,flag,a,ll,d)
! Solves       ( D   ) ( L' M') X = B
!              (   I ) (    I )
! where D is block diagonal of size q with blocks of size 1 and 2
! and L is unit lower triangular of size q.

    integer, intent (in) :: n ! Matrix order
    integer, intent (in) :: q ! Order of L (and D)
    integer, intent (in) :: nb ! The block size used for the blocked hybrid
!                                format.
    integer, intent (in) :: nrhs ! Number of right-hand sides
    integer, intent (in) :: ldb ! Leading extent of B
    real (wp), intent (inout) :: b(ldb,nrhs) ! Holds the right-hand
!                           sides on entry and is overwritten by the solution.
    integer, intent (out) :: flag
    integer(long), intent (in) :: ll ! Size of a
    real (wp), intent (in) :: a( ll )
!                               Holds ( L ) packed by columns.
!                                     ( M )
    real (wp), intent (in) :: d(2*q) !  Holds the inverse of D,
!              except that zero diagonal blocks (pivots) are not inverted.
!              Diagonal entries are in d(1:2*q-1:2) and entries below
!              the diagonal are in d(2:2*q-2:2).

    if (nb <= 1) then
       flag = -4
       return
    end if
    call ma64_solveD2_double(n,q,nrhs,b,ldb,flag,d)
    if (flag/=0) return
    call ma64_solveLT2_double(n,q,nb,nrhs,b,ldb,flag,a,ll)

  end subroutine ma64_solveDLT2_double


end module hsl_ma64_double


