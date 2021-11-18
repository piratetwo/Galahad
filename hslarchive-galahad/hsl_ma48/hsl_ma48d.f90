! COPYRIGHT (c) 2001 Council for the Central Laboratory
!               of the Research Councils
! Original date 18 October 2001
!
! Version 3.3.0
! For history see ChangeLog
!

!
! This file consists of four modules:
! hsl_ma48_ma50_internal_double (based on ma50)
! hsl_ma48_ma48_internal_double (based on ma48)
! hsl_ma48_ma51_internal_double (based on ma51)
! hsl_ma48_double (the user facing module)
!

module hsl_ma48_ma50_internal_double
   use hsl_zb01_double, zb01_real_info => zb01_info
   use hsl_zb01_integer, ZB01_integer => ZB01_resize1, &
       zb01_integer_info => zb01_info
   implicit none

! This module is based on ma50 2.0.0 (29 November 2006)
!
! 3 August 2010 Version 3.0.0.  Radical change to include some
!      INTEGER*64 integers (so that restriction on number of entries is
!      relaxed.  Also use of HSL_ZB01 to allow restarting with more
!      space allocated.

   private
   ! public for user
   public :: ma50ad, ma50bd, ma50cd
   ! public for unit tests
   public :: ma50dd, ma50ed, ma50fd, ma50gd, ma50hd, ma50id

   integer, parameter :: wp = kind(0.0d0)
   integer, parameter :: long = selected_int_kind(18) ! Long integer
   integer(long), parameter :: onel = 1
   integer :: stat

contains

SUBROUTINE ma50ad(m,n,ne,la,a,irn,ljcn,jcn,iq,cntl,icntl,ip,np, &
    jfirst,lenr,lastr, &
    nextr,ifirst,lenc,lastc,nextc,multiplier,info,rinfo)

! MA50A/AD chooses a pivot sequence using a Markowitz criterion with
!     threshold pivoting.

! If  the user requires a more convenient data interface then the MA48
!     package should be used. The MA48 subroutines call the MA50
!     subroutines after checking the user's input data and optionally
!     permute the matrix to block triangular form.

  integer m, n
  integer(long) :: ne, la, ljcn, oldlen
  real(wp), allocatable :: a(:)
  real(wp) :: cntl(10)
  real(wp), intent(in) :: multiplier
  integer, allocatable :: irn(:), jcn(:)
  integer, allocatable :: iw(:)
  integer icntl(20), np, jfirst(m), lenr(m), lastr(m), nextr(m), &
    ifirst(n), lenc(n), lastc(n), nextc(n)
  integer(long) :: info(15)
  integer(long) :: iq(n), ip(m)
  real(wp) :: rinfo(10)

! M is an integer variable that must be set to the number of rows.
!      It is not altered by the subroutine.
! N is an integer variable that must be set to the number of columns.
!      It is not altered by the subroutine.
! NE is an integer variable that must be set to the number of entries
!      in the input matrix. It is not altered by the subroutine.
! LA is an integer variable that must be set to the size of A, IRN, and
!      JCN. It is not altered by the subroutine.
! A is an array that holds the input matrix on entry and is used as
!      workspace.
! IRN  is an integer array.  Entries 1 to NE must be set to the
!      row indices of the corresponding entries in A.  IRN is used
!      as workspace and holds the row indices of the reduced matrix.
! JCN  is an integer array that need not be set by the user. It is
!      used to hold the column indices of entries in the reduced
!      matrix.
! IQ is an integer array of length N. On entry, it holds pointers
!      to column starts. During execution, IQ(j) holds the position of
!      the start of column j of the reduced matrix or -IQ(j) holds the
!      column index in the permuted matrix of column j. On exit, IQ(j)
!      holds the index of the column that is in position j of the
!      permuted matrix.
! CNTL must be set by the user as follows and is not altered.
!     CNTL(1)  Full matrix processing will be used if the density of
!       the reduced matrix is MIN(CNTL(1),1.0) or more.
!     CNTL(2) determines the balance between pivoting for sparsity and
!       for stability, values near zero emphasizing sparsity and values
!       near one emphasizing stability. Each pivot must have absolute
!       value at least CNTL(2) times the greatest absolute value in the
!       same column of the reduced matrix.
!     CNTL(3) If this is set to a positive value, any entry of the
!       reduced matrix whose modulus is less than CNTL(3) will be
!       dropped.
!     CNTL(4)  Any entry of the reduced matrix whose modulus is less
!       than or equal to CNTL(4) will be regarded as zero from the
!        point of view of rank.
! ICNTL must be set by the user as follows and is not altered.
!     ICNTL(1)  must be set to the stream number for error messages.
!       A value less than 1 suppresses output.
!     ICNTL(2) must be set to the stream number for diagnostic output.
!       A value less than 1 suppresses output.
!     ICNTL(3) must be set to control the amount of output:
!       0 None.
!       1 Error messages only.
!       2 Error and warning messages.
!       3 As 2, plus scalar parameters and a few entries of array
!         parameters on entry and exit.
!       4 As 3, plus all parameters on entry and exit.
!     ICNTL(4) If set to a positive value, the pivot search is limited
!       to ICNTL(4) columns (Zlatev strategy). This may result in
!       different fill-in and execution time. If ICNTL(4) is positive,
!       the workspace arrays LASTR and NEXTR are not referenced.
!     ICNTL(5) The block size to be used for full-matrix processing.
!     ICNTL(6) The last ICNTL(6) columns of A must be the last
!       ICNTL(6) columns of the permuted matrix. A value outside the
!       range 1 to N-1 is treated as zero.
!     ICNTL(7) If given the value 1, pivots are limited to
!       the main diagonal, which may lead to a premature switch to full
!       processing if no suitable diagonal entries are available.
!       If given the value 2, IFIRST must be set so that IFIRST(i) is
!       the column in position i of the permuted matrix and IP must
!       be set so that IP(i) < IP(j) if row i is recommended to
!       precede row j in the pivot sequence.
! IP is an integer array of length M that need not be set on entry
!      unless ICNTL(7)=2 (see ICNTL(7) for details of this case).
!      During execution, IP(i) holds the position of the start of row i
!      of the reduced matrix or -IP(i) holds the row index in the
!      permuted matrix of row i. Before exit, IP(i) is made positive.
! NP is an integer variable. It need not be set on entry. On exit,
!     it will be set to the number of columns to be processed in
!     packed storage.
! JFIRST is an integer workarray of length M. JFIRST(i) is the
!      first column of the reduced matrix to have i entries or is
!      zero if no column has i entries.
! LENR is an integer workarray of length M that is used to hold the
!      numbers of entries in the rows of the reduced matrix.
! LASTR is an integer workarray of length M, used only if ICNTL(4) = 0.
!      For rows in the reduced matrix, LASTR(i) indicates the previous
!      row to i with the same number of entries. LASTR(i) is zero if
!      no such row exists.
! NEXTR is an integer workarray of length M, used only if ICNTL(4) = 0
!      or ICNTL(7)=2. If ICNTL(4)=0, for rows in the reduced matrix,
!      NEXTR(i) indicates the next row to i with the same number of
!      entries; and if row i is the last in the chain, NEXTR is
!      equal to zero. If ICNTL(7)=2, NEXTR is a copy of the value of
!      IP on entry.
! IW is an integer array of length M used as workspace and is used to
!     assist the detection of duplicate entries and the sparse SAXPY
!     operations. It is reset to zero each time round the main loop.
! IFIRST is an integer array of length N, used only if ICNTL(4) = 0
!      or ICNTL(7)=2. If ICNTL(4) = 0, it is a workarray; IFIRST(i)
!      points to the first row of the reduced matrix to have i entries
!      or is zero if no row has i entries. If ICNTL(7)=2, IFIRST
!      must be set on entry (see ICNTL(7) for details of this case).
! LENC is an integer workarray of length N that is used to hold
!      the numbers of entries in the columns of the reduced matrix.
! LASTC is an integer workarray of length N.  For columns in the reduced
!      matrix, LASTC(j) indicates the previous column to j with the same
!      number of entries.  If column j is the first in the chain,
!      LASTC(j) is equal to zero.
! NEXTC is an integer workarray of length N.  For columns in the reduced
!      matrix, NEXTC(j) indicates the next column to j with the same
!      number of entries.  If column j is the last column in the chain,
!      NEXTC(j) is zero.
! INFO need not be set on entry. On exit, it holds the following:
!    INFO(1):
!       0  Successful entry.
!      -1  M < 1 or N < 1.
!      -2  NE < 1.
!      -4  Duplicated entries.
!      -5  Faulty column permutation in IFIRST when ICNTL(7)=2.
!      -6  ICNTL(4) not equal to 1 when ICNTL(7)=2.
!      +1  Rank deficient.
!      +2  Premature switch to full processing because of failure to
!          find a stable diagonal pivot (ICNTL(7)>=1 case only).
!      +3  Both of these warnings.
!    INFO(2) Number of compresses of the arrays.
!    INFO(3) Minimum LA recommended to analyse matrix.
!    INFO(4) Minimum LFACT required to factorize matrix.
!    INFO(5) Upper bound on the rank of the matrix.
!    INFO(6) Number of entries dropped from the data structure.
!    INFO(7) Number of rows processed in full storage.
! RINFO need not be set on entry. On exit, RINFO(1) holds the number of
!    floating-point operations needed for the factorization.

  integer idamax
  external idamax
  intrinsic abs, max, min

  real(wp) :: zero, one
  parameter (zero=0d0,one=1.0d0)

  type (zb01_real_info) :: zb01info
  type (zb01_integer_info) :: zb01intinfo
  real(wp) :: alen, amult, anew, asw, au, cost, cpiv
  integer len
  integer(long) :: eye, i, idummy, iend, ifill, ifir, ii, ij, &
    ijpos, iop, ipiv, ipos, isrch, i1, i2, j, jbeg, jend, jj, jlast, jmore, &
    jnew, jpiv, jpos, j1, j2, l, lc, lenpiv, lp, lr, ist
  integer(long) :: dispc, dispr, idrop
  real(wp) :: maxent
  integer(long) :: minc, mord, mp, msrch, nc, ndrop, nefact, nepr, nered, ne1, &
    nord, nord1, nr, nullc, nulli, nullj, nullr, pivbeg, pivcol, pivend, pivot
  real(wp) :: pivr, pivrat, u

! ALEN Real(LEN-1).
! AMULT Temporary variable used to store current multiplier.
! ANEW Temporary variable used to store value of fill-in.
! ASW Temporary variable used when swopping two real quantities.
! AU Temporary variable used in threshold test.
! COST Markowitz cost of current potential pivot.
! CPIV Markowitz cost of best pivot so far found.
! DISPC is the first free location in the column file.
! DISPR is the first free location in the row file.
! EYE Running relative position when processing pivot row.
! I Temporary variable holding row number. Also used as index in DO
!     loops used in initialization of arrays.
! IDROP Temporary variable used to accumulate number of entries dropped.
! IDUMMY DO index not referenced in the loop.
! IEND Position of end of pivot row.
! IFILL is the fill-in to the non-pivot column.
! IFIR Temporary variable holding first entry in chain.
! II Running position for current column.
! IJ Temporary variable holding row/column index.
! IJPOS Position of current pivot in A/IRN.
! IOP holds a running count of the number of rows with entries in both
!     the pivot and the non-pivot column.
! IPIV Row of the pivot.
! IPOS Temporary variable holding position in column file.
! ISRCH Temporary variable holding number of columns searched for pivot.
! I1 Position of the start of the current column.
! I2 Position of the end of the current column.
! J Temporary variable holding column number.
! JBEG Position of beginning of non-pivot column.
! JEND Position of end of non-pivot column.
! JJ Running position for current row.
! JLAST Last column acceptable as pivot.
! JMORE Temporary variable holding number of locations still needed
!     for fill-in in non-pivot column.
! JNEW Position of end of changed non-pivot column.
! JPIV Column of the pivot.
! JPOS Temporary variable holding position in row file.
! J1 Position of the start of the current row.
! J2 Position of the end of the current row.
! L Loop index.
! LC Temporary variable holding previous column in sequence.
! LEN Length of column or row.
! LENPIV Length of pivot column.
! LP Unit for error messages.
! LR Temporary variable holding previous row in sequence.
! MAXENT Temporary variable used to hold value of largest entry in
!    column.
! MINC Minimum number of entries of any row or column of the reduced
!     matrix, or in any column if ICNTL(4) > 0.
! MORD Number of rows ordered, excluding null rows.
! MP Unit for diagnostic messages.
! MSRCH Number of columns to be searched.
! NC Temporary variable holding next column in sequence.
! NDROP Number of entries dropped because of being in a column all of
!   whose entries are smaller than the pivot threshold.
! NEFACT Number of entries in factors.
! NEPR Number of entries in pivot row, excluding the pivot.
! NERED Number of entries in reduced matrix.
! NE1 Temporary variable used to hold number of entries in row/column
!     and to hold temporarily value of MINC.
! NORD Number of columns ordered, excluding null columns beyond JLAST.
! NORD1 Value of NORD at start of step.
! NR Temporary variable holding next row in sequence.
! NULLC Number of structurally zero columns found before any entries
!     dropped for being smaller than CNTL(3).
! NULLR Number of structurally zero rows found before any entries
!     dropped for being smaller than CNTL(3).
! NULLI Number of zero rows found.
! NULLJ Number of zero columns found beyond column JLAST.
! PIVBEG Position of beginning of pivot column.
! PIVCOL Temporary variable holding position in pivot column.
! PIVEND Position of end of pivot column.
! PIVOT Current step in Gaussian elimination.
! PIVR ratio of current pivot candidate to largest in its column.
! PIVRAT ratio of best pivot candidate to largest in its column.
! U Used to hold local copy of CNTL(2), changed if necessary so that it
!    is in range.

  lp = icntl(1)
  if (icntl(3)<=0) lp = 0
  mp = icntl(2)
  if (icntl(3)<=1) mp = 0
  info(1) = 0
  info(2) = 0
  info(3) = ne
  info(8) = ne
  info(4) = ne
  info(9) = ne
  info(5) = 0
  info(6) = 0
  info(7) = 0
  rinfo(1) = zero
  np = 0

! Make some simple checks
!!!
! Won't of course happen in the hsl_ma48 context
! if (m<1 .or. n<1) go to 690
  if (ne<1) go to 700
  if (la<ne) go to 710
!!    oldlen = la
!!    la = multiplier*la
!!write(6,*) 'oldlen,la',oldlen,la
!!    call ZB01_resize1(a,oldlen,la,zb01info)
!!write(6,*) 'oldlen,la',oldlen,la
!!      if (zb01info%flag .lt. 0) then
!!        info(1) = -2*n
!!        return
!!      endif
!!write(6,*) 'oldlen,la',oldlen,la
!!      call ZB01_integer(irn,oldlen,la,zb01intinfo)
!!write(6,*) 'oldlen,la',oldlen,la
!!      if (zb01intinfo%flag .lt. 0) then
!!        info(1) = -2*n
!!        return
!!      endif
!!  end if

  allocate(iw(m),stat=stat)
  if (stat .ne. 0) then
    info(1) = -2*n
    return
  endif

! Initial printing
  if (mp>0 .and. icntl(3)>2) then
    write (mp,'(/2(A,I6),A,I8,A,I8/A,1P,4E10.2/A,7I4)') &
      ' Entering MA50AD with M =', m, ' N =', n, ' NE =', ne, ' LA =', la, &
      ' CNTL =', (cntl(i),i=1,4), ' ICNTL =', (icntl(i),i=1,7)
    if (n==1 .or. icntl(3)>3) then
      do 10 j = 1, n - 1
        if (iq(j)<iq(j+1)) write (mp,'(A,I5,(T13,3(1P,E12.4,I5)))') ' Column', &
          j, (a(ii),irn(ii),ii=iq(j),iq(j+1)-1)
10    continue
      if (iq(n)<=ne) write (mp,'(A,I5,(T13,3(1P,E12.4,I5)))') ' Column', n, &
        (a(ii),irn(ii),ii=iq(n),ne)
    else
      if (iq(1)<iq(2)) write (mp,'(A,I5,(T13,3(1P,E12.4,I5)))') ' Column', 1, &
        (a(ii),irn(ii),ii=iq(1),iq(2)-1)
    end if
    if (icntl(7)==2) then
      write (mp,'(A,(T10,10(I7)))') ' IP = ', ip
      write (mp,'(A,(T10,10(I7)))') ' IFIRST = ', ifirst
    end if
  end if

! Initialization of counts etc.
  minc = 1
  nered = ne
  u = min(cntl(2),one)
  u = max(u,zero)
  msrch = icntl(4)
  if (msrch==0) msrch = n
  jlast = n - icntl(6)
  if (jlast<1 .or. jlast>n) jlast = n
  nulli = 0
  nullj = 0
  mord = 0
  nord = 0
  ndrop = 0
  nefact = 0
  do 20 i = 1, n - 1
    lenc(i) = iq(i+1) - iq(i)
20 continue
  lenc(n) = ne + 1 - iq(n)

  if (cntl(3)>zero) then
! Drop small entries
    nered = 0
    do 40 j = 1, n
      i = iq(j)
      iq(j) = nered + 1
      do 30 ii = i, i + lenc(j) - 1
        if (abs(a(ii))>=cntl(3)) then
          nered = nered + 1
          a(nered) = a(ii)
          irn(nered) = irn(ii)
        else
          info(6) = info(6) + 1
        end if
30    continue
      lenc(j) = nered + 1 - iq(j)
40  continue
  end if

  if (icntl(7)==2) then
! Column order specified - copy the row ordering array
    do 50 i = 1, m
      nextr(i) = ip(i)
50  continue
! Check ICNTL(4)
    if (icntl(4)/=1) go to 740
  end if

  dispr = nered + 1
  dispc = nered + 1

! Set up row oriented storage.
  do 60 i = 1, m
    iw(i) = 0
    lenr(i) = 0
    jfirst(i) = 0
60 continue
! Calculate row counts.
  do 70 ii = 1, nered
    i = irn(ii)
    lenr(i) = lenr(i) + 1
70 continue
! Set up row pointers so that IP(i) points to position after end
!     of row i in row file.
  ip(1) = lenr(1) + 1
  do 80 i = 2, m
    ip(i) = ip(i-1) + lenr(i)
80 continue
! Generate row file.
  do 100 j = 1, n
    ist = iq(j)
    do 90 ii = ist, ist + lenc(j) - 1
      i = irn(ii)
! Check for duplicate entry.
!!! Cannot happen when called from HSL_MA48
!     if (iw(i)==j) go to 720
      iw(i) = j
      ipos = ip(i) - 1
      jcn(ipos) = j
      ip(i) = ipos
90  continue
100 continue
  do 110 i = 1, m
    iw(i) = 0
110 continue

! Check for zero rows and (unless ICNTL(4) > 0), compute chains of rows
!    with equal numbers of entries.
  if (icntl(4)<=0) then
    do 120 i = 1, n
      ifirst(i) = 0
120 continue
    do 130 i = m, 1, -1
      ne1 = lenr(i)
      if (ne1>0) then
        ifir = ifirst(ne1)
        ifirst(ne1) = i
        lastr(i) = 0
        nextr(i) = ifir
        if (ifir>0) lastr(ifir) = i
      else
        ip(i) = -m + nulli
        nulli = nulli + 1
      end if
130 continue
  else
    do 140 i = m, 1, -1
      ne1 = lenr(i)
      if (ne1==0) then
        ip(i) = -m + nulli
        nulli = nulli + 1
      end if
140 continue
  end if
! Check for zero columns and compute chains of columns with equal
!   numbers of entries.
  do 150 j = n, 1, -1
    ne1 = lenc(j)
    if (ne1==0) then
      if (icntl(7)/=2) then
        if (j<=jlast) then
          nord = nord + 1
          iq(j) = -nord
          if (nord==jlast) then
! We have ordered the first N - ICNTL(6) columns.
            nord = nord + nullj
            jlast = n
            nullj = 0
          end if
        else
          nullj = nullj + 1
          iq(j) = -(jlast+nullj)
        end if
        lastc(j) = 0
        nextc(j) = 0
      end if
    else
      ifir = jfirst(ne1)
      jfirst(ne1) = j
      nextc(j) = ifir
      lastc(j) = 0
      if (ifir>0) lastc(ifir) = j
    end if
150 continue
  if (info(6)==0) then
    nullc = nord + nullj
    nullr = nulli
  end if


! **********************************************
! ****    Start of main elimination loop    ****
! **********************************************
  do 630 pivot = 1, n
! Check to see if reduced matrix should be considered as full.
    if (nered>=(min(cntl(1),one)*(n-nord))*(m-mord)) go to 640

    if (icntl(7)==2) then
! Column order specified - choose the pivot within the column
      ipiv = 0
      j = ifirst(pivot)
      if (j<1 .or. j>n) go to 730
      if (iq(j)<0) go to 730
      len = lenc(j)
      if (len<=0) go to 320
      alen = len - 1
      i1 = iq(j)
      i2 = i1 + len - 1
! Find largest entry in column
      ii = idamax(len,a(i1),1)
      maxent = abs(a(i1+ii-1))
! Is every entry in the column below the pivot threshold?
      if (maxent<=cntl(4)) go to 320
      au = max(maxent*u,cntl(4))
! Scan column for pivot
      do 160 ii = i1, i2
        if (abs(a(ii))<au) go to 160
! Candidate satisfies threshold criterion.
        i = irn(ii)
        if (ipiv/=0) then
          if (nextr(i)>=nextr(ipiv)) go to 160
        end if
        cpiv = alen*(lenr(i)-1)
        ijpos = ii
        ipiv = i
        jpiv = j
160   continue
      go to 330
    end if

! Find the least number of entries in a row or column (column only if
!   the Zlatev strategy is in use)
    len = minc
    do 170 minc = len, m - mord
      if (jfirst(minc)/=0) go to 180
      if (icntl(4)<=0) then
        if (ifirst(minc)/=0) go to 180
      end if
170 continue

! Find the next pivot or a column whose entries are all very small.
! CPIV is the Markowitz cost of the best pivot so far and PIVRAT is the
!      ratio of its absolute value to that of the largest entry in its
!      column.
180 cpiv = m
    cpiv = cpiv*n
    pivrat = zero
! Examine columns/rows in order of ascending count.
    isrch = 0
    do 300 len = minc, m - mord
      alen = len - 1
! Jump if Markowitz count cannot be bettered.
      if (cpiv<=alen**2 .and. icntl(4)<=0) go to 310
      ij = jfirst(len)
! Scan columns with LEN entries.
      do 220 idummy = 1, n
! If no more columns with LEN entries, exit loop.
        if (ij<=0) go to 230
        j = ij
        ij = nextc(j)
        if (j>jlast) go to 220
! Column J is now examined.
! First calculate multiplier threshold level.
        maxent = zero
        i1 = iq(j)
        i2 = i1 + len - 1
        ii = idamax(len,a(i1),1)
        maxent = abs(a(i1+ii-1))
! Exit loop if every entry in the column is below the pivot threshold.
        if (maxent<=cntl(4)) go to 320
        au = max(maxent*u,cntl(4))
! If diagonal pivoting requested, look for diagonal entry.
        if (icntl(7)==1) then
          do 190 ii = i1, i2
            if (irn(ii)==j) go to 200
190       continue
          go to 220
200       i1 = ii
          i2 = ii
        end if
! Scan column for possible pivots
        do 210 ii = i1, i2
          if (abs(a(ii))<au) go to 210
! Candidate satisfies threshold criterion.
          i = irn(ii)
          cost = alen*(lenr(i)-1)
          if (cost>cpiv) go to 210
          pivr = abs(a(ii))/maxent
          if (cost==cpiv) then
            if (pivr<=pivrat) go to 210
          end if
! Best pivot so far is found.
          cpiv = cost
          ijpos = ii
          ipiv = i
          jpiv = j
          if (cpiv<=alen**2 .and. icntl(4)<=0) go to 330
          pivrat = pivr
210     continue
! Increment number of columns searched.
        isrch = isrch + 1
! Jump if we have searched the number of columns stipulated and found a
!   pivot.
        if (isrch>=msrch) then
          if (pivrat>zero) go to 330
        end if
220   continue

! Rows with LEN entries now examined.
230   if (icntl(4)>0) go to 300
      if (cpiv<=alen*(alen+1)) go to 310
      if (len>n-nord) go to 300
      ij = ifirst(len)
      do 290 idummy = 1, m
        if (ij==0) go to 300
        i = ij
        ij = nextr(ij)
        j1 = ip(i)
        j2 = j1 + len - 1
! If diagonal pivoting requested, look for diagonal entry.
        if (icntl(7)==1) then
          do 240 jj = j1, j2
            if (jcn(jj)==i) go to 250
240       continue
          go to 290
250       j1 = jj
          j2 = jj
        end if
! Scan row I.
        do 280 jj = j1, j2
          j = jcn(jj)
          if (j>jlast) go to 280
          cost = alen*(lenc(j)-1)
          if (cost>=cpiv) go to 280
! Pivot has best Markowitz count so far. Now check its suitability
!     on numerical grounds by examining other entries in its column.
          i1 = iq(j)
          i2 = i1 + lenc(j) - 1
          ii = idamax(lenc(j),a(i1),1)
          maxent = abs(a(i1+ii-1))
          do 260 ii = i1, i2 - 1
            if (irn(ii)==i) go to 270
260       continue
270       jpos = ii
! Exit loop if every entry in the column is below the pivot threshold.
          if (maxent<=cntl(4)) go to 320
          if (abs(a(jpos))<maxent*u) go to 280
! Candidate satisfies threshold criterion.
          cpiv = cost
          ipiv = i
          jpiv = j
          ijpos = jpos
          pivrat = abs(a(jpos))/maxent
          if (cpiv<=alen*(alen+1)) go to 330
280     continue

290   continue

300 continue
310 if (pivrat>zero) go to 330
! No pivot found. Switch to full matrix processing.
    info(1) = info(1) + 2
    if (mp>0) write (mp,'(A/A)') &
      ' Warning message from MA50AD: no suitable diagonal pivot', &
      ' found, so switched to full matrix processing.'
    go to 640

! Every entry in the column is below the pivot threshold.
320 ipiv = 0
    jpiv = j

! The pivot has now been found in position (IPIV,JPIV) in location
!     IJPOS in column file or all entries of column JPIV are very small
!     (IPIV=0).
! Update row and column ordering arrays to correspond with removal
!     of the active part of the matrix. Also update NEFACT.
330 nefact = nefact + lenc(jpiv)
    pivbeg = iq(jpiv)
    pivend = pivbeg + lenc(jpiv) - 1
    nord = nord + 1
    nord1 = nord
    if (nord==jlast) then
! We have ordered the first N - ICNTL(6) columns.
      nord = nord + nullj
      jlast = n
      nullj = 0
    end if
    if (icntl(4)<=0) then
! Remove active rows from their row ordering chains.
      do 340 ii = pivbeg, pivend
        i = irn(ii)
        lr = lastr(i)
        nr = nextr(i)
        if (nr/=0) lastr(nr) = lr
        if (lr==0) then
          ne1 = lenr(i)
          ifirst(ne1) = nr
        else
          nextr(lr) = nr
        end if
340   continue
    end if
    if (ipiv>0) then
! NEPR is number of entries in strictly U part of pivot row.
      nepr = lenr(ipiv) - 1
      nefact = nefact + nepr
      rinfo(1) = rinfo(1) + cpiv*2 + lenr(ipiv)
      j1 = ip(ipiv)
! Remove active columns from their column ordering chains.
      do 350 jj = j1, j1 + nepr
        j = jcn(jj)
        lc = lastc(j)
        nc = nextc(j)
        if (nc/=0) lastc(nc) = lc
        if (lc==0) then
          ne1 = lenc(j)
          jfirst(ne1) = nc
        else
          nextc(lc) = nc
        end if
350   continue
! Move pivot to beginning of pivot column.
      if (pivbeg/=ijpos) then
        asw = a(pivbeg)
        a(pivbeg) = a(ijpos)
        a(ijpos) = asw
        irn(ijpos) = irn(pivbeg)
        irn(pivbeg) = ipiv
      end if
    else
      nepr = 0
      ne1 = lenc(jpiv)
      if (cntl(3)>zero) ndrop = ndrop + ne1
      if (ne1>0) then
! Remove column of small entries from its column ordering chain.
        lc = lastc(jpiv)
        nc = nextc(jpiv)
        if (nc/=0) lastc(nc) = lc
        if (lc==0) then
          jfirst(ne1) = nc
        else
          nextc(lc) = nc
        end if
      end if
    end if

! Set up IW array so that IW(i) holds the relative position of row i
!    entry from beginning of pivot column.
    do 360 ii = pivbeg + 1, pivend
      i = irn(ii)
      iw(i) = ii - pivbeg
360 continue
! LENPIV is length of strictly L part of pivot column.
    lenpiv = pivend - pivbeg

! Remove pivot column (including pivot) from row oriented file.
    do 390 ii = pivbeg, pivend
      i = irn(ii)
      lenr(i) = lenr(i) - 1
      j1 = ip(i)
! J2 is last position in old row.
      j2 = j1 + lenr(i)
      do 370 jj = j1, j2 - 1
        if (jcn(jj)==jpiv) go to 380
370   continue
380   jcn(jj) = jcn(j2)
      jcn(j2) = 0
390 continue

! For each active column, add the appropriate multiple of the pivot
!     column to it.
! We loop on the number of entries in the pivot row since the position
!     of this row may change because of compresses.
    do 600 eye = 1, nepr
      j = jcn(ip(ipiv)+eye-1)
! Search column J for entry to be eliminated, calculate multiplier,
!     and remove it from column file.
!  IDROP is the number of nonzero entries dropped from column J
!        because these fall beneath tolerance level.
      idrop = 0
      jbeg = iq(j)
      jend = jbeg + lenc(j) - 1
      do 400 ii = jbeg, jend - 1
        if (irn(ii)==ipiv) go to 410
400   continue
410   amult = -a(ii)/a(iq(jpiv))
      a(ii) = a(jend)
      irn(ii) = irn(jend)
      lenc(j) = lenc(j) - 1
      irn(jend) = 0
      jend = jend - 1
! Jump if pivot column is a singleton.
      if (lenpiv==0) go to 600
! Now perform necessary operations on rest of non-pivot column J.
      iop = 0
! Innermost loop.
!DIR$ IVDEP
      do 420 ii = jbeg, jend
        i = irn(ii)
        if (iw(i)>0) then
! Row i is involved in the pivot column.
          iop = iop + 1
          pivcol = iq(jpiv) + iw(i)
! Flag IW(I) to show that the operation has been done.
          iw(i) = -iw(i)
          a(ii) = a(ii) + amult*a(pivcol)
        end if
420   continue

      if (cntl(3)>zero) then
!  Run through non-pivot column compressing column so that entries less
!      than CNTL(3) are not stored. All entries less than CNTL(3) are
!      also removed from the row structure.
        jnew = jbeg
        do 450 ii = jbeg, jend
          if (abs(a(ii))>=cntl(3)) then
            a(jnew) = a(ii)
            irn(jnew) = irn(ii)
            jnew = jnew + 1
          else
!  Remove non-zero entry from row structure.
            i = irn(ii)
            j1 = ip(i)
            j2 = j1 + lenr(i) - 1
            do 430 jj = j1, j2 - 1
              if (jcn(jj)==j) go to 440
430         continue
440         jcn(jj) = jcn(j2)
            jcn(j2) = 0
            lenr(i) = lenr(i) - 1
          end if
450     continue
        do 460 ii = jnew, jend
          irn(ii) = 0
460     continue
        idrop = jend + 1 - jnew
        jend = jnew - 1
        lenc(j) = lenc(j) - idrop
        nered = nered - idrop
        info(6) = info(6) + idrop
      end if

! IFILL is fill-in left to do to non-pivot column J.
      ifill = lenpiv - iop
      nered = nered + ifill
      info(3) = max(info(3),nered + lenc(j))

! Treat no-fill case
      if (ifill==0) then
!DIR$ IVDEP
        do 470 ii = pivbeg + 1, pivend
          i = irn(ii)
          iw(i) = -iw(i)
470     continue
        go to 600
      end if

! See if there is room for fill-in at end of the column.
      do 480 ipos = jend + 1, min(jend + ifill,dispc-1)
        if (irn(ipos)/=0) go to 490
480   continue
      if (ipos==jend +ifill+1) go to 540
      if (jend +ifill+1<=la+1) then
        dispc = jend + ifill + 1
        go to 540
      end if
      ipos = la
      dispc = la + 1
! JMORE more spaces for fill-in are required.
490   jmore = jend + ifill - ipos + 1
! We now look in front of the column to see if there is space for
!     the rest of the fill-in.
      do 500 ipos = jbeg - 1, max(jbeg-jmore,1_long), -1
        if (irn(ipos)/=0) go to 510
500   continue
      ipos = ipos + 1
      if (ipos==jbeg-jmore) go to 520
! Column must be moved to the beginning of available storage.
510   if (dispc+lenc(j)+ifill>la+1) then
        info(2) = info(2) + 1
! Call compress and reset pointers
        call ma50dd(a,irn,iq,n,dispc,.true.)
        jbeg = iq(j)
        jend = jbeg + lenc(j) - 1
        pivbeg = iq(jpiv)
        pivend = pivbeg + lenc(jpiv) - 1
!!!! Call to HSL_ZB01
! Was  ... go to 705
        if (dispc+lenc(j)+ifill>la+1) then
          oldlen = la
          la = multiplier*real(la,kind=wp)
          call ZB01_resize1(a,oldlen,la,zb01info)
          if (zb01info%flag .lt. 0) then
!XXX In the test deck, size(a) can be larger than (new) la
!    This will provoke an error return of -1 from HSL_ZB01
!    It can also be caused by an allocation error in HSL_ZB01
            info(1) = -2*n
            return
          endif
          call ZB01_integer(irn,oldlen,la,zb01intinfo)
          if (zb01intinfo%flag .lt. 0) then
! It can be caused by an allocation error in HSL_ZB01
            info(1) = -2*n
            return
          endif
        endif
      end if
      ipos = dispc
      dispc = dispc + lenc(j) + ifill
! Move non-pivot column J.
520   iq(j) = ipos
      do 530 ii = jbeg, jend
        a(ipos) = a(ii)
        irn(ipos) = irn(ii)
        ipos = ipos + 1
        irn(ii) = 0
530   continue
      jbeg = iq(j)
      jend = ipos - 1
! Innermost fill-in loop which also resets IW.
! We know at this stage that there are IFILL positions free after JEND.
540   idrop = 0
      do 580 ii = pivbeg + 1, pivend
        i = irn(ii)
!       info(3) = max(info(3),nered +lenr(i)+1)
        info(8) = max(info(8),nered +lenr(i)+1)
        if (iw(i)<0) then
          iw(i) = -iw(i)
          go to 580
        end if
        anew = amult*a(ii)
        if (abs(anew)<cntl(3)) then
          idrop = idrop + 1
        else
          jend = jend + 1
          a(jend) = anew
          irn(jend) = i

! Put new entry in row file.
! iend is end of row i with the fill-in included
          iend = ip(i) + lenr(i)
          if (iend<dispr) then
! Jump if can copy new entry to position iend
            if (jcn(iend)==0) go to 560
          else
            if (dispr<=ljcn) then
              dispr = dispr + 1
              go to 560
            end if
          end if
          if (ip(i)>1) then
! Check if there is a free space before the row in the row file
            if (jcn(ip(i)-1)==0) then
! Copy row forward
              iend = iend - 1
              do 545 jj = ip(i), iend
                jcn(jj-1) = jcn(jj)
545           continue
              ip(i) = ip(i) - 1
              go to 560
            end if
          end if
! The full row must be copied to end of row file storage
          if (dispr+lenr(i)>ljcn) then
! Compress row file array.
            info(2) = info(2) + 1
            call ma50dd(a,jcn,ip,m,dispr,.false.)
!!!! Call to HSL_ZB01
! Was earlier go to 705
            if (dispr+lenr(i)>ljcn) then
              oldlen = ljcn
              ljcn = 2*ljcn
              call ZB01_integer(jcn,oldlen,ljcn,zb01intinfo)
              if (zb01intinfo%flag .lt. 0) then
! It can be caused by an allocation error in HSL_ZB01
                info(1) = -2*n
                return
              endif
            endif
          end if
! Copy row to first free position.
          j1 = ip(i)
          j2 = ip(i) + lenr(i) - 1
          ip(i) = dispr
          do 550 jj = j1, j2
            jcn(dispr) = jcn(jj)
            jcn(jj) = 0
            dispr = dispr + 1
550       continue
          iend = dispr
          dispr = iend + 1
! Fill-in put in position iend in row file
560       jcn(iend) = j
          lenr(i) = lenr(i) + 1
! End of adjustment to row file.
        end if
580   continue
      info(6) = info(6) + idrop
      nered = nered - idrop
      do 590 ii = 1, idrop
        irn(jend +ii) = 0
590   continue
      lenc(j) = lenc(j) + ifill - idrop
! End of scan of pivot row.
600 continue


! Remove pivot row from row oriented storage and update column
!     ordering arrays.  Remember that pivot row no longer includes
!     pivot.
    do 610 eye = 1, nepr
      jj = ip(ipiv) + eye - 1
      j = jcn(jj)
      jcn(jj) = 0
      ne1 = lenc(j)
      lastc(j) = 0
      if (ne1>0) then
        ifir = jfirst(ne1)
        jfirst(ne1) = j
        nextc(j) = ifir
        if (ifir/=0) lastc(ifir) = j
        minc = min(minc,ne1)
      else if (icntl(7)/=2) then
        if (info(6)==0) nullc = nullc + 1
        if (j<=jlast) then
          nord = nord + 1
          iq(j) = -nord
          if (nord==jlast) then
! We have ordered the first N - ICNTL(6) columns.
            nord = nord + nullj
            jlast = n
            nullj = 0
          end if
        else
          nullj = nullj + 1
          iq(j) = -(jlast+nullj)
        end if
      end if
610 continue
    nered = nered - nepr

! Restore IW and remove pivot column from column file.
!    Record the row permutation in IP(IPIV) and the column
!    permutation in IQ(JPIV), flagging them negative so that they
!    are not confused with real pointers in compress routine.
    if (ipiv/=0) then
      lenr(ipiv) = 0
      iw(ipiv) = 0
      irn(pivbeg) = 0
      mord = mord + 1
      pivbeg = pivbeg + 1
      ip(ipiv) = -mord
    end if
    nered = nered - lenpiv - 1
    do 620 ii = pivbeg, pivend
      i = irn(ii)
      iw(i) = 0
      irn(ii) = 0
      ne1 = lenr(i)
      if (ne1==0) then
        if (info(6)==0) nullr = nullr + 1
        ip(i) = -m + nulli
        nulli = nulli + 1
      else if (icntl(4)<=0) then
! Adjust row ordering arrays.
        ifir = ifirst(ne1)
        lastr(i) = 0
        nextr(i) = ifir
        ifirst(ne1) = i
        if (ifir/=0) lastr(ifir) = i
        minc = min(minc,ne1)
      end if
620 continue
    iq(jpiv) = -nord1
630 continue
! We may drop through this loop with NULLI nonzero.

! ********************************************
! ****    End of main elimination loop    ****
! ********************************************

! Complete the permutation vectors
640 info(5) = mord + min(m-mord-nulli,n-nord-nullj)
  do 650 l = 1, min(m-mord,n-nord)
    rinfo(1) = rinfo(1) + real(m - mord - l + 1,kind=wp) +  &
                          real(m-mord-l,kind=wp)*real(n-nord-l,kind=wp)*2.0d0
650 continue
  np = nord
! Estimate of reals for subsequent factorization
!! Because we can't easily decouple real/integer with BTF structure
!! We leave info(4) to be the maximum value
!!info(4) = 2 + nefact + m*2 + (n-nord)*(m-mord)
  info(4) = 2 + nefact + m*2 + max((n-nord)*(m-mord),n-nord +m-mord)
! Estimate of integers for subsequent factorization
  info(9) = 2 + nefact + m*2 + n-nord +m-mord
  info(6) = info(6) + ndrop
  info(7) = m - mord
  do 660 l = 1, m
    if (ip(l)<0) then
      ip(l) = -ip(l)
    else
      mord = mord + 1
      ip(l) = mord
    end if
660 continue
  do 670 l = 1, n
    if (iq(l)<0) then
      lastc(l) = -iq(l)
    else
      if (nord==jlast) nord = nord + nullj
      nord = nord + 1
      lastc(l) = nord
    end if
670 continue
! Store the inverse permutation
  do 680 l = 1, n
    iq(lastc(l)) = l
680 continue

! Test for rank deficiency
  if (info(5)<min(m,n)) info(1) = info(1) + 1

  if (mp>0 .and. icntl(3)>2) then
    write (mp,'(A,I6,A,F12.1/A,7I8)') ' Leaving MA50AD with NP =', np, &
      ' RINFO(1) =', rinfo(1), ' INFO =', (info(i),i=1,7)
    if (icntl(3)>3) then
      write (mp,'(A,(T6,10(I7)))') ' IP = ', ip
      write (mp,'(A,(T6,10(I7)))') ' IQ = ', iq
    end if
  end if

  go to 750

! Error conditions.
! Can't happen if call through HSL_MA48
!690 info(1) = -1
!  if (lp>0) write (lp,'(/A/(2(A,I8)))') ' **** Error return from MA50AD ****',&
!    ' M =', m, ' N =', n
!  go to 750
700 info(1) = -2
  if (lp>0) write (lp,'(/A/(A,I10))') ' **** Error return from MA50AD ****', &
    ' NE =', ne
  go to 750
!XXX Use of hsl_ZB01 prevents this problem
!!705 info(4) = nefact + nered
!!  info(6) = info(6) + ndrop
  710 info(1) = -3
    if (lp>0) write (lp,'(/A/A,I9,A,I9)') &
      ' **** Error return from MA50AD ****', &
      ' LA  must be increased from', la, ' to at least', info(3)
    go to 750
!!!
! By this time all duplicates should have been removed but am leaving
! this in in case MA50 is tested separately.
!720 info(1) = -4
!  if (lp>0) write (lp,'(/A/(3(A,I9)))') ' **** Error return from MA50AD ****',&
!    ' Entry in row', i, ' and column', j, ' duplicated'
!  go to 750
730 info(1) = -5
  if (lp>0) write (lp,'(/A/(3(A,I9)))') ' **** Error return from MA50AD ****', &
    ' Fault in component ', pivot, ' of column permutation given in IFIRST'
  go to 750
740 info(1) = -6
  if (lp>0) write (lp,'(/A/(3(A,I9)))') ' **** Error return from MA50AD ****', &
    ' ICNTL(4) = ', icntl(4), ' when ICNTL(6) = 2'
750 end  subroutine ma50ad


SUBROUTINE ma50bd(m,n,ne,job,aa,irna,iptra,cntl,icntl,ip,iq,np,lfact,fact, &
    lirnf,irnf,iptrl,iptru,info,rinfo)
! MA50B/BD factorizes the matrix in AA/IRNA/IPTRA as P L U Q where
!     P and Q are permutations, L is lower triangular, and U is unit
!     upper triangular. The prior information that it uses depends on
!     the value of the parameter JOB.

  integer m, n, job
  integer(long) ::  ne
  real(wp) :: aa(ne)
  integer(long) :: irna(ne)
  integer(long) :: iptra(n)
  real(wp) :: cntl(10)
  integer :: icntl(20), np
  integer(long) :: ip(m), iq(*), lfact, lirnf
!! At the moment, lfact = lirnf for this to work and especially for
!! the blocking (for BTF) in ma48bd to work.  We may change this later.
  real(wp) :: fact(lfact)
  integer(long) irnf(lirnf)
  integer(long) :: iptrl(n), iptru(n)
  real(wp), allocatable :: w(:)
  integer(long), allocatable :: iw(:)
  integer(long) :: info(15)
  real(wp) :: rinfo(10)

! M is an integer variable that must be set to the number of rows.
!      It is not altered by the subroutine.
! N is an integer variable that must be set to the number of columns.
!      It is not altered by the subroutine.
! NE is an integer variable that must be set to the number of entries
!      in the input matrix.  It is not altered by the subroutine.
! JOB is an integer variable that must be set to the value 1, 2, or 3.
!     If JOB is equal to 1 and any of the first NP recommended pivots
!      fails to satisfy the threshold pivot tolerance, the row is
!      interchanged with the earliest row in the recommended sequence
!      that does satisfy the tolerance. Normal row interchanges are
!      performed in the last N-NP columns.
!     If JOB is equal to 2, then M, N, NE, IRNA, IPTRA, IP, IQ,
!      LFACT, NP, IRNF, IPTRL, and IPTRU must be unchanged since a
!      JOB=1 entry for the same matrix pattern and no interchanges are
!      performed among the first NP pivots; if ICNTL(6) > 0, the first
!      N-ICNTL(6) columns of AA must also be unchanged.
!     If JOB is equal to 3, ICNTL(6) must be in the range 1 to N-1.
!      The effect is as for JOB=2 except that interchanges are
!      performed.
!     JOB is not altered by the subroutine.
! AA is an array that holds the entries of the matrix and
!      is not altered.
! IRNA is an integer array of length NE that must be set to hold the
!      row indices of the corresponding entries in AA. It is not
!      altered.
! IPTRA is an integer array that holds the positions of the starts of
!      the columns of AA. It is not altered by the subroutine.
! CNTL  must be set by the user as follows and is not altered.
!     CNTL(2) determines the balance between pivoting for sparsity and
!       for stability, values near zero emphasizing sparsity and values
!       near one emphasizing stability.
!     CNTL(3) If this is set to a positive value, any entry whose
!       modulus is less than CNTL(3) will be dropped from the factors.
!       The factorization will then require less storage but will be
!       inaccurate.
!     CNTL(4)  Any entry of the reduced matrix whose modulus is less
!       than or equal to CNTL(4) will be regarded as zero from the
!        point of view of rank.
! ICNTL must be set by the user as follows and is not altered.
!     ICNTL(1)  must be set to the stream number for error messages.
!       A value less than 1 suppresses output.
!     ICNTL(2) must be set to the stream number for diagnostic output.
!       A value less than 1 suppresses output.
!     ICNTL(3) must be set to control the amount of output:
!       0 None.
!       1 Error messages only.
!       2 Error and warning messages.
!       3 As 2, plus scalar parameters and a few entries of array
!         parameters on entry and exit.
!       4 As 2, plus all parameters on entry and exit.
!     ICNTL(5) The block size to be used for full-matrix processing.
!       If <=0, the BLAS1 version is used.
!       If =1, the BLAS2 version is used.
!     ICNTL(6) If N > ICNTL(6) > 0, only the columns of A that
!       correspond to the last ICNTL(6) columns of the permuted matrix
!       may change prior to an entry with JOB > 1.
!     ICNTL(8) If this has value 1, there is no attempt to compute a
!       recommended value for LFACT if it is too small.
! IP is an integer array. If JOB=1, it must be set so that IP(I) < IP(J)
!      if row I is recommended to precede row J in the pivot sequence.
!      If JOB>1, it need not be set. If JOB=1 or JOB=3, IP(I) is set
!      to -K when row I is chosen for pivot K and IP is eventually
!      reset to recommend the chosen pivot sequence to a subsequent
!      JOB=1 entry. If JOB=2, IP is not be referenced.
! IQ is an integer array that must be set so that either IQ(J) is the
!      column in position J in the pivot sequence, J=1,2,...,N,
!      or IQ(1)=0 and the columns are taken in natural order.
!      It is not altered by the subroutine.
! NP is an integer variable that holds the number of columns to be
!      processed in packed storage. It is not altered by the subroutine.
! LFACT is an integer variable set to the size of FACT and IRNF.
!      It is not altered by the subroutine.
! FACT is an array that need not be set on a JOB=1 entry and must be
!      unchanged since the previous entry if JOB>1. On return, FACT(1)
!      holds the value of CNTL(3) used, FACT(2) will holds the value
!      of CNTL(4) used, FACT(3:IPTRL(N)) holds the packed part of L/U
!      by columns, and the full part of L/U is held by columns
!      immediately afterwards. U has unit diagonal entries, which are
!      not stored. In each column of the packed part, the entries of
!      U precede the entries of L; also the diagonal entries of L
!      head each column of L and are reciprocated.
! IRNF is an integer array of length LFACT that need not be set on
!      a JOB=1 entry and must be unchanged since the previous entry
!      if JOB>1. On exit, IRNF(1) holds the number of dropped entries,
!      IRNF(2) holds the number of rows MF in full storage,
!      IRNF(3:IPTRL(N)) holds the row numbers of the packed part
!      of L/U, IRNF(IPTRL(N)+1:IPTRL(N)+MF) holds the row indices
!      of the full part of L/U, and IRNF(IPTRL(N)+MF+I), I=1,2,..,N-NP
!      holds the vector IPIV output by MA50GD.
!      If JOB=2, IRNF will be unaltered.
! IPTRL is an integer array that need not be set on a JOB=1 entry and
!     must be unchanged since the previous entry if JOB>1.
!     For J = 1,..., NP, IPTRL(J) holds the position in
!     FACT and IRNF of the end of column J of L.
!     For J = NP+1,..., N, IPTRL(J) is equal to IPTRU(J).
! IPTRU is an integer array that need not be set on a JOB=1 entry and
!     must be unchanged since the previous entry if JOB>1.
!     For J = 1,..., N, IPTRU(J) holds the position in
!     FACT and IRNF of the end of the packed part of column J of U.
! W is an array of length M used as workspace for holding
!      the expanded form of a sparse vector.
! IW is an integer array of length M+2*N used as workspace.
! INFO need not be set on entry. On exit, it holds the following:
!    INFO(1) A negative value will indicate an error return and a
!       positive value a warning. Possible nonzero values are:
!      -1  M < 1 or N < 1.
!      -2  NE < 0.
!      -3  Insufficient space.
!      -4  There are duplicated entries.
!      -5  JOB < 1, 3 when ICNTL(6)=0, or > 3.
!      -6  JOB = 2, but entries were dropped in the corresponding JOB=1
!          entry.
!      -7  NP < 0 or NP > N.
!     -(7+K) Pivot too small in column K when JOB=2.
!      +1  Rank deficient.
!    INFO(4) Minimum storage required to factorize matrix
!            if INFO(1) >= 0. Recommended value for LFACT
!            if ICNTL(8) = 0 and INFO(1) = -3.
!    INFO(5) Computed rank of the matrix.
!    INFO(6) Number of entries dropped from the data structure.
!    INFO(7) Number of rows processed in full storage.
! RINFO need not be set on entry. On exit, RINFO(1) holds the number of
!    floating-point operations performed.

  real(wp) :: zero, one
  parameter (zero=0.0_wp,one=1.0_wp)
  real(wp) :: amult, asw
  integer(long) :: begcol
  logical drop
  integer(long) :: endcol, eye, eye1, i, ia1, ia2, if1, if2, ii,   &
    il1, il2, ipiv, iqpiv, &
    iu1, iu2, isw, j, jdummy, jj, lp
  real(wp) :: maxent
  integer :: k, jlast, mf, mord, mp, nf, nullc
  integer(long) :: neu
  real(wp) :: pivlim
  integer :: rank
  real(wp) :: u
! AMULT Temporary variable used to store current multiplier.
! ASW Temporary variable used when swopping two real quantities.
! BEGCOL is pointer to beginning of section of column when pruning.
! DROP True if any entries dropped from current column.
! ENDCOL is pointer to end of section of column when pruning.
! EYE Running position for current column.
! EYE1 Position of the start of second current column.
! I Temporary variable holding row number. Also used as index in DO
!     loops used in initialization of arrays.
! IA1 Position of the start of the current column in AA.
! IA2 Position of the end of the current column in AA.
! IF1 Position of the start of the full submatrix.
! IF2 Position of the end of the full submatrix.
! II Running position for current column.
! IL1 Position of the first entry of the current column of L.
! IL2 Position of the last entry of the current column of L.
! IPIV Position of the pivot in FACT and IRNF.
! IQPIV Recommended position of the pivot row in the pivot sequence.
! IU1 Position of the start of current column of U.
! IU2 Position of the end of the current column of U.
! ISW Temporary variable used when swopping two integer quantities.
! J Temporary variable holding column number.
! JDUMMY DO index not referenced in the loop.
! JJ Running position for current column.
! JLAST The lesser of NP and the last column of A for which no new
!     factorization operations are needed.
! K Temporary variable holding the current pivot step in the elimination
! LP Unit for error messages.
! MAXENT Temporary variable used to hold value of largest entry in
!    column.
! MF Number of rows in full block.
! MORD Number of rows ordered.
! MP Unit for diagnostic messages.
! NEU Number of entries omitted from U and the full block in order to
!    calculate INFO(4) (0 unless INFO(1)=-3).
! NF Number of columns in full block.
! NULLC Number of columns found null before dropping any elements.
! PIVLIM Limit on pivot size.
! RANK Value returned by MA50E/ED or MA50F/FD
! U Used to hold local copy of CNTL(2), changed if necessary so that it
!    is in range.

  intrinsic abs, max, min
! LAPACK subroutine for triangular factorization.

  info(1) = 0
  info(4) = 0
  info(5) = 0
  info(6) = 0
  info(7) = 0
  info(9) = 0
  rinfo(1) = zero
  lp = icntl(1)
  mp = icntl(2)
  if (icntl(3)<=0) lp = 0
  if (icntl(3)<=1) mp = 0

! Check input values
!!!
! Won't of course happen in hsl_ma48 context
! if (m<1 .or. n<1) then
!   info(1) = -1
!   if (lp>0) write (lp,'(/A/A,I8,A,I8)') &
!     ' **** Error return from MA50BD ****', ' M =', m, ' N =', n
!   go to 550
! end if
  if (ne<=0) then
    info(1) = -2
    if (lp>0) write (lp,'(/A/A,I6)') ' **** Error return from MA50BD ****', &
      ' NE =', ne
    go to 550
  end if
  if (np<0 .or. np>n) then
    info(1) = -7
    if (lp>0) write (lp,'(/A/A,I8,A,I8)') &
      ' **** Error return from MA50BD ****', ' NP =', np, ' N =', n
    go to 550
  end if
  if (lfact<max(m*onel,ne+2)) then
    info(4) = max(m*onel,ne+2)
    go to 520
  end if
  if (job==1) then
  else if (job==2 .or. job==3) then
    if (irnf(1)/=0) then
      info(1) = -6
      if (lp>0) write (lp,'(/A/A,I1,A)') ' **** Error return from MA50BD ***', &
        ' Call with JOB=', job, &
        ' follows JOB=1 call in which entries were dropped'
      go to 550
    end if
  else
    info(1) = -5
    if (lp>0) write (lp,'(/A/A,I2)') ' **** Error return from MA50BD ****', &
      ' JOB =', job
    go to 550
  end if

! Allocate work arrays
  allocate(w(m),stat = stat)
  if (stat .ne. 0) then
    info(1) = -2*n
    return
  endif
  allocate(iw(m+2*n),stat = stat)
  if (stat .ne. 0) then
    info(1) = -2*n
    return
  endif

! Print input data
  if (mp>0) then
    if (icntl(3)>2) write (mp, &
      '(/2(A,I6),A,I8,A,I3/A,I8,A,I7/A,1P,4E10.2/A,7I8)') &
      ' Entering MA50BD with M =', m, ' N =', n, ' NE =', ne, ' JOB =', job, &
      ' LFACT =', lfact, ' NP =', np, ' CNTL =', (cntl(i),i=1,4), ' ICNTL =', &
      (icntl(i),i=1,7)
    if (icntl(3)>3) then
      write (mp,'(A,(T6,10(I7)))') ' IP = ', ip
      if (iq(1)>0) then
        write (mp,'(A,(T6,10(I7)))') ' IQ = ', (iq(j),j=1,n)
      else
        write (mp,'(A,(T6,I7))') ' IQ = ', iq(1)
      end if
      do 10 j = 1, n - 1
        if (iptra(j)<iptra(j+1)) write (mp,'(A,I5,(T13,3(1P,E12.4,I5)))') &
          ' Column', j, (aa(ii),irna(ii),ii=iptra(j),iptra(j+1)-1)
10    continue
      if (iptra(n)<=ne) write (mp,'(A,I5,(T13,3(1P,E12.4,I5)))') ' Column', n, &
        (aa(ii),irna(ii),ii=iptra(n),ne)
    end if
  end if

! Initializations.
  jlast = 0
  nullc = 0
  if (job>1 .and. icntl(6)>0 .and. icntl(6)<n) jlast = min(np,n-icntl(6))

  u = min(cntl(2),one)
  u = max(u,zero)
  do 20 i = 1, m
    iw(i+n) = 0
    w(i) = zero
20 continue
  mord = 0
  if1 = lfact + 1
  if2 = 0
  nf = n - np
  mf = 0
  il2 = 2
  if (jlast>0) il2 = iptrl(jlast)
  neu = 0

! Jump if JOB is equal to 2.
  if (job==2) go to 370

  if (job==3) then
! Reconstruct IP and set MORD
    do 30 j = 1, np
      ia1 = iptru(j) + 1
      if (ia1>iptrl(j)) go to 30
      if (j<=jlast) then
        mord = mord + 1
        ip(irnf(ia1)) = -j
      else
        ip(irnf(ia1)) = j
      end if
30  continue
    mf = irnf(2)
    ia1 = iptrl(n)
    do 40 j = 1, mf
      ip(irnf(ia1+j)) = np + j
40  continue
  end if

! Store copies of column ends ready for pruning
  do 50 k = 1, jlast
    iw(m+n+k) = iptrl(k)
50 continue

! Each pass through this main loop processes column K.
  do 310 k = jlast + 1, n
    drop = .false.
    if (k==np+1) then
! Set up data structure for full part.
      mf = m - mord
      if1 = lfact + 1 - mf
      ii = 0
      do 60 i = 1, m
        if (ip(i)>0) then
          iw(i+n) = n
          irnf(if1+ii) = i
          ii = ii + 1
          ip(i) = np + ii
        end if
60    continue
      if1 = lfact + 1 - max(mf*nf,mf+nf)
      if2 = if1 - 1 + mf*max(0,jlast-np)
    end if
    j = k
    if (iq(1)>0) j = iq(k)
    ia1 = iptra(j)
    ia2 = ne
    if (j/=n) ia2 = iptra(j+1) - 1
    iu1 = il2 + 1
    iu2 = iu1 - 1
    il1 = if1 - 1 + ia1 - ia2
    il2 = il1 - 1
    info(4) = max(info(4),neu+lfact-il1+iu2+m+1)
    info(9) = max(info(9),neu+lfact-il1+iu2+m+1)
    if (il1-iu2<=m) then
      if (info(1)/=-3) then
        info(1) = -3
        if (icntl(8)/=0) go to 480
! Get rid of U info.
        neu = il2 + lfact + 1 - mf - if1
        if1 = lfact + 1 - mf
        if2 = if1 - 1
        il2 = 0
        eye = 0
        do 80 j = 1, min(k-1,np)
          iu2 = iptru(j)
          iptru(j) = eye
          il2 = iptrl(j)
          neu = neu + iu2 - il2
          do 70 ii = iu2 + 1, il2
            eye = eye + 1
            irnf(eye) = irnf(ii)
            fact(eye) = fact(ii)
70        continue
          iptrl(j) = eye
          iw(m+n+j) = eye
80      continue
        iu1 = eye + 1
        iu2 = eye
        il1 = if1 - 1 + ia1 - ia2
        il2 = il1 - 1
      end if
! Quit if LFACT is much too small
      if (il1-iu2<=m) go to 480
    end if
! Load column K of AA into full vector W and into the back of IRNF.
! Check for duplicates.
    eye = il1
    do 90 ii = ia1, ia2
      i = irna(ii)
!!!
! Can't happen if called from HSL_MA48
!     if (iw(i+n)==-1) go to 540
      iw(i+n) = -1
      w(i) = aa(ii)
      irnf(eye) = i
      eye = eye + 1
90  continue
! Depth first search to find topological order for triangular solve
!     and structure of column K of L/U
! IW(J) is used to hold a pointer to next entry in column J
!     during the depth-first search at stage K, J = 1,..., N.
! IW(I+N) is set to K when row I has been processed, and to N for rows
!     of the full part once column NP has been passed. It is also
!     used for backtracking, a negative value being used to point to the
!     previous row in the chain.
! IW(M+N+I) is set to the position in FACT and IRNF of the end of the
!     active part of the column after pruning.  It is initially set to
!     IPTRL(I) and is flagged negative when column has been pruned.
! Set IPTRL temporarily for column K so that special code is
!     not required to process this column.
    iptrl(k) = eye - 1
    iw(m+n+k) = eye - 1
! IW(K) is set to beginning of original column K.
    iw(k) = il1
    j = k
! The outer loop of the depth-first search is executed once for column
!      K and twice for each entry in the upper-triangular part of column
!      K (once to initiate a search in the corresponding column and
!      once when the search in the column is finished).
    do 120 jdummy = 1, 2*k
! Look through column J of L (or column K of A). All the entries
!     are entries of the filled-in column K. Store new entries of the
!     lower triangle and continue until reaching an entry of the upper
!     triangle.
      do 100 ii = iw(j), abs(iw(m+n+j))
        i = irnf(ii)
! Jump if index I already encountered in column K or is in full part.
        if (iw(i+n)>=k) go to 100
        if (ip(i)<=0) go to 110
! Entry is in lower triangle. Flag it and store it in L.
        iw(i+n) = k
        il1 = il1 - 1
        irnf(il1) = i
100   continue
      if (j==k) go to 130
! Flag J, put its row index into U, and backtrack
      iu2 = iu2 + 1
      i = irnf(iptru(j)+1)
      irnf(iu2) = i
      j = -iw(i+n)
      iw(i+n) = k
      go to 120
! Entry in upper triangle.  Move search to corresponding column.
110   iw(i+n) = -j
      iw(j) = ii + 1
      j = -ip(i)
      iw(j) = iptru(j) + 2
120 continue
! Run through column K of U in the lexicographical order that was just
!     constructed, performing elimination operations.
130 do 150 ii = iu2, iu1, -1
      i = irnf(ii)
      j = -ip(i)
! Add multiple of column J of L to column K
      eye1 = iptru(j) + 1
      if (abs(w(i))<cntl(3)) go to 150
      amult = -w(i)*fact(eye1)
! Note we are storing negative multipliers
      w(i) = amult
      do 140 eye = eye1 + 1, iptrl(j)
        i = irnf(eye)
        w(i) = w(i) + amult*fact(eye)
140   continue
      rinfo(1) = rinfo(1) + one + 2*(iptrl(j)-eye1)
150 continue

! Unload reals of column of U and set pointer
    if (cntl(3)>zero) then
      eye = iu1
      do 160 ii = iu1, iu2
        i = irnf(ii)
        if (abs(w(i))<cntl(3)) then
          info(6) = info(6) + 1
        else
          irnf(eye) = -ip(i)
          fact(eye) = w(i)
          eye = eye + 1
        end if
        w(i) = zero
160   continue
      iu2 = eye - 1
    else
      do 170 ii = iu1, iu2
        i = irnf(ii)
        irnf(ii) = -ip(i)
        fact(ii) = w(i)
        w(i) = zero
170   continue
    end if
    if (info(1)==-3) then
      neu = neu + iu2 - iu1 + 1
      iu2 = iu1 - 1
    end if
    iptru(k) = iu2
    if (k<=np) then
! Find the largest entry in the column and drop any small entries
      maxent = zero
      if (cntl(3)>zero) then
        eye = il1
        do 180 ii = il1, il2
          i = irnf(ii)
          if (abs(w(i))<cntl(3)) then
            info(6) = info(6) + 1
            w(i) = zero
            drop = .true.
          else
            irnf(eye) = i
            eye = eye + 1
            maxent = max(abs(w(i)),maxent)
          end if
180     continue
        il2 = eye - 1
      else
        do 190 ii = il1, il2
          maxent = max(abs(w(irnf(ii))),maxent)
190     continue
      end if
! Unload column of L, performing pivoting and moving indexing
!      information.
      pivlim = u*maxent
      eye = iu2
      iqpiv = m + n
      if (il1>il2) nullc = nullc + 1
      do 200 ii = il1, il2
        i = irnf(ii)
        eye = eye + 1
        irnf(eye) = i
        fact(eye) = w(i)
        w(i) = zero
! Find position of pivot
        if (abs(fact(eye))>=pivlim) then
          if (abs(fact(eye))>cntl(4)) then
            if (ip(i)<iqpiv) then
              iqpiv = ip(i)
              ipiv = eye
            end if
          end if
        end if
200   continue
      il1 = iu2 + 1
      il2 = eye
      if (il1<=il2) then
! Column is not null
        if (iqpiv==m+n) then
! All entries in the column are too small to be pivotal. Drop them all.
          if (cntl(3)>zero) info(6) = info(6) + eye - iu2
          il2 = iu2
        else
          if (il1/=ipiv) then
! Move pivot to front of L
            asw = fact(ipiv)
            fact(ipiv) = fact(il1)
            fact(il1) = asw
            isw = irnf(il1)
            irnf(il1) = irnf(ipiv)
            irnf(ipiv) = isw
          end if
! Reciprocate pivot
          info(5) = info(5) + 1
          fact(il1) = one/fact(il1)
          rinfo(1) = rinfo(1) + one
! Record pivot row
          mord = mord + 1
          ip(irnf(il1)) = -k
        end if
      end if
    else
! Treat column as full
      il2 = iptru(k)
!DIR$ IVDEP
      do 210 ii = lfact - mf + 1, lfact
        i = irnf(ii)
        if2 = if2 + 1
        fact(if2) = w(i)
        w(i) = zero
210   continue
      if (info(1)==-3) if2 = if2 - mf
    end if
    iw(m+n+k) = il2
    iptrl(k) = il2
    if (drop) go to 310
! Scan columns involved in update of column K and remove trailing block.
    do 300 ii = iu1, iu2
      i = irnf(ii)
! Jump if column already pruned.
      if (iw(m+n+i)<0) go to 300
      begcol = iptru(i) + 2
      endcol = iptrl(i)
! Scan column to see if there is an entry in the current pivot row.
      if (k<=np) then
        do 220 jj = begcol, endcol
          if (ip(irnf(jj))==-k) go to 230
220     continue
        go to 300
      end if
! Sort the entries so that those in rows already pivoted (negative IP
!    values) precede the rest.
230   do 280 jdummy = begcol, endcol
        jj = begcol
        do 240 begcol = jj, endcol
          if (ip(irnf(begcol))>0) go to 250
240     continue
        go to 290
250     jj = endcol
        do 260 endcol = jj, begcol, -1
          if (ip(irnf(endcol))<0) go to 270
260     continue
        go to 290
270     asw = fact(begcol)
        fact(begcol) = fact(endcol)
        fact(endcol) = asw
        j = irnf(begcol)
        irnf(begcol) = irnf(endcol)
        irnf(endcol) = j
        begcol = begcol + 1
        endcol = endcol - 1
280   continue
290   iw(m+n+i) = -endcol
300 continue
310 continue
  if (n==np) then
! Set up data structure for the (null) full part.
    mf = m - mord
    if1 = lfact + 1 - mf
    ii = 0
    do 320 i = 1, m
      if (ip(i)>0) then
        iw(i+n) = n
        irnf(if1+ii) = i
        ii = ii + 1
        ip(i) = np + ii
      end if
320 continue
    if1 = lfact + 1 - max(mf*nf,mf+nf)
    if2 = if1 - 1 + mf*max(0,jlast-np)
  end if
  if (info(5)==min(m,n)) then
! Restore sign of IP
    do 330 i = 1, m
      ip(i) = abs(ip(i))
330 continue
  else
! Complete IP
    mord = np
    do 340 i = 1, m
      if (ip(i)<0) then
        ip(i) = -ip(i)
      else
        mord = mord + 1
        ip(i) = mord
      end if
340 continue
  end if
  irnf(1) = info(6)
  irnf(2) = mf
  info(7) = mf
  fact(1) = cntl(3)
  fact(2) = cntl(4)
  if (info(1)==-3) go to 520
! Move full part forward
  if2 = if2 - mf*nf
  do 350 ii = 1, mf*nf
    fact(il2+ii) = fact(if1-1+ii)
350 continue
  do 360 ii = 1, mf
    irnf(il2+ii) = irnf(lfact-mf+ii)
360 continue
  if1 = il2 + 1
  go to 440

! Fast factor (JOB = 2)
! Each pass through this main loop processes column K.
370 mf = irnf(2)
  if1 = iptrl(n) + 1
  if2 = if1 - 1
  do 430 k = jlast + 1, n
    j = k
    if (iq(1)>0) j = iq(k)
    ia1 = iptra(j)
    ia2 = ne
    if (j/=n) ia2 = iptra(j+1) - 1
    iu1 = il2 + 1
    iu2 = iptru(k)
    il1 = iu2 + 1
    il2 = iptrl(k)
! Load column K of A into full vector W
    do 380 ii = ia1, ia2
      w(irna(ii)) = aa(ii)
380 continue
! Run through column K of U in lexicographical order, performing
!      elimination operations.
    do 400 ii = iu2, iu1, -1
      j = irnf(ii)
      i = irnf(iptru(j)+1)
! Add multiple of column J of L to column K
      eye1 = iptru(j) + 1
      amult = -w(i)*fact(eye1)
! Note we are storing negative multipliers
      fact(ii) = amult
      w(i) = zero
      do 390 eye = eye1 + 1, iptrl(j)
        i = irnf(eye)
        w(i) = w(i) + amult*fact(eye)
390   continue
      rinfo(1) = rinfo(1) + one + 2*(iptrl(j)-eye1)
400 continue
    if (k<=np) then
      if (il1<=il2) then
! Load column of L.
!DIR$ IVDEP
        do 410 ii = il1, il2
          i = irnf(ii)
          fact(ii) = w(i)
          w(i) = zero
410     continue
! Test pivot. Note that this is the only numerical test when JOB = 2.
        if (abs(fact(il1))<=cntl(4)) then
          go to 530
        else
! Reciprocate pivot
          info(5) = info(5) + 1
          fact(il1) = one/fact(il1)
          rinfo(1) = rinfo(1) + one
        end if
      end if
    else
! Treat column as full
      do 420 ii = if1, if1 + mf - 1
        i = irnf(ii)
        if2 = if2 + 1
        fact(if2) = w(i)
        w(i) = zero
420   continue
    end if
430 continue
  info(4) = max(if1+mf+nf-1,if2)
  info(9) = max(if1+mf+nf-1,if2)

440 if (mf>0 .and. nf>0) then
! Factorize full block
    if (icntl(5)>1) call ma50gd(mf,nf,fact(if1),mf,icntl(5),cntl(4), &
      irnf(if1+mf),rank)
    if (icntl(5)==1) call ma50fd(mf,nf,fact(if1),mf,cntl(4),irnf(if1+mf),rank)
    if (icntl(5)<=0) call ma50ed(mf,nf,fact(if1),mf,cntl(4),irnf(if1+mf),rank)
    info(5) = info(5) + rank
    do 450 i = 1, min(mf,nf)
      rinfo(1) = rinfo(1) + real(mf - i + 1,kind=wp) + &
                            real(mf-i,kind=wp)*real(nf-i,kind=wp)*2.0d0
450 continue
  end if
  if (info(5)<min(m,n)) info(1) = 1
  if (mp>0 .and. icntl(3)>2) then
    write (mp,'(A,I6,A,F12.1/A,I3,A,4I8)') ' Leaving MA50BD with IRNF(2) =', &
      irnf(2), ' RINFO(1) =', rinfo(1), ' INFO(1) =', info(1), ' INFO(4:7) =', &
      (info(j),j=4,7)
    if (icntl(3)>3) then
      if (job/=2) write (mp,'(A,(T6,10(I7)))') ' IP = ', ip
      do 460 j = 1, n
        if (j>1) then
          if (iptrl(j-1)<iptru(j)) write (mp,'(A,I5,A,(T18,3(1P,E12.4,I5)))') &
            ' Column', j, ' of U', (fact(ii),irnf(ii),ii=iptrl(j-1)+1,iptru(j) &
            )
        end if
        if (iptru(j)<iptrl(j)) write (mp,'(A,I5,A,(T18,3(1P,E12.4,I5)))') &
          ' Column', j, ' of L', (fact(ii),irnf(ii),ii=iptru(j)+1,iptrl(j))
460   continue
      write (mp,'(A)') ' Full part'
      write (mp,'((6I12))') (irnf(if1+mf+j),j=0,nf-1)
      do 470 i = 0, mf - 1
        write (mp,'(I4,1P,6E12.4:/(4X,1P,6E12.4))') irnf(if1+i), &
          (fact(if1+i+j*mf),j=0,nf-1)
470   continue
    end if
  end if
  go to 550

! Error conditions
! LFACT is much too small or ICNTL(8)/=0. Patch up IP and quit.
480 do 490 i = 1, m
    iw(i) = 0
490 continue
  do 500 i = 1, m
    if (ip(i)>0) then
      iw(ip(i)) = i
    else
      ip(i) = -ip(i)
    end if
500 continue
  do 510 i = 1, m
    if (iw(i)>0) then
      ip(iw(i)) = k
      k = k + 1
    end if
510 continue
520 info(1) = -3
  if (lp>0) then
    write (lp,'(/A)') ' **** Error return from MA50BD **** '
    if (icntl(8)==0) then
      write (lp,'(A,I7,A,I7)') ' LFACT must be increased from', lfact, &
        ' to at least', info(4)
    else
      write (lp,'(A,I7)') ' LFACT must be increased from', lfact
    end if
  end if
  go to 550
530 info(1) = -(7+k)
  if (lp>0) write (lp,'(/A/A,I6,A)') ' **** Error return from MA50BD **** ', &
    ' Small pivot found in column', k, ' of the permuted matrix.'
!  go to 550
!!!
! Again this error can only happen if MA50 is tested separately
!540 info(1) = -4
!  if (lp>0) write (lp,'(/A/(3(A,I9)))') ' **** Error return from MA50BD ****',&
!    ' Entry in row', i, ' and column', j, ' duplicated'
550 end subroutine ma50bd

SUBROUTINE ma50cd(m,n,icntl,iq,np,trans,lfact,fact,irnf,iptrl,iptru,b,x,w)
! MA50C/CD uses the factorization produced by
!     MA50B/BD to solve A x = b or (A trans) x = b.

  integer m, n, icntl(20), np
  integer(long) :: iq(*)
  logical trans
  integer(long) :: lfact
  real(wp) :: fact(lfact)
  integer(long) irnf(lfact)
  integer(long) :: iptrl(n), iptru(n)
  real(wp) :: b(*), x(*), w(*)
!!!
! integer(long) :: info(15)

! M  is an integer variable set to the number of rows.
!     It is not altered by the subroutine.
! N  is an integer variable set to the number of columns.
!     It is not altered by the subroutine.
! ICNTL must be set by the user as follows and is not altered.
!     ICNTL(1)  must be set to the stream number for error messages.
!       A value less than 1 suppresses output.
!     ICNTL(2) must be set to the stream number for diagnostic output.
!       A value less than 1 suppresses output.
!     ICNTL(3) must be set to control the amount of output:
!       0 None.
!       1 Error messages only.
!       2 Error and warning messages.
!       3 As 2, plus scalar parameters and a few entries of array
!         parameters on entry and exit.
!       4 As 2, plus all parameters on entry and exit.
!     ICNTL(5) must be set to control the level of BLAS used:
!       0 Level 1 BLAS.
!      >0 Level 2 BLAS.
! IQ is an integer array holding the permutation Q.
!     It is not altered by the subroutine.
! NP is an integer variable that must be unchanged since calling
!     MA50B/BD. It holds the number of rows and columns in packed
!     storage. It is not altered by the subroutine.
! TRANS a logical variable thatmust be set to .TRUE. if (A trans)x = b
!     is to be solved and to .FALSE. if A x = b is to be solved.
!     TRANS is not altered by the subroutine.
! LFACT is an integer variable set to the size of FACT and IRNF.
!     It is not altered by the subroutine.
! FACT is an array that must be unchanged since calling MA50B/BD. It
!     holds the packed part of L/U by columns, and the full part of L/U
!     by columns. U has unit diagonal entries, which are not stored, and
!     the signs of the off-diagonal entries are inverted.  In the packed
!     part, the entries of U precede the entries of L; also the diagonal
!     entries of L head each column of L and are reciprocated.
!     FACT is not altered by the subroutine.
! IRNF is an integer array that must be unchanged since calling
!     MA50B/BD. It holds the row numbers of the packed part of L/U, and
!     the row numbers of the full part of L/U.
!     It is not altered by the subroutine.
! IPTRL is an integer array that must be unchanged since calling
!     MA50B/BD. For J = 1,..., NP, IPTRL(J) holds the position in
!     FACT and IRNF of the end of column J of L.
!     For J = NP+1,..., N, IPTRL(J) is equal to IPTRU(J).
!     It is not altered by the subroutine.
! IPTRU is an integer array that must be unchanged since calling
!     MA50B/BD. For J = 1,..., N, IPTRU(J) holds the position in
!     FACT and IRNF of the end of the packed part of column J of U.
!     It is not altered by the subroutine.
! B is an array that must be set to the vector b.
!     It is not altered.
! X is an array that need not be set on entry. On return, it holds the
!    solution x.
! W is a work array of length max(M,N).
! INFO need not be set on entry. On exit, it holds the following:
!    INFO(1) A nonzero value will indicate an error return. Possible
!      nonzero values are:
!      -1  M < 1 or N < 1

  real(wp) :: zero
  parameter (zero=0.0_wp)
  integer(long) :: ii, ia1, if1
  integer :: i, j, mf, mp, nf
  real(wp) :: prod
! I Temporary variable holding row number.
! II Position of the current entry in IRNF.
! IA1 Position of the start of the current row or column.
! IF1 Position of the start of the full part of U.
! J Temporary variable holding column number.
! LP Unit for error messages.
! MF Number of rows held in full format.
! MP Unit for diagnostic messages.
! NF Number of columns held in full format.
! PROD Temporary variable used to accumulate inner products.

!!!
! Now lp is set but never used
! lp = icntl(1)
  mp = icntl(2)
! if (icntl(3)<=0) lp = 0
  if (icntl(3)<=1) mp = 0

! Make some simple checks
!!!
! Can't happen if called from HSL_MA48
! if (m<1 .or. n<1) go to 250

  if (mp>0 .and. icntl(3)>2) write (mp,'(/2(A,I6),A,I4,A,L2)') &
    ' Entering MA50CD with M=', m, ' N =', n, ' NP =', np, ' TRANS =', trans
  if1 = iptrl(n) + 1
  mf = irnf(2)
  nf = n - np
  if (mp>0 .and. icntl(3)>2) write (mp,'(A,I5,A,I5)') &
    ' Size of full submatrix', mf, ' by', nf
  if (mp>0 .and. icntl(3)>3) then
    do 10 j = 1, n
      if (j>1) then
        if (iptrl(j-1)<iptru(j)) write (mp,'(A,I5,A,(T18,3(1P,E12.4,I5)))') &
          ' Column', j, ' of U', (fact(ii),irnf(ii),ii=iptrl(j-1)+1,iptru(j))
      end if
      if (iptru(j)<iptrl(j)) write (mp,'(A,I5,A,(T18,3(1P,E12.4,I5)))') &
        ' Column', j, ' of L', (fact(ii),irnf(ii),ii=iptru(j)+1,iptrl(j))
10  continue
    write (mp,'(A)') ' Full part'
    write (mp,'((6I12))') (irnf(if1+mf+j),j=0,nf-1)
    do 20 i = 0, mf - 1
      write (mp,'(I4,1P,6E12.4:/(4X,1P,6E12.4))') irnf(if1+i), &
        (fact(if1+i+j*mf),j=0,nf-1)
20  continue
  end if

  if (trans) then
    if (mp>0 .and. icntl(3)>3) write (mp,'(A4,5F10.4:/(4X,5F10.4))') ' B =', &
      (b(i),i=1,n)
    if (iq(1)>0) then
      do 30 i = 1, n
        w(i) = b(iq(i))
30    continue
    else
      do 40 i = 1, n
        w(i) = b(i)
40    continue
    end if
    do 50 i = 1, m
      x(i) = zero
50  continue
! Forward substitution through packed part of (U trans).
    do 70 i = 2, n
      prod = zero
      do 60 ii = iptrl(i-1) + 1, iptru(i)
        prod = prod + fact(ii)*w(irnf(ii))
60    continue
      w(i) = w(i) + prod
70  continue
! Backsubstitute through the full part of (PL) trans.
    do 80 i = 1, nf
      x(i) = w(np+i)
80  continue
    if (mf>0 .and. nf>0) then
      call ma50hd(trans,mf,nf,fact(if1),mf,irnf(if1+mf),x,icntl(5))
    else
      do 90 i = 1, mf
        x(i) = zero
90    continue
    end if
    do 100 i = mf, 1, -1
      j = irnf(if1+i-1)
      if (j/=i) x(j) = x(i)
100 continue
! Backsubstitute through the packed part of (PL) trans.
    do 120 i = np, 1, -1
      ia1 = iptru(i) + 1
      if (ia1>iptrl(i)) go to 120
      prod = zero
      do 110 ii = ia1 + 1, iptrl(i)
        prod = prod + fact(ii)*x(irnf(ii))
110   continue
      x(irnf(ia1)) = (w(i)-prod)*fact(ia1)
120 continue
    if (mp>0 .and. icntl(3)>3) write (mp,'(A/(4X,5F10.4))') &
      ' Leaving MA50CD with X =', (x(i),i=1,m)

  else
    if (mp>0 .and. icntl(3)>3) write (mp,'(A4,5F10.4:/(4X,5F10.4))') ' B =', &
      (b(i),i=1,m)
! Forward substitution through the packed part of PL
    do 130 i = 1, m
      w(i) = b(i)
130 continue
    do 150 i = 1, np
      ia1 = iptru(i) + 1
      if (ia1<=iptrl(i)) then
        x(i) = w(irnf(ia1))*fact(ia1)
        if (x(i)/=zero) then
!DIR$ IVDEP
          do 140 ii = ia1 + 1, iptrl(i)
            w(irnf(ii)) = w(irnf(ii)) - fact(ii)*x(i)
140       continue
        end if
      end if
150 continue
! Forward substitution through the full part of PL
    if (mf>0 .and. nf>0) then
      do 160 i = 1, mf
        w(i) = w(irnf(if1+i-1))
160   continue
      call ma50hd(trans,mf,nf,fact(if1),mf,irnf(if1+mf),w,icntl(5))
      do 170 i = 1, nf
        x(np+i) = w(i)
170   continue
    else
      do 180 i = 1, nf
        x(np+i) = zero
180   continue
    end if
! Back substitution through the packed part of U
    do 200 j = n, max(2,np+1), -1
      prod = x(j)
!DIR$ IVDEP
      do 190 ii = iptrl(j-1) + 1, iptru(j)
        x(irnf(ii)) = x(irnf(ii)) + fact(ii)*prod
190   continue
200 continue
    do 220 j = np, 2, -1
      ia1 = iptru(j)
      if (ia1>=iptrl(j)) then
        x(j) = zero
      else
        prod = x(j)
!DIR$ IVDEP
        do 210 ii = iptrl(j-1) + 1, ia1
          x(irnf(ii)) = x(irnf(ii)) + fact(ii)*prod
210     continue
      end if
220 continue
    if (np>=1 .and. iptru(1)>=iptrl(1)) x(1) = zero
    if (iq(1)>0) then
!         Permute X
      do 230 i = 1, n
        w(i) = x(i)
230   continue
      do 240 i = 1, n
        x(iq(i)) = w(i)
240   continue
    end if
    if (mp>0 .and. icntl(3)>3) write (mp,'(A/(4X,5F10.4))') &
      ' Leaving MA50CD with X =', (x(i),i=1,n)
  end if
  return
!!!
! Error condition.  Can only happen if MA50 is called directly.
!250 info(1) = -1
!  if (lp>0) write (lp,'(/A/2(A,I8))') ' **** Error return from MA50CD ****', &
!    ' M =', m, ' N =', n
end subroutine ma50cd

SUBROUTINE ma50dd(a,ind,iptr,n,disp,reals)
! This subroutine performs garbage collection on the arrays A and IND.
! DISP is the position in arrays A/IND immediately after the data
!     to be compressed.
!     On exit, DISP equals the position of the first entry
!     after the compressed part of A/IND.

  integer n
  integer(long) :: disp
  real(wp) :: a(:)
  integer(long) :: iptr(n)
  logical reals
  integer :: ind(:)
! Local variables.
  integer j
  integer(long) :: k, kn
! Set the first entry in each row(column) to the negative of the
!     row(column) and hold the column(row) index in the row(column)
!     pointer.  This enables the start of each row(column) to be
!     recognized in a subsequent scan.
  do 10 j = 1, n
    k = iptr(j)
    if (k>0) then
      iptr(j) = ind(k)
      ind(k) = -j
    end if
10 continue
  kn = 0
! Go through arrays compressing to the front so that there are no
!     zeros held in positions 1 to DISP-1 of IND.
!     Reset first entry of each row(column) and the pointer array IPTR.
  do 20 k = 1, disp - 1
    if (ind(k)==0) go to 20
    kn = kn + 1
    if (reals) a(kn) = a(k)
    if (ind(k)<=0) then
! First entry of row(column) has been located.
      j = -ind(k)
      ind(k) = iptr(j)
      iptr(j) = kn
    end if
    ind(kn) = ind(k)
20 continue
  disp = kn + 1
end subroutine ma50dd


SUBROUTINE ma50ed(m,n,a,lda,pivtol,ipiv,rank)
!*
  integer m, n, rank
  integer :: lda
  real(wp) :: pivtol

  integer(long) ipiv(n)
  real(wp) :: a(lda,n)


!  Purpose
!  =======

!  MA50ED computes an LU factorization of a general m-by-n matrix A.

!  The factorization has the form
!     A = P * L * U * Q
!  where P is a permutation matrix of order m, L is lower triangular
!  of order m with unit diagonal elements, U is upper trapezoidal of
!  order m * n, and Q is a permutation matrix of order n.

!  Row interchanges are used to ensure that the entries of L do not
!  exceed 1 in absolute value. Column interchanges are used to
!  ensure that the first r diagonal entries of U exceed PIVTOL in
!  absolute value. If r < m, the last (m-r) rows of U are zero.

!  This is the Level 1 BLAS version.

!  Arguments
!  =========

!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 1.

!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 1.

!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U*Q; the unit diagonal elements of L are not stored.

!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).

!  PIVTOL  (input) DOUBLE PRECISION
!          The pivot tolerance. Any entry with absolute value less
!          than or equal to PIVTOL is regarded as unsuitable to be a
!          pivot.
!*
!  IPIV    (output) INTEGER array, dimension (N)
!          The permutations; for 1 <= i <= RANK, row i of the
!          matrix was interchanged with row IPIV(i); for
!          RANK + 1 <= j <= N, column j of the
!          matrix was interchanged with column -IPIV(j).

!  RANK    (output) INTEGER
!          The computed rank of the matrix.

!  =====================================================================

  real(wp) :: one, zero
  parameter (one=1.0_wp,zero=0.0_wp)

  integer i, j, jp, k
  logical pivot
! I   Row index.
! J   Current column.
! JP  Pivot position.
! K   Main loop index.
! PIVOT True if there is a pivot in current column.

  integer idamax
  external idamax

  external daxpy, dscal, dswap
  intrinsic abs


  j = 1
  do 30 k = 1, n

!        Update elements in column J.
    do 10 i = 1, j - 1
      if (m>i) call daxpy(m-i,-a(i,j),a(i+1,i),1,a(i+1,j),1)
10  continue

!        Find pivot.
    if (j<=m) then
      jp = j - 1 + idamax(m-j+1,a(j,j),1)
      ipiv(j) = jp
      pivot = abs(a(jp,j)) > pivtol
    else
      pivot = .false.
    end if
    if (pivot) then

!           Apply row interchange to columns 1:N+J-K.
      if (jp/=j) call dswap(n+j-k,a(j,1),lda,a(jp,1),lda)

!           Compute elements J+1:M of J-th column.
      if (j<m) call dscal(m-j,one/a(j,j),a(j+1,j),1)

!           Update J
      j = j + 1

    else

      do 20 i = j, m
        a(i,j) = zero
20    continue
!           Apply column interchange and record it.
      if (k<n) call dswap(m,a(1,j),1,a(1,n-k+j),1)
      ipiv(n-k+j) = -j

    end if

30 continue

  rank = j - 1

!     End of MA50ED

end subroutine ma50ed


SUBROUTINE ma50fd(m,n,a,lda,pivtol,ipiv,rank)

!  -- This is a variant of the LAPACK routine DGETF2 --

  integer m, n, rank
  integer :: lda
  real(wp) :: pivtol

  integer(long) ipiv(n)
  real(wp) :: a(lda,n)


!  Purpose
!  =======

!  MA50FD computes an LU factorization of a general m-by-n matrix A.

!  The factorization has the form
!     A = P * L * U * Q
!  where P is a permutation matrix of order m, L is lower triangular
!  of order m with unit diagonal elements, U is upper trapezoidal of
!  order m * n, and Q is a permutation matrix of order n.

!  Row interchanges are used to ensure that the entries of L do not
!  exceed 1 in absolute value. Column interchanges are used to
!  ensure that the first r diagonal entries of U exceed PIVTOL in
!  absolute value. If r < m, the last (m-r) rows of U are zero.

!  This is the Level 2 BLAS version.

!  Arguments
!  =========

!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 1.

!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 1.

!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U*Q; the unit diagonal elements of L are not stored.

!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).

!  PIVTOL  (input) DOUBLE PRECISION
!          The pivot tolerance. Any entry with absolute value less
!          than or equal to PIVTOL is regarded as unsuitable to be a
!          pivot.
!*
!  IPIV    (output) INTEGER array, dimension (N)
!          The permutations; for 1 <= i <= RANK, row i of the
!          matrix was interchanged with row IPIV(i); for
!          RANK + 1 <= j <= N, column j of the
!          matrix was interchanged with column -IPIV(j).

!  RANK    (output) INTEGER
!          The computed rank of the matrix.

!  =====================================================================

  real(wp) :: one, zero
  parameter (one=1.0_wp,zero=0.0_wp)

  integer i, j, jp, k
  logical pivot
! I   Row index.
! J   Current column.
! JP  Pivot position.
! K   Main loop index.
! PIVOT True if there is a pivot in current column.

  integer idamax
  external idamax

  external dgemv, dscal, dswap

  intrinsic abs


  j = 1
  do 20 k = 1, n

    if (j<=m) then
!           Update diagonal and subdiagonal elements in column J.
      call dgemv('No transpose',m-j+1,j-1,-one,a(j,1),lda,a(1,j),1,one,a(j,j), &
        1)
!          Find pivot.
      jp = j - 1 + idamax(m-j+1,a(j,j),1)
      ipiv(j) = jp
      pivot = abs(a(jp,j)) > pivtol
    else
      pivot = .false.
    end if

    if (pivot) then

!           Apply row interchange to columns 1:N+J-K.
      if (jp/=j) call dswap(n+j-k,a(j,1),lda,a(jp,1),lda)

!           Compute elements J+1:M of J-th column.
      if (j<m) call dscal(m-j,one/a(j,j),a(j+1,j),1)

      if (j<n) then
!             Compute block row of U.
        call dgemv('Transpose',j-1,n-j,-one,a(1,j+1),lda,a(j,1),lda,one, &
          a(j,j+1),lda)
      end if

!           Update J
      j = j + 1

    else

      do 10 i = j, m
        a(i,j) = zero
10    continue
!           Apply column interchange and record it.
      if (k<n) call dswap(m,a(1,j),1,a(1,n-k+j),1)
      ipiv(n-k+j) = -j

    end if

20 continue

  rank = j - 1

!     End of MA50FD

end subroutine ma50fd


SUBROUTINE ma50gd(m,n,a,lda,nb,pivtol,ipiv,rank)

!  -- This is a variant of the LAPACK routine DGETRF --

  integer m, n, nb, rank
  integer :: lda
  real(wp) :: pivtol

  integer(long) ipiv(n)
  real(wp) :: a(lda,n)


!  Purpose
!  =======

!  MA50GD computes an LU factorization of a general m-by-n matrix A.

!  The factorization has the form
!     A = P * L * U * Q
!  where P is a permutation matrix of order m, L is lower triangular
!  of order m with unit diagonal elements, U is upper trapezoidal of
!  order m * n, and Q is a permutation matrix of order n.

!  Row interchanges are used to ensure that the entries of L do not
!  exceed 1 in absolute value. Column interchanges are used to
!  ensure that the first r diagonal entries of U exceed PIVTOL in
!  absolute value. If r < m, the last (m-r) rows of U are zero.

!  This is the Level 3 BLAS version.

!  Arguments
!  =========

!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 1.

!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 1.

!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U*Q; the unit diagonal elements of L are not stored.

!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).

!  NB      (input) INTEGER
!          The block size for BLAS3 processing.

!  PIVTOL  (input) DOUBLE PRECISION
!          The pivot tolerance. Any entry with absolute value less
!          than or equal to PIVTOL is regarded as unsuitable to be a
!          pivot.
!*
!  IPIV    (output) INTEGER array, dimension (N)
!          The permutations; for 1 <= i <= RANK, row i of the
!          matrix was interchanged with row IPIV(i); for
!          RANK + 1 <= j <= N, column j of the
!          matrix was interchanged with column -IPIV(j).

!  RANK    (output) INTEGER
!          The computed rank of the matrix.

!  =====================================================================

  real(wp) :: one, zero
  parameter (one=1.0_wp,zero=0.0_wp)

  integer i, j, jj, jp, j1, j2, k
  logical pivot
  real(wp) :: temp

! I   DO index for applying permutations.
! J   Current column.
! JJ  Column in which swaps occur.
! JP  Pivot position.
! J1  Column at start of current block.
! J2  Column at end of current block.
! K   Main loop index.
! PIVOT True if there is a pivot in current column.
! TEMP Temporary variable for swaps.

  external dgemm, dgemv, dswap, dscal, dtrsm, dtrsv

  integer idamax
  external idamax

  intrinsic abs, min


  j = 1
  j1 = 1
  j2 = min(n,nb)
  do 70 k = 1, n

    if (j<=m) then

!          Update diagonal and subdiagonal elements in column J.
      call dgemv('No transpose',m-j+1,j-j1,-one,a(j,j1),lda,a(j1,j),1,one, &
        a(j,j),1)

!          Find pivot.
      jp = j - 1 + idamax(m-j+1,a(j,j),1)
      ipiv(j) = jp
      pivot = abs(a(jp,j)) > pivtol
    else
      pivot = .false.
    end if

    if (pivot) then

!           Apply row interchange to columns J1:J2
      if (jp/=j) call dswap(j2-j1+1,a(j,j1),lda,a(jp,j1),lda)

!           Compute elements J+1:M of J-th column.
      if (j<m) call dscal(m-j,one/a(j,j),a(j+1,j),1)

      if (j+1<=j2) then
!             Compute row of U within current block
        call dgemv('Transpose',j-j1,j2-j,-one,a(j1,j+1),lda,a(j,j1),lda,one, &
          a(j,j+1),lda)
      end if

!           Update J
      j = j + 1

    else

      do 10 i = j, m
        a(i,j) = zero
10    continue

!           Record column interchange and revise J2 if necessary
      ipiv(n-k+j) = -j
!           Apply column interchange.
      if (k/=n) call dswap(m,a(1,j),1,a(1,n-k+j),1)
      if (n-k+j>j2) then
!              Apply operations to new column.
        do 20 i = j1, j - 1
          jp = ipiv(i)
          temp = a(i,j)
          a(i,j) = a(jp,j)
          a(jp,j) = temp
20      continue
        if (j>j1) call dtrsv('Lower','No transpose','Unit',j-j1,a(j1,j1),lda, &
          a(j1,j),1)
      else
        j2 = j2 - 1
      end if

    end if

    if (j>j2) then
!           Apply permutations to columns outside the block
      do 40 jj = 1, j1 - 1
        do 30 i = j1, j2
          jp = ipiv(i)
          temp = a(i,jj)
          a(i,jj) = a(jp,jj)
          a(jp,jj) = temp
30      continue
40    continue
      do 60 jj = j2 + 1, n - k + j - 1
        do 50 i = j1, j2
          jp = ipiv(i)
          temp = a(i,jj)
          a(i,jj) = a(jp,jj)
          a(jp,jj) = temp
50      continue
60    continue

      if (k/=n) then
!              Update the Schur complement
        call dtrsm('Left','Lower','No transpose','Unit',j2-j1+1,n-k,one, &
          a(j1,j1),lda,a(j1,j2+1),lda)
        if (m>j2) call dgemm('No transpose','No transpose',m-j2,n-k,j2-j1+1, &
          -one,a(j2+1,j1),lda,a(j1,j2+1),lda,one,a(j2+1,j2+1),lda)
      end if

      j1 = j2 + 1
      j2 = min(j2+nb,n-k+j-1)

    end if


70 continue
  rank = j - 1

!     End of MA50GD

end subroutine ma50gd

SUBROUTINE ma50hd(trans,m,n,a,lda,ipiv,b,icntl5)

!  -- This is a variant of the LAPACK routine DGETRS --
!     It handles the singular or rectangular case.

  logical trans
  integer m, n, icntl5
  integer :: lda
  integer(long) ipiv(n)
  real(wp) :: a(lda,n), b(*)


!  Purpose
!  =======

!  Solve a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general m by n matrix A using the LU factorization computed
!  by MA50ED, MA50FD, or MA50GD.

!  Arguments
!  =========

!  TRANS   (input) LOGICAL
!          Specifies the form of the system of equations.
!          = .FALSE. :  A * X = B  (No transpose)
!          = .TRUE.  :  A'* X = B  (Transpose)

!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 1.

!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 1.

!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by MA50GD.

!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).

!  IPIV    (input) INTEGER array, dimension (N)
!          The permutations; for 1 <= i <= RANK, row i of the
!          matrix was interchanged with row IPIV(i); for
!          RANK + 1 <= j <= N, column j of the
!          matrix was interchanged with column -IPIV(j).

!  B       (input/output) DOUBLE PRECISION array, size max(M,N)
!          On entry, the right hand side vectors B for the system of
!          linear equations.
!          On exit, the solution vectors, X.

!  ICNTL5  (input) INTEGER
!          0 for BLAS1 or >0 for BLAS2

!  =====================================================================

  integer i, k, rank
! I    Temporary variable.
! K    Temporary variable.
! RANK Rank of matrix.

  real(wp) :: zero
  parameter (zero=0.0_wp)
  real(wp) :: temp
  intrinsic min
  external daxpy, ddot, dtrsv
  real(wp) :: ddot

!   Find the rank
  rank = 0
  do 10 rank = min(m,n), 1, -1
    if (ipiv(rank)>0) go to 20
10 continue

20 if ( .not. trans) then

!        Solve A * X = B.

!        Apply row interchanges to the right hand side.
    do 30 i = 1, rank
      k = ipiv(i)
      temp = b(i)
      b(i) = b(k)
      b(k) = temp
30  continue

!        Solve L*X = B, overwriting B with X.
    if (icntl5>0) then
      if (rank>0) call dtrsv('L','NoTrans','Unit',rank,a,lda,b,1)
    else
      do 40 k = 1, rank - 1
        if (b(k)/=zero) call daxpy(rank-k,-b(k),a(k+1,k),1,b(k+1),1)
40    continue
    end if

!        Solve U*X = B, overwriting B with X.
    if (icntl5>0) then
      if (rank>0) call dtrsv('U','NoTrans','NonUnit',rank,a,lda,b,1)
    else
      do 50 k = rank, 2, -1
        if (b(k)/=zero) then
          b(k) = b(k)/a(k,k)
          call daxpy(k-1,-b(k),a(1,k),1,b(1),1)
        end if
50    continue
      if (rank>0) b(1) = b(1)/a(1,1)
    end if

!        Set singular part to zero
    do 60 k = rank + 1, n
      b(k) = zero
60  continue

!        Apply column interchanges to the right hand side.
    do 70 i = rank + 1, n
      k = -ipiv(i)
      temp = b(i)
      b(i) = b(k)
      b(k) = temp
70  continue

  else

!        Solve A' * X = B.

!        Apply column interchanges to the right hand side.
    do 80 i = n, rank + 1, -1
      k = -ipiv(i)
      temp = b(i)
      b(i) = b(k)
      b(k) = temp
80  continue

!        Solve U'*X = B, overwriting B with X.

    if (icntl5>0) then
      if (rank>0) call dtrsv('U','Trans','NonUnit',rank,a,lda,b,1)
    else
      if (rank>0) b(1) = b(1)/a(1,1)
      do 90 i = 2, rank
        temp = b(i) - ddot(i-1,a(1,i),1,b(1),1)
        b(i) = temp/a(i,i)
90    continue
    end if

!        Solve L'*X = B, overwriting B with X.
    if (icntl5>0) then
      if (rank>0) call dtrsv('L','Trans','Unit',rank,a,lda,b,1)
    else
      do 100 i = rank - 1, 1, -1
        b(i) = b(i) - ddot(rank-i,a(i+1,i),1,b(i+1),1)
100   continue
    end if

!        Set singular part to zero
    do 110 i = rank + 1, m
      b(i) = zero
110 continue

!        Apply row interchanges to the solution vectors.
    do 120 i = rank, 1, -1
      k = ipiv(i)
      temp = b(i)
      b(i) = b(k)
      b(k) = temp
120 continue
  end if

end subroutine ma50hd

SUBROUTINE ma50id(cntl,icntl)
! Set default values for the control arrays.

  real(wp) :: cntl(10)
  integer i, icntl(20)

  cntl(1) = 0.5_wp
  cntl(2) = 0.1_wp
  do 10 i = 3, 10
    cntl(i) = 0.0_wp
10 continue

  icntl(1) = 6
  icntl(2) = 6
  icntl(3) = 1
  icntl(4) = 3
  icntl(5) = 32
  do 20 i = 6, 20
    icntl(i) = 0
20 continue

end subroutine ma50id

end module hsl_ma48_ma50_internal_double

module hsl_ma48_ma48_internal_double
   use hsl_ma48_ma50_internal_double
   implicit none

! This module is based on ma48 2.1.0 (13th March 2007)
! 3 August 2010 Version 3.0.0.  Radical change to include some 
!      INTEGER*64 integers (so that restriction on number of entries is
!      relaxed.  Also use of HSL_ZB01 to allow restarting with more
!      space allocated.

   private
   public :: ma48ad, ma48bd, ma48cd
   public :: ma48id ! public for unit tests

   integer, parameter :: wp = kind(0.0d0)
   integer, parameter :: long = selected_int_kind(18) ! Long integer
   integer(long), parameter :: onel = 1
   integer :: stat

contains

    SUBROUTINE ma48ad(m,n,ne,job,la,a,irn,jcn,keep,cntl,icntl,info,rinfo, &
                      endcol)
! Given a sparse matrix, find its block upper triangular form, choose a
!     pivot sequence for each diagonal block, and prepare data
!     structures for actual factorization.

!     .. Arguments ..
      INTEGER m, n, job
      integer(long) :: ne, la
      real(wp) :: a(2*ne)
      real(wp), allocatable :: fact(:)
      INTEGER, allocatable ::  ifact(:), jfact(:)
      INTEGER(long) ::  irn(2*ne), jcn(ne)
      integer(long) :: keep(*)
      DOUBLE PRECISION cntl(10)
      DOUBLE PRECISION rinfo(10)
      integer, intent(in), optional :: endcol(n)
      INTEGER icntl(20)
      integer(long) :: info(20)
      integer, allocatable :: iw(:)
      real(wp) :: multiplier


! M must be set by the user to the number of rows.
!      It is not altered by the subroutine.  Restriction:  M > 0.
! N must be set by the user to the number of columns.
!      It is not altered by the subroutine.  Restriction:  N > 0.
! NE must be set by the user to the number of entries in the input
!      matrix. It is not altered by the subroutine. Restriction: NE > 0.
! JOB must be set by the user to 1 for automatic choice of pivots and
!      to 2 if the pivot sequence is specified in KEEP.  If JOB is set
!      to 3, then pivots are chosen automatically from the diagonal as
!      long as this is numerically feasible.
!      It is not altered by the subroutine.  Restriction: 1 <= JOB <= 3.
! LA must be set by the user to the size of A, IRN, and JCN.
!      It is not altered by the subroutine. Restriction LA >= 2*NE.
!      Normally a value of 3*NE will suffice.
! A must have its first NE elements set by the user to hold the matrix
!      entries. They may be in any order. If there is more than one for
!      a particular matrix position, they are accumulated. The first
!      NE entries are not altered by the subroutine. The rest is used
!      as workspace for the active submatrix.
! IRN  is an integer array. Entries 1 to NE must be set to the
!      row indices of the corresponding entries in A.  On return, the
!      leading part holds the row numbers of the permuted matrix, with
!      duplicates excluded. The permuted matrix is block upper
!      triangular with recommended pivots on its diagonal. The entries
!      of the block diagonal part are held contiguously from IRN(1)
!      and the entries of the off-diagonal part are held contiguously
!      backwards from IRN(NE) during the computation.  At the end the
!      off-diagonal blocks are held contiguously immediately after the
!      diagonal blocks.
! JCN  is an integer array. Entries 1 to NE must be set to the column
!      indices of the corresponding entries in A. On return, JCN(k)
!      holds the position in IRN that corresponds to the entry that was
!      input in A(k), k=1,NE.  Entries corresponding to out-of-range
!      indices are first set to 0 then, before exit, to NE+1.
!      If duplicates are found, the first entry of JCN is negated.
! KEEP is an integer array that need not be set on a JOB=1 or JOB=3
!      entry. On entry with JOB=2 and always
!      on a successful exit, KEEP(i) holds the position of row i in the
!      permuted matrix, I=1,M and KEEP(M+j) holds the index of the
!      column that is in position j of the permuted matrix, j=1,N.
!      The rest of the array need not be set on entry. On exit:
!        KEEP(IPTRD+j), IPTRD=M+3*N, holds the position in IRN of the
!          start of the block diagonal part of column j, J=1,N;
!        KEEP(IPTRD+N+1) holds the position that immediately follows
!          the end of the block diagonal part of column N;
!        KEEP(IPTRO+j),IPTRO=IPTRD+N+1, holds the position in IRN of
!          the start of the block off-diagonal part of column j, j=1,N;
!        KEEP(IPTRO+N+1) holds the position that immediately
!          follows the end of the block off-diagonal part of column N;
!          During the computation, the columns of the off-diagonal
!          blocks are held in reverse order and KEEP(IPTRO+N+1) points
!          to the position immediately before the off-diagonal block
!          storage.
!        KEEP(KBLOCK+3), KBLOCK=IPTRO+N+1, holds the number of
!          blocks NB in the block triangular form;
!        KEEP(NBLOCK+3*k), NBLOCK=IPTRO+N-1, holds the number of columns
!          in block k, k=1,NB; and
!        KEEP(MBLOCK+3*k), MBLOCK=IPTRO+N, is negative if block k
!          is triangular or holds the number of rows held in packed
!          storage when processing block k, k=1,NB.
!        KEEP(LBLOCK+k), LBLOCK=KBLOCK+3*NB is accessed only if ICNTL(8)
!          is not equal to 0 and will be set to the number of columns
!          in block k which are to be pivoted on last, k=1,NB.
! CNTL  is a real array of length 10 that must be set by the user
!       as follows and is not altered.
!     CNTL(1)  If this is set to a value less than or equal to one, full
!       matrix processing will be used by MA50A/AD when the density of
!       the reduced matrix reaches CNTL(1).
!     CNTL(2) determines the balance used by MA50A/AD and MA50B/BD
!       between pivoting for sparsity and for stability, values near
!       zero emphasizing sparsity and values near one emphasizing
!       stability.
!     CNTL(3) If this is set to a positive value, any entry whose
!       modulus is less than CNTL(3) will be dropped from the factors
!       calculated by MA50A/AD. The factorization will then require
!       less storage but will be inaccurate.
!     CNTL(4)  If this is set to a positive value, any entry whose
!       modulus is less than CNTL(4) will be regarded as zero from
!       the point of view of rank.
!     CNTL(5:10) are not accessed by this subroutine.
! ICNTL is an integer array of length 20 that must be set by the user
!       as follows and is not altered.
!     ICNTL(1)  must be set to the stream number for error messages.
!       A value less than 1 suppresses output.
!     ICNTL(2) must be set to the stream number for diagnostic output.
!       A value less than 1 suppresses output.
!     ICNTL(3) must be set to control the amount of output:
!       0 None.
!       1 Error messages only.
!       2 Error and warning messages.
!       3 As 2, plus scalar parameters and a few entries of array
!         parameters on entry and exit.
!       4 As 2, plus all parameters on entry and exit.
!     ICNTL(4)  If set to a positive value, the pivot search by MA50A/AD
!       is limited to ICNTL(4) columns. This may result in different
!       fill-in and execution time but could give faster execution.
!     ICNTL(5) The block size to be used for full-matrix processing.
!     ICNTL(6) is the minimum size for a block of the block triangular
!       form.
!     ICNTL(7) If not equal to 0, abort when structurally rank deficient
!       matrix found.
!     ICNTL(8) If set to a value other than zero and JOB = 1 or 3,
!       columns with IW flagged 0 are placed at
!       the end of their respective blocks in the first factorization
!       and the remaining columns are assumed unchanged in subsequent
!       factorizations. If set to a value other than zero and JOB = 2,
!       columns before the first 0 entry of IW are assumed unchanged in
!       subsequent factorizations.
!     ICNTL(9:20) are not accessed by this subroutine.
! INFO is an integer array of length 20 that need not be set on entry.
!   INFO(1)  On exit, a negative value  will indicate an error return
!      and a positive value will indicate a warning.
!      Possible nonzero values are:
!      -1  M < 1 or N < 1
!      -2  NE < 1
!      -3  Insufficient space
!      -4  ICNTL(7) not equal to 0 and the matrix is structurally rank
!          deficient.
!      -5  Faulty permutation input on JOB = 2 entry
!      -6  JOB out of range
!     -10  Failure in allocation
!      +1  Row or column number out of range (such entries are ignored)
!          and/or more than one entry for the same matrix position.
!      +2  Matrix rank deficient.  Estimated rank in INFO(5).
!          more than one entry for the same matrix position
!      +3  Combination of warnings +1 and +2.
!      +4  Premature switch to full code because of problems in choosing
!          diagonal pivots (JOB = 3 entry).
!      +5  Combination of warnings +1 and +4.
!      +6  Combination of warnings +2 and +4.
!      +7  Combination of warnings +1, +2 and +4.
!    INFO(2) Number of compresses of the files.
!    INFO(3) Minimum storage required to analyse matrix.
!    INFO(4) Minimum storage required to factorize matrix.
!    INFO(5) Upper bound on the rank of the matrix.
!    INFO(6) Number of entries dropped from the data structure.
!    INFO(7) Order of the largest nontriangular block on the diagonal
!            of the block triangular form.
!    INFO(8) Total of the orders of all the nontriangular blocks on the
!            diagonal of the block triangular form.
!    INFO(9) Total no. of entries in all the nontriangular blocks on
!            the diagonal of the block triangular form.
!    INFO(10) Structural rank.
!    INFO(11) Number of multiple entries in the input data.
!    INFO(12) Number of entries with out-of-range indices.
! RINFO is a real array that need not be set on entry. On exit,
!    RINFO(1) holds the number of floating-point operations needed for
!      the factorization.
! ENDCOL is an optional integer array of length N.
!      If ICNTL(8) is not 0, the ENDCOL must be set on entry so that
!      ENDCOL(i), i=1,N is zero for and only for the columns
!      designated to be at the end of the pivot sequence.

!     .. Local constants ..
      DOUBLE PRECISION zero
      PARAMETER (zero=0.0_wp)

!     .. Local variables ..
      DOUBLE PRECISION cntl5(10)
! These variables are integer because of calls to subroutines
      integer icntl5(20), nb, nc, nr, np, ndiag
      integer(long) :: info5(15)
      integer(long) :: eye, headc, i, ib, idummy, ips, &
        iptrd, iptro, isw, j, jay, jb, jfirst, &
        j1, j2, j3, k, kb, kblock, kd, kk, kl, ko, l, lastr, lastc, lblock
      integer(long), allocatable :: iptra(:),iptrp(:),ip(:),iq(:)
      integer, allocatable :: lastcol(:),lencol(:)
      LOGICAL ldup
      INTEGER lenc, lenp, lenr, lp, mblock, minblk, mp, nblock, &
        nextc, nextr, nxteye, ptrd, ptro
      integer(long) :: nza, nzb, nzd, ljfact
      DOUBLE PRECISION rinfo5(10), tol
! CNTL5 passed to MA50A/AD to correspond to dummy argument CNTL.
! EYE   position in arrays.
! HEADC displacement in IW. IW(HEADC+J) holds the position in A of the
!       head of the chain of entries for column J.
! I     row index.
! IB    displacement in IW. IW(IB+K) first holds the index in the
!       permuted matrix of the first column of the k-th block. Later
!       abs(IW(IB+K)) holds the number of columns in block K and is
!       negative for a triangular block.
! ICNTL5 passed to MA50A/AD to correspond to dummy argument ICNTL.
! IDUMMY do loop index not used within the loop.
! INFO5 passed to MA50A/AD to correspond to dummy argument INFO.
! IP    displacement in IW. IW(IP+I), I=1,M holds a row permutation.
! IPTRA displacement in IW. IW(IPTRA+J) holds the position in the
!       reordered matrix of the first entry of column J.
! IPTRD displacement in KEEP. See comment on KEEP.
! IPTRO displacement in KEEP. See comment on KEEP.
! IPTRP displacement in IW. IW(IP+I), I=1,N holds the column starts
!       of a permutation of the matrix, during the block triangular
!       calculation.
! IQ    displacement in IW. IW(IQ+J), J=1,N holds a column permutation.
! ISW   used for swopping two integers
! IW   is an integer array of length 6*M+3*N that is used as workspace.
! IW13  displacement in IW. IW(IW13+1) is the first entry of the
!       workarray IW of MC13DD.
! IW21  displacement in KEEP. KEEP(IW21+1) is the first entry of the
!       workarray IW of MC21AD.
! IW50  displacement in IW. IW(IW50+1) is the first entry of the
!       workarray IW of MA50A/AD.
! J     column index.
! JAY   column index.
! JB    block index.
! JFIRST displacement in IW. IW(JFIRST+1) is the first entry of the
!       workarray JFIRST of MA50A/AD.
! J1    first column index of a block.
! J2    last column index of a block.
! J3    last column index of a block ... used in printing
! K     running index for position in matrix.
! KB    block index.
! KBLOCK displacement in KEEP. See comment on KEEP.
! KD    running index for position in matrix.
! KK    running index for position in matrix.
! KL    last position in matrix of current column.
! KO    running index for position in matrix.
! L     length of a block.
! LBLOCK displacement in KEEP. See comment on KEEP.
! LASTR displacement in IW. IW(LASTR+I) holds the position of the last
!       entry encountered in row I.  Used for accumulating duplicates.
! LASTC displacement in IW. IW(LASTC+I) holds the index of the
!       last column encountered that had an entry in row I.
!       LASTC is also used to get workspace in KEEP for MA50AD.
! LDUP  logical flag used to indicate whether duplicates have been
!       found and recorded.
! LENC  displacement in IW. IW(LENC+J) holds the number of entries
!       in column J, excluding duplicates.
!       LENC is also used to get workspace in KEEP for MA50AD.
! LENP  displacement in IW. IW(LENP+J) holds the number of entries
!       in column J of the permuted matrix.
! LENR  displacement in IW. IW(LENR+1) is the first entry of the
!       workarray LENR of MA50A/AD.
! LP Unit for error messages.
! MBLOCK displacement in KEEP. See comment on KEEP.
! MINBLK Minimum size for a block.
! MP Unit for diagnostic messages.
! NB Number of blocks.
! NBLOCK displacement in KEEP. See comment on KEEP.
! NC Number of columns in block.
! NDIAG number entries placed on diagonal by MC21AD.
! NEXTC displacement in IW. IW(NEXTC+1) is the first entry of the
!       workarray NEXTC of MA50A/AD.
! NEXTR displacement in IW. IW(NEXTR+1) is the first entry of the
!       workarray NEXTR of MA50A/AD.
! NP Number of columns in packed storage.
! NR Number of rows in block.
! NXTEYE Next value of EYE.
! NZA   number of entries in the matrix, excluding duplicates.
! NZB   number of entries in the current diagonal block.
! NZD   total number of entries in the diagonal blocks.
! PTRD  displacement in IW. IW(PTRD+J) holds the position of the
!       first entry of the diagonal part of column J.
! PTRO  displacement in IW. IW(PTRO+J) holds the position of the
!       first entry of the off-diagonal part of column J.
! RINFO5 passed to MA50A/AD to correspond to dummy argument RINFO.
! TOL   pivot tolerance.  Entries less than this value are considered
!       as zero.

! MC13DD    finds permutation for block triangular form
! MC21AD    finds permutation for zero-free diagonal
! MA50A/AD factorizes non-triangular diagonal blocks
      INTRINSIC abs, max

      multiplier = cntl(6)

      allocate(iptra(n),iptrp(n),iq(n),ip(m),stat=stat)
      if (stat .ne. 0) then
        info(1) = -10
        return
      endif
! The +1 because pointers point to first free position but array used
! from next position?
      allocate(iw(4*max(m,n)+2*n+1),stat=stat)
      if (stat .ne. 0) then
        info(1) = -10
        deallocate(iptra,iptrp,iq,ip,stat=stat)
        return
      endif

! Assign data for columns to be operated on at end
      if (present(endcol)) iw(1:n) = endcol(1:n)
! Allocate storage for the factorization of the blocks
      ljfact = max(la/2,2*ne)
      allocate(fact(la),ifact(la),jfact(ljfact),stat=stat)
      if (stat .ne. 0) then
        info(1) = -10
        deallocate(iptra,iptrp,iq,ip,iw,stat=stat)
        return
      endif

! Initializations
      DO 53 i = 1, 20
        info(i) = 0
53    CONTINUE
      info(3) = ne*2
      info(14) = ne*2
      info(4) = ne
      info(10) = min(m,n)
      DO 56 i = 1, 5
        rinfo(i) = zero
56    CONTINUE
      DO 60 i = 1, 4
        cntl5(i) = cntl(i)
        icntl5(i) = icntl(i)
60    CONTINUE
! Switch off printing from MA50A/AD
      icntl5(3) = 0
! Set BLAS control
      icntl5(5) = icntl(5)
! Default control for restricted pivoting
      icntl5(6) = 0
! Set controls if pivot sequence provided
      icntl5(7) = 0
      IF (job==2) THEN
        icntl5(7) = 2
        icntl5(4) = 1
      END IF
! Set control for diagonal pivoting
      IF (job==3) icntl5(7) = 1

! Set pivot tolerance
      tol = max(zero,cntl(4))
! Set blcok size for BTF
      minblk = max(1,icntl(6))
! Set flag for scan for duplicates
      ldup = .FALSE.

! Simple data checks (some removed because done by hsl_ma48)
! These are executed by test deck but will not be executed by hsl_ma48
! Set printer streams
      lp = icntl(1)
      mp = icntl(2)
      IF (ne<=0) THEN
        info(1) = -2
        IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,9110) ne
        GO TO 530
      END IF
      IF (la<2*ne) THEN
        info(1) = -3
        info(3) = 2*ne
        IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,9120) la, info(3)
        GO TO 530
      END IF
      IF (job<1 .OR. job>3) THEN
        info(1) = -6
        IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,9130) job
        GO TO 530
      END IF
! Check input permutations
      IF (job==2) THEN
        DO 10 i = 1, max(m,n)
          iw(n+i) = 0
10      CONTINUE
! Check row permutation
        DO 20 i = 1, m
          j = keep(i)
          IF (j<1 .OR. j>m) GO TO 40
          IF (iw(n+j)==1) GO TO 40
          iw(n+j) = 1
20      CONTINUE
! Check column permutation
        DO 30 i = 1, n
          j = keep(m+i)
          IF (j<1 .OR. j>n) GO TO 40
          IF (iw(n+j)==2) GO TO 40
          iw(n+j) = 2
30      CONTINUE
        GO TO 50
40      info(1) = -5
        IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,9140)
        GO TO 530
      END IF

50    IF (mp>0 .AND. icntl(3)>2) THEN
        WRITE (mp,'(/A/A,I7,A,I6,A,I7,A,I2,A,I7/A,1P,4D12.4/A,4I8/A,3I8)') &
          ' Entering MA48A/AD with', ' M =', m, '     N =', n, '     NE =', &
          ne, '     JOB =', job, '     LA =', la, ' CNTL (1:4) =', &
          (cntl(i),i=1,4), ' ICNTL(1:4) = ', (icntl(i),i=1,4), &
          ' ICNTL(6:8) = ', (icntl(i),i=6,8)
        IF (icntl(3)>3) THEN
          WRITE (mp,9000) (a(k),irn(k),jcn(k),k=1,ne)
9000      FORMAT (' Entries:'/3(1P,D12.4,2I6))
        ELSE
          WRITE (mp,9000) (a(k),irn(k),jcn(k),k=1,min(onel*9,ne))
        END IF
        IF (job==2) THEN
          WRITE (mp,'(A)') ' Permutations input (JOB=2)'
          IF (icntl(3)>3) THEN
            WRITE (mp,9010) (keep(i),i=1,m)
9010        FORMAT (' Positions of original rows in the permuted ma', &
              'trix'/(10I6))
            WRITE (mp,9020) (keep(m+i),i=1,n)
9020        FORMAT (' Positions of columns of permuted matrix ','in',' or', &
              'iginal matrix '/(10I6))
          ELSE
            WRITE (mp,9010) (keep(i),i=1,min(10,m))
            WRITE (mp,9020) (keep(m+i),i=1,min(10,n))
          END IF
        END IF
        IF (icntl(8)/=0) THEN
          WRITE (mp,'(A,I6)') ' Value of IW entries on call with ICNTL(8) =', &
            icntl(8)
          IF (icntl(3)>3) THEN
            WRITE (mp,9030) (iw(i),i=1,n)
9030        FORMAT (10I6)
          ELSE
            WRITE (mp,9030) (iw(i),i=1,min(10,n))
          END IF
        END IF
      END IF


! Partition KEEP
      iptrd = m + 3*n
      iptro = iptrd + n + 1
      nblock = iptro + n - 1
      mblock = nblock + 1
      kblock = mblock + 1

! Partition IW
! We assume at this point that there might be data in IW(1:N)
      headc = n + 1
      lastc = headc + n

! Initialize header array for column links.
      DO 70 j = 1, n
        iw(headc+j) = 0
70    CONTINUE
! Overwrite JCN by column links, checking that rows and column numbers
!     are within range.
      DO 80 k = 1, ne
        i = irn(k)
        j = jcn(k)
        IF (i<1 .OR. i>m .OR. j<1 .OR. j>n) THEN
          info(12) = info(12) + 1
          IF (mp>0 .AND. info(12)<=10 .AND. icntl(3)>=2) WRITE (mp, &
            '(A,I7,A,2I6)') ' Message from MA48A/AD .. indices for entry ', k, &
            ' are', i, j
          jcn(k) = 0
        ELSE
          jcn(k) = iw(headc+j)
          iw(headc+j) = k
        END IF
80    CONTINUE


      IF (minblk>=n .OR. m/=n .OR. job>1) GO TO 190

! Obtain permutations to put matrix in block triangular form.

! Note that in this part of the code we know that M is equal to N.

! Partition IW for arrays used by MC21AD and MC13DD
      lenc = lastc + n
      ib = lenc
      ips = lenc + n
      lenp = ips + n

! Initialize array LASTC
      DO 90 i = 1, n
        iw(lastc+i) = 0
90    CONTINUE

! Run through the columns removing duplicates and generating copy of
!     structure by columns.
      ldup = .TRUE.
      k = 1
      DO 120 j = 1, n
        eye = iw(headc+j)
        iptra(j) = k
        DO 100 idummy = 1, ne
          IF (eye==0) GO TO 110
          i = irn(eye)
! Check for duplicates
          IF (iw(lastc+i)/=j) THEN
            iw(lastc+i) = j
            irn(ne+k) = i
            k = k + 1
          ELSE
! Record duplicate
            info(11) = info(11) + 1
            IF (mp>0 .AND. info(11)<=10 .AND. icntl(3)>=2) WRITE (mp, &
              '(A,I7,A,2I6)') &
              ' Message from MA48A/AD .. duplicate in position ', k, &
              ' with indices', i, j
          END IF
          eye = jcn(eye)
100     CONTINUE
110     iw(lenc+j) = k - iptra(j)
120   CONTINUE

! NZA is number of entries with duplicates and out-of-range indices
!     removed.
      nza = k - 1

!   Compute permutation for zero-free diagonal
      CALL mc21ad(n,irn(ne+1),nza,iptra,iw(lenc+1),iw(ips+1),ndiag,info)
      if (info(1) .eq. -10) go to 530
! IW(IP+1) ... col permutation ... new.to.old
      info(10) = ndiag
      IF (ndiag<n) THEN
        IF (icntl(7)/=0) THEN
          info(1) = -4
          IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,'(A,A/A,I7,A,I7)') &
            ' Error return from MA48A/AD because matrix structurally ', &
            ' singular', ' order is ', n, ' and structural rank', ndiag
          GO TO 530
        END IF
        GO TO 190
      END IF
! Set pointers and column counts for permuted matrix
!DIR$ IVDEP
      DO 130 j = 1, n
        jay = iw(ips+j)
        iptrp(j) = iptra(jay)
        iw(lenp+j) = iw(lenc+jay)
130   CONTINUE

! Use Tarjan's algorithm to obtain permutation to block triangular form.
! KEEP(M+1) .. col permutation ...new.to.old
      CALL mc13dd(n,irn(ne+1),nza,iptrp,iw(lenp+1),keep(m+1),iw(ib+1),nb,info)
      if (info(1) .eq. -10) go to 530
! IW(IB+1) ... col number in permuted matrix of start of blocks

! Merge adjacent 1x1 blocks into triangular blocks
! Change block pointers to block sizes
      DO 140 jb = 2, nb
        iw(ib+jb-1) = iw(ib+jb) - iw(ib+jb-1)
140   CONTINUE
      iw(ib+nb) = n + 1 - iw(ib+nb)
      IF (iw(ib+1)==1) iw(ib+1) = -1
      kb = 1
      DO 150 jb = 2, nb
        l = iw(ib+jb)
        IF (l==1 .AND. iw(ib+kb)<=0) THEN
          iw(ib+kb) = iw(ib+kb) - 1
        ELSE
          kb = kb + 1
          IF (l==1) THEN
            iw(ib+kb) = -1
          ELSE
            iw(ib+kb) = l
          END IF
        END IF
150   CONTINUE
      nb = kb
! Merge small blocks
      kb = 1
      DO 160 jb = 2, nb
        IF (abs(iw(ib+kb))<minblk) THEN
          iw(ib+kb) = abs(iw(ib+kb)) + abs(iw(ib+jb))
        ELSE
          kb = kb + 1
          iw(ib+kb) = iw(ib+jb)
        END IF
160   CONTINUE
      nb = kb
! Set size of blocks and triangular block flags
      DO 170 jb = 1, nb
        keep(nblock+3*jb) = abs(iw(ib+jb))
        keep(mblock+3*jb) = iw(ib+jb)
170   CONTINUE
! Record permutations.
! Set KEEP to position of column that is in position j of permuted
!     matrix .... new.to.old
      DO 180 j = 1, n
        keep(keep(m+j)) = j
        keep(m+j) = iw(ips+keep(m+j))
180   CONTINUE
      GO TO 220

! Set arrays for the case where block triangular form is not computed
190   nb = 1
      IF (job==1 .OR. job==3) THEN
        DO 200 i = 1, m
          keep(i) = i
200     CONTINUE
        DO 210 i = 1, n
          keep(m+i) = i
210     CONTINUE
      END IF
! Set number of columns in block
      keep(nblock+3) = n
      keep(mblock+3) = 0

! Control for forcing columns to end.
220   IF (icntl(8)/=0) THEN
        lblock = kblock + 3*nb
        IF (job==2) THEN
! Find first column that will be changed in subsequent factorizations.
          DO 230 i = 1, n
            IF (iw(i)==0) GO TO 240
230       CONTINUE
! Set LBLOCK to number of columns from and including first one flagged
240       keep(lblock+1) = n - i + 1
        ELSE
! Within each block move columns to the end
          j = 1
          DO 270 jb = 1, nb
            keep(lblock+jb) = 0
            j2 = j + keep(nblock+3*jb) - 1
            j1 = j2
! Jump if triangular block (leave columns in situ)
            IF (keep(mblock+3*jb)<0) GO TO 260
! Run through columns in block sorting those with IW flag to end
250         IF (j==j2) GO TO 260
            IF (iw(keep(m+j))==0) THEN
! Column is put at end
              keep(lblock+jb) = keep(lblock+jb) + 1
              isw = keep(m+j2)
              keep(m+j2) = keep(m+j)
              keep(m+j) = isw
              j2 = j2 - 1
            ELSE
              j = j + 1
            END IF
            GO TO 250
260         j = j1 + 1
270       CONTINUE
        END IF
      END IF

! Run through the columns to create block-ordered form, removing
!     duplicates, changing row indices, and holding map array in JCN.
! Grab space in IW for LASTR
      lastr = lastc + m
! Initialize LASTC
      DO 280 i = 1, m
        iw(lastc+i) = 0
280   CONTINUE
      keep(kblock+3) = nb
      k = 1
      kk = ne
      j2 = 0
      DO 310 jb = 1, nb
        j1 = j2 + 1
        j2 = j1 + keep(nblock+3*jb) - 1
! Run through columns in block JB
! LASTR is used for accumulation of duplicates
! LASTC is used to identify duplicates
        DO 300 jay = j1, j2
          j = keep(m+jay)
          eye = iw(headc+j)
          keep(iptrd+jay) = k
          keep(iptro+jay) = kk
          IF (keep(mblock+3*jb)<0) THEN
! Block is triangular.
! Reserve the leading position for the diagonal entry
            iw(lastc+jay) = jay
            a(ne+k) = zero
            irn(ne+k) = jay
            iw(lastr+jay) = k
            k = k + 1
          END IF
          DO 290 idummy = 1, ne
            IF (eye==0) GO TO 300
            nxteye = jcn(eye)
            i = keep(irn(eye))
            IF (iw(lastc+i)/=jay) THEN
! Entry encountered for the first time.
              iw(lastc+i) = jay
              IF ((i>=j1 .AND. i<=j2) .OR. (m/=n)) THEN
! Entry in diagonal block.
                a(ne+k) = a(eye)
                irn(ne+k) = i
                iw(lastr+i) = k
                jcn(eye) = k
                k = k + 1
              ELSE
! Entry in off-diagonal block.
                a(ne+kk) = a(eye)
                irn(ne+kk) = i
                iw(lastr+i) = kk
                jcn(eye) = kk
                kk = kk - 1
              END IF
            ELSE
! Entry has already been encountered.
              IF ( .NOT. ldup) THEN
! Duplicates should be counted here if they were not earlier
                info(11) = info(11) + 1
                IF (mp>0 .AND. info(11)<=10 .AND. icntl(3)>=2) WRITE (mp, &
                  '(A,I7,A,2I6)') &
                  ' Message from MA48A/AD .. duplicate in position ', eye, &
                  ' with indices', irn(eye), j
              END IF
              kl = iw(lastr+i)
              jcn(eye) = kl
              a(ne+kl) = a(ne+kl) + a(eye)
            END IF
            eye = nxteye
290       CONTINUE
300     CONTINUE
310   CONTINUE
      keep(iptrd+n+1) = k
      keep(iptro+n+1) = kk
      nzd = k - 1
      DO 320 i = 1, k - 1
        irn(i) = irn(ne+i)
320   CONTINUE
      DO 325 i = kk + 1, ne
        irn(i) = irn(ne+i)
325   CONTINUE
      DO 326 i = k, kk
        irn(i) = 0
326   CONTINUE

! Repartition IW (none of previous data is needed any more)
      ptrd = 0
      jfirst =  n
      lenr = m + n
      lastr = 2*m + n
      nextr = 3*m + n
      nextc = 4*m + n
      ptro = nextc
! Use scratch space in KEEP for work arrays in MA50AD.
! Use an allocatable work array for this now
!      lastc = m + n
!      lenc = lastc + n
       allocate(lastcol(n),lencol(n),stat=stat)
      if (stat .ne. 0) then
        info(1) = -10
        go to 530
      endif

! Call MA50A/AD block by block from the back.
      j1 = n + 1
! KB counts number of non-triangular blocks for calculation INFO(4)
      kb = 0
      DO 390 jb = nb, 1, -1
! NC is number of columns in block
        nc = keep(nblock+3*jb)
        j2 = j1 - 1
        j1 = j2 + 1 - nc
        IF (keep(mblock+3*jb)<0) THEN
! Action if triangular block
          DO 330 j = j1, j2
! Increment estimate of rank
            IF (abs(a(ne+keep(iptrd+j)))>tol) info(5) = info(5) + 1
            ip(j) = j
            iq(j) = j
330       CONTINUE
        ELSE
! Action if non-triangular block
! Copy block to allocatable work array, shifting indices.
          nzb = keep(iptrd+j2+1) - keep(iptrd+j1)
          kk = 0
          do k = keep(iptrd+j1), keep(iptrd+j2+1) - 1
!           irn(ne+k) = irn(k) - j1 + 1
            kk = kk + 1
            ifact(kk) = irn(k) - j1 + 1
!!!
!  In Version 3.0.0 fact(kk) was being set to a(k)
            fact(kk)  = a(ne+k)
          enddo
! Copy the column pointers, shifting them.
! K is position immediately before current block in IRN
          k = keep(iptrd+j1) - 1
          DO 350 j = j1, j2
            iq(j) = keep(iptrd+j) - k
350       CONTINUE
! NR is number of rows in block
          nr = nc
          IF (nb==1) nr = m
! Set permutations
! Set to identity because matrix has already been permuted to put
!   pivots on the diagonal.
          IF (job==2) THEN
            DO 360 j = j1, j1 + nr - 1
              ip(j) = j - j1 + 1
360         CONTINUE
            DO 370 j = j1, j2
              iw(ptrd+j) = j - j1 + 1
370         CONTINUE
          END IF
! Accumulate data on block triangular form
          info(7) = max(info(7),onel*nr)
          info(8) = info(8) + nc
          info(9) = info(9) + nzb
! Control for forcing columns to end.
          IF (icntl(8)/=0) icntl5(6) = keep(lblock+jb)
!         CALL ma50ad(nr,nc,nzb,la-ne-k,a(ne+k+1),irn(ne+k+1),jcn(ne+1), &
          CALL ma50ad(nr,nc,nzb,la,fact,ifact,ljfact,jfact, &
            iq(j1),cntl5,icntl5,ip(j1),np,iw(jfirst+1),iw(lenr+1), &
            iw(lastr+1),iw(nextr+1),iw(ptrd+j1),lencol(kb+1), &
            lastcol(kb+1),iw(nextc+1),multiplier,info5,rinfo5)
! Action if failure in allocate in ma50ad
          if (info5(1) .eq. -2*nr) then
            info(1) = -10
            go to 530
          endif
! KEEP(MBLOCK+3*JB) is number rows in packed storage
          keep(mblock+3*jb) = np
! Shift the permutations back
          DO 380 j = j1, j1 + nr - 1
            ip(j) = ip(j) + j1 - 1
380       CONTINUE
          DO 385 j = j1, j2
            iq(j) = iq(j) + j1 - 1
385       CONTINUE
! Adjust warning and error flag
          IF (info5(1)==1) THEN
            IF (info(1)==0 .OR. info(1)==4) info(1) = info(1) + 2
          END IF
          IF (info5(1)==2) THEN
            IF (info(1)==0 .OR. info(1)==2) info(1) = info(1) + 4
          END IF
          IF (info5(1)==3 .AND. info(1)>=0) info(1) = 6
!XXX not possible
!         IF (info5(1)==-3) info(1) = -3
! Accumulate number of garbage collections
          info(2) = info(2) + info5(2)
! Accumulate real space for subsequent analyse
          info(3) = max(info(3),ne+k+info5(3))
! Accumulate integer space for subsequent analyse
          info(14) = max(info(14),ne+k+info5(8))
! Accumulate estimate of rank
          info(5) = info(5) + info5(5)
! Accumulate number of dropped entries
          info(6) = info(6) + info5(6)
! Accumulate floating-point count
          rinfo(1) = rinfo(1) + rinfo5(1)
! Calculate data necessary for INFO(4)
! LENCOL(KB) holds the actual storage forecast for MA50BD entry.
          kb = kb + 1
          lencol(kb) = info5(4) - 2*nr
! LASTCOL(KB) holds the storage plus elbow room forecast
!         for MA50BD entry.
          lastcol(kb) = info5(4)
        END IF
390   CONTINUE

! Now calculate storage required for subsequent factorization.
      info(4) = ne*2
      k = ne
      do jb = kb, 1, -1
        info(4) = max(info(4),k+lastcol(jb))
        k = k + lencol(jb)
      enddo

      deallocate(lastcol,lencol,stat=stat)
      if (stat .ne. 0) then
        info(1) = -10
        go to 530
      endif

!XXX Cannot be executed .... remove?
!     IF (info(1)==-3) THEN
!       IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,9120) la, info(3)
!       GO TO 530
!     END IF


! Reorder IRN(1:NE) to the revised column permutation, storing in
!     IRN(NE+k) the new position for the entry that was in position k.
! Copy IRN data back.
      DO 410 k = 1, ne
        irn(ne+k) = irn(k)
410   CONTINUE
      DO 420 j = 1, n
        iw(ptrd+j) = keep(iptrd+j)
        iw(ptro+j) = keep(iptro+j+1) + 1
420   CONTINUE
      kd = 1
      ko = nzd + 1
      DO 450 j = 1, n
        keep(iptrd+j) = kd
        jay = iq(j)
        kl = nzd
        IF (jay/=n) kl = iw(ptrd+jay+1) - 1
        DO 430 kk = iw(ptrd+jay), kl
          irn(kd) = ip(irn(ne+kk))
          irn(ne+kk) = kd
          kd = kd + 1
430     CONTINUE
        keep(iptro+j) = ko
        kl = ne
        IF (jay/=1) kl = iw(ptro+jay-1) - 1
        DO 440 kk = iw(ptro+jay), kl
          irn(ko) = ip(irn(ne+kk))
          irn(ne+kk) = ko
          ko = ko + 1
440     CONTINUE
450   CONTINUE
      keep(iptro+n+1) = ko

! Compute the product permutations
      DO 460 i = 1, m
        keep(i) = ip(keep(i))
460   CONTINUE
      DO 465 i = 1, n
        iq(i) = keep(m+iq(i))
465   CONTINUE
      DO 470 i = 1, n
        keep(m+i) = iq(i)
470   CONTINUE

! Compute (product) map
! Set IRN(NE) to deal with out-of-range indices
      iw(1) = irn(ne)
      irn(ne) = ne
      DO 480 k = 1, ne
        jcn(k) = irn(ne+jcn(k))
480   CONTINUE
! Reset IRN(NE)
      irn(ne) = iw(1)

! Warning messages
      IF (info(11)>0 .OR. info(12)>0) THEN
        info(1) = info(1) + 1
        IF (mp>0 .AND. icntl(3)>=2) THEN
          IF (info(11)>0) WRITE (mp,9150) info(11)
          IF (info(12)>0) WRITE (mp,9160) info(12)
        END IF
        IF (info(11)>0) jcn(1) = -jcn(1)
      END IF

! Set rank estimate to min of structural and INFO(5) estimate.
      IF (info(10)<info(5)) THEN
        info(5) = info(10)
        IF (info(1)/=2 .AND. info(1)/=3 .AND. info(1)<6) info(1) = info(1) + 2
      END IF

      IF (mp>0 .AND. icntl(3)>2) THEN
        WRITE (mp,'(/A/A/12I6/A,F12.1)') ' Leaving MA48A/AD with', &
          ' INFO(1:12)  =', (info(i),i=1,12), ' RINFO(1) =', rinfo(1)
        WRITE (mp,'(A)') ' Permuted matrix by blocks (IRN)'
        kb = nb
        IF (icntl(3)==3) THEN
          WRITE (mp,'(A)') ' Only first column of up to 10 blocks printed'
          kb = min(10,nb)
        END IF
        WRITE (mp,'(A)') ' Diagonal blocks'
        j1 = 1
        DO 500 jb = 1, kb
          j2 = j1 + keep(nblock+3*jb) - 1
          IF (j1<=j2) WRITE (mp,'(A,I6)') ' Block', jb
          j3 = j2
          IF (icntl(3)==3) j3 = j1
          DO 490 j = j1, j3
            WRITE (mp,'(A,I5,(T13,10I6))') ' Column', j, &
              (irn(i),i=keep(iptrd+j),keep(iptrd+j+1)-1)
490       CONTINUE
          j1 = j2 + 1
500     CONTINUE
        IF (keep(iptro+n+1)>keep(iptrd+n+1)) THEN
          WRITE (mp,'(A)') ' Off-diagonal entries'
          j1 = 1
          DO 520 jb = 1, kb
            j2 = j1 + keep(nblock+3*jb) - 1
            j3 = j2
            IF (icntl(3)==3) j3 = j1
            DO 510 j = j1, j3
              IF (keep(iptro+j+1)>keep(iptro+j)) WRITE (mp, &
                '(A,I5,(T13,10I6))') ' Column', j, (irn(i),i=keep(iptro+j), &
                keep(iptro+j+1)-1)
510         CONTINUE
            j1 = j2 + 1
520       CONTINUE
        END IF
        IF (icntl(3)>3) THEN
          WRITE (mp,9040) (jcn(k),k=1,ne)
9040      FORMAT (' JCN (MAP) ='/(6X,10I6))
          WRITE (mp,'(A)') ' Contents of KEEP array'
          WRITE (mp,9010) (keep(i),i=1,m)
          WRITE (mp,9020) (keep(m+i),i=1,n)
          WRITE (mp,'(A)') ' Pointer information from KEEP'
          WRITE (mp,9050) (keep(iptrd+j),j=1,n+1)
9050      FORMAT (' IPTRD ='/(8X,10I6))
          WRITE (mp,9060) (keep(iptro+j),j=1,n+1)
9060      FORMAT (' IPTRO ='/(8X,10I6))
          WRITE (mp,'(A)') ' Block structure information from KEEP'
          WRITE (mp,9070) (keep(nblock+3*jb),jb=1,nb)
          WRITE (mp,9080) (keep(mblock+3*jb),jb=1,nb)
9070      FORMAT (' NBLOCK (order blocks) ='/(8X,10I6))
9080      FORMAT (' MBLOCK (triangular flag and number packed rows) ='/ &
            (8X,10I6))
9090      FORMAT (' LBLOCK (number of changed columns) ='/(8X,10I6))
          IF (icntl(8)/=0) WRITE (mp,9090) (keep(lblock+jb),jb=1,nb)
        ELSE
          WRITE (mp,9040) (jcn(k),k=1,min(10*onel,ne))
          WRITE (mp,'(A)') ' Contents of KEEP array'
          WRITE (mp,9010) (keep(i),i=1,min(10,m))
          WRITE (mp,9020) (keep(m+i),i=1,min(10,n))
          WRITE (mp,'(A)') ' Pointer information from KEEP'
          WRITE (mp,9050) (keep(iptrd+j),j=1,min(10,n+1))
          WRITE (mp,9060) (keep(iptro+j),j=1,min(10,n+1))
          WRITE (mp,'(A)') ' Block structure information from KEEP'
          WRITE (mp,9070) (keep(nblock+3*jb),jb=1,min(10,nb))
          WRITE (mp,9080) (keep(mblock+3*jb),jb=1,min(10,nb))
          IF (icntl(8)/=0) WRITE (mp,9090) (keep(lblock+jb),jb=1,min(10,nb))
        END IF
      END IF

! Deallocate arrays
530   deallocate(iptra,iptrp,iq,ip,stat=stat)
      deallocate(iw,stat=stat)
      deallocate(fact,ifact,jfact,stat=stat)
      if (stat .ne. 0) then
        info(1) = -10
        return
      endif
      RETURN

9110  FORMAT (' Error return from MA48A/AD because NE =',I10)
9120  FORMAT (' Error return from MA48A/AD because LA is',I10/' and ', &
        'must be at least',I10)
9130  FORMAT (' Error return from MA48A/AD because ','JOB = ',I10)
9140  FORMAT (' Error return from MA48A/AD because ','faulty permutati', &
        'ons input when JOB = 2')
9150  FORMAT (' Message from MA48A/AD ..',I8,' duplicates found')
9160  FORMAT (' Message from MA48A/AD ..',I8,' out-of-range indices fo','und')
    END subroutine ma48ad


    SUBROUTINE ma48bd(m,n,ne,job,la,a,irn,jcn,keep,cntl,icntl,info,rinfo)
! Factorize a sparse matrix, using data provided by MA48A/AD.

!     .. Arguments ..
      INTEGER m, n, job
      integer(long) :: ne, la
      DOUBLE PRECISION a(la)
      INTEGER(long) :: irn(la), jcn(ne)
      integer(long) :: keep(*)
      DOUBLE PRECISION cntl(10)
      INTEGER icntl(20)
      real(wp), allocatable :: w(:)
      integer(long), allocatable :: iw(:)
      INTEGER(long) :: info(20)
      DOUBLE PRECISION rinfo(10)

! M must be set by the user to the number of rows.
!      It is not altered by the subroutine. Restriction: M > 0.
! N must be set by the user to the number of columns.
!      It is not altered by the subroutine. Restriction: N > 0.
! NE must be set by the user to the number of entries in the input
!      matrix. It is not altered by the subroutine. Restriction: NE > 0.
! JOB  must be set by the user to 1 for a normal entry or to 2 for a
!      fast entry. JOB is not altered by the subroutine.
!      Restriction 1 <= JOB <= 2.
! LA must be set by the user to the size of A and IRN.
!      It is not altered by the subroutine. Restriction: LA > 2*NE.
! A    must be set by the user so that the first NE elements hold the
!      matrix entries in exactly the same order as were input to
!      MA48A/AD. On return, these elements hold the permuted matrix A
!      and the rest contains the factorizations of the block diagonal
!      part (excluding blocks that are triangular).
! IRN  The first KEEP(IPTRO+N+1)-1 (<=NE) entries must be as on return
!      from MA48A/AD and hold
!      the row numbers of the permuted matrix, with duplicates excluded.
!      The permuted matrix is block upper triangular with recommended
!      pivots on its diagonal. The rest of the array need not be set on
!      entry with JOB=1. It is used for the row indices of the
!      factorizations. On entry with JOB=2, it must be unchanged
!      since the entry with JOB=1 and is not altered.
! JCN  must be as on return from MA48A/AD. abs(JCN(k)) holds the
!      position in IRN that corresponds to the entry that was input in
!      A(k),k=1,NE.  It is not altered by the subroutine.
!      If JCN(1) is negative then there are duplicates.
! KEEP must be as on return from MA48A/AD or MA48B/BD.
!      KEEP(i) holds the position of row i in the permuted
!        matrix, I=1,M;
!      KEEP(N+j) holds the index of the  column that is in position j
!        of the permuted matrix, j=1,N;
!      KEEP(IPTRL+k), IPTRL=M+N, k=1,N, need not be set on an
!        entry with JOB=1. It holds pointers, relative to the start
!        of the block, to the columns of the
!        lower-triangular part of the factorization. On entry with
!        JOB=2, it must be unchanged since the entry with JOB=1 and is
!        not altered.
!      KEEP(IPTRU+k), IPTRU=IPTRL+N, k=1,N, need not be set on an
!        entry with JOB=1. It holds pointers, relative to the start
!        of the block, to the columns of the
!        upper-triangular part of the factorization. On an entry with
!        JOB=2, it must be unchanged since the entry with JOB=1 and
!        is not altered.
!      KEEP(IPTRD+j), IPTRD=IPTRU+N, holds the position in IRN of the
!        start of the block diagonal part of column j, J=1,N;
!      KEEP(IPTRD+N+1) holds the position that immediately follows
!        the end of the block diagonal part of column N;
!      KEEP(IPTRO+j),IPTRO=IPTRD+N+1, holds the position in IRN of
!        the start of the block off-diagonal part of column j, j=1,N;
!      KEEP(IPTRO+N+1) holds the position that immediately
!        follows the end of the block off-diagonal part of column N;
!      KEEP(KBLOCK+3), KBLOCK=IPTRO+N+1, holds the number of
!        blocks NB in the block triangular form;
!      KEEP(NBLOCK+3*k), NBLOCK=IPTRO+N-1, holds the number of columns
!        in block k, k=1,NB;
!      KEEP(MBLOCK+3*k), MBLOCK=IPTRO+N, is negative if block k
!        is triangular or holds the number of rows held in packed
!        storage when processing block k, k=1,NB;
!      KEEP(KBLOCK+3*k), k=2,NB need not be set on an entry with
!        JOB=1; on return it holds zero for a triangular block and for
!        a non-triangular block holds the position in A(NEWNE+1) and
!        IRN(NEWNE+1) of the start of the factorization of the block.
!        KEEP(LBLOCK+k), LBLOCK=KBLOCK+3*NB is accessed only if ICNTL(8)
!          is not 0 and holds the number of columns in
!          block k which are to be pivoted on last, k=1,NB.
!      On an entry with JOB=2, KEEP must be unchanged since the entry
!      with JOB=1 and is not altered.
! CNTL  must be set by the user as follows and is not altered.
!     CNTL(2) determines the balance used by MA50A/AD and MA50B/BD
!       between pivoting for sparsity and for stability, values near
!       zero emphasizing sparsity and values near one emphasizing
!       stability.
!     CNTL(3) If this is set to a positive value, any entry whose
!       modulus is less than CNTL(3) will be dropped from the factors
!       calculated by MA50A/AD. The factorization will then require
!       less storage but will be inaccurate.
!     CNTL(4)  If this is set to a positive value, any entry whose
!       modulus is less than CNTL(4) will be regarded as zero from
!       the point of view of rank.
!     CNTL(5:10) are not referenced by the subroutine.
! ICNTL must be set by the user as follows and is not altered.
!     ICNTL(1)  must be set to the stream number for error messages.
!       A value less than 1 suppresses output.
!     ICNTL(2) must be set to the stream number for diagnostic output.
!       A value less than 1 suppresses output.
!     ICNTL(3) must be set to control the amount of output:
!       0 None.
!       1 Error messages only.
!       2 Error and warning messages.
!       3 As 2, plus scalar parameters and a few entries of array
!         parameters on entry and exit.
!       4 As 2, plus all parameters on entry and exit.
!     ICNTL(4)  If set to a positive value, the pivot search by MA50A/AD
!       is limited to ICNTL(4) columns. This may result in different
!       fill-in and execution time.
!     ICNTL(5) The block size to be used for full-matrix processing.
!     ICNTL(6) is not referenced by the subroutine.
!     ICNTL(8) If set to a value other than zero and JOB = 2 or 3,
!       it is assumed that only the last KEEP(LBLOCK+JB), JB = 1,NB
!       columns have changed since the last factorization.
!     ICNTL(9) is not referenced by the subroutine.
!     ICNTL(10) has default value 0. If set to 1, there is an immediate
!       return from MA48B/BD if LA is too small, without continuing the
!       decomposition to compute the size necessary.
!     ICNTL(11) has default value 0. If set to 1 on a JOB=2 call to
!       MA48B/BD and the entries in one of the blocks on the diagonal
!       are unsuitable for the pivot sequence chosen on the previous
!       call, the block is refactorized as on a JOB=1 call.
!     ICNTL(12:20) are not referenced by the subroutine.

! W    is a workarray.
! IW   is a workarray.
! INFO(1) is an integer variable that need not be set on entry. On exit,
!      a nonzero value of INFO(1) will indicate an error return.
!      Possible nonzero values are:
!      -1  M < 1 or N < 1
!      -2  NE < 0
!      -3  Insufficient space
!      -6  JOB out of range or JOB = 2 or 3 after factorization in which
!          entries were dropped.
!      -7  On a call with JOB=2, the matrix entries are numerically
!          unsuitable
!      +2  Matrix is structurally rank deficient. Estimated rank in
!          INFO(5).
! INFO(4) is set to the minimum size required for a further
!      factorization on a matrix with the same pivot sequence.
! INFO(5) is set to computed rank of the matrix.
! INFO(6) is set to number of entries dropped from structure.
! RINFO need not be set on entry. On exit, RINFO(1) holds the number of
!    floating-point operations performed.

!     .. Local constants ..
      DOUBLE PRECISION zero
      PARAMETER (zero=0.0_wp)

!     .. Local variables ..
      DOUBLE PRECISION cntl5(10)
      integer icntl5(20), job5, nb, nr, nc, np
      integer(long) :: i, info5(15), iptrd, iptrl, iptro, iptru, iqb(1), &
        itry, j, jb, j1, j2, j3, k, kb, kblock, kk, lblock, lp, mblock, &
        mp, nblock, newne, nrf, nzb
      DOUBLE PRECISION rinfo5(10), tol
      LOGICAL trisng
! I Temporary DO index.
! ICNTL5 passed to MA50B/BD to correspond to dummy argument ICNTL.
! INFO5 passed to MA50B/BD to correspond to dummy argument INFO.
! IPTRD displacement in KEEP. See comment on KEEP.
! IPTRL displacement in KEEP. See comment on KEEP.
! IPTRO displacement in KEEP. See comment on KEEP.
! IPTRU displacement in KEEP. See comment on KEEP.
! IQB   passed to MA50B/BD to indicate that no column permutation is
!       required.
! ITRY  Loop index for retrying the factorization of a block.
! J     column index.
! JB    block index.
! JOB5  passed to MA50B/BD to correspond to dummy argument JOB
! J1    first column index of a block.
! J2    last column index of a block.
! J3    used in prints to hold last column index of block.
! K     running index for position in matrix.
! KB  number of blocks ... used in printing
! KBLOCK displacement in KEEP. See comment on KEEP.
! KK    running index for position in matrix.
! LBLOCK displacement in KEEP. See comment on KEEP.
! LP Unit for error messages.
! MBLOCK displacement in KEEP. See comment on KEEP.
! MP Unit for diagnostic messages.
! NBLOCK displacement in KEEP. See comment on KEEP.
! NB  number of blocks
! NC    number of columns in block
! NEWNE number of entries in original matrix with duplicates and
!       out-of-range indices removed
! NP Number of columns in packed storage.
! NR    number of rows in block
! NRF   is number of rows in full form for current block.
! NZB   number of entries in current diagonal block
! RINFO5 passed to MA50B/BD to correspond to dummy argument RINFO.
! TOL   pivot tolerance
! TRISNG is flag to indicate singularity in triangular block.

!     .. Externals ..
      INTRINSIC max

      lp = icntl(1)
      mp = icntl(2)
! Simple data checks
!!!
! Although these will not happen when called from HSL_MA48
! They are tested in internal calls from current test deck
      IF (m<=0 .OR. n<=0) THEN
        info(1) = -1
        IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,9160) m, n
        GO TO 240
      END IF
      IF (ne<=0) THEN
        info(1) = -2
        IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,9170) ne
        GO TO 240
      END IF
      IF (la<2*ne) THEN
        info(1) = -3
        info(4) = 2*ne
        IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,9180) la, info(4)
        GO TO 240
      END IF
      IF (job<1 .OR. job>3) THEN
        info(1) = -6
        IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,9190) job
        GO TO 240
      END IF
      info(1) = 0
      DO 10 i = 1, 4
        cntl5(i) = cntl(i)
        icntl5(i) = icntl(i)
10    CONTINUE
! Switch off printing from MA50B/BD
      icntl5(3) = 0
! Set BLAS control
      icntl5(5) = icntl(5)
! Initialize restricted pivoting control
      icntl5(6) = 0
      icntl5(7) = 0
      icntl5(8) = icntl(10)
      allocate(w(m),iw(2*m+2*n),stat=stat)
      if (stat .ne. 0) then
        info(1) = -10
        return
      endif
! Partition KEEP and A
      iptrl = m + n
      iptru = iptrl + n
      iptrd = iptru + n
      iptro = iptrd + n + 1
      nblock = iptro + n - 1
      mblock = nblock + 1
      kblock = mblock + 1
      nb = keep(kblock+3)
      lblock = kblock + 3*nb
      newne = keep(iptro+n+1) - 1

      IF (mp>0 .AND. icntl(3)>2) THEN
        WRITE (mp,'(/A/3(A,I8),A,I2,A,I8/A,1P,3D12.4/A,3I8/A,I8/A,I8)') &
          ' Entering MA48B/BD with', ' M =', m, '     N =', n, '     NE =', &
          ne, '     JOB =', job, '     LA =', la, ' CNTL (2:4) =', &
          (cntl(i),i=2,4), ' ICNTL(1:3) =', (icntl(i),i=1,3), ' ICNTL(5)   =', &
          icntl(5), ' ICNTL(8)   =', icntl(8)
        IF (icntl(3)>3) THEN
          WRITE (mp,9000) (a(k),k=1,ne)
        ELSE
          WRITE (mp,9000) (a(k),k=1,min(10*onel,ne))
        END IF
9000    FORMAT (' A ='/(4X,1P,5D12.4))
        WRITE (mp,'(A)') ' Indices for permuted matrix by blocks'
        kb = nb
        IF (icntl(3)==3) THEN
          WRITE (mp,'(A)') ' Only first column of up to 10 blocks printed'
          kb = min(10,nb)
        END IF
        WRITE (mp,'(A)') ' Diagonal blocks'
        j1 = 1
        DO 30 jb = 1, kb
          WRITE (mp,'(A,I6)') ' Block', jb
          j2 = j1 + keep(nblock+3*jb) - 1
          j3 = j2
          IF (icntl(3)==3) j3 = j1
          DO 20 j = j1, j3
            WRITE (mp,'(A,I5,(T13,10I6))') ' Column', j, &
              (irn(i),i=keep(iptrd+j),keep(iptrd+j+1)-1)
20        CONTINUE
          j1 = j2 + 1
30      CONTINUE
        IF (keep(iptro+n+1)>keep(iptrd+n+1)) THEN
          WRITE (mp,'(A)') ' Off-diagonal entries'
          j1 = 1
          DO 50 jb = 1, kb
            j2 = j1 + keep(nblock+3*jb) - 1
            j3 = j2
            IF (icntl(3)==3) j3 = j1
            DO 40 j = j1, j3
              IF (keep(iptro+j+1)>keep(iptro+j)) WRITE (mp, &
                '(A,I5,(T13,10I6))') ' Column', j, (irn(i),i=keep(iptro+j), &
                keep(iptro+j+1)-1)
40          CONTINUE
            j1 = j2 + 1
50        CONTINUE
        END IF
        IF (icntl(3)>3) THEN
          WRITE (mp,9010) (jcn(k),k=1,ne)
9010      FORMAT (' JCN (MAP) ='/(6X,10I6))
          WRITE (mp,'(A)') ' Contents of KEEP array'
          WRITE (mp,9020) (keep(i),i=1,m)
          WRITE (mp,9030) (keep(m+i),i=1,n)
9020      FORMAT (' Positions of original rows in the permuted matrix'/(10I6))
9030      FORMAT (' Positions of columns of permuted matrix ','in ', &
            'original matrix '/(10I6))
          WRITE (mp,'(A)') ' Pointer information from KEEP'
          WRITE (mp,9040) (keep(iptrd+j),j=1,n+1)
9040      FORMAT (' IPTRD ='/(8X,10I6))
          WRITE (mp,9050) (keep(iptro+j),j=1,n+1)
9050      FORMAT (' IPTRO ='/(8X,10I6))
          IF (job>1) THEN
            WRITE (mp,9060) (keep(iptrl+j),j=1,n)
9060        FORMAT (' IPTRL ='/(8X,10I6))
            WRITE (mp,9070) (keep(iptru+j),j=1,n)
9070        FORMAT (' IPTRU ='/(8X,10I6))
          END IF
          WRITE (mp,'(A)') ' Block structure information from KEEP'
          WRITE (mp,9080) (keep(nblock+3*jb),jb=1,nb)
          WRITE (mp,9090) (keep(mblock+3*jb),jb=1,nb)
9080      FORMAT (' NBLOCK (order blocks) ='/(8X,10I6))
9090      FORMAT (' MBLOCK (triangular flag and number packed rows) ='/ &
            (8X,10I6))
9100      FORMAT (' KBLOCK (position of beginning of block) ='/(8X,10I6))
9110      FORMAT (' LBLOCK (number of changed columns) ='/(8X,10I6))
          IF (job>1) THEN
            WRITE (mp,9100) (keep(kblock+3*jb),jb=1,nb)
            IF (icntl(8)/=0) WRITE (mp,9110) (keep(lblock+jb),jb=1,nb)
          END IF
        ELSE
          WRITE (mp,9010) (jcn(k),k=1,min(10*onel,ne))
          WRITE (mp,'(A)') ' Contents of KEEP array'
          WRITE (mp,9020) (keep(i),i=1,min(10,m))
          WRITE (mp,9030) (keep(m+i),i=1,min(10,n))
          WRITE (mp,'(A)') ' Pointer information from KEEP'
          WRITE (mp,9040) (keep(iptrd+j),j=1,min(10,n+1))
          WRITE (mp,9050) (keep(iptro+j),j=1,min(10,n+1))
          IF (job>1) THEN
            WRITE (mp,9060) (keep(iptrl+j),j=1,min(10,n))
            WRITE (mp,9070) (keep(iptru+j),j=1,min(10,n))
          END IF
          WRITE (mp,'(A)') ' Block structure information from KEEP'
          WRITE (mp,9080) (keep(nblock+3*jb),jb=1,min(10,nb))
          WRITE (mp,9090) (keep(mblock+3*jb),jb=1,min(10,nb))
          IF (job>1) THEN
            WRITE (mp,9100) (keep(kblock+3*jb),jb=1,min(10,nb))
            IF (icntl(8)/=0) WRITE (mp,9110) (keep(lblock+jb),jb=1,min(10,nb))
          END IF
        END IF
      END IF

! Initialize INFO and RINFO
      info(4) = ne
      info(5) = 0
      info(6) = 0
      rinfo(1) = zero
! Set pivot tolerance
      tol = max(zero,cntl(4))

! Use map to sort the new values into A.
! Mapping into first NEWNE locations in array A
      IF (jcn(1)>0) THEN
        DO 60 k = 1, ne
          a(ne+k) = a(k)
60      CONTINUE
!DIR$ IVDEP
        DO 70 k = 1, ne
          a(jcn(k)) = a(ne+k)
70      CONTINUE
      ELSE
! Duplicates exist
        DO 80 k = 1, ne
          a(ne+k) = a(k)
          a(k) = zero
80      CONTINUE
        a(-jcn(1)) = a(ne+1)
        DO 90 k = 2, ne
          kk = jcn(k)
          a(kk) = a(kk) + a(ne+k)
90      CONTINUE
      END IF

! Call MA50B/BD block by block.
      iqb(1) = 0
      kk = 0
      j2 = 0
      job5 = job
      DO 150 jb = 1, nb
        nc = keep(nblock+3*jb)
        j1 = j2 + 1
        j2 = j1 + nc - 1
        keep(kblock+3*jb) = 0
        IF (keep(mblock+3*jb)<0) THEN
! Action if triangular block
          trisng = .FALSE.
          DO 100 j = j1, j2
            IF (abs(a(keep(iptrd+j)))<=tol) trisng = .TRUE.
            keep(iptrl+j) = 0
            keep(iptru+j) = 0
100       CONTINUE
          IF ( .NOT. trisng) THEN
            info(5) = info(5) + nc
            GO TO 150
          END IF
! Block is singular. Treat as non-triangular.
          IF (job==2) THEN
            info(1) = -7
!!!
! The write statement is switched off by HSL_MA48 call so not tested
!           IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,'(A)') &
!             ' Error return from MA48B/BD with JOB=2 because ', &
!             ' the matrix is incompatible with expectations'
            GO TO 240
          ELSE
            keep(mblock+3*jb) = nc
          END IF
        END IF
! Action if non-triangular block
        DO 145 itry = 1, 2
! The second iteration of this loop is used only if JOB=2,
!   ICNTL(11)=1, and the call to MA50B with JOB=1 fails.
          nr = nc
          IF (nb==1) nr = m
          nzb = keep(iptrd+j2+1) - keep(iptrd+j1)
          IF (icntl(8)/=0) icntl5(6) = keep(lblock+jb)
! Shift the indices and pointers to local values.
          DO 110 k = keep(iptrd+j1), keep(iptrd+j2+1) - 1
            irn(k) = irn(k) - j1 + 1
110       CONTINUE
          k = keep(iptrd+j1) - 1
          DO 115 j = j1, j1 + nc - 1
            keep(iptrd+j) = keep(iptrd+j) - k
115       CONTINUE
          DO 120 j = j1, j1 + nr - 1
            iw(j) = j - j1 + 1
120       CONTINUE
          np = keep(mblock+3*jb)
!!!!!!!!!!!!!!!!!!! MA50BD !!!!!!!!!!!!!!!!!!!!!!!
          CALL ma50bd(nr,nc,nzb,job5,a(k+1),irn(k+1),keep(iptrd+j1),cntl5, &
            icntl5,iw(j1),iqb,np,la-newne-kk,a(newne+kk+1),  &
            la-newne-kk,irn(newne+kk+1), &
            keep(iptrl+j1),keep(iptru+j1),info5,rinfo5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Action if failure in allocate in ma50ad
          if (info5(1) .eq. -2*nr) then
            info(1) = -10
            GO TO 240
          endif
! Restore the indices and pointers
          DO 130 j = j1, j2
            keep(iptrd+j) = keep(iptrd+j) + k
130       CONTINUE
          DO 140 k = keep(iptrd+j1), keep(iptrd+j2+1) - 1
            irn(k) = irn(k) + j1 - 1
140       CONTINUE
! Set warning and error returns
!!!
!XXX This return cannot happen on hsl_ma48 call as "fast" is only
!    allowed if there are no dropped entries.
!         IF (info5(1)==-6) THEN
!           info(1) = -6
!           IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,'(A)') &
!             ' Error return from MA48B/BD with JOB greater than 1', &
!             ' and entries dropped during previous factorization'
!           GO TO 240
!         END IF
          IF (info5(1)<-7) THEN
            IF (icntl(11)==1 .AND. job==2) THEN
              job5 = 1
              IF (lp>0 .AND. icntl(3)>=2) WRITE (lp,'(A,2(A,I4))') &
                ' Warning from MA48B/BD. Switched from JOB=2 to JOB=1', &
                ' in block', jb, ' of ', nb
              info(1) = info(1) + 1
              GO TO 145
            END IF
            info(1) = -7
!!!
! However, this is actually tested in the internal MA48 call in test deck
            IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,'(A)') &
              ' Error return from MA48B/BD with JOB=2 because ', &
              ' the matrix is incompatible with expectations'
            GO TO 240
          ELSE
            GO TO 147
          END IF
145     CONTINUE
147     IF (info5(1)==-3) THEN
          info(1) = -3
          IF (icntl(10)==1) THEN
            keep(kblock+3) = nb
!!!
! However, this is actually tested in the internal MA48 call in test deck
            IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,'(A,2(A,I4))') &
              ' Error return from MA48B/BD because LA is too small.', &
              ' In block', jb, ' of ', nb
            GO TO 240
          END IF
        END IF
        IF (info(1)==-3) THEN
          info(4) = info(4) + info5(4)
          kk = 0
        ELSE
          info(4) = max(info(4),kk+newne+info5(4))
          nrf = irn(newne+kk+2)
! Set pointer to first entry in block JB in A(NEWNE+1), IRN(NEWNE+1)
          keep(kblock+3*jb) = kk + 1
          kk = kk + keep(iptrl+j2) + max((nc-keep(mblock+ &
            3*jb))*(nrf),(nc-keep(mblock+3*jb))+(nrf))
        END IF
! Is matrix rank deficient?
        IF (info5(1)==1) THEN
          IF (info(1)/=-3) info(1) = min(info(1)+2,3*onel)
        END IF
! Accumulate stats
        rinfo(1) = rinfo(1) + rinfo5(1)
        info(5) = info(5) + info5(5)
        info(6) = info(6) + info5(6)
150   CONTINUE

      info(4) = max(ne*2,info(4))
      keep(kblock+3) = nb

      IF (info(1)==-3) THEN
        IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,9180) la, info(4)
        GO TO 240
      END IF

      IF (mp>0 .AND. icntl(3)>2) THEN
        WRITE (mp,'(/A/A,I6/A,3I6/A,F12.1)') ' Leaving MA48B/BD with', &
          ' INFO(1)   = ', info(1), ' INFO(4:6) = ', (info(i),i=4,6), &
          ' RINFO(1)     =', rinfo(1)
        WRITE (mp,'(A)') ' Permuted matrix by blocks'
        kb = nb
        IF (icntl(3)==3) THEN
          WRITE (mp,'(A)') ' Only first column of up to 10 blocks printed'
          kb = min(10,nb)
        END IF
        WRITE (mp,'(A)') ' Diagonal blocks'
        j1 = 1
        DO 170 jb = 1, kb
          j2 = j1 + keep(nblock+3*jb) - 1
          IF (j1<=j2) WRITE (mp,'(A,I6)') ' Block', jb
          j3 = j2
          IF (icntl(3)==3) j3 = min(j1,j2)
          DO 160 j = j1, j3
            WRITE (mp,'(A,I5,(T13,3(1PD12.4,I5)))') ' Column', j, &
              (a(i),irn(i),i=keep(iptrd+j),keep(iptrd+j+1)-1)
160       CONTINUE
          j1 = j2 + 1
170     CONTINUE
        IF (keep(iptro+n+1)>keep(iptrd+n+1)) THEN
          WRITE (mp,'(A)') ' Off-diagonal entries'
          j1 = 1
          DO 190 jb = 1, kb
            j2 = j1 + keep(nblock+3*jb) - 1
            j3 = j2
            IF (icntl(3)==3) j3 = min(j1,j2)
            DO 180 j = j1, j3
              IF (keep(iptro+j+1)>keep(iptro+j)) WRITE (mp, &
                '(A,I5,(T13,3(1P,D12.4,I5)))') ' Column', j, &
                (a(i),irn(i),i=keep(iptro+j),keep(iptro+j+1)-1)
180         CONTINUE
            j1 = j2 + 1
190       CONTINUE
        END IF
        WRITE (mp,'(A)') ' Factorized matrix by blocks'
        j1 = 1
        DO 230 jb = 1, kb
          j2 = j1 + keep(nblock+3*jb) - 1
! Jump if triangular block
          IF (keep(mblock+3*jb)<0) GO TO 220
          nc = j2 - j1 + 1
          nr = nc
          IF (kb==1) nr = m
          WRITE (mp,'(A,I6)') ' Block', jb
          k = newne
          IF (jb>1) k = keep(kblock+3*jb) + newne - 1
          IF (keep(iptrl+j1)>keep(iptru+j1)) WRITE (mp, &
            '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column', j1, ' of L', &
            (a(k+i),irn(k+i),i=keep(iptru+j1)+1,keep(iptrl+j1))
          IF (icntl(3)==3) GO TO 210
          DO 200 j = j1 + 1, j2
            IF (keep(iptru+j)>keep(iptrl+j-1)) WRITE (mp, &
              '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column', j, ' of U', &
              (a(k+i),irn(k+i),i=keep(iptrl+j-1)+1,keep(iptru+j))
            IF (keep(iptru+j)<keep(iptrl+j)) WRITE (mp, &
              '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column', j, ' of L', &
              (a(k+i),irn(k+i),i=keep(iptru+j)+1,keep(iptrl+j))
200       CONTINUE
! Full blocks
210       WRITE (mp,'(A)') ' Full block'
          WRITE (mp,'(A)') ' Row indices'
          nrf = irn(k+2)
          k = k + keep(iptrl+j2)
          IF (icntl(3)>3) THEN
            WRITE (mp,9120) (irn(k+i),i=1,nrf)
            WRITE (mp,'(A)') ' Column pivoting information'
            WRITE (mp,9120) (irn(k+i),i=nrf+1,nrf+nc-keep(mblock+3*jb))
            WRITE (mp,'(A)') ' Reals by columns'
            WRITE (mp,9130) (a(k+i),i=1,(nrf)*(nc-keep(mblock+3*jb)))
9120        FORMAT (10I6)
9130        FORMAT (1P,5D12.4)
          ELSE
            WRITE (mp,9120) (irn(k+i),i=1,min(10*onel,nrf))
            WRITE (mp,'(A)') ' Column pivoting information'
            WRITE (mp,9120) (irn(k+i),i=nrf+1,nrf+min(10*onel,nc-keep(mblock+ &
              3*jb)))
            WRITE (mp,'(A)') ' Reals by columns'
            WRITE (mp,9130) &
               (a(k+i),i=1,min(10*onel,(nrf)*(nc-keep(mblock+3*jb))))
          END IF
220       j1 = j2 + 1
230     CONTINUE
        IF (job==1 .OR. job==3) THEN
          WRITE (mp,'(A)') ' Contents of KEEP array'
          IF (icntl(3)>3) THEN
            WRITE (mp,9020) (keep(i),i=1,m)
            WRITE (mp,9030) (keep(m+i),i=1,n)
            WRITE (mp,'(A)') ' Pointer information from KEEP'
            WRITE (mp,9140) (keep(iptrl+j),j=1,n)
9140        FORMAT (' IPTRL ='/(8X,10I6))
            WRITE (mp,9150) (keep(iptru+j),j=1,n)
9150        FORMAT (' IPTRU ='/(8X,10I6))
          ELSE
            WRITE (mp,9020) (keep(i),i=1,min(10,m))
            WRITE (mp,9030) (keep(m+i),i=1,min(10,n))
            WRITE (mp,'(A)') ' Pointer information from KEEP'
            WRITE (mp,9140) (keep(iptrl+j),j=1,min(10,n))
            WRITE (mp,9150) (keep(iptru+j),j=1,min(10,n))
          END IF
          WRITE (mp,'(A)') ' Block structure information from KEEP'
          WRITE (mp,9100) (keep(kblock+3*jb),jb=1,kb)
        END IF
      END IF

240   deallocate(w,iw,stat=stat)
      RETURN

9160  FORMAT (' Error return from MA48B/BD because M =',I10,' and N =',I10)
9170  FORMAT (' Error return from MA48B/BD because NE =',I10)
9180  FORMAT (' Error return from MA48B/BD because LA is',I10/' and mu', &
        'st be at least',I10)
9190  FORMAT (' Error return from MA48B/BD because ','JOB = ',I10)
    END  subroutine ma48bd


    SUBROUTINE ma48cd(m,n,trans,job,la,a,irn,keep,cntl,icntl,rhs,x,error, &
        info)

! Solve linear system, using data provided by MA48B/BD.

!     .. Arguments ..
      INTEGER m, n, job
      LOGICAL trans
      integer(long) :: la
      DOUBLE PRECISION a(la)
      INTEGER(long) irn(la)
      integer(long) :: keep(*)
      DOUBLE PRECISION cntl(10)
      INTEGER icntl(20)
      DOUBLE PRECISION rhs(*), x(*), error(3)
      real(wp), allocatable :: w(:)
      integer, allocatable :: iw(:)
      INTEGER(long) :: info(20)

! M must be set by the user to the number of rows in the matrix.
!      It is not altered by the subroutine. Restriction: M > 0.
! N must be set by the user to the number of columns in the matrix.
!      It is not altered by the subroutine. Restriction: N > 0.
! TRANS must be set by the user to indicate whether coefficient matrix
!      is A (TRANS=.FALSE.) or A transpose (TRANS=.TRUE.).
!      It is not altered by the subroutine.
! JOB  must be set by the user to control the solution phase.
!      Restriction: 1 <= JOB <= 4
!      Possible values are:
!       =1 Returns solution only.
!       =2 Returns solution and estimate of backward error.
!       =3 Performs iterative refinement and returns estimate of
!          backward error.
!       =4 As for 3 plus an estimate of the error in the solution.
!      It is not altered by the subroutine.
! LA must be set by the user to the size of arrays A and IRN.
!      It is not altered by the subroutine.
! A    must be set left unchanged from the last call to MA48B/BD.
!      It holds the original matrix in permuted form and the
!      factorization of the block diagonal part (excluding
!      any triangular blocks). It is not altered by the subroutine.
! IRN  must be set left unchanged from the last call to MA48B/BD.
!      It holds the row indices of the factorization in A and
!      the row indices of the original matrix in permuted form.
!      It is not altered by the subroutine.
! KEEP must be as on return from MA48B/BD.
!      It is not altered by the subroutine.
! CNTL  must be set by the user as follows and is not altered.
!     CNTL(I), I=1,4 not accessed by subroutine.
!     CNTL(5) is used in convergence test for termination of iterative
!     refinement.  Iteration stops if successive omegas do not decrease
!     by at least this factor.
! ICNTL must be set by the user as follows and is not altered.
!     ICNTL(1)  must be set to the stream number for error messages.
!       A value less than 1 suppresses output.
!     ICNTL(2) must be set to the stream number for diagnostic output.
!       A value less than 1 suppresses output.
!     ICNTL(3) must be set to control the amount of output:
!       0 None.
!       1 Error messages only.
!       2 Error and warning messages.
!       3 As 2, plus scalar parameters and a few entries of array
!         parameters on entry and exit.
!       4 As 2, plus all parameters on entry and exit.
!     ICNTL(I), I=4,8 not accessed by subroutine
!     ICNTL(9) is maximum number of iterations allowed in iterative
!         refinement.
! RHS  must be set by the user to contain the right-hand side of the
!      equations to be solved. If TRANS is .FALSE., RHS has length M;
!      otherwise, it has length N. It is used as workspace.
! X    need not be set by the user.  On exit, it will contain the
!      solution vector. If TRANS is .FALSE., X has length N; otherwise,
!      it has length M.
! ERROR is a real array of length 3 that need not be set by the user.
!      On exit with JOB >=2 ERROR(1) and ERROR(2) give estimates of
!      backward errors OMEGA1 and OMEGA2.  On exit with JOB = 4,
!      ERROR(3) holds an estimate of the infinity norm in the relative
!      error in the solution.
! IW   is a workarray.
! INFO need not be set on entry. On exit, it holds the following:
!    INFO(1) A nonzero value will indicate an error return. Possible
!      nonzero values are:
!      -1  M or N < 1
!      -6  JOB out of range
!      -8  Nonconvergence of iterative refinement
!      -9  Failure in MC71A/AD

!     .. Local constants ..
      real(wp) zero
      PARAMETER (zero=0.0_wp)

!     .. Local variables ..
      real(wp) cond(2), ctau, dxmax
      integer icntl5(20), kase, lp, mp, nb, nc, neq, nrf, nvar, iblock
      INTEGER(long) :: i, iptrd, iptrl, iptro, iptru, iqb(1), j, jb, jj, &
        j1, j2, j3, k, kb, kblock, kk
      INTEGER :: keep71(5)
      LOGICAL lcond(2)
      INTEGER(long) :: mblock, nblock, ne
      real(wp) oldomg(2), omega(2), om1, om2, tau
! COND  Arioli, Demmel, and Duff condition number
! CTAU is a real variable used to control the splitting of the equations
!     into two categories. This is constant in Arioli, Demmel, and Duff
!     theory.
! DXMAX max-norm of current solution estimate
! I     row index and DO loop variable
! ICNTL5 passed to MA50C/CD to correspond to dummy argument ICNTL.
! IPTRD displacement in KEEP. See comment on KEEP.
! IPTRL displacement in KEEP. See comment on KEEP.
! IPTRO displacement in KEEP. See comment on KEEP.
! IPTRU displacement in KEEP. See comment on KEEP.
! IQB   is passed to MA50C/CD to indicate that no column permutation
!       is used.
! J     column index
! JB    is current block.
! JJ    running index for column.
! J1    is index of first column in block.
! J2    is index of last column in block.
! J3    is index of last column in block used in printing diagnostics.
! K     iteration counter and DO loop variable
! KASE  control for norm calculating routine MC71A/AD
! KB    used in prints to hold number of blocks or less if ICNTL(3)
!       equal to 3.
! KBLOCK displacement in KEEP. See comment on KEEP.
! KEEP71 workspace to preserve MC71A/AD's locals between calls
! KK    DO loop variable
! LCOND LCOND(k) is set to .TRUE. if there are equations in category
!       k, k=1,2
! LP Unit for error messages.
! MBLOCK displacement in KEEP. See comment on KEEP.
! MP Unit for diagnostic messages.
! NB    number of diagonal blocks
! NBLOCK displacement in KEEP. See comment on KEEP.
! NC    is number of columns in packed form for current block.
! NE    is set to the displacement in A/ICN for the beginning of
!       information on the factors.
! NEQ   number of equations in system
! NRF   is number of rows in full form for current block.
! NVAR  number of variables in system
! OLDOMG value of omega from previous iteration. Kept in case it is
!       better.
! OMEGA backward error estimates.
! OM1   value of OMEGA(1)+OMEGA(2) from previous iteration
! OM2   value of OMEGA(1)+OMEGA(2) from current iteration
! TAU   from Arioli, Demmel, and Duff .. used to divide equations

!     .. Externals ..
      EXTERNAL mc71ad
      real(wp) eps
! EPS is the largest real such that 1+EPS is equal to 1.
! MA48D/DD solves block triangular system
! MA50C/CD solves non-blocked system
! MC71A/AD estimates norm of matrix
      INTRINSIC abs, max

! The amount allocated to w can be dependent on input parameters
      allocate(w(4*max(m,n)),iw(max(m,n)),stat=stat)
      if (stat .ne. 0) then
        info(1) = -10
        return
      endif
! Set streams for errors and warnings
      lp = icntl(1)
      mp = icntl(2)

! Simple data checks
      IF (n<=0 .OR. m<=0) THEN
        info(1) = -1
        IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,9140) m, n
        GO TO 380
      END IF
      IF (job>4 .OR. job<1) THEN
        info(1) = -6
        IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,9150) job
        GO TO 380
      END IF
      info(1) = 0

! Partition KEEP and A
      iptrl = m + n
      iptru = iptrl + n
      iptrd = iptru + n
      iptro = iptrd + n + 1
      nblock = iptro + n - 1
      mblock = nblock + 1
      kblock = mblock + 1
      nb = keep(kblock+3)

      ne = keep(iptro+n+1) - 1

      omega(1) = zero
      omega(2) = zero
      error(3) = zero
! Initialize EPS
      eps = epsilon(eps)
! CTAU ... 1000 eps (approx)
      ctau = 1000.*eps

      iqb(1) = 0

      DO 10 i = 1, 7
        icntl5(i) = 0
10    CONTINUE
! Set control for use of BLAS
      icntl5(5) = icntl(5)

      IF (mp>0 .AND. icntl(3)>2) THEN
        WRITE (mp,'(/A/3(A,I8),A,I2/A,L2,A,I7/A,1P,E12.4/A,3I6,2(/A,I6))') &
          ' Entering MA48C/CD with', ' M =', m, '     N =', n, '     LA =', &
          la, '      JOB =', job, '   TRANS =', trans, &
          '      No. of blocks =', nb, '   CNTL(5)    = ', cntl(5), &
          '   ICNTL(1:3) = ', (icntl(i),i=1,3), '   ICNTL(5)   = ', icntl(5), &
          '   ICNTL(9)   = ', icntl(9)
        WRITE (mp,'(A)') ' Permuted matrix by blocks'
        kb = nb
        IF (icntl(3)==3) THEN
          WRITE (mp,'(A)') ' Only first column of up to 10 blocks printed'
          kb = min(10,nb)
        END IF
        WRITE (mp,'(A)') ' Diagonal blocks'
        j1 = 1
        DO 30 jb = 1, kb
          j2 = j1 + keep(nblock+3*jb) - 1
          IF (j1<=j2) WRITE (mp,'(A,I6)') ' Block', jb
          j3 = j2
          IF (icntl(3)==3) j3 = min(j1,j2)
          DO 20 j = j1, j3
            WRITE (mp,'(A,I5,(T13,3(1P,E12.4,I5)))') ' Column', j, &
              (a(i),irn(i),i=keep(iptrd+j),keep(iptrd+j+1)-1)
20        CONTINUE
          j1 = j2 + 1
30      CONTINUE
        IF (keep(iptro+n+1)>keep(iptrd+n+1)) THEN
          WRITE (mp,'(A)') ' Off-diagonal entries'
          j1 = 1
          DO 50 jb = 1, kb
            j2 = j1 + keep(nblock+3*jb) - 1
            j3 = j2
            IF (icntl(3)==3) j3 = min(j1,j2)
            DO 40 j = j1, j3
              IF (keep(iptro+j+1)>keep(iptro+j)) WRITE (mp, &
                '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column', j, &
                (a(i),irn(i),i=keep(iptro+j),keep(iptro+j+1)-1)
40          CONTINUE
            j1 = j2 + 1
50        CONTINUE
        END IF
        WRITE (mp,'(A)') ' Factorized matrix by blocks'
        j1 = 1
        DO 90 jb = 1, kb
          j2 = j1 + keep(nblock+3*jb) - 1
! Jump if triangular block
          IF (keep(mblock+3*jb)<0) GO TO 80
          nc = j2 - j1 + 1
          WRITE (mp,'(A,I6)') ' Block', jb
          k = ne
          IF (jb>1) k = keep(kblock+3*jb) + ne - 1
          IF (keep(iptrl+j1)>keep(iptru+j1)) WRITE (mp, &
            '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column', j1, ' of L', &
            (a(k+i),irn(k+i),i=keep(iptru+j1)+1,keep(iptrl+j1))
          IF (icntl(3)==3) GO TO 70
          DO 60 j = j1 + 1, j2
            IF (keep(iptru+j)>keep(iptrl+j-1)) WRITE (mp, &
              '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column', j, ' of U', &
              (a(k+i),irn(k+i),i=keep(iptrl+j-1)+1,keep(iptru+j))
            IF (keep(iptru+j)<keep(iptrl+j)) WRITE (mp, &
              '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column', j, ' of L', &
              (a(k+i),irn(k+i),i=keep(iptru+j)+1,keep(iptrl+j))
60        CONTINUE
! Full blocks
70        WRITE (mp,'(A)') ' Full block'
          WRITE (mp,'(A)') ' Row indices'
          nrf = irn(k+2)
          k = k + keep(iptrl+j2)
          IF (icntl(3)>3) THEN
            WRITE (mp,9000) (irn(k+i),i=1,nrf)
            WRITE (mp,'(A)') ' Column pivoting information'
            WRITE (mp,9000) (irn(k+i),i=nrf+1,nrf+nc-keep(mblock+3*jb))
            WRITE (mp,'(A)') ' Reals by columns'
            WRITE (mp,9010) (a(k+i),i=1,(nrf)*(nc-keep(mblock+3*jb)))
9000        FORMAT (10I6)
9010        FORMAT (1P,5D12.4)
          ELSE
            WRITE (mp,9000) (irn(k+i),i=1,min(10,nrf))
            WRITE (mp,'(A)') ' Column pivoting information'
            WRITE (mp,9000) (irn(k+i),i=nrf+1,nrf+min(10*onel,nc-keep(mblock+ &
              3*jb)))
            WRITE (mp,'(A)') ' Reals by columns'
            WRITE (mp,9010) &
               (a(k+i),i=1,min(10*onel,(nrf)*(nc-keep(mblock+3*jb))))
          END IF
80        j1 = j2 + 1
90      CONTINUE
        IF (icntl(3)>3) THEN
          WRITE (mp,'(A)') ' Contents of KEEP array'
          WRITE (mp,9020) (keep(i),i=1,m)
          WRITE (mp,9030) (keep(m+i),i=1,n)
9020      FORMAT (' Positions of original rows in the permuted matrix'/(10I6))
9030      FORMAT (' Positions of columns of permuted matrix ','in or','ig', &
            'inal matrix '/(10I6))
          WRITE (mp,'(A)') ' Pointer information from KEEP'
          WRITE (mp,9040) (keep(iptrd+j),j=1,n+1)
9040      FORMAT (' IPTRD ='/(8X,10I6))
          WRITE (mp,9050) (keep(iptro+j),j=1,n+1)
9050      FORMAT (' IPTRO ='/(8X,10I6))
          WRITE (mp,9060) (keep(iptrl+j),j=1,n)
9060      FORMAT (' IPTRL ='/(8X,10I6))
          WRITE (mp,9070) (keep(iptru+j),j=1,n)
9070      FORMAT (' IPTRU ='/(8X,10I6))
          WRITE (mp,'(A)') ' Block structure information from KEEP'
          WRITE (mp,9080) (keep(nblock+3*jb),jb=1,nb)
          WRITE (mp,9090) (keep(mblock+3*jb),jb=1,nb)
          WRITE (mp,9100) (keep(kblock+3*jb),jb=1,nb)
9080      FORMAT (' NBLOCK (order blocks) ='/(8X,10I6))
9090      FORMAT (' MBLOCK (triangular flag and number packed rows) ='/ &
            (8X,10I6))
9100      FORMAT (' KBLOCK (position of beginning of block) ='/(8X,10I6))
          IF (trans) THEN
            WRITE (mp,9110) (rhs(i),i=1,n)
9110        FORMAT (' RHS =  ',1P,5D12.4/(8X,5D12.4))
          ELSE
            WRITE (mp,9110) (rhs(i),i=1,m)
          END IF
        ELSE
          WRITE (mp,'(A)') ' Contents of KEEP array'
          WRITE (mp,9020) (keep(i),i=1,min(10,m))
          WRITE (mp,9030) (keep(m+i),i=1,min(10,n))
          WRITE (mp,'(A)') ' Pointer information from KEEP'
          WRITE (mp,9040) (keep(iptrd+k),k=1,min(10,n+1))
          WRITE (mp,9050) (keep(iptro+k),k=1,min(10,n+1))
          WRITE (mp,9060) (keep(iptrl+j),j=1,min(10,n))
          WRITE (mp,9070) (keep(iptru+j),j=1,min(10,n))
          WRITE (mp,'(A)') ' Block structure information from KEEP'
          WRITE (mp,9080) (keep(nblock+3*jb),jb=1,kb)
          WRITE (mp,9090) (keep(mblock+3*jb),jb=1,kb)
          WRITE (mp,9100) (keep(kblock+3*jb),jb=1,kb)
          IF (trans) THEN
            WRITE (mp,9110) (rhs(i),i=1,min(10,n))
          ELSE
            WRITE (mp,9110) (rhs(i),i=1,min(10,m))
          END IF
        END IF
      END IF

! Apply global permutation to incoming right-hand side
!     and set NEQ and NVAR
      IF (trans) THEN
        neq = n
        nvar = m
        DO 100 i = 1, neq
          w(i) = rhs(keep(m+i))
100     CONTINUE
      ELSE
        neq = m
        nvar = n
        DO 110 i = 1, neq
          w(keep(i)) = rhs(i)
110     CONTINUE
      END IF


! Straight solution requested with no iterative refinement or error
!     estimate.
      IF (job==1) THEN

! Solve system using MA48D/DD or MA50C/CD.
        IF (nb==1 .AND. keep(mblock+3)>=0) THEN
          iblock = keep(mblock+3)
          CALL ma50cd(m,n,icntl5,iqb,iblock,trans,la-ne,a(ne+1), &
            irn(ne+1),keep(iptrl+1),keep(iptru+1),w,x,w(neq+1))
        ELSE
          CALL ma48dd(n,ne,la-ne,a(ne+1),a,irn(ne+1),irn,keep(iptrd+1), &
            keep(iptro+1),nb,keep(nblock+3),keep(iptrl+1),keep(iptru+1),w,x, &
            trans,icntl5,info)
          if (info(1) .eq. -10) go to 380
        END IF
        GO TO 340
      END IF

! Prepare for iterative refinement.

! Set initial estimate of solution to zero
      DO 120 i = 1, nvar
        x(i) = zero
120   CONTINUE
! Save permuted right-hand side for residual calculation
      DO 130 i = 1, neq
        rhs(i) = w(i)
130   CONTINUE


! Iterative refinement loop
! Initialize OM1 in case of problems with optimizing compiler
      om1 = zero
      DO 260 k = 1, icntl(9)

! Solve system using MA48D/DD or MA50C/CD.
        IF (nb==1 .AND. keep(mblock+3)>=0) THEN
          iblock = keep(mblock+3)
          CALL ma50cd(m,n,icntl5,iqb,iblock,trans,la-ne,a(ne+1), &
            irn(ne+1),keep(iptrl+1),keep(iptru+1),w,w(neq+1),w(m+n+1))
        ELSE
          CALL ma48dd(n,ne,la-ne,a(ne+1),a,irn(ne+1),irn,keep(iptrd+1), &
            keep(iptro+1),nb,keep(nblock+3),keep(iptrl+1),keep(iptru+1),w, &
            w(neq+1),trans,icntl5,info)
          if (info(1) .eq. -10) go to 380
        END IF
! Update solution
        DO 140 i = 1, nvar
          x(i) = x(i) + w(neq+i)
140     CONTINUE

! Calculate residual using information in A,IRN

        DO 150 i = 1, neq
! Residual  .. b-Ax
          w(i) = rhs(i)
! |A||x|
          w(neq+i) = zero
! Sum |a  |, j=1,N (= ||A  ||        )
!       ij               i.  infinity
          w(2*neq+i) = zero
150     CONTINUE
        IF (trans) THEN
          DO 180 j = 1, n
!DIR$ IVDEP
            DO 160 jj = keep(iptrd+j), keep(iptrd+j+1) - 1
              i = irn(jj)
              w(j) = w(j) - a(jj)*x(i)
              w(neq+j) = w(neq+j) + abs(a(jj)*x(i))
              w(2*neq+j) = w(2*neq+j) + abs(a(jj))
160         CONTINUE
!DIR$ IVDEP
            DO 170 jj = keep(iptro+j), keep(iptro+j+1) - 1
              i = irn(jj)
              w(j) = w(j) - a(jj)*x(i)
              w(neq+j) = w(neq+j) + abs(a(jj)*x(i))
              w(2*neq+j) = w(2*neq+j) + abs(a(jj))
170         CONTINUE
180       CONTINUE
        ELSE
          DO 210 j = 1, n
!DIR$ IVDEP
            DO 190 jj = keep(iptrd+j), keep(iptrd+j+1) - 1
              i = irn(jj)
              w(i) = w(i) - a(jj)*x(j)
              w(neq+i) = w(neq+i) + abs(a(jj)*x(j))
              w(2*neq+i) = w(2*neq+i) + abs(a(jj))
190         CONTINUE
!DIR$ IVDEP
            DO 200 jj = keep(iptro+j), keep(iptro+j+1) - 1
              i = irn(jj)
              w(i) = w(i) - a(jj)*x(j)
              w(neq+i) = w(neq+i) + abs(a(jj)*x(j))
              w(2*neq+i) = w(2*neq+i) + abs(a(jj))
200         CONTINUE
210       CONTINUE
        END IF
! Calculate max-norm of solution
        dxmax = zero
        DO 220 i = 1, nvar
          dxmax = max(dxmax,abs(x(i)))
220     CONTINUE
! Calculate also omega(1) and omega(2)
! tau is (||A  ||         ||x||   + |b| )*n*1000*epsilon
!            i.  infinity      max     i
        omega(1) = zero
        omega(2) = zero
        DO 230 i = 1, neq
          tau = (w(2*neq+i)*dxmax+abs(rhs(i)))*nvar*ctau
          IF ((w(neq+i)+abs(rhs(i)))>tau) THEN
! |Ax-b| /(|A||x| + |b|)
!       i               i
            omega(1) = max(omega(1),abs(w(i))/(w(neq+i)+abs(rhs(i))))
            iw(i) = 1
          ELSE
! TAU will be zero if all zero row in A, for example
            IF (tau>zero) THEN
! |Ax-b| /(|A||x| + ||A  ||        ||x||   )
!       i        i     i.  infinity     max
              omega(2) = max(omega(2),abs(w(i))/(w(neq+i)+w(2*neq+i)*dxmax))
            END IF
            iw(i) = 2
          END IF
230     CONTINUE

! Exit if iterative refinement not being performed
        IF (job==2) GO TO 340

!  Stop the calculations if the backward error is small

        om2 = omega(1) + omega(2)
! Jump if converged
! Statement changed because IBM SP held quantities in registers
!        IF ((OM2+ONE).LE.ONE) GO TO 270
        IF (om2<=eps) GO TO 270

!  Check the convergence.

        IF (k>1 .AND. om2>om1*cntl(5)) THEN
!  Stop if insufficient decrease in omega.
          IF (om2>om1) THEN
! Previous estimate was better ... reinstate it.
            omega(1) = oldomg(1)
            omega(2) = oldomg(2)
            DO 240 i = 1, nvar
              x(i) = w(3*neq+i)
240         CONTINUE
          END IF
          GO TO 270
        END IF
! Hold current estimate in case needed later
        DO 250 i = 1, nvar
          w(3*neq+i) = x(i)
250     CONTINUE
        oldomg(1) = omega(1)
        oldomg(2) = omega(2)
        om1 = om2
260   CONTINUE
! End of iterative refinement loop.
!!!
! Can only happen when IR fails ... have eyeballed this and
! tested it in practice.
      info(1) = -8
      IF (lp>0 .AND. icntl(3)>=1) WRITE (lp,9170) info(1), icntl(9)
      GO TO 340

270   IF (job<=3) GO TO 340
      IF (m/=n) GO TO 340

! Calculate condition numbers and estimate of the error.

!  Condition numbers obtained through use of norm estimation
!     routine MC71A/AD.

!  Initializations

      lcond(1) = .FALSE.
      lcond(2) = .FALSE.
      DO 280 i = 1, neq
        IF (iw(i)==1) THEN
          w(i) = w(neq+i) + abs(rhs(i))
! |A||x| + |b|
          w(neq+i) = zero
          lcond(1) = .TRUE.
        ELSE
! |A||x| + ||A  ||        ||x||
!             i.  infinity     max

          w(neq+i) = w(neq+i) + w(2*neq+i)*dxmax
          w(i) = zero
          lcond(2) = .TRUE.
        END IF
280   CONTINUE

!  Compute the estimate of COND

      kase = 0
      DO 330 k = 1, 2
        IF (lcond(k)) THEN
! MC71A/AD has its own built in limit to the number of iterations
!    allowed. It is this limit that will be used to terminate the
!    following loop.
          DO 310 kk = 1, 40
! MC71A/AD calculates norm of matrix
! We are calculating the infinity norm of INV(A).W
            CALL mc71ad(n,kase,w(3*neq+1),cond(k),rhs,iw,keep71)

!  KASE = 0........ Computation completed
!  KASE = 1........ W * INV(TRANSPOSE(A)) * Y
!  KASE = 2........ INV(A) * W * Y
!                   W is W/W(NEQ+1) .. Y is W(3*NEQ+1)

            IF (kase==0) GO TO 320
            IF (kase==1) THEN
! Solve system using MA48D/DD or MA50C/CD.
              IF (nb==1 .AND. keep(mblock+3)>=0) THEN
          iblock = keep(mblock+3)
                CALL ma50cd(m,n,icntl5,iqb,iblock, .NOT. trans,la-ne, &
                  a(ne+1),irn(ne+1),keep(iptrl+1),keep(iptru+1),w(3*neq+1), &
                  w(2*neq+1),rhs)
              ELSE
                CALL ma48dd(n,ne,la-ne,a(ne+1),a,irn(ne+1),irn,keep(iptrd+1), &
                  keep(iptro+1),nb,keep(nblock+3),keep(iptrl+1),keep(iptru+1), &
                  w(3*neq+1),w(2*neq+1), .NOT. trans,icntl5,info)
                 if (info(1) .eq. -10) go to 380
              END IF

              DO 290 i = 1, m
                w(3*neq+i) = w((k-1)*neq+i)*w(2*neq+i)
290           CONTINUE
            END IF
            IF (kase==2) THEN
              DO 300 i = 1, n
                w(2*neq+i) = w((k-1)*neq+i)*w(3*neq+i)
300           CONTINUE
! Solve system using MA48D/DD or MA50C/CD.
              IF (nb==1 .AND. keep(mblock+3)>=0) THEN
          iblock = keep(mblock+3)
                CALL ma50cd(m,n,icntl5,iqb,iblock,trans,la-ne,a(ne+1), &
                  irn(ne+1),keep(iptrl+1),keep(iptru+1),w(2*neq+1),w(3*neq+1), &
                  rhs)
              ELSE
                CALL ma48dd(n,ne,la-ne,a(ne+1),a,irn(ne+1),irn,keep(iptrd+1), &
                  keep(iptro+1),nb,keep(nblock+3),keep(iptrl+1),keep(iptru+1), &
                  w(2*neq+1),w(3*neq+1),trans,icntl5,info)
                 if (info(1) .eq. -10) go to 380
              END IF
            END IF
310       CONTINUE
!!!
! Can only happen when MC71 does not converge in 40 iterations.
! Have eyeballed this but never invoked it.
          info(1) = -9
          IF (lp/=0 .AND. icntl(3)>=1) WRITE (lp,9160)
          GO TO 340
320       IF (dxmax>zero) cond(k) = cond(k)/dxmax
          error(3) = error(3) + omega(k)*cond(k)
        END IF
330   CONTINUE

! Permute solution vector
340   DO 350 i = 1, nvar
        w(i) = x(i)
350   CONTINUE
      IF ( .NOT. trans) THEN
        DO 360 i = 1, nvar
          x(keep(m+i)) = w(i)
360     CONTINUE
      ELSE
        DO 370 i = 1, nvar
          x(i) = w(keep(i))
370     CONTINUE
      END IF
      IF (job>=2) THEN
        error(1) = omega(1)
        error(2) = omega(2)
      END IF

      IF (mp>0 .AND. icntl(3)>2) THEN
        WRITE (mp,'(/A,I6)') ' Leaving MA48C/CD with INFO(1) =', info(1)
        IF (job>1) THEN
          k = 2
          IF (job==4 .AND. info(1)/=-9) k = 3
          WRITE (mp,9120) (error(i),i=1,k)
9120      FORMAT (' ERROR =',1P,3D12.4)
        END IF
        IF (icntl(3)>3) THEN
          WRITE (mp,9130) (x(i),i=1,nvar)
9130      FORMAT (' X =    ',1P,5D12.4:/(8X,5D12.4))
        ELSE
          WRITE (mp,9130) (x(i),i=1,min(10,nvar))
        END IF
      END IF
380   deallocate(w,iw,stat=stat)
      RETURN

9140  FORMAT (' Error return from MA48C/CD because M =',I10,' and N =',I10)
9150  FORMAT (' Error return from MA48C/CD because ','JOB = ',I10)
9160  FORMAT (' Error return from MA48C/CD because of ','error in MC71', &
        'A/AD'/' ERROR(3) not calculated')
9170  FORMAT (' Error return from MA48C/CD because of ','nonconvergenc', &
        'e of iterative refinement'/' Error INFO(1) = ',I2,'  wit','h ICNTL', &
        '(9) = ',I10)
    END subroutine ma48cd


    SUBROUTINE ma48dd(n,ne,la,a,aa,irn,irna,iptrd,iptro,nb,iblock,iptrl,iptru, &
        rhs,x,trans,icntl5,info)
! Solve linear system, using data provided by MA48B/BD.

!     .. Arguments ..
      INTEGER n, nb
      integer(long) :: ne, la
      integer (long) :: info(20)
      real(wp) a(la), aa(ne)
      INTEGER(long) :: irn(la), irna(ne)
      integer(long) :: iptrd(n+1), iptro(n+1)
      integer(long) :: iblock(3,nb)
      integer(long) :: iptrl(n), iptru(n)
      real(wp) rhs(n), x(n)
      LOGICAL trans
      INTEGER icntl5(20)
      real(wp), allocatable :: w(:)

! N must be set by the user to the number of columns in the matrix.
!      It is not altered by the subroutine.
! NE must be set by the user to the size of arrays AA and IRNA.
!      It is not altered by the subroutine.
! LA must be set by the user to the size of arrays A and IRN.
!      It is not altered by the subroutine.
! A    must be set left unchanged from the last call to MA48B/BD.
!      It holds the factorization of the block diagonal part
!      (excluding triangular blocks).
!      It is not altered by the subroutine.
! AA   must be set left unchanged from the last call to MA48B/BD.
!      It holds the original matrix in permuted form.
!      It is not altered by the subroutine.
! IRN  must be set left unchanged from the last call to MA48B/BD.
!      It holds the row indices of the factorization in A.
!      It is not altered by the subroutine.
! IRNA must be set left unchanged from the last call to MA48B/BD.
!      It holds the row indices of the original matrix in permuted form.
!      It is not altered by the subroutine.
! IPTRD must be as on return from MA48B/BD. IPTRD(j) holds the position
!      in IRNA of the start of the block diagonal part of column j.
!      IPTRD(n+1) holds the position immediately after the end of
!      column n.
!      It is not altered by the subroutine.
! IPTRO must be as on return from MA48A/AD. IPTRO(j) holds the position
!      in IRNA of the start of the off block diagonal part of column j.
!      IPTRO(n+1) holds the position immediately after the end of
!      column n.
!      It is not altered by the subroutine.
! NB must be set by the user to the second dimension of IBLOCK.
!      It is not altered by the subroutine.
! IBLOCK must be as on return from MA48B/BD.  IBLOCK(1,k) holds the
!      order of diagonal block k, k=1,2,...  IBLOCK(2,k) holds the
!      number of rows held in packed storage when
!      processing block k, or a negative value for triangular blocks.
!      IBLOCK(3,1) holds the number of blocks and IBLOCK(3,k) holds the
!      position in the factorization of the start of block k, k =
!      2,3,...  IBLOCK is not altered by the subroutine.
! IPTRL must be unchanged since the last call to MA48B/BD.
!      It holds pointers to the
!      columns of the lower-triangular part of the factorization.
!      It is not altered by the subroutine.
! IPTRU must be unchanged since the last call to MA48B/BD.
!      It holds pointers to the
!      columns of the upper-triangular part of the factorization.
!      It is not altered by the subroutine.
! RHS  must be set by the user to contain the right-hand side of
!      the equations to be solved.
!      It is used as workspace.
! X    need not be set by the user.  On exit, it will contain the
!      solution vector.
! TRANS must be set by the user to indicate whether coefficient matrix
!      is A (TRANS=.FALSE.) or A transpose (TRANS=.TRUE.).
!      It is not altered by the subroutine.
! ICNTL5 passed to MA50C/CD to correspond to dummy argument ICNTL.
! W    is a workarray.


!     .. Local variables ..
      integer :: nc, np
      integer(long) :: i, iqb(1), j, jb, jj, j1, k1, k2, numb
! I     row index
! IFLAG error flag from MA50CD .. will never be set
! IQB   is passed to MA50C/CD to indicate that no column permutation
!       is used.
! J     column index
! JB    block index.
! JJ    running index for column.
! J1    position of beginning of block
! K1    index of first column in block
! K2    index of last column in block
! NC    number of columns in block
! NUMB  number of diagonal blocks

!     .. Externals ..

      allocate(w(n),stat=stat)
      if (stat .ne. 0) then
        info(1) = -10
        return
      endif

      iqb(1) = 0
      numb = iblock(3,1)
      IF ( .NOT. trans) THEN
! Solve system using block structure information and calls to MA50C
!     for each diagonal block.
! System is block upper triangular.
        k1 = n + 1
        DO 50 jb = numb, 1, -1
          nc = iblock(1,jb)
          k2 = k1 - 1
          k1 = k1 - nc
          IF (iblock(2,jb)<0) THEN
! Process triangular block
            DO 20 j = k2, k1, -1
              x(j) = rhs(j)/aa(iptrd(j))
!DIR$ IVDEP
              DO 10 jj = iptrd(j) + 1, iptrd(j+1) - 1
                i = irna(jj)
                rhs(i) = rhs(i) - aa(jj)*x(j)
10            CONTINUE
20          CONTINUE
          ELSE
            j1 = 1
            IF (jb>1) j1 = iblock(3,jb)
            np = iblock(2,jb)
            CALL ma50cd(nc,nc,icntl5,iqb,np,trans,la+1-j1,a(j1), &
              irn(j1),iptrl(k1),iptru(k1),rhs(k1),x(k1),w)
          END IF
          IF (jb==1) GO TO 50
! Substitution using off-diagonal block
          DO 40 j = k1, k2
!DIR$ IVDEP
            DO 30 jj = iptro(j), iptro(j+1) - 1
              i = irna(jj)
              rhs(i) = rhs(i) - aa(jj)*x(j)
30          CONTINUE
40        CONTINUE
50      CONTINUE
      ELSE
! Solve system using block structure information and calls to MA50C
!     for each diagonal block.
! System is block lower triangular.
        k2 = 0
        DO 100 jb = 1, numb
          nc = iblock(1,jb)
          k1 = k2 + 1
          k2 = k2 + nc
          IF (jb>1) THEN
! Substitution using off-diagonal block
            DO 70 j = k1, k2
              DO 60 jj = iptro(j), iptro(j+1) - 1
                i = irna(jj)
                rhs(j) = rhs(j) - aa(jj)*x(i)
60            CONTINUE
70          CONTINUE
          END IF
          IF (iblock(2,jb)<0) THEN
! Process triangular block
            DO 90 j = k1, k2
              DO 80 jj = iptrd(j) + 1, iptrd(j+1) - 1
                i = irna(jj)
                rhs(j) = rhs(j) - aa(jj)*x(i)
80            CONTINUE
              x(j) = rhs(j)/aa(iptrd(j))
90          CONTINUE
          ELSE
            j1 = 1
            IF (jb>1) j1 = iblock(3,jb)
            np = iblock(2,jb)
            CALL ma50cd(nc,nc,icntl5,iqb,np,trans,la+1-j1,a(j1), &
              irn(j1),iptrl(k1),iptru(k1),rhs(k1),x(k1),w)
          END IF
100     CONTINUE
      END IF

      deallocate(w,stat=stat)

      RETURN
    END subroutine ma48dd

    SUBROUTINE ma48id(cntl,icntl)
! Set default values for the control arrays.

      real(wp) cntl(10)
      INTEGER icntl(20)

! CNTL  is a real array of length 10.
!     CNTL(1)  If this is set to a value less than or equal to one, full
!       matrix processing will be used by MA50A/AD  when the density of
!       the reduced matrix reaches CNTL(1).
!     CNTL(2) determines the balance used by MA50A/AD and MA50B/BD
!       between pivoting for sparsity and for stability, values near
!       zero emphasizing sparsity and values near one emphasizing
!       stability.
!     CNTL(3) If this is set to a positive value, any entry whose
!       modulus is less than CNTL(3) will be dropped from the factors
!       calculated by MA50A/AD. The factorization will then require
!       less storage but will be inaccurate.
!     CNTL(4)  If this is set to a positive value, any entry whose
!       modulus is less than CNTL(4) will be regarded as zero from
!       the point of view of rank.
!     CNTL(5) is used in convergence test for termination of iterative
!       refinement.  Iteration stops if successive omegas do not
!       decrease by at least this factor.
!     CNTL(6) is set to the factor by which storage is increased if
!       arrays have to be reallocated.
!     CNTL(7:10) are not used.
! ICNTL is an integer array of length 20.
!     ICNTL(1)  must be set to the stream number for error messages.
!       A value less than 1 suppresses output.
!     ICNTL(2) must be set to the stream number for diagnostic output.
!       A value less than 1 suppresses output.
!     ICNTL(3) must be set to control the amount of output:
!       0 None.
!       1 Error messages only.
!       2 Error and warning messages.
!       3 As 2, plus scalar parameters and a few entries of array
!         parameters on entry and exit.
!       4 As 2, plus all parameters on entry and exit.
!     ICNTL(4)  If set to a positive value, the pivot search by MA50A/AD
!       is limited to ICNTL(4) columns. This may result in different
!       fill-in and execution time but could give faster execution.
!     ICNTL(5) The block size to be used for full-matrix processing.
!     ICNTL(6) is the minimum size for a block of the block triangular
!       form.
!     ICNTL(7) If not equal to 0, abort when structurally rank deficient
!       matrix found.
!     ICNTL(8) If set to a value other than zero and JOB = 1 or 3 on
!       entry to MA48A/AD, columns with IW flagged 0 are placed at
!       the end of their respective blocks in the first factorization
!       and the remaining columns are assumed unchanged in subsequent
!       factorizations. If set to a value other than zero and JOB = 2,
!       columns before the first 0 entry of IW are assumed unchanged in
!       subsequent factorizations.  On entry to MA48B/BD, the number
!       of columns that can change in each block is held in array KEEP.
!     ICNTL(9) is the limit on the number of iterative refinements
!        allowed by MA48C/CD
!     ICNTL(10) has default value 0. If set to 1, there is an immediate
!       return from MA48B/BD if LA is too small, without continuing the
!       decomposition to compute the size necessary.
!     ICNTL(11) has default value 0. If set to 1 on a JOB=2 call to
!       MA48B/BD and the entries in one of the blocks on the diagonal
!       are unsuitable for the pivot sequence chosen on the previous
!       call, the block is refactorized as on a JOB=1 call.
!     ICNTL(12:20) are not used.

      INTEGER i

      DO 10 i = 3, 10
        cntl(i) = 0.0D0
10    CONTINUE
      cntl(1) = 0.5D0
      cntl(2) = 0.1D0
      cntl(5) = 0.5D0
      cntl(6) = 2.0D0
      icntl(1) = 6
      icntl(2) = 6
      icntl(3) = 2
      icntl(4) = 3
      icntl(5) = 32
      icntl(6) = 1
      icntl(7) = 1
      icntl(8) = 0
      icntl(9) = 10
      DO 20 i = 10, 20
        icntl(i) = 0
20    CONTINUE

    END subroutine ma48id

! COPYRIGHT (c) 1977 AEA Technology
! Original date 8 Oct 1992
!######8/10/92 Toolpack tool decs employed.
!######8/10/92 D version created by name change only.
! 13/3/02 Cosmetic changes applied to reduce single/double differences
!
! 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE mc21ad(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,info)
!     .. Scalar Arguments ..
      INTEGER :: N,NUMNZ
      INTEGER(long) :: LICN,ICN(LICN)
!     ..
!     .. Array Arguments ..
      INTEGER IPERM(N),LENR(N)
      INTEGER(long) :: IP(N),info(20)
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
      CALL mc21bd(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,info)
      RETURN
!
      END SUBROUTINE
      SUBROUTINE mc21bd(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,info)
!   PR(I) IS THE PREVIOUS ROW TO I IN THE DEPTH FIRST SEARCH.
! IT IS USED AS A WORK ARRAY IN THE SORTING ALGORITHM.
!   ELEMENTS (IPERM(I),I) I=1, ... N  ARE NON-ZERO AT THE END OF THE
! ALGORITHM UNLESS N ASSIGNMENTS HAVE NOT BEEN MADE.  IN WHICH CASE
! (IPERM(I),I) WILL BE ZERO FOR N-NUMNZ ENTRIES.
!   CV(I) IS THE MOST RECENT ROW EXTENSION AT WHICH COLUMN I
! WAS VISITED.
!   ARP(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
! WHICH HAVE NOT BEEN SCANNED WHEN LOOKING FOR A CHEAP ASSIGNMENT.
!   OUT(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
! WHICH HAVE NOT BEEN SCANNED DURING ONE PASS THROUGH THE MAIN LOOP.

!   INITIALIZATION OF ARRAYS.
!     .. Scalar Arguments ..
      INTEGER :: N,NUMNZ
      INTEGER(long) :: LICN,ICN(LICN),info(20)
!     ..
!     .. Array Arguments ..
      INTEGER IPERM(N),LENR(N)
      INTEGER(long) :: IP(N)
!     ..
!     .. Local Scalars ..
      INTEGER I,IOUTK,J,J1,JORD,K,KK
      integer stat
      INTEGER(long) :: II,IN1,IN2

      integer,allocatable :: PR(:),ARP(:),CV(:),OUT(:)

      allocate(pr(n),arp(n),cv(n),out(n),stat=stat)
      if (stat .ne. 0) then
        info(1) = -10
        return
      endif
!     ..
!     .. Executable Statements ..
      DO 10 I = 1,N
        ARP(I) = LENR(I) - 1
        CV(I) = 0
        IPERM(I) = 0
   10 CONTINUE
      NUMNZ = 0


!   MAIN LOOP.
!   EACH PASS ROUND THIS LOOP EITHER RESULTS IN A NEW ASSIGNMENT
! OR GIVES A ROW WITH NO ASSIGNMENT.
      DO 100 JORD = 1,N
        J = JORD
        PR(J) = -1
        DO 70 K = 1,JORD
! LOOK FOR A CHEAP ASSIGNMENT
          IN1 = ARP(J)
          IF (IN1.LT.0) GO TO 30
          IN2 = IP(J) + LENR(J) - 1
          IN1 = IN2 - IN1
          DO 20 II = IN1,IN2
            I = ICN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
!   NO CHEAP ASSIGNMENT IN ROW.
          ARP(J) = -1
!   BEGIN LOOKING FOR ASSIGNMENT CHAIN STARTING WITH ROW J.
   30     CONTINUE
          OUT(J) = LENR(J) - 1
! INNER LOOP.  EXTENDS CHAIN BY ONE OR BACKTRACKS.
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENR(J) - 1
            IN1 = IN2 - IN1
! FORWARD SCAN.
            DO 40 II = IN1,IN2
              I = ICN(II)
              IF (CV(I).EQ.JORD) GO TO 40
!   COLUMN I HAS NOT YET BEEN ACCESSED DURING THIS PASS.
              J1 = J
              J = IPERM(I)
              CV(I) = JORD
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70

   40       CONTINUE

!   BACKTRACKING STEP.
   50       CONTINUE
            J = PR(J)
            IF (J.EQ.-1) GO TO 100
   60     CONTINUE

   70   CONTINUE

!   NEW ASSIGNMENT IS MADE.
   80   CONTINUE
        IPERM(I) = J
        ARP(J) = IN2 - II - 1
        NUMNZ = NUMNZ + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 100
          II = IP(J) + LENR(J) - OUT(J) - 2
          I = ICN(II)
          IPERM(I) = J
   90   CONTINUE

  100 CONTINUE

!   IF MATRIX IS STRUCTURALLY SINGULAR, WE NOW COMPLETE THE
! PERMUTATION IPERM.
      IF (NUMNZ.EQ.N) GO TO 200
      DO 110 I = 1,N
        ARP(I) = 0
  110 CONTINUE
      K = 0
      DO 130 I = 1,N
        IF (IPERM(I).NE.0) GO TO 120
        K = K + 1
        OUT(K) = I
        GO TO 130

  120   CONTINUE
        J = IPERM(I)
        ARP(J) = I
  130 CONTINUE
      K = 0
      DO 140 I = 1,N
        IF (ARP(I).NE.0) GO TO 140
        K = K + 1
        IOUTK = OUT(K)
        IPERM(IOUTK) = I
  140 CONTINUE
  
  200 deallocate(cv,out,pr,arp,stat=stat)

      RETURN

      END  subroutine
! COPYRIGHT (c) 1976 AEA Technology
! Original date 21 Jan 1993
!       Toolpack tool decs employed.
! Double version of MC13D (name change only)
! 10 August 2001 DOs terminated with CONTINUE
! 13/3/02 Cosmetic changes applied to reduce single/double differences
!
! 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE mc13dd(N,ICN,LICN,IP,LENR,IOR,IB,NUM,info)
!     .. Scalar Arguments ..
      INTEGER N,NUM
      INTEGER(long) :: LICN,ICN(LICN),info(20)

!     .. Array Arguments ..
      INTEGER IB(N),LENR(N)
! IOR is only long because of calling program (keep)
      INTEGER(long) :: IP(N),IOR(N)

!     .. Executable Statements ..
      CALL mc13ed(N,ICN,LICN,IP,LENR,IOR,IB,NUM,info)
      RETURN

      END SUBROUTINE
      SUBROUTINE mc13ed(N,ICN,LICN,IP,LENR,ARP,IB,NUM,info)
!
! ARP(I) IS ONE LESS THAN THE NUMBER OF UNSEARCHED EDGES LEAVING
!     NODE I.  AT THE END OF THE ALGORITHM IT IS SET TO A
!     PERMUTATION WHICH PUTS THE MATRIX IN BLOCK LOWER
!     TRIANGULAR FORM.
! IB(I) IS THE POSITION IN THE ORDERING OF THE START OF THE ITH
!     BLOCK.  IB(N+1-I) HOLDS THE NODE NUMBER OF THE ITH NODE
!     ON THE STACK.
! LOWL(I) IS THE SMALLEST STACK POSITION OF ANY NODE TO WHICH A PATH
!     FROM NODE I HAS BEEN FOUND.  IT IS SET TO N+1 WHEN NODE I
!     IS REMOVED FROM THE STACK.
! NUMB(I) IS THE POSITION OF NODE I IN THE STACK IF IT IS ON
!     IT, IS THE PERMUTED ORDER OF NODE I FOR THOSE NODES
!     WHOSE FINAL POSITION HAS BEEN FOUND AND IS OTHERWISE ZERO.
! PREV(I) IS THE NODE AT THE END OF THE PATH WHEN NODE I WAS
!     PLACED ON THE STACK.


!   ICNT IS THE NUMBER OF NODES WHOSE POSITIONS IN FINAL ORDERING HAVE
!     BEEN FOUND.
!     .. Scalar Arguments ..
      INTEGER N,NUM
      INTEGER(long) :: LICN,ICN(LICN),info(20)
!     ..
!     .. Array Arguments ..
      INTEGER IB(N),LENR(N)
      INTEGER(long) :: IP(N),ARP(N)
!     ..
!     .. Local Scalars ..
      INTEGER DUMMY,I,I1,I2,ICNT,II,ISN,IST,IST1,IV,IW,J,K,LCNT,NNM1,STP

      integer,allocatable :: lowl(:),numb(:),prev(:)
      integer :: stat
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MIN
!     ..
!     .. Executable Statements ..

      allocate(lowl(n),numb(n),prev(n),stat=stat)
      if (stat .ne. 0) then
        info(1) = -10
        return
      endif
      ICNT = 0
! NUM IS THE NUMBER OF BLOCKS THAT HAVE BEEN FOUND.
      NUM = 0
      NNM1 = N + N - 1

! INITIALIZATION OF ARRAYS.
      DO 20 J = 1,N
        NUMB(J) = 0
        ARP(J) = LENR(J) - 1
   20 CONTINUE


      DO 120 ISN = 1,N
! LOOK FOR A STARTING NODE
        IF (NUMB(ISN).NE.0) GO TO 120
        IV = ISN
! IST IS THE NUMBER OF NODES ON THE STACK ... IT IS THE STACK POINTER.
        IST = 1
! PUT NODE IV AT BEGINNING OF STACK.
        LOWL(IV) = 1
        NUMB(IV) = 1
        IB(N) = IV

! THE BODY OF THIS LOOP PUTS A NEW NODE ON THE STACK OR BACKTRACKS.
        DO 110 DUMMY = 1,NNM1
          I1 = ARP(IV)
! HAVE ALL EDGES LEAVING NODE IV BEEN SEARCHED.
          IF (I1.LT.0) GO TO 60
          I2 = IP(IV) + LENR(IV) - 1
          I1 = I2 - I1
!
! LOOK AT EDGES LEAVING NODE IV UNTIL ONE ENTERS A NEW NODE OR
!     ALL EDGES ARE EXHAUSTED.
          DO 50 II = I1,I2
            IW = ICN(II)
! HAS NODE IW BEEN ON STACK ALREADY.
            IF (NUMB(IW).EQ.0) GO TO 100
! UPDATE VALUE OF LOWL(IV) IF NECESSARY.
            LOWL(IV) = MIN(LOWL(IV),LOWL(IW))
   50     CONTINUE
!
! THERE ARE NO MORE EDGES LEAVING NODE IV.
          ARP(IV) = -1
! IS NODE IV THE ROOT OF A BLOCK.
   60     IF (LOWL(IV).LT.NUMB(IV)) GO TO 90

! ORDER NODES IN A BLOCK.
          NUM = NUM + 1
          IST1 = N + 1 - IST
          LCNT = ICNT + 1
! PEEL BLOCK OFF THE TOP OF THE STACK STARTING AT THE TOP AND
!     WORKING DOWN TO THE ROOT OF THE BLOCK.
          DO 70 STP = IST1,N
            IW = IB(STP)
            LOWL(IW) = N + 1
            ICNT = ICNT + 1
            NUMB(IW) = ICNT
            IF (IW.EQ.IV) GO TO 80
   70     CONTINUE
   80     IST = N - STP
          IB(NUM) = LCNT
! ARE THERE ANY NODES LEFT ON THE STACK.
          IF (IST.NE.0) GO TO 90
! HAVE ALL THE NODES BEEN ORDERED.
          IF (ICNT.LT.N) GO TO 120
          GO TO 130

! BACKTRACK TO PREVIOUS NODE ON PATH.
   90     IW = IV
          IV = PREV(IV)
! UPDATE VALUE OF LOWL(IV) IF NECESSARY.
          LOWL(IV) = MIN(LOWL(IV),LOWL(IW))
          GO TO 110

! PUT NEW NODE ON THE STACK.
  100     ARP(IV) = I2 - II - 1
          PREV(IW) = IV
          IV = IW
          IST = IST + 1
          LOWL(IV) = IST
          NUMB(IV) = IST
          K = N + 1 - IST
          IB(K) = IV
  110   CONTINUE

  120 CONTINUE


! PUT PERMUTATION IN THE REQUIRED FORM.
  130 DO 140 I = 1,N
        II = NUMB(I)
        ARP(II) = I
  140 CONTINUE


      deallocate(lowl,numb,prev,stat=stat)

      RETURN

      END subroutine
end module hsl_ma48_ma48_internal_double

module hsl_ma48_ma51_internal_double
   implicit none

! This module is based on ma51 1.0.0 (12 July 2004)

   private
   public :: ma51ad, ma51cd
   public :: ma51bd, ma51dd ! public for unit tests

   integer, parameter :: wp = kind(0.0d0)
   integer, parameter :: long = selected_int_kind(18) ! Long integer

contains

    SUBROUTINE ma51ad(m,n,la,irn,keep,rank,rows,cols,w)

!     This is for use in conjunction with the MA48 routines
!     for solving full sets of linear equations. It must follow
!     a call of MA48B/BD, and the values of the arguments
!     M,N,LA,IRN,KEEP must be unchanged.
!     It identifies which equations are ignored when solving Ax = b
!     and which solution components are always set to zero.
!     There are such equations and/or components in the singular or
!     rectangular case.


!     .. Arguments ..
      integer(long) :: la, keep(*), irn(la)
      integer :: m, n, rank, rows(m), cols(n), w(max(m,n))

! M must be set by the user to the number of rows in the matrix.
!      It is not altered by the subroutine. Restriction: M > 0.
! N must be set by the user to the number of columns in the matrix.
!      It is not altered by the subroutine. Restriction: N > 0.
! LA must be set by the user to the size of array IRN.
!      It is not altered by the subroutine.
! IRN  must be left unchanged from the last call to MA48B/BD.
!      It holds the row indices of the factorization in A and
!      the row indices of the original matrix in permuted form.
!      It is not altered by the subroutine.
! KEEP must be as on return from MA48B/BD.
!      It is not altered by the subroutine.
! RANK is an integer variable that need not be set on input.
!     On exit, it holds the rank computed by MA50.
! ROWS is an integer array that need not be set on input.
!      On exit, it holds a permutation. The indices of the rows that
!      are taken into account when solving Ax = b are ROWS(i),
!      i <= RANK.
! COLS is an integer array that need not be set on input.
!      On exit, it holds a permutation. The indices of the columns that
!      are taken into account when solving Ax = b are COLS(j),
!      j <= RANK.
! W is an integer workarray of length max(M,N).

!     .. Local variables ..
      integer :: i, iqb(1), j, jb, j1, k1, k2, &
        nb, nc, np, nr, rankb
      integer(long) :: iptrd, iptrl, iptro, iptru, &
        kblock, mblock, nblock, ne

      iqb(1) = 0

! Partition KEEP
      iptrl = m + n
      iptru = iptrl + n
      iptrd = iptru + n
      iptro = iptrd + n + 1
      nblock = iptro + n - 1
      mblock = nblock + 1
      kblock = mblock + 1
      nb = keep(kblock+3)
      ne = keep(iptro+n+1) - 1


      k2 = 0
      rank = 0
      DO 20 jb = 1, nb
        nc = keep(nblock+3*jb)
        nr = nc
        IF (nb==1) nr = m
        k1 = k2 + 1
        k2 = k2 + nc
        np = keep(mblock+3*jb)
        IF (np<0) THEN
! Triangular block
          DO 10 j = k1, k2
            rows(j) = 1
            cols(j) = 1
10        CONTINUE
          rank = rank + nc
        ELSE
          j1 = 1
          IF (jb>1) j1 = keep(kblock+3*jb)
          CALL ma51zd(nr,nc,iqb,np,la+1-ne-j1,irn(ne+j1),keep(iptrl+k1), &
            keep(iptru+k1),rankb,rows(k1),cols(k1),w)
          rank = rank + rankb
        END IF
20    CONTINUE

! Construct final row permutation
      DO 30 i = 1, m
        w(i) = rows(keep(i))
30    CONTINUE
      k1 = 1
      k2 = m
      DO 40 i = 1, m
        IF (w(i)==1) THEN
          rows(k1) = i
          k1 = k1 + 1
        ELSE
          rows(k2) = i
          k2 = k2 - 1
        END IF
40    CONTINUE

! Construct final column permutation
      DO 50 i = 1, n
        w(keep(m+i)) = cols(i)
50    CONTINUE
      k1 = 1
      k2 = n
      DO 60 i = 1, n
        IF (w(i)==1) THEN
          cols(k1) = i
          k1 = k1 + 1
        ELSE
          cols(k2) = i
          k2 = k2 - 1
        END IF
60    CONTINUE

    END subroutine ma51ad

    SUBROUTINE ma51bd(m,n,iq,np,lfact,irnf,iptrl,iptru,rank,rows,cols,w)

!  Purpose
!  =======
!     This is for use in conjunction with the MA50 routines
!     for solving full sets of linear equations. It must follow
!     a call of MA50B/BD, and the values of the arguments
!     M,N,IQ,NP,LFACT,IRNF,IPTRL,IPTRU must be unchanged.
!     It identifies which equations are ignored when solving Ax = b
!     and which solution components are always set to zero.
!     There are such equations and/or components in the singular or
!     rectangular case.

      integer  :: m, n, np, iq(*)
      integer(long) :: lfact, iptrl(n), iptru(n), irnf(lfact)
      integer :: rank, rows(m), cols(n), w(max(m,n))

! M  is an integer variable set to the number of rows.
!     It is not altered by the subroutine.
! N  is an integer variable set to the number of columns.
!     It is not altered by the subroutine.
! IQ is an integer array holding the permutation Q.
!     It is not altered by the subroutine.
! NP is an integer variable holding the number of rows and columns in
!     packed storage. It is not altered by the subroutine.
! LFACT is an integer variable set to the size of FACT and IRNF.
!     It is not altered by the subroutine.
! IRNF is an integer array holding the row numbers of the packed part
!     of L/U, and the row numbers of the full part of L/U.
!     It is not altered by the subroutine.
! IPTRL is an integer array. For J = 1,..., NP, IPTRL(J) holds the
!     position in FACT and IRNF of the end of column J of L.
!     For J = NP+1,..., N, IPTRL(J) is equal to IPTRU(J).
!     It is not altered by the subroutine.
! IPTRU is an integer array. For J = 1,..., N, IPTRU(J) holds the
!     position in FACT and IRNF of the end of the packed part of
!     column J of U. It is not altered by the subroutine.
! RANK is an integer variable that need not be set on input.
!     On exit, it holds the rank computed by MA50.
! ROWS is an integer array that need not be set on input.
!      On exit, it holds a permutation. The indices of the rows that
!      are taken into account when solving Ax = b are ROWS(i),
!      i <= RANK.
! COLS is an integer array that need not be set on input.
!      On exit, it holds a permutation. The indices of the columns that
!      are taken into account when solving Ax = b are COLS(j),
!      j <= RANK.
! W is an integer workarray of length max(M,N).


! Local variables
!  =========


      INTEGER i, k1, k2
! I Temporary variable holding row or column number.
! K1 Position of next row or column of the significant part.
! K2 Position of next row or column of the insignificant part.

      CALL ma51zd(m,n,iq,np,lfact,irnf,iptrl,iptru,rank,rows,cols,w)

! Construct final permutations
      DO 10 i = 1, m
        w(i) = rows(i)
10    CONTINUE
      k1 = 1
      k2 = m
      DO 20 i = 1, m
        IF (w(i)==1) THEN
          rows(k1) = i
          k1 = k1 + 1
        ELSE
          rows(k2) = i
          k2 = k2 - 1
        END IF
20    CONTINUE
      DO 30 i = 1, n
        w(i) = cols(i)
30    CONTINUE
      k1 = 1
      k2 = n
      DO 40 i = 1, n
        IF (w(i)==1) THEN
          cols(k1) = i
          k1 = k1 + 1
        ELSE
          cols(k2) = i
          k2 = k2 - 1
        END IF
40    CONTINUE
    END subroutine ma51bd

    SUBROUTINE ma51cd(m,n,la,a,irn,keep,sgndet,logdet,w)
!     This is for use in conjunction with the MA48 routines
!     for solving sparse sets of linear equations. It must follow
!     a call of MA48B/BD, and the values of the arguments
!     M,N,LA,A,IRN,KEEP must be unchanged.
!     It computes the determinant.

      INTEGER m, n
      integer(long) la
      real(wp) a(la)
      INTEGER sgndet
      INTEGER(long) irn(la)
      integer(long) keep(*)
      real(wp) logdet
      integer w(n)
      real(wp) zero
      PARAMETER (zero=0.0_wp)

! M must be set by the user to the number of rows in the matrix.
!      It is not altered by the subroutine.
! N must be set by the user to the number of columns in the matrix.
!      It is not altered by the subroutine.
! LA must be set by the user to the size of array IRN.
!      It is not altered by the subroutine.
! A must be left unchanged from the last call to MA48B/BD.
!      It holds the factorized matrix.
!      It is not altered by the subroutine.
! IRN  must be left unchanged from the last call to MA48B/BD.
!      It holds the row indices of the factorization in A and
!      the row indices of the original matrix in permuted form.
!      It is not altered by the subroutine.
! KEEP must be as on return from MA48B/BD.
!      It is not altered by the subroutine.
! SGNDET is an integer output variable that returns the
!      sign of the determinant or zero if the determinant is zero.
! LOGDET is a double precision output variable that returns the
!      LOG of the abslute value of the determinant or zero if the
!      determinant is zero.
! W is an integer work array.

!     .. Local variables ..
      INTEGER i, iqb(1), j, jb, j1, k, k1, k2, &
        l, nb, nc, np, nr, sgndt
      integer(long) :: iptrd, iptrl, iptro, iptru, &
        kblock, mblock, nblock, ne
      real(wp) piv, logdt

      iqb(1) = 0

! Partition KEEP
      iptrl = m + n
      iptru = iptrl + n
      iptrd = iptru + n
      iptro = iptrd + n + 1
      nblock = iptro + n - 1
      mblock = nblock + 1
      kblock = mblock + 1
      nb = keep(kblock+3)
      ne = keep(iptro+n+1) - 1

      logdet = zero
      sgndet = 1

! Check for rectangular case
      IF (m/=n) THEN
        sgndet = 0
        RETURN
      END IF

      DO 10 i = 1, n
        w(i) = 1
10    CONTINUE

! Compute sign of row permutation, leaving W=0
      DO 40 i = 1, n
        k = keep(i)
        IF (w(k)==1) THEN
          DO 30 j = 1, n
            l = keep(k)
            w(k) = 0
            IF (k==i) GO TO 40
            sgndet = -sgndet
            k = l
30        CONTINUE
        END IF
40    CONTINUE

! Compute sign of column permutation
      DO 60 i = 1, n
        k = keep(n+i)
        IF (w(k)==0) THEN
          DO 50 j = 1, n
            l = keep(n+k)
            w(k) = 1
            IF (k==i) GO TO 60
            sgndet = -sgndet
            k = l
50        CONTINUE
        END IF
60    CONTINUE

      k2 = 0
      DO 80 jb = 1, nb
        nc = keep(nblock+3*jb)
        nr = nc
        IF (nb==1) nr = m
        k1 = k2 + 1
        k2 = k2 + nc
        np = keep(mblock+3*jb)
        IF (np<0) THEN
! Triangular block
          DO 70 j = k1, k2
            piv = a(keep(iptrd+j))
            IF (piv<zero) THEN
              sgndet = -sgndet
              logdet = logdet + log(-piv)
            ELSE
              logdet = logdet + log(piv)
            END IF
70        CONTINUE
        ELSE
          j1 = 1
          IF (jb>1) j1 = keep(kblock+3*jb)
          CALL ma51dd(nr,nc,iqb,np,la+1-j1-ne,a(ne+j1),irn(ne+j1), &
            keep(iptrl+k1),keep(iptru+k1),sgndt,logdt,w)
          sgndet = sgndet*sgndt
          logdet = logdet + logdt
          IF (sgndet==0) THEN
            logdet = 0
            RETURN
          END IF
        END IF
80    CONTINUE

    END subroutine ma51cd

    SUBROUTINE ma51dd(m,n,iq,np,lfact,fact,irnf,iptrl,iptru,sgndet,logdet,w)
!  Computes the determinant of a factorized n-by-n matrix.
      INTEGER k, l, m, n, np, iq(*)
      integer(long) ::  lfact, irnf(lfact)
      real(wp) fact(lfact)
      INTEGER sgndet, w(n)
      integer(long) ::  iptrl(n), iptru(n)
      real(wp) logdet
! M is an integer variable that must be set to the number of rows.
!      It is not altered by the subroutine.
! N is an integer variable that must be set to the number of columns.
!      It is not altered by the subroutine.
! IQ is an integer array of length N holding the column permutation Q
!     or an array of length 1 with the value 0 if there is no column
!     permutation. It is not altered by the subroutine.
! NP is an integer variable that must be unchanged since calling
!     MA50B/BD. It holds the number of rows and columns in packed
!     storage. It is not altered by the subroutine.
! LFACT is an integer variable set to the size of FACT and IRNF.
!     It is not altered by the subroutine.
! FACT is an array that must be unchanged since calling MA50B/BD. It
!     holds the packed part of L/U by columns, and the full part of L/U
!     by columns. U has unit diagonal entries, which are not stored, and
!     the signs of the off-diagonal entries are inverted.  In the packed
!     part, the entries of U precede the entries of L; also the diagonal
!     entries of L head each column of L and are reciprocated.
!     FACT is not altered by the subroutine.
! IRNF is an integer array that must be unchanged since calling
!     MA50B/BD. It holds the row numbers of the packed part of L/U, and
!     the row numbers of the full part of L/U.
!     It is not altered by the subroutine.
! IPTRL is an integer array that must be unchanged since calling
!     MA50B/BD. For J = 1,..., NP, IPTRL(J) holds the position in
!     FACT and IRNF of the end of column J of L.
!     For J = NP+1,..., N, IPTRL(J) is equal to IPTRU(J).
!     It is not altered by the subroutine.
! IPTRU is an integer array that must be unchanged since calling
!     MA50B/BD. For J = 1,..., N, IPTRU(J) holds the position in
!     FACT and IRNF of the end of the packed part of column J of U.
!     It is not altered by the subroutine.
! SGNDET is an integer output variable that returns the
!      sign of the determinant or zero if the determinant is zero.
! LOGDET is a double precision output variable that returns the
!      LOG of the abslute value of the determinant or zero if the
!      determinant is zero.
! W is an integer work array.

      real(wp) a
      INTEGER i, j, mf, nf
      integer(long) :: ia1, if1
      real(wp) zero
      PARAMETER (zero=0.0_wp)
! A Pivot value.
! I Temporary variable.
! IA1 Position of the start of the current row or column.
! IF1 Position of the start of the full part of U.
! J Temporary variable.
! MF Number of rows held in full format.
! NF Number of columns held in full format.

      if1 = iptrl(n) + 1
      mf = irnf(2)
      nf = n - np
      logdet = 0
      sgndet = 1

! Find the determinant of the full part
      IF (mf>0) THEN
! Check for rank-deficient case
        IF (m/=n .OR. mf/=nf .OR. irnf(if1+mf+nf-1)<0) THEN
          sgndet = 0
          RETURN
        END IF
        CALL ma51xd(mf,fact(if1),mf,irnf(if1+mf),sgndet,logdet)
      END IF

! Find row permutation and compute LOGDET
      DO 10 i = 1, np
        ia1 = iptru(i) + 1
        a = fact(ia1)
        IF (a<zero) THEN
          sgndet = -sgndet
          logdet = logdet - log(-a)
        ELSE
          logdet = logdet - log(a)
        END IF
        w(i) = irnf(ia1)
10    CONTINUE
      DO 20 i = 1, mf
        w(np+i) = irnf(if1+i-1)
20    CONTINUE

! Compute sign of row permutation, leaving W=0
      DO 40 i = 1, n
        k = w(i)
        IF (k>0) THEN
          DO 30 j = 1, n
            l = w(k)
            w(k) = 0
            IF (k==i) GO TO 40
            sgndet = -sgndet
            k = l
30        CONTINUE
        END IF
40    CONTINUE
      IF (iq(1)<=0) RETURN

! Compute sign of column permutation
      DO 60 i = 1, n
        k = iq(i)
        IF (w(k)==0) THEN
          DO 50 j = 1, n
            l = iq(k)
            w(k) = 1
            IF (k==i) GO TO 60
            sgndet = -sgndet
            k = l
50        CONTINUE
        END IF
60    CONTINUE

    END subroutine ma51dd



    SUBROUTINE ma51xd(n,a,lda,ipiv,sgndet,logdet)
      IMPLICIT NONE
      INTEGER n
      integer :: lda
      real(wp) a(lda,n)
      INTEGER sgndet
      integer(long) ipiv(n)
      real(wp) logdet
!  Computes the determinant of a factorized n-by-n matrix A.

!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix of order n, L is lower triangular
!  of order n with unit diagonal elements, U is upper triangular of
!  order n, and Q is a permutation matrix of order n.

!  N       (input) INTEGER
!          Order of the matrix A.  N >= 1.

!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          Holds the factors L and U of the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.

!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= N.

!  IPIV    (input) INTEGER array, dimension (N)
!          The permutations; for 1 <= i <= N, row i of the
!          matrix was interchanged with row IPIV(i).

!  SGNDET  (output) INTEGER
!          The sign of the determinant.

!  LOGDET  (output) DOUBLE PRECISION
!          LOG of the determinant.

      INTEGER j
      real(wp) zero
      PARAMETER (zero=0.0_wp)
      INTRINSIC log

      sgndet = 1
      logdet = zero
      DO 30 j = 1, n
        IF (ipiv(j)/=j) sgndet = -sgndet
        IF (a(j,j)<zero) THEN
          sgndet = -sgndet
          logdet = logdet + log(-a(j,j))
        ELSE
          logdet = logdet + log(a(j,j))
        END IF
30    CONTINUE

    END subroutine ma51xd



    SUBROUTINE ma51yd(m,n,ipiv,rank,rows,cols)

!  Purpose
!  =======
!     This is for use in conjunction with the MA50 routines
!     for solving full sets of linear equations. It must follow
!     a call of MA50E/ED, MA50F/FD, or MA50G/GD.
!     It identifies which equations are ignored when solving Ax = b
!     and which solution components are always set to zero.
!     There are such equations and/or components in the singular or
!     rectangular case.

      INTEGER m, n, rank, rows(m), cols(n)
      integer (long) ipiv(n)

!  Arguments
!  =========

!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 1.

!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 1.

!  IPIV    (input) INTEGER array, dimension (N)
!          The permutations; for 1 <= i <= RANK, row i of the
!          matrix was interchanged with row IPIV(i); for
!          RANK + 1 <= j <= N, column j of the
!          matrix was interchanged with column -IPIV(j).

!  RANK    (output) INTEGER
!          On exit, the rank computed by MA50.

!  ROWS    (output) INTEGER array, dimension (M)
!          On exit, ROWS(i) = 0 if row i is ignored and
!          ROWS(i) = 1 otherwise.

!  COLS    (output) INTEGER array, dimension (N)
!          On exit, COLS(j) = 0 if component j is set to zero
!          always and COLS(j) = 1 otherwise.


! Local variables
!  =========

      INTEGER i, j, k
! I    Temporary variable.
! J    Temporary variable.
! K    Temporary variable.

! Find rank
      DO 10 k = 1, n
        IF (ipiv(k)<0) GO TO 20
        rows(k) = 1
        cols(k) = 1
10    CONTINUE
20    rank = k - 1

! Set singular parts to zero
      DO 30 k = rank + 1, m
        rows(k) = 0
30    CONTINUE
      DO 40 k = rank + 1, n
        cols(k) = 0
40    CONTINUE

! Apply interchanges.
      DO 50 i = rank + 1, n
        k = -ipiv(i)
        j = cols(i)
        cols(i) = cols(k)
        cols(k) = j
50    CONTINUE
      DO 60 i = rank, 1, -1
        k = ipiv(i)
        j = rows(i)
        rows(i) = rows(k)
        rows(k) = j
60    CONTINUE

    END subroutine ma51yd


    SUBROUTINE ma51zd(m,n,iq,np,lfact,irnf,iptrl,iptru,rank,rows,cols,w)

!  Purpose
!  =======
!     This is for use in conjunction with the MA50 routines
!     for solving full sets of linear equations. It must follow
!     a call of MA50B/BD, and the values of the arguments
!     M,N,IQ,NP,LFACT,IRNF,IPTRL,IPTRU must be unchanged.
!     It identifies which equations are ignored when solving Ax = b
!     and which solution components are always set to zero.
!     There are such equations and/or components in the singular or
!     rectangular case.

      integer  :: m, n, np, iq(*)
      integer(long) :: lfact, iptrl(n), iptru(n), irnf(lfact)
      integer :: rank, rows(m), cols(n), w(*)

! M  is an integer variable set to the number of rows.
!     It is not altered by the subroutine.
! N  is an integer variable set to the number of columns.
!     It is not altered by the subroutine.
! IQ is an integer array holding the permutation Q.
!     It is not altered by the subroutine.
! NP is an integer variable holding the number of rows and columns in
!     packed storage. It is not altered by the subroutine.
! LFACT is an integer variable set to the size of FACT and IRNF.
!     It is not altered by the subroutine.
! IRNF is an integer array holding the row numbers of the packed part
!     of L/U, and the row numbers of the full part of L/U.
!     It is not altered by the subroutine.
! IPTRL is an integer array. For J = 1,..., NP, IPTRL(J) holds the
!     position in FACT and IRNF of the end of column J of L.
!     For J = NP+1,..., N, IPTRL(J) is equal to IPTRU(J).
!     It is not altered by the subroutine.
! IPTRU is an integer array. For J = 1,..., N, IPTRU(J) holds the
!     position in FACT and IRNF of the end of the packed part of
!     column J of U. It is not altered by the subroutine.
! RANK is an integer variable that need not be set on input.
!     On exit, it holds the rank computed by MA50.
! ROWS is an integer array that need not be set on input.
!      On exit, ROWS(i) = 0 if row i is ignored and
!          ROWS(i) = 1 otherwise.
! COLS is an integer array that need not be set on input.
!      On exit, COLS(j) = 0 if component j is set to zero
!          always and COLS(j) = 1 otherwise.
! W is an integer workarray of length max(M,N).


! Local variables
!  =========


      INTEGER i, j, mf, nf
      integer(long) :: ia1, if1
! I Temporary variable holding row number.
! IA1 Position of the start of the current row or column.
! IF1 Position of the start of the full part of U.
! J Temporary variable holding column number.
! MF Number of rows held in full format.
! NF Number of columns held in full format.

      if1 = iptrl(n) + 1
      mf = irnf(2)
      nf = n - np

      IF (mf>0 .AND. nf>0) THEN
        CALL ma51yd(mf,nf,irnf(if1+mf),rank,w,cols(np+1))
      ELSE
        rank = 0
        DO 10 i = 1, mf
          w(i) = 0
10      CONTINUE
        DO 20 i = 1, nf
          cols(np+i) = 0
20      CONTINUE
      END IF

      DO 30 i = 1, m
        rows(i) = 1
30    CONTINUE


      DO 40 i = mf, 1, -1
        j = irnf(if1+i-1)
        rows(j) = w(i)
40    CONTINUE

      rank = rank + m - mf

      DO 220 j = np, 1, -1
        ia1 = iptru(j)
        IF (ia1>=iptrl(j)) THEN
          cols(j) = 0
        ELSE
          cols(j) = 1
        END IF
220   CONTINUE
      IF (iq(1)>0) THEN
! Permute COLS
        DO 230 i = 1, n
          w(i) = cols(i)
230     CONTINUE
        DO 240 i = 1, n
          cols(iq(i)) = w(i)
240     CONTINUE
      END IF

    END subroutine ma51zd

end module hsl_ma48_ma51_internal_double


module hsl_ma48_double
   use hsl_zd11_double
   use hsl_ma48_ma48_internal_double
   use hsl_ma48_ma51_internal_double
   implicit none

   private
   public :: ma48_factors,ma48_control,ma48_ainfo,ma48_finfo,ma48_sinfo, &
             ma48_initialize,ma48_analyse,ma48_factorize,ma48_solve, &
             ma48_finalize, ma48_get_perm,ma48_special_rows_and_cols, &
             ma48_determinant

   integer, parameter :: wp = kind(0.0d0)
   integer, parameter :: long = selected_int_kind(18) ! Long integer

   interface ma48_initialize
      module procedure ma48_initialize_double
   end interface

   interface ma48_analyse
      module procedure ma48_analyse_double
   end interface

   interface ma48_factorize
      module procedure ma48_factorize_double
   end interface

   interface ma48_solve
      module procedure ma48_solve_double
   end interface

   interface  ma48_finalize
      module procedure  ma48_finalize_double
   end interface

   interface ma48_get_perm
      module procedure ma48_get_perm_double
   end interface

   interface ma48_special_rows_and_cols
      module procedure ma48_special_rows_and_cols_double
   end interface

   interface ma48_determinant
      module procedure ma48_determinant_double
   end interface

   type ma48_factors
      integer(long), allocatable :: keep(:)
!! Must be long because used for mapping vectors
      integer(long), allocatable :: irn(:)
      integer(long), allocatable :: jcn(:)
      real(wp), allocatable :: val(:)
      integer :: m       ! Number of rows in matrix
      integer :: n       ! Number of columns in matrix
!! For the moment lirnreq = lareq because of BTF structure
      integer(long) :: lareq   ! Size of real array for further factorization
      integer(long) :: lirnreq ! Size of integer array for further 
                               ! factorization
      integer :: partial ! Number of columns kept to end for partial
                         ! factorization
      integer(long) :: ndrop   ! Number of entries dropped from data structure
      integer :: first   ! Flag to indicate whether it is first call to
                         ! factorize after analyse (set to 1 in analyse).
   end type ma48_factors

   type ma48_control
      real(wp) :: multiplier ! Factor by which arrays sizes are to be
                        ! increased if they are too small
      real(wp) :: u     ! Pivot threshold
      real(wp) :: switch ! Density for switch to full code
      real(wp) :: drop   ! Drop tolerance
      real(wp) :: tolerance ! anything less than this is considered zero
      real(wp) :: cgce  ! Ratio for required reduction using IR
      real(wp) :: reduce ! Only kept for compatibility with earlier
!                         version of hsl_ma48
      integer :: la     ! Only kept for compatibility with earlier
!                         version of hsl_ma48
      integer :: maxla  ! Only kept for compatibility with earlier
!                         version of hsl_ma48
      integer :: lp     ! Unit for error messages
      integer :: wp     ! Unit for warning messages
      integer :: mp     ! Unit for monitor output
      integer :: ldiag  ! Controls level of diagnostic output
      integer :: btf    ! Minimum block size for BTF ... >=N to avoid
      logical :: struct ! Control to abort if structurally singular
      integer :: maxit ! Maximum number of iterations
      integer :: factor_blocking ! Level 3 blocking in factorize
      integer :: solve_blas ! Switch for using Level 1 or 2 BLAS in solve.
      integer :: pivoting  ! Controls pivoting:
!                 Number of columns searched.  Zero for Markowitz
      logical :: diagonal_pivoting  ! Set to 0 for diagonal pivoting
      integer :: fill_in ! Initially fill_in * ne space allocated for factors
      logical :: switch_mode ! Whether to switch to slow when fast mode
!               given unsuitable pivot sequence.
   end type ma48_control

   type ma48_ainfo
      real(wp) :: ops = 0.0   ! Number of operations in elimination
      integer :: flag = 0  ! Flags success or failure case
      integer :: more = 0   ! More information on failure
      integer(long) :: lena_analyse  = 0! Size for analysis (main arrays)
      integer(long) :: lenj_analyse  = 0! Size for analysis (integer aux array)
      integer(long) :: len_analyse  = 0! Only kept for compatibility
!                      with earlier version of hsl_ma48. Set to maximum
!                      of lena_analyse and lenj_analyse
!! For the moment leni_factorize = lena_factorize because of BTF structure
      integer(long) :: lena_factorize  = 0 ! Size for factorize (real array)
      integer(long) :: leni_factorize  = 0 ! Size for factorize (integer array)
      integer(long) :: len_factorize  = 0! Only kept for compatibility
!                      with earlier version of hsl_ma48. Set to maximum
!                      of lena_factorize and leni_factorize
      integer :: ncmpa   = 0 ! Number of compresses in analyse
      integer :: rank   = 0  ! Estimated rank
      integer(long) :: drop   = 0  ! Number of entries dropped
      integer :: struc_rank  = 0! Structural rank of matrix
      integer(long) :: oor   = 0   ! Number of indices out-of-range
      integer(long) :: dup   = 0   ! Number of duplicates
      integer :: stat   = 0  ! STAT value after allocate failure
      integer :: lblock  = 0 ! Size largest non-triangular block
      integer :: sblock  = 0 ! Sum of orders of non-triangular blocks
      integer(long) :: tblock  = 0 ! Total entries in all non-triangular blocks
   end type ma48_ainfo

   type ma48_finfo
      real(wp) :: ops  = 0.0  ! Number of operations in elimination
      integer :: flag   = 0 ! Flags success or failure case
      integer :: more   = 0  ! More information on failure
      integer(long) :: size_factor   = 0! Number of words to hold factors
!! For the moment leni_factorize = lena_factorize because of BTF structure
      integer(long) :: lena_factorize  = 0 ! Size for factorize (real array)
      integer(long) :: leni_factorize  = 0 ! Size for factorize (integer array)
      integer(long) :: len_factorize  = 0! Only kept for compatibility
!                      with earlier version of hsl_ma48. Set to maximum
!                      of lena_factorize and leni_factorize
      integer(long) :: drop   = 0 ! Number of entries dropped
      integer :: rank   = 0  ! Estimated rank
      integer :: stat   = 0  ! STAT value after allocate failure
   end type ma48_finfo

   type ma48_sinfo
      integer :: flag   = 0 ! Flags success or failure case
      integer :: more   = 0  ! More information on failure
      integer :: stat   = 0  ! STAT value after allocate failure
   end type ma48_sinfo

contains

   SUBROUTINE ma48_initialize_double(factors,control)
      type(ma48_factors), intent(out), optional :: factors
      type(ma48_control), intent(out), optional :: control

      if (present(factors)) then
        factors%n = 0
        factors%first = 0
        factors%partial = 0
      end if
      if (present(control)) then
          control%switch = 0.5d0
          control%u      = 0.01d0
          control%drop   = 0.0d0
          control%tolerance = 0.0d0
          control%cgce      = 0.5d0
          control%lp = 6
          control%wp = 6
          control%mp = 6
          control%ldiag = 2
          control%pivoting = 3
          control%diagonal_pivoting = .false.
          control%fill_in = 3
          control%maxit = 10
          control%struct = .false.
          control%factor_blocking = 32
          control%solve_blas  = 2
          control%btf = 1
          control%multiplier = 2.0d0
          control%switch_mode = .false.
      end if
    end subroutine ma48_initialize_double

   SUBROUTINE ma48_analyse_double(matrix,factors,control,ainfo,finfo, &
                                  perm,endcol)
      type(zd11_type), Intent(in) :: matrix
      type(ma48_factors), intent(inout) :: factors
      type(ma48_control), intent(in) :: control
      type(ma48_ainfo), intent(out) :: ainfo
      type(ma48_finfo), intent(out), optional :: finfo
      integer, intent(in), optional :: perm(matrix%m+matrix%n) ! Input perm
      integer, intent(in), optional :: endcol(matrix%n) ! Define last cols

      integer, allocatable :: iwork(:)
      integer :: i,job,k,lkeep,m,n,stat,icntl(20)
      integer(long) :: la,ne,info(20)
      integer(long), parameter :: lone = 1
      real(wp):: rinfo(10),cntl(10)

      cntl(1) = control%switch
      cntl(2) = control%u
      cntl(3) = control%drop
      cntl(4) = control%tolerance
! cntl(6) set so that storage is always increased
      cntl(6) = max(1.2_wp,control%multiplier)
      icntl(1) = -1
      icntl(2) = -1
      if (control%ldiag.gt.2) icntl(2) = control%mp
      icntl(3) = control%ldiag
      icntl(4) = control%pivoting
      icntl(5) = control%factor_blocking
      icntl(6) = control%btf
      if (control%struct) then
! Abort run if matrix singular
        icntl(7) = 1
      else
        icntl(7) = 0
      end if

      m = matrix%m
      n = matrix%n
      ne = matrix%ne

      if(matrix%m .lt. 1) then
         ainfo%flag = -1
         ainfo%more = matrix%m
         if (control%ldiag>0 .and. control%lp>=0 ) &
         write (control%lp,'(/a,i3/a,i12)') &
         'Error return from MA48_ANALYSE with ainfo%flag = ',ainfo%flag, &
         'MATRIX%M has the value', matrix%m
       return
      end if

      if(matrix%n .lt. 1) then
         if (control%ldiag>0 .and. control%lp>=0 ) &
         ainfo%flag = -2
         ainfo%more = matrix%n
         write (control%lp,'(/a,i3/a,i12)') &
         'Error return from MA48_ANALYSE with ainfo%flag = ',ainfo%flag, &
         'MATRIX%N has the value', matrix%n
       return
      end if

      if(matrix%ne .lt. 0) then
         if (control%ldiag>0 .and. control%lp>=0 ) &
         ainfo%flag = -3
         ainfo%more = matrix%ne
         write (control%lp,'(/a,i3/a,i12)') &
         'Error return from MA48_ANALYSE with ainfo%flag = ',ainfo%flag, &
         'MATRIX%NE has the value', matrix%ne
       return
      end if

      ainfo%flag = 0
      ainfo%more = 0
      ainfo%stat = 0

      if(matrix%ne .eq. 0) then
        ainfo%ops = 0.0d0
        ainfo%rank = 0
        ainfo%drop = 0
        factors%ndrop = 0
        ainfo%oor = 0
        ainfo%dup = 0
        ainfo%lena_analyse = 0
        ainfo%len_analyse = 0
        ainfo%lena_factorize = 0
        ainfo%len_factorize = 0
        ainfo%struc_rank = 0
        ainfo%rank = 0
        ainfo%ncmpa = 0
        factors%first = 1
        factors%lareq = 0
        factors%lirnreq = 0
        factors%partial = 0
        factors%m = m
        factors%n = n
        if (control%struct) then
          ainfo%flag = -5
           if (control%ldiag>0 .and. control%lp>=0 ) &
                 write (control%lp,'(/a,i3/a,i5)') &
         'Error return from MA48_ANALYSE with ainfo%flag = ',ainfo%flag, &
                 'Matrix is structurally singular with rank ',ainfo%struc_rank
        else
          ainfo%flag = 4
          if (control%ldiag>1 .and. control%wp>=0) &
          write (control%wp,'(/a/a,i5)') &
              'Warning from MA48_ANALYSE: ', &
              'ainfo%flag is equal to ',ainfo%flag
        endif
        return
      endif

      lkeep = m+5*n+4*n/icntl(6)+7

      if(allocated(factors%keep)) then
!!! Problem with Fortran 95 extension
!        if(size(factors%keep,kind=long)/=lkeep) then
         if(size(factors%keep)*lone/=lkeep) then
            deallocate(factors%keep,stat=stat)
            allocate(factors%keep(lkeep),stat=stat)
            if (stat/=0) go to 100
         end if
      else
         allocate(factors%keep(lkeep),stat=stat)
         if (stat/=0) go to 100
      end if

      la = max(2*ne,(control%fill_in+0_long)*ne)

!! Can allocate to different lengths
      if(allocated(factors%val)) then
!!!! Same problem
!        if(la /= size(factors%val,kind=long)) then
         if(la /= size(factors%val)*lone) then
            deallocate(factors%irn,factors%jcn,factors%val,stat=stat)
            allocate(factors%irn(2*ne),factors%jcn(ne),factors%val(2*ne), &
               stat=stat)
            if (stat/=0) go to 100
         end if
      else
         allocate(factors%irn(2*ne),factors%jcn(ne),factors%val(2*ne),stat=stat)
         if (stat/=0) go to 100
      end if

      if (present(perm)) then
! Check permutation
         allocate (iwork(max(m,n)),stat=stat)
         if (stat/=0) go to 100

         iwork = 0
         do i = 1,m
           k = perm(i)
           if (k.gt.m .or. k.lt.0 .or. iwork(k) .ne. 0) then
             ainfo%flag = -6
             ainfo%more = i
             if (control%ldiag>0 .and. control%lp>=0 ) &
               write (control%lp,'(/a,i3/a,i12,a,i12)') &
         'Error return from MA48_ANALYSE with ainfo%flag = ',ainfo%flag, &
         'Invalid row permutation'
               deallocate (iwork,stat=stat)
             return
           end if
           iwork(k) = i
         end do

         iwork = 0
         do i = 1,n
           k = perm(m+i)
           if (k.gt.n .or. k.lt.0 .or. iwork(k) .ne. 0) then
             ainfo%flag = -6
             ainfo%more = i
             if (control%ldiag>0 .and. control%lp>=0 ) &
               write (control%lp,'(/a,i3/a,i12,a,i12)') &
         'Error return from MA48_ANALYSE with ainfo%flag = ',ainfo%flag, &
         'Invalid column permutation'
               deallocate (iwork,stat=stat)
             return
           end if
           iwork(k) = i
         end do

         deallocate (iwork,stat=stat)
         if (stat/=0) go to 100

         job = 2
      else
         job = 1
      end if

      if (control%diagonal_pivoting) job = 3

      icntl(8) = 0
      if (present(endcol)) then
        icntl(8) = 1
        factors%partial = count(endcol(1:n) == 0)
      end if

      if (present(perm)) factors%keep(1:m+n) = perm(1:m+n)

      factors%irn(1:ne) = matrix%row(1:ne)
      factors%jcn(1:ne) = matrix%col(1:ne)
      factors%val(1:ne) = matrix%val(1:ne)

!!!!!!!!!!!!!!!  Call to MA48AD !!!!!!!!!!!!!!!!!!!!!!
      if (present(endcol)) then
        call ma48ad(m,n,ne,job,la,factors%val,factors%irn,factors%jcn, &
                    factors%keep,cntl,icntl,info,rinfo,endcol)
      else
        call ma48ad(m,n,ne,job,la,factors%val,factors%irn,factors%jcn, &
                    factors%keep,cntl,icntl,info,rinfo)
      endif
! Jump if failure in allocations in MA48/MA50
       if (info(1).eq.-10) go to 100
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (info(1).eq.-4) then
        ainfo%flag = -5
        ainfo%struc_rank = info(10)

        if (control%ldiag>0 .and. control%lp>=0 ) &
              write (control%lp,'(/a,i3/a,i5)') &
         'Error return from MA48_ANALYSE with ainfo%flag = ',ainfo%flag, &
                 'Matrix is structurally singular with rank ',info(10)
        return
      endif

      if (info(1).gt.0) then
        if (info(1).eq.1 .and. info(11).gt.0) ainfo%flag = ainfo%flag + 2
        if (info(1).eq.1 .and. info(12).gt.0) ainfo%flag = ainfo%flag + 1
        if (info(1).eq.2) ainfo%flag = ainfo%flag + 4
        if (info(1).eq.4) ainfo%flag = ainfo%flag + 8
      endif

      if (ainfo%flag.gt.0) then
        if (control%ldiag>1 .and. control%wp>=0) &
          write (control%wp,'(/a/a,i5)') &
              'Warning from MA48_ANALYSE: ', &
              'ainfo%flag is equal to ',ainfo%flag
      end if

      factors%m = m
      factors%n = n
!! Both the same because of data structure in case of BTF
      factors%lareq = info(4)
      factors%lirnreq = info(4)
      factors%first = 1

      ainfo%ops    = rinfo(1)
      ainfo%rank   = info(5)
      ainfo%drop   = info(6)
      factors%ndrop   = info(6)
      ainfo%oor    = info(12)
      ainfo%dup    = info(11)
      ainfo%lena_analyse = info(3)
      ainfo%lenj_analyse = info(14)
      ainfo%len_analyse  = max(ainfo%lena_analyse,ainfo%lenj_analyse)
!! Both the same because of data structure in case of BTF
      ainfo%lena_factorize = info(4)
      ainfo%leni_factorize = info(4)
      ainfo%len_factorize  = max(ainfo%lena_factorize,ainfo%leni_factorize)
      ainfo%lblock = info(7)
      ainfo%sblock = info(8)
      ainfo%tblock = info(9)
      ainfo%struc_rank = info(10)
      ainfo%ncmpa  = info(2)

      if (present(finfo)) then
        call MA48_factorize(matrix,factors,control,finfo)
        if (finfo%flag<0) ainfo%flag = -7
      endif

      return

  100  ainfo%flag = -4
       ainfo%stat = stat
       if (control%ldiag>0 .and. control%lp>=0 ) &
         write (control%lp,'(/a,i3/a,i12)') &
         'Error return from MA48_ANALYSE with ainfo%flag = ',ainfo%flag, &
         'Allocate failed with STAT =',stat

   end subroutine ma48_analyse_double

   SUBROUTINE ma48_get_perm_double(factors,perm)
      type(ma48_factors), intent(in), optional :: factors
      integer, intent(out) :: perm(:)
      integer m,n

      m = factors%m
      n = factors%n

      perm(1:m+n) = factors%keep(1:m+n)

    end subroutine ma48_get_perm_double

   SUBROUTINE ma48_factorize_double(matrix,factors,control,finfo,fast,partial)
      type(zd11_type), intent(in) :: matrix
      type(ma48_factors), intent(inout) :: factors
      type(ma48_control), intent(in) :: control
      type(ma48_finfo), intent(out) :: finfo
      integer, optional, intent(in) :: fast,partial

      integer :: job,m,n
      integer(long) :: la,ne,info(20)
      integer stat  ! stat value in allocate statements

      integer(long), parameter :: lone = 1

      integer icntl(20)
      real(wp) multiplier,cntl(10),rinfo(10)

      integer (long), allocatable :: itemp(:)

! Check whether this has been preceded by a call to ma48_analyse
      if (factors%first .le. 0) then
         finfo%flag = -10
         write (control%lp,'(/a,i3/a,i12)') &
         'Error return from MA48_FACTORIZE with finfo%flag = ',finfo%flag, &
         'No prior call to MA48_ANALYSE'
         return
      endif

! Storage must be increased if ma48bd fails so effective value of 
!     multiplier is set greater than one
      multiplier = max(1.2_wp,control%multiplier)
      m = matrix%m
      n = matrix%n
      ne = matrix%ne

      if(factors%m/=matrix%m) then
       finfo%flag = -1
       finfo%more = matrix%m
       if (control%ldiag>0 .and. control%lp>=0 ) &
         write (control%lp,'(/a,i3/a,i12,a,i12)') &
         'Error return from MA48_FACTORIZE with finfo%flag = ',finfo%flag, &
         'MATRIX%M has the value', matrix%m,' instead of',factors%m
       return
      end if

      if(factors%n/=matrix%n) then
        finfo%flag = -2
        finfo%more = matrix%n
          if (control%ldiag>0 .and. control%lp>=0 ) &
            write (control%lp,'(/a,i3/a,i12,a,i12)') &
           'Error return from MA48_FACTORIZE with finfo%flag = ',finfo%flag, &
           'MATRIX%N has the value', matrix%n,' instead of',factors%n
        return
      end if

      if(matrix%ne .lt. 0) then
         finfo%flag = -3
         finfo%more = matrix%ne
         if (control%ldiag>0 .and. control%lp>=0 ) &
            write (control%lp,'(/a,i3/a,i12)') &
           'Error return from MA48_FACTORIZE with finfo%flag = ',finfo%flag, &
           'MATRIX%NE has the value', matrix%ne
       return
      end if

      finfo%flag = 0
      finfo%more = 0
      finfo%stat = 0

      if(matrix%ne .eq. 0) then
        finfo%ops = 0.0d0
        finfo%rank = 0
        finfo%drop = 0
        finfo%lena_factorize = 0
        finfo%leni_factorize = 0
        finfo%len_factorize = 0
        finfo%size_factor = 0
        factors%first = 2
        factors%lareq = 0
        factors%lirnreq = 0
        factors%partial = 0
        finfo%flag = 4
        if (control%ldiag>1 .and. control%wp>=0) &
          write (control%wp,'(/a/a,i5)') &
          'Warning from MA48_FACTORIZE: ', &
          'finfo%flag is equal to ',finfo%flag
        return
      endif

! Check length of main array

!!! Problem with Fortran extension
!     la = size(factors%val,kind=long)
      la = size(factors%val)*lone

      if (la<factors%lareq) then
            la = factors%lareq
            deallocate (factors%val,stat=stat)
            allocate (factors%val(la),stat=stat)
            if (stat/=0) go to 100
! Copy irn
            allocate (itemp(la),stat=stat)
            if (stat/=0) go to 100
            itemp(1:ne) = factors%irn(1:ne)
            deallocate (factors%irn,stat=stat)
            allocate (factors%irn(la),stat=stat)
            if (stat/=0) go to 100
            factors%irn(1:ne) = itemp(1:ne)
            deallocate (itemp,stat=stat)
      endif

      cntl(1) = control%switch
      cntl(2) = control%u
      cntl(3) = control%drop
      cntl(4) = control%tolerance
! Switch off error message printing from ma48bd
      icntl(1) = -1
      icntl(2) = -1
      if (control%ldiag.gt.2) icntl(2) = control%mp
      icntl(3) = control%ldiag
      icntl(4) = control%pivoting
      icntl(5) = control%factor_blocking
      icntl(6) = control%btf
! Immediate return from MA48B/BD if LA is too small.
! Will make this not operate
      icntl(10) = 1
! Switch to slow mode if necessary?
      icntl(11) = 0
      if (control%switch_mode) icntl(11) = 1

      job = 1
! Can only have job=2 if there was a previous factorize with no dropping
      if (present(fast) .and. factors%ndrop==0 .and. factors%first.gt.1) &
          job = 2
      icntl(8) = 0
      if (present(partial) .and. factors%ndrop==0) then
        if(factors%partial.gt.0) then ! Note: can't test on above line as
                                      ! may not have been set.
          job = 3
          icntl(8) = 1
        endif
      endif

! Copy matrix for factorization
 101  factors%val(1:ne) = matrix%val(1:ne)

!!!!!!!!!!!!!!!  Call to MA48BD !!!!!!!!!!!!!!!!!!!!!!
      call ma48bd(factors%m,factors%n,ne,job,la,factors%val, &
                  factors%irn,factors%jcn, &
                  factors%keep,cntl,icntl, &
                  info,rinfo)
! Jump if failure in allocations in MA48/MA50
      if (info(1).eq.-10) go to 100
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (info(1)==-3) then
! Increase size for factors and return with job = 1 call
         la = real(la,kind=wp)*multiplier
         deallocate (factors%val,stat=stat)
         allocate (factors%val(la),stat=stat)
         if (stat/=0) go to 100
! Copy irn
         allocate (itemp(la),stat=stat)
         if (stat/=0) go to 100
         itemp(1:ne) = factors%irn(1:ne)
         deallocate (factors%irn,stat=stat)
         allocate (factors%irn(la),stat=stat)
         if (stat/=0) go to 100
         factors%irn(1:ne) = itemp(1:ne)
         deallocate (itemp,stat=stat)

         job = 1
         goto 101

      end if

      if (info(1)==-7) then
         finfo%flag = -11
         if (control%ldiag>0 .and. control%lp>=0 ) &
            write (control%lp,'(/a,i3/a,i12)') &
            'Error return from MA48_FACTORIZE with finfo%flag = ',&
            finfo%flag,' Matrix entries unsuitable for fast factorization'
         return
      end if

      if (info(1).eq.2)    finfo%flag = finfo%flag + 4

      if (info(1)>=0) then
        finfo%lena_factorize = info(4)
        finfo%leni_factorize = info(4)
        finfo%len_factorize = info(4)
        factors%lareq = info(4)
        factors%lirnreq = info(4)
        finfo%rank   = info(5)
        finfo%drop   = info(6)
        finfo%ops    = rinfo(1)
        call nonzer(m,n,factors%keep,info)
        finfo%size_factor = info(2)+info(3)+info(4)
        factors%first = 2
      end if

      if (finfo%flag.gt.0) then
        if (control%ldiag>1 .and. control%wp>=0) &
          write (control%wp,'(/a/a,i5)') &
              'Warning from MA48_FACTORIZE: ', &
              'finfo%flag is equal to ',finfo%flag
      end if

      return

  100  finfo%flag = -4
       finfo%stat = stat
       if (control%ldiag>0 .and. control%lp>=0 ) &
         write (control%lp,'(/a,i3/a,i12)') &
         'Error return from MA48_FACTORIZE with finfo%flag = ',finfo%flag, &
         'Allocate failed with STAT =',stat

   end subroutine ma48_factorize_double

   SUBROUTINE ma48_solve_double(matrix,factors,rhs,x,control,sinfo,trans, &
                         resid,error)
      type(zd11_type), intent(in) :: matrix
      type(ma48_factors), intent(in) :: factors
      real(wp), intent(in) :: rhs(:)
      real(wp), intent(out) :: x(:)
      type(ma48_control), intent(in) :: control
      type(ma48_sinfo), intent(out) :: sinfo
      integer, optional, intent(in) :: trans
      real(wp), optional, intent(out) :: resid(2)
      real(wp), optional, intent(out) :: error
      integer icntl(20),job,m,n,stat
      integer(long) :: lasolve,info(20)
      integer(long), parameter :: lone = 1
      real(wp) cntl(10),err(3)
      logical trans48

      integer, allocatable :: iwork(:)
      real(wp), allocatable :: work(:)
      real(wp), allocatable :: rhswork(:)

! Check whether this has been preceded by a call to ma48_factorize
      if (factors%first .le. 1) then
         sinfo%flag = -10
         if (control%ldiag>0 .and. control%lp>=0 ) &
         write (control%lp,'(/a,i3/a,i12)') &
         'Error return from MA48_SOLVE with sinfo%flag = ',sinfo%flag, &
         'No prior call to MA48_FACTORIZE'
         return
      endif

      m = matrix%m
      n = matrix%n

      if(factors%m/=matrix%m) then
         sinfo%flag = -1
         sinfo%more = matrix%m
         if (control%ldiag>0 .and. control%lp>=0 ) &
           write (control%lp,'(/a,i3/a,i12,a,i12)') &
          'Error return from MA48_SOLVE with sinfo%flag = ',sinfo%flag, &
          'MATRIX%M has the value', matrix%m,' instead of',factors%m
       return
      end if

      if(factors%n/=matrix%n) then
         sinfo%flag = -2
         sinfo%more = matrix%n
         if (control%ldiag>0 .and. control%lp>=0 ) &
           write (control%lp,'(/a,i3/a,i12,a,i12)') &
          'Error return from MA48_SOLVE with sinfo%flag = ',sinfo%flag, &
          'MATRIX%N has the value', matrix%n,' instead of',factors%n
       return
      end if

      if(matrix%ne .lt. 0) then
         sinfo%flag = -3
         sinfo%more = matrix%ne
         if (control%ldiag>0 .and. control%lp>=0 ) &
         write (control%lp,'(/a,i3/a,i12)') &
          'Error return from MA48_SOLVE with sinfo%flag = ',sinfo%flag, &
          'MATRIX%NE has the value', matrix%ne
       return
      end if


      sinfo%flag = 0
      sinfo%more = 0
      sinfo%stat = 0

      if(matrix%ne .eq. 0) then
        x = 0.0d0
        if (present(resid)) resid = 0.0d0
        if (present(error)) error = 0.0d0
        return
      endif

      trans48 = present(trans)

      allocate (iwork(max(m,n)),work(max(3*m+n,3*n+m)),stat=stat)
      if (stat/=0) go to 100

      cntl(5)  = control%cgce
      icntl(1) = -1
      icntl(2) = -1
      if (control%ldiag.gt.2) icntl(2) = control%mp
      icntl(3) = control%ldiag
      icntl(5) = 0
      if (control%solve_blas.gt.1) icntl(5) = 2
      icntl(9) = control%maxit
      stat = 0

      job = 1
      if (control%maxit .gt. 0) then
        if (present(resid)) then
          if (present(error)) then
            job = 4
          else
            if (control%maxit == 1) then
              job = 2
            else
              job = 3
            end if
          end if
        end if
      end if

      if (trans48) then
        allocate (rhswork(n),stat=stat)
        if (stat/=0) go to 100
        rhswork = rhs(1:n)
      else
        allocate (rhswork(m),stat=stat)
        if (stat/=0) go to 100
        rhswork = rhs(1:m)
      endif

!!! Problem with kind = long is size statement
      lasolve = size(factors%val)*lone
      call ma48cd(factors%m,factors%n,trans48,job,  &
                  lasolve,factors%val, &
                  factors%irn,factors%keep,cntl,icntl,rhswork,   &
                  x,err,info)
! Jump if failure in allocations in MA48/MA50
       if (info(1).eq.-10) go to 100

      sinfo%flag = info(1)

!!!
! Only invoked if error return -8 or -9 from ma48cd
      if(sinfo%flag .lt. 0) then
         if (control%ldiag>0 .and. control%lp>=0 ) &
         write (control%lp,'(/a,i3)') &
          'Error return from MA48_SOLVE with sinfo%flag = ',sinfo%flag
       return
      end if

      if (job>1) resid(1) = err(1)
      if (job>1) resid(2) = err(2)
      if (job>3) error    = err(3)

      deallocate (iwork,work,rhswork,stat=stat)
      return

  100  sinfo%flag = -4
       sinfo%stat = stat
       if (control%ldiag>0 .and. control%lp>=0 ) &
         write (control%lp,'(/a,i3/a,i12)') &
         'Error return from MA48_SOLVE with sinfo%flag = ',sinfo%flag, &
         'Allocate failed with STAT =',stat

   end subroutine ma48_solve_double

   SUBROUTINE ma48_finalize_double(factors,control,info)
      type(ma48_factors), intent(inout) :: factors
      type(ma48_control), intent(in) :: control
      integer, intent(out) :: info

      info = 0

      if (allocated(factors%keep)) deallocate(factors%keep,stat=info)
      if (info==0 .and. allocated(factors%val))  deallocate &
                       (factors%irn,factors%jcn,factors%val,stat=info)
      if (info==0) return

      if (control%ldiag>0 .and. control%lp>=0 ) write (control%lp,'(/2a,i5)') &
         'Error return from MA48_finalize: ',&
         'deallocate failed with STAT =',info

    end subroutine ma48_finalize_double

      SUBROUTINE nonzer(M,N,KEEP,INFO)
!     NONZER should be called after MA48BD ..
! M,N,KEEP .. the same as for that call
! INFO(1)  Number of original entries excluding duplicates and out-of-range
! INFO(2)  Number of entries in off-diagonal blocks of block triangular form
! INFO(3)  Number of entries in triangular blocks on diagonal
! INFO(4)  Number of entries in L/U decomp of non-triangular blocks on diagonal
! INFO(3) and INFO(4) are separate because they are stored separately .. the
! INFO(3) entries are also used to calculate residuals.
      INTEGER M,N
      integer(long) :: keep(*)
      INTEGER(LONG) :: INFO(4)
! Declare internal variables
      INTEGER IPTRL,IPTRU,IPTRD,IPTRO,NBLOCK,MBLOCK,KBLOCK
      INTEGER NB,KB,JB,J1,J2,NC,NR,LC
! Pointers for KEEP array
      IPTRL = M+N
      IPTRU = IPTRL + N
      IPTRD = IPTRU + N
      IPTRO = IPTRD + N + 1
      NBLOCK = IPTRO + N - 1
      MBLOCK = NBLOCK + 1
      KBLOCK = MBLOCK + 1

      NB = KEEP(KBLOCK+3)

! Number of entries in original matrix with dups and o-o-r omitted
      INFO(1) = KEEP(IPTRO+N+1) - 1

! Number entries in off-diagonal blocks
      INFO(2) = KEEP(IPTRO+N+1) - KEEP(IPTRO+1)

! Number of entries in triangular diagonal blocks
      INFO(3) = 0
      J2 = 0
      LC = 0
      DO 100 JB = 1,NB
! NC is number of columns in block
        NC = KEEP(NBLOCK+3*JB)
        J1 = J2 + 1
        J2 = J1 + NC - 1
! MBLOCK negative flags triangular block
        IF (KEEP(MBLOCK+3*JB).LT.0) THEN
          INFO(3) = INFO(3) + KEEP(IPTRD+J2+1) - KEEP(IPTRD+J1)
        ELSE
! Recording information on last non-triangular block
          KB = JB
          NR = NC
          LC = J2
        ENDIF
  100 CONTINUE

! Compute number of entries in factors of non-triangular blocks
      IF (LC.EQ.0) THEN
! No non-triangular blocks
        INFO(4) = 0
      ELSE
        NC = NR
        IF (NB.EQ.1) THEN
! Only one block ...  matrix may be rectangular
          NR = M
          INFO(4) = 0
        ELSE
! Set to total number of entries in previous non-triangular blocks
          INFO(4) = KEEP(KBLOCK+3*KB) - 1
        ENDIF
! Add number of entries for last non-triangular block
! Sparse and full part
        INFO(4) = INFO(4) + KEEP(IPTRL+LC) +  &
                  MAX(((NC-KEEP(MBLOCK+3*KB))+(NR-KEEP(MBLOCK+3*KB))), &
                      ((NC-KEEP(MBLOCK+3*KB))*(NR-KEEP(MBLOCK+3*KB))))
      ENDIF

      end subroutine nonzer


  SUBROUTINE ma48_special_rows_and_cols_double(factors,rank,rows,cols,  &
                                               control,info)
      type(ma48_factors), intent(in) :: factors
      integer,intent(out) :: rank,info
      integer,intent(out),dimension(factors%m) :: rows
      integer,intent(out),dimension(factors%n) :: cols
      type(ma48_control), intent(in) :: control
      integer(long) :: la
      integer(long), parameter :: lone=1
      integer, allocatable :: iwork(:)

      allocate (iwork(max(factors%m,factors%n)),stat=info)
      if (info/=0) then
        info = -1

      if (control%ldiag>0 .and. control%lp>=0 ) write (control%lp,'(/2a,i5)') &
         'Error return from MA48_finalize: ',&
         'allocate failed with STAT =',info

        return
      end if
!!!!!
!     la = size(factors%val,kind=long)
      la = size(factors%val)*lone
      call ma51ad(factors%m,factors%n,la,factors%irn,factors%keep,rank, &
                  rows,cols,iwork)
      deallocate (iwork, stat=info)
   end subroutine ma48_special_rows_and_cols_double

   SUBROUTINE ma48_determinant_double(factors,sgndet,logdet,control,info)
      type(ma48_factors), intent(in) :: factors
      integer,intent(out) :: sgndet,info
      real(wp),intent(out) :: logdet
      type(ma48_control), intent(in) :: control
      integer(long) :: la
      integer(long), parameter :: lone = 1
      integer, allocatable :: iwork(:)
      allocate (iwork(factors%n),stat=info)
      if (info/=0) then
        info = -1

      if (control%ldiag>0 .and. control%lp>=0 ) write (control%lp,'(/2a,i5)') &
         'Error return from MA48_finalize: ',&
         'allocate failed with STAT =',info

        return
      end if
!!!!!
!     la = size(factors%val,kind=long)
      la = size(factors%val)*lone
      call ma51cd(factors%m,factors%n,la,factors%val,factors%irn, &
                  factors%keep,sgndet,logdet,iwork)
      deallocate (iwork, stat=info)
  end subroutine ma48_determinant_double

end module hsl_ma48_double
