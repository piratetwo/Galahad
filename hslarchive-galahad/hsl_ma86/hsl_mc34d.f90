! COPYRIGHT (c) 2009 Council for the Central Laboratory
!               of the Research Councils
! Original date 20 October 2009. Version 1.0.0.

! Fortran 95 version of the mc34 package.

! 18 May 2010 Version 1.1.0 -jhogg
!             Create hsl_mc34_integer
!             Change from logical Hermitian to integer sym_type to cope with
!             skew symmetric matrices as well as Hermitian and symmetric ones.

! to change precision:
!    change _double, kind(0.0d0)
! For complex version:
!    change real to complex
! For Hermitian case, take conjugate for upper triangle entries
! For integer version:
!    change real to integer, kind(0)

   module hsl_mc34_double
   implicit none
   private
   public mc34_expand

   integer, parameter :: wp = kind(0.0d0)

   interface mc34_expand
      module procedure mc34_expand_double
   end interface

   contains

      subroutine mc34_expand_double(n,row,ptr,iw,a,sym_type)

!  this subroutine generates the expanded structure for a
!   matrix a with a symmetric sparsity pattern given the structure 
!   for the lower triangular part.  diagonal entries need not be present.

      integer, intent(in) :: n  ! holds the order of a.

      integer, intent(inout) :: row(*) ! must be set by the user to
!       hold the row indices of the lower triangular part of a.
!       the entries of a single column must be
!       contiguous. the entries of column j must precede those of column
!       j+1, and there must be no wasted space between
!       columns. row indices within a column may be in any order.  on
!       exit, it will have the same meaning but will be changed to hold
!       the row indices of the entries in the expanded structure.  diagonal
!       entries need not be present. the new row indices added in the
!       upper triangular part will be in order for each column and will
!       precede the row indices for the lower triangular part which will
!       remain in the input order.

     integer, intent(inout) ::ptr(n+1)  !  must be set
!       by the user so that ptr(j) is the position in row
!       of the first entry in column j and
!       ptr(n+1) must be set to one more than the total number of
!       entries.  on exit, ptr(j) will have the same meaning but
!       will be changed to point to the position of the first entry of
!       column j in the expanded structure. the new value of
!       ptr(n+1) will be one greater than the number of entries in
!       the expanded structure.

     integer :: iw(n) ! workspace

     real(wp), optional, intent(inout) :: a(*) 
!       if present, a(1:ptr(n+1)-1) must be set by the user so that
!       a(k) holds the value of the entry in row(k). 
!       on exit, a will hold the values of the entries in the expanded 
!       structure corresponding to the output values of row.

     integer, optional, intent(in) :: sym_type 
!      if present with value 1, matrix is skew symmetric.
!      if present with value 2, matrix is hermitian.
!      otherwise matrix is symmetric.

      integer :: ckp1 ! used as running pointer
      integer :: i,i1,i2,ii,ipkp1,ipos
      integer :: j,jstart 
      integer :: lenk ! number of entries in col. j of original structure
      integer :: ndiag ! number diagonal entries present
      integer :: newtau ! number of entries in expanded storage
      integer :: oldtau ! number of entries in symmetric storage
      integer :: r_sym_type ! real sym_type value (used as argument is optional)

      oldtau = ptr(n+1) - 1
      iw(1:n) = 0

! iw(j) set to total number entries in col. j of expanded mx.
      ndiag = 0
      do j = 1,n
        i1 = ptr(j)
        i2 = ptr(j+1) - 1
        iw(j) = iw(j) + i2 - i1 + 1
        do ii = i1,i2
          i = row(ii)
          if (i /= j) then
            iw(i) = iw(i) + 1
          else
            ndiag = ndiag + 1
          end if
        end do
      end do

      newtau = 2*oldtau - ndiag
! ipkp1 points to position  after end of column being currently processed
      ipkp1 = oldtau + 1
! ckp1 points to position  after end of same column in expanded structure
      ckp1 = newtau + 1
! go through the array in the reverse order placing lower triangular
!     elements in  appropriate slots.
      do j = n,1,-1
        i1 = ptr(j)
        i2 = ipkp1
        lenk = i2 - i1
! jstart is running pointer to position in new structure
        jstart = ckp1
! set ikp1 for next column
        ipkp1 = i1
        i2 = i2 - 1
! run through columns in reverse order
! lower triangular part of col. moved to end of same column in expanded form
        if (present(a)) then
          do ii = i2,i1,-1
            jstart = jstart - 1
            a(jstart) = a(ii)
            row(jstart) = row(ii)
          end do
        else
          do ii = i2,i1,-1
            jstart = jstart - 1
            row(jstart) = row(ii)
          end do
        end if
! ptr is set to position of first entry in lower triangular part of
!     column j in expanded form
        ptr(j) = jstart
! set ckp1 for next column
        ckp1 = ckp1 - iw(j)
! reset iw(j) to number of entries in lower triangle of column.
        iw(j) = lenk
      end do
 
! again sweep through the columns in the reverse order, this
!     time when one is handling column j the upper triangular
!     elements a(j,i) are put in position.
        do j = n,1,-1
          i1 = ptr(j)
          i2 = ptr(j) + iw(j) - 1
! run down column in order
! note that i is always greater than or equal to j
          if (present(a)) then
            r_sym_type = 0 ! symmetric
            if(present(sym_type)) r_sym_type = sym_type
            select case(r_sym_type)
            case(1) ! skew symmetric
              do ii = i1,i2
                i = row(ii)
                if (i == j) cycle
                ptr(i) = ptr(i) - 1
                ipos = ptr(i)
                a(ipos) = -a(ii)
                row(ipos) = j
              end do
            case default ! symmetric or hermitian
              do ii = i1,i2
                i = row(ii)
                if (i == j) cycle
                ptr(i) = ptr(i) - 1
                ipos = ptr(i)
                a(ipos) = a(ii)
                row(ipos) = j
              end do
            end select
          else
            do ii = i1,i2
              i = row(ii)
              if (i == j) cycle
              ptr(i) = ptr(i) - 1
              ipos = ptr(i)
              row(ipos) = j
            end do
          end if
        end do
      ptr(n+1) = newtau + 1

      end subroutine mc34_expand_double
   end module hsl_mc34_double
