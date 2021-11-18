! COPYRIGHT (c) 2007 Science and Technology Facilities Council
! Original date 21st Sept 2007

!-*-*-*-*-*-  HSL_KB22_ long_integer MODULE -*-*-*-*-*-

!  Nick Gould and Philippe Toint, for GALAHAD productions
!  Copyright reserved
!  January 26th 1995

! 21st September 2007 Version 1.0.0. Version numbering added.

      MODULE HSL_KB22_long_integer

         IMPLICIT NONE

         PRIVATE
         PUBLIC :: KB22_build_heap, KB22_get_smallest

!-----------------------------
!   L o n g   i n t e g e r
!-----------------------------

         INTEGER, PARAMETER :: long = selected_int_kind( 18 )

         INTERFACE KB22_build_heap
             MODULE PROCEDURE KB22_build_heap_long_integer
         END INTERFACE
         INTERFACE KB22_get_smallest
             MODULE PROCEDURE KB22_get_smallest_long_integer
         END INTERFACE

      CONTAINS

!-*-*-*-  H S L _ K B 2 2 _ b u i l d _ h e a p   S U B R O U T I N E   -*-*-*

         SUBROUTINE KB22_build_heap_long_integer( n, A, inform, INDA )

!  Given an array A, elements A(1), ...., A(N), subroutine KB22_build_heap
!  re-arranges the elements to form a heap in which each parent has a smaller
!  value than either of its children.

!  Algorithm 232 of CACM (J. W. J. Williams): a combination of
!  ALGOL procedures SETHEAP and INHEAP

!  ------------------------- dummy arguments --------------------------
!
!  n      integer, which gives the number of values to be sorted.
!         n must be positive
!
!  A      long integer array of length n. On input, A must contain the
!         values which are to be sorted. On output, these values
!         will have been permuted so as to form a heap
!
!  inform integer, which informs the user of the success of KB22_build_heap.
!         If inform = 0 on exit, the heap has been formed.
!         If inform = 1 on exit, n was input with a value less than
!                       or equal to 0 and the heap has not been formed.
!
!  INDA   is an OPTIONAL integer array of length n. On input, INDA may
!         be used to hold indexing information (such as the original
!         order of the values) about A. On output, INDA will have been
!         permuted so that INDA(k) still refers to A(k).
!
!  ------------------ end of dummy arguments --------------------------

         INTEGER, INTENT( IN ) :: n
         INTEGER, INTENT( OUT ) :: inform
         INTEGER ( KIND = long ), INTENT( INOUT ), DIMENSION( n ) :: A
         INTEGER, INTENT( INOUT ), OPTIONAL, DIMENSION( n ) :: INDA

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

         INTEGER :: i, j, k, indin
         INTEGER ( KIND = long ) :: rin
         LOGICAL :: index

!  Add the elements to the heap one at a time

         index = PRESENT( INDA )
         IF ( n <= 0 ) THEN
            inform = 1
            RETURN
         ENDIF

         IF ( index ) THEN
           DO k = 2, n
              rin = A( k )
              indin = INDA( k )

!  The cycle may be repeated log2(k) times, but on average is repeated
!  merely twice

              i = k
              DO
                 IF ( i <= 1 ) EXIT
                 j = i / 2
                 IF ( A( j ) <= rin ) EXIT
                 A( i ) = A( j )
                 INDA( i ) = INDA( j )
                 i = j
              END DO
              A( i ) = rin
              INDA( i ) = indin
           END DO
         ELSE
           DO k = 2, n
              rin = A( k )

!  The cycle may be repeated log2(k) times, but on average is repeated
!  merely twice

              i = k
              DO
                 IF ( i <= 1 ) EXIT
                 j = i / 2
                 IF ( A( j ) <= rin ) EXIT
                 A( i ) = A( j )
                 i = j
              END DO
              A( i ) = rin
           END DO
         END IF
         inform = 0

         RETURN

!  End of subroutine KB22_build_heap

         END SUBROUTINE KB22_build_heap_long_integer

!-*-*-  H S L _ K B 2 2 _ g e t _ s m a l l e s t  S U B R O U T I N E   -*-*-

         SUBROUTINE KB22_get_smallest_long_integer( m, A, inform, INDA )

!  Given an array A, elements A(1), ...., A(m) forming a heap,
!  KB22_get_smallest assigns to rout the value of A(1), the smallest
!  member of the heap, and arranges the remaining members as elements
!  1 to m - 1 of A. rout is then placed in A(m)

!  Algorithm 232 of CACM (J. W. J. Williams): a combination of
!  ALGOL procedures OUTHEAP and SWOPHEAP

!  ------------------------- dummy arguments --------------------------
!
!  m      integer, which gives the number of values to be sorted.
!         m must be positive
!
!  A      long integer array of length m. On input, A must contain the values
!         which are to be sorted stored in a heap. On output, the smallest
!         value will have been moved into A(m) and the remaining values A(k),
!         k = 1,..., m-1 will have been restored to a heap
!
!  inform integer, which informs the user of the success of KB22_get_smallest.
!         If inform = 0 on exit, the smallest value has been found.
!         If inform = 1 on exit, m was input with a value less than
!                       or equal to 0 and the heap has not been formed
!
!  INDA   optional integer array of length m. On input, INDA may be used
!         to hold indexing information (see KB22_build_heap) about A.
!         On output, INDA will have been permuted so that INDA(k) still
!         refers to A(k). This argument is only permitted if it was present
!         when calling KB22_build_heap
!
!  ------------------ end of dummy arguments --------------------------

         INTEGER, INTENT( IN ) :: m
         INTEGER, INTENT( OUT ) :: inform
         INTEGER ( KIND = long ), INTENT( INOUT ), DIMENSION( m ) :: A
         INTEGER, INTENT( INOUT ), OPTIONAL, DIMENSION( m ) :: INDA

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

         INTEGER :: i, j, indin, indout
         INTEGER ( KIND = long ) :: rin, rout
         LOGICAL :: index

         index = PRESENT( INDA )

!  Add the element rin to the heap, extract and assign to rout
!  the value of the smallest member of the resulting set, and
!  leave the remaining elements in a heap of the original size.
!  In this process, elements 1 to n+1 of the array A may be disturbed

         IF ( m <= 0 ) THEN
            inform = 1
            RETURN
         ENDIF

         IF ( m > 1 ) THEN

           IF ( index ) THEN
              i = 1
              rout = A( 1 )
              indout = INDA( 1 )
              rin = A( m )
              indin = INDA( m )

!  Move from the top of the heap comparing the value of node i
!  with its two daughters. If node i is smallest, the heap has been
!  restored. If one of the children is smallest, promote this child
!  in the heap and move to the now vacated node.
!  This cycle may be repeated log2(m) times

              DO
                 j = i + i
                 IF ( j > m - 1 ) EXIT

!  Determine which of the two daughters is smallest

                 IF ( A( j + 1 ) < A( j ) ) j = j + 1

!  Determine if the smaller daughter is less than the value from node i

                 IF ( A( j ) >= rin ) EXIT
                 A( i ) = A( j )
                 INDA( i ) = INDA( j )
                 i = j
              END DO

!  The heap has been restored

              A( i ) = rin
              INDA( i ) = indin

!  Store the smallest value in the now vacated m-th position of the list

              A( m ) = rout
              INDA( m ) = indout
           ELSE
              i = 1
              rout = A( 1 )
              rin = A( m )

!  Move from the top of the heap comparing the value of node i
!  with its two daughters. If node i is smallest, the heap has been
!  restored. If one of the children is smallest, promote this child
!  in the heap and move to the now vacated node.
!  This cycle may be repeated log2(m) times

              DO
                 j = i + i
                 IF ( j > m - 1 ) EXIT

!  Determine which of the two daughters is smallest

                 IF ( A( j + 1 ) < A( j ) ) j = j + 1

!  Determine if the smaller daughter is less than the value from node i

                 IF ( A( j ) >= rin ) EXIT
                 A( i ) = A( j )
                 i = j
              END DO

!  The heap has been restored

              A( i ) = rin

!  Store the smallest value in the now vacated m-th position of the list

              A( m ) = rout
           END IF

         END IF
         inform = 0

         RETURN

!  End of subroutine KB22_get_smallest

         END SUBROUTINE KB22_get_smallest_long_integer

!  End of module HSL_KB22_long_integer

      END MODULE HSL_KB22_long_integer
