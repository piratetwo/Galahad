#include <fintrf.h>
      SUBROUTINE galmxCopyPtrToInteger( px, Y, n )
      mwPointer :: px
      mwSize :: n
      MWSIZE, DIMENSION( n ) :: Y

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copy INTEGER values from Matlab pointer array to Fortran array

!  Arguments

!  px  - Pointer to ir or jc array
!  Y   - Integer Fortran array
!  n   - number of elements to copy

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 4th 2007

      INTEGER :: dig

      dig = digits( dig ) + 1
      SELECT CASE ( dig )
      CASE ( 64 )
        CALL mxCopyPtrToInteger8( px, Y, n )
      CASE ( 16 )
        CALL mxCopyPtrToInteger2( px, Y, n )
      CASE ( 8 )
        CALL mxCopyPtrToInteger1( px, Y, n )
      CASE default
        CALL mxCopyPtrToInteger4( px, Y, n )
      END SELECT

      RETURN
      END SUBROUTINE galmxCopyPtrToInteger

      SUBROUTINE galmxCopyPtrToReal( px, Y, n )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      mwPointer :: px
      mwSize :: n
      REAL ( KIND = wp ), DIMENSION( n ) :: Y

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copy REAL values from Matlab pointer array to Fortran array

!  Arguments

!  px  - Pointer to ir or jc array
!  Y   - Real Fortran array
!  n   - number of elements to copy

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 4th 2007

      REAL ( KIND = wp ) :: ddig = 1.0_wp

      SELECT CASE ( digits( ddig ) )
      CASE ( 53 )
        CALL mxCopyPtrToReal8( px, Y, n )
      CASE default
        CALL mxCopyPtrToReal4( px, Y, n )
      END SELECT

      RETURN
      END SUBROUTINE galmxCopyPtrToReal

      SUBROUTINE galmxCopyIntegerToPtr( Y, px, n )
      mwPointer :: px
      mwSize :: n
      MWSIZE, DIMENSION( n ) :: Y

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copy INTEGER values from Fortran array to Matlab pointer array

!  Arguments

!  Y   - Integer Fortran array
!  px  - Pointer to ir or jc array
!  n   - number of elements to copy

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 4th 2007

      INTEGER :: dig

      dig = digits( dig ) + 1
      SELECT CASE ( dig )
      CASE ( 64 )
        CALL mxCopyInteger8ToPtr( Y, px, n )
      CASE ( 16 )
        CALL mxCopyInteger2ToPtr( Y, px, n )
      CASE ( 8 )
        CALL mxCopyInteger1ToPtr( Y, px, n )
      CASE default
        CALL mxCopyInteger4ToPtr( Y, px, n )
      END SELECT

      RETURN
      END SUBROUTINE galmxCopyIntegerToPtr

      SUBROUTINE galmxCopyRealToPtr( Y, px, n )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      mwPointer :: px
      mwSize :: n
      REAL ( KIND = wp ), DIMENSION( n ) :: Y

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copy REAL values from Fortran array to Matlab pointer array

!  Arguments

!  Y   - Real Fortran array
!  px  - Pointer to ir or jc array
!  n   - number of elements to copy

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 4th 2007

      REAL ( KIND = wp ) :: ddig = 1.0_wp

      SELECT CASE ( digits( ddig ) )
      CASE ( 53 )
        CALL mxCopyReal8ToPtr( Y, px, n )
      CASE default
        CALL mxCopyReal4ToPtr( Y, px, n )
      END SELECT

      RETURN
      END SUBROUTINE galmxCopyRealToPtr

      SUBROUTINE galmxCopyLogicalToPtr( Y, px, n )
      mwPointer :: px
      mwSize :: n
      LOGICAL, DIMENSION( n ) :: Y

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copy LOGICAL values from Fortran array to Matlab pointer array

!  Arguments

!  Y   - Logical Fortran array
!  px  - Pointer to ir or jc array
!  n   - number of elements to copy

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 4th 2007

      INTEGER :: dig
      INTEGER :: LY( n )

      DO i = 1, n
        IF ( Y( i ) ) THEN
          LY( i ) = 1
        ELSE
          LY( i ) = 0
        END IF
      END DO
      CALL galmxCopyIntegerToPtr( LY, px, n )

      RETURN
      END SUBROUTINE galmxCopyLogicalToPtr

      FUNCTION galmxCreateInteger( )
      mwPointer :: galmxCreateInteger

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Create an unpopulated INTEGER pointer

!  No arguments

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 4th 2007

      INTEGER :: dig
      mwSize :: m

      m = 1
      dig = digits( dig ) + 1
      SELECT CASE ( dig )
      CASE ( 64 )
        galmxCreateInteger = mxCreateNumericMatrix( 1, m,
     *                         mxClassIDFromClassName('int64'), 0 ) 
      CASE ( 16 )
        galmxCreateInteger = mxCreateNumericMatrix( 1, m,
     *                         mxClassIDFromClassName('int16'), 0 ) 
      CASE ( 8 )
        galmxCreateInteger = mxCreateNumericMatrix( 1, m,
     *                         mxClassIDFromClassName('int8'),  0 ) 
      CASE default
        galmxCreateInteger = mxCreateNumericMatrix( 1, m,
     *                         mxClassIDFromClassName('int32'), 0 ) 
      END SELECT

      RETURN
      END FUNCTION galmxCreateInteger

      FUNCTION galmxCreateReal( )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      mwPointer :: galmxCreateReal

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Create an unpopulated REAL pointer

!  No arguments

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 4th 2007


      REAL ( KIND = wp ) :: ddig = 1.0_wp
      mwSize :: m

      m = 1
      SELECT CASE ( digits( ddig ) )
      CASE ( 53 )
        galmxCreateReal = mxCreateNumericMatrix( 1, m,
     *                      mxClassIDFromClassName('double'), 0 ) 
      CASE default
        galmxCreateReal = mxCreateNumericMatrix( 1, m,
     *                      mxClassIDFromClassName('single'), 0 ) 
      END SELECT

      RETURN
      END FUNCTION galmxCreateReal

      FUNCTION galmxCreateIntegerMatrix( m, n )
      mwPointer :: galmxCreateIntegerMatrix
      mwSize :: m, n

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Create an unpopulated INTEGER pointer matrix

!  Arguments

!  m - number of rows
!  n - number of columns

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 4th 2007

      INTEGER :: dig

      dig = digits( dig ) + 1
      SELECT CASE ( dig )
      CASE ( 64 )
         galmxCreateIntegerMatrix = mxCreateNumericArray( 2, 
     *     (/ m, n /), mxClassIDFromClassName('int64'), 0 ) 
      CASE ( 16 )
        galmxCreateIntegerMatrix = mxCreateNumericArray( 2, 
     *     (/ m, n /), mxClassIDFromClassName('int16'), 0 ) 
      CASE ( 8 )
        galmxCreateIntegerMatrix = mxCreateNumericArray( 2, 
     *      (/ m, n /), mxClassIDFromClassName('int8'),  0 ) 
      CASE default
        galmxCreateIntegerMatrix = mxCreateNumericArray( 2, 
     *      (/ m, n /), mxClassIDFromClassName('int32'), 0 ) 
      END SELECT

      RETURN
      END FUNCTION galmxCreateIntegerMatrix

      FUNCTION galmxCreateRealMatrix( m, n )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      mwPointer :: galmxCreateRealMatrix
      mwSize :: m, n

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Create an unpopulated REAL pointer matrix

!  Arguments

!  m - number of rows
!  n - number of columns

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 4th 2007

      REAL ( KIND = wp ) :: ddig = 1.0_wp

      SELECT CASE ( digits( ddig ) )
      CASE ( 53 )
        galmxCreateRealMatrix = mxCreateNumericArray( 2, (/ m, n /),
     *                      mxClassIDFromClassName('double'), 0 ) 
      CASE default
        galmxCreateRealMatrix = mxCreateNumericArray( 2, (/ m, n /),
     *                      mxClassIDFromClassName('single'), 0 ) 
      END SELECT

      RETURN
      END FUNCTION galmxCreateRealMatrix


      SUBROUTINE galmxCreateStructureComponent( struct, name, pr,
     *                                          ninform, finform )
      mwPointer ::struct
      CHARACTER ( len = * ) :: name
      mwPointer :: pr
      mwSize :: ninform
      CHARACTER ( LEN = * ), DIMENSION( ninform ) :: finform

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Create a named INTEGER component of a structure

!  Arguments

!  struct - structure
!  name - name of component
!  pr - pointer to the structure
!  ninform - number of sub-components
!  finform - names of sub-components

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 24th 2007

      mwPointer :: galmxCreateInteger

      pr = mxCreateStructMatrix( 1, 1, ninform, finform )
      CALL mxSetField( struct, 1, name, pr )

      RETURN
      END SUBROUTINE galmxCreateStructureComponent

      SUBROUTINE galmxCreateIntegerComponent( struct, name, pr )
      mwPointer ::struct
      CHARACTER ( len = * ) :: name
      mwPointer :: pr

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Create a named INTEGER component of a structure

!  Arguments

!  struct - structure
!  name - name of component
!  pr - pointer to the structure

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 24th 2007

      mwPointer :: galmxCreateInteger

      pr = galmxCreateInteger( ) 
      CALL mxSetField( struct, 1, name, pr )

      RETURN
      END SUBROUTINE galmxCreateIntegerComponent

      SUBROUTINE galmxCreateRealComponent( struct, name, pr )
      mwPointer ::struct
      CHARACTER ( len = * ) :: name
      mwPointer :: pr

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Create a named REAL component of a structure

!  Arguments

!  struct - structure
!  name - name of component
!  pr - pointer to the structure

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 24th 2007

      mwPointer :: galmxCreateReal
 
      pr = galmxCreateReal( ) 
      CALL mxSetField( struct, 1, name, pr )

      RETURN
      END SUBROUTINE galmxCreateRealComponent

      SUBROUTINE galmxCreateCharacterComponent( struct, name, pr )
      mwPointer ::struct
      CHARACTER ( len = * ) :: name
      mwPointer :: pr

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Create a named CHARACTER component of a structure

!  Arguments

!  struct - structure
!  name - name of component
!  pr - pointer to the structure

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 24th 2007

      INTEGER, PARAMETER :: len_blank = 80
      CHARACTER ( len = len_blank ) :: blank
      blank = REPEAT( 'w', len_blank )

      pr = mxCreateCharMatrixFromStrings( 1, blank )
      CALL mxSetField( struct, 1, name, pr )

      RETURN
      END SUBROUTINE galmxCreateCharacterComponent

      SUBROUTINE galmxCreateLogicalComponent( struct, name, pr )
      mwPointer ::struct
      CHARACTER ( len = * ) :: name
      mwPointer :: pr

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Create a named LOGICAL component of a structure

!  Arguments

!  struct - structure
!  name - name of component
!  pr - pointer to the structure

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 24th 2007

! **NB** This is a bodge since Mex doesn't appear to handle Fortran logicals

      mwPointer :: galmxCreateInteger

      pr = galmxCreateInteger( ) 
      CALL mxSetField( struct, 1, name, pr )

      RETURN
      END SUBROUTINE galmxCreateLogicalComponent

      SUBROUTINE galmxSetCharacterComponent( struct, name, value )
      mwPointer ::struct
      CHARACTER ( len = * ) :: name
      CHARACTER ( len = * ) :: value

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Set a named CHARACTER component of a structure to a value

!  Arguments

!  struct - structure
!  name - name of component
!  pr - pointer to the structure
!  value - character value to be assigned

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 24th 2007

      mwPointer :: y_pr
      mwPointer :: mxCreateString

      CALL mxSetField( struct, 1, name, mxCreateString( value ) )

      RETURN
      END SUBROUTINE galmxSetCharacterComponent
