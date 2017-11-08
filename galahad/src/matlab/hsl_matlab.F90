#include <fintrf.h>

!  THIS VERSION: GALAHAD 2.1 - 25/07/2007 AT 11:00 GMT.

!-*-*-*-*-*-*-*-*- H S L _ M A T L A B   M O D U L E -*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1. July 25th 2007

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

    MODULE HSL_MATLAB

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: MATLAB_copy_from_ptr,                                          &
                MATLAB_copy_to_ptr,                                            &
                MATLAB_copy_single_to_ptr,                                     &
                MATLAB_create_integer,                                         &
                MATLAB_create_real,                                            &
                MATLAB_create_substructure,                                    &
                MATLAB_create_integer_component,                               &
                MATLAB_create_real_component,                                  &
                MATLAB_create_char_component,                                  &
                MATLAB_create_logical_component,                               &
                MATLAB_get_value

!--------------------
!   P r e c i s i o n
!--------------------

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!---------------------------------
!   I n t e r f a c e  B l o c k s
!---------------------------------

      INTERFACE MATLAB_copy_from_ptr
        MODULE PROCEDURE hslmxCopyPtrToInteger,                                &
                         hslmxCopyPtrToReal,                                   &
                         hslmxCopyPtrToRealArray,                              &
                         hslmxCopyPtrToRealMatrix                            
      END INTERFACE MATLAB_copy_from_ptr

      INTERFACE MATLAB_copy_to_ptr
        MODULE PROCEDURE hslmxCopyIntegerToPtr,                                &
                         hslmxCopyIntegerArrayToPtr,                           &
                         hslmxCopyRealToPtr,                                   &
                         hslmxCopyRealArrayToPtr,                              &
                         hslmxCopyRealMatrixToPtr,                             &
                         hslmxCopyLogicalToPtr,                                &
                         hslmxCopyLogicalArrayToPtr,                           &
                         hslmxSetCharacterComponent
      END INTERFACE MATLAB_copy_to_ptr

      INTERFACE MATLAB_copy_single_to_ptr
        MODULE PROCEDURE hslmxCopySingleToPtr,                                 &
                         hslmxCopySingleArrayToPtr,                            &
                         hslmxCopySingleMatrixToPtr
      END INTERFACE MATLAB_copy_single_to_ptr

      INTERFACE MATLAB_create_integer
        MODULE PROCEDURE hslmxCreateInteger,                                   &
                         hslmxCreateIntegerArray,                              &
                         hslmxCreateIntegerMatrix
      END INTERFACE MATLAB_create_integer

      INTERFACE MATLAB_create_real
        MODULE PROCEDURE hslmxCreateReal,                                      &
                         hslmxCreateRealArray,                                 &
                         hslmxCreateRealMatrix
      END INTERFACE MATLAB_create_real

      INTERFACE MATLAB_create_integer_component
        MODULE PROCEDURE hslmxCreateIntegerComponent,                          &
                         hslmxCreateIntegerArrayComponent,                     &
                         hslmxCreateIntegerMatrixComponent
      END INTERFACE MATLAB_create_integer_component

      INTERFACE MATLAB_create_real_component
        MODULE PROCEDURE hslmxCreateRealComponent,                             &
                         hslmxCreateRealArrayComponent,                        &
                         hslmxCreateRealMatrixComponent
      END INTERFACE MATLAB_create_real_component

      INTERFACE MATLAB_get_value
        MODULE PROCEDURE hslmxGetInteger,                                      &
                         hslmxGetReal,                                         &
                         hslmxGetLogical,                                      &
                         hslmxGetCharacter
      END INTERFACE MATLAB_get_value

    CONTAINS

!  -*-*-*-  M A T L A B _ c r e a t e _ s u b s t r u c t u r e  -*-*-*-*-

      SUBROUTINE MATLAB_create_substructure( struct, name, pr,                 &
                                              ninform, finform )
      mwPointer ::struct, pr
      CHARACTER ( len = * ) :: name
      mwSize :: ninform
      CHARACTER ( LEN = * ), DIMENSION( ninform ) :: finform

!  -----------------------------------------------------

!  Create a named sub-structure of a structure

!  Arguments

!  struct - existing pointer to the structure
!  name - name of component of the structure
!  pr - pointer to the sub-tructure
!  ninform - number of components of the substructure
!  finform - names of components of the substructure

!  ------------------------------------------------------

      mwPointer :: mxCreateStructMatrix

      pr = mxCreateStructMatrix( 1, 1, ninform, finform )
      CALL mxSetField( struct, 1, name, pr )

      RETURN
      END SUBROUTINE MATLAB_create_substructure


!  -*-*-*- h s l m x  C r e a t e  I n t e g e r  C o m p o n e n t -*-*-*-*-

      SUBROUTINE hslmxCreateIntegerComponent( struct, name, pr )
      mwPointer ::struct
      CHARACTER ( len = * ) :: name
      mwPointer :: pr

!  -----------------------------------------------

!  Create a named INTEGER component of a structure

!  Arguments

!  struct - structure
!  name - name of component
!  pr - pointer to the structure

!  -----------------------------------------------

      pr = hslmxCreateInteger( ) 
      CALL mxSetField( struct, 1, name, pr )

      RETURN
      END SUBROUTINE hslmxCreateIntegerComponent


!  -*- h s l m x  C r e a t e  I n t e g e r  A r r a y  C o m p o n e n t -*-

      SUBROUTINE hslmxCreateIntegerArrayComponent( struct, name, n, pr )
      mwPointer ::struct
      CHARACTER ( len = * ) :: name
      INTEGER :: n
      mwPointer :: pr

!  -----------------------------------------------------

!  Create a named INTEGER array component of a structure

!  Arguments

!  struct - structure
!  name - name of component
!  n - dimension of array
!  pr - pointer to the structure

!  -----------------------------------------------------

      pr = hslmxCreateIntegerArray( n ) 
      CALL mxSetField( struct, 1, name, pr )

      RETURN
      END SUBROUTINE hslmxCreateIntegerArrayComponent

!  -*- h s l m x  C r e a t e  I n t e g e r  M a t r i x  C o m p o n e n t -*-

      SUBROUTINE hslmxCreateIntegerMatrixComponent( struct, name, m, n, pr )
      mwPointer ::struct
      CHARACTER ( len = * ) :: name
      INTEGER :: m, n
      mwPointer :: pr

!  -----------------------------------------------------

!  Create a named INTEGER matrix component of a structure

!  Arguments

!  struct - structure
!  name - name of component
!  m - row dimension of matrix
!  n - column dimension of matrix
!  pr - pointer to the structure

!  -----------------------------------------------------

      pr = hslmxCreateIntegerMatrix( m, n ) 
      CALL mxSetField( struct, 1, name, pr )

      RETURN
      END SUBROUTINE hslmxCreateIntegerMatrixComponent


!  -*-*-*-*-*- h s l m x  C r e a t e  R e a l  C o m p o n e n t -*-*-*-*-*-*-

      SUBROUTINE hslmxCreateRealComponent( struct, name, pr )
      mwPointer ::struct
      CHARACTER ( len = * ) :: name
      mwPointer :: pr

!  --------------------------------------------

!  Create a named REAL component of a structure

!  Arguments

!  struct - structure
!  name - name of component
!  pr - pointer to the structure

!  --------------------------------------------

      pr = hslmxCreateReal( ) 
      CALL mxSetField( struct, 1, name, pr )

      RETURN
      END SUBROUTINE hslmxCreateRealComponent


!  -*-*-*- h s l m x  C r e a t e  R e a l  A r r a y  C o m p o n e n t -*-*-

      SUBROUTINE hslmxCreateRealArrayComponent( struct, name, n, pr )
      mwPointer ::struct
      CHARACTER ( len = * ) :: name
      INTEGER :: n
      mwPointer :: pr

!  --------------------------------------------------

!  Create a named REAL array component of a structure

!  Arguments

!  struct - structure
!  name - name of component
!  n - dimension of array
!  pr - pointer to the structure

!  --------------------------------------------------

      pr = hslmxCreateRealArray( n ) 
      CALL mxSetField( struct, 1, name, pr )

      RETURN
      END SUBROUTINE hslmxCreateRealArrayComponent


!  -*-*-*- h s l m x  C r e a t e  R e a l  M a t r i x  C o m p o n e n t -*-*-

      SUBROUTINE hslmxCreateRealMatrixComponent( struct, name, m, n, pr )
      mwPointer ::struct
      CHARACTER ( len = * ) :: name
      INTEGER :: m, n
      mwPointer :: pr

!  --------------------------------------------------

!  Create a named REAL matrix component of a structure

!  Arguments

!  struct - structure
!  name - name of component
!  m - row dimension of matrix
!  n - column dimension of matrix
!  pr - pointer to the structure

!  --------------------------------------------------

      pr = hslmxCreateRealMatrix( m, n ) 
      CALL mxSetField( struct, 1, name, pr )

      RETURN
      END SUBROUTINE hslmxCreateRealMatrixComponent


!  -*-*-*-  M A T L A B _ c r e a t e _ c h a r _ c o m p o n e n t  -*-*-*-

      SUBROUTINE MATLAB_create_char_component( struct, name, pr )
      mwPointer ::struct
      CHARACTER ( len = * ) :: name
      mwPointer :: pr

!  -------------------------------------------------

!  Create a named CHARACTER component of a structure

!  Arguments

!  struct - structure
!  name - name of component
!  pr - pointer to the structure

!  -------------------------------------------------

      INTEGER, PARAMETER :: len_blank = 80
      CHARACTER ( len = len_blank ) :: blank
      mwPointer :: mxCreateCharMatrixFromStrings

      blank = REPEAT( ' ', len_blank )

      pr = mxCreateCharMatrixFromStrings( 1, blank )
      CALL mxSetField( struct, 1, name, pr )

      RETURN
      END SUBROUTINE MATLAB_create_char_component


!  -*-*-  M A T L A B _ c r e a t e _ l o g i c a l _ c o m p o n e n t  -*-*-

      SUBROUTINE MATLAB_create_logical_component( struct, name, pr )
      mwPointer ::struct
      CHARACTER ( len = * ) :: name
      mwPointer :: pr

!  -----------------------------------------------------

!  Create a named LOGICAL component of a structure

!  Arguments

!  struct - structure
!  name - name of component
!  pr - pointer to the structure

! ** - NB - ** This is a bodge since Mex doesn't appear 
!              to handle Fortran logicals ** - NB - **

!  ----------------------------------------------------

      pr = hslmxCreateInteger( ) 
      CALL mxSetField( struct, 1, name, pr )

      RETURN
      END SUBROUTINE MATLAB_create_logical_component


!  -*-*-*-*-*-*-*-*-*-*- h s l m x  G e t  I n t e g e r  -*-*-*-*-*-*-*-*-*-*-

      SUBROUTINE hslmxGetInteger( ps, name, pc, value )
      mwPointer :: ps, pc
      CHARACTER ( LEN = * ) :: name
      INTEGER :: value

!  ---------------------------------------------------------

!  Obtain an integer value from the component of a structure

!  Arguments

!  ps - given pointer to the structure
!  name - name of component
!  pc - pointer to the component
!  value - value of the component

!  ---------------------------------------------------------

      mwPointer :: mxGetField
      REAL ( KIND = wp ) :: mxGetScalar 

      pc = mxGetField( ps, 1, name )
      value = INT( mxGetScalar( pc ) )

      RETURN
      END SUBROUTINE hslmxGetInteger


!  -*-*-*-*-*-*-*-*-*-*- h s l m x  G e t  R e a l  -*-*-*-*-*-*-*-*-*-*-*-

      SUBROUTINE hslmxGetReal( ps, name, pc, value )
      mwPointer :: ps, pc
      CHARACTER ( LEN = * ) :: name
      REAL ( KIND = wp ) :: value

!  -----------------------------------------------------

!  Obtain a real value from the component of a structure

!  Arguments

!  ps - given pointer to the structure
!  name - name of component
!  pc - pointer to the component
!  value - value of the component

!  -----------------------------------------------------

      mwPointer :: mxGetField
      REAL ( KIND = wp ) :: mxGetScalar 

      pc = mxGetField( ps, 1, name )
      value = mxGetScalar( pc )

      RETURN
      END SUBROUTINE hslmxGetReal


!  -*-*-*-*-*-*-*-*-*-*- h s l m x  G e t  L o g i c a l  -*-*-*-*-*-*-*-*-*-*-

      SUBROUTINE hslmxGetLogical( ps, name, pc, value )
      mwPointer :: ps, pc
      CHARACTER ( LEN = * ) :: name
      LOGICAL :: value

!  ---------------------------------------------------------

!  Obtain a logical value from the component of a structure

!  Arguments

!  ps - given pointer to the structure
!  name - name of component
!  pc - pointer to the component
!  value - value of the component

! ** - NB - ** This is a bodge since Mex doesn't appear 
!              to handle Fortran logicals ** - NB - **

!  ---------------------------------------------------------

      mwPointer :: mxGetField
      REAL ( KIND = wp ) :: mxGetScalar 

      pc = mxGetField( ps, 1, name )
      IF ( INT( mxGetScalar( pc ) ) == 1 ) THEN
        value = .TRUE.
      ELSE
        value = .FALSE.
      END IF

      RETURN
      END SUBROUTINE hslmxGetLogical


!  -*-*-*-*-*-*-*-*-*-*- h s l m x  G e t  C h a r a c t e r  -*-*-*-*-*-*-*-*-*-

      SUBROUTINE hslmxGetCharacter( ps, name, pc, value, len )
      mwPointer :: ps, pc
      CHARACTER ( LEN = * ) :: name
      CHARACTER ( LEN = * ) :: value
      INTEGER :: len

!  -----------------------------------------------------------------

!  Obtain a character string value from the component of a structure

!  Arguments

!  ps - given pointer to the structure
!  name - name of component
!  pc - pointer to the component
!  value - value of the character string
!  len - maximum length of the string

!  ------------------------------------------------------------------

      mwSize :: i, mxGetString
      mwPointer :: mxGetField

      pc = mxGetField( ps, 1, name )
      i = mxGetString( pc, value, len )

      RETURN
      END SUBROUTINE hslmxGetCharacter


!  -*-*-*-*-*-*- h s l m x  C o p y  P t r  T o  I n t e g e r  -*-*-*-*-*-*- 

      SUBROUTINE hslmxCopyPtrToInteger( px, Y, n )
      mwPointer :: px
      mwSize :: n
      MWSIZE, DIMENSION( n ) :: Y

!  --------------------------------------------------------------

!  Copy INTEGER values from Matlab pointer array to Fortran array

!  Arguments

!  px  - Pointer to ir or jc array
!  Y   - Integer Fortran array
!  n   - number of elements to copy

!  --------------------------------------------------------------

      INTEGER :: dig

      SELECT CASE ( digits( dig ) + 1 )
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
      END SUBROUTINE hslmxCopyPtrToInteger


!  -*-*-*-*-*-*- h s l m x  C o p y  P t r  T o  R e a l  -*-*-*-*-*-*- 

      SUBROUTINE hslmxCopyPtrToReal( px, Y )
      mwPointer :: px
      REAL ( KIND = wp ) :: Y

!  -----------------------------------------------------------

!  Copy REAL values from Matlab pointer array to Fortran array

!  Arguments

!  px  - Pointer to ir or jc array
!  Y   - Real Fortran array
!  n   - number of elements to copy

!  -----------------------------------------------------------

      REAL ( KIND = wp ) :: ddig = 1.0_wp

      SELECT CASE ( digits( ddig ) )
      CASE ( 53 )
        CALL mxCopyPtrToReal8( px, Y, 1 )
      CASE default
        CALL mxCopyPtrToReal4( px, Y, 1 )
      END SELECT

      RETURN
      END SUBROUTINE hslmxCopyPtrToReal


!  -*-*-*-*-*- h s l m x  C o p y  P t r  T o  R e a l  A r r a y   -*-*-*-*-*- 

      SUBROUTINE hslmxCopyPtrToRealArray( px, Y, n )
      mwPointer :: px
      mwSize :: n
      REAL ( KIND = wp ), DIMENSION( n ) :: Y

!  -----------------------------------------------------------

!  Copy REAL values from Matlab pointer array to Fortran array

!  Arguments

!  px  - Pointer to ir or jc array
!  Y   - Real Fortran array
!  n   - number of elements to copy

!  -----------------------------------------------------------

      REAL ( KIND = wp ) :: ddig = 1.0_wp

      SELECT CASE ( digits( ddig ) )
      CASE ( 53 )
        CALL mxCopyPtrToReal8( px, Y, n )
      CASE default
        CALL mxCopyPtrToReal4( px, Y, n )
      END SELECT

      RETURN
      END SUBROUTINE hslmxCopyPtrToRealArray


!  -*-*-*-*-*- h s l m x  C o p y  P t r  T o  R e a l  M a t r i x    -*-*-*-*- 

      SUBROUTINE hslmxCopyPtrToRealMatrix( px, Y, m, n )
      mwPointer :: px
      mwSize :: m, n
      REAL ( KIND = wp ), DIMENSION( m, n ) :: Y

!  -----------------------------------------------------------

!  Copy REAL values from Matlab pointer array to Fortran array

!  Arguments

!  px  - Pointer to ir or jc array
!  Y   - Real Fortran array
!  m   - row dimension
!  n   - column dimension
!
!  -----------------------------------------------------------

      REAL ( KIND = wp ) :: ddig = 1.0_wp

      SELECT CASE ( digits( ddig ) )
      CASE ( 53 )
        CALL mxCopyPtrToReal8( px, Y, m * n )
      CASE default
        CALL mxCopyPtrToReal4( px, Y, m * n )
      END SELECT

      RETURN
      END SUBROUTINE hslmxCopyPtrToRealMatrix


!  -*-*-*-*-*-  h s l m x  C o p y  I n t e g e r  T o  P t r  -*-*-*-*-*-*-*

      SUBROUTINE hslmxCopyIntegerToPtr( Y, px )
      mwPointer :: px
      MWSIZE :: Y

!  ---------------------------------------------------------

!  Copy INTEGER values from Fortran scalar to Matlab pointer

!  Arguments

!  Y   - Integer Fortran scalar
!  px  - Pointer to ir or jc array
!  n   - number of elements to copy

!  ---------------------------------------------------------

      INTEGER :: dig

      SELECT CASE ( digits( dig ) + 1 )
      CASE ( 64 )
        CALL mxCopyInteger8ToPtr( Y, px, 1 )
      CASE ( 16 )
        CALL mxCopyInteger2ToPtr( Y, px, 1 )
      CASE ( 8 )
        CALL mxCopyInteger1ToPtr( Y, px, 1 )
      CASE default
        CALL mxCopyInteger4ToPtr( Y, px, 1 )
      END SELECT

      RETURN
      END SUBROUTINE hslmxCopyIntegerToPtr


!  -*-*-*-  h s l m x  C o p y  I n t e g e r  A r r a y  T o  P t r  -*-*-*-*-

      SUBROUTINE hslmxCopyIntegerArrayToPtr( Y, px, n )
      mwPointer :: px
      mwSize :: n
      MWSIZE, DIMENSION( n ) :: Y

!  --------------------------------------------------------------

!  Copy INTEGER values from Fortran array to Matlab pointer array

!  Arguments

!  Y   - Integer Fortran array
!  px  - Pointer to ir or jc array
!  n   - number of elements to copy

!  --------------------------------------------------------------

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
      END SUBROUTINE hslmxCopyIntegerArrayToPtr


!  -*-*-*-*-  h s l m x  C o p y  R e a l  T o  P t r  -*-*-*-*-*-

      SUBROUTINE hslmxCopyRealToPtr( Y, px )
      mwPointer :: px
      REAL ( KIND = wp ) :: Y

!  -----------------------------------------------------

!  Copy REAL value from Fortran scalar to Matlab pointer

!  Arguments

!  Y   - Real Fortran scalar
!  px  - Pointer to ir or jc array

!  -----------------------------------------------------

      REAL ( KIND = wp ) :: ddig = 1.0_wp

      SELECT CASE ( digits( ddig ) )
      CASE ( 53 )
        CALL mxCopyReal8ToPtr( Y, px, 1 )
      CASE default
        CALL mxCopyReal4ToPtr( Y, px, 1 )
      END SELECT

      RETURN
      END SUBROUTINE hslmxCopyRealToPtr


!  -*-*-*-  h s l m x  C o p y  R e a l  A r r a y  T o  P t r  -*-*-*-*-

      SUBROUTINE hslmxCopyRealArrayToPtr( Y, px, n )
      mwPointer :: px
      mwSize :: n
      REAL ( KIND = wp ), DIMENSION( n ) :: Y

!  -----------------------------------------------------------

!  Copy REAL values from Fortran array to Matlab pointer array

!  Arguments

!  Y   - Real 1-D Fortran array
!  px  - Pointer to ir or jc array
!  n   - number of elements to copy

!  -----------------------------------------------------------

      REAL ( KIND = wp ) :: ddig = 1.0_wp

      SELECT CASE ( digits( ddig ) )
      CASE ( 53 )
        CALL mxCopyReal8ToPtr( Y, px, n )
      CASE default
        CALL mxCopyReal4ToPtr( Y, px, n )
      END SELECT

      RETURN
      END SUBROUTINE hslmxCopyRealArrayToPtr


!  -*-*-*-  h s l m x  C o p y  R e a l  M a t r i x  T o  P t r  -*-*-*-*-

      SUBROUTINE hslmxCopyRealMatrixToPtr( Y, px, m, n )
      mwPointer :: px
      mwSize :: m, n
      REAL ( KIND = wp ), DIMENSION( m, n ) :: Y

!  -----------------------------------------------------------

!  Copy REAL values from Fortran array to Matlab pointer array

!  Arguments

!  Y   - Real 2-D Fortran array
!  px  - Pointer to ir or jc array
!  m   - row dimension
!  n   - column dimension

!  -----------------------------------------------------------

      REAL ( KIND = wp ) :: ddig = 1.0_wp

      SELECT CASE ( digits( ddig ) )
      CASE ( 53 )
        CALL mxCopyReal8ToPtr( Y, px, m * n )
      CASE default
        CALL mxCopyReal4ToPtr( Y, px, m * n )
      END SELECT

      RETURN
      END SUBROUTINE hslmxCopyRealMatrixToPtr


!  -*-*-*-*-  h s l m x  C o p y  D o u b l e  T o  P t r  -*-*-*-*-*-

      SUBROUTINE hslmxCopySingleToPtr( Y, px )
      mwPointer :: px
      REAL :: Y

!  -----------------------------------------------------------------

!  Copy SINGLE PRECISION value from Fortran scalar to Matlab pointer

!  Arguments

!  Y   - Real Fortran scalar
!  px  - Pointer to ir or jc array

!  -----------------------------------------------------------------

      CALL mxCopyReal4ToPtr( Y, px, 1 )

      RETURN
      END SUBROUTINE hslmxCopySingleToPtr


!  -*-*-*-  h s l m x  C o p y  D o u b l e  A r r a y  T o  P t r  -*-*-*-*-

      SUBROUTINE hslmxCopySingleArrayToPtr( Y, px, n )
      mwPointer :: px
      mwSize :: n
      REAL, DIMENSION( n ) :: Y

!  ------------------------------------------------------------------------

!  Copy SINGLE PRECISION values from Fortran array to Matlab pointer array

!  Arguments

!  Y   - Single precision 1-D Fortran array
!  px  - Pointer to ir or jc array
!  n   - number of elements to copy

!  ------------------------------------------------------------------------

      CALL mxCopyReal4ToPtr( Y, px, n )

      RETURN
      END SUBROUTINE hslmxCopySingleArrayToPtr


!  -*-*-*-  h s l m x  C o p y  D o u b l e  M a t r i x  T o  P t r  -*-*-*-*-

      SUBROUTINE hslmxCopySingleMatrixToPtr( Y, px, m, n )
      mwPointer :: px
      mwSize :: m, n
      REAL, DIMENSION( m, n ) :: Y

!  ------------------------------------------------------------------------

!  Copy SINGLE PRECISION values from Fortran array to Matlab pointer array

!  Arguments

!  Y   - Single precision 2-D Fortran array
!  px  - Pointer to ir or jc array
!  m   - row dimension
!  n   - column dimension

!  ------------------------------------------------------------------------

      CALL mxCopyReal4ToPtr( Y, px, m * n )

      RETURN
      END SUBROUTINE hslmxCopySingleMatrixToPtr


!  -*-*-*-  h s l m x  C o p y  L o g i c a l  T o  P t r  -*-*-*-*-

      SUBROUTINE hslmxCopyLogicalToPtr( Y, px )
      mwPointer :: px
      LOGICAL :: Y

!  --------------------------------------------------------

!  Copy LOGICAL value from Fortran scalar to Matlab pointer

!  Arguments

!  Y   - Logical Fortran scalar
!  px  - Pointer to ir or jc array

!  --------------------------------------------------------

      INTEGER :: ly

      IF ( Y ) THEN
        ly = 1
      ELSE
        ly = 0
      END IF
      CALL hslmxCopyIntegerToPtr( ly, px )

      RETURN
      END SUBROUTINE hslmxCopyLogicalToPtr


!  -*-*-*-  h s l m x  C o p y  L o g i c a l  A r r a y  T o  P t r  -*-*-*-*-

      SUBROUTINE hslmxCopyLogicalArrayToPtr( Y, px, n )
      mwPointer :: px
      mwSize :: n
      LOGICAL, DIMENSION( n ) :: Y

!  --------------------------------------------------------------

!  Copy LOGICAL values from Fortran array to Matlab pointer array

!  Arguments

!  Y   - Logical Fortran array
!  px  - Pointer to ir or jc array
!  n   - number of elements to copy

!  --------------------------------------------------------------

      INTEGER :: i
      INTEGER, DIMENSION( n ) :: LY

      DO i = 1, n
        IF ( Y( i ) ) THEN
          LY( i ) = 1
        ELSE
          LY( i ) = 0
        END IF
      END DO
      CALL hslmxCopyIntegerArrayToPtr( LY, px, n )

      RETURN
      END SUBROUTINE hslmxCopyLogicalArrayToPtr


!  -*-*-*-  h s l m x  S e t  C h a r a c t e r  C o m p o n e n t  -*-*-*-*-

      SUBROUTINE hslmxSetCharacterComponent( struct, name, value )
      mwPointer ::struct
      CHARACTER ( len = * ) :: name
      CHARACTER ( len = * ) :: value

!  ---------------------------------------------------------

!  Set a named CHARACTER component of a structure to a value

!  Arguments

!  struct - structure
!  name - name of component
!  value - character value to be assigned

!  ---------------------------------------------------------

      mwPointer :: mxCreateString

      CALL mxSetField( struct, 1, name, mxCreateString( value ) )

      RETURN
      END SUBROUTINE hslmxSetCharacterComponent


!  -*-*-*-  h s l m x  C r e a t e  I n t e g e r  -*-*-*-*-

      FUNCTION hslmxCreateInteger( )
      mwPointer :: hslmxCreateInteger

!  -------------------------------------

!  Create an unpopulated INTEGER pointer

!  No arguments

!  -------------------------------------

      INTEGER :: dig
      INTEGER * 4 :: mxClassIDFromClassName
      mwPointer :: mxCreateNumericMatrix

      SELECT CASE ( digits( dig ) + 1 )
      CASE ( 64 )
        hslmxCreateInteger = mxCreateNumericMatrix( 1, 1,                      &
                               mxClassIDFromClassName('int64'), 0 ) 
      CASE ( 16 )
        hslmxCreateInteger = mxCreateNumericMatrix( 1, 1,                      &
                               mxClassIDFromClassName('int16'), 0 ) 
      CASE ( 8 )
        hslmxCreateInteger = mxCreateNumericMatrix( 1, 1,                      &
                               mxClassIDFromClassName('int8'),  0 ) 
      CASE default
        hslmxCreateInteger = mxCreateNumericMatrix( 1, 1,                      &
                               mxClassIDFromClassName('int32'), 0 ) 
      END SELECT

      RETURN
      END FUNCTION hslmxCreateInteger


!  -*-*-*-  h s l m x  C r e a t e  I n t e g e r  A r r a y -*-*-*-*-

      FUNCTION hslmxCreateIntegerArray( n )
      mwPointer :: hslmxCreateIntegerArray
      mwSize :: n

!  --------------------------------------------

!  Create an unpopulated INTEGER pointer array

!  Arguments

!  n - number of entries

!  --------------------------------------------

      INTEGER :: dig
      INTEGER * 4 :: mxClassIDFromClassName
      mwPointer :: mxCreateNumericArray

      SELECT CASE ( digits( dig ) + 1 )
      CASE ( 64 )
         hslmxCreateIntegerArray = mxCreateNumericArray( 2,                    &
           (/ n, 1 /), mxClassIDFromClassName('int64'), 0 ) 
      CASE ( 16 )
        hslmxCreateIntegerArray = mxCreateNumericArray( 2,                     & 
           (/ n, 1 /), mxClassIDFromClassName('int16'), 0 )
      CASE ( 8 )
        hslmxCreateIntegerArray = mxCreateNumericArray( 2,                     &
            (/ n, 1 /), mxClassIDFromClassName('int8'),  0 )
      CASE default
        hslmxCreateIntegerArray = mxCreateNumericArray( 2,                     &
            (/ n, 1 /), mxClassIDFromClassName('int32'), 0 )
      END SELECT

      RETURN
      END FUNCTION hslmxCreateIntegerArray

!  -*-*-*-  h s l m x  C r e a t e  I n t e g e r  M a t r i x -*-*-*-*-

      FUNCTION hslmxCreateIntegerMatrix( m, n )
      mwPointer :: hslmxCreateIntegerMatrix
      mwSize :: m, n

!  --------------------------------------------

!  Create an unpopulated INTEGER pointer array

!  Arguments

!  m - number of rows
!  n - number of columns

!  --------------------------------------------

      INTEGER :: dig
      INTEGER * 4 :: mxClassIDFromClassName
      mwPointer :: mxCreateNumericArray

      SELECT CASE ( digits( dig ) + 1 )
      CASE ( 64 )
         hslmxCreateIntegerMatrix = mxCreateNumericArray( 2,                    &
           (/ m, n /), mxClassIDFromClassName('int64'), 0 ) 
      CASE ( 16 )
        hslmxCreateIntegerMatrix = mxCreateNumericArray( 2,                     &
           (/ m, n /), mxClassIDFromClassName('int16'), 0 )
      CASE ( 8 )
        hslmxCreateIntegerMatrix = mxCreateNumericArray( 2,                     &
            (/ m, n /), mxClassIDFromClassName('int8'),  0 )
      CASE default
        hslmxCreateIntegerMatrix = mxCreateNumericArray( 2,                     &
            (/ m, n /), mxClassIDFromClassName('int32'), 0 )
      END SELECT

      RETURN
      END FUNCTION hslmxCreateIntegerMatrix


!  -*-*-*-  h s l m x  C r e a t e  R e a l   -*-*-*-*-

      FUNCTION hslmxCreateReal( )
      mwPointer :: hslmxCreateReal

!  ----------------------------------

!  Create an unpopulated REAL pointer

!  No arguments

!  ----------------------------------

      REAL ( KIND = wp ) :: ddig = 1.0_wp
      INTEGER * 4 :: mxClassIDFromClassName
      mwPointer :: mxCreateNumericMatrix

      SELECT CASE ( digits( ddig ) )
      CASE ( 53 )
        hslmxCreateReal = mxCreateNumericMatrix( 1, 1,                         &
                            mxClassIDFromClassName('double'), 0 ) 
      CASE default
        hslmxCreateReal = mxCreateNumericMatrix( 1, 1,                         &
                            mxClassIDFromClassName('single'), 0 ) 
      END SELECT

      RETURN
      END FUNCTION hslmxCreateReal


!  -*-*-*-  h s l m x  C r e a t e  R e a l  A r r a y   -*-*-*-*-

      FUNCTION hslmxCreateRealArray( n )
      mwPointer :: hslmxCreateRealArray
      mwSize :: n

!  -----------------------------------------

!  Create an unpopulated REAL pointer array

!  Arguments

!  n - number of entries

!  -----------------------------------------

      REAL ( KIND = wp ) :: ddig = 1.0_wp
      INTEGER * 4 :: mxClassIDFromClassName
      mwPointer :: mxCreateNumericArray

      SELECT CASE ( digits( ddig ) )
      CASE ( 53 )
        hslmxCreateRealArray = mxCreateNumericArray( 2, (/ n, 1 /),             &
                            mxClassIDFromClassName('double'), 0 ) 
      CASE default
        hslmxCreateRealArray = mxCreateNumericArray( 2, (/ n, 1 /),            &
                            mxClassIDFromClassName('single'), 0 ) 
      END SELECT

      RETURN
      END FUNCTION hslmxCreateRealArray


!  -*-*-*-  h s l m x  C r e a t e  R e a l  M a t r i x  -*-*-*-*-

      FUNCTION hslmxCreateRealMatrix( m, n )
      mwPointer :: hslmxCreateRealMatrix
      mwSize :: m, n

!  -----------------------------------------

!  Create an unpopulated REAL pointer matrix

!  Arguments

!  m - number of rows
!  n - number of columns

!  -----------------------------------------

      REAL ( KIND = wp ) :: ddig = 1.0_wp
      INTEGER * 4 :: mxClassIDFromClassName
      mwPointer :: mxCreateNumericArray

      SELECT CASE ( digits( ddig ) )
      CASE ( 53 )
        hslmxCreateRealMatrix = mxCreateNumericArray( 2, (/ m, n /),           &
                            mxClassIDFromClassName('double'), 0 ) 
      CASE default
        hslmxCreateRealMatrix = mxCreateNumericArray( 2, (/ m, n /),           &
                            mxClassIDFromClassName('single'), 0 ) 
      END SELECT

      RETURN
      END FUNCTION hslmxCreateRealMatrix

!-*-*-*-*-*-*- E N D  o f  H S L _ M A T L A B   M O D U L E -*-*-*-*-*-

    END MODULE HSL_MATLAB
