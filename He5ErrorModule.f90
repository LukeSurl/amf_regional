! $Id: He5ErrorModule.f90,v 1.1 2006/01/18 20:10:36 bmy Exp $
MODULE He5ErrorModule

  !========================================================================
  ! Module He5ErrorModule contains error check routines for the
  ! Fortran code that reads HDF5 and HDF-EOS5 data from disk.
  ! (bmy, 1/13/06)
  !
  ! Module Methods:
  ! ----------------------------------------------------------------------
  ! (1 ) He5AllocErr   : Prints an error message for allocating arrays
  ! (2 ) He5ErrMsg     : Prints an error message and halts execution
  ! (3 ) He5Msg        : Prints a message and flushes buffer
  ! (4 ) He5CheckValue : Checks a value for NaN or Infinity condition
  ! (5 ) ItIsNan       : Checks for NaN
  ! (6 ) ItIsFinite    : Checks for Infinity
  !
  ! NOTES:
  !========================================================================

  !-------------------------------
  ! PRIVATE / PUBLIC DECLARATIONS
  !-------------------------------

  ! Make everything PRIVATE ...
  PRIVATE

  ! ... except these routines
  PUBLIC :: He5AllocErr
  PUBLIC :: He5ErrMsg
  PUBLIC :: He5Msg
  PUBLIC :: He5CheckValue

  !-------------------------------
  ! MODULE ROUTINES
  !-------------------------------

CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE He5AllocErr( arrayName )

    !======================================================================
    ! Subroutine He5AllocErr halts program execution upon an allocation
    ! error. (bmy, 1/13/06)
    ! 
    ! Arguments as Input:
    ! ------------------------------------------------------------------
    ! (1) arrayName (CHARACTER) : Name of array with allocation error
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: arrayName

    !--------------------------
    ! He5AllocErr begins here!
    !--------------------------

    ! Print info
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    WRITE( 6, 100   ) TRIM( arrayName )
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    CALL FLUSH(6)

    ! Format string
100 FORMAT( 'Allocation error for array: ', a )

    ! Halt execution
    CALL EXIT(1)
      
  END SUBROUTINE He5AllocErr

!------------------------------------------------------------------------------

  SUBROUTINE He5ErrMsg( msg, loc )

    !======================================================================
    ! Subroutine He5ErrMsg halts program execution w/ an error message.
    ! (bmy, 1/13/06)
    ! 
    ! Arguments as Input:
    ! ------------------------------------------------------------------
    ! (1) msg (CHARACTER) : Error message to display
    ! (2) loc (CHARACTER) : Location where the error occurred
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: msg
    CHARACTER(LEN=*), INTENT(IN) :: loc

    !------------------------
    ! He5ErrMsg begins here!
    !------------------------

    ! Print info
    WRITE( 6, '(a)'               ) REPEAT( '=', 79 )
    WRITE( 6, '(a)'               ) TRIM( msg )
    WRITE( 6, '(''Stop in '', a)' ) TRIM( loc )
    WRITE( 6, '(a)'               ) REPEAT( '=', 79 )
    CALL FLUSH(6)

    ! Halt execution
    CALL EXIT(1)
      
  END SUBROUTINE He5ErrMsg

!------------------------------------------------------------------------------

  SUBROUTINE He5Msg( str )

    !======================================================================
    ! Subroutine He5Msg prints a string and flushes the output buffer.
    ! (bmy, 1/13/06)
    ! 
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) str (CHARACTER) : Message to display
    !======================================================================

    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: str

    !---------------------
    ! He5Msg begins here!
    !---------------------

    ! Print message
    WRITE( 6, '(a)' ) TRIM( str )
    CALL flush( 6 )
      
  END SUBROUTINE He5Msg

!-----------------------------------------------------------------------------

  SUBROUTINE He5CheckValue( value, name, loc )

    !======================================================================
    ! Subroutine He5CheckValue tests a value for IEEE NaN or Infinity.
    ! (bmy, 1/13/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) value (REAL*4   ) : value to be tested
    ! (2 ) name  (CHARACTER) : name of the variable 
    ! (3 ) loc   (INTEGER  ) : Grid box location (/i,j,l,t/)
    !======================================================================

    ! Arguments
    REAL*4,           INTENT(IN) :: value
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER,          INTENT(IN) :: loc(4)

    ! If VALUE is NaN, stop w/ error message
    IF ( itIsNaN( value ) ) THEN
       !$OMP CRITICAL
       WRITE( 6, 100 ) TRIM( name ), loc
 100   FORMAT( a, ' is NaN at grid box: ', 4i4, '!' )
       STOP
       !$OMP END CRITICAL
    ENDIF

    ! If VALUE is +/- Infinity, stop w/ error message
    IF ( .not. itIsFinite( value ) ) THEN
       !$OMP CRITICAL
       WRITE( 6, 110 ) TRIM( name ), loc
 110   FORMAT( a, ' is +/- Infinity at grid box: ', 4i4, '!' )
       STOP
       !$OMP END CRITICAL
    ENDIF

  END SUBROUTINE He5CheckValue

!-----------------------------------------------------------------------------

  FUNCTION itIsNan( value ) RESULT( itIsANaN )

    !===================================================================
    ! Subroutine itIsNaN tests a value for IEEE NaN on either SGI
    ! or Linux platforms (bmy, 1/13/06)
    !
    ! Arguments as Input:
    ! ------------------------------------------------------------------
    ! (1 ) value (REAL*4) : value to be tested
    !===================================================================

#include "define.h"

    ! Argument
    REAL*4, INTENT(IN) :: value
    LOGICAL            :: itIsANaN

    !----------------------
    ! itisNan begins here!
    !----------------------

#if defined( SGI32 ) || defined( SGI64 )

    ! Use SGI intrinsic function
    itIsANaN = IEEE_IS_NAN( value )   

#elif defined( ALTIX ) || defined( PC )

    ! Declare IS_NAN as an external function
    INTEGER, EXTERNAL  :: is_nan

    ! For LINUX or INTEL_FC compilers, use C routine "is_nan" to test if 
    ! VALUE is NaN.   VALUE must be cast to DBLE since "is_nan" only
    ! takes doubles.
    itIsANan = ( IS_NAN( DBLE( value ) ) /= 0 )

#endif

  END FUNCTION itIsNan

!-----------------------------------------------------------------------------

  FUNCTION itIsFinite( value ) RESULT( itIsAFinite )

    !===================================================================
    ! Subroutine itIsFinite tests a value for IEEE Finite on either 
    ! SGI or Linux platforms (bmy, 1/13/06)
    !
    ! Arguments as Input:
    ! ------------------------------------------------------------------
    ! (1 ) value (REAL*4) : value to be tested
    !===================================================================

#include "define.h"

    ! Arguments
    REAL*4, INTENT(IN) :: value
    LOGICAL            :: itIsAFinite

    !-------------------------
    ! itisFinite begins here!
    !-------------------------

#if defined( SGI32 ) || defined( SGI64 )

    ! Use SGI intrinsic function
    itIsAFinite = IEEE_FINITE( value )   

#elif defined( ALTIX ) || defined( PC )  

    ! Declare IS_NAN as an external function
    INTEGER, EXTERNAL  :: is_finite

    ! For LINUX or INTEL_FC compilers, use C routine "is_finite" to test if 
    ! VALUE is finite.   VALUE must be cast to DBLE since "is_finite" only
    ! takes doubles.
    itIsAFinite = ( IS_FINITE( DBLE( value ) ) /= 0 )

#endif

  END FUNCTION itIsFinite

!------------------------------------------------------------------------------

END MODULE He5ErrorModule
