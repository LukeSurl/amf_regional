C $Id: bpch2_mod.f,v 4.21 2001/10/29 15:48:34 bmy v421 $
      MODULE BPCH2_MOD
!
!******************************************************************************
!  Module BPCH2_MOD contains the routines used to read data from and write
!  data to binary punch file format (v. 2.0). (bmy, 6/28/00, 10/9/01)
!
!  Module Procedures:
!  ============================================================================
!  (1 ) BPCH2_HDR     : writes "top-of-file" header to a binary punch file
!  (2 ) BPCH2         : writes a data block to a binary punch file
!  (3 ) READ_BPCH2    : reads a data block from a binary punch file
!  (4 ) GET_MODELNAME : returns MODELNAME for the given met field data set
!  (5 ) GET_NAME_EXT  : returns file extension string for model name
!  (6 ) GET_RES_EXT   : returns file extension string for model resolution
!  (7 ) GET_TAU0_2A   : computes TAU0 from MONTH and YEAR
!  (8 ) GET_TAU0_6A   : computes TAU0 from MONTH, DAY, YEAR (, HOUR, MIN, SEC)
!  (9 ) JULDAY        : converts M/D/Y into astronomical Julian date
!  (10) MINT          : modified integer function, used by JULDAY
!  (11) CALDATE       : converts astronomical Julian date to YYYYMMDD, HHMMSS
!
!  Module Interfaces
!  ============================================================================
!  (1 ) GET_TAU0      : Overloads GET_TAU0_2A and GET_TAU0_6A
!
!  GEOS-CHEM modules referenced by bpch2_mod.f
!  ============================================================================
!  none
!
!  NOTES:
!  (1 ) Added routine GET_TAU0 (bmy, 7/20/00)
!  (2 ) Added years 1985-2001 for routine GET_TAU0 (bmy, 8/1/00)
!  (3 ) Use IOS /= 0 criterion to also check for EOF (bmy, 9/12/00)
!  (4 ) Removed obsolete code in "read_bpch2.f" (bmy, 12/18/00)
!  (5 ) Correct error for 1991 TAU values in GET_TAU0 (bnd, bmy, 1/4/01)
!  (6 ) BPCH2_MOD is now independent of any GEOS-CHEM size parameters.
!        (bmy, 4/18/01)
!  (7 ) Now have 2 versions of "GET_TAU0" overloaded by an interface.  The
!        original version takes 2 arguments (MONTH, YEAR).  The new version
!        takes 3 arguments (MONTH, DAY, YEAR). (bmy, 8/22/01)
!  (8 ) Updated comments (bmy, 9/4/01)
!  (9 ) Renamed GET_TAU0_3A to GET_TAU0_6A, and updated the GET_TAU0 
!        interface.  Also updated comments (bmy, 9/26/01)
!  (10) Now use special model name for GEOS-3 w/ 30 layers (bmy, 10/9/01)
!******************************************************************************
!
! Read LLPAR for GET_MODEL_NAME funtion
      USE CMN_SIZE_MOD,   ONLY: LLPAR

      IMPLICIT NONE

      ! Interface for current or new version of GET_TAU0
      INTERFACE GET_TAU0
         MODULE PROCEDURE GET_TAU0_2A, GET_TAU0_6A
      END INTERFACE

      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE BPCH2_HDR ( IUNIT, TITLE )
!
!******************************************************************************
!  Subroutine BPCH2_HDR writes a header at the top of the binary
!  punch file, version 2.0 (bmy, 5/27/99, 9/12/00).
!
!  Arguments as input:
!  ============================================================================
!  (1) IUNIT : INTEGER - logical unit number of binary punch file
!  (2) TITLE : CHAR*80 - description of data contained in binary punch file
!
!  NOTES:
!  (1 ) Added this routine to "bpch_mod.f" (bmy, 6/28/00)
!  (2 ) Use IOS /= 0 criterion to also check for EOF condition (bmy, 9/12/00)
!******************************************************************************
!
      ! Arguments
      INTEGER,           INTENT(IN) :: IUNIT
      CHARACTER(LEN=80), INTENT(IN) :: TITLE

      ! Local variable
      INTEGER                       :: IOS
      CHARACTER(LEN=40)             :: FTI = 'CTM bin 02'

      !=================================================================
      ! BPCH2_HDR begins here!
      !
      ! Write header information to binary punch file 
      ! Also be sure to trap I/O Error conditions
      !=================================================================
      WRITE ( IUNIT, IOSTAT=IOS ) FTI

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2_hdr:1' )

      WRITE ( IUNIT, IOSTAT=IOS ) TITLE

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2_hdr:2' )

      !=================================================================
      ! Return to calling program
      !=================================================================
      END SUBROUTINE BPCH2_HDR

!------------------------------------------------------------------------------

      SUBROUTINE BPCH2( IUNIT,     MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NTRACER,    
     &                  UNIT,      TAU0,      TAU1,     RESERVED,   
     &                  NI,        NJ,        NL,       IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY )
!
!******************************************************************************
!  Subroutine BPCH2 writes binary punch file (version 2.0) to disk.
!  Information about the model grid is also stored with each data block.
!  (bmy, 5/27/99, 9/12/00)
!
!  Arguments as input:
!  ============================================================================
!  (1    ) IUNIT      : INTEGER  - logical unit number of the file 
!  (2    ) MODELNAME  : CHAR*40  - Name of model used to create output
!  (3    ) LONRES     : REAL*4   - Longitude resolution of grid, in degrees
!  (4    ) LATRES     : REAL*4   - Latitude resolution of grid, in degrees
!  (4    ) HALFPOLAR  : INTEGER  - flag, =1 if model has half-polar boxes
!  (5    ) CENTER180  : INTEGER  - flag, =1 if model has lon center on 180 deg
!  (6    ) CATEGORY   : CHAR*40  - diagnostic category name
!  (7    ) NTRACER    : INTEGER  - number of tracer
!  (8    ) UNIT       : CHAR*40  - units of data
!  (9    ) TAU0       : REAL*8   - TAU at start of diagnostic interval
!  (10   ) TAU1       : REAL*8   - TAU at end   of diagnostic interval
!  (11   ) RESERVED   : CHAR*40  - Reserved for future use
!  (12-14) NI,NJ,NL   : INTEGER  - dimensions of ARRAY
!  (15   ) IFIRST     : INTEGER  - I-index of the first grid box
!  (16   ) JFIRST     : INTEGER  - J-index of the first grid box
!  (17   ) LFIRST     : INTEGER  - L-index of the first grid box
!  (18   ) ARRAY      : REAL*4   - data block to be written to the file
!
!  NOTES:
!  (1) Added indices to IOERROR calls (e.g. "bpch2:1", "bpch2:2", etc.) 
!       (bmy, 10/4/99)
!  (2) Added this routine to "bpch_mod.f" (bmy, 6/28/00)
!  (3) Use IOS /= 0 criterion to also check for EOF condition (bmy, 9/12/00)
!******************************************************************************
!  
      ! Arguments
      INTEGER,           INTENT(IN) :: IUNIT
      INTEGER,           INTENT(IN) :: NTRACER 
      INTEGER,           INTENT(IN) :: NI, NJ, NL 
      INTEGER,           INTENT(IN) :: IFIRST, JFIRST, LFIRST
      INTEGER,           INTENT(IN) :: HALFPOLAR, CENTER180
 
      REAL*4,            INTENT(IN) :: ARRAY( NI, NJ, NL )
      REAL*4,            INTENT(IN) :: LONRES, LATRES
      REAL*8,            INTENT(IN) :: TAU0,   TAU1

      CHARACTER(LEN=20), INTENT(IN) :: MODELNAME

      CHARACTER(LEN=40), INTENT(IN) :: CATEGORY
      CHARACTER(LEN=40), INTENT(IN) :: RESERVED
      CHARACTER(LEN=40), INTENT(IN) :: UNIT

      ! Local variables
      INTEGER                       :: I, J, L, NSKIP, IOS

      ! For computing NSKIP
      INTEGER, PARAMETER            :: BYTES_PER_NUMBER = 4
      INTEGER, PARAMETER            :: END_OF_RECORD    = 8

      !=================================================================
      ! BPCH2 begins here!!  
      !
      ! Compute the number of bytes to skip between the end of one 
      ! data block and the beginning of the next data header line
      !=================================================================
      NSKIP = ( BYTES_PER_NUMBER * ( NI * NJ * NL ) ) + END_OF_RECORD

      !=================================================================
      ! Write data block to binary punch file
      ! Check for I/O errors
      !=================================================================
      WRITE( IUNIT, IOSTAT=IOS ) 
     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2:1' )

      WRITE( IUNIT, IOSTAT = IOS ) 
     &     CATEGORY, NTRACER,  UNIT, TAU0,   TAU1,   RESERVED,
     &     NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &     NSKIP

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2:2' )

      WRITE( IUNIT, IOSTAT=IOS ) 
     &     ( ( ( ARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'bpch2:3' )

      !=================================================================
      ! Return to calling program      
      !=================================================================
      END SUBROUTINE BPCH2

!------------------------------------------------------------------------------

      SUBROUTINE READ_BPCH2( FILENAME, CATEGORY_IN, TRACER_IN, TAU0_IN, 
     &                       IX,       JX,          LX,        ARRAY ) 
!
!******************************************************************************
!  Subroutine READ_BPCH2 reads a binary punch file (v. 2.0) and extracts
!  a data block that matches the given category, tracer, and tau value.
!  (bmy, 12/10/99, 4/20/01)
!
!  Arguments as Input:
!  ============================================================================
!  (1  ) FILENAME    : (CHARACTER) String for input file name
!  (2  ) CATEGORY_IN : (CHARACTER) Category name for the desired data block
!  (3  ) TRACER_IN   : (INTEGER  ) Tracer number for which to extract data
!  (4  ) TAU0_IN     : (REAL*8   ) TAU value for which to extract data
!  (5-7) IX, JX, LX  : (INTEGER  ) Dimensions of ARRAY (see below) 
!
!  Arguments as Output:
!  ============================================================================
!  (8  ) ARRAY       : (REAL*4   ) Array to hold extracted data values
!
!  NOTES:
!  (1 ) Assumes that we are reading in a global-size data block.
!  (2 ) Trap all I/O errors with subroutine IOERROR.F.
!  (3 ) Now stop with an error message if no matches are found. (bmy, 3/9/00)
!  (4 ) Added this routine to "bpch_mod.f" (bmy, 6/28/00)
!  (5 ) Use IOS /= 0 criterion to also check for EOF condition (bmy, 9/12/00)
!  (6 ) TEMPARRAY now dimensioned to be of global size (bmy, 10/12/00) 
!  (7 ) Removed obsolete code from 10/12/00 (bmy, 12/18/00)
!  (8 ) Now make TEMPARRAY independent of CMN_SIZE parameters (bmy, 4/17/01)
!  (9 ) Removed old commented-out code (bmy, 4/20/01)
!******************************************************************************
!
      ! Arguments
      INTEGER,          INTENT(IN)  :: IX, JX, LX, TRACER_IN

      CHARACTER(LEN=*), INTENT(IN)  :: FILENAME, CATEGORY_IN 

      REAL*8,           INTENT(IN)  :: TAU0_IN
      REAL*4,           INTENT(OUT) :: ARRAY(IX, JX, LX)      

      ! Local variables
      INTEGER, PARAMETER :: IUNIT = 65
      LOGICAL            :: FOUND 

      INTEGER            :: I,  J,  L,  N,  IOS, M
      INTEGER            :: I1, I2, J1, J2, L1,  L2
      INTEGER            :: iii, jjj, kkk
      

      REAL*8             :: TIMEGAP !2015-11-16 LSURL
      ! Make TEMPARRAY big enough to for a 1x1 grid (bmy, 4/17/01)
      REAL*4             :: TEMPARRAY(360,181,70)

      ! For binary punch file, version 2.0
      INTEGER            :: NTRACER,   NSKIP
      INTEGER            :: HALFPOLAR, CENTER180
      INTEGER            :: NI,        NJ,        NL
      INTEGER            :: IFIRST,    JFIRST,    LFIRST

      REAL*4             :: LONRES,    LATRES

      REAL*8             :: ZTAU0,     ZTAU1

      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: CATEGORY
      CHARACTER(LEN=40)  :: UNIT     
      CHARACTER(LEN=40)  :: RESERVED
      CHARACTER(LEN=40)  :: FTI
      CHARACTER(LEN=80)  :: TITLE 

      !=================================================================
      ! READ_BPCH2 begins here!
      !  
      ! Initialize some variables
      !=================================================================
      print*,'Degub point 0.0'
      call flush(5)
      FOUND            = .FALSE.
      ARRAY(:,:,:)     = 0.0e0

 !     do iii = 1 , 360
 !           do jjj = 1 , 181
 !                do kkk = 1 , 70 ! innest is fastest
 !                      TEMPARRAY( iii , jjj , kkk ) = 0.0d0
 !                end do
 !           end do
 !     end do
      TEMPARRAY(:,:,:) = 0.0e0

      !======
      !Some debugging
      !print*,'entering READ_BPCH2'
      !call flush(5)
      !=================================================================
      ! Open binary punch file and read top-of-file header.
      ! Do some error checking to make sure the file is the right format.
      !=================================================================
      OPEN( IUNIT, FILE=TRIM( FILENAME ), STATUS='OLD',
     &             FORM='UNFORMATTED',    IOSTAT=IOS )

      IF ( IOS /= 0 ) CALL IOERROR( 
     &     IOS, IUNIT, 'read_bpch2:1' )

      ! Read file type identifier
      READ( IUNIT, IOSTAT=IOS ) FTI

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_bpch2:2' )

      IF ( TRIM( FTI ) /= 'CTM bin 02' ) THEN
         PRINT*, 'Input file is not in binary file format v. 2.0!'
         PRINT*, 'STOP in read_bpch2.f'
         STOP
      ENDIF
      
      !print*,'Read_bpch debug point 1.0'
      !call flush(5)
      
      ! Read top title
      READ( IUNIT, IOSTAT=IOS ) TITLE

      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_bpch2:3' )

      !=================================================================
      ! Read data from the binary punch file 
      !
      ! NOTE: IOS < 0 is end-of-file, IOS > 0 is error condition
      !=================================================================
      DO
         READ( IUNIT, IOSTAT=IOS ) 
     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
         
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IUNIT, 'read_bpch2:4' )

         READ( IUNIT, IOSTAT=IOS ) 
     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP

         !PRINT*, "Iteration"
 !!    &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
 !!    &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
 !!    &        NSKIP
         !call flush(5)

         IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_bpch2:5' )

         READ( IUNIT, IOSTAT=IOS ) 
     &        ( ( ( TEMPARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

         IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_bpch2:6' )
!! edit LSURL 2015-11-16. Allow TAU within a small tolerance range
         TIMEGAP = ABS(TAU0_IN - ZTAU0)
         ! Test for a match
         IF ( TRIM( CATEGORY_IN ) == TRIM( CATEGORY ) .and. 
     &        TRACER_IN           == NTRACER          .and.
     &        (TIMEGAP             .LE. 0.5)  ) THEN
!!     &        TAU0_IN             == ZTAU0 ) THEN
            FOUND = .TRUE.
            EXIT
         ENDIF

      ENDDO
      
      !print*,'Read_bpch debug point 2.0'
      !call flush(5)
      !=================================================================
      ! We have found a match!  Copy TEMPARRAY to ARRAY, taking into 
      ! account the starting positions (IFIRST, JFIRST, LFIRST) of 
      ! the data block.
      !=================================================================
      IF ( FOUND ) THEN 
         I1 = 1
         J1 = 1
        !I1 = IFIRST
        !J1 = JFIRST
         
         L1 = LFIRST
      
         I2 = NI - 1 !+ IFIRST
         J2 = NJ - 1 !+ JFIRST
         L2 = NL + LFIRST - 1
         

         !print*, I1
         !print*, I2
         !print*, J1
         !print*, J2
         !print*, NI
         !print*, NJ
         !call flush(5)
         ARRAY( I1:I2, J1:J2, L1:L2 ) = TEMPARRAY( 1:NI, 1:NJ, 1:NL )

         WRITE( 6, 100 ) ZTAU0, NTRACER
 100     FORMAT( 'READ_BPCH2: Found data for TAU = ', f10.2, 
     &           ' and tracer # ', i6 )
      ELSE
         WRITE( 6, 110 ) TRIM( FILENAME )
 110     FORMAT( 'READ_BPCH2: No matches found for file ', a )
         STOP
      ENDIF
      !print*,'Read_bpch debug point 3.0'
      !call flush(5)
      !=================================================================
      !  Close file and return to calling program  
      !=================================================================
      CLOSE( IUNIT )

      END SUBROUTINE READ_BPCH2

!------------------------------------------------------------------------------

      FUNCTION GET_MODELNAME() RESULT( MODELNAME )
!
!******************************************************************************
!  Function GET_MODELNAME returns the proper value of MODELNAME for GEOS-1,
!  GEOS-STRAT, GEOS-2, or GEOS-3 data.  MODELNAME is written to the binary
!  punch file and is used by the GAMAP package. (bmy, 6/22/00, 10/9/01)
!
!  NOTES:
!  (1 ) Now use special model name for GEOS-3 w/ 30 layers (bmy, 10/9/01)
!******************************************************************************
!

      ! MODELNAME holds the return value for the function
      CHARACTER(LEN=20) :: MODELNAME

#if   defined( GEOS_1 ) 
      MODELNAME = 'GEOS1'
     
#elif defined( GEOS_STRAT ) 
      MODELNAME = 'GEOS_STRAT'

#elif defined( GEOS_2 ) 
      MODELNAME = 'GEOS2'

#elif defined( GEOS_3 )
      !----------------------------------------------------------------- 
      ! Prior to 10/9/01:
      !MODELNAME = 'GEOS3'
      !----------------------------------------------------------------- 
      
      ! Write a special model name to the punch file for GAMAP
      ! if we are using regridded vertical resolution (bmy, 10/9/01)
      IF ( LLPAR == 30 ) THEN 
         MODELNAME = 'GEOS3_30L'
      ELSE
         MODELNAME = 'GEOS3'
      ENDIF

#elif defined( GEOS_4 )
      IF ( LLPAR == 30 ) THEN 
         MODELNAME = 'GEOS4_30L'
      ELSE
         MODELNAME = 'GEOS4'
      ENDIF
#elif defined( GEOS_5 )
      IF ( LLPAR == 47 ) THEN
         MODELNAME = 'GEOS5_47L'
      ELSE
         MODELNAME = 'GEOS5'
      ENDIF
#endif


      ! Return to calling program
      END FUNCTION GET_MODELNAME

!------------------------------------------------------------------------------

      FUNCTION GET_NAME_EXT() RESULT( NAME_EXT )
!
!******************************************************************************
!  Function GET_NAME_EXT returns the proper filename extension for CTM
!  model name (i.e. "geos1", "geoss", "geos2", or "geos3").  (bmy, 6/28/00)
!******************************************************************************
!
#     include "define.h"

      ! EXTENSION holds the return value for the function
      CHARACTER(LEN=5) :: NAME_EXT

#if   defined( GEOS_1 ) 
      NAME_EXT = 'geos1'
     
#elif defined( GEOS_STRAT ) 
      NAME_EXT = 'geoss'

#elif defined( GEOS_2 ) 
      NAME_EXT = 'geos2'

#elif defined( GEOS_3 )
      NAME_EXT = 'geos3'

#elif defined( GEOS_4 )
      NAME_EXT = 'geos4'
#elif defined( GEOS_5 )
      NAME_EXT = 'geos5'
#endif

      ! Return to calling program
      END FUNCTION GET_NAME_EXT

!------------------------------------------------------------------------------

      FUNCTION GET_RES_EXT() RESULT( RES_EXT )
!
!******************************************************************************
!  Function GET_RES_EXT returns the proper filename extension for
!  CTM grid resolution (i.e. "1x1", "2x25", "4x5").  (bmy, 6/28/00)
!******************************************************************************
!
#     include "define.h"

#if   defined( GRID4x5 )
      CHARACTER(LEN=3) :: RES_EXT
      RES_EXT = '4x5'
     
#elif defined( GRID2x25 ) 
      CHARACTER(LEN=4) :: RES_EXT
      RES_EXT = '2x25'

#elif defined( GRID1x1 ) 
      CHARACTER(LEN=3) :: RES_EXT
      RES_EXT = '1x1'

#elif defined( GRIDNEST )
      CHARACTER(LEN=9) :: RES_EXT
      RES_EXT = '025x03125'

#endif

      END FUNCTION GET_RES_EXT

!------------------------------------------------------------------------------

      FUNCTION GET_TAU0_2A( MONTH, YEAR ) RESULT( THIS_TAU0 )
!
!******************************************************************************
!  Function GET_TAU0_2A returns the corresponding TAU0 value for the first 
!  day of a given MONTH of a given YEAR.  This is necessary to index monthly 
!  mean binary punch files, which are used as input to GEOS-CHEM.  
!
!  NOTE: This function takes 2 arguments (MONTH, YEAR) and is the current 
!  default in the GEOS-CHEM code.  Overload this w/ an interface so that the 
!  user can also choose the version of GET_TAU0 w/ 3 arguments (MONTH, DAY,
!  YEAR), which will soon become the new standard. (bmy, 7/20/00, 8/22/01)
!
!  Arguments as Input:
!  ===========================================================================
!  (1) MONTH (INTEGER) : Month number (1-12)
!  (2) YEAR  (INTEGER) : 4-digit year number (e.g. 1985)
!
!  NOTES: 
!  (1 ) 1985 is the first year of the GEOS epoch.
!  (2 ) Add TAU0 values for years 1985-2001 (bmy, 8/1/00)
!  (3 ) Correct error for 1991 TAU values.  Also added 2002 and 2003.
!        (bnd, bmy, 1/4/01)
!*****************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN)           :: MONTH
      INTEGER, INTENT(IN), OPTIONAL :: YEAR

      ! Local variables
      INTEGER                       :: TMP_YEAR
      REAL*8                        :: TAU0(12)

      ! Return value
      REAL*8                        :: THIS_TAU0

      !=================================================================
      ! GET_TAU0_2A begins here!
      !=================================================================

      ! Error checking 
      IF ( MONTH < 1 .or. MONTH > 12 ) THEN
         WRITE( 6, '(a)' ) 'GET_TAU0: Invalid MONTH selection!'
         WRITE( 6, '(a)' ) 'STOP in GET_TAU0_VALUES ("bpch2_mod.f")!'
      ENDIF

      ! If YEAR is not passed, default to 1985 (first year of GEOS epoch)
      IF ( PRESENT( YEAR ) ) THEN
         TMP_YEAR = YEAR
      ELSE
         TMP_YEAR = 1985
      ENDIF
      
      ! CASE statement for year
      SELECT CASE ( TMP_YEAR )
         CASE ( 1985 )
            TAU0(:) = (/      0d0,    744d0,   1416d0,   2160d0, 
     &                     2880d0,   3624d0,   4344d0,   5088d0, 
     &                     5832d0,   6552d0,   7296d0,   8016d0 /)

         CASE ( 1986 )
            TAU0(:) = (/   8760d0,   9504d0,  10176d0,  10920d0, 
     &                    11640d0,  12384d0,  13104d0,  13848d0, 
     &                    14592d0,  15312d0,  16056d0,  16776d0 /)

         CASE ( 1987 )
            TAU0(:) = (/  17520d0,  18264d0,  18936d0,  19680d0,
     &                    20400d0,  21144d0,  21864d0,  22608d0,
     &                    23352d0,  24072d0,  24816d0,  25536d0 /)
 
         CASE ( 1988 )
            TAU0(:) = (/  26280d0,  27024d0,  27720d0,  28464d0,
     &                    29184d0,  29928d0,  30648d0,  31392d0,
     &                    32136d0,  32856d0,  33600d0,  34320d0 /)

         CASE ( 1989 )
            TAU0(:) = (/  35064d0,  35808d0,  36480d0,  37224d0,
     &                    37944d0,  38688d0,  39408d0,  40152d0,
     &                    40896d0,  41616d0,  42360d0,  43080d0 /)

         CASE ( 1990 )
            TAU0(:) = (/  43824d0,  44568d0,  45240d0,  45984d0,
     &                    46704d0,  47448d0,  48168d0,  48912d0,
     &                    49656d0,  50376d0,  51120d0,  51840d0 /)

         CASE ( 1991 )
            TAU0(:) = (/  52584d0,  53328d0,  54000d0,  54744d0,
     &                    55464d0,  56208d0,  56928d0,  57672d0,
     &                    58416d0,  59136d0,  59880d0,  60600d0 /)

         CASE ( 1992 )
            TAU0(:) = (/  61344d0,  62088d0,  62784d0,  63528d0,
     &                    64248d0,  64992d0,  65712d0,  66456d0,
     &                    67200d0,  67920d0,  68664d0,  69384d0 /)

         CASE ( 1993 )
            TAU0(:) = (/  70128d0,  70872d0,  71544d0,  72288d0,
     &                    73008d0,  73752d0,  74472d0,  75216d0,
     &                    75960d0,  76680d0,  77424d0,  78144d0 /)

         CASE ( 1994 ) 
            TAU0(:) = (/  78888d0,  79632d0,  80304d0,  81048d0,
     &                    81768d0,  82512d0,  83232d0,  83976d0,
     &                    84720d0,  85440d0,  86184d0,  86904d0 /)

         CASE ( 1995 )
            TAU0(:) = (/  87648d0,  88392d0,  89064d0,  89808d0,
     &                    90528d0,  91272d0,  91992d0,  92736d0,
     &                    93480d0,  94200d0,  94944d0,  95664d0 /)

         CASE ( 1996 )
            TAU0(:) = (/  96408d0,  97152d0,  97848d0,  98592d0,
     &                    99312d0, 100056d0, 100776d0, 101520d0,
     &                   102264d0, 102984d0, 103728d0, 104448d0 /)

         CASE ( 1997 )
            TAU0(:) = (/ 105192d0, 105936d0, 106608d0, 107352d0,
     &                   108072d0, 108816d0, 109536d0, 110280d0,
     &                   111024d0, 111744d0, 112488d0, 113208d0 /)

         CASE ( 1998 )
            TAU0(:) = (/ 113952d0, 114696d0, 115368d0, 116112d0,
     &                   116832d0, 117576d0, 118296d0, 119040d0,
     &                   119784d0, 120504d0, 121248d0, 121968d0 /)

         CASE ( 1999 )
            TAU0(:) = (/ 122712d0, 123456d0, 124128d0, 124872d0,
     &                   125592d0, 126336d0, 127056d0, 127800d0,
     &                   128544d0, 129264d0, 130008d0, 130728d0 /)

         CASE ( 2000 )
            TAU0(:) = (/ 131472d0, 132216d0, 132912d0, 133656d0,
     &                   134376d0, 135120d0, 135840d0, 136584d0,
     &                   137328d0, 138048d0, 138792d0, 139512d0 /)

         CASE ( 2001 )
            TAU0(:) = (/ 140256d0, 141000d0, 141672d0, 142416d0,
     &                   143136d0, 143880d0, 144600d0, 145344d0,
     &                   146088d0, 146808d0, 147552d0, 148272d0 /)
 
         CASE ( 2002 ) 
            TAU0(:) = (/ 149016d0, 149760d0, 150432d0, 151176d0,
     &                   151896d0, 152640d0, 153360d0, 154104d0,
     &                   154848d0, 155568d0, 156312d0, 157032d0 /)

         CASE ( 2003 ) 
            TAU0(:) = (/ 157776d0, 158520d0, 159192d0, 159936d0,
     &                   160656d0, 161400d0, 162120d0, 162864d0,
     &                   163608d0, 164328d0, 165072d0, 165792d0 /)

         CASE DEFAULT
            WRITE( 6, '(a)' ) 'GET_TAU0: Invalid YEAR selection!'
            WRITE( 6, '(a)' ) 'STOP in GET_TAU0 ("bpch2_mod.f")!'
            STOP

      END SELECT

      ! Select the TAU0 value for the given month
      THIS_TAU0 = TAU0( MONTH )

      ! Return to calling program
      END FUNCTION GET_TAU0_2A

!------------------------------------------------------------------------------

      FUNCTION GET_TAU0_6A( MONTH, DAY, YEAR, HOUR, MIN, SEC ) 
     &         RESULT( THIS_TAU0 )
!
!******************************************************************************
!  Function GET_TAU0_6A returns the corresponding TAU0 value for the first 
!  day of a given MONTH of a given YEAR.  This is necessary to index monthly 
!  mean binary punch files, which are used as input to GEOS-CHEM.
!  (bmy, 9/26/01) 
!
!  This function takes 3 mandatory arguments (MONTH, DAY, YEAR) and 3 
!  optional arguments (HOUR, MIN, SEC).  It is intended to replace the current 
!  2-argument version of GET_TAU0.  The advantage being that GET_TAU0_3A can 
!  compute a TAU0 for any date and time in the GEOS-CHEM epoch, rather than 
!  just the first day of each month.  Overload this w/ an interface so that 
!  the user can also choose the version of GET_TAU0 w/ 2 arguments 
!  (MONTH, YEAR), which is the prior version.
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) MONTH (INTEGER) : Month of year (1-12)
!  (2 ) DAY   (INTEGER) : Day of month (1-31)
!  (3 ) YEAR  (INTEGER) : 4-digit year number (e.g. 1985,2001)
!  (4 ) HOUR  (INTEGER) : OPTIONAL: Hour of day (0-24)
!  (5 ) MIN   (INTEGER) : OPTIONAL: Minute of hour (0-59)
!  (6 ) SEC   (INTEGER) : OPTIONAL: Seconds of minute (0-59)
!
!  NOTES: 
!  (1 ) 1985 is the first year of the GEOS epoch.
!  (2 ) Add TAU0 values for years 1985-2001 (bmy, 8/1/00)
!  (3 ) Correct error for 1991 TAU values.  Also added 2002 and 2003.
!        (bnd, bmy, 1/4/01)
!  (4 ) Updated comments  (bmy, 9/26/01)
!*****************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN)           :: MONTH
      INTEGER, INTENT(IN)           :: DAY
      INTEGER, INTENT(IN)           :: YEAR
      INTEGER, INTENT(IN), OPTIONAL :: HOUR
      INTEGER, INTENT(IN), OPTIONAL :: MIN
      INTEGER, INTENT(IN), OPTIONAL :: SEC

      ! Local variables
      INTEGER                       :: TMP_HOUR, TMP_MIN, TMP_SEC
      REAL*8                        :: DAYS

      ! Return value
      REAL*8                        :: THIS_TAU0
      
      !=================================================================
      ! GET_TAU0_3A begins here!
      !=================================================================

      ! Error checking 
      IF ( MONTH < 1 .or. MONTH > 12 ) THEN
         WRITE( 6, '(a)' ) 'GET_TAU0: Invalid MONTH selection!'
         WRITE( 6, '(a)' ) 'STOP in GET_TAU0 ("bpch2_mod.f")!'
         STOP
      ENDIF

      ! Error checking 
      IF ( DAY < 1 .or. DAY > 31 ) THEN
         write( *, * ) 'THIS IS THE DAY', DAY
         WRITE( 6, '(a)' ) 'GET_TAU0: Invalid DAY selection!'
         WRITE( 6, '(a)' ) 'STOP in GET_TAU0 ("bpch2_mod.f")!'
         STOP
      ENDIF

      ! If HOUR isn't passed, default to 0
      IF ( PRESENT( HOUR ) ) THEN
         TMP_HOUR = HOUR
      ELSE
         TMP_HOUR = 0
      ENDIF 

      ! If MIN isn't passed, default to 0
      IF ( PRESENT( MIN ) ) THEN
         TMP_MIN = MIN
      ELSE
         TMP_MIN = 0 
      ENDIF 

      ! If SEC isn't passed, default to 0
      IF ( PRESENT( SEC ) ) THEN
         TMP_SEC = SEC
      ELSE
         TMP_SEC = 0 
      ENDIF 

      ! Number of days since midnight on 1/1/1985
      THIS_TAU0 = JULDAY( YEAR, MONTH, DBLE( DAY ) ) - 2446066.5d0

      ! Multiply by 24 to get hours since 1/1/1985
      ! Also add in the hours elapsed since midnight on this date
      THIS_TAU0 = ( THIS_TAU0 * 24d0 ) + ( TMP_HOUR         ) + 
     &            ( TMP_MIN   / 60d0 ) + ( TMP_SEC / 3600d0 )

      ! Return to calling program
      END FUNCTION GET_TAU0_6A

!------------------------------------------------------------------------------

      FUNCTION JULDAY( YEAR, MONTH, DAY ) RESULT( JULIANDAY )
!
!******************************************************************************
!  Function JULDAY returns the astronomical Julian day.
!
!  Algorithm taken from "Practical Astronomy With Your Calculator",
!  Third Edition, by Peter Duffett-Smith, Cambridge UP, 1992.
! 
!  Arguments as Input:
!  ------------------------------------------------------------------
!  (1 ) YEAR  : (INTEGER) Current year
!  (2 ) MONTH : (INTEGER) Current month
!  (3 ) DAY   : (REAL*8 ) Current day (can be fractional, e.g. 17.25)
!
!  NOTES:
!  (1 ) JULDAY requires the external function MINT.F.
!  (2 ) JULDAY will compute the correct Julian day for any 
!        BC or AD date.
!  (3 ) For BC dates, subtract 1 from the year and append a minus 
!        sign.  For example, 1 BC is 0, 2 BC is -1, etc.  This is 
!        necessary for the algorithm.  
!******************************************************************************
!   
      ! Arguments
      INTEGER, INTENT(IN) :: year, month
      REAL*8              :: day,  julianDay
   
      ! Local variables
      INTEGER             :: year1, month1
      REAL*8              :: x1, a, b, c, d
      LOGICAL             :: isGregorian
   
      !==================================================================
      ! JULDAY begins here!
      !
      ! Follow algorithm from Peter Duffett-Smith (1992)
      !==================================================================
   
      ! Compute YEAR and MONTH1
      IF ( ( MONTH == 1 ) .OR. ( MONTH == 2 ) ) THEN
         YEAR1  = YEAR  - 1
         MONTH1 = MONTH + 12 
      ELSE
         YEAR1  = YEAR
         MONTH1 = MONTH
      ENDIF
   
      ! Compute the "A" term. 
      X1 = DBLE( YEAR ) / 100.0D0
      A  = MINT( X1 )
   
      ! The Gregorian calendar begins on 10 October 1582
      ! Any dates prior to this will be in the Julian calendar
      IF ( YEAR > 1582 ) THEN
         ISGREGORIAN = .TRUE.
      ELSE
         IF ( ( YEAR   == 1582 )  .AND. 
     &        ( MONTH1 >= 10   )  .AND. 
     &        ( DAY    >= 15.0 ) ) THEN 
            ISGREGORIAN = .TRUE.
         ELSE
            ISGREGORIAN = .FALSE.
         ENDIF
      ENDIF
            
      ! Compute the "B" term according to Gregorian or Julian calendar
      IF ( ISGREGORIAN ) THEN
         B = 2.0D0 - A + MINT( A / 4.0D0 )
      ELSE
         B = 0.0D0
      ENDIF
   
      ! Compute the "C" term for BC dates (YEAR1 <= 0 ) 
      ! or AD dates (YEAR1 > 0)
      IF ( YEAR1 < 0 ) THEN
         X1 = ( 365.25D0 * YEAR1 ) - 0.75D0
         C  = MINT( X1 )
      ELSE
         X1 = 365.25D0 * YEAR1
         C  = MINT( X1 ) 
      ENDIF
   
      ! Compute the "D" term    
      X1 = 30.6001D0 * DBLE( MONTH1 + 1 )
      D  = MINT( X1 )
   
      ! Add the terms to get the Julian Day number 
      JULIANDAY = B + C + D + DAY + 1720994.5D0
   
      END FUNCTION JULDAY
   
!------------------------------------------------------------------------------
   
      FUNCTION MINT( X ) RESULT ( VALUE )
!
!******************************************************************************
!  Function MINT is defined as follows:
! 
!  MINT = -INT( ABS( X ) ), X <  0
!  MINT =  INT( ABS( X ) ), X >= 0
! 
!  Arguments as Input:
!  ------------------------------------------------------------------
!  (1) X : (REAL*8) Argument for the function MINT
!
!  NOTES:
!  (1) MINT is primarily intended for use with routine JULDAY.
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: X
         
      ! Return value
      REAL*8             :: value
   
      !=================================================================
      ! MINT begins here!
      !=================================================================
      IF ( X < 0D0 ) THEN 
         VALUE = -INT( ABS( X ) )        
      ELSE
         VALUE =  INT( ABS( X ) )        
      ENDIF
   
      ! Return to calling program
      END FUNCTION MINT
   
!------------------------------------------------------------------------------
   
      SUBROUTINE CALDATE( JULDAY, NYMD, NHMS )
!
!******************************************************************************
!  Subroutine CALDATE converts an astronomical Julian day to 
!  the NYMD (e.g. YYYYMMDD) and NHMS (i.e. HHMMSS) format.
!   
!  Algorithm taken from "Practical Astronomy With Your Calculator",
!  Third Edition, by Peter Duffett-Smith, Cambridge UP, 1992.
!   
!  Arguments as Input:
!  ============================================================================
!  (1) JULDAY : REAL*8  : Astronomical julian day
!   
!  Arguments as Output:
!  ============================================================================
!  (1) NYMD   : INTEGER : YYYYMMDD corresponding to JDAY
!  (2) NHMS   : INTEGER : HHMMSS   corresponding to JDAY
!******************************************************************************
!
      ! Arguments
      REAL*8,  INTENT(IN)  :: JULDAY
      INTEGER, INTENT(OUT) :: NYMD, NHMS
    
      ! Local variables
      REAL*8               :: A, B, C, D, DAY, E, F 
      REAL*8               :: FDAY, G, I, J, JD, M, Y
    
      !=================================================================
      ! CALDATE begins here!
      ! See "Practical astronomy with your calculator", Peter Duffett-
      ! Smith 1992, for an explanation of the following algorithm.
      !=================================================================
      JD = JULDAY + 0.5D0
      I  = INT( JD )
      F  = JD - INT( I )
    
      IF ( I > 2299160D0 ) THEN
         A = INT( ( I - 1867216.25D0 ) / 36524.25 )
         B = I + 1 + A - INT( A / 4 )
      ELSE
         B = I
      ENDIF
    
      C = B + 1524D0
      
      D = INT( ( C - 122.1D0 ) / 365.25D0 )
    
      E = INT( 365.25D0 * D )
    
      G = INT( ( C - E ) / 30.6001D0 )
   
      ! DAY is the day number
      DAY  = C - E + F - INT( 30.6001D0 * G ) 
    
      ! FDAY is the fractional day number
      FDAY = DAY - INT( DAY )
      
      ! M is the month number
      IF ( G < 13.5D0 ) THEN
         M = G - 1D0
      ELSE
         M = G - 13D0
      ENDIF
    
      ! Y is the year number
      IF ( M > 2.5D0 ) THEN
         Y = D - 4716D0
      ELSE
         Y = D - 4715D0
      ENDIF
    
      ! NYMD is YYYYMMDD
      NYMD = ( INT( Y ) * 10000 ) + ( INT( M ) * 100 ) + INT( DAY )
      
      ! NHMS is HHMMSS
      NHMS = INT( FDAY * 24 ) * 10000 
    
      ! Return to calling program
      END SUBROUTINE CALDATE
   
!------------------------------------------------------------------------------

      END MODULE BPCH2_MOD
