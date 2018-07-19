! $Id: He5SwathModule.f90,v 1.2 2006/01/20 15:50:56 bmy Exp $
MODULE He5SwathModule

  !========================================================================
  ! Module He5SwathModule contains routines for reading data from swath
  ! data structures in HDF-EOS5 data files. (bmy, 1/20/06)
  !
  ! Module Variables:
  ! -----------------------------------------------------------------------
  ! (1 ) VERBOSE       (LOGICAL  ) : Flag for toggling verbose output
  ! (2 ) dataTypeName  (CHARACTER) : Array w/ names of HDF-EOS5 data types
  ! (3 ) saveFileName  (CHARACTER) : Shadow variable for filename
  ! (4 ) saveSwathName (CHARACTER) : Shadow variable for swath name
  !
  ! Module Routines:
  ! -----------------------------------------------------------------------
  ! (1 ) He5VerboseOutput          : Toggles verbose output for file I/O
  ! (2 ) He5FileOpen               : Opens HDF-EOS5 file; gets file ID #
  ! (3 ) He5FileClose              : Closes HDF-EOS5 file
  ! (4 ) He5SwathAttach            : Attaches to swath; gets swath ID #
  ! (5 ) He5SwathDetach            : Detaches from swath
  ! (6 ) He5SwathDimInfo           : Gets dimension names, types, sizes
  ! (7 ) He5SwathGeoFldInfo        : Gets info about swath geoloc fields
  ! (8 ) He5SwathDataFldInfo       : Gets info about swath data fields
  ! (9 ) He5SwathFldInfo           : Gets info about an individual field
  ! (10) He5SwathFillValue         : Gets missing data fill values
  ! (11) He5SwathAttrs             : Gets attributes from HDF-EOS5 swath
  ! (12) He5SwathReadData1dI2      : Reads 1-D INTEGER*2 data array 
  ! (13) He5SwathReadData1dI4      : Reads 1-D INTEGER*4 data array 
  ! (14) He5SwathReadData1dR4      : Reads 1-D REAL*4    data array 
  ! (15) He5SwathReadData1dR8      : Reads 1-D REAL*8    data array 
  ! (16) He5SwathReadData2dI2      : Reads 2-D INTEGER*2 data array 
  ! (17) He5SwathReadData2dI4      : Reads 2-D INTEGER*4 data array 
  ! (18) He5SwathReadData2dR4      : Reads 2-D REAL*4    data array 
  ! (19) He5SwathReadData2dR8      : Reads 2-D REAL*8    data array 
  ! (20) He5SwathReadData3dI2      : Reads 3-D INTEGER*2 data array 
  ! (21) He5SwathReadData3dI4      : Reads 3-D INTEGER*4 data array 
  ! (22) He5SwathReadData3dR4      : Reads 3-D REAL*4    data array 
  ! (23) He5SwathReadData3dR8      : Reads 3-D REAL*8    data array 
  ! (24) makeCharArrayFromCharList : Splits char list into char array
  ! (25) He5DataTypeName           : Returns data type name from type #
  !
  ! Module Interfaces:
  ! -----------------------------------------------------------------------
  ! (1 ) He5ReadSwathData          : Overloads the following routines:
  !                                   (a) He5SwathReadData1dI2
  !                                   (b) He5SwathReadData1dI4      
  !                                   (c) He5SwathReadData1dR4
  !                                   (d) He5SwathReadData1dR8      
  !                                   (e) He5SwathReadData2dI2      
  !                                   (f) He5SwathReadData2dI4      
  !                                   (g) He5SwathReadData2dR4      
  !                                   (h) He5SwathReadData2dR8      
  !                                   (i) He5SwathReadData3dI2      
  !                                   (j) He5SwathReadData3dI4      
  !                                   (k) He5SwathReadData3dR4      
  !                                   (l) He5SwathReadData3dR8      
  !
  ! Other Information:
  ! -----------------------------------------------------------------------
  ! (1 ) The data type HE5-INTEGER (represented by parameter HE5_INT in
  !       He5IncludeModule) is either:
  !       (a) INTEGER*4 (for 32-byte machines such as Linux Boxes)
  !       (b) INTEGER*8 (for 64-byte machines such as Altix & SGI Origin)
  ! (2 ) You must select your machine type in the file "He5Define.h".
  !       This will automatically set parameter HE5_INT accordingly.
  ! (3 ) Data arrays which are passed to the HDF-EOS5 library function
  !       HE5_SwRdFld must be dimensioned with values of type HE5_INT.
  !
  ! References:
  ! -----------------------------------------------------------------------
  ! (1 ) http://hdf.ncsa.uiuc.edu/HDF5/
  !         -- HDF5 home page
  ! (2 ) http://newsroom.gsfc.nasa.gov/sdptoolkit/toolkit.html
  !         -- ECS toolkit home page (home of HDF-EOS5)
  !
  ! NOTES:
  !========================================================================

  ! References to F90 modules
  USE He5ErrorModule
  USE He5IncludeModule

  ! Force explicit data types
  IMPLICIT NONE

  !------------------------------------------------------------------------
  ! PRIVATE / PUBLIC DECLARATIONS
  !------------------------------------------------------------------------

  ! Make everything PRIVATE ...
  PRIVATE

  ! ... except these routines
  PUBLIC :: He5VerboseOutput
  PUBLIC :: He5FileOpen
  PUBLIC :: He5FileClose
  PUBLIC :: He5SwathAttach
  PUBLIC :: He5SwathDetach
  PUBLIC :: He5SwathDimInfo
  PUBLIC :: He5SwathGeoFldInfo
  PUBLIC :: He5SwathDataFldInfo
  PUBLIC :: He5SwathFldInfo
  PUBLIC :: He5SwathFillValue
  PUBLIC :: He5SwathAttrs
  PUBLIC :: He5SwathReadData

  !------------------------------------------------------------------------
  ! MODULE VARIABLES 
  !------------------------------------------------------------------------
  LOGICAL                     :: VERBOSE      = .FALSE.
  CHARACTER(LEN=HE5_MAX_CHAR) :: dataTypeName(57)
  CHARACTER(LEN=HE5_MAX_CHAR) :: saveFileName
  CHARACTER(LEN=HE5_MAX_CHAR) :: saveSwathName
  
  !------------------------------------------------------------------------
  ! MODULE INTERFACES
  !------------------------------------------------------------------------
  INTERFACE He5SwathReadData
     MODULE PROCEDURE He5SwathReadData1dI2
     MODULE PROCEDURE He5SwathReadData1dI4
     MODULE PROCEDURE He5SwathReadData1dR4
     MODULE PROCEDURE He5SwathReadData1dR8
     MODULE PROCEDURE He5SwathReadData2dI2
     MODULE PROCEDURE He5SwathReadData2dI4
     MODULE PROCEDURE He5SwathReadData2dR4
     MODULE PROCEDURE He5SwathReadData2dR8
     MODULE PROCEDURE He5SwathReadData3dI2
     MODULE PROCEDURE He5SwathReadData3dI4
     MODULE PROCEDURE He5SwathReadData3dR4
     MODULE PROCEDURE He5SwathReadData3dR8
  END INTERFACE

CONTAINS

!------------------------------------------------------------------------------

  SUBROUTINE He5VerboseOutput( v )

    !======================================================================
    ! Subroutine He5VerboseOutput is used to trigger "extra" output from
    ! the routines in this module. (bmy, 1/4/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) v (LOGICAL) : TRUE or FALSE value
    !
    ! NOTES:
    !======================================================================

    ! Arguments 
    LOGICAL, INTENT(IN) :: v

    ! Set the value of verbose
    VERBOSE = v

  END SUBROUTINE He5VerboseOutput

!------------------------------------------------------------------------------

  SUBROUTINE He5FileOpen( fileName, fId )

    !======================================================================
    ! Subroutine He5FileOpen opens an HDF-EOS5 file and returns the
    ! file Id number. (bmy, 1/4/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) fileName (CHARACTER) : Name of HDF-EOS5 file to open
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (2) fId      (INTEGER  ) : HDF-EOS5 file ID number
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    CHARACTER(LEN=*), INTENT(IN)  :: fileName
    INTEGER,          INTENT(OUT) :: fId      

    ! Local variables
    CHARACTER(LEN=HE5_MAX_CHAR)   :: msg, loc

    ! HDF-EOS5 library routines
    INTEGER                       :: HE5_SwOpen

    !--------------------------
    ! He5FileOpen begins here!
    !--------------------------

    ! Store filename in a shadow variable
    saveFileName = fileName

    ! Open HDF-EOS5 file and get file ID #
    fId          = HE5_SwOpen( fileName, HE5F_ACC_RDONLY )

    ! Error check
    IF ( fId == FAILURE ) THEN
       msg = 'Error opening file ' // TRIM( saveFileName )
       loc = 'He5FileOpen ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF
 
    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( saveFileName )
       WRITE( 6, 110 ) fId 
100    FORMAT( '===> HDF-EOS5 file name : "', a, '"' )
110    FORMAT( '===> HDF-EOS5 file ID   : ', i10     )
    ENDIF

  END SUBROUTINE He5FileOpen

!------------------------------------------------------------------------------

  SUBROUTINE He5FileClose( fId )

    !======================================================================
    ! Subroutine He5FileClose closes an HDF-EOS5 file. (bmy, 1/4/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) fId (INTEGER) : HDF-EOS5 file ID number
    !
    ! NOTES:
    !======================================================================
    
    ! Arguments
    INTEGER, INTENT(IN)         :: fId      

    ! Local variables
    INTEGER                     :: status
    CHARACTER(LEN=HE5_MAX_CHAR) :: msg, loc

    ! HDF-EOS5 library routines
    INTEGER                     :: HE5_SwClose

    !---------------------------
    ! He5FileClose begins here!
    !---------------------------

    ! Get HDF-EOS5 file ID
    status = HE5_SwClose( fId )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error closing file ' // TRIM( savefileName )
       loc = 'He5FileClose ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( saveFileName )
100    FORMAT( '===> Closed file "', a, '"' )
    ENDIF

  END SUBROUTINE He5FileClose

!------------------------------------------------------------------------------

  SUBROUTINE He5SwathAttach( fId, swathName, sId )

    !======================================================================
    ! Subroutine He5SwathAttach attaches to an HDF-EOS5 swath data
    ! structure and returns the swath ID number. (bmy, 1/4/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) fId       (INTEGER)   : HDF-EOS5 file ID (see He5FileOpen)
    ! (2) swathName (CHARACTER) : Name of HDF-EOS5 swath to attach to
    ! 
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (3) sId       (INTEGER)   : HDF-EOS5 swath ID number
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: fId
    CHARACTER(LEN=*), INTENT(IN)  :: swathName
    INTEGER,          INTENT(OUT) :: sId  

    ! Local variables
    CHARACTER(LEN=HE5_MAX_CHAR)   :: msg, loc

    ! HDF-EOS5 library routines
    INTEGER                       :: HE5_SwAttach

    !-----------------------------
    ! He5SwathAttach begins here!
    !-----------------------------

    ! Save swathname in a shadow variable
    saveSwathName = swathName

    ! Attach to swath
    sId           = HE5_SwAttach( fId, swathName )

    ! Error check
    IF ( sId == FAILURE ) THEN
       msg = 'Error attaching to swath ' // TRIM( saveSwathName )
       loc = 'He5SwathAttach ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF
  
    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( saveSwathName )
       WRITE( 6, 110 ) sId 
100    FORMAT( '===> HDF-EOS5 swath name: "', a, '"' )
110    FORMAT( '===> HDF-EOS5 swath ID  : ', i10     )  
    ENDIF

  END SUBROUTINE He5SwathAttach

!------------------------------------------------------------------------------

  SUBROUTINE He5SwathDetach( sId )

    !======================================================================
    ! Subroutine He5SwathDetach detaches from an HDF-EOS5 swath
    ! data structure. (bmy, 1/4/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId (INTEGER)   : HDF-EOS5 swath ID number (see He5SwathAttach)
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER, INTENT(IN)         :: sId  

    ! Local variables
    INTEGER                     :: status
    CHARACTER(LEN=HE5_MAX_CHAR) :: msg, loc

    ! HDF-EOS5 library routines
    INTEGER                     :: HE5_SwDetach

    !-----------------------------
    ! He5SwathDetach begins here!
    !-----------------------------

    ! Detach from swath
    status = HE5_SwDetach( sId )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error detaching from swath ' // TRIM( saveSwathName )
       loc = 'He5SwathDetach ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( saveSwathName )
100    FORMAT( '===> Detached from swath "', a, '"' )
    ENDIF

  END SUBROUTINE He5SwathDetach

!------------------------------------------------------------------------------

  SUBROUTINE He5SwathDimInfo( sId, nDims, dims, dimNames )

    !======================================================================
    ! Subroutine He5SwathDetach detaches from an HDF-EOS5 swath
    ! data structure. (bmy, 1/4/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId (INTEGER)   : HDF-EOS5 swath ID number (see He5SwathAttach)
    !
    ! NOTES:
    !======================================================================
    
    ! Arguments
    INTEGER,                     INTENT(IN)  :: sId  
    INTEGER,                     INTENT(OUT) :: nDims  
    INTEGER,                     INTENT(OUT) :: dims(HE5_MAX_DIMS) 
    CHARACTER(LEN=HE5_MAX_CHAR), INTENT(OUT) :: dimNames(HE5_MAX_DIMS)
    
    ! Local variables
    INTEGER                                  :: N, C
    CHARACTER(LEN=HE5_MAX_CHAR)              :: msg, loc, dimList

    ! HDF-EOS5 library routines
    INTEGER                                  :: HE5_SwInqDims

    !------------------------------
    ! He5SwathDimInfo begins here!
    !------------------------------

    ! Initialize
    nDims       = 0
    dims(:)     = 0
    dimNames(:) = ''

    ! Get dimension info for this swath
    nDims       = HE5_SwInqDims( sId, dimList, dims )

    ! Make an array from the dimension list
    CALL makeCharArrayFromCharList( dimList, ',', dimNames )

    ! NOTE: Sometimes every other element of DIMS is zero.  
    ! I don't know why but we can just pack the array to be 
    ! on the safe side. (bmy, 1/4/06)
    C = 0
    DO N = 1, 2*nDims+1 
       IF ( dims(N) > 0  ) THEN
          C       = C + 1
          dims(C) = dims(N)
          IF ( N > 1 ) dims(N) = 0
       ENDIF
    ENDDO

    ! Error check
    IF ( nDims <= 0 ) THEN
       msg = 'Error getting dim info from swath ' // TRIM( saveSwathName )
       loc = 'He5SwathDetach ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) nDims

       DO N = 1, nDims 
          WRITE( 6, 110 ) TRIM( dimNames(N) ), dims(N)
       ENDDO

100    FORMAT( '===> There are ', i4, ' dimensions'      )
110    FORMAT( '===> ', a25,' is of size ', i10 )
    ENDIF

  END SUBROUTINE He5SwathDimInfo

!------------------------------------------------------------------------------

  SUBROUTINE He5SwathGeoFldInfo( sId, nGeo, geoRank, geoName, geoType )

    !======================================================================
    ! Subroutine He5SwathGeoFieldInfo obtains information about the 
    ! geolocation fields in the HDF-EOS5 swath data structure. 
    ! (bmy, 1/4/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId     (INTEGER  ) : HDF-EOS5 swath ID number 
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (2) nGeo    (INTEGER  ) : Number of geolocation fields
    ! (3) geoRank (INTEGER  ) : Number of dimensions for each geoloc field
    ! (4) geoName (CHARACTER) : Name of each geolocation field
    ! (5) geoType (CHARACTER) : Data type of each geolocation field
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,                     INTENT(IN)  :: sId
    INTEGER,                     INTENT(OUT) :: nGeo
    INTEGER,                     INTENT(OUT) :: geoRank(HE5_MAX_FLDS)
    CHARACTER(LEN=HE5_MAX_CHAR), INTENT(OUT) :: geoName(HE5_MAX_FLDS)
    CHARACTER(LEN=HE5_MAX_CHAR), INTENT(OUT) :: geoType(HE5_MAX_FLDS)

    ! Local variables
    INTEGER                                  :: N
    INTEGER                                  :: typeNum(HE5_MAX_FLDS)
    CHARACTER(LEN=HE5_MAX_CHAR)              :: msg, loc, geoList

    ! HDF-EOS5 library routines
    INTEGER                                  :: HE5_SwInqGFlds

    !---------------------------------
    ! He5SwathGeoFldInfo begins here!
    !---------------------------------

    ! Initialize
    nGeo       = 0
    geoRank(:) = 0
    geoName(:) = ''
    geoType(:) = ''

    ! Get number of geo fields and related info
    nGeo       = HE5_SwInqGFlds( sId, geoList, geoRank, typeNum )

    ! Error check
    IF ( nGeo <= 0 ) THEN
       msg = 'Error getting geolocation field information!'
       loc = 'He5SwathGeoFldInfo ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Separate list of field names into an array
    CALL makeCharArrayFromCharList( geoList, ',', geoName )

    ! Get HDF-EOS5 data type names for each data type number
    DO N = 1, nGeo
       geoType(N) = He5DataTypeName( typeNum(N) )
    ENDDO

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) nGeo

       DO N = 1, nGeo
          WRITE( 6, 110 ) TRIM( geoName(N) ), geoRank(N), TRIM( geoType(N) )
       ENDDO

100    FORMAT( '===> There are ', i4, ' Geolocation Fields' )
110    FORMAT( '===> ', a25, ' has ', i4 , ' dimensions and is ', a )
    ENDIF

  END SUBROUTINE He5SwathGeoFldInfo
    
!------------------------------------------------------------------------------

  SUBROUTINE He5SwathDataFldInfo( sId, nData, dataRank, dataName, dataType )

    !======================================================================
    ! Subroutine He5SwathDataFieldInfo obtains information about the 
    ! data fields in the HDF-EOS5 swath data structure. (bmy, 1/4/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId      (INTEGER  ) : HDF-EOS5 swath ID number 
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (2) nData    (INTEGER  ) : Number of data fields
    ! (3) dataRank (INTEGER  ) : Number of dimensions for each data field
    ! (4) dataName (CHARACTER) : Name of each data field
    ! (5) dataType (CHARACTER) : Data type of each data field
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,                     INTENT(IN)  :: sId
    INTEGER,                     INTENT(OUT) :: ndata
    INTEGER,                     INTENT(OUT) :: dataRank(HE5_MAX_FLDS)
    CHARACTER(LEN=HE5_MAX_CHAR), INTENT(OUT) :: dataName(HE5_MAX_FLDS)
    CHARACTER(LEN=HE5_MAX_CHAR), INTENT(OUT) :: dataType(HE5_MAX_FLDS)

    ! Local variables
    INTEGER                                  :: N
    INTEGER                                  :: typeNum(HE5_MAX_FLDS)
    CHARACTER(LEN=HE5_MAX_CHAR)              :: msg, loc, dataList

    ! HDF-EOS5 library routines
    INTEGER                                  :: HE5_SwInqDFlds

    !---------------------------------
    ! He5SwathGeoFldInfo begins here!
    !---------------------------------

    ! Initialize
    nData       = 0
    dataRank(:) = 0
    dataName(:) = ''
    dataType(:) = ''

    ! Get number of data fields and related info
    nData       = HE5_SwInqDFlds( sId, dataList, dataRank, typeNum )

    ! Error check
    IF ( nData <= 0 ) THEN
       msg = 'Error getting data field information!'
       loc = 'He5SwathDataFldInfo ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Separate list of field names into an array
    CALL makeCharArrayFromCharList( dataList, ',', dataName )

    ! Get HDF-EOS5 data type names for each data type number
    DO N = 1, nData
       dataType(N) = He5DataTypeName( typeNum(N) )
    ENDDO

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) nData

       DO N = 1, nData
          WRITE( 6, 110 ) TRIM( dataName(N) ), dataRank(N), TRIM( dataType(N) )
       ENDDO

100    FORMAT( '===> There are ', i6, ' Data Fields' )
110    FORMAT( '===> ', a40, ' has ', i4 , ' dimensions and is ', a )
    ENDIF

  END SUBROUTINE He5SwathDataFldInfo
    
!------------------------------------------------------------------------------

  SUBROUTINE He5SwathFldInfo( sId, name, typeName, rank, dims, dimNames )
    
    !======================================================================
    ! Subroutine He5SwathFldInfo obtains information about a particular
    ! data field in the HDF-EOS5 swath data structure. (bmy, 1/4/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId      (INTEGER  ) : HDF-EOS5 swath ID number 
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (2) nData    (INTEGER  ) : Number of data fields
    ! (3) typeName (CHARACTER) : Name of the data type for this field
    ! (4) rank     (INTEGER  ) : Number of dimensions for each data field
    ! (4) dims     (INTEGER  ) : Integer containing field dimensions
    ! (5) dimNames (CHARACTER) : Array containing names of each dimension
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,                     INTENT(IN)  :: sId
    CHARACTER(LEN=*),            INTENT(IN)  :: name
    INTEGER,                     INTENT(OUT) :: rank
    INTEGER,                     INTENT(OUT) :: dims(HE5_MAX_DIMS)
    CHARACTER(LEN=HE5_MAX_CHAR), INTENT(OUT) :: dimNames(HE5_MAX_DIMS)
    CHARACTER(LEN=HE5_MAX_CHAR), INTENT(OUT) :: typeName

    ! Local variables
    INTEGER                                  :: C,   N,   status,  typeNum
    CHARACTER(LEN=HE5_MAX_CHAR)              :: msg, loc, dimList, maxList

    ! HDF-EOS5 library routines
    INTEGER                                  :: HE5_SwFldInfo

    !------------------------------
    ! He5SwathFldInfo begins here!
    !------------------------------

    ! Initialize
    rank     = 0
    dims     = 0
    dimNames = ''

    ! Get number of data fields and related info
    status = HE5_SwFldInfo( sId, name, rank, dims, typeNum, dimList, maxList )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error getting info about field ' // TRIM( name )
       loc = 'He5SwathFldInfo ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Separate list of dimension names into an array
    CALL makeCharArrayFromCharList( dimList, ',', dimNames )

    ! Get HDF-EOS5 data type name
    typeName = He5DataTypeName( typeNum )

    ! NOTE: Sometimes every other element of DIMS is zero.  
    ! I don't know why but we can just pack the array to be 
    ! on the safe side. (bmy, 1/4/06)
    C = 0
    DO N = 1, HE5_MAX_DIMS
       IF ( dims(N) > 0 ) THEN
          C       = C + 1
          dims(C) = dims(N)
          IF ( N > 1 ) dims(N) = 0
       ENDIF
    ENDDO

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( name ), C, TRIM( typeName )
       WRITE( 6, 110 ) dims(1:C)
100    FORMAT( '===> ', a25 ' has ', i4 , ' dimensions and is ', a )
110    FORMAT( '===> ', 25x,' its dimensions are: ', 10i7 )
    ENDIF

  END SUBROUTINE He5SwathFldInfo

!------------------------------------------------------------------------------

  SUBROUTINE He5SwathFillValue( sId, dataName, dataFill )
    
    !======================================================================
    ! Subroutine He5SwathFillValue reads the missing data "fill" value
    ! for a field contained in the HDF-EOS5 swath data structure. 
    ! (bmy, 1/4/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId      (INTEGER  ) : HDF-EOS5 swath ID 
    ! (2) dataName (CHARACTER) : Name of the field to get fill value for
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------   
    ! (3) dataFill (REAL*4   ) : Fill value for missing data
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,                     INTENT(IN)  :: sId
    CHARACTER(LEN=HE5_MAX_CHAR), INTENT(IN)  :: dataName
    REAL*4,                      INTENT(OUT) :: dataFill

    ! Local variables
    INTEGER                                  :: status

    ! HDF-EOS5 library routines
    INTEGER                                  :: HE5_SwGetFill

    !--------------------------------
    ! He5SwathFillValue begins here!
    !--------------------------------

    ! Get the fill value
    status = HE5_SwGetFill( sId, TRIM( dataName ), dataFill )

    ! Set fill value to zero if it does not exist
    IF ( status == FAILURE ) dataFill = 0.0

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( dataName ), dataFill
100    FORMAT( '===> Fill value for ', a, ' is ', es13.6 )
    ENDIF

  END SUBROUTINE He5SwathFillValue
 
!-----------------------------------------------------------------------------

  SUBROUTINE He5SwathAttrs( sId, nAttrs, attrName, attrValue )
  
    !======================================================================
    ! Subroutine He5SwathAttrs returns the global attributes associated
    ! with the swath data structure. (bmy, 1/4/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1) sId       (INTEGER  ) : HDF-EOS5 swath ID number 
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (2) nAttrs    (INTEGER  ) : Number of swath attributes
    ! (3) attrName  (CHARACTER) : Array of attribute names
    ! (4) attrValue (CHARACTER) : Array of attribute values
    ! 
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,                     INTENT(IN)  :: sId
    INTEGER,                     INTENT(OUT) :: nAttrs
    CHARACTER(LEN=HE5_MAX_CHAR), INTENT(OUT) :: attrName(HE5_MAX_ATRS)
    REAL*4,                      INTENT(OUT) :: attrValue(HE5_MAX_ATRS)

    ! Local variables
    INTEGER                                  :: N, status, strBufSize
    REAL*4                                   :: value
    CHARACTER(LEN=HE5_MAX_CHAR)              :: attrList
    
    ! HDF_EOS5 library routines
    INTEGER                                  :: HE5_SwInqAttrs
    INTEGER                                  :: HE5_SwRdAttr

    !=====================================================================
    ! "gridGetAttributes" begins here!
    !=====================================================================

    ! Get list of attribute names
    nAttrs = HE5_SwInqAttrs( sId, attrList, strBufSize )

    ! Separate list into array
    CALL makeCharArrayFromCharList( attrList, ',', attrName )

    ! Get the data value for each attribute
    ! For each attribute
    DO N = 1, nAttrs
       status = HE5_SwRdAttr( sId, TRIM( attrName(N) ), attrValue )       
    ENDDO
    
    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) nAttrs

       DO N = 1, nAttrs
          WRITE( 6, 110 ) TRIM( attrName(N) ), attrValue(N)
       ENDDO

100    FORMAT( '===> There are ', i6, ' Swath Attributes' )
110    FORMAT( '===> ', a20, ' has value ', es13.6 )
    ENDIF

  END SUBROUTINE He5SwathAttrs

!-----------------------------------------------------------------------------

  SUBROUTINE He5SwathReadData1dI2( sId, fldName, nX, fldData )

    !======================================================================
    ! Routine He5SwathReadData1dI2 reads a 1-D INTEGER*2 data block from
    ! an HDF-EOS5 swath data structure.  This routine is included in the
    ! module interface He5SwathReadData. (bmy, 1/13/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS5 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE5-INTEGER) : 1st dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (4 ) fldData (INTEGER*2  ) : 1-D array w/ data from HDF-EOS5 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE5_INT), INTENT(IN)  :: nX
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    INTEGER*2,        INTENT(OUT) :: fldData(nX)

    ! Local variables
    INTEGER                       :: status
    INTEGER(HE5_INT)              :: start(1), stride(1), edge(1)
    CHARACTER(LEN=HE5_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS5 library routines
    INTEGER                       :: HE5_SwRdFld
    
    !-----------------------------------
    ! He5SwathReadData1dI2 begins here!
    !-----------------------------------

    ! Set up to read data for a given track
    start  = (/  0 /) 
    stride = (/  1 /)
    edge   = (/ nX /)

    ! Read data
    status = HE5_SwRdFld( sId, fldName, start, stride, edge, fldData )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData1dI2 ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF

  END SUBROUTINE He5SwathReadData1dI2

!-----------------------------------------------------------------------------

  SUBROUTINE He5SwathReadData1dI4( sId, fldName, nX, fldData )

    !======================================================================
    ! Routine He5SwathReadData1dI4 reads a 1-D INTEGER*4 data block from
    ! an HDF-EOS5 swath data structure.  This routine is included in the
    ! module interface He5SwathReadData. (bmy, 1/13/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS5 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE5-INTEGER) : 1st dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (4 ) fldData (INTEGER*2  ) : 1-D array w/ data from HDF-EOS5 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE5_INT), INTENT(IN)  :: nX
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    INTEGER*4,        INTENT(OUT) :: fldData(nX)

    ! Local variables
    INTEGER                       :: status
    INTEGER(HE5_INT)              :: start(1), stride(1), edge(1)
    CHARACTER(LEN=HE5_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS5 library routines
    INTEGER                       :: HE5_SwRdFld
    
    !-----------------------------------
    ! He5SwathReadData1dI2 begins here!
    !-----------------------------------

    ! Set up to read data for a given track
    start  = (/  0 /) 
    stride = (/  1 /)
    edge   = (/ nX /)

    ! Read data
    status = HE5_SwRdFld( sId, fldName, start, stride, edge, fldData )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData1dI4 ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF

  END SUBROUTINE He5SwathReadData1dI4

!-----------------------------------------------------------------------------

  SUBROUTINE He5SwathReadData1dR4( sId, fldName, nX, fldData )

    !======================================================================
    ! Routine He5SwathReadData1dR4 reads a 1-D REAL*4 data block from
    ! an HDF-EOS5 swath data structure.  This routine is included in the
    ! module interface He5SwathReadData. (bmy, 1/13/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS5 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE5-INTEGER) : 1st dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (4 ) fldData (REAL*4     ) : 1-D array w/ data from HDF-EOS5 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE5_INT), INTENT(IN)  :: nX
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    REAL*4,           INTENT(OUT) :: fldData(nX)

    ! Local variables
    INTEGER                       :: status
    INTEGER(HE5_INT)              :: start(1), stride(1), edge(1)
    CHARACTER(LEN=HE5_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS5 library routines
    INTEGER                       :: HE5_SwRdFld
    
    !-----------------------------------
    ! He5SwathReadData1dR4 begins here!
    !-----------------------------------

    ! Set up to read data for a given track
    start  = (/  0 /) 
    stride = (/  1 /)
    edge   = (/ nX /)

    ! Read data
    status = HE5_SwRdFld( sId, fldName, start, stride, edge, fldData )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData1dR4 ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF

  END SUBROUTINE He5SwathReadData1dR4

!------------------------------------------------------------------------------

  SUBROUTINE He5SwathReadData1dR8( sId, fldName, nX, fldData )

    !======================================================================
    ! Routine He5SwathReadData1dR8 reads a 1-D REAL*8 data block from
    ! an HDF-EOS5 swath data structure.  This routine is included in the
    ! module interface He5SwathReadData. (bmy, 1/13/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS5 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE5-INTEGER) : 1st dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (4 ) fldData (REAL*8     ) : 1-D array w/ data from HDF-EOS5 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)    :: sId
    INTEGER(HE5_INT), INTENT(IN)    :: nX
    CHARACTER(LEN=*), INTENT(IN)    :: fldName
    REAL*8,           INTENT(INOUT) :: fldData(nX)

    ! Local variables
    INTEGER                         :: status    
    INTEGER(HE5_INT)                :: start(1), stride(1), edge(1)
    CHARACTER(LEN=HE5_MAX_CHAR)     :: msg,      loc

    ! HDF-EOS5 library routines
    INTEGER                         :: HE5_SwRdFld
    
    !-----------------------------------
    ! He5SwathReadData1dR8 begins here!
    !-----------------------------------

    ! Set up dimension info
    start  = (/  0 /) 
    stride = (/  1 /)
    edge   = (/ nX /)

    ! Read data
    status = HE5_SwRdFld( sId, fldName, start, stride, edge, fldData )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData1dR8 ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF

  END SUBROUTINE He5SwathReadData1dR8

!------------------------------------------------------------------------------

  SUBROUTINE He5SwathReadData2dI2( sId, fldName, nX, nY, fldData )

    !======================================================================
    ! Routine He5SwathReadData2dI2 reads a 2-D INTEGER*2 data block from
    ! an HDF-EOS5 swath data structure.  This routine is included in the
    ! module interface He5SwathReadData. (bmy, 1/13/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS5 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE5-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE5-INTEGER) : 2nd dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (5 ) fldData (INTEGER*2  ) : 2-D array w/ data from HDF-EOS5 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE5_INT), INTENT(IN)  :: nX, nY
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    INTEGER*2,        INTENT(OUT) :: fldData(nX,nY)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE5_INT)              :: start(2), stride(2), edge(2)
    CHARACTER(LEN=HE5_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS5 library routines
    INTEGER                       :: HE5_SwRdFld
    
    !-----------------------------------
    ! He5SwathReadData2dR4 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0 /) 
    stride = (/  1,  1 /)
    edge   = (/ nX, nY /)
    
    ! Read data
    status = HE5_SwRdFld( sId, fldName, start, stride, edge, fldData )
    
    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData2dI2 ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF
       
  END SUBROUTINE He5SwathReadData2dI2

!------------------------------------------------------------------------------

  SUBROUTINE He5SwathReadData2dI4( sId, fldName, nX, nY, fldData )

    !======================================================================
    ! Routine He5SwathReadData2dI4 reads a 2-D INTEGER*4 data block from
    ! an HDF-EOS5 swath data structure.  This routine is included in the
    ! module interface He5SwathReadData. (bmy, 1/13/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS5 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE5-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE5-INTEGER) : 2nd dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (5 ) fldData (INTEGER*4  ) : 2-D array w/ data from HDF-EOS5 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE5_INT), INTENT(IN)  :: nX, nY
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    INTEGER*4,        INTENT(OUT) :: fldData(nX,nY)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE5_INT)              :: start(2), stride(2), edge(2)
    CHARACTER(LEN=HE5_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS5 library routines
    INTEGER                       :: HE5_SwRdFld
    
    !-----------------------------------
    ! He5SwathReadData2dR4 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0 /) 
    stride = (/  1,  1 /)
    edge   = (/ nX, nY /)
    
    ! Read data
    status = HE5_SwRdFld( sId, fldName, start, stride, edge, fldData )
    
    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData2dI4 ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF
       
  END SUBROUTINE He5SwathReadData2dI4

!------------------------------------------------------------------------------

  SUBROUTINE He5SwathReadData2dR4( sId, fldName, nX, nY, fldData )

    !======================================================================
    ! Routine He5SwathReadData2dR4 reads a 2-D REAL*4 data block from
    ! an HDF-EOS5 swath data structure.  This routine is included in the
    ! module interface He5SwathReadData. (bmy, 1/13/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS5 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE5-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE5-INTEGER) : 2nd dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (5 ) fldData (REAL*4     ) : 2-D array w/ data from HDF-EOS5 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE5_INT), INTENT(IN)  :: nX, nY
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    REAL*4,           INTENT(OUT) :: fldData(nX,nY)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE5_INT)              :: start(2), stride(2), edge(2)
    CHARACTER(LEN=HE5_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS5 library routines
    INTEGER                       :: HE5_SwRdFld
    
    !-----------------------------------
    ! He5SwathReadData2dR4 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0 /) 
    stride = (/  1,  1 /)
    edge   = (/ nX, nY /)
    
    ! Read data
    status = HE5_SwRdFld( sId, fldName, start, stride, edge, fldData )
    
    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData2dR4 ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF
       
  END SUBROUTINE He5SwathReadData2dR4

!-----------------------------------------------------------------------------

  SUBROUTINE He5SwathReadData2dR8( sId, fldName, nX, nY, fldData )

    !======================================================================
    ! Routine He5SwathReadData2dR8 reads a 2-D REAL*8 data block from
    ! an HDF-EOS5 swath data structure.  This routine is included in the
    ! module interface He5SwathReadData. (bmy, 1/13/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS5 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE5-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE5-INTEGER) : 2nd dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (5 ) fldData (REAL*8     ) : 2-D array w/ data from HDF-EOS5 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE5_INT), INTENT(IN)  :: nX, nY
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    REAL*8,           INTENT(OUT) :: fldData(nX,nY)

    ! Local variables
    INTEGER                       :: status  
    INTEGER(HE5_INT)              :: start(2), stride(2), count(2)
    CHARACTER(LEN=HE5_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS5 library routines
    INTEGER                       :: HE5_SwRdFld
    
    !-----------------------------------
    ! He5SwathReadData2dR8 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0 /) 
    stride = (/  1,  1 /)
    count  = (/ nX, nY /)

    ! Read data
    status = HE5_SwRdFld( sId, fldName, start, stride, count, fldData )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData2dR8 ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF

  END SUBROUTINE He5SwathReadData2dR8

!------------------------------------------------------------------------------

  SUBROUTINE He5SwathReadData3dI2( sId, fldName, nX, nY, nZ, fldData )

    !======================================================================
    ! Routine He5SwathReadData3dI2 reads a 3-D INTEGER*2 data block from
    ! an HDF-EOS5 swath data structure.  This routine is included in the
    ! module interface He5SwathReadData. (bmy, 1/13/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS5 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE5-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE5-INTEGER) : 2nd dimension of data array
    ! (5 ) nZ      (HE5-INTEGER) : 2nd dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (6 ) fldData (INTEGER*2  ) : 2-D array w/ data from HDF-EOS5 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE5_INT), INTENT(IN)  :: nX, nY, nZ
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    INTEGER*2,        INTENT(OUT) :: fldData(nX,nY,nZ)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE5_INT)              :: start(3), stride(3), edge(3)
    CHARACTER(LEN=HE5_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS5 library routines
    INTEGER                       :: HE5_SwRdFld
    
    !-----------------------------------
    ! He5SwathReadData3dI2 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0,  0 /) 
    stride = (/  1,  1,  1 /)
    edge   = (/ nX, nY, nZ /)
    
    ! Read data
    status = HE5_SwRdFld( sId, fldName, start, stride, edge, fldData )
    
    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData3dI2 ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF
       
  END SUBROUTINE He5SwathReadData3dI2

!------------------------------------------------------------------------------

  SUBROUTINE He5SwathReadData3dI4( sId, fldName, nX, nY, nZ, fldData )

    !======================================================================
    ! Routine He5SwathReadData3dI4 reads a 3-D INTEGER*4 data block from
    ! an HDF-EOS5 swath data structure.  This routine is included in the
    ! module interface He5SwathReadData. (bmy, 1/13/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS5 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE5-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE5-INTEGER) : 2nd dimension of data array
    ! (5 ) nZ      (HE5-INTEGER) : 2nd dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (6 ) fldData (INTEGER*4  ) : 2-D array w/ data from HDF-EOS5 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE5_INT), INTENT(IN)  :: nX, nY, nZ
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    INTEGER*4,        INTENT(OUT) :: fldData(nX,nY,nZ)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE5_INT)              :: start(3), stride(3), edge(3)
    CHARACTER(LEN=HE5_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS5 library routines
    INTEGER                       :: HE5_SwRdFld
    
    !-----------------------------------
    ! He5SwathReadData3dI2 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0,  0 /) 
    stride = (/  1,  1,  1 /)
    edge   = (/ nX, nY, nZ /)
    
    ! Read data
    status = HE5_SwRdFld( sId, fldName, start, stride, edge, fldData )
    
    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData3dI4 ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF
       
  END SUBROUTINE He5SwathReadData3dI4

!------------------------------------------------------------------------------

  SUBROUTINE He5SwathReadData3dR4( sId, fldName, nX, nY, nZ, fldData )

    !======================================================================
    ! Routine He5SwathReadData1dR4 reads a 2-D REAL*4 data block from
    ! an HDF-EOS5 swath data structure.  This routine is included in the
    ! module interface He5SwathReadData. (bmy, 1/13/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS5 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE5-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE5-INTEGER) : 2nd dimension of data array
    ! (5 ) nZ      (HE5-INTEGER) : 2nd dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (6 ) fldData (REAL*4     ) : 2-D array w/ data from HDF-EOS5 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE5_INT), INTENT(IN)  :: nX, nY, nZ
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    REAL*4,           INTENT(OUT) :: fldData(nX,nY,nZ)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE5_INT)              :: start(3), stride(3), edge(3)
    CHARACTER(LEN=HE5_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS5 library routines
    INTEGER                       :: HE5_SwRdFld
    
    !-----------------------------------
    ! He5SwathReadData3dR4 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0,  0 /) 
    stride = (/  1,  1,  1 /)
    edge   = (/ nX, nY, nZ /)
    
    ! Read data
    status = HE5_SwRdFld( sId, fldName, start, stride, edge, fldData )
    
    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData3dR4 ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF
       
  END SUBROUTINE He5SwathReadData3dR4

!-----------------------------------------------------------------------------

  SUBROUTINE He5SwathReadData3dR8( sId, fldName, nX, nY, nZ, fldData )

    !======================================================================
    ! Routine He5SwathReadData3dR8 reads a 3-D REAL*8 data block from
    ! an HDF-EOS5 swath data structure.  This routine is included in the
    ! module interface He5SwathReadData. (bmy, 1/13/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) sId     (INTEGER    ) : HDF-EOS5 Swath ID # 
    ! (2 ) fldName (CHARACTER  ) : Name of data array
    ! (3 ) nX      (HE5-INTEGER) : 1st dimension of data array
    ! (4 ) nY      (HE5-INTEGER) : 2nd dimension of data array
    ! (5 ) nZ      (HE5-INTEGER) : 3rd dimension of data array
    !
    ! Arguments as Output:
    ! ---------------------------------------------------------------------
    ! (6 ) fldData (REAL*8     ) : 3-D array w/ data from HDF-EOS5 file
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    INTEGER,          INTENT(IN)  :: sId
    INTEGER(HE5_INT), INTENT(IN)  :: nX, nY, nZ
    CHARACTER(LEN=*), INTENT(IN)  :: fldName
    REAL*8,           INTENT(OUT) :: fldData(nX,nY,nZ)

    ! Local variables
    INTEGER                       :: status   
    INTEGER(HE5_INT)              :: start(3), stride(3), edge(3)
    CHARACTER(LEN=HE5_MAX_CHAR)   :: msg,      loc

    ! HDF-EOS5 library routines
    INTEGER                       :: HE5_SwRdFld
    
    !-----------------------------------
    ! He5SwathReadData3dR8 begins here!
    !-----------------------------------

    ! Set up to read entire field
    start  = (/  0,  0,  0 /) 
    stride = (/  1,  1,  1 /)
    edge   = (/ nX, nY, nZ /)

    ! Read data
    status = HE5_SwRdFld( sId, fldName, start, stride, edge, fldData )

    ! Error check
    IF ( status == FAILURE ) THEN
       msg = 'Error reading data for ' // TRIM( fldName )
       loc = 'HdfSwathReadData3dR8 ("He5SwathModule.f90")'
       CALL He5ErrMsg( msg, loc )
    ENDIF

    ! Verbose output
    IF ( VERBOSE ) THEN
       WRITE( 6, 100 ) TRIM( fldName )
100    FORMAT( '===> Successfully read data for ', a )
    ENDIF

  END SUBROUTINE He5SwathReadData3dR8

!-----------------------------------------------------------------------------

  SUBROUTINE makeCharArrayFromCharList( list, separator, array )

    !=====================================================================
    ! Subroutine makeCharArrayFromCharList takes a comma-separated word 
    ! list, and places each word into a separate element of a character 
    ! array. (bmy, 1/4/06)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) list      (CHARACTER) : String with comma-separated words
    ! (2) separator (CHARACTER) : String for separator text
    !
    ! Arguments as output:
    ! --------------------------------------------------------------------
    ! (3) array     (CHARACTER) : Array of substrings
    !
    ! NOTES:
    !=====================================================================

    ! Arguments
    CHARACTER(LEN=HE5_MAX_CHAR), INTENT(IN)  :: list
    CHARACTER(LEN=1           ), INTENT(IN)  :: separator
    CHARACTER(LEN=HE5_MAX_CHAR), INTENT(OUT) :: array(:)

    ! local variables
    INTEGER                                  :: P, N, ind(255)
    CHARACTER(LEN=1)                         :: C

    !----------------------------------------
    ! makeCharArrayFromCharList begins here!
    !----------------------------------------

    ! Initialize
    N      = 1
    ind(:) = 0

    ! Find the positions of all the commas in LIST
    DO P = 1, LEN( list )

       ! Look at each character individually
       C = list(P:P)

       ! If a comma... 
       IF ( C == separator ) THEN 

          ! Increment comma
          N      = N + 1
          ind(N) = P 
       ENDIF
    ENDDO

    ! Add the position of the end of the string into IND
    ind(N+1) = LEN( list )

    ! Save text between the commas into ARRAY
    DO P = 1, N
       IF ( P == N ) THEN 
          array(P) = list( ind(P)+1:ind(P+1)   )
       ELSE
          array(P) = list( ind(P)+1:ind(P+1)-1 )
       ENDIF
    ENDDO

  END SUBROUTINE makeCharArrayFromCharList

!------------------------------------------------------------------------------

  FUNCTION He5DataTypeName( nType ) RESULT( typeStr )

    !=====================================================================
    ! Subroutine He5DataTypeName returns a descriptive string given a
    ! HDF-EOS5 data type number. (bmy, 1/4/06)
    !
    ! Arguments as Input:
    ! --------------------------------------------------------------------
    ! (1) nType     (INTEGER) : HDF-EOS number type 
    !
    ! NOTES:
    !=====================================================================

    ! Arguments
    INTEGER, INTENT(IN)         :: nType

    ! Local varaibles
    LOGICAL, SAVE               :: FIRST = .TRUE.
    CHARACTER(LEN=HE5_MAX_CHAR) :: typeStr

    !------------------------------
    ! He5DataTypeName begins here!
    !------------------------------

    ! First-time initialization
    IF ( FIRST ) THEN 
       dataTypeName(:)                  = ''
       dataTypeName(HE5T_NATIVE_SHORT ) = 'INTEGER*2'
       dataTypeName(HE5T_NATIVE_USHORT) = 'Unsigned INTEGER*2'
       dataTypeName(HE5T_NATIVE_LONG  ) = 'INTEGER*4'
       dataTypeName(HE5T_NATIVE_ULONG ) = 'Unsigned INTEGER*4'
       dataTypeName(HE5T_NATIVE_FLOAT ) = 'REAL*4'
       dataTypeName(HE5T_NATIVE_DOUBLE) = 'REAL*8'
       dataTypeName(HE5T_NATIVE_INT8  ) = 'INTEGER*8'
       dataTypeName(HE5T_NATIVE_HBOOL ) = 'LOGICAL'
       dataTypeName(HE5T_NATIVE_CHAR  ) = 'CHARACTER'
       dataTypeName(HE5T_CHARSTRING   ) = 'CHARACTER'
    ENDIF
    
    ! Return value
    typeStr = dataTypeName(nType)

  END FUNCTION He5DataTypeName

!------------------------------------------------------------------------------

END MODULE He5SwathModule









