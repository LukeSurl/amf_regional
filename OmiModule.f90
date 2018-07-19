! $Id$
MODULE OmiModule

  !========================================================================
  ! Module OmiModule contains variables and routines for reading OMI
  ! satellite swath data from HDF-EOS5 files (dbm, bmy, 1/20/06)
  !
  ! Module Var
  !========================================================================

  ! References to F90 modules
  Use He5ErrorModule
  Use He5IncludeModule
  USE He5SwathModule
  USE ParameterModule,  ONLY: NVALMAX

  ! Force explicit data types
  IMPLICIT NONE

  !------------------------------------------------------------------------
  ! PRIVATE / PUBLIC DECLARATIONS
  !------------------------------------------------------------------------

  ! Make everything Private by default, first subroutines contained are Public 
  PRIVATE
  PUBLIC :: OMI_GetAmfInput

  !------------------------------------------------------------------------
  ! MODULE VARIABLES
  !------------------------------------------------------------------------

  ! Fill value
  REAL*4                 :: dataFill

  ! Allocatable arrays
  INTEGER                :: as

  ! Cloud file geolocation fields. For NO2 these are in the Trace Gas File
  REAL*4,    ALLOCATABLE :: RelativeAzimuthAngle(:,:)

  ! Cloud file data fields. For NO2 these are in the Trace Gas File
  INTEGER*2,    ALLOCATABLE :: CloudPressure       (:,:)
  REAL*4,       ALLOCATABLE :: CloudFraction       (:,:)

  ! Data file geolocation fields
  REAL*8,    ALLOCATABLE :: OmiTime               (:)
  REAL*4,    ALLOCATABLE :: SolarZenithAngle   (:,:)
  REAL*4,    ALLOCATABLE :: ViewingZenithAngle (:,:)
  REAL*4,    ALLOCATABLE :: Latitude           (:,:)
  REAL*4,    ALLOCATABLE :: Longitude          (:,:)

  ! All Trace Gas file data fields
  INTEGER*2, ALLOCATABLE :: MainDataQualityFlag        (:,:)

  ! HCHO file data fields
  REAL*8,    ALLOCATABLE :: ColumnAmount               (:,:)
  REAL*8,    ALLOCATABLE :: ColumnUncertainty          (:,:)
  REAL*4,    ALLOCATABLE :: PixelCornerLatitudes       (:,:)
  REAL*4,    ALLOCATABLE :: PixelCornerLongitudes      (:,:)
  REAL*8,    ALLOCATABLE :: ColumnAmountDestriped      (:,:)
  REAL*8,    ALLOCATABLE :: FittingRMS                 (:,:)
  REAL*8,    ALLOCATABLE :: AirMassFactor              (:,:)
  INTEGER*2, ALLOCATABLE :: AirMassFactorDiagnosticFlag(:,:)

  ! NO2 file data fields
  INTEGER*2, ALLOCATABLE :: CloudFractionNO2        (:,:)
  REAL*4,    ALLOCATABLE :: ColumnAmountNO2Trop     (:,:)
  REAL*4,    ALLOCATABLE :: ColumnUncertaintyNO2Trop(:,:)
  REAL*4,    ALLOCATABLE :: AMFPolluted             (:,:)
  REAL*4,    ALLOCATABLE :: AMFUnpolluted           (:,:)
  REAL*4,    ALLOCATABLE :: ColumnAmountNO2BelowCloudTrop(:,:)
  REAL*4,    ALLOCATABLE :: ColumnUncertaintyNO2BelowCloudTrop(:,:)


  ! Saved variables for swath dimensions
  INTEGER*2, SAVE        :: OMI_dimX
  INTEGER*2, SAVE        :: OMI_dimY

  !------------------------------------------------------------------------
  ! MODULE ROUTINES
  !------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------
    

  SUBROUTINE OMI_GetAmfInput(satellite_file,TG_TYPE,CLDFILE,LAT_centers,LON_centers,SZA,SVA,AZM,PCLD,FCLD,MAXLINE,outPrefix,outSuffix)
    
    !======================================================================
    ! OMI_GetAmfInput extracts required inputs for AMF code from OMI file
    ! (GOB, 10/23/09)
    !
    ! NOTES: This routine was added to OmiModule.f90 so that satelliteIO.f90
    !        could be compiled without HDF libraries in future versions.
    !======================================================================

     ! Inputs
     !---------------------------------------------------------------------------------
     CHARACTER(LEN=255),intent(in)  :: satellite_file    ! The file to extract data from
     CHARACTER(LEN=255),intent(in)  :: CLDFILE           ! Cloud information in seperate file
     INTEGER, intent(in)            :: TG_TYPE           ! Trace Gas, 0=HCHO, 1=NO2, 2=SO2

     ! Outputs
     !---------------------------------------------------------------------------------
     !Variables Needed for all AMF calculations
      REAL, intent(out):: SZA(NVALMAX)                    ! Solar Zenith Angle
      REAL, intent(out):: SVA(NVALMAX)                    ! Satellite viewing angle
      REAL, intent(out):: AZM(NVALMAX)                    ! Relative azimuth angle, Not currently used
      REAL, intent(out):: LAT_centers(NVALMAX)            ! Latitude of pixel center
      REAL, intent(out):: LON_centers(NVALMAX)            ! Longitude of pixel center
      REAL, intent(out):: PCLD(NVALMAX)                   ! Cloud Pressure (hPa)
      REAL, intent(out):: FCLD(NVALMAX)                   ! Cloud Fraction (between 0 and 1)
      INTEGER, intent(out) :: MAXLINE                       ! The number of observations read in
      CHARACTER(LEN=200),intent(out) :: outSuffix(NVALMAX)  ! A string to pass on data from satellite file to AMF output
      CHARACTER(LEN=100),intent(out) :: outPrefix(NVALMAX)  ! A string to pass on data from satellite file to AMF output


     ! Other Variables
     !--------------------------------------------------------------------------------
      INTEGER:: X1,Y1
      INTEGER:: NSCAN,NPIX,LINE

     ! Start OMI_GetAmfInput
     !--------------------------------------------------------------------------------

! Read OMI cloud and Trace Gas Data files

      if (TG_TYPE .eq. 0) then
        !HCHO

         CALL HCHO_OmiReadCld (CLDFILE,'CloudFractionAndPressure')
         CALL HCHO_OmiReadData(satellite_file,'OMI Total Column Amount HCHO')

      endif

      if (TG_TYPE .eq. 1) then
        !NO2

         CALL NO2_OmiReadData( satellite_file, 'ColumnAmountNO2')
      endif


     ! Cast from INTEGER::*2 to INTEGER::*4
      X1      = OMI_dimX
      Y1      = OMI_dimY
      MAXLINE = X1 * Y1


      !Transfer OMI data from array to vector
       DO NSCAN = 1, Y1
       DO NPIX = 1, X1
 
         ! 1D array index
          LINE = ((NSCAN - 1) * X1 ) + NPIX

         !Read in required outputs
         !--------------------------------------------------------------
          SZA(LINE) =         SolarZenithAngle           (NPIX, NSCAN)
          SVA(LINE) =         ViewingZenithAngle         (NPIX, NSCAN)
          LAT_centers(LINE) = Latitude                   (NPIX, NSCAN)
          LON_centers(LINE) = Longitude                  (NPIX, NSCAN)
!          FLAG_SAO(LINE) =   MainDataQualityFlag        (NPIX, NSCAN)
          PCLD(LINE) =        1.0D0*CloudPressure        (NPIX, NSCAN)
          AZM(LINE) =         RelativeAzimuthAngle       (NPIX, NSCAN)

          if (TG_TYPE .eq. 0) then
             FCLD(LINE) =     CloudFraction              (NPIX, NSCAN)
          elseif (TG_TYPE .eq. 1) then
             FCLD(LINE) =     0.001D0*CloudFractionNO2   (NPIX, NSCAN)
          endif

         !Other optional outputs available to be
         ! passed along through outPrefix & outSufix
         !--------------------------------------------------------
          
 	  !OmiTime                    (NSCAN)

          !HCHO
 !           ColumnAmount               (NPIX, NSCAN)
 !           ColumnUncertainty          (NPIX, NSCAN)
 !           PixelCornerLongitudes      (NPIX, NSCAN)
 !           PixelCornerLongitudes    (NPIX+1, NSCAN)
 !           PixelCornerLongitudes  (NPIX+1, NSCAN+1)
 !           PixelCornerLongitudes    (NPIX, NSCAN+1)
 
 !           PixelCornerLatitudes       (NPIX, NSCAN)
 !           PixelCornerLatitudes     (NPIX+1, NSCAN)
 !           PixelCornerLatitudes   (NPIX+1, NSCAN+1)
 !           PixelCornerLatitudes     (NPIX, NSCAN+1)

 !           FittingRMS                 (NPIX, NSCAN)
 !           ColumnAmountDestriped      (NPIX, NSCAN)
 !           AirMassFactor              (NPIX, NSCAN)
 !           AirMassFactorDiagnosticFlag(NPIX, NSCAN)

          !NO2
 !           ColumnAmountNO2Trop      (NPIX, NSCAN)
 !           ColumnUncertaintyNO2Trop (NPIX, NSCAN)
 !           ColumnAmountNO2BelowCloudTrop (NPIX, NSCAN)
 !           ColumnUncertaintyNO2BelowCloudTrop(NPIX, NSCAN)
 !           AMFUnpolluted            (NPIX, NSCAN)
 !           AMFPolluted              (NPIX, NSCAN)

          write(outPrefix(LINE),'(2i6)') NSCAN, NPIX
       ENDDO
       ENDDO

       outSuffix(:) = ' suffix'

    end subroutine OMI_GetAmfInput


!------------------------------------------------------------------------------

  SUBROUTINE HCHO_OmiReadCld( fileName, swathName )

    !======================================================================
    ! HCHO_OmiReadCld reads cloud fields from an OMI HDF-EOS5 swath file.
    ! (bmy, 1/20/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) fileName  (CHARACTER) : Name of HDF-EOS5 file
    ! (2 ) swathName (CHARACTER) : Name of swath in HDF-EOS5 file:
    ! 
    ! NOTES:
    !======================================================================

    ! Arguments
    CHARACTER(LEN=*) , INTENT(IN) :: fileName, swathName
    
    ! Local variables
    INTEGER                       :: fId,     sId,     as  
    INTEGER                       :: fldRank, fldDims(HE5_MAX_DIMS)
    INTEGER(HE5_INT)              :: nX,      nY
    CHARACTER(LEN=HE5_MAX_CHAR)   :: fldType, name
    CHARACTER(LEN=HE5_MAX_CHAR)   :: fldDimNames(HE5_MAX_DIMS)
   
    !-----------------------------
    ! HCHO_OmiReadCld begins here!
    !-----------------------------

    ! Turn off verbose output
    CALL He5VerboseOutput( .FALSE. )

    ! Open HDF-EOS5 swath and get the file ID number
    CALL He5FileOpen( fileName, fId )

     ! Attach to swath and get swath ID number
    CALL He5SwathAttach( fId, swathName, sId )

    ! Echo info
    WRITE( 6, '(a)' ) REPEAT( '-', 79 )
    WRITE( 6, 100   ) TRIM( fileName )
    WRITE( 6, 110   ) TRIM( swathName )
100 FORMAT( 'HCHO_OmiReadCld: Reading ',       a )
110 FORMAT( 'HCHO_OmiReadCld: Swath name is ', a ) 

    !-----------------------------
    ! RelativeAzimuthAngle
    !-----------------------------

    ! Field name
    name = 'ViewingAzimuthAngle'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
    
    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)
    
    ! Allocate array
    ALLOCATE( RelativeAzimuthAngle( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, RelativeAzimuthAngle )
    
    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !-----------------------------
    ! CloudFraction
    !-----------------------------

    ! Field name
    name = 'CloudFraction'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
    
    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)
    
    ! Allocate array
    ALLOCATE( CloudFraction( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, CloudFraction )
    
    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !---------------------------
    ! CloudPressure
    !---------------------------

    ! Field name
    name = 'CloudPressure'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
    
    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)
    
    ! Allocate array
    ALLOCATE( CloudPressure( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, CloudPressure )

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !-----------------------------
    ! Cleanup and quit
    !-----------------------------

    ! Detach from swath
    CALL He5SwathDetach( sId )
    
    ! Close HDF-EOS5 file
    CALL He5FileClose( fId )
    
  END SUBROUTINE HCHO_OmiReadCld

!------------------------------------------------------------------------------

  SUBROUTINE HCHO_OmiReadData( fileName, swathName )

    !======================================================================
    ! HCHO_OmiReadData reads data fields from an OMI HDF-EOS5 swath file.
    ! (bmy, 1/20/06)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) fileName  (CHARACTER) : Name of HDF-EOS5 file
    ! (2 ) swathName (CHARACTER) : Name of swath in HDF-EOS5 file:
    ! 
    ! NOTES:
    !======================================================================

    ! Arguments
    CHARACTER(LEN=*) , INTENT(IN) :: fileName, swathName
    
    ! Local variables
    INTEGER                       :: fId,     sId,     as  
    INTEGER                       :: fldRank, fldDims(HE5_MAX_DIMS)
    INTEGER(HE5_INT)              :: nX,      nY
    CHARACTER(LEN=HE5_MAX_CHAR)   :: fldType, name
    CHARACTER(LEN=HE5_MAX_CHAR)   :: fldDimNames(HE5_MAX_DIMS)

    !-----------------------------
    ! HCHO_OmiReadData begins here!
    !-----------------------------

    ! Turn off verbose output
    CALL He5VerboseOutput( .FALSE. )

    ! Open HDF-EOS5 swath and get the file ID number
    CALL He5FileOpen( fileName, fId )

    ! Attach to swath and get swath ID number
    CALL He5SwathAttach( fId, swathName, sId )

    ! Echo info
    WRITE( 6, '(a)' ) REPEAT( '-', 79 )
    WRITE( 6, 100   ) TRIM( fileName )
    WRITE( 6, 110   ) TRIM( swathName )
100 FORMAT( 'HCHO_OmiReadData: Reading ',       a )
110 FORMAT( 'HCHO_OmiReadData: Swath name is ', a )

    !-----------------------------
    ! OmiTime
    !-----------------------------

    ! Field name
    name = 'Time'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
    
    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    
    ! Allocate array
    ALLOCATE( OmiTime( nX ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, OmiTime )

    ! Prevent FPE's when casting to REAL*4
    WHERE ( ABS(OmiTime) .GT. 1.0D30 ) 
       OmiTime = -1.0D30
    ENDWHERE

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !-----------------------------
    ! SolarZenithAngle
    !-----------------------------

    ! Field name
    name = 'SolarZenithAngle'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
    
    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)
    
    ! Allocate array
    ALLOCATE( SolarZenithAngle( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, SolarZenithAngle )

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )


    !-----------------------------
    ! ViewingZenithAngle
    !-----------------------------

    ! Field name
    name = 'ViewingZenithAngle'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
    
    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)
    
    ! Allocate array
    ALLOCATE( ViewingZenithAngle( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, ViewingZenithAngle )

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !-----------------------------
    ! Latitude
    !-----------------------------

    ! Field name
    name = 'Latitude'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
    
    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)
    
    ! Allocate array
    ALLOCATE( Latitude( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )
    Latitude = 0e0

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, Latitude )

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !-----------------------------
    ! Longitude
    !-----------------------------

    ! Field name
    name = 'Longitude'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
    
    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)
    
    ! Allocate array
    ALLOCATE( Longitude( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, Longitude )

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !---------------------------
    ! PixelCornerLatitudes
    !---------------------------

    ! Field name
    name = 'PixelCornerLatitudes'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
    
    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)
    
    ! Allocate array
    ALLOCATE( PixelCornerLatitudes( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, PixelCornerLatitudes )

    !---------------------------
    ! PixelCornerLongitudes
    !---------------------------

    ! Field name
    name = 'PixelCornerLongitudes'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
    
    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)
    
    ! Allocate array
    ALLOCATE( PixelCornerLongitudes( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, PixelCornerLongitudes )

    !---------------------------
    ! ColumnAmount
    !---------------------------

    ! Field name
    name = 'ColumnAmount'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
    
    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)
    
    ! Allocate array
    ALLOCATE( ColumnAmount( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, ColumnAmount )


    ! Prevent FPE's when casting to REAL*4
    WHERE ( ABS(ColumnAmount) .GT. 1.0D30 ) 
       ColumnAmount = -1.0D30
    ENDWHERE

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !---------------------------
    ! ColumnAmountDestriped
    !---------------------------

    ! Field name
    name = 'ColumnAmountDestriped'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
    
    ! Copy dims to saved variables for consistency check later
    OMI_dimX = fldDims(1)
    OMI_dimY = fldDims(2)


    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)
    
    ! Allocate array
    ALLOCATE( ColumnAmountDestriped( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, ColumnAmountDestriped )


    ! Prevent FPE's when casting to REAL*4
    WHERE ( ABS(ColumnAmountDestriped) .GT. 1.0D30 ) 
       ColumnAmountDestriped = -1.0D30
    ENDWHERE

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !---------------------------
    ! FittingRMS
    !---------------------------

    ! Field name
    name = 'FittingRMS'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
    
    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)
    
    ! Allocate array
    ALLOCATE( FittingRMS( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, FittingRMS )

    ! Prevent FPE's when casting to REAL*4
    WHERE ( ABS(FittingRMS) .GT. 1.0D30 ) 
       FittingRMS = -1.0D30
    ENDWHERE

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !---------------------------
    ! ColumnUncertainty
    !---------------------------

    ! Field name
    name = 'ColumnUncertainty'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
    
    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)
    
    ! Allocate array
    ALLOCATE( ColumnUncertainty( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, ColumnUncertainty )

    ! Prevent FPE's when casting to REAL*4
    WHERE ( ABS(ColumnUncertainty) .GT. 1.0D30 ) 
       ColumnUncertainty = -1.0D30
    ENDWHERE

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !---------------------------
    ! MainDataQualityFlag
    !---------------------------

    ! Field name
    name = 'MainDataQualityFlag'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
    
    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)
    
    ! Allocate array
    ALLOCATE( MainDataQualityFlag( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, MainDataQualityFlag )

    !-----------------------------
    ! AirMassFactor
    !-----------------------------

    ! Field name
    name = 'AirMassFactor'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
    
    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)
    
    ! Allocate array
    ALLOCATE( AirMassFactor( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, AirMassFactor )

    ! Prevent FPE's when casting to REAL*4
    WHERE ( ABS(AirMassFactor) .GT. 1.0D30 ) 
       AirMassFactor = -1.0D30
    ENDWHERE

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !-----------------------------
    ! AirMassFactorDiagnosticFlag
    !-----------------------------

    ! Field name
    name = 'AirMassFactorDiagnosticFlag'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )
    
    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)
    
    ! Allocate array
    ALLOCATE( AirMassFactorDiagnosticFlag( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, AirMassFactorDiagnosticFlag )

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )



    !---------------------------
    ! Cleanup and quit
    !---------------------------

    ! Detach from swath
    CALL He5SwathDetach( sId )
    
    ! Close HDF-EOS5 file
    CALL He5FileClose( fId )
    
 
  END SUBROUTINE HCHO_OmiReadData
    

!------------------------------------------------------------------------------

  SUBROUTINE NO2_OmiReadData( fileName, swathName )

    !======================================================================
    ! NO2_OmiReadData reads data fields from an OMI HDF-EOS5 swath file.
    ! Also reads Omi cloud data from this same file.
    !  Modified for NO2 by Gray O'Byrne (Jan, 2008)
    !
    ! Arguments as Input:
    ! ---------------------------------------------------------------------
    ! (1 ) fileName  (CHARACTER) : Name of HDF-EOS5 file
    ! (2 ) swathName (CHARACTER) : Name of swath in HDF-EOS5 file:
    !
    ! NOTES:
    !======================================================================

    ! Arguments
    CHARACTER(LEN=*) , INTENT(IN) :: fileName, swathName

    ! Local variables
    INTEGER                       :: fId,     sId,     as
    INTEGER                       :: fldRank, fldDims(HE5_MAX_DIMS)
    INTEGER(HE5_INT)              :: nX,      nY
    CHARACTER(LEN=HE5_MAX_CHAR)   :: fldType, name
    CHARACTER(LEN=HE5_MAX_CHAR)   :: fldDimNames(HE5_MAX_DIMS)

    !-----------------------------
    ! NO2_OmiReadData begins here!
    !-----------------------------

    ! Turn off verbose output
    CALL He5VerboseOutput( .FALSE. )

    ! Open HDF-EOS5 swath and get the file ID number
    CALL He5FileOpen( fileName, fId )

    ! Attach to swath and get swath ID number
    CALL He5SwathAttach( fId, swathName, sId )

    ! Echo info
    WRITE( 6, '(a)' ) REPEAT( '-', 79 )
    WRITE( 6, 100   ) TRIM( fileName )
    WRITE( 6, 110   ) TRIM( swathName )
100 FORMAT( 'NO2_OmiReadData: Reading ',       a )
110 FORMAT( 'NO2_OmiReadData: Swath name is ', a )

    !-----------------------------
    ! Cloud Data
    !-----------------------------
    ! RelativeAzimuthAngle
    !-----------------------------

    ! Field name
    name = 'ViewingAzimuthAngle'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )

    ! Copy dims to saved variables for consistency check later
    OMI_dimX = fldDims(1)
    OMI_dimY = fldDims(2)

    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)

    ! Allocate array
    ALLOCATE( RelativeAzimuthAngle( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, RelativeAzimuthAngle )

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !-----------------------------
    ! CloudFractionNO2
    !-----------------------------

    ! Field name
    name = 'CloudFraction'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )

    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)

    ! Allocate array
    ALLOCATE( CloudFractionNO2( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, CloudFractionNO2 )

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !---------------------------
    ! CloudPressure
    !---------------------------

    ! Field name
    name = 'CloudPressure'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )

    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)

    ! Allocate array
    ALLOCATE( CloudPressure( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, CloudPressure )

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !-----------------------------
    ! OmiTime
    !-----------------------------

    ! Field name
    name = 'Time'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )

    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)

    ! Allocate array
    ALLOCATE( OmiTime( nX ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, OmiTime )

    ! Prevent FPE's when casting to REAL*4
    WHERE ( ABS(OmiTime) .GT. 1.0D30 )
       OmiTime = -1.0D30
    ENDWHERE

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !-----------------------------
    ! SolarZenithAngle
    !-----------------------------

    ! Field name
    name = 'SolarZenithAngle'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )

    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)

    ! Allocate array
    ALLOCATE( SolarZenithAngle( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, SolarZenithAngle )

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )


    !-----------------------------
    ! ViewingZenithAngle
    !-----------------------------

    ! Field name
    name = 'ViewingZenithAngle'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )

    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)

    ! Allocate array
    ALLOCATE( ViewingZenithAngle( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, ViewingZenithAngle )

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !-----------------------------
    ! Latitude
    !-----------------------------

    ! Field name
    name = 'Latitude'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )

    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)

    ! Allocate array
    ALLOCATE( Latitude( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )
    Latitude = 0e0

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, Latitude )

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !-----------------------------
    ! Longitude
    !-----------------------------

    ! Field name
    name = 'Longitude'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )

    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)

    ! Allocate array
    ALLOCATE( Longitude( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, Longitude )

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !---------------------------
    ! ColumnAmountNO2 Trop
    !---------------------------

    ! Field name
    name = 'ColumnAmountNO2Trop'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )

    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)

    ! Allocate array
    ALLOCATE( ColumnAmountNO2Trop( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, ColumnAmountNO2Trop )

    ! Prevent FPE's when casting to REAL*4
    WHERE ( ABS(ColumnAmountNO2Trop) .GT. 1.0D30 )
       ColumnAmountNO2Trop = -1.0D30
    ENDWHERE

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !---------------------------
    ! ColumnAmountNO2 Below Cloud
    !---------------------------

    ! Field name
    name = 'ColumnAmountNO2BelowCloud'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )

    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)

    ! Allocate array
    ALLOCATE( ColumnAmountNO2BelowCloudTrop( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData(sId,name,nX,nY,ColumnAmountNO2BelowCloudTrop)

    ! Prevent FPE's when casting to REAL*4
    WHERE ( ABS(ColumnAmountNO2BelowCloudTrop) .GT. 1.0D30 )
       ColumnAmountNO2BelowCloudTrop = -1.0D30
    ENDWHERE

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )


    !---------------------------
    ! ColumnUncertaintyNO2
    !---------------------------

    ! Field name
    name = 'ColumnAmountNO2TropStd'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )

    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)

    ! Allocate array
    ALLOCATE( ColumnUncertaintyNO2Trop( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, ColumnUncertaintyNO2Trop )

    ! Prevent FPE's when casting to REAL*4
    WHERE ( ABS(ColumnUncertaintyNO2Trop) .GT. 1.0D30 )
       ColumnUncertaintyNO2Trop = -1.0D30
    ENDWHERE

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !---------------------------
    ! ColumnUncertaintyNO2 Below Cloud
    !---------------------------

    ! Field name
    name = 'ColumnAmountNO2BelowCloudStd'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )

    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)

    ! Allocate array
    ALLOCATE( ColumnUncertaintyNO2BelowCloudTrop( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData(sId,name,nX,nY,ColumnUncertaintyNO2BelowCloudTrop)

    ! Prevent FPE's when casting to REAL*4
    WHERE ( ABS(ColumnUncertaintyNO2BelowCloudTrop) .GT. 1.0D30 )
       ColumnUncertaintyNO2BelowCloudTrop = -1.0D30
    ENDWHERE

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )


    !---------------------------
    ! MainDataQualityFlag
    !---------------------------

    ! Field name
    name = 'vcdQualityFlags'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )

    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)

    ! Allocate array
    ALLOCATE( MainDataQualityFlag( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, MainDataQualityFlag )

    !-----------------------------
    ! AMFUnpolluted
    !-----------------------------

    ! Field name
    name = 'AMFUnpolluted'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )

    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)

    ! Allocate array
    ALLOCATE( AMFUnpolluted( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, AMFUnpolluted )

    ! Prevent FPE's when casting to REAL*4
    WHERE ( ABS(AMFUnpolluted) .GT. 1.0D30 )
       AMFUnpolluted = -1.0D30
    ENDWHERE

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !-----------------------------
    ! AMFPolluted
    !-----------------------------

    ! Field name
    name = 'AMFPolluted'

    ! Get field info (array dimensions are in fldDims)
    CALL He5SwathFldInfo( sId, name, fldType, fldRank, fldDims, fldDimNames )

    ! Copy dims to INTEGER*8 variables BEFORE array allocation!
    nX = fldDims(1)
    nY = fldDims(2)

    ! Allocate array
    ALLOCATE( AMFPolluted( nX, nY ), stat=as )
    IF ( as /= 0 ) CALL He5AllocErr( name )

    ! Read data
    CALL He5SwathReadData( sId, name, nX, nY, AMFPolluted )

    ! Prevent FPE's when casting to REAL*4
    WHERE ( ABS(AMFPolluted) .GT. 1.0D30 )
       AMFPolluted = -1.0D30
    ENDWHERE

    ! Get fill value (if necessary to strip out missing data values)
    !CALL He5SwathFillValue( sId, name, dataFill )

    !---------------------------
    ! Cleanup and quit
    !---------------------------

    ! Detach from swath
    CALL He5SwathDetach( sId )

    ! Close HDF-EOS5 file
    CALL He5FileClose( fId )

  END SUBROUTINE NO2_OmiReadData

!------------------------------------------------------------------------------

  SUBROUTINE OmiCleanup
    
  !========================================================================
  ! Subroutine OmiCleanup deallocates all module arrays (bmy, 1/20/06)
  !========================================================================
  IF (ALLOCATED(RelativeAzimuthAngle    )) DEALLOCATE(RelativeAzimuthAngle    )
  IF (ALLOCATED(CloudFractionNO2        )) DEALLOCATE(CloudFractionNO2        )
  IF (ALLOCATED(CloudFraction           )) DEALLOCATE(CloudFraction           )
  IF (ALLOCATED(CloudPressure           )) DEALLOCATE(CloudPressure           )
  IF (ALLOCATED(OmiTime                 )) DEALLOCATE(OmiTime                 )
  IF (ALLOCATED(SolarZenithAngle        )) DEALLOCATE(SolarZenithAngle        )
  IF (ALLOCATED(ViewingZenithAngle      )) DEALLOCATE(ViewingZenithAngle      )
  IF (ALLOCATED(Latitude                )) DEALLOCATE(Latitude                )
  IF (ALLOCATED(Longitude               )) DEALLOCATE(Longitude               )
  IF (ALLOCATED(PixelCornerLatitudes    )) DEALLOCATE(PixelCornerLatitudes    )
  IF (ALLOCATED(PixelCornerLongitudes   )) DEALLOCATE(PixelCornerLongitudes   )
  IF (ALLOCATED(ColumnAmount            )) DEALLOCATE(ColumnAmount            )
  IF (ALLOCATED(ColumnAmountNO2Trop     )) DEALLOCATE(ColumnAmountNO2Trop     )
  IF (ALLOCATED(ColumnAmountNO2BelowCloudTrop)) THEN
     DEALLOCATE(ColumnAmountNO2BelowCloudTrop )
  ENDIF
  
  IF (ALLOCATED(ColumnAmountDestriped   )) DEALLOCATE(ColumnAmountDestriped   )
  IF (ALLOCATED(FittingRMS              )) DEALLOCATE(FittingRMS              )
  IF (ALLOCATED(ColumnUncertainty       )) DEALLOCATE(ColumnUncertainty       )
  IF (ALLOCATED(ColumnUncertaintyNO2Trop)) DEALLOCATE(ColumnUncertaintyNO2Trop)
  IF (ALLOCATED(ColumnUncertaintyNO2BelowCloudTrop)) THEN
     DEALLOCATE(ColumnUncertaintyNO2BelowCloudTrop )
  ENDIF
  IF (ALLOCATED(MainDataQualityFlag     )) DEALLOCATE(MainDataQualityFlag     )
  IF (ALLOCATED(AirMassFactor           )) DEALLOCATE(AirMassFactor           )
  IF (ALLOCATED(AMFPolluted             )) DEALLOCATE(AMFPolluted             )
  IF (ALLOCATED(AMFUnpolluted           )) DEALLOCATE(AMFUnpolluted           )
  IF ( ALLOCATED( AirMassFactorDiagnosticFlag ) ) THEN
     DEALLOCATE( AirMassFactorDiagnosticFlag )
  ENDIF

  END SUBROUTINE OmiCleanup

!------------------------------------------------------------------------------

END MODULE OmiModule









