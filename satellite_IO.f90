module satellite_IO

      ! C Preprocessor #define statements for conditional compilation
#     include "define.h"

#if defined ( OMI )
   !The if statement is to avoid HDF libraries when not needed.
   !LSURL 2015-11-24 Because I've hacked a way to avoid HDF for OMI
   !I have commented out this line
   !USE OmiModule
#endif

   USE ParameterModule,  ONLY: NVALMAX, LUN

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: input_GOME,input_newGOME,input_SCIA,input_OMI

   INTEGER:: LINE
   INTEGER:: DAYNO

!-------------------------------------------------------------------------
! Satellite information not used for AMF calculation
!-------------------------------------------------------------------------

! Stellite instrument
!   INTEGER:: ULN(NVALMAX)                !unique line number
   INTEGER:: PIXEL(NVALMAX)                ! Pixel number
   INTEGER:: SCAN(NVALMAX)                 ! Scan position
   !REAL:: SC(NVALMAX)                    ! Slant column (mol/cm2)
   REAL:: DSC(NVALMAX)                   ! Fitting uncertainty (mol/cm2)
   REAL:: FITRMS(NVALMAX)                ! Fitting rms
   REAL:: DFCLD(NVALMAX)                 ! Cloud fraction error (%)
   REAL:: DPCLD(NVALMAX)                 ! Cloud pressure error (hPa)
   REAL:: I1(NVALMAX), J1(NVALMAX)       ! pixel edges
   REAL:: I2(NVALMAX), J2(NVALMAX) 
   REAL:: I3(NVALMAX), J3(NVALMAX)
   REAL:: I4(NVALMAX), J4(NVALMAX)
 
   CHARACTER(LEN=11), DIMENSION(NVALMAX)  :: TIME 


   CONTAINS
!=========================================================================
!     Satellite Input Routines
!=========================================================================
 
!--------------------------------
!     GOME
!--------------------------------

   subroutine input_GOME(satellite_file,TG_TYPE,LAT_centers,LON_centers,SZA,SVA,AZM,PCLD,FCLD,TCOT,SC,MAXLINE,outPrefix,outSuffix,flag)
 
      ! Inputs
      !---------------------------------------------------------------------------------
      CHARACTER(LEN=255),intent(in)  :: satellite_file    ! The file to extract data from
      INTEGER, intent(in)            :: TG_TYPE           ! Trace Gas, 0=HCHO, 1=NO2, 2=SO2

      ! Outputs
      !---------------------------------------------------------------------------------
      !Variables Needed for AMF calculation
      REAL, intent(out):: SZA(NVALMAX)                    ! Solar Zenith Angle
      REAL, intent(out):: SVA(NVALMAX)                    ! Satellite viewing angle
      REAL, intent(out):: AZM(NVALMAX)                    ! Relative azimuth angle, Not currently used
      REAL, intent(out):: LAT_centers(NVALMAX)            ! Latitude of pixel center
      REAL, intent(out):: LON_centers(NVALMAX)            ! Longitude of pixel center
      REAL, intent(out):: PCLD(NVALMAX)                   ! Cloud Pressure (hPa)
      REAL, intent(out):: FCLD(NVALMAX)                   ! Cloud Fraction (between 0 and 1)
      
      !observation
      REAL, intent(out):: SC(NVALMAX)  !observation mol/cm2

      !Specific outputs for this satellite
      REAL,intent(out):: TCOT(NVALMAX)                    ! Cloud optical thickness

      !Other
      INTEGER, intent(out) :: MAXLINE                       ! The number of observations read in
      CHARACTER(LEN=200),intent(out) :: outSuffix(NVALMAX)  ! A string to pass on data from satellite file to AMF output
      CHARACTER(LEN=100),intent(out) :: outPrefix(NVALMAX)  ! A string to pass on data from satellite file to AMF output
      REAL:: DATE
      LOGICAL, intent(out) :: flag(NVALMAX)       !Satellite flag, do not use measurement if true


      ! Start input_GOME
      !--------------------------------------------------------------------------------

      AZM(:) = 0      !Relative azimuth angle set to zero

      ! Read in slant columns
      OPEN ( FILE = satellite_file,                       &
        UNIT = LUN,                                       &
        FORM = 'FORMATTED',                               &
        STATUS = 'OLD')                                   

      DO LINE = 1, NVALMAX
         READ(LUN, '(i5,i2,15f7.2,3e11.3,i3)', END=900)   &
           PIXEL(LINE), SCAN(LINE),                       &
           SZA(LINE), SVA(LINE),                          &
           I1(LINE), J1(LINE),                            &
           I2(LINE), J2(LINE),                            &
           I3(LINE), J3(LINE),                            &
           I4(LINE), J4(LINE),                            &
           LAT_centers(LINE), LON_centers(LINE),          &
           FCLD(LINE), PCLD(LINE),                        &
           TCOT(LINE),                                    &
           SC(LINE), DSC(LINE), FITRMS(LINE),             &
           DATE
      ENDDO

900     CONTINUE

      TCOT(:)=TCOT(:)*0.25   !Scale optical thickness values to be comparable to ISCCP

      MAXLINE = LINE-1
      CLOSE ( LUN )

      outSuffix(:) = ' suffix'
      outPrefix(:) = 'prefix '
      flag(:)=.false.

   end subroutine input_GOME


!------------------------------
!  NEW GOME
!------------------------------

   subroutine input_newGOME(satellite_file,TG_TYPE,ULN,LAT_centers,LON_centers,SZA,SVA,PCLD,FCLD,MAXLINE,outPrefix,outSuffix,flag)

      ! Inputs
      !---------------------------------------------------------------------------------
      CHARACTER(LEN=255),intent(in)  :: satellite_file    ! The file to extract data from
      INTEGER, intent(in)            :: TG_TYPE           ! Trace Gas, 0=HCHO, 1=NO2, 2=SO2


      ! Outputs
      !---------------------------------------------------------------------------------
      !Variables Needed for all AMF calculations
      REAL, intent(out):: SZA(NVALMAX)                    ! Solar Zenith Angle
      REAL, intent(out):: SVA(NVALMAX)                    ! Satellite viewing angle
      !REAL, intent(out):: AZM(NVALMAX)                    ! Relative azimuth angle, Not currently used
      REAL, intent(out):: LAT_centers(NVALMAX)            ! Latitude of pixel center
      REAL, intent(out):: LON_centers(NVALMAX)            ! Longitude of pixel center
      REAL, intent(out):: PCLD(NVALMAX)                   ! Cloud Pressure (hPa)
      REAL, intent(out):: FCLD(NVALMAX)                   ! Cloud Fraction (between 0 and 1)

      !observation
      !REAL, intent(out):: SC(NVALMAX)  !observation mol/cm2

      !Specific outputs for this satellite
      !REAL,intent(out):: TCOT(NVALMAX)                    ! Cloud optical thickness

      !Other
      INTEGER, intent(out) :: ULN(NVALMAX)                  ! The unique line number
      INTEGER, intent(out) :: MAXLINE                       ! The number of observations read in
      CHARACTER(LEN=200),intent(out) :: outSuffix(NVALMAX)  ! A string to pass on data from satellite file to AMF output
      CHARACTER(LEN=100),intent(out) :: outPrefix(NVALMAX)  ! A string to pass on data from satellite file to AMF output
      LOGICAL, intent(out) :: flag(NVALMAX)       !Satellite flag, do not use measurement if true

      ! Start input_newGOME
      !--------------------------------------------------------------------------------

      !AZM(:) = 0      !Relative azimuth angle set to zero

      ! Read in slant columns
      OPEN ( FILE = satellite_file,                       &
        UNIT = LUN,                                       &
        FORM = 'FORMATTED',                               &
        STATUS = 'OLD')

      DO LINE = 1, NVALMAX
         READ(LUN, *, END=901) &
           ULN(LINE),SCAN(LINE), PIXEL(LINE),             &
           LAT_centers(LINE), LON_centers(LINE),        &
           SZA(LINE), SVA(LINE), FCLD(LINE),PCLD(LINE)                                                          
      ENDDO

 901  CONTINUE

      !TCOT(:)=TCOT(:)*0.25   !Scale optical thickness values to be comparable to ISCCP

      MAXLINE = LINE-1
      CLOSE ( LUN )

      outSuffix(:) = ' suffix'
      outPrefix(:) = 'prefix '
      flag(:)=.false.

   end subroutine input_newGOME


!---------------------------------
! SCIAMACHY
!---------------------------------

   subroutine input_SCIA(satellite_file,TG_TYPE,FRESCOv5,LAT_centers,LON_centers,SZA,SVA,AZM,PCLD,FCLD,FPS,SC,MAXLINE,outPrefix,outSuffix,flag)

      ! Inputs
      !---------------------------------------------------------------------------------
      CHARACTER(LEN=255),intent(in)  :: satellite_file    ! The file to extract data from
      INTEGER, intent(in) :: FRESCOv5 
      INTEGER, intent(in) :: TG_TYPE           ! Trace Gas, 0=HCHO, 1=NO2, 2=SO2

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
      !observation
      REAL, intent(out):: SC(NVALMAX)  !observation mol/cm2      
      !Specific outputs for this satellite
      REAL, intent(out):: FPS(NVALMAX)                    !Overwrite surface pressure with new fresco product surface pressure

      !Other
      INTEGER, intent(out) :: MAXLINE			    ! The number of observations read in
      CHARACTER(LEN=200),intent(out) :: outSuffix(NVALMAX)  ! A string to pass on data from satellite file to AMF output
      CHARACTER(LEN=100),intent(out) :: outPrefix(NVALMAX)  ! A string to pass on data from satellite file to AMF output
      LOGICAL, intent(out) :: flag(NVALMAX)       !Satellite flag, do not use measurement if true

      ! Start input_SCIA
      !--------------------------------------------------------------------------------

      ! Read in slant columns
      OPEN ( FILE = satellite_file,                       &
        UNIT = LUN,                                       &
        FORM = 'FORMATTED',                               &
        STATUS = 'OLD')
      print*,"IO debug point A"
      call flush(5)
      DO LINE = 1, NVALMAX
        IF (FRESCOv5 .EQ. 1) THEN    !FRESCO v5 available
           READ(LUN, '(i5,i3,15f8.3,3f9.3,3e11.3,a11)', END=902) &
             PIXEL(LINE), SCAN(LINE), 		                 &
             SZA(LINE), SVA(LINE),                        &
             AZM(LINE),                                   &
             I1(LINE), J1(LINE),                          &
             I2(LINE), J2(LINE),                          &
             I3(LINE), J3(LINE),                          &
             I4(LINE), J4(LINE),                          &
             LAT_centers(LINE), LON_centers(LINE),        &
             FCLD(LINE), DFCLD(LINE), PCLD(LINE),         &
             DPCLD(LINE), FPS(LINE), SC(LINE), DSC(LINE), FITRMS(LINE), &
             TIME(LINE)                                   
        ELSE             !older version of cloud product
          READ(LUN, '(i5,i3,14f8.3,f9.3,3e11.3,a11)', END=902) &
            PIXEL(LINE), SCAN(LINE),                      &
            SZA(LINE), SVA(LINE),                         &
            AZM(LINE),                                    &
            I1(LINE), J1(LINE),                           &
            I2(LINE), J2(LINE),                           &
            I3(LINE), J3(LINE),                           &
            I4(LINE), J4(LINE),                           &
            LAT_centers(LINE), LON_centers(LINE),         &
            FCLD(LINE), PCLD(LINE),                       &
            SC(LINE), DSC(LINE), FITRMS(LINE),            &
            TIME(LINE)
        ENDIF

        write(outPrefix(LINE),'(2i6)') PIXEL(LINE), SCAN(LINE)
      ENDDO

 902     CONTINUE

      MAXLINE = LINE-1
      
      CLOSE ( LUN )

      outSuffix(:) = ' suffix'
      flag(:)=.false.

   end subroutine input_SCIA


!---------------------------------------------
! OMI
!---------------------------------------------

   subroutine input_OMI(satellite_file,TG_TYPE,CLDFILE,LAT_centers,LON_centers,SZA,SVA,AZM,PCLD,FCLD,SC,MAXLINE,outPrefix,outSuffix,flag)

      !The following require HDF libraries. Comment out if not needed.

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
      !observation
      REAL, intent(out):: SC(NVALMAX)  !observation mol/cm2
      !Specific outputs for this satellite

      !Other
      INTEGER, intent(out) :: MAXLINE                       ! The number of observations read in
      CHARACTER(LEN=200),intent(out) :: outSuffix(NVALMAX)  ! A string to pass on data from satellite file to AMF output
      CHARACTER(LEN=100),intent(out) :: outPrefix(NVALMAX)  ! A string to pass on data from satellite file to AMF output
      LOGICAL, intent(out) :: flag(NVALMAX)       !Satellite flag, do not use measurement if true

      ! Start input_OMI
      !--------------------------------------------------------------------------------



! Read OMI cloud and Trace Gas Data files

       !Data read function in OmiModule.f90, hopefuly this will ease optional compilation
       ! with HDF libraries in future versions.

#if defined ( OMI )
! LSURL 2015-11-24 commenting out here as we've hacked an HDF-free solution
! to OMI
!       CALL OMI_GetAmfInput(satellite_file,TG_TYPE,CLDFILE,LAT_centers,LON_centers,SZA,SVA,AZM,PCLD,FCLD,MAXLINE,outPrefix,outSuffix)
#endif
       flag(:)=.false. 

   end subroutine input_OMI


end module satellite_IO
