      SUBROUTINE LUT_USE(RSPIN,RSAIN,
     &                   SZA,LOSZA,RELAZM,
     &                   AMFgeo,TRPAUSE,PCLD,
     &                   KHAT,KHATCLD,
     &                   ACLR,ACLD)

C------------------------------------------------------------------
C Sets up and calls look up table
C
C------------------------------------------------------------------

C  include files
      USE MAmfLut,      ONLY : ReadAmfLut, FitWindowCentre, GetAmf
      USE Mrweight,         ONLY : rweight
      USE CMN_SIZE_MOD		!Size Parameters

      IMPLICIT NONE


C  Define Variables

C  Inputs
      REAL*4,intent(in) :: RSPIN        !Surface Pressure
      REAL*4,intent(in) :: RSAIN        !Surface Albedo
      REAL*4,intent(in) :: PCLD        !Cloud Pressure
      REAL*4,intent(in) :: SZA          !Solar Zenith Angle
      REAL*4,intent(in) :: LOSZA        !Satellite Viewing Angle
      REAL*4,intent(in) :: RELAZM       !Relative Azimuth Angle

      REAL*4,intent(in) :: AMFgeo       !Geometric Air Mass Factor
      REAL*8,intent(in) :: TRPAUSE      !Model Tropopause Level

C  Outputs
      REAL*8,intent(out) :: KHAT(LLPAR), KHATCLD(LLPAR)    !Clear and Cloudy scattering weights
      REAL*8,intent(out) :: ACLD, ACLR                     !Cloudy and Clear Sky Albedos

C  Local Variables
      integer I, LEV      !Loop Indices

      REAL*4 RACLR, RACLD  !Clear and Cloud Sky Albedos as REAL*4

      REAL*4 TopPress, BottomPress, MidPress  !Model layer pressures
      REAL*4 AboveCloudFraction               !Fraction of model layer above cloud

      REAL*4 AMFClearRel, AMFCloudRel         ! = KHAT/AMFgeo 
      

!------------------------------------------------------------------
! Estimate radiation intensities for cloud and no cloud conditions
!   This is during AMF calculation when using Lidort
!------------------------------------------------------------------
               ! Cloudy
               CALL rweight( FitWindowCentre, PCLD,
     &              0.80, SZA, LOSZA,
     &              RELAZM, RACLD )

C               write(*,*) PCLD, SZA, LOSZA, RELAZM, RACLD

               ! Clear
               CALL rweight( FitWindowCentre,RSPIN,
     &              RSAIN, SZA, LOSZA,
     &              RELAZM, RACLR )

               ACLD=RACLD
               ACLR=RACLR

               IF ( ACLR .LE. 0.0 .OR. ACLD .LE. 0.0 ) THEN
                  PRINT*,'ERROR ass_Hx.f90: ACLR, ACLD = ',
     &                 ACLR, ACLD
                  STOP
               END IF



               DO LEV = 1, TRPAUSE-1

                  ! ====================
                  ! Get layer pressure
                  ! ====================
#if defined ( GEOS_4 ) || defined ( GEOS_5 )
                  BottomPress = ETA_A(LEV) +
     &                 ( ETA_B(LEV) * RSPIN )
                  TopPress    = ETA_A(LEV+1) +
     &                 ( ETA_B(LEV+1) * RSPIN )
                  MidPress    = 0.50 * ( BottomPress + TopPress )
#else
                  BottomPress = MSIGMAE(LEV) * ( RSPIN - 0.01 ) + 0.01
                  TopPress    = MSIGMAE(LEV+1) * ( RSPIN - 0.01 ) + 0.01
                  MidPress    = MSIGMA(LEV) * ( RSPIN - 0.01 ) + 0.01
#endif

                  ! ====================
                  ! Retrieve scattering weights
                  ! from lookup table
                  ! ====================

                  ! Clear sky case
                  CALL GetAmf( MidPress,RELAZM,LOSZA,
     &                 SZA,RSAIN, RSPIN,
     &                 AMFClearRel )

                  KHAT(LEV) = AMFClearRel * AMFgeo

                  ! Cloudy case
                  IF ( PCLD .GE. TopPress - 1.0E-5 ) THEN

                     IF ( PCLD .GE. BottomPress - 1.0E-5 ) THEN

                        ! Model level above the cloud top
                        CALL GetAmf( MidPress,RELAZM,LOSZA,
     &                       SZA,0.80,PCLD,
     &                       AMFCloudRel )

                        KHATCLD(LEV) = AMFCloudRel * AMFgeo

                     ELSE

                        ! Model level contains cloud top
                        AboveCloudFraction = ( PCLD - TopPress )
     &                       / (BottomPress-TopPress)

                        CALL GetAmf( MidPress,RELAZM, LOSZA,
     &                       SZA,0.80,MidPress + 1.0E-5,
     &                       AMFCloudRel )

                        KHATCLD(LEV) = AMFCloudRel * AMFgeo *
     &                       AboveCloudFraction

                     END IF
                  ELSE

                     ! model level below cloud top
                     KHATCLD(LEV) = 0.0
                  END IF

               ENDDO

C   Finish
      end
