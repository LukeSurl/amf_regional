      PROGRAM CONVERT2VERTICAL

C-------------------------------------------------------------------------
C This code will read the slant column and cloud observations, 
C find the corresponding 3-d 
C chemistry trace gas profile and surface pressure, read in the
C surface albedo information and use this information to compute 
C the air mass factor (AMF) used to convert the Satellite 
C slant columns to vertical columns.
C-------------------------------------------------------------------------

      ! C Preprocessor #define statements for conditional compilation
#     include "define.h"

      ! References to F90 modules
      USE BPCH2_MOD
      USE ParameterModule
      USE CMN_SIZE_MOD
      USE satellite_IO

#if defined(USELUT)
      USE MAmfLut
#endif

      IMPLICIT NONE

C  LIDORT input variables
#     include "include_s/LIDORT.PARS"
#     include "include_e/LIDORT_L.PARS"
#     include "include_s/LIDORT_CONTROL.VARS"
#     include "include_s/LIDORT_MODEL.VARS"
#     include "include_s/LIDORT_GEOPHYS.VARS"
#     include "include_e/LIDORT_L_CONTROL.VARS"


C=========================================================================
C
C   Variable Declarations
C
C=========================================================================

C-------------------------------------------------------------------------
C Inputs from run script
C-------------------------------------------------------------------------
      CHARACTER(LEN=255) :: SATFILE       ! File containing satellite geometry
      CHARACTER(LEN=255) :: CLDFILE       ! File containing cloud data
      INTEGER DAYNO                       ! Day of Month
      CHARACTER(LEN=255) :: TG_FILE       ! File containing trace gas profile (GEOS-Chem)
      CHARACTER(LEN=255) :: AEROSOL_FILE  ! File containing aerosol profile (GEOS-Chem)
      CHARACTER(LEN=255) :: TROPO_FILE    ! File containing tropopause heights (GEOS-Chem)
      CHARACTER(LEN=255) :: OUTFILE       ! Output file for AMFs
      CHARACTER(LEN=255) :: OUTFILE_WSCAT ! Output file for scattering weights
      CHARACTER(LEN=255) :: OUTFILE_Shape ! Output file for profile shapes
      CHARACTER(LEN=4  ) :: YR_C !Year as a character
      CHARACTER(LEN=2  ) :: MONTH_C !Month as a character
      CHARACTER(LEN=2  ) :: DAYNO_C !Day as a character
      
      
      
      
      INTEGER YR
      INTEGER MONTH
      INTEGER TG_TYPE           ! 0 = HCHO, 1 = NO2
      INTEGER FRESCOv5          ! 1 if True

C------------------------------------------------------------------------
C Loop Indices
C------------------------------------------------------------------------
      INTEGER I, J, L, J00, N
      INTEGER LINE
      INTEGER LEV

C-------------------------------------------------------------------------
C Binary punch file stuff
C-------------------------------------------------------------------------
      REAL*8             :: TAU_IN       ! Time stamp for GEOS-Chem files
      INTEGER            :: TRACER_IN    ! Tracer number in diag51.f (GEOS-Chem)
      CHARACTER(LEN=40)  :: CATEGORY_IN  ! (GEOS-Chem)

      INTEGER :: IND(NAER)  !Tracer numbers for each aerosol type (NAER set=5 in CMN_SIZE)
      DATA IND /4, 6, 9, 12, 15/  ! (4=dust,6=SO4,9=BC,12=OC,15=SSa,18=SSc)
      REAL*8  :: TG(IIPAR,JJPAR,LLPAR)         !Trace gas mixing ratio
      REAL*8  :: ARSL(IIPAR,JJPAR,LLPAR,NAER)  ! Aerosol optical depth
      REAL*8  :: ALBD(IIPAR_ALB,JJPAR_ALB)     !Surface reflectivity
      REAL*8  :: PS(IIPAR,JJPAR)               !Surface pressure
      REAL*8  :: TMP(IIPAR,JJPAR,LLPAR)        !Temperature
      REAL*8  :: FRCLND(IIPAR,JJPAR)           !Fraction of grid occupied by land
      REAL*8  :: TRPAUSE(IIPAR,JJPAR)          !Tropopause level
      REAL*8  :: GEOS_VC(NVALMAX) !Vertical Column in the GEOS-CHEM model

C Arrays to store temporary values read from BPCH files
      REAL*4  :: ARRAY(IIPAR,JJPAR,LLPAR) 
      REAL*4  :: ARRAY_ALB(IIPAR_ALB,JJPAR_ALB)


C-------------------------------------------------------------------------
C Observations
C-------------------------------------------------------------------------

C Satellite instrument
      INTEGER ULN(NVALMAX)                  ! Unique line number
      REAL*4 SZA(NVALMAX)                   ! Solar Zenith Angle
      REAL*4 LOSZA(NVALMAX)                 ! Satellite viewing angle
      REAL*4 RELAZM(NVALMAX)                ! Relative azimuth angle
      REAL*4 FCLD(NVALMAX)                  ! Cloud fraction (%)
      REAL*4 PCLD(NVALMAX)                  ! Cloud top pressure (hPa)
      REAL*4 OLATS(NVALMAX)                 ! Latitude of pixel center
      REAL*4 OLONS(NVALMAX)                 ! Longitude of pixel center
      REAL*4 TCOT(NVALMAX)                  ! Cloud optical thickness
      REAL*8 RSZA, RVZA, RCOT, RPCLD, RFCLD
      LOGICAL LAMBERTIAN_CLOUDS
      REAL FPS(NVALMAX)                  !Fresco Surface Pressure
      LOGICAL FLAG(NVALMAX)		 !Satellite flag, do not use measurement

      !observation
      REAL*4 SC(NVALMAX)  !observation mol/cm2

      ! To pass along extra satellite information to output:
      CHARACTER(LEN=200):: outSuffix(NVALMAX)  ! A string to pass on data from satellite file to AMF output
      CHARACTER(LEN=100):: outPrefix(NVALMAX)  ! A string to pass on data from satellite file to AMF output

      INTEGER MAXLINE		!Number of measurements read in


C Climatologies
      REAL*8 SPIN(NVALMAX)      ! Surface pressure
      REAL*8 SAIN(NVALMAX)      ! Surface albedo
      REAL*4 RSPIN, RSAIN
      REAL*8 AIRMIX(LLPAR)      ! Equal to 1 everywhere


C-------------------------------------------------------------------------
C AMF database
C-------------------------------------------------------------------------

      INTEGER IUSE, JUSE, IUSE_ALB, JUSE_ALB
      REAL*8 KHAT(NVALMAX,LLPAR)        ! Clear sky scattering weights
      REAL*8 ACLR		! Clear sky albedo
      REAL*8 ACLD		! Cloudy sky albedo
      REAL*8 KHATCLD(LLPAR)     ! Cloudy sky scattering weights
      REAL*8 S(NVALMAX,LLPAR)      			     	! Shape factors
      REAL*8 AIRSIGMA, TGSIGMA, TMPSIGMA, WEIGHT
      REAL*8 AMF(NVALMAX), AMFCLD(NVALMAX), AMFCLR(NVALMAX)	!Air Mass Factors
      REAL*8 TRAD, WCLR(NVALMAX), WCLD(NVALMAX)			!Clear vs Cloudy Weights
      REAL*8 DSIGMA(LLPAR)
      
      REAL*8 SCATCLD (LLPAR)            ! Cloud sky scattering weights
      REAL*8 SCATCLR (LLPAR)            ! Clear sky scattering weights
      REAL*8 WSCAT(NVALMAX,LLPAR)       ! Total Scattering Weight

      REAL*4 AMFgeo(NVALMAX)     !Geometric Air Max Factor
      REAL*4 htopatm, eps, xh    !Use to calculate AMFgeo

C------------------------------------------------------------------------
C Model Inputs
C------------------------------------------------------------------------

      REAL*8 TGMIX(LLPAR)       ! Trace gas mixing ratio
      REAL*8 AERPROF(LLPAR,NAER)  ! Aerosol profile


C-----------------------------------------
C   Some Physical constants
C-----------------------------------------
      REAL,PARAMETER    :: MA     =0.02897     ! Molecular weight of air kg/mol.
      REAL,PARAMETER    :: AVGRAD =6.022D23      ! Avagadro's number
      REAL,PARAMETER    :: G      =9.8065
      REAL,PARAMETER    :: EarthRadius = 6.3785e6
      REAL,PARAMETER    :: pi          = 3.1415927

C------------------------------------------
C Other
C------------------------------------------
      INTEGER STATUS_INPUTREAD      !Check lidort input success
      INTEGER SLEN                  !Length of filename
      LOGICAL :: FIRSTCALL=.TRUE.   !Signal for first call of AMF LUT subroutines


C=========================================================================
C
C  Begin Data Input
C
C=========================================================================
      print*,"It begins..."
      call flush(5)

C-------------------------------------------------------------------------
C User input 
C-------------------------------------------------------------------------
      !READ ( 5, '(A)' ) SATFILE  ! File containing slant columns
      call get_command_argument(1, SATFILE)
      
!      READ ( 5, '(A)' ) CLDFILE   ! File containing cloud data

!      READ ( 5, '(A)' ) TG_FILE  ! File containing trace gas profiles
      call get_command_argument(2, TG_FILE)
      
!      READ ( 5, '(A)' ) AEROSOL_FILE   ! File containing aerosol profiles

!      READ ( 5, '(A)' ) TROPO_FILE ! File containing tropopause height

 !     READ ( 5, '(A)' ) OUTFILE   ! Output file
      call get_command_argument(3, OUTFILE)
      
!      READ ( 5, '(A)' ) OUTFILE_WSCAT !Output file for scattering weights

!      READ ( 5, '(A)' ) OUTFILE_Shape !Output file for profile shapes

!      READ ( 5, *     ) YR
      call get_command_argument(4, YR_C)
      read(YR_C,'(i4)') YR !convert to int
      
!      READ ( 5, *     ) MONTH
      call get_command_argument(5, MONTH_C)
      read(MONTH_C,'(i2)') MONTH !convert to int
      
!      READ ( 5, *     ) DAYNO     ! Day of month
      call get_command_argument(6, DAYNO_C)
      read(DAYNO_C,'(i2)') DAYNO !convert to int      
      
      !READ ( 5, *     ) TG_TYPE
      
      !READ ( 5, *     ) FRESCOv5
      
      call get_command_argument(7, TG_TYPE_C_
      read(TG_TYPE_C,'(i2)') TG_FILE
      
      print*,"Input read..."
      call flush(5)

C     DUPLICATION OF INPUT - LUKE SURL 2016-08-24
      AEROSOL_FILE = TG_FILE
      TROPO_FILE   = TG_FILE
      
      !WE'll never change these, so just hard-code
      TG_TYPE = 0
      FRESCOv5 = 0 


#if defined ( USELUT )
!-------------------------------------------------------------------------
! On start: read air-mass factor lookup table
!-------------------------------------------------------------------------

      IF ( FIRSTCALL ) THEN
         CALL ReadAmfLut
         FIRSTCALL = .FALSE.
      ENDIF
#endif

C-------------------------------------------------------------------------
C Satellite-specific settings
C-------------------------------------------------------------------------

#if defined ( GOME )
      YR = YR + 1990
#endif

#if defined ( GOME ) || defined ( NEW_GOME )
      LAMBERTIAN_CLOUDS = .FALSE.
#else 
      LAMBERTIAN_CLOUDS = .TRUE.
#endif

C-------------------------------------------------------------------------
C Keep TAU_IN constant
C-------------------------------------------------------------------------
       TAU_IN = GET_TAU0_6A( MONTH, DAYNO, YR )

C-------------------------------------------------------------------------
C Read the trace gas mixing ratio from the binary punch file
C-------------------------------------------------------------------------
      print*,"Read in trace gas mixing ratios"
      call flush(5)
      IF (TG_TYPE .EQ. 0) THEN
         TRACER_IN   = 20	! HCHO
         CATEGORY_IN = 'IJ-AVG-$'
      ELSE
#if defined ( GOME )
         TRACER_IN   = 23       ! NO2
#else
         TRACER_IN   = 64 !25  	! NO2	
#endif
         CATEGORY_IN = 'IJ-AVG-$' !'TIME-SER'
      ENDIF
      print*,"BPCH2 details:"
      print*, TG_FILE
      print*, CATEGORY_IN
      print*, TRACER_IN
      print*, TAU_IN
      print*, IIPAR
      print*, JJPAR
      call flush(5)
      

      CALL READ_BPCH2( TG_FILE, CATEGORY_IN, TRACER_IN, TAU_IN,  
     &                 IIPAR,    JJPAR,       LLPAR,     ARRAY )
      WRITE(6,*) '=================================
     &================================='
      TG(:,:,:) = ARRAY(:,:,:)

C-------------------------------------------
C Surface Pressure
C-------------------------------------------
      print*,"Read in surface pressure"
      call flush(5)
      TRACER_IN = 1
      CATEGORY_IN = 'PEDGE-$'

      CALL READ_BPCH2( TG_FILE, CATEGORY_IN, TRACER_IN, TAU_IN,  
     &                 IIPAR,    JJPAR,       LLPAR,     ARRAY )
      WRITE(6,*) '=================================
     &================================='
      PS(:,:) = ARRAY(:,:,1)

C-------------------------------------------
C Temperature
C-------------------------------------------
      print*,"Read in temperature"
      call flush(5)
      TRACER_IN = 3
      CATEGORY_IN = 'DAO-3D-$'

      CALL READ_BPCH2( TG_FILE, CATEGORY_IN, TRACER_IN, TAU_IN,  
     &                 IIPAR,    JJPAR,       LLPAR,     ARRAY )
      WRITE(6,*) '=================================
     &================================='
      TMP(:,:,:) = ARRAY(:,:,:)

C-------------------------------------------------------------------------
C Read the aerosol profiles from the binary punch file
C Skip this step if using Look Up Table
C-------------------------------------------------------------------------
      !Skip if using look up table
#if defined ( NOUSELUT ) 
#if defined ( GOME )
      TAU_IN = GET_TAU0_2A( MONTH, YR )
#endif
      CATEGORY_IN = 'OD-MAP-$'
      print*,"read in aerosols"
      call flush(5)
      print*,"I think there are this many aersols: ",NAER
      call flush(5)
      !===hard edits to force 
      DO N = 1, NAER
        TRACER_IN   = IND(N)

        CALL READ_BPCH2( AEROSOL_FILE, CATEGORY_IN, TRACER_IN, TAU_IN,  
     &                   IIPAR,    JJPAR,       LLPAR,     ARRAY )
        WRITE(6,*) '=================================
     &================================='
        ARSL(:,:,:,N) = ARRAY(:,:,:)
      ENDDO

C     Add coarse sea-salt OD to accum sea-salt OD
      print*,"read in sea salt"
      call flush(5)
      TRACER_IN = 18
      CALL READ_BPCH2( AEROSOL_FILE, CATEGORY_IN, TRACER_IN, TAU_IN,  
     &                 IIPAR,    JJPAR,       LLPAR,     ARRAY )
      ARSL(:,:,:,NAER) = ARSL(:,:,:,NAER) + ARRAY(:,:,:)

C To turn off aerosol correction, uncomment the following line
C     ARSL(:,:,:,:) = 0.0D0

C-------------------------------------------------------------------------
C Read the aerosol parameters
C-------------------------------------------------------------------------
      print*,"Now the aerosol parameters"
      call flush(5)
      OPEN ( FILE = 'jv_spec.dat', UNIT = LUN,
     &     FORM = 'FORMATTED', STATUS = 'OLD')
      READ(LUN,*)  !Header
      DO N = 1, NAER+1
        READ(LUN,*) !Aerosol label
        DO I = 1, NWL
          READ(LUN,*) WL(I),QEXT(I,N),RAA(I,N),ASSA(I,N),
     &                (PHFCN(J,I,N),J=0,7)
        ENDDO
      ENDDO
      CLOSE ( LUN )
#endif   ! End section to skip if using look up table

C-------------------------------------------------------------------------
C Read in UV surface albedos 
C USE TOMS surface albedos (~360 nm) for HCHO
C and Satellite surface albedos (~440 nm) for NO2
C-------------------------------------------------------------------------
      IF (TG_TYPE .EQ. 0) THEN
#if   defined( GRID4x5 )
        TG_FILE = 
     &    '/group_workspaces/jasmin/geoschem/ExtData/&
     &GEOS_4x5/uvalbedo_200111/uvalbedo.geos.4x5'
#elif defined( GRID2x25 ) 
        TG_FILE = 'uvalbedo.geos.2x25'
#elif defined ( GRIDNEST )
        TG_FILE = 'uvalbedo.geos.025x03125'
#endif
      ELSE

      !GOB implementation of 1x1 albedos
#if   defined( GOME1x1 )
        TG_FILE = 'ALBEDO/ler440.generic.1x1'
#elif defined( GOME2x25 )
        TG_FILE = 'ALBEDO/ler440.geos.pwd
2x25'
#elif defined( GOME4x5 )
        TG_FILE = 'ALBEDO/ler440.geos.4x5'
#elif defined( OMI1x1 )
        TG_FILE = 'ALBEDO/OMLER440.generic.1x1'
#elif defined( OMI2x25 )
        TG_FILE = 'ALBEDO/OMLER440.geos.2x25'
#elif defined( OMI4x5 )
        TG_FILE = 'ALBEDO/OMLER440.geos.4x5'
#else
        TG_FILE = 'ALBEDO/ler440.geos.' // GET_RES_EXT()
#endif
      ENDIF
      print*, TG_FILE
      call flush(5)
      CALL READ_BPCH2( TG_FILE, 'UVALBEDO', 1, TAU_VAL( MONTH ),
     &                 IIPAR_ALB,     JJPAR_ALB,     1, ARRAY_ALB )
      ALBD(:,:) = ARRAY_ALB(:,:)

C-------------------------------------------
C Tropopause level
C-------------------------------------------
      print*,"Read in tropopause levels"
      call flush(5)
      CATEGORY_IN = 'TR-PAUSE'
C      TAU_IN = GET_TAU0_6A( MONTH, 1, YR )
      TAU_IN = GET_TAU0_6A(MONTH , DAYNO, YR )
      print*,"TAU_IN is ",TAU_IN
      TG_FILE = TRIM(TROPO_FILE)
      ! changing tracer number to 2 here., LSURL 2015-11-03
      ! and back to 1... LSURL 2016-02-02
      
      !Made a 3D read LSURL 2016-08-24
      CALL READ_BPCH2( TG_FILE, CATEGORY_IN, 15, TAU_IN,  
     &                 IIPAR,    JJPAR,   LLPAR,   ARRAY )
      WRITE(6,*) '============================
     &======================================'
      TRPAUSE(:,:) = ARRAY(:,:,1)

C-------------------------------------------------------------------------
C Read in slant columns
C-------------------------------------------------------------------------
        print*,"Read in slant columns"
        call flush(5)

C  Read LIDORT input file
C  ======================
        CALL LIDORT_V23E_INPUT
     &     ( 'amfgas_lidortinput.vza',
     &        'ACTIVEGAS_ERRORFILE', STATUS_INPUTREAD )
        IF ( STATUS_INPUTREAD .NE. LIDORT_SUCCESS ) THEN
          STOP' Problems encountered with input read'
        ENDIF

C  read USA ML atmosphere
        OPEN(1,FILE='usaml_ztp.dat',STATUS='OLD')
        DO J = 1, NUSAML
          J00 = NUSAML - J + 1
          READ(1,*)USAML_ZZZ(J00),USAML_TTT(J00),USAML_PPP(J00)
          USAML_LNP(J00) = DLOG(USAML_PPP(J00))
        ENDDO
        CLOSE(1)

C-------------------------------
C Read in Satellite Information
C-------------------------------
       print*,"Read satellite information"
       call flush(5)

#if defined ( GOME )
       CALL input_GOME(SATFILE,TG_TYPE,ULN,OLATS,OLONS,SZA,
     &                   LOSZA,RELAZM,PCLD,FCLD,TCOT,SC,
     &                   MAXLINE,outPrefix,outSuffix,flag)
#endif

#if defined ( NEW_GOME )
       CALL input_newGOME(SATFILE,TG_TYPE,OLATS,OLONS,SZA,
     &                   LOSZA,RELAZM,PCLD,FCLD,TCOT,SC,
     &                   MAXLINE,outPrefix,outSuffix,flag)
#endif

#if defined ( SCIA )
       print*,"It's time for input SCIA"
       call flush (5)
       CALL input_SCIA(SATFILE,TG_TYPE,FRESCOv5,OLATS,OLONS,SZA,
     &                   LOSZA,RELAZM,PCLD,FCLD,FPS,SC,
     &                   MAXLINE,outPrefix,outSuffix,flag)
#endif

#if defined ( OMI )
C      LSURL 2015-11-24 - change this to the NEW GOME reader (ASCII).
C      CALL input_OMI(SATFILE,TG_TYPE,CLDFILE,OLATS,OLONS,SZA,
C     &                   LOSZA,RELAZM,PCLD,FCLD,
C     &                   MAXLINE,outPrefix,outSuffix,flag)
       CALL input_newGOME(SATFILE,TG_TYPE,ULN,OLATS,OLONS,SZA,
     &                   LOSZA,PCLD,FCLD,
     &                   MAXLINE,outPrefix,outSuffix,flag)
     
#endif


      SLEN = INDEX ( SATFILE, ' ' ) - 1
      WRITE(6,*) 'Read in', MAXLINE, 'observations from ', 
     &     SATFILE(1:SLEN)


C==================================================
C
C  Begin AMF Calculation
C
C==================================================
      print*,"Begin AMF calculation"
      call flush(5)

!$omp parallel do
!$omp+ default(shared)
!$omp+ private(LINE)
!$omp+ firstprivate(TG_TYPE)
!$omp+ private(IUSE,JUSE,IUSE_ALB,JUSE_ALB)
!$omp+ private(AERPROF,MSIGMAE,MSIGMA, LEV)
!$omp+ private(AIRMIX,TGMIX,RSZA,RVZA)
!$omp+ private(RCOT,RPCLD,RSPIN,RSAIN)
!$omp+ private(htopatm,eps,xh)
!$omp+ private(KHATCLD,ACLR,ACLD)
!$omp+ private(DSIGMA,AIRSIGMA,TGSIGMA)
!$omp+ private(TMPSIGMA,WEIGHT)
!$omp+ private(SCATCLR,SCATCLD,RFCLD,TRAD)


      DO LINE = 1,MAXLINE
         if (MOD(LINE,1000) .eq. 0) then
            print*,LINE ," of ", MAXLINE
            call flush(5)
            !print*, "TG_TYPE is ", TG_TYPE
            !call flush(5) 
         end if
        !if (mod(line,1000)==1) then
        !    write(*,*) line
	!end if

#if defined ( USELUT )
         ! Skip AMF computation for data incompatible with LUT code

          if ( (SZA(line) .gt. 84)        !GOB Added so LUT code won't crash
     &      .or. (SZA(LINE) .lt. 1) ) then
            AMF   (LINE) = -1.0D30
            AMFCLR(LINE) = -1.0D30
            AMFCLD(LINE) = -1.0D30
            cycle
          end if
#endif
         !print*,OLATS(LINE)
         !print*,OLONS(LINE)
         !print*,JMIN
         !print*,JMAX
         !print*,IMIN
         !print*,IMAX
         !call flush(5)

         !Skip calculation if outside the nest
        if ((OLATS(LINE) .lt. JMIN).or.
     &      (OLATS(LINE) .gt. JMAX).or.
     &      (OLONS(LINE) .lt. IMIN).or.
     &      (OLONS(LINE) .gt. IMAX)) then
            AMF   (LINE) = -1.0D30
            AMFCLR(LINE) = -1.0D30
            AMFCLD(LINE) = -1.0D30
            !print*,"outside nest"
            !call flush(5)
            cycle
        end if
        !print*,"inside nest and SZA OK!"
        call flush(5) 

#if     defined( OMINEST )
        IUSE = 1 + (OLONS(LINE)-IMIN)/DISIZE
        JUSE = 1 + (OLATS(LINE)-JMIN)/DJSIZE
#else 


C     Determine the model grid cell the observation falls in
         IF (OLONS(LINE) .LT. IMAX) THEN
           IUSE = 1+(OLONS(LINE)-IMIN)/DISIZE
         ELSE
           IUSE = 1+(OLONS(LINE)-IMAX)/DISIZE
         ENDIF

         IF (OLATS(LINE) .LT. JMIN) THEN
           JUSE = 1
         ELSE 
           JUSE = 2+(OLATS(LINE)+JMAX)/DJSIZE
         ENDIF
#endif

        !GOB Implementation of albedos maps
#if     defined( OMI4x5 ) || defined( GOME4x5 ) || defined( OMI2x25 ) || defined( GOME2x25 )
         IF (OLONS(LINE) .LT. IMAX_ALB) THEN
           IUSE_ALB = 1+(OLONS(LINE)-IMIN_ALB)/DISIZE_ALB
         ELSE
           IUSE_ALB = 1+(OLONS(LINE)-IMAX_ALB)/DISIZE_ALB
         ENDIF

         IF (OLATS(LINE) .LT. JMIN_ALB) THEN
           JUSE_ALB = 1
         ELSE
           JUSE_ALB = 2+(OLATS(LINE)+JMAX_ALB)/DJSIZE_ALB
         ENDIF

#elif   defined( OMI1x1 ) || defined( GOME1x1 )
        !GOB Implementation of albedos maps
         IF (OLONS(LINE) .LT. IMAX_ALB) THEN
           IUSE_ALB = 1+(OLONS(LINE)-IMIN_ALB)/DISIZE_ALB
         ELSE
           IUSE_ALB = 1+(OLONS(LINE)-IMAX_ALB)/DISIZE_ALB
         ENDIF

         IF (OLATS(LINE) .EQ. JMAX_ALB) THEN
           JUSE_ALB = (OLATS(LINE)+JMAX_ALB)/DJSIZE_ALB
         ELSE
	   JUSE_ALB = 1+(OLATS(LINE)+JMAX_ALB)/DJSIZE_ALB
         ENDIF
#else
         !NESTED
         IUSE_ALB=IUSE;
         JUSE_ALB=JUSE;
#endif
C         print*,"AMF dugub point BB"
C         print*,JUSE
C         print*,JMIN
C         print*,JMAX
C         print*,DJSIZE
C         print*,OLATS(LINE)
         call flush(5)

            IF ( SUM(TG(IUSE,JUSE,:)) .GT. 0.0D0 ) THEN
             IF ( PCLD(LINE) .GT. 0.0 ) THEN

               IF (FRESCOv5 .EQ. 1) THEN
                  SPIN(LINE) = FPS(LINE)
               ELSE
                  SPIN(LINE)  = PS(IUSE,JUSE)
               ENDIF


C Albedos range from 0-1
               SAIN(LINE)   = ALBD(IUSE_ALB,JUSE_ALB)

C Aerosol profile
               AERPROF(:,:) = ARSL(IUSE,JUSE,:,:)

C               print*,"AMF degub point CC"
               call flush(5)
C------------------------------------------------------------------
C Calculate Eta Levels for GEOS_4 or GEOS_5
C------------------------------------------------------------------
#if   defined( GEOS_4 ) || defined( GEOS_5 )
               MSIGMAE(1) = ETA_A(1)/SPIN(LINE) + ETA_B(1)
               DO LEV = 1, LLPAR-1
                 MSIGMAE(LEV+1) = ETA_A(LEV+1)/SPIN(LINE) + ETA_B(LEV+1)
                 MSIGMA(LEV) = (MSIGMAE(LEV) + MSIGMAE(LEV+1) ) / 2.0d0
               ENDDO
               MSIGMAE(LLPAR+1) = 0.0d0
               MSIGMA(LLPAR)=(MSIGMAE(LLPAR)+MSIGMAE(LLPAR+1) )/2.0d0
#endif
C               print*,"AMF degub point DD"
               call flush(5)
       
C------------------------------------------------------------------
C Mixing ratio for air and Trace Gas
C------------------------------------------------------------------
               DO LEV = 1, LLPAR
                  AIRMIX(LEV) = 1.0D0
                  TGMIX(LEV) = TG(IUSE,JUSE,LEV)
               ENDDO

C               print*,"AMF degub point EE"
               call flush(5)
C Change from REAL*4 to REAL*8
               RSZA = SZA(LINE)
               RVZA = ABS(LOSZA(LINE))
C               RCOT = TCOT(LINE)  NOT using TCOT,
               RPCLD = PCLD(LINE)
               RSPIN = SPIN(LINE)
               RSAIN = SAIN(LINE)

               if (RPCLD .gt. RSPIN) then   !(GOB) move clouds up to surface
		   RPCLD = RSPIN
	       end if

C               print*,"AMF degub point FF"
               call flush(5)


!------------------------------------------------------------------
! Geometric AMF calculation
!
! From Carsten Leue, thesis:
!     AMFgeo = (1.0/cos(pi*vza/180.0)) + (1.0/cos(pi*sza/180.0))
!
! NOTE: the AMF lookup table contains AMF/AMFgeo, where
!         AMFgeo has exactly the same form as here!
!------------------------------------------------------------------

               htopatm = 60.0   ! km

               eps     = htopatm * 1000.0 / EarthRadius

               xh      = sqrt( cos(pi * SZA(LINE) / 180.0)**2 +
     &              eps**2 + 2.*eps )

               AMFgeo(LINE)  = 1.0 / cos( pi * LOSZA(LINE) / 180.0 ) +
     &              ( xh - cos( pi * SZA(LINE) / 180.0 ) ) / eps


!$OMP CRITICAL

#if defined ( USELUT )
C------------------------------------------------------------------
C Call to read weighting kernels from Look Up Table
C------------------------------------------------------------------
               CALL LUT_USE( RSPIN, RSAIN,
     &                       SZA(LINE), LOSZA(LINE), RELAZM(LINE),
     &                       AMFgeo(LINE), TRPAUSE(IUSE,JUSE),
     &                       PCLD(LINE),KHAT(LINE,:),KHATCLD,ACLR, ACLD)

#else

C------------------------------------------------------------------
C If not using look up table call lidort code
C Get weighting kernel for clear sky case
C------------------------------------------------------------------
C               print*,"AMF degub point GG"
C               call flush(5)               
               CALL RLIDORT(SPIN(LINE),SAIN(LINE),
     &              RSZA,TG_TYPE,0.0D0,900.0D0,
     &              TGMIX,KHAT(LINE,:),ACLR,RVZA,
     &              AERPROF,.FALSE.)

C------------------------------------------------------------------
C Get weighting kernel for cloudy case 
C------------------------------------------------------------------
C               print*,"AMF degub point HH"
C               call flush(5)
               CALL RLIDORT(SPIN(LINE),SAIN(LINE),
     &              RSZA,TG_TYPE,RCOT,RPCLD,
     &              TGMIX,KHATCLD,ACLD,RVZA,
     &              AERPROF,LAMBERTIAN_CLOUDS)

C              write(*,'(5f8.3)') ACLR, ACLD, RPCLD
#endif

!$OMP END CRITICAL

C-------------------------------------------------------------------
C Compute total column of air, Trace Gas, TG Weighted Mean Temperature
C------------------------------------------------------------------
               AIRSIGMA = 0.0D0
               DO LEV = 1, LLPAR
                  DSIGMA(LEV) = MSIGMAE(LEV) - 
     &                 MSIGMAE(LEV+1)
                  AIRSIGMA = DSIGMA(LEV)*AIRMIX(LEV) + AIRSIGMA
               ENDDO

                     
               TGSIGMA = 0.0D0
               TMPSIGMA = 0.0D0
               DO LEV = 1, TRPAUSE(IUSE,JUSE)-1
                  WEIGHT = DSIGMA(LEV)*TGMIX(LEV)
                  TGSIGMA = WEIGHT + TGSIGMA
                  TMPSIGMA = WEIGHT*TMP(IUSE,JUSE,LEV) + TMPSIGMA
               ENDDO


C------------------------------------------------------------------
C Compute S function
C------------------------------------------------------------------
               DO LEV = 1, TRPAUSE(IUSE,JUSE)-1
                  IF ( TGMIX(LEV) .GT. 0.0D0 ) THEN
                     S(LINE,LEV) = TGMIX(LEV)
                     S(LINE,LEV) = S(LINE,LEV)/
     &                    (TGSIGMA/AIRSIGMA)
                  ELSE
                     S(LINE,LEV) = 0.0D0
                  ENDIF
               ENDDO

C Check that S integrates to unity - use AMF variable to check this
               AMF(LINE) = 0.0D0
               DO LEV = 1, TRPAUSE(IUSE,JUSE)-1
                  AMF(LINE) = DSIGMA(LEV)*S(LINE,LEV) 
     &                 + AMF(LINE)
               ENDDO

C Allow for numerical error (especially when using REAL precision)
               IF ( AMF(LINE) .GT. 0.99 .AND. 
     &              AMF(LINE) .LT. 1.01) THEN
               ELSE
                  WRITE(6,*) 'Problem with model shape
     &                         integrating to unity'
                  WRITE(6,*) AMF(LINE)
                  !STOP
                  !let's not abort if this happens, just not do this AMF
                  AMF   (LINE) = -1.0D30
                  AMFCLR(LINE) = -1.0D30
                  AMFCLD(LINE) = -1.0D30
                  !print*,"outside nest"
                  !call flush(5)
                  cycle


               ENDIF

C------------------------------------------------------------------
C Compute air mass factor
C Satellite is more sensitive to NO2 at colder temps so scale KHAT
C   to account for this
C------------------------------------------------------------------
               AMFCLR(LINE) = 0.0D0
               AMFCLD(LINE) = 0.0D0
               
               DO LEV = 1, TRPAUSE(IUSE,JUSE)-1
                  IF (TG_TYPE .EQ. 1) THEN
#if defined ( GOME )
C Larger cross sections at 243K but spectral fit at 293K 
C Slope of 0.2 determined by comparing results for fits at each temp
                    WEIGHT = 1.0D0 + 0.2D0*(293.0D0-TMP(IUSE,JUSE,LEV))
     &                       /(293.0D0-243.0D0)
#else
C Smaller cross sections at 293K but spectral fit at 243K 
                    WEIGHT = 1.0D0 + 0.2D0*(243.0D0-TMP(IUSE,JUSE,LEV))
     &                       /(293.0D0-243.0D0)
#endif
                  ELSE
                    WEIGHT = 1.0D0 
                  ENDIF

                  AMFCLR(LINE) = (DSIGMA(LEV)*KHAT(LINE,LEV)*
     &                 S(LINE,LEV)*WEIGHT) + AMFCLR(LINE)
                   
                   
                  AMFCLD(LINE) = (DSIGMA(LEV)*KHATCLD(LEV)*
     &                 S(LINE,LEV)*WEIGHT) + AMFCLD(LINE)
                
                  !Scattering Weights for Clear and Cloudy Skies calculated
                  SCATCLR(LEV) = (DSIGMA(LEV)*KHAT(LINE,LEV)*WEIGHT)
  
                  SCATCLD(LEV) = (DSIGMA(LEV)*KHATCLD(LEV)*WEIGHT)

                  ENDDO

C Combine clear and cloudy Scattering weights
               RFCLD = FCLD(LINE)

#if defined ( GOME ) || defined ( SCIA ) || defined ( NEW_GOME )
               RFCLD = 0.01D0*RFCLD ! OMI CLDRR cloud fraction 0-1
#endif

               ! estimated observed radiance (normalization factor)
               TRAD = ACLR*(1.0D0-RFCLD) + ACLD*RFCLD

               ! Radiance fraction from clear part of pixel
               WCLR(LINE) = ACLR*(1.0D0-RFCLD)/TRAD

               ! Radiance fraction from cloudy part of pixel
               WCLD(LINE) = ACLD*RFCLD/TRAD

               ! Total AMF weighted by cloud fraction and
               ! radiation intensity
                AMF(LINE) = AMFCLR(LINE)*WCLR(LINE) 
     &              + AMFCLD(LINE)*WCLD(LINE)
               
              
               AMF(LINE)= 0.0D0 
               !Total Scattering weights calculated by weighting
               !cloud fraction and radiation intensity
               !AMF is calculated by integrating
               DO LEV = 1, TRPAUSE(IUSE,JUSE)-1
                  
                            
                  WSCAT(LINE,LEV)=(SCATCLD(LEV)*WCLD(LINE)
     &                 +SCATCLR(LEV)*WCLR(LINE))/AMFgeo(LINE)
               
                  
                  AMF(LINE)= (WSCAT(LINE,LEV)*S(LINE,LEV)*AMFgeo(LINE))+
     &                 AMF(LINE)
               ENDDO
  

               
        
C       write(*,'(5f8.3)') amfcld(line), amfclr(line), amf(line),
C     &                     wcld(line), wclr(line), AMFgeo(line)

C Compute vertical column
               GEOS_VC(LINE) = ((AVGRAD*(SPIN(LINE)-PTOP))/
     &              (G*MA) * TGSIGMA)/100.0D0


            ENDIF               ! PCLD is > 0.0
         ENDIF                  ! model TG is > 0.0
      ENDDO
!$omp end parallel do     

C-------------------------------------------------------------------------
C Open main output file
C-------------------------------------------------------------------------
      OPEN ( FILE = OUTFILE,
     &     UNIT = OLUN,
     &     STATUS = 'UNKNOWN',
     &     FORM = 'FORMATTED' )

      IF ( MAXLINE==0 ) THEN
         WRITE( 6, * ) '==No suitable observations found=='
      ELSE
         DO LINE = 1, MAXLINE
C     LUKE 2016-08-24 Write every line, we'll filter later.         
C           IF ((AMFCLD(LINE).LT.9.9D0) .AND. (AMFCLR(LINE).LT.9.9D0)
C     &     .AND. (AMFCLD(LINE).GE.0.0D0) .AND. (AMFCLR(LINE).GT.0.0D0)
C     &     .AND. (AMF(LINE).GT.0.0D0) .AND. (AMF(LINE).LT.9.9D0)) THEN
              
               WRITE(OLUN,'(a12,a1,I9,a1,E15.4,a1,E15.4,a1,a12)')
     &              
     &                 outPrefix(LINE),",",
     &                 ULN(LINE),",",
     &                 AMF(LINE),",",
     &                 GEOS_VC(LINE),",",
     &                 outSuffix(LINE)
C           ENDIF
        ENDDO
      ENDIF

      CLOSE ( OLUN ) ! Close output file


C   Don't ouput shape factor and scatering weights for OMI
C LUKE 2016-08-24 Don't write these useless files
C
C#if defined ( GOME ) || defined ( SCIA ) || defined ( NEW_GOME )
C
CC-------------------------------------------------------------------------
CC Open output file - for WSCAT
CC-------------------------------------------------------------------------
C      OPEN ( FILE = OUTFILE_WSCAT,
C     &     UNIT = OLUN,
C     &     STATUS = 'UNKNOWN',
C     &     FORM = 'FORMATTED' )
C
C      IF ( MAXLINE .ne. 0 ) THEN
C         DO LINE = 1, MAXLINE
C           IF ((AMFCLD(LINE).LT.9.9D0) .AND. (AMFCLR(LINE).LT.9.9D0)
C     &     .AND. (AMFCLD(LINE).GE.0.0D0) .AND. (AMFCLR(LINE).GT.0.0D0)
C     &     .AND. (AMF(LINE).GT.0.0D0) .AND. (AMF(LINE).LT.9.9D0)) THEN
C
C
C             WRITE(OLUN,'(a12, 47f9.5)')
C     &          outPrefix(LINE), KHAT(LINE,:)
C
C           ENDIF
C         ENDDO
C      ENDIF
C
C      CLOSE ( OLUN ) ! Close output file
C
CC-------------------------------------------------------------------------
CC Open output file - for shape factors
CC-------------------------------------------------------------------------
C      OPEN ( FILE = OUTFILE_Shape,
C     &     UNIT = OLUN,
C     &     STATUS = 'UNKNOWN',
C     &     FORM = 'FORMATTED' )
C
C      IF ( MAXLINE .ne. 0 ) THEN
C         DO LINE = 1, MAXLINE
C           IF ((AMFCLD(LINE).LT.9.9D0) .AND. (AMFCLR(LINE).LT.9.9D0)
C     &     .AND. (AMFCLD(LINE).GE.0.0D0) .AND. (AMFCLR(LINE).GT.0.0D0)
C     &     .AND. (AMF(LINE).GT.0.0D0) .AND. (AMF(LINE).LT.9.9D0)) THEN
C
C             WRITE(OLUN,'(a12, 30e15.8)')
C     &          outPrefix(LINE), S(LINE,:)
C
C           ENDIF
C         ENDDO
C      ENDIF
C
C      CLOSE ( OLUN ) ! Close output fil
C#endif

      STOP
      END
