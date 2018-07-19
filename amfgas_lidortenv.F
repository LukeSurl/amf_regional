      SUBROUTINE RLIDORT(SPIN,SAIN,SZA,TG_TYPE,CLDODIN,PCLD,
     &     ACTIVEGAS_VMR,MWGHTS,TOAALB,VZA,AERPROF,
     &     LAMBERTIAN_CLOUDS)

C
C Call LIDORT 
C Code originally written by Rob Spurr for HCHO, Spring 2000.
C Modified by Randall Martin 
C

        USE CMN_SIZE_MOD      ! Size parameters

        IMPLICIT NONE

C  include fileS of dimensions and numbers

#	include "include_s/LIDORT.PARS"
#	include "include_e/LIDORT_L.PARS"

C  input variables

#	include "include_s/LIDORT_CONTROL.VARS"
#	include "include_s/LIDORT_MODEL.VARS"
#	include "include_s/LIDORT_GEOPHYS.VARS"
#	include "include_e/LIDORT_L_CONTROL.VARS"

C  LIDORT result variables

#	include "include_s/LIDORT_RESULTS.VARS"
#	include "include_e/LIDORT_L_RESULTS.VARS"

C  input values

      INTEGER TG_TYPE           ! 1=NO2, 0=HCHO

      REAL*8 SPIN                 ! Input surface pressure
      REAL*8 SAIN                 ! Input surface albedo
      REAL*8 SZA                  ! Input solar zenith angle
      REAL*8 VZA                  ! Input viewing zenith angle
      REAL*8 CLDODIN              ! Input cloud optical depth
      REAL*8 PCLD                 ! Input cloud top pressure

      REAL*8 TGMIX(LLPAR)        ! Model profile at sigma coordinates 
      REAL*8 AERPROF(LLPAR,NAER) ! Aerosol profile at sigma coordinates 
      LOGICAL LAMBERTIAN_CLOUDS

C  output values
      REAL*8 MWGHTS(LLPAR)       ! Weights on model grid
      REAL*8 TOAALB              ! TOA albedo

C  dimensioning
      REAL*8 Z(MAXLAYER) ! Altitude
      REAL*8 DZ

C  toggle between HCHO and NO2
	CHARACTER*4	 TG

C  Maximum number of other gases

	INTEGER		MAX_OTHERGAS
	PARAMETER	( MAX_OTHERGAS = 2 )

C  help variables

	INTEGER		 I,UI,n,K,m,G,REV, L
	REAL*8           LAMDA, ACTIVEGAS_OPDEP,INTENS,SSA,
     &                   ACTIVEGAS_COLUMN,ACTIVEGAS_EFFXSEC
	INTEGER		 STATUS_INPUTCHECK
	INTEGER		 STATUS_CALCULATION
 	REAL*8		 e1, e2, AMF, SMKM
      PARAMETER       ( SMKM=1.D+05 )

C  local storage arrays for results

	REAL*8 VMRWF(MAXLAYER)

C  Local data 
C  -----------------------------------
        REAL*8 KHAT(MAXLAYER)        ! Weights on LIDORT grid
        REAL*8 PSURF                 ! Pressure of model atmosphere

C  Local data from prepared atmosphere
C  -----------------------------------

	INTEGER		 NGRID,N_OTHERGAS
	CHARACTER*4	 ACTIVEGAS_NAME, OTHERGAS_NAMES(MAX_OTHERGAS)
        REAL*8 FP
        PARAMETER       ( FP=7.24325825D+18 )

	REAL*8 ZGRID(MAXLAYER)
	REAL*8 DELTA(MAXLAYER)
	REAL*8 AERSCAT(MAXLAYER)
	REAL*8 AEREXTC(MAXLAYER)
	REAL*8 PHASEFCN(MAXLAYER,0:7)
	REAL*8 OTHERGAS_VMRS(MAXLAYER,MAX_OTHERGAS)
	REAL*8 ACTIVEGAS_CONC(MAXLAYER)
	REAL*8 CONC(MAXLAYER)
	REAL*8 ACTIVEGAS_VMR(MAXLAYER)
	REAL*8 PRESSURE(MAXLAYER)
	REAL*8 SIGMA(MAXLAYER)
	REAL*8 VMR(MAXLAYER)
	REAL*8 TEMPERATURE(MAXLAYER)
	COMMON / READATMOS /
     &      ZGRID, DELTA, PRESSURE, TEMPERATURE,
     &      ACTIVEGAS_CONC,
     &      OTHERGAS_VMRS, AERSCAT, AEREXTC,
     &      PHASEFCN, NGRID, N_OTHERGAS,
     &      ACTIVEGAS_NAME, OTHERGAS_NAMES
	SAVE / READATMOS /

#if   defined( GEOS_4 ) || defined (GEOS_5)
        MSIGMAE(1) = ETA_A(1)/SPIN + ETA_B(1)
        DO L = 1, LLPAR-1
          MSIGMAE(L+1) = ETA_A(L+1)/SPIN + ETA_B(L+1)
          MSIGMA(L) = (MSIGMAE(L) + MSIGMAE(L+1) ) / 2.0d0
        ENDDO
        MSIGMAE(LLPAR+1) = 0.0d0
        MSIGMA(LLPAR)=(MSIGMAE(LLPAR)+MSIGMAE(LLPAR+1) )/2.0d0
#endif
C        print*,"LIDORT point 001"
C        call flush(5)
C  ===============================================
      DO_LAMBERTIAN_ALBEDO = .TRUE.
      SUN0 = SZA         ! solar zenith angle 
      UI = 1     
      USER_ANGLES_INPUT(UI) = VZA !viewing angle

      IF ( LAMBERTIAN_CLOUDS ) THEN
        ALBEDO = 0.8D0
        CLDODIN = ZERO
        PSURF = PCLD
      ELSE
        ALBEDO = SAIN      ! surface albedo
        PSURF = SPIN
      ENDIF
      
C      print*,"LIDORT point 002"
C      call flush(5)
C  control input for the data base
C  ===============================

      IF (TG_TYPE .EQ. 0) THEN
         ACTIVEGAS_NAME   = 'HCHO'       
         LAMDA = 346.03766
         N_OTHERGAS = 2
         OTHERGAS_NAMES(1) = 'OZON'
         OTHERGAS_NAMES(2) = 'NO2 '
      ELSE
         ACTIVEGAS_NAME   = 'NO2 '          
         LAMDA = 438.19500D0
         N_OTHERGAS = 1
         OTHERGAS_NAMES(1) = 'OZON'
         OTHERGAS_NAMES(2) = '   '
      ENDIF

C      print*,"LIDORT point 003"
C      call flush(5)
C  get trace gases and ZPT from made-up atmosphere
C  ===============================================

       CALL MKATM
     I           (PSURF,CLDODIN,PCLD,ACTIVEGAS_VMR,LAMDA,
     I            AERPROF,
     O            ZGRID,NGRID,DELTA,
     O            PRESSURE,TEMPERATURE,OTHERGAS_VMRS,
     O            ACTIVEGAS_CONC,AERSCAT,AEREXTC,PHASEFCN)
      NLAYERS = NGRID

C Test output
c     OPEN ( FILE = 'test.out',
c    &     UNIT = 83, STATUS = 'UNKNOWN', FORM = 'FORMATTED' )
c     do n = 1,ngrid
c       write(83,*) zgrid(n),zvmr(n),pressure(n),temperature(n)
c     enddo
c     CLOSE ( 83 ) ! Close output file
C      print*,"LIDORT point 004"
C      call flush(5)
C  prepare LIDORT atmosphere
C  -------------------------

	CALL ULIDORT_PREP
     I      ( LAMDA, 
     O        ACTIVEGAS_OPDEP, ACTIVEGAS_COLUMN, ACTIVEGAS_EFFXSEC  )

C        print*,"LIDORT point 005"
C        call flush(5)

C  call to LIDORT
C  -------------------------
	CALL LIDORT_V23E_MASTER
     &                ( STATUS_INPUTCHECK, STATUS_CALCULATION )

C        print*,"LIDORT point 006"
C        call flush(5)

C  Check LIDORT output status

	IF ( STATUS_INPUTCHECK. NE. LIDORT_SUCCESS ) THEN
          write(*,*)'LIDORT error'
          DO N = 1,NLAYERS
            SSA = AERSCAT(N)/AEREXTC(N)
            write(*,*) zgrid(n), ssa, aerscat(n), aerextc(n)
          ENDDO
	  STOP' Problems encountered LIDORT input check'
	ENDIF
c       checked below ok
C       write(*,*) ACTIVEGAS_OPDEP, ACTIVEGAS_COLUMN, ACTIVEGAS_EFFXSEC 

C  save output
C  -----------
	 INTENS = INTENSITY(1,UI,UPIDX,1)
         TOAALB = FLUX(1,UPIDX)

C        numerical issues
         IF ( INTENS. LT. ZERO ) THEN
          INTENS = ZERO
         ENDIF

C-------------------------------------------------------
C Reverse order: surface is the 1st element
C-------------------------------------------------------

	 DO N = 1, NLAYERS
           IF ( LAMBERTIAN_CLOUDS ) THEN
             REV = (NLAYERS-N+5)
           ELSE
             REV = (NLAYERS-N+1)
           ENDIF
           Z(REV) = ZGRID(N)
	   VMRWF(REV) = ATMOSWF(1,N,1,UI,UPIDX,1)
           CONC(REV) = ACTIVEGAS_CONC(N)*SMKM
           SIGMA(REV) = (PRESSURE(N) - PTOP)/(SPIN - PTOP)
	 ENDDO

C        Set KHAT to zero below cloud
         IF ( LAMBERTIAN_CLOUDS ) THEN
           CONC(1:4) = ZERO
           SIGMA(4) = (PCLD - PTOP)/(SPIN - PTOP) - 1.0D-3
           SIGMA(3) = SIGMA(4) + 1.0D-3
           SIGMA(2) = SIGMA(3) + 1.0D-3
           SIGMA(1) = ONE
         ENDIF

C----------------------------------------------------------------------
C Compute $\frac{{\hat{w}}(z)}{I\overline{\sigma}}$, which is 
C the quantity to be interpolated. 
C----------------------------------------------------------------------
	 DO N = 1, NLAYERS-1
           DZ = Z(N) - Z(N+1)
           IF (INTENS*DZ*ACTIVEGAS_EFFXSEC*(CONC(N)+CONC(N+1)) 
     &         .EQ. ZERO) THEN
             KHAT(N) = ZERO
           ELSE 
             KHAT(N) = (VMRWF(N)+VMRWF(N+1))/
     &       (INTENS*DZ*ACTIVEGAS_EFFXSEC*(CONC(N)+CONC(N+1)))
           ENDIF
C          numerical issues
    	   IF ( KHAT(N). LT. ZERO ) THEN
            KHAT(N) = ZERO
  	   ENDIF
	 ENDDO

C----------------------------------------------------------------------
C Interpolate weights onto model sigma coordinates
C----------------------------------------------------------------------
       CALL LAGRANGE_INTERP_MASTER
     I  ( .TRUE., .TRUE.,
     I      SIGMA, KHAT, MAXLAYER, NGRID, 
     I      1, 1, 1.0D-12, LLPAR, LLPAR, MSIGMA, MWGHTS )

C Test interpolation
c     OPEN ( FILE = 'test1.out',
c    &     UNIT = 83, STATUS = 'UNKNOWN', FORM = 'FORMATTED' )
c     OPEN ( FILE = 'test2.out',
c    &     UNIT = 84, STATUS = 'UNKNOWN', FORM = 'FORMATTED' )
c     do i = 1,llpar
c       write(83,*) mwghts(i),msigma(i)
c     enddo
c     do i = 2,ngrid
c       write(84,*) khat(i),sigma(i)
c     enddo
c     CLOSE ( 83 ) ! Close output file
c     CLOSE ( 84 ) ! Close output file

	RETURN
	END

