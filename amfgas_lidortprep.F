
	SUBROUTINE ULIDORT_PREP
     &    ( LAMDA, 
     O      ACTIVEGAS_OPDEP, ACTIVEGAS_COLUMN, ACTIVEGAS_EFFXSEC )

        IMPLICIT NONE

C  include files of dimensions and numbers

#	include "include_s/LIDORT.PARS"
#	include "include_e/LIDORT_L.PARS"

C  input variables

#	include "include_s/LIDORT_CONTROL.VARS"
#	include "include_s/LIDORT_MODEL.VARS"
#	include "include_s/LIDORT_GEOPHYS.VARS"
#	include "include_e/LIDORT_L_CONTROL.VARS"
#	include "include_e/LIDORT_L_GEOPHYS.VARS"

C  module arguments

	REAL*8 LAMDA, ACTIVEGAS_OPDEP
	REAL*8 ACTIVEGAS_COLUMN, ACTIVEGAS_EFFXSEC

C  Maximum number of other gases

	INTEGER		MAX_OTHERGAS
	PARAMETER	( MAX_OTHERGAS = 2 )

C  Local data from prepared atmosphere
C  -----------------------------------

C  Atmospheric read variables

	INTEGER		 NGRID, N_OTHERGAS
	CHARACTER*4	 ACTIVEGAS_NAME, OTHERGAS_NAMES(MAX_OTHERGAS)
	REAL*8 ZGRID(MAXLAYER)
	REAL*8 DELTA(MAXLAYER)
	REAL*8 AERSCAT(MAXLAYER)
	REAL*8 AEREXTC(MAXLAYER)
	REAL*8 PHASEFCN(MAXLAYER,0:7)
	REAL*8 OTHERGAS_VMRS(MAXLAYER,MAX_OTHERGAS)
	REAL*8 ACTIVEGAS_CONC(MAXLAYER)
	REAL*8 PRESSURE(MAXLAYER)
	REAL*8 TEMPERATURE(MAXLAYER)
	COMMON / READATMOS /
     &      ZGRID, DELTA, PRESSURE, TEMPERATURE,
     &      ACTIVEGAS_CONC, 
     &      OTHERGAS_VMRS, AERSCAT, AEREXTC,
     &      PHASEFCN, NGRID, N_OTHERGAS,
     &      ACTIVEGAS_NAME, OTHERGAS_NAMES
	SAVE / READATMOS /

C  aerosol

	REAL*8 EXTCO_A(MAXLAYER)
	REAL*8 SCACO_A(MAXLAYER)
	REAL*8 PHASMOMS_INPUT_R(0:7,MAXLAYER)

C  preparation for a spherical atmosphere

	REAL*8 ZS_GRID(0:MAXLAYER)
	REAL*8 PS_GRID(0:MAXLAYER)
	REAL*8 TS_GRID(0:MAXLAYER)
	REAL*8 REARTH, NREFRAC

C  Spherical prepare output
C  fine gridding and slant path optical depths; solar zenith cosines

	INTEGER		 MAXFINELAYER
	PARAMETER	 ( MAXFINELAYER = 64 )
	INTEGER		 FINEGRID(MAXLAYER)
	REAL*8
     &      TAUTHICK_FINE_INPUT(MAXLAYER,MAXLAYER,0:MAXFINELAYER)
	REAL*8 TAUPATH_FINE_INPUT(MAXLAYER,0:MAXFINELAYER)
	REAL*8 MUPATH_FINE(MAXLAYER,0:MAXFINELAYER)

C  local parameters

	REAL*8 FP, DUKM_TO_CONC2,SMKM
	PARAMETER	( DUKM_TO_CONC2 = 2.68675D+16 )
	PARAMETER	( FP=7.24325825D+18 )
	PARAMETER	( SMKM=1.D+05 )

C  Constants for the Bates/Chance formula

        REAL*8  ABATES, BBATES, CBATES
        REAL*8  AKING, BKING, CKING
        PARAMETER	( ABATES = 3.9992662D-04 )
        PARAMETER	( BBATES = 1.0689770D-02 )
        PARAMETER	( CBATES = 6.6814090D-05 )
        PARAMETER	( AKING = 1.0469541 )
        PARAMETER	( BKING = 3.2503153D-04 )
        PARAMETER	( CKING = 3.8622851D-05 )

C  help variables

	INTEGER		 I, L, UI, J, N, M, G
	DOUBLE PRECISION
     &       SIG, SIG2, SIG4, RAYCOF, FKING, DEPOL, ABSCO_ACTIVEGAS,
     &       THICK, EXTCO, OTHERGAS_CONC,
     &       ACTIVEGAS, OTHERGAS_VOD(MAX_OTHERGAS),
     &       AIR_DENSITY, TEMP, ETA, X00, ABSGAS,
     &       SCACO_R, ABSCO_R, EXTCO_R, SIGMA
	REAL*8 TAU, DELT, DELTM, DELTM1
	REAL*8 EXTNCOF(MAXLAYER)
	logical		 aermask(maxlayer)

C  Fix inputs

	DO_LAYER_LINEARIZATION = .TRUE.
	DO_ALBEDO_LINEARIZATION = .TRUE.
	DO_SIMULATION_ONLY = .FALSE.

C  Aerosols
C  ---------------------
	DO I = 1, NGRID
	  EXTCO_a(i) = aerextc(i)
	  scaco_a(i) = aerscat(i)
	  aermask(i) = .true.
	  if (scaco_a(i) .lt. 1.0d-05 ) aermask(i) = .false.
	ENDDO

C  Rayleigh scattering coefficient + depolarization
C  ------------------------------------------------

C  Rayleigh scattering coefficient

	SIG = 1.0D+3 / LAMDA
	SIG2 = SIG * SIG
	SIG4 = SIG2 * SIG2
	RAYCOF = ABATES * SIG4 / ( ONE - BBATES*SIG2 - CBATES*SIG4 )
	RAYCOF = RAYCOF * 1.0D-24

C  Derivation of depolarisation factor D from King factors for air, as fitted
C  by K. Chance to data of Bates (1984).  Wavelength LAMDA must be in NM.
C        Fking = ( 6 + 3.Depol ) / ( 6 - 7.Depol ) 

        SIG = 1.0D+03 / LAMDA
        SIG2 = SIG * SIG
        SIG4 = SIG2 * SIG2
        FKING = AKING + SIG2 * BKING + SIG4 * CKING
        DEPOL = 6.0D0 * (FKING-1.0D0) / ( 3.0D0+7.0D0*FKING)

C  Set up optical depth, single scatter albedos and phase function moments
C  -----------------------------------------------------------------------

C  initialise

	TAUGRID_INPUT(0) = ZERO
	DO G = 1, N_OTHERGAS
	  OTHERGAS_VOD(G) = ZERO
	ENDDO
	ACTIVEGAS_OPDEP   = ZERO
	ACTIVEGAS_COLUMN  = ZERO
	ACTIVEGAS_EFFXSEC = ZERO

C  start layer loop

	DO I = 1, NLAYERS

C  air density

	  THICK = DELTA(I)
	  TEMP = TEMPERATURE(I)
	  AIR_DENSITY = PRESSURE(I)*FP/TEMP*SMKM

C  Othergas absorption 

	  ABSCO_R = ZERO
	  DO G = 1, N_OTHERGAS
            OTHERGAS_CONC   = OTHERGAS_VMRS(I,G) * AIR_DENSITY
	    IF ( OTHERGAS_NAMES(G) .EQ. 'OZON' ) THEN
              IF ( ACTIVEGAS_NAME .EQ. 'HCHO' ) THEN 
                CALL GET_O3XSEC ( LAMDA, TEMP, SIGMA, ETA )
              ELSE
	        CALL GET_O3XSEC_NO2 ( LAMDA, SIGMA )
	      ENDIF
            ELSE IF ((OTHERGAS_NAMES(G) .EQ. 'NO2') .AND. 
     &         (ACTIVEGAS_NAME .EQ. 'HCHO')) THEN
	       CALL GET_NO2XSEC ( LAMDA, TEMP, SIGMA )
	    ENDIF
	    ABSGAS = SIGMA * OTHERGAS_CONC
	    ABSCO_R = ABSCO_R + ABSGAS
	    OTHERGAS_VOD(G)  = OTHERGAS_VOD(G)  + ABSGAS  * THICK
	  ENDDO

C  ACTIVEGAS absorption

          ACTIVEGAS = ACTIVEGAS_CONC(I)*SMKM
          IF ( ACTIVEGAS_NAME .EQ. 'HCHO' ) THEN
            CALL GET_HCHOXSEC ( LAMDA, SIGMA )
          ENDIF
	  IF ( ACTIVEGAS_NAME .EQ. 'NO2 ' ) THEN
	    CALL GET_NO2XSEC_NO2 ( LAMDA, SIGMA )
	  ENDIF
	  ABSCO_ACTIVEGAS = SIGMA * ACTIVEGAS
	  ABSCO_R = ABSCO_R + ABSCO_ACTIVEGAS 
	  ACTIVEGAS_OPDEP = ACTIVEGAS_OPDEP + ABSCO_ACTIVEGAS * THICK
	  ACTIVEGAS_COLUMN = ACTIVEGAS_COLUMN + ACTIVEGAS * THICK

C  Molecular scattering and extinction

	  SCACO_R = RAYCOF * AIR_DENSITY
	  EXTCO_R = ABSCO_R + SCACO_R

C  add aerosol extinction

	  if ( aermask(i) ) then
	    EXTCO = EXTCO_R + EXTCO_A(I)
	  else
	    EXTCO = EXTCO_R
	  endif

C  save total extinction for later use in SPHERICAL PREPARE

	  EXTNCOF(I) = EXTCO

C  Aerosol single scatter albedo

	  if ( aermask(i) ) then
	    OMEGA_TOTAL_INPUT(I)  = (SCACO_A(I) + SCACO_R) / EXTCO
	  else
	    OMEGA_TOTAL_INPUT(I) = SCACO_R / EXTCO  !Rayleigh
	  endif
          IF (OMEGA_TOTAL_INPUT(I) .GT. 0.9999D0) 
     &      OMEGA_TOTAL_INPUT(I) = 0.9999D0

C  cumulative vertical optical depth

	  TAUGRID_INPUT(I) = TAUGRID_INPUT(I-1) + EXTCO*THICK

C  rayleigh phase function moments

	  PHASMOMS_TOTAL_INPUT(0,I) = ONE
	  PHASMOMS_INPUT_R(1,I) = ZERO
          PHASMOMS_INPUT_R(2,I) = (ONE-DEPOL)/(TWO+DEPOL)
	  DO L = 3, NMOMENTS
	    PHASMOMS_INPUT_R(L,I) = ZERO
	  ENDDO

C  aerosol phse function moments

	  if ( aermask(i) ) then
	    DO L = 1, NMOMENTS
	      PHASMOMS_TOTAL_INPUT(L,I) = 
     &        PHASEFCN(I,L)*SCACO_A(I)/(OMEGA_TOTAL_INPUT(I)*EXTCO) + 
     &        PHASMOMS_INPUT_R(L,I)*SCACO_R/(OMEGA_TOTAL_INPUT(I)*EXTCO)
	    ENDDO
	  else
	    DO L = 1, NMOMENTS
	      PHASMOMS_TOTAL_INPUT(L,I) = PHASMOMS_INPUT_R(L,I)
	    ENDDO
	  endif

C  linearization input - only with respect to ACTIVEGAS volume mixing ratio

	  IF ( DO_LAYER_LINEARIZATION ) THEN
	    LAYER_VARY_FLAG(I) = .TRUE.
	    LAYER_VARY_NUMBER(I) = 1
	    OMEGA_VARS_TOTAL_INPUT(1,I) = -ABSCO_ACTIVEGAS / EXTCO
            EXT_VARS_INPUT(1,I) = ABSCO_ACTIVEGAS / EXTCO
	  ELSE
	    LAYER_VARY_FLAG(I) = .FALSE.
	    LAYER_VARY_NUMBER(I) = 0
	    OMEGA_VARS_TOTAL_INPUT(1,I) = ZERO
            EXT_VARS_INPUT(1,I) = ZERO
	  ENDIF
	ENDDO

C  Quasi-spherical input

	NREFRAC = ZERO
	REARTH = 6371.0D0
	ZS_GRID(0) = ZGRID(1) + DELTA(1)*0.5
	DO J = 1, NLAYERS
	  ZS_GRID(J) = ZS_GRID(J-1) - DELTA(J)
	  FINEGRID(J) = MAXFINELAYER
	ENDDO

C  Prepare spherical attenuation (shell geometry)

	IF ( DO_QUASPHER_BEAM ) THEN

	  X00 = DCOS ( SUN0 * DEG_TO_RAD )
	  CALL SPHERICAL_PREPARE
     I     ( MAXLAYER, MAXFINELAYER, NLAYERS, FINEGRID,
     I       X00, REARTH, NREFRAC, DO_QSREFRAC_BEAM,
     I       ZS_GRID, TAUGRID_INPUT, PS_GRID, TS_GRID, EXTNCOF,
     O       TAUPATH_FINE_INPUT, MUPATH_FINE,
     O       TAUTHICK_INPUT, TAUTHICK_FINE_INPUT )

C  Prepare spherical attenuation (plane-parallel mock-up)

	ELSE

	  X00 = DCOS ( SUN0 * DEG_TO_RAD )
	  DO N = 1, NLAYERS
	    FINEGRID(N) = MAXFINELAYER
	    DELT = TAUGRID_INPUT(N) - TAUGRID_INPUT(N-1)
	    DELTM = DELT / DFLOAT(FINEGRID(N))
	    TAUTHICK_INPUT(N,N) = DELT/X00
	    DO J = 0, FINEGRID(N)
	      MUPATH_FINE(N,J) = X00
	      TAU = TAUGRID_INPUT(N-1) + DELTM*J
	      TAUTHICK_FINE_INPUT(N,N,J) = DELTM*J / X00
	      TAUPATH_FINE_INPUT(N,J) = TAU / X00
	      DO M = 1, N-1
	        DELTM1 = TAUGRID_INPUT(M) - TAUGRID_INPUT(M-1)
	        TAUTHICK_FINE_INPUT(N,M,J) = DELTM1 / X00
	        IF ( J.EQ.FINEGRID(N) ) TAUTHICK_INPUT(N,M) = DELTM1/X00
	      ENDDO
	    ENDDO
	  ENDDO

	ENDIF

C  prepare surface (Lambertian albedo)

	IF ( DO_LAMBERTIAN_ALBEDO ) THEN
	  DO I = 1, NSTREAMS
	    BIREFLEC_0(0,I) = ONE
	    EMISSIVITY(I) = ONE - ALBEDO
	    DO J = 1, NSTREAMS
	      BIREFLEC(0,I,J) = ONE
	    ENDDO
	  ENDDO
	  DO UI = 1, N_USER_STREAMS
	    USER_BIREFLEC_0(0,UI) = ONE
	    USER_EMISSIVITY(UI) = ONE - ALBEDO
	    DO J = 1, NSTREAMS
	      USER_BIREFLEC(0,UI,J) = ONE
	    ENDDO
	  ENDDO
	ENDIF

C  effective cross-section

	ACTIVEGAS_EFFXSEC = ACTIVEGAS_OPDEP / ACTIVEGAS_COLUMN

C  finish

	RETURN
	END
C

	SUBROUTINE SPHERICAL_PREPARE
     I     ( MAXLAYERS, MAXFINELEVELS, NLAYERS, FINEGRID,
     I       MU_TOA, REARTH, NREFRAC, DO_REFRAC,
     I       Z_GRID, TAU_GRID, P_GRID, T_GRID, EXT,
     O       TAUPATH_FINE, MUPATH_FINE,
     O       TAUTHICK_INPUT, TAUTHICK_FINE_INPUT )

C  Generate optical depths TAUPATH_FINE and zenith cosines MUPATH_FINE
C  for a curved ray-traced beam through a multilayer atmosphere.

C  Each layer is assumed to have constant optical properties, and the
C  layer gridding is input to the module. Values of Z_GRID etc. are
C  given at the layer boundaries, with the first value (index 0) at TOA.

C  The refractive geometry is assumed to start at the TOA level, where the
C  incoming beam zenith is X00. We also require the earth radius and the
C  refractive index constant NREFRAC at standard T&P (Born-Wolf approxmn.).

C  There is no refraction if the flag DO_REFRAC is not set. In this case
C  we do not require the pressue and temperature information. The calculation
C  will then run in a spehrical shell atmosphere ('geometric rays'). In
C  the refracting case, temperature will be interpolated linearlyl with height,
C  and pressure log-linearly. The refraction uses a simple Snell's law rule.

C  Output is specified at layer boundaries and at a number of intermediate
C  levels specified by the finelayer controlling array FINEGRID. This
C  assumes that each layer is divided into a number of equidistant (both
C  height and vertical optical depth!) points, which will vary from layer
C  to layer.

C  Module is stand-alone.

C  Input arguments
C  ===============

C  dimensioning

	INTEGER		 MAXLAYERS, MAXFINELEVELS

C  layer control

	INTEGER		 NLAYERS
	INTEGER		 FINEGRID(MAXLAYERS)

C  Refractive and geometry control inputs

	REAL*8 MU_TOA
	REAL*8 REARTH
	REAL*8 NREFRAC
	LOGICAL		 DO_REFRAC

C  input gridding

	REAL*8 Z_GRID(0:MAXLAYERS)
	REAL*8 P_GRID(0:MAXLAYERS)
	REAL*8 T_GRID(0:MAXLAYERS)
	REAL*8 TAU_GRID(0:MAXLAYERS)
	REAL*8 EXT(MAXLAYERS)

C  Output arguments
C  ================

	REAL*8 TAUPATH_FINE(MAXLAYERS,0:MAXFINELEVELS)
	REAL*8 MUPATH_FINE(MAXLAYERS,0:MAXFINELEVELS)

	REAL*8 TAUTHICK_INPUT(MAXLAYERS,MAXLAYERS)
	REAL*8
     &      TAUTHICK_FINE_INPUT(MAXLAYERS,MAXLAYERS,0:MAXFINELEVELS)

C  Local variables

	INTEGER		 N, J, NFINE, M, MAXDIM
	PARAMETER	 ( MAXDIM = 500 )

	REAL*8 GM_TOA, TH_TOA, TH0, TH1, GM0, GM1
	REAL*8 H(0:MAXDIM), DELZ(MAXDIM), TAUP
	REAL*8 Z_DIFF, Z_0, Z, CONST0
	REAL*8 X, XD, CUMDEP, S2, HF

C  Geometrical Case
C  ================

	GM_TOA = DSQRT ( 1.0D0 - MU_TOA * MU_TOA )
	TH_TOA = DASIN(GM_TOA)

	DO N = 0, NLAYERS
	  H(N) = Z_GRID(N) + REARTH
	ENDDO
	DO N = 1, NLAYERS
	  DELZ(N) = Z_GRID(N-1)-Z_GRID(N)
	ENDDO
	CONST0 = GM_TOA / H(0)
	CUMDEP = 0.0D0

C  zero output

	DO N = 1, NLAYERS
	  DO J = 0, FINEGRID(N)
	    TAUPATH_FINE(N,J) = 0.0D0
	  ENDDO
	ENDDO
	TAUPATH_FINE(1,0) = 0.0D0
	MUPATH_FINE(1,0)  = MU_TOA

C  zero first layer

	TAUTHICK_FINE_INPUT(1,1,0) = 0.0D0

	DO N = 1, NLAYERS
	  NFINE      = FINEGRID(N)
	  Z_DIFF   = DELZ(N) / DFLOAT (NFINE)
	  Z_0 = Z_GRID(N-1)
	  Z  = Z_0
	  DO M = 1, N-1
	    TAUTHICK_FINE_INPUT(N,M,0) =
     &            TAUTHICK_FINE_INPUT(N-1,N-1,NFINE)
	  ENDDO
	  DO J = 1, NFINE
	    Z = Z - Z_DIFF
	    X = Z_0 - Z
	    XD = X + CUMDEP
	    HF = H(0) - XD
	    GM0 = CONST0 * HF	  
	    TH0 = DASIN ( GM0 )
	    TAUP = 0.0D0
	    DO M = 1, N-1
	      GM1 = H(M-1) * GM0 / H(M)
	      TH1 = DASIN(GM1)
	      S2 = H(M-1) * DSIN(TH1-TH0) / GM1
	      TAUTHICK_FINE_INPUT(N,M,J) = EXT(M) * S2
	      IF ( J.EQ.NFINE ) TAUTHICK_INPUT(N,M) = EXT(M) * S2
	      TAUP = TAUP + TAUTHICK_FINE_INPUT(N,M,J)
	      TH0 = TH1
	      GM0 = GM1        
	    ENDDO
	    S2 = H(N-1)*DSIN(TH_TOA-TH0) / GM_TOA
	    TAUP = TAUP + EXT(N) * S2
	    TAUPATH_FINE(N,J) = TAUP
	    MUPATH_FINE(N,J)  = MU_TOA
	    TAUTHICK_FINE_INPUT(N,N,J) = EXT(N) * S2
	    IF ( J.EQ.NFINE ) TAUTHICK_INPUT(N,N) = EXT(N) * S2
	  ENDDO
	  CUMDEP = XD
	  IF (N.GT.1) THEN
	    TAUPATH_FINE(N,0) = TAUPATH_FINE(N-1,FINEGRID(N-1))
	    MUPATH_FINE(N,0)  = MU_TOA
	    DO M = 1, N-1
	      TAUTHICK_FINE_INPUT(N,M,0) = 
     &             TAUTHICK_FINE_INPUT(N-1,M,FINEGRID(M))
	    ENDDO
	  ENDIF

	ENDDO

C  Refractive geometry case
C  ========================

	IF ( DO_REFRAC ) THEN
	  NREFRAC = 0.0D0
	  P_GRID(1) = 0.0D0
	  T_GRID(1) = 0.0D0
	  PAUSE'placeholder on refraction'
	ENDIF

C  Finish

	END

