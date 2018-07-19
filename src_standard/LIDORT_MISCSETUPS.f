C ###############################################################
C #							    	#
C #                    THE LIDORT  MODEL			#
C #							        #
C #      (LInearized Discrete Ordinate Radiative Transfer)      #
C #	  --	     -	      -        - 	 -	        #
C #		  					        #
C ###############################################################

C ###############################################################
C #		  					        #
C #  Author :	Robert. J. D. Spurr			        #
C #		  					        #
C #  Address :	Harvard-Smithsonian Center for Astrophysics     #
C #		60 Garden Street			        #
C #	 	Cambridge, MA 02138, USA			#
C #		Tel: (617) 496 7819				#
C #		  					        #
C #  Email :      rspurr@cfa.harvard.edu			#
C #		  					        #
C #  Version :	  2.3					        #
C #  Release Date   January 2001				#
C #		  					        #
C ###############################################################

	SUBROUTINE LIDORT_MISCSETUPS

C  miscellaneous setup operations, Master routine

C  Delta-m scaling of input quantities

	CALL LIDORT_DELTAMSCALE_TOTAL

C  initialise single scatter albedo terms

	CALL LIDORT_SSALBINIT_TOTAL

C  Prepare quasi-spherical attenuation

	CALL LIDORT_QSPREP

C  Transmittances and Transmittance factors

	CALL LIDORT_PREPTRANS

C  Finish

	END  

C

	SUBROUTINE LIDORT_DELTAMSCALE_TOTAL

C  include files
C  -------------

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of set up variables (output to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  local variables

	DOUBLE PRECISION FDEL, FAC2, DNL1
	DOUBLE PRECISION TAU, TAUFAC, DNM1
	INTEGER		 N, N1, L, UT, UTA, NM1, K

C  Compute the layer optical thickness values for the input grid

	DO N = 1, NLAYERS
	  N1 = N - 1
	  DELTA(N) = TAUGRID_INPUT(N) - TAUGRID_INPUT(N1)
	ENDDO

C  DELTAM SCALING
C  ==============

	IF ( DO_DELTAM_SCALING ) THEN

	  TAUGRID(0) = ZERO
	  NM1 = NMOMENTS+1
	  DNM1 = DFLOAT(2*NM1+1)

C  Scaling for layer input
C  -----------------------

	  DO N = 1, NLAYERS

	    N1 = N - 1

C  overall truncation factor

	    FDEL = PHASMOMS_TOTAL_INPUT(NM1,N) / DNM1
	    TRUNC_FACTOR(N) = FDEL
	    FAC2 = ONE - FDEL
	    FAC1(N) = ONE - FDEL * OMEGA_TOTAL_INPUT(N)

C  Scale single-scatter albedos and phase function moments

	    PHASMOMS_TOTAL(0,N) = ONE
	    DO L = 1, NMOMENTS
	      DNL1 = DFLOAT(2*L + 1 )
	      PHASMOMS_TOTAL(L,N) =
     &            ( PHASMOMS_TOTAL_INPUT(L,N) - FDEL*DNL1 ) / FAC2
	    ENDDO
	    OMEGA_TOTAL(N) = OMEGA_TOTAL_INPUT(N) * FAC2 / FAC1(N)

C  scale optical depth grid

	    DELTA(N) = DELTA(N) * FAC1(N)
	    TAUGRID(N) = TAUGRID(N1) + DELTA(N)

	  ENDDO

C  Scaling for user-defined off-grid optical depths
C  ------------------------------------------------

C  off-grid values (on-grid values have already been scaled)

	  IF ( DO_USER_TAUS ) THEN
	    UT = 0
	    DO UTA = 1, N_OUT_USERTAUS
	      IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
	        UT = UT + 1
	        N = OFFGRID_UTAU_LAYERIDX(UT)
	        TAU = USER_TAUS_INPUT(UTA) - TAUGRID_INPUT(N-1)
	        TAUFAC = TAU * FAC1(N)
	        OFFGRID_UTAU_VALUES(UT) = TAUFAC
	      ENDIF
	    ENDDO
	  ENDIF

C  Scale layer path thickness values

	  DO N = 1, NLAYERS
	    DO K = 1, N
	      TAUTHICK(N,K) = TAUTHICK_INPUT(N,K) * FAC1(K)
	    ENDDO
	  ENDDO

C  NO DELTAM SCALING
C  =================

C  move input geophysical variables to Workspace quantities

	ELSE

	  TAUGRID(0) = ZERO
	  DO N = 1, NLAYERS
	    TRUNC_FACTOR(N) = ZERO
	    OMEGA_TOTAL(N) = OMEGA_TOTAL_INPUT(N)
	    DO L = 0, NMOMENTS
	      PHASMOMS_TOTAL(L,N) = PHASMOMS_TOTAL_INPUT(L,N)
	    ENDDO
	    TAUGRID(N) = TAUGRID_INPUT(N)
	  ENDDO

	  IF ( DO_USER_TAUS ) THEN
	    UT = 0
	    DO UTA = 1, N_OUT_USERTAUS
	      IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
	        UT = UT + 1
	        N = OFFGRID_UTAU_LAYERIDX(UT)
	        TAU = USER_TAUS_INPUT(UTA) - TAUGRID_INPUT(N-1)
	        OFFGRID_UTAU_VALUES(UT) = TAU
	      ENDIF
	    ENDDO
	  ENDIF

	  DO N = 1, NLAYERS
	    DO K = 1, N
	      TAUTHICK(N,K) = TAUTHICK_INPUT(N,K)
	    ENDDO
	  ENDDO

	ENDIF

C  Finish module

	END

C
    
	SUBROUTINE LIDORT_SSALBINIT_TOTAL

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of set up variables (output to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  local variables

	INTEGER		 N, L

C  phase moment-weighted OMEGA

 	DO N = 1, NLAYERS
	  DO L = 0, NMOMENTS
	    OMEGA_MOMS(N,L) = OMEGA_TOTAL(N)*PHASMOMS_TOTAL(L,N)
	  ENDDO
	ENDDO

C  Finish

	END

C

	SUBROUTINE LIDORT_QSPREP

C  include files
C  -------------

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of set up variables (local output to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  Local variables
C  ---------------

	INTEGER		 N, K
	DOUBLE PRECISION S_T_0, S_T_1, SEC0, TAU

C  plane-parallel case
C  -------------------

	IF ( .NOT. DO_QUASPHER_BEAM ) THEN

C  ## new code

	  SEC0 = ONE / X0
	  LAYER_PIS_CUTOFF = NLAYERS
	  DO N = 1, NLAYERS
	    IF  (N.LE.LAYER_PIS_CUTOFF ) THEN
	      IF ( TAUGRID(N) * SEC0 .GT. MAX_TAU_SPATH ) THEN
	        LAYER_PIS_CUTOFF = N
	      ENDIF
	      AVERAGE_SECANT(N) = SEC0
	      INITIAL_TRANS(N) = DEXP ( - TAUGRID(N-1) * SEC0 )
	      LOCAL_SZA(N) = X0
	    ELSE
	      AVERAGE_SECANT(N) = ZERO
	      INITIAL_TRANS(N)  = ZERO
	      LOCAL_SZA(N) = ZERO
	    ENDIF
	  ENDDO

C  ## old code
C	  SEC0 = ONE / X0
C	  DO N = 1, NLAYERS
C	    AVERAGE_SECANT(N) = SEC0
C	    INITIAL_TRANS(N) = DEXP ( - TAUGRID(N-1) * SEC0 )
C	    LOCAL_SZA(N) = X0
C	  ENDDO

	ELSE

C  pseudo-spherical case
C  ---------------------

C  Get the total spherical attenuation from layer thickness sums

	  TAUSLANT(0) = ZERO
	  DO N = 1, NLAYERS
	    TAU = ZERO
	    DO K = 1, N
	      TAU = TAU + TAUTHICK(N,K)
	    ENDDO
	    TAUSLANT(N) = TAU
	  ENDDO

C  set up the average secant formulation


C  ## new code

	  S_T_0 = ONE
	  LAYER_PIS_CUTOFF = NLAYERS
	  DO N = 1, NLAYERS
	    IF  (N.LE.LAYER_PIS_CUTOFF ) THEN
	      IF ( TAUSLANT(N) .GT. MAX_TAU_SPATH ) THEN
	        LAYER_PIS_CUTOFF = N
	      ELSE
	        S_T_1 = DEXP ( - TAUSLANT(N) )
	      ENDIF
	      AVERAGE_SECANT(N) = (TAUSLANT(N)-TAUSLANT(N-1))/DELTA(N)
	      INITIAL_TRANS(N)  = S_T_0
	      LOCAL_SZA(N) = X0
	      S_T_0 = S_T_1
	    ELSE
	      AVERAGE_SECANT(N) = ZERO
	      INITIAL_TRANS(N)  = ZERO
	    ENDIF
	  ENDDO

C  Set the Local solar zenith angles
C  Distinguish between the refractive and non-refractive cases.

	  IF ( DO_QSREFRAC_BEAM ) THEN
	    DO N = 1, NLAYERS
	      IF  (N.LE.LAYER_PIS_CUTOFF ) THEN
   	        LOCAL_SZA(N) = SUN_SZA_COSINES(N)
	      ELSE
 	        LOCAL_SZA(N) = ZERO
	      ENDIF
	    ENDDO
	  ELSE
	    DO N = 1, NLAYERS
	      IF  (N.LE.LAYER_PIS_CUTOFF ) THEN
   	        LOCAL_SZA(N) = X0
	      ELSE
 	        LOCAL_SZA(N) = ZERO
	      ENDIF
	    ENDDO
	  ENDIF


C  ## old code

C	  S_T_0 = ONE
C	  DO N = 1, NLAYERS
C	    S_T_1 = DEXP ( - TAUSLANT(N) )
C	    MB = DELTA(N) / DLOG ( S_T_0 / S_T_1 )
C	    SECBAR = ONE / MB
C	    AVERAGE_SECANT(N) = SECBAR
C	    INITIAL_TRANS(N)  = S_T_0
C	    LOCAL_SZA(N) = X0
C	    S_T_0 = S_T_1
C	  ENDDO

	ENDIF

C  finish

	END

C

	SUBROUTINE LIDORT_PREPTRANS

C  Prepare transmittances and transmittance factors

C  include files
C  -------------

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of setup variables (output stored here)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include file of Multiplier variables (output stored here)

	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  local variables
C  ---------------

	INTEGER		 N, UT, UM
	DOUBLE PRECISION WDEL, WX, WUDEL, UDEL, SPHER
	DOUBLE PRECISION XT, DIFF, SB, SECMUM, SU, SD
	DOUBLE PRECISION UX_DN, UX_UP, WDEL_UXUP

C  Transmittance factors for average secant stream
C  ===============================================

C  Whole layer Transmittance factors
C  ---------------------------------

C  layer transmittance 

	DO N = 1, NLAYERS
	  IF ( N. GT. LAYER_PIS_CUTOFF ) THEN
	    T_DELT_MUBAR(N) = ZERO
	  ELSE
	    SPHER = DELTA(N) * AVERAGE_SECANT(N)
	    IF ( SPHER .GT. MAX_TAU_SPATH ) THEN
	      T_DELT_MUBAR(N) = ZERO
	    ELSE
	      T_DELT_MUBAR(N) = DEXP ( - SPHER )
	    ENDIF
	  ENDIF
	ENDDO

C  Partial layer transmittance factors (for off-grid optical depths)
C  -----------------------------------------------------------------

	DO UT = 1, N_OFFGRID_USERTAUS
	  N     = OFFGRID_UTAU_LAYERIDX(UT)
	  IF ( N. GT. LAYER_PIS_CUTOFF ) THEN
	    T_UTDN_MUBAR(UT) = ZERO
	    T_UTUP_MUBAR(UT) = ZERO
	  ELSE
	    XT = OFFGRID_UTAU_VALUES(UT)
	    SPHER = XT * AVERAGE_SECANT(N)
	    IF ( SPHER .GT. MAX_TAU_SPATH ) THEN
	      T_UTDN_MUBAR(UT) = ZERO
	    ELSE
	      T_UTDN_MUBAR(UT) = DEXP ( - SPHER )
	    ENDIF
	    SPHER = ( DELTA(N) - XT ) * AVERAGE_SECANT(N)
	    IF ( SPHER .GT. MAX_TAU_SPATH ) THEN
	      T_UTUP_MUBAR(UT) = ZERO
	    ELSE
	      T_UTUP_MUBAR(UT) = DEXP ( - SPHER )
	    ENDIF
	  ENDIF
	ENDDO

C  Transmittances for User Streams
C  ===============================

C  return if not flagged

	IF ( .NOT. DO_USER_STREAMS ) RETURN

C  Initial transmittances divided by user streams
C  ----------------------------------------------

	DO N = 1, NLAYERS
	  DO UM = 1, N_USER_STREAMS
	    ITRANS_USERM(N,UM) = INITIAL_TRANS(N) * USER_SECANTS(UM)
	  ENDDO
	ENDDO

C  Whole Layer transmittances
C  --------------------------

	DO N = 1, NLAYERS
	  IF ( STERM_LAYERMASK_UP(N).OR.STERM_LAYERMASK_DN(N) ) THEN
	    DO UM = 1, N_USER_STREAMS
	      SPHER = DELTA(N) * USER_SECANTS(UM)
	      IF ( SPHER.GT.MAX_TAU_UPATH ) THEN
	        T_DELT_USERM(N,UM) = ZERO
	      ELSE
	        T_DELT_USERM(N,UM) = DEXP ( - SPHER )
	      ENDIF
	    ENDDO
	  ENDIF
	ENDDO

C  ## old code
C	DO N = 1, NLAYERS
C	  IF ( STERM_LAYERMASK_UP(N).OR.STERM_LAYERMASK_DN(N) ) THEN
C	    DO UM = 1, N_USER_STREAMS
C	      SPHER = DELTA(N) * USER_SECANTS(UM)
C	      HELP = DEXP ( - SPHER )
C	      T_DELT_USERM(N,UM) = HELP
C	    ENDDO
C	  ENDIF
C	ENDDO

C  Partial Layer transmittances for off-grid optical depths
C  --------------------------------------------------------

	DO UT = 1, N_OFFGRID_USERTAUS
	  N    = OFFGRID_UTAU_LAYERIDX(UT)
	  XT = OFFGRID_UTAU_VALUES(UT)
	  DO UM = 1, N_USER_STREAMS
	    SPHER = XT * USER_SECANTS(UM)
	    IF ( SPHER .GT. MAX_TAU_UPATH ) THEN
	      T_UTDN_USERM(UT,UM) = ZERO
	    ELSE
	      T_UTDN_USERM(UT,UM) = DEXP ( - SPHER )
	    ENDIF
	    SPHER = ( DELTA(N) - XT ) * USER_SECANTS(UM)
	    IF ( SPHER .GT. MAX_TAU_UPATH ) THEN
	      T_UTUP_USERM(UT,UM) = ZERO
	    ELSE
	      T_UTUP_USERM(UT,UM) = DEXP ( - SPHER )
	    ENDIF
	  ENDDO
	ENDDO

C  ## old code
C	DO UT = 1, N_OFFGRID_USERTAUS
C	  N    = OFFGRID_UTAU_LAYERIDX(UT)
C	  DO UM = 1, N_USER_STREAMS
C	    SPHER = OFFGRID_UTAU_VALUES(UT) * USER_SECANTS(UM)
C	    HELP = DEXP ( - SPHER )
C	    T_UTDN_USERM(UT,UM) = HELP
C	    T_UTUP_USERM(UT,UM) = T_DELT_USERM(N,UM) / HELP
C	  ENDDO
C	ENDDO

C  Transmittance Factors for User Streams
C  ======================================

C  L'Hopital's Rule flags for Downwelling EMULT
C  --------------------------------------------

	IF ( DO_DNWELLING ) THEN
	  DO N = 1, NLAYERS
	    SB = AVERAGE_SECANT(N)
	    DO UM = 1, N_USER_STREAMS
	      DIFF = DABS ( USER_SECANTS(UM) - SB )
	      IF ( DIFF .LT. HOPITAL_TOLERANCE ) THEN
	        EMULT_HOPRULE(N,UM) = .TRUE.
	      ELSE
	        EMULT_HOPRULE(N,UM) = .FALSE.
	      ENDIF
	    ENDDO
	  ENDDO
	ENDIF

C  sigma functions (all layers)
C  ----------------------------

	DO N = 1, NLAYERS
	  SB = AVERAGE_SECANT(N)
	  DO UM = 1, N_USER_STREAMS
	    SECMUM = USER_SECANTS(UM)
	    SIGMA_P(N,UM) = SB + SECMUM
	    SIGMA_M(N,UM) = SB - SECMUM
	  ENDDO
	ENDDO

C  upwelling External source function multipliers
C  ----------------------------------------------

	IF ( DO_UPWELLING ) THEN

C  whole layer

	  DO N = 1, NLAYERS
	    IF ( STERM_LAYERMASK_UP(N) ) THEN
	      IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	        DO UM = 1, N_USER_STREAMS
	          EMULT_UP(UM,N) = ZERO
	        ENDDO
	      ELSE
	        WDEL = T_DELT_MUBAR(N)
	        DO UM = 1, N_USER_STREAMS
	          WUDEL = WDEL * T_DELT_USERM(N,UM)
	          SU = ( ONE - WUDEL ) / SIGMA_P(N,UM)
	          EMULT_UP(UM,N) = ITRANS_USERM(N,UM) * SU
	        ENDDO
	      ENDIF
	    ENDIF
	  ENDDO

C  Partial layer

	  DO UT = 1, N_OFFGRID_USERTAUS
	    N  = OFFGRID_UTAU_LAYERIDX(UT)
	    IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	      DO UM = 1, N_USER_STREAMS
	        UT_EMULT_UP(UM,UT) = ZERO
	      ENDDO
	    ELSE
	      WX   = T_UTDN_MUBAR(UT)
	      WDEL = T_DELT_MUBAR(N)
	      DO UM = 1, N_USER_STREAMS
	        UX_UP = T_UTUP_USERM(UT,UM)
	        WDEL_UXUP = UX_UP * WDEL
	        SU = ( WX - WDEL_UXUP ) / SIGMA_P(N,UM)
	        UT_EMULT_UP(UM,UT) = ITRANS_USERM(N,UM) * SU
	      ENDDO
	    ENDIF
	  ENDDO

	ENDIF

C  downwelling External source function multipliers
C  ------------------------------------------------

C    .. Note use of L'Hopitals Rule

	IF ( DO_DNWELLING ) THEN

C  whole layer

	  DO N = 1, NLAYERS
	    IF ( STERM_LAYERMASK_DN(N) ) THEN
	      IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	        DO UM = 1, N_USER_STREAMS
	          EMULT_DN(UM,N) = ZERO
	        ENDDO
	      ELSE
	        WDEL = T_DELT_MUBAR(N)
	        DO UM = 1, N_USER_STREAMS
	          UDEL = T_DELT_USERM(N,UM)
	          IF ( EMULT_HOPRULE(N,UM) ) THEN
	            SD = DELTA(N) * UDEL
	          ELSE
	            SD = ( UDEL - WDEL ) / SIGMA_M(N,UM)
	          ENDIF 
	          EMULT_DN(UM,N) = ITRANS_USERM(N,UM) * SD
	        ENDDO
	      ENDIF
	    ENDIF
	  ENDDO

C  Partial layer

	  DO UT = 1, N_OFFGRID_USERTAUS
	    N  = OFFGRID_UTAU_LAYERIDX(UT)
	    IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	      DO UM = 1, N_USER_STREAMS
	        UT_EMULT_DN(UM,UT) = ZERO
	      ENDDO
	    ELSE
	      WX   = T_UTDN_MUBAR(UT)
	      DO UM = 1, N_USER_STREAMS
	        UX_DN = T_UTDN_USERM(UT,UM)
	        IF ( EMULT_HOPRULE(N,UM) ) THEN
	          SD = OFFGRID_UTAU_VALUES(UT) * UX_DN
	        ELSE
	          SD = ( UX_DN - WX ) / SIGMA_M(N,UM)
	        ENDIF
	        UT_EMULT_DN(UM,UT) = ITRANS_USERM(N,UM) * SD
	      ENDDO
	    ENDIF
	  ENDDO

	ENDIF

C  Finish

	END

C

	SUBROUTINE CHAPMAN_FUNCTION

C  This is the Chapman function calculation of the slant path
C  optical thickness array TAUTHICK_INPUT.

C  This module calculates TAUTHICK_INPUT internally inside LIDORT
C  saving you the job of doing it yourself.

C  You must specify the Earth_radius and the height grid in order
C  to make this work.

C  This is a straightforward geometrical calculation, and is only
C  valid for a NON-REFRACTIVE atmosphere. If you want refraction,
C  you must set TAUTHICK_INPUT yourself outsied of LIDORT!!

C  The non-refractive condition will be checked before this module
C  is called


C  Input files

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  Geophysical variables (output stored here)

	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'

C  Local variables

	INTEGER		 N, M
	DOUBLE PRECISION GM_TOA, MU_TOA, HELP1, HELP2
	DOUBLE PRECISION H(0:MAXLAYER), EXT(MAXLAYER) 
	DOUBLE PRECISION STH, CTH, DEL, DELT, DELS, S1, S0

C  get spherical optical depths
C  ----------------------------

C  Prepare spherical attenuation (shell geometry)

	IF ( DO_QUASPHER_BEAM ) THEN

	  MU_TOA = DCOS ( SUN0 * DEG_TO_RAD )
	  GM_TOA = DSQRT ( 1.0D0 - MU_TOA * MU_TOA )
	  DO N = 0, NLAYERS
	    H(N) = HEIGHT_GRID(N) + EARTH_RADIUS
	  ENDDO
	  DO N = 1, NLAYERS
	    DEL = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
	    DELT = TAUGRID_INPUT(N) - TAUGRID_INPUT(N-1)
	    EXT(N) = DELT/DEL
	  ENDDO

	  DO N = 1, NLAYERS
	    STH = GM_TOA * H(N)/H(0)
	    CTH = DSQRT ( ONE - STH * STH )
	    S0 = ZERO
	    HELP1 = H(0)*CTH
	    HELP2 = -H(0)*H(0)*STH*STH
	    DO M = 1, N
	      S1 = HELP1 - DSQRT(HELP2 + H(M)*H(M))
	      DELS = S1 - S0
	      TAUTHICK_INPUT(N,M) = EXT(M) * DELS
	      S0 = S1
	    ENDDO
	  ENDDO

C  Plane parallel

	ELSE

	  MU_TOA = DCOS ( SUN0 * DEG_TO_RAD )
	  DO N = 1, NLAYERS
	    DO M = 1, N
	      DELT = TAUGRID_INPUT(M) - TAUGRID_INPUT(M-1)
	      TAUTHICK_INPUT(N,M) = DELT / MU_TOA
	    ENDDO
	  ENDDO

	ENDIF

C  Finish

	END
