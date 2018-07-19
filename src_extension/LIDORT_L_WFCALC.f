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

	SUBROUTINE UPUSER_ATMOSWF
     I    ( DO_INCLUDE_ALBEDO,
     I      DO_INCLUDE_MVOUTPUT,
     I      DO_INCLUDE_DIRECTBEAM,
     I      R2, F1, FOURIER,
     I      LAYER_TO_VARY, K_PARAMETERS,
     I      USER_DIRECT_BEAM, DO_GMULT )

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include files of linearized setup and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Flux factor

	DOUBLE PRECISION F1

C  local control flags

	LOGICAL		 DO_INCLUDE_ALBEDO
	LOGICAL		 DO_INCLUDE_DIRECTBEAM
	LOGICAL		 DO_INCLUDE_MVOUTPUT

C  Linearization control

	INTEGER		 LAYER_TO_VARY, K_PARAMETERS

C  albedo factor

	DOUBLE PRECISION R2

C  Fourier number

	INTEGER		 FOURIER

C  reflected direct beam

	DOUBLE PRECISION USER_DIRECT_BEAM(MAX_USER_STREAMS)

C  Green function multiplier flag (Quadrature solutions only)

	LOGICAL		 DO_GMULT(MAX_OFFGRID_USERTAUS)

C  local variables
C  ---------------

	INTEGER		 N, NUT, NSTART, NUT_PREV, NLEVEL
	INTEGER		 UTA, UM, IUM, K, Q, NC, UT, AA, I, IQD
	DOUBLE PRECISION L_CUMUL_SOURCE(MAX_USER_STREAMS,MAX_PARAMETERS)
	DOUBLE PRECISION L_BOA_SOURCE(MAX_USER_STREAMS,MAX_PARAMETERS)
	DOUBLE PRECISION L_DIRBEAM(MAX_USER_STREAMS,MAX_PARAMETERS)
	DOUBLE PRECISION L_LAYER_SOURCE(MAX_USER_STREAMS,MAX_PARAMETERS)
	DOUBLE PRECISION
     &        L_MSCATUSE_LAYERSOURCE(MAX_USER_STREAMS,MAX_PARAMETERS)
	DOUBLE PRECISION VAR, TAU_DN, TAU_UP, H_UP, H_DN, L_FINAL_SOURCE
	LOGICAL		 DOCONS(MAXLAYER)

	K = LAYER_TO_VARY

C  Zero all Fourier components close to zenith

	IF ( DO_USER_STREAMS ) THEN
	  IF ( LOCAL_UM_START .GT. 1 ) THEN
	    DO UTA = 1, N_OUT_USERTAUS
	      DO Q = 1, K_PARAMETERS
	        DO UM = 1, LOCAL_UM_START - 1
	          IUM = USEROUTPUT_INDEX(UM)
	          ATMOSWF_F(Q,K,UTA,IUM,UPIDX) = ZERO
	        ENDDO
	      ENDDO
	    ENDDO
	  ENDIF
	ENDIF

C  Initialize recursion for user-defined stream angles only

	IF ( DO_USER_STREAMS ) THEN

	  CALL GET_L_BOASOURCE
     I    ( DO_INCLUDE_ALBEDO,
     I      DO_INCLUDE_DIRECTBEAM, R2, FOURIER, K, K_PARAMETERS,
     I      USER_DIRECT_BEAM, L_BOA_SOURCE, L_DIRBEAM )

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO Q = 1, K_PARAMETERS
	      L_CUMUL_SOURCE(UM,Q) = L_BOA_SOURCE(UM,Q)+L_DIRBEAM(UM,Q)
	    ENDDO
	  ENDDO

	  IF ( SAVE_LAYER_MSST ) THEN
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        ATMOSWF_MSCATBOA_STERM_F(Q,K,UM) = 
     &                     F1 * L_BOA_SOURCE(UM,Q)
	        ATMOSWF_DIRECTBOA_STERM_F(Q,K,UM) = 
     &                     F1 * L_DIRBEAM(UM,Q)
	      ENDDO
	    ENDDO
	  ENDIF

	ENDIF

C  initialise cumulative source term loop

 	NSTART = NLAYERS
	NUT_PREV = NSTART + 1
	DO N = 1, NLAYERS
	  DOCONS(N) = .TRUE.
	  L_HMULT_EXIST(N) = .FALSE.
	ENDDO

C  loop over all output optical depths
C  -----------------------------------

	DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

	  NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)
C    1. Get layer source terms
C    2. Find cumulative source term
C    3. Set multiple scatter source term (MSST) output if flagged

	  IF ( DO_USER_STREAMS ) THEN
	    NUT = NLEVEL + 1
	    DO N = NSTART, NUT, -1
	      NC = NLAYERS + 1 - N
	      CALL L_WHOLELAYER_STERM_UP
     &       ( FOURIER, N, K, K_PARAMETERS, DOCONS(N),
     O         L_LAYER_SOURCE, L_MSCATUSE_LAYERSOURCE )
	      IF ( N.EQ.K ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          DO Q = 1, K_PARAMETERS
	            L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)  +
     &                   T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q) +
     &                  L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_UP(UM,NC-1)
	          ENDDO
	        ENDDO
	      ELSE
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          DO Q = 1, K_PARAMETERS
	            L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)  +
     &                   T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q)
	          ENDDO
	        ENDDO
	      ENDIF
	      IF ( SAVE_LAYER_MSST ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          DO Q = 1, K_PARAMETERS
	            ATMOSWF_MSCATSTERM_F(Q,K,N,UM,UPIDX) =
     &                  L_MSCATUSE_LAYERSOURCE(UM,Q) * F1
	          ENDDO
	        ENDDO
	      ENDIF
	    ENDDO
	  ENDIF

C  Offgrid output
C  --------------

	  IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

	    UT = OFFGRID_UTAU_OUTINDEX(UTA)
	    N  = OFFGRID_UTAU_LAYERIDX(UT)

C  Linearized Transmittances (saved results).
  
	    IF ( DO_GMULT(UT) ) THEN
	      IF ( N .EQ. K ) THEN
	        TAU_DN = OFFGRID_UTAU_VALUES(UT)
	        TAU_UP = DELTA(N) - TAU_DN
	        DO AA = 1, NSTREAMS
	          H_DN = TAU_DN * T_UTDN_EIGEN(AA,UT)
	          H_UP = TAU_UP * T_UTUP_EIGEN(AA,UT)
	          DO Q = 1, K_PARAMETERS
	            VAR = - L_KEIGEN(AA,N,Q) - KEIGEN(AA,N) * VQ(N,Q)
	            L_T_UTDN_EIGEN(AA,UT,Q) = H_DN * VAR
	            L_T_UTUP_EIGEN(AA,UT,Q) = H_UP * VAR
	          ENDDO
	        ENDDO
	      ENDIF
	    ENDIF

C  Quadrature output at offgrid optical depths
C  ( Required if mean-value calculations are to be done)
C    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )
C  Copy into storage if quad output required

	    IF ( DO_QUAD_OUTPUT. OR. DO_INCLUDE_MVOUTPUT ) THEN
	      CALL QUADATMOSWF_OFFGRID_UP
     &         ( N, K, K_PARAMETERS, UTA, UT, F1, DO_GMULT(UT) )
	    ENDIF
	    IF ( DO_QUAD_OUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        IQD = QUADOUTPUT_INDEX(I)
	        DO Q = 1, K_PARAMETERS
 	          ATMOSWF_F(Q,K,UTA,IUM,UPIDX) = 
     &                    QUAD_ATMOSWF(Q,K,UTA,I,UPIDX)
	        ENDDO
	      ENDDO
	    ENDIF

C  User-defined stream output, add additional partial layer source term

	    IF ( DO_USER_STREAMS ) THEN
	      CALL L_PARTLAYER_STERM_UP
     I          ( FOURIER, N, UT, K, K_PARAMETERS, L_LAYER_SOURCE )
	      IF ( N.EQ.K ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          IUM = USEROUTPUT_INDEX(UM)
	          DO Q = 1, K_PARAMETERS
	            L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q)  +
     &                  T_UTUP_USERM(UT,UM)  * L_CUMUL_SOURCE(UM,Q) +
     &                L_T_UTUP_USERM(UT,UM,Q)*   CUMSOURCE_UP(UM,NC)
	            ATMOSWF_F(Q,K,UTA,IUM,UPIDX) = F1 * L_FINAL_SOURCE
	          ENDDO
	        ENDDO
	      ELSE
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          IUM = USEROUTPUT_INDEX(UM)
	          DO Q = 1, K_PARAMETERS
	            L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q)  +
     &                  T_UTUP_USERM(UT,UM)  * L_CUMUL_SOURCE(UM,Q)
	            ATMOSWF_F(Q,K,UTA,IUM,UPIDX) = F1 * L_FINAL_SOURCE
	          ENDDO
	        ENDDO
	      ENDIF
	    ENDIF

C  set multiplier flag

	    DO_GMULT(UT) = .FALSE.

C  Ongrid output
C  -------------

	  ELSE

C  Quadrature output at offgrid optical depths
C  ( Required if mean-value calculations are to be done)
C    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )
C  Copy into storage if quad output required

	    IF ( DO_QUAD_OUTPUT. OR. DO_INCLUDE_MVOUTPUT ) THEN
	      CALL QUADATMOSWF_LEVEL_UP
     *             ( NLEVEL, K, K_PARAMETERS, UTA, F1 )
	    ENDIF
	    IF ( DO_QUAD_OUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        IQD = QUADOUTPUT_INDEX(I)
	        DO Q = 1, K_PARAMETERS
 	          ATMOSWF_F(Q,K,UTA,IQD,UPIDX) = 
     &                    QUAD_ATMOSWF(Q,K,UTA,I,UPIDX)
	        ENDDO
	      ENDDO
	    ENDIF

C  User-defined stream output, just set to the cumulative source term

	    IF ( DO_USER_STREAMS ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        DO Q = 1, K_PARAMETERS
	          ATMOSWF_F(Q,K,UTA,IUM,UPIDX) =
     &                       F1 * L_CUMUL_SOURCE(UM,Q)
	        ENDDO
	      ENDDO
	    ENDIF

	  ENDIF

C  Check for updating the recursion

	  IF ( DO_USER_STREAMS ) THEN
	    IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
	    NUT_PREV = NUT
	  ENDIF

C  end loop over optical depth

	ENDDO

C  Finish

	END

C

	SUBROUTINE DNUSER_ATMOSWF
     I    ( FOURIER, F1, DO_GMULT, DO_INCLUDE_MVOUTPUT,
     I      LAYER_TO_VARY, K_PARAMETERS )

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include files of linearized setup and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments (flux factor constant)

	DOUBLE PRECISION F1
	INTEGER		 FOURIER

C  Linearization control

	INTEGER		 LAYER_TO_VARY, K_PARAMETERS

C  local MV output control

	LOGICAL		 DO_INCLUDE_MVOUTPUT

C  Green function multiplier flag (Quadrature solutions only)

	LOGICAL		 DO_GMULT(MAX_OFFGRID_USERTAUS)

C  local variables

	INTEGER		 N, NUT, NSTART, NUT_PREV, NLEVEL
	INTEGER		 UTA, UM, IUM, K, Q, NC, UT, AA, I, IQD
	DOUBLE PRECISION L_CUMUL_SOURCE(MAX_USER_STREAMS,MAX_PARAMETERS)
	DOUBLE PRECISION L_TOA_SOURCE(MAX_USER_STREAMS,MAX_PARAMETERS)
	DOUBLE PRECISION L_LAYER_SOURCE(MAX_USER_STREAMS,MAX_PARAMETERS)
	DOUBLE PRECISION
     &        L_MSCATUSE_LAYERSOURCE(MAX_USER_STREAMS,MAX_PARAMETERS)
	DOUBLE PRECISION VAR, TAU_DN, TAU_UP, H_UP, H_DN, L_FINAL_SOURCE
	LOGICAL		 DOCONS(MAXLAYER)

	K = LAYER_TO_VARY

C  Zero all Fourier components close to zenith

	IF ( DO_USER_STREAMS ) THEN
	  IF ( LOCAL_UM_START .GT. 1 ) THEN
	    DO UTA = 1, N_OUT_USERTAUS
	      DO Q = 1, K_PARAMETERS
	        DO UM = 1, LOCAL_UM_START - 1
	          IUM = USEROUTPUT_INDEX(UM)
	          ATMOSWF_F(Q,K,UTA,IUM,DNIDX) = ZERO
	        ENDDO
	      ENDDO
	    ENDDO
	  ENDIF
	ENDIF

C  Initialize recursion for user-defined stream angles only

	IF ( DO_USER_STREAMS ) THEN
	  CALL GET_L_TOASOURCE ( L_TOA_SOURCE, K_PARAMETERS )
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO Q = 1, K_PARAMETERS
	      L_CUMUL_SOURCE(UM,Q) = L_TOA_SOURCE(UM,Q)
	    ENDDO
	  ENDDO
	ENDIF

C  initialise cumulative source term loop

 	NSTART = 1
	NUT_PREV = NSTART - 1
	DO N = 1, NLAYERS
	  DOCONS(N) = .TRUE.
	  IF ( .NOT. DO_UPWELLING ) L_HMULT_EXIST(N) = .FALSE.
	ENDDO

C  loop over all output optical depths
C  -----------------------------------

	DO UTA = 1, N_OUT_USERTAUS

C  Layer index for given optical depth

	  NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)
C    1. Get layer source terms
C    2. Find cumulative source term
C    3. Set multiple scatter source term output if flagged

	  IF ( DO_USER_STREAMS ) THEN
	    NUT = NLEVEL
	    DO N = NSTART, NUT
	      NC = N
	      CALL L_WHOLELAYER_STERM_DN
     &       ( FOURIER, N, K, K_PARAMETERS, DOCONS(N),
     O         L_LAYER_SOURCE, L_MSCATUSE_LAYERSOURCE )
	      IF ( N.EQ.K ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          DO Q = 1, K_PARAMETERS
	            L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)  +
     &                   T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q) +
     &                  L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_DN(UM,NC-1)
	          ENDDO
	        ENDDO
	      ELSE
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          DO Q = 1, K_PARAMETERS
	            L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)  +
     &                   T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q)
	          ENDDO
	        ENDDO
	      ENDIF	
	      IF ( SAVE_LAYER_MSST ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          DO Q = 1, K_PARAMETERS
	            ATMOSWF_MSCATSTERM_F(Q,K,N,UM,DNIDX) =
     &                  L_MSCATUSE_LAYERSOURCE(UM,Q) * F1
	          ENDDO
	        ENDDO
	      ENDIF
	    ENDDO
	  ENDIF

C  Offgrid output
C  --------------

	  IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

	    UT = OFFGRID_UTAU_OUTINDEX(UTA)
	    N  = OFFGRID_UTAU_LAYERIDX(UT)

C  Linearized Transmittances (saved results).
  
	    IF ( DO_GMULT(UT) ) THEN
	      IF ( N .EQ. K ) THEN
	        TAU_DN = OFFGRID_UTAU_VALUES(UT)
	        TAU_UP = DELTA(N) - TAU_DN
	        DO AA = 1, NSTREAMS
	          H_DN = TAU_DN * T_UTDN_EIGEN(AA,UT)
	          H_UP = TAU_UP * T_UTUP_EIGEN(AA,UT)
	          DO Q = 1, K_PARAMETERS
	            VAR = - L_KEIGEN(AA,N,Q) - KEIGEN(AA,N) * VQ(N,Q)
	            L_T_UTDN_EIGEN(AA,UT,Q) = H_DN * VAR
	            L_T_UTUP_EIGEN(AA,UT,Q) = H_UP * VAR
	          ENDDO
	        ENDDO
	      ENDIF
	    ENDIF

C  Quadrature output at offgrid optical depths
C  ( Required if mean-value calculations are to be done)
C    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )
C  Copy into storage if quad output required

	    IF ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT ) THEN
	      CALL QUADATMOSWF_OFFGRID_DN
     &         ( N, K, K_PARAMETERS, UTA, UT, F1, DO_GMULT(UT) )
	    ENDIF
	    IF ( DO_QUAD_OUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        IQD = QUADOUTPUT_INDEX(I)
	        DO Q = 1, K_PARAMETERS
 	          ATMOSWF_F(Q,K,UTA,IUM,DNIDX) = 
     &                    QUAD_ATMOSWF(Q,K,UTA,I,DNIDX)
	        ENDDO
	      ENDDO
	    ENDIF

C  User-defined stream output, add additional partial layer source term

	    IF ( DO_USER_STREAMS ) THEN
	      CALL L_PARTLAYER_STERM_DN
     I          ( FOURIER, N, UT, K, K_PARAMETERS, L_LAYER_SOURCE )
	      IF ( N.EQ.K ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          IUM = USEROUTPUT_INDEX(UM)
	          DO Q = 1, K_PARAMETERS
	            L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q)  +
     &                  T_UTDN_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,Q) +
     &                L_T_UTDN_USERM(UT,UM,Q) *   CUMSOURCE_DN(UM,NC)
	          ATMOSWF_F(Q,K,UTA,IUM,DNIDX) = F1 * L_FINAL_SOURCE
	          ENDDO
	        ENDDO
	      ELSE
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          IUM = USEROUTPUT_INDEX(UM)
	          DO Q = 1, K_PARAMETERS
	            L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q)  +
     &                  T_UTDN_USERM(UT,UM)  * L_CUMUL_SOURCE(UM,Q)
	          ATMOSWF_F(Q,K,UTA,IUM,DNIDX) = F1 * L_FINAL_SOURCE
	          ENDDO
	        ENDDO
	      ENDIF
	    ENDIF

C  set multiplier flag

	    DO_GMULT(UT) = .FALSE.

C  Ongrid output
C  -------------

	  ELSE

C  Quadrature output at offgrid optical depths
C  ( Required if mean-value calculations are to be done)
C    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )
C  Copy into storage if quad output required

	    IF ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT ) THEN
	      CALL QUADATMOSWF_LEVEL_DN
     *             ( NLEVEL, K, K_PARAMETERS, UTA, F1 )
	    ENDIF
	    IF ( DO_QUAD_OUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        IQD = QUADOUTPUT_INDEX(I)
	        DO Q = 1, K_PARAMETERS
 	          ATMOSWF_F(Q,K,UTA,IQD,DNIDX) = 
     &                    QUAD_ATMOSWF(Q,K,UTA,I,DNIDX)
	        ENDDO
	      ENDDO
	    ENDIF

C  User-defined stream output, just set to the cumulative source term

	    IF ( DO_USER_STREAMS ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        DO Q = 1, K_PARAMETERS
	          ATMOSWF_F(Q,K,UTA,IUM,DNIDX) =
     &                       F1 * L_CUMUL_SOURCE(UM,Q)
	        ENDDO
	      ENDDO
	    ENDIF

	  ENDIF

C  Check for updating the recursion

	  IF ( DO_USER_STREAMS ) THEN
	    IF ( NUT. NE. NUT_PREV ) NSTART = NUT + 1
	    NUT_PREV = NUT
	  ENDIF

C  end loop over optical depth

	ENDDO

C  Finish

	END

C

	SUBROUTINE QUADATMOSWF_LEVEL_UP
     *             ( NL, K, K_PARAMETERS, UTA, F1 )

C  Upwelling weighting function Fourier components at level boundary NL
C  Quadrature angles only

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model variables (input)

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include files of linearized setup and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments (Flux factor, level)

	INTEGER		 NL, UTA
	DOUBLE PRECISION F1

C  linearization control

	INTEGER		 K, K_PARAMETERS

C  local variables

	INTEGER		 N, I, I1, AA, Q
	DOUBLE PRECISION SPAR, SHOM, HOM1, HOM2, HOM3, HOM4, HOM5

C  homogeneous and particular solution contributions SHOM and SPAR

C  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
C  looking at the perturbation field at the top of these layers. The
C  case where the level mask = NLAYERS is the upwelling perturbed fields
C  at the bottom of the atmosphere (treated separately).

	N = NL + 1

C  For the lowest level

	IF ( NL .EQ. NLAYERS ) THEN

C  If this is also the layer that is varying, extra contributions

	  IF ( K .EQ. NL ) THEN

	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          HOM1 = NCON_XVEC(I1,AA,NL,Q) * T_DELT_EIGEN(AA,NL)
	          HOM2 = LCON_XVEC(I1,AA,NL) * L_T_DELT_EIGEN(AA,NL,Q)
	          HOM3 =
     &              LCON(AA,NL)*T_DELT_EIGEN(AA,NL)*L_XPOS(I1,AA,NL,Q)
	          HOM4 = PCON_XVEC(I1,AA,NL,Q)
	          HOM5 = MCON(AA,NL) * L_XNEG(I1,AA,NL,Q)
	          SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
	        ENDDO
                SPAR = L_WLOWER(I1,NL,Q)
 	        QUAD_ATMOSWF(Q,K,UTA,I,UPIDX) = F1 * ( SPAR + SHOM )
	      ENDDO
	    ENDDO

C  non-varying lowest layer

	  ELSE

	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          HOM1 = NCON_XVEC(I1,AA,NL,Q) * T_DELT_EIGEN(AA,NL)
	          HOM2 = PCON_XVEC(I1,AA,NL,Q)
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
                SPAR = L_WLOWER(I1,NL,Q)
 	        QUAD_ATMOSWF(Q,K,UTA,I,UPIDX) = F1 * ( SPAR + SHOM )
	      ENDDO
	    ENDDO

	  ENDIF

C  For other levels in the atmosphere
C  ----------------------------------

	ELSE

C  If this is also the layer that is varying, extra contributions

	  IF ( K .EQ. N ) THEN

	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          HOM1 = NCON_XVEC(I1,AA,N,Q) 
	          HOM2 = LCON(AA,N) * L_XPOS(I1,AA,N,Q)
	          HOM3 = MCON_XVEC(I1,AA,N) * L_T_DELT_EIGEN(AA,N,Q)
	          HOM4 = MCON(AA,N)*T_DELT_EIGEN(AA,N)*L_XNEG(I1,AA,N,Q)
	          HOM5 = PCON_XVEC(I1,AA,N,Q) * T_DELT_EIGEN(AA,N)
	          SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
	        ENDDO
                SPAR = L_WUPPER(I1,N,Q)
 	        QUAD_ATMOSWF(Q,K,UTA,I,UPIDX) = F1 * ( SPAR + SHOM )
	      ENDDO
	    ENDDO

C  non-varying layer lower than the varying layer

	  ELSE IF ( K.LT.N ) THEN

	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          HOM1 = NCON_XVEC(I1,AA,N,Q)
	          HOM2 = PCON_XVEC(I1,AA,N,Q) * T_DELT_EIGEN(AA,N)
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
                SPAR = L_WUPPER(I1,N,Q)
 	        QUAD_ATMOSWF(Q,K,UTA,I,UPIDX) = F1 * ( SPAR + SHOM )
	      ENDDO
	    ENDDO

C  non-varying layer higher than the varying layer

	  ELSE IF ( K.GT.N ) THEN

	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          HOM1 = NCON_XVEC(I1,AA,N,Q)
	          HOM2 = PCON_XVEC(I1,AA,N,Q) * T_DELT_EIGEN(AA,N)
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
	        QUAD_ATMOSWF(Q,K,UTA,I,UPIDX) = F1 * SHOM
	      ENDDO
	    ENDDO

	  ENDIF

	ENDIF

C  Finish

	END

C

	SUBROUTINE QUADATMOSWF_LEVEL_DN
     *             ( NLEVEL, K, K_PARAMETERS, UTA, F1 )

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model variables (input)

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include files of linearized setup and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments (Flux factor, level)

	INTEGER		 NLEVEL, UTA
	DOUBLE PRECISION F1

C  linearization control

	INTEGER		 K, K_PARAMETERS

C  local variables

	INTEGER		 N, I, AA, Q
	DOUBLE PRECISION SPAR, SHOM, HOM1, HOM2, HOM3, HOM4, HOM5

C  homogeneous and particular solution contributions SHOM and SPAR

	N = NLEVEL

C  Downwelling weighting function at TOA ( or N = 0 ) is zero

	IF ( NLEVEL .EQ. 0 ) THEN

	  DO I = 1, NSTREAMS
	    DO Q = 1, K_PARAMETERS
 	      QUAD_ATMOSWF(Q,K,UTA,I,DNIDX) = ZERO
	    ENDDO
	  ENDDO

C  For other levels in the atmosphere
C  ----------------------------------

	ELSE

C  If this is also the layer that is varying, extra contributions

	  IF ( K .EQ. N ) THEN

	    DO I = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
	          HOM2 = LCON_XVEC(I,AA,N) * L_T_DELT_EIGEN(AA,N,Q)
	          HOM3 = LCON(AA,N)*T_DELT_EIGEN(AA,N)*L_XPOS(I,AA,N,Q)
	          HOM4 = PCON_XVEC(I,AA,N,Q)
	          HOM5 = MCON(AA,N) * L_XNEG(I,AA,N,Q)
	          SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
	        ENDDO
                SPAR = L_WLOWER(I,N,Q)
 	        QUAD_ATMOSWF(Q,K,UTA,I,DNIDX) = F1 * ( SPAR + SHOM )
	      ENDDO
	    ENDDO

C  non-varying layer lower than the varying layer

	  ELSE IF (K.LT.N) THEN

	    DO I = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
	          HOM2 = PCON_XVEC(I,AA,N,Q) 
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
                SPAR = L_WLOWER(I,N,Q)
 	        QUAD_ATMOSWF(Q,K,UTA,I,DNIDX) = F1 * ( SPAR + SHOM )
	      ENDDO
	    ENDDO

C  non-varying layer higher than the varying layer

	  ELSE IF (K.GT.N) THEN

	    DO I = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
	          HOM2 = PCON_XVEC(I,AA,N,Q) 
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
 	        QUAD_ATMOSWF(Q,K,UTA,I,DNIDX) = F1 * SHOM
	      ENDDO
	    ENDDO

	  ENDIF

	ENDIF

C  Finish

	END

C

	SUBROUTINE QUADATMOSWF_OFFGRID_UP
     &            ( N, K, K_PARAMETERS, UTA, UT, F1, DO_GMULT )

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments (Flux factor, layer and offgrid index)

	INTEGER		 N, K, K_PARAMETERS, UTA, UT
	DOUBLE PRECISION F1
	LOGICAL		 DO_GMULT

C  local variables

	INTEGER		 I, I1, AA, Q
	DOUBLE PRECISION SPAR, SHOM, HOM1, HOM2, PAR1, PAR2
	DOUBLE PRECISION H1, H2, H3, H4, H5, H6

C  For those optical depths at off-grid levels
C  -------------------------------------------

C  Classical solution

	IF ( DO_CLASSICAL_SOLUTION ) THEN

C  solution for N being the varying layer K

	  IF ( N. EQ. K ) THEN

	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      DO Q = 1, K_PARAMETERS
                SPAR = L_WUPPER(I1,N,Q) * T_UTDN_MUBAR(UT) +
     &                   WUPPER(I1,N)   * L_T_UTDN_MUBAR(UT,N,Q)
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          H1 = LCON_XVEC(I1,AA,N) * L_T_UTDN_EIGEN(AA,UT,Q)
	          H2 = NCON(AA,N,Q) * XPOS(I1,AA,N)
                  H3 = LCON(AA,N)   * L_XPOS(I1,AA,N,Q)
	          HOM1 = H1 + T_UTDN_EIGEN(AA,UT) * ( H2 + H3 )
	          H4 = MCON_XVEC(I1,AA,N) * L_T_UTUP_EIGEN(AA,UT,Q)
	          H5 = PCON(AA,N,Q) * XNEG(I1,AA,N)
                  H6 = MCON(AA,N)   * L_XNEG(I1,AA,N,Q)
	          HOM2 = H4 + T_UTUP_EIGEN(AA,UT) * ( H5 + H6 )
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
	        QUAD_ATMOSWF(Q,K,UTA,I,UPIDX) = F1 * ( SPAR + SHOM )
	      ENDDO
	    ENDDO

C  Solution for N > K (N below the layer that varies)

	  ELSE IF ( N. GT. K ) THEN

	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      DO Q = 1, K_PARAMETERS
                SPAR = L_WUPPER(I1,N,Q) *   T_UTDN_MUBAR(UT) +
     &                   WUPPER(I1,N)   * L_T_UTDN_MUBAR(UT,K,Q)
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          H2 = NCON(AA,N,Q) * XPOS(I1,AA,N)
	          HOM1 = T_UTDN_EIGEN(AA,UT) * H2
	          H5 = PCON(AA,N,Q) * XNEG(I1,AA,N)
	          HOM2 = T_UTUP_EIGEN(AA,UT) * H5
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
	        QUAD_ATMOSWF(Q,K,UTA,I,UPIDX) = F1 * ( SPAR + SHOM )
	      ENDDO
	    ENDDO

C  Solution for N < K (N above the layer that varies)

	  ELSE IF ( N. LT. K ) THEN

	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          H2 = NCON(AA,N,Q) * XPOS(I1,AA,N)
	          HOM1 = T_UTDN_EIGEN(AA,UT) * H2
	          H5 = PCON(AA,N,Q) * XNEG(I1,AA,N)
	          HOM2 = T_UTUP_EIGEN(AA,UT) * H5
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
	        QUAD_ATMOSWF(Q,K,UTA,I,UPIDX) = F1 * SHOM
	      ENDDO
	    ENDDO

	  ENDIF

C  Green function solution
C  =======================

	ELSE

C  get the multipliers

	  CALL L_QUAD_GFUNCMULT
     I      ( N, UT, K, K_PARAMETERS, DO_GMULT )

C  solution for N being the varying layer K

	  IF ( N. EQ. K ) THEN

	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = ZERO
	        DO AA = 1, NSTREAMS
	          PAR1 = L_XPOS(I,AA,N,Q)  * UT_GMULT_UP(AA,UT) +
     &                     XPOS(I,AA,N)    * L_UT_GMULT_UP(AA,UT,K,Q)
	          PAR2 = L_XPOS(I1,AA,N,Q)  * UT_GMULT_DN(AA,UT) +
     &                     XPOS(I1,AA,N)    * L_UT_GMULT_DN(AA,UT,K,Q)
	          SPAR = SPAR + PAR1 + PAR2
	        ENDDO
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          H1 = LCON_XVEC(I1,AA,N) * L_T_UTDN_EIGEN(AA,UT,Q)
	          H2 = NCON(AA,N,Q) * XPOS(I1,AA,N)
                  H3 = LCON(AA,N)   * L_XPOS(I1,AA,N,Q)
	          HOM1 = H1 + T_UTDN_EIGEN(AA,UT) * ( H2 + H3 )
	          H4 = MCON_XVEC(I1,AA,N) * L_T_UTUP_EIGEN(AA,UT,Q)
	          H5 = PCON(AA,N,Q) * XNEG(I1,AA,N)
                  H6 = MCON(AA,N)   * L_XNEG(I1,AA,N,Q)
	          HOM2 = H4 + T_UTUP_EIGEN(AA,UT) * ( H5 + H6 )
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
	        QUAD_ATMOSWF(Q,K,UTA,I,UPIDX) = F1 * ( SPAR + SHOM )
	      ENDDO
	    ENDDO

C  Solution for N > K (N below the layer that varies)

	  ELSE IF ( N. GT. K ) THEN

	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = ZERO
	        DO AA = 1, NSTREAMS
	          PAR1 = XPOS(I,AA,N)  * L_UT_GMULT_UP(AA,UT,K,Q)
	          PAR2 = XPOS(I1,AA,N) * L_UT_GMULT_DN(AA,UT,K,Q)
	          SPAR = SPAR + PAR1 + PAR2
	        ENDDO
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          H2 = NCON(AA,N,Q) * XPOS(I1,AA,N)
	          HOM1 = T_UTDN_EIGEN(AA,UT) * H2
	          H5 = PCON(AA,N,Q) * XNEG(I1,AA,N)
	          HOM2 = T_UTUP_EIGEN(AA,UT) * H5
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
	        QUAD_ATMOSWF(Q,K,UTA,I,UPIDX) = F1 * ( SPAR + SHOM )
	      ENDDO
	    ENDDO

C  Solution for N < K (N above the layer that varies)

	  ELSE IF ( N. LT. K ) THEN

	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          H2 = NCON(AA,N,Q) * XPOS(I1,AA,N)
	          HOM1 = T_UTDN_EIGEN(AA,UT) * H2
	          H5 = PCON(AA,N,Q) * XNEG(I1,AA,N)
	          HOM2 = T_UTUP_EIGEN(AA,UT) * H5
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
	        QUAD_ATMOSWF(Q,K,UTA,I,UPIDX) = F1 * SHOM
	      ENDDO
	    ENDDO

	  ENDIF

	ENDIF

C  Finish

	END

C

	SUBROUTINE QUADATMOSWF_OFFGRID_DN
     &            ( N, K, K_PARAMETERS, UTA, UT, F1, DO_GMULT )

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments (Flux factor, layer and offgrid index)

	INTEGER		 N, K, K_PARAMETERS, UTA, UT
	DOUBLE PRECISION F1
	LOGICAL		 DO_GMULT

C  local variables

	INTEGER		 I, I1, AA, Q
	DOUBLE PRECISION SPAR, SHOM, HOM1, HOM2, PAR1, PAR2
	DOUBLE PRECISION H1, H2, H3, H4, H5, H6

C  For those optical depths at off-grid levels
C  -------------------------------------------

C  Classical solution

	IF ( DO_CLASSICAL_SOLUTION ) THEN

C  solution for N being the varying layer K

	  IF ( N. EQ. K ) THEN

	    DO I = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
                SPAR = L_WUPPER(I,N,Q) *   T_UTDN_MUBAR(UT) +
     &                   WUPPER(I,N)   * L_T_UTDN_MUBAR(UT,N,Q)
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          H1 = LCON_XVEC(I,AA,N) * L_T_UTDN_EIGEN(AA,UT,Q)
	          H2 = NCON(AA,N,Q) * XPOS(I,AA,N)
                  H3 = LCON(AA,N)   * L_XPOS(I,AA,N,Q)
	          HOM1 = H1 + T_UTDN_EIGEN(AA,UT) * ( H2 + H3 )
	          H4 = MCON_XVEC(I,AA,N) * L_T_UTUP_EIGEN(AA,UT,Q)
	          H5 = PCON(AA,N,Q) * XNEG(I,AA,N)
                  H6 = MCON(AA,N)   * L_XNEG(I,AA,N,Q)
	          HOM2 = H4 + T_UTUP_EIGEN(AA,UT) * ( H5 + H6 )
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
	        QUAD_ATMOSWF(Q,K,UTA,I,DNIDX) = F1 * ( SPAR + SHOM )
	      ENDDO
	    ENDDO

C  Solution for N > K (N above the layer that varies)

	  ELSE IF ( N. GT. K ) THEN

	    DO I = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
                SPAR = L_WUPPER(I,N,Q) *   T_UTDN_MUBAR(UT) +
     &                   WUPPER(I,N)   * L_T_UTDN_MUBAR(UT,K,Q)
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          H2 = NCON(AA,N,Q) * XPOS(I,AA,N)
	          HOM1 = T_UTDN_EIGEN(AA,UT) * H2
	          H5 = PCON(AA,N,Q) * XNEG(I,AA,N)
	          HOM2 = T_UTUP_EIGEN(AA,UT) * H5
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
	        QUAD_ATMOSWF(Q,K,UTA,I,DNIDX) = F1 * ( SPAR + SHOM )
	      ENDDO
	    ENDDO

C  Solution for N < K (N below the layer that varies)

	  ELSE IF ( N. LT. K ) THEN

	    DO I = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          H2 = NCON(AA,N,Q) * XPOS(I,AA,N)
	          HOM1 = T_UTDN_EIGEN(AA,UT) * H2
	          H5 = PCON(AA,N,Q) * XNEG(I,AA,N)
	          HOM2 = T_UTUP_EIGEN(AA,UT) * H5
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
	        QUAD_ATMOSWF(Q,K,UTA,I,DNIDX) = F1 * SHOM
	      ENDDO
	    ENDDO

	  ENDIF

C  Green function solution
C  =======================

	ELSE

C  Get the multipliers if flagged

	  IF ( DO_GMULT ) THEN
	    CALL L_QUAD_GFUNCMULT
     I      ( N, UT, K, K_PARAMETERS, DO_GMULT )
	  ENDIF

C  solution for N being the varying layer K

	  IF ( N. EQ. K ) THEN

	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = ZERO
	        DO AA = 1, NSTREAMS
	          PAR1 = L_XPOS(I1,AA,N,Q)  * UT_GMULT_UP(AA,UT) +
     &                     XPOS(I1,AA,N)    * L_UT_GMULT_UP(AA,UT,K,Q)
	          PAR2 = L_XPOS(I,AA,N,Q)   * UT_GMULT_DN(AA,UT) +
     &                     XPOS(I,AA,N)     * L_UT_GMULT_DN(AA,UT,K,Q)
	          SPAR = SPAR + PAR1 + PAR2
	        ENDDO
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          H1 = LCON_XVEC(I,AA,N) * L_T_UTDN_EIGEN(AA,UT,Q)
	          H2 = NCON(AA,N,Q) * XPOS(I,AA,N)
                  H3 = LCON(AA,N)   * L_XPOS(I,AA,N,Q)
	          HOM1 = H1 + T_UTDN_EIGEN(AA,UT) * ( H2 + H3 )
	          H4 = MCON_XVEC(I,AA,N) * L_T_UTUP_EIGEN(AA,UT,Q)
	          H5 = PCON(AA,N,Q) * XNEG(I,AA,N)
                  H6 = MCON(AA,N)   * L_XNEG(I,AA,N,Q)
	          HOM2 = H4 + T_UTUP_EIGEN(AA,UT) * ( H5 + H6 )
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
	        QUAD_ATMOSWF(Q,K,UTA,I,DNIDX) = F1 * ( SPAR + SHOM )
	      ENDDO
	    ENDDO

C  Solution for N > K (N above the layer that varies)

	  ELSE IF ( N. GT. K ) THEN

	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = ZERO
	        DO AA = 1, NSTREAMS
	          PAR1 = XPOS(I1,AA,N) * L_UT_GMULT_UP(AA,UT,K,Q)
	          PAR2 = XPOS(I,AA,N)  * L_UT_GMULT_DN(AA,UT,K,Q)
	          SPAR = SPAR + PAR1 + PAR2
	        ENDDO
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          H2 = NCON(AA,N,Q) * XPOS(I,AA,N)
	          HOM1 = T_UTDN_EIGEN(AA,UT) * H2
	          H5 = PCON(AA,N,Q) * XNEG(I,AA,N)
	          HOM2 = T_UTUP_EIGEN(AA,UT) * H5
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
	        QUAD_ATMOSWF(Q,K,UTA,I,DNIDX) = F1 * ( SPAR + SHOM )
	      ENDDO
	    ENDDO

C  Solution for N < K (N below the layer that varies)

	  ELSE IF ( N. LT. K ) THEN

	    DO I = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          H2 = NCON(AA,N,Q) * XPOS(I,AA,N)
	          HOM1 = T_UTDN_EIGEN(AA,UT) * H2
	          H5 = PCON(AA,N,Q) * XNEG(I,AA,N)
	          HOM2 = T_UTUP_EIGEN(AA,UT) * H5
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
	        QUAD_ATMOSWF(Q,K,UTA,I,DNIDX) = F1 * SHOM
	      ENDDO
	    ENDDO

	  ENDIF

	ENDIF

C  Finish

	END

C

	SUBROUTINE LIDORT_MIFLUX_ATMOSWF ( K, K_PARAMETERS )

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include files of model and setup inputs

	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

C  Linearization control

	INTEGER		 K, K_PARAMETERS

C  local variables

	INTEGER		 I, IDIR, WDIR, UTA, UT, Q, N
	DOUBLE PRECISION SUM_MI, SUM_FX, FMU0
	DOUBLE PRECISION L_TRANS, L_DIRECT_FLUX, L_DIRECT_MEANI

C  mean intensity and flux
C  -----------------------

C  direction loop

	DO IDIR = 1, N_DIRECTIONS

	  WDIR = WHICH_DIRECTIONS(IDIR)

C  loop over all user-defined optical depths

	  DO UTA = 1, N_OUT_USERTAUS

C  loop over all parameters in the varying layer

	    DO Q = 1, K_PARAMETERS

C  integrations

	      SUM_MI = ZERO
	      SUM_FX = ZERO
	      DO I = 1, NSTREAMS
	        SUM_MI = SUM_MI + A(I)  * QUAD_ATMOSWF(Q,K,UTA,I,WDIR)
	        SUM_FX = SUM_FX + AX(I) * QUAD_ATMOSWF(Q,K,UTA,I,WDIR)
	      ENDDO
	      MINT_ATMOSWF(Q,K,UTA,WDIR) = SUM_MI * HALF
	      FLUX_ATMOSWF(Q,K,UTA,WDIR) = SUM_FX * PI2

C  end loops

	    ENDDO
	  ENDDO

C  For the downward direction, add the direct beam contributions

	  IF ( WDIR .EQ. DNIDX ) THEN

C  loop over all the output optical depths

	    DO UTA = 1, N_OUT_USERTAUS

C  For the offgrid values

	      IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
	        UT = OFFGRID_UTAU_OUTINDEX(UTA)
	        N  = OFFGRID_UTAU_LAYERIDX(UT)

C  Only contributions for layers above the PI cutoff

	        IF ( N .LE. LAYER_PIS_CUTOFF ) THEN
	          FMU0 = LOCAL_SZA(N) * FLUX_FACTOR
	          IF ( K.LE.N ) THEN
	            DO Q = 1, K_PARAMETERS
	              L_TRANS = L_INITIAL_TRANS(N,K,Q) * T_UTDN_MUBAR(UT) +
     &                      INITIAL_TRANS(N) * L_T_UTDN_MUBAR(UT,K,Q)
	              L_DIRECT_MEANI = FLUX_FACTOR * L_TRANS / PI4
		      MINT_ATMOSWF(Q,K,UTA,WDIR) =
     *                  MINT_ATMOSWF(Q,K,UTA,WDIR) + L_DIRECT_MEANI
	              L_DIRECT_FLUX  = FMU0 * L_TRANS
		      FLUX_ATMOSWF(Q,K,UTA,WDIR) =
     *                  FLUX_ATMOSWF(Q,K,UTA,WDIR) + L_DIRECT_FLUX
	            ENDDO
	          ENDIF
	        ENDIF

C  For the on-grid balues

	      ELSE
	        N = UTAU_LEVEL_MASK_DN(UTA)
	        IF ( N .LE. LAYER_PIS_CUTOFF ) THEN
	          IF ( N.GT.0 ) THEN
	            FMU0 = LOCAL_SZA(N) * FLUX_FACTOR
	            IF ( K.LE.N ) THEN
	              DO Q = 1, K_PARAMETERS
	                L_TRANS = L_INITIAL_TRANS(N,K,Q) * T_DELT_MUBAR(N) +
     &                         INITIAL_TRANS(N) * L_T_DELT_MUBAR(N,K,Q)
	              L_DIRECT_MEANI = FLUX_FACTOR * L_TRANS / PI4
		      MINT_ATMOSWF(Q,K,UTA,WDIR) =
     *                  MINT_ATMOSWF(Q,K,UTA,WDIR) + L_DIRECT_MEANI
	              L_DIRECT_FLUX  = FMU0 * L_TRANS
		      FLUX_ATMOSWF(Q,K,UTA,WDIR) =
     *                  FLUX_ATMOSWF(Q,K,UTA,WDIR) + L_DIRECT_FLUX
	              ENDDO
	            ENDIF
	          ENDIF
	        ENDIF
	      ENDIF

C  Close loops

      	    ENDDO
	  ENDIF   

C  end direction loop

	ENDDO

C  Finish

	END

C

  	SUBROUTINE GET_L_TOASOURCE ( L_TOA_SOURCE, K_PARAMETERS )

C  Include files
C  -------------

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model variables (input)

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  Subroutine arguments

	INTEGER		 K_PARAMETERS
	DOUBLE PRECISION L_TOA_SOURCE(MAX_USER_STREAMS,MAX_PARAMETERS)

C  local variables

	INTEGER		 UM, Q

C  initialise TOA source function

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  DO Q = 1, K_PARAMETERS
	    L_TOA_SOURCE(UM,Q) = ZERO
	  ENDDO
	ENDDO

C  Finish

	END

C

  	SUBROUTINE GET_L_BOASOURCE
     I    ( DO_INCLUDE_ALBEDO,
     I      DO_INCLUDE_DIRECTBEAM, R2, FOURIER, K, K_PARAMETERS,
     I      USER_DIRECT_BEAM, L_BOA_SOURCE, L_DIRBEAM)

C  Bottom of the atmosphere source term

C  Include files
C  -------------

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  Include file of Geophysical variables (input surface reflectivities)

	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'

C  include files of linearized setup and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  Subroutine input arguments
C  --------------------------

C  local control flags

	LOGICAL		 DO_INCLUDE_ALBEDO
	LOGICAL		 DO_INCLUDE_DIRECTBEAM

C  albedo and Fourier index

	DOUBLE PRECISION R2
	INTEGER		 FOURIER

C  reflected direct beam

	DOUBLE PRECISION USER_DIRECT_BEAM(MAX_USER_STREAMS)

C  linearization control

	INTEGER		 K, K_PARAMETERS

C  output

	DOUBLE PRECISION L_BOA_SOURCE(MAX_USER_STREAMS,MAX_PARAMETERS)
	DOUBLE PRECISION L_DIRBEAM(MAX_USER_STREAMS,MAX_PARAMETERS)

C  local variables
C  ---------------

	INTEGER		 M, N, J, I, UM, AA, Q
	DOUBLE PRECISION L_DOWN(MAXSTRM,MAX_PARAMETERS)
	DOUBLE PRECISION HELP(MAXSTRM,MAX_PARAMETERS)
	DOUBLE PRECISION REFLEC, L_BEAM, FAC
	DOUBLE PRECISION SPAR, SHOM, HOM1, HOM2, HOM3, HOM4, HOM5

C  Fourier number

	M = FOURIER

C  initialise boa source function

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  DO Q = 1, K_PARAMETERS
	    L_BOA_SOURCE(UM,Q) = ZERO
	    L_DIRBEAM(UM,Q) = ZERO
	  ENDDO
	ENDDO

	N = NLAYERS

C  reflectance from surface
C  ------------------------

	IF ( DO_INCLUDE_ALBEDO ) THEN

C  Downward intensity at computational angles (beam and homog)

C  If this is also the layer that is varying, extra contributions

	  IF ( K .EQ. N ) THEN

	    DO I = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
	          HOM2 = LCON_XVEC(I,AA,N) * L_T_DELT_EIGEN(AA,N,Q)
	          HOM3 = LCON(AA,N)*T_DELT_EIGEN(AA,N)*L_XPOS(I,AA,N,Q)
	          HOM4 = PCON_XVEC(I,AA,N,Q)
	          HOM5 = MCON(AA,N) * L_XNEG(I,AA,N,Q)
	          SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
	        ENDDO
                SPAR = L_WLOWER(I,N,Q)
 	        L_DOWN(I,Q) = SPAR + SHOM
	      ENDDO
	    ENDDO

C  non-varying layer lower than the varying layer

	  ELSE IF (K.LT.N) THEN

	    DO I = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
	          HOM2 = PCON_XVEC(I,AA,N,Q) 
	          SHOM = SHOM + HOM1 + HOM2
	        ENDDO
                SPAR = L_WLOWER(I,N,Q)
 	        L_DOWN(I,Q) = SPAR + SHOM
	      ENDDO
	    ENDDO

	  ENDIF

C  reflectance integrand  a(j).x(j).L_DOWN(-j)

	  DO Q = 1, K_PARAMETERS
	    DO I = 1, NSTREAMS
	      HELP(I,Q) = L_DOWN(I,Q) * AX(I)
	    ENDDO
	  ENDDO

C  reflected multiple scatter intensity at user defined-angles
C  -----------------------------------------------------------

C  .. integrate reflectance, same for all user-streams in Lambertian case

	  IF ( DO_LAMBERTIAN_ALBEDO ) THEN

	    DO Q = 1, K_PARAMETERS
	      REFLEC = ZERO
	      DO J = 1, NSTREAMS
	        REFLEC = REFLEC + HELP(J,Q)
	      ENDDO
	      REFLEC = R2 * REFLEC
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        L_BOA_SOURCE(UM,Q) = REFLEC
	      ENDDO
	    ENDDO

C  .. integrate with reflectance function at user angles (non-Lambertian)

	  ELSE

	    DO Q = 1, K_PARAMETERS
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        REFLEC = ZERO
	        DO J = 1, NSTREAMS
	          REFLEC = REFLEC + HELP(J,Q) * USER_BIREFLEC(M,UM,J)
	        ENDDO
	        L_BOA_SOURCE(UM,Q) = R2 * REFLEC
	      ENDDO
	    ENDDO

	  ENDIF
	  
C  Add direct beam if flagged

	  IF ( DO_INCLUDE_DIRECTBEAM ) THEN
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      FAC = - USER_DIRECT_BEAM(UM) * TAUTHICK(N,K) 
	      DO Q = 1, K_PARAMETERS
	        L_BEAM = VQ(K,Q) * FAC
	        L_DIRBEAM(UM,Q) = L_BEAM
	      ENDDO
	    ENDDO
	  ENDIF

	ENDIF

C  Finish

	END

C

	SUBROUTINE L_WHOLELAYER_STERM_UP
     I       ( FOURIER, GIVEN_LAYER,
     I         LAYER_TO_VARY, K_PARAMETERS,
     I         DO_CONSTANTS_SET,
     O         L_LAYERSOURCE, L_MSCATUSE_LAYERSOURCE )

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of multiplier and solution variables (input)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized multiplier and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  layer and Fourier number

	INTEGER		 GIVEN_LAYER, FOURIER

C  constants set flag

	LOGICAL		 DO_CONSTANTS_SET

C  output linearized layer source terms

	DOUBLE PRECISION
     &        L_LAYERSOURCE(MAX_USER_STREAMS,MAX_PARAMETERS),
     &        L_MSCATUSE_LAYERSOURCE(MAX_USER_STREAMS,MAX_PARAMETERS)

C  Linearization control

	INTEGER		 LAYER_TO_VARY, K_PARAMETERS

C  local variables

	INTEGER		 N, K, AA, UM, Q
	DOUBLE PRECISION SHOM, SFOR, SPAR, H1, H2, H3, H4, H5, H6

C  layer variables

	N = GIVEN_LAYER
	K = LAYER_TO_VARY

C  Save some calculation time. These quantities are always required.

	IF ( DO_CONSTANTS_SET ) THEN
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO AA = 1, NSTREAMS
	      LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
	      MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
	      DO Q = 1, K_PARAMETERS
	        NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XPOS(UM,AA,N)
	        PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XNEG(UM,AA,N)
  	      ENDDO
	    ENDDO
	  ENDDO
	ENDIF

C  Get the Multipliers
C  ===================

C  Get the external multiplier (results will be saved)
C  Only exists for N greater than or equal to K

	IF ( FOURIER .EQ. 0 ) THEN
	  IF ( N.GE.K ) THEN
	    CALL  L_WHOLELAYER_EMULT_UP ( N, K, K_PARAMETERS )
	  ENDIF
	ENDIF

C  Get the homogeneous solution multipliers (local output)
C  These are only required for the varying layer

	IF ( N .EQ. K ) THEN
	  CALL L_WHOLELAYER_HMULT_UP ( N, K_PARAMETERS )
	ENDIF

C  Get the Green's function solution multipliers 
C  Only exist for N greater than or equal to K

	IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN
	  IF ( N.GE.K ) THEN
	    CALL L_WHOLELAYER_GMULT_UP
     I       ( N, K, K_PARAMETERS )
	  ENDIF
	ENDIF

C  Homogeneous solutions
C  =====================

C  Special case when N = K
C  -----------------------

	IF ( N .EQ. K ) THEN

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO Q = 1, K_PARAMETERS
	      SHOM = ZERO
	      DO AA = 1, NSTREAMS
                H1 = LCON_UXVEC(UM,AA)   * L_HMULT_2(AA,UM,N,Q)
                H2 = NCON_UXVEC(UM,AA,Q) *   HMULT_2(AA,UM,N)
                H3 = LCON(AA,N)*L_U_XPOS(UM,AA,N,Q)*HMULT_2(AA,UM,N)
                H4 = MCON_UXVEC(UM,AA) * L_HMULT_1(AA,UM,N,Q)
                H5 = PCON_UXVEC(UM,AA,Q) *   HMULT_1(AA,UM,N)
                H6 =  MCON(AA,N)*L_U_XNEG(UM,AA,N,Q)*HMULT_1(AA,UM,N)
	        SHOM = SHOM + H1 + H2 + H3 + H4 + H5 + H6
	      ENDDO
	      L_LAYERSOURCE(UM,Q) = SHOM
	    ENDDO
	  ENDDO

C  Other cases when N not equal to K (only variation of Integ-Cons)
C  ---------------------------------

	ELSE

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO Q = 1, K_PARAMETERS
	      SHOM = ZERO
	      DO AA = 1, NSTREAMS
                H2 = NCON_UXVEC(UM,AA,Q) * HMULT_2(AA,UM,N)
                H5 = PCON_UXVEC(UM,AA,Q) * HMULT_1(AA,UM,N)
	        SHOM = SHOM + H2 + H5
	      ENDDO
	      L_LAYERSOURCE(UM,Q) = SHOM
	    ENDDO
	  ENDDO

	ENDIF

C  Particular and single scatter contributions - classical solution
C  ================================================================

	IF ( DO_CLASSICAL_SOLUTION ) THEN

C  Special case when N = K
C  -----------------------

	  IF ( N .EQ. K ) THEN

C  add particular solution

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = L_U_WPOS(UM,N,N,Q) * EMULT_UP(UM,N) +
     &                   U_WPOS2(UM,N) * L_EMULT_UP(UM,N,K,Q)
	        L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
	      ENDDO
	    ENDDO

C  Add single scatter term if flagged

	    IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        DO Q = 1, K_PARAMETERS
	          SFOR = L_U_WFORPOS(UM,N,Q) * EMULT_UP(UM,N) +
     &                   U_WPOS1(UM,N) * L_EMULT_UP(UM,N,K,Q)
	          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
	        ENDDO
	      ENDDO
	    ENDIF

C  Other cases when N > K
C  ----------------------

	  ELSE IF ( N .GT. K ) THEN

C  add particular solution

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = L_U_WPOS(UM,N,K,Q) *   EMULT_UP(UM,N) +
     &                      U_WPOS2(UM,N) * L_EMULT_UP(UM,N,K,Q)
	        L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
	      ENDDO
	    ENDDO

C  Add single scatter term if flagged

	    IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        DO Q = 1, K_PARAMETERS
	          SFOR = U_WPOS1(UM,N) * L_EMULT_UP(UM,N,K,Q)
	          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
	        ENDDO
	      ENDDO
	    ENDIF

	  ENDIF

C  Particular and single scatter contributions - Green's function solution
C  =======================================================================

	ELSE

C  Special case when N = K
C  -----------------------

	  IF ( N .EQ. K ) THEN

C  Add particular solution

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = ZERO
	        DO AA = 1, NSTREAMS
                  H1 =   U_XPOS(UM,AA,N)   * L_SGMULT_DN(AA,UM,Q)
                  H2 = L_U_XPOS(UM,AA,N,Q) *   SGMULT_UP_DN(AA,UM,N)
                  H3 =   U_XNEG(UM,AA,N)   * L_SGMULT_UP(AA,UM,Q)
                  H4 = L_U_XNEG(UM,AA,N,Q) *   SGMULT_UP_UP(AA,UM,N)
	          SPAR = SPAR + H1 + H2 + H3 + H4
	        ENDDO
	        L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
	      ENDDO
	    ENDDO

C  Add single scatter term if flagged

	    IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        DO Q = 1, K_PARAMETERS
	          SFOR = L_U_WFORPOS(UM,N,Q) * EMULT_UP(UM,N) +
     &                      U_WPOS1(UM,N)  * L_EMULT_UP(UM,N,K,Q)
	          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
	        ENDDO
	      ENDDO
	    ENDIF

C  Other cases when N > K
C  ----------------------

	  ELSE IF ( N .GT. K ) THEN

C  add particular solution

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = ZERO
	        DO AA = 1, NSTREAMS
                  H1 = U_XPOS(UM,AA,N) * L_SGMULT_DN(AA,UM,Q)
                  H3 = U_XNEG(UM,AA,N) * L_SGMULT_UP(AA,UM,Q)
	          SPAR = SPAR + H1 + H3
	        ENDDO
	        L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
	      ENDDO
	    ENDDO

C  Add single scatter term if flagged

	    IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        DO Q = 1, K_PARAMETERS
	          SFOR = U_WPOS1(UM,N)  * L_EMULT_UP(UM,N,K,Q)
	          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
	       ENDDO
	      ENDDO
	    ENDIF

	  ENDIF

	ENDIF

C  If operating in Ms-mode only, copy multiple scatter term for MSST
C  -----------------------------------------------------------------

	IF ( DO_MSMODE_LIDORT ) THEN
	  IF ( SAVE_LAYER_MSST ) THEN
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        L_MSCATUSE_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q)
	      ENDDO
	    ENDDO
	  ENDIF
	ENDIF

C  Finish

	END

C

	SUBROUTINE L_WHOLELAYER_STERM_DN
     I       ( FOURIER, GIVEN_LAYER,
     I         LAYER_TO_VARY, K_PARAMETERS,
     I         DO_CONSTANTS_SET,
     O         L_LAYERSOURCE, L_MSCATUSE_LAYERSOURCE )

C  source terms

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of multiplier and solution variables (input)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized multiplier and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  layer and Fourier number

	INTEGER		 GIVEN_LAYER, FOURIER

C  constants set flag

	LOGICAL		 DO_CONSTANTS_SET

C  output linearized layer source terms

	DOUBLE PRECISION
     &        L_LAYERSOURCE(MAX_USER_STREAMS,MAX_PARAMETERS),
     &        L_MSCATUSE_LAYERSOURCE(MAX_USER_STREAMS,MAX_PARAMETERS)

C  Linearization control

	INTEGER		 LAYER_TO_VARY, K_PARAMETERS

C  local variables

	INTEGER		 N, K, AA, UM, Q
	DOUBLE PRECISION SHOM, SFOR, SPAR, H1, H2, H3, H4, H5, H6

C  layer variables

	N = GIVEN_LAYER
	K = LAYER_TO_VARY

C  Save some calculation time. These quantities are always required.

	IF ( DO_CONSTANTS_SET ) THEN
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO AA = 1, NSTREAMS
	      LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
	      MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
	      DO Q = 1, K_PARAMETERS
	        NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XNEG(UM,AA,N)
	        PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XPOS(UM,AA,N)
  	      ENDDO
	    ENDDO
	  ENDDO
	ENDIF

C  Get the Multipliers
C  ===================

C  Get the external multiplier (results will be saved)
C  Only exists for N greater than or equal to K

	IF ( FOURIER .EQ. 0 ) THEN
	  IF ( N.GE.K ) THEN
	    CALL  L_WHOLELAYER_EMULT_DN ( N, K, K_PARAMETERS )
	  ENDIF
	ENDIF

C  Get the homogeneous solution multipliers (local output)
C  These are only required for the varying layer

	IF ( N .EQ. K ) THEN
	  CALL L_WHOLELAYER_HMULT_DN ( N, K_PARAMETERS )
	ENDIF

C  Get the Green's function solution multipliers 
C  Only exist for N greater than or equal to K

	IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN
	  IF ( N.GE.K ) THEN
	    CALL L_WHOLELAYER_GMULT_DN
     I       ( N, K, K_PARAMETERS )
	  ENDIF
	ENDIF

C  Homogeneous solutions
C  =====================

C  Special case when N = K
C  -----------------------

	IF ( N .EQ. K ) THEN

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO Q = 1, K_PARAMETERS
	      SHOM = ZERO
	      DO AA = 1, NSTREAMS
                H1 = LCON_UXVEC(UM,AA)   * L_HMULT_1(AA,UM,N,Q)
                H2 = NCON_UXVEC(UM,AA,Q) *   HMULT_1(AA,UM,N)
                H3 = LCON(AA,N)*L_U_XNEG(UM,AA,N,Q)*HMULT_1(AA,UM,N)
                H4 = MCON_UXVEC(UM,AA)   * L_HMULT_2(AA,UM,N,Q)
                H5 = PCON_UXVEC(UM,AA,Q) *   HMULT_2(AA,UM,N)
                H6 = MCON(AA,N)*L_U_XPOS(UM,AA,N,Q)*HMULT_2(AA,UM,N)
	        SHOM = SHOM + H1 + H2 + H3 + H4 + H5 + H6
	      ENDDO
	      L_LAYERSOURCE(UM,Q) = SHOM
	    ENDDO
	  ENDDO

C  Other cases when N not equal to K (only variation of Integ-Cons)
C  ---------------------------------

	ELSE

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO Q = 1, K_PARAMETERS
	      SHOM = ZERO
	      DO AA = 1, NSTREAMS
                H2 = NCON_UXVEC(UM,AA,Q) * HMULT_1(AA,UM,N)
                H5 = PCON_UXVEC(UM,AA,Q) * HMULT_2(AA,UM,N)
	        SHOM = SHOM + H2 + H5
	      ENDDO
	      L_LAYERSOURCE(UM,Q) = SHOM
	      L_MSCATUSE_LAYERSOURCE(UM,Q) = SHOM
	    ENDDO
	  ENDDO

	ENDIF

C  Particular and single scatter contributions - classical solution
C  ================================================================

	IF ( DO_CLASSICAL_SOLUTION ) THEN

C  Special case when N = K
C  -----------------------

	  IF ( N .EQ. K ) THEN

C  add particular solution

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = L_U_WNEG(UM,N,N,Q) * EMULT_DN(UM,N) +
     &                   U_WNEG2(UM,N) * L_EMULT_DN(UM,N,K,Q)
	        L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
	      ENDDO
	    ENDDO

C  Add single scatter term if flagged

	    IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        DO Q = 1, K_PARAMETERS
	          SFOR = L_U_WFORNEG(UM,N,Q) * EMULT_DN(UM,N) +
     &                     U_WNEG1(UM,N) * L_EMULT_DN(UM,N,K,Q)
	          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
	        ENDDO
	      ENDDO
	    ENDIF

C  Other cases when N > K
C  ----------------------

	  ELSE IF ( N .GT. K ) THEN

C  add particular solution

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = L_U_WNEG(UM,N,K,Q) *   EMULT_DN(UM,N) +
     &                      U_WNEG2(UM,N) * L_EMULT_DN(UM,N,K,Q)
	        L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
	      ENDDO
	    ENDDO

C  Add single scatter term if flagged

	    IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        DO Q = 1, K_PARAMETERS
	          SFOR = U_WNEG1(UM,N) * L_EMULT_DN(UM,N,K,Q)
 	          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
	        ENDDO
	      ENDDO
	    ENDIF

	  ENDIF

C  Particular and single scatter contributions - Green's function solution
C  =======================================================================

	ELSE

C  Special case when N = K
C  -----------------------

	  IF ( N .EQ. K ) THEN

C  Add particular solution

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = ZERO
	        DO AA = 1, NSTREAMS
                  H1 =   U_XNEG(UM,AA,N)   * L_SGMULT_DN(AA,UM,Q)
                  H2 = L_U_XNEG(UM,AA,N,Q) *   SGMULT_DN_DN(AA,UM,N)
                  H3 =   U_XPOS(UM,AA,N)   * L_SGMULT_UP(AA,UM,Q)
                  H4 = L_U_XPOS(UM,AA,N,Q) *   SGMULT_DN_UP(AA,UM,N)
	          SPAR = SPAR + H1 + H2 + H3 + H4
	        ENDDO
	        L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
	      ENDDO
	    ENDDO

C  Add single scatter term if flagged

	    IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        DO Q = 1, K_PARAMETERS
	          SFOR = L_U_WFORNEG(UM,N,Q) * EMULT_DN(UM,N) +
     &                      U_WNEG1(UM,N)  * L_EMULT_DN(UM,N,K,Q)
	        L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
	        ENDDO
	      ENDDO
	    ENDIF

C  Other cases when N > K
C  ----------------------

C  add particular solution

	  ELSE IF ( N .GT. K ) THEN

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = ZERO
	        DO AA = 1, NSTREAMS
                  H1 =   U_XNEG(UM,AA,N)   * L_SGMULT_DN(AA,UM,Q)
                  H3 =   U_XPOS(UM,AA,N)   * L_SGMULT_UP(AA,UM,Q)
	          SPAR = SPAR + H1 + H3
	        ENDDO
	        L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
	      ENDDO
	    ENDDO

C  Add single scatter term if flagged

	    IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        DO Q = 1, K_PARAMETERS
	          SFOR = U_WNEG1(UM,N)  * L_EMULT_DN(UM,N,K,Q)
	          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
	        ENDDO
	      ENDDO
	    ENDIF

	  ENDIF

	ENDIF

C  If operating in Ms-mode only, copy multiple scatter term for MSST
C  -----------------------------------------------------------------

	IF ( DO_MSMODE_LIDORT ) THEN
	  IF ( SAVE_LAYER_MSST ) THEN
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        L_MSCATUSE_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q)
	      ENDDO
	    ENDDO
	  ENDIF
	ENDIF

C  Finish

	END

C

	SUBROUTINE L_PARTLAYER_STERM_UP
     I       ( FOURIER, GIVEN_LAYER, OFFGRID_INDEX,
     I         LAYER_TO_VARY, K_PARAMETERS,
     O         L_LAYERSOURCE )

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of multiplier and solution variables (input)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized multiplier and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  layer and Fourier number, offgrid index

	INTEGER		 GIVEN_LAYER, FOURIER, OFFGRID_INDEX

C  output linearized layer source term

	DOUBLE PRECISION L_LAYERSOURCE(MAX_USER_STREAMS,MAX_PARAMETERS)

C  Linearization control

	INTEGER		 LAYER_TO_VARY, K_PARAMETERS

C  local variables

	INTEGER		 N, K, AA, UM, Q, UT
	DOUBLE PRECISION SHOM, SFOR, SPAR, H1, H2, H3, H4, H5, H6

C  layer variables

	N = GIVEN_LAYER
	K = LAYER_TO_VARY
	UT = OFFGRID_INDEX

C  Save some calculation time. These quantities are always required.

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  DO AA = 1, NSTREAMS
	    LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
	    MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
	    DO Q = 1, K_PARAMETERS
	      NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XPOS(UM,AA,N)
	      PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XNEG(UM,AA,N)
  	    ENDDO
	  ENDDO
	ENDDO

C  Multipliers
C  ===========

C  Get the external multiplier (results will be saved)
C  Only exists for N greater than or equal to K

	IF ( FOURIER .EQ. 0 ) THEN
	  IF ( N.GE.K ) THEN
	    CALL  L_PARTLAYER_EMULT_UP ( N, UT, K, K_PARAMETERS )
	  ENDIF
	ENDIF

C  Get the homogeneous solution multipliers (local output)
C  These are only required for the varying layer

	IF ( N .EQ. K ) THEN
	  CALL L_PARTLAYER_HMULT_UP ( N, UT, K_PARAMETERS )
	ENDIF

C  Get the Green's function solution multipliers
C  Only exist for N greater than or equal to K

	IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN
	  IF ( N.GE.K ) THEN
	    CALL L_PARTLAYER_GMULT_UP
     I       ( N, UT, K, K_PARAMETERS )
	  ENDIF
	ENDIF

C  Partial layer source function ( Homogeneous/constants variation )
C  =================================================================

C  Special case when N = K
C  -----------------------

	IF ( N .EQ. K ) THEN

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO Q = 1, K_PARAMETERS
	      SHOM = ZERO
	      DO AA = 1, NSTREAMS
                H1 = LCON_UXVEC(UM,AA)   * L_UT_HMULT_UP_DN(AA,UM,UT,Q)
                H2 = NCON_UXVEC(UM,AA,Q) *   UT_HMULT_UP_DN(AA,UM,UT)
                H3 = LCON(AA,N) * L_U_XPOS(UM,AA,N,Q)
                H3 = H3 * UT_HMULT_UP_DN(AA,UM,UT)
                H4 = MCON_UXVEC(UM,AA)   * L_UT_HMULT_UP_UP(AA,UM,UT,Q)
                H5 = PCON_UXVEC(UM,AA,Q) *   UT_HMULT_UP_UP(AA,UM,UT)
                H6 = MCON(AA,N) * L_U_XNEG(UM,AA,N,Q)
	        H6 = H6 * UT_HMULT_UP_UP(AA,UM,UT)
	        SHOM = SHOM + H1 + H2 + H3 + H4 + H5 + H6
	      ENDDO
	      L_LAYERSOURCE(UM,Q) = SHOM
	    ENDDO
	  ENDDO

C  Other cases when N not equal to K (only variation of Integ-Cons)
C  ---------------------------------

	ELSE 

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO Q = 1, K_PARAMETERS
	      SHOM = ZERO
	      DO AA = 1, NSTREAMS
                H2 = NCON_UXVEC(UM,AA,Q) * UT_HMULT_UP_DN(AA,UM,UT)
                H5 = PCON_UXVEC(UM,AA,Q) * UT_HMULT_UP_UP(AA,UM,UT)
	        SHOM = SHOM + H2 + H5
	      ENDDO
	      L_LAYERSOURCE(UM,Q) = SHOM
	    ENDDO
	  ENDDO

	ENDIF

C  Partial layer source function ( SS/Particular, Classical solution )
C  ===================================================================

	IF ( DO_CLASSICAL_SOLUTION ) THEN

C  Special case when N = K
C  -----------------------

	  IF ( N .EQ. K ) THEN

C  particular solution

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = L_U_WPOS(UM,N,N,Q) * UT_EMULT_UP(UM,UT) +
     &                   U_WPOS2(UM,N) * L_UT_EMULT_UP(UM,UT,K,Q)
	        L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
	      ENDDO
	    ENDDO

C  Add single scatter term if flagged

	    IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        DO Q = 1, K_PARAMETERS
	          SFOR = L_U_WFORPOS(UM,N,Q) * UT_EMULT_UP(UM,UT) +
     &                     U_WPOS1(UM,N) * L_UT_EMULT_UP(UM,UT,K,Q)
	          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
	        ENDDO
	      ENDDO
	    ENDIF

C  Other cases when N > K
C  -----------------------

	  ELSE IF ( N .GT. K ) THEN

C  particular solution

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = L_U_WPOS(UM,N,K,Q) *   UT_EMULT_UP(UM,UT) +
     &                   U_WPOS2(UM,N)    * L_UT_EMULT_UP(UM,UT,K,Q)
	        L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
	      ENDDO
	    ENDDO

C  Add single scatter term if flagged

	    IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        DO Q = 1, K_PARAMETERS
	          SFOR = U_WPOS1(UM,N) * L_UT_EMULT_UP(UM,UT,K,Q)
	          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
	        ENDDO
	      ENDDO
	    ENDIF

	  ENDIF

C  Partial layer source function ( SS/Particular, Green's solution )
C  =================================================================

	ELSE

C  Special case when N = K
C  -----------------------

	  IF ( N .EQ. K ) THEN

C  particular solution

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = ZERO
	        DO AA = 1, NSTREAMS
                  H1 =   U_XPOS(UM,AA,N)   *  L_SGMULT_DN(AA,UM,Q)
                  H2 = L_U_XPOS(UM,AA,N,Q) * UT_SGMULT_UP_DN(AA,UM,UT)
                  H3 =   U_XNEG(UM,AA,N)   *  L_SGMULT_UP(AA,UM,Q)
                  H4 = L_U_XNEG(UM,AA,N,Q) * UT_SGMULT_UP_UP(AA,UM,UT)
	          SPAR = SPAR + H1 + H2 + H3 + H4
	        ENDDO
	        L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
	      ENDDO
	    ENDDO

C  Add single scatter term if flagged

	    IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        DO Q = 1, K_PARAMETERS
	          SFOR = L_U_WFORPOS(UM,N,Q) * UT_EMULT_UP(UM,UT) +
     &                        U_WPOS1(UM,N)  * L_UT_EMULT_UP(UM,UT,K,Q)
	          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
	        ENDDO
	      ENDDO
	    ENDIF

C  Other cases when N > K
C  -----------------------

	  ELSE IF ( N .GT. K ) THEN

C  particular solution

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = ZERO
	        DO AA = 1, NSTREAMS
                  H1 = U_XPOS(UM,AA,N) * L_SGMULT_DN(AA,UM,Q)
                  H3 = U_XNEG(UM,AA,N) * L_SGMULT_UP(AA,UM,Q)
	          SPAR = SPAR + H1 + H3
	        ENDDO
	        L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
	      ENDDO
	    ENDDO

C  Add single scatter term if flagged

	    IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        DO Q = 1, K_PARAMETERS
	          SFOR = U_WPOS1(UM,N)  * L_UT_EMULT_UP(UM,UT,K,Q)
	          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
	        ENDDO
	      ENDDO
	    ENDIF

	  ENDIF

	ENDIF

C  Finish

	END

C

	SUBROUTINE L_PARTLAYER_STERM_DN
     I       ( FOURIER, GIVEN_LAYER, OFFGRID_INDEX,
     I         LAYER_TO_VARY, K_PARAMETERS,
     O         L_LAYERSOURCE )

C  source terms

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of multiplier and solution variables (input)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized multiplier and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  layer and Fourier number, offgrid index

	INTEGER		 GIVEN_LAYER, FOURIER, OFFGRID_INDEX

C  output linearized layer source term

	DOUBLE PRECISION L_LAYERSOURCE(MAX_USER_STREAMS,MAX_PARAMETERS)

C  Linearization control

	INTEGER		 LAYER_TO_VARY, K_PARAMETERS

C  local variables

	INTEGER		 N, K, AA, UM, Q, UT
	DOUBLE PRECISION SHOM, SFOR, SPAR, H1, H2, H3, H4, H5, H6

C  layer variables

	N = GIVEN_LAYER
	K = LAYER_TO_VARY
	UT = OFFGRID_INDEX

C  Save some calculation time. These quantities are always required.

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  DO AA = 1, NSTREAMS
	    LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
	    MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
	    DO Q = 1, K_PARAMETERS
	      NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XNEG(UM,AA,N)
	      PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XPOS(UM,AA,N)
  	    ENDDO
	  ENDDO
	ENDDO

C  Multipliers
C  ===========

C  Get the external multiplier (results will be saved)
C  Only exists for N greater than or equal to K

	IF ( FOURIER .EQ. 0 ) THEN
	  IF ( N.GE.K ) THEN
	    CALL  L_PARTLAYER_EMULT_DN ( N, UT, K, K_PARAMETERS )
	  ENDIF
	ENDIF

C  Get the homogeneous solution multipliers (local output)
C  These are only required for the varying layer

	IF ( N .EQ. K ) THEN
	  CALL L_PARTLAYER_HMULT_DN ( N, UT, K_PARAMETERS )
	ENDIF

C  Get the Green's function solution multipliers
C  Only exist for N greater than or equal to K

	IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN
	  IF ( N.GE.K ) THEN
	    CALL L_PARTLAYER_GMULT_DN
     I       ( N, UT, K, K_PARAMETERS )
	  ENDIF
	ENDIF

C  Partial layer source function ( Homogeneous/constants variation )
C  =================================================================

C  Special case when N = K
C  -----------------------

	IF ( N .EQ. K ) THEN

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO Q = 1, K_PARAMETERS
	      SHOM = ZERO
	      DO AA = 1, NSTREAMS
                H1 = LCON_UXVEC(UM,AA)   * L_UT_HMULT_DN_DN(AA,UM,UT,Q)
                H2 = NCON_UXVEC(UM,AA,Q) *   UT_HMULT_DN_DN(AA,UM,UT)
                H3 = LCON(AA,N) * L_U_XNEG(UM,AA,N,Q)
                H3 = H3 * UT_HMULT_DN_DN(AA,UM,UT)
                H4 = MCON_UXVEC(UM,AA)   * L_UT_HMULT_DN_UP(AA,UM,UT,Q)
                H5 = PCON_UXVEC(UM,AA,Q) *   UT_HMULT_DN_UP(AA,UM,UT)
                H6 = MCON(AA,N) * L_U_XPOS(UM,AA,N,Q)
	        H6 = H6 * UT_HMULT_DN_UP(AA,UM,UT)
	        SHOM = SHOM + H1 + H2 + H3 + H4 + H5 + H6
	      ENDDO
	      L_LAYERSOURCE(UM,Q) = SHOM + SFOR
	    ENDDO
	  ENDDO

C  Other cases when N not equal to K (only variation of Integ-Cons)
C  ---------------------------------

	ELSE 

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO Q = 1, K_PARAMETERS
	      SHOM = ZERO
	      DO AA = 1, NSTREAMS
                H2 = NCON_UXVEC(UM,AA,Q) * UT_HMULT_DN_DN(AA,UM,UT)
                H5 = PCON_UXVEC(UM,AA,Q) * UT_HMULT_DN_UP(AA,UM,UT)
	        SHOM = SHOM + H2 + H5
	      ENDDO
	      L_LAYERSOURCE(UM,Q) = SHOM
	    ENDDO
	  ENDDO

	ENDIF

C  Partial layer source function ( SS/Particular, Classical solution )
C  ===================================================================

	IF ( DO_CLASSICAL_SOLUTION ) THEN

C  Special case when N = K
C  -----------------------

	  IF ( N .EQ. K ) THEN

C  particular solution

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = L_U_WNEG(UM,N,N,Q) * UT_EMULT_DN(UM,UT) +
     &                   U_WNEG2(UM,N) * L_UT_EMULT_DN(UM,UT,K,Q)
	        L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
	      ENDDO
	    ENDDO

C  Add single scatter term if flagged

	    IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        DO Q = 1, K_PARAMETERS
	          SFOR = L_U_WFORNEG(UM,N,Q) * UT_EMULT_DN(UM,UT) +
     &                     U_WNEG1(UM,N) * L_UT_EMULT_DN(UM,UT,K,Q)
	          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
	        ENDDO
	      ENDDO
	    ENDIF

C  Other cases when N > K
C  -----------------------

	  ELSE IF ( N .GT. K ) THEN

C  particular solution

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = L_U_WNEG(UM,N,K,Q) *   UT_EMULT_DN(UM,UT) +
     &                   U_WNEG2(UM,N)    * L_UT_EMULT_DN(UM,UT,K,Q)
	        L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
	      ENDDO
	    ENDDO

C  Add single scatter term if flagged

	    IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        DO Q = 1, K_PARAMETERS
	          SFOR = U_WNEG1(UM,N) * L_UT_EMULT_DN(UM,UT,K,Q)
	          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
	        ENDDO
	      ENDDO
	    ENDIF

	  ENDIF

C  Partial layer source function ( SS/Particular, Green's solution )
C  =================================================================

	ELSE

C  Special case when N = K
C  -----------------------

	  IF ( N .EQ. K ) THEN

C  particular solution

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = ZERO
	        DO AA = 1, NSTREAMS
                  H1 =   U_XNEG(UM,AA,N)   *  L_SGMULT_DN(AA,UM,Q)
                  H2 = L_U_XNEG(UM,AA,N,Q) * UT_SGMULT_DN_DN(AA,UM,UT)
                  H3 =   U_XPOS(UM,AA,N)   *  L_SGMULT_UP(AA,UM,Q)
                  H4 = L_U_XPOS(UM,AA,N,Q) * UT_SGMULT_DN_UP(AA,UM,UT)
	          SPAR = SPAR + H1 + H2 + H3 + H4
	        ENDDO
	        L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
	      ENDDO
	    ENDDO

C  Add single scatter term if flagged

	    IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        DO Q = 1, K_PARAMETERS
	          SFOR = L_U_WFORNEG(UM,N,Q) * UT_EMULT_DN(UM,UT) +
     &                        U_WNEG1(UM,N)  * L_UT_EMULT_DN(UM,UT,K,Q)
	          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
	        ENDDO
	      ENDDO
	    ENDIF

C  Other cases when N > K
C  -----------------------

	  ELSE IF ( N .GT. K ) THEN

C  particular solution

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        SPAR = ZERO
	        DO AA = 1, NSTREAMS
                  H1 =   U_XNEG(UM,AA,N)   * L_SGMULT_DN(AA,UM,Q)
                  H3 =   U_XPOS(UM,AA,N)   * L_SGMULT_UP(AA,UM,Q)
	          SPAR = SPAR + H1 + H3
	        ENDDO
	        L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR
	      ENDDO
	    ENDDO

C  Add single scatter term if flagged

	    IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        DO Q = 1, K_PARAMETERS
	          SFOR = U_WNEG1(UM,N)  * L_UT_EMULT_DN(UM,UT,K,Q)
	          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
	        ENDDO
	      ENDDO
	    ENDIF

	  ENDIF

	ENDIF

C  Finish

	END

C

	SUBROUTINE UPUSER_ALBEDOWF
     I    ( DO_INCLUDE_SURFEMISS, DO_INCLUDE_MVOUTPUT,
     I      R2, F1, FOURIER, USER_DIRECT_BEAM )

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'
	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Flux factor

	DOUBLE PRECISION F1

C  local control flags

	LOGICAL		 DO_INCLUDE_SURFEMISS
	LOGICAL		 DO_INCLUDE_MVOUTPUT

C  albedo factor

	DOUBLE PRECISION R2

C  Fourier number

	INTEGER		 FOURIER

C  Direct beam

	DOUBLE PRECISION USER_DIRECT_BEAM(MAX_USER_STREAMS)

C  local variables
C  ---------------

	INTEGER		 N, NUT, NSTART, NUT_PREV, NLEVEL, NL
	INTEGER		 UTA, UM, IUM, IQD, I, I1, J, UT, AA
	DOUBLE PRECISION L_CUMUL_SOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION L_BOA_SOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION L_LAYER_SOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION NCON_HELP(MAX_USER_STREAMS,MAXSTRM)
	DOUBLE PRECISION PCON_HELP(MAX_USER_STREAMS,MAXSTRM)
	DOUBLE PRECISION HELP(MAXSTRM), REFLEC, L_DOWN
	DOUBLE PRECISION L_FINAL_SOURCE, H1, H2, SHOM
	LOGICAL		 DO_CONSTANTS_SET(MAXLAYER)

C  Zero all Fourier components close to zenith

	IF ( DO_USER_STREAMS ) THEN
	  IF ( LOCAL_UM_START .GT. 1 ) THEN
	    DO UTA = 1, N_OUT_USERTAUS
	      DO UM = 1, LOCAL_UM_START - 1
	        IUM = USEROUTPUT_INDEX(UM)
	        ALBEDOWF_F(UTA,IUM,UPIDX) = ZERO
	      ENDDO
	    ENDDO
	  ENDIF
	ENDIF

C  get the BOA source term
C  -----------------------

	IF ( DO_USER_STREAMS ) THEN

C  initialise boa source function

	  N = NLAYERS
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    L_BOA_SOURCE(UM) = ZERO
	  ENDDO

C  reflectance from surface
C  .. Downward intensity at computational angles (homog)
C  .. reflectance integrand  a(j).x(j).L_DOWN(-j)

	  DO I = 1, NSTREAMS
	    L_DOWN = ZERO
	    DO AA = 1, NSTREAMS
	      H1 = NCON_ALB(AA,N)*XPOS(I,AA,N)*T_DELT_EIGEN(AA,N)
	      H2 = PCON_ALB(AA,N)*XNEG(I,AA,N) 
	      L_DOWN = L_DOWN + H1 + H2
	    ENDDO
	    HELP(I) = L_DOWN * AX(I)
	  ENDDO

C  reflected multiple scatter intensity at user defined-angles
C  .. integrate reflectance, same for all user-streams in Lambertian case
C  .. integrate with reflectance function at user angles (non-Lambertian)

	  IF ( DO_LAMBERTIAN_ALBEDO ) THEN
	    REFLEC = ZERO
	    DO J = 1, NSTREAMS
	       REFLEC = REFLEC + HELP(J)
	    ENDDO
	    REFLEC = R2 * REFLEC
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      L_BOA_SOURCE(UM) = REFLEC
	    ENDDO
	  ELSE
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      REFLEC = ZERO
	      DO J = 1, NSTREAMS
	        REFLEC = REFLEC + HELP(J)*USER_BIREFLEC(FOURIER,UM,J)
	      ENDDO
	      L_BOA_SOURCE(UM) = R2 * REFLEC
	    ENDDO
	  ENDIF
	  
C  Add linearization term due to albedo variation 

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    L_BOA_SOURCE(UM) = L_BOA_SOURCE(UM) +
     &                  BOA_SOURCE(UM) + USER_DIRECT_BEAM(UM)
	  ENDDO

C  Add emissivity variation at user defined angles
C  (expression for emissivity variation follows from Kirchhoff's law)

	  IF ( DO_INCLUDE_SURFEMISS ) THEN
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      L_BOA_SOURCE(UM) = L_BOA_SOURCE(UM) - SURFBB
	    ENDDO
	  ENDIF

C  Initialize recursion for user-defined stream angles only

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    L_CUMUL_SOURCE(UM) = L_BOA_SOURCE(UM)
	  ENDDO

C  Save multiple scatter term output

	  IF ( SAVE_LAYER_MSST ) THEN
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      ALBEDOWF_MSCATBOA_STERM_F(UM) = 
     &           F1*L_BOA_SOURCE(UM)
	    ENDDO
	  ENDIF

	ENDIF

C  initialise cumulative source term loop

 	NSTART = NLAYERS
	NUT_PREV = NSTART + 1
	DO N = 1, NLAYERS
	  DO_CONSTANTS_SET(N) = .TRUE.
	ENDDO

C  loop over all output optical depths
C  -----------------------------------

	DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

	  NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)

	  IF ( DO_USER_STREAMS ) THEN
	    NUT = NLEVEL + 1
	    DO N = NSTART, NUT, -1

C  Save some calculation time. These quantities are always required.

	      IF ( DO_CONSTANTS_SET(N) ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          DO AA = 1, NSTREAMS
	            NCON_HELP(UM,AA) =  NCON_ALB(AA,N) * U_XPOS(UM,AA,N)
	            PCON_HELP(UM,AA) =  PCON_ALB(AA,N) * U_XNEG(UM,AA,N)
  	          ENDDO
	        ENDDO
 	      ENDIF

	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
                  H1 = NCON_HELP(UM,AA) * HMULT_2(AA,UM,N)
                  H2 = PCON_HELP(UM,AA) * HMULT_1(AA,UM,N)
	          SHOM = SHOM + H1 + H2
	        ENDDO
	        L_LAYER_SOURCE(UM) = SHOM
	      ENDDO

	      IF ( SAVE_LAYER_MSST ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          ALBEDOWF_MSCATSTERM_F(N,UM,UPIDX) =
     &                L_LAYER_SOURCE(UM) * F1
	        ENDDO
	      ENDIF

	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        L_CUMUL_SOURCE(UM) = L_LAYER_SOURCE(UM)  +
     &                   T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM)
	      ENDDO

	    ENDDO
	  ENDIF

C  Offgrid output
C  --------------

	  IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

	    UT = OFFGRID_UTAU_OUTINDEX(UTA)
	    N  = OFFGRID_UTAU_LAYERIDX(UT)

C  Quadrature output, mean value + flux output

	    IF ( DO_QUAD_OUTPUT. OR. DO_INCLUDE_MVOUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        I1 = I + NSTREAMS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          H1 = NCON_ALB(AA,N)*XPOS(I1,AA,N)*T_UTDN_EIGEN(AA,UT)
	          H2 = PCON_ALB(AA,N)*XNEG(I1,AA,N)*T_UTUP_EIGEN(AA,UT)
	          SHOM = SHOM + H1 + H2
	        ENDDO
	        QUAD_ALBEDO_WF(UTA,I,UPIDX) = F1 * SHOM
	      ENDDO
	    ENDIF

C  Copy results if quadrature output required (unlikely)

	    IF ( DO_QUAD_OUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        IQD = QUADOUTPUT_INDEX(I)
 	        ALBEDOWF_F(UTA,IQD,UPIDX) = QUAD_ALBEDO_WF(UTA,I,UPIDX)
	      ENDDO
	    ENDIF

C  User-defined stream output

	    IF ( DO_USER_STREAMS ) THEN

C  set help arrays

	      IF ( DO_CONSTANTS_SET(N) ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          DO AA = 1, NSTREAMS
	            NCON_HELP(UM,AA) =  NCON_ALB(AA,N) * U_XPOS(UM,AA,N)
	            PCON_HELP(UM,AA) =  PCON_ALB(AA,N) * U_XNEG(UM,AA,N)
  	          ENDDO
	        ENDDO
 	      ENDIF
	      DO_CONSTANTS_SET(N) = .FALSE.

C  add additional partial layer source term

	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
                  H1 = NCON_HELP(UM,AA) * UT_HMULT_UP_DN(AA,UM,UT)
                  H2 = PCON_HELP(UM,AA) * UT_HMULT_UP_UP(AA,UM,UT)
	          SHOM = SHOM + H1 + H2
	        ENDDO
	        L_LAYER_SOURCE(UM) = SHOM
	      ENDDO

C  assign final albedo weighting function

	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        L_FINAL_SOURCE = L_LAYER_SOURCE(UM)  +
     &              T_UTUP_USERM(UT,UM)  * L_CUMUL_SOURCE(UM)
	        ALBEDOWF_F(UTA,IUM,UPIDX) = F1 * L_FINAL_SOURCE
	      ENDDO

	    ENDIF

C  Ongrid output
C  -------------

	  ELSE

C  Quadrature output, mean value + flux output

	    IF ( DO_QUAD_OUTPUT. OR. DO_INCLUDE_MVOUTPUT ) THEN

C  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
C  looking at the perturbation field at the top of these layers. The
C  case where the level mask = NLAYERS is the upwelling perturbed fields
C  at the bottom of the atmosphere (treated separately).

	      NL = NLEVEL
	      N = NL + 1

C  For the lowest level

	      IF ( NLEVEL .EQ. NLAYERS ) THEN

	        DO I = 1, NSTREAMS
	          I1 = I + NSTREAMS
	          SHOM = ZERO
	          DO AA = 1, NSTREAMS
	            H1 = NCON_ALB(AA,NL)*XPOS(I1,AA,NL)*T_DELT_EIGEN(AA,NL)
	            H2 = PCON_ALB(AA,NL)*XNEG(I1,AA,NL)
	            SHOM = SHOM + H1 + H2
	          ENDDO
 	          QUAD_ALBEDO_WF(UTA,I,UPIDX) = F1 * SHOM
	        ENDDO	  

C  For other levels in the atmosphere

	      ELSE

	        DO I = 1, NSTREAMS
	          I1 = I + NSTREAMS
	          SHOM = ZERO
	          DO AA = 1, NSTREAMS
	            H1 = NCON_ALB(AA,N)*XPOS(I1,AA,N)
	            H2 =
     &                PCON_ALB(AA,N)*XNEG(I1,AA,N)*T_DELT_EIGEN(AA,N)
	            SHOM = SHOM + H1 + H2
	          ENDDO
 	          QUAD_ALBEDO_WF(UTA,I,UPIDX) = F1 * SHOM
	        ENDDO

	      ENDIF

	    ENDIF

C  Copy results if quadrature output required (unlikely)

	    IF ( DO_QUAD_OUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        IQD = QUADOUTPUT_INDEX(I)
 	        ALBEDOWF_F(UTA,IQD,UPIDX) = QUAD_ALBEDO_WF(UTA,I,UPIDX)
	      ENDDO
	    ENDIF

C  User-defined stream output, just set to the cumulative source term

	    IF ( DO_USER_STREAMS ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        ALBEDOWF_F(UTA,IUM,UPIDX) = F1 * L_CUMUL_SOURCE(UM)
	      ENDDO
	    ENDIF

	  ENDIF

C  Check for updating the recursion

	  IF ( DO_USER_STREAMS ) THEN
	    IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
	    NUT_PREV = NUT
	  ENDIF

C  end loop over optical depth

	ENDDO

C  Finish

	END

C

	SUBROUTINE DNUSER_ALBEDOWF ( DO_INCLUDE_MVOUTPUT, F1 )

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'
	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Flux factor

	DOUBLE PRECISION F1

C  local inclusion flag for mean-value ouptut

	LOGICAL		 DO_INCLUDE_MVOUTPUT

C  local variables
C  ---------------

	INTEGER		 N, NUT, NSTART, NUT_PREV, NLEVEL, NL
	INTEGER		 UTA, UM, IUM, IQD, I, UT, AA
	DOUBLE PRECISION L_CUMUL_SOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION L_TOA_SOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION L_LAYER_SOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION NCON_HELP(MAX_USER_STREAMS,MAXSTRM)
	DOUBLE PRECISION PCON_HELP(MAX_USER_STREAMS,MAXSTRM)
	DOUBLE PRECISION L_FINAL_SOURCE, H1, H2, SHOM
	LOGICAL		 DO_CONSTANTS_SET(MAXLAYER)

C  Zero all Fourier components close to zenith

	IF ( DO_USER_STREAMS ) THEN
	  IF ( LOCAL_UM_START .GT. 1 ) THEN
	    DO UTA = 1, N_OUT_USERTAUS
	      DO UM = 1, LOCAL_UM_START - 1
	        IUM = USEROUTPUT_INDEX(UM)
	        ALBEDOWF_F(UTA,IUM,DNIDX) = ZERO
	      ENDDO
	    ENDDO
	  ENDIF
	ENDIF

C  get the TOA source term
C  -----------------------

	IF ( DO_USER_STREAMS ) THEN

C  initialise boa source function

	  N = NLAYERS
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    L_TOA_SOURCE(UM) = ZERO
	  ENDDO

C  Initialize recursion for user-defined stream angles only

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    L_CUMUL_SOURCE(UM) = L_TOA_SOURCE(UM)
	  ENDDO

	ENDIF

C  initialise cumulative source term loop

 	NSTART = 1
	NUT_PREV = NSTART - 1
	DO N = 1, NLAYERS
	  DO_CONSTANTS_SET(N) = .TRUE.
	ENDDO

C  loop over all output optical depths
C  -----------------------------------

	DO UTA = 1, N_OUT_USERTAUS

C  Layer index for given optical depth

	  NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)

	  IF ( DO_USER_STREAMS ) THEN
	    NUT = NLEVEL
	    DO N = NSTART, NUT

C  Save some calculation time. These quantities are always required.

	      IF ( DO_CONSTANTS_SET(N) ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          DO AA = 1, NSTREAMS
	            NCON_HELP(UM,AA) =  NCON_ALB(AA,N) * U_XNEG(UM,AA,N)
	            PCON_HELP(UM,AA) =  PCON_ALB(AA,N) * U_XPOS(UM,AA,N)
  	          ENDDO
	        ENDDO
 	      ENDIF

	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
                  H1 = NCON_HELP(UM,AA) * HMULT_1(AA,UM,N)
                  H2 = PCON_HELP(UM,AA) * HMULT_2(AA,UM,N)
	          SHOM = SHOM + H1 + H2
	        ENDDO
	        L_LAYER_SOURCE(UM) = SHOM
	      ENDDO

	      IF ( SAVE_LAYER_MSST ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          ALBEDOWF_MSCATSTERM_F(N,UM,DNIDX) =
     &                L_LAYER_SOURCE(UM) * F1
	        ENDDO
	      ENDIF

	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        L_CUMUL_SOURCE(UM) = L_LAYER_SOURCE(UM)  +
     &                   T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM)
	      ENDDO

	    ENDDO
	  ENDIF

C  Offgrid output
C  --------------

	  IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

	    UT = OFFGRID_UTAU_OUTINDEX(UTA)
	    N  = OFFGRID_UTAU_LAYERIDX(UT)

C  Quadrature output, mean value + flux output

	    IF ( DO_QUAD_OUTPUT. OR. DO_INCLUDE_MVOUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          H1 = NCON_ALB(AA,N)*XPOS(I,AA,N)*T_UTDN_EIGEN(AA,UT)
	          H2 = PCON_ALB(AA,N)*XNEG(I,AA,N)*T_UTUP_EIGEN(AA,UT)
	          SHOM = SHOM + H1 + H2
	        ENDDO
	        QUAD_ALBEDO_WF(UTA,I,DNIDX) = F1 * SHOM
	      ENDDO
	    ENDIF

C  Copy results if quadrature output required (unlikely)

	    IF ( DO_QUAD_OUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        IQD = QUADOUTPUT_INDEX(I)
 	        ALBEDOWF_F(UTA,IQD,DNIDX) = QUAD_ALBEDO_WF(UTA,I,DNIDX)
	      ENDDO
	    ENDIF

C  User-defined stream output

	    IF ( DO_USER_STREAMS ) THEN

C  set help arrays

	      IF ( DO_CONSTANTS_SET(N) ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          DO AA = 1, NSTREAMS
	            NCON_HELP(UM,AA) =  NCON_ALB(AA,N) * U_XNEG(UM,AA,N)
	            PCON_HELP(UM,AA) =  PCON_ALB(AA,N) * U_XPOS(UM,AA,N)
  	          ENDDO
	        ENDDO
 	      ENDIF
	      DO_CONSTANTS_SET(N) = .FALSE.

C  add additional partial layer source term

	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
                  H1 = NCON_HELP(UM,AA) * UT_HMULT_DN_DN(AA,UM,UT)
                  H2 = PCON_HELP(UM,AA) * UT_HMULT_DN_UP(AA,UM,UT)
	          SHOM = SHOM + H1 + H2
	        ENDDO
	        L_LAYER_SOURCE(UM) = SHOM
	      ENDDO

C  assign final albedo weighting function

	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        L_FINAL_SOURCE = L_LAYER_SOURCE(UM)  +
     &              T_UTDN_USERM(UT,UM)  * L_CUMUL_SOURCE(UM)
	        ALBEDOWF_F(UTA,IUM,DNIDX) = F1 * L_FINAL_SOURCE
	      ENDDO

	    ENDIF

C  Ongrid output
C  -------------

	  ELSE

C  Quadrature output, mean value + flux output

	    IF ( DO_QUAD_OUTPUT. OR. DO_INCLUDE_MVOUTPUT ) THEN

C  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
C  looking at the perturbation field at the top of these layers. The
C  case where the level mask = NLAYERS is the upwelling perturbed fields
C  at the bottom of the atmosphere (treated separately).

	      NL = NLEVEL
	      N = NL

C  For the highest level

	      IF ( NLEVEL .EQ. 0 ) THEN

	        DO I = 1, NSTREAMS
 	          QUAD_ALBEDO_WF(UTA,I,DNIDX) = ZERO
	        ENDDO	  

C  For other levels in the atmosphere

	      ELSE

	        DO I = 1, NSTREAMS
	          SHOM = ZERO
	          DO AA = 1, NSTREAMS
	            H1 = NCON_ALB(AA,N)*XPOS(I,AA,N)*T_DELT_EIGEN(AA,N)
	            H2 = PCON_ALB(AA,N)*XNEG(I,AA,N)
	            SHOM = SHOM + H1 + H2
	          ENDDO
 	          QUAD_ALBEDO_WF(UTA,I,DNIDX) = F1 * SHOM
	        ENDDO

	      ENDIF

	    ENDIF

C  Copy results if quadrature output required (unlikely)

	    IF ( DO_QUAD_OUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        IQD = QUADOUTPUT_INDEX(I)
 	        ALBEDOWF_F(UTA,IQD,DNIDX) = QUAD_ALBEDO_WF(UTA,I,DNIDX)
	      ENDDO
	    ENDIF

C  User-defined stream output, just set to the cumulative source term

	    IF ( DO_USER_STREAMS ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        ALBEDOWF_F(UTA,IUM,DNIDX) = F1 * L_CUMUL_SOURCE(UM)
	      ENDDO
	    ENDIF

	  ENDIF

C  Check for updating the recursion

	  IF ( DO_USER_STREAMS ) THEN
	    IF ( NUT. NE. NUT_PREV ) NSTART = NUT + 1
	    NUT_PREV = NUT
	  ENDIF

C  end loop over optical depth

	ENDDO

C  Finish

	END

C

	SUBROUTINE LIDORT_MIFLUX_ALBEDOWF

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model inputs

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

C  local variables

	INTEGER		 I, IDIR, WDIR, UTA
	DOUBLE PRECISION SUM_MI, SUM_FX

C  mean intensity and flux
C  -----------------------

C  direction loop

	DO IDIR = 1, N_DIRECTIONS

	  WDIR = WHICH_DIRECTIONS(IDIR)

C  loop over all user-defined optical depths

	  DO UTA = 1, N_OUT_USERTAUS

C  integrations

	    SUM_MI = ZERO
	    SUM_FX = ZERO
	    DO I = 1, NSTREAMS
	      SUM_MI = SUM_MI + A(I)  * QUAD_ALBEDO_WF(UTA,I,WDIR)
	      SUM_FX = SUM_FX + AX(I) * QUAD_ALBEDO_WF(UTA,I,WDIR)
	    ENDDO
	    MINT_ALBEDOWF(UTA,WDIR) = SUM_MI * HALF
	    FLUX_ALBEDOWF(UTA,WDIR) = SUM_FX * PI2

C  end loops

	  ENDDO
	ENDDO

C  Finish

	END

C

	SUBROUTINE UPUSER_SURFBBWF
     I    ( DO_INCLUDE_ALBEDO, DO_INCLUDE_MVOUTPUT, R2, F1, FOURIER )

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'
	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Flux factor

	DOUBLE PRECISION F1

C  local control flags

	LOGICAL		 DO_INCLUDE_ALBEDO
	LOGICAL		 DO_INCLUDE_MVOUTPUT

C  albedo factor

	DOUBLE PRECISION R2

C  Fourier number

	INTEGER		 FOURIER

C  local variables
C  ---------------

	INTEGER		 N, NUT, NSTART, NUT_PREV, NLEVEL, NL
	INTEGER		 UTA, UM, IUM, IQD, I, I1, J, UT, AA
	DOUBLE PRECISION L_CUMUL_SOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION L_BOA_SOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION L_LAYER_SOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION NCON_HELP(MAX_USER_STREAMS,MAXSTRM)
	DOUBLE PRECISION PCON_HELP(MAX_USER_STREAMS,MAXSTRM)
	DOUBLE PRECISION HELP(MAXSTRM), REFLEC, L_DOWN
	DOUBLE PRECISION L_FINAL_SOURCE, H1, H2, SHOM
	LOGICAL		 DO_CONSTANTS_SET(MAXLAYER)

C  Zero all Fourier components close to zenith

	IF ( DO_USER_STREAMS ) THEN
	  IF ( LOCAL_UM_START .GT. 1 ) THEN
	    DO UTA = 1, N_OUT_USERTAUS
	      DO UM = 1, LOCAL_UM_START - 1
	        IUM = USEROUTPUT_INDEX(UM)
	        SURFBBWF(UTA,IUM,UPIDX) = ZERO
	      ENDDO
	    ENDDO
	  ENDIF
	ENDIF

C  get the BOA source term
C  -----------------------

	IF ( DO_USER_STREAMS ) THEN

C  initialise boa source function

	  N = NLAYERS
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    L_BOA_SOURCE(UM) = ZERO
	  ENDDO

C  Carry on only if there is an albedo

	  IF ( DO_INCLUDE_ALBEDO ) THEN

C  reflectance from surface
C  .. Downward intensity at computational angles (homog)
C  .. reflectance integrand  a(j).x(j).L_DOWN(-j)

	    DO I = 1, NSTREAMS
	      L_DOWN = ZERO
	      DO AA = 1, NSTREAMS
	        H1 = NCON_SURF(AA,N)*XPOS(I,AA,N)*T_DELT_EIGEN(AA,N)
	        H2 = PCON_SURF(AA,N)*XNEG(I,AA,N) 
	        L_DOWN = L_DOWN + H1 + H2
	      ENDDO
	      HELP(I) = L_DOWN * AX(I)
	    ENDDO

C  reflected multiple scatter intensity at user defined-angles
C  .. integrate reflectance, same for all user-streams in Lambertian case
C  .. integrate with reflectance function at user angles (non-Lambertian)

	    IF ( DO_LAMBERTIAN_ALBEDO ) THEN
	      REFLEC = ZERO
	      DO J = 1, NSTREAMS
	         REFLEC = REFLEC + HELP(J)
	      ENDDO
	      REFLEC = R2 * REFLEC
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        L_BOA_SOURCE(UM) = REFLEC
	      ENDDO
	    ELSE
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        REFLEC = ZERO
	        DO J = 1, NSTREAMS
	          REFLEC = REFLEC + HELP(J)*USER_BIREFLEC(FOURIER,UM,J)
	        ENDDO
	        L_BOA_SOURCE(UM) = R2 * REFLEC
	      ENDDO
	    ENDIF

	  ENDIF
	    
C  Add linearization term due to surf_bb variation 

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    L_BOA_SOURCE(UM) = L_BOA_SOURCE(UM) +
     &               SURFBB * USER_EMISSIVITY(UM)
	  ENDDO

C  Initialize recursion for user-defined stream angles only

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    L_CUMUL_SOURCE(UM) = L_BOA_SOURCE(UM)
	  ENDDO

C  Save multiple scatter term output

	  IF ( SAVE_LAYER_MSST ) THEN
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      SURFBBWF_BOA_STERM(UM) = F1*L_BOA_SOURCE(UM)
	    ENDDO
	  ENDIF

	ENDIF

C  initialise cumulative source term loop

 	NSTART = NLAYERS
	NUT_PREV = NSTART + 1
	DO N = 1, NLAYERS
	  DO_CONSTANTS_SET(N) = .TRUE.
	ENDDO

C  loop over all output optical depths
C  -----------------------------------

	DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

	  NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)

	  IF ( DO_USER_STREAMS ) THEN
	    NUT = NLEVEL + 1
	    DO N = NSTART, NUT, -1

C  Save some calculation time. These quantities are always required.

	      IF ( DO_CONSTANTS_SET(N) ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          DO AA = 1, NSTREAMS
	            NCON_HELP(UM,AA) = NCON_SURF(AA,N) * U_XPOS(UM,AA,N)
	            PCON_HELP(UM,AA) = PCON_SURF(AA,N) * U_XNEG(UM,AA,N)
  	          ENDDO
	        ENDDO
 	      ENDIF

	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
                  H1 = NCON_HELP(UM,AA) * HMULT_2(AA,UM,N)
                  H2 = PCON_HELP(UM,AA) * HMULT_1(AA,UM,N)
	          SHOM = SHOM + H1 + H2
	        ENDDO
	        L_LAYER_SOURCE(UM) = SHOM
	      ENDDO

	      IF ( SAVE_LAYER_MSST ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          SURFBBWF_MSCATSTERM(N,UM,UPIDX) =
     &                L_LAYER_SOURCE(UM) * F1
	        ENDDO
	      ENDIF

	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        L_CUMUL_SOURCE(UM) = L_LAYER_SOURCE(UM)  +
     &                   T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM)
	      ENDDO

	    ENDDO
	  ENDIF

C  Offgrid output
C  --------------

	  IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

	    UT = OFFGRID_UTAU_OUTINDEX(UTA)
	    N  = OFFGRID_UTAU_LAYERIDX(UT)

C  Quadrature output, mean value + flux output

	    IF ( DO_QUAD_OUTPUT. OR. DO_INCLUDE_MVOUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        I1 = I + NSTREAMS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          H1 = NCON_SURF(AA,N)*XPOS(I1,AA,N)*T_UTDN_EIGEN(AA,UT)
	          H2 = PCON_SURF(AA,N)*XNEG(I1,AA,N)*T_UTUP_EIGEN(AA,UT)
	          SHOM = SHOM + H1 + H2
	        ENDDO
	        QUAD_SURFBB_WF(UTA,I,UPIDX) = F1 * SHOM
	      ENDDO
	    ENDIF

C  Copy results if quadrature output required (unlikely)

	    IF ( DO_QUAD_OUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        IQD = QUADOUTPUT_INDEX(I)
 	        SURFBBWF(UTA,IQD,UPIDX) = QUAD_SURFBB_WF(UTA,I,UPIDX)
	      ENDDO
	    ENDIF

C  User-defined stream output

	    IF ( DO_USER_STREAMS ) THEN

C  set help arrays

	      IF ( DO_CONSTANTS_SET(N) ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          DO AA = 1, NSTREAMS
	            NCON_HELP(UM,AA) = NCON_SURF(AA,N) * U_XPOS(UM,AA,N)
	            PCON_HELP(UM,AA) = PCON_SURF(AA,N) * U_XNEG(UM,AA,N)
  	          ENDDO
	        ENDDO
 	      ENDIF
	      DO_CONSTANTS_SET(N) = .FALSE.

C  add additional partial layer source term

	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
                  H1 = NCON_HELP(UM,AA) * UT_HMULT_UP_DN(AA,UM,UT)
                  H2 = PCON_HELP(UM,AA) * UT_HMULT_UP_UP(AA,UM,UT)
	          SHOM = SHOM + H1 + H2
	        ENDDO
	        L_LAYER_SOURCE(UM) = SHOM
	      ENDDO

C  assign final albedo weighting function

	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        L_FINAL_SOURCE = L_LAYER_SOURCE(UM)  +
     &              T_UTUP_USERM(UT,UM)  * L_CUMUL_SOURCE(UM)
	        SURFBBWF(UTA,IUM,UPIDX) = F1 * L_FINAL_SOURCE
	      ENDDO

	    ENDIF

C  Ongrid output
C  -------------

	  ELSE

C  Quadrature output, mean value + flux output

	    IF ( DO_QUAD_OUTPUT. OR. DO_INCLUDE_MVOUTPUT ) THEN

C  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
C  looking at the perturbation field at the top of these layers. The
C  case where the level mask = NLAYERS is the upwelling perturbed fields
C  at the bottom of the atmosphere (treated separately).

	      NL = NLEVEL
	      N = NL + 1

C  For the lowest level

	      IF ( NLEVEL .EQ. NLAYERS ) THEN

	        DO I = 1, NSTREAMS
	          I1 = I + NSTREAMS
	          SHOM = ZERO
	          DO AA = 1, NSTREAMS
	            H1 =  NCON_SURF(AA,NL) *
     &                       XPOS(I1,AA,NL) * T_DELT_EIGEN(AA,NL)
	            H2 = PCON_SURF(AA,NL)*XNEG(I1,AA,NL)
	            SHOM = SHOM + H1 + H2
	          ENDDO
 	          QUAD_SURFBB_WF(UTA,IQD,UPIDX) = F1 * SHOM
	        ENDDO	  

C  For other levels in the atmosphere

	      ELSE

	        DO I = 1, NSTREAMS
	          I1 = I + NSTREAMS
	          SHOM = ZERO
	          DO AA = 1, NSTREAMS
	            H1 = NCON_SURF(AA,N)*XPOS(I1,AA,N)
	            H2 = PCON_SURF(AA,N) *
     &                      XNEG(I1,AA,N) * T_DELT_EIGEN(AA,N)
	            SHOM = SHOM + H1 + H2
	          ENDDO
 	          QUAD_SURFBB_WF(UTA,IQD,UPIDX) = F1 * SHOM
	        ENDDO

	      ENDIF

	    ENDIF

C  Copy results if quadrature output required (unlikely)

	    IF ( DO_QUAD_OUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        IQD = QUADOUTPUT_INDEX(I)
 	        SURFBBWF(UTA,IQD,UPIDX) = QUAD_SURFBB_WF(UTA,I,UPIDX)
	      ENDDO
	    ENDIF

C  User-defined stream output, just set to the cumulative source term

	    IF ( DO_USER_STREAMS ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        SURFBBWF(UTA,IUM,UPIDX) = F1 * L_CUMUL_SOURCE(UM)
	      ENDDO
	    ENDIF

	  ENDIF

C  Check for updating the recursion

	  IF ( DO_USER_STREAMS ) THEN
	    IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
	    NUT_PREV = NUT
	  ENDIF

C  end loop over optical depth

	ENDDO

C  Finish

	END

C

	SUBROUTINE DNUSER_SURFBBWF ( DO_INCLUDE_MVOUTPUT, F1 )

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'
	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Flux factor

	DOUBLE PRECISION F1

C  mean value output flag

	LOGICAL		 DO_INCLUDE_MVOUTPUT

C  local variables
C  ---------------

	INTEGER		 N, NUT, NSTART, NUT_PREV, NLEVEL, NL
	INTEGER		 UTA, UM, IUM, IQD, I, UT, AA
	DOUBLE PRECISION L_CUMUL_SOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION L_TOA_SOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION L_LAYER_SOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION NCON_HELP(MAX_USER_STREAMS,MAXSTRM)
	DOUBLE PRECISION PCON_HELP(MAX_USER_STREAMS,MAXSTRM)
	DOUBLE PRECISION L_FINAL_SOURCE, H1, H2, SHOM
	LOGICAL		 DO_CONSTANTS_SET(MAXLAYER)

C  Zero all Fourier components close to zenith

	IF ( DO_USER_STREAMS ) THEN
	  IF ( LOCAL_UM_START .GT. 1 ) THEN
	    DO UTA = 1, N_OUT_USERTAUS
	      DO UM = 1, LOCAL_UM_START - 1
	        IUM = USEROUTPUT_INDEX(UM)
	        SURFBBWF(UTA,IUM,DNIDX) = ZERO
	      ENDDO
	    ENDDO
	  ENDIF
	ENDIF

C  get the TOA source term
C  -----------------------

	IF ( DO_USER_STREAMS ) THEN

C  initialise boa source function

	  N = NLAYERS
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    L_TOA_SOURCE(UM) = ZERO
	  ENDDO

C  Initialize recursion for user-defined stream angles only

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    L_CUMUL_SOURCE(UM) = L_TOA_SOURCE(UM)
	  ENDDO

	ENDIF

C  initialise cumulative source term loop

 	NSTART = 1
	NUT_PREV = NSTART - 1
	DO N = 1, NLAYERS
	  DO_CONSTANTS_SET(N) = .TRUE.
	ENDDO

C  loop over all output optical depths
C  -----------------------------------

	DO UTA = 1, N_OUT_USERTAUS

C  Layer index for given optical depth

	  NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)

	  IF ( DO_USER_STREAMS ) THEN
	    NUT = NLEVEL
	    DO N = NSTART, NUT

C  Save some calculation time. These quantities are always required.

	      IF ( DO_CONSTANTS_SET(N) ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          DO AA = 1, NSTREAMS
	            NCON_HELP(UM,AA) = NCON_SURF(AA,N) * U_XNEG(UM,AA,N)
	            PCON_HELP(UM,AA) = PCON_SURF(AA,N) * U_XPOS(UM,AA,N)
  	          ENDDO
	        ENDDO
 	      ENDIF

	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
                  H1 = NCON_HELP(UM,AA) * HMULT_1(AA,UM,N)
                  H2 = PCON_HELP(UM,AA) * HMULT_2(AA,UM,N)
	          SHOM = SHOM + H1 + H2
	        ENDDO
	        L_LAYER_SOURCE(UM) = SHOM
	      ENDDO

	      IF ( SAVE_LAYER_MSST ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          SURFBBWF_MSCATSTERM(N,UM,DNIDX) =
     &                L_LAYER_SOURCE(UM) * F1
	        ENDDO
	      ENDIF

	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        L_CUMUL_SOURCE(UM) = L_LAYER_SOURCE(UM)  +
     &                   T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM)
	      ENDDO

	    ENDDO
	  ENDIF

C  Offgrid output
C  --------------

	  IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

	    UT = OFFGRID_UTAU_OUTINDEX(UTA)
	    N  = OFFGRID_UTAU_LAYERIDX(UT)

C  Quadrature output, mean value + flux output

	    IF ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
	          H1 = NCON_SURF(AA,N)*XPOS(I,AA,N)*T_UTDN_EIGEN(AA,UT)
	          H2 = PCON_SURF(AA,N)*XNEG(I,AA,N)*T_UTUP_EIGEN(AA,UT)
	          SHOM = SHOM + H1 + H2
	        ENDDO
	        QUAD_SURFBB_WF(UTA,I,DNIDX) = F1 * SHOM
	      ENDDO
	    ENDIF

C  Copy results if quadrature output required (unlikely)

	    IF ( DO_QUAD_OUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        IQD = QUADOUTPUT_INDEX(I)
 	        SURFBBWF(UTA,IQD,DNIDX) = QUAD_SURFBB_WF(UTA,I,DNIDX)
	      ENDDO
	    ENDIF

C  User-defined stream output

	    IF ( DO_USER_STREAMS ) THEN

C  set help arrays

	      IF ( DO_CONSTANTS_SET(N) ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          DO AA = 1, NSTREAMS
	            NCON_HELP(UM,AA) = NCON_SURF(AA,N) * U_XNEG(UM,AA,N)
	            PCON_HELP(UM,AA) = PCON_SURF(AA,N) * U_XPOS(UM,AA,N)
  	          ENDDO
	        ENDDO
 	      ENDIF
	      DO_CONSTANTS_SET(N) = .FALSE.

C  add additional partial layer source term

	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        SHOM = ZERO
	        DO AA = 1, NSTREAMS
                  H1 = NCON_HELP(UM,AA) * UT_HMULT_DN_DN(AA,UM,UT)
                  H2 = PCON_HELP(UM,AA) * UT_HMULT_DN_UP(AA,UM,UT)
	          SHOM = SHOM + H1 + H2
	        ENDDO
	        L_LAYER_SOURCE(UM) = SHOM
	      ENDDO

C  assign final albedo weighting function

	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        L_FINAL_SOURCE = L_LAYER_SOURCE(UM)  +
     &              T_UTDN_USERM(UT,UM)  * L_CUMUL_SOURCE(UM)
	        SURFBBWF(UTA,IUM,DNIDX) = F1 * L_FINAL_SOURCE
	      ENDDO

	    ENDIF

C  Ongrid output
C  -------------

	  ELSE

C  Quadrature output, mean value + flux output

	    IF ( DO_QUAD_OUTPUT. OR. DO_INCLUDE_MVOUTPUT ) THEN

C  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
C  looking at the perturbation field at the top of these layers. The
C  case where the level mask = NLAYERS is the upwelling perturbed fields
C  at the bottom of the atmosphere (treated separately).

	      NL = NLEVEL
	      N = NL

C  For the highest level

	      IF ( NLEVEL .EQ. 0 ) THEN

	        DO I = 1, NSTREAMS
 	          QUAD_SURFBB_WF(UTA,I,DNIDX) = ZERO
	        ENDDO	  

C  For other levels in the atmosphere

	      ELSE

	        DO I = 1, NSTREAMS
	          SHOM = ZERO
	          DO AA = 1, NSTREAMS
	            H1 = NCON_SURF(AA,N)*XPOS(I,AA,N)*T_DELT_EIGEN(AA,N)
	            H2 = PCON_SURF(AA,N)*XNEG(I,AA,N)
	            SHOM = SHOM + H1 + H2
	          ENDDO
 	          QUAD_SURFBB_WF(UTA,I,DNIDX) = F1 * SHOM
	        ENDDO

	      ENDIF

	    ENDIF

C  Copy results if quadrature output required (unlikely)

	    IF ( DO_QUAD_OUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        IQD = QUADOUTPUT_INDEX(I)
 	        SURFBBWF(UTA,IQD,DNIDX) = QUAD_SURFBB_WF(UTA,I,DNIDX)
	      ENDDO
	    ENDIF

C  User-defined stream output, just set to the cumulative source term

	    IF ( DO_USER_STREAMS ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        SURFBBWF(UTA,IUM,DNIDX) = F1 * L_CUMUL_SOURCE(UM)
	      ENDDO
	    ENDIF

	  ENDIF

C  Check for updating the recursion

	  IF ( DO_USER_STREAMS ) THEN
	    IF ( NUT. NE. NUT_PREV ) NSTART = NUT + 1
	    NUT_PREV = NUT
	  ENDIF

C  end loop over optical depth

	ENDDO

C  Finish

	END

C

	SUBROUTINE LIDORT_MIFLUX_SURFBBWF

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model inputs

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

C  local variables

	INTEGER		 I, IDIR, WDIR, UTA
	DOUBLE PRECISION SUM_MI, SUM_FX

C  mean intensity and flux
C  -----------------------

C  direction loop

	DO IDIR = 1, N_DIRECTIONS

	  WDIR = WHICH_DIRECTIONS(IDIR)

C  loop over all user-defined optical depths

	  DO UTA = 1, N_OUT_USERTAUS

C  integrations

	    SUM_MI = ZERO
	    SUM_FX = ZERO
	    DO I = 1, NSTREAMS
	      SUM_MI = SUM_MI + A(I)  * QUAD_SURFBB_WF(UTA,I,WDIR)
	      SUM_FX = SUM_FX + AX(I) * QUAD_SURFBB_WF(UTA,I,WDIR)
	    ENDDO
	    MINT_SURFBBWF(UTA,WDIR) = SUM_MI * HALF
	    FLUX_SURFBBWF(UTA,WDIR) = SUM_FX * PI2

C  end loops

	  ENDDO
	ENDDO

C  Finish

	END

