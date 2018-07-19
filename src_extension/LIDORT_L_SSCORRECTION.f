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

	SUBROUTINE LIDORT_ATMOSWF_SSCORRECTION
     &      ( SSFLUX, LAYER_TO_VARY, K_PARAMETERS)

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include files of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  Include file of single scatter saved variables

	INCLUDE '../include_s/LIDORT_SSCORRECTION.VARS'

C  Include files of model and control variables (input)

	INCLUDE '../include_e/LIDORT_L_CONTROL.VARS'
	INCLUDE '../include_e/LIDORT_L_GEOPHYS.VARS'

C  include files of linearized setup and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  Include file of results

	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

C  Flux factor

	DOUBLE PRECISION SSFLUX

C  Linearization control

	INTEGER		 LAYER_TO_VARY, K_PARAMETERS

C  local variables
C  ---------------

	INTEGER		 N, NUT, NSTART, NUT_PREV, NLEVEL
	INTEGER		 UT, L, UTA, UM, IUM, IA, NC, Q, NM1, K
	DOUBLE PRECISION HELP, TTT, UVAR, AVAR, FT1, FT2
	DOUBLE PRECISION L_FINAL_SOURCE, L_SS_LAYERSOURCE
	DOUBLE PRECISION MULTIPLIER, L_SSCORRECTION, BEFORE
	DOUBLE PRECISION VAR_TMS(MAX_PARAMETERS), LEGPOLY, L_PHASMOM

	DOUBLE PRECISION L_EXACTSCAT_UP
     S             (MAX_USER_STREAMS,MAX_USER_RELAZMS,MAX_PARAMETERS)
	DOUBLE PRECISION L_EXACTSCAT_DN
     S             (MAX_USER_STREAMS,MAX_USER_RELAZMS,MAX_PARAMETERS)
	DOUBLE PRECISION L_SS_CUMSOURCE
     S             (MAX_USER_STREAMS,MAX_USER_RELAZMS,MAX_PARAMETERS)

C  Set up operations
C  -----------------

	K = LAYER_TO_VARY

C  Create Linearized TMS factor for Layer K
C   ( Use UNSCALED inputs )

	TTT = TMS(K)
	IF ( DO_DELTAM_SCALING ) THEN
	  NM1 = NMOMENTS+1
	  FT1 = TRUNC_FACTOR(K) * TTT
	  FT2 = ONE + FT1
	  DO Q = 1, K_PARAMETERS
	    IF ( DO_PHASFUNC_VARIATION(K,Q) ) THEN
	      UVAR = OMEGA_VARS_TOTAL_INPUT(Q,K)
	      AVAR = UVAR + PHASMOM_VARS_TOTAL_INPUT(Q,NM1,K)
	      VAR_TMS(Q) = UVAR + FT1 * AVAR
	    ELSE
	      UVAR = OMEGA_VARS_TOTAL_INPUT(Q,K)
	      VAR_TMS(Q) = UVAR * FT2
	    ENDIF
	  ENDDO
	ELSE
	  DO Q = 1, K_PARAMETERS
	    UVAR = OMEGA_VARS_TOTAL_INPUT(Q,K)
	    VAR_TMS(Q) = UVAR
	  ENDDO
	ENDIF

C  ####################
C  #    UPWELLING     #
C  ####################

	IF ( DO_UPWELLING ) THEN

C  ===================================
C  Total phase function linearization (upwelling)
C  ===================================

	  IF ( STERM_LAYERMASK_UP(K)) THEN

C  Loop over output viewing directions

	    DO UM = 1, N_USER_STREAMS
	      DO IA = 1, N_USER_RELAZMS

C  Loop over varying parameters for layer K

	        DO Q = 1, K_PARAMETERS

C  Phase function moment variations
C    add TMS correction factor linearization

	          IF ( DO_PHASFUNC_VARIATION(K,Q) ) THEN

	            HELP = ZERO
	            DO L = 0, NMOMENTS_INPUT
	              LEGPOLY = SS_PLEG_UP(UM,IA,L)
	              L_PHASMOM = PHASMOM_VARS_TOTAL_INPUT(Q,L,K) *
     &                           PHASMOMS_TOTAL_INPUT(L,K)
	              HELP = HELP + L_PHASMOM*LEGPOLY
	            ENDDO
	            L_EXACTSCAT_UP(UM,IA,Q) = HELP * TMS(K) +
     &                       EXACTSCAT_UP(UM,IA,K) * VAR_TMS(Q)

	          ELSE

C  No phase function variations - just add TMS linearization

	            L_EXACTSCAT_UP(UM,IA,Q) = 
     &                       EXACTSCAT_UP(UM,IA,K) * VAR_TMS(Q)

	          ENDIF

C  end parameter loop and viewing directions loop

	        ENDDO
	      ENDDO
	    ENDDO

C  Only if layer K exists

	  ENDIF

C  ===================================
C  Upwelling single scatter recurrence
C  ===================================

C  initialize cumulative source term

	  DO UM = 1, N_USER_STREAMS
	    DO IA = 1, N_USER_RELAZMS
	      DO Q = 1, K_PARAMETERS
	        L_SS_CUMSOURCE(UM,IA,Q) = ZERO
	      ENDDO
	    ENDDO
	  ENDDO

C  initialise optical depth loop

 	  NSTART = NLAYERS
	  NUT_PREV = NSTART + 1

C  Main loop over all output optical depths
C  ========================================

	  DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

	    NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

C  finishing layer

	    NUT = NLEVEL + 1

C  Cumulative single scatter source terms to layer NUT
C  ---------------------------------------------------

C    1. Get layer source terms = Exact scattering * Multiplier
C    2. Loop over layers working upwards to NUT

	    DO N = NSTART, NUT, -1
	      NC = NLAYERS + 1 - N

C  If N = K (the layer that is varying)

	      IF ( N .EQ. K ) THEN

	        DO UM = 1, N_USER_STREAMS
	          DO IA = 1, N_USER_RELAZMS
	            DO Q = 1, K_PARAMETERS
	              L_SS_LAYERSOURCE = 
     &                  EXACTSCAT_UP(UM,IA,N) * L_EMULT_UP(UM,N,K,Q)
     &              + L_EXACTSCAT_UP(UM,IA,Q) *   EMULT_UP(UM,N)
	              L_SS_CUMSOURCE(UM,IA,Q) = L_SS_LAYERSOURCE
     &          +   T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(UM,IA,Q)
     &          + L_T_DELT_USERM(N,UM,Q) *   SS_CUMSOURCE_UP(UM,IA,NC-1)
	            ENDDO
	          ENDDO
	        ENDDO

C  Variations when N > K

	      ELSE IF ( N .GT. K ) THEN

	        DO UM = 1, N_USER_STREAMS	         
	          DO IA = 1, N_USER_RELAZMS
	            DO Q = 1, K_PARAMETERS
	              L_SS_LAYERSOURCE = 
     &                  EXACTSCAT_UP(UM,IA,N) * L_EMULT_UP(UM,N,K,Q)
	              L_SS_CUMSOURCE(UM,IA,Q) = L_SS_LAYERSOURCE
     &                  + T_DELT_USERM(N,UM) * L_SS_CUMSOURCE(UM,IA,Q)
	            ENDDO
	          ENDDO
	        ENDDO

C  Transmittance for N < K

	      ELSE

	        DO UM = 1, N_USER_STREAMS	         
	          DO IA = 1, N_USER_RELAZMS
	            DO Q = 1, K_PARAMETERS
	              L_SS_CUMSOURCE(UM,IA,Q) =
     &                    T_DELT_USERM(N,UM) * L_SS_CUMSOURCE(UM,IA,Q)
	            ENDDO
	          ENDDO
	        ENDDO

	      ENDIF

	    ENDDO

C  Offgrid output
C  --------------

C  Set final cumulative source and Correct the Weighting function

	    IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

	      UT = OFFGRID_UTAU_OUTINDEX(UTA)
	      N  = OFFGRID_UTAU_LAYERIDX(UT)

C  If N = K (the layer that is varying)
C    add Linearization of additional partial layer source term =
C        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)

	      IF ( N .EQ. K ) THEN

	        DO UM = 1, N_USER_STREAMS
	          IUM = USEROUTPUT_INDEX(UM)
	          DO IA = 1, N_USER_RELAZMS
	            DO Q = 1, K_PARAMETERS
	              L_SS_LAYERSOURCE = 
     &               EXACTSCAT_UP(UM,IA,N) * L_UT_EMULT_UP(UM,UT,N,Q)
     &           + L_EXACTSCAT_UP(UM,IA,Q) *   UT_EMULT_UP(UM,UT)
	              L_FINAL_SOURCE = L_SS_LAYERSOURCE
     &           + L_T_UTUP_USERM(UT,UM,Q) *   SS_CUMSOURCE_UP(UM,IA,NC)
     &           +   T_UTUP_USERM(UT,UM)   * L_SS_CUMSOURCE(UM,IA,Q)
	              L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
	              BEFORE = ATMOSWF(Q,K,UTA,IUM,UPIDX,IA)
	              ATMOSWF(Q,K,UTA,IUM,UPIDX,IA) =
     &                              BEFORE + L_SSCORRECTION
	            ENDDO
	          ENDDO
	        ENDDO

C  Variations when N > K
C    add Linearization of additional partial layer source term =
C         Exact_Scat * L_Multiplier(k)

	      ELSE IF ( N .GT. K ) THEN

	        DO UM = 1, N_USER_STREAMS	         
	          IUM = USEROUTPUT_INDEX(UM)
	          DO IA = 1, N_USER_RELAZMS
	            DO Q = 1, K_PARAMETERS
	              L_SS_LAYERSOURCE = 
     &               EXACTSCAT_UP(UM,IA,N) * L_UT_EMULT_UP(UM,UT,K,Q)
	              L_FINAL_SOURCE = L_SS_LAYERSOURCE
     &           +   T_UTUP_USERM(UT,UM)   * L_SS_CUMSOURCE(UM,IA,Q)
	              L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
	              BEFORE = ATMOSWF(Q,K,UTA,IUM,UPIDX,IA)
	              ATMOSWF(Q,K,UTA,IUM,UPIDX,IA) =
     &                              BEFORE + L_SSCORRECTION
	            ENDDO
	          ENDDO
	        ENDDO

C  Transmittance for the other layers

	      ELSE

	        DO UM = 1, N_USER_STREAMS	         
	          IUM = USEROUTPUT_INDEX(UM)
	          DO IA = 1, N_USER_RELAZMS
	            DO Q = 1, K_PARAMETERS
	              L_FINAL_SOURCE =
     &                  T_UTUP_USERM(UT,UM) * L_SS_CUMSOURCE(UM,IA,Q)
	              L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
	              BEFORE = ATMOSWF(Q,K,UTA,IUM,UPIDX,IA)
	              ATMOSWF(Q,K,UTA,IUM,UPIDX,IA) =
     &                              BEFORE + L_SSCORRECTION
	            ENDDO
	          ENDDO
	        ENDDO

	      ENDIF

C  Ongrid output
C  -------------

C  just set to the cumulative source term and Correct the Weighting Function

	    ELSE

	      DO UM = 1, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        DO IA = 1, N_USER_RELAZMS
	          DO Q = 1, K_PARAMETERS
	            L_FINAL_SOURCE = L_SS_CUMSOURCE(UM,IA,Q)
	            L_SSCORRECTION = SSFLUX * L_SS_CUMSOURCE(UM,IA,Q)
	            BEFORE = ATMOSWF(Q,K,UTA,IUM,UPIDX,IA)
	            ATMOSWF(Q,K,UTA,IUM,UPIDX,IA) =
     &                              BEFORE + L_SSCORRECTION
	          ENDDO
	        ENDDO
	      ENDDO

	    ENDIF

C  Check for updating the recursion

	    IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
	    NUT_PREV = NUT

C  end loop over optical depth

	  ENDDO

C  end Upwelling clause

	ENDIF

C  ####################
C  #   DOWNWELLING    #
C  ####################

	IF ( DO_DNWELLING ) THEN

C  ===================================
C  Total phase function linearization (upwelling)
C  ===================================

	  IF ( STERM_LAYERMASK_DN(K)) THEN

C  Loop over output viewing directions

	    DO UM = 1, N_USER_STREAMS
	      DO IA = 1, N_USER_RELAZMS

C  Loop over varying parameters for layer K

	        DO Q = 1, K_PARAMETERS

C  Phase function moment variations
C    add TMS correction factor linearization

	          IF ( DO_PHASFUNC_VARIATION(K,Q) ) THEN

	            HELP = ZERO
	            DO L = 0, NMOMENTS_INPUT
	              LEGPOLY = SS_PLEG_DN(UM,IA,L)
	              L_PHASMOM = PHASMOM_VARS_TOTAL_INPUT(Q,L,K) *
     &                             PHASMOMS_TOTAL_INPUT(L,K)
	              HELP = HELP + L_PHASMOM*LEGPOLY
	            ENDDO
	            L_EXACTSCAT_DN(UM,IA,Q) = HELP * TMS(K) +
     &                       EXACTSCAT_DN(UM,IA,K) * VAR_TMS(Q)

	          ELSE

C  No phase function variations - just add TMS linearization

	            L_EXACTSCAT_DN(UM,IA,Q) = 
     &                       EXACTSCAT_DN(UM,IA,K) * VAR_TMS(Q)

	          ENDIF

C  end parameter loop and viewing directions loop

	        ENDDO
	      ENDDO
	    ENDDO

C  Only if layer K exists

	  ENDIF

C  =====================================
C  Downwelling single scatter recurrence
C  =====================================

C  initialize cumulative source term

	  DO UM = 1, N_USER_STREAMS
	    DO IA = 1, N_USER_RELAZMS
	      DO Q = 1, K_PARAMETERS
	        L_SS_CUMSOURCE(UM,IA,Q) = ZERO
	      ENDDO
	    ENDDO
	  ENDDO

C  initialise optical depth loop

 	  NSTART = 1
	  NUT_PREV = NSTART - 1

C  Main loop over all output optical depths
C  ========================================

	  DO UTA = 1, N_OUT_USERTAUS

C  Layer index for given optical depth

	    NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

C  finishing layer

	    NUT = NLEVEL

C  Cumulative single scatter source terms to layer NUT
C  ---------------------------------------------------

C    1. Get layer source terms = Exact scattering * Multiplier
C    2. Loop over layers working upwards to NUT

	    DO N = NSTART, NUT
	      NC = N

C  If N = K (the layer that is varying)

	      IF ( N .EQ. K ) THEN

	        DO UM = 1, N_USER_STREAMS
	          DO IA = 1, N_USER_RELAZMS
	            DO Q = 1, K_PARAMETERS
	              L_SS_LAYERSOURCE = 
     &                  EXACTSCAT_DN(UM,IA,N) * L_EMULT_DN(UM,N,K,Q)
     &              + L_EXACTSCAT_DN(UM,IA,Q) *   EMULT_DN(UM,N)
	              L_SS_CUMSOURCE(UM,IA,Q) = L_SS_LAYERSOURCE
     &          +   T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(UM,IA,Q)
     &          + L_T_DELT_USERM(N,UM,Q) *   SS_CUMSOURCE_DN(UM,IA,NC-1)
	            ENDDO
	          ENDDO
	        ENDDO

C  Variations when N > K

	      ELSE IF ( N .GT. K ) THEN

	        DO UM = 1, N_USER_STREAMS	         
	          DO IA = 1, N_USER_RELAZMS
	            DO Q = 1, K_PARAMETERS
	              L_SS_LAYERSOURCE = 
     &                  EXACTSCAT_DN(UM,IA,N) * L_EMULT_DN(UM,N,K,Q)
	              L_SS_CUMSOURCE(UM,IA,Q) = L_SS_LAYERSOURCE
     &                  + T_DELT_USERM(N,UM) * L_SS_CUMSOURCE(UM,IA,Q)
	            ENDDO
	          ENDDO
	        ENDDO

C  Otherwise transmittance

	      ELSE

	        DO UM = 1, N_USER_STREAMS	         
	          DO IA = 1, N_USER_RELAZMS
	            DO Q = 1, K_PARAMETERS
	              L_SS_CUMSOURCE(UM,IA,Q) =
     &                  T_DELT_USERM(N,UM) * L_SS_CUMSOURCE(UM,IA,Q)
	            ENDDO
	          ENDDO
	        ENDDO

	      ENDIF

	    ENDDO

C  Offgrid output
C  --------------

C  add additional partial layer source term = Exact Scat * Multiplier
C  Set final cumulative source and Correct the intensity

	    IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

	      UT = OFFGRID_UTAU_OUTINDEX(UTA)
	      N  = OFFGRID_UTAU_LAYERIDX(UT)

C  If N = K (the layer that is varying)
C    add Linearization of additional partial layer source term =
C        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)

	      IF ( N .EQ. K ) THEN

	        DO UM = 1, N_USER_STREAMS
	          IUM = USEROUTPUT_INDEX(UM)
	          DO IA = 1, N_USER_RELAZMS
	            DO Q = 1, K_PARAMETERS
	              L_SS_LAYERSOURCE = 
     &               EXACTSCAT_DN(UM,IA,N) * L_UT_EMULT_DN(UM,UT,N,Q)
     &           + L_EXACTSCAT_DN(UM,IA,Q) *   UT_EMULT_DN(UM,UT)
	              L_FINAL_SOURCE = L_SS_LAYERSOURCE
     &           + L_T_UTDN_USERM(UT,UM,Q) *   SS_CUMSOURCE_DN(UM,IA,NC)
     &           +   T_UTDN_USERM(UT,UM)   * L_SS_CUMSOURCE(UM,IA,Q)
	              L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
	              BEFORE = ATMOSWF(Q,K,UTA,IUM,DNIDX,IA)
                      ATMOSWF(Q,K,UTA,IUM,DNIDX,IA) =
     &                        BEFORE + L_SSCORRECTION
	            ENDDO
	          ENDDO
	        ENDDO

C  Variations when N > K
C    add Linearization of additional partial layer source term =
C         Exact_Scat * L_Multiplier(k)

	      ELSE IF ( N .GT. K ) THEN

	        DO UM = 1, N_USER_STREAMS	         
	          IUM = USEROUTPUT_INDEX(UM)
	          DO IA = 1, N_USER_RELAZMS
	            DO Q = 1, K_PARAMETERS
	              L_SS_LAYERSOURCE = 
     &               EXACTSCAT_DN(UM,IA,N) * L_UT_EMULT_DN(UM,UT,K,Q)
	              L_FINAL_SOURCE = L_SS_LAYERSOURCE
     &           +   T_UTDN_USERM(UT,UM)   * L_SS_CUMSOURCE(UM,IA,Q)
	              L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
	              BEFORE = ATMOSWF(Q,K,UTA,IUM,DNIDX,IA)
                      ATMOSWF(Q,K,UTA,IUM,DNIDX,IA) =
     &                        BEFORE + L_SSCORRECTION
	            ENDDO
	          ENDDO
	        ENDDO

C  Otherwise transmittance

	      ELSE

	        DO UM = 1, N_USER_STREAMS	         
	          IUM = USEROUTPUT_INDEX(UM)
	          DO IA = 1, N_USER_RELAZMS
	            DO Q = 1, K_PARAMETERS
	              L_FINAL_SOURCE =
     &                 T_UTDN_USERM(UT,UM)   * L_SS_CUMSOURCE(UM,IA,Q)
	              L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
	              BEFORE = ATMOSWF(Q,K,UTA,IUM,DNIDX,IA)
                      ATMOSWF(Q,K,UTA,IUM,DNIDX,IA) =
     &                        BEFORE + L_SSCORRECTION
	            ENDDO
	          ENDDO
	        ENDDO

	      ENDIF

C  Ongrid output
C  -------------

C  just set to the cumulative source term and Correct the Weighting Function

	    ELSE

	      DO UM = 1, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        DO IA = 1, N_USER_RELAZMS
	          DO Q = 1, K_PARAMETERS
	            L_FINAL_SOURCE = L_SS_CUMSOURCE(UM,IA,Q)
	            L_SSCORRECTION = SSFLUX * L_SS_CUMSOURCE(UM,IA,Q)
	            BEFORE = ATMOSWF(Q,K,UTA,IUM,DNIDX,IA)
                    ATMOSWF(Q,K,UTA,IUM,DNIDX,IA) =
     &                        BEFORE + L_SSCORRECTION
	          ENDDO
	        ENDDO
	      ENDDO

	    ENDIF

C  Check for updating the recursion

	    IF ( NUT. NE. NUT_PREV ) NSTART = NUT + 1
	    NUT_PREV = NUT

C  end loop over optical depth

	  ENDDO

C  end Downwelling clause

	ENDIF

C  Finish

	END
