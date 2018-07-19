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

	SUBROUTINE L_WHOLELAYER_HMULT_UP ( N, N_PARAMETERS )

C  Include files
C  =============

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model variables (input)

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and number of parameters

	INTEGER		 N
	INTEGER		 N_PARAMETERS

C  Local variables
C  ---------------

	INTEGER		 UM, AA, Q
        DOUBLE PRECISION UDEL, ZDEL, L_UDEL, L_ZDEL
	DOUBLE PRECISION L_T2, L_T1, HOM1, HOM2, SM

C  existence flag

	L_HMULT_EXIST(N) = .TRUE.

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  UDEL = T_DELT_USERM(N,UM)
	  SM = USER_SECANTS(UM)
	  DO AA = 1, NSTREAMS
	    ZDEL  = T_DELT_EIGEN(AA,N)
	    DO Q = 1, N_PARAMETERS
	      L_ZDEL = L_T_DELT_EIGEN(AA,N,Q)
	      L_UDEL = L_T_DELT_USERM(N,UM,Q)
	      L_T2 = - ZDEL * L_UDEL - L_ZDEL * UDEL
	      L_T1 = L_ZDEL - L_UDEL
	      HOM1 =   L_KEIGEN(AA,N,Q) * HMULT_1(AA,UM,N) + SM * L_T1
	      HOM2 = - L_KEIGEN(AA,N,Q) * HMULT_2(AA,UM,N) + SM * L_T2
	      L_HMULT_1(AA,UM,N,Q) = ZETA_M(AA,UM,N) * HOM1
	      L_HMULT_2(AA,UM,N,Q) = ZETA_P(AA,UM,N) * HOM2
	    ENDDO
	  ENDDO
	ENDDO

C  Finish

	END

C

	SUBROUTINE L_WHOLELAYER_HMULT_DN ( N, N_PARAMETERS )

C  Include files
C  =============

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model variables (input)

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and number of parameters

	INTEGER		 N
	INTEGER		 N_PARAMETERS

C  Local variables
C  ---------------

	INTEGER		 UM, AA, Q
        DOUBLE PRECISION UDEL, ZDEL, L_UDEL, L_ZDEL
	DOUBLE PRECISION L_T2, L_T1, HOM1, HOM2, SM

C  existence flag, return if set

	IF ( L_HMULT_EXIST(N) ) RETURN

C  Get linearized multipliers

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  UDEL = T_DELT_USERM(N,UM)
	  SM = USER_SECANTS(UM)
	  DO AA = 1, NSTREAMS
	    ZDEL  = T_DELT_EIGEN(AA,N)
	    DO Q = 1, N_PARAMETERS
	      L_ZDEL = L_T_DELT_EIGEN(AA,N,Q)
	      L_UDEL = L_T_DELT_USERM(N,UM,Q)
	      L_T2 = - ZDEL * L_UDEL - L_ZDEL * UDEL
	      L_T1 = L_ZDEL - L_UDEL
	      HOM1 =   L_KEIGEN(AA,N,Q) * HMULT_1(AA,UM,N) + SM * L_T1
	      HOM2 = - L_KEIGEN(AA,N,Q) * HMULT_2(AA,UM,N) + SM * L_T2
	      L_HMULT_1(AA,UM,N,Q) = ZETA_M(AA,UM,N) * HOM1
	      L_HMULT_2(AA,UM,N,Q) = ZETA_P(AA,UM,N) * HOM2
	    ENDDO
	  ENDDO
	ENDDO

C  Finish

	END

C

	SUBROUTINE L_PARTLAYER_HMULT_UP  ( N, UT, N_PARAMETERS )

C  Include files
C  =============

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model variables (input)

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer indices and number of parameters

	INTEGER		 N, UT
	INTEGER		 N_PARAMETERS

C  Local variables
C  ---------------

	INTEGER		 UM, AA, Q
        DOUBLE PRECISION UX_UP, ZDEL, L_UX_UP, L_ZDEL
	DOUBLE PRECISION L_T2, L_T1, H1, H2, SM, FA

C  Get the multipliers

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  UX_UP = T_UTUP_USERM(UT,UM)
	  SM = USER_SECANTS(UM)
	  DO AA = 1, NSTREAMS
	    ZDEL  = T_DELT_EIGEN(AA,N)
	    DO Q = 1, N_PARAMETERS
	      FA = L_KEIGEN(AA,N,Q)
	      L_ZDEL  = L_T_DELT_EIGEN(AA,N,Q)
	      L_UX_UP = L_T_UTUP_USERM(UT,UM,Q)
	      L_T1 = - ZDEL * L_UX_UP - L_ZDEL * UX_UP
	      L_T1 = L_T_UTDN_EIGEN(AA,UT,Q) + L_T1
	      L_T2 = L_T_UTUP_EIGEN(AA,UT,Q) - L_UX_UP
	      H1 = - FA * UT_HMULT_UP_DN(AA,UM,UT) + SM*L_T1
	      H2 =   FA * UT_HMULT_UP_UP(AA,UM,UT) + SM*L_T2
	      L_UT_HMULT_UP_UP(AA,UM,UT,Q) = ZETA_M(AA,UM,N) * H2
	      L_UT_HMULT_UP_DN(AA,UM,UT,Q) = ZETA_P(AA,UM,N) * H1
	    ENDDO
	  ENDDO
	ENDDO

C  Finish

	END

C

	SUBROUTINE L_PARTLAYER_HMULT_DN ( N, UT, N_PARAMETERS )

C  Include files
C  =============

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model variables (input)

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer indices and number of parameters

	INTEGER		 N, UT
	INTEGER		 N_PARAMETERS

C  Local variables
C  ---------------

	INTEGER		 UM, AA, Q
        DOUBLE PRECISION UX_DN, ZDEL, L_UX_DN, L_ZDEL
	DOUBLE PRECISION L_T2, L_T1, H1, H2, SM, FA

C  Get the linearized multipliers

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  UX_DN = T_UTDN_USERM(UT,UM)
	  SM = USER_SECANTS(UM)
	  DO AA = 1, NSTREAMS
	    ZDEL  = T_DELT_EIGEN(AA,N)
	    DO Q = 1, N_PARAMETERS
	      FA = L_KEIGEN(AA,N,Q)
	      L_ZDEL  = L_T_DELT_EIGEN(AA,N,Q)
	      L_UX_DN = L_T_UTDN_USERM(UT,UM,Q)
	      L_T2 = - ZDEL * L_UX_DN - L_ZDEL * UX_DN
	      L_T2 = L_T_UTUP_EIGEN(AA,UT,Q) + L_T2
	      L_T1 = L_T_UTDN_EIGEN(AA,UT,Q) - L_UX_DN
	      H1 =   FA * UT_HMULT_DN_DN(AA,UM,UT) + SM*L_T1
	      H2 = - FA * UT_HMULT_DN_UP(AA,UM,UT) + SM*L_T2
	      L_UT_HMULT_DN_UP(AA,UM,UT,Q) = ZETA_P(AA,UM,N) * H2
	      L_UT_HMULT_DN_DN(AA,UM,UT,Q) = ZETA_M(AA,UM,N) * H1
	    ENDDO
	  ENDDO
	ENDDO

C  Finish

	END

C

	SUBROUTINE L_WHOLELAYER_EMULT_UP ( N, K, K_PARAMETERS )

C  Include files
C  =============

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index, varying layer index and number of parameters

	INTEGER		 N, K, K_PARAMETERS

C  local variables
C  ---------------

	DOUBLE PRECISION SU, V1, V2, WDEL, UDEL
	INTEGER		 UM, Q

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO UM = 1, N_USER_STREAMS
	    DO Q = 1, K_PARAMETERS
	      L_EMULT_UP(UM,N,K,Q) = ZERO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  transmittance factor

	WDEL = T_DELT_MUBAR(N)

C  For the quasi-spherical case
C  ----------------------------

	IF ( DO_QUASPHER_BEAM ) THEN

C  ..(a) If N = K, multiplier for due to variations in the layer N

	  IF ( K.EQ.N ) THEN

	    DO UM = 1, N_USER_STREAMS
	      UDEL = T_DELT_USERM(N,UM)
	      SU = - ITRANS_USERM(N,UM) / SIGMA_P(N,UM)
	      DO Q = 1, K_PARAMETERS
	        V1 = -L_AVERAGE_SECANT(N,K,Q) / SIGMA_P(N,UM)
	        V2 = WDEL * L_T_DELT_USERM(N,UM,Q) +
     *               UDEL * L_T_DELT_MUBAR(N,K,Q)
	        L_EMULT_UP(UM,N,K,Q) = EMULT_UP(UM,N) * V1 + SU * V2
	      ENDDO
	    ENDDO

C  ..(b) N > K, multiplier due to variations in a higher layer K

	  ELSE IF ( N.GT.K ) THEN

	    DO UM = 1, N_USER_STREAMS
	      UDEL = T_DELT_USERM(N,UM)
	      SU = - ITRANS_USERM(N,UM) / SIGMA_P(N,UM)
	      DO Q = 1, K_PARAMETERS
	        V1 = L_INITIAL_TRANS(N,K,Q) -
     &             ( L_AVERAGE_SECANT(N,K,Q) / SIGMA_P(N,UM) )
	        V2 =  UDEL * L_T_DELT_MUBAR(N,K,Q)
	        L_EMULT_UP(UM,N,K,Q) = EMULT_UP(UM,N) * V1 + SU * V2
	      ENDDO
	    ENDDO

	  ENDIF

C  For the plane-parallel case
C  ---------------------------

	ELSE

C  ..(a) If N = K, multiplier for due to variations in the layer N

	  IF ( K.EQ.N ) THEN

	    DO UM = 1, N_USER_STREAMS
	      UDEL = T_DELT_USERM(N,UM)
	      SU = - ITRANS_USERM(N,UM) / SIGMA_P(N,UM)
	      DO Q = 1, K_PARAMETERS
	        V2 = WDEL * L_T_DELT_USERM(N,UM,Q) +
     &               UDEL * L_T_DELT_MUBAR(N,K,Q)
	        L_EMULT_UP(UM,N,K,Q) =  SU * V2
	      ENDDO
	    ENDDO

C  ..(b) N > K, multiplier due to variations in a higher layer K

	  ELSE IF ( N.GT.K ) THEN

	    DO UM = 1, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        V1 = L_INITIAL_TRANS(N,K,Q)
	        L_EMULT_UP(UM,N,K,Q) = EMULT_UP(UM,N) * V1
	      ENDDO
	    ENDDO

	  ENDIF

C  End clause Qs versus PP

	ENDIF

C  Finish

	END

C

	SUBROUTINE L_WHOLELAYER_EMULT_DN ( N, K, K_PARAMETERS )

C  Include files
C  =============

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index, varying layer index and number of parameters

	INTEGER		 N, K, K_PARAMETERS

C  local variables
C  ---------------

	DOUBLE PRECISION SD, V1, V2
	INTEGER		 UM, Q

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO UM = 1, N_USER_STREAMS
	    DO Q = 1, K_PARAMETERS
	      L_EMULT_DN(UM,N,K,Q) = ZERO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  NOTE - use of L'Hopital's Rule is present in this module

C  For the quasi-spherical case
C  ----------------------------

	IF ( DO_QUASPHER_BEAM ) THEN

C  ..(a) If N = K, multiplier for due to variations in the layer N

	  IF ( K.EQ.N ) THEN

	    DO UM = 1, N_USER_STREAMS
	      IF ( EMULT_HOPRULE(N,UM) ) THEN
	        V1 = ONE - DELTA(N) * USER_SECANTS(UM)
	        V2 = - HALF * DELTA(N)
	        DO Q = 1, K_PARAMETERS
	          SD = V1 * VQ(N,Q) + V2 * L_AVERAGE_SECANT(N,K,Q)
	          L_EMULT_DN(UM,N,K,Q) = EMULT_DN(UM,N) * SD
	        ENDDO
	      ELSE
	        SD = ITRANS_USERM(N,UM) / SIGMA_M(N,UM)
	        DO Q = 1, K_PARAMETERS
	          V1 = - L_AVERAGE_SECANT(N,K,Q) / SIGMA_M(N,UM)
	          V2 = L_T_DELT_USERM(N,UM,Q) - L_T_DELT_MUBAR(N,K,Q)
	          L_EMULT_DN(UM,N,K,Q) = EMULT_DN(UM,N) * V1 + SD * V2
	        ENDDO
	      ENDIF
	    ENDDO

C  ..(b) N > K, multiplier due to variations in a higher layer K

	  ELSE IF ( N.GT.K ) THEN

	    DO UM = 1, N_USER_STREAMS
	      IF ( EMULT_HOPRULE(N,UM) ) THEN
	        V2 = - HALF * DELTA(N)
	        DO Q = 1, K_PARAMETERS
	          SD =      L_INITIAL_TRANS(N,K,Q) +
     &                 V2 * L_AVERAGE_SECANT(N,K,Q)
	          L_EMULT_DN(UM,N,K,Q) = EMULT_DN(UM,N) * SD
	        ENDDO
	      ELSE
	        SD = ITRANS_USERM(N,UM) / SIGMA_M(N,UM)
	        DO Q = 1, K_PARAMETERS
	          V1 = L_INITIAL_TRANS(N,K,Q) -
     &               ( L_AVERAGE_SECANT(N,K,Q) / SIGMA_M(N,UM) )
	          V2 = - L_T_DELT_MUBAR(N,K,Q)
	          L_EMULT_DN(UM,N,K,Q) = EMULT_DN(UM,N) * V1 + SD * V2
	        ENDDO
	      ENDIF
	    ENDDO

	  ENDIF

C  For the plane-parallel case
C  ---------------------------

	ELSE

C  ..(a) If N = K, multiplier for due to variations in the layer N

	  IF ( K.EQ.N ) THEN

	    DO UM = 1, N_USER_STREAMS
	      IF ( EMULT_HOPRULE(N,UM) ) THEN
	        V1 = ONE - DELTA(N) * USER_SECANTS(UM)
	        DO Q = 1, K_PARAMETERS
	          SD = V1 * VQ(N,Q)
	          L_EMULT_DN(UM,N,K,Q) = EMULT_DN(UM,N) * SD
	        ENDDO
	      ELSE
	        SD = ITRANS_USERM(N,UM) / SIGMA_M(N,UM)
	        DO Q = 1, K_PARAMETERS
	          V2 = L_T_DELT_USERM(N,UM,Q) - L_T_DELT_MUBAR(N,K,Q)
	          L_EMULT_DN(UM,N,K,Q) =  SD * V2
	        ENDDO
	      ENDIF
	    ENDDO

C  ..(b) N > K, multiplier due to variations in a higher layer K

	  ELSE IF ( N.GT.K ) THEN

	    DO UM = 1, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        V1 = L_INITIAL_TRANS(N,K,Q)
	        L_EMULT_DN(UM,N,K,Q) = EMULT_DN(UM,N) * V1
	      ENDDO
	    ENDDO

	  ENDIF

C  End clause Qs versus PP

	ENDIF

C  Finish

	END

C

	SUBROUTINE L_PARTLAYER_EMULT_UP ( N, UT, K, K_PARAMETERS )

C  Include files
C  =============

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given offgrid indices, varying layer index and number of parameters

	INTEGER		 N, UT, K, K_PARAMETERS

C  local variables
C  ---------------

	DOUBLE PRECISION SU, V1, V2, WDEL, UX_UP
	INTEGER		 UM, Q

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO UM = 1, N_USER_STREAMS
	    DO Q = 1, K_PARAMETERS
	      L_UT_EMULT_UP(UM,UT,K,Q) = ZERO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  transmittance factor

	WDEL = T_DELT_MUBAR(N)

C  For the quasi-spherical case
C  ----------------------------

	IF ( DO_QUASPHER_BEAM ) THEN

C  ..(a) If N = K, multiplier for due to variations in the layer N

	  IF ( K.EQ.N ) THEN

	    DO UM = 1, N_USER_STREAMS
	      UX_UP = T_UTUP_USERM(UT,UM)
	      SU = ITRANS_USERM(N,UM) / SIGMA_P(N,UM)
	      DO Q = 1, K_PARAMETERS
	        V1 = -L_AVERAGE_SECANT(N,K,Q) / SIGMA_P(N,UM)
	        V2 =         L_T_UTDN_MUBAR(UT,K,Q) -
     *               UX_UP * L_T_DELT_MUBAR(N,K,Q) -
     *               WDEL  * L_T_UTUP_USERM(UT,UM,Q)
	        L_UT_EMULT_UP(UM,UT,K,Q) =       SU * V2 +
     &                           UT_EMULT_UP(UM,UT) * V1
	      ENDDO
	    ENDDO

C  ..(b) N > K, multiplier due to variations in a higher layer K

	  ELSE IF ( N.GT.K ) THEN

	    DO UM = 1, N_USER_STREAMS
	      UX_UP = T_UTUP_USERM(UT,UM)
	      SU = ITRANS_USERM(N,UM) / SIGMA_P(N,UM)
	      DO Q = 1, K_PARAMETERS
	        V1 = L_INITIAL_TRANS(N,K,Q) -
     &             ( L_AVERAGE_SECANT(N,K,Q) / SIGMA_P(N,UM) )
	        V2 =         L_T_UTDN_MUBAR(UT,K,Q) -
     *               UX_UP * L_T_DELT_MUBAR(N,K,Q)
	        L_UT_EMULT_UP(UM,UT,K,Q) =       SU * V2 +
     &                           UT_EMULT_UP(UM,UT) * V1
	      ENDDO
	    ENDDO

	  ENDIF

C  For the plane-parallel case
C  ---------------------------

	ELSE

C  ..(a) If N = K, multiplier for due to variations in the layer N

	  IF ( K.EQ.N ) THEN

	    DO UM = 1, N_USER_STREAMS
	      UX_UP = T_UTUP_USERM(UT,UM)
	      SU = ITRANS_USERM(N,UM) / SIGMA_P(N,UM)
	      DO Q = 1, K_PARAMETERS
	        V2 =         L_T_UTDN_MUBAR(UT,K,Q) -
     *               UX_UP * L_T_DELT_MUBAR(N,K,Q) -
     *               WDEL  * L_T_UTUP_USERM(UT,UM,Q)
	        L_UT_EMULT_UP(UM,UT,K,Q) =       SU * V2
	      ENDDO
	    ENDDO

C  ..(b) N > K, multiplier due to variations in a higher layer K

	  ELSE IF ( N.GT.K ) THEN

	    DO UM = 1, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        V1 = L_INITIAL_TRANS(N,K,Q)
	        L_UT_EMULT_UP(UM,UT,K,Q) = UT_EMULT_UP(UM,UT) * V1
	      ENDDO
	    ENDDO

	  ENDIF

C  End clause Qs versus PP

	ENDIF

C  Finish

	END

C

	SUBROUTINE L_PARTLAYER_EMULT_DN ( N, UT, K, K_PARAMETERS )

C  Include files
C  =============

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given offgrid indices, varying layer index and number of parameters

	INTEGER		 N, UT, K, K_PARAMETERS

C  local variables
C  ---------------

	DOUBLE PRECISION SD, V1, V2
	INTEGER		 UM, Q

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO UM = 1, N_USER_STREAMS
	    DO Q = 1, K_PARAMETERS
	      L_UT_EMULT_DN(UM,UT,K,Q) = ZERO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  NOTE - use of L'Hopital's Rule is present in this module

C  For the quasi-spherical case
C  ----------------------------

	IF ( DO_QUASPHER_BEAM ) THEN

C  ..(a) If N = K, multiplier for due to variations in the layer N

	  IF ( K.EQ.N ) THEN

	    DO UM = 1, N_USER_STREAMS
	      IF ( EMULT_HOPRULE(N,UM) ) THEN
	        V1 = ONE - OFFGRID_UTAU_VALUES(UT) * USER_SECANTS(UM)
	        V2 = - HALF * OFFGRID_UTAU_VALUES(UT)
	        DO Q = 1, K_PARAMETERS
	          SD = V1 * VQ(N,Q) + V2 * L_AVERAGE_SECANT(N,K,Q)
	          L_UT_EMULT_DN(UM,UT,K,Q) = UT_EMULT_DN(UM,UT) * SD
	        ENDDO
	      ELSE
	        SD = ITRANS_USERM(N,UM) / SIGMA_M(N,UM)
	        DO Q = 1, K_PARAMETERS
	          V1 = - L_AVERAGE_SECANT(N,K,Q) / SIGMA_M(N,UM)
	          V2 = L_T_UTDN_USERM(UT,UM,Q) - L_T_UTDN_MUBAR(UT,K,Q)
	          L_UT_EMULT_DN(UM,UT,K,Q) =       SD * V2 +
     &                           UT_EMULT_DN(UM,UT) * V1
	        ENDDO
	      ENDIF
	    ENDDO

C  ..(b) N > K, multiplier due to variations in a higher layer K

	  ELSE IF ( N.GT.K ) THEN

	    DO UM = 1, N_USER_STREAMS
	      IF ( EMULT_HOPRULE(N,UM) ) THEN
	        V2 = - HALF * OFFGRID_UTAU_VALUES(UT)
	        DO Q = 1, K_PARAMETERS
	          SD =      L_INITIAL_TRANS(N,K,Q) +
     &                 V2 * L_AVERAGE_SECANT(N,K,Q)
	          L_UT_EMULT_DN(UM,UT,K,Q) = UT_EMULT_DN(UM,UT) * SD
	        ENDDO
	      ELSE
	        SD = ITRANS_USERM(N,UM) / SIGMA_M(N,UM)
	        DO Q = 1, K_PARAMETERS
	          V1 = L_INITIAL_TRANS(N,K,Q) -
     &               ( L_AVERAGE_SECANT(N,K,Q) / SIGMA_M(N,UM) )
	          V2 = - L_T_UTDN_MUBAR(UT,K,Q)
	          L_UT_EMULT_DN(UM,UT,K,Q) =       SD * V2 +
     &                           UT_EMULT_DN(UM,UT) * V1
	        ENDDO
	      ENDIF
	    ENDDO

	  ENDIF

C  For the plane-parallel case
C  ---------------------------

	ELSE

C  ..(a) If N = K, multiplier for due to variations in the layer N

	  IF ( K.EQ.N ) THEN

	    DO UM = 1, N_USER_STREAMS
	      IF ( EMULT_HOPRULE(N,UM) ) THEN
	        V1 = ONE - OFFGRID_UTAU_VALUES(UT) * USER_SECANTS(UM)
	        DO Q = 1, K_PARAMETERS
	          SD = V1 * VQ(N,Q)
	          L_UT_EMULT_DN(UM,UT,K,Q) = UT_EMULT_DN(UM,UT) * SD
	        ENDDO
	      ELSE
	        SD = ITRANS_USERM(N,UM) / SIGMA_M(N,UM)
	        DO Q = 1, K_PARAMETERS
	          V2 = L_T_UTDN_USERM(UT,UM,Q) - L_T_UTDN_MUBAR(UT,K,Q)
	          L_UT_EMULT_DN(UM,UT,K,Q) =  SD * V2 
	        ENDDO
	      ENDIF
	    ENDDO

C  ..(b) N > K, multiplier due to variations in a higher layer K

	  ELSE IF ( N.GT.K ) THEN

	    DO UM = 1, N_USER_STREAMS
	      DO Q = 1, K_PARAMETERS
	        V1 = L_INITIAL_TRANS(N,K,Q)
	        L_UT_EMULT_DN(UM,UT,K,Q) = UT_EMULT_DN(UM,UT) * V1
	      ENDDO
	    ENDDO

	  ENDIF

C  End clause Qs versus PP

	ENDIF

C  Finish

	END

C

	SUBROUTINE L_QUAD_GFUNCMULT
     &         ( N, UT, K, K_PARAMETERS, DO_GMULT )

C  Include files
C  =============

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and multiplier variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of linearized setup and solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'
	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'
	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given offgrid indices, varying layer index and number of parameters

	INTEGER		 N, UT, K, K_PARAMETERS
	LOGICAL		 DO_GMULT

C  local variables
C  ---------------

	INTEGER		 Q, AA
	DOUBLE PRECISION SD, SU, T0, TD, TU 
	DOUBLE PRECISION ZX_DN, ZX_UP, ZW, WX, WDEL
	DOUBLE PRECISION L_ZX_DN, L_ZX_UP, L_ZW, L_WX, L_WDEL

C  Return if already done

	IF ( .NOT. DO_GMULT ) RETURN

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO AA = 1, NSTREAMS
	    DO Q = 1, K_PARAMETERS
	      L_UT_GMULT_UP(AA,UT,K,Q) = ZERO
	      L_UT_GMULT_DN(AA,UT,K,Q) = ZERO
 	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  Layer constant terms

	WX    = T_UTDN_MUBAR(UT)
	WDEL  = T_DELT_MUBAR(N)

C  For the Quasi-spherical (average secant) multipliers
C  ====================================================

	IF ( DO_QUASPHER_BEAM ) THEN

C  If the layer containing UT is the same as the varying layer

	  IF ( N .EQ. K ) THEN

	    DO AA = 1, NSTREAMS

	      ZX_DN = T_UTDN_EIGEN(AA,UT)
	      ZX_UP = T_UTUP_EIGEN(AA,UT)
	      ZW    = WDEL * ZX_UP

	      DO Q = 1, K_PARAMETERS

	        L_WDEL = L_T_DELT_MUBAR(N,K,Q)
	        L_WX   = L_T_UTDN_MUBAR(UT,K,Q)
	        L_ZX_DN = L_T_UTDN_EIGEN(AA,UT,Q)
	        L_ZX_UP = L_T_UTUP_EIGEN(AA,UT,Q)
	        L_ZW    = WDEL * L_ZX_UP + L_WDEL * ZX_UP
	      
	        SD  = ( L_ZX_DN - L_WX ) / ( ZX_DN - WX )
	        TD = L_ATERM_SAVE(AA,N,Q) + L_GAMMA_M(AA,N,N,Q) + SD
	        L_UT_GMULT_DN(AA,UT,K,Q) = TD * UT_GMULT_DN(AA,UT)

	        SU  = ( L_WX    - L_ZW ) / ( WX    - ZW )
	        TU = L_BTERM_SAVE(AA,N,Q) + L_GAMMA_P(AA,N,N,Q) + SU
	        L_UT_GMULT_UP(AA,UT,K,Q) = TU * UT_GMULT_UP(AA,UT)

	      ENDDO

	    ENDDO

C  If the varying layer is above layer N

	  ELSE IF ( N .GT. K ) THEN

	    DO Q = 1, K_PARAMETERS

	      L_WDEL = L_T_DELT_MUBAR(N,K,Q)
	      L_WX   = L_T_UTDN_MUBAR(UT,K,Q)

	      DO AA = 1, NSTREAMS

	        ZX_DN = T_UTDN_EIGEN(AA,UT)
	        ZX_UP = T_UTUP_EIGEN(AA,UT)
	        ZW    = WDEL * ZX_UP
	        L_ZW  = L_WDEL * ZX_UP
	      
	        SD  = - L_WX / ( ZX_DN - WX )
	        TD = L_GAMMA_M(AA,N,K,Q) + SD + L_INITIAL_TRANS(N,K,Q)
	        L_UT_GMULT_DN(AA,UT,K,Q) = TD * UT_GMULT_DN(AA,UT)

	        SU  = ( L_WX  - L_ZW ) / ( WX - ZW )
	        TU = L_GAMMA_P(AA,N,K,Q) + SU + L_INITIAL_TRANS(N,K,Q)
	        L_UT_GMULT_UP(AA,UT,K,Q) = TU * UT_GMULT_UP(AA,UT)

	      ENDDO

	    ENDDO

	  ENDIF

C  Plane parallel case
C  ===================

	ELSE

C  If the layer containing UT is the same as the varying layer

	  IF ( N .EQ. K ) THEN

	    DO AA = 1, NSTREAMS

	      ZX_DN = T_UTDN_EIGEN(AA,UT)
	      ZX_UP = T_UTUP_EIGEN(AA,UT)
	      ZW    = WDEL * ZX_UP

	      DO Q = 1, K_PARAMETERS

	        L_WDEL = L_T_DELT_MUBAR(N,N,Q)
	        L_WX   = L_T_UTDN_MUBAR(UT,N,Q)
	        L_ZX_DN = L_T_UTDN_EIGEN(AA,UT,Q)
	        L_ZX_UP = L_T_UTUP_EIGEN(AA,UT,Q)
	        L_ZW    = WDEL * L_ZX_UP + L_WDEL * ZX_UP
	      
	        SD  = ( L_ZX_DN - L_WX ) / ( ZX_DN - WX )
	        TD = L_ATERM_SAVE(AA,N,Q) + L_GAMMA_M(AA,N,N,Q) + SD
	        L_UT_GMULT_DN(AA,UT,K,Q) = TD * UT_GMULT_DN(AA,UT)

	        SU  = ( L_WX - L_ZW ) / ( WX - ZW )
	        TU = L_BTERM_SAVE(AA,N,Q) + L_GAMMA_P(AA,N,N,Q) + SU
	        L_UT_GMULT_UP(AA,UT,K,Q) = TU * UT_GMULT_UP(AA,UT)

	      ENDDO

	    ENDDO

C  If the varying layer is above layer N

	  ELSE IF ( N .GT. K ) THEN

	    DO Q = 1, K_PARAMETERS
	      T0 = L_INITIAL_TRANS(N,K,Q)
	      DO AA = 1, NSTREAMS	      
	        L_UT_GMULT_DN(AA,UT,K,Q) = T0 * UT_GMULT_DN(AA,UT)
	        L_UT_GMULT_UP(AA,UT,K,Q) = T0 * UT_GMULT_UP(AA,UT)
	      ENDDO
	    ENDDO

	  ENDIF

	ENDIF

C  Finish

	END

C

	SUBROUTINE L_WHOLELAYER_GMULT_UP ( N, K, K_PARAMETERS )

C  Include files
C  =============

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

C  subroutine arguments
C  --------------------

C  Given layer index and number of parameters

	INTEGER		 N, K, K_PARAMETERS

C  Local variables
C  ---------------

	INTEGER		 UM, AA, Q
        DOUBLE PRECISION WDEL, L_WDEL, CONST, T1, L_SD, L_SU

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO AA = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        L_SGMULT_UP(AA,UM,Q) = ZERO
	        L_SGMULT_DN(AA,UM,Q) = ZERO
	      ENDDO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  transmittance

	WDEL  = T_DELT_MUBAR(N)
	CONST = INITIAL_TRANS(N)

C  For N = K, Quasi-spherical and plane parallel the same.

	IF ( N .EQ. K ) THEN

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO AA = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        L_WDEL = L_T_DELT_MUBAR(N,K,Q)
	        L_SD = CONST * L_HMULT_2(AA,UM,N,Q)
	        L_SD = L_SD - L_EMULT_UP(UM,N,N,Q)
	        T1 = L_ATERM_SAVE(AA,N,Q) + L_GAMMA_M(AA,N,N,Q)
	        T1 = T1 * SGMULT_UP_DN(AA,UM,N)
	        L_SGMULT_DN(AA,UM,Q) = AGM(AA,N) * L_SD + T1
	        L_SU = - CONST * ( L_HMULT_1(AA,UM,N,Q) * WDEL +
     &                               HMULT_1(AA,UM,N)   * L_WDEL )
	        L_SU =  L_SU + L_EMULT_UP(UM,N,N,Q)
	        T1 = L_BTERM_SAVE(AA,N,Q) + L_GAMMA_P(AA,N,N,Q)
	        T1 = T1 * SGMULT_UP_UP(AA,UM,N)
	        L_SGMULT_UP(AA,UM,Q) = BGP(AA,N) * L_SU + T1
	      ENDDO
	    ENDDO
	  ENDDO

C  Case where N > K

	ELSE IF ( N .GT. K ) THEN

C  Quasi-spherical

	  IF ( DO_QUASPHER_BEAM ) THEN

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO AA = 1, NSTREAMS
	        DO Q = 1, K_PARAMETERS
	          L_SD = HAQ(N,K,Q) * HMULT_2(AA,UM,N) -
     &                               L_EMULT_UP(UM,N,K,Q)
	          T1 = L_GAMMA_M(AA,N,K,Q) * SGMULT_UP_DN(AA,UM,N)
	          L_SGMULT_DN(AA,UM,Q) = AGM(AA,N) * L_SD + T1
	          L_SU = HBQ(N,K,Q) * HMULT_1(AA,UM,N) +
     &                               L_EMULT_UP(UM,N,K,Q)
	          T1 = L_GAMMA_P(AA,N,K,Q) * SGMULT_UP_UP(AA,UM,N)
	          L_SGMULT_UP(AA,UM,Q) = BGP(AA,N) * L_SU + T1
	        ENDDO
	      ENDDO
	    ENDDO

C  Plane parallel case

          ELSE

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO AA = 1, NSTREAMS
	        DO Q = 1, K_PARAMETERS
	          T1 = L_INITIAL_TRANS(N,K,Q)
	          L_SGMULT_DN(AA,UM,Q) = T1 * SGMULT_UP_DN(AA,UM,N)
	          L_SGMULT_UP(AA,UM,Q) = T1 * SGMULT_UP_UP(AA,UM,N)
	        ENDDO
	      ENDDO
	    ENDDO

	  ENDIF

	ENDIF

C  Finish

	END

C

	SUBROUTINE L_WHOLELAYER_GMULT_DN ( N, K, K_PARAMETERS )

C  Include files
C  =============

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

C  subroutine arguments
C  --------------------

C  Given layer index and number of parameters

	INTEGER		 N, K, K_PARAMETERS

C  Local variables
C  ---------------

	INTEGER		 UM, AA, Q
        DOUBLE PRECISION WDEL, L_WDEL, CONST, T1, L_SD, L_SU

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO AA = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        L_SGMULT_UP(AA,UM,Q) = ZERO
	        L_SGMULT_DN(AA,UM,Q) = ZERO
	      ENDDO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  transmittance

	WDEL = T_DELT_MUBAR(N)
	CONST = INITIAL_TRANS(N)

C  For N = K, Quasi-spherical and plane parallel the same.

	IF ( N .EQ. K ) THEN

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO AA = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        L_WDEL = L_T_DELT_MUBAR(N,K,Q)
	        L_SD = CONST * L_HMULT_1(AA,UM,N,Q)
	        L_SD = L_SD - L_EMULT_DN(UM,N,N,Q)
	        T1 = L_ATERM_SAVE(AA,N,Q) + L_GAMMA_M(AA,N,N,Q)
	        T1 = T1 * SGMULT_DN_DN(AA,UM,N)
	        L_SGMULT_DN(AA,UM,Q) = AGM(AA,N) * L_SD + T1
	        L_SU = - CONST * ( L_HMULT_2(AA,UM,N,Q) * WDEL +
     &                               HMULT_2(AA,UM,N)   * L_WDEL )
	        L_SU =  L_SU + L_EMULT_DN(UM,N,N,Q)
	        T1 = L_BTERM_SAVE(AA,N,Q) + L_GAMMA_P(AA,N,N,Q)
	        T1 = T1 * SGMULT_DN_UP(AA,UM,N)
	        L_SGMULT_UP(AA,UM,Q) = BGP(AA,N) * L_SU + T1
	      ENDDO
	    ENDDO
	  ENDDO

C  Case where N > K

	ELSE IF ( N .GT. K ) THEN

C  Quasi-spherical

	  IF ( DO_QUASPHER_BEAM ) THEN

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO AA = 1, NSTREAMS
	        DO Q = 1, K_PARAMETERS
	          L_SD = HAQ(N,K,Q) * HMULT_1(AA,UM,N) -
     &                             L_EMULT_DN(UM,N,K,Q)
	          T1 = L_GAMMA_M(AA,N,K,Q) * SGMULT_DN_DN(AA,UM,N)
	          L_SGMULT_DN(AA,UM,Q) = AGM(AA,N) * L_SD + T1
	          L_SU = HBQ(N,K,Q) * HMULT_2(AA,UM,N) +
     &                             L_EMULT_DN(UM,N,K,Q)
	          T1 = L_GAMMA_P(AA,N,K,Q) * SGMULT_DN_UP(AA,UM,N)
	          L_SGMULT_UP(AA,UM,Q) = BGP(AA,N) * L_SU + T1
	        ENDDO
	      ENDDO
	    ENDDO

C  Plane parallel case

          ELSE

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO AA = 1, NSTREAMS
	        DO Q = 1, K_PARAMETERS
	          T1 = L_INITIAL_TRANS(N,K,Q)
	          L_SGMULT_DN(AA,UM,Q) = T1 * SGMULT_DN_DN(AA,UM,N)
	          L_SGMULT_UP(AA,UM,Q) = T1 * SGMULT_DN_UP(AA,UM,N)
	        ENDDO
	      ENDDO
	    ENDDO

	  ENDIF

	ENDIF

C  Finish

	END

C

	SUBROUTINE L_PARTLAYER_GMULT_UP ( N, UT, K, K_PARAMETERS )

C  Include files
C  =============

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

C  subroutine arguments
C  --------------------

C  Given layer index and number of parameters

	INTEGER		 N, UT, K, K_PARAMETERS

C  Local variables
C  ---------------

	INTEGER		 UM, AA, Q
        DOUBLE PRECISION WDEL, L_WDEL, CONST, T1, L_SD, L_SU

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO AA = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        L_SGMULT_UP(AA,UM,Q) = ZERO
	        L_SGMULT_DN(AA,UM,Q) = ZERO
	      ENDDO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  transmittance

	WDEL = T_DELT_MUBAR(N)
	CONST = INITIAL_TRANS(N)

C  For N = K, Quasi-spherical and plane parallel the same.

	IF ( N .EQ. K ) THEN

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO AA = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        L_WDEL = L_T_DELT_MUBAR(N,K,Q)
	        L_SD = CONST * L_UT_HMULT_UP_DN(AA,UM,UT,Q)
	        L_SD = L_SD - L_UT_EMULT_UP(UM,UT,N,Q)
	        T1 = L_ATERM_SAVE(AA,N,Q) + L_GAMMA_M(AA,N,N,Q)
	        T1 = T1 * UT_SGMULT_UP_DN(AA,UM,UT)
	        L_SGMULT_DN(AA,UM,Q) = AGM(AA,N) * L_SD + T1
	        L_SU = - CONST * ( L_UT_HMULT_UP_UP(AA,UM,UT,Q) * WDEL +
     &                          UT_HMULT_UP_UP(AA,UM,UT)   * L_WDEL )
	        L_SU =  L_SU + L_UT_EMULT_UP(UM,UT,N,Q)
	        T1 = L_BTERM_SAVE(AA,N,Q) + L_GAMMA_P(AA,N,N,Q)
	        T1 = T1 * UT_SGMULT_UP_UP(AA,UM,UT)
	        L_SGMULT_UP(AA,UM,Q) = BGP(AA,N) * L_SU + T1
	      ENDDO
	    ENDDO
	  ENDDO

C  Case where N > K

	ELSE IF ( N .GT. K ) THEN

C  Quasi-spherical

	  IF ( DO_QUASPHER_BEAM ) THEN

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO AA = 1, NSTREAMS
	        DO Q = 1, K_PARAMETERS
	          L_SD = HAQ(N,K,Q) * UT_HMULT_UP_DN(AA,UM,UT) - 
     &                              L_UT_EMULT_UP(UM,UT,K,Q)
	          T1 = L_GAMMA_M(AA,N,K,Q) * UT_SGMULT_UP_DN(AA,UM,UT)
	          L_SGMULT_DN(AA,UM,Q) = AGM(AA,N) * L_SD + T1
	          L_SU = HBQ(N,K,Q) * UT_HMULT_UP_UP(AA,UM,UT) +
     &                              L_UT_EMULT_UP(UM,UT,K,Q)
	          T1 = L_GAMMA_P(AA,N,K,Q) * UT_SGMULT_UP_UP(AA,UM,UT)
	          L_SGMULT_UP(AA,UM,Q) = BGP(AA,N) * L_SU + T1
	        ENDDO
	      ENDDO
	    ENDDO

C  Plane parallel case

          ELSE

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO AA = 1, NSTREAMS
	        DO Q = 1, K_PARAMETERS
	          T1 = L_INITIAL_TRANS(N,K,Q)
	          L_SGMULT_DN(AA,UM,Q) = T1 * UT_SGMULT_UP_DN(AA,UM,UT)
	          L_SGMULT_UP(AA,UM,Q) = T1 * UT_SGMULT_UP_UP(AA,UM,UT)
	        ENDDO
	      ENDDO
	    ENDDO

	  ENDIF

	ENDIF

C  Finish

	END


C

	SUBROUTINE L_PARTLAYER_GMULT_DN ( N, UT, K, K_PARAMETERS )

C  Include files
C  =============

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

C  subroutine arguments
C  --------------------

C  Given layer index and number of parameters

	INTEGER		 N, UT, K, K_PARAMETERS

C  Local variables
C  ---------------

	INTEGER		 UM, AA, Q
        DOUBLE PRECISION WDEL, L_WDEL, CONST, T1, L_SD, L_SU

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO AA = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        L_SGMULT_UP(AA,UM,Q) = ZERO
	        L_SGMULT_DN(AA,UM,Q) = ZERO
	      ENDDO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  transmittance

	WDEL = T_DELT_MUBAR(N)
	CONST = INITIAL_TRANS(N)

C  For N = K, Quasi-spherical and plane parallel the same.

	IF ( N .EQ. K ) THEN

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO AA = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        L_WDEL = L_T_DELT_MUBAR(N,K,Q)
	        L_SD = CONST * L_UT_HMULT_DN_DN(AA,UM,UT,Q)
	        L_SD = L_SD - L_UT_EMULT_DN(UM,UT,N,Q)
	        T1 = L_ATERM_SAVE(AA,N,Q) + L_GAMMA_M(AA,N,N,Q)
	        T1 = T1 * UT_SGMULT_DN_DN(AA,UM,UT)
	        L_SGMULT_DN(AA,UM,Q) = AGM(AA,N) * L_SD + T1
	        L_SU = - CONST * ( L_UT_HMULT_DN_UP(AA,UM,UT,Q) * WDEL +
     &                        UT_HMULT_DN_UP(AA,UM,UT)   * L_WDEL )
	        L_SU =  L_SU + L_UT_EMULT_DN(UM,UT,N,Q)
	        T1 = L_BTERM_SAVE(AA,N,Q) + L_GAMMA_P(AA,N,N,Q)
	        T1 = T1 * UT_SGMULT_DN_UP(AA,UM,UT)
	        L_SGMULT_UP(AA,UM,Q) = BGP(AA,N) * L_SU + T1
	      ENDDO
	    ENDDO
	  ENDDO

C  Case where N > K

	ELSE IF ( N .GT. K ) THEN

C  Quasi-spherical

	  IF ( DO_QUASPHER_BEAM ) THEN

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO AA = 1, NSTREAMS
	        DO Q = 1, K_PARAMETERS
	          L_SD = HAQ(N,K,Q) * UT_HMULT_DN_DN(AA,UM,UT) - 
     &                              L_UT_EMULT_DN(UM,UT,K,Q)
	          T1 = L_GAMMA_M(AA,N,K,Q) * UT_SGMULT_DN_DN(AA,UM,UT)
	          L_SGMULT_DN(AA,UM,Q) = AGM(AA,N) * L_SD + T1
	          L_SU = HBQ(N,K,Q) * UT_HMULT_DN_UP(AA,UM,UT) +
     &                              L_UT_EMULT_DN(UM,UT,K,Q)
	          T1 = L_GAMMA_P(AA,N,K,Q) * UT_SGMULT_DN_UP(AA,UM,UT)
	          L_SGMULT_UP(AA,UM,Q) = BGP(AA,N) * L_SU + T1
	        ENDDO
	      ENDDO
	    ENDDO

C  Plane parallel case

          ELSE

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DO AA = 1, NSTREAMS
	        DO Q = 1, K_PARAMETERS
	          T1 = L_INITIAL_TRANS(N,K,Q)
	          L_SGMULT_DN(AA,UM,Q) = T1 * UT_SGMULT_DN_DN(AA,UM,UT)
	          L_SGMULT_UP(AA,UM,Q) = T1 * UT_SGMULT_DN_UP(AA,UM,UT)
	        ENDDO
	      ENDDO
	    ENDDO

	  ENDIF

	ENDIF

C  Finish

	END
