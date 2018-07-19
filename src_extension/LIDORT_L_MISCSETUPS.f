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

	SUBROUTINE LIDORT_L_MISCSETUPS

C  miscellaneous setup operations for linearized quantities

C  Deltam scaling of variational quantities

	CALL LIDORT_L_DELTAMSCALE_TOTAL

C  initialise single scatter albedo variational quantities

	CALL LIDORT_L_SSALBINIT_TOTAL

C  Transmittance factors variational quantities

	CALL LIDORT_L_PREPTRANS

C  Finish

	END

C

	SUBROUTINE LIDORT_L_DELTAMSCALE_TOTAL

C  include files
C  -------------

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  Additional include files

	INCLUDE '../include_e/LIDORT_L_CONTROL.VARS'
	INCLUDE '../include_e/LIDORT_L_GEOPHYS.VARS'

C  include files of set up variables (output to this module)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'

C  local variables
C  ---------------

	DOUBLE PRECISION BLD, DL, OF1, F, F1, UQ, EQ, FZM, T1
	DOUBLE PRECISION ZMQ, ZLQ, UZQ_SCALE, ZQ_SCALE, UQ_SCALE
	INTEGER		 N, Q, NM1, L
	LOGICAL		 LOOP

C  Set phasfunc linearization flag

	NM1 = NMOMENTS+1
	IF ( DO_LAYER_LINEARIZATION ) THEN
	  DO N = 1, NLAYERS
	    IF ( LAYER_VARY_FLAG(N) ) THEN
	      DO Q = 1, LAYER_VARY_NUMBER(N)
	        LOOP = .TRUE.
	        L = 0
	        DO WHILE (LOOP.AND.L.LT.NM1)
	          L = L + 1
	          LOOP = ( DABS(PHASMOM_VARS_TOTAL_INPUT(Q,L,N))
     &                    .LT.1000.0*SMALLNUM )
	        ENDDO
	        DO_PHASFUNC_VARIATION(N,Q) = .NOT.LOOP
	      ENDDO
	    ENDIF
	  ENDDO
	ENDIF

C  DELTAM SCALING
C  ==============

	IF ( DO_DELTAM_SCALING ) THEN

	  IF ( DO_LAYER_LINEARIZATION ) THEN

	    NM1 = NMOMENTS+1
	    DO N = 1, NLAYERS
	      IF ( LAYER_VARY_FLAG(N) ) THEN
	        OF1 = ( ONE - FAC1(N) ) / FAC1(N)
	        F   = TRUNC_FACTOR(N)
	        F1  = ONE - F
	        DO Q = 1, LAYER_VARY_NUMBER(N)

C  scale phase function linearization additionally (versions 2.2 +)

	          IF ( DO_PHASFUNC_VARIATION(N,Q) ) THEN

	            UQ = OMEGA_VARS_TOTAL_INPUT(Q,N)
	            EQ = EXT_VARS_INPUT(Q,N)
	            ZMQ = PHASMOM_VARS_TOTAL_INPUT(Q,NM1,N)
	            FZM = F * ZMQ
	            L_TRUNC_FACTOR(Q,N) = FZM
	            UZQ_SCALE = ( UQ + ZMQ ) * OF1
	            ZQ_SCALE = FZM / F1
	            OMEGA_VARS_TOTAL(Q,N) = UQ + UZQ_SCALE - ZQ_SCALE
	            EXT_VARS(Q,N)         = EQ - UZQ_SCALE
	            PHASMOM_VARS_TOTAL(Q,0,N) = ZERO
	            DO L = 1, NMOMENTS
	              ZLQ = PHASMOM_VARS_TOTAL_INPUT(Q,L,N)
	              DL = DFLOAT(2*L+1)
	              BLD = PHASMOMS_TOTAL_INPUT(L,N) / DL
	              T1 = ( BLD*ZLQ - FZM ) / ( BLD - F )
	              PHASMOM_VARS_TOTAL(Q,L,N) = T1 + ZQ_SCALE
	            ENDDO

C  No phase function linearization
C   Zero all PHASMOM_VARS_TOTAL quantities now;

	          ELSE

	            UQ = OMEGA_VARS_TOTAL_INPUT(Q,N)
	            EQ = EXT_VARS_INPUT(Q,N)
	            L_TRUNC_FACTOR(Q,N) = ZERO
	            UQ_SCALE = UQ * OF1
	            OMEGA_VARS_TOTAL(Q,N) = UQ + UQ_SCALE
	            EXT_VARS(Q,N)         = EQ - UQ_SCALE
	            DO L = 0, NMOMENTS
	              PHASMOM_VARS_TOTAL(Q,L,N) = ZERO
	            ENDDO

	          ENDIF

	        ENDDO
	      ENDIF
	    ENDDO

	  ENDIF

C  NO DELTAM SCALING
C  =================

C  move input geophysical variables to Workspace quantities

C  Version 2.2
C  New code to deal with phase function linearization

	ELSE

	  IF ( DO_LAYER_LINEARIZATION ) THEN

	    DO N = 1, NLAYERS
	      IF ( LAYER_VARY_FLAG(N) ) THEN
	        DO Q = 1, LAYER_VARY_NUMBER(N)
	          L_TRUNC_FACTOR(Q,N) = ZERO
	          EXT_VARS(Q,N)         = EXT_VARS_INPUT(Q,N)
	          OMEGA_VARS_TOTAL(Q,N) = OMEGA_VARS_TOTAL_INPUT(Q,N)
	          IF ( DO_PHASFUNC_VARIATION(N,Q) ) THEN
	            DO L = 0, MAXMOMENT
	              PHASMOM_VARS_TOTAL(Q,L,N) = 
     &                         PHASMOM_VARS_TOTAL_INPUT(Q,L,N)
	            ENDDO
	          ELSE
	            DO L = 0, MAXMOMENT
	              PHASMOM_VARS_TOTAL(Q,L,N) = ZERO
	            ENDDO
	          ENDIF
	        ENDDO
	      ENDIF
	    ENDDO

	  ENDIF

	ENDIF

C  Finish module

	END

C
    
	SUBROUTINE LIDORT_L_SSALBINIT_TOTAL

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'
	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  Additional include files

	INCLUDE '../include_e/LIDORT_L_CONTROL.VARS'
	INCLUDE '../include_e/LIDORT_L_GEOPHYS.VARS'

C  include files of set up variables (output to this module)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'

C  local variables

	INTEGER		 N, L, Q
	DOUBLE PRECISION VAR_L

C  phase moment-weighted OMEGA and linearizations
C  Including phase function linearization (version 2.2)

	IF ( DO_LAYER_LINEARIZATION ) THEN

	  DO N = 1, NLAYERS
  	    IF ( LAYER_VARY_FLAG(N) ) THEN
	      DO Q = 1, LAYER_VARY_NUMBER(N)
  	        DO L = 0, NMOMENTS
	          VAR_L = OMEGA_VARS_TOTAL(Q,N) +
     &                      PHASMOM_VARS_TOTAL(Q,L,N)
	          L_OMEGA_MOMS(N,L,Q) = OMEGA_MOMS(N,L) * VAR_L
	        ENDDO
	      ENDDO
	    ENDIF
	  ENDDO

	ENDIF

C  Finish

	END

C

	SUBROUTINE LIDORT_L_PREPTRANS

C  include files
C  -------------

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of set up variables (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  Additional include files

	INCLUDE '../include_e/LIDORT_L_GEOPHYS.VARS'

C  include files of linearized setup variables (output to this module)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'

C  local variables
C  ---------------

	INTEGER		 N, Q, UT, UM, K
	DOUBLE PRECISION VD, VU, TRANS, UX, TRANS_D, TRANS_U
	DOUBLE PRECISION WDEL, VAR, RHO, FAC, DELT, LAMDA
	DOUBLE PRECISION CONST, L_IT, L_WDEL

C  Ancillary quantities
C  --------------------

	DO N = 1, NLAYERS
	  DO Q = 1, LAYER_VARY_NUMBER(N)
	    VQ(N,Q) = EXT_VARS(Q,N)
	  ENDDO
	ENDDO

C  linearization of constants (relative variation)

	DO N = 1, NLAYERS
	  IF  (N.LE.LAYER_PIS_CUTOFF ) THEN
	    DO Q = 1, LAYER_VARY_NUMBER(N)
	      L_INITIAL_TRANS(N,N,Q) = ZERO
	    ENDDO
	    IF ( N .GT. 1 ) THEN
	      DO K = 1, N-1
	        DO Q = 1, LAYER_VARY_NUMBER(K)
	          L_INITIAL_TRANS(N,K,Q) = - VQ(K,Q) * TAUTHICK(N-1,K)
	        ENDDO
	      ENDDO
	    ENDIF
	  ELSE
	    DO K = 1, N
	      DO Q = 1, LAYER_VARY_NUMBER(K)
	        L_INITIAL_TRANS(N,K,Q) = ZERO
	      ENDDO
	    ENDDO
	  ENDIF
	ENDDO

C  linearization of LAMDA secants for quasi spherical case

	IF ( DO_QUASPHER_BEAM ) THEN

	  DO N = 1, NLAYERS

	    IF ( N .EQ. 1 ) THEN
	      DO Q = 1, LAYER_VARY_NUMBER(N)
	        L_AVERAGE_SECANT(N,N,Q) = ZERO
	      ENDDO
	    ELSE
	      IF  (N.LE.LAYER_PIS_CUTOFF ) THEN
	        DELT = DELTA(N)
	        LAMDA = AVERAGE_SECANT(N)
	        FAC = ( TAUTHICK(N,N) / DELT ) - LAMDA
	        DO Q = 1, LAYER_VARY_NUMBER(N)
	          L_AVERAGE_SECANT(N,N,Q) = VQ(N,Q) * FAC
	        ENDDO
	        DO K = 1, N-1
	          FAC = ( TAUTHICK(N,K) - TAUTHICK(N-1,K) ) / DELT
	          DO Q = 1, LAYER_VARY_NUMBER(K)
	            L_AVERAGE_SECANT(N,K,Q) = VQ(K,Q) * FAC
	          ENDDO
	        ENDDO
	      ELSE
	        DO K = 1, N
	          DO Q = 1, LAYER_VARY_NUMBER(K)
	            L_AVERAGE_SECANT(N,K,Q) = ZERO
	          ENDDO
	        ENDDO
	      ENDIF
	    ENDIF

	  ENDDO

	ENDIF

C  Whole layer Transmittance factors
C  ---------------------------------

	DO N = 1, NLAYERS

	  WDEL = T_DELT_MUBAR(N)
	  VAR = - DELTA(N) * WDEL
	  LAMDA = AVERAGE_SECANT(N)
	  FAC = VAR * AVERAGE_SECANT(N)

C  Quasi_spherical

	  IF ( DO_QUASPHER_BEAM ) THEN

	    IF ( N .EQ. 1 ) THEN
	      DO Q = 1, LAYER_VARY_NUMBER(N)
	        L_T_DELT_MUBAR(N,N,Q) = FAC * VQ(N,Q)
	      ENDDO
	    ELSE
	      IF  (N.LE.LAYER_PIS_CUTOFF ) THEN
	        DO Q = 1, LAYER_VARY_NUMBER(N)
	          RHO = L_AVERAGE_SECANT(N,N,Q)
	          L_T_DELT_MUBAR(N,N,Q) = VQ(N,Q) * FAC + VAR * RHO
	        ENDDO
	        DO K = 1, N-1
	          DO Q = 1, LAYER_VARY_NUMBER(K)
	            RHO = L_AVERAGE_SECANT(N,K,Q)
	            L_T_DELT_MUBAR(N,K,Q) = VAR * RHO
	          ENDDO
	        ENDDO
	      ELSE
	        DO K = 1, N
	          DO Q = 1, LAYER_VARY_NUMBER(K)
	            L_T_DELT_MUBAR(N,K,Q) = ZERO
	          ENDDO
	        ENDDO
	      ENDIF
	    ENDIF

C  Plane-parallel

	  ELSE

	    IF  (N.LE.LAYER_PIS_CUTOFF ) THEN
	      DO Q = 1, LAYER_VARY_NUMBER(N)
	        L_T_DELT_MUBAR(N,N,Q) = FAC * VQ(N,Q)
	      ENDDO
	      DO K = 1, N-1
	        DO Q = 1, LAYER_VARY_NUMBER(K)
	          L_T_DELT_MUBAR(N,K,Q) = ZERO
	        ENDDO
	      ENDDO
	    ELSE
	      DO K = 1, N
	        DO Q = 1, LAYER_VARY_NUMBER(K)
	          L_T_DELT_MUBAR(N,K,Q) = ZERO
	        ENDDO
	      ENDDO
	    ENDIF

	  ENDIF

	ENDDO

C  Partial layer transmittance factors (for off-grid optical depths)
C  -----------------------------------------------------------------

	DO UT = 1, N_OFFGRID_USERTAUS

	  N   = OFFGRID_UTAU_LAYERIDX(UT)
	  VAR = - OFFGRID_UTAU_VALUES(UT) * T_UTDN_MUBAR(UT)
	  FAC = VAR * AVERAGE_SECANT(N)

C  Quasi_spherical

	  IF ( DO_QUASPHER_BEAM ) THEN

	    IF ( N .EQ. 1 ) THEN
	      DO Q = 1, LAYER_VARY_NUMBER(N)
	        L_T_UTDN_MUBAR(UT,N,Q) = FAC * VQ(N,Q)
	      ENDDO
	    ELSE
	      IF ( N. LE. LAYER_PIS_CUTOFF ) THEN
	        DO Q = 1, LAYER_VARY_NUMBER(N)
	          RHO = L_AVERAGE_SECANT(N,N,Q)
	          L_T_UTDN_MUBAR(UT,N,Q) = VQ(N,Q) * FAC + VAR * RHO
	        ENDDO
	        DO K = 1, N-1
	          DO Q = 1, LAYER_VARY_NUMBER(K)
	            RHO = L_AVERAGE_SECANT(N,K,Q)
	            L_T_UTDN_MUBAR(UT,K,Q) = VAR * RHO
	          ENDDO
	        ENDDO
	      ELSE
	        DO K = 1, N
	          DO Q = 1, LAYER_VARY_NUMBER(K)
	            L_T_UTDN_MUBAR(UT,K,Q) = ZERO
	          ENDDO
	        ENDDO
	      ENDIF
	    ENDIF

C  Plane-parallel

	  ELSE

	    IF ( N. LE. LAYER_PIS_CUTOFF ) THEN
	      DO Q = 1, LAYER_VARY_NUMBER(N)
	        L_T_UTDN_MUBAR(UT,N,Q) = FAC * VQ(N,Q)
	      ENDDO
	      DO K = 1, N-1
	        DO Q = 1, LAYER_VARY_NUMBER(K)
	          L_T_UTDN_MUBAR(UT,K,Q) = ZERO
	        ENDDO
	      ENDDO
	    ELSE
	      DO K = 1, N
	        DO Q = 1, LAYER_VARY_NUMBER(K)
	          L_T_UTDN_MUBAR(UT,K,Q) = ZERO
	        ENDDO
	      ENDDO
	    ENDIF

	  ENDIF

	ENDDO

C  Linearization of Transmittance factors for User Streams
C  =======================================================

	IF ( .NOT. DO_USER_STREAMS  ) RETURN

C  help arrays (saves time later on)

	IF ( .NOT. DO_CLASSICAL_SOLUTION ) THEN
	  IF ( DO_QUASPHER_BEAM ) THEN
	    DO N = 1, NLAYERS
	      IF ( N .GT. 1 ) THEN
	        WDEL = T_DELT_MUBAR(N)
	        CONST = INITIAL_TRANS(N)
	        DO K = 1, N-1
	          DO Q = 1, LAYER_VARY_NUMBER(K)
	            L_WDEL = L_T_DELT_MUBAR(N,K,Q)
	            L_IT   = L_INITIAL_TRANS(N,K,Q)
	            HAQ(N,K,Q) = CONST * L_IT
	            HBQ(N,K,Q) = - CONST *  ( L_IT * WDEL + L_WDEL )
	          ENDDO
	        ENDDO
	      ENDIF
	    ENDDO
	  ENDIF
	ENDIF

C  Whole Layer transmittance factors
C  ---------------------------------

	DO N = 1, NLAYERS
	  IF ( LAYER_VARY_FLAG(N) ) THEN
	    DO UM = 1, N_USER_STREAMS
	      TRANS = T_DELT_USERM(N,UM) * USER_SECANTS(UM) * DELTA(N)
	      DO Q = 1, LAYER_VARY_NUMBER(N)
	        L_T_DELT_USERM(N,UM,Q) = - TRANS * VQ(N,Q)
	      ENDDO
	    ENDDO
	  ENDIF
	ENDDO

C  Partial Layer transmittance factors for off-grid optical depths
C  ---------------------------------------------------------------

	DO UT = 1, N_OFFGRID_USERTAUS
	  N  = OFFGRID_UTAU_LAYERIDX(UT)
	  UX = OFFGRID_UTAU_VALUES(UT)
	  IF ( LAYER_VARY_FLAG(N) ) THEN
	    DO UM = 1, N_USER_STREAMS
	      TRANS_D = T_UTDN_USERM(UT,UM) * USER_SECANTS(UM)
	      TRANS_U = T_UTUP_USERM(UT,UM) * USER_SECANTS(UM)
	      DO Q = 1, LAYER_VARY_NUMBER(N)
	        VD = EXT_VARS(Q,N) * UX
	        VU = EXT_VARS(Q,N) * ( DELTA(N) - UX )
	        L_T_UTDN_USERM(UT,UM,Q) = - TRANS_D * VD
	        L_T_UTUP_USERM(UT,UM,Q) = - TRANS_U * VU
	      ENDDO
	    ENDDO
	  ENDIF
	ENDDO

C  Finish

	END
