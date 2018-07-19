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

	SUBROUTINE LIDORT_ATMOSWF_COLUMN_SETUP
     I     ( DO_INCLUDE_ALBEDO, DO_INCLUDE_DIRECTBEAM,
     I       MODIFIED_BCL3,     MODIFIED_BCL4,
     I       LAYER_TO_VARY, N_PARAMETERS, FOURIER,
     I       R2, DIRECT_BEAM )

C  Standard Include files (all inputs)
C  ----------------------

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  include files of model variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of main solution and setup variables

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  Extended Include files
C  ----------------------

C  include files of linearized setup variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'

C  include files of linearized solution variables (output stored here)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  input arguments
C  ---------------

	INTEGER          N_PARAMETERS, LAYER_TO_VARY
	INTEGER		 FOURIER
	LOGICAL		 MODIFIED_BCL3, MODIFIED_BCL4
	LOGICAL		 DO_INCLUDE_ALBEDO,  DO_INCLUDE_DIRECTBEAM
	DOUBLE PRECISION DIRECT_BEAM(MAXSTRM), R2

C  local variables
C  ---------------

	INTEGER		 Q,AA,N,N1,I,I1,CM,C0
	DOUBLE PRECISION CPOS, CNEG, L_HOM, L_BEAM, FAC

	DOUBLE PRECISION R2_L_BEAM(MAXSTRM,MAX_PARAMETERS)
	DOUBLE PRECISION R2_L_HOMP(MAXSTRM,MAXSTRM,MAX_PARAMETERS)
 	DOUBLE PRECISION R2_L_HOMM(MAXSTRM,MAXSTRM,MAX_PARAMETERS)

	LOGICAL		 REGULAR_BCL3, REGULAR_BCL4

C  initialise
C  ----------

C  zero the results vectors

	DO I = 1, MAXTOTAL
	  DO Q = 1, MAX_PARAMETERS
	    COL2_WF(I,Q) = ZERO
	  ENDDO
	ENDDO

C  Get the linearized beam solution for the varying layer

	CALL L_BEAMSOL_LTOVARY ( LAYER_TO_VARY, N_PARAMETERS)

C  complete boundary condition flags

	REGULAR_BCL3 = .NOT.MODIFIED_BCL3
	REGULAR_BCL4 = .NOT.MODIFIED_BCL4

C  BCL1 or BCL3M - top of first layer (TOA), UPPER boundary condition
C  ------------------------------------------------------------------

	N = 1

C    If this layer is the one that is varied, use MODIFIED_BCL3 (BCL3M)

	IF ( MODIFIED_BCL3 ) THEN

C  .. contribution WVAR from beam solution variations
C  .. contribution HVAR homogeneous (eigenvalue) solution variations

	  DO I = 1, NSTREAMS
	    DO Q = 1, N_PARAMETERS
	      L_BEAM = - L_WUPPER(I,N,Q)
	      L_HOM  = ZERO
	      DO AA = 1, NSTREAMS
	        CPOS = L_XPOS(I,AA,N,Q)
	        CNEG = T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) + 
     &               L_T_DELT_EIGEN(AA,N,Q) * XNEG(I,AA,N)
	        L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
	      ENDDO
	      COL2_WF(I,Q) = L_BEAM - L_HOM
	    ENDDO
	  ENDDO

C  No variation case (BCL1)

	ELSE

	  DO I = 1, NSTREAMS
	    DO Q = 1, N_PARAMETERS
	      COL2_WF(I,Q) = ZERO
	    ENDDO
	  ENDDO

	ENDIF

C  BCL2 Intermediate levels between top layer and varying layer
C  ------------------------------------------------------------

C  [not required if top layer is varying, case MODIFIED_BCL3 above]

	IF ( REGULAR_BCL3 ) THEN

C  .. nothing varying in these layers

	  DO N = 2, LAYER_TO_VARY - 1
	    N1 = N - 1
            C0  = N1*NSTR2 - NSTREAMS
	    DO I = 1, NSTR2
	      CM = C0 + I
	      DO Q = 1, N_PARAMETERS
	        COL2_WF(CM,Q) = ZERO
	      ENDDO
	    ENDDO
	  ENDDO

	ENDIF

C  BCL3 - regular upper boundary condition for layer that is varying
C  -----------------------------------------------------------------

	IF ( REGULAR_BCL3 ) THEN

 	  N = LAYER_TO_VARY
	  N1  = N - 1
          C0  = N1*NSTR2 - NSTREAMS

C  .. contribution WVAR from beam solution variations
C  .. contribution HVAR homogeneous (eigenvalue) solution variations

	  DO I = 1, NSTR2
	    CM = C0 + I
	    DO Q = 1, N_PARAMETERS
	      L_BEAM  = + L_WUPPER(I,N,Q)
	      L_HOM = ZERO
	      DO AA = 1, NSTREAMS
	        CPOS = L_XPOS(I,AA,N,Q)
	        CNEG = T_DELT_EIGEN(AA,N)   * L_XNEG(I,AA,N,Q) + 
     &               L_T_DELT_EIGEN(AA,N,Q) * XNEG(I,AA,N)
	        L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
	      ENDDO
	      COL2_WF(CM,Q) = L_BEAM + L_HOM
	    ENDDO
	  ENDDO

	ENDIF
 
C  BCL4 - LOWER boundary condition for varying layer
C  -------------------------------------------------

C   special case when layer-to-vary = last (albedo) layer is treated
C   separately below under MODIFIED BCL4.

	IF ( REGULAR_BCL4 ) THEN

 	  N = LAYER_TO_VARY
	  N1 = N + 1
          C0  = N*NSTR2 - NSTREAMS

C  Get the linearized beam solution for the next layer

	  CALL L_BEAMSOL_NONVARY ( N1, LAYER_TO_VARY, N_PARAMETERS )

C  .. 2 contributions to WVAR from beam solution variations BEAM_V and BEAM_U 
C  .. contribution HVAR homogeneous (eigenvalue) solution variations

	  DO I = 1, NSTR2
	    CM = C0 + I
	    DO Q = 1, N_PARAMETERS
	      L_BEAM = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)	
	      L_HOM  = ZERO
	      DO AA = 1, NSTREAMS
	        CNEG = L_XNEG(I,AA,N,Q)
	        CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I,AA,N,Q) + 
     &               L_T_DELT_EIGEN(AA,N,Q) * XPOS(I,AA,N)
	        L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
	      ENDDO
	      COL2_WF(CM,Q) = L_BEAM - L_HOM
	    ENDDO
	  ENDDO

	ENDIF

C  BCL5 - Intermediate boundary conditions between varying layer & final layer
C  ---------------------------------------------------------------------------

	IF ( REGULAR_BCL4 ) THEN

	  DO N = LAYER_TO_VARY + 1, NLAYERS - 1

	    N1 = N + 1
            C0  = N*NSTR2 - NSTREAMS

C  Get the linearized beam solution for the next layer

	    CALL L_BEAMSOL_NONVARY ( N1, LAYER_TO_VARY, N_PARAMETERS)

C  .. contributions from beam solution (direct assign). No homog. variation

	    DO I = 1, NSTR2
	      CM = C0 + I
	      DO Q = 1, N_PARAMETERS
	        COL2_WF(CM,Q) = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)
	      ENDDO
	    ENDDO

C  end layer loop

	  ENDDO

C  end BCL5 boundary conditions

	ENDIF

C  Final layer - use BCL6 or BCL4M (last layer is varying)
C  -------------------------------------------------------

 	N = NLAYERS

C  Modified BCL4M Component loop

	IF ( MODIFIED_BCL4 ) THEN

C  get the linearized downward-reflected term

	  CALL L_INTEGRATED_DOWNREFLEC
     I     ( DO_INCLUDE_ALBEDO, MODIFIED_BCL4,
     I       FOURIER, R2, N_PARAMETERS,
     O       R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

C  Compute the solution

	  C0 = (N-1)*NSTR2 + NSTREAMS
	  DO I = 1, NSTREAMS
	    CM = C0 + I
	    I1 = I + NSTREAMS
	    DO Q = 1, N_PARAMETERS
	      L_BEAM = L_WLOWER(I1,N,Q) - R2_L_BEAM(I,Q)
	      L_HOM  = ZERO
	      DO AA = 1, NSTREAMS
	        CPOS = T_DELT_EIGEN(AA,N)   * L_XPOS(I1,AA,N,Q) + 
     &               L_T_DELT_EIGEN(AA,N,Q) * XPOS(I1,AA,N)
	        CPOS =        CPOS       - R2_L_HOMP(I,AA,Q)
	        CNEG = L_XNEG(I1,AA,N,Q) - R2_L_HOMM(I,AA,Q)
	        L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
	      ENDDO
	      COL2_WF(CM,Q) = - L_BEAM - L_HOM
	    ENDDO
	  ENDDO

C  ordinary BCL6 Component loop
 
	ELSE

C  get the linearized downward-reflected term

	  CALL L_INTEGRATED_DOWNREFLEC
     I     ( DO_INCLUDE_ALBEDO, MODIFIED_BCL4,
     I       FOURIER, R2, N_PARAMETERS,
     O       R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

C  Compute the solution

	  C0 = (N-1)*NSTR2 + NSTREAMS
	  DO I = 1, NSTREAMS
	    CM = C0 + I
	    I1 = I + NSTREAMS
	    DO Q = 1, N_PARAMETERS
	      L_BEAM = L_WLOWER(I1,N,Q) - R2_L_BEAM(I,Q)
	      COL2_WF(CM,Q) = - L_BEAM
	    ENDDO
	  ENDDO

	ENDIF

C  Add direct beam variation to Final boundary
C  -------------------------------------------

	IF ( DO_INCLUDE_DIRECTBEAM ) THEN
	  DO I = 1, NSTREAMS
	    CM = C0 + I
	    FAC = - DIRECT_BEAM(I) * TAUTHICK(N,LAYER_TO_VARY) 
	    DO Q = 1, N_PARAMETERS
	      L_BEAM = VQ(LAYER_TO_VARY,Q) * FAC
	      COL2_WF(CM,Q) = COL2_WF(CM,Q) + L_BEAM
	    ENDDO
	  ENDDO
	ENDIF

C  finish

	END

C

	SUBROUTINE L_INTEGRATED_DOWNREFLEC
     I     ( DO_INCLUDE_ALBEDO, MODIFIED_BCL4,
     I       FOURIER, R2, N_PARAMETERS,
     O       R2_L_BEAM, R2_L_HOMP, R2_L_HOMM )

C  Standard Include files (all inputs)
C  ----------------------

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  include files of Standard control, model & geophysical variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'

C  include files of Standard solution and setup variables

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  Extended Include files
C  ----------------------

C  include files of linearized setup variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'

C  include files of linearized solution variables (input)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  input arguments
C  ---------------

	INTEGER          N_PARAMETERS
	INTEGER		 FOURIER
	LOGICAL		 MODIFIED_BCL4
	LOGICAL		 DO_INCLUDE_ALBEDO
	DOUBLE PRECISION R2

C  Output arguments
C  ----------------

	DOUBLE PRECISION R2_L_BEAM(MAXSTRM,MAX_PARAMETERS)
	DOUBLE PRECISION R2_L_HOMP(MAXSTRM,MAXSTRM,MAX_PARAMETERS)
 	DOUBLE PRECISION R2_L_HOMM(MAXSTRM,MAXSTRM,MAX_PARAMETERS)

C  Local variables
C  ---------------


 	DOUBLE PRECISION PV_W(MAXSTRM,MAX_PARAMETERS)
 	DOUBLE PRECISION HV_P(MAXSTRM,MAXSTRM,MAX_PARAMETERS)
 	DOUBLE PRECISION HV_M(MAXSTRM,MAXSTRM,MAX_PARAMETERS)
	DOUBLE PRECISION HSP_U, HSM_U, H_P, H_M, HELP
	INTEGER		 AA, J, Q, N, I

C  Set up Auxiliary arrays
C  -----------------------

	IF ( DO_INCLUDE_ALBEDO ) THEN

	  N = NLAYERS

C    Modified boundary condition

	  IF ( MODIFIED_BCL4 ) THEN

	    DO J = 1, NSTREAMS
	      DO Q = 1, N_PARAMETERS
	        PV_W(J,Q) = L_WLOWER(J,N,Q) * AX(J)
	      ENDDO
	      DO AA = 1, NSTREAMS
	        DO Q = 1, N_PARAMETERS
	          HSP_U = L_XPOS(J,AA,N,Q) *   T_DELT_EIGEN(AA,N) +
     &                      XPOS(J,AA,N)   * L_T_DELT_EIGEN(AA,N,Q)
	          HSM_U = L_XNEG(J,AA,N,Q)
	          HV_P(J,AA,Q) = AX(J)*HSP_U
	          HV_M(J,AA,Q) = AX(J)*HSM_U
	        ENDDO
	      ENDDO
	    ENDDO

C     Ordinary boundary condition

	  ELSE

	    DO J = 1, NSTREAMS
	      DO Q = 1, N_PARAMETERS
	        PV_W(J,Q) = L_WLOWER(J,N,Q) * AX(J)
	      ENDDO
              DO AA = 1, NSTREAMS
	        DO Q = 1, N_PARAMETERS
	          HV_P(J,AA,Q) = ZERO
	          HV_M(J,AA,Q) = ZERO
	        ENDDO
	      ENDDO
	    ENDDO

	  ENDIF

C  Integrated Downward reflection (Calculation)
C  --------------------------------------------

C    Lambertian Reflection

	  IF ( DO_LAMBERTIAN_ALBEDO ) THEN

	    DO Q = 1, N_PARAMETERS
	      HELP = ZERO
	      DO J = 1, NSTREAMS
	        HELP = HELP + PV_W(J,Q)
	      ENDDO
	      DO I = 1, NSTREAMS
	        R2_L_BEAM(I,Q) = R2*HELP
	      ENDDO
	      DO AA = 1, NSTREAMS
	        H_P = ZERO
	        H_M = ZERO
	        DO J = 1, NSTREAMS
	          H_P = H_P + HV_P(J,AA,Q)
	          H_M = H_M + HV_M(J,AA,Q)
	        ENDDO
	        DO I = 1, NSTREAMS
	          R2_L_HOMP(I,AA,Q) = R2*H_P
	          R2_L_HOMM(I,AA,Q) = R2*H_M
	        ENDDO
	      ENDDO
	    ENDDO

C    Non-Lambertian Reflection (bidirectional)

	  ELSE

	    DO I = 1, NSTREAMS
	      DO Q = 1, N_PARAMETERS
	        HELP = ZERO
	        DO J = 1, NSTREAMS
	          HELP = HELP + PV_W(J,Q) * BIREFLEC(FOURIER,I,J)
	        ENDDO
	        R2_L_BEAM(I,Q) = R2*HELP
	        DO AA = 1, NSTREAMS
	          H_P = ZERO
	          H_M = ZERO
	          DO J = 1, NSTREAMS
	            H_P = H_P + HV_P(J,AA,Q)*BIREFLEC(FOURIER,I,J)
	            H_M = H_M + HV_M(J,AA,Q)*BIREFLEC(FOURIER,I,J)
	          ENDDO
	          R2_L_HOMP(I,AA,Q) = R2*H_P
	          R2_L_HOMM(I,AA,Q) = R2*H_M
	        ENDDO
	      ENDDO
	    ENDDO

	  ENDIF

C  no albedo - just zero the arrays (easier to code up!)
C  --------------------------------

	ELSE

	  DO I = 1, NSTREAMS
	    DO Q = 1, N_PARAMETERS
	      R2_L_BEAM(I,Q)    = ZERO
	      DO AA = 1, NSTREAMS
	        R2_L_HOMP(I,AA,Q) = ZERO
	        R2_L_HOMM(I,AA,Q) = ZERO
	      ENDDO
	    ENDDO
	  ENDDO

	ENDIF

C  Finish

	END

C

	SUBROUTINE L_BEAMSOL_LTOVARY ( N, N_PARAMETERS )

C  Linearization of Green's function beam particular integral, one layer only.
C  This is the Lyaer that contains the variation.

C  Include files
C  =============

C  input
C  -----

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_CONTROL.VARS'

C  include file of Set up stuff 

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include files of main solution variables

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include file of Multiplier coefficient variables

	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  output
C  ------

C  include file of main solution linearization variables (output)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  include file of setup linearization variables (output)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'

C  include file of linearized Multiplier coefficients (output)

	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index (input)

	INTEGER		 N

C  number of varying parameters (input)

	INTEGER		 N_PARAMETERS

C  Local variables
C  ---------------

C  help variables

	INTEGER		 AA, I, I1, Q
	DOUBLE PRECISION S_P_U, S_P_L, S_M_U, S_M_L, L_SD, L_SU 
	DOUBLE PRECISION CONST, WDEL, ZDEL, ZWDEL, T1, PQ, FQ
	DOUBLE PRECISION L_WDEL, L_ZDEL, L_ZWDEL, VAR1

C   ### NEW CODE ##################################################
C  No linearized particular solution beyond the cutoff layer.
C    [ Zero the boundary layer values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO Q = 1, N_PARAMETERS
	    DO I = 1, NSTR2
	      L_WUPPER(I,N,Q) = ZERO
	      L_WLOWER(I,N,Q) = ZERO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  Case 1. If N = LAYER_TO_VARY
C  ----------------------------

	CONST   = INITIAL_TRANS(N)
	WDEL    = T_DELT_MUBAR(N)

C  Classical solution
C  ------------------

	IF ( DO_CLASSICAL_SOLUTION ) THEN

C  ..(a) all cases

	  DO Q = 1, N_PARAMETERS
	    VAR1 = L_T_DELT_MUBAR(N,N,Q) * CONST
	    DO I = 1, NSTR2
	      L_WUPPER(I,N,Q) = CONST * L_WVEC(I,N,N,Q)
	      L_WLOWER(I,N,Q) = WDEL * L_WUPPER(I,N,Q) + VAR1*WVEC(I,N)
	    ENDDO
	  ENDDO

C  Green's function solution
C  -------------------------

	ELSE

C  Set up linearizations of GAMMA constants

C  ..(a) quasi-spherical for n > 1

	  IF ( DO_QUASPHER_BEAM .AND. N.GT.1 ) THEN

	    DO AA = 1, NSTREAMS
	      DO Q = 1, N_PARAMETERS
	        FQ = L_KEIGEN(AA,N,Q)
	        PQ = L_AVERAGE_SECANT(N,N,Q)
	        L_GAMMA_P(AA,N,N,Q) = - GAMMA_P(AA,N) * ( FQ + PQ )
	        L_GAMMA_M(AA,N,N,Q) =   GAMMA_M(AA,N) * ( FQ - PQ )
	      ENDDO
	    ENDDO

C  ..(b) plane-parallel or QS for n=1

	  ELSE

	    DO AA = 1, NSTREAMS
	      DO Q = 1, N_PARAMETERS
	        FQ = L_KEIGEN(AA,N,Q)
	        L_GAMMA_P(AA,N,N,Q) = - GAMMA_P(AA,N) * FQ 
	        L_GAMMA_M(AA,N,N,Q) =   GAMMA_M(AA,N) * FQ
	      ENDDO
	    ENDDO

	  ENDIF

C  Linearizations of optical depth integrations
C  Linearized Green function multipliers

 	  DO AA = 1, NSTREAMS
	    ZDEL  =   T_DELT_EIGEN(AA,N)
	    ZWDEL = ZDEL * WDEL
	    DO Q = 1, N_PARAMETERS
	      L_WDEL  = L_T_DELT_MUBAR(N,N,Q)
	      L_ZDEL  = L_T_DELT_EIGEN(AA,N,Q)
	      L_ZWDEL = ZDEL * L_WDEL + L_ZDEL * WDEL
	      L_SD  =  ( L_ZDEL - L_WDEL ) / ( ZDEL - WDEL )
	      L_SU  =        - L_ZWDEL     / ( ONE - ZWDEL )
	      T1 = L_ATERM_SAVE(AA,N,Q) + L_GAMMA_M(AA,N,N,Q)
	      L_GFUNC_DN(AA,Q) = GFUNC_DN(AA,N) * ( T1 + L_SD )
	      T1 = L_BTERM_SAVE(AA,N,Q) + L_GAMMA_P(AA,N,N,Q)
	      L_GFUNC_UP(AA,Q) = GFUNC_UP(AA,N) * ( T1 + L_SU )
	    ENDDO
	  ENDDO

C  Set linearized form of particular integral at boundaries

	  DO Q = 1, N_PARAMETERS
	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      S_P_U = ZERO
	      S_P_L = ZERO
	      S_M_U = ZERO
	      S_M_L = ZERO
	      DO AA = 1, NSTREAMS
	        S_P_U = S_P_U + L_GFUNC_UP(AA,Q) *   XPOS(I1,AA,N) +
     &                          GFUNC_UP(AA,N)   * L_XPOS(I1,AA,N,Q)
	        S_M_U = S_M_U + L_GFUNC_UP(AA,Q) *   XPOS(I,AA,N) +
     &                          GFUNC_UP(AA,N)   * L_XPOS(I,AA,N,Q)
	        S_P_L = S_P_L + L_GFUNC_DN(AA,Q) *   XPOS(I,AA,N) +
     &                          GFUNC_DN(AA,N)   * L_XPOS(I,AA,N,Q)
	        S_M_L = S_M_L + L_GFUNC_DN(AA,Q) *   XPOS(I1,AA,N) +
     &                          GFUNC_DN(AA,N)   * L_XPOS(I1,AA,N,Q)
	      ENDDO
	      L_WUPPER(I,N,Q) = S_P_U
	      L_WUPPER(I1,N,Q) = S_M_U
	      L_WLOWER(I1,N,Q) = S_M_L
	      L_WLOWER(I,N,Q) = S_P_L
	    ENDDO
	  ENDDO

C  end clause for choice of solution method

	ENDIF

C  Finish

	END

C

	SUBROUTINE L_BEAMSOL_NONVARY ( N, K, K_PARAMETERS )

C  Linearization of Green's function beam particular integral, one layer only.
C  This is for a layer that is non-varying.

C  Include files
C  =============

C  input
C  -----

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_CONTROL.VARS'

C  include file of Legendre stuff (input to this module)

	INCLUDE '../include_s/LIDORT_LEGENDRE.VARS'

C  include file of Set up stuff 

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include files of main solution variables

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include file of Multiplier coefficient variables

	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  output
C  ------

C  include file of main solution linearization variables (output)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  include file of setup linearization variables (output)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'

C  include file of linearized Multiplier coefficients (output)

	INCLUDE '../include_e/LIDORT_L_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer indices (inputs)

	INTEGER		 N, K

C  number of varying parameters (input)

	INTEGER		 K_PARAMETERS

C  Local variables
C  ---------------

C  help variables

	INTEGER		 AA, I, I1, Q
	DOUBLE PRECISION S_P_U, S_P_L, S_M_U, S_M_L, L_SD, L_SU 
	DOUBLE PRECISION CONST, WDEL, ZDEL, ZWDEL, T1
	DOUBLE PRECISION PQ, L_WDEL, L_ZWDEL
	DOUBLE PRECISION TRANS2, VAR1, VAR2, VAR_U

C   ### NEW CODE ##################################################
C  No linearized particular solution beyond the cutoff layer.
C    [ Zero the boundary layer values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO Q = 1, K_PARAMETERS
	    DO I = 1, NSTR2
	      L_WUPPER(I,N,Q) = ZERO
	      L_WLOWER(I,N,Q) = ZERO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  Case 2. If N NOT equal to LAYER_TO_VARY
C  ---------------------------------------

C  initialise

	CONST   = INITIAL_TRANS(N)
	WDEL    = T_DELT_MUBAR(N)

C  Classical solution
C  ------------------

	IF ( DO_CLASSICAL_SOLUTION ) THEN

	  TRANS2 = CONST * WDEL

C  ..(a) quasi-spherical for n > 1 (only gets done for this case anyway)

	  IF ( DO_QUASPHER_BEAM ) THEN

	    DO Q = 1, K_PARAMETERS
	      VAR1 = L_T_DELT_MUBAR(N,K,Q) * CONST
	      VAR2 = L_INITIAL_TRANS(N,K,Q)
	      DO I = 1, NSTR2
	        VAR_U = VAR2 * WVEC(I,N) + L_WVEC(I,N,K,Q)
	        L_WUPPER(I,N,Q) = CONST  * VAR_U
	        L_WLOWER(I,N,Q) = TRANS2 * VAR_U + VAR1 * WVEC(I,N)
	      ENDDO
	    ENDDO

C  ..(b) plane-parallel and qs for n = 1

	  ELSE

	    DO Q = 1, K_PARAMETERS
	      VAR1 = L_INITIAL_TRANS(N,K,Q) * CONST
	      DO I = 1, NSTR2
	        VAR_U = VAR2 * WVEC(I,N)
	        L_WUPPER(I,N,Q) = VAR1 * WVEC(I,N)
	        L_WLOWER(I,N,Q) = L_WUPPER(I,N,Q) * WDEL
	      ENDDO
	    ENDDO

	  ENDIF

C  Green's function solution
C  -------------------------

	ELSE

C  ..(a) quasi-spherical for n > 1 (only gets done for this case anyway)

	  IF ( DO_QUASPHER_BEAM ) THEN

C  Set up linearizations of GAMMA constants

	    DO AA = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        PQ = L_AVERAGE_SECANT(N,K,Q)
	        L_GAMMA_P(AA,N,K,Q) = - GAMMA_P(AA,N) * PQ
	        L_GAMMA_M(AA,N,K,Q) = - GAMMA_M(AA,N) * PQ
	      ENDDO
	    ENDDO

C  Linearizations of optical depth integrations
C  Linearized Green function multipliers

	    DO AA = 1, NSTREAMS
	      ZDEL  =   T_DELT_EIGEN(AA,N)
	      ZWDEL = ZDEL * WDEL
	      DO Q = 1, K_PARAMETERS
	        L_WDEL  = L_T_DELT_MUBAR(N,K,Q)
	        L_ZWDEL = ZDEL * L_WDEL
	        L_SD  = - L_WDEL  / ( ZDEL - WDEL )
	        L_SU  = - L_ZWDEL / ( ONE - ZWDEL )
	        T1 = L_INITIAL_TRANS(N,K,Q) + L_GAMMA_M(AA,N,K,Q)
	        L_GFUNC_DN(AA,Q) = GFUNC_DN(AA,N) * ( T1 + L_SD )
	        T1 = L_INITIAL_TRANS(N,K,Q) + L_GAMMA_P(AA,N,K,Q)
	        L_GFUNC_UP(AA,Q) = GFUNC_UP(AA,N) * ( T1 + L_SU )
	      ENDDO
	    ENDDO

C  ..(b) plane-parallel and qs for n = 1

	  ELSE

	    DO AA = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        T1 = L_INITIAL_TRANS(N,K,Q)
	        L_GFUNC_DN(AA,Q) = GFUNC_DN(AA,N) * T1
	        L_GFUNC_UP(AA,Q) = GFUNC_UP(AA,N) * T1
	      ENDDO
	    ENDDO

	  ENDIF

C  Set linearized form of particular integral

	  DO Q = 1, K_PARAMETERS
	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      S_P_U = ZERO
	      S_P_L = ZERO
	      S_M_U = ZERO
	      S_M_L = ZERO
	      DO AA = 1, NSTREAMS
	        S_P_U = S_P_U + L_GFUNC_UP(AA,Q) *   XPOS(I1,AA,N)
	        S_M_U = S_M_U + L_GFUNC_UP(AA,Q) *   XPOS(I,AA,N)
	        S_P_L = S_P_L + L_GFUNC_DN(AA,Q) *   XPOS(I,AA,N)
	        S_M_L = S_M_L + L_GFUNC_DN(AA,Q) *   XPOS(I1,AA,N)
	      ENDDO
	      L_WUPPER(I,N,Q) = S_P_U
	      L_WUPPER(I1,N,Q) = S_M_U
	      L_WLOWER(I1,N,Q) = S_M_L
	      L_WLOWER(I,N,Q) = S_P_L
	    ENDDO
	  ENDDO

C  End clause for choice of method

	ENDIF

C  Finish

	END

C

	SUBROUTINE LIDORT_ALBEDOWF_COLUMN_SETUP
     I   (  DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS, DIRECT_BEAM )

C  Standard Include files (all inputs)
C  ----------------------

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  include files of Standard geophysical and model variables

	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of Standard solution and setup variables

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  Extended Include files
C  ----------------------

C  include files of linearized solution variables (output stored here)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  inputs
C  ------

C  .. direct beam inputs

	LOGICAL		 DO_INCLUDE_DIRECTBEAM
	DOUBLE PRECISION DIRECT_BEAM(MAXSTRM)

C  .. surface emission inputs

	LOGICAL		 DO_INCLUDE_SURFEMISS

C  local variables
C  ---------------

	INTEGER		 N, I, C0, CM, AA, N1
	DOUBLE PRECISION HELP, EMISS_VAR, EMISS_INTENSITY

C  initialise

	DO I = 1, MAXTOTAL
	  COL2_WFALB(I,1) = ZERO
	ENDDO

C  boundary conditions not changed for first layer upper (TOA)

	N = 1
	DO I = 1, NSTREAMS
	  COL2_WFALB(I,1) = ZERO
	ENDDO

C  boundary conditions not changed for all intermediate levels

	DO N = 2, NLAYERS - 1
	  N1 = N - 1
          C0  = N1*NSTR2 - NSTREAMS
	  DO I = 1, NSTR2
	    CM = C0 + I
	    COL2_WFALB(CM,1) = ZERO
	  ENDDO
	ENDDO

C  Ground level boundary condition is changed

	N = NLAYERS
	C0 = (N-1)*NSTR2 + NSTREAMS
	DO I = 1, NSTREAMS
	  CM = C0 + I
	  HELP = R2_BEAM(I)
	  DO AA = 1, NSTREAMS
	    HELP = HELP + LCON(AA,N)*T_DELT_EIGEN(AA,N)*R2_HOMP(I,AA)
     &                  + MCON(AA,N)*R2_HOMM(I,AA)
	  ENDDO
	  COL2_WFALB(CM,1) = HELP
	ENDDO

C  Add direct beam variation of albedo

	IF ( DO_INCLUDE_DIRECTBEAM ) THEN
	  DO I = 1, NSTREAMS
	    CM = C0 + I
	    COL2_WFALB(CM,1) = COL2_WFALB(CM,1) + DIRECT_BEAM(I)
	  ENDDO
	ENDIF

C  include emissivity variation (checked 26 May)
C  (expression for emissivity variation follows from Kirchhoff's law)

	IF ( DO_INCLUDE_SURFEMISS ) THEN
	  HELP = SURFBB
	  DO I = 1, NSTREAMS
	    CM = C0 + I
	    EMISS_INTENSITY = EMISSIVITY(I) * HELP
	    EMISS_VAR = EMISS_INTENSITY * (ONE - (ONE/EMISSIVITY(I)))
	    COL2_WFALB(CM,1) = COL2_WFALB(CM,1) + EMISS_VAR
	  ENDDO
	ENDIF

C  Finish

	END

C

	SUBROUTINE LIDORT_SURFBBWF_COLUMN_SETUP

C  Standard Include files (all inputs)
C  ----------------------

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  include files of Standard geophysical and model variables

	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  Extended Include files
C  ----------------------

C  include files of linearized solution variables (output stored here)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  local variables
C  ---------------

	INTEGER		 N, I, C0, CM, N1
	DOUBLE PRECISION HELP, EMISS_VAR

C  boundary conditions not changed for first layer upper (TOA)

	N = 1
	DO I = 1, NSTREAMS
	  COL2_WFSBB(I,1) = ZERO
	ENDDO

C  boundary conditions not changed for all intermediate levels

	DO N = 2, NLAYERS - 1
	  N1 = N - 1
          C0  = N1*NSTR2 - NSTREAMS
	  DO I = 1, NSTR2
	    CM = C0 + I
	    COL2_WFSBB(CM,1) = ZERO
	  ENDDO
	ENDDO

C  Ground level boundary condition is changed by surface emission variation

	N = NLAYERS
	C0 = (N-1)*NSTR2 + NSTREAMS
	HELP = SURFBB
	DO I = 1, NSTREAMS
	  CM = C0 + I
	  EMISS_VAR = EMISSIVITY(I) * HELP
	  COL2_WFSBB(CM,1) = EMISS_VAR
	ENDDO

C  Finish

	END
