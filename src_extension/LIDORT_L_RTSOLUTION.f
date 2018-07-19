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

	SUBROUTINE LIDORT_L_RTSOLUTION
     I    ( LAYER, FOURIER, DOVARY, N_PARAMETERS,
     O      STATUS )

C  Linearizations of Discrete Ordinate Differential Equation Solutions

C  Include files
C  =============

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model and control input variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_CONTROL.VARS'

C  include file of setup variables (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include file of main solution variables (input to this module)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  Subroutine input arguments

	INTEGER		 FOURIER, LAYER, N_PARAMETERS
	LOGICAL		 DOVARY

C  output status

	INTEGER		 STATUS

C  Local variables
C  ===============

C  local variables for error tracing

	INTEGER		 STATUS_SUB, I
	CHARACTER*(70)	 MAIL, TRACE
	CHARACTER*3	 C2

	DOUBLE PRECISION QSAVE(MAXSTRM)

C  Start of code
C  =============

C  initialise status

	STATUS = LIDORT_SUCCESS

C  linearize homogeneous solutions (flagged)
C  ===============================

	IF ( DOVARY ) THEN

	  CALL L_HOMOGENEOUS_SOLUTION
     &       ( LAYER, FOURIER, N_PARAMETERS, STATUS_SUB )

	  IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
	    WRITE(C2,'(I2)')LAYER
	    MAIL = 'Error from L_HOMOGENEOUS_SOLUTION, layer'//C2
	    TRACE= 'Call in LIDORT_L_RTSOLUTION'
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
	    RETURN
	  ENDIF

	ENDIF

C  Linearize Particular integrals
C  ==============================

C  1. For the classical solution ( two parts )
C  -----------------------------

	IF ( DO_CLASSICAL_SOLUTION ) THEN

C  set up driving vector (using saved results)
C  This must be done regardless of whether layer N is varying or not.

	  DO I = 1, NSTREAMS
	    QSAVE(I) = QDIFVEC_SAVE(I) +
     &                 TWO * AVERAGE_SECANT(LAYER) * QVEC_SAVE(I)
	  ENDDO

C  Main linearization for the present layer N

	  IF ( DOVARY ) THEN

	    CALL L_CLASSICAL_BEAM_SOLUTION_1
     &        ( LAYER, FOURIER, QSAVE, N_PARAMETERS, STATUS_SUB )

	    IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
	      WRITE(C2,'(I2)')LAYER
	      MAIL =
     &           'Error from L_CLASSICAL_BEAM_SOLUTION_1, layer'//C2
	      TRACE= 'Call in LIDORT_L_RTSOLUTION'
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
	      RETURN
	    ENDIF

	  ENDIF

C  additional variations in this layer caused by other variational input
C  (quasi-spherical only)

	  IF ( DO_QUASPHER_BEAM ) THEN

	    CALL L_CLASSICAL_BEAM_SOLUTION_2
     &          ( LAYER, QSAVE, STATUS_SUB )

	    IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
	      WRITE(C2,'(I2)')LAYER
	      MAIL =
     &           'Error from L_CLASSICAL_BEAM_SOLUTION_2, layer'//C2
	      TRACE= 'Call in LIDORT_L_RTSOLUTION'
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
	      RETURN
	    ENDIF

	  ENDIF

C  2. For the Green function solution
C  ----------------------------------

	ELSE

	  IF ( DOVARY ) THEN
	    CALL L_GREENFUNC_BEAM_SOLUTION
     &         ( LAYER, FOURIER, N_PARAMETERS )
	  ENDIF

	ENDIF

C  Finish

	END

C

	SUBROUTINE L_HOMOGENEOUS_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, N_PARAMETERS,
     O      STATUS )

C  Linearization of the homogeneous solutions.
C  The actual solutions have already been found.

C  Include files
C  =============

C  input
C  -----

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of Legendre stuff (local to this module)

	INCLUDE '../include_s/LIDORT_LEGENDRE.VARS'

C  include file of setup variables (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include file of main solution variables (input to this module)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  output
C  ------

C  include file of main solution linearization variables (output)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  include file of setup linearization variables (output)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

	INTEGER		 GIVEN_LAYER
	INTEGER		 FOURIER

C  number of varying parameters (input)

	INTEGER		 N_PARAMETERS

C  output status

	INTEGER		 STATUS

C  Local variables
C  ---------------

C  local matrices for eigenvalue computation

	DOUBLE PRECISION
     &      HMAT(MAXSTRM_P1,MAXSTRM_P1),HVEC(MAXSTRM_P1,MAX_PARAMETERS),
     &      L_DIFVEC(MAXSTRM), L_EIGENVEC(MAXSTRM)

C  extra output/input for linear-algebra solver (LAPACK)

	INTEGER		 IPIV(MAXSTRM_P1), INFO

C  Miscellaneous local variables

	INTEGER		 I, J, I1, L, N, M, AA, K, Q, NSTRM_P1
 	DOUBLE PRECISION DP, DM, SUM, L_KVAL, LTQ, KSQD,
     &                   KVAL, HELP, XINV, FAC, TBAS

C  Start of code
C  -------------

C  initialise status

	STATUS = LIDORT_SUCCESS

C  Layer and Fourier number

	N = GIVEN_LAYER
	M = FOURIER
	NSTRM_P1 = NSTREAMS + 1

C  Linearize the Eignematrix
C  -------------------------

C  set up CP and CM (linearizations of SAB and DAB)

	DO I = 1, NSTREAMS
	  XINV = ONE/X(I)
	  DO J = 1, NSTREAMS
	    FAC = HALFA(J) * XINV
	    DO Q = 1, N_PARAMETERS
	      DP = ZERO
	      DM = ZERO
	      DO L = M, NMOMENTS
	        DP = DP + PLMI_PLMJ_P(I,J,L) * L_OMEGA_MOMS(N,L,Q)
	        DM = DM + PLMI_PLMJ_M(I,J,L) * L_OMEGA_MOMS(N,L,Q)
	      ENDDO
              L_SAB(I,J,Q) = FAC * ( DP + DM )
              L_DAB(I,J,Q) = FAC * ( DP - DM )
	    ENDDO
	  ENDDO
	ENDDO

C  set up linearized eigenmatrices (one for each varying parameter)

	DO Q = 1, N_PARAMETERS
	  DO I = 1, NSTREAMS
	    DO J = 1, NSTREAMS
	      SUM = ZERO
	      DO K = 1, NSTREAMS
	        SUM = SUM + L_DAB(I,K,Q)*SAB(K,J)+DAB(I,K)*L_SAB(K,J,Q)
	      ENDDO
	      L_EIGENMAT(I,J,Q) = SUM
	    ENDDO
	  ENDDO
        ENDDO

C  Matrix and column setup
C  -----------------------

C  Do this for each eigenvactor

	DO AA = 1, NSTREAMS

	  KVAL = KEIGEN(AA,N)
	  KSQD = KVAL * KVAL

C  initialise solution matrix HMAT (important to do this!)

	  DO I = 1, MAXSTRM_P1
	    DO J = 1, MAXSTRM_P1
	      HMAT(I,J) = ZERO
	    ENDDO
	  ENDDO

C  Determine solution matrix HMAT (A in the matrix equation A.x = B)

	  DO I = 1, NSTREAMS
	    DO J = 1, NSTREAMS
	      HMAT(I,J+1) = - EIGENMAT_SAVE(I,J)
	    ENDDO
	    HMAT(I,I+1) = HMAT(I,I+1) + KSQD
	    HMAT(I,1)   = TWO * KVAL * EIGENVEC_SAVE(I,AA)
	    HMAT(NSTRM_P1,I+1) = EIGENVEC_SAVE(I,AA)
	  ENDDO
	  HMAT(NSTRM_P1,1) = ZERO

C  solution column vectors (B in the matrix equation A.x = B).
C    (One for each parameter to be varied)

	  DO Q = 1, N_PARAMETERS
	    DO I = 1, NSTREAMS
	      HELP = ZERO
	      DO J = 1, NSTREAMS
	        HELP = HELP + L_EIGENMAT(I,J,Q) * EIGENVEC_SAVE(J,AA)
	      ENDDO
	      HVEC(I,Q) = HELP
	    ENDDO
	    HVEC(NSTRM_P1,Q) = ZERO
	  ENDDO

C  Solve matrix system for linearization factors
C  ---------------------------------------------

C  solve for the linearized sum-vector + eigenvalue
C  July 1 1999, test using LAPACK modules DGETRF and DEGTRS successful

C   .. LU-decomposition of the matrix HMAT using DGETRF

	  CALL DGETRF (NSTRM_P1,NSTRM_P1,HMAT,MAXSTRM_P1,IPIV,INFO)

	  IF ( INFO .NE. 0 ) THEN
	    CALL LIDORT_LAPACK_ERROR
     I      ( INFO, N, ' Homogeneous Linearization DGETRF', STATUS )
	    IF ( STATUS .EQ. LIDORT_SERIOUS ) RETURN
	  ENDIF

C   .. Back substitution for column vectors (one for each parameter varying)

	  CALL DGETRS ('N',NSTRM_P1,N_PARAMETERS,HMAT,
     &                      MAXSTRM_P1,IPIV,HVEC,MAXSTRM_P1,INFO)

	  IF ( INFO .NE. 0 ) THEN
	    CALL LIDORT_LAPACK_ERROR
     I      ( INFO, N, ' Homogeneous Linearization DGETRS', STATUS )
	    IF ( STATUS .EQ. LIDORT_SERIOUS ) RETURN
	  ENDIF

C  Assign linearization factors for each scatterer
C  -----------------------------------------------

C  Start loop over varying parameters

	  DO Q = 1, N_PARAMETERS

C  assign linearization for eigenvalue K

	    L_KVAL = HVEC(1,Q)
	    L_KEIGEN(AA,N,Q) = L_KVAL

C  linearization of the actual eigenvector

	    DO I = 1, NSTREAMS
	      L_EIGENVEC(I) = HVEC(I+1,Q)
	    ENDDO

C  linearized difference vector

	    DO I = 1, NSTREAMS
	      SUM = ZERO
	      DO K = 1, NSTREAMS
	        SUM = SUM - L_SAB(I,K,Q) * EIGENVEC_SAVE(K,AA)
     &                    - SAB(I,K)     * L_EIGENVEC(K)
	      ENDDO
	      HELP = ( SUM - L_KVAL * DIFVEC_SAVE(I,AA) ) / KVAL
	      L_DIFVEC(I) = HELP
	    ENDDO

C  assign linearization for 'positive' homogeneous solution vectors

	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      L_XPOS(I1,AA,N,Q) = HALF*(L_EIGENVEC(I)-L_DIFVEC(I))
	      L_XPOS(I,AA,N,Q)  = HALF*(L_EIGENVEC(I)+L_DIFVEC(I))
	    ENDDO

C  symmetry for linearized 'negative' homogeneous solution vectors

	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      L_XNEG(I1,AA,N,Q) = L_XPOS(I,AA,N,Q)
	      L_XNEG(I,AA,N,Q)  = L_XPOS(I1,AA,N,Q)
	    ENDDO

C  End parameter loop

	  ENDDO

C  End eigenvalue loop

	ENDDO

C  transmittance factors
C  ---------------------

C  linearization of transmittance factors

	DO AA = 1, NSTREAMS
	  TBAS = - T_DELT_EIGEN(AA,N) * DELTA(N)
	  DO Q = 1, N_PARAMETERS
	    LTQ = L_KEIGEN(AA,N,Q) + KEIGEN(AA,N) * VQ(N,Q)
	    L_T_DELT_EIGEN(AA,N,Q) = TBAS * LTQ
	  ENDDO
	ENDDO

C  Finish

	END

C

	SUBROUTINE L_CLASSICAL_BEAM_SOLUTION_1
     I    ( GIVEN_LAYER, FOURIER, QSAVE, N_PARAMETERS,
     O      STATUS )

C  linearized values of the classical particular solution.

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

C  output
C  ------

C  include file of main solution linearization variables (output)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  include file of setup linearization variables (output)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

	INTEGER		 GIVEN_LAYER
	INTEGER		 FOURIER

C  number of varying parameters (input)

	INTEGER		 N_PARAMETERS

C  Saved vector

	DOUBLE PRECISION QSAVE(MAXSTRM)

C  output status

	INTEGER		 STATUS

C  Local variables
C  ---------------

	DOUBLE PRECISION L_QVEC(MAXSTRM,MAX_PARAMETERS)
	DOUBLE PRECISION L_QSUMVEC(MAXSTRM,MAX_PARAMETERS)
	DOUBLE PRECISION L_QDIFVEC(MAXSTRM,MAX_PARAMETERS)

C  help variables

	INTEGER		 I, J, I1, L, N, M, Q, INFO
 	DOUBLE PRECISION TP, TM, TM1, TM2, TM3, XINV,
     &                   HELP, WSUM, WDIF, SECBAR

C  Start of code
C  -------------

C  Set layer, fourier

	N = GIVEN_LAYER
	M = FOURIER

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the boundary layer values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO I = 1, NSTR2
	    DO Q = 1, N_PARAMETERS
	      L_WVEC(I,N,N,Q) = ZERO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  solar zenith cosine for this layer

	SECBAR = AVERAGE_SECANT(N)

C  For each varying parameter

	DO Q = 1, N_PARAMETERS

C  Set up linearized sum and difference vectors for Beam source terms

	  DO I = 1, NSTREAMS
	    XINV = ONE / X(I)
	    TP = ZERO
	    TM = ZERO
	    DO L = M, NMOMENTS
	      TP = TP + PLMI_X0_P(I,L,N) * L_OMEGA_MOMS(N,L,Q)
	      TM = TM + PLMI_X0_M(I,L,N) * L_OMEGA_MOMS(N,L,Q)
	    ENDDO
	    L_QSUMVEC(I,Q) =  ( TP + TM ) * XINV
	    L_QDIFVEC(I,Q) =  ( TP - TM ) * XINV
	  ENDDO

C   setup linearized RHS vector
C  ( use results from the original solution )

	  DO I = 1, NSTREAMS
	    HELP = ZERO
	    DO J = 1, NSTREAMS
	      TM1 = L_EIGENMAT(I,J,Q) * QVEC_SAVE(J)
	      TM2 = DAB(I,J)     * L_QSUMVEC(J,Q)
	      TM3 = L_DAB(I,J,Q) * QSUMVEC_SAVE(J)
	      HELP = HELP - TM1 - TM2 - TM3
	    ENDDO
	    L_QVEC(I,Q)  = HELP + L_QDIFVEC(I,Q) * SECBAR
	  ENDDO

	ENDDO

C  additional terms for the quasi-spherical case
C  ( layers greater than one )

	IF ( DO_QUASPHER_BEAM ) THEN
	  IF ( N.GT.1 ) THEN
	    DO Q = 1, N_PARAMETERS
	      DO I = 1, NSTREAMS
	        HELP = QSAVE(I) * L_AVERAGE_SECANT(N,N,Q)
	        L_QVEC(I,Q) = L_QVEC(I,Q) + HELP
	      ENDDO
	    ENDDO
	  ENDIF
	ENDIF

C  Solve problem by back substitution for all of the RHS vectors
C  ( uses L_U decomposition of QMAT_SAVE from original solution )

	CALL DGETRS
     &        ('N',NSTREAMS,N_PARAMETERS,QMAT_SAVE,
     &          MAXSTRM,QPIVOT,L_QVEC,MAXSTRM,INFO)

	IF ( INFO .NE. 0 ) THEN
	  CALL LIDORT_LAPACK_ERROR
     I ( INFO, N, ' Beam P.I. linearization (N,N) (DGETRS)', STATUS )
	  IF ( STATUS .EQ. LIDORT_SERIOUS ) RETURN
	ENDIF

C  assign solutions for the quasi-spherical case, N > 1

	IF ( DO_QUASPHER_BEAM .AND. N.GT.1 ) THEN

	  DO I = 1, NSTREAMS
	    TM3 = - QDIF_SAVE(I) / SECBAR
	    DO Q = 1, N_PARAMETERS
	      I1 = I + NSTREAMS
	      HELP = ZERO
	      DO J = 1, NSTREAMS
	        TM1 = SAB(I,J)     * L_QVEC(J,Q)
	        TM2 = L_SAB(I,J,Q) * QVEC_SAVE(J)
	        HELP = HELP - TM1 - TM2
	      ENDDO
	      WSUM = L_QVEC(I,Q)
	      TM2 = ( HELP - L_QSUMVEC(I,Q) ) / SECBAR
	      WDIF = L_AVERAGE_SECANT(N,N,Q) * TM3 + TM2
	      L_WVEC(I,N,N,Q)  = HALF * ( WSUM + WDIF )
	      L_WVEC(I1,N,N,Q) = HALF * ( WSUM - WDIF )
	    ENDDO
	  ENDDO

C  assign solutions for plane/parallel & quasi-spherical case N = 1

	ELSE

	  DO I = 1, NSTREAMS
	    DO Q = 1, N_PARAMETERS
	      I1 = I + NSTREAMS
	      HELP = ZERO
	      DO J = 1, NSTREAMS
	        TM1 = SAB(I,J)     * L_QVEC(J,Q)
	        TM2 = L_SAB(I,J,Q) * QVEC_SAVE(J)
	        HELP = HELP - TM1 - TM2
	      ENDDO
	      WSUM = L_QVEC(I,Q)
	      WDIF = ( HELP - L_QSUMVEC(I,Q) ) / SECBAR
	      L_WVEC(I,N,N,Q)  = HALF * ( WSUM + WDIF )
	      L_WVEC(I1,N,N,Q) = HALF * ( WSUM - WDIF )
	    ENDDO
	  ENDDO

	ENDIF

C  Finish

	END

C

	SUBROUTINE L_CLASSICAL_BEAM_SOLUTION_2
     I    ( GIVEN_LAYER, QSAVE, 
     O      STATUS )

C  linearized values of the classical particular solution.

C  Include files
C  =============

C  input
C  -----

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include file of input variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of Set up stuff 

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include files of main solution variables

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include files of linearized geophysical variables

	INCLUDE '../include_e/LIDORT_L_GEOPHYS.VARS'

C  include file of setup linearization variables

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'

C  output
C  ------

C  include file of main solution linearization variables

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and driving vector (inputs)

	INTEGER		 GIVEN_LAYER
	DOUBLE PRECISION QSAVE(MAXSTRM)

C  output status

	INTEGER		 STATUS

C  Local variables
C  ---------------

	INTEGER		 I, J, I1, N, Q, INFO, K, K_PARAMETERS
	DOUBLE PRECISION L_QVEC(MAXSTRM,MAX_PARAMETERS)
 	DOUBLE PRECISION TM3, HELP, WSUM, WDIF, SECBAR

C  Start of code
C  -------------

C  Set layer, fourier

	N = GIVEN_LAYER

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the boundary layer values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO K = 1, N - 1
	    K_PARAMETERS = LAYER_VARY_NUMBER(K)
	    DO I = 1, NSTR2
	      DO Q = 1, K_PARAMETERS
	        L_WVEC(I,N,K,Q) = ZERO
	      ENDDO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  solar zenith cosine for this layer

	SECBAR = AVERAGE_SECANT(N)

C  Loop over layers above N

	DO K = 1, N - 1

C  If there is a varying layer

	  IF ( LAYER_VARY_FLAG(K) ) THEN

C  number of varying parameters for this layer

	    K_PARAMETERS = LAYER_VARY_NUMBER(K)

C  Set up vector

	    DO I = 1, NSTREAMS
	      DO Q = 1, K_PARAMETERS
	        L_QVEC(I,Q) = QSAVE(I) * L_AVERAGE_SECANT(N,K,Q)
	      ENDDO
	    ENDDO

C  Solve problem by back substitution for all of the RHS vectors
C  ( uses L_U decomposition of QMAT_SAVE from original solution )

	    CALL DGETRS
     &        ('N',NSTREAMS,K_PARAMETERS,QMAT_SAVE,
     &          MAXSTRM,QPIVOT,L_QVEC,MAXSTRM,INFO)

	    IF ( INFO .NE. 0 ) THEN
	      CALL LIDORT_LAPACK_ERROR
     I ( INFO, N, ' Beam P.I. linearization (N,K) (DGETRS)', STATUS )
	      IF ( STATUS .EQ. LIDORT_SERIOUS ) RETURN
	    ENDIF

C  assign linearized solutions for layer N due to variations in layer K

	    DO I = 1, NSTREAMS
	      TM3 = - QDIF_SAVE(I) / SECBAR
	      DO Q = 1, K_PARAMETERS
	        I1 = I + NSTREAMS
	        HELP = ZERO
	        DO J = 1, NSTREAMS
	          HELP = HELP - SAB(I,J) * L_QVEC(J,Q)
	        ENDDO
	        WSUM = L_QVEC(I,Q)
	        WDIF = L_AVERAGE_SECANT(N,K,Q) * TM3 + (HELP/SECBAR)
	        L_WVEC(I,N,K,Q)  = HALF * ( WSUM + WDIF )
	        L_WVEC(I1,N,K,Q) = HALF * ( WSUM - WDIF )
	      ENDDO
	    ENDDO

C  end K-layer loop

	  ENDIF
	ENDDO

C  Finish

	END

C

	SUBROUTINE L_GREENFUNC_BEAM_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, N_PARAMETERS )

C  Green's function beam particular integral, one layer only.

C  Include files
C  =============

C  input
C  -----

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of Legendre stuff (input to this module)

	INCLUDE '../include_s/LIDORT_LEGENDRE.VARS'

C  include file of Set up stuff 

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include files of main solution variables

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include file of setup linearization variables (input)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'

C  output
C  ------

C  include file of main solution linearization variables (output)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

	INTEGER		 GIVEN_LAYER
	INTEGER		 FOURIER

C  number of varying parameters (input)

	INTEGER		 N_PARAMETERS

C  Local variables
C  ---------------

C  help variables

	INTEGER		 AA, J, J1, L, I, I1, M, N, Q
	DOUBLE PRECISION NORM, T1, T2
	DOUBLE PRECISION L_SUM_LA, L_SUM_LB, L_NORM, L_ATERM, L_BTERM
	DOUBLE PRECISION TA1, TA2, TB1, TB2, L_TPA, L_TMA

C  linearizations of componenet variables

	DOUBLE PRECISION
     &     L_DMI(MAXSTRM,MAX_PARAMETERS),
     &     L_DPI(MAXSTRM,MAX_PARAMETERS),
     &     L_NORM_SAVED(MAXSTRM,MAX_PARAMETERS)

C  initialise indices

	N = GIVEN_LAYER
	M = FOURIER

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the boundary layer values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO AA = 1, NSTREAMS
	    DO Q = 1, N_PARAMETERS
	      L_ATERM_SAVE(AA,N,Q) = ZERO
	      L_BTERM_SAVE(AA,N,Q) = ZERO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  Form quantities independent of optical depth
C  ============================================

C  set up linearizations of help arrays (independent of eigenvector)

	DO Q = 1, N_PARAMETERS
	  DO I = 1, NSTREAMS
	    L_TPA = ZERO
	    L_TMA = ZERO
	    DO L = M, NMOMENTS
	      L_TPA = L_TPA + PLMI_X0_P(I,L,N) * L_OMEGA_MOMS(N,L,Q)
	      L_TMA = L_TMA + PLMI_X0_M(I,L,N) * L_OMEGA_MOMS(N,L,Q)
	    ENDDO
	    L_DPI(I,Q) = L_TPA
	    L_DMI(I,Q) = L_TMA
	  ENDDO
	ENDDO

C  Set up linearized Norm

	DO Q = 1, N_PARAMETERS
	  DO AA = 1, NSTREAMS
	    NORM = ZERO
	    DO J = 1, NSTREAMS
	      J1 = J + NSTREAMS
	      T1 = XPOS(J,AA,N)  * L_XPOS(J,AA,N,Q)
	      T2 = XPOS(J1,AA,N) * L_XPOS(J1,AA,N,Q)
	      NORM = NORM + AX(J)*(T1-T2)
	    ENDDO
	    L_NORM_SAVED(AA,Q) = TWO * NORM
	  ENDDO
	ENDDO

C  linearize quantities independent of TAU (L_ATERM_SAVE, L_BTERM_SAVE)

	DO AA = 1, NSTREAMS
	  DO Q = 1, N_PARAMETERS
	    L_SUM_LA = ZERO
	    L_SUM_LB = ZERO
	    DO I = 1, NSTREAMS
	      I1 = I + NSTREAMS
	      TA1 = L_DMI(I,Q)*XPOS(I1,AA,N) + DMI(I)*L_XPOS(I1,AA,N,Q)
	      TA2 = L_DPI(I,Q)*XPOS(I,AA,N)  + DPI(I)*L_XPOS(I,AA,N,Q)
	      L_SUM_LA  = L_SUM_LA + A(I) * ( TA1 + TA2 )
	      TB1 = L_DMI(I,Q)*XPOS(I,AA,N)  + DMI(I)*L_XPOS(I,AA,N,Q)
	      TB2 = L_DPI(I,Q)*XPOS(I1,AA,N) + DPI(I)*L_XPOS(I1,AA,N,Q)
	      L_SUM_LB  = L_SUM_LB + A(I) * ( TB1 + TB2 )
	    ENDDO
	    L_NORM = L_NORM_SAVED(AA,Q)
	    L_ATERM = ( L_SUM_LA / ATERM_SAVE(AA,N) ) - L_NORM 
	    L_BTERM = ( L_SUM_LB / BTERM_SAVE(AA,N) ) - L_NORM 
	    L_ATERM_SAVE(AA,N,Q) = L_ATERM / NORM_SAVED(AA)
	    L_BTERM_SAVE(AA,N,Q) = L_BTERM / NORM_SAVED(AA)
	  ENDDO
	ENDDO

C  Finish

	END

C

	SUBROUTINE LIDORT_L_USERSOLUTION
     I    ( LAYER, FOURIER, DOVARY, N_PARAMETERS )

C  Linearized Solution components and Multipliers for user-defined output.

C  Include files
C  =============

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'

C  Subroutine input arguments

	INTEGER		 FOURIER, LAYER, N_PARAMETERS
	LOGICAL		 DOVARY

C  Solutions at stream angles
C  ==========================

C  general flag

	IF  ( STERM_LAYERMASK_UP(LAYER) .OR.
     &        STERM_LAYERMASK_DN(LAYER) ) THEN

C  Homogeneous

	  IF ( DOVARY ) THEN
	    CALL L_USER_HOMOG_SOLUTION ( LAYER, FOURIER, N_PARAMETERS )
	  ENDIF

C  Particular integrals for the varying layer

	  IF ( DOVARY ) THEN
	    IF ( DO_CLASSICAL_SOLUTION ) THEN
	      CALL L_USER_CLASSICAL_SOLUTION_1
     &           ( LAYER, FOURIER, N_PARAMETERS )
	    ELSE
	      CALL L_USER_GREENFUNC_SOLUTION
     &           ( LAYER, FOURIER, N_PARAMETERS )
	    ENDIF
	  ENDIF

C  Particular integrals for the nonvayring layers above (classical q-s)

	  IF ( DO_CLASSICAL_SOLUTION ) THEN
	    IF ( DO_QUASPHER_BEAM .AND. LAYER.GT.1 ) THEN
	      CALL L_USER_CLASSICAL_SOLUTION_2 ( LAYER, FOURIER )
	    ENDIF
	  ENDIF

	ENDIF

C  Finish

	END

C

	SUBROUTINE L_USER_HOMOG_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, N_PARAMETERS )

C  Include files
C  =============

C  Standard files (all inputs)
C  --------------------------

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of Legendre stuff (input to this module)

	INCLUDE '../include_s/LIDORT_LEGENDRE.VARS'

C  include files of setup and solution stuff (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  Extended files
C  --------------

C  include file of dimensions

	INCLUDE '../include_e/LIDORT_L.PARS'

C  include files of Linearized setup stuff (input to this module)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'

C  include files of Linearized solution stuff (output to this module)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index, Fourier index, parameter number (inputs)

	INTEGER		 GIVEN_LAYER
	INTEGER		 FOURIER
	INTEGER		 N_PARAMETERS

C  Local variables
C  ---------------

	INTEGER		 UM, J, J1, L, N, M, AA, Q
	DOUBLE PRECISION L_U_HELP_P(MAXSTRM,0:MAXMOMENT)
	DOUBLE PRECISION L_U_HELP_M(MAXSTRM,0:MAXMOMENT)
        DOUBLE PRECISION SUM_NEG, SUM_POS, POS1, POS2, NEG1, NEG2
        DOUBLE PRECISION ULP, L_ULP

C  Layer and Fourier

	N = GIVEN_LAYER
	M = FOURIER

C  Eigenvector interpolation to user-defined angles
C  ------------------------------------------------

C  For each eigenvector and parameter

	DO Q = 1, N_PARAMETERS
          DO AA = 1, NSTREAMS

C  For each moment, do inner sum over computational angles
C  for the positive and negative linearized eigenvectors

	    DO L = M, NMOMENTS
              SUM_POS = ZERO
              SUM_NEG = ZERO
              DO  J = 1, NSTREAMS
	        J1 = J + NSTREAMS
	        POS1 = L_XPOS(J1,AA,N,Q) * WT_LEGP(J,L)
 	        POS2 = L_XPOS(J,AA,N,Q)  * WT_LEGM(J,L)
	        NEG1 = L_XNEG(J1,AA,N,Q) * WT_LEGP(J,L)
	        NEG2 = L_XNEG(J,AA,N,Q)  * WT_LEGM(J,L)
                SUM_POS = SUM_POS + POS1 + POS2
                SUM_NEG = SUM_NEG + NEG1 + NEG2
	      ENDDO
	      L_U_HELP_P(AA,L) = SUM_POS
 	      L_U_HELP_M(AA,L) = SUM_NEG
            ENDDO

C  Now sum over all harmonic contributions at each user-defined stream

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      SUM_POS = ZERO
	      SUM_NEG = ZERO
              DO L = M, NMOMENTS
	        ULP   = U_LEG_P(UM,L) * OMEGA_MOMS(N,L)
 	        L_ULP = U_LEG_P(UM,L) * L_OMEGA_MOMS(N,L,Q)
                SUM_POS = SUM_POS + L_U_HELP_P(AA,L) *   ULP +
     &                                U_HELP_P(AA,L) * L_ULP
                SUM_NEG = SUM_NEG + L_U_HELP_M(AA,L) *   ULP +
     &                                U_HELP_M(AA,L) * L_ULP
              ENDDO
	      L_U_XPOS(UM,AA,N,Q) = SUM_POS
	      L_U_XNEG(UM,AA,N,Q) = SUM_NEG
	    ENDDO

C  end eigenvector and parameter loop

	  ENDDO
	ENDDO

C  Finish

	END

C

	SUBROUTINE L_USER_CLASSICAL_SOLUTION_1
     I    ( GIVEN_LAYER, FOURIER, N_PARAMETERS )

C  Include files
C  =============

C  Standard files (all inputs)
C  --------------------------

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of model and control variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of Legendre stuff (input to this module)

	INCLUDE '../include_s/LIDORT_LEGENDRE.VARS'

C  include files of setup and solution stuff (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  Extended files
C  --------------

C  include file of dimensions

	INCLUDE '../include_e/LIDORT_L.PARS'

C  include files of Linearized setup stuff (input to this module)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'

C  include files of Linearized solution stuff (output to this module)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index, Fourier index, parameter number (inputs)

	INTEGER		 GIVEN_LAYER
	INTEGER		 FOURIER
	INTEGER		 N_PARAMETERS

C  Local variables
C  ---------------

	INTEGER		 UM, J, J1, L, N, M, Q
 	DOUBLE PRECISION HELP(0:MAXMOMENT), SUM, POS1, POS2
 	DOUBLE PRECISION HELP_FOR(0:MAXMOMENT)

C  Layer and Fourier

	N = GIVEN_LAYER
	M = FOURIER

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the boundary layer values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO Q = 1, N_PARAMETERS
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      L_U_WFORPOS(UM,N,Q) = ZERO
	      L_U_WFORNEG(UM,N,Q) = ZERO
	      L_U_WPOS(UM,N,N,Q) = ZERO
	      L_U_WNEG(UM,N,N,Q) = ZERO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  start parameter loop

	DO Q = 1, N_PARAMETERS

C  For each moment do inner sum over computational angles
C  taking care to use chain rule on the linearization

   	  DO L = M, NMOMENTS
            SUM = ZERO
            DO  J = 1, NSTREAMS
	      J1 = J + NSTREAMS
	      POS1 = L_WVEC(J1,N,N,Q) * WT_LEGP(J,L)
 	      POS2 = L_WVEC(J,N,N,Q)  * WT_LEGM(J,L)
              SUM = SUM + POS1 + POS2
	    ENDDO
	    HELP_FOR(L) = LEG0_M(L,N) * L_OMEGA_MOMS(N,L,Q)
Cold	    HELP(L) = ( W_HELP(L) + LEG0_M(L,N) ) * L_OMEGA_MOMS(N,L,Q)
Cold     &                          + SUM         *   OMEGA_MOMS(N,L)
	    HELP(L) =       SUM *   OMEGA_MOMS(N,L) +
     &                W_HELP(L) * L_OMEGA_MOMS(N,L,Q)
         ENDDO

C  Now sum over all harmonic contributions at each user-defined stream
C  Distinguish between upwelling and downwelling

	  IF ( DO_UPWELLING ) THEN
	    IF ( STERM_LAYERMASK_UP(N) ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
Cold            POS1 = ZERO
Cold            DO L = M, NMOMENTS
Cold              POS1 = POS1 + HELP(L)*U_LEG_P(UM,L)
Cold            ENDDO
Cold	        L_U_WPOS(UM,N,N,Q) = POS1
                POS1 = ZERO
                DO L = M, NMOMENTS
                  POS1 = POS1 + HELP_FOR(L)*U_LEG_P(UM,L)
                ENDDO
	        L_U_WFORPOS(UM,N,Q) = POS1
                POS1 = ZERO
                DO L = M, NMOMENTS
                  POS1 = POS1 + HELP(L)*U_LEG_P(UM,L)
                ENDDO
	        L_U_WPOS(UM,N,N,Q) = POS1
	      ENDDO
	    ENDIF
	  ENDIF

	  IF ( DO_DNWELLING ) THEN
	    IF ( STERM_LAYERMASK_DN(N) ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
Cold            POS2 = ZERO
Cold            DO L = M, NMOMENTS
Cold              POS2 = POS2 + HELP(L)*U_LEG_M(UM,L)
Cold            ENDDO
Cold            L_U_WNEG(UM,N,N,Q) = POS2
                POS2 = ZERO
                DO L = M, NMOMENTS
                  POS2 = POS2 + HELP_FOR(L)*U_LEG_M(UM,L)
                ENDDO
	        L_U_WFORNEG(UM,N,Q) = POS2
                POS2 = ZERO
                DO L = M, NMOMENTS
                  POS2 = POS2 + HELP(L)*U_LEG_M(UM,L)
                ENDDO
	        L_U_WNEG(UM,N,N,Q) =  POS2
	      ENDDO
	    ENDIF
	  ENDIF

C  end parameter loop

	ENDDO

C  Finish

	END

C

	SUBROUTINE L_USER_CLASSICAL_SOLUTION_2
     I    ( GIVEN_LAYER, FOURIER )

C  Include files
C  =============

C  Standard files (all inputs)
C  --------------------------

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include files of model and control variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of Legendre stuff (input to this module)

	INCLUDE '../include_s/LIDORT_LEGENDRE.VARS'

C  include files of setup stuff (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  Extended files
C  --------------

C  include files of Linearized geophysical stuff (input to this module)

	INCLUDE '../include_e/LIDORT_L_GEOPHYS.VARS'

C  include files of Linearized setup stuff (input to this module)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'

C  include files of Linearized solution stuff (output to this module)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index, Fourier index (inputs)

	INTEGER		 GIVEN_LAYER
	INTEGER		 FOURIER

C  Local variables
C  ---------------

	INTEGER		 UM, J, J1, L, N, M, Q, K, K_PARAMETERS
 	DOUBLE PRECISION HELP(0:MAXMOMENT), SUM, POS1, POS2

C  Layer and Fourier

	N = GIVEN_LAYER
	M = FOURIER

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the boundary layer values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO K = 1, N - 1
	    K_PARAMETERS = LAYER_VARY_NUMBER(K)
	    DO Q = 1, K_PARAMETERS
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        L_U_WPOS(UM,N,K,Q) = ZERO
	        L_U_WNEG(UM,N,K,Q) = ZERO
	      ENDDO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  start loop over all layers above N

	DO K = 1, N - 1

C  only do if layer K has some variation

	  IF ( LAYER_VARY_FLAG(K) ) THEN

	    K_PARAMETERS = LAYER_VARY_NUMBER(K)

C  start parameter loop

	    DO Q = 1, K_PARAMETERS

C  For each moment do inner sum over computational angles
C  taking care to use chain rule on the linearization

   	      DO L = M, NMOMENTS
                SUM = ZERO
                DO  J = 1, NSTREAMS
	          J1 = J + NSTREAMS
	          POS1 = L_WVEC(J1,N,K,Q) * WT_LEGP(J,L)
 	          POS2 = L_WVEC(J,N,K,Q)  * WT_LEGM(J,L)
                  SUM = SUM + POS1 + POS2
	        ENDDO
	        HELP(L) =  SUM * OMEGA_MOMS(N,L)
              ENDDO

C  Now sum over all harmonic contributions at each user-defined stream
C  Distinguish between upwelling and downwelling

	      IF ( DO_UPWELLING ) THEN
	        IF ( STERM_LAYERMASK_UP(N) ) THEN
	          DO UM = LOCAL_UM_START, N_USER_STREAMS
                    POS1 = ZERO
                    DO L = M, NMOMENTS
                      POS1 = POS1 + HELP(L)*U_LEG_P(UM,L)
                    ENDDO
	            L_U_WPOS(UM,N,K,Q) = POS1
	          ENDDO
	        ENDIF
	      ENDIF

	      IF ( DO_DNWELLING ) THEN
	        IF ( STERM_LAYERMASK_DN(N) ) THEN
	          DO UM = LOCAL_UM_START, N_USER_STREAMS
                    POS2 = ZERO
                    DO L = M, NMOMENTS
                      POS2 = POS2 + HELP(L)*U_LEG_M(UM,L)
                    ENDDO
	            L_U_WNEG(UM,N,K,Q) = POS2
	          ENDDO
	        ENDIF
	      ENDIF

C  end parameter loop

	    ENDDO

C  end K-layer loop

	  ENDIF
	ENDDO

C  Finish

	END

C

	SUBROUTINE L_USER_GREENFUNC_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, N_PARAMETERS )

C  Include files
C  =============

C  Standard files (all inputs)
C  --------------------------

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of Legendre stuff (input to this module)

	INCLUDE '../include_s/LIDORT_LEGENDRE.VARS'

C  Extended files
C  --------------

C  include file of dimensions

	INCLUDE '../include_e/LIDORT_L.PARS'

C  include files of Linearized setup stuff (input to this module)

	INCLUDE '../include_e/LIDORT_L_SETUPS.VARS'

C  include files of Linearized solution stuff (output to this module)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index, Fourier index, parameter number (inputs)

	INTEGER		 GIVEN_LAYER
	INTEGER		 FOURIER
	INTEGER		 N_PARAMETERS

C  Local variables
C  ---------------

	INTEGER		 UM, L, N, M, Q
 	DOUBLE PRECISION HELP(0:MAXMOMENT), POS1, POS2

C  Layer and Fourier

	N = GIVEN_LAYER
	M = FOURIER

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the boundary layer values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO Q = 1, N_PARAMETERS
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      L_U_WFORPOS(UM,N,Q) = ZERO
	      L_U_WFORNEG(UM,N,Q) = ZERO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  start parameter loop

	DO Q = 1, N_PARAMETERS

C  For each moment do inner sum over computational angles

   	  DO L = M, NMOMENTS
	    HELP(L) = LEG0_M(L,N)*L_OMEGA_MOMS(N,L,Q)
          ENDDO

C  Now sum over all harmonic contributions at each user-defined stream

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
            POS1 = ZERO
            POS2 = ZERO
            DO L = M, NMOMENTS
              POS1 = POS1 + HELP(L)*U_LEG_P(UM,L)
              POS2 = POS2 + HELP(L)*U_LEG_M(UM,L)
            ENDDO
Cold	    L_U_WPOS(UM,N,N,Q) = POS1
Cold	    L_U_WNEG(UM,N,N,Q) = POS2
	    L_U_WFORPOS(UM,N,Q) = POS1
	    L_U_WFORNEG(UM,N,Q) = POS2
	  ENDDO

C  end parameter loop

	ENDDO

C  Finish

	END
