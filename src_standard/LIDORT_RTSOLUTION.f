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

	SUBROUTINE LIDORT_RTSOLUTION
     I    ( LAYER, FOURIER,
     O      STATUS )

C  Discrete Ordinate Differential Equation Solutions

C  Include files
C  =============

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'

C  Subroutine input arguments

	INTEGER		 FOURIER, LAYER
C	LOGICAL		 DO_INCLUDE_THERMAL		! REMOVED

C  output status

	INTEGER		 STATUS

C  Local variables
C  ===============

C  local variables for error tracing

	INTEGER		 STATUS_SUB
	CHARACTER*(70)	 MAIL, TRACE
	CHARACTER*3	 C2

C  Start of code
C  =============

C  initialise status

	STATUS = LIDORT_SUCCESS

C  Homogeneous solution
C  ====================

C  eigenvalue problem

	CALL HOMOGENEOUS_SOLUTION ( LAYER, FOURIER, STATUS_SUB )

C  error flag

	IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
	  WRITE(C2,'(I2)')LAYER
	  MAIL = 'Error from HOMOGENEOUS_SOLUTION, layer'//C2
	  TRACE= 'Call in LIDORT_RTSOLUTION'
	  STATUS = LIDORT_SERIOUS
	  CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
	  RETURN
	ENDIF

C  Particular integrals
C  ====================

C  1. For the classical solutions
C  ------------------------------

	IF ( DO_CLASSICAL_SOLUTION ) THEN

C  beam solution

	  CALL CLASSICAL_BEAM_SOLUTION ( LAYER, FOURIER, STATUS_SUB )

C  error flag

	  IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
	    WRITE(C2,'(I2)')LAYER
	    MAIL = 'Error from CLASSICAL_BEAM_SOLUTION, layer'//C2
	    TRACE= 'Call in LIDORT_RTSOLUTION'
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
	    RETURN
	  ENDIF

C  thermal solution if flagged (REMOVED)
C	  IF ( DO_INCLUDE_THERMAL ) THEN
C  Disabled in this version
C	    CALL CLASSICAL_THERMAL_SOLUTION
C	  ENDIF

C  2. For the Green function solution
C  ----------------------------------

	ELSE

C  beam solution

	  CALL GREENFUNC_BEAM_SOLUTION (LAYER, FOURIER)

C  thermal solution if flagged (REMOVED)
C	  IF ( DO_INCLUDE_THERMAL ) THEN
C  Placeholder Placeholder Placeholder Placeholder Placeholder
C	    CALL GREENFUNC_THERMAL_SOLUTION
C  Placeholder Placeholder Placeholder Placeholder Placeholder
C	  ENDIF

	ENDIF


C  Finish

	END

C

	SUBROUTINE HOMOGENEOUS_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, 
     O      STATUS )

C  Numerical solution of Eigenproblem.

C  Include files
C  =============

C  input
C  -----

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of Legendre stuff (input to this module)

	INCLUDE '../include_s/LIDORT_LEGENDRE.VARS'

C  output
C  ------

C  include file of Set up stuff (output to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include files of main solution variables (output to this module)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

	INTEGER		 GIVEN_LAYER
	INTEGER		 FOURIER

C  output status

	INTEGER		 STATUS

C  Local variables
C  ---------------

C  local matrices for eigenvalue computation

	DOUBLE PRECISION
     &      EIGENMAT(MAXSTRM,MAXSTRM)

C  (output from Eigenpackage module ASYMTX)

	DOUBLE PRECISION KSQ(MAXSTRM), WK(MAXSTRM2)
	DOUBLE PRECISION EVEC(MAXSTRM,MAXSTRM)
	INTEGER		 IER
	LOGICAL		 ASYMTX_FAILURE
	CHARACTER*70	 MESSAGE

C  output/input for Eigenpackage DGEEV ( LAPACK)
C	DOUBLE PRECISION REAL_KSQ(MAXSTRM),IMAG_KSQ(MAXSTRM)
C	DOUBLE PRECISION LEFT_EVEC(MAXSTRM,MAXSTRM)
C	DOUBLE PRECISION RITE_EVEC(MAXSTRM,MAXSTRM)
C	INTEGER		 LWORK
C	PARAMETER	 (LWORK = 4*MAXSTRM)
C	DOUBLE PRECISION WORK(LWORK)

C  Miscellaneous local variables

	INTEGER		 I, J, I1, L, N, M, AA, K
 	DOUBLE PRECISION DP, DM, SUM, FAC, KVAL, NORM, XINV, HELP
    
C  local variables for error tracing

	CHARACTER*(70)	 MAIL, TRACE
	CHARACTER*3	 CN, CI

C  initialise status

	STATUS = LIDORT_SUCCESS

C  Layer and Fourier

	N = GIVEN_LAYER
	M = FOURIER

C  Construct Eigenmatrix
C  ---------------------

C  zero the Eigenmatrix

	DO I = 1, NSTREAMS
	  DO J = 1, NSTREAMS
	    EIGENMAT(I,J) = ZERO
	  ENDDO
	ENDDO

C  Develop Sum and Difference matrices

	DO I = 1, NSTREAMS
	  XINV = ONE/X(I)
	  DO J = 1, NSTREAMS
	    FAC = XINV * HALFA(J)
	    DP = ZERO
	    DM = ZERO
	    DO L = M, NMOMENTS
	      DP = DP + PLMI_PLMJ_P(I,J,L) * OMEGA_MOMS(N,L)
	      DM = DM + PLMI_PLMJ_M(I,J,L) * OMEGA_MOMS(N,L)
	    ENDDO
	    SAB(I,J) = FAC * ( DP + DM ) 
	    DAB(I,J) = FAC * ( DP - DM )
	  ENDDO
	  SAB(I,I) = SAB(I,I) - XINV
	  DAB(I,I) = DAB(I,I) - XINV
	ENDDO

C  Compute Eigenmatrix

	DO I = 1, NSTREAMS 
	  DO J = 1, NSTREAMS
	    SUM = ZERO
	    DO K = 1, NSTREAMS
	      SUM = SUM + DAB(I,K) * SAB(K,J)
	    ENDDO
	    EIGENMAT(I,J) = SUM
	  ENDDO
	ENDDO

C  save Eigenmatrix (original is destroyed by ASMTYX)

	DO I = 1, NSTREAMS 
	  DO J = 1, NSTREAMS
	    EIGENMAT_SAVE(I,J) = EIGENMAT(I,J)
	  ENDDO
	ENDDO

C  Eigensolution package
C  ---------------------

C  Let's see how we get on with the DISORT package
C  tested 19 May 1999. Gives same results as Analytical approach.

C  second test using DGEEV (LAPACK module) - Worked against ASYMTX.
C  However, DGEEV is 1.7 - 2.0 times slower because it must look for
C  complex roots. [ASYMTX was specially written to avoid this search].
C  Also first call to DGEEV is 20 times slower. Tested June 24, 1999.
C  Conclusion stick with ASYMTX for now.

	CALL  ASYMTX
     I         ( EIGENMAT, NSTREAMS, MAXSTRM, MAXSTRM,
     O           EVEC, KSQ, IER, WK,
     O           MESSAGE, ASYMTX_FAILURE )

C  error tracing

	IF ( ASYMTX_FAILURE .OR. IER.GT.0 ) THEN
	  WRITE(CI,'(I3)')IER
	  WRITE(CN,'(I3)')N
	  MAIL  = 'eigenvalue '//CI//' has not converged'
	  IF ( ASYMTX_FAILURE ) MAIL = MESSAGE
	  TRACE ='ASYMTX error, HOMOGENEOUS_SOLUTION, Layer='//CN
	  STATUS = LIDORT_SERIOUS
	  CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
	  RETURN
	ENDIF

C  second test using DGEEV (LAPACK module). Here is what to use.
C	  CALL DGEEV
C     I      ( 'N', 'V', NSTREAMS, EIGENMAT,
C     O        MAXSTRM, REAL_KSQ, IMAG_KSQ,
C     O        LEFT_EVEC, MAXSTRM, EVEC, MAXSTRM,
C     W        WORK, LWORK, LAPACK_INFO )

C  Find solution eigenvectors XPOS, XNEG for all eigenvalues
C  ---------------------------------------------------------

	DO AA = 1, NSTREAMS

C  Store positive values in output array for each layer

	  KVAL = DSQRT(KSQ(AA))
	  KEIGEN(AA,N) = KVAL

C  Normalize eigenvectors to 1

	  NORM = ZERO
	  DO I = 1, NSTREAMS
	    NORM = NORM + EVEC(I,AA)*EVEC(I,AA)
	  ENDDO
	  NORM = DSQRT(NORM)

C  Find normalized eigenvector EIGENVEC_SAVE

	  DO I = 1, NSTREAMS
	    EIGENVEC_SAVE(I,AA) = EVEC(I,AA)/NORM
	  ENDDO

C  Find difference eigenvector DIFVEC (Siewert's notation)

	  DO I = 1, NSTREAMS
	    SUM = ZERO
	    DO K = 1, NSTREAMS
	      SUM = SUM - SAB(I,K) * EIGENVEC_SAVE(K,AA)
	    ENDDO
	    DIFVEC_SAVE(I,AA) = SUM / KVAL
          ENDDO

C  assign original evectors; first N are "DOWN", last N are "UP" (streams)

	  DO I = 1, NSTREAMS
	    I1 = I + NSTREAMS
	    XPOS(I,AA,N)  =
     &            HALF * ( EIGENVEC_SAVE(I,AA) + DIFVEC_SAVE(I,AA) )
	    XPOS(I1,AA,N) =
     &            HALF * ( EIGENVEC_SAVE(I,AA) - DIFVEC_SAVE(I,AA) )
	  ENDDO

C  Use symmetry properties to set -ve eigenvectors

	  DO I = 1, NSTREAMS
	    I1 = I + NSTREAMS
	    XNEG(I1,AA,N) = XPOS(I,AA,N)
	    XNEG(I,AA,N)  = XPOS(I1,AA,N)
	  ENDDO

C  End eignestream loop

	ENDDO

C  transmittance factors
C  ---------------------

C  Eigenstream transmittance factors for this layer

	DO AA = 1, NSTREAMS
	  HELP = KEIGEN(AA,N)*DELTA(N)
	  IF ( HELP .GT. MAX_TAU_QPATH ) THEN
	    T_DELT_EIGEN(AA,N) = ZERO
	  ELSE
	    T_DELT_EIGEN(AA,N) = DEXP(-HELP)
	  ENDIF
	ENDDO

C  Finish

	END

C

	SUBROUTINE CLASSICAL_BEAM_SOLUTION
     I    ( GIVEN_LAYER, FOURIER, 
     O      STATUS )

C  This is the classical Chandrasekhar beam solution.
C  ( plane parallel or average secant only)
C  Linear Matrix algebra.

C  Include files
C  =============

C  input
C  -----

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of Legendre stuff (input to this module)

	INCLUDE '../include_s/LIDORT_LEGENDRE.VARS'

C  output
C  ------
C  include file of Set up stuff (output to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include files of main solution variables (output to this module)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

	INTEGER		 GIVEN_LAYER
	INTEGER		 FOURIER

C  output status

	INTEGER		 STATUS

C  Local variables
C  ---------------

C  help variables

	CHARACTER*70	 MAIL
	INTEGER		 I, J, I1, L, N, M, INFO
 	DOUBLE PRECISION TP, TM, INV_X0SQ, SECBAR, XINV,
     &                   HELP, WSUM, WDIF, TRANS1, TRANS2

C  Layer

	M = FOURIER
	N = GIVEN_LAYER

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the boundary layer values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO I = 1, NSTR2
	    WUPPER(I,N) = ZERO
	    WLOWER(I,N) = ZERO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  set local values

	SECBAR = AVERAGE_SECANT(N)
	INV_X0SQ = SECBAR * SECBAR

C  Initialise matrix and column vector for solution

	DO I = 1, MAXSTRM
	  DO J = 1, MAXSTRM
	    QMAT_SAVE(I,J) = ZERO
	  ENDDO
	  QVEC_SAVE(I) = ZERO
	ENDDO

C  Set up sum and difference vectors for Beam source terms
C  ( sum vector may be required again in linearization )

	DO I = 1, NSTREAMS
	  XINV = ONE / X(I)
	  TP = ZERO
	  TM = ZERO
	  DO L = M, NMOMENTS
	    TP = TP + PLMI_X0_P(I,L,N) * OMEGA_MOMS(N,L)
	    TM = TM + PLMI_X0_M(I,L,N) * OMEGA_MOMS(N,L)
	  ENDDO
	  QSUMVEC_SAVE(I) =  ( TP + TM ) * XINV
	  QDIFVEC_SAVE(I) =  ( TP - TM ) * XINV
	ENDDO

C  solution matrix for the reduced problem
C  ( matrix should be saved in the LU decomposition form)

	DO I = 1, NSTREAMS
	  DO J = 1, NSTREAMS
	    QMAT_SAVE(I,J) = EIGENMAT_SAVE(I,J)
	  ENDDO
	  QMAT_SAVE(I,I) = QMAT_SAVE(I,I) - INV_X0SQ
	ENDDO

C  RHS vector for the reduced problem
C  ( this vector will be the answer after the linear algebra solution,
C    and may be needed again if there is linearization )

	DO I = 1, NSTREAMS
	  HELP = ZERO
	  DO J = 1, NSTREAMS
	    HELP = HELP - DAB(I,J)*QSUMVEC_SAVE(J)
	  ENDDO
	  QVEC_SAVE(I) = HELP + QDIFVEC_SAVE(I) * SECBAR
	ENDDO

C  L-U decomposition of the solution matrix

	CALL DGETRF(NSTREAMS,NSTREAMS,QMAT_SAVE,MAXSTRM,QPIVOT,INFO)

	IF ( INFO .NE. 0 ) THEN
	  MAIL='Classical Beam solution LU decomposition (DGETRF)'
	  CALL LIDORT_LAPACK_ERROR ( INFO, N, MAIL, STATUS )
	  IF ( STATUS .EQ. LIDORT_SERIOUS ) RETURN
	ENDIF

C  Solution of reduced problem by back-substitution

	CALL DGETRS
     &    ('N',NSTREAMS,1,QMAT_SAVE,MAXSTRM,QPIVOT,
     &     QVEC_SAVE,MAXSTRM,INFO)

	IF ( INFO .NE. 0 ) THEN
	  MAIL='Classical Beam solution back-substitution (DGETRS)'
	  CALL LIDORT_LAPACK_ERROR ( INFO, N, MAIL, STATUS )
	  IF ( STATUS .EQ. LIDORT_SERIOUS ) RETURN
	ENDIF

C  Assigning beam particular integral solution vector

	DO I = 1, NSTREAMS
	  I1 = I + NSTREAMS
	  HELP = ZERO
	  DO J = 1, NSTREAMS
	    HELP = HELP - SAB(I,J)*QVEC_SAVE(J)
	  ENDDO
	  QDIF_SAVE(I) = ( HELP - QSUMVEC_SAVE(I) ) / SECBAR
	  WSUM = QVEC_SAVE(I)
	  WDIF = QDIF_SAVE(I)
	  WVEC(I,N)  = HALF * ( WSUM + WDIF )
	  WVEC(I1,N) = HALF * ( WSUM - WDIF )
	ENDDO

C  Values at the layer boundaries
C  (transmittance factors have been determined in SETUPS module)

	TRANS1 = INITIAL_TRANS(N)
	TRANS2 = T_DELT_MUBAR(N) * TRANS1
	DO I = 1, NSTR2
	  WUPPER(I,N) = WVEC(I,N)*TRANS1
	  WLOWER(I,N) = WVEC(I,N)*TRANS2
	ENDDO

C  Finish

	END

C

	SUBROUTINE GREENFUNC_BEAM_SOLUTION
     I    ( GIVEN_LAYER, FOURIER )

C  Green's function beam particular integral, one layer only.
C  Uses coefficient expansion of attenuation.

C  Include files
C  =============

C  input
C  -----

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model variables (input)

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of Legendre stuff (input)

	INCLUDE '../include_s/LIDORT_LEGENDRE.VARS'

C  output
C  ------

C  include file of setup stuff

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include file of solution variables

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include file of Multiplier coefficients

	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  ====================

C  Given layer index and Fourier number (inputs)

	INTEGER		 GIVEN_LAYER
	INTEGER		 FOURIER

C  local variables
C  ===============

	INTEGER		 AA, J, J1, L, I, I1, M, N
	DOUBLE PRECISION NORM, T1, T2
	DOUBLE PRECISION TPA, TMA, SUM_LA, SUM_LB
	DOUBLE PRECISION S_P_U, S_P_L, S_M_U, S_M_L

	DOUBLE PRECISION RHO_M, RHO_P
	DOUBLE PRECISION CONST, SECBAR, WDEL, ZDEL, ZWDEL

C  initialise indices

	N = GIVEN_LAYER
	M = FOURIER

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the boundary layer values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO I = 1, NSTR2
	    WUPPER(I,N) = ZERO
	    WLOWER(I,N) = ZERO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

	SECBAR = AVERAGE_SECANT(N)
	CONST  = INITIAL_TRANS(N)

C  Gamma constants

	DO AA = 1, NSTREAMS
	  RHO_M = SECBAR - KEIGEN(AA,N)
	  RHO_P = SECBAR + KEIGEN(AA,N)
	  GAMMA_P(AA,N) = ONE / RHO_P
	  GAMMA_M(AA,N) = ONE / RHO_M
	ENDDO

C  3. Optical depth integrations for the discrete ordinate solution
C     =============================================================

	WDEL    = T_DELT_MUBAR(N)
	DO AA = 1, NSTREAMS
	  ZDEL  = T_DELT_EIGEN(AA,N)
	  ZWDEL = ZDEL * WDEL
	  CFUNC(AA,N)  = ( ZDEL - WDEL ) * GAMMA_M(AA,N)
	  DFUNC(AA,N)  = ( ONE - ZWDEL ) * GAMMA_P(AA,N)
	ENDDO

C  4. Form quantities independent of optical depth
C     ============================================

C  set up help arrays (independent of eigenvector)

	DO I = 1, NSTREAMS
	  TPA = ZERO
	  TMA = ZERO
	  DO L = M, NMOMENTS
	    TPA = TPA + PLMI_X0_P(I,L,N) * OMEGA_MOMS(N,L)
	    TMA = TMA + PLMI_X0_M(I,L,N) * OMEGA_MOMS(N,L)
	  ENDDO
	  DPI(I) = TPA
	  DMI(I) = TMA
	ENDDO

C  Set up eigenfunction norms, and save them

	DO AA = 1, NSTREAMS
	  NORM = ZERO
	  DO J = 1, NSTREAMS
	    J1 = J + NSTREAMS
	    T1 = XPOS(J,AA,N)  * XPOS(J,AA,N)
	    T2 = XPOS(J1,AA,N) * XPOS(J1,AA,N)
	    NORM = NORM + AX(J)*(T1-T2)
	  ENDDO
	  NORM_SAVED(AA) = NORM
	ENDDO

C  For each eigenstream, get the terms ATERM_SAVE and BTERM_SAVE

	DO AA = 1, NSTREAMS
	  SUM_LA = ZERO
	  SUM_LB = ZERO
	  DO I = 1, NSTREAMS
	    I1 = I + NSTREAMS
	    TPA = A(I)*(DPI(I)*XPOS(I,AA,N)+DMI(I)*XPOS(I1,AA,N))
	    TMA = A(I)*(DMI(I)*XPOS(I,AA,N)+DPI(I)*XPOS(I1,AA,N))
	    SUM_LA  = SUM_LA + TPA
	    SUM_LB  = SUM_LB + TMA
	  ENDDO
	  ATERM_SAVE(AA,N) = SUM_LA / NORM_SAVED(AA)
	  BTERM_SAVE(AA,N) = SUM_LB / NORM_SAVED(AA)
	  AGM(AA,N)  = ATERM_SAVE(AA,N) * GAMMA_M(AA,N)
	  BGP(AA,N)  = BTERM_SAVE(AA,N) * GAMMA_P(AA,N)
	ENDDO

C  5. Green function multipliers
C     ==========================

C  For each eigenstream

	DO AA = 1, NSTREAMS
	  GFUNC_DN(AA,N) = CFUNC(AA,N) * ATERM_SAVE(AA,N) * CONST
	  GFUNC_UP(AA,N) = DFUNC(AA,N) * BTERM_SAVE(AA,N) * CONST
	ENDDO

C  6. Set particular integral from Green function expansion
C     =====================================================

C  particular integrals at lower and upper boundaries

	DO I = 1, NSTREAMS
	  I1 = I + NSTREAMS
	  S_P_U = ZERO
	  S_P_L = ZERO
	  S_M_U = ZERO
	  S_M_L = ZERO
	  DO AA = 1, NSTREAMS
	    S_P_U = S_P_U + GFUNC_UP(AA,N)*XPOS(I1,AA,N)
	    S_M_U = S_M_U + GFUNC_UP(AA,N)*XPOS(I,AA,N)
	    S_P_L = S_P_L + GFUNC_DN(AA,N)*XPOS(I,AA,N)
	    S_M_L = S_M_L + GFUNC_DN(AA,N)*XPOS(I1,AA,N)
	  ENDDO
	  WUPPER(I,N) = S_P_U
	  WUPPER(I1,N) = S_M_U
	  WLOWER(I1,N) = S_M_L
	  WLOWER(I,N) = S_P_L
	ENDDO

C  Finish

	END

C

	SUBROUTINE LIDORT_USERSOLUTION
     I    ( LAYER, FOURIER )

C  Standard Solution components and Multipliers for user-defined output.
C  Only compute the ones we need (not all layers required)

C  Include files
C  =============

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  Subroutine input arguments

	INTEGER		 FOURIER, LAYER

C  Solutions at stream angles
C  ==========================

C  general flag

	IF  ( STERM_LAYERMASK_UP(LAYER) .OR.
     &        STERM_LAYERMASK_DN(LAYER) ) THEN

C  Homogeneous

	  CALL USER_HOMOG_SOLUTION ( LAYER, FOURIER )

C  Particular integrals

	  IF ( DO_CLASSICAL_SOLUTION ) THEN
	    CALL USER_CLASSICAL_SOLUTION ( LAYER, FOURIER )
	  ELSE
	    CALL USER_GREENFUNC_SOLUTION ( LAYER, FOURIER )
	  ENDIF

	ENDIF

C  Finish

	END

C

	SUBROUTINE USER_HOMOG_SOLUTION
     I    ( GIVEN_LAYER, FOURIER )

C  Include files
C  =============

C  input
C  -----

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of Legendre stuff (input to this module)

	INCLUDE '../include_s/LIDORT_LEGENDRE.VARS'

C  include files of setup stuff (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include files of solution variables (output to this module)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include file of Multiplier coefficients

	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

	INTEGER		 GIVEN_LAYER
	INTEGER		 FOURIER

C  Local variables
C  ---------------

	INTEGER		 UM, J, J1, L, N, M, AA
        DOUBLE PRECISION SUM_NEG, SUM_POS, POS1, POS2, NEG1, NEG2
        DOUBLE PRECISION RHO_P, RHO_M, SECMUI, ULP

C  Layer and Fourier

	N = GIVEN_LAYER
	M = FOURIER

C  Zeta constants (always required)

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  SECMUI = USER_SECANTS(UM)
	  DO AA = 1, NSTREAMS
	    RHO_P = SECMUI + KEIGEN(AA,N)
	    RHO_M = SECMUI - KEIGEN(AA,N)
	    ZETA_P(AA,UM,N) = ONE / RHO_P
	    ZETA_M(AA,UM,N) = ONE / RHO_M
	  ENDDO
	ENDDO

C  Eigenvector interpolation to user-defined angles
C  ------------------------------------------------

C  For each eigenvector

        DO AA = 1, NSTREAMS

C  For each moment, do inner sum over computational angles
C  for the positive and negative eigenvectors

	  DO L = M, NMOMENTS
            SUM_POS = ZERO
            SUM_NEG = ZERO
            DO  J = 1, NSTREAMS
	      J1 = J + NSTREAMS
	      POS1 = XPOS(J1,AA,N) * WT_LEGP(J,L)
 	      POS2 = XPOS(J,AA,N)  * WT_LEGM(J,L)
	      NEG1 = XNEG(J1,AA,N) * WT_LEGP(J,L)
	      NEG2 = XNEG(J,AA,N)  * WT_LEGM(J,L)
              SUM_POS = SUM_POS + POS1 + POS2
              SUM_NEG = SUM_NEG + NEG1 + NEG2
	    ENDDO
	    U_HELP_P(AA,L) = SUM_POS
 	    U_HELP_M(AA,L) = SUM_NEG
          ENDDO

C  Now sum over all harmonic contributions at each user-defined stream

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    SUM_POS = ZERO
	    SUM_NEG = ZERO
            DO L = M, NMOMENTS
	      ULP = U_LEG_P(UM,L) * OMEGA_MOMS(N,L)
              SUM_POS = SUM_POS + U_HELP_P(AA,L) * ULP
              SUM_NEG = SUM_NEG + U_HELP_M(AA,L) * ULP
            ENDDO
	    U_XPOS(UM,AA,N) = SUM_POS
	    U_XNEG(UM,AA,N) = SUM_NEG
	  ENDDO

C  end eigenvector loop

	ENDDO

C  Finish

	END

C

	SUBROUTINE USER_CLASSICAL_SOLUTION
     I    ( GIVEN_LAYER, FOURIER )

C  Include files
C  =============

C  input
C  -----

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model and control variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of Legendre stuff (input to this module)

	INCLUDE '../include_s/LIDORT_LEGENDRE.VARS'

C  include files of setup stuff (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include files of solution variables (output to this module)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

	INTEGER		 GIVEN_LAYER
	INTEGER		 FOURIER

C  Local variables
C  ---------------

	INTEGER		 UM, J, J1, L, N, M
        DOUBLE PRECISION SUM, POS1, POS2
C 	DOUBLE PRECISION HELP(0:MAXMOMENT)
 	DOUBLE PRECISION HELP1(0:MAXMOMENT)
 	DOUBLE PRECISION HELP2(0:MAXMOMENT)

C  Layer and Fourier

	N = GIVEN_LAYER
	M = FOURIER

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the boundary layer values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    IF ( DO_UPWELLING ) THEN
	      U_WPOS1(UM,N) = ZERO
	      U_WPOS2(UM,N) = ZERO
	    ENDIF
	    IF ( DO_DNWELLING ) THEN
	      U_WNEG1(UM,N) = ZERO
	      U_WNEG2(UM,N) = ZERO
	    ENDIF
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  For each moment & particulate, do inner sum over computational angles
C  (required for both directions)

 	DO L = M, NMOMENTS
          SUM = ZERO
          DO  J = 1, NSTREAMS
	    J1 = J + NSTREAMS
	    POS1 = WVEC(J1,N) * WT_LEGP(J,L)
 	    POS2 = WVEC(J,N)  * WT_LEGM(J,L)
            SUM = SUM + POS1 + POS2
	  ENDDO
	  W_HELP(L) = SUM
C 	  HELP(L)   = ( W_HELP(L) + LEG0_M(L,N) ) * OMEGA_MOMS(N,L)
 	  HELP1(L)  = LEG0_M(L,N) * OMEGA_MOMS(N,L)
 	  HELP2(L)  = W_HELP(L)   * OMEGA_MOMS(N,L)
       ENDDO

C  Now sum over all harmonic contributions at each user-defined stream
C  Distinguish between upwelling and downwelling

	IF ( DO_UPWELLING ) THEN
	  IF ( STERM_LAYERMASK_UP(N) ) THEN
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
              POS1 = ZERO
              DO L = M, NMOMENTS
                POS1 = POS1 + HELP1(L)*U_LEG_P(UM,L)
              ENDDO
	      U_WPOS1(UM,N) = POS1
              POS1 = ZERO
              DO L = M, NMOMENTS
                POS1 = POS1 + HELP2(L)*U_LEG_P(UM,L)
              ENDDO
	      U_WPOS2(UM,N) = POS1
	    ENDDO
	  ENDIF
	ENDIF

	IF ( DO_DNWELLING ) THEN
	  IF ( STERM_LAYERMASK_DN(N) ) THEN
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
              POS2 = ZERO
              DO L = M, NMOMENTS
                POS2 = POS2 + HELP1(L)*U_LEG_M(UM,L)
              ENDDO
	      U_WNEG1(UM,N) = POS2
              POS2 = ZERO
              DO L = M, NMOMENTS
                POS2 = POS2 + HELP2(L)*U_LEG_M(UM,L)
              ENDDO
	      U_WNEG2(UM,N) = POS2
	    ENDDO
	  ENDIF
	ENDIF

C  Finish

	END

C

	SUBROUTINE USER_GREENFUNC_SOLUTION
     I    ( GIVEN_LAYER, FOURIER )

C  Include files
C  =============

C  input
C  -----

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of Legendre stuff (input to this module)

	INCLUDE '../include_s/LIDORT_LEGENDRE.VARS'

C  include files of setup stuff (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include files of solution variables (output to this module)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index and Fourier number (inputs)

	INTEGER		 GIVEN_LAYER
	INTEGER		 FOURIER

C  Local variables
C  ---------------

	INTEGER		 UM, L, N, M
        DOUBLE PRECISION POS1, POS2

C  Layer and Fourier

	N = GIVEN_LAYER
	M = FOURIER

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the boundary layer values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    U_WPOS1(UM,N) = ZERO
	    U_WNEG1(UM,N) = ZERO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  For each moment do inner sum over computational angles

 	DO L = M, NMOMENTS
	  W_HELP(L) = LEG0_M(L,N)*OMEGA_MOMS(N,L)
        ENDDO

C  Now sum over all harmonic contributions at each user-defined stream

	DO UM = LOCAL_UM_START, N_USER_STREAMS
          POS1 = ZERO
          POS2 = ZERO
          DO L = M, NMOMENTS
            POS1 = POS1 + W_HELP(L)*U_LEG_P(UM,L)
            POS2 = POS2 + W_HELP(L)*U_LEG_M(UM,L)
          ENDDO
	  U_WPOS1(UM,N) = POS1
	  U_WNEG1(UM,N) = POS2
	ENDDO

C  Finish

	END

C

	SUBROUTINE LIDORT_LEGENDRE_SETUP ( FOURIER )

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  include file of control and model and setup variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include file of Legendre polynomials and associated multipliers
C  ( output from this module)

	INCLUDE '../include_s/LIDORT_LEGENDRE.VARS'

C  input

	INTEGER		 FOURIER

C  local variables

	DOUBLE PRECISION CFPLG(0:MAXMOMENT), WT
	INTEGER		 M, I, J, L, LPM, N

C  Set integer M = Fourier number

	M = FOURIER

C  Legendre polynomials
C  --------------------

C  .. positive computational angle streams

	DO I = 1, NSTREAMS
	  CALL CFPLGARR (MAXMOMENT, NMOMENTS, M, X(I), CFPLG)
	  DO L = M, NMOMENTS
	    LEG_P(I,L) = CFPLG(L)
	  ENDDO
	ENDDO

C  .. negative streams by symmetry relations

	DO L = M, NMOMENTS
	  LPM = L + M
	  IF (MOD(LPM,2).EQ.0) THEN
	    DO I = 1, NSTREAMS
	      LEG_M(I,L) = LEG_P(I,L)
	    ENDDO
	  ELSE
	    DO I = 1, NSTREAMS
	      LEG_M(I,L) = - LEG_P(I,L)
	    ENDDO
	  ENDIF
	ENDDO

C  ..solar zenith angle values

	IF ( DO_QSREFRAC_BEAM ) THEN
	  DO N = 1, NLAYERS
	    CALL CFPLGARR (MAXMOMENT, NMOMENTS, M, LOCAL_SZA(N), LEG0_P)
	    DO L = M, NMOMENTS
	      LPM = L + M
	      IF (MOD(LPM,2).EQ.0) THEN
	        LEG0_M(L,N) = LEG0_P(L)
	      ELSE
	        LEG0_M(L,N) = - LEG0_P(L)
	      ENDIF    
	      DO I = 1, NSTREAMS
	        PLMI_X0_P(I,L,N) = LEG0_P(L)*LEG_P(I,L)
	        PLMI_X0_M(I,L,N) = LEG0_P(L)*LEG_M(I,L)
	      ENDDO
	    ENDDO
	  ENDDO
	ELSE
	  CALL CFPLGARR (MAXMOMENT, NMOMENTS, M, X0, LEG0_P)
	  DO L = M, NMOMENTS
	    LPM = L + M
	    IF (MOD(LPM,2).EQ.0) THEN
	      LEG0_M(L,1) = LEG0_P(L)
	    ELSE
	      LEG0_M(L,1) = - LEG0_P(L)
	    ENDIF
	    DO I = 1, NSTREAMS
	      PLMI_X0_P(I,L,1) = LEG0_P(L)*LEG_P(I,L)
	      PLMI_X0_M(I,L,1) = LEG0_P(L)*LEG_M(I,L)
	    ENDDO
	    DO N = 2, NLAYERS
	      LEG0_M(L,N) = LEG0_M(L,1)
	      DO I = 1, NSTREAMS
	        PLMI_X0_P(I,L,N) = PLMI_X0_P(I,L,1)
	        PLMI_X0_M(I,L,N) = PLMI_X0_M(I,L,1)
	      ENDDO
	    ENDDO
	  ENDDO
	ENDIF

C  set up phi_lm * pleg0 (Siewert notation)

C	DO L = M, NMOMENTS
C	  DO I = 1, NSTREAMS
C	    PLMI_X0_P(I,L) = LEG0_P(L)*LEG_P(I,L)
C	    PLMI_X0_M(I,L) = LEG0_P(L)*LEG_M(I,L)
C	  ENDDO
C	ENDDO

C  set up products of Associated Legendre polynomials
C  --------------------------------------------------

C  set up PLM(x).PMM(x)

	DO I = 1, NSTREAMS
	  DO J = 1, NSTREAMS
	    DO L = M, NMOMENTS
	      PLMI_PLMJ_P(I,J,L) = LEG_P(I,L) * LEG_P(J,L)
	      PLMI_PLMJ_M(I,J,L) = LEG_P(I,L) * LEG_M(J,L)
   	    ENDDO
	  ENDDO
	ENDDO

C  associated products

        DO  J = 1, NSTREAMS
	  WT = HALFA(J)
	  DO L = M, NMOMENTS
	    WT_LEGP(J,L) = LEG_P(J,L) * WT
	    WT_LEGM(J,L) = LEG_M(J,L) * WT
	  ENDDO
	ENDDO

C  Finish

	END

C

	SUBROUTINE LIDORT_USERLEGENDRE_SETUP ( FOURIER )

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  include file of control and model variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of Legendre polynomials and associated multipliers
C  ( output from this module)

	INCLUDE '../include_s/LIDORT_LEGENDRE.VARS'

C  input

	INTEGER		 FOURIER

C  local variables

	DOUBLE PRECISION CFPLG(0:MAXMOMENT)
	INTEGER		 M, L, UM, LPM

C  Set integer M = Fourier number

	M = FOURIER

C  Legendre polynomials defined on user stream angles
C  --------------------------------------------------

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  CALL CFPLGARR
     &       (  MAXMOMENT, NMOMENTS, M, USER_STREAMS(UM), CFPLG)
	  DO L = M, NMOMENTS
	    U_LEG_P(UM,L) = CFPLG(L)
	    LPM = L + M
	    IF (MOD(LPM,2).EQ.0) THEN
	      U_LEG_M(UM,L) = U_LEG_P(UM,L)
	    ELSE
	      U_LEG_M(UM,L) = - U_LEG_P(UM,L)
	    ENDIF  
	  ENDDO    
	ENDDO

C  Finish

	END
