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

	SUBROUTINE LIDORT_BVPROBLEM
     I     ( DO_INCLUDE_ALBEDO,
     I       DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS,
     I       FOURIER, R2, DIRECT_BEAM,
     O       STATUS )

C  Solves the boundary value problem.

C  include files
C  -------------

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of model input variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of main model variable (I/O to this module)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  subroutine arguments
C  --------------------

C  Albedo factor ( = 2xALBEDO for fourier = 0)

	DOUBLE PRECISION R2

C  direct beam contribution

	DOUBLE PRECISION DIRECT_BEAM(MAXSTRM)

C  Fourier index

	INTEGER		 FOURIER

C  local control flag

	LOGICAL		 DO_INCLUDE_ALBEDO
	LOGICAL		 DO_INCLUDE_SURFEMISS
	LOGICAL		 DO_INCLUDE_DIRECTBEAM
C	LOGICAL		 DO_INCLUDE_THERMAL		! REMOVED

C  status flag

	INTEGER		 STATUS

C  local variables
C  ---------------

	INTEGER		 I, C0, I1, LAY, AA
	CHARACTER*(70)	 MAIL, TRACE
	INTEGER		 INFO
	CHARACTER*3	 CI

C  Intialize status

	STATUS = LIDORT_SUCCESS

C  Set up boundary value problem
C  -----------------------------

C  Additional setups for the albedo layer

	CALL LIDORT_ALBEDO_SETUPS
     I  ( DO_INCLUDE_ALBEDO, FOURIER, R2 )

C  ########## Old code
C  set up boundary values matrix MAT2 (the "A" as in AX=B)
C	CALL LIDORT_BCMATRIX_SETUP
C     I     ( DO_INCLUDE_ALBEDO )
C  Compress solution matrix MAT2 to band storage form BANDMAT2
C	CALL LIDORT_BAND_STORAGE
C     &     ( NTOTAL, N_SUBDIAG, N_SUPDIAG, MAT2, BANDMAT2 )
C  #### end old code

C  initialize compression matrix

C	IF (FOURIER .EQ. 0 ) THEN
	  CALL LIDORT_BCMATRIX_INIT
C	ENDIF

C  set up boundary values matrix in compressed form (the "A" as in AX=B)

	CALL LIDORT_BCMATRIX_SETUP_NEW
     I     ( DO_INCLUDE_ALBEDO )

C  set up Column COL2 for solution vector (the "B" as in AX=B)

	CALL LIDORT_BCCOLUMN_SETUP
     I     ( DO_INCLUDE_ALBEDO,
     I       DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS,
     I       DIRECT_BEAM )

C  Solve the boundary problem for this Fourier component
C  -----------------------------------------------------

C  LAPACK LU-decomposition for band matrix

        CALL DGBTRF
     &     ( NTOTAL, NTOTAL, N_SUBDIAG, N_SUPDIAG,
     &       BANDMAT2, MAXBANDTOTAL, IPIVOT, INFO )

C  (Error tracing)

	IF ( INFO .GT. 0 ) THEN
	  WRITE(CI, '(I3)' ) INFO
	  MAIL  = 'Singular matrix, u(i,i)=0, for i = '//CI
	  TRACE = 'DGBTRF call in LIDORT_BVPROBLEM'
	  STATUS = LIDORT_SERIOUS
	  CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
	  RETURN
	ELSE IF ( INFO .LT. 0 ) THEN
	  WRITE(CI, '(I3)' ) INFO
	  MAIL  = 'argument i illegal value, for i = '//CI
	  TRACE = 'first DGBTRS call in LIDORT_FOURIER_MASTER'
	  STATUS = LIDORT_SERIOUS
	  CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
	  RETURN
	ENDIF

C  LAPACK substitution using RHS column vector COL2

        CALL DGBTRS
     &     ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1,
     &        BANDMAT2, MAXBANDTOTAL, IPIVOT,
     &        COL2, MAXTOTAL, INFO )

C  (error tracing)

	IF ( INFO .LT. 0 ) THEN
	  WRITE(CI, '(I3)' ) INFO
	  MAIL  = 'argument i illegal value, for i = '//CI
	  TRACE = 'DGBTRS call in LIDORT_BVPROBLEM'
	  STATUS = LIDORT_SERIOUS
	  CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
	  RETURN
	 ENDIF

C  Set integration constants LCON and MCON for -/+ eigensolutions, all layers
C  --------------------------------------------------------------------------

	DO LAY = 1, NLAYERS
	  C0 = (LAY-1)*NSTR2
	  DO I = 1, NSTREAMS
	    I1 = I+NSTREAMS
	    LCON(I,LAY) = COL2(C0+I,1)
	    MCON(I,LAY) = COL2(C0+I1,1)
	  ENDDO
	ENDDO

C  Associated quantities

	DO LAY = 1, NLAYERS
	  DO I = 1, NSTR2
	    DO AA = 1, NSTREAMS
	      LCON_XVEC(I,AA,LAY) = LCON(AA,LAY) * XPOS(I,AA,LAY)
	      MCON_XVEC(I,AA,LAY) = MCON(AA,LAY) * XNEG(I,AA,LAY)
	    ENDDO
	  ENDDO
	ENDDO

C  Finish

	END

C

	SUBROUTINE LIDORT_ALBEDO_SETUPS
     I  ( DO_INCLUDE_ALBEDO, FOURIER, R2 )

C  Additional sums for the final albedo-reflecting layer

C  include file of dimensiopns and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'

C  include file of solution stuff (i/o to this module)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  control input

	LOGICAL		 DO_INCLUDE_ALBEDO
C	LOGICAL		 DO_INCLUDE_THERMAL		! REMOVED
	INTEGER		 FOURIER
	DOUBLE PRECISION R2

C  local variables

	DOUBLE PRECISION AXBID(MAXSTRM,MAXSTRM), SUM, H_1, H_2
	INTEGER		 I, AA, J, LAY

C  No values if albedo flag not set
C  ================================

	IF ( .NOT. DO_INCLUDE_ALBEDO ) THEN
	  DO I = 1, NSTREAMS
	    R2_BEAM(I)    = ZERO
	    DO J = 1, NSTREAMS
 	      R2_HOMP(I,J) = ZERO
	      R2_HOMM(I,J) = ZERO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF

C  Integrate Downward streams of particular solutions
C  ==================================================

	LAY = NLAYERS

C  For Lambertian reflectance, all streams are the same

	IF ( DO_LAMBERTIAN_ALBEDO ) THEN

	  SUM = ZERO
	  DO J = 1, NSTREAMS
	    SUM = SUM + AX(J)*WLOWER(J,LAY)
	  ENDDO
	  DO I = 1, NSTREAMS
	    R2_BEAM(I) = R2*SUM
	  ENDDO

C  Non-Lambertian, bidirectional surface
C  ..  save AX*BIDIREC for later use

	ELSE

	  DO I = 1, NSTREAMS
	    SUM = ZERO
	    DO J = 1, NSTREAMS
	      AXBID(I,J) = AX(J)*BIREFLEC(FOURIER,I,J)
	      SUM = SUM + AXBID(I,J) * WLOWER(J,LAY)
	    ENDDO
	    R2_BEAM(I) = R2*SUM
	  ENDDO

	ENDIF    

C  Integrate Downward streams of Homogeneous solutions
C  ===================================================

C  For Lambertian reflectance, all streams are the same

	IF ( DO_LAMBERTIAN_ALBEDO ) THEN

	  DO AA = 1, NSTREAMS
	    DO I = 1, NSTREAMS
	      H_1 = ZERO
	      H_2 = ZERO
	      DO J = 1, NSTREAMS
	        H_1 = H_1 + AX(J)*XPOS(J,AA,LAY)
	        H_2 = H_2 + AX(J)*XNEG(J,AA,LAY)
	      ENDDO
	      R2_HOMP(I,AA) = H_1*R2
	      R2_HOMM(I,AA) = H_2*R2
	    ENDDO
	  ENDDO

C  Non-Lambertian, bidirectional surface

	ELSE

	  DO AA = 1, NSTREAMS
	    DO I = 1, NSTREAMS
	      H_1 = ZERO
	      H_2 = ZERO
	      DO J = 1, NSTREAMS
	        H_1 = H_1 + AXBID(I,J)*XPOS(J,AA,LAY)
	        H_2 = H_2 + AXBID(I,J)*XNEG(J,AA,LAY)
	      ENDDO
	      R2_HOMP(I,AA) = H_1*R2
	      R2_HOMM(I,AA) = H_2*R2
	    ENDDO
	  ENDDO

	ENDIF

C  Finish

	END

C

	SUBROUTINE LIDORT_BCMATRIX_INIT

C  Initialisee the compressed matrix

C  include files
C  -------------

C  include file of dimensiopns and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of main model variables (local to this module)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  Local variables
C  ---------------

	INTEGER	 NMIN(MAXTOTAL), NMAX(MAXTOTAL)
	INTEGER  I, J, N3, JS, JF, IS, LF, L, KALL, I1

c  compression row indices

	DO J = 1, N_SUPDIAG + 1
	  NMIN(J) = 1
	ENDDO
	DO J = N_SUPDIAG + 2, NTOTAL
	  NMIN(J) = J - N_SUPDIAG
	ENDDO
	DO J = 1, NTOTAL - N_SUBDIAG
	  NMAX(J) = J + N_SUBDIAG
	ENDDO
	DO J = NTOTAL - N_SUBDIAG + 1, NTOTAL
	  NMAX(J) = NTOTAL
	ENDDO

	KALL = N_SUBDIAG + N_SUPDIAG + 1
	DO I = 1, NTOTAL
	  DO J = 1, NTOTAL
	    IF ( (I.GE.NMIN(J)) .AND. (I.LE.NMAX(J)) ) THEN
	      BMAT_ROWMASK(I,J) = KALL + I - J
	    ENDIF
	  ENDDO
	ENDDO

C  compression matrix zeroing

	N3 = NSTR2 + NSTREAMS
	LF = NLAYERS - 2

C  upper band top

	JS = NSTR2 + 1
	JF = N3 - 1
	DO I = 1, NSTREAMS
	  DO J = JS, JF + I
	    BANDMAT2(BMAT_ROWMASK(I,J),J) = ZERO
	  ENDDO
	ENDDO

C  upper band

	DO L = 1, LF
	  IS = L*NSTR2 - NSTREAMS + 1
	  JS = IS + N3
	  JF = JS - 1
	  DO I = 1, NSTR2-1
	    I1 = I + IS
	    DO J = JS, JF + I
	      BANDMAT2(BMAT_ROWMASK(I1,J),J) = ZERO
	    ENDDO
	  ENDDO
	ENDDO

C  lower band

	DO L = 1, LF
	  IS = L*NSTR2 + NSTREAMS
	  JS = IS - N3 + 1
	  JF = IS - NSTREAMS
	  DO I = 1, NSTR2-1
	    I1 = I + IS
	    DO J = JS + I, JF
	      BANDMAT2(BMAT_ROWMASK(I1,J),J) = ZERO
	    ENDDO
	  ENDDO
	ENDDO

C  lower band bottom

	JS = LF * NSTR2 + 1
	IS = JS + N3 - 1
	JF = IS - NSTREAMS
	DO I = 1, NSTREAMS
	  I1 = I + IS
	  DO J = JS + I, JF
	    BANDMAT2(BMAT_ROWMASK(I1,J),J) = ZERO
	  ENDDO
	ENDDO

C  finish

	END

C

	SUBROUTINE LIDORT_BCMATRIX_SETUP_NEW
     I     (  DO_INCLUDE_ALBEDO )

C  Fills up the compressed matrix directly

C  include files
C  -------------

C  include file of dimensiopns and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of main model variables (local to this module)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  input arguments
C  ---------------

	LOGICAL		 DO_INCLUDE_ALBEDO

C  Local variables
C  ---------------

	INTEGER		 I,EP,EM,N,N1,I1,AA
	INTEGER		 C0,CE_OFFSET,CEM,CEP,CEM1,CEP1,CM,CP
	DOUBLE PRECISION XPNET, XMNET

C  top BC for layer 1: no downward diffuse radiation

	N = 1
	DO I = 1, NSTREAMS
	  DO EP = 1, NSTREAMS
	    EM = EP + NSTREAMS
	    BANDMAT2(BMAT_ROWMASK(I,EP),EP)  =
     &                   XPOS(I,EP,N)
	    BANDMAT2(BMAT_ROWMASK(I,EM),EM)  =
     &                   XNEG(I,EP,N)*T_DELT_EIGEN(EP,N)
	  ENDDO
	ENDDO

C  intermediate layer boundaries (will not be done if NLAYERS = 1 )

	C0 = - NSTREAMS
	DO N = 2, NLAYERS
	  N1 = N - 1
          C0   = C0 + NSTR2
	  CE_OFFSET = C0 - NSTREAMS
	  DO I = 1, NSTR2
	    CM = C0 + I
	    DO EP = 1, NSTREAMS
	      CEP = CE_OFFSET + EP
	      CEM = CEP + NSTREAMS
	      CEP1 = CEP + NSTR2
	      CEM1 = CEM + NSTR2
	      BANDMAT2(BMAT_ROWMASK(CM,CEP),CEP)   =
     &                   T_DELT_EIGEN(EP,N1)*XPOS(I,EP,N1)
	      BANDMAT2(BMAT_ROWMASK(CM,CEM),CEM)   =
     &                   XNEG(I,EP,N1)
	      BANDMAT2(BMAT_ROWMASK(CM,CEP1),CEP1) =
     &                   -XPOS(I,EP,N)
	      BANDMAT2(BMAT_ROWMASK(CM,CEM1),CEM1) =
     &                   -T_DELT_EIGEN(EP,N)*XNEG(I,EP,N)
	    ENDDO
	  ENDDO
	ENDDO

C  bottom BC (with albedo additions if flagged)

	N = NLAYERS
	C0 = C0 + NSTR2
	CE_OFFSET = C0 - NSTREAMS

	IF ( DO_INCLUDE_ALBEDO ) THEN
	  DO I = 1, NSTREAMS
	    CP = C0 + I
	    I1 = I + NSTREAMS
	    DO AA = 1, NSTREAMS
	      CEP = CE_OFFSET + AA
	      CEM = CEP + NSTREAMS
	      XPNET = XPOS(I1,AA,N) - R2_HOMP(I,AA)
	      XMNET = XNEG(I1,AA,N) - R2_HOMM(I,AA)
	      BANDMAT2(BMAT_ROWMASK(CP,CEP),CEP) =
     &                   T_DELT_EIGEN(AA,N) * XPNET
	      BANDMAT2(BMAT_ROWMASK(CP,CEM),CEM) =
     &                   XMNET
	    ENDDO
	  ENDDO
	ELSE
	  DO I = 1, NSTREAMS
	    CP = C0 + I
	    I1 = I + NSTREAMS
	    DO AA = 1, NSTREAMS
	      CEP = CE_OFFSET + AA
	      CEM = CEP + NSTREAMS
	      BANDMAT2(BMAT_ROWMASK(CP,CEP),CEP) =
     &                  T_DELT_EIGEN(AA,N) * XPOS(I1,AA,N)
	      BANDMAT2(BMAT_ROWMASK(CP,CEM),CEM) =
     &                  XNEG(I1,AA,N)
	    ENDDO
	  ENDDO
	ENDIF

C  finish

	END

C

	SUBROUTINE LIDORT_BCCOLUMN_SETUP
     I     ( DO_INCLUDE_ALBEDO,
     I       DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS,
     I       DIRECT_BEAM )

C  include files
C  -------------

C  include file of dimensiopns and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'

C  include files of solution variables (local to this module)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  input arguments
C  ---------------

	LOGICAL		 DO_INCLUDE_DIRECTBEAM
	LOGICAL		 DO_INCLUDE_ALBEDO, DO_INCLUDE_SURFEMISS
	DOUBLE PRECISION DIRECT_BEAM(MAXSTRM)
C	LOGICAL		 DO_INCLUDE_THERMAL		! REMOVED

C  Local variables
C  ---------------

	INTEGER		 I,LAY,LAY1,I1,C0,CM

C  zero column vector

	DO I = 1, NTOTAL
	  COL2(I,1) = ZERO
	ENDDO

C  Upper boundary for layer 1: no downward diffuse radiation
C  ---------------------------------------------------------

	LAY = 1
	DO I = 1, NSTREAMS
	  COL2(I,1)   = - WUPPER(I,LAY)
	ENDDO

C  intermediate layer boundaries (will not be done if NLAYERS = 1 )
C  -----------------------------

	DO LAY = 2, NLAYERS

	  LAY1 = LAY - 1
          C0 = LAY1*NSTR2 - NSTREAMS
	  DO I = 1, NSTR2
	    CM = C0 + I
	    COL2(CM,1) = WUPPER(I,LAY) - WLOWER(I,LAY1)
	  ENDDO

	ENDDO

C  lowest (surface) boundary with albedo (diffuse radiation terms only)
C  -------------------------------------

	LAY = NLAYERS
	C0 = (LAY-1)*NSTR2 + NSTREAMS

C  with non-zero albedo, include integrated downward reflectances

	IF ( DO_INCLUDE_ALBEDO ) THEN

	  DO I = 1, NSTREAMS
	    I1 = I + NSTREAMS
	    CM = C0 + I
	    COL2(CM,1) = - WLOWER(I1,LAY) + R2_BEAM(I)
	  ENDDO

C  no albedo, similar code excluding integrated reflectance

	ELSE

	  DO I = 1, NSTREAMS
	    I1 = I + NSTREAMS
	    CM = C0 + I
	    COL2(CM,1) = - WLOWER(I1,LAY)
	  ENDDO

	ENDIF

C  Add direct beam solution (only to final level)
C  ----------------------------------------------

	IF ( DO_INCLUDE_DIRECTBEAM ) THEN
	  IF ( DO_INCLUDE_ALBEDO ) THEN
	    DO I = 1, NSTREAMS
	      CM = C0 + I
	      COL2(CM,1) = COL2(CM,1) + DIRECT_BEAM(I)
	    ENDDO
	  ENDIF
	ENDIF

C  Add thermal emission of ground surface (only to final level)
C  ------------------------------------------------------------

	IF ( DO_INCLUDE_SURFEMISS ) THEN
	  DO I = 1, NSTREAMS
	    CM = C0 + I
	    COL2(CM,1) = COL2(CM,1) + SURFBB * EMISSIVITY(I)
	  ENDDO
	ENDIF

C  finish

	END

C#################################################################
C  the following is no longer required
C  (compressed band matrix is filled up directly)
C	SUBROUTINE LIDORT_BAND_STORAGE
C     I     ( N, KL, KU, MAT, BANDMAT )
C  Compression module for preparation of Band matrix required for
C  use of LAPACK modules.
C  include file of dimensions and numbers
C	INCLUDE '../include_s/LIDORT.PARS'
C  input
C	INTEGER		 N, KL, KU
C	DOUBLE PRECISION MAT(MAXTOTAL,MAXTOTAL)
C  output
C	DOUBLE PRECISION BANDMAT(MAXBANDTOTAL,MAXTOTAL)
C  local variables
C	INTEGER		 I,J, NMIN(MAXTOTAL), NMAX(MAXTOTAL), IROW
C  set up max/min
C	DO J=1,N
C	  NMIN(J) = MAX(J-KU,1)
C	  NMAX(J) = MIN(J+KL,N)
C	ENDDO
C  compress entries
C	DO I=1,N
C	  DO J=1,N
C	    IF ((I.GE.NMIN(J)).AND.(I.LE.NMAX(J))) THEN
C	      IROW = KL+KU+I-J+1
C	      BANDMAT(IROW,J)=MAT(I,J)
C	    END IF
C	  END DO
C	END DO
C  Finish
C	RETURN
C	END
C#################################################################
