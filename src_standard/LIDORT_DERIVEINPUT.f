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
C #  Version :	  2.2					        #
C #  Release Date   January 2001				#
C #		  					        #
C ###############################################################

	SUBROUTINE LIDORT_DERIVE_INPUT

C  Read, check and write to file of all control input

C  include files
C  -------------

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  include files with input variables to be checked out

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'

C  local variables
C  ---------------

	INTEGER		 INDEX_ANGLES ( MAX_OUT_STREAMS )
	DOUBLE PRECISION ALL_ANGLES ( MAX_OUT_STREAMS ), TAU
	INTEGER		 I, I1, UT, N, N1, UTA, NSTART, JINDEX, J
	LOGICAL		 LOOP
	INTEGER		 CONV_START_STREAMS

C  set additional numbers (derived input)
C  ======================

C  Mode of operation

	IF ( DO_FULLRAD_MODE ) THEN
	  IF ( DO_SSCORRECTION ) THEN
	    DO_MSMODE_LIDORT = .TRUE.
	  ELSE
	    DO_MSMODE_LIDORT = .FALSE.
	  ENDIF
	ELSE
	  DO_MSMODE_LIDORT = .TRUE.
	ENDIF
	    
C  Directional indices

	IF ( DO_UPWELLING .AND. DO_DNWELLING ) THEN
	  N_DIRECTIONS = 2
	  WHICH_DIRECTIONS(1) = UPIDX
	  WHICH_DIRECTIONS(2) = DNIDX
	ELSE
	  N_DIRECTIONS = 1
	  WHICH_DIRECTIONS(2) = 0
	  IF ( DO_UPWELLING ) THEN
	    WHICH_DIRECTIONS(1) = UPIDX
	  ELSE IF ( DO_DNWELLING) THEN
	    WHICH_DIRECTIONS(1) = DNIDX
	  ENDIF
	ENDIF

C  Number of moments

	IF ( DO_RAYLEIGH_ONLY ) THEN
	  NMOMENTS = 2
	ENDIF
	IF ( DO_ISOTROPIC_ONLY ) THEN
	  NMOMENTS = 0
	ENDIF
	IF ( .NOT.DO_RAYLEIGH_ONLY.AND..NOT.DO_ISOTROPIC_ONLY ) THEN
	  NMOMENTS = 2 * NSTREAMS - 1
	ENDIF

C  total quadratures (up and down)

	NSTR2 = 2*NSTREAMS

C  set auxiliary quantities

	DO I = 1, NSTREAMS
	  AX(I) = X(I)*A(I)
	  HALFA(I) = HALF * A(I)
	  XANG(I) = DACOS(X(I))/DEG_TO_RAD
	ENDDO

C  size of boundary value problem matrices and vectors

	NTOTAL = NLAYERS*NSTR2

C  number of sub and super diagonals in band matrix (boundary value problem)

	N_SUBDIAG = 3*NSTREAMS - 1
	N_SUPDIAG = 3*NSTREAMS - 1

C  solar zenith angle cosines

	X0 = DCOS ( SUN0 * DEG_TO_RAD )
	IF ( DO_QSREFRAC_BEAM ) THEN
	  DO N = 1, NLAYERS
	    SUN_SZA_COSINES(N) = DCOS ( SUNLOCAL_INPUT(N) * DEG_TO_RAD )
	  ENDDO
	ENDIF

C  Set the angle masks
C  -------------------

C  initialize

	DO I = 1, MAX_USER_STREAMS
	  USEROUTPUT_INDEX(I) = 0
	ENDDO
	DO I = 1, MAXSTRM
	  QUADOUTPUT_INDEX(I) = 0
	ENDDO

C  If quadrature output is required
C  --------------------------------

	IF ( DO_QUAD_OUTPUT ) THEN

C  .. For quadrature only, OUT_ANGLES are just the quadrature angles

	  IF ( .NOT. DO_USER_STREAMS ) THEN

C  .. number of streams for convergence and output

	    N_OUT_STREAMS  = NSTREAMS
	    CONV_START_STREAMS = 1

	    DO I = 1, N_OUT_STREAMS
	      I1 = NSTREAMS + 1 - I
	      QUADOUTPUT_INDEX(I) = I1
	      OUT_ANGLES(I) = XANG(I1)
	    ENDDO

C  For Quadrature AND user-defined stream angles

	  ELSE

C  .. number of streams for output

	    N_OUT_STREAMS  = NSTREAMS + N_USER_STREAMS

C   .. Set a dummy array ALL_ANGLES of all angles
C      ( order first with the user streams )

	    DO I = 1, N_USER_STREAMS
	      ALL_ANGLES(I) = USER_ANGLES_INPUT(I)
	    ENDDO
	    DO I = 1, NSTREAMS
	      I1 = I + N_USER_STREAMS
	      ALL_ANGLES(I1) = XANG(I)
	    ENDDO

C   .. Use an indexing routine INDEXX to sort ALL_ANGLES

	    CALL INDEXX(N_OUT_STREAMS, ALL_ANGLES, INDEX_ANGLES )

C   .. Set the output angles according to this indexing

	    DO I = 1, N_OUT_STREAMS
	      OUT_ANGLES(I) = ALL_ANGLES(INDEX_ANGLES(I))
	    ENDDO

C   .. Set the masks for quadrature and user-defined streams

	    DO I = 1, NSTREAMS
	      I1 = I + N_USER_STREAMS
	      LOOP = .TRUE.
	      J = 0
	      DO WHILE ( LOOP )
	        J = J + 1
	        JINDEX = INDEX_ANGLES(J)
	        IF ( JINDEX .EQ. I1 ) LOOP = .FALSE.
	      ENDDO
	      QUADOUTPUT_INDEX(I) = J
	    ENDDO

	    DO I = 1, N_USER_STREAMS
	      LOOP = .TRUE.
	      J = 0
	      DO WHILE ( LOOP )
	        J = J + 1
	        JINDEX = INDEX_ANGLES(J)
	        IF ( JINDEX .EQ. I ) LOOP = .FALSE.
	      ENDDO
	      USEROUTPUT_INDEX(I) = J
	    ENDDO

C  number of streams for convergence (uses zenith tolerance)

	    J = 0
	    DO I = 1, N_USER_STREAMS
	      IF ( DABS(ONE-DCOS(OUT_ANGLES(I)*DEG_TO_RAD)).LT.
     &                     ZENITH_TOLERANCE ) THEN
	        J = J + 1
	      ENDIF
	    ENDDO
	    CONV_START_STREAMS = J + 1

	  ENDIF

C  No Quadratures - only user-defined streams

	ELSE

	  IF ( DO_USER_STREAMS ) THEN

C  ..  number of streams for output

	    N_OUT_STREAMS  = N_USER_STREAMS

C   .. Index the user-defined angles and set the output anges & masks
C      ( Check already made that User-defined option is set)

	    DO I = 1, N_USER_STREAMS
	      ALL_ANGLES(I) = USER_ANGLES_INPUT(I)
	    ENDDO

	    IF ( N_OUT_STREAMS .EQ. 1 ) THEN
	      OUT_ANGLES(1) = ALL_ANGLES(1)
	      USEROUTPUT_INDEX(1) = 1
	    ELSE
	      CALL INDEXX
     &          ( N_OUT_STREAMS, USER_ANGLES_INPUT, INDEX_ANGLES )
	      DO I = 1, N_OUT_STREAMS
	        OUT_ANGLES(I) = ALL_ANGLES(INDEX_ANGLES(I))
	        LOOP = .TRUE.
	        J = 0
	        DO WHILE ( LOOP )
	          J = J + 1
	          JINDEX = INDEX_ANGLES(J)
	          IF ( JINDEX .EQ. I ) LOOP = .FALSE.
	        ENDDO
	        USEROUTPUT_INDEX(I) = J
	      ENDDO
	    ENDIF

C  number of streams for convergence (uses zenith tolerance)

	    J = 0
	    DO I = 1, N_USER_STREAMS
	      IF ( DABS(ONE-DCOS(OUT_ANGLES(I)*DEG_TO_RAD)).LT.
     &                     ZENITH_TOLERANCE ) THEN
	        J = J + 1
	      ENDIF
	    ENDDO
	    CONV_START_STREAMS = J + 1

	  ENDIF

	ENDIF

C  User stream cosines and secants

	IF ( DO_USER_STREAMS ) THEN
	  DO I = 1, N_USER_STREAMS
	    USER_STREAMS(I) = DCOS(DEG_TO_RAD*USER_ANGLES_INPUT(I))
	    USER_SECANTS(I) = ONE / USER_STREAMS(I)
	  ENDDO
	ENDIF

C  number of tests to be applied for convergence

	N_CONV_STREAMS = N_OUT_STREAMS - CONV_START_STREAMS + 1
	N_CONVTESTS = N_USER_RELAZMS * N_CONV_STREAMS * N_DIRECTIONS
	N_CONVTESTS = N_CONVTESTS * N_OUT_USERTAUS 

C  Sort out User optical depths
C  ----------------------------

C  ( These should be checked beforehand to lie within range of TAUGRID )

	IF ( DO_USER_TAUS ) THEN

C  Sort in ascending order

	  IF ( N_OUT_USERTAUS .GT. 1 ) THEN
	    CALL HPSORT(N_OUT_USERTAUS,USER_TAUS_INPUT)
	  ENDIF

C  mark all optical depths not equal to layer boundary values

	  NSTART = 0
	  UT = 0
	  DO UTA = 1, N_OUT_USERTAUS
	    TAU = USER_TAUS_INPUT(UTA)
	    OFFGRID_UTAU_OUTFLAG(UTA)  = .FALSE.
	    OFFGRID_UTAU_OUTINDEX(UTA) =   0
	    LOOP = .TRUE.
	    N = NSTART
	    DO WHILE (LOOP.AND.N.LT.NLAYERS)
	      N = N + 1
	      IF ( (TAU-TAUGRID_INPUT(N-1)).GT.SMALLNUM .AND.
     &             (TAUGRID_INPUT(N)-TAU).GT.SMALLNUM ) THEN
	        UT = UT + 1
	        OFFGRID_UTAU_OUTINDEX(UTA) =   UT
	        OFFGRID_UTAU_OUTFLAG(UTA)  = .TRUE.
	        OFFGRID_UTAU_LAYERIDX(UT) = N
	        UTAU_LEVEL_MASK_UP(UTA) = N
	        UTAU_LEVEL_MASK_DN(UTA) = N - 1
	        LOOP = .FALSE.
	        NSTART = N - 1
	      ENDIF
	    ENDDO
	  ENDDO
	  N_OFFGRID_USERTAUS = UT

C  set remaining optical depth mask (for values at layer boundaries)

	  NSTART = 0
	  DO UTA = 1, N_OUT_USERTAUS
	    IF ( .NOT. OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
	      TAU = USER_TAUS_INPUT(UTA)
	      LOOP = .TRUE.
	      N = NSTART
	      DO WHILE (LOOP.AND.N.LE.NLAYERS)
	        N1 = N
	        N = N + 1
	        IF ( DABS(TAU-TAUGRID_INPUT(N1)).LT.SMALLNUM ) THEN
	          LOOP = .FALSE.
	          UTAU_LEVEL_MASK_UP(UTA) = N1
	          UTAU_LEVEL_MASK_DN(UTA) = N1
	          NSTART = N1
	        ENDIF
	      ENDDO
	    ENDIF
	  ENDDO

C  Otherwise if the optical depths are set to be at all layer boundaries

	ELSE IF ( DO_LBOUND_TAUS ) THEN

C  assign

	  N_OFFGRID_USERTAUS = 0
	  DO UTA = 1, N_OUT_USERTAUS
	    USER_TAUS_INPUT(UTA) = TAUGRID_INPUT(LBOUND_TAUS_INPUT(UTA))
	    OFFGRID_UTAU_OUTFLAG(UTA) = .FALSE.
	  ENDDO

C  Sort in ascending order

	  IF ( N_OUT_USERTAUS .GT. 1 ) THEN
	    CALL HPSORT(N_OUT_USERTAUS,USER_TAUS_INPUT)
	  ENDIF

C  set optical depth mask (for values at layer boundaries)

	  NSTART = 0
	  DO UTA = 1, N_OUT_USERTAUS
	    TAU = USER_TAUS_INPUT(UTA)
	    LOOP = .TRUE.
	    N = NSTART
	    DO WHILE (LOOP.AND.N.LE.NLAYERS)
	      N1 = N
	      N = N + 1
	      IF ( DABS(TAU-TAUGRID_INPUT(N1)).LT.SMALLNUM ) THEN
	        LOOP = .FALSE.
	        UTAU_LEVEL_MASK_UP(UTA) = N1
	        UTAU_LEVEL_MASK_DN(UTA) = N1
	        NSTART = N1
	      ENDIF
	    ENDDO
	  ENDDO

	ENDIF

C  Set masking and number of layer source terms
C  --------------------------------------------

C   .. for upwelling

	IF ( DO_UPWELLING ) THEN
	  DO N = 1, NLAYERS
	    STERM_LAYERMASK_UP(N) = .FALSE.
	  ENDDO
	  UTA = 1
	  UT  = 1
	  IF ( .NOT. OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
	    N_LAYERSOURCE_UP = UTAU_LEVEL_MASK_UP(UTA) + 1
	    N_ALLLAYERS_UP   = N_LAYERSOURCE_UP
	  ELSE
	    N_LAYERSOURCE_UP = OFFGRID_UTAU_LAYERIDX(UT) + 1
	    N_ALLLAYERS_UP   = N_LAYERSOURCE_UP - 1
	  ENDIF
	  DO N = NLAYERS, N_ALLLAYERS_UP, -1
	    STERM_LAYERMASK_UP(N) = .TRUE.
	  ENDDO
	ENDIF

C   .. for downwelling

	IF ( DO_DNWELLING ) THEN
	  DO N = 1, NLAYERS
	    STERM_LAYERMASK_DN(N) = .FALSE.
	  ENDDO
	  UTA = N_OUT_USERTAUS
	  UT  = N_OFFGRID_USERTAUS
	  IF ( .NOT. OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
	    N_LAYERSOURCE_DN = UTAU_LEVEL_MASK_DN(UTA)
	    N_ALLLAYERS_DN   = N_LAYERSOURCE_DN
	  ELSE
	    N_LAYERSOURCE_DN = OFFGRID_UTAU_LAYERIDX(UT)
	    N_ALLLAYERS_DN   = N_LAYERSOURCE_DN
	  ENDIF
	  DO N = 1, N_ALLLAYERS_DN
	    STERM_LAYERMASK_DN(N) = .TRUE.
	  ENDDO
	ENDIF

C  Finish

	RETURN
	END
