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

	SUBROUTINE LIDORT_V23E_FOURIER
     I       ( FOURIER,
     O         DO_INCLUDE_ALBEDO, DO_INCLUDE_SURFEMISS, STATUS )

C  Complete Fourier component calculation for the Extended Code.

C  Standard include files
C  ----------------------

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  Include file of setup variables

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include files of solution variables

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include files of multiplier variables (output stored here)

	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  Extended Include files
C  ----------------------

C  include file of dimensions

	INCLUDE '../include_e/LIDORT_L.PARS'

C  include files of linearized control variables (input)

	INCLUDE '../include_e/LIDORT_L_CONTROL.VARS'

C  include files of linearized geophysical variables (input)

	INCLUDE '../include_e/LIDORT_L_GEOPHYS.VARS'

C  include files of linearized solution variables (output stored here)

	INCLUDE '../include_e/LIDORT_L_SOLUTION.VARS'


C  module arguments
C  ----------------

C  Input

	INTEGER		 FOURIER

C  Output

	LOGICAL		 DO_INCLUDE_ALBEDO
	LOGICAL		 DO_INCLUDE_SURFEMISS
	INTEGER		 STATUS

C  Local variables
C  ---------------

C  Green function multiplier flag (Quadrature solutions only)

	LOGICAL		 DO_GMULT(MAX_OFFGRID_USERTAUS)

C  local inclusion flags

	LOGICAL		 DO_INCLUDE_DIRECTBEAM
C	LOGICAL		 DO_INCLUDE_THERMAL		! REMOVED
	LOGICAL		 DO_INCLUDE_MVOUTPUT

C  Direct beam variables

	DOUBLE PRECISION ALB_X0_FLUX, ATTN
	DOUBLE PRECISION DIRECT_BEAM(MAXSTRM)
	DOUBLE PRECISION USER_DIRECT_BEAM(MAX_USER_STREAMS)

C  Other local help variables 

	INTEGER		 I, I1, UM, LAYER, Q, C0, N, UT, AA
	DOUBLE PRECISION FLUXFAC,DELFAC,R2,TAUP,SSFLUX

	INTEGER          N_PARAMETERS, LAYER_TO_VARY
	LOGICAL		 MODIFIED_BCL3, MODIFIED_BCL4

C  error tracing variables

	INTEGER          INFO
	CHARACTER*3	 CN, CI
	CHARACTER*70	 MAIL, TRACE
	INTEGER		 STATUS_SUB

C  ##############
C  initialization
C  ##############

C  module status

	STATUS = LIDORT_SUCCESS

C  Set local flags
C  ---------------

C  inclusion of thermal surface emission term, only for Fourier = 0

	DO_INCLUDE_SURFEMISS = .FALSE.
	IF ( DO_SURFACE_EMISSION ) THEN
	  IF ( FOURIER .EQ. 0 ) THEN
	    DO_INCLUDE_SURFEMISS = .TRUE.
	  ENDIF
	ENDIF

C  Albedo flag (for inclusion of some kind of reflecting boundary)

	DO_INCLUDE_ALBEDO = .TRUE.
	IF ( DO_LAMBERTIAN_ALBEDO ) THEN
	  IF ( FOURIER .NE. 0 ) DO_INCLUDE_ALBEDO = .FALSE.
	ENDIF

C  Direct beam flag (only if above albedo flag has been set)

	DO_INCLUDE_DIRECTBEAM = .FALSE.
	IF ( DO_DIRECT_BEAM ) THEN
	  IF ( DO_INCLUDE_ALBEDO ) THEN
	    DO_INCLUDE_DIRECTBEAM = .TRUE.
	  ENDIF
	ENDIF

C  Flux factor FLUXFAC

	DELFAC = TWO
	IF ( FOURIER .EQ. 0 ) DELFAC = ONE
	SSFLUX  = QUARTER*FLUX_FACTOR/PIE
	FLUXFAC = SSFLUX * DELFAC

C  albedo reflectance factor R2

	R2 = ALBEDO
	IF ( FOURIER .EQ. 0 ) R2 = TWO * ALBEDO

C  inclusion of thermal term (REMOVED)
C	DO_INCLUDE_THERMAL = .FALSE.
C	IF ( DO_THERMAL_EMISSION ) THEN
C	  IF ( FOURIER .EQ. 0 ) THEN
C	    DO_INCLUDE_THERMAL = .TRUE.
C	  ENDIF
C	ENDIF

C  inclusion of mean value output

	DO_INCLUDE_MVOUTPUT = .FALSE.
	IF ( DO_ADDITIONAL_MVOUT. OR. DO_MVOUT_ONLY ) THEN
	  IF ( FOURIER .EQ. 0 ) THEN
	    DO_INCLUDE_MVOUTPUT = .TRUE.
	  ENDIF
	ENDIF

C  Miscellaneous setups for Fourier = 0
C  ------------------------------------

	IF ( FOURIER .EQ. 0 ) THEN

C  standard setups

	  CALL LIDORT_MISCSETUPS

C  Additional setups for the linearizations

	  IF ( DO_LAYER_LINEARIZATION ) THEN
	    CALL LIDORT_L_MISCSETUPS
	  ENDIF

	ENDIF

C  Reflected Direct beam attenuation 
C  ---------------------------------

C  (division by FLUXFAC comes later)

	IF ( DO_INCLUDE_DIRECTBEAM ) THEN

C  Attenuation of beam

	  ALB_X0_FLUX = FOUR * ALBEDO * X0 * FLUX_FACTOR
	  IF ( DO_QUASPHER_BEAM ) THEN
	    TAUP = TAUSLANT(NLAYERS)
	  ELSE
	    TAUP = TAUGRID(NLAYERS)/X0
	  ENDIF
	  IF ( TAUP .GT. MAX_TAU_SPATH ) THEN
	    ATTN = ZERO
	  ELSE
	    ATTN = ALB_X0_FLUX * DEXP(-TAUP)
	  ENDIF

C  reflected along the quadrature directions

	  DO I = 1, NSTREAMS
	    DIRECT_BEAM(I) = ATTN * BIREFLEC_0(FOURIER,I)
	  ENDDO

C  reflected along user_stream directions

	  IF ( DO_USER_STREAMS ) THEN
	    DO UM = 1, N_USER_STREAMS
	      USER_DIRECT_BEAM(UM) = ATTN*USER_BIREFLEC_0(FOURIER,UM)
	    ENDDO
	  ENDIF

C  end direct beam calculation

	ENDIF

C  ###################################################
C  RT differential equation solutions + linearizations
C  ###################################################

C  Get Legendre polynomials for this Fourier component

	CALL LIDORT_LEGENDRE_SETUP ( FOURIER )
	IF ( DO_USER_STREAMS ) THEN
	  CALL LIDORT_USERLEGENDRE_SETUP ( FOURIER )
	ENDIF

C  Start layer loop

	DO LAYER = 1, NLAYERS

C  Standard Solutions
C  ------------------

C  Get Discrete ordinate solutions for this layer

	  CALL LIDORT_RTSOLUTION
     I    ( LAYER, FOURIER,
     O      STATUS_SUB )

C  .. error tracing

	  IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
	    MAIL = 'Error return from module LIDORT_RTSOLUTION'
	    TRACE= 'First Call in LIDORT_V23E_FOURIER'
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
	    RETURN
	  ENDIF

C  Get Post-processing ("user") solutions for this layer

	  CALL LIDORT_USERSOLUTION
     I    ( LAYER, FOURIER )

C  debug

C	write(99,*)'layer xpos',layer,FOURIER
C	write(99,'(a2,10(1pe14.6))')'EV',(keigen(k,layer),k=1,nstreams)
C	do i = 1, nstr2
C	write(99,'(i2,12(1pe14.6))')i,(xpos(i,k,layer),k=1,nstreams),
C     &     wupper(i,layer),wlower(i,layer)
C	enddo
C	write(99,*)'layer u_xpos',layer
CC	do i = 1, n_user_streams
C	write(99,'(i2,10(1pe14.6))')i,(u_xpos(i,k,layer),k=1,nstreams),
C     &     u_wpos(i,layer),u_wneg(i,layer)
C	enddo

C  Additional Solutions for Linearization
C  --------------------------------------

	  IF ( DO_LAYER_LINEARIZATION ) THEN

C  Linearizations of Discrete Ordinate solution

	    CALL LIDORT_L_RTSOLUTION
     I      ( LAYER, FOURIER,
     I        LAYER_VARY_FLAG(LAYER),
     I        LAYER_VARY_NUMBER(LAYER),
     O        STATUS_SUB )

C  .. error tracing

	    IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
	      MAIL = 'Error return from module LIDORT_L_RTSOLUTION'
	      TRACE= 'Second Call in LIDORT_V23E_FOURIER'
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
	      RETURN
	    ENDIF

C  Get Post-processing ("user") solutions for this layer

	    CALL LIDORT_L_USERSOLUTION
     I      ( LAYER, FOURIER,
     I        LAYER_VARY_FLAG(LAYER),
     I        LAYER_VARY_NUMBER(LAYER) )

C  End linearization control

	  ENDIF

C  end layer loop

	ENDDO

C  ######################
C  Intensity calculations
C  ######################

C  Standard boundary value problem
C  -------------------------------

	CALL LIDORT_BVPROBLEM
     I     ( DO_INCLUDE_ALBEDO,
     I       DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS,
     I       FOURIER, R2, DIRECT_BEAM,
     O       STATUS_SUB )

C  .. error tracing

	IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
	  MAIL = 'Error return from module LIDORT_BVPROBLEM'
	  TRACE= 'Third Call in LIDORT_V23E_FOURIER'
	  STATUS = LIDORT_SERIOUS
	  CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
	  RETURN
	ENDIF

C  Intensity
C  ---------

C  Do the original intensity calculations first, and some of linearizations
C  required BEFORE the perturbed boundary value problems are solved

C  initialise multiplier flag

	DO UT = 1, N_OFFGRID_USERTAUS
	  DO_GMULT(UT) = .TRUE.
	ENDDO

C  Upwelling Intensity

	IF ( DO_UPWELLING ) THEN
	  CALL UPUSER_INTENSITY
     I    ( DO_INCLUDE_ALBEDO,   DO_INCLUDE_SURFEMISS,
     I      DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTBEAM,
     I      R2, FLUXFAC, FOURIER,
     I      USER_DIRECT_BEAM, DO_GMULT )
	ENDIF

C  Downwelling Intensity

	IF ( DO_DNWELLING ) THEN
	  CALL DNUSER_INTENSITY
     I    ( DO_INCLUDE_MVOUTPUT, FLUXFAC, DO_GMULT )
	ENDIF

C  mean value output

	IF ( DO_INCLUDE_MVOUTPUT ) THEN
	  CALL LIDORT_MIFLUX_INTENSITY
	ENDIF

C  ###############################
C  Weighting Function calculations
C  ###############################

C  Finished if only intensity is required

	IF ( DO_SIMULATION_ONLY ) RETURN

C  Atmospheric weighting functions
C  ===============================

	IF ( DO_LAYER_LINEARIZATION ) THEN

C  Start over layers that will contain variations

	  DO LAYER_TO_VARY = 1, NLAYERS

C  flag for variation of layer

	    IF ( LAYER_VARY_FLAG(LAYER_TO_VARY) ) THEN

C  number of variations for this layer

	      N_PARAMETERS = LAYER_VARY_NUMBER(LAYER_TO_VARY)

C  Boundary condition flags for special cases

	      MODIFIED_BCL3 = ( LAYER_TO_VARY .EQ. 1 )
	      MODIFIED_BCL4 = ( LAYER_TO_VARY .EQ. NLAYERS )

C  Compute the main column B' where AX = B'

	      CALL LIDORT_ATMOSWF_COLUMN_SETUP
     I     ( DO_INCLUDE_ALBEDO, DO_INCLUDE_DIRECTBEAM,
     I       MODIFIED_BCL3,     MODIFIED_BCL4,
     I       LAYER_TO_VARY, N_PARAMETERS, FOURIER,
     I       R2, DIRECT_BEAM )

C  BV solution for perturbed integration constants
C  -----------------------------------------------

C  call to LAPACK solver routine for back substitution

              CALL DGBTRS
     &     ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_PARAMETERS,
     &        BANDMAT2, MAXBANDTOTAL, IPIVOT,
     &        COL2_WF, MAXTOTAL, INFO )

	      IF ( INFO .LT. 0 ) THEN
	        WRITE(CI, '(I3)' ) INFO
	        WRITE(CN, '(I3)' ) LAYER_TO_VARY
	        MAIL  = 'argument i illegal value, for i = '//CI
	        TRACE = 'Atmos_Wfs for layer '//CN//
     *                  '; DGBTRS call in LIDORT_V23E_FOURIER'
	        STATUS = LIDORT_SERIOUS
	        CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
	        RETURN
	      ENDIF

C  Set integration constants NCON and PCON for +/- eigensolutions

	      DO N = 1, NLAYERS
	        C0 = (N-1)*NSTR2
	        DO I = 1, NSTREAMS
	          I1 = I + NSTREAMS
	          DO Q = 1, N_PARAMETERS
	            NCON(I,N,Q) = COL2_WF(C0+I,Q)
	            PCON(I,N,Q) = COL2_WF(C0+I1,Q)
	          ENDDO
	        ENDDO
	      ENDDO

C  Associated quantities

	      DO N = 1, NLAYERS
	        DO I = 1, NSTR2
	          DO AA = 1, NSTREAMS
	            DO Q = 1, N_PARAMETERS
	              NCON_XVEC(I,AA,N,Q) = NCON(AA,N,Q) * XPOS(I,AA,N)
	              PCON_XVEC(I,AA,N,Q) = PCON(AA,N,Q) * XNEG(I,AA,N)
	            ENDDO
	          ENDDO
	        ENDDO
	      ENDDO

C  Get the weighting functions
C  ---------------------------

C  initialise multiplier flag

	      DO UT = 1, N_OFFGRID_USERTAUS
	        DO_GMULT(UT) = .TRUE.
	      ENDDO

C  Upwelling Atmospheric weighting functions

	      IF ( DO_UPWELLING ) THEN
	        CALL UPUSER_ATMOSWF
     I           ( DO_INCLUDE_ALBEDO,
     I             DO_INCLUDE_MVOUTPUT,
     I             DO_INCLUDE_DIRECTBEAM,
     I             R2, FLUXFAC, FOURIER,
     I             LAYER_TO_VARY, N_PARAMETERS,
     I             USER_DIRECT_BEAM, DO_GMULT )
	      ENDIF

C  Downwelling Atmospheric weighting functions

	      IF ( DO_DNWELLING ) THEN
	        CALL DNUSER_ATMOSWF
     I             ( FOURIER, FLUXFAC, DO_GMULT, DO_INCLUDE_MVOUTPUT,
     I               LAYER_TO_VARY, N_PARAMETERS )
	      ENDIF

C  mean value output

	       IF ( DO_INCLUDE_MVOUTPUT ) THEN
	         CALL LIDORT_MIFLUX_ATMOSWF
     I               ( LAYER_TO_VARY, N_PARAMETERS )
	       ENDIF

C  Finish loop over layers with variation

	    ENDIF
	  ENDDO
	ENDIF

C  Albedo weighting functions
C  ==========================

	IF ( DO_INCLUDE_ALBEDO ) THEN
	  IF ( DO_ALBEDO_LINEARIZATION ) THEN

C  BV solution for perturbed integration constants
C  -----------------------------------------------

C  Compute the main column B' where AX = B'

	    CALL LIDORT_ALBEDOWF_COLUMN_SETUP
     I   (  DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS, DIRECT_BEAM )

C  call to LAPACK solver routine (back substitution)

            CALL DGBTRS
     &     ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1,
     &        BANDMAT2, MAXBANDTOTAL, IPIVOT,
     &        COL2_WFALB, MAXTOTAL, INFO )

	    IF ( INFO .LT. 0 ) THEN
	      WRITE(CI, '(I3)' ) INFO
	      MAIL  = 'argument i illegal value, for i = '//CI
	      TRACE = 'Wfs for albedo'//
     *                '; DGBTRS call in LIDORT_V23E_FOURIER'
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
	      RETURN
	    ENDIF

C  Set all integration constants NCON_ALB and PCON_ALB for +/- eigensols

	    DO N = 1, NLAYERS
	      C0 = (N-1)*NSTR2
	      DO I = 1, NSTREAMS
	        I1 = I+NSTREAMS
	        NCON_ALB(I,N) = COL2_WFALB(C0+I,1)
	        PCON_ALB(I,N) = COL2_WFALB(C0+I1,1)
	      ENDDO
	    ENDDO

C  Get the weighting functions
C  ---------------------------

C  Upwelling Albedo weighting functions

	    IF ( DO_UPWELLING ) THEN
	      CALL UPUSER_ALBEDOWF
     I         ( DO_INCLUDE_SURFEMISS, DO_INCLUDE_MVOUTPUT,
     I           R2, FLUXFAC, FOURIER, USER_DIRECT_BEAM )
	    ENDIF

C  Downwelling Albedo weighting functions

	    IF ( DO_DNWELLING ) THEN
	        CALL DNUSER_ALBEDOWF ( DO_INCLUDE_MVOUTPUT, FLUXFAC )
	    ENDIF

C  mean value output

	     IF ( DO_INCLUDE_MVOUTPUT ) THEN
	       CALL LIDORT_MIFLUX_ALBEDOWF
	     ENDIF

C  end of albedo weighting functions

	  ENDIF
	ENDIF

C  SurfBB weighting functions
C  ==========================

	IF ( DO_INCLUDE_SURFEMISS ) THEN
	  IF ( DO_SURFBB_LINEARIZATION ) THEN

C  BV solution for perturbed integration constants
C  -----------------------------------------------

C  Compute the main column B' where AX = B'

	    CALL LIDORT_SURFBBWF_COLUMN_SETUP

C  call to LAPACK solver routine (back substitution)

            CALL DGBTRS
     &     ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1,
     &        BANDMAT2, MAXBANDTOTAL, IPIVOT,
     &        COL2_WFSBB, MAXTOTAL, INFO )

	    IF ( INFO .LT. 0 ) THEN
	      WRITE(CI, '(I3)' ) INFO
	      MAIL  = 'argument i illegal value, for i = '//CI
	      TRACE = 'Wfs for surfbb'//
     *                '; DGBTRS call in LIDORT_V23E_FOURIER'
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
	      RETURN
	    ENDIF

C  Set all integration constants NCON_SURF and PCON_SURF for +/- eigensols

	    DO N = 1, NLAYERS
	      C0 = (N-1)*NSTR2
	      DO I = 1, NSTREAMS
	        I1 = I+NSTREAMS
	        NCON_SURF(I,N) = COL2_WFSBB(C0+I,1)
	        PCON_SURF(I,N) = COL2_WFSBB(C0+I1,1)
	      ENDDO
	    ENDDO

C  Get the weighting functions
C  ---------------------------

C  Upwelling Surfbb weighting functions

	    IF ( DO_UPWELLING ) THEN
	      CALL UPUSER_SURFBBWF
     I           ( DO_INCLUDE_ALBEDO, DO_INCLUDE_MVOUTPUT,
     I             R2, FLUXFAC, FOURIER )
	    ENDIF

C  Downwelling Albedo weighting functions

	    IF ( DO_DNWELLING ) THEN
	        CALL DNUSER_SURFBBWF ( DO_INCLUDE_MVOUTPUT, FLUXFAC )
	    ENDIF

C  mean value output

	     IF ( DO_INCLUDE_MVOUTPUT ) THEN
	       CALL LIDORT_MIFLUX_SURFBBWF
	     ENDIF

C  end of surfbb Weighting function

	  ENDIF
	ENDIF

C  ######
C  finish
C  ######

	RETURN
	END
