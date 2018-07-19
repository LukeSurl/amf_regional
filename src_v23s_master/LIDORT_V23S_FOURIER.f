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

	SUBROUTINE LIDORT_V23S_FOURIER
     I       ( FOURIER,
     O         DO_INCLUDE_ALBEDO, DO_INCLUDE_SURFEMISS, STATUS )

C  Complete Fourier component calculation for the Standard Code.

C  include files
C  -------------

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  Include file of setup variables (output from SETUPS routine)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include file for debug

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

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

	INTEGER		 I, UM, LAYER, UT
	DOUBLE PRECISION FLUXFAC,DELFAC,R2,TAUP,SSFLUX

C  error tracing

	CHARACTER*(70)	 MAIL, TRACE
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

	  CALL LIDORT_MISCSETUPS

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

C  ##################################
C  RT differential equation solutions
C  ##################################

C  Get Legendre polynomials for this Fourier component

	CALL LIDORT_LEGENDRE_SETUP ( FOURIER )
	IF ( DO_USER_STREAMS ) THEN
	  CALL LIDORT_USERLEGENDRE_SETUP ( FOURIER )
	ENDIF

C  Start layer loop

	DO LAYER = 1, NLAYERS

C  Get Discrete ordinate solutions for this layer

	  CALL LIDORT_RTSOLUTION
     I    ( LAYER, FOURIER,
     O      STATUS_SUB )

C  .. error tracing

	  IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
	    MAIL = 'Error return from module LIDORT_RTSOLUTION'
	    TRACE= 'First Call in LIDORT_V23S_FOURIER'
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
	    RETURN
	  ENDIF

C  Get Post-processing ("user") solutions for this layer

	  CALL LIDORT_USERSOLUTION
     I    ( LAYER, FOURIER )

C  debug

C	if (layer.eq.1)write(99,*)
C     &  DO_QUASPHER_BEAM, DO_QSREFRAC_BEAM, DO_CLASSICAL_SOLUTION
C	if (layer.eq.1)write(99,*)X0
C	write(99,*)'layer xpos',layer
C	write(99,'(a2,8(1pe14.6))')'EV',(keigen(k,layer),k=1,nstreams)
C	do i = 1, nstr2
C	write(99,'(i2,10(1pe14.6))')i,(xpos(i,k,layer),k=1,nstreams),
C     &     wupper(i,layer),wlower(i,layer)
C	enddo
C	write(99,*)'layer u_xpos',layer
C	do i = 1, n_user_streams
C	write(99,'(i2,10(1pe14.6))')i,(u_xpos(i,k,layer),k=1,nstreams),
C     &     u_wpos(i,layer),u_wneg(i,layer)
C	enddo

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
	  TRACE= 'Second Call in LIDORT_V23S_FOURIER'
	  STATUS = LIDORT_SERIOUS
	  CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
	  RETURN
	ENDIF

C  Intensity
C  ---------

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

C  ######
C  finish
C  ######

	RETURN
	END
