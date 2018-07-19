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

	SUBROUTINE LIDORT_READINPUT ( STATUS )

C  Read all control inputs for LIDORT
C  ----------------------------------

C  include file of constants

	INCLUDE '../include_s/LIDORT.PARS'

C  include files to be filled out with input data

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  Module arguments (output status)

	INTEGER		STATUS

C  local variables

	CHARACTER*8	PREFIX
	PARAMETER	( PREFIX = 'LIDORT -' )
	LOGICAL		ERROR
        CHARACTER * 80  PAR_STR
        LOGICAL         GFINDPAR
        EXTERNAL        GFINDPAR
	INTEGER		I, FILUNIT, LEN_STRING
	CHARACTER*70	MAIL
	EXTERNAL	LEN_STRING

C  initialize status

	STATUS = LIDORT_SUCCESS
	ERROR  = .FALSE.

C  file unit

	FILUNIT = LIDORT_INUNIT

C  first read all variables in include file LIDORT_CONTROL.VARS
C  ============================================================

C  operation modes

        PAR_STR = 'Do full radiance calculation?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_FULLRAD_MODE
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Do single scatter correction?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_SSCORRECTION
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  multiple scatter source function output control

        PAR_STR = 'Output multiple scatter layer source functions?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) SAVE_LAYER_MSST
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  double convergence test

        PAR_STR = 'Perform double convergence test?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DOUBLE_CONV_TEST
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  particular solution control

        PAR_STR = 'Include direct beam?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_DIRECT_BEAM
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Use classical beam solution?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_CLASSICAL_SOLUTION
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Quasi-spherical treatment of direct beam?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_QUASPHER_BEAM
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Perform internal Chapman function calculation?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &     READ (FILUNIT,*,ERR=998) DO_CHAPMAN_FUNCTION
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

	IF ( DO_QUASPHER_BEAM ) THEN

          PAR_STR = 'Beam path with refractive atmosphere?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &       READ (FILUNIT,*,ERR=998) DO_QSREFRAC_BEAM
	  CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C          PAR_STR = 'Use coefficient expansion in beam solution?'
C          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
C     &       READ (FILUNIT,*,ERR=998) DO_USE_QSCOEFFS
C	  CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

	ENDIF

C  scatterers and phase function control

        PAR_STR='Rayleigh atmosphere only?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_RAYLEIGH_ONLY
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR='Isotropic atmosphere only?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_ISOTROPIC_ONLY
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Surface control

        PAR_STR = 'Lambertian albedo?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_LAMBERTIAN_ALBEDO
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Include surface emission?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_SURFACE_EMISSION
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  numerical control

        PAR_STR = 'No azimuth dependence in the solution?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_NO_AZIMUTH
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Use full -1 to +1 quadrature?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_FULL_QUADRATURE
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Compute all Fourier components?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_ALL_FOURIER
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Include Delta-M scaling?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_DELTAM_SCALING
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  mean-value output control

        PAR_STR = 'Generate mean value output additionally?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_ADDITIONAL_MVOUT
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Generate only mean value output?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_MVOUT_ONLY
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  directional output control

        PAR_STR = 'Upwelling output?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_UPWELLING
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Downwelling output?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_DNWELLING
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Output file control

        PAR_STR = 'Debug write?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_DEBUG_WRITE
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Input control write?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_WRITE_INPUT
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Input scenario write?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_WRITE_SCENARIO
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Results write?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_WRITE_RESULTS
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Results write IDL output?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_WRITE_RESULTS_IDL
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Fourier component output write?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_WRITE_FOURIER
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Stream angle and optical depth output control

        PAR_STR = 'Include quadrature angles in output?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_QUAD_OUTPUT
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'User-defined stream angles?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &   READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'User-defined optical depths?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &   READ (FILUNIT,*,ERR=998) DO_USER_TAUS
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Layer boundary optical depths?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &   READ (FILUNIT,*,ERR=998) DO_LBOUND_TAUS
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  output filenames

        PAR_STR = 'filename for main output'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,'(a)',ERR=998) RESULTS_WRITE_FILENAME
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'filename for main IDL output'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,'(a)',ERR=998) IDLRESULTS_WRITE_FILENAME
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'filename for input write'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,'(a)',ERR=998) INPUT_WRITE_FILENAME
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'filename for scenario write'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,'(a)',ERR=998) SCENARIO_WRITE_FILENAME
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Fourier output filename'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,'(a)',ERR=998) FOURIER_WRITE_FILENAME
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Now read variables in LIDORT_MODEL.VARS
C  =======================================

C  Stream/layer/moment (INTEGER input)

        PAR_STR = 'Number of half-space streams'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) NSTREAMS
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Number of atmospheric layers'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) NLAYERS
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Number of input Legendre moments'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) NMOMENTS_INPUT
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  All numbers are now checked against maximum diensions

	IF ( NSTREAMS .GT. MAXSTRM ) THEN
	  MAIL = 'Number of half-space streams > maximum dimension'
	  STATUS  = LIDORT_SERIOUS
	  CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
	ENDIF
	IF ( NLAYERS .GT. MAXLAYER ) THEN
	  MAIL = 'Number of layers > maximum dimension'
	  STATUS  = LIDORT_SERIOUS
	  CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
	ENDIF
	IF ( NMOMENTS_INPUT .GT. MAXMOMENT ) THEN
	  MAIL = 'Number of Legendre moments > maximum dimension'
	  STATUS  = LIDORT_SERIOUS
	  CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
	ENDIF

C  Flux/accuracy input

        PAR_STR = 'Solar flux constant'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) FLUX_FACTOR
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Fourier series convergence'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) LIDORT_ACCURACY
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Zenith tolerance level'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) ZENITH_TOLERANCE
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  beam model input

        PAR_STR = 'Solar zenith angle (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) SUN0
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  Earth radius (only for Chapman function calculation)

	IF ( DO_CHAPMAN_FUNCTION ) THEN
	  PAR_STR = 'Earth radius (km)'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &       READ (FILUNIT,*,ERR=998) EARTH_RADIUS
	  CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )
	ENDIF

C  User defined input
C      - Check numbers do not exceed specified dimensions

C  1. always need some azimuth angles if Azimuth flag set

	IF ( .NOT. DO_NO_AZIMUTH ) THEN

          PAR_STR = 'Number of user-defined relative azimuth angles'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
	  CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

	  IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
	    MAIL = 'Number of relative azimuths > maximum dimension'
	    STATUS  = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
	  ENDIF

          PAR_STR = 'User-defined relative azimuth angles'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
	    DO I = 1, N_USER_RELAZMS
	      READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
	    ENDDO
 	  ENDIF
	  CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

	ENDIF

C  2. User defined stream angles (should be positive)

	IF ( DO_USER_STREAMS ) THEN

          PAR_STR = 'Number of user-defined stream angles'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) N_USER_STREAMS
	  CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

	  IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
	    MAIL = 'Number of user streams > maximum dimension'
	    STATUS  = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
	  ENDIF

          PAR_STR = 'User-defined stream angles'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
	    DO I = 1, N_USER_STREAMS
	      READ (FILUNIT,*,ERR=998) USER_ANGLES_INPUT(I)
	    ENDDO
 	  ENDIF
	  CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

	ENDIF

C  3. User defined optical depths (should be positive)

	IF ( DO_USER_TAUS ) THEN

          PAR_STR = 'Number of user-defined optical depths'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &        READ (FILUNIT,*,ERR=998) N_OUT_USERTAUS
	  CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

	  IF ( N_OUT_USERTAUS .GT. MAX_OUT_USERTAUS ) THEN
	    MAIL=
     &      'Number of user-defined optical depths > maximum dimension'
	    STATUS  = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
	  ENDIF

          PAR_STR = 'User-defined optical depths'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
	    DO I = 1, N_OUT_USERTAUS
	      READ (FILUNIT,*,ERR=998) USER_TAUS_INPUT(I)
	    ENDDO
 	  ENDIF
	  CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C   ( override if layer boundaries are required )

	ELSE IF ( DO_LBOUND_TAUS ) THEN

          PAR_STR = 'Number of layer boundary optical depths'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &        READ (FILUNIT,*,ERR=998) N_OUT_USERTAUS
	  CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

	  IF ( N_OUT_USERTAUS .GT. MAX_OUT_USERTAUS ) THEN
	    MAIL=
     &    'Number of layer boundary optical depths > maximum dimension'
	    STATUS  = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
	  ENDIF

          PAR_STR = 'Which layer boundaries for optical depth'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
	    DO I = 1, N_OUT_USERTAUS
	      READ (FILUNIT,*,ERR=998) LBOUND_TAUS_INPUT(I)
	    ENDDO
 	  ENDIF
	  CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

	ENDIF

C  normal return

	RETURN

C  line read error - abort immediately

998	CONTINUE
	STATUS = LIDORT_SERIOUS
	MAIL = 'read failure for '//PAR_STR(1:LEN_STRING(PAR_STR))
	CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
	RETURN

C  Finish

	END

C

	SUBROUTINE LIDORT_INIT_CONTROL_VARS

C  Initialises all control inputs for LIDORT
C  -----------------------------------------

C  include file of constants

	INCLUDE '../include_s/LIDORT.PARS'

C  include files to be initialised

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'

	DO_FULLRAD_MODE = .FALSE.
	DO_SSCORRECTION = .FALSE.
	SAVE_LAYER_MSST = .FALSE.
	DOUBLE_CONV_TEST = .FALSE.

	DO_DIRECT_BEAM = .FALSE.
	DO_CLASSICAL_SOLUTION = .FALSE.

	DO_QUASPHER_BEAM = .FALSE.
	DO_QSREFRAC_BEAM = .FALSE.
	DO_CHAPMAN_FUNCTION = .FALSE.

	DO_RAYLEIGH_ONLY = .FALSE.
	DO_ISOTROPIC_ONLY = .FALSE.

	DO_LAMBERTIAN_ALBEDO = .FALSE.
	DO_SURFACE_EMISSION = .FALSE.

	DO_NO_AZIMUTH = .FALSE.
	DO_FULL_QUADRATURE = .FALSE.
	DO_ALL_FOURIER = .FALSE.
	DO_DELTAM_SCALING = .FALSE.

	DO_ADDITIONAL_MVOUT = .FALSE.
	DO_MVOUT_ONLY = .FALSE.

	DO_DEBUG_WRITE = .FALSE.
	DO_WRITE_INPUT = .FALSE.
	DO_WRITE_RESULTS = .FALSE.
	DO_WRITE_RESULTS_IDL = .FALSE.
	DO_WRITE_SCENARIO = .FALSE.
	DO_WRITE_FOURIER = .FALSE.

	DO_QUAD_OUTPUT = .FALSE.
	DO_USER_STREAMS = .FALSE.
	DO_USER_TAUS = .FALSE.

	DO_UPWELLING = .FALSE.
	DO_DNWELLING = .FALSE.

	END

C

	SUBROUTINE LIDORT_INIT_MODEL_VARS

C  Initialises all file-read model inputs for LIDORT
C  -------------------------------------------------

C  include file of constants

	INCLUDE '../include_s/LIDORT.PARS'

C  include files to be initialised

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  local variables

	INTEGER		 I

	FLUX_FACTOR = ZERO
	LIDORT_ACCURACY = ZERO
	SUN0 = ZERO
	EARTH_RADIUS = ZERO
	ZENITH_TOLERANCE = ZERO

	N_USER_STREAMS = 0
	N_USER_RELAZMS = 0
	N_OUT_USERTAUS = 0
	NSTREAMS = 0
	NLAYERS = 0
	NMOMENTS_INPUT = 0

	DO I = 1, MAX_USER_STREAMS
	  USER_ANGLES_INPUT(I) = ZERO
	ENDDO
	DO I = 1, MAX_USER_RELAZMS
	  USER_RELAZMS(I) = ZERO
	ENDDO
	DO I = 1, MAX_OUT_USERTAUS
	  USER_TAUS_INPUT(I) = ZERO
	ENDDO

	END

C

	SUBROUTINE FINDPAR_ERROR
     &     (  ERROR, PAR_STR, STATUS )

C  include file

	INCLUDE '../include_s/LIDORT.PARS'

C  subroutine arguments

	LOGICAL		ERROR
	CHARACTER*(*)	PAR_STR
	INTEGER		STATUS
	
C  local variables

	INTEGER		LEN_STRING
	CHARACTER*70	MAIL
	EXTERNAL	LEN_STRING

	IF ( ERROR ) THEN
	  STATUS = LIDORT_SERIOUS
	  MAIL = 'Cannot find string: '//PAR_STR(1:LEN_STRING(PAR_STR))
	  CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
	ENDIF

	RETURN
	END
