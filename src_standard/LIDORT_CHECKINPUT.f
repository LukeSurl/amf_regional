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
C #  Release Date   Janaury 2001				#
C #		  					        #
C ###############################################################

	SUBROUTINE LIDORT_CHECK_INPUT (STATUS)

C  check inputs (both file-read and derived)

	INCLUDE '../include_s/LIDORT.PARS'

C  include files with input variables to be checked out

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'

C  Module output

	INTEGER		 STATUS

C  local variables

	CHARACTER*70	MAIL
	CHARACTER*70	ACTION

	INTEGER		 I, L, N, UTA, NSTART, N1
	DOUBLE PRECISION TAU
	CHARACTER*2	 C2
	CHARACTER*3	 C3
	LOGICAL		 LOOP

C  Initialize output status

	STATUS = LIDORT_SUCCESS

C  Do some checking
C  ================

C  Check that the Chapman function will not be done for refractive geometry

	IF ( DO_CHAPMAN_FUNCTION ) THEN
	  IF ( DO_QUASPHER_BEAM .AND. DO_QSREFRAC_BEAM ) THEN
	    MAIL   =
     &  'Bad input: Cannot do Chapman function in refractive atmosphere'
	    ACTION =
     &    'Check DO_CHAPMAN_FUNCTION and DO_QSREFRAC_BEAM flags'
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	  ENDIF
	ENDIF

C  check consistency of mean value input control

	IF ( DO_ADDITIONAL_MVOUT ) THEN
	  IF ( DO_MVOUT_ONLY ) THEN
	    MAIL   = 'Bad input: Cannot have both mean-value flags set'
	    ACTION = 'Check DO_MVOUT_ONLY and DO_ADDITIONAL_MVOUT flags'
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	  ENDIF
	ELSE
	  IF ( DO_MVOUT_ONLY ) THEN
	    IF ( .NOT. DO_NO_AZIMUTH ) THEN
	      MAIL   =
     $          'Bad input: Mean-value option requires NO azimuth'
	      ACTION = 'Turn on the DO_NO_AZIMUTH flag'
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	    ENDIF
	    IF ( DO_USER_STREAMS ) THEN
	      MAIL   =
     &           'Bad input: Mean-value option needs quadratures only'
	      ACTION = 'Turn off DO_USER_STREAMS flag'
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	    ENDIF
	  ENDIF
	ENDIF

C  check consistency of multiple scatter source term output control

	IF ( SAVE_LAYER_MSST ) THEN
	  IF ( .NOT. DO_LBOUND_TAUS ) THEN
	    MAIL =
     & 'Bad input: MSCAT. source term - layer boundary TAU flag not set'
	    ACTION = 'Check DO_LBOUND_TAUS and DO_OUTPUT_MSCATSF flags'
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	  ENDIF
	  IF ( .NOT. DO_USER_STREAMS ) THEN
	    MAIL   =
     &   'Bad input: MSCAT. source term - user_streams flag not set'
	    ACTION = 'Check DO_USER_STREAMS and DO_OUTPUT_MSCATSF flags'
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	  ENDIF
	ENDIF

C  Check consistency of Quadrature output flag
C    - Quadrature output only with Uncorrected Full-radiance calculation.

	IF ( DO_QUAD_OUTPUT ) THEN
	  IF ( .NOT. DO_FULLRAD_MODE ) THEN
	    MAIL =
     &     'Bad input: Quadrature output not with Full radiance mode'
	    ACTION =
     &     'Set both DO_QUAD_OUTPUT and DO_FULLRAD_MODE flags'
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	  ELSE
	    IF ( DO_SSCORRECTION ) THEN
	      MAIL =
     &       'Bad input: Quadrature output not with SS correction'
	      ACTION =
     &       'Set both DO_QUAD_OUTPUT and DO_SSCORRECTION flags'
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	    ENDIF
	  ENDIF
	ENDIF

C  Check consistency of Rayleigh-only and Isotropic-only cases

	IF ( DO_RAYLEIGH_ONLY .AND. DO_ISOTROPIC_ONLY ) THEN
	  MAIL   =
     &   'Bad input: Isotropic_only & Rayleigh-only flags both set'
	  ACTION =
     &     'Check DO_RAYLEIGH_ONLY and DO_ISOTROPIC_ONLY flags'
	  STATUS = LIDORT_SERIOUS
	  CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	ENDIF

C  no Delta-M scaling with Rayleigh only

	IF ( DO_RAYLEIGH_ONLY ) THEN
	  IF ( DO_DELTAM_SCALING ) THEN
	    MAIL   =
     &       'Bad input: No delta-M scaling with Rayleigh-only output'
	    ACTION =
     &        'Check DO_RAYLEIGH_ONLY and DO_DELTAM_SCALING flags'
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	  ENDIF
	ENDIF

C  no Delta-M scaling with Isotropic only

	IF ( DO_ISOTROPIC_ONLY ) THEN
	  IF ( DO_DELTAM_SCALING ) THEN
	    MAIL   =
     &       'Bad input: No delta-M scaling with Isotropic_only output'
	    ACTION =
     &     'Check DO_ISOTROPIC_ONLY and DO_DELTAM_SCALING flags'
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	  ENDIF
	ENDIF

C  check optical depth flags

	IF ( DO_USER_TAUS .AND. DO_LBOUND_TAUS ) THEN
	  MAIL   = 'Bad input: both optical depth flags cannot be set'
	  ACTION = 'Check DO_USER_TAUS and DO_LBOUND_TAUS flags'
	  STATUS = LIDORT_SERIOUS
	  CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	ELSE IF ( .NOT. DO_USER_TAUS .AND. .NOT. DO_LBOUND_TAUS ) THEN
	  MAIL   = 'Bad input: Neither optical depth flags are set'
	  ACTION = 'Check DO_USER_TAUS and DO_LBOUND_TAUS flags'
	  STATUS = LIDORT_SERIOUS
	  CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	ENDIF
	  
C  check directional input

	IF ( .NOT.DO_UPWELLING .AND. .NOT. DO_DNWELLING ) THEN
	  MAIL   = 'Bad input: no directional input is set'
	  ACTION = 'Check DO_UPWELLING and DO_DNWELLING flags'
	  STATUS = LIDORT_SERIOUS
	  CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	ENDIF

C  check number of input Legendre moments
C  ======================================

C  (general scattering case)

	IF ( .NOT.DO_RAYLEIGH_ONLY.AND..NOT.DO_ISOTROPIC_ONLY ) THEN

	  IF ( DO_DELTAM_SCALING ) THEN
	    IF ( NMOMENTS_INPUT.LT.2*NSTREAMS ) THEN
	      MAIL =
     &     'Bad input: < 2.NSTREAMS moments for delta-M scaling'
	      ACTION = 
     &         'Re-set number of input moments when using delta-M'
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	    ENDIF
	  ELSE
	    IF ( NMOMENTS_INPUT.LT.2*NSTREAMS-1 ) THEN
	      MAIL =
     &      'Bad input: < 2.NSTREAMS - 1 moments without delta-M'
	      ACTION = 
     &         'Re-set number of input moments (no delta-M)'
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	    ENDIF
	  ENDIF

	ELSE

C  Checks for Rayleigh only option

	  IF ( DO_RAYLEIGH_ONLY ) THEN
	    IF ( NMOMENTS_INPUT.NE.2 ) THEN
	      MAIL   = 'Bad input: Rayleigh-only phase momemts not = 2'
	      ACTION = 'Set NMOMENTS_INPUT = 2 for Rayleigh-only'
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	    ENDIF
	  ENDIF

C  Checks for Isotropic only option

	  IF ( DO_ISOTROPIC_ONLY ) THEN
	    IF ( NMOMENTS_INPUT.NE.0 ) THEN
	      MAIL   =
     &       'Bad input: No phase function moments for isotropic'
	      ACTION = 'Set NMOMENTS_INPUT = 0 for Isotropic_only'
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	    ENDIF
	  ENDIF

	ENDIF

C  Check azimuth-only conditions
C  =============================

C  Check no-Azimuth flag

	IF ( .NOT.DO_NO_AZIMUTH ) THEN
	  IF ( DO_USER_STREAMS. AND. N_USER_STREAMS.EQ.1 ) THEN
	    IF ( USER_ANGLES_INPUT(1) .EQ. ZERO ) THEN
	       MAIL   =
     &    'Bad input: single zenith-sky output requires no azimuth'
	      ACTION =
     &       'Reset DO_NO_AZIMUTH flag (set to true)'
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	    ENDIF
	  ENDIF
	ENDIF

C  Checks for Isotropic only option

	IF ( DO_ISOTROPIC_ONLY ) THEN
 	  IF ( .NOT.DO_NO_AZIMUTH ) THEN
	    MAIL   =
     &       'Bad input: no azimuth dependence for isotropic_only'
	    ACTION = 'Turn on DO_NO_AZIMUTH for isotropic_only'
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	  ENDIF
	ENDIF

C  check viewing geometry input
C  ============================

C  Check earth radius (Chapman function only)

	IF ( DO_CHAPMAN_FUNCTION ) THEN
	  IF ( EARTH_RADIUS.LT.6320.0D0 .OR.
     &         EARTH_RADIUS.GT.6420.0D0 ) THEN
	    MAIL   = 'Bad input: Earh radius not in range [6320,6420]'
	    ACTION = 'Look at EARTH_RADIUS input'
	     STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	  ENDIF
 	ENDIF
	  
C  Check solar zenith angle input

	IF ( SUN0.LT.ZERO .OR.SUN0.GE.90.0 ) THEN
	  MAIL   = 'Bad input: Sun angle not in range [0,90)'
	  ACTION = 'Look at SUN0 input'
	  STATUS = LIDORT_SERIOUS
	  CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
 	ENDIF

C  Check zenith tolerance input

	IF ( ZENITH_TOLERANCE.LE.ZERO .OR.
     &       ZENITH_TOLERANCE.GT.0.001 ) THEN
	  MAIL   = 'Bad input: Zenith tolerance level out of bounds'
	  ACTION = 'Look at ZENITH_TOLERANCE input'
	  STATUS = LIDORT_SERIOUS
	  CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
 	ENDIF

C  Check relative azimuths

	LOOP = .TRUE.
	I = 0
 	DO WHILE (LOOP .AND. I.LT.N_USER_RELAZMS)
	  I = I + 1
	  IF ( USER_RELAZMS(I) .GT. 360.0   .OR.
     &         USER_RELAZMS(I) .LT. ZERO ) THEN
	    WRITE(C2,'(I2)')I
	    MAIL = 'Bad input: out-of-range azimuth angle, no. '//C2
	    ACTION = 'Look at azimuth angle input'
	    LOOP = .FALSE.
	    STATUS = LIDORT_SERIOUS
            CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	  ENDIF
	ENDDO

C  limits on user-defined options

	IF ( .NOT.DO_USER_STREAMS .AND..NOT.DO_QUAD_OUTPUT ) THEN
	  MAIL   = 'Bad input: No angular stream output allowed'
	  ACTION = 'Check DO_USER_STREAMS and DO_QUAD_OUTPUT flags'
	  STATUS = LIDORT_SERIOUS
	  CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	ENDIF

C  check user-defined stream angles (should always be [0,90])

	IF ( DO_USER_STREAMS ) THEN
	  LOOP = .TRUE.
	  I = 0
 	  DO WHILE (LOOP .AND. I.LT.N_USER_STREAMS)
	    I = I + 1
	    IF ( USER_ANGLES_INPUT(I) .GT. 90.0   .OR.
     &           USER_ANGLES_INPUT(I) .LT. ZERO ) THEN
	      WRITE(C2,'(I2)')I
	      MAIL   = 'Bad input: out-of-range user stream, no. '//C2
	      ACTION = 'Look at user-defined angle input'
	      LOOP = .FALSE.
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	    ENDIF
	  ENDDO
	ENDIF

C  Check height grid input (Chapman function only)

	IF ( DO_CHAPMAN_FUNCTION ) THEN
	  LOOP = .TRUE.
	  I = 0
 	  DO WHILE (LOOP .AND. I.LT.NLAYERS)
	    I = I + 1
	    IF ( HEIGHT_GRID(I-1).LE.HEIGHT_GRID(I) ) THEN
	      WRITE(C2,'(I2)')I
	      MAIL   =
     &'Bad input: Height-grid not monotonically decreasing; Layer '//C2
	      ACTION = 'Look at Height-grid input'
	      LOOP = .FALSE.
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	    ENDIF
	  ENDDO
	ENDIF

C  Check optical depth input
C  =========================

C  check user-defined optical depths (should always be within atmosphere!)

	IF ( DO_USER_TAUS ) THEN
	  LOOP = .TRUE.
	  I = 0
 	  DO WHILE (LOOP .AND. I.LT.N_OUT_USERTAUS)
	    I = I + 1
	    IF ( USER_TAUS_INPUT(I) .GT. TAUGRID_INPUT(NLAYERS))  THEN
	      WRITE(C2,'(I2)')I
	      MAIL   = 'Bad input: user optical depth too big, '//C2
	      ACTION = 'Check optical depth values in range'
	      LOOP = .FALSE.
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	    ENDIF
	  ENDDO
	ELSE IF ( DO_LBOUND_TAUS ) THEN
	  IF ( N_OUT_USERTAUS .GT. NLAYERS + 1 ) THEN
	    MAIL   = 'Bad input: Too many Layer boundary optical depths'
	    ACTION = 'Check N_OUT_USERTAUS = NLAYERS + 1'
	    LOOP = .FALSE.
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	  ENDIF
	  DO I = 1, N_OUT_USERTAUS
	    IF ( LBOUND_TAUS_INPUT(I) .GT. NLAYERS ) THEN
	      WRITE(C2,'(I2)')LBOUND_TAUS_INPUT(I)
	      MAIL   = 'Bad input: Layer boundary index > NLAYERS: '//C2
	      ACTION = 'Check LBOUND_TAUS_INPUT layer boundary indices'
	      LOOP = .FALSE.
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	    ENDIF
	  ENDDO
	ENDIF

C  check Number of offgrid user-defined optical depths

	IF ( DO_USER_TAUS ) THEN
	  NSTART = 0
	  DO UTA = 1, N_OUT_USERTAUS
	    TAU = USER_TAUS_INPUT(UTA)
	    LOOP = .TRUE.
	    N = 0
	    DO WHILE (LOOP.AND.N.LT.NLAYERS)
	      N = N + 1
	      IF ( DABS(TAU-TAUGRID_INPUT(N-1)).LT.SMALLNUM ) THEN
	        NSTART = NSTART + 1
	      ENDIF
	    ENDDO
	  ENDDO
	  N_OFFGRID_USERTAUS = N_OUT_USERTAUS - NSTART
	  IF ( N_OFFGRID_USERTAUS .GT. MAX_OFFGRID_USERTAUS ) THEN
	    WRITE(C2,'(I2)')N_OFFGRID_USERTAUS
	    MAIL   = 'Bad input: too many offgrid optical depths : '//C2
	    ACTION =
     &        'Check number of offgrid Tau does not exceed dimensioning'
	    LOOP = .FALSE.
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	  ENDIF
	ENDIF

C  check repetition of user-defined optical depth input

	IF ( DO_USER_TAUS ) THEN
	  UTA = 0
	  LOOP = .TRUE.
	  DO WHILE ( LOOP .AND. UTA .LT. N_OUT_USERTAUS )
	    UTA = UTA + 1
	    TAU = USER_TAUS_INPUT(UTA)
	    NSTART = 0
	    DO N = 1, N_OUT_USERTAUS
	      IF ( TAU .EQ. USER_TAUS_INPUT(N)) NSTART = NSTART + 1
	    ENDDO
	    IF ( NSTART .NE. 1 ) THEN
	      LOOP = .FALSE.
	      MAIL = 'Bad input: repetition of optical depth input value'
	      ACTION = 'Check user-defined optical depth input'
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	    ENDIF
	  ENDDO
	ENDIF
    
C  check repetition of layer boundary optical depth input

	IF ( DO_LBOUND_TAUS ) THEN
	  UTA = 0
	  LOOP = .TRUE.
	  DO WHILE (LOOP .AND. UTA .LT. N_OUT_USERTAUS )
	    UTA = UTA + 1
	    N1 = LBOUND_TAUS_INPUT(UTA)
	    NSTART = 0
	    DO N = 1, N_OUT_USERTAUS
	      IF ( LBOUND_TAUS_INPUT(N) .EQ. N1 ) NSTART = NSTART + 1
	    ENDDO
	    IF ( NSTART .GT. 1 ) THEN
	      LOOP = .FALSE.
	      MAIL = 'Bad input: repetition of optical depth input value'
	      ACTION = 'Check optical depth layer boundary input'
	      STATUS = LIDORT_SERIOUS
	      CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	    ENDIF
	  ENDDO
	ENDIF

C  check geophysical inputs
C  ========================

C  check single scatter albedos

	DO L = 1, NLAYERS
	  IF ( OMEGA_TOTAL_INPUT(L).GT.ONE-OMEGA_SMALLNUM ) THEN
	    WRITE(C3,'(I3)')L
	    MAIL = 'Bad input: SS-albedo too close to 1, layer '//C3
	    ACTION = 'Check SS-albedo input'
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	  ELSE IF ( OMEGA_TOTAL_INPUT(L).LT.OMEGA_SMALLNUM ) THEN
	    WRITE(C3,'(I3)')L
	    MAIL = 'Bad input: SS-albedo too close to 0, layer '//C3
	    ACTION = 'Check SS-albedo input'
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	  ENDIF
	ENDDO

C  Check first phase function moments

	DO L = 1, NLAYERS
	    IF ( PHASMOMS_TOTAL_INPUT(0,L).NE.ONE ) THEN
	    WRITE(C3,'(I3)')L
	    MAIL = '1st phase moment NOT = 1 for layer '//C3
	    ACTION = 'Check First phase function moment'
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
	  ENDIF
	ENDDO

C  Finish

	END
