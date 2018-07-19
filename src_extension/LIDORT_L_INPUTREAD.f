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

	SUBROUTINE LIDORT_L_INPUTREAD ( STATUS )

C  Read all linearization control inputs for LIDORT
C  ------------------------------------------------

C  include file of constants

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  include files to be filled out with input data

	INCLUDE '../include_e/LIDORT_L_CONTROL.VARS'

C  Module arguments (error filename,output status)

	INTEGER		STATUS

C  local variables

	CHARACTER*8	PREFIX
	PARAMETER	( PREFIX = 'LIDORT -' )
	LOGICAL		ERROR
        CHARACTER * 80  PAR_STR
        LOGICAL         GFINDPAR
        EXTERNAL        GFINDPAR
	INTEGER		FILUNIT, LEN_STRING, I
	CHARACTER*70	MAIL
	EXTERNAL	LEN_STRING

C  initialize status

	STATUS = LIDORT_SUCCESS
	ERROR = .FALSE.

C  file unit

	FILUNIT = LIDORT_INUNIT

C  read all variables in include file LIDORT_L_CONTROL.VARS

C  linearization control

        PAR_STR = 'Simulation only?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_SIMULATION_ONLY
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Variation for atmospheric layers?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_LAYER_LINEARIZATION
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Variation for surface albedo?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_ALBEDO_LINEARIZATION
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

        PAR_STR = 'Variation for surface emission?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) 
     &     READ (FILUNIT,*,ERR=998) DO_SURFBB_LINEARIZATION
	CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

C  layer weighting function names

	IF ( DO_LAYER_LINEARIZATION ) THEN

          PAR_STR = 'Total number of Layer weighting functions'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR))
     &    READ (FILUNIT,*,ERR=998) TOTAL_LAYERWF
	  CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

          PAR_STR = 'Layer weighting function names'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
	    DO I = 1, TOTAL_LAYERWF
	      READ (FILUNIT,'(a)',ERR=998) LAYERWF_NAMES(I)
	    ENDDO
	  ENDIF
	  CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS )

	ENDIF

C  set overall linearization flag

	DO_LINEARIZATION = ( DO_LAYER_LINEARIZATION  .OR.
     &                       DO_ALBEDO_LINEARIZATION .OR.
     &                       DO_SURFBB_LINEARIZATION )

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

	SUBROUTINE LIDORT_INIT_L_CONTROL_VARS

C  Initialises all control inputs for LIDORT
C  -----------------------------------------

C  include file of constants

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  include files to be filled out with input data

	INCLUDE '../include_e/LIDORT_L_CONTROL.VARS'

	DO_SIMULATION_ONLY        = .FALSE.
	DO_LINEARIZATION          = .FALSE.
	DO_LAYER_LINEARIZATION    = .FALSE.
	DO_ALBEDO_LINEARIZATION   = .FALSE.
	DO_SURFBB_LINEARIZATION   = .FALSE.

	END
