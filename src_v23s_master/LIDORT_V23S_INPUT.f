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

	SUBROUTINE LIDORT_V23S_INPUT ( FILNAM, EFILNAM, STATUS )

C  Read all control inputs for LIDORT
C  ----------------------------------

C  include file of constants

	INCLUDE '../include_s/LIDORT.PARS'

C  include file of CONTROL variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'

C  Module arguments (input filename, error filename,output status)

	CHARACTER*(*)	FILNAM
	CHARACTER*(*)	EFILNAM
	INTEGER		STATUS

C  local variables

	INTEGER		STATUS_SUB, FILUNIT, LEN_STRING
	CHARACTER*70	MAIL
	EXTERNAL	LEN_STRING

C  initialize status

	STATUS = LIDORT_SUCCESS

C  initialize error file

	LIDORT_ERROR_FILENAME = EFILNAM
	LIDORT_ERROR_INIT     = .TRUE.

C  initialize variables

	CALL LIDORT_INIT_CONTROL_VARS
	CALL LIDORT_INIT_MODEL_VARS

C  Open file

	FILUNIT = LIDORT_INUNIT
	OPEN(LIDORT_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

C  read standard inputs

	CALL LIDORT_READINPUT ( STATUS_SUB )
	IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
	  STATUS = LIDORT_SERIOUS
	  MAIL = 'Error from LIDORT_INPUTREAD'
	  CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
	  CLOSE(FILUNIT)
	  RETURN
	ENDIF

C  normal return

	CLOSE(FILUNIT)
	RETURN
	
C  Open file error

300 	CONTINUE
	STATUS = LIDORT_SERIOUS
	MAIL = 'openfile failure for '//FILNAM(1:LEN_STRING(FILNAM))
	CALL LIDORT_ERROR_TRACE ( MAIL, ' ', STATUS )
	CLOSE(FILUNIT)
	RETURN

C  Finish

	END

