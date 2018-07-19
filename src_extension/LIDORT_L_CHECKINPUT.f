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

	SUBROUTINE LIDORT_L_CHECK_INPUT (STATUS)

C  check inputs (both file-read and derived)

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  include files with input variables to be checked out

	INCLUDE '../include_e/LIDORT_L_CONTROL.VARS'

C  Module output

	INTEGER		 STATUS

C  local variables

	CHARACTER*70	MAIL
	CHARACTER*1	NOTRACE
	PARAMETER	( NOTRACE = ' ')

C  Initialize output status

	STATUS = LIDORT_SUCCESS

C  Do some checking
C  ================

C  Check formatting; should not exceed 10 entries for the number of
C  layer weighting functions. This is needed to avoid errors with
C  compilers not using variable formatting.

C	IF ( MAX_PARAMETERS .GT. 10 ) THEN
C	  MAIL='Max layer variations > 10, get formatting problems!'
C	  STATUS = LIDORT_SERIOUS
C	  CALL LIDORT_ERROR_TRACE ( MAIL, NOTRACE, STATUS )
C	ENDIF

C  Check something is being varied

	IF ( .NOT. DO_SIMULATION_ONLY ) THEN
	  IF ( .NOT. DO_LAYER_LINEARIZATION .AND.
     &         .NOT. DO_ALBEDO_LINEARIZATION. AND.
     &         .NOT. DO_SURFBB_LINEARIZATION  ) THEN
	    MAIL='Bad input: Nothing being varied or calculated'
	    STATUS = LIDORT_SERIOUS
	    CALL LIDORT_ERROR_TRACE ( MAIL, NOTRACE, STATUS )
	  ENDIF
	ENDIF

C  Finish

	END
