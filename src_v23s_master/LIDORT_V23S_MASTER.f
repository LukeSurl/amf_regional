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

	SUBROUTINE LIDORT_V23S_MASTER
     &     ( STATUS_INPUTCHECK, STATUS_CALCULATION )

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  subroutine arguments
C  --------------------

	INTEGER		 STATUS_INPUTCHECK
	INTEGER		 STATUS_CALCULATION

C  local variables
C  ---------------

	LOGICAL		 LOCAL_DO_NO_AZIMUTH,SAVE_DO_NO_AZIMUTH
	LOGICAL		 LOCAL_ITERATION
	LOGICAL		 DO_INCLUDE_ALBEDO,DO_INCLUDE_SURFEMISS

	INTEGER		 FOURIER, N_FOURIER, I, I1, UA, TESTCONV
	INTEGER		 LOCAL_N_USERAZM, STATUS_SUB
	INTEGER		 IUNIT, SUNIT, FUNIT, RUNIT

	DOUBLE PRECISION AZMFAC(MAX_USER_RELAZMS) 
	DOUBLE PRECISION X2(MAXSTRM2),A2(MAXSTRM2), SSFLUX

C  Local error handling

	CHARACTER*70	 MAIL, TRACE
	CHARACTER*3	 CF

C  initialize output status

	STATUS_CALCULATION = LIDORT_SUCCESS
	STATUS_INPUTCHECK  = LIDORT_SUCCESS

C  Check input
C  -----------

	CALL LIDORT_CHECK_INPUT ( STATUS_SUB )

	IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
	  STATUS_INPUTCHECK = LIDORT_SERIOUS
	  MAIL = ' LIDORT_CHECKINPUT failed'
	  TRACE = ' Called in LIDORT_V23S_MASTER '
	  CALL LIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS_INPUTCHECK )
	  CLOSE(LIDORT_ERRUNIT)
	  RETURN
	ENDIF

C  Chapman function calculation of TAUTHICK_INPUT
C  ----------------------------------------------

	IF ( DO_CHAPMAN_FUNCTION ) THEN
	  CALL CHAPMAN_FUNCTION
	ENDIF

C  write input variables
C  ---------------------

	IF ( DO_WRITE_INPUT ) THEN

C  open file

	  IUNIT = LIDORT_INUNIT
	  OPEN(IUNIT,FILE=INPUT_WRITE_FILENAME,STATUS='UNKNOWN')

C  standard input

	  CALL LIDORT_WRITEINPUT ( IUNIT )

C  close file

	  CLOSE(IUNIT)

	ENDIF

C  Set Quadrature abscissae and weights
C  ====================================

	IF ( DO_FULL_QUADRATURE ) THEN
	  CALL GAULEG(-ONE,ONE,X2,A2,NSTR2)
	  DO I = 1, NSTREAMS
	    I1 = I + NSTREAMS
	    X(I) = X2(I1)
	    A(I) = A2(I1)
	  ENDDO
	ELSE
	  CALL GAULEG(ZERO,ONE,X,A,NSTREAMS)
	ENDIF

C  Get derived inputs
C  ==================

	CALL LIDORT_DERIVE_INPUT

C  Initialise Fourier loop
C  =======================

C  Set Number of Fourier terms (NMOMENTS = Maximum).
C    ( Starting from 0 = Fundamental )

	SAVE_DO_NO_AZIMUTH  = DO_NO_AZIMUTH
	LOCAL_DO_NO_AZIMUTH = DO_NO_AZIMUTH

C  No azimuth dependency for following three cases

C	IF ( DO_TRANSMITTANCE_ONLY ) THEN
C	  LOCAL_DO_NO_AZIMUTH = .TRUE.
C	ENDIF

	IF ( DO_ISOTROPIC_ONLY  ) THEN
	  LOCAL_DO_NO_AZIMUTH = .TRUE.
	ENDIF

	IF ( DO_MVOUT_ONLY  ) THEN
	  LOCAL_DO_NO_AZIMUTH = .TRUE.
	ENDIF

C  set Fourier number (2 for Rayleigh only)

	IF ( LOCAL_DO_NO_AZIMUTH ) THEN
	  N_FOURIER = 0
	ELSE
	  IF ( DO_RAYLEIGH_ONLY  ) THEN
	    N_FOURIER = 2
	  ELSE
	    N_FOURIER = NMOMENTS
	  ENDIF
	ENDIF

C  re-set no-azimuth flag

	DO_NO_AZIMUTH = LOCAL_DO_NO_AZIMUTH

C  if there's no azimuth dependence, just do one value in azimuth loop

	IF ( DO_NO_AZIMUTH ) THEN
	  LOCAL_N_USERAZM = 1
 	ELSE
	  LOCAL_N_USERAZM = N_USER_RELAZMS
	ENDIF

C  Fourier loop
C  ============

	LOCAL_ITERATION = .TRUE.
	FOURIER = -1
	TESTCONV = 0

	DO WHILE ( LOCAL_ITERATION .AND. FOURIER.LT.N_FOURIER )

C  Fourier counter

	  FOURIER = FOURIER + 1

C  Local start of user-defined streams

	  IF ( FOURIER. EQ. 0 ) THEN
	    LOCAL_UM_START = 1
	  ELSE
	    LOCAL_UM_START = N_OUT_STREAMS - N_CONV_STREAMS + 1
	  ENDIF

C  azimuth cosine factor

	  IF ( FOURIER .GT. 0 ) THEN
	    DO UA = 1, LOCAL_N_USERAZM
	      AZMFAC(UA) = DCOS(DEG_TO_RAD*USER_RELAZMS(UA)*FOURIER)
	    ENDDO
	  ENDIF

C  Main call to Lidort Fourier module

C	  write(*,*)' ..calculating fourier component',FOURIER
	  CALL LIDORT_V23S_FOURIER
     I       (FOURIER,DO_INCLUDE_ALBEDO,DO_INCLUDE_SURFEMISS,
     O        STATUS_SUB)

C  error handling

	  IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
	    STATUS_CALCULATION = LIDORT_SERIOUS
	    WRITE(CF,'(I3)')FOURIER
	    MAIL = 'Error from LIDORT_V23S_FOURIER, Fourier = '//CF
	    TRACE = 'Called by LIDORT_V23S_MASTER'
	    CALL LIDORT_ERROR_TRACE
     &           ( MAIL, TRACE, STATUS_CALCULATION )
	    CLOSE(LIDORT_ERRUNIT)
	    RETURN
	  ENDIF

C  Convergence examination

	  CALL LIDORT_CONVERGE
     I      ( AZMFAC, FOURIER, LOCAL_N_USERAZM,
     O        TESTCONV, LOCAL_ITERATION )

C  Fourier output
C  --------------

	  IF ( DO_WRITE_FOURIER ) THEN

C  Open file if Fourier = 0

	    FUNIT = LIDORT_FUNIT
	    IF ( FOURIER .EQ. 0 ) THEN
              OPEN(FUNIT,FILE=FOURIER_WRITE_FILENAME,STATUS='UNKNOWN')
	    ENDIF

C  Write Standard Fourier output

	    CALL LIDORT_WRITEFOURIER ( FUNIT, FOURIER )

C  Close file if iteration has finished

	    IF ( .NOT.LOCAL_ITERATION ) CLOSE ( FUNIT )

	  ENDIF

C  end iteration loop

	ENDDO

C  restore no azimuth flag

	DO_NO_AZIMUTH = SAVE_DO_NO_AZIMUTH

C  Single scatter correction
C  =========================

	IF ( DO_SSCORRECTION ) THEN
	  SSFLUX  = QUARTER*FLUX_FACTOR/PIE
	  CALL LIDORT_SSCORRECTION ( SSFLUX ) 
	ENDIF

C  Major result output
C  ===================

C  Standard output
C  - - - - - - - -

	IF ( DO_WRITE_RESULTS ) THEN

C  Open file

	  RUNIT = LIDORT_RESUNIT
          OPEN(RUNIT,FILE=RESULTS_WRITE_FILENAME,STATUS='UNKNOWN')

C  Write standard results

	  CALL LIDORT_WRITERESULTS ( RUNIT )

C  close file

	  CLOSE(RUNIT)

	ENDIF

C  Write results (IDL output) (commented out)
C  - - - - - - - - - - - - - 

C	IF ( DO_WRITE_RESULTS_IDL ) THEN
C	  GUNIT = LIDORT_IDLUNIT
C          OPEN(GUNIT,FILE=IDLRESULTS_WRITE_FILENAME,STATUS='UNKNOWN')
C	  CALL LIDORT_IDLWRITERESULTS ( GUNIT )
C	  CLOSE(GUNIT)
C	ENDIF

C  Geophysical input (scenario) write
C  ----------------------------------

	IF ( DO_WRITE_SCENARIO ) THEN
	  SUNIT = LIDORT_SCENUNIT
	  OPEN(SUNIT,FILE=SCENARIO_WRITE_FILENAME,STATUS='UNKNOWN')
	  CALL LIDORT_WRITESCEN ( SUNIT )
	  CLOSE(SUNIT)
	ENDIF

C  close Error file if it was used

	IF ( .NOT. LIDORT_ERROR_INIT ) CLOSE(LIDORT_ERRUNIT)

C  Finish

	END
