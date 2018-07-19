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

	SUBROUTINE LIDORT_CONVERGE
     &      ( AZMFAC, FOURIER, LOCAL_N_USERAZM,
     &        TESTCONV, LOCAL_ITERATION )

C  convergence test on the intensity

C  Include files
C  -------------

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  Output results

	INCLUDE '../include_s/LIDORT_RESULTS.VARS'

C  Subroutine arguments
C  --------------------

C  input variables

	INTEGER		 FOURIER, LOCAL_N_USERAZM
	DOUBLE PRECISION AZMFAC(MAX_USER_RELAZMS)

C  modified/output variables

	LOGICAL		 LOCAL_ITERATION
	INTEGER		 TESTCONV

C  local variables
C  ---------------

	INTEGER		 COUNT, COUNT_A, I, IDIR, UT, UA, WDIR, N
	DOUBLE PRECISION TNEW, TOLD, TAZM

C  ###################
C  Fourier 0 component
C  ###################

	IF ( FOURIER.EQ.0 ) THEN

C  Copy Fourier component at all output angles and optical depths

	  DO IDIR = 1, N_DIRECTIONS
	    WDIR = WHICH_DIRECTIONS(IDIR)
	    DO UT = 1, N_OUT_USERTAUS
	      DO I = LOCAL_UM_START, N_OUT_STREAMS
	        DO UA = 1, LOCAL_N_USERAZM
	          INTENSITY(UT,I,WDIR,UA) = INTENSITY_F(UT,I,WDIR)
	        ENDDO
	      ENDDO
	    ENDDO
	  ENDDO

C  Multiple scatter source terms (Copy Fourier component)

	  IF ( SAVE_LAYER_MSST ) THEN
	    IF ( DO_DNWELLING ) THEN
	      DO N = 1, N_LAYERSOURCE_DN
	        DO I = LOCAL_UM_START, N_USER_STREAMS
	          DO UA = 1, LOCAL_N_USERAZM
	            MSCATSTERM(N,I,DNIDX,UA) = MSCATSTERM_F(N,I,DNIDX)
	          ENDDO
	        ENDDO
	      ENDDO
	    ENDIF
	    IF ( DO_UPWELLING ) THEN
	      DO I = LOCAL_UM_START, N_USER_STREAMS
	        DO UA = 1, LOCAL_N_USERAZM
	          MSCATBOA_SOURCETERM(I,UA) = MSCATBOA_SOURCETERM_F(I)
	          DIRECTBOA_SOURCETERM(I,UA) = DIRECTBOA_SOURCETERM_F(I)
	        ENDDO
	      ENDDO
	      DO N = N_LAYERSOURCE_UP, NLAYERS
	        DO I = LOCAL_UM_START, N_USER_STREAMS
	          DO UA = 1, LOCAL_N_USERAZM
	            MSCATSTERM(N,I,UPIDX,UA) = MSCATSTERM_F(N,I,UPIDX)
	          ENDDO
	        ENDDO
	      ENDDO
	    ENDIF
	  ENDIF

C  If no_azimuth, then set output and exit flag

	  IF ( DO_NO_AZIMUTH ) THEN
	    LOCAL_ITERATION = .FALSE.
	    RETURN
	  ENDIF

C  ######################
C  Fourier components > 0
C  ######################

	ELSE

C  No examination of convergence
C  -----------------------------

C  For Rayleigh atmosphere or if All Fourier components are required,
C     skip convergence test on intensity

	  IF ( DO_RAYLEIGH_ONLY. OR. DO_ALL_FOURIER ) THEN

C  For each azimuth, add Fourier component

	    DO UA = 1, LOCAL_N_USERAZM

C     - for direction, user optical depth, out stream

	      DO IDIR = 1, N_DIRECTIONS
	        WDIR = WHICH_DIRECTIONS(IDIR)
	        DO UT = 1, N_OUT_USERTAUS
	          DO I = LOCAL_UM_START, N_OUT_STREAMS
	            TOLD = INTENSITY(UT,I,WDIR,UA)
	            TAZM = AZMFAC(UA)*INTENSITY_F(UT,I,WDIR)
	            INTENSITY(UT,I,WDIR,UA) = TOLD + TAZM
	          ENDDO
	        ENDDO
	      ENDDO

	    ENDDO

C  Examine convergence on intensity only 
C  -------------------------------------

C  convergence test applied to ALL directions AND
C                              ALL stream values (except near zenith) AND
C                              ALL azimuths taken together
C                              ALL user optical depths

	  ELSE

C  Count number of occasions Fourier term addition is below accuracy level

	    COUNT = 0
	    DO UA = 1, N_USER_RELAZMS
	      COUNT_A = 0
	      DO IDIR = 1, N_DIRECTIONS
	        WDIR = WHICH_DIRECTIONS(IDIR)
	        DO UT = 1, N_OUT_USERTAUS
	          DO I = LOCAL_UM_START, N_OUT_STREAMS
	            TOLD = INTENSITY(UT,I,WDIR,UA)
	            TAZM = AZMFAC(UA)*INTENSITY_F(UT,I,WDIR)
	            TNEW = TOLD + TAZM
	            IF ( TAZM .NE. ZERO ) THEN
	              IF ( DABS(TAZM/TNEW) .LT. LIDORT_ACCURACY ) THEN
	                COUNT   = COUNT + 1
	                COUNT_A = COUNT_A + 1
	              ENDIF
	            ELSE
	              COUNT   = COUNT + 1
	              COUNT_A = COUNT_A + 1
	            ENDIF
	            INTENSITY(UT,I,WDIR,UA) = TNEW
	          ENDDO
	        ENDDO
	      ENDDO
	    ENDDO

C  set convergence counter TESTCONV (convergence requires twice)

	    IF ( COUNT .EQ. N_CONVTESTS ) THEN
	      TESTCONV = TESTCONV + 1
	      IF ( DOUBLE_CONV_TEST ) THEN
	        IF ( TESTCONV .EQ. 2 ) THEN
                  LOCAL_ITERATION = .FALSE.
	        ENDIF
	      ELSE
                LOCAL_ITERATION = .FALSE.
	      ENDIF
	      IF ( .NOT. LOCAL_ITERATION ) THEN
	        FOURIER_SAVED = FOURIER
	      ENDIF
	    ELSE
	      TESTCONV = 0
	    ENDIF

C  end convergence clause

	  ENDIF

C  For Rayleigh scattering alone, stop iteration after third harmonic

	  IF ( DO_RAYLEIGH_ONLY ) THEN
	    IF ( FOURIER .EQ. 2 ) THEN
	      LOCAL_ITERATION = .FALSE.
	      FOURIER_SAVED = FOURIER
	    ENDIF
	  ENDIF

C  For all Fourier, keep saveing the output number of Fourier terms

	  IF ( DO_ALL_FOURIER ) THEN
	    FOURIER_SAVED = FOURIER
	  ENDIF

C  Multiple scatter source terms (if flagged)

	  IF ( SAVE_LAYER_MSST ) THEN
	    DO UA = 1, LOCAL_N_USERAZM
	      IF ( DO_DNWELLING ) THEN
	        DO N = 1, N_LAYERSOURCE_DN
	          DO I = LOCAL_UM_START, N_USER_STREAMS
	            TOLD = MSCATSTERM(N,I,DNIDX,UA)
	            TAZM = AZMFAC(UA)*MSCATSTERM_F(N,I,DNIDX)
	            MSCATSTERM(N,I,DNIDX,UA) = TOLD + TAZM
	          ENDDO
	        ENDDO
	      ENDIF
	      IF ( DO_UPWELLING ) THEN
	        DO I = LOCAL_UM_START, N_USER_STREAMS
	          TOLD = MSCATBOA_SOURCETERM(I,UA)
	          TAZM = AZMFAC(UA)*MSCATBOA_SOURCETERM_F(I)
	          MSCATBOA_SOURCETERM(I,UA) = TOLD + TAZM
	          TOLD = DIRECTBOA_SOURCETERM(I,UA)
	          TAZM = AZMFAC(UA)*DIRECTBOA_SOURCETERM_F(I)
	          DIRECTBOA_SOURCETERM(I,UA) = TOLD + TAZM
	        ENDDO
	        DO N = N_LAYERSOURCE_UP, NLAYERS
	          DO I = LOCAL_UM_START, N_USER_STREAMS
	            TOLD = MSCATSTERM(N,I,UPIDX,UA)
	            TAZM = AZMFAC(UA)*MSCATSTERM_F(N,I,UPIDX)
	            MSCATSTERM(N,I,UPIDX,UA) = TOLD + TAZM
	          ENDDO
	        ENDDO
	      ENDIF
	    ENDDO
	  ENDIF

C  Finish iteration loop
   
	ENDIF

C  Finish

	END
