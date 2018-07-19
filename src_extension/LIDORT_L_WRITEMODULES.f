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

	SUBROUTINE LIDORT_L_WRITEFOURIER
     &     ( DO_INCLUDE_ALBEDO, FUNIT, FOURIER )

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_e/LIDORT_L_CONTROL.VARS'
	INCLUDE '../include_e/LIDORT_L_GEOPHYS.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

C  input variables

	INTEGER		 FOURIER, FUNIT
	LOGICAL		 DO_INCLUDE_ALBEDO
	
C  local variables

	INTEGER		 I, IDIR, UT, WDIR
	INTEGER		 NOF, BSIZE, NUM, NBLOCK, NB, NREM, Q, K

C  optical depth blocks

	BSIZE  = 5
	NBLOCK = N_OUT_USERTAUS / BSIZE
	NREM = 1
	IF (MOD(N_OUT_USERTAUS,BSIZE).EQ.0) NREM = 0

C  Atmospheric Weighting function output
C  #####################################

	IF ( DO_LAYER_LINEARIZATION ) THEN

C  overall header

	  WRITE(FUNIT,'(a)')' '
	  write(FUNIT,'(/a,i3/a/)')
     & 'Atmospheric Weighting function output for Fourier component',
     &         FOURIER,
     & '--------------------------------------------------------------'

	  WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of optical depths = ',N_OUT_USERTAUS
	  WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of output angles = ',N_OUT_STREAMS

C  detailed output

	  DO IDIR = 1, N_DIRECTIONS
	    WDIR = WHICH_DIRECTIONS(IDIR)

C  linearization control

	    DO K = 1, NLAYERS
	      IF ( LAYER_VARY_FLAG(K) ) THEN
	        DO Q = 1, LAYER_VARY_NUMBER(K)

C  direction headers

	          IF (WDIR .EQ. UPIDX ) THEN
	            WRITE(FUNIT,'(/A/A,I2/A,A31)')
     &     '  --> Upwelling Atmospheric Weighting functions',
     &     '        * for variations in layer = ',K,
     &     '        * with respect to : ',LAYERWF_NAMES(Q)
	          ELSE IF (WDIR .EQ. DNIDX ) THEN
	            WRITE(FUNIT,'(/A/A,I2/A,A31)')
     &     '  --> Downwelling Atmospheric Weighting functions',
     &     '        * for variations in layer = ',K,
     &     '        * with respect to : ',LAYERWF_NAMES(Q)
	          ENDIF

C  output block

	          DO NB = 1, NBLOCK + NREM
	            NOF = ( NB - 1 ) * BSIZE
	            IF ( NB .LE. NBLOCK ) NUM = BSIZE
	            IF ( NB .GT. NBLOCK ) NUM = N_OUT_USERTAUS-NOF

C  optical depth header

 	            WRITE(FUNIT,'(/a12,5(2x,1pe10.4,1x))')
     *         'Opt. depth->',(USER_TAUS_INPUT(NOF+UT),UT=1,NUM)
  	            WRITE(FUNIT,'(a7)')'  Angle'

C  write angle output

	            DO I = 1, N_OUT_STREAMS
	              WRITE(FUNIT,'(F9.5,3x,5(1PE13.5))')
     &               OUT_ANGLES(I),
     &                (ATMOSWF_F(Q,K,NOF+UT,I,WDIR),UT=1,NUM)
	            ENDDO
	          ENDDO

C  close control and loops

	        ENDDO
	      ENDIF
	    ENDDO
	  ENDDO

C  End atmospheric weighting function stuff

	ENDIF

C  Albedo Weighting function output
C  ################################

	IF ( DO_ALBEDO_LINEARIZATION.AND.DO_INCLUDE_ALBEDO ) THEN

C  only required for non-Lambertian cases

	  IF ( DO_LAMBERTIAN_ALBEDO ) RETURN

C  overall header

	  WRITE(FUNIT,'(a)')' '
	  write(FUNIT,'(/a,i3/a/)')
     & 'Albedo Weighting function output for Fourier component',
     &         FOURIER,
     & '----------------------------------------------------------'

	  WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of optical depths = ',N_OUT_USERTAUS
	  WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of output angles = ',N_OUT_STREAMS

C  detailed output

	  DO IDIR = 1, N_DIRECTIONS
	    WDIR = WHICH_DIRECTIONS(IDIR)

C  direction headers

	    IF (WDIR .EQ. UPIDX ) THEN
	      WRITE(FUNIT,'(/A)')
     &     '  --> Upwelling Albedo Weighting functions'
	    ELSE IF (WDIR .EQ. DNIDX ) THEN
	      WRITE(FUNIT,'(/A)')
     &     '  --> Downwelling Albedo Weighting functions'
	    ENDIF

C  output block

	    DO NB = 1, NBLOCK + NREM
	      NOF = ( NB - 1 ) * BSIZE
	      IF ( NB .LE. NBLOCK ) NUM = BSIZE
	      IF ( NB .GT. NBLOCK ) NUM = N_OUT_USERTAUS-NOF

C  optical depth header

 	      WRITE(FUNIT,'(/a12,5(2x,1pe10.4,1x))')
     *         'Opt. depth->',(USER_TAUS_INPUT(NOF+UT),UT=1,NUM)
  	      WRITE(FUNIT,'(a7)')'  Angle'

C  write angle output

	      DO I = 1, N_OUT_STREAMS
	        WRITE(FUNIT,'(F9.5,3x,5(1PE13.5))')
     &         OUT_ANGLES(I),(ALBEDOWF_F(NOF+UT,I,WDIR),UT=1,NUM)
	      ENDDO
	    ENDDO

	  ENDDO

C  End albedo weighting function stuff

	ENDIF

C  Finish

	END

C

	SUBROUTINE LIDORT_L_WRITERESULTS (RUN)

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_e/LIDORT_L_CONTROL.VARS'
	INCLUDE '../include_e/LIDORT_L_GEOPHYS.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

C  argument

	INTEGER		 RUN
	
C  local variables

	INTEGER		 I, UA, LOCAL_NUSERAZMS, IDIR, UT, WDIR
	INTEGER		 NOF, BSIZE, NUM, NBLOCK, NB, NREM, Q, K, N

C  optical depth blocks

	BSIZE  = 5
	NBLOCK = N_OUT_USERTAUS / BSIZE
	NREM = 1
	IF (MOD(N_OUT_USERTAUS,BSIZE).EQ.0) NREM = 0

C  Atmospheric Weighting function output
C  #####################################

	IF ( DO_LAYER_LINEARIZATION ) THEN

C  control point for avoiding weighting function output
	   
	  IF ( DO_MVOUT_ONLY ) GO TO 401

C  local number of azimuths

	  IF ( DO_NO_AZIMUTH ) THEN
	    LOCAL_NUSERAZMS = 1
	  ELSE
	    LOCAL_NUSERAZMS = N_USER_RELAZMS
	  ENDIF

C  overall header

	  write(RUN,'(/a/a/)')'Atmospheric Weighting function output',
     &                        '-------------------------------------'

C  start azimuth loop

	  DO UA = 1, LOCAL_NUSERAZMS

C  azimuth angle header

	    IF ( DO_NO_AZIMUTH ) THEN
	      write(RUN,FMT_CHAR)
     *       '** Results FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
	    ELSE
	      write(RUN,FMT_REAL)
     * '** RESULTS FOR RELATIVE AZIMUTH ANGLE (degs)=',USER_RELAZMS(UA)
	    ENDIF
	    WRITE(RUN,'(a)')' '
	    WRITE(RUN,FMT_INTEGER)
     &      'Total number of optical depths = ',N_OUT_USERTAUS
	    WRITE(RUN,FMT_INTEGER)
     &      'Total number of output angles = ',N_OUT_STREAMS

C  detailed output

	    DO IDIR = 1, N_DIRECTIONS
	      WDIR = WHICH_DIRECTIONS(IDIR)

C  linearization control

	      DO K = 1, NLAYERS
	        IF ( LAYER_VARY_FLAG(K) ) THEN
	          DO Q = 1, LAYER_VARY_NUMBER(K)

C  direction headers

	            IF (WDIR .EQ. UPIDX ) THEN
	              WRITE(RUN,'(/A/A,I2/A,A31)')
     &     '  --> Upwelling Atmospheric Weighting functions',
     &     '        * for variations in layer = ',K,
     &     '        * with respect to : ',LAYERWF_NAMES(Q)
	            ELSE IF (WDIR .EQ. DNIDX ) THEN
	              WRITE(RUN,'(/A/A,I2/A,A31)')
     &     '  --> Downwelling Atmospheric Weighting functions',
     &     '        * for variations in layer = ',K,
     &     '        * with respect to : ',LAYERWF_NAMES(Q)
	            ENDIF

C  output block

	            DO NB = 1, NBLOCK + NREM
	              NOF = ( NB - 1 ) * BSIZE
	              IF ( NB .LE. NBLOCK ) NUM = BSIZE
	              IF ( NB .GT. NBLOCK ) NUM = N_OUT_USERTAUS-NOF

C  optical depth header

 	              WRITE(RUN,'(/a12,5(2x,1pe10.4,1x))')
     *         'Opt. depth->',(USER_TAUS_INPUT(NOF+UT),UT=1,NUM)
  	              WRITE(RUN,'(a7)')'  Angle'

C  write angle output

	              DO I = 1, N_OUT_STREAMS
	                WRITE(RUN,'(F9.5,3x,5(1PE13.5))')
     &                 OUT_ANGLES(I),
     &                  (ATMOSWF(Q,K,NOF+UT,I,WDIR,UA),UT=1,NUM)
	              ENDDO
	            ENDDO

C  close control and loops

	          ENDDO
	        ENDIF
	      ENDDO
	    ENDDO
	  ENDDO

C  Mult-scat source term output
C  ----------------------------

	  IF ( .NOT. SAVE_LAYER_MSST ) GO TO 401

C  overall header

	  write(RUN,'(/a/a/)')
     &  'Mult-scat source term atmospheric weighting function output',
     &  '-----------------------------------------------------------'

C  start azimuth loop

	  DO UA = 1, LOCAL_NUSERAZMS

C  azimuth angle header

	    IF ( DO_NO_AZIMUTH ) THEN
	      write(RUN,FMT_CHAR)
     *       '** Results FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
	    ELSE
	      write(RUN,FMT_REAL)
     * '** RESULTS FOR RELATIVE AZIMUTH ANGLE (degs)=',USER_RELAZMS(UA)
	    ENDIF

C  detailed output

	    DO IDIR = 1, N_DIRECTIONS

	      WDIR = WHICH_DIRECTIONS(IDIR)

C  Upwelling
C  ~~~~~~~~~

	      IF (WDIR .EQ. UPIDX ) THEN

C  linearization control

	        DO K = 1, NLAYERS
	          IF ( LAYER_VARY_FLAG(K) ) THEN
	            DO Q = 1, LAYER_VARY_NUMBER(K)

C  direction headers

	              WRITE(RUN,'(/A/A,I2/A,I2)')
     &   '  --> Upwelling Mult-scat Source Term Weighting functions',
     &     '        * W.R.T parameters varying in layer = ',K,
     &     '        * Parameter number = ',Q

C  output

	              WRITE(RUN,'(/A,10(2x,F10.5,1X))')
     &     'angle->', (USER_ANGLES_INPUT(I),I=1,N_USER_STREAMS)
	              WRITE(RUN,'(A)')' layer'	
	              DO N = N_LAYERSOURCE_UP, NLAYERS
	                WRITE(RUN,'(1X,I3,3X,10(1PE13.5))')
     &	     N,(ATMOSWF_MSCATSTERM(Q,K,N,I,WDIR,UA),I=1,N_USER_STREAMS)
	              ENDDO
	              WRITE(RUN,'(/A/)')
     &   '  --> Upwelling diffuse BOA source term weighting functions'
	              WRITE(RUN,'(7X,10(1PE13.5))')
     &	         (ATMOSWF_MSCATBOA_STERM(Q,K,I,UA),I=1,N_USER_STREAMS)

	            ENDDO
	          ENDIF
	        ENDDO

C  Downwelling
C  ~~~~~~~~~~~

	      ELSE IF (WDIR .EQ. DNIDX ) THEN

	        DO K = 1, NLAYERS
	          IF ( LAYER_VARY_FLAG(K) ) THEN
	            DO Q = 1, LAYER_VARY_NUMBER(K)

	              WRITE(RUN,'(/A/A,I2/A,I2)')
     &   '  --> Downwelling Mult-scat Source Term Weighting functions',
     &     '        * W.R.T parameters varying in layer = ',K,
     &     '        * Parameter number = ',Q
	              WRITE(RUN,'(/A,10(2x,F10.5,1X))')
     &     'angle->', (USER_ANGLES_INPUT(I),I=1,N_USER_STREAMS)
	              WRITE(RUN,'(A)')' layer'	
	              DO N = 1, N_LAYERSOURCE_DN
	                WRITE(RUN,'(1X,I3,3X,10(1PE13.5))')
     &	     N,(ATMOSWF_MSCATSTERM(Q,K,N,I,WDIR,UA),I=1,N_USER_STREAMS)
	              ENDDO

	            ENDDO
	          ENDIF
	        ENDDO
	      ENDIF

C  Loops closed

	    ENDDO
	  ENDDO

C  Mean-value output
C  -----------------

401	  CONTINUE
	  IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN

	    write(RUN,'(/a/a)')
     &          'Mean value atmospheric weighting function output',
     &          '------------------------------------------------'

C  detailed output

	    DO IDIR = 1, N_DIRECTIONS
	      WDIR = WHICH_DIRECTIONS(IDIR)

C  linearization control

	      DO K = 1, NLAYERS
	        IF ( LAYER_VARY_FLAG(K) ) THEN
	          DO Q = 1, LAYER_VARY_NUMBER(K)

C  direction headers

	            IF (WDIR .EQ. UPIDX ) THEN
	              WRITE(RUN,'(/A/A,I2/A,I2/)')
     & '  --> Upwelling Mean intensity and Flux Weighting functions',
     &     '        * W.R.T parameters varying in layer = ',K,
     &     '        * Parameter number = ',Q
	            ELSE IF (WDIR .EQ. DNIDX ) THEN
	              WRITE(RUN,'(/A/A,I2/A,I2/)')
     & '  --> Downwelling Mean intensity and Flux Weighting functions',
     &     '        * W.R.T parameters varying in layer = ',K,
     &     '        * Parameter number = ',Q
	            ENDIF

C  optical depth header

	            WRITE(RUN,'(a/)')
     &      'optical depth     mean int. WF      flux WF'

C  optical depth loop

	            DO UT = 1, N_OUT_USERTAUS
	              WRITE(RUN,'(2x,F9.5,3x,2E17.7)')
     &                    USER_TAUS_INPUT(UT),
     &                    MINT_ATMOSWF(Q,K,UT,WDIR),
     &                    FLUX_ATMOSWF(Q,K,UT,WDIR)
	            ENDDO

C  close control and loops

	          ENDDO
	        ENDIF
	      ENDDO
	    ENDDO
	  ENDIF

C  End atmospheric weighting function stuff

	ENDIF

C  Albedo Weighting function output
C  ################################

	IF ( DO_ALBEDO_LINEARIZATION ) THEN

C  control point for avoiding weighting function output
	   
	  IF ( DO_MVOUT_ONLY ) GO TO 402

C  local number of azimuths

	  IF ( DO_NO_AZIMUTH ) THEN
	    LOCAL_NUSERAZMS = 1
	  ELSE
	    IF ( DO_LAMBERTIAN_ALBEDO ) THEN
	      LOCAL_NUSERAZMS = 1
	    ELSE
	      LOCAL_NUSERAZMS = N_USER_RELAZMS
	    ENDIF
	  ENDIF

C  overall header

	  write(RUN,'(/a/a/)')'Albedo Weighting function output',
     &                        '--------------------------------'

C  start azimuth loop

	  DO UA = 1, LOCAL_NUSERAZMS

C  azimuth angle header

	    IF ( DO_NO_AZIMUTH ) THEN
	      write(RUN,FMT_CHAR)
     *       '** Results FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
	    ELSE
	      IF ( DO_LAMBERTIAN_ALBEDO ) THEN
	        write(RUN,FMT_CHAR)
     *       '** Lambertian Albedo --> Results AZIMUTH-INDEPENDENT **'
	      ELSE
	        write(RUN,FMT_REAL)
     * '** RESULTS FOR RELATIVE AZIMUTH ANGLE (degs)=',USER_RELAZMS(UA)
	      ENDIF
	    ENDIF
	    WRITE(RUN,'(a)')' '
	    WRITE(RUN,FMT_INTEGER)
     &      'Total number of optical depths = ',N_OUT_USERTAUS
	    WRITE(RUN,FMT_INTEGER)
     &      'Total number of output angles = ',N_OUT_STREAMS

C  detailed output

	    DO IDIR = 1, N_DIRECTIONS
	      WDIR = WHICH_DIRECTIONS(IDIR)

C  direction headers

	      IF (WDIR .EQ. UPIDX ) THEN
	        WRITE(RUN,'(/A)')
     &     '  --> Upwelling Albedo Weighting functions'
	      ELSE IF (WDIR .EQ. DNIDX ) THEN
	        WRITE(RUN,'(/A)')
     &     '  --> Downwelling Albedo Weighting functions'
	      ENDIF

C  output block

	      DO NB = 1, NBLOCK + NREM
	        NOF = ( NB - 1 ) * BSIZE
	        IF ( NB .LE. NBLOCK ) NUM = BSIZE
	        IF ( NB .GT. NBLOCK ) NUM = N_OUT_USERTAUS-NOF

C  optical depth header

 	        WRITE(RUN,'(/a12,5(2x,1pe10.4,1x))')
     *         'Opt. depth->',(USER_TAUS_INPUT(NOF+UT),UT=1,NUM)
  	        WRITE(RUN,'(a7)')'  Angle'

C  write angle output

	        DO I = 1, N_OUT_STREAMS
	          WRITE(RUN,'(F9.5,3x,5(1PE13.5))')
     &         OUT_ANGLES(I),(ALBEDOWF(NOF+UT,I,WDIR,UA),UT=1,NUM)
	        ENDDO
	      ENDDO

C  close control and loops

	    ENDDO
	  ENDDO

C  Mult-scat source term output
C  ----------------------------

	  IF ( .NOT. SAVE_LAYER_MSST ) GO TO 402

C  overall header

	  write(RUN,'(/a/a/)')
     &       'Mult-scat source term albedo weighting function output',
     &       '------------------------------------------------------'

C  start azimuth loop

	  DO UA = 1, LOCAL_NUSERAZMS

C  azimuth angle header

	    IF ( DO_NO_AZIMUTH ) THEN
	      write(RUN,FMT_CHAR)
     *       '** Results FOR AZIMUTH-INDEPENDENT COMPONENT ONLY **'
	    ELSE
	      write(RUN,FMT_REAL)
     * '** RESULTS FOR RELATIVE AZIMUTH ANGLE (degs)=',USER_RELAZMS(UA)
	    ENDIF

C  detailed output

	    DO IDIR = 1, N_DIRECTIONS

	      WDIR = WHICH_DIRECTIONS(IDIR)

C  Upwelling
C  ~~~~~~~~~

	      IF (WDIR .EQ. UPIDX ) THEN

C  direction header and output

	        WRITE(RUN,'(/A)')
     &     '  --> Upwelling Mult-scat Source Term Weighting functions'
	        WRITE(RUN,'(/A,10(2x,F10.5,1X))')
     &     'angle->', (USER_ANGLES_INPUT(I),I=1,N_USER_STREAMS)
	        WRITE(RUN,'(A)')' layer'	
	        DO N = N_LAYERSOURCE_UP, NLAYERS
	          WRITE(RUN,'(1X,I3,3X,10(1PE13.5))')
     &	     N,(ALBEDOWF_MSCATSTERM(N,I,WDIR,UA),I=1,N_USER_STREAMS)
	        ENDDO
	        WRITE(RUN,'(/A/)')
     & '  --> Upwelling Mult-scat BOA Source Term Weighting function'
	        WRITE(RUN,'(7X,10(1PE13.5))')
     &	     (ALBEDOWF_MSCATBOA_STERM(I,UA),I=1,N_USER_STREAMS)


C  Downwelling
C  ~~~~~~~~~~~

	      ELSE IF (WDIR .EQ. DNIDX ) THEN

	        WRITE(RUN,'(/A)')
     &     '  --> Downwelling Mult-scat Source Term Weighting functions'
	        WRITE(RUN,'(/A,10(2x,F10.5,1X))')
     &     'angle->', (USER_ANGLES_INPUT(I),I=1,N_USER_STREAMS)
	        WRITE(RUN,'(A)')' layer'	
	        DO N = 1, N_LAYERSOURCE_DN
	          WRITE(RUN,'(1X,I3,3X,10(1PE13.5))')
     &	     N,(ALBEDOWF_MSCATSTERM(N,I,WDIR,UA),I=1,N_USER_STREAMS)
	        ENDDO

	      ENDIF
	    ENDDO
	  ENDDO

C  Mean-value output
C  -----------------

402	  CONTINUE
	  IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN

	  write(RUN,'(/a/a)')
     &          'Mean value albedo weighting function output',
     &          '-------------------------------------------'

C  detailed output

	    DO IDIR = 1, N_DIRECTIONS
	      WDIR = WHICH_DIRECTIONS(IDIR)

C  direction headers

	      IF (WDIR .EQ. UPIDX ) THEN
	        WRITE(RUN,'(/A/A/)')
     & '  --> Upwelling Mean intensity and Flux Weighting functions',
     &     '      * W.R.T albedo variation '
	      ELSE IF (WDIR .EQ. DNIDX ) THEN
	        WRITE(RUN,'(/A/A/)')
     & '  --> Downwelling Mean intensity and Flux Weighting functions',
     &     '      * W.R.T albedo variation '
	      ENDIF

C  optical depth header

	      WRITE(RUN,'(a/)')
     &      'optical depth     mean int. WF      flux WF'

C  optical depth loop

	      DO UT = 1, N_OUT_USERTAUS
	        WRITE(RUN,'(2x,F9.5,3x,2E17.7)')
     &            USER_TAUS_INPUT(UT),
     &            MINT_ALBEDOWF(UT,WDIR),
     &            FLUX_ALBEDOWF(UT,WDIR)
	      ENDDO

C  close control and loops

	    ENDDO
	  ENDIF

C  End albedo weighting function stuff

	ENDIF

C  Surfbb Weighting function output
C  ################################

	IF ( DO_SURFBB_LINEARIZATION ) THEN

C  control point for avoiding weighting function output
	   
	  IF ( DO_MVOUT_ONLY ) GO TO 403

C  overall header

	  write(RUN,'(/a/a/)')'Surfbb Weighting function output',
     &                        '--------------------------------'

C  headers

	  WRITE(RUN,FMT_CHAR)'** Results are AZIMUTH-INDEPENDENT **'
	  WRITE(RUN,'(a)')' '
	  WRITE(RUN,FMT_INTEGER)
     &      'Total number of optical depths = ',N_OUT_USERTAUS
	  WRITE(RUN,FMT_INTEGER)
     &      'Total number of output angles = ',N_OUT_STREAMS

C  detailed output

	  DO IDIR = 1, N_DIRECTIONS
	    WDIR = WHICH_DIRECTIONS(IDIR)

C  direction headers

	    IF (WDIR .EQ. UPIDX ) THEN
	      WRITE(RUN,'(/A)')
     &     '  --> Upwelling Surfbb Weighting functions'
	    ELSE IF (WDIR .EQ. DNIDX ) THEN
	      WRITE(RUN,'(/A)')
     &     '  --> Downwelling Surfbb Weighting functions'
	    ENDIF

C  output block

	    DO NB = 1, NBLOCK + 1
	      NOF = ( NB - 1 ) * BSIZE
	      IF ( NB .LE. NBLOCK ) NUM = BSIZE
	       IF ( NB .GT. NBLOCK ) NUM = N_OUT_USERTAUS-NOF

C  optical depth header

 	      WRITE(RUN,'(/a12,5(2x,1pe10.4,1x))')
     *         'Opt. depth->',(USER_TAUS_INPUT(NOF+UT),UT=1,NUM)
  	      WRITE(RUN,'(a7)')'  Angle'

C  write angle output

	      DO I = 1, N_OUT_STREAMS
	        WRITE(RUN,'(F9.5,3x,5(1PE13.5))')
     &         OUT_ANGLES(I),(SURFBBWF(NOF+UT,I,WDIR),UT=1,NUM)
	      ENDDO
	    ENDDO

	  ENDDO

C  Mult-scat source term output
C  ----------------------------

	  IF ( .NOT. SAVE_LAYER_MSST ) GO TO 403

C  overall header

	  write(RUN,'(/a/a/)')
     &      'Mult-scat source term Surfbb weighting function output',
     &      '------------------------------------------------------'

	  WRITE(RUN,FMT_CHAR)'** Results are AZIMUTH-INDEPENDENT **'

C  detailed output

	  DO IDIR = 1, N_DIRECTIONS
	    WDIR = WHICH_DIRECTIONS(IDIR)

C  Upwelling
C  ~~~~~~~~~

	    IF (WDIR .EQ. UPIDX ) THEN

	      WRITE(RUN,'(/A)')
     &     '  --> Upwelling Mult-scat Source Term Weighting functions'
	      WRITE(RUN,'(/A,10(2x,F10.5,1X))')
     &     'angle->', (USER_ANGLES_INPUT(I),I=1,N_USER_STREAMS)
	      WRITE(RUN,'(A)')' layer'	
	      DO N = N_LAYERSOURCE_UP, NLAYERS
	        WRITE(RUN,'(1X,I3,3X,10(1PE13.5))')
     &	     N,(SURFBBWF_MSCATSTERM(N,I,WDIR),I=1,N_USER_STREAMS)
	      ENDDO
	      WRITE(RUN,'(/A/)')
     &       '  --> Upwelling BOA Source Term Weighting function'
	      WRITE(RUN,'(7X,10(1PE13.5))')
     &	       (SURFBBWF_BOA_STERM(I),I=1,N_USER_STREAMS)

C  Downwelling
C  ~~~~~~~~~~~

	    ELSE IF (WDIR .EQ. DNIDX ) THEN

	      WRITE(RUN,'(/A)')
     &  '  --> Downwelling Mult-scat Source Term Weighting functions'
	      WRITE(RUN,'(/A,10(2x,F10.5,1X))')
     &     'angle->', (USER_ANGLES_INPUT(I),I=1,N_USER_STREAMS)
	      WRITE(RUN,'(A)')' layer'	
	      DO N = 1, N_LAYERSOURCE_DN
	        WRITE(RUN,'(1X,I3,3X,10(1PE13.5))')
     &	     N,(SURFBBWF_MSCATSTERM(N,I,WDIR),I=1,N_USER_STREAMS)
	      ENDDO

	    ENDIF
	  ENDDO

C  Mean-value output
C  -----------------

403	  CONTINUE

	  IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN

	  write(RUN,'(/a/a)')
     &          'Mean value surfbb weighting function output',
     &          '-------------------------------------------'

C  detailed output

	    DO IDIR = 1, N_DIRECTIONS
	      WDIR = WHICH_DIRECTIONS(IDIR)

C  direction headers

	      IF (WDIR .EQ. UPIDX ) THEN
	        WRITE(RUN,'(/A/A/)')
     & '  --> Upwelling Mean intensity and Flux Weighting functions',
     &     '      * W.R.T surfbb variation '
	      ELSE IF (WDIR .EQ. DNIDX ) THEN
	        WRITE(RUN,'(/A/A/)')
     & '  --> Downwelling Mean intensity and Flux Weighting functions',
     &     '      * W.R.T surfbb variation '
	      ENDIF

C  optical depth header

	      WRITE(RUN,'(a/)')
     &      'optical depth     mean int. WF      flux WF'

C  optical depth loop

	      DO UT = 1, N_OUT_USERTAUS
	        WRITE(RUN,'(2x,F9.5,3x,2E17.7)')
     &            USER_TAUS_INPUT(UT),
     &            MINT_SURFBBWF(UT,WDIR),
     &            FLUX_SURFBBWF(UT,WDIR)
	      ENDDO

C  close control and loops

	    ENDDO
	  ENDIF

C  End albedo weighting function stuff

	ENDIF

C  Finish

	END

C

	SUBROUTINE LIDORT_L_WRITEINPUT ( IUNIT )

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  include files with input variables to be checked out

	INCLUDE '../include_e/LIDORT_L_CONTROL.VARS'

C  argument

	INTEGER		IUNIT
	INTEGER		I

C  Linearization control

	WRITE(IUNIT, FMT_HEADING) 'Linearization control'

	IF ( DO_SIMULATION_ONLY ) THEN
	  WRITE(IUNIT, FMT_CHAR)
     &      ' Output will be intensity only (no weight. fns.)'
	ELSE
	  IF ( DO_LAYER_LINEARIZATION ) THEN
	    WRITE(IUNIT, FMT_CHAR)
     &    ' Atmospheric layer weighting functions will be output w.r.t:'
	    DO I = 1, TOTAL_LAYERWF
	      WRITE(IUNIT, FMT_CHAR)LAYERWF_NAMES(I)
	    ENDDO
	  ENDIF
	  WRITE(IUNIT,'(A)')' '
	  IF ( DO_ALBEDO_LINEARIZATION ) THEN
	    WRITE(IUNIT, FMT_CHAR)
     &      ' Albedo weighting functions will be output'
	  ENDIF
	  IF ( DO_SURFBB_LINEARIZATION ) THEN
	    WRITE(IUNIT, FMT_CHAR)
     &      ' Surface blackbody weighting functions will be output'
	  ENDIF
	ENDIF

C  Finish

	END

C

	SUBROUTINE LIDORT_L_WRITESCEN ( SUNIT )

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  include files with input variables to be checked out

	INCLUDE '../include_e/LIDORT_L_CONTROL.VARS'
	INCLUDE '../include_e/LIDORT_L_GEOPHYS.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'

C  Module input

	INTEGER		 SUNIT

C  local variables

	INTEGER		 N, V

C  set up

C   WF input
C   --------

C  layer variational quantities

	IF ( DO_LAYER_LINEARIZATION ) THEN

	  WRITE(SUNIT,FMT_SECTION)
     &     'Atmospheric input for Weighting function calculation'

	  WRITE(SUNIT, FMT_HEADING)
     &        'Single scattering albedo variations'
	  WRITE(SUNIT,'(a/)')
     &     'Layer varying | parameter numbers-->'
	  DO N = 1, NLAYERS
	    IF ( LAYER_VARY_FLAG(N) ) THEN
	       WRITE(SUNIT,'(I3,14x,5(1pe12.4))')
     &        N,(OMEGA_VARS_TOTAL_INPUT(V,N),V=1,LAYER_VARY_NUMBER(N))
	    ENDIF
	  ENDDO

	  WRITE(SUNIT, FMT_HEADING) 'Optical depth variations'
	  WRITE(SUNIT,'(a/)')
     &     'Layer varying | parameter numbers-->'
	  DO N = 1, NLAYERS
	    IF ( LAYER_VARY_FLAG(N) ) THEN
	      WRITE(SUNIT,'(I3,14x,5(1pe12.4))')
     &         N,(EXT_VARS_INPUT(V,N),V=1,LAYER_VARY_NUMBER(N))
	    ENDIF
	  ENDDO

	ENDIF

C  Finish

	END


C  Commented out IDL results write
C
	SUBROUTINE LIDORT_L_IDLWRITERESULTS (RUN)

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_e/LIDORT_L_CONTROL.VARS'
	INCLUDE '../include_e/LIDORT_L_GEOPHYS.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

	INTEGER		 RUN
	INTEGER		 I, UA, LOCAL_NUSERAZMS, IDIR, UT, WDIR
	INTEGER		 IOUT, J, Q, K
	DOUBLE PRECISION ANGLEMASK(MAX_OUT_STREAMS)

	LOCAL_NUSERAZMS = N_USER_RELAZMS
	IF ( DO_NO_AZIMUTH ) LOCAL_NUSERAZMS = 1

C  angle mask

	IF ( DO_QUAD_OUTPUT ) THEN
	  DO I = 1, NSTREAMS
	    IOUT = QUADOUTPUT_INDEX(I)
	    DO J = 1, N_OUT_STREAMS
	      IF (J.EQ.IOUT) ANGLEMASK(J) = ONE
	    ENDDO
	  ENDDO
	ENDIF
	IF ( DO_USER_STREAMS ) THEN
	  DO I = 1, N_USER_STREAMS
	    IOUT = USEROUTPUT_INDEX(I)
	    DO J = 1, N_OUT_STREAMS
	      IF (J.EQ.IOUT) ANGLEMASK(J) = ZERO
	    ENDDO
	  ENDDO
	ENDIF

C  IDL OUTPUT

	DO UA = 1, LOCAL_NUSERAZMS
C  azimuth angle header
	  IF ( DO_NO_AZIMUTH ) THEN
	    write(RUN,'(f9.3)')999.999
	  ELSE
	    write(RUN,'(f9.3)')USER_RELAZMS(UA)
	  ENDIF
C  idl output (debug only)
	  DO IDIR = 1, N_DIRECTIONS
	    WDIR = WHICH_DIRECTIONS(IDIR)
	    if (idir.eq.1)write(RUN,'(2i4)')N_OUT_USERTAUS,N_OUT_STREAMS
 	    WRITE(RUN,'(100(1x,1pe11.4,1x))')
     *       (USER_TAUS_INPUT(UT), UT = 1, N_OUT_USERTAUS)
	    DO I = 1, N_OUT_STREAMS
              WRITE(RUN,'(2f10.4)')OUT_ANGLES(I),ANGLEMASK(I)
	      DO K = 1, NLAYERS
	        IF ( LAYER_VARY_FLAG(K) ) THEN
	          DO Q = 1, LAYER_VARY_NUMBER(K)
	            WRITE(RUN,'(5(1PE13.5))')
     &               (ATMOSWF(Q,K,UT,I,WDIR,UA), UT = 1, N_OUT_USERTAUS)
	          ENDDO
	        ENDIF
	      ENDDO
	    ENDDO
	  ENDDO
	ENDDO

	END


