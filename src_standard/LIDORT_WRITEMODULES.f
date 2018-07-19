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

	SUBROUTINE LIDORT_WRITEFOURIER
     &     ( FUNIT, FOURIER )

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_RESULTS.VARS'

C  input variables

	INTEGER		 FOURIER, FUNIT
	
C  local variables

	INTEGER		 I, UT, IDIR, WDIR
	INTEGER		 NOF, BSIZE, NUM, NBLOCK, NB, NREM

C  optical depth blocks

	BSIZE  = 5
	NBLOCK = N_OUT_USERTAUS / BSIZE
	NREM = 1
	IF (MOD(N_OUT_USERTAUS,BSIZE).EQ.0) NREM = 0

C  write header

	WRITE(FUNIT,'(a)')' '
	write(FUNIT,'(/a,i3/a/)')
     &         'Intensity output for Fourier component',FOURIER,
     &         '------------------------------------------'

	WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of optical depths = ',N_OUT_USERTAUS
	WRITE(FUNIT,FMT_INTEGER)
     &      'Total number of output angles = ',N_OUT_STREAMS

C  detailed output

	DO IDIR = 1, N_DIRECTIONS

	  WDIR = WHICH_DIRECTIONS(IDIR)

C  direction header

	  IF (WDIR .EQ. UPIDX ) THEN
	    WRITE(FUNIT,'(/A)')
     &     '  --> Upwelling intensities all optical depths and angles'
	  ELSE IF (WDIR .EQ. DNIDX ) THEN
	      WRITE(FUNIT,'(/A)')
     &     '  --> Downwelling intensities all optical depths and angles'
	  ENDIF

	  DO NB = 1, NBLOCK + NREM
	    NOF = ( NB - 1 ) * BSIZE
	    IF ( NB .LE. NBLOCK ) NUM = BSIZE
	    IF ( NB .GT. NBLOCK ) NUM = N_OUT_USERTAUS - NOF

C  optical depth header

 	    WRITE(FUNIT,'(/a12,5(2x,1pe10.4,1x))')
     *         'Opt. depth->',(USER_TAUS_INPUT(NOF+UT),UT=1,NUM)
  	    WRITE(FUNIT,'(a7)')'  Angle'

C  write angle output

	    DO I = 1, N_OUT_STREAMS
	      WRITE(FUNIT,'(F9.5,3x,5(1PE13.5))')OUT_ANGLES(I),
     &              (INTENSITY_F(NOF+UT,I,WDIR),UT=1,NUM)
	    ENDDO
	  ENDDO

	ENDDO

C  Finish

	END

C

	SUBROUTINE LIDORT_WRITERESULTS (RUN)

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_RESULTS.VARS'

C  argument

	INTEGER		 RUN

C  local variables

	INTEGER		 I, UA, LOCAL_NUSERAZMS, IDIR, UT, WDIR, N
	INTEGER		 NOF, BSIZE, NUM, NBLOCK, NB, NREM

C  Beam Attenuation summary

	write(RUN,'(/a/a/)')'Beam Attenuation summary',
     &                     '------------------------'

	WRITE(RUN,'(a,F8.4)')'Solar zenith angle = (degrees)',SUN0

	IF ( DO_CLASSICAL_SOLUTION ) THEN
	  write(RUN,'(a)')
     &       'Classical solution of beam particular integral'
	ELSE
	  write(RUN,'(a)')
     &       'Green function solution of beam particular integral'
	ENDIF
	IF ( DO_QUASPHER_BEAM ) THEN
	  write(RUN,'(a)')
     &       'Quasi-spherical approximation to beam attenuation'
	  write(RUN,'(a)')
     &      ' * Average secant approximation was used (no coefficients)'
	ELSE
	  write(RUN,'(a)')
     &       'Plane parallel approximation to beam attenuation'
	ENDIF

C  Fourier Output summary

	write(RUN,'(/a/a/)')'Fourier Output summary',
     &                      '----------------------'

	IF ( DO_FULLRAD_MODE ) THEN
	  write(RUN,'(a)')
     &       'Full radiance calculation has been performed'
	  IF ( DO_SSCORRECTION ) THEN
	  write(RUN,'(a)')
     &  '  --> With the Nakajima-Tanaka TMS single scatter correction'
	  ELSE
	  write(RUN,'(a)')
     &  '  --> No single scatter correction has been applied'
	  ENDIF
	ELSE
	  write(RUN,'(a)')
     &  'ONLY Multiple-scatter radiance calculation has been performed'
	  IF ( SAVE_LAYER_MSST ) THEN
	    write(RUN,'(a)')
     &      '--> Layer multiple-scatter source terms were output'
    	  ENDIF
	ENDIF

	IF ( .NOT.DO_NO_AZIMUTH ) THEN
	  IF ( DOUBLE_CONV_TEST ) THEN
	    write(RUN,'(/a)')
     &   'Double convergence test was used for Fourier Azimuth Series'
	  ELSE
	    write(RUN,'(/a)')
     &   'Single convergence test was used for Fourier Azimuth Series'
	  ENDIF
	  write(RUN,'(a,F10.7/a,I3)')
     &   ' --> Accuracy level was pre-set at : ',LIDORT_ACCURACY,
     &   ' --> Number of Fourier terms used  : ',FOURIER_SAVED
	ELSE
	    write(RUN,'(/a)')
     &   'Azimuth independent output only (Fourier = 0)'
	ENDIF

C  optical depth blocks

	BSIZE  = 5
	NBLOCK = N_OUT_USERTAUS / BSIZE
	NREM = 1
	IF (MOD(N_OUT_USERTAUS,BSIZE).EQ.0) NREM = 0

C  control point for avoiding intensity output
	   
	IF ( DO_MVOUT_ONLY ) GO TO 400

C  intensity output
C  ----------------

C  local number of azimuths

	IF ( DO_NO_AZIMUTH ) THEN
	  LOCAL_NUSERAZMS = 1
	ELSE
	  LOCAL_NUSERAZMS = N_USER_RELAZMS
	ENDIF

C  overall header

	write(RUN,'(/a/a/)')'Intensity output',
     &                       '----------------'

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

C  direction header

	    IF (WDIR .EQ. UPIDX ) THEN
	      WRITE(RUN,'(/A)')
     &     '  --> Upwelling intensities all optical depths and angles'
	    ELSE IF (WDIR .EQ. DNIDX ) THEN
	      WRITE(RUN,'(/A)')
     &     '  --> Downwelling intensities all optical depths and angles'
	    ENDIF

	    DO NB = 1, NBLOCK + NREM
	      NOF = ( NB - 1 ) * BSIZE
	      IF ( NB .LE. NBLOCK ) NUM = BSIZE
	      IF ( NB .GT. NBLOCK ) NUM = N_OUT_USERTAUS - NOF

C  optical depth header

 	      WRITE(RUN,'(/a12,5(2x,1pe10.4,1x))')
     *         'Opt. depth->',(USER_TAUS_INPUT(NOF+UT),UT=1,NUM)
  	      WRITE(RUN,'(a7)')'  Angle'

C  write angle output

	      DO I = 1, N_OUT_STREAMS
	        WRITE(RUN,'(F9.5,3x,5(1PE13.5))')OUT_ANGLES(I),
     &              (INTENSITY(NOF+UT,I,WDIR,UA),UT=1,NUM)
	      ENDDO
	    ENDDO

C  Loops closed

	  ENDDO
	ENDDO

C  Multiple scatter layer source term output
C  -----------------------------------------

	IF ( .NOT. SAVE_LAYER_MSST ) GO TO 400

C  overall header

	write(RUN,'(/a/a/)')'Multiple scatter source term output',
     &                      '-----------------------------------'

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

	    IF (WDIR .EQ. UPIDX ) THEN
	      WRITE(RUN,'(/A)')
     &  '  --> Upwelling mult-scat layer source terms at user angles'
	      WRITE(RUN,'(/A,10(2x,F10.5,1X))')
     &     'angle->', (USER_ANGLES_INPUT(I),I=1,N_USER_STREAMS)
	      WRITE(RUN,'(A)')' layer'	
	      DO N = N_LAYERSOURCE_UP, NLAYERS
	          WRITE(RUN,'(1X,I3,3X,10(1PE13.5))')
     &	         N,(MSCATSTERM(N,I,WDIR,UA),I=1,N_USER_STREAMS)
	      ENDDO
	      WRITE(RUN,'(/A/)')
     &     '  --> Upwelling BOA source terms at same angles'
	      WRITE(RUN,'(7X,10(1PE13.5))')
     &	         (MSCATBOA_SOURCETERM(I,UA),I=1,N_USER_STREAMS)

	    ELSE IF (WDIR .EQ. DNIDX ) THEN
	      WRITE(RUN,'(/A)')
     &  '  --> Downwelling mult-scat layer source terms at user angles'
	      WRITE(RUN,'(/A,10(2x,F10.5,1X))')
     &     'angle->', (USER_ANGLES_INPUT(I),I=1,N_USER_STREAMS)
	      WRITE(RUN,'(A)')' layer'	
	      DO N = 1, N_LAYERSOURCE_DN
	          WRITE(RUN,'(1X,I3,3X,10(1PE13.5))')
     &	             N,(MSCATSTERM(N,I,WDIR,UA),I=1,N_USER_STREAMS)
	      ENDDO
	    ENDIF

C  Loops closed

	  ENDDO
	ENDDO

C  Mean-value output
C  -----------------

400	CONTINUE
	IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN

	  write(RUN,'(/a/a)')'Mean value output',
     &                         '-----------------'

C  detailed output

	  DO IDIR = 1, N_DIRECTIONS

	    WDIR = WHICH_DIRECTIONS(IDIR)

C  direction header

	    IF (WDIR .EQ. UPIDX ) THEN
	      WRITE(RUN,'(/A/)')
     &'  --> Upwelling mean intensities & fluxes, all optical depths'
	    ELSE IF (WDIR .EQ. DNIDX ) THEN
	      WRITE(RUN,'(/A/)')
     &'  --> Downwelling mean intensities & fluxes, all optical depths'
	    ENDIF

C  optical depth header

	    write(RUN,'(a/)')
     &      'optical depth     mean intensity       flux'

C  optical depth loop

	    DO UT = 1, N_OUT_USERTAUS
	      WRITE(RUN,'(2x,F9.5,3x,2E17.7)')
     &              USER_TAUS_INPUT(UT),
     &             MEAN_INTENSITY(UT,WDIR),FLUX(UT,WDIR)
	    ENDDO

C  end direction loop

	  ENDDO

	ENDIF

C  Finish

	END

C

	SUBROUTINE LIDORT_WRITEINPUT ( IUNIT )

C  Read, check and write to file of all control input

	INCLUDE '../include_s/LIDORT.PARS'

C  include files with input variables to be checked out

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  argument

	INTEGER		IUNIT

C  heading and version number

	WRITE(IUNIT, FMT_SECTION)
     &     ' Control input variables for run of LIDORT'
	WRITE(IUNIT, FMT_CHAR)
     &     ' LIDORT Version number = ',LIDORT_VERSION_NUMBER

C  general control

	WRITE(IUNIT, FMT_HEADING) ' control integers'
	IF ( DO_FULL_QUADRATURE ) THEN
	  WRITE(IUNIT, FMT_CHAR)
     &     ' Single Gaussian quadrature over [-1,1]'
	ELSE
	  WRITE(IUNIT, FMT_INTEGER)
     &   ' Double Gaussian quadrature over [-1,0] & [0,1]'
	ENDIF
	WRITE(IUNIT, FMT_INTEGER)
     &     '  ** Number of half space streams = ', NSTREAMS
	WRITE(IUNIT, FMT_INTEGER)
     &     '  ** Number of atmospheric layers = ', NLAYERS
	WRITE(IUNIT, FMT_INTEGER)
     &   '  ** Number of moments (input) = ',
     &             NMOMENTS_INPUT

C	IF ( DO_THERMAL_EMISSION ) THEN
C	  WRITE(IUNIT, FMT_INTEGER)
C     &  '  ** Number of thermal emission coefficients = ',
C     &       N_THERMAL_COEFFS
C	ENDIF

	WRITE(IUNIT, FMT_HEADING) ' flux/accuracy control'
	WRITE(IUNIT, FMT_REAL)
     &     ' Flux constant = ', FLUX_FACTOR
	IF ( .NOT.DO_NO_AZIMUTH ) THEN
	  WRITE(IUNIT, FMT_REAL)
     &  ' accuracy criterion (Fourier convergence) = ', LIDORT_ACCURACY
	ELSE
	  WRITE(IUNIT, FMT_CHAR)
     &  '  ** No Fourier series -- Azimuth-independent term only'
	ENDIF
	WRITE(IUNIT, FMT_REAL)
     &     ' Solar Zenith angle (degrees) = ', SUN0

C  RTE control

	WRITE(IUNIT, FMT_HEADING) ' RTE solution control'

	IF ( DO_DIRECT_BEAM ) THEN
	  WRITE(IUNIT, FMT_CHAR)
     &      ' Direct beam will be included in solution'
	ELSE
	  WRITE(IUNIT, FMT_CHAR)
     &      ' Direct beam will NOT be included in solution'
	ENDIF

C	IF ( DO_THERMAL_EMISSION ) THEN
C	  WRITE(IUNIT, FMT_CHAR)
C     &      ' Thermal Emission will be included in solution'
C	ELSE
C	  WRITE(IUNIT, FMT_CHAR)
C     &      ' NO Thermal emission in the solution'
C	ENDIF

	IF ( DO_SURFACE_EMISSION ) THEN
	  WRITE(IUNIT, FMT_CHAR)
     &      ' Surface Thermal Emission will be included in solution'
	ELSE
	  WRITE(IUNIT, FMT_CHAR)
     &      ' NO Surface Thermal emission in the solution'
	ENDIF

	IF ( DO_LAMBERTIAN_ALBEDO ) THEN
	  WRITE(IUNIT, FMT_CHAR)
     &      ' Surface will be treated as Lambertian'
	ELSE
	  WRITE(IUNIT, FMT_CHAR)
     &      ' Full Bidirectional surface input required'
	ENDIF

	IF ( DO_ISOTROPIC_ONLY ) THEN
	  WRITE(IUNIT, FMT_CHAR)' Medium is isotropic'
	ENDIF

	IF ( DO_RAYLEIGH_ONLY ) THEN
	  WRITE(IUNIT, FMT_CHAR)
     &      ' Medium has Rayleigh scattering only'
	ENDIF

	IF ( .NOT.DO_RAYLEIGH_ONLY.AND..NOT.DO_ISOTROPIC_ONLY ) THEN
	  WRITE(IUNIT, FMT_CHAR)
     &      ' Medium has general scattering phase function'
	ENDIF

C	IF ( DO_TRANSMITTANCE_ONLY ) THEN
C	  WRITE(IUNIT, FMT_CHAR)
C     &      ' RTE solution is Transmission-only (no scattering)'
C	ELSE
C	  WRITE(IUNIT, FMT_CHAR)
C     &      ' RTE solution is full multiple-scattering'
C	ENDIF

C  Beam particular integral control

	WRITE(IUNIT, FMT_HEADING) 'Beam solution control'

	IF ( DO_CLASSICAL_SOLUTION ) THEN
	  WRITE(IUNIT, FMT_CHAR)
     &  ' Beam solution determined by classical (Chandrasekhar) method'
	ELSE
	  WRITE(IUNIT, FMT_CHAR)
     &     ' Beam solution determined by Greens function method'
	ENDIF

	IF ( DO_QUASPHER_BEAM ) THEN
	  WRITE(IUNIT, FMT_CHAR)
     &      ' Beam solution in Quasi-spherical approximation'
	  WRITE(IUNIT, FMT_CHAR)
     &  '  ** Q-S attenuation will use average secant estimation'
	ELSE
	  WRITE(IUNIT, FMT_CHAR)
     &      ' Beam solution in Plane-parallel approximation'
	ENDIF

C  output control

	WRITE(IUNIT, FMT_HEADING) 'Output control'

	IF ( DO_ADDITIONAL_MVOUT ) THEN
	  WRITE(IUNIT, FMT_CHAR)
     &      ' Output of intensities AND fluxes & mean intensities'
	ELSE
	  IF ( DO_MVOUT_ONLY ) THEN
	    WRITE(IUNIT, FMT_CHAR)
     &      ' Output for fluxes & mean intensities ONLY'
	  ELSE
	    WRITE(IUNIT, FMT_CHAR)
     &      ' Output for intensities ONLY'
	  ENDIF
	ENDIF
	IF ( SAVE_LAYER_MSST ) THEN
	  WRITE(IUNIT, FMT_CHAR)
     &   ' Output of multiple-scatter layer source terms also given'
	ENDIF

	IF ( .NOT.DO_NO_AZIMUTH ) THEN
	  WRITE(IUNIT, FMT_INTEGER)
     &     ' Number of user-defined azimuth angles = ', N_USER_RELAZMS
	ENDIF

	IF ( DO_USER_TAUS ) THEN
	  WRITE(IUNIT, FMT_INTEGER)
     &       ' Total number of optical depth output levels = ',
     &        N_OUT_USERTAUS
	ENDIF

	IF ( DO_USER_STREAMS ) THEN
	  IF ( DO_QUAD_OUTPUT ) THEN
	    WRITE(IUNIT, FMT_CHAR)
     &       ' Stream angle output will include quadratures'
	  WRITE(IUNIT, FMT_INTEGER)
     &       '  ** Total Number of output stream angles = ',
     &                    NSTREAMS + N_USER_STREAMS
	  ELSE
	    WRITE(IUNIT, FMT_CHAR)
     &       ' Stream angle output for user-defined angles only'
	  WRITE(IUNIT, FMT_INTEGER)
     &       '  ** Total Number of output stream angles = ',
     &                    N_USER_STREAMS
	  ENDIF
	ELSE
	  IF ( .NOT. DO_MVOUT_ONLY ) THEN
	    WRITE(IUNIT, FMT_CHAR)
     &      ' Stream output at Quadrature angles only'
	    WRITE(IUNIT, FMT_INTEGER)
     &       '  ** Total Number of output stream angles = ',NSTREAMS
	  ENDIF
	ENDIF

C  Finish

	END

C

	SUBROUTINE LIDORT_WRITESCEN ( SUNIT )

C  write to file of all geophysical LIDORT input

	INCLUDE '../include_s/LIDORT.PARS'

C  include files with input variables to be checked out

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  Module input

	INTEGER		 SUNIT

C  local variables

	INTEGER		 N, L, UT

C  heading and version number

	WRITE(SUNIT,FMT_SECTION)
     &     ' Geophysical scenario variables for run of LIDORT'
	WRITE(SUNIT,FMT_CHAR)
     &     ' LIDORT Version number = ',LIDORT_VERSION_NUMBER

C   basic input
C   -----------

	WRITE(SUNIT,FMT_SECTION)
     &     ' Basic atmospheric input for intensity calculation'

C  Layer optical depths and single scatter albedos

	IF ( DO_DELTAM_SCALING ) THEN
	  WRITE(SUNIT, FMT_HEADING)
     &      '  scaled and unscaled inputs with Delta-M turned ON'
	  WRITE(SUNIT, FMT_HEADING)
     &      'Layer optical depths and single-scatter-albedos '
	  WRITE(SUNIT,'(a,a25/a,a25)')
     &       '      (unscaled) (scaled) ','  (unscaled)   (scaled)  ',
     &       'Layer  Op-depth  Op-depth ','  s.s albedo  s.s albedo '
	  WRITE(SUNIT,'(a)')' '
	  DO N = 1, NLAYERS
	    WRITE(SUNIT,'(I3,T6,2(f9.5,1x),2x,2(f9.5,3x))')
     &            N,TAUGRID_INPUT(N),TAUGRID(N),
     &            OMEGA_TOTAL_INPUT(N),OMEGA_TOTAL(N)
	  ENDDO
	ELSE
	  WRITE(SUNIT, FMT_HEADING)
     &      '  Unscaled inputs only with Delta-M turned OFF'
	  WRITE(SUNIT, FMT_HEADING)
     &    'Layer optical depths and single-scatter-albedos '
	  WRITE(SUNIT,'(a,a13/a,a13)')
     &        '      (unscaled) ','  (unscaled) ',
     &        'Layer  Op-depth  ','  s.s albedo '
	  WRITE(SUNIT,'(a)')' '
	  DO N = 1, NLAYERS
	    WRITE(SUNIT,'(I3,T6,f9.5,4x,f9.5)')
     &            N,TAUGRID_INPUT(N),OMEGA_TOTAL_INPUT(N)
	  ENDDO
	ENDIF

C  First 5 and final phase function moments

	IF ( .NOT.DO_RAYLEIGH_ONLY .AND. .NOT.DO_ISOTROPIC_ONLY ) THEN
	  IF ( DO_DELTAM_SCALING ) THEN
	    WRITE(SUNIT, FMT_HEADING)
     &         '1-5 and final SCALED TOTAL phase fn. moments'
	    WRITE(SUNIT,'(a,T53,a6,T60,a6,I3/)')
     &       'Layer | Moments 0 through 4-->',' ... ','Moment',NMOMENTS
	    DO N = 1, NLAYERS
	      WRITE(SUNIT,'(I3,5x,5f9.5,a5,f9.5)')
     &             N,(PHASMOMS_TOTAL(L,N),L=0,4),' ... ',
     &             PHASMOMS_TOTAL(NMOMENTS,N)
	    ENDDO
	  ELSE
	    WRITE(SUNIT, FMT_HEADING)
     &         '1-5 and final unscaled TOTAL phase fn. moments'
	    WRITE(SUNIT,'(a,T53,a6,T60,a6,I3/)')
     &       'Layer | Moments 0 through 4-->',' ... ','Moment',NMOMENTS
	    DO N = 1, NLAYERS
	      WRITE(SUNIT,'(I3,5x,5f9.5,a5,f9.5)')
     &             N,(PHASMOMS_TOTAL_INPUT(L,N),L=0,4),' ... ',
     &             PHASMOMS_TOTAL_INPUT(NMOMENTS,N)
	    ENDDO
	  ENDIF
	ENDIF

C  Commented out thermal expansion coefficients
C	IF ( DO_THERMAL_EMISSION ) THEN
C	  WRITE(SUNIT,FMT_HEADING)'thermal emission coefficients'
C	  WRITE(SUNIT,'(a/)')
C     &       'Layer | thermal emission expansion coeffs-->'
C	  DO N = 1, NLAYERS
C	    WRITE(SUNIT,'(I3,4x,10(f10.5))')
C     &          N,(THERMAL_COEFFS(N,S),S=1,N_THERMAL_COEFFS)
C	  ENDDO
C	ENDIF

C  surface emission

	WRITE(SUNIT,FMT_HEADING)'Surface reflecting property'
	IF ( DO_LAMBERTIAN_ALBEDO ) THEN

	  WRITE(SUNIT,FMT_REAL)'(Lambertian) surface albedo is',ALBEDO
	  WRITE(SUNIT,FMT_REAL)'(Lambertian) emissivity = ',ONE-ALBEDO

	ELSE

C  P L A C E H O L D E R

	ENDIF

	IF ( DO_SURFACE_EMISSION ) THEN
	  WRITE(SUNIT,FMT_REAL)'Surface blackbody function is', SURFBB
	ENDIF

C   Geometry input
C   --------------

	WRITE(SUNIT,FMT_SECTION)' Viewing geometry input'

	WRITE(SUNIT,FMT_HEADING)
     &         'Computational (quadrature) angles in the half space'
	WRITE(SUNIT,'(a/)')
     &           'Stream  |  Angle    |   Cosine   |   Weight'
	DO N = 1, NSTREAMS
	  WRITE(SUNIT,'(I3,5x,3(1x,f9.5,3x))')N,XANG(N),X(N),A(N)
	ENDDO

	WRITE(SUNIT,FMT_HEADING) 'User-defined stream angles'
	WRITE(SUNIT,'(a/)')'Stream  |  Angle    |   Cosine'
	DO N = 1, N_USER_STREAMS
	  WRITE(SUNIT,'(I3,5x,2(1x,f9.5,3x))')
     *              N,USER_ANGLES_INPUT(N),USER_STREAMS(N)
	ENDDO

	IF ( DO_NO_AZIMUTH ) THEN
	  WRITE(SUNIT,FMT_HEADING) 'No azimuth angles'
	ELSE
	  WRITE(SUNIT,FMT_HEADING)
     &          'User-defined relative azimuth angles'
	  WRITE(SUNIT,'(a/)')'Number |   Angle'
	  DO N = 1, N_USER_RELAZMS
	    WRITE(SUNIT,'(I3,5x,1x,f9.5)') N,USER_RELAZMS(N)
	  ENDDO
	ENDIF

C  Optical depth input
C  -------------------

	WRITE(SUNIT,FMT_SECTION)
     &       ' User-defined optical depths for output'
	WRITE(SUNIT,'(a/)')'Optical depth | Level/Layer of occurrence'
	UT = 0
	DO N = 1, N_OUT_USERTAUS
	  IF ( OFFGRID_UTAU_OUTFLAG(N) ) THEN
	    UT = UT + 1
	    WRITE(SUNIT,'(1x,f9.5,2x,a9,i3)')
     &         USER_TAUS_INPUT(N),'    layer',OFFGRID_UTAU_LAYERIDX(UT)
	  ELSE
	    WRITE(SUNIT,'(1x,f9.5,2x,a9,i3)')
     &         USER_TAUS_INPUT(N),'    level',UTAU_LEVEL_MASK_UP(N)
	  ENDIF
	ENDDO

C  Finish

	END

C  Commented out IDL results write
C
	SUBROUTINE LIDORT_IDLWRITERESULTS (RUN)

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_RESULTS.VARS'
	INTEGER		 RUN
	INTEGER		 I, UA, LOCAL_NUSERAZMS, IDIR, UT, WDIR
	INTEGER		 IOUT, J
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
	      WRITE(RUN,'(F9.5,f6.2,5(1PE13.5))')
     &             OUT_ANGLES(I),ANGLEMASK(I),
     &              (INTENSITY(UT,I,WDIR,UA),UT=1,N_OUT_USERTAUS)
	    ENDDO
	  ENDDO
	ENDDO

	END

