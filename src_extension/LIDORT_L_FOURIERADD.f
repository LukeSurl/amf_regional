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

	SUBROUTINE LIDORT_L_FOURIERADD
     &      ( AZMFAC, FOURIER, LOCAL_N_USERAZM )

C  Just upgrades the weighting function Fourier cosine series

C  include files of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'
	INCLUDE '../include_e/LIDORT_L.PARS'

C  Include files of input variables

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

	INCLUDE '../include_e/LIDORT_L_CONTROL.VARS'
	INCLUDE '../include_e/LIDORT_L_GEOPHYS.VARS'

C  Output results

	INCLUDE '../include_e/LIDORT_L_RESULTS.VARS'

C  input variables

	INTEGER		 FOURIER, LOCAL_N_USERAZM	
	DOUBLE PRECISION AZMFAC(MAX_USER_RELAZMS)

C  local variables

	INTEGER		 I, IDIR, UT, UA, Q, N, K, WDIR

C  ###################
C  Fourier 0 component
C  ###################

	IF ( FOURIER.EQ.0 ) THEN

C  atmospheric weighting functions
C  -------------------------------

	  IF ( DO_LAYER_LINEARIZATION ) THEN

C  Full output at all output angles

	    DO N = 1, NLAYERS
	      IF ( LAYER_VARY_FLAG(N) ) THEN
	        DO Q = 1, LAYER_VARY_NUMBER(N)
	          DO IDIR = 1, N_DIRECTIONS
	            WDIR = WHICH_DIRECTIONS(IDIR)
	            DO UT = 1, N_OUT_USERTAUS
	              DO I = 1, N_OUT_STREAMS
	                DO UA = 1, LOCAL_N_USERAZM
	                  ATMOSWF(Q,N,UT,I,WDIR,UA) =
     $                    ATMOSWF_F(Q,N,UT,I,WDIR)
	                ENDDO
	              ENDDO
	            ENDDO
	          ENDDO
	        ENDDO
	      ENDIF
	    ENDDO

C  M-SCAT source term weighting functions (copy Fourier component)

	    IF ( SAVE_LAYER_MSST ) THEN
	      IF ( DO_DNWELLING ) THEN
	        DO N = 1, NLAYERS
	          IF ( LAYER_VARY_FLAG(N) ) THEN
	            DO Q = 1, LAYER_VARY_NUMBER(N)
	              DO K = 1, N_LAYERSOURCE_DN
	                DO I = 1, N_USER_STREAMS
	                  DO UA = 1, LOCAL_N_USERAZM
	                    ATMOSWF_MSCATSTERM(Q,N,K,I,DNIDX,UA) =
     $                      ATMOSWF_MSCATSTERM_F(Q,N,K,I,DNIDX)
	                  ENDDO
	                ENDDO
	              ENDDO
	            ENDDO
	          ENDIF
	        ENDDO
	      ENDIF
	      IF ( DO_UPWELLING ) THEN
	        DO N = 1, NLAYERS
	          IF ( LAYER_VARY_FLAG(N) ) THEN
	            DO Q = 1, LAYER_VARY_NUMBER(N)
	              DO K = N_LAYERSOURCE_UP, NLAYERS
	                DO I = 1, N_USER_STREAMS
	                  DO UA = 1, LOCAL_N_USERAZM
	                    ATMOSWF_MSCATSTERM(Q,N,K,I,UPIDX,UA) =
     $                      ATMOSWF_MSCATSTERM_F(Q,N,K,I,UPIDX)
	                  ENDDO
	                ENDDO
	              ENDDO
	              DO I = 1, N_USER_STREAMS
	                DO UA = 1, LOCAL_N_USERAZM
	                  ATMOSWF_MSCATBOA_STERM(Q,N,I,UA) =
     $                    ATMOSWF_MSCATBOA_STERM_F(Q,N,I)
	                  ATMOSWF_DIRECTBOA_STERM(Q,N,I,UA) =
     $                    ATMOSWF_DIRECTBOA_STERM_F(Q,N,I)
	                ENDDO
	              ENDDO
	            ENDDO
	          ENDIF
	        ENDDO
	      ENDIF
	    ENDIF

C  end atmospheric WF clause

	  ENDIF

C  albedo weighting functions
C  --------------------------

	  IF ( DO_ALBEDO_LINEARIZATION ) THEN

C  Full output

	    DO IDIR = 1, N_DIRECTIONS
	      WDIR = WHICH_DIRECTIONS(IDIR)
	      DO UT = 1, N_OUT_USERTAUS
	        DO I = 1, N_OUT_STREAMS
	          DO UA = 1, LOCAL_N_USERAZM
	            ALBEDOWF(UT,I,WDIR,UA) = ALBEDOWF_F(UT,I,WDIR)
	          ENDDO
	        ENDDO
	      ENDDO
	    ENDDO

C  M-SCAT source term weighting functions (copy Fourier component)

	    IF ( SAVE_LAYER_MSST ) THEN
	      IF ( DO_DNWELLING ) THEN
	        DO K = 1, N_LAYERSOURCE_DN
	          DO I = 1, N_USER_STREAMS
	            DO UA = 1, LOCAL_N_USERAZM
C problematic line!
                      write(*,*) 'check'
                      ALBEDOWF_MSCATSTERM(K,I,DNIDX,UA) =
     $                ALBEDOWF_MSCATSTERM_F(K,I,DNIDX)
	            ENDDO
	          ENDDO
	        ENDDO
	      ENDIF
	      IF ( DO_UPWELLING ) THEN
	        DO K = N_LAYERSOURCE_UP, NLAYERS
	          DO I = 1, N_USER_STREAMS
	            DO UA = 1, LOCAL_N_USERAZM
	              ALBEDOWF_MSCATSTERM(K,I,UPIDX,UA) =
     $                ALBEDOWF_MSCATSTERM_F(K,I,UPIDX)
	            ENDDO
	          ENDDO
	        ENDDO
	        DO I = 1, N_USER_STREAMS
	          DO UA = 1, LOCAL_N_USERAZM
	            ALBEDOWF_MSCATBOA_STERM(I,UA) =
     $              ALBEDOWF_MSCATBOA_STERM_F(I)
	          ENDDO
	        ENDDO
	      ENDIF
	    ENDIF

C  End albedo WF clause

	  ENDIF

C  If no_azimuth, then exit
C  ------------------------

	  IF ( DO_NO_AZIMUTH ) RETURN

C  ######################
C  Fourier components > 0
C  ######################

	ELSE

C  atmospheric weighting functions
C  -------------------------------

	  IF ( DO_LAYER_LINEARIZATION ) THEN

C  Full output

	    DO N = 1, NLAYERS
	      IF ( LAYER_VARY_FLAG(N) ) THEN
	        DO Q = 1, LAYER_VARY_NUMBER(N)
	          DO UA = 1, LOCAL_N_USERAZM
	            DO IDIR = 1, N_DIRECTIONS
	              WDIR = WHICH_DIRECTIONS(IDIR)
	              DO UT = 1, N_OUT_USERTAUS
	                DO I = 1, N_OUT_STREAMS
	                  ATMOSWF(Q,N,UT,I,WDIR,UA) =
     $                    ATMOSWF(Q,N,UT,I,WDIR,UA) + 
     $                    ATMOSWF_F(Q,N,UT,I,WDIR)*AZMFAC(UA)
	                ENDDO
	              ENDDO
	            ENDDO
	          ENDDO
	        ENDDO
	      ENDIF
	    ENDDO

C  M_SCAT source term atmospheric weighting functions (add Fourier component)

	    IF ( SAVE_LAYER_MSST ) THEN
	      IF ( DO_DNWELLING ) THEN
	        DO N = 1, NLAYERS
	          IF ( LAYER_VARY_FLAG(N) ) THEN
	            DO Q = 1, LAYER_VARY_NUMBER(N)
	              DO UA = 1, LOCAL_N_USERAZM
	                DO K = 1, N_LAYERSOURCE_DN
	                  DO I = 1, N_USER_STREAMS
	                    ATMOSWF_MSCATSTERM(Q,N,K,I,DNIDX,UA) =
     $                  ATMOSWF_MSCATSTERM(Q,N,K,I,DNIDX,UA) + 
     $                  ATMOSWF_MSCATSTERM_F(Q,N,K,I,DNIDX)*AZMFAC(UA)
	                  ENDDO
	                ENDDO
	              ENDDO
	            ENDDO
	          ENDIF
	        ENDDO
	      ENDIF
	      IF ( DO_UPWELLING ) THEN
	        DO N = 1, NLAYERS
	          IF ( LAYER_VARY_FLAG(N) ) THEN
	            DO Q = 1, LAYER_VARY_NUMBER(N)
	              DO UA = 1, LOCAL_N_USERAZM
	                DO K = N_LAYERSOURCE_UP, NLAYERS
	                  DO I = 1, N_USER_STREAMS
	                    ATMOSWF_MSCATSTERM(Q,N,K,I,UPIDX,UA) =
     $                  ATMOSWF_MSCATSTERM(Q,N,K,I,UPIDX,UA) + 
     $                  ATMOSWF_MSCATSTERM_F(Q,N,K,I,UPIDX)*AZMFAC(UA)
	                  ENDDO
	                ENDDO
	                DO I = 1, N_USER_STREAMS
	                  ATMOSWF_MSCATBOA_STERM(Q,N,I,UA) =
     $                    ATMOSWF_MSCATBOA_STERM(Q,N,I,UA) + 
     $                    ATMOSWF_MSCATBOA_STERM_F(Q,N,I)*AZMFAC(UA)
	                  ATMOSWF_DIRECTBOA_STERM(Q,N,I,UA) =
     $                    ATMOSWF_DIRECTBOA_STERM(Q,N,I,UA) + 
     $                    ATMOSWF_DIRECTBOA_STERM_F(Q,N,I)*AZMFAC(UA)
	                ENDDO
	              ENDDO
	            ENDDO
	          ENDIF
	        ENDDO
	      ENDIF
	    ENDIF

C  End atmospheric WF clause

	  ENDIF

C  albedo weighting functions (non-Lambertian only)
C  ------------------------------------------------

	  IF ( DO_ALBEDO_LINEARIZATION ) THEN
	    IF ( .NOT. DO_LAMBERTIAN_ALBEDO ) THEN

C  Full output weighting functions

	      DO UA = 1, LOCAL_N_USERAZM
	        DO IDIR = 1, N_DIRECTIONS
	          WDIR = WHICH_DIRECTIONS(IDIR)
	          DO UT = 1, N_OUT_USERTAUS
	            DO I = 1, N_OUT_STREAMS
	              ALBEDOWF(UT,I,WDIR,UA) =
     &                ALBEDOWF(UT,I,WDIR,UA) +
     &                ALBEDOWF_F(UT,I,WDIR)*AZMFAC(UA)
	            ENDDO
	          ENDDO
	        ENDDO
	      ENDDO

C  M-SCAT source term atmospheric weighting functions (add Fourier component)

	      IF ( SAVE_LAYER_MSST ) THEN
	        IF ( DO_DNWELLING ) THEN
	          DO UA = 1, LOCAL_N_USERAZM
	            DO K = 1, N_LAYERSOURCE_DN
	              DO I = 1, N_USER_STREAMS
	                ALBEDOWF_MSCATSTERM(K,I,DNIDX,UA) =
     $                  ALBEDOWF_MSCATSTERM(K,I,DNIDX,UA) + 
     $                  ALBEDOWF_MSCATSTERM_F(K,I,DNIDX)*AZMFAC(UA)
	              ENDDO
	            ENDDO
	          ENDDO
	        ENDIF
	        IF ( DO_UPWELLING ) THEN
	          DO UA = 1, LOCAL_N_USERAZM
	            DO K = N_LAYERSOURCE_UP, NLAYERS
	              DO I = 1, N_USER_STREAMS
	                ALBEDOWF_MSCATSTERM(K,I,UPIDX,UA) =
     $                  ALBEDOWF_MSCATSTERM(K,I,UPIDX,UA) + 
     $                  ALBEDOWF_MSCATSTERM_F(K,I,UPIDX)*AZMFAC(UA)
	              ENDDO
	            ENDDO
	            DO I = 1, N_USER_STREAMS
	              ALBEDOWF_MSCATBOA_STERM(I,UA) =
     $                ALBEDOWF_MSCATBOA_STERM(I,UA) + 
     $                ALBEDOWF_MSCATBOA_STERM_F(I)*AZMFAC(UA)
	            ENDDO
	          ENDDO
	        ENDIF
	      ENDIF

C  End albedo WF clause

	    ENDIF
	  ENDIF

C  end Fourier clause

	ENDIF

C  Finish

	END
