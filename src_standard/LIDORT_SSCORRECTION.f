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

	SUBROUTINE LIDORT_SSCORRECTION (SSFLUX)

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'
	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  Include file of single scatter saved variables

	INCLUDE '../include_s/LIDORT_SSCORRECTION.VARS'

C  Include file of results

	INCLUDE '../include_s/LIDORT_RESULTS.VARS'

C  Argument

	DOUBLE PRECISION SSFLUX

C  local variables
C  ---------------

	INTEGER		 N, NUT, NSTART, NUT_PREV, NLEVEL
	INTEGER		 UT, L, UTA, UM, IUM, IA, NC
	DOUBLE PRECISION FINAL_SOURCE, HELP, SS_LAYERSOURCE
	DOUBLE PRECISION MULTIPLIER, SSCORRECTION, LEGPOLY, BEFORE
	DOUBLE PRECISION CTHETA, STHETA, CPHI(MAX_USER_RELAZMS), COSSCAT
	DOUBLE PRECISION SALPHA(MAX_USER_STREAMS)
	DOUBLE PRECISION CALPHA(MAX_USER_STREAMS)
	DOUBLE PRECISION DF1(0:MAXMOMENT),DF2(0:MAXMOMENT)

C  old code
C        DOUBLE PRECISION CFPLG (0:MAXMOMENT)

C  Set up operations
C  -----------------

C  Floating point numbers for Legendre polynomials

	DO L = 2, NMOMENTS_INPUT
	  HELP = DFLOAT(L)
	  DF1(L) = DFLOAT(2*L-1)/HELP
	  DF2(L) = DFLOAT(L-1)/HELP
	ENDDO

C  Cosines of relative azimuth

	DO IA = 1, N_USER_RELAZMS
	  CPHI(IA) = DCOS ( USER_RELAZMS(IA) * DEG_TO_RAD )
	ENDDO

C  Create TMS factors, these get stored

	DO N = 1, NLAYERS
	  IF ( DO_DELTAM_SCALING ) THEN
	    HELP   = ONE - TRUNC_FACTOR(N) * OMEGA_TOTAL_INPUT(N)
	    TMS(N) = OMEGA_TOTAL_INPUT(N) / HELP
	  ELSE
	    TMS(N) = OMEGA_TOTAL_INPUT(N)
	  ENDIF
	ENDDO

C  save some geometrical quantitie

	CTHETA = X0
	STHETA = DSQRT ( ONE - X0 * X0 )
	DO UM = 1, N_USER_STREAMS
	  CALPHA(UM) = USER_STREAMS(UM)
	  SALPHA(UM) = DSQRT ( ONE - CALPHA(UM) * CALPHA(UM) )
	ENDDO

C  ####################
C  #    UPWELLING     #
C  ####################

	IF ( DO_UPWELLING ) THEN

C  ===================================
C  exact scattering values (upwelling)
C  ===================================

C  Loop over output viewing directions

	  DO UM = 1, N_USER_STREAMS
	    DO IA = 1, N_USER_RELAZMS

C  cosins scatter angle (this is valid only for non-refracting atmosphere)

	      COSSCAT = - CTHETA * CALPHA(UM) +
     &                    STHETA * SALPHA(UM) * CPHI(IA)

C  Get the Legendre polynomials, save them

	      SS_PLEG_UP(UM,IA,0) = ONE
	      SS_PLEG_UP(UM,IA,1) = COSSCAT
	      DO L = 2, NMOMENTS_INPUT
	        SS_PLEG_UP(UM,IA,L) =
     &             DF1(L) * SS_PLEG_UP(UM,IA,L-1) * COSSCAT  -
     &             DF2(L) * SS_PLEG_UP(UM,IA,L-2)
	      ENDDO

C  Get the total phase function

	      DO N = 1, NLAYERS
	        IF ( STERM_LAYERMASK_UP(N)) THEN      
	          HELP = ZERO
	          DO L = 0, NMOMENTS_INPUT
	            LEGPOLY = SS_PLEG_UP(UM,IA,L)
	            HELP = HELP + PHASMOMS_TOTAL_INPUT(L,N)*LEGPOLY
	          ENDDO
	          EXACTSCAT_UP(UM,IA,N) = HELP
	        ENDIF
	      ENDDO

C  old code
C    	      CALL CFPLGARR(MAXMOMENT,NMOMENTS_INPUT,0,COSSCAT,CFPLG)
C	      DO N = 1, NLAYERS
C	        IF ( STERM_LAYERMASK_UP(N)) THEN      
C	          HELP = ZERO
C	          DO L = 0, NMOMENTS_INPUT
C	            HELP = HELP + PHASMOMS_TOTAL_INPUT(L,N)*CFPLG(L)
C	          ENDDO
C	          EXACTSCAT_UP(UM,IA,N) = HELP
C	        ENDIF
C	      ENDDO

C  end output viewing directions

	    ENDDO
	  ENDDO

C  add TMS correction factor

	  DO N = 1, NLAYERS
	    IF ( STERM_LAYERMASK_UP(N)) THEN
	      DO UM = 1, N_USER_STREAMS
	        DO IA = 1, N_USER_RELAZMS
	          EXACTSCAT_UP(UM,IA,N) = EXACTSCAT_UP(UM,IA,N) * TMS(N)
	        ENDDO
	      ENDDO
	    ENDIF
	  ENDDO

C  ===================================
C  Upwelling single scatter recurrence
C  ===================================

C  initialize cumulative source term

	  NC = 0
	  DO UM = 1, N_USER_STREAMS
	    DO IA = 1, N_USER_RELAZMS
	      SS_CUMSOURCE_UP(UM,IA,NC) = ZERO
	    ENDDO
	  ENDDO

C  initialise optical depth loop

 	  NSTART = NLAYERS
	  NUT_PREV = NSTART + 1

C  Main loop over all output optical depths
C  ========================================

	  DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

	    NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

C  finishing layer

	    NUT = NLEVEL + 1

C  Cumulative single scatter source terms to layer NUT
C  ---------------------------------------------------

C    1. Get layer source terms = Exact scattering * Multiplier
C    2. Loop over layers working upwards to NUT

	    DO N = NSTART, NUT, -1
	      NC = NLAYERS + 1 - N
	      DO UM = 1, N_USER_STREAMS
	        MULTIPLIER = EMULT_UP(UM,N)
	        DO IA = 1, N_USER_RELAZMS
	          SS_LAYERSOURCE = EXACTSCAT_UP(UM,IA,N) * MULTIPLIER
	          SS_CUMSOURCE_UP(UM,IA,NC) = SS_LAYERSOURCE +
     &              T_DELT_USERM(N,UM)*SS_CUMSOURCE_UP(UM,IA,NC-1)
	        ENDDO
	      ENDDO
	    ENDDO

C  Offgrid output
C  --------------

C  add additional partial layer source term = Exact Scat * Multiplier
C  Set final cumulative source and Correct the intensity

	    IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

	      UT = OFFGRID_UTAU_OUTINDEX(UTA)
	      N  = OFFGRID_UTAU_LAYERIDX(UT)
	      DO UM = 1, N_USER_STREAMS
	        MULTIPLIER = UT_EMULT_UP(UM,UT)
	        IUM = USEROUTPUT_INDEX(UM)
	        DO IA = 1, N_USER_RELAZMS
	          SS_LAYERSOURCE = EXACTSCAT_UP(UM,IA,N) * MULTIPLIER
	          FINAL_SOURCE = SS_LAYERSOURCE +
     &                 T_UTUP_USERM(UT,UM)*SS_CUMSOURCE_UP(UM,IA,NC)
	          SSCORRECTION = SSFLUX * FINAL_SOURCE
	          INTENSITY(UTA,IUM,UPIDX,IA) =
     &              INTENSITY(UTA,IUM,UPIDX,IA) + SSCORRECTION
	        ENDDO
	      ENDDO

C  Ongrid output
C  -------------

C  Set final cumulative source and Correct the intensity

	    ELSE

	      DO UM = 1, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        DO IA = 1, N_USER_RELAZMS
	          FINAL_SOURCE = SS_CUMSOURCE_UP(UM,IA,NC)
	          SSCORRECTION = SSFLUX * FINAL_SOURCE
	          BEFORE = INTENSITY(UTA,IUM,UPIDX,IA)
	          INTENSITY(UTA,IUM,UPIDX,IA) =
     &                          BEFORE + SSCORRECTION
	        ENDDO
	      ENDDO

	    ENDIF

C  Check for updating the recursion

	    IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
	    NUT_PREV = NUT

C  end loop over optical depth

	  ENDDO

C  end Upwelling clause

	ENDIF

C  ####################
C  #   DOWNWELLING    #
C  ####################

	IF ( DO_DNWELLING ) THEN

C  exact scattering values (downwelling)
C  -------------------------------------

C  Loop over output viewing directions

	  DO UM = 1, N_USER_STREAMS
	    DO IA = 1, N_USER_RELAZMS

C  cosins scatter angle (this is valid only for non-refracting atmosphere)

	      COSSCAT = + CTHETA * CALPHA(UM) +
     &                    STHETA * SALPHA(UM) * CPHI(IA)

C  Get the Legendre polynomials, save them

	      SS_PLEG_DN(UM,IA,0) = ONE
	      SS_PLEG_DN(UM,IA,1) = COSSCAT
	      DO L = 2, NMOMENTS_INPUT
	        SS_PLEG_DN(UM,IA,L) =
     &             DF1(L) * SS_PLEG_DN(UM,IA,L-1) * COSSCAT  -
     &             DF2(L) * SS_PLEG_DN(UM,IA,L-2)
	      ENDDO

C  Get the total phase function

	      DO N = 1, NLAYERS
	        IF ( STERM_LAYERMASK_DN(N)) THEN      
	          HELP = ZERO
	          DO L = 0, NMOMENTS_INPUT
	            LEGPOLY = SS_PLEG_DN(UM,IA,L)
	            HELP = HELP + PHASMOMS_TOTAL_INPUT(L,N)*LEGPOLY
	          ENDDO
	          EXACTSCAT_DN(UM,IA,N) = HELP
	        ENDIF
	      ENDDO

C  old code
C    	      CALL CFPLGARR(MAXMOMENT,NMOMENTS_INPUT,0,COSSCAT,CFPLG)
C	      DO N = 1, NLAYERS
C	        IF ( STERM_LAYERMASK_DN(N)) THEN      
C	          HELP = ZERO
C	          DO L = 0, NMOMENTS_INPUT
C	            HELP = HELP + PHASMOMS_TOTAL_INPUT(L,N)*CFPLG(L)
C	          ENDDO
C	          EXACTSCAT_DN(UM,IA,N) = HELP
C	        ENDIF
C	      ENDDO

C  end output viewing directions

	    ENDDO
	  ENDDO

C  add TMS correction factor

	  DO N = 1, NLAYERS
	    IF ( STERM_LAYERMASK_DN(N)) THEN
	      DO UM = 1, N_USER_STREAMS
	        DO IA = 1, N_USER_RELAZMS
	          EXACTSCAT_DN(UM,IA,N) = EXACTSCAT_DN(UM,IA,N) * TMS(N)
	        ENDDO
	      ENDDO
	    ENDIF
	  ENDDO

C  initialize cumulative source term

	  NC = 0
	  DO UM = 1, N_USER_STREAMS
	    DO IA = 1, N_USER_RELAZMS
	      SS_CUMSOURCE_DN(UM,IA,NC) = ZERO
	    ENDDO
	  ENDDO

C  initialise optical depth loop

 	  NSTART = 1
	  NUT_PREV = NSTART - 1

C  Main loop over all output optical depths
C  ========================================

	  DO UTA = 1, N_OUT_USERTAUS

C  Layer index for given optical depth

	    NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

C  finishing layer

	    NUT = NLEVEL

C  Cumulative single scatter source terms to layer NUT
C  ---------------------------------------------------

C    1. Get layer source terms = Exact scattering * Multiplier
C    2. Loop over layers working upwards to NUT

	    DO N = NSTART, NUT
	      NC = N
	      DO UM = 1, N_USER_STREAMS
	        MULTIPLIER = EMULT_DN(UM,N)
	        DO IA = 1, N_USER_RELAZMS
	          SS_LAYERSOURCE = EXACTSCAT_DN(UM,IA,N) * MULTIPLIER
	          SS_CUMSOURCE_DN(UM,IA,NC) = SS_LAYERSOURCE +
     &               T_DELT_USERM(N,UM)*SS_CUMSOURCE_DN(UM,IA,NC-1)
	        ENDDO
	      ENDDO
	    ENDDO

C  Offgrid output
C  --------------

C  add additional partial layer source term = Exact Scat * Multiplier
C  Set final cumulative source and Correct the intensity

	    IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

	      UT = OFFGRID_UTAU_OUTINDEX(UTA)
	      N  = OFFGRID_UTAU_LAYERIDX(UT)
	      DO UM = 1, N_USER_STREAMS
	        MULTIPLIER = UT_EMULT_DN(UM,UT)
	        IUM = USEROUTPUT_INDEX(UM)
	        DO IA = 1, N_USER_RELAZMS
	          SS_LAYERSOURCE = EXACTSCAT_DN(UM,IA,N) * MULTIPLIER
	          FINAL_SOURCE = SS_LAYERSOURCE +
     &               T_UTDN_USERM(UT,UM)*SS_CUMSOURCE_DN(UM,IA,NC)
	          SSCORRECTION = SSFLUX * FINAL_SOURCE
c	write(*,*)INTENSITY(UTA,IUM,DNIDX,IA),SSCORRECTION
	          INTENSITY(UTA,IUM,DNIDX,IA) =
     &              INTENSITY(UTA,IUM,DNIDX,IA) + SSCORRECTION
	        ENDDO
	      ENDDO

C  Ongrid output
C  -------------

C  Set final cumulative source and Correct the intensity

	    ELSE

	      DO UM = 1, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        DO IA = 1, N_USER_RELAZMS
	          FINAL_SOURCE = SS_CUMSOURCE_DN(UM,IA,NC)
	          SSCORRECTION = SSFLUX * FINAL_SOURCE
	          INTENSITY(UTA,IUM,DNIDX,IA) =
     &                INTENSITY(UTA,IUM,DNIDX,IA) + SSCORRECTION
c	write(*,*)INTENSITY(UTA,IUM,DNIDX,IA),SSCORRECTION
	        ENDDO
	      ENDDO

	    ENDIF

C  Check for updating the recursion

	    IF ( NUT. NE. NUT_PREV ) NSTART = NUT + 1
	    NUT_PREV = NUT

C  end loop over optical depth

	  ENDDO

C  end Downwelling clause

	ENDIF

C  Finish

	END
