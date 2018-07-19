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

	SUBROUTINE UPUSER_INTENSITY
     I    ( DO_INCLUDE_ALBEDO,   DO_INCLUDE_SURFEMISS,
     I      DO_INCLUDE_MVOUTPUT, DO_INCLUDE_DIRECTBEAM,
     I      R2, F1, FOURIER,
     I      USER_DIRECT_BEAM, DO_GMULT )

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_s/LIDORT_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Flux factor

	DOUBLE PRECISION F1

C  local control flags

	LOGICAL		 DO_INCLUDE_ALBEDO
	LOGICAL		 DO_INCLUDE_DIRECTBEAM
	LOGICAL		 DO_INCLUDE_SURFEMISS
	LOGICAL		 DO_INCLUDE_MVOUTPUT

C  albedo factor

	DOUBLE PRECISION R2

C  Fourier number

	INTEGER		 FOURIER

C  reflected direct beam

	DOUBLE PRECISION USER_DIRECT_BEAM(MAX_USER_STREAMS)

C  Green function multiplier flag (Quadrature solutions only)

	LOGICAL		 DO_GMULT(MAX_OFFGRID_USERTAUS)

C  local variables
C  ---------------

	INTEGER		 N, NUT, NSTART, NUT_PREV, NLEVEL, NC
	INTEGER		 UT, UTA, UM, IUM, AA, I, IQD
	DOUBLE PRECISION LAYER_SOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION MSCAT_LAYERSOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION FINAL_SOURCE, HELP, TAU_DN
	LOGICAL		 DOCONS(MAXLAYER)

C  Zero all Fourier components close to zenith

	IF ( DO_USER_STREAMS ) THEN
	  IF ( LOCAL_UM_START .GT. 1 ) THEN
	    DO UTA = 1, N_OUT_USERTAUS
	      DO UM = 1, LOCAL_UM_START - 1
	        IUM = USEROUTPUT_INDEX(UM)
	        INTENSITY_F(UTA,IUM,UPIDX) = ZERO
	      ENDDO
	    ENDDO
	  ENDIF
	ENDIF

C  Initialize recursion for user-defined stream angles only

	IF ( DO_USER_STREAMS ) THEN

	  CALL GET_BOASOURCE
     I    ( DO_INCLUDE_ALBEDO,   DO_INCLUDE_SURFEMISS,
     I      DO_INCLUDE_DIRECTBEAM, R2, FOURIER,
     I      USER_DIRECT_BEAM )

	  NC = 0
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    CUMSOURCE_UP(UM,NC) = BOA_SOURCE(UM) + DIRECT_BOA_SOURCE(UM)
	  ENDDO

C  Save separate terms if MSST output flagged

	  IF ( SAVE_LAYER_MSST ) THEN
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      MSCATBOA_SOURCETERM_F(UM)  = F1 * BOA_SOURCE(UM)
	      DIRECTBOA_SOURCETERM_F(UM) = F1 * DIRECT_BOA_SOURCE(UM)
	    ENDDO
	  ENDIF

	ENDIF

C  initialise cumulative source term loop

 	NSTART = NLAYERS
	NUT_PREV = NSTART + 1
	DO N = 1, NLAYERS
	  DOCONS(N) = .TRUE.
	  HMULT_EXIST(N) = .FALSE.
	ENDDO

C  loop over all output optical depths
C  -----------------------------------

	DO UTA = N_OUT_USERTAUS, 1, -1

C  Layer index for given optical depth

	  NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)
C    1. Get layer source terms
C    2. Find cumulative source term
C    3. Set multiple scatter source term (MSST) output if flagged

	  IF ( DO_USER_STREAMS ) THEN
	    NUT = NLEVEL + 1
	    DO N = NSTART, NUT, -1
	      NC = NLAYERS + 1 - N
	      CALL WHOLELAYER_STERM_UP
     &          ( N, LAYER_SOURCE, MSCAT_LAYERSOURCE, DOCONS(N) )
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        CUMSOURCE_UP(UM,NC) = LAYER_SOURCE(UM) +
     &                   T_DELT_USERM(N,UM)*CUMSOURCE_UP(UM,NC-1)
	      ENDDO
	      IF ( SAVE_LAYER_MSST ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          MSCATSTERM_F(N,UM,UPIDX) = MSCAT_LAYERSOURCE(UM) * F1
	        ENDDO
	      ENDIF
	    ENDDO
	  ENDIF

C  Offgrid output
C  --------------

	  IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

	    UT = OFFGRID_UTAU_OUTINDEX(UTA)
	    N  = OFFGRID_UTAU_LAYERIDX(UT)

C  Transmittances (saved results)
  
	    IF ( DO_GMULT(UT) ) THEN
	      TAU_DN = OFFGRID_UTAU_VALUES(UT)
	      DO AA = 1, NSTREAMS
	        HELP = KEIGEN(AA,N) * TAU_DN
	        IF ( HELP .GT. MAX_TAU_QPATH ) THEN
	          T_UTDN_EIGEN(AA,UT) = ZERO
	        ELSE
	          T_UTDN_EIGEN(AA,UT) = DEXP(-HELP)
	        ENDIF
	        HELP = KEIGEN(AA,N) * ( DELTA(N) - TAU_DN )
	        IF ( HELP .GT. MAX_TAU_QPATH ) THEN
	          T_UTUP_EIGEN(AA,UT) = ZERO
	        ELSE
	          T_UTUP_EIGEN(AA,UT) = DEXP(-HELP)
	        ENDIF
	      ENDDO
	    ENDIF

C  Quadrature output at offgrid optical depths
C  ( Required if mean-value calculations are to be done)
C    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )
C  Copy into storage if quad output required

	    IF ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT ) THEN
	      CALL QUADINTENS_OFFGRID_UP
     &               ( N, UTA, UT, F1, DO_GMULT(UT) )
	    ENDIF
	    IF ( DO_QUAD_OUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        IQD = QUADOUTPUT_INDEX(I)
 	        INTENSITY_F(UTA,IQD,UPIDX) = QUADINTENS(UTA,I,UPIDX)
	      ENDDO
	    ENDIF

C  User-defined stream output, add additional partial layer source term

	    IF ( DO_USER_STREAMS ) THEN
	      CALL PARTLAYER_STERM_UP
     &              ( N, UT, LAYER_SOURCE, DOCONS(N) )
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        FINAL_SOURCE = LAYER_SOURCE(UM) +
     &                       T_UTUP_USERM(UT,UM)*CUMSOURCE_UP(UM,NC)
	        INTENSITY_F(UTA,IUM,UPIDX) = F1 * FINAL_SOURCE
	      ENDDO
	    ENDIF

C  set multiplier flag

	    DO_GMULT(UT) = .FALSE.

C  Ongrid output
C  -------------

	  ELSE

C  Quadrature output at layer boundaries
C  ( Required if mean-value calculations are to be done)
C    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )
C  Copy into storage if quad output required

	    IF ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT ) THEN
	      CALL QUADINTENS_LEVEL_UP ( NLEVEL, UTA, F1 )
	    ENDIF
	    IF ( DO_QUAD_OUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        IQD = QUADOUTPUT_INDEX(I)
 	        INTENSITY_F(UTA,IQD,UPIDX) = QUADINTENS(UTA,I,UPIDX)
	      ENDDO
	    ENDIF

C  User-defined stream output, just set to the cumulative source term

	    IF ( DO_USER_STREAMS ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        INTENSITY_F(UTA,IUM,UPIDX) = F1 * CUMSOURCE_UP(UM,NC)
	      ENDDO
	    ENDIF

	  ENDIF

C  Check for updating the recursion

	  IF ( DO_USER_STREAMS ) THEN
	    IF ( NUT. NE. NUT_PREV ) NSTART = NUT - 1
	    NUT_PREV = NUT
	  ENDIF

C  end loop over optical depth

	ENDDO

C  Finish

	END

C

	SUBROUTINE DNUSER_INTENSITY
     &        ( DO_INCLUDE_MVOUTPUT, F1, DO_GMULT )

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of setup and solution variables (input)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_s/LIDORT_RESULTS.VARS'

C  Subroutine input argument (flux factor constant)

	DOUBLE PRECISION F1

C  Green function multiplier flag (Quadrature solutions only)

	LOGICAL		 DO_GMULT(MAX_OFFGRID_USERTAUS)

C  local inclusion flag

	LOGICAL		 DO_INCLUDE_MVOUTPUT

C  local variables

	INTEGER		 N, NUT, NSTART, NUT_PREV, NLEVEL
	INTEGER		 UT, UTA, UM, IUM, AA, NC, I, IQD
	DOUBLE PRECISION TOA_SOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION LAYER_SOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION MSCAT_LAYERSOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION FINAL_SOURCE, HELP, TAU_DN
	LOGICAL		 DOCONS(MAXLAYER)

C  Zero all Fourier components close to zenith

	IF ( DO_USER_STREAMS ) THEN
	  IF ( LOCAL_UM_START .GT. 1 ) THEN
	    DO UTA = 1, N_OUT_USERTAUS
	      DO UM = 1, LOCAL_UM_START - 1
	        IUM = USEROUTPUT_INDEX(UM)
	        INTENSITY_F(UTA,IUM,DNIDX) = ZERO
	      ENDDO
	    ENDDO
	  ENDIF
	ENDIF

C  Initialize recursion for user-defined stream angles only

	IF ( DO_USER_STREAMS ) THEN
	  CALL GET_TOASOURCE ( TOA_SOURCE )
	  NC = 0
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    CUMSOURCE_DN(UM,NC) = TOA_SOURCE(UM)
	  ENDDO
	ENDIF

C  initialise cumulative source term loop

 	NSTART = 1
	NUT_PREV = NSTART - 1
	DO N = 1, NLAYERS
	  DOCONS(N) = .TRUE.
	ENDDO

C  loop over all output optical depths
C  -----------------------------------

	DO UTA = 1, N_OUT_USERTAUS

C  Layer index for given optical depth

	  NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)
C    1. Get layer source terms
C    2. Find cumulative source term
C    3. Set multiple scatter source term output if flagged

	  IF ( DO_USER_STREAMS ) THEN
	    NUT = NLEVEL
	    DO N = NSTART, NUT
	      NC = N
	      CALL WHOLELAYER_STERM_DN
     &          ( N, LAYER_SOURCE, MSCAT_LAYERSOURCE, DOCONS(N) )
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        CUMSOURCE_DN(UM,NC) = LAYER_SOURCE(UM) +
     &                        T_DELT_USERM(N,UM)*CUMSOURCE_DN(UM,NC-1)
	      ENDDO
	      IF ( SAVE_LAYER_MSST ) THEN
	        DO UM = LOCAL_UM_START, N_USER_STREAMS
	          MSCATSTERM_F(N,UM,DNIDX) = MSCAT_LAYERSOURCE(UM) * F1
	        ENDDO
	      ENDIF
	    ENDDO
	  ENDIF

C  Offgrid output
C  --------------

	  IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN

	    UT = OFFGRID_UTAU_OUTINDEX(UTA)
	    N  = OFFGRID_UTAU_LAYERIDX(UT)

C  Transmittances (saved results)
  
	    IF ( DO_GMULT(UT) ) THEN
	      TAU_DN = OFFGRID_UTAU_VALUES(UT)
	      DO AA = 1, NSTREAMS
	        HELP = KEIGEN(AA,N) * TAU_DN
	        IF ( HELP .GT. MAX_TAU_QPATH ) THEN
	          T_UTDN_EIGEN(AA,UT) = ZERO
	        ELSE
	          T_UTDN_EIGEN(AA,UT) = DEXP(-HELP)
	        ENDIF
	        HELP = KEIGEN(AA,N) * ( DELTA(N) - TAU_DN )
	        IF ( HELP .GT. MAX_TAU_QPATH ) THEN
	          T_UTUP_EIGEN(AA,UT) = ZERO
	        ELSE
	          T_UTUP_EIGEN(AA,UT) = DEXP(-HELP)
	        ENDIF
	      ENDDO
	    ENDIF

C  Quadrature output at offgrid optical depths
C  ( Required if mean-value calculations are to be done)
C    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )
C  Copy into storage if quad output required

	    IF ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT ) THEN
	      CALL QUADINTENS_OFFGRID_DN
     &              ( N, UTA, UT, F1, DO_GMULT(UT) )
	    ENDIF
	    IF ( DO_QUAD_OUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        IQD = QUADOUTPUT_INDEX(I)
 	        INTENSITY_F(UTA,IQD,DNIDX) = QUADINTENS(UTA,I,DNIDX)
	      ENDDO
	    ENDIF

C  User-defined stream output, add additional partial layer source term

	    IF ( DO_USER_STREAMS ) THEN
	      CALL PARTLAYER_STERM_DN
     &              ( N, UT, LAYER_SOURCE, DOCONS(N) )
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        FINAL_SOURCE = LAYER_SOURCE(UM) +
     &                       T_UTDN_USERM(UT,UM)*CUMSOURCE_DN(UM,NC)
	        INTENSITY_F(UTA,IUM,DNIDX) = F1 * FINAL_SOURCE
	      ENDDO
	    ENDIF

C  set multiplier flag

	    DO_GMULT(UT) = .FALSE.

C  Ongrid output
C  -------------

	  ELSE

C  Quadrature output at layer boundaries
C  ( Required if mean-value calculations are to be done)
C    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )
C  Copy into storage if quad output required

	    IF ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT ) THEN
	      CALL QUADINTENS_LEVEL_DN ( NLEVEL, UTA, F1 )
	    ENDIF
	    IF ( DO_QUAD_OUTPUT ) THEN
	      DO I = 1, NSTREAMS
	        IQD = QUADOUTPUT_INDEX(I)
 	        INTENSITY_F(UTA,IQD,DNIDX) = QUADINTENS(UTA,I,DNIDX)
	      ENDDO
	    ENDIF

C  User-defined stream output, just set to the cumulative source term

	    IF ( DO_USER_STREAMS ) THEN
	      DO UM = LOCAL_UM_START, N_USER_STREAMS
	        IUM = USEROUTPUT_INDEX(UM)
	        INTENSITY_F(UTA,IUM,DNIDX) = F1 * CUMSOURCE_DN(UM,NC)
	      ENDDO
	    ENDIF

	  ENDIF

C  Check for updating the recursion

	  IF ( DO_USER_STREAMS ) THEN
	    IF ( NUT. NE. NUT_PREV ) NSTART = NUT + 1
	    NUT_PREV = NUT
	  ENDIF

C  end loop over optical depth

	ENDDO

C  Finish

	END

C

	SUBROUTINE QUADINTENS_LEVEL_UP ( NLEVEL, UTA, F1 )

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model variables (input)

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of main model variables (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_s/LIDORT_RESULTS.VARS'

C  Subroutine input arguments (Flux factor, level)

	INTEGER		 NLEVEL, UTA
	DOUBLE PRECISION F1

C  local variables

	INTEGER		 N, I, I1, AA
	DOUBLE PRECISION SPAR, SHOM, HOM1, HOM2

C  For those optical depths at layer boundaries
C  --------------------------------------------

C  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
C  looking at the intensity at the top of these layers. The
C  case where the level mask = NLAYERS is the upwelling intensity
C  at the bottom of the atmosphere (treated separately).

	N = NLEVEL + 1

C  homogeneous and particular solution contributions SHOM and SPAR

C  For the lowest level

	IF ( NLEVEL .EQ. NLAYERS ) THEN

	  DO I = 1, NSTREAMS
	    I1 = I + NSTREAMS
	    SHOM = ZERO
	    DO AA = 1, NSTREAMS
	      HOM1 = LCON_XVEC(I1,AA,NLEVEL) * T_DELT_EIGEN(AA,NLEVEL)
	      HOM2 = MCON_XVEC(I1,AA,NLEVEL)
	      SHOM = SHOM + HOM1 + HOM2
	    ENDDO
            SPAR = WLOWER(I1,NLEVEL)
 	    QUADINTENS(UTA,I,UPIDX) = F1 * ( SPAR + SHOM )
	  ENDDO

C  For other levels in the atmosphere

	ELSE

	  DO I = 1, NSTREAMS
	    I1 = I + NSTREAMS
	    SHOM = ZERO
	    DO AA = 1, NSTREAMS
	      HOM1 = LCON_XVEC(I1,AA,N)
	      HOM2 = MCON_XVEC(I1,AA,N) * T_DELT_EIGEN(AA,N)
	      SHOM = SHOM + HOM1 + HOM2
	    ENDDO
            SPAR = WUPPER(I1,N)
 	    QUADINTENS(UTA,I,UPIDX) = F1 * ( SPAR + SHOM )
	  ENDDO

	ENDIF

C  Finish

	END

C

	SUBROUTINE QUADINTENS_LEVEL_DN ( NLEVEL, UTA, F1 )

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model variables (input)

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of main model variables (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_s/LIDORT_RESULTS.VARS'

C  Subroutine input arguments (Flux factor, level)

	INTEGER		 NLEVEL, UTA
	DOUBLE PRECISION F1

C  local variables

	INTEGER		 N, I, AA
	DOUBLE PRECISION SPAR, SHOM, HOM1, HOM2

C  For those optical depths at layer boundaries
C  --------------------------------------------

	N = NLEVEL

C  Downwelling radiation at TOA ( or N = 0 ) is zero

	IF ( NLEVEL .EQ. 0 ) THEN

	  DO I = 1, NSTREAMS
	    QUADINTENS(UTA,I,DNIDX) = ZERO
	  ENDDO

C  Other levels

	ELSE

	  DO I = 1, NSTREAMS
            SPAR = WLOWER(I,N)
	    SHOM = ZERO
	    DO AA = 1, NSTREAMS
	      HOM1 = LCON_XVEC(I,AA,N) * T_DELT_EIGEN(AA,N)
	      HOM2 = MCON_XVEC(I,AA,N)
	      SHOM = SHOM + HOM1 + HOM2
	    ENDDO
	    QUADINTENS(UTA,I,DNIDX) = F1 * ( SPAR + SHOM )
	  ENDDO

	ENDIF

C  Finish

	END

C

	SUBROUTINE QUADINTENS_OFFGRID_UP ( N, UTA, UT, F1, DO_GMULT )

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of main model variables (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_s/LIDORT_RESULTS.VARS'

C  Subroutine input arguments (Flux factor, layer and offgrid index)

	INTEGER		 N, UTA, UT
	DOUBLE PRECISION F1
	LOGICAL		 DO_GMULT

C  local variables

	INTEGER		 I, I1, AA
	DOUBLE PRECISION SPAR, PAR1, PAR2, SHOM, HOM1, HOM2

C  For those optical depths at off-grid levels
C  -------------------------------------------

C  Classical solution

	IF ( DO_CLASSICAL_SOLUTION ) THEN

C  solution

	  DO I = 1, NSTREAMS
	    I1 = I + NSTREAMS
            SPAR = WUPPER(I1,N)*T_UTDN_MUBAR(UT)
	    SHOM = ZERO
	    DO AA = 1, NSTREAMS
	      HOM1 = LCON_XVEC(I1,AA,N) * T_UTDN_EIGEN(AA,UT)
	      HOM2 = MCON_XVEC(I1,AA,N) * T_UTUP_EIGEN(AA,UT)
	      SHOM = SHOM + HOM1 + HOM2
	    ENDDO
	    QUADINTENS(UTA,I,UPIDX) = F1 * ( SPAR + SHOM )
	  ENDDO

C  Green function solution
C  -- Particular solution SPAR first gets Green function multipliers

	ELSE

C  Multipliers

	  CALL QUAD_GFUNCMULT ( N, UT, DO_GMULT )

C  solution

	  DO I = 1, NSTREAMS
	    I1 = I + NSTREAMS
	    SPAR = ZERO
	    DO AA = 1, NSTREAMS
	      PAR1 = XPOS(I,AA,N)  * UT_GMULT_UP(AA,UT)
	      PAR2 = XPOS(I1,AA,N) * UT_GMULT_DN(AA,UT)
	      SPAR = SPAR + PAR1 + PAR2
	    ENDDO
	    SHOM = ZERO
	    DO AA = 1, NSTREAMS
	      HOM1 = LCON_XVEC(I1,AA,N) * T_UTDN_EIGEN(AA,UT)
	      HOM2 = MCON_XVEC(I1,AA,N) * T_UTUP_EIGEN(AA,UT)
	      SHOM = SHOM + HOM1 + HOM2
	    ENDDO
	    QUADINTENS(UTA,I,UPIDX) = F1 * ( SPAR + SHOM )
	  ENDDO

	ENDIF

C  Finish

	END

C

	SUBROUTINE QUADINTENS_OFFGRID_DN ( N, UTA, UT, F1, DO_GMULT )

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of main model variables (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_s/LIDORT_RESULTS.VARS'

C  Subroutine input arguments (Flux factor, layer and offgrid index)

	INTEGER		 N, UTA, UT
	DOUBLE PRECISION F1
	LOGICAL		 DO_GMULT

C  local variables

	INTEGER		 I, I1, AA
	DOUBLE PRECISION SPAR, PAR1, PAR2, SHOM, HOM1, HOM2

C  For those optical depths at off-grid levels
C  -------------------------------------------

C  Classical solution

	IF ( DO_CLASSICAL_SOLUTION ) THEN

C  solution

	  DO I = 1, NSTREAMS
            SPAR = WUPPER(I,N) * T_UTDN_MUBAR(UT)
	    SHOM = ZERO
	    DO AA = 1, NSTREAMS
	      HOM1 = LCON_XVEC(I,AA,N) * T_UTDN_EIGEN(AA,UT)
	      HOM2 = MCON_XVEC(I,AA,N) * T_UTUP_EIGEN(AA,UT)
	      SHOM = SHOM + HOM1 + HOM2
	    ENDDO
	    QUADINTENS(UTA,I,DNIDX) = F1 * ( SPAR + SHOM )
	  ENDDO

C  Green function solution
C  -- Check to see that multipliers have been calculated

	ELSE

C  Multipliers

	  IF ( DO_GMULT ) THEN
	    CALL QUAD_GFUNCMULT ( N, UT, DO_GMULT )
	  ENDIF

C  solution

	  DO I = 1, NSTREAMS
	    I1 = I + NSTREAMS
	    SPAR = ZERO
	    DO AA = 1, NSTREAMS
	      PAR1 = XPOS(I1,AA,N) * UT_GMULT_UP(AA,UT)
	      PAR2 = XPOS(I,AA,N)  * UT_GMULT_DN(AA,UT)
	      SPAR = SPAR + PAR1 + PAR2
	    ENDDO
	    SHOM = ZERO
	    DO AA = 1, NSTREAMS
	      HOM1 = LCON_XVEC(I,AA,N) * T_UTDN_EIGEN(AA,UT)
	      HOM2 = MCON_XVEC(I,AA,N) * T_UTUP_EIGEN(AA,UT)
	      SHOM = SHOM + HOM1 + HOM2
	    ENDDO
	    QUADINTENS(UTA,I,DNIDX) = F1 * ( SPAR + SHOM )
	  ENDDO

	ENDIF

C  Finish

	END

C

	SUBROUTINE LIDORT_MIFLUX_INTENSITY

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include files of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of main model variables (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include files of result variables (module output stored here)

	INCLUDE '../include_s/LIDORT_RESULTS.VARS'

C  local variables

	INTEGER		 I, IDIR, WDIR, UTA, UT, N
	DOUBLE PRECISION SUM_MI, SUM_FX, FMU0
	DOUBLE PRECISION DIRECT_TRANS, DIRECT_FLUX, DIRECT_MEANI

C  mean intensity and flux
C  -----------------------

C  direction loop

	DO IDIR = 1, N_DIRECTIONS

	  WDIR = WHICH_DIRECTIONS(IDIR)

C  loop over all user-defined optical depths

	  DO UTA = 1, N_OUT_USERTAUS

	    SUM_MI = ZERO
	    SUM_FX = ZERO
	    DO I = 1, NSTREAMS
	      SUM_MI = SUM_MI + A(I)  * QUADINTENS(UTA,I,WDIR)
	      SUM_FX = SUM_FX + AX(I) * QUADINTENS(UTA,I,WDIR)
	    ENDDO
	    MEAN_INTENSITY(UTA,WDIR) = SUM_MI * HALF
	    FLUX(UTA,WDIR)           = SUM_FX * PI2

	  ENDDO

C  For the downward direction, add the direct beam contributions

	  IF ( WDIR .EQ. DNIDX ) THEN

C  loop over all the output optical depths

	    DO UTA = 1, N_OUT_USERTAUS

C  For the offgrid values

	      IF ( OFFGRID_UTAU_OUTFLAG(UTA) ) THEN
	        UT = OFFGRID_UTAU_OUTINDEX(UTA)
	        N  = OFFGRID_UTAU_LAYERIDX(UT)

C  Only contributions for layers above the PI cutoff

	        IF ( N .LE. LAYER_PIS_CUTOFF ) THEN
	          DIRECT_TRANS = INITIAL_TRANS(N) * T_UTDN_MUBAR(UT)
	          DIRECT_MEANI = FLUX_FACTOR * DIRECT_TRANS / PI4
		  MEAN_INTENSITY(UTA,WDIR) =
     *                  MEAN_INTENSITY(UTA,WDIR) + DIRECT_MEANI
	          FMU0 = LOCAL_SZA(N) * FLUX_FACTOR
	          DIRECT_FLUX  = FMU0 * DIRECT_TRANS
	          FLUX(UTA,WDIR) = FLUX(UTA,WDIR) + DIRECT_FLUX
	        ENDIF

C  For the on-grid balues

	      ELSE
	        N = UTAU_LEVEL_MASK_DN(UTA)
	        IF ( N .LE. LAYER_PIS_CUTOFF ) THEN
	          IF ( N .EQ. 0 ) THEN
	            DIRECT_TRANS = ONE
	            FMU0 = LOCAL_SZA(1) * FLUX_FACTOR
 	          ELSE
	            DIRECT_TRANS = INITIAL_TRANS(N)*T_DELT_MUBAR(N)
	            FMU0 = LOCAL_SZA(N) * FLUX_FACTOR
	          ENDIF
	          DIRECT_MEANI = FLUX_FACTOR * DIRECT_TRANS / PI4
		  MEAN_INTENSITY(UTA,WDIR) =
     *                  MEAN_INTENSITY(UTA,WDIR) + DIRECT_MEANI
	          DIRECT_FLUX  = FMU0 * DIRECT_TRANS
	          FLUX(UTA,WDIR) = FLUX(UTA,WDIR) + DIRECT_FLUX
	        ENDIF
	      ENDIF

C  Close loops

      	    ENDDO
	  ENDIF   
	ENDDO

C  Finish

	END

C

  	SUBROUTINE GET_TOASOURCE ( TOA_SOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model variables (input)

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  Subroutine arguments

	DOUBLE PRECISION TOA_SOURCE(MAX_USER_STREAMS)

C  local variables

	INTEGER		 UM

C  initialise TOA source function

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  TOA_SOURCE(UM) = ZERO
	ENDDO

C  Finish

	END

C

  	SUBROUTINE GET_BOASOURCE
     I    ( DO_INCLUDE_ALBEDO,   DO_INCLUDE_SURFEMISS,
     I      DO_INCLUDE_DIRECTBEAM, R2, FOURIER,
     I      USER_DIRECT_BEAM )

C  Bottom of the atmosphere source term

C  Include files
C  -------------

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  Include file of Geophysical variables (input surface reflectivities)

	INCLUDE '../include_s/LIDORT_GEOPHYS.VARS'

C  include files of main solution and setup variables (input, output)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  Subroutine input arguments
C  --------------------------

C  local control flags

	LOGICAL		 DO_INCLUDE_ALBEDO
	LOGICAL		 DO_INCLUDE_DIRECTBEAM
	LOGICAL		 DO_INCLUDE_SURFEMISS

C  albedo and Fourier index

	DOUBLE PRECISION R2
	INTEGER		 FOURIER

C  reflected direct beam

	DOUBLE PRECISION USER_DIRECT_BEAM(MAX_USER_STREAMS)

C  local variables
C  ---------------

	INTEGER		 M, N, J, I, UM, AA
	DOUBLE PRECISION DOWN(MAXSTRM), HELP(MAXSTRM)
	DOUBLE PRECISION PAR, HOM, REFLEC

C  Fourier number

	M = FOURIER

C  initialise boa source function

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  BOA_SOURCE(UM) = ZERO
	  DIRECT_BOA_SOURCE(UM) = ZERO
	ENDDO

	N = NLAYERS

C  reflectance from surface
C  ------------------------

	IF ( DO_INCLUDE_ALBEDO ) THEN

C  Downward intensity at computational angles (beam and homog)

	  DO I = 1, NSTREAMS
	    PAR = WLOWER(I,N)
	    HOM = ZERO
	    DO AA = 1, NSTREAMS
	      HOM = HOM + LCON_XVEC(I,AA,N)*T_DELT_EIGEN(AA,N) +
     &                    MCON_XVEC(I,AA,N)
	    ENDDO
	    DOWN(I) = PAR + HOM
	  ENDDO

C  reflectance integrand  a(j).x(j).I(-j)

	  DO I = 1, NSTREAMS
	    HELP(I) = AX(I) * DOWN(I)
	  ENDDO

C  reflected multiple scatter intensity at user defined-angles
C  -----------------------------------------------------------

C  .. integrate reflectance, same for all user-streams in Lambertian case

	  IF ( DO_LAMBERTIAN_ALBEDO ) THEN

	    REFLEC = ZERO
	    DO J = 1, NSTREAMS
	      REFLEC = REFLEC + HELP(J)
	    ENDDO
	    REFLEC = R2 * REFLEC
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	     BOA_SOURCE(UM) = REFLEC
	    ENDDO

C  .. integrate with reflectance function at user angles (non-Lambertian)

	  ELSE

	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      REFLEC = ZERO
	      DO J = 1, NSTREAMS
	        REFLEC = REFLEC + HELP(J) * USER_BIREFLEC(M,UM,J)
	      ENDDO
	      BOA_SOURCE(UM) = R2 * REFLEC
	    ENDDO
	  ENDIF
	  
C  Add direct beam if flagged

	  IF ( DO_INCLUDE_DIRECTBEAM ) THEN
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      DIRECT_BOA_SOURCE(UM) = DIRECT_BOA_SOURCE(UM)
     &                      + USER_DIRECT_BEAM(UM)
	    ENDDO
	  ENDIF

	ENDIF

C  Add surface emission term if flagged

	IF ( DO_INCLUDE_SURFEMISS ) THEN
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DIRECT_BOA_SOURCE(UM) = DIRECT_BOA_SOURCE(UM) + 
     &                 SURFBB*USER_EMISSIVITY(UM)
	  ENDDO
	ENDIF

C  Finish

	END

C

	SUBROUTINE WHOLELAYER_STERM_UP
     &    ( N, LAYERSOURCE, MSCAT_LAYERSOURCE, DO_CONSTANTS_SET )

C  source terms

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of main solution and multiplier variables (input/output)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  arguments

	INTEGER		 N
	DOUBLE PRECISION LAYERSOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION MSCAT_LAYERSOURCE(MAX_USER_STREAMS)
	LOGICAL		 DO_CONSTANTS_SET

C  local variables

	INTEGER		 AA, UM
	DOUBLE PRECISION SPAR, SHOM, SFOR1, SFOR2

C  Save some calculation time

	IF ( DO_CONSTANTS_SET ) THEN
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO AA = 1, NSTREAMS
	      LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
	      MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
  	    ENDDO
	  ENDDO
	ENDIF
	DO_CONSTANTS_SET = .FALSE.

C  Get the homogeneous solution multipliers

	CALL WHOLELAYER_HMULT_UP ( N )

C  Whole layer source function ( Classical solution )

	IF ( DO_CLASSICAL_SOLUTION ) THEN

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    SFOR2 =  EMULT_UP(UM,N) * U_WPOS2(UM,N)
	    SHOM = ZERO
	    DO AA = 1, NSTREAMS
	      SHOM = SHOM + LCON_UXVEC(UM,AA)*HMULT_2(AA,UM,N) +
     &                      MCON_UXVEC(UM,AA)*HMULT_1(AA,UM,N)
	    ENDDO
	    LAYERSOURCE(UM) = SFOR2 + SHOM
	  ENDDO

C  Whole layer source function ( Green's function solution )

	ELSE

	  CALL WHOLELAYER_GMULT_UP ( N )
  
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    SHOM = ZERO
	    DO AA = 1, NSTREAMS
	      SHOM = SHOM + LCON_UXVEC(UM,AA)*HMULT_2(AA,UM,N) +
     &                      MCON_UXVEC(UM,AA)*HMULT_1(AA,UM,N)
	    ENDDO
	    SPAR = ZERO
	    DO AA = 1, NSTREAMS
	      SPAR = SPAR + U_XPOS(UM,AA,N)*SGMULT_UP_DN(AA,UM,N) +
     &                      U_XNEG(UM,AA,N)*SGMULT_UP_UP(AA,UM,N)
	    ENDDO
	    LAYERSOURCE(UM) = SPAR + SHOM
	  ENDDO

	ENDIF

C  Options

C  .. If operating in Ms-mode only, copy multiple scatter term

	IF ( DO_MSMODE_LIDORT ) THEN

	  IF ( SAVE_LAYER_MSST ) THEN
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      MSCAT_LAYERSOURCE(UM) = LAYERSOURCE(UM)
	    ENDDO
	  ENDIF

C  .. Full radiance mode, add single scatter part

	ELSE

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    SFOR1 = U_WPOS1(UM,N) * EMULT_UP(UM,N)
	    LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR1
	  ENDDO

	ENDIF

C  Finish

	END

C

	SUBROUTINE WHOLELAYER_STERM_DN
     &    ( N, LAYERSOURCE, MSCAT_LAYERSOURCE, DO_CONSTANTS_SET )

C  source terms

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of main solution and multiplier variables (input/output)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  arguments

	INTEGER		 N
	DOUBLE PRECISION LAYERSOURCE(MAX_USER_STREAMS)
	DOUBLE PRECISION MSCAT_LAYERSOURCE(MAX_USER_STREAMS)
	LOGICAL		 DO_CONSTANTS_SET

C  local variables

	INTEGER		 AA, UM
	DOUBLE PRECISION SPAR, SHOM, SFOR1, SFOR2

C  Save some calculation time

	IF ( DO_CONSTANTS_SET ) THEN
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO AA = 1, NSTREAMS
	      LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
	      MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
  	    ENDDO
	  ENDDO
	ENDIF

C  Get the homogeneous solution multipliers

	CALL WHOLELAYER_HMULT_DN ( N )

C  Whole layer source function ( Classical solution )

	IF ( DO_CLASSICAL_SOLUTION ) THEN

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    SFOR2 =  EMULT_DN(UM,N) * U_WNEG2(UM,N)
	    SHOM = ZERO
	    DO AA = 1, NSTREAMS
	      SHOM = SHOM + LCON_UXVEC(UM,AA)*HMULT_1(AA,UM,N) +
     &                      MCON_UXVEC(UM,AA)*HMULT_2(AA,UM,N)
	    ENDDO
	    LAYERSOURCE(UM) = SFOR2 + SHOM
	  ENDDO

C  Whole layer source function ( Green's function solution )

	ELSE

	  CALL WHOLELAYER_GMULT_DN ( N )
  
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    SHOM = ZERO
	    DO AA = 1, NSTREAMS
	      SHOM = SHOM + LCON_UXVEC(UM,AA)*HMULT_1(AA,UM,N) +
     &                      MCON_UXVEC(UM,AA)*HMULT_2(AA,UM,N)
	    ENDDO
	    SPAR = ZERO
	    DO AA = 1, NSTREAMS
	      SPAR = SPAR + U_XNEG(UM,AA,N)*SGMULT_DN_DN(AA,UM,N) +
     &                      U_XPOS(UM,AA,N)*SGMULT_DN_UP(AA,UM,N)
	    ENDDO
	    LAYERSOURCE(UM) = SPAR + SHOM
	  ENDDO

	ENDIF

C  Options

C  .. If operating in Ms-mode only, copy multiple scatter term

	IF ( DO_MSMODE_LIDORT ) THEN

	  IF ( SAVE_LAYER_MSST ) THEN
	    DO UM = LOCAL_UM_START, N_USER_STREAMS
	      MSCAT_LAYERSOURCE(UM) = LAYERSOURCE(UM)
	    ENDDO
	  ENDIF

C  .. Full radiance mode, add single scatter part

	ELSE

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    SFOR1 = U_WNEG1(UM,N) * EMULT_DN(UM,N)
	    LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR1
	  ENDDO

	ENDIF

C  Finish

	END

C

	SUBROUTINE PARTLAYER_STERM_UP
     &        ( N, UT, LAYERSOURCE, DO_CONSTANTS_SET )

C  source terms

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of main solution and multiplier variables (input/output)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  arguments

	INTEGER		 N, UT
	DOUBLE PRECISION LAYERSOURCE(MAX_USER_STREAMS)
	LOGICAL		 DO_CONSTANTS_SET

C  local variables

	INTEGER		 AA, UM
	DOUBLE PRECISION SPAR, SHOM, SFOR

C  Save some calculation time

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  DO AA = 1, NSTREAMS
	    LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
	    MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
  	  ENDDO
	ENDDO

C  Get the homogeneous solution multipliers

	CALL PARTLAYER_HMULT_UP ( N, UT )

C  Whole layer source function ( Classical solution )

	IF ( DO_CLASSICAL_SOLUTION ) THEN

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    SPAR =  UT_EMULT_UP(UM,UT) * U_WPOS2(UM,N)
	    SHOM = ZERO
	    DO AA = 1, NSTREAMS
	      SHOM = SHOM + LCON_UXVEC(UM,AA)*UT_HMULT_UP_DN(AA,UM,UT) +
     &                      MCON_UXVEC(UM,AA)*UT_HMULT_UP_UP(AA,UM,UT)
	    ENDDO
	    LAYERSOURCE(UM) = SPAR + SHOM
	  ENDDO

C  Whole layer source function ( Green's function solution )

	ELSE

	  CALL PARTLAYER_GMULT_UP ( N, UT )
  
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    SHOM = ZERO
	    DO AA = 1, NSTREAMS
	      SHOM = SHOM + LCON_UXVEC(UM,AA)*UT_HMULT_UP_DN(AA,UM,UT) +
     &                      MCON_UXVEC(UM,AA)*UT_HMULT_UP_UP(AA,UM,UT)
	    ENDDO
	    SPAR = ZERO
	    DO AA = 1, NSTREAMS
	      SPAR = SPAR + U_XPOS(UM,AA,N)*UT_SGMULT_UP_DN(AA,UM,UT) +
     &                      U_XNEG(UM,AA,N)*UT_SGMULT_UP_UP(AA,UM,UT)
	    ENDDO
	    LAYERSOURCE(UM) = SPAR + SHOM
	  ENDDO

	ENDIF

C  If NOT operating in Ms-mode only, add single scatter part

	IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    SFOR = U_WPOS1(UM,N) * UT_EMULT_UP(UM,UT)
	    LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR
	  ENDDO
	ENDIF
	    
C  Turn off constants_set

	DO_CONSTANTS_SET = .FALSE.

C  Finish

	END

C

	SUBROUTINE PARTLAYER_STERM_DN
     &      ( N, UT, LAYERSOURCE, DO_CONSTANTS_SET )

C  source terms

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model and control variables (input)

	INCLUDE '../include_s/LIDORT_CONTROL.VARS'
	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include files of main solution and multiplier variables (input/output)

	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'
	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  arguments

	INTEGER		 N, UT
	DOUBLE PRECISION LAYERSOURCE(MAX_USER_STREAMS)
	LOGICAL		 DO_CONSTANTS_SET

C  local variables

	INTEGER		 AA, UM
	DOUBLE PRECISION SPAR, SHOM, SFOR

C  Save some calculation time

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  DO AA = 1, NSTREAMS
	    LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
	    MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
  	  ENDDO
	ENDDO

C  Get the homogeneous solution multipliers

	CALL PARTLAYER_HMULT_DN ( N, UT )

C  Whole layer source function ( Classical solution )

	IF ( DO_CLASSICAL_SOLUTION ) THEN

	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    SPAR = UT_EMULT_DN(UM,UT) * U_WNEG2(UM,N)
	    SHOM = ZERO
	    DO AA = 1, NSTREAMS
	      SHOM = SHOM + LCON_UXVEC(UM,AA)*UT_HMULT_DN_DN(AA,UM,UT) +
     &                      MCON_UXVEC(UM,AA)*UT_HMULT_DN_UP(AA,UM,UT)
	    ENDDO
	    LAYERSOURCE(UM) = SPAR + SHOM
	  ENDDO

C  Whole layer source function ( Green's function solution )

	ELSE

	  CALL PARTLAYER_GMULT_DN ( N, UT )
  
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    SHOM = ZERO
	    DO AA = 1, NSTREAMS
	      SHOM = SHOM + LCON_UXVEC(UM,AA)*UT_HMULT_DN_DN(AA,UM,UT) +
     &                      MCON_UXVEC(UM,AA)*UT_HMULT_DN_UP(AA,UM,UT)
	    ENDDO
	    SPAR = ZERO
	    DO AA = 1, NSTREAMS
	      SPAR = SPAR + U_XNEG(UM,AA,N)*UT_SGMULT_DN_DN(AA,UM,UT) +
     &                      U_XPOS(UM,AA,N)*UT_SGMULT_DN_UP(AA,UM,UT)
	    ENDDO
	    LAYERSOURCE(UM) = SPAR + SHOM
	  ENDDO

	ENDIF

C  If NOT operating in Ms-mode only, add single scatter part

	IF ( .NOT. DO_MSMODE_LIDORT ) THEN
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    SFOR = U_WNEG1(UM,N) * UT_EMULT_DN(UM,UT)
	    LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR
	  ENDDO
	ENDIF

C  Turn off constants_set

	DO_CONSTANTS_SET = .FALSE.

C  Finish

	END
