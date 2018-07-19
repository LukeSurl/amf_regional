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

	SUBROUTINE WHOLELAYER_HMULT_UP ( N )

C  Whole layer INTEGRATED homogeneous solution MULTIPLIERS (Upwelling)

C  Include files
C  =============

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of setup stuff (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include files of Multiplier variables (output to this module)

	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index

	INTEGER		 N

C  Local variables
C  ---------------

	INTEGER		 UM, AA
        DOUBLE PRECISION UDEL, ZDEL, ZUDEL, SM, THETA_2, THETA_1

C  set the existence flag

	HMULT_EXIST(N) = .TRUE.

C  Get the multipliers

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  UDEL = T_DELT_USERM(N,UM)
	  SM = USER_SECANTS(UM)
	  DO AA = 1, NSTREAMS
	    ZDEL  = T_DELT_EIGEN(AA,N)
	    ZUDEL = ZDEL * UDEL
	    THETA_2 = ONE - ZUDEL
	    THETA_1 = ZDEL - UDEL
	    HMULT_1(AA,UM,N) = SM * THETA_1 * ZETA_M(AA,UM,N)
	    HMULT_2(AA,UM,N) = SM * THETA_2 * ZETA_P(AA,UM,N)
	  ENDDO
	ENDDO

C  Finish

	END

C

	SUBROUTINE WHOLELAYER_HMULT_DN ( N )

C  INTEGRATED homogeneous solution MULTIPLIERS (Downwelling)

C  Include files
C  =============

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of setup stuff (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include files of Multiplier variables (output to this module)

	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index

	INTEGER		 N

C  Local variables
C  ---------------

	INTEGER		 UM,  AA
        DOUBLE PRECISION UDEL, ZDEL, ZUDEL, SM, THETA_2, THETA_1

C  If this already exists, return (saves re-computation)

	IF ( HMULT_EXIST(N) ) RETURN

C  Get the multipliers

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  UDEL = T_DELT_USERM(N,UM)
	  SM = USER_SECANTS(UM)
	  DO AA = 1, NSTREAMS
	    ZDEL  = T_DELT_EIGEN(AA,N)
	    ZUDEL = ZDEL * UDEL
	    THETA_2 = ONE - ZUDEL
	    THETA_1 = ZDEL - UDEL
	    HMULT_1(AA,UM,N) = SM * THETA_1 * ZETA_M(AA,UM,N)
	    HMULT_2(AA,UM,N) = SM * THETA_2 * ZETA_P(AA,UM,N)
	  ENDDO
	ENDDO

C  Finish

	END

C

	SUBROUTINE PARTLAYER_HMULT_UP ( N, UT )

C  Partial layer INTEGRATED homogeneous solution MULTIPLIERS (Upwelling)

C  Include files
C  =============

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of setup stuff (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include files of Multiplier variables (output to this module)

	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer indices

	INTEGER		 N, UT

C  Local variables
C  ---------------

	INTEGER		 UM, AA
        DOUBLE PRECISION ZDEL, UX_UP, ZX_UP, ZX_DN
	DOUBLE PRECISION SM, THETA_DN, THETA_UP

C  Partial layer multipliers

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  UX_UP= T_UTUP_USERM(UT,UM)
	  SM = USER_SECANTS(UM)
	  DO AA = 1, NSTREAMS
	    ZDEL  = T_DELT_EIGEN(AA,N)
	    ZX_UP  = T_UTUP_EIGEN(AA,UT)
	    ZX_DN  = T_UTDN_EIGEN(AA,UT)
	    THETA_DN = ZX_DN - ZDEL*UX_UP
	    THETA_UP = ZX_UP - UX_UP
	    UT_HMULT_UP_DN(AA,UM,UT) = SM * THETA_DN * ZETA_P(AA,UM,N)
	    UT_HMULT_UP_UP(AA,UM,UT) = SM * THETA_UP * ZETA_M(AA,UM,N)
	  ENDDO
	ENDDO

C  Finish

	END

C

	SUBROUTINE PARTLAYER_HMULT_DN ( N, UT )

C  Partial layer INTEGRATED homogeneous solution MULTIPLIERS (Downwelling)

C  Include files
C  =============

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of setup stuff (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'

C  include files of Multiplier variables (output to this module)

	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer indices

	INTEGER		 N, UT

C  Local variables
C  ---------------

	INTEGER		 UM, AA
        DOUBLE PRECISION ZDEL, UX_DN, ZX_UP, ZX_DN
	DOUBLE PRECISION SM, THETA_DN, THETA_UP

C  Partial layer multipliers

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  UX_DN= T_UTDN_USERM(UT,UM)
	  SM = USER_SECANTS(UM)
	  DO AA = 1, NSTREAMS
	    ZDEL  = T_DELT_EIGEN(AA,N)
	    ZX_UP = T_UTUP_EIGEN(AA,UT)
	    ZX_DN = T_UTDN_EIGEN(AA,UT)
	    THETA_DN = ZX_DN - UX_DN
	    THETA_UP = ZX_UP - ZDEL*UX_DN
	    UT_HMULT_DN_DN(AA,UM,UT) = SM * THETA_DN * ZETA_M(AA,UM,N)
	    UT_HMULT_DN_UP(AA,UM,UT) = SM * THETA_UP * ZETA_P(AA,UM,N)
	  ENDDO
	ENDDO

C  Finish

	END

C

	SUBROUTINE QUAD_GFUNCMULT ( N, UT, DO_GMULT )

C  Include files
C  =============

C  input
C  -----

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of Legendre stuff (input to this module)

	INCLUDE '../include_s/LIDORT_LEGENDRE.VARS'

C  include file of setup and solution variables (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  output
C  ------

C  include files of Multiplier variables (output to this module)

	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

	INTEGER		 N, UT
	LOGICAL		 DO_GMULT

C  Local variables
C  ---------------

	INTEGER		 AA
	DOUBLE PRECISION SD, SU, WDEL 
	DOUBLE PRECISION ZX_DN, ZX_UP, ZW, WX, CONST

C  Return if already done

	IF ( .NOT. DO_GMULT ) RETURN

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO AA = 1, NSTREAMS
	    UT_GMULT_DN(AA,UT) = ZERO
	    UT_GMULT_UP(AA,UT) = ZERO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################
	  
C  Layer constant terms

	WX    = T_UTDN_MUBAR(UT)
	WDEL  = T_DELT_MUBAR(N)
	CONST = INITIAL_TRANS(N)

C  Tau integration without coefficients (average secant approximation)

	DO AA = 1, NSTREAMS
	  ZX_DN = T_UTDN_EIGEN(AA,UT)
	  ZX_UP = T_UTUP_EIGEN(AA,UT)
	  ZW    = WDEL * ZX_UP
	  SD =  ( ZX_DN - WX ) * GAMMA_M(AA,N)
	  SU =  ( WX    - ZW ) * GAMMA_P(AA,N)
	  UT_GMULT_DN(AA,UT) = SD * ATERM_SAVE(AA,N) * CONST
	  UT_GMULT_UP(AA,UT) = SU * BTERM_SAVE(AA,N) * CONST
	ENDDO

C  Finish

	END

C

	SUBROUTINE WHOLELAYER_GMULT_UP ( N )

C  INTEGRATED GREEN FUNCTION MULTIPLIER (Upwelling Whole layer)

C  Include files
C  =============

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of setup and solution stuff (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include files of Multiplier variables (output to this module)

	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index

	INTEGER		 N

C  Local variables
C  ---------------

	INTEGER		 UM, AA
	DOUBLE PRECISION WDEL, CONS, CONSWDEL, SD, SU, TD, TU

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO AA = 1, NSTREAMS
	      SGMULT_UP_DN(AA,UM,N) = ZERO
	      SGMULT_UP_UP(AA,UM,N) = ZERO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  Layer quantities

	WDEL = T_DELT_MUBAR(N)
	CONS = INITIAL_TRANS(N)
	CONSWDEL = - CONS * WDEL

C  Get the multipliers

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  DO AA = 1, NSTREAMS
	    TD = CONS     * HMULT_2(AA,UM,N)
	    TU = CONSWDEL * HMULT_1(AA,UM,N)
	    SD  = GAMMA_M(AA,N) * ( TD - EMULT_UP(UM,N) )
	    SU  = GAMMA_P(AA,N) * ( TU + EMULT_UP(UM,N) )
	    SGMULT_UP_DN(AA,UM,N) = SD * ATERM_SAVE(AA,N)
	    SGMULT_UP_UP(AA,UM,N) = SU * BTERM_SAVE(AA,N)
	  ENDDO
	ENDDO

C  Finish

	END

C

	SUBROUTINE WHOLELAYER_GMULT_DN ( N )

C  INTEGRATED GREEN FUNCTION MULTIPLIER (Downwelling Whole layer)

C  Include files
C  =============

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of setup and solution stuff (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include files of Multiplier variables (output to this module)

	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index

	INTEGER		 N

C  Local variables
C  ---------------

	INTEGER		 UM, AA
	DOUBLE PRECISION WDEL, CONS, CONSWDEL, SD, SU, TD, TU

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO AA = 1, NSTREAMS
	      SGMULT_DN_DN(AA,UM,N) = ZERO
	      SGMULT_DN_UP(AA,UM,N) = ZERO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  Layer quantities

	WDEL = T_DELT_MUBAR(N)
	CONS = INITIAL_TRANS(N)
	CONSWDEL = - CONS * WDEL

C  Get the basic multipliers

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  DO AA = 1, NSTREAMS
	    TD = CONS     * HMULT_1(AA,UM,N)
	    TU = CONSWDEL * HMULT_2(AA,UM,N)
	    SD  = GAMMA_M(AA,N) * ( TD - EMULT_DN(UM,N) )
	    SU  = GAMMA_P(AA,N) * ( TU + EMULT_DN(UM,N) )
	    SGMULT_DN_DN(AA,UM,N) = SD * ATERM_SAVE(AA,N)
	    SGMULT_DN_UP(AA,UM,N) = SU * BTERM_SAVE(AA,N)
	  ENDDO
	ENDDO

C  Finish

	END

C

	SUBROUTINE PARTLAYER_GMULT_UP ( N, UT )

C  INTEGRATED GREEN FUNCTION MULTIPLIER (Upwelling partial layer)

C  Include files
C  =============

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of setup and solution stuff (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include files of Multiplier variables (output to this module)

	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index, offgrid index

	INTEGER		 N, UT

C  Local variables
C  ---------------

	INTEGER		 UM, AA
	DOUBLE PRECISION WDEL, CONS, CONSWDEL, SD, SU, TD, TU

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO AA = 1, NSTREAMS
	      UT_SGMULT_UP_DN(AA,UM,UT) = ZERO
	      UT_SGMULT_UP_UP(AA,UM,UT) = ZERO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  Layer quantities

	WDEL = T_DELT_MUBAR(N)
	CONS = INITIAL_TRANS(N)
	CONSWDEL = - CONS * WDEL

C  Get the multipliers

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  DO AA = 1, NSTREAMS
	    TD =   CONS     * UT_HMULT_UP_DN(AA,UM,UT)
	    TU =   CONSWDEL * UT_HMULT_UP_UP(AA,UM,UT)
	    SD  = GAMMA_M(AA,N) * ( TD - UT_EMULT_UP(UM,UT) )
	    SU  = GAMMA_P(AA,N) * ( TU + UT_EMULT_UP(UM,UT) )
	    UT_SGMULT_UP_DN(AA,UM,UT) = SD * ATERM_SAVE(AA,N)
	    UT_SGMULT_UP_UP(AA,UM,UT) = SU * BTERM_SAVE(AA,N)
	  ENDDO
 	ENDDO

C  Finish

	END

C

	SUBROUTINE PARTLAYER_GMULT_DN ( N, UT )

C  INTEGRATED GREEN FUNCTION MULTIPLIER (Downwelling partial layer)

C  Include files
C  =============

C  include file of dimensions and numbers

	INCLUDE '../include_s/LIDORT.PARS'

C  Include file of model variables

	INCLUDE '../include_s/LIDORT_MODEL.VARS'

C  include file of setup and solution stuff (input to this module)

	INCLUDE '../include_s/LIDORT_SETUPS.VARS'
	INCLUDE '../include_s/LIDORT_SOLUTION.VARS'

C  include files of Multiplier variables (output to this module)

	INCLUDE '../include_s/LIDORT_MULTIPLIERS.VARS'

C  subroutine arguments
C  --------------------

C  Given layer index, offgrid index

	INTEGER		 N, UT

C  Local variables
C  ---------------

	INTEGER		 UM, AA
	DOUBLE PRECISION WDEL, CONS, CONSWDEL, SD, SU, TD, TU

C   ### NEW CODE ##################################################
C  No particular solution beyond the cutoff layer.
C    [ Zero the multiplier values and exit )
	IF ( N .GT. LAYER_PIS_CUTOFF ) THEN
	  DO UM = LOCAL_UM_START, N_USER_STREAMS
	    DO AA = 1, NSTREAMS
	      UT_SGMULT_DN_DN(AA,UM,UT) = ZERO
	      UT_SGMULT_DN_UP(AA,UM,UT) = ZERO
	    ENDDO
	  ENDDO
	  RETURN
	ENDIF
C   ### NEW CODE ##################################################

C  Layer quantities

	WDEL = T_DELT_MUBAR(N)
	CONS = INITIAL_TRANS(N)
	CONSWDEL = - CONS * WDEL

C  Get the multipliers

	DO UM = LOCAL_UM_START, N_USER_STREAMS
	  DO AA = 1, NSTREAMS
	    TD =   CONS     * UT_HMULT_DN_DN(AA,UM,UT)
	    TU =   CONSWDEL * UT_HMULT_DN_UP(AA,UM,UT)
	    SD  = GAMMA_M(AA,N) * ( TD - UT_EMULT_DN(UM,UT) )
	    SU  = GAMMA_P(AA,N) * ( TU + UT_EMULT_DN(UM,UT) )
	    UT_SGMULT_DN_DN(AA,UM,UT) = SD * ATERM_SAVE(AA,N)
	    UT_SGMULT_DN_UP(AA,UM,UT) = SU * BTERM_SAVE(AA,N)
	  ENDDO
 	ENDDO

C  Finish

	END

