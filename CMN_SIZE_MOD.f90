! $Id: CMN_SIZE,v 4.21 2001/10/29 15:48:34 bmy v421 $

module CMN_SIZE_MOD

      !=======================================================================
      ! CMN_SIZE: size parameters for GEOS-CHEM arrays (bmy, 3/16/01, 10/23/01)
      !           added size parameters for Albedo arrays (GOB, 10/28/09)
      !
      ! NOTES:
      ! (1 ) Now set LLTROP = 20 for GEOS-3 (bmy, 4/12/01)
      ! (2 ) Eliminated obsolete commented-out code (bmy, 4/20/01)
      ! (3 ) Now set MAXFAM = 12 for more P-L families (bmy, 6/28/01)  
      ! (4 ) Comment out {IJL}GCMPAR -- these are obosolete (bmy, 9/24/01)
      ! (5 ) Also set LLPAR = 30 for GEOS-3, will regrid online (bmy, 9/24/01) 
      ! (6 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
      !=======================================================================

      ! C Preprocessor #define statements for conditional compilation
#     include "define.h"

      !=======================================================================
      ! IM   - Longitudinal Dimension I=1,IM
      ! JM   - Latitudinal  Dimension J=1,JM
      ! I0   - Zero-PT start of Longitude grid (W. R. T. absolute global grid)
      ! J0   - Zero-PT start of Latitude  grid (W. R. T. absolute global grid)
      ! IMX  - Maximum Longitude Dimension
      ! JMX  - Maximum Latitude  Dimension
      !=======================================================================
      INTEGER ::      IM, JM, I0, J0, IMX, JMX
      COMMON /CMNSZ1/ IM, JM, I0, J0, IMX, JMX  ! needed for CTM

 ! ----------------
 ! Threadprivate declaration for OpenMP parallelization (dbm)
 !$OMP THREADPRIVATE( /CMNSZ1/ )
 ! ----------------

      !=======================================================================
      ! DISIZE = size (in degrees) of a longitude grid box
      ! DJSIZE = size (in degrees) of a latitude  grid box
      !=======================================================================
#if   defined( GRID4x5  ) 
      REAL*8, PARAMETER :: DISIZE = 5.0d0
      REAL*8, PARAMETER :: DJSIZE = 4.0d0

#elif defined( GRID2x25 )
      REAL*8, PARAMETER :: DISIZE = 2.5d0 
      REAL*8, PARAMETER :: DJSIZE = 2.0d0

#elif defined( GRID1x1 )
      REAL*8, PARAMETER :: DISIZE = 1.0d0 
      REAL*8, PARAMETER :: DJSIZE = 1.0d0

#elif defined( GRIDNEST )
      REAL*8, PARAMETER :: DISIZE = 0.3125d0
      REAL*8, PARAMETER :: DJSIZE = 0.25d0

#endif


      !=======================================================================
      ! DISIZE = size (in degrees) of a longitude grid box for ALBEDOS
      ! DJSIZE = size (in degrees) of a latitude  grid box for ALBEDOS
      !=======================================================================
#if   defined( OMI4x5 ) || defined( GOME4x5 )
      REAL*8, PARAMETER :: DISIZE_ALB = 5.0d0
      REAL*8, PARAMETER :: DJSIZE_ALB = 4.0d0

#elif defined( OMI2x25 ) || defined( GOME2x25 )
      REAL*8, PARAMETER :: DISIZE_ALB = 2.5d0
      REAL*8, PARAMETER :: DJSIZE_ALB = 2.0d0

#elif defined( OMI1x1 ) || defined( GOME1x1 )
      REAL*8, PARAMETER :: DISIZE_ALB = 1.0d0
      REAL*8, PARAMETER :: DJSIZE_ALB = 1.0d0

#elif defined( OMINEST )
      REAL*8, PARAMETER :: DISIZE_ALB = 0.3125d0
      REAL*8, PARAMETER :: DJSIZE_ALB = 0.25d0

#endif




      !=================================================================
      ! GRID PARAMETERS
      !
      ! IGLOB  = global longitude dimension
      ! JGLOB  = global latitude dimension
      ! LGLOB  = max number of sigma levels 
      ! IIPAR  = window longitude dimension
      ! JJPAR  = window latitude dimension
      ! LLPAR  = window vertical dimension
      ! LLTROP = number of tropospheric levels 
      ! PTOP   = model top pressure (mb)
      ! IMIN   = West edge of model
      ! IMAX   = East edge of model
      ! JMIN   = South edge of full model grid
      ! JMAX   = North edge of full model grid
      !
      ! NOTES: 
      ! (1) GEOS-CHEM is usually set up for a global run, thus,
      !     IIPAR=IGLOB, JJPAR=JGLOB, LLPAR=LGLOB (bmy, 4/9/99)
      !
      ! (2) IIPAR, JJPAR, LLPAR may be smaller than or equal to 
      !     IGLOB, JGLOB, LGLOB (bmy, 4/12/99)    
      !
      ! (3) PTOP is now correct (0.01 hPa) for GEOS-2, GEOS-3 grids
      !     (bmy, 10/2/00)
      !
      ! (4) Added IMIN, IMAX, JMIN, JMAX to this file (GOB, 10,28,09)
      !=================================================================

#if   defined( GRID4x5 )
      REAL, PARAMETER :: IMIN = -182.5 
      REAL, PARAMETER :: IMAX = 177.5 
      REAL, PARAMETER :: JMIN = -88.0 
      REAL, PARAMETER :: JMAX = 88.0 
    
#elif defined( GRID2x25 )
      REAL, PARAMETER :: IMIN = -181.25 
      REAL, PARAMETER :: IMAX = 178.75 
      REAL, PARAMETER :: JMIN = -89.0 
      REAL, PARAMETER :: JMAX = 89.0 

#elif defined( GRID1x1 )
      REAL, PARAMETER :: IMIN = -180.5 
      REAL, PARAMETER :: IMAX = 179.5 
      REAL, PARAMETER :: JMIN = -89.5 
      REAL, PARAMETER :: JMAX = 89.5 

#elif defined( GRIDNEST )
      REAL, PARAMETER :: IMIN = 64.84375
      REAL, PARAMETER :: IMAX = 99.84375
      REAL, PARAMETER :: JMIN = 2.125
      REAL, PARAMETER :: JMAX = 37.875

#endif

#if   defined( GEOS_1 ) && defined( GRID4x5 ) 
      INTEGER, PARAMETER :: IGLOB  = 72
      INTEGER, PARAMETER :: JGLOB  = 46
      INTEGER, PARAMETER :: LGLOB  = 20

      INTEGER, PARAMETER :: IIPAR  = IGLOB
      INTEGER, PARAMETER :: JJPAR  = JGLOB
      INTEGER, PARAMETER :: LLPAR  = LGLOB

      INTEGER, PARAMETER :: LLTROP = 16

      REAL*8,  PARAMETER :: PTOP   = 10d0

#elif defined( GEOS_1 ) && defined( GRID2x25 )
      INTEGER, PARAMETER :: IGLOB  = 144
      INTEGER, PARAMETER :: JGLOB  = 91
      INTEGER, PARAMETER :: LGLOB  = 20

      INTEGER, PARAMETER :: IIPAR  = IGLOB
      INTEGER, PARAMETER :: JJPAR  = JGLOB
      INTEGER, PARAMETER :: LLPAR  = LGLOB

      INTEGER, PARAMETER :: LLTROP = 16

      REAL*8,  PARAMETER :: PTOP   = 10d0

#elif defined( GEOS_STRAT ) && defined( GRID4x5 )
      INTEGER, PARAMETER :: IGLOB  = 72
      INTEGER, PARAMETER :: JGLOB  = 46
      INTEGER, PARAMETER :: LGLOB  = 26      

      INTEGER, PARAMETER :: IIPAR  = IGLOB
      INTEGER, PARAMETER :: JJPAR  = JGLOB
      INTEGER, PARAMETER :: LLPAR  = LGLOB

      INTEGER, PARAMETER :: LLTROP = 19

      REAL*8,  PARAMETER :: PTOP   = 0.1d0

#elif defined( GEOS_STRAT ) && defined( GRID2x25 )
      INTEGER, PARAMETER :: IGLOB  = 144
      INTEGER, PARAMETER :: JGLOB  = 91
      INTEGER, PARAMETER :: LGLOB  = 26            

      INTEGER, PARAMETER :: IIPAR  = IGLOB
      INTEGER, PARAMETER :: JJPAR  = JGLOB
      INTEGER, PARAMETER :: LLPAR  = LGLOB

      INTEGER, PARAMETER :: LLTROP = 19

      REAL*8,  PARAMETER :: PTOP   = 0.1d0

#elif defined( GEOS_3 ) && defined( GRID4x5 )
      INTEGER, PARAMETER :: IGLOB  = 72
      INTEGER, PARAMETER :: JGLOB  = 46
      INTEGER, PARAMETER :: LGLOB  = 48            

      INTEGER, PARAMETER :: IIPAR  = IGLOB
      INTEGER, PARAMETER :: JJPAR  = JGLOB
      !-----------------------------------------------------------------
      ! Uncomment this line for GEOS-3 4 x 5 w/ 48 sigma levels 
      ! Note if you do this, then you will need to update the sigma
      !  levels defined earlier in the code
      !INTEGER, PARAMETER :: LLPAR  = LGLOB      
      !-----------------------------------------------------------------
      ! Uncomment this line for GEOS-3 4 x 5 w/ 30 sigma levels 
      INTEGER, PARAMETER :: LLPAR  = 30
      !-----------------------------------------------------------------

      INTEGER, PARAMETER :: LLTROP = 20       

      REAL*8,  PARAMETER :: PTOP   = 0.01d0

#elif defined( GEOS_3 ) && defined( GRID2x25 )
      INTEGER, PARAMETER :: IGLOB  = 144
      INTEGER, PARAMETER :: JGLOB  = 91
      INTEGER, PARAMETER :: LGLOB  = 48            

      INTEGER, PARAMETER :: IIPAR  = IGLOB
      INTEGER, PARAMETER :: JJPAR  = JGLOB
      !-----------------------------------------------------------------
      ! Uncomment this line for GEOS-3 2 x 2.5 with 48 sigma levels
      !INTEGER, PARAMETER :: LLPAR  = LGLOB
      !-----------------------------------------------------------------
      ! Uncomment this line for GEOS-3 2 x 2.5 with 30 sigma levels
      INTEGER, PARAMETER :: LLPAR  = 30
      !-----------------------------------------------------------------

      INTEGER, PARAMETER :: LLTROP = 20       

      REAL*8,  PARAMETER :: PTOP   = 0.01d0

#elif defined( GEOS_3 ) && defined( GRID1x1 )
      INTEGER, PARAMETER :: IGLOB  = 360
      INTEGER, PARAMETER :: JGLOB  = 181
      INTEGER, PARAMETER :: LGLOB  = 48           

      INTEGER, PARAMETER :: IIPAR  = IGLOB
      INTEGER, PARAMETER :: JJPAR  = JGLOB
      !-----------------------------------------------------------------
      ! Uncomment this line for GEOS-3 1 x 1 with 48 sigma levels
      !INTEGER, PARAMETER :: LLPAR  = LGLOB
      !-----------------------------------------------------------------
      ! Uncomment this line for GEOS-3 1 x 1 with 30 sigma levels
      INTEGER, PARAMETER :: LLPAR  = 30
      !-----------------------------------------------------------------

      INTEGER, PARAMETER :: LLTROP = 20       

      REAL*8,  PARAMETER :: PTOP   = 0.01d0

#elif defined( GEOS_4 ) && defined( GRID4x5 )

      !--------------------
      ! GEOS-4 4 x 5
      !--------------------
      INTEGER, PARAMETER :: IGLOB  = 72
      INTEGER, PARAMETER :: JGLOB  = 46
      INTEGER, PARAMETER :: LGLOB  = 55

      INTEGER, PARAMETER :: IIPAR  = IGLOB
      INTEGER, PARAMETER :: JJPAR  = JGLOB

#if   defined( GRIDREDUCED )
      INTEGER, PARAMETER :: LLPAR  = 30     ! 30 levels
#else
      INTEGER, PARAMETER :: LLPAR  = LGLOB  ! 55 levels
#endif

      INTEGER, PARAMETER :: LLTROP = 17

      REAL*8,  PARAMETER :: PTOP   = 0.01d0

#elif defined( GEOS_4 ) && defined( GRID2x25 )

      !--------------------
      ! GEOS-4 2 x 2.5
      !--------------------
      INTEGER, PARAMETER :: IGLOB  = 144
      INTEGER, PARAMETER :: JGLOB  = 91
      INTEGER, PARAMETER :: LGLOB  = 55

      INTEGER, PARAMETER :: IIPAR  = IGLOB
      INTEGER, PARAMETER :: JJPAR  = JGLOB

#if   defined( GRIDREDUCED )
      INTEGER, PARAMETER :: LLPAR  = 30     ! 30 levels
#else
      INTEGER, PARAMETER :: LLPAR  = LGLOB  ! 55 levels
#endif

      INTEGER, PARAMETER :: LLTROP = 17

      REAL*8,  PARAMETER :: PTOP   = 0.01d0
 
#elif defined( GEOS_5 ) && defined( GRID2x25 )

      !--------------------
      ! GEOS-5 2 x 2.5
      !--------------------
      INTEGER, PARAMETER :: IGLOB      = 144
      INTEGER, PARAMETER :: JGLOB      = 91
      INTEGER, PARAMETER :: LGLOB      = 72

      INTEGER, PARAMETER :: IIPAR      = IGLOB
      INTEGER, PARAMETER :: JJPAR      = JGLOB

#elif defined( GEOS_5 ) && defined( GRIDNEST )

      !-------------------
      ! GEOS-5, INDIA 0.25x0.3125
      !-------------------
      INTEGER, PARAMETER :: IGLOB      = 112
      INTEGER, PARAMETER :: JGLOB      = 144
      INTEGER, PARAMETER :: LGLOB      = 72

      INTEGER, PARAMETER :: IIPAR      = IGLOB
      INTEGER, PARAMETER :: JJPAR      = JGLOB

!---------------------------------
! Prior to 2/2/07:
!#if   defined( GRID30LEV )
!---------------------------------
#if   defined( GRIDREDUCED )
      INTEGER, PARAMETER :: LLPAR      = 47     ! 47 levels
#else 
      INTEGER, PARAMETER :: LLPAR      = LGLOB  ! 72 levels
#endif

      ! Set this so that we read data in correctly for offline sims
      INTEGER, PARAMETER :: LLTROP_FIX = 38  
      INTEGER, PARAMETER :: LLTROP     = 38  

      REAL*8,  PARAMETER :: PTOP       = 0.01d0               

#endif

      REAL*8 MSIGMA(LLPAR) ! Sigma coordinates from model (centre)
      REAL*8 MSIGMAE(LLPAR+1) ! Sigma edges
      REAL*8 ETA_A(LLPAR+1) ! A coordinates from model (edges)
      REAL*8 ETA_B(LLPAR+1) ! B coordinates from model (edges)

#if   defined( GEOS_1  ) 
      DATA MSIGMA							&
       / .993936d0, .971301d0, .929925d0, .874137d0, .807833d0, 	&
         .734480d0, .657114d0, .578390d0, .500500d0, .424750d0,		&
         .352000d0, .283750d0, .222750d0, .172150d0, .132200d0, 	&
         .100050d0, .073000d0, .049750d0, .029000d0, .009500d0 /	

      DATA MSIGMAE							&
       / 1.0d0, .987871d0, .95473d0, .90512d0, .843153d0, .772512d0,	&
         .696448d0, .617779d0, .539000d0, .462000d0, .387500d0, 	&
         .316500d0, .251000d0, .194500d0, .149800d0, .114600d0,  	&
         .085500d0, .060500d0, .039000d0, .019000d0, .00000d0 /

#elif defined( GEOS_STRAT )
      DATA MSIGMA							&
        /     .993935d0, .971300d0, .929925d0, .875060d0, .812500d0,	&
              .745000d0, .674500d0, .604500d0, .536500d0, .471500d0,	&
              .410000d0, .352500d0, .301500d0, .257977d0, .220273d0,	&
              .187044d0, .157881d0, .132807d0, .111722d0, .094035d0,	&
              .079233d0, .066873d0, .056574d0, .044794d0, .028825d0,	&
              .009979d0 /

      DATA MSIGMAE							&
        / 1.0d0, .987871d0, .954730d0, .905120d0, .845000d0, .78d0,	&
              .710000d0, .639000d0, .570000d0, .503000d0, .440000d0,	&
              .380000d0, .325000d0, .278000d0, .237954d0, .202593d0,	&
              .171495d0, .144267d0, .121347d0, .102098d0, .085972d0,	&
              .072493d0, .061252d0, .051896d0, .037692d0, .019958d0,	&
              .000000d0 /

#elif defined( GEOS_3 )
      DATA MSIGMA							&
         / .998548d0, .994148d0, .986350d0, .974300d0, .956950d0,	&
           .933150d0, .901750d0, .861500d0, .811000d0, .750600d0,	&
           .682900d0, .610850d0, .537050d0, .463900d0, .393650d0, 	&
           .328275d0, .269500d0, .218295d0, .174820d0, .138840d0,	&
           .109790d0, .086690d0, .062045d0, .038605d0, .023990d0,	&
           .012710d0, .004785d0, .0016475d0, .00046d0, 			&
           7.75000d-05 /

      DATA MSIGMAE							&
         / 1.0d0, .997095d0, .9912d0, .9815d0, .9671d0, .9468d0,	&
           .919500d0, .884000d0, .839000d0, .783000d0, .718200d0,	&
           .647600d0, .574100d0, .500000d0, .427800d0, .359500d0,	&
           .297050d0, .241950d0, .194640d0, .155000d0, .122680d0,	&
           .096900d0, .076480d0, .047610d0, .029600d0, .018380d0,	&
           .000704d0, .002530d0, .000765d0, .000155d0, 			&
           .000000d0 /

#elif defined( GEOS_4 )
      DATA ETA_A							&
         / 0.000000d0, 0.000000d0, 12.704939d0, 35.465965d0, 		&
           66.098427d0, 101.671654d0, 138.744400d0, 173.403183d0, 	&
           198.737839d0, 215.417526d0, 223.884689d0, 224.362869d0, 	&
           216.864929d0, 201.192093d0, 176.929993d0, 150.393005d0, 	&
           127.837006d0, 108.663429d0, 92.365662d0, 78.512299d0, 	&
           56.387939d0, 40.175419d0, 28.367815d0, 19.791553d0,	 	&
           9.292943d0, 4.076567d0, 1.650792d0, 0.616779d0, 		&
           0.211349d0, 0.066000d0, 0.010000d0 /

      DATA ETA_B 							&
         / 1.000000d0, 0.985110d0, 0.943290d0, 0.867830d0, 		&
           0.764920d0, 0.642710d0, 0.510460d0, 0.378440d0, 		&
           0.270330d0, 0.183300d0, 0.115030d0, 0.063720d0, 		&
           0.028010d0, 0.006960d0, 0.000000d0, 0.000000d0, 		&
           0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 		&
           0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 		&
           0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 		&
           0.000000d0, 0.000000d0, 0.000000d0 /
#elif defined( GEOS_5 )
      DATA ETA_A							&
         / 0.000000d+00, 4.804826d-02, 6.593752d+00, 1.313480d+01,	&
             1.961311d+01, 2.609201d+01, 3.257081d+01, 3.898201d+01,	&
             4.533901d+01, 5.169611d+01, 5.805321d+01, 6.436264d+01,	&
             7.062198d+01, 7.883422d+01, 8.909992d+01, 9.936521d+01,	&
             1.091817d+02, 1.189586d+02, 1.286959d+02, 1.429100d+02,	&
             1.562600d+02, 1.696090d+02, 1.816190d+02, 1.930970d+02,	&
             2.032590d+02, 2.121500d+02, 2.187760d+02, 2.238980d+02,	&
             2.243630d+02, 2.168650d+02, 2.011920d+02, 1.769300d+02,	&
             1.503930d+02, 1.278370d+02, 1.086630d+02, 9.236572d+01,	&
             7.851231d+01, 5.638791d+01, 4.017541d+01, 2.836781d+01,	&
             1.979160d+01, 9.292942d+00, 4.076571d+00, 1.650790d+00,	&
             6.167791d-01, 2.113490d-01, 6.600001d-02, 1.000000d-02 /	

      DATA ETA_B 							&
         / 1.000000d+00, 9.849520d-01, 9.634060d-01, 9.418650d-01,	&
             9.203870d-01, 8.989080d-01, 8.774290d-01, 8.560180d-01,	&
             8.346609d-01, 8.133039d-01, 7.919469d-01, 7.706375d-01,	&
             7.493782d-01, 7.211660d-01, 6.858999d-01, 6.506349d-01,	&
             6.158184d-01, 5.810415d-01, 5.463042d-01, 4.945902d-01,	&
             4.437402d-01, 3.928911d-01, 3.433811d-01, 2.944031d-01,	&
             2.467411d-01, 2.003501d-01, 1.562241d-01, 1.136021d-01,	&
             6.372006d-02, 2.801004d-02, 6.960025d-03, 8.175413d-09,	&
             0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,	&
             0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,	&
             0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,	&
             0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00 /	
#endif

!-----------------------------------------------------------------
!  !Define grid for surface reflectivities
!-----------------------------------------------------------------
#if defined ( GOME1x1 ) || defined ( OMI1x1 )
      INTEGER, PARAMETER :: IIPAR_ALB=360  !(GOB)
      INTEGER, PARAMETER :: JJPAR_ALB=180
#else
      INTEGER, PARAMETER :: IIPAR_ALB=IIPAR
      INTEGER, PARAMETER :: JJPAR_ALB=JJPAR
#endif

#if   defined( OMI4x5 ) || defined ( GOME4x5 )
      REAL, PARAMETER :: IMIN_ALB = -182.5 
      REAL, PARAMETER :: IMAX_ALB = 177.5 
      REAL, PARAMETER :: JMIN_ALB = -88.0
      REAL, PARAMETER :: JMAX_ALB = 88.0 

#elif defined( OMI2x25 ) || defined( GOME2x25 )
      REAL, PARAMETER :: IMIN_ALB = -181.25 
      REAL, PARAMETER :: IMAX_ALB = 178.75 
      REAL, PARAMETER :: JMIN_ALB = -89.0 
      REAL, PARAMETER :: JMAX_ALB = 89.0 

#elif defined( OMI1x1 ) || defined( GOME1x1)
      REAL, PARAMETER :: IMIN_ALB = -180 
      REAL, PARAMETER :: IMAX_ALB = 180 
      REAL, PARAMETER :: JMIN_ALB = -90 
      REAL, PARAMETER :: JMAX_ALB = 90 

#elif defined( OMINEST )
      REAL, PARAMETER :: IMIN_ALB = 64.84375
      REAL, PARAMETER :: IMAX_ALB = 100.15625
      REAL, PARAMETER :: JMIN_ALB = 1.875
      REAL, PARAMETER :: JMAX_ALB = 38.125
#endif






      !=================================================================
      ! TRACER & EMISSION SPECIES PARAMETERS
      !
      ! NNPAR   = max number of tracers
      ! NEMPARA = max number of anthropogenic emission species
      ! NEMPARB = max number of biogenic      emission species
      !
      !=================================================================
#if   defined( SMALLCHEM )
      INTEGER, PARAMETER :: NNPAR   = 6
      INTEGER, PARAMETER :: NEMPARA = 10
      INTEGER, PARAMETER :: NEMPARB = 1

#elif defined( FULLCHEM  )
      INTEGER, PARAMETER :: NNPAR   = 24
      INTEGER, PARAMETER :: NEMPARA = 10
      INTEGER, PARAMETER :: NEMPARB = 1

#elif defined( LGEOSCO ) && !defined( FULLCHEM ) && !defined( SMALLCHEM )
      INTEGER, PARAMETER :: NNPAR   = 10
      INTEGER, PARAMETER :: NEMPARA = 10
      INTEGER, PARAMETER :: NEMPARB = 1
      
#endif

      !=================================================================
      ! OTHER PARAMETERS 
      !=================================================================

      ! NVEGTYPE - Maximum number of surface types: 74 olson
      ! NTYPE    - Maximum number of veg types in a 4x5 box
      ! NPOLY    - Number of coefficients for polynomial fits
      INTEGER, PARAMETER :: NVEGTYPE = 74
      INTEGER, PARAMETER :: NTYPE    = 15
      INTEGER, PARAMETER :: NPOLY    = 20

      ! NAIR     - Maximum number (km) for aircraft NOx emissions   
      ! LAIREMS  - Maximum number of layers for aircraft emissions AIREMIS
      INTEGER, PARAMETER :: NAIR    = 20
      INTEGER, PARAMETER :: LAIREMS = LLTROP

      ! NDUST -- number of mineral dust categories rvm, 08/17/00
      INTEGER, PARAMETER :: NDUST = 7

      ! NAER -- number of other aerosol categories rvm, 11/14/01
      INTEGER, PARAMETER :: NAER = 5

      ! NRH -- number of relative humidity bins rvm, 12/14/01
      INTEGER, PARAMETER :: NRH = 5

      ! NNSTA = max number of time series stations (in inptr.ctm)
      INTEGER, PARAMETER :: NNSTA = 800

      ! MAXIJ - Maximum number of 1st level grid boxes
      INTEGER, PARAMETER :: MAXIJ = IIPAR * JJPAR

      ! MAXDEP - Maximum number of depositing species for drydep
      INTEGER, PARAMETER :: MAXDEP = 25

      ! Now define NUMDEP in DRYDEP_SETUP.F, but put it into
      ! a common block here (bmy, 4/9/99)
      INTEGER            :: NUMDEP
      COMMON /BMYDEP1/      NUMDEP
      
 ! ----------------
 ! Threadprivate declaration for OpenMP parallelization (dbm)
 !$OMP THREADPRIVATE( /BMYDEP1/ )
 ! ----------------

      ! LLCONVM - Max number of layers for convection
      INTEGER, PARAMETER :: LLCONVM = LLPAR - 1

      ! NOXLEVELS = Number of levels of anthro NOx emission 
      !             (e.g. surface and 100m)
      ! NOXEXTENT = Highest sigma level that receives anthro NOx emission 
      INTEGER, PARAMETER :: NOXLEVELS = 2
      INTEGER, PARAMETER :: NOXEXTENT = 2 

      ! MAXFAM -- Max number of families for prod and loss output
      INTEGER, PARAMETER :: MAXFAM = 17

      ! TAU values at the first day of each month for 1985
      ! needed to index binary punch files
      REAL*8               :: TAU_VAL(12) = (/   0d0,  744d0, 1416d0,	&
                                             2160d0, 2880d0, 3624d0,	&
                                             4344d0, 5088d0, 5832d0,	&
                                             6552d0, 7296d0, 8016d0/)

!   Mid-latitude atmosphere (Z, P, T) supplied by K. Chance
!   US standard atmosphere for 45N, July in Geometric height

        INTEGER          NUSAML
        PARAMETER        ( NUSAML = 105 )
        REAL*8 USAML_ZZZ(NUSAML)
        REAL*8 USAML_TTT(NUSAML)
        REAL*8 USAML_PPP(NUSAML)
        REAL*8 USAML_LNP(NUSAML)

        COMMON / USAML / USAML_ZZZ, USAML_TTT, USAML_PPP, USAML_LNP

        SAVE   / USAML /
 ! ----------------
 ! Threadprivate declaration for OpenMP parallelization (dbm)
 !$OMP THREADPRIVATE( /USAML/ )
 ! ----------------

!   Aerosol parameters
        INTEGER         NWL
        PARAMETER       ( NWL = 7 )
        REAL*8 WL(NWL)          ! Wavelength
        REAL*8 QEXT(NWL,NAER+1)	!Extinction efficiency
        REAL*8 RAA(NWL,NAER+1)	!Aerosol effective radius
        REAL*8 ASSA(NWL,NAER+1)	!Aerosol single scattering albedo
        REAL*8 PHFCN(0:7,NWL,NAER+1)	!Aerosol phase function

        COMMON / ARSL_PARAMS /  WL, QEXT, RAA, ASSA, PHFCN
        SAVE   / ARSL_PARAMS /
 ! ----------------
 ! Threadprivate declaration for OpenMP parallelization (dbm)
 !$OMP THREADPRIVATE( /ARSL_PARAMS/ )
 !----------------

END MODULE CMN_SIZE_MOD
