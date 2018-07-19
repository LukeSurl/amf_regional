! $Id: define.h,v 4.21 2001/10/29 15:48:34 bmy v421 $
!
!==============================================================================
! Undefine all "switches" so that they cannot be accidentally reset  
!==============================================================================
#undef GEOS_1
#undef GEOS_STRAT
#undef GEOS_2
#undef GEOS_3
#undef GEOS_4
#undef GEOS_5
#undef GRIDREDUCED
#undef GRID2x25  
#undef GRID4x5
#undef GRID1x1   
#undef GRIDNEST
#undef FULLCHEM  
#undef SMALLCHEM 
#undef LGEOSCO
#undef LFASTJ
#undef LSLOWJ
#undef GOME
#undef SCIA
#undef NEW_GOME
#undef OMI
#undef USELUT
#undef NOUSELUT
#undef OMI1x01
#undef GOME1x1
#undef OMI2x25
#undef OMI4x5
#undef OMINEST
#undef GOME2x25
#undef GOME4x5


!==============================================================================
! Define the necessary "switches" for the AMF code. 
! Give each switch its own name as a value, since this will prevent
! the C-preprocessor overwriting the name everywhere in the code.
!==============================================================================
!-----------------------------------------------------
! Pick the satellite
!-----------------------------------------------------

!#define GOME 'GOME'
!#define SCIA 'SCIA'
!#define NEW_GOME 'NEW_GOME'
#define OMI 'OMI'

!-----------------------------------------------------
! Pick use of look up table
!-----------------------------------------------------

#define NOUSELUT   'NOUSELUT'
!#define USELUT 'USELUT'

!-----------------------------------------------------
! GEOS-Chem switches
!-----------------------------------------------------

!#define GEOS_1      'GEOS_1'       
!#define GEOS_STRAT  'GEOS_STRAT'
!#define GEOS_2      'GEOS_2'
!#define GEOS_3      'GEOS_3'
!#define GEOS_4      'GEOS_4'
#define GEOS_5      'GEOS_5'

!#define GRID1x1     'GRID1x1'
!#define GRID2x25    'GRID2x25'
!#define GRID4x5     'GRID4x5'
#define GRIDNEST    'GRIDNEST'
#define GRIDREDUCED   'GRIDREDUCED'

!------------------------------------------------------------------
! Choose Albedos for NO2 (GOB)
!  Default dataset is GOME 5 year minimum surface albedos (Koelemeijer 2003)
!  Also available is OMI 3 year surface LER (Kleipool 2008)
!
!  Default resolution is the above GRID resolution
!  But others can be chosen
!------------------------------------------------------------------
!#define OMI1x1     'OMI1x1'
!#define OMI2x25    'OMI2x25'
!#define OMI4x5     'OMI4x5'
#define OMINEST    'OMINEST'

!#define GOME1x1    'GOME1x1'
!#define GOME2x25   'GOME2x25'
!#define GOME4x5    'GOME4x5'


!-----------------------------------------------------
! Pick the machine that you are running on
!-----------------------------------------------------

#define ALTIX 'ALTIX'
!#define PC    'PC'
!#define SGI32 'SGI32'
!#define SGI64 'SGI64'

!-----------------------------------------------------
! Define integer size for certain HDF-EOS5 functions
!
! On Altix, you need to pass 64-bit integers 
! (aka INTEGER*8) to some HDF routines
!-----------------------------------------------------

#if defined( ALTIX )

! 64-bit Altix
#define NEED_INT_64 'NEED_INT_64'

#elif defined( PC )

! 32-bit linux box
#define NEED_INT_32 'NEED_INT_32'

#elif defined( SGI32 )

! 32-bit SGI 
#define NEED_INT_32 'NEED_INT_32'

#elif defined( SGI64 )

! 64-bit SGI
#define NEED_INT_32 'NEED_INT_32'

#endif



#if !defined( GEOS_1 ) && !defined( GEOS_STRAT ) && !defined( GEOS_2 ) && !defined( GEOS_3 ) && !defined (GEOS_4) &&  !defined (GEOS_5)
#error "ERROR: GEOS_1, GEOS_STRAT, GEOS_2, and GEOS_3"
#error "are ALL undefined in header file define.h"
#endif

#if !defined( GRID2x25 ) && !defined( GRID4x5 ) && !defined( GRID1x1 ) && !defined( GRIDNEST )
#error "ERROR: GRID2x25, GRID4x5, and GRID1x1, and GRIDNEST"
#error "are ALL undefined in header file define.h"
#endif
