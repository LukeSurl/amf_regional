      SUBROUTINE MKATM(SPIN,CLDOD,CTP,GVMR,LAMDA,AERPROF,
     &                 ZGRID,NGRID,DELTA,
     &                 PRESS,TEMPERATURE,OTHERGAS_VMR,NDENS,
     &                 AERSCA,AEREXT,PHASEFCN)

C----------------------------------------------------------------
C Constructs the atmospheres used for the LIDORT runs.
C
C----------------------------------------------------------------

C  include files
      USE CMN_SIZE_MOD        ! Size parameters


      IMPLICIT NONE

C  input values

      REAL*8 SPIN                 ! Input surface pressure
      REAL*8 CLDOD                ! Input cloud optical depth
      REAL*8 CTP                  ! Input cloud top pressure
      REAL*8 GVMR(LLPAR)          ! Input Mixing Ratio on GEOS Grid
      REAL*8 AERPROF(LLPAR,NAER)  ! Input aerosol profile on GEOS Grid
      REAL*8 LAMDA                ! Wavelength


C  dimensioning

	INTEGER	MAXGRID
	PARAMETER	 (MAXGRID = 101)

C  Output
C  ------

C  output height grid from ground to TOA_HEIGHT
C    NGRID levels counting downward from top

	INTEGER		 NGRID
	REAL*8 ZGRID(MAXGRID), ZTEMPO(MAXGRID)
	REAL*8 DELTA(MAXGRID)
        REAL*8 NDENS(MAXGRID)   !Trace gas number density profile 
	REAL*8 AEREXT(MAXGRID)  !Total aerosol and cloud absorption
	REAL*8 AERSCA(MAXGRID)  !Total aerosol and cloud scattering
	REAL*8 PHASEFCN(MAXGRID,0:7) !Weighted average of arsl types

C  output Pressure and temperatures

	REAL*8 PRESS(MAXGRID)
	REAL*8 TEMPERATURE(MAXGRID)

C  output Trace gas quantities

	REAL*8 O3VMR(MAXGRID)
 	REAL*8 NO2VMR(MAXGRID)
        REAL*8 LVMR(MAXGRID)    !Trace gas profile on LIDORT grid

C  Maximum number of other gases

        INTEGER         MAX_OTHERGAS
        PARAMETER       ( MAX_OTHERGAS = 2 )
        REAL*8 OTHERGAS_VMR(MAXGRID,MAX_OTHERGAS)
	INTEGER		LOW, HIGH

C  USA data :
C  ----------

C   USA altitude levels                                               
C   USA v.m.r for O3  : Mid-latitude Spring/Summer                    
C   USA v.m.r for NO2                                                 

	INTEGER		 I
	REAL*8 USA_ZZZ(46), USA_VMR_O3(46), USA_VMR_NO2(46)
        DATA (USA_VMR_O3(I), I = 1, 46) / 
     $  4.0000D-01,  7.0000D-01,  7.5000D-01,  5.7000D-01,  2.0000D-01,
     $  1.9000D-01,  4.0000D-01,  8.0000D-01,  1.3000D+00,  1.8000D+00,
     $  2.8000D+00,  3.5000D+00,  4.5000D+00,  5.9000D+00,  7.5500D+00,
     $  8.7000D+00,  8.9000D+00,  8.1000D+00,  7.0000D+00,  6.0000D+00,
     $  4.8000D+00,  4.0000D+00,  3.4000D+00,  2.9000D+00,  2.4000D+00,
     $  2.0000D+00,  1.5000D+00,  1.0000D+00,  7.0000D-01,  6.0000D-01,
     $  5.0000D-01,  4.4000D-01,  3.0000D-01,  2.2300D-01,  1.7930D-01,
     $  1.3040D-01,  1.1110D-01,  9.1260D-02,  7.7640D-02,  6.4080D-02,
     $  5.5120D-02,  4.8210D-02,  4.2220D-02,  3.6940D-02,  3.3370D-02,
     $  3.0170D-02/
        DATA (USA_ZZZ(I), I = 1, 46) / 
     $  1.0000D+02,  9.5000D+01,  9.0000D+01,  8.5000D+01,  8.0000D+01,
     $  7.5000D+01,  7.0000D+01,  6.5000D+01,  6.0000D+01,  5.5000D+01,
     $  5.0000D+01,  4.7500D+01,  4.5000D+01,  4.2500D+01,  4.0000D+01,
     $  3.7500D+01,  3.5000D+01,  3.2500D+01,  3.0000D+01,  2.7500D+01,
     $  2.5000D+01,  2.4000D+01,  2.3000D+01,  2.2000D+01,  2.1000D+01,
     $  2.0000D+01,  1.9000D+01,  1.8000D+01,  1.7000D+01,  1.6000D+01,
     $  1.5000D+01,  1.4000D+01,  1.3000D+01,  1.2000D+01,  1.1000D+01,
     $  1.0000D+01,  9.0000D+00,  8.0000D+00,  7.0000D+00,  6.0000D+00,
     $  5.0000D+00,  4.0000D+00,  3.0000D+00,  2.0000D+00,  1.0000D+00,
     $  0.0000D+00/
        DATA (USA_VMR_NO2(I), I = 1, 46) / 
     $  1.7000D-04,  1.7600D-04,  1.8300D-04,  1.9200D-04,  2.0200D-04,
     $  2.1500D-04,  2.3100D-04,  2.5300D-04,  2.8500D-04,  3.3900D-04,
     $  4.4300D-04,  6.6600D-04,  1.1500D-03,  2.1700D-03,  4.0300D-03,
     $  6.2600D-03,  7.2800D-03,  7.2100D-03,  6.1600D-03,  4.8100D-03,
     $  3.7400D-03,  3.0600D-03,  2.5800D-03,  2.1600D-03,  1.7600D-03,
     $  1.3900D-03,  1.0600D-03,  7.7100D-04,  5.1900D-04,  3.1900D-04,
     $  1.7100D-04,  7.4800D-05,  4.4500D-05,  3.1500D-05,  2.6200D-05,
     $  2.3800D-05,  2.3200D-05,  2.3000D-05,  2.3000D-05,  2.3000D-05,
     $  2.3000D-05,  2.3000D-05,  2.3000D-05,  2.3000D-05,  2.3000D-05,
     $  2.3000D-05/

C  local variables
C  ---------------

	INTEGER		 J, K, L, M, N, cidx_low, cidx_hig
	REAL*8 TOA_HEIGHT, GROUND, TAILOFF, SMALL, ZN, SSA
	REAL*8 G, L0, L1, AOD, CBP, MPRESS(LLPAR), XG, DP
	PARAMETER	 ( SMALL = 1.0D-12 )
        REAL*8 TEMP(LLPAR)      !Aerosol profile on GEOS grid
        REAL*8 TAER(MAXGRID)    !Aerosol profile on LIDORT grid
        REAL*8 LAER(MAXGRID,NAER+1)  !Aerosol profile on LIDORT grid
	REAL*8 TAERSCA
	REAL*8 TPHASEFCN(MAXGRID,0:7)
	REAL*8 FP, P, AIR_DENSITY, DPK, ZLOWER, LEV_COT
	PARAMETER	( FP=7.24325825D+18 )

#if   defined( GEOS_4 ) || defined ( GEOS_5 )
        MSIGMAE(1) = ETA_A(1)/SPIN + ETA_B(1)
        DO L = 1, LLPAR-1
          MSIGMAE(L+1) = ETA_A(L+1)/SPIN + ETA_B(L+1)
          MSIGMA(L) = (MSIGMAE(L) + MSIGMAE(L+1) ) / 2.0d0
        ENDDO
        MSIGMAE(LLPAR+1) = 0.0d0
        MSIGMA(LLPAR)=(MSIGMAE(LLPAR)+MSIGMAE(LLPAR+1) )/2.0d0
#endif

C  Get lowest height
	DPK = DLOG(SPIN)
	CALL LAGRANGE_INTERP_MASTER
     I  ( .TRUE., .TRUE.,  
     I      USAML_LNP, USAML_ZZZ, NUSAML, NUSAML,
     I      1, 1, SMALL, 1, 1, DPK, ZLOWER )

	ZN = ZLOWER

C  set grid problem

        SSA = 0.9999D0
	DO I = 1, MAXGRID
	  DELTA(I) = 0.0D0
	  ZGRID(I) = 0.0D0
	  PRESS(I) = 0.0D0
	  TEMPERATURE(I) = 0.0D0
 	  NO2VMR(I) = 0.0D0
	  O3VMR(I) = 0.0D0
          NDENS(I) = 0.0D0
          AEREXT(I) = 1.0D-6
          AERSCA(I) = AEREXT(I)*SSA
          PHASEFCN(I,0) = 1.0D0
          DO J = 1,7
            PHASEFCN(I,J) = 0.0D0
          ENDDO
	ENDDO
	TOA_HEIGHT = 65.0D0
	I = 1
	ZTEMPO(1) = ZN
	DO WHILE (ZTEMPO(I).LT.17.5D0)
	  I = I + 1
	  ZTEMPO(I) = ZTEMPO(I-1) + 0.25D0
	ENDDO
	DO WHILE (ZTEMPO(I).LT.18.0D0)
	  I = I + 1
	  ZTEMPO(I) = ZTEMPO(I-1) + 0.5D0
	ENDDO
	DO WHILE (ZTEMPO(I).LT.20.0D0)
	  I = I + 1
	  ZTEMPO(I) = ZTEMPO(I-1) + 1.0D0
	ENDDO
	DO WHILE (ZTEMPO(I).LT.25.0D0)
	  I = I + 1
	  ZTEMPO(I) = ZTEMPO(I-1) + 2.5D0
	ENDDO
	DO WHILE (ZTEMPO(I).LT.55.0D0)
	  I = I + 1
	  ZTEMPO(I) = ZTEMPO(I-1) + 5.0D0
	ENDDO
        ZTEMPO(I+1) = TOA_HEIGHT
	NGRID = I + 1

C  get midpoints and layer thickness values

	DO I = 1, NGRID
	  ZGRID(I) = ZTEMPO(NGRID-I+1)
	ENDDO
	DO I = 1, NGRID
	  ZTEMPO(I) = ZGRID(I)
	ENDDO
	DO I = 1, NGRID - 1
	  DELTA(I) = ZTEMPO(I) - ZTEMPO(I+1)
	  ZGRID(I) = 0.5D0*(ZTEMPO(I)+ZTEMPO(I+1))
	ENDDO
	ZGRID(NGRID) = 0.0D0
	NGRID = NGRID - 1

C  get NO2, O3  from USA

 	CALL LAGRANGE_INTERP_MASTER
     I  ( .TRUE., .TRUE.,  
     I      USA_ZZZ, USA_VMR_NO2, 46, 46,
     I      1, 1, SMALL, MAXGRID, NGRID, ZGRID, NO2VMR )

	CALL LAGRANGE_INTERP_MASTER
     I  ( .TRUE., .TRUE.,  
     I      USA_ZZZ, USA_VMR_O3, 46, 46,
     I      1, 1, SMALL, MAXGRID, NGRID, ZGRID, O3VMR )

         DO I = 1, NGRID
           OTHERGAS_VMR(I,1) = 1.0D-06*O3VMR(I)
           OTHERGAS_VMR(I,2) = 1.0D-06*NO2VMR(I)
         ENDDO

C  get pressure and temperature from ML atmosphere

	CALL LAGRANGE_INTERP_MASTER
     I  ( .TRUE., .TRUE.,  
     I      USAML_ZZZ, USAML_TTT, NUSAML, NUSAML,
     I      1, 1, SMALL, MAXGRID, NGRID, ZGRID, TEMPERATURE )

	CALL LAGRANGE_INTERP_MASTER
     I  ( .TRUE., .FALSE.,  
     I      USAML_ZZZ, USAML_PPP, NUSAML, NUSAML,
     I      1, 1, SMALL, MAXGRID, NGRID, ZGRID, PRESS )

C----------------------------------------------------------------------
C  put trace gas profile onto lidort grid
C----------------------------------------------------------------------

C Compute pressure levels of trace gas profile
      DO I = 1,LLPAR
         MPRESS(I) = MSIGMA(I)*(SPIN - PTOP) + PTOP
      ENDDO

C Interpolate mixing ratio profile onto lidort grid
       CALL LAGRANGE_INTERP_MASTER
     I  ( .TRUE., .TRUE.,  
     I      MPRESS, GVMR, LLPAR, LLPAR,
     I      1, 1, SMALL, MAXGRID, NGRID, PRESS, LVMR )

        DO I = 1, NGRID
C          IF (LVMR(I) .LT. 0.0D0)
C    &         LVMR(I) = 1.0D-15 !Numerical Issues
           AIR_DENSITY = PRESS(I)*FP/TEMPERATURE(I)
           NDENS(I) = LVMR(I) * AIR_DENSITY
        ENDDO

C Insert stratosphere for NO2
        IF (LAMDA .GT. 400.0D0) THEN    
          I = 1
          DO WHILE (ZGRID(I) .GT. 19.0D0) 
             AIR_DENSITY = PRESS(I)*FP/TEMPERATURE(I)
             NDENS(I) = 1.0D-6*NO2VMR(I) * AIR_DENSITY
             I = I + 1
          ENDDO
        ENDIF

C  aerosol and clouds
C  ==================

C Interpolate aerosol profiles onto lidort grid, preserving the 
C   total column optical thickness
C Since the GEOS-CHEM aerosol profiles are optical thickness for
C   a given layer, first divide by the thickness of the layer
       DO N = 1, NAER
         DO I = 1, LLPAR
           DP = SPIN*(MSIGMAE(I)-MSIGMAE(I+1))
           TEMP(I) = AERPROF(I,N) / DP
         ENDDO

         CALL LAGRANGE_INTERP_MASTER
     I    ( .TRUE., .TRUE.,  
     I        MPRESS, TEMP, LLPAR, LLPAR,
     I        1, 1, SMALL, MAXGRID, NGRID, PRESS, TAER )
         
         LAER(1,N) = 0.0D0   !Top of lidort model
         DO I = 2, NGRID
           DP = PRESS(I) - PRESS(I-1)
           LAER(I,N) = TAER(I) * DP
         ENDDO
       ENDDO

       laer(1:ngrid,naer+1) = 0.0d0

       IF (CLDOD .GT. 0.0D0) THEN
C Determine cloud profile (adapted from code by tpk)
C Assume optical depth of 8 per 100 hPa
C If given CTP and CLDOD such that thickness would be greater than
C   8 per 100 hPa, set to 10 per 100 hPa to prevent unrealistically 
C   thick clouds in the boundary layer
       CBP = MIN(CTP+12.5D0*CLDOD,SPIN)

       cidx_low = MAXVAL(MINLOC(press(1:ngrid), 
     $      MASK=(press(1:ngrid) >= ctp) ))
       cidx_hig = MAXVAL(MAXLOC(press(1:ngrid), 
     $      MASK=(press(1:ngrid) <= cbp) ))

       if (press(ngrid) .le. ctp) then
         cidx_low = cidx_hig
       else
         lev_cot = cldod/ABS(press(cidx_hig)-press(cidx_low)) * 100.0D0
         cidx_low = MIN(cidx_low + 1, cidx_hig)
       endif

       IF (cidx_low .LT. cidx_hig) THEN 
         DO i = cidx_low, cidx_hig
          laer(i,naer+1) = lev_cot*(press(i)-press(i-1))/100.0d0
         ENDDO
       ELSE 
         LEV_COT = MIN(CLDOD,10.0D0)
         LAER(CIDX_LOW,NAER+1) = LEV_COT
       ENDIF
       ENDIF


C Interpolate aerosol parameters to wavelength LAMDA
       IF (LAMDA .GT. 400.0D0) THEN    
         LOW = 2
         HIGH = 3
       ELSE
         LOW = 1
         HIGH = 2
       ENDIF
       XG = (LAMDA-WL(LOW))/(WL(HIGH)-WL(LOW))

C      Get total extinction 
       DO N = 1, NAER + 1
         SSA = ASSA(LOW,N) + XG*(ASSA(HIGH,N)-ASSA(LOW,N))
         DO I = 2, NGRID
            IF (LAER(I,N) .GT. 1.0D-10) THEN
              ! aerext has units of km^-1
	      AEREXT(I) = AEREXT(I) + LAER(I,N)/(ZGRID(I-1)-ZGRID(I)) 	
	      AERSCA(I) = AERSCA(I) + 
     &                    SSA * LAER(I,N)/(ZGRID(I-1)-ZGRID(I)) 	
	    ENDIF
	  ENDDO
        ENDDO

       DO N = 1, NAER + 1
         SSA = ASSA(LOW,N) + XG*(ASSA(HIGH,N)-ASSA(LOW,N))
         DO I = 2, NGRID
            IF (LAER(I,N) .GT. 1.0D-10) THEN
	      TAERSCA = SSA * LAER(I,N)/(ZGRID(I-1)-ZGRID(I))

              DO M = 1,7
                TPHASEFCN(I,M) = PHFCN(M,LOW,N) + 
     &            XG*(PHFCN(M,HIGH,N)-PHFCN(M,LOW,N)) 

C Phase function is the weighted average as required by LIDORT
                PHASEFCN(I,M) = PHASEFCN(I,M) + 
     &             TPHASEFCN(I,M) * TAERSCA / AERSCA(I)
	      ENDDO
	    ENDIF
	  ENDDO
        ENDDO

C Test interpolation
c     OPEN ( FILE = 'test1.out',
c    &     UNIT = 83, STATUS = 'UNKNOWN', FORM = 'FORMATTED' )
c     OPEN ( FILE = 'test2.out',
c    &     UNIT = 84, STATUS = 'UNKNOWN', FORM = 'FORMATTED' )
c     do i = 1,llpar
c       write(83,*) gvmr(i),mpress(i)
c     enddo
c     do i = 2,ngrid
c       write(84,*) lvmr(i),press(i),aerext(i)*(zgrid(i-1)-zgrid(i))
c     enddo
c     CLOSE ( 83 ) ! Close output file
c     CLOSE ( 84 ) ! Close output file

C  Finish
	RETURN
	END
C ******************************************************************
C
      SUBROUTINE LINTP2 (N1,X1,Y1,N2,X2,Y2)

      INTEGER           N1,N2
      REAL*8  X1(N1),Y1(N1),X2(N2),Y2(N2)
      INTEGER           INV, N2P, J, I2, I1, J2, I
      REAL*8  X, YP, XP, S

      IF(N1.LE.0.OR.N2.LE.0.OR.N1.GT.32767.OR.N2.GT.32767)GOTO 9
C
C --- Checking the condition for inverting arrays ---
C
      INV=0
      IF (N1.EQ.1.OR.N2.EQ.1) GOTO 500
      IF ((X1(N1)-X1(1))*(X2(N2)-X2(1))) 300,300,500
C
C --- Inversion of new grid ---
C
  300 CONTINUE
      INV=1
      N2P=N2/2
      DO 301 J=1,N2P
      XP=X2(J)
      X2(J)=X2(N2-J+1)
      X2(N2-J+1)=XP
  301 CONTINUE
C
C --- Main block ---
C
  500 IF (N1.EQ.1) GOTO 7
      S=DSIGN(1.0D0,X1(N1)-X1(1))
      I2=1
      I1=2
    1 IF((X2(I2)-X1(1))*S) 2,2,3
    2 Y2(I2)=Y1(1)
      I2=I2+1
      IF(I2.GT.N2) GOTO 999
      GOTO 1
    3 IF((X2(I2)-X1(I1))*S) 4,4,5
    4 X=(X2(I2)-X1(I1-1))/(X1(I1)-X1(I1-1))
      Y2(I2)=Y1(I1)*X+Y1(I1-1)*(1.0D0-X)
      I2=I2+1
      IF(I2.GT.N2) GOTO 999
      GOTO 3
    5 I1=I1+1
      IF(I1.LE.N1) GOTO 3
      DO 6 J2=I2,N2
    6 Y2(J2)=Y1(N1)
      GOTO 999
    7 DO 8 I=1,N2
    8 Y2(I)=Y1(1)
  999 CONTINUE
C
C --- Checking the condition for back inversion ---
C
      IF (INV.NE.1) GOTO 1000
C
C --- Back ersion ---
C
      DO 302 J=1,N2P
      XP=X2(J)
      YP=Y2(J)
      X2(J)=X2(N2-J+1)
      Y2(J)=Y2(N2-J+1)
      X2(N2-J+1)=XP
      Y2(N2-J+1)=YP
  302 CONTINUE
C
C --- Exit block ---
C
 1000 RETURN
    9 PRINT 100,N1,N2
  100 FORMAT(10X,'  Error in subroutine LINTP2 !',2I15)
      STOP '  Error in subroutine LINTP2 !'
      END

        SUBROUTINE LAGRANGE_INTERP_MASTER
     &     ( L, DL, ZFIX, WFIX, MAXFIX, NFIX, N2, M2, FLOOR_VALUE,
     &       MAXOUT, NOUT, ZOUT, WOUT )

C  Four-point lagrange extrapolation routine (extracted from Lowtran).

C  Given grid ZFIX of actual size NFIX and external dimension MAXFIX, and
C  array of interpolates WFIX, interpolate onto grid ZOUT.

C  Given array ZOUT of actual size NOUT and external dimension MAXOUT,
C  this routine finds extrapolated points WOUT by trawling through
C  array ZFIX to find correct four points round each entry ZOUT(j), then 
C  calculating extrapolated values WOUT(j,n) from knowledge of four nearby
C  values of WFIX.  If Z-values coincide exactly, then interpolation is
C  skipped entirely, and answer is a straight copy from WFIX.

C  The second index N2 (total of M2 possible values) allows simultaneous 
C  interpolation (for example of concentrations). Logical arrays L & DL
C  (indexed by N) control whether interpolation is required for quantity
C  and whether interpolation is linear or log linear, respectively.

C  All data values in WFIX are assumed greater than or equal to zero.
C  Resulting WOUT values smaller than limiting value FLOOR_VALUE are
C  set to zero (input argument).

C  The LOWTRAN7 version assumed the values of ZFIX were in ascending
C  order.  This version has been modified to work with ZFIX either
C  ascending or descending---as long as it's monotonic.

C  Output array ZOUT must be monotonic in the same sense as ZFIX.

C  There are no internal dimension arguments.  There are no external
C  references - the module is entirely generic.

C       ----------------------------------------------------------------
C       Input
C       ----------------------------------------------------------------

C  MAXFIX is external dimension of WFIX; actual size NFIX

        INTEGER          N2, M2, NFIX, MAXFIX
        REAL*8 ZFIX (MAXFIX), WFIX (MAXFIX,M2)

C  L is an array of logical flags: L(I) is .TRUE. if we want W(I).
C  DL is an array of logical flags: DL(I) is .TRUE. if the
C  interpolation for variable I is linear. L(I) must be true.

        LOGICAL           L(M2), DL(M2)

C  Array of output Z-values and dimensioning

        INTEGER          NOUT, MAXOUT
        REAL*8 ZOUT (MAXOUT)

C  limiting value

	REAL*8 FLOOR_VALUE

C       ----------------------------------------------------------------
C       Output
C       ----------------------------------------------------------------

        REAL*8 WOUT (MAXOUT,M2)
	
C       ----------------------------------------------------------------
C       Local variables
C       ----------------------------------------------------------------

        LOGICAL          LEXACT, ASCENDING
        LOGICAL          EXTRAPOLATE_1, EXTRAPOLATE_4
        INTEGER          IX(4), ISTART_PREVIOUS, K
        REAL*8 DX(4), ZX(4)

C       ----------------------------------------------------------------
C       initialise
C       ----------------------------------------------------------------

C  determine ascending flag

        ASCENDING = ( ZFIX ( 1 ) .LT. ZFIX ( NFIX ) )

C       ----------------------------------------------------------------
C       Main loop over output abscissae
C       ----------------------------------------------------------------

	DO K = 1, NOUT

C  Call to get pointers and setup factors

	  ISTART_PREVIOUS = 999

	  CALL LAGRANGE_INTERP_1
     I    ( ZOUT(K), ZFIX, MAXFIX, NFIX, ASCENDING, 
     O      LEXACT, EXTRAPOLATE_1, EXTRAPOLATE_4,
     O      ISTART_PREVIOUS, IX, ZX, DX )

C  carry out interpolation

	  CALL LAGRANGE_INTERP_2
     I      ( IX, ZX, DX, L, DL, N2, M2, FLOOR_VALUE,
     I        LEXACT, EXTRAPOLATE_1, EXTRAPOLATE_4,
     I        WFIX, MAXFIX, ZOUT(K), K, MAXOUT,
     O        WOUT )

	ENDDO

C       ----------------------------------------------------------------
C       Finish
C       ----------------------------------------------------------------

	END

C

        SUBROUTINE LAGRANGE_INTERP_1
     I    ( Z, ZFIX, MAXFIX, NFIX, ASCENDING,
     O      LEXACT, EXTRAPOLATE_1, EXTRAPOLATE_4,
     O      ISTART_PREVIOUS, IX, ZX, DX )

C  Trawl routine. Given starting point ISTART, trawl through array ZFIX
C  to find correct four points ZX(4) around Z, indexed by pointers IX(4).
C  Output IX, ZX and additional array DX of weight factors.

C       ----------------------------------------------------------------
C       Input
C       ----------------------------------------------------------------

C  MAXFIX is for dimensioning; NFIX is the actual number of entries

        INTEGER          NFIX, MAXFIX
        REAL*8 ZFIX (MAXFIX), Z
	LOGICAL		 ASCENDING

C       ----------------------------------------------------------------
C       Output (ISTART_PREVIOUS is modified argument)
C       ----------------------------------------------------------------

        LOGICAL          LEXACT, EXTRAPOLATE_1, EXTRAPOLATE_4
        INTEGER          IX(4), ISTART_PREVIOUS
        REAL*8 DX(4), ZX(4)

C       ----------------------------------------------------------------
C       Local variables
C       ----------------------------------------------------------------

	LOGICAL		LOOP
        INTEGER         IM, ISTART, LOCAL_ISTART_PREV

C       ----------------------------------------------------------------
C       Trawling  code
C       ----------------------------------------------------------------

C  ASCENDING when .TRUE. means that the values in ZFIX are monotonically
C   ascending; they're monotonically descending if false.

        LEXACT = .FALSE.
        EXTRAPOLATE_1 = .FALSE.
        EXTRAPOLATE_4 = .FALSE.
	IF ( ISTART_PREVIOUS .EQ. 999 ) THEN
	  LOCAL_ISTART_PREV = 1
	ELSE
	  LOCAL_ISTART_PREV = ISTART_PREVIOUS
	ENDIF

C  Trawl through array ZFIX to find the value nearest to Z.
C  LEXACT when .TRUE. means Z exactly matches ZFIX (IX(3)), so that there's
C   no need to do any interpolation at all.

        IF ( ASCENDING ) THEN
	  LOOP = .TRUE.
	  IM = LOCAL_ISTART_PREV - 1
	  DO WHILE ( LOOP .AND. IM .LT. NFIX )
	    IM = IM + 1
            IF ( ZFIX(IM) .EQ. Z ) THEN
              IX(3) = IM
              LEXACT = .TRUE.
	      LOOP = .FALSE.
            ELSE IF ( ZFIX(IM) .GT. Z ) THEN
              IX(3) = IM
 	      LOOP = .FALSE.
            END IF
          ENDDO
        ELSE
	  LOOP = .TRUE.
	  IM = LOCAL_ISTART_PREV - 1
	  DO WHILE ( LOOP .AND. IM .LT. NFIX )
	    IM = IM + 1
            IF ( ZFIX(IM) .EQ. Z ) THEN
              IX(3) = IM
              LEXACT = .TRUE.
	      LOOP = .FALSE.
            ELSE IF ( ZFIX(IM) .LT. Z ) THEN
              IX(3) = IM
 	      LOOP = .FALSE.
            END IF
          ENDDO
        END IF
        IF (LOOP) IX(3) = NFIX

C       ----------------------------------------------------------------
C       Exactness - set starter.
C       ----------------------------------------------------------------

	IF ( LEXACT ) THEN

	  ISTART = IX(3)

C       ----------------------------------------------------------------
C       Lagrange coefficient determination (if not exact)
C       ----------------------------------------------------------------

        ELSE IF ( .NOT. LEXACT ) THEN

C  Endpoint correction: we have to take care of cases where Z is very
C near, or beyond, one of the ends of ZFIX.
C  If we're beyond the large-Z end of ZFIX, we extrapolate to get the
C point beyond the end we'll need.  If we're beyond the small-Z end,
C we just use the four points near the end.  This follows what LOWTRAN
C was doing; no doubt the idea is that at the top of the atmosphere,
C with small concentrations, the extrapolation is OK, but at the
C bottom, with larger concentrations, it could blow up.  But the code
C here has been written so that the checks on ASCENDING in this
C block may be removed, to get extrapolation at both ends,
C and all will work.

C         --------------------------------------------------------------
C         Lagrange coefficient pointers (if not exact)
C         --------------------------------------------------------------

C  general pointers

          IX(1) = IX(3) - 2
          IX(2) = IX(3) - 1
          IX(4) = IX(3) + 1
	  ISTART = IX(2)

C  Lower endpoint special case

          IF ( IX(1) .LT. 1 ) THEN

	    ISTART = 1
            IF ( .NOT. ASCENDING  .AND.  Z .GT. ZFIX(1) ) THEN
              EXTRAPOLATE_1 = .TRUE.     ! Lower endpoint extrapolation
              IX(2) = 1
              IX(3) = 2
              IX(4) = 3
            ELSE                         ! Lower endpoint correction
              IX(1) = 1
              IX(2) = 2
              IX(3) = 3
              IX(4) = 4
            END IF

C  Upper endpoint special case

          ELSE IF ( IX(4) .GT. NFIX ) THEN

	    ISTART = NFIX
            IF ( ASCENDING  .AND.  Z .GT. ZFIX(NFIX) ) THEN
              EXTRAPOLATE_4 = .TRUE.     ! Upper endpoint extrapolation
              IX(1) = NFIX - 2
              IX(2) = NFIX - 1
              IX(3) = NFIX
            ELSE                         ! Upper endpoint correction
              IX(1) = NFIX - 3
              IX(2) = NFIX - 2
              IX(3) = NFIX - 1
              IX(4) = NFIX
            END IF

          END IF

C         --------------------------------------------------------------
C         Lagrange coefficient determination (if required)
C         --------------------------------------------------------------

	  IF ( ISTART .NE. ISTART_PREVIOUS ) THEN

C  Abscissae

            ZX(2) = ZFIX (IX(2))
            ZX(3) = ZFIX (IX(3))
            IF ( EXTRAPOLATE_1 ) THEN
              ZX(1) = ZX(2) + 2.0 * ( Z - ZX(2) )
            ELSE
              ZX(1) = ZFIX(IX(1))
            END IF
            IF ( EXTRAPOLATE_4 ) THEN
              ZX(4) = ZX(3) + 2.0 * ( Z - ZX(3) )
            ELSE
              ZX(4) = ZFIX (IX(4))
            END IF

C  Set up interpolation coefficients

            DX(1) = (ZX(1)-ZX(2)) * (ZX(1)-ZX(3)) * (ZX(1)-ZX(4))
            DX(2) = (ZX(2)-ZX(3)) * (ZX(2)-ZX(4)) * (ZX(2)-ZX(1))
            DX(3) = (ZX(3)-ZX(4)) * (ZX(3)-ZX(1)) * (ZX(3)-ZX(2))
            DX(4) = (ZX(4)-ZX(1)) * (ZX(4)-ZX(2)) * (ZX(4)-ZX(3))

	  ENDIF

C       ----------------------------------------------------------------
C  	Finish (set starter pointer)
C       ----------------------------------------------------------------

        END IF

	ISTART_PREVIOUS = ISTART

	END


	SUBROUTINE LAGRANGE_INTERP_2
     I   ( IX, ZX, DX, L, DL, N2, M2, FLOOR_VALUE,
     I     LEXACT, EXTRAPOLATE_1, EXTRAPOLATE_4,
     I     WFIX, MAXFIX, Z, KOUT, MAXOUT,
     O     WOUT )
C       ----------------------------------------------------------------
C  	Module Input variables
C       ----------------------------------------------------------------

C  pointers to interpolation ordinates

	INTEGER		 IX(4)

C  Lagrangian abscissae and weight factors

	REAL*8 ZX(4), DX(4)

C  Exactness flag

        LOGICAL          LEXACT

C  end-point extrapolation flags

        LOGICAL          EXTRAPOLATE_1, EXTRAPOLATE_4

C  MAXFIX is  for declaring the first dimension of WFIX;

        INTEGER          N2, M2, MAXFIX
        REAL*8 WFIX (MAXFIX,M2)

C  Interpolation point

        REAL*8 Z
	INTEGER		 KOUT, MAXOUT

C  L is an array of logical flags: L(I) is .TRUE. if we want W(I) returned.

        LOGICAL          L(M2)

C  DL is an array of logical flags: DL(I) is .TRUE. if the
C  interpolation for variable I is linear. L(I) must be true first.

        LOGICAL          DL(M2)

C  Limiting value

	REAL*8 FLOOR_VALUE
	
C       ----------------------------------------------------------------
C  	Module Output variables
C       ----------------------------------------------------------------

        REAL*8 WOUT(MAXOUT,M2)

C       ----------------------------------------------------------------
C  	Local variables
C       ----------------------------------------------------------------

	INTEGER		 K
        REAL*8 A1, A2, A3, A4, X1, X2, X3, X4
        REAL*8 ARG, BIGEXP_VALUE
	PARAMETER	 ( BIGEXP_VALUE = 37.0 )

C  Initialise
C  ----------

	DO K = 1, N2
	  WOUT(KOUT,K) = 0.0
	ENDDO

C  Exact values
C  ------------

        IF ( LEXACT ) THEN

          DO K = 1, N2
	    IF ( L(K) ) THEN
              WOUT(KOUT,K) = WFIX (IX(3),K)
            ENDIF
	  ENDDO

C  Interpolation calculation
C  -------------------------

	ELSE

C  Lagrangian weights

          A1 = ((Z-ZX(2)) * (Z-ZX(3)) * (Z-ZX(4))) / DX(1)
          A2 = ((Z-ZX(3)) * (Z-ZX(4)) * (Z-ZX(1))) / DX(2)
          A3 = ((Z-ZX(4)) * (Z-ZX(1)) * (Z-ZX(2))) / DX(3)
          A4 = ((Z-ZX(1)) * (Z-ZX(2)) * (Z-ZX(3))) / DX(4)

C  For each quantity (second index N)

          DO K = 1, N2

	    IF ( L(K) ) THEN

C  .. set up ordinates X1, X2, X3, X4

              X2 = WFIX(IX(2),K)
              X3 = WFIX(IX(3),K)

              IF ( EXTRAPOLATE_1 ) THEN
                X1 = X2 + 2.0 * ( X2 - X3 )
              ELSE
                X1 = WFIX(IX(1),K)
              END IF

              IF ( EXTRAPOLATE_4 ) THEN
                X4 = X3 + 2.0 * ( X3 - X2 )
              ELSE
                X4 = WFIX(IX(4),K)
              END IF

C  .. Interpolation options

C  	* Lagrangian abandoned to linear if any X < 0
C  	* Lagrangian with linear values
C       * Lagrangian with Log-linear values (pressures!)

              IF (X1.LE.0.0.OR.X2.LE.0.0.OR.X3.LE.0.0.OR.X4.LE.0.0)THEN
                WOUT(KOUT,K) = (X3*(Z-ZX(2))-X2*(Z-ZX(3)))/(ZX(3)-ZX(2))
              ELSE IF ( DL(K) ) THEN
                WOUT(KOUT,K) = A1*X1 + A2*X2 + A3*X3 + A4*X4
              ELSE IF ( .NOT. DL(K) ) THEN
                X1 = DLOG(X1)
                X2 = DLOG(X2)
                X3 = DLOG(X3)
                X4 = DLOG(X4)
                ARG = A1*X1 + A2*X2 + A3*X3 + A4*X4
                ARG = MIN ( ARG, BIGEXP_VALUE )
                WOUT(KOUT,K) = DEXP ( ARG )
              END IF

            END IF

C  Floor value exceeded - set to zero

            IF (DABS(WOUT(KOUT,K)).LT.FLOOR_VALUE)WOUT(KOUT,K)=0.0

          ENDDO

	ENDIF

C  Finish

        END

