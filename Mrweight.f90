module Mrweight
!---------------------------------------------------------------------------
! Contains routine to compute the intensity at the top of the atmosphere
! for a Rayleigh scattering atmosphere with a Lambertian surface.
! 
! call rweight(wavelx,p_sfc,a_sfc,th0,th,dphix,refl)
!
!             Pepijn Veefkind, Henk Eskes, KNMI, 1999
!---------------------------------------------------------------------------
  implicit none

  private

  public :: rweight
 
contains


  subroutine csalbr(xtau,xalb)
    implicit none
    real,intent(in)  :: xtau
    real,intent(out) :: xalb
    xalb=0
    xalb=(3*xtau-fintexp3(xtau)*(4+2*xtau)+2*exp(-xtau))
    xalb=xalb/(4.+3*xtau)
  end subroutine csalbr

  function fintexp3(xtau)
    implicit none
    real,intent(in)  :: xtau
    real :: fintexp3
    real :: xx
    xx=(exp(-xtau)*(1.-xtau)+xtau*xtau*fint1exp(xtau))/2.
    fintexp3=xx
  end function fintexp3

  function fint1exp(xtau)
! accuracy 2e-07... for 0<xtau<1
    implicit none
    real,intent(in) :: xtau
    real :: fint1exp
    real :: xx,xftau
    real,dimension(0:5),parameter :: &
           a=(/-.57721566, 0.99999193,-0.24991055, &
               0.05519968,-0.00976004, 0.00107857 /)
    integer :: i
    xx=a(0)
    xftau=1.
    do i=1,5
       xftau=xftau*xtau
       xx=xx+a(i)*xftau
    enddo
    fint1exp=xx-log(xtau)
  end function fint1exp



  subroutine chand (xphi,xmuv,xmus,xtau,xrray)
!
! input parameters: xphi,xmus,xmuv,xtau
!   xphi: azimuthal difference between sun and observation (xphi=0,
!         in backscattering) and expressed in degree (0.:360.)
!   xmus: cosine of the sun zenith angle
!   xmuv: cosine of the observation zenith angle
!   xtau: molecular optical depth
! output parameter: 
!   xrray : molecular reflectance (0.:1.)
! constant : xdep: depolarization factor (0.0279)
!
    real,intent(in)  :: xphi,xmuv,xmus,xtau
    real,intent(out) :: xrray
!
    real :: xdep,pl(10)
    real :: fs0,fs1,fs2
    real :: as0(10),as1(2),as2(2)
    real :: fac,pi,phios,xcosf1,xcosf2
    real :: xcosf3,xbeta2,xfd,xph1,xph2,xph3,xitm, xp1, xp2, xp3
    real :: cfonc1,cfonc2,cfonc3,xlntau,xitot1,xitot2,xitot3
    integer   :: i
!
    as0 = (/.33243832,-6.777104e-02,.16285370 &
           ,1.577425e-03,-.30924818,-1.240906e-02,-.10324388 &
           ,3.241678e-02,.11493334,-3.503695e-02/)
    as1 = (/.19666292, -5.439061e-02/)
    as2 = (/.14545937, -2.910845e-02/)
!
!	pi=3.1415927
    fac=pi/180.
    phios=180.-xphi
    xcosf1=1.
    xcosf2=cos(phios*fac)
    xcosf3=cos(2*phios*fac)
    xbeta2=0.5
    xdep=0.0279
    xfd=xdep/(2-xdep)
    xfd=(1-xfd)/(1+2*xfd)
    xph1=1+(3*xmus*xmus-1)*(3*xmuv*xmuv-1)*xfd/8.
    xph2=-xmus*xmuv*sqrt(1-xmus*xmus)*sqrt(1-xmuv*xmuv)
    xph2=xph2*xfd*xbeta2*1.5
    xph3=(1-xmus*xmus)*(1-xmuv*xmuv)
    xph3=xph3*xfd*xbeta2*0.375
    xitm=(1-exp(-xtau*(1/xmus+1/xmuv)))*xmus/(4*(xmus+xmuv))
    xp1=xph1*xitm
    xp2=xph2*xitm
    xp3=xph3*xitm
    xitm=(1-exp(-xtau/xmus))*(1-exp(-xtau/xmuv))
    cfonc1=xph1*xitm
    cfonc2=xph2*xitm
    cfonc3=xph3*xitm
    xlntau=log(xtau)
    pl(1)=1.
    pl(2)=xlntau
    pl(3)=xmus+xmuv
    pl(4)=xlntau*pl(3)
    pl(5)=xmus*xmuv
    pl(6)=xlntau*pl(5)
    pl(7)=xmus*xmus+xmuv*xmuv
    pl(8)=xlntau*pl(7)
    pl(9)=xmus*xmus*xmuv*xmuv
    pl(10)=xlntau*pl(9)
    fs0=0.
    do i=1,10
       fs0=fs0+pl(i)*as0(i)
    enddo
    fs1=pl(1)*as1(1)+pl(2)*as1(2)
    fs2=pl(1)*as2(1)+pl(2)*as2(2)
    xitot1=xp1+cfonc1*fs0*xmus
    xitot2=xp2+cfonc2*fs1*xmus
    xitot3=xp3+cfonc3*fs2*xmus
    xrray=xitot1*xcosf1
    xrray=xrray+xitot2*xcosf2*2
    xrray=xrray+xitot3*xcosf3*2
    xrray=xrray/xmus
!
  end subroutine chand



  subroutine rweight(wavelx,p_sfc,a_sfc,th0,th,dphix,refl)
!-----------------------------------------------------------------------
! name       :weight
! description:computes radiance weight as a function of cloud fraction
! date       :18 August 1999
! version    :0.0
! authors    : J.P. Veefkind
!-----------------------------------------------------------------------
!  Variable declaration
!  Input:
!      Wavelength [nm]           : wavel0
!      Surface pressure [hPa]    : p_sfc 
!      Surface albedo            : a_sfc   
!      Solar zenith angle        : th0
!      Viewing zenith angle      : th
!      Relative azimuth angle    : dphix
!  Output:
!      total reflectance         : refl
!-----------------------------------------------------------------------
    implicit none
!
    real,intent(in)  :: wavelx, p_sfc, a_sfc, th0, th, dphix
    real,intent(out) :: refl
!
    real :: refl_a, t_mu, t_mu0, a_sph, a
    real :: mu, mu0, pi, dphi, wavel
    real :: tau_sfc, tau
!  
    pi=acos(-1.)
      
!-----------------------------------------------------------------------
!  Main
!-----------------------------------------------------------------------      
! Conversion of some parameters
! wavelength in micrometers:
    wavel=wavelx*1e-3
! compute mu and mu0
    mu0=cos(th0/180.*pi) 
    mu =cos(th /180.*pi)  
!
! use symmetry of dphi
    dphi=abs(dphix)
    if( dphi > 180.0 ) dphi=360.0-dphi
!
! redefine dphi
    dphi=180.-dphi
!
! Compute the rayleigh optical depth for surface and cloud top pressure
! approximation formula by Hansen and Travis      
    tau=0.008569 * wavel**(-4.) * &
               (1. + 0.0113 * wavel**(-2.) + 0.00013 * wavel**(-4.))
     
    tau_sfc= p_sfc / 1013.25 * tau

    tau=tau_sfc
    a=a_sfc
	
! compute the transmission, see Vermote and Tanre JQSRT 47, pp 305, 1992
    t_mu0=( (2./3.+mu0)+(2./3.-mu0)*exp(-tau/mu0) )/( (4./3.)+tau )
    t_mu =( (2./3.+mu )+(2./3.-mu )*exp(-tau/mu ) )/( (4./3.)+tau )	
	
! compute sherical albedo
    call csalbr(tau, a_sph)
	
! compute atm reflectance
    call chand (dphi,mu,mu0,tau,refl_a)

! compute total reflectance
    refl=refl_a + t_mu0 * t_mu* a / (1. - a_sph * a)

	
! Write result of RT to screen
!    print*,'dphix, dphi ',dphix, dphi
!    print*,'mu0, mu ',mu0, mu
!    write(*,'(a,f9.3)') 'Transmission mu0          ',t_mu0	
!    write(*,'(a,f9.3)') 'Transmission mu0          ',t_mu
!    write(*,'(a,f9.3)') 'Spherical albedo          ',a_sph
!    write(*,'(a,f9.3)') 'Atmospheric reflectance   ',refl_a
!    write(*,'(a,f9.3)') 'Reflectance               ',refl

  end subroutine rweight


end module Mrweight
