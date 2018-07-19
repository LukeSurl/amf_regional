module MAmfLut
  !
  ! Module to load the height-dependent air-mass factor lookup table
  ! and to compute (interpolate) the AMFs
  !
  !  subroutine ReadAmfLut
  !  subroutine GetAmf( pres,azi,view,sza,albedo,spres,amf ) 
  !
  implicit none

  private

  public :: ReadAmfLut, GetAmf, FitWindowCentre
   
  ! -------------------------------------------------------------
  ! filename of the air-mass factor lookup table
  ! -------------------------------------------------------------         
  character(160),parameter :: &
       AmfPath = '/as/home/ctm/kfb/fordylan'
       
  character(16),parameter :: amf_file = 'hcho_amf_lut.hdf'      
  ! -------------------------------------------------------------
  ! memory block containing the lookup table for air-mass factors
  ! -------------------------------------------------------------
  integer,parameter :: amfmpres=24
  integer,parameter :: amfmazi=10
  integer,parameter :: amfmview=9
  integer,parameter :: amfmsza=13
  integer,parameter :: amfmalbedo=10
  integer,parameter :: amfmspres=10
  real,dimension(amfmpres,amfmazi,amfmview,amfmsza,amfmalbedo,amfmspres) :: amflut
  real,dimension(amfmpres)  :: amfpress   ! pressure levels (TM3)
  real,dimension(amfmazi)   :: amfazis    ! azimuth angles
  real,dimension(amfmview)  :: amfviews   ! viewing angles
  real,dimension(amfmsza)   :: amfszas    ! solar zenith angles
  real,dimension(amfmalbedo):: amfalbedos ! albedo's
  real,dimension(amfmspres) :: amfspress  ! surface pressures 

  ! --------------------------------------------------------------
  ! FitWindowCentre: central wavelength of fit window
  ! --------------------------------------------------------------

  real,parameter  :: FitWindowCentre=348.0


contains


  subroutine ReadAmfLut
    !=======================================================================
    !
    ! ReadAmfLut:    read lookup table for AMF calculations
    !                updated HDF format lookup table (may 2004)
    !
    !                             Henk Eskes, KNMI, nov 1999
    !                             Folkert Boersma, KNMI, may 2004
    !=======================================================================
    implicit none
  
    ! include HDF parameters
    include 'hdf.f90'

    integer :: status
    integer :: sfstart, sffattr, sfgainfo, sfrattr, sfselect, sfrdata, sfend
    integer :: num_type, cnt
    integer :: sd_id = -1    ! identifier for scia HDF file
    integer :: attr_index, sds_id, sds_index
  
    integer, dimension(6) :: start, stride, edge
  
    character(len=18) :: buff_ver
    character(len=50) :: buff_wvl
    integer :: wvl

    start(1:6)  = 0 
    stride(1:6) = 1 
    edge(1:6)   = (/ amfmpres,amfmazi,amfmview,amfmsza,amfmalbedo,amfmspres /) 
        
  !integer :: status
  !integer :: hopen, vfstart
  !integer :: file_id = -1    ! identifier for scia HDF file

    !--------------------------------------------------------------
    ! Open HDF file for reading and initialize s interface
    ! --------------------------------------------------------------
    sd_id = sfstart(amf_file, dfacc_read)    
    if (sd_id == -1) then
       write(*,*) 'ReadAmfLut: Failed to open: ',amf_file
       stop
    else
       write(*,*) 'ReadAmfLut: File ', amf_file,' opened'
    end if
    attr_index = sffattr(sd_id,'Wavelength:')  
    status = sfrattr(sd_id, attr_index, buff_wvl)
    if (status == -1) then
       write(*,*) 'ReadAmfLut: Failed to read wavelength attribute'
    else   
       write(*,*) 'Wavelength attribute: ',buff_wvl
    end if
    
    attr_index = sffattr(sd_id,'Pressure levels:')  
    status = sfrattr(sd_id, attr_index, amfpress)
    if (status == -1) then
       write(*,*) 'ReadAmfLut: Failed to read pressure levels'
    else   
       write(*,*) 'Pressure levels: ',amfpress
    end if
    
    attr_index = sffattr(sd_id,'mu_0 (cosine sza) levels: ')  
    status = sfrattr(sd_id, attr_index, amfszas)
    if (status == -1) then
       write(*,*) 'ReadAmfLut: Failed to read mu_0 levels'
    else   
       write(*,*) 'mu_0 levels: ',amfszas
    end if
    
    attr_index = sffattr(sd_id,'mu (cosine vza) levels: ')  
    status = sfrattr(sd_id, attr_index, amfviews)
    if (status == -1) then
       write(*,*) 'ReadAmfLut: Failed to read mu levels'
    else   
       write(*,*) 'mu levels: ',amfviews
    end if
    
    attr_index = sffattr(sd_id,'rel. azimuth angles: ')  
    status = sfrattr(sd_id, attr_index, amfazis)
    if (status == -1) then
       write(*,*) 'ReadAmfLut: Failed to read rel. azimuth angles'
    else   
       write(*,*) 'rel. azimuth angles: ',amfazis
    end if
    
    attr_index = sffattr(sd_id,'albedo levels: ')  
    status = sfrattr(sd_id, attr_index, amfalbedos)
    if (status == -1) then
       write(*,*) 'ReadAmfLut: Failed to read albedo levels'
    else   
       write(*,*) 'albedo levels: ',amfalbedos
    end if
    
    attr_index = sffattr(sd_id,'surface pressure levels: ')  
    status = sfrattr(sd_id, attr_index, amfspress)
    if (status == -1) then
       write(*,*) 'ReadAmfLut: Failed to read surface pressure levels'
    else   
       write(*,*) 'surface pressure levels: ',amfspress
    end if

    !---------------------------------------------------------------
    ! Select HDF scientific data set
    ! --------------------------------------------------------------
    sds_index = 0
    sds_id = sfselect(sd_id, sds_index)
    if (sds_id == -1) then
       write(*,*) 'ReadAmfLut: Failed to select sds'
       stop
    else
       write(*,*) 'ReadAmfLut: Succeeded in selecting sds'
    end if
    
    !---------------------------------------------------------------
    ! Read in the data set
    ! --------------------------------------------------------------
    status = sfrdata(sds_id, start, stride, edge, amflut)
    if (status == -1) then
       write(*,*) 'ReadAmfLut: Failed to read hdf data set'
       stop
    else
       write(*,*) 'ReadAmfLut: Succeeded in reading hdf data set'
       write(*,*) 'ReadAmfLut: amf (min,max) = ', &
            '(',minval(amflut),',',maxval(amflut),')'
    end if
    
    !---------------------------------------------------------------
    ! Close HDF file 
    ! --------------------------------------------------------------
    status =sfend(sd_id)
    if (status == -1) then
       write(*,*) 'ReadAmfLut: Failed to close: ',amf_file
       stop
    else
       write(*,*) 'ReadAmfLut: File ', amf_file,' closed'
    end if

  end subroutine ReadAmfLut

  subroutine GetAmf( pres,azi1,view1,sza1,albedo,spres,amf1 ) 
    !=======================================================================
    !
    ! GetAmf:    Return AMF (divided by the geometrical AMF) 
    !            interpolated from the lookup table
    !
    ! pres   : model pressure level in hPa
    ! azi1   : azimuth angle (degrees) [-360..360]
    ! view1  : viewing angle (degrees) [-60..60]
    ! sza1   : solar zenith angle (degrees) [-90..90]
    ! albedo : ground albedo el. [0..1]
    ! spres  : surface pressure in hPa
    !
    ! CHANGE:
    !          In GetAmf the model predicted pressure 'pres' is 
    !          now replaced bu amfpres(1) if the pres > amfpress(1)
    !          and extrapolation would occur. This is done because
    !          on January 28, 2003, the ECMWF surface pressure
    !          of 1047 hPa exceeds the highest amfpres of 1045.94
    !          and the program would otherwise stop. KFB - 3/11/2003
    !           
    !          Identical for spres, since on January 28, a surface
    !          pressure of 1053.89 hPa was obtained. KFB - 3/11/2003
    !
    !                             Henk Eskes, KNMI, nov 1999
    !=======================================================================
    implicit none
    real,intent(in)   :: pres,spres,view1,azi1,sza1,albedo
    real,intent(out)  :: amf1
    !
    real              :: pres_tmp, spres_tmp
    integer           :: ipres,ispres,iview,iazi,isza,ialbedo
    integer           :: i,i1,i2,i3,i4,i5,i6
    real,dimension(2) :: xpres,xspres,xview,xazi,xsza,xalbedo
    real              :: azi,view,sza
    real,parameter    :: pi = 3.1415927
    !
    if ( pres > spres ) then
       write(*,*)'WARNING GetAmf: pressure > surface pressure'
    end if
    if ( pres < 0.0 ) then
       write(*,*)'ERROR GetAmf: negative pressure'
       stop
    end if
    if ( pres > 1100.0 ) then
       write(*,*)'ERROR GetAmf: unrealistic high pressure'
       stop
    end if
    if ( spres < 0.0 ) then
       write(*,*)'ERROR GetAmf: negative surface pressure'
       stop
    end if
    if ( spres > 1100.0 ) then
       write(*,*)'ERROR GetAmf: unrealistic high surfacepressure'
       stop
    end if
    !
    if ( azi1 > 360.1 .or. azi1 < -360.1 ) then
       write(*,*)'ERROR GetAmf: unrealistic large azimuth angle'
       stop
    end if
    azi=abs(azi1)
    if( azi > 180.0 ) azi=360.0-azi

    if ( view1 > 180.1 .or. view1 < -180.1 ) then
       write(*,*)'ERROR GetAmf: unrealistic large viewing angle'
       stop
    end if
    view=cos(pi*view1/180.0)

    if ( sza1 > 180.1 .or. sza1 < -0.001 ) then
       write(*,*)'ERROR GetAmf: unrealistic large solar angle'
       stop
    end if

    sza=cos(pi*sza1/180.0)
    !
    ! Determine reference index (ipres,iview,iazi,isza,ialbedo)

    ! Pressure: (decreasing)

    pres_tmp = pres
    if ( pres < amfpress(amfmpres) ) pres_tmp = amfpress(amfmpres)
    if ( pres > amfpress(1) ) then
       pres_tmp = amfpress(1)       
       write(*,*)'Model pressure: ',pres, &
            ' AMF first pressure level: ',amfpress(1)
       write(*,*)'WARNING: Overwriting Model pressure with AMF', &
            ' first pressure level'
       write(*,*)'for getting lowest level amf...'
    end if
    do i=1,amfmpres-1
       ipres = i
       if( amfpress(ipres+1) < pres_tmp ) exit 
    end do
    xpres(1)=(pres_tmp-amfpress(ipres+1))/(amfpress(ipres)-amfpress(ipres+1))
    xpres(2)=1.0-xpres(1)       

    ! Azimuth angle: (increasing)

    if( (azi < amfazis(1)) .or. (azi > amfazis(amfmazi)) )then
       write(*,*) "ERROR GetAmf: azimuth angle out of range"
       stop
    end if
    do i=1,amfmazi-1
       iazi = i
       if( amfazis(iazi+1) > azi ) exit 
    end do
    xazi(1)=(amfazis(iazi+1)-azi)/(amfazis(iazi+1)-amfazis(iazi))
    xazi(2)=1.0-xazi(1)

    ! Viewing angle mu: (increasing)

    if( (view < amfviews(1)) .or. (view > amfviews(amfmview)) )then
       write(*,*) "ERROR GetAmf: viewing angle out of range"
       write(*,*) "ERROR GetAmf: allowed cos(vza) range ",amfviews(1), &
            "-",amfviews(amfmview)
       stop
    end if
    do i=1,amfmview-1
       iview = i
       if( amfviews(iview+1) > view ) exit 
    end do
    xview(1)=(amfviews(iview+1)-view)/(amfviews(iview+1)-amfviews(iview))
    xview(2)=1.0-xview(1)

    ! Solar zenith angle: (decreasing)

    if( (sza < amfszas(amfmsza)) .or. (sza > amfszas(1)) )then
       write(*,*) "ERROR GetAmf: solar-zenith angle out of range"
       write(*,*) "ERROR GetAmf: sza:",sza
       write(*,*) "ERROR GetAmf: allowed sza range:", amfszas(amfmsza), &
            "-", amfszas(1)
       stop   
    end if
    do i=1,amfmsza-1
       isza = i
       if( amfszas(isza+1) < sza ) exit 
    end do
    xsza(1)=(sza-amfszas(isza+1))/(amfszas(isza)-amfszas(isza+1))
    xsza(2)=1.0-xsza(1)

    ! Ground albedo: (increasing)

    if( (albedo < amfalbedos(1)) .or. (albedo > amfalbedos(amfmalbedo)) )then
       write(*,*) "ERROR GetAmf: albedo out of range"
       stop
    end if
    do i=1,amfmalbedo-1
       ialbedo = i
       if( amfalbedos(ialbedo+1) > albedo ) exit 
    end do
    xalbedo(1)=(amfalbedos(ialbedo+1)-albedo)/ &
         (amfalbedos(ialbedo+1)-amfalbedos(ialbedo))
    xalbedo(2)=1.0-xalbedo(1)

    ! Surface pressure: (decreasing)

    if( (spres < amfspress(amfmspres)) )then
       write(*,*) amfspress(amfmspres),amfspress(1)
       write(*,*) "ERROR GetAmf: surface pressure out of range"
       write(*,*) "ERROR GetAmf: spres = ",spres
       stop
    end if
    if( spres > amfspress(1) ) then
       !write(*,*)'WARNING: Overwriting surface pressure with AMF first surface pressure level'
       !write(*,*)'for getting lowest level amf...'
       spres_tmp = amfpress(1)
       
       do i=1,amfmspres-1
          ispres = i
          if( amfspress(ispres+1) < spres_tmp ) exit 
       end do
       
       xspres(1)=(spres_tmp-amfspress(ispres+1))/ &
            (amfspress(ispres)-amfspress(ispres+1))
       xspres(2)=1.0-xspres(1)       
    else
       do i=1,amfmspres-1
          ispres = i
          if( amfspress(ispres+1) < spres ) exit 
       end do
       xspres(1)=(spres-amfspress(ispres+1))/ &
            (amfspress(ispres)-amfspress(ispres+1))
       xspres(2)=1.0-xspres(1)
    end if
    !
    ! Linear interpolation
    !
    amf1=0.0
    do i1=1,2
       do i2=1,2
          do i3=1,2
             do i4=1,2
                do i5=1,2
                   do i6=1,2
                      amf1=amf1+xpres(i1)*xazi(i2)*xview(i3)* &
                           xsza(i4)*xalbedo(i5)*xspres(i6)* &
                           amflut(ipres+i1-1,iazi+i2-1,iview+i3-1, &
                           isza+i4-1,ialbedo+i5-1,ispres+i6-1)
                   end do
                end do
             end do
          end do
       end do
    end do
    !
  end subroutine GetAmf


end module MAmfLut
