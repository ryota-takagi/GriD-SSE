module m_dc3d_driver
!<
!! Driver module for DC3D (Okada 1992)
!!
!!   Fault parameter
!!     lon     (degree)
!!     lat     (degree)
!!     dep     (km)
!!     strike  (degree)
!!     dip     (degree)
!!     length  (km)
!!     width   (km)
!!     slip    (mm)
!!     rake    (degree)
!!   Station
!!     stlo    (degree)
!!     stla    (degree)
!!     stdp    (km)
!!   Output
!!     u(3)    Displacement (mm); component 1=E, 2=N, 3=U
!!     e(3,3)  (strain)
!!
!>

!  use m_dc3d

  implicit none
  private

  public :: dc3d__ul  ! Upper Left corner (single station)
  public :: dc3d__ulr ! Upper Left corner (multi stations relative to reference stations)
  public :: dc3d__cm  ! Center Middle point (single station)
  public :: dc3d__cmr ! Center Middle point (multi stations relative to reference stations)

  real(8), parameter :: pid = acos(-1.0d0)/180.0d0

contains

  !-----------------------------------------------------------------------------
  subroutine dc3d__ul(alpha,lon,lat,dep,strike,dip,length,width,slip,rake,&
                      stlo,stla,stdp,u,e,is_strain)
  !-----------------------------------------------------------------------------
    real(8),  intent(in)  :: alpha,lon,lat,dep,strike,dip,length,width,slip,rake
    real(8),  intent(in)  :: stlo,stla,stdp
    real(8),  intent(out) :: u(3),e(3,3)
    logical,  intent(in)  :: is_strain
    real(8)               :: disl1,disl2
    real(8)               :: g(3,3),x,y
    real(8)               :: ux,uy,uz
    real(8)               :: uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
    real(8)               :: e1(3,3),e2(3,3)
    integer               :: iret
    call dc3d__sr2disl(slip,rake,disl1,disl2)
    call dc3d__ll2fxy(lon,lat,strike,stlo,stla,x,y,g)
    call dc3d(alpha,x,y,-stdp,dep,dip,0.0d0,length,-width,0.0d0,disl1,disl2,0.0d0,&
              ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret)
    call dc3d__uf2g(ux,uy,uz,g,u)
    if(is_strain)then
      call dc3d__ef2g(uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,g,e)
    endif
  end subroutine dc3d__ul

  !-----------------------------------------------------------------------------
  subroutine dc3d__cm(alpha,lon,lat,dep,strike,dip,length,width,slip,rake,&
                      stlo,stla,stdp,u,e,is_strain)
  !-----------------------------------------------------------------------------
    real(8),  intent(in)  :: alpha,lon,lat,dep,strike,dip,length,width,slip,rake
    real(8),  intent(in)  :: stlo,stla,stdp
    real(8),  intent(out) :: u(3),e(3,3)
    logical,  intent(in)  :: is_strain
    real(8)               :: disl1,disl2
    real(8)               :: g(3,3),x,y
    real(8)               :: ux,uy,uz
    real(8)               :: uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
    real(8)               :: e1(3,3),e2(3,3)
    integer               :: iret
    call dc3d__sr2disl(slip,rake,disl1,disl2)
    call dc3d__ll2fxy(lon,lat,strike,stlo,stla,x,y,g)
    call dc3d(alpha,x,y,-stdp,dep,dip,-length*0.5d0,length*0.5d0,-width*0.5d0,width*0.5d0,disl1,disl2,0.0d0,&
              ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret)
    call dc3d__uf2g(ux,uy,uz,g,u)
    if(is_strain)then
      call dc3d__ef2g(uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,g,e)
    endif
  end subroutine dc3d__cm

  !-----------------------------------------------------------------------------
  subroutine dc3d__ulr(alpha,lon,lat,dep,strike,dip,length,width,slip,rake,&
                       ns,stlo,stla,stdp,iuse,u,e,nr,rlo,rla,rdp,is_strain)
  !-----------------------------------------------------------------------------
    real(8),  intent(in)  :: alpha,lon,lat,dep,strike,dip,rake,length,width,slip
    integer,  intent(in)  :: ns
    real(8),  intent(in)  :: stlo(ns),stla(ns),stdp(ns)
    integer,  intent(in)  :: iuse(ns)
    real(8),  intent(out) :: u(:,:),e(:,:,:)
    integer,  intent(in)  :: nr
    real(8),  intent(in)  :: rlo(nr),rla(nr),rdp(nr)
    logical,  intent(in)  :: is_strain
    real(8)               :: disl1,disl2
    real(8)               :: g(3,3),x,y
    real(8)               :: ux,uy,uz
    real(8)               :: uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
    real(8)               :: e1(3,3),e2(3,3)
    integer               :: iret
    real(8)               :: ur(3),ur1(3)
    integer               :: ir,is
    call dc3d__sr2disl(slip,rake,disl1,disl2)
    call dc3d__g2f(strike,g)
    ! Reference stations
    ur(1:3) = 0.0d0
    do ir=1,nr
      call dc3d__ll2fxy1(lon,lat,g,rlo(ir),rla(ir),x,y)
      call dc3d(alpha,x,y,-rdp(ir),dep,dip,0.0d0,length,-width,0.0d0,disl1,disl2,0.0d0,&
                ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret)
      call dc3d__uf2g(ux,uy,uz,g,ur1)
      ur(1:3) = ur(1:3) + ur1(1:3)
    enddo
    if(nr>=1) ur(1:3) = ur(1:3)/dble(nr)
    ! Observation stations
    !$omp parallel private(x,y,ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret)
    !$omp do
    do is=1,ns
      if(iuse(is)==0)then
        u(1:3,is)     = 0.0d0
        e(1:3,1:3,is) = 0.0d0
        cycle
      endif
      call dc3d__ll2fxy1(lon,lat,g,stlo(is),stla(is),x,y)
      call dc3d(alpha,x,y,-stdp(is),dep,dip,0.0d0,length,-width,0.0d0,disl1,disl2,0.0d0,&
                ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret)
      call dc3d__uf2g(ux,uy,uz,g,u(1:3,is))
      u(1:3,is) = u(1:3,is) - ur(1:3)
      if(is_strain)then
        call dc3d__ef2g(uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,g,e(1:3,1:3,is))
      endif
    enddo
    !$omp end do
    !$omp end parallel
  end subroutine dc3d__ulr

  !-----------------------------------------------------------------------------
  subroutine dc3d__cmr(alpha,lon,lat,dep,strike,dip,length,width,slip,rake,&
                           ns,stlo,stla,stdp,iuse,u,e,nr,rlo,rla,rdp,is_strain)
  !-----------------------------------------------------------------------------
    real(8),  intent(in)  :: alpha,lon,lat,dep,strike,dip,length,width,slip,rake
    integer,  intent(in)  :: ns
    real(8),  intent(in)  :: stlo(ns),stla(ns),stdp(ns)
    integer,  intent(in)  :: iuse(ns)
    real(8),  intent(out) :: u(:,:),e(:,:,:)
    integer,  intent(in)  :: nr
    real(8),  intent(in)  :: rlo(nr),rla(nr),rdp(nr)
    logical,  intent(in)  :: is_strain
    real(8)               :: disl1,disl2
    real(8)               :: g(3,3),x,y
    real(8)               :: ux,uy,uz
    real(8)               :: uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
    real(8)               :: e1(3,3),e2(3,3)
    integer               :: iret
    real(8)               :: ur(3),ur1(3)
    integer               :: ir,is
    call dc3d__sr2disl(slip,rake,disl1,disl2)
    call dc3d__g2f(strike,g)
    ! Reference stations
    ur(1:3) = 0.0d0
    do ir=1,nr
      call dc3d__ll2fxy1(lon,lat,g,rlo(ir),rla(ir),x,y)
      call dc3d(alpha,x,y,-rdp(ir),dep,dip,-length*0.5d0,length*0.5d0,-width*0.5d0,width*0.5d0,disl1,disl2,0.0d0,&
                ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret)
      call dc3d__uf2g(ux,uy,uz,g,ur1)
      ur(1:3) = ur(1:3) + ur1(1:3)
    enddo
    if(nr>=1) ur(1:3) = ur(1:3)/dble(nr)
    ! Observation stations
    !$omp parallel private(x,y,ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret)
    !$omp do
    do is=1,ns
      if(iuse(is)==0)then
        u(1:3,is)     = 0.0d0
        e(1:3,1:3,is) = 0.0d0
        cycle
      endif
      call dc3d__ll2fxy1(lon,lat,g,stlo(is),stla(is),x,y)
      call dc3d(alpha,x,y,-stdp(is),dep,dip,-length*0.5d0,length*0.5d0,-width*0.5d0,width*0.5d0,disl1,disl2,0.0d0,&
                ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret)
      call dc3d__uf2g(ux,uy,uz,g,u(1:3,is))
      u(1:3,is) = u(1:3,is) - ur(1:3)
      if(is_strain)then
        call dc3d__ef2g(uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,g,e(1:3,1:3,is))
      endif
    enddo
    !$omp end do
    !$omp end parallel
  end subroutine dc3d__cmr

  !-----------------------------------------------------------------------------
  !! private subroutine
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  subroutine dc3d__sr2disl(slip,rake,disl1,disl2)
  !<
  !! Convert slip/rake to DISL1(strike slip)/DISL2(dip slip)
  !>
  !-----------------------------------------------------------------------------
    real(8),  intent(in)  :: slip,rake
    real(8),  intent(out) :: disl1,disl2
    disl1 = slip*cos(rake*pid)
    disl2 = slip*sin(rake*pid)
  end subroutine dc3d__sr2disl

  !-----------------------------------------------------------------------------
  subroutine dc3d__g2f(strike,g)
  !<
  !! Rotation matrix to convert XY coordinate to fault coordinate
  !!  XY coordinate
  !!    - X: East
  !!    - Y: North
  !!    - Z: Up
  !!  Fault coordinate
  !!    - X: Strike
  !!    - Y: Updip
  !!    - Z: Up
  !>
  !-----------------------------------------------------------------------------
    real(8),  intent(in)  :: strike
    real(8),  intent(out) :: g(3,3)
    real(8)               :: sn,cs
    sn = sin(strike*pid)
    cs = cos(strike*pid)
    g(1,1) =  sn
    g(1,2) =  cs
    g(2,1) = -cs
    g(2,2) =  sn
    g(1,3) = 0.0d0
    g(2,3) = 0.0d0
    g(3,1) = 0.0d0
    g(3,2) = 0.0d0
    g(3,3) = 1.0d0
  end subroutine dc3d__g2f

  !-----------------------------------------------------------------------------
  subroutine dc3d__ll2fxy(lon,lat,strike,stlo,stla,x,y,g)
  !<
  !! Convert lon lat to fault coordinate
  !!  Fault coordinate
  !!    - X: Strike
  !!    - Y: Updip
  !!    - Z: Up
  !>
  !-----------------------------------------------------------------------------
    real(8),  intent(in)  :: lon,lat,strike
    real(8),  intent(in)  :: stlo,stla
    real(8),  intent(out) :: x,y,g(3,3)
    real(8)               :: x1,y1
    call dc3d__g2f(strike,g)
    call dc3d__pltxy(stla,stlo,x1,y1,0,lat,lon,0)
    x = g(1,1)*x1 + g(1,2)*y1
    y = g(2,1)*x1 + g(2,2)*y1
  end subroutine dc3d__ll2fxy

  !-----------------------------------------------------------------------------
  subroutine dc3d__ll2fxy1(lon,lat,g,stlo,stla,x,y)
  !<
  !! Convert lon lat to fault coordinate
  !>
  !-----------------------------------------------------------------------------
    real(8),  intent(in)  :: lon,lat,g(3,3)
    real(8),  intent(in)  :: stlo,stla
    real(8),  intent(out) :: x,y
    real(8)               :: x1,y1
    call dc3d__pltxy(stla,stlo,x1,y1,0,lat,lon,0)
    x = g(1,1)*x1 + g(1,2)*y1
    y = g(2,1)*x1 + g(2,2)*y1
  end subroutine dc3d__ll2fxy1

  !-----------------------------------------------------------------------------
  subroutine dc3d__uf2g(ux,uy,uz,g,u)
  !<
  !! Convert displacement from fault coordinate to XY coordinate
  !>
  !-----------------------------------------------------------------------------
    real(8),  intent(in)  :: ux,uy,uz,g(3,3)
    real(8),  intent(out) :: u(3)
    u(1) = g(1,1)*ux + g(2,1)*uy
    u(2) = g(1,2)*ux + g(2,2)*uy
    u(3) = uz
  end subroutine dc3d__uf2g

  !-----------------------------------------------------------------------------
  subroutine dc3d__ef2g(uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,g,e)
  !<
  !! Convert strain from fault coordinate to XY coordinate
  !>
  !-----------------------------------------------------------------------------
    real(8),  intent(in)  :: uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,g(3,3)
    real(8),  intent(out) :: e(3,3)
    real(8)               :: e1(3,3),e2(3,3)
    integer               :: j1,j2,l
    e1(1,1) = uxx*1.0d-6
    e1(1,2) = (uxy+uyx)*0.5d-6
    e1(1,3) = (uzx+uxz)*0.5d-6
    e1(2,1) = e1(1,2)
    e1(2,2) = uyy*1.0d-6
    e1(2,3) = (uyz+uzy)*0.5d-6
    e1(3,1) = e1(1,3)
    e1(3,2) = e1(2,3)
    e1(3,3) = uzz*1.0d-6
    do j1=1,3
      do j2=1,3
        e2(j2,j1) = 0.0d0
        do l=1,3
          e2(j2,j1) = e2(j2,j1) + g(l,j2)*e1(l,j1)
        enddo
      enddo
    enddo
    do j1=1,3
      do j2=1,3
        e(j2,j1) = 0.0d0
        do l=1,3
          e(j2,j1) = e(j2,j1) + e2(j2,l)*g(l,j1)
        enddo
      enddo
    enddo
  end subroutine dc3d__ef2g

  !-----------------------------------------------------------------------------
  subroutine dc3d__pltxy(alat,along,x,y,ind,alat0,alng0,icord)
  !<
  !! coordinate transformation between geodetic and local cartesian
  !!   pltxy transforms (x,y) to (alat,along) if ind.eq.1
  !!   pltxy transforms (alat,along) to (x,y) if ind.eq.0
  !!   when icord.ne.0  pltxy makes no change in transformation  between
  !!                    (x,y) and (alat,along).
  !<
  !-----------------------------------------------------------------------------
    implicit none
    integer :: ind,icord
    real(8) :: alat,along,x,y,alat0,alng0
    real(8) :: a,e2,e12,d,rd
    real(8) :: rlat,slat,clat,v2,al,ph1,rph1,rph2,r,an,c1,c2
    real(8) :: rlato,slato,tphi1,cphi1
    ! GRS80
    real(8), parameter :: R_EARTH = 6371.0d0
    real(8), parameter :: A_EARTH = 6378.137d0
    real(8), parameter :: F_EARTH = 1.0d0/298.257222101d0
    real(8), parameter :: B_EARTH  = A_EARTH*(1.0d0-F_EARTH)
    real(8), parameter :: E_EARTH  = sqrt(F_EARTH*(2.0d0-F_EARTH))
    real(8), parameter :: E2_EARTH = E_EARTH**2
    real(8), parameter :: D2R = acos(-1.0d0)/180.0d0
    real(8), parameter :: R2D = 1.0d0/D2R
    a   = A_EARTH
    e2  = E2_EARTH
    e12 = E2_EARTH/(1.0d0-E2_EARTH)
    d   = R2D
    rd  = D2R
    if(ind==0)then
      if(icord/=0) then
        x = alat
        y = along
        return
      end if
      rlat = rd*alat
      slat = dsin(rlat)
      clat = dcos(rlat)
      v2   = 1.0d0 + e12*clat**2
      al   = along-alng0
      ph1  = alat + (v2*al**2*slat*clat)/(2.0d0*d )
      rph1 = ph1*rd
      rph2 = (ph1 + alat0)*0.5d0*rd
      r    = a*(1.0d0-e2)/dsqrt((1.0d0-e2*dsin(rph2)**2)**3)
      an   = a/dsqrt(1.0d0-e2*dsin(rph1)**2)
      c1   = d/r
      c2   = d/an
      y    = (ph1-alat0)/c1
      x    = (al*clat)/c2+(al**3*clat*dcos(2.0d0*rlat))/(6.0d0*c2*d**2)
      return
    elseif(ind==1)then
      if(icord/=0) then
        alat  = x
        along = y
        return
      end if
      rlato = alat0*rd
      slato = dsin(rlato)
      r     = a*(1.0d0-e2)/dsqrt((1.0d0-e2*slato**2)**3)
      an    = a/dsqrt(1.0d0-e2*slato**2)
      v2    = 1.0d0 + e12*dcos(rlato)**2
      c1    = d/r
      c2    = d/an
      ph1   = alat0+c1*y
      rph1  = ph1*rd
      tphi1 = dtan(rph1)
      cphi1 = dcos(rph1)
      alat  = ph1-(c2*x)**2*v2*tphi1/(2.0d0*d)
      along = alng0+c2*x/cphi1-(c2*x)**3*(1.0d0+2.0d0*tphi1**2)/(6.0d0*d**2*cphi1)
      return
    else
      return
    endif
  end subroutine dc3d__pltxy

end module m_dc3d_driver
