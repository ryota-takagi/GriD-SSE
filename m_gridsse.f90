module m_gridsse
!<
!! GriD-SSE module for timeseries fitting to estimate a rectangler fault
!!
!! 2018/12/XX Takagi et al. (2019, JGR)
!!
!!   Fault parameters
!!     flt(1): longitude
!!     flt(2): latitude
!!     flt(3): depth
!!     flt(4): strike
!!     flt(5): dip
!!     flt(6): length
!!     flt(7): width
!!     flt(8): slip
!!     flt(9): rake
!>

  use m_dc3d_driver
  use m_rtrend
  use m_inv
  use m_sort

  implicit none
  public  :: gridsse__mh87     ! Nonliner inversion (modified from Matsu'ura & Hasegawa 1987)
  private :: gridsse__partial
  private :: gridsse__ramp

  real(8), parameter :: alpha = 2.0d0/3.0d0

contains

  !-----------------------------------------------------------------------------
  subroutine gridsse__mh87(ns,stlo,stla,stdp,iuse,nr,rlo,rla,rdp,  &
                           ndmax,nd,td,ud,eud,t1,t2,               &
                           flt,eflt,dflt,wv,nitr,ak,               &
                           nsvr,rss,res,vr,rssall,resall,vrall)
  !-----------------------------------------------------------------------------
    ! Input & output parameters
    integer,  intent(in)    :: ns           ! number of stations
    real(8),  intent(in)    :: stlo(ns)     ! station longitude
    real(8),  intent(in)    :: stla(ns)     ! station latitude
    real(8),  intent(in)    :: stdp(ns)     ! station depth
    integer,  intent(in)    :: iuse(ns)     ! station flag (1: use, 0:NOT use)
    integer,  intent(in)    :: nr           ! number of reference stations
    real(8),  intent(in)    :: rlo(nr)      ! reference station longitude
    real(8),  intent(in)    :: rla(nr)      ! reference station latitude
    real(8),  intent(in)    :: rdp(nr)      ! reference station depth
    integer,  intent(in)    :: ndmax        ! maximum number of data points
    integer,  intent(in)    :: nd(ns)       ! number of data points
    real(8),  intent(in)    :: td(:,:)      ! td(nd(is),ns): time (day)
    real(8),  intent(in)    :: ud(:,:,:)    ! ud(nd(is),nc,ns): displacement (m)
    real(8),  intent(in)    :: eud(:,:,:)   ! eud(nd(is),nc,ns): S.D. of displacement (m)
    real(8),  intent(in)    :: t1           ! onset time of ramp function
    real(8),  intent(in)    :: t2           ! end time of ramp function
    real(8),  intent(inout) :: flt(9)       ! fault parameters (initial model for input)
    real(8),  intent(inout) :: eflt(9)      ! S.D. of fault paramters (a priori model error for input)
    real(8),  intent(in)    :: dflt(9)      ! infinitesimal difference of model parameters for partial derivative
    real(8),  intent(in)    :: wv           ! weighting facter for vertical component
    integer,  intent(in)    :: nitr         ! number for iteration
    real(8),  intent(in)    :: ak(nitr)     ! scaling parameter for update
    integer,  intent(in)    :: nsvr         ! number of stations for calculating variance reduction
    real(8),  intent(out)   :: rss          ! residual sum of square for nsvr stations
    real(8),  intent(out)   :: res          ! weighted residual sum of square sum[w*(d-c)**2]/sum[w] for nsvr stations
    real(8),  intent(out)   :: vr           ! variance reduction for nsvr stations
    real(8),  intent(out)   :: rssall       ! rss for all statinos used
    real(8),  intent(out)   :: resall       ! res for all stations used
    real(8),  intent(out)   :: vrall        ! vr for all stations used
    ! Data preparation
    integer                 :: nu,nc,nm,ni
    integer                 :: isu(ns)
    integer                 :: iflt(9)
    real(8)                 :: wud(ndmax,3,ns),ud1(ndmax),wsum(3,ns)
    real(8)                 :: ft(ndmax,ns)
    real(8)                 :: dc(3,ns),cc(3,ns),ac(3,ns)
    real(8)                 :: u(3,ns),du(3,ns,9)
    ! Nonlinear inversion
    real(8),  allocatable   :: mx(:),m0(:),Wm(:,:)
    real(8),  allocatable   :: G(:,:),H(:,:),E(:,:)
    real(8),  allocatable   :: b(:),c(:),r(:)
    ! VR calculation
    integer                 :: idx(ns)
    real(8)                 :: uabs(ns)
    real(8)                 :: var,wsm
    ! Cost function
    real(8)                 :: S(0:nitr)
    real(8),  allocatable   :: wm1(:),mx1(:)
    ! Zero infinitesimal difference to skip partial derivative
    real(8),  parameter     :: dflt0(9) = (/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
    ! Loop
    integer                 :: is,iu,ic,im,id,itr,i

    ! Vertical component
    if(wv<=0.0d0) nc = 2    ! use only 2 horizontal components
    if(wv>0.0d0)  nc = 3    ! use both 2 horizontal & vertical components

    ! Stations used
    iu = 0
    do is=1,ns
      if(iuse(is)==1)then
        iu      = iu + 1
        isu(iu) = is
      endif
    enddo
    nu = iu

    ! number of total data points
    ni = nu*nc

    ! Weighting factor that depends on data variance
    !$omp parallel do
    do is=1,ns
      if(iuse(is)==0)then
        wsum(1:nc,is) = 0.0d0
        cycle
      endif
      wud(1:nd(is),1:nc,is) = 1.0d0/eud(1:nd(is),1:nc,is)**2
      do ic=1,nc
        wsum(ic,is) = sum(wud(1:nd(is),ic,is))
      enddo
    enddo
    !$omp end parallel do

    ! Slip time function
    !$omp parallel do
    do is=1,ns
      if(iuse(is)==0)cycle
      call gridsse__ramp(nd(is),td(1:nd(is),is),t1,t2,ft(1:nd(is),is))
      call rtrend__rtrend(nd(is),ft(1:nd(is),is),td(1:nd(is),is))
    enddo
    !$omp end parallel do

    ! Cross and auto correlations between slip time function and data
    !$omp parallel do private(ud1)
    do is=1,ns
      if(iuse(is)==0)then
        dc(1:nc,is) = 0.0d0
        cc(1:nc,is) = 0.0d0
        ac(1:nc,is) = 0.0d0
        cycle
      endif
      do ic=1,nc
        ud1(1:nd(is)) = ud(1:nd(is),ic,is)
        call rtrend__rtrend(nd(is),ud1(1:nd(is)),td(1:nd(is),is))
        dc(ic,is) = sum( ud1(1:nd(is))   * ud1(1:nd(is))    * wud(1:nd(is),ic,is) )
        cc(ic,is) = sum( ft(1:nd(is),is) * ud1(1:nd(is))    * wud(1:nd(is),ic,is) )
        ac(ic,is) = sum( ft(1:nd(is),is) * ft(1:nd(is),is)  * wud(1:nd(is),ic,is) )
      enddo
      if(nc==3)then
        cc(3,is) = cc(3,is)*wv
        ac(3,is) = ac(3,is)*wv
      endif
    enddo
    !$omp end parallel do

    ! The number of model parameters (if model error of a parameter = 0, the mode parameter is fixed.)
    nm = 0
    do im=1,9
      if(eflt(im)/=0.0d0)then
        nm       = nm + 1
        iflt(nm) = im
      endif
    enddo

    ! Allocation memory
    allocate(mx(nm),m0(nm),Wm(nm,nm))
    allocate(G(ni,nm),H(ni,nm),E(nm,nm))
    allocate(b(nm),c(nm),r(nm))
    allocate(wm1(nm),mx1(nm))

    ! Initial value & model variance
    Wm(1:nm,1:nm) = 0.0d0
    do im=1,nm
      mx(im)    = flt(iflt(im))
      m0(im)    = flt(iflt(im))
      Wm(im,im) = 1.0d0/eflt(iflt(im))**2
      wm1(im)   = Wm(im,im)
    enddo

    ! Iteration
    do itr=1,nitr

      ! Provisional value & partial derivative
      do im=1,nm
        flt(iflt(im)) = mx(im)
      enddo
      call gridsse__partial(alpha,flt,dflt,ns,stlo,stla,stdp,iuse,u,du,nr,rlo,rla,rdp)

      ! Initial cost function
      if(itr==1)then
        S(0) = sum(dc(1:nc,1:ns))                    &
             + sum(u(1:nc,1:ns)**2*ac(1:nc,1:ns))    &
             - 2.0d0*sum(u(1:nc,1:ns)*cc(1:nc,1:ns)) &
             + sum(wm1(1:nm)*(m0(1:nm)-mx(1:nm))**2)
      endif

      ! Green's function
      !$omp parallel do private(id)
      do im=1,nm
        b(im) = 0.0d0
        c(im) = 0.0d0
        do iu=1,nu
          do ic=1,nc
            id = iu + nu*(ic-1)
            G(id,im) = du(ic,isu(iu),iflt(im))                        ! G : partial derivative
            H(id,im) = G(id,im)*ac(ic,isu(iu))                        ! H = (Rt*Wd*R)*G
            b(im)    = b(im) + G(id,im)*cc(ic,isu(iu))                ! b = (Gt*Rt*Wd)*d
            c(im)    = c(im) + G(id,im)*ac(ic,isu(iu))*u(ic,isu(iu))  ! c = (Gt*Rt*Wd)*f
          enddo
        enddo
      enddo
      !$omp end parallel do

      ! E = At*Wd*A + Wm
      E(1:nm,1:nm) = Wm(1:nm,1:nm)
      call DGEMM("T","N",nm,nm,ni,1.0d0,G,ni,H,ni,1.0d0,E,nm)
      ! E = inv(At*Wd*A+Wm)
      call inv(nm,E)
      ! r = (b-c) + Wm*(x0-x(k))
      do im=1,nm
        r(im) = (b(im)-c(im)) + Wm(im,im)*(m0(im)-mx(im))
      enddo

      ! Optimize the correction term to reduce the cost function
      mx1(1:nm) = mx(1:nm)
      do i=1,20
        mx(1:nm) = mx1(1:nm)
        call DGEMM("N","N",nm,1,nm,ak(itr)/2.0d0**(i-1),E,nm,r,nm,1.0d0,mx,nm)
        flt(iflt(1:nm)) = mx(1:nm)
        call gridsse__partial(alpha,flt,dflt0,ns,stlo,stla,stdp,iuse,u,du,nr,rlo,rla,rdp)
        S(itr) = sum(dc(1:nc,1:ns))                    &
               + sum(u(1:nc,1:ns)**2*ac(1:nc,1:ns))    &
               - 2.0d0*sum(u(1:nc,1:ns)*cc(1:nc,1:ns)) &
               + sum(wm1(1:nm)*(m0(1:nm)-mx(1:nm))**2)
        if(S(itr)<=S(itr-1)) exit
      enddo

    enddo

    ! Final value
    do im=1,nm
      flt(iflt(im))  = mx(im)
      eflt(iflt(im)) = sqrt(E(im,im))
    enddo

    ! Residual & VR for all stations used
    call gridsse__partial(alpha,flt,dflt0,ns,stlo,stla,stdp,iuse,u,du,nr,rlo,rla,rdp)
    rssall = sum(dc(1:nc,1:ns))                    &
           + sum(u(1:nc,1:ns)**2*ac(1:nc,1:ns))    &
           - 2.0d0*sum(u(1:nc,1:ns)*cc(1:nc,1:ns))
    vrall  = (1.0d0-rssall/sum(dc(1:nc,1:ns)))*100.0d0
    resall = rssall/sum(wsum(1:nc,1:ns))

    ! Residual & VR for top nsvr with the largest model displacements
    do is=1,ns
      idx(is) = is
      if(iuse(is)==0)then
        uabs(is) = 0.0d0
      else
        uabs(is) = sqrt( sum( u(1:nc,is)**2 ) )
      endif
    enddo
    call sort__descend(ns,uabs,idx)
    rss = 0.0d0
    var = 0.0d0
    wsm = 0.0d0
    do is=1,nsvr
      rss = rss + sum(dc(1:nc,idx(is)))                       &
                + sum(u(1:nc,idx(is))**2*ac(1:nc,idx(is)))    &
                - 2.0d0*sum(u(1:nc,idx(is))*cc(1:nc,idx(is)))
      var = var + sum(dc(1:nc,idx(is)))
      wsm = wsm + sum(wsum(1:nc,idx(is)))
    enddo
    vr  = (1.0d0-rss/var)*100.0d0
    res = rss/wsm

    deallocate(mx,m0,Wm)
    deallocate(G,H,E)
    deallocate(b,c,r)

  end subroutine gridsse__mh87


  !-----------------------------------------------------------------------------
  subroutine gridsse__partial(alpha,flt,dflt,ns,stlo,stla,stdp,iuse,u,du,nr,rlo,rla,rdp)
  !<
  !! Partial derivatives with respective to model parameters
  !!
  !! if dflt = 0 for a model paramer, the partial derivative for the parameter is not computed.
  !>
  !-----------------------------------------------------------------------------
    real(8),  intent(in)  :: alpha        ! "alpha" in DC3D by Okada(1992)
    real(8),  intent(in)  :: flt(9)       ! fault parameters
    real(8),  intent(in)  :: dflt(9)      ! infinitesimal difference for partial derivative
    integer,  intent(in)  :: ns           ! number of stations
    real(8),  intent(in)  :: stlo(ns)     ! longitude
    real(8),  intent(in)  :: stla(ns)     ! latitude
    real(8),  intent(in)  :: stdp(ns)     ! depth
    integer,  intent(in)  :: iuse(ns)     ! station flag (iuse=1:use, 0:not use)
    real(8),  intent(out) :: u(:,:)       ! displacement u(3,ns)
    real(8),  intent(out) :: du(:,:,:)    ! partial derivative du(3,ns,9)
    integer,  intent(in)  :: nr           ! number of reference stations
    real(8),  intent(in)  :: rlo(nr)      ! longitude of reference stations
    real(8),  intent(in)  :: rla(nr)      ! latigude of reference stations
    real(8),  intent(in)  :: rdp(nr)      ! depth of reference stations
    real(8)               :: f1(9),f2(9),u1(3,ns),u2(3,ns),e(3,3,ns)
    integer               :: jf,is,ic
    call dc3d__cmr(alpha,flt(1),flt(2),flt(3),flt(4),flt(5),flt(6),flt(7),flt(8),flt(9),&
                   ns,stlo,stla,stdp,iuse,u(1:3,1:ns),e,nr,rlo,rla,rdp,.false.)
    do jf=1,9
      if(dflt(jf)/=0.0d0)then
        f1(1:9) = flt(1:9)
        f2(1:9) = flt(1:9)
        f1(jf)  = f1(jf) - dflt(jf)
        f2(jf)  = f2(jf) + dflt(jf)
        call dc3d__cmr(alpha,f1(1),f1(2),f1(3),f1(4),f1(5),f1(6),f1(7),f1(8),f1(9),&
                       ns,stlo,stla,stdp,iuse,u1,e,nr,rlo,rla,rdp,.false.)
        call dc3d__cmr(alpha,f2(1),f2(2),f2(3),f2(4),f2(5),f2(6),f2(7),f2(8),f2(9),&
                       ns,stlo,stla,stdp,iuse,u2,e,nr,rlo,rla,rdp,.false.)
        do is=1,ns
          if(iuse(is)==0)cycle
          do ic=1,3
            du(ic,is,jf) = (u2(ic,is)-u1(ic,is))/(2.0d0*dflt(jf))
          enddo
        enddo
      endif
    enddo
  end subroutine gridsse__partial

  !-----------------------------------------------------------------------------
  subroutine gridsse__ramp(n,t,t1,t2,ramp)
  !<
  !!              _____    1
  !!             /.    .
  !!        ____/ .    .   0
  !!           t1 t2
  !>
  !-----------------------------------------------------------------------------
    integer,  intent(in)  :: n
    real(8),  intent(in)  :: t(:),t1,t2
    real(8),  intent(out) :: ramp(:)
    real(8)               :: a
    integer               :: i
    a = 1.0d0/(t2-t1)
    do i=1,n
      if(t(i)<t1)then
        ramp(i) = 0.0d0
      elseif(t(i)>=t2)then
        ramp(i) = 1.0d0
      else
        ramp(i) = a*(t(i)-t1)
      endif
    enddo
  end subroutine gridsse__ramp

end module m_gridsse
