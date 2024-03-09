MODULE module_bl_mynn_1d



  use module_bl_mynn_common,only:                        &
        cp        , cpv       , cliq       , cice      , &
        p608      , ep_2      , ep_3       , gtr       , &
        grav      , g_inv     , karman     , p1000mb   , &
        rcp       , r_d       , r_v        , rk        , &
        rvovrd    , svp1      , svp2       , svp3      , &
        xlf       , xlv       , xls        , xlscp     , &
        xlvcp     , tv0       , tv1        , tref      , &
        zero      , half      , one        , two       , &
        onethird  , twothirds , tkmin      , t0c       , &
        tice      , kind_phys


  IMPLICIT NONE

!===================================================================
! From here on, these are MYNN-specific parameters:
! The parameters below depend on stability functions of module_sf_mynn.
  real(kind_phys), parameter :: cphm_st=5.0, cphm_unst=16.0, &
                                cphh_st=5.0, cphh_unst=16.0

! Closure constants
  real(kind_phys), parameter ::  &
       &pr  =  0.74,             &
       &g1  =  0.235,            &  ! NN2009 = 0.235
       &b1  = 24.0,              &
       &b2  = 15.0,              &  ! CKmod     NN2009
       &c2  =  0.729,            &  ! 0.729, & !0.75, &
       &c3  =  0.340,            &  ! 0.340, & !0.352, &
       &c4  =  0.0,              &
       &c5  =  0.2,              &
       &a1  = b1*( 1.0-3.0*g1 )/6.0, &
!       &c1  = g1 -1.0/( 3.0*a1*b1**(1.0/3.0) ), &
       &c1  = g1 -1.0/( 3.0*a1*2.88449914061481660), &
       &a2  = a1*( g1-c1 )/( g1*pr ), &
       &g2  = b2/b1*( 1.0-c3 ) +2.0*a1/b1*( 3.0-2.0*c2 )

  real(kind_phys), parameter ::  &
       &cc2 =  1.0-c2,           &
       &cc3 =  1.0-c3,           &
       &e1c =  3.0*a2*b2*cc3,    &
       &e2c =  9.0*a1*a2*cc2,    &
       &e3c =  9.0*a2*a2*cc2*( 1.0-c5 ), &
       &e4c = 12.0*a1*a2*cc2,    &
       &e5c =  6.0*a1*a1

! Constants for min tke in elt integration (qmin), max z/L in els (zmax), 
! and factor for eddy viscosity for TKE (Kq = Sqfac*Km):
  real(kind_phys), parameter :: qmin=0.0, zmax=1.0, Sqfac=3.0
! Note that the following mixing-length constants are now specified in mym_length
!      &cns=3.5, alp1=0.23, alp2=0.3, alp3=3.0, alp4=10.0, alp5=0.2

  real(kind_phys), parameter :: gpw=5./3., qcgmin=1.e-8, qkemin=1.e-12
  real(kind_phys), parameter :: tliq = 269. !all hydrometeors are liquid when T > tliq

! Constants for cloud PDF (mym_condensation)
  real(kind_phys), parameter :: rr2=0.7071068, rrp=0.3989423

  !>Use Canuto/Kitamura mod (remove Ric and negative TKE) (1:yes, 0:no)
  !!For more info, see Canuto et al. (2008 JAS) and Kitamura (Journal of the 
  !!Meteorological Society of Japan, Vol. 88, No. 5, pp. 857-864, 2010).
  !!Note that this change required further modification of other parameters
  !!above (c2, c3). If you want to remove this option, set c2 and c3 constants 
  !!(above) back to NN2009 values (see commented out lines next to the
  !!parameters above). This only removes the negative TKE problem
  !!but does not necessarily improve performance - neutral impact.
  real(kind_phys), parameter :: CKmod=1.

  !Option to activate environmental subsidence in mass-flux scheme
  logical, parameter :: env_subs = .false.


CONTAINS

!=======================================================================
!     SUBROUTINE  mym_initialize:
!
!     Input variables:
!       iniflag         : <>0; turbulent quantities will be initialized
!                         = 0; turbulent quantities have been already
!                              given, i.e., they will not be initialized
!       nx, nz          : Dimension sizes of the
!                         x and z directions, respectively
!       tref            : Reference temperature                      (K)
!       dz(nz)          : Vertical grid spacings                     (m)
!                         # dz(nz)=dz(nz-1)
!       zw(nz+1)        : Heights of the walls of the grid boxes     (m)
!                         # zw(1)=0.0 and zw(k)=zw(k-1)+dz(k-1)
!       exner(nx,nz)    : Exner function at zw*h+zg             (J/kg K)
!                         defined by c_p*( p_basic/1000hPa )^kappa
!                         This is usually computed by integrating
!                         d(pi0)/dz = -h*g/tref.
!       rmo(nx)         : Inverse of the Obukhov length         (m^(-1))
!       flt, flq(nx)    : Turbulent fluxes of potential temperature and
!                         total water, respectively:
!                                    flt=-u_*Theta_*             (K m/s)
!                                    flq=-u_*qw_*            (kg/kg m/s)
!       ust(nx)         : Friction velocity                        (m/s)
!       pmz(nx)         : phi_m-zeta at z1*h+z0, where z1 (=0.5*dz(1))
!                         is the first grid point above the surafce, z0
!                         the roughness length and zeta=(z1*h+z0)*rmo
!       phh(nx)         : phi_h at z1*h+z0
!       u, v(nx,nz)     : Components of the horizontal wind        (m/s)
!       thl(nx,nz)      : Liquid water potential temperature
!                                                                    (K)
!       qw(nx,nz)       : Total water content Q_w                (kg/kg)
!
!     Output variables:
!       ql(nx,nz)       : Liquid water content                   (kg/kg)
!       vt, vq(nx,nz)   : Functions for computing the buoyancy flux
!       qke(nx,nz)      : Twice the turbulent kinetic energy q^2
!                                                              (m^2/s^2)
!       tsq(nx,nz)      : Variance of Theta_l                      (K^2)
!       qsq(nx,nz)      : Variance of Q_w
!       cov(nx,nz)      : Covariance of Theta_l and Q_w              (K)
!       el(nx,nz)       : Master length scale L                      (m)
!                         defined on the walls of the grid boxes
!
!     Work arrays:        see subroutine mym_level2
!       pd?(nx,nz,ny) : Half of the production terms at Level 2
!                         defined on the walls of the grid boxes
!       qkw(nx,nz,ny) : q on the walls of the grid boxes         (m/s)
!
!     # As to dtl, ...gh, see subroutine mym_turbulence.
!
!-------------------------------------------------------------------

!>\ingroup gsd_mynn_edmf
!! This subroutine initializes the mixing length, TKE, \f$\theta^{'2}\f$,
!! \f$q^{'2}\f$, and \f$\theta^{'}q^{'}\f$.
!!\section gen_mym_ini GSD MYNN-EDMF mym_initialize General Algorithm
!> @{
  SUBROUTINE  mym_initialize (                                &
       &            kts,kte,xland,                            &
       &            dz, dx, zw,                               &
       &            u, v, thl, qw,                            &
!       &            ust, rmo, pmz, phh, flt, flq,             &
       &            zi, theta, thetav, sh, sm,                &
       &            ust, rmo, el,                             &
       &            Qke, Tsq, Qsq, Cov, Psig_bl, cldfra_bl1D, &
       &            bl_mynn_mixlength,                        &
       &            edmf_w1,edmf_a1,                          &
       &            INITIALIZE_QKE,                           &
       &            spp_pbl,rstoch_col)
!
!-------------------------------------------------------------------

    integer, intent(in)           :: kts,kte
    integer, intent(in)           :: bl_mynn_mixlength
    logical, intent(in)           :: INITIALIZE_QKE
!    real(kind_phys), intent(in)   :: ust, rmo, pmz, phh, flt, flq
    real(kind_phys), intent(in)   :: rmo, Psig_bl, xland
    real(kind_phys), intent(in)   :: dx, ust, zi
    real(kind_phys), dimension(kts:kte),   intent(in) :: dz
    real(kind_phys), dimension(kts:kte+1), intent(in) :: zw
    real(kind_phys), dimension(kts:kte),   intent(in) :: u,v,thl,&
         &qw,cldfra_bl1D,edmf_w1,edmf_a1
    real(kind_phys), dimension(kts:kte),   intent(out) :: tsq,qsq,cov
    real(kind_phys), dimension(kts:kte),   intent(inout) :: el,qke
    real(kind_phys), dimension(kts:kte) ::                       &
         &ql,pdk,pdt,pdq,pdc,dtl,dqw,dtv,                        &
         &gm,gh,sm,sh,qkw,vt,vq
    integer :: k,l,lmax
    real(kind_phys):: phm,vkz,elq,elv,b1l,b2l,pmz=1.,phh=1.,     &
         &flt=0.,fltv=0.,flq=0.,tmpq
    real(kind_phys), dimension(kts:kte) :: theta,thetav
    real(kind_phys), dimension(kts:kte) :: rstoch_col
    integer ::spp_pbl

!> - At first ql, vt and vq are set to zero.
    DO k = kts,kte
       ql(k) = 0.0
       vt(k) = 0.0
       vq(k) = 0.0
    END DO
!
!> - Call mym_level2() to calculate the stability functions at level 2.
    CALL mym_level2 ( kts,kte,                      &
         &            dz,                           &
         &            u, v, thl, thetav, qw,        &
         &            ql, vt, vq,                   &
         &            dtl, dqw, dtv, gm, gh, sm, sh )
!
!   **  Preliminary setting  **

    el (kts) = 0.0
    IF (INITIALIZE_QKE) THEN
       !qke(kts) = ust**2 * ( b1*pmz )**(2.0/3.0)
       qke(kts) = 1.5 * ust**2 * ( b1*pmz )**(2.0/3.0)
       DO k = kts+1,kte
          !qke(k) = 0.0
          !linearly taper off towards top of pbl
          qke(k)=qke(kts)*MAX((ust*700. - zw(k))/(MAX(ust,0.01)*700.), 0.01)
       ENDDO
    ENDIF
!
    phm      = phh*b2 / ( b1*pmz )**(1.0/3.0)
    tsq(kts) = phm*( flt/ust )**2
    qsq(kts) = phm*( flq/ust )**2
    cov(kts) = phm*( flt/ust )*( flq/ust )
!
    DO k = kts+1,kte
       vkz = karman*zw(k)
       el (k) = vkz/( 1.0 + vkz/100.0 )
!       qke(k) = 0.0
!
       tsq(k) = 0.0
       qsq(k) = 0.0
       cov(k) = 0.0
    END DO
!
!   **  Initialization with an iterative manner          **
!   **  lmax is the iteration count. This is arbitrary.  **
    lmax = 5
!
    DO l = 1,lmax
!
!> - call mym_length() to calculate the master length scale.
       CALL mym_length (                          &
            &            kts,kte,xland,           &
            &            dz, dx, zw,              &
            &            rmo, flt, fltv, flq,     &
            &            vt, vq,                  &
            &            u, v, qke,               &
            &            dtv,                     &
            &            el,                      &
            &            zi,theta,                &
            &            qkw,Psig_bl,cldfra_bl1D, &
            &            bl_mynn_mixlength,       &
            &            edmf_w1,edmf_a1          )
!
       DO k = kts+1,kte
          elq = el(k)*qkw(k)
          pdk(k) = elq*( sm(k)*gm(k) + &
               &         sh(k)*gh(k) )
          pdt(k) = elq*  sh(k)*dtl(k)**2
          pdq(k) = elq*  sh(k)*dqw(k)**2
          pdc(k) = elq*  sh(k)*dtl(k)*dqw(k)
       END DO
!
!   **  Strictly, vkz*h(i,j) -> karman*( 0.5*dz(1)*h(i,j)+z0 )  **
       vkz = karman*0.5*dz(kts)
       elv = 0.5*( el(kts+1)+el(kts) ) /  vkz
       IF (INITIALIZE_QKE)THEN
          !qke(kts) = ust**2 * ( b1*pmz*elv    )**(2.0/3.0)
          qke(kts) = 1.0 * MAX(ust,0.02)**2 * ( b1*pmz*elv    )**(2.0/3.0)
       ENDIF

       phm      = phh*b2 / ( b1*pmz/elv**2 )**(1.0/3.0)
       tsq(kts) = phm*( flt/ust )**2
       qsq(kts) = phm*( flq/ust )**2
       cov(kts) = phm*( flt/ust )*( flq/ust )

       DO k = kts+1,kte-1
          b1l = b1*0.25*( el(k+1)+el(k) )
          !tmpq=MAX(b1l*( pdk(k+1)+pdk(k) ),qkemin)
          !add MIN to limit unreasonable QKE
          tmpq=MIN(MAX(b1l*( pdk(k+1)+pdk(k) ),qkemin),125.)
!          PRINT *,'tmpqqqqq',tmpq,pdk(k+1),pdk(k)
          IF (INITIALIZE_QKE)THEN
             qke(k) = tmpq**twothirds
          ENDIF

          IF ( qke(k) .LE. 0.0 ) THEN
             b2l = 0.0
          ELSE
             b2l = b2*( b1l/b1 ) / SQRT( qke(k) )
          END IF

          tsq(k) = b2l*( pdt(k+1)+pdt(k) )
          qsq(k) = b2l*( pdq(k+1)+pdq(k) )
          cov(k) = b2l*( pdc(k+1)+pdc(k) )
       END DO

    END DO

!!    qke(kts)=qke(kts+1)
!!    tsq(kts)=tsq(kts+1)
!!    qsq(kts)=qsq(kts+1)
!!    cov(kts)=cov(kts+1)

    IF (INITIALIZE_QKE)THEN
       qke(kts)=0.5*(qke(kts)+qke(kts+1))
       qke(kte)=qke(kte-1)
    ENDIF
    tsq(kte)=tsq(kte-1)
    qsq(kte)=qsq(kte-1)
    cov(kte)=cov(kte-1)

!
!    RETURN

  END SUBROUTINE mym_initialize
!> @}

!
! ==================================================================
!     SUBROUTINE  mym_level2:
!
!     Input variables:    see subroutine mym_initialize
!
!     Output variables:
!       dtl(nx,nz,ny) : Vertical gradient of Theta_l             (K/m)
!       dqw(nx,nz,ny) : Vertical gradient of Q_w
!       dtv(nx,nz,ny) : Vertical gradient of Theta_V             (K/m)
!       gm (nx,nz,ny) : G_M divided by L^2/q^2                (s^(-2))
!       gh (nx,nz,ny) : G_H divided by L^2/q^2                (s^(-2))
!       sm (nx,nz,ny) : Stability function for momentum, at Level 2
!       sh (nx,nz,ny) : Stability function for heat, at Level 2
!
!       These are defined on the walls of the grid boxes.
!

!>\ingroup gsd_mynn_edmf
!! This subroutine calculates the level 2, non-dimensional wind shear
!! \f$G_M\f$ and vertical temperature gradient \f$G_H\f$ as well as
!! the level 2 stability funcitons \f$S_h\f$ and \f$S_m\f$.
!!\param kts    horizontal dimension
!!\param kte    vertical dimension
!!\param dz     vertical grid spacings (\f$m\f$)
!!\param u      west-east component of the horizontal wind (\f$m s^{-1}\f$)
!!\param v      south-north component of the horizontal wind (\f$m s^{-1}\f$)
!!\param thl    liquid water potential temperature
!!\param qw     total water content \f$Q_w\f$
!!\param ql     liquid water content (\f$kg kg^{-1}\f$)
!!\param vt
!!\param vq
!!\param dtl     vertical gradient of \f$\theta_l\f$ (\f$K m^{-1}\f$)
!!\param dqw     vertical gradient of \f$Q_w\f$
!!\param dtv     vertical gradient of \f$\theta_V\f$ (\f$K m^{-1}\f$)
!!\param gm      \f$G_M\f$ divided by \f$L^{2}/q^{2}\f$ (\f$s^{-2}\f$)
!!\param gh      \f$G_H\f$ divided by \f$L^{2}/q^{2}\f$ (\f$s^{-2}\f$)
!!\param sm      stability function for momentum, at Level 2
!!\param sh      stability function for heat, at Level 2
!!\section gen_mym_level2 GSD MYNN-EDMF mym_level2 General Algorithm
!! @ {
  SUBROUTINE  mym_level2 (kts,kte,                &
       &            dz,                           &
       &            u, v, thl, thetav, qw,        &
       &            ql, vt, vq,                   &
       &            dtl, dqw, dtv, gm, gh, sm, sh )
!
!-------------------------------------------------------------------

    integer, intent(in)   :: kts,kte

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    real(kind_phys), dimension(kts:kte), intent(in)  :: dz
    real(kind_phys), dimension(kts:kte), intent(in)  :: u,v, &
         &thl,qw,ql,vt,vq,thetav
    real(kind_phys), dimension(kts:kte), intent(out) ::      &
         &dtl,dqw,dtv,gm,gh,sm,sh

    integer :: k

    real(kind_phys):: rfc,f1,f2,rf1,rf2,smc,shc,             &
         &ri1,ri2,ri3,ri4,duz,dtz,dqz,vtt,vqq,dtq,dzk,       &
         &afk,abk,ri,rf

    real(kind_phys)::   a2fac

!    ev  = 2.5e6
!    tv0 = 0.61*tref
!    tv1 = 1.61*tref
!    gtr = 9.81/tref
!
    rfc = g1/( g1+g2 )
    f1  = b1*( g1-c1 ) +3.0*a2*( 1.0    -c2 )*( 1.0-c5 ) &
    &                   +2.0*a1*( 3.0-2.0*c2 )
    f2  = b1*( g1+g2 ) -3.0*a1*( 1.0    -c2 )
    rf1 = b1*( g1-c1 )/f1
    rf2 = b1*  g1     /f2
    smc = a1 /a2*  f1/f2
    shc = 3.0*a2*( g1+g2 )
!
    ri1 = 0.5/smc
    ri2 = rf1*smc
    ri3 = 4.0*rf2*smc -2.0*ri2
    ri4 = ri2**2
!
    DO k = kts+1,kte
       dzk = 0.5  *( dz(k)+dz(k-1) )
       afk = dz(k)/( dz(k)+dz(k-1) )
       abk = 1.0 -afk
       duz = ( u(k)-u(k-1) )**2 +( v(k)-v(k-1) )**2
       duz =   duz                    /dzk**2
       dtz = ( thl(k)-thl(k-1) )/( dzk )
       dqz = ( qw(k)-qw(k-1) )/( dzk )
!
       vtt =  1.0 +vt(k)*abk +vt(k-1)*afk  ! Beta-theta in NN09, Eq. 39
       vqq =  tv0 +vq(k)*abk +vq(k-1)*afk  ! Beta-q
       dtq =  vtt*dtz +vqq*dqz
       !Alternatively, use theta-v without the SGS clouds
       !dtq = ( thetav(k)-thetav(k-1) )/( dzk )
!
       dtl(k) =  dtz
       dqw(k) =  dqz
       dtv(k) =  dtq
!?      dtv(i,j,k) =  dtz +tv0*dqz
!?   :              +( xlv/pi0(i,j,k)-tv1 )
!?   :              *( ql(i,j,k)-ql(i,j,k-1) )/( dzk*h(i,j) )
!
       gm (k) =  duz
       gh (k) = -dtq*gtr
!
!   **  Gradient Richardson number  **
       ri = -gh(k)/MAX( duz, 1.0e-10 )

    !a2fac is needed for the Canuto/Kitamura mod
    IF (CKmod .eq. 1) THEN
       a2fac = 1./(1. + MAX(ri,0.0))
    ELSE
       a2fac = 1.
    ENDIF

       rfc = g1/( g1+g2 )
       f1  = b1*( g1-c1 ) +3.0*a2*a2fac *( 1.0    -c2 )*( 1.0-c5 ) &
    &                     +2.0*a1*( 3.0-2.0*c2 )
       f2  = b1*( g1+g2 ) -3.0*a1*( 1.0    -c2 )
       rf1 = b1*( g1-c1 )/f1
       rf2 = b1*  g1     /f2
       smc = a1 /(a2*a2fac)*  f1/f2
       shc = 3.0*(a2*a2fac)*( g1+g2 )

       ri1 = 0.5/smc
       ri2 = rf1*smc
       ri3 = 4.0*rf2*smc -2.0*ri2
       ri4 = ri2**2

!   **  Flux Richardson number  **
       rf = MIN( ri1*( ri + ri2-SQRT(ri**2 - ri3*ri + ri4) ), rfc )
!
       sh (k) = shc*( rfc-rf )/( 1.0-rf )
       sm (k) = smc*( rf1-rf )/( rf2-rf ) * sh(k)
    END DO
!
!    RETURN

#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mym_level2
!! @}

! ==================================================================
!     SUBROUTINE  mym_length:
!
!     Input variables:    see subroutine mym_initialize
!
!     Output variables:   see subroutine mym_initialize
!
!     Work arrays:
!       elt(nx,ny)      : Length scale depending on the PBL depth    (m)
!       vsc(nx,ny)      : Velocity scale q_c                       (m/s)
!                         at first, used for computing elt
!
!     NOTE: the mixing lengths are meant to be calculated at the full-
!           sigmal levels (or interfaces beween the model layers).
!
!>\ingroup gsd_mynn_edmf
!! This subroutine calculates the mixing lengths.
  SUBROUTINE  mym_length (                     &
    &            kts,kte,xland,                &
    &            dz, dx, zw,                   &
    &            rmo, flt, fltv, flq,          &
    &            vt, vq,                       &
    &            u1, v1, qke,                  &
    &            dtv,                          &
    &            el,                           &
    &            zi, theta, qkw,               &
    &            Psig_bl, cldfra_bl1D,         &
    &            bl_mynn_mixlength,            &
    &            edmf_w1,edmf_a1               )

!-------------------------------------------------------------------

    integer, intent(in)   :: kts,kte

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    integer, intent(in)   :: bl_mynn_mixlength
    real(kind_phys), dimension(kts:kte), intent(in)   :: dz
    real(kind_phys), dimension(kts:kte+1), intent(in) :: zw
    real(kind_phys), intent(in) :: rmo,flt,fltv,flq,Psig_bl,xland
    real(kind_phys), intent(in) :: dx,zi
    real(kind_phys), dimension(kts:kte), intent(in)   :: u1,v1,  &
         &qke,vt,vq,cldfra_bl1D,edmf_w1,edmf_a1
    real(kind_phys), dimension(kts:kte), intent(out)  :: qkw, el
    real(kind_phys), dimension(kts:kte), intent(in)   :: dtv
    real(kind_phys):: elt,vsc
    real(kind_phys), dimension(kts:kte), intent(in) :: theta
    real(kind_phys), dimension(kts:kte) :: qtke,elBLmin,elBLavg,thetaw
    real(kind_phys):: wt,wt2,zi2,h1,h2,hs,elBLmin0,elBLavg0,cldavg

    ! THE FOLLOWING CONSTANTS ARE IMPORTANT FOR REGULATING THE
    ! MIXING LENGTHS:
    real(kind_phys):: cns,   &   !< for surface layer (els) in stable conditions
            alp1,            &   !< for turbulent length scale (elt)
            alp2,            &   !< for buoyancy length scale (elb)
            alp3,            &   !< for buoyancy enhancement factor of elb
            alp4,            &   !< for surface layer (els) in unstable conditions
            alp5,            &   !< for BouLac mixing length or above PBLH
            alp6                 !< for mass-flux/

    !THE FOLLOWING LIMITS DO NOT DIRECTLY AFFECT THE ACTUAL PBLH.
    !THEY ONLY IMPOSE LIMITS ON THE CALCULATION OF THE MIXING LENGTH
    !SCALES SO THAT THE BOULAC MIXING LENGTH (IN FREE ATMOS) DOES
    !NOT ENCROACH UPON THE BOUNDARY LAYER MIXING LENGTH (els, elb & elt).
    real(kind_phys), parameter :: minzi = 300.  !< min mixed-layer height
    real(kind_phys), parameter :: maxdz = 750.  !< max (half) transition layer depth
                                     !! =0.3*2500 m PBLH, so the transition
                                     !! layer stops growing for PBLHs > 2.5 km.
    real(kind_phys), parameter :: mindz = 300.  !< 300  !min (half) transition layer depth

    !SURFACE LAYER LENGTH SCALE MODS TO REDUCE IMPACT IN UPPER BOUNDARY LAYER
    real(kind_phys), parameter :: ZSLH = 100. !< Max height correlated to surface conditions (m)
    real(kind_phys), parameter :: CSL = 2.    !< CSL = constant of proportionality to L O(1)


    integer :: i,j,k
    real(kind_phys):: afk,abk,zwk,zwk1,dzk,qdz,vflx,bv,tau_cloud,     &
            & wstar,elb,els,elf,el_stab,el_mf,el_stab_mf,elb_mf,      &
            & PBLH_PLUS_ENT,Uonset,Ugrid,wt_u,el_les
    real(kind_phys), parameter :: ctau = 1000. !constant for tau_cloud

!    tv0 = 0.61*tref
!    gtr = 9.81/tref

    SELECT CASE(bl_mynn_mixlength)

      CASE (0) ! ORIGINAL MYNN MIXING LENGTH + BouLac

        cns  = 2.7
        alp1 = 0.23
        alp2 = 1.0
        alp3 = 5.0
        alp4 = 100.
        alp5 = 0.3

        ! Impose limits on the height integration for elt and the transition layer depth
        zi2  = MIN(10000.,zw(kte-2))  !originally integrated to model top, not just 10 km.
        h1=MAX(0.3*zi2,mindz)
        h1=MIN(h1,maxdz)         ! 1/2 transition layer depth
        h2=h1/2.0                ! 1/4 transition layer depth

        qkw(kts) = SQRT(MAX(qke(kts),1.0e-10))
        DO k = kts+1,kte
           afk = dz(k)/( dz(k)+dz(k-1) )
           abk = 1.0 -afk
           qkw(k) = SQRT(MAX(qke(k)*abk+qke(k-1)*afk,1.0e-3))
        END DO

        elt = 1.0e-5
        vsc = 1.0e-5

        !   **  Strictly, zwk*h(i,j) -> ( zwk*h(i,j)+z0 )  **
        k = kts+1
        zwk = zw(k)
        DO WHILE (zwk .LE. zi2+h1)
           dzk = 0.5*( dz(k)+dz(k-1) )
           qdz = MAX( qkw(k)-qmin, 0.03 )*dzk
           elt = elt +qdz*zwk
           vsc = vsc +qdz
           k   = k+1
           zwk = zw(k)
        END DO

        elt =  alp1*elt/vsc
        vflx = ( vt(kts)+1.0 )*flt +( vq(kts)+tv0 )*flq
        vsc = ( gtr*elt*MAX( vflx, 0.0 ) )**(1.0/3.0)

        !   **  Strictly, el(i,k=1) is not zero.  **
        el(kts) = 0.0
        zwk1    = zw(kts+1)

        DO k = kts+1,kte
           zwk = zw(k)              !full-sigma levels

           !   **  Length scale limited by the buoyancy effect  **
           IF ( dtv(k) .GT. 0.0 ) THEN
              bv  = SQRT( gtr*dtv(k) )
              elb = alp2*qkw(k) / bv &
                  &       *( 1.0 + alp3/alp2*&
                  &SQRT( vsc/( bv*elt ) ) )
              elf = alp2 * qkw(k)/bv

           ELSE
              elb = 1.0e10
              elf = elb
           ENDIF

           !   **  Length scale in the surface layer  **
           IF ( rmo .GT. 0.0 ) THEN
              els  = karman*zwk/(1.0+cns*MIN( zwk*rmo, zmax ))
           ELSE
              els  =  karman*zwk*( 1.0 - alp4* zwk*rmo )**0.2
           END IF

           !   ** HARMONC AVERGING OF MIXING LENGTH SCALES:
           !       el(k) =      MIN(elb/( elb/elt+elb/els+1.0 ),elf)
           !       el(k) =      elb/( elb/elt+elb/els+1.0 )

           wt=.5*TANH((zwk - (zi2+h1))/h2) + .5

           el(k) = MIN(elb/( elb/elt+elb/els+1.0 ),elf)

        END DO

      CASE (1) !NONLOCAL (using BouLac) FORM OF MIXING LENGTH

        ugrid = sqrt(u1(kts)**2 + v1(kts)**2)
        uonset= 15.
        wt_u = (1.0 - min(max(ugrid - uonset, 0.0)/30.0, 0.5))
        cns  = 2.7 !was 3.5
        alp1 = 0.23
        alp2 = 0.3
        alp3 = 2.5 * wt_u !taper off bouyancy enhancement in shear-driven pbls
        alp4 = 5.0
        alp5 = 0.3
        alp6 = 50.

        ! Impose limits on the height integration for elt and the transition layer depth
        zi2=MAX(zi,300.) !minzi)
        h1=MAX(0.3*zi2,300.)
        h1=MIN(h1,600.)          ! 1/2 transition layer depth
        h2=h1/2.0                ! 1/4 transition layer depth

        qtke(kts)=MAX(0.5*qke(kts), 0.01) !tke at full sigma levels
        thetaw(kts)=theta(kts)            !theta at full-sigma levels
        qkw(kts) = SQRT(MAX(qke(kts),1.0e-10))

        DO k = kts+1,kte
           afk = dz(k)/( dz(k)+dz(k-1) )
           abk = 1.0 -afk
           qkw(k) = SQRT(MAX(qke(k)*abk+qke(k-1)*afk,1.0e-3))
           qtke(k) = 0.5*(qkw(k)**2)     ! q -> TKE
           thetaw(k)= theta(k)*abk + theta(k-1)*afk
        END DO

        elt = 1.0e-5
        vsc = 1.0e-5

        !   **  Strictly, zwk*h(i,j) -> ( zwk*h(i,j)+z0 )  **
        k = kts+1
        zwk = zw(k)
        DO WHILE (zwk .LE. zi2+h1)
           dzk = 0.5*( dz(k)+dz(k-1) )
           qdz = min(max( qkw(k)-qmin, 0.03 ), 30.0)*dzk
           elt = elt +qdz*zwk
           vsc = vsc +qdz
           k   = k+1
           zwk = zw(k)
        END DO

        elt = MIN( MAX( alp1*elt/vsc, 10.), 400.)
        !avoid use of buoyancy flux functions which are ill-defined at the surface
        !vflx = ( vt(kts)+1.0 )*flt + ( vq(kts)+tv0 )*flq
        vflx = fltv
        vsc = ( gtr*elt*MAX( vflx, 0.0 ) )**onethird

        !   **  Strictly, el(i,j,1) is not zero.  **
        el(kts) = 0.0
        zwk1    = zw(kts+1)              !full-sigma levels

        ! COMPUTE BouLac mixing length
        CALL boulac_length(kts,kte,zw,dz,qtke,thetaw,elBLmin,elBLavg)

        DO k = kts+1,kte
           zwk = zw(k)              !full-sigma levels

           !   **  Length scale limited by the buoyancy effect  **
           IF ( dtv(k) .GT. 0.0 ) THEN
              bv  = max( sqrt( gtr*dtv(k) ), 0.0001)
              elb = MAX(alp2*qkw(k),                          &
                  &    alp6*edmf_a1(k-1)*edmf_w1(k-1)) / bv   &
                  &  *( 1.0 + alp3*SQRT( vsc/(bv*elt) ) )
              elb = MIN(elb, zwk)
              elf = 1.0 * qkw(k)/bv
              elBLavg(k) = MAX(elBLavg(k), alp6*edmf_a1(k-1)*edmf_w1(k-1)/bv)
           ELSE
              elb = 1.0e10
              elf = elb
           ENDIF

           !   **  Length scale in the surface layer  **
           IF ( rmo .GT. 0.0 ) THEN
              els  = karman*zwk/(1.0+cns*MIN( zwk*rmo, zmax ))
           ELSE
              els  =  karman*zwk*( 1.0 - alp4* zwk*rmo )**0.2
           END IF

           !   ** NOW BLEND THE MIXING LENGTH SCALES:
           wt=.5*TANH((zwk - (zi2+h1))/h2) + .5

           !add blending to use BouLac mixing length in free atmos;
           !defined relative to the PBLH (zi) + transition layer (h1)
           !el(k) = MIN(elb/( elb/elt+elb/els+1.0 ),elf)
           !try squared-blending - but take out elb (makes it underdiffusive)
           !el(k) = SQRT( els**2/(1. + (els**2/elt**2) +(els**2/elb**2)))
           el(k) = sqrt( els**2/(1. + (els**2/elt**2)))
           el(k) = min(el(k), elb)
           el(k) = MIN (el(k), elf)
           el(k) = el(k)*(1.-wt) + alp5*elBLavg(k)*wt

           ! include scale-awareness, except for original MYNN
           el(k) = el(k)*Psig_bl

         END DO

      CASE (2) !Local (mostly) mixing length formulation

        Uonset = 3.5 + dz(kts)*0.1
        Ugrid  = sqrt(u1(kts)**2 + v1(kts)**2)
        cns  = 3.5 !JOE-test  * (1.0 - MIN(MAX(Ugrid - Uonset, 0.0)/10.0, 1.0))
        alp1 = 0.22
        alp2 = 0.30
        alp3 = 2.0
        alp4 = 5.0
        alp5 = alp2 !like alp2, but for free atmosphere
        alp6 = 50.0 !used for MF mixing length

        ! Impose limits on the height integration for elt and the transition layer depth
        !zi2=MAX(zi,minzi)
        zi2=MAX(zi,    300.)
        !h1=MAX(0.3*zi2,mindz)
        !h1=MIN(h1,maxdz)         ! 1/2 transition layer depth
        h1=MAX(0.3*zi2,300.)
        h1=MIN(h1,600.)
        h2=h1*0.5                ! 1/4 transition layer depth

        qtke(kts)=MAX(0.5*qke(kts),0.01) !tke at full sigma levels
        qkw(kts) = SQRT(MAX(qke(kts),1.0e-4))

        DO k = kts+1,kte
           afk = dz(k)/( dz(k)+dz(k-1) )
           abk = 1.0 -afk
           qkw(k) = SQRT(MAX(qke(k)*abk+qke(k-1)*afk,1.0e-3))
           qtke(k) = 0.5*qkw(k)**2  ! qkw -> TKE
        END DO

        elt = 1.0e-5
        vsc = 1.0e-5

        !   **  Strictly, zwk*h(i,j) -> ( zwk*h(i,j)+z0 )  **
        PBLH_PLUS_ENT = MAX(zi+h1, 100.)
        k = kts+1
        zwk = zw(k)
        DO WHILE (zwk .LE. PBLH_PLUS_ENT)
           dzk = 0.5*( dz(k)+dz(k-1) )
           qdz = min(max( qkw(k)-qmin, 0.03 ), 30.0)*dzk
           elt = elt +qdz*zwk
           vsc = vsc +qdz
           k   = k+1
           zwk = zw(k)
        END DO

        elt = MIN( MAX(alp1*elt/vsc, 10.), 400.)
        !avoid use of buoyancy flux functions which are ill-defined at the surface
        !vflx = ( vt(kts)+1.0 )*flt +( vq(kts)+tv0 )*flq
        vflx = fltv
        vsc = ( gtr*elt*MAX( vflx, 0.0 ) )**onethird

        !   **  Strictly, el(i,j,1) is not zero.  **
        el(kts) = 0.0
        zwk1    = zw(kts+1)

        DO k = kts+1,kte
           zwk = zw(k)              !full-sigma levels
           dzk = 0.5*( dz(k)+dz(k-1) )
           cldavg = 0.5*(cldfra_bl1D(k-1)+cldfra_bl1D(k))

           !   **  Length scale limited by the buoyancy effect  **
           IF ( dtv(k) .GT. 0.0 ) THEN
              !impose min value on bv
              bv  = MAX( SQRT( gtr*dtv(k) ), 0.001)
              !elb_mf = alp2*qkw(k) / bv  &
              elb_mf = MAX(alp2*qkw(k),                    &
                  &    alp6*edmf_a1(k-1)*edmf_w1(k-1)) / bv    &
                  &  *( 1.0 + alp3*SQRT( vsc/( bv*elt ) ) )
              elb = MIN(MAX(alp5*qkw(k), alp6*edmf_a1(k)*edmf_w1(k))/bv, zwk)

              !tau_cloud = MIN(MAX(0.5*zi/((gtr*zi*MAX(vflx,1.0e-4))**onethird),30.),150.)
              wstar = 1.25*(gtr*zi*MAX(vflx,1.0e-4))**onethird
              tau_cloud = MIN(MAX(ctau * wstar/grav, 30.), 150.)
              !minimize influence of surface heat flux on tau far away from the PBLH.
              wt=.5*TANH((zwk - (zi2+h1))/h2) + .5
              tau_cloud = tau_cloud*(1.-wt) + 50.*wt
              elf = MIN(MAX(tau_cloud*SQRT(MIN(qtke(k),40.)), &
                  &         alp6*edmf_a1(k)*edmf_w1(k)/bv), zwk)

              !IF (zwk > zi .AND. elf > 400.) THEN
              !   ! COMPUTE BouLac mixing length
              !   !CALL boulac_length0(k,kts,kte,zw,dz,qtke,thetaw,elBLmin0,elBLavg0)
              !   !elf = alp5*elBLavg0
              !   elf = MIN(MAX(50.*SQRT(qtke(k)), 400.), zwk)
              !ENDIF

           ELSE
              ! use version in development for RAP/HRRR 2016
              ! JAYMES-
              ! tau_cloud is an eddy turnover timescale;
              ! see Teixeira and Cheinet (2004), Eq. 1, and
              ! Cheinet and Teixeira (2003), Eq. 7.  The
              ! coefficient 0.5 is tuneable. Expression in
              ! denominator is identical to vsc (a convective
              ! velocity scale), except that elt is relpaced
              ! by zi, and zero is replaced by 1.0e-4 to
              ! prevent division by zero.
              !tau_cloud = MIN(MAX(0.5*zi/((gtr*zi*MAX(vflx,1.0e-4))**onethird),50.),150.)
              wstar = 1.25*(gtr*zi*MAX(vflx,1.0e-4))**onethird
              tau_cloud = MIN(MAX(ctau * wstar/grav, 50.), 200.)
              !minimize influence of surface heat flux on tau far away from the PBLH.
              wt=.5*TANH((zwk - (zi2+h1))/h2) + .5
              !tau_cloud = tau_cloud*(1.-wt) + 50.*wt
              tau_cloud = tau_cloud*(1.-wt) + MAX(100.,dzk*0.25)*wt

              elb = MIN(tau_cloud*SQRT(MIN(qtke(k),40.)), zwk)
              !elf = elb
              elf = elb !/(1. + (elb/800.))  !bound free-atmos mixing length to < 800 m.
              elb_mf = elb
         END IF
         elf    = elf/(1. + (elf/800.))  !bound free-atmos mixing length to < 800 m.
         elb_mf = MAX(elb_mf, 0.01) !to avoid divide-by-zero below

         !   **  Length scale in the surface layer  **
         IF ( rmo .GT. 0.0 ) THEN
            els  = karman*zwk/(1.0+cns*MIN( zwk*rmo, zmax ))
         ELSE
            els  =  karman*zwk*( 1.0 - alp4* zwk*rmo )**0.2
         END IF

         !   ** NOW BLEND THE MIXING LENGTH SCALES:
         wt=.5*TANH((zwk - (zi2+h1))/h2) + .5

         !try squared-blending
         el(k) = SQRT( els**2/(1. + (els**2/elt**2) +(els**2/elb_mf**2)))
         el(k) = el(k)*(1.-wt) + elf*wt

         ! include scale-awareness. For now, use simple asymptotic kz -> 12 m (should be ~dz).
         el_les= MIN(els/(1. + (els/12.)), elb_mf)
         el(k) = el(k)*Psig_bl + (1.-Psig_bl)*el_les

       END DO

    END SELECT


#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mym_length


! ==================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine was taken from the BouLac scheme in WRF-ARW
!! and modified for integration into the MYNN PBL scheme.
!! WHILE loops were added to reduce the computational expense.
!! This subroutine computes the length scales up and down
!! and then computes the min, average of the up/down
!! length scales, and also considers the distance to the
!! surface.
  SUBROUTINE boulac_length(kts,kte,zw,dz,qtke,theta,lb1,lb2)
!      dlu = the distance a parcel can be lifted upwards give a finite
!            amount of TKE.
!      dld = the distance a parcel can be displaced downwards given a
!            finite amount of TKE.
!      lb1 = the minimum of the length up and length down
!      lb2 = the average of the length up and length down
!-------------------------------------------------------------------

     integer, intent(in) :: kts,kte
     real(kind_phys), dimension(kts:kte),   intent(in) :: qtke,dz,theta
     real(kind_phys), dimension(kts:kte),   intent(out):: lb1,lb2
     real(kind_phys), dimension(kts:kte+1), intent(in) :: zw

     !LOCAL VARS
     integer :: iz, izz, found
     real(kind_phys), dimension(kts:kte) :: dlu,dld
     real(kind_phys), parameter :: Lmax=2000.  !soft limit
     real(kind_phys):: dzt, zup, beta, zup_inf, bbb, tl, zdo, zdo_sup, zzz

     !print*,"IN MYNN-BouLac",kts, kte

     do iz=kts,kte

        !----------------------------------
        ! FIND DISTANCE UPWARD
        !----------------------------------
        zup=0.
        dlu(iz)=zw(kte+1)-zw(iz)-dz(iz)*0.5
        zzz=0.
        zup_inf=0.
        beta=gtr           !Buoyancy coefficient (g/tref)

        !print*,"FINDING Dup, k=",iz," zw=",zw(iz)

        if (iz .lt. kte) then      !cant integrate upwards from highest level

          found = 0
          izz=iz
          DO WHILE (found .EQ. 0)

            if (izz .lt. kte) then
              dzt=dz(izz)                    ! layer depth above
              zup=zup-beta*theta(iz)*dzt     ! initial PE the parcel has at iz
              !print*,"  ",iz,izz,theta(izz),dz(izz)
              zup=zup+beta*(theta(izz+1)+theta(izz))*dzt*0.5 ! PE gained by lifting a parcel to izz+1
              zzz=zzz+dzt                   ! depth of layer iz to izz+1
              !print*,"  PE=",zup," TKE=",qtke(iz)," z=",zw(izz)
              if (qtke(iz).lt.zup .and. qtke(iz).ge.zup_inf) then
                 bbb=(theta(izz+1)-theta(izz))/dzt
                 if (bbb .ne. 0.) then
                    !fractional distance up into the layer where TKE becomes < PE
                    tl=(-beta*(theta(izz)-theta(iz)) + &
                      & sqrt( max(0.,(beta*(theta(izz)-theta(iz)))**2 + &
                      &       2.*bbb*beta*(qtke(iz)-zup_inf))))/bbb/beta
                 else
                    if (theta(izz) .ne. theta(iz))then
                       tl=(qtke(iz)-zup_inf)/(beta*(theta(izz)-theta(iz)))
                    else
                       tl=0.
                    endif
                 endif
                 dlu(iz)=zzz-dzt+tl
                 !print*,"  FOUND Dup:",dlu(iz)," z=",zw(izz)," tl=",tl
                 found =1
              endif
              zup_inf=zup
              izz=izz+1
             ELSE
              found = 1
            ENDIF

          ENDDO

        endif

        !----------------------------------
        ! FIND DISTANCE DOWN
        !----------------------------------
        zdo=0.
        zdo_sup=0.
        dld(iz)=zw(iz)
        zzz=0.

        !print*,"FINDING Ddown, k=",iz," zwk=",zw(iz)
        if (iz .gt. kts) then  !cant integrate downwards from lowest level

          found = 0
          izz=iz
          DO WHILE (found .EQ. 0)

            if (izz .gt. kts) then
              dzt=dz(izz-1)
              zdo=zdo+beta*theta(iz)*dzt
              !print*,"  ",iz,izz,theta(izz),dz(izz-1)
              zdo=zdo-beta*(theta(izz-1)+theta(izz))*dzt*0.5
              zzz=zzz+dzt
              !print*,"  PE=",zdo," TKE=",qtke(iz)," z=",zw(izz)
              if (qtke(iz).lt.zdo .and. qtke(iz).ge.zdo_sup) then
                 bbb=(theta(izz)-theta(izz-1))/dzt
                 if (bbb .ne. 0.) then
                    tl=(beta*(theta(izz)-theta(iz))+ &
                      & sqrt( max(0.,(beta*(theta(izz)-theta(iz)))**2 + &
                      &       2.*bbb*beta*(qtke(iz)-zdo_sup))))/bbb/beta
                 else
                    if (theta(izz) .ne. theta(iz)) then
                       tl=(qtke(iz)-zdo_sup)/(beta*(theta(izz)-theta(iz)))
                    else
                       tl=0.
                    endif
                 endif
                 dld(iz)=zzz-dzt+tl
                 !print*,"  FOUND Ddown:",dld(iz)," z=",zw(izz)," tl=",tl
                 found = 1
              endif
              zdo_sup=zdo
              izz=izz-1
            ELSE
              found = 1
            ENDIF
          ENDDO

        endif

        !----------------------------------
        ! GET MINIMUM (OR AVERAGE)
        !----------------------------------
        !The surface layer length scale can exceed z for large z/L,
        !so keep maximum distance down > z.
        dld(iz) = min(dld(iz),zw(iz+1))!not used in PBL anyway, only free atmos
        lb1(iz) = min(dlu(iz),dld(iz))     !minimum
        !JOE-fight floating point errors
        dlu(iz)=MAX(0.1,MIN(dlu(iz),1000.))
        dld(iz)=MAX(0.1,MIN(dld(iz),1000.))
        lb2(iz) = sqrt(dlu(iz)*dld(iz))    !average - biased towards smallest
        !lb2(iz) = 0.5*(dlu(iz)+dld(iz))   !average

        !Apply soft limit (only impacts very large lb; lb=100 by 5%, lb=500 by 20%).
        lb1(iz) = lb1(iz)/(1. + (lb1(iz)/Lmax))
        lb2(iz) = lb2(iz)/(1. + (lb2(iz)/Lmax))

        if (iz .eq. kte) then
           lb1(kte) = lb1(kte-1)
           lb2(kte) = lb2(kte-1)
        endif
        !print*,"IN MYNN-BouLac",kts, kte,lb1(iz)
        !print*,"IN MYNN-BouLac",iz,dld(iz),dlu(iz)

     ENDDO

  END SUBROUTINE boulac_length
!
! ==================================================================
!     SUBROUTINE  mym_turbulence:
!
!     Input variables:    see subroutine mym_initialize
!       closure        : closure level (2.5, 2.6, or 3.0)
!
!     # ql, vt, vq, qke, tsq, qsq and cov are changed to input variables.
!
!     Output variables:   see subroutine mym_initialize
!       dfm(nx,nz,ny) : Diffusivity coefficient for momentum,
!                         divided by dz (not dz*h(i,j))            (m/s)
!       dfh(nx,nz,ny) : Diffusivity coefficient for heat,
!                         divided by dz (not dz*h(i,j))            (m/s)
!       dfq(nx,nz,ny) : Diffusivity coefficient for q^2,
!                         divided by dz (not dz*h(i,j))            (m/s)
!       tcd(nx,nz,ny)   : Countergradient diffusion term for Theta_l
!                                                                  (K/s)
!       qcd(nx,nz,ny)   : Countergradient diffusion term for Q_w
!                                                              (kg/kg s)
!       pd?(nx,nz,ny) : Half of the production terms
!
!       Only tcd and qcd are defined at the center of the grid boxes
!
!     # DO NOT forget that tcd and qcd are added on the right-hand side
!       of the equations for Theta_l and Q_w, respectively.
!
!     Work arrays:        see subroutine mym_initialize and level2
!
!     # dtl, dqw, dtv, gm and gh are allowed to share storage units with
!       dfm, dfh, dfq, tcd and qcd, respectively, for saving memory.
!
!>\ingroup gsd_mynn_edmf
!! This subroutine calculates the vertical diffusivity coefficients and the
!! production terms for the turbulent quantities.
!>\section gen_mym_turbulence GSD mym_turbulence General Algorithm
!! Two subroutines mym_level2() and mym_length() are called within this
!!subrouine to collect variable to carry out successive calculations:
!! - mym_level2() calculates the level 2 nondimensional wind shear \f$G_M\f$
!! and vertical temperature gradient \f$G_H\f$ as well as the level 2 stability
!! functions \f$S_h\f$ and \f$S_m\f$.
!! - mym_length() calculates the mixing lengths.
!! - The stability criteria from Helfand and Labraga (1989) are applied.
!! - The stability functions for level 2.5 or level 3.0 are calculated.
!! - If level 3.0 is used, counter-gradient terms are calculated.
!! - Production terms of TKE,\f$\theta^{'2}\f$,\f$q^{'2}\f$, and \f$\theta^{'}q^{'}\f$
!! are calculated.
!! - Eddy diffusivity \f$K_h\f$ and eddy viscosity \f$K_m\f$ are calculated.
!! - TKE budget terms are calculated (if the namelist parameter \p tke_budget
!! is set to True)
  SUBROUTINE  mym_turbulence (                                &
    &            kts,kte,                                     &
    &            xland,closure,                               &
    &            dz, dx, zw,                                  &
    &            u, v, thl, thetav, ql, qw,                   &
    &            qke, tsq, qsq, cov,                          &
    &            vt, vq,                                      &
    &            rmo, flt, fltv, flq,                         &
    &            zi,theta,                                    &
    &            sh, sm,                                      &
    &            El,                                          &
    &            Dfm, Dfh, Dfq, Tcd, Qcd, Pdk, Pdt, Pdq, Pdc, &
    &		 qWT1D,qSHEAR1D,qBUOY1D,qDISS1D,              &
    &            tke_budget,                                  &
    &            Psig_bl,Psig_shcu,cldfra_bl1D,               &
    &            bl_mynn_mixlength,                           &
    &            edmf_w1,edmf_a1,                             &
    &            TKEprodTD,                                   &
    &            spp_pbl,rstoch_col,debug_code                )

!-------------------------------------------------------------------

    integer, intent(in)   :: kts,kte
    logical, intent(in)   :: debug_code

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    integer, intent(in)               :: bl_mynn_mixlength,tke_budget
    real(kind_phys), intent(in)       :: closure
    real(kind_phys), dimension(kts:kte),   intent(in) :: dz
    real(kind_phys), dimension(kts:kte+1), intent(in) :: zw
    real(kind_phys), intent(in)       :: rmo,flt,fltv,flq,                 &
         &Psig_bl,Psig_shcu,xland,dx,zi
    real(kind_phys), dimension(kts:kte), intent(in) :: u,v,thl,thetav,qw,  &
         &ql,vt,vq,qke,tsq,qsq,cov,cldfra_bl1D,edmf_w1,edmf_a1,            &
         &TKEprodTD

    real(kind_phys), dimension(kts:kte), intent(out) :: dfm,dfh,dfq,       &
         &pdk,pdt,pdq,pdc,tcd,qcd,el

    real(kind_phys), dimension(kts:kte), intent(inout) ::                  &
         qWT1D,qSHEAR1D,qBUOY1D,qDISS1D
    real(kind_phys):: q3sq_old,dlsq1,qWTP_old,qWTP_new
    real(kind_phys):: dudz,dvdz,dTdz,upwp,vpwp,Tpwp

    real(kind_phys), dimension(kts:kte) :: qkw,dtl,dqw,dtv,gm,gh,sm,sh

    integer :: k
!    real(kind_phys):: cc2,cc3,e1c,e2c,e3c,e4c,e5c
    real(kind_phys):: e6c,dzk,afk,abk,vtt,vqq,                             &
         &cw25,clow,cupp,gamt,gamq,smd,gamv,elq,elh

    real(kind_phys):: cldavg
    real(kind_phys), dimension(kts:kte), intent(in) :: theta

    real(kind_phys)::  a2fac, duz, ri !JOE-Canuto/Kitamura mod

    real:: auh,aum,adh,adm,aeh,aem,Req,Rsl,Rsl2,                           &
           gmelq,sm20,sh20,sm25max,sh25max,sm25min,sh25min,                &
           sm_pbl,sh_pbl,zi2,wt,slht,wtpr

    DOUBLE PRECISION  q2sq, t2sq, r2sq, c2sq, elsq, gmel, ghel
    DOUBLE PRECISION  q3sq, t3sq, r3sq, c3sq, dlsq, qdiv
    DOUBLE PRECISION  e1, e2, e3, e4, enum, eden, wden

!   Stochastic
    integer,         intent(in)                   ::    spp_pbl
    real(kind_phys), dimension(kts:kte)           ::    rstoch_col
    real(kind_phys):: Prnum, shb
    real(kind_phys), parameter :: Prlimit = 5.0

!
!    tv0 = 0.61*tref
!    gtr = 9.81/tref
!
!    cc2 =  1.0-c2
!    cc3 =  1.0-c3
!    e1c =  3.0*a2*b2*cc3
!    e2c =  9.0*a1*a2*cc2
!    e3c =  9.0*a2*a2*cc2*( 1.0-c5 )
!    e4c = 12.0*a1*a2*cc2
!    e5c =  6.0*a1*a1
!

    CALL mym_level2 (kts,kte,                   &
    &            dz,                            &
    &            u, v, thl, thetav, qw,         &
    &            ql, vt, vq,                    &
    &            dtl, dqw, dtv, gm, gh, sm, sh  )
!
    CALL mym_length (                           &
    &            kts,kte,xland,                 &
    &            dz, dx, zw,                    &
    &            rmo, flt, fltv, flq,           &
    &            vt, vq,                        &
    &            u, v, qke,                     &
    &            dtv,                           &
    &            el,                            &
    &            zi,theta,                      &
    &            qkw,Psig_bl,cldfra_bl1D,       &
    &            bl_mynn_mixlength,             &
    &            edmf_w1,edmf_a1                )
!

    DO k = kts+1,kte
       dzk = 0.5  *( dz(k)+dz(k-1) )
       afk = dz(k)/( dz(k)+dz(k-1) )
       abk = 1.0 -afk
       elsq = el (k)**2
       q3sq = qkw(k)**2
       q2sq = b1*elsq*( sm(k)*gm(k)+sh(k)*gh(k) )

       sh20 = MAX(sh(k), 1e-5)
       sm20 = MAX(sm(k), 1e-5)
       sh(k)= MAX(sh(k), 1e-5)

       !Canuto/Kitamura mod
       duz = ( u(k)-u(k-1) )**2 +( v(k)-v(k-1) )**2
       duz =   duz                    /dzk**2
       !   **  Gradient Richardson number  **
       ri = -gh(k)/MAX( duz, 1.0e-10 )
       IF (CKmod .eq. 1) THEN
          a2fac = 1./(1. + MAX(ri,0.0))
       ELSE
          a2fac = 1.
       ENDIF
       !end Canuto/Kitamura mod

       !level 2.0 Prandtl number
       !Prnum = MIN(sm20/sh20, 4.0)
       !The form of Zilitinkevich et al. (2006) but modified
       !half-way towards Esau and Grachev (2007, Wind Eng)
       !Prnum = MIN(0.76 + 3.0*MAX(ri,0.0), Prlimit)
       Prnum = MIN(0.76 + 4.0*MAX(ri,0.0), Prlimit)
       !Prnum = MIN(0.76 + 5.0*MAX(ri,0.0), Prlimit)
!
!  Modified: Dec/22/2005, from here, (dlsq -> elsq)
       gmel = gm (k)*elsq
       ghel = gh (k)*elsq
!  Modified: Dec/22/2005, up to here

       ! Level 2.0 debug prints
       IF ( debug_code ) THEN
         IF (sh(k)<0.0 .OR. sm(k)<0.0) THEN
           print*,"MYNN; mym_turbulence 2.0; sh=",sh(k)," k=",k
           print*," gm=",gm(k)," gh=",gh(k)," sm=",sm(k)
           print*," q2sq=",q2sq," q3sq=",q3sq," q3/q2=",q3sq/q2sq
           print*," qke=",qke(k)," el=",el(k)," ri=",ri
           print*," PBLH=",zi," u=",u(k)," v=",v(k)
         ENDIF
       ENDIF

!     **  Since qkw is set to more than 0.0, q3sq > 0.0.  **

!     new stability criteria in level 2.5 (as well as level 3) - little/no impact
!     **  Limitation on q, instead of L/q  **
       dlsq =  elsq
       IF ( q3sq/dlsq .LT. -gh(k) ) q3sq = -dlsq*gh(k)

       IF ( q3sq .LT. q2sq ) THEN
          !Apply Helfand & Labraga mod
          qdiv = SQRT( q3sq/q2sq )   !HL89: (1-alfa)
!
          !Use level 2.5 stability functions
          !e1   = q3sq - e1c*ghel*a2fac
          !e2   = q3sq - e2c*ghel*a2fac
          !e3   = e1   + e3c*ghel*a2fac**2
          !e4   = e1   - e4c*ghel*a2fac
          !eden = e2*e4 + e3*e5c*gmel
          !eden = MAX( eden, 1.0d-20 )
          !sm(k) = q3sq*a1*( e3-3.0*c1*e4       )/eden
          !!JOE-Canuto/Kitamura mod
          !!sh(k) = q3sq*a2*( e2+3.0*c1*e5c*gmel )/eden
          !sh(k) = q3sq*(a2*a2fac)*( e2+3.0*c1*e5c*gmel )/eden
          !sm(k) = Prnum*sh(k)
          !sm(k) = sm(k) * qdiv

          !Use level 2.0 functions as in original MYNN
          sh(k) = sh(k) * qdiv
          sm(k) = sm(k) * qdiv
        !  !sm_pbl = sm(k) * qdiv
        !
        !  !Or, use the simple Pr relationship
        !  sm(k) = Prnum*sh(k)
        !
        !  !or blend them:
        !  zi2   = MAX(zi, 300.)
        !  wt    =.5*TANH((zw(k) - zi2)/200.) + .5
        !  sm(k) = sm_pbl*(1.-wt) + sm(k)*wt

          !Recalculate terms for later use
          !JOE-Canuto/Kitamura mod
          !e1   = q3sq - e1c*ghel * qdiv**2
          !e2   = q3sq - e2c*ghel * qdiv**2
          !e3   = e1   + e3c*ghel * qdiv**2
          !e4   = e1   - e4c*ghel * qdiv**2
          e1   = q3sq - e1c*ghel*a2fac * qdiv**2
          e2   = q3sq - e2c*ghel*a2fac * qdiv**2
          e3   = e1   + e3c*ghel*a2fac**2 * qdiv**2
          e4   = e1   - e4c*ghel*a2fac * qdiv**2
          eden = e2*e4 + e3*e5c*gmel * qdiv**2
          eden = MAX( eden, 1.0d-20 )
          !!JOE-Canuto/Kitamura mod
          !!sh(k) = q3sq*a2*( e2+3.0*c1*e5c*gmel )/eden  - retro 5
          !sh(k) = q3sq*(a2*a2fac)*( e2+3.0*c1*e5c*gmel )/eden
          !sm(k) = Prnum*sh(k)
       ELSE
          !JOE-Canuto/Kitamura mod
          !e1   = q3sq - e1c*ghel
          !e2   = q3sq - e2c*ghel
          !e3   = e1   + e3c*ghel
          !e4   = e1   - e4c*ghel
          e1   = q3sq - e1c*ghel*a2fac
          e2   = q3sq - e2c*ghel*a2fac
          e3   = e1   + e3c*ghel*a2fac**2
          e4   = e1   - e4c*ghel*a2fac
          eden = e2*e4 + e3*e5c*gmel
          eden = MAX( eden, 1.0d-20 )

          qdiv = 1.0
          !Use level 2.5 stability functions
          sm(k) = q3sq*a1*( e3-3.0*c1*e4       )/eden
        !  sm_pbl = q3sq*a1*( e3-3.0*c1*e4       )/eden
          !!JOE-Canuto/Kitamura mod
          !!sh(k) = q3sq*a2*( e2+3.0*c1*e5c*gmel )/eden
          sh(k) = q3sq*(a2*a2fac)*( e2+3.0*c1*e5c*gmel )/eden
        !  sm(k) = Prnum*sh(k)

        !  !or blend them:
        !  zi2   = MAX(zi, 300.)
        !  wt    = .5*TANH((zw(k) - zi2)/200.) + .5
        !  sm(k) = sm_pbl*(1.-wt) + sm(k)*wt
       END IF !end Helfand & Labraga check

       !Impose broad limits on Sh and Sm:
       gmelq    = MAX(gmel/q3sq, 1e-8)
       sm25max  = 4.  !MIN(sm20*3.0, SQRT(.1936/gmelq))
       sh25max  = 4.  !MIN(sh20*3.0, 0.76*b2)
       sm25min  = 0.0 !MAX(sm20*0.1, 1e-6)
       sh25min  = 0.0 !MAX(sh20*0.1, 1e-6)

       !JOE: Level 2.5 debug prints
       ! HL88 , lev2.5 criteria from eqs. 3.17, 3.19, & 3.20
       IF ( debug_code ) THEN
         IF ((sh(k)<sh25min .OR. sm(k)<sm25min .OR. &
              sh(k)>sh25max .OR. sm(k)>sm25max) ) THEN
           print*,"In mym_turbulence 2.5: k=",k
           print*," sm=",sm(k)," sh=",sh(k)
           print*," ri=",ri," Pr=",sm(k)/MAX(sh(k),1e-8)
           print*," gm=",gm(k)," gh=",gh(k)
           print*," q2sq=",q2sq," q3sq=",q3sq, q3sq/q2sq
           print*," qke=",qke(k)," el=",el(k)
           print*," PBLH=",zi," u=",u(k)," v=",v(k)
           print*," SMnum=",q3sq*a1*( e3-3.0*c1*e4)," SMdenom=",eden
           print*," SHnum=",q3sq*(a2*a2fac)*( e2+3.0*c1*e5c*gmel ),&
                  " SHdenom=",eden
         ENDIF
       ENDIF

       !Enforce constraints for level 2.5 functions
       IF ( sh(k) > sh25max ) sh(k) = sh25max
       IF ( sh(k) < sh25min ) sh(k) = sh25min
       !IF ( sm(k) > sm25max ) sm(k) = sm25max
       !IF ( sm(k) < sm25min ) sm(k) = sm25min
       !sm(k) = Prnum*sh(k)

       !surface layer PR
       !slht  = zi*0.1
       !wtpr  = min( max( (slht - zw(k))/slht, 0.0), 1.0) ! 1 at z=0, 0 above sfc layer
       !Prlim = 1.0*wtpr + (1.0 - wtpr)*Prlimit
       !Prlim = 2.0*wtpr + (1.0 - wtpr)*Prlimit
       !sm(k) = MIN(sm(k), Prlim*Sh(k))
       !Pending more testing, keep same Pr limit in sfc layer
       shb   = max(sh(k), 0.002)
       sm(k) = MIN(sm(k), Prlimit*shb)

!   **  Level 3 : start  **
       IF ( closure .GE. 3.0 ) THEN
          t2sq = qdiv*b2*elsq*sh(k)*dtl(k)**2
          r2sq = qdiv*b2*elsq*sh(k)*dqw(k)**2
          c2sq = qdiv*b2*elsq*sh(k)*dtl(k)*dqw(k)
          t3sq = MAX( tsq(k)*abk+tsq(k-1)*afk, 0.0 )
          r3sq = MAX( qsq(k)*abk+qsq(k-1)*afk, 0.0 )
          c3sq =      cov(k)*abk+cov(k-1)*afk

!  Modified: Dec/22/2005, from here
          c3sq = SIGN( MIN( ABS(c3sq), SQRT(t3sq*r3sq) ), c3sq )
!
          vtt  = 1.0 +vt(k)*abk +vt(k-1)*afk
          vqq  = tv0 +vq(k)*abk +vq(k-1)*afk

          t2sq = vtt*t2sq +vqq*c2sq
          r2sq = vtt*c2sq +vqq*r2sq
          c2sq = MAX( vtt*t2sq+vqq*r2sq, 0.0d0 )
          t3sq = vtt*t3sq +vqq*c3sq
          r3sq = vtt*c3sq +vqq*r3sq
          c3sq = MAX( vtt*t3sq+vqq*r3sq, 0.0d0 )
!
          cw25 = e1*( e2 + 3.0*c1*e5c*gmel*qdiv**2 )/( 3.0*eden )
!
!     **  Limitation on q, instead of L/q  **
          dlsq =  elsq
          IF ( q3sq/dlsq .LT. -gh(k) ) q3sq = -dlsq*gh(k)
!
!     **  Limitation on c3sq (0.12 =< cw =< 0.76) **
          ! Use Janjic's (2001; p 13-17) methodology (eqs 4.11-414 and 5.7-5.10)
          ! to calculate an exact limit for c3sq:
          auh = 27.*a1*((a2*a2fac)**2)*b2*(gtr)**2
          aum = 54.*(a1**2)*(a2*a2fac)*b2*c1*(gtr)
          adh = 9.*a1*((a2*a2fac)**2)*(12.*a1 + 3.*b2)*(gtr)**2
          adm = 18.*(a1**2)*(a2*a2fac)*(b2 - 3.*(a2*a2fac))*(gtr)

          aeh = (9.*a1*((a2*a2fac)**2)*b1 +9.*a1*((a2*a2fac)**2)* &
                (12.*a1 + 3.*b2))*(gtr)
          aem = 3.*a1*(a2*a2fac)*b1*(3.*(a2*a2fac) + 3.*b2*c1 + &
                (18.*a1*c1 - b2)) + &
                (18.)*(a1**2)*(a2*a2fac)*(b2 - 3.*(a2*a2fac))

          Req = -aeh/aem
          Rsl = (auh + aum*Req)/(3.*adh + 3.*adm*Req)
          !For now, use default values, since tests showed little/no sensitivity
          Rsl = .12             !lower limit
          Rsl2= 1.0 - 2.*Rsl    !upper limit
          !IF (k==2)print*,"Dynamic limit RSL=",Rsl
          !IF (Rsl < 0.10 .OR. Rsl > 0.18) THEN
          !   print*,'--- ERROR: MYNN: Dynamic Cw '// &
          !        'limit exceeds reasonable limits'
          !   print*," MYNN: Dynamic Cw limit needs attention=",Rsl
          !ENDIF

          !JOE-Canuto/Kitamura mod
          !e2   = q3sq - e2c*ghel * qdiv**2
          !e3   = q3sq + e3c*ghel * qdiv**2
          !e4   = q3sq - e4c*ghel * qdiv**2
          e2   = q3sq - e2c*ghel*a2fac * qdiv**2
          e3   = q3sq + e3c*ghel*a2fac**2 * qdiv**2
          e4   = q3sq - e4c*ghel*a2fac * qdiv**2
          eden = e2*e4  + e3 *e5c*gmel * qdiv**2

          !JOE-Canuto/Kitamura mod
          !wden = cc3*gtr**2 * dlsq**2/elsq * qdiv**2 &
          !     &        *( e2*e4c - e3c*e5c*gmel * qdiv**2 )
          wden = cc3*gtr**2 * dlsq**2/elsq * qdiv**2 &
               &        *( e2*e4c*a2fac - e3c*e5c*gmel*a2fac**2 * qdiv**2 )

          IF ( wden .NE. 0.0 ) THEN
             !JOE: test dynamic limits
             clow = q3sq*( 0.12-cw25 )*eden/wden
             cupp = q3sq*( 0.76-cw25 )*eden/wden
             !clow = q3sq*( Rsl -cw25 )*eden/wden
             !cupp = q3sq*( Rsl2-cw25 )*eden/wden
!
             IF ( wden .GT. 0.0 ) THEN
                c3sq  = MIN( MAX( c3sq, c2sq+clow ), c2sq+cupp )
             ELSE
                c3sq  = MAX( MIN( c3sq, c2sq+clow ), c2sq+cupp )
             END IF
          END IF
!
          e1   = e2 + e5c*gmel * qdiv**2
          eden = MAX( eden, 1.0d-20 )
!  Modified: Dec/22/2005, up to here

          !JOE-Canuto/Kitamura mod
          !e6c  = 3.0*a2*cc3*gtr * dlsq/elsq
          e6c  = 3.0*(a2*a2fac)*cc3*gtr * dlsq/elsq

          !============================
          !     **  for Gamma_theta  **
          !!          enum = qdiv*e6c*( t3sq-t2sq )
          IF ( t2sq .GE. 0.0 ) THEN
             enum = MAX( qdiv*e6c*( t3sq-t2sq ), 0.0d0 )
          ELSE
             enum = MIN( qdiv*e6c*( t3sq-t2sq ), 0.0d0 )
          ENDIF
          gamt =-e1  *enum    /eden

          !============================
          !     **  for Gamma_q  **
          !!          enum = qdiv*e6c*( r3sq-r2sq )
          IF ( r2sq .GE. 0.0 ) THEN
             enum = MAX( qdiv*e6c*( r3sq-r2sq ), 0.0d0 )
          ELSE
             enum = MIN( qdiv*e6c*( r3sq-r2sq ), 0.0d0 )
          ENDIF
          gamq =-e1  *enum    /eden

          !============================
          !     **  for Sm' and Sh'd(Theta_V)/dz  **
          !!          enum = qdiv*e6c*( c3sq-c2sq )
          enum = MAX( qdiv*e6c*( c3sq-c2sq ), 0.0d0)

          !JOE-Canuto/Kitamura mod
          !smd  = dlsq*enum*gtr/eden * qdiv**2 * (e3c+e4c)*a1/a2
          smd  = dlsq*enum*gtr/eden * qdiv**2 * (e3c*a2fac**2 + &
               & e4c*a2fac)*a1/(a2*a2fac)

          gamv = e1  *enum*gtr/eden
          sm(k) = sm(k) +smd

          !============================
          !     **  For elh (see below), qdiv at Level 3 is reset to 1.0.  **
          qdiv = 1.0

          ! Level 3 debug prints
          IF ( debug_code ) THEN
            IF (sh(k)<-0.3 .OR. sm(k)<-0.3 .OR. &
              qke(k) < -0.1 .or. ABS(smd) .gt. 2.0) THEN
              print*," MYNN; mym_turbulence3.0; sh=",sh(k)," k=",k
              print*," gm=",gm(k)," gh=",gh(k)," sm=",sm(k)
              print*," q2sq=",q2sq," q3sq=",q3sq," q3/q2=",q3sq/q2sq
              print*," qke=",qke(k)," el=",el(k)," ri=",ri
              print*," PBLH=",zi," u=",u(k)," v=",v(k)
            ENDIF
          ENDIF

!   **  Level 3 : end  **

       ELSE
!     **  At Level 2.5, qdiv is not reset.  **
          gamt = 0.0
          gamq = 0.0
          gamv = 0.0
       END IF
!
!      Add min background stability function (diffusivity) within model levels
!      with active plumes and clouds.
       cldavg = 0.5*(cldfra_bl1D(k-1) + cldfra_bl1D(k))
       IF (edmf_a1(k) > 0.001 .OR. cldavg > 0.02) THEN
           ! for mass-flux columns
           sm(k) = MAX(sm(k), 0.03*MIN(10.*edmf_a1(k)*edmf_w1(k),1.0) )
           sh(k) = MAX(sh(k), 0.03*MIN(10.*edmf_a1(k)*edmf_w1(k),1.0) )
           ! for clouds
           sm(k) = MAX(sm(k), 0.05*MIN(cldavg,1.0) )
           sh(k) = MAX(sh(k), 0.05*MIN(cldavg,1.0) )
       ENDIF
!
       elq = el(k)*qkw(k)
       elh = elq*qdiv

       ! Production of TKE (pdk), T-variance (pdt),
       ! q-variance (pdq), and covariance (pdc)
       pdk(k) = elq*( sm(k)*gm(k)                &
            &        +sh(k)*gh(k)+gamv ) +       &
            &    0.5*TKEprodTD(k)        ! xmchen
       pdt(k) = elh*( sh(k)*dtl(k)+gamt )*dtl(k)
       pdq(k) = elh*( sh(k)*dqw(k)+gamq )*dqw(k)
       pdc(k) = elh*( sh(k)*dtl(k)+gamt )        &
            &   *dqw(k)*0.5                      &
            & + elh*( sh(k)*dqw(k)+gamq )*dtl(k)*0.5

       ! Contergradient terms
       tcd(k) = elq*gamt
       qcd(k) = elq*gamq

       ! Eddy Diffusivity/Viscosity divided by dz
       dfm(k) = elq*sm(k) / dzk
       dfh(k) = elq*sh(k) / dzk
!  Modified: Dec/22/2005, from here
!   **  In sub.mym_predict, dfq for the TKE and scalar variance **
!   **  are set to 3.0*dfm and 1.0*dfm, respectively. (Sqfac)   **
       dfq(k) =     dfm(k)
!  Modified: Dec/22/2005, up to here

   IF (tke_budget .eq. 1) THEN
       !TKE BUDGET
!       dudz = ( u(k)-u(k-1) )/dzk
!       dvdz = ( v(k)-v(k-1) )/dzk
!       dTdz = ( thl(k)-thl(k-1) )/dzk

!       upwp = -elq*sm(k)*dudz
!       vpwp = -elq*sm(k)*dvdz
!       Tpwp = -elq*sh(k)*dTdz
!       Tpwp = SIGN(MAX(ABS(Tpwp),1.E-6),Tpwp)


!!  TKE budget  (Puhales, 2020, WRF 4.2.1)  << EOB

       !!!Shear Term
       !!!qSHEAR1D(k)=-(upwp*dudz + vpwp*dvdz)
       qSHEAR1D(k) = elq*sm(k)*gm(k) !staggered

       !!!Buoyancy Term
       !!!qBUOY1D(k)=grav*Tpwp/thl(k)
       !qBUOY1D(k)= elq*(sh(k)*gh(k) + gamv)
       !qBUOY1D(k) = elq*(sh(k)*(-dTdz*grav/thl(k)) + gamv) !! ORIGINAL CODE

       !! Buoyncy term takes the TKEprodTD(k) production now
       qBUOY1D(k) = elq*(sh(k)*gh(k)+gamv)+0.5*TKEprodTD(k) ! xmchen

       !!!Dissipation Term (now it evaluated in mym_predict)
       !qDISS1D(k) = (q3sq**(3./2.))/(b1*MAX(el(k),1.)) !! ORIGINAL CODE

       !! >> EOB
    ENDIF

    END DO
!

    dfm(kts) = 0.0
    dfh(kts) = 0.0
    dfq(kts) = 0.0
    tcd(kts) = 0.0
    qcd(kts) = 0.0

    tcd(kte) = 0.0
    qcd(kte) = 0.0

!
    DO k = kts,kte-1
       dzk = dz(k)
       tcd(k) = ( tcd(k+1)-tcd(k) )/( dzk )
       qcd(k) = ( qcd(k+1)-qcd(k) )/( dzk )
    END DO
!
    if (spp_pbl==1) then
       DO k = kts,kte
          dfm(k)= dfm(k) + dfm(k)* rstoch_col(k) * 1.5 * MAX(exp(-MAX(zw(k)-8000.,0.0)/2000.),0.001)
          dfh(k)= dfh(k) + dfh(k)* rstoch_col(k) * 1.5 * MAX(exp(-MAX(zw(k)-8000.,0.0)/2000.),0.001)
       END DO
    endif

!    RETURN
#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mym_turbulence


! ==================================================================
!>\ingroup gsd_mynn_edmf
  SUBROUTINE retrieve_exchange_coeffs(kts,kte,&
       &dfm,dfh,dz,K_m,K_h)

!-------------------------------------------------------------------

    integer , intent(in) :: kts,kte

    real(kind_phys), dimension(KtS:KtE), intent(in) :: dz,dfm,dfh

    real(kind_phys), dimension(KtS:KtE), intent(out) :: K_m, K_h


    integer :: k
    real(kind_phys):: dzk

    K_m(kts)=0.
    K_h(kts)=0.

    DO k=kts+1,kte
       dzk = 0.5  *( dz(k)+dz(k-1) )
       K_m(k)=dfm(k)*dzk
       K_h(k)=dfh(k)*dzk
    ENDDO

  END SUBROUTINE retrieve_exchange_coeffs

! ==================================================================
!>\ingroup gsd_mynn_edmf
  SUBROUTINE tridiag(n,a,b,c,d)

!! to solve system of linear eqs on tridiagonal matrix n times n
!! after Peaceman and Rachford, 1955
!! a,b,c,d - are vectors of order n
!! a,b,c - are coefficients on the LHS
!! d - is initially RHS on the output becomes a solution vector

!-------------------------------------------------------------------

    integer, intent(in):: n
    real(kind_phys), dimension(n), intent(in) :: a,b
    real(kind_phys), dimension(n), intent(inout) :: c,d

    integer :: i
    real(kind_phys):: p
    real(kind_phys), dimension(n) :: q

    c(n)=0.
    q(1)=-c(1)/b(1)
    d(1)=d(1)/b(1)

    DO i=2,n
       p=1./(b(i)+a(i)*q(i-1))
       q(i)=-c(i)*p
       d(i)=(d(i)-a(i)*d(i-1))*p
    ENDDO

    DO i=n-1,1,-1
       d(i)=d(i)+q(i)*d(i+1)
    ENDDO

  END SUBROUTINE tridiag

! ==================================================================
!>\ingroup gsd_mynn_edmf
      subroutine tridiag2(n,a,b,c,d,x)
      implicit none
!      a - sub-diagonal (means it is the diagonal below the main diagonal)
!      b - the main diagonal
!      c - sup-diagonal (means it is the diagonal above the main diagonal)
!      d - right part
!      x - the answer
!      n - number of unknowns (levels)

        integer,intent(in) :: n
        real(kind_phys), dimension(n), intent(in) :: a,b,c,d
        real(kind_phys), dimension(n), intent(out):: x
        real(kind_phys), dimension(n)  :: cp,dp
        real(kind_phys):: m
        integer :: i

        ! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
        ! solve for vectors c-prime and d-prime
        do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i))/m
        enddo
        ! initialize x
        x(n) = dp(n)
        ! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
           x(i) = dp(i)-cp(i)*x(i+1)
        end do

    end subroutine tridiag2
! ==================================================================
!>\ingroup gsd_mynn_edmf
       subroutine tridiag3(kte,a,b,c,d,x)

!ccccccccccccccccccccccccccccccc
! Aim: Inversion and resolution of a tridiagonal matrix
!          A X = D
! Input:
!  a(*) lower diagonal (Ai,i-1)
!  b(*) principal diagonal (Ai,i)
!  c(*) upper diagonal (Ai,i+1)
!  d
! Output
!  x     results
!ccccccccccccccccccccccccccccccc

       implicit none
        integer,intent(in)   :: kte
        integer, parameter   :: kts=1
        real(kind_phys), dimension(kte) :: a,b,c,d
        real(kind_phys), dimension(kte), intent(out) :: x
        integer :: in

!       integer kms,kme,kts,kte,in
!       real(kind_phys)a(kms:kme,3),c(kms:kme),x(kms:kme)

        do in=kte-1,kts,-1
         d(in)=d(in)-c(in)*d(in+1)/b(in+1)
         b(in)=b(in)-c(in)*a(in+1)/b(in+1)
        enddo

        do in=kts+1,kte
         d(in)=d(in)-a(in)*d(in-1)/b(in-1)
        enddo

        do in=kts,kte
         x(in)=d(in)/b(in)
        enddo

        return
        end subroutine tridiag3

! ==================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine calculates hybrid diagnotic boundary-layer height (PBLH).
!!
!! NOTES ON THE PBLH FORMULATION: The 1.5-theta-increase method defines
!!PBL heights as the level at.
!!which the potential temperature first exceeds the minimum potential.
!!temperature within the boundary layer by 1.5 K. When applied to.
!!observed temperatures, this method has been shown to produce PBL-
!!height estimates that are unbiased relative to profiler-based.
!!estimates (Nielsen-Gammon et al. 2008 \cite Nielsen_Gammon_2008).
!! However, their study did not
!!include LLJs. Banta and Pichugina (2008) \cite Pichugina_2008  show that a TKE-based.
!!threshold is a good estimate of the PBL height in LLJs. Therefore,
!!a hybrid definition is implemented that uses both methods, weighting
!!the TKE-method more during stable conditions (PBLH < 400 m).
!!A variable tke threshold (TKEeps) is used since no hard-wired
!!value could be found to work best in all conditions.
!>\section gen_get_pblh  GSD get_pblh General Algorithm
!> @{
  SUBROUTINE GET_PBLH(KTS,KTE,zi,thetav1D,qke1D,zw1D,dz1D,landsea,kzi)

    !---------------------------------------------------------------
    !             NOTES ON THE PBLH FORMULATION
    !
    !The 1.5-theta-increase method defines PBL heights as the level at
    !which the potential temperature first exceeds the minimum potential
    !temperature within the boundary layer by 1.5 K. When applied to
    !observed temperatures, this method has been shown to produce PBL-
    !height estimates that are unbiased relative to profiler-based
    !estimates (Nielsen-Gammon et al. 2008). However, their study did not
    !include LLJs. Banta and Pichugina (2008) show that a TKE-based
    !threshold is a good estimate of the PBL height in LLJs. Therefore,
    !a hybrid definition is implemented that uses both methods, weighting
    !the TKE-method more during stable conditions (PBLH < 400 m).
    !A variable tke threshold (TKEeps) is used since no hard-wired
    !value could be found to work best in all conditions.
    !---------------------------------------------------------------

    integer,intent(in) :: KTS,KTE

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    real(kind_phys), intent(out) :: zi
    real(kind_phys), intent(in) :: landsea
    real(kind_phys), dimension(kts:kte), intent(in) :: thetav1D, qke1D, dz1D
    real(kind_phys), dimension(kts:kte+1), intent(in) :: zw1D
    !LOCAL VARS
    real(kind_phys)::  PBLH_TKE,qtke,qtkem1,wt,maxqke,TKEeps,minthv
    real(kind_phys):: delt_thv   !delta theta-v; dependent on land/sea point
    real(kind_phys), parameter :: sbl_lim  = 200. !upper limit of stable BL height (m).
    real(kind_phys), parameter :: sbl_damp = 400. !transition length for blending (m).
    integer :: I,J,K,kthv,ktke,kzi

    !Initialize KPBL (kzi)
    kzi = 2

    !> - FIND MIN THETAV IN THE LOWEST 200 M AGL
    k = kts+1
    kthv = 1
    minthv = 9.E9
    DO WHILE (zw1D(k) .LE. 200.)
    !DO k=kts+1,kte-1
       IF (minthv > thetav1D(k)) then
           minthv = thetav1D(k)
           kthv = k
       ENDIF
       k = k+1
       !IF (zw1D(k) .GT. sbl_lim) exit
    ENDDO

    !> - FIND THETAV-BASED PBLH (BEST FOR DAYTIME).
    zi=0.
    k = kthv+1
    IF((landsea-1.5).GE.0)THEN
        ! WATER
        delt_thv = 1.0
    ELSE
        ! LAND
        delt_thv = 1.25
    ENDIF

    zi=0.
    k = kthv+1
!    DO WHILE (zi .EQ. 0.)
    DO k=kts+1,kte-1
       IF (thetav1D(k) .GE. (minthv + delt_thv))THEN
          zi = zw1D(k) - dz1D(k-1)* &
             & MIN((thetav1D(k)-(minthv + delt_thv))/ &
             & MAX(thetav1D(k)-thetav1D(k-1),1E-6),1.0)
       ENDIF
       !k = k+1
       IF (k .EQ. kte-1) zi = zw1D(kts+1) !EXIT SAFEGUARD
       IF (zi .NE. 0.0) exit
    ENDDO
    !print*,"IN GET_PBLH:",thsfc,zi

    !> - FOR STABLE BOUNDARY LAYERS, USE TKE METHOD TO COMPLEMENT THE
    !! THETAV-BASED DEFINITION (WHEN THE THETA-V BASED PBLH IS BELOW ~0.5 KM).
    !!THE TANH WEIGHTING FUNCTION WILL MAKE THE TKE-BASED DEFINITION NEGLIGIBLE
    !!WHEN THE THETA-V-BASED DEFINITION IS ABOVE ~1 KM.
    ktke = 1
    maxqke = MAX(Qke1D(kts),0.)
    !Use 5% of tke max (Kosovic and Curry, 2000; JAS)
    !TKEeps = maxtke/20. = maxqke/40.
    TKEeps = maxqke/40.
    TKEeps = MAX(TKEeps,0.02) !0.025)
    PBLH_TKE=0.

    k = ktke+1
!    DO WHILE (PBLH_TKE .EQ. 0.)
    DO k=kts+1,kte-1
       !QKE CAN BE NEGATIVE (IF CKmod == 0)... MAKE TKE NON-NEGATIVE.
       qtke  =MAX(Qke1D(k)/2.,0.)      ! maximum TKE
       qtkem1=MAX(Qke1D(k-1)/2.,0.)
       IF (qtke .LE. TKEeps) THEN
           PBLH_TKE = zw1D(k) - dz1D(k-1)* &
             & MIN((TKEeps-qtke)/MAX(qtkem1-qtke, 1E-6), 1.0)
           !IN CASE OF NEAR ZERO TKE, SET PBLH = LOWEST LEVEL.
           PBLH_TKE = MAX(PBLH_TKE,zw1D(kts+1))
           !print *,"PBLH_TKE:",i,PBLH_TKE, Qke1D(k)/2., zw1D(kts+1)
       ENDIF
       !k = k+1
       IF (k .EQ. kte-1) PBLH_TKE = zw1D(kts+1) !EXIT SAFEGUARD
       IF (PBLH_TKE .NE. 0.) exit
    ENDDO

    !> - With TKE advection turned on, the TKE-based PBLH can be very large
    !! in grid points with convective precipitation (> 8 km!),
    !! so an artificial limit is imposed to not let PBLH_TKE exceed the
    !!theta_v-based PBL height +/- 350 m.
    !!This has no impact on 98-99% of the domain, but is the simplest patch
    !!that adequately addresses these extremely large PBLHs.
    PBLH_TKE = MIN(PBLH_TKE,zi+350.)
    PBLH_TKE = MAX(PBLH_TKE,MAX(zi-350.,10.))

    wt=.5*TANH((zi - sbl_lim)/sbl_damp) + .5
    IF (maxqke <= 0.05) THEN
       !Cold pool situation - default to theta_v-based def
    ELSE
       !BLEND THE TWO PBLH TYPES HERE:
       zi=PBLH_TKE*(1.-wt) + zi*wt
    ENDIF

    !Compute KPBL (kzi)
    DO k=kts+1,kte-1
       IF ( zw1D(k) >= zi) THEN
          kzi = k-1
          exit
       ENDIF
    ENDDO

#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE GET_PBLH
!> @}

! ==================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine is the Dynamic Multi-Plume (DMP) Mass-Flux Scheme.
!!
!! dmp_mf() calculates the nonlocal turbulent transport from the dynamic
!! multiplume mass-flux scheme as well as the shallow-cumulus component of
!! the subgrid clouds. Note that this mass-flux scheme is called when the
!! namelist paramter \p bl_mynn_edmf is set to 1 (recommended).
!!
!! Much thanks to Kay Suslj of NASA-JPL for contributing the original version
!! of this mass-flux scheme. Considerable changes have been made from it's
!! original form. Some additions include:
!!  -# scale-aware tapering as dx -> 0
!!  -# transport of TKE (extra namelist option)
!!  -# Chaboureau-Bechtold cloud fraction & coupling to radiation (when icloud_bl > 0)
!!  -# some extra limits for numerical stability
!!
!! This scheme remains under development, so consider it experimental code.
!!
  SUBROUTINE DMP_mf(                            &
                 & kts,kte,dt,zw,dz,p,rho,      &
                 & momentum_opt,                &
                 & tke_opt,                     &
                 & scalar_opt,                  &
                 & u,v,w,th,thl,thv,tk,         &
                 & qt,qv,qc,qke,                &
                 & qnc,qni,qnwfa,qnifa,qnbca,   &
                 & exner,vt,vq,sgm,             &
                 & ust,flt,fltv,flq,flqv,       &
                 & pblh,kpbl,dx,landsea,ts,     &
            ! outputs - updraft properties
                 & edmf_a,edmf_w,               &
                 & edmf_qt,edmf_thl,            &
                 & edmf_ent,edmf_qc,            &
            ! outputs - variables needed for solver
                 & s_aw,s_awthl,s_awqt,         &
                 & s_awqv,s_awqc,               &
                 & s_awu,s_awv,s_awqke,         &
                 & s_awqnc,s_awqni,             &
                 & s_awqnwfa,s_awqnifa,         &
                 & s_awqnbca,                   &
                 & sub_thl,sub_sqv,             &
                 & sub_u,sub_v,                 &
                 & det_thl,det_sqv,det_sqc,     &
                 & det_u,det_v,                 &
            ! chem/smoke
                 & nchem,chem1,s_awchem,        &
                 & mix_chem,                    &
            ! in/outputs - subgrid scale clouds
                 & qc_bl1d,cldfra_bl1d,         &
                 & qc_bl1D_old,cldfra_bl1D_old, &
            ! inputs - flags for moist arrays
                 & F_QC,F_QI,                   &
                 & F_QNC,F_QNI,                 &
                 & F_QNWFA,F_QNIFA,F_QNBCA,     &
                 & Psig_shcu,                   &
            ! output info
                 & maxwidth,ktop,maxmf,ztop,    &
            ! inputs for stochastic perturbations
                 & spp_pbl,rstoch_col   )

  ! inputs:
     integer, intent(in) :: KTS,KTE,KPBL,momentum_opt,tke_opt,scalar_opt
#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

! Stochastic
     integer,  intent(in)                 :: spp_pbl
     real(kind_phys), dimension(kts:kte)  :: rstoch_col

     real(kind_phys),dimension(kts:kte), intent(in) ::                 &
          &U,V,W,TH,THL,TK,QT,QV,QC,                                   &
          &exner,dz,THV,P,rho,qke,qnc,qni,qnwfa,qnifa,qnbca
     real(kind_phys),dimension(kts:kte+1), intent(in) :: zw    !height at full-sigma
     real(kind_phys), intent(in) :: flt,fltv,flq,flqv,Psig_shcu,       &
          &landsea,ts,dx,dt,ust,pblh
     logical, optional :: F_QC,F_QI,F_QNC,F_QNI,F_QNWFA,F_QNIFA,F_QNBCA

  ! outputs - updraft properties
     real(kind_phys),dimension(kts:kte), intent(out) :: edmf_a,edmf_w, &
                      & edmf_qt,edmf_thl,edmf_ent,edmf_qc
     !add one local edmf variable:
     real(kind_phys),dimension(kts:kte) :: edmf_th
  ! output
     integer, intent(out) :: ktop
     real(kind_phys), intent(out) :: maxmf,ztop,maxwidth
  ! outputs - variables needed for solver
     real(kind_phys),dimension(kts:kte+1) :: s_aw,                     & !sum ai*rho*wis_awphi
          &s_awthl,s_awqt,s_awqv,s_awqc,s_awqnc,s_awqni,               &
          &s_awqnwfa,s_awqnifa,s_awqnbca,s_awu,s_awv,                  &
          &s_awqke,s_aw2

     real(kind_phys),dimension(kts:kte), intent(inout) ::              &
          &qc_bl1d,cldfra_bl1d,qc_bl1d_old,cldfra_bl1d_old

    integer, parameter :: nup=8, debug_mf=0
    real(kind_phys)    :: nup2

  !------------- local variables -------------------
  ! updraft properties defined on interfaces (k=1 is the top of the
  ! first model layer
     real(kind_phys),dimension(kts:kte+1,1:NUP) ::                     &
          &UPW,UPTHL,UPQT,UPQC,UPQV,                                   &
          &UPA,UPU,UPV,UPTHV,UPQKE,UPQNC,                              &
          &UPQNI,UPQNWFA,UPQNIFA,UPQNBCA
  ! entrainment variables
     real(kind_phys),dimension(kts:kte,1:NUP) :: ENT,ENTf
     integer,dimension(kts:kte,1:NUP)         :: ENTi
  ! internal variables
     integer :: K,I,k50
     real(kind_phys):: fltv2,wstar,qstar,thstar,sigmaW,sigmaQT,        &
          &sigmaTH,z0,pwmin,pwmax,wmin,wmax,wlv,Psig_w,maxw,maxqc,wpbl
     real(kind_phys):: B,QTn,THLn,THVn,QCn,Un,Vn,QKEn,QNCn,QNIn,       &
          &  QNWFAn,QNIFAn,QNBCAn,                                     &
          &  Wn2,Wn,EntEXP,EntEXM,EntW,BCOEFF,THVkm1,THVk,Pk,rho_int

  ! w parameters
     real(kind_phys), parameter ::                                     &
          &Wa=2./3.,                                                   &
          &Wb=0.002,                                                   &
          &Wc=1.5

  ! Lateral entrainment parameters ( L0=100 and ENT0=0.1) were taken from
  ! Suselj et al (2013, jas). Note that Suselj et al (2014,waf) use L0=200 and ENT0=0.2.
     real(kind_phys),parameter ::                                      &
          & L0=100.,                                                   &
          & ENT0=0.1

  ! Parameters/variables for regulating plumes:
     real(kind_phys), parameter :: Atot = 0.10 ! Maximum total fractional area of all updrafts
     real(kind_phys), parameter :: lmax = 1000.! diameter of largest plume (absolute maximum, can be smaller)
     real(kind_phys), parameter :: lmin = 300. ! diameter of smallest plume (absolute minimum, can be larger)
     real(kind_phys), parameter :: dlmin = 0.  ! delta increase in the diameter of smallest plume (large fltv)
     real(kind_phys)            :: minwidth    ! actual width of smallest plume
     real(kind_phys)            :: dl          ! variable increment of plume size
     real(kind_phys), parameter :: dcut = 1.2  ! max diameter of plume to parameterize relative to dx (km)
     real(kind_phys)::  d     != -2.3 to -1.7  ;=-1.9 in Neggers paper; power law exponent for number density (N=Cl^d).
          ! Note that changing d to -2.0 makes each size plume equally contribute to the total coverage of all plumes.
          ! Note that changing d to -1.7 doubles the area coverage of the largest plumes relative to the smallest plumes.
     real(kind_phys):: cn,c,l,n,an2,hux,wspd_pbl,cloud_base,width_flx

  ! chem/smoke
     integer, intent(in) :: nchem
     real(kind_phys),dimension(:, :) :: chem1
     real(kind_phys),dimension(kts:kte+1, nchem) :: s_awchem
     real(kind_phys),dimension(nchem) :: chemn
     real(kind_phys),dimension(kts:kte+1,1:NUP, nchem) :: UPCHEM
     integer :: ic
     real(kind_phys),dimension(kts:kte+1, nchem) :: edmf_chem
     logical, intent(in) :: mix_chem

  !JOE: add declaration of ERF
   real(kind_phys):: ERF

   logical :: superadiabatic

  ! VARIABLES FOR CHABOUREAU-BECHTOLD CLOUD FRACTION
   real(kind_phys),dimension(kts:kte), intent(inout) :: vt, vq, sgm
   real(kind_phys):: sigq,xl,rsl,cpm,a,qmq,mf_cf,Aup,Q1,diffqt,qsat_tk,&
           Fng,qww,alpha,beta,bb,f,pt,t,q2p,b9,satvp,rhgrid,           &
           Ac_mf,Ac_strat,qc_mf
   real(kind_phys), parameter :: cf_thresh = 0.5 ! only overwrite stratus CF less than this value

  ! Variables for plume interpolation/saturation check
   real(kind_phys),dimension(kts:kte) :: exneri,dzi,rhoz
   real(kind_phys):: THp, QTp, QCp, QCs, esat, qsl
   real(kind_phys):: csigma,acfac,ac_wsp

   !plume overshoot
   integer :: overshoot
   real(kind_phys):: bvf, Frz, dzp

   !Flux limiter: not let mass-flux of heat between k=1&2 exceed (fluxportion)*(surface heat flux).
   !This limiter makes adjustments to the entire column.
   real(kind_phys):: adjustment, flx1
   real(kind_phys), parameter :: fluxportion=0.75 ! set liberally, so has minimal impact. Note that
                                       ! 0.5 starts to have a noticeable impact
                                       ! over land (decrease maxMF by 10-20%), but no impact over water.

   !Subsidence
   real(kind_phys),dimension(kts:kte) :: sub_thl,sub_sqv,sub_u,sub_v,  &  !tendencies due to subsidence
                      det_thl,det_sqv,det_sqc,det_u,det_v,             &  !tendencied due to detrainment
                 envm_a,envm_w,envm_thl,envm_sqv,envm_sqc,             &
                                       envm_u,envm_v  !environmental variables defined at middle of layer
   real(kind_phys),dimension(kts:kte+1) ::  envi_a,envi_w        !environmental variables defined at model interface
   real(kind_phys):: temp,sublim,qc_ent,qv_ent,qt_ent,thl_ent,detrate, &
           detrateUV,oow,exc_fac,aratio,detturb,qc_grid,qc_sgs,        &
           qc_plume,exc_heat,exc_moist,tk_int,tvs
   real(kind_phys), parameter :: Cdet   = 1./45.
   real(kind_phys), parameter :: dzpmax = 300. !limit dz used in detrainment - can be excessing in thick layers
   !parameter "Csub" determines the propotion of upward vertical velocity that contributes to
   !environmenatal subsidence. Some portion is expected to be compensated by downdrafts instead of
   !gentle environmental subsidence. 1.0 assumes all upward vertical velocity in the mass-flux scheme
   !is compensated by "gentle" environmental subsidence.
   real(kind_phys), parameter :: Csub=0.25

   !Factor for the pressure gradient effects on momentum transport
   real(kind_phys), parameter :: pgfac = 0.00  ! Zhang and Wu showed 0.4 is more appropriate for lower troposphere
   real(kind_phys):: Uk,Ukm1,Vk,Vkm1,dxsa

! check the inputs
!     print *,'dt',dt
!     print *,'dz',dz
!     print *,'u',u
!     print *,'v',v
!     print *,'thl',thl
!     print *,'qt',qt
!     print *,'ust',ust
!     print *,'flt',flt
!     print *,'flq',flq
!     print *,'pblh',pblh

! Initialize individual updraft properties
  UPW=0.
  UPTHL=0.
  UPTHV=0.
  UPQT=0.
  UPA=0.
  UPU=0.
  UPV=0.
  UPQC=0.
  UPQV=0.
  UPQKE=0.
  UPQNC=0.
  UPQNI=0.
  UPQNWFA=0.
  UPQNIFA=0.
  UPQNBCA=0.
  if ( mix_chem ) then
     UPCHEM(kts:kte+1,1:NUP,1:nchem)=0.0
  endif

  ENT=0.001
! Initialize mean updraft properties
  edmf_a  =0.
  edmf_w  =0.
  edmf_qt =0.
  edmf_thl=0.
  edmf_ent=0.
  edmf_qc =0.
  if ( mix_chem ) then
     edmf_chem(kts:kte+1,1:nchem) = 0.0
  endif

! Initialize the variables needed for implicit solver
  s_aw=0.
  s_awthl=0.
  s_awqt=0.
  s_awqv=0.
  s_awqc=0.
  s_awu=0.
  s_awv=0.
  s_awqke=0.
  s_awqnc=0.
  s_awqni=0.
  s_awqnwfa=0.
  s_awqnifa=0.
  s_awqnbca=0.
  if ( mix_chem ) then
     s_awchem(kts:kte+1,1:nchem) = 0.0
  endif

! Initialize explicit tendencies for subsidence & detrainment
  sub_thl = 0.
  sub_sqv = 0.
  sub_u   = 0.
  sub_v   = 0.
  det_thl = 0.
  det_sqv = 0.
  det_sqc = 0.
  det_u   = 0.
  det_v   = 0.
  nup2    = nup !start with nup, but set to zero if activation criteria fails

  ! Taper off MF scheme when significant resolved-scale motions
  ! are present This function needs to be asymetric...
  maxw    = 0.0
  cloud_base  = 9000.0
  do k=1,kte-1
     if (zw(k) > pblh + 500.) exit

     wpbl = w(k)
     if (w(k) < 0.)wpbl = 2.*w(k)
     maxw = max(maxw,abs(wpbl))

     !Find highest k-level below 50m AGL
     if (ZW(k)<=50.)k50=k

     !Search for cloud base
     qc_sgs = max(qc(k), qc_bl1d(k))
     if (qc_sgs> 1E-5 .and. (cldfra_bl1d(k) .ge. 0.5) .and. cloud_base == 9000.0) then
       cloud_base = 0.5*(ZW(k)+ZW(k+1))
     endif
  enddo

  !do nothing for small w (< 1 m/s), but linearly taper off for w > 1.0 m/s
  maxw = max(0.,maxw - 1.0)
  Psig_w = max(0.0, 1.0 - maxw)
  Psig_w = min(Psig_w, Psig_shcu)

  !Completely shut off MF scheme for strong resolved-scale vertical velocities.
  fltv2 = fltv
  if(Psig_w == 0.0 .and. fltv > 0.0) fltv2 = -1.*fltv

  ! If surface buoyancy is positive we do integration, otherwise no.
  ! Also, ensure that it is at least slightly superadiabatic up through 50 m
  superadiabatic = .false.
  if ((landsea-1.5).ge.0) then
     hux = -0.001   ! WATER  ! dT/dz must be < - 0.1 K per 100 m.
  else
     hux = -0.005  ! LAND    ! dT/dz must be < - 0.5 K per 100 m.
  endif
  tvs = ts*(1.0+p608*qv(kts))
  do k=1,max(1,k50-1) !use "-1" because k50 used interface heights (zw).
    if (k == 1) then
      if ((thv(k)-tvs)/(0.5*dz(k)) < hux) then
        superadiabatic = .true.
      else
        superadiabatic = .false.
        exit
      endif
    else
      if ((thv(k)-thv(k-1))/(0.5*(dz(k)+dz(k-1))) < hux) then
        superadiabatic = .true.
      else
        superadiabatic = .false.
        exit
      endif
    endif
  enddo

  ! Determine the numer of updrafts/plumes in the grid column:
  ! Some of these criteria may be a little redundant but useful for bullet-proofing.
  !   (1) largest plume = 1.2 * dx.
  !   (2) Apply a scale-break, assuming no plumes with diameter larger than 1.1*PBLH can exist.
  !   (3) max plume size beneath clouds deck approx = 0.5 * cloud_base.
  !   (4) add wspd-dependent limit, when plume model breaks down. (hurricanes)
  !   (5) limit to reduce max plume sizes in weakly forced conditions. This is only
  !       meant to "soften" the activation of the mass-flux scheme.
  ! Criteria (1)
    maxwidth = min(dx*dcut, lmax)
  !Criteria (2)
    maxwidth = min(maxwidth, 1.1_kind_phys*PBLH)
  ! Criteria (3)
    if ((landsea-1.5) .lt. 0) then  !land
       maxwidth = MIN(maxwidth, 0.5_kind_phys*cloud_base)
    else                            !water
       maxwidth = MIN(maxwidth, 0.9_kind_phys*cloud_base)
    endif
  ! Criteria (4)
    wspd_pbl=SQRT(MAX(u(kts)**2 + v(kts)**2, 0.01_kind_phys))
    !Note: area fraction (acfac) is modified below
  ! Criteria (5) - only a function of flt (not fltv)
    if ((landsea-1.5).LT.0) then  !land
      width_flx = MAX(MIN(1000.*(0.6*tanh((fltv - 0.040)/0.04) + .5),1000._kind_phys), 0._kind_phys)
    else                          !water
      width_flx = MAX(MIN(1000.*(0.6*tanh((fltv - 0.007)/0.02) + .5),1000._kind_phys), 0._kind_phys)
    endif
    maxwidth = MIN(maxwidth, width_flx)
    minwidth = lmin
    !allow min plume size to increase in large flux conditions (eddy diffusivity should be
    !large enough to handle the representation of small plumes).
    if (maxwidth .ge. (lmax - 1.0) .and. fltv .gt. 0.2)minwidth = lmin + dlmin*min((fltv-0.2)/0.3, 1._kind_phys)

    if (maxwidth .le. minwidth) then ! deactivate MF component
       nup2 = 0
       maxwidth = 0.0
    endif

  ! Initialize values for 2d output fields:
    ktop = 0
    ztop = 0.0
    maxmf= 0.0

!Begin plume processing if passes criteria
if ( fltv2 > 0.002 .AND. (maxwidth > minwidth) .AND. superadiabatic) then

    ! Find coef C for number size density N
    cn = 0.
    d  =-1.9  !set d to value suggested by Neggers 2015 (JAMES).
    dl = (maxwidth - minwidth)/real(nup-1,kind=kind_phys)
    do i=1,NUP
       ! diameter of plume
       l = minwidth + dl*real(i-1)
       cn = cn + l**d * (l*l)/(dx*dx) * dl  ! sum fractional area of each plume
    enddo
    C = Atot/cn   !Normalize C according to the defined total fraction (Atot)

    ! Make updraft area (UPA) a function of the buoyancy flux
    if ((landsea-1.5).LT.0) then  !land
       acfac = .5*tanh((fltv2 - 0.02)/0.05) + .5
    else                          !water
       acfac = .5*tanh((fltv2 - 0.01)/0.03) + .5
    endif
    !add a windspeed-dependent adjustment to acfac that tapers off
    !the mass-flux scheme linearly above sfc wind speeds of 10 m/s.
    !Note: this effect may be better represented by an increase in
    !entrainment rate for high wind consitions (more ambient turbulence).
    if (wspd_pbl .le. 10.) then
       ac_wsp = 1.0
    else
       ac_wsp = 1.0 - min((wspd_pbl - 10.0)/15., 1.0)
    endif
    acfac  = acfac * ac_wsp

    ! Find the portion of the total fraction (Atot) of each plume size:
    An2 = 0.
    do i=1,NUP
       ! diameter of plume
       l  = minwidth + dl*real(i-1)
       N  = C*l**d                          ! number density of plume n
       UPA(1,i) = N*l*l/(dx*dx) * dl        ! fractional area of plume n

       UPA(1,i) = UPA(1,i)*acfac
       An2 = An2 + UPA(1,i)                 ! total fractional area of all plumes
       !print*," plume size=",l,"; area=",UPA(1,I),"; total=",An2
    end do

    ! set initial conditions for updrafts
    z0=50.
    pwmin=0.1       ! was 0.5
    pwmax=0.4       ! was 3.0

    wstar=max(1.E-2,(gtr*fltv2*pblh)**(onethird))
    qstar=max(flq,1.0E-5)/wstar
    thstar=flt/wstar

    if ((landsea-1.5) .ge. 0) then
       csigma = 1.34   ! WATER
    else
       csigma = 1.34   ! LAND
    endif

    if (env_subs) then
       exc_fac = 0.0
    else
       if ((landsea-1.5).GE.0) then
         !water: increase factor to compensate for decreased pwmin/pwmax
         exc_fac = 0.58*4.0
       else
         !land: no need to increase factor - already sufficiently large superadiabatic layers
         exc_fac = 0.58
       endif
    endif
    !decrease excess for large wind speeds
    exc_fac = exc_fac * ac_wsp

    !Note: sigmaW is typically about 0.5*wstar
    sigmaW =csigma*wstar*(z0/pblh)**(onethird)*(1 - 0.8*z0/pblh)
    sigmaQT=csigma*qstar*(z0/pblh)**(onethird)
    sigmaTH=csigma*thstar*(z0/pblh)**(onethird)

    !Note: Given the pwmin & pwmax set above, these max/mins are
    !      rarely exceeded.
    wmin=MIN(sigmaW*pwmin,0.1)
    wmax=MIN(sigmaW*pwmax,0.5)

    !SPECIFY SURFACE UPDRAFT PROPERTIES AT MODEL INTERFACE BETWEEN K = 1 & 2
    do i=1,NUP
       wlv=wmin+(wmax-wmin)/NUP2*(i-1)

       !SURFACE UPDRAFT VERTICAL VELOCITY
       UPW(1,I)=wmin + real(i)/real(NUP)*(wmax-wmin)
       UPU(1,I)=(U(KTS)*DZ(KTS+1)+U(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
       UPV(1,I)=(V(KTS)*DZ(KTS+1)+V(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
       UPQC(1,I)=0.0
       !UPQC(1,I)=(QC(KTS)*DZ(KTS+1)+QC(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))

       exc_heat = exc_fac*UPW(1,I)*sigmaTH/sigmaW
       UPTHV(1,I)=(THV(KTS)*DZ(KTS+1)+THV(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1)) &
           &     + exc_heat
       UPTHL(1,I)=(THL(KTS)*DZ(KTS+1)+THL(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1)) &
           &     + exc_heat

       !calculate exc_moist by use of surface fluxes
       exc_moist=exc_fac*UPW(1,I)*sigmaQT/sigmaW
       UPQT(1,I)=(QT(KTS)*DZ(KTS+1)+QT(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))&
            &     +exc_moist

       UPQKE(1,I)=(QKE(KTS)*DZ(KTS+1)+QKE(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
       UPQNC(1,I)=(QNC(KTS)*DZ(KTS+1)+QNC(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
       UPQNI(1,I)=(QNI(KTS)*DZ(KTS+1)+QNI(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
       UPQNWFA(1,I)=(QNWFA(KTS)*DZ(KTS+1)+QNWFA(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
       UPQNIFA(1,I)=(QNIFA(KTS)*DZ(KTS+1)+QNIFA(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
       UPQNBCA(1,I)=(QNBCA(KTS)*DZ(KTS+1)+QNBCA(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
    enddo

    if ( mix_chem ) then
      do i=1,NUP
        do ic = 1,nchem
          UPCHEM(1,i,ic)=(chem1(KTS,ic)*DZ(KTS+1)+chem1(KTS+1,ic)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
        enddo
      enddo
    endif

    !Initialize environmental variables which can be modified by detrainment
    envm_thl(kts:kte)=THL(kts:kte)
    envm_sqv(kts:kte)=QV(kts:kte)
    envm_sqc(kts:kte)=QC(kts:kte)
    envm_u(kts:kte)=U(kts:kte)
    envm_v(kts:kte)=V(kts:kte)
    do k=kts,kte-1
       rhoz(k)  = (rho(k)*dz(k+1)+rho(k+1)*dz(k))/(dz(k+1)+dz(k))
    enddo
    rhoz(kte) = rho(kte)

    !dxsa is scale-adaptive factor governing the pressure-gradient term of the momentum transport
    dxsa = 1. - MIN(MAX((12000.0-dx)/(12000.0-3000.0), 0.), 1.)

    ! do integration  updraft
    do i=1,NUP
       QCn = 0.
       overshoot = 0
       l  = minwidth + dl*real(i-1)            ! diameter of plume
       do k=kts+1,kte-1
          !Entrainment from Tian and Kuang (2016)
          !ENT(k,i) = 0.35/(MIN(MAX(UPW(K-1,I),0.75),1.9)*l)
          wmin = 0.3 + l*0.0005 !* MAX(pblh-ZW(k+1), 0.0)/pblh
          ENT(k,i) = 0.33/(MIN(MAX(UPW(K-1,I),wmin),0.9)*l)

          !Entrainment from Negggers (2015, JAMES)
          !ENT(k,i) = 0.02*l**-0.35 - 0.0009
          !ENT(k,i) = 0.04*l**-0.50 - 0.0009   !more plume diversity
          !ENT(k,i) = 0.04*l**-0.495 - 0.0009  !"neg1+"

          !Minimum background entrainment
          ENT(k,i) = max(ENT(k,i),0.0003)
          !ENT(k,i) = max(ENT(k,i),0.05/ZW(k))  !not needed for Tian and Kuang

          !increase entrainment for plumes extending very high.
          IF(ZW(k) >= MIN(pblh+1500., 4000.))THEN
            ENT(k,i)=ENT(k,i) + (ZW(k)-MIN(pblh+1500.,4000.))*5.0E-6
          ENDIF

          !SPP
          ENT(k,i) = ENT(k,i) * (1.0 - rstoch_col(k))

          ENT(k,i) = min(ENT(k,i),0.9/(ZW(k+1)-ZW(k)))

          ! Define environment U & V at the model interface levels
          Uk  =(U(k)*DZ(k+1)+U(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
          Ukm1=(U(k-1)*DZ(k)+U(k)*DZ(k-1))/(DZ(k-1)+DZ(k))
          Vk  =(V(k)*DZ(k+1)+V(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
          Vkm1=(V(k-1)*DZ(k)+V(k)*DZ(k-1))/(DZ(k-1)+DZ(k))

          ! Linear entrainment:
          EntExp= ENT(K,I)*(ZW(k+1)-ZW(k))
          EntExm= EntExp*0.3333    !reduce entrainment for momentum
          QTn =UPQT(k-1,I) *(1.-EntExp) + QT(k)*EntExp
          THLn=UPTHL(k-1,I)*(1.-EntExp) + THL(k)*EntExp
          Un  =UPU(k-1,I)  *(1.-EntExm) + U(k)*EntExm + dxsa*pgfac*(Uk - Ukm1)
          Vn  =UPV(k-1,I)  *(1.-EntExm) + V(k)*EntExm + dxsa*pgfac*(Vk - Vkm1)
          QKEn=UPQKE(k-1,I)*(1.-EntExp) + QKE(k)*EntExp
          QNCn=UPQNC(k-1,I)*(1.-EntExp) + QNC(k)*EntExp
          QNIn=UPQNI(k-1,I)*(1.-EntExp) + QNI(k)*EntExp
          QNWFAn=UPQNWFA(k-1,I)*(1.-EntExp) + QNWFA(k)*EntExp
          QNIFAn=UPQNIFA(k-1,I)*(1.-EntExp) + QNIFA(k)*EntExp
          QNBCAn=UPQNBCA(k-1,I)*(1.-EntExp) + QNBCA(k)*EntExp

          !capture the updated qc, qt & thl modified by entranment alone,
          !since they will be modified later if condensation occurs.
          qc_ent  = QCn
          qt_ent  = QTn
          thl_ent = THLn

          ! Exponential Entrainment:
          !EntExp= exp(-ENT(K,I)*(ZW(k)-ZW(k-1)))
          !QTn =QT(K) *(1-EntExp)+UPQT(K-1,I)*EntExp
          !THLn=THL(K)*(1-EntExp)+UPTHL(K-1,I)*EntExp
          !Un  =U(K)  *(1-EntExp)+UPU(K-1,I)*EntExp
          !Vn  =V(K)  *(1-EntExp)+UPV(K-1,I)*EntExp
          !QKEn=QKE(k)*(1-EntExp)+UPQKE(K-1,I)*EntExp

          if ( mix_chem ) then
            do ic = 1,nchem
              ! Exponential Entrainment:
              !chemn(ic) = chem(k,ic)*(1-EntExp)+UPCHEM(K-1,I,ic)*EntExp
              ! Linear entrainment:
              chemn(ic)=UPCHEM(k-1,i,ic)*(1.-EntExp) + chem1(k,ic)*EntExp
            enddo
          endif

          ! Define pressure at model interface
          Pk    =(P(k)*DZ(k+1)+P(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
          ! Compute plume properties thvn and qcn
          call condensation_edmf(QTn,THLn,Pk,ZW(k+1),THVn,QCn)

          ! Define environment THV at the model interface levels
          THVk  =(THV(k)*DZ(k+1)+THV(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
          THVkm1=(THV(k-1)*DZ(k)+THV(k)*DZ(k-1))/(DZ(k-1)+DZ(k))

!          B=g*(0.5*(THVn+UPTHV(k-1,I))/THV(k-1) - 1.0)
          B=grav*(THVn/THVk - 1.0)
          IF(B>0.)THEN
            BCOEFF = 0.15        !w typically stays < 2.5, so doesnt hit the limits nearly as much
          ELSE
            BCOEFF = 0.2 !0.33
          ENDIF

          ! Original StEM with exponential entrainment
          !EntW=exp(-2.*(Wb+Wc*ENT(K,I))*(ZW(k)-ZW(k-1)))
          !Wn2=UPW(K-1,I)**2*EntW + (1.-EntW)*0.5*Wa*B/(Wb+Wc*ENT(K,I))
          ! Original StEM with linear entrainment
          !Wn2=UPW(K-1,I)**2*(1.-EntExp) + EntExp*0.5*Wa*B/(Wb+Wc*ENT(K,I))
          !Wn2=MAX(Wn2,0.0)
          !WA: TEMF form
!          IF (B>0.0 .AND. UPW(K-1,I) < 0.2 ) THEN
          IF (UPW(K-1,I) < 0.2 ) THEN
             Wn = UPW(K-1,I) + (-2. * ENT(K,I) * UPW(K-1,I) + BCOEFF*B / MAX(UPW(K-1,I),0.2)) * MIN(ZW(k)-ZW(k-1), 250.)
          ELSE
             Wn = UPW(K-1,I) + (-2. * ENT(K,I) * UPW(K-1,I) + BCOEFF*B / UPW(K-1,I)) * MIN(ZW(k)-ZW(k-1), 250.)
          ENDIF
          !Do not allow a parcel to accelerate more than 1.25 m/s over 200 m.
          !Add max increase of 2.0 m/s for coarse vertical resolution.
          IF(Wn > UPW(K-1,I) + MIN(1.25*(ZW(k)-ZW(k-1))/200., 2.0) ) THEN
             Wn = UPW(K-1,I) + MIN(1.25*(ZW(k)-ZW(k-1))/200., 2.0)
          ENDIF
          !Add symmetrical max decrease in w
          IF(Wn < UPW(K-1,I) - MIN(1.25*(ZW(k)-ZW(k-1))/200., 2.0) ) THEN
             Wn = UPW(K-1,I) - MIN(1.25*(ZW(k)-ZW(k-1))/200., 2.0)
          ENDIF
          Wn = MIN(MAX(Wn,0.0), 3.0)

          !Check to make sure that the plume made it up at least one level.
          !if it failed, then set nup2=0 and exit the mass-flux portion.
          IF (k==kts+1 .AND. Wn == 0.) THEN
             NUP2=0
             exit
          ENDIF

          IF (debug_mf == 1) THEN
            IF (Wn .GE. 3.0) THEN
              ! surface values
              print *," **** SUSPICIOUSLY LARGE W:"
              print *,' QCn:',QCn,' ENT=',ENT(k,i),' Nup2=',Nup2
              print *,'pblh:',pblh,' Wn:',Wn,' UPW(k-1)=',UPW(K-1,I)
              print *,'K=',k,' B=',B,' dz=',ZW(k)-ZW(k-1)
            ENDIF
          ENDIF

          !Allow strongly forced plumes to overshoot if KE is sufficient
          IF (Wn <= 0.0 .AND. overshoot == 0) THEN
             overshoot = 1
             IF ( THVk-THVkm1 .GT. 0.0 ) THEN
                bvf = SQRT( gtr*(THVk-THVkm1)/dz(k) )
                !vertical Froude number
                Frz = UPW(K-1,I)/(bvf*dz(k))
                !IF ( Frz >= 0.5 ) Wn =  MIN(Frz,1.0)*UPW(K-1,I)
                dzp = dz(k)*MAX(MIN(Frz,1.0),0.0) ! portion of highest layer the plume penetrates
             ENDIF
          ELSE
             dzp = dz(k)
          ENDIF

          !minimize the plume penetratration in stratocu-topped PBL
          !IF (fltv2 < 0.06) THEN
          !   IF(ZW(k+1) >= pblh-200. .AND. qc(k) > 1e-5 .AND. I > 4) Wn=0.
          !ENDIF

          !Modify environment variables (representative of the model layer - envm*)
          !following the updraft dynamical detrainment of Asai and Kasahara (1967, JAS).
          !Reminder: w is limited to be non-negative (above)
          aratio   = MIN(UPA(K-1,I)/(1.-UPA(K-1,I)), 0.5) !limit should never get hit
          detturb  = 0.00008
          oow      = -0.060/MAX(1.0,(0.5*(Wn+UPW(K-1,I))))   !coef for dynamical detrainment rate
          detrate  = MIN(MAX(oow*(Wn-UPW(K-1,I))/dz(k), detturb), .0002) ! dynamical detrainment rate (m^-1)
          detrateUV= MIN(MAX(oow*(Wn-UPW(K-1,I))/dz(k), detturb), .0001) ! dynamical detrainment rate (m^-1)
          envm_thl(k)=envm_thl(k) + (0.5*(thl_ent + UPTHL(K-1,I)) - thl(k))*detrate*aratio*MIN(dzp,dzpmax)
          qv_ent = 0.5*(MAX(qt_ent-qc_ent,0.) + MAX(UPQT(K-1,I)-UPQC(K-1,I),0.))
          envm_sqv(k)=envm_sqv(k) + (qv_ent-QV(K))*detrate*aratio*MIN(dzp,dzpmax)
          IF (UPQC(K-1,I) > 1E-8) THEN
             IF (QC(K) > 1E-6) THEN
                qc_grid = QC(K)
             ELSE
                qc_grid = cldfra_bl1d(k)*qc_bl1d(K)
             ENDIF
             envm_sqc(k)=envm_sqc(k) + MAX(UPA(K-1,I)*0.5*(QCn + UPQC(K-1,I)) - qc_grid, 0.0)*detrate*aratio*MIN(dzp,dzpmax)
          ENDIF
          envm_u(k)  =envm_u(k)   + (0.5*(Un + UPU(K-1,I)) - U(K))*detrateUV*aratio*MIN(dzp,dzpmax)
          envm_v(k)  =envm_v(k)   + (0.5*(Vn + UPV(K-1,I)) - V(K))*detrateUV*aratio*MIN(dzp,dzpmax)

          IF (Wn > 0.) THEN
             !Update plume variables at current k index
             UPW(K,I)=Wn  !sqrt(Wn2)
             UPTHV(K,I)=THVn
             UPTHL(K,I)=THLn
             UPQT(K,I)=QTn
             UPQC(K,I)=QCn
             UPU(K,I)=Un
             UPV(K,I)=Vn
             UPQKE(K,I)=QKEn
             UPQNC(K,I)=QNCn
             UPQNI(K,I)=QNIn
             UPQNWFA(K,I)=QNWFAn
             UPQNIFA(K,I)=QNIFAn
             UPQNBCA(K,I)=QNBCAn
             UPA(K,I)=UPA(K-1,I)
             IF ( mix_chem ) THEN
               do ic = 1,nchem
                 UPCHEM(k,I,ic) = chemn(ic)
               enddo
             ENDIF
             ktop = MAX(ktop,k)
          ELSE
             exit  !exit k-loop
          END IF
       ENDDO

       IF (debug_mf == 1) THEN
          IF (MAXVAL(UPW(:,I)) > 10.0 .OR. MINVAL(UPA(:,I)) < 0.0 .OR. &
              MAXVAL(UPA(:,I)) > Atot .OR. NUP2 > 10) THEN
             ! surface values
             print *,'flq:',flq,' fltv:',fltv2,' Nup2=',Nup2
             print *,'pblh:',pblh,' wstar:',wstar,' ktop=',ktop
             print *,'sigmaW=',sigmaW,' sigmaTH=',sigmaTH,' sigmaQT=',sigmaQT
             ! means
             print *,'u:',u
             print *,'v:',v
             print *,'thl:',thl
             print *,'UPA:',UPA(:,I)
             print *,'UPW:',UPW(:,I)
             print *,'UPTHL:',UPTHL(:,I)
             print *,'UPQT:',UPQT(:,I)
             print *,'ENT:',ENT(:,I)
          ENDIF
       ENDIF
    ENDDO
ELSE
    !At least one of the conditions was not met for activating the MF scheme.
    NUP2=0.
END IF !end criteria check for mass-flux scheme

ktop=MIN(ktop,KTE-1)
IF (ktop == 0) THEN
   ztop = 0.0
ELSE
   ztop=zw(ktop)
ENDIF

IF (nup2 > 0) THEN
   !Calculate the fluxes for each variable
   !All s_aw* variable are == 0 at k=1
   DO i=1,NUP
      DO k=KTS,KTE-1
        s_aw(k+1)   = s_aw(k+1)    + rhoz(k)*UPA(K,i)*UPW(K,i)*Psig_w
        s_awthl(k+1)= s_awthl(k+1) + rhoz(k)*UPA(K,i)*UPW(K,i)*UPTHL(K,i)*Psig_w
        s_awqt(k+1) = s_awqt(k+1)  + rhoz(k)*UPA(K,i)*UPW(K,i)*UPQT(K,i)*Psig_w
        !to conform to grid mean properties, move qc to qv in grid mean
        !saturated layers, so total water fluxes are preserved but
        !negative qc fluxes in unsaturated layers is reduced.
!        if (qc(k) > 1e-12 .or. qc(k+1) > 1e-12) then
          qc_plume = UPQC(K,i)
!        else
!          qc_plume = 0.0
!        endif
        s_awqc(k+1) = s_awqc(k+1)  + rhoz(k)*UPA(K,i)*UPW(K,i)*qc_plume*Psig_w
        s_awqv(k+1) = s_awqt(k+1)  - s_awqc(k+1)
      ENDDO
   ENDDO
   !momentum
   if (momentum_opt > 0) then
      do i=1,nup
         do k=kts,kte-1
            s_awu(k+1)  = s_awu(k+1)   + rhoz(k)*UPA(K,i)*UPW(K,i)*UPU(K,i)*Psig_w
            s_awv(k+1)  = s_awv(k+1)   + rhoz(k)*UPA(K,i)*UPW(K,i)*UPV(K,i)*Psig_w
         enddo
      enddo
   endif
   !tke
   if (tke_opt > 0) then
      do i=1,nup
         do k=kts,kte-1
            s_awqke(k+1)= s_awqke(k+1) + rhoz(k)*UPA(K,i)*UPW(K,i)*UPQKE(K,i)*Psig_w
         enddo
      enddo
   endif
   !chem
   if ( mix_chem ) then
      do k=kts,kte
         do i=1,nup
            do ic = 1,nchem
              s_awchem(k+1,ic) = s_awchem(k+1,ic) + rhoz(k)*UPA(K,i)*UPW(K,i)*UPCHEM(K,i,ic)*Psig_w
            enddo
         enddo
      enddo
   endif

   if (scalar_opt > 0) then
      do k=kts,kte
         do I=1,nup
            s_awqnc(k+1)  = s_awqnc(K+1)   + rhoz(k)*UPA(K,i)*UPW(K,i)*UPQNC(K,i)*Psig_w
            s_awqni(k+1)  = s_awqni(K+1)   + rhoz(k)*UPA(K,i)*UPW(K,i)*UPQNI(K,i)*Psig_w
            s_awqnwfa(k+1)= s_awqnwfa(K+1) + rhoz(k)*UPA(K,i)*UPW(K,i)*UPQNWFA(K,i)*Psig_w
            s_awqnifa(k+1)= s_awqnifa(K+1) + rhoz(k)*UPA(K,i)*UPW(K,i)*UPQNIFA(K,i)*Psig_w
            s_awqnbca(k+1)= s_awqnbca(K+1) + rhoz(k)*UPA(K,i)*UPW(K,i)*UPQNBCA(K,i)*Psig_w
         enddo
      enddo
   endif

   !Flux limiter: Check ratio of heat flux at top of first model layer
   !and at the surface. Make sure estimated flux out of the top of the
   !layer is < fluxportion*surface_heat_flux
   IF (s_aw(kts+1) /= 0.) THEN
       dzi(kts) = 0.5*(DZ(kts)+DZ(kts+1)) !dz centered at model interface
       flx1   = MAX(s_aw(kts+1)*(TH(kts)-TH(kts+1))/dzi(kts),1.0e-5)
   ELSE
       flx1 = 0.0
       !print*,"ERROR: s_aw(kts+1) == 0, NUP=",NUP," NUP2=",NUP2,&
       !       " superadiabatic=",superadiabatic," KTOP=",KTOP
   ENDIF
   adjustment=1.0
   !Print*,"Flux limiter in MYNN-EDMF, adjustment=",fluxportion*flt/dz(kts)/flx1
   !Print*,"flt/dz=",flt/dz(kts)," flx1=",flx1," s_aw(kts+1)=",s_aw(kts+1)
   IF (flx1 > fluxportion*flt/dz(kts) .AND. flx1>0.0) THEN
       adjustment= fluxportion*flt/dz(kts)/flx1
       s_aw      = s_aw*adjustment
       s_awthl   = s_awthl*adjustment
       s_awqt    = s_awqt*adjustment
       s_awqc    = s_awqc*adjustment
       s_awqv    = s_awqv*adjustment
       s_awqnc   = s_awqnc*adjustment
       s_awqni   = s_awqni*adjustment
       s_awqnwfa = s_awqnwfa*adjustment
       s_awqnifa = s_awqnifa*adjustment
       s_awqnbca = s_awqnbca*adjustment
       IF (momentum_opt > 0) THEN
          s_awu  = s_awu*adjustment
          s_awv  = s_awv*adjustment
       ENDIF
       IF (tke_opt > 0) THEN
          s_awqke= s_awqke*adjustment
       ENDIF
       IF ( mix_chem ) THEN
          s_awchem = s_awchem*adjustment
       ENDIF
       UPA = UPA*adjustment
   ENDIF
   !Print*,"adjustment=",adjustment," fluxportion=",fluxportion," flt=",flt

   !Calculate mean updraft properties for output:
   !all edmf_* variables at k=1 correspond to the interface at top of first model layer
   do k=kts,kte-1
      do I=1,nup
         edmf_a(K)  =edmf_a(K)  +UPA(K,i)
         edmf_w(K)  =edmf_w(K)  +rhoz(k)*UPA(K,i)*UPW(K,i)
         edmf_qt(K) =edmf_qt(K) +rhoz(k)*UPA(K,i)*UPQT(K,i)
         edmf_thl(K)=edmf_thl(K)+rhoz(k)*UPA(K,i)*UPTHL(K,i)
         edmf_ent(K)=edmf_ent(K)+rhoz(k)*UPA(K,i)*ENT(K,i)
         edmf_qc(K) =edmf_qc(K) +rhoz(k)*UPA(K,i)*UPQC(K,i)
      enddo
   enddo
   do k=kts,kte-1
      !Note that only edmf_a is multiplied by Psig_w. This takes care of the
      !scale-awareness of the subsidence below:
      if (edmf_a(k)>0.) then
         edmf_w(k)=edmf_w(k)/edmf_a(k)
         edmf_qt(k)=edmf_qt(k)/edmf_a(k)
         edmf_thl(k)=edmf_thl(k)/edmf_a(k)
         edmf_ent(k)=edmf_ent(k)/edmf_a(k)
         edmf_qc(k)=edmf_qc(k)/edmf_a(k)
         edmf_a(k)=edmf_a(k)*Psig_w
         !FIND MAXIMUM MASS-FLUX IN THE COLUMN:
         if(edmf_a(k)*edmf_w(k) > maxmf) maxmf = edmf_a(k)*edmf_w(k)
      endif
   enddo ! end k

   !smoke/chem
   if ( mix_chem ) then
      do k=kts,kte-1
        do I=1,nup
          do ic = 1,nchem
            edmf_chem(k,ic) = edmf_chem(k,ic) + rhoz(k)*UPA(K,I)*UPCHEM(k,i,ic)
          enddo
        enddo
      enddo
      do k=kts,kte-1
        if (edmf_a(k)>0.) then
          do ic = 1,nchem
            edmf_chem(k,ic) = edmf_chem(k,ic)/edmf_a(k)
          enddo
        endif
      enddo ! end k
   endif

   !Calculate the effects environmental subsidence.
   !All envi_*variables are valid at the interfaces, like the edmf_* variables
   IF (env_subs) THEN
       DO k=kts+1,kte-1
          !First, smooth the profiles of w & a, since sharp vertical gradients
          !in plume variables are not likely extended to env variables
          !Note1: w is treated as negative further below
          !Note2: both w & a will be transformed into env variables further below
          envi_w(k) = onethird*(edmf_w(k-1)+edmf_w(k)+edmf_w(k+1))
          envi_a(k) = onethird*(edmf_a(k-1)+edmf_a(k)+edmf_a(k+1))*adjustment
       ENDDO
       !define env variables at k=1 (top of first model layer)
       envi_w(kts) = edmf_w(kts)
       envi_a(kts) = edmf_a(kts)
       !define env variables at k=kte
       envi_w(kte) = 0.0
       envi_a(kte) = edmf_a(kte)
       !define env variables at k=kte+1
       envi_w(kte+1) = 0.0
       envi_a(kte+1) = edmf_a(kte)
       !Add limiter for very long time steps (i.e. dt > 300 s)
       !Note that this is not a robust check - only for violations in
       !   the first model level.
       IF (envi_w(kts) > 0.9*DZ(kts)/dt) THEN
          sublim = 0.9*DZ(kts)/dt/envi_w(kts)
       ELSE
          sublim = 1.0
       ENDIF
       !Transform w & a into env variables
       DO k=kts,kte
          temp=envi_a(k)
          envi_a(k)=1.0-temp
          envi_w(k)=csub*sublim*envi_w(k)*temp/(1.-temp)
       ENDDO
       !calculate tendencies from subsidence and detrainment valid at the middle of
       !each model layer. The lowest model layer uses an assumes w=0 at the surface.
       dzi(kts)    = 0.5*(dz(kts)+dz(kts+1))
       sub_thl(kts)= 0.5*envi_w(kts)*envi_a(kts)*                               &
                     (rho(kts+1)*thl(kts+1)-rho(kts)*thl(kts))/dzi(kts)/rhoz(k)
       sub_sqv(kts)= 0.5*envi_w(kts)*envi_a(kts)*                               &
                     (rho(kts+1)*qv(kts+1)-rho(kts)*qv(kts))/dzi(kts)/rhoz(k)
       DO k=kts+1,kte-1
          dzi(k)    = 0.5*(dz(k)+dz(k+1))
          sub_thl(k)= 0.5*(envi_w(k)+envi_w(k-1))*0.5*(envi_a(k)+envi_a(k-1)) * &
                      (rho(k+1)*thl(k+1)-rho(k)*thl(k))/dzi(k)/rhoz(k)
          sub_sqv(k)= 0.5*(envi_w(k)+envi_w(k-1))*0.5*(envi_a(k)+envi_a(k-1)) * &
                      (rho(k+1)*qv(k+1)-rho(k)*qv(k))/dzi(k)/rhoz(k)
       ENDDO

       DO k=KTS,KTE-1
          det_thl(k)=Cdet*(envm_thl(k)-thl(k))*envi_a(k)*Psig_w
          det_sqv(k)=Cdet*(envm_sqv(k)-qv(k))*envi_a(k)*Psig_w
          det_sqc(k)=Cdet*(envm_sqc(k)-qc(k))*envi_a(k)*Psig_w
       ENDDO

       IF (momentum_opt > 0) THEN
         sub_u(kts)=0.5*envi_w(kts)*envi_a(kts)*                               &
                    (rho(kts+1)*u(kts+1)-rho(kts)*u(kts))/dzi(kts)/rhoz(k)
         sub_v(kts)=0.5*envi_w(kts)*envi_a(kts)*                               &
                    (rho(kts+1)*v(kts+1)-rho(kts)*v(kts))/dzi(kts)/rhoz(k)
         DO k=kts+1,kte-1
            sub_u(k)=0.5*(envi_w(k)+envi_w(k-1))*0.5*(envi_a(k)+envi_a(k-1)) * &
                      (rho(k+1)*u(k+1)-rho(k)*u(k))/dzi(k)/rhoz(k)
            sub_v(k)=0.5*(envi_w(k)+envi_w(k-1))*0.5*(envi_a(k)+envi_a(k-1)) * &
                      (rho(k+1)*v(k+1)-rho(k)*v(k))/dzi(k)/rhoz(k)
         ENDDO

         DO k=KTS,KTE-1
           det_u(k) = Cdet*(envm_u(k)-u(k))*envi_a(k)*Psig_w
           det_v(k) = Cdet*(envm_v(k)-v(k))*envi_a(k)*Psig_w
         ENDDO
       ENDIF
   ENDIF !end subsidence/env detranment

   !First, compute exner, plume theta, and dz centered at interface
   !Here, k=1 is the top of the first model layer. These values do not
   !need to be defined at k=kte (unused level).
   DO K=KTS,KTE-1
       exneri(k) = (exner(k)*dz(k+1)+exner(k+1)*dz(k))/(dz(k+1)+dz(k))
       edmf_th(k)= edmf_thl(k) + xlvcp/exneri(k)*edmf_qc(K)
       dzi(k)    = 0.5*(dz(k)+dz(k+1))
   ENDDO

!JOE: ADD CLDFRA_bl1d, qc_bl1d. Note that they have already been defined in
!     mym_condensation. Here, a shallow-cu component is added, but no cumulus
!     clouds can be added at k=1 (start loop at k=2).
   do k=kts+1,kte-2
      if (k > KTOP) exit
         if(0.5*(edmf_qc(k)+edmf_qc(k-1))>0.0 .and. (cldfra_bl1d(k) < cf_thresh))THEN
            !interpolate plume quantities to mass levels
            Aup = (edmf_a(k)*dzi(k-1)+edmf_a(k-1)*dzi(k))/(dzi(k-1)+dzi(k))
            THp = (edmf_th(k)*dzi(k-1)+edmf_th(k-1)*dzi(k))/(dzi(k-1)+dzi(k))
            QTp = (edmf_qt(k)*dzi(k-1)+edmf_qt(k-1)*dzi(k))/(dzi(k-1)+dzi(k))
            !convert TH to T
!            t = THp*exner(k)
            !SATURATED VAPOR PRESSURE
            esat = esat_blend(tk(k))
            !SATURATED SPECIFIC HUMIDITY
            qsl=ep_2*esat/max(1.e-7,(p(k)-ep_3*esat))

            !condensed liquid in the plume on mass levels
            if (edmf_qc(k)>0.0 .and. edmf_qc(k-1)>0.0) then
              QCp = (edmf_qc(k)*dzi(k-1)+edmf_qc(k-1)*dzi(k))/(dzi(k-1)+dzi(k))
            else
              QCp = max(edmf_qc(k),edmf_qc(k-1))
            endif

            !COMPUTE CLDFRA & QC_BL FROM MASS-FLUX SCHEME and recompute vt & vq
            xl = xl_blend(tk(k))                ! obtain blended heat capacity
            qsat_tk = qsat_blend(tk(k),p(k))    ! get saturation water vapor mixing ratio
                                                !   at t and p
            rsl = xl*qsat_tk / (r_v*tk(k)**2)   ! slope of C-C curve at t (abs temp)
                                                ! CB02, Eqn. 4
            cpm = cp + qt(k)*cpv                ! CB02, sec. 2, para. 1
            a   = 1./(1. + xl*rsl/cpm)          ! CB02 variable "a"
            b9  = a*rsl                         ! CB02 variable "b"

            q2p  = xlvcp/exner(k)
            pt = thl(k) +q2p*QCp*Aup ! potential temp (env + plume)
            bb = b9*tk(k)/pt ! bb is "b9" in BCMT95.  Their "b9" differs from
                           ! "b9" in CB02 by a factor
                           ! of T/theta.  Strictly, b9 above is formulated in
                           ! terms of sat. mixing ratio, but bb in BCMT95 is
                           ! cast in terms of sat. specific humidity.  The
                           ! conversion is neglected here.
            qww   = 1.+0.61*qt(k)
            alpha = 0.61*pt
            beta  = pt*xl/(tk(k)*cp) - 1.61*pt
            !Buoyancy flux terms have been moved to the end of this section...

            !Now calculate convective component of the cloud fraction:
            if (a > 0.0) then
               f = MIN(1.0/a, 4.0)              ! f is vertical profile scaling function (CB2005)
            else
               f = 1.0
            endif

            !CB form:
            !sigq = 3.5E-3 * Aup * 0.5*(edmf_w(k)+edmf_w(k-1)) * f  ! convective component of sigma (CB2005)
            !sigq = SQRT(sigq**2 + sgm(k)**2)    ! combined conv + stratus components
            !Per S.DeRoode 2009?
            !sigq = 5. * Aup * (QTp - qt(k))
            sigq = 10. * Aup * (QTp - qt(k))
            !constrain sigq wrt saturation:
            sigq = max(sigq, qsat_tk*0.02 )
            sigq = min(sigq, qsat_tk*0.25 )

            qmq = a * (qt(k) - qsat_tk)           ! saturation deficit/excess;
            Q1  = qmq/sigq                        !   the numerator of Q1

            if ((landsea-1.5).GE.0) then      ! WATER
               !modified form from LES
               !mf_cf = min(max(0.5 + 0.36 * atan(1.20*(Q1+0.2)),0.01),0.6)
               !Original CB
               mf_cf = min(max(0.5 + 0.36 * atan(1.55*Q1),0.01),0.6)
               mf_cf = max(mf_cf, 1.2 * Aup)
               mf_cf = min(mf_cf, 5.0 * Aup)
            else                              ! LAND
               !LES form
               !mf_cf = min(max(0.5 + 0.36 * atan(1.20*(Q1+0.4)),0.01),0.6)
               !Original CB
               mf_cf = min(max(0.5 + 0.36 * atan(1.55*Q1),0.01),0.6)
               mf_cf = max(mf_cf, 1.8 * Aup)
               mf_cf = min(mf_cf, 5.0 * Aup)
            endif

            !IF ( debug_code ) THEN
            !   print*,"In MYNN, StEM edmf"
            !   print*,"  CB: env qt=",qt(k)," qsat=",qsat_tk
            !   print*,"  k=",k," satdef=",QTp - qsat_tk," sgm=",sgm(k)
            !   print*,"  CB: sigq=",sigq," qmq=",qmq," tk=",tk(k)
            !   print*,"  CB: mf_cf=",mf_cf," cldfra_bl=",cldfra_bl1d(k)," edmf_a=",edmf_a(k)
            !ENDIF

            ! Update cloud fractions and specific humidities in grid cells
            ! where the mass-flux scheme is active. The specific humidities
            ! are converted to grid means (not in-cloud quantities).
            if ((landsea-1.5).GE.0) then     ! water
               if (QCp * Aup > 5e-5) then
                  qc_bl1d(k) = 1.86 * (QCp * Aup) - 2.2e-5
               else
                  qc_bl1d(k) = 1.18 * (QCp * Aup)
               endif
               cldfra_bl1d(k) = mf_cf
               Ac_mf          = mf_cf
            else                             ! land
               if (QCp * Aup > 5e-5) then
                  qc_bl1d(k) = 1.86 * (QCp * Aup) - 2.2e-5
               else
                  qc_bl1d(k) = 1.18 * (QCp * Aup)
               endif
               cldfra_bl1d(k) = mf_cf
               Ac_mf          = mf_cf
            endif

            !Now recalculate the terms for the buoyancy flux for mass-flux clouds:
            !See mym_condensation for details on these formulations.
            !Use Bechtold and Siebesma (1998) piecewise estimation of Fng with
            !limits ,since they really should be recalculated after all the other changes...:
            !Only overwrite vt & vq in non-stratus condition
            !if ((landsea-1.5).GE.0) then      ! WATER
               Q1=max(Q1,-2.25)
            !else
            !   Q1=max(Q1,-2.0)
            !endif

            if (Q1 .ge. 1.0) then
               Fng = 1.0
            elseif (Q1 .ge. -1.7 .and. Q1 .lt. 1.0) then
               Fng = EXP(-0.4*(Q1-1.0))
            elseif (Q1 .ge. -2.5 .and. Q1 .lt. -1.7) then
               Fng = 3.0 + EXP(-3.8*(Q1+1.7))
            else
               Fng = min(23.9 + EXP(-1.6*(Q1+2.5)), 60.)
            endif

            !link the buoyancy flux function to active clouds only (c*Aup):
            vt(k) = qww   - (1.5*Aup)*beta*bb*Fng - 1.
            vq(k) = alpha + (1.5*Aup)*beta*a*Fng  - tv0
         endif !check for (qc in plume) .and. (cldfra_bl < threshold)
      enddo !k-loop

ENDIF  !end nup2 > 0

!modify output (negative: dry plume, positive: moist plume)
if (ktop > 0) then
   maxqc = maxval(edmf_qc(1:ktop))
   if ( maxqc < 1.E-8) maxmf = -1.0*maxmf
endif

!
! debugging
!
if (edmf_w(1) > 4.0) then
! surface values
    print *,'flq:',flq,' fltv:',fltv2
    print *,'pblh:',pblh,' wstar:',wstar
    print *,'sigmaW=',sigmaW,' sigmaTH=',sigmaTH,' sigmaQT=',sigmaQT
! means
!   print *,'u:',u
!   print *,'v:',v
!   print *,'thl:',thl
!   print *,'thv:',thv
!   print *,'qt:',qt
!   print *,'p:',p

! updrafts
! DO I=1,NUP2
!   print *,'up:A',i
!   print *,UPA(:,i)
!   print *,'up:W',i
!   print*,UPW(:,i)
!   print *,'up:thv',i
!   print *,UPTHV(:,i)
!   print *,'up:thl',i
!   print *,UPTHL(:,i)
!   print *,'up:qt',i
!   print *,UPQT(:,i)
!   print *,'up:tQC',i
!   print *,UPQC(:,i)
!   print *,'up:ent',i
!   print *,ENT(:,i)
! ENDDO

! mean updrafts
   print *,' edmf_a',edmf_a(1:14)
   print *,' edmf_w',edmf_w(1:14)
   print *,' edmf_qt:',edmf_qt(1:14)
   print *,' edmf_thl:',edmf_thl(1:14)

ENDIF !END Debugging


#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

END SUBROUTINE DMP_MF
!=================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine
subroutine condensation_edmf(QT,THL,P,zagl,THV,QC)
!
! zero or one condensation for edmf: calculates THV and QC
!
real(kind_phys),intent(in)   :: QT,THL,P,zagl
real(kind_phys),intent(out)  :: THV
real(kind_phys),intent(inout):: QC

integer :: niter,i
real(kind_phys):: diff,exn,t,th,qs,qcold

! constants used from module_model_constants.F
! p1000mb
! rcp ... Rd/cp
! xlv ... latent heat for water (2.5e6)
! cp
! rvord .. r_v/r_d  (1.6)

! number of iterations
  niter=50
! minimum difference (usually converges in < 8 iterations with diff = 2e-5)
  diff=1.e-6

  EXN=(P/p1000mb)**rcp
  !QC=0.  !better first guess QC is incoming from lower level, do not set to zero
  do i=1,NITER
     T=EXN*THL + xlvcp*QC
     QS=qsat_blend(T,P)
     QCOLD=QC
     QC=0.5*QC + 0.5*MAX((QT-QS),0.)
     if (abs(QC-QCOLD)<Diff) exit
  enddo

  T=EXN*THL + xlvcp*QC
  QS=qsat_blend(T,P)
  QC=max(QT-QS,0.)

  !Do not allow saturation below 100 m
  if(zagl < 100.)QC=0.

  !THV=(THL+xlv/cp*QC).*(1+(1-rvovrd)*(QT-QC)-QC);
  THV=(THL+xlvcp*QC)*(1.+QT*(rvovrd-1.)-rvovrd*QC)

!  IF (QC > 0.0) THEN
!    PRINT*,"EDMF SAT, p:",p," iterations:",i
!    PRINT*," T=",T," THL=",THL," THV=",THV
!    PRINT*," QS=",QS," QT=",QT," QC=",QC,"ratio=",qc/qs
!  ENDIF

  !THIS BASICALLY GIVE THE SAME RESULT AS THE PREVIOUS LINE
  !TH = THL + xlv/cp/EXN*QC
  !THV= TH*(1. + p608*QT)

  !print *,'t,p,qt,qs,qc'
  !print *,t,p,qt,qs,qc


end subroutine condensation_edmf

!===============================================================

subroutine condensation_edmf_r(QT,THL,P,zagl,THV,QC)
!
! zero or one condensation for edmf: calculates THL and QC
! similar to condensation_edmf but with different inputs
!
real(kind_phys),intent(in)   :: QT,THV,P,zagl
real(kind_phys),intent(out)  :: THL, QC

integer :: niter,i
real(kind_phys):: diff,exn,t,th,qs,qcold

! number of iterations
  niter=50
! minimum difference
  diff=2.e-5

  EXN=(P/p1000mb)**rcp
  ! assume first that th = thv
  T = THV*EXN
  !QS = qsat_blend(T,P)
  !QC = QS - QT

  QC=0.

  do i=1,NITER
     QCOLD = QC
     T = EXN*THV/(1.+QT*(rvovrd-1.)-rvovrd*QC)
     QS=qsat_blend(T,P)
     QC= MAX((QT-QS),0.)
     if (abs(QC-QCOLD)<Diff) exit
  enddo
  THL = (T - xlv/cp*QC)/EXN

end subroutine condensation_edmf_r

!===============================================================

!>\ingroup gsd_mynn_edmf
!! \author JAYMES- added 22 Apr 2015
!! This function calculates saturation vapor pressure.  Separate ice and liquid functions
!! are used (identical to those in module_mp_thompson.F, v3.6). Then, the
!! final returned value is a temperature-dependant "blend". Because the final
!! value is "phase-aware", this formulation may be preferred for use throughout
!! the module (replacing "svp").
  FUNCTION esat_blend(t)

      IMPLICIT NONE

      real(kind_phys), intent(in):: t
      real(kind_phys):: esat_blend,XC,ESL,ESI,chi
      !liquid
      real(kind_phys), parameter:: J0= .611583699E03
      real(kind_phys), parameter:: J1= .444606896E02
      real(kind_phys), parameter:: J2= .143177157E01
      real(kind_phys), parameter:: J3= .264224321E-1
      real(kind_phys), parameter:: J4= .299291081E-3
      real(kind_phys), parameter:: J5= .203154182E-5
      real(kind_phys), parameter:: J6= .702620698E-8
      real(kind_phys), parameter:: J7= .379534310E-11
      real(kind_phys), parameter:: J8=-.321582393E-13
      !ice
      real(kind_phys), parameter:: K0= .609868993E03
      real(kind_phys), parameter:: K1= .499320233E02
      real(kind_phys), parameter:: K2= .184672631E01
      real(kind_phys), parameter:: K3= .402737184E-1
      real(kind_phys), parameter:: K4= .565392987E-3
      real(kind_phys), parameter:: K5= .521693933E-5
      real(kind_phys), parameter:: K6= .307839583E-7
      real(kind_phys), parameter:: K7= .105785160E-9
      real(kind_phys), parameter:: K8= .161444444E-12

      XC=MAX(-80.,t - t0c) !note t0c = 273.15, tice is set in module mynn_common to 240

! For 240 < t < 268.16 K, the vapor pressures are "blended" as a function of temperature,
! using the approach similar to Chaboureau and Bechtold (2002), JAS, p. 2363.  The resulting
! values are returned from the function.
      IF (t .GE. (t0c-6.)) THEN
          esat_blend = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
      ELSE IF (t .LE. tice) THEN
          esat_blend = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
      ELSE
          ESL = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
          ESI = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          chi = ((t0c-6.) - t)/((t0c-6.) - tice)
          esat_blend = (1.-chi)*ESL  + chi*ESI
      END IF

  END FUNCTION esat_blend

! ====================================================================

!>\ingroup gsd_mynn_edmf
!! This function extends function "esat" and returns a "blended"
!! saturation mixing ratio. Tice currently set to 240 K, t0c = 273.15 K.
!!\author JAYMES
  FUNCTION qsat_blend(t, P)

      IMPLICIT NONE

      real(kind_phys), intent(in):: t, P
      real(kind_phys):: qsat_blend,XC,ESL,ESI,RSLF,RSIF,chi
      !liquid
      real(kind_phys), parameter:: J0= .611583699E03
      real(kind_phys), parameter:: J1= .444606896E02
      real(kind_phys), parameter:: J2= .143177157E01
      real(kind_phys), parameter:: J3= .264224321E-1
      real(kind_phys), parameter:: J4= .299291081E-3
      real(kind_phys), parameter:: J5= .203154182E-5
      real(kind_phys), parameter:: J6= .702620698E-8
      real(kind_phys), parameter:: J7= .379534310E-11
      real(kind_phys), parameter:: J8=-.321582393E-13
      !ice
      real(kind_phys), parameter:: K0= .609868993E03
      real(kind_phys), parameter:: K1= .499320233E02
      real(kind_phys), parameter:: K2= .184672631E01
      real(kind_phys), parameter:: K3= .402737184E-1
      real(kind_phys), parameter:: K4= .565392987E-3
      real(kind_phys), parameter:: K5= .521693933E-5
      real(kind_phys), parameter:: K6= .307839583E-7
      real(kind_phys), parameter:: K7= .105785160E-9
      real(kind_phys), parameter:: K8= .161444444E-12

      XC=MAX(-80.,t - t0c)

      IF (t .GE. (t0c-6.)) THEN
          ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
          ESL  = min(ESL, P*0.15) ! Even with P=1050mb and T=55C, the sat. vap. pres only contributes to ~15% of total pres.
          qsat_blend = 0.622*ESL/max(P-ESL, 1e-5)
      ELSE IF (t .LE. tice) THEN
          ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          ESI  = min(ESI, P*0.15)
          qsat_blend = 0.622*ESI/max(P-ESI, 1e-5)
      ELSE
          ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
          ESL  = min(ESL, P*0.15)
          ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          ESI  = min(ESI, P*0.15)
          RSLF = 0.622*ESL/max(P-ESL, 1e-5)
          RSIF = 0.622*ESI/max(P-ESI, 1e-5)
!          chi  = (268.16-t)/(268.16-240.)
          chi  = ((t0c-6.) - t)/((t0c-6.) - tice)
         qsat_blend = (1.-chi)*RSLF + chi*RSIF
      END IF

  END FUNCTION qsat_blend

! ===================================================================

!>\ingroup gsd_mynn_edmf
!! This function interpolates the latent heats of vaporization and sublimation into
!! a single, temperature-dependent, "blended" value, following
!! Chaboureau and Bechtold (2002) \cite Chaboureau_2002, Appendix.
!!\author JAYMES
  FUNCTION xl_blend(t)

      IMPLICIT NONE

      real(kind_phys), intent(in):: t
      real(kind_phys):: xl_blend,xlvt,xlst,chi
      !note: t0c = 273.15, tice is set in mynn_common

      IF (t .GE. t0c) THEN
          xl_blend = xlv + (cpv-cliq)*(t-t0c)  !vaporization/condensation
      ELSE IF (t .LE. tice) THEN
          xl_blend = xls + (cpv-cice)*(t-t0c)  !sublimation/deposition
      ELSE
          xlvt = xlv + (cpv-cliq)*(t-t0c)  !vaporization/condensation
          xlst = xls + (cpv-cice)*(t-t0c)  !sublimation/deposition
!          chi  = (273.16-t)/(273.16-240.)
          chi  = (t0c - t)/(t0c - tice)
          xl_blend = (1.-chi)*xlvt + chi*xlst     !blended
      END IF

  END FUNCTION xl_blend

! ===================================================================

END MODULE module_bl_mynn_1d
