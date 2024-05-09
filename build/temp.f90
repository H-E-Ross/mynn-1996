  SUBROUTINE#[?4mmf(                            &
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
                 & edmf_a,edmf_w,               &
                 & edmf_qt,edmf_thl,            &
                 & edmf_ent,edmf_qc,            &
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
                 & nchem,chem1,s_awchem,        &
                 & mix_chem,                    &
                 & qc_bl1d,cldfra_bl1d,         &
                 & qc_bl1D_old,cldfra_bl1D_old, &
                 & F_QC,F_QI,                   &
                 & F_QNC,F_QNI,                 &
                 & F_QNWFA,F_QNIFA,F_QNBCA,     &
                 & Psig_shcu,                   &
                 & maxwidth,ktop,maxmf,ztop,    &
                 & spp_pbl,rstoch_col   )

     integer, intent(in) :: KTS,KTE,KPBL,momentum_opt,tke_opt,scalar_opt
#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

     integer,  intent(in)                 :: spp_pbl
     real(kind_phys), dimension(kts:kte)  :: rstoch_col

     real(kind_phys),dimension(kts:kte), intent(in) ::                 &
          &U,V,W,TH,THL,TK,QT,QV,QC,                                   &
          &exner,dz,THV,P,rho,qke,qnc,qni,qnwfa,qnifa,qnbca
     real(kind_phys),dimension(kts:kte+1), intent(in) :: zw    
     real(kind_phys), intent(in) :: flt,fltv,flq,flqv,Psig_shcu,       &
          &landsea,ts,dx,dt,ust,pblh
     logical, optional :: F_QC,F_QI,F_QNC,F_QNI,F_QNWFA,F_QNIFA,F_QNBCA
     real(kind_phys),dimension(kts:kte), intent(out) :: edmf_a,edmf_w, &
                      & edmf_qt,edmf_thl,edmf_ent,edmf_qc
     real(kind_phys),dimension(kts:kte) :: edmf_th
     integer, intent(out) :: ktop
     real(kind_phys), intent(out) :: maxmf,ztop,maxwidth
     real(kind_phys),dimension(kts:kte+1) :: s_aw,                     & 
          &s_awthl,s_awqt,s_awqv,s_awqc,s_awqnc,s_awqni,               &
          &s_awqnwfa,s_awqnifa,s_awqnbca,s_awu,s_awv,                  &
          &s_awqke,s_aw2

     real(kind_phys),dimension(kts:kte), intent(inout) ::              &
          &qc_bl1d,cldfra_bl1d,qc_bl1d_old,cldfra_bl1d_old

    integer, parameter :: nup=8, debug_mf=0
    real(kind_phys)    :: nup2
     real(kind_phys),dimension(kts:kte+1,1:NUP) ::                     &
          &UPW,UPTHL,UPQT,UPQC,UPQV,                                   &
          &UPA,UPU,UPV,UPTHV,UPQKE,UPQNC,                              &
          &UPQNI,UPQNWFA,UPQNIFA,UPQNBCA
     real(kind_phys),dimension(kts:kte,1:NUP) :: ENT,ENTf
     integer,dimension(kts:kte,1:NUP)         :: ENTi
     integer :: K,I,k50
     real(kind_phys):: fltv2,wstar,qstar,thstar,sigmaW,sigmaQT,        &
          &sigmaTH,z0,pwmin,pwmax,wmin,wmax,wlv,Psig_w,maxw,maxqc,wpbl
     real(kind_phys):: B,QTn,THLn,THVn,QCn,Un,Vn,QKEn,QNCn,QNIn,       &
          &  QNWFAn,QNIFAn,QNBCAn,                                     &
          &  Wn2,Wn,EntEXP,EntEXM,EntW,BCOEFF,THVkm1,THVk,Pk,rho_int
     real(kind_phys), parameter ::                                     &
          &Wa=2./3.,                                                   &
          &Wb=0.002,                                                   &
          &Wc=1.5

     real(kind_phys),parameter ::                                      &
          & L0=100.,                                                   &
          & ENT0=0.1

     real(kind_phys), parameter :: Atot = 0.10 
     real(kind_phys), parameter :: lmax = 1000.
     real(kind_phys), parameter :: lmin = 300. 
     real(kind_phys), parameter :: dlmin = 0.  
     real(kind_phys)            :: minwidth    
     real(kind_phys)            :: dl          
     real(kind_phys), parameter :: dcut = 1.2  
     real(kind_phys)::  d     
     real(kind_phys):: cn,c,l,n,an2,hux,wspd_pbl,cloud_base,width_flx

     integer, intent(in) :: nchem
     real(kind_phys),dimension(:, :) :: chem1
     real(kind_phys),dimension(kts:kte+1, nchem) :: s_awchem
     real(kind_phys),dimension(nchem) :: chemn
     real(kind_phys),dimension(kts:kte+1,1:NUP, nchem) :: UPCHEM
     integer :: ic
     real(kind_phys),dimension(kts:kte+1, nchem) :: edmf_chem
     logical, intent(in) :: mix_chem

   real(kind_phys):: ERF

   logical :: superadiabatic

   real(kind_phys),dimension(kts:kte), intent(inout) :: vt, vq, sgm
   real(kind_phys):: sigq,xl,rsl,cpm,a,qmq,mf_cf,Aup,Q1,diffqt,qsat_tk,&
           Fng,qww,alpha,beta,bb,f,pt,t,q2p,b9,satvp,rhgrid,           &
           Ac_mf,Ac_strat,qc_mf
   real(kind_phys), parameter :: cf_thresh = 0.5 

   real(kind_phys),dimension(kts:kte) :: exneri,dzi,rhoz
   real(kind_phys):: THp, QTp, QCp, QCs, esat, qsl
   real(kind_phys):: csigma,acfac,ac_wsp

   integer :: overshoot
   real(kind_phys):: bvf, Frz, dzp

   real(kind_phys):: adjustment, flx1
   real(kind_phys), parameter :: fluxportion=0.75 
   real(kind_phys),dimension(kts:kte) :: sub_thl,sub_sqv,sub_u,sub_v,  &  
                      det_thl,det_sqv,det_sqc,det_u,det_v,             &  
                 envm_a,envm_w,envm_thl,envm_sqv,envm_sqc,             &
                                       envm_u,envm_v  
   real(kind_phys),dimension(kts:kte+1) ::  envi_a,envi_w        
   real(kind_phys):: temp,sublim,qc_ent,qv_ent,qt_ent,thl_ent,detrate, &
           detrateUV,oow,exc_fac,aratio,detturb,qc_grid,qc_sgs,        &
           qc_plume,exc_heat,exc_moist,tk_int,tvs
   real(kind_phys), parameter :: Cdet   = 1./45.
   real(kind_phys), parameter :: dzpmax = 300. 
   real(kind_phys), parameter :: Csub=0.25

   real(kind_phys), parameter :: pgfac = 0.00  
   real(kind_phys):: Uk,Ukm1,Vk,Vkm1,dxsa

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
  edmf_a  =0.
  edmf_w  =0.
  edmf_qt =0.
  edmf_thl=0.
  edmf_ent=0.
  edmf_qc =0.
  if ( mix_chem ) then
     edmf_chem(kts:kte+1,1:nchem) = 0.0
  endif

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

  sub_thl = 0.
  sub_sqv = 0.
  sub_u   = 0.
  sub_v   = 0.
  det_thl = 0.
  det_sqv = 0.
  det_sqc = 0.
  det_u   = 0.
  det_v   = 0.
  nup2    = nup 

  maxw    = 0.0
  cloud_base  = 9000.0
  do k=1,kte-1
     if (zw(k) > pblh + 500.) exit

     wpbl = w(k)
     if (w(k) < 0.)wpbl = 2.*w(k)
     maxw = max(maxw,abs(wpbl))

     if (ZW(k)<=50.)k50=k

     qc_sgs = max(qc(k), qc_bl1d(k))
     if (qc_sgs> 1E-5 .and. (cldfra_bl1d(k) .ge. 0.5) .and. cloud_base == 9000.0) then
       cloud_base = 0.5*(ZW(k)+ZW(k+1))
     endif
  enddo

  maxw = max(0.,maxw - 1.0)
  Psig_w = max(0.0, 1.0 - maxw)
  Psig_w = min(Psig_w, Psig_shcu)

  fltv2 = fltv
  if(Psig_w == 0.0 .and. fltv > 0.0) fltv2 = -1.*fltv

  superadiabatic = .false.
  if ((landsea-1.5).ge.0) then
     hux = -0.001   
  else
     hux = -0.005  
  endif
  tvs = ts*(1.0+p608*qv(kts))
  do k=1,max(1,k50-1) 
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
    maxwidth = min(dx*dcut, lmax)
    maxwidth = min(maxwidth, 1.1_kind_phys*PBLH)
    if ((landsea-1.5) .lt. 0) then  
       maxwidth = MIN(maxwidth, 0.5_kind_phys*cloud_base)
    else                            
       maxwidth = MIN(maxwidth, 0.9_kind_phys*cloud_base)
    endif
    wspd_pbl=SQRT(MAX(u(kts)**2 + v(kts)**2, 0.01_kind_phys))
    if ((landsea-1.5).LT.0) then  
      width_flx = MAX(MIN(1000.*(0.6*tanh((fltv - 0.040)/0.04) + .5),1000._kind_phys), 0._kind_phys)
    else                          
      width_flx = MAX(MIN(1000.*(0.6*tanh((fltv - 0.007)/0.02) + .5),1000._kind_phys), 0._kind_phys)
    endif
    maxwidth = MIN(maxwidth, width_flx)
    minwidth = lmin
    if (maxwidth .ge. (lmax - 1.0) .and. fltv .gt. 0.2)minwidth = lmin + dlmin*min((fltv-0.2)/0.3, 1._kind_phys)
    if (maxwidth .le. minwidth) then 
       nup2 = 0
       maxwidth = 0.0
    endif

    ktop = 0
    ztop = 0.0
    maxmf= 0.0

if ( fltv2 > 0.002 .AND. (maxwidth > minwidth) .AND. superadiabatic) then

    cn = 0.
    d  =-1.9  
    dl = (maxwidth - minwidth)/real(nup-1,kind=kind_phys)
    do i=1,NUP
       l = minwidth + dl*real(i-1)
       cn = cn + l**d * (l*l)/(dx*dx) * dl 
    enddo
    C = Atot/cn   

    if ((landsea-1.5).LT.0) then  
       acfac = .5*tanh((fltv2 - 0.02)/0.05) + .5
    else                          
       acfac = .5*tanh((fltv2 - 0.01)/0.03) + .5
    endif
    if (wspd_pbl .le. 10.) then
       ac_wsp = 1.0
    else
       ac_wsp = 1.0 - min((wspd_pbl - 10.0)/15., 1.0)
    endif
    acfac  = acfac * ac_wsp

    An2 = 0.
    do i=1,NUP
       l  = minwidth + dl*real(i-1)
       N  = C*l**d                          
       UPA(1,i) = N*l*l/(dx*dx) * dl        

       UPA(1,i) = UPA(1,i)*acfac
       An2 = An2 + UPA(1,i)                 
    end do

    z0=50.
    pwmin=0.1       
    pwmax=0.4       

    wstar=max(1.E-2,(gtr*fltv2*pblh)**(onethird))
    qstar=max(flq,1.0E-5)/wstar
    thstar=flt/wstar

    if ((landsea-1.5) .ge. 0) then
       csigma = 1.34   
    else
       csigma = 1.34   
    endif

    if (env_subs) then
       exc_fac = 0.0
    else
       if ((landsea-1.5).GE.0) then
         exc_fac = 0.58*4.0
       else
         exc_fac = 0.58
       endif
    endif
    exc_fac = exc_fac * ac_wsp

    sigmaW =csigma*wstar*(z0/pblh)**(onethird)*(1 - 0.8*z0/pblh)
    sigmaQT=csigma*qstar*(z0/pblh)**(onethird)
    sigmaTH=csigma*thstar*(z0/pblh)**(onethird)

    wmin=MIN(sigmaW*pwmin,0.1)
    wmax=MIN(sigmaW*pwmax,0.5)

    do i=1,NUP
       wlv=wmin+(wmax-wmin)/NUP2*(i-1)

       UPW(1,I)=wmin + real(i)/real(NUP)*(wmax-wmin)
       UPU(1,I)=(U(KTS)*DZ(KTS+1)+U(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
       UPV(1,I)=(V(KTS)*DZ(KTS+1)+V(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1))
       UPQC(1,I)=0.0

       exc_heat = exc_fac*UPW(1,I)*sigmaTH/sigmaW
       UPTHV(1,I)=(THV(KTS)*DZ(KTS+1)+THV(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1)) &
           &     + exc_heat
       UPTHL(1,I)=(THL(KTS)*DZ(KTS+1)+THL(KTS+1)*DZ(KTS))/(DZ(KTS)+DZ(KTS+1)) &
           &     + exc_heat

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

    envm_thl(kts:kte)=THL(kts:kte)
    envm_sqv(kts:kte)=QV(kts:kte)
    envm_sqc(kts:kte)=QC(kts:kte)
    envm_u(kts:kte)=U(kts:kte)
    envm_v(kts:kte)=V(kts:kte)
    do k=kts,kte-1
       rhoz(k)  = (rho(k)*dz(k+1)+rho(k+1)*dz(k))/(dz(k+1)+dz(k))
    enddo
    rhoz(kte) = rho(kte)

    dxsa = 1. - MIN(MAX((12000.0-dx)/(12000.0-3000.0), 0.), 1.)

    do i=1,NUP
       QCn = 0.
       overshoot = 0
       l  = minwidth + dl*real(i-1)            
       do k=kts+1,kte-1
          wmin = 0.3 + l*0.0005 
          ENT(k,i) = 0.33/(MIN(MAX(UPW(K-1,I),wmin),0.9)*l)

          ENT(k,i) = max(ENT(k,i),0.0003)

          IF(ZW(k) >= MIN(pblh+1500., 4000.))THEN
            ENT(k,i)=ENT(k,i) + (ZW(k)-MIN(pblh+1500.,4000.))*5.0E-6
          ENDIF

          ENT(k,i) = ENT(k,i) * (1.0 – rstoch_col(k))
          ENT(k,i) = min(ENT(k,i),0.9/(ZW(k+1)-ZW(k)))

          Uk  =(U(k)*DZ(k+1)+U(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
          Ukm1=(U(k-1)*DZ(k)+U(k)*DZ(k-1))/(DZ(k-1)+DZ(k))
          Vk  =(V(k)*DZ(k+1)+V(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
          Vkm1=(V(k-1)*DZ(k)+V(k)*DZ(k-1))/(DZ(k-1)+DZ(k))

          EntExp= ENT(K,I)*(ZW(k+1)-ZW(k))
          EntExm= EntExp*0.3333    
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

          qc_ent  = QCn
          qt_ent  = QTn
          thl_ent = THLn


          if ( mix_chem ) then
            do ic = 1,nchem
              chemn(ic)=UPCHEM(k-1,i,ic)*(1.-EntExp) + chem1(k,ic)*EntExp
            enddo
          endif

          Pk    =(P(k)*DZ(k+1)+P(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
          call condensation_edmf(QTn,THLn,Pk,ZW(k+1),THVn,QCn)

          THVk  =(THV(k)*DZ(k+1)+THV(k+1)*DZ(k))/(DZ(k+1)+DZ(k))
          THVkm1=(THV(k-1)*DZ(k)+THV(k)*DZ(k-1))/(DZ(k-1)+DZ(k))

          B=grav*(THVn/THVk - 1.0)
          IF(B>0.)THEN
            BCOEFF = 0.15        
          ELSE
            BCOEFF = 0.2 
          ENDIF

          IF (UPW(K-1,I) < 0.2 ) THEN
             Wn = UPW(K-1,I) + (-2. * ENT(K,I) * UPW(K-1,I) + BCOEFF*B / MAX(UPW(K-1,I),0.2)) * MIN(ZW(k)-ZW(k-1), 250.)
          ELSE
             Wn = UPW(K-1,I) + (-2. * ENT(K,I) * UPW(K-1,I) + BCOEFF*B / UPW(K-1,I)) * MIN(ZW(k)-ZW(k-1), 250.)
          ENDIF
          IF(Wn > UPW(K-1,I) + MIN(1.25*(ZW(k)-ZW(k-1))/200., 2.0) ) THEN
             Wn = UPW(K-1,I) + MIN(1.25*(ZW(k)-ZW(k-1))/200., 2.0)
          ENDIF
          IF(Wn < UPW(K-1,I) - MIN(1.25*(ZW(k)-ZW(k-1))/200., 2.0) ) THEN
             Wn = UPW(K-1,I) - MIN(1.25*(ZW(k)-ZW(k-1))/200., 2.0)
          ENDIF
          Wn = MIN(MAX(Wn,0.0), 3.0)

          IF (k==kts+1 .AND. Wn == 0.) THEN
             NUP2=0
             exit
          ENDIF

          IF (debug_mf == 1) THEN
            IF (Wn .GE. 3.0) THEN
              print *," **** SUSPICIOUSLY LARGE W:"
              print *,' QCn:',QCn,' ENT=',ENT(k,i),' Nup2=',Nup2
              print *,'pblh:',pblh,' Wn:',Wn,' UPW(k-1)=',UPW(K-1,I)
              print *,'K=',k,' B=',B,' dz=',ZW(k)-ZW(k-1)
            ENDIF
          ENDIF

          IF (Wn <= 0.0 .AND. overshoot == 0) THEN
             overshoot = 1
             IF ( THVk-THVkm1 .GT. 0.0 ) THEN
                bvf = SQRT( gtr*(THVk-THVkm1)/dz(k) )
                Frz = UPW(K-1,I)/(bvf*dz(k))
                dzp = dz(k)*MAX(MIN(Frz,1.0),0.0) 
             ENDIF
          ELSE
             dzp = dz(k)
          ENDIF

          aratio   = MIN(UPA(K-1,I)/(1.-UPA(K-1,I)), 0.5) 
          detturb  = 0.00008
          oow      = -0.060/MAX(1.0,(0.5*(Wn+UPW(K-1,I))))   
          detrate  = MIN(MAX(oow*(Wn-UPW(K-1,I))/dz(k), detturb), .0002) 
          detrateUV= MIN(MAX(oow*(Wn-UPW(K-1,I))/dz(k), detturb), .0001) 
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
             UPW(K,I)=Wn  
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
             exit  
          END IF
       ENDDO

       IF (debug_mf == 1) THEN
          IF (MAXVAL(UPW(:,I)) > 10.0 .OR. MINVAL(UPA(:,I)) < 0.0 .OR. &
              MAXVAL(UPA(:,I)) > Atot .OR. NUP2 > 10) THEN
             print *,'flq:',flq,' fltv:',fltv2,' Nup2=',Nup2
             print *,'pblh:',pblh,' wstar:',wstar,' ktop=',ktop
             print *,'sigmaW=',sigmaW,' sigmaTH=',sigmaTH,' sigmaQT=',sigmaQT
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
    NUP2=0.
END IF 

ktop=MIN(ktop,KTE-1)
IF (ktop == 0) THEN
   ztop = 0.0
ELSE
   ztop=zw(ktop)
ENDIF

IF (nup2 > 0) THEN
   DO i=1,NUP
      DO k=KTS,KTE-1
        s_aw(k+1)   = s_aw(k+1)    + rhoz(k)*UPA(K,i)*UPW(K,i)*Psig_w
        s_awthl(k+1)= s_awthl(k+1) + rhoz(k)*UPA(K,i)*UPW(K,i)*UPTHL(K,i)*Psig_w
        s_awqt(k+1) = s_awqt(k+1)  + rhoz(k)*UPA(K,i)*UPW(K,i)*UPQT(K,i)*Psig_w
          qc_plume = UPQC(K,i)
        s_awqc(k+1) = s_awqc(k+1)  + rhoz(k)*UPA(K,i)*UPW(K,i)*qc_plume*Psig_w
        s_awqv(k+1) = s_awqt(k+1)  - s_awqc(k+1)
      ENDDO
   ENDDO
   if (momentum_opt > 0) then
      do i=1,nup
         do k=kts,kte-1
            s_awu(k+1)  = s_awu(k+1)   + rhoz(k)*UPA(K,i)*UPW(K,i)*UPU(K,i)*Psig_w
            s_awv(k+1)  = s_awv(k+1)   + rhoz(k)*UPA(K,i)*UPW(K,i)*UPV(K,i)*Psig_w
         enddo
      enddo
         endif
   if (tke_opt > 0) then
      do i=1,nup
         do k=kts,kte-1
            s_awqke(k+1)= s_awqke(k+1) + rhoz(k)*UPA(K,i)*UPW(K,i)*UPQKE(K,i)*Psig_w
         enddo
      enddo
   endif
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

   IF (s_aw(kts+1) /= 0.) THEN
       dzi(kts) = 0.5*(DZ(kts)+DZ(kts+1)) 
       flx1   = MAX(s_aw(kts+1)*(TH(kts)-TH(kts+1))/dzi(kts),1.0e-5)
   ELSE
       flx1 = 0.0
   ENDIF
   adjustment=1.0
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
      if (edmf_a(k)>0.) then
         edmf_w(k)=edmf_w(k)/edmf_a(k)
         edmf_qt(k)=edmf_qt(k)/edmf_a(k)
         edmf_thl(k)=edmf_thl(k)/edmf_a(k)
         edmf_ent(k)=edmf_ent(k)/edmf_a(k)
         edmf_qc(k)=edmf_qc(k)/edmf_a(k)
         edmf_a(k)=edmf_a(k)*Psig_w
         if(edmf_a(k)*edmf_w(k) > maxmf) maxmf = edmf_a(k)*edmf_w(k)
      endif
   enddo

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
      enddo 
   endif
           IF (env_subs) THEN
       DO k=kts+1,kte-1
          envi_w(k) = onethird*(edmf_w(k-1)+edmf_w(k)+edmf_w(k+1))
          envi_a(k) = onethird*(edmf_a(k-1)+edmf_a(k)+edmf_a(k+1))*adjustment
       ENDDO
       envi_w(kts) = edmf_w(kts)
       envi_a(kts) = edmf_a(kts)
       envi_w(kte) = 0.0
       envi_a(kte) = edmf_a(kte)
       envi_w(kte+1) = 0.0
       envi_a(kte+1) = edmf_a(kte)
       IF (envi_w(kts) > 0.9*DZ(kts)/dt) THEN
          sublim = 0.9*DZ(kts)/dt/envi_w(kts)
       ELSE
          sublim = 1.0
       ENDIF
       DO k=kts,kte
          temp=envi_a(k)
          envi_a(k)=1.0-temp
          envi_w(k)=csub*sublim*envi_w(k)*temp/(1.-temp)
       ENDDO
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
   ENDIF 

   DO K=KTS,KTE-1
       exneri(k) = (exner(k)*dz(k+1)+exner(k+1)*dz(k))/(dz(k+1)+dz(k))
       edmf_th(k)= edmf_thl(k) + xlvcp/exneri(k)*edmf_qc(K)
       dzi(k)    = 0.5*(dz(k)+dz(k+1))
   ENDDO

   do k=kts+1,kte-2
      if (k > KTOP) exit
         if(0.5*(edmf_qc(k)+edmf_qc(k-1))>0.0 .and. (cldfra_bl1d(k) < cf_thresh))THEN
            Aup = (edmf_a(k)*dzi(k-1)+edmf_a(k-1)*dzi(k))/(dzi(k-1)+dzi(k))
            THp = (edmf_th(k)*dzi(k-1)+edmf_th(k-1)*dzi(k))/(dzi(k-1)+dzi(k))
            QTp = (edmf_qt(k)*dzi(k-1)+edmf_qt(k-1)*dzi(k))/(dzi(k-1)+dzi(k))
            esat = esat_blend(tk(k))
            qsl=ep_2*esat/max(1.e-7,(p(k)-ep_3*esat))
            if (edmf_qc(k)>0.0 .and. edmf_qc(k-1)>0.0) then
              QCp = (edmf_qc(k)*dzi(k-1)+edmf_qc(k-1)*dzi(k))/(dzi(k-1)+dzi(k))
            else
              QCp = max(edmf_qc(k),edmf_qc(k-1))
            endif
            xl = xl_blend(tk(k))                
            qsat_tk = qsat_blend(tk(k),p(k))    
            rsl = xl*qsat_tk / (r_v*tk(k)**2)   
            cpm = cp + qt(k)*cpv                
            a   = 1./(1. + xl*rsl/cpm)          
            b9  = a*rsl                         
            q2p  = xlvcp/exner(k)
            pt = thl(k) +q2p*QCp*Aup 
            bb = b9*tk(k)/pt 
            qww   = 1.+0.61*qt(k)
            alpha = 0.61*pt
            beta  = pt*xl/(tk(k)*cp) - 1.61*pt
            if (a > 0.0) then
               f = MIN(1.0/a, 4.0)              
            else
               f = 1.0
            endif
            sigq = 10. * Aup * (QTp - qt(k))
            sigq = max(sigq, qsat_tk*0.02 )
            sigq = min(sigq, qsat_tk*0.25 )

            qmq = a * (qt(k) - qsat_tk)           
            Q1  = qmq/sigq                        

            if ((landsea-1.5).GE.0) then      
               mf_cf = min(max(0.5 + 0.36 * atan(1.55*Q1),0.01),0.6)
               mf_cf = max(mf_cf, 1.2 * Aup)
               mf_cf = min(mf_cf, 5.0 * Aup)
            else                              
               mf_cf = min(max(0.5 + 0.36 * atan(1.55*Q1),0.01),0.6)
               mf_cf = max(mf_cf, 1.8 * Aup)
               mf_cf = min(mf_cf, 5.0 * Aup)
            endif
            if ((landsea-1.5).GE.0) then     
               if (QCp * Aup > 5e-5) then
                  qc_bl1d(k) = 1.86 * (QCp * Aup) - 2.2e-5
               else
                  qc_bl1d(k) = 1.18 * (QCp * Aup)
               endif
               cldfra_bl1d(k) = mf_cf
               Ac_mf          = mf_cf
            else                       
               if (QCp * Aup > 5e-5) then
                  qc_bl1d(k) = 1.86 * (QCp * Aup) - 2.2e-5
               else
                  qc_bl1d(k) = 1.18 * (QCp * Aup)
               endif
               cldfra_bl1d(k) = mf_cf
               Ac_mf          = mf_cf
            endif
               Q1=max(Q1,-2.25)
            if (Q1 .ge. 1.0) then
               Fng = 1.0
            elseif (Q1 .ge. -1.7 .and. Q1 .lt. 1.0) then
               Fng = EXP(-0.4*(Q1-1.0))
            elseif (Q1 .ge. -2.5 .and. Q1 .lt. -1.7) then
               Fng = 3.0 + EXP(-3.8*(Q1+1.7))
            else
               Fng = min(23.9 + EXP(-1.6*(Q1+2.5)), 60.)
            endif
            vt(k) = qww   - (1.5*Aup)*beta*bb*Fng - 1.
            vq(k) = alpha + (1.5*Aup)*beta*a*Fng  - tv0
         endif 
      enddo 
ENDIF
if (ktop > 0) then
   maxqc = maxval(edmf_qc(1:ktop))
   if ( maxqc < 1.E-8) maxmf = -1.0*maxmf
endif
if (edmf_w(1) > 4.0) then
    print *,'flq:',flq,' fltv:',fltv2
    print *,'pblh:',pblh,' wstar:',wstar
    print *,'sigmaW=',sigmaW,' sigmaTH=',sigmaTH,' sigmaQT=',sigmaQT
   print *,' edmf_a',edmf_a(1:14)
   print *,' edmf_w',edmf_w(1:14)
   print *,' edmf_qt:',edmf_qt(1:14)
   print *,' edmf_thl:',edmf_thl(1:14)
ENDIF
#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif
END SUBROUTINE DMP_MF
