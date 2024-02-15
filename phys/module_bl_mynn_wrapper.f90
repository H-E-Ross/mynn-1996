







      MODULE module_bl_mynn_wrapper

      use module_bl_mynn_common

      contains




      subroutine mynnedmf_wrapper_init (              &
        &  RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,RQCBLTEN,&
        &  RQIBLTEN,QKE,                              &
        &  restart,allowed_to_read,                   &
        &  P_QC,P_QI,PARAM_FIRST_SCALAR,              &
        &  IDS,IDE,JDS,JDE,KDS,KDE,                   &
        &  IMS,IME,JMS,JME,KMS,KME,                   &
        &  ITS,ITE,JTS,JTE,KTS,KTE                    )

        implicit none
        
        LOGICAL,INTENT(IN) :: ALLOWED_TO_READ,RESTART

        INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE,  &
             &                IMS,IME,JMS,JME,KMS,KME,  &
             &                ITS,ITE,JTS,JTE,KTS,KTE


        REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(INOUT) :: &
             &RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,                 &
             &RQCBLTEN,RQIBLTEN,QKE

        INTEGER,  intent(in) :: P_QC,P_QI,PARAM_FIRST_SCALAR

        INTEGER :: I,J,K,ITF,JTF,KTF

        JTF=MIN0(JTE,JDE-1)
        KTF=MIN0(KTE,KDE-1)
        ITF=MIN0(ITE,IDE-1)

        IF (.NOT.RESTART) THEN
         DO J=JTS,JTF
          DO K=KTS,KTF
           DO I=ITS,ITF
              RUBLTEN(i,k,j)=0.
              RVBLTEN(i,k,j)=0.
              RTHBLTEN(i,k,j)=0.
              RQVBLTEN(i,k,j)=0.
              if( p_qc >= param_first_scalar ) RQCBLTEN(i,k,j)=0.
              if( p_qi >= param_first_scalar ) RQIBLTEN(i,k,j)=0.
           ENDDO
          ENDDO
         ENDDO
        ENDIF

      end subroutine mynnedmf_wrapper_init

      subroutine mynnedmf_wrapper_finalize ()
      end subroutine mynnedmf_wrapper_finalize





SUBROUTINE mynnedmf_wrapper_run(                 &
     &  initflag,restart,cycling,                &
     &  delt,dz,dxc,znt,                         &
     &  u,v,w,th,                                &
     &  qv,qc,qi,qs,qnc,qni,qnwfa,qnifa,qnbca,   &

     &  p,exner,rho,t3d,                         &
     &  xland,ts,qsfc,ps,                        &
     &  ust,ch,hfx,qfx,rmol,wspd,                &
     &  uoce,voce,                               &
     &  qke,qke_adv,sh3d,sm3d,                   &


     &  Tsq,Qsq,Cov,                             &
     &  rublten,rvblten,rthblten,                &
     &  rqvblten,rqcblten,rqiblten,rqsblten,     &
     &  rqncblten,rqniblten,                     &
     &  rqnwfablten,rqnifablten,rqnbcablten,     &

     &  exch_h,exch_m,pblh,kpbl,el_pbl,          &
     &  dqke,qwt,qshear,qbuoy,qdiss,             &
     &  qc_bl,qi_bl,cldfra_bl,                   &
     &  edmf_a,edmf_w,edmf_qt,                   &
     &  edmf_thl,edmf_ent,edmf_qc,               &
     &  sub_thl3d,sub_sqv3d,                     &
     &  det_thl3d,det_sqv3d,                     &
     &  maxwidth,maxMF,ztop_plume,ktop_plume,    &
     &  rthraten,                                &
     &  tke_budget,         bl_mynn_tkeadvect,   &
     &  bl_mynn_cloudpdf,   bl_mynn_mixlength,   &
     &  icloud_bl,          bl_mynn_edmf,        &
     &  bl_mynn_edmf_mom,   bl_mynn_edmf_tke,    &
     &  bl_mynn_cloudmix,   bl_mynn_mixqt,       &
     &  bl_mynn_output,     bl_mynn_closure,     &
     &  bl_mynn_mixscalars,                      &
     &  spp_pbl,pattern_spp_pbl,                 &
     &  flag_qc,flag_qi,flag_qs,                 &
     &  flag_qnc,flag_qni,                       &
     &  flag_qnwfa,flag_qnifa,flag_qnbca,        &
     &  ids,ide,jds,jde,kds,kde,                 &
     &  ims,ime,jms,jme,kms,kme,                 &
     &  its,ite,jts,jte,kts,kte                  )


     use module_bl_mynn_driver, only: mynn_bl_driver


     implicit none


     
     
     
     logical, parameter ::                                  &
     &       mix_chem   =.false.,                           &
     &       enh_mix    =.false.,                           &
     &       rrfs_sd    =.false.,                           &
     &       smoke_dbg  =.false.
     integer, parameter :: nchem=2, ndvel=2, kdvel=1,       &
             num_vert_mix = 1


     logical, intent(in) ::                                 &
     &       bl_mynn_tkeadvect,                             &
     &       cycling
      integer, intent(in) ::                                &
     &       bl_mynn_cloudpdf,                              &
     &       bl_mynn_mixlength,                             &
     &       icloud_bl,                                     &
     &       bl_mynn_edmf,                                  &
     &       bl_mynn_edmf_mom,                              &
     &       bl_mynn_edmf_tke,                              &
     &       bl_mynn_cloudmix,                              &
     &       bl_mynn_mixqt,                                 &
     &       bl_mynn_output,                                &
     &       bl_mynn_mixscalars,                            &
     &       spp_pbl,                                       &
     &       tke_budget
      real(kind_phys), intent(in) ::                        &
     &       bl_mynn_closure

      logical, intent(in) ::                                &
     &       FLAG_QI, FLAG_QNI, FLAG_QC, FLAG_QNC,          &
     &       FLAG_QS, FLAG_QNWFA, FLAG_QNIFA, FLAG_QNBCA
      logical, parameter :: FLAG_OZONE = .false.


      REAL(kind_phys),    intent(in) :: delt, dxc
      LOGICAL, intent(in) :: restart
      INTEGER :: i, j, k, itf, jtf, ktf, n
      INTEGER, intent(in) :: initflag,                      &
     &           IDS,IDE,JDS,JDE,KDS,KDE,                   &
     &           IMS,IME,JMS,JME,KMS,KME,                   &
     &           ITS,ITE,JTS,JTE,KTS,KTE


      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(in) ::               &
     &       u,v,w,t3d,th,rho,exner,p,dz
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(inout) ::            &
     &       rublten,rvblten,rthblten,                                                 &
     &       rqvblten,rqcblten,rqiblten,rqsblten,                                      &
     &       rqncblten,rqniblten,                                                      &
     &       rqnwfablten,rqnifablten,rqnbcablten 
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(inout) ::            &
     &       qke, qke_adv, el_pbl, sh3d, sm3d
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(inout) ::            &
     &       Tsq, Qsq, Cov, exch_h, exch_m
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(in) :: rthraten


      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), optional, intent(in)  ::    &
     &       pattern_spp_pbl
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), optional, intent(inout) ::  &
     &       qc_bl, qi_bl, cldfra_bl
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), optional, intent(inout) ::  &
     &       edmf_a,edmf_w,edmf_qt,                                                    &
     &       edmf_thl,edmf_ent,edmf_qc,                                                &
     &       sub_thl3d,sub_sqv3d,det_thl3d,det_sqv3d
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), optional, intent(inout) ::  &
     &       dqke,qWT,qSHEAR,qBUOY,qDISS
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), optional, intent(inout) ::  &
     &       qv,qc,qi,qs,qnc,qni,qnwfa,qnifa,qnbca


      real(kind_phys), allocatable, dimension(:,:) ::                                  &
     &       qc_bl2d, qi_bl2d, cldfra_bl2d, pattern_spp_pbl2d
      real(kind_phys), allocatable, dimension(:,:) ::                                  &
     &       edmf_a2d,edmf_w2d,edmf_qt2d,                                              &
     &       edmf_thl2d,edmf_ent2d,edmf_qc2d,                                          &
     &       sub_thl2d,sub_sqv2d,det_thl2d,det_sqv2d
      real(kind_phys), allocatable, dimension(:,:) ::                                  &
     &       dqke2d,qWT2d,qSHEAR2d,qBUOY2d,qDISS2d
      real(kind_phys), allocatable, dimension(:,:) ::                                  &
     &       qc2d,qi2d,qs2d,qnc2d,qni2d,qnwfa2d,qnifa2d,qnbca2d


      real(kind_phys), dimension(ims:ime,kms:kme,nchem) :: chem
      real(kind_phys), dimension(ims:ime,ndvel)         :: vd
      real(kind_phys), dimension(ims:ime) :: frp_mean, emis_ant_no


      real(kind_phys), dimension(ims:ime,jms:jme), intent(in) ::                       &
     &        xland,ts,qsfc,ps,ch
      real(kind_phys), dimension(ims:ime,jms:jme), intent(inout) ::                    &
     &        znt,pblh,maxwidth,maxmf,ztop_plume,rmol,hfx,qfx,ust,wspd,                &
     &        uoce,voce
      integer, dimension(ims:ime,jms:jme), intent(inout) ::                            &
     &        kpbl,ktop_plume


      real(kind_phys), dimension(ims:ime,kms:kme)         :: delp,sqv,sqc,sqi,sqs,ikzero
      real(kind_phys), dimension(ims:ime)                 :: dx
      logical, parameter                                  :: debug = .false.
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme) :: ozone,rO3blten

   
   
   
   

   jtf=MIN0(JTE,JDE-1)
   ktf=MIN0(KTE,KDE-1)
   itf=MIN0(ITE,IDE-1)

   
   ozone=0.0
   rO3blten=0.0
   ikzero=0.0

   
   if (icloud_bl > 0) then
      allocate(qc_bl2d(ims:ime,kms:kme))
      allocate(qi_bl2d(ims:ime,kms:kme))
      allocate(cldfra_bl2d(ims:ime,kms:kme))
      qc_bl2d=0.0
      qi_bl2d=0.0
      cldfra_bl2d=0.0
   endif
   if (spp_pbl > 0) then
      allocate(pattern_spp_pbl2d(ims:ime,kms:kme))
   endif
   if (bl_mynn_output > 0) then
      allocate(edmf_a2d(ims:ime,kms:kme))
      allocate(edmf_w2d(ims:ime,kms:kme))
      allocate(edmf_qt2d(ims:ime,kms:kme))
      allocate(edmf_thl2d(ims:ime,kms:kme))
      allocate(edmf_ent2d(ims:ime,kms:kme))
      allocate(edmf_qc2d(ims:ime,kms:kme))
      allocate(sub_thl2d(ims:ime,kms:kme))
      allocate(sub_sqv2d(ims:ime,kms:kme))
      allocate(det_thl2d(ims:ime,kms:kme))
      allocate(det_sqv2d(ims:ime,kms:kme))
   endif
   if (tke_budget .eq. 1) then
      allocate(dqke2d(ims:ime,kms:kme))
      allocate(qWT2d(ims:ime,kms:kme))
      allocate(qSHEAR2d(ims:ime,kms:kme))
      allocate(qBUOY2d(ims:ime,kms:kme))
      allocate(qDISS2d(ims:ime,kms:kme))
      dqke2d  =0.0
      qWT2d   =0.0
      qSHEAR2d=0.0
      qBUOY2d =0.0
      qDISS2d =0.0
   endif
   if (flag_qc) then
      allocate(qc2d(ims:ime,kms:kme))
      qc2d=0.0
   endif
   if (flag_qi) then
      allocate(qi2d(ims:ime,kms:kme))
      qi2d=0.0
   endif
   if (flag_qs) then
      allocate(qs2d(ims:ime,kms:kme))
      qs2d=0.0
   endif
   if (flag_qnc) then
      allocate(qnc2d(ims:ime,kms:kme))
      qnc2d=0.0
   endif
   if (flag_qni) then
      allocate(qni2d(ims:ime,kms:kme))
      qni2d=0.0
   endif
   if (flag_qnwfa) then
      allocate(qnwfa2d(ims:ime,kms:kme))
      qnwfa2d=0.0
   endif
   if (flag_qnifa) then
      allocate(qnifa2d(ims:ime,kms:kme))
      qnifa2d=0.0
   endif
   if (flag_qnbca) then
      allocate(qnbca2d(ims:ime,kms:kme))
      qnbca2d=0.0
   endif
   
   
   
   do j = jts,jtf

      
      if (icloud_bl > 0) then
         do k=kts,ktf
         do i=its,itf
            qc_bl2d(i,k)     = qc_bl(i,k,j)
            qi_bl2d(i,k)     = qi_bl(i,k,j)
            cldfra_bl2d(i,k) = cldfra_bl(i,k,j)
         enddo
         enddo
      endif

      
      if (spp_pbl > 0) then
         do k=kts,ktf
         do i=its,itf
            pattern_spp_pbl2d(i,k) = pattern_spp_pbl(i,k,j)
         enddo
         enddo
      endif

      
      if (flag_qc) then
         do k=kts,ktf
         do i=its,itf
            qc2d(i,k) = qc(i,k,j)
         enddo
         enddo
      endif
      if (flag_qi) then
         do k=kts,ktf
         do i=its,itf
            qi2d(i,k) = qi(i,k,j)
         enddo
         enddo
      endif
      if (flag_qs) then
         do k=kts,ktf
         do i=its,itf
            qs2d(i,k) = qs(i,k,j)
         enddo
         enddo
      endif
      if (flag_qnc) then
         do k=kts,ktf
         do i=its,itf
            qnc2d(i,k) = qnc(i,k,j)
         enddo
         enddo
      endif
      if (flag_qni) then
         do k=kts,ktf
         do i=its,itf
            qni2d(i,k) = qni(i,k,j)
         enddo
         enddo
      endif
      if (flag_qnwfa) then
         do k=kts,ktf
         do i=its,itf
            qnwfa2d(i,k) = qnwfa(i,k,j)
         enddo
         enddo
      endif
      if (flag_qnifa) then
         do k=kts,ktf
         do i=its,itf
            qnifa2d(i,k) = qnifa(i,k,j)
         enddo
         enddo
      endif
      if (flag_qnbca) then
         do k=kts,ktf
         do i=its,itf
            qnbca2d(i,k) = qnbca(i,k,j)
         enddo
         enddo
      endif

      chem     = 0.0
      vd       = 0.0
      frp_mean = 0.0
      emis_ant_no = 0.0

      
      
      do i=its,itf
         dx(i)=dxc
      enddo








      
      do k=kts,ktf
         do i=its,itf
             sqv(i,k)=qv(i,k,j)/(1.0 + qv(i,k,j))
             sqc(i,k)=qc2d(i,k)/(1.0 + qv(i,k,j))
         enddo
      enddo
      if (flag_qi) then
         do k=kts,ktf
           do i=its,itf
             sqi(i,k)=qi2d(i,k)/(1.0 + qv(i,k,j))
           enddo
         enddo
      else
         sqi(:,:)=0.0
      endif
      if (flag_qs) then
         do k=kts,ktf
           do i=its,itf
             sqs(i,k)=qs2d(i,k)/(1.0 + qv(i,k,j))
           enddo
         enddo
      else
         sqs(:,:)=0.0
      endif

      if (debug) then
         print*
         write(0,*)"===CALLING mynn_bl_driver; input:"
         print*,"tke_budget=",tke_budget
         print*,"bl_mynn_tkeadvect=",bl_mynn_tkeadvect
         print*,"bl_mynn_cloudpdf=",bl_mynn_cloudpdf
         print*,"bl_mynn_mixlength=",bl_mynn_mixlength
         print*,"bl_mynn_edmf=",bl_mynn_edmf
         print*,"bl_mynn_edmf_mom=",bl_mynn_edmf_mom
         print*,"bl_mynn_edmf_tke=",bl_mynn_edmf_tke
         print*,"bl_mynn_cloudmix=",bl_mynn_cloudmix
         print*,"bl_mynn_mixqt=",bl_mynn_mixqt
         print*,"icloud_bl=",icloud_bl
         print*,"T:",t3d(its,1,j),t3d(its,2,j),t3d(its,kte,j)
         print*,"TH:",th(its,1,j),th(its,2,j),th(its,kte,j)
         print*,"rho:",rho(its,1,j),rho(its,2,j),rho(its,kte,j)
         print*,"exner:",exner(its,1,j),exner(its,2,j),exner(its,kte,j)
         print*,"p:",p(its,1,j),p(its,2,j),p(its,kte,j)
         print*,"dz:",dz(its,1,j),dz(its,2,j),dz(its,kte,j)
         print*,"u:",u(its,1,j),u(its,2,j),u(its,kte,j)
         print*,"v:",v(its,1,j),v(its,2,j),v(its,kte,j)
         print*,"sqv:",sqv(its,1),sqv(its,2),sqv(its,kte)
         print*,"sqc:",sqc(its,1),sqc(its,2),sqc(its,kte)
         print*,"sqi:",sqi(its,1),sqi(its,2),sqi(its,kte)
         print*,"rmol:",rmol(its,j)," ust:",ust(its,j)
         print*,"dx=",dx(its),"initflag=",initflag
         print*,"Thetasurf:",ts(its,j)
         print*,"HFX:",hfx(its,j)," qfx",qfx(its,j)
         print*,"qsfc:",qsfc(its,j)," ps:",ps(its,j)
         print*,"wspd:",wspd(its,j)
         print*,"znt:",znt(its,j)," delt=",delt
         print*,"ite=",ite," kte=",kte
         print*,"PBLH=",pblh(its,j)," KPBL=",KPBL(its,j)," xland=",xland(its,j)
         print*," ch=",ch(its,j)
         print*,"qke:",qke(its,1,j),qke(its,2,j),qke(its,kte,j)
         print*,"el_pbl:",el_pbl(its,1,j),el_pbl(its,2,j),el_pbl(its,kte,j)
         print*,"Sh3d:",Sh3d(its,1,j),sh3d(its,2,j),sh3d(its,kte,j)
         print*,"max cf_bl:",maxval(cldfra_bl(its,:,j))
      endif


              CALL  mynn_bl_driver(                                    &
     &        initflag=initflag,restart=restart,cycling=cycling,       &
     &        delt=delt,dz=dz(:,:,j),dx=dx,znt=znt(:,j),               &
     &        u=u(:,:,j),v=v(:,:,j),w=w(:,:,j),                        &
     &        th=th(:,:,j),sqv3D=sqv,sqc3D=sqc,                        &
     &        sqi3D=sqi,sqs3D=sqs,qnc=qnc2d,qni=qni2d,                 &
     &        qnwfa=qnwfa2d,qnifa=qnifa2d,qnbca=qnbca2d,               &
     &        ozone=ozone(:,:,j),                                      &
     &        p=p(:,:,j),exner=exner(:,:,j),rho=rho(:,:,j),            &
     &        T3D=t3d(:,:,j),xland=xland(:,j),                         &
     &        ts=ts(:,j),qsfc=qsfc(:,j),ps=ps(:,j),                    &
     &        ust=ust(:,j),ch=ch(:,j),hfx=hfx(:,j),qfx=qfx(:,j),       &
     &        rmol=rmol(:,j),wspd=wspd(:,j),                           &
     &        uoce=uoce(:,j),voce=voce(:,j),                           & 
     &        qke=QKE(:,:,j),qke_adv=qke_adv(:,:,j),                   & 
     &        sh3d=Sh3d(:,:,j),sm3d=Sm3d(:,:,j),                       & 
     &        nchem=nchem,kdvel=kdvel,ndvel=ndvel,                     & 
     &        Chem3d=chem,Vdep=vd,                                     &
     &        FRP=frp_mean,EMIS_ANT_NO=emis_ant_no,                    &
     &        mix_chem=mix_chem,enh_mix=enh_mix,                       &
     &        rrfs_sd=rrfs_sd,smoke_dbg=smoke_dbg,                     & 
     &        tsq=tsq(:,:,j),qsq=qsq(:,:,j),cov=cov(:,:,j),            & 
     &        RUBLTEN=RUBLTEN(:,:,j),RVBLTEN=RVBLTEN(:,:,j),           & 
     &        RTHBLTEN=RTHBLTEN(:,:,j),RQVBLTEN=RQVBLTEN(:,:,j),       & 
     &        RQCBLTEN=rqcblten(:,:,j),RQIBLTEN=rqiblten(:,:,j),       & 
     &        RQNCBLTEN=rqncblten(:,:,j),RQNIBLTEN=rqniblten(:,:,j),   & 
     &        RQSBLTEN=ikzero,                                         & 
     &        RQNWFABLTEN=RQNWFABLTEN(:,:,j),                          & 
     &        RQNIFABLTEN=RQNIFABLTEN(:,:,j),                          & 
     &        RQNBCABLTEN=RQNBCABLTEN(:,:,j),                          & 
     &        dozone=rO3blten(:,:,j),                                  & 
     &        EXCH_H=exch_h(:,:,j),EXCH_M=exch_m(:,:,j),               & 
     &        pblh=pblh(:,j),KPBL=KPBL(:,j),                           & 
     &        el_pbl=el_pbl(:,:,j),                                    & 
     &        dqke=dqke2d,qWT=qWT2d,qSHEAR=qSHEAR2d,                   & 
     &        qBUOY=qBUOY2d,qDISS=qDISS2d,                             & 
     &        qc_bl=qc_bl2d,qi_bl=qi_bl2d,cldfra_bl=cldfra_bl2d,       & 
     &        bl_mynn_tkeadvect=bl_mynn_tkeadvect,                     & 
     &        tke_budget=tke_budget,                                   & 
     &        bl_mynn_cloudpdf=bl_mynn_cloudpdf,                       & 
     &        bl_mynn_mixlength=bl_mynn_mixlength,                     & 
     &        icloud_bl=icloud_bl,                                     & 
     &        closure=bl_mynn_closure,bl_mynn_edmf=bl_mynn_edmf,       & 
     &        bl_mynn_edmf_mom=bl_mynn_edmf_mom,                       & 
     &        bl_mynn_edmf_tke=bl_mynn_edmf_tke,                       & 
     &        bl_mynn_mixscalars=bl_mynn_mixscalars,                   & 
     &        bl_mynn_output=bl_mynn_output,                           & 
     &        bl_mynn_cloudmix=bl_mynn_cloudmix,                       & 
     &        bl_mynn_mixqt=bl_mynn_mixqt,                             & 
     &        edmf_a=edmf_a2d,edmf_w=edmf_w2d,                         & 
     &        edmf_qt=edmf_qt2d,edmf_thl=edmf_thl2d,                   & 
     &        edmf_ent=edmf_ent2d,edmf_qc=edmf_qc2d,                   & 
     &        sub_thl3D=sub_thl2d,sub_sqv3D=sub_sqv2d,                 & 
     &        det_thl3D=det_thl2d,det_sqv3D=det_sqv2d,                 & 
     &        maxwidth=maxwidth(:,j),maxMF=maxMF(:,j),                 & 
     &        ztop_plume=ztop_plume(:,j),ktop_plume=ktop_plume(:,j),   & 
     &        spp_pbl=spp_pbl,pattern_spp_pbl=pattern_spp_pbl2d,       & 
     &        RTHRATEN=rthraten(:,:,j),                                & 
     &        FLAG_QI=flag_qi,FLAG_QNI=flag_qni,FLAG_QS=flag_qs,       & 
     &        FLAG_QC=flag_qc,FLAG_QNC=flag_qnc,                       & 
     &        FLAG_QNWFA=FLAG_QNWFA,FLAG_QNIFA=FLAG_QNIFA,             & 
     &        FLAG_QNBCA=FLAG_QNBCA,FLAG_OZONE=flag_ozone,             & 
     &        IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde,         & 
     &        IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme,         & 
     &        ITS=its,ITE=itf,JTS=jts,JTE=jtf,KTS=kts,KTE=kte)           


     
      do k=kts,ktf
        do i=its,itf
           RQVBLTEN(i,k,j) = RQVBLTEN(i,k,j)/(1.0 - sqv(i,k))
           RQCBLTEN(i,k,j) = RQCBLTEN(i,k,j)/(1.0 - sqv(i,k))
           RQIBLTEN(i,k,j) = RQIBLTEN(i,k,j)/(1.0 - sqv(i,k))
        enddo
      enddo
      if (.false.) then 
        do k=kts,ktf
          do i=its,itf
            RQSBLTEN(i,k,j) = RQSBLTEN(i,k,j)/(1.0 - sqv(i,k))
          enddo
        enddo
      endif

     
      if (icloud_bl > 0) then
         do k=kts,ktf
         do i=its,itf
            qc_bl(i,k,j)     = qc_bl2d(i,k)/(1.0 - sqv(i,k))
            qi_bl(i,k,j)     = qi_bl2d(i,k)/(1.0 - sqv(i,k))
            cldfra_bl(i,k,j) = cldfra_bl2d(i,k)
         enddo
         enddo
      endif

      if (tke_budget .eq. 1) then
         do k=kts,ktf
         do i=its,itf
            dqke(i,k,j)     = dqke2d(i,k)
            qwt(i,k,j)      = qwt2d(i,k)
            qshear(i,k,j)   = qshear2d(i,k)
            qbuoy(i,k,j)    = qbuoy2d(i,k)
            qdiss(i,k,j)    = qdiss2d(i,k)
         enddo
         enddo
      endif

      if (bl_mynn_output > 0) then
         do k=kts,ktf
         do i=its,itf
            edmf_a(i,k,j)     = edmf_a2d(i,k)
            edmf_w(i,k,j)     = edmf_w2d(i,k)
            edmf_qt(i,k,j)    = edmf_qt2d(i,k)
            edmf_thl(i,k,j)   = edmf_thl2d(i,k)
            edmf_ent(i,k,j)   = edmf_ent2d(i,k)
            edmf_qc(i,k,j)    = edmf_qc2d(i,k)
            sub_thl3d(i,k,j)  = sub_thl2d(i,k)
            sub_sqv3d(i,k,j)  = sub_sqv2d(i,k)
            det_thl3d(i,k,j)  = det_thl2d(i,k)
            det_sqv3d(i,k,j)  = det_sqv2d(i,k)
         enddo
         enddo
      endif

      if (debug) then
          print*
          print*,"===Finished with mynn_bl_driver; output:"
          print*,"T:",t3d(its,1,j),t3d(its,2,j),t3d(its,kte,j)
          print*,"TH:",th(its,1,j),th(its,2,j),th(its,kte,j)
          print*,"rho:",rho(its,1,j),rho(its,2,j),rho(its,kte,j)
          print*,"exner:",exner(its,1,j),exner(its,2,j),exner(its,kte,j)
          print*,"p:",p(its,1,j),p(its,2,j),p(its,kte,j)
          print*,"dz:",dz(its,1,j),dz(its,2,j),dz(its,kte,j)
          print*,"u:",u(its,1,j),u(its,2,j),u(its,kte,j)
          print*,"v:",v(its,1,j),v(its,2,j),v(its,kte,j)
          print*,"sqv:",sqv(its,1),sqv(its,2),sqv(its,kte)
          print*,"sqc:",sqc(its,1),sqc(its,2),sqc(its,kte)
          print*,"sqi:",sqi(its,1),sqi(its,2),sqi(its,kte)
          print*,"rmol:",rmol(its,j)," ust:",ust(its,j)
          print*,"dx(its,j)=",dx(its),"initflag=",initflag
          print*,"Thetasurf:",ts(its,j)
          print*,"HFX:",hfx(its,j)," qfx",qfx(its,j)
          print*,"qsfc:",qsfc(its,j)," ps:",ps(its,j)
          print*,"wspd:",wspd(its,j)
          print*,"znt:",znt(its,j)," delt=",delt
          print*,"im=",ite," kte=",kte
          print*,"PBLH=",pblh(its,j)," KPBL=",KPBL(its,j)," xland=",xland(its,j)
          print*,"ch=",ch(its,j)
          print*,"qke:",qke(its,1,j),qke(its,2,j),qke(its,kte,j)
          print*,"el_pbl:",el_pbl(its,1,j),el_pbl(its,2,j),el_pbl(its,kte,j)
          print*,"Sh3d:",Sh3d(its,1,j),sh3d(its,2,j),sh3d(its,kte,j)
          print*,"exch_h:",exch_h(its,1,j),exch_h(its,2,j),exch_h(its,kte,j)
          print*,"exch_m:",exch_m(its,1,j),exch_m(its,2,j),exch_m(its,kte,j)
          print*,"max cf_bl:",maxval(cldfra_bl(its,:,j))
          print*,"max qc_bl:",maxval(qc_bl(its,:,j))
          print*,"dtdt:",rthblten(its,1,j),rthblten(its,2,j),rthblten(its,kte,j)
          print*,"dudt:",rublten(its,1,j),rublten(its,2,j),rublten(its,kte,j)
          print*,"dvdt:",rvblten(its,1,j),rvblten(its,2,j),rvblten(its,kte,j)
          print*,"dqdt:",rqvblten(its,1,j),rqvblten(its,2,j),rqvblten(its,kte,j)
          print*,"ztop_plume:",ztop_plume(its,j)," maxmf:",maxmf(its,j)
          print*
       endif

   enddo  

   
   if (bl_mynn_output > 0) then
      deallocate(edmf_a2d)
      deallocate(edmf_w2d)
      deallocate(edmf_qt2d)
      deallocate(edmf_thl2d)
      deallocate(edmf_ent2d)
      deallocate(edmf_qc2d)
      deallocate(sub_thl2d)
      deallocate(sub_sqv2d)
      deallocate(det_thl2d)
      deallocate(det_sqv2d)
   endif
   if (tke_budget .eq. 1) then
      deallocate(dqke2d)
      deallocate(qwt2d)
      deallocate(qshear2d)
      deallocate(qbuoy2d)
      deallocate(qdiss2d)
   endif
   if (icloud_bl > 0) then
      deallocate(qc_bl2d)
      deallocate(qi_bl2d)
      deallocate(cldfra_bl2d)
   endif
   if (flag_qc)   deallocate(qc2d)
   if (flag_qi)   deallocate(qi2d)
   if (flag_qs)   deallocate(qs2d)
   if (flag_qnc)  deallocate(qnc2d)
   if (flag_qni)  deallocate(qni2d)
   if (flag_qnwfa)deallocate(qnwfa2d)
   if (flag_qnifa)deallocate(qnifa2d)
   if (flag_qnbca)deallocate(qnbca2d)
   if (spp_pbl > 0) then
      deallocate(pattern_spp_pbl2d)
   endif



  CONTAINS


  SUBROUTINE moisture_check2(kte, delt, dp, exner, &
                             qv, qc, qi, th        )
  
  
  
  
  
  
  
  
  
  

    implicit none
    integer,  intent(in)     :: kte
    real, intent(in)     :: delt
    real, dimension(kte), intent(in)     :: dp
    real, dimension(kte), intent(in)     :: exner
    real, dimension(kte), intent(inout)  :: qv, qc, qi, th
    integer   k
    real ::  dqc2, dqi2, dqv2, sum, aa, dum
    real, parameter :: qvmin1= 1e-8,    & 
                       qvmin = 1e-20,   & 
                       qcmin = 0.0,     &
                       qimin = 0.0

    do k = kte, 1, -1  
       dqc2 = max(0.0, qcmin-qc(k)) 
       dqi2 = max(0.0, qimin-qi(k)) 

       
       qc(k)  = qc(k)  +  dqc2
       qi(k)  = qi(k)  +  dqi2
       qv(k)  = qv(k)  -  dqc2 - dqi2
       
       
       
       
       th(k)  = th(k)  +  xlvcp*dqc2 + &
                          xlscp*dqi2

       
       if (k .eq. 1) then
          dqv2   = max(0.0, qvmin1-qv(k)) 
          qv(k)  = qv(k)  + dqv2
          qv(k)  = max(qv(k),qvmin1)
          dqv2   = 0.0
       else
          dqv2   = max(0.0, qvmin-qv(k)) 
          qv(k)  = qv(k)  + dqv2
          qv(k-1)= qv(k-1)  - dqv2*dp(k)/dp(k-1)
          qv(k)  = max(qv(k),qvmin)
       endif
       qc(k) = max(qc(k),qcmin)
       qi(k) = max(qi(k),qimin)
    end do



    
    if( dqv2 .gt. 1.e-20 ) then
        sum = 0.0
        do k = 1, kte
           if( qv(k) .gt. 2.0*qvmin ) sum = sum + qv(k)*dp(k)
        enddo
        aa = dqv2*dp(1)/max(1.e-20,sum)
        if( aa .lt. 0.5 ) then
            do k = 1, kte
               if( qv(k) .gt. 2.0*qvmin ) then
                   dum    = aa*qv(k)
                   qv(k)  = qv(k) - dum
               endif
            enddo
        else
        

        endif
    endif

    return

  END SUBROUTINE moisture_check2

  END SUBROUTINE mynnedmf_wrapper_run



END MODULE module_bl_mynn_wrapper
