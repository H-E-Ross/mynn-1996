!>\file module_bl_mynn.F90
!! This file contains the entity of MYNN-EDMF PBL scheme.
! **********************************************************************
! *   An improved Mellor-Yamada turbulence closure model               *
! *                                                                    *
! *      Original author: M. Nakanishi (N.D.A), naka@nda.ac.jp         *
! *      Translated into F90 and implemented in WRF-ARW by:            *
! *                       Mariusz Pagowski (NOAA-GSL)                  *
! *      Subsequently developed by:                                    *
! *                 Joseph Olson, Jaymes Kenyon (NOAA/GSL),            *
! *                 Wayne Angevine (NOAA/CSL), Kay Suselj (NASA/JPL),  *
! *                 Franciano Puhales (UFSM), Laura Fowler (NCAR),     *
! *                 Elynn Wu (UCSD), and Jordan Schnell (NOAA/GSL)     *
! *                                                                    *
! *   Contents:                                                        *
! *                                                                    *
! *   mynn_bl_driver - main subroutine which calls all other routines  *
! *   --------------                                                   *
! *     1. mym_initialize  (to be called once initially)               *
! *        gives the closure constants and initializes the turbulent   *
! *        quantities.                                                 *
! *     2. get_pblh                                                    *
! *        Calculates the boundary layer height                        *
! *     3. scale_aware                                                 *
! *        Calculates scale-adaptive tapering functions                *
! *     4. mym_condensation                                            *
! *        determines the liquid water content and the cloud fraction  *
! *        diagnostically.                                             *
! *     5. dmp_mf                                                      *
! *        Calls the (nonlocal) mass-flux component                    *
! *     6. ddmf_jpl                                                    *
! *        Calls the downdraft mass-flux component                     *
! *    (-) mym_level2      (called in the other subroutines)           *
! *        calculates the stability functions at Level 2.              *
! *    (-) mym_length      (called in the other subroutines)           *
! *        calculates the master length scale.                         *
! *     7. mym_turbulence                                              *
! *        calculates the vertical diffusivity coefficients and the    *
! *        production terms for the turbulent quantities.              *
! *     8. mym_predict                                                 *
! *        predicts the turbulent quantities at the next step.         *
! *                                                                    *
! *             call mym_initialize                                    *
! *                  |                                                 *
! *                  |<----------------+                               *
! *                  |                 |                               *
! *             call get_pblh          |                               *
! *             call scale_aware       |                               *
! *             call mym_condensation  |                               *
! *             call dmp_mf            |                               *
! *             call ddmf_jpl          |                               *
! *             call mym_turbulence    |                               *
! *             call mym_predict       |                               *
! *                  |                 |                               *
! *                  |-----------------+                               *
! *                  |                                                 *
! *                 end                                                *
! *                                                                    *
! *   Variables worthy of special mention:                             *
! *     tref   : Reference temperature                                 *
! *     thl    : Liquid water potential temperature                    *
! *     qw     : Total water (water vapor+liquid water) content        *
! *     ql     : Liquid water content                                  *
! *     vt, vq : Functions for computing the buoyancy flux             *
! *     qke    : 2 * TKE                                               *
! *     el     : mixing length                                         *
! *                                                                    *
! *     If the water contents are unnecessary, e.g., in the case of    *
! *     ocean models, thl is the potential temperature and qw, ql, vt  *
! *     and vq are all zero.                                           *
! *                                                                    *
! *   Grid arrangement:                                                *
! *             k+1 +---------+                                        *
! *                 |         |     i = 1 - nx                         *
! *             (k) |    *    |     k = 1 - nz                         *
! *                 |         |                                        *
! *              k  +---------+                                        *
! *                 i   (i)  i+1                                       *
! *                                                                    *
! *     All the predicted variables are defined at the center (*) of   *
! *     the grid boxes. The diffusivity coefficients and two of their  *
! *     components (el and stability functions sh & sm) are, however,  *
! *     defined on the walls of the grid boxes.                        *
! *     # Upper boundary values are given at k=nz.                     *
! *                                                                    *
! *   References:                                                      *
! *     1. Nakanishi, M., 2001:                                        *
! *        Boundary-Layer Meteor., 99, 349-378.                        *
! *     2. Nakanishi, M. and H. Niino, 2004:                           *
! *        Boundary-Layer Meteor., 112, 1-31.                          *
! *     3. Nakanishi, M. and H. Niino, 2006:                           *
! *        Boundary-Layer Meteor., 119, 397-407.                       *
! *     4. Nakanishi, M. and H. Niino, 2009:                           *
! *        Jour. Meteor. Soc. Japan, 87, 895-912.                      *
! *     5. Olson J. and coauthors, 2019: A description of the          *
! *        MYNN-EDMF scheme and coupling to other components in        *
! *        WRF-ARW. NOAA Tech. Memo. OAR GSD, 61, 37 pp.,              *
! *        https://doi.org/10.25923/n9wm-be49.                         * 
! *     6. Puhales, Franciano S. and coauthors, 2020: Turbulent        *
! *        Kinetic Energy Budget for MYNN-EDMF PBL Scheme in WRF model.*
! *        Universidade Federal de Santa Maria Technical Note. 9 pp.   *
! **********************************************************************
! ==================================================================
! Notes on original implementation into WRF-ARW
! changes to original code:
! 1. code is 1D (in z)
! 2. option to advect TKE, but not the covariances and variances
! 3. Cranck-Nicholson replaced with the implicit scheme
! 4. removed terrain-dependent grid since input in WRF in actual
!    distances in z[m]
! 5. cosmetic changes to adhere to WRF standard (remove common blocks,
!            intent etc)
!-------------------------------------------------------------------
! Further modifications post-implementation
!
! 1. Addition of BouLac mixing length in the free atmosphere.
! 2. Changed the turbulent mixing length to be integrated from the
!    surface to the top of the BL + a transition layer depth.
! v3.4.1:    Option to use Kitamura/Canuto modification which removes 
!            the critical Richardson number and negative TKE (default).
!            Hybrid PBL height diagnostic, which blends a theta-v-based
!            definition in neutral/convective BL and a TKE-based definition
!            in stable conditions.
!            TKE budget output option
! v3.5.0:    TKE advection option (bl_mynn_tkeadvect)
! v3.5.1:    Fog deposition related changes.
! v3.6.0:    Removed fog deposition from the calculation of tendencies
!            Added mixing of qc, qi, qni
!            Added output for wstar, delta, TKE_PBL, & KPBL for correct 
!                   coupling to shcu schemes  
! v3.8.0:    Added subgrid scale cloud output for coupling to radiation
!            schemes (activated by setting icloud_bl =1 in phys namelist).
!            Added WRF_DEBUG prints (at level 3000)
!            Added Tripoli and Cotton (1981) correction.
!            Added namelist option bl_mynn_cloudmix to test effect of mixing
!                cloud species (default = 1: on). 
!            Added mass-flux option (bl_mynn_edmf, = 1 for DMP mass-flux, 0: off).
!                Related options: 
!                 bl_mynn_edmf_mom = 1 : activate momentum transport in MF scheme
!                 bl_mynn_edmf_tke = 1 : activate TKE transport in MF scheme
!            Added mixing length option (bl_mynn_mixlength, see notes below)
!            Added more sophisticated saturation checks, following Thompson scheme
!            Added new cloud PDF option (bl_mynn_cloudpdf = 2) from Chaboureau
!                and Bechtold (2002, JAS, with mods) 
!            Added capability to mix chemical species when env variable
!                WRF_CHEM = 1, thanks to Wayne Angevine.
!            Added scale-aware mixing length, following Junshi Ito's work
!                Ito et al. (2015, BLM).
! v3.9.0    Improvement to the mass-flux scheme (dynamic number of plumes,
!                better plume/cloud depth, significant speed up, better cloud
!                fraction). 
!            Added Stochastic Parameter Perturbation (SPP) implementation.
!            Many miscellaneous tweaks to the mixing lengths and stratus
!                component of the subgrid clouds.
! v.4.0      Removed or added alternatives to WRF-specific functions/modules
!                for the sake of portability to other models.
!                the sake of portability to other models.
!            Further refinement of mass-flux scheme from SCM experiments with
!                Wayne Angevine: switch to linear entrainment and back to
!                Simpson and Wiggert-type w-equation.
!            Addition of TKE production due to radiation cooling at top of 
!                clouds (proto-version); not activated by default.
!            Some code rewrites to move if-thens out of loops in an attempt to
!                improve computational efficiency.
!            New tridiagonal solver, which is supposedly 14% faster and more
!                conservative. Impact seems very small.
!            Many miscellaneous tweaks to the mixing lengths and stratus
!                component of the subgrid-scale (SGS) clouds.
! v4.1       Big improvements in downward SW radiation due to revision of subgrid clouds
!                - better cloud fraction and subgrid scale mixing ratios.
!                - may experience a small cool bias during the daytime now that high 
!                  SW-down bias is greatly reduced...
!            Some tweaks to increase the turbulent mixing during the daytime for
!                bl_mynn_mixlength option 2 to alleviate cool bias (very small impact).
!            Improved ensemble spread from changes to SPP in MYNN
!                - now perturbing eddy diffusivity and eddy viscosity directly
!                - now perturbing background rh (in SGS cloud calc only)
!                - now perturbing entrainment rates in mass-flux scheme
!            Added IF checks (within IFDEFS) to protect mixchem code from being used
!                when HRRR smoke is used (no impact on regular non-wrf chem use)
!            Important bug fix for wrf chem when transporting chemical species in MF scheme
!            Removed 2nd mass-flux scheme (no only bl_mynn_edmf = 1, no option 2)
!            Removed unused stochastic code for mass-flux scheme
!            Changed mass-flux scheme to be integrated on interface levels instead of
!                mass levels - impact is small
!            Added option to mix 2nd moments in MYNN as opposed to the scalar_pblmix option.
!                - activated with bl_mynn_mixscalars = 1; this sets scalar_pblmix = 0
!                - added tridagonal solver used in scalar_pblmix option to duplicate tendencies
!                - this alone changes the interface call considerably from v4.0.
!            Slight revision to TKE production due to radiation cooling at top of clouds
!            Added the non-Guassian buoyancy flux function of Bechtold and Siebesma (1998, JAS).
!                - improves TKE in SGS clouds
!            Added heating due to dissipation of TKE (small impact, maybe + 0.1 C daytime PBL temp)
!            Misc changes made for FV3/MPAS compatibility
! v4.2       A series of small tweaks to help reduce a cold bias in the PBL:
!                - slight increase in diffusion in convective conditions
!                - relaxed criteria for mass-flux activation/strength
!                - added capability to cycle TKE for continuity in hourly updating HRRR
!                - added effects of compensational environmental subsidence in mass-flux scheme,
!                  which resulted in tweaks to detrainment rates.
!            Bug fix for diagnostic-decay of SGS clouds - noticed by Greg Thompson. This has
!                a very small, but primarily  positive, impact on SW-down biases.
!            Tweak to calculation of KPBL - urged by Laura Fowler - to make more intuitive.
!            Tweak to temperature range of blending for saturation check (water to ice). This
!                slightly reduces excessive SGS clouds in polar region. No impact warm clouds. 
!            Added namelist option bl_mynn_output (0 or 1) to suppress or activate the
!                allocation and output of 10 3D variables. Most people will want this
!                set to 0 (default) to save memory and disk space.
!            Added new array qi_bl as opposed to using qc_bl for both SGS qc and qi. This
!                gives us more control of the magnitudes which can be confounded by using
!                a single array. As a results, many subroutines needed to be modified,
!                especially mym_condensation.
!            Added the blending of the stratus component of the SGS clouds to the mass-flux
!                clouds to account for situations where stratus and cumulus may exist in the
!                grid cell.
!            Misc small-impact bugfixes:
!                1) dz was incorrectly indexed in mym_condensation
!                2) configurations with icloud_bl = 0 were using uninitialized arrays
! v4.5 / CCPP
!            This version includes many modifications that proved valuable in the global
!            framework and removes some key lingering bugs in the mixing of chemical species.
!            TKE Budget output fixed (Puhales, 2020-12)
!            New option for stability function: (Puhales, 2020-12)
!                bl_mynn_stfunc = 0 (original, Kansas-type function, Paulson, 1970 )
!                bl_mynn_stfunc = 1 (expanded range, same as used for Jimenez et al (MWR)
!                see the Technical Note for this implementation (small impact).
!            Improved conservation of momentum and higher-order moments.
!            Important bug fixes for mixing of chemical species.
!            Addition of pressure-gradient effects on updraft momentum transport.
!            Addition of bl_mynn_closure option = 2.5, 2.6, or 3.0
!            Addition of higher-order moments for sigma when using 
!                bl_mynn_cloudpdf = 2 (Chab-Becht).
!            Removed WRF_CHEM dependencies.
!            Many miscellaneous tweaks.
! v4.6 / CCPP
!            Some code optimization. Removed many conditions from loops. Redesigned the mass-
!                flux scheme to use 8 plumes instead of a variable n plumes. This results in
!                the removal of the output variable "nudprafts" and adds maxwidth and ztop_plume.
!            Revision option bl_mynn_cloudpdf = 2, which now ensures cloud fractions for all
!                optically relevant mixing ratios (tip from Greg Thompson). Also, added flexibility
!                for tuning near-surface cloud fractions to remove excess fog/low ceilings.
!            Now outputs all SGS cloud mixing ratios as grid-mean values, not in-cloud. This 
!                results in a change in the pre-radiation code to no longer multiply mixing ratios
!                by cloud fractions.
!            Bug fix for the momentum transport.
!            Lots of code cleanup: removal of test code, comments, changing text case, etc.
!            Many misc tuning/tweaks.
!
! Many of these changes are now documented in references listed above.
!====================================================================

MODULE module_bl_mynn_driver

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

  use module_bl_mynn_1d
  use iso_c_binding !, only: c_double, c_int, c_float

  IMPLICIT NONE

!!===============================================INTERFACE_TO_C++====================================================

  interface
    SUBROUTINE  mym_predict_cc (kts,kte,                                &
       &            closure,                                            &
       &            delt,                                               &
       &            dz,                                                 &
       &            ust, flt, flq, pmz, phh,                            &
       &            el,  dfq, rho,                                      &
       &            pdk, pdt, pdq, pdc,                                 &
       &            qke, tsq, qsq, cov,                                 &
       &            s_aw,s_awqke,bl_mynn_edmf_tke,                      &
       &            qWT1D, qDISS1D,tke_budget,                          &
       &            xlvcp,xlscp) bind(c)
      import c_int,c_float

      real(kind=c_float), intent(in)    :: closure
      real(kind=c_float), intent(in)  :: xlvcp, xlscp
      integer(kind=c_int), intent(in) :: bl_mynn_edmf_tke,tke_budget
      real(kind=c_float), dimension(kts:kte), intent(in) :: dz, dfq, el, rho
      real(kind=c_float), dimension(kts:kte), intent(inout) :: pdk, pdt, pdq, pdc
      real(kind=c_float), intent(in)    :: flt, flq, pmz, phh
      real(kind=c_float), intent(in)    :: ust, delt
      real(kind=c_float), dimension(kts:kte), intent(inout) :: qke,tsq, qsq, cov
! WA 8/3/15
      real(kind=c_float), dimension(kts:kte+1), intent(inout) :: s_awqke,s_aw

      !!  TKE budget  (Puhales, 2020, WRF 4.2.1)  << EOB
      real(kind=c_float), dimension(kts:kte), intent(out) :: qWT1D, qDISS1D
    END SUBROUTINE mym_predict_cc
   end interface
      
  interface
  SUBROUTINE mynn_tendencies_cc(kts,kte,delt,dz,rho,              &
       &        u,v,tk,qv,psfc,p,thl,sqv,sqc,sqw,ust,flt,flq,flqv,   &
       &        flqc,wspd,uoce,voce,tcd,qcd,dfm,dfh,Du,Dv,        &
       &        Dth,diss_heat,s_aw,s_awthl,s_awqt,s_awqv,s_awqc,  &
       &        s_awu,s_awv,sub_thl,sub_sqv,sub_u,sub_v,    &
       &        det_thl,det_sqv,det_sqc,det_u,det_v,FLAG_QC_I,    &
       &        bl_mynn_cloudmix,bl_mynn_mixqt,bl_mynn_edmf,      &
       &        bl_mynn_edmf_mom,debug_code,r_d, p608, ep_2,      &
       &        ep_3,tv0,xlv,xlvcp) bind(c)
     import c_int,c_float
    integer(kind=c_int), intent(in) :: kts,kte
    integer(kind=c_int), intent(in) :: debug_code, FLAG_QC_I
    integer(kind=c_int), intent(in) :: bl_mynn_cloudmix,bl_mynn_mixqt,  &
                           bl_mynn_edmf,bl_mynn_edmf_mom               
    real(kind=c_float), dimension(kts:kte+1), intent(in) :: s_awthl,    &
         &s_aw,s_awqt,s_awqv,s_awqc,s_awu,s_awv 
    real(kind=c_float) :: ep_2,ep_3,tv0,xlv,xlvcp
! tendencies from mass-flux environmental subsidence and detrainment
    real(kind=c_float), dimension(kts:kte), intent(in) :: sub_thl,sub_sqv,   &
         &sub_u,sub_v,det_thl,det_sqv,det_sqc,det_u,det_v
    real(kind=c_float), dimension(kts:kte), intent(in) ::u,v,tk,rho,p,dz,   &
    & qv,tcd,qcd,diss_heat
    real(kind=c_float), dimension(kts:kte), intent(inout) :: thl,sqw,sqv,sqc,&
         &dfm,dfh
    real(kind=c_float), dimension(kts:kte), intent(inout) ::Du,Dv,Dth
    real(kind=c_float), intent(in) :: flt,flq,flqv,flqc,uoce,voce
    real(kind=c_float), intent(in) :: ust,delt,psfc,wspd,r_d, p608
  END SUBROUTINE mynn_tendencies_cc
  end interface

interface
    SUBROUTINE mynn_mix_chem_cc (kts,kte,i,        &
     &          delt,dz,pblh,                      &
     &          nchem, kdvel, ndvel,               &
     &          chem1, vd1,                        &
     &          rho,                               &
     &          flt, tcd, qcd,                     &
     &          dfh,                               &
     &          s_aw, s_awchem,                    &
     &          emis_ant_no, frp, rrfs_sd,         &
     &          enh_mix                            ) bind(c)
     import c_int,c_float
    integer(kind=c_int), intent(in) :: kts,kte,i
    real(kind=c_float), dimension(kts:kte), intent(in)    :: dfh,dz,tcd,qcd
    real(kind=c_float), dimension(kts:kte), intent(inout) :: rho
    real(kind=c_float), intent(in)    :: flt
    real(kind=c_float), intent(in)    :: delt,pblh
    integer(kind=c_int), intent(in) :: nchem, kdvel, ndvel
    real(kind=c_float), dimension( kts:kte+1), intent(in) :: s_aw
    real(kind=c_float), dimension( kts:kte, nchem ), intent(inout) :: chem1
    real(kind=c_float), dimension( kts:kte+1,nchem), intent(in) :: s_awchem
    real(kind=c_float), dimension( ndvel ), intent(in) :: vd1
    real(kind=c_float), intent(in) :: emis_ant_no,frp
    integer(kind=c_int), intent(in) :: rrfs_sd, enh_mix
    END SUBROUTINE mynn_mix_chem_cc
  end interface

  interface
  SUBROUTINE moisture_check_cc(kte, delt, dp, exner,   &
                            qv, qc, qi, qs, th,     &
                            dqv, dqc, dqi, dqs, dth ) bind (c)
  !! This subroutine was adopted from the CAM-UW ShCu scheme and
  !! adapted for use here.
  !!
  !! If qc < qcmin, qi < qimin, or qv < qvmin happens in any layer,
  !! force them to be larger than minimum value by (1) condensating
  !! water vapor into liquid or ice, and (2) by transporting water vapor
  !! from the very lower layer.
  !!
  !! We then update the final state variables and tendencies associated
  !! with this correction. If any condensation happens, update theta too.
  !! Note that (qv,qc,qi,th) are the final state variables after
  !! applying corresponding input tendencies and corrective tendencies.
    import c_int, c_float
    integer(kind=c_int),         intent(in)     :: kte
    real(kind=c_float), intent(in)     :: delt
    real(kind=c_float), dimension(kte), intent(in)     :: dp, exner
    real(kind=c_float), dimension(kte), intent(inout)  :: qv, qc, qi, qs, th
    real(kind=c_float), dimension(kte), intent(inout)  :: dqv, dqc, dqi, dqs, dth
  END SUBROUTINE moisture_check_cc
  end interface

interface
    SUBROUTINE  mym_condensation_cc (kts,kte,   &
    &            dx, dz, zw, xland,        &
    &            thl, qw, qv, qc, qi, qs,  &
    &            p,exner,                  &
    &            tsq, qsq, cov,            &
    &            Sh, el, bl_mynn_cloudpdf, &
    &            qc_bl1D, qi_bl1D,         &
    &            cldfra_bl1D,              &
    &            PBLH1,HFX1,               &
    &            Vt, Vq, th, sgm, rmo,     &
    &            spp_pbl,rstoch_col,ep_2,  &
    &            ep_3,xlv,r_d,xlvcp,p608,  &
    &            tv0, cpv, r_v, cice,cliq, & 
    &            cp,xls,rcp) bind(c)
!-------------------------------------------------------------------
    import c_int, c_float

    integer(kind=c_int), intent(in) :: kts,kte,bl_mynn_cloudpdf,spp_pbl
    real(kind=c_float), intent(in)    :: ep_2, ep_3,xlv, r_d, xlvcp, &
                  p608, tv0,cpv,r_v,cice,cliq,cp,xls 
    real(kind=c_float), intent(in)      :: HFX1,rmo,xland
    real(kind=c_float), intent(in)      :: dx,pblh1
    real(kind=c_float), dimension(kts:kte), intent(in) :: dz
    real(kind=c_float), dimension(kts:kte+1), intent(in) :: zw
    real(kind=c_float), dimension(kts:kte), intent(in) :: p,exner,thl,qw,   &
         &qv,qc,qi,qs,tsq,qsq,cov,th,rstoch_col
    real(kind=c_float), dimension(kts:kte), intent(inout) :: vt,vq,sgm
    real(kind=c_float), dimension(kts:kte) :: alp,a,bet,b,ql,q1,RH
    real(kind=c_float), dimension(kts:kte), intent(out) :: qc_bl1D,qi_bl1D, &
         &cldfra_bl1D
    real(kind=c_float), dimension(kts:kte), intent(in) :: Sh,el
 END SUBROUTINE mym_condensation_cc
 end interface 

interface
 SUBROUTINE topdown_cloudrad_cc(kts,kte,                         &
               &dz1,zw,fltv,xland,kpbl,PBLH,                  &
               &sqc,sqi,sqw,thl,th1,ex1,p1,rho1,thetav,       &
               &cldfra_bl1D,rthraten,                         &
               &maxKHtopdown,KHtopdown,TKEprodTD) bind(c)
    import c_int, c_float
    real(kind=c_float)                    :: maxKHtopdown
    real(kind=c_float),dimension(kts:kte) :: KHtopdown,TKEprodTD
    integer(kind=c_int),intent(in)        :: kte,kts
    real(kind=c_float), dimension(kts:kte), intent(in)   :: dz1,sqc,sqi,sqw,&
          thl,th1,ex1,p1,rho1,thetav,cldfra_bl1D
    real(kind=c_float), dimension(kts:kte), intent(in)   :: rthraten
    real(kind=c_float), dimension(kts:kte+1), intent(in) :: zw
    real(kind=c_float), intent(in)                       :: pblh,fltv
    real(kind=c_float), intent(in)                       :: xland
    integer(kind=c_int), intent(in)                      :: kpbl
  END SUBROUTINE topdown_cloudrad_cc
 end interface

  interface
  SUBROUTINE DDMF_JPL_cc(kts,kte,dt,zw,dz,p,         &
              &u,v,th,thl,thv,tk,qt,qv,qc,           &
              &rho,exner,ust,wthl,wqt,pblh,kpbl,     &
              &edmf_a_dd,edmf_w_dd, edmf_qt_dd,      &
              &edmf_thl_dd,edmf_ent_dd,edmf_qc_dd,   &
              &sd_aw,sd_awthl,sd_awqt,               &
              &sd_awqv,sd_awqc,sd_awu,sd_awv,        &
              &sd_awqke,qc_bl1d,cldfra_bl1d,         &
              &rthraten,svp1,grav,onethird,          &
              & p1000mb,rcp,xlvcp                   ) bind(c)
        import c_int, c_float
        integer(kind=c_int), intent(in) :: kts,kte,kpbl
        real(kind=c_float), dimension(kts:kte), intent(in) ::            &
            u,v,th,thl,tk,qt,qv,qc,thv,p,rho,exner,dz
        real(kind=c_float), dimension(kts:kte), intent(in) :: rthraten
        real(kind=c_float), dimension(kts:kte), intent(in) :: edmf_a_dd,&
                edmf_w_dd,edmf_thl_dd,edmf_ent_dd,edmf_qc_dd,edmf_qt_dd 
        real(kind=c_float), dimension(kts:kte), intent(in) ::sd_aw, &
                sd_awthl,sd_awqt,sd_awqv,sd_awqc,sd_awu,sd_awv,sd_awqke
        real(kind=c_float), dimension(kts:kte), intent(in) :: qc_bl1d,cldfra_bl1d 
        ! zw .. heights of the downdraft levels (edges of boxes)
        real(kind=c_float), dimension(kts:kte+1), intent(in) :: zw
        real(kind=c_float), intent(in)  :: dt,ust,pblh,wthl,wqt
        real(kind=c_float), intent(in)  :: svp1,grav,onethird,p1000mb,rcp,xlvcp
      END SUBROUTINE DDMF_JPL_cc
  end interface

  interface
    SUBROUTINE scale_aware_cc(dx,PBL1,Psig_bl,Psig_shcu) bind(c)
        import c_float
        real(kind=c_float), intent(in)  :: dx,pbl1
        real(kind=c_float), intent(out) :: Psig_bl,Psig_shcu
    END SUBROUTINE scale_aware_cc
  end interface

  interface
    SUBROUTINE GET_PBLH_cc(KTS,KTE,zi,thetav1D,qke1D,zw1D, &
      dz1D,landsea,kzi) bind(c)
      import c_float, c_int
      integer(kind=c_int),intent(in) :: KTS,KTE
      real(kind=c_float), intent(in) :: landsea
      real(kind=c_float), dimension(kts:kte), intent(in) :: thetav1D, qke1D, dz1D
      real(kind=c_float), dimension(kts:kte+1), intent(in) :: zw1D
    END SUBROUTINE GET_PBLH_cc
  end interface

 interface
   SUBROUTINE retrieve_exchange_coeffs_cc(kts,kte,&
       &dfm,dfh,dz,K_m,K_h) bind(c)
        import c_float, c_int
        integer(kind=c_int) , intent(in) :: kts,kte
        real(kind=c_float), dimension(KtS:KtE), intent(in) :: dz,dfm,dfh
    END SUBROUTINE retrieve_exchange_coeffs_cc
 end interface


!!===============================================INTERFACES_END====================================================

  !>Use Ito et al. (2015, BLM) scale-aware (0: no, 1: yes). Note that this also has impacts
  !!on the cloud PDF and mass-flux scheme, using LES-derived similarity function.
  real(kind_phys), parameter :: scaleaware=1.

  !>Of the following the options, use one OR the other, not both.
  !>Adding top-down diffusion driven by cloud-top radiative cooling
  integer, parameter :: bl_mynn_topdown = 0
  !>Option to activate downdrafts, from Elynn Wu (0: deactive, 1: active)
  integer, parameter :: bl_mynn_edmf_dd = 0

  !>Option to activate heating due to dissipation of TKE (to activate, set to 1.0)
  integer, parameter :: dheat_opt = 1

  !Option to switch flux-profile relationship for surface (from Puhales et al. 2020)
  !0: use original Dyer-Hicks, 1: use Cheng-Brustaert and Blended COARE
  integer, parameter :: bl_mynn_stfunc = 1

  !option to print out more stuff for debugging purposes
  logical, parameter :: debug_code = .false.
  integer, parameter :: debug_code_i = 0
  integer, parameter :: idbg = 23 !specific i-point to write out

  ! Used in WRF-ARW module_physics_init.F
  integer :: mynn_level


CONTAINS

! ==================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine is the GSD MYNN-EDNF PBL driver routine,which
!! encompassed the majority of the subroutines that comprise the 
!! procedures that ultimately solve for tendencies of 
!! \f$U, V, \theta, q_v, q_c, and q_i\f$.
!!\section gen_mynn_bl_driver GSD mynn_bl_driver General Algorithm
!> @{
  SUBROUTINE mynn_bl_driver(            &
       &initflag,restart,cycling,       &
       &delt,dz,dx,znt,                 &
       &u,v,w,th,sqv3d,sqc3d,sqi3d,     &
       &sqs3d,qnc,qni,                  &
       &qnwfa,qnifa,qnbca,ozone,        &
       &p,exner,rho,t3d,                &
       &xland,ts,qsfc,ps,               &
       &ust,ch,hfx,qfx,rmol,wspd,       &
       &uoce,voce,                      & !ocean current
       &qke,qke_adv,                    &
       &sh3d,sm3d,                      &
       &nchem,kdvel,ndvel,              & !smoke/chem variables
       &chem3d,vdep,                    &
       &frp,emis_ant_no,                &
       &mix_chem,enh_mix,               & !note: these arrays/flags are still under development
       &rrfs_sd,smoke_dbg,              & !end smoke/chem variables
       &tsq,qsq,cov,                    &
       &rublten,rvblten,rthblten,       &
       &rqvblten,rqcblten,rqiblten,     &
       &rqncblten,rqniblten,rqsblten,   &
       &rqnwfablten,rqnifablten,        &
       &rqnbcablten,dozone,             &
       &exch_h,exch_m,                  &
       &pblh,kpbl,                      & 
       &el_pbl,                         &
       &dqke,qwt,qshear,qbuoy,qdiss,    &
       &qc_bl,qi_bl,cldfra_bl,          &
       &bl_mynn_tkeadvect,              &
       &tke_budget,                     &
       &bl_mynn_cloudpdf,               &
       &bl_mynn_mixlength,              &
       &icloud_bl,                      &
       &closure,                        &
       &bl_mynn_edmf,                   &
       &bl_mynn_edmf_mom,               &
       &bl_mynn_edmf_tke,               &
       &bl_mynn_mixscalars,             &
       &bl_mynn_output,                 &
       &bl_mynn_cloudmix,bl_mynn_mixqt, &
       &edmf_a,edmf_w,edmf_qt,          &
       &edmf_thl,edmf_ent,edmf_qc,      &
       &sub_thl3D,sub_sqv3D,            &
       &det_thl3D,det_sqv3D,            &
       &maxwidth,maxMF,ztop_plume,      &
       &ktop_plume,                     &
       &spp_pbl,pattern_spp_pbl,        &
       &rthraten,                       &
       &FLAG_QC,FLAG_QI,FLAG_QNC,       &
       &FLAG_QNI,FLAG_QS,               &
       &FLAG_QNWFA,FLAG_QNIFA,          &
       &FLAG_QNBCA,FLAG_OZONE,          &
       &IDS,IDE,JDS,JDE,KDS,KDE,        &
       &IMS,IME,JMS,JME,KMS,KME,        &
       &ITS,ITE,JTS,JTE,KTS,KTE         )
    
!-------------------------------------------------------------------

    integer, intent(in) :: initflag
    !INPUT NAMELIST OPTIONS:
    logical, intent(in) :: restart,cycling
    integer, intent(in) :: tke_budget
    integer, intent(in) :: bl_mynn_cloudpdf
    integer, intent(in) :: bl_mynn_mixlength
    integer, intent(in) :: bl_mynn_edmf
    logical, intent(in) :: bl_mynn_tkeadvect
    integer, intent(in) :: bl_mynn_edmf_mom
    integer, intent(in) :: bl_mynn_edmf_tke
    integer, intent(in) :: bl_mynn_mixscalars
    integer, intent(in) :: bl_mynn_output
    integer, intent(in) :: bl_mynn_cloudmix
    integer, intent(in) :: bl_mynn_mixqt
    integer, intent(in) :: icloud_bl
    real(kind_phys), intent(in) :: closure

    logical, intent(in) :: FLAG_QI,FLAG_QNI,FLAG_QC,FLAG_QNC,&
                           FLAG_QNWFA,FLAG_QNIFA,FLAG_QNBCA, &
                           FLAG_OZONE,FLAG_QS

    logical, intent(in) :: mix_chem,enh_mix,rrfs_sd,smoke_dbg

    integer, intent(in) ::                                   &
                         & IDS,IDE,JDS,JDE,KDS,KDE           &
                         &,IMS,IME,JMS,JME,KMS,KME           &
                         &,ITS,ITE,JTS,JTE,KTS,KTE

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

! initflag > 0  for TRUE
! else        for FALSE
!       closure       : <= 2.5;  Level 2.5
!                  2.5< and <3;  Level 2.6
!                        =   3;  Level 3

! SGT: Changed this to use assumed shape arrays (dimension(:,:,:)) with no "optional" arguments
!      to prevent a crash on Cheyenne. Do not change it back without testing if the code runs
!      on Cheyenne with the GNU compiler.
    integer :: enh_mix_in 
    real(kind_phys), intent(in) :: delt
    real(kind_phys), dimension(ims:ime),   intent(in) :: dx
    real(kind_phys), dimension(ims:ime,kms:kme), intent(in) :: dz,              &
         &u,v,w,th,sqv3D,p,exner,rho,T3D
    real(kind_phys), dimension(ims:ime,kms:kme), optional, intent(in) ::        &
         &sqc3D,sqi3D,sqs3D,qni,qnc,qnwfa,qnifa,qnbca
    real(kind_phys), dimension(ims:ime,kms:kme), optional,intent(in):: ozone
    real(kind_phys), dimension(ims:ime),   intent(in):: ust,                    &
         &ch,qsfc,ps,wspd
    real(kind_phys), dimension(ims:ime,kms:kme), intent(inout) ::               &
         &Qke,Tsq,Qsq,Cov,qke_adv
    real(kind_phys), dimension(ims:ime,kms:kme), intent(inout) ::               &
         &rublten,rvblten,rthblten,rqvblten,rqcblten,                           &
         &rqiblten,rqsblten,rqniblten,rqncblten,                                &
         &rqnwfablten,rqnifablten,rqnbcablten
    real(kind_phys), dimension(ims:ime,kms:kme), intent(inout) :: dozone
    real(kind_phys), dimension(ims:ime,kms:kme), intent(in)    :: rthraten

    real(kind_phys), dimension(ims:ime,kms:kme), intent(out)   :: exch_h,exch_m
    real(kind_phys), dimension(ims:ime),   intent(in)    :: xland,              &
         &ts,znt,hfx,qfx,uoce,voce

   !These 10 arrays are only allocated when bl_mynn_output > 0
   real(kind_phys), dimension(ims:ime,kms:kme), optional, intent(inout) ::      &
         & edmf_a,edmf_w,edmf_qt,edmf_thl,edmf_ent,edmf_qc,                     &
         & sub_thl3D,sub_sqv3D,det_thl3D,det_sqv3D

!   real, dimension(ims:ime,kms:kme)   ::                                       &
!         & edmf_a_dd,edmf_w_dd,edmf_qt_dd,edmf_thl_dd,edmf_ent_dd,edmf_qc_dd

    real(kind_phys), dimension(ims:ime), intent(inout) :: pblh
    real(kind_phys), dimension(ims:ime), intent(inout) :: rmol

    real(kind_phys), dimension(ims:ime) :: psig_bl,psig_shcu

    integer,dimension(ims:ime),intent(inout) ::                                 &
         &KPBL,ktop_plume

    real(kind_phys), dimension(ims:ime), intent(out) ::                         &
         &maxmf,maxwidth,ztop_plume

    real(kind_phys), dimension(ims:ime,kms:kme), intent(inout) :: el_pbl

    real(kind_phys), dimension(ims:ime,kms:kme), optional, intent(inout) ::     &
         &qWT,qSHEAR,qBUOY,qDISS,dqke
    ! 3D budget arrays are not allocated when tke_budget == 0
    ! 1D (local) budget arrays are used for passing between subroutines.
    real(kind_phys), dimension(kts:kte) ::                                      &
         &qwt1,qshear1,qbuoy1,qdiss1,dqke1,diss_heat

    real(kind_phys), dimension(ims:ime,kms:kme), intent(out) :: Sh3D,Sm3D

    real(kind_phys), dimension(ims:ime,kms:kme), intent(inout) ::               &
         &qc_bl,qi_bl,cldfra_bl
    real(kind_phys), dimension(kts:kte) :: qc_bl1D,qi_bl1D,                     &
         &cldfra_bl1D,qc_bl1D_old,qi_bl1D_old,cldfra_bl1D_old

! smoke/chemical arrays
    integer, intent(IN   ) ::   nchem, kdvel, ndvel
    real(kind_phys), dimension(ims:ime,kms:kme,nchem), optional, intent(inout) :: chem3d
    real(kind_phys), dimension(ims:ime, ndvel), optional,  intent(in)    :: vdep
    real(kind_phys), dimension(ims:ime),   optional,   intent(in)    :: frp,EMIS_ANT_NO
    !local
    real(kind_phys), dimension(kts:kte  ,nchem)      :: chem1
    real(kind_phys), dimension(kts:kte+1,nchem)      :: s_awchem1
    real(kind_phys), dimension(ndvel)                :: vd1
    integer :: ic

!local vars
    integer :: ITF,JTF,KTF, IMD,JMD
    integer :: i,j,k,kproblem
    real(kind_phys), dimension(kts:kte) ::                  &
         &thl,tl,qv1,qc1,qi1,qs1,sqw,                       &
         &el, dfm, dfh, dfq, tcd, qcd, pdk, pdt, pdq, pdc,  &
         &vt, vq, sgm, kzero
    real(kind_phys), dimension(kts:kte) ::                  &
         &thetav,sh,sm,u1,v1,w1,p1,                         &
         &ex1,dz1,th1,tk1,rho1,qke1,tsq1,qsq1,cov1,         &
         &sqv,sqi,sqc,sqs,                                  &
         &du1,dv1,dth1,dqv1,dqc1,dqi1,dqs1,ozone1,          &
         &k_m1,k_h1,qni1,dqni1,qnc1,dqnc1,qnwfa1,qnifa1,    &
         &qnbca1,dqnwfa1,dqnifa1,dqnbca1,dozone1

    !mass-flux variables
    real(kind_phys), dimension(kts:kte) ::                  &
         &dth1mf,dqv1mf,dqc1mf,du1mf,dv1mf
    real(kind_phys), dimension(kts:kte) ::                  &
         &edmf_a1,edmf_w1,edmf_qt1,edmf_thl1,               &
         &edmf_ent1,edmf_qc1
    real(kind_phys), dimension(kts:kte) ::                  &
         &edmf_a_dd1,edmf_w_dd1,edmf_qt_dd1,edmf_thl_dd1,   &
         &edmf_ent_dd1,edmf_qc_dd1
    real(kind_phys), dimension(kts:kte) ::                  &
         &sub_thl,sub_sqv,sub_u,sub_v,                      &
         &det_thl,det_sqv,det_sqc,det_u,det_v
    real(kind_phys), dimension(kts:kte+1) ::                &
         &s_aw1,s_awthl1,s_awqt1,                           &
         &s_awqv1,s_awqc1,s_awu1,s_awv1,s_awqke1,           &
         &s_awqnc1,s_awqni1,s_awqnwfa1,s_awqnifa1,          &
         &s_awqnbca1
    real(kind_phys), dimension(kts:kte+1) ::                &
         &sd_aw1,sd_awthl1,sd_awqt1,                        &
         &sd_awqv1,sd_awqc1,sd_awu1,sd_awv1,sd_awqke1

    real(kind_phys), dimension(kts:kte+1) :: zw
    real(kind_phys) :: cpm,sqcg,flt,fltv,flq,flqv,flqc,     &
         &pmz,phh,exnerg,zet,phi_m,                         &
         &afk,abk,ts_decay, qc_bl2, qi_bl2,                 &
         &th_sfc,wsp

    !top-down diffusion
    real(kind_phys), dimension(ITS:ITE) :: maxKHtopdown
    real(kind_phys), dimension(kts:kte) :: KHtopdown,TKEprodTD

    logical :: INITIALIZE_QKE,problem

    ! Stochastic fields 
    integer,  intent(in)                         :: spp_pbl
    real(kind_phys), dimension(ims:ime,kms:kme), optional, intent(in)  :: pattern_spp_pbl
    real(kind_phys), dimension(kts:kte)          :: rstoch_col


    ! Substepping TKE
    integer :: nsub
    real(kind_phys) :: delt2

    integer :: FLAG_QC_I
    
    if (FLAG_QC .eqv. .true.) then
      FLAG_QC_I=1
    else
      FLAG_QC_I=0
    endif


    write(*,*) "Calling mynn"
   !! CALL helloworld()
    if (debug_code) then !check incoming values
      do i=its,ite
        problem = .false.
        do k=kts,kte
          wsp  = sqrt(u(i,k)**2 + v(i,k)**2)
          if (abs(hfx(i)) > 1200. .or. abs(qfx(i)) > 0.001 .or.         &
              wsp > 200. .or. t3d(i,k) > 360. .or. t3d(i,k) < 160. .or. &
              sqv3d(i,k)< 0.0 .or. sqc3d(i,k)< 0.0 ) then
             kproblem = k
             problem = .true.
             print*,"Incoming problem at: i=",i," k=1"
             print*," QFX=",qfx(i)," HFX=",hfx(i)
             print*," wsp=",wsp," T=",t3d(i,k)
             print*," qv=",sqv3d(i,k)," qc=",sqc3d(i,k)
             print*," u*=",ust(i)," wspd=",wspd(i)
             print*," xland=",xland(i)," ts=",ts(i)
             print*," z/L=",0.5*dz(i,1)*rmol(i)," ps=",ps(i)
             print*," znt=",znt(i)," dx=",dx(i)
          endif
        enddo
        if (problem) then
          print*,"===tk:",t3d(i,max(kproblem-3,1):min(kproblem+3,kte))
          print*,"===qv:",sqv3d(i,max(kproblem-3,1):min(kproblem+3,kte))
          print*,"===qc:",sqc3d(i,max(kproblem-3,1):min(kproblem+3,kte))
          print*,"===qi:",sqi3d(i,max(kproblem-3,1):min(kproblem+3,kte))
          print*,"====u:",u(i,max(kproblem-3,1):min(kproblem+3,kte))
          print*,"====v:",v(i,max(kproblem-3,1):min(kproblem+3,kte))
        endif
      enddo
    endif

!***  Begin debugging
    IMD=(IMS+IME)/2
    JMD=(JMS+JME)/2
!***  End debugging

    JTF=JTE
    ITF=ITE
    KTF=KTE

    IF (bl_mynn_output > 0) THEN !research mode
       edmf_a(its:ite,kts:kte)=0.
       edmf_w(its:ite,kts:kte)=0.
       edmf_qt(its:ite,kts:kte)=0.
       edmf_thl(its:ite,kts:kte)=0.
       edmf_ent(its:ite,kts:kte)=0.
       edmf_qc(its:ite,kts:kte)=0.
       sub_thl3D(its:ite,kts:kte)=0.
       sub_sqv3D(its:ite,kts:kte)=0.
       det_thl3D(its:ite,kts:kte)=0.
       det_sqv3D(its:ite,kts:kte)=0.

       !edmf_a_dd(its:ite,kts:kte)=0.
       !edmf_w_dd(its:ite,kts:kte)=0.
       !edmf_qt_dd(its:ite,kts:kte)=0.
       !edmf_thl_dd(its:ite,kts:kte)=0.
       !edmf_ent_dd(its:ite,kts:kte)=0.
       !edmf_qc_dd(its:ite,kts:kte)=0.
    ENDIF
    ktop_plume(its:ite)=0   !int
    ztop_plume(its:ite)=0.
    maxwidth(its:ite)=0.
    maxmf(its:ite)=0.
    maxKHtopdown(its:ite)=0.
    kzero(kts:kte)=0.

    ! DH* CHECK HOW MUCH OF THIS INIT IF-BLOCK IS ACTUALLY NEEDED FOR RESTARTS
!> - Within the MYNN-EDMF, there is a dependecy check for the first time step,
!! If true, a three-dimensional initialization loop is entered. Within this loop,
!! several arrays are initialized and k-oriented (vertical) subroutines are called 
!! at every i and j point, corresponding to the x- and y- directions, respectively.  
    IF (initflag > 0 .and. .not.restart) THEN

       !Test to see if we want to initialize qke
       IF ( (restart .or. cycling)) THEN
          IF (MAXVAL(QKE(its:ite,kts)) < 0.0002) THEN
             INITIALIZE_QKE = .TRUE.
             !print*,"QKE is too small, must initialize"
          ELSE
             INITIALIZE_QKE = .FALSE.
             !print*,"Using background QKE, will not initialize"
          ENDIF
       ELSE ! not cycling or restarting:
          INITIALIZE_QKE = .TRUE.
          !print*,"not restart nor cycling, must initialize QKE"
       ENDIF
 
       if (.not.restart .or. .not.cycling) THEN
         Sh3D(its:ite,kts:kte)=0.
         Sm3D(its:ite,kts:kte)=0.
         el_pbl(its:ite,kts:kte)=0.
         tsq(its:ite,kts:kte)=0.
         qsq(its:ite,kts:kte)=0.
         cov(its:ite,kts:kte)=0.
         cldfra_bl(its:ite,kts:kte)=0.
         qc_bl(its:ite,kts:kte)=0.
         qke(its:ite,kts:kte)=0.
       else
         qc_bl1D(kts:kte)=0.0
         qi_bl1D(kts:kte)=0.0
         cldfra_bl1D(kts:kte)=0.0
       end if
       dqc1(kts:kte)=0.0
       dqi1(kts:kte)=0.0
       dqni1(kts:kte)=0.0
       dqnc1(kts:kte)=0.0
       dqnwfa1(kts:kte)=0.0
       dqnifa1(kts:kte)=0.0
       dqnbca1(kts:kte)=0.0
       dozone1(kts:kte)=0.0
       qc_bl1D_old(kts:kte)=0.0
       cldfra_bl1D_old(kts:kte)=0.0
       edmf_a1(kts:kte)=0.0
       edmf_w1(kts:kte)=0.0
       edmf_qc1(kts:kte)=0.0
       edmf_a_dd1(kts:kte)=0.0
       edmf_w_dd1(kts:kte)=0.0
       edmf_qc_dd1(kts:kte)=0.0
       sgm(kts:kte)=0.0
       vt(kts:kte)=0.0
       vq(kts:kte)=0.0

       DO k=KTS,KTE
          DO i=ITS,ITF
             exch_m(i,k)=0.
             exch_h(i,k)=0.
          ENDDO
       ENDDO

       IF (tke_budget .eq. 1) THEN
          DO k=KTS,KTE
             DO i=ITS,ITF
                qWT(i,k)=0.
                qSHEAR(i,k)=0.
                qBUOY(i,k)=0.
                qDISS(i,k)=0.
                dqke(i,k)=0.
             ENDDO
          ENDDO
       ENDIF

       DO i=ITS,ITF !this is where we go from 2d to 1d
          if (FLAG_QI ) then
             sqi(:)=sqi3D(i,:)
          else
             sqi = 0.0
          endif
          if (FLAG_QS ) then
             sqs(:)=sqs3D(i,:)
          else
             sqs = 0.0
          endif
          if (icloud_bl > 0) then
             cldfra_bl1d(:)=cldfra_bl(i,:)
             qc_bl1d(:)=qc_bl(i,:)
             qi_bl1d(:)=qi_bl(i,:)
          endif

          do k=KTS,KTE !KTF
                dz1(k)=dz(i,k)
                u1(k) = u(i,k)
                v1(k) = v(i,k)
                w1(k) = w(i,k)
                th1(k)=th(i,k)
                tk1(k)=T3D(i,k)
                ex1(k)=exner(i,k)
                rho1(k)=rho(i,k)
                sqc(k)=sqc3D(i,k) !/(1.+qv(i,k))
                sqv(k)=sqv3D(i,k) !/(1.+qv(i,k))
                thetav(k)=th(i,k)*(1.+p608*sqv(k))
                !keep snow out for now - increases ceiling bias
                sqw(k)=sqv(k)+sqc(k)+sqi(k)!+sqs(k)
                thl(k)=th1(k) - xlvcp/ex1(k)*sqc(k) &
                    &         - xlscp/ex1(k)*(sqi(k))!+sqs(k))
                !Use form from Tripoli and Cotton (1981) with their
                !suggested min temperature to improve accuracy.
                !thl(k)=th(i,k)*(1.- xlvcp/MAX(tk1(k),TKmin)*sqc(k) &
                !    &               - xlscp/MAX(tk1(k),TKmin)*sqi(k))

                IF (k==kts) THEN
                   zw(k)=0.
                ELSE
                   zw(k)=zw(k-1)+dz(i,k-1)
                ENDIF
                IF (INITIALIZE_QKE) THEN
                   !Initialize tke for initial PBLH calc only - using 
                   !simple PBLH form of Koracin and Berkowicz (1988, BLM)
                   !to linearly taper off tke towards top of PBL.
                   qke1(k)=5.*ust(i) * MAX((ust(i)*700. - zw(k))/(MAX(ust(i),0.01)*700.), 0.01)
                ELSE
                   qke1(k)=qke(i,k)
                ENDIF
                el(k)=el_pbl(i,k)
                sh(k)=Sh3D(i,k)
                sm(k)=Sm3D(i,k)
                tsq1(k)=tsq(i,k)
                qsq1(k)=qsq(i,k)
                cov1(k)=cov(i,k)
                if (spp_pbl==1) then
                    rstoch_col(k)=pattern_spp_pbl(i,k)
                else
                    rstoch_col(k)=0.0
                endif

             ENDDO

             zw(kte+1)=zw(kte)+dz(i,kte)

!>  - Call get_pblh() to calculate hybrid (\f$\theta_{v}-TKE\f$) PBL height.
             CALL GET_PBLH_cc(KTS,KTE,PBLH(i),thetav,&
               &  Qke1,zw,dz1,xland(i),KPBL(i))
             
!>  - Call scale_aware() to calculate similarity functions for scale-adaptive control
!! (\f$P_{\sigma-PBL}\f$ and \f$P_{\sigma-shcu}\f$).
             IF (scaleaware > 0.) THEN
                CALL scale_aware_cc(dx(i),PBLH(i),Psig_bl(i),Psig_shcu(i))
             ELSE
                Psig_bl(i)=1.0
                Psig_shcu(i)=1.0
             ENDIF

             ! DH* CHECK IF WE CAN DO WITHOUT CALLING THIS ROUTINE FOR RESTARTS
!>  - Call mym_initialize() to initializes the mixing length, TKE, \f$\theta^{'2}\f$,
!! \f$q^{'2}\f$, and \f$\theta^{'}q^{'}\f$. These variables are calculated after 
!! obtaining prerequisite variables by calling the following subroutines from 
!! within mym_initialize(): mym_level2() and mym_length().
             CALL mym_initialize (                & 
                  &kts,kte,xland(i),              &
                  &dz1, dx(i), zw,                &
                  &u1, v1, thl, sqv,              &
                  &PBLH(i), th1, thetav, sh, sm,  &
                  &ust(i), rmol(i),               &
                  &el, Qke1, Tsq1, Qsq1, Cov1,    &
                  &Psig_bl(i), cldfra_bl1D,       &
                  &bl_mynn_mixlength,             &
                  &edmf_w1,edmf_a1,               &
                  &INITIALIZE_QKE,                &
                  &spp_pbl,rstoch_col             )

             IF (.not.restart) THEN
                !UPDATE 3D VARIABLES
                DO k=KTS,KTE !KTF
                   el_pbl(i,k)=el(k)
                   sh3d(i,k)=sh(k)
                   sm3d(i,k)=sm(k)
                   qke(i,k)=qke1(k)
                   tsq(i,k)=tsq1(k)
                   qsq(i,k)=qsq1(k)
                   cov(i,k)=cov1(k)
                ENDDO
                !initialize qke_adv array if using advection
                IF (bl_mynn_tkeadvect) THEN
                   DO k=KTS,KTE
                      qke_adv(i,k)=qke1(k)
                   ENDDO
                ENDIF
             ENDIF

!***  Begin debugging
!             IF(I==IMD .AND. J==JMD)THEN
!               PRINT*,"MYNN DRIVER INIT: k=",1," sh=",sh(k)
!               PRINT*," sqw=",sqw(k)," thl=",thl(k)," k_m=",exch_m(i,k)
!               PRINT*," xland=",xland(i)," rmol=",rmol(i)," ust=",ust(i)
!               PRINT*," qke=",qke(i,k)," el=",el_pbl(i,k)," tsq=",Tsq(i,k)
!               PRINT*," PBLH=",PBLH(i)," u=",u(i,k)," v=",v(i,k)
!             ENDIF
!***  End debugging

       ENDDO !end i-loop

    ENDIF ! end initflag

!> - After initializing all required variables, the regular procedures 
!! performed at every time step are ready for execution.
    !ACF- copy qke_adv array into qke if using advection
    IF (bl_mynn_tkeadvect) THEN
       qke=qke_adv
    ENDIF

    DO i=ITS,ITF
       !Initialize some arrays
       if (tke_budget .eq. 1) then
          dqke(i,:)=qke(i,:)
       endif
       if (FLAG_QI ) then
          sqi(:)=sqi3D(i,:)
       else
          sqi = 0.0
       endif
       if (FLAG_QS ) then
          sqs(:)=sqs3D(i,:)
       else
          sqs = 0.0
       endif
       if (icloud_bl > 0) then
          CLDFRA_BL1D(:)=CLDFRA_BL(i,:)
          QC_BL1D(:)   =QC_BL(i,:)
          QI_BL1D(:)   =QI_BL(i,:)
          cldfra_bl1D_old(:)=cldfra_bl(i,:)
          qc_bl1D_old(:)=qc_bl(i,:)
          qi_bl1D_old(:)=qi_bl(i,:)
       else
          CLDFRA_BL1D  =0.0
          QC_BL1D      =0.0
          QI_BL1D      =0.0
          cldfra_bl1D_old=0.0
          qc_bl1D_old  =0.0
          qi_bl1D_old  =0.0
       endif
       dz1(kts:kte)    =dz(i,kts:kte)
       u1(kts:kte)     =u(i,kts:kte)
       v1(kts:kte)     =v(i,kts:kte)
       w1(kts:kte)     =w(i,kts:kte)
       th1(kts:kte)    =th(i,kts:kte)
       tk1(kts:kte)    =T3D(i,kts:kte)
       p1(kts:kte)     =p(i,kts:kte)
       ex1(kts:kte)    =exner(i,kts:kte)
       rho1(kts:kte)   =rho(i,kts:kte)
       sqv(kts:kte)    =sqv3D(i,kts:kte) !/(1.+qv(i,kts:kte))
       sqc(kts:kte)    =sqc3D(i,kts:kte) !/(1.+qv(i,kts:kte))
       qv1(kts:kte)    =sqv(kts:kte)/(1.-sqv(kts:kte))
       qc1(kts:kte)    =sqc(kts:kte)/(1.-sqv(kts:kte))
       qi1(kts:kte)    =sqi(kts:kte)/(1.-sqv(kts:kte))
       qs1(kts:kte)    =sqs(kts:kte)/(1.-sqv(kts:kte))
       dqc1(kts:kte)   =0.0
       dqi1(kts:kte)   =0.0
       dqs1(kts:kte)   =0.0
       dqni1(kts:kte)  =0.0
       dqnc1(kts:kte)  =0.0
       dqnwfa1(kts:kte)=0.0
       dqnifa1(kts:kte)=0.0
       dqnbca1(kts:kte)=0.0
       dozone1(kts:kte)=0.0
       IF (FLAG_QNI ) THEN
          qni1(kts:kte)=qni(i,kts:kte)
       ELSE
          qni1(kts:kte)=0.0
       ENDIF
       IF (FLAG_QNC ) THEN
          qnc1(kts:kte)=qnc(i,kts:kte)
       ELSE
          qnc1(kts:kte)=0.0
       ENDIF
       IF (FLAG_QNWFA ) THEN
          qnwfa1(kts:kte)=qnwfa(i,kts:kte)
       ELSE
          qnwfa1(kts:kte)=0.0
       ENDIF
       IF (FLAG_QNIFA ) THEN
          qnifa1(kts:kte)=qnifa(i,kts:kte)
       ELSE
          qnifa1(kts:kte)=0.0
       ENDIF
       IF (FLAG_QNBCA ) THEN
          qnbca1(kts:kte)=qnbca(i,kts:kte)
       ELSE
          qnbca1(kts:kte)=0.0
       ENDIF
       IF (FLAG_OZONE ) THEN
          ozone1(kts:kte)=ozone(i,kts:kte)
       ELSE
          ozone1(kts:kte)=0.0
       ENDIF
       el(kts:kte)  =el_pbl(i,kts:kte)
       qke1(kts:kte)=qke(i,kts:kte)
       sh(kts:kte)  =sh3d(i,kts:kte)
       sm(kts:kte)  =sm3d(i,kts:kte)
       tsq1(kts:kte)=tsq(i,kts:kte)
       qsq1(kts:kte)=qsq(i,kts:kte)
       cov1(kts:kte)=cov(i,kts:kte)
       if (spp_pbl==1) then
          rstoch_col(kts:kte)=pattern_spp_pbl(i,kts:kte)
       else
          rstoch_col(kts:kte)=0.0
       endif
       !edmf
       edmf_a1    =0.0
       edmf_w1    =0.0
       edmf_qc1   =0.0
       s_aw1      =0.0
       s_awthl1   =0.0
       s_awqt1    =0.0
       s_awqv1    =0.0
       s_awqc1    =0.0
       s_awu1     =0.0
       s_awv1     =0.0
       s_awqke1   =0.0
       s_awqnc1   =0.0
       s_awqni1   =0.0
       s_awqnwfa1 =0.0
       s_awqnifa1 =0.0
       s_awqnbca1 =0.0
       ![EWDD]
       edmf_a_dd1 =0.0
       edmf_w_dd1 =0.0
       edmf_qc_dd1=0.0
       sd_aw1     =0.0
       sd_awthl1  =0.0
       sd_awqt1   =0.0
       sd_awqv1   =0.0
       sd_awqc1   =0.0
       sd_awu1    =0.0
       sd_awv1    =0.0
       sd_awqke1  =0.0
       sub_thl    =0.0
       sub_sqv    =0.0
       sub_u      =0.0
       sub_v      =0.0
       det_thl    =0.0
       det_sqv    =0.0
       det_sqc    =0.0
       det_u      =0.0
       det_v      =0.0

       do k = kts,kte
          if (k==kts) then
             zw(k)=0.
          else
             zw(k)=zw(k-1)+dz(i,k-1)
          endif
          !keep snow out for now - increases ceiling bias
          sqw(k)= sqv(k)+sqc(k)+sqi(k)!+sqs(k)
          thl(k)= th1(k) - xlvcp/ex1(k)*sqc(k) &
               &         - xlscp/ex1(k)*(sqi(k))!+sqs(k))
          !Use form from Tripoli and Cotton (1981) with their
          !suggested min temperature to improve accuracy.
          !thl(k)=th(i,k)*(1.- xlvcp/MAX(tk1(k),TKmin)*sqc(k) &
          !    &               - xlscp/MAX(tk1(k),TKmin)*sqi(k))
          thetav(k)=th1(k)*(1.+p608*sqv(k))
       enddo ! end k
       zw(kte+1)=zw(kte)+dz(i,kte)

       !initialize smoke/chem arrays (if used):
       if ( mix_chem ) then
          do ic = 1,ndvel
             vd1(ic) = vdep(i,ic) ! dry deposition velocity
          enddo
          do k = kts,kte
             do ic = 1,nchem
                chem1(k,ic) = chem3d(i,k,ic)
             enddo
          enddo
       else
          do ic = 1,ndvel
             vd1(ic) = 0. ! dry deposition velocity
          enddo
          do k = kts,kte
             do ic = 1,nchem
                chem1(k,ic) = 0.
             enddo
          enddo
       endif
       s_awchem1(kts:kte+1,1:nchem)   = 0.0

!>  - Call get_pblh() to calculate the hybrid \f$\theta_{v}-TKE\f$
!! PBL height diagnostic.
       CALL GET_PBLH_cc(KTS,KTE,PBLH(i),thetav,&
       & Qke1,zw,dz1,xland(i),KPBL(i))

!>  - Call scale_aware() to calculate the similarity functions,
!! \f$P_{\sigma-PBL}\f$ and \f$P_{\sigma-shcu}\f$, to control 
!! the scale-adaptive behaviour for the local and nonlocal 
!! components, respectively.
       if (scaleaware > 0.) then
          call SCALE_AWARE_CC(dx(i),PBLH(i),Psig_bl(i),Psig_shcu(i))
       else
          Psig_bl(i)=1.0
          Psig_shcu(i)=1.0
       endif

       sqcg= 0.0   !ill-defined variable; qcg has been removed
       cpm=cp*(1.+0.84*qv1(kts))
       exnerg=(ps(i)/p1000mb)**rcp

       !-----------------------------------------------------
       !ORIGINAL CODE
       !flt = hfx(i)/( rho(i,kts)*cpm ) &
       ! +xlvcp*ch(i)*(sqc(kts)/exner(i,kts) -sqcg/exnerg)
       !flq = qfx(i)/  rho(i,kts)       &
       !    -ch(i)*(sqc(kts)   -sqcg )
       !-----------------------------------------------------
       flqv = qfx(i)/rho1(kts)
       flqc = 0.0 !currently no sea-spray fluxes, fog settling handled elsewhere
       th_sfc = ts(i)/ex1(kts)

       ! TURBULENT FLUX FOR TKE BOUNDARY CONDITIONS
       flq =flqv+flqc                   !! LATENT
       flt =hfx(i)/(rho1(kts)*cpm )-xlvcp*flqc/ex1(kts)  !! Temperature flux
       fltv=flt + flqv*p608*th_sfc      !! Virtual temperature flux

       ! Update 1/L using updated sfc heat flux and friction velocity
       rmol(i) = -karman*gtr*fltv/max(ust(i)**3,1.0e-6)
       zet = 0.5*dz(i,kts)*rmol(i)
       zet = MAX(zet, -20.)
       zet = MIN(zet,  20.)
       !if(i.eq.idbg)print*,"updated z/L=",zet
       if (bl_mynn_stfunc == 0) then
          !Original Kansas-type stability functions
          if ( zet >= 0.0 ) then
             pmz = 1.0 + (cphm_st-1.0) * zet
             phh = 1.0 +  cphh_st      * zet
          else
             pmz = 1.0/    (1.0-cphm_unst*zet)**0.25 - zet
             phh = 1.0/SQRT(1.0-cphh_unst*zet)
          end if
       else
          !Updated stability functions (Puhales, 2020)
          phi_m = phim(zet)
          pmz   = phi_m - zet
          phh   = phih(zet)
       end if

!>  - Call mym_condensation_cc() to calculate the nonconvective component
!! of the subgrid cloud fraction and mixing ratio as well as the functions
!! used to calculate the buoyancy flux. Different cloud PDFs can be
!! selected by use of the namelist parameter \p bl_mynn_cloudpdf.

       call mym_condensation_cc (kts,kte,                   &
            &dx(i),dz1,zw,xland(i),                      &
            &thl,sqw,sqv,sqc,sqi,sqs,                    &
            &p1,ex1,tsq1,qsq1,cov1,                      &
            &Sh,el,bl_mynn_cloudpdf,                     &
            &qc_bl1D,qi_bl1D,cldfra_bl1D,                &
            &PBLH(i),HFX(i),                             &
            &Vt, Vq, th1, sgm, rmol(i),                  &
            &spp_pbl,rstoch_col,ep_2,ep_3,xlv,r_d,xlvcp, &
            &p608,tv0,r_v,cice,cliq,cpv,cp,xls,rcp)

!>  - Add TKE source driven by cloud top cooling
!!  Calculate the buoyancy production of TKE from cloud-top cooling when
!! \p bl_mynn_topdown =1.
       if (bl_mynn_topdown.eq.1) then
          call topdown_cloudrad_cc(kts,kte,dz1,zw,fltv,     &
               &xland(i),kpbl(i),PBLH(i),                &
               &sqc,sqi,sqw,thl,th1,ex1,p1,rho1,thetav,  &
               &cldfra_bl1D,rthraten(i,:),               &
               &maxKHtopdown(i),KHtopdown,TKEprodTD      )
       else
          maxKHtopdown(i)  = 0.0
          KHtopdown(kts:kte) = 0.0
          TKEprodTD(kts:kte) = 0.0
       endif

       if (bl_mynn_edmf > 0) then
          !PRINT*,"Calling DMP Mass-Flux: i= ",i
          call DMP_mf(                                   &
               &kts,kte,delt,zw,dz1,p1,rho1,             &
               &bl_mynn_edmf_mom,                        &
               &bl_mynn_edmf_tke,                        &
               &bl_mynn_mixscalars,                      &
               &u1,v1,w1,th1,thl,thetav,tk1,             &
               &sqw,sqv,sqc,qke1,                        &
               &qnc1,qni1,qnwfa1,qnifa1,qnbca1,          &
               &ex1,Vt,Vq,sgm,                           &
               &ust(i),flt,fltv,flq,flqv,                &
               &PBLH(i),KPBL(i),DX(i),                   &
               &xland(i),th_sfc,                         &
            ! now outputs - tendencies
            ! &,dth1mf,dqv1mf,dqc1mf,du1mf,dv1mf         &
            ! outputs - updraft properties
               &edmf_a1,edmf_w1,edmf_qt1,                &
               &edmf_thl1,edmf_ent1,edmf_qc1,            &
            ! for the solver
               &s_aw1,s_awthl1,s_awqt1,                  &
               &s_awqv1,s_awqc1,                         &
               &s_awu1,s_awv1,s_awqke1,                  &
               &s_awqnc1,s_awqni1,                       &
               &s_awqnwfa1,s_awqnifa1,s_awqnbca1,        &
               &sub_thl,sub_sqv,                         &
               &sub_u,sub_v,                             &
               &det_thl,det_sqv,det_sqc,                 &
               &det_u,det_v,                             &
            ! chem/smoke mixing
               &nchem,chem1,s_awchem1,                   &
               &mix_chem,                                &
               &qc_bl1D,cldfra_bl1D,                     &
               &qc_bl1D_old,cldfra_bl1D_old,             &
               &FLAG_QC,FLAG_QI,                         &
               &FLAG_QNC,FLAG_QNI,                       &
               &FLAG_QNWFA,FLAG_QNIFA,FLAG_QNBCA,        &
               &Psig_shcu(i),                            &
               &maxwidth(i),ktop_plume(i),               &
               &maxmf(i),ztop_plume(i),                  &
               &spp_pbl,rstoch_col                       )
       endif

       if (bl_mynn_edmf_dd == 1) then
          call DDMF_JPL_cc(kts,kte,delt,zw,dz1,p1,       &
               &u1,v1,th1,thl,thetav,tk1,                &
               &sqw,sqv,sqc,rho1,ex1,                    &
               &ust(i),flt,flq,                          &
               &PBLH(i),KPBL(i),                         &
               &edmf_a_dd1,edmf_w_dd1,edmf_qt_dd1,       &
               &edmf_thl_dd1,edmf_ent_dd1,               &
               &edmf_qc_dd1,                             &
               &sd_aw1,sd_awthl1,sd_awqt1,               &
               &sd_awqv1,sd_awqc1,sd_awu1,sd_awv1,       &
               &sd_awqke1,qc_bl1d,cldfra_bl1d,           &
               &rthraten(i,:),svp1,grav,onethird,        &
               &p1000mb,rcp, xlvcp                       )
       endif

       !Capability to substep the eddy-diffusivity portion
       !do nsub = 1,2
       delt2 = delt !*0.5    !only works if topdown=0

       call mym_turbulence(                              & 
               &kts,kte,xland(i),closure,                &
               &dz1, DX(i), zw,                          &
               &u1, v1, thl, thetav, sqc, sqw,           &
               &qke1, tsq1, qsq1, cov1,                  &
               &vt, vq,                                  &
               &rmol(i), flt, fltv, flq,                 &
               &PBLH(i),th1,                             &
               &Sh,Sm,el,                                &
               &Dfm,Dfh,Dfq,                             &
               &Tcd,Qcd,Pdk,                             &
               &Pdt,Pdq,Pdc,                             &
               &qWT1,qSHEAR1,qBUOY1,qDISS1,              &
               &tke_budget,                              &
               &Psig_bl(i),Psig_shcu(i),                 &
               &cldfra_bl1D,bl_mynn_mixlength,           &
               &edmf_w1,edmf_a1,                         &
               &TKEprodTD,                               &
               &spp_pbl,rstoch_col,debug_code            )

!>  - Call mym_predict() to solve TKE and 
!! \f$\theta^{'2}, q^{'2}, and \theta^{'}q^{'}\f$
!! for the following time step.
       call mym_predict_cc(kts,kte,closure,                 &
               &delt2, dz1,                              &
               &ust(i), flt, flq, pmz, phh,              &
               &el, dfq, rho1, pdk, pdt, pdq, pdc,       &
               &Qke1, Tsq1, Qsq1, Cov1,                  &
               &s_aw1, s_awqke1, bl_mynn_edmf_tke,       &
               &qWT1, qDISS1, tke_budget,                &
               xlvcp,xlscp                               )

       if (dheat_opt > 0) then
          do k=kts,kte-1
             ! Set max dissipative heating rate to 7.2 K per hour
             diss_heat(k) = MIN(MAX(1.0*(qke1(k)**1.5)/(b1*MAX(0.5*(el(k)+el(k+1)),1.))/cp, 0.0),0.002)
             ! Limit heating above 100 mb:
             diss_heat(k) = diss_heat(k) * exp(-10000./MAX(p1(k),1.)) 
          enddo
          diss_heat(kte) = 0.
       else
          diss_heat(1:kte) = 0.
       endif

!>  - Call mynn_tendencies() to solve for tendencies of 
!! \f$U, V, \theta, q_{v}, q_{c}, and q_{i}\f$.
       call mynn_tendencies_cc(kts,kte,delt, dz1, rho1,    &
               &u1,v1,tk1,qv1,ps(i),p1,thl,sqv,sqc,sqw,ust(i),flt,   &
               &flq,flqv,flqc,wspd(i),uoce(i),voce(i),tcd, & 
               &qcd,dfm,dfh,Du1,Dv1,Dth1,diss_heat,        &
               ! mass flux components
               &s_aw1,s_awthl1,s_awqt1,s_awqv1,s_awqc1,    &
               &s_awu1,s_awv1,sub_thl,sub_sqv,sub_u,sub_v, &
               &det_thl,det_sqv,det_sqc,det_u,det_v,       &      
               &FLAG_QC_I,bl_mynn_cloudmix,bl_mynn_mixqt,  &
               &bl_mynn_edmf,bl_mynn_edmf_mom,             &
               &debug_code_i,r_d,p608,ep_2,ep_3,tv0,xlv,   &
               &xlvcp)


       if ( mix_chem ) then
          if (enh_mix) then
                  enh_mix_in=1
          else
                  enh_mix_in=0
          endif
          if ( rrfs_sd ) then 
             call mynn_mix_chem_cc(kts,kte,i,               &
                  &delt, dz1, pblh(i),                   &
                  &nchem, kdvel, ndvel,                  &
                  &chem1, vd1,                           &
                  &rho1,flt,                             &
                  &tcd, qcd,                             &
                  &dfh,                                  &
                  &s_aw1,s_awchem1,                      &
                  &emis_ant_no(i),                       &
                  &frp(i), 1,                            & !//rrfs_sd passed as integer
                  &enh_mix_in                            )
          else
             call mynn_mix_chem_cc(kts,kte,i,               &
                  &delt, dz1, pblh(i),                   &
                  &nchem, kdvel, ndvel,                  &
                  &chem1, vd1,                           &
                  &rho1,flt,                             &
                  &tcd, qcd,                             &
                  &dfh,                                  &
                  &s_aw1,s_awchem1,                      &
                  &zero,                                 &
                  &zero, 0,                              &
                  &enh_mix_in                            )
          endif
          do ic = 1,nchem
             do k = kts,kte
                chem3d(i,k,ic) = max(1.e-12, chem1(k,ic))
             enddo
          enddo
       endif
 
       call retrieve_exchange_coeffs(kts,kte,            &
            &dfm, dfh, dz1, K_m1, K_h1                   )

       !UPDATE 3D ARRAYS
       exch_m(i,kts:kte)  =k_m1(kts:kte)
       exch_h(i,kts:kte)  =k_h1(kts:kte)
       rublten(i,kts:kte) =du1(kts:kte)
       rvblten(i,kts:kte) =dv1(kts:kte)
       rthblten(i,kts:kte)=dth1(kts:kte)
       rqvblten(i,kts:kte)=dqv1(kts:kte)
       if (bl_mynn_cloudmix > 0) then
          if (flag_qc) rqcblten(i,kts:kte)=dqc1(kts:kte)
          if (flag_qi) rqiblten(i,kts:kte)=dqi1(kts:kte)
          if (flag_qs) rqsblten(i,kts:kte)=dqs1(kts:kte)
       else
          if (flag_qc) rqcblten(i,:)=0.
          if (flag_qi) rqiblten(i,:)=0.
          if (flag_qs) rqsblten(i,:)=0.
       endif
       if (bl_mynn_cloudmix > 0 .and. bl_mynn_mixscalars > 0) then
          if (flag_qnc) rqncblten(i,kts:kte)    =dqnc1(kts:kte)
          if (flag_qni) rqniblten(i,kts:kte)    =dqni1(kts:kte)
          if (flag_qnwfa) rqnwfablten(i,kts:kte)=dqnwfa1(kts:kte)
          if (flag_qnifa) rqnifablten(i,kts:kte)=dqnifa1(kts:kte)
          if (flag_qnbca) rqnbcablten(i,kts:kte)=dqnbca1(kts:kte)
       else
          if (flag_qnc) rqncblten(i,:)    =0.
          if (flag_qni) rqniblten(i,:)    =0.
          if (flag_qnwfa) rqnwfablten(i,:)=0.
          if (flag_qnifa) rqnifablten(i,:)=0.
          if (flag_qnbca) rqnbcablten(i,:)=0.
       endif
       dozone(i,kts:kte)=dozone1(kts:kte)
       if (icloud_bl > 0) then
          qc_bl(i,kts:kte)    =qc_bl1D(kts:kte)
          qi_bl(i,kts:kte)    =qi_bl1D(kts:kte)
          cldfra_bl(i,kts:kte)=cldfra_bl1D(kts:kte)
       endif
       el_pbl(i,kts:kte)=el(kts:kte)
       qke(i,kts:kte)   =qke1(kts:kte)
       tsq(i,kts:kte)   =tsq1(kts:kte)
       qsq(i,kts:kte)   =qsq1(kts:kte)
       cov(i,kts:kte)   =cov1(kts:kte)
       sh3d(i,kts:kte)  =sh(kts:kte)
       sm3d(i,kts:kte)  =sm(kts:kte)

       if (tke_budget .eq. 1) then
          !! TKE budget is now given in m**2/s**-3 (Puhales, 2020)
          !! Lower boundary condtions (using similarity relationships such as the prognostic equation for Qke)
          k=kts
          qSHEAR1(k)    =4.*(ust(i)**3*phi_m/(karman*dz(i,k)))-qSHEAR1(k+1) !! staggered
          qBUOY1(k)     =4.*(-ust(i)**3*zet/(karman*dz(i,k)))-qBUOY1(k+1) !! staggered
          !! unstaggering SHEAR and BUOY and trasfering all TKE budget to 3D array               
          do k = kts,kte-1
             qSHEAR(i,k)=0.5*(qSHEAR1(k)+qSHEAR1(k+1)) !!! unstaggering in z
             qBUOY(i,k) =0.5*(qBUOY1(k)+qBUOY1(k+1)) !!! unstaggering in z
             qWT(i,k)   =qWT1(k)
             qDISS(i,k) =qDISS1(k)
             dqke(i,k)  =(qke1(k)-dqke(i,k))*0.5/delt
          enddo
          !! Upper boundary conditions               
          k=kte
          qSHEAR(i,k)   =0.
          qBUOY(i,k)    =0.
          qWT(i,k)      =0.
          qDISS(i,k)    =0.
          dqke(i,k)     =0.
       endif

       !update updraft/downdraft properties
       if (bl_mynn_output > 0) then !research mode == 1
          if (bl_mynn_edmf > 0) then
             edmf_a(i,kts:kte)   =edmf_a1(kts:kte)
             edmf_w(i,kts:kte)   =edmf_w1(kts:kte)
             edmf_qt(i,kts:kte)  =edmf_qt1(kts:kte)
             edmf_thl(i,kts:kte) =edmf_thl1(kts:kte)
             edmf_ent(i,kts:kte) =edmf_ent1(kts:kte)
             edmf_qc(i,kts:kte)  =edmf_qc1(kts:kte)
             sub_thl3D(i,kts:kte)=sub_thl(kts:kte)
             sub_sqv3D(i,kts:kte)=sub_sqv(kts:kte)
             det_thl3D(i,kts:kte)=det_thl(kts:kte)
             det_sqv3D(i,kts:kte)=det_sqv(kts:kte)
          endif
          !if (bl_mynn_edmf_dd > 0) THEN
          !   edmf_a_dd(i,kts:kte)  =edmf_a_dd1(kts:kte)
          !   edmf_w_dd(i,kts:kte)  =edmf_w_dd1(kts:kte)
          !   edmf_qt_dd(i,kts:kte) =edmf_qt_dd1(kts:kte)
          !   edmf_thl_dd(i,kts:kte)=edmf_thl_dd1(kts:kte)
          !   edmf_ent_dd(i,kts:kte)=edmf_ent_dd1(kts:kte)
          !   edmf_qc_dd(i,kts:kte) =edmf_qc_dd1(kts:kte)
          !endif
       endif

       !***  Begin debug prints
       if ( debug_code .and. (i .eq. idbg)) THEN
          if ( ABS(QFX(i))>.001)print*,&
             "SUSPICIOUS VALUES AT: i=",i," QFX=",QFX(i)
          if ( ABS(HFX(i))>1100.)print*,&
             "SUSPICIOUS VALUES AT: i=",i," HFX=",HFX(i)
          do k = kts,kte
             IF ( sh(k) < 0. .OR. sh(k)> 200.)print*,&
                "SUSPICIOUS VALUES AT: i,k=",i,k," sh=",sh(k)
             IF ( ABS(vt(k)) > 2.0 )print*,&
                "SUSPICIOUS VALUES AT: i,k=",i,k," vt=",vt(k)
             IF ( ABS(vq(k)) > 7000.)print*,&
                "SUSPICIOUS VALUES AT: i,k=",i,k," vq=",vq(k)
             IF ( qke(i,k) < -1. .OR. qke(i,k)> 200.)print*,&
                "SUSPICIOUS VALUES AT: i,k=",i,k," qke=",qke(i,k)
             IF ( el_pbl(i,k) < 0. .OR. el_pbl(i,k)> 1500.)print*,&
                "SUSPICIOUS VALUES AT: i,k=",i,k," el_pbl=",el_pbl(i,k)
             IF ( exch_m(i,k) < 0. .OR. exch_m(i,k)> 2000.)print*,&
                "SUSPICIOUS VALUES AT: i,k=",i,k," exxch_m=",exch_m(i,k)
             IF (icloud_bl > 0) then
                IF ( cldfra_bl(i,k) < 0.0 .OR. cldfra_bl(i,k)> 1.)THEN
                   PRINT*,"SUSPICIOUS VALUES: CLDFRA_BL=",cldfra_bl(i,k)," qc_bl=",QC_BL(i,k)
                ENDIF
             ENDIF

             !IF (I==IMD .AND. J==JMD) THEN
             !   PRINT*,"MYNN DRIVER END: k=",k," sh=",sh(k)
             !   PRINT*," sqw=",sqw(k)," thl=",thl(k)," exch_m=",exch_m(i,k)
             !   PRINT*," xland=",xland(i)," rmol=",rmol(i)," ust=",ust(i)
             !   PRINT*," qke=",qke(i,k)," el=",el_pbl(i,k)," tsq=",tsq(i,k)
             !   PRINT*," PBLH=",PBLH(i)," u=",u(i,k)," v=",v(i,k)
             !   PRINT*," vq=",vq(k)," vt=",vt(k)
             !ENDIF
          enddo !end-k
       endif

    enddo !end i-loop

!ACF copy qke into qke_adv if using advection
    IF (bl_mynn_tkeadvect) THEN
       qke_adv=qke
    ENDIF
!ACF-end

#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mynn_bl_driver


 FUNCTION phim(zet)
     ! New stability function parameters for momentum (Puhales, 2020, WRF 4.2.1)
     ! The forms in unstable conditions (z/L < 0) use Grachev et al. (2000), which are a blend of 
     ! the classical (Kansas) forms (i.e., Paulson 1970, Dyer and Hicks 1970), valid for weakly 
     ! unstable conditions (-1 < z/L < 0). The stability functions for stable conditions use an
     ! updated form taken from Cheng and Brutsaert (2005), which extends the validity into very
     ! stable conditions [z/L ~ O(10)].
      IMPLICIT NONE

      real(kind_phys), intent(in):: zet
      real(kind_phys):: dummy_0,dummy_1,dummy_11,dummy_2,dummy_22,dummy_3,dummy_33,dummy_4,dummy_44,dummy_psi
      real(kind_phys), parameter :: am_st=6.1, bm_st=2.5, rbm_st=1./bm_st
      real(kind_phys), parameter :: ah_st=5.3, bh_st=1.1, rbh_st=1./bh_st
      real(kind_phys), parameter :: am_unst=10., ah_unst=34.
      real(kind_phys):: phi_m,phim

      if ( zet >= 0.0 ) then
         dummy_0=1+zet**bm_st
         dummy_1=zet+dummy_0**(rbm_st)
         dummy_11=1+dummy_0**(rbm_st-1)*zet**(bm_st-1)
         dummy_2=(-am_st/dummy_1)*dummy_11
         phi_m = 1-zet*dummy_2
      else
         dummy_0 = (1.0-cphm_unst*zet)**0.25
         phi_m = 1./dummy_0
         dummy_psi = 2.*log(0.5*(1.+dummy_0))+log(0.5*(1.+dummy_0**2))-2.*atan(dummy_0)+1.570796

         dummy_0=(1.-am_unst*zet)          ! parentesis arg
         dummy_1=dummy_0**0.333333         ! y
         dummy_11=-0.33333*am_unst*dummy_0**(-0.6666667) ! dy/dzet
         dummy_2 = 0.33333*(dummy_1**2.+dummy_1+1.)    ! f
         dummy_22 = 0.3333*dummy_11*(2.*dummy_1+1.)    ! df/dzet
         dummy_3 = 0.57735*(2.*dummy_1+1.) ! g
         dummy_33 = 1.1547*dummy_11        ! dg/dzet
         dummy_4 = 1.5*log(dummy_2)-1.73205*atan(dummy_3)+1.813799364 !psic
         dummy_44 = (1.5/dummy_2)*dummy_22-1.73205*dummy_33/(1.+dummy_3**2)! dpsic/dzet

         dummy_0 = zet**2
         dummy_1 = 1./(1.+dummy_0) ! denon
         dummy_11 = 2.*zet         ! denon/dzet
         dummy_2 = ((1-phi_m)/zet+dummy_11*dummy_4+dummy_0*dummy_44)*dummy_1
         dummy_22 = -dummy_11*(dummy_psi+dummy_0*dummy_4)*dummy_1**2

         phi_m = 1.-zet*(dummy_2+dummy_22)
      end if

      !phim = phi_m - zet
      phim = phi_m

  END FUNCTION phim
! ===================================================================

  FUNCTION phih(zet)
    ! New stability function parameters for heat (Puhales, 2020, WRF 4.2.1)
    ! The forms in unstable conditions (z/L < 0) use Grachev et al. (2000), which are a blend of
    ! the classical (Kansas) forms (i.e., Paulson 1970, Dyer and Hicks 1970), valid for weakly
    ! unstable conditions (-1 < z/L < 0). The stability functions for stable conditions use an
    ! updated form taken from Cheng and Brutsaert (2005), which extends the validity into very
    ! stable conditions [z/L ~ O(10)].
      IMPLICIT NONE

      real(kind_phys), intent(in):: zet
      real(kind_phys):: dummy_0,dummy_1,dummy_11,dummy_2,dummy_22,dummy_3,dummy_33,dummy_4,dummy_44,dummy_psi
      real(kind_phys), parameter :: am_st=6.1, bm_st=2.5, rbm_st=1./bm_st
      real(kind_phys), parameter :: ah_st=5.3, bh_st=1.1, rbh_st=1./bh_st
      real(kind_phys), parameter :: am_unst=10., ah_unst=34.
      real(kind_phys):: phh,phih

      if ( zet >= 0.0 ) then
         dummy_0=1+zet**bh_st
         dummy_1=zet+dummy_0**(rbh_st)
         dummy_11=1+dummy_0**(rbh_st-1)*zet**(bh_st-1)
         dummy_2=(-ah_st/dummy_1)*dummy_11
         phih = 1-zet*dummy_2
      else
         dummy_0 = (1.0-cphh_unst*zet)**0.5
         phh = 1./dummy_0
         dummy_psi = 2.*log(0.5*(1.+dummy_0))

         dummy_0=(1.-ah_unst*zet)          ! parentesis arg
         dummy_1=dummy_0**0.333333         ! y
         dummy_11=-0.33333*ah_unst*dummy_0**(-0.6666667) ! dy/dzet
         dummy_2 = 0.33333*(dummy_1**2.+dummy_1+1.)    ! f
         dummy_22 = 0.3333*dummy_11*(2.*dummy_1+1.)    ! df/dzet
         dummy_3 = 0.57735*(2.*dummy_1+1.) ! g
         dummy_33 = 1.1547*dummy_11        ! dg/dzet
         dummy_4 = 1.5*log(dummy_2)-1.73205*atan(dummy_3)+1.813799364 !psic
         dummy_44 = (1.5/dummy_2)*dummy_22-1.73205*dummy_33/(1.+dummy_3**2)! dpsic/dzet

         dummy_0 = zet**2
         dummy_1 = 1./(1.+dummy_0)         ! denon
         dummy_11 = 2.*zet                 ! ddenon/dzet
         dummy_2 = ((1-phh)/zet+dummy_11*dummy_4+dummy_0*dummy_44)*dummy_1
         dummy_22 = -dummy_11*(dummy_psi+dummy_0*dummy_4)*dummy_1**2

         phih = 1.-zet*(dummy_2+dummy_22)
      end if

END FUNCTION phih

! ==================================================================
! ===================================================================
! ===================================================================

END MODULE module_bl_mynn_driver
