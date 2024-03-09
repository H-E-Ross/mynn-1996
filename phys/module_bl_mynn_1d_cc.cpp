#include <algorithm> 
#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <limits>

extern "C" void mynn_tendencies_cc(int kts, int kte, float delt, float* dz, float* rho, float* u, float* v, float* tk, float* qv, float* psfc, float* p, float* thl, float* sqv, float* sqc, float* sqw, float* ust, float flt, float flq, float flqv, float flqc, float wspd, float uoce, float voce, float* tcd, float* qcd, float* dfm, float* dfh, float* Du, float* Dv, float* Dth, float* diss_heat, float* s_aw, float* s_awthl, float* s_awqt, float* s_awqv, float* s_awqc, float* s_awu, float* s_awv, float* sd_aw, float* sd_awthl, float* sd_awqt, float* sd_awqv, float* sd_awqc, float* sd_awu, float* sd_awv, float* sub_thl, float* sub_sqv, float* sub_u, float* sub_v, float* det_thl, float* det_sqv, float* det_sqc, float* det_u, float* det_v, int FLAG_QC, int bl_mynn_cloudmix, int bl_mynn_mixqt, int bl_mynn_edmf_mom, int debug_code, float r_d, float p608, float ep_2,float ep_3,float tv0,float xlv,float xlvcp);

extern "C" void mym_predict_cc(int kts, int kte, float closure, float delt, float* dz, float* ust, float flt, float flq, float pmz, float phh, float* el, float* dfq, float* rho, float* pdk, float* pdt, float* pdq, float* pdc, float* qke, float* tsq, float* qsq, float* cov, float* s_aw, float* s_awqke, int bl_mynn_edmf_tke, int tke_budget, float xlvcp, float xlscp, float karman);

extern "C" void mynn_mix_chem_cc(int kts, int kte, int i,float delt, std::vector<float>& dz, float pblh, int nchem, int kdvel, int ndvel,std::vector<std::vector<float>>& chem1, std::vector<float>& vd1, std::vector<float>& rho,float flt, std::vector<float>& tcd, std::vector<float>& qcd, std::vector<float>& dfh,std::vector<float>& s_aw, std::vector<std::vector<float>>& s_awchem, float emis_ant_no, float frp, int rrfs_sd, int enh_mix); 

extern "C" void moisture_check_cc(int kte, float delt, float* dp, float* exner,float* qv, float* qc, float* qi, float* qs, float* th,float* dqv, float* dqc, float* dqi, float* dqs, float* dth,float dqv2, float xlvcp, float xlscp); 

extern "C" void mym_condensation_cc(int kts, int kte, float dx, float dz[], float zw[], float xland,float thl[], float qw[], float qv[], float qc[], float qi[], float qs[],float p[], float exner[], float tsq[], float qsq[], float cov[], float Sh[], float el[], int bl_mynn_cloudpdf,float qc_bl1D[], float qi_bl1D[], float cldfra_bl1D[], float PBLH1, float HFX1,float Vt[], float Vq[], float th[], float sgm[], float rmo[],int spp_pbl, float rstoch_col[], float ep_2, float ep_3, float xlv, float r_d, float xlvcp, float p608, float tv0, float cpv,float r_v, float cice, float cliq, float cp, float xls, float rcp); 


extern "C" void topdown_cloudrad_cc(int kts, int kte, const std::vector<float>& dz1, const std::vector<float>& zw, float fltv, float xland, int kpbl, float PBLH, const std::vector<float>& sqc, const std::vector<float>& sqi, const std::vector<float>& sqw, const std::vector<float>& thl, const std::vector<float>& th1, const std::vector<float>& ex1, const std::vector<float>& p1, const std::vector<float>& rho1, const std::vector<float>& thetav, const std::vector<float>& cldfra_bl1D, const std::vector<float>& rthraten, float& maxKHtopdown, std::vector<float>& KHtopdown, std::vector<float>& TKEprodTD);

extern "C" void DDMF_JPL_cc(int kts, int kte, float dt, std::vector<float> zw, std::vector<float> dz, std::vector<float> p,std::vector<float> u, std::vector<float> v, std::vector<float> th, std::vector<float> thl, std::vector<float> thv, std::vector<float> tk, std::vector<float> qt, std::vector<float> qv, std::vector<float> qc, std::vector<float> rho, std::vector<float> exner, float ust, float wthl, float wqt, float pblh, int kpbl,std::vector<float>& edmf_a_dd, std::vector<float>& edmf_w_dd, std::vector<float>& edmf_qt_dd,std::vector<float>& edmf_thl_dd, std::vector<float>& edmf_ent_dd, std::vector<float>& edmf_qc_dd,std::vector<float>& sd_aw, std::vector<float>& sd_awthl, std::vector<float>& sd_awqt,std::vector<float>& sd_awqv, std::vector<float>& sd_awqc, std::vector<float>& sd_awu,std::vector<float>& sd_awv, std::vector<float>& sd_awqke,std::vector<float> qc_bl1d, std::vector<float> cldfra_bl1d,std::vector<float> rthraten,float svp1, float grav,float onethird,float p1000mb,float rcp,float xlvcp);

extern "C" void scale_aware_cc(float dx, float pbl1, float& Psig_bl, float& Psig_shcu); 

extern "C" void GET_PBLH_cc(int KTS, int KTE, float& zi, float landsea, const std::vector<float>& thetav1D, const std::vector<float>& qke1D, const std::vector<float>& zw1D, const std::vector<float>& dz1D, int& kzi);

extern "C" void retrieve_exchange_coeffs(int kts, int kte, const std::vector<float>& dfm, const std::vector<float>& dfh, const std::vector<float>& dz, std::vector<float>& K_m, std::vector<float>& K_h);

	
//----------------------------------------CONTSTANTS-------------------------------------------

// Constants
const float NO_threshold = 10.0;     // For anthropogenic sources
const float frp_threshold = 10.0;    // Increased the frp threshold to enhance mixing over big fires
const float pblh_threshold = 100.0;

const float t0c = 273.15; // Assuming t0c is 273.15
const float tice = 240.0; // Assuming tice is 240 based on the comment

// Assuming float corresponds to float precision
const float cphm_st = 5.0, cphm_unst = 16.0,
                 cphh_st = 5.0, cphh_unst = 16.0;

// Closure constants
constexpr float pr = 0.74,
                 g1 = 0.235, // NN2009 = 0.235
                 b1 = 24.0,
                 b2 = 15.0, // CKmod     NN2009
                 c2 = 0.729, // 0.729, //0.75,
                 c3 = 0.340, // 0.340, //0.352,
                 c4 = 0.0,
                 c5 = 0.2,
                 a1 = b1 * (1.0 - 3.0 * g1) / 6.0,
                 c1 = g1 - 1.0 / (3.0 * a1 * std::pow(2.88449914061481660, 1.0/3.0)),
                 a2 = a1 * (g1 - c1) / (g1 * pr),
                 g2 = b2 / b1 * (1.0 - c3) + 2.0 * a1 / b1 * (3.0 - 2.0 * c2);

constexpr float cc2 = 1.0 - c2,
                 cc3 = 1.0 - c3,
                 e1c = 3.0 * a2 * b2 * cc3,
                 e2c = 9.0 * a1 * a2 * cc2,
                 e3c = 9.0 * a2 * a2 * cc2 * (1.0 - c5),
                 e4c = 12.0 * a1 * a2 * cc2,
                 e5c = 6.0 * a1 * a1;

// Constants for min tke in elt integration (qmin), max z/L in els (zmax),
// and factor for eddy viscosity for TKE (Kq = Sqfac*Km):
constexpr float qmin = 0.0, zmax = 1.0, Sqfac = 3.0;

constexpr float gpw = 5.0 / 3.0, qcgmin = 1e-8, qkemin = 1e-12;
constexpr float tliq = 269.0; // all hydrometeors are liquid when T > tliq

// Constants for cloud PDF (mym_condensation)
constexpr float rr2 = 0.7071068, rrp = 0.3989423;

// Use Canuto/Kitamura mod (remove Ric and negative TKE) (1:yes, 0:no)
constexpr float CKmod = 1.0;

// Option to activate environmental subsidence in mass-flux scheme
constexpr bool env_subs = false;

//---------------------------------------------------------------------------------------------
float vsc = 1.0e-5;
float elt = 1.0e-5;

float esat_blend_cc(float t) {

    // Constants for liquid
    const float J0 = .611583699E03;
    const float J1 = .444606896E02;
    const float J2 = .143177157E01;
    const float J3 = .264224321E-1;
    const float J4 = .299291081E-3;
    const float J5 = .203154182E-5;
    const float J6 = .702620698E-8;
    const float J7 = .379534310E-11;
    const float J8 = -.321582393E-13;

    // Constants for ice
    const float K0 = .609868993E03;
    const float K1 = .499320233E02;
    const float K2 = .184672631E01;
    const float K3 = .402737184E-1;
    const float K4 = .565392987E-3;
    const float K5 = .521693933E-5;
    const float K6 = .307839583E-7;
    const float K7 = .105785160E-9;
    const float K8 = .161444444E-12;

    float XC = std::max(-80.0f, t - t0c);
    float esat_blend_cc;

    if (t >= (t0c - 6.0)) {
        esat_blend_cc = J0 + XC * (J1 + XC * (J2 + XC * (J3 + XC * (J4 + XC * (J5 + XC * (J6 + XC * (J7 + XC * J8)))))));
    } else if (t <= tice) {
        esat_blend_cc = K0 + XC * (K1 + XC * (K2 + XC * (K3 + XC * (K4 + XC * (K5 + XC * (K6 + XC * (K7 + XC * K8)))))));
    } else {
        float ESL = J0 + XC * (J1 + XC * (J2 + XC * (J3 + XC * (J4 + XC * (J5 + XC * (J6 + XC * (J7 + XC * J8)))))));
        float ESI = K0 + XC * (K1 + XC * (K2 + XC * (K3 + XC * (K4 + XC * (K5 + XC * (K6 + XC * (K7 + XC * K8)))))));
        float chi = ((t0c - 6.0) - t) / ((t0c - 6.0) - tice);
        esat_blend_cc = (1.0 - chi) * ESL + chi * ESI;
    }

    return esat_blend_cc;
}


float qsat_blend_cc(float t, float P) {
    // Constants for liquid
    const float J0 = .611583699E03;
    const float J1 = .444606896E02;
    const float J2 = .143177157E01;
    const float J3 = .264224321E-1;
    const float J4 = .299291081E-3;
    const float J5 = .203154182E-5;
    const float J6 = .702620698E-8;
    const float J7 = .379534310E-11;
    const float J8 = -.321582393E-13;
    // Constants for ice
    const float K0 = .609868993E03;
    const float K1 = .499320233E02;
    const float K2 = .184672631E01;
    const float K3 = .402737184E-1;
    const float K4 = .565392987E-3;
    const float K5 = .521693933E-5;
    const float K6 = .307839583E-7;
    const float K7 = .105785160E-9;
    const float K8 = .161444444E-12; 
    // Temperature thresholds
    const float t0c = 0.0; // Assuming 0 for t0c (temperature in Celsius)
    const float tice = -273.15; // Assuming -273.15 for tice (absolute zero, could be different)
    float XC = std::max(-80.0f, t - t0c);
    float qsat_blend_cc, ESL, ESI, RSLF, RSIF, chi;

    if (t >= (t0c - 6.0)) {
        ESL = J0 + XC * (J1 + XC * (J2 + XC * (J3 + XC * (J4 + XC * (J5 + XC * (J6 + XC * (J7 + XC * J8)))))));
        ESL = std::min(ESL, P * 0.15f);
        qsat_blend_cc = 0.622 * ESL / std::max(P - ESL, 1e-5f);
    } else if (t <= tice) {
        ESI = K0 + XC * (K1 + XC * (K2 + XC * (K3 + XC * (K4 + XC * (K5 + XC * (K6 + XC * (K7 + XC * K8)))))));
        ESI = std::min(ESI, P * 0.15f);
        qsat_blend_cc = 0.622 * ESI / std::max(P - ESI, 1e-5f);
    } else {
        ESL = J0 + XC * (J1 + XC * (J2 + XC * (J3 + XC * (J4 + XC * (J5 + XC * (J6 + XC * (J7 + XC * J8)))))));
        ESL = std::min(ESL, P * 0.15f);
        ESI = K0 + XC * (K1 + XC * (K2 + XC * (K3 + XC * (K4 + XC * (K5 + XC * (K6 + XC * (K7 + XC * K8)))))));
        ESI = std::min(ESI, P * 0.15f);
        RSLF = 0.622 * ESL / std::max(P - ESL, 1e-5f);
        RSIF = 0.622 * ESI / std::max(P - ESI, 1e-5f);
        chi = ((t0c - 6.0) - t) / ((t0c - 6.0) - tice);
        qsat_blend_cc = (1.0 - chi) * RSLF + chi * RSIF;
    }
    return qsat_blend_cc;
}


float xl_blend_cc(float t,float xlv, float xls, float cpv, float cliq, float cice) {
    float xl_blend_cc, xlvt, xlst, chi;
    // t0c = 273.15, tice is set elsewhere
    if (t >= t0c) {
        xl_blend_cc = xlv + (cpv - cliq) * (t - t0c); // vaporization/condensation
    } else if (t <= tice) {
        xl_blend_cc = xls + (cpv - cice) * (t - t0c); // sublimation/deposition
    } else {
        xlvt = xlv + (cpv - cliq) * (t - t0c); // vaporization/condensation
        xlst = xls + (cpv - cice) * (t - t0c); // sublimation/deposition
        chi = (t0c - t) / (t0c - tice);
        xl_blend_cc = (1. - chi) * xlvt + chi * xlst; // blended
    }
    return xl_blend_cc;
}

void condensation_edmf(float QT, float THL, float P, float zagl, float& THV, float& QC, float p1000mb, float rcp, float xlvcp, float rvovrd) {
    const int niter = 50;
    const float diff = 1.e-6;
    float EXN = std::pow((P / p1000mb), rcp);
    // QC is assumed to be initialized before calling this function
    for (int i = 0; i < niter; ++i) {
        float T = EXN * THL + xlvcp * QC;
        float QS = qsat_blend_cc(T, P);
        float QCOLD = QC;
        QC = 0.5 * QC + 0.5 * std::max((QT - QS), 0.0f);
        if (std::abs(QC - QCOLD) < diff) break;
    }
    float T = EXN * THL + xlvcp * QC;
    float QS = qsat_blend_cc(T, P);
    QC = std::max(QT - QS, 0.0f);
    // Do not allow saturation below 100 m
    if (zagl < 100.0) QC = 0.0;
    THV = (THL + xlvcp * QC) * (1.0 + QT * (rvovrd - 1.0) - rvovrd * QC);
}

// Function to solve system of linear equations on tridiagonal matrix n times n
// after Peaceman and Rachford, 1955
// a, b, c, d - are std::vectors of order n
// a, b, c - are coefficients on the LHS
// d - is initially RHS on the output becomes a solution std::vector
void tridiag_cc(int n, const std::vector<float>& a, const std::vector<float>& b, std::vector<float>& c, std::vector<float>& d) {
    std::vector<float> q(n);
    c[n-1] = 0.0;
    q[0] = -c[0] / b[0];
    d[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        float p = 1.0 / (b[i] + a[i] * q[i - 1]);
        q[i] = -c[i] * p;
        d[i] = (d[i] - a[i] * d[i - 1]) * p;
    }

    for (int i = n - 2; i >= 0; --i) {
        d[i] = d[i] + q[i] * d[i + 1];
    }
}

void tridiag2_cc(int n, float* a, float* b, float* c, float* d, float* x) {
    float* cp = new float[n];
    float* dp = new float[n];
    float m;

    // Initialize c-prime and d-prime
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    // Solve for std::vectors c-prime and d-prime
    for (int i = 1; i < n; ++i) {
        m = b[i] - cp[i - 1] * a[i];
        cp[i] = c[i] / m;
        dp[i] = (d[i] - dp[i - 1] * a[i]) / m;
    }

    // Initialize x
    x[n - 1] = dp[n - 1];

    // Solve for x from the std::vectors c-prime and d-prime
    for (int i = n - 2; i >= 0; --i) {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }

    delete[] cp;
    delete[] dp;
}

// Function to perform tridiagonal matrix algorithm
void tridiag3_cc(int kte, std::vector<float>& a, std::vector<float>& b, std::vector<float>& c, std::vector<float>& d, std::vector<float>& x) {
    // Inversion and resolution of a tridiagonal matrix A X = D
    // a - lower diagonal (Ai,i-1)
    // b - principal diagonal (Ai,i)
    // c - upper diagonal (Ai,i+1)
    // d - right-hand side std::vector
    // x - solution std::vector

    for (int in = kte - 1; in >= 1; --in) {
        d[in] = d[in] - c[in] * d[in + 1] / b[in + 1];
        b[in] = b[in] - c[in] * a[in + 1] / b[in + 1];
    }
    for (int in = 1 + 1; in < kte; ++in) {
        d[in] = d[in] - a[in] * d[in - 1] / b[in - 1];
    }
    for (int in = 1; in < kte; ++in) {
        x[in] = d[in] / b[in];
    }
}


//
// ==================================================================
//     SUBROUTINE  mym_level2:
//
//     Input variables:    see subroutine mym_initialize
//
//     Output variables:
//       dtl(nx,nz,ny) : Vertical gradient of Theta_l             (K/m)
//       dqw(nx,nz,ny) : Vertical gradient of Q_w
//       dtv(nx,nz,ny) : Vertical gradient of Theta_V             (K/m)
//       gm (nx,nz,ny) : G_M divided by L^2/q^2                (s^(-2))
//       gh (nx,nz,ny) : G_H divided by L^2/q^2                (s^(-2))
//       sm (nx,nz,ny) : Stability function for momentum, at Level 2
//       sh (nx,nz,ny) : Stability function for heat, at Level 2
//
//       These are defined on the walls of the grid boxes.
//

//>\ingroup gsd_mynn_edmf
// This subroutine calculates the level 2, non-dimensional wind shear
// \f$G_M\f$ and vertical temperature gradient \f$G_H\f$ as well as
// the level 2 stability funcitons \f$S_h\f$ and \f$S_m\f$.
//\param kts    horizontal dimension
//\param kte    vertical dimension
//\param dz     vertical grid spacings (\f$m\f$)
//\param u      west-east component of the horizontal wind (\f$m s^{-1}\f$)
//\param v      south-north component of the horizontal wind (\f$m s^{-1}\f$)
//\param thl    liquid water potential temperature
//\param qw     total water content \f$Q_w\f$
//\param ql     liquid water content (\f$kg kg^{-1}\f$)
//\param Vt
//\param Vq
//\param dtl     vertical gradient of \f$\theta_l\f$ (\f$K m^{-1}\f$)
//\param dqw     vertical gradient of \f$Q_w\f$
//\param dtv     vertical gradient of \f$\theta_V\f$ (\f$K m^{-1}\f$)
//\param gm      \f$G_M\f$ divided by \f$L^{2}/q^{2}\f$ (\f$s^{-2}\f$)
//\param gh      \f$G_H\f$ divided by \f$L^{2}/q^{2}\f$ (\f$s^{-2}\f$)
//\param sm      stability function for momentum, at Level 2
//\param sh      stability function for heat, at Level 2
//\section gen_mym_level2 GSD MYNN-EDMF mym_level2 General Algorithm
// @ {

void mym_level2(int kts, int kte, float* dz, float* u, float* v,
                float* thl, float* thetav, float* qw, float* ql,
                float* Vt, float* Vq, float* dtl, float* dqw,
                float* dtv, float* gm, float* gh, float* sm, float* sh, 
		float tv0, float gtr, float sqw2) {
    float rfc, f1, f2, rf1, rf2, smc, shc, ri1, ri2, ri3, ri4, duz, dtz, dqz, Vtt, Vqq, dtq, dzk, afk, abk, ri, rf;
    float a2fac;

    rfc = g1 / (g1 + g2);
    f1 = b1 * (g1 - c1) + 3.0 * a2 * (1.0 - c2) * (1.0 - c5) + 2.0 * a1 * (3.0 - 2.0 * c2);
    f2 = b1 * (g1 + g2) - 3.0 * a1 * (1.0 - c2);
    rf1 = b1 * (g1 - c1) / f1;
    rf2 = b1 * g1 / f2;
    smc = a1 / a2 * f1 / f2;
    shc = 3.0 * a2 * (g1 + g2);

    ri1 = 0.5 / smc;
    ri2 = rf1 * smc;
    ri3 = 4.0 * rf2 * smc - 2.0 * ri2;
    ri4 = ri2 * ri2;

    for (int k = kts + 1; k <= kte; ++k) {
        dzk = 0.5 * (dz[k] + dz[k - 1]);
        afk = dz[k] / (dz[k] + dz[k - 1]);
        abk = 1.0 - afk;
        duz = std::pow(u[k] - u[k - 1], 2) + std::pow(v[k] - v[k - 1], 2);
        duz = duz / std::pow(dzk, 2);
        dtz = (thl[k] - thl[k - 1]) / dzk;
        dqz = (qw[k] - qw[k - 1]) / dzk;

        Vtt = 1.0 + Vt[k] * abk + Vt[k - 1] * afk; // Beta-theta in NN09, Eq. 39
        Vqq = tv0 + Vq[k] * abk + Vq[k - 1] * afk; // Beta-q
        dtq = Vtt * dtz + Vqq * dqz;
        // Alternatively, use theta-v without the SGS clouds
        // dtq = (thetav[k] - thetav[k - 1]) / dzk;

        dtl[k] = dtz;
        dqw[k] = dqz;
        dtv[k] = dtq;

        gm[k] = duz;
        gh[k] = -dtq * gtr;

        // Gradient Richardson number
        ri = -gh[k] / std::max(duz, 1.0e-10f);
        // a2fac is needed for the Canuto/Kitamura mod
        if (CKmod == 1) {
            a2fac = 1.0f / (1.0f + std::max(ri, 0.0f));
        } else {
            a2fac = 1.0f;
        }
        rfc = g1 / (g1 + g2);
        f1 = b1 * (g1 - c1) + 3.0 * a2 * a2fac * (1.0 - c2) * (1.0 - c5) + 2.0 * a1 * (3.0 - 2.0 * c2);
        f2 = b1 * (g1 + g2) - 3.0 * a1 * (1.0 - c2);
        rf1 = b1 * (g1 - c1) / f1;
        rf2 = b1 * g1 / f2;
        smc = a1 / (a2 * a2fac) * f1 / f2;
        shc = 3.0 * (a2 * a2fac) * (g1 + g2);
        ri1 = 0.5 / smc;
        ri2 = rf1 * smc;
        ri3 = 4.0 * rf2 * smc - 2.0 * ri2;
        ri4 = ri2 * ri2;

        // Flux Richardson number
        rf = std::min(ri1 * (ri + ri2 - std::sqrt(ri * ri - ri3 * ri + ri4)), rfc);

        sh[k] = shc * (rfc - rf) / (1.0 - rf);
        sm[k] = smc * (rf1 - rf) / (rf2 - rf) * sh[k];
    }
}

// @}

// ==================================================================
//     SUBROUTINE  mym_length:
//
//     Input variables:    see subroutine mym_initialize
//
//     Output variables:   see subroutine mym_initialize
//
//     Work arrays:
//       elt(nx,ny)      : Length scale depending on the PBL depth    (m)
//       vsc(nx,ny)      : Velocity scale q_c                       (m/s)
//                         at first, used for computing elt
//
//     NOTE: the mixing lengths are meant to be calculated at the full-
//           sigmal levels (or interfaces beween the model layers).
//
//>\ingroup gsd_mynn_edmf
// This subroutine calculates the mixing lengths.
void mym_length(int kts, int kte, float xland, float* dz, float* zw, float rmo, float flt, float fltv, float flq, float* Vt, float* Vq, float* u1, float* v1, float* qke, float* dtv, float* el, float zi, float* theta, float Psig_bl, float* cldfra_bl1D, int bl_mynn_mixlength, float* edmf_w1, float* edmf_a1) {
    float cns, alp1, alp2, alp3, alp4, alp5, alp6;
    float minzi = 300.0;
    float maxdz = 750.0;
    float mindz = 300.0;
    float ZSLH = 100.0;
    float CSL = 2.0;
    float afk, abk, zwk, zwk1, dzk, qdz, vflx, bv, tau_cloud, wstar, elb, els, elf, el_stab, el_mf, el_stab_mf, elb_mf, PBLH_PLUS_ENT, Uonset, Ugrid, wt_u, el_les;
    float ctau = 1000.0;
    float karman = 0.4;
    float gtr = 9.81 / theta[0];
    float qmin = 1.0e-10;
    float tv0 = 0.61 * theta[0];
    float onethird = 1.0 / 3.0;
    float grav = 9.81;
    float* qkw = new float[kte-kts+1];
    float qtke[kte+1];
    float thetaw[kte+1];
    float elBLmin[kte+1];
    float elBLavg[kte+1];
    float h1, h2, hs, elBLmin0, elBLavg0, cldavg;
    int i, j, k;

    switch(bl_mynn_mixlength) {
        case 0:
            cns = 2.7;
            alp1 = 0.23;
            alp2 = 1.0;
            alp3 = 5.0;
            alp4 = 100.0;
            alp5 = 0.3;
            zi = std::min(10000.0f, float(zw[kte-2]));
            h1 = std::max(0.3f * float(zi), float(mindz));
            h1 = std::min(float(h1), float(maxdz));
            h2 = h1 / 2.0;
            qkw[kts] = std::sqrt(std::max(float(qke[kts]), 1.0e-10f));
            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0 - afk;
                qkw[k] = std::sqrt(std::max(float(qke[k] * abk + qke[k-1] * afk), 1.0e-3f));
            }
            k = kts + 1;
            zwk = zw[k];
            while (zwk <= zi + h1) {
                dzk = 0.5 * (dz[k] + dz[k-1]);
                qdz = std::max(float(qkw[k] - qmin), 0.03f) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
            }
            elt = alp1 * elt / vsc;
            vflx = (Vt[kts] + 1.0) * flt + (Vq[kts] + tv0) * flq;
            vsc = std::pow(gtr * elt * std::max(float(vflx), 0.0f), onethird);
            el[kts] = 0.0;
            zwk1 = zw[kts+1];
            for (k = kts+1; k <= kte; k++) {
                zwk = zw[k];
                if (dtv[k] > 0.0) {
                    bv = std::sqrt(gtr * dtv[k]);
                    elb = alp2 * qkw[k] / bv * (1.0 + alp3 / alp2 * std::sqrt(vsc / (bv * elt)));
                    elf = alp2 * qkw[k] / bv;
                } else {
                    elb = 1.0e10;
                    elf = elb;
                }
                if (rmo > 0.0) {
                    els = karman * zwk / (1.0 + cns * std::min(zwk * rmo, ZSLH));
                } else {
                    els = karman * zwk * std::pow(1.0 - alp4 * zwk * rmo, 0.2);
                }
                float wt = 0.5 * std::tanh((zwk - (zi + h1)) / h2) + 0.5;
                el[k] = std::min(elb / (elb / elt + elb / els + 1.0f), elf);
            }
            break;
        case 1:
            Ugrid = std::sqrt(u1[kts] * u1[kts] + v1[kts] * v1[kts]);
            Uonset = 15.0;
            wt_u = (1.0 - std::min(std::max(float(Ugrid - Uonset), 0.0f) / 30.0f, 0.5f));
            cns = 2.7;
            alp1 = 0.23;
            alp2 = 0.3;
            alp3 = 2.5 * wt_u;
            alp4 = 5.0;
            alp5 = 0.3;
            alp6 = 50.0;
            zi = std::max(float(zi), float(minzi));
            h1 = std::max(float(0.3f * zi), 300.0f);
            h1 = std::min(float(h1), 600.0f);
            h2 = h1 / 2.0;
            qtke[kts] = std::max(float(0.5f * qke[kts]), 0.01f);
            thetaw[kts] = theta[kts];
            qkw[kts] = std::sqrt(std::max(float(qke[kts]), 1.0e-10f));
            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0 - afk;
                qkw[k] = std::sqrt(std::max(float(qke[k] * abk + qke[k-1] * afk), 1.0e-3f));
                qtke[k] = 0.5 * qkw[k] * qkw[k];
                thetaw[k] = theta[k] * abk + theta[k-1] * afk;
            }
            k = kts + 1;
            zwk = zw[k];
            while (zwk <= zi + h1) {
                dzk = 0.5 * (dz[k] + dz[k-1]);
                qdz = std::min(std::max(float(qkw[k] - qmin), 0.03f), 30.0f) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
            }
            elt = std::min(std::max(float(alp1 * elt / vsc), 10.0f), 400.0f);
            vflx = fltv;
            vsc = std::pow(gtr * elt * std::max(float(vflx), 0.0f), onethird);
            el[kts] = 0.0;
            zwk1 = zw[kts+1];
            for (k = kts+1; k <= kte; k++) {
                zwk = zw[k];
                if (dtv[k] > 0.0) {
                    bv = std::max(float(std::sqrt(gtr * dtv[k])), 0.0001f);
                    elb = std::max(float(alp2 * qkw[k]), float(alp6 * edmf_a1[k-1] * edmf_w1[k-1])) / bv * (1.0 + alp3 * std::sqrt(vsc / (bv * elt)));
                    elb = std::min(elb, zwk);
                    elf = 1.0 * qkw[k] / bv;
                    elBLavg[k] = std::max(float(elBLavg[k]), float(alp6 * edmf_a1[k-1] * edmf_w1[k-1] / bv));
                } else {
                    elb = 1.0e10;
                    elf = elb;
                }
                if (rmo > 0.0) {
                    els = karman * zwk / (1.0 + cns * std::min(zwk * rmo, ZSLH));
                } else {
                    els = karman * zwk * std::pow(1.0 - alp4 * zwk * rmo, 0.2);
                }
                float wt = 0.5 * std::tanh((zwk - (zi + h1)) / h2) + 0.5;
                el[k] = std::min(elb / (elb / elt + elb / els + 1.0f), elf);
                el[k] = el[k] * Psig_bl;
            }
            break;
        case 2:
            Uonset = 3.5 + dz[kts] * 0.1;
            Ugrid = std::sqrt(u1[kts] * u1[kts] + v1[kts] * v1[kts]);
            cns = 3.5;
            alp1 = 0.22;
            alp2 = 0.30;
            alp3 = 2.0;
            alp4 = 5.0;
            alp5 = alp2;
            alp6 = 50.0;
            zi = std::max(float(zi), float(minzi));
            h1 = std::max(float(0.3 * zi), 300.0f);
            h1 = std::min(float(h1), 600.0f);
            h2 = h1 * 0.5;
            qtke[kts] = std::max(float(0.5 * qke[kts]), 0.01f);
            qkw[kts] = std::sqrt(std::max(float(qke[kts]), 1.0e-4f));
            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0 - afk;
                qkw[k] = std::sqrt(std::max(float(qke[k] * abk + qke[k-1] * afk), 1.0e-3f));
                qtke[k] = 0.5 * qkw[k] * qkw[k];
            }
            k = kts + 1;
            zwk = zw[k];
            while (zwk <= zi + h1) {
                dzk = 0.5 * (dz[k] + dz[k-1]);
                qdz = std::min(std::max(float(qkw[k] - qmin), 0.03f), 30.0f) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
            }
            elt = std::min(std::max(float(alp1 * elt / vsc), 10.0f), 400.0f);
            vflx = fltv;
            vsc = std::pow(gtr * elt * std::max(float(vflx), 0.0f), onethird);
            el[kts] = 0.0;
            zwk1 = zw[kts+1];
            for (k = kts+1; k <= kte; k++) {
                zwk = zw[k];
                dzk = 0.5 * (dz[k] + dz[k-1]);
                cldavg = 0.5 * (cldfra_bl1D[k-1] + cldfra_bl1D[k]);
                if (dtv[k] > 0.0) {
                    bv = std::max(float(std::sqrt(gtr * dtv[k])), 0.001f);
                    elb_mf = std::max(float(alp2 * qkw[k]), float(alp6 * edmf_a1[k-1] * edmf_w1[k-1]) / bv * (1.0f + alp3 * std::sqrt(vsc / (bv * elt))));
                    elb = std::min(std::max(float(alp5 * qkw[k]), float(alp6 * edmf_a1[k] * edmf_w1[k]) / bv), float(zwk));
                    wstar = 1.25 * std::pow(gtr * zi * std::max(float(vflx), 1.0e-4f), onethird);
                    tau_cloud = std::min(std::max(float(ctau * wstar / grav), 30.0f), 150.0f);
                    float wt = 0.5 * std::tanh((zwk - (zi + h1)) / h2) + 0.5;
                    tau_cloud = tau_cloud * (1.0 - wt) + 50.0 * wt;
                    elf = std::min(std::max(float(tau_cloud * std::sqrt(std::min(float(qtke[k]), 40.0f))), float(alp6 * edmf_a1[k] * edmf_w1[k] / bv)), float(zwk));
                } else {
                    wstar = 1.25 * std::pow(gtr * zi * std::max(float(vflx), 1.0e-4f), onethird);
                    tau_cloud = std::min(std::max(float(ctau * wstar / grav), 50.0f), 200.0f);
                    float wt = 0.5 * std::tanh((zwk - (zi + h1)) / h2) + 0.5;
                    tau_cloud = tau_cloud * (1.0 - wt) + std::max(100.0f, dzk * 0.25f) * wt;
                    elb = std::min(tau_cloud * std::sqrt(std::min(qtke[k], 40.0f)), zwk);
                    elf = elb;
                    elb_mf = elb;
                }
                elf = elf / (1.0 + (elf / 800.0));
                elb_mf = std::max(float(elb_mf), 0.01f);
                if (rmo > 0.0) {
                    els = karman * zwk / (1.0 + cns * std::min(zwk * rmo, ZSLH));
                } else {
                    els = karman * zwk * std::pow(1.0 - alp4 * zwk * rmo, 0.2);
                }
                float wt = 0.5 * std::tanh((zwk - (zi + h1)) / h2) + 0.5;
                el[k] = std::sqrt(els * els / (1.0 + (els * els / elt * elt) + (els * els / elb_mf * elb_mf)));
                el[k] = el[k] * (1.0 - wt) + elf * wt;
                el[k] = el[k] * Psig_bl + (1.0 - Psig_bl) * el_les;
            }
            break;
    }
}




// called from driver 
void moisture_check_cc(int kte, float delt, float* dp, float* exner,
                    float* qv, float* qc, float* qi, float* qs, float* th,
                    float* dqv, float* dqc, float* dqi, float* dqs, float* dth, 
		    float dqv2, float xlvcp, float xlscp) {

    // Constants (assuming xlvcp and xlscp are defined elsewhere)
    const float qvmin = 1e-20, qcmin = 0.0, qimin = 0.0;

    for (int k = kte; k >= 1; --k) { // From the top to the surface
        float dqc2 = std::max(0.0f, qcmin - qc[k-1]); // Adjusting for 1-based indexing
        float dqi2 = std::max(0.0f, qimin - qi[k-1]);
        float dqs2 = std::max(0.0f, qimin - qs[k-1]);

        // Fix tendencies
        dqc[k-1] += dqc2 / delt;
        dqi[k-1] += dqi2 / delt;
        dqs[k-1] += dqs2 / delt;
        dqv[k-1] -= (dqc2 + dqi2 + dqs2) / delt;
        dth[k-1] += xlvcp / exner[k-1] * (dqc2 / delt) + xlscp / exner[k-1] * ((dqi2 + dqs2) / delt);

        // Update species
        qc[k-1] += dqc2;
        qi[k-1] += dqi2;
        qs[k-1] += dqs2;
        qv[k-1] -= dqc2 + dqi2 + dqs2;
        th[k-1] += xlvcp / exner[k-1] * dqc2 + xlscp / exner[k-1] * (dqi2 + dqs2);

        // Then fix qv
        float dqv2 = std::max(0.0f, qvmin - qv[k-1]);
        dqv[k-1] += dqv2 / delt;
        qv[k-1] += dqv2;
        if (k != 1) {
            qv[k-2] -= dqv2 * dp[k-1] / dp[k-2]; // Adjusting for 1-based indexing
            dqv[k-2] -= dqv2 * dp[k-1] / dp[k-2] / delt;
        }
        qv[k-1] = std::max(float(qv[k-1]), float(qvmin));
        qc[k-1] = std::max(float(qc[k-1]), float(qcmin));
        qi[k-1] = std::max(float(qi[k-1]), float(qimin));
        qs[k-1] = std::max(float(qs[k-1]), float(qimin));
    }

        float sum = 0.0;
    float aa, dum;

    // Only execute if dqv2 > 1.e-20, which indicates adjustment was made at the top layer
    if(dqv2 > 1e-20) {
        for (int k = 1; k <= kte; ++k) { // Loop through all layers
            if (qv[k-1] > 2.0 * qvmin) {
                sum += qv[k-1] * dp[k-1];
            }
        }

        aa = dqv2 * dp[0] / std::max(1.e-20f, sum); // Adjust for 1-based indexing with dp[0]

        if (aa < 0.5) {
            for (int k = 1; k <= kte; ++k) { // Loop through all layers again
                if (qv[k-1] > 2.0 * qvmin) {
                    dum = aa * qv[k-1];
                    qv[k-1] -= dum;
                    dqv[k-1] -= dum / delt;
                }
            }
        } else {
            // For testing purposes only (not yet found in any output):
            // std::cout << "Full moisture conservation is impossible" << std::endl;
        }
    }

}

/*
! ==================================================================
!     SUBROUTINE  mym_predict:
!
!     Input variables:    see subroutine mym_initialize and turbulence
!       qke(nx,nz,ny) : qke at (n)th time level
!       tsq, ...cov     : ditto
!
!     Output variables:
!       qke(nx,nz,ny) : qke at (n+1)th time level
!       tsq, ...cov     : ditto
!
!     Work arrays:
!       qkw(nx,nz,ny)   : q at the center of the grid boxes        (m/s)
!       bp (nx,nz,ny)   : = 1/2*F,     see below
!       rp (nx,nz,ny)   : = P-1/2*F*Q, see below
!
!     # The equation for a turbulent quantity Q can be expressed as
!          dQ/dt + Ah + Av = Dh + Dv + P - F*Q,                      (1)
!       where A is the advection, D the diffusion, P the production,
!       F*Q the dissipation and h and v denote horizontal and vertical,
!       respectively. If Q is q^2, F is 2q/B_1L.
!       Using the Crank-Nicholson scheme for Av, Dv and F*Q, a finite
!       difference equation is written as
!          Q{n+1} - Q{n} = dt  *( Dh{n}   - Ah{n}   + P{n} )
!                        + dt/2*( Dv{n}   - Av{n}   - F*Q{n}   )
!                        + dt/2*( Dv{n+1} - Av{n+1} - F*Q{n+1} ),    (2)
!       where n denotes the time level.
!       When the advection and diffusion terms are discretized as
!          dt/2*( Dv - Av ) = a(k)Q(k+1) - b(k)Q(k) + c(k)Q(k-1),    (3)
!       Eq.(2) can be rewritten as
!          - a(k)Q(k+1) + [ 1 + b(k) + dt/2*F ]Q(k) - c(k)Q(k-1)
!                 = Q{n} + dt  *( Dh{n}   - Ah{n}   + P{n} )
!                        + dt/2*( Dv{n}   - Av{n}   - F*Q{n}   ),    (4)
!       where Q on the left-hand side is at (n+1)th time level.
!
!       In this subroutine, a(k), b(k) and c(k) are obtained from
!       subprogram coefvu and are passed to subprogram tinteg via
!       common. 1/2*F and P-1/2*F*Q are stored in bp and rp,
!       respectively. Subprogram tinteg solves Eq.(4).
!
!       Modify this subroutine according to your numerical integration
!       scheme (program).
!
!-------------------------------------------------------------------
!>\ingroup gsd_mynn_edmf
!! This subroutine predicts the turbulent quantities at the next step.
*/
void mym_predict_cc(int kts, int kte, float closure, float delt, float* dz, float* ust, float flt, float flq, float pmz, float phh, float* el, float* dfq, float* rho, float* pdk, float* pdt, float* pdq, float* pdc, float* qke, float* tsq, float* qsq, float* cov, float* s_aw, float* s_awqke, int bl_mynn_edmf_tke, int tke_budget, float xlvcp, float xlscp, float karman) {
    float vkz, pdk1, phm, pdt1, pdq1, pdc1, b1l, b2l, onoff;
    float* dtz = new float[kte-kts+1];
    float* a = new float[kte-kts+1];
    float* b = new float[kte-kts+1];
    float* c = new float[kte-kts+1];
    float* d = new float[kte-kts+1];
    float* x = new float[kte-kts+1];
    float* rhoinv = new float[kte-kts+1];
    float* rhoz = new float[kte-kts+2];
    float* kqdz = new float[kte-kts+2];
    float* kmdz = new float[kte-kts+2];
    float* qkw = new float[kte-kts+1];
    float* bp = new float[kte-kts+1];
    float* rp = new float[kte-kts+1];
    float* df3q = new float[kte-kts+1];
    float* tke_up = new float[kte-kts+1];
    float* dzinv = new float[kte-kts+1];
    
    // REGULATE THE MOMENTUM MIXING FROM THE MASS-FLUX SCHEME (on or off)
    if (bl_mynn_edmf_tke == 0) {
        onoff = 0.0;
    } else {
        onoff = 1.0;
    }
    
    // Calculate vkz
    vkz = karman * 0.5 * dz[kts];
    
    // Calculate df3q and dtz
    for (int k = kts; k <= kte; k++) {
        qkw[k] = sqrt(std::max(qke[k], 0.0f));
        df3q[k] = Sqfac * dfq[k];
        dtz[k] = delt / dz[k];
    }
    
    // Prepare "constants" for diffusion equation
    rhoz[kts] = rho[kts];
    rhoinv[kts] = 1.0 / rho[kts];
    kqdz[kts] = rhoz[kts] * df3q[kts];
    kmdz[kts] = rhoz[kts] * dfq[kts];
    for (int k = kts+1; k <= kte; k++) {
        rhoz[k] = (rho[k] * dz[k-1] + rho[k-1] * dz[k]) / (dz[k-1] + dz[k]);
        rhoz[k] = std::max(rhoz[k], 1E-4f);
        rhoinv[k] = 1.0 / std::max(rho[k], 1E-4f);
        kqdz[k] = rhoz[k] * df3q[k];
        kmdz[k] = rhoz[k] * dfq[k];
    }
    rhoz[kte+1] = rhoz[kte];
    kqdz[kte+1] = rhoz[kte+1] * df3q[kte];
    kmdz[kte+1] = rhoz[kte+1] * dfq[kte];
    
    // Calculate pdk1, phm, pdt1, pdq1, pdc1
    pdk1 = 2.0f * pow(*ust, 3.0f) * pmz / vkz;
    phm = 2.0f / *ust * phh / vkz;
    pdt1 = phm * pow(flt, 2);
    pdq1 = phm * pow(flq, 2);
    pdc1 = phm * flt * flq;
    
    // Calculate pdk, pdt, pdq, pdc
    pdk[kts] = pdk1 - pdk[kts+1];
    pdt[kts] = pdt[kts+1];
    pdq[kts] = pdq[kts+1];
    pdc[kts] = pdc[kts+1];
    
    // Prediction of twice the turbulent kinetic energy
    for (int k = kts; k <= kte-1; k++) {
        b1l = b1 * 0.5 * (el[k+1] + el[k]);
        bp[k] = 2.0 * qkw[k] / b1l;
        rp[k] = pdk[k+1] + pdk[k];
    }
    for (int k = kts; k <= kte-1; k++) {
        a[k] = -dtz[k] * df3q[k] + 0.5 * dtz[k] * s_aw[k] * onoff;
        b[k] = 1.0 + dtz[k] * (df3q[k] + df3q[k+1]) + 0.5 * dtz[k] * (s_aw[k] - s_aw[k+1]) * onoff + bp[k] * delt;
        c[k] = -dtz[k] * df3q[k+1] - 0.5 * dtz[k] * s_aw[k+1] * onoff;
        d[k] = rp[k] * delt + qke[k] + dtz[k] * (s_awqke[k] - s_awqke[k+1]) * onoff;
    }
    a[kte] = 0.0;
    b[kte] = 1.0;
    c[kte] = 0.0;
    d[kte] = qke[kte];
    tridiag2_cc(kte, a, b, c, d, x);
    for (int k = kts; k <= kte; k++) {
        qke[k] = std::max(x[k], 1e-4f);
        qke[k] = std::min(qke[k], 150.0f);
    }
    
    // TKE budget
    if (tke_budget == 1) {
        float* qWT1D = new float[kte-kts+1];
        float* qDISS1D = new float[kte-kts+1];
        float* tke_up = new float[kte-kts+1];
        float* dzinv = new float[kte-kts+1];
        
        // TKE Vertical transport
	for (int k=kts; k <=kte; k++) 
	{
		tke_up[k] = 0.5f * qke[k];
                dzinv[k] = 1.0f / dz[k];
	}
        qWT1D[kts] = dzinv[kts] * ((kqdz[kts+1] * (tke_up[kts+1] - tke_up[kts])) - (kqdz[kts] * tke_up[kts])) + 0.5 * rhoinv[kts] * (s_aw[kts+1] * tke_up[kts+1] + ((s_aw[kts+1] - s_aw[kts]) * tke_up[kts]) + (s_awqke[kts] - s_awqke[kts+1])) * onoff;
        for (int k = kts+1; k <= kte-1; k++) {
            qWT1D[k] = dzinv[k] * ((kqdz[k+1] * (tke_up[k+1] - tke_up[k])) - (kqdz[k] * (tke_up[k] - tke_up[k-1]))) + 0.5 * rhoinv[k] * (s_aw[k+1] * tke_up[k+1] + ((s_aw[k+1] - s_aw[k]) * tke_up[k]) - (s_aw[k] * tke_up[k-1]) + (s_awqke[k] - s_awqke[k+1])) * onoff;
        }
        qWT1D[kte] = dzinv[kte] * (-(kqdz[kte] * (tke_up[kte] - tke_up[kte-1]))) + 0.5 * rhoinv[kte] * (-(s_aw[kte] * tke_up[kte]) - (s_aw[kte] * tke_up[kte-1]) + s_awqke[kte]) * onoff;
        
        // TKE dissipation rate
	for (int k=kts; k <=kte; k++) 
	{
		qDISS1D[k] = bp[k] * tke_up[k];
	}
    }
    
    if (closure > 2.5) {
        // Prediction of the moisture variance
        for (int k = kts; k <= kte-1; k++) {
            b2l = b2 * 0.5 * (el[k+1] + el[k]);
            bp[k] = 2.0 * qkw[k] / b2l;
            rp[k] = pdq[k+1] + pdq[k];
        }
        for (int k = kts; k <= kte-1; k++) {
            a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
            b[k] = 1.0 + dtz[k] * (kmdz[k] + kmdz[k+1]) * rhoinv[k] + bp[k] * delt;
            c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k];
            d[k] = rp[k] * delt + qsq[k];
        }
        a[kte] = -1.0;
        b[kte] = 1.0;
        c[kte] = 0.0;
        d[kte] = 0.0;
        tridiag2_cc(kte, a, b, c, d, x);
        for (int k = kts; k <= kte; k++) {
            qsq[k] = std::max(x[k], 1e-17f);
        }
    } else {
        // Level 2.5 - use level 2 diagnostic
        for (int k = kts; k <= kte-1; k++) {
            if (qkw[k] <= 0.0) {
                b2l = 0.0;
            } else {
                b2l = b2 * 0.25 * (el[k+1] + el[k]) / qkw[k];
            }
            qsq[k] = b2l * (pdq[k+1] + pdq[k]);
        }
        qsq[kte] = qsq[kte-1];
    }
    
    if (closure >= 3.0) {
        // Prediction of the temperature variance
        for (int k = kts; k <= kte-1; k++) {
            b2l = b2 * 0.5 * (el[k+1] + el[k]);
            bp[k] = 2.0 * qkw[k] / b2l;
            rp[k] = pdt[k+1] + pdt[k];
        }
        for (int k = kts; k <= kte-1; k++) {
            a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
            b[k] = 1.0 + dtz[k] * (kmdz[k] + kmdz[k+1]) * rhoinv[k] + bp[k] * delt;
            c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k];
            d[k] = rp[k] * delt + tsq[k];
        }
        a[kte] = 0.0;
        b[kte] = 1.0;
        c[kte] = 0.0;
        d[kte] = tsq[kte];
        tridiag2_cc(kte, a, b, c, d, x);
        for (int k = kts; k <= kte; k++) {
            tsq[k] = x[k];
        }
        
        // Prediction of the temperature-moisture covariance
        for (int k = kts; k <= kte-1; k++) {
            b2l = b2 * 0.5 * (el[k+1] + el[k]);
            bp[k] = 2.0 * qkw[k] / b2l;
            rp[k] = pdc[k+1] + pdc[k];
        }
        for (int k = kts; k <= kte-1; k++) {
            a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
            b[k] = 1.0 + dtz[k] * (kmdz[k] + kmdz[k+1]) * rhoinv[k] + bp[k] * delt;
            c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k];
            d[k] = rp[k] * delt + cov[k];
        }
        a[kte] = 0.0;
        b[kte] = 1.0;
        c[kte] = 0.0;
        d[kte] = 0.0;
        tridiag2_cc(kte, a, b, c, d, x);
        for (int k = kts; k <= kte; k++) {
            cov[k] = x[k];
        }
    } else {
        // Not level 3 - default to level 2 diagnostic
        for (int k = kts; k <= kte-1; k++) {
            if (qkw[k] <= 0.0) {
                b2l = 0.0;
            } else {
                b2l = b2 * 0.25 * (el[k+1] + el[k]) / qkw[k];
            }
            tsq[k] = b2l * (pdt[k+1] + pdt[k]);
            cov[k] = b2l * (pdc[k+1] + pdc[k]);
        }
        tsq[kte] = tsq[kte-1];
        cov[kte] = cov[kte-1];
    }
    
    delete[] dtz;
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    delete[] x;
    delete[] rhoinv;
    delete[] rhoz;
    delete[] kqdz;
    delete[] kmdz;
    delete[] qkw;
    delete[] bp;
    delete[] rp;
    delete[] df3q;
    delete[] tke_up;
    delete[] dzinv;
}

void mynn_mix_chem_cc(int kts, int kte, int i,
                   float delt, std::vector<float>& dz, float pblh,
                   int nchem, int kdvel, int ndvel,
                   std::vector<std::vector<float>>& chem1, std::vector<float>& vd1,
                   std::vector<float>& rho,
                   float flt, std::vector<float>& tcd, std::vector<float>& qcd,
                   std::vector<float>& dfh,
                   std::vector<float>& s_aw, std::vector<std::vector<float>>& s_awchem,
                   float emis_ant_no, float frp, int rrfs_sd, int enh_mix) {

    // Local vars
    std::vector<float> dtz(kte - kts + 1);
    std::vector<float> a(kte - kts + 1), b(kte - kts + 1), c(kte - kts + 1), d(kte - kts + 1), x(kte - kts + 1);
    float dztop = 0.5 * (dz[kte - 1] + dz[kte - 2]);
    for (int k = kts; k <= kte; ++k) {
        dtz[k - kts] = delt / dz[k - 1];
    }
    // Prepare "constants" for diffusion equation.
    std::vector<float> rhoz(kte - kts + 2), khdz(kte - kts + 2), rhoinv(kte - kts + 1);
    rhoz[0] = rho[kts - 1];
    rhoinv[0] = 1.0 / rho[kts - 1];
    khdz[0] = rhoz[0] * dfh[kts - 1];
    for (int k = kts + 1; k <= kte; ++k) {
        rhoz[k - kts] = (rho[k - 1] * dz[k - 2] + rho[k - 2] * dz[k - 1]) / (dz[k - 2] + dz[k - 1]);
        rhoz[k - kts] = std::max(float(rhoz[k - kts]), 1e-4f);
        rhoinv[k - kts] = 1.0 / std::max(float(rho[k - 1]), 1e-4f);
        float dzk = 0.5 * (dz[k - 1] + dz[k - 2]);
        khdz[k - kts] = rhoz[k - kts] * dfh[k - 1];
    }
    rhoz[kte - kts + 1] = rhoz[kte - kts];
    khdz[kte - kts + 1] = rhoz[kte - kts + 1] * dfh[kte - 1];
    // Stability criteria for mf
    for (int k = kts + 1; k <= kte - 1; ++k) {
        khdz[k - kts] = std::max(float(khdz[k - kts]), float(0.5 * s_aw[k - kts]));
        khdz[k - kts] = std::max(float(khdz[k - kts]), float(-0.5 * (s_aw[k - kts] - s_aw[k - kts + 1])));
    }
    // Enhanced mixing over fires
    if (rrfs_sd==1 && enh_mix==1) {
        for (int k = kts + 1; k <= kte - 1; ++k) {
            float khdz_old = khdz[k - kts];
            float khdz_back = pblh * 0.15 / dz[k - 1];
            // Modify based on anthropogenic emissions of NO and FRP
            if (pblh < pblh_threshold) {
                if (emis_ant_no > NO_threshold) {
                    khdz[k - kts] = std::max(1.1f * float(khdz[k - kts]), float(std::sqrt((emis_ant_no / NO_threshold)) / dz[k - 1] * rhoz[k - kts]));
                }
                if (frp > frp_threshold) {
                    int kmaxfire = std::ceil(std::log(frp));
                    khdz[k - kts] = std::max(float(1.1 * khdz[k - kts]), float((1.0 - k / (kmaxfire * 2.0)) * (std::pow(std::log(frp), 2.0) - 2.0 * std::log(frp)) / dz[k - 1] * rhoz[k - kts]));
                }
            }
        }
    }
    // Mixing of chemical species
    for (int ic = 0; ic < nchem; ++ic) {
        int k = kts;
        a[0] = -dtz[0] * khdz[0] * rhoinv[0];
        b[0] = 1.0 + dtz[0] * (khdz[1] + khdz[0]) * rhoinv[0] - 0.5 * dtz[0] * rhoinv[0] * s_aw[1];
        c[0] = -dtz[0] * khdz[1] * rhoinv[0] - 0.5 * dtz[0] * rhoinv[0] * s_aw[1];
        d[0] = chem1[k - 1][ic] - dtz[0] * vd1[ic] * chem1[k - 1][ic] - dtz[0] * rhoinv[0] * s_awchem[1][ic];
        for (k = kts + 1; k <= kte - 1; ++k) {
            a[k - kts] = -dtz[k - kts] * khdz[k - kts] * rhoinv[k - kts] + 0.5 * dtz[k - kts] * rhoinv[k - kts] * s_aw[k - kts];
            b[k - kts] = 1.0 + dtz[k - kts] * (khdz[k - kts] + khdz[k - kts + 1]) * rhoinv[k - kts] + 0.5 * dtz[k - kts] * rhoinv[k - kts] * (s_aw[k - kts] - s_aw[k - kts + 1]);
            c[k - kts] = -dtz[k - kts] * khdz[k - kts + 1] * rhoinv[k - kts] - 0.5 * dtz[k - kts] * rhoinv[k - kts] * s_aw[k - kts + 1];
            d[k - kts] = chem1[k - 1][ic] + dtz[k - kts] * rhoinv[k - kts] * (s_awchem[k - kts][ic] - s_awchem[k - kts + 1][ic]);
        }
        // Prescribed value at top
        a[kte - kts] = 0.0;
        b[kte - kts] = 1.0;
        c[kte - kts] = 0.0;
        d[kte - kts] = chem1[kte - 1][ic];
        tridiag3_cc(kte, a, b, c, d, x);
        for (k = kts; k <= kte; ++k) {
            chem1[k - 1][ic] = x[k - kts];
        }
    }
}


// ==================================================================
//>\ingroup gsd_mynn_edmf
// This subroutine solves for tendencies of U, V, \f$\theta\f$, qv,
// qc, and qi
void mynn_tendencies_cc(int kts, int kte, float delt, float* dz, float* rho, float* u, float* v, float* tk, float* qv, float* psfc, float* p, float* thl, float* sqv, float* sqc, float* sqw, float* ust, float flt, float flq, float flqv, float flqc, float wspd, float uoce, float voce, float* tcd, float* qcd, float* dfm, float* dfh, float* Du, float* Dv, float* Dth, float* diss_heat, float* s_aw, float* s_awthl, float* s_awqt, float* s_awqv, float* s_awqc, float* s_awu, float* s_awv, float* sd_aw, float* sd_awthl, float* sd_awqt, float* sd_awqv, float* sd_awqc, float* sd_awu, float* sd_awv, float* sub_thl, float* sub_sqv, float* sub_u, float* sub_v, float* det_thl, float* det_sqv, float* det_sqc, float* det_u, float* det_v, int FLAG_QC, int bl_mynn_cloudmix, int bl_mynn_mixqt, int bl_mynn_edmf_mom, int debug_code, float r_d, float p608, float ep_2,float ep_3,float tv0,float xlv,float xlvcp) {
    float nonloc = 1.0;
    float dztop = 0.5 * (dz[kte] + dz[kte-1]);
    float onoff = (bl_mynn_edmf_mom == 0) ? 0.0 : 1.0;
    float rhosfc = *psfc / (r_d * (tk[kts] + p608 * qv[kts]));
    float* dtz = new float[kte+1];
    float* dfhc = new float[kte]; 
    float* dfmc = new float[kte];
    float* delp = new float[kte]; 
    float* sqv2 = new float[kte]; 
    float* sqc2 = new float[kte];
    float* rhoz = new float[kte+1];
    float* khdz = new float[kte+1];
    float* kmdz = new float[kte+1];
    float* rhoinv = new float[kte+1];
    float* sqw2 = new float[kte+1];
    float* a = new float[kte+1];
    float* b = new float[kte+1]; 
    float* c = new float[kte+1];
    float* d = new float[kte+1];
    float* x = new float[kte+1];
    float qvflux;
    float ust_v=*ust;
    
    dtz[kts] = delt / dz[kts];
    rhoz[kts] = rho[kts];
    rhoinv[kts] = 1.0 / rho[kts];
    khdz[kts] = rhoz[kts] * dfh[kts];
    kmdz[kts] = rhoz[kts] * dfm[kts];
    delp[kts] = *psfc - (p[kts+1] * dz[kts] + p[kts] * dz[kts+1]) / (dz[kts] + dz[kts+1]);
    
    for (int k = kts+1; k <= kte; k++) {
        dtz[k] = delt / dz[k];
        rhoz[k] = (rho[k] * dz[k-1] + rho[k-1] * dz[k]) / (dz[k-1] + dz[k]);
        rhoz[k] = std::max(float(rhoz[k]), 1e-4f);
        rhoinv[k] = 1.0 / std::max(float(rho[k]), 1e-4f);
        float dzk = 0.5 * (dz[k] + dz[k-1]);
        khdz[k] = rhoz[k] * dfh[k];
        kmdz[k] = rhoz[k] * dfm[k];
    }
    
    for (int k = kts+1; k <= kte-1; k++) {
        delp[k] = (p[k] * dz[k-1] + p[k-1] * dz[k]) / (dz[k] + dz[k-1]) - (p[k+1] * dz[k] + p[k] * dz[k+1]) / (dz[k] + dz[k+1]);
    }
    delp[kte] = delp[kte-1];
    rhoz[kte+1] = rhoz[kte];
    khdz[kte+1] = rhoz[kte+1] * dfh[kte];
    kmdz[kte+1] = rhoz[kte+1] * dfm[kte];
    
    for (int k = kts+1; k <= kte-1; k++) {
        khdz[k] = std::max(float(khdz[k]), float(0.5 * s_aw[k]));
        khdz[k] = std::max(float(khdz[k]), float(-0.5 * (s_aw[k] - s_aw[k+1])));
        kmdz[k] = std::max(float(kmdz[k]), float(0.5 * s_aw[k]));
        kmdz[k] = std::max(float(kmdz[k]), float(-0.5 * (s_aw[k] - s_aw[k+1])));
    }
    
    float ustdrag = std::min(ust_v * ust_v, 0.99f) / wspd;
    float ustdiff = std::min(ust_v * ust_v, 0.01f) / wspd;
    
    float dth[kte+1];
    for (int k = kts; k <= kte; k++) {
        dth[k] = 0.0;
    }
    
    int k = kts;
    a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
    b[k] = 1.0 + dtz[k] * (kmdz[k+1] + rhosfc * ust_v * ust_v / wspd) * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
    c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
    d[k] = u[k] + dtz[k] * uoce * ust_v * ust_v / wspd - dtz[k] * rhoinv[k] * s_awu[k+1] * onoff + dtz[k] * rhoinv[k] * sd_awu[k+1] * onoff + sub_u[k] * delt + det_u[k] * delt;
    
    for (int k = kts+1; k <= kte-1; k++) {
        a[k] = -dtz[k] * kmdz[k] * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * s_aw[k] * onoff + 0.5 * dtz[k] * rhoinv[k] * sd_aw[k] * onoff;
        b[k] = 1.0 + dtz[k] * (kmdz[k] + kmdz[k+1]) * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) * onoff + 0.5 * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]) * onoff;
        c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
        d[k] = u[k] + dtz[k] * rhoinv[k] * (s_awu[k] - s_awu[k+1]) * onoff - dtz[k] * rhoinv[k] * (sd_awu[k] - sd_awu[k+1]) * onoff + sub_u[k] * delt + det_u[k] * delt;
    }
    
    a[kte] = 0.0;
    b[kte] = 1.0;
    c[kte] = 0.0;
    d[kte] = u[kte];
    tridiag2_cc(kte, a, b, c, d, x);
    for (int k = kts; k <= kte; k++) {
        Du[k] = (x[k] - u[k]) / delt;
    }
    
    k = kts;
    a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
    b[k] = 1.0 + dtz[k] * (kmdz[k+1] + rhosfc * ust_v * ust_v / wspd) * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
    c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
    d[k] = v[k] + dtz[k] * voce * ust_v * ust_v / wspd - dtz[k] * rhoinv[k] * s_awv[k+1] * onoff + dtz[k] * rhoinv[k] * sd_awv[k+1] * onoff + sub_v[k] * delt + det_v[k] * delt;
    
    for (int k = kts+1; k <= kte-1; k++) {
        a[k] = -dtz[k] * kmdz[k] * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * s_aw[k] * onoff + 0.5 * dtz[k] * rhoinv[k] * sd_aw[k] * onoff;
        b[k] = 1.0 + dtz[k] * (kmdz[k] + kmdz[k+1]) * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) * onoff + 0.5 * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]) * onoff;
        c[k] = -dtz[k] * kmdz[k+1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] * onoff - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1] * onoff;
        d[k] = v[k] + dtz[k] * rhoinv[k] * (s_awv[k] - s_awv[k+1]) * onoff - dtz[k] * rhoinv[k] * (sd_awv[k] - sd_awv[k+1]) * onoff + sub_v[k] * delt + det_v[k] * delt;
    }
    
    a[kte] = 0.0;
    b[kte] = 1.0;
    c[kte] = 0.0;
    d[kte] = v[kte];
    tridiag2_cc(kte, a, b, c, d, x);
    for (int k = kts; k <= kte; k++) {
        Dv[k] = (x[k] - v[k]) / delt;
    }
    
    k = kts;
    a[k] = -dtz[k] * khdz[k] * rhoinv[k];
    b[k] = 1.0 + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1];
    c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1];
    d[k] = thl[k] + dtz[k] * rhosfc * flt * rhoinv[k] + tcd[k] * delt - dtz[k] * rhoinv[k] * s_awthl[k+1] - dtz[k] * rhoinv[k] * sd_awthl[k+1] + diss_heat[k] * delt + sub_thl[k] * delt + det_thl[k] * delt;
    
    for (int k = kts+1; k <= kte-1; k++) {
        a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * s_aw[k] + 0.5 * dtz[k] * rhoinv[k] * sd_aw[k];
        b[k] = 1.0 + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) + 0.5 * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]);
        c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1];
        d[k] = thl[k] + tcd[k] * delt + dtz[k] * rhoinv[k] * (s_awthl[k] - s_awthl[k+1]) + dtz[k] * rhoinv[k] * (sd_awthl[k] - sd_awthl[k+1]) + diss_heat[k] * delt + sub_thl[k] * delt + det_thl[k] * delt;
    }
    
    a[kte] = 0.0;
    b[kte] = 1.0;
    c[kte] = 0.0;
    d[kte] = thl[kte];
    tridiag2_cc(kte, a, b, c, d, x);
    for (int k = kts; k <= kte; k++) {
        Dth[k] = (x[k] - thl[k]) / delt;
    }
    
    if (bl_mynn_mixqt > 0) {
        k = kts;
        a[k] = -dtz[k] * khdz[k] * rhoinv[k];
        b[k] = 1.0 + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1];
        c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1];
        d[k] = sqw[k] + dtz[k] * rhosfc * flq * rhoinv[k] + qcd[k] * delt - dtz[k] * rhoinv[k] * s_awqt[k+1] - dtz[k] * rhoinv[k] * sd_awqt[k+1];
        
        for (int k = kts+1; k <= kte-1; k++) {
            a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * s_aw[k] + 0.5 * dtz[k] * rhoinv[k] * sd_aw[k];
            b[k] = 1.0 + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) + 0.5 * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]);
            c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1];
            d[k] = sqw[k] + qcd[k] * delt + dtz[k] * rhoinv[k] * (s_awqt[k] - s_awqt[k+1]) + dtz[k] * rhoinv[k] * (sd_awqt[k] - sd_awqt[k+1]);
        }
        
        a[kte] = 0.0;
        b[kte] = 1.0;
        c[kte] = 0.0;
        d[kte] = sqw[kte];
        tridiag2_cc(kte, a, b, c, d, sqw2);
    } else {
        for (int k = kts; k <= kte; k++) {
            sqw2[k] = sqw[k];
        }
    }
    
    if (bl_mynn_mixqt == 0) {
        if (bl_mynn_cloudmix > 0 && FLAG_QC==1) {
            k = kts;
            a[k] = -dtz[k] * khdz[k] * rhoinv[k];
            b[k] = 1.0 + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1];
            c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1];
            d[k] = sqc[k] + dtz[k] * rhosfc * flqc * rhoinv[k] + qcd[k] * delt - dtz[k] * rhoinv[k] * s_awqc[k+1] - dtz[k] * rhoinv[k] * sd_awqc[k+1] + det_sqc[k] * delt;
            
            for (int k = kts+1; k <= kte-1; k++) {
                a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * s_aw[k] + 0.5 * dtz[k] * rhoinv[k] * sd_aw[k];
                b[k] = 1.0 + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) + 0.5 * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]);
                c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1];
                d[k] = sqc[k] + qcd[k] * delt + dtz[k] * rhoinv[k] * (s_awqc[k] - s_awqc[k+1]) + dtz[k] * rhoinv[k] * (sd_awqc[k] - sd_awqc[k+1]) + det_sqc[k] * delt;
            }
            
            a[kte] = 0.0;
            b[kte] = 1.0;
            c[kte] = 0.0;
            d[kte] = sqc[kte];
            tridiag2_cc(kte, a, b, c, d, sqc2);
        } else {
            for (int k = kts; k <= kte; k++) {
                sqc2[k] = sqc[k];
            }
        }
    }
    
    if (bl_mynn_mixqt == 0) {
        k = kts;
        qvflux = flqv;
        if (qvflux < 0.0) {
            qvflux = std::max(float(qvflux), (std::min(0.9f * float(sqv[kts]) - 1e-8f, 0.0f) / dtz[kts]));
        }
        a[k] = -dtz[k] * khdz[k] * rhoinv[k];
        b[k] = 1.0 + dtz[k] * (khdz[k+1] + khdz[k]) * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1];
        c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1];
        d[k] = sqv[k] + dtz[k] * rhosfc * qvflux * rhoinv[k] + qcd[k] * delt - dtz[k] * rhoinv[k] * s_awqv[k+1] - dtz[k] * rhoinv[k] * sd_awqv[k+1] + sub_sqv[k] * delt + det_sqv[k] * delt;
        
        for (int k = kts+1; k <= kte-1; k++) {
            a[k] = -dtz[k] * khdz[k] * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * s_aw[k] + 0.5 * dtz[k] * rhoinv[k] * sd_aw[k];
            b[k] = 1.0 + dtz[k] * (khdz[k] + khdz[k+1]) * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * (s_aw[k] - s_aw[k+1]) + 0.5 * dtz[k] * rhoinv[k] * (sd_aw[k] - sd_aw[k+1]);
            c[k] = -dtz[k] * khdz[k+1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw[k+1] - 0.5 * dtz[k] * rhoinv[k] * sd_aw[k+1];
            d[k] = sqv[k] + qcd[k] * delt + dtz[k] * rhoinv[k] * (s_awqv[k] - s_awqv[k+1]) + dtz[k] * rhoinv[k] * (sd_awqv[k] - sd_awqv[k+1]) + sub_sqv[k] * delt + det_sqv[k] * delt;
        }
        
        a[kte] = 0.0;
        b[kte] = 1.0;
        c[kte] = 0.0;
        d[kte] = sqv[kte];
        tridiag2_cc(kte, a, b, c, d, sqv2);
    } else {
        for (int k = kts; k <= kte; k++) {
            sqv2[k] = sqv[k];
        }
    }

    delete [] dtz;
    delete [] dfhc;
    delete [] dfmc;
    delete [] delp;
    delete [] sqv2;
    delete [] sqc2;
    delete [] rhoz;
    delete [] khdz;
    delete [] kmdz;
    delete [] rhoinv;
    delete [] sqw2;
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] d;
    delete [] x;

}



void mym_condensation_cc(int kts, int kte, float dx, float dz[], float zw[], float xland,float thl[], 
		float qw[], float qv[], float qc[], float qi[], float qs[],float p[], float exner[], 
		float tsq[], float qsq[], float cov[],float Sh[], float el[], int bl_mynn_cloudpdf,
                float qc_bl1D[], float qi_bl1D[], float cldfra_bl1D[],float PBLH1, float HFX1,
                float Vt[], float Vq[], float th[], float sgm[], float rmo[],int spp_pbl, float rstoch_col[], 
		float ep_2, float ep_3, float xlv, float r_d, float xlvcp, float p608, float tv0, float cpv, 
		float r_v, float cice, float cliq, float cp, float xls, float rcp) {
    int k;
    float t3sq, r3sq, c3sq;
    float qsl, esat, qsat, dqsl, cld0, q1k, qlk, eq1, qll, q2p, pt, rac, qt, t, xl, rsl, cpm, Fng, qww, alpha, beta, bb, ls, wt, wt2, qpct, cld_factor, fac_damp, liq_frac, ql_ice, ql_water, qmq, qsat_tk, q1_rh, rh_hack, dzm1, zsl, maxqc;
    const float qpct_sfc = 0.025;
    const float qpct_pbl = 0.030;
    const float qpct_trp = 0.040;
    const float rhcrit = 0.83;
    const float rhmax = 1.02;
    float erf;
    float dth, dtl, dqw, dzk, els;
    float zagl, damp, PBLH2;
    float cfmax;
    float theta1, theta2, ht1, ht2;
    float qw_pert;
    int k_tropo;
//real(float), dimension(kts:kte) :: alp,a,bet,b,ql,q1,RH
    float* alp = new float[kte-kts]; 
    float* a = new float[kte-kts]; 
    float* bet = new float[kte-kts]; 
    float* b = new float[kte-kts]; 
    float* ql = new float[kte-kts]; 
    float* q1 = new float[kte-kts]; 
    float* rh = new float[kte-kts]; 

    // Obtain an estimate for the tropopause height (k)
    for (k = kte - 3; k >= kts; k--) {
        theta1 = th[k];
        theta2 = th[k + 2];
        ht1 = 44307.692 * (1.0 - pow(p[k] / 101325.0, 0.190));
        ht2 = 44307.692 * (1.0 - pow(p[k + 2] / 101325.0, 0.190));
        if ((((theta2 - theta1) / (ht2 - ht1)) < 10.0 / 1500.0) && (ht1 < 19000.0) && (ht1 > 4000.0)) {
            break;
        }
    }
    k_tropo = std::max(kts + 2, k + 2);
    zagl = 0.0;

    switch (bl_mynn_cloudpdf) {
        case 0: // ORIGINAL MYNN PARTIAL-CONDENSATION SCHEME
            for (k = kts; k < kte; k++) {
                t = th[k] * exner[k];
                esat = esat_blend_cc(t);
                qsl = ep_2 * esat / std::max(1e-4f, (p[k] - ep_3 * esat));
                dqsl = qsl * ep_2 * xlv / (r_d * pow(t, 2));
                alp[k] = 1.0 / (1.0 + dqsl * xlvcp);
                bet[k] = dqsl * exner[k];
                t3sq = std::max(tsq[k], 0.0f);
                r3sq = std::max(qsq[k], 0.0f);
                c3sq = cov[k];
                c3sq = std::copysign(std::min(std::abs(c3sq), std::sqrt(t3sq * r3sq)), c3sq);
                r3sq = r3sq + bet[k] * bet[k] * t3sq - 2.0 * bet[k] * c3sq;
                qmq = qw[k] - qsl;
                sgm[k] = std::sqrt(std::max(r3sq, 1.0e-10f));
                q1[k] = qmq / sgm[k];
                cldfra_bl1D[k] = 0.5 * (1.0 + std::erf(q1[k] * rr2));
                q1k = q1[k];
                eq1 = rrp * std::exp(-0.5 * q1k * q1k);
                qll = std::max(cldfra_bl1D[k] * q1k + eq1, 0.0f);
                ql[k] = alp[k] * sgm[k] * qll;
                liq_frac = std::min(1.0f, std::max(0.0f, (t - 240.0f) / 29.0f));
                qc_bl1D[k] = liq_frac * ql[k];
                qi_bl1D[k] = (1.0 - liq_frac) * ql[k];
                q2p = xlvcp / exner[k];
                pt = thl[k] + q2p * ql[k];
                qt = 1.0 + p608 * qw[k] - (1.0 + p608) * (qc_bl1D[k] + qi_bl1D[k]) * cldfra_bl1D[k];
                rac = alp[k] * (cldfra_bl1D[k] - qll * eq1) * (q2p * qt - (1.0 + p608) * pt);
                Vt[k] = qt - 1.0 - rac * bet[k];
                Vq[k] = p608 * pt - tv0 + rac;
            }
            break;
        case 1:
        case -1: // ALTERNATIVE FORM (Nakanishi & Niino 2004 BLM, eq. B6, and Kuwano-Yoshida et al. 2010 QJRMS, eq. 7)
            for (k = kts; k < kte; k++) {
                t = th[k] * exner[k];
                esat = esat_blend_cc(t);
                qsl = ep_2 * esat / std::max(1e-4f, (p[k] - ep_3 * esat));
                dqsl = qsl * ep_2 * xlv / (r_d * pow(t, 2));
                alp[k] = 1.0 / (1.0 + dqsl * xlvcp);
                bet[k] = dqsl * exner[k];
                if (k == kts) {
                    dzk = 0.5 * dz[k];
                } else {
                    dzk = dz[k];
                }
                dth = 0.5 * (thl[k + 1] + thl[k]) - 0.5 * (thl[k] + thl[std::max(k - 1, kts)]);
                dqw = 0.5 * (qw[k + 1] + qw[k]) - 0.5 * (qw[k] + qw[std::max(k - 1, kts)]);
                sgm[k] = std::sqrt(std::max(float((pow(alp[k], 2) * std::max(pow(el[k], 2), 0.1) * b2 * std::max(Sh[k], 0.03f)) / 4.0f * pow((dqw / dzk - bet[k] * (dth / dzk)), 2)), 1.0e-10f));
                qmq = qw[k] - qsl;
                q1[k] = qmq / sgm[k];
                cldfra_bl1D[k] = 0.5 * (1.0 + std::erf(q1[k] * rr2));
                q1k = q1[k];
                eq1 = rrp * std::exp(-0.5 * q1k * q1k);
                qll = std::max(cldfra_bl1D[k] * q1k + eq1, 0.0f);
                ql[k] = alp[k] * sgm[k] * qll;
                liq_frac = std::min(1.0f, std::max(0.0f, (t - 240.0f) / 29.0f));
                qc_bl1D[k] = liq_frac * ql[k];
                qi_bl1D[k] = (1.0 - liq_frac) * ql[k];
                q2p = xlvcp / exner[k];
                pt = thl[k] + q2p * ql[k];
                qt = 1.0 + p608 * qw[k] - (1.0 + p608) * (qc_bl1D[k] + qi_bl1D[k]) * cldfra_bl1D[k];
                rac = alp[k] * (cldfra_bl1D[k] - qll * eq1) * (q2p * qt - (1.0 + p608) * pt);
                Vt[k] = qt - 1.0 - rac * bet[k];
                Vq[k] = p608 * pt - tv0 + rac;
            }
            break;
        case 2:
        case -2: // Diagnostic statistical scheme of Chaboureau and Bechtold (2002), JAS
                 // but with use of higher-order moments to estimate sigma
            PBLH2 = std::max(10.0f, PBLH1);
            zagl = 0.0f;
            dzm1 = 0.0f;
            for (k = kts; k < kte; k++) {
                zagl += 0.5f * (dz[k] + dzm1);
                dzm1 = dz[k];
                t = th[k] * exner[k];
                xl = xl_blend_cc(t,xlv,xls,cpv,cliq,cice);
                qsat_tk = qsat_blend_cc(t, p[k]);
                rh[k] = std::max(std::min(rhmax, qw[k] / std::max(1e-10f, qsat_tk)), 0.001f);
                dqsl = qsat_tk * ep_2 * xlv / (r_d * pow(t, 2));
                alp[k] = 1.0f / (1.0f + dqsl * xlvcp);
                bet[k] = dqsl * exner[k];
                rsl = xl * qsat_tk / (r_v * pow(t, 2));
                cpm = cp + qw[k] * cpv;
                a[k] = 1.0f / (1.0f + xl * rsl / cpm);
                b[k] = a[k] * rsl;
                qw_pert = qw[k] + qw[k] * 0.5f * rstoch_col[k] * spp_pbl;
                qmq = qw_pert - qsat_tk;
                r3sq = std::max(qsq[k], 0.0f);
                sgm[k] = std::sqrt(r3sq);
                sgm[k] = std::min(sgm[k], qsat_tk * 0.666f);
                wt = std::max(500.0f - std::max(dz[k] - 100.0f, 0.0f), 0.0f) / 500.0f;
                sgm[k] += sgm[k] * 0.2f * (1.0f - wt);
                qpct = qpct_pbl * wt + qpct_trp * (1.0f - wt);
                qpct = std::min(qpct, std::max(qpct_sfc, qpct_pbl * zagl / 500.0f));
                sgm[k] = std::max(sgm[k], qsat_tk * qpct);
                q1[k] = qmq / sgm[k];
                q1k = q1[k];
                eq1 = rrp * std::exp(-0.5f * q1k * q1k);
                qll = std::max(cldfra_bl1D[k] * q1k + eq1, 0.0f);
                ql[k] = alp[k] * sgm[k] * qll;
                liq_frac = std::min(1.0f, std::max(0.0f, (t - tice) / (tliq - tice)));
                qc_bl1D[k] = liq_frac * ql[k];
                qi_bl1D[k] = (1.0f - liq_frac) * ql[k];
                if (cldfra_bl1D[k] < 0.001f) {
                    ql_ice = 0.0f;
                    ql_water = 0.0f;
                    cldfra_bl1D[k] = 0.0f;
                }
                if ((qi[k] + qs[k]) > 1e-9f && zagl > PBLH2) {
                    rh_hack = std::min(rhmax, rhcrit + wt2 * 0.045f * (9.0f + std::log10(qi[k] + qs[k])));
                    rh[k] = std::max(rh[k], rh_hack);
                    q1_rh = -3.0f + 3.0f * (rh[k] - rhcrit) / (1.0f - rhcrit);
                    q1[k] = std::max(q1_rh, q1[k]);
                }
                if (qc[k] > 1e-6f && zagl > PBLH2) {
                    rh_hack = std::min(rhmax, rhcrit + wt2 * 0.08f * (6.0f + std::log10(qc[k])));
                    rh[k] = std::max(rh[k], rh_hack);
                    q1_rh = -3.0f + 3.0f * (rh[k] - rhcrit) / (1.0f - rhcrit);
                    q1[k] = std::max(q1_rh, q1[k]);
                }
                q1k = q1[k];
                cldfra_bl1D[k] = std::max(0.0f, std::min(1.0f, 0.5f + 0.36f * std::atan(1.8f * (q1[k] + 0.2f))));
                maxqc = std::max(qw[k] - qsat_tk, 0.0f);
                if (q1k < 0.0) {
                    ql_water = sgm[k] * std::exp(1.2f * q1k - 1.0f);
                    ql_ice = sgm[k] * std::exp(1.2f * q1k - 1.0f);
                } else if (q1k > 2.0) {
                    ql_water = std::min(sgm[k] * q1k, maxqc);
                    ql_ice = sgm[k] * q1k;
                } else {
                    ql_water = std::min(float(sgm[k] * (std::exp(-1.0f) + 0.66f * q1k + 0.086f * pow(q1k, 2.0f))), maxqc);
                    ql_ice = sgm[k] * (std::exp(-1.0) + 0.66f * q1k + 0.086f * pow(q1k, 2));
                }
                if (cldfra_bl1D[k] < 0.001f) {
                    ql_ice = 0.0f;
                    ql_water = 0.0f;
                    cldfra_bl1D[k] = 0.0f;
                }
                liq_frac = std::min(1.0f, std::max(0.0f, (t - tice) / (tliq - tice)));
                qc_bl1D[k] = liq_frac * ql_water;
                qi_bl1D[k] = (1.0f - liq_frac) * ql_ice;
                if (k >= k_tropo) {
                    cldfra_bl1D[k] = 0.0f;
                    qc_bl1D[k] = 0.0f;
                    qi_bl1D[k] = 0.0f;
                }
                q2p = xlvcp / exner[k];
                pt = thl[k] + q2p * ql[k];
                qt = 1.0f + p608 * qw[k] - (1.0f + p608) * (qc_bl1D[k] + qi_bl1D[k]) * cldfra_bl1D[k];
                rac = alp[k] * (cldfra_bl1D[k] - qll * eq1) * (q2p * qt - (1.0f + p608) * pt);
                Vt[k] = qt - 1.0f - rac * bet[k];
                Vq[k] = p608 * pt - tv0 + rac;
                fac_damp = std::min(zagl * 0.0025f, 1.0f);
                cld_factor = 1.0 + fac_damp * std::min(std::pow(std::max(0.0f, (rh[k] - 0.92f)) / 0.145f, 2.0f), 0.37f);
                cldfra_bl1D[k] = std::min(1.0f, cld_factor * cldfra_bl1D[k]);
            }
            break;
    }

    ql[kte] = ql[kte - 1];
    Vt[kte] = Vt[kte - 1];
    Vq[kte] = Vq[kte - 1];
    qc_bl1D[kte] = 0.0;
    qi_bl1D[kte] = 0.0;
    cldfra_bl1D[kte] = 0.0;
}


//===============================================================
// ===================================================================
// This is the downdraft mass flux scheme - analogus to edmf_JPL but
// flipped updraft to downdraft. This scheme is currently only tested
// for Stratocumulus cloud conditions. For a detailed desctiption of the
// model, see paper.
void DDMF_JPL_cc(int kts, int kte, float dt, std::vector<float> zw, std::vector<float> dz, std::vector<float> p,
              std::vector<float> u, std::vector<float> v, std::vector<float> th, std::vector<float> thl, std::vector<float> thv, 
	      std::vector<float> tk,std::vector<float> qt, std::vector<float> qv, std::vector<float> qc, std::vector<float> 
	      rho, std::vector<float> exner,float ust, float wthl, float wqt, float pblh, int kpbl,
              std::vector<float>& edmf_a_dd, std::vector<float>& edmf_w_dd, std::vector<float>& edmf_qt_dd,
              std::vector<float>& edmf_thl_dd, std::vector<float>& edmf_ent_dd, std::vector<float>& edmf_qc_dd,
              std::vector<float>& sd_aw, std::vector<float>& sd_awthl, std::vector<float>& sd_awqt,
              std::vector<float>& sd_awqv, std::vector<float>& sd_awqc, std::vector<float>& sd_awu,
              std::vector<float>& sd_awv, std::vector<float>& sd_awqke,
              std::vector<float> qc_bl1d, std::vector<float> cldfra_bl1d,
              std::vector<float> rthraten, float svp1, float grav, float onethird, float p1000mb, 
	      float rcp, float xlvcp, float cp, float rvovrd ) {
    int ndown = 5;
    std::vector<int> DD_initK(ndown);
    std::vector<float> randNum(ndown);
    std::vector<std::vector<float>> DOWNW(kte + 1, std::vector<float>(ndown));
    std::vector<std::vector<float>> DOWNTHL(kte + 1, std::vector<float>(ndown));
    std::vector<std::vector<float>> DOWNQT(kte + 1, std::vector<float>(ndown));
    std::vector<std::vector<float>> DOWNQC(kte + 1, std::vector<float>(ndown));
    std::vector<std::vector<float>> DOWNA(kte + 1, std::vector<float>(ndown));
    std::vector<std::vector<float>> DOWNU(kte + 1, std::vector<float>(ndown));
    std::vector<std::vector<float>> DOWNV(kte + 1, std::vector<float>(ndown));
    std::vector<std::vector<float>> DOWNTHV(kte + 1, std::vector<float>(ndown));
    std::vector<std::vector<float>> ENT(kte + 1, std::vector<float>(ndown));
    std::vector<std::vector<int>> ENTi(kte + 1, std::vector<int>(ndown));
    int k, I, ki, kminrad, qlTop, p700_ind, qlBase;
    float wthv, wstar, qstar, thstar, sigmaW, sigmaQT, sigmaTH, z0, pwmin, pwmax, wmin, wmax, wlv, wtv, went, mindownw;
    float B, QTn, THLn, THVn, QCn, Un, Vn, QKEn, Wn2, Wn, THVk, Pk, EntEXP, EntW, beta_dm, EntExp_M, rho_int;
    float jump_thetav, jump_qt, jump_thetal, refTHL, refTHV, refQT;
    float minrad, zminrad, radflux, F0, wst_rad, wst_dd;
    bool cloudflg;
    float sigq, xl, rsl, cpm, a, mf_cf, diffqt, Fng, qww, alpha, beta, bb, f, pt, t, q2p, b9, satvp, rhgrid;
    float Wa = 1.0, Wb = 1.5, Z00 = 100.0, BCOEFF = 0.2;
    float L0 = 80, ENT0 = 0.2;
    float dp, dl, Adn;
    int debug_mf = 0;
    dl = (1000.0 - 500.0) / ndown;
    pwmin = -3.0;
    pwmax = -1.0;
    DOWNW.assign(kte + 1, std::vector<float>(ndown, 0.0));
    DOWNTHL.assign(kte + 1, std::vector<float>(ndown, 0.0));
    DOWNTHV.assign(kte + 1, std::vector<float>(ndown, 0.0));
    DOWNQT.assign(kte + 1, std::vector<float>(ndown, 0.0));
    DOWNQC.assign(kte + 1, std::vector<float>(ndown, 0.0));
    DOWNA.assign(kte + 1, std::vector<float>(ndown, 0.0));
    DOWNU.assign(kte + 1, std::vector<float>(ndown, 0.0));
    DOWNV.assign(kte + 1, std::vector<float>(ndown, 0.0));
    ENT.assign(kte + 1, std::vector<float>(ndown, 0.0));
    DD_initK.assign(ndown, 0);
    edmf_a_dd.assign(kte + 1, 0.0);
    edmf_w_dd.assign(kte + 1, 0.0);
    edmf_qt_dd.assign(kte + 1, 0.0);
    edmf_thl_dd.assign(kte + 1, 0.0);
    edmf_ent_dd.assign(kte + 1, 0.0);
    edmf_qc_dd.assign(kte + 1, 0.0);
    sd_aw.assign(kte + 1, 0.0);
    sd_awthl.assign(kte + 1, 0.0);
    sd_awqt.assign(kte + 1, 0.0);
    sd_awqv.assign(kte + 1, 0.0);
    sd_awqc.assign(kte + 1, 0.0);
    sd_awu.assign(kte + 1, 0.0);
    sd_awv.assign(kte + 1, 0.0);
    sd_awqke.assign(kte + 1, 0.0);
    for (int i = 0; i < ndown; i++) {
        DD_initK[i] = qlTop;
    }
    F0 = 0.0;
    for (int i = 0; i < qlTop; i++) {
        radflux = rthraten[i] * exner[i];
        radflux = radflux * cp / grav * (p[i] - p[i + 1]);
        if (radflux < 0.0) {
            F0 = abs(radflux) + F0;
        }
    }
    F0 = std::max(F0, 1.0f);
    Adn = std::min(0.05 + F0 * 0.001, 0.3);
    cloudflg = false;
    minrad = 100.0;
    kminrad = kpbl;
    zminrad = pblh;
    qlTop = 1;
    qlBase = 1;
    wthv = wthl + svp1 * wqt;
    for (int i = std::max(3, kpbl - 2); i <= kpbl + 3; i++) {
        if (qc[i] > 1.0e-6 && cldfra_bl1d[i] > 0.5) {
            cloudflg = true;
            qlTop = i;
        }
    }
    for (int i = qlTop; i >= kts; i--) {
        if (qc[i] > 1E-6) {
            qlBase = i;
        }
    }
    qlBase = (qlTop + qlBase) / 2;
    for (int i = 0; i < ndown; i++) {
        DD_initK[i] = qlTop;
    }
    if (cloudflg) {

    p700_ind = 0;
    float min_value = p[0];
    for (int i = 1; i < p.size(); ++i) {
	float pval=abs(p[i]-70000.0);
        if (pval < min_value) {
            p700_ind = i;
        }
    }

        //p700_ind = minloc(abs(p - 70000.0), 1.0f);
        jump_thetav = thv[p700_ind] - thv[1] - (thv[p700_ind] - thv[qlTop + 3]) / (zw[p700_ind] - zw[qlTop + 3]) * (zw[p700_ind] - zw[qlTop]);
        jump_qt = qc[p700_ind] + qv[p700_ind] - qc[1] - qv[1];
        jump_thetal = thl[p700_ind] - thl[1] - (thl[p700_ind] - thl[qlTop + 3]) / (zw[p700_ind] - zw[qlTop + 3]) * (zw[p700_ind] - zw[qlTop]);
        refTHL = thl[qlTop];
        refTHV = thv[qlTop];
        refQT = qt[qlTop];
        wst_rad = pow(grav * zw[qlTop] * F0 / (refTHL * rho[qlTop] * cp), 0.333);
        wst_rad = std::max(wst_rad, 0.1f);
        wstar = std::max(0.0, pow(grav / thv[1] * wthv * pblh, onethird));
        went = thv[1] / (grav * jump_thetav * zw[qlTop]) * (0.15 * (pow(wstar, 3) + 5 * pow(ust, 3)) + 0.35 * pow(wst_rad, 3));
        qstar = abs(went * jump_qt / wst_rad);
        thstar = F0 / (rho[qlTop] * cp * wst_rad) - went * jump_thetav / wst_rad;
        wst_dd = pow(0.15 * (pow(wstar, 3) + 5 * pow(ust, 3)) + 0.35 * pow(wst_rad, 3), 0.333);
	std::cout << "qstar=" << qstar << " thstar=" << thstar << " wst_dd=" << wst_dd << std::endl;
	std::cout << "F0=" << F0 << " wst_rad=" << wst_rad << " jump_thv=" << jump_thetav << std::endl;
	std::cout << "entrainment velocity=" << went << std::endl;
        sigmaW = 0.2 * wst_dd;
        sigmaQT = 40 * qstar;
        sigmaTH = 1.0 * thstar;
        wmin = sigmaW * pwmin;
        wmax = sigmaW * pwmax;
        for (int i = 0; i < ndown; i++) {
            ki = DD_initK[i];
            wlv = wmin + (wmax - wmin) / ndown * (i - 1);
            wtv = wmin + (wmax - wmin) / ndown * i;
            DOWNW[ki][i] = wlv;
            DOWNA[ki][i] = Adn / ndown;
            DOWNU[ki][i] = (u[ki - 1] * dz[ki] + u[ki] * dz[ki - 1]) / (dz[ki] + dz[ki - 1]);
            DOWNV[ki][i] = (v[ki - 1] * dz[ki] + v[ki] * dz[ki - 1]) / (dz[ki] + dz[ki - 1]);
            refTHL = (thl[ki - 1] * dz[ki] + thl[ki] * dz[ki - 1]) / (dz[ki] + dz[ki - 1]);
            refTHV = (thv[ki - 1] * dz[ki] + thv[ki] * dz[ki - 1]) / (dz[ki] + dz[ki - 1]);
            refQT = (qt[ki - 1] * dz[ki] + qt[ki] * dz[ki - 1]) / (dz[ki] + dz[ki - 1]);
            DOWNQC[ki][i] = (qc[ki - 1] * dz[ki] + qc[ki] * dz[ki - 1]) / (dz[ki] + dz[ki - 1]);
            DOWNQT[ki][i] = refQT;
            DOWNTHV[ki][i] = refTHV + 0.01 * DOWNW[ki][i] * sigmaTH / sigmaW;
            DOWNTHL[ki][i] = refTHL + 0.01 * DOWNW[ki][i] * sigmaTH / sigmaW;
        
        for (int k = DD_initK[i] - 1; k >= kts + 1; k--) {
            wmin = 0.3 + dp * 0.0005;
            ENT[k][i] = 0.33 / (std::min(std::max(-1.0f * DOWNW[k + 1][i], wmin), 0.9f) * dp);
            EntEXP = ENT[k][i] * dz[k];
            EntExp_M = ENT[k][i] * 0.333 * dz[k];
            QTn = DOWNQT[k + 1][i] * (1.0 - EntEXP) + qt[k] * EntEXP;
            THLn = DOWNTHL[k + 1][i] * (1.0 - EntEXP) + thl[k] * EntEXP;
            Un = DOWNU[k + 1][i] * (1.0 - EntEXP) + u[k] * EntExp_M;
            Vn = DOWNV[k + 1][i] * (1.0 - EntEXP) + v[k] * EntExp_M;
            Pk = (p[k - 1] * dz[k] + p[k] * dz[k - 1]) / (dz[k] + dz[k - 1]);
            condensation_edmf(QTn, THLn, Pk, zw[k], THVn, QCn,p1000mb,rcp,xlvcp,rvovrd);
            THVk = (thv[k - 1] * dz[k] + thv[k] * dz[k - 1]) / (dz[k] + dz[k - 1]);
            B = grav * (THVn / THVk - 1.0);
            EntW = EntEXP;
            mindownw = std::min(DOWNW[k + 1][i], -0.2f);
            Wn = DOWNW[k + 1][i] + (-2.0 * ENT[k][i] * DOWNW[k + 1][i] - BCOEFF * B / mindownw) * std::min(dz[k], 250.0f);
            if (Wn < DOWNW[k + 1][i] - std::min(1.25 * dz[k] / 200.0, -2.0)) {
                Wn = DOWNW[k + 1][i] - std::min(1.25 * dz[k] / 200.0, -2.0);
            }
            if (Wn > DOWNW[k + 1][i] + std::min(1.25 * dz[k] / 200.0, 2.0)) {
                Wn = DOWNW[k + 1][i] + std::min(1.25 * dz[k] / 200.0, 2.0);
            }
            Wn = std::max(std::min(Wn, 0.0f), -3.0f);
            if (Wn < 0.0) {
                DOWNW[k][i] = Wn;
                DOWNTHV[k][i] = THVn;
                DOWNTHL[k][i] = THLn;
                DOWNQT[k][i] = QTn;
                DOWNQC[k][i] = QCn;
                DOWNU[k][i] = Un;
                DOWNV[k][i] = Vn;
                DOWNA[k][i] = DOWNA[k + 1][i];
            } 
	    else {
                if (DD_initK[i] - k < 2) {
                    DOWNW.assign(kte + 1, std::vector<float>(ndown, 0.0));
                    DOWNTHV.assign(kte + 1, std::vector<float>(ndown, 0.0));
                    DOWNTHL.assign(kte + 1, std::vector<float>(ndown, 0.0));
                    DOWNQT.assign(kte + 1, std::vector<float>(ndown, 0.0));
                    DOWNQC.assign(kte + 1, std::vector<float>(ndown, 0.0));
                    DOWNU.assign(kte + 1, std::vector<float>(ndown, 0.0));
                    DOWNV.assign(kte + 1, std::vector<float>(ndown, 0.0));
                    }
                break;
                }
            }
        }
    }
    DOWNW[0].assign(ndown, 0.0);
    DOWNA[0].assign(ndown, 0.0);
    for (int k = qlTop; k <= kts; k++) {
        for (int i = 0; i < ndown; i++) {
            edmf_a_dd[k] += DOWNA[k - 1][i];
            edmf_w_dd[k] += DOWNA[k - 1][i] * DOWNW[k - 1][i];
            edmf_qt_dd[k] += DOWNA[k - 1][i] * DOWNQT[k - 1][i];
            edmf_thl_dd[k] += DOWNA[k - 1][i] * DOWNTHL[k - 1][i];
            edmf_ent_dd[k] += DOWNA[k - 1][i] * ENT[k - 1][i];
            edmf_qc_dd[k] += DOWNA[k - 1][i] * DOWNQC[k - 1][i];
        }
        if (edmf_a_dd[k] > 0.0) {
            edmf_w_dd[k] /= edmf_a_dd[k];
            edmf_qt_dd[k] /= edmf_a_dd[k];
            edmf_thl_dd[k] /= edmf_a_dd[k];
            edmf_ent_dd[k] /= edmf_a_dd[k];
            edmf_qc_dd[k] /= edmf_a_dd[k];
        }
    }
    for (int k = kts; k <= qlTop; k++) {
        rho_int = (rho[k] * dz[k + 1] + rho[k + 1] * dz[k]) / (dz[k + 1] + dz[k]);
        for (int i = 0; i < ndown; i++) {
            sd_aw[k] += rho_int * DOWNA[k][i] * DOWNW[k][i];
            sd_awthl[k] += rho_int * DOWNA[k][i] * DOWNW[k][i] * DOWNTHL[k][i];
            sd_awqt[k] += rho_int * DOWNA[k][i] * DOWNW[k][i] * DOWNQT[k][i];
            sd_awqc[k] += rho_int * DOWNA[k][i] * DOWNW[k][i] * DOWNQC[k][i];
            sd_awu[k] += rho_int * DOWNA[k][i] * DOWNW[k][i] * DOWNU[k][i];
            sd_awv[k] += rho_int * DOWNA[k][i] * DOWNW[k][i] * DOWNV[k][i];
        }
        sd_awqv[k] = sd_awqt[k] - sd_awqc[k];
    }
}

// Assuming float is equivalent to float or float. Adjust as necessary.
void topdown_cloudrad_cc(int kts, int kte, const std::vector<float>& dz1, const std::vector<float>& zw, float fltv, float xland, int kpbl, float PBLH, const std::vector<float>& sqc, const std::vector<float>& sqi, const std::vector<float>& sqw, const std::vector<float>& thl, const std::vector<float>& th1, const std::vector<float>& ex1, const std::vector<float>& p1, const std::vector<float>& rho1, const std::vector<float>& thetav, const std::vector<float>& cldfra_bl1D, const std::vector<float>& rthraten, float& maxKHtopdown, std::vector<float>& KHtopdown, std::vector<float>& TKEprodTD) {
    // Constants
    const float pfac = 2.0, zfmin = 0.01, phifac = 8.0;
    const float grav = 9.81, cp = 1004.0, xlv = 2.5e6, xlvcp = xlv / cp, r_d = 287.0, ep_2 = 0.622, p608 = 0.608, karman = 0.4;
    const float twothirds = 2.0 / 3.0, onethird = 1.0 / 3.0;

    // Local variables
    std::vector<float> zfac(kte - kts + 1), wscalek2(kte - kts + 1), zfacent(kte - kts + 1);
    float bfx0, wm2 = 0, wm3, bfxpbl, dthvx, tmp1;
    float temps, templ, zl1, wstar3_2;
    float ent_eff, radsum, radflux, we, rcldb, rvls, minrad = 100., zminrad;
    int k, kk, kminrad = kpbl;
    bool cloudflg = false;

    KHtopdown.assign(kte - kts + 1, 0.0);
    TKEprodTD.assign(kte - kts + 1, 0.0);
    maxKHtopdown = 0.0;

    // Check for stratocumulus-topped boundary layers
    for (kk = std::max(1, kpbl - 2); kk <= kpbl + 3; ++kk) {
        if (sqc[kk - kts] > 1.e-6 || sqi[kk - kts] > 1.e-6 || cldfra_bl1D[kk - kts] > 0.5) {
            cloudflg = true;
        }
        if (rthraten[kk - kts] < minrad) {
            minrad = rthraten[kk - kts];
            kminrad = kk;
            zminrad = zw[kk - kts] + 0.5 * dz1[kk - kts];
        }
    }
    if (std::max(kminrad, kpbl) < 2) cloudflg = false;

    if (cloudflg) {
        zl1 = dz1[0]; // Assuming kts is 1-based index in Fortran and adjusted to 0-based in C++
        k = std::max(kpbl - 1, kminrad - 1);
        templ = thl[k - kts] * ex1[k - kts];
        rvls = 100. * 6.112 * std::exp(17.67 * (templ - 273.16) / (templ - 29.65)) * (ep_2 / p1[k + 1 - kts]);
        temps = templ + (sqw[k - kts] - rvls) / (cp / xlv + ep_2 * xlv * rvls / (r_d * std::pow(templ, 2)));
        rvls = 100. * 6.112 * std::exp(17.67 * (temps - 273.15) / (temps - 29.65)) * (ep_2 / p1[k + 1 - kts]);
        rcldb = std::max(sqw[k - kts] - rvls, 0.0f);
        dthvx = (thl[k + 2 - kts] + th1[k + 2 - kts] * p608 * sqw[k + 2 - kts]) - (thl[k - kts] + th1[k - kts] * p608 * sqw[k - kts]);
        dthvx = std::max(dthvx, 0.1f);
        tmp1 = xlvcp * rcldb / (ex1[k - kts] * dthvx);
        ent_eff = 0.2 + 0.2 * 8. * tmp1;
        radsum = 0.0;
        for (kk = std::max(1, kpbl - 3); kk <= kpbl + 3; ++kk) {
            radflux = rthraten[kk - kts] * ex1[kk - kts]; // Converts theta/s to temp/s
            radflux = radflux * cp / grav * (p1[kk - kts] - p1[kk + 1 - kts]); // Converts temp/s to W/m^2
            if (radflux < 0.0) radsum = std::abs(radflux) + radsum;
        }
        if ((xland - 1.5) >= 0) { // WATER
            radsum = std::min(radsum, 90.0f);
            bfx0 = std::max(radsum / rho1[k - kts] / cp, 0.0f);
        } else { // LAND
            radsum = std::min(0.25 * radsum, 30.0); // Practically turn off over land
            bfx0 = std::max(radsum / rho1[k - kts] / cp - std::max(fltv, 0.0f), 0.0f);
        }
        wm3 = grav / thetav[k - kts] * bfx0 * std::min(PBLH, 1500.f); // This is wstar3
        wm2 = wm2 + std::pow(wm3, twothirds);
        bfxpbl = -ent_eff * bfx0;
        dthvx = std::max(thetav[k + 1 - kts] - thetav[k - kts], 0.1f);
        we = std::max(bfxpbl / dthvx, -std::sqrt(std::pow(wm3, twothirds)));
        for (kk = kts; kk <= kpbl + 3; ++kk) {
            zfac[kk - kts] = std::min(std::max((1.f - (zw[kk + 1 - kts] - zl1) / (zminrad - zl1)), zfmin), 1.0f);
            zfacent[kk - kts] = 10. * std::max((zminrad - zw[kk + 1 - kts]) / zminrad, 0.0f) * std::pow((1. - zfac[kk - kts]), 3);
            wscalek2[kk - kts] = std::pow((phifac * karman * wm3 * (zfac[kk - kts])), onethird);
            KHtopdown[kk - kts] = wscalek2[kk - kts] * karman * (zminrad - zw[kk + 1 - kts]) * std::pow((1. - zfac[kk - kts]), 3); // pfac
            KHtopdown[kk - kts] = std::max(KHtopdown[kk - kts], 0.0f);
            TKEprodTD[kk - kts] = 2. * ent_eff * wm3 / std::max(PBLH, 100.f) * zfacent[kk - kts];
            TKEprodTD[kk - kts] = std::max(TKEprodTD[kk - kts], 0.0f);
        }
    }
    maxKHtopdown = *std::max_element(KHtopdown.begin(), KHtopdown.end());
}

void scale_aware_cc(float dx, float pbl1, float& Psig_bl, float& Psig_shcu) {
    float dxdh;
    Psig_bl = 1.0f;
    Psig_shcu = 1.0f;
    dxdh = std::max(2.5f * dx, 10.0f) / std::min(pbl1, 3000.0f);
    
    // New form to preserve parameterized mixing - only down 5% at dx = 750 m
    Psig_bl = ((dxdh * dxdh) + 0.106f * std::pow(dxdh, 0.667f)) / ((dxdh * dxdh) + 0.066f * std::pow(dxdh, 0.667f) + 0.071f);
    
    // Assume a 500 m cloud depth for shallow-cu clouds
    dxdh = std::max(2.5f * dx, 10.0f) / std::min(pbl1 + 500.0f, 3500.0f);
    
    // Hyeyum Hailey Shin and Song-You Hong 2013, TKE in entrainment zone
    Psig_shcu = ((dxdh * dxdh) + 0.145f * std::pow(dxdh, 0.667f)) / ((dxdh * dxdh) + 0.172f * std::pow(dxdh, 0.667f) + 0.170f);
    
    // Clamping Psig_bl and Psig_shcu to [0, 1]
    Psig_bl = std::max(0.0f, std::min(Psig_bl, 1.0f));
    Psig_shcu = std::max(0.0f, std::min(Psig_shcu, 1.0f));
}

// ==================================================================
//>\ingroup gsd_mynn_edmf
// This subroutine calculates hybrid diagnotic boundary-layer height (PBLH).
//
// NOTES ON THE PBLH FORMULATION: The 1.5-theta-increase method defines
//PBL heights as the level at.
//which the potential temperature first exceeds the minimum potential.
//temperature within the boundary layer by 1.5 K. When applied to.
//observed temperatures, this method has been shown to produce PBL-
//height estimates that are unbiased relative to profiler-based.
//estimates (Nielsen-Gammon et al. 2008 \cite Nielsen_Gammon_2008).
// However, their study did not
//include LLJs. Banta and Pichugina (2008) \cite Pichugina_2008  show that a TKE-based.
//threshold is a good estimate of the PBL height in LLJs. Therefore,
//a hybrid definition is implemented that uses both methods, weighting
//the TKE-method more during stable conditions (PBLH < 400 m).
//A variable tke threshold (TKEeps) is used since no hard-wired
//value could be found to work best in all conditions.
//>\section gen_get_pblh  GSD get_pblh General Algorithm
//> @{
void GET_PBLH_CC(int KTS, int KTE, float& zi, float landsea, const std::vector<float>& thetav1D, const std::vector<float>& qke1D, const std::vector<float>& zw1D, const std::vector<float>& dz1D, int& kzi) {
    // Constants
    const float sbl_lim = 200.0;
    const float sbl_damp = 400.0;

    // Local variables
    float PBLH_TKE, qtke, qtkem1, maxqke, TKEeps, minthv, delt_thv;
    int kthv, ktke;

    // Initialize KPBL (kzi)
    kzi = 2;

    // Find min thetav in the lowest 200 m AGL
    kthv = 1;
    minthv = 9E9;
    for (int k = KTS + 1; k <= KTE && zw1D[k - KTS] <= 200.; ++k) {
        if (minthv > thetav1D[k - KTS]) {
            minthv = thetav1D[k - KTS];
            kthv = k;
        }
    }

    // Find thetav-based PBLH (best for daytime)
    zi = 0.0;
    delt_thv = (landsea - 1.5) >= 0 ? 1.0 : 1.25;

    for (int k = kthv + 1; k < KTE; ++k) {
        if (thetav1D[k - KTS] >= (minthv + delt_thv)) {
            zi = zw1D[k - KTS] - dz1D[k - 1 - KTS] * std::min((thetav1D[k - KTS] - (minthv + delt_thv)) / std::max(thetav1D[k - KTS] - thetav1D[k - 1 - KTS], 1.e-6f), 1.0f);
            break;
        }
        if (k == KTE - 1) zi = zw1D[KTS + 1 - KTS]; // Exit safeguard
    }

    // For stable boundary layers, use TKE method to complement the thetav-based definition
    PBLH_TKE = 0.0;
    maxqke = std::max(qke1D[KTS - KTS], 0.0f);
    TKEeps = maxqke / 40.0;
    TKEeps = std::max(TKEeps, 0.02f);

    for (int k = KTS + 1; k < KTE; ++k) {
        qtke = std::max(qke1D[k - KTS] / 2.0, 0.0);
        qtkem1 = std::max(qke1D[k - 1 - KTS] / 2.0, 0.0);
        if (qtke <= TKEeps) {
            PBLH_TKE = zw1D[k - KTS] - dz1D[k - 1 - KTS] * std::min((TKEeps - qtke) / std::max(qtkem1 - qtke, 1.0e-6f), 1.0f);
            PBLH_TKE = std::max(PBLH_TKE, zw1D[KTS + 1 - KTS]);
            break;
        }
        if (k == KTE - 1) PBLH_TKE = zw1D[KTS + 1 - KTS]; // Exit safeguard
    }

    // Limit PBLH_TKE to not exceed the thetav-based PBL height +/- 350 m
    PBLH_TKE = std::min(PBLH_TKE, zi + 350.0f);
    PBLH_TKE = std::max(PBLH_TKE, std::max(zi - 350.0f, 10.0f));

    float wt = 0.5 * std::tanh((zi - sbl_lim) / sbl_damp) + 0.5;
    if (maxqke > 0.05) {
        zi = PBLH_TKE * (1.0 - wt) + zi * wt;
    }

    // Compute KPBL (kzi)
    for (int k = KTS + 1; k < KTE; ++k) {
        if (zw1D[k - KTS] >= zi) {
            kzi = k - 1;
            break;
        }
    }
}

void retrieve_exchange_coeffs(int kts, int kte, const std::vector<float>& dfm, const std::vector<float>& dfh, const std::vector<float>& dz, std::vector<float>& K_m, std::vector<float>& K_h) {
    float dzk;
    K_m[kts] = 0.0;
    K_h[kts] = 0.0;
    for (int k = kts + 1; k <= kte; ++k) {
        dzk = 0.5 * (dz[k] + dz[k - 1]);
        K_m[k] = dfm[k] * dzk;
        K_h[k] = dfh[k] * dzk;
    }
}



void mym_level2(int kts, int kte, const std::vector<float>& dz,
                const std::vector<float>& u, const std::vector<float>& v, 
                const std::vector<float>& thl, const std::vector<float>& thetav, 
                const std::vector<float>& qw, const std::vector<float>& ql, 
                const std::vector<float>& vt, const std::vector<float>& vq, 
                std::vector<float>& dtl, std::vector<float>& dqw, 
                std::vector<float>& dtv, std::vector<float>& gm, 
                std::vector<float>& gh, std::vector<float>& sm, 
                std::vector<float>& sh, float tv0, float gtr) {
    // Constants need to be defined or included from the original Fortran code
    // For example:
    // float g1, g2, b1, c1, c2, c5, a1, a2, tv0, gtr;
    // bool CKmod; // Assuming CKmod is a boolean flag

    float rfc, f1, f2, rf1, rf2, smc, shc, ri1, ri2, ri3, ri4, duz, dtz, dqz, vtt, vqq, dtq, dzk, afk, abk, ri, rf;
    float a2fac;

    rfc = g1 / (g1 + g2);
    // The rest of the calculations for f1, f2, rf1, etc. need the actual values of the constants used

    for (int k = kts + 1; k <= kte; ++k) {
        dzk = 0.5 * (dz[k] + dz[k - 1]);
        afk = dz[k] / (dz[k] + dz[k - 1]);
        abk = 1.0 - afk;
        duz = std::pow(u[k] - u[k - 1], 2) + std::pow(v[k] - v[k - 1], 2);
        duz = duz / std::pow(dzk, 2);
        dtz = (thl[k] - thl[k - 1]) / dzk;
        dqz = (qw[k] - qw[k - 1]) / dzk;

        vtt = 1.0 + vt[k] * abk + vt[k - 1] * afk; // Beta-theta
        vqq = tv0 + vq[k] * abk + vq[k - 1] * afk; // Beta-q
        dtq = vtt * dtz + vqq * dqz;

        dtl[k] = dtz;
        dqw[k] = dqz;
        dtv[k] = dtq;

        gm[k] = duz;
        gh[k] = -dtq * gtr;

        ri = -gh[k] / std::max(duz, 1.0e-10f);

        if (CKmod) {
            a2fac = 1.0 / (1.0 + std::max(ri, 0.0f));
        } else {
            a2fac = 1.0;
        }

        // Recalculate rfc, f1, f2, rf1, rf2, smc, shc with a2fac considered
        // The calculations need the actual values of the constants used

        rf = std::min(float(ri1 * (ri + ri2 - std::sqrt(std::pow(ri, 2) - ri3 * ri + ri4))), rfc);

        sh[k] = shc * (rfc - rf) / (1.0 - rf);
        sm[k] = smc * (rf1 - rf) / (rf2 - rf) * sh[k];
    }
}




