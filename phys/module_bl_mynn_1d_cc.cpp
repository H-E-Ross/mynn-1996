#include <algorithm> 
#include <iostream>
#include <vector>
#include <cmath>
#include <functional>

extern "C" void mym_predict(int kts, int kte, float closure, float delt, float* dz, float* ust, float flt, float flq, float pmz, float phh, float* el, float* dfq, float* rho, float* pdk, float* pdt, float* pdq, float* pdc, float* qke, float* tsq, float* qsq, float* cov, float* s_aw, float* s_awqke, int bl_mynn_edmf_tke, int tke_budget, float xlvcp, float xlscp, float karman);

extern "C" void mynn_mix_chem_cc(int kts, int kte, int i,float delt, std::vector<float>& dz, float pblh, int nchem, int kdvel, int ndvel,std::vector<std::vector<float>>& chem1, std::vector<float>& vd1, std::vector<float>& rho,float flt, std::vector<float>& tcd, std::vector<float>& qcd, std::vector<float>& dfh,std::vector<float>& s_aw, std::vector<std::vector<float>>& s_awchem, float emis_ant_no, float frp, int rrfs_sd, int enh_mix); 

extern "C" void mynn_tendencies(int kts, int kte, int i, float delt, float* dz, float* rho, float* u, float* v, float* th, float* tk, float* qv, float* qc, float* qi, float* qs, float* qnc, float* qni, float* psfc, float* p, float* exner, float* thl, float* sqv, float* sqc, float* sqi, float* sqs, float* sqw, float* qnwfa, float* qnifa, float* qnbca, float* ozone, float* ust, float flt, float flq, float flqv, float flqc, float wspd, float uoce, float voce, float* tsq, float* qsq, float* cov, float* tcd, float* qcd, float* dfm, float* dfh, float* dfq, float* Du, float* Dv, float* Dth, float* Dqv, float* Dqc, float* Dqi, float* Dqs, float* Dqnc, float* Dqni, float* Dqnwfa, float* Dqnifa, float* Dqnbca, float* Dozone, float* diss_heat, float* s_aw, float* s_awthl, float* s_awqt, float* s_awqv, float* s_awqc, float* s_awu, float* s_awv, float* s_awqnc, float* s_awqni, float* s_awqnwfa, float* s_awqnifa, float* s_awqnbca, float* sd_aw, float* sd_awthl, float* sd_awqt, float* sd_awqv, float* sd_awqc, float* sd_awu, float* sd_awv, float* sub_thl, float* sub_sqv, float* sub_u, float* sub_v, float* det_thl, float* det_sqv, float* det_sqc, float* det_u, float* det_v, bool FLAG_QC, bool FLAG_QI, bool FLAG_QNC, bool FLAG_QNI, bool FLAG_QS, bool FLAG_QNWFA, bool FLAG_QNIFA, bool FLAG_QNBCA, bool FLAG_OZONE, float* cldfra_bl1d, int bl_mynn_cloudmix, int bl_mynn_mixqt, int bl_mynn_edmf, int bl_mynn_edmf_mom, int bl_mynn_mixscalars, bool debug_code, float* rhoinv, float* sqw2, float R_d, float p608);

using std::bind;

//----------------------------------------CONTSTANTS-------------------------------------------

// Constants
const float NO_threshold = 10.0;     // For anthropogenic sources
const float frp_threshold = 10.0;    // Increased the frp threshold to enhance mixing over big fires
const float pblh_threshold = 100.0;

// Assuming float corresponds to float precision
constexpr float cphm_st = 5.0, cphm_unst = 16.0,
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



// Interally used


// Function to solve system of linear equations on tridiagonal matrix n times n
// after Peaceman and Rachford, 1955
// a, b, c, d - are vectors of order n
// a, b, c - are coefficients on the LHS
// d - is initially RHS on the output becomes a solution vector
void tridiag_c(int n, const std::vector<float>& a, const std::vector<float>& b, std::vector<float>& c, std::vector<float>& d) {
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

void tridiag2_c(int n, float* a, float* b, float* c, float* d, float* x) {
    float* cp = new float[n];
    float* dp = new float[n];
    float m;

    // Initialize c-prime and d-prime
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    // Solve for vectors c-prime and d-prime
    for (int i = 1; i < n; ++i) {
        m = b[i] - cp[i - 1] * a[i];
        cp[i] = c[i] / m;
        dp[i] = (d[i] - dp[i - 1] * a[i]) / m;
    }

    // Initialize x
    x[n - 1] = dp[n - 1];

    // Solve for x from the vectors c-prime and d-prime
    for (int i = n - 2; i >= 0; --i) {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }

    delete[] cp;
    delete[] dp;
}

// Function to perform tridiagonal matrix algorithm
void tridiag3_c(int kte, std::vector<float>& a, std::vector<float>& b, std::vector<float>& c, std::vector<float>& d, std::vector<float>& x) {
    // Inversion and resolution of a tridiagonal matrix A X = D
    // a - lower diagonal (Ai,i-1)
    // b - principal diagonal (Ai,i)
    // c - upper diagonal (Ai,i+1)
    // d - right-hand side vector
    // x - solution vector

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
//\param vt
//\param vq
//\param dtl     vertical gradient of \f$\theta_l\f$ (\f$K m^{-1}\f$)
//\param dqw     vertical gradient of \f$Q_w\f$
//\param dtv     vertical gradient of \f$\theta_V\f$ (\f$K m^{-1}\f$)
//\param gm      \f$G_M\f$ divided by \f$L^{2}/q^{2}\f$ (\f$s^{-2}\f$)
//\param gh      \f$G_H\f$ divided by \f$L^{2}/q^{2}\f$ (\f$s^{-2}\f$)
//\param sm      stability function for momentum, at Level 2
//\param sh      stability function for heat, at Level 2
//\section gen_mym_level2 GSD MYNN-EDMF mym_level2 General Algorithm
// @ {
#include <algorithm> // for std::max and std::min
#include <cmath>     // for std::sqrt

void mym_level2(int kts, int kte, float* dz, float* u, float* v,
                float* thl, float* thetav, float* qw, float* ql,
                float* vt, float* vq, float* dtl, float* dqw,
                float* dtv, float* gm, float* gh, float* sm, float* sh, 
		float tv0, float gtr, float sqw2) {
    float rfc, f1, f2, rf1, rf2, smc, shc, ri1, ri2, ri3, ri4, duz, dtz, dqz, vtt, vqq, dtq, dzk, afk, abk, ri, rf;
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

        vtt = 1.0 + vt[k] * abk + vt[k - 1] * afk; // Beta-theta in NN09, Eq. 39
        vqq = tv0 + vq[k] * abk + vq[k - 1] * afk; // Beta-q
        dtq = vtt * dtz + vqq * dqz;
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
void mym_length(int kts, int kte, float xland, float* dz, float* zw, float rmo, float flt, float fltv, float flq, float* vt, float* vq, float* u1, float* v1, float* qke, float* dtv, float* el, float zi, float* theta, float Psig_bl, float* cldfra_bl1D, int bl_mynn_mixlength, float* edmf_w1, float* edmf_a1) {
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
            zi = std::min(10000.0, double(zw[kte-2]));
            h1 = std::max(0.3 * double(zi), double(mindz));
            h1 = std::min(double(h1), double(maxdz));
            h2 = h1 / 2.0;
            qkw[kts] = std::sqrt(std::max(double(qke[kts]), 1.0e-10));
            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0 - afk;
                qkw[k] = std::sqrt(std::max(double(qke[k] * abk + qke[k-1] * afk), 1.0e-3));
            }
            k = kts + 1;
            zwk = zw[k];
            while (zwk <= zi + h1) {
                dzk = 0.5 * (dz[k] + dz[k-1]);
                qdz = std::max(double(qkw[k] - qmin), 0.03) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
            }
            elt = alp1 * elt / vsc;
            vflx = (vt[kts] + 1.0) * flt + (vq[kts] + tv0) * flq;
            vsc = std::pow(gtr * elt * std::max(double(vflx), 0.0), onethird);
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
            wt_u = (1.0 - std::min(std::max(double(Ugrid - Uonset), 0.0) / 30.0, 0.5));
            cns = 2.7;
            alp1 = 0.23;
            alp2 = 0.3;
            alp3 = 2.5 * wt_u;
            alp4 = 5.0;
            alp5 = 0.3;
            alp6 = 50.0;
            zi = std::max(double(zi), double(minzi));
            h1 = std::max(double(0.3 * zi), 300.0);
            h1 = std::min(double(h1), 600.0);
            h2 = h1 / 2.0;
            qtke[kts] = std::max(double(0.5 * qke[kts]), 0.01);
            thetaw[kts] = theta[kts];
            qkw[kts] = std::sqrt(std::max(double(qke[kts]), 1.0e-10));
            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0 - afk;
                qkw[k] = std::sqrt(std::max(double(qke[k] * abk + qke[k-1] * afk), 1.0e-3));
                qtke[k] = 0.5 * qkw[k] * qkw[k];
                thetaw[k] = theta[k] * abk + theta[k-1] * afk;
            }
            k = kts + 1;
            zwk = zw[k];
            while (zwk <= zi + h1) {
                dzk = 0.5 * (dz[k] + dz[k-1]);
                qdz = std::min(std::max(double(qkw[k] - qmin), 0.03), 30.0) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
            }
            elt = std::min(std::max(double(alp1 * elt / vsc), 10.0), 400.0);
            vflx = fltv;
            vsc = std::pow(gtr * elt * std::max(double(vflx), 0.0), onethird);
            el[kts] = 0.0;
            zwk1 = zw[kts+1];
            for (k = kts+1; k <= kte; k++) {
                zwk = zw[k];
                if (dtv[k] > 0.0) {
                    bv = std::max(double(std::sqrt(gtr * dtv[k])), 0.0001);
                    elb = std::max(double(alp2 * qkw[k]), double(alp6 * edmf_a1[k-1] * edmf_w1[k-1])) / bv * (1.0 + alp3 * std::sqrt(vsc / (bv * elt)));
                    elb = std::min(elb, zwk);
                    elf = 1.0 * qkw[k] / bv;
                    elBLavg[k] = std::max(double(elBLavg[k]), double(alp6 * edmf_a1[k-1] * edmf_w1[k-1] / bv));
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
            zi = std::max(double(zi), double(minzi));
            h1 = std::max(double(0.3 * zi), 300.0);
            h1 = std::min(double(h1), 600.0);
            h2 = h1 * 0.5;
            qtke[kts] = std::max(double(0.5 * qke[kts]), 0.01);
            qkw[kts] = std::sqrt(std::max(double(qke[kts]), 1.0e-4));
            for (k = kts+1; k <= kte; k++) {
                afk = dz[k] / (dz[k] + dz[k-1]);
                abk = 1.0 - afk;
                qkw[k] = std::sqrt(std::max(double(qke[k] * abk + qke[k-1] * afk), 1.0e-3));
                qtke[k] = 0.5 * qkw[k] * qkw[k];
            }
            k = kts + 1;
            zwk = zw[k];
            while (zwk <= zi + h1) {
                dzk = 0.5 * (dz[k] + dz[k-1]);
                qdz = std::min(std::max(double(qkw[k] - qmin), 0.03), 30.0) * dzk;
                elt = elt + qdz * zwk;
                vsc = vsc + qdz;
                k = k + 1;
                zwk = zw[k];
            }
            elt = std::min(std::max(double(alp1 * elt / vsc), 10.0), 400.0);
            vflx = fltv;
            vsc = std::pow(gtr * elt * std::max(double(vflx), 0.0), onethird);
            el[kts] = 0.0;
            zwk1 = zw[kts+1];
            for (k = kts+1; k <= kte; k++) {
                zwk = zw[k];
                dzk = 0.5 * (dz[k] + dz[k-1]);
                cldavg = 0.5 * (cldfra_bl1D[k-1] + cldfra_bl1D[k]);
                if (dtv[k] > 0.0) {
                    bv = std::max(double(std::sqrt(gtr * dtv[k])), 0.001);
                    elb_mf = std::max(double(alp2 * qkw[k]), double(alp6 * edmf_a1[k-1] * edmf_w1[k-1]) / bv * (1.0 + alp3 * std::sqrt(vsc / (bv * elt))));
                    elb = std::min(std::max(double(alp5 * qkw[k]), double(alp6 * edmf_a1[k] * edmf_w1[k]) / bv), double(zwk));
                    wstar = 1.25 * std::pow(gtr * zi * std::max(double(vflx), 1.0e-4), onethird);
                    tau_cloud = std::min(std::max(double(ctau * wstar / grav), 30.0), 150.0);
                    float wt = 0.5 * std::tanh((zwk - (zi + h1)) / h2) + 0.5;
                    tau_cloud = tau_cloud * (1.0 - wt) + 50.0 * wt;
                    elf = std::min(std::max(double(tau_cloud * std::sqrt(std::min(double(qtke[k]), 40.0))), double(alp6 * edmf_a1[k] * edmf_w1[k] / bv)), double(zwk));
                } else {
                    wstar = 1.25 * std::pow(gtr * zi * std::max(double(vflx), 1.0e-4), onethird);
                    tau_cloud = std::min(std::max(double(ctau * wstar / grav), 50.0), 200.0);
                    float wt = 0.5 * std::tanh((zwk - (zi + h1)) / h2) + 0.5;
                    tau_cloud = tau_cloud * (1.0 - wt) + std::max(100.0f, dzk * 0.25f) * wt;
                    elb = std::min(tau_cloud * std::sqrt(std::min(qtke[k], 40.0f)), zwk);
                    elf = elb;
                    elb_mf = elb;
                }
                elf = elf / (1.0 + (elf / 800.0));
                elb_mf = std::max(double(elb_mf), 0.01);
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
void moisture_check_c(int kte, float delt, float* dp, float* exner,
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
        qv[k-1] = std::max(double(qv[k-1]), double(qvmin));
        qc[k-1] = std::max(double(qc[k-1]), double(qcmin));
        qi[k-1] = std::max(double(qi[k-1]), double(qimin));
        qs[k-1] = std::max(double(qs[k-1]), double(qimin));
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
void mym_predict_c(int kts, int kte, float closure, float delt, float* dz, float* ust, float flt, float flq, float pmz, float phh, float* el, float* dfq, float* rho, float* pdk, float* pdt, float* pdq, float* pdc, float* qke, float* tsq, float* qsq, float* cov, float* s_aw, float* s_awqke, int bl_mynn_edmf_tke, int tke_budget, float xlvcp, float xlscp, float karman) {
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
    tridiag2_c(kte, a, b, c, d, x);
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
        tridiag2_c(kte, a, b, c, d, x);
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
        tridiag2_c(kte, a, b, c, d, x);
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
        tridiag2_c(kte, a, b, c, d, x);
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
        rhoz[k - kts] = std::max(double(rhoz[k - kts]), 1E-4);
        rhoinv[k - kts] = 1.0 / std::max(double(rho[k - 1]), 1E-4);
        float dzk = 0.5 * (dz[k - 1] + dz[k - 2]);
        khdz[k - kts] = rhoz[k - kts] * dfh[k - 1];
    }
    rhoz[kte - kts + 1] = rhoz[kte - kts];
    khdz[kte - kts + 1] = rhoz[kte - kts + 1] * dfh[kte - 1];
    // Stability criteria for mf
    for (int k = kts + 1; k <= kte - 1; ++k) {
        khdz[k - kts] = std::max(double(khdz[k - kts]), double(0.5 * s_aw[k - kts]));
        khdz[k - kts] = std::max(double(khdz[k - kts]), double(-0.5 * (s_aw[k - kts] - s_aw[k - kts + 1])));
    }
    // Enhanced mixing over fires
    if (rrfs_sd==1 && enh_mix==1) {
        for (int k = kts + 1; k <= kte - 1; ++k) {
            float khdz_old = khdz[k - kts];
            float khdz_back = pblh * 0.15 / dz[k - 1];
            // Modify based on anthropogenic emissions of NO and FRP
            if (pblh < pblh_threshold) {
                if (emis_ant_no > NO_threshold) {
                    khdz[k - kts] = std::max(1.1 * double(khdz[k - kts]), double(std::sqrt((emis_ant_no / NO_threshold)) / dz[k - 1] * rhoz[k - kts]));
                }
                if (frp > frp_threshold) {
                    int kmaxfire = std::ceil(std::log(frp));
                    khdz[k - kts] = std::max(double(1.1 * khdz[k - kts]), double((1.0 - k / (kmaxfire * 2.0)) * (std::pow(std::log(frp), 2.0) - 2.0 * std::log(frp)) / dz[k - 1] * rhoz[k - kts]));
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
        tridiag3_c(kte, a, b, c, d, x);
        for (k = kts; k <= kte; ++k) {
            chem1[k - 1][ic] = x[k - kts];
        }
    }
}


// ==================================================================
//>\ingroup gsd_mynn_edmf
// This subroutine solves for tendencies of U, V, \f$\theta\f$, qv,
// qc, and qi
void mynn_tendencies(int kts, int kte, int i, float delt, float* dz, float* rho, float* u, float* v, float* th, float* tk, float* qv, float* qc, float* qi, float* qs, float* qnc, float* qni, float* psfc, float* p, float* exner, float* thl, float* sqv, float* sqc, float* sqi, float* sqs, float* sqw, float* qnwfa, float* qnifa, float* qnbca, float* ozone, float* ust, float flt, float flq, float flqv, float flqc, float wspd, float uoce, float voce, float* tsq, float* qsq, float* cov, float* tcd, float* qcd, float* dfm, float* dfh, float* dfq, float* Du, float* Dv, float* Dth, float* Dqv, float* Dqc, float* Dqi, float* Dqs, float* Dqnc, float* Dqni, float* Dqnwfa, float* Dqnifa, float* Dqnbca, float* Dozone, float* diss_heat, float* s_aw, float* s_awthl, float* s_awqt, float* s_awqv, float* s_awqc, float* s_awu, float* s_awv, float* s_awqnc, float* s_awqni, float* s_awqnwfa, float* s_awqnifa, float* s_awqnbca, float* sd_aw, float* sd_awthl, float* sd_awqt, float* sd_awqv, float* sd_awqc, float* sd_awu, float* sd_awv, float* sub_thl, float* sub_sqv, float* sub_u, float* sub_v, float* det_thl, float* det_sqv, float* det_sqc, float* det_u, float* det_v, bool FLAG_QC, bool FLAG_QI, bool FLAG_QNC, bool FLAG_QNI, bool FLAG_QS, bool FLAG_QNWFA, bool FLAG_QNIFA, bool FLAG_QNBCA, bool FLAG_OZONE, float* cldfra_bl1d, int bl_mynn_cloudmix, int bl_mynn_mixqt, int bl_mynn_edmf, int bl_mynn_edmf_mom, int bl_mynn_mixscalars, bool debug_code, float* rhoinv, float* sqw2, float R_d, float p608) {
    float nonloc = 1.0;
    float dztop = 0.5 * (dz[kte] + dz[kte-1]);
    float onoff = (bl_mynn_edmf_mom == 0) ? 0.0 : 1.0;
    float rhosfc = *psfc / (R_d * (tk[kts] + p608 * qv[kts]));
    float dtz[kte+1], dfhc[kte], dfmc[kte], delp[kte], sqv2[kte], sqc2[kte];
    float rhoz[kte+1], khdz[kte+1], kmdz[kte+1];
    float a[kte+1], b[kte+1], c[kte+1], d[kte+1], x[kte+1];
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
        rhoz[k] = std::max(double(rhoz[k]), 1E-4);
        rhoinv[k] = 1.0 / std::max(double(rho[k]), 1E-4);
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
        khdz[k] = std::max(double(khdz[k]), double(0.5 * s_aw[k]));
        khdz[k] = std::max(double(khdz[k]), double(-0.5 * (s_aw[k] - s_aw[k+1])));
        kmdz[k] = std::max(double(kmdz[k]), double(0.5 * s_aw[k]));
        kmdz[k] = std::max(double(kmdz[k]), double(-0.5 * (s_aw[k] - s_aw[k+1])));
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
    tridiag2_c(kte, a, b, c, d, x);
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
    tridiag2_c(kte, a, b, c, d, x);
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
    tridiag2_c(kte, a, b, c, d, x);
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
        tridiag2_c(kte, a, b, c, d, sqw2);
    } else {
        for (int k = kts; k <= kte; k++) {
            sqw2[k] = sqw[k];
        }
    }
    
    if (bl_mynn_mixqt == 0) {
        if (bl_mynn_cloudmix > 0 && FLAG_QC) {
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
            tridiag2_c(kte, a, b, c, d, sqc2);
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
            qvflux = std::max(double(qvflux), (std::min(0.9 * double(sqv[kts]) - 1e-8, 0.0) / dtz[kts]));
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
        tridiag2_c(kte, a, b, c, d, sqv2);
    } else {
        for (int k = kts; k <= kte; k++) {
            sqv2[k] = sqv[k];
        }
    }
}



