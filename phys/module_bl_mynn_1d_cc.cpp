#include <algorithm> 
#include <functional>

extern "C" void mym_predict(int kts, int kte, float closure, float delt, float* dz, float* ust, float flt, float flq, float pmz, float phh, float* el, float* dfq, float* rho, float* pdk, float* pdt, float* pdq, float* pdc, float* qke, float* tsq, float* qsq, float* cov, float* s_aw, float* s_awqke, int bl_mynn_edmf_tke, int tke_budget, float xlvcp, float xlscp, float karman);


using std::bind;

//----------------------------------------CONTSTANTS-------------------------------------------
#include <cmath> // For pow function, if needed

// Assuming kind_phys corresponds to float precision
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



// Interally used

#include <iostream>

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
        qv[k-1] = std::max(qv[k-1], qvmin);
        qc[k-1] = std::max(qc[k-1], qcmin);
        qi[k-1] = std::max(qi[k-1], qimin);
        qs[k-1] = std::max(qs[k-1], qimin);
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
void mym_predict(int kts, int kte, float closure, float delt, float* dz, float* ust, float flt, float flq, float pmz, float phh, float* el, float* dfq, float* rho, float* pdk, float* pdt, float* pdq, float* pdc, float* qke, float* tsq, float* qsq, float* cov, float* s_aw, float* s_awqke, int bl_mynn_edmf_tke, int tke_budget, float xlvcp, float xlscp, float karman) {
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


