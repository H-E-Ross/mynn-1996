#include <algorithm> 
#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <limits>

extern "C" void mym_predict_cc(int kts, int kte, 
                 double closure, 
                 double delt,
                 double dz[],
                 double ust, double flt, double flq, double pmz, double phh,
                 double el[], double dfq[], double rho[],
                 double pdk[], double pdt[], double pdq[], double pdc[],
                 double qke[], double tsq[], double qsq[], double cov[],
                 double s_aw1[], double s_awqke1[], int bl_mynn_edmf_tke,
                 double qWT1[], double qDISS1[], int tke_budget, double xlvcp, 
		 double xlscp, double karman, double Sqfac, double b1, double b2,
                 double qkemin);


// ===============================INTERNALLY USED ===============================

void tridiag2(int kte, const double* a, const double* b, const double* c, const double* d, double* x) {
    // a - sub-diagonal (means it is the diagonal below the main diagonal)
    // b - the main diagonal
    // c - sup-diagonal (means it is the diagonal above the main diagonal)
    // d - right part
    // x - the answer
    // kte - number of unknowns (levels)

    std::vector<double> cp(kte);
    std::vector<double> dp(kte);
    double m;

    // initialize c-prime and d-prime
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    // solve for vectors c-prime and d-prime
    for (int k = 1; k < kte; ++k) {
        m = b[k] - cp[k - 1] * a[k];
        cp[k] = c[k] / m;
        dp[k] = (d[k] - dp[k - 1] * a[k]) / m;
    }

    // initialize x
    x[kte - 1] = dp[kte - 1];

    // solve for x from the vectors c-prime and d-prime
    for (int k = kte - 2; k >= 0; --k) {
        x[k] = dp[k] - cp[k] * x[k + 1];
    }
}

// ==============================================================================
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
void mym_predict_cc(int kts, int kte,
                 double closure,
                 double delt,
                 double dz[],
                 double ust, double flt, double flq, double pmz, double phh,
                 double el[], double dfq[], double rho[],
                 double pdk[], double pdt[], double pdq[], double pdc[],
                 double qke[], double tsq[], double qsq[], double cov[],
                 double s_aw1[], double s_awqke1[], int bl_mynn_edmf_tke,
                 double qWT1[], double qDISS1[], int tke_budget, double xlvcp, 
		 double xlscp, double karman,double Sqfac, double b1, double b2,
		 double qkemin)
{
    // Constants
    const double zero = 0.0;
    const double one = 1.0;

    // Variables
    int k;
    double vkz, pdk1, phm, pdt1, pdq1, pdc1, b1l, b2l, onoff;
    std::vector<double> dtz(kte - kts + 1), a(kte - kts + 1), b(kte - kts + 1), c(kte - kts + 1), d(kte - kts + 1), x(kte - kts + 1);
    std::vector<double> qkw(kte - kts + 1), bp(kte - kts + 1), rp(kte - kts + 1), df3q(kte - kts + 1);
    std::vector<double> rhoinv(kte - kts + 1), rhoz(kte - kts + 2), kqdz(kte - kts + 1), kmdz(kte - kts + 1);

    // REGULATE THE MOMENTUM MIXING FROM THE MASS-FLUX SCHEME (on or off)
    onoff = (bl_mynn_edmf_tke == 0) ? zero : one;

    vkz = karman * 0.5 * dz[kts];

    for (k = kts; k <= kte; ++k) {
        qkw[k] = std::sqrt(std::max(qke[k], zero));
        df3q[k] = Sqfac * dfq[k];
        dtz[k] = delt / dz[k];
    }

    // Prepare "constants" for diffusion equation.
    rhoz[kts] = rho[kts];
    rhoinv[kts] = 1.0 / rho[kts];
    kqdz[kts] = rhoz[kts] * df3q[kts];
    kmdz[kts] = rhoz[kts] * dfq[kts];

    for (k = kts + 1; k <= kte; ++k) {
        rhoz[k] = (rho[k] * dz[k - 1] + rho[k - 1] * dz[k]) / (dz[k - 1] + dz[k]);
        rhoz[k] = std::max(rhoz[k], 1E-4);
        rhoinv[k] = 1.0 / std::max(rho[k], 1E-4);
        kqdz[k] = rhoz[k] * df3q[k]; // for TKE
        kmdz[k] = rhoz[k] * dfq[k];  // for T'2, q'2, and T'q'
    }
    rhoz[kte + 1] = rhoz[kte];
    kqdz[kte + 1] = rhoz[kte + 1] * df3q[kte];
    kmdz[kte + 1] = rhoz[kte + 1] * dfq[kte];

    // Stability criteria for mf
    for (k = kts + 1; k <= kte - 1; ++k) {
        kqdz[k] = std::max(kqdz[k], 0.5 * s_aw1[k]);
        kqdz[k] = std::max(kqdz[k], -0.5 * (s_aw1[k] - s_aw1[k + 1]));
        kmdz[k] = std::max(kmdz[k], 0.5 * s_aw1[k]);
        kmdz[k] = std::max(kmdz[k], -0.5 * (s_aw1[k] - s_aw1[k + 1]));
    }

    pdk1 = 2.0 * ust * ust * ust * pmz / vkz;
    phm = 2.0 / ust * phh / vkz;
    pdt1 = phm * flt * flt;
    pdq1 = phm * flq * flq;
    pdc1 = phm * flt * flq;

    pdk[kts] = pdk1 - pdk[kts + 1];
    pdt[kts] = pdt[kts + 1];
    pdq[kts] = pdq[kts + 1];
    pdc[kts] = pdc[kts + 1];

    for (k = kts; k <= kte - 1; ++k) {
        b1l = b1 * 0.5 * (el[k + 1] + el[k]);
        bp[k] = 2.0 * qkw[k] / b1l;
        rp[k] = pdk[k + 1] + pdk[k];
    }

    for (k = kts; k <= kte - 1; ++k) {
        a[k] = -dtz[k] * kqdz[k] * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * s_aw1[k] * onoff;
        b[k] = 1.0 + dtz[k] * (kqdz[k] + kqdz[k + 1]) * rhoinv[k] + 0.5 * dtz[k] * rhoinv[k] * (s_aw1[k] - s_aw1[k + 1]) * onoff + bp[k] * delt;
        c[k] = -dtz[k] * kqdz[k + 1] * rhoinv[k] - 0.5 * dtz[k] * rhoinv[k] * s_aw1[k + 1] * onoff;
        d[k] = rp[k] * delt + qke[k] + dtz[k] * rhoinv[k] * (s_awqke1[k] - s_awqke1[k + 1]) * onoff;
    }

    a[kte] = 0.0;
    b[kte] = 1.0;
    c[kte] = 0.0;
    d[kte] = qke[kte];
    x[kte] = 0.0;

    // Call to tridiagonal solver
    tridiag2(kte, a.data(), b.data(), c.data(), d.data(), x.data());

    for (k = kts; k <= kte; ++k) {
        qke[k] = std::max(x[k], qkemin);
        qke[k] = std::min(qke[k], 150.0);
    }

    // TKE budget
    if (tke_budget == 1) {
        // TKE Vertical transport
        std::vector<double> tke_up(kte - kts + 1), dzinv(kte - kts + 1);
        dzinv[k] = 1.0 / dz[k];
        k = kts;
        qWT1[k] = dzinv[k] * ((kqdz[k + 1] * (tke_up[k + 1] - tke_up[k]) - kqdz[k] * tke_up[k]) +
            0.5 * rhoinv[k] * (s_aw1[k + 1] * tke_up[k + 1] + (s_aw1[k + 1] - s_aw1[k]) * tke_up[k] +
            (s_awqke1[k] - s_awqke1[k + 1])) * onoff);

        for (k = kts + 1; k <= kte - 1; ++k) {
            qWT1[k] = dzinv[k] * ((kqdz[k + 1] * (tke_up[k + 1] - tke_up[k]) - kqdz[k] * (tke_up[k] - tke_up[k - 1])) +
                0.5 * rhoinv[k] * (s_aw1[k + 1] * tke_up[k + 1] + (s_aw1[k + 1] - s_aw1[k]) * tke_up[k] -
                s_aw1[k] * tke_up[k - 1] + (s_awqke1[k] - s_awqke1[k + 1])) * onoff);
        }
        k = kte;
        qWT1[k] = dzinv[k] * (-kqdz[k] * (tke_up[k] - tke_up[k - 1]) +
            0.5 * rhoinv[k] * (-s_aw1[k] * tke_up[k] - s_aw1[k] * tke_up[k - 1] + s_awqke1[k]) * onoff);
        for (int i = kts + 1; i <= kte - 1; ++i) {
            qDISS1[i] = bp[i] * tke_up[i]; // TKE dissipation rate
	}	
    }

    // Prediction of the moisture variance
    if (closure > 2.5) {
        for (k = kts; k <= kte - 1; ++k) {
            b2l = b2 * 0.5 * (el[k + 1] + el[k]);
            bp[k] = 2.0 * qkw[k] / b2l;
            rp[k] = pdq[k + 1] + pdq[k];
        }

        for (k = kts; k <= kte - 1; ++k) {
            a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
            b[k] = 1.0 + dtz[k] * (kmdz[k] + kmdz[k + 1]) * rhoinv[k] + bp[k] * delt;
            c[k] = -dtz[k] * kmdz[k + 1] * rhoinv[k];
            d[k] = rp[k] * delt + qsq[k];
        }

        a[kte] = -1.0;
        b[kte] = 1.0;
        c[kte] = 0.0;
        d[kte] = 0.0;

        tridiag2(kte, a.data(), b.data(), c.data(), d.data(), x.data());

        for (k = kts; k <= kte; ++k) {
            qsq[k] = std::max(x[k], 1e-17);
        }
    } else {
        for (k = kts; k <= kte - 1; ++k) {
            if (qkw[k] <= zero) {
                b2l = zero;
            } else {
                b2l = b2 * 0.25 * (el[k + 1] + el[k]) / qkw[k];
            }
            qsq[k] = b2l * (pdq[k + 1] + pdq[k]);
        }
        qsq[kte] = qsq[kte - 1];
    }

    // Prediction of the temperature variance
    if (closure >= 3.0) {
        for (k = kts; k <= kte - 1; ++k) {
            b2l = b2 * 0.5 * (el[k + 1] + el[k]);
            bp[k] = 2.0 * qkw[k] / b2l;
            rp[k] = pdt[k + 1] + pdt[k];
        }

        for (k = kts; k <= kte - 1; ++k) {
            a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
            b[k] = 1.0 + dtz[k] * (kmdz[k] + kmdz[k + 1]) * rhoinv[k] + bp[k] * delt;
            c[k] = -dtz[k] * kmdz[k + 1] * rhoinv[k];
            d[k] = rp[k] * delt + tsq[k];
        }

        a[kte] = -1.0;
        b[kte] = 1.0;
        c[kte] = 0.0;
        d[kte] = 0.0;

        tridiag2(kte, a.data(), b.data(), c.data(), d.data(), x.data());

        for (k = kts; k <= kte; ++k) {
            tsq[k] = x[k];
        }

        // Prediction of the temperature-moisture covariance
        for (k = kts; k <= kte - 1; ++k) {
            b2l = b2 * 0.5 * (el[k + 1] + el[k]);
            bp[k] = 2.0 * qkw[k] / b2l;
            rp[k] = pdc[k + 1] + pdc[k];
        }

        for (k = kts; k <= kte - 1; ++k) {
            a[k] = -dtz[k] * kmdz[k] * rhoinv[k];
            b[k] = 1.0 + dtz[k] * (kmdz[k] + kmdz[k + 1]) * rhoinv[k] + bp[k] * delt;
            c[k] = -dtz[k] * kmdz[k + 1] * rhoinv[k];
            d[k] = rp[k] * delt + cov[k];
        }

        a[kte] = -1.0;
        b[kte] = 1.0;
        c[kte] = 0.0;
        d[kte] = 0.0;

        tridiag2(kte, a.data(), b.data(), c.data(), d.data(), x.data());

        for (k = kts; k <= kte; ++k) {
            cov[k] = x[k];
        }
    } else {
        for (k = kts; k <= kte - 1; ++k) {
            if (qkw[k] <= zero) {
                b2l = zero;
            } else {
                b2l = b2 * 0.25 * (el[k + 1] + el[k]) / qkw[k];
            }
            tsq[k] = b2l * (pdt[k + 1] + pdt[k]);
            cov[k] = b2l * (pdc[k + 1] + pdc[k]);
        }
        tsq[kte] = tsq[kte - 1];
        cov[kte] = cov[kte - 1];
    }
}


