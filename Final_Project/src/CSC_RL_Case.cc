#include "CSC_RL_Case.hh"
#include "cmath"

// ============================================================
// Create CSC-RL parameters
// ============================================================
CSC_RL_Parameters make_csc_rl_parameters()
{
    const ODE::real_type pi = std::acos(-1.0);

    CSC_RL_Parameters p;

    // ------------------------------------------------------------
    // Base parameters
    // ------------------------------------------------------------
    const ODE::real_type Vbase = 380.0;
    const ODE::real_type Sbase = 60e3;
    const ODE::real_type wbase = 2.0 * pi * 50.0;

    const ODE::real_type Zbase = Vbase * Vbase / Sbase;
    const ODE::real_type Lbase = Zbase / wbase;
    const ODE::real_type Cbase = 1.0 / (Zbase * wbase);

    const ODE::real_type Tsw     = 1.0 / 10e3;
    const ODE::real_type delta_i = 2.0;

    // ------------------------------------------------------------
    // PLL design
    // ------------------------------------------------------------
    const ODE::real_type f0_Hz = 50.0;
    const ODE::real_type w0_pll = 2.0 * pi * f0_Hz;

    const ODE::real_type Vdq_nom = Vbase;

    const ODE::real_type fn_Hz = 4.0;
    const ODE::real_type zeta  = 1.0 / std::sqrt(2.0);
    const ODE::real_type wn_pll = 2.0 * pi * fn_Hz;

    const ODE::real_type tau_pll = 1.0 / (zeta * wn_pll);
    const ODE::real_type ts_pll  = 4.6 * tau_pll;

    const ODE::real_type kp_pll = 9.2 / ts_pll;
    const ODE::real_type ki_pll = 21.16 / (ts_pll * ts_pll * zeta * zeta);

    p.w0_pll  = w0_pll;
    p.kp_pll  = kp_pll;
    p.ki_pll  = ki_pll;
    p.Vdq_nom = Vdq_nom;

    // ------------------------------------------------------------
    // Grid parameters
    // ------------------------------------------------------------
    p.Lg = 0.02 * Lbase;
    p.Cg = 0.05 * Cbase;
    p.Rg = 0.01 * Zbase;

    // ------------------------------------------------------------
    // CSC parameters
    // ------------------------------------------------------------
    p.Rf  = 0.01 * Zbase;
    p.Lf  = 0.05 * Lbase;
    p.Cf  = 0.10 * Cbase;
    p.Ldc = (Vbase * Tsw) / delta_i;
    p.RL  = 0.6;

    // ------------------------------------------------------------
    // Grid voltage in dq frame
    // ------------------------------------------------------------
    p.egd = Vbase;
    p.egq = 0.0;

    // ------------------------------------------------------------
    // Outer loop controller parameters
    // ------------------------------------------------------------
    const ODE::real_type fno = 5.0;
    const ODE::real_type wno = 2.0 * pi * fno;
    const ODE::real_type zetao = 1.0 / std::sqrt(2.0);
    p.kpO = 2.0 * zetao * wno * p.Ldc;
    p.kiO = wno * wno * p.Ldc;

    p.V2_min = 1.0;

    // ------------------------------------------------------------
    // Inner loop controller parameters
    // ------------------------------------------------------------
    const ODE::real_type ts2    = 10e-3;
    const ODE::real_type alpha2 = std::log(9.0) / ts2;
    p.ki2 = alpha2 * p.Rf;
    p.kp2 = alpha2 * p.Lf;

    const ODE::real_type ts1    = 0.05e-3;
    const ODE::real_type alpha1 = std::log(9.0) / ts1;
    p.ki1 = alpha1 * p.Rf;
    p.kp1 = alpha1 * p.Lf;

    p.Idc_min = 20.0;

    // ------------------------------------------------------------
    // Initial conditions
    // ------------------------------------------------------------
    p.x0.resize(NSTATES);
    p.x0.setZero();

    // AC source / PCC
    p.x0(IGD) = 16.0;
    p.x0(IGQ) = 0.0;
    p.x0(ED)  = Vdq_nom;
    p.x0(EQ)  = 0.0;

    // CSC
    p.x0(ID)   = 16.0;
    p.x0(IQ)   = 0.0;
    p.x0(VD)   = Vdq_nom;
    p.x0(VQ)   = 0.0;
    p.x0(ISTK) = 250.0;

    // PLL
    p.x0(THETA_HAT) = 0.0;
    p.x0(XI_PLL)    = 0.0;

    // Inner loop integrators
    p.x0(XI_UCD) = 0.0;
    p.x0(XI_UCQ) = 0.0;
    p.x0(XI_ISD) = 0.0;
    p.x0(XI_ISQ) = 0.0;

    // Outer loop integrator
    p.x0(XI_IDC2) = 0.0;

    return p;
}

void configure_references(CSC_RL_Parameters& p)
{
    p.idc_ref.times  = {0.0, 0.5, 2.0, 3.0, 4.0};
    p.idc_ref.values = {200.0, 250.0, 200.0, 150.0, 100.0};

    p.Q_ref.times  = {0.0, 1.7, 3.4};
    p.Q_ref.values = {30000.0, -30000.0, 0.0};
}