#include "ODE4.hh"
#include "CSC_RL.hh"

#include <cmath>

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

    p.omega = 2.0 * pi * 50.0;

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
    p.igd0 = 16.0;
    p.igq0 = 0.0;

    p.ed0 = Vdq_nom;
    p.eq0 = 0.0;

    p.id0 = 16.0;
    p.iq0 = 0.0;

    p.vd0 = Vdq_nom;
    p.vq0 = 0.0;

    p.istk0 = 250.0;

    p.theta_hat0 = 0.0;
    p.xi_pll0    = 0.0;

    p.xi_ucd0 = 0.0;
    p.xi_ucq0 = 0.0;
    p.xi_isd0 = 0.0;
    p.xi_isq0 = 0.0;

    p.xi_idc2_0 = 0.0;

    return p;
}

int main()
{
    try
    {
        std::system("mkdir results 2>nul");

        CSC_RL_Parameters p = make_csc_rl_parameters();

        // ------------------------------------------------------------
        // Simulation settings
        // ------------------------------------------------------------
        p.t0 = 0.0;
        p.tf = 5.0;

        const ODE::real_type h = 1e-5;
        const ODE::integer save_every = 20;

        // ------------------------------------------------------------
        // Reference signals
        // ------------------------------------------------------------
        p.idc_ref.times  = {0.0, 1.0, 2.0, 3.0, 4.0};
        p.idc_ref.values = {200.0, 250.0, 200.0, 150.0, 100.0};

        p.Q_ref.values  = {20000.0, -20000, 0};
        p.Q_ref.times = {0.0, 1.7, 3.4};

        // ------------------------------------------------------------
        // Build and solve problem
        // ------------------------------------------------------------
        CSC_RL problem(p);

        ODE::solve_rk4(problem, h, "results/csc_rl.csv", save_every);

        std::cout << "Simulation completed successfully.\n";
        std::cout << "Output file: results/csc_rl.csv\n";

        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}