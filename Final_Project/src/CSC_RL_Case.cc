#include "CSC_RL_Case.hh"
#include "cmath"

/**
 * @file CSC_RL_Case.cc
 * @brief Default nominal parameter set for the CSC-RL system.
 *
 * @details
 * This file defines the base quantities, grid parameters, PWM-CSC filter
 * parameters, DC-side load parameters, PLL gains, controller gains, reference
 * profiles, and initial conditions used by the default simulation case.
 *
 * The electrical base quantities are
 *
 * \f[
 * Z_b=\frac{V_b^2}{S_b},
 * \qquad
 * L_b=\frac{Z_b}{\omega_b},
 * \qquad
 * C_b=\frac{1}{Z_b\omega_b}.
 * \f]
 *
 * The DC-link inductor is selected from a current-ripple approximation:
 *
 * \f[
 * L_{dc}=\frac{V_bT_{sw}}{\Delta i}.
 * \f]
 */

// ============================================================
// Default parameter set for the CSC-RL case.
//
// This file defines the nominal electrical bases, grid parameters,
// CSC filter parameters, DC-side parameters, PLL gains, controller
// gains, initial conditions, and default reference profiles.
// ============================================================

/**
 * @brief Creates the nominal CSC-RL simulation case.
 *
 * @return Parameter set containing model parameters, controller gains,
 * references, and initial conditions.
 */

CSC_RL_Parameters make_csc_rl_case()
{
    CSC_RL_Parameters p;

    const ODE::real_type pi = std::acos(-1.0);  // Compute pi

    // ------------------------------------------------------------
    // Base parameters
    // ------------------------------------------------------------
    const ODE::real_type Vbase = 380.0;                     // Base voltage [V]
    const ODE::real_type Sbase = 60e3;                      // Base power [VA]
    const ODE::real_type wbase = 2.0 * pi * 50.0;           // Base angular frequency [rad/s]
    const ODE::real_type Zbase = Vbase * Vbase / Sbase;     // Base impedance [ohm]
    const ODE::real_type Lbase = Zbase / wbase;             // Base inductance [H]
    const ODE::real_type Cbase = 1.0 / (Zbase * wbase);     // Base capacitance [F]
    const ODE::real_type Tsw     = 1.0 / 10e3;              // Switching time [s]
    const ODE::real_type delta_i = 2.0;                     // Ripple of DC current [A]

    // ------------------------------------------------------------
    // PLL design
    // ------------------------------------------------------------
    const ODE::real_type f0_Hz = 50.0;                                        // Nominal grid frequency [Hz]
    const ODE::real_type w0_pll = 2.0 * pi * f0_Hz;                           // Nominal PLL angular frequency [rad/s]
    const ODE::real_type fn_Hz = 4.0;                       
    const ODE::real_type zeta  = 1.0 / std::sqrt(2.0);      
    const ODE::real_type wn_pll = 2.0 * pi * fn_Hz;         
    const ODE::real_type tau_pll = 1.0 / (zeta * wn_pll);   
    const ODE::real_type ts_pll  = 4.6 * tau_pll;           
    const ODE::real_type kp_pll = 9.2 / ts_pll;                               // PLL proportional gain.
    const ODE::real_type ki_pll = 21.16 / (ts_pll * ts_pll * zeta * zeta);    // PLL integral gain.
    p.w0_pll  = w0_pll;
    p.kp_pll  = kp_pll;
    p.ki_pll  = ki_pll;
    p.Vdq_nom = Vbase;

    // Default simulation time interval.
    p.t0 =  0.0;
    p.tf = 10.0;

    // ------------------------------------------------------------
    // Grid-side parameters expressed from base values.
    // ------------------------------------------------------------
    p.Lg = 0.02 * Lbase;                   // Inductor grid
    p.Cg = 0.05 * Cbase;                   // Capacitor grid
    p.Rg = 0.01 * Zbase;                   // Resistance grid
    p.egd = Vbase;                         // Grid voltage vector aligned with d-axis
    p.egq = 0.0;                           // 

    // ------------------------------------------------------------
    // CSC parameters
    // ------------------------------------------------------------
    p.Rf  = 0.01 * Zbase;                  // Resistance filter AC
    p.Lf  = 0.05 * Lbase;                  // Inductance filter AC
    p.Cf  = 0.10 * Cbase;                  // Capacitor  filter AC
    p.Ldc = (Vbase * Tsw) / delta_i;       // Inductor filter DC
    p.RL  = 0.6;                           // Resistance load

    // ------------------------------------------------------------
    // Outer loop controller parameters
    // ------------------------------------------------------------
    const ODE::real_type fno = 5.0;
    const ODE::real_type wno = 2.0 * pi * fno;
    const ODE::real_type zetao = 1.0 / std::sqrt(2.0);
    // PI gains for the outer DC current control loop.
    p.kpO = 2.0 * zetao * wno * p.Ldc;
    p.kiO = wno * wno * p.Ldc;

    p.V2_min = 1.0;             // For avoid divition by zero
   
    // ------------------------------------------------------------
    // Inner loop controller parameters
    // ------------------------------------------------------------
    // Stage-2 current-loop gains based on desired settling time.
    const ODE::real_type ts2    = 10e-3;
    const ODE::real_type alpha2 = std::log(9.0) / ts2;
    p.ki2 = alpha2 * p.Rf;
    p.kp2 = alpha2 * p.Lf;
    // Stage-1 voltage-loop gains based on desired fast settling time.
    const ODE::real_type ts1    = 0.05e-3;
    const ODE::real_type alpha1 = std::log(9.0) / ts1;
    p.ki1 = alpha1 * p.Rf;
    p.kp1 = alpha1 * p.Lf;

    p.Idc_min = 20.0;          // For avoid divition by zero

    // ------------------------------------------------------------
    // Initial conditions
    // ------------------------------------------------------------

    // Allocate and initialize the state vector.
    p.x0.resize(NSTATES);
    p.x0.setZero();

    // AC source / PCC
    p.x0(IGD) = 16.0;
    p.x0(IGQ) = 0.0;
    p.x0(ED)  = Vbase;
    p.x0(EQ)  = 0.0;

    // CSC
    p.x0(ID)   = 16.0;
    p.x0(IQ)   = 0.0;
    p.x0(VD)   = Vbase;
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

    // Default dc current and reactive power reference profile.
    p.idc_ref.times  = {0.0, 0.5, 2.0, 3.0, 4.0};
    p.idc_ref.values = {200.0, 250.0, 200.0, 150.0, 100.0};
    p.Q_ref.times  = {0.0, 1.7, 3.4};
    p.Q_ref.values = {30000.0, -30000.0, 0.0};

    return p;
}