#pragma once
#include "ODE4.hh"

// ============================================================
// Generic stair signal
//
//
// Example:
//
// times  = {0.0, 1.0, 2.0, 3.0}
// values = {16.0, 50.0, 100.0, 80.0}
//
// Then:
//
// 0.0 <= t < 1.0  -> 16.0
// 1.0 <= t < 2.0  -> 50.0
// 2.0 <= t < 3.0  -> 100.0
// t >= 3.0        -> 80.0
//
// If times and values are empty, the signal uses constant_value.
// ============================================================
struct StairSignal
{

    ODE::real_type constant_value = 0.0;    // Default value if no stair is provided
    std::vector<ODE::real_type> times;      // Switching times of the stair signal
    std::vector<ODE::real_type> values;     // Values of the stair signal at the corresponding times

    ODE::real_type value(ODE::real_type t) const;
    void validate(const std::string& name) const;
};

struct CSC_RL_Parameters
{
    // ------------------------------------------------------------
    // PLL parameters
    // ------------------------------------------------------------
    ODE::real_type w0_pll = 0.0;      // Nominal angular frequency [rad/s]
    ODE::real_type kp_pll = 0.0;      // PLL proportional gain
    ODE::real_type ki_pll = 0.0;      // PLL integral gain
    ODE::real_type Vdq_nom = 1.0;      // Phase peak voltage used for normalization

    ODE::real_type theta_hat0; // Initial angle of the PLL's internal oscillator
    ODE::real_type xi_pll0;    // Initial frequency deviation of the PLL's internal oscillator

    // ------------------------------------------------------------
    // CSC AC filter parameters
    // ------------------------------------------------------------
    ODE::real_type Rf = 0.0;
    ODE::real_type Lf = 0.0;
    ODE::real_type Cf = 0.0;

    // ------------------------------------------------------------
    // DC-side parameters
    // ------------------------------------------------------------
    ODE::real_type Ldc;
    ODE::real_type RL;  // Load resistance on DC side

    // Grid parameters
    ODE::real_type Lg;
    ODE::real_type Cg;
    ODE::real_type Rg;
    ODE::real_type egd;
    ODE::real_type egq;

    // ------------------------------------------------------------
    // Initial conditions
    // x = [igd, igq, ed, eq, id, iq, vd, vq, istk theta_hat xi_pll]^T
    // ------------------------------------------------------------

    // PCC states
    ODE::real_type igd0;
    ODE::real_type igq0;
    ODE::real_type ed0;
    ODE::real_type eq0;

    // CSC states
    ODE::real_type id0;
    ODE::real_type iq0;
    ODE::real_type vd0;
    ODE::real_type vq0;
    ODE::real_type istk0;

    // ------------------------------------------------------------
    // Inner loop controller parameters
    //
    // Stage 2: current error -> capacitor voltage reference
    // Stage 1: capacitor voltage error -> modulation index
    // ------------------------------------------------------------
    ODE::real_type kp1;
    ODE::real_type ki1;

    ODE::real_type kp2;
    ODE::real_type ki2;

    ODE::real_type Idc_min = 1.0;

    // ------------------------------------------------------------
    // Inner loop initial conditions
    // x_inner = [xi_ucd, xi_ucq, xi_isd, xi_isq]^T
    // ------------------------------------------------------------
    ODE::real_type xi_ucd0;
    ODE::real_type xi_ucq0;
    ODE::real_type xi_isd0;
    ODE::real_type xi_isq0;

    // ------------------------------------------------------------
    // Outer loop initial conditions
    // ------------------------------------------------------------
    ODE::real_type kpO;      // Outer-loop proportional gain
    ODE::real_type kiO;      // Outer-loop integral gain
    ODE::real_type V2_min = 1.0;

    StairSignal idc_ref;     // DC current reference [A]
    StairSignal Q_ref;       // Reactive power reference [var]

    ODE::real_type xi_idc2_0;

    // ------------------------------------------------------------
    // Simulation time
    // ------------------------------------------------------------
    ODE::real_type t0;
    ODE::real_type tf;
};

enum StateIndex
{
    IGD = 0, IGQ, ED, EQ,               // AC Source
    ID, IQ, VD, VQ, ISTK,               // CSC
    THETA_HAT, XI_PLL,                  // PLL  
    XI_UCD, XI_UCQ, XI_ISD, XI_ISQ,     // Inner
    XI_IDC2,                            // Outer
    NSTATES
};

enum OutputIndex
{
    IA = 0, IB, IC, VA, VB, VC,       // Voltage and current in abc
    ID_R, IQ_R, IDC_REF, Q_REF, PREF, // References
    P_OUT, Q_OUT,                     // Power
    NOUTPUTS
};

struct ControlOutput
{
    ODE::real_type e_pll = 0.0;
    ODE::real_type w_hat = 0.0;

    ODE::real_type idc_ref = 0.0;
    ODE::real_type Qref = 0.0;
    ODE::real_type Pref = 0.0;

    ODE::real_type idr = 0.0;
    ODE::real_type iqr = 0.0;

    ODE::real_type Eid = 0.0;
    ODE::real_type Eiq = 0.0;
    ODE::real_type Eidc2 = 0.0;

    ODE::real_type vdr = 0.0;
    ODE::real_type vqr = 0.0;

    ODE::real_type Evd = 0.0;
    ODE::real_type Evq = 0.0;

    ODE::real_type md = 0.0;
    ODE::real_type mq = 0.0;
};

class CSC_RL : public ODE::ODE_Problem_base
{
private:
    CSC_RL_Parameters p_;
    void validate_parameters() const;
    ControlOutput compute_control(ODE::real_type t, const ODE::vec_type& x) const;

public:
    explicit CSC_RL(const CSC_RL_Parameters& parameters);

    ODE::integer n() const override;
    ODE::real_type t0() const override;
    ODE::real_type tf() const override;
    ODE::real_type initial_condition(ODE::integer i) const override;
    std::string state_name(ODE::integer i) const override;

    void rhs(ODE::real_type t, const ODE::vec_type& x, ODE::vec_type& dxdt) const override;

    ODE::integer n_outputs() const override; // Number of algebraic outputs
    std::string output_name(ODE::integer i) const override; // Name of algebraic output i
    ODE::real_type output(ODE::real_type t, const ODE::vec_type& x, ODE::integer i) const override; // value of each algebraic output
};
