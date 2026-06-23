#pragma once
#include "ODE4.hh"

#include <string>

class CSC_RL_Model;

struct CSC_RL_Parameters
{
    // Simulation time
    ODE::real_type t0 = 0.0;
    ODE::real_type tf = 0.0;

    // Grid voltage in dq frame
    ODE::real_type egd = 0.0;
    ODE::real_type egq = 0.0;

    // Grid parameters
    ODE::real_type Lg = 0.0;
    ODE::real_type Cg = 0.0;
    ODE::real_type Rg = 0.0;

    // CSC AC filter parameters
    ODE::real_type Rf = 0.0;
    ODE::real_type Lf = 0.0;
    ODE::real_type Cf = 0.0;

    // DC-side parameters
    ODE::real_type Ldc = 0.0;
    ODE::real_type RL = 0.0;

    // PLL parameters
    ODE::real_type w0_pll = 0.0;
    ODE::real_type kp_pll = 0.0;
    ODE::real_type ki_pll = 0.0;
    ODE::real_type Vdq_nom = 1.0;

    // Inner loop controller
    ODE::real_type kp1 = 0.0;
    ODE::real_type ki1 = 0.0;
    ODE::real_type kp2 = 0.0;
    ODE::real_type ki2 = 0.0;
    ODE::real_type Idc_min = 1.0;

    // Outer loop controller
    ODE::real_type kpO = 0.0;
    ODE::real_type kiO = 0.0;
    ODE::real_type V2_min = 1.0;

    // References
    ODE::StairSignal idc_ref;
    ODE::StairSignal Q_ref;

    // Initial conditions
    ODE::vec_type x0;
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

struct PLLOutput
{
    ODE::real_type e_pll = 0.0;
    ODE::real_type w_hat = 0.0;
};

struct ReferenceOutput
{
    ODE::real_type idc_ref = 0.0;
    ODE::real_type Qref = 0.0;
};

struct OuterLoopOutput
{
    ODE::real_type Eidc2 = 0.0;
    ODE::real_type Pref = 0.0;
};

struct InnerLoopOutput
{
    ODE::real_type idr = 0.0;
    ODE::real_type iqr = 0.0;
    // Stage 2
    ODE::real_type Eid = 0.0;
    ODE::real_type Eiq = 0.0;
    ODE::real_type vdr = 0.0;
    ODE::real_type vqr = 0.0;
    // Stage 1
    ODE::real_type Evd = 0.0;
    ODE::real_type Evq = 0.0;
    ODE::real_type md = 0.0;
    ODE::real_type mq = 0.0;
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
    const CSC_RL_Model& model_;
    CSC_RL_Parameters p_;
    void validate_parameters() const;

public:

    CSC_RL(const CSC_RL_Model& model, const CSC_RL_Parameters& parameters);
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
    ODE::vec_type outputs(ODE::real_type t, const ODE::vec_type& x) const override;
};
