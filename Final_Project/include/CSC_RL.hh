#pragma once

/**
 * @file CSC_RL.hh
 * @brief ODE wrapper and public data structures for the CSC-RL model.
 *
 * @details
 * This file declares the parameter container, state indexes, algebraic output
 * indexes, controller-output structures, and the CSC_RL class.
 *
 * The CSC_RL class adapts the physical model implemented in CSC_RL_Model.hh
 * to the generic ODE::ODE_Problem_base interface required by the RK4 solver.
 */

#include "ODE4.hh"

#include <string>

// ============================================================
// CSC-RL ODE problem interface.
//
// The actual model equations are implemented in CSC_RL_Model.hh.
// This file only declares the problem wrapper used by the ODE solver.
// ============================================================

class CSC_RL_Model;                 // Forward declaration of the model class to avoid including CSC_RL_Model.hh here.


/**
 * @brief Complete parameter set for the CSC-RL simulation.
 *
 * @details
 * This structure stores all quantities required by the model:
 *
 * - simulation time interval,
 * - PLL gains,
 * - grid-side parameters,
 * - PWM-CSC filter parameters,
 * - DC-side load parameters,
 * - controller gains,
 * - external stair references,
 * - and initial conditions.
 *
 * The state vector stored in x0 must have NSTATES entries.
 */

// All physical parameters, controller gains,
// references, and initial conditions of the CSC-RL system.
struct CSC_RL_Parameters
{
    // Simulation time
    ODE::real_type t0 = 0.0;        // Initial simulation time [s]
    ODE::real_type tf = 0.0;        // Final simulation time [s]

    // PLL parameters
    ODE::real_type w0_pll = 0.0;    // Nominal angular frequency [rad/s]
    ODE::real_type kp_pll = 0.0;    // PLL proportional gain
    ODE::real_type ki_pll = 0.0;    // PLL integral gain
    ODE::real_type Vdq_nom = 1.0;   // Nominal dq voltage used for PLL normalization [V]

    // Grid parameters
    ODE::real_type Lg = 0.0;        // Grid inductance [H]
    ODE::real_type Cg = 0.0;        // Grid/PCC capacitance [F]
    ODE::real_type Rg = 0.0;        // Grid resistance [ohm]
    ODE::real_type egd = 0.0;       // d-axis component of grid voltage [V]
    ODE::real_type egq = 0.0;       // q-axis component of grid voltage [V]

    // CSC parameters
    ODE::real_type Rf = 0.0;        // AC-side filter resistance [ohm]
    ODE::real_type Lf = 0.0;        // AC-side filter inductance [H]
    ODE::real_type Cf = 0.0;        // AC-side filter capacitance [F]
    ODE::real_type Ldc = 0.0;       // DC-side filter inductance [H]
    ODE::real_type RL = 0.0;        // Equivalent DC-side load resistance [ohm]

    // Inner loop controller
    ODE::real_type kp1 = 0.0;       // Stage-1 proportional gain
    ODE::real_type ki1 = 0.0;       // Stage-1 integral gain
    ODE::real_type kp2 = 0.0;       // Stage-2 proportional gain
    ODE::real_type ki2 = 0.0;       // Stage-2 integral gain
    ODE::real_type Idc_min = 1.0;   // Minimum DC current threshold for numerical protection [A]

    // Outer loop controller
    ODE::real_type kpO = 0.0;       // Outer-loop proportional gain
    ODE::real_type kiO = 0.0;       // Outer-loop integral gain
    ODE::real_type V2_min = 1.0;    // Minimum squared-voltage threshold for numerical protection [V^2]

    // References
    ODE::StairSignal idc_ref;       // Piecewise-constant DC current reference [A]
    ODE::StairSignal Q_ref;         // Piecewise-constant reactive power reference [var]

    // Initial conditions
    ODE::vec_type x0;               // Initial state vector

};

/**
 * @brief Index map for the dynamic state vector.
 *
 * @details
 * The state vector is ordered as
 *
 * \f[
 * x =
 * \begin{bmatrix}
 * i_{gd} &
 * i_{gq} &
 * e_d &
 * e_q &
 * i_d &
 * i_q &
 * v_d &
 * v_q &
 * i_{stk} &
 * \hat{\theta} &
 * \xi_{\mathrm{pll}} &
 * \xi_{\mathrm{ucd}} &
 * \xi_{\mathrm{ucq}} &
 * \xi_{\mathrm{isd}} &
 * \xi_{\mathrm{isq}} &
 * \xi_{\mathrm{idc2}}
 * \end{bmatrix}^{\top}.
 * \f]
 */

enum StateIndex
{
    IGD = 0, IGQ, ED, EQ,               // AC source and PCC states in dq coordinates
    ID, IQ, VD, VQ, ISTK,               // CSC AC-side states and DC stack current
    THETA_HAT, XI_PLL,                  // PLL estimated angle and PLL integral state
    XI_UCD, XI_UCQ, XI_ISD, XI_ISQ,     // Inner-loop integral states
    XI_IDC2,                            // Outer-loop integral state
    NSTATES                             // Total number of dynamic states
};

/**
 * @brief Index map for algebraic outputs exported to the CSV file.
 *
 * @details
 * These outputs are not dynamic states. They are reconstructed from the state
 * vector and the controller algebraic variables. They include abc currents,
 * abc voltages, references, active power, and reactive power.
 */

enum OutputIndex
{
    IA = 0, IB, IC, VA, VB, VC,         // Three-phase currents and voltages in abc coordinates
    ID_R, IQ_R, IDC_REF, Q_REF, PREF,   // Current references, DC current reference, reactive power reference, and active power reference
    P_OUT, Q_OUT,                       // Active and reactive power outputs
    NOUTPUTS                            // Total number of algebraic outputs
};

/**
 * @brief PLL algebraic output variables.
 */

struct PLLOutput
{
    ODE::real_type e_pll = 0.0;         // Normalized PLL error
    ODE::real_type w_hat = 0.0;         // Estimated angular frequency [rad/s]
};

/**
 * @brief External reference values evaluated at a given simulation time.
 */

struct ReferenceOutput
{
    ODE::real_type idc_ref = 0.0;       // DC current reference evaluated at the current time [A]
    ODE::real_type Qref = 0.0;          // Reactive power reference evaluated at the current time [var]
};

/**
 * @brief Outer DC-current loop algebraic output variables.
 */

struct OuterLoopOutput
{
    ODE::real_type Eidc2 = 0.0;         // DC current tracking error
    ODE::real_type Pref = 0.0;          // Active power reference generated by the outer loop [W]
};

/**
 * @brief Inner current and voltage control algebraic output variables.
 */

struct InnerLoopOutput
{
    ODE::real_type idr = 0.0;           // d-axis current reference [A]
    ODE::real_type iqr = 0.0;           // q-axis current reference [A]

    // Stage 2
    ODE::real_type Eid = 0.0;           // d-axis current error
    ODE::real_type Eiq = 0.0;           // q-axis current error
    ODE::real_type vdr = 0.0;           // d-axis voltage reference generated by stage 2 [V]
    ODE::real_type vqr = 0.0;           // q-axis voltage reference generated by stage 2 [V]

    // Stage 1
    ODE::real_type Evd = 0.0;           // d-axis voltage error
    ODE::real_type Evq = 0.0;           // q-axis voltage error
    ODE::real_type md = 0.0;            // d-axis modulation signal
    ODE::real_type mq = 0.0;            // q-axis modulation signal

};

/**
 * @brief Complete set of controller algebraic variables.
 */

struct ControlOutput
{
    ODE::real_type e_pll = 0.0;         // Normalized PLL error
    ODE::real_type w_hat = 0.0;         // Estimated angular frequency [rad/s]

    ODE::real_type idc_ref = 0.0;       // DC current reference [A]
    ODE::real_type Qref = 0.0;          // Reactive power reference [var]
    ODE::real_type Pref = 0.0;          // Active power reference [W]

    ODE::real_type idr = 0.0;           // d-axis current reference [A]
    ODE::real_type iqr = 0.0;           // q-axis current reference [A]

    ODE::real_type Eid = 0.0;           // d-axis current error
    ODE::real_type Eiq = 0.0;           // q-axis current error
    ODE::real_type Eidc2 = 0.0;         // DC current tracking error

    ODE::real_type vdr = 0.0;           // d-axis voltage reference [V]
    ODE::real_type vqr = 0.0;           // q-axis voltage reference [V]

    ODE::real_type Evd = 0.0;           // d-axis voltage error
    ODE::real_type Evq = 0.0;           // q-axis voltage error

    ODE::real_type md = 0.0;            // d-axis modulation signal
    ODE::real_type mq = 0.0;            // q-axis modulation signal

};

/**
 * @brief ODE wrapper for the CSC-RL model.
 *
 * @details
 * This class exposes the CSC-RL model through the 
 * ODE::ODE_Problem_base interface. It does not implement the model equations
 * directly; instead, it delegates the dynamic equations and control laws to
 * CSC_RL_Model.
 *
 * Its main responsibilities are:
 *
 * - validating parameters,
 * - providing state names and initial conditions,
 * - evaluating the right-hand side,
 * - reconstructing algebraic outputs for CSV export.
 */

class CSC_RL : public ODE::ODE_Problem_base
{
    private:
        const CSC_RL_Model& model_;         // Reference to the model object containing equations and control laws
        CSC_RL_Parameters p_;               // Local copy of the CSC-RL parameter set
        void validate_parameters() const;   // Validate physical parameters, references, and initial conditions

    public:

        CSC_RL(const CSC_RL_Model& model, const CSC_RL_Parameters& parameters);  // Constructor using an external model object
        explicit CSC_RL(const CSC_RL_Parameters& parameters);                    // Constructor using a default internal model object

        ODE::integer n() const override;                                          // Return the number of dynamic states
        ODE::real_type t0() const override;                                       // Return the initial simulation time
        ODE::real_type tf() const override;                                       // Return the final simulation time
        ODE::real_type initial_condition(ODE::integer i) const override;          // Return the initial condition of state i
        std::string state_name(ODE::integer i) const override;                    // Return the name of state i for CSV output

        void rhs(ODE::real_type t, const ODE::vec_type& x, ODE::vec_type& dxdt) const override; // Compute state derivatives dxdt = f(t, x)

        ODE::integer n_outputs() const override;                                  // Return the number of algebraic outputs
        std::string output_name(ODE::integer i) const override;                   // Return the name of algebraic output i
        ODE::real_type output(ODE::real_type t, const ODE::vec_type& x, ODE::integer i) const override; // Return the value of algebraic output i
        ODE::vec_type outputs(ODE::real_type t, const ODE::vec_type& x) const override; // Return all algebraic outputs evaluated at time t and state x

};
