#include "CSC_RL.hh"
#include "CSC_RL_Model.hh"

#include <array>
#include <algorithm>
#include <cmath>
#include <stdexcept>

/**
 * @file CSC_RL.cc
 * @brief Implementation of the CSC_RL ODE wrapper.
 *
 * @details
 * This file connects the mathematical CSC-RL model to the generic RK4 solver.
 * It also reconstructs abc quantities from dq variables and computes algebraic
 * output variables exported to the CSV file.
 */

// ============================================================
// CSC-RL ODE problem interface.
//
// The actual model equations are implemented in CSC_RL_Model.hh.
// This file only declares the problem wrapper used by the ODE solver.
// ============================================================

namespace
{
    constexpr ODE::real_type sqrt_3    = 1.7320508075688772935;     // sqr(3), used in the inverse dq-to-abc transformation.
    constexpr ODE::real_type dq0_scale = 0.8164965809277260327;     // Power-invariant dq0-to-abc scaling factor: sqrt(2/3).

    /**
     * @brief Reconstructs one abc component from dq variables.
     *
     * @details
     * The implementation assumes balanced zero-sequence-free quantities and uses
     * the power-invariant inverse dq transformation. For a dq vector
     * \f$(d,q)\f$, the abc components are reconstructed using the PLL angle
     * \f$\hat{\theta}\f$.
     *
     * @param d d-axis component.
     * @param q q-axis component.
     * @param theta Electrical angle used by the inverse transformation.
     * @param component Phase index: 0 for phase a, 1 for phase b, 2 for phase c.
     * @return Selected abc phase component.
     */

    ODE::real_type abc_output_component(ODE::real_type d,           // d-axis component
                                        ODE::real_type q,           // q-axis component
                                        ODE::real_type theta,       // angle used for transformation
                                        ODE::integer component)     // 0 -> a, 1 -> b, 2 -> c
    {
        const ODE::real_type ct = std::cos(theta);                  // Precompute cosine of the electrical angle.
        const ODE::real_type st = std::sin(theta);                  // Precompute sine of the electrical angle.

        // Assumes zero-sequence-free balanced quantities.
        switch (component)  // Select the requested abc phase.
        {
        case 0:             // Phase-a component.
            return dq0_scale * (ct * d - st * q);
        case 1:             // Phase-b component.
            return dq0_scale * ((-0.5 * ct + 0.5 * sqrt_3 * st) * d - (-0.5 * st - 0.5 * sqrt_3 * ct) * q);
        case 2:             // Phase-c component.
            return dq0_scale * ((-0.5 * ct - 0.5 * sqrt_3 * st) * d - (-0.5 * st + 0.5 * sqrt_3 * ct) * q);
        default:            // Reject invalid phase indices.
            throw std::out_of_range("CSC_RL: invalid abc component index");
        }
    }
}

namespace
{

    /**
     * @brief Provides a default static CSC_RL_Model instance.
     *
     * @return Reference to an internal stateless CSC_RL_Model object.
     */
    const CSC_RL_Model& default_csc_rl_model()
    {
        static CSC_RL_Model model;
        return model;
    }
}


CSC_RL::CSC_RL(const CSC_RL_Model& model,
               const CSC_RL_Parameters& parameters)
    : model_(model)
    , p_(parameters)
{
    validate_parameters();
}

CSC_RL::CSC_RL(const CSC_RL_Parameters& parameters)
    : model_(default_csc_rl_model())
    , p_(parameters)
{
    validate_parameters();
}

// Check physical and numerical consistency of the CSC-RL parameter set.
void CSC_RL::validate_parameters() const
{
    if (p_.tf < p_.t0)      throw std::runtime_error("CSC_RL: tf must be greater than or equal to t0");
    if (p_.Idc_min == 0.0)  throw std::runtime_error("CSC_RL: Idc_min must be different to zero");
    if (p_.V2_min == 0.0)   throw std::runtime_error("CSC_RL: V2_min must be different to zero");
    if (p_.Vdq_nom == 0.0)  throw std::runtime_error("CSC_RL: Vdq_nom must be different to zero");
    if (p_.x0.size() != NSTATES) throw std::runtime_error("CSC_RL: x0 has wrong size");

    // Validate the generic stair references
    p_.idc_ref.validate("CSC_RL: idc_ref");
    p_.Q_ref.validate("CSC_RL: Q_ref");
}

ODE::integer CSC_RL::n() const {return NSTATES;}        // Return the number of dynamic states.
ODE::real_type CSC_RL::t0() const {return p_.t0;}       // Return the initial simulation time.
ODE::real_type CSC_RL::tf() const {return p_.tf;}       // Return the final simulation time.

ODE::real_type CSC_RL::initial_condition(ODE::integer i) const                          // Return the initial condition of state i.
{
    if (i < 0 || i >= NSTATES)  throw std::out_of_range("CSC_RL: invalid state index"); // Reject invalid state indices.

    return p_.x0(i);
}

// Compute the state derivatives by delegating the calculation to CSC_RL_Model.
void CSC_RL::rhs(ODE::real_type t, const ODE::vec_type& x, ODE::vec_type& dxdt) const
{
    model_.rhs(p_, t, x, dxdt);
}

ODE::integer CSC_RL::n_outputs() const {return NOUTPUTS;}   // Return the number of algebraic outputs exported to the CSV file.

std::string CSC_RL::output_name(ODE::integer i) const       // Return the CSV header name of algebraic output i.
{
    // Fixed list of algebraic output names.
    static const std::array<std::string, NOUTPUTS> names = 
    {
        "ia", "ib", "ic",
        "va", "vb", "vc",
        "id_r", "iq_r",
        "idc_ref", "Q_ref", "Pref",
        "P", "Q"
    };

    if (i >= 0 && i < NOUTPUTS)     return names[i];            // Return the requested name if the output index is valid.
    throw std::out_of_range("CSC_RL: invalid output index");    // Reject invalid output indices.
}

std::string CSC_RL::state_name(ODE::integer i) const            // Return the CSV header name of state i.
{
    static const std::array<std::string, NSTATES> names = 
    {
        "igd", "igq", "ed", "eq",                           // AC Source 
        "id", "iq", "vd", "vq", "istk",                     // CSC
        "theta_hat", "xi_pll",                              // PLL
        "xi_ucd", "xi_ucq", "xi_isd", "xi_isq",             // Inner
        "xi_idc2"                                           // Outer
    };

    if (i >= 0 && i < NSTATES)      return names[i];        // Return the requested name if the state index is valid.
    throw std::out_of_range("CSC_RL: invalid state index"); // Reject invalid state indices.
}

/**
 * @brief Computes all algebraic outputs exported to the CSV file.
 *
 * @details
 * The function reconstructs abc currents and voltages from dq variables and
 * evaluates the control algebraic quantities. The active and reactive powers
 * are computed as
 *
 * \f[
 * P = e_d i_d + e_q i_q,
 * \f]
 *
 * \f[
 * Q = e_q i_d - e_d i_q.
 * \f]
 *
 * @param t Current simulation time.
 * @param x Current state vector.
 * @return Vector containing all algebraic outputs.
 */

// Compute all algebraic outputs exported to the CSV file.
ODE::vec_type CSC_RL::outputs(ODE::real_type t, const ODE::vec_type& x) const
{
    ODE::vec_type y(NOUTPUTS);                          // Allocate the algebraic output vector.

    const ODE::real_type ed = x(ED);                    // PCC voltage d-axis.
    const ODE::real_type eq = x(EQ);                    // PCC voltage q-axis.
    const ODE::real_type id = x(ID);                    // CSC AC-side current d-axis.
    const ODE::real_type iq = x(IQ);                    // CSC AC-side current q-axis.
    const ODE::real_type vd = x(VD);                    // CSC AC-side voltage q-axis.
    const ODE::real_type vq = x(VQ);                    // CSC AC-side voltage q-axis.
    const ODE::real_type theta = x(THETA_HAT);          // Electrical angle estimated by the PLL.

    const auto input = model_.input_from_parameters<ODE::real_type>(p_, t);     // Evaluate time-dependent model inputs and references.
    const auto c = model_.control(p_, x, input);                                // Compute control variables from the current state and references.

    // Reconstruct three-phase CSC currents from dq components.
    y(IA) = abc_output_component(id, iq, theta, 0);
    y(IB) = abc_output_component(id, iq, theta, 1);
    y(IC) = abc_output_component(id, iq, theta, 2);

    // Reconstruct three-phase CSC voltages from dq components.
    y(VA) = abc_output_component(vd, vq, theta, 0);
    y(VB) = abc_output_component(vd, vq, theta, 1);
    y(VC) = abc_output_component(vd, vq, theta, 2);

    y(ID_R) = c.idr;                    // d-axis current reference generated by controller
    y(IQ_R) = c.iqr;                    // q-axis current reference generated by controller

    y(IDC_REF) = c.idc_ref;             // DC current reference (external)
    y(Q_REF)   = c.Qref;                // Reactive power reference (external)
    y(PREF)    = c.Pref;                // Active power reference

    y(P_OUT) = ed * id + eq * iq;       // Compute active power
    y(Q_OUT) = eq * id - ed * iq;       // Compute reactive power

    return y;
}

// Return one algebraic output selected by index.
ODE::real_type CSC_RL::output(ODE::real_type t, const ODE::vec_type& x, ODE::integer i) const
{
    if (i < 0 || i >= NOUTPUTS)     throw std::out_of_range("CSC_RL: invalid output index");    // Reject invalid output indices.
    
    return outputs(t, x)(i);
}
