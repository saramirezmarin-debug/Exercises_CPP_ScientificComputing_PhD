#include "CSC_RL.hh"
#include "CSC_RL_Model.hh"

#include <array>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace
{
    constexpr ODE::real_type sqrt_3    = 1.7320508075688772935;
    constexpr ODE::real_type dq0_scale = 0.8164965809277260327; // sqrt(2/3)

    ODE::real_type abc_output_component(ODE::real_type d,
                                        ODE::real_type q,
                                        ODE::real_type theta,
                                        ODE::integer component)
    {
        const ODE::real_type ct = std::cos(theta);
        const ODE::real_type st = std::sin(theta);

        switch (component)
        {
        case 0:
            return dq0_scale * (ct * d - st * q);
        case 1:
            return dq0_scale * ((-0.5 * ct + 0.5 * sqrt_3 * st) * d - (-0.5 * st - 0.5 * sqrt_3 * ct) * q);
        case 2:
            return dq0_scale * ((-0.5 * ct - 0.5 * sqrt_3 * st) * d - (-0.5 * st + 0.5 * sqrt_3 * ct) * q);
        default:
            throw std::out_of_range("CSC_RL: invalid abc component index");
        }
    }
}

namespace
{
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

void CSC_RL::validate_parameters() const
{
    if (p_.Lg <= 0.0)       throw std::runtime_error("CSC_RL: Lg must be positive");
    if (p_.Cg <= 0.0)       throw std::runtime_error("CSC_RL: Cg must be positive");
    if (p_.Lf <= 0.0)       throw std::runtime_error("CSC_RL: Lf must be positive");
    if (p_.Cf <= 0.0)       throw std::runtime_error("CSC_RL: Cf must be positive");
    if (p_.Ldc <= 0.0)      throw std::runtime_error("CSC_RL: Ldc must be positive");
    if (p_.RL < 0.0)        throw std::runtime_error("CSC_RL: RL must be non-negative");
    if (p_.tf < p_.t0)      throw std::runtime_error("CSC_RL: tf must be greater than or equal to t0");
    if (p_.Idc_min <= 0.0)  throw std::runtime_error("CSC_RL: Idc_min must be positive");
    if (p_.V2_min <= 0.0)   throw std::runtime_error("CSC_RL: V2_min must be positive");
    if (p_.Vdq_nom == 0.0)  throw std::runtime_error("CSC_RL: Vdq_nom must be different to zero");

    if (p_.x0.size() != NSTATES) throw std::runtime_error("CSC_RL: x0 has wrong size");


    // Validate the generic stair references
    p_.idc_ref.validate("CSC_RL: idc_ref");
    p_.Q_ref.validate("CSC_RL: Q_ref");
}

ODE::integer CSC_RL::n() const {return NSTATES;}
ODE::real_type CSC_RL::t0() const {return p_.t0;}
ODE::real_type CSC_RL::tf() const {return p_.tf;}

ODE::real_type CSC_RL::initial_condition(ODE::integer i) const
{
    if (i < 0 || i >= NSTATES)
    {
        throw std::out_of_range("CSC_RL: invalid state index");
    }

    return p_.x0(i);
}

void CSC_RL::rhs(ODE::real_type t, const ODE::vec_type& x, ODE::vec_type& dxdt) const
{
    model_.rhs(p_, t, x, dxdt);
}

// ------------------------------------------------------------
// Number of algebraic outputs exported to the CSV.
// ------------------------------------------------------------
ODE::integer CSC_RL::n_outputs() const {return NOUTPUTS;}

// ------------------------------------------------------------
// Names of the algebraic outputs.
// ------------------------------------------------------------
std::string CSC_RL::output_name(ODE::integer i) const
{
    static const std::array<std::string, NOUTPUTS> names = 
    {
        "ia", "ib", "ic",
        "va", "vb", "vc",
        "id_r", "iq_r",
        "idc_ref", "Q_ref", "Pref",
        "P", "Q"
    };

    if (i >= 0 && i < NOUTPUTS)     return names[i];

    throw std::out_of_range("CSC_RL: invalid output index");
}

std::string CSC_RL::state_name(ODE::integer i) const
{
    static const std::array<std::string, NSTATES> names = 
    {
        "igd", "igq", "ed", "eq",               // AC Source 
        "id", "iq", "vd", "vq", "istk",         // CSC
        "theta_hat", "xi_pll",                  // PLL
        "xi_ucd", "xi_ucq", "xi_isd", "xi_isq", // Inner
        "xi_idc2"                               // Outer
    };

    if (i >= 0 && i < NSTATES)      return names[i];

    throw std::out_of_range("CSC_RL: invalid state index");
}

ODE::vec_type CSC_RL::outputs(ODE::real_type t, const ODE::vec_type& x) const
{
    ODE::vec_type y(NOUTPUTS);

    const ODE::real_type ed = x(ED);
    const ODE::real_type eq = x(EQ);
    const ODE::real_type id = x(ID);
    const ODE::real_type iq = x(IQ);
    const ODE::real_type vd = x(VD);
    const ODE::real_type vq = x(VQ);
    const ODE::real_type theta = x(THETA_HAT);

    // const ControlOutput c = compute_control(t, x);
    const auto input = model_.input_from_parameters<ODE::real_type>(p_, t);
    const auto c = model_.control(p_, x, input);

    y(IA) = abc_output_component(id, iq, theta, 0);
    y(IB) = abc_output_component(id, iq, theta, 1);
    y(IC) = abc_output_component(id, iq, theta, 2);

    y(VA) = abc_output_component(vd, vq, theta, 0);
    y(VB) = abc_output_component(vd, vq, theta, 1);
    y(VC) = abc_output_component(vd, vq, theta, 2);

    y(ID_R) = c.idr;
    y(IQ_R) = c.iqr;

    y(IDC_REF) = c.idc_ref;
    y(Q_REF)   = c.Qref;
    y(PREF)    = c.Pref;

    y(P_OUT) = ed * id + eq * iq;
    y(Q_OUT) = eq * id - ed * iq;

    return y;
}

// ------------------------------------------------------------
//   iabc = dq0_to_abc([id, iq, 0], theta_hat)
//   eabc = dq0_to_abc([ed, eq, 0], theta_hat)
// ------------------------------------------------------------
ODE::real_type CSC_RL::output(ODE::real_type t, const ODE::vec_type& x, ODE::integer i) const
{
    if (i < 0 || i >= NOUTPUTS)
    {
        throw std::out_of_range("CSC_RL: invalid output index");
    }

    return outputs(t, x)(i);
}
