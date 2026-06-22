#include "CSC_RL.hh"

namespace
{
    constexpr ODE::real_type sqrt_3 = 1.7320508075688772935;
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
            return dq0_scale * ((-0.5 * ct + 0.5 * sqrt_3 * st) * d
                              - (-0.5 * st - 0.5 * sqrt_3 * ct) * q);
        case 2:
            return dq0_scale * ((-0.5 * ct - 0.5 * sqrt_3 * st) * d
                              - (-0.5 * st + 0.5 * sqrt_3 * ct) * q);
        default:
            throw std::out_of_range("CSC_RL: invalid abc component index");
        }
    }
}

CSC_RL::CSC_RL(const CSC_RL_Parameters& parameters) : p_(parameters)
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

    // Validate the generic stair references
    // p_.id_r.validate("CSC_RL: id_r");
    // p_.iq_r.validate("CSC_RL: iq_r");
    p_.idc_ref.validate("CSC_RL: idc_ref");
    p_.Q_ref.validate("CSC_RL: Q_ref");
}

ODE::integer CSC_RL::n() const
{
    return 16;
}

ODE::real_type CSC_RL::t0() const
{
    return p_.t0;
}

ODE::real_type CSC_RL::tf() const
{
    return p_.tf;
}

ODE::real_type CSC_RL::initial_condition(ODE::integer i) const
{
    switch (i)
    {
    case 0:  return p_.igd0;
    case 1:  return p_.igq0;
    case 2:  return p_.ed0;
    case 3:  return p_.eq0;
    case 4:  return p_.id0;
    case 5:  return p_.iq0;
    case 6:  return p_.vd0;
    case 7:  return p_.vq0;
    case 8:  return p_.istk0;
    case 9:  return p_.theta_hat0;
    case 10: return p_.xi_pll0;
    case 11: return p_.xi_ucd0;
    case 12: return p_.xi_ucq0;
    case 13: return p_.xi_isd0;
    case 14: return p_.xi_isq0;
    case 15: return p_.xi_idc2_0;
    default: throw std::out_of_range("CSC_RL: invalid state index");
    }
}

void CSC_RL::rhs(ODE::real_type t,
                 const ODE::vec_type& x,
                 ODE::vec_type& dxdt) const
{
    (void)t;

    // ------------------------------------------------------------
    // State vector
    // x = [igd, igq, ed, eq, id, iq, vd, vq, istk, theta_hat, xi_pll, xi_ucd, xi_ucq, xi_isd, xi_isq, xi_idc2]^T
    // ------------------------------------------------------------

    // PCC states
    const ODE::real_type igd  = x(0);
    const ODE::real_type igq  = x(1);
    const ODE::real_type ed   = x(2);
    const ODE::real_type eq   = x(3);

    // CSC states
    const ODE::real_type id   = x(4);
    const ODE::real_type iq   = x(5);
    const ODE::real_type vd   = x(6);
    const ODE::real_type vq   = x(7);
    const ODE::real_type istk = x(8);

    const ODE::real_type xi_pll    = x(10);

    // Inner loop controller states
    const ODE::real_type xi_ucd = x(11);
    const ODE::real_type xi_ucq = x(12);
    const ODE::real_type xi_isd = x(13);
    const ODE::real_type xi_isq = x(14);

    // Outer loop controller state
    const ODE::real_type xi_idc2 = x(15);

    // PLL
    const ODE::real_type e_pll = eq / p_.Vdq_nom; 
    const ODE::real_type w_hat = p_.w0_pll + p_.kp_pll * e_pll + p_.ki_pll * xi_pll;

    // ------------------------------------------------------------
    // OUTER LOOP CONTROLLER
    // ------------------------------------------------------------
    const ODE::real_type idc_ref = p_.idc_ref.value(t);
    const ODE::real_type Qref   = p_.Q_ref.value(t);

    const ODE::real_type Eidc2 = idc_ref*idc_ref - istk*istk; // Error in DC current squared
    const ODE::real_type Pref = p_.kpO * Eidc2 + p_.kiO * xi_idc2; // Active power reference 

    const ODE::real_type V2 = ed*ed + eq*eq; // PCC voltage magnitude square
    const ODE::real_type V2_eff = std::max(V2, p_.V2_min);

    const ODE::real_type idr = (Pref * ed + Qref * eq) / V2_eff;
    const ODE::real_type iqr = (Pref * eq - Qref * ed) / V2_eff;

    // ------------------------------------------------------------
    // INNER LOOP CONTROLLER
    // ------------------------------------------------------------

    // ------------------------------------------------------------
    // Stage 2: current reference -> capacitor voltage reference
    // ------------------------------------------------------------
    const ODE::real_type Eid = idr - id;
    const ODE::real_type Eiq = iqr - iq;
    const ODE::real_type vdr = ed + p_.Lf * w_hat * iq - (p_.kp2 * Eid + p_.ki2 * xi_isd); 
    const ODE::real_type vqr = eq - p_.Lf * w_hat * id - (p_.kp2 * Eiq + p_.ki2 * xi_isq);

    // ------------------------------------------------------------
    // Stage 1: capacitor voltage reference -> modulation
    // ------------------------------------------------------------
    const ODE::real_type Evd = vdr - vd;
    const ODE::real_type Evq = vqr - vq;

    const ODE::real_type Idc_eff = std::max(std::abs(istk), p_.Idc_min); // Avoid division by zero or very small values
    ODE::real_type md = 0;
    ODE::real_type mq = 0;

    if (std::abs(istk) >= p_.Idc_min)
    {
        md = (id + p_.Cf * w_hat * iq - (p_.kp1 * Evd + p_.ki1 * xi_ucd)) / Idc_eff;
        mq = (iq - p_.Cf * w_hat * id - (p_.kp1 * Evq + p_.ki1 * xi_ucq)) / Idc_eff;
    }
    
    // ------------------------------------------------------------
    // PCC
    // ------------------------------------------------------------
    dxdt(0) = (p_.egd - p_.Rg * igd - ed + w_hat * p_.Lg * igq) / p_.Lg;
    dxdt(1) = (p_.egq - p_.Rg * igq - eq - w_hat * p_.Lg * igd) / p_.Lg;
    dxdt(2) = (igd - id + w_hat * p_.Cg * eq) / p_.Cg;
    dxdt(3) = (igq - iq - w_hat * p_.Cg * ed) / p_.Cg;

    // ------------------------------------------------------------
    // CSC
    // ------------------------------------------------------------
    dxdt(4) = (ed - vd - p_.Rf * id + w_hat * p_.Lf * iq) / p_.Lf;
    dxdt(5) = (eq - vq - p_.Rf * iq - w_hat * p_.Lf * id) / p_.Lf;
    dxdt(6) = (id - md * istk + w_hat * p_.Cf * vq) / p_.Cf;
    dxdt(7) = (iq - mq * istk - w_hat * p_.Cf * vd) / p_.Cf;
    dxdt(8) = (md * vd + mq * vq - p_.RL * istk) / p_.Ldc;

    // ------------------------------------------------------------
    // PLL states
    // ------------------------------------------------------------
    dxdt(9)  = w_hat;
    dxdt(10) = e_pll;

    // ------------------------------------------------------------
    // Inner loop controller states
    // ------------------------------------------------------------
    dxdt(11) = Evd;
    dxdt(12) = Evq; 
    dxdt(13) = Eid;
    dxdt(14) = Eiq;

    // ------------------------------------------------------------
    // Outer loop controller states
    // ------------------------------------------------------------
    dxdt(15) = Eidc2;

}

// ------------------------------------------------------------
// Number of algebraic outputs exported to the CSV.
//
// Outputs:
//   0 -> ia
//   1 -> ib
//   2 -> ic
//   3 -> va
//   4 -> vb
//   5 -> vc
//   6 -> id_r
//   7 -> iq_r
//   8 -> idc_ref
//   9 -> Qref
//   10 -> Pref
//   11 -> P
//   12 -> Q
// ------------------------------------------------------------
ODE::integer CSC_RL::n_outputs() const
{
    return 13;
}

// ------------------------------------------------------------
// Names of the algebraic outputs.
// ------------------------------------------------------------
std::string CSC_RL::output_name(ODE::integer i) const
{
    static const std::string names[] = {
        "ia", "ib", "ic",
        "va", "vb", "vc",
        "id_r", "iq_r",
        "idc_ref", "Q_ref", "Pref",
        "P", "Q"
    };

    if (i >= 0 && i < n_outputs()) return names[i];

    throw std::out_of_range("CSC_RL: invalid output index");
}

// ------------------------------------------------------------
//   iabc = dq0_to_abc([id, iq, 0], theta_hat)
//   eabc = dq0_to_abc([ed, eq, 0], theta_hat)
// ------------------------------------------------------------
ODE::real_type CSC_RL::output(ODE::real_type t,
                              const ODE::vec_type& x,
                              ODE::integer i) const
{
    const ODE::real_type ed = x(2);
    const ODE::real_type eq = x(3);
    const ODE::real_type id = x(4);
    const ODE::real_type iq = x(5);

    if (i >= 0 && i <= 2)
    {
        return abc_output_component(id, iq, x(9), i);
    }

    if (i >= 3 && i <= 5)
    {
        return abc_output_component(x(6), x(7), x(9), i - 3);
    }

    if (i == 11) return ed * id + eq * iq;
    if (i == 12) return eq * id - ed * iq;

    const ODE::real_type idc_ref = p_.idc_ref.value(t);
    if (i == 8) return idc_ref;

    const ODE::real_type Qref = p_.Q_ref.value(t);
    if (i == 9) return Qref;

    const ODE::real_type istk = x(8);
    const ODE::real_type Pref =
        p_.kpO * (idc_ref * idc_ref - istk * istk) + p_.kiO * x(15);
    if (i == 10) return Pref;

    const ODE::real_type V2_eff = std::max(ed * ed + eq * eq, p_.V2_min);
    if (i == 6) return (Pref * ed + Qref * eq) / V2_eff;
    if (i == 7) return (Pref * eq - Qref * ed) / V2_eff;

    throw std::out_of_range("CSC_RL: invalid output index");
}
