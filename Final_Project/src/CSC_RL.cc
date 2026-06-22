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
    case IGD:  return p_.igd0;
    case IGQ:  return p_.igq0;
    case ED:  return p_.ed0;
    case EQ:  return p_.eq0;
    case ID:  return p_.id0;
    case IQ:  return p_.iq0;
    case VD:  return p_.vd0;
    case VQ:  return p_.vq0;
    case ISTK:  return p_.istk0;
    case THETA_HAT:  return p_.theta_hat0;
    case XI_PLL: return p_.xi_pll0;
    case XI_ECD: return p_.xi_ucd0;
    case XI_ECQ: return p_.xi_ucq0;
    case XI_ISD: return p_.xi_isd0;
    case XI_ISQ: return p_.xi_isq0;
    case XI_IDC2: return p_.xi_idc2_0;
    default: throw std::out_of_range("CSC_RL: invalid state index");
    }
}

ControlOutput CSC_RL::compute_control(ODE::real_type t,
                                      const ODE::vec_type& x) const
{
    ControlOutput c;

    const ODE::real_type ed   = x(ED);
    const ODE::real_type eq   = x(EQ);
    const ODE::real_type id   = x(ID);
    const ODE::real_type iq   = x(IQ);
    const ODE::real_type vd   = x(VD);
    const ODE::real_type vq   = x(VQ);
    const ODE::real_type istk = x(ISTK);

    const ODE::real_type xi_pll  = x(XI_PLL);
    const ODE::real_type xi_ucd  = x(XI_ECD);
    const ODE::real_type xi_ucq  = x(XI_ECQ);
    const ODE::real_type xi_isd  = x(XI_ISD);
    const ODE::real_type xi_isq  = x(XI_ISQ);
    const ODE::real_type xi_idc2 = x(XI_IDC2);

    // ------------------------------------------------------------
    // PLL
    // ------------------------------------------------------------
    c.e_pll = eq / p_.Vdq_nom;
    c.w_hat = p_.w0_pll + p_.kp_pll * c.e_pll + p_.ki_pll * xi_pll;

    // ------------------------------------------------------------
    // Outer loop: i_dc reference -> active power reference
    // ------------------------------------------------------------
    c.idc_ref = p_.idc_ref.value(t);
    c.Qref    = p_.Q_ref.value(t);

    c.Eidc2 = c.idc_ref * c.idc_ref - istk * istk;
    c.Pref  = p_.kpO * c.Eidc2 + p_.kiO * xi_idc2;

    // ------------------------------------------------------------
    // P/Q reference -> current references
    // ------------------------------------------------------------
    const ODE::real_type V2_eff =
        std::max(ed * ed + eq * eq, p_.V2_min);

    c.idr = (c.Pref * ed + c.Qref * eq) / V2_eff;
    c.iqr = (c.Pref * eq - c.Qref * ed) / V2_eff;

    // ------------------------------------------------------------
    // Inner loop stage 2: current error -> voltage reference
    // ------------------------------------------------------------
    c.Eid = c.idr - id;
    c.Eiq = c.iqr - iq;

    c.vdr = ed + p_.Lf * c.w_hat * iq
          - (p_.kp2 * c.Eid + p_.ki2 * xi_isd);

    c.vqr = eq - p_.Lf * c.w_hat * id
          - (p_.kp2 * c.Eiq + p_.ki2 * xi_isq);

    // ------------------------------------------------------------
    // Inner loop stage 1: voltage error -> modulation index
    // ------------------------------------------------------------
    c.Evd = c.vdr - vd;
    c.Evq = c.vqr - vq;

    const ODE::real_type Idc_eff =
        std::max(std::abs(istk), p_.Idc_min);

    c.md = 0.0;
    c.mq = 0.0;

    if (std::abs(istk) >= p_.Idc_min)
    {
        c.md = (id + p_.Cf * c.w_hat * iq
             - (p_.kp1 * c.Evd + p_.ki1 * xi_ucd)) / Idc_eff;

        c.mq = (iq - p_.Cf * c.w_hat * id
             - (p_.kp1 * c.Evq + p_.ki1 * xi_ucq)) / Idc_eff;
    }

    return c;
}

void CSC_RL::rhs(ODE::real_type t,
                 const ODE::vec_type& x,
                 ODE::vec_type& dxdt) const
{
    const ODE::real_type igd  = x(IGD);
    const ODE::real_type igq  = x(IGQ);
    const ODE::real_type ed   = x(ED);
    const ODE::real_type eq   = x(EQ);

    const ODE::real_type id   = x(ID);
    const ODE::real_type iq   = x(IQ);
    const ODE::real_type vd   = x(VD);
    const ODE::real_type vq   = x(VQ);
    const ODE::real_type istk = x(ISTK);

    const ControlOutput c = compute_control(t, x);

    // ------------------------------------------------------------
    // PCC dynamics
    // ------------------------------------------------------------
    dxdt(IGD) = (p_.egd - p_.Rg * igd - ed
              + c.w_hat * p_.Lg * igq) / p_.Lg;

    dxdt(IGQ) = (p_.egq - p_.Rg * igq - eq
              - c.w_hat * p_.Lg * igd) / p_.Lg;

    dxdt(ED) = (igd - id
             + c.w_hat * p_.Cg * eq) / p_.Cg;

    dxdt(EQ) = (igq - iq
             - c.w_hat * p_.Cg * ed) / p_.Cg;

    // ------------------------------------------------------------
    // CSC AC-side and DC-side dynamics
    // ------------------------------------------------------------
    dxdt(ID) = (ed - vd - p_.Rf * id
             + c.w_hat * p_.Lf * iq) / p_.Lf;

    dxdt(IQ) = (eq - vq - p_.Rf * iq
             - c.w_hat * p_.Lf * id) / p_.Lf;

    dxdt(VD) = (id - c.md * istk
             + c.w_hat * p_.Cf * vq) / p_.Cf;

    dxdt(VQ) = (iq - c.mq * istk
             - c.w_hat * p_.Cf * vd) / p_.Cf;

    dxdt(ISTK) = (c.md * vd + c.mq * vq
               - p_.RL * istk) / p_.Ldc;

    // ------------------------------------------------------------
    // PLL states
    // ------------------------------------------------------------
    dxdt(THETA_HAT) = c.w_hat;
    dxdt(XI_PLL)    = c.e_pll;

    // ------------------------------------------------------------
    // Controller integrator states
    // ------------------------------------------------------------
    dxdt(XI_ECD)  = c.Evd;
    dxdt(XI_ECQ)  = c.Evq;
    dxdt(XI_ISD)  = c.Eid;
    dxdt(XI_ISQ)  = c.Eiq;
    dxdt(XI_IDC2) = c.Eidc2;
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
    return NOUTPUTS;
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

    if (i >= 0 && i < NOUTPUTS) 
    {
        return names[i];
    }

    throw std::out_of_range("CSC_RL: invalid output index");
}

std::string CSC_RL::state_name(ODE::integer i) const
{
    static const std::string names[] = {
        "igd",
        "igq",
        "ed",
        "eq",
        "id",
        "iq",
        "vd",
        "vq",
        "istk",
        "theta_hat",
        "xi_pll",
        "xi_ucd",
        "xi_ucq",
        "xi_isd",
        "xi_isq",
        "xi_idc2"
    };

    if (i >= 0 && i < NSTATES)
    {
        return names[i];
    }

    throw std::out_of_range("CSC_RL: invalid state index");
}

// ------------------------------------------------------------
//   iabc = dq0_to_abc([id, iq, 0], theta_hat)
//   eabc = dq0_to_abc([ed, eq, 0], theta_hat)
// ------------------------------------------------------------
ODE::real_type CSC_RL::output(ODE::real_type t,
                              const ODE::vec_type& x,
                              ODE::integer i) const
{
    const ODE::real_type ed = x(ED);
    const ODE::real_type eq = x(EQ);
    const ODE::real_type id = x(ID);
    const ODE::real_type iq = x(IQ);

    if (i >= 0 && i <= 2)
    {
        return abc_output_component(id, iq, x(THETA_HAT), i);
    }

    if (i >= 3 && i <= 5)
    {
        return abc_output_component(x(VD), x(VQ), x(THETA_HAT), i - 3);
    }

    if (i == P_OUT) return ed * id + eq * iq;
    if (i == Q_OUT) return eq * id - ed * iq;

    const ODE::real_type idc_ref = p_.idc_ref.value(t);
    if (i == IDC_REF) return idc_ref;

    const ODE::real_type Qref = p_.Q_ref.value(t);
    if (i == Q_REF) return Qref;

    const ODE::real_type istk = x(ISTK);
    const ODE::real_type Pref =
        p_.kpO * (idc_ref * idc_ref - istk * istk) + p_.kiO * x(XI_IDC2);
    if (i == PREF) return Pref;

    const ODE::real_type V2_eff = std::max(ed * ed + eq * eq, p_.V2_min);
    if (i == ID_R) return (Pref * ed + Qref * eq) / V2_eff;
    if (i == IQ_R) return (Pref * eq - Qref * ed) / V2_eff;

    throw std::out_of_range("CSC_RL: invalid output index");
}
