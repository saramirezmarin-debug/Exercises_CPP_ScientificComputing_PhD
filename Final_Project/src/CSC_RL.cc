#include "CSC_RL.hh"

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

CSC_RL::CSC_RL(const CSC_RL_Parameters& parameters) : p_(parameters) {validate_parameters();}

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

PLLOutput CSC_RL::compute_pll(const ODE::vec_type& x) const
{
    PLLOutput pll;

    const ODE::real_type eq = x(EQ);
    const ODE::real_type xi_pll = x(XI_PLL);

    pll.e_pll = eq / p_.Vdq_nom;
    pll.w_hat = p_.w0_pll + p_.kp_pll * pll.e_pll + p_.ki_pll * xi_pll;

    return pll;
}

ReferenceOutput CSC_RL::compute_references(ODE::real_type t) const
{
    ReferenceOutput ref;
    
    ref.idc_ref = p_.idc_ref.value(t);
    ref.Qref = p_.Q_ref.value(t);

    return ref;
}

OuterLoopOutput CSC_RL::compute_outer_loop(const ODE::vec_type& x, const ReferenceOutput& ref) const
{
    OuterLoopOutput outer;

    const ODE::real_type istk = x(ISTK);
    const ODE::real_type xi_idc2 = x(XI_IDC2);

    outer.Eidc2 = ref.idc_ref * ref.idc_ref - istk * istk;
    outer.Pref = p_.kpO * outer.Eidc2 + p_.kiO * xi_idc2;

    return outer;
}

InnerLoopOutput CSC_RL::compute_inner_loop(const ODE::vec_type& x, const ReferenceOutput& ref, const PLLOutput& pll, const OuterLoopOutput& outer) const
{
    InnerLoopOutput inner;

    const ODE::real_type ed = x(ED);
    const ODE::real_type eq = x(EQ);
    const ODE::real_type id = x(ID);
    const ODE::real_type iq = x(IQ);
    const ODE::real_type vd = x(VD);
    const ODE::real_type vq = x(VQ);
    const ODE::real_type istk = x(ISTK);

    const ODE::real_type xi_isd = x(XI_ISD);
    const ODE::real_type xi_isq = x(XI_ISQ);
    const ODE::real_type xi_ucd = x(XI_UCD);
    const ODE::real_type xi_ucq = x(XI_UCQ);

    const ODE::real_type V2_eff = std::max(ed * ed + eq * eq, p_.V2_min);
    const ODE::real_type Idc_eff = std::max(std::abs(istk), p_.Idc_min);

    inner.idr = (ed * outer.Pref + eq * ref.Qref) / V2_eff;
    inner.iqr = (eq * outer.Pref - ed * ref.Qref) / V2_eff;

    inner.Eid = inner.idr - id;
    inner.Eiq = inner.iqr - iq;

    inner.vdr = ed + pll.w_hat * p_.Lf * iq - p_.Rf * id - p_.kp2 * inner.Eid - p_.ki2 * xi_isd;
    inner.vqr = eq - pll.w_hat * p_.Lf * id - p_.Rf * iq - p_.kp2 * inner.Eiq - p_.ki2 * xi_isq;

    inner.Evd = inner.vdr - vd;
    inner.Evq = inner.vqr - vq;

    inner.md = (id + pll.w_hat * p_.Cf * vq - p_.kp1 * inner.Evd - p_.ki1 * xi_ucd) / Idc_eff;
    inner.mq = (iq - pll.w_hat * p_.Cf * vd - p_.kp1 * inner.Evq - p_.ki1 * xi_ucq) / Idc_eff;

    return inner;
}

ControlOutput CSC_RL::compute_control(ODE::real_type t, const ODE::vec_type& x) const
{
    ControlOutput c;

    const PLLOutput pll = compute_pll(x);
    const ReferenceOutput ref = compute_references(t);
    const OuterLoopOutput outer = compute_outer_loop(x, ref);
    const InnerLoopOutput inner = compute_inner_loop(x, ref, pll, outer);

    // PLL
    c.e_pll = pll.e_pll;
    c.w_hat = pll.w_hat;

    // References
    c.idc_ref = ref.idc_ref;
    c.Qref = ref.Qref;

    // Outer
    c.Eidc2 = outer.Eidc2;
    c.Pref = outer.Pref;

    // Inner
    c.idr = inner.idr;
    c.iqr = inner.iqr;

    c.Eid = inner.Eid;
    c.Eiq = inner.Eiq;
    c.vdr = inner.vdr;
    c.vqr = inner.vqr;

    c.Evd = inner.Evd;
    c.Evq = inner.Evq;
    c.md  = inner.md;
    c.mq  = inner.mq;

    return c;
}

void CSC_RL::rhs(ODE::real_type t, const ODE::vec_type& x, ODE::vec_type& dxdt) const
{
    // ------------------------------------------------------------
    // Read states
    // ------------------------------------------------------------
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
    dxdt(IGD) = (p_.egd - p_.Rg * igd - ed + c.w_hat * p_.Lg * igq) / p_.Lg;
    dxdt(IGQ) = (p_.egq - p_.Rg * igq - eq - c.w_hat * p_.Lg * igd) / p_.Lg;
    dxdt(ED) = (igd - id + c.w_hat * p_.Cg * eq) / p_.Cg;
    dxdt(EQ) = (igq - iq - c.w_hat * p_.Cg * ed) / p_.Cg;

    // ------------------------------------------------------------
    // CSC AC-side and DC-side dynamics
    // ------------------------------------------------------------
    dxdt(ID) = (ed - vd - p_.Rf * id + c.w_hat * p_.Lf * iq) / p_.Lf;
    dxdt(IQ) = (eq - vq - p_.Rf * iq - c.w_hat * p_.Lf * id) / p_.Lf;
    dxdt(VD) = (id - c.md * istk + c.w_hat * p_.Cf * vq) / p_.Cf;
    dxdt(VQ) = (iq - c.mq * istk - c.w_hat * p_.Cf * vd) / p_.Cf;
    dxdt(ISTK) = (c.md * vd + c.mq * vq - p_.RL * istk) / p_.Ldc;

    // ------------------------------------------------------------
    // PLL states
    // ------------------------------------------------------------
    dxdt(THETA_HAT) = c.w_hat;
    dxdt(XI_PLL)    = c.e_pll;

    // ------------------------------------------------------------
    // Controller integrator states
    // ------------------------------------------------------------
    dxdt(XI_UCD)  = c.Evd;
    dxdt(XI_UCQ)  = c.Evq;
    dxdt(XI_ISD)  = c.Eid;
    dxdt(XI_ISQ)  = c.Eiq;
    dxdt(XI_IDC2) = c.Eidc2;
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

    const ControlOutput c = compute_control(t, x);

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
