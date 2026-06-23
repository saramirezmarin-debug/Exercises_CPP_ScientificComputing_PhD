#pragma once

#include "CSC_RL.hh"
#include "Equilibrium.hh"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdexcept>

class CSC_RL_Model
{
public:

    template <typename T>
    struct Input
    {
        T idc_ref;
        T Q_ref;
    };

    template <typename T>
    struct Control
    {
        T e_pll;
        T w_hat;

        T idc_ref;
        T Qref;
        T Pref;

        T idr;
        T iqr;

        T Eid;
        T Eiq;
        T Eidc2;

        T vdr;
        T vqr;

        T Evd;
        T Evq;

        T md;
        T mq;
    };

    template <typename T>
    Input<T> input_from_parameters(const CSC_RL_Parameters& p,
                                   ODE::real_type t) const
    {
        Input<T> input;

        input.idc_ref = T(p.idc_ref.value(t));
        input.Q_ref   = T(p.Q_ref.value(t));

        return input;
    }

    template <typename Vec>
    Control<typename Vec::Scalar>
    control(const CSC_RL_Parameters& p,
            const Vec& x,
            const Input<typename Vec::Scalar>& input) const
    {
        using T = typename Vec::Scalar;

        Control<T> c;

        const T ed   = x(ED);
        const T eq   = x(EQ);
        const T id   = x(ID);
        const T iq   = x(IQ);
        const T vd   = x(VD);
        const T vq   = x(VQ);
        const T istk = x(ISTK);

        const T xi_pll  = x(XI_PLL);
        const T xi_ucd  = x(XI_UCD);
        const T xi_ucq  = x(XI_UCQ);
        const T xi_isd  = x(XI_ISD);
        const T xi_isq  = x(XI_ISQ);
        const T xi_idc2 = x(XI_IDC2);

        // ------------------------------------------------------------
        // PLL
        // ------------------------------------------------------------
        c.e_pll = eq / p.Vdq_nom;
        c.w_hat = p.w0_pll + p.kp_pll * c.e_pll + p.ki_pll * xi_pll;

        // ------------------------------------------------------------
        // References
        // ------------------------------------------------------------
        c.idc_ref = input.idc_ref;
        c.Qref    = input.Q_ref;

        // ------------------------------------------------------------
        // Outer loop
        // ------------------------------------------------------------
        c.Eidc2 = c.idc_ref * c.idc_ref - istk * istk;
        c.Pref  = p.kpO * c.Eidc2 + p.kiO * xi_idc2;

        // ------------------------------------------------------------
        // AC current references
        // ------------------------------------------------------------
        const T V2 = ed * ed + eq * eq;

        c.idr = (ed * c.Pref + eq * c.Qref) / V2;
        c.iqr = (eq * c.Pref - ed * c.Qref) / V2;

        c.Eid = c.idr - id;
        c.Eiq = c.iqr - iq;

        // ------------------------------------------------------------
        // Current controller
        // ------------------------------------------------------------
        c.vdr = ed + c.w_hat * p.Lf * iq - p.Rf * id - p.kp2 * c.Eid - p.ki2 * xi_isd;
        c.vqr = eq - c.w_hat * p.Lf * id - p.Rf * iq - p.kp2 * c.Eiq - p.ki2 * xi_isq;
        c.Evd = c.vdr - vd;
        c.Evq = c.vqr - vq;

        // ------------------------------------------------------------
        // Voltage controller
        // ------------------------------------------------------------
        c.md = (id + c.w_hat * p.Cf * vq - p.kp1 * c.Evd - p.ki1 * xi_ucd) / istk;
        c.mq = (iq - c.w_hat * p.Cf * vd - p.kp1 * c.Evq - p.ki1 * xi_ucq) / istk;

        return c;
    }

    template <typename Vec>
    void rhs(const CSC_RL_Parameters& p,
             ODE::real_type t,
             const Vec& x,
             Vec& dxdt,
             const Input<typename Vec::Scalar>& input) const
    {
        using T = typename Vec::Scalar;

        dxdt.resize(NSTATES);

        const T igd  = x(IGD);
        const T igq  = x(IGQ);
        const T ed   = x(ED);
        const T eq   = x(EQ);

        const T id   = x(ID);
        const T iq   = x(IQ);
        const T vd   = x(VD);
        const T vq   = x(VQ);
        const T istk = x(ISTK);

        const Control<T> c = control(p, x, input);

        // ------------------------------------------------------------
        // PCC dynamics
        // ------------------------------------------------------------
        dxdt(IGD) = (p.egd - p.Rg * igd - ed + c.w_hat * p.Lg * igq) / p.Lg;
        dxdt(IGQ) = (p.egq - p.Rg * igq - eq - c.w_hat * p.Lg * igd) / p.Lg;

        dxdt(ED) = (igd - id + c.w_hat * p.Cg * eq) / p.Cg;
        dxdt(EQ) = (igq - iq - c.w_hat * p.Cg * ed) / p.Cg;

        // ------------------------------------------------------------
        // CSC dynamics
        // ------------------------------------------------------------
        dxdt(ID) = (ed - vd - p.Rf * id + c.w_hat * p.Lf * iq) / p.Lf;
        dxdt(IQ) = (eq - vq - p.Rf * iq - c.w_hat * p.Lf * id) / p.Lf;

        dxdt(VD) = (id - c.md * istk + c.w_hat * p.Cf * vq) / p.Cf;
        dxdt(VQ) = (iq - c.mq * istk - c.w_hat * p.Cf * vd) / p.Cf;

        dxdt(ISTK) = (c.md * vd + c.mq * vq - p.RL * istk) / p.Ldc;

        // ------------------------------------------------------------
        // PLL
        // ------------------------------------------------------------
        dxdt(THETA_HAT) = c.w_hat;
        dxdt(XI_PLL)    = c.e_pll;

        // ------------------------------------------------------------
        // Integrators
        // ------------------------------------------------------------
        dxdt(XI_UCD)  = c.Evd;
        dxdt(XI_UCQ)  = c.Evq;
        dxdt(XI_ISD)  = c.Eid;
        dxdt(XI_ISQ)  = c.Eiq;
        dxdt(XI_IDC2) = c.Eidc2;
    }

    template <typename Vec>
    void rhs(const CSC_RL_Parameters& p,
             ODE::real_type t,
             const Vec& x,
             Vec& dxdt) const
    {
        const auto input = input_from_parameters<typename Vec::Scalar>(p, t);
        rhs(p, t, x, dxdt, input);
    }
};

namespace CSC_RL_Equilibrium_Detail
{
    enum EquilibriumIndex
    {
        Z_IGD = 0,
        Z_IGQ,
        Z_ED,
        Z_EQ,
        Z_ID,
        Z_IQ,
        Z_VD,
        Z_VQ,
        Z_ISTK,
        Z_XI_PLL,
        Z_XI_UCD,
        Z_XI_UCQ,
        Z_XI_ISD,
        Z_XI_ISQ,
        Z_XI_IDC2,
        Z_N
    };
}

inline ODE::vec_type make_csc_rl_equilibrium_initial_guess(const CSC_RL_Parameters& p,
                                                           const EquilibriumReferences& ref)
{
    using namespace CSC_RL_Equilibrium_Detail;

    ODE::vec_type z0;
    z0.resize(Z_N);
    z0.setZero();

    const ODE::real_type ed0 = p.Vdq_nom;
    const ODE::real_type eq0 = 0.0;

    const ODE::real_type P0 = p.RL * ref.idc_ref * ref.idc_ref;

    const ODE::real_type id0 = P0 / ed0;
    const ODE::real_type iq0 = -ref.Q_ref / ed0;

    z0(Z_IGD) = id0;
    z0(Z_IGQ) = iq0;

    z0(Z_ED) = ed0;
    z0(Z_EQ) = eq0;

    z0(Z_ID) = id0;
    z0(Z_IQ) = iq0;

    z0(Z_VD) = ed0;
    z0(Z_VQ) = eq0;

    z0(Z_ISTK) = ref.idc_ref;

    z0(Z_XI_PLL)  = 0.0;
    z0(Z_XI_UCD)  = 0.0;
    z0(Z_XI_UCQ)  = 0.0;
    z0(Z_XI_ISD)  = 0.0;
    z0(Z_XI_ISQ)  = 0.0;
    z0(Z_XI_IDC2) = 0.0;

    return z0;
}

class CSC_RL_EquilibriumProblem : public AD::NewtonProblemAD<ODE::real_type>
{
public:
    using real = ODE::real_type;

    using AD::NewtonProblemAD<real>::F;
    using AD::NewtonProblemAD<real>::JF;

    static constexpr int Z_N = CSC_RL_Equilibrium_Detail::Z_N;

    CSC_RL_EquilibriumProblem(const CSC_RL_Model& model,
                              const CSC_RL_Parameters& p,
                              const EquilibriumReferences& ref)
        : model_(model)
        , p_(p)
        , ref_(ref)
    {
        this->_x0 = make_csc_rl_equilibrium_initial_guess(p_, ref_);
    }

    d_dual_vec F(d_dual_vec const& X) const override
    {
        using DualReal = AD::Dual<real>;
        using namespace CSC_RL_Equilibrium_Detail;

        d_dual_vec x;
        d_dual_vec dxdt;
        d_dual_vec res;

        x.resize(NSTATES);
        dxdt.resize(NSTATES);
        res.resize(Z_N);

        x.setZero();
        dxdt.setZero();
        res.setZero();

        // ------------------------------------------------------------
        // Reduced equilibrium vector -> full CSC_RL state vector
        // ------------------------------------------------------------
        x(IGD)  = X(Z_IGD);
        x(IGQ)  = X(Z_IGQ);
        x(ED)   = X(Z_ED);
        x(EQ)   = X(Z_EQ);

        x(ID)   = X(Z_ID);
        x(IQ)   = X(Z_IQ);
        x(VD)   = X(Z_VD);
        x(VQ)   = X(Z_VQ);
        x(ISTK) = X(Z_ISTK);

        // theta_hat is not part of the equilibrium unknowns.
        x(THETA_HAT) = DualReal(p_.x0(THETA_HAT));

        x(XI_PLL)  = X(Z_XI_PLL);
        x(XI_UCD)  = X(Z_XI_UCD);
        x(XI_UCQ)  = X(Z_XI_UCQ);
        x(XI_ISD)  = X(Z_XI_ISD);
        x(XI_ISQ)  = X(Z_XI_ISQ);
        x(XI_IDC2) = X(Z_XI_IDC2);

        // ------------------------------------------------------------
        // Fixed references for equilibrium
        // ------------------------------------------------------------
        typename CSC_RL_Model::template Input<DualReal> input;

        input.idc_ref = DualReal(ref_.idc_ref);
        input.Q_ref   = DualReal(ref_.Q_ref);

        // ------------------------------------------------------------
        // Evaluate the same model used by ODE4
        // ------------------------------------------------------------
        model_.rhs(p_, p_.t0, x, dxdt, input);

        // ------------------------------------------------------------
        // Equilibrium residuals: dx/dt = 0
        // theta_hat is excluded because dtheta_hat/dt = omega_hat.
        // ------------------------------------------------------------
        res(Z_IGD) = dxdt(IGD);
        res(Z_IGQ) = dxdt(IGQ);
        res(Z_ED)  = dxdt(ED);
        res(Z_EQ)  = dxdt(EQ);

        res(Z_ID)   = dxdt(ID);
        res(Z_IQ)   = dxdt(IQ);
        res(Z_VD)   = dxdt(VD);
        res(Z_VQ)   = dxdt(VQ);
        res(Z_ISTK) = dxdt(ISTK);

        res(Z_XI_PLL)  = dxdt(XI_PLL);
        res(Z_XI_UCD)  = dxdt(XI_UCD);
        res(Z_XI_UCQ)  = dxdt(XI_UCQ);
        res(Z_XI_ISD)  = dxdt(XI_ISD);
        res(Z_XI_ISQ)  = dxdt(XI_ISQ);
        res(Z_XI_IDC2) = dxdt(XI_IDC2);

        return res;
    }

private:
    const CSC_RL_Model& model_;
    CSC_RL_Parameters p_;
    EquilibriumReferences ref_;
};

inline ODE::vec_type expand_csc_rl_equilibrium_to_full_state(const ODE::vec_type& z,
                                                             const CSC_RL_Parameters& p)
{
    using namespace CSC_RL_Equilibrium_Detail;

    if (z.size() != Z_N)
    {
        throw std::runtime_error("CSC_RL equilibrium: invalid z size");
    }

    ODE::vec_type x0;
    x0.resize(NSTATES);
    x0.setZero();

    x0(IGD) = z(Z_IGD);
    x0(IGQ) = z(Z_IGQ);

    x0(ED) = z(Z_ED);
    x0(EQ) = z(Z_EQ);

    x0(ID) = z(Z_ID);
    x0(IQ) = z(Z_IQ);

    x0(VD) = z(Z_VD);
    x0(VQ) = z(Z_VQ);

    x0(ISTK) = z(Z_ISTK);

    x0(THETA_HAT) = p.x0(THETA_HAT);

    x0(XI_PLL)  = z(Z_XI_PLL);
    x0(XI_UCD)  = z(Z_XI_UCD);
    x0(XI_UCQ)  = z(Z_XI_UCQ);
    x0(XI_ISD)  = z(Z_XI_ISD);
    x0(XI_ISQ)  = z(Z_XI_ISQ);
    x0(XI_IDC2) = z(Z_XI_IDC2);

    return x0;
}

inline void apply_csc_rl_equilibrium_to_parameters(CSC_RL_Parameters& p,
                                                   const ODE::vec_type& x0_full)
{
    if (x0_full.size() != NSTATES)
    {
        throw std::runtime_error("CSC_RL equilibrium: invalid full state size");
    }

    p.x0 = x0_full;
}

inline bool compute_csc_rl_equilibrium(const CSC_RL_Model& model,
                                       CSC_RL_Parameters& p,
                                       const EquilibriumOptions& options)
{
    EquilibriumReferences eq_ref;

    if (options.use_reference_override)
    {
        eq_ref = options.reference;
    }
    else
    {
        eq_ref.idc_ref = p.idc_ref.value(p.t0);
        eq_ref.Q_ref   = p.Q_ref.value(p.t0);
    }

    if (options.verbose)
    {
        std::cout << "Computing CSC_RL equilibrium point...\n";
        std::cout << "idc_ref = " << eq_ref.idc_ref << " A\n";
        std::cout << "Q_ref   = " << eq_ref.Q_ref << " var\n";
    }

    CSC_RL_EquilibriumProblem eq_problem(model, p, eq_ref);

    ODE::vec_type z_star;
    const bool ok = solve_equilibrium_problem(eq_problem, options, z_star);

    if (!ok)
    {
        return false;
    }

    const ODE::vec_type x0_eq =
        expand_csc_rl_equilibrium_to_full_state(z_star, p);

    apply_csc_rl_equilibrium_to_parameters(p, x0_eq);

    return true;
}

inline bool compute_csc_rl_equilibrium(const CSC_RL_Model& model,
                                       CSC_RL_Parameters& p,
                                       bool verbose = false)
{
    EquilibriumOptions options;
    options.verbose = verbose;

    return compute_csc_rl_equilibrium(model, p, options);
}

inline bool compute_csc_rl_equilibrium(CSC_RL_Parameters& p,
                                       const EquilibriumOptions& options)
{
    static CSC_RL_Model model;

    return compute_csc_rl_equilibrium(model, p, options);
}

inline bool compute_csc_rl_equilibrium(CSC_RL_Parameters& p,
                                       bool verbose = false)
{
    EquilibriumOptions options;
    options.verbose = verbose;

    return compute_csc_rl_equilibrium(p, options);
}