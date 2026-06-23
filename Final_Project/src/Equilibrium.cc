#include "Equilibrium.hh"

#include <iostream>
#include <cmath>
#include <stdexcept>

namespace
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

    template <typename T> T square(T const& x) {return x * x;}
}

CSCEquilibriumProblem::CSCEquilibriumProblem(const CSC_RL_Parameters& p, const EquilibriumReferences& ref)
    : p_(p)
    , ref_(ref)
{
    _x0 = make_equilibrium_initial_guess(p_, ref_);
}

CSCEquilibriumProblem::d_dual_vec
CSCEquilibriumProblem::F(d_dual_vec const& X) const
{
    d_dual_vec res;
    res.resize(Z_N);

    auto const& igd  = X(Z_IGD);
    auto const& igq  = X(Z_IGQ);
    auto const& ed   = X(Z_ED);
    auto const& eq   = X(Z_EQ);
    auto const& id   = X(Z_ID);
    auto const& iq   = X(Z_IQ);
    auto const& vd   = X(Z_VD);
    auto const& vq   = X(Z_VQ);
    auto const& istk = X(Z_ISTK);

    auto const& xi_pll  = X(Z_XI_PLL);
    auto const& xi_ucd  = X(Z_XI_UCD);
    auto const& xi_ucq  = X(Z_XI_UCQ);
    auto const& xi_isd  = X(Z_XI_ISD);
    auto const& xi_isq  = X(Z_XI_ISQ);
    auto const& xi_idc2 = X(Z_XI_IDC2);

    // ------------------------------------------------------------
    // PLL
    // ------------------------------------------------------------
    auto const e_pll = eq / p_.Vdq_nom;
    auto const w_hat = p_.w0_pll + p_.kp_pll * e_pll + p_.ki_pll * xi_pll;

    // ------------------------------------------------------------
    // Outer DC loop
    // ------------------------------------------------------------
    auto const idc_ref = ref_.idc_ref;
    auto const Q_ref   = ref_.Q_ref;

    auto const E_idc2 = square(idc_ref) - square(istk);
    auto const P_ref = p_.kpO * E_idc2 + p_.kiO * xi_idc2;

    // ------------------------------------------------------------
    // AC current references
    // ------------------------------------------------------------
    auto const V2 = square(ed) + square(eq);

    auto const id_ref = (ed * P_ref + eq * Q_ref) / V2;
    auto const iq_ref = (eq * P_ref - ed * Q_ref) / V2;

    auto const E_id = id_ref - id;
    auto const E_iq = iq_ref - iq;

    // ------------------------------------------------------------
    // Current controller
    // ------------------------------------------------------------
    auto const vd_ref = ed + w_hat * p_.Lf * iq - p_.Rf * id - p_.kp2 * E_id - p_.ki2 * xi_isd;
    auto const vq_ref = eq - w_hat * p_.Lf * id - p_.Rf * iq - p_.kp2 * E_iq - p_.ki2 * xi_isq;

    auto const E_vd = vd_ref - vd;
    auto const E_vq = vq_ref - vq;

    // ------------------------------------------------------------
    // Voltage controller
    // ------------------------------------------------------------
    auto const md = (id + w_hat * p_.Cf * vq - p_.kp1 * E_vd - p_.ki1 * xi_ucd) / istk;
    auto const mq = (iq - w_hat * p_.Cf * vd - p_.kp1 * E_vq - p_.ki1 * xi_ucq) / istk;

    // ------------------------------------------------------------
    // Differential equations imposed as zero
    // ------------------------------------------------------------
    res(Z_IGD) = (p_.egd - p_.Rg * igd - ed + w_hat * p_.Lg * igq) / p_.Lg;
    res(Z_IGQ) = (p_.egq - p_.Rg * igq - eq - w_hat * p_.Lg * igd) / p_.Lg;

    res(Z_ED) = (igd - id + w_hat * p_.Cg * eq) / p_.Cg;
    res(Z_EQ) = (igq - iq - w_hat * p_.Cg * ed) / p_.Cg;

    res(Z_ID) = (ed - vd - p_.Rf * id + w_hat * p_.Lf * iq) / p_.Lf;
    res(Z_IQ) = (eq - vq - p_.Rf * iq - w_hat * p_.Lf * id) / p_.Lf;
    res(Z_VD) = (id - md * istk + w_hat * p_.Cf * vq) / p_.Cf;
    res(Z_VQ) = (iq - mq * istk - w_hat * p_.Cf * vd) / p_.Cf;
    res(Z_ISTK) = (md * vd + mq * vq - p_.RL * istk) / p_.Ldc;

    // ------------------------------------------------------------
    // Integrator residuals
    // ------------------------------------------------------------
    res(Z_XI_PLL)  = e_pll;
    res(Z_XI_UCD)  = E_vd;
    res(Z_XI_UCQ)  = E_vq;
    res(Z_XI_ISD)  = E_id;
    res(Z_XI_ISQ)  = E_iq;
    res(Z_XI_IDC2) = E_idc2;

    return res;
}

ODE::vec_type make_equilibrium_initial_guess(const CSC_RL_Parameters& p, const EquilibriumReferences& ref)
{
    ODE::vec_type z0;
    z0.resize(CSCEquilibriumProblem::Z_N);
    z0.setZero();

    const ODE::real_type ed0 = p.Vdq_nom;
    const ODE::real_type eq0 = 0.0;

    const ODE::real_type P0 =
        p.RL * ref.idc_ref * ref.idc_ref;

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

ODE::vec_type expand_equilibrium_to_full_state(const ODE::vec_type& z, const CSC_RL_Parameters& p)
{
    if (z.size() != CSCEquilibriumProblem::Z_N)
    {
        throw std::runtime_error("Equilibrium: invalid z size");
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

void apply_equilibrium_to_parameters(CSC_RL_Parameters& p,const ODE::vec_type& x0_full)
{
    if (x0_full.size() != NSTATES)
    {
        throw std::runtime_error("Equilibrium: invalid full state size");
    }

    p.x0 = x0_full;
}

bool compute_csc_rl_equilibrium(CSC_RL_Parameters& p,
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
        std::cout << "Computing equilibrium point...\n";
        std::cout << "idc_ref = " << eq_ref.idc_ref << " A\n";
        std::cout << "Q_ref   = " << eq_ref.Q_ref << " var\n";
    }

    CSCEquilibriumProblem eq_problem(p, eq_ref);

    AD::NewtonOptions newton_options;
    newton_options.verbose      = options.newton_verbose;
    newton_options.tolerance    = options.tolerance;
    newton_options.max_iter     = options.max_iter;
    newton_options.max_sub_iter = options.max_sub_iter;
    newton_options.damp_factor  = options.damp_factor;

    AD::NewtonSolver<ODE::real_type> eq_solver(&eq_problem, newton_options);

    if (options.verbose)
    {
        std::cout << "Starting Newton solver...\n";
    }

    eq_solver.solve();

    if (options.verbose)
    {
        std::cout << "Newton solver finished.\n";
    }

    if (!eq_solver.converged())
    {
        if (options.verbose)
        {
            std::cerr << "Equilibrium Newton failed.\n";
            std::cerr << "Simulation will continue with manual initial conditions.\n";
        }

        return false;
    }

    if (options.verbose)
    {
        std::cout << "Equilibrium converged.\n";
        std::cout << "z_star = "
                  << eq_solver.solution().transpose()
                  << '\n';
    }

    const ODE::vec_type x0_eq =
        expand_equilibrium_to_full_state(eq_solver.solution(), p);

    apply_equilibrium_to_parameters(p, x0_eq);

    return true;
}

bool compute_csc_rl_equilibrium(CSC_RL_Parameters& p, bool verbose)
{
    EquilibriumOptions options;
    options.verbose = verbose;

    return compute_csc_rl_equilibrium(p, options);
}