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

CSCEquilibriumProblem::CSCEquilibriumProblem(const CSC_RL_Model& model, const CSC_RL_Parameters& p, const EquilibriumReferences& ref)
    : model_(model)
    , p_(p)
    , ref_(ref)
{
    _x0 = make_equilibrium_initial_guess(p_, ref_);
}

CSCEquilibriumProblem::d_dual_vec
CSCEquilibriumProblem::F(d_dual_vec const& X) const
{
    using DualReal = AD::Dual<real>;

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
    // Reduced equilibrium vector -> full state vector
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

    // theta_hat is not part of the equilibrium unknowns
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
    // Evaluate the same model used by ODE
    // ------------------------------------------------------------
    model_.rhs(p_, p_.t0, x, dxdt, input);

    // ------------------------------------------------------------
    // Equilibrium residuals: dx/dt = 0
    // theta_hat is excluded because dtheta/dt = omega_hat
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

bool compute_csc_rl_equilibrium(const CSC_RL_Model& model, CSC_RL_Parameters& p, const EquilibriumOptions& options)
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

    CSCEquilibriumProblem eq_problem(model, p, eq_ref);

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

bool compute_csc_rl_equilibrium(const CSC_RL_Model& model, CSC_RL_Parameters& p, bool verbose)
{
    EquilibriumOptions options;
    options.verbose = verbose;

    return compute_csc_rl_equilibrium(model, p, options);
}

bool compute_csc_rl_equilibrium(CSC_RL_Parameters& p, const EquilibriumOptions& options)
{
    static CSC_RL_Model model;

    return compute_csc_rl_equilibrium(model, p, options);
}
