#pragma once



#include "CSC_RL.hh"
#include "Newton.hh"

struct EquilibriumReferences
{
    ODE::real_type idc_ref = 0.0;
    ODE::real_type Q_ref   = 0.0;
};

struct EquilibriumOptions
{
    bool verbose = false;

    int max_iter = 50;
    int max_sub_iter = 20;
    ODE::real_type damp_factor = 0.8;
    ODE::real_type tolerance = 1e-10;

    bool newton_verbose = false;

    bool use_reference_override = false;
    EquilibriumReferences reference;
};

class CSCEquilibriumProblem : public AD::NewtonProblemAD<ODE::real_type>
{
public:
    using real = ODE::real_type;

    using AD::NewtonProblemAD<real>::F;
    using AD::NewtonProblemAD<real>::JF;

    static constexpr int Z_N = 15;
    CSCEquilibriumProblem(const CSC_RL_Parameters& p, const EquilibriumReferences& ref);

    d_dual_vec F(d_dual_vec const& X) const override;

private:
    CSC_RL_Parameters p_;
    EquilibriumReferences ref_;
};

ODE::vec_type make_equilibrium_initial_guess(const CSC_RL_Parameters& p, const EquilibriumReferences& ref);
ODE::vec_type expand_equilibrium_to_full_state(const ODE::vec_type& z, const CSC_RL_Parameters& p);
void apply_equilibrium_to_parameters(CSC_RL_Parameters& p, const ODE::vec_type& x0_full);
bool compute_csc_rl_equilibrium(CSC_RL_Parameters& p, const EquilibriumOptions& options);
bool compute_csc_rl_equilibrium(CSC_RL_Parameters& p, bool verbose = false);


