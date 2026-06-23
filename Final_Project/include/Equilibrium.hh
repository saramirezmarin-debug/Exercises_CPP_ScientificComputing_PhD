#pragma once

#ifndef EQUILIBRIUM_HH
#define EQUILIBRIUM_HH

#include "CSC_RL.hh"
#include "Newton.hh"

struct EquilibriumReferences
{
    ODE::real_type idc_ref = 0.0;
    ODE::real_type Q_ref   = 0.0;
};

class CSCEquilibriumProblem : public AD::NewtonProblemAD<ODE::real_type>
{
public:
    using real = ODE::real_type;

    using AD::NewtonProblemAD<real>::F;
    using AD::NewtonProblemAD<real>::JF;

    static constexpr int Z_N = 15;

    CSCEquilibriumProblem(const CSC_RL_Parameters& p,
                          const EquilibriumReferences& ref);

    d_dual_vec F(d_dual_vec const& X) const override;

private:
    CSC_RL_Parameters p_;
    EquilibriumReferences ref_;
};

ODE::vec_type make_equilibrium_initial_guess(const CSC_RL_Parameters& p,
                                             const EquilibriumReferences& ref);

ODE::vec_type expand_equilibrium_to_full_state(const ODE::vec_type& z,
                                               const CSC_RL_Parameters& p);

void apply_equilibrium_to_parameters(CSC_RL_Parameters& p,
                                     const ODE::vec_type& x0_full);

#endif