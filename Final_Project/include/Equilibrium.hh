#pragma once

#include "ODE4.hh"
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

bool solve_equilibrium_problem(AD::NewtonProblemAD<ODE::real_type>& eq_problem, const EquilibriumOptions& options, ODE::vec_type& z_star);