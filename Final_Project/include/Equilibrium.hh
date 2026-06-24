#pragma once

#include "ODE4.hh"
#include "Newton.hh"

// ============================================================
// Equilibrium solver interface.
//
// This file defines the options and reference containers used
// to solve steady-state problems through a Newton solver with
// automatic differentiation.
// ============================================================

/**
 * @file Equilibrium.hh
 * @brief Generic interface for Newton-based equilibrium calculation.
 *
 * @details
 * This file defines the reference and option structures used to compute a
 * steady-state operating point before running the time-domain simulation.
 *
 * A model-specific equilibrium problem must derive from
 * AD::NewtonProblemAD<ODE::real_type> and define the nonlinear residual:
 *
 * \f[
 * F(z)=0.
 * \f]
 */


 /**
 * @brief References imposed during equilibrium calculation.
 */

struct EquilibriumReferences
{
    ODE::real_type idc_ref = 0.0;
    ODE::real_type Q_ref   = 0.0;
};

/**
 * @brief Options for the damped Newton equilibrium solver.
 */

struct EquilibriumOptions
{
    bool verbose = false;

    int max_iter = 50;                       // Maximum number of Newton iterations.
    int max_sub_iter = 20;                   // Maximum number of damping attempts per Newton iteration.
    ODE::real_type damp_factor = 0.8;        // Damping factor used when a Newton step is rejected.
    ODE::real_type tolerance = 1e-10;        // Residual tolerance for declaring convergence.
    bool newton_verbose = false;             // Enable detailed Newton iteration messages.
    bool use_reference_override = false;     // If true, use manually provided references for the equilibrium point.
    EquilibriumReferences reference;         // Optional manually specified equilibrium references.
};

/**
 * @brief Solves a generic automatic-differentiation Newton equilibrium problem.
 *
 * @param eq_problem Nonlinear problem providing residual and Jacobian.
 * @param options Equilibrium solver options.
 * @param z_star Output vector containing the converged equilibrium solution.
 * @return true if Newton converged, false otherwise.
 */

// Solve a generic AD-based Newton equilibrium problem.
bool solve_equilibrium_problem(AD::NewtonProblemAD<ODE::real_type>& eq_problem, const EquilibriumOptions& options, ODE::vec_type& z_star);