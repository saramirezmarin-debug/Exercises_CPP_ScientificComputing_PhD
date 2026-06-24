#include "Equilibrium.hh"

#include <iostream>
#include <iomanip>


// ============================================================
// Equilibrium solver implementation.
// This file wraps the Newton solver.
// ============================================================

/**
 * @file Equilibrium.cc
 * @brief Implementation of the generic equilibrium solver wrapper.
 *
 * @details
 * This file translates EquilibriumOptions into AD::NewtonOptions, runs the
 * damped Newton solver, checks convergence, and returns the converged
 * equilibrium vector.
 */

 /**
 * @brief Runs the Newton solver for a generic equilibrium problem.
 *
 * @details
 * The function returns false if the nonlinear residual does not reach the
 * requested tolerance. In that case, the simulation can continue using manual
 * initial conditions.
 *
 * @param eq_problem Nonlinear AD-based equilibrium problem.
 * @param options Solver configuration.
 * @param z_star Converged solution if Newton succeeds.
 * @return true if the equilibrium was found.
 */

bool solve_equilibrium_problem(
    AD::NewtonProblemAD<ODE::real_type>& eq_problem, // Newton problem
    const EquilibriumOptions& options,               // Solver configuration
    ODE::vec_type& z_star)                           // Find solution
{
    // Translate equilibrium options into Newton solver options.
    AD::NewtonOptions newton_options;
    newton_options.verbose      = options.newton_verbose;
    newton_options.tolerance    = options.tolerance;
    newton_options.max_iter     = options.max_iter;
    newton_options.max_sub_iter = options.max_sub_iter;
    newton_options.damp_factor  = options.damp_factor;

    AD::NewtonSolver<ODE::real_type> eq_solver(&eq_problem, newton_options);


    // Run the Newton iterations.
    if (options.verbose)    std::cout << "Starting Newton solver...\n";
    eq_solver.solve();
    if (options.verbose)    std::cout << "Newton solver finished.\n\n";

    // Check whether Newton reached the requested tolerance.
    if (!eq_solver.converged())
    {
        if (options.verbose)
        {
            std::cerr << "Equilibrium Newton failed.\n";
            std::cerr << "Simulation will continue with manual initial conditions.\n";
        }

        return false;
    }

    z_star = eq_solver.solution();      // Store the converged equilibrium vector.

    const char* z_names[] =
    {
        "igd",
        "igq",
        "ed",
        "eq",
        "id",
        "iq",
        "vd",
        "vq",
        "istk",
        "xi_pll",
        "xi_ucd",
        "xi_ucq",
        "xi_isd",
        "xi_isq",
        "xi_idc2"
    };


    // Print the converged equilibrium point 
    if (options.verbose)
    {
        std::cout << "Equilibrium converged.\n";
        std::cout << "z_star:\n";

        for (Eigen::Index i = 0; i < z_star.size(); ++i)
        {
            std::cout << "  "
                    << std::setw(8) << z_names[i]
                    << " = "
                    << z_star(i)
                    << '\n';
        }
    }

    return true;
}