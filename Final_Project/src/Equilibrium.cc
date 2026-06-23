#include "Equilibrium.hh"

#include <iostream>

bool solve_equilibrium_problem(AD::NewtonProblemAD<ODE::real_type>& eq_problem, const EquilibriumOptions& options, ODE::vec_type& z_star)
{
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

    z_star = eq_solver.solution();

    if (options.verbose)
    {
        std::cout << "Equilibrium converged.\n";
        std::cout << "z_star = " << z_star.transpose() << '\n';
    }

    return true;
}