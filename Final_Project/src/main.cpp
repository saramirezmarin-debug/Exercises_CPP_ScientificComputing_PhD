#include "ODE4.hh"
#include "CSC_RL.hh"
#include "CSC_RL_Case.hh"
#include "Equilibrium.hh"
#include "SimulationConfig.hh"



int main()
{
    try
    {
        // ------------------------------------------------------------
        // Simulation settings
        // ------------------------------------------------------------

        CSC_RL_Parameters p = make_csc_rl_parameters();

        p.t0 = 0.0;
        p.tf = 5;

        configure_references(p);

        SimulationConfig sim;
        sim.h = 1e-5;
        sim.save_every = 20;
        sim.output_file = "results/csc_rl.csv";
        
        // ------------------------------------------------------------
        // Calculate equilibrium point
        // ------------------------------------------------------------

        EquilibriumReferences eq_ref;

        eq_ref.idc_ref = p.idc_ref.value(p.t0);
        eq_ref.Q_ref   = p.Q_ref.value(p.t0);

        std::cout << "Computing equilibrium point...\n";
        std::cout << "idc_ref = " << eq_ref.idc_ref << " A\n";
        std::cout << "Q_ref   = " << eq_ref.Q_ref << " var\n";

        CSCEquilibriumProblem eq_problem(p, eq_ref);

        AD::NewtonOptions newton_options;
        newton_options.verbose = true;
        newton_options.tolerance = 1e-10;
        newton_options.max_iter = 50;
        newton_options.max_sub_iter = 20;
        newton_options.damp_factor = 0.8;
        AD::NewtonSolver<ODE::real_type> eq_solver(&eq_problem, newton_options);

        std::cout << "Starting Newton solver...\n";
        eq_solver.solve();
        std::cout << "Newton solver finished.\n";

        if (eq_solver.converged())
        {
            std::cout << "Equilibrium converged.\n";
            std::cout << "z_star = "
                      << eq_solver.solution().transpose()
                      << '\n';

            const ODE::vec_type x0_eq =
                expand_equilibrium_to_full_state(eq_solver.solution(), p);

            apply_equilibrium_to_parameters(p, x0_eq);
        }
        else
        {
            std::cerr << "Equilibrium Newton failed.\n";
            std::cerr << "Simulation will continue with manual initial conditions.\n";
        }

        // ------------------------------------------------------------
        // Build and solve problem
        // ------------------------------------------------------------
        CSC_RL problem(p);

        ODE::solve_rk4(problem, sim.h, "results/csc_rl.csv", sim.save_every);

        std::cout << "Simulation completed successfully.\n";
        std::cout << "Output file: results/csc_rl.csv\n";

        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}