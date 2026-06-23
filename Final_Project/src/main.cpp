#include "ODE4.hh"
#include "CSC_RL.hh"
#include "CSC_RL_Case.hh"
#include "Equilibrium.hh"
#include "CSC_RL_Model.hh"


int main()
{
    try
    {
        // ------------------------------------------------------------
        // Simulation settings
        // ------------------------------------------------------------

        CSC_RL_Model model;
        CSC_RL_Parameters p = make_csc_rl_case();

        p.t0 = 0.0;
        p.tf = 1.0;
        p.idc_ref.times  = {0.0, 0.8, 1.6};
        p.idc_ref.values = {250.0, 300.0, 350.0};
        p.Q_ref.times    = {0.0, 0.8, 1.6};
        p.Q_ref.values   = {0.0, 30000.0, -30000.0};

        ODE::SolverOptions sim;
        sim.h = 1e-5;
        sim.save_every = 20;
        sim.output_file = "results/csc_rl.csv";
        
        // ------------------------------------------------------------
        // Calculate equilibrium point
        // ------------------------------------------------------------
        EquilibriumOptions eq_opt;
        eq_opt.verbose = true;
        eq_opt.tolerance = 1e-8;
        eq_opt.max_iter = 80;
        eq_opt.max_sub_iter = 30;
        eq_opt.damp_factor = 0.6;
        eq_opt.newton_verbose = true;

        compute_csc_rl_equilibrium(model, p, eq_opt);

        // ------------------------------------------------------------
        // Build and solve problem
        // ------------------------------------------------------------
        CSC_RL problem(model, p);

        ODE::solve_rk4(problem, sim);

        std::cout << "Simulation completed successfully.\n";
        std::cout << "Output file: " << sim.output_file << "\n";

        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}