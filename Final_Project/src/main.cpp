#include "ODE4.hh"                     // Solver RK4
#include "CSC_RL.hh"                   // CSC-RL ODE problem wrapper.
#include "CSC_RL_Case.hh"              // Factory function for the default CSC-RL simulation case.
#include "Equilibrium.hh"              // Equilibrium solver
#include "CSC_RL_Model.hh"             // CSC-RL model equations and control laws

#include <iostream>
#include <exception>

/**
 * @file main.cpp
 * @brief Main executable for the CSC-RL RK4 simulation.
 *
 * @details
 * This file defines the complete simulation workflow:
 *
 * - create the CSC-RL model,
 * - load the default parameter set,
 * - modify simulation references,
 * - configure the RK4 solver,
 * - compute a Newton equilibrium point,
 * - initialize the system from that equilibrium point,
 * - run the time-domain simulation,
 * - export the results to CSV.
 */


// =============================================================================
// Main simulation program for the Current source converter with resistive laod.
//
// This file:
//   1. Builds the CSC-RL model and its parameter set.
//   2. Configures simulation references and solver options.
//   3. Computes a consistent equilibrium initial condition.
//   4. Runs the RK4 time-domain simulation.
//   5. Exports the results to a CSV file.
// =============================================================================

/**
 * @brief Program entry point.
 *
 * @details
 * The program first computes a consistent steady-state operating point using
 * Newton's method. If Newton converges, the equilibrium point replaces p.x0 and
 * the RK4 simulation starts from this state. This avoids artificial start-up
 * transients caused by inconsistent manual initial conditions.
 */


int main()
{
    try
    {
        // ------------------------------------------------------------
        // Model, parameter set, references, and RK4 solver settings.
        // ------------------------------------------------------------

        CSC_RL_Model model;
        CSC_RL_Parameters p = make_csc_rl_case();  // Build the default parameter set for the CSC-RL case.

        // Simulation time interval.
        p.t0 = 0.0;
        p.tf = 5.0;
        // DC current reference
        p.idc_ref.times  = {0.0, 0.5, 2.0, 3.0, 4.0};
        p.idc_ref.values = {200.0, 250.0, 200.0, 150.0, 100.0};
        // Reactive power reference
        p.Q_ref.times  = {0.0, 1.7, 3.4};
        p.Q_ref.values = {30000.0, -30000.0, 0.0};




        ODE::SolverOptions sim;                   // Create RK4 solver configuration.
        sim.h = 1e-5;                             // Fixed RK4 integration step size [s].    
        sim.save_every = 20;                      // Save one sample every 20 integration steps to reduce CSV size.
        sim.output_file = "results/csc_rl.csv";   // CSV file where simulation results will be written.
        
        // ------------------------------------------------------------
        // Calculate equilibrium point
        // ------------------------------------------------------------
        EquilibriumOptions eq_opt;                // Equilibrium solver options.
        eq_opt.tolerance = 1e-8;                  // Residual tolerance for Newton convergence.
        eq_opt.max_iter = 80;                     // Maximum number of Newton iterations.
        eq_opt.max_sub_iter = 30;                 // Maximum number of damping attempts per Newton iteration.
        eq_opt.damp_factor = 0.6;                 // Step reduction factor used by the damped Newton method.
        eq_opt.verbose = true;
        eq_opt.newton_verbose = true;


        const bool eq_ok = compute_csc_rl_equilibrium(model, p, eq_opt);  // Compute and apply a consistent steady-state initial condition.
        if (!eq_ok)
        {
            std::cerr << "\nWarning: equilibrium was not found. Using manual initial conditions.\n\n";
        }

        // ------------------------------------------------------------
        // Build and solve problem
        // ------------------------------------------------------------
        CSC_RL problem(model, p);                  // Wrap the CSC-RL model as a generic ODE problem.
        ODE::solve_rk4(problem, sim);              // Solve the initial value problem using fixed-step RK4.
        std::cout << "\nSimulation completed successfully.\n";
        std::cout << "Output file: " << sim.output_file << "\n";

        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}