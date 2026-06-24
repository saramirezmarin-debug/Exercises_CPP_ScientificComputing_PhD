#pragma once

/**
 * @file ODE4.hh
 * @brief Generic ODE interface and fixed-step fourth-order Runge-Kutta solver.
 *
 * @details
 * This file defines the generic interface used by the RK4 solver. Any dynamic
 * system must inherit from ODE::ODE_Problem_base and implement the state
 * dimension, initial condition, time limits, and right-hand side:
 *
 * \f[
 * \frac{dx}{dt} = f(t,x).
 * \f]
 *
 */

// Eigen
#include <Eigen/Dense>

// Standard
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <iomanip> 

// ============================================================
// ODE solver interface and RK4 declarations.
// ============================================================

/**
 * @brief Namespace containing ODE types, signal utilities, and RK4 solvers.
 */

namespace ODE
{
    
    using real_type = double;                                                   // Floating point type
    using integer = int;                                                        // Integer type
    using mat_type = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;  // Dynamic dense matrix
    using vec_type = Eigen::Matrix<real_type, Eigen::Dynamic, 1>;               // Dynamic column vector

    /**
     * @brief Configuration options for the RK4 solver.
     *
     * @details
     * These options define the numerical integration step, the output sampling
     * frequency, and the CSV file used to store simulation results.
     */

    // Configuration options for the fixed-step RK4 solver.
    struct SolverOptions
    {
        real_type h = 1e-5;                                 // Fixed integration step size [s].
        integer save_every = 1;                             // Save one sample every save_every integration steps.
        std::string output_file = "results/output.csv";     // Output CSV file path.
    };

    /**
     * @brief Piecewise-constant stair signal.
     *
     * @details
     * This structure is used to represent external references such as the DC
     * current reference and the reactive-power reference. Given a list of switching
     * times \f$t_k\f$ and values \f$v_k\f$, the signal holds the last active value:
     *
     * \f[
     * s(t)=v_k,
     * \qquad
     * t_k \leq t < t_{k+1}.
     * \f]
     *
     * If no profile is provided, constant_value is returned.
     */

    // Stair signal used for references or inputs.
    struct StairSignal
    {
        real_type constant_value = 0.0;                 // Value used when no time/value profile is provided.
        std::vector<real_type> times;                   // Switching times of the stair signal.
        std::vector<real_type> values;                  // Signal values associated with each switching time.
        real_type value(real_type t) const;             // Return the signal value at time t.
        void validate(const std::string& name) const;   // Check that the signal definition is consistent.
    };

    /**
     * @brief Abstract base class for initial-value ODE problems.
     *
     * @details
     * A derived class must provide the state dimension, initial condition, time
     * interval, and right-hand side of the dynamic system:
     *
     * \f[
     * \frac{dx}{dt}=f(t,x),
     * \qquad
     * x(t_0)=x_0.
     * \f]
     *
     */
    class ODE_Problem_base
    {
    public:

        ODE_Problem_base() = default;                               // Default constructor.
        virtual ~ODE_Problem_base() = default;                      // Virtual destructor.
        virtual integer n() const = 0;                              // Return the number of dynamic states.
        virtual real_type initial_condition(integer i) const = 0;   // Return the initial condition of state i.

        virtual std::string state_name([[maybe_unused]] integer i) const
        {
            return "x" + std::to_string(i + 1);                     // Return the name of state i for CSV headers.
        }

        virtual real_type t0() const = 0;       // Return the initial simulation time.
        virtual real_type tf() const = 0;       // Return the final simulation time.


        virtual void rhs(real_type t, const vec_type& x, vec_type& dxdt) const = 0;                 // Compute the right-hand side dxdt = f(t, x).

        // Return the exact solution of state i at time t, if available.
        virtual real_type exact([[maybe_unused]] real_type t, [[maybe_unused]] integer i) const     
        {
            return 0;
        }

        // Indicate whether an exact analytical solution is available.
        virtual bool has_exact_solution() const 
        {
            return false;
        }

        // ----------------------------------------------------------------
        // Algebraic outputs 
        // ----------------------------------------------------------------

        // Return the number of algebraic outputs to be exported.
        virtual integer n_outputs() const
        {
            return 0;
        }

        // Return the name of algebraic output i for CSV headers.
        virtual std::string output_name([[maybe_unused]] integer i) const
        {
            return "";
        }

        // Return algebraic output i evaluated at time t and state x.
        virtual real_type output([[maybe_unused]] real_type t, [[maybe_unused]] const vec_type& x, [[maybe_unused]] integer i) const
        {
            return 0.0;
        }

        // Return all algebraic outputs evaluated at time t and state x.
        virtual vec_type outputs(real_type t, const vec_type& x) const;
    };

    /**
     * @brief Performs one classical fourth-order Runge-Kutta step.
     *
     * @details
     * Given the current state \f$x_k\f$ at time \f$t_k\f$, the RK4 stages are
     *
     * \f[
     * k_1=f(t_k,x_k),
     * \f]
     *
     * \f[
     * k_2=f\left(t_k+\frac{h}{2},x_k+\frac{h}{2}k_1\right),
     * \f]
     *
     * \f[
     * k_3=f\left(t_k+\frac{h}{2},x_k+\frac{h}{2}k_2\right),
     * \f]
     *
     * \f[
     * k_4=f(t_k+h,x_k+h k_3).
     * \f]
     *
     * The next state is
     *
     * \f[
     * x_{k+1}
     * =
     * x_k
     * +
     * \frac{h}{6}
     * \left(k_1+2k_2+2k_3+k_4\right).
     * \f]
     *
     * @param problem ODE problem to integrate.
     * @param t Current simulation time.
     * @param x Current state vector.
     * @param h Fixed integration step.
     * @return State vector after one RK4 step.
     */
    vec_type rk4_step(const ODE_Problem_base& problem, real_type t, const vec_type& x, real_type h);

    /**
     * @brief Solves an initial-value problem using fixed-step RK4 and exports CSV results.
     *
     * @param problem ODE problem to solve.
     * @param h Fixed integration step size.
     * @param filename Output CSV file path.
     * @param save_every Saves one sample every @p save_every RK4 steps.
     */
    // Solve the full initial value problem and write results to a CSV file.
    void solve_rk4(const ODE_Problem_base& problem, real_type h, const std::string& filename, integer save_every = 1);

    /**
     * @brief Solves an initial-value problem using a SolverOptions structure.
     *
     * @param problem ODE problem to solve.
     * @param options RK4 solver configuration.
     */
    // Convenience overload using a SolverOptions structure.
    void solve_rk4(const ODE_Problem_base& problem, const SolverOptions& options);
}
