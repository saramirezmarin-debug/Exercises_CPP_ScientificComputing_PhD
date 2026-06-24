#include "ODE4.hh"

/**
 * @file ODE4.cc
 * @brief Implementation of the generic RK4 solver and stair-signal utilities.
 *
 * @details
 * This file implements the fixed-step RK4 algorithm, CSV export, algebraic
 * output handling, and the piecewise-constant reference signal evaluation.
 */

// ============================================================
// ODE solver interface and RK4 declarations.
// ============================================================

namespace ODE
{
    namespace
    {
        
        /**
         * @brief Internal allocation-free RK4 step.
         *
         * @details
         * This helper receives all temporary vectors by reference to avoid repeated
         * memory allocation inside the time-integration loop.
         *
         * @param problem ODE problem to integrate.
         * @param t Current simulation time.
         * @param x Current state vector.
         * @param h Integration step.
         * @param x_next Output vector for the next state.
         * @param k1 First RK4 stage.
         * @param k2 Second RK4 stage.
         * @param k3 Third RK4 stage.
         * @param k4 Fourth RK4 stage.
         * @param xtmp Temporary state vector.
         */
        
        void rk4_step_inplace(const ODE_Problem_base& problem,                     // Internal RK4 step implementation
                            real_type t,
                            const vec_type& x,
                            real_type h,
                            vec_type& x_next,
                            vec_type& k1,
                            vec_type& k2,
                            vec_type& k3,
                            vec_type& k4,
                            vec_type& xtmp)
        {
            problem.rhs(t, x, k1);                                                  // Compute k1 = f(t, x)
            xtmp.noalias() = x + 0.5 * h * k1;
            problem.rhs(t + 0.5 * h, xtmp, k2);                                     // Compute k2 = f(t+h/2, x+h*k1/2)
            xtmp.noalias() = x + 0.5 * h * k2;
            problem.rhs(t + 0.5 * h, xtmp, k3);                                     // Compute k3 = f(t+h/2, x+h*k2/2)
            xtmp.noalias() = x + h * k3;
            problem.rhs(t + h, xtmp, k4);                                           // Compute k4 = f(t+h, x+h*k3)
            x_next.noalias() = x + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);     // Compute x_{k+1}
        }
    }

    // Evaluate the piecewise-constant stair signal at time t.
    real_type StairSignal::value(real_type t) const
    {
        if (times.empty())  return constant_value;                  // If no switching times are provided, use the constant value.
        if (t < times.front())  return values.front();              // Before the first switching time, hold the first specified value.
        auto it = std::upper_bound(times.begin(), times.end(), t);  // Find the first switching time strictly greater than t.
        const std::size_t index = static_cast<std::size_t>(std::distance(times.begin(), it) - 1);   // Select the last interval whose switching time is less than or equal to t.

        return values[index];                                       // Return the value associated with the active stair interval.
    }

    // Validate the consistency of the stair signal definition.
    void StairSignal::validate(const std::string& name) const
    {
        // Empty time and value vectors define a constant signal.
        if (times.empty() && values.empty()) return;    
        // Times and values must be either both empty or both non-empty.
        if (times.empty() || values.empty()) throw std::runtime_error(name + ": times and values must both be empty or both be non-empty");
        // Throw a descriptive error if the signal is partially defined.
        if (times.size() != values.size())   throw std::runtime_error(name + ": times and values must have the same size");

        for (std::size_t i = 0; i < times.size(); ++i)  // Check every time/value pair.
        {
            // Switching times must be finite numbers.
            if (!std::isfinite(times[i]))          throw std::runtime_error(name + ": all times must be finite");
            // Signal values must be finite numbers.
            if (!std::isfinite(values[i]))         throw std::runtime_error(name + ": all values must be finite");
            // Switching times must be strictly increasing.
            if (i > 0 && times[i] <= times[i - 1]) throw std::runtime_error(name + ": times must be strictly increasing");
        }
    }

    // RK4 one-step function that allocates its own temporary vectors.
    vec_type rk4_step(const ODE_Problem_base& problem, real_type t, const vec_type& x, real_type h)
    {
        if (h <= 0.0)   throw std::runtime_error("rk4_step: step size must be positive");       // Check that step size is positive
        const integer n = problem.n();                                                          // Get the number of dynamic states.
        if (x.size() != n)  throw std::runtime_error("rk4_step: state vector has wrong size");  // Check that the input state vector has the correct dimension.
        vec_type x_next(n), k1(n), k2(n), k3(n), k4(n), xtmp(n);                                // Allocate temporary vectors for the RK4 stages.
        rk4_step_inplace(problem, t, x, h, x_next, k1, k2, k3, k4, xtmp);                       // Perform RK4 step

        return x_next;                                                                          // Return the next state
    }

    // Solve the initial value problem using fixed-step RK4 and export to CSV.
    void solve_rk4(const ODE_Problem_base& problem, real_type h, const std::string& filename, integer save_every)
    {
        if (h <= 0.0)   throw std::runtime_error("solve_rk4: step size must be positive");          // Check that the integration step is positive.
        if (save_every <= 0)    throw std::runtime_error("solve_rk4: save_every must be positive"); // Check that the output saving interval is positive.
        
        const integer n = problem.n();                          // Get the number of dynamic states.
        const integer n_outputs = problem.n_outputs();          // Get number of algebraic outputs
        const bool has_exact = problem.has_exact_solution();    // Check if exact solution is available
        const real_type t0 = problem.t0();                      // Get initial time
        const real_type tf = problem.tf();                      // Get final time

        if (n <= 0) throw std::runtime_error("solve_rk4: problem dimension must be positive");              // The problem must contain at least one dynamic state.
        if (tf < t0) throw std::runtime_error("solve_rk4: final time must be greater than initial time");   // The final time must not be smaller than the initial time.
        
        vec_type x(n);                                  // State vector
        vec_type x_next(n);                             // Next state vector (for RK4 step)
        vec_type k1(n), k2(n), k3(n), k4(n), xtmp(n);   // Temporary vectors for RK4 stages
        vec_type x_exact(n);                            // Vector for exact solution (if available)
 
        // Initialize the state vector from the problem initial conditions.
        for (integer i = 0; i < n; ++i) 
        {
            x(i) = problem.initial_condition(i);
        }

        std::ofstream csv(filename);                                                    // Open the output CSV file.
        if (!csv)   throw std::runtime_error("solve_rk4: could not open output file");  // Stop if the output file cannot be opened.
        csv << std::setprecision(15);                                                   // Use high numerical precision when writing floating-point values.

        csv << "t";                                 // First CSV column: simulation time.
        
        for (integer i = 0; i < n; ++i)             
        {
            csv << "," << problem.state_name(i);    // Write dynamic state names to the CSV header.
        }
        
        for (integer i = 0; i < n_outputs; ++i)
        {
            csv << "," << problem.output_name(i);   // write algebraic output names if available
        }

        if (has_exact)                              // Add exact solution and error columns when an exact solution is available.
        {
            for (integer i = 0; i < n; ++i)
            {
                csv << "," << problem.state_name(i) << "_exact";    // Write exact solution column names.
            }
            for (integer i = 0; i < n; ++i)
            {
                csv << ",error_x" << i + 1;                         // Error columns
            }
        }
        csv << "\n";                                                // Finish the CSV header row.


        real_type t = t0;           // Initialize simulation time.
        integer step = 0;           // Initialize step counter.

        while (t <= tf + 1e-12)     // Time integration loop with a small tolerance for floating-point roundoff.
        {
            const bool save_this_step = (step % save_every == 0) || (t >= tf);  // Save this sample according to save_every, always including the final time.

            // Write the current simulation sample to the CSV file.
            if (save_this_step)
            {
                csv << t;                           // Write current simulation time.
                for (integer i = 0; i < n; ++i)
                {
                    csv << "," << x(i);             // Write current state values.
                }

                const vec_type y = problem.outputs(t, x);   // Evaluate algebraic outputs at the current time and state.

                if (y.size() != n_outputs)  throw std::runtime_error("solve_rk4: invalid output vector size");  // Check that the output vector has the expected dimension.
                
                for (integer i = 0; i < n_outputs; ++i)
                {
                    csv << "," << y(i);             // Write algebraic output values.
                }

                // If available, write exact solution and absolute error values.
                if (has_exact)
                {
                    for (integer i = 0; i < n; ++i)
                    {
                        x_exact(i) = problem.exact(t, i);
                        csv << "," << x_exact(i);                   // Evaluate and write the exact solution of state i.
                    }

                    for (integer i = 0; i < n; ++i)
                    {
                        csv << "," << std::abs(x(i) - x_exact(i));  // Write the absolute integration error for state i.
                    }
                }

                csv << "\n";                                        // Finish the current CSV row.
            }

            // Stop after writing the final-time sample.
            if (t >= tf)
            {
                break;
            }

            const real_type h_step = std::min(h, tf - t);   // Reduce the last step if needed so that the simulation ends exactly at tf.

            rk4_step_inplace(problem, t, x, h_step, x_next, k1, k2, k3, k4, xtmp);  // Advance the state by one RK4 step.

            x.swap(x_next);     // Swap current and next states without copying vector data.
            t += h_step;        // Advance simulation time 
            ++step;             // Step counter
        }
    }

    // Convenience overload that forwards SolverOptions to the main solver.
    void solve_rk4(const ODE_Problem_base& problem, const SolverOptions& options)
    {
        solve_rk4(problem, options.h, options.output_file, options.save_every);
    }


    // Evaluate all algebraic outputs by calling output(t, x, i).
    vec_type ODE_Problem_base::outputs(real_type t, const vec_type& x) const    
    {
        vec_type y(n_outputs());

        for (integer i = 0; i < n_outputs(); ++i)
        {
            y(i) = output(t, x, i);     // Fill the output vector by evaluating each scalar output.
        }

        return y;
    }
}