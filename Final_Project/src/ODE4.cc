#include "ODE4.hh"

namespace ODE
{
    namespace
    {
        
        void rk4_step_inplace(const ODE_Problem_base& problem,
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
            problem.rhs(t, x, k1); // Compute k1 = f(t, x)
            xtmp.noalias() = x + 0.5 * h * k1;
            problem.rhs(t + 0.5 * h, xtmp, k2); // Compute k2 = f(t+h/2, x+h*k1/2)

            xtmp.noalias() = x + 0.5 * h * k2;
            problem.rhs(t + 0.5 * h, xtmp, k3); // Compute k3 = f(t+h/2, x+h*k2/2)

            xtmp.noalias() = x + h * k3;
            problem.rhs(t + h, xtmp, k4); // Compute k4 = f(t+h, x+h*k3)

            x_next.noalias() = x + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4); // Compute x_{k+1}
        }
    }

    real_type StairSignal::value(real_type t) const
    {
        if (times.empty())
        {
            return constant_value;
        }

        if (t < times.front())
        {
            return values.front();
        }

        auto it = std::upper_bound(times.begin(), times.end(), t);
        const std::size_t index =
            static_cast<std::size_t>(std::distance(times.begin(), it) - 1);

        return values[index];
    }

    void StairSignal::validate(const std::string& name) const
    {
        if (times.empty() && values.empty()) return;
        if (times.empty() || values.empty()) throw std::runtime_error(name + ": times and values must both be empty or both be non-empty");
        if (times.size() != values.size())   throw std::runtime_error(name + ": times and values must have the same size");

        for (std::size_t i = 0; i < times.size(); ++i)
        {
            if (!std::isfinite(times[i]))          throw std::runtime_error(name + ": all times must be finite");
            if (!std::isfinite(values[i]))         throw std::runtime_error(name + ": all values must be finite");
            if (i > 0 && times[i] <= times[i - 1]) throw std::runtime_error(name + ": times must be strictly increasing");
        }
    }

    vec_type rk4_step(const ODE_Problem_base& problem, real_type t, const vec_type& x, real_type h)
    {

        if (h <= 0.0) // Check that step size is positive
        {
            throw std::runtime_error("rk4_step: step size must be positive");
        }

        const integer n = problem.n(); // Get problem dimension

        if (x.size() != n) // Check that input state vector has correct size
        {
            throw std::runtime_error("rk4_step: state vector has wrong size");
        }

        vec_type x_next(n), k1(n), k2(n), k3(n), k4(n), xtmp(n); // Temporary vectors for RK4 stages

        rk4_step_inplace(problem, t, x, h, x_next, k1, k2, k3, k4, xtmp); // Perform RK4 step

        return x_next; // Return the next state
    }

    void solve_rk4(const ODE_Problem_base& problem,
                   real_type h,
                   const std::string& filename,
                   integer save_every)
    {
        if (h <= 0.0) // Check that step size is positive
        {
            throw std::runtime_error("solve_rk4: step size must be positive");
        }

        if (save_every <= 0)
        {
            throw std::runtime_error("solve_rk4: save_every must be positive");
        }

        const integer n = problem.n(); // Get problem dimension
        const integer n_outputs = problem.n_outputs(); // Get number of algebraic outputs
        const bool has_exact = problem.has_exact_solution(); // Check if exact solution is available
        const real_type t0 = problem.t0(); // Get initial time
        const real_type tf = problem.tf(); // Get final time

        if (n <= 0) // Check that problem dimension is positive
        {
            throw std::runtime_error("solve_rk4: problem dimension must be positive");
        }

        if (tf < t0) // Check that final time is greater than initial time
        {
            throw std::runtime_error("solve_rk4: final time must be greater than initial time");
        }

        vec_type x(n); // State vector
        vec_type x_next(n); // Next state vector (for RK4 step)

        vec_type k1(n), k2(n), k3(n), k4(n), xtmp(n); // Temporary vectors for RK4 stages
        vec_type x_exact(n); // Vector for exact solution (if available)
 
        for (integer i = 0; i < n; ++i) // Initialize state vector with initial conditions
        {
            x(i) = problem.initial_condition(i);
        }

        // CSV file for writing results
        std::ofstream csv(filename); 
        if (!csv)
        {
            throw std::runtime_error("solve_rk4: could not open output file");
        }
        csv << std::setprecision(15);

        ///                   Header                 ///
        csv << "t"; // Time column
        for (integer i = 0; i < n; ++i) // State columns
        {
            csv << "," << problem.state_name(i);
        }

        // write algebraic output names if available
        // t, x1, x2, ..., ia, ib, ic, ...
        for (integer i = 0; i < n_outputs; ++i)
        {
            csv << "," << problem.output_name(i);
        }

        if (has_exact)
        {
            for (integer i = 0; i < n; ++i)
            {
                csv << "," << problem.state_name(i);
            }

            for (integer i = 0; i < n; ++i)
            {
                csv << ",error_x" << i + 1; // Error columns
            }
        }
        csv << "\n";


        real_type t = t0;
        integer step = 0;

        while (t <= tf + 1e-12)
        {
            const bool save_this_step = (step % save_every == 0) || (t >= tf);

            if (save_this_step)
            {
                csv << t;

                for (integer i = 0; i < n; ++i)
                {
                    csv << "," << x(i);
                }

                 // Algebraic outputs
                const vec_type y = problem.outputs(t, x);

                if (y.size() != n_outputs)
                {
                    throw std::runtime_error("solve_rk4: invalid output vector size");
                }

                for (integer i = 0; i < n_outputs; ++i)
                {
                    csv << "," << y(i);
                }

                if (has_exact)
                {
                    for (integer i = 0; i < n; ++i)
                    {
                        x_exact(i) = problem.exact(t, i);
                        csv << "," << x_exact(i);
                    }

                    for (integer i = 0; i < n; ++i)
                    {
                        csv << "," << std::abs(x(i) - x_exact(i));
                    }
                }

                csv << "\n";
            }

            if (t >= tf)
            {
                break;
            }

            const real_type h_step = std::min(h, tf - t);

            rk4_step_inplace(problem,
                            t,
                            x,
                            h_step,
                            x_next,
                            k1,
                            k2,
                            k3,
                            k4,
                            xtmp);

            x.swap(x_next);
            t += h_step;
            ++step;
        }
    }

    void solve_rk4(const ODE_Problem_base& problem,
               const SolverOptions& options)
    {
        solve_rk4(problem,
                options.h,
                options.output_file,
                options.save_every);
    }

    vec_type ODE_Problem_base::outputs(real_type t, const vec_type& x) const
    {
        vec_type y(n_outputs());

        for (integer i = 0; i < n_outputs(); ++i)
        {
            y(i) = output(t, x, i);
        }

        return y;
    }
}