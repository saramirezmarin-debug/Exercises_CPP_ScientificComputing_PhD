#pragma once

// Eigen
#include <Eigen/Dense>

// Standard
#include <iostream>
#include <algorithm>
#include <string>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <iomanip> 

namespace ODE
{
    /// Floating point type
    using real_type = double;

    /// Integer type
    using integer = int;

    /// Dynamic dense matrix
    using mat_type = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

    /// Dynamic column vector
    using vec_type = Eigen::Matrix<real_type, Eigen::Dynamic, 1>;

    /**
     * @brief Abstract interface for initial value problems of the form
     *        dx/dt = f(t,x).
     */
    class ODE_Problem_base
    {
    public:
        /// Default constructor.
        ODE_Problem_base() = default;

        /// Virtual destructor.
        virtual ~ODE_Problem_base() = default;

        /// @return number of states.
        virtual integer n() const = 0;

        /**
         * @brief Initial condition.
         * @param i state index.
         * @return initial value of x_i(t0).
         */
        virtual real_type initial_condition(integer i) const = 0;

        virtual std::string state_name([[maybe_unused]] integer i) const
        {
            return "x" + std::to_string(i + 1);
        }

        /**
         * @brief Initial time.
         */
        virtual real_type t0() const = 0;

        /**
         * @brief Final simulation time.
         */
        virtual real_type tf() const = 0;

        /**
         * @brief Right-hand side of the ODE system.
         *
         * Computes:
         * dxdt = f(t,x)
         *
         * @param t current time.
         * @param x current state vector.
         * @param dxdt output derivative vector.
         */
        virtual void rhs(real_type t, const vec_type& x, vec_type& dxdt) const = 0;

        /**
         * @brief Exact solution, when available.
         * @param t evaluation time.
         * @param i state index.
         * @return exact value of x_i(t).
         */
        virtual real_type exact([[maybe_unused]] real_type t,
                                [[maybe_unused]] integer i) const
        {
            return 0;
        }

        /**
         * @brief Whether an exact solution is available.
         */
        virtual bool has_exact_solution() const
        {
            return false;
        }

        // ----------------------------------------------------------------
        // Optional algebraic outputs 
        // ----------------------------------------------------------------
        virtual integer n_outputs() const
        {
            return 0;
        }
        // name of the i-th algebraic output.
        virtual std::string output_name([[maybe_unused]] integer i) const
        {
            return "";
        }
        // value of the i-th algebraic output at time t and state x.
        virtual real_type output([[maybe_unused]] real_type t,
                                 [[maybe_unused]] const vec_type& x,
                                 [[maybe_unused]] integer i) const
        {
            return 0.0;
        }
    };

    /**
     * @brief One Runge-Kutta 4 integration step.
     *
     * Computes x_{k+1} from x_k using:
     *
     * k1 = f(t, x)
     * k2 = f(t+h/2, x+h*k1/2)
     * k3 = f(t+h/2, x+h*k2/2)
     * k4 = f(t+h,   x+h*k3)
     *
     * x_{k+1} = x_k + h/6 * (k1 + 2k2 + 2k3 + k4)
     */
    vec_type rk4_step(const ODE_Problem_base& problem,
                      real_type t,
                      const vec_type& x,
                      real_type h);

    /**
     * @brief Solve an initial value problem using RK4 and export results to CSV.
     */
    void solve_rk4(const ODE_Problem_base& problem,
                   real_type h,
                   const std::string& filename,
                   integer save_every = 1);
}
