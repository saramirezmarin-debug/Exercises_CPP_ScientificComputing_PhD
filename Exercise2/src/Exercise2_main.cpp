#include "Exercise2.hh"

#include <cstdlib> // for std::atoi -> convert string to integer
#include <limits> // for std::numeric_limits -> get the largest representable value of a type
#include <fstream> // for std::ofstream -> write to files

namespace Exercise2
{

    /**
     * @brief Test problem used in the assignment text.
     *
     * Sistema continuo:
     *
     * - \f$-u'' + 2u - v = f_1\f$
     * - \f$-v'' - u + 2v = f_2\f$
     *
     * with exact solution:
     *
     * - \f$u(x) = \sin(\pi x)\f$
     * - \f$v(x) = \sin(2\pi x)\f$
     *
     * The boundary conditions are homogeneous, but they are still expressed
     * through `exact(a,·)` and `exact(b,·)` so that the code stays consistent
     * with more general cases.
     */
    class BVP_Problem : public BVP_Problem_base
    {
        /// Number of coupled equations.
        integer _n = 2;

        /// π stored in double precision.
        static constexpr real_type pi = 3.14159265358979323846;

        /// Constant coupling matrix `A`.
        Eigen::Matrix<real_type, 2, 2> _A;

    public:
        /// Build the coupling matrix of the problem.
        BVP_Problem()
        {
            _A << 2, -1,
                -1, 2;
        }

        /// @copydoc BVP_Problem_base::n
        integer n() const override { return _n; }

        /// @copydoc BVP_Problem_base::a
        real_type a() const override { return 0; }

        /// @copydoc BVP_Problem_base::b
        real_type b() const override { return 1; }

        /// @copydoc BVP_Problem_base::A
        real_type A([[maybe_unused]] real_type x, integer i, integer j) const override
        {
            return _A(i, j);
        }

        /// @copydoc BVP_Problem_base::F
        real_type F(real_type x, integer i) const override
        {
            switch (i)
            {
            case 0:
                return (pi * pi + 2) * sin(pi * x) - sin(2 * pi * x);
            case 1:
                return (4 * pi * pi + 2) * sin(2 * pi * x) - sin(pi * x);
            }
            throw std::out_of_range("BVP_Problem::F: invalid component index");
        }

        /// @copydoc BVP_Problem_base::left_BC
        real_type left_BC(integer i) const override
        {
            return exact(a(), i);
        }

        /// @copydoc BVP_Problem_base::right_BC
        real_type right_BC(integer i) const override
        {
            return exact(b(), i);
        }

        /// @copydoc BVP_Problem_base::has_exact_solution
        bool has_exact_solution() const override
        {
            return true;
        }

        /// @copydoc BVP_Problem_base::exact
        real_type exact(real_type x, integer i) const override
        {
            switch (i)
            {
            case 0:
                return sin(pi * x);
            case 1:
                return sin(2 * pi * x);
            }
            throw std::out_of_range("BVP_Problem::exact: invalid component index");
        }
    };

    /**
     * @brief Result of one numerical experiment on a given grid.
     */
    struct Experiment_Result
    {
        /// Number of interior grid points in the discrete problem.
        integer N = 0;

        /// Grid size `h = (b-a)/(N+1)`.
        real_type h = 0;

        /// Maximum nodal error.
        real_type max_error = 0;

        /// Order estimate: `log₂(E_old / E_new)`.
        real_type order = std::numeric_limits<real_type>::quiet_NaN();

        /// Time spent solving the linear system.
        real_type cpu_ms = 0;
    };

    /**
     * @brief Run a complete experiment for a given grid size `N`.
     *
     * Steps:
     * - build the uniform grid;
     * - assemble the block tridiagonal matrix;
     * - inject boundary contributions into the right-hand side;
     * - solve the linear system;
     * - compute the error, when the exact solution is available.
     *
     * @param problem continuous problem to discretize.
     * @param N number of interior nodes.
     * @return all relevant data produced by the experiment.
     */
    static Experiment_Result
    run_experiment(BVP_Problem_base const &problem, integer N)
    {
        integer const n = problem.n();

        /// Left endpoint of the search interval.
        real_type const xa = problem.a();

        /// Right endpoint of the search interval.
        real_type const xb = problem.b();

        /// Uniform grid size of the centered finite-difference discretization.
        real_type const h = (xb - xa) / real_type(N + 1);

        /// Factor `1/h²` appearing in all finite-difference coefficients.
        real_type const inv_h2 = 1.0 / (h * h);

        /// Final linear system in block tridiagonal format.
        Block_Tridiagonal system;
        system.setup(n, N);

        for (integer blk = 0; blk < N; ++blk)
        {
            /// Interior node `x_i = a + (i+1)h`.
            real_type const x = xa + (blk + 1) * h;

            for (integer i = 0; i < n; ++i)
            {
                // Discrete RHS: F(x_i) plus the Dirichlet contributions moved
                // from the left/right boundary whenever the node touches a boundary.
                system.b(blk, i) = problem.F(x, i);
                if (blk == 0)
                    system.b(blk, i) += problem.left_BC(i) * inv_h2;
                if (blk == N - 1)
                    system.b(blk, i) += problem.right_BC(i) * inv_h2;

                // Diagonal block:
                //   D_i = A(x_i) + (2/h²) I.
                for (integer j = 0; j < n; ++j)
                {
                    system.D(blk, i, j) = problem.A(x, i, j);
                }
                system.D(blk, i, i) += 2 * inv_h2;
            }
        }

        // Off-diagonal blocks:
        //   L_i = U_i = -(1/h²) I.
        for (integer blk = 0; blk < N - 1; ++blk)
        {
            for (integer i = 0; i < n; ++i)
            {
                for (integer j = 0; j < n; ++j)
                {
                    system.L(blk, i, j) = 0;
                    system.U(blk, i, j) = 0;
                }
                system.L(blk, i, i) = -inv_h2;
                system.U(blk, i, i) = -inv_h2;
            }
        }

        tic();
        system.solve();
        real_type const elapsed_ms = toc();

        // When available, evaluate
        //   E∞ = max_i || y_i^num - y_i^exact ||₂.
        real_type max_error = std::numeric_limits<real_type>::quiet_NaN();
        if (problem.has_exact_solution())
        {
            max_error = 0;
            for (integer blk = 0; blk < N; ++blk)
            {
                real_type const x = xa + (blk + 1) * h;
                real_type err2 = 0;
                for (integer i = 0; i < n; ++i)
                {
                    real_type const diff = system.x(blk, i) - problem.exact(x, i);
                    err2 += diff * diff;
                }
                max_error = std::max(max_error, std::sqrt(err2));
            }
        }

        return {N, h, max_error, std::numeric_limits<real_type>::quiet_NaN(), elapsed_ms};
    }

}

int main(int argc, char **argv)
{
    SetConsoleOutputCP(CP_UTF8);
    using namespace Exercise2;

    /// Concrete problem used in the experiments.
    BVP_Problem problem;

    /// List of `N` values used to assess convergence and timings.
    std::vector<integer> grid_sizes;

    // Store in csv file the convergence history for plotting.
    std::ofstream csv("output/convergence.csv");

    if(!csv)
    {
        fmt::print(fg(fmt::color::red), "Error: could not open output/convergence.csv for writing\n");
        return 1;
    }

    csv << "N,h,max_error,order,cpu_ms\n";

    // If command-line arguments are provided, interpret them as grid sizes;
    // otherwise use the sequence requested in the assignment text.
    if (argc > 1)
    {
        for (int i = 1; i < argc; ++i)
        {
            integer N = std::atoi(argv[i]);
            if (N <= 0)
            {
                fmt::print(fg(fmt::color::red), "Invalid grid size '{}'\n", argv[i]);
                return 1;
            }
            grid_sizes.push_back(N);
        }
    }
    else
    {
        grid_sizes = {5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480};
    }

    fmt::print("Coupled BVP solved with centered finite differences and block Thomas-LDU\n");
    fmt::print("Problem size per node: {}\n\n", problem.n());

    // Unicode table for a cleaner textual report.
    fmt::print("┌──────────┬────────────────┬──────────────────┬──────────────┬────────────┐\n");
    fmt::print("│ {:>8} │ {:>14} │ {:>16} │ {:>12} │ {:>10} │\n", "N", "h", "max error", "order", "time [ms]");
    fmt::print("├──────────┼────────────────┼──────────────────┼──────────────┼────────────┤\n");

    /// Error from the previous run, used for the asymptotic order estimate.
    real_type prev_error = std::numeric_limits<real_type>::quiet_NaN();
    integer const repeats = 20; // number of runs to average timings and mitigate outliers
    for (integer N : grid_sizes)
    {
        // Run various experiments and collect the mean error and CPU time.
        Experiment_Result best_res;
        real_type total_time = 0;

        for(int k = 0; k < repeats; ++k)
        {
            Experiment_Result tmp = run_experiment(problem, N);
            total_time += tmp.cpu_ms;
            if (k == 0)
            {
                best_res = tmp;
            }
        }

        best_res.cpu_ms = total_time / repeats;

        Experiment_Result res = best_res;
        if (std::isfinite(prev_error) && std::isfinite(res.max_error))
        {
            res.order = std::log2(prev_error / res.max_error);
        }

        fmt::print(
            "│ {:8d} │ {:14.6e} │ {:16.8e} │ {:>12} │ {:10.3f} │\n",
            res.N,
            res.h,
            res.max_error,
            std::isfinite(res.order) ? fmt::format("{:.4f}", res.order) : std::string("-"),
            res.cpu_ms);

        prev_error = res.max_error;

        csv << res.N << ","
            << res.h << ","
            << res.max_error << "," 
            << (std::isfinite(res.order) ? res.order : std::numeric_limits<real_type>::quiet_NaN()) << ","
            << res.cpu_ms << "\n";
    }

    fmt::print("└──────────┴────────────────┴──────────────────┴──────────────┴────────────┘\n");

    fmt::print("\nResults exported to output/convergence.csv\n");

    return 0;
}