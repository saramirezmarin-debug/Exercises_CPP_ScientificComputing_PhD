#pragma once

// Old way to do pragma once
#ifndef EXERCISE2_HH
#define EXERCISE2_HH

// {fmt}
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
#include <fmt/chrono.h>

// Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <Eigen/Cholesky>
#include <Eigen/LU>

// Standard
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <windows.h>

/**
 * @file Exercise2.hh
 * @brief Data structures and interfaces for the Exercise2 problem.
 *
 * This file contains:
 * - Eigen-based numeric aliases;
 * - a block tridiagonal linear solver;
 * - an abstract base class for problems of the form
 *   \f$-Y''(x) + A(x)Y(x) = F(x)\f$.
 *
 * Unicode symbols are intentionally used in the comments to make the
 * relation between mathematics and implementation clearer.
 */

namespace Exercise2
{
    /// Floating point type used throughout the project.
    using real_type = double;

    /// Integer type used for indices and sizes.
    using integer = int;

    /// Dynamic dense matrix of real scalars
    using mat_type = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

    /// Dynamic columns vector of real scalars
    using vec_type = Eigen::Matrix<real_type, Eigen::Dynamic, 1>;

    /**
     * @brief Start a lightweight global timer.
     *
     * It is meant for simple local timing measurements:
     * `tic(); ...; toc();`
     */
    void tic();

    /**
     * @brief Return the elapsed time in milliseconds since the last call to `tic()`.
     *
     * `tic(); ...; toc();`
     */
    real_type toc(); // return elapsed time in milliseconds

    /**
     * @brief Container for a block tridiagonal linear system.
     *
     * The matrix has the structure
     *
     * \f[
     * \begin{bmatrix}
     * D_0 & U_0 \\
     * L_0 & D_1 & U_1 \\
     *     & L_1 & D_2 & \ddots \\
     *     &     & \ddots & \ddots & U_{m-2} \\
     *     &     &        & L_{m-2} & D_{m-1}
     * \end{bmatrix}
     * \f]
     *
     * where each block has size \f$n \times n\f$.
     *
     * The implementation stores separately:
     * - the diagonal blocks `D_bloks`,
     * - the lower off-diagonal blocks `L_bloks`,
     * - the upper off-diagonal blocks `U_bloks`,
     * - the right-hand side and the solution.
     *
     * The system is solved with a block Thomas algorithm,
     * interpreted as an `LDU` factorization.
     */

    class Block_Tridiagonal
    {

        /// Size of each square block.
        integer _n = 0;

        /// Total number of blocks along the diagonal.
        integer _n_blk = 0;

        /// Diagonal blocks, stored side by side: `n × (n·n_blk)`.
        mat_type D_bloks;

        /// Lower off-diagonal blocks, stored side by side: `n × (n·(n_blk-1))`.
        mat_type L_bloks;

        /// Upper off-diagonal blocks, stored side by side: `n × (n·(n_blk-1))`.
        mat_type U_bloks;

        /// Global right-hand side, stored as the concatenation of blocks `b_i`.
        vec_type B_bloks;

        /// Global solution, stored as the concatenation of blocks `x_i`.
        vec_type X_bloks;

    public:
        /// Default constructor
        Block_Tridiagonal() = default;

        /**
         * @brief Allocate the data structure.
         * @param n size of each square block.
         * @param n_blk number of diagonal blocks.
         *
         * After `setup`, all arrays are reset to zero.
         */
        void setup(integer n, integer n_blk);

        /// @return block size.
        integer block_size() const { return _n; }

        /// @return number of diagonal blocks.
        integer num_blocks() const { return _n_blk; }

        /**
         * @brief Solve the linear system with the block Thomas algorithm.
         *
         * The procedure performs:
         * - forward elimination on the modified diagonal blocks;
         * - local dense solves on the blocks;
         * - backward substitution.
         *
         * @throw std::runtime_error if the structure is not initialized or
         *         if one of the modified diagonal blocks is singular.
         */
        void solve();

        /// @brief Extract a copy of diagonal block `blk`.
        mat_type diagonal_block(integer blk) const { return D_bloks.middleCols(blk * _n, _n); }

        /// @brief Extract a copy of lower block `blk`.
        mat_type lower_block(integer blk) const { return L_bloks.middleCols(blk * _n, _n); }

        /// @brief Extract a copy of upper block `blk`.
        mat_type upper_block(integer blk) const { return U_bloks.middleCols(blk * _n, _n); }

        /// @brief Extract a copy of right-hand side block `blk`.
        vec_type rhs_block(integer blk) const { return B_bloks.segment(blk * _n, _n); }

        /// @brief Extract a copy of solution block `blk`.
        vec_type solution_block(integer blk) const { return X_bloks.segment(blk * _n, _n); }

        /// @name Pointwise coefficient access
        /// @{
        real_type D(integer blk, integer i, integer j) const { return D_bloks(i, blk * _n + j); }
        real_type &D(integer blk, integer i, integer j) { return D_bloks(i, blk * _n + j); }

        real_type L(integer blk, integer i, integer j) const { return L_bloks(i, blk * _n + j); }
        real_type &L(integer blk, integer i, integer j) { return L_bloks(i, blk * _n + j); }

        real_type U(integer blk, integer i, integer j) const { return U_bloks(i, blk * _n + j); }
        real_type &U(integer blk, integer i, integer j) { return U_bloks(i, blk * _n + j); }

        real_type b(integer blk, integer i) const { return B_bloks(blk * _n + i); }
        real_type &b(integer blk, integer i) { return B_bloks(blk * _n + i); }
        /// @}

        /// @brief Read-only access to component `i` of solution block `blk`.
        real_type x(integer blk, integer i) const { return X_bloks(blk * _n + i); }
    };

    /**
     * @brief Abstract interface for vector boundary value problems of the form
     *        \f$-Y''(x) + A(x)Y(x) = F(x)\f$.
     *
     * Each concrete problem must define:
     * - the dimension of the unknown vector `Y`,
     * - the coefficient matrix `A(x)`,
     * - the forcing term `F(x)`,
     * - the left/right boundary conditions,
     * - the interval `[a,b]` where the solution is sought.
     *
     * The exact solution is optional:
     * - by default it is not available (`has_exact_solution() = false`);
     * - when available, it can be used for error and convergence estimates.
     */
    class BVP_Problem_base
    {
    public:
        /// Default constructor.
        BVP_Problem_base() = default;

        /// @return number of components of the unknown vector `Y(x)`.
        virtual integer n() const = 0;

        /// @return left endpoint `a` of the interval `[a,b]`.
        virtual real_type a() const = 0;

        /// @return right endpoint `b` of the interval `[a,b]`.
        virtual real_type b() const = 0;

        /**
         * @brief Matrix coefficient `A(x)` of the continuous problem.
         * @param x evaluation point.
         * @param i row index.
         * @param j column index.
         * @return value of entry `(i,j)` of `A(x)`.
         */
        virtual real_type A(real_type x, integer i, integer j) const = 0;

        /**
         * @brief Forcing term `F(x)`.
         * @param x evaluation point.
         * @param i requested component of `F(x)`.
         * @return value of component `i`.
         */
        virtual real_type F(real_type x, integer i) const = 0;

        /// @brief Value of component `i` at the left boundary `x = a`.
        virtual real_type left_BC(integer i) const = 0;

        /// @brief Value of component `i` at the right boundary `x = b`.
        virtual real_type right_BC(integer i) const = 0;

        /// @return `true` if an analytical exact solution is available.
        virtual bool has_exact_solution() const { return false; }

        /**
         * @brief Analytical exact solution, when available.
         * @param x evaluation point.
         * @param i requested component.
         * @return exact value of component `i`.
         *
         * The default is `0`, a neutral fallback when the exact solution
         * is not known.
         */
        virtual real_type exact([[maybe_unused]] real_type x, [[maybe_unused]] integer i) const { return 0; }
    };
}

#endif