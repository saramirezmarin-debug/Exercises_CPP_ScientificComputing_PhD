#pragma once

#ifndef BVP_HH
#define BVP_HH

#include <iostream>
#include <fstream>

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
#include <complex>
#include <chrono>
#include <random>
#include <cmath>
#include <cstdint>

namespace BVP 
{
    // alias of double
    using real_type = double;

    //timing
    void tic();
    real_type toc(); //return elapsed time in milliseconds

    using integer = int;

    // alias for Eigen

    // define mat_type as a dynamic matrix of double
    using mat_type = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;
    using vec_type = Eigen::Matrix<real_type, Eigen::Dynamic, 1>;

    class Block_Tridiagonal
    {
        //
        // / 
        // |    D  U
        // |    L  D  U
        // |       L  D  U
        // |          L  D
        //
        //
        //

        integer _n      = 0; // size of the block
        integer _n_blk  = 0; // number of blocks

        // storage of tridiagonal system and the solution
        mat_type D_blocks; // n x (n*b_blk) matrix
        mat_type L_blocks; // n x (n*b_blk - 1) matrix
        mat_type U_blocks; // n x (n*b_blk - 1) matrix  

        vec_type B_blocks; // (n*b_blk) vector
        vec_type X_blocks; // (n*b_blk) vector
        
    public:

        Block_Tridiagonal() {}


        // allocate memory for the blocks
        void setup(integer n, integer n_blk);

        // acces the elements
        real_type   D(integer blk, integer i, integer j) const;
        real_type & D(integer blk, integer i, integer j);

        real_type   L(integer blk, integer i, integer j) const;
        real_type & L(integer blk, integer i, integer j);

        real_type   U(integer blk, integer i, integer j) const;
        real_type & U(integer blk, integer i, integer j);

        real_type   b(integer blk, integer i) const;
        real_type & b(integer blk, integer i);

        //access the solution vector
        real_type   x(integer blk, integer i) const;
        real_type & x(integer blk, integer i);
    
    };

    class BVP_Problem {


    };
}


#endif