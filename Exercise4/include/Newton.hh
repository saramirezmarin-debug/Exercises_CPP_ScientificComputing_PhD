#pragma once

#ifndef NEWTON_HH
#define NEWTON_HH

#include "dual.hh"

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

// ---------------------------------------------------------------------------
// Type aliases
// ---------------------------------------------------------------------------
using Vector2d = Eigen::Vector2d;
using Vector3d = Eigen::Vector3d;
using Vector4d = Eigen::Vector4d;
using VectorXd = Eigen::VectorXd;

using Matrix2d = Eigen::Matrix2d;
using Matrix3d = Eigen::Matrix3d;
using Matrix4d = Eigen::Matrix4d;
using MatrixXd = Eigen::MatrixXd;

using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;


namespace AD {

  template <typename T>
  class NewtonProblem {
  public:
    using d_vec   = Eigen::Matrix<T,Eigen::Dynamic,1>;
    using d_mat   = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;
    using integer = int;

  protected:
    d_vec _x0; // initial point
  
  public:
    NewtonProblem() = default;
    
    integer dim() const { return _x0.size(); }
  
    // compute F(X) for the problem
    virtual d_vec F( d_vec const & X ) const = 0;

    // compute Jacobian JF(X) for the problem
    virtual d_mat JF( d_vec const & X ) const = 0;

    d_vec initial_point() const { return _x0; }
  
  };

  template <typename T>
  class NewtonProblemAD : public NewtonProblem<T> {

  public:
    using d_dual_vec = Eigen::Matrix<Dual<T>,Eigen::Dynamic,1>;
    using integer    = typename NewtonProblem<T>::integer;

  public:
  
    using typename NewtonProblem<T>::d_vec;
    using typename NewtonProblem<T>::d_mat;
  
    NewtonProblemAD() = default;

    // compute F(X) for the problem
    virtual d_dual_vec F( d_dual_vec const & X ) const = 0;
  
    // compute F(X) for the problem
    d_vec F( d_vec const & X ) const override {
      d_dual_vec Xd, dual_res;
      d_vec      res;
      integer N = X.size();
      Xd.resize( X.size() );
      Xd       = X.template cast<Dual<T>>();   // convert from double to dual
      dual_res = this->F( Xd );
      res.resize(N);
      for ( integer i = 0; i < N; ++i )
        res(i) = dual_res(i).value(); // convert from dual to double
      return res;
    }

    // compute Jacobian JF(X) for the problem
    d_mat JF( d_vec const & X ) const override {
      d_dual_vec Xd, column;
      d_mat      JF;
      integer N = X.rows();
      JF.resize( N, N ); // Jacobian allocation
      Xd.resize( N );
      column.resize( N );
      Xd = X.template cast<Dual<T>>();   // convert from double to dual
      Dual<T> D(0,1); // 0 + 1 * epsilon
      for ( integer i = 0; i < N; ++i ) {
        // set x(i) = x(i) + epsilon --> compute derivative orespect to x(i)
        Xd(i) += D;
        column.noalias() = this->F( Xd );
        Xd(i) -= D;
        for ( integer j = 0; j < N; ++j ) {
          JF( j, i ) = column(j).dual();
        }
      }
      return JF;
    }
  };
  
  template <typename T>
  class NewtonSolver {
    using integer = typename NewtonProblem<T>::integer;
    using d_vec   = typename NewtonProblem<T>::d_vec;
    using d_mat   = typename NewtonProblem<T>::d_mat;
  
    NewtonProblem<T> const * _problem = nullptr;
    
    integer _max_iter     = 100;
    integer _max_sub_iter = 50;
    T       _damp_factor  = 0.5;
    T       _tolerance    = 1e-9;
    
    bool    _converged = false;
    d_vec   _x; // solution
  public:

    NewtonSolver() = delete;
    explicit NewtonSolver( NewtonProblem<T> const * p ) : _problem( p ) {}
    
    bool converged() const { return _converged; }
    d_vec solution() const { return _x; }
    
    void
    solve() {
      // initialize the loop
      d_vec x0, x1, d0, d1, F0, F1;
      d_mat JF0;
      integer N = _problem->dim();
      x0.resize(N);
      x1.resize(N);
      F0.resize(N);
      F1.resize(N);
      JF0.resize(N,N);
      d0.resize(N);
      d1.resize(N);

      _converged = false;
      x0         = _problem->initial_point();
      for ( integer iter = 0; iter < _max_iter; ++iter ) {
        F0 = _problem->F( x0 );
        // check || F(x0) || for temination
        T F0_norm = F0.template lpNorm<Eigen::Infinity>();
        fmt::print( "{} ||F|| = {}\n", iter, F0_norm );
        _converged = F0_norm <= _tolerance;
        if ( _converged )
        {
            _x = x0;
            return;
        }   
        

        // compute direction
        JF0 = _problem->JF( x0 );
        d0  = -JF0.fullPivLu().solve(F0);
        T d0norm2 = d0.dot(d0);
        // damp loop
        T alpha = 1;
        bool ok = false;
        for ( integer i = 0; !ok && i < _max_sub_iter; ++i ) {
          x1 = x0 + alpha * d0;
          // check Deufhard condition for x1
          F1 = _problem->F( x1 );
          d1 = -JF0.fullPivLu().solve(F1);
          T d1norm2 = d1.dot(d1);
          ok = d1norm2 <= d0norm2*(1-alpha/2);
          if ( !ok ) alpha *= _damp_factor;
        }
        if ( !ok ) return; // failed Newton
        x0 = x1; // step accepted
        //_x = x0; // update solution
      }
    }
  };
}

#endif