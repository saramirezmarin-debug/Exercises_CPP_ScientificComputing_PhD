#include "Exercise2.hh"

namespace Exercise2 {

  /// Local alias for the high-resolution clock used in the micro-benchmark.
  using clock_type = std::chrono::high_resolution_clock;

  /// Global start time used by the `tic`/`toc` pair.
  static clock_type::time_point m_start;

  /**
   * @brief Save the current time as the initial reference.
   *
   * The implementation is intentionally minimal: it is only used to
   * measure the cost of the linear solve in the test program.
   */
  void tic() {
    m_start = clock_type::now();
  }

  /**
   * @brief Return the time elapsed since `tic()` in milliseconds.
   * @return elapsed time in milliseconds.
   */
  real_type toc() {
    auto now = clock_type::now();
    return std::chrono::duration<real_type,std::milli>(now- m_start).count();
  }
  
  /**
   * @brief Allocate and zero the block tridiagonal structure.
   *
   * Memory is organized in "stripes":
   * - diagonal blocks are placed side by side in `D_bloks`,
   * - the same idea is used for lower and upper off-diagonal blocks.
   *
   * This makes both pointwise access and extraction of whole blocks
   * through `middleCols` straightforward.
   */
  void Block_Tridiagonal::setup( integer n, integer n_blk ) {
    if ( n <= 0 )     throw std::invalid_argument("Block_Tridiagonal::setup: block size must be positive");
    if ( n_blk <= 0 ) throw std::invalid_argument("Block_Tridiagonal::setup: number of blocks must be positive");
    // Store the internal dimensions
    _n     = n;
    _n_blk = n_blk;

    // Resize the arrays
    D_bloks.resize( n, n*n_blk );
    L_bloks.resize( n, n*std::max<integer>(n_blk-1,0) );
    U_bloks.resize( n, n*std::max<integer>(n_blk-1,0) );
    B_bloks.resize( n*n_blk );
    X_bloks.resize( n*n_blk );

    // Reset all entries to zero.
    D_bloks.setZero();
    L_bloks.setZero();
    U_bloks.setZero();
    B_bloks.setZero();
    X_bloks.setZero();
  }
    
  /**
   * @brief Solve the system with block Thomas ⇔ `LDU` factorization.
   *
   * Numerical scheme:
   * - step 1: recursive construction of the modified diagonal blocks;
   * - step 2: solve for the modified upper blocks `Û_i`;
   * - step 3: forward elimination on the right-hand side;
   * - step 4: backward substitution on the solution.
   *
   * No explicit block inverse is ever formed:
   * each occurrence of `D^{-1}` is implemented as a local dense solve.
   */
  void Block_Tridiagonal::solve() {
    if ( _n <= 0 || _n_blk <= 0 ) {
      throw std::runtime_error("Block_Tridiagonal::solve: system not initialized");
    }

    // Using LU with full pivoting -> A = LU or PAQ = LU with pivotings P and Q.
    /// `FullPivLU` is robust and simple for small dense blocks.
    using LU_type = Eigen::FullPivLU<mat_type>;

    /// Modified upper blocks `Û_i` of the factorization.
    std::vector<mat_type> U_hat( std::max<integer>(_n_blk - 1, 0) );

    /// Modified right-hand side after forward elimination.
    std::vector<vec_type> rhs_hat( _n_blk );

    /// First modified diagonal block: it coincides with `D₀`.
    mat_type D_hat = diagonal_block(0); 
    LU_type  lu( D_hat );
    if ( !lu.isInvertible() ) {
      throw std::runtime_error("Block_Tridiagonal::solve: singular first diagonal block");
    }

    // First step of the recurrence:
    //   Û₀ = D̂₀⁻¹ U₀,   r̂₀ = D̂₀⁻¹ b₀.
    if ( _n_blk > 1 ) U_hat[0] = lu.solve( upper_block(0) );
    rhs_hat[0] = lu.solve( rhs_block(0) );

    // Forward elimination on the remaining blocks.
    for ( integer blk = 1; blk < _n_blk; ++blk ) {
      // D̂ᵢ = Dᵢ - Lᵢ₋₁ Ûᵢ₋₁.
      D_hat = diagonal_block(blk) - lower_block(blk-1) * U_hat[blk-1];
      lu.compute( D_hat );
      if ( !lu.isInvertible() ) {
        throw std::runtime_error(
          fmt::format("Block_Tridiagonal::solve: singular diagonal block {}", blk)
        );
      }
      // Compute the modified upper block only when it is actually needed.
      if ( blk < _n_blk-1 ) U_hat[blk] = lu.solve( upper_block(blk) );

      // r̂ᵢ = D̂ᵢ⁻¹ ( bᵢ - Lᵢ₋₁ r̂ᵢ₋₁ ).
      rhs_hat[blk] = lu.solve( rhs_block(blk) - lower_block(blk-1) * rhs_hat[blk-1] );
    }

    // Backward substitution:
    //   x_{m-1} = r̂_{m-1},
    //   x_i     = r̂_i - Û_i x_{i+1}.
    X_bloks.segment( (_n_blk-1)*_n, _n ) = rhs_hat[_n_blk-1];
    for ( integer blk = _n_blk-2; blk >= 0; --blk ) {
      X_bloks.segment( blk*_n, _n ) = rhs_hat[blk] - U_hat[blk] * X_bloks.segment( (blk+1)*_n, _n );
      if ( blk == 0 ) break;
    }
  }
}


