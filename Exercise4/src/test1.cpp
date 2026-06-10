#include "Newton.hh"

using real = double;

class Problem : public AD::NewtonProblemAD<real> {

public:

  Problem() {
    _x0.resize(3);
    _x0(0) = -1.0; // x1
    _x0(1) =  0.0; // x2
    _x0(2) =  0.0; // x3
  }

  using AD::NewtonProblemAD<real>::F;
  using AD::NewtonProblemAD<real>::JF;

  d_dual_vec F( d_dual_vec const & X ) const override {
    d_dual_vec res;
    res.resize(3);
    auto const & x1 = X(0);
    auto const & x2 = X(1);
    auto const & x3 = X(2);

    constexpr real pi = 3.1415926535897932384626433832795;

    using std::atan;
    using std::atan2;
    using std::sqrt;

    auto theta = (1.0 / (2.0 * pi)) * atan2(x2, x1);


    res(0) = 10.0 * (x3 - 10.0*theta);
    res(1) = 10.0 * (sqrt(x1*x1 + x2*x2) - 1.0);
    res(2) = x3;

    return res;
  }
};

int main() {

  Problem problem;
  
  Problem::d_vec ip = problem.initial_point();
  
  Problem::d_vec F  = problem.F( ip );
  Problem::d_mat JF = problem.JF( ip );
  
  AD::NewtonSolver solver( &problem );
  
  //std::cout << "F = " << F << '\n';
  //std::cout << "JF\n" << JF << '\n';
  
  solver.solve(); // compute solution
  
  if ( solver.converged() ) {
    std::cout << "CONVERGED!!!!\n";
    std::cout << "x = " << solver.solution().transpose() << '\n' ;
  } else {
    std::cout << "FAILED!!!!\n";
  }

  return 0;
}