#include "ODE4.hh"

class Oscillator2D : public ODE::ODE_Problem_base
{
public:
    ODE::integer n() const override
    {
        return 2;
    }

    ODE::real_type t0() const override
    {
        return 0.0;
    }

    ODE::real_type tf() const override
    {
        return 20.0;
    }

    ODE::real_type initial_condition(ODE::integer i) const override
    {
        if (i == 0) return 1.0; // x1(0)
        if (i == 1) return 0.0; // x2(0)

        throw std::out_of_range("Invalid state index");
    }

    void rhs(ODE::real_type t,
             const ODE::vec_type& x,
             ODE::vec_type& dxdt) const override
    {
        (void)t;

        dxdt(0) =  x(1);   // dx1/dt = x2
        dxdt(1) = -x(0);  // dx2/dt = -x1
    }

    bool has_exact_solution() const override
    {
        return true;
    }

    ODE::real_type exact(ODE::real_type t, ODE::integer i) const override
    {
        if (i == 0) return std::cos(t);
        if (i == 1) return -std::sin(t);

        throw std::out_of_range("Invalid state index");
    }
};

int main()
{
    Oscillator2D problem;

    const ODE::real_type h = 0.01; // Time step size

    ODE::solve_rk4(problem, h, "results/oscillator2d.csv");

    return 0;
}