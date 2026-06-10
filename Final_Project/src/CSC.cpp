#include "CSC.hh"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

// ------------------------------------------------------------
// Dynamic system:
// dx/dt = f(t, x)
// Example: dx/dt = -x
// ------------------------------------------------------------
double rhs(double t, double x)
{
    (void)t; // t is not used in this example
    return -x;
}

// ------------------------------------------------------------
// One RK4 integration step
// ------------------------------------------------------------
double rk4_step(double t, double x, double h)
{
    double k1 = rhs(t, x);
    double k2 = rhs(t + 0.5 * h, x + 0.5 * h * k1);
    double k3 = rhs(t + 0.5 * h, x + 0.5 * h * k2);
    double k4 = rhs(t + h,       x + h * k3);

    return x + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

int main()
{
    // Simulation parameters
    double t0 = 0.0;
    double tf = 5.0;
    double h  = 0.01;

    // Initial condition
    double x = 1.0;
    double t = t0;

    std::ofstream csv("results/CSC_rk4.csv");

    if(!csv)
    {
        fmt::print(fg(fmt::color::red), "Error: could not open output/CSC_rk4.csv for writing\n");
        return 1;
    }

    // CSV header
    csv << "t,x_rk4,x_exact,error\n";


    while (t <= tf)
    {
        double x_exact = std::exp(-t);
        double error = std::abs(x - x_exact);

        csv << std::setprecision(15)
             << t << ","
             << x << ","
             << x_exact << ","
             << error << "\n";

        x = rk4_step(t, x, h);
        t += h;
    }
    

    csv.close();

    fmt::print("Simulation finished.\n");
    fmt::print("Results saved in CSC_rk4.csv\n");

    return 0;
}