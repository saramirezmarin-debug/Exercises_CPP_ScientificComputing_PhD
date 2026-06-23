#pragma once

#include "ODE4.hh"
#include <string>

struct SimulationConfig
{
    ODE::real_type h = 1e-5;
    ODE::integer save_every = 1;
    std::string output_file = "results/output.csv";
};