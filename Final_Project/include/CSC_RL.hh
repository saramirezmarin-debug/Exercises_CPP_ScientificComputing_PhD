#pragma once

#include "ODE4.hh"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

// ============================================================
// Generic stair signal
//
//
// Example:
//
// times  = {0.0, 1.0, 2.0, 3.0}
// values = {16.0, 50.0, 100.0, 80.0}
//
// Then:
//
// 0.0 <= t < 1.0  -> 16.0
// 1.0 <= t < 2.0  -> 50.0
// 2.0 <= t < 3.0  -> 100.0
// t >= 3.0        -> 80.0
//
// If times and values are empty, the signal uses constant_value.
// ============================================================
struct StairSignal
{

    ODE::real_type constant_value = 0.0;    // Default value if no stair is provided
    std::vector<ODE::real_type> times;      // Switching times of the stair signal
    std::vector<ODE::real_type> values;     // Values of the stair signal at the corresponding times

    ODE::real_type value(ODE::real_type t) const
    {
        if (times.empty() && values.empty())
        {
            return constant_value; // Return constant value if no stair is defined
        }

        if (times.empty() || values.empty())
        {
            throw std::runtime_error("StairSignal: times and values must both be non-empty or both empty");
        }

        // Ensure times and values have the same size
        if (times.size() != values.size())
        {
            throw std::runtime_error("StairSignal: times and values must have the same size");
        }

        // If t is before the first switching time,
        // use the first defined value.
        if (t < times.front())
        {
            return values.front();
        }

        // upper_bound gives the first switching time strictly greater than t.
        // The active value is therefore one index before that.
        auto it = std::upper_bound(times.begin(), times.end(), t);

        const std::size_t index =
            static_cast<std::size_t>(std::distance(times.begin(), it) - 1);

        return values[index];
    }

    void validate(const std::string& name) const
    {
        // Empty vectors are allowed.
        // In that case the signal behaves as a constant.
        if (times.empty() && values.empty())
        {
            return;
        }

        // It is invalid to define only times or only values.
        if (times.empty() || values.empty())
        {
            throw std::runtime_error(name + ": times and values must both be empty or both be non-empty");
        }

        // Each switching time must have one associated value.
        if (times.size() != values.size())
        {
            throw std::runtime_error(name + ": times and values must have the same size");
        }

        // The stair logic assumes that times are sorted.
        if (!std::is_sorted(times.begin(), times.end()))
        {
            throw std::runtime_error(name + ": times must be sorted in ascending order");
        }
    }
};

struct CSC_RL_Parameters
{
    // ------------------------------------------------------------
    // PLL parameters
    // ------------------------------------------------------------
    ODE::real_type w0_pll;      // Nominal angular frequency [rad/s]
    ODE::real_type kp_pll;      // PLL proportional gain
    ODE::real_type ki_pll;      // PLL integral gain
    ODE::real_type Vdq_nom;      // Phase peak voltage used for normalization

    ODE::real_type theta_hat0; // Initial angle of the PLL's internal oscillator
    ODE::real_type xi_pll0;    // Initial frequency deviation of the PLL's internal oscillator

    // ------------------------------------------------------------
    // CSC AC filter parameters
    // ------------------------------------------------------------
    ODE::real_type Rf;
    ODE::real_type Lf;
    ODE::real_type Cf;

    // ------------------------------------------------------------
    // DC-side parameters
    // ------------------------------------------------------------
    ODE::real_type Ldc;
    ODE::real_type RL;  // Load resistance on DC side

    // ------------------------------------------------------------
    // Synchronous reference frame
    // ------------------------------------------------------------
    ODE::real_type omega;

    // Grid parameters
    ODE::real_type Lg;
    ODE::real_type Cg;
    ODE::real_type Rg;
    ODE::real_type egd;
    ODE::real_type egq;

    // ------------------------------------------------------------
    // Initial conditions
    // x = [igd, igq, ed, eq, id, iq, vd, vq, istk theta_hat xi_pll]^T
    // ------------------------------------------------------------

    // PCC states
    ODE::real_type igd0;
    ODE::real_type igq0;
    ODE::real_type ed0;
    ODE::real_type eq0;

    // CSC states
    ODE::real_type id0;
    ODE::real_type iq0;
    ODE::real_type vd0;
    ODE::real_type vq0;
    ODE::real_type istk0;

    // ------------------------------------------------------------
    // Inner loop controller parameters
    //
    // Stage 2: current error -> capacitor voltage reference
    // Stage 1: capacitor voltage error -> modulation index
    // ------------------------------------------------------------
    ODE::real_type kp1;
    ODE::real_type ki1;

    ODE::real_type kp2;
    ODE::real_type ki2;

    ODE::real_type Idc_min;

    // // Current references
    // StairSignal id_r;
    // StairSignal iq_r;

    // ------------------------------------------------------------
    // Inner loop initial conditions
    // x_inner = [xi_ucd, xi_ucq, xi_isd, xi_isq]^T
    // ------------------------------------------------------------
    ODE::real_type xi_ucd0;
    ODE::real_type xi_ucq0;
    ODE::real_type xi_isd0;
    ODE::real_type xi_isq0;

    // ------------------------------------------------------------
    // Outer loop initial conditions
    // ------------------------------------------------------------
    ODE::real_type kpO;      // Outer-loop proportional gain
    ODE::real_type kiO;      // Outer-loop integral gain
    ODE::real_type V2_min;

    StairSignal idc_ref;     // DC current reference [A]
    StairSignal Q_ref;       // Reactive power reference [var]

    ODE::real_type xi_idc2_0;

    // ------------------------------------------------------------
    // Simulation time
    // ------------------------------------------------------------
    ODE::real_type t0;
    ODE::real_type tf;
};

class CSC_RL : public ODE::ODE_Problem_base
{
private:
    CSC_RL_Parameters p_;

    void validate_parameters() const;

public:
    explicit CSC_RL(const CSC_RL_Parameters& parameters);

    ODE::integer n() const override;

    ODE::real_type t0() const override;

    ODE::real_type tf() const override;

    ODE::real_type initial_condition(ODE::integer i) const override;

    void rhs(ODE::real_type t,
             const ODE::vec_type& x,
             ODE::vec_type& dxdt) const override;

    ODE::integer n_outputs() const override; // Number of algebraic outputs
    std::string output_name(ODE::integer i) const override; // Name of algebraic output i
    ODE::real_type output(ODE::real_type t, const ODE::vec_type& x, ODE::integer i) const override; // value of each algebraic output
};
