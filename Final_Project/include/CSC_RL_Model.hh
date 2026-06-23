#pragma once

#include "CSC_RL.hh"
#include "dual.hh"

#include <cmath>
#include <algorithm>

class CSC_RL_Model
{
public:

    template <typename T>
    struct Input
    {
        T idc_ref;
        T Q_ref;
    };

    template <typename T>
    struct Control
    {
        T e_pll;
        T w_hat;

        T idc_ref;
        T Qref;
        T Pref;

        T idr;
        T iqr;

        T Eid;
        T Eiq;
        T Eidc2;

        T vdr;
        T vqr;

        T Evd;
        T Evq;

        T md;
        T mq;
    };

    template <typename T>
    Input<T> input_from_parameters(const CSC_RL_Parameters& p,
                                   ODE::real_type t) const
    {
        Input<T> input;

        input.idc_ref = T(p.idc_ref.value(t));
        input.Q_ref   = T(p.Q_ref.value(t));

        return input;
    }

    template <typename Vec>
    Control<typename Vec::Scalar>
    control(const CSC_RL_Parameters& p,
            const Vec& x,
            const Input<typename Vec::Scalar>& input) const
    {
        using T = typename Vec::Scalar;

        Control<T> c;

        const T ed   = x(ED);
        const T eq   = x(EQ);
        const T id   = x(ID);
        const T iq   = x(IQ);
        const T vd   = x(VD);
        const T vq   = x(VQ);
        const T istk = x(ISTK);

        const T xi_pll  = x(XI_PLL);
        const T xi_ucd  = x(XI_UCD);
        const T xi_ucq  = x(XI_UCQ);
        const T xi_isd  = x(XI_ISD);
        const T xi_isq  = x(XI_ISQ);
        const T xi_idc2 = x(XI_IDC2);

        // ------------------------------------------------------------
        // PLL
        // ------------------------------------------------------------
        c.e_pll = eq / p.Vdq_nom;
        c.w_hat = p.w0_pll + p.kp_pll * c.e_pll + p.ki_pll * xi_pll;

        // ------------------------------------------------------------
        // References
        // ------------------------------------------------------------
        c.idc_ref = input.idc_ref;
        c.Qref    = input.Q_ref;

        // ------------------------------------------------------------
        // Outer loop
        // ------------------------------------------------------------
        c.Eidc2 = c.idc_ref * c.idc_ref - istk * istk;
        c.Pref  = p.kpO * c.Eidc2 + p.kiO * xi_idc2;

        // ------------------------------------------------------------
        // AC current references
        // ------------------------------------------------------------
        const T V2 = ed * ed + eq * eq;

        c.idr = (ed * c.Pref + eq * c.Qref) / V2;
        c.iqr = (eq * c.Pref - ed * c.Qref) / V2;

        c.Eid = c.idr - id;
        c.Eiq = c.iqr - iq;

        // ------------------------------------------------------------
        // Current controller
        // ------------------------------------------------------------
        c.vdr = ed + c.w_hat * p.Lf * iq - p.Rf * id - p.kp2 * c.Eid - p.ki2 * xi_isd;
        c.vqr = eq - c.w_hat * p.Lf * id - p.Rf * iq - p.kp2 * c.Eiq - p.ki2 * xi_isq;
        c.Evd = c.vdr - vd;
        c.Evq = c.vqr - vq;

        // ------------------------------------------------------------
        // Voltage controller
        // ------------------------------------------------------------
        c.md = (id + c.w_hat * p.Cf * vq - p.kp1 * c.Evd - p.ki1 * xi_ucd) / istk;
        c.mq = (iq - c.w_hat * p.Cf * vd - p.kp1 * c.Evq - p.ki1 * xi_ucq) / istk;

        return c;
    }

    template <typename Vec>
    void rhs(const CSC_RL_Parameters& p,
             ODE::real_type t,
             const Vec& x,
             Vec& dxdt,
             const Input<typename Vec::Scalar>& input) const
    {
        using T = typename Vec::Scalar;

        dxdt.resize(NSTATES);

        const T igd  = x(IGD);
        const T igq  = x(IGQ);
        const T ed   = x(ED);
        const T eq   = x(EQ);

        const T id   = x(ID);
        const T iq   = x(IQ);
        const T vd   = x(VD);
        const T vq   = x(VQ);
        const T istk = x(ISTK);

        const Control<T> c = control(p, x, input);

        // ------------------------------------------------------------
        // PCC dynamics
        // ------------------------------------------------------------
        dxdt(IGD) = (p.egd - p.Rg * igd - ed + c.w_hat * p.Lg * igq) / p.Lg;
        dxdt(IGQ) = (p.egq - p.Rg * igq - eq - c.w_hat * p.Lg * igd) / p.Lg;

        dxdt(ED) = (igd - id + c.w_hat * p.Cg * eq) / p.Cg;
        dxdt(EQ) = (igq - iq - c.w_hat * p.Cg * ed) / p.Cg;

        // ------------------------------------------------------------
        // CSC dynamics
        // ------------------------------------------------------------
        dxdt(ID) = (ed - vd - p.Rf * id + c.w_hat * p.Lf * iq) / p.Lf;
        dxdt(IQ) = (eq - vq - p.Rf * iq - c.w_hat * p.Lf * id) / p.Lf;

        dxdt(VD) = (id - c.md * istk + c.w_hat * p.Cf * vq) / p.Cf;
        dxdt(VQ) = (iq - c.mq * istk - c.w_hat * p.Cf * vd) / p.Cf;

        dxdt(ISTK) = (c.md * vd + c.mq * vq - p.RL * istk) / p.Ldc;

        // ------------------------------------------------------------
        // PLL
        // ------------------------------------------------------------
        dxdt(THETA_HAT) = c.w_hat;
        dxdt(XI_PLL)    = c.e_pll;

        // ------------------------------------------------------------
        // Integrators
        // ------------------------------------------------------------
        dxdt(XI_UCD)  = c.Evd;
        dxdt(XI_UCQ)  = c.Evq;
        dxdt(XI_ISD)  = c.Eid;
        dxdt(XI_ISQ)  = c.Eiq;
        dxdt(XI_IDC2) = c.Eidc2;
    }

    template <typename Vec>
    void rhs(const CSC_RL_Parameters& p,
             ODE::real_type t,
             const Vec& x,
             Vec& dxdt) const
    {
        const auto input = input_from_parameters<typename Vec::Scalar>(p, t);
        rhs(p, t, x, dxdt, input);
    }
};