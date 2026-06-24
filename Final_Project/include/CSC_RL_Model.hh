#pragma once

#include "CSC_RL.hh"
#include "Equilibrium.hh"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdexcept>

/**
 * @file CSC_RL_Model.hh
 * @brief Dynamic model, control laws, and equilibrium formulation for the CSC-RL system.
 *
 * @details
 * This file contains the mathematical core of the simulation. The class
 * CSC_RL_Model evaluates the right-hand side of the nonlinear ODE system used
 * by the RK4 solver and by the Newton-based equilibrium calculation.
 *
 * The model represents a grid-connected PWM-CSC with:
 *
 * - grid-side equivalent impedance and PCC capacitance,
 * - AC-side filter,
 * - DC-link inductor,
 * - equivalent DC-side resistive load,
 * - PLL,
 * - outer DC-current loop,
 * - inner current loop,
 * - inner capacitor-voltage loop.
 *
 */

 /**
 * @brief Evaluator of the CSC-RL equations and control laws.
 *
 * @details
 * All physical parameters and controller gains are supplied through
 * CSC_RL_Parameters. The model itself stores no internal state.
 *
 */

class CSC_RL_Model
{
public:

    /**
     * @brief External inputs and references evaluated at a given time.
     *
     * @tparam T Scalar type. It can be ODE::real_type or AD::Dual<ODE::real_type>.
     */

    template <typename T> struct Input
    {
        T idc_ref;
        T Q_ref;
    };

    /**
     * @brief Algebraic variables produced by the PLL and control system.
     *
     * @details
     * This structure stores intermediate control variables. These variables are
     * computed from the current state and references, but they are not dynamic
     * states themselves.
     *
     * The PLL frequency estimate is
     *
     * \f[
     * \hat{\omega}
     * =
     * \omega_0
     * +
     * k_{p,\mathrm{pll}} e_{\mathrm{pll}}
     * +
     * k_{i,\mathrm{pll}} \xi_{\mathrm{pll}}.
     * \f]
     *
     * The outer loop generates the active-power reference:
     *
     * \f[
     * P_{\mathrm{ref}}
     * =
     * k_{pO}
     * \left(i_{dc,\mathrm{ref}}^2-i_{stk}^2\right)
     * +
     * k_{iO}\xi_{\mathrm{idc2}}.
     * \f]
     *
     * @tparam T Scalar type used for the algebraic variables.
     */

    template <typename T> struct Control
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

    /**
     * @brief Evaluates stair references from the parameter structure.
     *
     * @param p Complete CSC-RL parameter set.
     * @param t Current simulation time.
     * @return Input<T> DC-current and reactive-power references at time @p t.
     *
     * @tparam T Scalar type used in the returned input structure.
     */

    template <typename T>
    Input<T> input_from_parameters(const CSC_RL_Parameters& p,
                                   ODE::real_type t) const
    {
        Input<T> input;

        input.idc_ref = T(p.idc_ref.value(t));
        input.Q_ref   = T(p.Q_ref.value(t));

        return input;
    }

    /**
     * @brief Computes PLL, outer-loop, current-loop, and voltage-loop variables.
     *
     * @details
     * The PLL error is
     *
     * \f[
     * e_{\mathrm{pll}}=\frac{e_q}{V_{\mathrm{dq,nom}}}.
     * \f]
     *
     * The estimated angular frequency is
     *
     * \f[
     * \hat{\omega}
     * =
     * \omega_0
     * +
     * k_{p,\mathrm{pll}} e_{\mathrm{pll}}
     * +
     * k_{i,\mathrm{pll}} \xi_{\mathrm{pll}}.
     * \f]
     *
     * The outer DC-current loop uses
     *
     * \f[
     * E_{\mathrm{idc2}}
     * =
     * i_{dc,\mathrm{ref}}^2
     * -
     * i_{stk}^2,
     * \f]
     *
     * \f[
     * P_{\mathrm{ref}}
     * =
     * k_{pO}E_{\mathrm{idc2}}
     * +
     * k_{iO}\xi_{\mathrm{idc2}}.
     * \f]
     *
     * The AC current references are obtained from the synchronous-frame power
     * equations:
     *
     * \f[
     * i_d^*
     * =
     * \frac{
     * e_dP_{\mathrm{ref}}+e_qQ_{\mathrm{ref}}
     * }{
     * e_d^2+e_q^2
     * },
     * \f]
     *
     * \f[
     * i_q^*
     * =
     * \frac{
     * e_qP_{\mathrm{ref}}-e_dQ_{\mathrm{ref}}
     * }{
     * e_d^2+e_q^2
     * }.
     * \f]
     *
     * @param p Complete CSC-RL parameter set.
     * @param x Current state vector.
     * @param input External references evaluated at the current time.
     * @return Control<typename Vec::Scalar> Controller algebraic variables.
     *
     */

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

    /**
     * @brief Evaluates the nonlinear right-hand side of the CSC-RL model.
     *
     * @details
     * This function evaluates
     *
     * \f[
     * \frac{dx}{dt}=f(t,x).
     * \f]
     *
     * The grid-side current equations are
     *
     * \f[
     * \frac{d i_{gd}}{dt}
     * =
     * \frac{
     * e_{gd}
     * -
     * R_g i_{gd}
     * -
     * e_d
     * +
     * \hat{\omega}L_g i_{gq}
     * }{L_g},
     * \f]
     *
     * \f[
     * \frac{d i_{gq}}{dt}
     * =
     * \frac{
     * e_{gq}
     * -
     * R_g i_{gq}
     * -
     * e_q
     * -
     * \hat{\omega}L_g i_{gd}
     * }{L_g}.
     * \f]
     *
     * The PCC voltage equations are
     *
     * \f[
     * \frac{d e_d}{dt}
     * =
     * \frac{
     * i_{gd}
     * -
     * i_d
     * +
     * \hat{\omega}C_g e_q
     * }{C_g},
     * \f]
     *
     * \f[
     * \frac{d e_q}{dt}
     * =
     * \frac{
     * i_{gq}
     * -
     * i_q
     * -
     * \hat{\omega}C_g e_d
     * }{C_g}.
     * \f]
     *
     * The PWM-CSC equations are
     *
     * \f[
     * \frac{d i_d}{dt}
     * =
     * \frac{
     * e_d
     * -
     * v_d
     * -
     * R_f i_d
     * +
     * \hat{\omega}L_f i_q
     * }{L_f},
     * \f]
     *
     * \f[
     * \frac{d i_q}{dt}
     * =
     * \frac{
     * e_q
     * -
     * v_q
     * -
     * R_f i_q
     * -
     * \hat{\omega}L_f i_d
     * }{L_f},
     * \f]
     *
     * \f[
     * \frac{d v_d}{dt}
     * =
     * \frac{
     * i_d
     * -
     * m_d i_{stk}
     * +
     * \hat{\omega}C_f v_q
     * }{C_f},
     * \f]
     *
     * \f[
     * \frac{d v_q}{dt}
     * =
     * \frac{
     * i_q
     * -
     * m_q i_{stk}
     * -
     * \hat{\omega}C_f v_d
     * }{C_f},
     * \f]
     *
     * \f[
     * \frac{d i_{stk}}{dt}
     * =
     * \frac{
     * m_dv_d
     * +
     * m_qv_q
     * -
     * R_L i_{stk}
     * }{L_{dc}}.
     * \f]
     *
     * @param p Complete CSC-RL parameter set.
     * @param t Current simulation time.
     * @param x Current state vector.
     * @param dxdt Output vector for state derivatives.
     * @param input References used during this evaluation.
     *
     * @tparam Vec Eigen vector type.
     */

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

    /**
     * @brief Evaluates the right-hand side using references stored in the parameter set.
     *
     * @details
     * This overload first evaluates the stair references at time @p t and then
     * calls the full rhs overload.
     *
     * @param p Complete CSC-RL parameter set.
     * @param t Current simulation time.
     * @param x Current state vector.
     * @param dxdt Output vector for state derivatives.
     *
     * @tparam Vec Eigen vector type.
     */

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

/**
 * @brief Internal namespace for the reduced equilibrium vector indexing.
 *
 * @details
 * The equilibrium vector excludes the PLL angle because
 *
 * \f[
 * \frac{d\hat{\theta}}{dt}=\hat{\omega}\neq 0
 * \f]
 *
 * during synchronous steady-state operation.
 */

namespace CSC_RL_Equilibrium_Detail
{
    enum EquilibriumIndex
    {
        Z_IGD = 0,
        Z_IGQ,
        Z_ED,
        Z_EQ,
        Z_ID,
        Z_IQ,
        Z_VD,
        Z_VQ,
        Z_ISTK,
        Z_XI_PLL,
        Z_XI_UCD,
        Z_XI_UCQ,
        Z_XI_ISD,
        Z_XI_ISQ,
        Z_XI_IDC2,
        Z_N
    };
}

/**
 * @brief Builds the initial guess for the reduced Newton equilibrium problem.
 *
 * @details
 * The initial active power is approximated from the resistive DC-side load:
 *
 * \f[
 * P_0=R_L i_{dc,\mathrm{ref}}^2.
 * \f]
 *
 * Assuming the voltage vector is initially aligned with the d-axis,
 *
 * \f[
 * e_d^0=V_{\mathrm{dq,nom}},
 * \qquad
 * e_q^0=0,
 * \f]
 *
 * the initial current estimates are
 *
 * \f[
 * i_d^0=\frac{P_0}{e_d^0},
 * \qquad
 * i_q^0=-\frac{Q_{\mathrm{ref}}}{e_d^0}.
 * \f]
 *
 * @param p Complete CSC-RL parameter set.
 * @param ref References imposed at the equilibrium point.
 * @return Initial reduced vector for Newton.
 */

inline ODE::vec_type make_csc_rl_equilibrium_initial_guess(const CSC_RL_Parameters& p,
                                                           const EquilibriumReferences& ref)
{
    using namespace CSC_RL_Equilibrium_Detail;

    ODE::vec_type z0;
    z0.resize(Z_N);
    z0.setZero();

    const ODE::real_type ed0 = p.Vdq_nom;
    const ODE::real_type eq0 = 0.0;

    const ODE::real_type P0 = p.RL * ref.idc_ref * ref.idc_ref;

    const ODE::real_type id0 = P0 / ed0;
    const ODE::real_type iq0 = -ref.Q_ref / ed0;

    z0(Z_IGD) = id0;
    z0(Z_IGQ) = iq0;

    z0(Z_ED) = ed0;
    z0(Z_EQ) = eq0;

    z0(Z_ID) = id0;
    z0(Z_IQ) = iq0;

    z0(Z_VD) = ed0;
    z0(Z_VQ) = eq0;

    z0(Z_ISTK) = ref.idc_ref;

    z0(Z_XI_PLL)  = 0.0;
    z0(Z_XI_UCD)  = 0.0;
    z0(Z_XI_UCQ)  = 0.0;
    z0(Z_XI_ISD)  = 0.0;
    z0(Z_XI_ISQ)  = 0.0;
    z0(Z_XI_IDC2) = 0.0;

    return z0;
}

/**
 * @brief Builds the initial guess for the reduced Newton equilibrium problem.
 *
 * @details
 * The initial active power is approximated from the resistive DC-side load:
 *
 * \f[
 * P_0=R_L i_{dc,\mathrm{ref}}^2.
 * \f]
 *
 * Assuming the voltage vector is initially aligned with the d-axis,
 *
 * \f[
 * e_d^0=V_{\mathrm{dq,nom}},
 * \qquad
 * e_q^0=0,
 * \f]
 *
 * the initial current estimates are
 *
 * \f[
 * i_d^0=\frac{P_0}{e_d^0},
 * \qquad
 * i_q^0=-\frac{Q_{\mathrm{ref}}}{e_d^0}.
 * \f]
 *
 * @param p Complete CSC-RL parameter set.
 * @param ref References imposed at the equilibrium point.
 * @return Initial reduced vector for Newton.
 */

class CSC_RL_EquilibriumProblem : public AD::NewtonProblemAD<ODE::real_type>
{
public:
    using real = ODE::real_type;

    using AD::NewtonProblemAD<real>::F;
    using AD::NewtonProblemAD<real>::JF;

    static constexpr int Z_N = CSC_RL_Equilibrium_Detail::Z_N;

    /**
     * @brief Newton problem used to compute a CSC-RL steady-state operating point.
     *
     * @details
     * This class defines the residual
     *
     * \f[
     * F(z)=0,
     * \f]
     *
     * where \f$z\f$ is the reduced state vector. The reduced vector excludes the
     * PLL angle \f$\hat{\theta}\f$ but includes the remaining physical and
     * controller states.
     *
     * The reduced vector is expanded to the full CSC-RL state vector, the same
     * rhs() function used by the RK4 solver is evaluated, and the residual is
     * formed by enforcing zero derivatives for the selected states.
     */

    CSC_RL_EquilibriumProblem(const CSC_RL_Model& model,
                              const CSC_RL_Parameters& p,
                              const EquilibriumReferences& ref)
        : model_(model)
        , p_(p)
        , ref_(ref)
    {
        this->_x0 = make_csc_rl_equilibrium_initial_guess(p_, ref_);
    }

    /**
     * @brief Evaluates the equilibrium residual using dual numbers.
     *
     * @details
     * Because the residual is evaluated with AD::Dual numbers, the base class can
     * construct the Jacobian automatically by injecting one dual perturbation per
     * unknown.
     *
     * @param X Reduced equilibrium vector represented with dual-number scalars.
     * @return Residual vector with dual components carrying derivative information.
     */

    d_dual_vec F(d_dual_vec const& X) const override
    {
        using DualReal = AD::Dual<real>;
        using namespace CSC_RL_Equilibrium_Detail;

        d_dual_vec x;
        d_dual_vec dxdt;
        d_dual_vec res;

        x.resize(NSTATES);
        dxdt.resize(NSTATES);
        res.resize(Z_N);

        x.setZero();
        dxdt.setZero();
        res.setZero();

        // ------------------------------------------------------------
        // Reduced equilibrium vector -> full CSC_RL state vector
        // ------------------------------------------------------------
        x(IGD)  = X(Z_IGD);
        x(IGQ)  = X(Z_IGQ);
        x(ED)   = X(Z_ED);
        x(EQ)   = X(Z_EQ);

        x(ID)   = X(Z_ID);
        x(IQ)   = X(Z_IQ);
        x(VD)   = X(Z_VD);
        x(VQ)   = X(Z_VQ);
        x(ISTK) = X(Z_ISTK);

        // theta_hat is not part of the equilibrium unknowns.
        x(THETA_HAT) = DualReal(p_.x0(THETA_HAT));

        x(XI_PLL)  = X(Z_XI_PLL);
        x(XI_UCD)  = X(Z_XI_UCD);
        x(XI_UCQ)  = X(Z_XI_UCQ);
        x(XI_ISD)  = X(Z_XI_ISD);
        x(XI_ISQ)  = X(Z_XI_ISQ);
        x(XI_IDC2) = X(Z_XI_IDC2);

        // ------------------------------------------------------------
        // Fixed references for equilibrium
        // ------------------------------------------------------------
        typename CSC_RL_Model::template Input<DualReal> input;

        input.idc_ref = DualReal(ref_.idc_ref);
        input.Q_ref   = DualReal(ref_.Q_ref);

        // ------------------------------------------------------------
        // Evaluate the same model used by ODE4
        // ------------------------------------------------------------
        model_.rhs(p_, p_.t0, x, dxdt, input);

        // ------------------------------------------------------------
        // Equilibrium residuals: dx/dt = 0
        // theta_hat is excluded because dtheta_hat/dt = omega_hat.
        // ------------------------------------------------------------
        res(Z_IGD) = dxdt(IGD);
        res(Z_IGQ) = dxdt(IGQ);
        res(Z_ED)  = dxdt(ED);
        res(Z_EQ)  = dxdt(EQ);

        res(Z_ID)   = dxdt(ID);
        res(Z_IQ)   = dxdt(IQ);
        res(Z_VD)   = dxdt(VD);
        res(Z_VQ)   = dxdt(VQ);
        res(Z_ISTK) = dxdt(ISTK);

        res(Z_XI_PLL)  = dxdt(XI_PLL);
        res(Z_XI_UCD)  = dxdt(XI_UCD);
        res(Z_XI_UCQ)  = dxdt(XI_UCQ);
        res(Z_XI_ISD)  = dxdt(XI_ISD);
        res(Z_XI_ISQ)  = dxdt(XI_ISQ);
        res(Z_XI_IDC2) = dxdt(XI_IDC2);

        return res;
    }

private:
    const CSC_RL_Model& model_;
    CSC_RL_Parameters p_;
    EquilibriumReferences ref_;
};

/**
 * @brief Expands the reduced Newton solution into the full state vector.
 *
 * @param z Reduced Newton solution.
 * @param p Parameter set used to recover the fixed PLL angle.
 * @return Full CSC-RL state vector with NSTATES entries.
 */

inline ODE::vec_type expand_csc_rl_equilibrium_to_full_state(const ODE::vec_type& z,
                                                             const CSC_RL_Parameters& p)
{
    using namespace CSC_RL_Equilibrium_Detail;

    if (z.size() != Z_N)
    {
        throw std::runtime_error("CSC_RL equilibrium: invalid z size");
    }

    ODE::vec_type x0;
    x0.resize(NSTATES);
    x0.setZero();

    x0(IGD) = z(Z_IGD);
    x0(IGQ) = z(Z_IGQ);

    x0(ED) = z(Z_ED);
    x0(EQ) = z(Z_EQ);

    x0(ID) = z(Z_ID);
    x0(IQ) = z(Z_IQ);

    x0(VD) = z(Z_VD);
    x0(VQ) = z(Z_VQ);

    x0(ISTK) = z(Z_ISTK);

    x0(THETA_HAT) = p.x0(THETA_HAT);

    x0(XI_PLL)  = z(Z_XI_PLL);
    x0(XI_UCD)  = z(Z_XI_UCD);
    x0(XI_UCQ)  = z(Z_XI_UCQ);
    x0(XI_ISD)  = z(Z_XI_ISD);
    x0(XI_ISQ)  = z(Z_XI_ISQ);
    x0(XI_IDC2) = z(Z_XI_IDC2);

    return x0;
}

/**
 * @brief Stores a full equilibrium state into the parameter structure.
 *
 * @details
 * If Newton converges, this function replaces CSC_RL_Parameters::x0 with the
 * equilibrium point. Therefore, the RK4 simulation starts from a consistent
 * operating point and avoids artificial start-up dynamics.
 *
 * @param p Parameter set to update.
 * @param x0_full Full equilibrium state.
 */

inline void apply_csc_rl_equilibrium_to_parameters(CSC_RL_Parameters& p,
                                                   const ODE::vec_type& x0_full)
{
    if (x0_full.size() != NSTATES)
    {
        throw std::runtime_error("CSC_RL equilibrium: invalid full state size");
    }

    p.x0 = x0_full;
}

inline bool compute_csc_rl_equilibrium(const CSC_RL_Model& model,
                                       CSC_RL_Parameters& p,
                                       const EquilibriumOptions& options)
{
    EquilibriumReferences eq_ref;

    if (options.use_reference_override)
    {
        eq_ref = options.reference;
    }
    else
    {
        eq_ref.idc_ref = p.idc_ref.value(p.t0);
        eq_ref.Q_ref   = p.Q_ref.value(p.t0);
    }

    if (options.verbose)
    {
        std::cout << "Computing CSC_RL equilibrium point...\n";
        std::cout << "idc_ref = " << eq_ref.idc_ref << " A\n";
        std::cout << "Q_ref   = " << eq_ref.Q_ref << " var\n\n";
    }

    CSC_RL_EquilibriumProblem eq_problem(model, p, eq_ref);

    ODE::vec_type z_star;
    const bool ok = solve_equilibrium_problem(eq_problem, options, z_star);

    if (!ok)
    {
        return false;
    }

    const ODE::vec_type x0_eq =
        expand_csc_rl_equilibrium_to_full_state(z_star, p);

    apply_csc_rl_equilibrium_to_parameters(p, x0_eq);

    return true;
}

inline bool compute_csc_rl_equilibrium(const CSC_RL_Model& model,
                                       CSC_RL_Parameters& p,
                                       bool verbose = false)
{
    EquilibriumOptions options;
    options.verbose = verbose;

    return compute_csc_rl_equilibrium(model, p, options);
}

inline bool compute_csc_rl_equilibrium(CSC_RL_Parameters& p,
                                       const EquilibriumOptions& options)
{
    static CSC_RL_Model model;

    return compute_csc_rl_equilibrium(model, p, options);
}

inline bool compute_csc_rl_equilibrium(CSC_RL_Parameters& p,
                                       bool verbose = false)
{
    EquilibriumOptions options;
    options.verbose = verbose;

    return compute_csc_rl_equilibrium(p, options);
}