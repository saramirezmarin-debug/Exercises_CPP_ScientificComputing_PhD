#pragma once

#include "CSC_RL.hh"

// ============================================================
// Declaration of the default CSC-RL case factory.
//
// The function declared here returns a fully initialized
// CSC_RL_Parameters structure containing physical parameters,
// controller gains, references, and initial conditions.
// ============================================================

/**
 * @file CSC_RL_Case.hh
 * @brief Declaration of the default CSC-RL simulation case.
 *
 * @details
 * The function declared here returns a fully initialized
 * CSC_RL_Parameters structure containing nominal electrical parameters,
 * controller gains, references, and initial conditions.
 */

 /**
 * @brief Builds the default CSC-RL parameter set.
 *
 * @return Complete CSC_RL_Parameters structure ready to be modified or passed
 *         directly to the CSC_RL ODE wrapper.
 */

CSC_RL_Parameters make_csc_rl_case(); // Build and return the default parameter set for the CSC-RL model.