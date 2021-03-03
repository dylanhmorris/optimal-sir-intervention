#!/usr/bin/env python3

####################################################
# filename: early_late_parameter_sweep.py
# author: Dylan Morris <dhmorris@princeton.edu>
#
# description: calculate parameter sweeps
# for cost of being early vs. late with R0, tau, gamma
####################################################
import sys
import numpy as np
import pandas as pd

import parameters as params
import InterventionSIR as sir
import optimize_interventions as oi

covid_sir = sir.InterventionSIR(
    b_func = None,
    R0 = params.R0_default,
    gamma = params.gamma_default,
    inits = params.inits_sweep_small)
covid_sir.reset()
tau = params.tau_default

t_sim_max = 1000


def sweep_tau_early_late(
        taus,
        R0s,
        offset,
        sir_model,
        strategy = "maintain-suppress-time",
        t_sim_max = 1000):

    sir_model.inits = params.inits_sweep_small
    
    results = np.zeros(len(taus) * len(R0s)).reshape(
        len(R0s), -1)

    sir_model.b_func = sir.Intervention(
        strategy = strategy)
    
    for col, tau in enumerate(taus):
        print(col)
        for row, R0 in enumerate(R0s):
            sir_model.b_func.tau = tau
            sir_model.R0 = R0
            S_i_expected = 0
            
            if strategy == "maintain-suppress-time":
                S_i_expected, f = oi.calc_Sf_opt(
                    sir_model.R0,
                    sir_model.gamma * tau)
                I_i_expected = sir_model.I_of_S(S_i_expected)
                sir_model.b_func.S_i_expected = S_i_expected
                sir_model.b_func.I_i_expected = I_i_expected
                sir_model.b_func.f = f
            elif strategy == "fixed":
                S_i_expected, sigma = oi.calc_Sb_opt(
                    sir_model.R0,
                    sir_model.gamma,
                    tau)
                sir_model.b_func.sigma = sigma
            elif strategy == "full-suppression":
                S_i_expected = oi.calc_S_var_opt(
                    sir_model.R0,
                    sir_model.gamma * tau,
                    0)
                sir_model.b_func.sigma = 0

            t_i_opt = sir_model.t_of_S(S_i_expected)[0]

            ## numerically solve the sir
            ## with our b function
            ## starting too early 
            sir_model.b_func.t_i = t_i_opt - offset
            sir_model.reset()
            sir_model.integrate(t_sim_max)
            early = sir_model.get_I_max()
            
            ## numerically solve the sir
            ## with our b function
            ## starting too late 
            sir_model.b_func.t_i = t_i_opt + offset
            sir_model.reset()
            sir_model.integrate(t_sim_max)
            late = sir_model.get_I_max()
            results[row, col] = late - early
        pass
    return results


def sweep_gamma_early_late(
        gammas,
        R0s,
        offset,
        tau,
        sir_model,
        strategy = "maintain-suppress-time",
        t_sim_normed = 500 / 14):

    results = np.zeros(len(gammas) * len(R0s)).reshape(
        len(R0s), -1)
    sir_model.inits = params.inits_sweep_small

    if not sir.check_gamma(np.max(gammas),
                           sir_model,
                           -offset,
                           strategy,
                           tau = tau,
                           verbose = True):
        raise ValueError("t_i + most negative offset "
                         "is negative!")
    
    sir_model.b_func = sir.Intervention(
        strategy = strategy,
        tau = tau)

    for col, gamma in enumerate(gammas):
        print(col)
        for row, R0 in enumerate(R0s):
            sir_model.gamma = gamma
            sir_model.R0 = R0
            S_i_expected = 0
        
            if strategy == "maintain-suppress-time":
                S_i_expected, f = oi.calc_Sf_opt(
                    sir_model.R0,
                    sir_model.gamma * sir_model.b_func.tau)
                I_i_expected = sir_model.I_of_S(S_i_expected)
                sir_model.b_func.S_i_expected = S_i_expected
                sir_model.b_func.I_i_expected = I_i_expected
                sir_model.b_func.f = f
            elif strategy == "fixed":
                S_i_expected, sigma = oi.calc_Sb_opt(
                    sir_model.R0,
                    sir_model.gamma,
                    sir_model.b_func.tau)
                sir_model.b_func.sigma = sigma
            elif strategy == "full-suppression":
                S_i_expected = oi.calc_S_var_opt(
                    sir_model.R0,
                    sir_model.gamma * sir_model.b_func.tau,
                    0)
                sir_model.b_func.sigma = 0

            t_i_opt = sir_model.t_of_S(S_i_expected)[0]

            ## numerically solve the sir
            ## with our b function
            ## starting too early 
            sir_model.b_func.t_i = t_i_opt - offset
            sir_model.reset()
            sir_model.integrate(t_sim_max)
            early = sir_model.get_I_max()
            
            ## numerically solve the sir
            ## with our b function
            ## starting too late 
            sir_model.b_func.t_i = t_i_opt + offset
            sir_model.reset()
            sir_model.integrate(t_sim_max)
            late = sir_model.get_I_max()
            results[row, col] = late - early
        pass
    return results

def main(var, strategy, outpath):

    with open(outpath, 'w') as out:
        ## overwrite file, so make knows
        ## we're working on it
        pass

    taus = params.heatmap_taus
    R0s = params.heatmap_R0s
    gammas = params.heatmap_gammas

    if strategy == "fixed":
        taus = params.fixed_heatmap_taus
        R0s = params.fixed_heatmap_R0s
        gammas = params.fixed_heatmap_gammas
    offset = 7

    print("calculating early/late parameter sweep for "
          "R0 and {}...\n".format(var))
    if var == "tau":
        result = pd.DataFrame(
            sweep_tau_early_late(
                taus,
                R0s,
                offset,
                covid_sir,
                strategy = strategy))
        result.columns = taus
    elif var == "gamma":
        result = pd.DataFrame(
            sweep_gamma_early_late(
                gammas,
                R0s,
                offset,
                tau,
                covid_sir,
                strategy = strategy))
        result.columns = gammas
    print("saving results to "
          "{}\n".format(outpath))
    result.index = R0s
    result.to_csv(outpath)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("USAGE ./{} <variable to sweep> "
              "<strategy to use> <outpath>"
              "\n".format(sys.argv[0]))
    else:
        main(sys.argv[1],
             sys.argv[2],
             sys.argv[3])

