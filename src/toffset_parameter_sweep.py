#!/usr/bin/env python3

####################################################
# filename: toffset_parameter_sweep.py
# author: Dylan Morris <dhmorris@princeton.edu>
#
# description: calculate parameter sweeps
# for toffset with tau, gamma
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

def sweep_tau_offset(
        taus,
        offsets,
        sir_model,
        strategy = "mc-time",
        t_sim_max = 1000):

    results = np.zeros(len(taus) * len(offsets)).reshape(
        len(offsets), -1)

    sir_model.b_func = sir.Intervention(
        strategy = strategy)
    
    for col, tau in enumerate(taus):
        sir_model.b_func.tau = tau
        S_i_expected = 0
        
        if strategy == "mc-time":
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

        for row, offset_val in enumerate(offsets):
            sir_model.b_func.t_i = t_i_opt + offset_val

            ## numerically solve the sir
            ## with our mistimed b function
            sir_model.reset()
            sir_model.integrate(t_sim_max)
            results[row, col] = sir_model.get_I_max()
            
    return results


def sweep_gamma_offset(
        gammas,
        offsets,
        tau,
        sir_model,
        strategy = "mc-time",
        t_sim_normed = 400 / 14):

    results = np.zeros(len(gammas) * len(offsets)).reshape(
        len(offsets), -1)

    if not sir.check_gamma(np.max(gammas),
                           sir_model,
                           np.min(offsets),
                           strategy,
                           tau = tau,
                           verbose = True):
        raise ValueError("t_i + most negative offset "
                         "is negative!")
    
    sir_model.b_func = sir.Intervention(
        strategy = strategy,
        tau = tau)

    for col, gamma in enumerate(gammas):
        sir_model.gamma = gamma
        S_i_expected = 0
        
        if strategy == "mc-time":
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

        ## calculate intended intervention time
        t_i_opt = sir_model.t_of_S(S_i_expected)[0]
        
        for row, offset_val in enumerate(offsets):
            sir_model.b_func.t_i = t_i_opt + offset_val

            ## solve the sir model
            ## with our intervention b_func
            sir_model.reset()
            sir_model.integrate(t_sim_normed / gamma)
            results[row, col] = sir_model.get_I_max()
    return results

def main(var, strategy, outpath):
    with open(outpath, 'w') as out:
        ## overwrite file, so make knows
        ## we're working on it
        pass
    
    taus = params.heatmap_taus
    gammas = params.heatmap_gammas
    offsets = params.heatmap_offsets

    print("calculating parameter sweep for "
          "t_offset and {}...\n".format(var))
    if var == "tau":
        result = pd.DataFrame(
            sweep_tau_offset(
                taus,
                offsets,
                covid_sir,
                strategy = strategy))
        result.columns = taus
    elif var == "gamma":
        result = pd.DataFrame(
            sweep_gamma_offset(
                gammas,
                offsets,
                tau,
                covid_sir,
                strategy = strategy))
        result.columns = gammas
    print("saving results to "
          "{}\n".format(outpath))
    result.index = offsets
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

