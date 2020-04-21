#!/usr/bin/env python3

####################################################
# filename: calc_fixed_optima.py
# author: Dylan Morris <dhmorris@princeton.edu>
#
# description: calculate optimal interventions
# for the fixed strategy over a range of taus
# and save to a file (fixed is slower than others
# because it requires numerical integration)
####################################################

import sys
import parameters as params
import pandas as pd
import optimize_interventions as oi 
import csv
import numpy as np

def main(outpath = None, taus = None,
         R0 = params.R0_default,
         gamma = params.gamma_default,
         I0 = params.inits_default[1],
         Rec0 = params.inits_default[2]):
    np.random.seed(8544)

    with open(outpath, "w") as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        header_row = [
            "tau",
            "t_i",
            "sigma",
            "Imax",
            "S_i",
            "R_e_i"]
        writer.writerow(header_row)
    
    for tau in taus:
        with open(outpath, "a") as csvfile:
            writer = csv.writer(csvfile, delimiter=',')

            S_i, sigma = oi.calc_Sb_opt(
                R0,
                gamma,
                tau,
                n_max_tries = 15)
            t_i = oi.t_of_S(S_i, R0, gamma,
                            I0 = I0,
                            Rec0 = Rec0)[0]
            
            R_e_i = R0 * S_i * sigma

            Imax = oi.Imax_of_S_i_b(
                S_i,
                sigma,
                R0,
                gamma,
                tau)

            writer.writerow([tau, t_i, sigma, Imax, S_i, R_e_i])


if __name__ == "__main__":
    if len(sys.argv) < 1:
        print("USAGE: ./{} <outpath>".format(sys.argv[0]))
    else:
        main(outpath = sys.argv[1],
             taus = params.taus_sweep)
