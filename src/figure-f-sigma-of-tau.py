#!/usr/bin/env python3

####################################################
# filename: figure-f-sigma-of-tau.py
# author: Dylan Morris <dhmorris@princeton.edu>
#
# description: figure showing how f and sigma vary
# with duration tau
####################################################
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from matplotlib.gridspec import GridSpec

import optimize_interventions as oi
import plotting_style as ps
import parameters as params

def setup_figure():
    width = 10
    height = 5

    figure = plt.figure(figsize = (width, height))

    plot_positions = [
        {"name": "f-of-tau",
         "row": 1,
         "col": 1,
         "sharex": None,
         "sharey": None},
        
        {"name": "tau-crash-of-R0",
         "row": 1,
         "col": 2,
         "sharex": None,
         "sharey": None}
    ]
    plots = ps.setup_multipanel(figure,
                                plot_positions,
                                letter_loc = (-0.025, 1.1))
    
    return(figure, plots)



def main(fixed_path,
         outpath):
    fig, plots = setup_figure()

    vline_darkness = 0.8
    vline_lw = mpl.rcParams["lines.linewidth"] * 0.7

    Sfs = [oi.calc_Sf_opt(
        params.R0_default,
        params.gamma_default * tau)
           for tau in params.taus_sweep]

    R0s = np.linspace(1.001, 20, 100)
    tau_crashes = [oi.tau_crash(R0) / params.gamma_default
                   for R0 in R0s]
    fixed_dat = pd.read_csv(fixed_path)
    plots["f-of-tau"].axvline(
        oi.tau_crash(params.R0_default) / params.gamma_default,
        color = ps.full_suppression_cmap(
            vline_darkness),
        lw = vline_lw,
        linestyle = "dotted")
    
    plots["f-of-tau"].axvline(
        fixed_dat[fixed_dat["R_e_i"] > 1]["tau"].min(),
        color = ps.fixed_cmap(
            vline_darkness),
        lw = vline_lw,
        linestyle = "dotted")
    
    plots["f-of-tau"].plot(
        params.taus_sweep,
        [item[1] for item in Sfs],
        color = ps.opt_cmap(0.8))

    plots["f-of-tau"].plot(
        fixed_dat["tau"],
        fixed_dat["sigma"],
        color = ps.fixed_cmap(0.8))

    plots["tau-crash-of-R0"].plot(
        R0s,
        tau_crashes,
        color = ps.opt_cmap(0.8))
    
    plots["f-of-tau"].set_ylim(bottom = -0.1, top = 1.1)
    plots["f-of-tau"].set_yticks([0, 0.25, 0.5, 0.75, 1])
    plots["f-of-tau"].set_xticks(np.arange(0, 90, 20))
    plots["f-of-tau"].set_ylabel("$f, \sigma$")
    plots["f-of-tau"].set_xlabel("$\\tau$")
    plots["tau-crash-of-R0"].set_xlabel("$\mathcal{R}_0$")

    plots["tau-crash-of-R0"].set_ylabel("$\\tau_{\mathrm{crit}}$")
    plots["tau-crash-of-R0"].set_ylim(bottom = -5, top = 55)
    plots["tau-crash-of-R0"].set_yticks(np.arange(0, 60, 15))
    plots["tau-crash-of-R0"].set_xticks(np.arange(0, 22, 5))
    for plot in plots.values():
        plot.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))


    fig.tight_layout()
    fig.savefig(outpath)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("USAGE ./{} "
              "<path to fixed strategies> "
              "<outpath>\n".format(sys.argv[0]))
    else:
        main(sys.argv[1],
             sys.argv[2])
