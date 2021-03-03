#!/usr/bin/env python3

####################################################
# filename: figure-f-sigma-of-tau.py
# author: Dylan Morris <dhmorris@princeton.edu>
#
# description: figure showing how f and sigma vary
# with duration tau, and example timecourses
####################################################
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from matplotlib.gridspec import GridSpec

import optimize_interventions as oi
import InterventionSIR as sir
import plotting_style as ps
import parameters as params

def setup_figure():
    width = 10
    height = 10

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
         "sharey": None},

        {"name": "maintain-suppress-time-timecourse",
         "row": 2,
         "col": 1,
         "sharex": None,
         "sharey": None},

        {"name": "fixed-timecourse",
         "row": 2,
         "col": 2,
         "sharex": "maintain-suppress-time-timecourse",
         "sharey": "maintain-suppress-time-timecourse"}
    ]
    
    plots = ps.setup_multipanel(figure,
                                plot_positions,
                                letter_loc = (-0.025, 1.1))
    
    return(figure, plots)



def main(fixed_path,
         outpath):
    fixed_dat = pd.read_csv(fixed_path)

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


    tau_crit = fixed_dat[fixed_dat["R_e_i"] > 1]["tau"].min()

    taus = [tau_crit, 80, 160]
    tau_crit_lab = "$\\tau_1 \\approx$ {:.0f}".format(tau_crit)
    print(tau_crit_lab)
    tau_labels = [tau_crit_lab] + taus[1:]
    
    plots["f-of-tau"].axvline(
        oi.tau_crash(params.R0_default) / params.gamma_default,
        color = ps.full_suppression_cmap(
            vline_darkness),
        lw = vline_lw,
        linestyle = "dotted")
    
    plots["f-of-tau"].axvline(
        tau_crit,
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


    covid_sir = sir.InterventionSIR(
        b_func = sir.Intervention(),
        R0 = params.R0_default,
        gamma = params.gamma_default,
        inits = params.inits_default)
    covid_sir.reset()

    t_sim_max = 360
    t_plot_max = 300
    null_color = "black"
    null_style = "dashed"
    null_dashes = (3, 3)
    timecourse_alpha = 0.9
    
    color_vals = np.linspace(0.8, 0.3, len(taus))
    leg_cmap = plt.cm.Greys
    handlelines = []
    for color_val, tau, tau_label in zip(color_vals, taus, tau_labels): 
        handlelines.append(
            mpl.lines.Line2D([], [],
                             color = leg_cmap(color_val),
                             label = tau_label))

    ## calculate and plot non intervention results
    null_time, null_result = covid_sir.integrate_null(t_sim_max)

    ## calculate and plot intervention results
    for name, cmap in zip(
            ["maintain-suppress-time", "fixed"],
            [ps.opt_cmap, ps.fixed_cmap]):


        ## set intervention strategy
        covid_sir.b_func.strategy = name
        plotname = name + "-timecourse"
        plots[plotname].plot(null_time,
                             null_result[:, 1],
                             color = null_color,
                             ls = null_style,
                             dashes = null_dashes)

        ## iterate over intervention durations
        for i_tau, tau in enumerate(taus):
            covid_sir.b_func.tau = tau
            color_val = color_vals[i_tau]
            S_i_expected = 0
            print("optimizing strategy for {} "
                  "with tau = {}".format(name, tau))
            if name == "maintain-suppress-time":
                S_i_expected, f = oi.calc_Sf_opt(
                    covid_sir.R0,
                    covid_sir.gamma * tau)
                I_i_expected = covid_sir.I_of_S(S_i_expected)
                covid_sir.b_func.S_i_expected = S_i_expected
                covid_sir.b_func.I_i_expected = I_i_expected
                covid_sir.b_func.f = f

            elif name == "fixed":
                S_i_expected, sigma = oi.calc_Sb_opt(
                    covid_sir.R0,
                    covid_sir.gamma,
                    tau)
                covid_sir.b_func.sigma = sigma
                pass
            
            t_i_opt = covid_sir.t_of_S(S_i_expected)[0]
            covid_sir.b_func.t_i = t_i_opt
            
            covid_sir.reset()
            covid_sir.integrate(t_sim_max)

            plots[plotname].plot(
                covid_sir.time_ts,
                covid_sir.state_ts[:, 1],
                color = cmap(color_val),
                alpha = timecourse_alpha)
            pass
        pass
    pass


    plots["maintain-suppress-time-timecourse"].set_xlim(0, 300)
    plots["maintain-suppress-time-timecourse"].set_ylim(ps.ymin_timecourse,
                                         ps.ymax_timecourse)
    plots["maintain-suppress-time-timecourse"].set_yticks(np.arange(0, 0.4, 0.1))
    plots["maintain-suppress-time-timecourse"].set_xticks(np.arange(0, 280,
                                                     90))

    plots["f-of-tau"].set_ylim(bottom = -0.1, top = 1.1)
    plots["f-of-tau"].set_yticks([0, 0.25, 0.5, 0.75, 1])
    plots["f-of-tau"].set_xticks(np.arange(0, 90, 20))
    plots["f-of-tau"].set_ylabel("$f, \sigma$")
    plots["f-of-tau"].set_xlabel("$\\tau$ (days)")
    plots["tau-crash-of-R0"].set_xlabel("$\mathcal{R}_0$")

    plots["tau-crash-of-R0"].set_ylabel("$\\tau_{\mathrm{crit}}$")
    plots["tau-crash-of-R0"].set_ylim(bottom = -5, top = 55)
    plots["tau-crash-of-R0"].set_yticks(np.arange(0, 60, 15))
    plots["tau-crash-of-R0"].set_xticks(np.arange(0, 22, 5))
    for plot in plots.values():
        plot.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        
    plots["maintain-suppress-time-timecourse"].legend(title = "$\\tau$ (days)",
                                       handles = handlelines)
    plots["maintain-suppress-time-timecourse"].set_ylabel("prevalence $I(t)$")


    plots["maintain-suppress-time-timecourse"].set_title("Optimal")
    plots["fixed-timecourse"].set_title("Fixed control")
    plots["maintain-suppress-time-timecourse"].set_xlabel("time (days)")
    plots["fixed-timecourse"].set_xlabel("time (days)")

    fig.tight_layout()
    fig.align_ylabels()
    fig.savefig(outpath)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("USAGE ./{} "
              "<path to fixed strategies> "
              "<outpath>\n".format(sys.argv[0]))
    else:
        main(sys.argv[1],
             sys.argv[2])
