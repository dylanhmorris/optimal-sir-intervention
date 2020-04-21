#!/usr/bin/env python3

####################################################
# filename: figure-interventions.py
# author: Dylan Morris <dhmorris@princeton.edu>
#
# description: plots figure showing example
# interventions, how they are achieved, and
# their effects
####################################################
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib as mpl

import matplotlib.lines as mlines
from matplotlib.ticker import FormatStrFormatter

import InterventionSIR as sir
import plotting_style as ps
import parameters as params
import optimize_interventions as oi

np.random.seed(8544)

fig_cols = 3
fig_rows = 4
width = 17
height = 8 

def setup_figure():    
    figure = plt.figure(
        figsize = (width, height))

    opt_col, fixed_col, full_suppression_col = 1, 2, 3
    ti_col, imax_col = 4, 4
    imax_row, ti_row= 1, 2
    timecourse_row, b_row = 1, 2

    
    plot_positions = [
        {"name": "mc-time-timecourse",
         "row": timecourse_row,
         "col": opt_col,
         "sharex": None,
         "sharey": None},
        
        {"name": "fixed-timecourse",
         "row": timecourse_row,
         "col": fixed_col,
         "sharex": "mc-time-timecourse",
         "sharey": "mc-time-timecourse"},
        
        {"name": "full-suppression-timecourse",
         "row": timecourse_row,
         "col": full_suppression_col,
         "sharex": "mc-time-timecourse",
         "sharey": "mc-time-timecourse"},

        {"name": "mc-time-b",
         "row": b_row,
         "col": opt_col,
         "sharex": "mc-time-timecourse",
         "sharey": None},
        
        {"name": "fixed-b",
         "row": b_row,
         "col": fixed_col,
         "sharex": "mc-time-b",
         "sharey": "mc-time-b"},
        
        {"name": "full-suppression-b",
         "row": b_row,
         "col": full_suppression_col,
         "sharex": "mc-time-b",
         "sharey": "mc-time-b"},

        {"name": "imax",
         "row": imax_row,
         "col": imax_col,
         "sharex": None,
         "sharey": "mc-time-timecourse"},

        {"name": "ti",
         "row": ti_row,
         "col": ti_col,
         "sharex": "imax",
         "sharey": None}
    ]

    plots = ps.setup_multipanel(figure, plot_positions)
    
    return figure, plots


def main(path_to_fixed,
         outpath):
    
    ## read in data
    fixed = pd.read_csv(path_to_fixed)

    ## set up numerics
    taus = params.taus_figure_interventions
    
    covid_sir = sir.InterventionSIR(
        b_func = sir.Intervention(),
        R0 = params.R0_default,
        gamma = params.gamma_default,
        inits = params.inits_default)
    covid_sir.reset()

    t_sim_max = 360
    t_plot_max = 300

    ## set up figure
    fig, plots = setup_figure()
    null_color = "black"
    null_style = "dashed"
    null_dashes = (3, 3)
    timecourse_alpha = 0.9
    tau_plot_darkness = 0.8
    tau_crash_lw = mpl.rcParams["lines.linewidth"] * 0.7
    
    ## set up legend
    color_vals = np.linspace(0.8, 0.3, len(taus))
    leg_cmap = plt.cm.Greys
    handlelines = []
    for color_val, tau in zip(color_vals, taus): 
        handlelines.append(
            mlines.Line2D([], [],
                          color = leg_cmap(color_val),
                          label = tau))

    ## calculate and plot non intervention results
    null_time, null_result = covid_sir.integrate_null(t_sim_max)

    ## plot tau crash below plots of tau
    tau_crash = oi.tau_crash(covid_sir.R0) / covid_sir.gamma
    plots["imax"].axvline(tau_crash,
                          color = ps.full_suppression_cmap(
                              tau_plot_darkness),
                          lw = tau_crash_lw,
                          linestyle = "dotted")
    
    plots["ti"].axvline(tau_crash,
                        color = ps.full_suppression_cmap(
                            tau_plot_darkness),
                        lw = tau_crash_lw,
                        linestyle = "dotted")

    for plotname in ["mc-time-timecourse",
                     "fixed-timecourse",
                     "full-suppression-timecourse"]:
        plots[plotname].plot(null_time,
                             null_result[:, 1],
                             color = null_color,
                             ls = null_style,
                             dashes = null_dashes)

    ## calculate and plot intervention results
    for name, cmap in zip(
            ["mc-time", "fixed", "full-suppression"],
            [ps.opt_cmap, ps.fixed_cmap, ps.full_suppression_cmap]):
        solids = [] # save solid b lines for plotting later

        ## set intervention strategy
        covid_sir.b_func.strategy = name

        ## iterate over intervention durations
        for i_tau, tau in enumerate(taus):
            covid_sir.b_func.tau = tau
            color_val = color_vals[i_tau]
            S_i_expected = 0
            print("optimizing strategy for {} "
                  "with tau = {}".format(name, tau))
            if name == "mc-time":
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
                
            elif name == "full-suppression":
                S_i_expected = oi.calc_S_var_opt(
                    covid_sir.R0,
                    covid_sir.gamma * tau,
                    0)
                covid_sir.b_func.sigma = 0
            
            t_i_opt = covid_sir.t_of_S(S_i_expected)[0]
            covid_sir.b_func.t_i = t_i_opt
            
            covid_sir.reset()
            covid_sir.integrate(t_sim_max)

            plotname_a = name + "-timecourse"
            plots[plotname_a].plot(
                covid_sir.time_ts,
                covid_sir.state_ts[:, 1],
                color = cmap(color_val),
                alpha = timecourse_alpha)
            
            plotname_b = name + "-b"
            times = np.array(covid_sir.time_ts)

            b_of_t =  np.array([
                covid_sir.b_func(time,
                                 covid_sir.R0 * covid_sir.gamma,
                                 covid_sir.gamma,
                                 S,
                                 I)
                for time, S, I in zip(
                        times,
                        covid_sir.state_ts[:, 0],
                        covid_sir.state_ts[:, 1])]).astype("float")

            disc = ((np.abs(np.diff(b_of_t, append = [1])) > 0.01) |
                    (np.abs(np.diff(b_of_t, prepend = [1])) > 0.01))
            not_disc = np.logical_not(disc)
            sol_y = np.array(b_of_t)
            dash_y = np.array(b_of_t)
            dash_y[not_disc] = np.nan
            sol_y[disc] = np.nan
            solids.append(sol_y)
            # save to plot later,
            # so that solid lines pass over
            # dashed lines
            
            plots[plotname_b].plot(
                times,
                dash_y,
                color = cmap(color_val),
                linestyle = "dashed",
                lw = height / 4,
                dashes = (2, 2),
                alpha = timecourse_alpha)
            
        for i_tau, tau in enumerate(taus):
            color_val = color_vals[i_tau]
            sol_y = solids[i_tau]
            plotname_b = name + "-b"
            times = np.array(covid_sir.time_ts)
            plots[plotname_b].plot(
                times,
                sol_y,
                color = cmap(color_val),
                linestyle = "solid",
                solid_capstyle = "butt",
                alpha = timecourse_alpha)


        ## plot optima as a function of tau
        taus_sweep = params.taus_sweep
        print("calculating and plotting t_i and I_max "
              "for strategy {}...\n".format(name))
        if name == "mc-time":
           Sfs = [oi.calc_Sf_opt(
               covid_sir.R0,
               covid_sir.gamma * tau)
                  for tau in taus_sweep]
           t_is = [covid_sir.t_of_S(Sf[0])
                   for Sf in Sfs]
           Imaxes = [oi.I_of_S(S, covid_sir.R0) for S, f in Sfs]
           taus_plot = taus_sweep
        elif name == "fixed":
            t_is = fixed["t_i"]
            Imaxes = fixed["Imax"]
            taus_plot = fixed["tau"]
                       
        elif name == "full-suppression":
            S_is = [oi.calc_S_var_opt(
                covid_sir.R0,
                covid_sir.gamma * tau,
                0)
                    for tau in taus_sweep]

            t_is = [covid_sir.t_of_S(S)
                    for S in S_is]
            Imaxes = [oi.I_max_opt_of_S_i(S, 0,
                                         covid_sir.R0,
                                         covid_sir.gamma * tau)
                      for S, tau in zip(S_is, taus_sweep)]
            taus_plot = taus_sweep
            
        plotname_Imax = "imax"
        plotname_t_i =  "ti"
        plots[plotname_Imax].plot(taus_plot,
                                  Imaxes,
                                  color = cmap(tau_plot_darkness),
                                  alpha = timecourse_alpha)
        plots[plotname_t_i].plot(taus_plot,
                                 t_is,
                                 color = cmap(tau_plot_darkness),
                                 alpha = timecourse_alpha)

    ## labeling and styling
    # ylim for 0 to 1 plots
    ymin = -0.1
    ymax = 1.1
    plots["mc-time-timecourse"].legend(title = "$\\tau$ (days)",
                                       handles = handlelines)

    plots["mc-time-timecourse"].set_ylabel("prevalence $I(t)$")
    plots["mc-time-b"].set_ylabel("intervention $b(t)$")
    plots["imax"].set_ylabel("peak $I^{\max}$",
                             rotation = 270,
                             labelpad = 30)
    plots["ti"].set_xlabel("duration $\\tau$ (days)")
    plots["ti"].set_ylabel("timing $t^{\mathrm{opt}}_i$",
                           rotation = 270,
                           labelpad = 30)
    plots["fixed-b"].set_xlabel("time (days)")

    plots["mc-time-timecourse"].set_xlim(0, 300)
    plots["mc-time-timecourse"].set_ylim(ps.ymin_timecourse,
                                         ps.ymax_timecourse)
    plots["mc-time-timecourse"].set_yticks(np.arange(0, 0.4, 0.1))
    plots["mc-time-timecourse"].set_xticks(np.arange(0, 280,
                                                     90))
    plots["mc-time-b"].set_yticks(np.arange(0, 1.1, 0.25))
    plots["mc-time-b"].set_ylim(ymin, ymax)
    
    plots["imax"].set_xticks(np.arange(0, 120, 30))
    plots["imax"].set_xlim(0, 90)

    plots["ti"].set_yticks(np.arange(75, 100, 5))
    plots["ti"].set_ylim(73, 97)


    ## add some titles
    plots["mc-time-timecourse"].set_title("Optimal")
    plots["fixed-timecourse"].set_title("Fixed control")
    plots["full-suppression-timecourse"].set_title("Full suppression")
    plots["full-suppression-timecourse"].set_title("Full suppression")
    plots["imax"].set_title("Effect of $\\tau$")
    
    for plotname, plot in plots.items():
        if not plot.is_first_col() and not plotname in ["ti", "imax"]:
            plt.setp(plot.get_yticklabels(), visible = False)
        else:
            plot.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        if plotname in ["ti", "imax"]:
            plot.yaxis.set_label_position("right")
            plot.yaxis.tick_right()

    fig.tight_layout()
    fig.savefig(outpath)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("USAGE ./{} <path to fixed b results> "
              "<outpath>\n".format(sys.argv[0]))
    else:
        main(sys.argv[1],
             sys.argv[2])
