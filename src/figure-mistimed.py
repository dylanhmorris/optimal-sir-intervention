#!/usr/bin/env python3

####################################################
# filename: figure-mistimed.py
# author: Dylan Morris <dhmorris@princeton.edu>
#
# description: plots figure showing consequences
# of mistiming an optimized intervention
####################################################
import sys
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib as mpl
import numpy as np

import parameters as params
import plotting_style as ps
import InterventionSIR as sir
import optimize_interventions as oi

fig_cols = 3
fig_rows = 3
width = 15
height = 12


tau = params.tau_default
offsets_timecourse = params.offsets_default
sc_reduction = params.sc_reduction

cmaps = ps.cmap_dict 

covid_sir = sir.InterventionSIR(
    b_func = sir.Intervention(),
    R0 = params.R0_default,
    gamma = params.gamma_default,
    inits = params.inits_default)
covid_sir.reset()

covid_sc_sir = sir.InterventionSIR(
    b_func = sir.Intervention(),
    R0 = params.R0_default * sc_reduction,
    gamma = params.gamma_default,
    inits = params.inits_default)
covid_sir.reset()


t_sim_max = 300
t_plot_max = 270

null_time, null_result = covid_sir.integrate_null(t_sim_max)
null_alone_peak =np.max(null_result[:, 1])

null_sc_time, null_sc_result = covid_sc_sir.integrate_null(t_sim_max)
null_sc_peak = np.max(null_sc_result[:, 1])
print(null_sc_peak)

def setup_figure():    
    figure = plt.figure(figsize = (width, height))

    opt_row, fixed_row, full_suppression_row = 1, 2, 3
    alone_col, sc_col, vs_sc_col = 1, 2, 3 

    
    
    plot_positions = [
        {"name": "maintain-suppress-time-alone",
         "row": opt_row,
         "col": alone_col,
         "sharex": None,
         "sharey": None},
        
        {"name": "fixed-alone",
         "row": fixed_row,
         "col": alone_col,
         "sharex": "maintain-suppress-time-alone",
         "sharey": "maintain-suppress-time-alone"},
        
        {"name": "full-suppression-alone",
         "row": full_suppression_row,
         "col": alone_col,
         "sharex": "maintain-suppress-time-alone",
         "sharey": "maintain-suppress-time-alone"},

        {"name": "maintain-suppress-time-sc",
         "row": opt_row,
         "col": sc_col,
         "sharex": "maintain-suppress-time-alone",
         "sharey": "maintain-suppress-time-alone"},
        
        {"name": "fixed-sc",
         "row": fixed_row,
         "col": sc_col,
         "sharex": "maintain-suppress-time-sc",
         "sharey": "maintain-suppress-time-sc"},
        
        {"name": "full-suppression-sc",
         "row": full_suppression_row,
         "col": sc_col,
         "sharex": "maintain-suppress-time-sc",
         "sharey": "maintain-suppress-time-sc"},

        {"name": "maintain-suppress-time-vs-sc",
         "row": opt_row,
         "col": vs_sc_col,
         "sharex": None,
         "sharey": "maintain-suppress-time-alone"},
        
        {"name": "fixed-vs-sc",
         "row": fixed_row,
         "col": vs_sc_col,
         "sharex": "maintain-suppress-time-vs-sc",
         "sharey": "maintain-suppress-time-vs-sc"},
        
        {"name": "full-suppression-vs-sc",
         "row": full_suppression_row,
         "col": vs_sc_col,
         "sharex": "maintain-suppress-time-vs-sc",
         "sharey": "maintain-suppress-time-vs-sc"}
    ]

    plots = ps.setup_multipanel(figure,
                                plot_positions,
                                letter_loc = (-0.025, 1.1))
    
    return figure, plots


def calc_offset_sweep(
        tau,
        sir_model,
        t_sim_max,
        offset = params.offsets_for_sweep,
        offsets = None,
        strategy = "maintain-suppress-time",
        fineness = 25):
    results = []

    sir_model.b_func = sir.Intervention(
        strategy = strategy,
        tau = tau)
    
    if offsets is None:
        offsets = np.linspace(-offset,
                              offset,
                              fineness)
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

    print("calculating I_max given offsets...\n")
    for offset_val in offsets:
        sir_model.b_func.t_i = t_i_opt + offset_val
        sir_model.reset()
        sir_model.integrate(t_sim_max)
        results.append(sir_model.get_I_max())
    return (offsets, np.array(results))

           

def main(outpath,
         offsets_timecourse = offsets_timecourse):
    fig, plots = setup_figure()
    null_color = "black"
    null_sc_color = "gray"
    null_style = "dashed"
    null_dashes = (2.5, 2)
    timecourse_alpha = 0.9
    offset_darkness = 0.8
    handlelines = []
    sc_handlelines = []
    leg_cmap = plt.cm.Greys
    offset_labels = ["week late", "week early", "optimal time"]

    sc_style = "solid"
    
    ## set up legend
    color_vals = {
        "week early": 0.25,
        "optimal time": 0.5,
        "week late": 0.75}

    handlelines.append(
        mlines.Line2D([], [],
                      ls = null_style,
                      dashes = (0.75, 0.75),
                      color = null_color,
                      label = "no intervention"))

    for offset_type in ["week early",
                        "optimal time",
                        "week late"]:
        
        color_val = color_vals.get(offset_type, 0)
        handlelines.append(
            mlines.Line2D([], [],
                          color = leg_cmap(color_val),
                          label = offset_type))

    
    sc_handlelines.append(
        mlines.Line2D([], [],
                      ls = null_style,
                      dashes = (0.75, 0.75),
                      color = null_sc_color,
                      label = "sustained control\nonly"))

    ## plot dashed no control / sc only lines 
    for plotname in ["maintain-suppress-time-alone",
                     "fixed-alone",
                     "full-suppression-alone"]:
        plots[plotname].plot(null_time,
                             null_result[:, 1],
                             color = null_color,
                             ls = null_style,
                             dashes = null_dashes)

    for plotname in ["maintain-suppress-time-sc",
                     "fixed-sc",
                     "full-suppression-sc"]:
        plots[plotname].plot(null_sc_time,
                             null_sc_result[:, 1],
                             color = null_sc_color,
                             ls = null_style,
                             dashes = null_dashes,
                             label = "sustained control\nonly")

    for sir_model, suffix in zip([covid_sir, covid_sc_sir],
                                 ["-alone", "-sc"]):
        sir_model.b_func.tau = tau
        for strategy in ["maintain-suppress-time", "fixed", "full-suppression"]:
            sir_model.b_func.strategy = strategy
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
            for offset_val, offset_label in zip(offsets_timecourse,
                                                offset_labels):
                color_val = color_vals.get(offset_label, 0)
                cmap = cmaps[strategy]
                ## perform integration
                sir_model.b_func.t_i = t_i_opt + offset_val
                sir_model.reset()
                sir_model.integrate(t_sim_max)

                ## plot results
                plotname = strategy + suffix
                print("plotting to {}...\n".format(plotname))                
                plots[plotname].plot(
                    sir_model.time_ts,
                    sir_model.state_ts[:, 1],
                    color = cmap(color_val),
                    alpha = timecourse_alpha)
            pass
        pass
    
    plots["maintain-suppress-time-alone"].set_xlim(0, 300)
    plots["maintain-suppress-time-alone"].legend(
        handles = handlelines,
        frameon = True,
        fancybox = True,
        framealpha = 1,
        labelspacing = 0.25,
        handlelength = 0.75,
        loc = (0.55, 0.55))
    plots["maintain-suppress-time-sc"].legend(
        handles = sc_handlelines,
        frameon = True,
        fancybox = True,
        framealpha = 1,
        labelspacing = 0.25,
        handlelength = 0.75,
        loc = (0.55, 0.75))


    sweep_results = {}
    for strategy in ["maintain-suppress-time", "fixed", "full-suppression"]:
        plotname_vs_sc = strategy + "-vs-sc"
        alone_varname = strategy + "-sweep-alone"
        sc_varname = strategy + "-sweep-sc"
        cmap = cmaps[strategy]
        offsets, alone_sweep = calc_offset_sweep(
            tau,
            covid_sir,
            1000,
            strategy = strategy)

        offsets, sc_sweep = calc_offset_sweep(
            tau,
            covid_sc_sir,
            1000,
            strategy = strategy)

        sweep_results[sc_varname] = sc_sweep
        sweep_results[alone_varname] = alone_sweep
        plots[plotname_vs_sc].axhline(
            null_sc_peak,
            linestyle = "dashed",
            dashes = null_dashes,
            color = null_sc_color,
            lw = 3)
        plots[plotname_vs_sc].axhline(
            null_alone_peak,
            linestyle = "dashed",
            dashes = null_dashes,
            color = "black",
            lw = 3)
        
        plots[plotname_vs_sc].plot(offsets,
                            alone_sweep,
                            color = cmap(offset_darkness))
        plots[plotname_vs_sc].plot(
            offsets,
            sc_sweep,
            linestyle = "dotted",
            color = cmap(offset_darkness))
        plots[plotname_vs_sc].set_xticks(np.arange(
            -params.offsets_for_sweep,
            params.offsets_for_sweep + 7,
            7))
        plots[plotname_vs_sc].yaxis.set_label_position("right")
        plots[plotname_vs_sc].yaxis.tick_right()
        plots[plotname_vs_sc].set_xlim(-params.offsets_for_sweep,
                                       params.offsets_for_sweep)


    plots["maintain-suppress-time-alone"].set_xlim(0, 300)
    plots["maintain-suppress-time-alone"].set_ylim(ps.ymin_timecourse,
                                    ps.ymax_timecourse)
    plots["maintain-suppress-time-alone"].set_yticks(np.arange(0, 0.4, 0.1))
    plots["maintain-suppress-time-alone"].set_xticks(np.arange(0, 280, 90))

    for plotname, plot in plots.items():
        if "-vs-sc" not in plotname:
            plot.label_outer()
        elif not plot.is_last_row():
            plot.tick_params(axis = "x",
                             labelbottom=False)

    plots["fixed-alone"].set_ylabel("prevalence $I(t)$")

    plots["full-suppression-alone"].set_xlabel("time (days)")
    plots["full-suppression-sc"].set_xlabel("time (days)")
    plots["full-suppression-vs-sc"].set_xlabel("offset from $t_i^{\mathrm{opt}}$ (days)")
    plots["fixed-vs-sc"].set_ylabel("peak prevalence $I^{\max}$",
                                    rotation = 270,
                                    labelpad = 30)

    plots["maintain-suppress-time-alone"].set_title("Time-limited")
    plots["maintain-suppress-time-sc"].set_title("Sustained")
    plots["maintain-suppress-time-vs-sc"].set_title("Costs of mistiming")

    fig.tight_layout()
    fig.align_ylabels()
    fig.savefig(outpath)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("USAGE ./{} <outpath>\n".format(sys.argv[0]))
    else:
        main(sys.argv[1])
