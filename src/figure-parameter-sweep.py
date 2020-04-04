#!/usr/bin/env python3

####################################################
# filename: figure-parameter-sweep.py
# author: Dylan Morris <dhmorris@princeton.edu>
#
# description: heatmaps showing how our results
# vary with parameters
####################################################
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
from matplotlib.gridspec import GridSpec

import parameters as params
import plotting_style as ps
import InterventionSIR as sir
import optimize_interventions as oi


def setup_figure():
    width = 16
    height = 15

    figure = plt.figure(figsize = (width, height))

    offset_row, early_late_row= 0, 1
    tau_col, gamma_col, cbar_col = 0, 1, 2
    gs = GridSpec(2, 3, width_ratios=[1, 1, 0.05], height_ratios=[1, 1])
    
    plot_positions = [
        {"name": "offset-vs-tau",
         "grid_position": np.s_[offset_row, tau_col],
         "sharex": None,
         "sharey": None},
        
        {"name": "offset-vs-gamma",
         "grid_position": np.s_[offset_row, gamma_col],
         "sharex": None,
         "sharey": "offset-vs-tau"},
        
        {"name": "early-late-tau",
         "grid_position": np.s_[early_late_row,
                                tau_col],
         "sharex": "offset-vs-tau",
         "sharey": None},
        
        {"name": "early-late-gamma",
         "grid_position": np.s_[early_late_row,
                                gamma_col],
        "sharex": "offset-vs-gamma",
         "sharey": "early-late-tau"},
        
        {"name": "offset-cbar",
         "grid_position": np.s_[offset_row,
                                cbar_col],
         "include_letter": False,
         "sharex": None,
         "sharey": None},
        
        {"name": "early-late-cbar",
         "grid_position": np.s_[early_late_row,
                                cbar_col],
         "include_letter": False,
         "sharex": None,
         "sharey": None}
    ]

    plots = ps.setup_multipanel(figure,
                                plot_positions,
                                gridspec = gs,
                                letter_loc = (-0.025, 1.1))
    
    return(figure, plots)



def main(tau_offset_path,
         gamma_offset_path,
         tau_early_late_path,
         gamma_early_late_path,
         outpath):
    fig, plots = setup_figure()

    offset_cmap = plt.cm.YlGnBu
    early_late_cmap = plt.cm.magma_r
    
    tau_offset_sweep = pd.read_csv(
        tau_offset_path, header = 0,
        index_col = 0)
    gamma_offset_sweep = pd.read_csv(
        gamma_offset_path,
        header = 0,
        index_col = 0)
    tau_early_late_sweep = pd.read_csv(
        tau_early_late_path, header = 0,
        index_col = 0)
    gamma_early_late_sweep = pd.read_csv(
        gamma_early_late_path,
        header = 0,
        index_col = 0)

    null_peak = oi.I_max(1, 0, params.R0_default)

    min_tau_o = float(tau_offset_sweep.columns.min())
    max_tau_o = float(tau_offset_sweep.columns.max())
    min_offset_t_o = tau_offset_sweep.index.min()
    max_offset_t_o = tau_offset_sweep.index.max()
    
    extent_tau_offset = [
        min_tau_o,
        max_tau_o,
        max_offset_t_o,
        min_offset_t_o]

    min_gamma_o = float(gamma_offset_sweep.columns.min())
    max_gamma_o = float(gamma_offset_sweep.columns.max())
    min_offset_g_o = gamma_offset_sweep.index.min()
    max_offset_g_o = gamma_offset_sweep.index.max()
    
    extent_gamma_offset = [
        min_gamma_o,
        max_gamma_o,
        max_offset_g_o,
        min_offset_g_o]

    min_tau_el = float(tau_early_late_sweep.columns.min())
    max_tau_el = float(tau_early_late_sweep.columns.max())
    min_R0_t_el = tau_early_late_sweep.index.min()
    max_R0_t_el = tau_early_late_sweep.index.max()
    
    extent_tau_early_late = [
        min_tau_el,
        max_tau_el,
        max_R0_t_el,
        min_R0_t_el]

    min_gamma_el = float(gamma_early_late_sweep.columns.min())
    max_gamma_el = float(gamma_early_late_sweep.columns.max())
    min_R0_g_el = gamma_early_late_sweep.index.min()
    max_R0_g_el = gamma_early_late_sweep.index.max()
    
    extent_gamma_early_late = [
        min_gamma_el,
        max_gamma_el,
        max_R0_g_el,
        min_R0_g_el]
  
    plots["offset-vs-tau"].imshow(
        tau_offset_sweep,
        cmap = offset_cmap,
        extent = extent_tau_offset,
        aspect = "auto",
        vmin = 0,
        vmax = null_peak)
    plots["offset-vs-tau"].invert_yaxis()

    os_gamma = plots["offset-vs-gamma"].imshow(
        gamma_offset_sweep,
        cmap = offset_cmap,
        extent = extent_gamma_offset,
        aspect = "auto",
        vmin = 0,
        vmax = null_peak)
    plots["offset-vs-gamma"].invert_yaxis()

    max_diff = max(
        np.max(np.array(tau_early_late_sweep)),
        np.max(np.array(gamma_early_late_sweep)))    
    print(max_diff)
    min_diff = 0
    
    plots["early-late-tau"].imshow(
        tau_early_late_sweep,
        cmap = early_late_cmap,
        extent = extent_tau_early_late,
        vmin = min_diff,
        vmax = max_diff,
        aspect = "auto")
    plots["early-late-tau"].invert_yaxis()

    el_gamma = plots["early-late-gamma"].imshow(
        gamma_early_late_sweep,
        cmap = early_late_cmap,
        extent = extent_gamma_early_late,
        vmin = min_diff,
        vmax = max_diff,
        aspect = "auto")
    plots["early-late-gamma"].invert_yaxis()
    
    
    plots["offset-vs-tau"].set_ylabel("offset from $t_i^{\mathrm{opt}}$"
                                      " (days)")

    plots["early-late-tau"].set_ylabel("basic reproduction number "
                                       "$\mathcal{R}_0$")

    plots["early-late-gamma"].set_xticks(
        np.arange(0, max_gamma_el + 0.1, 0.3))
    plots["offset-vs-gamma"].set_xticks(
        np.arange(0, max_gamma_el + 0.1, 0.1))
    plots["early-late-gamma"].xaxis.set_major_formatter(
        FormatStrFormatter('$%g$'))
    plots["early-late-gamma"].set_xlabel(
        "recovery rate $\gamma$")
    
    plots["offset-vs-tau"].set_yticks(np.arange(-10, 11, 5))
    plots["offset-vs-gamma"].set_yticks(np.arange(-10, 11, 5))
    plots["early-late-tau"].set_yticks(np.arange(1, 6, 1))
    plots["early-late-gamma"].set_yticks(np.arange(1, 6, 1))

    plots["offset-vs-tau"].set_xticks(
        np.arange(0, max_tau_el + 10, 25))
    plots["early-late-tau"].set_xticks(
        np.arange(0, max_tau_el + 10, 30))
    
    plots["early-late-tau"].set_xlabel(
        "duration $\\tau$")
    cb_os = fig.colorbar(
        os_gamma,
        cax = plots["offset-cbar"])
    cb_el = fig.colorbar(
        el_gamma,
        cax = plots["early-late-cbar"])
    cb_os.outline.set_visible(False)
    cb_el.outline.set_visible(False)
    plots["offset-cbar"].set_ylabel("peak $I^{\max}$",
                                    rotation = 270,
                                    labelpad = 40)
    plots["early-late-cbar"].set_ylabel("Asymmetry ($I^{\max}$ "
                                        "late $- I^{\max}$ early)",
                                        rotation = 270,
                                        labelpad = 40)
    for plotname, plot in plots.items():
        if "cbar" not in plotname:
            plot.label_outer()
        plot.grid(b = False)
    
    fig.tight_layout()
    fig.savefig(outpath)


if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("USAGE ./{} "
              "<path to offset tau sweep> "
              "<path to offset gamma sweep> "
              "<path to early/late tau sweep> "
              "<path to early/late gamma sweep> "
              "<outpath>\n".format(sys.argv[0]))
    else:
        main(sys.argv[1],
             sys.argv[2],
             sys.argv[3],
             sys.argv[4],
             sys.argv[5])
