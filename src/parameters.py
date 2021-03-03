#!/usr/bin/env python3

####################################################
# filename: parameters.py
# author: Dylan Morris <dhmorris@princeton.edu>
#
# description: sets parameters for running
# various numerical analyses of the Morris/Rossine
# model of optimal and mistimed intervention
####################################################

import numpy as np

# default epidemiological parameters
R0_default = 3
gamma_default = 1/14

# default intervention parameters
tau_default = 28
taus_figure_interventions = [14, 28, 56]
taus_sweep = np.linspace(0.5, 90, 90)
sc_reduction = 0.75 # sustained control reduction in R0
offsets_default = [7, -7, 0] # mistiming by +/- a week
offsets_for_sweep = 14 # for sweep in the heatmaps

## parameter sweeps for heatmaps
heatmap_fineness = 100
heatmap_taus = np.linspace(1, 90, heatmap_fineness)
heatmap_gammas = np.linspace(0.01, 0.7, heatmap_fineness)
heatmap_R0s = np.linspace(1.1, 5, heatmap_fineness)
heatmap_offsets = np.linspace(-10,
                              10,
                              heatmap_fineness)

fixed_heatmap_fineness = 25
fixed_heatmap_taus = np.linspace(2, 90, fixed_heatmap_fineness)
fixed_heatmap_gammas = np.linspace(0.001, 0.7, fixed_heatmap_fineness)
fixed_heatmap_R0s = np.linspace(1.01, 5, fixed_heatmap_fineness)
fixed_heatmap_offsets = np.linspace(-offsets_for_sweep,
                                    offsets_for_sweep,
                                    fixed_heatmap_fineness)

# default initial values
I_init_default = 1e-6
R_init_default = 0
S_init_default = 1 - I_init_default - R_init_default
inits_default = np.array([S_init_default,
                          I_init_default,
                          R_init_default])

# sweep initial values
inits_sweep_small = np.array([1-1e-8, 1e-8, 0])
