#!/usr/bin/env python3

####################################################
# filename: plotting_style.py
# author: Dylan Morris <dhmorris@princeton.edu>
#
# description: sets plot styling for the
# Morris/Rossine model of
# mistimed interventions 
####################################################

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('seaborn-whitegrid')
mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = "x-large"
mpl.rcParams['xtick.labelsize'] = "x-large"
mpl.rcParams['ytick.labelsize'] = "x-large"
mpl.rcParams['axes.titlesize'] = "large"
mpl.rcParams['text.usetex'] = True
mpl.rcParams['axes.formatter.use_mathtext'] = True
mpl.rcParams['axes.formatter.limits'] = ((-3, 3))
mpl.rcParams['axes.grid'] = True
mpl.rcParams['legend.frameon'] = True
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.framealpha'] = 1
mpl.rcParams['legend.title_fontsize'] = "small"
mpl.rcParams['legend.fontsize'] = "small"
mpl.rcParams['lines.linewidth'] = 5
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.left'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.bottom'] = False
## custom colors
ourblue = "#0571a3"
ourred = "#c20a1f"

## colormaps for different model conditions
opt_cmap = plt.cm.Greens
fixed_cmap = plt.cm.Blues
full_suppression_cmap = plt.cm.Reds
cmap_dict = {
    "optimal": opt_cmap,
    "maintain-suppress-time": opt_cmap,
    "fixed": fixed_cmap,
    "full-suppression": full_suppression_cmap}

## multipanel parameters
default_letter_loc = (0, 1.125)
default_letter_size = "large"

# default limits
ymax_timecourse = 0.33
ymin_timecourse = -0.03
## useful plotting functions
def add_bounding_subplot(figure, position=None):
    if position is None:
        position = 111
    if ((type(position) is tuple) or
        (type(position) is list)):
        ax = figure.add_subplot(
            position[0], position[1], position[2])
    else:
        ax = figure.add_subplot(position)
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.grid(b=False)
    ax.patch.set_alpha(0)
    ax.tick_params(labelcolor='none',
                   grid_alpha=0,
                   top = False,
                   bottom = False,
                   left = False,
                   right = False)
    for label in ax.get_xticklabels():
        label = ''
    for label in ax.get_yticklabels():
        label = ''
    ax.set_zorder(0)
    return ax


def get_position(position_dict,
                 n_cols):
    pos = ((position_dict["row"] - 1) *
           n_cols + position_dict["col"])
    return pos

def setup_multipanel(fig, plot_positions,
                     gridspec=None,
                     add_letters=True,
                     letter_loc=None,
                     verbose=False,
                     letters=None,
                     upper=False):

    if gridspec is None:
        n_rows = np.max([plot_dict["row"] for plot_dict in plot_positions])
        n_cols = np.max([plot_dict["col"] for plot_dict in plot_positions])
        
    plots = {}

    if letters is None:
        letters = ["a", "b", "c", "d", "e", "f",
                   "g", "h", "i", "j", "k", "l",
                   "m", "n", "o", "p", "q", "r",
                   "s", "t", "u", "v", "w", "x",
                   "y", "z", "aa", "bb", "cc"]

    if upper:
        letters = [let.upper() for let in letters]

    for plot_no, plot_dict in enumerate(plot_positions):

        sharex = plot_dict["sharex"]
        sharey = plot_dict["sharey"]

        if sharex is not None:
            plot_sharex = plots[sharex]
        else:
            plot_sharex = sharex
            
        if sharey is not None:
            plot_sharey = plots[sharey]
        else:
            plot_sharey = sharey

        if gridspec is None:
            position_id = get_position(plot_dict,
                                       n_cols)
            if verbose:
                print("Adding plot {} at position {}"
                      "".format(plot_dict['name'],
                                    position_id))
            ax = fig.add_subplot(n_rows,
                                 n_cols,
                                 position_id,
                                 sharex=plot_sharex,
                                 sharey=plot_sharey)
        else:
            grid_pos = plot_dict["grid_position"]
            
            if verbose:
                print("Adding plot {} at position {}"
                      "".format(plot_dict['name'],
                                    grid_pos))

            ax = fig.add_subplot(gridspec[grid_pos],
                                 sharex=plot_sharex,
                                 sharey=plot_sharey)
        if add_letters:
            letter_locations = plot_dict.get("letter_loc", letter_loc)
            if letter_locations is None:
                letter_x, letter_y = default_letter_loc
                # handle spanny letters equivalently
                # by default
                if gridspec is not None:
                    (_, _,
                     row_start,
                     row_stop,
                     col_start,
                     col_stop) = ax.get_subplotspec().get_rows_columns()
                    colspan = 1 + col_stop - col_start
                    rowspan = 1 + row_start - row_stop
                    x_anchor = (
                        0 if abs(1 - letter_x) > abs(0 - letter_x) else 1)
                    y_anchor = (
                        0 if abs(1 - letter_y) > abs(0 - letter_y) else 1)
                    l_x_rel = letter_x - x_anchor
                    l_y_rel = letter_y - y_anchor
                    letter_x = x_anchor + l_x_rel / colspan
                    letter_y = y_anchor + l_y_rel / rowspan                  
            elif type(letter_locations) is tuple:
                letter_x, letter_y = letter_locations
            if plot_dict.get("include_letter", True):
                ax.text(letter_x,
                        letter_y,
                        letters[plot_no],
                        transform = ax.transAxes,
                        fontsize = default_letter_size,
                        fontweight = "bold",
                        va = "top")
        
        plots[plot_dict["name"]] = ax

    return plots
