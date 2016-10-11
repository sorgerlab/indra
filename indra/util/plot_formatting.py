from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import matplotlib

fontsize=7

def set_fig_params():
    matplotlib.rcParams['font.sans-serif'] = 'Arial'
    matplotlib.rcParams['text.usetex'] = True
    matplotlib.rcParams['text.latex.preamble'] = [
            '\\usepackage{helvet}',
            '\\usepackage{sansmath}',
            '\\sansmath',
            '\\usepackage{underscore}',]

def format_axis(ax, label_padding=2, tick_padding=0, yticks_position='left'):
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position(yticks_position)
    ax.yaxis.set_tick_params(which='both', direction='out', labelsize=fontsize,
                             pad=tick_padding, length=2, width=0.5)
    ax.xaxis.set_tick_params(which='both', direction='out', labelsize=fontsize,
                             pad=tick_padding, length=2, width=0.5)
    ax.xaxis.labelpad = label_padding
    ax.yaxis.labelpad = label_padding
    ax.xaxis.label.set_size(fontsize)
    ax.yaxis.label.set_size(fontsize)

