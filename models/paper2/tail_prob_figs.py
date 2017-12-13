import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt
from indra.util import plot_formatting as pf

plt.ion()
pf.set_fig_params()


def plot_tail_prob(loc, scale=0.2):
    x = np.linspace(-2, 2, 200)
    plt.figure(figsize=(3, 2), dpi=150)
    yvals = norm.pdf(x, loc=loc, scale=scale)
    plt.plot(x, yvals, linestyle='none')
    ax = plt.gca()
    ax.set_yticks([])
    #ax.spines['left'].set_position(('axes', 0))
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    #ax.spines['top'].set_position(('axes', 0))
    #ax.spines['left'].set_smart_bounds(True)
    #ax.spines['bottom'].set_smart_bounds(True)
    ymax = np.max(yvals)*1.05
    ax.set_ylim(0, ymax)
    ax.fill_between(x, 0, yvals, where=x>0, facecolor='red')
    ax.fill_between(x, 0, yvals, where=x<0, facecolor='blue')
    pf.format_axis(ax)
    plt.xlabel('Log(Fold-Change)')
    plt.ylabel('Probability')
    plt.vlines(0, 0, ymax, color='k', linewidth=1)
    plt.vlines(loc, 0, np.max(yvals), color='gray', linewidth=1)
    plt.subplots_adjust(bottom=0.15)
    fontdict = {'size': 6}
    p_inc = norm.sf(0, loc=loc, scale=scale)
    p_dec = norm.cdf(0, loc=loc, scale=scale)
    plt.text(1, ymax, 'P(increase$|$D) = %.2f' % p_inc, color='red',
             fontdict=fontdict, horizontalalignment='center')
    plt.text(-1, ymax, 'P(decrease$|$D) = %.2f' % p_dec, color='blue',
             fontdict=fontdict, horizontalalignment='center')
    plt.savefig('tail_prob_loc_%.2f.pdf' % loc)

if __name__ == '__main__':
    plt.close('all')
    plot_tail_prob(0)
    plot_tail_prob(0.15)
    plot_tail_prob(1)
    plot_tail_prob(-1)

