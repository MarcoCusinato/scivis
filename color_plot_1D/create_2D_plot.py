import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker

def prepare_plot(xlim, ylim, xlogscale, ylogscale, figsize):
    fig = plt.figure(figsize=figsize)
    ax = fig.subplot_mosaic(
                    """A.a""",
                    width_ratios = [0.92, 0.03, 0.05])
    if xlim is not None:
        ax['A'].set_xlim(xlim)
    if ylim is not None:
        ax['A'].set_ylim(ylim)
    if xlogscale:
        ax['A'].set_xscale('log')
    if ylogscale:
        ax['A'].set_yscale('log')
    return fig, ax

def set_labels(xlabel, ylabel, cb_label, ax, cbar):
    ax['A'].set_xlabel(xlabel)
    ax['A'].set_ylabel(ylabel)
    cbar.set_label(cb_label)

def levels(vmin, vmax, nlevels, logscale):
    if logscale:
        levels = 10 ** np.linspace(vmin, vmax, nlevels)
        loc = ticker.LogLocator
    else:
        levels = np.linspace(vmin, vmax, nlevels)
        loc = ticker.LinearLocator
    return levels, loc

def plot_quantities(x, y, quantity, fig, ax, vmin, vmax, nlevels, logscale, cmap):
    X, Y = np.meshgrid(x, y)
    lev, locator = levels(vmin, vmax, nlevels, logscale)
    pcm = ax['A'].contourf(X, Y, quantity, levels=lev, locator=locator(), cmap=cmap)
    cbar = fig.colorbar(pcm, cax=ax["a"])
    try:
        pcm = ax['A'].contourf(X, Y, quantity, levels=lev, locator=locator(), cmap=cmap)
        cbar = fig.colorbar(pcm, cax=ax["a"])
    except:
        pcm = ax['A'].contourf(X, Y, quantity.T, levels=lev, locator=locator(), cmap=cmap)
        cbar = fig.colorbar(pcm, cax=ax["a"])
    return cbar, ax