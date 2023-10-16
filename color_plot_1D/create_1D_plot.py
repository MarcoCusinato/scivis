import matplotlib.pyplot as plt
from scidata.units.units import units
import numpy as np
from matplotlib import ticker

u = units()


def create_averaged_quantities(sim, quantity, rxlim):
    file_list = sim.file_list_hdf()
    indices, time = sim.get_PNS_radius(indices=True, ret_time=True)
    time = u.convert_to_ms(time)
    if rxlim is not None:
        t_index = np.argmax(time > rxlim) + 1
        indices = indices[:t_index]
        file_list = file_list[:t_index]
    radius = u.convert_to_km(sim.cell.radius(sim.ghost))
    quantity_av = np.zeros((radius.size, time.size))
    quantity_method = getattr(sim, quantity)
    for (file, i) in zip(file_list, indices):
        data = sim.open_h5(file)
        quantity_data_tmp = quantity_method(data)
        if sim.dim == 1:
            quantity_av[:, i] = quantity_data_tmp
        elif sim.dim == 2:
            quantity_av[:, i] = quantity_data_tmp.mean(axis = 0)
        elif sim.dim == 3:
            quantity_av[:, i] = quantity_data_tmp.mean(axis = (0,1))
        sim.close_h5(data)
    return time, quantity_av, radius

def find_xlim(time, right_lim):
    lxlim = -10
    if right_lim is None:
        rxlim = time[-1]
    else:
        rxlim = right_lim
    return lxlim, rxlim

def set_xlimits_labels(ax, time, right_lim):
    lxlim, rxlim = find_xlim(time, right_lim)
    if rxlim >= 1000:
        time = u.convert_to_s(time)
        ax.set_xlabel(r't-t$_b$ [s]')
        ax.set_xlim(u.convert_to_s(lxlim), u.convert_to_s(rxlim))
    else:
        ax.set_xlabel(r't-t$_b$ [ms]')
        ax.set_xlim(lxlim, rxlim)

def set_ylimits_labels(ax, radius, top_lim):
    bottom_lim = radius[0]
    if top_lim is None:
        top_lim = radius[-1]
    ax.set_ylim(bottom_lim, top_lim)
    ax.set_yscale('log')
    ax.set_ylabel('R [km]')

def linspace_levels(quant, physical_zero, nvmin, nvmax):
    vmin = np.floor(np.nanmin(quant))
    vmax = np.ceil(np.nanmax(quant))
    if vmin < 0 and physical_zero:
        quant[quant < 0] = 0
        vmin = np.floor(np.nanmin(quant[np.nonzero(quant)]))
    if nvmin is not None:
        vmin = nvmin
    if nvmax is not None:
        vmax = nvmax
    
    if vmax - vmin > 50:
        nlevels = int(vmax - vmin)
    else:
        nlevels = 50
    return np.linspace(vmin, vmax, nlevels), ticker.LinearLocator

def logspaced_levels(quant, physical_zero, nvmin, nvmax):
    if physical_zero:
        quant[quant < 0] = 0
    if np.nanmin(quant) < 0 and (nvmin is None):
        return linspace_levels(quant, nvmin, nvmax)
    if nvmin is None:
        vmin = np.floor(np.log10(np.nanmin(quant[np.nonzero(quant)])))
    else:
        vmin = np.floor(np.log10(nvmin))

    if nvmax is None:
        vmax = np.floor(np.log10(np.nanmax(quant))) + 1
    else:
        vmax = np.ceil(np.log10(nvmax))
    
    nlevels = int(vmax - vmin)
    levels = 10 ** np.linspace(vmin, vmax, nlevels*2 + 1)
    return levels, ticker.LogLocator


def set_up_plot(sim, quantity, rxlim=None, tylim=None, log=False, cbar_label = None, cmap='inferno',
                physical_zero=False, vmin=None, vmax=None):
    time, quant, radius = create_averaged_quantities(sim, quantity, rxlim)
    fig = plt.figure(figsize=(12.5, 10))
    ax = fig.subplot_mosaic(
                    """A.a""",
                    width_ratios = [0.92, 0.03, 0.05],
                    gridspec_kw={'wspace': 0})
    set_xlimits_labels(ax['A'], time, rxlim)
    set_ylimits_labels(ax['A'], radius, tylim)
    X, Y = np.meshgrid(time, radius)

    if log:
        levels, locator = logspaced_levels(quant, physical_zero, vmin, vmax)
    else:
        levels, locator = linspace_levels(quant, physical_zero, vmin, vmax)
  

    pcm = ax['A'].contourf(X, Y, quant, levels=levels, locator=locator(), cmap=cmap)
    cbar = fig.colorbar(pcm, cax=ax["a"])
    if cbar_label is None:
        cbar.set_label(quantity)
    else:
        cbar.set_label(cbar_label)
    return fig






