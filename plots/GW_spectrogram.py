import numpy as np
import matplotlib.pyplot as plt
from scidata.quantities.quantities import SimulationAnalysis
from scidata.units.units import units
import matplotlib.colors as colors
import os
import argparse
u = units()

## Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("--sim-name", required=True, help="Name of the simulation to analyse")
parser.add_argument("--sim-path", required=False, default=None, help="Path of the selected simulation")
parser.add_argument("--rxlim", required=False, default=None, type=float, help="Right x-axis limit in ms, if empty full GW is plotted")
parser.add_argument("--name", required=True, help="Plot name")
parser.add_argument("--save-path", required=True, help="Plot path")


args = parser.parse_args()

def find_xlim(time, right_lim):
    lxlim = -10
    if right_lim is None:
        rxlim = time[-1]
        rindex = time.shape[0]
    else:
        rxlim = right_lim
        rindex = np.argmax(time >= rxlim)
    return lxlim, rxlim, rindex

def set_ylimits(strain, ax):
    limits = np.max(np.abs(strain)) * 1.05
    if limits >= 100:
        limits = 100
    ax.set_ylim(-limits, limits)

def plot_spectrogram(fig, ax, ax_cb, time, freq, Zxx):
    pcm = ax.pcolormesh(time, freq, Zxx, shading='gouraud', cmap='inferno',
                             norm=colors.LogNorm(vmin=1e-24, vmax=1e-19))
    cbar = fig.colorbar(pcm, cax=ax_cb)
    cbar.set_label(r'$\frac{\mathrm{dE_{GW}}}{\mathrm{df}}$ [B$\cdot$HZ$^{-1}$]')

## Get the data
sim = SimulationAnalysis(args.sim_name, None, args.sim_path)
time, freq, Zxx = sim.GW_spectrogram(GW_norm=1/3.086e+22)
time = u.convert_to_ms(time)
freq /= 1000
GWs = sim.GW_Amplitudes()
GWs[:, 0] = u.convert_to_ms(GWs[:,0])
lxlim, rxlim, rindex = find_xlim(GWs[:,0], args.rxlim)

if sim.dim == 1:
    pass
elif sim.dim == 2:
    fig = plt.figure(figsize=(12.5, 10))
    ax = fig.subplot_mosaic(
        """A..
        B.b
        """,
        height_ratios = [1, 2],
        #withd ratio between columns
        width_ratios = [0.92, 0.03, 0.05],
        gridspec_kw={'hspace':0,
                    'wspace': 0}
    )

    ax["A"].sharex(ax["B"])
    
    
    if rxlim >= 1000:
        time = u.convert_to_s(time)
        GWs[:, 0] = u.convert_to_s(GWs[:,0])
        ax["B"].set_xlabel(r't-t$_b$ [s]')
        ax["B"].set_xlim(u.convert_to_s(lxlim), u.convert_to_s(rxlim))
    else:
        ax["B"].set_xlabel(r't-t$_b$ [ms]')
        ax["B"].set_xlim(lxlim, rxlim)
    
    set_ylimits(GWs[:rindex, 1], ax["A"])

    ax["A"].set_ylabel(r'$\mathcal{D}h_+$ [cm]')
    ax["B"].set_ylabel(r'$f$ [$10^3$Hz]')
    
    ax["B"].set_ylim(0,2)
    ax["A"].plot(GWs[:,0], GWs[:,1])
    plot_spectrogram(fig, ax['B'], ax['b'], time, freq, Zxx)
    

elif sim.dim == 3:
    fig = plt.figure(figsize=(25, 20))
    ax = fig.subplot_mosaic(
        """A...B..
           C.c.D.d
           .......
           E...F..
           G.g.H.h
        """,
        height_ratios = [1, 2, 0.3, 1, 2],
        #withd ratio between columns
        width_ratios = [0.92, 0.03, 0.04, 0.2, 0.92, 0.03, 0.04],
        gridspec_kw={'hspace':0,
                    'wspace': 0}
    )
    #share X axis
    ax["B"].sharex(ax["A"])
    ax["C"].sharex(ax["A"])
    ax["D"].sharex(ax["A"])
    ax["E"].sharex(ax["A"])
    ax["F"].sharex(ax["A"])
    ax["G"].sharex(ax["A"])
    ax["H"].sharex(ax["A"])
    
    #x limits and labels
    if rxlim >= 1000:
        time = u.convert_to_s(time)
        GWs[:, 0] = u.convert_to_s(GWs[:,0])
        ax["C"].set_xlabel(r't-t$_b$ [s]')
        ax["D"].set_xlabel(r't-t$_b$ [s]')
        ax["G"].set_xlabel(r't-t$_b$ [s]')
        ax["H"].set_xlabel(r't-t$_b$ [s]')
        ax["A"].set_xlim(u.convert_to_s(lxlim), u.convert_to_s(rxlim))
    else:
        ax["C"].set_xlabel(r't-t$_b$ [ms]')
        ax["D"].set_xlabel(r't-t$_b$ [ms]')
        ax["G"].set_xlabel(r't-t$_b$ [ms]')
        ax["H"].set_xlabel(r't-t$_b$ [ms]')
        ax["A"].set_xlim(lxlim, rxlim)

    #y labels GWs
    ax["A"].set_ylabel(r'$\mathcal{D}h_{\mathrm{eq, +}}$ [cm]')
    ax["B"].set_ylabel(r'$\mathcal{D}h_{\mathrm{pol, +}}$ [cm]')
    ax["E"].set_ylabel(r'$\mathcal{D}h_{\mathrm{eq, x}}$ [cm]')
    ax["F"].set_ylabel(r'$\mathcal{D}h_{\mathrm{pol, x}}$ [cm]')

    #y label spectro
    ax["C"].set_ylabel(r'$f$ [$10^3$Hz]')
    ax["D"].set_ylabel(r'$f$ [$10^3$Hz]')
    ax["G"].set_ylabel(r'$f$ [$10^3$Hz]')
    ax["H"].set_ylabel(r'$f$ [$10^3$Hz]')

    # y limits GWs
    set_ylimits(GWs[:rindex, 1], ax["A"])    
    set_ylimits(GWs[:rindex, 2], ax["B"])    
    set_ylimits(GWs[:rindex, 3], ax["E"])    
    set_ylimits(GWs[:rindex, 4], ax["F"])

    # y limits spectrogram
    ax["C"].set_ylim(0,2)
    ax["D"].set_ylim(0,2)
    ax["G"].set_ylim(0,2)
    ax["H"].set_ylim(0,2)

    ## Plot GWs
    ax["A"].plot(GWs[:,0], GWs[:,1], lw=2)
    ax["B"].plot(GWs[:,0], GWs[:,2], lw=2)
    ax["E"].plot(GWs[:,0], GWs[:,3], lw=2)
    ax["F"].plot(GWs[:,0], GWs[:,4], lw=2)

    ## Plot spectra
    plot_spectrogram(fig, ax['C'], ax['c'], time[:,0], freq[:,0], Zxx[...,0])
    plot_spectrogram(fig, ax['D'], ax['d'], time[:,1], freq[:,1], Zxx[...,1])
    plot_spectrogram(fig, ax['G'], ax['g'], time[:,2], freq[:,2], Zxx[...,2])
    plot_spectrogram(fig, ax['H'], ax['h'], time[:,3], freq[:,3], Zxx[...,3])


fig.savefig(os.path.join(args.save_path, args.name))

