import matplotlib.pyplot as plt
from scidata.quantities.quantities import SimulationAnalysis
from scidata.units.units import units
import argparse
import numpy as np
import os

u = units()

## Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("--sim-name", required=True, help="Name of the simulation to analyse")
parser.add_argument("--sim-path", required=False, default=None, help="Path of the selected simulation")
parser.add_argument("--rxlim", required=False, default=None, type=float, help="Right x-axis limit in ms, if empty full GW is plotted")
parser.add_argument("--name", required=True, help="Plot name")
parser.add_argument("--save-path", required=True, help="Plot path")

args = parser.parse_args()

## Select simulation

sim = SimulationAnalysis(args.sim_name, None, args.sim_path)
dim = sim.dim

## Get the GWs
GWs = sim.GW_Amplitudes()
GWs[:, 0] = u.convert_to_ms(GWs[:, 0])

def set_ylimits(strain, ax):
    limits = np.max(np.abs(strain)) * 1.05
    if limits >= 100:
        limits = 100
    ax.set_ylim(-limits, limits)

def find_xlim(strain, right_lim):
    lxlim = -5
    if right_lim is None:
        rxlim = strain[-1, 0]
        rindex = strain[:, 0].shape[0]
    else:
        rxlim = right_lim
        rindex = np.argmax(strain[:,0] >= rxlim)
    return lxlim, rxlim, rindex


if dim == 1:
    pass
elif dim == 2:
    ## Set the figure
    fig, ax = plt.subplots(1, 1, figsize=(12.5, 10))
    lxlim, rxlim, rindex = find_xlim(GWs, args.rxlim)

    ## Find xlimits
    if rxlim >= 1000:
        GWs[:, 0] = u.convert_to_s(GWs[:, 0])
        ax.set_xlabel(r't-t$_b$ [s]')
        ax.set_xlim(u.convert_to_s(lxlim), u.convert_to_s(rxlim))
    else:
        ax.set_xlabel(r't-t$_b$ [ms]')
        ax.set_xlim(lxlim, rxlim)

    ## Find ylimits
    ax.set_ylabel(r'$\mathcal{D}h_+$ [cm]')
    set_ylimits(GWs[:rindex, 1], ax)
    ax.plot(GWs[:,0], GWs[:,1])
    fig.savefig(os.path.join(args.save_path, args.name))
elif dim == 3:
    ## Setting up the plot
    fig, ax = plt.subplots(2, 2, figsize=(25, 20), sharex=True)
    fig.subplots_adjust(hspace=0)
    #Find xlimits
    lxlim, rxlim, rindex = find_xlim(GWs, args.rxlim)
    if rxlim >= 1000:
        GWs[:, 0] = u.convert_to_s(GWs[:, 0])
        ax[0, 0].set_xlim(u.convert_to_s(lxlim), u.convert_to_s(rxlim))
        measure = '[s]'
    else:
        ax[0,0].set_xlim(lxlim, rxlim)
        measure = '[ms]'
    ## Set the labels
    ax[1,0].set_xlabel('t-t$_b$ ' + measure)
    ax[1,1].set_xlabel('t-t$_b$ ' + measure)
    ax[0,0].set_ylabel(r'$\mathcal{D}h_{+, eq}$ [cm]')
    ax[1,0].set_ylabel(r'$\mathcal{D}h_{+, pol}$ [cm]')
    ax[0,1].set_ylabel(r'$\mathcal{D}h_{x, eq}$ [cm]')
    ax[1,1].set_ylabel(r'$\mathcal{D}h_{x, pol}$ [cm]')
    ## Setting ylimits
    set_ylimits(GWs[:rindex, 1], ax[0, 0])
    set_ylimits(GWs[:rindex, 2], ax[1, 0])
    set_ylimits(GWs[:rindex, 3], ax[0, 1])
    set_ylimits(GWs[:rindex, 4], ax[1, 1])
    #Plotting things
    ax[0,0].plot(GWs[:,0], GWs[:,1], lw=3)
    ax[1,0].plot(GWs[:,0], GWs[:,2], lw=3)
    ax[0,1].plot(GWs[:,0], GWs[:,3], lw=3)
    ax[1,1].plot(GWs[:,0], GWs[:,4], lw=3)
    fig.savefig(os.path.join(args.save_path, args.name))
    



