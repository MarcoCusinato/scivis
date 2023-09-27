import numpy as np
import matplotlib.pyplot as plt
from scidata.quantities.quantities import SimulationAnalysis
from scidata.units.units import units
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

## Get the data

sim = SimulationAnalysis(args.sim_name, None, args.sim_path)

mean_ene = sim.neutrino_mean_energies_grey()
neutrino = sim.integrated_neutrino_lum()
mean_ene[:,0] = u.convert_to_ms(mean_ene[:,0])

## Design the plot

fig, ax = plt.subplots(1, 2, sharex=True, figsize=(25, 10))


lxlim, rxlim, rindex = find_xlim(mean_ene[:, 0], args.rxlim)
## Find xlimits
if rxlim >= 1000:
    mean_ene[:, 0] = u.convert_to_s(mean_ene[:, 0])
    ax[0].set_xlabel(r't-t$_b$ [s]')
    ax[1].set_xlabel(r't-t$_b$ [s]')
    ax[0].set_xlim(u.convert_to_s(lxlim), u.convert_to_s(rxlim))
else:
    ax[0].set_xlabel(r't-t$_b$ [ms]')
    ax[1].set_xlabel(r't-t$_b$ [ms]')
    ax[0].set_xlim(lxlim, rxlim)

#Set ylabels
ax[0].set_ylabel(r'L$_{\nu}$ [$10^{53}$erg$\cdot$s$^{-1}$]')

ax[1].set_ylabel(r'E$_{\nu}$ [MeV]')

#Set ylims
ax[0].set_ylim(0,1)
ax[1].set_ylim(0,30)

## Plotting 
ax[0].plot(mean_ene[:,0], neutrino[:,1] * 1e-53, c='#ef476f', lw=3)
ax[0].plot(mean_ene[:,0], neutrino[:,2] * 1e-53, c='#ffd166', lw=3)
ax[0].plot(mean_ene[:,0], neutrino[:,3] * 1e-53, c='#26547c', lw=3)
ax[1].plot(mean_ene[:,0], mean_ene[:,1], label=r'$\nu_e$', c='#ef476f', lw=3)
ax[1].plot(mean_ene[:,0], mean_ene[:,2], label=r'$\overline{\nu}_e$', c='#ffd166', lw=3)
ax[1].plot(mean_ene[:,0], mean_ene[:,3], label=r'$\nu_x$', c='#26547c', lw=3)
ax[1].legend()


fig.savefig(os.path.join(args.save_path, args.name))


