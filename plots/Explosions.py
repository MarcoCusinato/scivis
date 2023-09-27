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

time, ej_mass, expl_ene = sim.get_masses_and_energies(time=True, ejected_mass=True, explosion_energy=True)
time = u.convert_to_ms(time)

## Design the plot

fig, ax = plt.subplots(2, 1, sharex=True, figsize=(25, 10))
fig.subplots_adjust(hspace=0)

lxlim, rxlim, rindex = find_xlim(time, args.rxlim)
## Find xlimits
if rxlim >= 1000:
    time = u.convert_to_s(time)
    ax[1].set_xlabel(r't-t$_b$ [s]')
    ax[0].set_xlim(u.convert_to_s(lxlim), u.convert_to_s(rxlim))
else:
    ax[1].set_xlabel(r't-t$_b$ [ms]')
    ax[0].set_xlim(lxlim, rxlim)

ax[0].set_ylabel('E$_{expl}$ [erg]')
ax[1].set_ylabel('M$_{ej}$ [M$_\odot$]')

ax[0].set_yscale('log')
ax[1].set_yscale('log')
ax[0].set_ylim(1e47, 1e52)
ax[1].set_ylim(1e-3, 1e1)

##Plot quantities
ax[0].plot(time, expl_ene, lw=3)
ax[1].plot(time, ej_mass, lw=3)

fig.savefig(os.path.join(args.save_path, args.name))


