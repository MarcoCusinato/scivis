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

gain_radius = sim.get_gain_radius(min_max_average=True)
gain_mass, gain_heating, \
    mass_accr = sim.get_masses_and_energies(gain_mass=True,
                                            gain_heating=True,
                                            mass_accretion_500km=True)
gain_radius[:,0] = u.convert_to_ms(gain_radius[:,0])
gain_radius[:,3] = u.convert_to_km(gain_radius[:,3])

## Design the plot

fig, ax = plt.subplots(2, 2, sharex=True, figsize=(25, 20))
fig.subplots_adjust(hspace=0)

lxlim, rxlim, rindex = find_xlim(gain_radius[:, 0], args.rxlim)
## Find xlimits
if rxlim >= 1000:
    gain_radius[:, 0] = u.convert_to_s(gain_radius[:, 0])
    ax[1,0].set_xlabel(r't-t$_b$ [s]')
    ax[1,1].set_xlabel(r't-t$_b$ [s]')
    ax[1,0].set_xlim(u.convert_to_s(lxlim), u.convert_to_s(rxlim))
else:
    ax[1,0].set_xlabel(r't-t$_b$ [ms]')
    ax[1,1].set_xlabel(r't-t$_b$ [ms]')
    ax[1,0].set_xlim(lxlim, rxlim)


##Set ylabels
ax[0,0].set_ylabel('M$_{gain}$ [M$_\odot$]')
ax[0,1].set_ylabel('M$_{accretion,\,500}$ [$\dot{M}_\odot$]')
ax[1,0].set_ylabel('R$_{gain}$ [km]')
ax[1,1].set_ylabel(r'Q$_{\nu}$/M$_{gain}$ [$10^{20}$ erg/g]')
##Set ylimits
ax[0,0].set_ylim(1e-3, 1)
ax[0,0].set_yscale('log')
ax[1,0].set_ylim(0, 130)
ax[0,1].set_ylim(0, 8)
ax[1,1].set_ylim(-1, 8)


##Plot quantities
ax[0,0].plot(gain_radius[:,0], gain_mass, lw=3)
ax[1,0].plot(gain_radius[:,0], gain_radius[:, 3], lw=3)
ax[0,1].plot(gain_radius[:,0], mass_accr, lw=3)
ax[1,1].plot(gain_radius[:,0], gain_heating / u.convert_to_grams(gain_mass) * 1e-20, lw=3)

fig.savefig(os.path.join(args.save_path, args.name))


