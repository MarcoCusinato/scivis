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

shock_radius = sim.get_shock_radius(min_max_average=True)
shock_radius[:,0] = u.convert_to_ms(shock_radius[:, 0])
shock_radius[:,3] = u.convert_to_km(shock_radius[:, 3])

## Design the plot

fig, ax = plt.subplots(1, 1, sharex=True, figsize=(12.5, 10))

lxlim, rxlim, rindex = find_xlim(shock_radius[:, 0], args.rxlim)
## Find xlimits
if rxlim >= 1000:
    shock_radius[:, 0] = u.convert_to_s(shock_radius[:, 0])
    ax.set_xlabel(r't-t$_b$ [s]')
    ax.set_xlim(u.convert_to_s(lxlim), u.convert_to_s(rxlim))
else:
    ax.set_xlabel(r't-t$_b$ [ms]')
    ax.set_xlim(lxlim, rxlim)

ax.set_ylabel('R$_{shock}$ [km]')
ax.set_yscale('log')

##Plot quantities
ax.plot(shock_radius[:,0], shock_radius[:, 3], lw=3)

fig.savefig(os.path.join(args.save_path, args.name))


