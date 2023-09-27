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
nue, nua, nux = sim.get_neutrino_sphere_radius(min_max_average_nue=True, min_max_average_nua=True, min_max_average_nux=True)
nue[:,0] = u.convert_to_ms(nue[:,0])
nue[:,3] = u.convert_to_km(nue[:,3])
nua[:,3] = u.convert_to_km(nua[:,3])
nux[:,3] = u.convert_to_km(nux[:,3])

## Design the plot

fig, ax = plt.subplots(1, 1, sharex=True, figsize=(12.5, 10))

lxlim, rxlim, rindex = find_xlim(nue[:,0], args.rxlim)
## Find xlimits
if rxlim >= 1000:
    nue[:, 0] = u.convert_to_s(nue[:, 0])
    ax.set_xlabel(r't-t$_b$ [s]')
    ax.set_xlim(u.convert_to_s(lxlim), u.convert_to_s(rxlim))
else:
    ax.set_xlabel(r't-t$_b$ [ms]')
    ax.set_xlim(lxlim, rxlim)

#Set ylabels
ax.set_ylabel(r'R$_{\nu}$ [km]')


#Set ylims
ax.set_ylim(0,120)


## Plotting 
ax.plot(nue[:,0], nue[:,3], c='#ef476f', lw=3, label = r'$\nu_\mathrm{e}$')
ax.plot(nue[:,0], nua[:,3], c='#ffd166', lw=3, label = r'$\overline{\nu}_\mathrm{e}$')
ax.plot(nue[:,0], nux[:,3], c='#26547c', lw=3, label = r'$\nu_\mathrm{x}$')

ax.legend()


fig.savefig(os.path.join(args.save_path, args.name))


