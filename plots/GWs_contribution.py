import matplotlib.pyplot as plt
from scidata.quantities.quantities import SimulationAnalysis
import numpy as np
from scidata.units.units import units
import argparse
import os

u = units()

parser = argparse.ArgumentParser()

parser.add_argument("--sim-name", required=True, help="Name of the simulation to analyse")
parser.add_argument("--sim-path", required=False, default=None, help="Path of the selected simulation")
parser.add_argument("--rxlim", required=False, default=None, type=float, help="Right x-axis limit in ms, if empty full GW is plotted")
parser.add_argument("--name", required=True, help="Plot name")
parser.add_argument("--save-path", required=True, help="Plot path")

args = parser.parse_args()
## GET THE SIMULATION
sim = SimulationAnalysis(args.sim_name, None, args.sim_path)
## GET THE GRID 
radius, time, AE220 = sim.AE220()
radius = u.convert_to_km(radius)
time = u.convert_to_ms(time)
R, T = np.meshgrid(radius, time)

## SUM UP THE CONTRIBUTIONS
rad_20 = np.argmax(radius>=20) + 1
GWs_contr = np.zeros((len(time), 3))
GWs_20km = np.sum(AE220[:rad_20, :], axis=0)
GWs_other = np.sum(AE220[rad_20:, :], axis=0)

## FIND THE GW STRAIN
GWs = sim.GW_Amplitudes()
GWs[:, 0] = u.convert_to_ms(GWs[:,0])

## STRAIN LIMITS
if args.rxlim is not None:
    GWs_lim = np.max( np.abs( GWs[:np.argmax( GWs[:, 0] >= args.rxlim), 1] ) )
else:
    GWs_lim = np.max( np.abs( GWs[:,1] ))
if GWs_lim > 100:
    GWs_lim = 100

## CREATE THE PLOT
fig = plt.figure(figsize=(5.7, 7.3), layout='tight')
axd = fig.subplot_mosaic(
    """A..
       B.b""",
    #height ratio between rows
    height_ratios = [1, 3],
    #withd ratio between columns
    width_ratios = [1, 0.05, 0.07],
    gridspec_kw={
        "wspace": 0,
        "hspace": 0,
    }
)

## FIX THE LABELS AND LIMITS
axd["A"].tick_params(top=True, labeltop=True,
                        bottom=False, labelbottom=False,
                        left=True, labelleft=True,
                        right=False, labelright=False)
## GWs PLOT
axd["A"].xaxis.set_label_position('top')
axd["A"].plot(time, GWs_other, label='$\leq$ 20 km')
axd["A"].plot(time, GWs_20km, label = '$\geq$ 20 km')
axd["A"].set_xlim(-5, args.rxlim)
axd["A"].set_ylim(-GWs_lim, GWs_lim)
axd["A"].set_xlabel('t-t$_b$ [ms]')
axd["A"].set_ylabel('$\mathcal{D}h_+$')

## 2D PLOT
axd["B"].set_xscale('log')
axd["B"].set_ylim(-5, args.rxlim)
axd["B"].set_ylabel('t-t$_b$ [ms]')
axd["B"].set_xlabel('radius [km]')
level_boundaries = 4

## PLOT THE QUANTITIES
## GWs
axd["A"].plot(GWs[:,0], GWs[:, 1], alpha=.6, label='GWs')
axd["A"].legend(loc='lower right')

## 2D PLOT
levels = np.linspace(-level_boundaries, level_boundaries, 25)
pcm = axd["B"].contourf(R, T, AE220.T, levels=levels, cmap='seismic', vmin=-level_boundaries, vmax=level_boundaries)
fig.colorbar(pcm, cax = axd["b"])

fig.savefig(os.path.join(args.save_path, args.name), dpi=300)
