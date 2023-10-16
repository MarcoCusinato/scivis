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
parser.add_argument("--maxr", required=False, default=None, type=float, help="Right radius limit in km, if empty full spacial profile is plotted")
parser.add_argument("--maxt", required=False, default=None, type=float, help="Right time limit in ms, if empty full sim time is plotted")
parser.add_argument('--disable-log-scale-r', action='store_true', default=False, help="Makes the y-axis scale linear")
parser.add_argument('--BV_freq_lim', nargs='+', required=False, type=float, default=[1e-3, 1e9], help="Limits of the BV frequency, default 1e-3, 1e9")
parser.add_argument("--name", required=True, help="Plot name")
parser.add_argument("--save-path", required=True, help="Plot path")


args = parser.parse_args()

   

## Get the data
sim = SimulationAnalysis(args.sim_name, None, args.sim_path)
time, radius, BV_freq = sim.BV_frequency_profile()
radius = u.convert_to_km(radius)
time = u.convert_to_ms(time)
fig = plt.figure(figsize=(12.5, 10))
ax = fig.subplot_mosaic(
        """A.a""",
        width_ratios = [0.97, 0.01, 0.02],
        gridspec_kw={'hspace':0} )


if args.maxr is None:
    maxr = radius[-1]
else:
    maxr = args.maxr

if args.maxt is None:
    maxt = time[-1]
else:
    maxt = args.maxt


if maxt >= 1000:
    ax['A'].set_xlabel(r't-t$_b$ [s]')
    ax['A'].set_xlim(u.convert_to_s(-5), u.convert_to_s(maxt))
    time = u.convert_to_s(time)
else:
    ax['A'].set_xlabel(r't-t$_b$ [ms]')
    ax['A'].set_xlim(-5, maxt)

ax['A'].set_ylim(radius[0], maxr)
if not args.disable_log_scale_r:
    ax['A'].set_yscale('log')

ax['A'].set_ylabel('R [km]')

pcm = ax['A'].pcolormesh(time, radius, np.abs(BV_freq), shading='gouraud', cmap='viridis',
                        norm=colors.LogNorm(vmin=args.BV_freq_lim[0], vmax=args.BV_freq_lim[1]))

cbar = fig.colorbar(pcm, cax=ax['a'])
cbar.set_label(r'$|\omega_{BV}^2|$ [Hz$^2$]')



fig.savefig(os.path.join(args.save_path, args.name))

