from scidata.quantities.quantities import SimulationAnalysis
from scivis.color_plot_1D.create_1D_plot import set_up_plot
import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--sim-name", required=True, help="Name of the simulation to analyse")
parser.add_argument("--sim-path", required=False, default=None, help="Path of the selected simulation")
parser.add_argument("--rxlim", required=False, default=None, type=float, help="Right x-axis limit in ms, if empty full time profile is plotted")
parser.add_argument("--tylim", required=False, default=None, type=float, help="Top limit of the y-axis in km, if empty full radial profile is plotted")
parser.add_argument("--name", required=True, help="Plot name")
parser.add_argument("--save-path", required=True, help="Plot path")
parser.add_argument("--vmin", required=False, default=None, type=float, help="Lower value of the colormap")
parser.add_argument("--vmax", required=False, default=None, type=float, help="Upper value of the colormap")

args = parser.parse_args()

sim = SimulationAnalysis(args.sim_name, None, args.sim_path)
fig = set_up_plot(sim, 'Ye', args.rxlim, args.tylim, False, r'Y$_\mathrm{e}$',
                  'viridis', True, args.vmin, args.vmax)
fig.savefig(os.path.join(args.save_path, args.name))