from scivis.color_plot_1D.create_2D_plot import prepare_plot, set_labels, plot_quantities, plt
from scidata.quantities.quantities import SimulationAnalysis
from scidata.units.units import units
import argparse, os

u = units()

parser = argparse.ArgumentParser()
parser.add_argument("--sim-name", required=True, help="Name of the simulation to analyse")
parser.add_argument("--sim-path", required=False, default=None, help="Path of the selected simulation")
parser.add_argument("--rxlim", required=False, default=None, type=float, help="Right x-axis limit in ms, if empty full GW is plotted")
parser.add_argument("--name", required=True, help="Plot name")
parser.add_argument("--save-path", required=True, help="Plot path")
parser.add_argument("--logscale", required=False, action='store_true')
parser.add_argument("--time", required=False, default=None, type=float, help="Time to plot after the bounce in ms, if empty last step is plotted.")

args = parser.parse_args()

sim = SimulationAnalysis(args.sim_name, None, args.sim_path)
time = []
if args.time is None:
    file_list = reversed(sim.file_list_hdf())
    for file in file_list:
        data = sim.open_h5(file)
        time.append(sim.time(data))
        unph_cells = sim.errors(data)
        if unph_cells.sum() > 0:
            sim.close_h5(data)
            break
        sim.close_h5(data)
    
else:
    file = sim.find_file_from_time(args.time, True)



radius = u.convert_to_km(sim.cell.radius(sim.ghost))
theta = sim.cell.theta(sim.ghost)
if args.rxlim is not None:
    rxlim = args.rxlim
else:
    rxlim = radius[-1]

fig, ax = prepare_plot((radius[0], rxlim), (theta[0], theta[-1]), args.logscale, False, (12.5, 10))
cbar, ax = plot_quantities(radius, theta, unph_cells, fig, ax, 0, 1, 3, False, 'gray_r')
ax["A"].set_title(r't$_\mathrm{end}$-t = ' + str((time[0] - time[-1])[0]))
set_labels('r [km]', r'$\theta$ [rad]', 'umphysical cells', ax, cbar)
fig.savefig(os.path.join(args.save_path, args.name))
#plt.show()
