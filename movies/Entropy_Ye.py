import matplotlib.pyplot as plt
import numpy as np
from scidata.quantities.quantities import SimulationAnalysis
from scidata.units.units import units
from scidata.grid.grid import grid
import matplotlib.ticker as tick
import os
import cv2
import argparse
import shutil
parser = argparse.ArgumentParser()
parser.add_argument('--sim', type=str, required=True)
parser.add_argument('--sim-path', type=str, required=False, default=None)
parser.add_argument('--rxlim', type=float, default=None)
parser.add_argument('--keep-frames', action='store_true', default=False)
parser.add_argument('--name', required=True, type=str, help='Name of the movie without extension')
parser.add_argument('--save-path', required=True, type=str, help='path where to save the movie')
args = parser.parse_args()

u = units()

def prepare_plot(fig, xlim_grid, y_lim_GWs, rxlimGWs):
    axd = fig.subplot_mosaic(
        """..AA..
           b.BC.c""",
        #height ratio between rows
        height_ratios = [1, 3.5],
        #withd ratio between columns
        width_ratios = [0.07, 0.17, 0.9, 0.9, 0.17, 0.07],
        gridspec_kw={
            "wspace": 0,
            "hspace": 0,
        }
        )
    ## fix labels and ticks
    ## Labels
    axd["A"].set_xlabel("t-t$_b$ [ms]")
    axd["B"].set_xlabel("x [km]")
    axd["C"].set_xlabel("x [km]")
    axd["A"].set_ylabel(r"$\mathcal{D}h_+$ [cm]")
    axd["B"].set_ylabel("z [km]")
    axd["C"].set_ylabel("z [km]")
    ## Ticks
    axd["A"].tick_params(top=True, labeltop=True,
                        bottom=False, labelbottom=False,
                        left=True, labelleft=True,
                        right=False, labelright=False)
    axd["B"].tick_params(top=False, labeltop=False,
                        bottom=True, labelbottom=True,
                        left=True, labelleft=True,
                        right=False, labelright=False)
    axd["C"].tick_params(top=False, labeltop=False,
                        bottom=True, labelbottom=True,
                        left=False, labelleft=False,
                        right=True, labelright=True)
    axd["b"].tick_params(top=False, labeltop=False,
                        bottom=False, labelbottom=False,
                        left=True, labelleft=True,
                        right=False, labelright=False)
    axd["c"].tick_params(top=False, labeltop=False,
                        bottom=False, labelbottom=False,
                        left=False, labelleft=False,
                        right=True, labelright=True)
    axd["A"].xaxis.set_label_position('top')
    axd["B"].xaxis.set_label_position('bottom')
    axd["B"].xaxis.set_label_position('bottom')
    axd["A"].yaxis.set_label_position('left')
    axd["B"].yaxis.set_label_position('left')
    axd["C"].yaxis.set_label_position('right')

    #axis
    axd["C"].set_xlim(0, xlim_grid)
    axd["C"].set_ylim(-xlim_grid, xlim_grid)
    axd["B"].set_ylim(-xlim_grid, xlim_grid)
    axd["B"].set_xlim(xlim_grid, 0)
    axd["A"].set_xlim(-5, rxlimGWs)
    axd["A"].set_ylim(-y_lim_GWs, y_lim_GWs)

    return axd

def plot_quantities(fig, axd, X, Y, entropy, Ye, Bfield_lines, Vx, Vy, GWs, time, points_skip):
    ## Number of levels in contourf
    levels_Ye = np.linspace(0, 0.5, 25)
    levels_S = np.linspace(1.5, 12, 25)
    ## Plotting the gravitational waves
    axd["A"].plot(GWs[:,0], GWs[:, 1])
    axd["A"].axvline(time, lw=1, ls='--', c='black')
    ## Plotting the Ye
    pcm_ye = axd["B"].contourf(X, Y, Ye.T, levels=levels_Ye,
                               vmin = 0.0, vmax = 0.5, antialiased=True,
                               cmap ='gist_rainbow_r')
    cbar_ye = fig.colorbar(pcm_ye, cax=axd["b"])
    cbar_ye.set_label("Y$_e$")
    cbar_ye.ax.yaxis.set_ticks_position('left')
    cbar_ye.ax.yaxis.set_label_position('left')
    cbar_ye.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    if Bfield_lines is not None:
        bfield = axd["B"].contour(X, Y, Bfield_lines.T, 50, colors = 'black')
        plt.clabel(bfield, inline = True, fontsize=8, fmt='%.2e')
    ## Plotting the entropy
    pcm_s = axd["C"].contourf(X, Y, entropy.T, levels=levels_S,
                              vmin = 1.5, vmax = 12, antialiased=True,
                              cmap = 'gist_rainbow_r')
    axd["C"].quiver(X[points_skip], Y[points_skip], Vx.T[points_skip], Vy.T[points_skip])
    cbar_s = fig.colorbar(pcm_s, cax=axd["c"])
    cbar_s.set_label("S [k$_B$/bry]")

def plot_radii(axd, R_inner, R_PNS, angle):
    X_inner = R_inner * np.sin(angle)
    Y_inner = R_inner * np.cos(angle)
    X_PNS = R_PNS * np.sin(angle)
    Y_PNS = R_PNS * np.cos(angle)
    axd["B"].plot(X_inner, Y_inner, lw=1, ls='--', c='#A52A2A')
    axd["C"].plot(X_inner, Y_inner, lw=1, ls='--', c='#A52A2A')
    axd["B"].plot(X_PNS, Y_PNS, lw=1, ls='--', c='#FFE4C4')
    axd["C"].plot(X_PNS, Y_PNS, lw=1, ls='--', c='#FFE4C4')


##  SIMULATION PARAMETER AND SAVING FOLDER
sim_name = args.sim
sim_path = args.sim_path
output_folder = args.save_path
frame_folder = os.path.join(output_folder, 'frames')
if not os.path.exists(frame_folder):
    os.mkdir(frame_folder)

## INITIALIZE THE SIMULATION
sim = SimulationAnalysis(sim_name, simulation_folder_path=args.sim_path)
# GET THE GRID
sim.ghost.update_ghost_cells(t_l=3, t_r=3)
radius = u.convert_to_km(sim.cell.radius(sim.ghost))
theta = sim.cell.theta(sim.ghost)
gr = grid(sim.dim, radius, theta)
X, Y = gr.cartesian_grid()
sim.ghost.restore_default()
## GET THE RADII
Innercore_radius, indices, time, g_cells = sim.get_innercore_radius(innercore_radius = True, indices = True,
                                                            ret_time = True, ghost_cells = True)
try:
    PNS_radius = sim.get_neutrino_sphere_radius(radius_nua=True)
except:
    print('PNS')
    PNS_radius = sim.get_PNS_radius(PNS_radius = True)
PNS_radius = u.convert_to_km(PNS_radius)
Innercore_radius = u.convert_to_km(Innercore_radius)
time = u.convert_to_ms(time)

## FIND THE STARTING INDEX
index_start = np.argmax(time >= -6)
index_end = np.argmax(time >= 100) + 1

## CUT THE ARRAYS
indices = indices[index_start:index_end]

## DEFINE THE SKIPPING 
skip = (slice(None, None, 4), slice(None, None, 4))

## GET THE GWS
GW_strain = sim.GW_Amplitudes()
GW_strain[:,0] = u.convert_to_ms(GW_strain[:,0])


## DEFINE THE PLOT LIMITS
plot_xlim = 100 ## km
if args.rxlim is None:
    GWs_xlim = GW_strain[-1, 0]
    GWs_ylim = np.ceil( np.max( np.abs( GW_strain[:, 1] ) ) + 10 )
else:
    GWs_xlim = args.rxlim
    GWs_ylim = np.ceil( np.max( np.abs( GW_strain[:np.argmax(GW_strain[:, 0] >= 100), 1] ) ) + 10 )
if GWs_ylim >= 100:
    GWs_ylim = 100

## CREATE THE FIGURE
figure = plt.figure()
frame_name = 1

## FILE LIST

file_list = sim.file_list_hdf()
frame_list = []

## CYCLE OVER THE NEEDED FRAMES
for index in indices:
    ## OPEN THE DATA
    ##UPDATE THE GHOST CELLS
    sim.ghost.update_ghost_cells(t_l=3, t_r=3)
    data = sim.open_h5(file_list[index])
    S = sim.entropy(data)
    ye = sim.Ye(data)
    try:
        streamlines = sim.stream_function(data, 'xz')
    except:
        None
    Vr = sim.radial_velocity(data)
    Vtheta = sim.theta_velocity(data)
    sim.close_h5(data)
    vx = Vr * np.sin(theta)[:, None] +  Vtheta * np.cos(theta)[:, None]
    vy = Vr * np.cos(theta)[:, None] -  Vtheta * np.sin(theta)[:, None]
    ## PREPARE THE FIGURE
    figure.clear()
    axs = prepare_plot(figure, plot_xlim, GWs_ylim, GWs_xlim)
    ## PLOT THE QUANTITIES
    plot_quantities(figure, axs, X, Y, S, ye, streamlines, vx, vy, GW_strain, time[index], skip)
    ## PLOT THE RADII
    sim.ghost.restore_default()
    plot_radii(axs, sim.ghost.remove_ghost_cells_radii(Innercore_radius[..., index], sim.dim, t_l=4, t_r=4),
            sim.ghost.remove_ghost_cells_radii(PNS_radius[..., index], sim.dim, t_l=4, t_r=4), theta)
    ## SAVE THE FIGURE
    figure.set_size_inches(10, 15)
    save_path = os.path.join(frame_folder, str(frame_name) + '.png')
    frame_list.append(save_path)
    figure.savefig(save_path, bbox_inches='tight', dpi=300)
    frame_name += 1


frame = cv2.imread(frame_list[0])
height, width, channels = frame.shape

# Define the video codec and create a VideoWriter object
fourcc = cv2.VideoWriter_fourcc(*'vp80')
video_writer = cv2.VideoWriter(os.path.join(output_folder, args.name + '.webm'), fourcc, 10.0, (width, height))

# Write each frame to the video
for frame_name in frame_list:
    frame = cv2.imread(frame_name)
    video_writer.write(frame)

# Release the video writer and print a success message
video_writer.release()
if not args.keep_frames:
    shutil.rmtree(frame_folder)

