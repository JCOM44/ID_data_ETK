import matplotlib.pyplot as plt
import os
import numpy as np
from kuibit.grid_data import UniformGrid
import matplotlib.animation as animation
import matplotlib.cm as cm
from kuibit.simdir import SimDir
from tqdm import tqdm

# Constants
sims_dir = "/home/jolivera/simulations"
saveplotdir = '/home/jolivera/python/plots'
GR_dir = sims_dir + "/BDCollN/output-0000/output_directory"

G = 6.673e-11       # m^3/(kg s^2)
c = 299792458       # m/s
M_sol = 1.98892e30  # kg
M_to_ms = 1. / (1000 * M_sol * G / (c * c * c))
M_to_density = c**5 / (G**3 * M_sol**2)  # kg/m^3

cmap = cm.gist_rainbow
fig_width, fig_height = 6, 5

# Updated update function that avoids using the deprecated `collections` attribute
def update_plot(it_n, time_steps, quant_toplot, scene, grid):
    time = time_steps[it_n]
    fig.suptitle(f'time = {time / M_to_ms:.3f} ms')

    # Recompute new data to plot
    ploted = quant_toplot.read_on_grid(it_nbr[it_n], grid)

    # Clear the previous contour plot completely
    for artist in scene[0]:
        artist.remove()  # Remove previous contours

    # Plot the new data
    new_contour_plot = plt.contourf(*ploted.coordinates_meshgrid(), np.log10(ploted.data_xyz))
    scene[0] = new_contour_plot.collections  # Store the new contour artists

    return new_contour_plot.collections  # Return the updated artists for blitting

# Load data from simulation directory
testdir = GR_dir
os.chdir(testdir)
testsim = SimDir(testdir)

rhosfef = testsim.gf.xy.fields.rho
it_nbr = rhosfef.available_iterations
time_steps = rhosfef.available_times

print(f"Available iterations: {len(it_nbr)}")
print(f"Available timesteps: {len(time_steps)}")

# Create initial grid and plot
grid = UniformGrid([100, 100], x0=[-15, -10], x1=[15, 10])

rho0_center = rhosfef.read_on_grid(0, grid)
fig = plt.figure(figsize=(fig_width, fig_height))
plt.ioff()  # Turn interactive plotting off

cf = plt.contourf(*rho0_center.coordinates_meshgrid(), np.log10(rho0_center.data_xyz))
plt.colorbar(cf)
scene = [cf.collections]  # Store the initial plot's artists

# Animation settings
fps = 20
frames = np.arange(0, len(time_steps), 1)

# Progress bar for animation frames
with tqdm(total=len(frames), desc="Animating", unit="frame") as pbar:
    anim = animation.FuncAnimation(
        fig,
        update_plot,
        frames=frames,
        fargs=(time_steps, rhosfef, scene, grid),
        interval=1000 / fps,
        blit=False,  # Blitting might be tricky with contour plots
        repeat=False
    )

    for _ in frames:
        pbar.update(1)

# Progress bar for saving animation
total_frames_to_save = len(frames)

def progress_callback(frame_idx, total_frames):
    pbar_save.update(1)

# Make a second progress bar for saving the animation
with tqdm(total=total_frames_to_save, desc="Saving animation", unit="frame") as pbar_save:
    anim.save(f'{saveplotdir}/test.mp4', writer='ffmpeg', fps=fps, progress_callback=progress_callback)