import matplotlib.pyplot as plt
import os, sys
import numpy as np

import matplotlib.animation as animation
import matplotlib.cm as cm

from tqdm import tqdm
import zipfile
import plot_info as pinf


#sys.path.append("/home/jmeneses/ET_files/ID_data_ETK/plotting_scripts/kuibit/")

from kuibit.simdir import SimDir
from kuibit.grid_data import UniformGrid
#Input paramets

if len(sys.argv) < 5:
    print("Usage: python3 2dxy.py computer_name sim_name var fps")
    sys.exit(1)

# Retrieve the output directory from the command-line arguments
current_computer = sys.argv[1]
sim_name = sys.argv[2]
var_name = sys.argv[3]
fps = int(sys.argv[4])

home_dir, sim_dir = pinf.IDcomputer(current_computer)

current_dir = os.getcwd()
saveplotdir = current_dir+"/vids"


def update_plot(it_n, time_steps, rhosfef, phifef, scene, grid):
    time = time_steps[it_n]
    fig.suptitle(f'time = {time / M_to_ms:.3f} ms')
    
    # Update rhosfef plot (left subplot)
    rho_ploted = rhosfef.read_on_grid(it_nbr[it_n], grid)
    for artist in scene[0]:
        artist.remove()  # Remove previous contours
    new_rho_contour = ax1.contourf(*rho_ploted.coordinates_meshgrid(), np.log10(rho_ploted.data_xyz))
    scene[0] = new_rho_contour.collections
    
    # Update phifef plot (right subplot)
    phi_ploted = phifef.read_on_grid(it_nbr[it_n], grid)
    for artist in scene[1]:
        artist.remove()  # Remove previous contours
    new_phi_contour = ax2.contourf(*phi_ploted.coordinates_meshgrid(), phi_ploted.data_xyz)
    scene[1] = new_phi_contour.collections
    
    # Update subplot titles
    ax1.set_title(f'rhosfef at t = {time / M_to_ms:.3f} ms')
    ax2.set_title(f'phifef at t = {time / M_to_ms:.3f} ms')
    
    return scene[0] + scene[1]  # Return all updated artists for blitting

# Constants

outdir = sim_dir+sim_name+"/output-0000/output_directory"

G = 6.673e-11       # m^3/(kg s^2)
c = 299792458       # m/s
M_sol = 1.98892e30  # kg
M_to_ms = 1. / (1000 * M_sol * G / (c * c * c))
M_to_density = c**5 / (G**3 * M_sol**2)  # kg/m^3

cmap = cm.gist_rainbow
fig_width, fig_height = 6, 5

# Load data from simulation directory
testdir = outdir
os.chdir(testdir)

print("creating simulation object...")
testsim = SimDir(testdir)
print("Generated Sim object." )

print("Retrieving 2d field data..." )
rhosfef = testsim.gf.xy.fields.rho


print("Obtaining file information...")
it_nbr = rhosfef.available_iterations
time_steps = rhosfef.available_times

phifef = testsim.gf.xy.fields.phi1
print(f"Available iterations: {len(it_nbr)}")
print(f"Available timesteps: {len(time_steps)}")

# Create directories to store the frames and output
frames_dir = os.path.join(saveplotdir, 'frames')
os.makedirs(frames_dir, exist_ok=True)

# Create initial grid and plot
grid = UniformGrid([100, 100], x0=[-10, -10], x1=[10, 10])
rho0_center = rhosfef.read_on_grid(0, grid)
phi0_center = phifef.read_on_grid(0, grid)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fig_width * 2, fig_height))
plt.ioff()  # Turn interactive plotting off


# First subplot for rhosfef
cf1 = ax1.contourf(*rho0_center.coordinates_meshgrid(), np.log10(rho0_center.data_xyz))
plt.colorbar(cf1, ax=ax1)
ax1.set_title('rhosfef')

# Second subplot for phifef
cf2 = ax2.contourf(*phi0_center.coordinates_meshgrid(), phi0_center.data_xyz)
plt.colorbar(cf2, ax=ax2)
ax2.set_title('phifef')

#scene = [cf.collections]  # Store the initial plot
scene = [cf1.collections, cf2.collections]  # Store the initial plots

# Animation settings
frames = np.arange(0, len(time_steps), 1)

# Function to save each frame
def save_frame(frame):
    update_plot(frame, time_steps, rhosfef, phifef, scene, grid)
    plt.savefig(os.path.join(frames_dir, f'frame_{frame:04d}.png'))
    pbar_save.update(1)

# Progress bar for animation creation and saving
total_frames = len(frames)
with tqdm(total=total_frames, desc="Creating and saving animation", unit="frame") as pbar_save:
    anim = animation.FuncAnimation(
        fig,
        update_plot,
        frames=frames,
        fargs=(time_steps, rhosfef, phifef, scene, grid),
        interval=1000 / fps,
        blit=False,
        repeat=False
    )
    
    # Save as GIF
    gif_path = os.path.join(saveplotdir, f"{sim_name}.gif")
    anim.save(gif_path, writer='pillow', fps=fps, progress_callback=lambda i, n: pbar_save.update(1))
    print(f"GIF saved to {gif_path}")
   
# Reset the progress bar for frame saving
    pbar_save.reset()
 
    # Save individual frames
    for i in range(total_frames):
        save_frame(i)

# Zip the frames
zip_filename = os.path.join(saveplotdir, f'{sim_name}_frames.zip')
with zipfile.ZipFile(zip_filename, 'w') as zipf:
    for root, _, files in os.walk(frames_dir):
        for file in tqdm(files, desc="Zipping frames", unit="file"):
            zipf.write(os.path.join(root, file), 
                       os.path.relpath(os.path.join(root, file), 
                                       os.path.join(frames_dir, '..')))

print(f"Frames saved and zipped in {zip_filename}")

# Optionally, remove the individual frame files to save space
import shutil
shutil.rmtree(frames_dir)
print(f"Individual frame files removed from {frames_dir}")





