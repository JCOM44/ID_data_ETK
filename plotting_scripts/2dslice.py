#################################################################
# code for creating a 2d slice in xy and xz plane. 
#
# 
#    Usage:
#     python3 2dslice.py computer_name sim_name var fps
#        var : choose var (rho, phi, etc)
#
#    Output:
#    {sim_name}_{var_name}_frames.zip containing all the 
#       plots to then make a video
#################################################################


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
    print("Usage: python3 2dslice.py computer_name sim_name var fps")
    sys.exit(1)

# Retrieve the output directory from the command-line arguments
current_computer = sys.argv[1]
sim_name = sys.argv[2]
var_name = sys.argv[3]
fps = int(sys.argv[4])

home_dir, sim_dir = pinf.IDcomputer(current_computer)

current_dir = os.getcwd()
saveplotdir = current_dir+"/vids"


def update_plot(it_n, time_steps, rhoxy, rhoxz, scene, grid):
    time = time_steps[it_n]
    fig.suptitle(f'time = {time / M_to_ms:.3f} ms')
    
    # Update rhoxy plot (left subplot)
    rhoxy_ploted = rhoxy.read_on_grid(it_nbr[it_n], grid)
    for artist in scene[0]:
        artist.remove()  # Remove previous contours
    gridp_xy_in_km = np.array(rhoxy_ploted.coordinates_meshgrid()) * 1.477 # transform to km 
    if var_name == "rho":
      new_rhoxy_contour = ax1.contourf(*gridp_xy_in_km, np.log10(rhoxy_ploted.data_xyz))
    elif var_name == "phi":
      new_rhoxy_contour = ax1.contourf(*gridp_xy_in_km, rhoxy_ploted.data_xyz)
    scene[0] = new_rhoxy_contour.collections
    
    # Update rhoxz plot (right subplot)
    rhoxz_ploted = rhoxz.read_on_grid(it_nbr[it_n], grid)
    for artist in scene[1]:
        artist.remove()  # Remove previous contours
    gridp_xz_in_km = np.array(rhoxz_ploted.coordinates_meshgrid()) * 1.477 # transform to km 
    if var_name == "rho":
      new_rhoxz_contour = ax2.contourf(*gridp_xz_in_km, np.log10(rhoxz_ploted.data_xyz))
    elif var_name == "phi":
      new_rhoxz_contour = ax2.contourf(*gridp_xz_in_km, rhoxz_ploted.data_xyz)
    scene[1] = new_rhoxz_contour.collections
    
    # Update subplot titles
    ax1.set_title(f'{var_name}xy at t = {time / M_to_ms:.3f} ms')
    ax2.set_title(f'{var_name}xz at t = {time / M_to_ms:.3f} ms')
    
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

if var_name == "rho":
        rhoxy = testsim.gf.xy.fields.rho
        rhoxz = testsim.gf.xz.fields.rho
elif var_name == "phi":
        rhoxy = testsim.gf.xy.fields.phi1
        rhoxz = testsim.gf.xz.fields.phi1

print("Obtaining file information...")
it_nbr = rhoxy.available_iterations
time_steps = rhoxy.available_times


print(f"Available iterations: {len(it_nbr)}")
print(f"Available timesteps: {len(time_steps)}")

# Create directories to store the frames and output
frames_dir = os.path.join(saveplotdir, 'frames')
os.makedirs(frames_dir, exist_ok=True)

# Create initial grid and plot
grid = UniformGrid([100, 100], x0=[-10, -10], x1=[10, 10])
rho0xy_center = rhoxy.read_on_grid(0, grid)
rho0xz_center = rhoxz.read_on_grid(0, grid)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fig_width * 2, fig_height))
plt.ioff()  # Turn interactive plotting off


# First subplot for rhoxy
grid_xy_in_km = np.array(rho0xy_center.coordinates_meshgrid()) * 1.477 # transform to km 
if var_name == "rho":
      cf1 = ax1.contourf(*grid_xy_in_km, np.log10(rho0xy_center.data_xyz))
elif var_name == "phi":
      cf1 = ax1.contourf(*grid_xy_in_km, rho0xy_center.data_xyz)

plt.colorbar(cf1, ax=ax1)
ax1.set_title(f'{var_name}_xy')
ax1.set_xlabel('x (km)')
ax1.set_ylabel('y (km)')

# Second subplot for rhoxz
grid_xz_in_km = np.array(rho0xz_center.coordinates_meshgrid()) * 1.477 # transform to km 
if var_name == "rho":
      cf2 = ax2.contourf(*grid_xz_in_km, np.log10(rho0xz_center.data_xyz))
elif var_name == "phi":
      cf2 = ax2.contourf(*grid_xz_in_km, rho0xz_center.data_xyz)
plt.colorbar(cf2, ax=ax2)
ax2.set_title(f'{var_name}_xz')
ax2.set_xlabel('x (km)')
ax2.set_ylabel('z (km)')

#scene = [cf.collections]  # Store the initial plot
scene = [cf1.collections, cf2.collections]  # Store the initial plots

# Animation settings
#frames = np.array([1]) 
frames = np.arange(0, len(time_steps), 1)

# Function to save each frame
def save_frame(frame):
    update_plot(frame, time_steps, rhoxy, rhoxz, scene, grid)
    plt.savefig(os.path.join(frames_dir, f'frame_{frame:04d}.png'))
    pbar_save.update(1)

# Progress bar for animation creation and saving
total_frames = len(frames)
with tqdm(total=total_frames, desc="Creating and saving animation", unit="frame") as pbar_save:
    anim = animation.FuncAnimation(
        fig,
        update_plot,
        frames=frames,
        fargs=(time_steps, rhoxy, rhoxz, scene, grid),
        interval=1000 / fps,
        blit=False,
        repeat=False
    )
    
    # Save as GIF
    gif_path = os.path.join(saveplotdir, f"{sim_name}_{var_name}.gif")
    anim.save(gif_path, writer='pillow', fps=fps, progress_callback=lambda i, n: pbar_save.update(1))
    print(f"GIF saved to {gif_path}")
   
# Reset the progress bar for frame saving
    pbar_save.reset()
 
    # Save individual frames
    for i in range(total_frames):
        save_frame(i)
    #save_frame(37)
# Zip the frames
zip_filename = os.path.join(saveplotdir, f'{sim_name}_{var_name}_frames.zip')
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





