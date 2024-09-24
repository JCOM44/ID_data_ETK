import matplotlib.pyplot as plt
import os, sys
import numpy as np
from kuibit import simdir as sd
from kuibit.grid_data import UniformGrid
import matplotlib.animation as animation
import matplotlib.cm as cm
from tqdm import tqdm
import plot_info as pinf

# Constants

G = 6.673e-11       # m^3/(kg s^2)
c = 299792458       # m/s
M_sol = 1.98892e30  # kg
M_to_ms = 1. / (1000 * M_sol * G / (c * c * c))
M_to_density = c**5 / (G**3 * M_sol**2)  # kg/m^3



if len(sys.argv) < 3:
    print("Usage: python3 2dxy.py computer_name sim_name var fps")
    sys.exit(1)

current_computer = sys.argv[1]
sim_name = sys.argv[2]


home_dir, sim_dir = pinf.IDcomputer(current_computer)
dirGw= sim_dir + sim_name 
current_dir = os.getcwd()
saveplotdir = current_dir+"/pdfplots/"

sim = sd.SimDir(dirGw)
psi4 = sim.gws

print(f"Available extraction radii: {psi4.radii[:]}")

radius = psi4.radii[-1] # r = 300
psi4 = (psi4[radius])[(2,2)]
fig = plt.figure()
plt.plot(psi4.t, np.real(psi4.y),  label = r"$\Psi_4^{2,2}$ Real Part at r ="+str(radius))
plt.xlabel(r"$t/M_{sun}$")
plt.ylabel(r"$\Psi_4$")
plt.legend()
pinf.apply_second_xaxis(plt.gca())

figf_name = f"{saveplotdir}{sim_name}_Psi4_l2_m2_r{radius}"
fig.savefig(figf_name+".png")
fig.savefig(figf_name+".pdf")
