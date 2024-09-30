import matplotlib.pyplot as plt
import os, sys
import numpy as np
from kuibit import simdir as sd
from kuibit.grid_data import UniformGrid
import matplotlib.animation as animation
import matplotlib.cm as cm
from tqdm import tqdm
import plot_info as pinf
import subprocess

# Constants

G = 6.673e-11       # m^3/(kg s^2)
c = 299792458       # m/s
M_sol = 1.98892e30  # kg
M_to_ms = 1. / (1000 * M_sol * G / (c * c * c))
M_to_density = c**5 / (G**3 * M_sol**2)  # kg/m^3



if len(sys.argv) < 6:
    print("Usage: python3 GW.py computer_name sim_name l_mode m_mode STT_phi? l_mode_STT m_mode_STT")
    sys.exit(1)

current_computer = sys.argv[1]
sim_name = sys.argv[2]
l_mode = int(sys.argv[3])
m_mode = int(sys.argv[4])
STT_phi = sys.argv[5]
if STT_phi:
   l_mode_STT = sys.argv[6]
   m_mode_STT = sys.argv[7]

home_dir, sim_dir = pinf.IDcomputer(current_computer)
dirGw= sim_dir + sim_name 
current_dir = os.getcwd()
saveplotdir = current_dir+"/pdfplots/"
savedatadir = current_dir+"/plots/"


sim = sd.SimDir(dirGw)
psi4 = sim.gws

r_values = psi4.radii[:]

print(f"{len(r_values)} available extraction radii: {r_values}")

for radius in r_values:
   psi4_in_r = (psi4[radius])[(l_mode,m_mode)]
   
   print(f"Extracting GW at r={radius}")

   t_var = psi4_in_r.t
   var = np.real(psi4_in_r.y)
   dataf_name = f"{savedatadir}{sim_name}_Psi4_l{l_mode}_m{m_mode}_r{radius}.txt"

   np.savetxt(dataf_name, np.column_stack((t_var, var)), header='t {real(Psi4) 2,2}', comments='', fmt='%.18e')
   print(f"Saving file as {dataf_name}")
   
   multipole = f"({l_mode},{m_mode})"

   fig = plt.figure()
   plt.plot(t_var, var,  label = fr"$\Psi_4^{multipole}$ Real Part at r ="+str(radius))
   plt.xlabel(r"$t/M_{sun}$")
   plt.ylabel(r"$\Psi_4$")
   plt.legend()
   pinf.apply_second_xaxis(plt.gca())
	
   figf_name = f"{saveplotdir}{sim_name}_Psi4_l{l_mode}_m{m_mode}_r{radius}"
   #fig.savefig(figf_name+".png")
   fig.savefig(figf_name+".pdf")

   if STT_phi: 
       print(f"Extracting scalar  GW at r={radius}")
       phi_filename =f"mp_phi_l{l_mode_STT}_m{m_mode_STT}_r{radius}0"
       data_phi = np.loadtxt(f"{dirGw}/output-0000/output_directory/{phi_filename}.asc",dtype=np.float64)
       t_phi = data_phi[:,0]   
       phi   = data_phi[:,1]   
   
       np.savetxt(phi_filename+".txt", np.column_stack((t_var, var)), header=f't phi {l_mode_STT},{m_mode_STT}', comments='', fmt='%.18e')
       print(f"Saving file as {phi_filename}")
       fig = plt.figure()
       plt.plot(t_phi, phi,  label = rf"$\varphi^{ {l_mode_STT} , {m_mode_STT} }$ at r ="+str(radius))
       plt.xlabel(r"$t/M_{sun}$")
       plt.ylabel(r"$\varphi$")
       plt.legend()
       pinf.apply_second_xaxis(plt.gca())
	
    #   fig.savefig(saveplotdir+sim_name+phi_filename+".png")
       fig.savefig(saveplotdir+sim_name+phi_filename+".pdf")
 
subprocess.run("git add plots/*",shell=True)
subprocess.run("git add pdfplots/*",shell=True)
subprocess.run(f'git commit -m "GW {sim_name}" ',shell= True)
subprocess.run("git push",shell=True)
