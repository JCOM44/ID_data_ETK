#################################################################
# GW extraction code. 
#
# Extracts Psi4_l_m and if in STT, also extracts phi_l_m
# 
#    Usage:
# python3 GW.py computer_name sim_name l_mode m_mode STT_phi? l_mode_STT m_mode_STT
#  STT_phi? : boolean 1 or 0
#
#    Output:
#    {savedatadir}{sim_name}_Psi4_l{l_mode}_m{m_mode}_r{radius}.txt             :  t_ret[M], psi4_lm[M] 
#    {saveplotdir}{sim_name}_Psi4_l{l_mode}_m{m_mode}_r{radius}.pdf 
#    {savedatadir}{sim_name}_mp_phi_l{l_mode_STT}_m{m_mode_STT}_r{radius}0.txt  : t_ret[M], phi_lm[M] 
#    {savedatadir}{sim_name}_mp_phi_l{l_mode_STT}_m{m_mode_STT}_r{radius}0.pdf 
#################################################################

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


################################################
 # Set command line arguments
################################################

if len(sys.argv) < 7 and len(sys.argv) > 5:
    print("Usage: python3 GW.py computer_name sim_name l_mode m_mode STT_phi? l_mode_STT m_mode_STT")
    sys.exit(1)



current_computer = sys.argv[1]
sim_name = sys.argv[2]
l_mode = int(sys.argv[3])
m_mode = int(sys.argv[4])
STT_phi = False

if len(sys.argv)>5:
   STT_phi = sys.argv[5]
   l_mode_STT = sys.argv[6]
   m_mode_STT = sys.argv[7]
else:
   print("This is just GR")


################################################
#       set directories
################################################

home_dir, sim_dir = pinf.IDcomputer(current_computer)
dirGw       = sim_dir + sim_name 
out_number  = pinf.get_out_number(dirGw) 
current_dir = os.getcwd()
saveplotdir = current_dir+"/pdfplots/"
savedatadir = current_dir+"/plots/"



################################################
# Extract GW
################################################

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
# Calculate the difference8
   t_var_adjusted = t_var - radius
# Boolean mask to select only positive values of t_var - radius
   positive_mask = t_var_adjusted > 0
# Use the mask to filter the t_var and var arrays
   t_var_filtered = t_var_adjusted[positive_mask]
   var_filtered = var[positive_mask]

################################################
# Plot Psi4
################################################

   fig = plt.figure()
   plt.plot(t_var_filtered,var_filtered,  label = fr"$\Psi_4^{multipole}$ Real Part at r ="+str(radius))
   plt.xlabel(r"$t_{ret}[M]$")
   plt.ylabel(r"$\Psi_4$")
   plt.legend()
   pinf.apply_second_xaxis(plt.gca())
	
   figf_name = f"{saveplotdir}{sim_name}_Psi4_l{l_mode}_m{m_mode}_r{radius}"
   #fig.savefig(figf_name+".png")
   fig.savefig(figf_name+".pdf")



################################################
#  Extract scalar field 
################################################

   if STT_phi:
       t_var = []
       var   = [] 
       seen_t_values = set()  # Set to track unique t_var values


       phi_filename =f"mp_phi_l{l_mode_STT}_m{m_mode_STT}_r{radius}0"
       for i in range(int(out_number)):
          print(f"Extracting scalar at r={radius}")
          data_phi = np.loadtxt(f"{dirGw}/output-{i:04d}/output_directory/{phi_filename}.asc",dtype=np.float64)
          t_phi = data_phi[:,0]   
          phi   = data_phi[:,1]  

# Calculate the retarded time t-r_ext
          t_phi_adjusted = t_phi - radius
# Boolean mask to select only positive values of t_var - radius, so that t_0 = 0
          positive_mask = t_phi_adjusted > 0
# Use the mask to filter the t_var and var arrays
          t_var_filtered = t_phi_adjusted[positive_mask]
          var_filtered = phi[positive_mask]

 # Append only unique values
          for t, v in zip(t_var_filtered, var_filtered):
            if t not in seen_t_values:
                seen_t_values.add(t)
                t_var.append(t)
                var.append(v)

       
       new_phi_filename = f"{savedatadir}{sim_name}_{phi_filename}.txt"
       np.savetxt(new_phi_filename, np.column_stack((t_var, var)), header=f't phi {l_mode_STT},{m_mode_STT}', comments='', fmt='%.18e')
       print(f"Saving file as {new_phi_filename}")

################################################
#  Plot scalar field 
################################################

       fig = plt.figure()
       plt.plot(t_var, var,  label = rf"$\varphi^{ {l_mode_STT} , {m_mode_STT} }$ at r ="+str(radius))
       plt.xlabel(r"$t_{ret} [M]$")
       plt.ylabel(r"$\varphi$")
       plt.legend()
       pinf.apply_second_xaxis(plt.gca())
	
       fig.savefig(saveplotdir+sim_name+phi_filename+".pdf")
 
################################################
#  Upload to git repo 
################################################

subprocess.run("git add plots/*",shell=True)
subprocess.run("git add pdfplots/*",shell=True)
subprocess.run(f'git commit -m "GW {sim_name}" ',shell= True)
subprocess.run("git push",shell=True)
