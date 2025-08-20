import  subprocess 
import shutil 
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import plot_info as pinf
from kuibit.simdir import SimDir                 #generates an environment which shows available plots 

# constants, in SI
G = 6.673e-11       # m^3/(kg s^2)
c = 299792458       # m/s
M_sol = 1.98892e30  # kg
# convertion factors
M_to_ms = 1./(1000*M_sol*G/(c*c*c))
M_to_density = c**5 / (G**3 * M_sol**2) # kg/m^3


if len(sys.argv) < 4:
    print("Usage: python3 horizon.py computer_name sim_name pxrint_qlm")
    sys.exit(1)

current_computer = sys.argv[1]
sim_name = sys.argv[2]
print_qlm = sys.argv[3]

home_dir, sim_dir = pinf.IDcomputer(current_computer)
data_dir= sim_dir + sim_name
current_dir = os.getcwd()
saveplotdir = current_dir+"/plots/"

# Number of columns correspond to those given in the file BH_diagnostics.ah1.gp: 

cols_to_extract = [1, 7, 25, 26, 27]

horizonfile = f"{saveplotdir}{sim_name}_horizon.txt"
horizon_data = [] 
out_number  = pinf.get_out_number(data_dir) 

for i in range(int(out_number)):
    input_file  = f"{data_dir}/output-{i:04d}/output_directory/BH_diagnostics.ah1.gp"
    if not os.path.exists(input_file):
        print(f"Warning: File {input_file} not found, skipping.")
        continue
    
    try:
        data = np.loadtxt(input_file, comments="#", usecols=cols_to_extract)
        horizon_data.append(data)
        print(f"Read {input_file} successfully.")
        print(horizon_data)
    except Exception as e:
        print(f"Error reading {input_file}: {e}")

horizon_data1 = np.vstack(horizon_data)

# Save to new file
np.savetxt(horizonfile, horizon_data1, fmt="%.18e", header="time mean_radius area m_irreducible areal_radius")

subprocess.run(f"git add plots/*",shell=True)
subprocess.run(f'git commit -m "horizon {sim_name}" ',shell= True)
subprocess.run("git push",shell=True)
