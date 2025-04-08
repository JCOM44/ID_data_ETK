
import  subprocess 
import shutil 
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import plot_info as pinf

#define plot info directory 

title_coll = {
    "rho": "Central rest-mass density",
    "lapse": "Central lapse",
    "phi": "Central scalar field"
    
}

# constants, in SI
G = 6.673e-11       # m^3/(kg s^2)
c = 299792458       # m/s
M_sol = 1.98892e30  # kg
# convertion factors
M_to_ms = 1./(1000*M_sol*G/(c*c*c))
M_to_density = c**5 / (G**3 * M_sol**2) # kg/m^3



if len(sys.argv) < 4:
    print("""Usage: python3 1d_tshot.py computer_name sim_name coordinate
          computer_name: Specify computer (to add paths)
          sim_name: name of the ETK simulation 
          coordinate: x, y, z
           """)
    sys.exit(1)

current_computer = sys.argv[1]
sim_name = sys.argv[2]
coordinate = sys.argv[3]


home_dir, sim_dir = pinf.IDcomputer(current_computer)
direc= sim_dir + sim_name+"/output-0000/output_directory" 
current_dir = os.getcwd()
saveplotdir = current_dir+"/vids/"
tmp_dir = saveplotdir+"/"+sim_name


j=0
thorns = ["hydrobase","admbase","scalarbase"]
quantities = ["rho","lapse","phi"]
ax_lims = [(-0.0001,0.006),(-0.01,1),(-0.01,0.03)]

tk1,xk1,rl1,rln1,datax = pinf.get_info(thorns[0],quantities[0],direc,0.0,coordinate)

if not os.path.isdir(tmp_dir):
     os.mkdir(tmp_dir)

for itd in range(0,len(tk1),10):
  if itd>len(tk1):
     print("Finally")
     break

# Plot the data
  fig, axs = plt.subplots(1, 3, figsize=(18, 8))

  for thorn, quantity, ax,lim_ax in zip(thorns,quantities, axs,ax_lims):
    tk1,xk1,rl1,rln1,datax = pinf.get_info(thorn,quantity,direc,0.0,coordinate)

    print(f"Getting 1d-{coordinate} slice at t = {tk1[itd]}")
    x_index = datax[:,8] == tk1[itd]

    if coordinate == "x": 
       f_x_ti = np.vstack(  (datax[x_index,8],  datax[x_index,9]  , datax[x_index,12]  ))
    if coordinate == "y": 
       f_x_ti = np.vstack(  (datax[x_index,8],  datax[x_index,10]  , datax[x_index,12]  ))
    if coordinate == "z": 
       f_x_ti = np.vstack(  (datax[x_index,8],  datax[x_index,11]  , datax[x_index,12]  ))

    #print(tk1)

    tj = (f_x_ti[0]).tolist()
    xj = (f_x_ti[1]).tolist()
    f_xi_tj = (f_x_ti[2]).tolist()

    # Convert lists back to numpy arrays for sorting
    xj = np.array(xj)
    f_xi_tj = np.array(f_xi_tj)

# Sort the arrays based on xj
    sorted_indices = np.argsort(xj)

# Reorder both xj and f_xi_tj based on sorted indices
    xj_sorted = xj[sorted_indices]
    f_xi_tj_sorted = f_xi_tj[sorted_indices]


# Save sorted data
    
    var_file = f"{current_dir}/plots/{sim_name}_{quantity}_{tj[0]:.4f}.txt"
    np.savetxt(var_file, np.column_stack((xj_sorted,f_xi_tj_sorted)),header=f"{coordinate} {quantity}", fmt ="%.13f")

# Show the plot


os.chdir(current_dir)
subprocess.run(f"git add plots/*",shell=True)
subprocess.run(f'git commit -m "1d slice frames {sim_name}" ',shell= True)
subprocess.run("git push",shell=True)


