
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
    print("""Usage: python3 1d_slice_mov.py computer_name sim_name adjust_axs?
          computer_name: Specify computer (to add paths)
          sim_name: name of the ETK simulation 
          adjust_axs: 1 to keep axis in plot fixed, 0 to auto adjust
           """)
    sys.exit(1)

current_computer = sys.argv[1]
sim_name = sys.argv[2]
adj_axs = sys.argv[3]

home_dir, sim_dir = pinf.IDcomputer(current_computer)
direc= sim_dir + sim_name+"/output-0000/output_directory" 
current_dir = os.getcwd()
saveplotdir = current_dir+"/vids/"
tmp_dir = saveplotdir+"/"+sim_name


j=0
thorns = ["hydrobase","admbase","scalarbase"]
quantities = ["rho","lapse","phi"]
ax_lims = [(-0.0001,0.006),(-0.01,1),(-0.03,0.4)]

tk1,xk1,rl1,rln1,datax = pinf.get_info(thorns[0],quantities[0],direc,0.0)

if not os.path.isdir(tmp_dir):
     os.mkdir(tmp_dir)

for itd in range(0,len(tk1),1):

# Plot the data
  fig, axs = plt.subplots(1, 3, figsize=(18, 8))

  for thorn, quantity, ax,lim_ax in zip(thorns,quantities, axs,ax_lims):
    tk1,xk1,rl1,rln1,datax = pinf.get_info(thorn,quantity,direc,0.0)

    print("Getting 1d-x slice at t = {}".format(tk1[itd]))
    x_index = datax[:,8] == tk1[itd]

    f_x_ti = np.vstack(  (datax[x_index,8],  datax[x_index,9]  , datax[x_index,12]  ))

    print(tk1)

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

# Plot with sorted data
    ax.plot(xj_sorted, f_xi_tj_sorted, label=rf'BD $k_0-0.6$  t={tj[0]/M_to_ms:.4f} ms', color='red')

    ax.set_xlabel(r'X [$M$]', fontsize=14)
    if adj_axs:
       ax.set_ylim(lim_ax)
    ax.set_xlim(-35,35)
    #ax.set_ylabel(r'$\{}_c(t)/\{}_c(0)-1$'.format(quantity,quantity), fontsize=14)
    ax.legend()
    ax.set_title(title_coll[quantity], fontsize=16)
    ax.grid(True)
    pinf.apply_second_xaxis_distance(ax)

# Show the plot
  j=j+1
  plt.tight_layout()
  plt.savefig(f"{tmp_dir}/p_{j}.png")


os.chdir(current_dir)
subprocess.run(f"git add vids/*",shell=True)
subprocess.run(f'git commit -m "1d slice frames {sim_name}" ',shell= True)
subprocess.run("git push",shell=True)


