
import  subprocess 
import shutil 
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import plot_info as pinf

#directories
sims_dir = "/home/jolivera/simulations"
python_dir= '/home/jolivera/python'
saveplotdir= '/home/jolivera/python/plots'
datadir = "/home/jolivera/ET_files/ID_data_ETK/plotting_scripts/plots/"
datfile_dir = '/home/jolivera/python/dat_files/'

direc = sims_dir+"/BDColl/output-0000/output_directory"   #BD no A mid


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

j=0
thorns = ["hydrobase","admbase","scalarbase"]
quantities = ["rho","lapse","phi"]
ax_lims = [(-0.01,0.005),(-0.01,1),(-0.01,0.2)]

tk1,xk1,rl1,rln1,datax = pinf.get_info(thorns[0],quantities[0],direc)

tmp_dir = python_dir+"/tmp/"
os.mkdir(tmp_dir)

for itd in range(0,len(tk1),10):
#for itd in range():

# Plot the data
  fig, axs = plt.subplots(1, 3, figsize=(18, 8))

  for thorn, quantity, ax,lim_ax in zip(thorns,quantities, axs,ax_lims):
    tk1,xk1,rl1,rln1,datax = pinf.get_info(thorn,quantity,direc)

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
    ax.plot(xj_sorted, f_xi_tj_sorted, label=rf'BD $k_0-0.6$  t={tj[0]/M_to_ms:.4f}', color='red')

    ax.set_xlabel(r'X [$M$]', fontsize=14)
    ax.set_ylim(lim_ax)
    #ax.set_ylabel(r'$\{}_c(t)/\{}_c(0)-1$'.format(quantity,quantity), fontsize=14)
    ax.legend()
    ax.set_title(title_coll[quantity], fontsize=16)
    ax.grid(True)
    pinf.apply_second_xaxis_distance(ax)

# Show the plot
  j=j+1
  plt.tight_layout()
  plt.savefig(f"{tmp_dir}p_{j}.png")
  #plt.show()


os.chdir(tmp_dir)

ffcommand = "ffmpeg -framerate 20 -i p_%d.png -c:v libx264 -pix_fmt yuv420p out.mp4" 
mvff = f"mv out.mp4 {python_dir}"

subprocess.run(ffcommand,shell=True)
subprocess.run(mvff,shell=True)

os.chdir(python_dir)

shutil.rmtree(tmp_dir)

