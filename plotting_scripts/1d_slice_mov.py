#################################################################
#  Extract 1d slice along with star radius. 
#
# 
#    Usage:
#   python3 1d_slice_mov.py computer_name sim_name coordinate adjust_axs? make_movie       
#          computer_name: Specify computer (to add paths)       
#          coordinate: x, y, z                                                                                               
#          adjust_axs: 1 to keep axis in plot fixed, 0 to auto adjust
#          make_movie: r, full, movie
#    Output:
#    {savedatadir}{sim_name}_Psi4_l{l_mode}_m{m_mode}_r{radius}.txt             :  t_ret[M], psi4_lm[M] 
#    {saveplotdir}{sim_name}_Psi4_l{l_mode}_m{m_mode}_r{radius}.pdf 
#################################################################

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


rho_atm = 1e-7 # density threshold to extract surface radius


################################################
 # Set command line arguments
################################################

if len(sys.argv) < 6:
    print("""Usage: python3 1d_slice_mov.py computer_name sim_name coordinate adjust_axs? make_movie
          computer_name: Specify computer (to add paths)
          sim_name: name of the ETK simulation 
          coordinate: x, y, z
          adjust_axs: 1 to keep axis in plot fixed, 0 to auto adjust
          make_movie: r, full, movie
           """)
    sys.exit(1)

# set args

current_computer = sys.argv[1]
sim_name = sys.argv[2]
coordinate = sys.argv[3]
adj_axs = sys.argv[4]
make_movie  = sys.argv[5]


home_dir, sim_dir = pinf.IDcomputer(current_computer)
direc= sim_dir + sim_name+"/output-0000/output_directory" 
current_dir = os.getcwd()
saveplotdir = current_dir+"/vids/"
tmp_dir = saveplotdir+"/"+sim_name



t_re_x = []
re_x   = []  # array to save the radius

j=0


################################################
 # Only get r_e or r_p
################################################
if make_movie == "r":
   thorn = "hydrobase"
   quantity = "rho"

   tk1,xk1,rl1,rln1,datax = pinf.get_info(thorn,quantity,direc,0.0,coordinate)
 
   for itd in range(0,len(tk1),1):
       x_j, f_xj_ti = pinf.get_1d_slice(tk1,xk1,datax,itd,coordinate) 

     # Extract radius of star
       x_surface = pinf.get_star_surface(x_j,f_xj_ti,rho_atm)
       print(f"At t={tk1[itd]}, {coordinate}_surface= {x_surface}")
       
       t_re_x.append(tk1[itd])
       re_x.append(x_surface)

   re_file = f"{current_dir}/plots/re{coordinate}_{sim_name}.txt"
   np.savetxt(re_file, np.column_stack((t_re_x,re_x)),header=f"t re_{coordinate}", fmt ="%.13f")


   os.chdir(current_dir)
   subprocess.run(f"git add plots/*",shell=True)
   subprocess.run(f'git commit -m "r_{coordinate} {sim_name}" ',shell= True)
   subprocess.run("git push",shell=True)


###################################################
 # make movie of 1d slice evolution at different t
###################################################
   
elif make_movie == "yes":
   thorns = ["hydrobase","admbase","scalarbase"]
   quantities = ["rho","lapse","phi"]
   ax_lims = [(-0.0001,0.006),(-0.01,1),(-0.01,0.03)]

   tk1,xk1,rl1,rln1,datax = pinf.get_info(thorns[0],quantities[0],direc,0.0,coordinate)
   if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)
   for itd in range(0,len(tk1),1):
# Plot the data
      fig, axs = plt.subplots(1, len(quantities), figsize=(18, 8))
      for thorn, quantity, ax,lim_ax in zip(thorns,quantities, axs,ax_lims):
         tk1,xk1,rl1,rln1,datax = pinf.get_info(thorns[0],quantities[0],direc,0.0,coordinate)
         x_j, f_xj_ti = pinf.get_1d_slice(tk1,xk1,datax,itd,coordinate) 
         if quantity=="rho":
            x_surface = pinf.get_star_surface(x_j,f_xj_ti,rho_atm)
            print(f"At t={tk1[itd]}, {coordinate}_surface= {x_surface}")
       
            t_re_x.append(tk1[itd])
            re_x.append(x_surface)
# Plot with sorted data
         ax.plot(x_j, f_xj_ti, label=rf'BD $k_0-0.1$, $m_\varphi=1.0e-2$,  t={tk1[itd]/M_to_ms:.4f} ms', color='red')
         ax.set_xlabel(rf'{coordinate} [$M$]', fontsize=14)
         if adj_axs:
            ax.set_ylim(lim_ax)
         ax.set_xlim(-35,35)
         ax.legend()
         ax.set_title(title_coll[quantity], fontsize=16)
         ax.grid(True)
         pinf.apply_second_xaxis_distance(ax)
      j=j+1
      plt.tight_layout()
      plt.savefig(f"{tmp_dir}/p_{j}.png")

################################################
 # only get  1d slice evolution data
################################################
else:
   thorns = ["hydrobase","admbase","scalarbase"]
   quantities = ["rho","lapse","phi"]

   for thorn, quantity in zip(thorns,quantities):
      tk1,xk1,rl1,rln1,datax = pinf.get_info(thorn,quantity,direc,0.0,coordinate)
      for itd in range(0,len(tk1),1):
         x_j, f_xj_ti = pinf.get_1d_slice(tk1,xk1,datax,itd,coordinate) 
         var_file = f"{current_dir}/plots/plots_{sim_name}/{sim_name}_{quantity}_{tk1[itd]:.4f}.txt"
         np.savetxt(var_file, np.column_stack((x_j,f_xj_ti)),header=f"{coordinate} {quantity}", fmt ="%.13f")

     # Extract radius of star
         if quantity=="rho":
            x_surface = pinf.get_star_surface(x_j,f_xj_ti,rho_atm)
            print(f"At t={tk1[itd]}, {coordinate}_surface= {x_surface}")
       
            t_re_x.append(tk1[itd])
            re_x.append(x_surface)


   re_file = f"{current_dir}/plots/re{coordinate}_{sim_name}.txt"
   np.savetxt(re_file, np.column_stack((t_re_x,re_x)),header=f"t re_{coordinate}", fmt ="%.13f")


   os.chdir(current_dir)
   subprocess.run(f"git add plots/*",shell=True)
   subprocess.run(f'git commit -m "1d slice frames {sim_name}" ',shell= True)
   subprocess.run("git push",shell=True)


