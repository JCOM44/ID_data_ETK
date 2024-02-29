
import os, sys
import numpy as np
#sys.path.append("/home/jmeneses/ET_files/ID_data_ETK/plotting_scripts/packaging/src")
#sys.path.append("/home/jmeneses/ET_files/ID_data_ETK/plotting_scripts/matplotlib/lib")

#import matplotlib.pyplot as plt
#import matplotlib.ticker as mticker
#from matplotlib import rc
import re 
import plot_info as pinf


#def fx_plot()

# Check if the output directory is provided as a command-line argument
if len(sys.argv) < 4:
    print("Usage: python3 plotter.py computer_name sim_name output_number")
    sys.exit(1)

# Retrieve the output directory from the command-line arguments
current_computer = sys.argv[1]
sim_name = sys.argv[2]
out_number = sys.argv[3]

if re.search(r"wicky",current_computer):
    home_dir='/home/jolivera'  
    sim_dir= home_dir+'/simulations/'
    print("Computer identified as zwicky.")
elif re.search(r"iscovere",current_computer):
    home_dir='/home/jmeneses'
    sim_dir='/discofs/bg-phys-02/ETK/simulations/'
    print("Computer identified as Discoverer.")
elif re.search(r"inac",current_computer):
    home_dir='/home/tu/tu_tu/tu_pelol01' 
    sim_dir= '/beegfs/work/workspace/ws/tu_tupelol01-NS_SF-0/sims/'
    print("Computer identified as BinaC.") 
else:
    home_dir="NULL"
    print("Computer not recognized")
print("\nSetting central directory:       {}".format(home_dir))
print("Setting simulation directory:    {}".format(sim_dir))


current_dir = os.getcwd()
plot_dir = current_dir+"/plots"

#Define thorns 
var_list = ["rho","phi"]
thorn_list = ["hydrobase","scalarbase"]


for j in range(len(thorn_list)):
     
    #Define arrays to store data
    t_var   = []
    var     = []
    thorn = thorn_list[j]
    quantity = var_list[j]
    
    for i in range(int(out_number)):
                
        outdir = sim_dir+sim_name+"/output-000"+str(i)+"/output_directory"
        # Now you can use the output_dir variable in your program
        print("Looking for simulation directory:", outdir)

        t1,x1,rl1,rl_n1,datax1 = pinf.get_info(thorn,quantity,outdir)
        t_var_temp,var_temp=pinf.fx_timeseries(t1,x1,datax1)

        t_var.extend(t_var_temp)
        var.extend(var_temp)

    os.chdir(plot_dir)
    np.savetxt('{}_{}.txt'.format(quantity,sim_name), np.column_stack((t_var, var)), header='t {}'.format(quantity), comments='', fmt='%f')




