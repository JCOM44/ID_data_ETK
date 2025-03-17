
import os, sys
import numpy as np
import subprocess
#sys.path.append("/home/jmeneses/ET_files/ID_data_ETK/plotting_scripts/packaging/src")
#sys.path.append("/home/jmeneses/ET_files/ID_data_ETK/plotting_scripts/matplotlib/lib")

#import matplotlib.pyplot as plt
#import matplotlib.ticker as mticker
#from matplotlib import rc
import re 
import plot_info as pinf


#def fx_plot()

# Check if the output directory is provided as a command-line argument
if len(sys.argv) < 5:
    print("Usage: python3 plotter.py computer_name sim_name var_set  output_number")
    sys.exit(1)

# Retrieve the output directory from the command-line arguments
current_computer = sys.argv[1]
sim_name = sys.argv[2]
var_set = sys.argv[3]
out_number = sys.argv[4]



home_dir, sim_dir =  pinf.IDcomputer(current_computer)

current_dir = os.getcwd()
plot_dir = current_dir+"/plots"


# sET to 1 to print something else
debug_on = 0
old_format = 0


#Define thorns 
if re.search(r"imple",var_set):
    var_list = ["rho","phi","ham"]
    thorn_list = ["hydrobase","scalarbase","jbssn"]

elif re.search(r"GR",var_set):
    var_list = ["rho"]
    thorn_list = ["hydrobase"]

elif re.search(r"ollaps",var_set):
    var_list = ["rho","phi","lapse"]
    thorn_list = ["hydrobase","scalarbase","admbase"]

elif (debug_on ==1):

    var_list = ["phi_x_derivative","phi_rhs_total","kphi_rhs_total"]
    thorn_list = ["scalarevolve","scalarevolve","scalarevolve"]

else:
    var_list = ["rho","phi","lapse","kphi","shift","ml_ham"]
    thorn_list = ["hydrobase","scalarbase","admbase","scalarbase","admbase","ml_admconstraints"]


for j in range(len(thorn_list)):

    #Define arrays to store data
    t_var   = []
    var     = []
    thorn = thorn_list[j]
    quantity = var_list[j]
    print("Getting data for {}-{}".format(thorn,quantity)) 
    savefilename = f"{quantity}_{sim_name}.txt"
    
    for i in range(int(out_number)):
        
        if (old_format == 1):
            outdir = sim_dir+sim_name+"/output-000"+str(i)+f"/{sim_name}"
        else:
            outdir = sim_dir+sim_name+"/output-000"+str(i)+"/output_directory"
	
        t0, file_exist = pinf.check_file(thorn,quantity,f"{plot_dir}/{savefilename}")

        # Now you can use the output_dir variable in your program
        print("Looking for simulation directory:", outdir)

        t1,x1,rl1,rl_n1,datax1 = pinf.get_info(thorn,quantity,outdir,t0)
        if t0<t1[-1]:
                t_var_temp,var_temp=pinf.fx_timeseries(t1,x1,datax1)
                t_var.extend(t_var_temp)
                var.extend(var_temp)

    os.chdir(plot_dir)
    if file_exist:
        with open(savefilename,'a') as f:
              np.savetxt(f,np.column_stack((t_var, var)), fmt='%.18e')
    else:
        np.savetxt(savefilename, np.column_stack((t_var, var)), header='t {}'.format(quantity), comments='', fmt='%.18e')
    print(f"Saving file as {savefilename}.txt")

subprocess.run("git add *",shell=True)
subprocess.run(f'git commit -m "{sim_name}" ',shell= True)
subprocess.run("git push",shell=True)

