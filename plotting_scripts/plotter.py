#################################################################
# Timeseries extraction code. 
#
# Extracts timeseries for a given quantity. 
# 
#    Usage:
# python3 plotter.py computer_name sim_name var_set coordinate 
#    var_set  : 
#               GR  (densitty, ham)
#               simple ( GR + phi)
#               collapse ( simple + lapse) 
#               all  ( collapse + shift and  kphi)
#    Output:
#    {quantity}_{sim_name}.txt   : t [M], vars[code units] 
#################################################################


import os, sys
import numpy as np
import subprocess
import re 
import plot_info as pinf


# Check if the output directory is provided as a command-line argument
if len(sys.argv) < 5:
    print("Usage: python3 plotter.py computer_name sim_name var_set coordinate")
    sys.exit(1)

# Retrieve the output directory from the command-line arguments
current_computer = sys.argv[1]
sim_name = sys.argv[2]
var_set = sys.argv[3]
coordinate  = sys.argv[4]

home_dir, sim_dir =  pinf.IDcomputer(current_computer)
out_number  = pinf.get_out_number(sim_dir+sim_name) 

current_dir = os.getcwd()
plot_dir = current_dir+"/plots"


# sET to 1 to print something else
debug_on = 0
old_format = 0


#Define thorns 
if re.search(r"imple",var_set):
    var_list = ["rho","phi","ham"]
    thorn_list = ["hydrobase","scalarbase","jbssn"]
    redux_list = ["norm2"]

elif re.search(r"GRcoll",var_set):
    var_list = ["rho","ham","lapse"]
    thorn_list = ["hydrobase","jbssn","admbase"]

elif re.search(r"GR",var_set):
    var_list = ["rho","ham"]
    thorn_list = ["hydrobase","jbssn"]

elif re.search(r"ollaps",var_set):
    var_list = ["rho","phi","lapse","ham"]
    thorn_list = ["hydrobase","scalarbase","admbase","jbssn"]
    redux_list = ["norm2"]

else:
    var_list = ["rho","phi","lapse","kphi","shift","ham"]
    thorn_list = ["hydrobase","scalarbase","admbase","scalarbase","admbase","jbssn"]


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

        if quantity=="ham":
            t_var_temp,var_temp=pinf.redux_timeseries(thorn,quantity,outdir,"norm2")                
        else:
           t1,x1,rl1,rl_n1,datax1 = pinf.get_info(thorn,quantity,outdir,t0)
           if t0<t1[-1]:
                t_var_temp,var_temp=pinf.fx_timeseries(t1,x1,datax1)
        t_var.extend(t_var_temp)
        var.extend(var_temp)

    os.chdir(plot_dir)
    if file_exist:
        if quantity=="ham":
           np.savetxt(savefilename, np.column_stack((t_var, var)), header='t {}'.format(quantity), comments='', fmt='%.18e')
        else:
           with open(savefilename,'a') as f:
              np.savetxt(f,np.column_stack((t_var, var)), fmt='%.18e')
    else:
        np.savetxt(savefilename, np.column_stack((t_var, var)), header='t {}'.format(quantity), comments='', fmt='%.18e')
    print(f"Saving file as {savefilename}.txt")

subprocess.run("git add *",shell=True)
subprocess.run(f'git commit -m "{sim_name}" ',shell= True)
subprocess.run("git push",shell=True)

print("COMPLETED SUCCESFULLY.")

