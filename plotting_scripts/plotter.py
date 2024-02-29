import matplotlib.pyplot as plt
import os, sys
import numpy as np
import matplotlib.ticker as mticker
from matplotlib import rc
import re 
import plot_info as pinf


#def fx_plot()

# Check if the output directory is provided as a command-line argument
if len(sys.argv) < 4:
    print("Usage: python3 plotter.py home_dir sim_name output_number")
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



outdir = sim_dir+sim_name+"/output-000"+out_number+"/output_directory"

print(outdir)


binac_dir = "/home/jolivera/BINAC_output"
discoverer_dir = "/home/jolivera/Discoverer_data"
python_dir = "/home/jolivera/python"
ns_disc_k2 = discoverer_dir+"/NS_SF_k2"
gr_dir = binac_dir+"/SF_NS/tests/decoupling/NS_GR"








# Now you can use the output_dir variable in your program
print("Looking for simulation directory:", outdir)

#output_dir = ns_disc_k2
thorn = "hydrobase"
quantity = 'rho'

t1,x1,rl1,rl_n1,datax1 = pinf.get_info(thorn,quantity,outdir)
t_rho,rho=pinf.fx_timeseries(t1,x1,datax1)

thorn = "scalarbase"
quantity = 'phi'

os.getcwd()

t1,x1,rl1,rl_n1,datax1 = pinf.get_info(thorn,quantity,outdir)
t_phi,phi=pinf.fx_timeseries(t1,x1,datax1)
