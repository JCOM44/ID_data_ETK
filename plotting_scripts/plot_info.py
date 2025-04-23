#import matplotlib.pyplot as plt
import os, sys, re
import numpy as np
#import matplotlib.ticker as mticker
#from matplotlib import rc



# Plot basic stuff


################################################
 # Define constants and conversion factors
################################################
# constants, in SI
G = 6.673e-11       # m^3/(kg s^2)
c = 299792458       # m/s
M_sol = 1.98892e30  # kg
# convertion factors
M_to_ms = 1./(1000*M_sol*G/(c*c*c))
M_to_density = c**5 / (G**3 * M_sol**2) # kg/m^3



  
################################################
# Function to ID computer.
# input  : 
#        current_copmputer : str, the name of the computer (or something like it)
# output :
#        home_dir          : str, home dir directory
#        sim_dir           : str, simulation directory  
################################################

def IDcomputer(current_computer):

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
	    sim_dir= '/beegfs/work/workspace/ws/tu_pelol01-NS_JBSSN-0/simulations/'
	    print("Computer identified as BinaC.") 
	else:
	    home_dir="NULL"
	    print("Computer not recognized")
	print("\nSetting central directory:       {}".format(home_dir))
	print("Setting simulation directory:    {}".format(sim_dir))
	return home_dir, sim_dir

################################################
# Function to identify the number of output dirs. 
# input  : 
#        dir_name: directory name
#	
# output :
#        out_number: output number
################################################

def get_out_number(output_dir):
# Reg exp to match folders like "output-0000", "output-0012", etc.
   pattern = re.compile(r'^output-\d+$')
# put all dirs with the pattern in a list
   matching_dirs = [
      name for name in os.listdir(output_dir)
      if os.path.isdir(os.path.join(output_dir, name)) and pattern.match(name)
   ]
   return len(matching_dirs)


################################################
# Function to set_tick size in plots.
# input  : 
#        ax : current ax
#	
# output :
################################################

def set_tick_sizes(ax, major, minor):
    for l in ax.get_xticklines() + ax.get_yticklines():
        l.set_markersize(major)
    for tick in ax.xaxis.get_minor_ticks() + ax.yaxis.get_minor_ticks():
        tick.tick1line.set_markersize(minor)
        tick.tick2line.set_markersize(minor)
    ax.xaxis.LABELPAD=10.
    ax.xaxis.OFFSETTEXTPAD=10.

################################################
# Function to check if file already exists to avoid reading everything again
# input  : 
#        thorn     : str, thorn name, e.j. "hydrobase"
#        quantity  : str, quantity name, e.j. "rho" 
#        file_path : str, path where the file should be
# output :
#        returns the max t value that has been calcualted in the file, 
#        if file doesnt exist returns zero to start there 
################################################

def check_file(thorn,quantity,file_path):
	print("Check if file exists...")
	if os.path.exists(file_path):
		data = np.loadtxt(file_path,skiprows=1)
		return np.max(data[:,0]), True  #returns value of t0
	else:
		print("File does not exist (yet)")
		return 0.0, False
		

################################################
# Function to get info on the simulation properties for a given quantity
# input  : 
#        thorn      : str, thorn name, e.j. "hydrobase"
#        quantity   : str, quantity name, e.j. "rho" 
#        folder     : str, path where the file should be
#        t0         : float, value where we begin
#        coordinate : str, choose coord (x,y,z)
# output :
#        t          : array, values in t higher than t0
#        x_p        : array, points in chosen coordinate (x,y,z) 
#        rl         : array, refinement levels (0,1,2,3,....)
#        rl_n       : int  , number of refinement levels   #maybe this is not necessary 
#        data_x     : array, f(t,xp)
################################################


def get_info(thorn,quantity, folder,t0,coordinate="x"):
        print("Looking for files in the folder: {}".format(folder))
        os.chdir(folder)
        filex = f"{thorn}-{quantity}.{coordinate}.asc"
        print("Opening file: {}.....".format(filex))
        datax = np.loadtxt(filex, comments='#')
        print("Reading file....")
        
        # Number of iterations
        it = np.unique(datax[:, 0])
        it_n = len(it)
        print("Number of iterations:", it_n)
        
        # Time values
        t = np.unique(datax[:, 8])
        # X values
        if coordinate == "x": 
           x_p = np.unique(datax[:, 9])
        if coordinate == "y": 
           x_p = np.unique(datax[:, 10])
        if coordinate == "z": 
           x_p = np.unique(datax[:, 11])
        # Refinement levels	
        rl = np.unique(datax[:, 2])
        rl_n = len(rl)
        print("Total number of refinement levels:", rl_n)

        if t0<t[-1]:
            t=t[t>t0]
            t_n = len(t)
            print("Number of different time values:", t_n)

        # Points
            x_p_n = len(x_p)
            print("Total number of points:", x_p_n)


            points_per_rl = []
            rl_max_point = []
            for i in range(rl_n):
                x_in_rl = np.unique(datax[datax[:, 2] == rl[i], 9])
                points_in_rl = len(x_in_rl)
                print("Number of points in refinement level", i, ":", points_in_rl)
                rl_max_point.append(np.max(x_in_rl))
                points_per_rl.append(points_in_rl)
       # rl_max_point.append(0.0)
        
        return t,x_p,rl,rl_n,datax


################################################
# Function to get TIMESRIES for a specific quantity, needs as input the output generated by get_info
# input  : 
#        t           : array,  t array coming from function get_info(..)
#        x_p         : array,  x_p array coming from function get_info(..) 
#        datax       : array,  data_x array coming from get_info()..
#        ixd         : int  ,  index in x, choose value in x_p 
#        coordinate  : str  ,  choose coord (x,y,z)
# output :
#        time_values : array, values in t higher than t0 
#        f_xt_values : array, f(t,x=x[ixd])
################################################

def fx_timeseries(t,x_p,datax,ixd=0, coordinate="x"):     #index value of x as input
    #create output lists
    print(os.getcwd())
    t_n = len(t)
    time_values = []
    f_xt_values = []
    if ixd==0:
        print(f"Calculating timeseries for {coordinate} = 0.0")
    else:
        print(f"Calculating timeseries for {coordinate} = {x_p[ixd]}")
    print(f"Starting at  t = {t[0]}")
 # create filter for time steps
    for j in range(t_n): 
        t_index = datax[:,8] == t[j]
 #create new array only with x and data columns for given t
        f_x_ti = np.vstack(  (datax[t_index,8],  datax[t_index,9]  , datax[t_index,12]  ))  #now x=f_x_ti[0][:] and f(x)=f_x_ti[1][:]
 #create filter for space points
        if ixd==0:
           x_index = f_x_ti[1][:] == 0.0
        else:
           x_index = f_x_ti[1][:] == x_p[ixd]

 # save t, x and f(x,t) in a list (use lists to improve efficiency when extending)
        tj = (f_x_ti[0][x_index]).tolist()
        f_xi_tj = (f_x_ti[2][x_index]).tolist()

 #append values
        time_values.extend(tj)
        f_xt_values.extend(f_xi_tj)
        if(j==np.round(1/8*t_n) or j==np.round(1/4*t_n)) or j==np.round(3/8*t_n) or  j==np.round(1/2*t_n) or j==np.round(5/8*t_n) or j==np.round(3/4*t_n) or j==np.round(7/8*t_n):
                print("Progress: {} %".format(j/t_n *100))
    print("Done...!")
    return time_values,f_xt_values



################################################
# Function to get TIMESRIES for a specific reduced quantity (usually max, min, norm2, etc)
# input  : 
#        thorn      : str, thorn name, e.j. "hydrobase"
#        quantity   : str, quantity name, e.j. "rho" 
#        folder     : str, path where the file should be
#        redux      : str, reduction (maximum,minimum,norm2..)
# output :
#        time_values : array, values in t higher than t0 
#        f_xt_values : array, f(t)
################################################

def redux_timeseries(thorn,quantity,folder,redux):
        os.chdir(folder)
        filex = f"{thorn}-{quantity}.{redux}.asc"
        print("Opening file: {}.....".format(filex))
        datax = np.loadtxt(filex, comments='#')
        print("Reading file....")
        time_values =  datax[:,1]
        f_values    =  datax[:,2]
        return time_values, f_values


################################################
# Function to get 1D slice for a specific quantity, needs as input the output generated by get_info
# input  : 
#        tk           : array,  t array coming from function get_info(..)
#        x_p          : array,  x_p array coming from function get_info(..) 
#        datax        : array,  data_x array coming from get_info()..
#        itd          : int  ,  index in t, choose value in tk 
#        coordinate   : str  ,  choose coord (x,y,z)
# output :
#        xj_sorted    : array, values in coordinate 
#        f_xi_tj_sort : array, f(t=t[ixd],x)
################################################

def get_1d_slice(tk1, xk1, datax, itd, coordinate):

    print(f"Getting 1d-{coordinate} slice at t = {tk1[itd]}")

    t_index = datax[:,8] == tk1[itd] # get all values at fixed time t_i = tk1[itd]

       # get data  as t,coordinate,f(t,coordinate) 
    if coordinate == "x": 
       f_x_ti = np.vstack(  (datax[t_index,8],  datax[t_index,9]  , datax[t_index,12]  ))
    if coordinate == "y": 
       f_x_ti = np.vstack(  (datax[t_index,8],  datax[t_index,10]  , datax[t_index,12]  ))
    if coordinate == "z": 
       f_x_ti = np.vstack(  (datax[t_index,8],  datax[t_index,11]  , datax[t_index,12]  ))

       # split into t_i, x_j, f(x_j,t_i) 
    tj = (f_x_ti[0]).tolist()    # t_i should be all the same
    xj = (f_x_ti[1]).tolist()    # array of {x,y,z} values
    f_xi_tj = (f_x_ti[2]).tolist()

    # Convert lists back to numpy arrays for sorting
    xj = np.array(xj)
    f_xi_tj = np.array(f_xi_tj)
    # Sort the arrays based on xj
    sorted_indices = np.argsort(xj)
    # Reorder both xj and f_xi_tj based on sorted indices
    xj_sorted = xj[sorted_indices]
    f_xi_tj_sorted = f_xi_tj[sorted_indices]

    return xj_sorted, f_xi_tj_sorted


################################################
# Function to get the position of the star surface
# for the chosen coordinate from  1d slice data
# input  : 
#        x_p          : array,  x_p array 
#        f_xi_tj      : array,  array coming from 1d slice f(t=t[idt],x)
#        rho_atm      : int  ,  index in t, choose value in tk 
# output :
#        x_surface    : float, coordinate value at surface in km 
################################################

def get_star_surface(xp,f_xi_tj,rho_atm):
       condition = (f_xi_tj < rho_atm) & (xp >= 0) # identify where density is smaller than rho_atm
       indices   = np.where(condition)[0]                
       if len(indices)>0:
          surface_index   = indices[0]
          x_surface       = xp[surface_index]*1.477
       else:
          x_surface       = 0.0				  # no star found found
       return x_surface 

################################################
# Function to make plot of timeseries rho(t) and phi(t)
################################################


def plot_funcs(t_rho,rho,t_phi,phi,plot_dir=os.getcwd()):
    fig, ax = plt.subplots(1,2,figsize=(16, 10))

    fig.subplots_adjust(top=0.85, bottom=0.16, left=0.11,right=0.97)

    #ax = fig.add_subplot(1,1,1)
    xlim = (0,1)
    ylim = (-0.003,0)
    
    # Plot density

    ax[0].plot(t_rho,rho, linestyle='-.', color='blue', label = r"$DEF \, beta = -5$")                   #BINAC


    # Plot scalar field

    ax[1].plot(t_phi,phi, linestyle='-', color='blue', label = r"$BIN, \, k_0 = 0.001$")


    #ax.grid(True)

    # plot properties
    #ax.set_xlim(xlim)
    #ax.set_ylim(ylim) 

    #plt.title('PSD', fontsize=fontsize)

    ax[0].set_xlabel(r't [M]')
    ax[0].xaxis.set_major_locator(mticker.MaxNLocator(7))
    ax[0].xaxis.set_minor_locator(mticker.MaxNLocator(14))
    ax[0].xaxis.grid(False)
    ax[1].set_xlabel(r't [M]')
    ax[1].xaxis.set_major_locator(mticker.MaxNLocator(7))
    ax[1].xaxis.set_minor_locator(mticker.MaxNLocator(14))
    ax[1].xaxis.grid(False)
    ax2 = ax[0].twiny()
    ax2.set_xlabel(r't [ms]')
    ax2.set_xlim((ax[0].get_xlim()[0]/M_to_ms, ax[0].get_xlim()[1]/M_to_ms))
    ax2.xaxis.set_major_locator(mticker.MaxNLocator(7))
    ax2.xaxis.set_minor_locator(mticker.MaxNLocator(14))
    ax3 = ax[1].twiny()
    ax3.set_xlabel(r't [ms]')
    ax3.set_xlim((ax[1].get_xlim()[0]/M_to_ms, ax[1].get_xlim()[1]/M_to_ms))
    ax3.xaxis.set_major_locator(mticker.MaxNLocator(7))
    ax3.xaxis.set_minor_locator(mticker.MaxNLocator(14))
    ax[0].set_ylabel(r'$\rho_c(t)$')
    ax[0].yaxis.set_major_locator(mticker.MaxNLocator(8))
    ax[0].yaxis.set_minor_locator(mticker.MaxNLocator(10))
    ax[0].yaxis.grid(True)
    set_tick_sizes(ax[0], 8, 8)
    ax[1].legend(loc='upper right', fontsize=fontsize)
    ax[1].set_ylabel(r'$\phi_c(t)$')
    ax[1].yaxis.set_major_locator(mticker.MaxNLocator(8))
    ax[1].yaxis.set_minor_locator(mticker.MaxNLocator(10))
    ax[1].yaxis.grid(True)
    set_tick_sizes(ax[1], 8, 8)
    #ax[1].legend(loc='upper right', fontsize=fontsize)
    ax[0].set_title(r'Central density')
    ax[1].set_title(r'Central scalar field')



    name = "plot1"
    #os.chdir(plots_dir)
    #fig.savefig("phi_min.jpeg")
    os.chdir(plot_dir)
    plt.savefig("{}.png".format(name))




# 2D plotting

def save_frame(frame):
    plt.savefig(os.path.join(frames_dir, f'frame_{frame:04d}.png'))
    pbar_save.update(1)


################################################
# Function to add a second time axis in
# timeseries, this on is in ms
################################################

def apply_second_xaxis(ax):
   ax2=ax.twiny()
   ax2.set_xlabel(r't [ms] ')
   ax2.set_xlim((ax.get_xlim()[0]/M_to_ms, ax.get_xlim()[1]/M_to_ms))

################################################
# Function to add a second coordinate axis in
# spatial profiles, this one is in km 
################################################

def apply_second_xaxis_distance(ax):
   ax2=ax.twiny()
   ax2.set_xlabel(r'x [km] ')
   ax2.set_xlim((ax.get_xlim()[0]*1.477, ax.get_xlim()[1]*1.477))



