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

ah=1
dx=0.02
sdd= SimDir(data_dir)
sdd_hor = sdd.horizons
sdd_hor.available_apparent_horizons

print(f"QuasiLocalMeasures horizons available: {sdd_hor.available_qlm_horizons}")



if ah not in sdd_hor.available_apparent_horizons:
            raise ValueError(f"Apparent horizons {ah} is not available")

horizon = sdd_hor.get_apparent_horizon(ah).ah

current_horizon = sdd_hor.get_apparent_horizon(ah)
time_found = current_horizon.ah.cctk_iteration.t
print(f"Horizon found at {time_found[0] / M_to_ms} ms")

if print_qlm==1:
        horizon_qlm = sdd_hor.get_qlm_horizon(0)

        time_qlm = time_found

        irr_mass = horizon_qlm["irreducible_mass"](time_qlm)

        print(
                    f"""\
        QuasiLocalMeasures:
        Time:                     {time_qlm:4.5f}
        Irreducible mass:         {irr_mass:4.5f}
        Christodoulou mass:       {horizon_qlm['mass'](time_qlm):4.5f}
        Angular momentum:         {horizon_qlm['spin'](time_qlm):4.5f}"""
        )

plt.ylabel(f"Radius of horizon [km]")
plt.xlabel("Time")
plt.plot(horizon.mean_radius*1.477, label="Mean radius")
#plt.plot(horizon.min_radius*1.477, label="Min radius")
#plt.plot(horizon.max_radius*1.477, label="Max radius")
plt.legend()

if dx:
    #logger.debug("Adding resolution y axis")
    plt.twinx()
    plt.plot(horizon.mean_radius/dx)
    #plt.plot(horizon.min_radius/dx )
   # plt.plot(horizon.max_radius/dx )
    plt.ylabel("Number of points on radius")

pinf.apply_second_xaxis(plt.gca())

horizonfile = f"{saveplotdir}{sim_name}_horizon.txt"

print(horizon.mean_radius)

np.savetxt(horizonfile,np.column_stack((horizon.mean_radius)),fmt="%.13f")

plt.savefig(f"{saveplotdir}horizon_{sim_name}.pdf")
subprocess.run(f"git add plots/*",shell=True)
subprocess.run(f'git commit -m "horizon {sim_name}" ',shell= True)
subprocess.run("git push",shell=True)
