#!/bin/bash
#SBATCH --account=bg-phys-02
#SBATCH --qos=bg-phys-02

#SBATCH --job-name=py_Job
#SBATCH --time=12:00:00
#SBATCH --output=out.out
#SBATCH --error=err.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cn


# 3.9 works with all my crap, unlikekly that they change latest in 3.9
module purge
module load python/3/3.9/latest
module load ffmpeg/6/6.0.0-amd64-static 


# source env path 
#ENV_PATH="/discofs/bg-phys-02/ETK/python/ETKanalysis/bin/activate"
ENV_PATH="/valhalla/projects/bg-phys-02/JoseC/py_envs/ETKanalysis/bin/activate"
source "$ENV_PATH"

# examples:
#   -- To get Gravitational waves
#                                             l m sf_modes l_sf  m_sf
#  sbatch pyjob.sh GW.py server_name sim_name 2 0   1 	    1     1
#
#   -- To get spatial profiles                          coordinate   adjust_plot   mode: r, full, movie
#  sbatch pyjob.sh 1d_slice_mov.py server_name sim_name   x             1             r
#
#   -- To get timeseries                            mode(simple, GR, collapse, all)     
#   sbatch pyjob.sh plotter.py server_name sim_name   simple                 
#
#   -- To get 2D slice
#   

run_script() {
    script=$1
    shift
    echo "Running $script : $@"
    python3 $script "$@"
}


# get parameters from a file, specially useful if we want to run all scripts
if [[ "$1" == "--from-file" ]]; then
    param_file="$2"
    while read -r line; do
        [[ "$line" =~ ^#.*$ || -z "$line" ]] && continue  # skip comments and empty lines
        run_script $line
    done < "$param_file"
else
    run_script "$@"
fi
