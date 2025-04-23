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

# 3.9 works with all my crap, unlikekly that they change latest in 3.9
module purge
module load python/3/3.9/latest

# source env path 
ENV_PATH="/discofs/bg-phys-02/ETK/python/ETKanalysis/bin/activate"
source "$ENV_PATH"

# example: 
#  sbatch pyjob.sh GW.py server1 sim_name 2 0 1 1 1
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
