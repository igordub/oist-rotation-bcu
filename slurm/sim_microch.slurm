#!/bin/bash

# Time flag has been commented out
# Supply time when submitting the job as it depends 
# on the number of simualtion steps requested
# Example: sbatch --time-min=10 sim_microch.slurm src/microch_30x1.py 12000

#SBATCH --partition=short                       # Partition (short, compute)
#SBATCH --job-name=sim_microch                  # Job name
#SBATCH --mail-type=END,FAIL                    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=igors.dubanevics@oist.jp    # Where to send mail	
#SBATCH --ntasks=1                              # Run on a single CPU
#SBATCH --cpus-per-task=1                       # CPUs per tasks
#SBATCH --mem=1gb                               # Job memory request
##SBATCH --time=00:10:00                        # Time limit hrs:min:sec
#SBATCH --output=logs/sim_microch_%j.log        # Standard output and error log

SCRIPT_PATH=$1
NUM_OF_STEPS=$2
module purge; pwd; hostname; date
echo

echo "Activate Conda environment"
source ~/apps/miniconda/etc/profile.d/conda.sh
conda activate cellmodeller

echo "Running CellModeller on a single CPU core"
python3 ~/apps/cellmodeller/Scripts/batch.py "$SCRIPT_PATH" 0 0 "$NUM_OF_STEPS"

echo "Job finished at"
date
