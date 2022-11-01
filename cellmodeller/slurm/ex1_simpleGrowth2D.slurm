#!/bin/bash
#SBATCH --partition=short                       # Partition (short, compute)
#SBATCH --job-name=cellmod_ex1_2D               # Job name
#SBATCH --mail-type=END,FAIL                    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=igors.dubanevics@oist.jp    # Where to send mail	
#SBATCH --ntasks=1                              # Run on a single CPU
#SBATCH --cpus-per-task=1                       # CPUs per tasks
#SBATCH --mem=1gb                               # Job memory request
#SBATCH --time=00:15:00                         # Time limit hrs:min:sec
#SBATCH --output=logs/cellmod_ex1_2D_%j.log     # Standard output and error log
module purge; pwd; hostname; date
echo

echo "Activate Conda environment"
source apps/miniconda/etc/profile.d/conda.sh
conda activate cellmodeller

echo "Running CellModeller on a single CPU core"
python3 github/CellModeller/Scripts/batch.py github/CellModeller/Examples/ex1_simpleGrowth2D.py 0 0

echo "Job finished at"
date