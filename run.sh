#!/bin/bash -x
#SBATCH --account=171206131221
#SBATCH --job-name="sc"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --time=1:00:00
#SBATCH -o output.txt
#SBATCH -e errlog.txt

cd $SLURM_SUBMIT_DIR

srun -n 16 python do_it_all.py -b 1 -o "testout" -f "../*/*cif"
