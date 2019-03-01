#!/bin/bash -x
#SBATCH --account=171206131221
##SBATCH -p debug
#SBATCH --job-name="sc"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
##SBATCH --ntasks=64
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --time=1:00:00
#SBATCH -o output.txt
#SBATCH -e errlog.txt

cd $SLURM_SUBMIT_DIR

module purge
module load PrgEnv/python/gcc/2.7.11
module load PrgEnv/mpi/openmpi/gcc/3.0.0
export PYTHONPATH=/u/cb/cr/felixtherrien/p2env/lib/python2.7/site-packages:$PYTHONPATH
source /u/cb/cr/felixtherrien/p2env/bin/activate

srun -n 16 python do_it_all.py -b 1 -o "testout" -f "../*/*cif"