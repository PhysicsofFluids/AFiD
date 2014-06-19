#!/bin/bash
#n this job an hybrid MPI-OpenMP program is started.
# This job will run on two nodes, 8 processes on each node, 
# each process with 4 OpenMP threads making a total of 64 threads.
#
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --job-name=Ra1e8_pencil
#
# Job will take at most 1 hour, 10 minutes and 20 seconds wallclock time
#
#SBATCH -t 24:00:00
#
# Here the shell script starts. 
# go to the working directory:
#
cd $SLURM_SUBMIT_DIR
#
# and call the MPI program
# 
export OMP_NUM_THREADS=1
srun ./boutnp

