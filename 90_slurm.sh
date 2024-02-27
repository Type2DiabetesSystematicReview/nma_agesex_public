#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=none   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --job-name=reghba1c    # some descriptive job name of your choice
#SBATCH --output=%x-%j.out      # output file name will contain job name + job ID
#SBATCH --error=%x-%j.err       # error file name will contain job name + job ID
#SBATCH --partition=nodes        # which partition to use, default on MARS is “nodes"
#SBATCH --time=0-01:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=20G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --nodes=1               # number of nodes to allocate, default is 1
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=4       # number of processor cores to be assigned for each task, default is 1, increase for multi-threading runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node

############# LOADING MODULES (optional) #############
#module load apps/xxx
#module load libs/xxx
module load libs/gsl
module load apps/R/4.3.0/gcc-13.1.0+lapack-3.10.0+blas-3.10.0+x11

############# MY CODE ############ specify network, whether common or independent estimates and fixed or random
Rscript 06b_fit_hba1c_reg.R 1
