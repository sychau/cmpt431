#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --partition=slow
#SBATCH --mem=5G

srun ./heat_transfer_parallel --iCX 0.15 --iCY 0.1 --tSteps 1000 --gSize 4000 --mTemp 600
