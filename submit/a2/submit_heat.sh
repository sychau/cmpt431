#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00
#SBATCH --mem=1G
#SBATCH --partition=slow

srun ./heat_transfer_parallel --nThreads 4 --gSize 1000 --iCX 0.1 --iCY 0.15 --mTemp 600 --tSteps 4000