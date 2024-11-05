#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=slow
#SBATCH --mem=5G

srun ./curve_area_parallel --coeffA 1.2 --coeffB 0.8 --nPoints 4000000000 --rSeed 129
