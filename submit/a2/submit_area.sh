#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00
#SBATCH --mem=1G
#SBATCH --partition=slow

srun ./curve_area_parallel --nThreads 8 --nPoints 1000000000 --coeffA 3.2 --coeffB 7.1 --rSeed 37