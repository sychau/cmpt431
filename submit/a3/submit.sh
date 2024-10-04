#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00
#SBATCH --mem=4G
#SBATCH --partition=slow

srun page_rank_push --nThreads 8 --nIterations 10 --inputFile input_graphs/roadNet-CA