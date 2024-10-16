#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00
#SBATCH --mem=4G
#SBATCH --partition=slow

srun page_rank_pull_parallel --nThreads 8 --nIterations 20 --inputFile input_graphs/test_25M_50M --strategy 4 --granularity 10