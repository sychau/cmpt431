#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00
#SBATCH --mem=4G
#SBATCH --partition=slow

srun page_rank_push_parallel_atomic --nThreads 8 --nIterations 20 --inputFile /scratch/input_graphs/test_25M_50M --strategy 4 --granularity 10
srun page_rank_push_parallel_atomic --nThreads 8 --nIterations 20 --inputFile /scratch/input_graphs/test_25M_50M --strategy 4 --granularity 100
srun page_rank_push_parallel_atomic --nThreads 8 --nIterations 20 --inputFile /scratch/input_graphs/test_25M_50M --strategy 4 --granularity 1000
srun page_rank_push_parallel_atomic --nThreads 8 --nIterations 20 --inputFile /scratch/input_graphs/test_25M_50M --strategy 4 --granularity 2000