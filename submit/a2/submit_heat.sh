#!/bin/bash
#

executable=$1
nThreads=$3
nIterations=$5
inputFiles=$7

# Display the arguments
echo "Argument 1: $arg1"
echo "Argument 2: $arg2"

#SBATCH --cpus-per-task=$nThreads
#SBATCH --time=10:00
#SBATCH --mem=4G
#SBATCH --partition=slow

srun ./page_rank_pull_parallel --nThreads $nThreads --nIterations $nIterations --inputFile $inputFile