#!/bin/bash

# Assign command-line arguments to variables
executable=$1
nThreads=$3
nIterations=$5
inputFiles=$7

# Display the arguments
echo "Executable: $executable"
echo "Number of Threads: $nThreads"
echo "Number of Iterations: $nIterations"
echo "Input Files: $inputFiles"

# Submit job with dynamic parameters using sbatch
sbatch --cpus-per-task=$nThreads \
       --time=10:00 \
       --mem=4G \
       --partition=slow \
       --wrap="srun $executable --nThreads $nThreads --nIterations $nIterations --inputFile $inputFiles"
