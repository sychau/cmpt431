#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00
#SBATCH --mem=1G
#SBATCH --partition=slow

srun python /home/$USER/cmpt431/a2/test_scripts/submission_validator.pyc \
     --tarPath=/home/$USER/cmpt431/a2/assignment2.tar.gz
