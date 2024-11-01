#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00
#SBATCH --mem=1G
#SBATCH --partition=slow

srun python /home/$USER/cmpt431/a2/test_scripts/curve_area_tester.pyc \ 
     --execPath=/home/$USER/cmpt431/a2/curve_area_parallel  \ 
     --scriptPath=/home/$USER/cmpt431/a2/test_scripts/curve_area_evaluator.pyc