#!/bin/bash
#SBATCH --partition=normal                      # "normal" partition
#SBATCH --nodes=1                               # we should only need one node for this!
#SBATCH --cpus-per-task=24                      # 8 CPUs will probably suffice.
#SBATCH --mem-per-cpu=16GB                      # we'll need an assload of memory for dealing with these huge matrices.
#SBATCH --time=0-04:00:00                       # Automatically set experiments to run for five days.
#SBATCH --mail-type=ALL
#SBATCH --output=./profiles/hopper/%x.%j.out
#SBATCH --error=./profiles/hopper/%x.%j.err

#SBATCH --mail-user=apizzime@gmu.edu


module load gnu10
module load python
export PYTHONPATH=~/lib/python:$PYTHONPATH

python profile.models.IC.py $1 $2 $3 $4 $5 24 1
