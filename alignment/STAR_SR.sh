#! /bin/bash


### run STAR array ###
FILES=($(ls -1 *_1.txt.gz))

# get size of array
NUMFASTQ=${#FILES[@]}

# now submit to SLURM
if [ $NUMFASTQ -ge 0 ]; then
	sbatch --array=1-$NUMFASTQ STAR_SR.sbatch
fi
