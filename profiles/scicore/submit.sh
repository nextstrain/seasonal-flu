#!/bin/sh

#SBATCH --output=log/%j.out                 # where to store the output ( %j is the JOBID )
#SBATCH --error=log/%j.err                  # where to store error messages

# activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
export PATH="/scicore/home/neher/GROUP/bin:/scicore/home/neher/neher/julia-1.7.2/bin:/scicore/home/neher/neher/.julia/bin:"$PATH
conda activate nextstrain

#Test
export AUGUR_MINIFY_JSON=1
export AUGUR_RECURSION_LIMIT=10000

{exec_job}


