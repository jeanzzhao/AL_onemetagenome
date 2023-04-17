#!/bin/bash -login
#SBATCH -p med2                # use 'med2' partition for medium priority
#SBATCH -J grist               # name for job
#SBATCH -c 11                   # 4 core, Tessa: use 11
#SBATCH -t 10-00:00:00             # ask for 24h, Tessa: use 4-:00:00:00, 4 days
#SBATCH --mem=200G             # memory (2000 mb = 2gb)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zyzhao@ucdavis.edu
#SBATCH -e jobs/grist.j%j.err
#SBATCH -o jobs/grist.j%j.out

# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# activate your desired conda environment
conda activate grist

# fail on weird errors
set -e
set -x

### YOUR COMMANDS GO HERE ### Tessa: -
genome-grist run conf_ERR.yml summarize_gather summarize_mapping --cores 11 --unlock
genome-grist run conf_ERR.yml summarize_gather summarize_mapping --cores 11 --keep-going --rerun-incomplete

# Print out values of the current jobs SLURM environment variables
env | grep SLURM

# Print out final statistics about resource use before job exits
scontrol show job ${SLURM_JOB_ID}
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
