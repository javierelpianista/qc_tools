#!/bin/bash

#SBATCH --job-name=rdx
#SBATCH --mem=18gb
#SBATCH --tasks-per-node=8
#SBATCH --time=2:00:00
#SBATCH --output=%x-%j.out
#SBATCH --exclusive

echo "Compute node: $HOSTNAME"
echo "Job ID: $SLURM_JOBID"
echo "Started at:"
date

module purge
module load gcc/8.2.0 openmpi/3.1.3-intel

JOBNAME=rdx

LOCALDIR=/mnt/kolos2-io1/home/jgarcia/qc_tools/examples/rdx/test1
SCRATCHDIR=/scratch/local/jgarcia/$SLURM_JOB_ID
mkdir -p $SCRATCHDIR
cp -r $LOCALDIR/* $SCRATCHDIR/
cd $SCRATCHDIR

/data/sw/SAPT2020.1/bin/sapt-fastdf rdx
cat *.timer
rm -r $SCRATCHDIR

echo "Done at:"
hostname
date
