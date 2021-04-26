#!/bin/bash
#SBATCH -A m2865
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=16
#SBATCH -c 2
#SBATCH --gpus-per-node=8

# Despite what its name suggests, --gpus-per-task in the examples below only counts the number of GPUs to allocate to the job; it does not enforce any binding or affinity of GPUs to CPUs or tasks.

export SLURM_CPU_BIND="cores"
srun -n 16 /global/cscratch1/sd/gguidi/diBELLA.2D/build_release/./dibella -i /global/cscratch1/sd/gguidi/diBELLA.2D/ecoli_hifi.fasta -k 31 --idxmap ecoli.hifi.idxmap -c 191028 --alph dna --af ecoli.hifi.result -s 1 -O 100000 --afreq 100000 --xa 15
