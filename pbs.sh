#!/bin/sh
#PBS -l nodes=1:ppn=32
#PBS -l walltime=00:05:00
#PBS -q express
#PBS -M nils.kohl@studium.uni-erlangen.de -m abe
#PBS -N SiWiR-ex04
#PBS -o $PBS_JOBNAME.out -e $PBS_JOBNAME.err

. /etc/profile.d/modules.sh
module load openmpi/1.6.5-ib
module load gcc/4.8.2

cd /home/stud/me31kove/Documents/siwir_git/siwir_ex04
make clean
make


mpirun -np 4 ./heat 50 50 10 0.0001 10 0.02 0.5 0.5 0

