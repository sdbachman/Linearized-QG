#!/bin/bash
### Job Name
#PBS -N QG
#PBS -l walltime=120:00:00
#PBS -q verylong
#PBS -l nodes=1:ppn=16
REPL1

### Run the executable

NPROCS=`wc -l < $PBS_NODEFILE`

echo "Node File:"
echo "----------"
cat  "$PBS_NODEFILE"
echo ""

REPL2
#mpiexec -v -n 36 -display-map --hostfile $PBS_NODEFILE /home/bachman/.linuxbrew/bin/python3 /home/bachman/dedalus/examples/ivp/layered_qg/qg_2.py 

python3 combine_output.py
#/usr/local/bin/qsub start_dedalus.sh

