#!/bin/bash
# export OMP_NUM_THREADS=1
nohup mpiexec -n 4 ./debugcode > out.debugcode.txt &