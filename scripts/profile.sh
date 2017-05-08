#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q c1normal

LIBRARY_PATH=$LIBRARY_PATH:/home/hal9000/library/lib
export LIBRARY_PATH
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/hal9000/library/lib
export LD_LIBRARY_PATH

cd ../src/cc

#GMON_OUT_PREFIX=gmon.out
#export GMON_OUT_PREFIX

./pmc.prf -l 0.95 -i ../R/readydata/$inst.data -c ../R/readydata/$inst.class -o ../R/results/$inst.yuu.results -a 0.05

gprof ./pmc.prf gmon.out > $inst.yuu.prof
