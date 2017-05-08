#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q c1normal
LIBRARY_PATH=$LIBRARY_PATH:/home/hal9000/library/lib
export LIBRARY_PATH
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/hal9000/library/lib
export LD_LIBRARY_PATH

# Get variables
# This wacky line is required for whatever reasons.
#. /fefs/opt/x86_64/intel/parallel_studio_xe_2017/impi/2017.2.174/bin64/mpivars.sh
#source /home/hal9000/.bashrc

# TODO: choose which data to tes
cd ../src/cc

#echo "./pmc -l -i ../R/readydata/$inst.data -c ../R/readydata/$inst.class -o ../R/results/$inst.yuu.results -a 0.05"
./pmc.opt -l $thre -i ../R/readydata/$inst.data -c ../R/readydata/$inst.class -o ../../results/$inst.yuu-$thre.a5.results -a 0.05
#./pmc    -i ../R/readydata/$inst.data -c ../R/readydata/$inst.class -o ../R/results/$inst.mahito.results -a 0.05


