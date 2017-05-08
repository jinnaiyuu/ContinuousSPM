#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q c1normal

# Get variables
# This wacky line is required for whatever reasons.
#. /fefs/opt/x86_64/intel/parallel_studio_xe_2017/impi/2017.2.174/bin64/mpivars.sh
#source /home/hal9000/.bashrc
LIBRARY_PATH=$LIBRARY_PATH:/home/hal9000/library/lib
export LIBRARY_PATH
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/hal9000/library/lib
export LD_LIBRARY_PATH

# TODO: choose which data to test
cd ../src/cc
INSTDIR=../../instances
#echo "./pmc -l -i ../R/readydata/$inst.data -c ../R/readydata/$inst.class -o ../R/results/$inst.yuu.results -a 0.05"
#./pmc -l -i ../R/readydata/$inst.data -c ../R/readydata/$inst.class -o ../R/results/$inst.yuu.results -a 0.05
./pmc.opt    -i $INSTDIR/$inst.data -c $INSTDIR/$inst.class -o ../../results/$inst.mahito.a5.patterns -a 0.05


