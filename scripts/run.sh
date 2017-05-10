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
INSTDIR=../../instances

if [ $method = "mahito" ]
then
    ./pmc.opt -i $INSTDIR/$inst.data -c $INSTDIR/$inst.class -o ../../results/$inst.${method}.a5.patterns -a 0.05
else
    thre=`echo $method | awk -F '-' '{print $2}'`
    ./pmc.opt -l $thre -i $INSTDIR/$inst.data -c $INSTDIR/$inst.class -o ../../results/$inst.${method}.a5.patterns -a 0.05    
fi


