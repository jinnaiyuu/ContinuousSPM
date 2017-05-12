#!/bin/bash

testsummary () {
    vsamples="100 1000 10000 100000"
    vfeatures=20
#    vsamples=200
#    vfeatures="10 50 100"
#    vfeatures="10"
    rs="0.1 0.5"
#    rs="0.1"
    jids=""
    results=""
    method="jinnai-0.95"
    summary=../results/summary/${method}.a5
    for samples in $vsamples
    do
	for features in $vfeatures
	do
	    for r0 in $rs
	    do
		instance="synth_${samples}_${features}_${r0}"
		RESULT=../results/$instance.${method}.a5.stat
		id=`qsub -v inst=$instance,method=${method} -e "../results.$instance.${method}.a5.e" -o "$RESULT"  run.sh | awk '{print $3}'`
		jids="$jids,$id"
		results="${results} ${RESULT}"
		done
	done
    done
    jids=`echo $jids | sed -r 's/^.//'`
#    results=`echo $results | sed -r 's/^.{2}//'`
    echo "jids=$jids"
    echo "results=$results"
    qsub -v results="$results",summary="$summary" \
	-hold_jid $jids summary.sh
}
testsummary
exit


simple () {
    samples=100
    features=20
    r0=0.5
    instance="synth_${samples}_${features}_$r0"
    qsub -v inst=$instance -e "./$instance.yuu.stats.e" -o "./$instance.yuu.stats.o" run.sh
}

profsimple () {
    samples=200
    features=50
    r0=0.5
    instance="synth_${samples}_${features}_$r0"
    qsub -v inst=$instance -e "./$instance.yuu.stats.e" -o "./$instance.yuu.stats.o" profile.sh
}

#profsimple
#exit

varingsamplesexperiment () {
    vsamples="100 1000 10000 100000"
    features=20
    r0=0.1
    for samples in $vsamples
    do
	instance="synth_${samples}_${features}_$r0"
	qsub -v inst=$instance -e "./$instance.yuu.stats.e" -o "./$instance.yuu.stats.o"  run-mahito.sh
    done
}

varingfeaturesexperiment () {
    samples=200
    vfeatures="10 50 100 200"
    rs="0.1 0.5"
    thres="0.95"
#    thres="0.95 0.9"
    for features in $vfeatures
    do
	for r in $rs
	do
	    for thre in $thres
	    do
		instance="synth_${samples}_${features}_$r0"
#		qsub -v inst=$instance,thre=$thre -e "./$instance.yuu-$thre.stats.e" -o "./$instance.yuu-$thre.stats.o" run.sh
		qsub -v inst=$instance -e "./$instance.mahito.stats.e" -o "./$instance.mahito.stats.o"  run-mahito.sh
	    done
	done
    done
}

#simple
#simple
#varingsamplesexperiment
varingfeaturesexperiment
