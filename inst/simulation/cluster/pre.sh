#! /bin/bash

work=`pwd`
total=1000
node=5
core=10
nn=5

let sp=$(($total/$nn))

for P1 in 200 400; do
    for P2 in M2 M5; do
	for P3 in TRUE FALSE; do
	    for P4 in TRUE; do
		if [ ! -d dat-${P1}-${P2}-${P3}-${P4} ]; then
		    mkdir dat-${P1}-${P2}-${P3}-${P4}
		    cp submit.sh dat-${P1}-${P2}-${P3}-${P4}
		fi
		cd $work/dat-${P1}-${P2}-${P3}-${P4}
		ln -s ../setup.R .
		sed "s/NODE/$node/; s/CORE/$core/; s/SIMJOB/job-${P1}-${P2}-${P3}-${P4}/;" ../submit.sh > submit.sh
		sed "s/P1/$P1/; s/P2/$P2/; s/P3/$P3/; s/P4/$P4/; s/nn/$nn/;" ../base.R > run-${P1}-${P2}-${P3}-${P4}.R
		for i in $( seq 1 $sp);
		do
		    let i=i-1
		    echo "Rscript --slave --no-restore --no-save run-${P1}-${P2}-${P3}-${P4}.R ${i} > out${i}.Rout" >> job-${P1}-${P2}-${P3}-${P4}
		done
		sbatch submit.sh
		sleep 1
		cd ..
	    done
	done
    done
done
