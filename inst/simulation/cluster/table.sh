#! /bin/bash

work=`pwd`

for P1 in 400 200; do
    for P2 in M6 M2 M3 M4 M5; do
	for P3 in TRUE FALSE; do
	    cd $work/dat-${P1}-${P2}-${P3}
	    echo `pwd`
	    cat output-${P1}-${P2}-${P3}-* > gany-${P1}-${P2}-${P3}
	    rm *.Rout
	    cd $work
	    cp $work/dat-${P1}-${P2}-${P3}/gany-* $work/output/
	done
    done
done
