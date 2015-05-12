#!/bin/bash

theta=$1
rho=$2
numSites=$3
P=$4
filename=$5

ms 2 1 -t $theta -r $rho $numSites | sed -e '/positions/!d' -e 's/positions: //' | tr ' ' '\n' | sed '/^$/d' > $filename
R --slave -e "cat(sort(runif(rpois(1, ${theta}*${P}))), sep = \"\n\")" >> $filename
