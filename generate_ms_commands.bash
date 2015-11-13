# without asex

numReps=5

unstructuredCommand="scrm 2 100 -t 30000 -r 18000 30000000 | ./ms2psmcfa.py"

for rep in $(seq 1 $numReps)
do
    unstructuredFilename="constantN_no_asex_rep${rep}.psmc"
    echo "$unstructuredCommand > $unstructuredFilename"
done

for asexT in 0.1 0.25 0.5
do
    asexCommand="scrm 2 100 -t 30000 -r 18000 30000000 | ./add_asex.py - ${asexT} | ./ms2psmcfa.py"
    for rep in $(seq 1 $numReps)
    do
        asexFilename="constantN_asex${asexT}_rep${rep}.psmc"
        echo "$asexCommand > $asexFilename"
    done
done
