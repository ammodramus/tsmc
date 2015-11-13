# without asex

numReps=5


initialT="0.5"


for rep in $(seq 1 $numReps)
do
    unstructuredFilename="constantN_no_asex_rep${rep}.psmc"
    unstructuredFilenameResults="results_constantN_no_asex_rep${rep}.txt"
    unstructuredCommand="./psmc ${unstructuredFilename} > $unstructuredFilenameResults"
    echo $unstructuredCommand
    unstructuredFilenameResultsWithT="results_constantN_no_asex_with_divergence_rep${rep}.txt"
    unstructuredCommandWithT="./psmc ${unstructuredFilename} -T $initialT > $unstructuredFilenameResultsWithT"
    echo $unstructuredCommandWithT
done

for asexT in 0.1 0.25 0.5
do
    for rep in $(seq 1 $numReps)
    do
        asexFilename="constantN_asex${asexT}_rep${rep}.psmc"
        asexResults="results_constantN_asex${asexT}_rep${rep}"
        asexResultsWithT="results_constantN_asex${asexT}_with_divergence_rep${rep}"
        command="./psmc ${asexFilename} > $asexResults"
        echo $command
        commandWithT="./psmc ${asexFilename} -T $initialT > $asexResultsWithT"
        echo $commandWithT
    done
done
