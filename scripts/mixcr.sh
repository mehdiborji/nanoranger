#!/bin/bash

echo 'filename = ' $2
echo 'save_prefix = ' $1
echo 'species = ' $3
echo 'nthreads = ' $4

mixcr align -r $1_report.txt -s $3 -t $4 -OsaveOriginalReads=true  $2 $1.vdjca -f

mixcr assemble -t $4 -r $1_report.txt $1.vdjca $1.clns -f \
-OclusteringFilter.specificMutationProbability=5E-2 \
-OaddReadsCountOnClustering=true

mixcr exportClones -p fullImputed  -topChains -chains $1.clns $1_clones.txt -f

mixcr assemble -t $4 -r $1_report.txt -a $1.vdjca $1.clna -f \
-OclusteringFilter.specificMutationProbability=5E-2 \
-OaddReadsCountOnClustering=true

mixcr exportAlignments -descrsR1 -cloneIdWithMappingType -cloneId -topChains -chains $1.clna $1_cloneID.txt -f

pigz -f $1_cloneID.txt
pigz -f $1_clones.txt
#rm $1.vdjca $1.clna $1.clns
