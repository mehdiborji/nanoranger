#!/bin/bash

module load java/jdk-1.8u112

echo 'filename = ' $2
echo 'save_prefix = ' $1
echo 'species = ' $3
echo 'nthreads = ' $4

#mixcr align -r $1_report_asmbl.txt -s $3 -t $4 $2 $1_asmbl.vdjca -f

mixcr assemble -t $4 -r $1_report_asmbl.txt -a $1_asmbl.vdjca $1_asmbl.clna -f \
-OclusteringFilter.specificMutationProbability=5E-2 \
-OaddReadsCountOnClustering=true

mixcr assembleContigs -t $4 -r $1_report_asmbl.txt $1_asmbl.clna $1_asmbl.clns

mixcr exportClones -p fullImputed $1_asmbl.clns $1_clones_asmbl.txt -f

pigz -f $1_clones_asmbl.txt
#rm $1.vdjca $1.clna $1.clns
