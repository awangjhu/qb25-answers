#!/bin/bash

#Question 1
#How many features (lines)?
wc ce11_genes.bed
#Output: 53935  323610 2200094 ce11_genes.bed
#Answer: 53935

#How many features per chr? e.g. chrI, chrII
cmdb@QuantBio-13 unix-pythin-scripts % cut -f 1 ce11_genes.bed | uniq -c 
#Answer:
#5460 chrI
#12 chrM
#9057 chrV
#6840 chrX
#6299 chrII
#21418 chrIV
#4849 chrIII

#How many features per strand? e.g. +, -
unix-pythin-scripts % cut -f 6 ce11_genes.bed | sort | uniq -c 
#Answer:
#26626 -
#27309 +

#Question 3
#Which three SMTSDs (Tissue Site Detail) have the most samples?
cut -f 7 GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt | sort | uniq -c |sort -r | head -n 3 
#Answer
#3288 Whole Blood
#1132 Muscle - Skeletal
#867 Lung
#How many lines have “RNA”?
cut -f 12 GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt | grep RNA | wc
#Answer
#20016  117213  947910