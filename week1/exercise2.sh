#!/usr/bin/env bash

#Bash script with commands to calculate unique genes
#To find the total in each assembly:
bedtools sort -i hg19-kc.bed > hg19-kc.sorted.bed
wc -l hg19-kc.sorted.bed #Output: 80269 hg19-kc.sorted.bed
bedtools sort -i hg16-kc.bed > hg16-kc.sorted.bed
wc -l hg16-kc.sorted.bed #Output: 21364 hg16-kc.sorted.bed
#In hg19 but not in hg16
bedtools intersect -v -a hg19-kc.sorted.bed -b hg16-kc.sorted.bed > hg19_unique_vs_hg16.bed
wc -l hg19_unique_vs_hg16.bed #Output:  42717 hg19_unique_vs_hg16.bed
#In hg16 but not in hg19
bedtools intersect -v -a hg16-kc.sorted.bed -b hg19-kc.sorted.bed > hg16_unique_vs_hg19.bed
wc -l hg16_unique_vs_hg19.bed #Output: 3460 hg16_unique_vs_hg19.bed

#Questions:
#How many genes are in hg19?
#80269 genes are in hg19
#How many genes are in hg19 but not in hg16?
#42717 genes are in hg19 but not in hg16
#Why are some genes in hg19 but not in hg16?
#There are significantly more genes in hg19 than 16 because it was the newly updated and, thus, has more annotations.

#How many genes are in hg16?
#21364 genes are in hg16
#How many genes are in hg16 but not in hg19?
#3460 genes are in hg16 but not in hg19
#Why are some genes in hg16 but not in hg19?
#There are some 3460 genes in hg16 but not in hg19 likely because they use different systems to annotate and compare sequences, so some genes may not line up correctly.