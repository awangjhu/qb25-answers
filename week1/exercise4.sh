#!/usr/bin/env bash

#Confirming that the snps-chr1.bed has 1091148 lines
wc snps-chr1.bed #1091148 6546888 42505620 chr1snp.bed

#Determining which gene has the most SNPS
bedtools sort -i snps-chr1.bed > snps-chr1.sorted.bed
bedtools sort -i hg19-kc.bed   > hg19-kc.sorted.bed
bedtools intersect -c -a hg19-kc.sorted.bed -b snps-chr1.sorted.bed > hg19-kc_snpcounts_chr1.bed
sort -k5,5nr hg19-kc_snpcounts_chr1.bed | head -n1 | tee top_gene_chr1.tsv 
#Output: chr1	245912648	246670581	ENST00000490107.6_7	5445
#Describing the gene:
#Synthetic name: ENST00000490107.6_7
#Human readout: SMYD3
#Position: g19 chr1:245,912,649-246,670,581
#Size: 757,933
#Exon count: 12

#Determining which SNPs lie within vs outside of a gene
bedtools sample -i snps-chr1.sorted.bed -n 20 -seed 42 > snps-chr1-20.bed
bedtools sort -i snps-chr1-20.bed > snps-chr1-20.sorted.bed
bedtools closest -d -t first -a snps-chr1-20.sorted.bed -b hg19-kc.sorted.bed > snp20_closest.tsv

#How many SNPS are inside of a gene
awk '$NF==0' snp20_closest.tsv | wc -l #Output: 15
#15/20 SNPs are inside of the gene

#Distance for SNPs outside genes
awk '$NF>0{ if(!n++){min=$NF;max=$NF} if($NF<min)min=$NF; if($NF>max)max=$NF }
     END{ if(!n) print "All inside"; else printf "%d–%d\n", min, max }' snp20_closest.tsv
#The range of distances for the ones outside of a gene are 1664–22944