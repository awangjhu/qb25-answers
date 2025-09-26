#!/bin/bash

tar -xzvf BYxRM_bam.tar.gz #extract the tar file
#Index each bam file
samtools index A01_09.bam
samtools index A01_11.bam 
samtools index A01_23.bam 
samtools index A01_24.bam
samtools index A01_27.bam 
samtools index A01_31.bam 
samtools index A01_35.bam
samtools index A01_39.bam 
samtools index A01_62.bam
samtools index A01_63.bam
#Read-count loop
for f in *.bam; do 
n=$(samtools view -c -F 260 "$f") 
printf "%s\t%s\n" "$f" "$n" >> read_counts.txt 
done
#Output
head read_counts.txt:
file	primary_mapped_reads
A01_09.bam	669520
A01_11.bam	656196
A01_23.bam	708697
A01_24.bam	797335
A01_27.bam	602368
A01_31.bam	610325
A01_35.bam	803495
A01_39.bam	713671
A01_62.bam	816598

# run FreeBayes to discover variants 
freebayes -f sacCer3.fa -L bamListFile.txt --genotype-qualities -p 1> unfiltered.vcf

# the resulting VCF file is unfiltered, meaning that it contains low-confidence calls and also has some quirky formatting, so the following steps use a software suite called vcflib to clean up the VCF

# filter the variants based on their quality score and remove sites where any sample had missing data
vcffilter -f "QUAL > 20" -f "AN > 9" unfiltered.vcf > filtered.vcf

# FreeBayes has a quirk where it sometimes records haplotypes rather than individual variants; we want to override this behavior
vcfallelicprimitives -kg filtered.vcf > decomposed.vcf

# in very rare cases, a single site may have more than two alleles detected in your sample; while these cases may be interesting, they may also reflect technical errors and also pose a challenge for parsing the data, so we remove them
vcfbreakmulti decomposed.vcf > biallelic.vcf