#!/bin/bash

bedtools makewindows -g hg19-main.chrom.sizes -w 1000000 > hg19-1mb.bed
bedtools intersect -c -a hg19-kc.bed -b hg19-1mb.bed