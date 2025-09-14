#!/usr/bin/env bash

#Construct a bedtools command to test where there is any overlap between 1_Active and 12_Repressed in a given condition (aka mutually exclusive)

#Test is there is any overlap between 1_Active and 12_Repressed in a given condition
bedtools intersect -u -a nhek-active.bed -b nhek-repressed.bed > nhek_active_ and_repressed.bed
wc -l nhek_active_and_repressed.bed  #Output: 0 nhek_active_ and_repressed.bed

#Bedtool command to find regions active in NHEK and NHLF, and one to find regions that are active in NHEK but not active in NHLF
#Active in both
bedtools intersect -u -a nhek-active.bed -b nhlf-active.bed > active_in_both.bed
wc -l active_in_both.bed #Output: 11608 active_in_both.bed
#Active in NHEK but not active in NHLF
bedtools intersect -v -a nhek-active.bed -b nhlf-active.bed > active_only_nhek.bed
wc -l active_only_nhek.bed #Output: 2405 active_only_nhek.bed

# How many features are output by the first command? by the second command?
#Answer: 11608 output in active_in_both.bed, 2405 in active_only_nhek.bed

#Do these two numbers add up to the original number of lines in nhek-active.bed?
#Yes, these numbers add up to the original 14013 in the original bed file

#Bedtool intersect commands to see the effect of using the arguments -f 1, -F 1, and -f 1 -F 1 when comparing -a nhek-active.bed -b nhlf-active.bed
bedtools intersect -u -f 1 -a nhek-active.bed -b nhlf-active.bed > f1.bed
wc -l f1.bed #Output: 4821 f1.bed
head -n1 f1.bed #Output: chr1	25558413	25559413	1_Active_Promoter	0	.	25558413	25559413

bedtools intersect -u -F 1 -a nhek-active.bed -b nhlf-active.bed > F1.bed
wc -l F1.bed #Output: 6401 F1.bed 
head -n1 F1.bed #Output: chr1	19922613	19924613	1_Active_Promoter	0	.	19922613	19924613

bedtools intersect -u -f 1 -F 1 -a nhek-active.bed -b nhlf-active.bed > f1F1.bed
wc -l f1F1.bed #Output: 1409 f1F1.bed
head -n1 f1F1.bed #Output: chr1	1051137	1051537	1_Active_Promoter	0	.	1051137	1051537

#How does the relationship between the NHEK and NHLF chromatin state change as you alter the overlap parameter?
#As I alter the overlap parameter to -f 1, the loci where the NHEK active promoter seems contained within NHLF active promoter. When looking at the active promoter of NHEK, it sits inside of NHLF but the boundaries do not extend pas NHLF active promoter. This means there are fewer overlaps and show a shared core promoter segment.
#When altering the overlap parameter to -F 1, the loci where NHLF active promoter is seems to be contained inside of NHEK active promoter.
#Finally, when altering the overlap to -f 1 -F 1, there seem to be identical segments retained with the same boundaries.

#Construct three bedtools intersect commands to identify the following types of regions.
#Describe the chromatin state across all nine conditions.
# Active in NHEK, Active in NHLF
bedtools intersect -u -a nhek-active.bed -b nhlf-active.bed > active_active.bed 
head -n1 active_active.bed #Output: chr1	19922613	19924613	1_Active_Promoter	0	.	19922613	19924613
#Description: Mostly active

# Active in NHEK, Repressed in NHLF.
bedtools intersect -u -a nhek-active.bed -b nhlf-repressed.bed > active_repressed.bed
head -n1 active_repressed.bed #Output: chr1	1981140	1981540	1_Active_Promoter	0	.	1981140	1981540
#Description: Mostly active, a few repressed, with HUVEC having a poised promoter.

# Repressed in NHEK, Repressed in NHLF
bedtools intersect -u -a nhek-repressed.bed -b nhlf-repressed.bed > repressed_repressed.bed
head -n1 repressed_repressed.bed #Output: chr1	11534013	11538613	12_Repressed	0	.	115340111538613
#Description: Very few active, mostly repressed and poised.


