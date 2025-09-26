#!/usr/bin/env python3

#3.1. based on the IGV, I think A01_09, A01_11, A01_23, A01_35, A01_39 are BY lab strains while the rest are likely RM wine samples since these have significantly more mismatches in the chromosome region.
#3.2
import sys

vcf_file = sys.argv[1] if len(sys.argv) > 1 else "biallelic.vcf"

sample_ids = ["A01_62", "A01_39", "A01_63", "A01_35", "A01_31",
              "A01_27", "A01_24", "A01_23", "A01_11", "A01_09"]

with open(vcf_file) as fh, open("gt_long.txt", "w") as out:
    out.write("sample\tchrom\tpos\tgenotype\n")
    for line in fh:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        chrom = fields[0]
        pos   = fields[1]
        sample_cols = fields[9:]
        for s_id, s_col in zip(sample_ids, sample_cols):
            gt = s_col.split(":")[0]
            if gt == "0":
                out.write(f"{s_id}\t{chrom}\t{pos}\t0\n")
            elif gt == "1":
                out.write(f"{s_id}\t{chrom}\t{pos}\t1\n")
