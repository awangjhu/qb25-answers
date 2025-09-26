#!/usr/bin/env python3

for line in open("/Users/cmdb/qb25-answers/week3/biallelic.vcf"):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t")
    fmt = f[8].split(":")
    if "DP" not in fmt: continue
    dp_idx = fmt.index("DP")
    for s in f[9:]:
        vals = s.split(":")
        if dp_idx < len(vals) and vals[dp_idx] not in (".",""):
            print(vals[dp_idx])