#!/usr/bin/env python3

import sys

if len(sys.argv) != 2:
    sys.exit(1)

sam_file = open(sys.argv[1])

my_dict = {} #this is for chromsome count
nm_dict = {} #for mismatch count 

for line in sam_file:
    if line.startswith("@"):
        continue #this will skip for SAM headers
    fields = line.rstrip("\n").split("\t") 
    if len(fields) < 11:
        continue
    
    rname = fields[2]
    my_dict[rname] = my_dict.get(rname, 0) + 1

    for f in fields[11:]: #count by mismatch
        if f.startswith("NM:i:"):
            nm = int(f[5:])
            nm_dict[nm] = nm_dict.get(nm, 0) + 1
            break

for chrom in my_dict.keys():
    print(f"{chrom}\t{my_dict[chrom]}")

print("---")

for nm in sorted(nm_dict.keys()):
    print(f"{nm}\t{nm_dict[nm]}") #nm count