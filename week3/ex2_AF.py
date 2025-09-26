#!/usr/bin/env python3

for line in open("/Users/cmdb/qb25-answers/week3/biallelic.vcf"):
    if line.startswith('#'):
        continue
    fields = line.rstrip('\n').split('\t')
    info_items = fields[7].split(';')
    for item in info_items:
        if item.startswith('AF='):
            af = item.split('=', 1)[1].split(',')[0]  # first AF value
            print(af)
            break