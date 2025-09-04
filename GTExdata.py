#!/usr/bin/env python3

import sys

my_dict = {}
my_file = open("/Users/cmdb/qb25-answers/unix-python-scripts/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct")

_ = my_file.readline()
_ = my_file.readline()
header = my_file.readline().rstrip().split("\t")
data = my_file.readline().rstrip().split("\t")
#split and list field variables
for line in range(len(data)):
    my_dict[header[line]] = data[line]
my_file.close()
#create a dictionary and loop header and data as a value

GTEx_file_2 = open("/Users/cmdb/qb25-answers/unix-python-scripts/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
for line_2 in GTEx_file_2:
    line_2 = line_2.strip('\n').split("\t")
    SAMPID = line_2[0]
    if SAMPID in my_dict:
        print(f"SAMPID: {SAMPID}")
        print(f"Expression: {my_dict.get(SAMPID)}")
        print(f"SMTSD: {line_2[6]}")

print("Brain, Adrenal, Thyroid")
GTEx_file_2.close() 