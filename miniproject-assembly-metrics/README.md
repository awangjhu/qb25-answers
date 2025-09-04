# Caenorhabditis remanei
## Four genomic assemblies
### Python Script
#!/usr/bin/env python3

import sys
import fasta

my_file = open(sys.argv[1])
contigs = fasta.FASTAReader(my_file)
counter = 0
length = 0 
Contig_list = []

for ident, sequence in contigs:
    counter += 1
    length += len(sequence)
    Contig_list.append(len(sequence))
Contig_list.sort(reverse = True)
Average_length = length / counter

C_length = 0
for i in Contig_list:
    C_length += i
    if C_length >= 0.5*(length):
        print(i)
        break

print(f"Total number of contigs: {counter}")
print(f"Total length: {length}")
print (f"Average length: {Average_length} ")
print (f"N50 : {i} ")

### Instructions
1. cd ~/qb25-answers/miniproject-assembly-metrics
2. wget https://raw.githubusercontent.com/bxlab/cmdb-quantbio/refs/heads/main/resources/code/fasta.py
3. open fasta.py
4. wget "link fa.gz from database" 
i.e. wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/caenorhabditis_remanei/PRJNA248909/caenorhabditis_remanei.PRJNA248909.WBPS19.genomic.fa.gz 
5. gunzip "file name"
i.e. gunzip caenorhabditis_remanei.PRJNA248909.WBPS19.genomic.fa.gz 
6. ./assembly-metrics.py "link fa.gz from database"
i.e. ./assembly-metrics.py caenorhabditis_remanei/PRJNA248909/caenorhabditis_remanei.PRJNA248909.WBPS19.genomic.fa.gz 

## Contig Results
### Contig 1
["https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/caenorhabditis_remanei/PRJNA248909/caenorhabditis_remanei.PRJNA248909.WBPS19.genomic.fa.gz]

Results:
Total number of contigs: 1591
Total length: 118549266
Average length: 74512.42363293526 
N50 : 1522088

### Contig 2
[https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/caenorhabditis_remanei/PRJNA248911/caenorhabditis_remanei.PRJNA248911.WBPS19.genomic.fa.gz]

Results: 
Total number of contigs: 912
Total length: 124541912
Average length: 136559.11403508772 
N50 : 1765890 

### Contig 3
[https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/caenorhabditis_remanei/PRJNA53967/caenorhabditis_remanei.PRJNA53967.WBPS19.genomic.fa.gz]

Results:
Total number of contigs: 3670
Total length: 145442736
Average length: 39630.17329700272 
N50 : 435512

### Contig 4
[https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/caenorhabditis_remanei/PRJNA577507/caenorhabditis_remanei.PRJNA577507.WBPS19.genomic.fa.gz]

Results
Total number of contigs: 187
Total length: 130480874
Average length: 697758.6844919786 
N50 : 21501900