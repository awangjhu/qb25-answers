#!/usr/bin/env python3

def generate_debruijn_edges(reads, k):
    #suggested by Chatgpt to automatically duplicate edges
    edges = set()
    # Pseudocode
    for read in reads:
        for i in range(len(read) - k):
            kmer1 = read[i:i+k]
            kmer2 = read[i+1:i+1+k]
            edges.add((kmer1, kmer2))
    
    return edges

# write edge to text file, one edge per line
def write_edges_file(edges, filename='debruijn_edges.txt'):
    with open(filename, 'w') as f:
        for kmer1, kmer2 in sorted(edges):
            f.write(f'{kmer1} -> {kmer2}\n')
    
    print(f"Edges at {filename}")

def write_graphviz_file(edges, filename='debruijn_graph.dot'):
    with open(filename, 'w') as f:
        # Chatgpt suggested: Digraph creates a directed graph
        f.write('digraph DeBruijn {\n')
        f.write('    rankdir=LR;\n')  #left to right
        f.write('    node [shape=circle];\n') #circular nodes
        
        for kmer1, kmer2 in sorted(edges):
            f.write(f'    "{kmer1}" -> "{kmer2}";\n')
        
        f.write('}\n')
    
    print(f"Graphviz file written to {filename}")
    print(f"\nTo generate the graph image, run:")
    print(f"dot -Tpng {filename} -o ex2_digraph.png")

reads = ['ATTCA', 'ATTGA', 'CATTG', 'CTTAT', 'GATTG', 'TATTT', 
         'TCATT', 'TCTTA', 'TGATT', 'TTATT', 'TTCAT', 'TTCTT', 'TTGAT']

# Generate de Bruijn graph with k=3
k = 3
edges = generate_debruijn_edges(reads, k)

# Print stats (Chtagpt helped me with the string formatting)
print(f"Number of reads: {len(reads)}")
print(f"K-mer size: {k}")
print(f"Number of unique edges: {len(edges)}")
print()
write_edges_file(edges)
print()
write_graphviz_file(edges)
print()
print("All edges:")
for kmer1, kmer2 in sorted(edges):
    print(f"{kmer1} -> {kmer2}")