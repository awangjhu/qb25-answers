#!/usr/bin/env python3

import sys
from fasta import readFASTA

# 1.1
def main():
    if len(sys.argv) != 5:
        print("<fasta_file> <scoring_matrix> <gap_penalty> <output_file>")
        sys.exit(1)
    fasta_file = sys.argv[1]
    scoring_matrix_file = sys.argv[2]
    gap_penalty = float(sys.argv[3])
    output_file = sys.argv[4]
    # Read sequences
    input_sequences = readFASTA(open(fasta_file))
    seq1_id, sequence1 = input_sequences[0]
    seq2_id, sequence2 = input_sequences[1]
    # Read scoring matrix
    scoring_matrix = read_scoring_matrix(scoring_matrix_file)
    # Alignment
    aligned_seq1, aligned_seq2, alignment_score = needleman_wunsch(
        sequence1, sequence2, scoring_matrix, gap_penalty)
    # Calculate statistics
    gaps_seq1, gaps_seq2, identity_seq1, identity_seq2 = calculate_statistics(
        aligned_seq1, aligned_seq2)
    # Generate identity alignment 
    # (| for matches, space for mismatches)
    identity_alignment = ''
    for i in range(len(aligned_seq1)):
        if aligned_seq1[i] == aligned_seq2[i]:
            identity_alignment += '|'
        else:
            identity_alignment += ' '
    # 1.5
    # Output
    with open(output_file, 'w') as f:
        f.write(f"Sequence 1: {seq1_id}\n")
        f.write(f"Sequence 2: {seq2_id}\n")
        f.write(f"Alignment Score: {alignment_score}\n")
        f.write(f"Gaps in Sequence 1: {gaps_seq1}\n")
        f.write(f"Gaps in Sequence 2: {gaps_seq2}\n")
        f.write(f"Sequence 1 Identity: {identity_seq1:.2f}%\n")
        f.write(f"Sequence 2 Identity: {identity_seq2:.2f}%\n\n")
        # Write alignment 
        for i in range(0, len(identity_alignment), 100):
            f.write(aligned_seq1[i:i+100] + '\n')
            f.write(identity_alignment[i:i+100] + '\n')
            f.write(aligned_seq2[i:i+100] + '\n\n\n')
    # Print statistics 
    print(f"Gaps in Sequence 1: {gaps_seq1}")
    print(f"Gaps in Sequence 2: {gaps_seq2}")
    print(f"Sequence 1 Identity: {identity_seq1:.2f}%")
    print(f"Sequence 2 Identity: {identity_seq2:.2f}%")
    print(f"Alignment Score: {alignment_score}")
# Chatgpt recommend me format the scoring matrix function like this
# This parses a scoring matrix and creates a dictionary where the keys are tuples of
#character pairs and values are scores
def read_scoring_matrix(filepath):
    scoring_matrix = {}
    with open(filepath) as f:
        lines = f.readlines()
        # first line has characters
        chars = lines[0].strip().split()
        # the rest uses scores
        for i, line in enumerate(lines[1:]):
            values = line.strip().split()
            char1 = values[0]
            for j, score in enumerate(values[1:]):
                char2 = chars[j]
                scoring_matrix[(char1, char2)] = float(score)
    return scoring_matrix
def needleman_wunsch(seq1, seq2, scoring_matrix, gap_penalty):
# Needleman-Wunsch alignment
    m = len(seq1)
    n = len(seq2)
    # 1.2
    # I got help using Chatgpt to create the comprehension syntax. The matrix
    # with (m + 1) rows and (n +1 ) columns are initialized to 0
    F = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    traceback = [['' for _ in range(n + 1)] for _ in range(m + 1)]
    for i in range(m + 1):
        F[i][0] = i * gap_penalty
        if i > 0:
            traceback[i][0] = 'up'  # gap in seq2
    for j in range(n + 1):
        F[0][j] = j * gap_penalty
        if j > 0:
            traceback[0][j] = 'left'  # gap in seq1
    traceback[0][0] = 'done'
    # 1.3
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            char1 = seq1[i - 1]
            char2 = seq2[j - 1]
            # Calculation
            match_score = F[i - 1][j - 1] + scoring_matrix.get((char1, char2), 0)
            gap_seq1 = F[i][j - 1] + gap_penalty  # gap in seq1
            gap_seq2 = F[i - 1][j] + gap_penalty  # gap in seq2 
            if match_score >= gap_seq1 and match_score >= gap_seq2:
                F[i][j] = match_score
                traceback[i][j] = 'diag'
            elif gap_seq1 >= gap_seq2:
                F[i][j] = gap_seq1
                traceback[i][j] = 'left'
            else:
                F[i][j] = gap_seq2
                traceback[i][j] = 'up'
    # 1.4 
    # Traceback is used to find most optimal alignment
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = m, n
    while i > 0 or j > 0:
        direction = traceback[i][j]
        if direction == 'diag':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif direction == 'left':
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1
        elif direction == 'up':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else: 
            break
    # Reversed sequences
    # .join(reversed()) used for reversing and joining all in one step. This allows
    # me to reverse the list and convert it to string - recommended by Chatgpt
    aligned_seq1 = ''.join(reversed(aligned_seq1))
    aligned_seq2 = ''.join(reversed(aligned_seq2))
    alignment_score = F[m][n]
    return aligned_seq1, aligned_seq2, alignment_score
def calculate_statistics(aligned_seq1, aligned_seq2):
    gaps_seq1 = aligned_seq1.count('-')
    gaps_seq2 = aligned_seq2.count('-')
    # Calculate matches
    # This counts positions for sequences that have the same character
    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
    # Percent identity = (matches / length without gaps) * 100
    len_seq1_no_gaps = len(aligned_seq1) - gaps_seq1
    len_seq2_no_gaps = len(aligned_seq2) - gaps_seq2
    identity_seq1 = (matches / len_seq1_no_gaps * 100) if len_seq1_no_gaps > 0 else 0
    identity_seq2 = (matches / len_seq2_no_gaps * 100) if len_seq2_no_gaps > 0 else 0 
    return gaps_seq1, gaps_seq2, identity_seq1, identity_seq2
if __name__ == "__main__":
    main()