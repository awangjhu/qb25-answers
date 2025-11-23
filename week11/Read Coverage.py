#!/usr/bin/env python3

import numpy as np
from scipy.stats import poisson, norm

def simulate_coverage(genome_size, read_length, coverage):
    num_reads = int((genome_size * coverage) / read_length)   
    # I was suggested to use dtype=int to save memory according to Chatgpt 
    genome_coverage = np.zeros(genome_size, dtype=int)
    # simulate reads (from the psuedocode)
    for i in range(num_reads):
        # random start position
        start_pos = np.random.randint(0, genome_size - read_length + 1)
        end_pos = start_pos + read_length
        # Chatgpt helped with NumPy slice syntax for incrementing ranges
        genome_coverage[start_pos:end_pos] += 1
    return genome_coverage, num_reads

def calculate_distribution_estimates(max_coverage, coverage_lambda):
    xs = list(range(0, max_coverage + 1))
    poisson_estimates = poisson.pmf(xs, mu=coverage_lambda)
    mean = coverage_lambda
    stddev = np.sqrt(coverage_lambda)
    normal_estimates = norm.pdf(xs, loc=mean, scale=stddev)
    
    return xs, poisson_estimates, normal_estimates

def save_coverage_data(genome_coverage, coverage_level, num_reads):
    filename = f'coverage_{coverage_level}x.txt'
    np.savetxt(filename, genome_coverage, fmt='%d')
    
    # Calculate statistics
    max_coverage = np.max(genome_coverage)
    zero_coverage = np.sum(genome_coverage == 0)
    zero_percent = (zero_coverage / len(genome_coverage)) * 100
    xs, poisson_est, normal_est = calculate_distribution_estimates(
        max_coverage, coverage_level)
    #Chatgpt helped with file writing pattern    
    dist_filename = f'distributions_{coverage_level}x.txt'
    with open(dist_filename, 'w') as f:
        f.write('coverage\tpoisson\tnormal\n')
        for x, p, n in zip(xs, poisson_est, normal_est):
            f.write(f'{x}\t{p}\t{n}\n')
    
    # Print statistics
    print(f"\n{coverage_level}x Coverage Results:")
    print(f"Number of reads: {num_reads}")
    print(f"Max coverage: {max_coverage}")
    print(f"Positions with 0x coverage: {zero_coverage} ({zero_percent:.4f}%)")
    print(f"Expected 0x coverage (Poisson): {poisson.pmf(0, coverage_level) * 100:.4f}%")
    
    return zero_percent

if __name__ == "__main__":
    genome_size = 1_000_000
    read_length = 100
    
    # Simulate 3x coverage
    print("Simulating 3x coverage...")
    coverage_3x, num_reads_3x = simulate_coverage(genome_size, read_length, 3)
    save_coverage_data(coverage_3x, 3, num_reads_3x)
    
    # Simulate 10x coverage
    print("\nSimulating 10x coverage...")
    coverage_10x, num_reads_10x = simulate_coverage(genome_size, read_length, 10)
    save_coverage_data(coverage_10x, 10, num_reads_10x)
    
    # Simulate 30x coverage
    print("\nSimulating 30x coverage...")
    coverage_30x, num_reads_30x = simulate_coverage(genome_size, read_length, 30)
    save_coverage_data(coverage_30x, 30, num_reads_30x)