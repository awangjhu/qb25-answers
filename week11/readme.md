*Question 1.1: How many 100bp reads are needed to sequence a 1Mbp genome to 3x coverage?
Answer: Number of reads = (10^6 * 3)/100 = 30,000 reads

*Question 1.4: In your simulation, how much of the genome has not been sequenced (has 0x coverage)? How well does this match Poisson expectations? How well does the normal distribution fit the data?
Answer: In the 3x coverage plot, the bar at position 0 show a frequency of 50,000 which means around 5% of the genome has 0x coverage. 
This matches well with the Poisson distribution, there are a positions unsequenced but overall its a good match. The normal distribution fits fine but its not perfect at 3x coverage because the normal distributiona and Poisson do not overlap and are identical.

*Question 1.5: In your simulation, how much of the genome has not been sequenced (has 0x coverage)?
How well does this match Poisson expectations? How well does the normal distribution fit the data?
Answer: The 10x plot at position 0 is small. The frequency is around 40 or 50, which means around 0.005% of the genome has 0x coverage. This matches well with the Poisson expectations since we have few positiosn that are unsequenced. Arguably, the normal distribution fits much better at 10x coverage that 3x coverage since the normal distribution and Poisson lines practically overlap in the distribution.

*Question 1.6: *Question 1.5: In your simulation, how much of the genome has not been sequenced (has 0x coverage)? How well does this match Poisson expectations? How well does the normal distribution fit the data?
Answer: Position 0 in the 30x plot is not visable, which means that practically 0% of the genome has 0 coverage. The Poison prediction matches well because the probability of any position having 0 coverage is 0. The normal distribution perfectly overlaps with the Poisson at 30x coverage which means that the normal distribution fits well with the data.

*Question 2.4:
Answer: dot -Tpng debruijn_graph.dot -o ex2_digraph.png

*Question 2.5:
Answer: Based on my path, one possible sequence that would produce these reads would be: ATTGATTCATTATTTCTTA

*Question 2.6:
Answer: To accurately reconstruct the full genome, we would first need higher converage because higher coverage makes sure that the genome is fully represented. Longer reads would also be helpful for sorting out ambiguities in repetitive sequences and with the path. I would also argue coverage depth information is important to determining the correct path. Additionally, coverage depth information is also quite important because knowing how many times each edge should be traversed can easily help determine the correct path. Perhaps it would alos be helpful to know the starting and endpoints. I would argue this is what it would take to accurately reconstruct the genome.