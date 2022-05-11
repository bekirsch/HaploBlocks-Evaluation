import sys
import argparse
import os
from os import path
import msprime, pyslim
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="Input-file")

args = parser.parse_args()

input=args.input

# Load the SLiM -.trees file
ts = pyslim.load(input)

# Recapitate, see https://pyslim.readthedocs.io/en/latest/tutorial.html# for further info
recap = ts.recapitate(recombination_rate=1e-8, Ne=10000)

# Simulate neutral mutaions, see https://pyslim.readthedocs.io/en/latest/tutorial.html# for further info; Save the resulting treesequence for future reference
mutated = pyslim.SlimTreeSequence(msprime.mutate(recap, rate=2.5e-8, keep=True))
mutated.dump(input + ".recap.mutate")

# Find individuals alive at present and sample 1000 randomly
alive_inds = mutated.individuals_alive_at(0)
keep_indivs = np.random.choice(alive_inds, 1000, replace=False)
keep_nodes = []
for i in keep_indivs:
   keep_nodes.extend(mutated.individual(i).nodes)

# Simplify the treesequence for the sampled individuals
simplified = mutated.simplify(keep_nodes)

# Get all the nodes (chromosomes) alive at present (aka, time 0)
nodes = simplified.individuals_alive_at(0)

# Determine all samples that are in a population
samples_of_pop = dict()
for u in nodes:
    pop = simplified.individual(u).population
    samples_of_pop.setdefault(pop, []).append(u)

# Create new sample IDs
sample_names = []
sampled_nodes = []
for sam in samples_of_pop[pop]:
	sampled_nodes.append(sam)
	new_sam_name = 'Ind_{}'.format(sam)
	sample_names.append(new_sam_name)

# Save as VCF; position_transform="legacy" prevents multiple SNPs from occoring at one site due to rounding
with open(input + ".vcf", "wt") as vcf_file:
	simplified.write_vcf(output=vcf_file,
        individuals=sampled_nodes,
        individual_names=sample_names,
        position_transform="legacy",
        contig_id='1')
