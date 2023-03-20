import sys
import argparse
import os
from os import path
import gzip
import msprime, pyslim
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="Input-file")

args = parser.parse_args()

input=args.input

# Load the .trees file
orig_ts = pyslim.load(input) # no simplify!

alive = orig_ts.individuals_alive_at(0)

print(f"There are {len(alive)} individuals alive from the final generation.")

num_alive = [0 for _ in range(orig_ts.num_populations)]
for i in alive:
   ind = orig_ts.individual(i)
   num_alive[ind.population] += 1

for pop, num in enumerate(num_alive):
   print(f"Number of individuals in population {pop}: {num}")


pop_configs = [msprime.PopulationConfiguration(initial_size=7310)
               for _ in range(orig_ts.num_populations)]
rts = orig_ts.recapitate(population_configurations=pop_configs,
                         migration_matrix=[[0.0, 0.0, 0.0, 0.0],
                                           [0.0, 0.0, 0.000025, 0.0000078],
                                           [0.0, 0.000025, 0.0, 0.0000311],
                                           [0.0, 0.0000078, 0.0000311, 0.0]],
                         recombination_rate=1e-8,
                         random_seed=4)

mutated = pyslim.SlimTreeSequence(msprime.mutate(rts, rate=2.38e-8, keep=True))


pop_indivs = [[], [], [], []]
pop_nodes = [[], [], [], []]
for i in mutated.individuals_alive_at(0):
   ind = mutated.individual(i)
   pop_indivs[ind.population].append(i)
   pop_nodes[ind.population].extend(ind.nodes)

p2 = mutated.simplify(pop_nodes[2])

alive_inds = p2.individuals_alive_at(0)
keep_indivs = np.random.choice(alive_inds, 1000, replace=False)
keep_nodes = []
for i in keep_indivs:
   keep_nodes.extend(p2.individual(i).nodes)

simplified = p2.simplify(keep_nodes)


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
	
with open(input + ".vcf", "wt") as vcf_file:
	simplified.write_vcf(output=vcf_file,
	individuals=sampled_nodes,
	individual_names=sample_names,
	position_transform="legacy",
	contig_id='1')



