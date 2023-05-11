import os
import msprime
import demes
import sys

if len(sys.argv) != 3:
    print("Usage: python simulate.py <demes_file_path> <output_file_path>")
    sys.exit(1)

demes_file_path = sys.argv[1]
if not os.path.exists(demes_file_path):
    print(f"File {demes_file_path} not found.")
    sys.exit(1)
output_file_path = sys.argv[2]

mutation_rate = 1e-8 #set mutation rate
sq_len = 1e5 #length of chromosome (change to 250e6)
n_individuals = 100*2 #remember diploid
rate_recomb = 1e-8 #average recombination rate in humans
static_seed = 5555 #seed 5555 to eliminate randomness in simulation


graph = demes.load(demes_file_path)
demography = msprime.Demography.from_demes(graph)
ts = msprime.sim_ancestry({"YRI": n_individuals,"CEU": n_individuals,"CHB": n_individuals,"KAR": n_individuals,"ADMIX": n_individuals},
 sequence_length = sq_len, demography = demography, random_seed = static_seed, recombination_rate = rate_recomb)
mts = msprime.sim_mutations(ts, rate = mutation_rate, random_seed = static_seed)

mts.write_vcf(sys.stdout)

#saving the the ancestery of every mutation as a tsv file
with open(output_file_path, "w") as f:
    # Header row to the file
    f.write("Position\tAncestral_population\n")
    
    # Iterate over the variants and write the data to the file
    for v in mts.variants():
        site = v.site
        pos = site.position
        mut = site.mutations
        for mutation in mut:
            node_num = mutation.node
            pop_id = mts.population(mts.nodes_population[node_num]).metadata['name']
        f.write(f"{pos}\t{pop_id}\n")  # Write the data to the file

