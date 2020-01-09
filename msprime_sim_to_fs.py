import math
import argparse
import msprime
import os

parser = argparse.ArgumentParser()
parser.add_argument("-out", help="Path and filename to save the tree sequence to e.g. Documents/msprimetest NB directory must already exist", required=True)
parser.add_argument("-path", help="Path to save FINESTRUCTURE input files to e.g. Documents/msprime2fstest", required=True)
args = parser.parse_args()

# initial population sizes:
N_modern = 100000
N_steppe = 20000
N_whg = 10000
N_ehg = 10000
N_neo = 50000
N_chg = 10000
N_A = 5000  # Ancestor of WHG and EHG
N_B = 5000  # Ancestor of CHG and Neolithic farmers

# Time of events
T_modern = 170
T_steppe = 200
T_neo = 250
T_near_east = 800
T_europe = 500
T_basal = 1500

# Migration
hg_mig_rate = 2.5e-5  # Unsure
migration_matrix = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
nhaps = [20, 40, 40, 20, 20]
times = [0, 200, 180, 200, 200]
popnames = ["modern", "neolithic", "steppe", "WHG", "EHG"]
populations = [0, 0, 1, 2, 3]

# Populations: 0=modern/neolithic_farmers/B,1=steppe/CHG,2=WHG/A, 3=EHG
population_configurations = [
    msprime.PopulationConfiguration(initial_size=N_modern),
    msprime.PopulationConfiguration(initial_size=N_steppe),
    msprime.PopulationConfiguration(initial_size=N_whg),
    msprime.PopulationConfiguration(initial_size=N_ehg)
]

Modern_formation = [
    msprime.MassMigration(time=T_modern, source=0, dest=1, proportion=0.5),
    msprime.PopulationParametersChange(time=T_modern, initial_size=N_neo, population=0)
]

Steppe_formation = [
    msprime.MassMigration(time=T_steppe, source=1, dest=3, proportion=0.5),
    msprime.PopulationParametersChange(time=T_steppe, initial_size=N_chg, population=1),
    msprime.MigrationRateChange(time=T_steppe, rate=hg_mig_rate, matrix_index=(2, 3)),
    msprime.MigrationRateChange(time=T_steppe, rate=hg_mig_rate, matrix_index=(3, 2))
]

European_neolithic = [msprime.MassMigration(time=T_neo, source=0, dest=2, proportion=1.0 / 3.0)]

HG_split = [
    msprime.MassMigration(time=T_europe, source=3, dest=2, proportion=1),
    msprime.MigrationRateChange(time=T_europe, rate=0),
    msprime.PopulationParametersChange(time=T_europe, initial_size=N_A, population=2)
]

Near_east_split = [
    msprime.MassMigration(time=T_near_east, source=1, dest=0, proportion=1),
    msprime.PopulationParametersChange(time=T_near_east, initial_size=N_B, population=0)
]

Basal_split = [msprime.MassMigration(time=T_basal, source=2, dest=0, proportion=1)]

demographic_events = Modern_formation + Steppe_formation + European_neolithic + HG_split + Near_east_split + Basal_split

# Define samples
samples = []
indnames = []
indpopnames = []
for i, p in enumerate(populations):
    sample = [msprime.Sample(time=times[i], population=p)]
    samples = samples + sample * nhaps[i]
    tmp = [popnames[i] + "_" + str(j) for j in range(math.floor(nhaps[i] / 2))]
    for n in tmp:
        indnames.append(n)
    for n in tmp:
        indpopnames.append(popnames[i])

# modern_sample = [msprime.Sample(time=0, population=0)]
# neolithic_sample = [msprime.Sample(time=200, population=0)]
# steppe_sample = [msprime.Sample(time=180, population=1)]
# WHG_sample = [msprime.Sample(time=300, population=2)]
# EHG_sample = [msprime.Sample(time=300, population=3)]

# samples = modern_sample*10 + neolithic_sample*20 + steppe_sample*20 + WHG_sample*10 + EHG_sample*10 #Can multiply each by the number of samples needed for each population

# Debugging the demography
dd = msprime.DemographyDebugger(population_configurations=population_configurations, migration_matrix=migration_matrix,
                                demographic_events=demographic_events)
dd.print_history()

# Simulate chromosome 3 only
tree_sequence = msprime.simulate(
    length=198295000,
    recombination_rate=1e-8,
    mutation_rate=1.25e-8,
    population_configurations=population_configurations,
    demographic_events=demographic_events,
    samples=samples
)

tree_sequence.dump(
    path="/Users/williambarrie/" + args.out)  # Saving tree.sequence object to desired directory
tree = tree_sequence.first()
print(tree.draw(format="unicode"))  # Plot the first tree in the tree sequence

# Use tree_sequence to create finestructure input files

path = ("/Users/williambarrie/" + args.path)
if not os.path.exists(path):
    os.makedirs(path)

# pop_ids file
with open(os.path.join(path, "pop_ids"), "w") as file:
    for i in range(len(indnames)):
        file.write(indnames[i] + " " + indpopnames[i] + " 1\n")

# phase_file
list_of_inds = []
for sample in tree_sequence.samples():
    list_of_inds.append(sample)
number_of_inds = len(list_of_inds)

listofsites = []
for site in tree_sequence.sites():
    listofsites.append(math.floor(site.position))
numsnps = len(listofsites)

#Ensure all SNP sites are unique
uniquelistofsites = [listofsites[0]]
for i in range(1, len(listofsites)):
    if listofsites[i] <= listofsites[i - 1]:
        uniquelistofsites.append(listofsites[i - 1] + 1)
    else:
        uniquelistofsites.append(listofsites[i])

#Run again in case there are three the same in a row
uniquelistofsites2 = [uniquelistofsites[0]]
for i in range(1, len(uniquelistofsites)):
    if uniquelistofsites[i] <= uniquelistofsites[i - 1]:
        uniquelistofsites2.append(uniquelistofsites[i - 1] + 1)
    else:
        uniquelistofsites2.append(uniquelistofsites[i])

#Test if all sites are unique, otherwise error message
if uniquelistofsites2 != set(uniquelistofsites2):
    print("Error: some SNPs have same location - need to edit")

with open(os.path.join(path, "phasefile"), "w") as file:
    file.write(str(number_of_inds) + "\n")
    file.write(str(numsnps) + "\n")
    file.write("P")
    for site in uniquelistofsites2:
        file.write(" " + str(site))
    file.write("\n")
    for hap in tree_sequence.haplotypes():
        file.write(str(hap) + "\n")

# recomb_file
with open(os.path.join(path, "recombfile"), "w") as file:
    file.write("start.pos " + "recom.rate.perbp" + "\n")
    for site in uniquelistofsites2:
        file.write(str(site) + " " + "1e-7" + "\n")
