import os
import argparse
import math
from scripts.run_ms_prime_v3_for_fs import RunMsPrime


class RunMsWithFsOutput:

    def __init__(self, nhaps=[200, 102, 12, 162, 14, 86, 74, 20, 24],
                 sample_times=[0, 149, 150, 200, 175, 300, 275, 251, 250],
                 hg_mig_rate=2e-3, length=198295559, mutation_rate=1.25e-8,
                 popnames=["modern", "bronze", "baa", "neolithic", "yam", "WHG", "EHG", "ana", "CHG"],
                 populations=[0, 0, 4, 0, 1, 2, 3, 0, 1], infile = "/willerslev/datasets/hapmapRecomb/2011-01_phaseII_B37/genetic_map_GRCh37_chr3.txt"
                 ):
        self.nhaps = nhaps
        self.sample_times = sample_times
        self.hg_mig_rate = hg_mig_rate
        self.length = length
        self.mutation_rate = mutation_rate
        self.popnames = popnames
        self.populations = populations

    def run(self):
        #run simulation
        tree_sequence = RunMsPrime(nhaps=self.nhaps, sample_times=self.sample_times, hg_mig_rate=self.hg_mig_rate,
                 length=self.length, mutation_rate=self.mutation_rate,
                 popnames=self.popnames, populations=self.populations).run_model()

        # args for script
        parser = argparse.ArgumentParser()
        parser.add_argument("-out_ms",
                            help="Path and filename to save the tree sequence to e.g. Documents/msprimetest NB directory must already exist",
                            required=True)
        parser.add_argument("-path_fs",
                            help="Path to save FINESTRUCTURE input files to e.g. Documents/msprime2fstest",
                            required=True)
        parser.add_argument("-seed",
                            help="If this is None, automatically generated",
                            required=False)
        args = parser.parse_args()
        tree_sequence.dump(path=args.out_ms)
        path = args.path_fs
        if not os.path.exists(path):
            os.makedirs(path)

        indnames = []
        indpopnames = []
        for i, p in enumerate(self.populations):
            tmp = [self.popnames[i] + "_" + str(j) for j in range(math.floor(self.nhaps[i] / 2))]
            for n in tmp:
                indnames.append(n)
                indpopnames.append(self.popnames[i])

        # pop_ids file
        with open(os.path.join(path, "pop_ids"), "w") as file:
            for i in range(len(indnames)):
                file.write(indnames[i] + " " + indpopnames[i] + " 1\n")

        # phase_file
        list_of_inds = [sample for sample in tree_sequence.samples()]
        number_of_inds = len(list_of_inds)
        listofsites = [math.floor(site.position) for site in tree_sequence.sites()]
        numsnps = len(listofsites)
        # Ensure all site positions are unique
        while len(listofsites) != len(set(listofsites)):
            for i in range(1, len(listofsites)):
                if listofsites[i] <= listofsites[i - 1]:
                    listofsites[i] = (listofsites[i - 1] + 1)
        # Test if all site positions are unique, otherwise error message
        if len(listofsites) != len(set(listofsites)):
            print("Error: some SNPs have same location - need to edit")
        with open(os.path.join(path, "phasefile"), "w") as file:
            file.write(str(number_of_inds) + "\n")
            file.write(str(numsnps) + "\n")
            file.write("P")
            for site in listofsites:
                file.write(" " + str(site))
            file.write("\n")
            for hap in tree_sequence.haplotypes():
                file.write(str(hap) + "\n")

        # recomb_file
        # Load recomb map
        infile = "/willerslev/datasets/hapmapRecomb/2011-01_phaseII_B37/genetic_map_GRCh37_chr3.txt"
        recomb_map = msprime.RecombinationMap.read_hapmap(infile)

        with open(os.path.join(path, "recombfile"), "w") as file:
            file.write("start.pos " + "recom.rate.perbp" + "\n")
            for site in listofsites:
                file.write(str(site) + " " + "1e-7" + "\n")
