import os
import argparse
import math
from run_ms_prime_v4_for_fs import RunMsPrime


class RunMsWithFsOutput:

    def __init__(self, nhaps=[2, 0, 0, 360, 30, 96, 84, 10, 26],
                 sample_times=[0, 149, 150, 187, 158, 266, 234, 314, 289],
                 hg_mig_rate=2e-3, mutation_rate=1.25e-8,
                 popnames=["modern", "bronze", "baa", "farmer", "yam", "WHG", "EHG", "farmeranatolian", "CHG"],
                 populations=[0, 0, 4, 0, 1, 2, 3, 0, 1],
                 ):
        self.nhaps = nhaps
        self.sample_times = sample_times
        self.hg_mig_rate = hg_mig_rate
        self.mutation_rate = mutation_rate
        self.popnames = popnames
        self.populations = populations

    def run(self):
        #  Args for script
        parser = argparse.ArgumentParser()
        parser.add_argument("-path_out",
                            help="Path to save tree sequence and FINESTRUCTURE input files to",
                            required=True)
        parser.add_argument("-chromosome",
                            help="Which chromosome to simulate",
                            required=True)
        parser.add_argument("-seed",
                            help="If this is None, automatically generated",
                            required=False)
        args = parser.parse_args()

        #  Run simulation
        tree_sequence = RunMsPrime(nhaps=self.nhaps, sample_times=self.sample_times, hg_mig_rate=self.hg_mig_rate,
                                   mutation_rate=self.mutation_rate, popnames=self.popnames,
                                   populations=self.populations, chromosome=args.chromosome).run_model()

        #  Make output directory and save tree sequence there
        path = args.path_out
        if not os.path.exists(path):
            os.makedirs(path)
        tree_sequence.dump(args.path_out + "/" + str(args.chromosome) + ".tree_sequence.ts")

        indnames = []
        indpopnames = []
        for i, p in enumerate(self.populations):
            tmp = [self.popnames[i] + str(j) for j in range(math.floor(self.nhaps[i] / 2))]
            for n in tmp:
                indnames.append(n)
                indpopnames.append(self.popnames[i])

        #  pop_ids file
        with open(os.path.join(path, "pop_ids"), "w") as file:
            for i in range(len(indnames)):
                file.write(indnames[i] + " " + indpopnames[i] + " 1\n")

        #  phase_file
        list_of_inds = [sample for sample in tree_sequence.samples()]
        number_of_inds = len(list_of_inds)
        listofsites = [math.floor(site.position) for site in tree_sequence.sites()]
        numsnps = len(listofsites)
        #  Ensure all site positions are unique
        while len(listofsites) != len(set(listofsites)):
            for i in range(1, len(listofsites)):
                if listofsites[i] <= listofsites[i - 1]:
                    listofsites[i] = (listofsites[i - 1] + 1)
        #  Test if all site positions are unique, otherwise error message
        if len(listofsites) != len(set(listofsites)):
            print("Error: some SNPs have same location - need to edit")
        with open(os.path.join(path, str(args.chromosome) + ".phasefile"), "w") as file:
            file.write(str(number_of_inds) + "\n")
            file.write(str(numsnps) + "\n")
            file.write("P")
            for site in listofsites:
                file.write(" " + str(site))
            file.write("\n")
            for hap in tree_sequence.haplotypes():
                file.write(str(hap) + "\n")

        # recomb_file
        with open(os.path.join(path, str(args.chromosome) + ".recombfile"), "w") as file:
            file.write("start.pos " + "recom.rate.perbp" + "\n")
            for site in listofsites:
                file.write(str(site) + " " + "1e-7" + "\n")
