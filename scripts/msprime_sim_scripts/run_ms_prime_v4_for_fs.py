#  Updated msprime model to reflect samples in mesoneo painting of UKBiobank
#  Generation times: assume 30 years average per generation here

import math
import msprime


class RunMsPrime:

    def __init__(self, nhaps=[0, 0, 0, 360, 30, 96, 84, 10, 26],
                 sample_times=[0, 149, 150, 187, 158, 266, 234, 314, 289],
                 hg_mig_rate=2e-3, length=198295559, recombination_rate=1e-8, mutation_rate=1.25e-8,
                 popnames=["modern", "bronze", "baa", "farmer", "yam", "WHG", "EHG", "ana", "CHG"],
                 populations=[0, 0, 4, 0, 1, 2, 3, 0, 1]
                 ):
        self.nhaps = nhaps
        self.sample_times = sample_times
        self.hg_mig_rate = hg_mig_rate
        self.length = length
        self.recombination_rate = recombination_rate
        self.mutation_rate = mutation_rate
        self.popnames = popnames
        self.populations = populations

        infile = "/Users/williambarrie/Desktop/genetic_map_GRCh37_chr3.txt"
        self.recomb_map = msprime.RecombinationMap.read_hapmap(infile)


    def run_model(self):

        # initial population sizes:
        N_bronze = 50000
        N_Yam = 20000
        N_baa = 10000
        N_whg = 10000
        N_ehg = 10000
        N_neo = 50000
        N_chg = 10000
        N_A = 5000  # Ancestor of WHG and EHG
        N_B = 5000  # Ancestor of CHG and Neolithic farmers

        # Time of events
        T_bronze = 150
        T_Yam = 200
        T_neo = 250
        T_baa = 275
        T_near_east = 800
        T_europe = 500
        T_basal = 1500

        # Growth rate and initial population size for present day from bronze age
        r_EU = 0.067
        N_present = N_bronze / math.exp(-r_EU * T_bronze)

        #Populations: 0=present/bronze/neolithic_farmers/Ana/B,1=Yam/CHG,2=WHG/A, 3=EHG, 4=BAA
        population_configurations = [
            msprime.PopulationConfiguration(initial_size=N_present, growth_rate=r_EU),
            msprime.PopulationConfiguration(initial_size=N_Yam),
            msprime.PopulationConfiguration(initial_size=N_whg),
            msprime.PopulationConfiguration(initial_size=N_ehg),
            msprime.PopulationConfiguration(initial_size=N_baa)
        ]
        bronze_formation = [
            msprime.MassMigration(time=T_bronze, source=0, dest=1, proportion=0.5),
            msprime.PopulationParametersChange(time=T_bronze, initial_size=N_neo, growth_rate=0, population=0)
        ]
        yam_formation = [
            msprime.MassMigration(time=T_Yam, source=1, dest=3, proportion=0.5),
            msprime.PopulationParametersChange(time=T_Yam, initial_size=N_chg, population=1),
            msprime.MigrationRateChange(time=T_Yam, rate=self.hg_mig_rate, matrix_index=(2, 3)),
            msprime.MigrationRateChange(time=T_Yam, rate=self.hg_mig_rate, matrix_index=(3, 2))
        ]
        european_neolithic = [msprime.MassMigration(time=T_neo, source=0, dest=2, proportion=1.0 / 4.0)]
        baa_formation = [msprime.MassMigration(time=T_baa, source=4, dest=1, proportion=1.0 / 4.0)]
        ana_split = [msprime.MassMigration(time=276, source=4, dest=0, proportion=1)]
        hg_split = [
            msprime.MassMigration(time=T_europe, source=3, dest=2, proportion=1),
            msprime.MigrationRateChange(time=T_europe, rate=0),
            msprime.PopulationParametersChange(time=T_europe, initial_size=N_A, population=2)
        ]
        near_east_split = [
            msprime.MassMigration(time=T_near_east, source=1, dest=0, proportion=1),
            msprime.PopulationParametersChange(time=T_near_east, initial_size=N_B, population=0)
        ]
        basal_split = [msprime.MassMigration(time=T_basal, source=2, dest=0, proportion=1)]
        demographic_events = bronze_formation + yam_formation + european_neolithic + baa_formation + ana_split + hg_split + near_east_split + basal_split

        # Define samples
        samples = []
        for i, p in enumerate(self.populations):
            sample = [msprime.Sample(time=self.sample_times[i], population=p)]
            samples = samples + sample * self.nhaps[i]

        # Debugging the demography
        migration_matrix = [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]
        dd = msprime.DemographyDebugger(population_configurations=population_configurations,
                                        migration_matrix=migration_matrix, demographic_events=demographic_events)
        dd.print_history()

        # Simulate chromosome 3 only
        tree_sequence = msprime.simulate(
            # length=self.length,
            # recombination_rate=self.recombination_rate,
            mutation_rate=self.mutation_rate,
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            samples=samples,
            recombination_map=self.recomb_map
        )
        return tree_sequence
