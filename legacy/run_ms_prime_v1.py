import msprime


class RunMsPrime:

    def __init__(self, nhaps=[10, 20, 20, 10, 10], sample_times=[0, 200, 180, 200, 200], hg_mig_rate=2.5e-5,
                 length=198295000, recombination_rate=1e-8, mutation_rate=1.25e-8,
                 popnames=["modern", "neolithic", "steppe", "WHG", "EHG"], populations=[0, 0, 1, 2, 3]):
        self.nhaps = nhaps
        self.sample_times = sample_times
        self.hg_mig_rate = hg_mig_rate
        self.length = length
        self.recombination_rate = recombination_rate
        self.mutation_rate = mutation_rate
        self.popnames = popnames
        self.populations = populations

    def run_model(self):
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
            msprime.MigrationRateChange(time=T_steppe, rate=self.hg_mig_rate, matrix_index=(2, 3)),
            msprime.MigrationRateChange(time=T_steppe, rate=self.hg_mig_rate, matrix_index=(3, 2))
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
        for i, p in enumerate(self.populations):
            sample = [msprime.Sample(time=self.sample_times[i], population=p)]
            samples = samples + sample * self.nhaps[i]

        # Simulate chromosome 3 only
        tree_sequence = msprime.simulate(
            length=self.length,
            recombination_rate=self.recombination_rate,
            mutation_rate=self.mutation_rate,
            population_configurations=population_configurations,
            demographic_events=demographic_events,
            samples=samples
        )
        return tree_sequence
