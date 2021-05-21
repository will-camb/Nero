import pandas as pd
import argparse
import seaborn as sns
import gzip

sns.set(color_codes=True)

parser = argparse.ArgumentParser()
parser.add_argument("-chr",
                    help="The chromosome that the sites of interest are on",
                    required=True)
parser.add_argument("-sites_file",
                    help="The positions of the sites for which average minor allele ancestry is to be calculated, "
                         "one per row",
                    required=True)
args = parser.parse_args()

sites_list = pd.read_csv(args.sites_file, header=None)[0].tolist()
sitesDf = pd.read_csv(args.sites_file, header=None)
idfile = pd.read_csv("/willerslev/ukbiobank/ordered_all_pop_ids_mapped", header=None, sep=" ", index_col=0)
# idfile = pd.read_csv("/willerslev/ukbiobank/painting_results_split/split_24001-48000/ordered_all_pop_ids_mapped", header=None, sep=" ", index_col=0)
cols = pd.read_csv("/willerslev/ukbiobank/painting_results_aggregate/copyprobs_per_anc/WHG." +
                   str(args.chr) + ".master_all_copyprobsperlocus.txt.gz", sep=" ", nrows=0).columns.tolist()
# cols = pd.read_csv("/willerslev/ukbiobank/painting_results_split/split_24001-48000/temp.WHG." +
#                    str(args.chr) + ".master_all_copyprobsperlocus.txt.gz", sep=" ", nrows=0).columns.tolist()
del cols[0]
cols = [int(x) for x in cols]
cols_reversed = cols[::-1]

phaselist = [cols.index(x) for x in sites_list]  # Get index of columns from phasefile to include
# phase = pd.read_csv("/willerslev/ukbiobank/painting_results_aggregate/phasefiles_transformed/" + str(args.chr) + ".merged.phase.gz", header=None, sep=" ", dtype='int8', usecols=phaselist)
with open("/willerslev/ukbiobank/painting_results_aggregate/phasefiles_transformed/" + str(args.chr) + ".merged.phase.gz", 'rb') as fd:
    gzip_fd = gzip.GzipFile(fileobj=fd)
    phase = pd.read_csv(gzip_fd, header=None, sep=" ", dtype='int8', usecols=phaselist)
# phase = pd.read_csv("/willerslev/ukbiobank/painting_results_split/split_24001-48000/phasefiles/transformed." + str(args.chr) +
#                     ".merged.phase.gz", header=None, sep=" ", dtype='int8', usecols=phaselist)
phase.columns = sites_list
phase.index = [val for val in idfile.index.unique().tolist() for _ in (0, 1)]
haps = list()
for h in range(int(phase.shape[0] / 2)):
    haps.extend([1, 2])

wrong_right_map = pd.DataFrame(cols)
wrong_right_map.loc[:, 1] = cols_reversed
wrong_right_map.columns = ['Wrong', 'Right']
mapped_positions = pd.merge(sitesDf, wrong_right_map, left_on=0, right_on='Right')  # Check this
cols_mapped = mapped_positions['Wrong'].tolist()
types_dict = {'0': str}
types_dict.update({str(col): 'int8' for col in cols_mapped if col not in types_dict})

results_list = []
results_dict = {}
# for anc in ["CHG", "EHG", "Farmer", "African", "EastAsian", "WHG", "Yamnaya"]:
for anc in ["Farmer", "Yamnaya"]:
    copyprobs = pd.read_csv("/willerslev/ukbiobank/painting_results_aggregate/copyprobs_per_anc/" + str(anc) + "." +
                            str(args.chr) + ".master_all_copyprobsperlocus.txt.gz",
                            sep=" ", dtype=types_dict, usecols=types_dict.keys())
    # copyprobs = pd.read_csv("/willerslev/ukbiobank/painting_results_split/split_24001-48000/temp." + str(anc) + "." +
    #                         str(args.chr) + ".master_all_copyprobsperlocus.txt.gz",
    #                         sep=" ", dtype=types_dict, usecols=types_dict.keys())
    copyprobs.set_index("0", inplace=True)
    copyprobs.columns = sites_list
    copyprobs["haps"] = haps
    for i in sites_list:
        phase_i = pd.DataFrame(zip(phase.index, haps, phase[i]))
        anc_copyprobs_temp_i = pd.DataFrame(zip(copyprobs.index, haps, copyprobs[i]))
        df = pd.merge(anc_copyprobs_temp_i, phase_i, left_on=[0, 1], right_on=[0, 1]).set_index(0)
        df.columns = ['haps', anc, 'phase']
        if df.phase.sum() < 408884:
        # if df.phase.sum() < 24000:
            #  1 is minor allele
            minor = 1
        else:
            #  1 is major allele
            minor = 0
        mean_minor_painting = df.loc[df['phase'] == minor][anc].mean()
        results_list.append(mean_minor_painting)
    results_dict[anc] = results_list

final_results = pd.DataFrame.from_dict(results_dict, orient='index')
final_results.columns = sites_list
final_results.to_csv("leave_1_out.results.csv", index=False)
