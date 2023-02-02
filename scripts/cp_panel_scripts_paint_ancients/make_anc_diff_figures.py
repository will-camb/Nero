import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import random
from scipy.stats import zscore
sns.set_style("dark")
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-chr",
                    help="the chromosome, NB must be named in form $chr.master_all_copyprobsperlocus.txt",
                    required=True)
args = parser.parse_args()

name2id = pd.read_csv("name2id", header=None, sep=" ")
final_metadata = pd.read_csv("final_metadata.csv")
final_metadata_id = pd.merge(name2id, final_metadata, left_on=0, right_on='fam')
NNLS_combined_metadata = pd.read_csv("NNLS_combined_metadata.csv")

denmark_medieval = NNLS_combined_metadata[(NNLS_combined_metadata.country == 'Denmark') & (NNLS_combined_metadata.ageAverage <= 850) & (NNLS_combined_metadata.ageAverage > 100)]
denmark_medieval_ids = denmark_medieval['Unnamed: 0'].tolist()

denmark_bronze_iron_viking = NNLS_combined_metadata[NNLS_combined_metadata.groupLabel.isin(['Denmark_BronzeAge', 'Denmark_IronAge', 'Denmark_VikingAge'])]
denmark_bronze_iron_viking_ids = denmark_bronze_iron_viking['Unnamed: 0'].tolist()

europe_mesolithic = NNLS_combined_metadata[(NNLS_combined_metadata.groupLabel.str.contains('Mesolithic')) & (NNLS_combined_metadata.region.str.contains('Europe'))]
europe_mesolithic_ids = europe_mesolithic['Unnamed: 0'].tolist()

europe_neolithic = NNLS_combined_metadata[(NNLS_combined_metadata.groupLabel.str.contains('Neolithic')) & (NNLS_combined_metadata.region.str.contains('Europe'))]
europe_neolithic_ids = europe_neolithic['Unnamed: 0'].tolist()

anc_copyprobs = pd.read_csv(f"tranformed_copyprobs/{args.chr}.master_all_copyprobsperlocus.txt.gz", index_col=0)
CHG = anc_copyprobs.iloc[0::10, :]
EHG = anc_copyprobs.iloc[1::10, :]
FarmerAnatolian = anc_copyprobs.iloc[2::10, :]
FarmerEarly = anc_copyprobs.iloc[3::10, :]
FarmerLate = anc_copyprobs.iloc[4::10, :]
FarmerMiddle = anc_copyprobs.iloc[5::10, :]
EastAsian = anc_copyprobs.iloc[6::10, :]
African = anc_copyprobs.iloc[7::10, :]
WHG = anc_copyprobs.iloc[8::10, :]
Yamnaya = anc_copyprobs.iloc[9::10, :]
Farmer = FarmerAnatolian + FarmerEarly + FarmerLate + FarmerMiddle
Steppe = Yamnaya + EastAsian

difference_per_anc = pd.DataFrame(index=Yamnaya.columns)
early_ids = europe_neolithic_ids
late_ids = denmark_medieval_ids

# for ancestry_str in ['African', 'CHG', 'EHG', 'EastAsian', 'FarmerAnatolian', 'FarmerEarly','FarmerLate', 'FarmerMiddle', 'WHG', 'Yamnaya']:
# for ancestry_str in ['African', 'CHG', 'EHG', 'EastAsian', 'Farmer', 'WHG', 'Yamnaya', 'Steppe']:
for ancestry_str in ['CHG', 'EHG','Farmer', 'WHG', 'Steppe']:
    ancestry = eval(ancestry_str)
    late_mean = ancestry.loc[ancestry.index.isin(late_ids)].mean()
    early_mean = ancestry.loc[ancestry.index.isin(early_ids)].mean()
    merged = pd.DataFrame(zip(early_mean, late_mean), index=early_mean.index, columns=['early', 'late'])
    merged = merged.apply(zscore) # z-score for each ancestry average
    merged[ancestry_str] = merged['late'] - merged['early']
    difference_per_anc[ancestry_str] = merged[ancestry_str]
difference_per_anc.reset_index(inplace=True)
difference_per_anc['index'] = difference_per_anc['index'].astype(float)

# Get position labels
df = pd.read_csv("../angelos_snps/ms_snps_final_discovery_0.7_combined_with_effect_rs.csv")
ms_positions = list(df.loc[(df['proxy_chr']==int(args.chr))][['finemapped_position','Patsopoulos_effect_rs']].itertuples(index=False))

# Plot
fig,ax = plt.subplots(figsize=[35,10])
sns.lineplot(x='index', y='value', hue='variable',
             data=pd.melt(difference_per_anc, ['index']), ax=ax)
ax.axhline(0, color='black')
plt.ylim(-7.5, 7.5)
plt.ylabel(f'Change in ancestry (z-score)')
plt.xlabel(f'Chromosome {args.chr}')
# plt.title("Europe Mesolithic - Denmark Medieval", fontsize=40)

for position,effect_rs in ms_positions:
    ax.axvline(position, color='black', alpha=0.1)
    plt.text(position+500,3.5+random.uniform(0,1.5), effect_rs, fontsize=15)
plt.savefig(f"Ancestry_difference_figures/EuropeNeo-DenmarkMed_chr{args.chr}.png")
plt.ylabel(f'Change in ancestry')
plt.xlabel(f'Chromosome {args.chr}')

print(f"Done chr {args.chr}!")