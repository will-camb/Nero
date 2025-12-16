import pandas as pd


for name in ['smallpox', 'HIV', 'plague', 'measles', 'leprosy', 'B_recurrentis', 'T_pallidum', 'L_interrogans', 'M_tuberculosis', 'P_vivax']:
	df = pd.read_csv(f'{name}_analysis/{name}_hla_binding_results.csv')
	df_i = df.loc[df['hla_class']=='I']
	df_i = df_i.loc[df_i['%Rank_EL'] < 0.05] # 'strong' < 0.5
	df_i = df_i.groupby('hla_allele')['protein_id'].nunique().sort_values()
	df_i.to_csv(f"analysed_results/{name}_class_i_hla_binding_results_analysed.tsv", sep="\t")

	df_ii = df.loc[df['hla_class']=='II']
	df_ii = df_ii.loc[df_ii['%Rank_EL'] < 0.1] # 'strong' < 1
	df_ii = df_ii.groupby('hla_allele')['protein_id'].nunique().sort_values()
	df_ii.to_csv(f"analysed_results/{name}_class_ii_hla_binding_results_analysed.tsv", sep="\t")