import pandas as pd

for name in ['smallpox', 'HIV', 'plague', 'measles', 'leprosy', 'B_recurrentis', 'T_pallidum', 'L_interrogans', 'M_tuberculosis', 'P_vivax', 'EBV', 'JCV', 'Salmonella_typhi', 'HSV-2', 'MCPyV', 'VZV']:
	df = pd.read_csv(f'{name}_analysis/{name}_hla_binding_results.csv')
	
	df_i = df.loc[df['hla_class']=='I']
	results_df = pd.DataFrame(columns=['hla_allele', 'Weak binders (%)', 'Strong binders (%)'])
	for HLA in df_i['hla_allele'].unique().tolist():
		df_i_HLA = df_i.loc[df_i['hla_allele'] == HLA]
		n_total = df_i_HLA.shape[0]
		n_strong = df_i_HLA.loc[df_i_HLA['%Rank_EL'] <= 0.5].shape[0]
		percent_strong = (n_strong / n_total)*100
		n_weak = df_i_HLA.loc[df_i_HLA['%Rank_EL'] <= 2].shape[0]
		percent_weak = (n_weak / n_total)*100
		results_df.loc[len(results_df)] = [HLA, percent_weak, percent_strong]
	results_df.sort_values('Strong binders (%)', ascending=False).to_csv(f"analysed_results/{name}_class_i_hla_binding_results_analysed.tsv", sep="\t", index=False)

	df_ii = df.loc[df['hla_class']=='II']
	results_df = pd.DataFrame(columns=['hla_allele', 'Weak binders (%)', 'Strong binders (%)'])
	for HLA in df_ii['hla_allele'].unique().tolist():
		df_ii_HLA = df_ii.loc[df_ii['hla_allele'] == HLA]
		n_total = df_ii_HLA.shape[0]
		n_strong = df_ii_HLA.loc[df_ii_HLA['%Rank_EL'] <= 1].shape[0]
		percent_strong = (n_strong / n_total)*100
		n_weak = df_ii_HLA.loc[df_ii_HLA['%Rank_EL'] <= 5].shape[0]
		percent_weak = (n_weak / n_total)*100
		results_df.loc[len(results_df)] = [HLA, percent_weak, percent_strong]
	results_df.sort_values('Strong binders (%)', ascending=False).to_csv(f"analysed_results/{name}_class_ii_hla_binding_results_analysed.tsv", sep="\t", index=False)



	# df_i = df_i.loc[df_i['%Rank_EL'] <= 0.05] # 'strong' < 0.5; 'weak' < 2
	# df_i = df_i.groupby('hla_allele')['protein_id'].nunique().sort_values()
	# df_i.to_csv(f"analysed_results/{name}_class_i_hla_binding_results_analysed.tsv", sep="\t")

	# df_ii = df.loc[df['hla_class']=='II']
	# df_ii = df_ii.loc[df_ii['%Rank_EL'] < 0.1] # 'strong' < 1; 'weak' < 5
	# df_ii = df_ii.groupby('hla_allele')['protein_id'].nunique().sort_values()
	# df_ii.to_csv(f"analysed_results/{name}_class_ii_hla_binding_results_analysed.tsv", sep="\t")