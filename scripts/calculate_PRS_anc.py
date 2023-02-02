def calculate_PRS_anc(data, ancestry):
    results_df = pd.DataFrame(columns=['phenotype', 'ancestry', 'PRS'])
    maf_beta = pd.DataFrame()
    bootstrap = data
    dict_of_PRS = dict()
    maf_df_all = pd.DataFrame()
    for anc in anc_list:
        df1 = df.loc[(df[1] == phen) & (df[8] == 0) & (df[2] == anc)]
        maf_list = list()
        for i in range(1, 23):
            try:
                maf_list += df1.loc[(df1[3] == i)][9].item()
            except ValueError:
                continue
        beta_position_map['maf'] = maf_list
        maf_beta = pd.merge(bootstrap, beta_position_map)
        maf_beta['meangen'] = maf_beta['maf'] * maf_beta['beta']
        dict_of_PRS[anc] = maf_beta['meangen'].sum()
        maf_df_all[anc] = maf_beta['maf']
    maf_df_all['beta'] = maf_beta['beta']
    maf_df_all['mean_freq'] = maf_df_all[anc_list].mean(axis=1)
    maf_df_all['score'] = maf_df_all['mean_freq']*(1-maf_df_all['mean_freq']) * (maf_df_all['beta']**2)
    mean_PRS = statistics.mean(list(dict_of_PRS.values()))
    varmean = maf_df_all['score'].sum()
    for anc in anc_list:
        meangen = dict_of_PRS[anc] - mean_PRS
        meangenvec = 2 * meangen / (math.sqrt(4*varmean))
        results_df = results_df.append({'phenotype': phen, 'ancestry': anc, 'PRS': meangenvec}, ignore_index=True)
    return results_df.loc[results_df['ancestry'] == ancestry]['PRS'].item()

