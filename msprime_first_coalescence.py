import msprime
import math
import tskit
import pandas as pd
import numpy as np
import argparse


def nan_all_but_min(df):
    arr = df.values
    idx = np.argmin(arr, axis=1)
    newarr = np.full_like(arr, np.nan, dtype='float')
    newarr[np.arange(arr.shape[0]), idx] = arr[np.arange(arr.shape[0]), idx]
    df = pd.DataFrame(newarr, columns=df.columns, index=df.index)
    return df


parser = argparse.ArgumentParser()
parser.add_argument("-input_ts",
                    help="Path and filename of msprime output (tree_sequence)",
                    required=True)
parser.add_argument("-out",
                    help="Path to save msprime first coalescence csv to",
                    required=False)
args = parser.parse_args()

ts = tskit.load(args.input_ts)

sample_ids = [i.id for i in ts.nodes() if i.flags==msprime.NODE_IS_SAMPLE]

df = pd.DataFrame(columns = (x for x in range(len(sample_ids))))
dictofdata = dict()
index_list = list()

intervallist2 = list()
for tree in ts.aslist():
    for sample_id1 in sample_ids:
        index_list.append('tree_' + str(tree.index) + '_HAP_' + str(sample_id1))
        for sample_id2 in sample_ids:
            if sample_id1 == sample_id2:
                data = math.inf
                dictofdata[sample_id2] = data
            else:
                data = tree.mrca(sample_id1, sample_id2)
                dictofdata[sample_id2] = data
        df = df.append(dictofdata, ignore_index=True)

df.index = index_list

# Produce df with Nan values for all but the nearest neighbour mrca
df = nan_all_but_min(df=df)

# Produce df with mrca, nearest_neighbour and tree interval for each hap for each tree
df3 = pd.DataFrame(index=index_list)
list_of_mrca = list()
list_of_nearest_neighbour = list()

for index, row in df.iterrows():
    for x in df.columns:
        if str(row[x]) != 'nan':
            list_of_mrca.append(row[x])
            list_of_nearest_neighbour.append('hap' + str(x))
df3['mrca'] = list_of_mrca
df3['nearest_neighbour'] = list_of_nearest_neighbour

intervallist2 = list()
for tree in ts.aslist():
    for sample_id1 in sample_ids:
        intervallist2.append(tree.interval)
df3['interval'] = intervallist2


# Add hap_id and tree columns to df based on existing index
split_list = list(df3.index)
split_hap_list = list()
split_tree_list = list()
for x in split_list:
    y = x.split('_')
    split_hap_list.append(y[2] + '_' + y[3])
    split_tree_list.append(y[0] + '_' + y[1])

df3['hap_id'] = split_hap_list
df3['tree'] = split_tree_list

# Make df with one row per variant per hap_id
variant_list = list()
hap_ids_list = list()
df4 = pd.DataFrame()
for sample_id in sample_ids:
    for variant in ts.variants():
        variant_list.append(math.floor(variant.site.position))
        hap_ids_list.append('HAP_' + str(sample_id))
df4['pos'] = variant_list
df4['hap_id'] = hap_ids_list

nearest_neighbour_list = list()
tree_list = list()
for pos,hap_id1 in zip(df4['pos'], df4['hap_id']):
    for nearest_neighbour,interval,hap_id2,tree in zip(df3['nearest_neighbour'],df3['interval'],df3['hap_id'],df3['tree']):
        if pos >= interval[0] and pos<= interval[1] and hap_id1 == hap_id2:
            nearest_neighbour_list.append(nearest_neighbour)
            tree_list.append(tree)

df4['ancestry'] = nearest_neighbour_list
df4['tree'] = tree_list

df4.to_csv(args.out + '/df_msprime_final_ancestry.csv')