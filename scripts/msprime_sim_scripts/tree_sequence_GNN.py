import math
import tskit
import msprime
import argparse
import pandas as pd
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("-input_ts",
                    help="tree sequence input file",
                    required=True)
parser.add_argument("-ts_number_modern_haps",
                    help="the number of modern haplotypes that are sampled in the tree sequence",
                    required=True)
parser.add_argument("-GNN__modern_haps",
                    help="the number of modern haplotypes that we want to record the GNNs for. Must start at 0",
                    required=True)
parser.add_argument("-pop_ids",
                    help="the pop ids file from the tree sequence - with pop ids in the same order as the samples",
                    required=True)
args = parser.parse_args()
ts = tskit.load(args.input_ts)

sample_ids = [i.id for i in ts.nodes() if i.flags == msprime.NODE_IS_SAMPLE]
modern_ids = sample_ids[0:args.GNN_modern_haps]
ancient_ids = sample_ids[args.ts_number_modern_haps:]

df = pd.DataFrame(columns=(x for x in range(len(ancient_ids))))
dictofdata = dict()
index_list = list()

for tree in ts.aslist():
    for modern_id in modern_ids:
        index_list.append('HAP_' + str(modern_id) + '_tree_' + str(tree.index))
        for ancient_id in ancient_ids:
            dictofdata[ancient_id] = tree.mrca(modern_id, ancient_id)
        df = df.append(dictofdata, ignore_index=True)

df.index = index_list

# produces df with Nan values for all but the nearest neighbours
df = df.where(df.values == df.min(axis=1)[:,None])

# Produces df with mrca, nearest_neighbour and tree interval for each hap for each tree
df2 = pd.DataFrame(index=index_list)
set_of_mrca = {}
dict_of_nearest_neighbour = dict()

for index, row in df.iterrows():
    list_of_nearest_neighbour = list()
    list_of_mrca = list()
    for x in df.columns:
        if str(row[x]) != 'nan':
            set_of_mrca.update(row[x])
            list_of_nearest_neighbour.append('hap' + str(x))
    dict_of_mrca['mrca'] = set_of_mrca
    dict_of_nearest_neighbour['GNNs'] = list_of_nearest_neighbour
    df2 = df2.append(dict_of_mrca, ignore_index=True)
    df2 = df2.append(dict_of_nearest_neighbour, ignore_index=True)

list_of_intervals = list()
for tree in ts.aslist():
    for sample_id in modern_ids + ancient_ids:
        list_of_intervals.append(tree.interval)
df2['interval'] = list_of_intervals

split_list = list(df2.index)
split_hap_list = list()
split_tree_list = list()
for x in split_list:
    y = x.split('_')
    split_hap_list.append(y[2] + '_' + y[3])
    split_tree_list.append(y[0] + '_' + y[1])

df2['hap_id'] = split_hap_list
df2['tree'] = split_tree_list

variant_list = list()
hap_ids_list = list()
df3 = pd.DataFrame()
for sample_id in sample_ids:
    for variant in ts.variants():
        variant_list.append(math.floor(variant.site.position))
        hap_ids_list.append('HAP_' + str(sample_id))
df3['pos'] = variant_list
df3['hap_id'] = hap_ids_list

nearest_neighbour_list = list()
tree_list = list()
for pos,hap_id1 in zip(df3['pos'], df3['hap_id']):
    for nearest_neighbour,interval,hap_id2,tree in zip(df2['nearest_neighbour'],df2['interval'],df2['hap_id'],df2['tree']):
        if pos >= interval[0] and pos<= interval[1] and hap_id1 == hap_id2:
            nearest_neighbour_list.append(nearest_neighbour)
            tree_list.append(tree)

df3['ancestry'] = nearest_neighbour_list
df3['tree'] = tree_list

