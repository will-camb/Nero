#!/bin/bash
chrlist=`seq 1 22`
for chr in $chrlist; do
    length=`head -n 1 $chr.master_all_copyprobsperlocus.txt`
    awk -F ',' -v nsnps="$length" 'length($3) == nsnps' $chr.master_all_copyprobsperlocus.txt > temp.$chr.master_all_copyprobsperlocus.txt
    cut -f1 -d"," temp.$chr.master_all_copyprobsperlocus.txt | uniq -c > temp.{$chr}.counts
done

python3 << END
import pandas as pd
import numpy as np

ordered_all_pop_ids_mapped = pd.read_csv("ordered_all_pop_ids_mapped", header=None, sep=" ")
ordered_all_pop_ids_mapped_list = ordered_all_pop_ids_mapped[0].tolist()

to_paint = list()
for i in range(1,23):
    painted = pd.read_csv("temp."+str(i)+".counts", sep=" ", header=None, error_bad_lines=False)
    painted = painted[[5,6]]
    painted[5] = painted[5].astype(str)
    painted_successfully=painted[painted[5]=='20']
    painted_successfully_list = painted_successfully[6].tolist()
    to_paint_temp = np.setdiff1d(ordered_all_pop_ids_mapped_list, painted_successfully_list)
    to_paint.extend([s for s in to_paint_temp if "UKBB" in s])
    print("Done assessment of chr " + str(i) + " and results written to to_paint; now processing chr " + str(i+1))
to_paint = list(set(to_paint))
chunkcounts_list = pd.read_csv("ordered_all_pop_ids_mapped.allchr.chunkcounts.out", header=None, sep=" ", error_bad_lines=False)[0].tolist()
chunklengths_list = pd.read_csv("ordered_all_pop_ids_mapped.allchr.chunklengths.out", header=None, sep=" ", error_bad_lines=False)[0].tolist()
mutationprobs_list = pd.read_csv("ordered_all_pop_ids_mapped.allchr.mutationprobs.out", header=None, sep=" ", error_bad_lines=False)[0].tolist()
regionchunkcounts_list = pd.read_csv("ordered_all_pop_ids_mapped.allchr.regionchunkcounts.out", header=None, sep=" ", error_bad_lines=False)[0].tolist()
regionsquaredchunkcounts_list = pd.read_csv("ordered_all_pop_ids_mapped.allchr.regionsquaredchunkcounts.out", header=None, sep=" ", error_bad_lines=False)[0].tolist()
to_paint.extend(np.setdiff1d(ordered_all_pop_ids_mapped_list, chunkcounts_list))
to_paint.extend(np.setdiff1d(ordered_all_pop_ids_mapped_list, chunklengths_list))
to_paint.extend(np.setdiff1d(ordered_all_pop_ids_mapped_list, mutationprobs_list))
to_paint.extend(np.setdiff1d(ordered_all_pop_ids_mapped_list, regionchunkcounts_list))
to_paint.extend(np.setdiff1d(ordered_all_pop_ids_mapped_list, regionsquaredchunkcounts_list))
to_paint = list(set(to_paint))
commands = pd.read_csv("paintvspanel1by1_commands.txt", header=None, sep=" ")
commands = commands[commands[2].isin(to_paint)]
commands.to_csv("repaint.paintvspanel1by1_commands.txt", sep=" ", header=False, index=False)
print("now run repaint.paintvspanel1by1_commands.txt e.g. cat repaint.paintvspanel1by1_commands.txt | parallel")
END

# rm -r temp.*

commands = pd.read_csv("paintvspanel1by1_commands.txt", header=None, sep=" ")
commands = commands[commands[2].isin(to_paint)]
commands.to_csv("repaint.paintvspanel1by1_commands.txt", sep=" ", header=False, index=False)
print("now run repaint.paintvspanel1by1_commands.txt e.g. cat repaint.paintvspanel1by1_commands.txt | parallel")

for chr in $chrlist; do
    tmpfile=$(mktemp)
    cp file "$tmpfile" &&
    awk '...some program here...' "$tmpfile" >file
    rm "$tmpfile"