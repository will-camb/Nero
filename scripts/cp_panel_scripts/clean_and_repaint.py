import subprocess
import pandas as pd
import numpy as np

# bashCommand1 = "length=`head -n 1 dummy`"
# process = subprocess.Popen(bashCommand1.split(), stdout=subprocess.PIPE)
# output, error = process.communicate()
# awk '{ print length }' dummy_temp

# chrlist=`seq 1 22`
# for chr in $chrlist; do
#     length=`head -n 1 $chr.master_all_copyprobsperlocus.txt`
#     awk -F ',' -v nhaps="$length" 'length($3) == nhaps' $chr.master_all_copyprobsperlocus.txt > temp.$chr.master_all_copyprobsperlocus.txt
# done

def get_inds_to_paint(input):
    "Input a n.master_all_copyprobsperlocus.txt and write individuals to repaint to a list to_paint"
    painted = pd.read_csv(input, usecols=[0], header=None, skiprows=1)
    counts = painted[0].value_counts()
    counts.to_csv("counts")
    counts = pd.read_csv("counts", header=None)
    painted_successfully = counts[counts[1] == 20]
    painted_successfully_list = painted_successfully[0].tolist()
    to_paint_temp = np.setdiff1d(ordered_all_pop_ids_mapped_list, painted_successfully_list)
    to_paint.extend([s for s in to_paint_temp if "UKBB" in s])

ordered_all_pop_ids_mapped = pd.read_csv("ordered_all_pop_ids_mapped", header=None, sep=" ")
ordered_all_pop_ids_mapped_list = ordered_all_pop_ids_mapped[0].tolist()

to_paint = list()
for i in range(1,3):
    get_inds_to_paint("temp." + str(i) + ".master_all_copyprobsperlocus.txt")
to_paint = list(set(to_paint))

# for chr in $chrlist; do
#     rm temp.$chr.master_all_copyprobsperlocus.txt
# done

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
commands = commands[2].isin(to_paint)
commands.to_csv("repaint.paintvspanel1by1_commands.txt", sep=" ", header=False, index=False)
