#!/bin/bash
#cd $PBS_O_WORKDIR
module load anaconda3/4.4.0

if [ "$#" -ne "5" ] ; then
    echo "Usage: run_PRS_calculator.sh <manifest> <file_name> <copyprobs_directory> <phasefile_directory>"
    echo "<manifest>: the manifest csv file"
    echo "<file_name>: filename from manifest to be calculated"
    echo "<copyprobs_directory>: Directory of all_copyprobsperlocus files; should be named in form n.master_all_copyprobsperlocus.txt.gz"
    echo "<phasefile_directory>: Directory of phasefiles; should be named in form chr#.merged.phase"
    echo "<idfile>: Path to idfile: ordered_all_pop_ids_mapped ('individualID popID 1')"
    exit 0
fi

manifest="$1"
file_name="$2"
copyprobs_directory="$3"
phasefile_directory="$4"
idfile="$5"
echo "*** Extracting results file from Dropbox ***"
#`awk -F ',' -v file_name="$file_name" '$5 == file_name { print $6 }' $manifest`
MY_FILE=$file_name
NEW_EXT=${MY_FILE/bgz/gz}
#mv $file_name "$NEW_EXT"
echo "*** Success! Now making temp copyprobs files per ancestry ***"
chrlist=`seq 1 22`
for chr in $chrlist; do
awk 'NR % 10 == 1' $3/$chr.master_all_copyprobsperlocus.txt > temp.CHG.$chr.master_all_copyprobsperlocus.txt
awk 'NR % 10 == 2' $3/$chr.master_all_copyprobsperlocus.txt > temp.EHG.$chr.master_all_copyprobsperlocus.txt
awk 'NR % 10 == 3' $3/$chr.master_all_copyprobsperlocus.txt > temp.FarmerAnatolian.$chr.master_all_copyprobsperlocus.txt
awk 'NR % 10 == 4' $3/$chr.master_all_copyprobsperlocus.txt > temp.FarmerEarly.$chr.master_all_copyprobsperlocus.txt
awk 'NR % 10 == 5' $3/$chr.master_all_copyprobsperlocus.txt > temp.FarmerLate.$chr.master_all_copyprobsperlocus.txt
awk 'NR % 10 == 6' $3/$chr.master_all_copyprobsperlocus.txt > temp.FarmerMiddle.$chr.master_all_copyprobsperlocus.txt
awk 'NR % 10 == 7' $3/$chr.master_all_copyprobsperlocus.txt > temp.EastAsian.$chr.master_all_copyprobsperlocus.txt
awk 'NR % 10 == 8' $3/$chr.master_all_copyprobsperlocus.txt > temp.African.$chr.master_all_copyprobsperlocus.txt
awk 'NR % 10 == 9' $3/$chr.master_all_copyprobsperlocus.txt > temp.WHG.$chr.master_all_copyprobsperlocus.txt
awk 'NR % 10 == 0' $3/$chr.master_all_copyprobsperlocus.txt > temp.Yamnaya.$chr.master_all_copyprobsperlocus.txt
done
python3 << END
import pandas as pd
for i in range(1,23):
  copyprobsFA = pd.read_csv("temp.FarmerAnatolian."+str(i)+".master_all_copyprobsperlocus.txt",
                              header=None,
                              index_col=0)
  copyprobsFE = pd.read_csv("temp.FarmerEarly."+str(i)+".master_all_copyprobsperlocus.txt",
                              header=None,
                              index_col=0)
  copyprobsFL = pd.read_csv("temp.FarmerLate."+str(i)+".master_all_copyprobsperlocus.txt",
                              header=None,
                              index_col=0)
  copyprobsFM = pd.read_csv("temp.FarmerMiddle."+str(i)+".master_all_copyprobsperlocus.txt",
                              header=None,
                              index_col=0)
  copyprobsFarmer = copyprobsFA[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric) + \
                              copyprobsFE[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric) + \
                              copyprobsFL[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric) + \
                              copyprobsFM[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric)
  copyprobsFarmer.to_csv("temp.Farmer."+str(i)+".master_all_copyprobsperlocus.txt", header=False)
  for anc in ["CHG", "EHG", "Farmer", "WHG", "Yamnaya"]:
    print("Processing " + str(anc) + " chromosome " + str(i))
    anc_copyprobs = pd.read_csv("temp." + str(anc) + "." + str(i) + ".master_all_copyprobsperlocus.txt", header=None, index_col=0)
    if anc != 'Farmer':
      anc_copyprobs = anc_copyprobs[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric)
    positions = pd.read_csv("$phasefile_directory" + str(i) + ".merged.phase",
                            skiprows=2,
                            nrows=1,
                            sep=" ",
                            header=None).T.drop(0)
    anc_copyprobs.columns = positions[0].tolist()
    anc_copyprobs.to_csv("temp." + str(anc) + "." + str(i) + ".master_all_copyprobsperlocus.txt", sep=" ")
    print(str(anc) + " done!")
END
echo "*** Success! Now running PRS calculator: ***"
echo "python3 PRS_calculator.py -phasefile_directory $phasefile_directory -idfile $idfile -file_name $NEW_EXT"
`python3 PRS_calculator.py -phasefile_directory $phasefile_directory -idfile $idfile -file_name $NEW_EXT`
echo "*** All done, now tidying up temp files ***"
#rm $NEW_EXT
#rm temp*