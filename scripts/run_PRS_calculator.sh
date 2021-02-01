#!/bin/bash
#cd $PBS_O_WORKDIR
module load anaconda3/4.4.0
#source /willerslev/software/venv_python3.6/bin/activate

if [ "$#" -ne "5" ] ; then
    echo "Usage: run_PRS_calculator.sh <manifest> <file_name> <copyprobs_directory> <phasefile_directory>"
    echo "<manifest>: the manifest csv file"
    echo "<phenotype_file>: file containing filename from manifest to be calculated, 1 per row eg E4_OBESITY.gwas.imputed_v3.both_sexes.tsv.bgz"
    echo "<copyprobs_directory>: Directory of all_copyprobsperlocus files; should be named in form n.master_all_copyprobsperlocus.txt.gz"
    echo "<phasefile_directory>: Directory of phasefiles; should be named in form chr#.merged.phase"
    echo "<idfile>: Path to idfile: ordered_all_pop_ids_mapped ('individualID popID 1')"
    exit 0
fi

manifest="$1"
phenotype_file="$2"
#file_name=$2
copyprobs_directory="$3"
phasefile_directory="$4"
idfile="$5"

#Loop through lines and download results
#while IFS= read -r line; do
#echo "Extracting results from Dropbox for $line"
#`awk -F ',' -v $line="$line" '$5 == line { print $6 }' $manifest`
#MY_FILE=$line
#NEW_EXT=${MY_FILE/bgz/gz}
#mv $line "$NEW_EXT"
#done < $phenotype_file

echo "*** Extracting results file from Dropbox ***"
#`awk -F ',' -v file_name="$file_name" '$5 == file_name { print $6 }' $manifest`
#MY_FILE=$file_name
#NEW_EXT=${MY_FILE/bgz/gz}
#mv $file_name "$NEW_EXT"
echo "*** Success! Now making temp copyprobs files per ancestry ***"
#chrlist=`seq 1 22`
chrlist=`seq 1 21`
#for chr in $chrlist; do
#zcat $3/$chr.master_all_copyprobsperlocus.txt.gz | awk 'NR % 10 == 1' - > temp.CHG.$chr.master_all_copyprobsperlocus.txt
#gzip temp.CHG.$chr.master_all_copyprobsperlocus.txt
#zcat $3/$chr.master_all_copyprobsperlocus.txt.gz | awk 'NR % 10 == 2' - > temp.EHG.$chr.master_all_copyprobsperlocus.txt
#gzip temp.EHG.$chr.master_all_copyprobsperlocus.txt
#zcat $3/$chr.master_all_copyprobsperlocus.txt.gz | awk 'NR % 10 == 3' - > temp.FarmerAnatolian.$chr.master_all_copyprobsperlocus.txt
#gzip temp.FarmerAnatolian.$chr.master_all_copyprobsperlocus.txt
#zcat $3/$chr.master_all_copyprobsperlocus.txt.gz | awk 'NR % 10 == 4' - > temp.FarmerEarly.$chr.master_all_copyprobsperlocus.txt
#gzip temp.FarmerEarly.$chr.master_all_copyprobsperlocus.txt
#zcat $3/$chr.master_all_copyprobsperlocus.txt.gz | awk 'NR % 10 == 5' - > temp.FarmerLate.$chr.master_all_copyprobsperlocus.txt
#gzip temp.FarmerLate.$chr.master_all_copyprobsperlocus.txt
#zcat $3/$chr.master_all_copyprobsperlocus.txt.gz | awk 'NR % 10 == 6' - > temp.FarmerMiddle.$chr.master_all_copyprobsperlocus.txt
#gzip temp.FarmerMiddle.$chr.master_all_copyprobsperlocus.txt
#zcat $3/$chr.master_all_copyprobsperlocus.txt.gz | awk 'NR % 10 == 7' - > temp.EastAsian.$chr.master_all_copyprobsperlocus.txt
#gzip temp.EastAsian.$chr.master_all_copyprobsperlocus.txt
#zcat $3/$chr.master_all_copyprobsperlocus.txt.gz | awk 'NR % 10 == 8' - > temp.African.$chr.master_all_copyprobsperlocus.txt
#gzip temp.African.$chr.master_all_copyprobsperlocus.txt
#zcat $3/$chr.master_all_copyprobsperlocus.txt.gz | awk 'NR % 10 == 9' - > temp.WHG.$chr.master_all_copyprobsperlocus.txt
#gzip temp.WHG.$chr.master_all_copyprobsperlocus.txt
#zcat $3/$chr.master_all_copyprobsperlocus.txt.gz | awk 'NR % 10 == 0' - > temp.Yamnaya.$chr.master_all_copyprobsperlocus.txt
#gzip temp.Yamnaya.$chr.master_all_copyprobsperlocus.txt
#done
echo "*** Success! Now making farmer aggregate, and reformatting temp files in Python ***"
python3 << END
import pandas as pd
for i in range(1,23):
  print("Reading Farmer chr" + str(i))
  copyprobsFA = pd.read_csv("temp.FarmerAnatolian."+str(i)+".master_all_copyprobsperlocus.txt.gz",
                              header=None,
                              index_col=0)
  copyprobsFE = pd.read_csv("temp.FarmerEarly."+str(i)+".master_all_copyprobsperlocus.txt.gz",
                              header=None,
                              index_col=0)
  copyprobsFL = pd.read_csv("temp.FarmerLate."+str(i)+".master_all_copyprobsperlocus.txt.gz",
                              header=None,
                              index_col=0)
  copyprobsFM = pd.read_csv("temp.FarmerMiddle."+str(i)+".master_all_copyprobsperlocus.txt.gz",
                              header=None,
                              index_col=0)
  copyprobsFarmer = copyprobsFA[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric) + \
                              copyprobsFE[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric) + \
                              copyprobsFL[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric) + \
                              copyprobsFM[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric)
  copyprobsFarmer.to_csv("temp.Farmer."+str(i)+".master_all_copyprobsperlocus.txt", header=False)
  print("Farmer done, now converting other ancestries to correct format")
  for anc in ["CHG", "EHG", "Farmer", "African", "EastAsian", "WHG", "Yamnaya"]:
    print("Processing " + str(anc) + " chromosome " + str(i))
    if anc == 'Farmer':
      anc_copyprobs = pd.read_csv("temp." + str(anc) + "." + str(i) + ".master_all_copyprobsperlocus.txt", header=None, index_col=0)
    else:
      anc_copyprobs = pd.read_csv("temp." + str(anc) + "." + str(i) + ".master_all_copyprobsperlocus.txt.gz", header=None, index_col=0)
      anc_copyprobs = anc_copyprobs[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric)
    positions = pd.read_csv("$phasefile_directory/" + str(i) + ".merged.phase",
                            skiprows=2,
                            nrows=1,
                            sep=" ",
                            header=None).T.drop(0)
    anc_copyprobs.columns = positions[0].tolist()
    anc_copyprobs.to_csv("temp." + str(anc) + "." + str(i) + ".master_all_copyprobsperlocus.txt", sep=" ")
    print(str(anc) + " done!")
END
echo "*** Success! Now gzipping temp files ***"
for chr in $chrlist; do
gzip --force temp.CHG.$chr.master_all_copyprobsperlocus.txt
gzip --force temp.EHG.$chr.master_all_copyprobsperlocus.txt
gzip --force temp.Farmer.$chr.master_all_copyprobsperlocus.txt
gzip --force temp.EastAsian.$chr.master_all_copyprobsperlocus.txt
gzip --force temp.African.$chr.master_all_copyprobsperlocus.txt
gzip --force temp.WHG.$chr.master_all_copyprobsperlocus.txt
gzip --force temp.Yamnaya.$chr.master_all_copyprobsperlocus.txt
done
echo "*** Success! Now running PRS calculator: ***"
#while IFS= read -r line; do
#echo "Writing PRS_calculator commands for $line"
##MY_FILE=$line
##NEW_EXT=${MY_FILE/bgz/gz}
#NEW_EXT=$line
#echo "python3 PRS_calculator.py -phasefile_directory $phasefile_directory -idfile $idfile -file_name $NEW_EXT" >> PRS_calculator_commands
#done < "$phenotype_file"

echo "Now running commands in PRS_calculator_commands in parallel"
#cat PRS_calculator_commands | parallel

#echo "python3 PRS_calculator.py -phasefile_directory $phasefile_directory -idfile $idfile -file_name $NEW_EXT"
#`python3 PRS_calculator.py -phasefile_directory $phasefile_directory -idfile $idfile -file_name $NEW_EXT`
echo "*** All done, now tidying up temp files ***"
#rm $NEW_EXT
#rm temp*
#rm PRS_calculator_commands