import pandas as pd 
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument("-copyprobsperlocus_location",
                    help="location of indivual copyprobsperlocus.out output from painting",
                    required=True)
args = parser.parse_args()

print('This is your specified location for copyprobsperlocus.out: ' + args.copyprobsperlocus_location)

copyprobsDF = pd.read_csv(args.copyprobsperlocus_location, delim_whitespace=True)

#Add haplotype and individual ID column
Hap_start_sites = copyprobsDF.loc[df.pos == 'HAP'].index.tolist()
copyprobsDF['ID'] = ''
copyprobsDF.iloc[:Hap_start_sites[1]].ID = str(copyprobsDF.iloc[Hap_start_sites[0]][0]) + '_' + str(copyprobsDF.iloc[Hap_start_sites[0]][1]) + '_' + str(copyprobsDF.iloc[Hap_start_sites[0]][2])
copyprobsDF.iloc[Hap_start_sites[1]:].ID = str(copyprobsDF.iloc[Hap_start_sites[1]][0]) + '_' + str(copyprobsDF.iloc[Hap_start_sites[1]][1])+ '_' + str(copyprobsDF.iloc[Hap_start_sites[1]][2])

#Eliminate rows containing Hap number
rowstoignore = list(copyprobsDF[copyprobsDF['pos'].str.contains("HAP") == True].index)
introws = list()
for i in list(copyprobsDF.index):
    if i not in rowstoignore:
        introws.append(i)
copyprobs_modDF = copyprobsDF.iloc[introws]

#Flip output to have first haplotype at the top
copyprobsDF = copyprobsDF[::-1].reset_index(drop=True)

#Save output
copyprobsDF.to_csv(args.copyprobsperlocus_location)
