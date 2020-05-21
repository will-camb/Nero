#!/bin/bash
## Copyright William Barrie, 2020
## All rights reserved 
## Example usage 'bash vcf_to_cp.sh mesoneo 3 22'

if [ "$#" -ne "3" ] ; then
    echo "Usage: vcf_to_cp.sh <dir> <firstchr> <lastchr>"
    echo "	<dir>: The location for the output"
    echo "	<firstchr>: The first chromosome in the sequence to be analysed"
    echo "	<lastchr>: The last chromosome in the sequence to be analysed i.e. we analyse all chromosomes between firstchr and lastchr inclusive"
    exit 0
fi
set -e

verbose=TRUE
vcfdir="/willerslev/scratch/neo/impute_20200221/03022020"
genmapdir="/willerslev/datasets/hapmapRecomb/2011-01_phaseII_B37"
dir="$1"
chrlist=`seq $2 $3`
plink=/willerslev/software/plink1.9/plink
plink2chromopainter=/willerslev/software/fs_4.0.1/plink2chromopainter.pl
convertrecfile=/willerslev/software/fs_4.0.1/convertrecfile.pl

echo "These are the chromosomes you've specified: $chrlist"
###################

# Get the list of sites that we have
echo "Extracting site list of the panel data into $dir/sitelist" 
mkdir -p $dir/sitelist
for chr in $chrlist; do
echo "Processing Chromosome $chr"
zcat ${vcfdir}/1000G.chr${chr}.maf1pm.best.03022020.vcf.gz | grep -v "^#" | cut -f3 | sed 's/:/_/' > $dir/sitelist/sitelist.chr${chr}.txt
echo "Done extracting site list!"
done

## Convert vcf to ped/map/fam
echo "Converting vcf files to plink format using plink1.9"
mkdir -p $dir/plinkformat
for chr in $chrlist; do
echo "Processing Chromosome $chr"
$plink --vcf ${vcfdir}/1000G.chr$chr.maf1pm.best.03022020.vcf.gz --recode12 --double-id --out $dir/plinkformat/processed.chr${chr}
done
echo "Done converting vcf to plink for all chromosomes!"

## Convert ped/map to cp input files
mkdir -p $dir/cpinput
echo "Converting plink files to cp phase file"
for chr in $chrlist; do
echo "Processing chromosome $chr"
$plink2chromopainter -p $dir/plinkformat/processed.chr${chr}.ped -m $dir/plinkformat/processed.chr${chr}.map -o $dir/cpinput/phasefile.chr$chr
done
echo "Done converting plink to cp input for all chromosomes!"

## Make cp recombfile using human recombination map 
echo "Making recomb file using provided genetic maps in directory $genmapdir"
for chr in $chrlist; do
echo "Processing chromosome $chr"
echo "Making recomb file recomb.chr${chr}"
$convertrecfile -M hapmap -U cM/MB $dir/cpinput/phasefile.chr$chr $genmapdir/genetic_map_GRCh37_chr${chr}.txt $dir/cpinput/recomb.chr${chr}
done
echo "Done making all recomb files!"

## Make cp idfile using 
echo "Making ids file $dir/cpinput/pop_ids"
for chr in $chrlist; do
echo "Processing Chromosome $chr"
bcftools query -l ${vcfdir}/1000G.chr${chr}.maf1pm.best.03022020.vcf.gz > $dir/cpinput/pop_ids
done 
echo "Done making cp input pop_ids file!"
#############################################################
