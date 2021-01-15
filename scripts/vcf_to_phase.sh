#!/bin/bash
cd $PBS_O_WORKDIR
module load pbwt/20200610
module load tools finestructure/4.1.1
module load perl/5.20.1

if [ "$#" -ne "2" ] ; then
    echo "Usage: vcf_to_phase.sh <Individual ID> <Chromosome number>"
    echo "<vcf file>: The name of the vcf file; should be without .vcf.gz"
    echo "<Chromosome number>: The chromosome to be analysed"
    exit 0
fi

vcf="$1"
chr="$2"
pbwt -readVcfGT haps.$chr/$vcf.vcf.gz -writeImputeHapsG haps.$chr/$vcf.haps
impute2chromopainter.pl haps.$chr/$vcf.haps phase.$chr/$vcf.phase
rm haps.$chr/$vcf.haps
tail -n 2 phase.$chr/$vcf.phase > phase.$chr/$vcf.reduced.phase
rm phase.$chr/$vcf.phase