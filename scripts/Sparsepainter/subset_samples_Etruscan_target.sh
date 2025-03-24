#!/bin/bash
module load perl/5.38.0
module load gsl/2.5
module load bcftools/1.20
module load qctool/2.2.0
module load plink/2.0.0

# Etruscans other ancients
ref_dir="/maps/projects/lundbeck/scratch/vrb229/project_etrus/data/240123_etrus/merge/vcf"
ref_samples="ref_samples.txt"
target_samples="target_samples.txt" # Etruscan and Roman samples
output_dir="genetic_data"
refname2id="refname2id"

mkdir -p $output_dir

process_chromosome() {
    chr=$1
    ref="$chr.240123_etrus.glimpse.map_1240_info0.5.vcf.gz"
    bcftools reheader -s $refname2id -o $output_dir/"$ref" $ref_dir/"$ref"
    ref=$output_dir/$ref
    tabix "$ref"
    /projects/lundbeck/people/jdn321/software/pbwt/./pbwt -readVcfGT "$ref" -selectSamples $target_samples -writeVcf $output_dir/target."$chr".vcf
    /projects/lundbeck/people/jdn321/software/pbwt/./pbwt -readVcfGT "$ref" -selectSamples $ref_samples -writeVcf $output_dir/ref."$chr".vcf
    bgzip $output_dir/target."$chr".vcf
    bgzip $output_dir/ref."$chr".vcf
    tabix $output_dir/target."$chr".vcf.gz
    tabix $output_dir/ref."$chr".vcf.gz
    rm "$ref"
}

export -f process_chromosome
export output_dir ref_dir target_dir refname2id keep_samples target_samples ref_samples

# Use GNU Parallel to run the function for each chromosome
parallel process_chromosome ::: {1..22}

bcftools query -l $output_dir/target.22.vcf.gz > target_samples_ordered