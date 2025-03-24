#!/bin/bash
module load perl/5.38.0
module load gsl/2.5
module load bcftools/1.20
module load qctool/2.2.0
module load plink/2.0.0

genetic_data_dir="/maps/projects/lundbeck/scratch/vrb229/Project_IronAge/data/240627_impute_iadk/merge/vcf/"
ref_samples="ref_samples.txt"
target_samples="target_samples.txt"
output_dir="genetic_data"

mkdir -p $output_dir

process_chromosome() {
    chr=$1
    file="$chr.240627_impute_iadk.glimpse.map_1240_info0.5.vcf.gz"
    /projects/lundbeck/people/jdn321/software/pbwt/./pbwt -readVcfGT "$genetic_data_dir/$file" -selectSamples $target_samples -writeVcf $output_dir/target."$chr".vcf
    /projects/lundbeck/people/jdn321/software/pbwt/./pbwt -readVcfGT "$genetic_data_dir/$file" -selectSamples $ref_samples -writeVcf $output_dir/ref."$chr".vcf
    bgzip $output_dir/target."$chr".vcf
    bgzip $output_dir/ref."$chr".vcf
    tabix $output_dir/target."$chr".vcf.gz
    tabix $output_dir/ref."$chr".vcf.gz
}

export -f process_chromosome
export output_dir genetic_data_dir target_samples ref_samples

# Use GNU Parallel to run the function for each chromosome
parallel process_chromosome ::: {1..22}

bcftools query -l $output_dir/target.22.vcf.gz > target_samples_ordered