#!/bin/bash
module load perl/5.38.0
module load gsl/2.5
module load bcftools/1.20
module load qctool/2.2.0
module load plink/2.0.0

# 1000g and IA
ref_dir="/maps/projects/lundbeck/scratch/vrb229/Project_IronAge/data/230905_impute_iadk/merge/vcf/"
target_dir="/projects/lundbeck/data/1000genomes_2015_nature"
ref_samples="ref_samples.txt"  # Update
target_samples="1000g_samples"
output_dir="genetic_data"
refname2id="refname2id"
bcftools query -l $target_dir/22.1000g.freeze9.umich.GRCh37.snps.biallelic.pass.vcf.gz > 1000g_samples

mkdir -p $output_dir

process_chromosome() {
    chr=$1
    ref="$chr.230905_impute_iadk.glimpse.map_1240_info0.5.vcf.gz"
    target="$chr.1000g.freeze9.umich.GRCh37.snps.biallelic.pass.vcf.gz"
    bcftools reheader -s $refname2id -o $output_dir/"$ref" $ref_dir/"$ref"
    ref=$output_dir/$ref
    target=$target_dir/$target
    qctool -g "$ref" -g "$target" -og $output_dir/ref.target."$chr".vcf
    rm "$ref"
    bgzip $output_dir/ref.target."$chr".vcf
    tabix $output_dir/ref.target."$chr".vcf.gz
    /projects/lundbeck/people/jdn321/software/pbwt/./pbwt -readVcfGT $output_dir/ref.target."$chr".vcf.gz -selectSamples $target_samples -writeVcf $output_dir/target."$chr".vcf
    /projects/lundbeck/people/jdn321/software/pbwt/./pbwt -readVcfGT $output_dir/ref.target."$chr".vcf.gz -selectSamples $ref_samples -writeVcf $output_dir/ref."$chr".vcf
    bgzip $output_dir/target."$chr".vcf
    bgzip $output_dir/ref."$chr".vcf
    tabix $output_dir/target."$chr".vcf.gz
    tabix $output_dir/ref."$chr".vcf.gz
    rm $output_dir/ref.target."$chr".vcf.gz $output_dir/ref.target."$chr".vcf.gz.tbi
}

export -f process_chromosome
export output_dir ref_dir target_dir refname2id keep_samples target_samples ref_samples

# Use GNU Parallel to run the function for each chromosome
parallel process_chromosome ::: {1..22}