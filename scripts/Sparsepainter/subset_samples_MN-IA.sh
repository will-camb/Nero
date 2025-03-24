#!/bin/bash
module load perl/5.38.0
module load gsl/2.5
module load bcftools/1.20
module load qctool/2.2.0
module load plink/2.0.0

# IA and MN
# Samples come from same directory so no need to merge
ref_dir="/maps/projects/lundbeck/scratch/vrb229/Project_IronAge/data/230905_impute_iadk/merge/vcf/"
ref_samples="mesoneo_samples.txt"
refname2id="mesoneo_name2id"
target_samples="ref_samples.txt"
targetname2id="refname2id"
output_dir="/datasets/ukb-AUDIT/iron_age/continental_ref_April24/genetic_data"

mkdir -p $output_dir

# Re-header using ref_name2id, select samples, save
# Re-header using mesoneo_name2id, select samples, save
process_chromosome() {
    chr=$1
    ref="$chr.230905_impute_iadk.glimpse.map_1240_info0.5.vcf.gz"
    bcftools reheader -s $refname2id -o $output_dir/"$ref" $ref_dir/"$ref"
    /projects/lundbeck/people/jdn321/software/pbwt/./pbwt -readVcfGT $output_dir/"$ref" -selectSamples $ref_samples -writeVcf $output_dir/ref."$chr".vcf
    bgzip $output_dir/ref."$chr".vcf
    tabix $output_dir/ref."$chr".vcf.gz
    rm $output_dir/"$ref"
    bcftools reheader -s $targetname2id -o $output_dir/"$ref" $ref_dir/"$ref"
    /projects/lundbeck/people/jdn321/software/pbwt/./pbwt -readVcfGT $output_dir/"$ref" -selectSamples $target_samples -writeVcf $output_dir/target."$chr".vcf
    bgzip $output_dir/target."$chr".vcf
    tabix $output_dir/target."$chr".vcf.gz
    rm $output_dir/"$ref"
}

export -f process_chromosome
export output_dir ref_dir target_dir refname2id targetname2id target_samples ref_samples # Update this if needed

# Use GNU Parallel to run the function for each chromosome
parallel process_chromosome ::: {1..22}


#UK Biobank
#target_dir="/datasets/ukb-AUDIT/haplotype_vcf"
#target_sample_file="/datasets/ukb-AUDIT/ukb58935_imp_chr1_v3_s487296.sample"
#keep_samples="/datasets/ukb-AUDIT/nonbritish/merged_vcfs/keep_samples"
#    target="$target_dir/ukb_imp_chr${chr}_v3.bgen"
#    /projects/lundbeck/people/jdn321/software/pbwt/./pbwt -readVcfGT $target -readSamples $target_sample_file -selectSamples $target_sample_file -writeVcf $output_dir/target.filtered.$chr.vcf.gz
#    #    willerslev/software/pbwt/pbwt -readVcfGT /willerslev/ukbiobank/haplotype_vcf/haplotype_chr10.vcf.gz -writeVcf /willerslev/ukbiobank/nonbritish/merged_vcfs/ukbb.haplotype_chr10.GT.vcf
#    /willerslev/software/pbwt/pbwt -readVcfGT /willerslev/ukbiobank/haplotype_vcf/haplotype_chr10.vcf.gz -writeVcf /willerslev/ukbiobank/nonbritish/merged_vcfs/ukbb.haplotype_chr10.GT.vcf
#    bcftools annotate ukbb.haplotype_chr1.GT.vcf.gz --rename-chrs /willerslev/ukbiobank/merged_UKBB_ref/pbwt/merged_vcfs/chr_rename_maps/1.map -o ukbb.haplotype_chr1.GT.chrinfo.vcf
#    bcftools reheader --samples UKBB_samples ukbb.haplotype_chr1.GT.chrinfo.vcf.gz -o ukbb.haplotype_chr1.GT.chrinfo.samples.vcf.gz