#!/bin/bash
module load perl/5.38.0
module load gsl/2.5
module load bcftools/1.20

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_vcf output_vcf"
    exit 1
fi

# input_vcf="/datasets/ukb-AUDIT/iron_age/continental_ref_April24/genetic_data/target.21.vcf"
# output_vcf="/datasets/ukb-AUDIT/iron_age/continental_ref_April24/genetic_data/test.target.21.vcf"
input_vcf="$1"
output_vcf="$2"
temp_dir=$(mktemp -d)
cd "$temp_dir" || exit

if [[ $input_vcf == *.vcf.gz ]]; then
    gunzip "$input_vcf"
fi

# Extract sample names
bcftools query -l "$input_vcf" > samples.txt

# Function to process each sample
process_sample() {
    sample=$1
    input_vcf=$2
    temp_dir=$3

    # Split each sample into two haplotypes
    bcftools view -s "$sample" "$input_vcf" | \
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT]\n' | \
    awk -v sample="$sample" '{split($6, g, "|"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t.\t.\t.\tGT\t"g[1]}' > "${temp_dir}"/"${sample}"_hap1.txt

    bcftools view -s "$sample" "$input_vcf" | \
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT]\n' | \
    awk -v sample="$sample" '{split($6, g, "|"); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t.\t.\t.\tGT\t"g[2]}' > "${temp_dir}"/"${sample}"_hap2.txt

    # Convert text files to VCF format and compress with bgzip
    (echo '##fileformat=VCFv4.2';
     echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">';
     echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'"${sample}"_hap1;
     cat "${temp_dir}"/"${sample}"_hap1.txt) | bgzip -c > "${temp_dir}"/"${sample}"_hap1.vcf.gz

    (echo '##fileformat=VCFv4.2';
     echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">';
     echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'"${sample}"_hap2;
     cat "${temp_dir}"/"${sample}"_hap2.txt) | bgzip -c > "${temp_dir}"/"${sample}"_hap2.vcf.gz

    # Index the compressed VCF files
    tabix -p vcf "${temp_dir}"/"${sample}"_hap1.vcf.gz
    tabix -p vcf "${temp_dir}"/"${sample}"_hap2.vcf.gz
}

export -f process_sample

# Use parallel to process all samples concurrently
parallel --jobs 8 process_sample {} "$input_vcf" "$temp_dir" :::: samples.txt

# Merge all haploid VCF files into one final VCF
hap_files=$(ls "${temp_dir}"/*_hap*.vcf.gz)
bcftools merge "$hap_files" -Oz -o "$output_vcf"
gzip "$output_vcf"

# Clean up
cd - || exit
rm -r "$temp_dir"

