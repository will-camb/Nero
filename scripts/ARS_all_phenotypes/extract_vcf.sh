#!/bin/bash
conda activate PRS # run 'source' not 'bash' when executing script
module load bgen/1.1.4
module load perl/5.38.0
module load gsl/2.5
module load bcftools/1.16

if [ "$#" -ne "1" ] ; then
    echo "Usage: source extract_vcf.sh <SNP_file>"
    echo "<SNP_file>: A space separated file with SNP (rsID) and CHR columns"
    exit 0
fi

SNP_file="$1"
bgen_dir="/maps/datasets/ukb-AUDIT/imputation_bgen"

header=$(head -n 1 "$SNP_file")
snp_col=$(echo "$header" | tr ' ' '\n' | grep -n -w "SNP" | cut -d: -f1)
chr_col=$(echo "$header" | tr ' ' '\n' | grep -n -w "CHR" | cut -d: -f1)

# Process each chromosome separately
for chr in {1..22}; do
    echo "Processing chromosome $chr..."

    # Extract rsIDs for the current chromosome
    rsIDs=$(awk -v chr_col=$chr_col -v snp_col=$snp_col -v chr="$chr" 'NR > 1 && $chr_col == chr {print $snp_col}' "$SNP_file")

    # Check if there are any rsIDs for this chromosome
    if [ -z "$rsIDs" ]; then
        echo "No rsIDs found for chromosome $chr, skipping..."
        continue
    fi
    rsid_file="chr${chr}_rsids.txt"
    echo "$rsIDs" > "$rsid_file"

    # Define the input BGEN file for the current chromosome
    bgen_file="$bgen_dir/ukb_imp_chr${chr}_v3.bgen"

    # Define intermediate and final output files
    intermediate_vcf="chr${chr}.vcf"
    gzipped_vcf="${intermediate_vcf}.gz"
    final_vcf="chr${chr}.GT.vcf"
    final_gzipped_vcf="${final_vcf}.gz"

    # Extract SNPs from BGEN file to VCF format
    bgenix -g "$bgen_file" -incl-rsids "$rsid_file" -vcf > "$intermediate_vcf"

    # Compress and index the intermediate VCF file
    bgzip "$intermediate_vcf"
    tabix "$gzipped_vcf"

    # Use pbwt to read and write VCF files
    pbwt -readVcfGT "$gzipped_vcf" -writeVcf "$final_vcf"

    # Compress and index the final VCF file
    bgzip "$final_vcf"
    tabix "$final_gzipped_vcf"
    rm "$rsid_file"
done

# Combine all the final gzipped VCF files into one
echo "Combining all chromosomes into one VCF file..."
bcftools concat -Oz -o "combined.vcf.gz" chr*.GT.vcf.gz
tabix "combined.vcf.gz"

rm -f chr*.vcf.gz chr*.vcf.gz.tbi chr*.GT.vcf.gz chr*.GT.vcf.gz.tbi

