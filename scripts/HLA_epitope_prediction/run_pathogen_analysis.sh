#!/bin/bash
# Batch runner for HLA epitope prediction across multiple pathogens
# Uses improved incremental processing pipeline

set -e  # Exit on error

# Configuration
EMAIL="wb275@cam.ac.uk"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${SCRIPT_DIR}/hla_epitope_predictor.py"

# Tool paths (modify if needed)
NETMHCPAN_PATH="/home/jdn321/software/netMHCpan-4.2/netMHCpan"
NETMHCIIPAN_PATH="/home/jdn321/software/netMHCIIpan-4.3/netMHCIIpan"
SIGNALP_MODELS="/home/jdn321/software/signalp6_fast/signalp-6-package/models/"

# Allele files (use defaults if not specified)
CLASS_I_ALLELES="${SCRIPT_DIR}/alleles_class_i_default.txt"
CLASS_II_ALLELES="${SCRIPT_DIR}/alleles_class_ii_default.txt"

# Number of parallel jobs (reduced from 80 to prevent memory issues)
# For large proteomes (>3000 proteins), use 20-30
# For small proteomes (<500 proteins), can use 40-60
N_JOBS=30

# Peptide generation parameters (conservative defaults for ~10x speedup)
# Class I: 9-mers only (instead of 8,9,10,11)
# Step size: 3 (instead of 1 = fully overlapping)
CLASS_I_LENGTHS="9"
CLASS_I_STEP=3
CLASS_II_LENGTH=15
CLASS_II_STEP=3

echo "============================================"
echo "HLA Epitope Prediction Pipeline"
echo "============================================"
echo "Email: $EMAIL"
echo "Parallel jobs: $N_JOBS"
echo "Class I: lengths=$CLASS_I_LENGTHS, step=$CLASS_I_STEP"
echo "Class II: length=$CLASS_II_LENGTH, step=$CLASS_II_STEP"
echo "Class I alleles: $CLASS_I_ALLELES"
echo "Class II alleles: $CLASS_II_ALLELES"
echo ""

# Function to run analysis for one organism
run_analysis() {
    local organism_query="$1"
    local output_dir="$2"
    local organism_name="$3"
    local n_jobs="${4:-$N_JOBS}"  # Use default if not provided

    echo ""
    echo "========================================"
    echo "Analyzing: $organism_name"
    echo "Output: $output_dir"
    echo "Parallel jobs: $n_jobs"
    echo "========================================"

    python3 "$SCRIPT" \
        --email "$EMAIL" \
        --organism "$organism_query" \
        --output-dir "$output_dir" \
        --name "$organism_name" \
        --class-i-alleles "$CLASS_I_ALLELES" \
        --class-ii-alleles "$CLASS_II_ALLELES" \
        --class-i-lengths $CLASS_I_LENGTHS \
        --class-i-step "$CLASS_I_STEP" \
        --class-ii-length "$CLASS_II_LENGTH" \
        --class-ii-step "$CLASS_II_STEP" \
        --netmhcpan-path "$NETMHCPAN_PATH" \
        --netmhciipan-path "$NETMHCIIPAN_PATH" \
        --signalp-model-dir "$SIGNALP_MODELS" \
        --n-jobs "$n_jobs"

    if [ $? -eq 0 ]; then
        echo "✓ $organism_name completed successfully"
    else
        echo "✗ $organism_name failed"
        return 1
    fi
}

# Run analyses for multiple pathogens using REFERENCE STRAINS
# Syntax: run_analysis "NCBI query" "output_dir" "name" [n_jobs]

echo "=========================================="
echo "VIRAL PATHOGENS"
echo "=========================================="

# DNA viruses - Poxviridae
run_analysis "Variola virus strain Bangladesh-1975[Organism]" \
             "./smallpox_analysis" \
             "smallpox" \
             40

# RNA viruses - Paramyxoviridae
run_analysis "Measles morbillivirus strain Edmonston[Organism]" \
             "./measles_analysis" \
             "measles" \
             40

run_analysis "Mumps orthorubulavirus strain Enders[Organism]" \
             "./mumps_analysis" \
             "mumps" \
             40

# RNA viruses - Togaviridae
run_analysis "Rubella virus strain RA27/3[Organism]" \
             "./rubella_analysis" \
             "rubella" \
             40

# RNA viruses - Retroviridae
run_analysis "Human immunodeficiency virus 1 HXB2[Organism]" \
             "./HIV_analysis" \
             "HIV" \
             40

# DNA viruses - Herpesviridae
run_analysis "Human gammaherpesvirus 4 strain B95-8[Organism]" \
             "./EBV_analysis" \
             "EBV" \
             40

run_analysis "Human alphaherpesvirus 2 strain HG52[Organism]" \
             "./HSV2_analysis" \
             "HSV2" \
             40

run_analysis "Human alphaherpesvirus 3 strain Dumas[Organism]" \
             "./VZV_analysis" \
             "VZV" \
             40

# DNA viruses - Polyomaviridae
run_analysis "Merkel cell polyomavirus strain R17b[Organism]" \
             "./MCPyV_analysis" \
             "MCPyV" \
             40

run_analysis "JC polyomavirus strain Mad1[Organism]" \
             "./JCV_analysis" \
             "JCV" \
             40

echo ""
echo "=========================================="
echo "BACTERIAL PATHOGENS"
echo "=========================================="

# Plague
run_analysis "Yersinia pestis CO92[Organism]" \
             "./plague_analysis" \
             "plague" \
             30

# Tuberculosis
run_analysis "Mycobacterium tuberculosis H37Rv[Organism]" \
             "./M_tuberculosis_analysis" \
             "M_tuberculosis" \
             30

# Leprosy
run_analysis "Mycobacterium leprae TN[Organism]" \
             "./leprosy_analysis" \
             "leprosy" \
             30

# Typhoid
run_analysis "Salmonella enterica subsp. enterica serovar Typhi str. CT18[Organism]" \
             "./Salmonella_typhi_analysis" \
             "Salmonella_typhi" \
             30

# Syphilis
run_analysis "Treponema pallidum subsp. pallidum str. Nichols[Organism]" \
             "./T_pallidum_analysis" \
             "T_pallidum" \
             30

# Streptococcus pyogenes
run_analysis "Streptococcus pyogenes M1 GAS[Organism]" \
             "./Streptococcus_pyogenes_analysis" \
             "Streptococcus_pyogenes" \
             30

# Streptococcus pneumoniae
run_analysis "Streptococcus pneumoniae TIGR4[Organism]" \
             "./Streptococcus_pneumoniae_analysis" \
             "Streptococcus_pneumoniae" \
             30

# Borrelia recurrentis
run_analysis "Borrelia recurrentis A1[Organism]" \
             "./B_recurrentis_analysis" \
             "B_recurrentis" \
             30

# Brucella melitensis
run_analysis "Brucella melitensis 16M[Organism]" \
             "./Brucella_melitensis_analysis" \
             "Brucella_melitensis" \
             30

# Listeria monocytogenes
run_analysis "Listeria monocytogenes EGD-e[Organism]" \
             "./Listeria_monocytogenes_analysis" \
             "Listeria_monocytogenes" \
             30

# Leptospira interrogans
run_analysis "Leptospira interrogans serovar Copenhageni str. Fiocruz L1-130[Organism]" \
             "./L_interrogans_analysis" \
             "L_interrogans" \
             20

echo ""
echo "=========================================="
echo "PARASITIC PATHOGENS"
echo "=========================================="

# Plasmodium vivax
run_analysis "Plasmodium vivax Sal-1[Organism]" \
             "./P_vivax_analysis" \
             "P_vivax" \
             20

# Plasmodium falciparum
run_analysis "Plasmodium falciparum 3D7[Organism]" \
             "./P_falciparum_analysis" \
             "P_falciparum" \
             20

echo ""
echo "============================================"
echo "All analyses complete!"
echo "============================================"
echo ""
echo "Next steps:"
echo "1. Review results in each *_analysis/ directory"
echo "2. Check *_hla_predictions.csv for binding predictions"
echo "3. Run downstream analysis (e.g., aggregate by allele)"
echo ""
