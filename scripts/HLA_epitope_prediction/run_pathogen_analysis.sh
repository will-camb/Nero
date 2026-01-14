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

# Run analyses for multiple pathogens using RefSeq reference sequences
# Syntax: run_analysis "NCBI query" "output_dir" "name" [n_jobs]
#
# Query format: Uses RefSeq[Filter] to get NCBI-curated reference genomes
# This ensures high-quality, non-redundant sequences without strain name issues
#
# Naming convention: Use common disease/pathogen names for easy interpretation
# Format: lowercase with underscores (e.g., hepatitis_B, strep_pyogenes)

echo "=========================================="
echo "VIRAL PATHOGENS"
echo "=========================================="

# Poxviridae
run_analysis "Variola virus[Organism] AND RefSeq[Filter]" \
             "./smallpox_analysis" \
             "smallpox" \
             50

# Paramyxoviridae
run_analysis "Measles morbillivirus[Organism] AND RefSeq[Filter]" \
             "./measles_analysis" \
             "measles" \
             50

run_analysis "Mumps orthorubulavirus[Organism] AND RefSeq[Filter]" \
             "./mumps_analysis" \
             "mumps" \
             50

# Togaviridae
run_analysis "Rubella virus[Organism] AND RefSeq[Filter]" \
             "./rubella_analysis" \
             "rubella" \
             50

# Retroviridae
run_analysis "Human immunodeficiency virus 1[Organism] AND RefSeq[Filter]" \
             "./HIV_analysis" \
             "HIV" \
             50

# Herpesviridae
run_analysis "Human gammaherpesvirus 4[Organism] AND RefSeq[Filter]" \
             "./EBV_analysis" \
             "EBV" \
             50

run_analysis "Human alphaherpesvirus 2[Organism] AND RefSeq[Filter]" \
             "./herpes_simplex_2_analysis" \
             "herpes_simplex_2" \
             50

run_analysis "Human alphaherpesvirus 3[Organism] AND RefSeq[Filter]" \
             "./varicella_zoster_analysis" \
             "varicella_zoster" \
             50

# Polyomaviridae
run_analysis "Merkel cell polyomavirus[Organism] AND RefSeq[Filter]" \
             "./merkel_polyomavirus_analysis" \
             "merkel_polyomavirus" \
             50

run_analysis "JC polyomavirus[Organism] AND RefSeq[Filter]" \
             "./JC_polyomavirus_analysis" \
             "JC_polyomavirus" \
             50

# Hepadnaviridae
run_analysis "Hepatitis B virus[Organism] AND RefSeq[Filter]" \
             "./hepatitis_B_analysis" \
             "hepatitis_B" \
             50

# Flaviviridae
run_analysis "Hepatitis C virus[Organism] AND RefSeq[Filter]" \
             "./hepatitis_C_analysis" \
             "hepatitis_C" \
             50

echo ""
echo "=========================================="
echo "BACTERIAL PATHOGENS"
echo "=========================================="

# Yersinia pestis - Plague
run_analysis "Yersinia pestis CO92[Organism] AND RefSeq[Filter]" \
             "./plague_analysis" \
             "plague" \
             30

# Mycobacterium tuberculosis - Tuberculosis
run_analysis "Mycobacterium tuberculosis H37Rv[Organism] AND RefSeq[Filter]" \
             "./tuberculosis_analysis" \
             "tuberculosis" \
             30

# Mycobacterium leprae - Leprosy
run_analysis "Mycobacterium leprae[Organism] AND RefSeq[Filter]" \
             "./leprosy_analysis" \
             "leprosy" \
             30

# Salmonella typhi - Typhoid fever
run_analysis "Salmonella enterica serovar Typhi CT18[Organism] AND RefSeq[Filter]" \
             "./typhoid_analysis" \
             "typhoid" \
             30

# Treponema pallidum - Syphilis
run_analysis "Treponema pallidum[Organism] AND RefSeq[Filter]" \
             "./syphilis_analysis" \
             "syphilis" \
             30

# Treponema pallidum pertenue - Yaws
run_analysis "Treponema pallidum pertenue[Organism] AND RefSeq[Filter]" \
             "./yaws_analysis" \
             "yaws" \
             30

# Streptococcus pyogenes - Group A Strep (scarlet fever, necrotizing fasciitis)
run_analysis "Streptococcus pyogenes[Organism] AND RefSeq[Filter]" \
             "./strep_pyogenes_analysis" \
             "strep_pyogenes" \
             30

# Streptococcus pneumoniae - Pneumococcus (pneumonia, meningitis)
run_analysis "Streptococcus pneumoniae[Organism] AND RefSeq[Filter]" \
             "./pneumococcus_analysis" \
             "pneumococcus" \
             30

# Borrelia recurrentis - Relapsing fever
run_analysis "Borrelia recurrentis[Organism] AND RefSeq[Filter]" \
             "./relapsing_fever_analysis" \
             "relapsing_fever" \
             30

# Brucella melitensis - Brucellosis
run_analysis "Brucella melitensis[Organism] AND RefSeq[Filter]" \
             "./brucellosis_analysis" \
             "brucellosis" \
             30

# Listeria monocytogenes - Listeriosis
run_analysis "Listeria monocytogenes[Organism] AND RefSeq[Filter]" \
             "./listeriosis_analysis" \
             "listeriosis" \
             30

# Leptospira interrogans - Leptospirosis
run_analysis "Leptospira interrogans Copenhageni[Organism] AND RefSeq[Filter]" \
             "./leptospirosis_analysis" \
             "leptospirosis" \
             20

# Vibrio cholerae - Cholera
run_analysis "Vibrio cholerae[Organism] AND RefSeq[Filter]" \
             "./cholera_analysis" \
             "cholera" \
             30

# Corynebacterium diphtheriae - Diphtheria
run_analysis "Corynebacterium diphtheriae[Organism] AND RefSeq[Filter]" \
             "./diphtheria_analysis" \
             "diphtheria" \
             30

# Bacillus anthracis - Anthrax
run_analysis "Bacillus anthracis[Organism] AND RefSeq[Filter]" \
             "./anthrax_analysis" \
             "anthrax" \
             30

# Rickettsia prowazekii - Epidemic typhus
run_analysis "Rickettsia prowazekii[Organism] AND RefSeq[Filter]" \
             "./typhus_analysis" \
             "typhus" \
             30

echo ""
echo "=========================================="
echo "PARASITIC PATHOGENS"
echo "=========================================="

# Plasmodium vivax - Malaria (tertian)
run_analysis "Plasmodium vivax PvP01[Organism] AND RefSeq[Filter]" \
             "./malaria_vivax_analysis" \
             "malaria_vivax" \
             20

# Plasmodium falciparum - Malaria (falciparum)
run_analysis "Plasmodium falciparum 3D7[Organism] AND RefSeq[Filter]" \
             "./malaria_falciparum_analysis" \
             "malaria_falciparum" \
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
