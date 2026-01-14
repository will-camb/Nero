#!/bin/bash
# Example usage of the Pathogen HLA Class I & II Binding Analysis Pipeline

# Install dependencies
# echo "Installing Python dependencies..."
# pip install -r requirements.txt --break-system-packages

# Make the main script executable
chmod +x pathogen_hla_binding_class_i_ii.py

# Example 1: Analyze smallpox (Variola virus) - both Class I and II
echo "Running analysis for Variola virus (smallpox)..."
python pathogen_hla_binding_class_i_ii.py \
    --email "wb275@cam.ac.uk" \
    --organism "Variola virus[Organism] AND RefSeq[Filter]" \
    --output-dir "./smallpox_analysis" \
    --name "smallpox" \
    --n-jobs 80

echo "Running analysis for Yersinia pestis CO92 (plague)..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Yersinia pestis CO92[Organism] AND RefSeq[Filter]" \
   --output-dir "plague_analysis" \
   --name "plague" \
   --n-jobs 80

echo "Running analysis for HIV..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "HIV[Organism] AND RefSeq[Filter]" \
   --output-dir "HIV_analysis" \
   --name "HIV" \
   --n-jobs 80

echo "Running analysis for measles..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "measles[Organism] AND RefSeq[Filter]" \
   --output-dir "measles_analysis" \
   --name "measles" \
   --n-jobs 80

echo "Running analysis for leprosy (Mycobacterium leprae)..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Mycobacterium leprae[Organism] AND RefSeq[Filter]" \
   --output-dir "leprosy_analysis" \
   --name "leprosy" \
   --n-jobs 80

echo "Running analysis for Borrelia recurrentis..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Borrelia recurrentis[Organism] AND RefSeq[Filter]" \
   --output-dir "B_recurrentis_analysis" \
   --name "B_recurrentis" \
   --n-jobs 80

echo "Running analysis for Leptospira interrogans..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Leptospira interrogans Copenhageni[Organism] AND RefSeq[Filter]" \
   --output-dir "L_interrogans_analysis" \
   --name "L_interrogans" \
   --n-jobs 30

echo "Running analysis for Treponema pallidum..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Treponema pallidum[Organism] AND RefSeq[Filter]" \
   --output-dir "T_pallidum_analysis" \
   --name "T_pallidum" \
   --n-jobs 40

echo "Running analysis for Mycobacterium tuberculosis..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Mycobacterium tuberculosis H37Rv[Organism] AND RefSeq[Filter]" \
   --output-dir "M_tuberculosis_analysis" \
   --name "M_tuberculosis" \
   --n-jobs 30

echo "Running analysis for Plasmodium vivax..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Plasmodium vivax[Organism] AND RefSeq[Filter]" \
   --output-dir "P_vivax_analysis" \
   --name "P_vivax" \
   --n-jobs 30

echo "Running analysis for EBV..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Epstein–Barr virus[Organism] AND RefSeq[Filter]" \
   --output-dir "EBV_analysis" \
   --name "EBV" \
   --n-jobs 50

echo "Running analysis for Salmonella typhi CT18..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Salmonella typhi CT18[Organism] AND RefSeq[Filter]" \
   --output-dir "Salmonella_typhi_analysis" \
   --name "Salmonella_typhi" \
   --n-jobs 50

echo "Running analysis for HSV-2..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Herpes simplex virus 2[Organism] AND RefSeq[Filter]" \
   --output-dir "HSV-2_analysis" \
   --name "HSV-2" \
   --n-jobs 50

echo "Running analysis for Merkel cell polyomavirus..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Merkel cell polyomavirus[Organism] AND RefSeq[Filter]" \
   --output-dir "MCPyV_analysis" \
   --name "MCPyV" \
   --n-jobs 50

echo "Running analysis for JC virus (Human polyomavirus 2)..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Human polyomavirus 2[Organism] AND RefSeq[Filter]" \
   --output-dir "JCV_analysis" \
   --name "JCV" \
   --n-jobs 50

echo "Running analysis for Varicella zoster virus..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Varicella zoster virus[Organism] AND RefSeq[Filter]" \
   --output-dir "VZV_analysis" \
   --name "VZV" \
   --n-jobs 50

echo "Running analysis for Mumps virus..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Mumps virus[Organism] AND RefSeq[Filter]" \
   --output-dir "Mumps_analysis" \
   --name "Mumps" \
   --n-jobs 50

echo "Running analysis for Rubella virus..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Rubella virus[Organism] AND RefSeq[Filter]" \
   --output-dir "Rubella_analysis" \
   --name "Rubella" \
   --n-jobs 50

echo "Running analysis for Plasmodium falciparum..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Plasmodium falciparum 3D7[Organism] AND RefSeq[Filter]" \
   --output-dir "P_falciparum_analysis" \
   --name "P_falciparum" \
   --n-jobs 40

echo "Running analysis for Streptococcus pyogenes..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Streptococcus pyogenes[Organism] AND RefSeq[Filter]" \
   --output-dir "S_pyogenes_analysis" \
   --name "S_pyogenes" \
   --n-jobs 50

echo "Running analysis for Hepatitis B virus..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Hepatitis B virus[Organism] AND RefSeq[Filter]" \
   --output-dir "HBV_analysis" \
   --name "HBV" \
   --n-jobs 50

echo "Running analysis for Hepatitis C virus..."
python3 pathogen_hla_binding_class_i_ii.py \
   --email "wb275@cam.ac.uk" \
   --organism "Hepatitis C virus[Organism] AND RefSeq[Filter]" \
   --output-dir "HCV_analysis" \
   --name "HCV" \
   --n-jobs 50

•••••••••••••••••••••••••••••••••••••••••••••••••


dandycomp05fl: falciparum

echo "Analysis pipeline setup complete!"
echo ""
echo "IMPORTANT NOTES:"
echo "1. Replace 'your.email@example.com' with your actual email address"
echo "2. Ensure SignalP 6.0+ is installed and in your PATH"
echo "3. Verify paths are correct:"
echo "   - netMHCIIpan: /home/jdn321/software/netMHCIIpan-4.3/netMHCIIpan"
echo "   - netMHCpan: /home/jdn321/software/netMHCpan-4.2/netMHCpan"
echo "4. The analysis may take considerable time depending on proteome size"
echo "5. Class I analysis generates 8-11mer peptides, Class II generates 15mers"
echo ""
echo "Output files:"
echo "- organism_proteome.fasta: Full proteome from NCBI"
echo "- organism_secreted.fasta: Filtered secreted proteins"
echo "- signalp/: SignalP analysis results"
echo "- netmhciipan_class_ii/: netMHCIIpan prediction files"
echo "- netmhcpan_class_i/: netMHCpan prediction files"
echo "- organism_hla_binding_results.csv: Final combined results table"
echo ""
echo "Results file columns:"
echo "- protein_id: Source protein identifier"
echo "- peptide_sequence: Peptide sequence"
echo "- peptide_length: Length of peptide (8-11 for Class I, 15 for Class II)"
echo "- position: Position in source protein"
echo "- hla_class: 'I' or 'II'"
echo "- hla_allele: HLA allele name"
echo "- Score_EL: Binding score"
echo "- %Rank_EL: Percentile rank (lower = stronger binding)"
