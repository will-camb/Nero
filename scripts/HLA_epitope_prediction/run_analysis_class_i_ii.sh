#!/bin/bash
# Example usage of the Pathogen HLA Class I & II Binding Analysis Pipeline

# Install dependencies
echo "Installing Python dependencies..."
pip install -r requirements.txt --break-system-packages

# Make the main script executable
chmod +x pathogen_hla_binding.py

# Example 1: Analyze smallpox (Variola virus) - both Class I and II
echo "Running analysis for Variola virus (smallpox) - Class I & II..."
python pathogen_hla_binding.py \
    --email "your.email@example.com" \
    --organism "Variola virus[Organism] AND RefSeq[Filter]" \
    --output-dir "./smallpox_analysis" \
    --name "smallpox" \
    --n-jobs 20

# Example 2: Analyze another pathogen (e.g., Monkeypox)
echo "Example command for Monkeypox virus..."
echo "python pathogen_hla_binding.py \\"
echo "    --email \"your.email@example.com\" \\"
echo "    --organism \"Monkeypox virus[Organism] AND RefSeq[Filter]\" \\"
echo "    --output-dir \"./monkeypox_analysis\" \\"
echo "    --name \"monkeypox\" \\"
echo "    --n-jobs 20"

# Example 3: Analyze specific bacterial strain
echo "Example command for Y. pestis specific strain..."
echo "python pathogen_hla_binding.py \\"
echo "    --email \"your.email@example.com\" \\"
echo "    --organism \"Yersinia pestis CO92[Organism] AND RefSeq[Filter]\" \\"
echo "    --output-dir \"./ypestis_analysis\" \\"
echo "    --name \"ypestis_co92\" \\"
echo "    --n-jobs 20"

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
echo "- organism_hla_binding_results.xlsx: Final combined results table"
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
