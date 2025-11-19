#!/bin/bash
# Example usage of the Pathogen HLA Binding Analysis Pipeline

# Install dependencies
echo "Installing Python dependencies..."
pip install -r requirements.txt --break-system-packages

# Make the main script executable
chmod +x pathogen_hla_binding.py

# Example 1: Analyze smallpox (Variola virus) - default
echo "Running analysis for Variola virus (smallpox)..."
python pathogen_hla_binding.py \
    --email "your.email@example.com" \
    --organism "Variola virus[Organism]" \
    --output-dir "./smallpox_analysis" \
    --name "smallpox"

# Example 2: Analyze another pathogen (e.g., Monkeypox)
echo "Example command for Monkeypox virus..."
echo "python pathogen_hla_binding.py \\"
echo "    --email \"your.email@example.com\" \\"
echo "    --organism \"Monkeypox virus[Organism]\" \\"
echo "    --output-dir \"./monkeypox_analysis\" \\"
echo "    --name \"monkeypox\""

# Example 3: Analyze specific strain
echo "Example command for specific strain..."
echo "python pathogen_hla_binding.py \\"
echo "    --email \"your.email@example.com\" \\"
echo "    --organism \"Variola virus India-1967[Organism]\" \\"
echo "    --output-dir \"./variola_india_analysis\" \\"
echo "    --name \"variola_india_1967\""

echo "Analysis pipeline setup complete!"
echo ""
echo "IMPORTANT NOTES:"
echo "1. Replace 'your.email@example.com' with your actual email address"
echo "2. Ensure SignalP 6.0+ is installed and in your PATH"
echo "3. Verify netMHCIIpan path is correct: /home/jdn321/software/netMHCIIpan-4.3/netMHCIIpan"
echo "4. The analysis may take considerable time depending on proteome size"
echo ""
echo "Output files:"
echo "- organism_proteome.fasta: Full proteome from NCBI"
echo "- organism_secreted.fasta: Filtered secreted proteins"
echo "- signalp/: SignalP analysis results"
echo "- netmhciipan/: netMHCIIpan prediction files"
echo "- organism_hla_binding_results.xlsx: Final results table"
