#!/usr/bin/env python3
"""
Step 4: Split target samples into batches for parallel painting
"""

import sys
from pathlib import Path

def split_target_phase_file(phase_file, batch_size, output_dir):
    """
    Split target phase file into batches.
    
    SparsePainter phase format:
    Line 1: Number of samples
    Line 2: Number of SNPs
    Line 3: Position marker (P)
    Lines 4+: SNP data
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"  Reading target phase file: {phase_file}")
    
    with open(phase_file, 'r') as f:
        # Read header
        n_samples = int(f.readline().strip())
        n_snps_line = f.readline().strip()
        p_line = f.readline().strip()
        
        # Calculate number of batches
        n_batches = (n_samples + batch_size - 1) // batch_size
        print(f"  Total samples: {n_samples}")
        print(f"  Batch size: {batch_size}")
        print(f"  Number of batches: {n_batches}")
        
        # Read all SNP lines
        snp_lines = f.readlines()
    
    # Split into batches
    # Each haplotype line in phase file represents one haplotype
    # So for N samples, there are 2N lines (2 haplotypes per sample)
    
    # For simplicity, we'll split based on sample IDs
    # Need to re-parse to get sample structure
    
    # Actually, phase files have haplotype data interleaved
    # Let's read the sample IDs from the VCF instead
    pass

def split_targets_simple(chr_num, layer, batch_size, base_dir):
    """
    Simplified approach: create batch files pointing to sections of targets.
    We'll use PBWT's selectSamples feature to create actual batch files.
    """
    layer_dir = Path(base_dir) / "data" / f"layer{layer}"
    target_phase = layer_dir / "target" / f"chr{chr_num}.target.phase"
    target_vcf = layer_dir / "target" / f"chr{chr_num}.target.vcf.gz"
    batch_dir = layer_dir / "target" / "batches" / f"chr{chr_num}"
    batch_dir.mkdir(parents=True, exist_ok=True)
    
    # Get sample list from VCF
    import subprocess
    result = subprocess.run(
        ['bcftools', 'query', '-l', str(target_vcf)],
        capture_output=True, text=True
    )
    samples = result.stdout.strip().split('\n')
    
    n_samples = len(samples)
    n_batches = (n_samples + batch_size - 1) // batch_size
    
    print(f"  Chromosome {chr_num}: {n_samples} samples -> {n_batches} batches")
    
    # Create sample list files for each batch
    for batch_idx in range(n_batches):
        start_idx = batch_idx * batch_size
        end_idx = min((batch_idx + 1) * batch_size, n_samples)
        batch_samples = samples[start_idx:end_idx]
        
        # Write sample list
        batch_samples_file = batch_dir / f"batch{batch_idx+1}_samples.txt"
        with open(batch_samples_file, 'w') as f:
            f.write('\n'.join(batch_samples))
        
        print(f"    Batch {batch_idx+1}: {len(batch_samples)} samples")
    
    # Write batch count file
    with open(batch_dir / "n_batches.txt", 'w') as f:
        f.write(str(n_batches))
    
    return n_batches

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: python split_target_batches.py <layer> <chr> <batch_size>")
        sys.exit(1)
    
    layer = sys.argv[1]
    chr_num = sys.argv[2]
    batch_size = int(sys.argv[3])
    
    # Get base directory (two levels up from script)
    script_dir = Path(__file__).parent
    base_dir = script_dir.parent.parent
    
    # Parse config to get output_base
    config_file = base_dir / "config" / "master_config.yaml"
    with open(config_file) as f:
        for line in f:
            if 'output_base:' in line:
                output_base = line.split(':', 1)[1].strip().strip('"')
                break
    
    print(f"Step 4: Splitting targets into batches for chr{chr_num}")
    n_batches = split_targets_simple(chr_num, layer, batch_size, output_base)
    print(f"  âœ“ Created {n_batches} batches for chromosome {chr_num}")
