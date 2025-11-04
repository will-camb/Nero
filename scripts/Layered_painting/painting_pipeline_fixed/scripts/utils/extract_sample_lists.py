#!/usr/bin/env python3
"""
Extract sample lists from reference population file.
Creates individual files for each population and combined files.
"""

import sys
from pathlib import Path
from collections import defaultdict

def process_ref_popfile(input_file, output_dir):
    """
    Process reference population file and create sample lists.
    
    Input format: sample_id \t population
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read popfile and group by population
    pop_samples = defaultdict(list)
    all_samples = []
    
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) != 2:
                print(f"Warning: Skipping malformed line: {line}")
                continue
            
            sample_id, population = parts
            pop_samples[population].append(sample_id)
            all_samples.append(sample_id)
    
    print(f"Processed {len(all_samples)} total samples across {len(pop_samples)} populations")
    
    # Write individual population files
    for pop, samples in pop_samples.items():
        # Use the full population name (everything before the last dot-separated component if it looks like a classifier)
        # e.g., "IA.Denmark.1575-900BP.XpcCWC.0_1_3" -> "IA.Denmark.1575-900BP"
        # Take everything up to the part that has underscores (classifier codes)
        parts = pop.split('.')
        # Find where the classifier codes start (parts with underscores or just numbers)
        pop_name_parts = []
        for part in parts:
            # if '_' in part or part.isdigit():
            #    break
            pop_name_parts.append(part)

        safe_pop = '.'.join(pop_name_parts) if pop_name_parts else pop

        output_file = output_dir / f"{safe_pop}.txt"
        with open(output_file, 'w') as f:
            f.write('\n'.join(samples))

        print(f"  {output_file}: {len(samples)} samples")
    
    # Write combined reference samples file (all samples, one per line)
    ref_samples_file = output_dir / "ref_samples.txt"
    with open(ref_samples_file, 'w') as f:
        f.write('\n'.join(all_samples))
    print(f"\n  {ref_samples_file}: {len(all_samples)} samples (all references)")
    
    # Also keep the original format popfile for SparsePainter
    # (already exists as ref_popfile.txt)
    
    return pop_samples, all_samples

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python extract_sample_lists.py <ref_popfile> <output_dir>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    
    process_ref_popfile(input_file, output_dir)
