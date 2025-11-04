#!/usr/bin/env python3
"""
Process ancient sample metadata into standardized format.

This script reads the provided ancient sample metadata TSV and creates
a cleaned, standardized version for use in the painting pipeline.
"""

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path


def parse_age_average(row):
    """
    Extract mean age in BP from the ageAverage column.
    If ageAverage is not available, calculate from ageHigh and ageLow.
    """
    if pd.notna(row['ageAverage']):
        return row['ageAverage']
    elif pd.notna(row['ageHigh']) and pd.notna(row['ageLow']):
        return (row['ageHigh'] + row['ageLow']) / 2
    else:
        return np.nan


def determine_quality_tier(row):
    """
    Assign quality tier based on coverage and flags.
    Tiers: high, medium, low, exclude
    """
    flag = str(row['flag'])
    coverage = row['coverage']
    
    # Exclude samples with critical quality issues
    exclude_flags = [
        'QUESTIONABLE_CRITICAL', 'FAILED', 'CONTAM', 
        'lowcov', 'contam', 'QUESTIONABLE'
    ]
    
    if any(exc in flag for exc in exclude_flags):
        return 'exclude'
    
    # High quality: coverage >= 1.0, no quality flags
    if pd.notna(coverage) and coverage >= 1.0 and (flag == '0' or flag == 'nan'):
        return 'high'
    
    # Medium quality: coverage >= 0.5
    if pd.notna(coverage) and coverage >= 0.5:
        return 'medium'
    
    # Low quality: coverage < 0.5 but not excluded
    return 'low'


def process_ancient_metadata(input_file, output_file):
    """
    Process ancient sample metadata into standardized format.
    
    Args:
        input_file: Path to input TSV
        output_file: Path to output TSV
    """
    print(f"Reading ancient sample metadata from {input_file}...")
    df = pd.read_csv(input_file, sep='\t', low_memory=False)
    
    print(f"Loaded {len(df)} samples")
    
    # Extract and calculate key fields
    processed = pd.DataFrame({
        'sample_id': df['sampleId'],
        'pop_id': df['popId'],
        'group_label': df['groupLabel'],
        'site': df['site'],
        'country': df['country'],
        'region': df['region'],
        'age_average_bp': df.apply(parse_age_average, axis=1),
        'age_low_bp': df['ageLow'],
        'age_high_bp': df['ageHigh'],
        'dating_source': df['datingSource'],
        'coverage': df['coverage'],
        'sex': df['sex'],
        'flag': df['flag'],
        'quality_tier': df.apply(determine_quality_tier, axis=1),
        'latitude': df['latitude'],
        'longitude': df['longitude'],
        'data_source': df['dataSource'],
        'callset': df['callset']
    })
    
    # Replace empty strings with NaN
    processed = processed.replace('', np.nan)
    
    # Print summary statistics
    print("\n=== Sample Summary ===")
    print(f"Total samples: {len(processed)}")
    print(f"\nQuality tiers:")
    print(processed['quality_tier'].value_counts())
    
    print(f"\nTop 10 population groups:")
    print(processed['group_label'].value_counts().head(10))
    
    print(f"\nCoverage distribution:")
    print(processed['coverage'].describe())
    
    print(f"\nAge distribution (BP):")
    print(processed['age_average_bp'].describe())
    
    # Save processed metadata
    processed.to_csv(output_file, sep='\t', index=False)
    print(f"\nProcessed metadata saved to {output_file}")
    
    # Also create a high-quality subset
    high_qual = processed[processed['quality_tier'].isin(['high', 'medium'])]
    high_qual_file = output_file.replace('.tsv', '_highquality.tsv')
    high_qual.to_csv(high_qual_file, sep='\t', index=False)
    print(f"High-quality subset ({len(high_qual)} samples) saved to {high_qual_file}")
    
    return processed


def create_population_sample_lists(metadata_df, output_dir):
    """
    Create sample ID lists for common population groups.
    This is a helper to quickly create initial sample lists.
    
    Args:
        metadata_df: Processed metadata DataFrame
        output_dir: Directory to save sample lists
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\nCreating sample lists in {output_dir}...")
    
    # Filter to high/medium quality only
    df = metadata_df[metadata_df['quality_tier'].isin(['high', 'medium'])]
    
    # Example: Create lists for common populations
    populations = df['group_label'].value_counts().head(20).index
    
    for pop in populations:
        if pd.isna(pop):
            continue
        samples = df[df['group_label'] == pop]['sample_id'].tolist()
        if len(samples) > 0:
            # Sanitize filename
            safe_name = pop.replace(' ', '_').replace('/', '_').replace('.', '_')
            filename = output_dir / f"{safe_name}.txt"
            with open(filename, 'w') as f:
                f.write('\n'.join(samples))
            print(f"  Created {filename} with {len(samples)} samples")


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python process_ancient_metadata.py <input_tsv> <output_tsv>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Process metadata
    df = process_ancient_metadata(input_file, output_file)
    
    # Optionally create sample lists for common populations
    if len(sys.argv) > 3:
        sample_lists_dir = sys.argv[3]
        create_population_sample_lists(df, sample_lists_dir)
