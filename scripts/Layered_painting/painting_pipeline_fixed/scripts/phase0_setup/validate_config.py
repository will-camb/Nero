#!/usr/bin/env python3
"""
Validate pipeline configuration and check that all required files exist.
"""

import yaml
import sys
from pathlib import Path
import re


def load_config(config_file):
    """Load YAML configuration file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)


def expand_path(path_template, config):
    """Expand path template with environment variables."""
    # Replace ${variable} with values from config
    if path_template is None:
        return None
    
    path_str = str(path_template)
    # Replace ${output_base} etc.
    for key, value in config.get('paths', {}).items():
        path_str = path_str.replace(f"${{{key}}}", str(value))
    
    return Path(path_str)


def check_file_pattern(pattern, name, required=True):
    """Check if a file pattern is valid (doesn't check if files exist)."""
    if pattern is None:
        if required:
            print(f"  âŒ ERROR: {name} is not specified")
            return False
        else:
            print(f"  âš ï¸  WARNING: {name} is not specified (optional)")
            return True
    
    # Check if pattern contains chromosome placeholder
    if '{' in pattern or 'chr' in pattern.lower():
        print(f"  âœ“ {name}: {pattern} (pattern)")
    else:
        print(f"  âœ“ {name}: {pattern}")
    
    return True


def validate_master_config(config_file):
    """Validate master configuration file."""
    print("=" * 60)
    print("Validating Master Configuration")
    print("=" * 60)
    
    try:
        config = load_config(config_file)
    except Exception as e:
        print(f"âŒ ERROR: Could not load {config_file}: {e}")
        return False
    
    valid = True
    
    # Check required sections
    required_sections = ['paths', 'global_filters', 'sparsepainter', 'hpc']
    for section in required_sections:
        if section not in config:
            print(f"âŒ ERROR: Missing required section: {section}")
            valid = False
    
    # Check paths
    print("\nðŸ“ Checking paths:")
    paths = config.get('paths', {})
    
    # Check data input paths
    check_file_pattern(paths.get('ancient_vcf'), "Ancient VCF", required=True)
    check_file_pattern(paths.get('ukb_bgen'), "UKB BGEN/VCF", required=False)
    check_file_pattern(paths.get('kg_vcf'), "1000 Genomes VCF", required=False)
    check_file_pattern(paths.get('genetic_maps'), "Genetic maps", required=True)
    
    # Check software paths
    sparsepainter = paths.get('sparsepainter')
    pbwt = paths.get('pbwt')
    
    if not sparsepainter:
        print(f"  âŒ ERROR: SparsePainter path not specified")
        valid = False
    else:
        print(f"  âœ“ SparsePainter: {sparsepainter}")
    
    if not pbwt:
        print(f"  âŒ ERROR: PBWT path not specified")
        valid = False
    else:
        print(f"  âœ“ PBWT: {pbwt}")
    
    # Check ancient metadata
    metadata_path = paths.get('ancient_metadata')
    if metadata_path:
        full_path = Path(config_file).parent.parent / metadata_path
        if full_path.exists():
            print(f"  âœ“ Ancient metadata: {metadata_path} (exists)")
        else:
            print(f"  âš ï¸  WARNING: Ancient metadata file not found: {full_path}")
    
    # Check filters
    print("\nðŸ” Checking global filters:")
    filters = config.get('global_filters', {})
    info_thresh = filters.get('info_threshold')
    maf = filters.get('maf_global')
    
    if info_thresh is not None:
        print(f"  âœ“ INFO threshold: {info_thresh}")
        if info_thresh < 0 or info_thresh > 1:
            print(f"    âš ï¸  WARNING: INFO threshold should be between 0 and 1")
    
    if maf is not None:
        print(f"  âœ“ Global MAF: {maf}")
        if maf < 0 or maf > 0.5:
            print(f"    âš ï¸  WARNING: MAF should be between 0 and 0.5")
    
    # Check HPC settings
    print("\nðŸ’» Checking HPC settings:")
    hpc = config.get('hpc', {})
    scheduler = hpc.get('scheduler')
    if scheduler not in ['slurm', 'pbs']:
        print(f"  âš ï¸  WARNING: Scheduler '{scheduler}' not recognized (should be 'slurm' or 'pbs')")
    else:
        print(f"  âœ“ Scheduler: {scheduler}")
    
    print(f"  âœ“ Partition: {hpc.get('partition')}")
    print(f"  âœ“ Time limit: {hpc.get('time_limit')}")
    print(f"  âœ“ Memory: {hpc.get('mem_per_job')}")
    
    return valid


def validate_layer_config(layer_config_file, master_config_file):
    """Validate layer configuration file."""
    print("\n" + "=" * 60)
    print(f"Validating Layer Configuration: {layer_config_file}")
    print("=" * 60)
    
    try:
        layer_config = load_config(layer_config_file)
        master_config = load_config(master_config_file)
    except Exception as e:
        print(f"âŒ ERROR: Could not load config: {e}")
        return False
    
    valid = True
    
    # Check required fields
    layer_id = layer_config.get('layer_id')
    layer_name = layer_config.get('layer_name')
    
    if layer_id is None:
        print("âŒ ERROR: layer_id not specified")
        valid = False
    else:
        print(f"âœ“ Layer ID: {layer_id}")
    
    if layer_name is None:
        print("âŒ ERROR: layer_name not specified")
        valid = False
    else:
        print(f"âœ“ Layer name: {layer_name}")
    
    # Check reference populations
    print("\nðŸ‘¥ Checking reference populations:")
    ref_pops = layer_config.get('reference_populations', [])
    
    if len(ref_pops) == 0:
        print("âŒ ERROR: No reference populations specified")
        valid = False
    else:
        print(f"âœ“ Found {len(ref_pops)} reference populations:")
        
        config_dir = Path(layer_config_file).parent.parent
        
        for i, pop in enumerate(ref_pops, 1):
            pop_name = pop.get('name')
            sample_file = pop.get('sample_file')
            time_range = pop.get('time_range_bp')
            
            print(f"\n  Population {i}: {pop_name}")
            
            if not sample_file:
                print(f"    âŒ ERROR: sample_file not specified")
                valid = False
            else:
                full_path = config_dir / sample_file
                if full_path.exists():
                    # Count samples in file
                    with open(full_path) as f:
                        n_samples = sum(1 for line in f if line.strip())
                    print(f"    âœ“ Sample file: {sample_file} ({n_samples} samples)")
                else:
                    print(f"    âš ï¸  WARNING: Sample file not found: {full_path}")
                    print(f"       You will need to create this file before running the pipeline")
            
            if time_range:
                print(f"    âœ“ Time range: {time_range[0]}-{time_range[1]} BP")
            else:
                print(f"    âš ï¸  WARNING: time_range_bp not specified")
    
    # Check filters
    print("\nðŸ” Checking layer-specific filters:")
    filters = layer_config.get('filters', {})
    if filters:
        print(f"  âœ“ Reference max missing: {filters.get('ref_max_missing')}")
        print(f"  âœ“ Target max missing: {filters.get('target_max_missing')}")
        print(f"  âœ“ MAF threshold: {filters.get('maf_threshold')}")
    
    return valid


def validate_all_layers(config_dir):
    """Validate all layer configuration files."""
    config_path = Path(config_dir)
    master_config = config_path / 'master_config.yaml'
    
    if not master_config.exists():
        print(f"âŒ ERROR: Master config not found: {master_config}")
        return False
    
    # Validate master config first
    if not validate_master_config(master_config):
        print("\nâŒ Master configuration validation FAILED")
        return False
    
    print("\nâœ“ Master configuration validation PASSED")
    
    # Find and validate all layer configs
    layers_dir = config_path / 'layers'
    if not layers_dir.exists():
        print(f"\nâš ï¸  WARNING: Layers directory not found: {layers_dir}")
        return True
    
    layer_configs = sorted(layers_dir.glob('layer*_config.yaml'))
    
    if len(layer_configs) == 0:
        print(f"\nâš ï¸  WARNING: No layer configurations found in {layers_dir}")
        return True
    
    all_valid = True
    for layer_config in layer_configs:
        if not validate_layer_config(layer_config, master_config):
            all_valid = False
            print(f"\nâŒ Layer configuration validation FAILED: {layer_config.name}")
        else:
            print(f"\nâœ“ Layer configuration validation PASSED: {layer_config.name}")
    
    return all_valid


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python validate_config.py <config_directory>")
        print("Example: python validate_config.py /path/to/painting_pipeline/config")
        sys.exit(1)
    
    config_dir = sys.argv[1]
    
    print("ðŸ”§ Pipeline Configuration Validator")
    print("=" * 60)
    
    if validate_all_layers(config_dir):
        print("\n" + "=" * 60)
        print("âœ… ALL VALIDATIONS PASSED")
        print("=" * 60)
        sys.exit(0)
    else:
        print("\n" + "=" * 60)
        print("âŒ VALIDATION FAILED - Please fix the errors above")
        print("=" * 60)
        sys.exit(1)
