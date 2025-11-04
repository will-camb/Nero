# Temporal Ancestry Painting Pipeline

A complete pipeline for multi-layer chromosome painting of modern samples against ancient DNA reference panels.

## What's Included

This package contains all the scripts needed to run the painting pipeline from start to finish.

### Directory Structure

```
painting_pipeline_fixed/
├── run_pipeline.sh                    # Master orchestrator script
├── config/
│   ├── master_config.yaml            # Global configuration
│   ├── layers/
│   │   ├── layer1_config.yaml        # Layer 1 (Migration Period, 900-2800 BP)
│   │   ├── layer2_config.yaml        # Layer 2 (Bronze Age, 2800-4500 BP)
│   │   ├── layer3_config.yaml        # Layer 3 (Neolithic, 4500-8000 BP)
│   │   ├── layer4_config.yaml        # Layer 4 (Mesolithic, 8000-12000 BP)
│   │   └── layer5_config.yaml        # Layer 5 (Archaic, >50000 BP)
│   └── samples/
│       ├── layer1/
│       │   └── ref_popfile.txt       # Reference populations (311 samples)
│       ├── layer2/                   # Bronze Age samples (placeholder)
│       ├── layer3/                   # Neolithic samples (placeholder)
│       ├── layer4/                   # Mesolithic samples (placeholder)
│       └── layer5/                   # Archaic samples (placeholder)
├── scripts/
│   ├── phase0_setup/
│   │   └── validate_config.py        # Configuration validation
│   ├── phase1_data_prep/
│   │   └── run_all.sh                # Data preparation (FIXED VERSION)
│   ├── phase2_layer_painting/
│   │   ├── run_layer.sh              # Layer orchestrator
│   │   ├── 01_prepare_ref_target.sh  # Extract ref/target VCFs
│   │   ├── 02_vcf_to_phase.sh        # Convert to phase format
│   │   ├── 03_make_genetic_maps.sh   # Generate genetic maps
│   │   ├── 04_split_target_batches.py # Split into batches
│   │   ├── 05_estimate_lambda.sh     # Estimate lambda parameter
│   │   ├── 06_paint_ref_vs_ref.sh    # Paint reference vs reference
│   │   ├── 07_submit_paint_targets.sh # Paint target batches
│   │   └── 08_merge_results.sh       # Merge batch results (FIXED VERSION)
│   └── utils/
│       ├── extract_sample_lists.py   # Extract sample IDs from popfile (FIXED VERSION)
│       └── process_ancient_metadata.py # Process ancient sample metadata
└── logs/                              # Log files (created at runtime)
```

## What Was Fixed

### Version 1.2 (Current)

#### Critical Bug Fixes

1. **Syntax Error in run_pipeline.sh**
   - **Problem**: Python syntax (`else:`) used instead of bash syntax
   - **Impact**: Script would crash when running postprocess command
   - **Fixed**: Changed to correct bash syntax (`else`)

2. **Broken Python Script Execution in 08_merge_results.sh**
   - **Problem**: Heredoc Python scripts weren't executing properly
   - **Impact**: Result merging would fail silently
   - **Fixed**: Now creates temporary Python scripts, executes them, and cleans up

3. **Sample List Naming Mismatch in extract_sample_lists.py**
   - **Problem**: Created files like `IA.Denmark.txt` but configs expected `IA.Denmark.1575-900BP.txt`
   - **Impact**: Sample files wouldn't match layer configuration expectations
   - **Fixed**: Now preserves full population names including time ranges

4. **Path Inconsistencies**
   - **Problem**: Genetic map path used `{i}` instead of `{chrom}`
   - **Fixed**: Standardized to use `{chrom}` throughout

5. **Redundant Module Loads**
   - **Problem**: Modules loaded twice in run_layer.sh
   - **Fixed**: Removed duplicate module load commands

#### Configuration Updates

- **Test Chromosome**: Changed from chr6 to chr22
- **Lambda Parameter**: Updated default from 25.7926 to 50 (configurable, will parse from SparsePainter in future)
- **Number of Layers**: Updated from 1 to 5 to reflect all available configurations

#### New Features

- **Complete Layer Configurations**: Added placeholder configs for Layers 2-4
  - Layer 2: Bronze Age (2800-4500 BP)
  - Layer 3: Neolithic (4500-8000 BP)
  - Layer 4: Mesolithic (8000-12000 BP)
- **Sample Directories**: Created directories with documentation for all layers
- **Phase 3 Handling**: Removed from main workflow (not yet implemented), kept as optional command

### Version 1.1

#### Phase 1 Data Preparation (`run_all.sh`)

**Problem**: Used non-existent `--max-missing-count` option in bcftools

**Solution**: Replaced with correct filter expression:
```bash
# OLD (BROKEN):
bcftools view --max-missing-count $(calculation) ...

# NEW (FIXED):
bcftools view -i "INFO/INFO>$INFO_THRESH && MAF>$MAF_THRESH && F_MISSING<$MAX_MISSING" ...
```

Now uses `F_MISSING` (fraction of missing genotypes) which is the standard bcftools field.

## Quick Start

### 1. Extract and Setup

```bash
# Extract to your working directory
cd /datasets/ukb-AUDIT
tar -xzf painting_pipeline_fixed.tar.gz
cd painting_pipeline_fixed

# Make scripts executable
chmod +x run_pipeline.sh
chmod +x scripts/phase1_data_prep/run_all.sh
chmod +x scripts/phase2_layer_painting/*.sh
```

### 2. Configure Paths

Edit `config/master_config.yaml` to set your paths:

```yaml
paths:
  ancient_vcf: "/projects/lundbeck/scratch/vrb229/project_sea/data/250601_impute_sea/merge/vcf/{chrom}.250601_impute_sea.glimpse.vcf.gz"
  ukb_haplotype: "/datasets/ukb-AUDIT/haplotype_vcf/haplotype_chr{chrom}.vcf.gz"
  ukb_chr_rename_map: "/datasets/ukb-AUDIT/merged_UKBB_ref/pbwt/merged_vcfs/chr_rename_maps/{chrom}.map"
  ukb_samples_file: "/datasets/ukb-AUDIT/merged_UKBB_ref/pbwt/merged_vcfs/UKBB_samples"
  output_base: "/datasets/ukb-AUDIT/layered_painting"

chromosomes: [22]  # Start with chr22 for testing
```

### 3. Create UKB Subset

```bash
# Create subset for testing (10K samples)
head -n 10000 /datasets/ukb-AUDIT/samples_wb > config/samples/ukb_subset.txt
```

### 4. Validate Configuration

```bash
python3 scripts/phase0_setup/validate_config.py config/
```

Expected output: `✅ ALL VALIDATIONS PASSED`

### 5. Run Phase 1: Data Preparation

```bash
./run_pipeline.sh prep
```

**Duration**: ~20-40 minutes for chr22 (smaller than chr6)

**What it does**:
- Preprocesses UKB VCF (adds chr info, sample names)
- Filters ancient VCF (INFO>0.5, MAF>0.01, F_MISSING<0.2)
- Finds common SNP sites between datasets
- Extracts reference sample lists

### 6. Run Phase 2: Layer Painting

```bash
./run_pipeline.sh layer 1
```

**Duration**: ~16-32 hours for 10K samples (chr22, smaller than chr6)

**What it does**:
- Converts VCFs to phase format
- Generates genetic maps
- Splits targets into batches (1000 each)
- Estimates lambda parameter
- Paints reference vs reference
- Paints all target batches (parallel SLURM jobs)

### 7. Monitor Jobs

```bash
# Watch job queue
watch -n 30 'squeue -u $USER'

# Check logs
tail -f logs/paint_layer1_chr22_batch*.out
```

### 8. Merge Results

After all painting jobs complete:

```bash
bash scripts/phase2_layer_painting/08_merge_results.sh 1
```

**Note**: Phase 3 (post-processing) is planned for a future release. The `postprocess` command will show you what's coming.

## Configuration Files

### Master Config (`config/master_config.yaml`)

Global settings for all layers:
- File paths (VCFs, tools, output directory)
- Quality filters (INFO, MAF, missingness)
- HPC settings (partition, memory, time limits)
- Chromosomes to process

### Layer Configs (`config/layers/layer{N}_config.yaml`)

Layer-specific settings:
- Reference populations
- Time period
- Layer-specific filters
- Painting parameters

## Key Features

### ✓ Fixed bcftools Filter
- Uses correct `F_MISSING` field
- Single-pass filtering (more efficient)
- Properly handles INFO, MAF, and missingness thresholds

### ✓ Checkpoint/Resume
- All scripts check for existing output files
- Automatically skips completed steps
- Safe to re-run after failures

### ✓ Layer Independence
- Each layer runs independently
- Can process in any order
- Results stored separately

### ✓ Scalability
- Batch processing for large datasets
- SLURM array jobs for parallelization
- Configurable batch sizes

### ✓ Quality Control
- Configuration validation before running
- Progress reporting throughout
- Detailed logging

## Commands Reference

```bash
# Validate configuration
./run_pipeline.sh validate

# Run Phase 1 (data prep)
./run_pipeline.sh prep

# Run specific layer
./run_pipeline.sh layer 1

# Run all layers
./run_pipeline.sh all-layers

# Monitor jobs
squeue -u $USER
watch -n 30 'squeue -u $USER'

# Check specific log
tail -f logs/paint_layer1_chr22_batch5.out

# Merge results (after painting completes)
bash scripts/phase2_layer_painting/08_merge_results.sh 1

# Post-processing (not yet implemented)
./run_pipeline.sh postprocess
```

## Expected Outputs

After successful run for Layer 1 on chr22:

```
/datasets/ukb-AUDIT/layered_painting/
├── data/
│   ├── processed/
│   │   ├── ukb/chr22.ukb.vcf.gz
│   │   ├── ancient/chr22.ancient.filtered.vcf.gz
│   │   └── common_sites/chr22_common_sites.txt
│   └── layer1/
│       ├── ref/chr22.ref.vcf.gz
│       ├── target/chr22.target.vcf.gz
│       ├── ref/chr22.ref.phase
│       ├── target/chr22.target.phase
│       ├── genetic_maps/chr22.map
│       ├── lambda/chr22_lambda.txt
│       ├── ref_vs_ref/chr22_refvsref.chunklengths.out
│       └── target/batches/chr22/batch_*/chr22_batch*.out
└── results/
    └── layer1/
        ├── raw_probabilities/chr22_prob.txt.gz
        └── chunklengths/chr22_chunklength.txt.gz
```

Each additional layer (2-5) will create similar structures in `data/layer{N}/` and `results/layer{N}/`.

## Resource Requirements

### Phase 1 (per chromosome)
- CPU: 1-2 cores
- Memory: ~8GB
- Time: ~30-60 minutes
- Storage: ~15GB

### Phase 2 (per layer, per chromosome)
- Setup (steps 1-6): 4-8 cores, ~16GB RAM, ~1-2 hours
- Painting (step 7): 1 core per batch, ~8GB RAM, ~2-4 hours per batch
- Total with 10K samples: ~20-40 hours (parallelized)
- Storage: ~50-100GB per layer

## Troubleshooting

### Module Load Failures
**Problem**: Required modules not found  
**Solution**: Check which modules are available on your system
```bash
module avail bcftools
module avail python
module avail gsl
```

### Path Errors
**Problem**: File not found errors  
**Solution**: Verify all paths in `master_config.yaml` are correct
```bash
ls /projects/lundbeck/scratch/vrb229/project_sea/data/250601_impute_sea/merge/vcf/22.*
ls /datasets/ukb-AUDIT/haplotype_vcf/haplotype_chr22.vcf.gz
ls /datasets/ukb-AUDIT/merged_UKBB_ref/pbwt/merged_vcfs/chr_rename_maps/22.map
```

### No Common Sites Found
**Problem**: Common sites file is empty
**Solution**: Check that both VCFs have matching chromosome names
```bash
bcftools query -f '%CHROM\n' ancient.vcf.gz | head
bcftools query -f '%CHROM\n' ukb.vcf.gz | head
```

### SLURM Job Failures
**Problem**: Painting jobs fail
**Solution**: Check logs for specific errors
```bash
ls logs/paint_layer1_chr22_batch*.err
tail logs/paint_layer1_chr22_batch1.err
```

### Out of Memory
**Problem**: Jobs killed due to memory  
**Solution**: Reduce batch size in `master_config.yaml`
```yaml
sparsepainter:
  parallelization:
    target_batch_size: 500  # Reduce from 1000
```

## Reference Populations (Layer 1)

Your reference panel includes 311 ancient samples from 6 populations:

| Population | Time Period | Count |
|------------|-------------|-------|
| IA.Denmark.1575-900BP | Danish Vikings | ~63 |
| England.Anglo-Saxon.1450-1150BP | Anglo-Saxons | ~95 |
| IA.France.2400-900BP | Continental Celts | ~58 |
| IA.Norway.1850-1250BP | Norwegian Vikings | ~29 |
| IA.ScotlandIreland.2800-900BP | Insular Celts | ~45 |
| IA.Sweden.1575-900BP | Swedish Vikings | ~20 |

## Available Layers

The pipeline now includes configurations for all 5 temporal layers:

| Layer | Period | Time Range | Status |
|-------|--------|------------|--------|
| 1 | Migration Period | 900-2800 BP | Ready (samples provided) |
| 2 | Bronze Age | 2800-4500 BP | Placeholder (needs samples) |
| 3 | Neolithic | 4500-8000 BP | Placeholder (needs samples) |
| 4 | Mesolithic | 8000-12000 BP | Placeholder (needs samples) |
| 5 | Archaic | >50000 BP | Placeholder (needs samples) |

For layers 2-5, you need to:
1. Create `ref_popfile.txt` in the layer's sample directory
2. Run `extract_sample_lists.py` to generate individual population files
3. See the README.md in each layer's sample directory for details

## Next Steps

After successful test run with chr22:

1. **Add samples for other layers**: Create sample files for Layers 2-5 using your ancient VCF data

2. **Expand to all chromosomes**: Edit `master_config.yaml` to include all autosomes
   ```yaml
   chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
   ```

3. **Scale up targets**: Use full UKB sample set instead of subset

4. **Implement Phase 3**: Post-processing scripts for database creation and analysis (coming in v2.0)

## Support

For questions or issues:
1. Check the troubleshooting section above
2. Review log files in the `logs/` directory
3. Verify your configuration with `./run_pipeline.sh validate`

## Version

**Pipeline Version**: 1.2
**Release Date**: 2024
**Changes**:
- Fixed critical bugs (syntax errors, Python execution, sample naming)
- Added complete layer configurations (Layers 2-4)
- Updated to test on chr22 instead of chr6
- Improved Phase 3 handling
- Updated lambda default to 50

See CHANGELOG.md for full version history.

## License

Custom pipeline for temporal ancestry analysis research.
