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
│   │   ├── layer1_config.yaml        # Layer 1 (Migration Period)
│   │   └── layer5_config.yaml        # Layer 5 (Archaic)
│   └── samples/
│       └── layer1/
│           └── ref_popfile.txt       # Reference populations (311 samples)
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
│   │   └── 08_merge_results.sh       # Merge batch results
│   └── utils/
│       ├── extract_sample_lists.py   # Extract sample IDs from popfile
│       └── process_ancient_metadata.py # Process ancient sample metadata
└── logs/                              # Log files (created at runtime)
```

## What Was Fixed

### Phase 1 Data Preparation (`run_all.sh`)

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
  ukb_chr_rename_map: "/datasets/ukb-AUDIT/haplotype_vcf/ukb_chr{chrom}_rename.txt"
  ukb_samples_file: "/datasets/ukb-AUDIT/samples_wb"
  output_base: "/datasets/ukb-AUDIT/layered_painting"
  
data:
  chromosomes: [6]  # Start with chr6 for testing
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

**Duration**: ~30-60 minutes for chr6

**What it does**:
- Preprocesses UKB VCF (adds chr info, sample names)
- Filters ancient VCF (INFO>0.5, MAF>0.01, F_MISSING<0.2)
- Finds common SNP sites between datasets
- Extracts reference sample lists

### 6. Run Phase 2: Layer Painting

```bash
./run_pipeline.sh layer 1
```

**Duration**: ~24-48 hours for 10K samples (chr6)

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
tail -f logs/paint_layer1_chr6_batch*.out
```

### 8. Merge Results

After all painting jobs complete:

```bash
bash scripts/phase2_layer_painting/08_merge_results.sh 1
```

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
tail -f logs/paint_layer1_chr6_batch5.out

# Merge results (after painting completes)
bash scripts/phase2_layer_painting/08_merge_results.sh 1
```

## Expected Outputs

After successful run:

```
/datasets/ukb-AUDIT/layered_painting/
├── data/
│   ├── processed/
│   │   ├── ukb/chr6.ukb.vcf.gz
│   │   ├── ancient/chr6.ancient.filtered.vcf.gz
│   │   └── common_sites/chr6_common_sites.txt
│   └── layer1/
│       ├── ref/chr6.ref.vcf.gz
│       ├── target/chr6.target.vcf.gz
│       ├── ref/chr6.ref.phase
│       ├── target/chr6.target.phase
│       ├── genetic_maps/chr6.map
│       ├── lambda/chr6_lambda.txt
│       ├── ref_vs_ref/chr6_refvsref.chunklengths.out
│       └── target/batches/chr6/batch_*/chr6_batch*.out
└── results/
    └── layer1/
        ├── raw_probabilities/chr6_prob.txt.gz
        ├── chunklengths/chr6_chunklengths.txt.gz
        └── chunkcounts/chr6_chunkcounts.txt.gz
```

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
ls /projects/lundbeck/scratch/vrb229/project_sea/data/250601_impute_sea/merge/vcf/6.*
ls /datasets/ukb-AUDIT/haplotype_vcf/haplotype_chr6.vcf.gz
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
ls logs/paint_layer1_chr6_batch*.err
tail logs/paint_layer1_chr6_batch1.err
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

## Next Steps

After successful test run with chr6:

1. **Expand to all chromosomes**: Edit `master_config.yaml` to include all autosomes
   ```yaml
   chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
   ```

2. **Add more layers**: Create configs for Layers 2-4 (Bronze Age, Neolithic, Mesolithic)

3. **Scale up targets**: Use full UKB sample set instead of subset

4. **Implement Phase 3**: Post-processing scripts for database creation and analysis

## Support

For questions or issues:
1. Check the troubleshooting section above
2. Review log files in the `logs/` directory
3. Verify your configuration with `./run_pipeline.sh validate`

## Version

**Pipeline Version**: 1.1 (Fixed)  
**Date**: 2024  
**Change**: Fixed bcftools filter in Phase 1 data preparation

## License

Custom pipeline for temporal ancestry analysis research.
