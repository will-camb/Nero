# Quick Start Guide - Pipeline v1.2

## Prerequisites

- Access to HPC cluster with SLURM scheduler
- Required modules: `python/3.9.5`, `bcftools/1.20`, `gsl/2.5`, `perl/5.38.0`
- Write access to output directory
- Ancient DNA VCF files
- UK Biobank VCF files

## What's New in v1.2

- ✅ Fixed critical bugs (syntax errors, Python execution, sample naming)
- ✅ All 5 temporal layers now configured (Layers 2-4 need sample files)
- ✅ Testing on chr22 (faster than chr6)
- ✅ Lambda default updated to 50
- ✅ Phase 3 properly marked as not yet implemented

## 5-Minute Setup

### Step 1: Extract (30 seconds)

```bash
cd /datasets/ukb-AUDIT
tar -xzf painting_pipeline_fixed.tar.gz
cd painting_pipeline_fixed
chmod +x run_pipeline.sh scripts/*/*.sh
```

### Step 2: Configure (2 minutes)

Edit `config/master_config.yaml`:

```yaml
paths:
  output_base: "/datasets/ukb-AUDIT/layered_painting"  # Your output directory

chromosomes: [22]  # Test with chr22 first (faster than chr6)
```

All other paths should already be correct for your setup.

### Step 3: Create Sample List (30 seconds)

```bash
head -n 10000 /datasets/ukb-AUDIT/samples_wb > config/samples/ukb_subset.txt
```

### Step 4: Validate (30 seconds)

```bash
python3 scripts/phase0_setup/validate_config.py config/
```

Should see: `✅ ALL VALIDATIONS PASSED`

### Step 5: Run (1 minute to start)

```bash
# Start Phase 1
./run_pipeline.sh prep
```

Then wait ~20-40 minutes for Phase 1 to complete (chr22 is smaller than chr6).

## What Happens Next?

### Phase 1: Data Preparation (~20-40 min for chr22)

The pipeline will:
1. ✅ Preprocess UKB VCF (add chr info, sample names)
2. ✅ Filter ancient VCF with **FIXED** bcftools command
3. ✅ Find common SNP sites
4. ✅ Extract sample lists with **FIXED** naming (preserves full population names)

**You'll know it worked when you see:**
```
==================================================
Phase 1 Complete!
==================================================

Next steps:
  1. Review the common sites counts
  2. Create a UKB subset file: config/samples/ukb_subset.txt
  3. Run Phase 2: ./run_pipeline.sh layer 1
```

### Phase 2: Layer Painting (~16-32 hours for chr22)

```bash
./run_pipeline.sh layer 1
```

The pipeline will:
1. Convert VCFs to phase format
2. Generate genetic maps
3. Split targets into batches
4. Estimate lambda (default: 50)
5. Paint ref vs ref
6. Submit painting jobs (parallel)
7. Wait for completion

**Note**: Chr22 is smaller than chr6, so painting is faster.

**Monitor with:**
```bash
watch -n 30 'squeue -u $USER'
```

### Final Step: Merge Results (~10-30 min)

```bash
bash scripts/phase2_layer_painting/08_merge_results.sh 1
```

## Key Changes from Previous Versions

### ✅ Version 1.2 Fixes

**Critical bug fixes:**
1. Fixed syntax error in `run_pipeline.sh` (Python `else:` → bash `else`)
2. Fixed broken Python script execution in `08_merge_results.sh` (now uses temp files)
3. Fixed sample naming in `extract_sample_lists.py` (preserves full population names with time ranges)
4. Fixed path inconsistencies (genetic maps now use `{chrom}` consistently)
5. Removed redundant module loads

**Configuration updates:**
- Test chromosome: chr6 → chr22 (faster testing)
- Lambda default: 25.7926 → 50
- All 5 layers now have configurations

### ✅ Version 1.1 Fix: bcftools Filter

**OLD (broken):**
```bash
bcftools view --max-missing-count $(calculation) ...
# Error: unknown option '--max-missing-count'
```

**NEW (working):**
```bash
bcftools view -i "INFO/INFO>0.5 && MAF>0.01 && F_MISSING<0.2" ...
# Uses correct F_MISSING field
```

## Expected Results

After Phase 1:
```
data/processed/ukb/chr22.ukb.vcf.gz                    ✓
data/processed/ancient/chr22.ancient.filtered.vcf.gz   ✓
data/processed/common_sites/chr22_common_sites.txt     ✓
```

After Phase 2:
```
results/layer1/raw_probabilities/chr22_prob.txt.gz     ✓
results/layer1/chunklengths/chr22_chunklength.txt.gz   ✓
```

## Common Issues & Solutions

| Issue | Quick Fix |
|-------|-----------|
| Module not found | `module avail bcftools` to check available versions |
| Path not found | Double-check paths in `master_config.yaml` |
| No common sites | Verify chromosome names match in both VCFs |
| Job fails | Check `logs/paint_layer1_chr22_batch*.err` |
| Out of memory | Reduce `target_batch_size` from 1000 to 500 |

## Tips for Success

1. **Start small**: Test with chr22 and 10K samples first (faster than chr6)
2. **Check logs**: Look in `logs/` directory if anything fails
3. **Use screen/tmux**: For long-running commands
4. **Monitor jobs**: Use `squeue -u $USER` frequently
5. **Be patient**: Phase 2 takes 16-32 hours for 10K samples on chr22

## After Test Success

Once chr22 test completes successfully:

1. Add samples for other layers (2-5):
   - See `config/samples/layer{2,3,4,5}/README.md` for instructions
   - Create `ref_popfile.txt` for each layer
   - Run `extract_sample_lists.py` to generate population files

2. Scale to all chromosomes:
   ```yaml
   chromosomes: [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
   ```

3. Use full UKB dataset:
   ```bash
   # Remove or modify ukb_subset.txt constraint
   ```

4. Run all layers for complete temporal depth analysis

## Need Help?

1. Check `README.md` for detailed documentation
2. Review troubleshooting section
3. Examine log files in `logs/` directory
4. Validate config: `./run_pipeline.sh validate`

## Time Budget

- **Setup**: 5 minutes
- **Phase 1 (chr22)**: 20-40 minutes
- **Phase 2 (chr22, 10K samples)**: 16-32 hours
- **Merge**: 10-20 minutes
- **Total**: ~1 day for first test (chr22 is faster than chr6)

## Success Criteria

✅ Phase 1 completes without errors  
✅ Common sites file has >10,000 SNPs  
✅ All Phase 2 jobs complete successfully  
✅ Results files created in `results/layer1/`  
✅ Can load and inspect probability files  

You're ready to go! 🚀
