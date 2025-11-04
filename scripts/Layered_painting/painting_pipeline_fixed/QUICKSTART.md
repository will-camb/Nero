# Quick Start Guide - Fixed Pipeline

## Prerequisites

- Access to HPC cluster with SLURM scheduler
- Required modules: `python/3.9.5`, `bcftools/1.20`, `gsl/2.5`, `perl/5.38.0`
- Write access to output directory
- Ancient DNA VCF files
- UK Biobank VCF files

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

data:
  chromosomes: [6]  # Test with chr6 first
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

Then wait ~30-60 minutes for Phase 1 to complete.

## What Happens Next?

### Phase 1: Data Preparation (~30-60 min)

The pipeline will:
1. ✅ Preprocess UKB VCF (add chr info, sample names)
2. ✅ Filter ancient VCF with **FIXED** bcftools command
3. ✅ Find common SNP sites
4. ✅ Extract sample lists

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

### Phase 2: Layer Painting (~24-48 hours)

```bash
./run_pipeline.sh layer 1
```

The pipeline will:
1. Convert VCFs to phase format
2. Generate genetic maps
3. Split targets into batches
4. Estimate lambda
5. Paint ref vs ref
6. Submit painting jobs (parallel)
7. Wait for completion

**Monitor with:**
```bash
watch -n 30 'squeue -u $USER'
```

### Final Step: Merge Results (~10-30 min)

```bash
bash scripts/phase2_layer_painting/08_merge_results.sh 1
```

## Key Changes from Previous Version

### ✅ Fixed: bcftools Filter

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
data/processed/ukb/chr6.ukb.vcf.gz                    ✓
data/processed/ancient/chr6.ancient.filtered.vcf.gz   ✓
data/processed/common_sites/chr6_common_sites.txt     ✓
```

After Phase 2:
```
results/layer1/raw_probabilities/chr6_prob.txt.gz     ✓
results/layer1/chunklengths/chr6_chunklengths.txt.gz  ✓
results/layer1/chunkcounts/chr6_chunkcounts.txt.gz    ✓
```

## Common Issues & Solutions

| Issue | Quick Fix |
|-------|-----------|
| Module not found | `module avail bcftools` to check available versions |
| Path not found | Double-check paths in `master_config.yaml` |
| No common sites | Verify chromosome names match in both VCFs |
| Job fails | Check `logs/paint_layer1_chr6_batch*.err` |
| Out of memory | Reduce `target_batch_size` from 1000 to 500 |

## Tips for Success

1. **Start small**: Test with chr6 and 10K samples first
2. **Check logs**: Look in `logs/` directory if anything fails
3. **Use screen/tmux**: For long-running commands
4. **Monitor jobs**: Use `squeue -u $USER` frequently
5. **Be patient**: Phase 2 takes 24-48 hours for 10K samples

## After Test Success

Once chr6 test completes successfully:

1. Scale to all chromosomes:
   ```yaml
   chromosomes: [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
   ```

2. Use full UKB dataset:
   ```bash
   # Remove or modify ukb_subset.txt constraint
   ```

3. Add more layers (2-5) for temporal depth

## Need Help?

1. Check `README.md` for detailed documentation
2. Review troubleshooting section
3. Examine log files in `logs/` directory
4. Validate config: `./run_pipeline.sh validate`

## Time Budget

- **Setup**: 5 minutes
- **Phase 1 (chr6)**: 30-60 minutes
- **Phase 2 (chr6, 10K samples)**: 24-48 hours
- **Merge**: 10-30 minutes
- **Total**: ~1-2 days for first test

## Success Criteria

✅ Phase 1 completes without errors  
✅ Common sites file has >10,000 SNPs  
✅ All Phase 2 jobs complete successfully  
✅ Results files created in `results/layer1/`  
✅ Can load and inspect probability files  

You're ready to go! 🚀
