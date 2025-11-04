# Changelog

## Version 1.1 - Fixed (Current Release)

### 🔧 Bug Fixes

#### Critical: Phase 1 Data Preparation
- **Fixed bcftools filter syntax error** in `scripts/phase1_data_prep/run_all.sh`
  - **Problem**: Used non-existent `--max-missing-count` option
  - **Impact**: Phase 1 would fail with "unknown option" error
  - **Solution**: Replaced with correct `F_MISSING` field in filter expression
  
**Before (Lines 95-104):**
```bash
bcftools view \
    $ANCIENT_VCF_CHR \
    -i "INFO/INFO>$INFO_THRESH" \
    -Ou | \
bcftools view \
    --min-af $MAF_THRESH:minor \
    -Ou | \
bcftools view \
    --max-missing-count $(echo "$MAX_MISSING * $(bcftools query -l $ANCIENT_VCF_CHR | wc -l)" | bc | cut -d. -f1) \
    -Oz -o $ANCIENT_FILTERED
```

**After (Lines 95-99):**
```bash
bcftools view \
    $ANCIENT_VCF_CHR \
    -i "INFO/INFO>$INFO_THRESH && MAF>$MAF_THRESH && F_MISSING<$MAX_MISSING" \
    -Oz -o $ANCIENT_FILTERED
```

**Benefits:**
- ✅ Uses correct bcftools syntax
- ✅ Single-pass filtering (more efficient)
- ✅ Clearer filter logic
- ✅ No intermediate piping needed

### 📦 Package Improvements

- Consolidated all essential scripts into single downloadable package
- Removed obsolete/duplicate files
- Added comprehensive documentation:
  - `README.md` - Complete pipeline documentation
  - `QUICKSTART.md` - 5-minute setup guide
  - `CHANGELOG.md` - This file

### 📁 File Organization

**Included files:**
```
painting_pipeline_fixed/
├── run_pipeline.sh                    # Master orchestrator
├── README.md                          # Full documentation
├── QUICKSTART.md                      # Quick setup guide
├── CHANGELOG.md                       # Version history
├── config/
│   ├── master_config.yaml            # Global config
│   ├── layers/
│   │   ├── layer1_config.yaml        # Layer 1 config
│   │   └── layer5_config.yaml        # Layer 5 config
│   └── samples/
│       └── layer1/
│           └── ref_popfile.txt       # Reference populations
└── scripts/
    ├── phase0_setup/
    │   └── validate_config.py        # Config validation
    ├── phase1_data_prep/
    │   └── run_all.sh                # **FIXED VERSION**
    ├── phase2_layer_painting/
    │   ├── run_layer.sh              # Layer orchestrator
    │   ├── 01_prepare_ref_target.sh
    │   ├── 02_vcf_to_phase.sh
    │   ├── 03_make_genetic_maps.sh
    │   ├── 04_split_target_batches.py
    │   ├── 05_estimate_lambda.sh
    │   ├── 06_paint_ref_vs_ref.sh
    │   ├── 07_submit_paint_targets.sh
    │   └── 08_merge_results.sh
    └── utils/
        ├── extract_sample_lists.py
        └── process_ancient_metadata.py
```

**Excluded files:**
- Obsolete documentation files
- Test/development scripts
- Duplicate versions
- Notes and planning documents

---

## Version 1.0 - Initial Release

### Features

#### Phase 0: Setup & Validation
- Configuration system (master + per-layer configs)
- Configuration validation script
- Sample list extraction utilities
- Ancient metadata processing

#### Phase 1: Data Preparation
- UKB VCF preprocessing (chromosome info, sample names)
- Ancient VCF quality filtering
- Common site identification
- Sample list management

#### Phase 2: Layer Painting
- Complete 8-step painting workflow
- VCF to phase format conversion
- Genetic map generation
- Target batch splitting
- Lambda parameter estimation
- Reference vs reference painting
- Target painting with parallelization
- Result merging

#### Configuration
- YAML-based configuration system
- Global and per-layer settings
- Flexible population definitions
- Temporal filtering logic
- HPC/SLURM integration

#### Documentation
- Comprehensive setup guides
- Test run instructions
- Interactive mode guide
- Troubleshooting documentation

### Known Issues (Fixed in 1.1)

- ❌ bcftools `--max-missing-count` option doesn't exist
- ❌ Phase 1 would fail with "unknown option" error
- ❌ Required manual editing of script to fix

---

## Migration Guide: 1.0 → 1.1

If you were using version 1.0, here's how to upgrade:

### Option 1: Fresh Install (Recommended)
```bash
# Extract new version
cd /datasets/ukb-AUDIT
tar -xzf painting_pipeline_fixed.tar.gz
cd painting_pipeline_fixed

# Copy your old configs (if customized)
cp ../painting_pipeline/config/master_config.yaml config/
cp ../painting_pipeline/config/samples/ukb_subset.txt config/samples/

# Validate and run
./run_pipeline.sh validate
```

### Option 2: In-Place Update
```bash
# Just replace the broken script
cd /datasets/ukb-AUDIT/painting_pipeline
cp /path/to/new/run_all.sh scripts/phase1_data_prep/run_all.sh
```

### What to Check After Update

1. **Verify the fix:**
   ```bash
   grep "F_MISSING" scripts/phase1_data_prep/run_all.sh
   # Should show: -i "INFO/INFO>$INFO_THRESH && MAF>$MAF_THRESH && F_MISSING<$MAX_MISSING"
   ```

2. **Test Phase 1:**
   ```bash
   ./run_pipeline.sh prep
   # Should complete without bcftools errors
   ```

3. **Check outputs:**
   ```bash
   ls -lh data/processed/ancient/chr*.ancient.filtered.vcf.gz
   # Should have reasonable file sizes
   ```

---

## Future Plans

### Version 1.2 (Planned)
- [ ] Implement Phase 3: Post-processing scripts
- [ ] Add database creation for ancestry paths
- [ ] VCF annotation with ancestry information
- [ ] Summary statistics generation

### Version 1.3 (Planned)
- [ ] Implement Phase 4: Analysis tools
- [ ] Ancestry path query scripts
- [ ] CLUES export functionality
- [ ] Visualization tools

### Version 2.0 (Future)
- [ ] Python-based configuration management
- [ ] Better error handling and logging
- [ ] Automated quality control checks
- [ ] Web-based result viewer

---

## Support

For questions or issues with this release:
1. Check the `README.md` for detailed documentation
2. Review the `QUICKSTART.md` for common setup issues
3. Examine log files in the `logs/` directory
4. Validate configuration with `./run_pipeline.sh validate`

---

**Release Date**: 2024  
**Compatibility**: Tested with bcftools 1.20, Python 3.9.5  
**Status**: Production-ready
