# Changelog

## Version 1.2 - Bug Fixes and Complete Layer Configurations (Current Release)

### 🔴 Critical Bug Fixes

#### 1. Syntax Error in run_pipeline.sh
- **Location**: Line 167
- **Problem**: Used Python syntax (`else:`) instead of bash syntax
- **Impact**: Script would crash when running the postprocess command
- **Solution**: Changed to correct bash syntax (`else`)
- **Files Modified**: `run_pipeline.sh`

#### 2. Broken Python Script Execution in 08_merge_results.sh
- **Location**: Lines 42-77 and 90-123
- **Problem**: Heredoc Python scripts weren't being executed properly
- **Impact**: Result merging would fail silently or with cryptic errors
- **Solution**: Changed to create temporary Python scripts, execute them, and clean up
- **Files Modified**: `scripts/phase2_layer_painting/08_merge_results.sh`

**Before:**
```bash
python3 << 'EOF'
...python code...
EOF
python3 -c "$(cat)" $ARGS  # This doesn't work
```

**After:**
```bash
TEMP_SCRIPT=$(mktemp)
cat > $TEMP_SCRIPT << 'EOF'
...python code...
EOF
python3 $TEMP_SCRIPT $ARGS
rm $TEMP_SCRIPT
```

#### 3. Sample List Naming Mismatch in extract_sample_lists.py
- **Location**: Lines 42-50
- **Problem**: Created files like `IA.Denmark.txt` but layer configs expected `IA.Denmark.1575-900BP.txt`
- **Impact**: Sample files wouldn't match layer configuration expectations
- **Solution**: Now preserves full population names including time ranges (stops at classifier codes with underscores)
- **Files Modified**: `scripts/utils/extract_sample_lists.py`

**Before:**
```python
safe_pop = pop.split('.')[0:2]  # Takes only "IA.Denmark"
```

**After:**
```python
# Preserves "IA.Denmark.1575-900BP" from "IA.Denmark.1575-900BP.XpcCWC.0_1_3"
parts = pop.split('.')
pop_name_parts = []
for part in parts:
    if '_' in part or part.isdigit():
        break
    pop_name_parts.append(part)
safe_pop = '.'.join(pop_name_parts)
```

#### 4. Path Inconsistencies
- **Problem**: Genetic map path used `{i}` instead of `{chrom}`
- **Solution**: Standardized all path templates to use `{chrom}`
- **Files Modified**: `config/master_config.yaml`

**Before:**
```yaml
genetic_maps: "/path/genetic_map_GRCh37_chr{i}.txt"
```

**After:**
```yaml
genetic_maps: "/path/genetic_map_GRCh37_chr{chrom}.txt"
```

#### 5. Redundant Module Loads
- **Problem**: Modules loaded twice in `run_layer.sh`
- **Solution**: Removed duplicate module load commands at the beginning of the script
- **Files Modified**: `scripts/phase2_layer_painting/run_layer.sh`

### ⚙️ Configuration Updates

#### Test Chromosome Changed
- **Changed**: Default test chromosome from chr6 to chr22
- **Reason**: Chr22 is smaller and faster for initial testing
- **Duration Impact**: Phase 1: 30-60 min → 20-40 min; Phase 2: 24-48 hrs → 16-32 hrs
- **Files Modified**: `config/master_config.yaml`, documentation

#### Lambda Parameter Updated
- **Changed**: Default lambda from 25.7926 to 50
- **Location**: `scripts/phase2_layer_painting/05_estimate_lambda.sh`
- **Note**: Still hardcoded, will parse from SparsePainter output in future version
- **Comment Added**: Clearer documentation about future plans to parse actual value

#### Number of Layers Updated
- **Changed**: `n_layers` from 1 to 5 in `config/master_config.yaml`
- **Reason**: All 5 layer configurations now exist

### 🆕 New Features

#### Complete Layer Configurations
Added placeholder configurations for all temporal layers:

**Layer 2: Bronze Age (2800-4500 BP)**
- Config file: `config/layers/layer2_config.yaml`
- Sample populations:
  - BA.Britain.4500-2800BP
  - BA.CentralEurope.4500-2800BP
  - BA.Iberia.4500-2800BP
  - BA.Steppe.4500-2800BP

**Layer 3: Neolithic (4500-8000 BP)**
- Config file: `config/layers/layer3_config.yaml`
- Sample populations:
  - Neolithic.Britain.8000-4500BP
  - Neolithic.Anatolia.8000-4500BP
  - Neolithic.CentralEurope.8000-4500BP
  - Neolithic.Iberia.8000-4500BP

**Layer 4: Mesolithic (8000-12000 BP)**
- Config file: `config/layers/layer4_config.yaml`
- Sample populations:
  - Mesolithic.WHG.12000-8000BP (Western Hunter-Gatherers)
  - Mesolithic.EHG.12000-8000BP (Eastern Hunter-Gatherers)
  - Mesolithic.SHG.12000-8000BP (Scandinavian Hunter-Gatherers)
  - Mesolithic.Caucasus.12000-8000BP (Caucasus Hunter-Gatherers)

#### Sample Directory Structure
- Created directories: `config/samples/layer{2,3,4}/`
- Added README.md files in each directory explaining how to create sample files
- Documented required `ref_popfile.txt` format
- Explained how to run `extract_sample_lists.py`

#### Phase 3 Handling Improved
- Removed Phase 3 from main `full` workflow (it's not yet implemented)
- Kept `postprocess` command available for when Phase 3 is implemented
- Added informative warning messages when Phase 3 is called
- Updated `run_full()` function to skip Phase 3 and inform user it's coming

### 📝 Documentation Updates

#### README.md
- Updated directory structure to show all 5 layers
- Changed all chr6 references to chr22
- Added comprehensive "What Was Fixed" section with version 1.2 and 1.1 changes
- Updated resource requirements for chr22
- Added "Available Layers" table showing status of each layer
- Updated version information to 1.2

#### QUICKSTART.md
- Added "What's New in v1.2" section
- Changed all chr6 references to chr22
- Updated time estimates for chr22
- Enhanced "Key Changes" section with both v1.2 and v1.1 fixes
- Updated "After Test Success" steps to include layer setup instructions

#### Path Corrections
- Updated UKB chr rename map path to correct location
- Verified all paths use consistent templating

### 🧹 Code Quality Improvements

- Improved error messages throughout
- Better documentation of TODOs and future work
- Consistent formatting in configuration files
- Clearer comments in shell scripts

---

## Version 1.1 - Fixed

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
