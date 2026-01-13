# HLA Epitope Prediction Pipeline - Improvements

## Summary of Changes

This refactored version addresses memory issues, enables incremental processing, and improves code quality while maintaining the same core functionality.

---

## Problems Fixed

### 1. **Memory and File Completion Issues**

**Root causes identified:**

- **Excessive parallelization**: Running 80 parallel jobs with 234 alleles creates:
  - 80+ concurrent netMHC processes, each loading ~500MB-1GB models
  - Hundreds of temporary files simultaneously
  - File descriptor exhaustion

- **No cleanup**: Temporary directories accumulated across runs

- **All-or-nothing processing**: If one allele failed, entire analysis had to restart

**Solutions implemented:**

- **Incremental processing**: Each allele's results saved separately and tracked
- **Resume capability**: Can resume from where analysis stopped
- **Lower default parallelization**: Reduced from 80 to 20-40 jobs depending on proteome size
- **Progress tracking**: JSON file tracks completed alleles (`completed_alleles.json`)
- **Better error handling**: Partial failures don't invalidate entire analysis

### 2. **Inability to Add New Alleles**

**Problem**: To analyze new alleles, you had to re-run all 234 alleles (~hours of compute)

**Solution**:
- Results stored per-allele in `allele_results/` directory
- Script automatically skips already-completed alleles
- Adding new alleles only processes the new ones

**Example workflow:**
```bash
# Initial run with 45 Class I alleles
python3 hla_epitope_predictor.py --organism "Variola virus" --email user@example.com

# Later, add 5 new alleles to alleles_class_i_default.txt
# Re-run: only the 5 new alleles are processed
python3 hla_epitope_predictor.py --organism "Variola virus" --email user@example.com
```

### 3. **Code Quality Issues**

**Problems**:
- Hardcoded alleles in Python code
- Monolithic 1100-line class
- No separation of concerns
- Poor naming conventions
- Commented dead code
- Magic numbers throughout

**Solutions**:
- **Modular design**: Separate classes for each responsibility:
  - `ProteomeDownloader`: NCBI data retrieval
  - `SignalPFilter`: SignalP filtering
  - `PeptideGenerator`: Peptide generation
  - `HLAPredictor`: Base class for predictions
  - `NetMHCpanPredictor`, `NetMHCIIpanPredictor`: Specific predictors
  - `ResultParser`: Parse prediction outputs
  - `HLAEpitopePipeline`: Orchestrates full workflow

- **External configuration**: Alleles in text files (`alleles_class_i_default.txt`, `alleles_class_ii_default.txt`)

- **Better naming**:
  - `analyze_organism()` → `run_full_analysis()`
  - `get_netmhciipan_alleles()` → `load_alleles_from_file()`
  - Clear, descriptive variable names

- **Type hints**: All functions have type annotations

- **Documentation**: Docstrings for all classes and methods

---

## New Features

### 1. **Configurable Alleles**

Alleles now specified in simple text files:

```text
# alleles_class_i_default.txt
HLA-A01:01
HLA-A02:01
# Add more...
```

Can use custom allele files:
```bash
python3 hla_epitope_predictor.py \
    --class-i-alleles my_alleles.txt \
    --class-ii-alleles my_class_ii.txt
```

### 2. **Incremental Processing**

Results cached per-allele in `output_dir/class_i/allele_results/` and `output_dir/class_ii/allele_results/`.

Progress tracked in `completed_alleles.json`:
```json
[
  "HLA-A01:01",
  "HLA-A02:01",
  "HLA-B07:02"
]
```

### 3. **Better Memory Management**

- Reduced default parallelization (20-40 instead of 80)
- Per-allele result files instead of monolithic temp directory
- Clearer guidelines for adjusting `--n-jobs` based on proteome size

### 4. **Improved Error Handling**

- Continues on individual allele failures
- Marks successful alleles as complete
- Reports which alleles failed at the end
- Can resume from failures

### 5. **Configuration Object**

Clean configuration via `AnalysisConfig` dataclass:
```python
config = AnalysisConfig(
    email="user@example.com",
    organism_query="Variola virus[Organism]",
    output_dir="./analysis",
    organism_name="smallpox",
    netmhcpan_path="/path/to/netMHCpan",
    netmhciipan_path="/path/to/netMHCIIpan",
    signalp_model_dir="/path/to/signalp/models",
    n_parallel_jobs=20
)
```

---

## Code Structure Comparison

### Before (Monolithic)
```
pathogen_hla_binding_class_i_ii.py (1119 lines)
├── PathogenHLAAnalyzer class (does everything)
│   ├── download_proteome()
│   ├── run_signalp()
│   ├── filter_secreted_proteins()
│   ├── generate_peptides()
│   ├── generate_class_i_peptides()
│   ├── get_netmhciipan_alleles() [189 hardcoded alleles]
│   ├── get_netmhcpan_alleles() [45 hardcoded alleles]
│   ├── run_netmhciipan_batch()
│   ├── run_netmhcpan_batch()
│   ├── parse_netmhciipan_results()
│   ├── parse_netmhcpan_results()
│   └── analyze_organism()
└── main()
```

### After (Modular)
```
hla_epitope_predictor.py (950 lines, better organized)
├── AnalysisConfig (dataclass)
├── ProteomeDownloader
│   └── download()
├── SignalPFilter
│   ├── run()
│   └── filter_fasta()
├── PeptideGenerator
│   ├── generate_peptides()
│   └── generate_variable_length_peptides()
├── HLAPredictor (base class)
│   ├── get_completed_alleles()
│   ├── mark_allele_complete()
│   └── run_predictions()
├── NetMHCpanPredictor (Class I)
├── NetMHCIIpanPredictor (Class II)
├── ResultParser
│   └── parse_netmhc_results()
├── HLAEpitopePipeline
│   ├── run_class_i_analysis()
│   ├── run_class_ii_analysis()
│   └── run_full_analysis()
└── main()

alleles_class_i_default.txt (45 alleles)
alleles_class_ii_default.txt (189 alleles)
```

---

## Performance Improvements

### Memory Usage
- **Before**: 80 jobs × 1GB per job = ~80GB peak memory
- **After**: 20-40 jobs × 1GB per job = ~20-40GB peak memory

### Resilience
- **Before**: One failure = restart entire analysis (hours wasted)
- **After**: Resume from last completed allele (minutes)

### Adding Alleles
- **Before**: Re-run all 234 alleles = ~2-6 hours
- **After**: Run only new alleles = ~5-15 minutes per allele

---

## Usage Examples

### Basic Usage (with defaults)
```bash
python3 hla_epitope_predictor.py \
    --email user@example.com \
    --organism "Variola virus[Organism] AND RefSeq[Filter]" \
    --output-dir ./smallpox_analysis
```

### Custom Alleles
```bash
python3 hla_epitope_predictor.py \
    --email user@example.com \
    --organism "HIV[Organism]" \
    --class-i-alleles my_alleles.txt \
    --class-ii-alleles my_class_ii_alleles.txt
```

### Adding New Alleles (Incremental)
```bash
# Add new alleles to alleles_class_i_default.txt
echo "HLA-A24:03" >> alleles_class_i_default.txt
echo "HLA-B35:01" >> alleles_class_i_default.txt

# Re-run: only new alleles processed
python3 hla_epitope_predictor.py \
    --email user@example.com \
    --organism "Variola virus[Organism]" \
    --output-dir ./smallpox_analysis  # Same directory
```

### Adjust Parallelization
```bash
# Large proteome (e.g., P. vivax ~5000 proteins)
python3 hla_epitope_predictor.py \
    --organism "Plasmodium vivax[Organism]" \
    --n-jobs 20

# Small proteome (e.g., HIV ~20 proteins)
python3 hla_epitope_predictor.py \
    --organism "HIV[Organism]" \
    --n-jobs 60
```

### Batch Processing
```bash
# Use improved batch runner
chmod +x run_pathogen_analysis.sh
./run_pathogen_analysis.sh
```

---

## Output Structure

```
output_dir/
├── organism_name_proteome.fasta          # Downloaded proteome
├── organism_name_secreted.fasta          # Filtered for Class II
├── signalp/                               # SignalP results
│   ├── proteins_short_ids.fasta
│   └── prediction_results.txt
├── class_i/                               # Class I analysis
│   ├── peptides.txt                       # Generated peptides
│   ├── completed_alleles.json             # Progress tracking
│   ├── allele_results/                    # Individual allele results
│   │   ├── HLA-A01_01.txt
│   │   ├── HLA-A02_01.txt
│   │   └── ...
│   └── combined_results.txt               # All alleles combined
├── class_ii/                              # Class II analysis
│   ├── peptides.txt
│   ├── completed_alleles.json
│   ├── allele_results/
│   │   ├── DRB1_0101.txt
│   │   └── ...
│   └── combined_results.txt
└── organism_name_hla_predictions.csv      # Final combined results
```

---

## Analysis Recommendations

### 1. **Avoid Memory Issues**

**Guideline**: `n_jobs` × 1GB < Available RAM

| Proteome Size | Proteins | Recommended n_jobs |
|---------------|----------|-------------------|
| Small         | < 500    | 40-60             |
| Medium        | 500-2000 | 30-40             |
| Large         | > 2000   | 20-30             |

### 2. **Monitor Progress**

Check completed alleles:
```bash
cat output_dir/class_i/completed_alleles.json
cat output_dir/class_ii/completed_alleles.json
```

### 3. **Resume Failed Runs**

Simply re-run the same command - completed alleles are skipped automatically.

---

## Remaining Analysis Considerations

### 1. **SignalP Organism Parameter**

**Current**: Uses `organism_type="other"` for all pathogens

**Consideration**: SignalP6 has organism-specific models:
- `gram+` for Gram-positive bacteria
- `gram-` for Gram-negative bacteria
- `eukarya` for eukaryotic parasites

**Recommendation**: Add organism-specific logic:
```python
if "Plasmodium" in organism_name or "vivax" in organism_name:
    organism_type = "eukarya"
elif "Mycobacterium" in organism_name:
    organism_type = "gram+"
elif "Yersinia" in organism_name or "Leptospira" in organism_name:
    organism_type = "gram-"
else:
    organism_type = "other"
```

### 2. **Class II Filtering Rationale**

**Current approach**:
- Class I: ALL proteins (intracellular processing)
- Class II: SignalP-filtered only (extracellular processing)

**Consideration**: This is biologically appropriate BUT:
- Some Class II epitopes come from endocytosed/phagocytosed proteins
- Pathogens that live intracellularly (e.g., M. tuberculosis) may present non-secreted proteins via Class II

**Recommendation**: Consider analyzing all proteins for Class II as well, or provide option.

### 3. **Peptide Length Flexibility**

**Current**: Fixed 15-mers for Class II

**Consideration**: netMHCIIpan can handle 15-19mers. Variable length may capture more epitopes.

**Recommendation**: Add option for Class II peptide lengths:
```python
class_ii_peptide_lengths: List[int] = [15, 16, 17, 18, 19]
```

### 4. **Binding Threshold Selection**

**Current analysis script** (analyse_results.py):
- Class I: `rank < 0.05`
- Class II: `rank < 0.1`

**Consideration**: These are very stringent thresholds. Common thresholds:
- Strong binders: rank < 0.5
- Weak binders: rank < 2.0

**Recommendation**: Document threshold rationale and provide options.

### 5. **Statistical Analysis Gaps**

**Missing**:
- Population coverage calculation (allele frequency × binding)
- Epitope conservation analysis
- Protein immunogenicity scoring
- Comparison across pathogens

**Recommendation**: Add downstream analysis module for these metrics.

---

## Migration Guide

To switch from old to new code:

### 1. **Backup existing results**
```bash
cp -r existing_analysis existing_analysis_backup
```

### 2. **Create allele files** (or use defaults)
```bash
# Use provided defaults
cp alleles_class_i_default.txt my_alleles_i.txt
cp alleles_class_ii_default.txt my_alleles_ii.txt
```

### 3. **Run new script**
```bash
python3 hla_epitope_predictor.py \
    --email your@email.com \
    --organism "Organism[Organism] AND RefSeq[Filter]" \
    --output-dir ./new_analysis \
    --class-i-alleles my_alleles_i.txt \
    --class-ii-alleles my_alleles_ii.txt
```

### 4. **Verify results match**
The prediction columns should be identical. Output filename changed from `_hla_binding_results.csv` to `_hla_predictions.csv`.

---

## Future Enhancements

1. **Distributed computing**: Spread allele predictions across multiple nodes
2. **Database backend**: Store results in SQLite/PostgreSQL for complex queries
3. **Web interface**: GUI for configuring and monitoring analyses
4. **Automatic allele validation**: Check if alleles are supported by prediction tools
5. **Peptide clustering**: Group similar peptides to reduce redundancy
6. **MHC Class I pathway modeling**: Predict proteasomal cleavage, TAP transport
7. **Immunogenicity scoring**: Integrate T-cell recognition models

---

## Questions?

For issues or questions:
1. Check the log output for specific error messages
2. Review `completed_alleles.json` to see progress
3. Verify tool paths are correct
4. Ensure sufficient memory for `n_jobs` setting
