# HLA Epitope Prediction - Performance Tuning Guide

## Quick Performance Reference

### Conservative (Default) - Recommended for most users
**~10x speedup over comprehensive mode**

```bash
python3 hla_epitope_predictor.py \
    --email user@example.com \
    --organism "Organism[Organism]" \
    --class-i-lengths 9 \
    --class-i-step 3 \
    --class-ii-step 3
```

**Coverage**: ~70-80% of epitopes
**Use case**: Production analysis, vaccine design

---

### Comprehensive - Maximum sensitivity
**Baseline speed (no optimization)**

```bash
python3 hla_epitope_predictor.py \
    --email user@example.com \
    --organism "Organism[Organism]" \
    --class-i-lengths 8 9 10 11 \
    --class-i-step 1 \
    --class-ii-step 1
```

**Coverage**: ~100% of epitopes
**Use case**: Research, complete epitope mapping

---

### Aggressive - Fast screening
**~30x speedup over comprehensive mode**

```bash
python3 hla_epitope_predictor.py \
    --email user@example.com \
    --organism "Organism[Organism]" \
    --class-i-lengths 9 \
    --class-i-step 5 \
    --class-ii-step 5
```

**Coverage**: ~50-60% of epitopes
**Use case**: Rapid pathogen screening, initial assessment

---

### Ultra-Fast - Non-overlapping peptides
**~100x speedup over comprehensive mode**

```bash
python3 hla_epitope_predictor.py \
    --email user@example.com \
    --organism "Organism[Organism]" \
    --class-i-lengths 9 \
    --class-i-step 9 \
    --class-ii-step 15
```

**Coverage**: ~30-40% of epitopes
**Use case**: Very large proteomes, exploratory analysis

---

## Performance Impact Table

| Configuration | Class I Peptides/Protein | Class II Peptides/Protein | Total Speedup | Epitope Coverage |
|--------------|-------------------------|--------------------------|---------------|------------------|
| **Comprehensive** (8-11mers, step=1) | 1,166 | 286 | 1x | ~100% |
| **Conservative** (9mers, step=3) | 97 | 95 | ~10-12x | ~70-80% |
| **Aggressive** (9mers, step=5) | 59 | 58 | ~20-25x | ~50-60% |
| **Ultra-Fast** (9mers, step=9/15) | 33 | 20 | ~70-100x | ~30-40% |

*Assumes 300 amino acid average protein length*

---

## Parameter Guide

### Class I Peptide Lengths (`--class-i-lengths`)

**Options**: 8, 9, 10, 11 (or any combination)

| Configuration | Rationale | Coverage |
|--------------|-----------|----------|
| `9` only | Most common epitope length (~60-70% of epitopes) | Good |
| `9 10` | Covers ~85% of epitopes | Better |
| `8 9 10 11` | All common lengths | Best |

**Recommendation**: Use `9` for most analyses. Add other lengths if specifically studying rare epitopes.

---

### Step Size (`--class-i-step`, `--class-ii-step`)

**Definition**: Number of amino acids to shift the sliding window

| Step | Overlap | Peptides Generated | Coverage | Use Case |
|------|---------|-------------------|----------|----------|
| 1 | High (length-1 AA) | Maximum | ~100% | Complete mapping |
| 3 | Medium (length-3 AA) | ~3x fewer | ~70-80% | Standard analysis |
| 5 | Low (length-5 AA) | ~5x fewer | ~50-60% | Fast screening |
| 9+ | None (or minimal) | ~9-15x fewer | ~30-40% | Exploratory |

**Recommendation**:
- **Step=3** for most analyses (good balance)
- **Step=1** if you need comprehensive epitope mapping
- **Step=5+** for very large proteomes or initial screening

---

## Memory Considerations

Adjust `--n-jobs` based on proteome size and available RAM:

| Proteome Size | n_jobs (conservative) | n_jobs (comprehensive) | Expected Memory |
|---------------|----------------------|------------------------|-----------------|
| Small (<500 proteins) | 40-60 | 20-30 | ~40-60 GB |
| Medium (500-2000) | 30-40 | 15-20 | ~30-40 GB |
| Large (>2000) | 20-30 | 10-15 | ~20-30 GB |

**Rule of thumb**: Each parallel job needs ~1 GB RAM

---

## Example Workflows

### Workflow 1: Standard Vaccine Design
```bash
# Conservative mode with all alleles
python3 hla_epitope_predictor.py \
    --email user@example.com \
    --organism "Variola virus[Organism]" \
    --class-i-lengths 9 \
    --class-i-step 3 \
    --class-ii-step 3 \
    --n-jobs 40
```

### Workflow 2: Large Proteome Screening
```bash
# Aggressive mode for P. vivax (~5000 proteins)
python3 hla_epitope_predictor.py \
    --email user@example.com \
    --organism "Plasmodium vivax[Organism]" \
    --class-i-lengths 9 \
    --class-i-step 5 \
    --class-ii-step 5 \
    --n-jobs 20
```

### Workflow 3: Comprehensive Research Analysis
```bash
# Full coverage, high quality
python3 hla_epitope_predictor.py \
    --email user@example.com \
    --organism "HIV[Organism]" \
    --class-i-lengths 8 9 10 11 \
    --class-i-step 1 \
    --class-ii-step 1 \
    --n-jobs 30
```

### Workflow 4: Multi-Pathogen Screening
```bash
# Ultra-fast mode for comparing many pathogens
for org in "Measles" "Mumps" "Rubella"; do
    python3 hla_epitope_predictor.py \
        --email user@example.com \
        --organism "${org}[Organism]" \
        --class-i-step 9 \
        --class-ii-step 15 \
        --n-jobs 40 \
        --output-dir "./${org}_analysis"
done
```

---

## Changing Settings Mid-Analysis

The incremental processing feature allows you to adjust parameters and re-run:

### Example: Start fast, then refine

```bash
# Step 1: Quick screening with ultra-fast mode
python3 hla_epitope_predictor.py \
    --organism "Organism[Organism]" \
    --class-i-step 9 \
    --class-ii-step 15 \
    --output-dir ./screening

# Step 2: Identify interesting proteins from results
# (e.g., proteins with many hits)

# Step 3: Re-run with comprehensive mode on full proteome
# Since alleles are already completed, this adds new peptides
python3 hla_epitope_predictor.py \
    --organism "Organism[Organism]" \
    --class-i-lengths 8 9 10 11 \
    --class-i-step 1 \
    --class-ii-step 1 \
    --output-dir ./comprehensive
```

**Note**: Changing peptide parameters requires a new analysis (different peptides). Incremental processing only applies to adding new alleles with the same peptide settings.

---

## Performance vs. Biology Trade-offs

### When to use Conservative (step=3)?
✅ Standard vaccine design
✅ Epitope discovery for known pathogens
✅ Population-level coverage analysis
✅ Most research applications

### When to use Comprehensive (step=1)?
✅ Novel pathogen analysis
✅ Studying epitope structure/binding
✅ Identifying all possible epitopes
✅ Publication-quality complete mapping

### When to use Aggressive (step=5+)?
✅ Initial pathogen screening
✅ Very large proteomes (>5000 proteins)
✅ Comparative analyses across many organisms
✅ Exploratory research

---

## Validation

To assess if you're missing important epitopes, compare conservative vs. comprehensive:

```bash
# Run both modes
python3 hla_epitope_predictor.py --organism "X" --class-i-step 3 --output-dir ./conservative
python3 hla_epitope_predictor.py --organism "X" --class-i-step 1 --output-dir ./comprehensive

# Compare strong binders
python3 analyze_predictions.py --input ./conservative/X_hla_predictions.csv
python3 analyze_predictions.py --input ./comprehensive/X_hla_predictions.csv

# If conservative captures >70% of strong binders, it's adequate for most purposes
```

---

## Quick Decision Tree

```
Are you working with a large proteome (>2000 proteins)?
├─ YES → Use aggressive mode (step=5) or conservative (step=3) with low n-jobs (20)
└─ NO → Continue

Is this a novel/unstudied pathogen?
├─ YES → Use comprehensive mode (step=1) for complete mapping
└─ NO → Continue

Do you need publication-quality complete epitope maps?
├─ YES → Use comprehensive mode (step=1, all lengths)
└─ NO → Continue

Is this for vaccine design or standard analysis?
└─ YES → Use conservative mode (step=3, 9-mers) ← DEFAULT
```

---

## Monitoring Performance

Watch the log output for peptide counts:

```
INFO: Generating 9-mer peptides (step=3)...
INFO: Generated 12,450 peptides               ← Conservative: ~12k
INFO: Running predictions for 45 alleles
```

Compare to comprehensive mode:
```
INFO: Generating peptides of lengths: [8, 9, 10, 11] (step=1)...
INFO: Generated 145,890 peptides              ← Comprehensive: ~146k
INFO: Running predictions for 45 alleles
```

**If peptide count is >100k per analysis**, consider using more aggressive settings.
