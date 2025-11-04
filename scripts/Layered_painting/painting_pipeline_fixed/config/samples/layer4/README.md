# Layer 4 Sample Files (Mesolithic)

This directory should contain sample list files for Mesolithic reference populations.

## Required Files

You need to create a `ref_popfile.txt` with the format:
```
sample_id1    population_name
sample_id2    population_name
...
```

Where population names match those in `layer4_config.yaml`:
- Mesolithic.WHG.12000-8000BP
- Mesolithic.EHG.12000-8000BP
- Mesolithic.SHG.12000-8000BP
- Mesolithic.Caucasus.12000-8000BP

## How to Generate

Run the extraction script after creating ref_popfile.txt:
```bash
python3 ../../scripts/utils/extract_sample_lists.py ref_popfile.txt ./
```

This will create individual files for each population.
