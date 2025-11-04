# Layer 3 Sample Files (Neolithic)

This directory should contain sample list files for Neolithic reference populations.

## Required Files

You need to create a `ref_popfile.txt` with the format:
```
sample_id1    population_name
sample_id2    population_name
...
```

Where population names match those in `layer3_config.yaml`:
- Neolithic.Britain.8000-4500BP
- Neolithic.Anatolia.8000-4500BP
- Neolithic.CentralEurope.8000-4500BP
- Neolithic.Iberia.8000-4500BP

## How to Generate

Run the extraction script after creating ref_popfile.txt:
```bash
python3 ../../scripts/utils/extract_sample_lists.py ref_popfile.txt ./
```

This will create individual files for each population.
