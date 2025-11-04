# Layer 2 Sample Files (Bronze Age)

This directory should contain sample list files for Bronze Age reference populations.

## Required Files

You need to create a `ref_popfile.txt` with the format:
```
sample_id1    population_name
sample_id2    population_name
...
```

Where population names match those in `layer2_config.yaml`:
- BA.Britain.4500-2800BP
- BA.CentralEurope.4500-2800BP
- BA.Iberia.4500-2800BP
- BA.Steppe.4500-2800BP

## How to Generate

Run the extraction script after creating ref_popfile.txt:
```bash
python3 ../../scripts/utils/extract_sample_lists.py ref_popfile.txt ./
```

This will create individual files for each population.
