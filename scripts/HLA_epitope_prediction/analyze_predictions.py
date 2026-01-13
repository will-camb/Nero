#!/usr/bin/env python3
"""
Analyze HLA epitope predictions to identify strong binders.

This script aggregates predictions by allele and counts unique
proteins with strong binding affinity.
"""

import sys
import argparse
from pathlib import Path
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def analyze_predictions(
    input_file: Path,
    output_dir: Path,
    class_i_rank_threshold: float = 0.5,
    class_ii_rank_threshold: float = 1.0,
    organism_name: str = None
):
    """
    Analyze HLA predictions and generate summary statistics.

    Args:
        input_file: Path to combined predictions CSV
        output_dir: Directory for output files
        class_i_rank_threshold: %Rank threshold for Class I strong binders
        class_ii_rank_threshold: %Rank threshold for Class II strong binders
        organism_name: Optional organism name for output files
    """
    logger.info(f"Loading predictions from {input_file}")

    try:
        df = pd.read_csv(input_file)
    except FileNotFoundError:
        logger.error(f"File not found: {input_file}")
        return False
    except pd.errors.EmptyDataError:
        logger.error(f"File is empty: {input_file}")
        return False

    # Validate required columns
    required_cols = ['hla_class', 'hla_allele', 'protein_id', 'rank']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        logger.error(f"Missing required columns: {missing_cols}")
        logger.info(f"Available columns: {df.columns.tolist()}")
        return False

    logger.info(f"Loaded {len(df)} predictions")
    logger.info(f"Classes found: {df['hla_class'].unique()}")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine organism name from input file if not provided
    if organism_name is None:
        organism_name = input_file.stem.replace('_hla_predictions', '')

    # Analyze Class I
    logger.info("\n=== Class I Analysis ===")
    df_i = df[df['hla_class'] == 'I'].copy()

    if len(df_i) > 0:
        logger.info(f"Total Class I predictions: {len(df_i)}")

        # Filter strong binders
        df_i_strong = df_i[df_i['rank'] < class_i_rank_threshold]
        logger.info(f"Strong binders (rank < {class_i_rank_threshold}): {len(df_i_strong)}")

        if len(df_i_strong) > 0:
            # Count unique proteins per allele
            allele_counts = (df_i_strong
                           .groupby('hla_allele')['protein_id']
                           .nunique()
                           .sort_values(ascending=False))

            output_file = output_dir / f"{organism_name}_class_i_strong_binders.tsv"
            allele_counts.to_csv(output_file, sep='\t', header=['unique_proteins'])
            logger.info(f"Saved Class I summary to {output_file}")

            # Print top 10
            logger.info("\nTop 10 Class I alleles by unique proteins:")
            for allele, count in allele_counts.head(10).items():
                logger.info(f"  {allele}: {count} proteins")
        else:
            logger.warning("No strong Class I binders found")
    else:
        logger.warning("No Class I predictions found")

    # Analyze Class II
    logger.info("\n=== Class II Analysis ===")
    df_ii = df[df['hla_class'] == 'II'].copy()

    if len(df_ii) > 0:
        logger.info(f"Total Class II predictions: {len(df_ii)}")

        # Filter strong binders
        df_ii_strong = df_ii[df_ii['rank'] < class_ii_rank_threshold]
        logger.info(f"Strong binders (rank < {class_ii_rank_threshold}): {len(df_ii_strong)}")

        if len(df_ii_strong) > 0:
            # Count unique proteins per allele
            allele_counts = (df_ii_strong
                           .groupby('hla_allele')['protein_id']
                           .nunique()
                           .sort_values(ascending=False))

            output_file = output_dir / f"{organism_name}_class_ii_strong_binders.tsv"
            allele_counts.to_csv(output_file, sep='\t', header=['unique_proteins'])
            logger.info(f"Saved Class II summary to {output_file}")

            # Print top 10
            logger.info("\nTop 10 Class II alleles by unique proteins:")
            for allele, count in allele_counts.head(10).items():
                logger.info(f"  {allele}: {count} proteins")
        else:
            logger.warning("No strong Class II binders found")
    else:
        logger.warning("No Class II predictions found")

    # Generate combined summary
    logger.info("\n=== Combined Summary ===")

    summary_data = []

    if len(df_i) > 0:
        df_i_strong = df_i[df_i['rank'] < class_i_rank_threshold]
        summary_data.append({
            'class': 'I',
            'total_predictions': len(df_i),
            'strong_binders': len(df_i_strong),
            'unique_proteins': df_i['protein_id'].nunique(),
            'proteins_with_strong_binders': df_i_strong['protein_id'].nunique() if len(df_i_strong) > 0 else 0,
            'n_alleles': df_i['hla_allele'].nunique(),
            'rank_threshold': class_i_rank_threshold
        })

    if len(df_ii) > 0:
        df_ii_strong = df_ii[df_ii['rank'] < class_ii_rank_threshold]
        summary_data.append({
            'class': 'II',
            'total_predictions': len(df_ii),
            'strong_binders': len(df_ii_strong),
            'unique_proteins': df_ii['protein_id'].nunique(),
            'proteins_with_strong_binders': df_ii_strong['protein_id'].nunique() if len(df_ii_strong) > 0 else 0,
            'n_alleles': df_ii['hla_allele'].nunique(),
            'rank_threshold': class_ii_rank_threshold
        })

    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_file = output_dir / f"{organism_name}_summary.tsv"
        summary_df.to_csv(summary_file, sep='\t', index=False)
        logger.info(f"\nSaved summary to {summary_file}")
        logger.info("\n" + summary_df.to_string(index=False))

    return True


def batch_analyze(analysis_dirs: list, output_dir: Path,
                 class_i_threshold: float = 0.5,
                 class_ii_threshold: float = 1.0):
    """
    Analyze multiple organism predictions in batch.

    Args:
        analysis_dirs: List of analysis directories
        output_dir: Output directory for aggregated results
        class_i_threshold: Class I rank threshold
        class_ii_threshold: Class II rank threshold
    """
    logger.info(f"Batch analysis of {len(analysis_dirs)} organisms")
    output_dir.mkdir(parents=True, exist_ok=True)

    for analysis_dir in analysis_dirs:
        analysis_path = Path(analysis_dir)

        if not analysis_path.exists():
            logger.warning(f"Directory not found: {analysis_path}")
            continue

        # Find predictions file
        predictions_files = list(analysis_path.glob("*_hla_predictions.csv"))

        if not predictions_files:
            # Try old naming convention
            predictions_files = list(analysis_path.glob("*_hla_binding_results.csv"))

        if not predictions_files:
            logger.warning(f"No predictions file found in {analysis_path}")
            continue

        predictions_file = predictions_files[0]
        organism_name = analysis_path.name.replace('_analysis', '')

        logger.info(f"\n{'='*60}")
        logger.info(f"Analyzing: {organism_name}")
        logger.info(f"{'='*60}")

        analyze_predictions(
            predictions_file,
            output_dir,
            class_i_threshold,
            class_ii_threshold,
            organism_name
        )


def main():
    parser = argparse.ArgumentParser(
        description="Analyze HLA epitope predictions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze single organism
  python analyze_predictions.py --input smallpox_analysis/smallpox_hla_predictions.csv

  # Batch analyze multiple organisms
  python analyze_predictions.py --batch \\
      smallpox_analysis HIV_analysis plague_analysis

  # Custom thresholds
  python analyze_predictions.py --input results.csv \\
      --class-i-threshold 0.05 --class-ii-threshold 0.1
        """
    )

    # Input mode
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '--input', type=Path,
        help='Single predictions CSV file'
    )
    input_group.add_argument(
        '--batch', nargs='+',
        help='Multiple analysis directories to process'
    )

    # Thresholds
    parser.add_argument(
        '--class-i-threshold', type=float, default=0.5,
        help='%%Rank threshold for Class I strong binders (default: 0.5)'
    )
    parser.add_argument(
        '--class-ii-threshold', type=float, default=1.0,
        help='%%Rank threshold for Class II strong binders (default: 1.0)'
    )

    # Output
    parser.add_argument(
        '--output-dir', type=Path, default='./analyzed_results',
        help='Output directory (default: ./analyzed_results)'
    )

    args = parser.parse_args()

    try:
        if args.input:
            # Single file analysis
            success = analyze_predictions(
                args.input,
                args.output_dir,
                args.class_i_threshold,
                args.class_ii_threshold
            )
            sys.exit(0 if success else 1)

        elif args.batch:
            # Batch analysis
            batch_analyze(
                args.batch,
                args.output_dir,
                args.class_i_threshold,
                args.class_ii_threshold
            )
            sys.exit(0)

    except Exception as e:
        logger.error(f"Analysis failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
