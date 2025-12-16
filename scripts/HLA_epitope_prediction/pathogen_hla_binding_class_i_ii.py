#!/usr/bin/env python3
"""
Pathogen HLA Class I & II Binding Analysis Pipeline


Author: Generated for systematic pathogen immunogenomics analysis
"""

import os
import sys
import subprocess
import tempfile
import pandas as pd
from pathlib import Path
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse
import logging
from typing import List, Dict, Tuple
import re

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class PathogenHLAAnalyzer:
    def __init__(self, email: str, netmhciipan_path: str = "/home/jdn321/software/netMHCIIpan-4.3/netMHCIIpan", 
                 netmhcpan_path: str = "/home/jdn321/software/netMHCpan-4.2/netMHCpan"):
        """
        Initialize the analyzer with required parameters.
        
        Args:
            email: Email for NCBI Entrez API
            netmhciipan_path: Path to netMHCIIpan executable (Class II)
            netmhcpan_path: Path to netMHCpan executable (Class I)
        """
        self.email = email
        self.netmhciipan_path = netmhciipan_path
        self.netmhcpan_path = netmhcpan_path
        Entrez.email = email
        
        # Verify both tools are available
        if not os.path.exists(self.netmhciipan_path):
            raise FileNotFoundError(f"netMHCIIpan not found at {self.netmhciipan_path}")
        if not os.path.exists(self.netmhcpan_path):
            raise FileNotFoundError(f"netMHCpan not found at {self.netmhcpan_path}")
    
    def download_proteome(self, organism_query: str, output_file: str) -> str:
        """
        Download proteome from NCBI for specified organism.
        
        Args:
            organism_query: NCBI search query (e.g., "Variola virus[Organism]")
            output_file: Path to save the proteome FASTA file
            
        Returns:
            Path to downloaded proteome file
        """
        logger.info(f"Searching for proteins: {organism_query}")
        
        # Search for protein sequences
        search_handle = Entrez.esearch(
            db="protein", 
            term=organism_query,
            retmax=15000  # Adjust if needed
        )
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        protein_ids = search_results["IdList"]
        logger.info(f"Found {len(protein_ids)} protein sequences")
        
        if not protein_ids:
            raise ValueError(f"No proteins found for query: {organism_query}")
        
        # Download sequences in batches
        batch_size = 100
        all_records = []
        
        for i in range(0, len(protein_ids), batch_size):
            batch_ids = protein_ids[i:i+batch_size]
            logger.info(f"Downloading batch {i//batch_size + 1}/{(len(protein_ids)-1)//batch_size + 1}")
            
            fetch_handle = Entrez.efetch(
                db="protein",
                id=batch_ids,
                rettype="fasta",
                retmode="text"
            )
            
            batch_records = list(SeqIO.parse(fetch_handle, "fasta"))
            all_records.extend(batch_records)
            fetch_handle.close()
        
        # Save to file
        with open(output_file, "w") as f:
            SeqIO.write(all_records, f, "fasta")
        
        logger.info(f"Downloaded {len(all_records)} sequences to {output_file}")
        return output_file
    
    def create_short_id_fasta(self, input_fasta: str, output_fasta: str) -> Dict[str, str]:
        """
        Create a FASTA file with shortened IDs and return mapping.
        
        Returns:
            Dictionary mapping short_id -> original_id
        """
        id_mapping = {}
        short_records = []
        
        with open(input_fasta, 'r') as f:
            for i, record in enumerate(SeqIO.parse(f, "fasta")):
                short_id = f"protein_{i+1:06d}"  # e.g., protein_000001
                id_mapping[short_id] = record.id
                
                # Create new record with short ID
                short_record = SeqRecord(
                    record.seq,
                    id=short_id,
                    description=""  # Remove long description too
                )
                short_records.append(short_record)
        
        with open(output_fasta, 'w') as f:
            SeqIO.write(short_records, f, "fasta")
        
        logger.info(f"Created shortened FASTA with {len(short_records)} sequences")
        return id_mapping
    
    def run_signalp(self, input_fasta: str, output_dir: str) -> Tuple[str, List[str]]:
        """
        Run SignalP to identify secreted proteins.
        """
        logger.info("Running SignalP analysis...")
        
        # Create shortened FASTA to avoid filename issues
        short_fasta = os.path.join(output_dir, "proteins_short_ids.fasta")
        id_mapping = self.create_short_id_fasta(input_fasta, short_fasta)
        
        # Run SignalP on shortened file
        cmd = [
            "signalp6",
            "--fastafile", short_fasta,  # Use shortened file
            "--output_dir", output_dir,
            "--format", "none",
            "--mode", "fast",
            "--organism", "eukarya",
            "--model_dir", "/home/jdn321/software/signalp6_fast/signalp-6-package/models/"
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.info("SignalP completed successfully")
        except subprocess.CalledProcessError as e:
            logger.error(f"SignalP failed: {e}")
            raise
        
        # Parse results and map back to original IDs
        secreted_ids = []
        prediction_file = os.path.join(output_dir, "prediction_results.txt")
        
        if os.path.exists(prediction_file):
            with open(prediction_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 2 and parts[1] == "SP":
                        short_id = parts[0]
                        original_id = id_mapping.get(short_id, short_id)
                        secreted_ids.append(original_id)  # Use original ID
        
        logger.info(f"Found {len(secreted_ids)} secreted proteins")
        return prediction_file, secreted_ids

    def skip_signalp_get_all_proteins(self, input_fasta: str) -> List[str]:
        """
        Return all protein IDs without SignalP filtering (for Class I analysis).
        """
        logger.info("Getting all proteins for Class I analysis (no filtering)...")
        
        all_protein_ids = []
        with open(input_fasta, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                protein_id = record.id.split('|')[1] if '|' in record.id else record.id
                all_protein_ids.append(protein_id)
        
        logger.info(f"Found {len(all_protein_ids)} proteins for Class I analysis")
        return all_protein_ids

    def filter_secreted_proteins(self, proteome_file: str, secreted_ids: List[str], output_file: str) -> str:
        """
        Filter proteome to keep only secreted proteins.
        
        Args:
            proteome_file: Path to full proteome FASTA
            secreted_ids: List of secreted protein IDs
            output_file: Path to save filtered FASTA
            
        Returns:
            Path to filtered FASTA file
        """
        logger.info("Filtering for secreted proteins...")
        
        secreted_records = []
        with open(proteome_file, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                # Extract protein ID (handle different FASTA header formats)
                protein_id = record.id.split('|')[1] if '|' in record.id else record.id
                
                if protein_id in secreted_ids:
                    secreted_records.append(record)
        
        with open(output_file, 'w') as f:
            SeqIO.write(secreted_records, f, "fasta")
        
        logger.info(f"Filtered to {len(secreted_records)} secreted proteins")
        return output_file
    
    def generate_peptides(self, fasta_file: str, peptide_length: int = 15) -> List[Dict]:
        """
        Generate overlapping peptides from protein sequences.
        
        Args:
            fasta_file: Path to FASTA file
            peptide_length: Length of peptides to generate
            
        Returns:
            List of dictionaries with protein_id and peptide_sequence
        """
        logger.info(f"Generating {peptide_length}-mer peptides...")
        
        peptides = []
        with open(fasta_file, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                protein_id = record.id
                sequence = str(record.seq)
                
                # Generate overlapping peptides
                for i in range(len(sequence) - peptide_length + 1):
                    peptide = sequence[i:i + peptide_length]
                    # Skip peptides with ambiguous amino acids
                    if 'X' not in peptide and '*' not in peptide:
                        peptides.append({
                            'protein_id': protein_id,
                            'peptide_sequence': peptide,
                            'position': i + 1
                        })
        
        logger.info(f"Generated {len(peptides)} class II peptides")
        return peptides
    
    def generate_class_i_peptides(self, fasta_file: str, peptide_lengths: List[int] = [8, 9, 10, 11]) -> List[Dict]:
        """
        Generate overlapping peptides for HLA Class I (8-11mers).
        
        Args:
            fasta_file: Path to FASTA file
            peptide_lengths: List of peptide lengths to generate
            
        Returns:
            List of dictionaries with protein_id, peptide_sequence, length, position
        """
        logger.info(f"Generating Class I peptides (lengths: {peptide_lengths})...")
        
        peptides = []
        with open(fasta_file, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                protein_id = record.id
                sequence = str(record.seq)
                
                # Generate peptides of different lengths
                for length in peptide_lengths:
                    for i in range(len(sequence) - length + 1):
                        peptide = sequence[i:i + length]
                        # Skip peptides with ambiguous amino acids
                        if 'X' not in peptide and '*' not in peptide:
                            peptides.append({
                                'protein_id': protein_id,
                                'peptide_sequence': peptide,
                                'peptide_length': length,
                                'position': i + 1
                            })
        
        logger.info(f"Generated {len(peptides)} Class I peptides")
        return peptides

    def get_netmhciipan_alleles(self) -> List[str]:
        """Generate alleles in correct format for netMHCIIpan."""
        logger.info("Using predefined HLA class II allele list...")
        alleles = []
        
        # DR alleles
        dr_alleles = ['DRB59901', # all that we tested for selection with estimated_HLA_frequency>0.01
         # 'DRB11301',
         # 'DRB10701',
         # 'DRB10801',
         # 'DRB10301',
         # 'DRB10101',
         # 'DRB11601',
         # 'DRB50202',
         # 'DRB11401',
         # 'DRB11501',
         # 'DRB50101',
         # 'DRB10102',
         # 'DRB11201',
         # 'DRB11302',
         # 'DRB30301',
         # 'DRB11303',
         # 'DRB30101',
         # 'DRB39901',
         # 'DRB40101',
         # 'DRB40103',
         'DRB49901']
        
        # All combinations of alleles that we tested for selection with estimated_HLA_frequency>0.01
        dp_dq_combinations = ['HLA-DPA10103-DPB10101',
         # 'HLA-DPA10201-DPB10101',
         # 'HLA-DPA10202-DPB10101',
         # 'HLA-DQA10101-DPB10101',
         # 'HLA-DQA10102-DPB10101',
         # 'HLA-DQA10103-DPB10101',
         # 'HLA-DQA10201-DPB10101',
         # 'HLA-DQA10301-DPB10101',
         # 'HLA-DQA10401-DPB10101',
         # 'HLA-DQA10501-DPB10101',
         # 'HLA-DPA10103-DPB10401',
         # 'HLA-DPA10201-DPB10401',
         # 'HLA-DPA10202-DPB10401',
         # 'HLA-DQA10101-DPB10401',
         # 'HLA-DQA10102-DPB10401',
         # 'HLA-DQA10103-DPB10401',
         # 'HLA-DQA10201-DPB10401',
         # 'HLA-DQA10301-DPB10401',
         # 'HLA-DQA10401-DPB10401',
         # 'HLA-DQA10501-DPB10401',
         # 'HLA-DPA10103-DPB10402',
         # 'HLA-DPA10201-DPB10402',
         # 'HLA-DPA10202-DPB10402',
         # 'HLA-DQA10101-DPB10402',
         # 'HLA-DQA10102-DPB10402',
         # 'HLA-DQA10103-DPB10402',
         # 'HLA-DQA10201-DPB10402',
         # 'HLA-DQA10301-DPB10402',
         # 'HLA-DQA10401-DPB10402',
         # 'HLA-DQA10501-DPB10402',
         # 'HLA-DPA10103-DPB10501',
         # 'HLA-DPA10201-DPB10501',
         # 'HLA-DPA10202-DPB10501',
         # 'HLA-DQA10101-DPB10501',
         # 'HLA-DQA10102-DPB10501',
         # 'HLA-DQA10103-DPB10501',
         # 'HLA-DQA10201-DPB10501',
         # 'HLA-DQA10301-DPB10501',
         # 'HLA-DQA10401-DPB10501',
         # 'HLA-DQA10501-DPB10501',
         # 'HLA-DPA10103-DPB11001',
         # 'HLA-DPA10201-DPB11001',
         # 'HLA-DPA10202-DPB11001',
         # 'HLA-DQA10101-DPB11001',
         # 'HLA-DQA10102-DPB11001',
         # 'HLA-DQA10103-DPB11001',
         # 'HLA-DQA10201-DPB11001',
         # 'HLA-DQA10301-DPB11001',
         # 'HLA-DQA10401-DPB11001',
         # 'HLA-DQA10501-DPB11001',
         # 'HLA-DPA10103-DPB11101',
         # 'HLA-DPA10201-DPB11101',
         # 'HLA-DPA10202-DPB11101',
         # 'HLA-DQA10101-DPB11101',
         # 'HLA-DQA10102-DPB11101',
         # 'HLA-DQA10103-DPB11101',
         # 'HLA-DQA10201-DPB11101',
         # 'HLA-DQA10301-DPB11101',
         # 'HLA-DQA10401-DPB11101',
         # 'HLA-DQA10501-DPB11101',
         # 'HLA-DPA10103-DPB11301',
         # 'HLA-DPA10201-DPB11301',
         # 'HLA-DPA10202-DPB11301',
         # 'HLA-DQA10101-DPB11301',
         # 'HLA-DQA10102-DPB11301',
         # 'HLA-DQA10103-DPB11301',
         # 'HLA-DQA10201-DPB11301',
         # 'HLA-DQA10301-DPB11301',
         # 'HLA-DQA10401-DPB11301',
         # 'HLA-DQA10501-DPB11301',
         # 'HLA-DPA10103-DPB11701',
         # 'HLA-DPA10201-DPB11701',
         # 'HLA-DPA10202-DPB11701',
         # 'HLA-DQA10101-DPB11701',
         # 'HLA-DQA10102-DPB11701',
         # 'HLA-DQA10103-DPB11701',
         # 'HLA-DQA10201-DPB11701',
         # 'HLA-DQA10301-DPB11701',
         # 'HLA-DQA10401-DPB11701',
         # 'HLA-DQA10501-DPB11701',
         # 'HLA-DPA10103-DQB10501',
         # 'HLA-DPA10201-DQB10501',
         # 'HLA-DPA10202-DQB10501',
         # 'HLA-DQA10101-DQB10501',
         # 'HLA-DQA10102-DQB10501',
         # 'HLA-DQA10103-DQB10501',
         # 'HLA-DQA10201-DQB10501',
         # 'HLA-DQA10301-DQB10501',
         # 'HLA-DQA10401-DQB10501',
         # 'HLA-DQA10501-DQB10501',
         # 'HLA-DPA10103-DQB10603',
         # 'HLA-DPA10201-DQB10603',
         # 'HLA-DPA10202-DQB10603',
         # 'HLA-DQA10101-DQB10603',
         # 'HLA-DQA10102-DQB10603',
         # 'HLA-DQA10103-DQB10603',
         # 'HLA-DQA10201-DQB10603',
         # 'HLA-DQA10301-DQB10603',
         # 'HLA-DQA10401-DQB10603',
         # 'HLA-DQA10501-DQB10603',
         # 'HLA-DPA10103-DQB10202',
         # 'HLA-DPA10201-DQB10202',
         # 'HLA-DPA10202-DQB10202',
         # 'HLA-DQA10101-DQB10202',
         # 'HLA-DQA10102-DQB10202',
         # 'HLA-DQA10103-DQB10202',
         # 'HLA-DQA10201-DQB10202',
         # 'HLA-DQA10301-DQB10202',
         # 'HLA-DQA10401-DQB10202',
         # 'HLA-DQA10501-DQB10202',
         # 'HLA-DPA10103-DQB10402',
         # 'HLA-DPA10201-DQB10402',
         # 'HLA-DPA10202-DQB10402',
         # 'HLA-DQA10101-DQB10402',
         # 'HLA-DQA10102-DQB10402',
         # 'HLA-DQA10103-DQB10402',
         # 'HLA-DQA10201-DQB10402',
         # 'HLA-DQA10301-DQB10402',
         # 'HLA-DQA10401-DQB10402',
         # 'HLA-DQA10501-DQB10402',
         # 'HLA-DPA10103-DQB10201',
         # 'HLA-DPA10201-DQB10201',
         # 'HLA-DPA10202-DQB10201',
         # 'HLA-DQA10101-DQB10201',
         # 'HLA-DQA10102-DQB10201',
         # 'HLA-DQA10103-DQB10201',
         # 'HLA-DQA10201-DQB10201',
         # 'HLA-DQA10301-DQB10201',
         # 'HLA-DQA10401-DQB10201',
         # 'HLA-DQA10501-DQB10201',
         # 'HLA-DPA10103-DQB10301',
         # 'HLA-DPA10201-DQB10301',
         # 'HLA-DPA10202-DQB10301',
         # 'HLA-DQA10101-DQB10301',
         # 'HLA-DQA10102-DQB10301',
         # 'HLA-DQA10103-DQB10301',
         # 'HLA-DQA10201-DQB10301',
         # 'HLA-DQA10301-DQB10301',
         # 'HLA-DQA10401-DQB10301',
         # 'HLA-DQA10501-DQB10301',
         # 'HLA-DPA10103-DQB10302',
         # 'HLA-DPA10201-DQB10302',
         # 'HLA-DPA10202-DQB10302',
         # 'HLA-DQA10101-DQB10302',
         # 'HLA-DQA10102-DQB10302',
         # 'HLA-DQA10103-DQB10302',
         # 'HLA-DQA10201-DQB10302',
         # 'HLA-DQA10301-DQB10302',
         # 'HLA-DQA10401-DQB10302',
         # 'HLA-DQA10501-DQB10302',
         # 'HLA-DPA10103-DQB10303',
         # 'HLA-DPA10201-DQB10303',
         # 'HLA-DPA10202-DQB10303',
         # 'HLA-DQA10101-DQB10303',
         # 'HLA-DQA10102-DQB10303',
         # 'HLA-DQA10103-DQB10303',
         # 'HLA-DQA10201-DQB10303',
         # 'HLA-DQA10301-DQB10303',
         # 'HLA-DQA10401-DQB10303',
         # 'HLA-DQA10501-DQB10303',
         # 'HLA-DPA10103-DQB10502',
         # 'HLA-DPA10201-DQB10502',
         # 'HLA-DPA10202-DQB10502',
         # 'HLA-DQA10101-DQB10502',
         # 'HLA-DQA10102-DQB10502',
         # 'HLA-DQA10103-DQB10502',
         # 'HLA-DQA10201-DQB10502',
         # 'HLA-DQA10301-DQB10502',
         # 'HLA-DQA10401-DQB10502',
         # 'HLA-DQA10501-DQB10502',
         # 'HLA-DPA10103-DQB10503',
         # 'HLA-DPA10201-DQB10503',
         # 'HLA-DPA10202-DQB10503',
         # 'HLA-DQA10101-DQB10503',
         # 'HLA-DQA10102-DQB10503',
         # 'HLA-DQA10103-DQB10503',
         # 'HLA-DQA10201-DQB10503',
         # 'HLA-DQA10301-DQB10503',
         # 'HLA-DQA10401-DQB10503',
         # 'HLA-DQA10501-DQB10503',
         # 'HLA-DPA10103-DQB10602',
         # 'HLA-DPA10201-DQB10602',
         # 'HLA-DPA10202-DQB10602',
         # 'HLA-DQA10101-DQB10602',
         # 'HLA-DQA10102-DQB10602',
         # 'HLA-DQA10103-DQB10602',
         # 'HLA-DQA10201-DQB10602',
         # 'HLA-DQA10301-DQB10602',
         # 'HLA-DQA10401-DQB10602',
         # 'HLA-DQA10501-DQB10602',
         # 'HLA-DPA10103-DQB10604',
         # 'HLA-DPA10201-DQB10604',
         # 'HLA-DPA10202-DQB10604',
         # 'HLA-DQA10101-DQB10604',
         # 'HLA-DQA10102-DQB10604',
         # 'HLA-DQA10103-DQB10604',
         # 'HLA-DQA10201-DQB10604',
         # 'HLA-DQA10301-DQB10604',
         # 'HLA-DQA10401-DQB10604',
         # 'HLA-DQA10501-DQB10604',
         # 'HLA-DPA10103-DQB10609',
         # 'HLA-DPA10201-DQB10609',
         # 'HLA-DPA10202-DQB10609',
         # 'HLA-DQA10101-DQB10609',
         # 'HLA-DQA10102-DQB10609',
         # 'HLA-DQA10103-DQB10609',
         # 'HLA-DQA10201-DQB10609',
         # 'HLA-DQA10301-DQB10609',
         # 'HLA-DQA10401-DQB10609',
         'HLA-DQA10501-DQB10609']

        alleles = dr_alleles + dp_dq_combinations

        logger.info(f"Using {len(alleles)} HLA class II alleles")
        return alleles
    
    def get_netmhcpan_alleles(self) -> List[str]:
        """
        Get list of HLA Class I alleles for netMHCpan.
        """
        logger.info("Using predefined HLA Class I allele list...")
        
        # Common HLA Class I alleles
        alleles = [ # all that we tested for selection with estimated_HLA_frequency>0.01
         'HLA-A01:01',
         # 'HLA-A02:01',
         # 'HLA-A03:01',
         # 'HLA-A11:01',
         # 'HLA-A23:01',
         # 'HLA-A24:02',
         # 'HLA-A25:01',
         # 'HLA-A29:02',
         # 'HLA-A30:01',
         # 'HLA-A31:01',
         # 'HLA-A32:01',
         # 'HLA-A68:01',
         # 'HLA-B07:02',
         # 'HLA-C07:02',
         # 'HLA-B08:01',
         # 'HLA-B13:02',
         # 'HLA-B14:01',
         # 'HLA-B14:02',
         # 'HLA-B15:01',
         # 'HLA-B18:01',
         # 'HLA-B27:05',
         # 'HLA-B37:01',
         # 'HLA-B38:01',
         # 'HLA-B40:01',
         # 'HLA-C03:04',
         # 'HLA-B44:02',
         # 'HLA-B44:03',
         # 'HLA-B49:01',
         # 'HLA-B50:01',
         # 'HLA-B51:01',
         # 'HLA-B55:01',
         # 'HLA-B57:01',
         # 'HLA-C01:02',
         # 'HLA-C02:02',
         # 'HLA-C03:03',
         # 'HLA-C04:01',
         # 'HLA-C05:01',
         # 'HLA-C06:02',
         # 'HLA-C07:01',
         # 'HLA-C07:04',
         # 'HLA-C08:02',
         # 'HLA-C12:03',
         # 'HLA-C14:02',
         # 'HLA-C15:02',
         'HLA-C16:01']
        
        logger.info(f"Using {len(alleles)} HLA Class I alleles")
        return alleles

    def run_netmhciipan_batch(self, peptides: List[Dict], alleles: List[str], output_dir: str, n_jobs: int = 30) -> str:
        """
        Run netMHCIIpan predictions in parallel batches.
        """
        logger.info(f"Running netMHCIIpan predictions for {len(peptides)} peptides against {len(alleles)} alleles...")
        logger.info(f"Using {n_jobs} parallel jobs")
        
        # Create peptide file
        peptide_file = os.path.join(output_dir, "peptides.txt")
        with open(peptide_file, 'w') as f:
            for peptide in peptides:
                f.write(f"{peptide['peptide_sequence']}\n")
        
        # Create commands file for parallel execution
        commands_file = os.path.join(output_dir, "netmhciipan_commands.txt")
        temp_dir = os.path.join(output_dir, "temp_results")
        os.makedirs(temp_dir, exist_ok=True)
        
        with open(commands_file, 'w') as f:
            for i, allele in enumerate(alleles):
                output_file = os.path.join(temp_dir, f"result_{allele}_{i:03d}.txt")
                cmd = f"{self.netmhciipan_path} -f {peptide_file} -a {allele} -inptype 1 -xls -xlsfile {output_file}"
                f.write(f"{cmd}\n")
        
        # Run in parallel using GNU parallel
        logger.info("Running parallel netMHCIIpan jobs...")
        parallel_cmd = f"cat {commands_file} | parallel -j {n_jobs}"
        
        try:
            result = subprocess.run(parallel_cmd, shell=True, capture_output=True, text=True, check=True)
            logger.info("Parallel netMHCIIpan completed successfully")
        except subprocess.CalledProcessError as e:
            logger.error(f"Parallel netMHCIIpan failed: {e}")
            raise
        
        # Combine all result files
        combined_output = os.path.join(output_dir, "netmhciipan_results.txt")
        self._combine_netmhciipan_results(temp_dir, combined_output)
        
        return combined_output

    def _combine_netmhciipan_results(self, temp_dir: str, output_file: str):
        """
        Combine individual netMHCIIpan result files into one.
        """
        logger.info("Combining parallel results...")
        
        with open(output_file, 'w') as outf:
            for result_file in sorted(os.listdir(temp_dir)):
                if not result_file.endswith('.txt'):
                    continue
                    
                file_path = os.path.join(temp_dir, result_file)
                try:
                    with open(file_path, 'r') as inf:
                        # Copy all content from each file
                        outf.write(inf.read())
                        outf.write('\n')  # Add separator between files
                        
                except Exception as e:
                    logger.warning(f"Failed to read {result_file}: {e}")
        
        logger.info(f"Combined {len(os.listdir(temp_dir))} result files")
    
    def run_netmhcpan_batch(self, peptides: List[Dict], alleles: List[str], output_dir: str, n_jobs: int = 30) -> str:
        """
        Run netMHCpan (Class I) predictions in parallel batches.
        """
        logger.info(f"Running netMHCpan (Class I) predictions for {len(peptides)} peptides against {len(alleles)} alleles...")
        logger.info(f"Using {n_jobs} parallel jobs")
        
        # Create peptide file
        peptide_file = os.path.join(output_dir, "class_i_peptides.txt")
        with open(peptide_file, 'w') as f:
            for peptide in peptides:
                f.write(f"{peptide['peptide_sequence']}\n")
        
        # Create commands file for parallel execution
        commands_file = os.path.join(output_dir, "netmhcpan_commands.txt")
        temp_dir = os.path.join(output_dir, "temp_results")
        os.makedirs(temp_dir, exist_ok=True)
        
        with open(commands_file, 'w') as f:
            for i, allele in enumerate(alleles):
                output_file = os.path.join(temp_dir, f"class_i_result_{allele}_{i:03d}.txt")
                cmd = f"{self.netmhcpan_path} -f {peptide_file} -a {allele} -inptype 1 -xls -xlsfile {output_file}"
                f.write(f"{cmd}\n")
        
        # Run in parallel
        logger.info("Running parallel netMHCpan jobs...")
        parallel_cmd = f"cat {commands_file} | parallel -j {n_jobs}"
        
        try:
            result = subprocess.run(parallel_cmd, shell=True, capture_output=True, text=True, check=True)
            logger.info("Parallel netMHCpan completed successfully")
        except subprocess.CalledProcessError as e:
            logger.error(f"Parallel netMHCpan failed: {e}")
            raise
        
        # Combine results
        combined_output = os.path.join(output_dir, "netmhcpan_results.txt")
        self._combine_netmhciipan_results(temp_dir, combined_output)
        
        return combined_output
    
    def parse_netmhciipan_results(self, results_file: str, peptides: List[Dict]) -> pd.DataFrame:
        """
        Parse netMHCIIpan results and combine with peptide metadata.
        """
        logger.info("Parsing netMHCIIpan results...")
        
        # Create peptide lookup
        peptide_lookup = {p['peptide_sequence']: p for p in peptides}
        
        all_results = []
        
        # Read the combined results file
        try:
            with open(results_file, 'r') as f:
                current_allele = None
                
                for line in f:
                    line = line.strip()
                    
                    # Skip empty lines
                    if not line:
                        continue
                    
                    # Check if this line contains an allele name (header line)
                    # Must start with HLA- or be just the allele name, and not contain tabs
                    if ('\t' not in line and 
                        (line.startswith('HLA-') or 
                         line.startswith('DRB') or 
                         line.startswith('DQA') or 
                         line.startswith('DPA')) and
                        'Pos' not in line):
                        current_allele = line.strip()
                        continue
                    
                    # Skip header lines
                    if line.startswith('Pos') or 'Peptide ID' in line:
                        continue
                    
                    # Parse data lines
                    parts = line.split('\t')
                    if len(parts) >= 7 and parts[0].isdigit():
                        try:
                            position = int(parts[0])
                            peptide_sequence = parts[1]  # This is actually the peptide sequence
                            score = float(parts[6]) if parts[6] != 'NA' else None
                            rank = float(parts[7]) if parts[7] != 'NA' else None
                            
                            # Get protein metadata
                            peptide_info = peptide_lookup.get(peptide_sequence, {})
                            protein_id = peptide_info.get('protein_id', 'Unknown')
                            protein_position = peptide_info.get('position', 0)
                            
                            all_results.append({
                                'protein_id': protein_id,
                                'peptide_sequence': peptide_sequence,
                                'position': protein_position,
                                'hla_allele': current_allele,
                                'Score_EL': score,  # Note: this might be a different score type
                                '%Rank_EL': rank
                            })
                            
                        except ValueError as e:
                            logger.warning(f"Could not parse line: {line} - {e}")
                            continue
                            
        except Exception as e:
            logger.error(f"Failed to parse netMHCIIpan results: {e}")
            raise
        
        if not all_results:
            raise ValueError("No valid results found in netMHCIIpan output")
        
        df = pd.DataFrame(all_results)
        logger.info(f"Parsed {len(df)} prediction results")
        return df
    
    def parse_netmhcpan_results(self, results_file: str, peptides: List[Dict]) -> pd.DataFrame:
        """
        Parse netMHCpan (Class I) results.
        """
        logger.info("Parsing netMHCpan (Class I) results...")
        
        peptide_lookup = {p['peptide_sequence']: p for p in peptides}
        all_results = []
        
        try:
            with open(results_file, 'r') as f:
                current_allele = None
                
                for line in f:
                    line = line.strip()
                    
                    # Skip empty lines
                    if not line:
                        continue
                    
                    # Check if this line contains an allele name (header line)
                    # Must start with HLA- and not contain tabs
                    if ('\t' not in line and 
                        line.startswith('HLA-') and
                        'Pos' not in line):
                        current_allele = line.strip()
                        continue
                    
                    # Skip header lines
                    if line.startswith('Pos') or 'Peptide ID' in line:
                        continue
                    
                    # Parse data lines
                    parts = line.split('\t')
                    if len(parts) >= 7 and parts[0].isdigit():
                        try:
                            position = int(parts[0])
                            peptide_sequence = parts[1]  # Peptide sequence is in column 1
                            score = float(parts[5]) if parts[5] != 'NA' else None  # Score in column 4
                            rank = float(parts[6]) if parts[6] != 'NA' else None   # Rank in column 5
                            
                            # Skip if peptide_seq is invalid
                            if pd.isna(peptide_sequence) or not isinstance(peptide_sequence, str):
                                continue
                                
                            # Get protein metadata
                            peptide_info = peptide_lookup.get(peptide_sequence, {})
                            protein_id = peptide_info.get('protein_id', 'Unknown')
                            protein_position = peptide_info.get('position', 0)
                            peptide_length = peptide_info.get('peptide_length', len(peptide_sequence))
                            
                            all_results.append({
                                'protein_id': protein_id,
                                'peptide_sequence': peptide_sequence,
                                'peptide_length': peptide_length,
                                'position': protein_position,
                                'hla_allele': current_allele,
                                'Score_EL': score,
                                '%Rank_EL': rank
                            })
                            
                        except ValueError as e:
                            logger.warning(f"Could not parse line: {line} - {e}")
                            continue
                            
        except Exception as e:
            logger.error(f"Failed to parse netMHCpan results: {e}")
            raise
        
        if not all_results:
            raise ValueError("No valid results found in netMHCpan output")
        
        df = pd.DataFrame(all_results)
        logger.info(f"Parsed {len(df)} Class I prediction results")
        return df
    
    # def analyze_organism(self, organism_query: str, output_dir: str, organism_name: str = None, n_jobs: int = 20) -> str:
    #     """
    #     Complete analysis pipeline for both Class I and Class II.
        
    #     Args:
    #         organism_query: NCBI search query
    #         output_dir: Output directory
    #         organism_name: Name for output files (optional)
    #         n_jobs: Number of parallel jobs
            
    #     Returns:
    #         Path to final Excel results file
    #     """
    #     if organism_name is None:
    #         organism_name = re.sub(r'[^\w\-_\.]', '_', organism_query)
        
    #     os.makedirs(output_dir, exist_ok=True)
        
    #     # Steps 1-3: Same as before (download, SignalP, filter)
    #     logger.info("=== Starting Pathogen HLA Analysis Pipeline ===")
        
    #     # Step 1: Download proteome
    #     proteome_file = os.path.join(output_dir, f"{organism_name}_proteome.fasta")
    #     self.download_proteome(organism_query, proteome_file)
        
    #     # # Step 2: Run SignalP
    #     # signalp_dir = os.path.join(output_dir, "signalp")
    #     # os.makedirs(signalp_dir, exist_ok=True)
    #     # _, secreted_ids = self.run_signalp(proteome_file, signalp_dir)
        
    #     # # Step 3: Filter for secreted proteins
    #     # secreted_file = os.path.join(output_dir, f"{organism_name}_secreted.fasta")
    #     # self.filter_secreted_proteins(proteome_file, secreted_ids, secreted_file)

    #     # Step 2: Skip filtering - use all proteins
    #     logger.info("Using all proteins for analysis (no SignalP filtering)")
    #     _, all_protein_ids = self.run_signalp(proteome_file, output_dir)

    #     # Step 3: Copy proteome file (no filtering needed)
    #     secreted_file = os.path.join(output_dir, f"{organism_name}_all_proteins.fasta")
    #     import shutil
    #     shutil.copy2(proteome_file, secreted_file)
        
    #     # Step 4a: Run Class II analysis
    #     logger.info("=== Running HLA Class II Analysis ===")
    #     class_ii_peptides = self.generate_peptides(secreted_file, 15)  # 15mers for Class II
        
    #     if not class_ii_peptides:
    #         raise ValueError("No Class II peptides generated. Check secreted protein filtering.")
        
    #     class_ii_alleles = self.get_netmhciipan_alleles()
        
    #     netmhciipan_dir = os.path.join(output_dir, "netmhciipan_class_ii")
    #     os.makedirs(netmhciipan_dir, exist_ok=True)
    #     class_ii_results_file = self.run_netmhciipan_batch(class_ii_peptides, class_ii_alleles, netmhciipan_dir, n_jobs)
    #     class_ii_df = self.parse_netmhciipan_results(class_ii_results_file, class_ii_peptides)
    #     class_ii_df['hla_class'] = 'II'
        
    #     # Step 4b: Run Class I analysis  
    #     logger.info("=== Running HLA Class I Analysis ===")
    #     class_i_peptides = self.generate_class_i_peptides(secreted_file, [8, 9, 10, 11])
        
    #     if not class_i_peptides:
    #         raise ValueError("No Class I peptides generated. Check secreted protein filtering.")
        
    #     class_i_alleles = self.get_netmhcpan_alleles()
        
    #     netmhcpan_dir = os.path.join(output_dir, "netmhcpan_class_i")
    #     os.makedirs(netmhcpan_dir, exist_ok=True)
    #     class_i_results_file = self.run_netmhcpan_batch(class_i_peptides, class_i_alleles, netmhcpan_dir, n_jobs)
    #     class_i_df = self.parse_netmhcpan_results(class_i_results_file, class_i_peptides)
    #     class_i_df['hla_class'] = 'I'
        
    #     # Step 5: Combine results
    #     logger.info("=== Combining Results ===")
        
    #     # Ensure both dataframes have consistent columns
    #     # Add peptide_length column for Class II if missing
    #     if 'peptide_length' not in class_ii_df.columns:
    #         class_ii_df['peptide_length'] = class_ii_df['peptide_sequence'].str.len()
        
    #     # Combine the dataframes
    #     combined_df = pd.concat([class_ii_df, class_i_df], ignore_index=True)
        
    #     # Reorder columns for consistency
    #     columns = ['protein_id', 'peptide_sequence', 'peptide_length', 'position', 
    #                'hla_class', 'hla_allele', 'Score_EL', '%Rank_EL']
        
    #     # Only keep columns that exist
    #     available_columns = [col for col in columns if col in combined_df.columns]
    #     combined_df = combined_df[available_columns]
        
    #     # Step 6: Save combined results
    #     excel_file = os.path.join(output_dir, f"{organism_name}_hla_binding_results.csv")
    #     combined_df.to_csv(excel_file, index=False)
        
    #     logger.info("=== Analysis Complete ===")
    #     logger.info(f"Class I results: {len(class_i_df)} predictions")
    #     logger.info(f"Class II results: {len(class_ii_df)} predictions") 
    #     logger.info(f"Total results: {len(combined_df)} predictions")
    #     logger.info(f"Combined results saved to: {excel_file}")
        
    #     return excel_file

    def print_analysis_instructions():
        """
        Print instructions for interpreting Class I and Class II results.
        """
        print("\n" + "="*80)
        print("IMPORTANT: HOW TO INTERPRET CLASS I vs CLASS II RESULTS")
        print("="*80)
        
        print("\n🔬 ANALYSIS APPROACH USED:")
        print("• Class I: ALL proteins analyzed (comprehensive intracellular targets)")
        print("• Class II: SignalP-filtered proteins only (extracellular targets)")
        print("• This reflects biological antigen processing pathways!")
        
        print("\n🔬 BIOLOGICAL DIFFERENCES:")
        print("• Class I (CD8+ T cells): Cytotoxic response, kills infected cells")
        print("• Class II (CD4+ T cells): Helper response, coordinates immune response")
        print("• Both are needed for comprehensive immunity!")
        
        print("\n📊 ANALYSIS RECOMMENDATIONS:")
        print("• Use %Rank_EL as primary metric (more reliable than IC50)")
        print("• Analyze each class SEPARATELY - don't directly compare scores")
        
        print("\n💡 SUGGESTED WORKFLOW:")
        print("1. Filter strong binders by class:")
        print("   class_i_strong = df[(df['hla_class']=='I') & (df['%Rank_EL']<0.5)]")
        print("   class_ii_strong = df[(df['hla_class']=='II') & (df['%Rank_EL']<1)]")
        
        print("\n2. Identify proteins with broad HLA coverage:")
        print("   coverage_i = class_i_strong.groupby('protein_id')['hla_allele'].nunique()")
        print("   coverage_ii = class_ii_strong.groupby('protein_id')['hla_allele'].nunique()")
        
        print("\n3. Find proteins with BOTH Class I and II epitopes (best vaccine targets):")
        print("   dual_targets = set(coverage_i.index) & set(coverage_ii.index)")
        
        print("\n⚠️  IMPORTANT NOTES:")
        print("• Class I dataset is larger (all proteins) - expect more predictions")
        print("• Class II dataset is smaller (secreted only) - more focused results")
        print("• This size difference is EXPECTED and biologically appropriate")
        
        print("="*80 + "\n")

    def analyze_organism(self, organism_query: str, output_dir: str, organism_name: str = None, n_jobs: int = 20) -> str:
        """
        Complete analysis pipeline for both Class I (all proteins) and Class II (SignalP filtered).
        """
        if organism_name is None:
            organism_name = re.sub(r'[^\w\-_\.]', '_', organism_query)
        
        os.makedirs(output_dir, exist_ok=True)
        
        logger.info("=== Starting Pathogen HLA Analysis Pipeline ===")
        
        # Step 1: Download proteome
        proteome_file = os.path.join(output_dir, f"{organism_name}_proteome.fasta")
        self.download_proteome(organism_query, proteome_file)
        
        # Step 2a: Run Class I analysis (ALL PROTEINS)
        logger.info("=== Running HLA Class I Analysis (All Proteins) ===")
        all_protein_ids = self.skip_signalp_get_all_proteins(proteome_file)
        
        # Create file with all proteins for Class I
        class_i_file = os.path.join(output_dir, f"{organism_name}_all_proteins.fasta")
        import shutil
        shutil.copy2(proteome_file, class_i_file)
        
        # Generate Class I peptides from all proteins
        class_i_peptides = self.generate_class_i_peptides(class_i_file, [8, 9, 10, 11])
        
        if not class_i_peptides:
            raise ValueError("No Class I peptides generated.")
        
        class_i_alleles = self.get_netmhcpan_alleles()
        
        netmhcpan_dir = os.path.join(output_dir, "netmhcpan_class_i")
        os.makedirs(netmhcpan_dir, exist_ok=True)
        class_i_results_file = self.run_netmhcpan_batch(class_i_peptides, class_i_alleles, netmhcpan_dir, n_jobs)
        class_i_df = self.parse_netmhcpan_results(class_i_results_file, class_i_peptides)
        class_i_df['hla_class'] = 'I'
        
        # Step 2b: Run Class II analysis (SIGNALP FILTERED)
        logger.info("=== Running HLA Class II Analysis (SignalP Filtered) ===")
        
        # Run SignalP for Class II
        signalp_dir = os.path.join(output_dir, "signalp")
        os.makedirs(signalp_dir, exist_ok=True)
        _, secreted_ids = self.run_signalp(proteome_file, signalp_dir)
        
        # Filter for secreted proteins for Class II
        secreted_file = os.path.join(output_dir, f"{organism_name}_secreted.fasta")
        self.filter_secreted_proteins(proteome_file, secreted_ids, secreted_file)
        
        # Generate Class II peptides from secreted proteins only
        class_ii_peptides = self.generate_peptides(secreted_file, 15)  # 15mers for Class II
        
        if not class_ii_peptides:
            raise ValueError("No Class II peptides generated. Check secreted protein filtering.")
        
        class_ii_alleles = self.get_netmhciipan_alleles()
        
        netmhciipan_dir = os.path.join(output_dir, "netmhciipan_class_ii")
        os.makedirs(netmhciipan_dir, exist_ok=True)
        class_ii_results_file = self.run_netmhciipan_batch(class_ii_peptides, class_ii_alleles, netmhciipan_dir, n_jobs)
        class_ii_df = self.parse_netmhciipan_results(class_ii_results_file, class_ii_peptides)
        class_ii_df['hla_class'] = 'II'
        
        # Step 3: Combine results
        logger.info("=== Combining Results ===")
        
        # Ensure both dataframes have consistent columns
        if 'peptide_length' not in class_ii_df.columns:
            class_ii_df['peptide_length'] = class_ii_df['peptide_sequence'].str.len()
        
        # Combine the dataframes
        combined_df = pd.concat([class_i_df, class_ii_df], ignore_index=True)
        
        # Reorder columns for consistency
        columns = ['protein_id', 'peptide_sequence', 'peptide_length', 'position', 
                   'hla_class', 'hla_allele', 'Score_EL', '%Rank_EL']
        
        available_columns = [col for col in columns if col in combined_df.columns]
        combined_df = combined_df[available_columns]
        
        # Step 4: Save combined results
        excel_file = os.path.join(output_dir, f"{organism_name}_hla_binding_results.csv")
        combined_df.to_csv(excel_file, index=False)
        
        logger.info("=== Analysis Complete ===")
        logger.info(f"Class I results (all proteins): {len(class_i_df)} predictions")
        logger.info(f"Class II results (secreted only): {len(class_ii_df)} predictions") 
        logger.info(f"Total results: {len(combined_df)} predictions")
        logger.info(f"Combined results saved to: {excel_file}")
        
        # Print analysis instructions
        self.print_analysis_instructions()
        
        return excel_file

def main():
    parser = argparse.ArgumentParser(description="Pathogen HLA Class I & II Binding Analysis")
    parser.add_argument("--email", required=True, help="Email for NCBI Entrez API")
    parser.add_argument("--organism", default="Variola virus[Organism]", 
                       help="NCBI organism query (default: Variola virus[Organism])")
    parser.add_argument("--output-dir", default="./hla_analysis", 
                       help="Output directory (default: ./hla_analysis)")
    parser.add_argument("--name", help="Name for output files (auto-generated if not provided)")
    parser.add_argument("--netmhciipan-path", default="/home/jdn321/software/netMHCIIpan-4.3/netMHCIIpan",
                       help="Path to netMHCIIpan executable (Class II)")
    parser.add_argument("--netmhcpan-path", default="/home/jdn321/software/netMHCpan-4.2/netMHCpan",
                       help="Path to netMHCpan executable (Class I)")
    parser.add_argument("--n-jobs", type=int, default=20,
                       help="Number of parallel jobs (default: 20)")
    
    args = parser.parse_args()
    
    try:
        analyzer = PathogenHLAAnalyzer(args.email, args.netmhciipan_path, args.netmhcpan_path)
        result_file = analyzer.analyze_organism(args.organism, args.output_dir, args.name, args.n_jobs)
        print(f"Analysis completed successfully!")
        print(f"Results saved to: {result_file}")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
