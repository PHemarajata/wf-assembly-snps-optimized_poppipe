#!/usr/bin/env python3
"""
Validate PopPUNK cluster file against FASTA directory to check for sample name matching issues.
"""

import pandas as pd
import argparse
import glob
import os
from pathlib import Path

def create_sample_variants(sample_name):
    """Create multiple variants of a sample name for flexible matching."""
    variants = set()
    
    # 1. Original name
    variants.add(sample_name)
    
    # 2. Remove common suffixes
    clean_name = sample_name.replace('-SPAdes', '').replace('_contigs', '').replace('_genomic', '')
    variants.add(clean_name)
    
    # 3. Extract base sample ID (before underscore)
    if '_' in clean_name:
        base_name = clean_name.split('_')[0]
        if len(base_name) > 3:
            variants.add(base_name)
    
    # 4. Extract base sample ID (before dash)
    if '-' in clean_name:
        base_name = clean_name.split('-')[0]
        if len(base_name) > 3:
            variants.add(base_name)
    
    # 5. For GCA/GCF accessions, handle versions
    if clean_name.startswith(('GCA_', 'GCF_')):
        parts = clean_name.split('_')
        if len(parts) >= 2:
            accession = '_'.join(parts[:2])  # GCA_002587385
            variants.add(accession)
            if '.' in parts[1]:
                number_part = parts[1].split('.')[0]
                variants.add(f"{parts[0]}_{number_part}")
    
    return variants

def main():
    parser = argparse.ArgumentParser(description='Validate PopPUNK cluster file against FASTA directory')
    parser.add_argument('--poppunk_csv', required=True, help='PopPUNK cluster assignment CSV file')
    parser.add_argument('--fasta_dir', required=True, help='Directory containing FASTA files')
    parser.add_argument('--output', default='validation_report.txt', help='Output validation report')
    
    args = parser.parse_args()
    
    # Read PopPUNK results
    try:
        df = pd.read_csv(args.poppunk_csv)
        print(f"✅ Loaded PopPUNK results: {len(df)} samples")
    except Exception as e:
        print(f"❌ Error reading PopPUNK CSV: {e}")
        return 1
    
    # Find FASTA files
    fasta_extensions = ['.fasta', '.fas', '.fna', '.fsa', '.fa', '.fasta.gz', '.fas.gz', '.fna.gz', '.fsa.gz', '.fa.gz']
    fasta_files = []
    
    for ext in fasta_extensions:
        pattern = os.path.join(args.fasta_dir, f'*{ext}')
        fasta_files.extend(glob.glob(pattern))
    
    print(f"✅ Found {len(fasta_files)} FASTA files")
    
    # Create sample name to file mapping
    sample_to_file = {}
    file_variants = {}
    
    for file_path in fasta_files:
        file_name = os.path.basename(file_path)
        # Extract sample name from filename (remove extensions)
        sample_name = Path(file_name).stem
        
        # Remove common genomic file extensions
        for ext in ['.genomic', '.fasta', '.fas', '.fna', '.fsa', '.fa']:
            if sample_name.endswith(ext):
                sample_name = sample_name[:-len(ext)]
        
        # Create variants
        variants = create_sample_variants(sample_name)
        
        # Store all variants
        for variant in variants:
            if variant:  # Only non-empty variants
                sample_to_file[variant] = file_path
                if variant not in file_variants:
                    file_variants[variant] = []
                file_variants[variant].append(file_path)
    
    print(f"✅ Created {len(sample_to_file)} sample name variants")
    
    # Check matching
    matched_samples = []
    unmatched_samples = []
    
    for _, row in df.iterrows():
        sample_name = row['Taxon']
        if sample_name in sample_to_file:
            matched_samples.append((sample_name, sample_to_file[sample_name]))
        else:
            unmatched_samples.append(sample_name)
    
    # Generate report
    with open(args.output, 'w') as f:
        f.write("PopPUNK Cluster Validation Report\n")
        f.write("=" * 50 + "\n\n")
        
        f.write(f"PopPUNK file: {args.poppunk_csv}\n")
        f.write(f"FASTA directory: {args.fasta_dir}\n")
        f.write(f"Total PopPUNK samples: {len(df)}\n")
        f.write(f"Total FASTA files: {len(fasta_files)}\n")
        f.write(f"Matched samples: {len(matched_samples)}\n")
        f.write(f"Unmatched samples: {len(unmatched_samples)}\n")
        f.write(f"Match rate: {len(matched_samples)/len(df)*100:.1f}%\n\n")
        
        if matched_samples:
            f.write("MATCHED SAMPLES:\n")
            f.write("-" * 30 + "\n")
            for sample, file_path in matched_samples[:20]:  # Show first 20
                f.write(f"{sample} -> {os.path.basename(file_path)}\n")
            if len(matched_samples) > 20:
                f.write(f"... and {len(matched_samples) - 20} more\n")
            f.write("\n")
        
        if unmatched_samples:
            f.write("UNMATCHED SAMPLES:\n")
            f.write("-" * 30 + "\n")
            for sample in unmatched_samples[:20]:  # Show first 20
                f.write(f"{sample}\n")
            if len(unmatched_samples) > 20:
                f.write(f"... and {len(unmatched_samples) - 20} more\n")
            f.write("\n")
        
        f.write("AVAILABLE FASTA FILES:\n")
        f.write("-" * 30 + "\n")
        for file_path in fasta_files[:20]:  # Show first 20
            f.write(f"{os.path.basename(file_path)}\n")
        if len(fasta_files) > 20:
            f.write(f"... and {len(fasta_files) - 20} more\n")
    
    # Print summary
    print(f"\n📊 VALIDATION SUMMARY:")
    print(f"   Total samples: {len(df)}")
    print(f"   Matched: {len(matched_samples)} ({len(matched_samples)/len(df)*100:.1f}%)")
    print(f"   Unmatched: {len(unmatched_samples)} ({len(unmatched_samples)/len(df)*100:.1f}%)")
    print(f"   Report saved: {args.output}")
    
    if len(matched_samples) == 0:
        print("❌ No samples matched! Check sample naming conventions.")
        return 1
    elif len(unmatched_samples) > len(matched_samples):
        print("⚠️  More unmatched than matched samples. Review naming conventions.")
        return 1
    else:
        print("✅ Validation successful!")
        return 0

if __name__ == "__main__":
    exit(main())