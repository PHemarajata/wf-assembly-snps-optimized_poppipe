# PopPUNK Compatibility Analysis

## Quick Answer: ✅ Your filenames should work fine!

After analyzing your file examples and the pipeline's PopPUNK matching logic, **your current filenames are actually quite good** and should work with the clustering workflow without renaming.

## Analysis Results

### Your Filename Quality: ✅ EXCELLENT
- **0 out of 10** sample files have problematic characters
- No spaces, excessive dots, or shell-problematic characters
- Version numbers (like `.1`, `.2`) are handled correctly
- Length and format are appropriate

### PopPUNK Matching Compatibility: ✅ HIGH
The pipeline uses **flexible matching** with multiple strategies:

1. **Direct matching**: Exact sample name match
2. **Suffix removal**: Removes `-SPAdes`, `_contigs`, `_genomic` automatically  
3. **Base name extraction**: Extracts core identifiers (e.g., `IP-0001-8` from `IP-0001-8_S1_L001-SPAdes`)
4. **Version handling**: Handles GCA/GCF accession versions properly
5. **Fuzzy matching**: Partial string matching as fallback

## How PopPUNK Matching Works

When you provide a PopPUNK cluster file, the pipeline will match samples using these patterns:

| Your FASTA File | PopPUNK Sample Name | Match Method |
|----------------|-------------------|--------------|
| `Burkholderia_pseudomallei_28P_GCF_025847415.1_Viet_Nam.fasta` | `Burkholderia_pseudomallei_28P_GCF_025847415.1_Viet_Nam` | Direct match ✅ |
| `Burkholderia_pseudomallei_28P_GCF_025847415.1_Viet_Nam.fasta` | `28P` | Base name extraction ✅ |
| `IP-0001-8_S1_L001-SPAdes.fasta` | `IP-0001-8_S1_L001-SPAdes` | Direct match ✅ |
| `IP-0001-8_S1_L001-SPAdes.fasta` | `IP-0001-8` | Base name extraction ✅ |
| `GCA_900595585.2.fasta` | `GCA_900595585.2` | Direct match ✅ |
| `GCA_900595585.2.fasta` | `GCA_900595585` | Version handling ✅ |
| `SRR12527871_contigs.fasta` | `SRR12527871_contigs` | Direct match ✅ |
| `SRR12527871_contigs.fasta` | `SRR12527871` | Suffix removal ✅ |

## Recommendations

### 🎯 **Primary Recommendation: Try as-is first**
Your filenames look good! Try running the pipeline with your current files first.

### 🛡️ **Backup Plan: Safe renaming available**
If you encounter any issues, I've provided a `safe_rename_fasta_files.py` script that:
- Makes **minimal changes** to preserve PopPUNK compatibility
- Only fixes truly problematic characters
- Maintains all version numbers and identifiers
- Preserves sample name matching

### 📋 **PopPUNK Cluster File Format**
Make sure your PopPUNK cluster file has this format:
```csv
sample,cluster
Burkholderia_pseudomallei_28P_GCF_025847415.1_Viet_Nam,1
IP-0001-8_S1_L001-SPAdes,2
GCA_900595585.2,1
...
```

## Pipeline Improvements Made

The repository includes fixes to the filename parsing logic:

1. **Better extension detection** in `modules/local/infile_handling_unix/main.nf`
2. **Improved sample ID generation** in `subworkflows/local/input_check.nf`
3. **Flexible PopPUNK matching** in `modules/local/parse_poppunk_clusters/main.nf`

## Conclusion

✅ **Your filenames should work fine with PopPUNK clustering!**

The pipeline has robust matching logic that handles your naming patterns well. The fixes I made to the parsing modules should resolve the original filename recognition issues you encountered.

**Next steps:**
1. Download the cleaned repository
2. Try running with your current filenames
3. Use the safe renaming script only if needed