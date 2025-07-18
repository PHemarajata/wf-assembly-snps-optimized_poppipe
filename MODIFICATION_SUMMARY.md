# PopPIPE-bp Integration Modifications Summary

This document summarizes all modifications made to integrate PopPIPE-bp outputs with the wf-assembly-snps pipeline.

## New Files Created

### 1. Modules

#### `modules/local/parse_poppipe_output/main.nf`
- **Purpose**: Parse PopPIPE-bp output directory and cluster assignments
- **Functionality**:
  - Finds and reads cluster assignment CSV files
  - Maps sample names to assembly files with fuzzy matching
  - Creates cluster-specific file lists
  - Generates cluster summary and metadata

#### `modules/local/tree_graft/main.nf`
- **Purpose**: Implement PopPIPE-style tree grafting algorithm
- **Functionality**:
  - Loads cluster-specific ML trees from SNP analysis
  - Loads PopPIPE NJ trees for branch length scaling
  - Scales ML tree branches based on NJ tree total lengths
  - Grafts scaled ML trees onto backbone tree structure
  - Outputs integrated phylogenetic tree

### 2. Workflows

#### `workflows/poppipe_snp_analysis.nf`
- **Purpose**: Main workflow for PopPIPE-bp integration
- **Functionality**:
  - Orchestrates the entire PopPIPE integration process
  - Parses PopPIPE outputs
  - Runs cluster-based SNP analysis
  - Performs tree grafting
  - Generates comprehensive outputs

### 3. Documentation

#### `POPPIPE_INTEGRATION.md`
- Comprehensive guide for using the PopPIPE integration
- Usage examples and parameter descriptions
- Troubleshooting guide
- Expected input/output formats

#### `test_poppipe_integration.nf`
- Test script to validate the integration
- Creates mock data for testing
- Validates parsing functionality

#### `MODIFICATION_SUMMARY.md` (this file)
- Summary of all changes made

## Modified Files

### 1. Core Pipeline Files

#### `main.nf`
- **Changes**:
  - Added import for `POPPIPE_SNP_ANALYSIS` workflow
  - Modified main workflow logic to check for `params.poppipe_output`
  - Added conditional execution: PopPIPE → PopPUNK clusters → Standard

#### `conf/params.config`
- **Changes**:
  - Added `poppipe_output` parameter (default: null)
  - Added `enable_tree_grafting` parameter (default: true)
  - Existing `min_cluster_size` parameter now applies to PopPIPE clusters

#### `nextflow_schema.json`
- **Changes**:
  - Added schema definitions for new parameters:
    - `poppipe_output`: Path to PopPIPE-bp output directory
    - `min_cluster_size`: Minimum cluster size for analysis
    - `enable_tree_grafting`: Boolean to enable/disable tree grafting

#### `README.md`
- **Changes**:
  - Updated introduction to mention PopPIPE integration
  - Added PopPIPE usage example
  - Added feature list highlighting new capabilities

## Integration Architecture

### Workflow Decision Tree

```
main.nf
├── params.poppipe_output? 
│   ├── YES → POPPIPE_SNP_ANALYSIS()
│   └── NO → params.poppunk_clusters?
│       ├── YES → ASSEMBLY_SNPS_CLUSTERED()
│       └── NO → ASSEMBLY_SNPS()
```

### PopPIPE Integration Flow

```
1. PARSE_POPPIPE_OUTPUT
   ├── Read cluster CSV files
   ├── Map samples to assembly files
   ├── Create cluster file lists
   └── Generate metadata

2. CLUSTER_SNP_ANALYSIS (per cluster)
   ├── Core genome alignment (ParSNP)
   ├── SNP distance calculation
   ├── Recombination detection (Gubbins)
   ├── Phylogenetic tree construction
   └── Distance matrix generation

3. TREE_GRAFT
   ├── Load cluster ML trees
   ├── Load PopPIPE NJ trees
   ├── Scale branch lengths
   ├── Graft onto backbone
   └── Output integrated tree
```

## Key Features Implemented

### 1. Flexible Input Parsing
- **TSV File Input (Recommended)**: Direct use of PopPIPE rfiles.txt or similar TSV files with sample-to-file mapping
- **Directory Input**: Supports various PopPIPE output directory structures
- Handles different CSV column naming conventions
- Fuzzy matching for sample name variations (directory input only)
- Robust error handling for missing files

### 2. PopPIPE Algorithm Implementation
- Faithful implementation of the tree grafting algorithm from PopPIPE-bp
- Branch length scaling based on NJ tree total lengths
- Proper handling of backbone trees and cluster positioning
- Midpoint rooting of final integrated tree

### 3. Scalable Cluster Processing
- Parallel processing of clusters
- Configurable minimum cluster size
- Quality control checks per cluster
- Comprehensive output collection

### 4. Comprehensive Output Generation
- Cluster-specific SNP analysis results
- Integrated phylogenetic tree
- Distance matrices and alignments
- Summary reports and logs
- Quality control information

## Usage Scenarios

### 1. Basic PopPIPE Integration
```bash
nextflow run main.nf \
  -profile docker \
  --poppipe_output poppipe_results/ \
  --input assemblies/ \
  --outdir results/
```

### 2. With Custom Parameters
```bash
nextflow run main.nf \
  -profile docker \
  --poppipe_output poppipe_results/ \
  --input assemblies/ \
  --ref reference.fasta \
  --min_cluster_size 5 \
  --outdir results/
```

## Expected Input Structure

### PopPIPE-bp Output Directory
```
poppipe_output/
├── combined_clusters.csv          # Required: cluster assignments
├── strain/                        # Optional: strain-specific data
│   ├── cluster_1/
│   │   ├── njtree.nwk            # For branch scaling
│   │   └── ...
│   └── other/
│       └── njtree.nwk            # Backbone tree
└── full_tree.nwk                 # Optional: overall tree
```

### Assembly Directory
```
assemblies/
├── sample1.fasta
├── sample2.fasta
├── sample3.fasta
└── ...
```

## Output Structure

```
results/
├── grafted_tree.nwk              # Main integrated phylogeny
├── cluster_summary.tsv           # Cluster information
├── tree_graft_log.txt           # Grafting process log
├── cluster_1/                   # Per-cluster results
│   ├── distance_matrix.tsv
│   ├── phylogenetic_tree.nwk
│   └── ...
├── cluster_2/
│   └── ...
└── Summaries/
    └── Summary.QC_File_Checks.tsv
```

## Testing and Validation

### Test Script
- `test_poppipe_integration.nf` provides validation of core functionality
- Creates mock PopPIPE output and assembly data
- Tests parsing and basic workflow execution

### Quality Control
- File size and format validation
- Sample name matching verification
- Cluster size filtering
- Tree grafting success logging

## Benefits of Integration

1. **Scalability**: Process large datasets by analyzing clusters separately
2. **Accuracy**: Detailed SNP analysis within homogeneous clusters
3. **Integration**: Combine results using established PopPIPE algorithms
4. **Flexibility**: Works with various PopPIPE output configurations
5. **Compatibility**: Maintains existing pipeline functionality

## Future Enhancements

Potential areas for future development:
1. Support for additional SNP calling methods
2. Integration with other clustering algorithms
3. Enhanced visualization outputs
4. Performance optimizations for very large datasets
5. Additional tree grafting algorithms

## Dependencies

### New Dependencies Added
- `ete3`: For tree manipulation and grafting
- `pandas`: For CSV parsing and data manipulation
- `json`: For metadata handling

### Existing Dependencies Used
- All existing wf-assembly-snps dependencies
- ParSNP for SNP calling
- Gubbins for recombination detection

## Compatibility

- **Nextflow**: DSL2 compatible
- **Containers**: Docker and Singularity support
- **Platforms**: Linux, macOS (with appropriate containers)
- **PopPIPE versions**: Compatible with PopPIPE-bp output formats

This integration maintains full backward compatibility with the existing wf-assembly-snps pipeline while adding powerful new capabilities for PopPIPE-bp users.