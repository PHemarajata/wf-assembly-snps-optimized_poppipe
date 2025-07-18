# Relative Path Support in TSV Files

## Overview

The pipeline now supports relative paths in TSV files (such as `rfiles.txt` or `combined_rfile.txt`). This feature allows you to organize your project files in a more flexible way and makes your analysis more portable.

## How It Works

When the pipeline reads a TSV file containing sample-to-file mappings, it now:

1. **Detects the TSV file location**: The pipeline identifies the directory containing the TSV file
2. **Resolves relative paths**: Any relative paths in the TSV file are resolved relative to the TSV file's directory
3. **Normalizes paths**: Paths are normalized to handle `..` and `.` components correctly
4. **Validates files**: The pipeline verifies that all resolved file paths exist

## Example Project Structure

```
my_project/
├── data/
│   ├── rfiles.txt              # TSV file with sample mappings
│   └── metadata.csv
├── assemblies/
│   ├── batch1/
│   │   ├── sample1.fasta
│   │   └── sample2.fasta
│   └── batch2/
│       ├── sample3.fasta
│       └── sample4.fasta
└── poppipe_output/
    ├── combined_clusters.csv
    └── strain/
        ├── cluster_1/
        └── cluster_2/
```

## TSV File Format with Relative Paths

In the above structure, your `data/rfiles.txt` file can use relative paths:

```tsv
sample1	../assemblies/batch1/sample1.fasta
sample2	../assemblies/batch1/sample2.fasta
sample3	../assemblies/batch2/sample3.fasta
sample4	../assemblies/batch2/sample4.fasta
```

## Path Resolution Examples

| TSV File Location | Path in TSV | Resolved Path |
|-------------------|-------------|---------------|
| `/project/data/rfiles.txt` | `../assemblies/sample1.fasta` | `/project/assemblies/sample1.fasta` |
| `/project/data/rfiles.txt` | `./local_file.fasta` | `/project/data/local_file.fasta` |
| `/project/data/rfiles.txt` | `/absolute/path/file.fasta` | `/absolute/path/file.fasta` |

## Benefits

1. **Portability**: Projects can be moved between systems without updating file paths
2. **Organization**: Keep TSV files separate from assembly files while maintaining references
3. **Flexibility**: Mix relative and absolute paths as needed
4. **Clarity**: Relative paths make project structure more obvious

## Usage

### Command Line Example

```bash
nextflow run main.nf \
  -profile docker \
  --poppipe_output /path/to/poppipe_output \
  --input /path/to/data/rfiles.txt \
  --outdir results
```

### TSV File Requirements

- Tab-separated format
- First column: sample name
- Second column: file path (can be relative or absolute)
- Comments (lines starting with `#`) are ignored
- Empty lines are ignored

### Supported File Extensions

The pipeline will look for assembly files with these extensions:
- `.fasta`, `.fas`, `.fna`, `.fsa`, `.fa`
- Gzip compressed versions: `.fasta.gz`, `.fas.gz`, etc.

## Testing

You can test the relative path functionality using the provided test script:

```bash
nextflow run test_relative_paths.nf
```

This will create a test project structure and verify that relative paths are resolved correctly.

## Troubleshooting

### Common Issues

1. **File not found errors**: Check that relative paths are correct relative to the TSV file location
2. **Path resolution**: Use `../` to go up directories, `./` for current directory
3. **Mixed separators**: Use forward slashes `/` even on Windows (Nextflow handles conversion)

### Debug Information

The pipeline provides detailed logging showing:
- TSV file directory location
- Original paths from TSV file
- Resolved absolute paths
- File existence verification

Check the process logs for detailed path resolution information.