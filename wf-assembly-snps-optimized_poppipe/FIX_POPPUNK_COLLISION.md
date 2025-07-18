# Fix for PopPUNK File Name Collision Error

## Problem
The `PARSE_POPPUNK_CLUSTERS` process was failing with a "file name collision" error because Nextflow was trying to stage multiple files with the same basename into the same working directory.

## Error Message
```
Process `ASSEMBLY_SNPS_CLUSTERED:PARSE_POPPUNK_CLUSTERS` input file name collision -- There are multiple input files for each of the following file names: [long list of files]
```

## Root Cause
The issue occurred because:
1. The process was receiving all input files via `path input_files` 
2. Nextflow stages all `path` inputs into the process working directory
3. If multiple files have the same basename (even from different directories), Nextflow sees this as a collision
4. This can happen when you have duplicate filenames or when the same files are being processed multiple times

## Solution
Modified the `PARSE_POPPUNK_CLUSTERS` process to:

1. **Changed input type**: From `path input_files` to `val input_file_list`
2. **Pass file paths as strings**: Instead of staging files, we now pass file paths as a list of strings
3. **Added file validation**: The script now checks if files exist and filters for FASTA files
4. **Eliminated staging conflicts**: No files are staged, so no collision can occur

## Changes Made

### In `modules/local/parse_poppunk_clusters/main.nf`:
- Changed input from `path input_files` to `val input_file_list`
- Updated script to work with file paths instead of staged files
- Added file existence validation

### In `workflows/assembly_snps_clustered.nf`:
- Modified the channel to pass file paths as strings using `.toString()`
- Added `.unique()` to remove any duplicate paths
- Updated process call to use the new input format

## Benefits
1. **Eliminates file staging collisions**
2. **Handles duplicate filenames gracefully**
3. **More efficient** (no unnecessary file staging)
4. **More robust** (validates file existence)

## Testing
After applying this fix, the pipeline should be able to handle:
- Files with duplicate basenames from different directories
- Large numbers of input files
- Mixed file naming patterns

The PopPUNK clustering workflow should now run without the collision error.