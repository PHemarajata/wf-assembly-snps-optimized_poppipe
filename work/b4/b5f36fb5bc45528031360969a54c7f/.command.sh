#!/bin/bash -euo pipefail
set -x  # Enable debug mode
# Define bash functions inline
msg() {
  echo "[$(date '+%Y-%b-%d %a %H:%M:%S')] $@"
}

verify_minimum_file_size() {
  # Boolean test to ensure the filepath is a file, is non-zero size, and
  #  is at least the minimum specified size (in Bytes).
  # Parameters: $1=filename, $2=file description, $3=minimum size
  local file_path="${1}"
  local file_desc="${2}"
  local min_size="${3}"

  msg "DEBUG: Checking file: ${file_path}"
  msg "DEBUG: File exists: $(test -f "${file_path}" && echo 'YES' || echo 'NO')"
  msg "DEBUG: File size: $(ls -lh "${file_path}" 2>/dev/null | awk '{print $5}' || echo 'UNKNOWN')"
  msg "DEBUG: Minimum required: ${min_size}"

  if [ -f "${file_path}" ]; then
    if [ -s "${file_path}" ]; then
      # Use a more reliable method to check file size
      local find_result
      find_result=$(find -L "${file_path}" -type f -size +"${min_size}" 2>/dev/null)
      msg "DEBUG: Find result: '${find_result}'"

      if [ -n "${find_result}" ]; then
        msg "DEBUG: File ${file_path} passed size check"
        return 0
      else
        msg "ERROR: ${file_desc} file $(basename ${file_path}) present but too small (less than ${min_size})" >&2
        return 1
      fi
    else
      msg "ERROR: ${file_desc} file $(basename ${file_path}) present but empty" >&2
      return 1
    fi
  else
    msg "ERROR: ${file_desc} file $(basename ${file_path}) absent" >&2
    return 1
  fi
}

# Rename input files to prefix and move to inputfiles dir
mkdir inputfiles
cp ERS013358.fasta inputfiles/"1_ERS013358.fasta"

# gunzip all files that end in .{gz,Gz,GZ,gZ}
find -L inputfiles/ -type f -name '*.[gG][zZ]' -exec gunzip -f {} +

# Filter out small inputfiles
msg "Checking input file sizes for sample 1_ERS013358.."
echo -e "Sample name	QC step	Outcome (Pass/Fail)" > "1_ERS013358.Initial_Input_File.tsv"

file_count=0
for file in inputfiles/*; do
  ((file_count++))
  msg "Processing file ${file_count}: $(basename ${file})"

  # Check if file exists before processing
  if [ ! -f "${file}" ]; then
    msg "WARNING: File ${file} does not exist, skipping"
    continue
  fi

  msg "DEBUG: About to check file size for: ${file}"
  if verify_minimum_file_size "${file}" 'Input' "1k"; then
    msg "DEBUG: File ${file} passed size check"
    echo -e "$(basename ${file%%.*})	Input File	PASS"         >> "1_ERS013358.Initial_Input_File.tsv"
  else
    msg "DEBUG: File ${file} failed size check"
    echo -e "$(basename ${file%%.*})	Input File	FAIL"         >> "1_ERS013358.Initial_Input_File.tsv"
    rm "${file}"
  fi
done

msg "Completed processing ${file_count} files for sample 1_ERS013358"

cat <<-END_VERSIONS > versions.yml
"POPPIPE_SNP_ANALYSIS:CLUSTER_SNP_ANALYSIS:INFILE_HANDLING_UNIX":
  ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\n')
END_VERSIONS
