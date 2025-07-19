#!/bin/bash -euo pipefail
# Define bash functions inline
msg() {
  echo "[$(date '+%Y-%b-%d %a %H:%M:%S')] $@"
}

verify_minimum_file_size() {
  # Boolean test to ensure the filepath is a file, is non-zero size, and
  #  is at least the minimum specified size (in Bytes).
  # Parameters: $1=filename, $2=file description, $3=minimum size
  if [ -f "${1}" ]; then
    if [ -s "${1}" ]; then
      if [[ $(find -L "${1}" -type f -size +"${3}") ]]; then
        return 0
      else
        msg "ERROR: ${2} file $(basename ${1}) present but too small (less than ${3})" >&2
        return 1
      fi
    else
      msg "ERROR: ${2} file $(basename ${1}) present but empty" >&2
      return 1
    fi
  else
    msg "ERROR: ${2} file $(basename ${1}) absent" >&2
    return 1
  fi
}

# Rename input files to prefix and move to inputfiles dir
mkdir inputfiles
cp IP-0162.fasta inputfiles/"1_IP-0162.fasta"

# gunzip all files that end in .{gz,Gz,GZ,gZ}
find -L inputfiles/ -type f -name '*.[gG][zZ]' -exec gunzip -f {} +

# Filter out small inputfiles
msg "Checking input file sizes for sample 1_IP-0162.."
echo -e "Sample name	QC step	Outcome (Pass/Fail)" > "1_IP-0162.Initial_Input_File.tsv"

file_count=0
for file in inputfiles/*; do
  # Check if the glob matched any files
  if [ ! -e "${file}" ]; then
    msg "No files found in inputfiles directory"
    break
  fi

  ((file_count++))

  # Safely get basename
  local filename
  filename=$(basename "${file}" 2>/dev/null || echo "unknown_file")
  msg "Processing file ${file_count}: ${filename}"

  if verify_minimum_file_size "${file}" 'Input' "1k"; then
    local base_name
    base_name=$(basename "${file%%.*}" 2>/dev/null || echo "unknown")
    echo -e "${base_name}	Input File	PASS"         >> "1_IP-0162.Initial_Input_File.tsv"
  else
    local base_name
    base_name=$(basename "${file%%.*}" 2>/dev/null || echo "unknown")
    echo -e "${base_name}	Input File	FAIL"         >> "1_IP-0162.Initial_Input_File.tsv"
    rm "${file}"
  fi
done

msg "Completed processing ${file_count} files for sample 1_IP-0162"

cat <<-END_VERSIONS > versions.yml
"POPPIPE_SNP_ANALYSIS:CLUSTER_SNP_ANALYSIS:INFILE_HANDLING_UNIX":
  ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\n')
END_VERSIONS
