process INFILE_HANDLING_UNIX {

    tag( "${meta.id}" )
    container "ubuntu:jammy"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${meta.id}.Initial_Input_File.tsv"), emit: qc_filecheck
    path("inputfiles/*")                                      , emit: input_files
    path(".command.{out,err}")
    path("versions.yml")                                      , emit: versions

    shell:
    // Rename files with meta.id (has spaces and periods removed)
    gzip_compressed = input.toString().contains('.gz') ? '.gz' : ''
    // Improved file extension extraction to handle complex filenames
    def filename_no_gz = input.toString().replaceAll(/(?i)\.gz$/, '')
    def known_extensions = ['fasta', 'fas', 'fna', 'fsa', 'fa']
    file_extension = known_extensions.find { ext -> 
        filename_no_gz.toLowerCase().endsWith('.' + ext) 
    } ?: filename_no_gz.split('\\.')[-1]
    '''
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
    cp !{input} inputfiles/"!{meta.id}.!{file_extension}!{gzip_compressed}"

    # gunzip all files that end in .{gz,Gz,GZ,gZ}
    find -L inputfiles/ -type f -name '*.[gG][zZ]' -exec gunzip -f {} +

    # Filter out small inputfiles
    msg "Checking input file sizes for sample !{meta.id}.."
    echo -e "Sample name\tQC step\tOutcome (Pass/Fail)" > "!{meta.id}.Initial_Input_File.tsv"
    
    file_count=0
    for file in inputfiles/*; do
      ((file_count++))
      msg "Processing file ${file_count}: $(basename ${file})"
      
      if verify_minimum_file_size "${file}" 'Input' "!{params.min_input_filesize}"; then
        echo -e "$(basename ${file%%.*})\tInput File\tPASS" \
        >> "!{meta.id}.Initial_Input_File.tsv"
      else
        echo -e "$(basename ${file%%.*})\tInput File\tFAIL" \
        >> "!{meta.id}.Initial_Input_File.tsv"
        rm ${file}
      fi
    done
    
    msg "Completed processing ${file_count} files for sample !{meta.id}"

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
      ubuntu: $(awk -F ' ' '{print $2,$3}' /etc/issue | tr -d '\\n')
    END_VERSIONS
    '''
}
