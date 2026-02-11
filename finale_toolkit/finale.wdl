version 1.0

workflow finale_end_motifs {
  input {
    File bam
    File bam_idx
    File reference_2bit
    Int cpu
    String? sample_name

    Int kmer_length
  }
  #Including a sample name is optionl. If there is a sample name, it will use that. If not, it will use the filename on the bam file itself for the sample name
  String resolved_sample_name=select_first([sample_name, basename(bam, ".bam")])

  #calculate the disk size so that we don't have any disk space problems on the instance. Conservitavly asking for 150% of the size of the bam file being processed
  Int disk_size = ceil(size(bam, "GB")*1.5)
  

  call end_motifs {
    input:
      bam = bam,
      bam_idx = bam_idx,
      resolved_sample_name = resolved_sample_name,
      disk_size = disk_size,
      cpu = cpu,
      reference_2bit = reference_2bit,
      kmer_length = kmer_length
  }

  output {
    File end_motif_tsv = end_motifs.end_motif_tsv
  }
}

task end_motifs {
  input {
    File bam
    File bam_idx
    String resolved_sample_name
    Int cpu
    Int disk_size
    File reference_2bit
    Int kmer_length
  }

  command <<<
    set -euo pipefail

    #transfering these variables to bash so that we can use them in output fields
   
    bam_idx="~{bam_idx}"
    resolved_sample_name="~{resolved_sample_name}"

    bam_base="$(basename "~{bam}")"        # e.g. tumor.bam
    expected_idx="${bam_base}.bai"       # tumor.bam.bai

    # If the expected index file isn't already present, copy it there.
    if [[ ! -f "$expected_idx" ]]; then
    cp -f "$bam_idx" "$expected_idx"
    fi

    echo "These are the local files"
    pwd
    ls

    finaletoolkit end-motifs "~{bam}" "~{reference_2bit}"\
        -q 30 \
        -k ~{kmer_length} \
        -o "${resolved_sample_name}.end_motif.tsv" \
        -w ~{cpu} \
        -v

    
    echo "These are all the final files"
    ls

  >>>

  runtime {
    docker: "dx://project-J2zPy600ZV7y1F26fbzv4vk9:/workflows/docker_images/finale_toolkit_0.11.tar.gz"
    docker: "dx://project-J2zPy600ZV7y1F26fbzv4vk9:/workflows/docker_images/fianle_toolkit_0.11.tar.gz"
    docker: "dx://project-J2zPy600ZV7y1F26fbzv4vk9:/workflows/docker_images/samtools_qc.tar.gz"

    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
  }

  output {
    File end_motif_tsv = "${resolved_sample_name}.end_motif.tsv"
  }
}
