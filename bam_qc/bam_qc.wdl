version 1.0

workflow bam_qc {
  input {
    File bam
    File bam_idx
    Int cpu
  }

  String sample_name = basename(bam, ".bam")

  Int disk_size = ceil(size(bam, "GB")*1.5)
  

  call samtools_qc {
    input:
      bam = bam,
      bam_idx = bam_idx,
      sample_name = sample_name,
      disk_size = disk_size,
      cpu = cpu
  }

  output {
    File flagstat = samtools_qc.flagstat
    File idxstats = samtools_qc.idxstats
    File stats    = samtools_qc.stats

    File? acgt_cycles_data  = samtools_qc.acgt_cycles_data
    File? acgt_cycles_graph = samtools_qc.acgt_cycles_graph

    File? coverage_data  = samtools_qc.coverage_data
    File? coverage_graph = samtools_qc.coverage_graph

    File? gc_content_data  = samtools_qc.gc_content_data
    File? gc_content_graph = samtools_qc.gc_content_graph

    File? gc_depth_data  = samtools_qc.gc_depth_data
    File? gc_depth_graph = samtools_qc.gc_depth_graph

    File? quals_hm_data  = samtools_qc.quals_hm_data
    File? quals_hm_graph = samtools_qc.quals_hm_graph

    File? quals_data  = samtools_qc.quals_data
    File? quals_graph = samtools_qc.quals_graph
    File? quals2_data = samtools_qc.quals2_data
    File? quals2_graph = samtools_qc.quals2_graph
    File? quals3_data = samtools_qc.quals3_data
    File? quals3_graph = samtools_qc.quals3_graph

    File? read_length_hist      = samtools_qc.read_length_hist
    File? read_length_raw_hist  = samtools_qc.read_length_raw_hist
    File? read_length_data      = samtools_qc.read_length_data
    File? read_length_graph     = samtools_qc.read_length_graph
  }
}

task samtools_qc {
  input {
    File bam
    File bam_idx
    String sample_name
    Int cpu
    Int disk_size
  }

  command <<<
    set -euo pipefail

    bam="~{bam}"
    bam_idx="~{bam_idx}"

    bam_base="$(basename "$bam")"        # e.g. tumor.bam
    expected_idx="${bam_base}.bai"       # tumor.bam.bai

    # If the expected index file isn't already present, copy it there.
    if [[ ! -f "$expected_idx" ]]; then
    cp -f "$bam_idx" "$expected_idx"
    fi

    echo "These are the local files"
    pwd
    ls

    sample_name="~{sample_name}"

    flagstat_file="${sample_name}.flagstat.txt"
    idxstats_file="${sample_name}.idx_stats.txt"
    bamstats_file="${sample_name}.bam_stats.txt"

    touch "$flagstat_file" "$idxstats_file" "$bamstats_file"

    # -O tsv is not supported on all samtools versions; remove unless you have confirmed it.
    samtools flagstat "$bam" > "$flagstat_file"

    samtools idxstats "$bam" > "$idxstats_file"
    samtools stats "$bam" > "$bamstats_file"

    plot-bamstats -p "$sample_name" "$bamstats_file"

    echo "These are all the final files"
    ls

  >>>

  runtime {
    docker: "dx://project-J2zPy600ZV7y1F26fbzv4vk9:/workflows/docker_images/samtools_qc.tar.gz"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
  }

  output {
    File flagstat = "${sample_name}.flagstat.txt"
    File idxstats = "${sample_name}.idx_stats.txt"
    File stats    = "${sample_name}.bam_stats.txt"

    # plot-bamstats outputs (directory-based). These names are what plot-bamstats typically uses.
    File? acgt_cycles_data  = "${sample_name}-acgt-cycles.gp"
    File? acgt_cycles_graph = "${sample_name}-acgt-cycles.png"

    File? coverage_data  = "${sample_name}-coverage.gp"
    File? coverage_graph = "${sample_name}-coverage.png"

    File? gc_content_data  = "${sample_name}-gc-content.gp"
    File? gc_content_graph = "${sample_name}-gc-content.png"

    File? gc_depth_data  = "${sample_name}-gc-depth.gp"
    File? gc_depth_graph = "${sample_name}-gc-depth.png"

    File? quals_hm_data  = "${sample_name}-quals-hm.gp"
    File? quals_hm_graph = "${sample_name}-quals-hm.png"

    File? quals_data  = "${sample_name}-quals.gp"
    File? quals_graph = "${sample_name}-quals.png"
    File? quals2_data = "${sample_name}-quals2.gp"
    File? quals2_graph = "${sample_name}-quals2.png"
    File? quals3_data = "${sample_name}-quals3.gp"
    File? quals3_graph = "${sample_name}-quals3.png"

    File? read_length_hist     = "${sample_name}-read_len_histdata.txt"
    File? read_length_raw_hist = "${sample_name}-read_len_rawdata.txt"
    File? read_length_data     = "${sample_name}-read_length.gp"
    File? read_length_graph    = "${sample_name}-read_length.png"
  }
}
