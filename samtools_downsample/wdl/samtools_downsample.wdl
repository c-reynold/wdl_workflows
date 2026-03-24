version 1.0
task samtools_downsample {
  input {
    File bam
    File bai
    String sample_name
    Int cpu
    Int disk_size

    Float current_coverage
    Float desired_coverage
  }

    Float downsample_factor = desired_coverage/current_coverage

    Int desired_coverage_round = round(desired_coverage)
  command <<<
    
    samtools view -@ $(nproc) -s ~{downsample_factor} -b "~{bam}" > "~{sample_name}_depth_~{desired_coverage_round}x.bam"
    samtools index "~{sample_name}_depth_~{desired_coverage_round}x.bam" "~{sample_name}_depth_~{desired_coverage_round}x.bam.bai"

  >>>

  runtime {
    docker: "dx://project-J2zPy600ZV7y1F26fbzv4vk9:/workflows/docker_images/samtools_qc.tar.gz" # for miniwdl testing, can use "biocontainers/samtools:v1.9-4-deb_cv1"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
  }

  output {
    File bam_out = "~{sample_name}_depth_~{desired_coverage_round}x.bam"
    File bai_out = "~{sample_name}_depth_~{desired_coverage_round}x.bam.bai"
  }
}
