#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob  # so empty dir = no iteration



# api_key="$dx_api_key"

dx whoami
dx select --name "liquid_hrd"



workflow_id="workflow-J6p7vB80ZV7gGqBq7qQP6b7k"
output_dir="liquid_hrd:/bam_files/aligned_to_MRD_hg19/trim_umi+2"
json_directory="parabricks_fq2bam/fastq2bam_sample_jsons"


for file in "$json_directory"/*.json; do
  sample_name=$(basename "$file" .inputs.json)
  sample_dir="$output_dir/$sample_name"
 

  job_id=$(dx run "$workflow_id" \
    -f "$file" \
    --destination "$sample_dir" \
    --priority high \
    --brief \
    --name "fq2bam_${sample_name}" \
    --instance-type align_bam=mem2_ssd2_gpu1_v2_x64 \
    -y)

  echo "Launched job $job_id for sample $sample_name"

done


