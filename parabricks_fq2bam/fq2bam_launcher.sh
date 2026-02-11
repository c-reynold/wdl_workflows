#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob  # so empty dir = no iteration

message_to_add_to_readme="These samples had their UMIs trimmed and two additional bases trimmed and then were aligned with parabricks fastq to bam using MRDs hg19 genome"


# api_key="$dx_api_key"

dx whoami
dx select --name "liquid_hrd"



workflow_id="workflow-J65gg2Q0ZV7ZXb4ZG5zvvZ97"
output_dir="liquid_hrd:/bam_files/aligned_to_MRD_hg19/$(date +%Y%m%d_%H%M%S)-umi+2"
json_directory="parabricks_fq2bam/fastq2bam_sample_jsons"

dx mkdir -p "$output_dir"

readme="readme.txt"
{
  echo "Parabricks Fastq-to-bam run with workflow ID: $workflow_id"
  echo "Launched on: $(date)"
  echo "$message_to_add_to_readme"
  echo "Samples and DNA nexus IDs:"
} > "$readme"

for file in "$json_directory"/*.json; do
  sample_name=$(basename "$file" .inputs.json)
  sample_dir="$output_dir/$sample_name"
  dx mkdir -p "$sample_dir"

  job_id=$(dx run "$workflow_id" \
    -f "$file" \
    --destination "$sample_dir" \
    --priority high \
    --brief \
    -y)

  dx upload "$file" --path "${sample_dir}/$(basename "$file")"

  echo "$sample_name: $job_id" | tee -a "$readme"
done

dx upload "$readme" --path "$output_dir/readme.txt"
rm "$readme"

