#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob  # so empty dir = no iteration

message_to_add_to_readme="These are QC metrics made using samtools flagstat, stats, and idxstats commands."



# api_key="$dx_api_key"

dx whoami
dx select --name "liquid_hrd"



workflow_id="workflow-J617j1j0ZV7gB4gjJpqK671k"
output_dir="liquid_hrd:/bam_files/aligned_to_MRD_hg19/qc_metrics/$(date +%Y%m%d_%H%M%S)"
json_directory="bam_qc/input_jsons"

dx mkdir -p "$output_dir"

readme="readme.txt"
{
  echo "samtools qc based on samtools v1.13: $workflow_id"
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

