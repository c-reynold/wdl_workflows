#!/usr/bin/env bash

miniwdl run "samtools_downsample/wdl/samtools_downsample.wdl" -i "samtools_downsample/mini_wdl/input_jsons/ichor_cna_ichorcna_inputs.json" -d "samtools_downsample/mini_wdl/outputs" --verbose
