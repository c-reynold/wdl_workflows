#!/usr/bin/env bash

project="project-J2zPy600ZV7y1F26fbzv4vk9"
folder=/workflows/bam_utilities/samtools_downsample/


dxcompile compile "samtools_downsample/wdl/samtools_downsample.wdl" -project "$project" -folder "$folder" -archive