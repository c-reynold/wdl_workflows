version 1.0
##########################################################################
# WORKFLOW: FQ2BAM
##########################################################################


workflow fq2bam {
    
    input {
        Array[File] r1_fastqs
        Array[File] r2_fastqs
        Int trim_length
        File reference_genome
    }

String sample_name = sub(basename(r1_fastqs[0]), "_.*$", "")
String bam = sample_name + ".bam"
String bam_idx = sample_name + ".bam.bai"

Int n = length(r1_fastqs)


    scatter (i in range(n)) {
        call trim_umi {
            input:
                r1_fastq = r1_fastqs[i],
                r2_fastq = r2_fastqs[i],
                trim_length = trim_length
}

    }

    call align_bam {
        input:
            r1_fastqs = trim_umi.trimmed_r1_fastq,
            r2_fastqs = trim_umi.trimmed_r2_fastq,
            reference_tarball = reference_genome,
            bam = bam,
            bam_idx = bam_idx
    }

    output {
        File bam_out = align_bam.bam_out
        File bam_idx_out = align_bam.bam_idx_out

    }
}

#######################################################################################################
#  TASK: TRIM FASTQ
#######################################################################################################

task trim_umi {

    input {

        File r1_fastq
        File r2_fastq
        Int trim_length

    }

    String r1_name = basename(r1_fastq, ".fastq.gz")
    String r2_name = basename(r2_fastq, ".fastq.gz")

    command <<<
       
        set -euo pipefail

        trimmomatic PE -phred33 ~{r1_fastq} ~{r2_fastq} ~{r1_name}.trimmed.fastq.gz ~{r1_name}.unpaired.fastq.gz ~{r2_name}.trimmed.fastq.gz ~{r2_name}.unpaired.fastq.gz \
        HEADCROP:~{trim_length}

    >>>

    runtime {
        docker: "quay.io/biocontainers/trimmomatic:0.40--hdfd78af_0"
        dx_instance_type: "mem1_ssd2_v2_x8"
    }

    output {
        
        File trimmed_r1_fastq = "~{r1_name}.trimmed.fastq.gz"
        File trimmed_r2_fastq = "~{r2_name}.trimmed.fastq.gz"
         
    }
}





#################################################################################################################
#  TASK: ALIGN BAM
#################################################################################################################

task align_bam {

    runtime {
    
    docker: "nvcr.io/nvidia/clara/clara-parabricks:4.2.0-1"
    dx_instance_type: "mem2_ssd2_gpu1_v2_x64"
    
    }

    input {

        Array[File] r1_fastqs
        Array[File] r2_fastqs
        File reference_tarball
        String bam
        String bam_idx

    }

    command <<<

        #!/usr/bin/env bash

        set -euo pipefail

        apt-get update
        apt-get install -y samtools


        R1_LINES="~{sep='\n' r1_fastqs}"
        R2_LINES="~{sep='\n' r2_fastqs}"

        mapfile -t r1_arr <<< "$R1_LINES"
        mapfile -t r2_arr <<< "$R2_LINES"


        if [[ "${#r1_arr[@]}" -ne "${#r2_arr[@]}" ]]; then
            echo "ERROR: R1 and R2 arrays differ in length" >&2
            exit 1
        fi

        #This loop builds the in_fq_args array input to fq2bam
        in_fq_args=()

        for ((i=0; i<${#r1_arr[@]}; i++)); do
            fq_base="$(basename "${r1_arr[i]}")"

            sample_name="${fq_base%%_S*}"

            #creates an arbitrary string for flowcell
            flowcell="flowcell"

            # # lane digits after _L (robust)
            # sample_lane="$(echo "$fq_base" | sed -n 's/.*_L\([0-9]\+\)_.*/\1/p')"
            # if [[ -z "$sample_lane" ]]; then
            #     echo "ERROR: could not parse lane from filename: $fq_base" >&2
            #     exit 1
            # fi

            # # optional: derive flowcell from header if not provided elsewhere
            # flowcell="$(zcat "${r1_arr[i]}" | head -n 1 | cut -d: -f3 | sed 's/^@//' || true)"

            rg="@RG\tID:${sample_name}.${i}\tPL:Illumina\tPU:${flowcell}.${i}\tSM:${sample_name}\tLB:default"

            in_fq_args+=( --in-fq "${r1_arr[i]}" "${r2_arr[i]}" "$rg" )

        done

        echo "Built fq2bam args:" >&2
        printf '  %q' "${in_fq_args[@]}" >&2
        echo >&2
        
        checkpoint=0    #Checkpoint initial value for progress and debugging


        #makes a directory to hold the reference files that are extracted later
        mkdir -p reference

        echo "checkpoint $checkpoint"
        echo "Directories created: reference and output"
        checkpoint=$((checkpoint+1))    


        #Extracts the reference tarballs
        tar -zxvf ~{reference_tarball} -C reference/

        echo "checkpoint $checkpoint"
        echo "Reference genome and index files extracted"
        checkpoint=$((checkpoint+1))    


        ref_fasta=$(find reference/ -type f \( -name '*.fa' -o -name '*.fasta' \) | head -n 1) #derives the reference fasta file name to be passed to --ref in fq2bam
        
       

        # echo "checkpoint $checkpoint"
        # echo "Input FASTQ arguments constructed for fq2bam"
        # checkpoint=$((checkpoint+1))
    

        pbrun fq2bam \
        --ref "$ref_fasta" \
        "${in_fq_args[@]}" \
        --bwa-options="-Y" \
        --out-bam "~{bam}" \
        --low-memory
        
        echo "checkpoint $checkpoint"
        echo "Alignment completed and BAM file generated"
        checkpoint=$((checkpoint+1))

        samtools index "~{bam}" "~{bam_idx}"

        echo "checkpoint $checkpoint"
        echo "BAM file indexed"
        checkpoint=$((checkpoint+1))
    >>>

    output {
        File bam_out = "~{bam}"
        File bam_idx_out = "~{bam_idx}"
    }

}