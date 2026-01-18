version 1.0
##########################################################################
# WORKFLOW: FQ2BAM
##########################################################################

struct fastq_pair {
    File r1_fastq
    File r2_fastq
}

workflow fq2bam {
    
    input {
        Array[fastq_pair] fastq_pairs
        Int umi_length
        File reference_genome
        File bwa_index_tarball
        File gatk_index_tarball
    }

    scatter (p in fastq_pairs) {
        call trim_umi {
            input:
                untrimmed_fastqs = p,
                umi_length = umi_length
        }
    }

    call align_bam {
        input:
            r1_fastqs = trim_umi.r1_fastqs,
            r2_fastqs = trim_umi.r2_fastqs,
            reference_tarball = reference_genome,
            bwa_index_tarball = bwa_index_tarball,
            gatk_index_tarball = gatk_index_tarball
    }

    output {
        File bam = align_bam.output_bam

    }
}

#######################################################################################################
#  TASK: TRIM FASTQ
#######################################################################################################

task trim_umi {

    input {

        fastq_pair untrimmed_fastqs
        Int umi_length

    }

    String r1_name = basename(untrimmed_fastqs.r1_fastq, ".fastq.gz")
    String r2_name = basename(untrimmed_fastqs.r2_fastq, ".fastq.gz")

    command <<<
       
        set -euo pipefail

        trimmomatic PE ~{untrimmed_fastqs.r1_fastq} ~{untrimmed_fastqs.r2_fastq} ~{r1_name}.trimmed.fastq.gz ~{r1_name}.unpaired.fastq.gz ~{r2_name}.trimmed.fastq.gz ~{r2_name}.unpaired.fastq.gz \
        HEADCROP:~{umi_length}

    >>>

    runtime {
        docker: "quay.io/biocontainers/trimmomatic:0.40--hdfd78af_0"
        dx_instance_type: "mem1_ssd1_v2_x8"
    }

    output {
        
        File r1_fastqs = "~{r1_name}.trimmed.fastq.gz"
        File r2_fastqs = "~{r2_name}.trimmed.fastq.gz"
         
    }
}





#################################################################################################################
#  TASK: ALIGN BAM
#################################################################################################################

task align_bam {

    runtime {
    
    docker: "nvcr.io/nvidia/clara/clara-parabricks:4.2.0-1"
    dx_instance_type: "mem2_ssd1_gpu_x48"
    
    }

    input {

        Array[File] r1_fastqs
        Array[File] r2_fastqs
        File reference_tarball
        File bwa_index_tarball
        File gatk_index_tarball

    }

    command <<<

        #!/usr/bin/env bash

        set -euo pipefail

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
        in_fq_args+=( --in-fq "${r1_arr[i]}" "${r2_arr[i]}" )
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


        #Extracts the three reference tarballs (can probably change this to use one tarballed reference with all the needed indicies in it)
        tar -zxvf ~{reference_tarball} -C reference/
        tar -zxvf ~{bwa_index_tarball} -C reference/
        tar -zxvf ~{gatk_index_tarball} -C reference/

        echo "checkpoint $checkpoint"
        echo "Reference genome and index files extracted"
        checkpoint=$((checkpoint+1))    


        ref_fasta=$(find reference/ -name '*.fa' | head -n 1)  #derives the reference fasta file name to be passed to --ref in fq2bam

        sample_name=$(basename "${r1_arr[0]}" | awk -F"_" '{print $1}')  #Derives the sample name from the first trimmed R1 FASTQ file, assuming the naming convention is SampleName_Lane_Read.fastq.gz
        
       

        # echo "checkpoint $checkpoint"
        # echo "Input FASTQ arguments constructed for fq2bam"
        # checkpoint=$((checkpoint+1))
    

        pbrun fq2bam \
        --ref "$ref_fasta" \
        "${in_fq_args[@]}" \
        --bwa-options="-Y" \
        --out-bam "output.bam" \
        --low-memory
        
        echo "checkpoint $checkpoint"
        echo "Alignment completed and BAM file generated"
        checkpoint=$((checkpoint+1))

    >>>

    output {
        File output_bam = "output.bam"
    }

}








