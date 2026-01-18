version 1.0
##########################################################################
# WORKFLOW: FQ2BAM
##########################################################################
workflow fq2bam {
    
    input {
        Array[File] r1_fastqs
        Array[File] r2_fastqs
        Int umi_length
        File reference_genome
        File bwa_index_tarball
        File gatk_index_tarball
    }

    scatter (i in range(length(r1_fastqs))) {
        call trim_umi {
            input:
                r1_file = r1_fastqs[i],
                r2_file = r2_fastqs[i],
                umi_length = umi_length
        }
    }

    call align_bam {
        input:
            trimmed_r1_fastq = trim_umi.trimmed_r1_fastq,
            trimmed_r2_fastq = trim_umi.trimmed_r2_fastq,
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

        File r1_file
        File r2_file
        Int umi_length

    }

    String r1_name = basename(r1_file, ".fastq.gz")
    String r2_name = basename(r2_file, ".fastq.gz")

    command <<<
       
        set -euo pipefail

        trimmomatic PE ~{r1_file} ~{r2_file} ~{r1_name}.trimmed.fastq.gz ~{r1_name}.unpaired.fastq.gz ~{r2_name}.trimmed.fastq.gz ~{r2_name}.unpaired.fastq.gz \
        HEADCROP:~{umi_length}

    >>>

    runtime {
        docker: "quay.io/biocontainers/trimmomatic:0.40--hdfd78af_0"
        dx_instance_type: "mem1_ssd1_v2_x8"
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
    dx_instance_type: "mem2_ssd1_gpu_x48"
    
    }

    input {

        Array[File] trimmed_r1_fastq
        Array[File] trimmed_r2_fastq
        File reference_tarball
        File bwa_index_tarball
        File gatk_index_tarball

    }


    command <<<

        set -euo pipefail

        R1_LIST=$'~{sep="\n" trimmed_r1_fastq}'
        R2_LIST=$'~{sep="\n" trimmed_r2_fastq}'

        mapfile -t r1_fastqs <<< "$R1_LIST"
        mapfile -t r2_fastqs <<< "$R2_LIST"


        #This checks to make sure that the number of fastq files in R1 and R2 arrays match
        if [[ "${#r1_fastqs[@]}" -ne "${#r2_fastqs[@]}" ]]; then
        echo "ERROR: R1 and R2 arrays differ in length" >&2
        exit 1
        fi

        #This loop builds the in_fq_args array input to fq2bam
        in_fq_args=()
        for ((i=0; i<${#r1_fastqs[@]}; i++)); do
        in_fq_args+=( --in-fq "${r1_fastqs[i]}" "${r2_fastqs[i]}" )
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

        sample_name=$(basename "~{trimmed_r1_fastq[0]}" | awk -F"_" '{print $1}')  #Derives the sample name from the first trimmed R1 FASTQ file, assuming the naming convention is SampleName_Lane_Read.fastq.gz
        
       

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








