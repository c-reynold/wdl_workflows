version 1.0

import "mutect2.wdl" as mutect2

workflow pb_mutect2 {
    input {
        File genome_tarball
        File tumor_bam
        File? normal_bam
        File? tumor_bai
        File? normal_bai
        String sample_name
    }

    call pb_m2 {
        input:
            tumor_bam = tumor_bam,
            normal_bam = normal_bam,
            tumor_bai = tumor_bai,
            normal_bai = normal_bai,
            genome_tarball = genome_tarball,
    }

    call mutect2.filter_calls {
        input:
            genome_tarball = genome_tarball,
            vcf = pb_m2.somatic_vcf,
            vcf_idx = pb_m2.somatic_vcf_idx,
            stats = pb_m2.stats,
            sample_name = sample_name,
    }

    output {
        File filtered_vcf = filter_calls.filtered_vcf
        File filtering_stats = filter_calls.filtering_stats
    }
}

task pb_m2 {
    input {
        File tumor_bam
        File? normal_bam
        File? tumor_bai
        File? normal_bai
        File genome_tarball
    }
    Boolean tumor_bai_defined = defined(tumor_bai)
    Boolean normal_bai_defined = defined(normal_bai)
    Boolean normal_bam_defined = defined(normal_bam)
    command <<<
        LD_LIBRARY_PATH='/usr/local/cuda/compat:/usr/local/cuda/lib64'

        set -exo pipefail

        # Show parabricks mutectcaller options
        pbrun mutectcaller -h

        mkdir /tmp/genome/
        tar -zxvf ~{genome_tarball} -C /tmp/genome/ &

        # Determine if we need to re-write the ReadGroup for the Normal bam. Parabricks mutectcaller
        # requires the "ID:" fields in the tumor bam be different than those fields in the normal bam.
        # Currently MRD_alignment uses the LANE as the ID, so tumor and normal have overlapping IDs.
        # If that's the case we re-write the header to have a unique ID
        touch normal_ids_unsorted.txt
        touch tumor_ids_unsorted.txt
        for line in $(samtools view -H ~{normal_bam} | grep "^@RG") ; do
            if [[ "$line" =~ ^"ID:" ]]; then
                echo $line >> normal_ids_unsorted.txt
            fi
        done
        for line in $(samtools view -H ~{tumor_bam} | grep "^@RG") ; do
            if [[ "$line" =~ ^"ID:" ]]; then
                echo $line >> tumor_ids_unsorted.txt
            fi
        done
        if [ ~{normal_bam_defined} = true ]; then
            sort normal_ids_unsorted.txt > normal_ids_sorted.txt
            sort tumor_ids_unsorted.txt > tumor_ids_sorted.txt
            num_in_common=$(comm -12 tumor_ids_sorted.txt normal_ids_sorted.txt | wc -l)
            # If the tumor and normal bams have any ID: tags in common then re-write the normal_bam ID
            if (( num_in_common > 0 )); then
                echo "*** Tumor and Normal BAMs share ReadGroup ID(s) in common, need to re-write: $(comm -12 tumor_ids_sorted.txt normal_ids_sorted.txt | tr '\n' ' ') ***"
                # Save the first ReadGroup's entire header line to a file
                samtools view -H ~{normal_bam} | grep "^@RG" | head -1 > RG.txt

                # Replace (in-place) the first RG's ID: tag with ID:NEW_NORMAL_ to create a new, unique ReadGroup identifier
                sed -i -e 's/ID:/ID:NEW_NORMAL_/g' RG.txt
                echo "RE-WRITING NORMAL RG AS: $(cat RG.txt)"

                ### EXAMPLE COMMAND:
                ### samtools addreplacerg -r "@RG\tID:005\tLB:default\tPL:Illumina\tPU:HJLL2DMXY\tSM:HJLL2DMXY-00000303-MEX-22" -o HJLL2DMXY-00000303-MEX-22.updated.bam HJLL2DMXY-00000303-MEX-22.bam
                # Use the new RG header line (with ID:NEW_NORMAL_) to replace the RG header in the normal bam
                mv ~{normal_bam} ~{normal_bam}.original.bam
                samtools addreplacerg -@ 15 -r "`cat RG.txt`" -o ~{normal_bam} ~{normal_bam}.original.bam
                samtools index -@ 15 ~{normal_bam}
            else
                echo "ReadGroup IDs are unique between normal and tumor bams, no need to re-write RG IDs"
            fi
        else
            sort tumor_ids_unsorted.txt > tumor_ids_sorted.txt
            echo "Only tumor bam exists."
        fi

        # Index tumor bam if no index is present. If an index IS submitted,
        # touch the file to ensure it was modified more recently than the bam.
        if [ ~{tumor_bai_defined} = false ]; then
            samtools index -@ 15 ~{tumor_bam} &
        else
            touch ~{tumor_bai}
        fi
        # For the normal bam, the index may have been created as part of ReadGroup re-writing so
        # also check if it exists, not just if it was submitted to the task
        if [ ~{normal_bam_defined} = true ]; then
            if [ ~{normal_bai_defined} = false ] && ! [ -f ~{normal_bam}.bai  ]; then
                samtools index -@ 15 ~{normal_bam} &
            else
                # Here normal_bai may not be defined, but an index was created, so use the bam name to find it
                touch ~{normal_bam}.bai
            fi
        else
            echo "Only tumor bam exists."
        fi

        # Wait for any indexing to finish, as well as the genome untarring
        bg_pids=$(jobs -p)
        # Check if the background PID variable is empty
        if [ -z "$bg_pids" ]; then
            echo "No Jobs in the Background, no wait necessary..."
        else
            echo "Calling wait for background jobs with PIDs: $(echo $bg_pids)"
            wait $(echo $bg_pids)
        fi

        #### Set exome options - NOT currently available in parabricks mutect caller: "--disable-adaptive-pruning"
        ## The following are currently supported original mutectcaller options:
        ##     -pcr-indel-model <NONE, HOSTILE, AGGRESSIVE, CONSERVATIVE> (e.g. --mutectcaller-options="-pcr-indel-model HOSTILE")

        # For PB mutect, Sample names for both the normal and tumor samples need to be specified
        # Function to generate the tumor and normal mutect name args from each bam
        generate_name_arg() {
            local bam_path=$1
            # e.g. "--tumor-name" or "--normal-name"
            local arg_prefix=$2

            sample_names=()
            headers=$(samtools view -H $bam_path | grep "^@RG")
            for field in $headers; do
                if [[ "$field" =~ ^"SM:" ]]; then
                    sample_names+=(${field:3})
                fi
            done
            unique_sample_names=$(echo ${sample_names[@]} | tr ' ' '\n' | sort -u | tr '\n' ' ')
            sample_name_arg=""
            for name in $unique_sample_names; do
                sample_name_arg+="$arg_prefix $name "
            done

            # echo the result
            echo $sample_name_arg
        }

        if [ ~{normal_bam_defined} = true ]; then
            normal_name_arg=$(generate_name_arg "~{normal_bam}" "--normal-name")
        else
            normal_name_arg=$(generate_name_arg "" "--normal-name")
        fi

        tumor_name_arg=$(generate_name_arg "~{tumor_bam}" "--tumor-name")

        ## --mutect-low-memory is due to T4 GPUs having not quite enough GPU memory - should be able to remove with T5 GPUs
        ## --x3 just shows underlying call to mutectcaller with full args
        ## --run-partition splits up the available GPUs to run the calling in parallel for performance
        ref_fasta=`find /tmp/genome/ -name '*.fa' | head -n 1`
        
        if [ ~{normal_bam_defined} = true ]; then
            pbrun mutectcaller \
                --ref $ref_fasta \
                --in-tumor-bam ~{tumor_bam} \
                --in-normal-bam ~{normal_bam} \
                $tumor_name_arg \
                $normal_name_arg \
                --out-vcf somatic.vcf \
                --x3 --mutect-low-memory --run-partition
        else
            pbrun mutectcaller \
                --ref $ref_fasta \
                --in-tumor-bam ~{tumor_bam} \
                $tumor_name_arg \
                --out-vcf somatic.vcf \
                --x3 --mutect-low-memory --run-partition
        fi

        # bgzip the results, tabix index them, and rename the stats file to what is expected (.vcf.gz.stats)
        bgzip -l 9 somatic.vcf
        tabix -p vcf somatic.vcf.gz
        mv somatic.vcf.stats somatic.vcf.gz.stats
    >>>
    parameter_meta {
        genome_tarball: "stream"
    }
    runtime {
        # original_parabricks_docker_image: "nvcr.io/nvidia/clara/clara-parabricks:4.2.0-1"
        docker: "dx://myelin_mrd_images:/images/clara-parabricks-4.2.0-1-with-tabix.tar.gz"
        dx_instance_type: "mem2_ssd1_gpu_x48"
    }
    output {
        File somatic_vcf = "somatic.vcf.gz"
        File somatic_vcf_idx = "somatic.vcf.gz.tbi"
        File stats = "somatic.vcf.gz.stats"
    }
}
