import os
import pandas as pd
import json
import dxpy as dx
import sys
from pathlib import Path

#RD_Analysis project: project-GPgxkj8070pyZXjk58KbQffV
#liquid_hrd project:project-J2zPy600ZV7y1F26fbzv4vk9

#These should be the project in which the fastq files are located.
fastq_project_id = "project-GPgxkj8070pyZXjk58KbQffV"

ref_project_id = "project-J2zPy600ZV7y1F26fbzv4vk9"
reference_genome = "file-GpkQV9Q07fG36gzq1v1K8FjB"

sample_sheet = "parabricks_fq2bam/sample_sheets/3-16-26_sample_sheet_trim2.csv"
output_dir = "parabricks_fq2bam/fastq2bam_sample_jsons"

#Whether or not to use the File ID to look up the DNA nexus filename to validate that the filenames for 2 paried reads differ only by "R1" and "R2"
#It checks for the number of differences (only 1 difference allowed) and confirms that the difference is a "1" or a "2", no other differences are allowed.
#Set this to false if the DNA nexus filenames are not consistent between reads. This should be a rare and extenuating circumstance if this is not the case.
#Fastq filenames really should not be different between read 1 and read 2. If they are, then whatever is generating fastq files needs to be fixed or changed.
#This check is slow with many fastq
filename_check = True

Path(output_dir).rmdir
Path(output_dir).mkdir(parents=True, exist_ok=True)





os.makedirs(output_dir, exist_ok=True)

#CUSTOM FUNCTION DEFINITIONS
def dxlink(file_id: str, project_id: str) -> dict:
    return {"$dnanexus_link": {"project": project_id, "id": file_id}}

def difference(str1, str2):
    """Compares two strings and returns a list of character differences.
    Output is [character position, [char_from_str1, char_from_str2]]
    
    """
    list1=list(str1)
    list2=list(str2)

    largest_string =  max(len(list1), len(list2))

    diffs = []
    for i in range(largest_string):
        c1 = list1[i] if i < len(list1) else None
        c2 = list2[i] if i < len(list2) else None

        if c1 != c2:
            diffs.append([i, [c1, c2]])

    return diffs

samples = pd.read_csv(sample_sheet)

#check for duplicated files and/or file ids

concatonated_samples = samples["sample_id"].astype(str) + "_" + samples["lane"].astype(str) + "_" + samples["read"].astype(str)

duplicate_samples_list = samples[concatonated_samples.duplicated(keep=False)]
duplicate_file_id_list = samples[samples["file_id"].duplicated(keep=False)]
empty_file_id_list = samples[samples["file_id"].isna()]


#no duplicate samples can exist for the pd.pivot that occurs later. 
if len(duplicate_samples_list) == 0:
    print("No duplicate samples found. Proceeding to check File IDs...")
else:
    print("Duplicate samples with identical read and lane information found. Please remove one or both of these duplicates from the sample sheet:")
    print(duplicate_samples_list.to_string())
    sys.exit(1)

if len(duplicate_file_id_list) == 0:
    print("No duplicate DNA Nexus File IDs identified. Continuing...")
else:
    print("Duplicate DNA Nexus File IDs identified on the sample sheet. Please remove one or both of these duplicates from the sample sheet:")
    print(duplicate_file_id_list.to_string())
    sys.exit(1)

if len(empty_file_id_list) == 0:
    print("All samples contain a DNA Nexus File ID. Continuing...")
else:
    print("The following samples are missing a DNA Nexus File ID:")
    print(empty_file_id_list.to_string())
    sys.exit(1)




#The sample dataframe is pivoted to make it easier to assign R1 and R2 groups correctly. This puts R1 and R2 in the same row of the dataframe.
#In order for this pivot to work, there cannot be any duplicate samples with the same read and lane information.
#In the case that flowcells need to be combined, this means that they should be combined at the .bam level, not the alignment level.
#In a future version of this code, it would be quite easy to add an additional column to the sample sheet for "flowcell", which would allow the combining of samples from flowcells
#prior to alignment instead of after.
#The parabricks FQ2BAM program requires read 1 and read 2 from the same lanes to be present in the input args.
#This is the reason behind all of this checking. Perhaps in the future, this validation/verification can be done in an initial step of the fq2bam.wdl
#For now, we are using file names to get read and  
pivoted_samples = samples.pivot(
    index=["sample_id", "lane", "trim_length"],
    columns="read",
    values="file_id"
    ).reset_index()

pivoted_samples.rename(columns={1:"R1", 2:"R2"}, inplace=True)


#This block checks to make sure that there are not any R1 or R2 reads missing. I.e., a sample that has read 1 for a lane, but not a read 2, or vice versa
missing_mask = pivoted_samples[["R1", "R2"]].isna().any(axis=1)

if missing_mask.any():
    bad_rows = pivoted_samples.loc[
        missing_mask,
        ["sample_id", "lane", "trim_length", "R1", "R2"]
    ]

    print("The following sample/lane combinations are missing R1 or R2:")
    print(bad_rows.to_string(index=False))

    raise ValueError(f"{len(bad_rows)} sample/lane rows are missing R1 or R2 after pivot")

print("All samples have an R1 and R2 associated with them. Continuing...")

#This section is for double checking DNA nexus filenames by directly looking them up with the dx.describe command
if filename_check:

    #This variable ticks up with errors and if there are 1 or more errors it prints the number of errors before exiting
    filename_read_check = 0

    #For each row in the pivoted frame, we will now check that:
    #The length of the filenames in DNA nexus are exactly the same for R1/R2 pairs

    for _, row in pivoted_samples.iterrows():
        
        #pull file IDs from the dataframe from the csv file
        file_id_read_1=(row["R1"])
        file_id_read_2=(row["R2"])

        #Get the filename on DNA nexus from the fileID
        filename_r1 = dx.describe(file_id_read_1)["name"] #type: ignore
        filename_r2 = dx.describe(file_id_read_2)["name"] #type: ignore

        print(f"Using {file_id_read_1} and {file_id_read_2} checking {filename_r1} and {filename_r2}")

        diffs= difference(filename_r1,filename_r2)

        if len(diffs) != 1:
            print(f"{filename_r1} and {filename_r2} differ at {len(diffs)} positions (expected exactly 1).")
            filename_read_check = filename_read_check + 1
            continue

        pos, (r1, r2) = diffs[0]

        if {r1,r2} != {"1","2"}:
            print(f"The difference between {filename_r1} and {filename_r2} is not expected. Expecting difference only in the read number in the DNA Nexus File names")
            print(f"Instead, {filename_r1} has a '{r1}' while {filename_r2} has a '{r2}' at at character {pos}")
            filename_read_check = filename_read_check + 1
            continue

    if filename_read_check > 0:
        print("There were one or more errors when validating DNA Nexus file names")
        print(f"Total errors: {filename_read_check}")
        print("Exiting...")
        sys.exit(1)

grouped_frames = pivoted_samples.groupby("sample_id")




#This block is what actually builds the sample sheets
for sample_id, g in grouped_frames:

    if g["trim_length"].nunique() != 1:
        raise ValueError(f"Sample {sample_id} has multiple trim_length values: {g['trim_length'].unique()}")

    # Sort lanes so R1/R2 arrays align by lane order (L001, L002, ...)
    g = g.sort_values("lane")

    r1_list = [dxlink(fid,fastq_project_id) for fid in g["R1"].tolist()]
    r2_list = [dxlink(fid,fastq_project_id) for fid in g["R2"].tolist()]

    # If trim_length is per-sample, take the first (or validate uniqueness upstream)
    umi = int(g["trim_length"].iloc[0])

    inputs = {
        "stage-common.r1_fastqs": r1_list,
        "stage-common.r2_fastqs": r2_list,
        "stage-common.trim_length": umi,
        "stage-common.reference_genome": dxlink(reference_genome, ref_project_id)
    }

    print(f"All checks for {sample_id} were successful. Generating input.json... ")

    out_path = os.path.join(output_dir, f"{sample_id}.inputs.json")
    with open(out_path, "w") as f:
        json.dump(inputs, f, indent=2)

    print(f"{sample_id} input .json file created at {output_dir}/{sample_id}.inputs.json")



   