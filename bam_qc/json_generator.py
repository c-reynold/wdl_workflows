import json
import pandas as pd
import dxpy as dx
import os

project_id = "project-J2zPy600ZV7y1F26fbzv4vk9"
output_dir = "bam_qc/input_jsons"
sample_sheet = "bam_qc/bam_qc_sample_sheet.csv"

number_cpu=8

def dxlink(file_id: str) -> dict:
    return {"$dnanexus_link": {"project": project_id, "id": file_id}}


df=pd.read_csv(sample_sheet)



for _,row in df.iterrows():
    bam_file=row["bam_id"]
    bai_file=row["bai_id"]

    bam_dx_name = dx.describe(bam_file)["name"].split(".")
    bai_dx_name = dx.describe(bai_file)["name"].split(".")

    if bam_dx_name[0] == bai_dx_name[0]:
        print(f"bam file {bam_dx_name[0]} (file-id = {bam_file}) and bai index file {bai_dx_name[0]} (file-id = {bai_dx_name} file names share the same root. Continuing...")
    else:
        print(f"bam file {bam_dx_name[0]} (file-id = {bam_file}) and bai index file {bai_dx_name[0]} (file-id = {bai_dx_name} file names DO NOT share the same root.")
        raise ValueError('A very specific bad thing happened.')


    


    sample_id = bam_dx_name[0]

    dict={

        "stage-common.bam": {
            "$dnanexus_link": {
            "id": bam_file,
            "project": project_id
            }
        },
        "stage-common.bam_idx": {
            "$dnanexus_link": {
            "id": bai_file,
            "project": project_id
            }
        },
        "stage-common.cpu": number_cpu

    }

    out_path = os.path.join(output_dir, f"{sample_id}.inputs.json")
    with open(out_path, "w") as f:
        json.dump(dict, f, indent=2)

print("Job completed. Thank you for using Christopher's amazing .json sample sheet generator!")