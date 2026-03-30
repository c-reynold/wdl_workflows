import dxpy as dx
import pandas as pd
import subprocess

app="applet-J71Qx600ZV7ZqBK6fyx1520J"
project_id="project-J2zPy600ZV7y1F26fbzv4vk9"
app_name=dx.describe(app)["name"] # type: ignore

#These are specifics for 
disk_size=300 #This is the disk size in GB that will be requested for the instance on DNA nexus
cpu=4 




priority="high" #can be {low, normal, high}

sample_csv="samtools_downsample/launch/sample_sheets/samtools_downsample_.csv"
dx_parent_dir=f"{project_id}:/workflows/bam_utilities/samtools_downsample/results/downsample_to_1.3X_cfdna_bams"





df=pd.read_csv(sample_csv, index_col=None)

print(df)

for _,row in df.iterrows():

    sample_name=row["sample_name"]
    bam_id=row["bam_file_id"]
    bai_id=row["bai_file_id"]
    current_coverage=row["current_coverage"]
    desired_coverage=int(row["desired_coverage"])


    #raises an error if the desired coverage is larger than current coverage
    if desired_coverage > current_coverage:
        raise ValueError(f"Desired coverage of {desired_coverage} is larger than current coverage of {current_coverage}")

    dx_run_cmd=["dx", "run", app,
                   f"-isample_name={sample_name}", 
                   f"-ibam={bam_id}",
                   f"-ibai={bai_id}",
                   f"-icpu={cpu}",
                   f"-idisk_size={disk_size}",
                   f"-icurrent_coverage={current_coverage}",
                   f"-idesired_coverage={desired_coverage}",
                   "--name", f"{app_name}_{sample_name}_{desired_coverage}x", 
                   "--destination", f"{dx_parent_dir}/{sample_name}",
                   "--yes", "--priority", f"{priority}"
    ]

    print(dx_run_cmd)
    bash_run = ' '.join(str(x) for x in dx_run_cmd)
    print(f"this is the bash command:{bash_run}")


    job=subprocess.run(dx_run_cmd, capture_output=True, text=True, check=True)
    print(dx_run_cmd)

  






 

