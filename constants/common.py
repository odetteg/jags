# All imports
import os
from pathlib import Path
import pandas as pd
import re
from dataclasses import dataclass, field

# All scripts
base_dir = Path(__file__).resolve().parent.parent
results_dir = os.path.join(base_dir, "results")
data_dir = os.path.join(results_dir, "data")
fastqc_dir=os.path.join(results_dir, "fastqc")
cmds_dir=os.path.join(base_dir, "temp")
download_script = os.path.join(base_dir, "workflows", "scripts", "download.sh")
# Assuming that all sample names are stored in a txt or csv file
samples_df = pd.read_csv(os.path.join(base_dir, "links.txt"), header=None)
sample_names = []
read = ["R1", "R2"]
ref_genome = "will/obtain/later"
def is_url(link):
    url_pattern=re.compile(r'^(https?|ftp)://[^\s/$.?#].[^\s]*$', re.IGNORECASE)
    return bool(url_pattern.match(link))
def getInputFiles(exts=["fastq","fastq.gz"]):
    list_of_files = []
    for file in samples_df.iloc[:,0].to_list():
        if any(file.endswith(ext) for ext in exts):
            if "_1" in file:
                file = file.replace("_1", "_R1").split("/")[-1]
            elif "_2" in file:
                file = file.replace("_2", "_R2").split("/")[-1]
            list_of_files.append(os.path.join(data_dir, file))
        else:
            list_of_files.append(os.path.join(data_dir, f"{file}_R1.fastq"))
            list_of_files.append(os.path.join(data_dir, f"{file}_R2.fastq"))
    return list_of_files
for id in samples_df.iloc[:, 0].to_list():
    if is_url(id) is True:
        sample_id = id.split("/")[-1]
        sample_names.append(sample_id)
    else:
        sample_names.append(id)
file_exts = set()

for file in getInputFiles():
    r1_ext = file.split("/")[-1]
    file_exts.add(r1_ext)
    r2_ext = r1_ext.replace("_R1", "_R2")
    file_exts.add(r2_ext)
exts_=list(file_exts)

accession_nums = set()               
for id in exts_:
    samples = id.split("_")[-2]
    accession_nums.add(samples)
ids_ = list(accession_nums)

def paired_end():
    paired_reads = []
    for read in getInputFiles():
        r1_ext = read.split("/")[-1]
        if "_R1" in r1_ext:
            r1 = os.path.join(data_dir, r1_ext)
            r2_ext = r1_ext.replace("_R1", "_R2")
            r2 = os.path.join(data_dir, r2_ext)
            if os.path.exists(r1) and os.path.exists(r2):
                paired_reads.append((r1, r2))
            else:
                print(f"You might be missing a pair for {read}")
    return paired_reads

def bwa_cmds(ref: str):
    map_file = os.path.join(cmds_dir, "map.txt")
    with open(map_file, 'w') as f:
        for r1, r2 in paired_end():
            out_=os.path.join(results_dir, "bwa", f"{os.path.basename(r1).split("_R1")[-2]}_aligned.sam" )
            bwa_cmds_ = f"bwa mem {ref} {r1} {r2} {out_}\n"
            f.write(bwa_cmds_)
    return map_file

def fastp_cmds():
    fastp_cmds_file = os.path.join(base_dir, "temp", "fastp_cmds.txt")
    with open(fastp_cmds_file, 'w') as f:
        for r1, r2 in paired_end():
                r1_out = os.path.join(base_dir, "results", "fastp", os.path.basename(r1))
                r2_out = os.path.join(base_dir, "results", "fastp", os.path.basename(r2))
                json=os.path.join(base_dir, "results", "fastp", "fastp.json")
                html=os.path.join(base_dir, "results", "fastp", "fastp.html")
                fastp_cmd = f"fastp -i {r1} -o {r1_out} -I {r2} -O {r2_out} -j {json} -h {html}\n"
                f.write(fastp_cmd)
    return fastp_cmds_file