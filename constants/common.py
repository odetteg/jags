# All imports
import os
from pathlib import Path
import pandas as pd
import re
from dataclasses import dataclass, field
from typing import List, Tuple
import urllib.request
import yaml

# All scripts
base_dir = Path(__file__).resolve().parent.parent
results_dir = os.path.join(base_dir, "results")
data_dir = os.path.join(results_dir, "data")
fastqc_dir=os.path.join(results_dir, "fastqc")
cmds_dir=os.path.join(base_dir, "temp")
ref_genome_dir = os.path.join(results_dir, "ref")
download_script = os.path.join(base_dir, "workflows", "scripts", "download.sh")
# Assuming that all sample names are stored in a txt or csv file
samples_df = pd.read_csv(os.path.join(base_dir, "links.txt"), header=None)
sample_names = []
read = ["R1", "R2"]
ref_exts = ["fa", "fna.gz", "fasta", "fna", "fasta.gz"]  
with open("config/config.yaml", 'r') as f:
    config = yaml.safe_load(f)
bwa_mem_params = config["bwa"]["bwa_mem"]["extra_params"]
def is_url(link):
    url_pattern=re.compile(r'^(https?|ftp)://[^\s/$.?#].[^\s]*$', re.IGNORECASE)
    return bool(url_pattern.match(link))
def is_ref_file(file: str, ref_exts = ref_exts) -> bool:
    for ext in ref_exts:
        if file.endswith(f".{ext}"):
            return True
    return False
#Major issue to resolve here
def getInputFiles(exts=["fastq","fastq.gz"]):
    list_of_files = []
    for file in samples_df.iloc[:,0].to_list():
        if any(file.endswith(ext) for ext in exts):
            if "_1" in file:
                file = file.replace("_1", "_R1").split("/")[-1]
            elif "_2" in file:
                file = file.replace("_2", "_R2").split("/")[-1]
            list_of_files.append(os.path.join(data_dir, file))
        elif not is_ref_file(file):
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
ids_unfiltered = list(accession_nums)

ids_ = [id for id in ids_unfiltered if not id.endswith(tuple(ref_exts))]
# def download_ref_genome(exts: list[str] = ref_exts):
#     for link in samples_df.iloc[:,0].to_list():
#         if any(link.endswith(ext) for ext in exts):
#             ref_file_name = os.path.join(ref_genome_dir, os.path.basename(link))
#             try:
#                 urllib.request.urlretrieve(link, ref_file_name)
#                 print(f"Ref file {os.path.basename(link)} succesfullfy downloaded and stored in the {ref_genome_dir}")
#             except Exception as e:
#                 print(f"Enountered a problem while downloading the {os.path.basename(link)}. {e}")
                        
def get_ref_genome(exts: str = ref_exts) -> str | None:
    ref_list = []
    for link in samples_df.iloc[:,0].to_list():
        if any(link.endswith(ext) for ext in exts):
            _genome_name = os.path.join(ref_genome_dir, os.path.basename(link))
            ref_list.append(_genome_name)
    if len(ref_list) >1:
        ref_genome = sorted(ref_list)[0]
    elif len(ref_list) == 1:
        ref_genome = ref_list[0]
    else:
        ref_genome = None
    return ref_genome  
def get_ref_genome_name(link_url: str = get_ref_genome(), ref_exts=ref_exts):
    for ext in ref_exts:
        base_name = os.path.basename(link_url)
        if ext in link_url:
            return base_name.split(f".{ext}")[0], ext 
ref_name = get_ref_genome_name()[-2]
ext_in_ref = get_ref_genome_name()[-1]
@dataclass
class PairedFiles:
    data_dir: str = data_dir

    def paired_end(self) -> List[Tuple[str, str]]:
        paired_reads = []
        for read in getInputFiles():
            r1_ext = read.split("/")[-1]
            if "_R1" in r1_ext:
                r1 = os.path.join(self.data_dir, r1_ext)
                r2_ext = r1_ext.replace("_R1", "_R2")
                r2 = os.path.join(self.data_dir, r2_ext)
                paired_reads.append((r1, r2))
        return paired_reads

def bwa_cmds(ref: str, paired: PairedFiles = PairedFiles()) -> tuple[str, str]:
    fastp_map_file = os.path.join(cmds_dir, "fp_map.sh")
    trimmomatic_map_file = os.path.join(cmds_dir, "tm_map.sh")
    with open(fastp_map_file, 'w') as fp, open(trimmomatic_map_file, 'w') as tp:
            fp.write(f"mkdir -p {results_dir}/bwa/fastp/sam\n")
            fp.write(f"mkdir -p {results_dir}/bwa/fastp/sam\n")
            for r1, r2 in paired.paired_end():
                fp_r1, fp_r2 = r1.replace("data", "fastp"), r2.replace("data", "fastp")
                tp_r1, tp_r2 = r1.replace("data", "fastp"), r2.replace("data", "fastp")
                sample_name = os.path.basename(r1).split("_R1")[0]
                fp_out_ = os.path.join(results_dir, "bwa", "fastp", "sam", f"fastp_trimmed_{sample_name}_aligned.sam")
                tm_out_ = os.path.join(results_dir, "bwa", "trimmomatic", "sam", f"matic_trimmed_{sample_name}_aligned.sam")
                fp_bwa_cmds_ = f"bwa mem {bwa_mem_params} {ref} {fp_r1} {fp_r2} > {fp_out_}\n"
                tm_bwa_cmds_ = f"bwa mem {bwa_mem_params} {ref} {tp_r1} {tp_r2} > {tm_out_}\n"
                fp.write(fp_bwa_cmds_)
                tp.write(tm_bwa_cmds_)
            return fastp_map_file, trimmomatic_map_file

def fastp_cmds(paired: PairedFiles = PairedFiles()) -> str:
    fastp_cmds_file = os.path.join(base_dir, "temp", "fastp.txt")
    with open(fastp_cmds_file, 'w') as f:
        for r1, r2 in paired.paired_end():
            r1_out = os.path.join(base_dir, "results", "fastp", os.path.basename(r1))
            r2_out = os.path.join(base_dir, "results", "fastp", os.path.basename(r2))
            json = os.path.join(base_dir, "results", "fastp", "fastp.json")
            html = os.path.join(base_dir, "results", "fastp", "fastp.html")
            fastp_cmd = f"fastp -i {r1} -o {r1_out} -I {r2} -O {r2_out} -j {json} -h {html}\n"
            f.write(fastp_cmd)
    return fastp_cmds_file