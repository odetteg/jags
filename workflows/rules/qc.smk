import sys
from pathlib import Path
from constants.common import *

rule laod_data:
    output:
        *getInputFiles(),
        get_ref_genome()
    log:
        str(base_dir) + "/logs/download.log"
    script:
        download_script
rule fastqc:
    input:
        *getInputFiles()
    output:
        fqs=expand("results/fastqc/{sample}_{read}_fastqc.{format}", read=read, format=["html", "zip"], sample=ids_),
        fqc_out=directory(fastqc_dir)
    message:
        "Running FASTQC with on {input}"
    log:
        str(base_dir) + "/logs/fastqc.log"
    shell:
        """
        mkdir -p {output.fqc_out}
        fastqc -o {output.fqc_out} -f fastq {input} >> {log} 2>&1
        """


rule get_fastp_cmds:
    output:
        temp(cmds_dir + "/fastp.txt")
    run:
        fastp_cmds()
rule fastp:
    input:
        *getInputFiles(),
        fastp_cmds_= os.path.join(cmds_dir, "fastp.txt")
    output:
        expand(str(base_dir) + "/results/fastp/{sample}", sample=exts_),
        expand(str(base_dir) + "/results/fastp/{rpt}", rpt=["fastp.json", "fastp.html"])
    message:
        "Running FASTP on {input}"
    log:
        str(base_dir) + "/logs/bwa.log"
    shell:
        """
        while read -r cmd;do
        $cmd
        done < {input.fastp_cmds_} >> {log} 2>&1
        """
