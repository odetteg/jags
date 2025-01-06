
configfile: "config/config.yaml"
include: "workflows/rules/qc.smk"
include:  "workflows/rules/map.smk"
rule all:
    input:
        expand("results/fastqc/{sample}_{read}_fastqc.{format}", read=read, format=["html", "zip"], sample=ids_),
        expand(str(base_dir) + "/results/fastp/{sample}", sample=exts_),
        expand(str(base_dir) + "/results/fastp/{rpt}", rpt=["fastp.json", "fastp.html"]),