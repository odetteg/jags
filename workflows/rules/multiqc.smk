import sys
from pathlib import Path
from constants.common import *

rule multiqc:
    input:
        expand("results/fastqc/{sample}_{read}_fastqc.{format}", read=read, format=["html", "zip"], sample=ids_),
        fastqc_dir,
        expand(str(base_dir) + "/results/fastp/{sample}", sample=exts_),
        expand(str(base_dir) + "/results/fastp/{rpt}", rpt=["fastp.json", "fastp.html"]),
        expand(
            str(results_dir) + "/bwa/fastp/bam/fastp_trimmed_{sample}.sorted.bam", sample=ids_
        ),
        expand(
            str(results_dir) + "/variants/bcftools/fastp_{sample}.filtered.vcf", sample=ids_
        ),
        expand(
            str(results_dir) + "/variants/bcftools/fastp_{sample}.{type}.stats.txt", sample=ids_, type=["raw", "filtered"]
        )
    output:
        directory(multiqc_dir),
        report(str(results_dir) + "/multiqc/multiqc_report.html", caption="report/multiqc.rst", category="Quality control"),

    log:
        "logs/multiqc.log",
    params:
        outdir = multiqc_dir,
    shell:
        """
        multiqc --force --outdir {params.outdir} {input} >> {log} 2>&1
        """