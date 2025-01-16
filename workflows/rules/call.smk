import sys
from pathlib import Path
from constants.common import *

rule call:
    input:
        fastp_bam = expand(
            str(results_dir) + "/bwa/fastp/bam/fastp_trimmed_{sample}.sorted.bam", sample=ids_
        ),
        fastp_bai = expand(
            str(results_dir) + "/bwa/fastp/bam/fastp_trimmed_{sample}.sorted.bam.bai", sample=ids_
        ),
        ref_file = get_ref_genome(),
    output:
        raw_vcfs= expand(
            str(results_dir) + "/variants/bcftools/fastp_{sample}.raw.vcf", sample=ids_
        )
    message:
        "Running bcftools..."
    shell:
        """
        set -x
        for bam_file in {input.fastp_bam}; do
        base_name=$(basename $bam_file | sed 's/^fastp_trimmed_//' | sed 's/\\.sorted.bam$//')
        bai={results_dir}/bwa/fastp/bam/fastp_trimmed_${{base_name}}.sorted.bam.bai
        vcf_file={results_dir}/variants/bcftools/fastp_${{base_name}}.raw.vcf
        bcftools mpileup -Ov -f {input.ref_file} $bam_file | \
        bcftools call -mv -Ov -o $vcf_file
        done
        """
# rule filter_variants:
#     input:
#         raw_vcfs= expand(
#             str(results_dir) + "/variants/bcftools/fastp_{sample}.raw.vcf", sample=ids_
#         )
#     output:
#         filtered_vcfs= expand(
#             str(results_dir) + "/variants/bcftools/fastp_{sample}.filtered.vcf", sample=ids_
#         )
#     shell:
#         """
        
#         """