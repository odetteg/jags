import sys
from pathlib import Path
from constants.common import *

rule get_map_cmds:
    output:
        str(cmds_dir) + "/fp_map.sh"
    run:
        bwa_cmds(ref=get_ref_genome())
rule index_fa:
    input:
        ref_file = get_ref_genome()
    output:
        index_out_ = expand(
            str(ref_genome_dir) + "/{ref}.{ext_in_ref}.{ext}",
            ref=ref_name,
            ext=["amb", "ann", "bwt", "pac", "sa"],
            ext_in_ref=ext_in_ref
        ),
    shell:
        """
        bwa index {input.ref_file}
        """
rule map:
    input:
        ref_file = get_ref_genome(),
        fp_map_cmds =  str(cmds_dir) + "/fp_map.sh",
        fastp_filtered = expand(
            str(results_dir + "/fastp/{sample}"), sample=exts_
        ),

        index_files = expand(
            str(ref_genome_dir) + "/{ref}.{ext_in_ref}.{ext}",
            ref=ref_name,
            ext=["amb", "ann", "bwt", "pac", "sa"],
            ext_in_ref=ext_in_ref
        ),
    output:
        fastp_bam = expand(
            str(results_dir) + "/bwa/fastp/bam/fastp_trimmed_{sample}.sorted.bam", sample=ids_
        ),
        fastp_bai = expand(
            str(results_dir) + "/bwa/fastp/bam/fastp_trimmed_{sample}.sorted.bam.bai", sample=ids_
        )
    message:
        "Mapping reads..."
    shell:
        """
        set -x 
        chmod 775 {input.fp_map_cmds}
        {input.fp_map_cmds}
        for bam_fi in {output.fastp_bam}; do
        fp_sam="${{bam_fi/bam/sam}}"
        fp_sam="${{fp_sam%.sorted.bam}}_aligned.sam"
        unsorted_bam="${{bam_fi%.sorted.bam}}.bam"
        samtools view -Sb $fp_sam -o $unsorted_bam
        done
        for bam_file in {output.fastp_bam};do
        bam_in="${{bam_file%.sorted.bam}}.bam"
        samtools sort $bam_in -o $bam_file && samtools index $bam_file
        done
        """
    