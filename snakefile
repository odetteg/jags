
configfile: "config/config.yaml"
include: "workflows/rules/qc.smk"
include:  "workflows/rules/map.smk"
include:  "workflows/rules/call.smk"
rule all:
    input:
        expand(
            str(results_dir) + "/variants/bcftools/fastp_{sample}.raw.vcf", sample=ids_
        ),
        