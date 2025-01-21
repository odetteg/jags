
configfile: "config/config.yaml"
include: "workflows/rules/qc.smk"
include:  "workflows/rules/map.smk"
include:  "workflows/rules/call.smk"
include:  "workflows/rules/multiqc.smk"
rule all:
    input:
        multiqc_dir,
        report(
            str(results_dir) + "/multiqc/multiqc_report.html", caption="report/multiqc.rst", category="Quality control"),
