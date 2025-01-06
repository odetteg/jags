#!/bin/bash

declare -A snakemake_input=( )
declare -A snakemake_output=( [0]="/Users/sg/Documents/micromake/results/data/ERR246975_R2.fastq" )
declare -A snakemake_params=( )
declare -A snakemake_wildcards=( [0]="ERR246975_R2.fastq" [sample]="ERR246975_R2.fastq" )
declare -A snakemake_resources=( [0]="1" [_cores]="1" [1]="1" [_nodes]="1" [2]="/var/folders/fk/ftbf2rrs5lv2cvwbp68vdv5r0000gn/T" [tmpdir]="/var/folders/fk/ftbf2rrs5lv2cvwbp68vdv5r0000gn/T" )
declare -A snakemake_log=( )
declare -A snakemake_config=( [snakefile]="snakefile" [cores]="8" )
declare -A snakemake=( [_params_types]="{}" [threads]="1" [rule]="laod_data" [bench_iteration]="None" [scriptdir]="/Users/sg/Documents/micromake/workflows/scripts" )
set -u
set -o pipefail
accessions_file="links.txt"
data_dir="results/data"
mkdir -p $data_dir

# Downloading files when given a list of accession numbers
accessions=()
while IFS= read -r accession_no || [[ -n "$accession_no" ]]; do
    accessions+=("$accession_no")
done < "$accessions_file"

for id in ${accessions[@]}; do
    if [[ "$id" == "SRR"* ||  "$id" == "ERR"* ]]; then
        fasterq-dump $id --split-files -O $data_dir
    elif [[ "$id" == "ftp://"* ]]; then 
        out_dir=$(basename "$id")
        wget -nc -P "$data_dir" "$id"
    else
        echo "Invalid accession number"
        exit 1
    fi
done
 
for fq in results/data/*.fastq*; do
    if [[ "$fq" == *"_1.fastq" && "$fq" != *"_R1.fastq" ]]; then
        renamed_r1="${fq/_1.fastq/_R1.fastq}"
        mv "$fq" "$renamed_r1"
    elif [[ "$fq" == *"_1.fastq.gz" && "$fq" != *"_R1.fastq.gz" ]]; then
        renamed_r1="${fq/_1.fastq.gz/_R1.fastq.gz}"
        mv "$fq" "$renamed_r1"
    elif [[ "$fq" == *"_2.fastq" && "$fq" != *"_R2.fastq" ]]; then
        renamed_r2="${fq/_2.fastq/_R2.fastq}"
        mv "$fq" "$renamed_r2"
    elif [[ "$fq" == *"_2.fastq.gz" && "$fq" != *"_R2.fastq.gz" ]]; then
        renamed_r2="${fq/_2.fastq.gz/_R2.fastq.gz}"
        mv "$fq" "$renamed_r2"
    fi
done
