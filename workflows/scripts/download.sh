#!/bin/bash
set -u
set -o pipefail
accessions_file="links.txt"
data_dir="results/data"
ref_dir="results/ref"

# Downloading files when given a list of accession numbers
accessions=()
while IFS= read -r accession_no || [[ -n "$accession_no" ]]; do
    accessions+=("$accession_no")
done < "$accessions_file"

for id in ${accessions[@]}; do
    if [[ "$id" == "SRR"* ||  "$id" == "ERR"* ]]; then
        fasterq-dump $id --split-files -O $data_dir
    elif [[ "$id" == "ftp://"* ]]; then 
        wget -nc -P "$data_dir" "$id"
    elif [[ "$id" == *.fna.gz* || "$id" == *.fa.gz || "$id" == *.fasta* || "$id" == *fasta.gz ]]; then
        wget -nc -P "$ref_dir" "$id"
        ref_file=$(basename "$id")
        ref_path="$ref_dir/$ref_file"
        if file "$ref_path" | grep -q "BGZF"; then
        echo "file in the correct format"
        else 
        base_name="$ref_dir/$(basename "$id" .gz)"
        gunzip "$ref_path"
        bgzip "$base_name"
        fi
    else
        echo "Skipping download for $id"
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
