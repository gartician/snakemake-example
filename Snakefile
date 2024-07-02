import os
import glob
import pandas as pd

"""
Teaching example to download CUT&Tag datasets, rename files, and align to a reference genome with bowtie2
"""

configfile: 'config/config.yaml'

# list all SRA accession IDs to download
sra_examples = [
    "SRR8383512",
    "SRR8383513",
    "SRR8383514",
    "SRR8383515",
    "SRR8383516",
    "SRR8383517",
    "SRR8435051",
    "SRR8435052" 
]

# list corresponding sample names
sample_examples = [
    "K562_1_H3K4me1",
    "K562_2_H3K4me1",
    "K562_1_H3K4me2",
    "K562_2_H3K4me2",
    "K562_1_H3K4me3",
    "K562_2_H3K4me3",
    "K562_1_IgG",
    "K562_2_IgG",
]

# the above lists can be read through a configuration table, yaml file, or more!

# link SRA accession IDs to a sample name in two different ways
sra_samples_dict = {}
for sra, sample in zip(sra_examples, sample_examples):
    sra_samples_dict[sample] = sra

print("\nExample SRA-sample dictionary")
for i, j in sra_samples_dict.items():
    print(f"{i}: {j}")

sra_samples_df = pd.DataFrame({
    "sra": sra_examples,
    "sample": sample_examples
})

print("\nExample SRA-sample dataframe")
print(sra_samples_df)

rule all:
    input:
        # download FASTQs
        expand("data/raw/{sra}_1.fastq.gz", sra = sra_examples),
        expand("data/raw/{sra}_2.fastq.gz", sra = sra_examples),

        # rename files (the following examples are equivalent ways to express output files)
        expand("data/rename/{sample}_R1.fastq.gz",sample = sra_samples_dict.keys()),
        expand("data/rename/{sample}_R2.fastq.gz",sample = sra_samples_dict.keys()),
        # expand("data/rename/{sample}_{read}.fastq.gz", sample = sra_samples_dict.keys(), read = ["R1", "R2"]),
        # expand("data/rename/{sample}_{read}.fastq.gz", sample = sra_samples_df["sample"], read = ["R1", "R2"])
        
        # align reads to genome
        expand("data/bowtie2/{sample}.bam", sample = sample_examples)

rule download:
    output:
        "data/raw/{sra}_1.fastq.gz",
        "data/raw/{sra}_2.fastq.gz"
    conda:
        "workflow/envs/sra-tools.yaml"
    shell:
        "fasterq-dump -O data/raw -p --split-files {wildcards.sra}; "
        "gzip -f data/raw/{wildcards.sra}_1.fastq; "
        "gzip -f data/raw/{wildcards.sra}_2.fastq "

def get_sra(wildcards):
    sra_id = sra_samples_dict[wildcards.sample]
    R1 = f"data/raw/{sra_id}_1.fastq.gz"
    R2 = f"data/raw/{sra_id}_2.fastq.gz"
    return(R1, R2)

rule rename:
    input:
        get_sra
    output:
        R1 = "data/rename/{sample}_R1.fastq.gz",
        R2 = "data/rename/{sample}_R2.fastq.gz"
    shell:
        "cp {input[0]} {output.R1}; "
        "cp {input[1]} {output.R2} "

rule subset:
    input:
        R1 = "data/rename/{sample}_R1.fastq.gz",
        R2 = "data/rename/{sample}_R2.fastq.gz"
    output:
        R1 = "data/subset/{sample}_R1.fastq.gz",
        R2 = "data/subset/{sample}_R2.fastq.gz"
    params:
        rows = 1_000_000
    shell:
        "zcat {input.R1} | sed -n '1,{params.rows}p' | gzip -f > {output[0]}; "
        "zcat {input.R2} | sed -n '1,{params.rows}p' | gzip -f > {output[1]} "
# something wrong with the head function so I'm using sed for now

rule bowtie2:
    input:
        R1 = "data/subset/{sample}_R1.fastq.gz",
        R2 = "data/subset/{sample}_R2.fastq.gz"
    output:
        "data/bowtie2/{sample}.bam"
    params:
        genome = config["genome"]
    threads: 2
    conda:
        "workflow/envs/bowtie2.yaml"
    shell:
        "bowtie2 "
        "--local "
        "--very-sensitive-local "
        "--no-unal "
        "--no-mixed "
        "--threads {threads} "
        "--no-discordant "
        "--phred33 "
        "-I 10 "
        "-X 700 "
        "-x {params.genome} "
        "-1 {input.R1} -2 {input.R2} | samtools view -Sbh - > {output}"