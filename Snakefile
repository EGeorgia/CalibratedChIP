"""
Snakemake workflow to process Calibrated ChIP-seq data.
Author: Emily Georgiades
Hughes Group, University of Oxford
March 2020
"""

configfile: "config.yaml"
ruleorder:
    uniqAlignGen > uniqAlignSpi > indexing_genome > indexing_spikein > counting_gen > counting_spi

#IDS, = glob_wildcards("data/outputBams/{id}.bam")

rule all: # contains all final outputs as inputs
    input:
        expand("data/outputFiles/{sample}_readcount.txt",sample=config["samples"]),
        expand("data/outputFiles/{sample}_{spikein}_readcount.txt",sample=config["samples"],spikein=config["spikein"]),
        expand("data/outputFiles/{sample}_{genome}_readcount.txt",sample=config["samples"],genome=config["genome"])
rule bamming:
    input:
        READ1="data/samples/{sample}_READ1.fastq.gz",
        READ2="data/samples/{sample}_READ2.fastq.gz"
    params:
        bt2=config["bt2Prefix"]
    output:
        temp("data/outputBams/{sample}_UniqMapped.bam")
    conda:
        "envs/mapping.yaml"
    shell:
        "bowtie2 -x data/bowtie2/{params.bt2} -1 {input.READ1} -2 {input.READ2} -p 56 --no-mixed --no-discordant |"
        "grep -v XS: - |"
        "samtools view -b -h -S -F4 > {output}"

rule sorting:
    input:
        "data/outputBams/{sample}_UniqMapped.bam"
    output:
        temp("data/outputBams/{sample}_UniqMapped_sorted.bam")
    conda:
        "envs/mapping.yaml"
    shell:
        "sambamba sort  -t 56 -m 60G -o {output} {input}"

rule remove_duplicates:
    input:
        "data/outputBams/{sample}_UniqMapped_sorted.bam"
    output:
        "data/outputBams/{sample}_UniqMapped_sorted_rmdup.bam"
    conda:
        "envs/mapping.yaml"
    shell:
        "sambamba markdup -r -t 56 {input} {output}"

rule uniqAlignGen:
    input:
        "data/outputBams/{sample}_UniqMapped_sorted_rmdup.bam"
    params:
        genome=config["genome"],
        spikein=config["spikein"]
    output:
        gen="data/outputBams/{sample}_UniqMapped_sorted_rmdup_{genome}.bam"
    shell:
        "samtools view -h {input} |"
        "grep -v {params.spikein} | sed s/{params.genome}\_chr/chr/g |"
        "samtools view -bhS - > {output.gen}"

rule uniqAlignSpi:
    input:
        "data/outputBams/{sample}_UniqMapped_sorted_rmdup.bam"
    params:
        genome=config["genome"],
        spikein=config["spikein"]
    output:
        spi="data/outputBams/{sample}_UniqMapped_sorted_rmdup_{spikein}.bam"
    shell:
        "samtools view -h {input} |"
        "grep -v {params.genome} | sed s/{params.spikein}\_chr/chr/g |"
        "samtools view -bhS - > {output.spi}"

rule indexing_genome:
    input:
        "data/outputBams/{sample}_UniqMapped_sorted_rmdup_{genome}.bam"
    params:
        genome=config["genome"],
        spikein=config["spikein"]
    output:
        temp("data/outputBams/{sample}_UniqMapped_sorted_rmdup_{genome}.bam.bai")
    shell:
        "sambamba index -t 56 {input}"

rule indexing_spikein:
    input:
        inspikein="data/outputBams/{sample}_UniqMapped_sorted_rmdup_{spikein}.bam"
    params:
        genome=config["genome"],
        spikein=config["spikein"]
    output:
        temp("data/outputBams/{sample}_UniqMapped_sorted_rmdup_{spikein}.bam.bai")
    shell:
        "sambamba index -t 56 {input}"

rule counting:
    input:
        "data/outputBams/{sample}_UniqMapped_sorted_rmdup.bam"
    output:
        "data/outputFiles/{sample}_readcount.txt"
    shell:
        "sambamba view -c -t 50 {input} >> {output}"

rule counting_gen:
    input:
        "data/outputBams/{sample}_UniqMapped_sorted_rmdup_{genome}.bam"
    params:
        genome=config["genome"],
        spikein=config["spikein"]
    output:
        "data/outputFiles/{sample}_{genome}_readcount.txt"
    shell:
        "sambamba view -c -t 50 {input} >> {output}"

rule counting_spi:
    input:
        "data/outputBams/{sample}_UniqMapped_sorted_rmdup_{spikein}.bam"
    params:
        genome=config["genome"],
        spikein=config["spikein"]
    output:
        "data/outputFiles/{sample}_{spikein}_readcount.txt"
    shell:
        "sambamba view -c -t 50 {input} >> {output}"
