#fastq_to_bam_paired.snakefile
#Anna Hoge
#Template made October 18th, 2019
#Ha Lab
#Fred Hutchinson Cancer Research Center

"""
#before running snakemake, do in tmux terminal:
ml snakemake/5.2.4-foss-2016b-Python-3.6.6
ml BWA/0.7.17-foss-2016b
ml SAMtools/1.9-foss-2016b
ml java/jdk1.8.0_31
ml GATK/4.1.0.0-foss-2016b-Python-3.6.6
ml picard/2.18.29-Java

#command to run snakemake (remove -np at end when done validating):
snakemake -s fastq_to_bam_paired.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 100 -np
"""

configfile: "config/samples.yaml"
configfile: "config/config.yaml"

rule all:
  input:
    expand("results/{samples}/{samples}_unsorted.bam", samples=config["samples"]),
    expand("results/{samples}/{samples}_sorted.bam", samples=config["samples"]),
    expand("results/{samples}/{samples}_dups_marked.bam", samples=config["samples"]),
    expand("results/{samples}/{samples}_marked_dup_metrics.txt", samples=config["samples"]),
    expand("results/{samples}/{samples}_recalibration_data.table", samples=config["samples"]),
    expand("results/{samples}/{samples}_recalibrated.bam", samples=config["samples"]),
    expand("results/{samples}/{samples}_recalibrated.bam.bai", samples=config["samples"])

rule map_to_reference:
    input:
        fastq1=lambda wildcards: config["samples"][wildcards.samples][0],
        fastq2=lambda wildcards: config["samples"][wildcards.samples][1]
    output:
        "results/{samples}/{samples}_unsorted.bam"
    params:
        reference_genome=config["reference_genome"],
        bwa=config["bwa"],
        samtools=config["samtools"],
        bwa_threads=config["bwa_threads"]
    log:
        "logs/map_to_reference/{samples}_map_to_reference.txt"
    shell:
        "{params.bwa} mem -t {params.bwa_threads} -M \
        -R '@RG\\tID:no_id\\tLB:no_library\\tPL:no_platform\\tPU:no_unit\\tSM:{wildcards.samples}' \
        {params.reference_genome} \
        {input.fastq1} {input.fastq2} | {params.samtools} view -b - > {output}"


rule sort_by_coord:
    input:
        "results/{samples}/{samples}_unsorted.bam"
    output:
        "results/{samples}/{samples}_sorted.bam"
    params:
        samtools=config["samtools"]
    log:
        "logs/sort_by_coord/{samples}_sort_by_coord.txt"
    shell:
        "{params.samtools} sort -o {output} {input}"


rule mark_dups:
    input:
        "results/{samples}/{samples}_sorted.bam"
    output:
        bam="results/{samples}/{samples}_dups_marked.bam",
        metrics="results/{samples}/{samples}_marked_dup_metrics.txt"
    params:
        java=config["java"],
        picard_jar = config["picard_jar"]
    log:
        "logs/mark_dups/{samples}_mark_dups.txt"
    shell:
        "{params.java} -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx72G \
        -jar {params.picard_jar} MarkDuplicates \
        I={input} \
        O={output.bam} \
        M={output.metrics} \
        TMP_DIR={wildcards.samples}_tmp"


rule build_recalibrator_model:
    input:
        "results/{samples}/{samples}_dups_marked.bam"
    output:
        "results/{samples}/{samples}_recalibration_data.table"
    params:
        gatk=config["gatk"],
        reference_genome=config["reference_genome"],
        known_polymorphic_sites1=config["known_polymorphic_sites1"],
        known_polymorphic_sites2=config["known_polymorphic_sites2"],
        base_recalibrator_gap_open_penalty=config["base_recalibrator_gap_open_penalty"]
    log:
        "logs/build_recalibrator_model/{samples}_build_recalibrator_model.txt"
    shell:
        "{params.gatk} BaseRecalibrator \
        -I {input} \
        -R {params.reference_genome} \
        --known-sites {params.known_polymorphic_sites1} \
        --known-sites {params.known_polymorphic_sites2} \
        --bqsr-baq-gap-open-penalty {params.base_recalibrator_gap_open_penalty} \
        -O {output}"


rule apply_recalibration:
    input:
        bam="results/{samples}/{samples}_dups_marked.bam",
        model="results/{samples}/{samples}_recalibration_data.table"
    output:
        bam="results/{samples}/{samples}_recalibrated.bam",
        index="results/{samples}/{samples}_recalibrated.bai"
    params:
        gatk=config["gatk"],
        reference_genome=config["reference_genome"]
    log:
        "logs/apply_recalibration/{samples}_apply_recalibration.txt"
    shell:
        "{params.gatk} ApplyBQSR \
        -R {params.reference_genome} \
        -I {input.bam} \
        --bqsr-recal-file {input.model} \
        -O {output.bam}"


rule rename_index_files:
    input:
        "results/{samples}/{samples}_recalibrated.bai"
    output:
        "results/{samples}/{samples}_recalibrated.bam.bai"
    log:
        "logs/rename_index_files/{samples}_rename_index_files.txt"
    shell:
        "mv {input} {output}"
