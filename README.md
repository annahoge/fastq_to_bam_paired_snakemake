# fastq_to_bam_paired_snakemake
A snakemake to convert paired fastq files to analysis-ready bams, following GATK best practices.

This repository contains all the folders needed to run the snakemake.  results/ and logs/cluster/ currently contain placeholder files because GitHub does not allow empty folders, but after copying this repository structure to the server, these placeholder files can be removed from the server.

To run the snakemake, update config/samples.yaml and config/config.yaml as per the directions in those files, and then follow the instructions in fastq_to_bam_paired.snakefile.  The jobs will be launched to the cluster at the Fred Hutchinson Cancer Research Center.
