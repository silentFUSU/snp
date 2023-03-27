#! /usr/bin/env bash
## Snakefile
####################

ref_data="/storage/zhangyanxiaoLab/suzhuojie/ref_data/"
bam_data="/storage/zhangyanxiaoLab/suzhuojie/software/varCA/data/"
result="/storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/result/merge_data/before_doubletfinder/gatk/"
raw_data="/storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/data/rawdata/bam/"
Samples=['S04']
subset_bam="/storage/zhangyanxiaoLab/suzhuojie/software/subset-bam_linux"
chr="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
rule all:
    input:
        expand(raw_data+"normal_chromosome/{sample}.bam",sample=Samples),
        expand(result+"recal_table/{sample}_recal_data.table1",sample=Samples),
        expand(raw_data+"fixed/{sample}_marked_fixed.bam",sample=Samples),
        expand(result+"recal_table/{sample}_recal_data.table2",sample=Samples),
        expand(result+"plot/{sample}_AnalyzeCovariates.pdf",sample=Samples),
        expand(raw_data+"fixed_cancer/{sample}_marked_fixed_cancer.bam",sample=Samples),
        expand(result+"vcf/{sample}_tumor.vcf.gz",sample=Samples)        
        
rule normalchromosome:
    input:
        bam_data+"{sample}.bam"
    output:
        raw_data+"normal_chromosome/{sample}.bam"
    shell:
        "samtools view -b {input} {chr} > {output}"

rule BaseRecalibrator:
    input:
        raw_data+"normal_chromosome/{sample}.bam"
    output:
        result+"recal_table/{sample}_recal_data.table1"
    shell: 
        "gatk --java-options '-Xmx20G' BaseRecalibrator "
        "-I {input} "
        "-R {ref_data}hg38.fa "
        "--known-sites {ref_data}gatk/1000G_phase1.snps.high_confidence.hg38.vcf.gz "
        "--known-sites {ref_data}gatk/dbsnp_146.hg38.vcf.gz "
        "--known-sites {ref_data}gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz "
        "-O {output}"

rule ApplyBQSR:
    input:
        bam=bam_data+"{sample}.bam",
        table=result+"recal_table/{sample}_recal_data.table1"
    output:
        raw_data+"fixed/{sample}_marked_fixed.bam"
    shell:
        "gatk  --java-options '-Xmx20G' ApplyBQSR " 
        "-R {ref_data}hg38.fa "
        "-I {input.bam} "
        "--bqsr-recal-file  {input.table} "
        "-O {output}"

rule BaseRecalibrator2:
    input:
        raw_data+"fixed/{sample}_marked_fixed.bam"
    output:
        result+"recal_table/{sample}_recal_data.table2"
    shell: 
        "gatk --java-options '-Xmx20G' BaseRecalibrator "
        "-I {input} "
        "-R {ref_data}hg38.fa "
        "--known-sites {ref_data}gatk/1000G_phase1.snps.high_confidence.hg38.vcf.gz "
        "--known-sites {ref_data}gatk/dbsnp_146.hg38.vcf.gz "
        "--known-sites {ref_data}gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz "
        "-O {output}"

rule AnalyzeCovariates:
    input:
        before=result+"recal_table/{sample}_recal_data.table1",
        after=result+"recal_table/{sample}_recal_data.table2"
    output:
        result+"plot/{sample}_AnalyzeCovariates.pdf"
    shell:
        "gatk --java-options '-Xmx20G' AnalyzeCovariates " 
        "-before {input.before} "
        "-after {input.after} "
        "-plots {output}"

rule subset_tumor:
    input:
        bam= raw_data+"fixed/{sample}_marked_fixed.bam",
        barcode= result+"barcode_cancer/{sample}.txt"
    output:
        raw_data+"fixed_cancer/{sample}_marked_fixed_cancer.bam"
    shell:
        "{subset_bam} --bam {input.bam} --bam-tag CB --cell-barcodes {input.barcode} cores 4 --log-level info --out-bam {output}"

rule mutect2:
    input:
        bam=raw_data+"fixed_cancer/{sample}_marked_fixed_cancer.bam"
    output:
        result+"vcf/{sample}_tumor.vcf.gz"
    shell:
        "gatk --java-options '-Xmx20G'  Mutect2 "
        "-R {ref_data}hg38.fa "
        "-I {input} "
        "-O {output}"
