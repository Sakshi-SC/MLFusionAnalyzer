# Snakefile

with open("paired_ids.txt", "r") as f:
    PE_SAMPLES = [line.strip() for line in f if line.strip()]

with open("single_ids.txt", "r") as f:
    SE_SAMPLES = [line.strip() for line in f if line.strip()]

rule all:
    input:
        expand("trimmed/{sample}_1_val_1.fq.gz", sample=PE_SAMPLES) if PE_SAMPLES else [],
        expand("trimmed/{sample}_2_val_2.fq.gz", sample=PE_SAMPLES) if PE_SAMPLES else [],
        expand("trimmed/{sample}_val.fq.gz", sample=SE_SAMPLES) if SE_SAMPLES else [],
         expand("output/{sample}/pe_star-fusion.fusion_predictions.abridged.coding_effect.tsv", sample=PE_SAMPLES),
        expand("output/{sample}/se_star-fusion.fusion_predictions.abridged.coding_effect.tsv", sample=SE_SAMPLES)

rule download_sra:
    output:
        sra="raw/{sample}.sra"
    log:
        "logs/fastq_dump_{sample}.log"
    conda:
        "fusion.yaml"
    shell:
        """
        prefetch {wildcards.sample} -o {output.sra} --verbose &> {log} && \
        vdb-validate {output.sra} &>> {log}
        """

rule sra_to_fastq_pe:
    input:
        "raw/{sample}.sra"
    output:
        r1="fastq/{sample}_1.fastq",
        r2="fastq/{sample}_2.fastq"
    log:
        "logs/fastq_dump_pe_{sample}.log"
    conda:
        "fusion.yaml"
    threads: 4
    shell:
        """
        fastq-dump --split-files --outdir fastq/ {input} &> {log}
        """

rule sra_to_fastq_se:
    input:
        "raw/{sample}.sra"
    output:
        "fastq/{sample}.fastq"
    log:
        "logs/fastq_dump_se_{sample}.log"
    conda:
        "fusion.yaml"
    threads: 2
    shell:
        """
        fastq-dump --outdir fastq/ {input} &> {log}
        """

rule trim_galore_pe:
    input:
        r1="fastq/{sample}_1.fastq",
        r2="fastq/{sample}_2.fastq"
    output:
        r1_trim="trimmed/{sample}_1_val_1.fq.gz",
        r2_trim="trimmed/{sample}_2_val_2.fq.gz"
    log:
        "logs/trim_galore_pe_{sample}.log"
    conda:
        "fusion.yaml"
    threads: 8  # Increased from 2 to 4 (16 cores available)
    shell:
        """
        trim_galore --paired --fastqc --gzip --output_dir trimmed {input.r1} {input.r2} &> {log} && \
        [ -f {output.r1_trim} ] && [ -f {output.r2_trim} ] || \
        (echo "Error: One or both output files missing" >> {log} && exit 1)
        """

rule trim_galore_se:
    input:
        "fastq/{sample}.fastq"
    output:
        "trimmed/{sample}_val.fq.gz"
    log:
        "logs/trim_galore_se_{sample}.log"
    conda:
        "fusion.yaml"
    threads: 1
    shell:
        """
        trim_galore --fastqc --gzip --output_dir trimmed {input} &> {log} && \
        [ -f trimmed/{wildcards.sample}_trimmed.fq.gz ] && mv trimmed/{wildcards.sample}_trimmed.fq.gz {output} || \
        (echo "Error: Expected output trimmed/{wildcards.sample}_trimmed.fq.gz not found" >> {log} && exit 1)
        """

rule star_fusion_pe:
    input:
        pe_r1="trimmed/{sample}_1_val_1.fq.gz",
        pe_r2="trimmed/{sample}_2_val_2.fq.gz"
    output:
        "output/{sample}/pe_star-fusion.fusion_predictions.abridged.coding_effect.tsv"
    log:
        "logs/star_fusion_{sample}.log"
    conda:
        "fusion.yaml"
    threads: 8
    shell:
        """
        mkdir -p output/{wildcards.sample} && \
        gunzip -c {input.pe_r1} > trimmed/{wildcards.sample}_1_val_1.fq && \
        gunzip -c {input.pe_r2} > trimmed/{wildcards.sample}_2_val_2.fq && \
        STAR-Fusion \
            --left_fq trimmed/{wildcards.sample}_1_val_1.fq \
            --right_fq trimmed/{wildcards.sample}_2_val_2.fq \
            --genome_lib_dir arabidopsis_lib \
            -O output/{wildcards.sample} \
            #--FusionInspector validate \
            --examine_coding_effect \
            --denovo_reconstruct &> {log} && \
        mv output/{wildcards.sample}/star-fusion.fusion_predictions.abridged.coding_effect.tsv output/{wildcards.sample}/pe_star-fusion.fusion_predictions.abridged.coding_effect.tsv && \
rm trimmed/{wildcards.sample}_1_val_1.fq trimmed/{wildcards.sample}_2_val_2.fq
        """

rule star_fusion_se:
    input:
        "trimmed/{sample}_val.fq.gz"
    output:
        "output/{sample}/se_star-fusion.fusion_predictions.abridged.coding_effect.tsv"
    log:
        "logs/star_fusion_se_{sample}.log"
    conda:
        "fusion.yaml"
    threads: 8
    shell:
        """
        mkdir -p output/{wildcards.sample} && \
        gunzip -c {input} > trimmed/{wildcards.sample}_val.fq && \
        STAR-Fusion \
            --left_fq trimmed/{wildcards.sample}_val.fq \
            --genome_lib_dir arabidopsis_lib \
            -O output/{wildcards.sample} \
            #--FusionInspector validate \
            --examine_coding_effect \
            --denovo_reconstruct &> {log} && \
        mv output/{wildcards.sample}/star-fusion.fusion_predictions.abridged.coding_effect.tsv output/{wildcards.sample}/se_star-fusion.fusion_predictions.abridged.coding_effect.tsv && \
rm trimmed/{wildcards.sample}_val.fq
        """

ruleorder: sra_to_fastq_pe > sra_to_fastq_se # trim_galore_pe > trim_galore_se > star_fusion_pe > star_fusion_se
