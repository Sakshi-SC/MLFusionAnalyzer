# Snakefile

with open("paired_ids.txt", "r") as f:
    PE_SAMPLES = [line.strip() for line in f if line.strip()]

with open("single_ids.txt", "r") as f:
    SE_SAMPLES = [line.strip() for line in f if line.strip()]

rule all:
    input:
        # Paired-end final outputs
        #expand("pe_output/{sample}/star-fusion.fusion_predictions.abridged.tsv", sample=PE_SAMPLES),
        # Single-end final outputs
        #expand("se_output/{sample}/star-fusion.fusion_predictions.abridged.tsv", sample=SE_SAMPLES),
    input:
        expand("FUSION_OUTPUT/{sample}/.done_pe", sample=PE_SAMPLES),
        expand("FUSION_OUTPUT/{sample}/.done_se", sample=SE_SAMPLES)
    output:
        "Fusion_output/tempo_fusion_data_annotated.tsv"
    script:
        "2.py"

# --- STAGE 1: Download SRA files ---
rule download_sra:
    output:
        "raw/{sample}.sra"
    log:
        "logs/download_sra_{sample}.log"
    conda:
        "fusion.yaml"
    shell:
        """
        prefetch {wildcards.sample} -o {output} --verbose &> {log} && \
        vdb-validate {output} &>> {log}
        """

# --- STAGE 2: Convert SRA to FASTQ ---
rule sra_to_fastq_pe:
    input:
        "raw/{sample}.sra"
    output:
        r1="pe_fastq/{sample}_1.fastq",
        r2="pe_fastq/{sample}_2.fastq"
    log:
        "logs/fastq_dump_pe_{sample}.log"
    conda:
        "fusion.yaml"
    threads: 4
    shell:
        """
        mkdir -p pe_fastq && \
        fastq-dump --split-files --outdir pe_fastq {input} &> {log}
        """

rule sra_to_fastq_se:
    input:
        "raw/{sample}.sra"
    output:
        "se_fastq/{sample}.fastq"
    log:
        "logs/fastq_dump_se_{sample}.log"
    conda:
        "fusion.yaml"
    threads: 2
    shell:
        """
        mkdir -p se_fastq && \
        fastq-dump --outdir se_fastq {input} &> {log}
        """

# --- STAGE 3: Quality Trimming ---
rule trim_galore_pe:
    input:
        r1="pe_fastq/{sample}_1.fastq",
        r2="pe_fastq/{sample}_2.fastq"
    output:
        r1="pe_trimmed/{sample}_1_val_1.fq.gz",
        r2="pe_trimmed/{sample}_2_val_2.fq.gz"
    log:
        "logs/trim_galore_pe_{sample}.log"
    conda:
        "fusion.yaml"
    threads: 8
    shell:
        """
        mkdir -p pe_trimmed && \
        trim_galore --paired --gzip --output_dir pe_trimmed {input.r1} {input.r2} &> {log}
        """

rule trim_galore_se:
    input:
        "se_fastq/{sample}.fastq"
    output:
        "se_trimmed/{sample}_trimmed.fq.gz"
    log:
        "logs/trim_galore_se_{sample}.log"
    conda:
        "fusion.yaml"
    threads: 4
    shell:
        """
        mkdir -p se_trimmed && \
        trim_galore --gzip --output_dir se_trimmed {input} &> {log}
        """

# --- STAGE 4: STAR-Fusion Analysis ---
rule star_fusion_pe:
    input:
        r1="pe_trimmed/{sample}_1_val_1.fq.gz",
        r2="pe_trimmed/{sample}_2_val_2.fq.gz"
    output:
        "pe_output/{sample}/star-fusion.fusion_predictions.abridged.tsv"
    log:
        "logs/star_fusion_pe_{sample}.log"
    conda:
        "fusion.yaml"
    threads: 8
    shell:
        """
        mkdir -p pe_output/{wildcards.sample} && \
        gunzip -c {input.r1} > pe_trimmed/{wildcards.sample}_1_val_1.fq && \
        gunzip -c {input.r2} > pe_trimmed/{wildcards.sample}_2_val_2.fq && \
        STAR-Fusion \
            --left_fq pe_trimmed/{wildcards.sample}_1_val_1.fq \
            --right_fq pe_trimmed/{wildcards.sample}_2_val_2.fq \
            --genome_lib_dir arabidopsis_lib \
            --examine_coding_effect \
            -O pe_output/{wildcards.sample} &> {log} && \
        rm pe_trimmed/{wildcards.sample}_1_val_1.fq pe_trimmed/{wildcards.sample}_2_val_2.fq
        """

rule star_fusion_se:
    input:
        "se_trimmed/{sample}_trimmed.fq.gz"
    output:
        "se_output/{sample}/star-fusion.fusion_predictions.abridged.tsv"
    log:
        "logs/star_fusion_se_{sample}.log"
    conda:
        "fusion.yaml"
    threads: 8
    shell:
        """
        mkdir -p se_output/{wildcards.sample} && \
        gunzip -c {input} > se_trimmed/{wildcards.sample}_val.fq && \
        STAR-Fusion \
            --left_fq se_trimmed/{wildcards.sample}_val.fq \
            --genome_lib_dir arabidopsis_lib \
            --examine_coding_effect \
            -O se_output/{wildcards.sample} &> {log} && \
#STAGE 5: Move STAR-Fusion Output to Unified Directory ---
rule move_star_fusion_pe:
    input:
        "pe_output/{sample}/star-fusion.fusion_predictions.abridged.tsv"
    output:
        touch("FUSION_OUTPUT/{sample}/.done_pe")
    run:
        import shutil, os
        src_dir = f"pe_output/{wildcards.sample}"
        dest_dir = f"FUSION_OUTPUT/{wildcards.sample}"
        if not os.path.exists("FUSION_OUTPUT"):
            os.makedirs("FUSION_OUTPUT")
        shutil.move(src_dir, dest_dir)

rule move_star_fusion_se:
    input:
        "se_output/{sample}/star-fusion.fusion_predictions.abridged.tsv"
    output:
        touch("FUSION_OUTPUT/{sample}/.done_se")
    run:
        import shutil, os
        src_dir = f"se_output/{wildcards.sample}"
        dest_dir = f"FUSION_OUTPUT/{wildcards.sample}"
        if not os.path.exists("FUSION_OUTPUT"):
            os.makedirs("FUSION_OUTPUT")
        shutil.move(src_dir, dest_dir)
        """


# --- Execution Order Control ---
ruleorder:
ruleorder: download_sra > sra_to_fastq_pe
ruleorder: sra_to_fastq_pe > trim_galore_pe
ruleorder: trim_galore_pe > star_fusion_pe

ruleorder: download_sra > sra_to_fastq_se
ruleorder: sra_to_fastq_se > trim_galore_se
ruleorder: trim_galore_se > star_fusion_s
