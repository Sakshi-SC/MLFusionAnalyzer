# MLFusionAnalyzer
Validating Fusion Transcripts with Genomic Precision: A Machine Learning Pipeline for Plant Data


# **Overview**

Gene fusions can be identified using high-throughput techniques such as RNA sequencing (RNA-Seq) and whole genome sequencing (WGS). RNA-Seq offers insights into gene and fusion transcript expression and is often more cost-effective. However, RNA-Seq-based fusion detection tools tend to produce false positives due to read misalignment and cannot accurately identify genomic breakpoints within intronic regions.

WGS, on the other hand, can detect precise genomic breakpoints and structural variations at the DNA level. But WGS alone cannot distinguish between expressed and non-expressed fusions, often leading to the detection of non-functional genomic rearrangements.

Combining RNA-Seq and WGS provides a more comprehensive and accurate strategy for fusion detection. When chimeric transcripts identified by RNA-Seq are supported by discordant read pairs in WGS data, it provides strong evidence for true fusion events. Additionally, pinpointing exact genomic breakpoints further validates the biological relevance of these events.

However, a major challenge remains: validating predicted fusion events when matched RNA-Seq and WGS data are not available for the same biological sample or genotype—particularly common in plant genomics, where often only RNA-Seq data is accessible.

# **Our Contribution**

Here, we propose a novel and fully automated fusion gene validation pipeline built using the workflow management system **Snakemake**. Our pipeline integrates advanced machine learning (ML) and deep learning (DL) architecture to validate predicted fusion transcripts at the genomic level. Leveraging models such as **XGBoost**, **LightGBM**, and **Convolutional Neural Networks (CNNs)** **Deep neural network (DNN) models**, the pipeline is designed to handle the complexity and heterogeneity of genomic data across a wide range of plant species. Unlike prior work—such as the pipeline by Hafstað et al., which focuses on matched RNA-Seq and WGS data from human cancers—our framework is broadly adaptable and suitable for diverse plant genomic datasets, making it a valuable tool for plant genomics and transcriptomics research.

