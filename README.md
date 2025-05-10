# MLFusionRoots
Validating Fusion Transcripts with Genomic Precision: A Machine Learning Pipeline for Plant Data


Overview
Gene fusions can be identified using high-throughput techniques such as RNA sequencing (RNA-Seq) and whole genome sequencing (WGS). RNA-Seq offers insights into gene and fusion transcript expression and is often more cost-effective. However, RNA-Seq-based fusion detection tools tend to produce false positives due to read misalignment and cannot accurately identify genomic breakpoints within intronic regions.

WGS, on the other hand, can detect precise genomic breakpoints and structural variations at the DNA level. But WGS alone cannot distinguish between expressed and non-expressed fusions, often leading to the detection of non-functional genomic rearrangements.

Combining RNA-Seq and WGS provides a more comprehensive and accurate strategy for fusion detection. When chimeric transcripts identified by RNA-Seq are supported by discordant read pairs in WGS data, it provides strong evidence for true fusion events. Additionally, pinpointing exact genomic breakpoints further validates the biological relevance of these events.

However, a major challenge remains: validating predicted fusion events when matched RNA-Seq and WGS data are not available for the same biological sample or genotype—particularly common in plant genomics, where often only RNA-Seq data is accessible
