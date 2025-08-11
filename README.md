# amplicall
A variant calling pipeline optimized for ONT bacterial gene amplicons

This project was developed for the purpose of benchmarking variant callers on ONT bacterial gene amplicons, specifically the *amuc_1412* gene in *A. muciniphila*.

## Software Versions

| Software | Version | Purpose |
|----------|---------|---------|
| NanoFilt | v2.8.0 | Quality filtering for long-read sequencing data |
| NanoPlot | v1.41.0 | Quality visualization and statistics for long-read data |
| Minimap2 | v2.28 | Sequence alignment for long-read sequencing |
| Samtools | v1.20 | SAM/BAM file manipulation and processing |
| Freebayes | v1.3.7 | Bayesian genetic variant detector |
| LoFreq | v2.1.5 | Low-frequency variant caller |
| Clair3 | v1.0.10 | Deep learning-based variant caller |
| BCFtools | v1.22 | VCF/BCF file manipulation and analysis |
| RTG-tools | v3.13 | Variant calling evaluation and comparison tools |

## Script Descriptions
The scripts should be executed in the following order:

| Script Name | Function | Input | Output | Execution Location |
|-------------|----------|-------|---------|-----------------|
| `alignment.sh` | Filters raw reads and aligns them to a reference sequence | Raw FASTQ files | Sorted and indexed BAM files, quality reports | HPC |
| `freebayes.sh` | Calls variants using Freebayes | BAM files | VCF files | HPC |
| `lofreq.sh` | Calls variants using LoFreq | BAM files | VCF files | Local PC |
| `clair3.sh` | Calls variants using Clair3 | BAM files | VCF files | HPC |
| `filter.sh` | Filters variant calls based on QUAL scores | VCF files | Filtered VCF files | HPC |
| `format-fb.sh` | Formats Freebayes VCF files | Filtered VCF files | Reformatted and compressed VCF files | Local PC |
| `format-lf.sh` | Formats LoFreq VCF files | Filtered VCF files | Reformatted and compressed VCF files | Local PC |
| `format-c3.sh` | Formats Clair3 VCF files | Filtered VCF files | Reformatted and compressed VCF files | Local PC |
| `vcfeval.sh` |  Evaluates the performance of variant callers | Compressed VCF files | Evaluation results | Local PC |
| `performance_plots.R` | Generates plots comparing the perfromance of variant callers | Evaluation results | Performance plots | Local PC |

