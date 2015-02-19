# ngs-slider
My homebrew pipeline for taking next-gen sequencing data and finding decent variants!

Designed to follow the Best Practices page of the GATK website - [link](https://www.broadinstitute.org/gatk/guide/best-practices?bpm=DNAseq)




When given a fastq file in bioinformatics, the goal is to run through a reasonable number of tools and eventually produce a vcf file. This program allows you to start from any one of those steps (provided you meet the naming conventions) and slide through the rest of the pipeline.

This tool is a combination of two scripts and a series of bioinformatic tools.
ngs-slide.sh - used to run the slide, try running it for an explanation
ngs-config.sh - a config file for your environment, includes hardware and software params

Good luck, I'll post more here as time goes on.

### Running it from the command line
The idea is to choose your starting tool/step, and it will slide on from there. If you're at the very beginning, it assumes that you have two gzipped fastq files from a single sample that have been properly demultiplexed (eg sampleA-r1.fq.gz and sampleA-r2.fq.gz). For that you would type:
```ngs-slide.sh sampleA bwa-mem```

Here are the other options

| Functions              | Expected <sample> |
|------------------------|-------------------------- |
| bwa-mem                | -r1.fq.gz + -r2.fq.gz |
| samtools-fixmate       | -1-raw.sam |
| sambamba-sort          | -2-fm.bam |
| picard-mark-duplicates | -3-s.bam |
| samtools-index         | -4-md.bam |
| gatk-rtc               | -4-md.bam + -4-md.bam.bai |
| gatk-ir                | -4-md.bam + -4-md.bam.bai  + -4.intervals |
| gatk-br                | -5-ir.bam + -5-ir.bam.bai |
| gatk-pr                | -5-ir.bam + -5-ir.bam.bai  + -5-ir.bam.table |
| gatk-callable          | -ready.bam + -ready.bam.bai |
| gatk-hc                | -ready.bam + -ready.bam.bai |
| gatk-snp-select        | -raw.vcf |
| gatk-snp-filt          | -raw-snp.vcf |
| gatk-ind-select        | -raw.vcf      + -mark-snp.vcf |
| gatk-ind-filt          | -raw-ind.vcf  + -mark-snp.vcf |
| gatk-variant-combine   | -mark-ind.vcf + -mark-snp.vcf |
| gatk-variant-pass      | -mark.vcf |
| gatk-variant-eval      | -hard.vcf |

### Installation

### Required Tools
