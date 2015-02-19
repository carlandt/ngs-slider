#!/bin/bash
####
##  The system related information needed by my ngs suite of tools
##  Written by Tristan M. Carland, PhD (Avera, JCVI)
####


####
## Ensembl human_g1k_v37_decoy reference related files (Grch37)
refdir="/data/database/GATK/Grch37"
ref="$refdir/ref.fasta"

# note that the 142 version has been problematic, 138 works so far
dbsnp_vcf="$refdir/dbsnp_138.b37.vcf"
dbsnp_bed="$refdir/dbsnp_138.b37.bed"
indel_vcf="$refdir/Mills_and_1000G_gold_standard.indels.b37.vcf"
pharmGKB_bed="$refdir/genes-pharmgkb.bed"

temp="/mnt/tristan/temp-$sample"

cosmic="/data/storage/b37/b37_cosmic_v54_120711.vcf"

care_list_bed="$refdir/clia-list.bed"

####
##
job_thr=6
job_ram=12G

node_thr=8
node_ram=16G

####
## Hg19
#refdir="/data/database/GATK/hg19"
#ref="/data/database/GATK/hg19/ucsc.hg19.fasta"
#dbsnp="/data/database/GATK/hg19/dbsnp_137.hg19.vcf"
#indel="/data/database/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf"


####
## Paths to tools, unnecessary with proper environmental variables
tools="/opt/software"
mutect="$tools/mutect/1.1.4/muTect-1.1.4.jar"
gatk="$tools/gatk/3.3-0/GenomeAnalysisTK-3.3-0.jar"
MarkDup="$tools/picard/picard-tools-1.121/MarkDuplicates.jar"
PicardAlignMetrics="$tools/picard/picard-tools-1.121/CollectAlignmentSummaryMetrics.jar"


####
## Function used to check termination variables and generate timestamps and error messages
finishCheck ()
{
  if [ "$1" -ne "0" ]
    then
      echo "NGS-CMD - $2 Failed" 1>&2
      echo 'NGS-TIMESTAMP' `date`
      exit 1
    else
      echo "NGS-CMD - $2 Completed Successfully" 1>&2
      echo 'NGS-TIMESTAMP' `date`
  fi
}
