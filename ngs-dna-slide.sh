#!/bin/bash
####
##  The complete NGS pipeline for sequence to variants with haplotypecaller and QC metrics
##   Made to take a given step and progress from there, good for faults and stops
##  Written by Tristan M. Carland, PhD
####
if [ $# -lt 3 ]; then
  echo ""
  echo " Tristan's NGS Scripts - Full NGS Slide v1.0 (GATK 3.3-0)"
  echo ""
  echo "  usage: `basename $0` <config> <sample> <starting-sub-routine> [optional-bed-file]"
  echo ""
  echo "   Functions                Expected <sample>"
  echo "     bwa-mem                  -r1.fq.gz     + -r2.fq.gz"
  echo "     sambamba-sort            -fm-md.bam"
  echo "     samtools-index           -s.bam"
  echo "     gatk-rtc                 -s.bam        + -s.bam.bai"
  echo "     gatk-ir                  -s.bam        + -s.bam.bai     + -s.intervals"
  echo "     gatk-br                  -ir.bam       + -ir.bam.bai"
  echo "     gatk-pr                  -ir.bam       + -ir.bam.bai    + -ir.bam.table"
  echo "     gatk-callable            -ready.bam    + -ready.bam.bai"
  echo "     gatk-hc                  -ready.bam    + -ready.bam.bai"
  echo "     gatk-snp-select          -raw.vcf"
  echo "     gatk-snp-filt            -raw-snp.vcf"
  echo "     gatk-ind-select          -raw.vcf      + -mark-snp.vcf"
  echo "     gatk-ind-filt            -raw-ind.vcf  + -mark-snp.vcf"
  echo "     gatk-variant-combine     -mark-ind.vcf + -mark-snp.vcf"
  echo "     gatk-variant-pass        -mark.vcf"
  echo "     gatk-variant-eval        -hard.vcf"
  echo ""
  exit 1
fi

##
# Parsing of command line input variables and various bundled library files
config=$1
sample=$2
functn=$3
tmp=temp-$sample

if($bedfile ne ""); then
	gatk_regions=$4
else
	gatk_regions=" -L $gatk_regions"
fi

##
# Most of the variables, as well as the finishCheck function, are saved in a config file
#  in particular, the config file now saves the number of threads (job_thr) and job_ram settings
source $config
echo "NGS-CONFIG - Sourcing $config for configuration"
echo "NGS-CONFIG - Using Reference $ref"
echo ""

# ## Starting timestamp
echo 'NGS-TIMESTAMP' `date`

main () # called from the bottom so that the other functions are known
{
	# takes the chosen function and actually calls that function, nice and simple
	case "$functn" in
		bwa-mem)					bwa-mem ;;
		sambamba-sort)				sambamba-sort ;;
		samtools-index)				samtools-index ;;
		gatk-rtc)					gatk-rtc ;;
		gatk-ir)					gatk-ir ;;
		gatk-br)					gatk-br ;;
		gatk-pr)					gatk-pr ;;
		gatk-callable)				gatk-callable ;;
		gatk-hc)					gatk-hc ;;
		gatk-snp-select)			gatk-snp-select ;;
		gatk-snp-filt)				gatk-snp-filt ;;
		gatk-ind-select)			gatk-ind-select ;;
		gatk-ind-filt)				gatk-ind-filt ;;
		gatk-variant-combine)		gatk-variant-combine ;;
		gatk-variant-pass)			gatk-variant-pass ;;
		gatk-variant-eval)			gatk-variant-eval ;;
	esac
}

##
# problem of note - with the piping, there isn't proper error propogation
bwa-mem ()
{
	echo "NGS-CMD - Begin bwa mem pipe" # includes basic readgroup info from sample name	
	bwa mem -t $job_thr -M \
	  -R '@RG\tID:'$sample'\tSM:'$sample'\tPL:illumina\tPU:ns500\tLB:1234' \
	  $ref $sample-r1.fq.gz $sample-r2.fq.gz \
	  | samblaster -r \
	  | samtools fixmate - -O bam $sample-fm-md.bam

	finishCheck $? "bwa mem pipe"

	sambamba-sort # call the next tool
}

sambamba-sort ()
{
	echo "NGS-CMD - Begin sambamba sort" # creates many many temporary files
	sambamba sort -m $job_ram -t $job_thr -p -o $sample-s.bam $sample-fm-md.bam

	finishCheck $? "sambamba sort"

	## delete more intermediate files
	rm $sample-fm-md.bam

	samtools-index # call the next tool
}

samtools-index ()
{
	echo "NGS-CMD - Begin samtools indexing"
	samtools index $sample-s.bam

	finishCheck $? "samtools indexing"

	gatk-rtc # call the next tool
}

gatk-rtc ()
{
	echo "NGS-CMD - Start GATK RealignerTargetCreator" # known should include known indels
	java -Xmx$job_ram -jar $gatk \
	  -T RealignerTargetCreator -nt $job_thr \
	  -R $ref \
	  --known $dbsnp_vcf \
	  --known $indel_vcf \
	  -I $sample-s.bam \
	  -o $sample-s.intervals

	finishCheck $? "GATK RealignerTargetCreator"

	gatk-ir # call the next tool
}

gatk-ir ()
{
	echo "NGS-CMD - Start GATK IndelRealigner" # knownAlleles should include known indels
	java -Xmx$job_ram -jar $gatk \
	  -T IndelRealigner \
	  -R $ref \
	  --knownAlleles $dbsnp_vcf \
	  --knownAlleles $indel_vcf \
	  -filterNoBases \
	  -filterMBQ \
	  -filterRNC \
	  -targetIntervals $sample-s.intervals \
	  -I $sample-s.bam \
	  -o $sample-ir.bam

	finishCheck $? "GATK IndelRealigner"

	rm $sample-s.bam
	rm $sample-s.bam.bai

	gatk-br # call the next tool
}

gatk-br ()
{
	echo "NGS-CMD - Start GATK BaseRecalibrator" # knownSites should include any known variants
	java -Xmx$job_ram -jar $gatk \
	  -T BaseRecalibrator -nct $job_thr \
	  -R $ref \
	  --knownSites $dbsnp_vcf \
	  --knownSites $indel_vcf \
	  -I $sample-ir.bam \
	  -o $sample-ir.table

	finishCheck $? "GATK BaseRecalibrator"

	gatk-pr # call the next tool
}

gatk-pr ()
{
	echo "NGS-CMD - Start GATK PrintReads" # uses info from the previous steps
	java -Xmx$job_ram -jar $gatk \
	  -T PrintReads -nct $job_thr \
	  -R $ref \
	  --BQSR $sample-ir.table \
	  -I $sample-ir.bam \
	  -o $sample-ready.bam

	finishCheck $? "GATK PrintReads"

	rm $sample-ir.bam
	rm $sample-ir.bai
	rm $sample-s.intervals
	rm $sample-ir.table

	gatk-hc # call the next tool
}

# gatk-callable ()
# {
# 	echo "NGS-CMD - Start GATK Callable Loci tool" # can we call certain loci
# 	echo $clia_list_bed
# 	java -jar $gatk \
# 	  -T CallableLoci \
# 	  -R $ref \
# 	  -I $sample-ready.bam \
# 	  -L $care_list_bed \
# 	  -summary qc-$sample-5-call.sum \
# 	  -o qc-$sample-call.bed

# 	finishCheck $? "GATK Callable Loci AveraList"

# 	rm qc-$sample-call.bed

# 	gatk-hc # call the next tool
# }

gatk-hc ()
{
	echo "NGS-CMD - Start GATK HaplotypeCaller" # the genotyper
	java -Xmx$job_ram -jar $gatk \
	 -T HaplotypeCaller -nct $job_thr \
	 -R $ref \
	 --dbsnp $dbsnp_vcf \
	 -I $sample-ready.bam \
	 -stand_call_conf 30 \
	 -stand_emit_conf 10 \
	 -o $sample-raw.vcf

	finishCheck $? "GATK HaplotypeCaller"

	gatk-snp-select # call the next tool
}

gatk-snp-select ()
{
	echo "NGS-CMD - Begin SNP selection"
	java -Xmx2g -jar $gatk \
	 -T SelectVariants \
	 -R $ref \
	 -V $sample-raw.vcf \
	 -selectType SNP \
	 -o $sample-raw-snp.vcf

	finishCheck $? "SNP selection"

	gatk-snp-filt # call the next tool
}

gatk-snp-filt ()
{
	echo "NGS-CMD - Begin SNP filtration"
	java -Xmx2g -jar $gatk \
	 -T VariantFiltration \
	 -R $ref \
	 -V $sample-raw-snp.vcf \
	 --filterExpression "(vc.hasAttribute('QD') && QD < 2.0) || (vc.hasAttribute('FS') && FS > 60.0) || (vc.hasAttribute('MQ') && MQ < 40.0) || (vc.hasAttribute('MappingQualityRankSum') && MappingQualityRankSum < -12.5) || (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)" \
	 --filterName "vda-hard14" \
	 -o $sample-mark-snp.vcf

	finishCheck $? "SNP filtration"

	rm $sample-raw-snp.vcf
	rm $sample-raw-snp.vcf.idx

	gatk-ind-select # call the next tool
}

gatk-ind-select ()
{
	echo "NGS-CMD - Begin InDel selection" # Extract the Indels from the call se
	java -Xmx2g -jar $gatk \
	 -T SelectVariants \
	 -R $ref \
	 -V $sample-raw.vcf \
	 -selectType INDEL \
	 -o $sample-raw-indel.vcf

	finishCheck $? "InDel selection"

	gatk-ind-filt # call the next tool
}

gatk-ind-filt ()
{
	echo "NGS-CMD - Begin InDel filtering" # Apply the filter to the Indel call set
	java -Xmx2g -jar $gatk \
	 -T VariantFiltration \
	 -R $ref \
	 -V $sample-raw-indel.vcf \
	 --filterExpression "(vc.hasAttribute('QD') && QD < 2.0) || (vc.hasAttribute('FS') && FS > 200.0) || (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0)" \
	 --filterName "vda-hard14" \
	 -o $sample-mark-indel.vcf

	finishCheck $? "InDel filtration"

	rm $sample-raw-indel.vcf
	rm $sample-raw-indel.vcf.idx

	gatk-variant-combine # call the next tool
}

gatk-variant-combine ()
{
	echo "NGS-CMD - Start Combining Variants" # Merge the filtered sets
	java -Xmx2g -jar $gatk \
	 -T CombineVariants \
	 -R $ref \
	 --variant:snp $sample-mark-snp.vcf \
	 --variant:ind $sample-mark-indel.vcf \
	 --genotypemergeoption PRIORITIZE \
	 --rod_priority_list snp,ind \
	 -o $sample-mark.vcf

	finishCheck $? "Combine Variants"

	rm $sample-raw.vcf
	rm $sample-raw.vcf.idx
	rm $sample-mark-snp.vcf
	rm $sample-mark-snp.vcf.idx
	rm $sample-mark-indel.vcf
	rm $sample-mark-indel.vcf.idx

	gatk-variant-pass # call the next tool
}

gatk-variant-pass ()
{
	echo "NGS-CMD - Start PASS Collection" # Apply a simple PASS/FAIL to get the final set
	java -Xmx2g -jar $gatk \
	 -T SelectVariants \
	 -R $ref \
	 -V $sample-mark.vcf \
	 --excludeFiltered \
	 --excludeNonVariants \
	 -o $sample-hard.vcf

	finishCheck $? "PASS Collection"

	echo " To clean up all intermediate files, just delete the temp-$2-filt directory"

	gatk-variant-eval # call the next tool
}

gatk-variant-eval ()
{
	echo "NGS-CMD - Variant Evaluation"
	java -Xmx4g -jar $gatk \
	 -T VariantEval -nt $job_thr \
	 -R $ref \
	 --dbsnp $dbsnp_vcf \
	 --goldStandard $indel_vcf \
	 --doNotUseAllStandardModules \
	 --doNotUseAllStandardStratifications \
	 --evalModule VariantSummary \
	 --evalModule CountVariants \
	 --evalModule IndelSummary \
	 --eval:$sample-mark $sample-mark.vcf \
	 -o qc-$sample-varEval.csv

	finishCheck $? "Variant Evaluation"
}

main # call the main function - this way it knows of the other functions

exit 0 # end shell scripts with an exit for easier commenting at the en
##############################################################################################
# Notes #
#########

### A note on parallelism in GATK ###
http://gatkforums.broadinstitute.org/discussion/1975/recommendations-for-parallelizing-gatk-tools

nt  = --job_threads - the number of data threads sent to the processor
nct = -- num_cpu_threads_per_data_thread - the number of CPU threads allocated to each data thread
- not that this changes a whole lot...


Available stratification modules:
 to use, add -ST or --stratIntervals
(Standard modules are starred) 
	AlleleCount 
	AlleleFrequency 
	CompRod* 
	Contig 
	CpG 
	Degeneracy 
	EvalRod* 
	Filter 
	FunctionalClass 
	IndelSize 
	IntervalStratification 
	JexlExpression* 
	Novelty* 
	OneBPIndel 
	Sample 
	SnpEffPositionModifier 
	TandemRepeat 
	VariantType 
 
Available evaluation modules:
 to use, add -EV or --evalModule
(Standard modules are starred) 
	CompOverlap* 
	CountVariants* 
	IndelLengthHistogram* 
	IndelSummary* 
	MendelianViolationEvaluator 
	MultiallelicSummary* 
	PrintMissingComp 
	ThetaVariantEvaluator 
	TiTvVariantEvaluator* 
	ValidationReport* 
	VariantSummary* 
