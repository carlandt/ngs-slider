#!/bin/bash
####
##  The complete NGS pipeline for sequence to variants with haplotypecaller and QC metrics
##   Made to take a given step and progress from there, good for faults and stops
##   Also takes a config file as input, see bottom of script for template
##  Written by Tristan M. Carland, PhD (Avera, JCVI)
####
if [ $# -lt 2 ]; then
  echo ""
  echo " Tristan's NGS Scripts - Full NGS Slide v1.2 (GATK 3.3-0)"
  echo ""
  echo "  usage: `basename $0` <sample> <starting-sub-routine>"
  echo ""
  echo "   Functions                Expected <sample>"
  echo "     bwa-mem                  -r1.fq.gz     + -r2.fq.gz"
  echo "     samtools-fixmate         -1-raw.sam"
  echo "     sambamba-sort            -2-fm.bam"
  echo "     picard-mark-duplicates   -3-s.bam"
  echo "     samtools-index           -4-md.bam"
  echo "     gatk-rtc                 -4-md.bam     + -4-md.bam.bai"
  echo "     gatk-ir                  -4-md.bam     + -4-md.bam.bai  + -4.intervals"
  echo "     gatk-br                  -5-ir.bam     + -5-ir.bam.bai"
  echo "     gatk-pr                  -5-ir.bam     + -5-ir.bam.bai  + -5-ir.bam.table"
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
sample=$1
functn=$2
tmp=temp-$sample

##
# Most of the variables, as well as the finishCheck function, are saved in a config file
#  in particular, the config file now saves the number of threads (num_thr) and max_ram settings
source ngs-config-bbox.sh

# ## Starting timestamp
echo 'NGS-TIMESTAMP' `date`

main () # called from the bottom so that the other functions are known
{
	# takes the chosen function and actually calls that function, nice and simple
	case "$functn" in
		bwa-mem)								bwa-mem ;;
		samtools-fixmate)				samtools-fixmate ;;
		sambamba-sort)					sambamba-sort ;;
		picard-mark-duplicates)	picard-mark-duplicates ;;
		samtools-index)					samtools-index ;;
		gatk-rtc)								gatk-rtc ;;
		gatk-ir)								gatk-ir ;;
		gatk-br)								gatk-br ;;
		gatk-pr)								gatk-pr ;;
		gatk-callable)					gatk-callable ;;
		gatk-hc)								gatk-hc ;;
		gatk-snp-select)				gatk-snp-select ;;
		gatk-snp-filt)					gatk-snp-filt ;;
		gatk-ind-select)				gatk-ind-select ;;
		gatk-ind-filt)					gatk-ind-filt ;;
		gatk-variant-combine)		gatk-variant-combine ;;
		gatk-variant-pass)			gatk-variant-pass ;;
		gatk-variant-eval)			gatk-variant-eval ;;
	esac
}

bwa-mem ()
{
	echo "NGS-CMD - Begin bwa mem mapping" # includes basic readgroup info from sample name	
	bwa mem -t $num_thr -M \
	  -R '@RG\tID:'$sample'\tSM:'$sample'\tPL:illumina\tPU:ns500\tLB:1234' \
	  $ref $sample-r1.fq.gz $sample-r2.fq.gz > $sample-1-raw.sam

	finishCheck $? "bwa mem mapping"

	samtools-fixmate # call the next tool
}

samtools-fixmate ()
{
	echo "NGS-CMD - Begin samtools fixmate" # fixes matepair info and converts sam to bam
	samtools fixmate $sample-1-raw.sam -O bam $sample-2-fm.bam

	finishCheck $? "samtools fixmate"

	## now we may delete the intermediate sam file
	rm $sample-1-raw.sam

	sambamba-sort # call the next tool
}

sambamba-sort ()
{
	echo "NGS-CMD - Begin sambamba sort" # creates many many temporary files
	sambamba sort -m $max_ram -t $num_thr -p -o $sample-3-s.bam $sample-2-fm.bam

	finishCheck $? "sambamba sort"

	## delete more intermediate files
	rm $sample-2-fm.bam

	picard-mark-duplicates # call the next tool
}

picard-mark-duplicates ()
{
	echo "NGS-CMD - Picard Mark Duplicates" # Also returns QC metrics
	java -jar $MarkDup \
	 TMP_DIR=$temp \
	 VALIDATION_STRINGENCY=LENIENT \
	 INPUT=$sample-3-s.bam \
	 METRICS_FILE=qc-$sample-4-md.metrics \
	 REMOVE_DUPLICATES=true \
	 OUTPUT=$sample-4-md.bam

	finishCheck $? "Picard Mark Duplicates"

	# delete intermediate files and make a copy of the QC metrics
	#rm $sample-3-s.bam
	#rm $sample-3-s.bam.bai

	samtools-index # call the next tool
}

samtools-index ()
{
	echo "NGS-CMD - Begin samtools indexing"
	samtools index $sample-4-md.bam

	finishCheck $? "samtools indexing"

	gatk-rtc # call the next tool
}

gatk-rtc ()
{
	echo "NGS-CMD - Start GATK RealignerTargetCreator" # known should include known indels
	java -Xmx$max_ram -jar $gatk \
	  -T RealignerTargetCreator -nt $num_thr \
	  -R $ref \
	  --known $dbsnp_vcf \
	  --known $indel_vcf \
	  -I $sample-4-md.bam \
	  -o $sample-4.intervals

	finishCheck $? "GATK RealignerTargetCreator"

	gatk-ir # call the next tool
}

gatk-ir ()
{
	echo "NGS-CMD - Start GATK IndelRealigner" # knownAlleles should include known indels
	java -Xmx$max_ram -jar $gatk \
	  -T IndelRealigner \
	  -R $ref \
	  --knownAlleles $dbsnp_vcf \
	  --knownAlleles $indel_vcf \
	  -filterNoBases \
	  -filterMBQ \
	  -filterRNC \
	  -targetIntervals $sample-4.intervals \
	  -I $sample-4-md.bam \
	  -o $sample-5-ir.bam

	finishCheck $? "GATK IndelRealigner"

	rm $sample-4-md.bam
	rm $sample-4-md.bam.bai

	gatk-br # call the next tool
}

gatk-br ()
{
	echo "NGS-CMD - Start GATK BaseRecalibrator" # knownSites should include any known variants
	java -Xmx$max_ram -jar $gatk \
	  -T BaseRecalibrator -nct $num_thr \
	  -R $ref \
	  --knownSites $dbsnp_vcf \
	  --knownSites $indel_vcf \
	  -I $sample-5-ir.bam \
	  -o $sample-5-ir.table

	finishCheck $? "GATK BaseRecalibrator"

	gatk-pr # call the next tool
}

gatk-pr ()
{
	echo "NGS-CMD - Start GATK PrintReads" # uses info from the previous steps
	java -Xmx$max_ram -jar $gatk \
	  -T PrintReads -nct $num_thr \
	  -R $ref \
	  --BQSR $sample-5-ir.table \
	  -I $sample-5-ir.bam \
	  -o $sample-ready.bam

	finishCheck $? "GATK PrintReads"

	rm $sample-5-ir.bam
	rm $sample-5-ir.bai
	rm $sample-4.intervals
	rm $sample-5-ir.table

	gatk-callable # call the next tool
}

gatk-callable ()
{
	echo "NGS-CMD - Start GATK Callable Loci tool" # can we call certain loci
	echo $clia_list_bed
	java -jar $gatk \
	  -T CallableLoci \
	  -R $ref \
	  -I $sample-ready.bam \
	  -L $care_list_bed \
	  -summary qc-$sample-5-call.sum \
	  -o qc-$sample-call.bed

	finishCheck $? "GATK Callable Loci AveraList"

	rm qc-$sample-call.bed

	gatk-hc # call the next tool
}

gatk-hc ()
{
	echo "NGS-CMD - Start GATK HaplotypeCaller" # the genotyper
	java -Xmx$max_ram -jar $gatk \
	 -T HaplotypeCaller -nct $num_thr \
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
	 -T VariantEval -nt $num_thr \
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

nt  = --num_threads - the number of data threads sent to the processor
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