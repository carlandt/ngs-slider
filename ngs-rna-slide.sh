#!/bin/bash
####
##  The complete NGS pipeline for RNA sequence to variants with haplotypecaller and QC metrics
##   
##  Written by Tristan M. Carland, PhD
####
if [ $# -lt 3 ]; then
  echo ""
  echo " Tristan's NGS Scripts - Full NGS Slide v1.0 (GATK 3.3-0)"
  echo ""
  echo "  usage: `basename $0` <config> <sample> <starting-sub-routine>"
  echo ""
  echo "   Functions                Expected <sample=S>"
  echo "     star-index-ref           (ref genome)"
  echo "     star-map-first-pass      S-r1.fq.gz   + S-r2.fq.gz      + temp-ref1-S"
  echo "     star-index-sample        temp-map1-S  + temp-ref1-S"
  echo "     star-map-second-pass     temp-ref2-S"
  echo "     sam-dup-mate-sort        temp-map1-S/Aligned.out.sam"
  echo "     samtools-index           -s.bam"
  echo "     gatk-rtc                 -s.bam       + -s.bam.bai"
  echo "     gatk-ir                  -s.bam       + -s.bam.bai      + -s.intervals"
  echo "     gatk-br                  -ir.bam      + -ir.bam.bai"
  echo "     gatk-pr                  -ir.bam      + -ir.bam.bai     + -ir.bam.table"
  echo "     gatk-callable            -ready.bam   + -ready.bam.bai"
  echo "     gatk-hc                  -ready.bam   + -ready.bam.bai"
  echo "     gatk-filt                -xxx.vcf"
  echo "     gatk-variant-eval        -xx.vcf"
  echo ""
  exit 1
fi

##
# Parsing of command line input variables and various bundled library files
config=$1
sample=$2
functn=$3
tmp=temp-$sample

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
		star-index-ref)			star-index-ref ;;
		star-map-first-pass)	star-map-first-pass ;;
		star-index-sample)		star-index-sample ;;
		star-map-second-pass)	star-map-second-pass ;;
		sam-dup-mate-sort)		sam-dup-mate-sort ;;
		samtools-index)			samtools-index ;;
		gatk-rtc)				gatk-rtc ;;
		gatk-ir)				gatk-ir ;;
		gatk-br)				gatk-br ;;
		gatk-pr)				gatk-pr ;;
		gatk-callable)			gatk-callable ;;
		gatk-hc)				gatk-hc ;;
		gatk-snp-select)		gatk-snp-select ;;
		gatk-snp-filt)			gatk-snp-filt ;;
		gatk-ind-select)		gatk-ind-select ;;
		gatk-ind-filt)			gatk-ind-filt ;;
		gatk-variant-combine)	gatk-variant-combine ;;
		gatk-variant-pass)		gatk-variant-pass ;;
		gatk-variant-eval)		gatk-variant-eval ;;
	esac
}

##
# 
star-index-ref ()
{
	echo "NGS-CMD - Begin STAR index of genome" # 
	mkdir -p temp-ref1-$sample

	STAR --runThreadN 10 \
	--runMode genomeGenerate \
	--genomeDir temp-ref1-$sample \
	--genomeFastaFiles $ref \
	--sjdbGTFfile $refgtf \
	--sjdbOverhang 100

	finishCheck $? "STAR genome index"

	star-map-first-pass # call the next tool
}

##
#
star-map-first-pass ()
{
	echo "NGS-CMD - Begin STAR map 1pass" # 

	mkdir -p temp-map1-$sample
	cd temp-map1-$sample

	echo 'pwd   :'`pwd`
	echo 'ls    :'`ls`
	echo 'ls ../:'`ls ../`

	STAR --runThreadN 10 \
	--genomeDir ../temp-ref1-$sample \
	--readFilesIn ../$sample-r1.fq.gz ../$sample-r2.fq.gz \
	--readFilesCommand gzip -dc

	finishCheck $? "STAR map 1pass"

	cd ..

	# call the next tool
	star-index-sample
}

star-index-sample ()
{
	echo "NGS-CMD - Begin STAR index of sample" # 
	mkdir -p temp-ref2-$sample

	STAR --runThreadN 10 \
	--runMode genomeGenerate \
	--genomeDir temp-ref2-$sample \
	--genomeFastaFiles $ref \
	--sjdbFileChrStartEnd temp-map1-$sample/SJ.out.tab \
	--sjdbOverhang 75 \
	--sjdbGTFfile $refgtf

	finishCheck $? "STAR sample index"

	# call the next tool
	star-map-second-pass
}
#
# --outSAMtype BAM SortedByCoordinate \
#
star-map-second-pass ()
{
	echo "NGS-CMD - Begin STAR map 2pass" # 
	mkdir -p temp-map2-$sample
	cd temp-map2-$sample

	STAR --runThreadN 10 \
	--genomeDir ../temp-ref2-$sample \
	--readFilesIn ../$sample-r1.fq.gz ../$sample-r2.fq.gz \
	--readFilesCommand gzip -dc \
	--outSAMstrandField intronMotif \
	--outFilterIntronMotifs RemoveNoncanonical \
	--outSAMattrRGline ID:$sample SM:$sample PL:illumina PU:hs LB:1234 \
	--outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped Fastx

	finishCheck $? "STAR map 2pass"

	# call the next tool
	picard-markdup
}

sam-dup-mate-sort ()
{
	tool="Samblaster MarkDup + Samtools Fixmate + Sambamba Sort"
	echo "NGS-CMD - Begin $tool"
	
	samblaster -r -i temp-map2-$sample/Aligned.out.sam \
	| samtools fixmate - -O bam \
	| sambamba sort -m $job_ram -t $job_thr -p -o $sample-md.bam -

	finishCheck $? $tool

	# call the next tool
	samtools-index
}

picard-markdup ()
{
	tool="Picard MarkDup"
	echo "NGS-CMD - Begin $tool"

	java -jar $MarkDup \
	I=temp-map2-$sample/Aligned.out.sam \
	O=$sample-md.bam \
	M=$sample-md.metrics \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT

	finishCheck $? $tool

	gatk-sncr # call the next tool
}

gatk-sncr ()
{
	echo "NGS-CMD - Start GATK SplitNCigarReads" # 
	java -Xmx$job_ram -jar $gatk \
	  -T SplitNCigarReads \
	  -R $ref \
	  -I $sample-md.bam \
	  -rf ReassignOneMappingQuality \
	  -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS \
	  -o $sample-md-c.bam

	finishCheck $? "GATK SplitNCigarReads"

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
	  -I $sample-md-c.bam \
	  -o $sample-md-c.intervals

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
	  -targetIntervals $sample-md-c.intervals \
	  -I $sample-md-c.bam \
	  -o $sample-ir.bam

	finishCheck $? "GATK IndelRealigner"

	rm $sample-md-c.bam
	rm $sample-md-c.bam.bai

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
	rm $sample-md-c.intervals
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
	 -dontUseSoftClippedBases \
	 -stand_call_conf 20 \
	 -stand_emit_conf 20 \
	 -o $sample-raw.vcf

	finishCheck $? "GATK HaplotypeCaller"

	gatk-filt # call the next tool
}

gatk-filt ()
{
	echo "NGS-CMD - Begin variant filtration"
	java -Xmx$job_ram -jar $gatk \
	 -T VariantFiltration \
	 -R $ref \
	 -V $sample-raw.vcf \
	 -window 35 -cluster 3 \
	 -filterName FS -filter "FS > 30.0" \
	 -filterName QD -filter "QD < 2.0" \
	 -o $sample-mark.vcf

	finishCheck $? "variant filtration"

	gatk-variant-eval # call the next tool
}

# gatk-variant-pass ()
# {
# 	echo "NGS-CMD - Start PASS Collection" # Apply a simple PASS/FAIL to get the final set
# 	java -Xmx2g -jar $gatk \
# 	 -T SelectVariants \
# 	 -R $ref \
# 	 -V $sample-mark.vcf \
# 	 --excludeFiltered \
# 	 --excludeNonVariants \
# 	 -o $sample-hard.vcf

# 	finishCheck $? "PASS Collection"

# 	echo " To clean up all intermediate files, just delete the temp-$2-filt directory"

# 	gatk-variant-eval # call the next tool
# }

gatk-variant-eval ()
{
	echo "NGS-CMD - Variant Evaluation"
	java -Xmx$job_ram -jar $gatk \
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
