#!/bin/bash
####
##  Wrapper script for running MuTect, conforms to ngs- package design
##  Written by Tristan M. Carland, PhD
####
if [ $# -lt 4 ]; then
  echo ""
  echo " Tristan's NGS Scripts - MuTect (MuTect 1.1.4, built upon GATK v2.2-25)"
  echo ""
  echo "  usage: `basename $0` <config> <tumor-sample> <normal-sample> <out-base>"
  echo ""

  exit 1
fi

##
# Parsing of command line input variables and various bundled library files
source $1 # the config file in $1 contains the location of mutect and reference files
tumor_bam=$2
normy_bam=$3
out_base=$4

# ## Starting timestamp
echo 'NGS-TIMESTAMP' `date`

##
# MuTect! - see http://www.broadinstitute.org/cancer/cga/mutect_run
#   --intervals <intervals_to_process>   --coverage_file <coverage.wig.txt>
#  - requires java 1.6 for some damned reason...
mutect ()
{
  echo "NGS-CMD - Begin MuTect 1.1.4"
  /usr/lib/jvm/java-6-openjdk-amd64/jre/bin/java -Xmx$job_ram -jar $mutect \
  --analysis_type MuTect \
  --reference_sequence $ref \
  --cosmic $cosmic \
  --dbsnp $dbsnp_vcf \
  --input_file:tumor $tumor_bam \
  --input_file:normal $normy_bam \
  --out $out_base-raw.csv \
  --vcf $out_base-raw.vcf

  finishCheck $? "bwa mem pipe"
}

echo 'NGS-TIMESTAMP' `date`

mutect # call the mutect function - this way it's saved as a function we can move

exit 0 # end shell scripts with an exit for easier commenting at the end 
###########################################################################
# Notes #
#########


---------------------------------------------------------------------------------
The Genome Analysis Toolkit (GATK) v2.2-25-g2a68eab, Compiled 2012/11/08 10:30:02
Copyright (c) 2010 The Broad Institute
For support and documentation go to http://www.broadinstitute.org/gatk
---------------------------------------------------------------------------------
---------------------------------------------------------------------------------
usage: java -jar muTect-1.1.4.jar -T <analysis_type> [-args <arg_file>] [-I <input_file>] [-rbs <read_buffer_size>] [-et 
       <phone_home>] [-K <gatk_key>] [-tag <tag>] [-rf <read_filter>] [-L <intervals>] [-XL <excludeIntervals>] [-isr 
       <interval_set_rule>] [-im <interval_merging>] [-ip <interval_padding>] [-R <reference_sequence>] [-ndrs] 
       [--disableRandomization] [-maxRuntime <maxRuntime>] [-maxRuntimeUnits <maxRuntimeUnits>] [-dt <downsampling_type>] 
       [-dfrac <downsample_to_fraction>] [-dcov <downsample_to_coverage>] [-baq <baq>] [-baqGOP <baqGapOpenPenalty>] [-PF 
       <performanceLog>] [-OQ] [-BQSR <BQSR>] [-DIQ] [-EOQ] [-preserveQ <preserve_qscores_less_than>] [-DBQ 
       <defaultBaseQualities>] [-S <validation_strictness>] [-rpr] [-kpr] [-U <unsafe>] [-nt <num_threads>] [-nct 
       <num_cpu_threads_per_data_thread>] [-mte] [-bfh <num_bam_file_handles>] [-rgbl <read_group_black_list>] [-ped 
       <pedigree>] [-pedString <pedigreeString>] [-pedValidationType <pedigreeValidationType>] [-l <logging_level>] [-log 
       <log_to_file>] [-h] [-filterMBQ] [--noop] [--tumor_sample_name <tumor_sample_name>] [--bam_tumor_sample_name 
       <bam_tumor_sample_name>] [--normal_sample_name <normal_sample_name>] [--force_output] [--force_alleles] 
       [--only_passing_calls] [--initial_tumor_lod <initial_tumor_lod>] [--tumor_lod <tumor_lod>] [--fraction_contamination 
       <fraction_contamination>] [--minimum_mutation_cell_fraction <minimum_mutation_cell_fraction>] [--normal_lod 
       <normal_lod>] [--dbsnp_normal_lod <dbsnp_normal_lod>] [--somatic_classification_normal_power_threshold 
       <somatic_classification_normal_power_threshold>] [--minimum_normal_allele_fraction <minimum_normal_allele_fraction>] 
       [--tumor_f_pretest <tumor_f_pretest>] [--min_qscore <min_qscore>] [--gap_events_threshold <gap_events_threshold>] 
       [--heavily_clipped_read_fraction <heavily_clipped_read_fraction>] [--clipping_bias_pvalue_threshold 
       <clipping_bias_pvalue_threshold>] [--fraction_mapq0_threshold <fraction_mapq0_threshold>] [--pir_median_threshold 
       <pir_median_threshold>] [--pir_mad_threshold <pir_mad_threshold>] [--required_maximum_alt_allele_mapping_quality_score 
       <required_maximum_alt_allele_mapping_quality_score>] [--max_alt_alleles_in_normal_count 
       <max_alt_alleles_in_normal_count>] [--max_alt_alleles_in_normal_qscore_sum <max_alt_alleles_in_normal_qscore_sum>] 
       [--max_alt_allele_in_normal_fraction <max_alt_allele_in_normal_fraction>] [--power_constant_qscore 
       <power_constant_qscore>] [--absolute_copy_number_data <absolute_copy_number_data>] [--power_constant_af 
       <power_constant_af>] [-o <out>] [-vcf <vcf>] [-dbsnp <dbsnp>] [-cosmic <cosmic>] [-cov <coverage_file>] [-cov_q20 
       <coverage_20_q20_file>] [-pow <power_file>] [-tdf <tumor_depth_file>] [-ndf <normal_depth_file>]

 -T,--analysis_type <analysis_type>                                         Type of analysis to run
 -args,--arg_file <arg_file>                                                Reads arguments from the specified file
 -I,--input_file <input_file>                                               SAM or BAM file(s)
 -rf,--read_filter <read_filter>                                            Specify filtration criteria to apply to each 
                                                                            read individually
 -L,--intervals <intervals>                                                 One or more genomic intervals over which to 
                                                                            operate. Can be explicitly specified on the 
                                                                            command line or in a file (including a rod 
                                                                            file)
 -XL,--excludeIntervals <excludeIntervals>                                  One or more genomic intervals to exclude 
                                                                            from processing. Can be explicitly specified 
                                                                            on the command line or in a file (including 
                                                                            a rod file)
 -isr,--interval_set_rule <interval_set_rule>                               Indicates the set merging approach the 
                                                                            interval parser should use to combine the 
                                                                            various -L or -XL inputs (UNION|
                                                                            INTERSECTION)
 -im,--interval_merging <interval_merging>                                  Indicates the interval merging rule we 
                                                                            should use for abutting intervals (ALL|
                                                                            OVERLAPPING_ONLY)
 -ip,--interval_padding <interval_padding>                                  Indicates how many basepairs of padding to 
                                                                            include around each of the intervals 
                                                                            specified with the -L/--intervals argument
 -R,--reference_sequence <reference_sequence>                               Reference sequence file
 -ndrs,--nonDeterministicRandomSeed                                         Makes the GATK behave non deterministically, 
                                                                            that is, the random numbers generated will 
                                                                            be different in every run
 --disableRandomization                                                     Completely eliminates randomization from 
                                                                            nondeterministic methods. To be used mostly 
                                                                            in the testing framework where dynamic 
                                                                            parallelism can result in differing numbers 
                                                                            of calls to the generator.
 -dt,--downsampling_type <downsampling_type>                                Type of reads downsampling to employ at a 
                                                                            given locus.  Reads will be selected 
                                                                            randomly to be removed from the pile based 
                                                                            on the method described here (NONE|ALL_READS|
                                                                            BY_SAMPLE)
 -dfrac,--downsample_to_fraction <downsample_to_fraction>                   Fraction [0.0-1.0] of reads to downsample to
 -dcov,--downsample_to_coverage <downsample_to_coverage>                    Coverage [integer] to downsample to at any 
                                                                            given locus; note that downsampled reads are 
                                                                            randomly selected from all possible reads at 
                                                                            a locus
 -baq,--baq <baq>                                                           Type of BAQ calculation to apply in the 
                                                                            engine (OFF|CALCULATE_AS_NECESSARY|
                                                                            RECALCULATE)
 -baqGOP,--baqGapOpenPenalty <baqGapOpenPenalty>                            BAQ gap open penalty (Phred Scaled). 
                                                                             Default value is 40.  30 is perhaps better 
                                                                            for whole genome call sets
 -PF,--performanceLog <performanceLog>                                      If provided, a GATK runtime performance log 
                                                                            will be written to this file
 -OQ,--useOriginalQualities                                                 If set, use the original base quality scores 
                                                                            from the OQ tag when present instead of the 
                                                                            standard scores
 -BQSR,--BQSR <BQSR>                                                        The input covariates table file which 
                                                                            enables on-the-fly base quality score 
                                                                            recalibration
 -DIQ,--disable_indel_quals                                                 If true, disables printing of base insertion 
                                                                            and base deletion tags (with -BQSR)
 -EOQ,--emit_original_quals                                                 If true, enables printing of the OQ tag with 
                                                                            the original base qualities (with -BQSR)
 -preserveQ,--preserve_qscores_less_than <preserve_qscores_less_than>       Bases with quality scores less than this 
                                                                            threshold won't be recalibrated (with -BQSR)
 -DBQ,--defaultBaseQualities <defaultBaseQualities>                         If reads are missing some or all base 
                                                                            quality scores, this value will be used for 
                                                                            all base quality scores
 -S,--validation_strictness <validation_strictness>                         How strict should we be with validation 
                                                                            (STRICT|LENIENT|SILENT)
 -rpr,--remove_program_records                                              Should we override the Walker's default and 
                                                                            remove program records from the SAM header
 -kpr,--keep_program_records                                                Should we override the Walker's default and 
                                                                            keep program records from the SAM header
 -U,--unsafe <unsafe>                                                       If set, enables unsafe operations: nothing 
                                                                            will be checked at runtime.  For expert 
                                                                            users only who know what they are doing.  We 
                                                                            do not support usage of this argument. 
                                                                            (ALLOW_UNINDEXED_BAM|
                                                                            ALLOW_UNSET_BAM_SORT_ORDER|
                                                                            NO_READ_ORDER_VERIFICATION|
                                                                            ALLOW_SEQ_DICT_INCOMPATIBILITY|
                                                                            LENIENT_VCF_PROCESSING|ALL)
 -nt,--num_threads <num_threads>                                            How many data threads should be allocated to 
                                                                            running this analysis.
 -nct,--num_cpu_threads_per_data_thread <num_cpu_threads_per_data_thread>   How many CPU threads should be allocated per 
                                                                            data thread to running this analysis?
 -mte,--monitorThreadEfficiency                                             Enable GATK threading efficiency monitoring
 -bfh,--num_bam_file_handles <num_bam_file_handles>                         The total number of BAM file handles to keep 
                                                                            open simultaneously
 -rgbl,--read_group_black_list <read_group_black_list>                      Filters out read groups matching 
                                                                            <TAG>:<STRING> or a .txt file containing the 
                                                                            filter strings one per line.
 -ped,--pedigree <pedigree>                                                 Pedigree files for samples
 -pedString,--pedigreeString <pedigreeString>                               Pedigree string for samples
 -pedValidationType,--pedigreeValidationType <pedigreeValidationType>       How strict should we be in validating the 
                                                                            pedigree information? (STRICT|SILENT)
 -l,--logging_level <logging_level>                                         Set the minimum level of logging, i.e. 
                                                                            setting INFO get's you INFO up to FATAL, 
                                                                            setting ERROR gets you ERROR and FATAL level 
                                                                            logging.
 -log,--log_to_file <log_to_file>                                           Set the logging location
 -h,--help                                                                  Generate this help message

Arguments for MalformedReadFilter:
 -filterMBQ,--filter_mismatching_base_and_quals   if a read has mismatching number of bases and base qualities, filter 
                                                  out the read instead of blowing up.

Arguments for MuTect:
 --noop                                                                                   used for debugging, basically 
                                                                                          exit as soon as we get the 
                                                                                          reads
 --tumor_sample_name <tumor_sample_name>                                                  name to use for tumor in 
                                                                                          output files
 --bam_tumor_sample_name <bam_tumor_sample_name>                                          if the tumor bam contains 
                                                                                          multiple samples, only use 
                                                                                          read groups with SM equal to 
                                                                                          this value
 --normal_sample_name <normal_sample_name>                                                name to use for normal in 
                                                                                          output files
 --force_output                                                                           force output for each site
 --force_alleles                                                                          force output for all alleles 
                                                                                          at each site
 --only_passing_calls                                                                     only emit passing calls
 --initial_tumor_lod <initial_tumor_lod>                                                  Initial LOD threshold for 
                                                                                          calling tumor variant
 --tumor_lod <tumor_lod>                                                                  LOD threshold for calling 
                                                                                          tumor variant
 --fraction_contamination <fraction_contamination>                                        estimate of fraction (0-1) of 
                                                                                          physical contamination with 
                                                                                          other unrelated samples
 --minimum_mutation_cell_fraction <minimum_mutation_cell_fraction>                        minimum fraction of cells 
                                                                                          which are presumed to have a 
                                                                                          mutation, used to handle 
                                                                                          non-clonality and 
                                                                                          contamination
 --normal_lod <normal_lod>                                                                LOD threshold for calling 
                                                                                          normal non-germline
 --dbsnp_normal_lod <dbsnp_normal_lod>                                                    LOD threshold for calling 
                                                                                          normal non-variant at dbsnp 
                                                                                          sites
 --somatic_classification_normal_power_threshold                                          Power threshold for normal to 
<somatic_classification_normal_power_threshold>                                           determine germline vs variant
 --minimum_normal_allele_fraction <minimum_normal_allele_fraction>                        minimum allele fraction to be 
                                                                                          considered in normal, useful 
                                                                                          for normal sample contaminated 
                                                                                          with tumor
 --tumor_f_pretest <tumor_f_pretest>                                                      for computational efficiency, 
                                                                                          reject sites with allelic 
                                                                                          fraction below this threshold
 --min_qscore <min_qscore>                                                                threshold for minimum base 
                                                                                          quality score
 --gap_events_threshold <gap_events_threshold>                                            how many gapped events 
                                                                                          (ins/del) are allowed in 
                                                                                          proximity to this candidate
 --heavily_clipped_read_fraction <heavily_clipped_read_fraction>                          if this fraction or more of 
                                                                                          the bases in a read are 
                                                                                          soft/hard clipped, do not use 
                                                                                          this read for mutation calling
 --clipping_bias_pvalue_threshold <clipping_bias_pvalue_threshold>                        pvalue threshold for fishers 
                                                                                          exact test of clipping bias in 
                                                                                          mutant reads vs ref reads
 --fraction_mapq0_threshold <fraction_mapq0_threshold>                                    threshold for determining if 
                                                                                          there is relatedness between 
                                                                                          the alt and ref allele read 
                                                                                          piles
 --pir_median_threshold <pir_median_threshold>                                            threshold for clustered read 
                                                                                          position artifact median
 --pir_mad_threshold <pir_mad_threshold>                                                  threshold for clustered read 
                                                                                          position artifact MAD
 --required_maximum_alt_allele_mapping_quality_score                                      required minimum value for 
<required_maximum_alt_allele_mapping_quality_score>                                       tumor alt allele maximum 
                                                                                          mapping quality score
 --max_alt_alleles_in_normal_count <max_alt_alleles_in_normal_count>                      threshold for maximum 
                                                                                          alternate allele counts in 
                                                                                          normal
 --max_alt_alleles_in_normal_qscore_sum <max_alt_alleles_in_normal_qscore_sum>            threshold for maximum 
                                                                                          alternate allele quality score 
                                                                                          sum in normal
 --max_alt_allele_in_normal_fraction <max_alt_allele_in_normal_fraction>                  threshold for maximum 
                                                                                          alternate allele fraction in 
                                                                                          normal
 --power_constant_qscore <power_constant_qscore>                                          Phred scale quality score 
                                                                                          constant to use in power 
                                                                                          calculations
 --absolute_copy_number_data <absolute_copy_number_data>                                  Absolute Copy Number Data, as 
                                                                                          defined by Absolute, to use in 
                                                                                          power calculations
 --power_constant_af <power_constant_af>                                                  Allelic fraction constant to 
                                                                                          use in power calculations
 -o,--out <out>                                                                           Call-stats output
 -vcf,--vcf <vcf>                                                                         VCF output of mutation 
                                                                                          candidates
 -dbsnp,--dbsnp <dbsnp>                                                                   VCF file of DBSNP information
 -cosmic,--cosmic <cosmic>                                                                VCF file of COSMIC sites
 -cov,--coverage_file <coverage_file>                                                     write out coverage in WIGGLE 
                                                                                          format to this file
 -cov_q20,--coverage_20_q20_file <coverage_20_q20_file>                                   write out 20x of Q20 coverage 
                                                                                          in WIGGLE format to this file
 -pow,--power_file <power_file>                                                           write out power in WIGGLE 
                                                                                          format to this file
 -tdf,--tumor_depth_file <tumor_depth_file>                                               write out tumor read depth in 
                                                                                          WIGGLE format to this file
 -ndf,--normal_depth_file <normal_depth_file>                                             write out normal read depth in 
                                                                                          WIGGLE format to this file

Available Reference Ordered Data types:
         Name       FeatureType   Documentation
         BCF2    VariantContext   http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_utils_codecs_bcf2_BCF2Codec.html
          BED        BEDFeature   http://www.broadinstitute.org/gatk/gatkdocs/org_broad_tribble_bed_BEDCodec.html
EXAMPLEBINARY           Feature   http://www.broadinstitute.org/gatk/gatkdocs/org_broad_tribble_example_ExampleBinaryCodec.html
     GELITEXT   GeliTextFeature   http://www.broadinstitute.org/gatk/gatkdocs/org_broad_tribble_gelitext_GeliTextCodec.html
     OLDDBSNP   OldDbSNPFeature   http://www.broadinstitute.org/gatk/gatkdocs/org_broad_tribble_dbsnp_OldDbSNPCodec.html
          VCF    VariantContext   http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_utils_codecs_vcf_VCFCodec.html

For a full description of this walker, see its GATKdocs at:
http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_cga_tools_gatk_walkers_cancer_mutect_MuTect.html

