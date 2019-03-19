// Bpipe pipeline
// Harriet Dashnow 22 March 2015
// Takes mapped exome bams. Simulates pooled exomes.


BASE="/group/bioi1/harrietd/pooled-parent/pooled_simulation"

// Set a good location for storing large temp files here
TMPDIR="$BASE/tmpdata"

// Tools
TOOLS="$BASE/tools"
PICARD="/group/bioi1/harrietd/src/picard-2.18.11/picard.jar"
GATK="/group/bioi1/harrietd/src/gatk-4.0.11.0/gatk"

// Scripts
MEDIAN_AWK="/group/bioi1/harrietd/git/micro-genotyper-long/repeat_genotyper_bpipes/median.awk"

// Ref files
EXCLUDE="$BASE/CS_excluded.bed"
REFBASE="/group/bioi1/shared/genomes/hg19/gatk" 
REF="$REFBASE/gatk.ucsc.hg19.fasta"
DBSNP="$REFBASE/dbsnp_138.hg19.vcf"
GOLD_STANDARD_INDELS="$REFBASE/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
GNOMAD="/group/bioi1/shared/genomes/hg19/gnomadData/gnomad_r2.0.2"
VEP_SYN="/group/bioi1/shared/genomes/GRCh37/vep_chr_synonyms.txt"

// Configuration files
VCFANNO_CONFIG="/group/bioi1/harrietd/git/pooled_simulation/gnomadExomeGenomeSHORT_self.toml"

// Bed files
combined_bed="combined_target.bed"
INT_BED="/group/bioi1/harrietd/pooled-parent/filter/intersected/intersection_Nextera_SureSelectClinical.bed"

PLATFORM='illumina'
//threads=8 //for BWA
