// Bpipe pipeline
// Harriet Dashnow 22 March 2015
// Takes mapped exome bams. Simulates pooled exomes.


BASE="/group/bioi1/harrietd/pooled-parent/pooled_simulation"

// Tools
TOOLS="$BASE/tools"
PICARD="/usr/local/installed/picard/2.0.1/picard.jar"
GATK="/usr/local/installed/gatk/3.8"

// Scripts
MEDIAN_AWK="/group/bioi1/harrietd/git/micro-genotyper-long/repeat_genotyper_bpipes/median.awk"

// Ref files
EXCLUDE="$BASE/CS_excluded.bed"
REFBASE="/group/bioi1/shared/genomes/hg19/gatk" 
REF="$REFBASE/gatk.ucsc.hg19.fasta"
DBSNP="$REFBASE/dbsnp_138.hg19.vcf"
GOLD_STANDARD_INDELS="$REFBASE/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
GNOMAD="/group/bioi1/shared/genomes/hg19/gnomadData/gnomad_r2.0.2"

VCFANNO_CONFIG="/group/bioi1/harrietd/git/pooled_simulation/gnomadExomeGenomeSHORT_self.toml"

combined_bed="combined_target.bed"

PLATFORM='illumina'
threads=8 //for BWA
