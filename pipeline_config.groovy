// Bpipe pipeline
// Harriet Dashnow 22 March 2015
// Takes mapped exome bams. Simulates pooled exomes.


BASE="/group/bioi1/harrietd/pooled-parent/pooled_simulation"

// Tools
TOOLS="$BASE/tools"
PICARD="/usr/local/installed/picard/2.0.1/picard.jar"
GATK="/usr/local/installed/gatk/3.8"

// Ref files
EXOME_TARGET="/group/bioi1/shared/genomes/hg19/exome_targets/nextera_rapid_capture_exome/target_region.bed"
EXCLUDE="$BASE/CS_excluded.bed"
REFBASE="/group/bioi1/shared/genomes/hg19/gatk" 
REF="$REFBASE/gatk.ucsc.hg19.fasta"
DBSNP="$REFBASE/dbsnp_138.hg19.vcf"
GOLD_STANDARD_INDELS="$REFBASE/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"

combined_bed="combined_target.bed"
