library('ggplot2')
library('dplyr')

# Read in recall/variant recovery data
HC_vars = read.csv("data/pooled_sim_variants.csv")
#FB_vars = read.csv("data/pooled_sim_variants_freebayes.csv")
FB_vars = read.csv("data/random_pools2-freebayes-pooled_sim_compare.csv")
#FBqual_vars = read.csv("data/pooled_sim_variants_freebayes_highqual.csv")
FBqual_vars = read.csv("data/random_pools1-freebayes_qual-pooled_sim_compare.csv")
HCqual_vars = read.csv("data/pooled_sim_variants_haplotypecaller_qual.csv")
# FBqual_ind_vars = read.csv("data/pooled_sim_variants_freebayes_highqual-ind.csv")

# group by pool and calcuate stats
HC_vars %>% group_by(pool) %>% 
  summarise(n_recovered = sum(recovered), total = length(recovered), 
            variantcaller = "GATK Haplotype Caller") -> HC_poolstats
FB_vars %>% group_by(pool) %>% 
  summarise(n_recovered = sum(recovered), total = length(recovered), 
            variantcaller = "FreeBayes") -> FB_poolstats
FBqual_vars %>% group_by(pool) %>% 
  summarise(n_recovered = sum(recovered), total = length(recovered), 
            variantcaller = "FreeBayes QUAL>20") -> FBqual_poolstats
HCqual_vars %>% group_by(pool) %>% 
  summarise(n_recovered = sum(recovered), total = length(recovered), 
            variantcaller = "GATK Haplotype Caller QUAL>20") -> HCqual_poolstats
# FBqual_ind_vars %>% group_by(pool) %>% 
#   summarise(n_recovered = sum(recovered), total = length(recovered), 
#             variantcaller = "FreeBayes individuals QUAL>20, pool unfiltered") -> FBqual_ind_poolstats

poolstats = rbind(HC_poolstats, FB_poolstats, FBqual_poolstats, HCqual_poolstats)#, FBqual_ind_poolstats)
poolstats$percent = 100*poolstats$n_recovered/poolstats$total


# Read in false positive data
HC_falsepos = read.csv("data/pooled_sim_variants_haplotypecaller_falsepos.csv")
FB_falsepos = read.csv("data/pooled_sim_variants_freebayes_falsepos.csv")
FBqual_falsepos = read.csv("data/pooled_sim_variants_freebayes_highqual_falsepos.csv")
HCqual_falsepos = read.csv("data/pooled_sim_variants_haplotypecaller_qual_falsepos.csv")
# FBqual_ind_falsepos = read.csv("data/pooled_sim_variants_freebayes_highqual-ind_falsepos.csv")

# group by pool and calcuate stats
HC_falsepos %>% group_by(pool) %>% 
  summarise(n_falsepos = sum(false_positive), total = length(false_positive),
            variantcaller = "GATK Haplotype Caller") -> HC_falsepos_stats
FB_falsepos %>% group_by(pool) %>% 
  summarise(n_falsepos = sum(false_positive), total = length(false_positive),
            variantcaller = "FreeBayes") -> FB_falsepos_stats
FBqual_falsepos %>% group_by(pool) %>% 
  summarise(n_falsepos = sum(false_positive), total = length(false_positive),
            variantcaller = "FreeBayes QUAL>20") -> FBqual_falsepos_stats
HCqual_falsepos %>% group_by(pool) %>% 
  summarise(n_falsepos = sum(false_positive), total = length(false_positive),
            variantcaller = "GATK Haplotype Caller QUAL>20") -> HCqual_falsepos_stats
# FBqual_ind_falsepos %>% group_by(pool) %>% 
#   summarise(n_falsepos = sum(false_positive), total = length(false_positive),
#             variantcaller = "FreeBayes individuals QUAL>20, pool unfiltered") -> FBqual_ind_falsepos_stats

falsepos_stats = rbind(HC_falsepos_stats, FB_falsepos_stats, FBqual_falsepos_stats, HCqual_falsepos_stats)#, FBqual_ind_falsepos_stats)
falsepos_stats$percent_fp = 100*falsepos_stats$n_falsepos/falsepos_stats$total

#Recall just for HaplotypeCaller
ggplot(data=poolstats[poolstats$variantcaller == "GATK Haplotype Caller",], 
       aes(x=pool, y=percent)) + 
  geom_line() + 
  geom_point() +
  labs(x="Number of samples in pool", y="Recall %")

# Plot recall and false positives as a percentage
ggplot(data=poolstats, aes(x=pool, y=percent, colour = variantcaller)) + 
  geom_line() + 
  geom_point() +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x="Number of samples in pool", y="Recall %") +
  geom_point(size=2) + theme_classic(base_size = 20)

ggplot(data=falsepos_stats, aes(x=pool, y=percent_fp, colour = variantcaller)) + 
  geom_line() + 
  geom_point() +
  labs(x="Number of samples in pool", y="False positives %") +
  geom_point(size=2) + theme_classic(base_size = 20)

# Plot absolute number of recalls and false positives
ggplot(data=poolstats, aes(x=pool, y=n_recovered, colour = variantcaller)) + 
  geom_line() + 
  geom_point() +
  labs(x="Number of samples in pool", y="Recall")

ggplot(data=falsepos_stats, aes(x=pool, y=n_falsepos, colour = variantcaller)) + 
  geom_line() + 
  geom_point() +
  labs(x="Number of samples in pool", y="False positives")


# Group by number of non-reference alleles at locus for recall
group_nonref_alleles = function(variants, variantcaller) {
  variants %>% group_by(pool, nonref_allele_count_truth) %>% 
    summarise(n_recovered = sum(recovered), total = length(recovered)) -> vars_alleles
  vars_alleles$percent_recovered = vars_alleles$n_recovered/vars_alleles$total*100
  vars_alleles$pool_id = paste0('Pool_', vars_alleles$pool)
  vars_alleles$variantcaller = variantcaller
  return(vars_alleles)
}

HC_vars_alleles = group_nonref_alleles(HC_vars, "GATK Haplotype Caller")
FB_vars_alleles = group_nonref_alleles(FB_vars, "Freebayes")
FBqual_vars_alleles = group_nonref_alleles(FBqual_vars, "FreeBayes QUAL>20")
HCqual_vars_alleles = group_nonref_alleles(HCqual_vars, "GATK Haplotype Caller QUAL>20")
#FBqual_ind_vars_alleles = group_nonref_alleles(FBqual_ind_vars, "FreeBayes individuals QUAL>20, pool unfiltered")

vars_alleles = rbind(HC_vars_alleles, FB_vars_alleles, FBqual_vars_alleles, 
                     HCqual_vars_alleles)#, FBqual_ind_vars_alleles)

ggplot(data=vars_alleles[vars_alleles$nonref_allele_count_truth >0,], aes(x=pool, y=percent_recovered, 
                              colour=factor(nonref_allele_count_truth))) + 
  geom_point() + geom_line() + facet_wrap(~variantcaller) +
  labs(x="Number of samples in pool", y="Recall %") +
  geom_point(size=2) + theme_classic(base_size = 16)


ggplot(data=vars_alleles[vars_alleles$nonref_allele_count_truth == 1,], 
       aes(x=pool, y=percent_recovered, colour = variantcaller)) + 
  geom_line() + 
  geom_point() +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x="Number of samples in pool", y="Recall %", 
       title = "Variant loci with one allele present in the pool") +
  geom_point(size=2) + theme_classic(base_size = 20)

ggplot(data=vars_alleles[vars_alleles$nonref_allele_count_truth == 1,], 
       aes(x=pool, y=n_recovered, colour = variantcaller)) + 
  geom_line() + 
  geom_point() +
  labs(x="Number of samples in pool", y="Recall", title = "Variant loci with one allele present in the pool")


# Group by number of non-reference alleles at locus for false positives
group_nonref_alleles_falsepos = function(variants, variantcaller) {
  variants %>% group_by(pool, nonref_allele_count_obs) %>% 
    summarise(n_falsepos = sum(false_positive), total = length(false_positive)) -> vars_alleles
  vars_alleles$percent_falsepos = vars_alleles$n_falsepos/vars_alleles$total*100
  vars_alleles$pool_id = paste0('Pool_', vars_alleles$pool)
  vars_alleles$variantcaller = variantcaller
  return(vars_alleles)
}

HC_falsepos_alleles = group_nonref_alleles_falsepos(HC_falsepos, "GATK Haplotype Caller")
FB_falsepos_alleles = group_nonref_alleles_falsepos(FB_falsepos, "Freebayes")
FBqual_falsepos_alleles = group_nonref_alleles_falsepos(FBqual_falsepos, "FreeBayes QUAL>20")
HCqual_falsepos_alleles = group_nonref_alleles_falsepos(HCqual_falsepos, "GATK Haplotype Caller QUAL>20")
#FBqual_ind_falsepos_alleles = group_nonref_alleles_falsepos(FBqual_ind_falsepos, "FreeBayes individuals QUAL>20, pool unfiltered")

falsepos_alleles = rbind(HC_falsepos_alleles, FB_falsepos_alleles, FBqual_falsepos_alleles, 
                     HCqual_falsepos_alleles)#, FBqual_ind_falsepos_alleles)

ggplot(data=falsepos_alleles, aes(x=nonref_allele_count_obs, y=percent_falsepos, 
                              colour=pool_id)) + 
  geom_point() + geom_line() + facet_wrap(~variantcaller)

ggplot(data=falsepos_alleles[falsepos_alleles$nonref_allele_count_obs == 1,], 
       aes(x=pool, y=percent_falsepos, colour = variantcaller)) + 
  geom_line() + 
  geom_point() +
  labs(x="Number of samples in pool", y="False postives %", title = "Variant loci with one allele present in the pool")

ggplot(data=falsepos_alleles[falsepos_alleles$nonref_allele_count_obs == 1,], 
       aes(x=pool, y=n_falsepos, colour = variantcaller)) + 
  geom_line() + 
  geom_point() +
  labs(x="Number of samples in pool", y="False positives", title = "Variant loci with one allele present in the pool")


# HC_vars$variantcaller = "GATK Haplotype Caller"
# FBqual_vars$variantcaller = "FreeBayes QUAL>20"
# vars_recovered = rbind(HC_vars, FBqual_vars)



HC_falsepos %>% group_by(pool, nonref_allele_count_obs) %>% 
  summarise(n_false_positive = sum(false_positive), total = length(false_positive)) -> HC_falsepos_alleles
HC_falsepos_alleles$percent_false_positive = HC_falsepos_alleles$n_false_positive/HC_falsepos_alleles$total
HC_falsepos_alleles$pool_id = paste0('Pool_', HC_falsepos_alleles$pool)
HC_falsepos_alleles$variantcaller = "GATK Haplotype Caller"

FBqual_falsepos %>% group_by(pool, nonref_allele_count_obs) %>% 
  summarise(n_false_positive = sum(false_positive), total = length(false_positive)) -> FBqual_falsepos_alleles
FBqual_falsepos_alleles$percent_false_positive = FBqual_falsepos_alleles$n_false_positive/FBqual_falsepos_alleles$total
FBqual_falsepos_alleles$pool_id = paste0('Pool_', FBqual_falsepos_alleles$pool)
FBqual_falsepos_alleles$variantcaller = "FreeBayes QUAL>20"

falsepos_alleles = rbind(HC_falsepos_alleles, FBqual_falsepos_alleles)

ggplot(data=falsepos_alleles, aes(x=nonref_allele_count_obs, y=percent_false_positive, 
                              colour=pool_id)) + 
  geom_point() + geom_line() + facet_grid(~variantcaller, scales = 'free')
