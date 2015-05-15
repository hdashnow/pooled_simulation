# Cost estimates
# Assumes minimum submission of 12 samples and ~150X

# Cost of the library prep
lib.prep.cost = 450 
# $350 per 4 Gb of sequencing. 8 Gb standard for 150X
seq.Gb.cost = 350/4 # cost per Gb of sequencing
# cost per average 1X of exome coverage (assuming standard exome is 8Gb and 150X - it's actually 153X)
seq.cov.cost = seq.Gb.cost*8/ 150


pool_recovery <- data.frame(pool.size = c(2,4,6,8,10), 
                            prop.recovered = c(0.947957966, 0.898628636, 0.568131842, 
                                               0.307364446, 0.171492227) )
par(mar=c(5,7,4,2))

pdf("plots/pooled_genotyping_simulation_allCHR.pdf")
plot(pool_recovery, type = 'b', xlab="Number of individuals in pool", 
     ylab="Proportion of loci called \n (in pool compared to same individuals sequenced separately)", 
     main="Pooled genotyping simulation")
dev.off()


#--------------------------------------------------------

# Now downsample the reads, but maintain the same ploidy (6 samples in this case)

pool_recovery <- data.frame(downsampling = c(1,2,3,4,5,6), 
                            percent.recovered = c( 92, 87.8, 83.0, 77.9, 75.6, 68.8) )

pool_recovery$prop.recovered = pool_recovery$percent.recovered/100
pool_recovery$downsampling = 1/pool_recovery$downsampling
pool_recovery <- subset(pool_recovery, select = -c(percent.recovered) )

pdf("plots/pooled_genotyping_simulation_allCHR_varyDepth.pdf")
par(mar=c(5,7,4,2))
plot(pool_recovery, type = 'b', 
     xlab="Downsampling: proportion of reads relative to orignal sequencing depth", 
     ylab="Proportion of loci called \n (in pool compared to same individuals sequenced separately)", 
     main="Pooled genotyping simulation\nAll pools have 6 samples but varying depth")
dev.off()

#--------------------------------------------------------

#Replotted with sequencing depth:
pool_recovery <- data.frame(downsampling = c(1,2,3,4,5,6), 
                            percent.recovered = c( 92, 87.8, 83.0, 77.9, 75.6, 68.8) )

pool_recovery$prop.recovered = pool_recovery$percent.recovered/100
pool_recovery$coverage = 1/pool_recovery$downsampling * 153*6 # Average mean coverage for these samples is 153
pool_recovery <- subset(pool_recovery, select = -c(percent.recovered,downsampling) )

pdf("plots/pooled_genotyping_simulation_allCHR_varyDepth2.pdf")
par(mar=c(5,7,4,2))
plot(x=pool_recovery$coverage, y=pool_recovery$prop.recovered, type = 'b',
     xlab="Mean sequencing depth",
     ylab="Proportion of loci called \n (in pool compared to same individuals sequenced separately)",
     main="Pooled genotyping simulation\nAll pools have 6 samples but varying depth")
dev.off()

#--------------------------------------------------------

#Fit a log curve:

xdata = pool_recovery$coverage
ydata = pool_recovery$prop.recovered

# some starting values
p1 = 0.1
p2 = 0.2

# do the fit
fit = nls(ydata ~ p1 + p2*log(xdata), start=list(p1=p1,p2=p2))

new = data.frame(xdata = seq(min(xdata),max(xdata),len=200))
lines(new$xdata,predict(fit,newdata=new))


#--------------------------------------------------------

pool_recovery <- data.frame(downsampling = c(1,2,3,4,5,6), percent.recovered = c( 92, 87.8, 83.0, 77.9, 75.6, 68.8) )

pool_recovery$prop.recovered = pool_recovery$percent.recovered/100
pool_recovery$coverage = 1/pool_recovery$downsampling * 153*6 # Average mean coverage for these samples is 153
pool_recovery <- subset(pool_recovery, select = -c(percent.recovered,downsampling) )

par(mar=c(5,7,4,2))
plot(x=pool_recovery$coverage, y=pool_recovery$prop.recovered, type = 'p',
     xlab="Mean sequencing depth",
     ylab="Proportion of loci called \n (in pool compared to same individuals sequenced separately)",
     main="Pooled genotyping simulation\nAll pools have 6 samples but varying depth")

curve(2*log(x)+2)

xdata = pool_recovery$coverage
ydata = pool_recovery$prop.recovered

# some starting values
p1 = 0.1
p2 = 0.2

# do the fit
fit = nls(ydata ~ p1 + p2*log(xdata), start=list(p1=p1,p2=p2))

new = data.frame(xdata = seq(min(xdata),max(xdata),len=200))
lines(new$xdata,predict(fit,newdata=new))

pool_recovery <- data.frame(pool.size = c(2,4,6,8,10),
                            prop.recovered = c(0.947957966, 0.898628636, 0.568131842, 0.307364446, 0.171492227) )
par(mar=c(5,7,4,2))
plot(pool_recovery, type = 'p', xlab="Number of individuals in pool",
     ylab="Proportion of loci called \n (in pool compared to same individuals sequenced separately)",
     main="Pooled genotyping simulation")

xdata = pool_recovery$pool.size
ydata = pool_recovery$prop.recovered

# some starting values
p1 = 0.1
p2 = 0.2

# do the fit
fit = nls(ydata ~ p1 + p2*log(xdata), start=list(p1=p1,p2=p2))
fit = nls(ydata ~ p1 + p2*xdata, start=list(p1=p1,p2=p2))

new = data.frame(xdata = seq(min(xdata),max(xdata),len=200))
lines(new$xdata,predict(fit,newdata=new))