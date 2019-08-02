## simurg: Paper figures ##

# Load simurg
library(simurg)

# Set test env
wd <- tempdir()
tgz <- system.file('extdata', 'ref_tutorial.tar.gz', package = 'simurg')
untar(tarfile = tgz, exdir = wd)
ref <- list.files(path = wd, pattern = 'ref_tutorial[.]fasta$', full.names = TRUE)
set.seed(555)

# Set parameters

norg <- 20
ne <- 1e11
C <- 100
u <-  1e-8; theta <- 2 * ne * u
v <-  1e-11; rho <- 2 * ne * v
mu <- 5e-12




#############
# FIGURE 1A #
#############
# Mean genetic distances vs. Normalized cophenetic distances

pg <- simpg(ref = ref,
            norg = norg,
            ne = ne,
            C = C,
            u = u,
            v = v,
            mu = mu,
            dir_out = paste(wd, 'dir_out', sep = '/'),
            replace = TRUE,
            force = TRUE,
            verbose = FALSE)

spg <- summary(pg)
dis <- spg$Evo_dist

pdf('Figure1A.pdf')

par(mar = c(5,6,1,1))

plot(x = dis$norm_cophenetic,
     y = dis$mean_gene_dist,
     xlim = c(0, 1), xlab = 'Normalized cophenetic distances',
     ylim = c(0, 1), ylab = 'Mean genetic distances',
     cex.lab = 1.3, pch = 19, las = 1,cex = 1.7, col =' grey',
     xaxt='n',yaxt='n')

axis(1, at = seq(0 ,1 , 0.2), las = 1, lwd.ticks = 2, cex.axis = 1.4)
axis(2, at = seq(0 ,1 , 0.2), las = 1, lwd.ticks = 2, cex.axis = 1.4)

abline(0, 1, lwd = 2,lty = 2)

box(lwd = 2)

dev.off()




#############
# FIGURE 1B #
#############
# Gene family frequency spectrum

# Replicates
reps <- 1:10

# Matrix to store results
resu <- matrix(0L, nrow = length(reps), ncol = norg,
               dimnames = list(paste0('rep', reps), 1:norg))

for (i in reps){

  pg <- simpg(ref = ref,
              norg = norg,
              ne = ne,
              C = C,
              u = u,
              v = v,
              mu = 5e-12,
              dir_out = paste(wd, 'dir_out', sep = '/'),
              replace = TRUE,
              force = TRUE,
              verbose = FALSE)
  spg <- summary(pg)
  Ofreq <- spg$Gene_family_frequency
  resu[i, names(Ofreq)] <- Ofreq

}

# Load function to calculate gene family frequency spectrum:
source('https://raw.githubusercontent.com/rec3141/pangenome/master/f-pangenome.R')
# Testing the Infinitely Many Genes Model of the Bacterial Pangenome. R. Eric
# Collins and Paul G. Higgs. Molecular Biology and Evolution, Volume 29, Issue
# 11, November 2012, Pages 3413–3425, https://doi.org/10.1093/molbev/mss163
Efreq <- f.coalescent.spec(x = c(rho, theta, C), ng = norg)
mx <- max(c(resu, Efreq))

pdf('Figure1B.pdf')

par(mar = c(5,7,1,1))

boxplot(resu, xlab = 'Number of genomes', ylab= 'Gene family frequency',
        ylim = c(0L, mx*1.1),
        boxcol = 'grey60', medcol  = 'grey60',
        whiskcol = 'grey60', staplecol='grey60', lwd = 3,
        outcol = 'grey60', cex.lab = 1.3, xaxt = 'n', yaxt = 'n')

axis(1, at = seq(0 ,20, 5), las = 1, lwd.ticks = 2, cex.axis = 1.4 )
axis(2, at = seq(0 ,3500 , 500), las = 1, lwd.ticks = 2, cex.axis = 1.4)

lines(Efreq, lwd = 2, lty = 2)

box(lwd = 2)

dev.off()
