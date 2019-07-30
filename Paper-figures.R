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

tiff('Figure1.tiff')
plot(x = dis$norm_cophenetic, 
     y = dis$mean_gene_dist, 
     xlim = c(0, 1), xlab = 'Normalized cophenetic distances',
     ylim = c(0, 1), ylab = 'Mean genetic distances', 
     cex.lab=1.3)
abline(0, 1)
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

tiff('Figure2.tiff')
boxplot(resu, xlab = 'Genomes', ylab= 'Gene family frequency',
        ylim = c(0L, mx*1.1),
        boxcol = 'grey50', medcol  = 'grey50', 
        whiskcol = 'grey50', staplecol='grey50', 
        outcol='grey50', cex.lab=1.3, xaxt = 'n')
axis(1, at = c(1, 5, 10, 15, 20))
lines(Efreq)
dev.off()
