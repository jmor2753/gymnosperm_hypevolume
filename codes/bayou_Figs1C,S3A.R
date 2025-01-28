### HPC Cluster ###

rm(list = ls())

set.seed(1)

setwd("~/gymnosperm_niches")

#install.packages("fitdistrplus", lib = "/uoa/home/c02fc23/R/x86_64-pc-linux-gnu-library/4.2")

library(phytools, lib.loc = "/uoa/home/c02fc23/R/x86_64-pc-linux-gnu-library/4.2")
library(dplyr, lib.loc = "/uoa/home/c02fc23/R/x86_64-pc-linux-gnu-library/4.2")
library(bayou, lib.loc = "/uoa/home/c02fc23/R/x86_64-pc-linux-gnu-library/4.2")

tr_gen <- read.tree("data/tr_gen.tre")
tr_fam <- read.tree("data/tr_fam.tre")

load(file = "data/hypervol_df_gen_50.RData")
load(file = "data/hypervol_df_fam_50.RData")
load(file = "data/hypervol_df_gen_out.RData")
load(file = "data/hypervol_df_fam_out.RData")

## summary dimensions
summary_persistence <- function(df) {
  output <- c()
  
  if ("family" %in% colnames(df)) {
    ## for loop for pairwise interactions
    for (i in 1:nrow(df)) {
      outdata <- df$persistenceHomology[[i]]$persHomol %>%
        mutate(surv = death - birth) %>%
        group_by(dimension) %>%
        summarise(avg_surv = mean(surv),
                  sd_surv = sd(surv),
                  max_surv = max(surv)) %>%
        mutate(family = df$family[i])
      output <- as.data.frame(bind_rows(output, outdata))
    }
    output <- output[c(5, 1, 2, 3, 4)]
  } else {
    for (i in 1:nrow(df)) {
      outdata <- df$persistenceHomology[[i]]$persHomol %>%
        mutate(surv = death - birth) %>%
        group_by(dimension) %>%
        summarise(avg_surv = mean(surv),
                  sd_surv = sd(surv),
                  max_surv = max(surv)) %>%
        mutate(genus = df$genus[i])
      output <- as.data.frame(bind_rows(output, outdata))
    }
    output <- output[c(5, 1, 2, 3, 4)]
  }
  return(as.data.frame(output))
}

## calculating average survival of connected components
persistence_gen_dt <- summary_persistence(rbind(hypervol_df_gen_50, hypervol_df_gen_out)) 
persistence_fam_dt <- summary_persistence(rbind(hypervol_df_fam_50, hypervol_df_fam_out)) 

volume_hyper_gen <-  rbind(hypervol_df_gen_50, hypervol_df_gen_out) %>%
  mutate(hv_vol = purrr::map(hypervol, function(.x) {
    .x@Volume
  })) %>%
  dplyr::select(genus, hv_vol) %>%
  tidyr::unnest(hv_vol) %>%
  right_join(., persistence_gen_dt) 

volume_hyper_fam <-  rbind(hypervol_df_fam_50, hypervol_df_fam_out) %>%
  mutate(hv_vol = purrr::map(hypervol, function(.x) {
    .x@Volume
  })) %>%
  dplyr::select(family, hv_vol) %>%
  tidyr::unnest(hv_vol) %>%
  right_join(., persistence_fam_dt) 

tr_gen <- reorder(tr_gen, "postorder")
tr_fam <- reorder(tr_fam, "postorder")

hyp_gen <- volume_hyper_gen$hv_vol
names(hyp_gen) <- volume_hyper_gen$genus
hyp_gen <- hyp_gen[!duplicated(names(hyp_gen))]
hyp_gen <- log(hyp_gen)

hyp_fam <- volume_hyper_fam$hv_vol
names(hyp_fam) <- volume_hyper_fam$family
hyp_fam <- hyp_fam[!duplicated(names(hyp_fam))]
hyp_fam <- log(hyp_fam)

priorOU_gen <- make.prior(tr_gen, 
  dists = list(dalpha = "dhalfcauchy", 
               dsig2 = "dhalfcauchy", 
               dk = "cdpois", dtheta = "dnorm"),
  param = list(dalpha = list(scale = 0.1), 
               dsig2 = list(scale = 0.1),
               dk = list(lambda = 10, kmax = 50), 
               dsb = list(bmax = 1, prob = 1), 
                dtheta = list(mean = mean(hyp_gen), sd = 1.5*sd(hyp_gen)))
  )

priorOU_fam <- make.prior(tr_fam, 
  dists = list(dalpha = "dhalfcauchy", 
               dsig2 = "dhalfcauchy", 
               dk = "cdpois", dtheta = "dnorm"),
  param = list(dalpha = list(scale = 0.1), 
               dsig2 = list(scale = 0.1),
               dk = list(lambda = 10, kmax = 50), 
               dsb = list(bmax = 1, prob = 1), 
                dtheta = list(mean = mean(hyp_fam), sd = 1.5*sd(hyp_fam)))
  )

mcmcOU_gen <- bayou.makeMCMC(tr_gen, hyp_gen, prior = priorOU_gen, 
                             outname = "modelOU_gen001", plot.freq = NULL) 
mcmcOU_fam <- bayou.makeMCMC(tr_fam, hyp_fam, prior = priorOU_fam, 
                             outname = "modelOU_fam001", plot.freq = NULL) 

mcmcOU_gen$run(1000000) # Run the MCMC
mcmcOU_fam$run(1000000) # Run the MCMC

save(mcmcOU_gen, file = "data/mcmcOU_gen.RData")
save(mcmcOU_fam, file = "data/mcmcOU_fam.RData")

chainOU_gen <- mcmcOU_gen$load()
chainOU_fam <- mcmcOU_fam$load()

save(chainOU_gen, file = "data/chainOU_gen.RData")
save(chainOU_fam, file = "data/chainOU_fam.RData")

pdf("bayou_conv_gen.pdf", height = 15)
plot(chainOU_gen, auto.layout = T)
dev.off()

pdf("bayou_conv_fam.pdf", height = 15)
plot(chainOU_fam, auto.layout = T)
dev.off()

shiftsums_gen <- shiftSummaries(chainOU_gen, mcmcOU_gen)
shiftsums_fam <- shiftSummaries(chainOU_fam, mcmcOU_fam)

pdf("Fig1C4.pdf")
par(mar = c(0, 0, 0, 0))
plotSimmap.mcmc(chainOU_gen, burnin = 0.3, pp.cutoff = 0.3, cex = 0.7)
par(mar = c(0, 3, 0, 0))
plotBranchHeatMap(tr_gen, chainOU_gen, "theta", burnin = 0.3, pal = cm.colors, cex = 0.7)
par(mar = c(4, 4, 4, 4))
phenogram.density(tr_gen, hyp_gen, burnin = 0.3, chainOU_gen, pp.cutoff = 0.3)
dev.off()

pdf("FigS3A2.pdf")
par(mar = c(0, 0, 0, 0))
plotSimmap.mcmc(chainOU_fam, burnin = 0.3, pp.cutoff = 0.3, cex = 0.7)
par(mar = c(0, 3, 0, 0))
plotBranchHeatMap(tr_fam, chainOU_fam, "theta", burnin = 0.3, pal = cm.colors, cex = 0.7)
par(mar = c(4, 4, 4, 4))
phenogram.density(tr_fam, hyp_fam, burnin = 0.3, chainOU_fam, pp.cutoff = 0.3)
dev.off()

### local computer ###

rm(list = ls())

set.seed(1)

setwd("~/Documents/lab/gymnosperm_niches")

library(phytools)
library(dplyr)
library(geiger)

source("codes/functions.R")

load("data/chainOU1_gen.RData")
load("data/chainOU2_gen.RData")
load("data/chainOU3_gen.RData")
load("data/chainOU1_fam.RData")
load("data/chainOU2_fam.RData")
load("data/chainOU3_fam.RData")

load("data/mcmcOU1_gen.RData")
load("data/mcmcOU2_gen.RData")
load("data/mcmcOU3_gen.RData")
load("data/mcmcOU1_fam.RData")
load("data/mcmcOU2_fam.RData")
load("data/mcmcOU3_fam.RData")

shiftsums1_gen <- shiftSummaries(chainOU1_gen, mcmcOU1_gen)
shiftsums2_gen <- shiftSummaries(chainOU2_gen, mcmcOU2_gen)
shiftsums3_gen <- shiftSummaries(chainOU3_gen, mcmcOU3_gen)
shiftsums1_fam <- shiftSummaries(chainOU1_fam, mcmcOU1_fam)
shiftsums2_fam <- shiftSummaries(chainOU2_fam, mcmcOU2_fam)
shiftsums3_fam <- shiftSummaries(chainOU3_fam, mcmcOU3_fam)

plotSimmap.mcmc(chainOU_gen, burnin = 0.3, pp.cutoff = 0.3)


pdf("run1.pdf")
par(mfrow = c(2,2))
plotSimmap.mcmc(chainOU, burnin = 0.3, pp.cutoff = 0.3)

plotBranchHeatMap(tr_gen, chainOU, "theta", burnin = 0.3, pal = cm.colors)
phenogram.density(tr_gen, hyp_gen, burnin = 0.3, chainOU, pp.cutoff = 0.3)
dev.off()











## phylogenetic signal centroid
load(file = "data/env_df_gen_dist.RData")
load(file = "data/env_df_fam_dist.RData")

dist_cent_gen <- do.call(c, env_df_gen_dist$dist_cent)
names(dist_cent_gen) <- env_df_gen_dist$genus

dist_cent_fam <- do.call(c, env_df_fam_dist$dist_cent)
names(dist_cent_fam) <- env_df_fam_dist$family

phylosig(keep.tip(tr_gen, names(dist_cent_gen)),
         log(dist_cent_gen), method = "K", test = T, nsim = 1000) 
phylosig(keep.tip(tr_gen, names(dist_cent_gen)),
         log(dist_cent_gen), method = "lambda", test = T, nsim = 1000) 

phylosig(keep.tip(tr_fam, names(dist_cent_fam)),
         log(dist_cent_fam), method = "K", test = T, nsim = 1000) 
phylosig(keep.tip(tr_fam, names(dist_cent_fam)),
         log(dist_cent_fam), method = "lambda", test = T, nsim = 1000) 

plot_centroid_gen <- do.call(rbind, env_df_gen_dist$centroid)
plot_centroid_fam <- do.call(rbind, env_df_fam_dist$centroid)

mean_all <- c(1.150726e-13, -5.372683e-15, -3.000657e-15)

pdf("figures/FigS1D_gen.pdf", height = 5, width = 7)

ggplot(plot_centroid_gen, aes(PC1, PC2, col = PC3)) +
  geom_point(size = 3) +
  geom_point(aes(x = mean_all[1], y = mean_all[2], col = mean_all[3]), shape = 8, size = 5) + 
  scale_color_viridis_c() +
  theme_bw(base_size = 19) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 9))

dev.off()

pdf("figures/FigS1E_fam.pdf", height = 5, width = 7)

ggplot(plot_centroid_fam, aes(PC1, PC2, col = PC3)) +
  geom_point(size = 3) +
  geom_point(aes(x = mean_all[1], y = mean_all[2], col = mean_all[3]), shape = 8, size = 5) + 
  scale_color_viridis_c() +
  theme_bw(base_size = 19) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 9))

dev.off()
