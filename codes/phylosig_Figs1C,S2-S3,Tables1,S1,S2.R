rm(list = ls())

setwd("~/Documents/lab/gymnosperm_niches")

library(dplyr)
library(purrr)
library(TDAstats)
library(ggplot2)
library(patchwork)
library(phytools)
library(ggtree)
library(deeptime)
library(TDA)
library(plotrix)
library(viridis)

tr_gen <- read.tree("data/tr_gen.tre")
tr_fam <- read.tree("data/tr_fam.tre")

load(file = "data/cluster/gen/env_df_gen_hypervolume_50.RData")
load(file = "data/cluster/fam/env_df_fam_hypervolume_50.RData")

load(file = "data/cluster/gen/env_df_gen_hypervolume_out.RData")
load(file = "data/cluster/fam/env_df_fam_hypervolume_out.RData")

# persistence homology calculations
subsample_dist <- function(df, size, d, ...) {
  
  subsample_dist_demo <- function(x, size, d, d.max = NULL, replacement = FALSE, 
                                  latlong = FALSE, echo = FALSE) {
    
    if (missing(x))
      stop("Must define a spatial object")
    
    if (missing(d))
      stop("Must define minimum separation distance")
    
    if (!is.null(d.max)) {
      if (d.max <= d) stop("Maximum distance must be larger than minimum")
    }
    
    if (!any(class(x)[1] == c("SpatialPointsDataFrame", "SpatialPolygonsDataFrame")))
      stop("x must be sp class polygons or points")
    
    if (latlong == TRUE) {
      message("geographic projection distances must be in kilometers")
    }
    
    if (size >= nrow(x))
      stop("subsample size must be smaller than population")
    
    rs <- sample(1:nrow(x), 1)
    
    s <- x[rs, ]
    
    if (replacement == FALSE) {
      x <- x[-rs, ]
    }
    
    deval = TRUE
    
    for (i in 2:size) {
      nsamp = 0
      
      while (deval == TRUE) {
        rs <- sample(1:nrow(x), 1)
        pts.dist = sp::spDists(s, x[rs, ], longlat = latlong)
        
        if (is.null(d.max)) {
          deval <- any(pts.dist < d, na.rm = TRUE)
        } else {
          deval <- any(pts.dist < d, na.rm = TRUE) | any(pts.dist > d.max, 
                                                         na.rm = TRUE)
        }
        
        nsamp = nsamp + 1
        
        if (echo)
          cat("Sample iteration=", nsamp, "\n")
        
        if (nsamp == nrow(x))
          break
      }
      
      if (echo) {
        cat("\n", "Min distance for", i, "=", min(pts.dist, na.rm = TRUE), "\n")
        cat(" Max distance for", i, "=", max(pts.dist, na.rm = TRUE), "\n")
      }
      
      if (nsamp == nrow(x)) {
        message(paste0("Warning: sampling cannot converge at n=", size, 
                       " returning n=", nrow(s)))
        return(s)
      }
      
      deval = TRUE
      
      s <- rbind(s, x[rs, ])
      
      if (replacement == FALSE) {
        x <- x[-rs, ]
      }
      
    }
    
    return(s)
    
  }
  
  sub.meuse <- c()
  
  for(i in 1:ncol(df)){
    
    for(j in 1:ncol(df)) {
      
      if(i == j | j > i){
        next
      } else {
        testsp <- df
        
        sp::coordinates(testsp) <- c(i, j)
        
        sub <- subsample_dist_demo(testsp, size = size/ncol(df), d = d) ##d is for the minimum distance between the points sampled
        
        sub.meuse <- bind_rows(sub.meuse, as.data.frame(sub))
      }
      
    }
    
  }
  
  return(sub.meuse)
  
}

calc_ph <- function(data, tree) {

  ph <- data %>%
  mutate(randomHyper = purrr::map(hypervol, function(.x) {
    out <- apply(.x@RandomPoints, 2, base::scale)
  })) %>%
  mutate(persistenceHomology = purrr::map(randomHyper, function(.x) {
    set.seed(145) 
    subset_hyper <- as.data.frame(.x) 
    subset_hyper <- subsample_dist(subset_hyper, size = 300, d = 0.2) 
    PH_calc <- calculate_homology(subset_hyper, dim = 2)
    return(list(hyper_subset = subset_hyper,
                persHomol = as.data.frame(PH_calc)))
  }))

  ph_df <- ph %>%
    mutate(hull = purrr::map(randomHyper, function(.x) {
      hullout <- with(as.data.frame(.x), chull(PC1, PC2))
      hullout
    }), out = map2(randomHyper, hull, ~ .x[.y,, drop = FALSE]))

  if ("family" %in% colnames(ph_df)) {
    index <- which(ph_df$family %in% tree$tip.label)
    ph_df <- ph_df[index,]
  } else {
    index <- which(ph_df$genus %in% tree$tip.label)
    ph_df <- ph_df[index,]
  }

  return(ph_df)

}

hypervol_df_gen_50 <- calc_ph(env_df_gen_hypervolume_50, tr_gen)
#save(hypervol_df_gen_50, file = "data/hypervol_df_gen_50.RData")

hypervol_df_fam_50 <- calc_ph(env_df_fam_hypervolume_50, tr_fam)
#save(hypervol_df_fam_50, file = "data/hypervol_df_fam_50.RData")

hypervol_df_gen_out <- calc_ph(env_df_gen_hypervolume_out, tr_gen)
#save(hypervol_df_gen_out, file = "data/hypervol_df_gen_out.RData")

hypervol_df_fam_out <- calc_ph(env_df_fam_hypervolume_out, tr_fam)
#save(hypervol_df_fam_out, file = "data/hypervol_df_fam_out.RData")

#load(file = "data/hypervol_df_gen_50.RData")
#load(file = "data/hypervol_df_fam_50.RData")

#load(file = "data/hypervol_df_gen_out.RData")
#load(file = "data/hypervol_df_fam_out.RData")

## figure 1c

hypervol_df_gen_50_2 <- hypervol_df_gen_50 %>% 
  filter(genus %in% c("Ginkgo", "Juniperus", "Pinus", "Taxus")) %>%
  arrange(genus)

plot_points_gen <- list()
for (i in 1:nrow(hypervol_df_gen_50_2)) {

  plot_points_gen[[i]] <- ggplot(as.data.frame(hypervol_df_gen_50_2$out[[i]]), 
                             aes(PC1, PC2, col = PC3)) +
    geom_polygon(alpha = 0.05, fill = "royalblue2", col = NA) +
    geom_point(data = as.data.frame(hypervol_df_gen_50_2$randomHyper[[i]]),
               aes(PC1, PC2, col = PC3), inherit.aes = FALSE, alpha = 0.07, 
               size = 0.2) +
    geom_point(data = hypervol_df_gen_50_2$persistenceHomology[[i]]$hyper_subset,
               aes(PC1, PC2, col = PC3), inherit.aes = FALSE, alpha = 1, 
               size = 1.5) +
    ggtitle(paste0(hypervol_df_gen_50_2$genus[i])) +
    scale_color_viridis_c() +
    theme_bw(base_size = 16) +
    theme(panel.grid = element_blank(),
          #axis.text = element_text(size = 8),
          #axis.title = element_text(size = 9),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 9))
}

# plots of the persistence diagrams
plot_persistence_gen <- list()
for (i in 1:nrow(hypervol_df_gen_50_2)) {
  plot_persistence_gen[[i]] <- ggplot(data = hypervol_df_gen_50_2$persistenceHomology[[i]]$persHomol,
                                  aes(x = birth, 
                                      y = death, 
                                      shape = as.factor(dimension), 
                                      color = as.factor(dimension))) +
    geom_point(alpha = 0.9) +
    geom_abline(slope = 1, size = 0.2) +
    xlab("Birth") +
    ylab("Death") +
    theme_linedraw(base_size = 16) +
    theme(panel.grid = element_blank()) +
    scale_color_manual("Dimension", 
                       values = c("darkseagreen3", "royalblue2", "hotpink")) +
    scale_shape_manual("Dimension", values = c(15, 17, 19))
}

final_plot_gen <- c(rbind(plot_points_gen, plot_persistence_gen))

pdf("figures/Fig1C3.pdf", height = 7, width = 19)
wrap_plots(final_plot_gen) + plot_layout(ncol = 4)
dev.off()

## figure 1c - phylogeny

tr_gen_ladd <- ladderize(tr_gen, right = FALSE)
tr_gen_ladd <- ips::fixNodes(tr_gen)

## plotting the phylogenies
pdf("figures/Fig1C1.pdf", width = 8, height = 10)
revts(ggtree(tr_gen_ladd) + geom_tiplab(size = 3)) +
  coord_geo(xlim = c(-450, 100), ylim = c(0, 75), neg = TRUE, abbrv = F, 
            dat = "period") +
  scale_x_continuous(breaks = seq(-450, 0, 25), labels = abs(seq(-450, 0, 25))) +
  theme_tree2()
dev.off()

## function to estimate the distance between two persistence diagrams
homology_distance <- function(x, y, dim = 2, ...) {
  output <- c()
  for (i in 0:dim) {
    dist <- TDA::bottleneck(as.matrix(x), as.matrix(y), i)
    output <- cbind(output, dist)
    colnames(output)[i + 1] <- paste(i)
  }
  return(as.data.frame(output))
}

## pairwise distance of the persistence (bottleneck distances)
pairwise_persistence <- function(df) {
  output <- c()
  
  if ("family" %in% colnames(df)) {
    # for loop for pairwise interactions
    for (i in 1:nrow(df)) {
      for (j in 1:nrow(df)) {
        if (i <= j) { # excluding self-comparison and double-comparison
          next
        } else {
          pairwise_dist_out <- homology_distance(df$persistenceHomology[[i]]$persHomol,
                                                 df$persistenceHomology[[j]]$persHomol)
          fam1 <- df$family[j]
          fam2 <- df$family[i]
            
          outdata <- data.frame(family1 = fam1, family2 = fam2, dist = pairwise_dist_out)
          output <- as.data.frame(bind_rows(output, outdata))
        }
      }
    }

    output_final <- data.frame(family1 = output$family2, 
                               family2 = output$family1) %>%
      bind_cols(., output[3:5]) %>%
      bind_rows(., output)
  } else {
    # for loop for pairwise interactions
    for (i in 1:nrow(df)) {
      for (j in 1:nrow(df)) {
        if (i <= j) { # excluding self-comparison and double-comparison
          next
        } else {
          pairwise_dist_out <- homology_distance(df$persistenceHomology[[i]]$persHomol,
                                                 df$persistenceHomology[[j]]$persHomol)
          gen1 <- df$genus[j]
          gen2 <- df$genus[i]
            
          outdata <- data.frame(genus1 = gen1, genus2 = gen2, dist = pairwise_dist_out)
          output <- as.data.frame(bind_rows(output, outdata))
        }
      }
    }

    output_final <- data.frame(genus1 = output$genus2, 
                               genus2 = output$genus1) %>%
      bind_cols(., output[3:5]) %>%
      bind_rows(., output)
  }
  
  return(as.data.frame(output_final))
}

pairw_persist_gen <- pairwise_persistence(rbind(hypervol_df_gen_50, hypervol_df_gen_out))
colnames(pairw_persist_gen) <- c("genus1", "genus2", "dist0", "dist1", "dist2")
#save(pairw_persist_gen, file = "data/pairw_persist_gen.RData")

pairw_persist_gen$genus1 <- factor(pairw_persist_gen$genus1,
  levels = c((hypervol_df_gen_50 %>% arrange(genus))$genus,
             (hypervol_df_gen_out %>% arrange(genus))$genus))

pairw_persist_gen$genus2 <- factor(pairw_persist_gen$genus2,
 levels = c((hypervol_df_gen_50 %>% arrange(genus))$genus,
            (hypervol_df_gen_out %>% arrange(genus))$genus))

plot_gen_dist2 <- ggplot(pairw_persist_gen, aes(genus1, genus2, fill = dist2)) +
  geom_tile() +
  scale_fill_viridis_c("Distance\n(Dim 2)", option = "mako") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.0)) +
  xlab("Genus 1") + 
  ylab("Genus 2")

pdf("figures/FigS2E.pdf", height = 10, width = 10)
plot_gen_dist2
dev.off()

pairw_persist_fam <- pairwise_persistence(rbind(hypervol_df_fam_50, hypervol_df_fam_out))
colnames(pairw_persist_fam) <- c("famus1", "famus2", "dist0", "dist1", "dist2")
#save(pairw_persist_fam, file = "data/pairw_persist_fam.RData")

pairw_persist_fam$family1 <- factor(pairw_persist_fam$famus1,
                                   levels = c((hypervol_df_fam_50 %>% arrange(family))$family,
                                              (hypervol_df_fam_out %>% arrange(family))$family))

pairw_persist_fam$family2 <- factor(pairw_persist_fam$famus2,
                                   levels = c((hypervol_df_fam_50 %>% arrange(family))$family,
                                              (hypervol_df_fam_out %>% arrange(family))$family))

plot_fam_dist2 <- ggplot(pairw_persist_fam, aes(family1, family2, fill = dist2)) +
  geom_tile() +
  scale_fill_viridis_c("Distance\n(Dim 2)", option = "mako") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.0)) +
  xlab("famus 1") + 
  ylab("famus 2")

pdf("figures/FigS3F.pdf", height = 10, width = 10)
plot_fam_dist2
dev.off()

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

## calculate phylogenetic signal
calculate_physig <- function(tree, data, level, variable, method) {
  
  output <- c()
  inner <- c()

  sub_dat <- as.data.frame(subset(data, dimension == 2, select = variable))[, 1]
  names(sub_dat) <- as.data.frame(subset(data, dimension == 2, select = level))[, 1]
  
  if (method == "K") {
    signalphy <- phylosig(tree, log(sub_dat), method = "K", test = T, 
                          nsim = 1000) 
    inner <- bind_rows(inner, data.frame(dim = 2,
                                         estimate = signalphy$K,
                                         pval = signalphy$P))
  } else {
    signalphy <- phylosig(tree, log(sub_dat), method = "lambda", test = T, 
                          niter = 1000) 
    inner <- bind_rows(inner, data.frame(dim = 2,
                                         estimate = signalphy$lambda,
                                         logL0 = signalphy$logL0,
                                         logL = signalphy$logL,
                                         pval = signalphy$P))
  }
  
  output <- bind_rows(output, inner)
  return(output)
}

hypervolume_sig_K_gen <- calculate_physig(tr_gen, volume_hyper_gen, 
                                          "genus", "hv_vol", "K") %>%
  mutate(type = "Hypervolume")
hypervolume_sig_L_gen <- calculate_physig(tr_gen, volume_hyper_gen,
                                          "genus", "hv_vol", "lambda") %>%
  mutate(type = "Hypervolume")

avgdist_sig_K_gen <- calculate_physig(tr_gen, volume_hyper_gen, 
                                      "genus", "avg_surv", "K") %>%
  mutate(type = "Mean dist")
avgdist_sig_L_gen <- calculate_physig(tr_gen, volume_hyper_gen,
                                      "genus", "avg_surv", "lambda") %>%
  mutate(type = "Mean dist")

sddist_sig_K_gen <- calculate_physig(tr_gen, volume_hyper_gen, 
                                     "genus", "sd_surv", "K") %>%
  mutate(type = "SD dist")
sddist_sig_L_gen <- calculate_physig(tr_gen, volume_hyper_gen,
                                     "genus", "sd_surv", "lambda") %>%
  mutate(type = "SD dist")

maxdist_sig_K_gen <- calculate_physig(tr_gen, volume_hyper_gen, 
                                      "genus", "max_surv", "K") %>%
  mutate(type = "Max dist")
maxdist_sig_L_gen <- calculate_physig(tr_gen, volume_hyper_gen,
                                      "genus", "max_surv", "lambda") %>%
  mutate(type = "Max dist")

final_signal_gen <- bind_rows(cbind(method = "K", hypervolume_sig_K_gen),
                              cbind(method = "L", hypervolume_sig_L_gen),
                              cbind(method = "K", avgdist_sig_K_gen),
                              cbind(method = "L", avgdist_sig_L_gen),
                              cbind(method = "K", sddist_sig_K_gen),
                              cbind(method = "L", sddist_sig_L_gen),
                              cbind(method = "K", maxdist_sig_K_gen),
                              cbind(method = "L", maxdist_sig_L_gen))

table_physig_gen <- final_signal_gen %>%
  group_by(method, type) %>%
  dplyr::select(-dim) 

write.csv(table_physig_gen, "tables/Table1_unformatted.csv")

hypervolume_sig_K_fam <- calculate_physig(tr_fam, volume_hyper_fam,
                                          "family", "hv_vol", "K") %>%
  mutate(type = "Hypervolume")
hypervolume_sig_L_fam <- calculate_physig(tr_fam, volume_hyper_fam, 
                                          "family", "hv_vol", "lambda") %>%
  mutate(type = "Hypervolume")

avgdist_sig_K_fam <- calculate_physig(tr_fam, volume_hyper_fam,
                                      "family", "avg_surv", "K") %>%
  mutate(type = "Mean dist")
avgdist_sig_L_fam <- calculate_physig(tr_fam, volume_hyper_fam, 
                                      "family", "avg_surv", "lambda") %>%
  mutate(type = "Mean dist")

sddist_sig_K_fam <- calculate_physig(tr_fam, volume_hyper_fam,
                                     "family", "sd_surv", "K") %>%
  mutate(type = "SD dist")
sddist_sig_L_fam <- calculate_physig(tr_fam, volume_hyper_fam, 
                                     "family", "sd_surv", "lambda") %>%
  mutate(type = "SD dist")

maxdist_sig_K_fam <- calculate_physig(tr_fam, volume_hyper_fam,
                                      "family", "max_surv", "K") %>%
  mutate(type = "Max dist")
maxdist_sig_L_fam <- calculate_physig(tr_fam, volume_hyper_fam, 
                                      "family", "max_surv", "lambda") %>%
  mutate(type = "Max dist")

final_signal_fam <- bind_rows(cbind(method = "K", hypervolume_sig_K_fam),
                              cbind(method = "L", hypervolume_sig_L_fam),
                              cbind(method = "K", avgdist_sig_K_fam),
                              cbind(method = "L", avgdist_sig_L_fam),
                              cbind(method = "K", sddist_sig_K_fam),
                              cbind(method = "L", sddist_sig_L_fam),
                              cbind(method = "K", maxdist_sig_K_fam),
                              cbind(method = "L", maxdist_sig_L_fam))

table_physig_fam <- final_signal_fam %>%
  group_by(method, type) %>%
  dplyr::select(-dim) 

write.csv(table_physig_fam, "tables/TableS1_unformatted.csv")

## figure 1C, S2A-C, S3A-D

hyp_gen <- log(volume_hyper_gen$hv_vol[volume_hyper_gen$dimension == 2])
names(hyp_gen) <- volume_hyper_gen$genus[volume_hyper_gen$dimension == 2]
hyp_gen <- hyp_gen[tr_gen_ladd$tip.label]

avg_gen <- log(volume_hyper_gen$avg_surv[volume_hyper_gen$dimension == 2])
names(avg_gen) <- volume_hyper_gen$genus[volume_hyper_gen$dimension == 2]
avg_gen <- avg_gen[tr_gen_ladd$tip.label]

max_gen <- log(volume_hyper_gen$max_surv[volume_hyper_gen$dimension == 2])
names(max_gen) <- volume_hyper_gen$genus[volume_hyper_gen$dimension == 2]
max_gen <- max_gen[tr_gen_ladd$tip.label]

sd_gen <- log(volume_hyper_gen$sd_surv[volume_hyper_gen$dimension == 2])
names(sd_gen) <- volume_hyper_gen$genus[volume_hyper_gen$dimension == 2]
sd_gen <- sd_gen[tr_gen_ladd$tip.label]

hyp_gen_contMap <- contMap(tr_gen_ladd, hyp_gen, plot = F, res = 200)
hyp_gen_contMap <- setMap(hyp_gen_contMap, turbo(100))
hyp_gen_contMap$tree <- ladderize.simmap(hyp_gen_contMap$tree, right = F)

avg_gen_contMap <- contMap(tr_gen_ladd, avg_gen, plot = F, res = 200)
avg_gen_contMap <- setMap(avg_gen_contMap, turbo(100))
avg_gen_contMap$tree <- ladderize.simmap(avg_gen_contMap$tree, right = F)

max_gen_contMap <- contMap(tr_gen_ladd, max_gen, plot = F, res = 200)
max_gen_contMap <- setMap(max_gen_contMap, turbo(100))
max_gen_contMap$tree <- ladderize.simmap(max_gen_contMap$tree, right = F)

sd_gen_contMap <- contMap(tr_gen_ladd, sd_gen, plot = F, res = 200)
sd_gen_contMap <- setMap(sd_gen_contMap, turbo(100))
sd_gen_contMap$tree <- ladderize.simmap(sd_gen_contMap$tree, right = F)

pdf("figures/Fig1C2.pdf", height = 15, width = 12)

plot(hyp_gen_contMap, fsize = c(0.7, 0.8), leg.txt = "Hypervolume",
     outline = F, legend = 300)

dev.off()

pdf("figures/FigS2A-C.pdf", height = 15, width = 12)

layout(matrix(1:4, ncol = 2))

plot(avg_gen_contMap, fsize = c(0.7, 0.8), leg.txt = "Mean distance (dim 2)",
     outline = F, legend = 300)

plot(max_gen_contMap, fsize = c(0.7, 0.8), leg.txt = "Max distance (dim 2)",
     outline = F, legend = 300)

plot(sd_gen_contMap, fsize = c(0.7, 0.8), leg.txt = "SD distance (dim 2)",
     outline = F, legend = 300)

dev.off()

hyp_fam <- log(volume_hyper_fam$hv_vol[volume_hyper_fam$dimension == 2])
names(hyp_fam) <- volume_hyper_fam$family[volume_hyper_fam$dimension == 2]
hyp_fam <- hyp_fam[tr_fam_ladd$tip.label]

avg_fam <- log(volume_hyper_fam$avg_surv[volume_hyper_fam$dimension == 2])
names(avg_fam) <- volume_hyper_fam$family[volume_hyper_fam$dimension == 2]
avg_fam <- avg_fam[tr_fam_ladd$tip.label]

max_fam <- log(volume_hyper_fam$max_surv[volume_hyper_fam$dimension == 2])
names(max_fam) <- volume_hyper_fam$family[volume_hyper_fam$dimension == 2]
max_fam <- max_fam[tr_fam_ladd$tip.label]

sd_fam <- log(volume_hyper_fam$sd_surv[volume_hyper_fam$dimension == 2])
names(sd_fam) <- volume_hyper_fam$family[volume_hyper_fam$dimension == 2]
sd_fam <- sd_fam[tr_fam_ladd$tip.label]

hyp_fam_contMap <- contMap(tr_fam_ladd, hyp_fam, plot = F, res = 200)
hyp_fam_contMap <- setMap(hyp_fam_contMap, turbo(100))
hyp_fam_contMap$tree <- ladderize.simmap(hyp_fam_contMap$tree, right = F)

avg_fam_contMap <- contMap(tr_fam_ladd, avg_fam, plot = F, res = 200)
avg_fam_contMap <- setMap(avg_fam_contMap, turbo(100))
avg_fam_contMap$tree <- ladderize.simmap(avg_fam_contMap$tree, right = F)

max_fam_contMap <- contMap(tr_fam_ladd, max_fam, plot = F, res = 200)
max_fam_contMap <- setMap(max_fam_contMap, turbo(100))
max_fam_contMap$tree <- ladderize.simmap(max_fam_contMap$tree, right = F)

sd_fam_contMap <- contMap(tr_fam_ladd, sd_fam, plot = F, res = 200)
sd_fam_contMap <- setMap(sd_fam_contMap, turbo(100))
sd_fam_contMap$tree <- ladderize.simmap(sd_fam_contMap$tree, right = F)

pdf("figures/FigS3A-D.pdf", height = 5, width = 14)

layout(matrix(1:4, ncol = 4))

plot(hyp_fam_contMap, fsize = c(0.7, 0.8), leg.txt = "Hypervolume",
     outline = F, legend = 300)

plot(avg_fam_contMap, fsize = c(0.7, 0.8), leg.txt = "Mean distance (dim 2)",
     outline = F, legend = 300)

plot(max_fam_contMap, fsize = c(0.7, 0.8), leg.txt = "Max distance (dim 2)",
     outline = F, legend = 300)

plot(sd_fam_contMap, fsize = c(0.7, 0.8), leg.txt = "SD distance (dim 2)",
     outline = F, legend = 300)

dev.off()

## phylogenetic signal centroid
load(file = "data/env_df_gen_dist.RData")
load(file = "data/env_df_fam_dist.RData")

dist_cent_gen <- do.call(c, env_df_gen_dist$dist_cent)
names(dist_cent_gen) <- env_df_gen_dist$genus

dist_cent_fam <- do.call(c, env_df_fam_dist$dist_cent)
names(dist_cent_fam) <- env_df_fam_dist$family

cent_gen_sig_K <- phylosig(keep.tip(tr_gen, names(dist_cent_gen)),
                           log(dist_cent_gen), method = "K", test = T, 
                           nsim = 1000) 
cent_gen_sig_L <- phylosig(keep.tip(tr_gen, names(dist_cent_gen)),
                           log(dist_cent_gen), method = "lambda", test = T, 
                           nsim = 1000) 

cent_fam_sig_K <- phylosig(keep.tip(tr_fam, names(dist_cent_fam)),
                           log(dist_cent_fam), method = "K", test = T, 
                           nsim = 1000) 
cent_fam_sig_L <- phylosig(keep.tip(tr_fam, names(dist_cent_fam)),
                           log(dist_cent_fam), method = "lambda", test = T, 
                           nsim = 1000) 

cent_sig <- data.frame(
  level = c("genus", "family"),
  L = c(round(cent_gen_sig_L$lambda, 3), round(cent_fam_sig_L$lambda, 3)),
  logL0 = c(round(cent_gen_sig_L$logL0, 3), round(cent_fam_sig_L$logL0, 3)),
  logL = c(round(cent_gen_sig_L$logL, 3), round(cent_fam_sig_L$logL, 3)),
  p = c(round(cent_gen_sig_L$P, 3), round(cent_fam_sig_L$P, 3)),
  K = c(round(cent_gen_sig_K$K, 3), round(cent_fam_sig_K$K, 3)),
  pK = c(round(cent_gen_sig_K$P, 3), round(cent_fam_sig_K$P, 3)))

write.csv(cent_sig, "tables/TableS2_unformatted.csv")

plot_centroid_gen <- do.call(rbind, env_df_gen_dist$centroid)
plot_centroid_fam <- do.call(rbind, env_df_fam_dist$centroid)

mean_all <- c(1.150726e-13, -5.372683e-15, -3.000657e-15)

pdf("figures/FigS2D.pdf", height = 5, width = 7)

ggplot(plot_centroid_gen, aes(PC1, PC2, col = PC3)) +
  geom_point(size = 3) +
  geom_point(aes(x = mean_all[1], y = mean_all[2], col = mean_all[3]), shape = 8, size = 5) + 
  scale_color_viridis_c() +
  theme_bw(base_size = 19) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 9))

dev.off()

pdf("figures/FigS3E.pdf", height = 5, width = 7)

ggplot(plot_centroid_fam, aes(PC1, PC2, col = PC3)) +
  geom_point(size = 3) +
  geom_point(aes(x = mean_all[1], y = mean_all[2], col = mean_all[3]), shape = 8, size = 5) + 
  scale_color_viridis_c() +
  theme_bw(base_size = 19) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 9))

dev.off()
