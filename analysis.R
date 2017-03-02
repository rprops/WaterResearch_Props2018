library("Phenoflow")
library("dplyr")
library("ggplot2")
library("grid")
library("gridExtra")
library("egg")
library("qdapRegex")
library("scales")

# # Bin the data
# 
# # Import original data
# original_data <- read.flowSet(path = "original_data")
# 
# # We then extract the total analysis time for each sample in seconds
# analysis_length <- c()
# for(i in 1:length(original_data)) analysis_length[i] <- as.numeric(original_data[[i]]@description$`#ACQUISITIONTIMEMILLI`)/1000
# analysis_length <- data.frame(time=analysis_length); rownames(analysis_length) <- flowCore::sampleNames(original_data)
# 
# # Choose the size of your time-gate (in seconds)
# time_interval <- c(10, 30, 60, 120, 180, 240, 300, 600,
#                    1200, 1800, 2400, 3000, 3600)
# time.step = 0.1
# 
# # Bin the data and export into new directories
# setwd("binned_data/")
# for(i in time_interval){
#   time_discretization(x = original_data, analysis.length = analysis_length, create=TRUE,
#                                               time.interval = i, height = c(0,200000000),
#                                             start = 0, time.step = 0.1)
# }
# setwd("..")


# Import binned data
stability_FCS <- read.flowSet(path = "binned_data/Stability_tap1_60")
bacteria_FCS <- read.flowSet(path = "binned_data/Bacteria_Run_60")
mixed_FCS <- read.flowSet(path = "binned_data/Mixed_Run_60")

# Transform parameters by asinh transformation and select primary parameters of interest
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")
stability_FCS <- stability_FCS[,param]
bacteria_FCS <- bacteria_FCS[,param]
mixed_FCS <- mixed_FCS[,param]
stability_FCS <- transform(stability_FCS,`FL1-H`=asinh(`FL1-H`), `SSC-H`=asinh(`SSC-H`), 
                           `FL3-H`=asinh(`FL3-H`), `FSC-H`=asinh(`FSC-H`))
bacteria_FCS <- transform(bacteria_FCS,`FL1-H`=asinh(`FL1-H`), `SSC-H`=asinh(`SSC-H`), 
                           `FL3-H`=asinh(`FL3-H`), `FSC-H`=asinh(`FSC-H`))
mixed_FCS <- transform(mixed_FCS,`FL1-H`=asinh(`FL1-H`), `SSC-H`=asinh(`SSC-H`), 
                           `FL3-H`=asinh(`FL3-H`), `FSC-H`=asinh(`FSC-H`))


# Create a PolygonGate for extracting the single-cell information
sqrcut1 <- matrix(c(8.75,8.75,15,15,3,8,15,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")
### Creating a rectangle gate, set correct threshold here for FL1
sqrcut2 <- matrix(c(asinh(20000),asinh(20000),20,20,
                    0,20,20,0),ncol=2, nrow=4)
colnames(sqrcut2) <- c("FL1-H","FL3-H")
rGate_HNA <- polygonGate(.gate=sqrcut2, filterId = "HNA bacteria")

# Check if gate is correct for all data
xyplot(`FL3-H` ~ `FL1-H`, data=stability_FCS[1], filter=polyGate1,
       scales=list(y=list(limits=c(0,15)),
                   x=list(limits=c(6,15))),
       axis = axis.default, nbin=125, par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

xyplot(`FL3-H` ~ `FL1-H`, data=bacteria_FCS[1], filter=polyGate1,
       scales=list(y=list(limits=c(0,15)),
                   x=list(limits=c(6,15))),
       axis = axis.default, nbin=125, par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

xyplot(`FL3-H` ~ `FL1-H`, data=mixed_FCS[1], filter=polyGate1,
       scales=list(y=list(limits=c(0,15)),
                   x=list(limits=c(6,15))),
       axis = axis.default, nbin=125, par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

# Proceed only with the bacterial information for fingerprinting
stability_FCS <- Subset(stability_FCS, polyGate1)
bacteria_FCS <- Subset(bacteria_FCS, polyGate1)
mixed_FCS <- Subset(mixed_FCS, polyGate1)

# Extract the cell counts
# 1. Stability data
a <- flowCore::filter(stability_FCS, rGate_HNA)
HNACount <- summary(a);HNACount <- toTable(HNACount)
s <- flowCore::filter(stability_FCS, polyGate1)
TotalCount <- summary(s);TotalCount <- toTable(TotalCount)
vol <- c()
for(i in 1:length(stability_FCS)) vol[i] <- as.numeric(stability_FCS[[i]]@description$`$VOL`)/1000
counts_stability <- data.frame(Sample = flowCore::sampleNames(stability_FCS),
                               Total_cells = TotalCount$true/vol,
                               HNA_cells = HNACount$true/vol,
                               LNA_cells = (TotalCount$true - HNACount$true)/vol)
# 2. Bacteria data
a <- flowCore::filter(bacteria_FCS, rGate_HNA)
HNACount <- summary(a);HNACount <- toTable(HNACount)
s <- flowCore::filter(bacteria_FCS, polyGate1)
TotalCount <- summary(s);TotalCount <- toTable(TotalCount)
vol <- c()
for(i in 1:length(bacteria_FCS)) vol[i] <- as.numeric(bacteria_FCS[[i]]@description$`$VOL`)/1000
counts_bacteria <- data.frame(Sample = flowCore::sampleNames(bacteria_FCS),
                               Total_cells = TotalCount$true/vol,
                               HNA_cells = HNACount$true/vol,
                               LNA_cells = (TotalCount$true - HNACount$true)/vol)
# 3. Mixed data
a <- flowCore::filter(mixed_FCS, rGate_HNA)
HNACount <- summary(a);HNACount <- toTable(HNACount)
s <- flowCore::filter(mixed_FCS, polyGate1)
TotalCount <- summary(s);TotalCount <- toTable(TotalCount)
vol <- c()
for(i in 1:length(mixed_FCS)) vol[i] <- as.numeric(mixed_FCS[[i]]@description$`$VOL`)/1000
counts_mixed <- data.frame(Sample = flowCore::sampleNames(mixed_FCS),
                               Total_cells = TotalCount$true/vol,
                               HNA_cells = HNACount$true/vol,
                               LNA_cells = (TotalCount$true - HNACount$true)/vol)

# Normalize the signals based on the maximum FL1-H fluorescence
summary <- fsApply(x=stability_FCS ,FUN=function(x) apply(x,2,max),use.exprs=TRUE)
stability_FCS <- stability_FCS[!is.infinite(summary[,1])]
summary <- summary[!is.infinite(summary[,1]),]
max = max(summary[,1])
mytrans <- function(x) x/max
stability_FCS <- transform(stability_FCS ,`FL1-H`=mytrans(`FL1-H`),
                                  `FL3-H`=mytrans(`FL3-H`), 
                                  `SSC-H`=mytrans(`SSC-H`),
                                  `FSC-H`=mytrans(`FSC-H`))

summary <- fsApply(x=bacteria_FCS ,FUN=function(x) apply(x,2,max),use.exprs=TRUE)
bacteria_FCS <- bacteria_FCS[!is.infinite(summary[,1])]
summary <- summary[!is.infinite(summary[,1]),]
max = max(summary[,1])
mytrans <- function(x) x/max
bacteria_FCS <- transform(bacteria_FCS ,`FL1-H`=mytrans(`FL1-H`),
                           `FL3-H`=mytrans(`FL3-H`), 
                           `SSC-H`=mytrans(`SSC-H`),
                           `FSC-H`=mytrans(`FSC-H`))

summary <- fsApply(x=mixed_FCS ,FUN=function(x) apply(x,2,max),use.exprs=TRUE)
mixed_FCS <- mixed_FCS[!is.infinite(summary[,1])]
max = max(summary[,1])
mytrans <- function(x) x/max
mixed_FCS <- transform(mixed_FCS ,`FL1-H`=mytrans(`FL1-H`),
                           `FL3-H`=mytrans(`FL3-H`), 
                           `SSC-H`=mytrans(`SSC-H`),
                           `FSC-H`=mytrans(`FSC-H`))

# Run phenotypic diversity analysis
diversity_stability <- Diversity_rf(stability_FCS, R = 3, param = param, d = 3)
diversity_bacteria <- Diversity_rf(bacteria_FCS, R = 3, param = param, d = 3)
diversity_mixed <- Diversity_rf(mixed_FCS, R = 3, param = param, d = 3)

# Merge count and phenotypic diversity data in one file
results_stability <- left_join(counts_stability, diversity_stability, by = c("Sample" = "Sample_names"))
results_stability <- results_stability[results_stability$Total_cells > 0, ]
results_bacteria <- left_join(counts_bacteria, diversity_bacteria, by = c("Sample" = "Sample_names"))
results_bacteria <- results_bacteria[results_bacteria$Total_cells > 0, ]
results_mixed <- left_join(counts_mixed, diversity_mixed, by = c("Sample" = "Sample_names"))
results_mixed <- results_mixed[results_mixed$Total_cells > 0, ]

# Add time points
meta_stability <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(stability_FCS),"_"), rbind)))
meta_bacteria <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(bacteria_FCS),"_"), rbind)))
meta_mixed <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(mixed_FCS),"_"), rbind)))
results_stability <- data.frame(results_stability, Time = as.numeric(as.character(meta_stability$X1)))
results_bacteria <- data.frame(results_bacteria, Time = as.numeric(as.character(meta_bacteria$X1)))
results_mixed <- data.frame(results_mixed, Time = as.numeric(as.character(meta_mixed$X1)))

# Calculate mean and sd of HNA, count and diversity
mean_stab_HNA <- mean(results_stability$HNA_cells)
mean_stab_count <- mean(results_stability$Total_cells)
mean_stab_D2 <- mean(results_stability$D2)
sd_stab_HNA <- sd(results_stability$HNA_cells)
sd_stab_count <- sd(results_stability$Total_cells)
sd_stab_D2 <- sd(results_stability$D2)

# Add column to results indicating difference between stability run and observed values
results_stability <- data.frame(results_stability, diff_HNA = abs(results_stability$HNA_cells-mean_stab_HNA)/sd_stab_HNA, 
                                diff_count = abs(results_stability$Total_cells-mean_stab_count)/sd_stab_count, 
                                diff_D2 = abs(results_stability$D2-mean_stab_D2)/sd_stab_D2)
results_bacteria <- data.frame(results_bacteria, diff_HNA = abs(results_bacteria$HNA_cells-mean_stab_HNA)/sd_stab_HNA, 
                               diff_count = abs(results_bacteria$Total_cells-mean_stab_count)/sd_stab_count, 
                               diff_D2 = abs(results_bacteria$D2-mean_stab_D2)/sd_stab_D2)
results_mixed <- data.frame(results_mixed, diff_HNA = abs(results_mixed$HNA_cells-mean_stab_HNA)/sd_stab_HNA, 
                            diff_count = abs(results_mixed$Total_cells-mean_stab_count)/sd_stab_count, 
                            diff_D2 = abs(results_mixed$D2-mean_stab_D2)/sd_stab_D2)
# Import clustered data
# These were generated on a 64x64 fingerprint by performing PCA on the bins
clusters_bacteria <- read.csv("clustering_data/Bacteria_run_silhouette.csv", stringsAsFactors = FALSE)
clusters_mixed <- read.csv("clustering_data/Mixed_run_silhouette.csv", stringsAsFactors = FALSE)

# formate names correctly
clusters_bacteria$X <- rm_between(clusters_bacteria$X, "'", "'", extract=TRUE)
clusters_mixed$X <- rm_between(clusters_mixed$X, "'", "'", extract=TRUE)
order_bacteria <- as.numeric(gsub(clusters_bacteria$X, pattern = "_.*", replacement=""))
order_mixed <- as.numeric(gsub(clusters_mixed$X, pattern = "_.*", replacement=""))
clusters_mixed$X <- gsub(clusters_mixed$X, pattern = "mixed", replacement="Mixed_Run")
clusters_bacteria$X <- gsub(clusters_bacteria$X, pattern = "Run1", replacement="Run")

# Sort rows of cluster data
clusters_mixed <- clusters_mixed[order(order_mixed),]
clusters_bacteria <- clusters_bacteria[order(order_bacteria),]

# Replace sample numbers by correct ones
clusters_bacteria$X <- paste(seq(2:nrow(clusters_bacteria)),"_60_", 
                             gsub(clusters_bacteria$X, pattern = ".*_60_", replacement=""), sep="")
clusters_mixed$X <- paste(seq(2:nrow(clusters_mixed)),"_60_", 
                             gsub(clusters_mixed$X, pattern = ".*_60_", replacement=""), sep="")
# Only take cluster allocation and sample name
clusters_bacteria <- data.frame(Sample = clusters_bacteria$X, cluster_alloc = clusters_bacteria$Cluster.prediction)
clusters_mixed <- data.frame(Sample = clusters_mixed$X, cluster_alloc = clusters_mixed$Cluster.prediction)

# Merge data with other results
results_bacteria <- left_join(results_bacteria, clusters_bacteria, by = "Sample")
results_mixed <- left_join(results_mixed, clusters_mixed, by = "Sample")

# Create stability plots
original_data <- read.flowSet(path = "original_data")
FL1_stability <- data.frame(FL1 = asinh(exprs(original_data[[4]])[,9]),
                            Time = exprs(original_data[[4]])[,14])

p_FL1 <- ggplot(FL1_stability, aes(x = Time, y = FL1))+
  geom_density_2d(n=100, alpha = 0.5)+ 
  stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE, n = 100)+
  ggplot2::scale_fill_distiller(palette="RdBu", na.value="white",
                                trans = "sqrt")+
  labs(y = "Green fluorescence intensity (FL1-H)")+
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  ylim(8,13)+
  guides(fill=FALSE)
  
p_density <-  ggplot(results_stability, aes(x = as.numeric(Time), y = Total_cells, fill = diff_count))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu")+
  theme_bw()+
  ylim(0,100)+
  labs(y="Total cell density (cells/µL)", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  geom_line(color="black", alpha = 0.9)

p_HNA <-  ggplot(results_stability, aes(x = as.numeric(Time), y = 100*HNA_cells/Total_cells, fill = diff_HNA))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu")+
  theme_bw()+
  ylim(0,100)+
  labs(y="% HNA cells", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  geom_line(color="black", alpha = 0.9)

p_diversity <-  ggplot(results_stability, aes(x = as.numeric(Time), y = D2, fill = diff_D2))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu")+
  theme_bw()+
  ylim(1000,3350)+
  labs(y="Phenotypic diversity", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  geom_line(color="black", alpha = 0.9)+
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.01)

g1 <- ggplotGrob(p_FL1)
g2 <- ggplotGrob(p_density)
g3 <- ggplotGrob(p_HNA)
g4 <- ggplotGrob(p_diversity)
fg2 <- gtable_frame(g2)
fg3 <- gtable_frame(g3)
fg4 <- gtable_frame(g4)
fg234 <- gtable_frame(rbind(fg2, fg3, fg4))
fg1 <- gtable_frame(g1)
grid.newpage()
combined <- rbind(fg1, fg234)
grid.draw(combined)


# Create bacteria contamination plots
FL1_bacteria <- data.frame(FL1 = asinh(exprs(original_data[[1]])[,9]),
                            Time = exprs(original_data[[1]])[,14])

p_FL1 <- ggplot(FL1_bacteria, aes(x = Time, y = FL1))+
  geom_density_2d(n=100, alpha = 0.5)+ 
  stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE, n = 100)+
  ggplot2::scale_fill_distiller(palette="RdBu", na.value="white",
                                trans = "sqrt")+
  labs(y = "Green fluorescence intensity (FL1-H)")+
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  ylim(8,13)+
  guides(fill=FALSE)

p_density <-  ggplot(results_bacteria, aes(x = as.numeric(Time), y = Total_cells, fill = diff_count))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu", limits = c(1,6), breaks=c(1,3,5), oob=squish)+
  theme_bw()+
  ylim(0,300)+
  labs(y="Total cell density (cells/µL)", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  geom_line(color="black", alpha = 0.9)

# p_HNA <-  ggplot(results_bacteria, aes(x = as.numeric(Time), y = 100*HNA_cells/Total_cells, fill = diff_HNA))+
#   geom_point(shape=21, size = 4, alpha = 0.5)+
#   scale_fill_distiller(palette="RdBu")+
#   theme_bw()+
#   ylim(0,100)+
#   labs(y="% HNA cells", fill = "Deviation (s.d.)")+
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
#   geom_line(color="black", alpha = 0.9)

p_diversity <-  ggplot(results_bacteria, aes(x = as.numeric(Time), y = D2, fill = diff_D2))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu", limits = c(1,6), breaks=c(1,3,5), oob=squish)+
  theme_bw()+
  ylim(1000,3350)+
  labs(y="Phenotypic diversity", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  geom_line(color="black", alpha = 0.9)+
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.01)

p_cluster <-  ggplot(results_bacteria, aes(x = as.numeric(Time), y = cluster_alloc))+
  geom_point(shape=21, size = 4, alpha = 0.5, aes(fill = factor(cluster_alloc)))+
  scale_fill_manual(values = c("blue","red","red","red"))+
  theme_bw()+
  labs(y="Cluster allocation", x = "Time (min.)")+
  ylim(0,3)+
  geom_line(color="black", alpha = 0.9)+
  guides(fill = FALSE)

g1 <- ggplotGrob(p_FL1)
g2 <- ggplotGrob(p_density)
g3 <- ggplotGrob(p_diversity)
g4 <- ggplotGrob(p_cluster)
fg2 <- gtable_frame(g2)
fg3 <- gtable_frame(g3)
fg4 <- gtable_frame(g4)
fg234 <- gtable_frame(rbind(fg2, fg3, fg4))
fg1 <- gtable_frame(g1)
grid.newpage()
combined <- rbind(fg1, fg234)
grid.draw(combined)


# Create mixed contamination plots
FL1_mixed <- data.frame(FL1 = asinh(exprs(original_data[[2]])[,9]),
                           Time = exprs(original_data[[2]])[,14])

p_FL1 <- ggplot(FL1_mixed, aes(x = Time, y = FL1))+
  geom_density_2d(n=100, alpha = 0.5)+ 
  stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE, n = 100)+
  ggplot2::scale_fill_distiller(palette="RdBu", na.value="white",
                                trans = "sqrt")+
  labs(y = "Green fluorescence intensity (FL1-H)", x = "")+
  theme_bw()+ 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  ylim(8,13)+
  guides(fill=FALSE)

p_density <-  ggplot(results_mixed, aes(x = as.numeric(Time), y = Total_cells, fill = diff_count))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu", limits = c(1,3), breaks=c(1,2,3), oob=squish)+
  theme_bw()+
  ylim(0,600)+
  labs(y="Total cell density (cells/µL)", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  geom_line(color="black", alpha = 0.9)

# p_HNA <-  ggplot(results_bacteria, aes(x = as.numeric(Time), y = 100*HNA_cells/Total_cells, fill = diff_HNA))+
#   geom_point(shape=21, size = 4, alpha = 0.5)+
#   scale_fill_distiller(palette="RdBu")+
#   theme_bw()+
#   ylim(0,100)+
#   labs(y="% HNA cells", fill = "Deviation (s.d.)")+
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
#   geom_line(color="black", alpha = 0.9)

p_diversity <-  ggplot(results_mixed, aes(x = as.numeric(Time), y = D2, fill = diff_D2))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu", limits = c(1,3), breaks=c(1,2,3), oob=squish)+
  theme_bw()+
  ylim(1000,3350)+
  labs(y="Phenotypic diversity", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  geom_line(color="black", alpha = 0.9)+
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.01)

p_cluster <-  ggplot(results_mixed, aes(x = as.numeric(Time), y = cluster_alloc))+
  geom_point(shape=21, size = 4, alpha = 0.5, aes(fill = factor(cluster_alloc)))+
  scale_fill_manual(values = c("blue","red","red","red","red"))+
  theme_bw()+
  labs(y="Cluster allocation", x = "Time (min.)")+
  ylim(0,6)+
  geom_line(color="black", alpha = 0.9)+
  guides(fill = FALSE)

g1 <- ggplotGrob(p_FL1)
g2 <- ggplotGrob(p_density)
g3 <- ggplotGrob(p_diversity)
g4 <- ggplotGrob(p_cluster)
fg2 <- gtable_frame(g2)
fg3 <- gtable_frame(g3)
fg4 <- gtable_frame(g4)
fg234 <- gtable_frame(rbind(fg2, fg3, fg4))
fg1 <- gtable_frame(g1, height = unit(0.5,"null"))
grid.newpage()
combined <- rbind(fg1, fg234)
png("test_fig.png", res=500, height = 10, width = 6, units="in")
grid.arrange(combined)
dev.off()
