library("Phenoflow")
library("dplyr")
library("ggplot2")
library("grid")
library("gridExtra")
library("egg")
library("qdapRegex")
library("scales")
library("cluster")
library("factoextra")
source("functions.R")
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

# Remove last sample as this is a null observation (no or little cells).
stability_FCS <- stability_FCS[-58]
bacteria_FCS <- bacteria_FCS[-86]
mixed_FCS <- mixed_FCS[-74]

# Transform parameters by asinh transformation and select primary parameters of interest
param=c("FL1-H", "FL3-H", "SSC-H", "FSC-H")
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
sqrcut1 <- matrix(c(8.5,8.5,16,16,3,7.25,16,3),ncol=2, nrow=4)
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

# Create fingerprints for cluster analysis
fp_bacteria <- flowBasis(bacteria_FCS, param=param, nbin = 128, bw = 0.01, normalize = function(x) x)
fp_mixed <- flowBasis(mixed_FCS, param=param, nbin = 128, bw = 0.01, normalize = function(x) x)

# For exporting features
# write.table(t(fp_bacteria@basis),file="bacteria_60_128.txt",col.names=FALSE)
# write.table(t(fp_mixed@basis),file="mixed_60_128.txt",col.names=FALSE)
# write.table(rownames(fp_bacteria@basis), file="bacteria_names_60.txt",row.names=FALSE,col.names=FALSE)
# write.table(rownames(fp_mixed@basis), file="mixed_names_60.txt",row.names=FALSE,col.names=FALSE)

# Perform PCA to reduce number of features in fingerprint
pca_bacteria <- prcomp(fp_bacteria@basis)
pca_mixed <- prcomp(fp_mixed@basis)

# Only retain PC which explain x% of the variance
thresh <- 0.9
nr_pc_bacteria <- min(which((cumsum(vegan::eigenvals(pca_bacteria)/sum(vegan::eigenvals(pca_bacteria)))>thresh)==TRUE))
nr_pc_mixed <- min(which((cumsum(vegan::eigenvals(pca_mixed)/sum(vegan::eigenvals(pca_mixed)))>thresh)==TRUE))

pc_cluster_bacteria <- pca_bacteria$x[, 1:nr_pc_bacteria]
pc_cluster_mixed <- pca_mixed$x[, 1:nr_pc_mixed]

# Evaluate number of robust clusters by means of silhouette index

tmp.si <- c()
for(i in 2:(nrow(pc_cluster_bacteria)-1)){
  tmp.si[i] <- pam(pc_cluster_bacteria, k=i)$silinfo$avg.width
}
nr_clusters_bacteria <- which(tmp.si == max(tmp.si, na.rm = TRUE))

tmp.si <- c()
for(i in 2:(nrow(pc_cluster_mixed)-1)){
  tmp.si[i] <- pam(pc_cluster_mixed, k=i)$silinfo$avg.width
}
nr_clusters_mixed <- which(tmp.si == max(tmp.si, na.rm = TRUE))


# Cluster samples and export cluster labels
clusters_bacteria <- pam(pc_cluster_bacteria, k=nr_clusters_bacteria)
clusters_mixed <- pam(pc_cluster_mixed, k=nr_clusters_mixed)

# # Extract cluster labels
cluster_labels_bacteria <- data.frame(Sample = names(clusters_bacteria$clustering),
                                      cluster_label = clusters_bacteria$clustering)
cluster_labels_mixed <- data.frame(Sample = names(clusters_mixed$clustering),
                                      cluster_label = clusters_mixed$clustering)

# Merge count and phenotypic diversity data in one file
results_stability <- left_join(counts_stability, diversity_stability, by = c("Sample" = "Sample_names"))
results_stability <- results_stability[results_stability$Total_cells > 0, ]
results_bacteria <- left_join(counts_bacteria, diversity_bacteria, by = c("Sample" = "Sample_names"))
results_bacteria <- results_bacteria[results_bacteria$Total_cells > 0, ]
results_mixed <- left_join(counts_mixed, diversity_mixed, by = c("Sample" = "Sample_names"))
results_mixed <- results_mixed[results_mixed$Total_cells > 0, ]

# Merge results with cluster labels
results_bacteria <- left_join(results_bacteria, cluster_labels_bacteria, by = "Sample")
results_mixed <- left_join(results_mixed, cluster_labels_mixed, by = "Sample")

# Add time points
meta_stability <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(stability_FCS),"_"), rbind)))
meta_bacteria <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(bacteria_FCS),"_"), rbind)))
meta_mixed <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(mixed_FCS),"_"), rbind)))
results_stability <- data.frame(results_stability, Time = as.numeric(as.character(meta_stability$X1)))
results_bacteria <- data.frame(results_bacteria, Time = as.numeric(as.character(meta_bacteria$X1)))
results_mixed <- data.frame(results_mixed, Time = as.numeric(as.character(meta_mixed$X1)))

# Calculate mean and sd of HNA, count and diversity based on first 15 minutes of run time
# For the bacteria and mixed contaminations. Consider the entire run for the stability
# experiment.
mean_stab_HNA_s <- mean(results_stability$HNA_cells/results_stability$Total_cells)
mean_stab_count_s <- mean(results_stability$Total_cells)
mean_stab_D2_s <- mean(results_stability$D2)
sd_stab_HNA_s <- sd(results_stability$HNA_cells/results_stability$Total_cells)
sd_stab_count_s <- sd(results_stability$Total_cells)
sd_stab_D2_s <- sd(results_stability$D2)

mean_stab_HNA_b <- mean(results_bacteria$HNA_cells[1:15]/results_bacteria$Total_cells[1:15])
mean_stab_count_b <- mean(results_bacteria$Total_cells[1:15])
mean_stab_D2_b <- mean(results_bacteria$D2[1:15])
sd_stab_HNA_b <- sd(results_bacteria$HNA_cells[1:15]/results_bacteria$Total_cells[1:15])
sd_stab_count_b <- sd(results_bacteria$Total_cells[1:15])
sd_stab_D2_b <- sd(results_bacteria$D2[1:15])

mean_stab_HNA_m <- mean(results_mixed$HNA_cells[1:15]/results_mixed$Total_cells[1:15])
mean_stab_count_m <- mean(results_mixed$Total_cells[1:15])
mean_stab_D2_m <- mean(results_mixed$D2[1:15])
sd_stab_HNA_m <- sd(results_mixed$HNA_cells[1:15]/results_mixed$Total_cells[1:15])
sd_stab_count_m <- sd(results_mixed$Total_cells[1:15])
sd_stab_D2_m <- sd(results_mixed$D2[1:15])

# Add column to results indicating difference between stability run and observed values
results_stability <- data.frame(results_stability, diff_HNA = abs(results_stability$HNA_cells/results_stability$Total_cells-mean_stab_HNA_s)/sd_stab_HNA_s, 
                            diff_count = abs(results_stability$Total_cells-mean_stab_count_s)/sd_stab_count_s, 
                            diff_D2 = abs(results_stability$D2-mean_stab_D2_s)/sd_stab_D2_s)
results_bacteria <- data.frame(results_bacteria, diff_HNA = abs(results_bacteria$HNA_cells/results_bacteria$Total_cells-mean_stab_HNA_b)/sd_stab_HNA_b, 
                               diff_count = abs(results_bacteria$Total_cells-mean_stab_count_b)/sd_stab_count_b, 
                               diff_D2 = abs(results_bacteria$D2-mean_stab_D2_b)/sd_stab_D2_b)
results_mixed <- data.frame(results_mixed, diff_HNA = abs(results_mixed$HNA_cells/results_mixed$Total_cells-mean_stab_HNA_m)/sd_stab_HNA_m, 
                            diff_count = abs(results_mixed$Total_cells-mean_stab_count_m)/sd_stab_count_m, 
                            diff_D2 = abs(results_mixed$D2-mean_stab_D2_m)/sd_stab_D2_m)
# # Import clustered data
# # These were generated on a 64x64 fingerprint by performing PCA on the bins
# clusters_bacteria <- read.csv("clustering_data/Bacteria_run_silhouette.csv", stringsAsFactors = FALSE)
# clusters_mixed <- read.csv("clustering_data/Mixed_run_silhouette.csv", stringsAsFactors = FALSE)
# 
# # formate names correctly
# clusters_bacteria$X <- rm_between(clusters_bacteria$X, "'", "'", extract=TRUE)
# clusters_mixed$X <- rm_between(clusters_mixed$X, "'", "'", extract=TRUE)
# order_bacteria <- as.numeric(gsub(clusters_bacteria$X, pattern = "_.*", replacement=""))
# order_mixed <- as.numeric(gsub(clusters_mixed$X, pattern = "_.*", replacement=""))
# clusters_mixed$X <- gsub(clusters_mixed$X, pattern = "mixed", replacement="Mixed_Run")
# clusters_bacteria$X <- gsub(clusters_bacteria$X, pattern = "Run1", replacement="Run")
# 
# # Sort rows of cluster data
# clusters_mixed <- clusters_mixed[order(order_mixed),]
# clusters_bacteria <- clusters_bacteria[order(order_bacteria),]
# 
# # Replace sample numbers by correct ones
# clusters_bacteria$X <- paste(seq(2:nrow(clusters_bacteria)),"_60_", 
#                              gsub(clusters_bacteria$X, pattern = ".*_60_", replacement=""), sep="")
# clusters_mixed$X <- paste(seq(2:nrow(clusters_mixed)),"_60_", 
#                              gsub(clusters_mixed$X, pattern = ".*_60_", replacement=""), sep="")
# # Only take cluster allocation and sample name
# clusters_bacteria <- data.frame(Sample = clusters_bacteria$X, cluster_alloc = clusters_bacteria$Cluster.prediction)
# clusters_mixed <- data.frame(Sample = clusters_mixed$X, cluster_alloc = clusters_mixed$Cluster.prediction)
# 
# # Merge data with other results
# results_bacteria <- left_join(results_bacteria, clusters_bacteria, by = "Sample")
# results_mixed <- left_join(results_mixed, clusters_mixed, by = "Sample")

# Create stability plots
original_data <- read.flowSet(path = "original_data")
FL1_stability <- data.frame(FL1 = asinh(exprs(original_data[[4]])[,9]),
                            Time = exprs(original_data[[4]])[,14])

p_FL1 <- ggplot(FL1_stability, aes(x = Time, y = FL1))+
  stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE, n = 100)+
  ggplot2::scale_fill_distiller(palette="Greens", na.value="white",
                                direction = 1, type ="div")+
  labs(y = "")+
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0,size=18))+
  ggtitle("(A) Green fluorescence intensity (FL1-H)")+
  ylim(8,13)+
  guides(fill=FALSE)+
  xlim(0,max(results_stability$Time*60/0.1))+
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80)*60/0.1)

p_density <- ggplot(results_stability, aes(x = as.numeric(Time), y = Total_cells, fill = diff_count))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu", limits = c(0,3), breaks=c(0,1,2,3), oob=squish)+
  theme_bw()+
  ylim(0,100)+
  labs(y="", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0, size=18))+
  ggtitle(bquote("(B) Total cell density (cells µL"^{-1}*")") )+
  geom_line(color="black", alpha = 0.9)+
  xlim(0,max(results_stability$Time))+
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80))

p_diversity <-  ggplot(results_stability, aes(x = as.numeric(Time), y = D2, fill = diff_D2))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu", limits = c(0,3), breaks=c(0,1,2,3), oob=squish)+
  theme_bw()+
  ylim(900,3350)+
  labs(y="", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text=element_text(size=14),
        plot.title = element_text(hjust = 0, size=18))+
  labs(y="", x = "Time (min.)")+
  ggtitle("(C) Phenotypic diversity index (a.u.)")+
  geom_line(color="black", alpha = 0.9)+
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.01)+
  xlim(0,max(results_stability$Time))+
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80))

p_HNA <- ggplot(results_stability, aes(x = as.numeric(Time), y = 100*HNA_cells/Total_cells, fill = diff_HNA))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu", limits = c(0,3), breaks=c(0,1,2,3), oob=squish)+
  theme_bw()+
  ylim(0,100)+
  labs(y="", fill = "Deviation (s.d.)", x = "Time (min.)")+
  ggtitle("(D) % HNA cells")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), plot.title = element_text(hjust = 0, size=18))+
  geom_line(color="black", alpha = 0.9)+
  xlim(0,max(results_stability$Time))+
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80))

png("Fig1.png", res=500, height = 10, width = 10, units="in")
ggarrange(p_FL1, p_density, p_diversity, p_HNA, ncol=1)
dev.off()

# Create bacteria contamination plots
FL1_bacteria <- data.frame(FL1 = asinh(exprs(original_data[[1]])[,9]),
                            Time = exprs(original_data[[1]])[,14])

p_FL1_b <- ggplot(FL1_bacteria, aes(x = Time, y = FL1))+
  stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE, n = 100)+
  ggplot2::scale_fill_distiller(palette="Greens", na.value="white",
                                direction = 1, type ="div")+
  labs(y = "")+
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0,size=18))+
  ggtitle("(A) Green fluorescence intensity (FL1-H)")+
  ylim(8,13)+
  guides(fill=FALSE)+
  xlim(0,max(results_bacteria$Time*60/0.1))+
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80)*60/0.1)

p_density_b <- ggplot(results_bacteria, aes(x = as.numeric(Time), y = Total_cells, fill = diff_count))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu", limits = c(0,3), breaks=c(0,1,2,3), oob=squish)+
  theme_bw()+
  ylim(0,250)+
  labs(y="", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0, size=18))+
  ggtitle(bquote("(B) Total cell density (cells µL"^{-1}*")") )+
  geom_line(color="black", alpha = 0.9)+
  xlim(0,max(results_bacteria$Time))+
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80))

# p_HNA_b <-  ggplot(results_bacteria, aes(x = as.numeric(Time), y = 100*HNA_cells/Total_cells, fill = diff_HNA))+
#   geom_point(shape=21, size = 4, alpha = 0.5)+
#   scale_fill_distiller(palette="RdBu")+
#   theme_bw()+
#   ylim(0,100)+
#   labs(y="% HNA cells", fill = "Deviation (s.d.)")+
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
#   geom_line(color="black", alpha = 0.9)+
#   xlim(0,max(results_bacteria$Time))

p_diversity_b <-  ggplot(results_bacteria, aes(x = as.numeric(Time), y = D2, fill = diff_D2))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu", limits = c(0,3), breaks=c(0,1,2,3), oob=squish)+
  theme_bw()+
  ylim(900,3350)+
  labs(y="", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0, size=18))+
  ggtitle("(C) Phenotypic diversity index (a.u.)")+
  geom_line(color="black", alpha = 0.9)+
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.01)+
  xlim(0,max(results_bacteria$Time))+
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80))

p_cluster_b <- ggplot(results_bacteria, aes(x = as.numeric(Time), y = cluster_label))+
  geom_point(shape=21, size = 4, alpha = 0.5, aes(fill = factor(cluster_label)))+
  scale_fill_manual(values = c("#2166AC","#33A02C","#6A3D9A"))+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), plot.title = element_text(hjust = 0, size=18))+
  labs(y="", x = "Time (min.)")+
  ggtitle("(D) Phenotypic community type")+
  ylim(0,3.5)+
  geom_line(color="black", alpha = 0.9)+
  guides(fill = FALSE)+
  xlim(0,max(results_bacteria$Time))+
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80))

png("Fig2.png", res=500, height = 10, width = 10, units="in")
ggarrange(p_FL1_b, p_density_b, p_diversity_b, p_cluster_b, ncol=1)
dev.off()

# Create mixed contamination plots
FL1_mixed <- data.frame(FL1 = asinh(exprs(original_data[[2]])[,9]),
                           Time = exprs(original_data[[2]])[,14])

p_FL1_mixed <- ggplot(FL1_mixed, aes(x = Time, y = FL1))+
  stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE, n = 100)+
  ggplot2::scale_fill_distiller(palette="Greens", na.value="white",
                               direction = 1, type ="div")+
  labs(y = "")+
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0,size=18))+
  ggtitle("(A) Green fluorescence intensity (FL1-H)")+
  ylim(8, 13)+
  guides(fill=FALSE)+
  xlim(0, 80*60/0.1)+
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80)*60/0.1)

p_density_mixed <-  ggplot(results_mixed, aes(x = as.numeric(Time), y = Total_cells, fill = diff_count))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu", limits = c(0,3), breaks=c(0,1,2,3), oob=squish)+
  theme_bw()+
  ylim(0,250)+
  labs(y="", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0, size=18))+
  ggtitle(bquote("(B) Total cell density (cells µL"^{-1}*")") )+
  geom_line(color="black", alpha = 0.9)+
  xlim(0, 80)

# p_HNA_mixed <-  ggplot(results_mixed, aes(x = as.numeric(Time), y = 100*HNA_cells/Total_cells, fill = diff_HNA))+
#   geom_point(shape=21, size = 4, alpha = 0.5)+
#   scale_fill_distiller(palette="RdBu")+
#   theme_bw()+
#   ylim(0,100)+
#   labs(y="% HNA cells", fill = "Deviation (s.d.)")+
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
#   geom_line(color="black", alpha = 0.9)+
#   xlim(0,max(results_mixed$Time))

p_diversity_mixed <-  ggplot(results_mixed, aes(x = as.numeric(Time), y = D2, fill = diff_D2))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu", limits = c(0,3), breaks=c(0,1,2,3), oob=squish)+
  theme_bw()+
  ylim(900,3350)+
  labs(y="", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0, size=18))+
  ggtitle("(C) Phenotypic diversity index (a.u.)")+
  geom_line(color="black", alpha = 0.9)+
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.01)+
  xlim(0,80)

p_cluster_mixed <- ggplot(results_mixed, aes(x = as.numeric(Time), y = cluster_label))+
  geom_point(shape=21, size = 4, alpha = 0.5, aes(fill = factor(cluster_label)))+
  scale_fill_manual(values = c("#2166AC","#33A02C","#E31A1C","#FF7F00","#6A3D9A"))+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), plot.title = element_text(hjust = 0, size=18))+
  labs(y="", x = "Time (min.)")+
  ggtitle("(D) Phenotypic community type")+
  ylim(0,5.5)+
  geom_line(color="black", alpha = 0.9)+
  guides(fill = FALSE)+
  xlim(0,80)

png("Fig3.png", res=500, height = 10, width = 10, units="in")
ggarrange(p_FL1_mixed, p_density_mixed, p_diversity_mixed, p_cluster_mixed, ncol=1)
dev.off()