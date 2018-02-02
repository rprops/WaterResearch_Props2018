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
library("flowClean")
library("nlme")
library("mgcv")
source("functions.R")
my.settings <- list(
  strip.background=list(col="transparent"),
  strip.border=list(col="transparent", cex=5),
  gate=list(col="black", fill="lightblue", alpha=0.2,border=NA,lwd=2),
  panel.background=list(col="lightgray"),
  background=list(col="white"))
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

set.seed(777)
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

summary <- fsApply(x=mixed_FCS, FUN=function(x) apply(x,2,max),use.exprs=TRUE)
mixed_FCS <- mixed_FCS[!is.infinite(summary[,1])]
max = max(summary[,1])
mytrans <- function(x) x/max
mixed_FCS <- transform(mixed_FCS ,`FL1-H`=mytrans(`FL1-H`),
                           `FL3-H`=mytrans(`FL3-H`), 
                           `SSC-H`=mytrans(`SSC-H`),
                           `FSC-H`=mytrans(`FSC-H`))

# Run phenotypic diversity analysis
diversity_stability <- Diversity_rf(stability_FCS, R = 100, param = param, d = 3)
diversity_bacteria <- Diversity_rf(bacteria_FCS, R = 100, param = param, d = 3)
diversity_mixed <- Diversity_rf(mixed_FCS, R = 100, param = param, d = 3)

# Create fingerprints for cluster analysis
fp_bacteria <- flowBasis(bacteria_FCS, param=param, nbin = 128, bw = 0.01, normalize = function(x) x)
fp_mixed <- flowBasis(mixed_FCS, param=param, nbin = 128, bw = 0.01, normalize = function(x) x)

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

# Extract cluster labels
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

mean_stab_HNA_b <- mean(results_bacteria$HNA_cells[results_bacteria$Time<16]/results_bacteria$Total_cells[results_bacteria$Time<16])
mean_stab_count_b <- mean(results_bacteria$Total_cells[results_bacteria$Time<16])
mean_stab_D2_b <- mean(results_bacteria$D2[results_bacteria$Time<16])
sd_stab_HNA_b <- sd(results_bacteria$HNA_cells[results_bacteria$Time<16]/results_bacteria$Total_cells[results_bacteria$Time<16])
sd_stab_count_b <- sd(results_bacteria$Total_cells[results_bacteria$Time<16])
sd_stab_D2_b <- sd(results_bacteria$D2[results_bacteria$Time<16])

mean_stab_HNA_m <- mean(results_mixed$HNA_cells[results_mixed$Time<16]/results_mixed$Total_cells[results_mixed$Time<16])
mean_stab_count_m <- mean(results_mixed$Total_cells[results_mixed$Time<16])
mean_stab_D2_m <- mean(results_mixed$D2[results_mixed$Time<16])
sd_stab_HNA_m <- sd(results_mixed$HNA_cells[results_mixed$Time<16]/results_mixed$Total_cells[results_mixed$Time<16])
sd_stab_count_m <- sd(results_mixed$Total_cells[results_mixed$Time<16])
sd_stab_D2_m <- sd(results_mixed$D2[results_mixed$Time<16])

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
  labs(y="", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0, size=18))+
  ggtitle(bquote("(B) Total cell density (cells µL"^{-1}*")") )+
  geom_line(color="black", alpha = 0.9)+
  xlim(0,max(results_stability$Time))+
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80))+
  scale_y_continuous(breaks=c(0, 25, 50, 75), limits = c(0,75))

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
  labs(y="", fill = "Deviation (s.d.)", x = "Time (min.)")+
  ggtitle("(D) % HNA cells")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), plot.title = element_text(hjust = 0, size=18))+
  geom_line(color="black", alpha = 0.9)+
  xlim(0,max(results_stability$Time))+
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80))+
  scale_y_continuous(breaks=c(0, 25, 50, 75), limits = c(0,75)) 

png("Fig2.png", res=500, height = 10, width = 10, units="in")
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
  ggtitle("Green fluorescence intensity (FL1-H)")+
  ylim(8,13)+
  guides(fill=FALSE)+
  xlim(0,max(results_bacteria$Time*60/0.1))+
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80)*60/0.1)

p_density_b <- ggplot(results_bacteria, aes(x = as.numeric(Time), y = Total_cells, fill = diff_count))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu", limits = c(0,3), breaks=c(0,1,2,3), oob=squish)+
  theme_bw()+
  ylim(0,225)+
  labs(y="", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0, size=18))+
  ggtitle(bquote("(A) Total cell concentration (cells µL"^{-1}*")") )+
  geom_line(color="black", alpha = 0.9)+
  xlim(0,max(results_bacteria$Time))+
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80))+
  guides(fill = FALSE)

p_HNA_b <-  ggplot(results_bacteria, aes(x = as.numeric(Time), y = 100*HNA_cells/Total_cells, fill = diff_HNA))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu", limits = c(0,3), breaks=c(0,1,2,3), oob=squish)+
  theme_bw()+
  ylim(0,100)+
  labs(y="", fill = "Deviation (s.d.)",x = "Time (min.)")+
  theme(axis.text.x = element_text(size=14),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0, size=18),
        axis.title=element_text(size=16))+
  ggtitle(bquote("(B) % HNA cells") )+
  geom_line(color="black", alpha = 0.9)+
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80))+
  # xlim(0,max(results_bacteria$Time))+
  guides(fill = FALSE)

p_diversity_b <-  ggplot(results_bacteria, aes(x = as.numeric(Time), y = D2, fill = diff_D2))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu", limits = c(0,3), breaks=c(0,1,2,3), oob=squish)+
  theme_bw()+
  ylim(900,3000)+
  labs(y="", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0, size=18))+
  ggtitle("(C) Phenotypic diversity index (a.u.)")+
  geom_line(color="black", alpha = 0.9)+
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.01)+
  xlim(0,max(results_bacteria$Time))+
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80))+
  guides(fill = FALSE)

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

png("Fig3.png", res=450, height = 8, width = 13, units="in")
ggarrange(p_density_b, p_diversity_b,
          p_HNA_b,  p_cluster_b,
          ncol=2)
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
  ylim(0,225)+
  labs(y="", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0, size=18))+
  ggtitle(bquote("(A) Total cell concentration (cells µL"^{-1}*")") )+
  geom_line(color="black", alpha = 0.9)+
  xlim(0, 80)+
  guides(fill = FALSE)

p_HNA_mixed <-  ggplot(results_mixed, aes(x = as.numeric(Time), y = 100*HNA_cells/Total_cells, fill = diff_HNA))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu", limits = c(0,3), breaks=c(0,1,2,3), oob=squish)+
  theme_bw()+
  ylim(0,100)+
  labs(y="", fill = "Deviation (s.d.)",x = "Time (min.)")+
  theme(axis.text.x = element_text(size=14),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0, size=18),
        axis.title=element_text(size=16))+
  ggtitle(bquote("(B) % HNA cells") )+
  geom_line(color="black", alpha = 0.9)+
  xlim(0, 80)+
  guides(fill = FALSE)

p_diversity_mixed <-  ggplot(results_mixed, aes(x = as.numeric(Time), y = D2, fill = diff_D2))+
  geom_point(shape=21, size = 4, alpha = 0.5)+
  scale_fill_distiller(palette="RdBu", limits = c(0,3), breaks=c(0,1,2,3), oob=squish)+
  theme_bw()+
  labs(y="", fill = "Deviation (s.d.)")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y=element_text(size=14),
        plot.title = element_text(hjust = 0, size=18))+
  ggtitle("(C) Phenotypic diversity index (a.u.)")+
  geom_line(color="black", alpha = 0.9)+
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.01)+
  xlim(0,80)+
  ylim(900,3000)+
  guides(fill = FALSE)

p_cluster_mixed <- ggplot(results_mixed, aes(x = as.numeric(Time), y = cluster_label))+
  geom_point(shape=21, size = 4, alpha = 0.5, aes(fill = factor(cluster_label)))+
  scale_fill_manual(values = c("#2166AC","#33A02C","#E31A1C","#FF7F00","#6A3D9A"))+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), plot.title = element_text(hjust = 0, size=18)
        #, panel.background = element_rect(fill = "lightblue",
        #                                 colour = "lightblue",
        #                                 size = 0.5, linetype = "solid"),
        # panel.grid.major = element_line(size = 0.5, linetype = 'solid',
        #                                 colour = "white"), 
        # panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
        #                                 colour = "white")
        )+
  labs(y="", x = "Time (min.)")+
  ggtitle("(D) Phenotypic community type")+
  geom_line(color="black", alpha = 0.9)+
  guides(fill = FALSE)+
  xlim(0,80)+
  scale_y_continuous(breaks=c(0:5), limits = c(0,5.5))

png("Fig4.png", res=450, height = 8, width = 13, units="in")
ggarrange(p_density_mixed, p_diversity_mixed, 
          p_HNA_mixed, p_cluster_mixed, 
          ncol=2)
dev.off()

# Phase 2: Evaluate temporal resolution required for robust estimate of
# cells / HNA / phenotypic diversity

# Import data
dirs <- list.dirs("binned_data", recursive = FALSE)
dirs <- dirs[grep(dirs, pattern = "Stability")]
dirs <- dirs[1:26]

# Merge all flowData
for(i in dirs){
  flowSum <- read.flowSet(path = i)
  param=c("FL1-H", "FL3-H", "SSC-H", "FSC-H")
  flowSum <- flowSum[,param]
  # Transform parameters
  flowSum <- transform(flowSum,`FL1-H`=asinh(`FL1-H`), `SSC-H`=asinh(`SSC-H`), 
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
  
  # Extract counts
  a <- flowCore::filter(flowSum, rGate_HNA)
  HNACount <- summary(a);HNACount <- toTable(HNACount)
  s <- flowCore::filter(flowSum, polyGate1)
  TotalCount <- summary(s);TotalCount <- toTable(TotalCount)
  vol <- c()
  for(j in 1:length(flowSum)){
    vol[j] <- as.numeric(flowSum[[j]]@description$`$VOL`)/1000
  }
  counts_flowSum <- data.frame(Sample = flowCore::sampleNames(flowSum),
                               Total_cells = TotalCount$true,
                               HNA_cells = HNACount$true,
                               LNA_cells = (TotalCount$true - HNACount$true),
                               volume = vol)
  
  # Normalize parameters
  summary <- fsApply(x=flowSum ,FUN=function(x) apply(x,2,max),use.exprs=TRUE)
  flowSum <- flowSum[!is.infinite(summary[,1])]
  max = max(summary[,1])
  mytrans <- function(x) x/17
  flowSum <- transform(flowSum ,`FL1-H`=mytrans(`FL1-H`),
                       `FL3-H`=mytrans(`FL3-H`), 
                       `SSC-H`=mytrans(`SSC-H`),
                       `FSC-H`=mytrans(`FSC-H`))
  
  # Run phenotypic diversity analysis
  diversity_flowSum <- Diversity_rf(flowSum, R = 3, R.b = 100, param = param, d = 3)
  
  # Extract metadata
  metadata_flowSum <- data.frame(Sample_names = flowCore::sampleNames(flowSum), 
                                 do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowSum),"_"), rbind)))[,1:5]
  colnames(metadata_flowSum)[2:5] <- c("Time", "Resolution", "Experiment", "Replicate")
  metadata_flowSum$Replicate <- gsub(metadata_flowSum$Replicate, pattern = ".fcs", replacement = "") 
  
  # Merge dataframes
  results_flowSum <- left_join(diversity_flowSum, counts_flowSum, by = c("Sample_names" = "Sample"))
  results_flowSum <- left_join(results_flowSum, metadata_flowSum, by = "Sample_names")
  
  if(i == dirs[1]) results_final <- results_flowSum else results_final <- rbind(results_final, results_flowSum)
}

# Filter out some binning sizes
results_final <- results_final[results_final$Resolution != '3600' & results_final$Resolution != '3000' & results_final$Resolution !='2400',]
results_final$Resolution <- droplevels(results_final$Resolution)

# Filter out outliers at 10 second resolution
results_final_river <- filter(results_final, Total_cells/volume > 100 & Replicate == "river")
results_final_tap <- filter(results_final, Replicate == "tap1")

# Order resolution factor level
results_final_river$Resolution <- factor(results_final_river$Resolution, 
                                         levels = c("10", "30", "60", "120", "180", "240", "300", "600", "1200", "1800"))
results_final_tap$Resolution <- factor(results_final_tap$Resolution, 
                                         levels = c("10", "30", "60", "120", "180", "240", "300", "600", "1200", "1800"))
# Calculate coefficient of variation for each temporal resolution

# remove outliers by replacing the ones larger than 1.5*IQR by the 95% and 5% quantiles // called capping
x <- results_final_river$Total_cells/results_final_river$volume
qnt <- quantile(x, probs=c(.25, .75), na.rm = T)
caps <- quantile(x, probs=c(.05, .95), na.rm = T)
H <- 1.5 * IQR(x, na.rm = T)
x[x < (qnt[1] - H)] <- caps[1]
x[x > (qnt[2] + H)] <- caps[2]

results_final_river$Total_cells <- x*results_final_river$volume
results_final_river <- results_final_river %>% 
  group_by(Resolution) %>% 
  mutate(CV_dens = 100*sd(Total_cells/volume)/mean(Total_cells/volume),
         CV_D2 = 100*sd(D2)/mean(D2))

CV_river <- do.call(rbind,by(results_final_river[,c(13,15,16,17)], INDICES = factor(results_final_river$Resolution), 
                 FUN = unique))

CV_river <- data.frame(CV_river, Total_cells = do.call(rbind,by(results_final_river[,c(8,9)], INDICES = factor(results_final_river$Resolution), 
                                                            FUN = colMeans))[,1])

results_final_tap <- results_final_tap %>% 
  group_by(Resolution) %>% 
  mutate(CV_dens = 100*sd(Total_cells/volume)/mean(Total_cells/volume),
         CV_D2 = 100*sd(D2)/mean(D2))

CV_tap <- do.call(rbind,by(results_final_tap[,c(13,15,16,17)], INDICES = factor(results_final_tap$Resolution), 
                             FUN = unique))
CV_tap <- data.frame(CV_tap, Total_cells = do.call(rbind,by(results_final_tap[,c(8,9)], INDICES = factor(results_final_tap$Resolution), 
                                                       FUN = colMeans))[,1])

# Combine both
CV_total <- rbind(CV_river, CV_tap)

# Make plots
p_stab_density_tap <- ggplot(data = results_final_tap, aes(x = Total_cells, y = Total_cells/volume, fill = Resolution))+
  geom_point(shape = 21, size = 4, alpha = 0.5)+
  scale_fill_brewer(palette = "Paired")+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), plot.title = element_text(hjust = 0, size=18),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),legend.title=element_text(size=16))+
  # ylim(10,50)+
  geom_smooth(fill="gray", color = "black", span = 1)+
  annotation_logticks(sides = "b")+
  labs(y = bquote("Total cell density (cells µL"^{-1}*")"), x = "Cells measured", fill = "Resolution (s)")+
  ggtitle("(A)")

p_stab_diversity_tap <- ggplot(data = results_final_tap, aes(x = Total_cells, y = D2, fill = Resolution))+
  geom_point(shape = 21, size = 4, alpha = 0.5)+
  scale_fill_brewer(palette = "Paired")+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), plot.title = element_text(hjust = 0, size=18),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),legend.title=element_text(size=16))+  ylim(1750,2750)+
  geom_smooth(fill="gray", color = "black", span = 1)+
  annotation_logticks(sides = "b")+
  labs(y = bquote("Phenotypic diversity (a.u.)"), x = "Cells measured", fill = "Resolution (s)")+
  ggtitle("(B)")

p_stab_density_river <- ggplot(data = results_final_river, aes(x = Total_cells, y = Total_cells/volume, fill = Resolution))+
  geom_point(shape = 21, size = 4, alpha = 0.5)+
  scale_fill_brewer(palette = "Paired")+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), plot.title = element_text(hjust = 0, size=18),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),legend.title=element_text(size=16))+  # ylim(10,50)+
  geom_smooth(fill="gray", color = "black", span = 1)+
  annotation_logticks(sides = "b")+
  labs(y = bquote("Total cell density (cells µL"^{-1}*")"), x = "Cells measured", fill = "Resolution (s)")+
  ggtitle("(C)")+
  ylim(200,500)

p_stab_diversity_river <- ggplot(data = results_final_river, aes(x = Total_cells, y = D2, fill = Resolution))+
  geom_point(shape = 21, size = 4, alpha = 0.5)+
  scale_fill_brewer(palette = "Paired")+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), plot.title = element_text(hjust = 0, size=18),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),legend.title=element_text(size=16))+
  ylim(3000,4000)+
  geom_smooth(fill="gray", color = "black", span = 1)+
  annotation_logticks(sides = "b")+
  labs(y = bquote("Phenotypic diversity (a.u.)"), x = "Cells measured", fill = "Resolution (s)")+
  ggtitle("(D)")


png("Fig5.png", res=500, height = 10, width = 10, units="in")
grid_arrange_shared_legend(p_stab_density_tap, p_stab_diversity_tap, p_stab_density_river, p_stab_diversity_river, ncol = 2, nrow = 2)
dev.off()

# Plots for coefficient of variation
CV_total$Replicate <- gsub(CV_total$Replicate, pattern = "tap1", replacement = "Tap water")
CV_total$Replicate <- gsub(CV_total$Replicate, pattern = "river", replacement = "River water")

p_CV_D2 <- ggplot(data = CV_total, aes(x = Total_cells, y = CV_D2, fill = Resolution))+
  geom_point(size = 4, aes(fill = Resolution, shape = Replicate), alpha = 0.5)+
  scale_fill_brewer(palette = "Paired")+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw()+
  scale_shape_manual(values=c(21,22))+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14), plot.title = element_text(hjust = 0, size=18),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),legend.title=element_text(size=16))+
  # geom_smooth(fill="gray", color = "black", span = 1,formula=y~x)+
  geom_smooth(method="lm",color="black",fill="gray",formula=y~x)+
  annotation_logticks(sides = "b")+
  labs(y = "CV on phenotypic diversity (%)", x = "Cells measured", fill = "Resolution (s)", shape = "")+
  ggtitle("(A)")+
  guides(fill = guide_legend(override.aes = list(shape = 21)), ncol=1)+
  ylim(0,5)+
  geom_hline(yintercept = 5, linetype = 2)


p_CV_dens <- ggplot(data = CV_total, aes(x = Total_cells, y = CV_dens, fill = Resolution))+
  geom_point(size = 4, aes(fill = Resolution,  shape = Replicate), alpha = 0.5)+
  scale_fill_brewer(palette = "Paired")+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw()+
  scale_shape_manual(values=c(21,22))+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14), plot.title = element_text(hjust = 0, size=18),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=15),legend.title=element_text(size=16))+
  # geom_smooth(fill="gray", color = "black", span = 1,formula=y~x)+
  geom_smooth(color="black",fill="gray",formula=y~x)+
  annotation_logticks(sides = "b")+
  labs(y = "CV on cell density (%)", x = "Cells measured", fill = "Resolution (s)", shape = "")+
  ggtitle("(B)")+
  guides(fill = guide_legend(override.aes = list(shape = 21)), ncol=1)+
  ylim(0,15)+
  geom_hline(yintercept = 5, linetype = 2)

png("Fig6_run2.png", res=500, height = 5, width = 10, units="in")
grid_arrange_shared_legend(p_CV_D2, p_CV_dens, ncol = 2, nrow = 1)
dev.off()

# Autocorrelation check and gam analysis for tap water
results_final_tap_60 <- results_stability

### GAMs for stability inference on phenotypic diversity index 

gam.D2 <- gamm(D2~s(Time, k=60), data=results_final_tap_60, 
               correlation=corCAR1(form =~Time, value = 0.5))
plot(gam.D2$gam,residuals=TRUE)
plot(residuals.gam(gam.D2$gam), type = "l")
acf(residuals.gam(gam.D2$gam))

# Check for normality in model residuals
car::qqPlot(residuals.gam(gam.D2$gam))

plot(D2~Time, data = results_final_tap_60, ylim = c(1000,3000))
points(predict(gam.D2$gam), x = results_final_tap_60$Time, col ="red")

# Check for significant autocorrelation
Box.test(residuals.gam(gam.D2$gam), lag=15, type="Ljung-Box")

# Check for significance of slope
anova(gam.D2$gam)

### GAMs for stability inference on total cell density 
gam.TC <- gamm(Total_cells~s(Time, k=60), data=results_final_tap_60, 
               correlation=corCAR1(form =~Time, value = 0.5))
plot(gam.TC$gam,residuals=TRUE)
plot(residuals.gam(gam.TC$gam), type = "l")
acf(residuals.gam(gam.TC$gam))

qqPlot(residuals.gam(gam.TC$gam))


plot(Total_cells~Time, data = results_final_tap_60)
points(predict(gam.TC$gam), x = results_final_tap_60$Time, col ="red")

# Check for significant autocorrelation
Box.test(residuals.gam(gam.TC$gam), lag=15, type="Ljung-Box")

# Check for significance of slope
anova(gam.TC$gam)


# Autocorrelation check and gam analysis for river water
results_final_river_60 <- results_final_river[results_final_river$Resolution ==60,]
results_final_river_60$Time <- as.numeric(as.character(results_final_river_60$Time))

### GAMs for stability inference on phenotypic diversity index 

gam.D2 <- gamm(D2~s(Time, k=60), data=results_final_river_60, 
               correlation=corCAR1(form =~Time, value = 0.5))
plot(gam.D2$gam,residuals=TRUE)
plot(residuals.gam(gam.D2$gam), type = "l")
acf(residuals.gam(gam.D2$gam))

# Check for normality in model residuals
qqPlot(residuals.gam(gam.D2$gam))

plot(D2~Time, data = results_final_river_60, ylim = c(1000,5000))
points(predict(gam.D2$gam), x = results_final_river_60$Time, col ="red")

# Check for significant autocorrelation
Box.test(residuals.gam(gam.D2$gam), lag=15, type="Ljung-Box")

# Check for significance of slope
anova(gam.D2$gam)

### GAMs for stability inference on total cell density 
gam.TC <- gamm(Total_cells/volume~s(Time, k=60), data=results_final_river_60, 
               correlation=corCAR1(form =~Time, value = 0.5))
plot(gam.TC$gam,residuals=TRUE)
plot(residuals.gam(gam.TC$gam), type = "l")

# Check for normality in model residuals
qqPlot(residuals.gam(gam.TC$gam))

plot(Total_cells/volume~Time, data = results_final_river_60, ylim = c(0,320))
points(predict(gam.TC$gam), x = results_final_river_60$Time, col ="red")

# Check for significant autocorrelation
acf(residuals.gam(gam.TC$gam))
Box.test(residuals.gam(gam.TC$gam), lag=15, type="Ljung-Box")

# Check for significance of slope
anova(gam.TC$gam)



### Fig S1

fs_S1 <- transform(original_data[1],`FL1-H`=asinh(`FL1-H`), `SSC-H`=asinh(`SSC-H`), 
                                 `FL3-H`=asinh(`FL3-H`), `FSC-H`=asinh(`FSC-H`))

sqrcut1 <- matrix(c(asinh(20000),asinh(20000),16,16,
                    3,9.55,16,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
rGate_HNA <- polygonGate(.gate=sqrcut1, filterId = "HNA")
sqrcut1 <- matrix(c(8.5,8.5,asinh(20000),asinh(20000),
                    3,7.25,9.55,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
rGate_LNA <- polygonGate(.gate=sqrcut1, filterId = "LNA")


sqrcut1 <- matrix(c(8.5,8.5,16,16,
                    3,7.25,16,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")
### Creating a rectangle gate, set correct threshold here for FL1
sqrcut2 <- matrix(c(asinh(20000),asinh(20000),20,20,
                    0,20,20,0),ncol=2, nrow=4)


filters <- filters(list(rGate_LNA, rGate_HNA))
flist <- list(filters)
names(flist) <- flowCore::sampleNames(fs_S1)

png("FigS1.png", res=500, height = 5, width = 6, units="in")
print(xyplot(`FL3-H`~`FL1-H`, data=fs_S1,
             filter=flist,
             xbins=500,nbin=128, par.strip.text=list(col="black", font=3,cex=1), 
             smooth=FALSE, xlim = c(7.5,15), ylim = c(2.5,15), xlab=list(label="Green fluorescence intensity (FL1-H)",cex=1.35),ylab=list(label="Red fluorescence intensity (FL3-H)",cex=1.35),
             par.settings=my.settings,
             scales=list(x=list(at=seq(from=0, to=15, by=2.5),cex=1),
                         y=list(at=seq(from=0, to=15, by=2.5),cex=1)), layout=c(1,1),
             margin=TRUE,
             binTrans="log",
             strip=strip.custom(factor.levels="")
)
)
dev.off()


### Calculate coefficient of variation for each bin in each experiment and visualize
### this fingerprint for Figure 6
cv_fp_bacteria <- apply(fp_bacteria@basis[, 1:16384], 2, function(x) sd(x))
cv_fp_mixed <- apply(fp_mixed@basis[, 1:16384], 2, function(x) sd(x))
dummy_fp <- fp_bacteria@basis[1, 1:16384]
dummy_fp[dummy_fp > 0.00000001] <- NA

nbin <- 128
Y <- c()
for (i in 1:nbin) Y <- c(Y, rep(i, 128))
X = rep(1:nbin, nbin)
df_cv_bacteria <- data.frame(dev = cv_fp_bacteria, X = rep(1:nbin, nbin), 
           Y = Y)
df_cv_mixed <- data.frame(dev = cv_fp_mixed, X = rep(1:nbin, nbin), 
                             Y = Y)
df_dummy <- data.frame(dev = dummy_fp, X = rep(1:nbin, nbin), 
                       Y = Y)
colnames(df_cv_mixed)[2:3] <- c("FL1-H", "FL3-H")
colnames(df_cv_bacteria)[2:3] <- c("FL1-H", "FL3-H")
colnames(df_dummy)[2:3] <- c("FL1-H", "FL3-H")
df_dummy <- df_dummy[is.na(df_dummy$dev), ]

# Filter out low SD values
df_cv_bacteria <- df_cv_bacteria %>% dplyr::filter(dev > 5)
df_cv_mixed <- df_cv_mixed %>% filter(dev > 5)

v_bact <- ggplot2::ggplot(df_cv_bacteria, ggplot2::aes(`FL1-H`, `FL3-H`, z = dev))+
  ggplot2::geom_point(data = df_dummy, colour="gray",size = 1, shape = 22)+
  ggplot2::geom_tile(ggplot2::aes(fill=dev))+
  ggplot2::scale_fill_distiller(palette="RdBu", na.value="white", trans = "log",
                                labels = seq(15, 60, 15), breaks = seq(15, 60, 15)) + 
  # ggplot2::stat_contour(ggplot2::aes(fill=..level..), geom="polygon", binwidth=0.1)+
  ggplot2::theme_bw()+
  # ggplot2::geom_contour(color = "white", alpha = 1)+
  ylim(25,100)+
  xlim(40,110)+
  theme(axis.title=element_text(size=16), strip.text.x=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2))
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )

print(v_bact)


v_mixed <- ggplot2::ggplot(df_cv_mixed, ggplot2::aes(`FL1-H`, `FL3-H`, z = dev))+
  ggplot2::geom_tile(ggplot2::aes(fill=dev)) +
  # ggplot2::geom_point(size = 1, ggplot2::aes(fill=cor), shape =22,
  # color = "transparent")+
  ggplot2::scale_fill_distiller(palette="RdBu", na.value="white", trans = "log",
                                labels = seq(15, 60, 15), breaks = seq(15, 60, 15)) + 
  # ggplot2::stat_contour(ggplot2::aes(fill=..level..), geom="polygon", binwidth=0.1)+
  ggplot2::theme_bw()+
  # ggplot2::geom_contour(color = "white", alpha = 1)+
  ylim(0,100)+
  xlim(40,110)+
  theme(axis.title=element_text(size=16), strip.text.x=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2))
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )


print(v_mixed)

### Extract labels of samples belonging to different "robust" community types
### For natural community contaminations
lab_mixed_c0 <- results_mixed$cluster_label == 1 & results_mixed$Time < 10
lab_mixed_c1 <- results_mixed$cluster_label == 2
lab_mixed_c2 <- results_mixed$cluster_label == 3
lab_mixed_c3 <- results_mixed$cluster_label == 4
lab_mixed_c4 <- results_mixed$cluster_label == 5


### Calculate contrasts between control and contaminations for all 5 contaminations
fp_mixed_c1 <- fp_contrasts(fp_mixed, comp2 = lab_mixed_c0, 
                            comp1 = lab_mixed_c1)
fp_mixed_c2 <- fp_contrasts(fp_mixed, comp2 = lab_mixed_c0, 
                            comp1 = lab_mixed_c2)
fp_mixed_c3 <- fp_contrasts(fp_mixed, comp2 = lab_mixed_c0, 
                            comp1 = lab_mixed_c3)
fp_mixed_c4 <- fp_contrasts(fp_mixed, comp2 = lab_mixed_c0, 
                            comp1 = lab_mixed_c4)

### Merge all contrasts
fp_merged <- data.frame(rbind(fp_mixed_c1, fp_mixed_c2,
                              fp_mixed_c3, fp_mixed_c4), 
                        Contamination = c(rep("Evian", nrow(fp_mixed_c1)),
                                            rep("Pond", nrow(fp_mixed_c2)),
                                            rep("River-1", nrow(fp_mixed_c3)),
                                            rep("River-2", nrow(fp_mixed_c4))))

### Visualize different contrasts
v_c_all <- ggplot(fp_merged, aes(`FL1.H`, `FL3.H`, z = Density))+
  geom_tile(aes(fill=Density)) + 
  geom_point(colour="gray", alpha=0.7)+
  scale_fill_distiller(expression(Delta~"Density"), palette="RdYlBu", na.value="white", limits = c(-0.65, 0.65)) + 
  stat_contour(aes(fill=..level..), geom="polygon", binwidth=0.1)+
  theme_bw()+
  facet_grid(~Contamination)+
  geom_contour(color = "white", alpha = 1)+
  labs(x="Green fluorescence intensity (a.u.)", y="Red fluorescence intensity (a.u.)")+
  theme(axis.title=element_text(size=16), strip.text.x=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2))
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )

png("Fig6a.png", res=500, height = 5, width = 10, units="in")
print(v_c_all)
dev.off()

### Visualize different contrasts
# idx <- duplicated(interaction(fp_merged$FL1.H, fp_merged$FL3.H)) | duplicated(interaction(fp_merged$FL1.H, fp_merged$FL3.H), fromLast = TRUE)
# fp_merged$Contamination <- as.character(fp_merged$Contamination)
# fp_merged$Contamination[idx] <- "Shared"


v_c_cat <- fp_merged %>% dplyr::filter(Density>=0) %>% 
  ggplot(aes(`FL1.H`, `FL3.H`, z = Density))+
  geom_tile(aes(fill=Contamination), alpha=0.7) + 
  geom_point(colour="gray", alpha=0.4)+
  # scale_fill_distiller(palette="RdBu", na.value="white") + 
  scale_fill_manual(values = c("#33A02C","#E31A1C","#FF7F00","#6A3D9A", "Grey"))+
  # stat_contour(aes(fill=..level..), geom="polygon", binwidth=0.1)+
  theme_bw()+
  # geom_contour(color = "white", alpha = 1)+
  labs(x="Green fluorescence intensity (a.u.)", y="Red fluorescence intensity (a.u.)")+
  theme(axis.title=element_text(size=16), strip.text.x=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2))
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )

png("Fig6b.png", res=500, height = 5, width = 7, units="in")
print(v_c_cat)
dev.off()

