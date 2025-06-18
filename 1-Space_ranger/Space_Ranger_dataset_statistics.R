.libPaths(c("/home/christophe.lepriol/NeuroDev_ADD/R/r_3.6.0", .libPaths()))
library(ggplot2)
library(ggforce)
library(ggpubr)
library(RColorBrewer)


##############
# Parameters #
##############
work_dir <- "/home/christophe.lepriol/NeuroDev_ADD/spatial_transcriptomics/projects/30-EpiReg"

# reference genomes
## one single reference genome
# genome_name <- "NCBIRefSeq108_Ensembl105GTF"
# genomes <- c("NCBIRefSeq108_Ensembl105GTF")
## reference genome comparison
genomes <- c("NCBIRefSeq106_Ensembl104GTF", "NCBIRefSeq106_NCBIRefSeq106GTF", "NCBIRefSeq108_NCBIRefSeq108GTF", "NCBIRefSeq108_Ensembl105GTF")
# samples
samples <- c("A_L1_S1", "A_L2_S5", "A_L3_S9", "A_L4_S13", "B_L1_S2", "B_L2_S6", "B_L3_S10", "B_L4_S14", "C_L1_S3", "C_L2_S7", "C_L3_S11", "C_L4_S15", "D_L1_S4", "D_L2_S8", "D_L3_S12", "D_L4_S16")
design_df <- data.frame(sample=samples, reset=rep(c(rep("none", 2), rep("reset", 2)), 4), condition=rep(c(rep("pilo", 4), rep("NaCl", 4)), 2), time=factor(c(rep(c(5, 10), 4), rep(c(20, 40, 40, 20), 2))))
# config filename
trimming_filtering_config_basename <- "10-TSO_polyA_R1hardtrim1_ov5_n2_min20"
# output directory
## one single reference genome
# output_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/%s/%s/10-Dataset", work_dir, trimming_filtering_config_basename, genome_name)
## reference genome comparison
output_dir <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/%s/10-Reference_comparison/10-Dataset", work_dir, trimming_filtering_config_basename)

plot_theme <- theme_bw() +
  theme(panel.border=element_rect(color="grey50")) +
  theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.title=element_text(size=12)) +
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=10))

# plot colors
plot_colors <- c("Ensembl 104" = brewer.pal(12, "Paired")[5], 
                 "Ensembl 105" = brewer.pal(12, "Paired")[6], 
                 "RefSeq 106" = brewer.pal(12, "Paired")[1], 
                 "RefSeq 108" = brewer.pal(12, "Paired")[2])

############
# Analysis #
############
# statistics data frame
all_samples_metrics_summary_df <- data.frame()
genome_vector <- c()
sample_vector <- c()
for (genome_name in genomes) {
  for (one_sample in samples) {
    # Space Ranger metrics summary csv file
    space_ranger_metrics_summary_file <- sprintf("%s/10-ST_analysis/10-Space_Ranger/output/%s/%s/00-Samples/%s/10-Pipeline/outs/metrics_summary.csv", work_dir, trimming_filtering_config_basename, genome_name, one_sample)
    space_ranger_metrics_summary_df <- read.csv(space_ranger_metrics_summary_file, quote="")
    all_samples_metrics_summary_df <- rbind(all_samples_metrics_summary_df, space_ranger_metrics_summary_df)
    genome_vector <- c(genome_vector, genome_name)
    sample_vector <- c(sample_vector, one_sample)
  }
}
## remove 'Sample.ID' column
all_samples_metrics_summary_df <- all_samples_metrics_summary_df[, ! colnames(all_samples_metrics_summary_df) == "Sample.ID"]
all_samples_metrics_summary_df <- cbind(genome=genome_vector, sample=sample_vector, all_samples_metrics_summary_df)
write.csv(all_samples_metrics_summary_df, file=sprintf("%s/all_samples_Space_Ranger_metrics_summary.csv", output_dir), quote=FALSE, row.names=FALSE)

# build data frame for ggplot
sample_vector <- time_vector <- condition_vector <- reset_vector <- genome_vector <- stat_vector <- value_vector <- c()
for (genome_name in genomes) {
  for (one_sample in samples) {
    sample_vector <- c(sample_vector, rep(one_sample, length(colnames(all_samples_metrics_summary_df))-2))
    time_vector <- c(time_vector, rep(as.character(design_df[which(design_df$sample==one_sample), "time"]), length(colnames(all_samples_metrics_summary_df))-2))
    condition_vector <- c(condition_vector, rep(as.character(design_df[which(design_df$sample==one_sample), "condition"]), length(colnames(all_samples_metrics_summary_df))-2))
    reset_vector <- c(reset_vector, rep(as.character(design_df[which(design_df$sample==one_sample), "reset"]), length(colnames(all_samples_metrics_summary_df))-2))
    genome_vector <- c(genome_vector, rep(genome_name, length(colnames(all_samples_metrics_summary_df))-2))
    stat_vector <- c(stat_vector, colnames(all_samples_metrics_summary_df)[3:length(colnames(all_samples_metrics_summary_df))])
    value_vector <- c(value_vector, as.numeric(all_samples_metrics_summary_df[which(all_samples_metrics_summary_df$genome==genome_name & all_samples_metrics_summary_df$sample==one_sample), 3:length(colnames(all_samples_metrics_summary_df))]))
  }
}
df2ggplot <- data.frame(sample=sample_vector, time=time_vector, condition=condition_vector, reset=reset_vector, genome=genome_vector, stat=stat_vector, value=value_vector)
levels(df2ggplot$genome)[levels(df2ggplot$genome) == "NCBIRefSeq106_Ensembl104GTF"] <- "Ensembl 104"
levels(df2ggplot$genome)[levels(df2ggplot$genome) == "NCBIRefSeq106_NCBIRefSeq106GTF"] <- "RefSeq 106"
levels(df2ggplot$genome)[levels(df2ggplot$genome) == "NCBIRefSeq108_NCBIRefSeq108GTF"] <- "RefSeq 108"
levels(df2ggplot$genome)[levels(df2ggplot$genome) == "NCBIRefSeq108_Ensembl105GTF"] <- "Ensembl 105"
write.csv(df2ggplot, file=sprintf("%s/all_samples_Space_Ranger_metrics_summary_df2ggplot.csv", output_dir), quote=FALSE, row.names=FALSE)
#### aggregated mean statistics
df2ggplot_mean <- aggregate(df2ggplot$value, by=list(genome=df2ggplot$genome, stat=df2ggplot$stat), FUN=mean)
colnames(df2ggplot_mean)[which(colnames(df2ggplot_mean) == "x")] <- "mean_value"
write.csv(df2ggplot_mean, file=sprintf("%s/all_samples_Space_Ranger_metrics_summary_df2ggplot_mean.csv", output_dir), quote=FALSE, row.names=FALSE)

# plots
## statistics boxplots and barplots
print("statistics boxplots and barplots")
### all stats
print("all stats")
pdf(sprintf("%s/all_samples_Space_Ranger_metrics_summary.pdf", output_dir))
for (one_stat in colnames(all_samples_metrics_summary_df)[3:length(colnames(all_samples_metrics_summary_df))]) {
  p <- ggplot(df2ggplot[which(df2ggplot$stat==one_stat),], aes(x=genome, y=value, fill=genome)) +
    geom_boxplot() +
    scale_fill_manual(values=plot_colors) +
    labs(title=sprintf("%s", gsub("\\.", " ", one_stat)), x="Reference genome", y="Value", fill="Reference genome") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.text.x=element_text(angle=30, hjust=1)) +
    theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.title=element_text(size=8)) +
    theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=6))
  print(p)
  
  p <- ggplot(df2ggplot[which(df2ggplot$stat==one_stat),], aes(x=sample, y=value, fill=genome)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=plot_colors) +
    labs(title=sprintf("%s", gsub("\\.", " ", one_stat)), x="Sample", y="Value", fill="Reference genome") +
    theme_bw() +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.title=element_text(size=8)) +
    theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=6))
  print(p)
}
dev.off()

# reset effect analysis
## order sample levels according to time, condition and reset
reset_samples <- c()
for (one_time in levels(design_df$time)) {
  for (one_condition in levels(design_df$condition)) {
    reset_samples <- c(reset_samples, c(as.character(design_df[which(design_df$time==one_time & design_df$condition==one_condition & design_df$reset=="none"), "sample"]), as.character(design_df[which(design_df$time==one_time & design_df$condition==one_condition & design_df$reset=="reset"), "sample"])))
  }
}
df2ggplot$sample <- factor(df2ggplot$sample, levels=reset_samples)
df2ggplot$time <- factor(df2ggplot$time, levels=c("5", "10", "20", "40"))
levels(df2ggplot$time) <- sprintf("J%s", levels(df2ggplot$time))
levels(df2ggplot$reset) <- c("sans", "avec")


### all stats
print("all stats")
levels(df2ggplot$genome)[levels(df2ggplot$genome) == "NCBIRefSeq106_Ensembl104GTF"] <- "Ensembl 104"
levels(df2ggplot$genome)[levels(df2ggplot$genome) == "NCBIRefSeq106_NCBIRefSeq106GTF"] <- "RefSeq 106"
levels(df2ggplot$genome)[levels(df2ggplot$genome) == "NCBIRefSeq108_NCBIRefSeq108GTF"] <- "RefSeq 108"
levels(df2ggplot$genome)[levels(df2ggplot$genome) == "NCBIRefSeq108_Ensembl105GTF"] <- "Ensembl 105"
pdf(sprintf("%s/all_samples_Space_Ranger_metrics_summary_reset.pdf", output_dir))
for (one_stat in colnames(all_samples_metrics_summary_df)[3:length(colnames(all_samples_metrics_summary_df))]) {
  one_stat_df2ggplot <- df2ggplot[which(df2ggplot$stat==one_stat),]
  one_stat_df2ggplot <- one_stat_df2ggplot[order(one_stat_df2ggplot$reset, one_stat_df2ggplot$time, one_stat_df2ggplot$condition),]
  
  p <- ggplot(one_stat_df2ggplot, aes(x=sample, y=value, fill=reset)) +
    geom_bar(stat="identity", position=position_dodge()) +
    labs(title=sprintf("%s", gsub("\\.", " ", one_stat)), x="Echantillon", y="Valeur", fill="Reset") +
    theme_bw() +
    plot_theme +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  print(p)
  
  p <- ggplot(one_stat_df2ggplot, aes(x=time, y=value, fill=condition, group=sample)) +
    geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8)) +
    labs(title=sprintf("%s", gsub("\\.", " ", one_stat)), x="Echantillon", y="Valeur", fill="Condition") +
    theme_bw() +
    plot_theme
  print(p)
  
  none_values <- one_stat_df2ggplot[which(one_stat_df2ggplot$reset=="sans"), "value"]
  reset_values <- one_stat_df2ggplot[which(one_stat_df2ggplot$reset=="avec"), "value"]
  differences <- none_values - reset_values
  shapiro_wilk_test_result <- shapiro.test(differences)
  wilcoxon_test_result <- wilcox.test(none_values, reset_values, paired=TRUE)
  t_test_result <- t.test(value ~ reset, data=one_stat_df2ggplot, paired=TRUE)
  # --> identical to: t.test(none_values, reset_values, paired=TRUE)
  p <- ggplot(one_stat_df2ggplot, aes(x=reset, y=value, fill=reset)) +
    geom_boxplot() +
    labs(title=sprintf("%s\nShapiro-Wilk test p-value=%.2f\nWilcoxon signed-rank test p-value=%.2f\npaired t-test p-value=%.2f", gsub("\\.", " ", one_stat), shapiro_wilk_test_result$p.value, wilcoxon_test_result$p.value, t_test_result$p.value), x="Reset", y="Valeur") +
    theme_bw() +
    plot_theme +
    theme(legend.position="none")
  print(p)
  
  p <- ggpaired(one_stat_df2ggplot, x="reset", y="value", color="reset", line.color="gray", line.size=0.4) +
    stat_compare_means(paired=TRUE) +
    labs(title=sprintf("%s", gsub("\\.", " ", one_stat)), x="Reset", y="Valeur") +
    theme(legend.position="none")
  print(p)
  
}
dev.off()


stats2plot <- c("Fraction.of.Spots.Under.Tissue", "Mean.Reads.Under.Tissue.per.Spot", "Median.Genes.per.Spot", "Total.Genes.Detected", "Median.UMI.Counts.per.Spot")
pdf(sprintf("%s/all_samples_Space_Ranger_metrics_summary_reset_2.pdf", output_dir))
for (one_stat in stats2plot) {
  one_stat_df2ggplot <- df2ggplot[which(df2ggplot$stat==one_stat),]
  one_stat_df2ggplot <- one_stat_df2ggplot[order(one_stat_df2ggplot$reset, one_stat_df2ggplot$time, one_stat_df2ggplot$condition),]
  
  p <- ggplot(one_stat_df2ggplot, aes(x=time, y=value, fill=condition, group=sample)) +
    geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8)) +
    labs(title=sprintf("%s", gsub("\\.", " ", one_stat)), x="Echantillon", y="Valeur", fill="Condition") +
    theme_bw() +
    plot_theme
  print(p)
  
  p <- ggpaired(one_stat_df2ggplot, x="reset", y="value", color="reset", line.color="gray", line.size=0.4) +
    stat_compare_means(paired=TRUE) +
    labs(title=sprintf("%s", gsub("\\.", " ", one_stat)), x="Reset", y="Valeur") +
    theme(legend.position="none")
  print(p)
  
}
dev.off()
