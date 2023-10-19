############################################################################
###This script analyzes the ShortBRED profiling of DOME sample resistomes###
############################################################################

###Load the necessary packages
library(reshape2)
library(vegan)
library(ggplot2)
library(labdsv)
library(Maaslin2)
library(tidyverse)
library(tibble)
library(ggpubr)
library(rstatix)
library(stats)
library(dplyr)
library(NatParksPalettes)
library(rstatix)
library(coin)
library(readr)
library(Hmisc)
library(ggbeeswarm)
library(pheatmap)

###Creating necessary functions
`%notin%` <- Negate(`%in%`)

season_names <- list(
  'SPRING 1' = 'Spring 2019',
  'SUMMER 1' = 'Summer 2019',
  'FALL 1' = 'Fall 2019',
  'WINTER 1' = 'Winter 2019',
  'SPRING 2' = 'Spring 2020'
)

season_labeller <- function(variable, value) {
  return(season_names[value])
}

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

position_nudgedodge <- function(x = 0, y = 0, width = 0.75) {
  ggproto(NULL, PositionNudgedodge,
          x = x,
          y = y,
          width = width
  )
}

PositionNudgedodge <- ggproto("PositionNudgedodge", PositionDodge,
                              x = 0,
                              y = 0,
                              width = 0.3,
                              
                              setup_params = function(self, data) {
                                l <- ggproto_parent(PositionDodge,self)$setup_params(data)
                                append(l, list(x = self$x, y = self$y))
                              },
                              
                              compute_layer = function(self, data, params, layout) {
                                d <- ggproto_parent(PositionNudge,self)$compute_layer(data,params,layout)
                                d <- ggproto_parent(PositionDodge,self)$compute_layer(d,params,layout)
                                d
                              }
)

###Read in the ShortBRED combined txt file
shortbred_combined <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_shortbred_reads/ShortBRED_combined.txt", header = TRUE)
shortbred_counts_wide <- dcast(shortbred_combined, Sample ~ Family, value.var = "Count")
#Two of the cow samples (19-C0205-W001fc01 & 19-C0416-F001fc01) have been clustering with human samples when looking at taxa (shotgun and 16S) and ARGs
#These two samples have likley been misplabeled. Removing these two from further analysis
shortbred_counts_wide <- shortbred_counts_wide[which(shortbred_counts_wide$Sample != "19-C0205-W001fc01" & shortbred_counts_wide$Sample != "19-C0416-F001fc01"),]
rownames(shortbred_counts_wide) <- shortbred_counts_wide$Sample
shortbred_counts_wide <- shortbred_counts_wide[,2:9294]
t_shortbred_counts_wide <- as.data.frame(t(shortbred_counts_wide))
#write.csv(t_shortbred_counts_wide, file = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_shortbred_reads/ShortBRED_combined_wide_transposed_filtered.csv")
#write.csv(shortbred_counts_wide, file = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_shortbred_reads/ShortBRED_combined_wide_filtered.csv")


###Read in the metadata
samples_metadata <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/220504_DOME_SamplesMasterList_v10.csv", header = TRUE)
samples_metadata <- samples_metadata[which(samples_metadata$Sample_Type == "Fecal"),1:16]
samples_metadata <- samples_metadata[which(samples_metadata$Sample_ID != "19-C0205-W001fc01" & samples_metadata$Sample_ID != "19-C0416-F001fc01"),]

###Read in the mapping file containing the abx class information for the ShortBRED hits
ARG_mapping <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/ShortBRED/shortBRED-identify/220529/Mapping/all_annot_AbxClass.txt", header = FALSE)
colnames(ARG_mapping)[1] <- "Gene"
colnames(ARG_mapping)[2] <- "AbxClass"
shortbred_AbxClass <- merge(ARG_mapping, t_shortbred_counts_wide, by.x = "Gene", by.y = "row.names")
shortbred_AbxClass <- shortbred_AbxClass[2:dim(shortbred_AbxClass)[2]]
shortbred_AbxClass <- shortbred_AbxClass %>% group_by(AbxClass) %>% summarise_all(sum)
shortbred_AbxClass <- as.data.frame(shortbred_AbxClass)
rownames(shortbred_AbxClass) <- shortbred_AbxClass$AbxClass
shortbred_AbxClass <- shortbred_AbxClass[2:dim(shortbred_AbxClass)[2]]
shortbred_AbxClass <- as.data.frame(t(shortbred_AbxClass))
shortbred_AbxClass_subjects <- merge(samples_metadata[,c(1,3)], shortbred_AbxClass, by.x = "Sample_ID", by.y = "row.names")
shortbred_AbxClass_subjects <- shortbred_AbxClass_subjects[2:dim(shortbred_AbxClass_subjects)[2]]
shortbred_AbxClass_subjects <- shortbred_AbxClass_subjects %>% group_by(SID) %>% summarise_all(mean)
shortbred_AbxClass_subjects <- as.data.frame(shortbred_AbxClass_subjects)
rownames(shortbred_AbxClass_subjects) <- shortbred_AbxClass_subjects$SID
shortbred_AbxClass_subjects <- shortbred_AbxClass_subjects[2:dim(shortbred_AbxClass_subjects)[2]]

###Calculating the significance in the difference of ARG burden (total RPKM) across subjects in each season
ARGburden <- as.data.frame(rowSums(shortbred_AbxClass))
ARGburden$Sample <- rownames(ARGburden)
colnames(ARGburden)[1] <- "TotalRPKM"
ARGburden <- merge(ARGburden, samples_metadata[,c(3,12,13,16)], by.x = "Sample", by.y = "Sample_ID")
ARGburden$SeasonFull <- paste(ARGburden$Season, ARGburden$Year)
ARGburden_MeansCompared <- compare_means(TotalRPKM ~ Group, data = ARGburden[which(ARGburden$SeasonFull != "SPRING 2"),], group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")
ARGburden_MeansCompared <- ARGburden_MeansCompared[which(ARGburden_MeansCompared$p.adj < 0.05),]

ARGburden$SeasonFull <- factor(ARGburden$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))
ARGburden$Group <- factor(ARGburden$Group, levels = c("C", "W", "D"))

p <- ggplot(data = ARGburden, aes(x = Group, y = TotalRPKM, color = Group)) + 
  geom_boxplot(aes(fill = Group), color = "black", outlier.shape = NA, width = 0.4, size = 1, alpha = 0.5) + 
  geom_flat_violin(aes(fill = Group), position = position_nudge(x = 0.3), alpha = 0.7, scale = "width", width = 0.45) +
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.5) +
  theme_linedraw() + 
  ylim(NA, 17000) +
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text = element_blank(), strip.text = element_text(face = "bold", size = 20),
        legend.position = "none", axis.title = element_blank()) +
  scale_color_manual(values = c("#873e23", "#18678d", "#626262")) + 
  scale_fill_manual(values = c("#873e23", "#18678d", "#626262"))
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/Shotgun_ARGburdenRPKM_v1.svg", p, bg = "transparent", width = 15, height = 10)


###Plotting abx class abundances
AbxClass_plot <- shortbred_AbxClass[,order(colSums(shortbred_AbxClass, na.rm = TRUE),decreasing=TRUE)]

top10_AbxClasses <- colnames(AbxClass_plot)[1:12]
top10_AbxClasses <- top10_AbxClasses[top10_AbxClasses != "OTHER"]
top10_AbxClasses <- top10_AbxClasses[top10_AbxClasses != "EFFLUX"]

AbxClass_plot <- data.frame(AbxClass_plot[,colnames(AbxClass_plot) %in% top10_AbxClasses],Others=rowSums(AbxClass_plot[,!colnames(AbxClass_plot) %in% top10_AbxClasses]))
AbxClass_plot <- melt(as.matrix(AbxClass_plot))
colnames(AbxClass_plot)[1] <- "Sample"
colnames(AbxClass_plot)[2] <- "AbxClass"
colnames(AbxClass_plot)[3] <- "Value"

AbxClass_plot <- merge(AbxClass_plot, samples_metadata[,c(3,12,13,16)], by.x = "Sample", by.y = "Sample_ID")
AbxClass_plot$SeasonFull <- paste(AbxClass_plot$Season, AbxClass_plot$Year)
AbxClass_plot <- AbxClass_plot[,c(2,3,5,7)]
AbxClass_plot <- AbxClass_plot %>% group_by(AbxClass, Group, SeasonFull) %>% summarise(Value = mean(Value))
AbxClass_plot$Group <- factor(AbxClass_plot$Group, levels = c("C", "W", "D"))
AbxClass_plot$SeasonFull <- factor(AbxClass_plot$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))

ggplot(data = AbxClass_plot, aes(x = Group, y = Value, fill = AbxClass)) + geom_bar(stat = "identity") +
  theme_test() + 
  scale_fill_manual(values = natparks.pals("DeathValley", 11)) + 
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) +
  ylab("RPKM") +
  scale_x_discrete(labels = c("Cow", "Farmer", "Office"))

ggplot(data = AbxClass_plot, aes(x = Group, y = Value, fill = AbxClass)) + geom_bar(stat = "identity") +
  theme_test() + 
  scale_fill_manual(values = c(natparks.pals("DeathValley", 11)[0:10], "#7d7d7d")) + 
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) +
  ylab("RPKM") +
  scale_x_discrete(labels = c("Cow", "Farmer", "Office"))

ggplot(data = AbxClass_plot, aes(x = Group, y = Value, fill = AbxClass)) + geom_bar(stat = "identity") +
  theme_test() + 
  scale_fill_manual(values = c(natparks.pals("DeathValley", 10), "#bcbcbc")) + 
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) +
  ylab("RPKM") +
  scale_x_discrete(labels = c("Cow", "Farmer", "Office"))
ggplot(data = AbxClass_plot, aes(x = Group, y = Value, fill = AbxClass)) + 
  geom_bar(stat = "identity", color = "white", size = 0.1) +
  theme_linedraw() + 
  scale_fill_manual(values = c(natparks.pals("DeathValley", 10), "#848484")) + 
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) +
  ylab("RPKM") +
  scale_x_discrete(labels = c("Cow", "Farmer", "Office")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), 
        panel.grid.minor.y = element_blank(), strip.text = element_text(face = "bold", size = 20))
ggplot(data = AbxClass_plot, aes(x = Group, y = Value, fill = AbxClass)) + 
  geom_bar(stat = "identity", color = "white", size = 0.1) +
  theme_linedraw() + 
  scale_fill_manual(values = c("#8C2B0E", "#B25422", "#D8813B", "#FEB359", "#9F7E59", "#233F6C", "#435F90", "#5B4C64", "#81565F", "#B47E83", "#848484")) + 
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) +
  ylab("RPKM") +
  scale_x_discrete(labels = c("Cow", "Farmer", "Office")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), 
        panel.grid.minor.y = element_blank(), strip.text = element_text(face = "bold", size = 20))
p <- ggplot(data = AbxClass_plot, aes(x = Group, y = Value, fill = AbxClass)) + 
  geom_bar(stat = "identity", color = "white", size = 0.1) +
  theme_linedraw() + 
  scale_fill_manual(values = c("#8C2B0E", "#B25422", "#D8813B", "#FEB359", "#9F7E59", "#233F6C", "#435F90", "#5B4C64", "#81565F", "#B47E83", "#848484")) + 
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) +
  labs(fill = "Antibiotic class") +
  scale_x_discrete(labels = c("Cow", "Farmer", "Office")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), 
        panel.grid.minor.y = element_blank(), strip.text = element_text(face = "bold", size = 20),
        axis.title = element_blank(), axis.text = element_blank(), legend.position = "none")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Shotgun_Resistome_RPKM_v4.svg", p, bg = "transparent", width = 15, height = 10)

###Get the ARG richness for all samples and plot
ARG_richness <- as.data.frame(specnumber(shortbred_counts_wide))
colnames(ARG_richness)[1] <- "ARG_richness"
ARG_richness$Sample <- rownames(ARG_richness)
rownames(ARG_richness) <- 1:710

ARG_richness <- merge(ARG_richness, samples_metadata, by.x = "Sample", by.y = "Sample_ID")
ARG_richness$SeasonFull <- paste(ARG_richness$Season, ARG_richness$Year)
ARG_richness$SeasonFull <- factor(ARG_richness$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))
ARG_richness$GroupSeasonFull <- paste(ARG_richness$Group, ARG_richness$Season, ARG_richness$Year)
ARG_richness$Group <- factor(ARG_richness$Group, levels = c("C", "W", "D"))

#richness_means_compared <- compare_means(ARG_richness ~ GroupSeasonFull, data = ARG_richness, method = "wilcox.test", p.adjust.method = "BH")
#compare_means(ARG_richness ~ Group, data = ARG_richness[which(ARG_richness$SeasonFull == "SPRING 1"),], method = "wilcox.test", p.adjust.method = "BH")
#compare_means(ARG_richness ~ Group, data = ARG_richness[which(ARG_richness$SeasonFull == "SUMMER 1"),], method = "wilcox.test", p.adjust.method = "BH")
#compare_means(ARG_richness ~ Group, data = ARG_richness[which(ARG_richness$SeasonFull == "FALL 1"),], method = "wilcox.test", p.adjust.method = "BH")
#compare_means(ARG_richness ~ Group, data = ARG_richness[which(ARG_richness$SeasonFull == "WINTER 1"),], method = "wilcox.test", p.adjust.method = "BH")
#compare_means(ARG_richness ~ Group, data = ARG_richness[which(ARG_richness$SeasonFull == "SPRING 2"),], method = "wilcox.test", p.adjust.method = "BH")
richness_means_compared <- compare_means(ARG_richness ~ Group, data = ARG_richness[which(ARG_richness$SeasonFull != "SPRING 2"),], group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")
richness_means_compared <- richness_means_compared[which(richness_means_compared$p.adj < 0.05),]
pval_positions <- ARG_richness %>%
  group_by(SeasonFull) %>%
  t_test(ARG_richness ~ Group)
pval_positions <- pval_positions %>%
  add_xy_position(x = "SeasonFull", dodge = 0.8)
pval_positions <- pval_positions[, c('SeasonFull', 'group1', 'group2', 'y.position', 'groups', 'x', 'xmin', 'xmax')]
richness_means_compared <- merge(richness_means_compared, pval_positions, by = c('SeasonFull', 'group1', 'group2'))

ggplot(data = ARG_richness, aes(x = SeasonFull, y = ARG_richness, color = Group)) + 
  geom_boxplot() + 
  theme_bw() + 
  ylim(NA, 900) + 
  stat_pvalue_manual(richness_means_compared, label = "p.adj") + 
  ylab("ARG richness") + 
  xlab("Season") + 
  labs(title = "ARG richness of subject types across seasons", subtitle = "Wilcoxon rank sum test. Benjamini-Hochberg correction.") +
  scale_color_discrete(name = "Subject", labels = c("Cows", "Office workers", "Farmers"))

p <- ggplot(data = ARG_richness, aes(x = Group, y = ARG_richness, color = Group)) + 
  geom_boxplot(aes(fill = Group), color = "black", outlier.shape = NA, width = 0.4, size = 1, alpha = 0.5) + 
  geom_flat_violin(aes(fill = Group), position = position_nudge(x = 0.3), alpha = 0.7, scale = "width", width = 0.45) +
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.5) +
  theme_linedraw() + 
  ylim(NA, 750) +
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text = element_blank(), strip.text = element_text(face = "bold", size = 20),
        legend.position = "none", axis.title = element_blank()) +
  scale_color_manual(values = c("#873e23", "#18678d", "#626262")) + 
  scale_fill_manual(values = c("#873e23", "#18678d", "#626262"))
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Shotgun_ARGrichness_v1.svg", p, bg = "transparent", width = 15, height = 10)

###Look at the richness of each group (cows, office workers, farmers) individually
cows_richness_means_compared <- compare_means(ARG_richness ~ SeasonFull, data = ARG_richness[which(ARG_richness$Group == "C"),], method = "wilcox.test", p.adjust.method = "BH")
cows_richness_means_compared <- cows_richness_means_compared[which(cows_richness_means_compared$p.adj < 0.05),]

ggplot(data = ARG_richness[which(ARG_richness$Group == "C"),], aes(x = SeasonFull, y = ARG_richness)) + 
  geom_boxplot(color = c("#F8766D")) + 
  theme_bw() + 
  ylab("ARG richness") +
  xlab("Season") +
  stat_pvalue_manual(cows_richness_means_compared, y.position = 700, step.increase = 0.1, label = "p.adj") +
  labs(title = "ARG richness of cows across seasons", subtitle = "Wilcoxon rank sum test. Benjamini-Hochberg correction.") +
  theme(legend.position = "none")


office_richness_means_compared <- compare_means(ARG_richness ~ SeasonFull, data = ARG_richness[which(ARG_richness$Group == "D"),], method = "wilcox.test", p.adjust.method = "BH") #nothing was significant after adjustment

ggplot(data = ARG_richness[which(ARG_richness$Group == "D"),], aes(x = SeasonFull, y = ARG_richness)) + 
  geom_boxplot(color = c("#00BA38")) + 
  theme_bw() + 
  ylab("ARG richness") +
  xlab("Season") +
  labs(title = "ARG richness of office workers across seasons", subtitle = "Wilcoxon rank sum test. Benjamini-Hochberg correction.") +
  theme(legend.position = "none")


farmers_richness_means_compared <- compare_means(ARG_richness ~ SeasonFull, data = ARG_richness[which(ARG_richness$Group == "W"),], method = "wilcox.test", p.adjust.method = "BH")
farmers_richness_means_compared <- farmers_richness_means_compared[which(farmers_richness_means_compared$p.adj < 0.05),]

ggplot(data = ARG_richness[which(ARG_richness$Group == "W"),], aes(x = SeasonFull, y = ARG_richness)) + 
  geom_boxplot(fill = "#18678d", color = "black", outlier.shape = NA, width = 0.4, size = 1, alpha = 0.5) +
  geom_flat_violin(fill = "#18678d", position = position_nudge(x = 0.3), alpha = 0.7, scale = "width", width = 0.45) +
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.5, color = "#18678d") +
  theme_linedraw() + 
  ylab("ARG richness") +
  xlab("Season") +
  stat_pvalue_manual(farmers_richness_means_compared, y.position = 700, step.increase = 0.1, label = "p.adj") +
  labs(title = "ARG richness of farmers across seasons", subtitle = "Wilcoxon rank sum test. Benjamini-Hochberg correction.") +
  theme(legend.position = "none")
p <- ggplot(data = ARG_richness[which(ARG_richness$Group == "W"),], aes(x = SeasonFull, y = ARG_richness)) + 
  geom_boxplot(fill = "#18678d", color = "black", outlier.shape = NA, width = 0.5, size = 1, alpha = 0.5, position = position_nudge(x = -.2)) +
  geom_flat_violin(fill = "#18678d", position = position_nudge(x = 0.11), alpha = 0.7, scale = "width", width = 0.5) +
  geom_point(position = position_dodge2nudge(width = 0.4, x = -.2), alpha = 0.5, color = "#18678d", size = 4) +
  theme_linedraw() + 
  scale_y_continuous(limits = c(NA, 760), breaks = seq(200, 800, 200)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text = element_blank(), strip.text = element_text(face = "bold", size = 20),
        legend.position = "none", axis.title = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Shotgun_ARGrichness_Farmers_v1.svg", p, bg = "transparent", width = 10, height = 10)

###Look at the Bray-Curtis dissimilarity of sample resistomes
ARG_BrayCurtis <- vegdist(shortbred_counts_wide, method = "bray")
pcoa_BC <- pco(ARG_BrayCurtis, k = 4)
pcoa_BC_scree <- as.data.frame(pcoa_BC$eig*100/sum(pcoa_BC$eig))
pcoa_BC <- as.data.frame(pcoa_BC$points)
pcoa_BC <- merge(pcoa_BC, samples_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_BC)[1] <- "Sample"
pcoa_BC$SeasonFull <- paste(pcoa_BC$Season, pcoa_BC$Year)
pcoa_BC$SeasonFull <- factor(pcoa_BC$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))
pcoa_BC$Group <- factor(pcoa_BC$Group, levels = c("C", "W", "D"))

ggplot(data = pcoa_BC, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull)) +
  theme_test()
p <- ggplot(data = pcoa_BC, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull), alpha = 0.75, size = 4) +
  theme_test() + stat_ellipse(aes(color = Group), size = 2.5) +
  scale_color_manual(values = c("#873e23", "#18678d", "#626262"), labels = c("Cows", "Farmers", "Controls")) + xlab("PCoA 1 (47.6%)") + ylab("PCoA 2 (5.0%)") +
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"), axis.text = element_text(), axis.title = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 20), legend.position = "none") +
  scale_shape_discrete(name = "Season", labels = c("Spring 2019", "Summer 2019", "Fall 2019", "Winter 2019", "Spring 2020"))
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Shotgun_Resistome_PCOA_v1.svg", p, bg = "transparent", width = 15, height = 10)

p <- ggplot(data = pcoa_BC[which(pcoa_BC$SeasonFull != "SPRING 2"),], aes(x = V1, y = V2)) + geom_point(aes(color = Group), alpha = 0.75, size = 6) +
  theme_test() + stat_ellipse(aes(color = Group), size = 2.5) +
  scale_color_manual(values = c("#873e23", "#18678d", "#626262"), labels = c("Cows", "Farmers", "Controls")) + xlab("PCoA 1 (47.6%)") + ylab("PCoA 2 (5.0%)") +
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"), axis.text = element_blank(), axis.title = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 20), legend.position = "none") +
  scale_shape_manual(values = c(15, 15, 15, 15), name = "Season", labels = c("Spring 2019", "Summer 2019", "Fall 2019", "Winter 2019"))
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Shotgun_Resistome_PCOA_v3.svg", p, bg = "transparent", width = 11, height = 9.8)

###Looking at Bray-Curtis dissimilarity within each of the seasons separately
metadata_SPRING1 <- samples_metadata[which(samples_metadata$Season == "SPRING" & samples_metadata$Year == "1"),]
metadata_SUMMER1 <- samples_metadata[which(samples_metadata$Season == "SUMMER" & samples_metadata$Year == "1"),]
metadata_FALL1 <- samples_metadata[which(samples_metadata$Season == "FALL" & samples_metadata$Year == "1"),]
metadata_WINTER1 <- samples_metadata[which(samples_metadata$Season == "WINTER" & samples_metadata$Year == "1"),]
metadata_SPRING2 <- samples_metadata[which(samples_metadata$Season == "SPRING" & samples_metadata$Year == "2"),]

ARG_BC_SPRING1 <- vegdist(shortbred_counts_wide[which(row.names(shortbred_counts_wide) %in% metadata_SPRING1$Sample_ID),], index = "bray")
pcoa_BC_SPRING1 <- pco(ARG_BC_SPRING1, k = 4)
pcoa_BC_SPRING1 <- as.data.frame(pcoa_BC_SPRING1$points)
pcoa_BC_SPRING1 <- merge(pcoa_BC_SPRING1, metadata_SPRING1, by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_BC_SPRING1)[1] <- "Sample"
pcoa_BC_SPRING1$SeasonFull <- paste(pcoa_BC_SPRING1$Season, pcoa_BC_SPRING1$Year)

ggplot(data = pcoa_BC_SPRING1, aes(x = V1, y = V2)) + geom_point(aes(color = Group))


ARG_BC_SUMMER1 <- vegdist(shortbred_counts_wide[which(row.names(shortbred_counts_wide) %in% metadata_SUMMER1$Sample_ID),], index = "bray")
pcoa_BC_SUMMER1 <- pco(ARG_BC_SUMMER1, k = 4)
pcoa_BC_SUMMER1 <- as.data.frame(pcoa_BC_SUMMER1$points)
pcoa_BC_SUMMER1 <- merge(pcoa_BC_SUMMER1, metadata_SUMMER1, by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_BC_SUMMER1)[1] <- "Sample"

ggplot(data = pcoa_BC_SUMMER1, aes(x = V1, y = V2)) + geom_point(aes(color = Group))


ARG_BC_FALL1 <- vegdist(shortbred_counts_wide[which(row.names(shortbred_counts_wide) %in% metadata_FALL1$Sample_ID),], index = "bray")
pcoa_BC_FALL1 <- pco(ARG_BC_FALL1, k = 4)
pcoa_BC_FALL1 <- as.data.frame(pcoa_BC_FALL1$points)
pcoa_BC_FALL1 <- merge(pcoa_BC_FALL1, metadata_FALL1, by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_BC_FALL1)[1] <- "Sample"

ggplot(data = pcoa_BC_FALL1, aes(x = V1, y = V2)) + geom_point(aes(color = Group))


ARG_BC_WINTER1 <- vegdist(shortbred_counts_wide[which(row.names(shortbred_counts_wide) %in% metadata_WINTER1$Sample_ID),], index = "bray")
pcoa_BC_WINTER1 <- pco(ARG_BC_WINTER1, k = 4)
pcoa_BC_WINTER1 <- as.data.frame(pcoa_BC_WINTER1$points)
pcoa_BC_WINTER1 <- merge(pcoa_BC_WINTER1, metadata_WINTER1, by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_BC_WINTER1)[1] <- "Sample"

ggplot(data = pcoa_BC_WINTER1, aes(x = V1, y = V2)) + geom_point(aes(color = Group))


ARG_BC_SPRING2 <- vegdist(shortbred_counts_wide[which(row.names(shortbred_counts_wide) %in% metadata_SPRING2$Sample_ID),], index = "bray")
pcoa_BC_SPRING2 <- pco(ARG_BC_SPRING2, k = 4)
pcoa_BC_SPRING2 <- as.data.frame(pcoa_BC_SPRING2$points)
pcoa_BC_SPRING2 <- merge(pcoa_BC_SPRING2, metadata_SPRING2, by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_BC_SPRING2)[1] <- "Sample"

ggplot(data = pcoa_BC_SPRING2, aes(x = V1, y = V2)) + geom_point(aes(color = Group))

###Getting the pairwise Bray-Curtis distances for sample resistomes
ARG_BrayCurtis_pairwise <- as.matrix(ARG_BrayCurtis)
ARG_BrayCurtis_pairwise <- melt(ARG_BrayCurtis_pairwise)
ARG_BrayCurtis_pairwise <- ARG_BrayCurtis_pairwise[!duplicated(apply(ARG_BrayCurtis_pairwise, 1, function(x) paste(sort(x), collapse = ""))),]
ARG_BrayCurtis_pairwise <- ARG_BrayCurtis_pairwise[which(ARG_BrayCurtis_pairwise$X1 != ARG_BrayCurtis_pairwise$X2),]
ARG_BrayCurtis_pairwise <- merge(ARG_BrayCurtis_pairwise, samples_metadata[,c(1,3,12,13,15,16)], by.x = "X1", by.y = "Sample_ID")
colnames(ARG_BrayCurtis_pairwise)[4] <- "SID1"
colnames(ARG_BrayCurtis_pairwise)[5] <- "Season 1"
colnames(ARG_BrayCurtis_pairwise)[6] <- "Group 1"
colnames(ARG_BrayCurtis_pairwise)[7] <- "Site 1"
colnames(ARG_BrayCurtis_pairwise)[8] <- "Year 1"
ARG_BrayCurtis_pairwise <- merge(ARG_BrayCurtis_pairwise, samples_metadata[,c(1,3,12,13,15,16)], by.x = "X2", by.y = "Sample_ID")
colnames(ARG_BrayCurtis_pairwise)[9] <- "SID2"
colnames(ARG_BrayCurtis_pairwise)[10] <- "Season 2"
colnames(ARG_BrayCurtis_pairwise)[11] <- "Group 2"
colnames(ARG_BrayCurtis_pairwise)[12] <- "Site 2"
colnames(ARG_BrayCurtis_pairwise)[13] <- "Year 2"

ARG_BrayCurtis_pairwise_filtered <- ARG_BrayCurtis_pairwise[which(ARG_BrayCurtis_pairwise$`Season 1` == ARG_BrayCurtis_pairwise$`Season 2` & ARG_BrayCurtis_pairwise$`Year 1` == ARG_BrayCurtis_pairwise$`Year 2`),]
ARG_BrayCurtis_pairwise_filtered <- ARG_BrayCurtis_pairwise_filtered[which(ARG_BrayCurtis_pairwise_filtered$`Group 1` != ARG_BrayCurtis_pairwise_filtered$`Group 2`),]
ARG_BrayCurtis_pairwise_filtered$Groups <- paste(ARG_BrayCurtis_pairwise_filtered$`Group 1`, ARG_BrayCurtis_pairwise_filtered$`Group 2`)
ARG_BrayCurtis_pairwise_filtered <- ARG_BrayCurtis_pairwise_filtered[which(ARG_BrayCurtis_pairwise_filtered$Groups != "D W" & ARG_BrayCurtis_pairwise_filtered$Groups != "W D"),]
ARG_BrayCurtis_pairwise_filtered$Sites[ARG_BrayCurtis_pairwise_filtered$`Site 1` == ARG_BrayCurtis_pairwise_filtered$`Site 2`] <- "Same site"
ARG_BrayCurtis_pairwise_filtered$Sites[ARG_BrayCurtis_pairwise_filtered$`Site 1` != ARG_BrayCurtis_pairwise_filtered$`Site 2`] <- "Different sites"
ARG_BrayCurtis_pairwise_filtered$SeasonFull <- paste(ARG_BrayCurtis_pairwise_filtered$`Season 1`, ARG_BrayCurtis_pairwise_filtered$`Year 1`)

ARG_BrayCurtis_pairwise_filtered <- ARG_BrayCurtis_pairwise_filtered %>%
  group_by(SID1, Sites, SeasonFull, `Group 1`) %>%
  summarise(
    value = mean(value)
  )

ARG_BrayCurtis_pairwise_mean <- ARG_BrayCurtis_pairwise_filtered[c(1,2,4,5)] %>%
  group_by(SID1, Sites, `Group 1`) %>%
  summarise(
    value = mean(value)
  )

ARG_BrayCurtis_pairwise_filtered$GroupSites <- paste(ARG_BrayCurtis_pairwise_filtered$`Group 1`, ARG_BrayCurtis_pairwise_filtered$Sites)
ARG_BrayCurtis_pairwise_filtered$SeasonFull <- factor(ARG_BrayCurtis_pairwise_filtered$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))

ARG_BrayCurtis_pairwise_mean$GroupSites <- paste(ARG_BrayCurtis_pairwise_mean$`Group 1`, ARG_BrayCurtis_pairwise_mean$Sites)
ARG_BrayCurtis_pairwise_mean$`Group 1` <- factor(ARG_BrayCurtis_pairwise_mean$`Group 1`, levels = c("W", "D"))

BC_pairwise_compared <- compare_means(value ~ GroupSites, data = ARG_BrayCurtis_pairwise_filtered, group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")
BC_pairwise_compared <- BC_pairwise_compared[which(BC_pairwise_compared$p.adj < 0.05),]

compare_means(value ~ GroupSites, data = ARG_BrayCurtis_pairwise_mean, method = "wilcox.test", p.adjust.method = "BH")

ggplot(data = ARG_BrayCurtis_pairwise_filtered, aes(x = Sites, y = value, color = `Group 1`)) + geom_boxplot() + 
  facet_wrap(~ARG_BrayCurtis_pairwise_filtered$SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_test() + ylab("Bray-Curtis distance") +
  ylim(0,1.2) +
  labs(title = "Average (within subject) ARG Bray-Curtis distances to cow samples", subtitle = "Wilcoxon rank sum test. Benjamini-Hochberg correction.") +
  scale_color_manual(values=c("#00BA38", "#619CFF"), name = "", labels = c("Office worker", "Farmer"))
  
ggplot(data = ARG_BrayCurtis_pairwise_mean, aes(x = Sites, y = value, color = `Group 1`)) + geom_boxplot() +
  theme_test() + ylab("Bray-Curtis distance") +
  labs(title = "Average (within subject) ARG Bray-Curtis distances to cow samples", subtitle = "Wilcoxon rank sum test. Benjamini-Hochberg correction.") +
  scale_color_manual(values=c("#00BA38", "#619CFF"), name = "")
ggplot(data = ARG_BrayCurtis_pairwise_mean, aes(x = Sites, y = value)) + geom_boxplot(aes(fill = `Group 1`), color = "black", outlier.shape = NA, width = 0.5, size = 1, alpha = 0.5, position = position_dodge(width = 0.9)) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.9, jitter.height = 0, jitter.width = 0.45), alpha = 0.5, aes(group = `Group 1`, color = `Group 1`)) +
  geom_flat_violin(aes(fill = `Group 1`), position = position_nudgedodge(x = 0.15, width = 0.9), alpha = 0.7, scale = "width", width = 1.2) +
  theme_linedraw() +
  ylim(0.6, 1) +
  scale_fill_manual(values=c("#18678d", "#626262"), name = "", labels = c("Farmers", "Controls")) +
  scale_color_manual(values=c("#18678d", "#626262"), name = "", labels = c("Farmers", "Controls")) +
  scale_x_discrete(labels = c("Different\nsites", "Same\nsite")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text = element_text(),
        legend.position = "right", axis.title.y = element_text(), axis.title.x = element_text())
p <- ggplot(data = ARG_BrayCurtis_pairwise_mean, aes(x = Sites, y = value)) + geom_boxplot(aes(fill = `Group 1`), color = "black", outlier.shape = NA, width = 0.5, size = 1, alpha = 0.5, position = position_dodge(width = 0.9)) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.9, jitter.height = 0, jitter.width = 0.45), alpha = 0.5, aes(group = `Group 1`, color = `Group 1`)) +
  geom_flat_violin(aes(fill = `Group 1`), position = position_nudgedodge(x = 0.15, width = 0.9), alpha = 0.7, scale = "width", width = 1.2) +
  theme_linedraw() +
  ylim(0.6, 1) +
  scale_fill_manual(values=c("#18678d", "#626262"), name = "",) +
  scale_color_manual(values=c("#18678d", "#626262"), name = "",) +
  scale_x_discrete(labels = c("Different\nsites", "Same\nsite")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text = element_blank(),
        legend.position = "none", axis.title = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/Shortbred_PairwiseBC_v1.svg", p, bg = "transparent", width = 6, height = 10)

######Look at the Bray-Curtis dissimilarity of resistomes of farmers and office workers
human_metadata <- samples_metadata[which(samples_metadata$Group != "C"),]
human_metadata <- human_metadata[order(human_metadata$Sample_ID),]
shortbred_counts_wide_humans <- shortbred_counts_wide[which(row.names(shortbred_counts_wide) %in% human_metadata$Sample_ID),]
shortbred_counts_wide_humans <- shortbred_counts_wide_humans[order(row.names(shortbred_counts_wide_humans)),]
human_ARG_BrayCurtis <- vegdist(shortbred_counts_wide_humans, index = "bray")
human_pcoa_BC <- pco(human_ARG_BrayCurtis, k = 4)
human_pcoa_BC <- as.data.frame(human_pcoa_BC$points)
human_pcoa_BC <- merge(human_pcoa_BC, human_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(human_pcoa_BC)[1] <- "Sample"
human_pcoa_BC$SeasonFull <- paste(human_pcoa_BC$Season, human_pcoa_BC$Year)

ggplot(data = human_pcoa_BC, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull)) + stat_ellipse(aes(color = Group)) + 
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test()

#With adonis, make sure that the row names of the distance matrix match the Sample_IDs in the metadata
adonis2(shortbred_counts_wide_humans ~ human_metadata$Group,
       permutations = 9999, method = "bray")
#p < 0.001

###Looking at the Bray-Curtis dissimilarity of resistomes of farmers and office workers within each season
human_metadata_SPRING1 <- samples_metadata[which(samples_metadata$Season == "SPRING" & samples_metadata$Year == "1" & samples_metadata$Group != "C"),]
human_metadata_SUMMER1 <- samples_metadata[which(samples_metadata$Season == "SUMMER" & samples_metadata$Year == "1" & samples_metadata$Group != "C"),]
human_metadata_FALL1 <- samples_metadata[which(samples_metadata$Season == "FALL" & samples_metadata$Year == "1" & samples_metadata$Group != "C"),]
human_metadata_WINTER1 <- samples_metadata[which(samples_metadata$Season == "WINTER" & samples_metadata$Year == "1" & samples_metadata$Group != "C"),]
human_metadata_SPRING2 <- samples_metadata[which(samples_metadata$Season == "SPRING" & samples_metadata$Year == "2" & samples_metadata$Group != "C"),]

shortbred_counts_wide_humans_SPRING1 <- shortbred_counts_wide[which(row.names(shortbred_counts_wide) %in% human_metadata_SPRING1$Sample_ID),]
ARG_BC_SPRING1_human <- vegdist(shortbred_counts_wide_humans_SPRING1, index = "bray")
pcoa_BC_SPRING1_human <- pco(ARG_BC_SPRING1_human, k = 4)
pcoa_BC_SPRING1_human_scree <- as.data.frame(100*pcoa_BC_SPRING1_human$eig/sum(pcoa_BC_SPRING1_human$eig))
pcoa_BC_SPRING1_human <- as.data.frame(pcoa_BC_SPRING1_human$points)
pcoa_BC_SPRING1_human <- merge(pcoa_BC_SPRING1_human, human_metadata_SPRING1, by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_BC_SPRING1_human)[1] <- "Sample"
pcoa_BC_SPRING1_human$Group <- factor(pcoa_BC_SPRING1_human$Group, levels = c("W", "D"))

human_metadata_SPRING1 <- human_metadata_SPRING1[order(human_metadata_SPRING1$Sample_ID),]

adonis2(shortbred_counts_wide_humans_SPRING1 ~ human_metadata_SPRING1$Group,
       permutations = 9999, method = "bray")

ggplot(data = pcoa_BC_SPRING1_human, aes(x = V1, y = V2)) + geom_point(aes(color = Group)) +
  theme_test() +
  labs(title = "PCOA (Bray-Curtis) of human samples from SPRING 1", subtitle = "PERMANOVA. Farmer vs. Office worker. BH correction. P = 0.063.")
p <- ggplot(data = pcoa_BC_SPRING1_human, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull), alpha = 1, size = 5, shape = 16) +
  theme_test() + stat_ellipse(aes(color = Group), size = 2.5, alpha = 0.75) +
  scale_color_manual(values = c("#18678d", "#626262")) + xlab("PCoA 1 (19.5%)") + ylab("PCoA 2 (14.1%)") +
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"), axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 20),
        legend.position = "none")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/Resistome_Human_PCOA_Spring2019_v1.svg", p, bg = "transparent", width = 7, height = 7)


shortbred_counts_wide_humans_SUMMER1 <- shortbred_counts_wide[which(row.names(shortbred_counts_wide) %in% human_metadata_SUMMER1$Sample_ID),]
ARG_BC_SUMMER1_human <- vegdist(shortbred_counts_wide_humans_SUMMER1, index = "bray")
pcoa_BC_SUMMER1_human <- pco(ARG_BC_SUMMER1_human, k = 4)
pcoa_BC_SUMMER1_human_scree <- as.data.frame(100*pcoa_BC_SUMMER1_human$eig/sum(pcoa_BC_SUMMER1_human$eig))
pcoa_BC_SUMMER1_human <- as.data.frame(pcoa_BC_SUMMER1_human$points)
pcoa_BC_SUMMER1_human <- merge(pcoa_BC_SUMMER1_human, human_metadata_SUMMER1, by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_BC_SUMMER1_human)[1] <- "Sample"
pcoa_BC_SUMMER1_human$Group <- factor(pcoa_BC_SUMMER1_human$Group, levels = c("W", "D"))

human_metadata_SUMMER1 <- human_metadata_SUMMER1[order(human_metadata_SUMMER1$Sample_ID),]

adonis2(shortbred_counts_wide_humans_SUMMER1 ~ human_metadata_SUMMER1$Group,
        permutations = 9999, method = "bray")

ggplot(data = pcoa_BC_SUMMER1_human, aes(x = V1, y = V2)) + geom_point(aes(color = Group)) + 
  theme_test() +
  labs(title = "PCOA (Bray-Curtis) of human samples from SUMMER 1", subtitle = "PERMANOVA. Farmer vs. Office worker. BH correction. P = 0.253.")
p <- ggplot(data = pcoa_BC_SUMMER1_human, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull), alpha = 1, size = 5, shape = 17) +
  theme_test() + stat_ellipse(aes(color = Group), size = 2.5, alpha = 0.75) +
  scale_color_manual(values = c("#18678d", "#626262")) + xlab("PCoA 1 (15.6%)") + ylab("PCoA 2 (9.6%)") +
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"), axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 20),
        legend.position = "none")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/Resistome_Human_PCOA_Summer2019_v1.svg", p, bg = "transparent", width = 7, height = 7)


shortbred_counts_wide_humans_FALL1 <- shortbred_counts_wide[which(row.names(shortbred_counts_wide) %in% human_metadata_FALL1$Sample_ID),]
ARG_BC_FALL1_human <- vegdist(shortbred_counts_wide_humans_FALL1, index = "bray")
pcoa_BC_FALL1_human <- pco(ARG_BC_FALL1_human, k = 4)
pcoa_BC_FALL1_human_scree <- as.data.frame(100*pcoa_BC_FALL1_human$eig/sum(pcoa_BC_FALL1_human$eig))
pcoa_BC_FALL1_human <- as.data.frame(pcoa_BC_FALL1_human$points)
pcoa_BC_FALL1_human <- merge(pcoa_BC_FALL1_human, human_metadata_FALL1, by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_BC_FALL1_human)[1] <- "Sample"
pcoa_BC_FALL1_human$Group <- factor(pcoa_BC_FALL1_human$Group, levels = c("W", "D"))

human_metadata_FALL1 <- human_metadata_FALL1[order(human_metadata_FALL1$Sample_ID),]

adonis2(shortbred_counts_wide_humans_FALL1 ~ human_metadata_FALL1$Group,
        permutations = 9999, method = "bray")

ggplot(data = pcoa_BC_FALL1_human, aes(x = V1, y = V2)) + geom_point(aes(color = Group)) + 
  theme_test() +
  labs(title = "PCOA (Bray-Curtis) of human samples from FALL 1", subtitle = "PERMANOVA. Farmer vs. Office worker. BH correction. P = 0.0055.", x = "PCOA 1", y = "PCOA 2") +
  stat_ellipse(aes(color = Group)) +
  scale_color_manual(values=c("#00BA38", "#619CFF"), name = "", labels = c("Office worker", "Farmer")) + 
  theme(axis.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 12, face = "bold"), 
        legend.background = element_rect(fill = "lightgray"), legend.key = element_rect(fill = "lightgray"),
        legend.position = c(0.91, 0.93))
p <- ggplot(data = pcoa_BC_FALL1_human, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull), alpha = 1, size = 5, shape = 15) +
  theme_test() + stat_ellipse(aes(color = Group), size = 2.5, alpha = 0.75) +
  scale_color_manual(values = c("#18678d", "#626262")) + xlab("PCoA 1 (12.5%)") + ylab("PCoA 2 (9.0%)") +
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"), axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 20),
        legend.position = "none")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/Resistome_Human_PCOA_Fall2019_v1.svg", p, bg = "transparent", width = 7, height = 7)


shortbred_counts_wide_humans_WINTER1 <- shortbred_counts_wide[which(row.names(shortbred_counts_wide) %in% human_metadata_WINTER1$Sample_ID),]
ARG_BC_WINTER1_human <- vegdist(shortbred_counts_wide_humans_WINTER1, index = "bray")
pcoa_BC_WINTER1_human <- pco(ARG_BC_WINTER1_human, k = 4)
pcoa_BC_WINTER1_human_scree <- as.data.frame(100*pcoa_BC_WINTER1_human$eig/sum(pcoa_BC_WINTER1_human$eig))
pcoa_BC_WINTER1_human <- as.data.frame(pcoa_BC_WINTER1_human$points)
pcoa_BC_WINTER1_human <- merge(pcoa_BC_WINTER1_human, human_metadata_WINTER1, by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_BC_WINTER1_human)[1] <- "Sample"
pcoa_BC_WINTER1_human$Group <- factor(pcoa_BC_WINTER1_human$Group, levels = c("W", "D"))

human_metadata_WINTER1 <- human_metadata_WINTER1[order(human_metadata_WINTER1$Sample_ID),]

adonis2(shortbred_counts_wide_humans_WINTER1 ~ human_metadata_WINTER1$Group,
        permutations = 9999, method = "bray")

ggplot(data = pcoa_BC_WINTER1_human, aes(x = V1, y = V2)) + geom_point(aes(color = Group)) + 
  theme_test() +
  labs(title = "PCOA (Bray-Curtis) of human samples from WINTER 1", subtitle = "PERMANOVA. Farmer vs. Office worker. BH correction. P = 0.063.")
p <- ggplot(data = pcoa_BC_WINTER1_human, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull), alpha = 1, size = 5, shape = 3) +
  theme_test() + stat_ellipse(aes(color = Group), size = 2.5, alpha = 0.75) +
  scale_color_manual(values = c("#18678d", "#626262")) + xlab("PCoA 1 (11.5%)") + ylab("PCoA 2 (9.2%)") +
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"), axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 20),
        legend.position = "none")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/Resistome_Human_PCOA_Winter2019_v1.svg", p, bg = "transparent", width = 7, height = 7)


shortbred_counts_wide_humans_SPRING2 <- shortbred_counts_wide[which(row.names(shortbred_counts_wide) %in% human_metadata_SPRING2$Sample_ID),]
ARG_BC_SPRING2_human <- vegdist(shortbred_counts_wide_humans_SPRING2, index = "bray")
pcoa_BC_SPRING2_human <- pco(ARG_BC_SPRING2_human, k = 4)
pcoa_BC_SPRING2_human_scree <- as.data.frame(100*pcoa_BC_SPRING2_human$eig/sum(pcoa_BC_SPRING2_human$eig))
pcoa_BC_SPRING2_human <- as.data.frame(pcoa_BC_SPRING2_human$points)
pcoa_BC_SPRING2_human <- merge(pcoa_BC_SPRING2_human, human_metadata_SPRING2, by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_BC_SPRING2_human)[1] <- "Sample"
pcoa_BC_SPRING2_human$Group <- factor(pcoa_BC_SPRING2_human$Group, levels = c("W", "D"))

human_metadata_SPRING2 <- human_metadata_SPRING2[order(human_metadata_SPRING2$Sample_ID),]

adonis2(shortbred_counts_wide_humans_SPRING2 ~ human_metadata_SPRING2$Group,
        permutations = 9999, method = "bray")

ggplot(data = pcoa_BC_SPRING2_human, aes(x = V1, y = V2)) + geom_point(aes(color = Group)) + 
  theme_test() +
  labs(title = "PCOA (Bray-Curtis) of human samples from SPRING 2", subtitle = "PERMANOVA. Farmer vs. Office worker. BH correction. P = 0.4176.")
p <- ggplot(data = pcoa_BC_SPRING2_human, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull), alpha = 1, size = 5, shape = 7) +
  theme_test() + stat_ellipse(aes(color = Group), size = 2.5, alpha = 0.75) +
  scale_color_manual(values = c("#18678d", "#626262")) + xlab("PCoA 1 (14.7%)") + ylab("PCoA 2 (11.7%)") +
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"), axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 20),
        legend.position = "none")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/Resistome_Human_PCOA_Spring2020_v1.svg", p, bg = "transparent", width = 7, height = 7)

p.adjust(c(0.038, 0.2022, 0.0011, 0.0378, 0.4176), method = "BH")
p.adjust(c(0.038, 0.2022, 0.0011, 0.0378), method = "BH")

###Looking at Jaccard dissimilarity of the samples
ARG_Jaccard <- vegdist(shortbred_counts_wide, method = "jaccard")
pcoa_Jacc <- pco(ARG_Jaccard, k = 4)
pcoa_Jacc <- as.data.frame(pcoa_Jacc$points)
pcoa_Jacc <- merge(pcoa_Jacc, samples_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_Jacc)[1] <- "Sample"
pcoa_Jacc$SeasonFull <- paste(pcoa_Jacc$Season, pcoa_Jacc$Year)

ggplot(data = pcoa_Jacc, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull)) +
  theme_test()

###Getting the pairwise Bray-Curtis distances for samples
ARG_Jaccard_pairwise <- as.matrix(ARG_Jaccard)
ARG_Jaccard_pairwise <- melt(ARG_Jaccard_pairwise)
ARG_Jaccard_pairwise <- ARG_Jaccard_pairwise
ARG_Jaccard_pairwise <- ARG_Jaccard_pairwise[!duplicated(apply(ARG_Jaccard_pairwise, 1, function(x) paste(sort(x), collapse = ""))),]
ARG_Jaccard_pairwise <- ARG_Jaccard_pairwise[which(ARG_Jaccard_pairwise$X1 != ARG_Jaccard_pairwise$X2),]
ARG_Jaccard_pairwise <- merge(ARG_Jaccard_pairwise, samples_metadata[,c(1,3,12,13,15,16)], by.x = "X1", by.y = "Sample_ID")
colnames(ARG_Jaccard_pairwise)[4] <- "SID1"
colnames(ARG_Jaccard_pairwise)[5] <- "Season 1"
colnames(ARG_Jaccard_pairwise)[6] <- "Group 1"
colnames(ARG_Jaccard_pairwise)[7] <- "Site 1"
colnames(ARG_Jaccard_pairwise)[8] <- "Year 1"
ARG_Jaccard_pairwise <- merge(ARG_Jaccard_pairwise, samples_metadata[,c(1,3,12,13,15,16)], by.x = "X2", by.y = "Sample_ID")
colnames(ARG_Jaccard_pairwise)[9] <- "SID2"
colnames(ARG_Jaccard_pairwise)[10] <- "Season 2"
colnames(ARG_Jaccard_pairwise)[11] <- "Group 2"
colnames(ARG_Jaccard_pairwise)[12] <- "Site 2"
colnames(ARG_Jaccard_pairwise)[13] <- "Year 2"

ARG_Jaccard_pairwise_filtered <- ARG_Jaccard_pairwise[which(ARG_Jaccard_pairwise$`Season 1` == ARG_Jaccard_pairwise$`Season 2` & ARG_Jaccard_pairwise$`Year 1` == ARG_Jaccard_pairwise$`Year 2`),]
ARG_Jaccard_pairwise_filtered <- ARG_Jaccard_pairwise_filtered[which(ARG_Jaccard_pairwise_filtered$`Group 1` != ARG_Jaccard_pairwise_filtered$`Group 2`),]
ARG_Jaccard_pairwise_filtered$Groups <- paste(ARG_Jaccard_pairwise_filtered$`Group 1`, ARG_Jaccard_pairwise_filtered$`Group 2`)
ARG_Jaccard_pairwise_filtered <- ARG_Jaccard_pairwise_filtered[which(ARG_Jaccard_pairwise_filtered$Groups != "D W" & ARG_Jaccard_pairwise_filtered$Groups != "W D"),]
ARG_Jaccard_pairwise_filtered$Sites[ARG_Jaccard_pairwise_filtered$`Site 1` == ARG_Jaccard_pairwise_filtered$`Site 2`] <- "Same site"
ARG_Jaccard_pairwise_filtered$Sites[ARG_Jaccard_pairwise_filtered$`Site 1` != ARG_Jaccard_pairwise_filtered$`Site 2`] <- "Different sites"
ARG_Jaccard_pairwise_filtered$SeasonFull <- paste(ARG_Jaccard_pairwise_filtered$`Season 1`, ARG_Jaccard_pairwise_filtered$`Year 1`)

ARG_Jaccard_pairwise_filtered <- ARG_Jaccard_pairwise_filtered %>%
  group_by(SID1, Sites, SeasonFull, `Group 1`) %>%
  summarise(
    value = mean(value)
  )

ARG_Jaccard_pairwise_filtered$GroupSites <- paste(ARG_Jaccard_pairwise_filtered$`Group 1`, ARG_Jaccard_pairwise_filtered$Sites)
ARG_Jaccard_pairwise_filtered$SeasonFull <- factor(ARG_Jaccard_pairwise_filtered$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))

Jacc_pairwise_compared <- compare_means(value ~ GroupSites, data = ARG_Jaccard_pairwise_filtered, group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")
Jacc_pairwise_compared <- Jacc_pairwise_compared[which(Jacc_pairwise_compared$p.adj < 0.05),]

ggplot(data = ARG_Jaccard_pairwise_filtered, aes(x = Sites, y = value, color = `Group 1`)) + geom_boxplot() + 
  facet_wrap(~ARG_Jaccard_pairwise_filtered$SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_test() + ylab("Jaccard distance") +
  ylim(0,1.2) +
  labs(title = "ARG Jaccard distances to cow samples", subtitle = "Wilcoxon rank sum test. Benjamini-Hochberg correction.") +
  scale_color_manual(values=c("#00BA38", "#619CFF"), name = "", labels = c("Office worker", "Farmer"))


###Running Maaslin2 to check for significant associations of ARGs and abx classes with farmers relative to office workers
maaslin_human_metadata <- human_metadata %>% dplyr::select(3,1, 12:16)
maaslin_human_metadata$SeasonFull <- paste(maaslin_human_metadata$Season, maaslin_human_metadata$Year)
row.names(maaslin_human_metadata) <- maaslin_human_metadata$Sample_ID
maaslin_shortbred_counts_wide_humans <- shortbred_counts_wide_humans[,colSums(shortbred_counts_wide_humans) > 0]
#maaslin_shortbred_counts_wide_humans <- maaslin_shortbred_counts_wide_humans[,colSums(maaslin_shortbred_counts_wide_humans != 0) > 150]
#maaslin_shortbred_counts_wide_humans <- rownames_to_column(maaslin_shortbred_counts_wide_humans, "Sample_ID")

Maaslin2(
  input_data = maaslin_shortbred_counts_wide_humans,
  input_metadata = maaslin_human_metadata,
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/ARG_farmerVSoffice",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)

Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_human_metadata,
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/AbxClasses_farmerVSoffice",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)

###Running Maaslin2 to check for significant associations of ARGs and abx classes with farmers within each season
maaslin_human_metadata_SPRING1 <- maaslin_human_metadata[which(maaslin_human_metadata$SeasonFull == "SPRING 1"),]
row.names(maaslin_human_metadata_SPRING1) <- maaslin_human_metadata_SPRING1$Sample_ID

#ARG
Maaslin2(
  input_data = maaslin_shortbred_counts_wide_humans,
  input_metadata = maaslin_human_metadata_SPRING1,
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/ARG_farmerVSoffice_SPRING1",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)


maaslin_human_metadata_SUMMER1 <- maaslin_human_metadata[which(maaslin_human_metadata$SeasonFull == "SUMMER 1"),]
row.names(maaslin_human_metadata_SUMMER1) <- maaslin_human_metadata_SUMMER1$Sample_ID

Maaslin2(
  input_data = maaslin_shortbred_counts_wide_humans,
  input_metadata = maaslin_human_metadata_SUMMER1,
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/ARG_farmerVSoffice_SUMMER1",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)


maaslin_human_metadata_FALL1 <- maaslin_human_metadata[which(maaslin_human_metadata$SeasonFull == "FALL 1"),]
row.names(maaslin_human_metadata_FALL1) <- maaslin_human_metadata_FALL1$Sample_ID

Maaslin2(
  input_data = maaslin_shortbred_counts_wide_humans,
  input_metadata = maaslin_human_metadata_FALL1,
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/ARG_farmerVSoffice_FALL1",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)


maaslin_human_metadata_WINTER1 <- maaslin_human_metadata[which(maaslin_human_metadata$SeasonFull == "WINTER 1"),]
row.names(maaslin_human_metadata_WINTER1) <- maaslin_human_metadata_WINTER1$Sample_ID

Maaslin2(
  input_data = maaslin_shortbred_counts_wide_humans,
  input_metadata = maaslin_human_metadata_WINTER1,
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/ARG_farmerVSoffice_WINTER1",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)


maaslin_human_metadata_SPRING2 <- maaslin_human_metadata[which(maaslin_human_metadata$SeasonFull == "SPRING 2"),]
row.names(maaslin_human_metadata_SPRING2) <- maaslin_human_metadata_SPRING2$Sample_ID

Maaslin2(
  input_data = maaslin_shortbred_counts_wide_humans,
  input_metadata = maaslin_human_metadata_SPRING2,
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/ARG_farmerVSoffice_SPRING2",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

#Abx classes
Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_human_metadata[which(maaslin_human_metadata$SeasonFull == "SPRING 1"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/AbxClasses_farmerVSoffice_SPRING1",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_human_metadata[which(maaslin_human_metadata$SeasonFull == "SUMMER 1"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/AbxClasses_farmerVSoffice_SUMMER1",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_human_metadata[which(maaslin_human_metadata$SeasonFull == "FALL 1"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/AbxClasses_farmerVSoffice_FALL1",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_human_metadata[which(maaslin_human_metadata$SeasonFull == "WINTER 1"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/AbxClasses_farmerVSoffice_WINTER1",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_human_metadata[which(maaslin_human_metadata$SeasonFull == "SPRING 2"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/AbxClasses_farmerVSoffice_SPRING2",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

###I identified the taxa of the shared strains between cow and farmer guts
###For the genera that were also enriched in farmer noses relative to office worker noses (16S analysis) ran a corelation analysis against ARGs
###Running Maaslin2 to check for significant associations of only significantly correlating ARGs (spearman coef > 0.5, pval < 0.05) with farmers relative to office workers
sharedGenera_ARGs <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/Lindsey's analysis/BManalysis/221221_Cow_Genus-ARG_correlation_SharedGenera_v2.csv",
                              header = TRUE)

Maaslin2(
  input_data = shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% rownames(maaslin_human_metadata)), c(sharedGenera_ARGs$SharedTaxaCorrelGenes)],
  input_metadata = maaslin_human_metadata,
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/SharedTaxaARG_farmerVSoffice",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)


###Running Maaslin2 to check for significant associations of ARGs and abx classes with cows relative to office workers in general and within seasons
metadata_CowsOffice <- samples_metadata[which(samples_metadata$Group != "W"),]
maaslin_CD_metadata <- metadata_CowsOffice %>% dplyr::select(3,1, 12,13, 14, 15, 16)
row.names(maaslin_CD_metadata) <- maaslin_CD_metadata$Sample_ID
maaslin_CD_metadata$SeasonFull <- paste(maaslin_CD_metadata$Season, maaslin_CD_metadata$Year)

#ARGs
Maaslin2(
  input_data = shortbred_counts_wide,
  input_metadata = maaslin_CD_metadata,
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/ARG_cowVSoffice",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)

Maaslin2(
  input_data = shortbred_counts_wide,
  input_metadata = maaslin_CD_metadata[which(maaslin_CD_metadata$Season == "SPRING" & maaslin_CD_metadata$Year == 1),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/ARG_cowVSoffice_SPRING1",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

Maaslin2(
  input_data = shortbred_counts_wide,
  input_metadata = maaslin_CD_metadata[which(maaslin_CD_metadata$Season == "SUMMER" & maaslin_CD_metadata$Year == 1),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/ARG_cowVSoffice_SUMMER1",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

Maaslin2(
  input_data = shortbred_counts_wide,
  input_metadata = maaslin_CD_metadata[which(maaslin_CD_metadata$Season == "FALL" & maaslin_CD_metadata$Year == 1),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/ARG_cowVSoffice_FALL1",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

Maaslin2(
  input_data = shortbred_counts_wide,
  input_metadata = maaslin_CD_metadata[which(maaslin_CD_metadata$Season == "WINTER" & maaslin_CD_metadata$Year == 1),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/ARG_cowVSoffice_WINTER1",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

Maaslin2(
  input_data = shortbred_counts_wide,
  input_metadata = maaslin_CD_metadata[which(maaslin_CD_metadata$Season == "SPRING" & maaslin_CD_metadata$Year == 2),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/ARG_cowVSoffice_SPRING2",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

#Abx classes
Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_CD_metadata,
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/AbxClass_cowVSoffice",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)

Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_CD_metadata[which(maaslin_CD_metadata$Season == "SPRING" & maaslin_CD_metadata$Year == 1),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/AbxClass_cowVSoffice_SPRING1",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_CD_metadata[which(maaslin_CD_metadata$Season == "SUMMER" & maaslin_CD_metadata$Year == 1),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/AbxClass_cowVSoffice_SUMMER1",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_CD_metadata[which(maaslin_CD_metadata$Season == "FALL" & maaslin_CD_metadata$Year == 1),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/AbxClass_cowVSoffice_FALL1",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_CD_metadata[which(maaslin_CD_metadata$Season == "WINTER" & maaslin_CD_metadata$Year == 1),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/AbxClass_cowVSoffice_WINTER1",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_CD_metadata[which(maaslin_CD_metadata$Season == "SPRING" & maaslin_CD_metadata$Year == 2),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/AbxClass_cowVSoffice_SPRING2",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

###Identified ARGs that are enriched in cows relative to humans
###Running Maaslin2 to check for significant associations of these ARGs with farmers relative to office workers
#cowARGs <- read_csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/ARG_cowVSoffice/221222_SigEnrichedGenes_v3.csv")
cowARGs <- read_csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/ARG_cowVShuman/230217_SigEnrichedGenes_v3.csv")

Maaslin2(
  input_data = shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% rownames(maaslin_human_metadata)), c(cowARGs$OGname)],
  input_metadata = maaslin_human_metadata,
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/CowEnrichedARGs_farmerVSoffice",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)

###Looking at the abx classes of the ARGs enriched in cows relative to office workers
cowARGs <- merge(cowARGs, ARG_mapping, by.x = "OGname", by.y = "Gene")

shortbred_cowARGs <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %notin% rownames(maaslin_human_metadata)), c(cowARGs$OGname)]
prevalence_cowARGs <- as.data.frame(colSums(shortbred_cowARGs != 0))
colnames(prevalence_cowARGs)[1] <- "Prevalence_Count"

cowARGs <- merge(cowARGs, prevalence_cowARGs, by.x = "OGname", by.y = "row.names")
cowARGs$Prevalence_Frac <- cowARGs$Prevalence_Count/dim(shortbred_cowARGs)[1]

cowARGs_maaslin <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/ARG_cowVShuman/significant_results.tsv", header = TRUE, "\t")
cowARGs_maaslin <- cowARGs_maaslin[which(cowARGs_maaslin$coef < 0 & cowARGs_maaslin$qval < 0.05),]
cowARGs_maaslin <- cowARGs_maaslin[,c(1,4)]
cowARGs_maaslin$coef <- cowARGs_maaslin$coef * -1

cowARGs <- merge(cowARGs, cowARGs_maaslin, by.x = "M2name", by.y = "feature")
#write_csv(cowARGs, file = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/ARG_cowVSoffice/221222_SigEnrichedGenes_v4.csv")

colors_ARGclasses <- data_frame(ARGclass = c("TETRACYCLINE", "TRIMETHOPRIM", "GLYCOPEPTIDE", "BETA-LACTAM", "LINCOSAMIDE", "MACROLIDE", "PHENICOL", "QUINOLONE", "AMINOGLYCOSIDE", "FOSFOMYCIN", "Other"),
                                Color = c(c(natparks.pals("DeathValley", 11))[1:10], "#bcbcbc"))

cowARGs$FormattedAbxClass <- ""

for (i in seq(1, length(cowARGs$FormattedAbxClass), 1)) {
  if (cowARGs$AbxClass[i] %in% colors_ARGclasses$ARGclass) {
    cowARGs$FormattedAbxClass[i] <- cowARGs$AbxClass[i]
  } else {
    cowARGs$FormattedAbxClass[i] <- "Other"
  }
}

cowARGs <- merge(cowARGs, colors_ARGclasses, by.x = "FormattedAbxClass", by.y = "ARGclass")

cowARGs$FormattedAbxClass <- factor(cowARGs$FormattedAbxClass, levels = c("TETRACYCLINE", "TRIMETHOPRIM", "GLYCOPEPTIDE", "BETA-LACTAM", "LINCOSAMIDE", "MACROLIDE", "PHENICOL", "QUINOLONE", "AMINOGLYCOSIDE", "FOSFOMYCIN", "Other"))

ggplot(cowARGs) + geom_point(aes(x = coef, y = Prevalence_Frac*100, group = FormattedAbxClass, color = FormattedAbxClass)) + coord_trans(y = "log10") +
  scale_color_manual(name = "", values = unique(cowARGs$Color), labels = unique(cowARGs$FormattedAbxClass)) +
  theme_test() + theme(legend.position = c(0.85, 0.25))
ggplot(cowARGs, aes(x = coef, y = Prevalence_Frac*100)) + coord_trans(y = "log10") + geom_jitter()
ggplot(cowARGs, aes(x = coef, y = Prevalence_Frac*100)) + coord_trans(y = "log10") + geom_point()

ggplot(cowARGs) + geom_point(aes(x = coef, y = Prevalence_Frac*100, color = FormattedAbxClass)) +
  scale_color_manual(values = c(natparks.pals("DeathValley", 10), "#848484")) +
  theme_linedraw() + ylim(0, NA) +
  theme(legend.position = c(0.85, 0.25))
ggplot(cowARGs) + geom_point(aes(x = coef, y = Prevalence_Frac*100, fill = FormattedAbxClass), size = 4, alpha = 1, pch = 21, color = "black") +
  scale_fill_manual(values = c(natparks.pals("DeathValley", 10), "#848484")) +
  theme_linedraw() + ylim(0, NA) +
  labs(color = "Antibiotic class") +
  theme(legend.position = "none", 
        panel.grid.major.x = element_line(size= 0.2, color = "gray", linetype = "dashed"), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(size= 0.2, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        axis.text = element_blank(), axis.title = element_blank())
p <- ggplot(cowARGs) + geom_point(aes(x = coef, y = Prevalence_Frac*100, fill = FormattedAbxClass), size = 4, alpha = 1, pch = 21, color = "black") +
  scale_fill_manual(values = c("#8C2B0E", "#B25422", "#D8813B", "#FEB359", "#9F7E59", "#233F6C", "#435F90", "#5B4C64", "#81565F", "#B47E83", "#848484")) +
  theme_linedraw() + ylim(0, NA) +
  labs(color = "Antibiotic class") +
  theme(legend.position = "none", 
        panel.grid.major.x = element_line(size= 0.2, color = "gray", linetype = "dashed"), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(size= 0.2, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        axis.text = element_blank(), axis.title = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Shortbred_Maaslin2_CowARGs_v4.svg", p, bg = "transparent", width = 10, height = 10)

ggplot(cowARGs[which(cowARGs$Prevalence_Frac < 0.5 & cowARGs$coef < 2.5),]) + geom_point(aes(x = coef, y = Prevalence_Frac*100, fill = FormattedAbxClass), size = 4, alpha = 0.8, pch = 21, color = "black") +
  scale_fill_manual(values = c(natparks.pals("DeathValley", 10), "#848484")) +
  theme_linedraw() +
  labs(color = "Antibiotic class") +
  scale_x_continuous(breaks = seq(0, 3, by = 1)) +
  scale_y_continuous(breaks = seq(0, 60, by = 20)) +
  theme(legend.position = "none", 
        panel.grid.major.x = element_line(size= 0.2, color = "gray", linetype = "dashed"), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(size= 0.2, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        axis.title = element_blank())
p <- ggplot(cowARGs[which(cowARGs$Prevalence_Frac < 0.4 & cowARGs$coef < 2.1),]) + geom_point(aes(x = coef, y = Prevalence_Frac*100, fill = FormattedAbxClass), size = 4, alpha = 0.8, pch = 21, color = "black") +
  scale_fill_manual(values = c("#8C2B0E", "#B25422", "#D8813B", "#FEB359", "#9F7E59", "#233F6C", "#435F90", "#5B4C64", "#81565F", "#B47E83", "#848484")) +
  theme_linedraw() +
  labs(color = "Antibiotic class") +
  scale_x_continuous(breaks = seq(0.5, 2.5, by = 1)) +
  scale_y_continuous(breaks = seq(0, 60, by = 20)) +
  theme(legend.position = "none", 
        panel.grid.major.x = element_line(size= 0.2, color = "gray", linetype = "dashed"), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(size= 0.2, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        axis.title = element_blank(), axis.text = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Shortbred_Maaslin2_CowARGs_inlet_v3.svg", p, bg = "transparent", width = 4, height = 4)

cowARGs$Count <- 1

cowARGs_ClassTotals <- cowARGs[,c(1,9)] %>% group_by(FormattedAbxClass) %>% summarise_all(sum)
cowARGs_ClassTotals <- merge(cowARGs_ClassTotals, colors_ARGclasses, by.x = "FormattedAbxClass", by.y = "ARGclass")
cowARGs_ClassTotals$FormattedAbxClass <- factor(cowARGs_ClassTotals$FormattedAbxClass, levels = c("TETRACYCLINE", "TRIMETHOPRIM", "GLYCOPEPTIDE", "BETA-LACTAM", "LINCOSAMIDE", "MACROLIDE", "PHENICOL", "QUINOLONE", "AMINOGLYCOSIDE", "FOSFOMYCIN", "Other"))
#cowARGs_ClassTotals$Color <- factor(cowARGs_ClassTotals$Color, levels = c("#8C2B0E", "#AE5020", "#D07735", "#F2A450", "#9F7E59", "#132F5B", "#2F4B7A", "#4A5982", "#60485B", "#865A63", "#bcbcbc"))

cowARGs_ClassTotals <- rbind(cowARGs_ClassTotals[10,], cowARGs_ClassTotals[11,], cowARGs_ClassTotals[4,], cowARGs_ClassTotals[2,], cowARGs_ClassTotals[5,], cowARGs_ClassTotals[6,],
      cowARGs_ClassTotals[8,], cowARGs_ClassTotals[9,], cowARGs_ClassTotals[1,], cowARGs_ClassTotals[3,], cowARGs_ClassTotals[7,])

ggplot(cowARGs_ClassTotals, aes(x = "", y = Count, fill = FormattedAbxClass), fill = "transparent") + 
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0, direction = -1) + theme_void() +
  scale_fill_manual(name = "", values = unique(cowARGs_ClassTotals$Color), labels = unique(cowARGs_ClassTotals$FormattedAbxClass)) +
  theme(legend.position = "none")
ggplot(cowARGs_ClassTotals, aes(x = "", y = Count, fill = FormattedAbxClass), fill = "transparent") + 
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0, direction = -1) + theme_void() +
  scale_fill_manual(values = c(natparks.pals("DeathValley", 10), "#848484")) + 
  theme(legend.position = "none")
p <- ggplot(cowARGs_ClassTotals, aes(x = "", y = Count, fill = FormattedAbxClass), fill = "transparent") + 
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0, direction = -1) + theme_void() +
  scale_fill_manual(values = c("#8C2B0E", "#B25422", "#D8813B", "#FEB359", "#9F7E59", "#233F6C", "#435F90", "#5B4C64", "#81565F", "#B47E83", "#848484")) + 
  theme(legend.position = "none")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Shortbred_Maaslin2_CowARGs_Breakdown_v2.svg", p, bg = "transparent", width = 5, height = 5)

cowARGs$ARGsource <- ""

for (i in seq(1, length(cowARGs$FormattedAbxClass), 1)) {
  if (strsplit(cowARGs$OGname[i], split = character(0))[[1]][1] == "D") {
    cowARGs$ARGsource[i] <- "FMG"
  } else {
    cowARGs$ARGsource[i] <- "Non-FMG"
  }
}

cowARG_ClassSourceTotal <- cowARGs[,c(1,9,10)] %>% group_by(FormattedAbxClass, ARGsource) %>% summarise_all(sum)
cowARG_ClassSourceTotal$Color <- ""

for (i in seq(1, dim(cowARG_ClassSourceTotal)[1], 1)) {
  if (cowARG_ClassSourceTotal$ARGsource[i] == "FMG") {
    cowARG_ClassSourceTotal$Color[i] <- "#151514"
  } else {
    cowARG_ClassSourceTotal$Color[i] <- "#f9f9f9"
  }
}

cowARG_ClassSourceTotal$ClassSource <- paste(cowARG_ClassSourceTotal$FormattedAbxClass, cowARG_ClassSourceTotal$ARGsource)
cowARG_ClassSourceTotal$ClassSource <- factor(cowARG_ClassSourceTotal$ClassSource, levels = c("TETRACYCLINE FMG", "TETRACYCLINE Non-FMG", "TRIMETHOPRIM FMG",
                                                                                              "GLYCOPEPTIDE FMG", "BETA-LACTAM FMG", "BETA-LACTAM Non-FMG",
                                                                                              "LINCOSAMIDE FMG", "LINCOSAMIDE Non-FMG", "MACROLIDE FMG",
                                                                                              "MACROLIDE Non-FMG", "PHENICOL FMG", "QUINOLONE FMG",
                                                                                              "AMINOGLYCOSIDE FMG", "FOSFOMYCIN FMG", "Other FMG"))

p <- ggplot(cowARG_ClassSourceTotal, aes(x = "", y = Count, fill = ClassSource)) + 
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0, direction = -1) + theme_void() +
  scale_fill_manual(name = "", values = cowARG_ClassSourceTotal$Color, labels = cowARG_ClassSourceTotal$ARGsource) + 
  theme(legend.position = "none")
p
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Shortbred_Maaslin2_CowARGs_Source_v1.svg", p, bg = "transparent", width = 5, height = 5)

###Running Maaslin2 to check for significant associations of ARGs and abx classes with cows relative to humans (farmers and office workers)
maaslin_metadata <- samples_metadata %>% dplyr::select(3, 1, 12:16)
maaslin_metadata$Origin <- ""

for (i in seq(1, length(maaslin_metadata$Sample_ID), 1)) {
  if (maaslin_metadata$Group[i] == "C") {
    maaslin_metadata$Origin[i] <- "C"
  } else {
    maaslin_metadata$Origin[i] <- "H"
  }
}

maaslin_metadata$SeasonFull <- paste(maaslin_metadata$Season, maaslin_metadata$Year)

row.names(maaslin_metadata) <- maaslin_metadata$Sample_ID

Maaslin2(
  input_data = shortbred_counts_wide,
  input_metadata = maaslin_metadata,
  fixed_effects = "Origin",
  random_effects = c("SID", "SeasonFull", "Site"),
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/ARG_cowVShuman"
)

Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_metadata,
  fixed_effects = "Origin",
  random_effects = c("SID", "SeasonFull", "Site"),
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/AbxClass_cowVShuman"
)

M2_AbxClass_cowVShuman <- read_delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/AbxClass_cowVShuman/significant_results.tsv",
                                     "\t", escape_double = FALSE, trim_ws = TRUE)
M2_AbxClass_cowVShuman$shape <- ifelse(M2_AbxClass_cowVShuman$coef > 0, 1, -1)
M2_AbxClass_cowVShuman <- subset(M2_AbxClass_cowVShuman, qval < 0.05)
M2_AbxClass_cowVShuman$CowCoef <- M2_AbxClass_cowVShuman$coef * -1
ggplot(M2_AbxClass_cowVShuman, aes(x=coef, y=reorder(feature, coef), color=coef)) +
  geom_point(aes(shape=factor(shape)), size=2)+
  labs(x="Relative to Cow Gut Resistome",
       y="Antibiotic class", shape="Shape")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_vline(xintercept=0, colour='black',size=1, linetype="dotted")+
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0)) + 
  scale_colour_gradient2()

ggplot(M2_AbxClass_cowVShuman[which(M2_AbxClass_cowVShuman$feature %in% unique(AbxClass_plot$AbxClass)),], aes(x=CowCoef, y=reorder(feature, CowCoef))) +
  geom_point(color = "#873e23", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=CowCoef + stderr, xmin =CowCoef - stderr, height = 0), color = "#873e23", size = 0.7) +
  geom_vline(xintercept=0, colour='black',size=2, linetype="dotted") +
  labs(x="Relative to Human Gut Resistome",
       y="Antibiotic class") +
  theme_linedraw()
p <- ggplot(M2_AbxClass_cowVShuman[which(M2_AbxClass_cowVShuman$feature %in% unique(AbxClass_plot$AbxClass)),], aes(x=CowCoef, y=reorder(feature, CowCoef))) +
  geom_point(color = "#873e23", shape = 16, size = 9, alpha = 0.8) +
  geom_errorbarh(aes(xmax=CowCoef + stderr, xmin =CowCoef - stderr, height = 0), color = "#873e23", size = 2) +
  geom_vline(xintercept=0, colour='black',size=2, linetype="dotted") +
  labs(x="Relative to Human Gut Resistome",
       y="Antibiotic class") +
  theme_linedraw() +
  theme(axis.title = element_blank(), axis.text = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Shortbred_Maaslin2_AbxClasses_CowHuman_v4.svg", p, bg = "transparent", width = 5, height = 8)


###Running Maaslin2 to check for significant associations of ARGs and abx classes with cows/farmers relative to office workers within collection site 4
###Site 4 has only 1 farmer and 1 office worker. Can't do Maaslin on this
maaslin_site4_metadata <- samples_metadata[which(samples_metadata$Site == "4"),]
maaslin_site4_metadata <- maaslin_site4_metadata %>% select(3, 1, 12:16)
maaslin_site4_metadata$SeasonFull <- paste(maaslin_site4_metadata$Season, maaslin_site4_metadata$Year)
row.names(maaslin_site4_metadata) <- maaslin_site4_metadata$Sample_ID

###Running Maaslin2 to check for significant associations of ARGs and abx classes with cows/farmers relative to office workers within collection site 5
###Site 5 has only 1 farmer and 1 office worker. Can't do Maaslin on this
maaslin_site5_metadata <- samples_metadata[which(samples_metadata$Site == "5"),]
maaslin_site5_metadata <- maaslin_site5_metadata %>% select(3, 1, 12:16)
maaslin_site5_metadata$SeasonFull <- paste(maaslin_site5_metadata$Season, maaslin_site5_metadata$Year)
row.names(maaslin_site5_metadata) <- maaslin_site5_metadata$Sample_ID

###Running Maaslin2 to check for significant associations of ARGs and abx classes with cows/farmers relative to office workers within collection site 20
###Site 20 has only 1 farmer and 1 office worker. Can't do Maaslin on this
maaslin_site20_metadata <- samples_metadata[which(samples_metadata$Site == "20"),]
maaslin_site20_metadata <- maaslin_site20_metadata %>% select(3, 1, 12:16)
maaslin_site20_metadata$SeasonFull <- paste(maaslin_site20_metadata$Season, maaslin_site20_metadata$Year)
row.names(maaslin_site20_metadata) <- maaslin_site20_metadata$Sample_ID

###Running Maaslin2 to check for significant associations of ARGs and abx classes with cows/farmers relative to office workers within collection site 37
###Site 37 has only 1 farmer and 1 office worker. Can't do Maaslin on this
maaslin_site37_metadata <- samples_metadata[which(samples_metadata$Site == "37"),]
maaslin_site37_metadata <- maaslin_site37_metadata %>% select(3, 1, 12:16)
maaslin_site37_metadata$SeasonFull <- paste(maaslin_site37_metadata$Season, maaslin_site37_metadata$Year)
row.names(maaslin_site37_metadata) <- maaslin_site37_metadata$Sample_ID

###Running Maaslin2 to check for significant associations of ARGs and abx classes with cows/farmers relative to office workers within collection site 2
maaslin_site2_metadata <- samples_metadata[which(samples_metadata$Site == "2"),]
maaslin_site2_metadata <- maaslin_site2_metadata %>% select(3, 1, 12:16)
maaslin_site2_metadata$SeasonFull <- paste(maaslin_site2_metadata$Season, maaslin_site2_metadata$Year)
row.names(maaslin_site2_metadata) <- maaslin_site2_metadata$Sample_ID

Maaslin2(
  input_data = shortbred_counts_wide,
  input_metadata = maaslin_site2_metadata[which(maaslin_site2_metadata$Group != "W"),],
  fixed_effects = "Group",
  random_effects = c("SID", "SeasonFull"),
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/site2_ARG_cowVSoffice"
)

Maaslin2(
  input_data = shortbred_counts_wide,
  input_metadata = maaslin_site2_metadata[which(maaslin_site2_metadata$Group != "C"),],
  fixed_effects = "Group",
  random_effects = c("SID", "SeasonFull"),
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/site2_ARG_farmerVSoffice"
)

Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_site2_metadata[which(maaslin_site2_metadata$Group != "W"),],
  fixed_effects = "Group",
  random_effects = c("SID", "SeasonFull"),
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/site2_AbxClasses_cowVSoffice"
)

Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_site2_metadata[which(maaslin_site2_metadata$Group != "C"),],
  fixed_effects = "Group",
  random_effects = c("SID", "SeasonFull"),
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/site2_AbxClasses_farmerVSoffice"
)

Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_site2_metadata[which(maaslin_site2_metadata$Group != "D"),],
  fixed_effects = "Group",
  random_effects = c("SID", "SeasonFull"),
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/site2_AbxClasses_cowVSfarmer"
)

###Running Maaslin2 to check for significant associations of ARGs and abx classes with cows/farmers relative to office workers within collection site 10
maaslin_site10_metadata <- samples_metadata[which(samples_metadata$Site == "10"),]
maaslin_site10_metadata <- maaslin_site10_metadata %>% select(3, 1, 12:16)
maaslin_site10_metadata$SeasonFull <- paste(maaslin_site10_metadata$Season, maaslin_site10_metadata$Year)
row.names(maaslin_site10_metadata) <- maaslin_site10_metadata$Sample_ID

Maaslin2(
  input_data = shortbred_counts_wide,
  input_metadata = maaslin_site10_metadata[which(maaslin_site10_metadata$Group != "W"),],
  fixed_effects = "Group",
  random_effects = c("SID", "SeasonFull"),
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/site10_ARG_cowVSoffice"
)

Maaslin2(
  input_data = shortbred_counts_wide,
  input_metadata = maaslin_site10_metadata[which(maaslin_site10_metadata$Group != "C"),],
  fixed_effects = "Group",
  random_effects = c("SID", "SeasonFull"),
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/site10_ARG_farmerVSoffice"
)

Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_site10_metadata[which(maaslin_site10_metadata$Group != "W"),],
  fixed_effects = "Group",
  random_effects = c("SID", "SeasonFull"),
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/site10_AbxClasses_cowVSoffice"
)

Maaslin2(
  input_data = shortbred_AbxClass,
  input_metadata = maaslin_site10_metadata[which(maaslin_site10_metadata$Group != "C"),],
  fixed_effects = "Group",
  random_effects = c("SID", "SeasonFull"),
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED_Maaslin2/site10_AbxClasses_farmerVSoffice"
)

###Looking for ARGs that are uniquely present in cows and farmers (not in office workers) of site 2
shortbred_site2_cows <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% rownames(maaslin_site2_metadata[which(maaslin_site2_metadata$Group == "C"),])),]
shortbred_site2_cows <- shortbred_site2_cows[,which(colSums(shortbred_site2_cows) > 0)]

shortbred_site2_office <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% rownames(maaslin_site2_metadata[which(maaslin_site2_metadata$Group == "D"),])),]
shortbred_site2_office <- shortbred_site2_office[,which(colSums(shortbred_site2_office) > 0)]

shortbred_site2_farmer <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% rownames(maaslin_site2_metadata[which(maaslin_site2_metadata$Group == "W"),])),]
shortbred_site2_farmer <- shortbred_site2_farmer[,which(colSums(shortbred_site2_farmer) > 0)]

site2_cowANDfarmer_ARGs <- colnames(shortbred_site2_cows)
site2_cowANDfarmer_ARGs <- site2_cowANDfarmer_ARGs[which(site2_cowANDfarmer_ARGs %notin% colnames(shortbred_site2_office))]
site2_cowANDfarmer_ARGs <- site2_cowANDfarmer_ARGs[which(site2_cowANDfarmer_ARGs %in% colnames(shortbred_site2_farmer))]

site2_cowANDoffice_ARGs <- colnames(shortbred_site2_cows)
site2_cowANDoffice_ARGs <- site2_cowANDoffice_ARGs[which(site2_cowANDoffice_ARGs %notin% colnames(shortbred_site2_farmer))]
site2_cowANDoffice_ARGs <- site2_cowANDoffice_ARGs[which(site2_cowANDoffice_ARGs %in% colnames(shortbred_site2_office))]

site2_cowANDfarmerANDoffice_ARGs <- colnames(shortbred_site2_cows)
site2_cowANDfarmerANDoffice_ARGs <- site2_cowANDfarmerANDoffice_ARGs[which(site2_cowANDfarmerANDoffice_ARGs %in% colnames(shortbred_site2_office))]
site2_cowANDfarmerANDoffice_ARGs <- site2_cowANDfarmerANDoffice_ARGs[which(site2_cowANDfarmerANDoffice_ARGs %in% colnames(shortbred_site2_farmer))]

###Looking for ARGs that are uniquely present in cows and farmers (not in office workers) of site 10
shortbred_site10_cows <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% rownames(maaslin_site10_metadata[which(maaslin_site10_metadata$Group == "C"),])),]
shortbred_site10_cows <- shortbred_site10_cows[,which(colSums(shortbred_site10_cows) > 0)]

shortbred_site10_office <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% rownames(maaslin_site10_metadata[which(maaslin_site10_metadata$Group == "D"),])),]
shortbred_site10_office <- shortbred_site10_office[,which(colSums(shortbred_site10_office) > 0)]

shortbred_site10_farmer <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% rownames(maaslin_site10_metadata[which(maaslin_site10_metadata$Group == "W"),])),]
shortbred_site10_farmer <- shortbred_site10_farmer[,which(colSums(shortbred_site10_farmer) > 0)]

site10_cowANDfarmer_ARGs <- colnames(shortbred_site10_cows)
site10_cowANDfarmer_ARGs <- site10_cowANDfarmer_ARGs[which(site10_cowANDfarmer_ARGs %notin% colnames(shortbred_site10_office))]
site10_cowANDfarmer_ARGs <- site10_cowANDfarmer_ARGs[which(site10_cowANDfarmer_ARGs %in% colnames(shortbred_site10_farmer))]

site10_cowANDoffice_ARGs <- colnames(shortbred_site10_cows)
site10_cowANDoffice_ARGs <- site10_cowANDoffice_ARGs[which(site10_cowANDoffice_ARGs %notin% colnames(shortbred_site10_farmer))]
site10_cowANDoffice_ARGs <- site10_cowANDoffice_ARGs[which(site10_cowANDoffice_ARGs %in% colnames(shortbred_site10_office))]

site10_cowANDfarmerANDoffice_ARGs <- colnames(shortbred_site10_cows)
site10_cowANDfarmerANDoffice_ARGs <- site10_cowANDfarmerANDoffice_ARGs[which(site10_cowANDfarmerANDoffice_ARGs %in% colnames(shortbred_site10_office))]
site10_cowANDfarmerANDoffice_ARGs <- site10_cowANDfarmerANDoffice_ARGs[which(site10_cowANDfarmerANDoffice_ARGs %in% colnames(shortbred_site10_farmer))]

###Looking for ARGs that are uniquely present in cows and farmers (not in office workers) of site 23
shortbred_site23_cows <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% samples_metadata[which(samples_metadata$Site == "23" & samples_metadata$Group == "C"),]$Sample_ID),]
shortbred_site23_cows <- shortbred_site23_cows[,which(colSums(shortbred_site23_cows) > 0)]

shortbred_site23_office <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% samples_metadata[which(samples_metadata$Site == "23" & samples_metadata$Group == "D"),]$Sample_ID),]
shortbred_site23_office <- shortbred_site23_office[,which(colSums(shortbred_site23_office) > 0)]

shortbred_site23_farmer <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% samples_metadata[which(samples_metadata$Site == "23" & samples_metadata$Group == "W"),]$Sample_ID),]
shortbred_site23_farmer <- shortbred_site23_farmer[,which(colSums(shortbred_site23_farmer) > 0)]

site23_cowANDfarmer_ARGs <- colnames(shortbred_site23_cows)
site23_cowANDfarmer_ARGs <- site23_cowANDfarmer_ARGs[which(site23_cowANDfarmer_ARGs %notin% colnames(shortbred_site23_office))]
site23_cowANDfarmer_ARGs <- site23_cowANDfarmer_ARGs[which(site23_cowANDfarmer_ARGs %in% colnames(shortbred_site23_farmer))]

site23_cowANDoffice_ARGs <- colnames(shortbred_site23_cows)
site23_cowANDoffice_ARGs <- site23_cowANDoffice_ARGs[which(site23_cowANDoffice_ARGs %notin% colnames(shortbred_site23_farmer))]
site23_cowANDoffice_ARGs <- site23_cowANDoffice_ARGs[which(site23_cowANDoffice_ARGs %in% colnames(shortbred_site23_office))]

site23_cowANDfarmerANDoffice_ARGs <- colnames(shortbred_site23_cows)
site23_cowANDfarmerANDoffice_ARGs <- site23_cowANDfarmerANDoffice_ARGs[which(site23_cowANDfarmerANDoffice_ARGs %in% colnames(shortbred_site23_office))]
site23_cowANDfarmerANDoffice_ARGs <- site23_cowANDfarmerANDoffice_ARGs[which(site23_cowANDfarmerANDoffice_ARGs %in% colnames(shortbred_site23_farmer))]

###Looking for ARGs that are uniquely present in cows and farmers (not in office workers) of site 26
shortbred_site26_cows <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% samples_metadata[which(samples_metadata$Site == "26" & samples_metadata$Group == "C"),]$Sample_ID),]
shortbred_site26_cows <- shortbred_site26_cows[,which(colSums(shortbred_site26_cows) > 0)]

shortbred_site26_office <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% samples_metadata[which(samples_metadata$Site == "26" & samples_metadata$Group == "D"),]$Sample_ID),]
shortbred_site26_office <- shortbred_site26_office[,which(colSums(shortbred_site26_office) > 0)]

shortbred_site26_farmer <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% samples_metadata[which(samples_metadata$Site == "26" & samples_metadata$Group == "W"),]$Sample_ID),]
shortbred_site26_farmer <- shortbred_site26_farmer[,which(colSums(shortbred_site26_farmer) > 0)]

site26_cowANDfarmer_ARGs <- colnames(shortbred_site26_cows)
site26_cowANDfarmer_ARGs <- site26_cowANDfarmer_ARGs[which(site26_cowANDfarmer_ARGs %notin% colnames(shortbred_site26_office))]
site26_cowANDfarmer_ARGs <- site26_cowANDfarmer_ARGs[which(site26_cowANDfarmer_ARGs %in% colnames(shortbred_site26_farmer))]

site26_cowANDoffice_ARGs <- colnames(shortbred_site26_cows)
site26_cowANDoffice_ARGs <- site26_cowANDoffice_ARGs[which(site26_cowANDoffice_ARGs %notin% colnames(shortbred_site26_farmer))]
site26_cowANDoffice_ARGs <- site26_cowANDoffice_ARGs[which(site26_cowANDoffice_ARGs %in% colnames(shortbred_site26_office))]

site26_cowANDfarmerANDoffice_ARGs <- colnames(shortbred_site26_cows)
site26_cowANDfarmerANDoffice_ARGs <- site26_cowANDfarmerANDoffice_ARGs[which(site26_cowANDfarmerANDoffice_ARGs %in% colnames(shortbred_site26_office))]
site26_cowANDfarmerANDoffice_ARGs <- site26_cowANDfarmerANDoffice_ARGs[which(site26_cowANDfarmerANDoffice_ARGs %in% colnames(shortbred_site26_farmer))]

##### Looking for ARGs that are uniquely present in cows and farmers (not in office workers) of site 16
shortbred_site16_cows <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% samples_metadata[which(samples_metadata$Site == "16" & samples_metadata$Group == "C"),]$Sample_ID),]
shortbred_site16_cows <- shortbred_site16_cows[,which(colSums(shortbred_site16_cows) > 0)]

shortbred_site16_office <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% samples_metadata[which(samples_metadata$Site == "16" & samples_metadata$Group == "D"),]$Sample_ID),]
shortbred_site16_office <- shortbred_site16_office[,which(colSums(shortbred_site16_office) > 0)]

shortbred_site16_farmer <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% samples_metadata[which(samples_metadata$Site == "16" & samples_metadata$Group == "W"),]$Sample_ID),]
shortbred_site16_farmer <- shortbred_site16_farmer[,which(colSums(shortbred_site16_farmer) > 0)]

site16_cowANDfarmer_ARGs <- colnames(shortbred_site16_cows)
site16_cowANDfarmer_ARGs <- site16_cowANDfarmer_ARGs[which(site16_cowANDfarmer_ARGs %notin% colnames(shortbred_site16_office))]
site16_cowANDfarmer_ARGs <- site16_cowANDfarmer_ARGs[which(site16_cowANDfarmer_ARGs %in% colnames(shortbred_site16_farmer))]

site16_cowANDoffice_ARGs <- colnames(shortbred_site16_cows)
site16_cowANDoffice_ARGs <- site16_cowANDoffice_ARGs[which(site16_cowANDoffice_ARGs %notin% colnames(shortbred_site16_farmer))]
site16_cowANDoffice_ARGs <- site16_cowANDoffice_ARGs[which(site16_cowANDoffice_ARGs %in% colnames(shortbred_site16_office))]

site16_cowANDfarmerANDoffice_ARGs <- colnames(shortbred_site16_cows)
site16_cowANDfarmerANDoffice_ARGs <- site16_cowANDfarmerANDoffice_ARGs[which(site16_cowANDfarmerANDoffice_ARGs %in% colnames(shortbred_site16_office))]
site16_cowANDfarmerANDoffice_ARGs <- site16_cowANDfarmerANDoffice_ARGs[which(site16_cowANDfarmerANDoffice_ARGs %in% colnames(shortbred_site16_farmer))]

###Looking for ARGs that are uniquely present in cows and farmers (not in office workers) of site 19
shortbred_site19_cows <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% samples_metadata[which(samples_metadata$Site == "19" & samples_metadata$Group == "C"),]$Sample_ID),]
shortbred_site19_cows <- shortbred_site19_cows[,which(colSums(shortbred_site19_cows) > 0)]

shortbred_site19_office <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% samples_metadata[which(samples_metadata$Site == "19" & samples_metadata$Group == "D"),]$Sample_ID),]
shortbred_site19_office <- shortbred_site19_office[,which(colSums(shortbred_site19_office) > 0)]

shortbred_site19_farmer <- shortbred_counts_wide[which(rownames(shortbred_counts_wide) %in% samples_metadata[which(samples_metadata$Site == "19" & samples_metadata$Group == "W"),]$Sample_ID),]
shortbred_site19_farmer <- shortbred_site19_farmer[,which(colSums(shortbred_site19_farmer) > 0)]

site19_cowANDfarmer_ARGs <- colnames(shortbred_site19_cows)
site19_cowANDfarmer_ARGs <- site19_cowANDfarmer_ARGs[which(site19_cowANDfarmer_ARGs %notin% colnames(shortbred_site19_office))]
site19_cowANDfarmer_ARGs <- site19_cowANDfarmer_ARGs[which(site19_cowANDfarmer_ARGs %in% colnames(shortbred_site19_farmer))]

site19_cowANDoffice_ARGs <- colnames(shortbred_site16_cows)
site19_cowANDoffice_ARGs <- site19_cowANDoffice_ARGs[which(site19_cowANDoffice_ARGs %notin% colnames(shortbred_site19_farmer))]
site19_cowANDoffice_ARGs <- site19_cowANDoffice_ARGs[which(site19_cowANDoffice_ARGs %in% colnames(shortbred_site19_office))]

site19_cowANDfarmerANDoffice_ARGs <- colnames(shortbred_site19_cows)
site19_cowANDfarmerANDoffice_ARGs <- site19_cowANDfarmerANDoffice_ARGs[which(site19_cowANDfarmerANDoffice_ARGs %in% colnames(shortbred_site19_office))]
site19_cowANDfarmerANDoffice_ARGs <- site19_cowANDfarmerANDoffice_ARGs[which(site19_cowANDfarmerANDoffice_ARGs %in% colnames(shortbred_site19_farmer))]

###Looking at the shortBRED calls for the Prevotella MAGs
Prevotella_shortbred <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_shortbred_Prevotella/combined_shortBRED_Prevotella_edited.txt", header = TRUE)
Prevotella_shortbred$Count <- 1
colnames(Prevotella_shortbred)[2] <- "SID"

MAG_metadata <- Prevotella_shortbred[,c(1,2)]
MAG_metadata <- MAG_metadata %>% group_by(MAG) %>% summarise_all(unique)
MAG_metadata <- merge(MAG_metadata, samples_metadata[,c(1,13:15)] %>% group_by(SID) %>% summarise_all(unique), by.x = "SID", by.y = "SID")

Prevotella_shortbred_ARGs <- Prevotella_shortbred[,c(1,4,5)]
Prevotella_shortbred_ARGs <- dcast(Prevotella_shortbred_ARGs, MAG ~ ARG, value.var = "Count", fun.aggregate = sum)
rownames(Prevotella_shortbred_ARGs) <- Prevotella_shortbred_ARGs$MAG
Prevotella_shortbred_ARGs <- Prevotella_shortbred_ARGs[,2:34]
ARG_Prevotella_Jaccard <- vegdist(Prevotella_shortbred_ARGs, method = "jaccard")
pcoa_Prevotella_jacc <- pco(ARG_Prevotella_Jaccard, k = 4)
pcoa_Prevotella_jacc <- as.data.frame(pcoa_Prevotella_jacc$points)
pcoa_Prevotella_jacc <- merge(pcoa_Prevotella_jacc, MAG_metadata, by.x = "row.names", by.y = "MAG")
colnames(pcoa_Prevotella_jacc)[1] <- "MAG"

ggplot(data = pcoa_Prevotella_jacc, aes(x = V1, y = V2, color = Group)) + geom_point()

pairwise_ARG_Prevotella_Jaccard <- as.matrix(ARG_Prevotella_Jaccard)
pairwise_ARG_Prevotella_Jaccard <- melt(pairwise_ARG_Prevotella_Jaccard)
pairwise_ARG_Prevotella_Jaccard <- pairwise_ARG_Prevotella_Jaccard[!duplicated(apply(pairwise_ARG_Prevotella_Jaccard, 1, function(x) paste(sort(x), collapse = ""))),]
pairwise_ARG_Prevotella_Jaccard <- pairwise_ARG_Prevotella_Jaccard[which(pairwise_ARG_Prevotella_Jaccard$Var1 != pairwise_ARG_Prevotella_Jaccard$Var2),]
colnames(pairwise_ARG_Prevotella_Jaccard)[1] <- "MAG1"
colnames(pairwise_ARG_Prevotella_Jaccard)[2] <- "MAG2"
pairwise_ARG_Prevotella_Jaccard <- merge(pairwise_ARG_Prevotella_Jaccard, MAG_metadata, by.x = "MAG1", by.y = "MAG")
colnames(pairwise_ARG_Prevotella_Jaccard)[4] <- "SID1"
colnames(pairwise_ARG_Prevotella_Jaccard)[5] <- "Group1"
colnames(pairwise_ARG_Prevotella_Jaccard)[6] <- "Subject1"
colnames(pairwise_ARG_Prevotella_Jaccard)[7] <- "Site1"
pairwise_ARG_Prevotella_Jaccard <- merge(pairwise_ARG_Prevotella_Jaccard, MAG_metadata, by.x = "MAG2", by.y = "MAG")
colnames(pairwise_ARG_Prevotella_Jaccard)[8] <- "SID2"
colnames(pairwise_ARG_Prevotella_Jaccard)[9] <- "Group2"
colnames(pairwise_ARG_Prevotella_Jaccard)[10] <- "Subject2"
colnames(pairwise_ARG_Prevotella_Jaccard)[11] <- "Site2"
pairwise_ARG_Prevotella_Jaccard <- pairwise_ARG_Prevotella_Jaccard[which(pairwise_ARG_Prevotella_Jaccard$Subject1 != pairwise_ARG_Prevotella_Jaccard$Subject2),]
pairwise_ARG_Prevotella_Jaccard$GroupPair <- ""

for (i in seq(1, length(pairwise_ARG_Prevotella_Jaccard$GroupPair), 1)) {
  pairwise_ARG_Prevotella_Jaccard$GroupPair[i] <- paste(sort(c(pairwise_ARG_Prevotella_Jaccard$Group1[i], pairwise_ARG_Prevotella_Jaccard$Group2[i]))[1], sort(c(pairwise_ARG_Prevotella_Jaccard$Group1[i], pairwise_ARG_Prevotella_Jaccard$Group2[i]))[2])
}

ggplot(pairwise_ARG_Prevotella_Jaccard[which(pairwise_ARG_Prevotella_Jaccard$Site1 == pairwise_ARG_Prevotella_Jaccard$Site2),], aes(x = GroupPair, y = value)) + geom_violin() + geom_jitter()
pairwise_ARG_Prevotella_Jaccard[which(pairwise_ARG_Prevotella_Jaccard$value != 1),c(12,3)] %>% group_by(GroupPair) %>% summarise(value = mean(value))

Prevotella_shortbred_abxClasses <- Prevotella_shortbred[,c(1,3,5)]
Prevotella_shortbred_abxClasses <- dcast(Prevotella_shortbred_abxClasses, MAG ~ AbxClass, value.var = "Count", fun.aggregate = sum)
rownames(Prevotella_shortbred_abxClasses) <- Prevotella_shortbred_abxClasses$MAG
Prevotella_shortbred_abxClasses <- Prevotella_shortbred_abxClasses[,2:11]
AbxClass_Prevotella_Jaccard <- vegdist(Prevotella_shortbred_abxClasses, method = "jaccard")
pcoa_AbxClass_Prev_jaccard <- pco(AbxClass_Prevotella_Jaccard, k = 4)
pcoa_AbxClass_Prev_jaccard <- as.data.frame(pcoa_AbxClass_Prev_jaccard$points)
pcoa_AbxClass_Prev_jaccard <- merge(pcoa_AbxClass_Prev_jaccard, MAG_metadata, by.x = "row.names", by.y = "MAG")
colnames(pcoa_AbxClass_Prev_jaccard)[1] <- "MAG"

ggplot(data = pcoa_AbxClass_Prev_jaccard, aes(x = V1, y = V2, color = Group)) + geom_point()

###Using Lindsey's script, looked at the genes that are correlated with the shared genera (strain sharing) found to be enriched in farmer noses relative to office workers
###Here, plotting a heatmap of those genera and those ARGs
cow_sharedGenera_corARG <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/Lindsey's analysis/BManalysis/230114_Cow_Genus-ARG_correlation_SharedGenesGenera_v1.csv")
cow_sharedGenera_corARG_03 <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/Lindsey's analysis/BManalysis/230308_Cow_Genus-ARG_correlation_SharedGenera_0.3_v1.csv")

#The ARG names are off, so laoding in this dataframe to convert the ARG names back to their original names
cow_sharedGenera_ARGs <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/Lindsey's analysis/BManalysis/221221_Cow_Genus-ARG_correlation_SharedGenera_v3.csv")

cow_sharedGenera_corARG <- merge(cow_sharedGenera_corARG, cow_sharedGenera_ARGs, by.x = "gene", by.y = "M2name")
cow_sharedGenera_corARG <- merge(cow_sharedGenera_corARG, ARG_mapping, by.x = "Ogname", by.y = "Gene")

cow_sharedGenera_corARG$CowEnrich <- ""

for (i in seq(1, length(cow_sharedGenera_corARG$CowEnrich), 1)) {
  if (cow_sharedGenera_corARG$Ogname[i] %in% cowARGs$OGname) {
    cow_sharedGenera_corARG$CowEnrich[i] <- "Yes"
  } else {
    cow_sharedGenera_corARG$CowEnrich[i] <- "No"
  }
}

cow_sharedGenera_corARG$genus <- factor(cow_sharedGenera_corARG$genus, levels = c("Prevotellaceae.UCG.003",
                                                                                  "Romboutsia",
                                                                                  "Turicibacter",
                                                                                  "Rikenellaceae.RC9.gut.group",
                                                                                  "Prevotellaceae.UCG.004",
                                                                                  "Treponema", 
                                                                                  "Bifidobacterium"))

cow_sharedGenera_corARG$Ogname <- factor(cow_sharedGenera_corARG$Ogname, levels = c(
  "DOME_05090|P5-TE_98|gene-1|+|188-688|complete",
  "DOME_04541|P7-TE_16|gene-1|-|111-2114|complete",
  "DLFM12_07041|DFP09SR2_TMP_354|gene-2|+|176-970|complete",
  "DLFM12_04398|DFP06SR1_TMP_30|gene-2|+|457-942|complete",
  "DLFM04_03505|AmoxConc-TwinA-Time1_TMP_12|gene-2|+|1014-1502|complete",
  "494496960|WP_007286427_1|1|1|erm(Q)|erm(Q)|target_modification|2|MACROLIDE|MACROLIDE|23S_rRNA_(adenine(2058)-N(6))-methyltransferase_Erm(Q)",
  "489575244|WP_003479690_1|1|1|tetA(P)|tetA(P)|efflux|2|TETRACYCLINE|TETRACYCLINE|tetracycline_efflux_MFS_transporter_TetA(P)",
  "DOME_00805|P6-CH_34|gene-1|-|55-405|complete",
  "DOME_05695|P5-CH_57|gene-2|-|438-899|complete",
  "DLFM12_02266|DFP03SR1_TMP_127|gene-2|+|637-1134|complete",
  "DOME_06223|P8-TR-SX_134|gene-2|+|1308-1826|complete",
  "DLFM12_00762|DFP02SR1_DCS_331|gene-2|-|526-1653|complete",
  "DOME_01131|P1-TE_215|gene-1|+|119-1828|complete",
  "DOME_06687|P6-CX_137|gene-2|+|304-984|complete",
  "DOME_00554|P7-TR-SX_131|gene-2|-|666-1808|complete",
  "DOME_04760|P5-PI_307|gene-3|+|2439-2822|complete",
  "DOME_06197|P8-TR-SX_94|gene-2|+|505-1182|complete",
  "DOME_02841|P4-GE_66|gene-2|-|631-1422|complete",
  "DLFM02_01841|F15_DCS_39|gene-5|-|2911-3189|complete",
  "DOME_06516|P8-PE_61|gene-1|-|27-1598|complete",
  "DLFM11_02621|P21SR5_CTX_14|gene-2|+|518-1321|complete",
  "DLFM07_02181|01G-000_PEN_1|gene-1|+|134-676|complete",
  "DOME_00905|P5-CX_24|gene-3|-|2118-3083|complete",
  "DOME_05912|P3-CH_16|gene-1|+|85-651|complete",
  "DOME_06203|P8-TR-SX_103|gene-2|-|244-735|complete",
  "DOME_01056|P1-TE_94|gene-1|-|21-713|complete",
  "DOME_00269|P8-TE_49|gene-4|+|1722-3059|complete",
  "DOME_00905|P5-CX_24|gene-4|-|3186-3680|complete",
  "DOME_00905|P5-CX_24|gene-1|+|532-1143|complete"
))

for (i in c(
  "DOME_05090|P5-TE_98|gene-1|+|188-688|complete",
  "DOME_04541|P7-TE_16|gene-1|-|111-2114|complete",
  "DLFM12_07041|DFP09SR2_TMP_354|gene-2|+|176-970|complete",
  "DLFM12_04398|DFP06SR1_TMP_30|gene-2|+|457-942|complete",
  "DLFM04_03505|AmoxConc-TwinA-Time1_TMP_12|gene-2|+|1014-1502|complete",
  "494496960|WP_007286427_1|1|1|erm(Q)|erm(Q)|target_modification|2|MACROLIDE|MACROLIDE|23S_rRNA_(adenine(2058)-N(6))-methyltransferase_Erm(Q)",
  "489575244|WP_003479690_1|1|1|tetA(P)|tetA(P)|efflux|2|TETRACYCLINE|TETRACYCLINE|tetracycline_efflux_MFS_transporter_TetA(P)",
  "DOME_00805|P6-CH_34|gene-1|-|55-405|complete",
  "DOME_05695|P5-CH_57|gene-2|-|438-899|complete",
  "DLFM12_02266|DFP03SR1_TMP_127|gene-2|+|637-1134|complete",
  "DOME_06223|P8-TR-SX_134|gene-2|+|1308-1826|complete",
  "DLFM12_00762|DFP02SR1_DCS_331|gene-2|-|526-1653|complete",
  "DOME_01131|P1-TE_215|gene-1|+|119-1828|complete",
  "DOME_06687|P6-CX_137|gene-2|+|304-984|complete",
  "DOME_00554|P7-TR-SX_131|gene-2|-|666-1808|complete",
  "DOME_04760|P5-PI_307|gene-3|+|2439-2822|complete",
  "DOME_06197|P8-TR-SX_94|gene-2|+|505-1182|complete",
  "DOME_02841|P4-GE_66|gene-2|-|631-1422|complete",
  "DLFM02_01841|F15_DCS_39|gene-5|-|2911-3189|complete",
  "DOME_06516|P8-PE_61|gene-1|-|27-1598|complete",
  "DLFM11_02621|P21SR5_CTX_14|gene-2|+|518-1321|complete",
  "DLFM07_02181|01G-000_PEN_1|gene-1|+|134-676|complete",
  "DOME_00905|P5-CX_24|gene-3|-|2118-3083|complete",
  "DOME_05912|P3-CH_16|gene-1|+|85-651|complete",
  "DOME_06203|P8-TR-SX_103|gene-2|-|244-735|complete",
  "DOME_01056|P1-TE_94|gene-1|-|21-713|complete",
  "DOME_00269|P8-TE_49|gene-4|+|1722-3059|complete",
  "DOME_00905|P5-CX_24|gene-4|-|3186-3680|complete",
  "DOME_00905|P5-CX_24|gene-1|+|532-1143|complete"
)) {
  print(paste(i, ":", ARG_mapping[which(ARG_mapping$Gene == i),2]))
}

cow_sharedGenera_corARG_filtered <- cow_sharedGenera_corARG[which(cow_sharedGenera_corARG$BH_corrected_pvalue < 0.05 & cow_sharedGenera_corARG$spearman_correlation_value > 0.5),]

cow_sharedGenera_corARG_filtered <- dcast(cow_sharedGenera_corARG_filtered[,c(1,3,6)], genus ~ Ogname, value.var = "spearman_correlation_value")
rownames(cow_sharedGenera_corARG_filtered) <- cow_sharedGenera_corARG_filtered$genus
cow_sharedGenera_corARG_filtered <- cow_sharedGenera_corARG_filtered[,c(2:30)]
#cow_sharedGenera_corARG_filtered[is.na(cow_sharedGenera_corARG_filtered)] <- 0
#test <- sapply(cow_sharedGenera_corARG_filtered, as.numeric)
#rownames(test) <- rownames(cow_sharedGenera_corARG_filtered)
#cow_sharedGenera_corARG_filtered <- as.data.frame(test)
cow_sharedGenera_corARG_filtered <- cow_sharedGenera_corARG_filtered[order(nrow(cow_sharedGenera_corARG_filtered):1),]
cow_sharedGenera_corARG_filtered_annotation <- as.data.frame(colnames(cow_sharedGenera_corARG_filtered))
colnames(cow_sharedGenera_corARG_filtered_annotation)[1] <- "ARG"
cow_sharedGenera_corARG_filtered_annotation <- merge(cow_sharedGenera_corARG_filtered_annotation, ARG_mapping, by.x = "ARG", by.y = "Gene")
cow_sharedGenera_corARG_filtered_annotation$Enriched <- ""
for (i in seq(1,nrow(cow_sharedGenera_corARG_filtered_annotation))) {
  if (cow_sharedGenera_corARG_filtered_annotation$ARG[i] %in% cowARGs$OGname)  {
    cow_sharedGenera_corARG_filtered_annotation$Enriched[i] <- "Yes"
  } else {
    cow_sharedGenera_corARG_filtered_annotation$Enriched[i] <- "No"
  }
}
rownames(cow_sharedGenera_corARG_filtered_annotation) <- cow_sharedGenera_corARG_filtered_annotation$ARG
cow_sharedGenera_corARG_filtered_annotation <- cow_sharedGenera_corARG_filtered_annotation[,c(3,2)]

my_colors = list(
  AbxClass = c("TETRACYCLINE" = "#8C2B0E", "MACROLIDE" = "#233F6C", "TRIMETHOPRIM" = "#B25422", "AMINOGLYCOSIDE" = "#81565F", "BETA-LACTAM" = "#FEB359", "GLYCOPEPTIDE" = "#D8813B",
               "QUINOLONE" = "#5B4C64", "FOSFOMYCIN" = "#B47E83", "LINCOSAMIDE" = "#9F7E59",
               "OTHER" = "#848484", "EFFLUX" = "#848484", "MULTIDRUG" = "#848484", "STREPTOTHRICIN" = "#848484", "QUATERNARY AMMONIUM" = "#848484"),
  Enriched = c("Yes" = "#151514", "No" = "#f9f9f9")
)

p <- pheatmap(as.matrix(cow_sharedGenera_corARG_filtered), na_col = "white", cluster_rows = FALSE, cluster_cols = FALSE, 
         show_colnames = FALSE, color = rev(natparks.pals("Arches2", type = c("continuous"))),
         annotation_col = cow_sharedGenera_corARG_filtered_annotation, annotation_colors = my_colors)
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/SpearCorr_SharedGenera_CowARGs_v1.svg", p, bg = "transparent", width = 15, height = 20)

cow_sharedGenera_corARG_filtered <- cow_sharedGenera_corARG[which(cow_sharedGenera_corARG$BH_corrected_pvalue < 0.05 & cow_sharedGenera_corARG$spearman_correlation_value > 0.5),]

ggplot() + geom_tile(data = cow_sharedGenera_corARG[which(cow_sharedGenera_corARG$genus != "Prevotellaceae.UCG.003"),], aes(x = Ogname, y = genus), fill = "white", color = "black") + 
  geom_tile(data = cow_sharedGenera_corARG_filtered, aes(x = Ogname, y = genus, fill = spearman_correlation_value), color = "black") + theme_test() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_gradientn(colors = rev(natparks.pals("Arches2")))
ggplot() + geom_tile(data = cow_sharedGenera_corARG[which(cow_sharedGenera_corARG$genus != "Prevotellaceae.UCG.003"),], aes(x = Ogname, y = genus), fill = "white", color = "black") + 
  geom_tile(data = cow_sharedGenera_corARG_filtered, aes(x = Ogname, y = genus, fill = spearman_correlation_value), color = "black") + theme_test() + 
  theme(axis.text.x = element_blank(), panel.border = element_rect(size = 0)) + scale_fill_gradientn(colors = rev(natparks.pals("Arches2")))
p <- ggplot() + geom_tile(data = cow_sharedGenera_corARG[which(cow_sharedGenera_corARG$genus != "Prevotellaceae.UCG.003"),], aes(x = Ogname, y = genus), fill = "white", color = "black") + 
  geom_tile(data = cow_sharedGenera_corARG_filtered, aes(x = Ogname, y = genus, fill = spearman_correlation_value), color = "black") + theme_test() + 
  theme(axis.text = element_text(angle = 90, vjust = 0.5, hjust=1), panel.border = element_rect(size = 0), axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "none") + 
  scale_fill_gradientn(colors = rev(natparks.pals("Arches2")))
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/SpearCorr_SharedGenera_CowARGs_v4.svg", p, bg = "transparent", width = 20, height = 30)


cow_allSharedGenera_corARG <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/Lindsey's analysis/BManalysis/230405_Cow_Genus-ARG_correlation_AllSharedGenera_v1.csv",
                                       header = TRUE)

cow_allSharedGenera_corARG_genes <- as.data.frame(unique(cow_allSharedGenera_corARG$gene))
colnames(cow_allSharedGenera_corARG_genes)[1] <- "Gene"
cow_allSharedGenera_corARG_genes <- merge(cow_allSharedGenera_corARG_genes, ARG_mapping, by = "Gene")
write_csv(cow_allSharedGenera_corARG_genes, file = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/Lindsey's analysis/BManalysis/230424_Cow_Genus-ARG_correlation_AllSharedGenera_Genes_v1.csv")

cow_allSharedGenera_corARG_heatmap <- dcast(cow_allSharedGenera_corARG[,c(1,2,5)], genus ~ gene, value.var = "spearman_correlation_value")
rownames(cow_allSharedGenera_corARG_heatmap) <- cow_allSharedGenera_corARG_heatmap$genus
cow_allSharedGenera_corARG_heatmap <- cow_allSharedGenera_corARG_heatmap[,c(2:191)]
cow_allSharedGenera_corARG_heatmap[is.na(cow_allSharedGenera_corARG_heatmap)] <- 0

pheatmap(as.matrix(cow_allSharedGenera_corARG_heatmap), na_col = "white", cluster_rows = FALSE, cluster_cols = FALSE, 
         show_colnames = FALSE, color = rev(natparks.pals("Arches2", type = c("continuous"))))
pheatmap(as.matrix(cow_allSharedGenera_corARG_heatmap), na_col = "white", 
         show_colnames = TRUE, color = rev(natparks.pals("Arches2", type = c("continuous"))))

test <- hclust(dist(t(cow_allSharedGenera_corARG_heatmap), method = "euclidian"), method = "complete")
#write_delim(as.data.frame(test$order), file = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/Lindsey's analysis/BManalysis/230406_Cow_Genus-ARG_correlation_heatmap_GeneOrder.txt")
#write_delim(as.data.frame(test$labels), file = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/Lindsey's analysis/BManalysis/230406_Cow_Genus-ARG_correlation_heatmap_Genes.txt")

cow_allSharedGenera_corARG_heatmap <- dcast(cow_allSharedGenera_corARG[,c(1,2,5)], genus ~ gene, value.var = "spearman_correlation_value")
rownames(cow_allSharedGenera_corARG_heatmap) <- cow_allSharedGenera_corARG_heatmap$genus
cow_allSharedGenera_corARG_heatmap <- cow_allSharedGenera_corARG_heatmap[,c(2:191)]
cow_allSharedGenera_corARG_heatmap <- cow_allSharedGenera_corARG_heatmap[c("Holdemanella", "Bifidobacterium", "Romboutsia", "Turicibacter", "Treponema", "Rikenellaceae.RC9.gut.group", "Prevotellaceae.UCG.004", "Agathobacter", "Blautia"),
                                   c(
                                     "DOME_00905|P5-CX_24|gene-1|+|532-1143|complete",
                                     "DOME_00905|P5-CX_24|gene-4|-|3186-3680|complete",
                                     "DOME_00269|P8-TE_49|gene-4|+|1722-3059|complete",
                                     "DOME_01056|P1-TE_94|gene-1|-|21-713|complete",
                                     "DOME_06203|P8-TR-SX_103|gene-2|-|244-735|complete",
                                     "DOME_00905|P5-CX_24|gene-3|-|2118-3083|complete",
                                     "DOME_05912|P3-CH_16|gene-1|+|85-651|complete",
                                     "DLFM07_02181|01G-000_PEN_1|gene-1|+|134-676|complete",
                                     "DLFM02_01841|F15_DCS_39|gene-5|-|2911-3189|complete",
                                     "DOME_02841|P4-GE_66|gene-2|-|631-1422|complete",
                                     "DOME_00554|P7-TR-SX_131|gene-2|-|666-1808|complete",
                                     "DOME_06687|P6-CX_137|gene-2|+|304-984|complete",
                                     "DOME_04760|P5-PI_307|gene-3|+|2439-2822|complete",
                                     "DOME_06197|P8-TR-SX_94|gene-2|+|505-1182|complete",
                                     "DLFM11_02621|P21SR5_CTX_14|gene-2|+|518-1321|complete",
                                     "DOME_06516|P8-PE_61|gene-1|-|27-1598|complete",
                                     "DLFM12_03383|DFP05SR1_TET_65|gene-2|-|360-854|complete",
                                     "DLFM12_02359|DFP04SR1_DCS_26|gene-2|-|496-1488|complete",
                                     "DOME_00742|P9-TE_114|gene-2|+|2371-2964|complete",
                                     "DOME_01214|P2-TE_38|gene-4|+|3881-4573|complete",
                                     "DLFM02_01543|F16_DCS_128|gene-3|-|740-1834|complete",
                                     "DLFM11_03135|P08SR2_DCS_170|gene-2|-|498-1658|complete",
                                     "DLFM11_00523|P05SR4_AZM_56|gene-2|-|99-863|complete",
                                     "DLFM02_03451|F33_TMP_94|gene-1|+|96-1121|complete",
                                     "DLFM02_03356|F33_SXT_52|gene-2|-|827-1675|complete",
                                     "DLFM02_00147|F30_DCS_95|gene-1|+|366-1925|complete",
                                     "DLFM02_01416|F18_DCS_41|gene-2|-|165-1340|complete",
                                     "DOME_03008|P12-CH_24|gene-2|+|434-886|complete",
                                     "DOME_02235|P9-CX_96|gene-3|-|765-1142|complete",
                                     "DLFM02_00865|F22_SXT_16|gene-1|+|17-1042|complete",
                                     "DLFM11_03897|P05SR2_FOX_103|gene-2|-|743-1180|complete",
                                     "DLFM11_05079|P17SR3_PEN_200|gene-2|+|419-1048|complete",
                                     "DLFM11_03053|P05SR2_DCS_161|gene-2|-|91-735|complete",
                                     "DLFM02_01128|F19_SXT_59|gene-1|+|21-815|complete",
                                     "DLFM11_03038|P05SR2_DCS_124|gene-2|+|78-1256|complete",
                                     "DLFM02_02671|F06_GEN_17|gene-2|+|1026-2033|complete",
                                     "DLFM11_03651|P03SR2_FOX_207|gene-2|-|359-841|complete",
                                     "DLFM02_02639|F06_SXT_86|gene-3|-|2098-2946|complete",
                                     "gb|BAE78082_1|ARO_3003550|mdtP__Escherichia_coli",
                                     "DOME_05581|P11-TR-SX_5|gene-1|-|105-824|complete",
                                     "DOME_04074|P2-CH_47|gene-3|+|2317-2937|complete",
                                     "DOME_01915|P13-CX_67|gene-2|+|1035-1610|complete",
                                     "DOME_01724|P9-TR-SX_404|gene-2|-|354-1118|complete",
                                     "DOME_00382|P10-TE_151|gene-6|-|5648-6130|complete",
                                     "DOME_00303|P10-TE_35|gene-2|-|164-1051|complete",
                                     "DLFM13_01139|CBP02SR2_SXT_38|gene-1|+|465-1634|complete",
                                     "DLFM13_00362|CBP02SR2_CAZ_50|gene-3|-|1391-2305|complete",
                                     "DLFM13_00268|CBP01SR2_TZP_10|gene-2|-|151-1113|complete",
                                     "DLFM12_05667|DFP08SR2_TET_78|gene-2|-|224-706|complete",
                                     "DLFM11_07448|P02SR2_TMP_233|gene-4|+|1928-2506|complete",
                                     "DLFM11_07392|P02SR2_TMP_45|gene-3|+|1302-1790|complete",
                                     "DLFM11_05602|P05SR2_SXT_22|gene-1|+|162-1010|complete",
                                     "DLFM11_02252|P15SR5_CHL_38|gene-2|-|242-832|complete",
                                     "DLFM11_01230|P17SR4_AZM_142|gene-5|-|2103-2762|complete",
                                     "DLFM11_00136|P12SR3_ATM_121|gene-2|+|404-1585|complete",
                                     "DLFM11_00001|P01SR2_ATM_6|gene-2|-|583-1317|complete",
                                     "DLFM07_01138|03F-000_CHL_13|gene-4|-|2369-2821|complete",
                                     "DLFM04_03710|AmoxConc-Mom_CAZ_47|gene-5|-|3622-4332|complete",
                                     "DLFM04_03632|AmoxConc-Mom_TMP_89|gene-2|+|852-1448|complete",
                                     "DLFM04_00330|Control-TwinB-Time2_CST_6|gene-2|+|1361-1957|complete",
                                     "DLFM04_00122|Control-TwinB-Time3_FEP_8|gene-2|+|194-1138|complete",
                                     "DLFM03_03175|S08_FOX_1|gene-4|-|3048-3995|complete",
                                     "DLFM02_03488|F33_TMP_88|gene-2|-|232-723|complete",
                                     "DLFM02_03399|F33_SXT_5|gene-2|+|1124-1627|complete",
                                     "DLFM02_02956|F05_DCS_47|gene-2|+|580-1641|complete",
                                     "DLFM02_02926|F05_DCS_86|gene-2|+|524-1852|complete",
                                     "DLFM02_02643|F06_DCS_145|gene-2|-|886-1986|complete",
                                     "DLFM02_02426|F07_SXT_30|gene-3|+|1507-1776|complete",
                                     "DLFM02_02287|F09_SXT_1|gene-2|+|362-1210|complete",
                                     "DLFM02_01901|F15_TMP_20|gene-2|+|243-734|complete",
                                     "DLFM02_01268|F19_DCS_106|gene-2|+|907-1641|complete",
                                     "DLFM02_01029|F21_SXT_3|gene-3|+|1841-2335|complete",
                                     "DLFM02_00785|F22_DCS_127|gene-2|-|367-960|complete",
                                     "DLFM02_00720|F23_DCS_58|gene-2|+|320-1276|complete",
                                     "488251003|WP_002322211_1|1|1|dfrF|dfrF|insensitivity|2|TRIMETHOPRIM|TRIMETHOPRIM|trimethoprim-resistant_dihydrofolate_reductase_DfrF",
                                     "DLFM02_00337|F26_DCS_55|gene-1|+|55-1104|complete",
                                     "DOME_03508|P7-TR_239|gene-1|+|81-917|complete",
                                     "DLFM10_03704|F19_TET_14|gene-3|-|1943-3046|complete",
                                     "164457634|BAF96541_1|1|1|nleA|nleA||1|||type_III_secretion_system_effector_NleA",
                                     "25987939|AAN76093_1|1|1|iroN|iroN||1|||siderophore_salmochelin_receptor_IroN",
                                     "DLFM07_00940|03G-000_TET_22|gene-1|+|213-1565|complete",
                                     "DOME_06029|P3-CH_282|gene-2|-|1322-1885|complete",
                                     "488223297|WP_002294505_1|1|1|ant(6)-Ia|ant(6)-Ia|nucleotidyltransferase|2|STREPTOMYCIN|AMINOGLYCOSIDE|aminoglycoside_nucleotidyltransferase_ANT(6)-Ia",
                                     "DOME_06849|P9-CP_6|gene-6|+|3675-5030|complete",
                                     "DOME_06516|P8-PE_61|gene-3|+|2038-3480|complete",
                                     "DOME_06450|P8-TR-SX_525|gene-3|+|1286-1771|complete",
                                     "DOME_05339|P9-PI_24|gene-1|-|71-973|complete",
                                     "DOME_03636|P13-PE_23|gene-2|-|582-1838|complete",
                                     "DOME_03218|P12-TR_310|gene-1|+|58-537|complete",
                                     "DOME_03133|P12-TR_168|gene-1|+|112-906|complete",
                                     "DOME_03132|P12-TR_167|gene-2|+|174-656|complete",
                                     "DOME_02993|P12-CH_6|gene-1|+|34-1200|complete",
                                     "DOME_02791|P12-TR-SX_351|gene-1|+|400-1821|complete",
                                     "DOME_02744|P12-TR-SX_269|gene-1|+|17-1558|complete",
                                     "DOME_02666|P12-TR-SX_129|gene-2|-|527-1027|complete",
                                     "DOME_02646|P12-TR-SX_94|gene-5|+|1948-2799|complete",
                                     "DOME_02417|P10-CX_32|gene-3|-|2399-3226|complete",
                                     "DOME_02405|P10-CX_15|gene-2|+|596-988|complete",
                                     "DOME_02203|P9-CX_58|gene-3|+|1285-1725|complete",
                                     "DOME_02167|P9-CX_12|gene-3|+|731-2329|complete",
                                     "DOME_01689|P9-TR-SX_339|gene-4|+|1316-2056|complete",
                                     "DOME_01689|P9-TR-SX_339|gene-2|+|196-1047|complete",
                                     "DOME_00712|P9-TE_65|gene-8|-|7710-8483|complete",
                                     "DOME_00710|P9-TE_63|gene-5|+|3139-4305|complete",
                                     "DOME_00700|P9-TE_48|gene-3|+|1211-2632|complete",
                                     "DOME_00305|P10-TE_40|gene-3|+|1398-2285|complete",
                                     "DLFM13_04112|CBP14SR2_CAZ_6|gene-2|-|1426-2529|complete",
                                     "DLFM13_03336|CBP10SR2_TET_52|gene-2|-|618-1325|complete",
                                     "DLFM13_02729|CBP08SR2_CAZ_30|gene-2|-|841-1386|complete",
                                     "DLFM12_03562|DFP05SR1_TET_413|gene-2|+|468-959|complete",
                                     "DLFM11_08277|P14SR3_TZP_142|gene-2|-|798-1646|complete",
                                     "DLFM11_07937|P12SR3_TZP_4|gene-2|-|55-831|complete",
                                     "DLFM11_06989|P15SR4_TIC_240|gene-2|-|994-1437|complete",
                                     "DLFM11_06924|P15SR4_TIC_77|gene-2|+|492-1238|complete",
                                     "DLFM11_05694|P12SR3_SXT_35|gene-2|+|776-1840|complete",
                                     "DLFM11_05533|P02SR5_SXT_55|gene-3|-|1193-2041|complete",
                                     "DLFM11_04959|P15SR3_PEN_230|gene-1|-|34-525|complete",
                                     "DLFM11_04403|P22SR1_FOX_18|gene-2|+|375-980|complete",
                                     "DLFM11_04079|P08SR2_FOX_285|gene-2|-|1013-1696|complete",
                                     "DLFM11_04012|P08SR2_FOX_108|gene-2|-|181-600|complete",
                                     "DLFM11_03292|P22SR1_DCS_89|gene-2|+|292-1464|complete",
                                     "DLFM11_03147|P08SR2_DCS_199|gene-2|+|285-1508|complete",
                                     "DLFM11_02855|P02SR2_DCS_302|gene-2|-|261-797|complete",
                                     "DLFM11_02450|P11SR3_CST_10|gene-2|-|233-370|complete",
                                     "DLFM11_02447|P11SR3_CST_6|gene-2|-|376-1551|complete",
                                     "DLFM11_02362|P23SR5_CHL_22|gene-2|+|352-912|complete",
                                     "DLFM11_02135|P10SR3_CHL_51|gene-2|+|367-1542|complete",
                                     "DLFM11_01874|P17SR3_CAZ_60|gene-1|-|85-567|complete",
                                     "DLFM11_01418|P19SR4_AZM_313|gene-1|-|126-890|complete",
                                     "DLFM11_00738|P11SR4_AZM_95|gene-2|+|505-1395|complete",
                                     "DLFM11_00138|P12SR3_ATM_125|gene-3|+|721-1419|complete",
                                     "DLFM08_00520|B4_TMP_12|gene-2|+|1243-2073|complete",
                                     "DLFM07_01560|03C-000_CHL_8|gene-3|-|1394-3187|complete",
                                     "DLFM07_01364|03D-000_PEN_88|gene-2|-|724-1593|complete",
                                     "DLFM06_01668|1-C6_FOX_19|gene-3|-|1239-1748|complete",
                                     "DLFM05_00050|F23_TET_1|gene-1|-|69-1415|complete",
                                     "DLFM04_03148|AmoxDisc-Mom_DCS_160|gene-2|+|319-1356|complete",
                                     "DLFM04_03113|AmoxDisc-Mom_DCS_231|gene-2|+|1443-2498|complete",
                                     "DLFM04_03101|AmoxDisc-Mom_DCS_251|gene-1|+|33-1085|complete",
                                     "DLFM04_02975|AmoxDisc-Mom_SXT_22|gene-6|-|1972-2451|complete",
                                     "DLFM02_03572|F33_DCS_55|gene-2|+|226-1263|complete",
                                     "DLFM02_03498|F33_TMP_23|gene-2|+|279-770|complete",
                                     "DLFM02_03297|F04_DCS_179|gene-1|+|432-1499|complete",
                                     "DLFM02_03231|F04_TMP_87|gene-2|-|287-781|complete",
                                     "DLFM02_03120|F04_DCS_57|gene-3|+|1346-2530|complete",
                                     "DLFM02_03101|F04_DCS_173|gene-3|+|843-1916|complete",
                                     "DLFM02_02844|F06_SXT_49|gene-3|-|1949-2797|complete",
                                     "DLFM02_02812|F06_TMP_86|gene-2|+|311-805|complete",
                                     "DLFM02_02756|F06_SXT_28|gene-4|-|1974-2822|complete",
                                     "DLFM02_02710|F06_TMP_158|gene-4|+|1628-2341|complete",
                                     "DLFM02_02639|F06_SXT_86|gene-2|-|1526-2047|complete",
                                     "DLFM02_02562|F06_TMP_149|gene-2|+|1180-1974|complete",
                                     "DLFM02_02314|F09_SXT_46|gene-3|+|1119-1607|complete",
                                     "DLFM02_02261|F09_CHL_1|gene-2|-|501-824|complete",
                                     "DLFM02_02250|F09_DCS_3|gene-2|-|222-1112|complete",
                                     "DLFM02_02186|F11_SXT_62|gene-2|+|258-749|complete",
                                     "DLFM02_02171|F11_SXT_38|gene-3|-|1768-2616|complete",
                                     "DLFM02_02104|F12_DCS_38|gene-2|-|328-1380|complete",
                                     "DLFM02_01961|F15_DCS_15|gene-2|+|79-1134|complete",
                                     "DLFM02_01930|F15_DCS_26|gene-2|-|550-1587|complete",
                                     "DLFM02_01920|F15_DCS_32|gene-1|+|187-1236|complete",
                                     "DLFM02_01919|F15_DCS_78|gene-1|+|55-1110|complete",
                                     "DLFM02_01859|F15_DCS_35|gene-1|-|454-1515|complete",
                                     "DLFM02_01802|F16_DCS_25|gene-2|+|465-1517|complete",
                                     "DLFM02_01510|F18_DCS_153|gene-2|-|347-1396|complete",
                                     "DLFM02_01377|F18_DCS_49|gene-2|+|238-1302|complete",
                                     "DLFM02_01346|F19_DCS_210|gene-3|+|1544-1906|complete",
                                     "DLFM02_01303|F19_DCS_130|gene-2|-|412-1407|complete",
                                     "DLFM02_01287|F19_DCS_22|gene-2|+|81-1115|complete",
                                     "DLFM02_01128|F19_SXT_59|gene-2|+|846-1361|complete",
                                     "DLFM02_00910|F21_DCS_87|gene-2|+|417-1487|complete",
                                     "DLFM02_00807|F22_TMP_43|gene-2|-|604-1095|complete",
                                     "DLFM02_00797|F22_DCS_116|gene-2|+|232-1293|complete",
                                     "DLFM02_00686|F24_DCS_2|gene-1|-|4-1185|complete",
                                     "DLFM02_00647|F24_TMP_54|gene-4|+|1435-1929|complete",
                                     "DLFM02_00267|F28_TMP_9|gene-2|-|650-1138|complete",
                                     "DLFM02_00647|F24_TMP_54|gene-3|+|1007-1414|complete",
                                     "DOME_01131|P1-TE_215|gene-1|+|119-1828|complete",
                                     "DLFM12_00762|DFP02SR1_DCS_331|gene-2|-|526-1653|complete",
                                     "DOME_06223|P8-TR-SX_134|gene-2|+|1308-1826|complete",
                                     "DLFM12_02266|DFP03SR1_TMP_127|gene-2|+|637-1134|complete",
                                     "DOME_00805|P6-CH_34|gene-1|-|55-405|complete",
                                     "DOME_05695|P5-CH_57|gene-2|-|438-899|complete",
                                     "DLFM12_04398|DFP06SR1_TMP_30|gene-2|+|457-942|complete",
                                     "DLFM12_07041|DFP09SR2_TMP_354|gene-2|+|176-970|complete",
                                     "DOME_04541|P7-TE_16|gene-1|-|111-2114|complete",
                                     "DOME_05090|P5-TE_98|gene-1|+|188-688|complete",
                                     "489575244|WP_003479690_1|1|1|tetA(P)|tetA(P)|efflux|2|TETRACYCLINE|TETRACYCLINE|tetracycline_efflux_MFS_transporter_TetA(P)",
                                     "494496960|WP_007286427_1|1|1|erm(Q)|erm(Q)|target_modification|2|MACROLIDE|MACROLIDE|23S_rRNA_(adenine(2058)-N(6))-methyltransferase_Erm(Q)",
                                     "DLFM04_03505|AmoxConc-TwinA-Time1_TMP_12|gene-2|+|1014-1502|complete"
                                   )]

cow_allSharedGenera_corARG_heatmap_annotation <- as.data.frame(colnames(cow_allSharedGenera_corARG_heatmap))
colnames(cow_allSharedGenera_corARG_heatmap_annotation)[1] <- "ARG"
cow_allSharedGenera_corARG_heatmap_annotation <- merge(cow_allSharedGenera_corARG_heatmap_annotation, ARG_mapping, by.x = "ARG", by.y = "Gene")
cow_allSharedGenera_corARG_heatmap_annotation$Enriched <- ""
for (i in seq(1,nrow(cow_allSharedGenera_corARG_heatmap_annotation))) {
  if (cow_allSharedGenera_corARG_heatmap_annotation$ARG[i] %in% cowARGs$OGname)  {
    cow_allSharedGenera_corARG_heatmap_annotation$Enriched[i] <- "Yes"
  } else {
    cow_allSharedGenera_corARG_heatmap_annotation$Enriched[i] <- "No"
  }
}
rownames(cow_allSharedGenera_corARG_heatmap_annotation) <- cow_allSharedGenera_corARG_heatmap_annotation$ARG
cow_allSharedGenera_corARG_heatmap_annotation <- cow_allSharedGenera_corARG_heatmap_annotation[,c(3,2)]

my_colors = list(
  AbxClass = c("TETRACYCLINE" = "#8C2B0E", "MACROLIDE" = "#233F6C", "TRIMETHOPRIM" = "#B25422", "AMINOGLYCOSIDE" = "#81565F", "BETA-LACTAM" = "#FEB359", "GLYCOPEPTIDE" = "#D8813B",
               "QUINOLONE" = "#5B4C64", "FOSFOMYCIN" = "#B47E83", "LINCOSAMIDE" = "#9F7E59", "PHENICOL" = "#435F90",
               "OTHER" = "#848484", "EFFLUX" = "#848484", "MULTIDRUG" = "#848484", "STREPTOTHRICIN" = "#848484", "QUATERNARY AMMONIUM" = "#848484","PHENICOL/QUINOLONE" = "#848484",
               "BACITRACIN" = "#848484", "STREPTOGRAMIN" = "#848484", "ARSENIC" = "#848484", "AVILAMYCIN" = "#848484", "LINCOSAMIDE/STREPTOGRAMIN" = "#848484", "BLEOMYCIN" = "#848484", 
               "PHENICOL/OXAZOLIDINONE" = "#848484"),
  Enriched = c("Yes" = "#151514", "No" = "#f9f9f9")
)

pheatmap(as.matrix(cow_allSharedGenera_corARG_heatmap), na_col = "white", cluster_rows = FALSE, cluster_cols = FALSE, 
         show_colnames = FALSE, color = rev(natparks.pals("Arches2", type = c("continuous"))))
pheatmap(as.matrix(t(cow_allSharedGenera_corARG_heatmap)), na_col = "white", cluster_rows = FALSE, cluster_cols = FALSE, 
         show_rownames = FALSE, color = rev(natparks.pals("Arches2", type = c("continuous"))))
pheatmap(as.matrix(cow_allSharedGenera_corARG_heatmap), na_col = "white", cluster_rows = FALSE, cluster_cols = FALSE, 
         show_colnames = FALSE, color = rev(natparks.pals("Arches2", type = c("continuous"))),
         annotation_col = cow_allSharedGenera_corARG_heatmap_annotation, annotation_colors = my_colors)

###Looking at the rel abundance in humans of only ARGs present in cows, ARGs enriched in cows, ARGs associated in cows with shared taxa, or ARGs found shared between cow and farmer MAGs
shortbred_counts_wide_cows <- shortbred_counts_wide[which(row.names(shortbred_counts_wide) %notin% human_metadata$Sample_ID),]
shortbred_counts_wide_cows <- shortbred_counts_wide_cows[,which(colSums(shortbred_counts_wide_cows) > 0)]
shortbred_wide_humans_cowPresARGs <- shortbred_counts_wide_humans[,which(colnames(shortbred_counts_wide_humans) %in% colnames(shortbred_counts_wide_cows))]
human_cowPresARG_BrayCurtis <- vegdist(shortbred_wide_humans_cowPresARGs, index = "bray")
human_cowPresARG_pcoa_BC <- pco(human_cowPresARG_BrayCurtis, k = 4)
human_cowPresARG_pcoa_BC <- as.data.frame(human_cowPresARG_pcoa_BC$points)
human_cowPresARG_pcoa_BC <- merge(human_cowPresARG_pcoa_BC, human_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(human_cowPresARG_pcoa_BC)[1] <- "Sample"
human_cowPresARG_pcoa_BC$SeasonFull <- paste(human_cowPresARG_pcoa_BC$Season, human_cowPresARG_pcoa_BC$Year)

ggplot(data = human_cowPresARG_pcoa_BC[which(human_cowPresARG_pcoa_BC$SeasonFull == "SPRING 2"),], aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull)) + stat_ellipse(aes(color = Group)) + 
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test()

human_totalARG_cowPresARG <- as.data.frame(rowSums(shortbred_wide_humans_cowPresARGs))
human_totalARG_cowPresARG <- merge(human_totalARG_cowPresARG, human_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(human_totalARG_cowPresARG)[1] <- "Sample"
colnames(human_totalARG_cowPresARG)[2] <- "TotalRPKM"
human_totalARG_cowPresARG$SeasonFull <- paste(human_totalARG_cowPresARG$Season, human_totalARG_cowPresARG$Year)
human_totalARG_cowPresARG$SeasonFull <- factor(human_totalARG_cowPresARG$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))

compare_means(TotalRPKM ~ Group, data = human_totalARG_cowPresARG, group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")

ggplot(human_totalARG_cowPresARG, aes(x = Group, y = TotalRPKM)) + geom_boxplot() + geom_point(aes(color = Group)) +
  facet_wrap(~SeasonFull, nrow = 1, labeller = season_labeller) +
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test()

human_richness_cowPresARG <- as.data.frame(specnumber(shortbred_wide_humans_cowPresARGs))
human_richness_cowPresARG <- merge(human_richness_cowPresARG, human_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(human_richness_cowPresARG)[1] <- "Sample"
colnames(human_richness_cowPresARG)[2] <- "Richness"
human_richness_cowPresARG$SeasonFull <- paste(human_richness_cowPresARG$Season, human_richness_cowPresARG$Year)
human_richness_cowPresARG$SeasonFull <- factor(human_richness_cowPresARG$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))

compare_means(Richness ~ Group, data = human_richness_cowPresARG, group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")

ggplot(human_richness_cowPresARG, aes(x = Group, y = Richness)) + geom_boxplot() + geom_point(aes(color = Group)) +
  facet_wrap(~SeasonFull, nrow = 1, labeller = season_labeller) +
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test()

#----------------------------------------------------------------------------------------------------------------------------------

shortbred_wide_humans_cowARGs <- shortbred_counts_wide_humans[,which(colnames(shortbred_counts_wide_humans) %in% cowARGs$OGname)]
human_cowARGs_BrayCurtis <- vegdist(shortbred_wide_humans_cowARGs, index = "bray")
human_cowARGs_pcoa_BC <- pco(human_cowARGs_BrayCurtis, k = 4)
human_cowARGs_pcoa_BC <- as.data.frame(human_cowARGs_pcoa_BC$points)
human_cowARGs_pcoa_BC <- merge(human_cowARGs_pcoa_BC, human_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(human_cowARGs_pcoa_BC)[1] <- "Sample"
human_cowARGs_pcoa_BC$SeasonFull <- paste(human_cowARGs_pcoa_BC$Season, human_cowARGs_pcoa_BC$Year)

ggplot(data = human_cowARGs_pcoa_BC[which(human_cowARGs_pcoa_BC$SeasonFull == "WINTER 1"),], aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull)) + stat_ellipse(aes(color = Group)) + 
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test()

human_totalARG_cowARGs <- as.data.frame(rowSums(shortbred_wide_humans_cowARGs))
human_totalARG_cowARGs <- merge(human_totalARG_cowARGs, human_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(human_totalARG_cowARGs)[1] <- "Sample"
colnames(human_totalARG_cowARGs)[2] <- "TotalRPKM"
human_totalARG_cowARGs$SeasonFull <- paste(human_totalARG_cowARGs$Season, human_totalARG_cowARGs$Year)
human_totalARG_cowARGs$SeasonFull <- factor(human_totalARG_cowARGs$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))

compare_means(TotalRPKM ~ Group, data = human_totalARG_cowARGs, group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")

ggplot(human_totalARG_cowARGs, aes(x = Group, y = TotalRPKM)) + geom_boxplot() + geom_point(aes(color = Group)) +
  facet_wrap(~SeasonFull, nrow = 1, labeller = season_labeller) +
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test()

human_richness_cowARGs <- as.data.frame(specnumber(shortbred_wide_humans_cowARGs))
human_richness_cowARGs <- merge(human_richness_cowARGs, human_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(human_richness_cowARGs)[1] <- "Sample"
colnames(human_richness_cowARGs)[2] <- "Richness"
human_richness_cowARGs$SeasonFull <- paste(human_richness_cowARGs$Season, human_richness_cowARGs$Year)
human_richness_cowARGs$SeasonFull <- factor(human_richness_cowARGs$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))

compare_means(Richness ~ Group, data = human_richness_cowARGs, group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")

ggplot(human_richness_cowARGs, aes(x = Group, y = Richness)) + geom_boxplot() + geom_point(aes(color = Group)) +
  facet_wrap(~SeasonFull, nrow = 1, labeller = season_labeller) +
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test()

#----------------------------------------------------------------------------------------------------------------------------------

shortbred_wide_humans_cowSharedGenera <- shortbred_counts_wide_humans[,which(colnames(shortbred_counts_wide_humans) %in% cow_sharedGenera_corARG_filtered$Ogname)]
#shortbred_wide_humans_cowSharedGenera <- shortbred_counts_wide_humans[,which(colnames(shortbred_counts_wide_humans) %in% cow_sharedGenera_corARG_filtered[which(cow_sharedGenera_corARG_filtered$genus == "Bifidobacterium"),]$Ogname)]
#shortbred_wide_humans_cowSharedGenera <- shortbred_wide_humans_cowSharedGenera[which(rowSums(shortbred_wide_humans_cowSharedGenera) > 0),]
human_SharedGenera_BrayCurtis <- vegdist(shortbred_wide_humans_cowSharedGenera, index = "bray")
human_SharedGenera_pcoa_BC <- pco(human_SharedGenera_BrayCurtis, k = 4)
human_SharedGenera_pcoa_BC <- as.data.frame(human_SharedGenera_pcoa_BC$points)
human_SharedGenera_pcoa_BC <- merge(human_SharedGenera_pcoa_BC, human_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(human_SharedGenera_pcoa_BC)[1] <- "Sample"
human_SharedGenera_pcoa_BC$SeasonFull <- paste(human_SharedGenera_pcoa_BC$Season, human_SharedGenera_pcoa_BC$Year)

ggplot(data = human_SharedGenera_pcoa_BC[which(human_SharedGenera_pcoa_BC$SeasonFull == "SPRING 1"),], aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull)) + stat_ellipse(aes(color = Group)) + 
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test()

#Looking at how many times each of these genes are present in different groups (W and D)
shortbred_humans_cowSharedGenera_GroupPresence <- shortbred_wide_humans_cowSharedGenera
shortbred_humans_cowSharedGenera_GroupPresence[shortbred_humans_cowSharedGenera_GroupPresence != 0] <- 1
shortbred_humans_cowSharedGenera_GroupPresence <- merge(shortbred_humans_cowSharedGenera_GroupPresence, human_metadata[,c(3,12,13,16)], by.x = "row.names", by.y = "Sample_ID")
shortbred_humans_cowSharedGenera_GroupPresence$SeasonFull <- paste(shortbred_humans_cowSharedGenera_GroupPresence$Season, shortbred_humans_cowSharedGenera_GroupPresence$Year)
shortbred_humans_cowSharedGenera_GroupPresence <- shortbred_humans_cowSharedGenera_GroupPresence[which(shortbred_humans_cowSharedGenera_GroupPresence$SeasonFull == "SPRING 1"),]
shortbred_humans_cowSharedGenera_GroupPresence <- shortbred_humans_cowSharedGenera_GroupPresence[,c(32, 2:30)]
shortbred_humans_cowSharedGenera_GroupPresence <- shortbred_humans_cowSharedGenera_GroupPresence %>% group_by(Group) %>% summarise_all(sum)

human_totalARG_SharedGenera <- as.data.frame(rowSums(shortbred_wide_humans_cowSharedGenera))
human_totalARG_SharedGenera <- merge(human_totalARG_SharedGenera, human_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(human_totalARG_SharedGenera)[1] <- "Sample"
colnames(human_totalARG_SharedGenera)[2] <- "TotalRPKM"
human_totalARG_SharedGenera$SeasonFull <- paste(human_totalARG_SharedGenera$Season, human_totalARG_SharedGenera$Year)
human_totalARG_SharedGenera$SeasonFull <- factor(human_totalARG_SharedGenera$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))

compare_means(TotalRPKM ~ Group, data = human_totalARG_SharedGenera, group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")

ggplot(human_totalARG_SharedGenera, aes(x = Group, y = TotalRPKM)) + geom_boxplot() + geom_point(aes(color = Group)) +
  facet_wrap(~SeasonFull, nrow = 1, labeller = season_labeller) +
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test()

human_richness_SharedGenera <- as.data.frame(specnumber(shortbred_wide_humans_cowSharedGenera))
human_richness_SharedGenera <- merge(human_richness_SharedGenera, human_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(human_richness_SharedGenera)[1] <- "Sample"
colnames(human_richness_SharedGenera)[2] <- "Richness"
human_richness_SharedGenera$SeasonFull <- paste(human_richness_SharedGenera$Season, human_richness_SharedGenera$Year)
human_richness_SharedGenera$SeasonFull <- factor(human_richness_SharedGenera$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))
human_richness_SharedGenera$Group <- factor(human_richness_SharedGenera$Group, levels = c("W", "D"))
human_richness_SharedGenera_SIDaverage <- human_richness_SharedGenera[,c(2,3,14)] %>% group_by(SID, Group) %>% summarise_all(mean)

compare_means(Richness ~ Group, data = human_richness_SharedGenera[which(human_richness_SharedGenera$SeasonFull != "SPRING 2"),], group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")
compare_means(Richness ~ Group, data = human_richness_SharedGenera_SIDaverage, method = "wilcox.test", p.adjust.method = "BH")

ggplot(human_richness_SharedGenera, aes(x = Group, y = Richness)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(color = Group), position = position_jitter(width = 0.3, height = 0.1)) +
  facet_wrap(~SeasonFull, nrow = 1, labeller = season_labeller) +
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test() + labs(title = "Richness of ARGs found to be correlated with shared genera") + ylab("ARG richness")

ggplot(human_richness_SharedGenera, aes(x = Group, y = Richness)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(color = Group), position = position_dodge2(width = 0.4), alpha = 0.75) +
  facet_wrap(~SeasonFull, nrow = 1, labeller = season_labeller) +
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test() + labs(title = "Richness of ARGs found to be correlated with shared genera") + ylab("ARG richness")

ggplot(human_richness_SharedGenera, aes(x = Group, y = Richness)) + geom_boxplot(outlier.shape = NA) + geom_beeswarm(aes(color = Group), dodge.width = 0.9, cex = 3) +
  facet_wrap(~SeasonFull, nrow = 1, labeller = season_labeller) +
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test() + labs(title = "Richness of ARGs found to be correlated with shared genera") + ylab("ARG richness")

ggplot(human_richness_SharedGenera[which(human_richness_SharedGenera$Season == "SPRING"),], aes(x = Group, y = Richness, color = Group)) + geom_boxplot(aes(fill = Group), color = "black", outlier.shape = NA, width = 0.5, size = 1, alpha = 0.5, position = position_nudge(x = -0.25)) +
  geom_flat_violin(aes(fill = Group), position = position_nudge(x = 0.15), alpha = 0.7, scale = "width", width = 0.55) +
  geom_point(position = position_dodgenudge(width = 0.4, x = -0.25), alpha = 0.8) +
  facet_wrap(~SeasonFull, nrow = 1, labeller = season_labeller) +
  theme_linedraw() + labs(title = "Richness of ARGs found to be correlated with shared genera") + ylab("ARG richness") +
  scale_x_discrete(labels = c("Farm.", "Cont.")) + scale_color_manual(values = c("#18678d", "#626262")) + 
  scale_fill_manual(values = c("#18678d", "#626262")) + xlab(element_blank())
p <- ggplot(human_richness_SharedGenera[which(human_richness_SharedGenera$Season == "SPRING"),], aes(x = Group, y = Richness, color = Group)) + geom_boxplot(aes(fill = Group), color = "black", outlier.shape = NA, width = 0.6, size = 1.2, alpha = 0.5, position = position_nudge(x = -0.23)) +
  geom_flat_violin(aes(fill = Group), position = position_nudge(x = 0.12), alpha = 0.7, scale = "width", width = 0.52) +
  geom_point(position = position_dodge2nudge(width = 0.5, x = -0.23), alpha = 0.6, size = 4) +
  facet_wrap(~SeasonFull, nrow = 1, labeller = season_labeller) +
  theme_linedraw() +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 24, 5)) +
  scale_x_discrete(labels = c("Farm.", "Cont.")) + scale_color_manual(values = c("#18678d", "#626262")) + 
  scale_fill_manual(values = c("#18678d", "#626262")) + xlab(element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "bold", size = 25),
        legend.position = "none", axis.title.y = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/SharedGenera_ARGRichness_v1.svg", p, bg = "transparent", width = 5, height = 6.5)

p <- ggplot(human_richness_SharedGenera[which(human_richness_SharedGenera$SeasonFull == "SPRING 1"),], aes(x = Group, y = Richness, color = Group)) + geom_boxplot(aes(fill = Group), color = "black", outlier.shape = NA, width = 0.6, size = 1.2, alpha = 0.5, position = position_nudge(x = -0.1)) +
  geom_flat_violin(aes(fill = Group), position = position_nudge(x = 0.25), alpha = 0.7, scale = "width", width = 0.5) +
  geom_point(position = position_dodge2nudge(width = 0.5, x = -0.1), alpha = 0.6, size = 4) +
  theme_linedraw() +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 24, 5)) +
  scale_x_discrete(labels = c("Farm.", "Cont.")) + scale_color_manual(values = c("#18678d", "#626262")) + 
  scale_fill_manual(values = c("#18678d", "#626262")) + xlab(element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "bold", size = 25),
        legend.position = "none", axis.title.y = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/SharedGenera_ARGRichness_v2.svg", p, bg = "transparent", width = 3, height = 6.5)

#----------------------------------------------------------------------------------------------------------------------------------

shortbred_wide_humans_cowAllSharedGenera <- shortbred_counts_wide_humans[,which(colnames(shortbred_counts_wide_humans) %in% cow_allSharedGenera_corARG$gene)]

human_totalARG_AllSharedGenera <- as.data.frame(rowSums(shortbred_wide_humans_cowAllSharedGenera))
human_totalARG_AllSharedGenera <- merge(human_totalARG_AllSharedGenera, human_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(human_totalARG_AllSharedGenera)[1] <- "Sample"
colnames(human_totalARG_AllSharedGenera)[2] <- "TotalRPKM"
human_totalARG_AllSharedGenera$SeasonFull <- paste(human_totalARG_AllSharedGenera$Season, human_totalARG_AllSharedGenera$Year)
human_totalARG_AllSharedGenera$SeasonFull <- factor(human_totalARG_AllSharedGenera$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))
human_totalARG_AllSharedGenera_SIDAverage <- human_totalARG_AllSharedGenera[,c(2,3,14)] %>% group_by(SID, Group) %>% summarise_all(mean)

compare_means(TotalRPKM ~ Group, data = human_totalARG_AllSharedGenera[which(human_totalARG_AllSharedGenera$SeasonFull != "SPRING 2"),], group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")
compare_means(TotalRPKM ~ Group, data = human_totalARG_AllSharedGenera_SIDAverage, method = "wilcox.test", p.adjust.method = "BH")

ggplot(human_totalARG_AllSharedGenera, aes(x = Group, y = TotalRPKM)) + geom_boxplot() + geom_point(aes(color = Group)) +
  facet_wrap(~SeasonFull, nrow = 1, labeller = season_labeller) +
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test()
ggplot(human_totalARG_AllSharedGenera_SIDAverage, aes(x = Group, y = TotalRPKM)) + geom_boxplot() + geom_point(aes(color = Group)) +
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test()


human_richness_AllSharedGenera <- as.data.frame(specnumber(shortbred_wide_humans_cowAllSharedGenera))
human_richness_AllSharedGenera <- merge(human_richness_AllSharedGenera, human_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(human_richness_AllSharedGenera)[1] <- "Sample"
colnames(human_richness_AllSharedGenera)[2] <- "Richness"
human_richness_AllSharedGenera$SeasonFull <- paste(human_richness_AllSharedGenera$Season, human_richness_AllSharedGenera$Year)
human_richness_AllSharedGenera$SeasonFull <- factor(human_richness_AllSharedGenera$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))
human_richness_AllSharedGenera$Group <- factor(human_richness_AllSharedGenera$Group, levels = c("W", "D"))
human_richness_AllSharedGenera_SIDaverage <- human_richness_AllSharedGenera[,c(2,3,14)] %>% group_by(SID, Group) %>% summarise_all(mean)

compare_means(Richness ~ Group, data = human_richness_AllSharedGenera[which(human_richness_SharedGenera$SeasonFull != "SPRING 2"),], group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")
compare_means(Richness ~ Group, data = human_richness_AllSharedGenera_SIDaverage, method = "wilcox.test", p.adjust.method = "BH")

ggplot(human_richness_AllSharedGenera, aes(x = Group, y = Richness)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(color = Group), position = position_jitter(width = 0.3, height = 0.1)) +
  facet_wrap(~SeasonFull, nrow = 1, labeller = season_labeller) +
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test() + labs(title = "Richness of ARGs found to be correlated with all hared genera") + ylab("ARG richness")

ggplot(human_richness_AllSharedGenera_SIDaverage, aes(x = Group, y = Richness)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(color = Group), position = position_jitter(width = 0.3, height = 0.1)) +
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test() + labs(title = "Richness of ARGs found to be correlated with all shared genera") + ylab("ARG richness")
p <- ggplot(human_richness_AllSharedGenera_SIDaverage, aes(x = Group, y = Richness, color = Group)) + geom_boxplot(aes(fill = Group), color = "black", outlier.shape = NA, width = 0.6, size = 1.2, alpha = 0.5, position = position_nudge(x = -0.1)) +
  geom_flat_violin(aes(fill = Group), position = position_nudge(x = 0.25), alpha = 0.7, scale = "width", width = 0.5) +
  geom_point(position = position_dodge2nudge(width = 0.5, x = -0.1), alpha = 0.6, size = 4) +
  theme_linedraw() +
  scale_y_continuous(limits = c(NA, 120), breaks = seq(0, 120, 25)) +
  scale_x_discrete(labels = c("Farm.", "Cont.")) + scale_color_manual(values = c("#18678d", "#626262")) + 
  scale_fill_manual(values = c("#18678d", "#626262")) + xlab(element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text = element_blank(), strip.text = element_text(face = "bold", size = 25),
        legend.position = "none", axis.title.y = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/SharedGenera_ARGRichness_v4.svg", p, bg = "transparent", width = 5, height = 6.5)


#----------------------------------------------------------------------------------------------------------------------------------

shortbred_wide_humans_cowSharedGenera_03 <- shortbred_counts_wide_humans[,which(colnames(shortbred_counts_wide_humans) %in% cow_sharedGenera_corARG_03$cow_merged_SharedGenera.gene)]

human_totalARG_SharedGenera_03 <- as.data.frame(rowSums(shortbred_wide_humans_cowSharedGenera_03))
human_totalARG_SharedGenera_03 <- merge(human_totalARG_SharedGenera_03, human_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(human_totalARG_SharedGenera_03)[1] <- "Sample"
colnames(human_totalARG_SharedGenera_03)[2] <- "TotalRPKM"
human_totalARG_SharedGenera_03$SeasonFull <- paste(human_totalARG_SharedGenera_03$Season, human_totalARG_SharedGenera_03$Year)
human_totalARG_SharedGenera_03$SeasonFull <- factor(human_totalARG_SharedGenera_03$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))

compare_means(TotalRPKM ~ Group, data = human_totalARG_SharedGenera_03, group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")

ggplot(human_totalARG_SharedGenera_03, aes(x = Group, y = TotalRPKM)) + geom_boxplot() + geom_point(aes(color = Group)) +
  facet_wrap(~SeasonFull, nrow = 1, labeller = season_labeller) +
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test()


human_richness_SharedGenera_03 <- as.data.frame(specnumber(shortbred_wide_humans_cowSharedGenera_03))
human_richness_SharedGenera_03 <- merge(human_richness_SharedGenera_03, human_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(human_richness_SharedGenera_03)[1] <- "Sample"
colnames(human_richness_SharedGenera_03)[2] <- "Richness"
human_richness_SharedGenera_03$SeasonFull <- paste(human_richness_SharedGenera_03$Season, human_richness_SharedGenera_03$Year)
human_richness_SharedGenera_03$SeasonFull <- factor(human_richness_SharedGenera_03$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))
human_richness_SharedGenera_03$Group <- factor(human_richness_SharedGenera_03$Group, levels = c("W", "D"))
human_richness_SharedGenera_03_SIDaverage <- human_richness_SharedGenera_03[,c(2,3,14)] %>% group_by(SID, Group) %>% summarise_all(mean)

compare_means(Richness ~ Group, data = human_richness_SharedGenera_03, group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")
compare_means(Richness ~ Group, data = human_richness_SharedGenera_03_SIDaverage, method = "wilcox.test", p.adjust.method = "BH")

ggplot(human_richness_SharedGenera_03, aes(x = Group, y = Richness)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(color = Group), position = position_jitter(width = 0.3, height = 0.1)) +
  facet_wrap(~SeasonFull, nrow = 1, labeller = season_labeller) +
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test() + labs(title = "Richness of ARGs found to be correlated with shared genera") + ylab("ARG richness")

#----------------------------------------------------------------------------------------------------------------------------------

sharedARGs <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d22_BLASTn_HQMQ/CW_100cov_sharedARGs.txt", check.names = FALSE, header = FALSE)
colnames(sharedARGs)[1] <- "ARG"

shortbred_wide_humans_sharedARGs <- shortbred_counts_wide_humans[,which(colnames(shortbred_counts_wide_humans) %in% sharedARGs$ARG)]

human_richness_sharedARGs <- as.data.frame(specnumber(shortbred_wide_humans_sharedARGs))
human_richness_sharedARGs <- merge(human_richness_sharedARGs, human_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(human_richness_sharedARGs)[1] <- "Sample"
colnames(human_richness_sharedARGs)[2] <- "Richness"
human_richness_sharedARGs$SeasonFull <- paste(human_richness_sharedARGs$Season, human_richness_sharedARGs$Year)
human_richness_sharedARGs$SeasonFull <- factor(human_richness_sharedARGs$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))
human_richness_sharedARGs$Group <- factor(human_richness_sharedARGs$Group, levels = c("W", "D"))
human_richness_SharedGenera_SIDaverage <- human_richness_sharedARGs[,c(2,3,14)] %>% group_by(SID, Group) %>% summarise_all(mean)

compare_means(Richness ~ Group, data = human_richness_sharedARGs, group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")


sharedARGs_enrichedCows <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d22_BLASTn_HQMQ/CW_100cov_sharedARGs_EnrichedInCows.txt", check.names = FALSE, header = FALSE)
colnames(sharedARGs_enrichedCows)[1] <- "ARG"

shortbred_wide_humans_sharedARGs_enrichedCows <- shortbred_counts_wide_humans[,which(colnames(shortbred_counts_wide_humans) %in% sharedARGs_enrichedCows$ARG)]

human_richness_sharedARGs_enrichedCows <- as.data.frame(specnumber(shortbred_wide_humans_sharedARGs_enrichedCows))
human_richness_sharedARGs_enrichedCows <- merge(human_richness_sharedARGs_enrichedCows, human_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(human_richness_sharedARGs_enrichedCows)[1] <- "Sample"
colnames(human_richness_sharedARGs_enrichedCows)[2] <- "Richness"
human_richness_sharedARGs_enrichedCows$SeasonFull <- paste(human_richness_sharedARGs_enrichedCows$Season, human_richness_sharedARGs_enrichedCows$Year)
human_richness_sharedARGs_enrichedCows$SeasonFull <- factor(human_richness_sharedARGs_enrichedCows$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))
human_richness_sharedARGs_enrichedCows$Group <- factor(human_richness_sharedARGs_enrichedCows$Group, levels = c("W", "D"))
human_richness_SharedGenera_SIDaverage <- human_richness_sharedARGs_enrichedCows[,c(2,3,14)] %>% group_by(SID, Group) %>% summarise_all(mean)

compare_means(Richness ~ Group, data = human_richness_sharedARGs_enrichedCows, group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")


human_totalARG_sharedARGs_enrichedCows <- as.data.frame(rowSums(shortbred_wide_humans_sharedARGs_enrichedCows))
human_totalARG_sharedARGs_enrichedCows <- merge(human_totalARG_sharedARGs_enrichedCows, human_metadata, by.x = "row.names", by.y = "Sample_ID")
colnames(human_totalARG_sharedARGs_enrichedCows)[1] <- "Sample"
colnames(human_totalARG_sharedARGs_enrichedCows)[2] <- "TotalRPKM"
human_totalARG_sharedARGs_enrichedCows$SeasonFull <- paste(human_totalARG_sharedARGs_enrichedCows$Season, human_totalARG_sharedARGs_enrichedCows$Year)
human_totalARG_sharedARGs_enrichedCows$SeasonFull <- factor(human_totalARG_sharedARGs_enrichedCows$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))

compare_means(TotalRPKM ~ Group, data = human_totalARG_sharedARGs_enrichedCows, group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")

ggplot(human_totalARG_sharedARGs_enrichedCows, aes(x = Group, y = TotalRPKM)) + geom_boxplot() + geom_point(aes(color = Group)) +
  facet_wrap(~SeasonFull, nrow = 1, labeller = season_labeller) +
  scale_color_manual(values=c("#00BA38", "#619CFF")) + 
  theme_test()

###Looking at the correlation of the cow ARGs
shortbred_wide_cowARGs <- shortbred_counts_wide[samples_metadata[which(samples_metadata$Group == "C"),"Sample_ID"],]
shortbred_wide_cowARGs <- shortbred_wide_cowARGs[,colSums(shortbred_wide_cowARGs) > 0]
cowARGs_corr <- rcorr(as.matrix(shortbred_wide_cowARGs))
cowARGs_corr <- flattenCorrMatrix(cowARGs_corr$r, cowARGs_corr$P)
colnames(cowARGs_corr)[1] <- "ARG1"
colnames(cowARGs_corr)[2] <- "ARG2"
colnames(cowARGs_corr)[3] <- "R"
colnames(cowARGs_corr)[4] <- "Pval"
cowARGs_corr$Padj <- p.adjust(cowARGs_corr$Pval, method = "BH")
cowARGs_corr <- merge(cowARGs_corr, ARG_mapping, by.x = "ARG1", by.y = "Gene")
colnames(cowARGs_corr)[6] <- "AbxClass1"
cowARGs_corr <- merge(cowARGs_corr, ARG_mapping, by.x = "ARG2", by.y = "Gene")
colnames(cowARGs_corr)[7] <- "AbxClass2"
cowARGs_corr_filt <- cowARGs_corr[which(cowARGs_corr$Padj < 0.05 & cowARGs_corr$R > 0.7),]
write.csv(cowARGs_corr_filt, 
          file = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_AnalyzeShortBRED/PearsonCorr/230201_AllCowARGsPearsonCorr_07_005BH.csv",
          row.names = FALSE)

cowARGs_corr_filt <- cowARGs_corr[which(cowARGs_corr$Padj < 0.05 & cowARGs_corr$R > 0.5),]
cowARGs_corr_filt <- dcast(cowARGs_corr_filt[,c(2,1,3)], ARG1 ~ ARG2, value.var = "R")
rownames(cowARGs_corr_filt) <- cowARGs_corr_filt$ARG1
cowARGs_corr_filt <- cowARGs_corr_filt[,c(2:64)]
cowARGs_corr_filt[is.na(cowARGs_corr_filt)] <- 0
test <- sapply(cowARGs_corr_filt, as.numeric)
rownames(test) <- rownames(cowARGs_corr_filt)
cowARGs_corr_filt <- test

heatmap(cowARGs_corr_filt)


dim(cowARGs_corr[which(cowARGs_corr$Padj < 0.05 & cowARGs_corr$R > 0.3),])

###Looking at Lnu(AN2) and Mef(EN2) genes
shortbred_LnuMef <- shortbred_combined[which(shortbred_combined$Family %in% c("490437755|WP_004308783_1|1|1|lnu(AN2)|lnu(AN2)|nucleotidyltransferase|2|LINCOSAMIDE|LINCOSAMIDE|lincosamide_nucleotidyltransferase_Lnu(AN2)", "1028097032|WP_063853729_1|1|1|mef(En2)|mef(En2)|efflux|2|MACROLIDE|MACROLIDE|macrolide_efflux_MFS_transporter_Mef(En2)")),]
shortbred_LnuMef <- merge(shortbred_LnuMef, samples_metadata, by.x = "Sample", by.y = "Sample_ID")
shortbred_LnuMef$SeasonFull <- paste(shortbred_LnuMef$Season, shortbred_LnuMef$Year)
shortbred_LnuMef$SeasonFull <- factor(shortbred_LnuMef$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))
shortbred_LnuMef$Group <- factor(shortbred_LnuMef$Group, levels = c("C", "W", "D"))
shortbred_LnuMef_meanSID <- shortbred_LnuMef[which(shortbred_LnuMef$Group != "C" & shortbred_LnuMef$SeasonFull != "SPRING 2"),c(1, 3, 6, 17)]
shortbred_LnuMef_meanSID <- shortbred_LnuMef_meanSID %>% group_by(Sample, Group, SID) %>% summarise(Count = sum(Count))
shortbred_LnuMef_meanSID <- shortbred_LnuMef_meanSID[,c(2,3,4)] %>% group_by(Group, SID) %>% summarise(Count = mean(Count))

compare_means(Count ~ Group,
              data = shortbred_LnuMef_meanSID, method = "t.test", p.adjust.method = "BH")

test <- compare_means(Count ~ Group, 
                      data = shortbred_LnuMef[which(shortbred_LnuMef$Family == "1028097032|WP_063853729_1|1|1|mef(En2)|mef(En2)|efflux|2|MACROLIDE|MACROLIDE|macrolide_efflux_MFS_transporter_Mef(En2)" & shortbred_LnuMef$Group != "C" & shortbred_LnuMef$SeasonFull != "SPRING 2"),],
                      group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")
test <- test[which(test$p.adj < 0.05),]

ggplot(data = shortbred_LnuMef[which(shortbred_LnuMef$Family == "490437755|WP_004308783_1|1|1|lnu(AN2)|lnu(AN2)|nucleotidyltransferase|2|LINCOSAMIDE|LINCOSAMIDE|lincosamide_nucleotidyltransferase_Lnu(AN2)"),],
       aes(x = Group, y = Count, color = Group)) +
  geom_boxplot() +
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) +
  theme_linedraw()

ggplot(data = shortbred_LnuMef_meanSID, aes(x = Group, y = Count, color = Group)) + 
  geom_boxplot()

###Looking at the taxa (MAGs) that carry the Lnu(AN2) and Mef(EN2) genes
MAG_subjects <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d07_bin/dastool/allbins_noUnbin/MAG_subjects.csv")
MAG_taxa <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d07_bin/dastool/HQandMQ_MAGs/joint_classify.csv")

shortbred_HQMQ <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d12_shortbred_HQMQ/joint_HQMQ_output.txt", check.names = FALSE)
shortbred_HQMQ_LnuMef <- shortbred_HQMQ[which(shortbred_HQMQ$Family %in% c("490437755|WP_004308783_1|1|1|lnu(AN2)|lnu(AN2)|nucleotidyltransferase|2|LINCOSAMIDE|LINCOSAMIDE|lincosamide_nucleotidyltransferase_Lnu(AN2)", "1028097032|WP_063853729_1|1|1|mef(En2)|mef(En2)|efflux|2|MACROLIDE|MACROLIDE|macrolide_efflux_MFS_transporter_Mef(En2)")),]
shortbred_HQMQ_LnuMef <- merge(shortbred_HQMQ_LnuMef, MAG_subjects, by = "MAG")
shortbred_HQMQ_LnuMef <- merge(shortbred_HQMQ_LnuMef, MAG_taxa, by = "MAG")

unique(shortbred_HQMQ_LnuMef[which(shortbred_HQMQ_LnuMef$Group == "C"),c("Taxon")])
unique(shortbred_HQMQ_LnuMef[which(shortbred_HQMQ_LnuMef$Group == "W"),c("Taxon")])
unique(shortbred_HQMQ_LnuMef[which(shortbred_HQMQ_LnuMef$Group == "D"),c("Taxon")])
