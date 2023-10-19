###This script is to analyze the 16S data generate by Rhiannon

###Load the packages
library(ggplot2)
library(vegan)
library(dplyr)
library(labdsv)
library(reshape2)
library(stringr)
library(rstatix)
library(ggpubr)
library(Maaslin2)
library(NatParksPalettes)
library(PupillometryR)
library(ggpattern)
library(ggsci)
library(ggpp)
library(svglite)

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

group_names <- list(
  'C' = 'Cows',
  'W' = 'Farmers',
  'D' = 'Controls'
)

group_labeller <- function(variable, value) {
  return(group_names[value])
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

###Read in the metadata
samples_metadata <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/220504_DOME_SamplesMasterList_v10.csv", header = TRUE)
#Two of the cow samples (19-C0205-W001fc01 & 19-C0416-F001fc01) have been clustering with human samples when looking at taxa (shotgun and 16S) and ARGs
#These two samples have likley been misplabeled. Removing these two from further analysis
#Another cow sample (19-C0210-S001fc01) is clustering with human ones, so removing that as well
#Sample 19-C0327-V301ns01 also seem to be clustering with the human groups. Removing these as well
samples_metadata <- samples_metadata[which(samples_metadata$Sample_ID != "19-C0205-W001fc01" & samples_metadata$Sample_ID != "19-C0416-F001fc01"),]
samples_metadata <- samples_metadata[which(samples_metadata$Sample_ID != "19-C0210-S001fc01" & samples_metadata$Sample_ID != "19-C0327-V301ns01"),]
samples_metadata <- samples_metadata[,1:16]

###Read in the csv file with the genera read counts and generating relative abundance matrix with 0.1% abundance cutoff
genus_abund <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/ASV_tables/230209_Genus_editedSampleID.csv", header = TRUE, row.names = 1)
genus_abund[is.na(genus_abund)] <- 0
genus_abund <- genus_abund[which(rowSums(genus_abund) > 0),]
genus_abund <- genus_abund/rowSums(genus_abund)
genus_abund[is.na(genus_abund)] <- 0
genus_abund[genus_abund < 0.001] <- 0
genus_abund <- genus_abund/rowSums(genus_abund)
genus_abund[is.na(genus_abund)] <- 0
genus_abund <- genus_abund[which(rownames(genus_abund) %in% samples_metadata$Sample_ID),]
genus_abund <- genus_abund[, colSums(genus_abund) > 0]
genus_abund_ns <- genus_abund[which(rownames(genus_abund) %in% samples_metadata[which(samples_metadata$Sample_Type == "Nasal"),]$Sample_ID),]
genus_abund_ns <- genus_abund_ns[, colSums(genus_abund_ns) > 0]
genus_abund_ns_2019 <- genus_abund[which(rownames(genus_abund) %in% samples_metadata[which(samples_metadata$Sample_Type == "Nasal" & samples_metadata$Year != 2),]$Sample_ID),]
genus_abund_ns_2019 <- genus_abund_ns_2019[, colSums(genus_abund_ns_2019) > 0]
genus_abund_fc <- genus_abund[which(rownames(genus_abund) %in% samples_metadata[which(samples_metadata$Sample_Type == "Fecal"),]$Sample_ID),]
genus_abund_fc <- genus_abund_fc[, colSums(genus_abund_fc) > 0]
genus_abund_fc_2019 <- genus_abund[which(rownames(genus_abund) %in% samples_metadata[which(samples_metadata$Sample_Type == "Fecal" & samples_metadata$Year != 2),]$Sample_ID),]
genus_abund_fc_2019 <- genus_abund_fc_2019[, colSums(genus_abund_fc_2019) > 0]

###Read in the csv file with the family read counts and generating relative abundance matrix with 0.1% abundance cutoff
family_abund <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/ASV_tables/230209_Family_editedSamplesID.csv", header = TRUE, row.names = 1)
family_abund[is.na(family_abund)] <- 0
family_abund <- family_abund[which(rowSums(family_abund) > 0),]
family_abund <- family_abund/rowSums(family_abund)
family_abund[is.na(family_abund)] <- 0
family_abund[family_abund < 0.001] <- 0
family_abund <- family_abund/rowSums(family_abund)
family_abund[is.na(family_abund)] <- 0
family_abund <- family_abund[which(rownames(family_abund) %in% samples_metadata$Sample_ID),]
family_abund <- family_abund[, colSums(family_abund) > 0]
family_abund_ns <- family_abund[which(rownames(family_abund) %in% samples_metadata[which(samples_metadata$Sample_Type == "Nasal"),]$Sample_ID),]
family_abund_ns <- family_abund_ns[, colSums(family_abund_ns) > 0]
family_abund_fc <- family_abund[which(rownames(family_abund) %in% samples_metadata[which(samples_metadata$Sample_Type == "Fecal"),]$Sample_ID),]
family_abund_fc <- family_abund_fc[, colSums(family_abund_fc) > 0]

###Get the genus alpha diversity for all samples
richness <- as.data.frame(specnumber(genus_abund))
richness$Sample <- rownames(richness)
colnames(richness)[1] <-"Richness"

shannon_div <- as.data.frame(diversity(genus_abund, index = "shannon"))
shannon_div$Sample <- row.names(shannon_div)
colnames(shannon_div)[1] <- "ShannonDiv"

alpha_div <- merge(richness, shannon_div, by = "Sample")
alpha_div <- merge(alpha_div, samples_metadata[, c(1,2,3,12,13,16)], by.x = "Sample", by.y = "Sample_ID")
alpha_div$SeasonFull <- paste(alpha_div$Season, alpha_div$Year)
alpha_div$SeasonFull <- factor(alpha_div$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))
alpha_div$Group <- factor(alpha_div$Group, levels = c("C", "W", "D"))

test <- compare_means(ShannonDiv ~ Group, data = alpha_div[which(alpha_div$Sample_Type == "Nasal" & alpha_div$SeasonFull != "SPRING 2"),], group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")
test <- test[which(test$p.adj < 0.05),]

ggplot(alpha_div[which(alpha_div$Sample_Type == "Fecal"),], aes(x = Group, y = Richness, color = Group)) + geom_boxplot() +
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_test() +
  theme(legend.position = "none") + ylab("Genus Richness (0.1% threshold)") + labs(title = "Fecal samples. 16S")
p <- ggplot(alpha_div[which(alpha_div$Sample_Type == "Fecal"),], aes(x = Group, y = Richness, color = Group)) + geom_boxplot(aes(fill = Group), color = "black", outlier.shape = NA, width = 0.4, size = 1, alpha = 0.5) + 
  geom_flat_violin(aes(fill = Group), position = position_nudge(x = 0.3), alpha = 0.7, scale = "width", width = 0.45) +
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.8) +
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_linedraw() + 
  scale_y_continuous(limits = c(NA, 105), breaks = seq(0, 100, 20)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(linewidth = 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "bold", size = 20),
        legend.position = "none", axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Cows", "Farm.", "Cont.")) + scale_color_manual(values = c("#873e23", "#18678d", "#626262")) + 
  scale_fill_manual(values = c("#873e23", "#18678d", "#626262")) + xlab(element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/16S_Fecal_Richness_v1.svg", p, bg = "transparent", width = 15, height = 10)

ggplot(alpha_div[which(alpha_div$Sample_Type == "Nasal"),], aes(x = Group, y = Richness, color = Group)) + geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.8) +
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_linedraw() + ylim(NA, 175) +
  theme(legend.position = "none") + ylab("Genus Richness (0.1% threshold)") + labs(title = "Nasal samples. 16S") +
  scale_x_discrete(labels = c("Cows", "Farmers", "Controls")) + scale_color_manual(values = c("#873e23", "#18678d", "#626262"))
p <- ggplot(alpha_div[which(alpha_div$Sample_Type == "Nasal"),], aes(x = Group, y = Richness, color = Group)) + geom_boxplot(aes(fill = Group), color = "black", outlier.shape = NA, width = 0.4, size = 1, alpha = 0.5) + 
  geom_flat_violin(aes(fill = Group), position = position_nudge(x = 0.3), alpha = 0.7, scale = "width", width = 0.45) +
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.8) +
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_linedraw() + ylim(NA, 160) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(linewidth = 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "bold", size = 20),
        legend.position = "none", axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Cows", "Farm.", "Cont.")) + scale_color_manual(values = c("#873e23", "#18678d", "#626262")) + 
  scale_fill_manual(values = c("#873e23", "#18678d", "#626262")) + xlab(element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/16S_Nasal_Richness_v1.svg", p, bg = "transparent", width = 15, height = 10)

p <- ggplot(alpha_div[which(alpha_div$Sample_Type == "Fecal"),], aes(x = Group, y = ShannonDiv, color = Group)) + geom_boxplot(aes(fill = Group), color = "black", outlier.shape = NA, width = 0.4, size = 1, alpha = 0.5) +
  geom_flat_violin(aes(fill = Group), position = position_nudge(x = 0.3), alpha = 0.7, scale = "width", width = 0.45) +
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.8) +
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_linedraw() + ylim(NA, 4.2) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16, face = "bold"), strip.text = element_text(face = "bold", size = 20),
        legend.position = "none", axis.title.y = element_text(face = "bold", size = 18)) + ylab("Shannon diversity (genera)") +
  scale_x_discrete(labels = c("Cows", "Farm.", "Cont.")) + scale_color_manual(values = c("#873e23", "#18678d", "#626262")) + 
  scale_fill_manual(values = c("#873e23", "#18678d", "#626262")) + xlab(element_blank())
p <- ggplot(alpha_div[which(alpha_div$Sample_Type == "Fecal"),], aes(x = Group, y = ShannonDiv, color = Group)) + geom_boxplot(aes(fill = Group), color = "black", outlier.shape = NA, width = 0.4, size = 1, alpha = 0.5) +
  geom_flat_violin(aes(fill = Group), position = position_nudge(x = 0.3), alpha = 0.7, scale = "width", width = 0.45) +
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.8) +
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_linedraw() + ylim(NA, 4.2) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "bold", size = 20),
        legend.position = "none", axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Cows", "Farm.", "Cont.")) + scale_color_manual(values = c("#873e23", "#18678d", "#626262")) + 
  scale_fill_manual(values = c("#873e23", "#18678d", "#626262")) + xlab(element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/16S_Fecal_ShannonDiv_v2.svg", p, bg = "transparent", width = 15, height = 10)
#  scale_x_discrete(labels = c("Cow", "Office", "Farmer"))

ggplot(alpha_div[which(alpha_div$Sample_Type == "Nasal"),], aes(x = Group, y = ShannonDiv, color = Group)) + geom_boxplot(aes(fill = Group), color = "black", outlier.shape = NA, width = 0.4, size = 1, alpha = 0.5) + 
  geom_flat_violin(aes(fill = Group), position = position_nudge(x = 0.3), alpha = 0.7, scale = "width", width = 0.45) +
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.8) +
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_linedraw() + ylim(NA, 5) + 
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16, face = "bold"), strip.text = element_text(face = "bold", size = 20),
        legend.position = "none", axis.title.y = element_text(face = "bold", size = 18)) + ylab("Shannon diversity (genera)") +
  scale_x_discrete(labels = c("Cows", "Farm.", "Cont.")) + scale_color_manual(values = c("#873e23", "#18678d", "#626262")) + 
  scale_fill_manual(values = c("#873e23", "#18678d", "#626262")) + xlab(element_blank())
p <- ggplot(alpha_div[which(alpha_div$Sample_Type == "Nasal"),], aes(x = Group, y = ShannonDiv, color = Group)) + geom_boxplot(aes(fill = Group), color = "black", outlier.shape = NA, width = 0.4, size = 1, alpha = 0.5) + 
  geom_flat_violin(aes(fill = Group), position = position_nudge(x = 0.3), alpha = 0.7, scale = "width", width = 0.45) +
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.8) +
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_linedraw() + ylim(NA, 5) + 
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "bold", size = 20),
        legend.position = "none", axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Cows", "Farm.", "Cont.")) + scale_color_manual(values = c("#873e23", "#18678d", "#626262")) + 
  scale_fill_manual(values = c("#873e23", "#18678d", "#626262")) + xlab(element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/16S_Nasal_ShannonDiv_v2.svg", p, bg = "transparent", width = 15, height = 10)

###Getting the stacked bar plots
genus_abund_ns_top20 <- genus_abund_ns[,order(colSums(genus_abund_ns, na.rm = TRUE), decreasing = TRUE)]
top20genera_ns <- colnames(genus_abund_ns_top20)[1:20]
genus_abund_ns_top20 <- data.frame(genus_abund_ns_top20[,colnames(genus_abund_ns_top20) %in% top20genera_ns],Others=rowSums(genus_abund_ns_top20[,!colnames(genus_abund_ns_top20) %in% top20genera_ns]))
genus_abund_ns_top20 <- melt(as.matrix(genus_abund_ns_top20))
colnames(genus_abund_ns_top20)[1] <- "Sample"
colnames(genus_abund_ns_top20)[2] <- "Genus"
colnames(genus_abund_ns_top20)[3] <- "RelAbund"
genus_abund_ns_top20 <- merge(genus_abund_ns_top20, samples_metadata[, c(3,12,13,16)], by.x = "Sample", by.y = "Sample_ID")
genus_abund_ns_top20$SeasonFull <- paste(genus_abund_ns_top20$Season, genus_abund_ns_top20$Year)
genus_abund_ns_top20 <- genus_abund_ns_top20[, c(5,7,2,3)] %>% group_by(Group, SeasonFull, Genus) %>% summarise(RelAbund = mean(RelAbund))
genus_abund_ns_top20$SeasonFull <- factor(genus_abund_ns_top20$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))
genus_abund_ns_top20$Group <- factor(genus_abund_ns_top20$Group, levels = c("C", "W", "D"))
genus_abund_ns_top20$Genus <- factor(genus_abund_ns_top20$Genus, levels = c(top20genera_ns, "Others"))

ggplot(data = genus_abund_ns_top20, aes(x = Group, y = RelAbund, fill = Genus)) + geom_bar(stat = "identity", color = "white", size = 0.1) +
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) +
  theme_linedraw() + scale_fill_manual(values = c(c(natparks.pals("Triglav", 20)), "#7d7d7d")) +
  scale_x_discrete(labels = c("Cows", "Farm.", "Cont.")) + 
  guides(fill = guide_legend(ncol = 1)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16, face = "bold"), strip.text = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 18), axis.title.x = element_blank()) + ylab("Relative abundance")
ggplot(data = genus_abund_ns_top20, aes(x = Group, y = RelAbund, fill = Genus)) + geom_bar(stat = "identity", color = "white", size = 0.1) +
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) +
  theme_linedraw() + scale_fill_manual(values = c("#2574b0", "#b5c3e1", "#ff7d0c", "#ffbd77", "#2f9c2c", "#94e08a", "#d3262d", "#fe9a9b", "#9867be", "#c3b1da",
                                                  "#8d5749", "#c09f92", "#e079c1", "#f9b4d0", "#532A34", "#c4c9cf", "#b7bc2a", "#ddda8a", "#1abfcd", "#9cdae6", "#7d7d7d")) +
  scale_x_discrete(labels = c("Cows", "Farm.", "Cont.")) + 
  guides(fill = guide_legend(ncol = 1)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "bold", size = 20),
        axis.title.y = element_blank(), axis.title.x = element_blank())
p <- ggplot(data = genus_abund_ns_top20, aes(x = Group, y = RelAbund, fill = Genus)) + geom_bar(stat = "identity", color = "white", size = 0.1) +
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) +
  theme_linedraw() + scale_fill_manual(values = c(c(natparks.pals("Triglav", 20)), "#7d7d7d")) +
  scale_x_discrete(labels = c("Cows", "Farm.", "Cont.")) + 
  guides(fill = guide_legend(ncol = 1)) + scale_y_continuous(position = "right") +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "bold", size = 20),
        axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = "none")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/16S_Nasal_RelAbund_v2.svg", p, bg = "transparent", width = 15, height = 10)

ggplot(data = genus_abund_ns_top20, aes(x = Group, y = RelAbund, fill = Genus)) + geom_bar(stat = "identity", color = "white", size = 0.1) +
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) +
  theme_linedraw() + scale_fill_manual(values = c(sample(sample(natparks.pals("Cuyahoga", 20))), "#7d7d7d"))

ggplot(data = genus_abund_ns_top20, aes(x = Group, y = RelAbund, fill = Genus, pattern = Genus)) + geom_bar_pattern(stat = "identity") +
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) +
  theme_linedraw() + scale_fill_manual(values = c(natparks.pals("Cuyahoga", 20), "#7d7d7d")) +
  scale_pattern_manual(values = c(Staphylococcus = "none", Corynebacterium = "stripe", Anaerococcus = "none", Weissella = "stripe", Streptococcus = "none",
                                  Moraxella = "stripe", Peptoniphilus = "none", Finegoldia = "stripe", Romboutsia = "none", Mannheimia = "stripe",
                                  Acetobacter = "none", Psychrobacter = "stripe", UCG.005 = "none", Rothia = "stripe", Acinetobacter = "none",
                                  Lentilactobacillus = "stripe", Latilactobacillus = "none", Aerococcus = "stripe", Jeotgalicoccus = "none", Bacillus = "stripe",
                                  Others = "none"))


###Calculate beta diversity for all samples
all_Genus_BrayCurtis <- vegdist(genus_abund, method = "bray")
pcoa_genusAll_BC <- pco(all_Genus_BrayCurtis, k = 4)
pcoa_genusAll_BC <- as.data.frame(pcoa_genusAll_BC$points)
pcoa_genusAll_BC <- merge(pcoa_genusAll_BC, samples_metadata[, c(1,2,3,12,13,16)], by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_genusAll_BC)[1] <- "Sample"
pcoa_genusAll_BC$SeasonFull <- paste(pcoa_genusAll_BC$Season, pcoa_genusAll_BC$Year)

ggplot(data = pcoa_genusAll_BC, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = Sample_Type)) +
  theme_test()


fc_Genus_BrayCurtis <- vegdist(genus_abund_fc, method = "bray")
pcoa_FCgenus_BC <- pco(fc_Genus_BrayCurtis, k = 4)
pcoa_FCgenus_scree <- as.data.frame(pcoa_FCgenus_BC$eig*100/sum(pcoa_FCgenus_BC$eig))
pcoa_FCgenus_BC <- as.data.frame(pcoa_FCgenus_BC$points)
pcoa_FCgenus_BC <- merge(pcoa_FCgenus_BC, samples_metadata[, c(1,2,3,12,13,16)], by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_FCgenus_BC)[1] <- "Sample"
pcoa_FCgenus_BC$SeasonFull <- paste(pcoa_FCgenus_BC$Season, pcoa_FCgenus_BC$Year)
pcoa_FCgenus_BC$Season <- factor(pcoa_FCgenus_BC$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))
pcoa_FCgenus_BC$Group <- factor(pcoa_FCgenus_BC$Group, levels = c("C", "W", "D"))


ggplot(data = pcoa_FCgenus_BC, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull), alpha = 0.75, size = 3) +
  theme_test() + stat_ellipse(aes(color = Group), size = 2.5) +
  scale_color_manual(values = c("#873e23", "#18678d", "#626262"), labels = c("Cows", "Farmers", "Controls")) + xlab("PCoA 1 (55.0%)") + ylab("PCoA 2 (5.9%)") +
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"), axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 20)) +
  scale_shape_discrete(name = "Season", labels = c("Spring 2019", "Summer 2019", "Fall 2019", "Winter 2019", "Spring 2020"))
p <- ggplot(data = pcoa_FCgenus_BC, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull), alpha = 0.75, size = 4) +
  theme_test() + stat_ellipse(aes(color = Group), size = 2.5) +
  scale_color_manual(values = c("#873e23", "#18678d", "#626262"), labels = c("Cows", "Farmers", "Controls")) + xlab("PCoA 1 (55.0%)") + ylab("PCoA 2 (5.9%)") +
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"), axis.text = element_blank(), axis.title = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 20), legend.position = "none") +
  scale_shape_discrete(name = "Season", labels = c("Spring 2019", "Summer 2019", "Fall 2019", "Winter 2019", "Spring 2020"))
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/16S_Fecal_PCOA_v2.svg", p, bg = "transparent", width = 19.5, height = 9.9)

p <- ggplot(data = pcoa_FCgenus_BC[which(pcoa_FCgenus_BC$SeasonFull != "SPRING 2"),], aes(x = V1, y = V2)) + geom_point(aes(color = Group), alpha = 0.75, size = 6) +
  theme_test() + stat_ellipse(aes(color = Group), size = 2.5) +
  scale_color_manual(values = c("#873e23", "#18678d", "#626262"), labels = c("Cows", "Farmers", "Controls")) + xlab("PCoA 1 (55.0%)") + ylab("PCoA 2 (5.9%)") +
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"), axis.text = element_blank(), axis.title = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 20), legend.position = "none") +
  scale_shape_manual(values = c(15, 15, 15, 15), name = "Season", labels = c("Spring 2019", "Summer 2019", "Fall 2019", "Winter 2019"))
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/16S_Fecal_PCOA_v4.svg", p, bg = "transparent", width = 12, height = 10)

ns_Genus_BrayCurtis <- vegdist(genus_abund_ns, method = "bray")
pcoa_NSgenus_BC <- pco(ns_Genus_BrayCurtis, k = 4)
pcoa_NSgenus_scree <- as.data.frame(pcoa_NSgenus_BC$eig*100/sum(pcoa_NSgenus_BC$eig))
pcoa_NSgenus_BC <- as.data.frame(pcoa_NSgenus_BC$points)
pcoa_NSgenus_BC <- merge(pcoa_NSgenus_BC, samples_metadata[, c(1,2,3,12,13,16)], by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_NSgenus_BC)[1] <- "Sample"
pcoa_NSgenus_BC$SeasonFull <- paste(pcoa_NSgenus_BC$Season, pcoa_NSgenus_BC$Year)
pcoa_NSgenus_BC$SeasonFull <- factor(pcoa_NSgenus_BC$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))
pcoa_NSgenus_BC$Group <- factor(pcoa_NSgenus_BC$Group, levels = c("C", "W", "D"))
#pcoa_NSgenus_BC_MeanSeasonFull <- pcoa_NSgenus_BC[,c(9,11,2,3,4,5)] %>% group_by(Group, SeasonFull) %>% summarise_all(mean)

ggplot(data = pcoa_NSgenus_BC, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull), alpha = 0.75, size = 3) +
  theme_test() + stat_ellipse(aes(color = Group), size = 2.5) +
  scale_color_manual(values = c("#873e23", "#18678d", "#626262"), labels = c("Cows", "Farmers", "Controls")) + xlab("PCoA 1 (15.2%)") + ylab("PCoA 2 (10.6%)") +
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"), axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 20)) +
  scale_shape_discrete(name = "Season", labels = c("Spring 2019", "Summer 2019", "Fall 2019", "Winter 2019", "Spring 2020"))
p <- ggplot(data = pcoa_NSgenus_BC, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull), alpha = 0.75, size = 4) +
  theme_test() + stat_ellipse(aes(color = Group), size = 2.5) +
  scale_color_manual(values = c("#873e23", "#18678d", "#626262"), labels = c("Cows", "Farmers", "Controls")) + xlab("PCoA 1 (15.2%)") + ylab("PCoA 2 (10.6%)") +
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"), axis.text = element_blank(), axis.title = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 20), legend.position = "none") +
  scale_shape_discrete(name = "Season", labels = c("Spring 2019", "Summer 2019", "Fall 2019", "Winter 2019", "Spring 2020"))
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/16S_Nasal_PCOA_v3.svg", p, bg = "transparent", width = 14.34, height = 10)

p <- ggplot(data = pcoa_NSgenus_BC[which(pcoa_NSgenus_BC$SeasonFull != "SPRING 2"),], aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull), alpha = 0.75, size = 4) +
  theme_test() + stat_ellipse(aes(color = Group), size = 2.5) +
  scale_color_manual(values = c("#873e23", "#18678d", "#626262"), labels = c("Cows", "Farmers", "Controls")) + xlab("PCoA 1 (15.2%)") + ylab("PCoA 2 (10.6%)") +
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"), axis.text = element_blank(), axis.title = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 20), legend.position = "none") +
  scale_shape_manual(values = c(15, 15, 15, 15), name = "Season", labels = c("Spring 2019", "Summer 2019", "Fall 2019", "Winter 2019"))
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/16S_Nasal_PCOA_v6.svg", p, bg = "transparent", width = 10, height = 10)

ns_Genus_BrayCurtis_2019 <- vegdist(genus_abund_ns_2019, method = "bray")
pcoa_NSgenus_BC_2019 <- pco(ns_Genus_BrayCurtis_2019, k = 4)
pcoa_NSgenus_2019_scree <- as.data.frame(pcoa_NSgenus_BC_2019$eig*100/sum(pcoa_NSgenus_BC_2019$eig))
pcoa_NSgenus_BC_2019 <- as.data.frame(pcoa_NSgenus_BC_2019$points)
pcoa_NSgenus_BC_2019 <- merge(pcoa_NSgenus_BC_2019, samples_metadata[, c(1,2,3,12,13,16)], by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_NSgenus_BC_2019)[1] <- "Sample"
pcoa_NSgenus_BC_2019$SeasonFull <- paste(pcoa_NSgenus_BC_2019$Season, pcoa_NSgenus_BC_2019$Year)
pcoa_NSgenus_BC_2019$SeasonFull <- factor(pcoa_NSgenus_BC_2019$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1"))
pcoa_NSgenus_BC_2019$Group <- factor(pcoa_NSgenus_BC_2019$Group, levels = c("C", "W", "D"))

p <- ggplot(data = pcoa_NSgenus_BC_2019, aes(x = V1, y = V2)) + geom_point(aes(color = Group), alpha = 0.75, size = 6) +
  theme_test() + stat_ellipse(aes(color = Group), size = 2.5) +
  scale_color_manual(values = c("#873e23", "#18678d", "#626262"), labels = c("Cows", "Farmers", "Controls")) + xlab("PCoA 1 (16.1%)") + ylab("PCoA 2 (11.1%)") +
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"), axis.text = element_blank(), axis.title = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 20), legend.position = "none") +
  scale_shape_manual(values = c(15, 15, 15, 15), name = "Season", labels = c("Spring 2019", "Summer 2019", "Fall 2019", "Winter 2019"))
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/16S_Nasal_PCOA_v8.svg", p, bg = "transparent", width = 10, height = 10)

adonis_genus_abund_ns_2019 <- merge(genus_abund_ns_2019, samples_metadata[,c(3,13)], by.x = "row.names", by.y = "Sample_ID")
row.names(adonis_genus_abund_ns_2019) <- adonis_genus_abund_ns_2019$Row.names
adonis_genus_abund_ns_2019 <- adonis_genus_abund_ns_2019[,c(-1)]

adonis_genus_abund_ns_2019 <- adonis_genus_abund_ns_2019[which(adonis_genus_abund_ns_2019$Group != "D"),]

adonis2(adonis_genus_abund_ns_2019[, c(-1169)] ~ adonis_genus_abund_ns_2019$Group,
        permutations = 9999, method = "bray")

adonis_genus_abund_fc_2019 <- merge(genus_abund_fc_2019, samples_metadata[,c(3,13)], by.x = "row.names", by.y = "Sample_ID")
row.names(adonis_genus_abund_fc_2019) <- adonis_genus_abund_fc_2019$Row.names
adonis_genus_abund_fc_2019 <- adonis_genus_abund_fc_2019[,c(-1)]

adonis_genus_abund_fc_2019 <- adonis_genus_abund_fc_2019[which(adonis_genus_abund_fc_2019$Group != "C"),]

adonis2(adonis_genus_abund_fc_2019[, c(-322)] ~ adonis_genus_abund_fc_2019$Group,
        permutations = 9999, method = "bray")

###Looking at the human beta diversity
samples_metadata_human <- samples_metadata[which(samples_metadata$Group != "C" & samples_metadata$Sample_Type == "Fecal" & samples_metadata$Season == "WINTER" & samples_metadata$Year == "1"),]
human_genus_fc <- genus_abund_fc[which(rownames(genus_abund_fc) %in% samples_metadata_human$Sample_ID),]
human_genus_fc <- human_genus_fc[, colSums(human_genus_fc) > 0]
samples_metadata_human <- samples_metadata_human[which(samples_metadata_human$Sample_ID %in% row.names(human_genus_fc)),]
samples_metadata_human <- samples_metadata_human[order(samples_metadata_human$Sample_ID),]
human_genus_fc <- human_genus_fc[order(row.names(human_genus_fc)),]

fc_human_Genus_Bray_Curtis <- vegdist(human_genus_fc, method = "bray")
pcoa_human_FCGenus_BC <- pco(fc_human_Genus_Bray_Curtis, k = 4)
pcoa_human_FCGenus_BC_scree <- as.data.frame(100*pcoa_human_FCGenus_BC$eig/sum(pcoa_human_FCGenus_BC$eig))
pcoa_human_FCGenus_BC <- as.data.frame(pcoa_human_FCGenus_BC$points)
pcoa_human_FCGenus_BC <- merge(pcoa_human_FCGenus_BC, samples_metadata[, c(1,2,3,12,13,16)], by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_human_FCGenus_BC)[1] <- "Sample"
pcoa_human_FCGenus_BC$SeasonFull <- paste(pcoa_human_FCGenus_BC$Season, pcoa_human_FCGenus_BC$Year)
pcoa_human_FCGenus_BC$Group <- factor(pcoa_human_FCGenus_BC$Group, levels = c("W", "D"))

ggplot(data = pcoa_human_FCGenus_BC, aes(x = V1, y = V2)) + geom_point(aes(color = Group, fill = Group), pch = 7, size = 5) +
  theme_test() + labs(title = "Fecal samples PCOA. Bray-Curtis") + stat_ellipse(aes(color = Group))
p <- ggplot(data = pcoa_human_FCGenus_BC, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull), alpha = 1, size = 5, shape = 16) +
  theme_test() + stat_ellipse(aes(color = Group), size = 2.5, alpha = 0.75) +
  scale_color_manual(values = c("#18678d", "#626262"), labels = c("Farmers", "Controls")) + xlab("PCoA 1 (16.3%)") + ylab("PCoA 2 (11.4%)") +
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"), axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 20),
        legend.position = "none")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/16S_Fecal_Human_PCOA_Fall2019_v2.svg", p, bg = "transparent", width = 7, height = 7)

adonis2(human_genus_fc ~ samples_metadata_human$Group,
        permutations = 9999, method = "bray")

p.adjust(c(0.0897, 0.0589, 0.0691, 0.0154, 0.4978), method = "BH")
p.adjust(c(0.0897, 0.0589, 0.0691, 0.0154), method = "BH")

###Looking at pairwise BC distances
fc_BrayCurtis_pairwise <- as.matrix(fc_Genus_BrayCurtis)
fc_BrayCurtis_pairwise <- melt(fc_BrayCurtis_pairwise)
fc_BrayCurtis_pairwise <- fc_BrayCurtis_pairwise[!duplicated(apply(fc_BrayCurtis_pairwise, 1, function(x) paste(sort(x), collapse = ""))),]
fc_BrayCurtis_pairwise <- fc_BrayCurtis_pairwise[which(fc_BrayCurtis_pairwise$X1 != fc_BrayCurtis_pairwise$X2),]
fc_BrayCurtis_pairwise <- merge(fc_BrayCurtis_pairwise, samples_metadata[,c(1,3,12,13,15,16)], by.x = "X1", by.y = "Sample_ID")
colnames(fc_BrayCurtis_pairwise)[4] <- "SID1"
colnames(fc_BrayCurtis_pairwise)[5] <- "Season1"
colnames(fc_BrayCurtis_pairwise)[6] <- "Group1"
colnames(fc_BrayCurtis_pairwise)[7] <- "Site1"
colnames(fc_BrayCurtis_pairwise)[8] <- "Year1"
fc_BrayCurtis_pairwise <- merge(fc_BrayCurtis_pairwise, samples_metadata[,c(1,3,12,13,15,16)], by.x = "X2", by.y = "Sample_ID")
colnames(fc_BrayCurtis_pairwise)[9] <- "SID2"
colnames(fc_BrayCurtis_pairwise)[10] <- "Season2"
colnames(fc_BrayCurtis_pairwise)[11] <- "Group2"
colnames(fc_BrayCurtis_pairwise)[12] <- "Site2"
colnames(fc_BrayCurtis_pairwise)[13] <- "Year2"

fc_BrayCurtis_pairwise <- fc_BrayCurtis_pairwise[which(fc_BrayCurtis_pairwise$Season1 == fc_BrayCurtis_pairwise$Season2 & fc_BrayCurtis_pairwise$Year1 == fc_BrayCurtis_pairwise$Year2),]
fc_BrayCurtis_pairwise <- fc_BrayCurtis_pairwise[which(fc_BrayCurtis_pairwise$Group1 != fc_BrayCurtis_pairwise$Group2),]
fc_BrayCurtis_pairwise$Groups <- ""

for (i in seq(1,dim(fc_BrayCurtis_pairwise)[1])) {
  fc_BrayCurtis_pairwise$Groups[i] = paste(str_sort(c(fc_BrayCurtis_pairwise$Group1[i], fc_BrayCurtis_pairwise$Group2[i]))[1], str_sort(c(fc_BrayCurtis_pairwise$Group1[i], fc_BrayCurtis_pairwise$Group2[i]))[2])
}

fc_BrayCurtis_pairwise <- fc_BrayCurtis_pairwise[which(fc_BrayCurtis_pairwise$Groups != "D W"),]
fc_BrayCurtis_pairwise$Sites[fc_BrayCurtis_pairwise$Site1 == fc_BrayCurtis_pairwise$Site2] <- "Same site"
fc_BrayCurtis_pairwise$Sites[fc_BrayCurtis_pairwise$Site1 != fc_BrayCurtis_pairwise$Site2] <- "Different sites"
fc_BrayCurtis_pairwise$SeasonFull <- paste(fc_BrayCurtis_pairwise$Season1, fc_BrayCurtis_pairwise$Year1)
fc_BrayCurtis_pairwise$SeasonFull <- factor(fc_BrayCurtis_pairwise$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))

fc_BrayCurtis_pairwise$HumanGroup <- ""

for (i in seq(1,dim(fc_BrayCurtis_pairwise)[1])) {
  fc_BrayCurtis_pairwise$HumanGroup[i] = str_sort(c(fc_BrayCurtis_pairwise$Group1[i], fc_BrayCurtis_pairwise$Group2[i]))[2]
}

fc_BrayCurtis_pairwise$HumanSID <- ""

for (i in seq(1,dim(fc_BrayCurtis_pairwise)[1])) {
  fc_BrayCurtis_pairwise$HumanSID[i] = str_sort(c(fc_BrayCurtis_pairwise$SID1[i], fc_BrayCurtis_pairwise$SID2[i]))[2]
}

fc_BrayCurtis_pairwise <- fc_BrayCurtis_pairwise %>%
  group_by(HumanSID, Sites, SeasonFull, HumanGroup) %>%
  summarise(
    value = mean(value)
  )

fc_BrayCurtis_pairwise$GroupSites <- paste(fc_BrayCurtis_pairwise$HumanGroup, fc_BrayCurtis_pairwise$Sites)
fc_BrayCurtis_pairwise$HumanGroup <- factor(fc_BrayCurtis_pairwise$HumanGroup, levels = c("W", "D"))

fc_BC_pairwise_compared <- compare_means(value ~ GroupSites, data = fc_BrayCurtis_pairwise[which(fc_BrayCurtis_pairwise$SeasonFull != "SPRING 2"),], group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")
fc_BC_pairwise_compared <- fc_BC_pairwise_compared[which(fc_BC_pairwise_compared$p.adj < 0.05),]
compare_means(value ~ GroupSites, data = fc_BrayCurtis_pairwise[which(fc_BrayCurtis_pairwise$SeasonFull != "SPRING 2"),], method = "wilcox.test", p.adjust.method = "BH")

ggplot(data = fc_BrayCurtis_pairwise, aes(x = Sites, y = value, color = HumanGroup)) + geom_boxplot() + 
  facet_wrap(~fc_BrayCurtis_pairwise$SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_test() + ylab("Bray-Curtis distance") +
  ylim(0,1.2) +
  labs(title = "Fecal genus bray-Curtis distances to cow samples", subtitle = "Wilcoxon rank sum test. Benjamini-Hochberg correction.")
ggplot(data = fc_BrayCurtis_pairwise, aes(x = Sites, y = value)) + geom_boxplot(aes(fill = HumanGroup), color = "black", outlier.shape = NA, width = 0.3, size = 1, alpha = 0.5, position = position_dodge(width = 0.95)) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.95, jitter.height = 0, jitter.width = 0.25), alpha = 0.5, aes(group = HumanGroup, color = HumanGroup)) +
  geom_flat_violin(aes(fill = HumanGroup), position = position_nudgedodge(x = 0.12, width = 0.95), alpha = 0.7, scale = "width", width = 1.6) +
  facet_wrap(~fc_BrayCurtis_pairwise$SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_linedraw() +
  scale_fill_manual(values=c("#18678d", "#626262"), name = "", labels = c("Farmers", "Controls")) +
  scale_color_manual(values=c("#18678d", "#626262"), name = "", labels = c("Farmers", "Controls")) +
  scale_x_discrete(labels = c("Different\nsites", "Same\nsite")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_text(), axis.text.x = element_text(), strip.text = element_text(face = "bold", size =20),
        legend.position = "right", axis.title.y = element_blank(), axis.title.x = element_blank())

p <- ggplot(data = fc_BrayCurtis_pairwise, aes(x = Sites, y = value)) + geom_boxplot(aes(fill = HumanGroup), color = "black", outlier.shape = NA, width = 0.3, size = 1, alpha = 0.5, position = position_dodge(width = 0.95)) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.95, jitter.height = 0, jitter.width = 0.25), alpha = 0.5, aes(group = HumanGroup, color = HumanGroup)) +
  geom_flat_violin(aes(fill = HumanGroup), position = position_nudgedodge(x = 0.12, width = 0.95), alpha = 0.7, scale = "width", width = 1.6) +
  facet_wrap(~fc_BrayCurtis_pairwise$SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_linedraw() +
  scale_fill_manual(values=c("#18678d", "#626262"), name = "", labels = c("Farmers", "Controls")) +
  scale_color_manual(values=c("#18678d", "#626262"), name = "", labels = c("Farmers", "Controls")) +
  scale_x_discrete(labels = c("Different\nsites", "Same\nsite")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "bold", size =20),
        legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/16S_Fecal_PairwiseBC_v1.svg", p, bg = "transparent", width = 15, height = 10)

ggplot(data = fc_BrayCurtis_pairwise[which(fc_BrayCurtis_pairwise$SeasonFull != "SPRING 2"),], aes(x = Sites, y = value)) + geom_boxplot(aes(fill = HumanGroup), color = "black", outlier.shape = NA, width = 0.3, size = 1, alpha = 0.5, position = position_dodge(width = 0.95)) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.95, jitter.height = 0, jitter.width = 0.25), alpha = 0.5, aes(group = HumanGroup, color = HumanGroup)) +
  geom_flat_violin(aes(fill = HumanGroup), position = position_nudgedodge(x = 0.12, width = 0.95), alpha = 0.7, scale = "width", width = 1.6) +
  theme_linedraw() +
  scale_fill_manual(values=c("#18678d", "#626262"), name = "", labels = c("Farmers", "Controls")) +
  scale_color_manual(values=c("#18678d", "#626262"), name = "", labels = c("Farmers", "Controls")) +
  scale_x_discrete(labels = c("Different\nsites", "Same\nsite")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "bold", size =20),
        legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank())

fc_BrayCurtis_pairwise_mean <- fc_BrayCurtis_pairwise[fc_BrayCurtis_pairwise$SeasonFull != "SPRING 2", c(-3, -6)] %>%
  group_by(HumanSID, Sites, HumanGroup) %>%
  summarise(
    value = mean(value)
  )

fc_BrayCurtis_pairwise_mean$GroupSites <- paste(fc_BrayCurtis_pairwise_mean$HumanGroup, fc_BrayCurtis_pairwise_mean$Sites)
fc_BrayCurtis_pairwise_mean$HumanGroup <- factor(fc_BrayCurtis_pairwise_mean$HumanGroup, levels = c("W", "D"))

fc_BC_pairwise_compared <- compare_means(value ~ GroupSites, data = fc_BrayCurtis_pairwise_mean[,c(4,5)], method = "wilcox.test", p.adjust.method = "BH")
fc_BC_pairwise_compared <- fc_BC_pairwise_compared[which(fc_BC_pairwise_compared$p.adj < 0.05),]

p <- ggplot(data = fc_BrayCurtis_pairwise_mean, aes(x = Sites, y = value)) + geom_boxplot(aes(fill = HumanGroup), color = "black", outlier.shape = NA, width = 0.5, size = 1, alpha = 0.5, position = position_dodge(width = 0.9)) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.9, jitter.height = 0, jitter.width = 0.45), alpha = 0.5, aes(group = HumanGroup, color = HumanGroup)) +
  geom_flat_violin(aes(fill = HumanGroup), position = position_nudgedodge(x = 0.15, width = 0.9), alpha = 0.7, scale = "width", width = 1.2) +
  theme_linedraw() +
  scale_fill_manual(values=c("#18678d", "#626262"), name = "", labels = c("Farmers", "Controls")) +
  scale_color_manual(values=c("#18678d", "#626262"), name = "", labels = c("Farmers", "Controls")) +
  scale_x_discrete(labels = c("Different\nsites", "Same\nsite")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text = element_blank(),
        legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/16S_Fecal_PairwiseBC_v4.svg", p, bg = "transparent", width = 6, height = 10)


ns_BrayCurtis_pairwise <- as.matrix(ns_Genus_BrayCurtis)
ns_BrayCurtis_pairwise <- melt(ns_BrayCurtis_pairwise)
ns_BrayCurtis_pairwise <- ns_BrayCurtis_pairwise[!duplicated(apply(ns_BrayCurtis_pairwise, 1, function(x) paste(sort(x), collapse = ""))),]
ns_BrayCurtis_pairwise <- ns_BrayCurtis_pairwise[which(ns_BrayCurtis_pairwise$X1 != ns_BrayCurtis_pairwise$X2),]
ns_BrayCurtis_pairwise <- merge(ns_BrayCurtis_pairwise, samples_metadata[,c(1,3,12,13,15,16)], by.x = "X1", by.y = "Sample_ID")
colnames(ns_BrayCurtis_pairwise)[4] <- "SID1"
colnames(ns_BrayCurtis_pairwise)[5] <- "Season1"
colnames(ns_BrayCurtis_pairwise)[6] <- "Group1"
colnames(ns_BrayCurtis_pairwise)[7] <- "Site1"
colnames(ns_BrayCurtis_pairwise)[8] <- "Year1"
ns_BrayCurtis_pairwise <- merge(ns_BrayCurtis_pairwise, samples_metadata[,c(1,3,12,13,15,16)], by.x = "X2", by.y = "Sample_ID")
colnames(ns_BrayCurtis_pairwise)[9] <- "SID2"
colnames(ns_BrayCurtis_pairwise)[10] <- "Season2"
colnames(ns_BrayCurtis_pairwise)[11] <- "Group2"
colnames(ns_BrayCurtis_pairwise)[12] <- "Site2"
colnames(ns_BrayCurtis_pairwise)[13] <- "Year2"

ns_BrayCurtis_pairwise <- ns_BrayCurtis_pairwise[which(ns_BrayCurtis_pairwise$Season1 == ns_BrayCurtis_pairwise$Season2 & ns_BrayCurtis_pairwise$Year1 == ns_BrayCurtis_pairwise$Year2),]
ns_BrayCurtis_pairwise <- ns_BrayCurtis_pairwise[which(ns_BrayCurtis_pairwise$Group1 != ns_BrayCurtis_pairwise$Group2),]
ns_BrayCurtis_pairwise$Groups <- ""

for (i in seq(1,dim(ns_BrayCurtis_pairwise)[1])) {
  ns_BrayCurtis_pairwise$Groups[i] = paste(str_sort(c(ns_BrayCurtis_pairwise$Group1[i], ns_BrayCurtis_pairwise$Group2[i]))[1], str_sort(c(ns_BrayCurtis_pairwise$Group1[i], ns_BrayCurtis_pairwise$Group2[i]))[2])
}

ns_BrayCurtis_pairwise <- ns_BrayCurtis_pairwise[which(ns_BrayCurtis_pairwise$Groups != "D W"),]
ns_BrayCurtis_pairwise$Sites[ns_BrayCurtis_pairwise$Site1 == ns_BrayCurtis_pairwise$Site2] <- "Same site"
ns_BrayCurtis_pairwise$Sites[ns_BrayCurtis_pairwise$Site1 != ns_BrayCurtis_pairwise$Site2] <- "Different sites"
ns_BrayCurtis_pairwise$SeasonFull <- paste(ns_BrayCurtis_pairwise$Season1, ns_BrayCurtis_pairwise$Year1)
ns_BrayCurtis_pairwise$SeasonFull <- factor(ns_BrayCurtis_pairwise$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))

ns_BrayCurtis_pairwise$HumanGroup <- ""

for (i in seq(1,dim(ns_BrayCurtis_pairwise)[1])) {
  ns_BrayCurtis_pairwise$HumanGroup[i] = str_sort(c(ns_BrayCurtis_pairwise$Group1[i], ns_BrayCurtis_pairwise$Group2[i]))[2]
}

ns_BrayCurtis_pairwise$HumanSID <- ""

for (i in seq(1,dim(ns_BrayCurtis_pairwise)[1])) {
  ns_BrayCurtis_pairwise$HumanSID[i] = str_sort(c(ns_BrayCurtis_pairwise$SID1[i], ns_BrayCurtis_pairwise$SID2[i]))[2]
}

ns_BrayCurtis_pairwise <- ns_BrayCurtis_pairwise %>%
  group_by(HumanSID, Sites, SeasonFull, HumanGroup) %>%
  summarise(
    value = mean(value)
  )

ns_BrayCurtis_pairwise$GroupSites <- paste(ns_BrayCurtis_pairwise$HumanGroup, ns_BrayCurtis_pairwise$Sites)
ns_BrayCurtis_pairwise$HumanGroup <- factor(ns_BrayCurtis_pairwise$HumanGroup, levels = c("W", "D"))

ns_BC_pairwise_compared <- compare_means(value ~ GroupSites, data = ns_BrayCurtis_pairwise[which(ns_BrayCurtis_pairwise$SeasonFull != "SPRING 2"),c(3,5,6)], group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")
ns_BC_pairwise_compared <- ns_BC_pairwise_compared[which(ns_BC_pairwise_compared$p.adj < 0.05),]
compare_means(value ~ GroupSites, data = ns_BrayCurtis_pairwise[which(ns_BrayCurtis_pairwise$SeasonFull != "SPRING 2"),c(3,5,6)],method = "wilcox.test", p.adjust.method = "BH")

ggplot(data = ns_BrayCurtis_pairwise, aes(x = Sites, y = value, color = HumanGroup)) + geom_boxplot() + 
  facet_wrap(~ns_BrayCurtis_pairwise$SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_test() + ylab("Bray-Curtis distance") +
  ylim(0,1.2) +
  labs(title = "Nasal genus bray-Curtis distances to cow samples", subtitle = "Wilcoxon rank sum test. Benjamini-Hochberg correction.") +
  scale_color_manual(values=c("#00BA38", "#619CFF"), name = "", labels = c("Office worker", "Farmer"))

p <- ggplot(data = ns_BrayCurtis_pairwise, aes(x = Sites, y = value)) + geom_boxplot(aes(fill = HumanGroup), color = "black", outlier.shape = NA, width = 0.3, size = 1, alpha = 0.5, position = position_dodge(width = 0.95)) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.95, jitter.height = 0, jitter.width = 0.25), alpha = 0.5, aes(group = HumanGroup, color = HumanGroup)) +
  geom_flat_violin(aes(fill = HumanGroup), position = position_nudgedodge(x = 0.12, width = 0.95), alpha = 0.7, scale = "width", width = 1.6) +
  facet_wrap(~ns_BrayCurtis_pairwise$SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_linedraw() +
  ylim(0.25,1.2) +
  scale_fill_manual(values=c("#18678d", "#626262"), name = "", labels = c("Farmers", "Controls")) +
  scale_color_manual(values=c("#18678d", "#626262"), name = "", labels = c("Farmers", "Controls")) +
  scale_x_discrete(labels = c("Different\nsites", "Same\nsite")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "bold", size =20),
        legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/16S_Nasal_PairwiseBC_v2.svg", p, bg = "transparent", width = 15, height = 10)

p <- ggplot(data = ns_BrayCurtis_pairwise[which(ns_BrayCurtis_pairwise$SeasonFull != "SPRING 2"),], aes(x = Sites, y = value)) + geom_boxplot(aes(fill = HumanGroup), color = "black", outlier.shape = NA, width = 0.5, size = 1, alpha = 0.5, position = position_dodge(width = 0.9)) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.9, jitter.height = 0, jitter.width = 0.45), alpha = 0.5, aes(group = HumanGroup, color = HumanGroup)) +
  geom_flat_violin(aes(fill = HumanGroup), position = position_nudgedodge(x = 0.15, width = 0.9), alpha = 0.7, scale = "width", width = 1.2) +
  theme_linedraw() +
  ylim(0.25,1.2) +
  scale_fill_manual(values=c("#18678d", "#626262"), name = "", labels = c("Farmers", "Controls")) +
  scale_color_manual(values=c("#18678d", "#626262"), name = "", labels = c("Farmers", "Controls")) +
  scale_x_discrete(labels = c("Different\nsites", "Same\nsite")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "bold", size =20),
        legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/16S_Nasal_PairwiseBC_v3.svg", p, bg = "transparent", width = 5, height = 10)

ns_BrayCurtis_pairwise_mean <- ns_BrayCurtis_pairwise[which(ns_BrayCurtis_pairwise$SeasonFull != "SPRING 2"), c(-3, -6)] %>%
  group_by(HumanSID, Sites, HumanGroup) %>%
  summarise(
    value = mean(value)
  )

ns_BrayCurtis_pairwise_mean$GroupSites <- paste(ns_BrayCurtis_pairwise_mean$HumanGroup, ns_BrayCurtis_pairwise_mean$Sites)
ns_BrayCurtis_pairwise_mean$HumanGroup <- factor(ns_BrayCurtis_pairwise_mean$HumanGroup, levels = c("W", "D"))

ns_BC_pairwise_compared <- compare_means(value ~ GroupSites, data = ns_BrayCurtis_pairwise_mean[,c(4,5)], method = "wilcox.test", p.adjust.method = "BH")
ns_BC_pairwise_compared <- ns_BC_pairwise_compared[which(ns_BC_pairwise_compared$p.adj < 0.05),]

p <- ggplot(data = ns_BrayCurtis_pairwise_mean, aes(x = Sites, y = value)) + geom_boxplot(aes(fill = HumanGroup), color = "black", outlier.shape = NA, width = 0.5, size = 1, alpha = 0.5, position = position_dodge(width = 0.9)) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.9, jitter.height = 0, jitter.width = 0.45), alpha = 0.5, aes(group = HumanGroup, color = HumanGroup)) +
  geom_flat_violin(aes(fill = HumanGroup), position = position_nudgedodge(x = 0.15, width = 0.9), alpha = 0.7, scale = "width", width = 1.2) +
  theme_linedraw() +
  ylim(0.25,1.2) +
  scale_fill_manual(values=c("#18678d", "#626262"), name = "", labels = c("Farmers", "Controls")) +
  scale_color_manual(values=c("#18678d", "#626262"), name = "", labels = c("Farmers", "Controls")) +
  scale_x_discrete(labels = c("Different\nsites", "Same\nsite")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text = element_blank(),
        legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/16S_Nasal_PairwiseBC_v6.svg", p, bg = "transparent", width = 6, height = 10)

###Looking at the total relative abundances in the gut of the genera enriched in the farmer noses
genus_abund_fc_FarmNasalEnrichedGenus <- as.data.frame(rowSums(genus_abund_fc[colnames(genus_abund_fc) %in% M2_ns_genus_famerVSoffice[which(M2_ns_genus_famerVSoffice$coef > 0),"feature"]]))
colnames(genus_abund_fc_FarmNasalEnrichedGenus)[1] <- "TotRelAbund"
genus_abund_fc_FarmNasalEnrichedGenus <- merge(genus_abund_fc_FarmNasalEnrichedGenus, samples_metadata[,c(1,3,12,13,15,16)], by.x = "row.names", by.y = "Sample_ID")
ggplot(data = genus_abund_fc_FarmNasalEnrichedGenus[which(genus_abund_fc_FarmNasalEnrichedGenus$Group != "C" & genus_abund_fc_FarmNasalEnrichedGenus$Season == "SUMMER" & genus_abund_fc_FarmNasalEnrichedGenus$Year == "1"),]) + 
  geom_boxplot(aes(x = Group, y = TotRelAbund))

###Running Maaslin2 to check for significant associations of genera and families with farmers relative to office workers
metadata_maaslin <- samples_metadata
rownames(metadata_maaslin) <- metadata_maaslin$Sample_ID
metadata_maaslin$SeasonFull <- paste(metadata_maaslin$Season, metadata_maaslin$Year)

Maaslin2(
  input_data = genus_abund_fc,
  input_metadata = metadata_maaslin[which(metadata_maaslin$Sample_Type == "Fecal" & metadata_maaslin$Group != "C"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/fc_Genus_farmerVSoffice",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)

Maaslin2(
  input_data = genus_abund_fc,
  input_metadata = metadata_maaslin[which(metadata_maaslin$Sample_Type == "Fecal" & metadata_maaslin$Group != "C" & metadata_maaslin$SeasonFull == "SPRING 2"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/fc_Genus_farmerVSoffice_SPRING2",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

Maaslin2(
  input_data = genus_abund_fc[colnames(genus_abund_fc) %in% M2_ns_genus_famerVSoffice[which(M2_ns_genus_famerVSoffice$coef > 0),"feature"]],
  input_metadata = metadata_maaslin[which(metadata_maaslin$Sample_Type == "Fecal" & metadata_maaslin$Group != "C"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/fc_FarmerNasalEnrichedGenus_farmerVSoffice",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)

Maaslin2(
  input_data = genus_abund_fc[colnames(genus_abund_fc) %in% M2_fc_genus_cowVSoffice[which(M2_fc_genus_cowVSoffice$coef > 0),"feature"]],
  input_metadata = metadata_maaslin[which(metadata_maaslin$Sample_Type == "Fecal" & metadata_maaslin$Group != "C"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/fc_CowGutlEnrichedGenus_farmerVSoffice",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)

Maaslin2(
  input_data = genus_abund_fc[colnames(genus_abund_fc) %in% M2_fc_genus_cowVSoffice[which(M2_fc_genus_cowVSoffice$coef > 0),"feature"]],
  input_metadata = metadata_maaslin[which(metadata_maaslin$Sample_Type == "Fecal" & metadata_maaslin$Group != "C" & metadata_maaslin$SeasonFull == "SPRING 1"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/fc_CowGutlEnrichedGenus_farmerVSoffice_SPRING1",
  fixed_effects = c("Group"),
  random_effects = c("Site")
)

Maaslin2(
  input_data = genus_abund_fc[colnames(genus_abund_fc) %in% M2_ns_genus_famerVSoffice[which(M2_ns_genus_famerVSoffice$coef > 0),"feature"]],
  input_metadata = metadata_maaslin[which(metadata_maaslin$Sample_Type == "Fecal" & metadata_maaslin$Group != "C" & metadata_maaslin$SeasonFull == "SPRING 1"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/fc_FarmerNasalEnrichedGenus_farmerVSoffice_SPRING1",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)


Maaslin2(
  input_data = family_abund_fc,
  input_metadata = metadata_maaslin[which(metadata_maaslin$Sample_Type == "Fecal" & metadata_maaslin$Group != "C"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/fc_Family_farmerVSoffice",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)

M2_fc_family_famerVSoffice <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/fc_Family_farmerVSoffice/significant_results.tsv", header = TRUE, "\t")
M2_fc_family_famerVSoffice$shape <- ifelse(M2_fc_family_famerVSoffice$coef > 0, 1, -1)
M2_fc_family_famerVSoffice <- subset(M2_fc_family_famerVSoffice, qval < 0.05)
ggplot(M2_fc_family_famerVSoffice, aes(x=coef, y=reorder(feature, coef), color=shape))+
  geom_point(size=2)+
  labs(x="Relative to Office Worker Gut Microbiome",
       y="Family", shape="Shape")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_vline(xintercept=0, colour='black',size=1, linetype="dotted")+
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0)) + 
  scale_colour_gradient2() + labs(title = "Maaslin2 on gut microbiota")

Maaslin2(
  input_data = genus_abund_fc,
  input_metadata = metadata_maaslin[which(metadata_maaslin$Sample_Type == "Fecal" & metadata_maaslin$Group != "W"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/fc_Genus_cowVSoffice",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)

M2_fc_genus_cowVSoffice <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/fc_Genus_cowVSoffice/significant_results.tsv", header = TRUE, "\t")
M2_fc_genus_cowVSoffice$coef = -1 * M2_fc_genus_cowVSoffice$coef
M2_fc_genus_cowVSoffice$value <- "C"
M2_fc_genus_cowVSoffice$shape <- ifelse(M2_fc_genus_cowVSoffice$coef > 0, 1, -1)
M2_fc_genus_cowVSoffice <- subset(M2_fc_genus_cowVSoffice, qval < 0.05)
ggplot(M2_fc_genus_cowVSoffice, aes(x=coef, y=reorder(feature, coef), color=shape))+
  geom_point(size=2)+
  labs(x="Relative to Office Worker Gut Microbiome",
       y="Genus", shape="Shape")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_vline(xintercept=0, colour='black',size=1, linetype="dotted")+
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0)) + 
  scale_colour_gradient2() + labs(title = "Maaslin2 on gut microbiota of cows and office workers")
p <- ggplot(M2_fc_genus_cowVSoffice, aes(x=coef, y=reorder(feature, coef)))+
  geom_point(color = "#873e23", size = 5, alpha = 0.7)+
  geom_vline(xintercept=0, colour='black',size=1.5, linetype="dotted")+
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0), color = "#873e23", size = 1.5, alpha = 1) + 
  theme_linedraw()+
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 14)) +
  labs(x="Coefficient (relative to controls)", y="Family")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/16S_Fecal_Maaslin2_Cow-vs-Office_Genus_v4.svg", p, bg = "transparent", width = 10, height = 30)


Maaslin2(
  input_data = family_abund_fc,
  input_metadata = metadata_maaslin[which(metadata_maaslin$Sample_Type == "Fecal" & metadata_maaslin$Group != "W"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/fc_Family_cowVSoffice",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)

M2_fc_family_cowVSoffice <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/fc_Family_cowVSoffice/significant_results.tsv", header = TRUE, "\t")
M2_fc_family_cowVSoffice$coef = -1 * M2_fc_family_cowVSoffice$coef
M2_fc_family_cowVSoffice$value <- "C"
M2_fc_family_cowVSoffice$shape <- ifelse(M2_fc_family_cowVSoffice$coef > 0, 1, -1)
M2_fc_family_cowVSoffice <- subset(M2_fc_family_cowVSoffice, qval < 0.05)
ggplot(M2_fc_family_cowVSoffice, aes(x=coef, y=reorder(feature, coef), color=shape))+
  geom_point(size=2)+
  labs(x="Relative to Office Worker Gut Microbiome",
       y="Family", shape="Shape")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_vline(xintercept=0, colour='black',size=1, linetype="dotted")+
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0)) + 
  scale_colour_gradient2() + labs(title = "Maaslin2 on gut microbiota of cows and office workers")
p <- ggplot(M2_fc_family_cowVSoffice, aes(x=coef, y=reorder(feature, coef)))+
  geom_point(color = "#873e23", size = 5, alpha = 0.7)+
  geom_vline(xintercept=0, colour='black',size=1.5, linetype="dotted")+
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0), color = "#873e23", size = 1.5, alpha = 1) + 
  theme_linedraw()+
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 14)) +
  labs(x="Coefficient (relative to controls)", y="Family")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/16S_Fecal_Maaslin2_Cow-vs-Office_Family_v4.svg", p, bg = "transparent", width = 10, height = 30)


Maaslin2(
  input_data = genus_abund_ns,
  input_metadata = metadata_maaslin[which(metadata_maaslin$Sample_Type == "Nasal" & metadata_maaslin$Group != "C"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/ns_Genus_farmerVSoffice",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)

M2_ns_genus_famerVSoffice <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/ns_Genus_farmerVSoffice/significant_results.tsv", header = TRUE, "\t")
M2_ns_genus_famerVSoffice$shape <- ifelse(M2_ns_genus_famerVSoffice$coef > 0, 1, -1)
M2_ns_genus_famerVSoffice <- subset(M2_ns_genus_famerVSoffice, qval < 0.05)
ggplot(M2_ns_genus_famerVSoffice, aes(x=coef, y=reorder(feature, coef), color=shape))+
  geom_point(size=2)+
  labs(x="Relative to Office Worker Nasal Microbiome",
       y="Genus", shape="Shape")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_vline(xintercept=0, colour='black',size=1, linetype="dotted")+
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0)) + 
  scale_colour_gradient2() + labs(title = "Maaslin2 on nasal microbiota of farmers and office workers")

dim(M2_ns_genus_famerVSoffice[which(M2_ns_genus_famerVSoffice$coef > 0),])

Maaslin2(
  input_data = genus_abund_ns,
  input_metadata = metadata_maaslin[which(metadata_maaslin$Sample_Type == "Nasal" & metadata_maaslin$Group != "W"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/ns_Genus_cowVSoffice",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)

M2_ns_genus_cowVSoffice <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/ns_Genus_cowVSoffice/significant_results.tsv", header = TRUE, "\t")
M2_ns_genus_cowVSoffice$coef = -1 * M2_ns_genus_cowVSoffice$coef
M2_ns_genus_cowVSoffice$value <- "C"
M2_ns_genus_cowVSoffice$shape <- ifelse(M2_ns_genus_cowVSoffice$coef > 0, 1, -1)
M2_ns_genus_cowVSoffice <- subset(M2_ns_genus_cowVSoffice, qval < 0.05)
ggplot(M2_ns_genus_cowVSoffice, aes(x=coef, y=reorder(feature, coef), color=shape))+
  geom_point(size=2)+
  labs(x="Relative to Office Worker Nasal Microbiome",
       y="Genus", shape="Shape")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_vline(xintercept=0, colour='black',size=1, linetype="dotted")+
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0)) + 
  scale_colour_gradient2() + labs(title = "Maaslin2 on nasal microbiota of cows and office workers")
p <- ggplot(M2_ns_genus_cowVSoffice, aes(x=coef, y=reorder(feature, coef)))+
  geom_point(color = "#873e23", size = 5, alpha = 0.7)+
  geom_vline(xintercept=0, colour='black',size=1.5, linetype="dotted")+
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0), color = "#873e23", size = 1.5, alpha = 1) + 
  theme_linedraw()+
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 14)) +
  labs(x="Coefficient (relative to controls)", y="Genus")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/16S_Nasal_Maaslin2_Cow-vs-Office_Genus_v4.svg", p, bg = "transparent", width = 10, height = 30)

M2_ns_genus_farmerANDcowVSoffice <- M2_ns_genus_famerVSoffice[,c("feature", "value", "coef", "stderr")]
colnames(M2_ns_genus_farmerANDcowVSoffice)[1] <- "Genus"
colnames(M2_ns_genus_farmerANDcowVSoffice)[2] <- "Group_W"
colnames(M2_ns_genus_farmerANDcowVSoffice)[3] <- "Coef_W"
colnames(M2_ns_genus_farmerANDcowVSoffice)[4] <- "StdErr_W"

M2_ns_genus_farmerANDcowVSoffice <- merge(M2_ns_genus_farmerANDcowVSoffice, M2_ns_genus_cowVSoffice[which(M2_ns_genus_cowVSoffice$feature %in% M2_ns_genus_farmerANDcowVSoffice$Genus),c("feature", "value", "coef", "stderr")], by.x = "Genus", by.y = "feature", all = TRUE)
colnames(M2_ns_genus_farmerANDcowVSoffice)[5] <- "Group_C"
colnames(M2_ns_genus_farmerANDcowVSoffice)[6] <- "Coef_C"
colnames(M2_ns_genus_farmerANDcowVSoffice)[7] <- "StdErr_C"

M2_ns_genus_farmerANDcowVSoffice$MeanCoef <- rowMeans(M2_ns_genus_farmerANDcowVSoffice[,c("Coef_W", "Coef_C")], na.rm = TRUE)

ggplot() + 
  geom_point(data = M2_ns_genus_farmerANDcowVSoffice, aes(x = Coef_W, y = reorder(Genus, MeanCoef)), color = "#18678d", shape = 16, size = 6, alpha = 0.8) +
  geom_errorbarh(data = M2_ns_genus_farmerANDcowVSoffice, aes(xmax = Coef_W + StdErr_W, xmin = Coef_W - StdErr_W, y = Genus, height = 0), color = "#18678d", size = 1.8, alpha = 0.7) + 
  geom_point(data = M2_ns_genus_farmerANDcowVSoffice[which(!is.na(M2_ns_genus_farmerANDcowVSoffice$Group_C)),], aes(x = Coef_C, y = Genus), color = "#873e23", size = 6, alpha = 0.8) + 
  geom_errorbarh(data = M2_ns_genus_farmerANDcowVSoffice[which(!is.na(M2_ns_genus_farmerANDcowVSoffice$Group_C)),], aes(xmax = Coef_C + StdErr_C, xmin = Coef_C - StdErr_C, y = Genus, height = 0), color = "#873e23", size = 1.8, alpha = 0.7) +
  theme_linedraw()+
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 20, face = "bold"), axis.text.y = element_text(size = 14)) +
  geom_vline(xintercept=0, colour='black',size=1, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Genus")
p <- ggplot() + 
  geom_point(data = M2_ns_genus_farmerANDcowVSoffice, aes(x = Coef_W, y = reorder(Genus, MeanCoef)), color = "#18678d", shape = 16, size = 6, alpha = 0.8) +
  geom_errorbarh(data = M2_ns_genus_farmerANDcowVSoffice, aes(xmax = Coef_W + StdErr_W, xmin = Coef_W - StdErr_W, y = Genus, height = 0), color = "#18678d", size = 1.8, alpha = 0.7) + 
  geom_point(data = M2_ns_genus_farmerANDcowVSoffice[which(!is.na(M2_ns_genus_farmerANDcowVSoffice$Group_C)),], aes(x = Coef_C, y = Genus), color = "#873e23", size = 6, alpha = 0.8) + 
  geom_errorbarh(data = M2_ns_genus_farmerANDcowVSoffice[which(!is.na(M2_ns_genus_farmerANDcowVSoffice$Group_C)),], aes(xmax = Coef_C + StdErr_C, xmin = Coef_C - StdErr_C, y = Genus, height = 0), color = "#873e23", size = 1.8, alpha = 0.7) +
  theme_linedraw()+
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 14)) +
  geom_vline(xintercept=0, colour='black',size=1, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Genus")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/16S_Nasal_Maaslin2_FarmerCow_Genus_v2.svg", p, bg = "transparent", width = 15, height = 20)

dim(M2_ns_genus_farmerANDcowVSoffice[which(M2_ns_genus_farmerANDcowVSoffice$Coef_C >0),])

ggplot(data = M2_ns_genus_farmerANDcowVSoffice[which(M2_ns_genus_farmerANDcowVSoffice$Coef_C >0),], aes(x = Coef_C, y = Coef_W)) + geom_point()

cor.test(x = M2_ns_genus_farmerANDcowVSoffice[which(M2_ns_genus_farmerANDcowVSoffice$Coef_C >0), c("Coef_C")], y = M2_ns_genus_farmerANDcowVSoffice[which(M2_ns_genus_farmerANDcowVSoffice$Coef_C >0), c("Coef_W")])

Maaslin2(
  input_data = family_abund_ns,
  input_metadata = metadata_maaslin[which(metadata_maaslin$Sample_Type == "Nasal" & metadata_maaslin$Group != "C"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/ns_Family_farmerVSoffice",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)

M2_ns_family_famerVSoffice <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/ns_Family_farmerVSoffice/significant_results.tsv", header = TRUE, "\t")
M2_ns_family_famerVSoffice$shape <- ifelse(M2_ns_family_famerVSoffice$coef > 0, 1, -1)
M2_ns_family_famerVSoffice <- subset(M2_ns_family_famerVSoffice, qval < 0.05)
ggplot(M2_ns_family_famerVSoffice, aes(x=coef, y=reorder(feature, coef), color=shape))+
  geom_point(size=2)+
  labs(x="Relative to Office Worker Nasal Microbiome",
       y="Family", shape="Shape")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_vline(xintercept=0, colour='black',size=1, linetype="dotted")+
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0)) + 
  scale_colour_gradient2() + labs(title = "Maaslin2 on nasal microbiota")

dim(M2_ns_family_famerVSoffice[which(M2_ns_family_famerVSoffice$coef > 0),])

Maaslin2(
  input_data = family_abund_ns,
  input_metadata = metadata_maaslin[which(metadata_maaslin$Sample_Type == "Nasal" & metadata_maaslin$Group != "W"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/ns_Family_cowVSoffice",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)

M2_ns_family_cowVSoffice <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/ns_Family_cowVSoffice/significant_results.tsv", header = TRUE, "\t")
M2_ns_family_cowVSoffice$coef = -1 * M2_ns_family_cowVSoffice$coef
M2_ns_family_cowVSoffice$value <- "C"
M2_ns_family_cowVSoffice$shape <- ifelse(M2_ns_family_cowVSoffice$coef > 0, 1, -1)
M2_ns_family_cowVSoffice <- subset(M2_ns_family_cowVSoffice, qval < 0.05)
ggplot(M2_ns_family_cowVSoffice, aes(x=coef, y=reorder(feature, coef), color=shape))+
  geom_point(size=2)+
  labs(x="Relative to Office Worker Nasal Microbiome",
       y="Family", shape="Shape")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_vline(xintercept=0, colour='black',size=1, linetype="dotted")+
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0)) + 
  scale_colour_gradient2() + labs(title = "Maaslin2 on nasal microbiota of cows and office workers")
p <- ggplot(M2_ns_family_cowVSoffice, aes(x=coef, y=reorder(feature, coef)))+
  geom_point(color = "#873e23", size = 5, alpha = 0.7)+
  geom_vline(xintercept=0, colour='black',size=1.5, linetype="dotted")+
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0), color = "#873e23", size = 1.5, alpha = 1) + 
  theme_linedraw()+
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 14)) +
  labs(x="Coefficient (relative to controls)", y="Genus")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/16S_Nasal_Maaslin2_Cow-vs-Office_Family_v4.svg", p, bg = "transparent", width = 10, height = 30)

dim(M2_ns_genus_cowVSoffice[which(M2_ns_genus_cowVSoffice$coef > 0),])

M2_ns_family_farmerANDcowVSoffice <- M2_ns_family_famerVSoffice[,c("feature", "value", "coef", "stderr")]
colnames(M2_ns_family_farmerANDcowVSoffice)[1] <- "Family"
colnames(M2_ns_family_farmerANDcowVSoffice)[2] <- "Group_W"
colnames(M2_ns_family_farmerANDcowVSoffice)[3] <- "Coef_W"
colnames(M2_ns_family_farmerANDcowVSoffice)[4] <- "StdErr_W"

M2_ns_family_farmerANDcowVSoffice <- merge(M2_ns_family_farmerANDcowVSoffice, M2_ns_family_cowVSoffice[which(M2_ns_family_cowVSoffice$feature %in% M2_ns_family_farmerANDcowVSoffice$Family),c("feature", "value", "coef", "stderr")], by.x = "Family", by.y = "feature", all = TRUE)
colnames(M2_ns_family_farmerANDcowVSoffice)[5] <- "Group_C"
colnames(M2_ns_family_farmerANDcowVSoffice)[6] <- "Coef_C"
colnames(M2_ns_family_farmerANDcowVSoffice)[7] <- "StdErr_C"

M2_ns_family_farmerANDcowVSoffice$MeanCoef <- rowMeans(M2_ns_family_farmerANDcowVSoffice[,c("Coef_W", "Coef_C")], na.rm = TRUE)

ggplot() + 
  geom_point(data = M2_ns_family_farmerANDcowVSoffice, aes(x = Coef_W, y = reorder(Family, MeanCoef)), color = "#18678d", shape = 16, size = 6, alpha = 0.8) +
  geom_errorbarh(data = M2_ns_family_farmerANDcowVSoffice, aes(xmax = Coef_W + StdErr_W, xmin = Coef_W - StdErr_W, y = Family, height = 0), color = "#18678d", size = 1.8, alpha = 0.7) + 
  geom_point(data = M2_ns_family_farmerANDcowVSoffice[which(!is.na(M2_ns_family_farmerANDcowVSoffice$Group_C)),], aes(x = Coef_C, y = Family), color = "#873e23", size = 6, alpha = 0.8) + 
  geom_errorbarh(data = M2_ns_family_farmerANDcowVSoffice[which(!is.na(M2_ns_family_farmerANDcowVSoffice$Group_C)),], aes(xmax = Coef_C + StdErr_C, xmin = Coef_C - StdErr_C, y = Family, height = 0), color = "#873e23", size = 1.8, alpha = 0.7) +
  theme_linedraw()+
  theme(axis.text.x = element_text(), axis.title = element_blank(), axis.text.y = element_text(size = 14)) +
  geom_vline(xintercept=0, colour='black',size=1, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Family")
p <- ggplot() + 
  geom_point(data = M2_ns_family_farmerANDcowVSoffice, aes(x = Coef_W, y = reorder(Family, MeanCoef)), color = "#18678d", shape = 16, size = 6, alpha = 0.8) +
  geom_errorbarh(data = M2_ns_family_farmerANDcowVSoffice, aes(xmax = Coef_W + StdErr_W, xmin = Coef_W - StdErr_W, y = Family, height = 0), color = "#18678d", size = 1.8, alpha = 0.7) + 
  geom_point(data = M2_ns_family_farmerANDcowVSoffice[which(!is.na(M2_ns_family_farmerANDcowVSoffice$Group_C)),], aes(x = Coef_C, y = Family), color = "#873e23", size = 6, alpha = 0.8) + 
  geom_errorbarh(data = M2_ns_family_farmerANDcowVSoffice[which(!is.na(M2_ns_family_farmerANDcowVSoffice$Group_C)),], aes(xmax = Coef_C + StdErr_C, xmin = Coef_C - StdErr_C, y = Family, height = 0), color = "#873e23", size = 1.8, alpha = 0.7) +
  theme_linedraw()+
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 14)) +
  geom_vline(xintercept=0, colour='black',size=1, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Family")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/16S_Nasal_Maaslin2_FarmerCow_Family_v1.svg", p, bg = "transparent", width = 15, height = 20)

dim(M2_ns_family_farmerANDcowVSoffice[which(M2_ns_family_farmerANDcowVSoffice$Coef_C >0),])

###Running Maaslin2 to check for the association of the genera found among shared strains and enriched in farmer noses in the human guts
shared_genera <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d17_inStrain_MAGs/221221_SharedStrains_FarmerNoseEnrichedGenera.csv")

Maaslin2(
  input_data = genus_abund_fc[,c(shared_genera$Genus)],
  input_metadata = metadata_maaslin[which(metadata_maaslin$Sample_Type == "Fecal" & metadata_maaslin$Group != "C"),],
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/16S/221219_DOME_16S_RV/Maaslin2/fc_SharedGenera_farmerVSoffice",
  fixed_effects = c("Group"),
  random_effects = c("SID", "SeasonFull", "Site")
)
