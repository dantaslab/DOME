#########################################################################################
###This script analyzes the Metaphlan4 taxonomic profiles of the DOME  gut microbiomes###
#########################################################################################

###Load the necessary packages
library(ggplot2)
library(reshape2)
library(vegan)
library(Maaslin2)
library(dplyr)

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


###Read in the metadata
samples_metadata <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/220504_DOME_SamplesMasterList_v10.csv", header = TRUE)
samples_metadata <- samples_metadata[which(samples_metadata$Sample_Type == "Fecal"),1:16]
#Two of the cow samples (19-C0205-W001fc01 & 19-C0416-F001fc01) have been clustering with human samples when looking at taxa (shotgun and 16S) and ARGs
#These two samples have likley been misplabeled. Removing these two from further analysis
samples_metadata <- samples_metadata[which(samples_metadata$Sample_ID != "19-C0205-W001fc01" & samples_metadata$Sample_ID != "19-C0416-F001fc01"),]

###Read in the Metaphlan4 species output file
species <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d23_metaphlan4/Extracted/species/Merged_profile_species.txt")
species <- dcast(species, Sample ~ Species, value.var = "Abundance")
species <- species[which(species$Sample != "19-C0205-W001fc01" & species$Sample != "19-C0416-F001fc01"),]
species[is.na(species)] <- 0
species_abundfilt_01 <- species
species_abundfilt_01[species_abundfilt_01 < 0.1] <- 0
rownames(species_abundfilt_01) <- species_abundfilt_01$Sample
species_abundfilt_01 <- species_abundfilt_01[2:dim(species_abundfilt_01)[2]]
species_abundfilt_01 <- species_abundfilt_01 / rowSums(species_abundfilt_01)
species_abundfilt_01 <- species_abundfilt_01[, colSums(species_abundfilt_01) > 0]
species_human <- species_abundfilt_01[which(rownames(species_abundfilt_01) %in% samples_metadata[which(samples_metadata$Group != "C"),]$Sample_ID),]
species_human <- species_human[, colSums(species_human) > 0]

###Read in the Metaphlan4 genus output file
genus <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d23_metaphlan4/Extracted/genus/Merged_profile_genus.txt")
genus <- dcast(genus, Sample ~ Genus, value.var = "Abundance")
genus <- genus[which(genus$Sample != "19-C0205-W001fc01" & genus$Sample != "19-C0416-F001fc01"),]
genus[is.na(genus)] <- 0
genus_abundfilt_01 <- genus
genus_abundfilt_01[genus_abundfilt_01 < 0.1] <- 0
rownames(genus_abundfilt_01) <- genus_abundfilt_01$Sample
genus_abundfilt_01 <- genus_abundfilt_01[2:dim(genus_abundfilt_01)[2]]
genus_abundfilt_01 <- genus_abundfilt_01 / rowSums(genus_abundfilt_01)
genus_abundfilt_01 <- genus_abundfilt_01[, colSums(genus_abundfilt_01) > 0]
genus_human <- genus_abundfilt_01[which(rownames(genus_abundfilt_01) %in% samples_metadata[which(samples_metadata$Group != "C"),]$Sample_ID),]
genus_human <- genus_human[, colSums(genus_human) > 0]

###Read in the Metaphlan4 family output file
family <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d23_metaphlan4/Extracted/family/Merged_profile_family.txt")
family <- dcast(family, Sample ~ Family, value.var = "Abundance")
family <- family[which(family$Sample != "19-C0205-W001fc01" & family$Sample != "19-C0416-F001fc01"),]
family[is.na(family)] <- 0
family_abundfilt_01 <- family
family_abundfilt_01[family_abundfilt_01 < 0.1] <- 0
rownames(family_abundfilt_01) <- family_abundfilt_01$Sample
family_abundfilt_01 <- family_abundfilt_01[2:dim(family_abundfilt_01)[2]]
family_abundfilt_01 <- family_abundfilt_01 / rowSums(family_abundfilt_01)
family_abundfilt_01 <- family_abundfilt_01[, colSums(family_abundfilt_01) > 0]
family_human <- family_abundfilt_01[which(rownames(family_abundfilt_01) %in% samples_metadata[which(samples_metadata$Group != "C"),]$Sample_ID),]
family_human <- family_human[, colSums(family_human) > 0]

###Get the alpha diversity for all samples
richness <- as.data.frame(specnumber(genus_abundfilt_01))
richness$Sample <- row.names(richness)
colnames(richness)[1] <-"Richness"

shannon_div <- as.data.frame(diversity(genus_abundfilt_01, index = "shannon"))
shannon_div$Sample <- row.names(shannon_div)
colnames(shannon_div)[1] <- "ShannonDiv"

alpha_div <- merge(richness, shannon_div, by = "Sample")
alpha_div <- merge(alpha_div, samples_metadata[, c(1,3,12,13,16)], by.x = "Sample", by.y = "Sample_ID")
alpha_div$SeasonFull <- paste(alpha_div$Season, alpha_div$Year)
alpha_div$SeasonFull <- factor(alpha_div$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))
alpha_div$Group <- factor(alpha_div$Group, levels = c("C", "W", "D"))

ggplot(alpha_div, aes(x = Group, y = Richness, color = Group)) + geom_boxplot() + 
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_test() +
  theme(legend.position = "none") + ylab("Genus Richness (0.1% threshold)") +
  scale_x_discrete(labels = c("Cow", "Farmer", "Office"))
ggplot(alpha_div[which(alpha_div$Season == "SPRING" & alpha_div$Year == 1),], aes(x = Group, y = Richness, color = Group)) + geom_boxplot(aes(fill = Group), color = "black", outlier.shape = NA, width = 0.4, size = 1, alpha = 0.5) + 
  geom_flat_violin(aes(fill = Group), position = position_nudge(x = 0.24), alpha = 0.7, scale = "width", width = 0.45) +
  geom_point(position = position_dodge2(width = 0.4), alpha = 0.8) +
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_linedraw() + ylab("Genus richness") +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_text(size = 15), axis.text.x = element_text(face = "bold", size = 20), strip.text = element_blank(),
        legend.position = "none", axis.title.y = element_text(face = "bold", size = 20)) +
  scale_x_discrete(labels = c("Cows", "Farmer", "Non-farmer")) + scale_color_manual(values = c("#873e23", "#18678d", "#626262")) + 
  scale_fill_manual(values = c("#873e23", "#18678d", "#626262")) + xlab(element_blank())
compare_means(Richness ~ Group, data = alpha_div[which(alpha_div$Season == "SPRING" & alpha_div$Year == 1),], method = "wilcox.test", p.adjust.method = "BH")


ggplot(alpha_div, aes(x = Group, y = ShannonDiv, color = Group)) + geom_boxplot() + 
  facet_wrap(~ SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_test() +
  theme(legend.position = "none") + ylab("Shannon Diversity (0.1% threshold)") +
  scale_x_discrete(labels = c("Cow", "Farmer", "Office"))

###Looking only at the genera that were identified during strain sharing
genus_human_SharedStrainGenera <- as.data.frame(rowSums(genus_human[,c("g__Bifidobacterium", "g__Romboutsia", "g__Treponema", "g__Prevotella", "g__Turicibacter", "g__Holdemanella", "g__Streptococcus", "g__Escherichia", "g__Blautia")]))
genus_human_SharedStrainGenera <- as.data.frame(genus_human[1:285,c("g__Prevotella")])
rownames(genus_human_SharedStrainGenera) <- rownames(genus_human)
colnames(genus_human_SharedStrainGenera)[1] <- "TotRelAbund"
genus_human_SharedStrainGenera <- merge(genus_human_SharedStrainGenera, samples_metadata[, c(1,3,12,13,16)], by.x = "row.names", by.y = "Sample_ID")
ggplot(data = genus_human_SharedStrainGenera) + 
  geom_boxplot(aes(x = Group, y = TotRelAbund))
compare_means(TotRelAbund ~ Group, data = genus_human_SharedStrainGenera, method = "wilcox.test")

###Calculate beta diversity for all samples
Species_BrayCurtis <- vegdist(species_abundfilt_01, method = "bray")
pcoa_species_BC <- pco(Species_BrayCurtis, k = 4)
pcoa_species_BC <- as.data.frame(pcoa_species_BC$points)
pcoa_species_BC <- merge(pcoa_species_BC, samples_metadata[, c(1,3,12,13,16)], by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_species_BC)[1] <- "Sample"
pcoa_species_BC$SeasonFull <- paste(pcoa_species_BC$Season, pcoa_species_BC$Year)

ggplot(data = pcoa_species_BC, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull)) +
  theme_test()


###Running Maaslin2 to get the species and genera associated with farmers relative to office workers
maaslin_metadata <- samples_metadata %>% dplyr::select(3, 1, 12:16)
maaslin_metadata$SeasonFull <- paste(maaslin_metadata$Season, maaslin_metadata$Year)
rownames(maaslin_metadata) <- maaslin_metadata$Sample_ID

Maaslin2(
  input_data = species_abundfilt_01,
  input_metadata = maaslin_metadata[which(maaslin_metadata$Group != "C"),],
  fixed_effects = "Group",
  random_effects = c("SID", "SeasonFull", "Site"),
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d23_metaphlan4/Maaslin2/Species_farmerVSoffice"
)

#M2_species_famerVSoffice <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d23_metaphlan4/Maaslin2/Species_farmerVSoffice/significant_results.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
M2_species_famerVSoffice <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d23_metaphlan4/Maaslin2/Species_farmerVSoffice/significant_results.tsv", header = TRUE, "\t")
M2_species_famerVSoffice$shape <- ifelse(M2_species_famerVSoffice$coef > 0, 1, -1)
M2_species_famerVSoffice <- subset(M2_species_famerVSoffice, qval < 0.05)
p <- ggplot(M2_species_famerVSoffice, aes(x=coef, y=reorder(feature, coef), color=coef)) +
  geom_point(color = "#18678d", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0), color = "#18678d", size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 20, face = "bold"), axis.text.y = element_text(size = 14)) +
  geom_vline(xintercept=0, colour='black',size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Species")
p <- ggplot(M2_species_famerVSoffice, aes(x=coef, y=reorder(feature, coef), color=coef)) +
  geom_point(color = "#18678d", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0), color = "#18678d", size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_blank()) +
  geom_vline(xintercept=0, colour='black',size=2, linetype="dotted")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Shotgun_Maaslin2_FarmerOffice_Species_v3.svg", p, bg = "transparent", width = 5, height = 10)


Maaslin2(
  input_data = genus_abundfilt_01,
  input_metadata = maaslin_metadata[which(maaslin_metadata$Group != "C"),],
  fixed_effects = "Group",
  random_effects = c("SID", "SeasonFull", "Site"),
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d23_metaphlan4/Maaslin2/Genus_farmerVSoffice"
)

#M2_genus_famerVSoffice <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d23_metaphlan4/Maaslin2/Genus_farmerVSoffice/significant_results.tsv", header = TRUE, "\t", escape_double = FALSE, trim_ws = TRUE)
M2_genus_famerVSoffice <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d23_metaphlan4/Maaslin2/Genus_farmerVSoffice/significant_results.tsv", header = TRUE, "\t")
M2_genus_famerVSoffice$shape <- ifelse(M2_genus_famerVSoffice$coef > 0, 1, -1)
M2_genus_famerVSoffice <- subset(M2_genus_famerVSoffice, qval < 0.05)
ggplot(M2_genus_famerVSoffice, aes(x=coef, y=reorder(feature, coef), color=coef))+
  geom_point(aes(shape=factor(shape)), size=2)+
  labs(x="Relative to Office Worker Gut Microbiome",
       y="Genus", shape="Shape")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  geom_vline(xintercept=0, colour='black',size=1, linetype="dotted")+
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0)) + 
  scale_colour_gradient2()
p <- ggplot(M2_genus_famerVSoffice, aes(x=coef, y=reorder(feature, coef), color=coef)) +
  geom_point(color = "#18678d", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0), color = "#18678d", size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_blank()) +
  geom_vline(xintercept=0, colour='black',size=2, linetype="dotted")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/Shotgun_Maaslin2_FarmerOffice_Genus_v2.svg", p, bg = "transparent", width = 5, height = 10)

Maaslin2(
  input_data = family_abundfilt_01,
  input_metadata = maaslin_metadata[which(maaslin_metadata$Group != "C"),],
  fixed_effects = "Group",
  random_effects = c("SID", "SeasonFull", "Site"),
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d23_metaphlan4/Maaslin2/Family_farmerVSoffice"
)

###Running Maaslin2 to get the species associated with cows relative to office workers
Maaslin2(
  input_data = species_abundfilt_01,
  input_metadata = maaslin_metadata[which(maaslin_metadata$Group != "W"),],
  fixed_effects = "Group",
  random_effects = c("SID", "SeasonFull", "Site"),
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d23_metaphlan4/Maaslin2/Species_cowVSoffice"
)

###Looking at the Prevotella abundance in humans
genus_human_Prevotella <- genus_human
genus_human_Prevotella$Sample <- rownames(genus_human_Prevotella)
genus_human_Prevotella <- genus_human_Prevotella[,c("Sample", "g__Prevotella")]
genus_human_Prevotella <- merge(genus_human_Prevotella, samples_metadata[, c(1,3,12,13,16)], by.x = "Sample", by.y = "Sample_ID")
genus_human_Prevotella$SeasonFull <- paste(genus_human_Prevotella$Season, genus_human_Prevotella$Year)
genus_human_Prevotella$SeasonFull <- factor(genus_human_Prevotella$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))
genus_human_Prevotella$Group <- factor(genus_human_Prevotella$Group, levels = c("W", "D"))

ggplot(data = genus_human_Prevotella, aes(x = Group, y = g__Prevotella)) + 
  geom_boxplot(aes(fill = Group), color = "black", outlier.shape = NA, size = 2, alpha = 0.5) +
  geom_point(aes(fill = Group), position = position_dodge2(width = 0.5), alpha = 0.6, size = 3, pch = 21, color = "black") +
  theme_linedraw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16, face = "bold"), strip.text = element_text(face = "bold", size = 20),
        legend.position = "none", axis.title.y = element_text(face = "bold", size = 18)) +
  scale_x_discrete(labels = c("Farmers", "Controls")) + scale_color_manual(values = c("#18678d", "#626262")) + 
  scale_fill_manual(values = c("#18678d", "#626262")) + xlab(element_blank())
p <- ggplot(data = genus_human_Prevotella, aes(x = Group, y = g__Prevotella)) + 
  geom_boxplot(aes(fill = Group), color = "black", outlier.shape = NA, size = 2, alpha = 0.5) +
  geom_point(aes(fill = Group), position = position_dodge2(width = 0.6), alpha = 0.6, size = 3, pch = 21, color = "black") +
  theme_linedraw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text = element_blank(),
        legend.position = "none", axis.title = element_blank(), panel.border = element_rect(linewidth = 1.5, color = "black")) +
  scale_color_manual(values = c("#18678d", "#626262")) + 
  scale_fill_manual(values = c("#18678d", "#626262")) + xlab(element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Shotgun_Maaslin2_FarmerOffice_Prevotella_v1.svg", p, bg = "transparent", width = 5, height = 5)

###Calculating beta diversity for human samples only at the species level
human_spec_BrayCrutis <- vegdist(species_human, method = "bray")
pcoa_human_spec_BC <- pco(human_spec_BrayCrutis, k = 4)
pcoa_human_spec_BC <- as.data.frame(pcoa_human_spec_BC$points)
pcoa_human_spec_BC <- merge(pcoa_human_spec_BC, samples_metadata[, c(1,3,12,13,16)], by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_human_spec_BC)[1] <- "Sample"
pcoa_human_spec_BC$SeasonFull <- paste(pcoa_human_spec_BC$Season, pcoa_human_spec_BC$Year)

ggplot(data = pcoa_human_spec_BC, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull)) +
  theme_test() + stat_ellipse(aes(color = Group)) + labs(title = "Species PCOA (BC) of farmer and office worker guts")

###Looking at the species-level communities within each of the season
samples_metadata_season <- samples_metadata[which(samples_metadata$Group != "C" & samples_metadata$Sample_Type == "Fecal" & samples_metadata$Season == "WINTER" & samples_metadata$Year == "1"),]

species_human_season <- species_human[rownames(species_human) %in% samples_metadata_season$Sample_ID,]
human_spec_BC_season <- vegdist(species_human_season, method = "bray")
pcoa_human_spec_BC_season <- pco(human_spec_BC_season, k = 4)
pcoa_human_spec_BC_season_scree <- as.data.frame(100*pcoa_human_spec_BC_season$eig/sum(pcoa_human_spec_BC_season$eig))
pcoa_human_spec_BC_season <- as.data.frame(pcoa_human_spec_BC_season$points)
pcoa_human_spec_BC_season <- merge(pcoa_human_spec_BC_season, samples_metadata[, c(1,3,12,13,16)], by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_human_spec_BC_season)[1] <- "Sample"
pcoa_human_spec_BC_season$Group <- factor(pcoa_human_spec_BC_season$Group, levels = c("W", "D"))

p <- ggplot(data = pcoa_human_spec_BC_season, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull), alpha = 1, size = 5, shape = 16) +
  theme_test() + stat_ellipse(aes(color = Group), size = 2.5, alpha = 0.75) +
  scale_color_manual(values = c("#18678d", "#626262"), labels = c("Farmers", "Controls")) + xlab("PCoA 1 (10.7%)") + ylab("PCoA 2 (9.4%)") +
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"), axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 20),
        legend.position = "none")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/Shotgun_Species_Human_PCOA_Winter2019_v1.svg", p, bg = "transparent", width = 7, height = 7)

species_human_season <- merge(species_human_season, samples_metadata_season[, c(3, 13)], by.x = "row.names", by.y = "Sample_ID")
row.names(species_human_season) <- species_human_season$Row.names
species_human_season <- species_human_season[,c(-1)]
adonis2(species_human_season[,c(-773)] ~ species_human_season$Group,
        permutations = 9999, method = "bray")

p.adjust(c(0.0269, 0.081, 0.0196, 0.0088), method = "BH")

#samples_metadata_season <- samples_metadata_season[order(samples_metadata_season$Sample_ID),]
#species_human_season <- species_human[rownames(species_human) %in% samples_metadata_season$Sample_ID,]
#species_human_season <- species_human_season[, colSums(species_human_season) > 0]
#species_human_season <- species_human_season[order(rownames(species_human_season)),]

#adonis2(species_human_season ~ samples_metadata_season$Group,
#        permutations = 9999, method = "bray")

#p.adjust(c(0.0296, 0.0769, 0.0207, 0.0073, 0.2368), method = "BH")
#0.04933333 0.09612500 0.04933333 0.03650000 0.23680000

genus_human_season <- genus_human[rownames(genus_human) %in% samples_metadata_season$Sample_ID,]

human_genus_BC_season <- vegdist(genus_human_season, method = "bray")
pcoa_human_genus_BC_season <- pco(human_genus_BC_season, k = 4)
pcoa_human_genus_BC_season_scree <- as.data.frame(100*pcoa_human_genus_BC_season$eig/sum(pcoa_human_genus_BC_season$eig))
pcoa_human_genus_BC_season <- as.data.frame(pcoa_human_genus_BC_season$points)
pcoa_human_genus_BC_season <- merge(pcoa_human_genus_BC_season, samples_metadata[, c(1,3,12,13,16)], by.x = "row.names", by.y = "Sample_ID")
colnames(pcoa_human_genus_BC_season)[1] <- "Sample"
pcoa_human_genus_BC_season$Group <- factor(pcoa_human_genus_BC_season$Group, levels = c("W", "D"))

p <- ggplot(data = pcoa_human_genus_BC_season, aes(x = V1, y = V2)) + geom_point(aes(color = Group, shape = SeasonFull), alpha = 1, size = 5, shape = 16) +
  theme_test() + stat_ellipse(aes(color = Group), size = 2.5, alpha = 0.75) +
  scale_color_manual(values = c("#18678d", "#626262"), labels = c("Farmers", "Controls")) + xlab("PCoA 1 (16.0%)") + ylab("PCoA 2 (11.4%)") +
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"), axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 20),
        legend.position = "none")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/Shotgun_Genus_Human_PCOA_Winter2019_v1.svg", p, bg = "transparent", width = 7, height = 7)

genus_human_season <- merge(genus_human_season, samples_metadata_season[, c(3, 13)], by.x = "row.names", by.y = "Sample_ID")
row.names(genus_human_season) <- genus_human_season$Row.names
genus_human_season <- genus_human_season[,c(-1)]
adonis2(genus_human_season[,c(-419)] ~ genus_human_season$Group,
        permutations = 9999, method = "bray")

p.adjust(c(0.1082, 0.1836, 0.041, 0.0283), method = "BH")
#0.1442667 0.1836000 0.0820000 0.0820000



###Comparing the pairwise BC distances based on species (Metaphlan4) and kmers (Simka)
human_spec_BC_pairwise <- as.matrix(human_spec_BrayCrutis)
human_spec_BC_pairwise <- melt(human_spec_BC_pairwise)
human_spec_BC_pairwise <- human_spec_BC_pairwise[!duplicated(apply(human_spec_BC_pairwise, 1, function(x) paste(sort(x), collapse = ""))),]
colnames(human_spec_BC_pairwise)[1] <- "Sample1"
colnames(human_spec_BC_pairwise)[2] <- "Sample2"
colnames(human_spec_BC_pairwise)[3] <- "Spec_BC"
human_spec_BC_pairwise <- human_spec_BC_pairwise[which(human_spec_BC_pairwise$Sample1 != human_spec_BC_pairwise$Sample2),]
human_spec_BC_pairwise$Spec_BC <- as.numeric(human_spec_BC_pairwise$Spec_BC)

kmer_16 <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d21_Simka/k16/220930_k16_BC.csv", header = FALSE)
names(kmer_16) <- kmer_16[1,]
kmer_16 <- kmer_16[-1,]
colnames(kmer_16)[1] <- "Sample"
rownames(kmer_16) <- kmer_16$Sample
kmer_16 <- kmer_16[,-1]
kmer_16 <- melt(as.matrix(kmer_16))
colnames(kmer_16)[1] <- "Sample1"
colnames(kmer_16)[2] <- "Sample2"
colnames(kmer_16)[3] <- "kmer_16"
kmer_16$kmer_16 <- as.numeric(kmer_16$kmer_16)

kmer_15 <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d21_Simka/k15/220929_k15_BC.csv", header = FALSE)
names(kmer_15) <- kmer_15[1,]
kmer_15 <- kmer_15[-1,]
colnames(kmer_15)[1] <- "Sample"
rownames(kmer_15) <- kmer_15$Sample
kmer_15 <- kmer_15[,-1]
kmer_15 <- melt(as.matrix(kmer_15))
colnames(kmer_15)[1] <- "Sample1"
colnames(kmer_15)[2] <- "Sample2"
colnames(kmer_15)[3] <- "kmer_15"
kmer_15$kmer_15 <- as.numeric(kmer_15$kmer_15)

kmer_17 <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d21_Simka/k17/220928_k17_BC.csv", header = FALSE)
names(kmer_17) <- kmer_17[1,]
kmer_17 <- kmer_17[-1,]
colnames(kmer_17)[1] <- "Sample"
rownames(kmer_17) <- kmer_17$Sample
kmer_17 <- kmer_17[,-1]
kmer_17 <- melt(as.matrix(kmer_17))
colnames(kmer_17)[1] <- "Sample1"
colnames(kmer_17)[2] <- "Sample2"
colnames(kmer_17)[3] <- "kmer_17"
kmer_17$kmer_17 <- as.numeric(kmer_17$kmer_17)

kmer_18 <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d21_Simka/k18/221013_k18_BC.csv", header = FALSE)
names(kmer_18) <- kmer_18[1,]
kmer_18 <- kmer_18[-1,]
colnames(kmer_18)[1] <- "Sample"
rownames(kmer_18) <- kmer_18$Sample
kmer_18 <- kmer_18[,-1]
kmer_18 <- melt(as.matrix(kmer_18))
colnames(kmer_18)[1] <- "Sample1"
colnames(kmer_18)[2] <- "Sample2"
colnames(kmer_18)[3] <- "kmer_18"
kmer_18$kmer_18 <- as.numeric(kmer_18$kmer_18)

kmer_19 <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d21_Simka/k19/221013_k19_BC.csv", header = FALSE)
names(kmer_19) <- kmer_19[1,]
kmer_19 <- kmer_19[-1,]
colnames(kmer_19)[1] <- "Sample"
rownames(kmer_19) <- kmer_19$Sample
kmer_19 <- kmer_19[,-1]
kmer_19 <- melt(as.matrix(kmer_19))
colnames(kmer_19)[1] <- "Sample1"
colnames(kmer_19)[2] <- "Sample2"
colnames(kmer_19)[3] <- "kmer_19"
kmer_19$kmer_19 <- as.numeric(kmer_19$kmer_19)

kmer_20 <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d21_Simka/k20/221013_k20_BC.csv", header = FALSE)
names(kmer_20) <- kmer_20[1,]
kmer_20 <- kmer_20[-1,]
colnames(kmer_20)[1] <- "Sample"
rownames(kmer_20) <- kmer_20$Sample
kmer_20 <- kmer_20[,-1]
kmer_20 <- melt(as.matrix(kmer_20))
colnames(kmer_20)[1] <- "Sample1"
colnames(kmer_20)[2] <- "Sample2"
colnames(kmer_20)[3] <- "kmer_20"
kmer_20$kmer_20 <- as.numeric(kmer_20$kmer_20)

kmer_21 <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d21_Simka/k21/220927_k21_BC.csv", header = FALSE)
names(kmer_21) <- kmer_21[1,]
kmer_21 <- kmer_21[-1,]
colnames(kmer_21)[1] <- "Sample"
rownames(kmer_21) <- kmer_21$Sample
kmer_21 <- kmer_21[,-1]
kmer_21 <- melt(as.matrix(kmer_21))
colnames(kmer_21)[1] <- "Sample1"
colnames(kmer_21)[2] <- "Sample2"
colnames(kmer_21)[3] <- "kmer_21"
kmer_21$kmer_21 <- as.numeric(kmer_21$kmer_21)

human_spec_BC_pairwise <- merge(human_spec_BC_pairwise, kmer_15, by = c("Sample1", "Sample2"))
human_spec_BC_pairwise <- merge(human_spec_BC_pairwise, kmer_16, by = c("Sample1", "Sample2"))
human_spec_BC_pairwise <- merge(human_spec_BC_pairwise, kmer_17, by = c("Sample1", "Sample2"))
human_spec_BC_pairwise <- merge(human_spec_BC_pairwise, kmer_18, by = c("Sample1", "Sample2"))
human_spec_BC_pairwise <- merge(human_spec_BC_pairwise, kmer_19, by = c("Sample1", "Sample2"))
human_spec_BC_pairwise <- merge(human_spec_BC_pairwise, kmer_20, by = c("Sample1", "Sample2"))
human_spec_BC_pairwise <- merge(human_spec_BC_pairwise, kmer_21, by = c("Sample1", "Sample2"))


ggplot(data = human_spec_BC_pairwise, aes(x = Spec_BC, y = kmer_15)) + geom_point() + geom_smooth(method = "lm", se = FALSE) + theme_test()
ggplot(data = human_spec_BC_pairwise, aes(x = Spec_BC, y = kmer_16)) + geom_point() + geom_smooth(method = "lm", se = FALSE) + theme_test()
ggplot(data = human_spec_BC_pairwise, aes(x = Spec_BC, y = kmer_17)) + geom_point() + geom_smooth(method = "lm", se = FALSE) + theme_test()
ggplot(data = human_spec_BC_pairwise, aes(x = Spec_BC, y = kmer_18)) + geom_point() + geom_smooth(method = "lm", se = FALSE) + theme_test()
ggplot(data = human_spec_BC_pairwise, aes(x = Spec_BC, y = kmer_19)) + geom_point() + geom_smooth(method = "lm", se = FALSE) + theme_test()
ggplot(data = human_spec_BC_pairwise, aes(x = Spec_BC, y = kmer_20)) + geom_point() + geom_smooth(method = "lm", se = FALSE) + theme_test()
ggplot(data = human_spec_BC_pairwise, aes(x = Spec_BC, y = kmer_21)) + geom_point() + geom_smooth(method = "lm", se = FALSE) + theme_test()

summary(lm(human_spec_BC_pairwise$Spec_BC ~ human_spec_BC_pairwise$kmer_15)) #R2 = 0.7819
summary(lm(human_spec_BC_pairwise$Spec_BC ~ human_spec_BC_pairwise$kmer_16)) #R2 = 0.8725
summary(lm(human_spec_BC_pairwise$Spec_BC ~ human_spec_BC_pairwise$kmer_17)) #R2 = 0.8819
summary(lm(human_spec_BC_pairwise$Spec_BC ~ human_spec_BC_pairwise$kmer_18)) #R2 = 0.8793
summary(lm(human_spec_BC_pairwise$Spec_BC ~ human_spec_BC_pairwise$kmer_19)) #R2 = 0.8771
summary(lm(human_spec_BC_pairwise$Spec_BC ~ human_spec_BC_pairwise$kmer_20)) #R2 = 0.8754
summary(lm(human_spec_BC_pairwise$Spec_BC ~ human_spec_BC_pairwise$kmer_21)) #R2 = 0.8737

kmer_R2 <- data.frame(kmer = c(15, 16, 17, 18, 19, 20, 21),
           R2 = c(0.7819, 0.8725, 0.8819, 0.8793, 0.8771, 0.8754, 0.8737))

ggplot(data = kmer_R2, aes(x = kmer, y = R2)) + geom_point() + geom_line() + theme_bw() + scale_x_continuous(name = "kmer size", breaks = seq(15,21,1))
