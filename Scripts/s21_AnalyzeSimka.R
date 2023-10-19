######################################################################################
###This script analyzes the kmer-based beta-diversity profiling of DOME metagenomes###
######################################################################################

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
#library(coin)

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

###Read in the BC distance matrix based on kmer counts from Simka
kmer_BC <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d21_Simka/k17/220928_k17_BC.csv", header = FALSE)
names(kmer_BC) <- kmer_BC[1,]
kmer_BC <- kmer_BC[-1,]
colnames(kmer_BC)[1] <- "Sample"
kmer_BC <- kmer_BC[which(kmer_BC$Sample != "19-C0205-W001fc01" & kmer_BC$Sample != "19-C0416-F001fc01"),]
rownames(kmer_BC) <- kmer_BC$Sample
kmer_BC <- kmer_BC[,-1]
kmer_BC <- kmer_BC[,which(colnames(kmer_BC) %in% rownames(kmer_BC))]

###Getting the pairwise Bray-Curtis distances for sample kmers
kmer_BC_pairwise <- melt(as.matrix(kmer_BC))
colnames(kmer_BC_pairwise)[1] <- "Sample1"
colnames(kmer_BC_pairwise)[2] <- "Sample2"
colnames(kmer_BC_pairwise)[3] <- "BCdist"
#kmer_BC_pairwise <- kmer_BC_pairwise[which(kmer_BC_pairwise$Sample1 %in% samples_metadata[which(samples_metadata$Group != "C"),]$Sample_ID),]
kmer_BC_pairwise <- kmer_BC_pairwise[!duplicated(apply(kmer_BC_pairwise, 1, function(x) paste(sort(x), collapse = ""))),]
kmer_BC_pairwise <- kmer_BC_pairwise[which(kmer_BC_pairwise$Sample1 != kmer_BC_pairwise$Sample2),]
kmer_BC_pairwise <- merge(kmer_BC_pairwise, samples_metadata[,c(1,3,12,13,15,16)], by.x = "Sample1", by.y = "Sample_ID")
colnames(kmer_BC_pairwise)[4] <- "SID1"
colnames(kmer_BC_pairwise)[5] <- "Season 1"
colnames(kmer_BC_pairwise)[6] <- "Group 1"
colnames(kmer_BC_pairwise)[7] <- "Site 1"
colnames(kmer_BC_pairwise)[8] <- "Year 1"
kmer_BC_pairwise <- merge(kmer_BC_pairwise, samples_metadata[,c(1,3,12,13,15,16)], by.x = "Sample2", by.y = "Sample_ID")
colnames(kmer_BC_pairwise)[9] <- "SID2"
colnames(kmer_BC_pairwise)[10] <- "Season 2"
colnames(kmer_BC_pairwise)[11] <- "Group 2"
colnames(kmer_BC_pairwise)[12] <- "Site 2"
colnames(kmer_BC_pairwise)[13] <- "Year 2"

kmer_BC_pairwise$Species1[kmer_BC_pairwise$`Group 1` == "C"] <- "C"
kmer_BC_pairwise$Species1[kmer_BC_pairwise$`Group 1` != "C"] <- "H"
kmer_BC_pairwise$Species2[kmer_BC_pairwise$`Group 2` == "C"] <- "C"
kmer_BC_pairwise$Species2[kmer_BC_pairwise$`Group 2` != "C"] <- "H"
kmer_BC_pairwise$SpeciesJoint <- paste(kmer_BC_pairwise$Species1, kmer_BC_pairwise$Species2)
kmer_BC_pairwise$SpeciesRelat[kmer_BC_pairwise$SpeciesJoint == "C C"] <- "Cow-Cow"
kmer_BC_pairwise$SpeciesRelat[kmer_BC_pairwise$SpeciesJoint == "C H"] <- "Cow-Human"
kmer_BC_pairwise$SpeciesRelat[kmer_BC_pairwise$SpeciesJoint == "H C"] <- "Cow-Human"
kmer_BC_pairwise$SpeciesRelat[kmer_BC_pairwise$SpeciesJoint == "H H"] <- "Human-Human"
kmer_BC_pairwise$BCdist <- as.numeric(kmer_BC_pairwise$BCdist)

ggplot(data = kmer_BC_pairwise, aes(x = kmer_BC_pairwise$SpeciesRelat, y = kmer_BC_pairwise$BCdist, color = kmer_BC_pairwise$SpeciesRelat)) + geom_boxplot() + theme_test() +
  ylab("Bray-Curtis distance (kmer = 17)") + xlab("Host pairs")

kmer_BC_pairwise_filtered <- kmer_BC_pairwise[which(kmer_BC_pairwise$`Season 1` == kmer_BC_pairwise$`Season 2` & kmer_BC_pairwise$`Year 1` == kmer_BC_pairwise$`Year 2`),]
kmer_BC_pairwise_filtered <- kmer_BC_pairwise_filtered[which(kmer_BC_pairwise_filtered$`Group 1` != kmer_BC_pairwise_filtered$`Group 2`),]
kmer_BC_pairwise_filtered <- kmer_BC_pairwise_filtered[,c(1:13)]
kmer_BC_pairwise_filtered$Groups <- paste(kmer_BC_pairwise_filtered$`Group 1`, kmer_BC_pairwise_filtered$`Group 2`)
kmer_BC_pairwise_filtered <- kmer_BC_pairwise_filtered[which(kmer_BC_pairwise_filtered$Groups != "D W" & kmer_BC_pairwise_filtered$Groups != "W D"),]
kmer_BC_pairwise_filtered$Sites[kmer_BC_pairwise_filtered$`Site 1` == kmer_BC_pairwise_filtered$`Site 2`] <- "Same site"
kmer_BC_pairwise_filtered$Sites[kmer_BC_pairwise_filtered$`Site 1` != kmer_BC_pairwise_filtered$`Site 2`] <- "Different sites"
kmer_BC_pairwise_filtered$SeasonFull <- paste(kmer_BC_pairwise_filtered$`Season 1`, kmer_BC_pairwise_filtered$`Year 1`)
kmer_BC_pairwise_filtered$HumanGroup <- ""

for (i in seq(1,dim(kmer_BC_pairwise_filtered)[1])) {
  if (kmer_BC_pairwise_filtered$`Group 1`[i] == "C") {
    kmer_BC_pairwise_filtered$HumanGroup[i] = kmer_BC_pairwise_filtered$`Group 2`[i]
  } else {
    kmer_BC_pairwise_filtered$HumanGroup[i] = kmer_BC_pairwise_filtered$`Group 1`[i]
  }
}

kmer_BC_pairwise_filtered$HumanSID <- ""

for (i in seq(1,dim(kmer_BC_pairwise_filtered)[1])) {
  if (kmer_BC_pairwise_filtered$`Group 1`[i] == "C") {
    kmer_BC_pairwise_filtered$HumanSID[i] = kmer_BC_pairwise_filtered$SID2[i]
  } else {
    kmer_BC_pairwise_filtered$HumanSID[i] = kmer_BC_pairwise_filtered$SID1[i]
  }
}

kmer_BC_pairwise_filtered <- kmer_BC_pairwise_filtered[c(3,15:18)]

kmer_BC_pairwise_filtered <- kmer_BC_pairwise_filtered %>%
  group_by(HumanGroup, HumanSID, SeasonFull, Sites) %>%
  summarise(
    BCdist = mean(BCdist)
  )

kmer_BC_pairwise_filtered$GroupSites <- paste(kmer_BC_pairwise_filtered$HumanGroup, kmer_BC_pairwise_filtered$Sites)
kmer_BC_pairwise_filtered$SeasonFull <- factor(kmer_BC_pairwise_filtered$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))

BC_kmer_pairwise_compared <- compare_means(BCdist ~ GroupSites, data = kmer_BC_pairwise_filtered[,c(3,5,6)], group.by = "SeasonFull", method = "wilcox.test", p.adjust.method = "BH")
BC_kmer_pairwise_compared <- BC_kmer_pairwise_compared[which(BC_kmer_pairwise_compared$p.adj < 0.05),]

wilcox_effsize(data = kmer_BC_pairwise_filtered[which(kmer_BC_pairwise_filtered$SeasonFull == "WINTER 1"),], BCdist ~ GroupSites)

ggplot(data = kmer_BC_pairwise_filtered, aes(x = Sites, y = (BCdist), color = HumanGroup)) + geom_boxplot() + 
  facet_wrap(~kmer_BC_pairwise_filtered$SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_test() + ylab("Bray-Curtis distance") +
  ylim(0.75,1) +
  labs(title = "kmer Bray-Curtis distances to cow samples", subtitle = "Wilcoxon rank sum test. Benjamini-Hochberg correction.") +
  scale_color_manual(values=c("#00BA38", "#619CFF"), name = "", labels = c("Office worker", "Farmer"))

kmer_BC_summarized <- kmer_BC_pairwise_filtered

kmer_BC_summarized <- kmer_BC_summarized %>%
  group_by(HumanGroup, Sites, SeasonFull) %>%
  summarise(
    BCdist_mean = mean(BCdist),
    BCdist_median = median(BCdist)
  )