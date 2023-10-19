#########################################################
###This script analyzes the output of inStrain-Compare###
#########################################################

###Load the necessary packages
library(ggplot2)
library(ggpubr)
library(stats)

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

###Creating a dataframe with all unique pairwise combinations of samples
all_samples <- samples_metadata$Sample_ID
sample_pairs <- as.data.frame(apply(combn(all_samples,2),2,paste,collapse='_'))
colnames(sample_pairs)[1] <- "Sample_pairs"
sample_pairs <- as.data.frame(do.call('rbind', strsplit(as.character(sample_pairs$Sample_pairs),'_',fixed=TRUE)))
colnames(sample_pairs)[1] <- "sample1"
colnames(sample_pairs)[2] <- "sample2"
sample_pairs <- merge(sample_pairs, samples_metadata[,c(1,3,12,13,15,16)], by.x = "sample1", by.y = "Sample_ID")
colnames(sample_pairs)[3] <- "SID1"
colnames(sample_pairs)[4] <- "Season1"
colnames(sample_pairs)[5] <- "Group1"
colnames(sample_pairs)[6] <- "Site1"
colnames(sample_pairs)[7] <- "Year1"
sample_pairs <- merge(sample_pairs, samples_metadata[,c(1,3,12,13,15,16)], by.x = "sample2", by.y = "Sample_ID")
colnames(sample_pairs)[8] <- "SID2"
colnames(sample_pairs)[9] <- "Season2"
colnames(sample_pairs)[10] <- "Group2"
colnames(sample_pairs)[11] <- "Site2"
colnames(sample_pairs)[12] <- "Year2"

sample_pairs$SubjectPair[sample_pairs$SID1 == sample_pairs$SID2] <- "Same subject"
sample_pairs$SubjectPair[sample_pairs$SID1 != sample_pairs$SID2] <- "Different subjects"

sample_pairs$Sites[sample_pairs$Site1 == sample_pairs$Site2] <- "Same site"
sample_pairs$Sites[sample_pairs$Site1 != sample_pairs$Site2] <- "Different sites"

sample_pairs$SubTypePair[sample_pairs$Group1 == "C" & sample_pairs$Group2 == "W"] <- "Farmer-Cow"
sample_pairs$SubTypePair[sample_pairs$Group1 == "W" & sample_pairs$Group2 == "C"] <- "Farmer-Cow"
sample_pairs$SubTypePair[sample_pairs$Group1 == "C" & sample_pairs$Group2 == "D"] <- "Office-Cow"
sample_pairs$SubTypePair[sample_pairs$Group1 == "D" & sample_pairs$Group2 == "C"] <- "Office-Cow"

###Read in the inStrain-Compare results
inStrainCompare_all <- read.csv("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d17_inStrain_MAGs/221014_jointISCompare_genomeWide_compare.csv", header = TRUE)
#Two of the cow samples (19-C0205-W001fc01 & 19-C0416-F001fc01) have been clustering with human samples when looking at taxa (shotgun and 16S) and ARGs
#These two samples have likley been misplabeled. Removing these two from further analysis
inStrainCompare_all <- inStrainCompare_all[which(inStrainCompare_all$sample1 != "19-C0205-W001fc01" & inStrainCompare_all$sample1 != "19-C0416-F001fc01" & inStrainCompare_all$sample2 != "19-C0205-W001fc01" & inStrainCompare_all$sample2 != "19-C0416-F001fc01"),]

###Looking at the distribution of popANI values
ISC_metadata <- inStrainCompare_all[,c(1,2,3,8,10)]
ISC_metadata <- merge(ISC_metadata, samples_metadata[,c(1,3,12,13,15,16)], by.x = "sample1", by.y = "Sample_ID")
colnames(ISC_metadata)[6] <- "SID1"
colnames(ISC_metadata)[7] <- "Season1"
colnames(ISC_metadata)[8] <- "Group1"
colnames(ISC_metadata)[9] <- "Site1"
colnames(ISC_metadata)[10] <- "Year1"
ISC_metadata <- merge(ISC_metadata, samples_metadata[,c(1,3,12,13,15,16)], by.x = "sample2", by.y = "Sample_ID")
colnames(ISC_metadata)[11] <- "SID2"
colnames(ISC_metadata)[12] <- "Season2"
colnames(ISC_metadata)[13] <- "Group2"
colnames(ISC_metadata)[14] <- "Site2"
colnames(ISC_metadata)[15] <- "Year2"

ISC_metadata$SubjectPair[ISC_metadata$SID1 == ISC_metadata$SID2] <- "Same subject"
ISC_metadata$SubjectPair[ISC_metadata$SID1 != ISC_metadata$SID2] <- "Different subjects"

ISC_metadata$Sites[ISC_metadata$Site1 == ISC_metadata$Site2] <- "Same site"
ISC_metadata$Sites[ISC_metadata$Site1 != ISC_metadata$Site2] <- "Different sites"

ISC_metadata$popANI <- as.numeric(ISC_metadata$popANI)

compare_means(popANI ~ SubjectPair, data = ISC_metadata[,c(4,16)], method = "wilcox.test")
compare_means(popANI ~ Sites, data = ISC_metadata[,c(4,17)], method = "wilcox.test")

p <- ggplot(ISC_metadata, aes(x = SubjectPair, y = popANI)) + geom_boxplot() + theme_test() + ylim(NA, 1.005) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"))
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/inStrain_popANI_SubjectPairs_v1.svg", p, bg = "transparent", width = 2.5, height = 5)

p <- ggplot(ISC_metadata, aes(x = Sites, y = popANI)) + geom_boxplot() + theme_test() + ylim(NA, 1.005) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size= 0.3, color = "gray", linetype = "dashed"), panel.grid.minor.y = element_blank(), 
        axis.text.y = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"))
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/inStrain_popANI_SitePairs_v1.svg", p, bg = "transparent", width = 2.5, height = 5)

###Looking at instances of cross-subject sharing of strains
ISC_subjectFilt_99999 <- ISC_metadata[which(ISC_metadata$SID1 != ISC_metadata$SID2 & ISC_metadata$popANI >= 0.99999 & ISC_metadata$percent_compared >= 0.5),]
#ISC_subjectFilt_99999$SubTypePair[ISC_subjectFilt_99999$Group1 == "C" & ISC_subjectFilt_99999$Group2 == "C"] <- "Cow-Cow"
#ISC_subjectFilt_99999$SubTypePair[ISC_subjectFilt_99999$Group1 == "W" & ISC_subjectFilt_99999$Group2 == "W"] <- "Farmer-Farmer"
#ISC_subjectFilt_99999$SubTypePair[ISC_subjectFilt_99999$Group1 == "D" & ISC_subjectFilt_99999$Group2 == "D"] <- "Office-Office"
#ISC_subjectFilt_99999$SubTypePair[ISC_subjectFilt_99999$Group1 == "W" & ISC_subjectFilt_99999$Group2 == "D"] <- "Office-Farmer"
#ISC_subjectFilt_99999$SubTypePair[ISC_subjectFilt_99999$Group1 == "D" & ISC_subjectFilt_99999$Group2 == "W"] <- "Office-Farmer"
ISC_subjectFilt_99999$SubTypePair[ISC_subjectFilt_99999$Group1 == "C" & ISC_subjectFilt_99999$Group2 == "W"] <- "Farmer-Cow"
ISC_subjectFilt_99999$SubTypePair[ISC_subjectFilt_99999$Group1 == "W" & ISC_subjectFilt_99999$Group2 == "C"] <- "Farmer-Cow"
ISC_subjectFilt_99999$SubTypePair[ISC_subjectFilt_99999$Group1 == "C" & ISC_subjectFilt_99999$Group2 == "D"] <- "Office-Cow"
ISC_subjectFilt_99999$SubTypePair[ISC_subjectFilt_99999$Group1 == "D" & ISC_subjectFilt_99999$Group2 == "C"] <- "Office-Cow"
ISC_subjectFilt_99999 <- ISC_subjectFilt_99999[which(ISC_subjectFilt_99999$SubTypePair != ""),]
dim(ISC_subjectFilt_99999[which(ISC_subjectFilt_99999$SubTypePair == "Farmer-Cow"),])[1]

ISC_subjectFilt_998 <- ISC_metadata[which(ISC_metadata$SID1 != ISC_metadata$SID2 & ISC_metadata$popANI >= 0.998),]
ISC_subjectFilt_998$SubTypePair[ISC_subjectFilt_998$Group1 == "C" & ISC_subjectFilt_998$Group2 == "W"] <- "Farmer-Cow"
ISC_subjectFilt_998$SubTypePair[ISC_subjectFilt_998$Group1 == "W" & ISC_subjectFilt_998$Group2 == "C"] <- "Farmer-Cow"
ISC_subjectFilt_998$SubTypePair[ISC_subjectFilt_998$Group1 == "C" & ISC_subjectFilt_998$Group2 == "D"] <- "Office-Cow"
ISC_subjectFilt_998$SubTypePair[ISC_subjectFilt_998$Group1 == "D" & ISC_subjectFilt_998$Group2 == "C"] <- "Office-Cow"
ISC_subjectFilt_998 <- ISC_subjectFilt_998[which(ISC_subjectFilt_998$SubTypePair != ""),]
dim(ISC_subjectFilt_998[which(ISC_subjectFilt_998$SubTypePair == "Farmer-Cow"),])[1]

ISC_subjectFilt_995 <- ISC_metadata[which(ISC_metadata$SID1 != ISC_metadata$SID2 & ISC_metadata$popANI >= 0.995 & ISC_metadata$percent_compared >= 0.5),]
ISC_subjectFilt_995$SubTypePair[ISC_subjectFilt_995$Group1 == "C" & ISC_subjectFilt_995$Group2 == "W"] <- "Farmer-Cow"
ISC_subjectFilt_995$SubTypePair[ISC_subjectFilt_995$Group1 == "W" & ISC_subjectFilt_995$Group2 == "C"] <- "Farmer-Cow"
ISC_subjectFilt_995$SubTypePair[ISC_subjectFilt_995$Group1 == "C" & ISC_subjectFilt_995$Group2 == "D"] <- "Office-Cow"
ISC_subjectFilt_995$SubTypePair[ISC_subjectFilt_995$Group1 == "D" & ISC_subjectFilt_995$Group2 == "C"] <- "Office-Cow"
ISC_subjectFilt_995 <- ISC_subjectFilt_995[which(ISC_subjectFilt_995$SubTypePair != ""),]
#ISC_subjectFilt_995 <- ISC_subjectFilt_995[which(ISC_subjectFilt_995$Sites == "Same site"),]
#write.csv(ISC_subjectFilt_995, file = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d17_inStrain_MAGs/221118_jointISCompare_genomeWide_compare_ANI995PercComp50filter.csv", row.names = FALSE)
dim(ISC_subjectFilt_995[which(ISC_subjectFilt_995$SubTypePair == "Office-Cow"),])[1]
dim(sample_pairs[which(sample_pairs$SubTypePair == "Farmer-Cow"),])[1]
unique(ISC_subjectFilt_995$genome)
dim(ISC_subjectFilt_995[which(ISC_subjectFilt_995$SubTypePair == "Office-Cow" & ISC_subjectFilt_995$Sites == "Different sites"),])[1]
dim(sample_pairs[which(sample_pairs$SubTypePair == "Office-Cow" & sample_pairs$Sites == "Different sites"),])[1]

ggplot(ISC_subjectFilt_99999, aes(x = SubTypePair, y = popANI)) + geom_boxplot()

###Looking at instances of within-site strain sharing between farmers and cows
ISC_SameFarm_99999 <- ISC_metadata[which(ISC_metadata$popANI >= 0.99999 & ISC_metadata$percent_compared >= 0.5 & ISC_metadata$Sites == "Same site"),]
ISC_SameFarm_99999$SubTypePair[ISC_SameFarm_99999$Group1 == "C" & ISC_SameFarm_99999$Group2 == "W"] <- "Farmer-Cow"
ISC_SameFarm_99999$SubTypePair[ISC_SameFarm_99999$Group1 == "W" & ISC_SameFarm_99999$Group2 == "C"] <- "Farmer-Cow"
ISC_SameFarm_99999$SubTypePair[ISC_SameFarm_99999$Group1 == "C" & ISC_SameFarm_99999$Group2 == "C"] <- "Cow-Cow"
ISC_SameFarm_99999$SubTypePair[ISC_SameFarm_99999$Group1 == "W" & ISC_SameFarm_99999$Group2 == "W"] <- "Farmer-Farmer"
ISC_SameFarm_99999 <- ISC_SameFarm_99999[which(ISC_SameFarm_99999$SubTypePair != ""),]

###Plotting the results of permutation tests for the enrichment of lineage-sharing events
farm_office_SiteInd <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d17_inStrain_MAGs/Permutations/230124_Farmer-Office_SiteInd_FarmerCounts.txt", header = FALSE)
ggplot() + geom_histogram(aes(x = as.integer(farm_office_SiteInd[3:10002,1])), bins = 26, fill = "#18678d") +
  geom_vline(xintercept = 410, color = "black", size = 2) + theme_linedraw() + xlab("Number of Farmer-Cow lineage-sharing events") +
  labs(title = "Permutation (10,000) of the farmer-cow lineage-sharing events indepentdent of site", subtitle = "z-score = 2.87. P = 0.0027 (BH adjusted)")
p <- ggplot() + geom_histogram(aes(x = as.integer(farm_office_SiteInd[3:10002,1])), bins = 26, fill = "#18678d") +
  geom_vline(xintercept = 410, color = "black", size = 2) + theme_linedraw() +
  theme(axis.text = element_blank(), axis.title = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/inStrain_Farmer-Office_SiteInd_v1.svg", p, bg = "transparent", width = 5, height = 5)

farm_office_SameSite <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d17_inStrain_MAGs/Permutations/230124_Farmer-Office_SameSite_FarmerCounts.txt", header = FALSE)
ggplot() + geom_histogram(aes(x = as.integer(farm_office_SameSite[3:10002,1])), bins = 30, fill = "#18678d") + 
  geom_vline(xintercept = 46, color = "black", size = 2) + theme_linedraw() + xlab("Number of Farmer-Cow lineage-sharing events") +
  labs(title = "Permutation (10,000) of the farmer-cow lineage-sharing events within the same sites", subtitle = "z-score = 3.74. P = 0.00018")
p <- ggplot() + geom_histogram(aes(x = as.integer(farm_office_SameSite[3:10002,1])), bins = 30, fill = "#18678d") + 
  geom_vline(xintercept = 46, color = "black", size = 2) + theme_linedraw() +
  theme(axis.title = element_blank(), axis.text = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/inStrain_Farmer-Office_SameSite_v1.svg", p, bg = "transparent", width = 5, height = 5)

farm_office_DiffSite <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d17_inStrain_MAGs/Permutations/230124_Farmer-Office_DiffSite_FarmerCounts.txt", header = FALSE)
ggplot() + geom_histogram(aes(x = as.integer(farm_office_DiffSite[3:10002,1])), bins = 31) + geom_vline(xintercept = 364, color = "red", size = 1.5) + theme_bw() + xlab("Number of Farmer-Cow lineage-sharing events") + 
  labs(title = "Permutation (10,000) of the farmer-cow lineage-sharing events within different sites", subtitle = "z-score = 1.90. P = 0.029")
p <- ggplot() + geom_histogram(aes(x = as.integer(farm_office_DiffSite[3:10002,1])), bins = 31, fill = "#18678d") + 
  geom_vline(xintercept = 364, color = "black", size = 2) + theme_linedraw() +
  theme(axis.title = element_blank(), axis.text = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/inStrain_Farmer-Office_DiffSite_v1.svg", p, bg = "transparent", width = 5, height = 5)

farm_farm_Sites <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d17_inStrain_MAGs/Permutations/230124_Farmer-Farmer_Sites_SameSiteCounts.txt", header = FALSE)
ggplot() + geom_histogram(aes(x = as.integer(farm_farm_Sites[3:10002,1])), bins = 23) + geom_vline(xintercept = 46, color = "red", size = 1.5) + theme_bw() + xlab("Number of same-site Farmer-Cow lineage-sharing events") + 
  labs(title = "Permutation (10,000) of same-site farmer-cow lineage-sharing events", subtitle = "z-score = 9.49. P = 4.6e-21")
p <- ggplot() + geom_histogram(aes(x = as.integer(farm_farm_Sites[3:10002,1])), bins = 23, fill = "#18678d") + 
  geom_vline(xintercept = 46, color = "black", size = 2) + theme_bw() + theme_linedraw() +
  theme(axis.title = element_blank(), axis.text = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/inStrain_Farmer-Farmer_Sites_v1.svg", p, bg = "transparent", width = 5, height = 5)

pnorm(2.8690451664290686, lower.tail = FALSE) #SiteInd
pnorm(3.740881235432306, lower.tail = FALSE) #SameSite
pnorm(1.9027242060600873, lower.tail = FALSE) #DiffSite
pnorm(9.490562843444174, lower.tail = FALSE) #Sites

p.adjust(c(0.002058565, 9.168808e-05, 0.02853827, 1.148951e-21), method = "BH")

###Plotting the results of permutation tests for the enrichment of human-cow pairs or human subject where shared lineages were identified
farm_office_SubjectPairs <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d17_inStrain_MAGs/Permutations/231008_Farmer-Office_SubjectPairs_FarmerCowPairCounts.txt", header = FALSE)
ggplot() + geom_histogram(aes(x = as.integer(farm_office_SubjectPairs[3:10002,1])), bins = 25, fill = "#18678d") + 
  geom_vline(xintercept = 204, color = "black", size = 2) + theme_linedraw() +
  theme(axis.title = element_blank(), axis.text = element_text())
p <- ggplot() + geom_histogram(aes(x = as.integer(farm_office_SubjectPairs[3:10002,1])), bins = 25, fill = "#18678d") + 
  geom_vline(xintercept = 204, color = "black", size = 2) + theme_linedraw() +
  theme(axis.title = element_blank(), axis.text = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/inStrain_Farmer-Office_SubjectPairs_v1.svg", p, bg = "transparent", width = 5, height = 5)

farm_office_Subjects <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d17_inStrain_MAGs/Permutations/231008_Farmer-Office_Subjects_FarmerCounts.txt", header = FALSE)
ggplot() + geom_histogram(aes(x = as.integer(farm_office_Subjects[3:10002,1])), bins = 25, fill = "#18678d") + 
  geom_vline(xintercept = 27, color = "black", size = 2) + theme_linedraw() +
  theme(axis.title = element_blank(), axis.text = element_text())
p <- ggplot() + geom_histogram(aes(x = as.integer(farm_office_Subjects[3:10002,1])), bins = 25, fill = "#18678d") + 
  geom_vline(xintercept = 27, color = "black", size = 2) + theme_linedraw() +
  theme(axis.title = element_blank(), axis.text = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Supplements/inStrain_Farmer-Office_Subjects_v1.svg", p, bg = "transparent", width = 5, height = 5)

pnorm(2.070396111496948, lower.tail = FALSE)
pnorm(1.255582231816633, lower.tail = FALSE)

p.adjust(c(0.01920763, 0.1046337), method = "BH")


