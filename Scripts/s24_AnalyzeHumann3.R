################################################################
###This script analyzes the HUMANN3 profiles of human samples###
################################################################

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


###Read in the Metaphlan4 species output file
pathways <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d24_humann3/norm/Merged_pathabundance.txt")
colnames(pathways)[3] <- "Abundance_CPM"
pathways <- dcast(pathways, Sample ~ Pathway, value.var = "Abundance_CPM")
pathways[is.na(pathways)] <- 0
rownames(pathways) <- pathways$Sample
pathways <- pathways[,2:475]
pathways_filtered <- pathways[,colnames(pathways) %notin% c("UNMAPPED", "UNINTEGRATED")]

###Running Maaslin2 to check for associations of pathways with suject types
maaslin_metadata <- samples_metadata %>% dplyr::select(3, 1, 12:16)
maaslin_metadata$SeasonFull <- paste(maaslin_metadata$Season, maaslin_metadata$Year)
rownames(maaslin_metadata) <- maaslin_metadata$Sample_ID

Maaslin2(
  input_data = pathways_filtered,
  input_metadata = maaslin_metadata[which(maaslin_metadata$Group != "C"),],
  fixed_effects = "Group",
  random_effects = c("SID", "SeasonFull", "Site"),
  output = "/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d24_humann3/Maaslin2/Pathways_farmerVSoffice"
)

M2_pathways_famerVSoffice <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d24_humann3/Maaslin2/Pathways_farmerVSoffice/significant_results.tsv", header = TRUE, "\t")
M2_pathways_famerVSoffice$shape <- ifelse(M2_pathways_famerVSoffice$coef > 0, 1, -1)
M2_pathways_famerVSoffice <- subset(M2_pathways_famerVSoffice, qval < 0.05)
p <- ggplot(M2_pathways_famerVSoffice, aes(x=coef, y=reorder(feature, coef))) +
  geom_point(color = "#18678d", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0), color = "#18678d", size = 0.7) +
  theme_linedraw() +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 20, face = "bold"), axis.text.y = element_text(size = 14)) +
  geom_vline(xintercept=0, colour='black',size=2, linetype="dotted") +
  labs(x="Coefficient (relative to controls)", y="Pathway")
p <- ggplot(M2_pathways_famerVSoffice, aes(x=coef, y=reorder(feature, coef))) +
  geom_point(color = "#18678d", shape = 16, size = 3, alpha = 0.8) +
  geom_errorbarh(aes(xmax=coef + stderr, xmin =coef - stderr, height = 0), color = "#18678d", size = 0.7) +
  theme_linedraw()+
  theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.text.y = element_blank()) +
  geom_vline(xintercept=0, colour='black',size=2, linetype="dotted")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/Shotgun_Maaslin2_FarmerOffice_Pathways_v3.svg", p, bg = "transparent", width = 5, height = 10)
