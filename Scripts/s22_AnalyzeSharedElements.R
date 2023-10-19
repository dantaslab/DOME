########################################################################################
###This script analyzes the fragment sharing events between DOME metagenomes and MAGs###
########################################################################################

###Load the necessary packages
library(reshape2)
library(ggplot2)
library(circlize)
library(dplyr)
library(NatParksPalettes)

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
samples_metadata <- samples_metadata[which(samples_metadata$Sample_ID != "19-C0205-W001fc01" & samples_metadata$Sample_ID != "19-C0416-F001fc01"),]

###Read in the file with counts of sharead elements between samples
sharedElements_count <- read_delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d22_BLASTn_samples/BLASTn/221004_ShareElementsCounts_condensed.txt")

###Look at the number of shared elements between different subject types
sharedElements_filtered <- merge(sharedElements_count, samples_metadata[,c(1,3,12,13,15,16)], by.x = "Sample1", by.y = "Sample_ID")
colnames(sharedElements_filtered)[4] <- "Subject1"
colnames(sharedElements_filtered)[5] <- "Season1"
colnames(sharedElements_filtered)[6] <- "Group1"
colnames(sharedElements_filtered)[7] <- "Site1"
colnames(sharedElements_filtered)[8] <- "Year1"
sharedElements_filtered <- merge(sharedElements_filtered, samples_metadata[,c(1,3,12,13,15,16)], by.x = "Sample2", by.y = "Sample_ID")
colnames(sharedElements_filtered)[9] <- "Subject2"
colnames(sharedElements_filtered)[10] <- "Season2"
colnames(sharedElements_filtered)[11] <- "Group2"
colnames(sharedElements_filtered)[12] <- "Site2"
colnames(sharedElements_filtered)[13] <- "Year2"
sharedElements_filtered <- sharedElements_filtered[which(sharedElements_filtered$Subject1 != sharedElements_filtered$Subject2),]
sharedElements_filtered$Groups <- ""

for (i in seq(1, dim(sharedElements_filtered)[1])) {
  sharedElements_filtered$Groups[i] <- paste(sort(c(sharedElements_filtered$Group1[i], sharedElements_filtered$Group2[i]))[1], sort(c(sharedElements_filtered$Group1[i], sharedElements_filtered$Group2[i]))[2])
}

ggplot(sharedElements_filtered, aes(x = Groups, y = log(Count))) + geom_boxplot() +
  theme_test() + ylab("log( number of shared DNA segments between sample pairs )")

sharedElements_CowHuman <- sharedElements_filtered[which(sharedElements_filtered$Groups == "C D" | sharedElements_filtered$Groups == "C W"),]
sharedElements_CowHuman$Sites[sharedElements_CowHuman$Site1 == sharedElements_CowHuman$Site2] <- "Same site"
sharedElements_CowHuman$Sites[sharedElements_CowHuman$Site1 != sharedElements_CowHuman$Site2] <- "Different sites"
sharedElements_CowHuman <- sharedElements_CowHuman[which(sharedElements_CowHuman$Season1 == sharedElements_CowHuman$Season2 & sharedElements_CowHuman$Year1 == sharedElements_CowHuman$Year2),]
sharedElements_CowHuman$SeasonFull <- paste(sharedElements_CowHuman$Season1, sharedElements_CowHuman$Year1)
sharedElements_CowHuman$HumanGroup <- ""

for (i in seq(1, dim(sharedElements_CowHuman)[1])) {
 sharedElements_CowHuman$HumanGroup[i] <- sort(c(sharedElements_CowHuman$Group1[i], sharedElements_CowHuman$Group2[i]))[2] 
}

sharedElements_CowHuman$HumanSID <- ""

for (i in seq(1, dim(sharedElements_CowHuman)[1])) {
  if (sharedElements_CowHuman$Group1[i] == sharedElements_CowHuman$HumanGroup[i]) {
    sharedElements_CowHuman$HumanSID[i] <- sharedElements_CowHuman$Subject1[i]
  } else {
    sharedElements_CowHuman$HumanSID[i] <- sharedElements_CowHuman$Subject2[i]
  }
}

sharedElements_CowHuman <- sharedElements_CowHuman[,c(15:18,3)]

sharedElements_CowHuman <- sharedElements_CowHuman %>%
  group_by(HumanSID, HumanGroup, Sites, SeasonFull) %>%
  summarise(
    Count = mean(Count)
  )

sharedElements_CowHuman$GroupSites <- paste(sharedElements_CowHuman$HumanGroup, sharedElements_CowHuman$Sites)
sharedElements_CowHuman$SeasonFull <- factor(sharedElements_CowHuman$SeasonFull, levels = c("SPRING 1", "SUMMER 1", "FALL 1", "WINTER 1", "SPRING 2"))

ggplot(data = sharedElements_CowHuman, aes(x = Sites, y = log(Count), color = HumanGroup)) + geom_boxplot() + 
  facet_wrap(~sharedElements_CowHuman$SeasonFull, nrow = 1, labeller = season_labeller) + 
  theme_test() + ylab("log( number of shared DNA segments between sample pairs )") +
  scale_color_manual(values=c("#00BA38", "#619CFF"), name = "", labels = c("Office worker", "Farmer"))


###Read in the file with shared elements between MAGs
MAGs_sharedElements <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d22_BLASTn_HQMQ/joint_BLASTn_allHQMQMAGs_nodup_annot_filt.txt")
dim(MAGs_sharedElements[which(MAGs_sharedElements$MAGqual1 == "HQ" & MAGs_sharedElements$MAGqual2 == "HQ" & MAGs_sharedElements$HostPair == "C-W"),])
dim(MAGs_sharedElements[which(MAGs_sharedElements$HostPair == "C-D"),])

###Creating a chord diagram to show the sharing of MGEs between cow and farmer MAGs
CW_hits_top10families <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d22_BLASTn_HQMQ/ChordDiagram/CW_hitsBLASTn_top10families_ARG.txt")
CW_hits_top10families <- CW_hits_top10families %>% group_by(W_family_alt, C_family_alt, W_family, C_family, Families, ARG) %>% summarise_all(sum)
CW_hits_top10families$FamARG <- paste(CW_hits_top10families$Families, CW_hits_top10families$ARG)
CW_hits_top10families$zIndex <- ""
CW_hits_top10families$zIndex[CW_hits_top10families$ARG == "Present"] <- 1
CW_hits_top10families$zIndex[CW_hits_top10families$ARG == "None"] <- 2

colors_FamARG <- data_frame(FamARG = c("Same family None", "Same family Present", "Different families None", "Different families Present"),
                            Color = c("#75B9FA", "#1A3D82", "#E49362", "#832B0F"),
                            Transparency = c(0.7, 0, 0.7, 0))

CW_hits_top10families <- merge(CW_hits_top10families, colors_FamARG, by = "FamARG")

colors_families <- c("W-f__Lachnospiraceae" = "#082544", "W-f__Bifidobacteriaceae" = "#11395D", "W-f__Bacteroidaceae" = "#1B4E76", "W-f__Acutalibacteraceae" = "#3C5A82", "W-f__Ruminococcaceae" = "#646288", "W-f__Eggerthellaceae" = "#8F5C79", "W-f__UBA1381" = "#BC4A53", "W-f__Treponemataceae" = "#E04D3E", "W-f__Turicibacteraceae" = "#E9945E", "W-f__Erysipelatoclostridiaceae" = "#F2DC7E", "W-Other" = "#bcbcbc",
                     "C-f__Lachnospiraceae" = "#082544", "C-f__Bifidobacteriaceae" = "#11395D", "C-f__Bacteroidaceae" = "#1B4E76", "C-f__Acutalibacteraceae" = "#3C5A82", "C-f__Ruminococcaceae" = "#646288", "C-f__Eggerthellaceae" = "#8F5C79", "C-f__UBA1381" = "#BC4A53", "C-f__Treponemataceae" = "#E04D3E", "C-f__Turicibacteraceae" = "#E9945E", "C-f__Erysipelatoclostridiaceae" = "#F2DC7E", "C-Other" = "#bcbcbc")

colors_families <- c("W-f__Lachnospiraceae" = "#082544", "W-f__Bifidobacteriaceae" = "#154268", "W-f__Bacteroidaceae" = "#355880", "W-f__Acutalibacteraceae" = "#6F648A", "W-f__Ruminococcaceae" = "#985871", "W-f__Eggerthellaceae" = "#C3474D", "W-f__UBA1381" = "#E15542", "W-f__Treponemataceae" = "#E99860", "W-f__Turicibacteraceae" = "#EDBA6F", "W-f__Erysipelatoclostridiaceae" = "#F2DC7E", "W-Other" = "#bcbcbc",
                     "C-f__Lachnospiraceae" = "#082544", "C-f__Bifidobacteriaceae" = "#11395D", "C-f__Bacteroidaceae" = "#1B4E76", "C-f__Acutalibacteraceae" = "#3C5A82", "C-f__Ruminococcaceae" = "#646288", "C-f__Eggerthellaceae" = "#8F5C79", "C-f__UBA1381" = "#BC4A53", "C-f__Treponemataceae" = "#E04D3E", "C-f__Turicibacteraceae" = "#E9945E", "C-f__Erysipelatoclostridiaceae" = "#F2DC7E", "C-Other" = "#bcbcbc")

colors_families <- c("W-f__Lachnospiraceae" = "#082544", "W-f__Bifidobacteriaceae" = "#154268", "W-f__Bacteroidaceae" = "#355880", "W-f__Acutalibacteraceae" = "#6F648A", "W-f__Ruminococcaceae" = "#985871", "W-f__Eggerthellaceae" = "#C3474D", "W-f__UBA1381" = "#E15542", "W-f__Treponemataceae" = "#FD8700", "W-f__Turicibacteraceae" = "#FEB424", "W-f__Erysipelatoclostridiaceae" = "#F2DC7E", "W-Other" = "#bcbcbc",
                     "C-f__Lachnospiraceae" = "#082544", "C-f__Bifidobacteriaceae" = "#154268", "C-f__Bacteroidaceae" = "#355880", "C-f__Acutalibacteraceae" = "#6F648A", "C-f__Ruminococcaceae" = "#985871", "C-f__Eggerthellaceae" = "#C3474D", "C-f__UBA1381" = "#E15542", "C-f__Treponemataceae" = "#FD8700", "C-f__Turicibacteraceae" = "#FEB424", "C-f__Erysipelatoclostridiaceae" = "#F2DC7E", "C-Other" = "#bcbcbc")
#colors_families <- c("W-f__Acutalibacteraceae" = "#082544", "W-f__Bacteroidaceae" = "#11395D", "W-f__Bifidobacteriaceae" = "#1B4E76", "W-f__Eggerthellaceae" = "#3C5A82", "W-f__Erysipelatoclostridiaceae" = "#646288", "W-f__Lachnospiraceae" = "#8F5C79", "W-f__Ruminococcaceae" = "#BC4A53", "W-f__Treponemataceae" = "#E04D3E", "W-f__Turicibacteraceae" = "#E9945E", "W-f__UBA1381" = "#F2DC7E", "W-Other" = "#bcbcbc",
#                     "C-f__Acutalibacteraceae" = "#082544", "C-f__Bacteroidaceae" = "#11395D", "C-f__Bifidobacteriaceae" = "#1B4E76", "C-f__Eggerthellaceae" = "#3C5A82", "C-f__Erysipelatoclostridiaceae" = "#646288", "C-f__Lachnospiraceae" = "#8F5C79", "C-f__Ruminococcaceae" = "#BC4A53", "C-f__Treponemataceae" = "#E04D3E", "C-f__Turicibacteraceae" = "#E9945E", "C-f__UBA1381" = "#F2DC7E", "C-Other" = "#bcbcbc")

chordDiagram(CW_hits_top10families[c(2,3,8)], grid.col = colors_families, col = CW_hits_top10families$Color, transparency = CW_hits_top10families$Transparency,
             link.zindex = CW_hits_top10families$zIndex,
             order = c("W-f__Lachnospiraceae", "W-f__Bifidobacteriaceae", "W-f__Bacteroidaceae", "W-f__Acutalibacteraceae", "W-f__Ruminococcaceae", "W-f__Eggerthellaceae", "W-f__UBA1381", "W-f__Treponemataceae", "W-f__Turicibacteraceae", "W-f__Erysipelatoclostridiaceae", "W-Other",
                       "C-Other", "C-f__Erysipelatoclostridiaceae", "C-f__Turicibacteraceae", "C-f__Treponemataceae", "C-f__UBA1381", "C-f__Eggerthellaceae", "C-f__Ruminococcaceae", "C-f__Acutalibacteraceae", "C-f__Bacteroidaceae", "C-f__Bifidobacteriaceae", "C-f__Lachnospiraceae"),
             annotationTrack = c("axis", "grid"))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 1)
}, bg.border = NA)

svg("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/SharedFragments_ChordDiagram_v1.svg")
chordDiagram(CW_hits_top10families[c(2,3,8)], grid.col = colors_families, col = CW_hits_top10families$Color, transparency = CW_hits_top10families$Transparency,
             link.zindex = CW_hits_top10families$zIndex,
             order = c("W-Other", "W-f__Erysipelatoclostridiaceae", "W-f__Turicibacteraceae", "W-f__Treponemataceae", "W-f__UBA1381", "W-f__Eggerthellaceae", "W-f__Ruminococcaceae", "W-f__Acutalibacteraceae", "W-f__Bacteroidaceae", "W-f__Bifidobacteriaceae", "W-f__Lachnospiraceae",
                       "C-f__Lachnospiraceae", "C-f__Bifidobacteriaceae", "C-f__Bacteroidaceae", "C-f__Acutalibacteraceae", "C-f__Ruminococcaceae", "C-f__Eggerthellaceae", "C-f__UBA1381", "C-f__Treponemataceae", "C-f__Turicibacteraceae", "C-f__Erysipelatoclostridiaceae", "C-Other"),
             annotationTrack = c("axis", "grid"))
dev.off()

###Creating a chord diagram to show the sharing of MGE-encoded ARGs between cow and farmer MAGs
CW_ARGhits_top10fam <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d22_BLASTn_HQMQ/ChordDiagram/CW_hitsBLASTn_top10families_onlyARG.txt")
CW_ARGhits_top10fam <- CW_ARGhits_top10fam %>% group_by(W_family_alt, C_family_alt, W_family, C_family, Families, ARG) %>% summarise_all(sum)

colors_ARGclasses <- data_frame(ARGclass = c("TETRACYCLINE", "TRIMETHOPRIM", "GLYCOPEPTIDE", "BETA.LACTAM", "LINCOSAMIDE", "MACROLIDE", "PHENICOL", "QUINOLONE", "AMINOGLYCOSIDE", "FOSFOMYCIN", "Other"),
                                Color = c(c(natparks.pals("DeathValley", 11))[1:10], "#bcbcbc"))

colors_ARGclasses <- data_frame(ARGclass = c("TETRACYCLINE", "TRIMETHOPRIM", "GLYCOPEPTIDE", "BETA.LACTAM", "LINCOSAMIDE", "MACROLIDE", "PHENICOL", "QUINOLONE", "AMINOGLYCOSIDE", "FOSFOMYCIN", "Other"),
                                Color = c("#8C2B0E", "#B25422", "#D8813B", "#FEB359", "#9F7E59", "#233F6C", "#435F90", "#5B4C64", "#81565F", "#B47E83", "#848484"))

CW_ARGhits_top10fam <- merge(CW_ARGhits_top10fam, colors_ARGclasses, by.x = "ARG", by.y = "ARGclass")

chordDiagram(CW_ARGhits_top10fam[c(2,3,7)], grid.col = colors_families, col = CW_ARGhits_top10fam$Color, transparency = 0.3,
             order = c("W-f__Lachnospiraceae", "W-f__Bifidobacteriaceae", "W-f__Bacteroidaceae", "W-f__Acutalibacteraceae", "W-f__Ruminococcaceae", "W-f__Erysipelatoclostridiaceae", "W-Other",
                       "C-Other", "C-f__Erysipelatoclostridiaceae", "C-f__Turicibacteraceae", "C-f__Treponemataceae", "C-f__UBA1381", "C-f__Acutalibacteraceae", "C-f__Bacteroidaceae", "C-f__Bifidobacteriaceae", "C-f__Lachnospiraceae"),
             annotationTrack = c("axis", "grid"))

unique(CW_ARGhits_top10fam$C_family_alt)
c(c(natparks.pals("DeathValley", 11))[1:10],"TEST")
c(natparks.pals("DeathValley", 11))
natparks.pals("DeathValley", 11)

CW_top3ARGhits_top10fam <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d22_BLASTn_HQMQ/ChordDiagram/CW_hitsBLASTn_top10families_onlyARGtop3.txt")
CW_top3ARGhits_top10fam <- CW_top3ARGhits_top10fam %>% group_by(W_family_alt, C_family_alt, W_family, C_family, Families, ARG) %>% summarise_all(sum)
CW_top3ARGhits_top10fam <- merge(CW_top3ARGhits_top10fam, colors_ARGclasses, by.x = "ARG", by.y = "ARGclass")

chordDiagram(CW_top3ARGhits_top10fam[c(2,3,7)], grid.col = colors_families, col = CW_top3ARGhits_top10fam$Color, transparency = 0.3,
             order = c("W-f__Lachnospiraceae", "W-f__Bifidobacteriaceae", "W-f__Bacteroidaceae", "W-f__Acutalibacteraceae", "W-f__Ruminococcaceae", "W-f__Erysipelatoclostridiaceae", "W-Other",
                       "C-Other", "C-f__Erysipelatoclostridiaceae", "C-f__Turicibacteraceae", "C-f__Treponemataceae", "C-f__UBA1381", "C-f__Acutalibacteraceae", "C-f__Bacteroidaceae", "C-f__Bifidobacteriaceae", "C-f__Lachnospiraceae"),
             annotationTrack = c("axis", "grid"))

svg("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/SharedFragments_ChordDiagram_top3ARGs_v1.svg")
chordDiagram(CW_top3ARGhits_top10fam[c(2,3,7)], grid.col = colors_families, col = CW_top3ARGhits_top10fam$Color, transparency = 0.3,
             order = c("W-Other", "W-f__Erysipelatoclostridiaceae", "W-f__Ruminococcaceae", "W-f__Acutalibacteraceae", "W-f__Bacteroidaceae", "W-f__Bifidobacteriaceae", "W-f__Lachnospiraceae",
                       "C-f__Lachnospiraceae", "C-f__Bifidobacteriaceae", "C-f__Bacteroidaceae", "C-f__Acutalibacteraceae", "C-f__UBA1381", "C-f__Treponemataceae",  "C-f__Turicibacteraceae", "C-f__Erysipelatoclostridiaceae", "C-Other"),
             annotationTrack = c("axis", "grid"))
dev.off()

CW_top3CowEnrichedARGhits_top10fam <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d22_BLASTn_HQMQ/ChordDiagram/CW_hitsBLASTn_top10families_onlyCowEnrichedARGtop3.txt")
CW_top3CowEnrichedARGhits_top10fam <- CW_top3CowEnrichedARGhits_top10fam %>% group_by(W_family_alt, C_family_alt, W_family, C_family, Families, ARG) %>% summarise_all(sum)
CW_top3CowEnrichedARGhits_top10fam <- merge(CW_top3CowEnrichedARGhits_top10fam, colors_ARGclasses, by.x = "ARG", by.y = "ARGclass")

chordDiagram(CW_top3CowEnrichedARGhits_top10fam[c(2,3,7)], grid.col = colors_families, col = CW_top3CowEnrichedARGhits_top10fam$Color, transparency = 0.3,
             order = c("W-f__Lachnospiraceae", "W-f__Bifidobacteriaceae", "W-f__Bacteroidaceae", "W-f__Acutalibacteraceae", "W-f__Ruminococcaceae", "W-f__Erysipelatoclostridiaceae", "W-Other",
                       "C-Other", "C-f__Erysipelatoclostridiaceae", "C-f__Turicibacteraceae", "C-f__Treponemataceae", "C-f__UBA1381", "C-f__Acutalibacteraceae", "C-f__Bacteroidaceae", "C-f__Bifidobacteriaceae", "C-f__Lachnospiraceae"),
             annotationTrack = c("axis", "grid"))

svg("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/SharedFragments_ChordDiagram_top3CowEnrichedARGs_v1.svg")
chordDiagram(CW_top3CowEnrichedARGhits_top10fam[c(2,3,7)], grid.col = colors_families, col = CW_top3CowEnrichedARGhits_top10fam$Color, transparency = 0.3,
             order = c("W-Other", "W-f__Erysipelatoclostridiaceae", "W-f__Ruminococcaceae", "W-f__Acutalibacteraceae", "W-f__Bacteroidaceae", "W-f__Bifidobacteriaceae", "W-f__Lachnospiraceae",
                       "C-f__Lachnospiraceae", "C-f__Bifidobacteriaceae", "C-f__Bacteroidaceae", "C-f__Acutalibacteraceae", "C-f__UBA1381", "C-f__Treponemataceae", "C-f__Erysipelatoclostridiaceae", "C-Other"),
             annotationTrack = c("axis", "grid"))
dev.off()
