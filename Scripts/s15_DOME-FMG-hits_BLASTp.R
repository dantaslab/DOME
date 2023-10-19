###This script is to plot the 2d density plot of the blast results of DOME FMG hits against NCBI AMRProt and CARD references to show varibility
###Also to plot the taxonomic origins of the top NCBI nr database hits of FMG hits

library(ggplot2)
library(dplyr)

blast_sum <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d15_BLAST/DOME_FMG-NCBI_CARD-blastp-BLOSUM45_summary.txt",
                        header = TRUE)

mapping <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d15_BLAST/FMGhits_mapping_v1.txt",
                      header = TRUE)

FMG_origins <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d15_BLAST/FMG_RefSeq_origins.txt",
                          header = TRUE)

blast_sum <- merge(blast_sum, mapping, by.x = "Hit", by.y = "Hit")

ggplot(data = blast_sum, aes(x = Coverage, y = Identity)) +
  geom_hex(bins = 100) + scale_fill_continuous(type = "viridis") + theme_test() + xlab("Coverage (%)") + ylab("Identity (%)")

ggplot(data = blast_sum[which(blast_sum$Origin == "C"),], aes(x = Coverage, y = Identity)) +
  geom_hex(bins = 100) + scale_fill_continuous(type = "viridis") + theme_test() + xlab("Coverage (%)") + ylab("Identity (%)")

ggplot(data = blast_sum[which(blast_sum$Origin == "W"),], aes(x = Coverage, y = Identity)) +
  geom_hex(bins = 100) + scale_fill_continuous(type = "viridis") + theme_test() + xlab("Coverage (%)") + ylab("Identity (%)")

ggplot(data = blast_sum[which(blast_sum$Origin == "D"),], aes(x = Coverage, y = Identity)) +
  geom_hex(bins = 100) + scale_fill_continuous(type = "viridis") + theme_test() + xlab("Coverage (%)") + ylab("Identity (%)")

ggplot(data = blast_sum, aes(x = Coverage, y = Identity)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) + scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position='none'
  )

ggplot(data = blast_sum, aes(x = Coverage, y = Identity, color = Origin)) + geom_point() + geom_jitter()

ggplot(data = blast_sum, aes(x = Origin, y = Coverage)) + geom_boxplot()

FMG_origins <- FMG_origins[,c(3,6)]
FMG_origins$Count <- 1

top_taxa <- FMG_origins %>% select(RefSeq.origin, Count) %>% group_by(RefSeq.origin) %>% summarise(Count = sum(Count))
top_taxa <- top_taxa[order(-top_taxa$Count),]
top_taxa <- top_taxa[c(1:9, 11, 12, 13, 15, 16, 18:24),]

FMG_origins <- FMG_origins %>% group_by(Source, RefSeq.origin) %>% summarise(Count = sum(Count))
FMG_origins$TopTaxa <- ""
for (i in seq(1, 496, 1)) {
  if (FMG_origins$RefSeq.origin[i] %in% top_taxa$RefSeq.origin) {
    FMG_origins$TopTaxa[i] <- FMG_origins$RefSeq.origin[i]
  } else {
    FMG_origins$TopTaxa[i] <- "Others"
  }
}

FMG_origins$TopTaxa <- factor(FMG_origins$TopTaxa, levels = c("Others", rev(top_taxa$RefSeq.origin)))
FMG_origins$Source <- factor(FMG_origins$Source, levels = c("D", "W", "C"))

ggplot(FMG_origins[which(FMG_origins$TopTaxa != "Others"),], aes(fill = Source, y = Count, x = TopTaxa)) + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.98) + coord_flip() + theme_linedraw() + 
  xlab("Best BLAST hit (NCBI nr protein database)") + ylab("Unique resistance-conferring ORFs") +
  scale_fill_manual(values = c("#626262", "#18678d", "#873e23")) +
  theme(panel.border = element_rect(linewidth = 0.75, color = "black"), panel.grid.major.y = element_line(size = 0.15), 
        panel.grid.major.x = element_line(size = 0.10), panel.grid.minor = element_blank())
p <- ggplot(FMG_origins[which(FMG_origins$TopTaxa != "Others"),], aes(fill = Source, y = Count, x = TopTaxa)) + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.98) + coord_flip() + theme_linedraw() + 
  xlab("Best BLAST hit (NCBI nr protein database)") + ylab("Unique resistance-conferring ORFs") +
  scale_fill_manual(values = c("#626262", "#18678d", "#873e23")) +
  theme(panel.border = element_rect(linewidth = 0.75, color = "black"), panel.grid.major.y = element_line(size = 0.15), 
        panel.grid.major.x = element_line(size = 0.10), panel.grid.minor = element_blank(),
        axis.title = element_blank(), legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank())
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/FMG_BLASTp_Sources_v2.svg", p, bg = "transparent", width = 5, height = 10)

