###This script plots the identities of the FMG hits with the NCBI and CARD/AMRProt references

library(ggplot2)
library(ggExtra)
library(reshape2)

identities <- read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d19_needle/220818_FMG-NCBI-ARG_identities.txt")
identities.m <- melt(identities, id.vars = "Source", measure.vars = c("NCBI", "ARG"))

identities$Source <- factor(identities$Source, levels = c("C", "W", "D"))

ggplot(data = identities, aes(x = ARG, y = NCBI, color = Source)) + geom_point(alpha = 0.7, size = 3) + ylim(0,100) +
  theme_test() + 
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("#873e23", "#18678d", "#626262"))
p <- ggplot(data = identities, aes(x = ARG, y = NCBI, color = Source)) + geom_point(alpha = 0.5, size = 3) + ylim(0,100) +
  theme_test() + 
  theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), panel.border = element_rect(linewidth = 0.75, color = "black")) +
  scale_color_manual(values = c("#873e23", "#18678d", "#626262"))
p <- ggMarginal(p, groupColour = TRUE, groupFill = TRUE, type = "boxplot")
ggsave("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Manuscript/Figures/FMG_CARD-NCBI_Identities_v1.svg", p, bg = "transparent", width = 10, height = 10)

ggplot(data = identities.m, aes(x = variable, y = value, color = Source)) + geom_boxplot() + theme_test() + ylab("Percent amino acid identity") + xlab("Database")

test <- compare_means(value ~ Source, data = identities.m, group.by = "variable", method = "wilcox.test", p.adjust.method = "BH")
test <- test[which(test$p.adj < 0.05),]

median(identities$ARG)
median(identities$NCBI)

median(identities[which(identities$Source == "C"), c("NCBI")])
median(identities[which(identities$Source == "W"), c("NCBI")])
median(identities[which(identities$Source == "D"), c("NCBI")])
