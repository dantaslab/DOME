###This script is used to visualize the distribution of the PlasClass plasmid probabilities for all assembled contigs

library(ggplot2)

data = read.delim("/Users/bejanmahmud/Library/CloudStorage/Box-Box/Bejan's files/Projects/DOME/DOME/Analysis/Metagenome/220221_DOME_fc_all/d11_plasclass/CombinedOutput.txt", header = FALSE)

colnames(data)[1] <- "Contig"
colnames(data)[2] <- "Probability"

ggplot(data, aes(Probability)) + geom_histogram(binwidth = 0.01)