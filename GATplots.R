#!/usr/bin/env Rscript

# GAT CpG and genic annotation enrichment plots
# Takes the output of GAT for hg38 CpG and genic annotation enrichments and makes a ggplot2
# Ben Laufer

# Install -----------------------------------------------------------------

if (!requireNamespace(c("ggplot2", "reshape2"), quietly = TRUE))
  install.packages(c("ggplot2", "reshape2"), repos = "https://cloud.r-project.org")
library(ggplot2)
library(reshape2)

# Functions ---------------------------------------------------------------

#' gg_color_hue
#' @description Generate ggplot2 style colors
#' @param n Number of samples
#' @return Character string of colors
#' @references \url{https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette}
#' @importFrom glue glue
#' @export gg_color_hue
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# CpG ---------------------------------------------------------------------

# Load
GAT <- read.delim("Gat_CpG_hyper_hypo_results.tsv")

# Nice names
GAT$annotation <- as.character(GAT$annotation)
GAT$annotation[GAT$annotation == "hg38_cpg_inter"] <- "Open Sea"
GAT$annotation[GAT$annotation == "hg38_cpg_shelves"] <- "CpG Shelf"
GAT$annotation[GAT$annotation == "hg38_cpg_islands"] <- "CpG Island"
GAT$annotation[GAT$annotation == "hg38_cpg_shores"] <- "CpG Shore"
GAT$annotation <- as.factor(GAT$annotation)

GAT.1 <- GAT

# Fix Fold Change
index <- GAT$fold < 1
GAT$fold[index] <- -(1/GAT$fold[index])

# Order
GAT <- dcast(GAT, annotation ~ track, value.var = "fold")
order <- c("CpG Island", "CpG Shore", "CpG Shelf", "Open Sea")
GAT$annotation <- factor(GAT$annotation,levels = order)
GAT <- GAT[match(order, GAT$annotation),]
row.names(GAT) <- c(1:nrow(GAT))
GAT$annotation <- factor(GAT$annotation, levels = rev(order))

# Colors for Custom Order
n = nrow(GAT)
cols = gg_color_hue(n)
names(cols) <- order

# Melt
GAT <- melt(GAT)

# Signifcance
GAT.1 <- dcast(GAT.1, annotation ~ track, value.var = "qvalue")
order <- c("CpG Island", "CpG Shore", "CpG Shelf", "Open Sea")
GAT.1$annotation <- factor(GAT.1$annotation,levels = order)
GAT.1 <- GAT.1[match(order, GAT.1$annotation),]
row.names(GAT.1) <- c(1:nrow(GAT.1))
GAT.1$annotation <- factor(GAT.1$annotation, levels = rev(order))
GAT.1 <- melt(GAT.1)
GAT.1$signif <- ifelse(GAT.1$value < 0.05,1,0)
GAT.1 <- GAT.1[-3]
GAT <- merge(GAT, GAT.1, by =c("annotation", "variable"))

# Plot
CpG <- ggplot(data = GAT, aes(annotation, y = value, fill = annotation)) +
  geom_bar(stat = "identity") + coord_flip() + labs(y = "Fold Enrichment", x = element_blank()) + theme_classic() + 
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16),
        strip.text = element_text(size = 16), legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  scale_y_continuous(limits = c(-4,4), breaks = c(-4,-3,-2,-1,0,1,2,3,4)) +
  scale_fill_manual(values = cols, breaks = order, name = "Annotation") + facet_grid(~variable) +
  geom_hline(yintercept = 0) +
  geom_text(data = GAT[(GAT$signif == 1 & GAT$value > 0), ], label = "*", size = 8, show.legend = FALSE, nudge_y = 0.5, nudge_x = -0.05) +
  geom_text(data = GAT[(GAT$signif == 1 & GAT$value < 0), ], label = "*", size = 8, show.legend = FALSE, nudge_y = -0.5, nudge_x = -0.05)

ggsave("Gat_CpG_hyper_hypo.pdf",
       plot = CpG,
       device = NULL,
       width = 11,
       height = 6)

# Genic  ------------------------------------------------------------------

# Load
GAT2 <- read.delim("Gat_genic_hyper_hypo_results.tsv")

# Nice names
GAT2$annotation <- as.character(GAT2$annotation)
GAT2$annotation[GAT2$annotation == "enhancers_merged.bed"] <- "Enhancers"
GAT2$annotation[GAT2$annotation == "intergenic_merged.bed"] <- "Intergenic"
GAT2$annotation[GAT2$annotation == "introns_merged.bed"] <- "Intron"
GAT2$annotation[GAT2$annotation == "onetofivekb_merged.bed"] <- "1-5 kb Upstream"
GAT2$annotation[GAT2$annotation == "boundaries_merged.bed"] <- "Exon/Intron Boundaries"
GAT2$annotation[GAT2$annotation == "threeUTRs_merged.bed"] <- "3' UTR"
GAT2$annotation[GAT2$annotation == "exons_merged.bed"] <- "Exon"
GAT2$annotation[GAT2$annotation == "fiveUTRs_merged.bed"] <- "5' UTR"
GAT2$annotation[GAT2$annotation == "promoters_merged.bed"] <- "Promoter"
GAT2$annotation <- as.factor(GAT2$annotation)

GAT2.1 <- GAT2

# Fix Fold Change
index1 <- GAT2$fold < 1
GAT2$fold[index1] <- -(1/GAT2$fold[index1])

# Order
GAT2 <- dcast(GAT2, annotation ~ track, value.var = "fold")
order <- c("Enhancers", "1-5 kb Upstream", "Promoter", "5' UTR", "Exon", "Exon/Intron Boundaries", "Intron", "3' UTR", "Intergenic")
GAT2$annotation <- factor(GAT2$annotation,levels=order)
GAT2 <- GAT2[match(order, GAT2$annotation),]
row.names(GAT2) <- c(1:nrow(GAT2))
GAT2$annotation <- factor(GAT2$annotation, levels = rev(order))

# Colors for Custom Order
n = nrow(GAT2)
cols = gg_color_hue(n)
names(cols) <- order

# Melt
GAT2 <- melt(GAT2)

# Signifcance
GAT2.1 <- dcast(GAT2.1, annotation ~ track, value.var = "qvalue")
order <- c("Enhancers", "1-5 kb Upstream", "Promoter", "5' UTR", "Exon", "Exon/Intron Boundaries", "Intron", "3' UTR", "Intergenic")
GAT2.1$annotation <- factor(GAT2.1$annotation,levels=order)
GAT2.1 <- GAT2.1[match(order, GAT2.1$annotation),]
row.names(GAT2.1) <- c(1:nrow(GAT2.1))
GAT2.1$annotation <- factor(GAT2.1$annotation, levels = rev(order))
GAT2.1 <- melt(GAT2.1)
GAT2.1$signif <- ifelse(GAT2.1$value < 0.05,1,0)
GAT2.1 <- GAT2.1[-3]
GAT2 <- merge(GAT2, GAT2.1, by =c("annotation", "variable"))

# Plot
gene <- ggplot(data=GAT2, aes(x=annotation, y =value, fill = annotation)) +
  geom_bar(stat="identity") + coord_flip() + labs(y = "Fold Enrichment", x= element_blank()) + theme_classic() + 
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16),
        strip.text = element_text(size = 16), legend.text = element_text(size = 14),
        legend.title=element_text(size = 14)) +
  scale_y_continuous(limits = c(-3.5,3.5), breaks = c(-3,-2,-1,0,1,2,3)) + 
  scale_fill_manual(values = cols, breaks = order, name = "Annotation") + facet_grid(~variable) +
  geom_hline(yintercept=0) +
  geom_text(data = GAT2[(GAT2$signif ==1 & GAT2$value >0), ], label = "*", size=8, show.legend = FALSE, nudge_y = 0.5, nudge_x = -0.09) +
  geom_text(data = GAT2[(GAT2$signif ==1 & GAT2$value <0), ], label = "*", size=8, show.legend = FALSE, nudge_y = -0.5, nudge_x = -0.09)

ggsave("Gat_genic_hyper_hypo.pdf",
       plot = gene,
       device = NULL,
       width = 11,
       height = 6)

