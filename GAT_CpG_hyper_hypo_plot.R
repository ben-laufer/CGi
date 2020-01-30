# GAT CpG Annotation Enrichment Plot
# Ben Laufer

# Takes the output of GAT for hg38 CpG annotation enrichments and makes a ggplot2

if (!requireNamespace(c("ggplot2", "reshape2"), quietly = TRUE))
  install.packages(c("ggplot2", "reshape2"))
library(ggplot2)
library(reshape2)

# Load
GAT <- read.delim("GAT_CpG_hyper_hypo.tsv")

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
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
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

ggsave("GAT_CpG_hyper_hypo.pdf",
       plot = CpG,
       device = NULL,
       width = 11,
       height = 6)
