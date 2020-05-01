#!/usr/bin/env Rscript

# GAT CpG and genic annotation enrichment plots
# Takes the output of GAT for hg38 CpG and genic annotation enrichments and makes a ggplot2
# Ben Laufer

# Packages ----------------------------------------------------------------

if (!requireNamespace(c("tidyverse", "wesanderson"), quietly = TRUE))
  install.packages("tidyverse", repos = "https://cloud.r-project.org")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
options(readr.num_columns = 0)

# CpG ---------------------------------------------------------------------

cat("\n","Tidying CpG annotations...")
GAT <- readr::read_tsv("Gat_CpG_hyper_hypo_results.tsv") %>%
  dplyr::mutate(annotation = dplyr::recode_factor(annotation,
                                                  "hg38_cpg_inter" = "Open Sea",
                                                  "hg38_cpg_shelves" = "CpG Shelf",
                                                  "hg38_cpg_shores" = "CpG Shore",
                                                  "hg38_cpg_islands" = "CpG Island")
  ) %>%
  dplyr::mutate(fold = dplyr::case_when(fold < 1 ~ -1/fold,
                                        fold >= 1 ~ fold)
  ) %>%
  dplyr::mutate(signif = dplyr::case_when(qvalue <= 0.05 ~ 1,
                                          qvalue > 0.05 ~ 0)
  ) %>%
  dplyr::select(annotation,
                track,
                fold,
                signif)
cat("Done.")

cat("\n", "Plotting...")
(ggplot(data = GAT,
        aes(annotation,
            y = fold,
            fill = annotation)
        ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(y = "Fold Enrichment",
         x = element_blank()) +
    theme_classic() + 
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)
    ) +
    scale_y_continuous(limits = c(-4,4),
                       breaks = c(-4,-3,-2,-1,0,1,2,3,4)
    ) +
    scale_fill_manual(values = c("forestgreen", "goldenrod2", "dodgerblue", "blue3"),
                      breaks = GAT$annotation %>%
                        levels() %>%
                        rev(),
                      name = "Annotation") +
    facet_grid(~track) +
    geom_hline(yintercept = 0) +
    geom_text(data = GAT[(GAT$signif == 1 & GAT$fold > 0), ],
              label = "*",
              size = 8,
              show.legend = FALSE,
              nudge_y = 0.5,
              nudge_x = -0.05) +
    geom_text(data = GAT[(GAT$signif == 1 & GAT$fold < 0), ],
              label = "*",
              size = 8,
              show.legend = FALSE,
              nudge_y = -0.5,
              nudge_x = -0.05)
) %>%  
  ggsave("Gat_CpG_hyper_hypo.pdf",
         plot = .,
         width = 11,
         height = 6)
cat("Done.")

# Genic  ------------------------------------------------------------------

cat("\n","\n", "Tidying genic annotations...")
GAT2 <- readr::read_tsv("Gat_genic_hyper_hypo_results.tsv") %>% 
  dplyr::mutate(annotation = dplyr::recode_factor(annotation,
                                                  "intergenic_merged.bed" = "Intergenic",
                                                  "threeUTRs_merged.bed" = "3' UTR",
                                                  "introns_merged.bed" = "Intron",
                                                  "boundaries_merged.bed" = "Exon/Intron Boundaries",
                                                  "exons_merged.bed" = "Exon",
                                                  "fiveUTRs_merged.bed" = "5' UTR",
                                                  "promoters_merged.bed" = "Promoter",
                                                  "onetofivekb_merged.bed" = "1-5 kb Upstream",
                                                  "enhancers_merged.bed" = "Enhancers")
  ) %>%
  dplyr::mutate(fold = dplyr::case_when(fold < 1 ~ -1/fold,
                                        fold >= 1 ~ fold)
  ) %>%
  dplyr::mutate(signif = dplyr::case_when(qvalue <= 0.05 ~ 1,
                                          qvalue > 0.05 ~ 0)
  ) %>%
  dplyr::select(annotation,
                track,
                fold,
                signif)
cat("Done.")

cat("\n", "Plotting...")
(ggplot(data = GAT2,
        aes(x = annotation,
            y = fold,
            fill = annotation)
        ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(y = "Fold Enrichment",
         x = element_blank()
    ) +
    theme_classic() + 
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)
    ) +
    scale_y_continuous(limits = c(-3.5,3.5),
                       breaks = c(-3,-2,-1,0,1,2,3)
    ) + 
    scale_fill_manual(values =  GAT2$annotation %>%
                        nlevels() %>% 
                        wesanderson::wes_palette("Zissou1", n = ., type = "continuous") %>%
                        rev(),
                      breaks = GAT2$annotation %>%
                        levels() %>%
                        rev(),
                      name = "Annotation") +
    facet_grid(~track) +
    geom_hline(yintercept = 0) +
    geom_text(data = GAT2[(GAT2$signif == 1 & GAT2$fold > 0), ],
              label = "*",
              size = 8,
              show.legend = FALSE,
              nudge_y = 0.5,
              nudge_x = -0.09) +
    geom_text(data = GAT2[(GAT2$signif == 1 & GAT2$fold < 0), ],
              label = "*",
              size = 8,
              show.legend = FALSE,
              nudge_y = -0.5,
              nudge_x = -0.09)
) %>% 
  ggsave("Gat_genic_hyper_hypo.pdf",
         plot = ., 
         width = 11,
         height = 6)
cat("Done.", "\n", "\n")
