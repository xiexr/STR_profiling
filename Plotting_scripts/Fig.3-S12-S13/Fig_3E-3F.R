rm(list = ls())
setwd('E:/eSTR_enrichment')
library(magrittr)
library(dplyr)
library(readr)
library(ggplot2)
library(cowplot)
eSTR_feature = read_tsv("eSTRgat.out") %>%
  mutate(annotation = case_when(
    annotation == "CDS" ~ "CDS",
    annotation == "intron" ~ "Intron",
    annotation == "intergenic" ~ "Intergenic",
    annotation == "three_prime_UTR" ~ "3'UTR",
    annotation == "five_prime_UTR" ~ "5'UTR",
    annotation == "upstream" ~ "Promoter",
    annotation == "downstream" ~ "Downstream"
  ))

eSTR_feature = mutate(eSTR_feature,fdr_bh = p.adjust(pvalue, method = "BH"))

eNip = read_tsv('Allleaves.out')%>%
  mutate(fdr_bh = p.adjust(pvalue, method = "BH"))
plot_enrich <- function(dat=NULL) {
  dat %>%
    mutate(sig = case_when(
      fdr_bh >= 0.05 ~ "n.s.",
      fold > 1 ~ "pos",
      TRUE ~ "neg"
    )) %>%
    mutate(sig = factor(sig, levels = c("pos", "neg", "n.s."))) %>%
    ggplot(aes(x = fold, y = reorder(annotation, fold), color = sig)) +
    geom_point() +
    geom_pointrange(aes(xmin = observed/CI95low, xmax = observed/CI95high),size=1,lwd=1) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("#ca3226", "#476fa9", "gray")) +
    #scale_color_manual(values = c("pos" = "#ca3226", "neg" = "#476fa9", "n.s." = "gray"))+
    labs(x = "Fold enrichment/depletion", y = NULL,
         color = "pvalue") +
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())+
    theme(axis.line.x=element_line(linetype=1,color="black",linewidth=.5))+
    theme(axis.line.y=element_line(linetype=1,color="black",linewidth=.5))+
    theme(axis.ticks.x=element_line(color="black",size=.5,lineend = 1))+
    theme(axis.ticks.y=element_line(color="black",size=.5,lineend = 10))+
    theme(legend.position = "none")}
 

p1 <- plot_enrich(eSTR_feature)
ggsave(p1, filename = "eSTR_feature.pdf", height = 4, width = 4)

p2 <- plot_enrich(eNip)
ggsave(p2, filename = "eNip_leaves.pdf", height = 4, width = 4)

p.nip <- plot_grid(p1, p2, labels = c("", ""), nrow = 1, rel_heights = c(110, 110),hjust = -3)
ggsave(p.nip, filename = "combined_plot.pdf",dpi=300, height = 4, width = 7)

