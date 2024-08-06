rm(list = ls())
library(ggbeeswarm)
library(ggrepel)
setwd('E:/experiment/STR2/02.STR_features')

theme_cust <- theme_bw() + 
  theme(plot.title = ggplot2::element_text(size=10,  color = "black"),
        legend.text = ggplot2::element_text(size=10,  color = "black"),
        legend.title =  ggplot2::element_text(size=10,  color = "black"),
        axis.title =  ggplot2::element_text(size=10,  color = "black"),
        axis.text =  ggplot2::element_text(size=10,  color = "black"),
        strip.text = ggplot2::element_text(size=10, vjust = 1,  color = "black"),
        strip.background = ggplot2::element_blank(), 
        panel.grid = ggplot2::element_blank(),
        text = ggplot2::element_text(family="Helvetica"))

period_size_color <- alpha(c("1" = "#5B9BD5", "2" = "#ED7D31", "3" = "black",
                             "4" = "#FFC000", "5" = "#4472C4", "6" = "#70AD47"), alpha = 1)
fisher_enrich_func <- function(pt_deg){
  sis <- pt_deg$sig_n  
  silp <- pt_deg$n - sis 
  fis <- pt_deg$sig_total - sis   
  filp <- pt_deg$total - pt_deg$sig_total - silp  
  ftp <- fisher.test(matrix(c(sis,silp,fis,filp), 2, 2), alternative='greater')
  return(ftp$p.value)
}

table_s1_poly<- read.csv('result.txt',sep = '\t')

enrichMotif_polySTR_region_stats <- table_s1_poly %>%
  dplyr::group_by(gfeature, motif_geno) %>%
  dplyr::count(name = "sig_n") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(gfeature) %>%
  dplyr::mutate(sig_total = sum(sig_n)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(motif_geno) %>%
  dplyr::mutate(n = sum(sig_n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total = sum(sig_n)) %>%
  dplyr::group_by(gfeature, motif_geno) %>%
  dplyr::do(data.frame(x = fisher_enrich_func(.))) %>%
  dplyr::rename(fisherp = x) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(fisherp_adj = p.adjust(fisherp, method = "bonferroni")) %>%
  dplyr::mutate(group_factor = "gfeature, motif_geno",
                group_factor_catogory = paste0(gfeature, ", ", motif_geno),
                method = "one-sided Fisher's Exact test",
                padjustment = "BF")


enrichMotif_polySTR_region <- enrichMotif_polySTR_region_stats %>% 
  dplyr::filter(fisherp_adj<0.05)   %>%  
  dplyr::mutate(logp=-log10(fisherp_adj)) %>% 
  dplyr::mutate(motif_length=nchar(motif_geno))%>% 
  dplyr::arrange(  motif_length )  %>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::top_n(3, logp) 


enrichMotif_polySTR_region$gfeatures<- factor(enrichMotif_polySTR_region$gfeature,levels = c("IGR","upstream","five_prime_UTR","CDS","intron","three_prime_UTR","downstream"))



fig_2d <- ggplot(enrichMotif_polySTR_region,
                 aes(y=gfeatures,x=logp,color=factor(motif_length)))+
  ggbeeswarm::geom_beeswarm(cex=2,groupOnX = F,size=3) +
  theme_cust+
  theme( legend.position = "none") +
  scale_color_manual(values=period_size_color) +
  ylab("Genomic\nfeatures")+
  xlab(expression(-log[10](italic(p)))) 


fig_2d
ggsave('outputwithtitle.png',width = 7,height = 6)
