line_compare_generate <- function (transcriptomics,proteomics) {
  stopifnot(all(transcriptomics$patient_id %in% proteomics$patient_id))

}

library(tidyverse)
library(scales)

rna <- read_tsv("/home/rstudio/staRfish/data_test/data_mrna_seq_v2_rsem.txt")
rpot <- read_tsv("/home/rstudio/staRfish/data_test/data_protein_quantification.txt") %>%
  separate(Composite.Element.REF,into=c("gene","gene_2")) %>% select(!gene_2) %>% na.omit() %>%
  mutate(avg_expr=rowMeans(pick(!contains("gene")),to=c((-1),1)))
rna_2<- rna[3:ncol(rna)][colnames(rna[3:ncol(rna)]) %in% colnames(rpot[3:ncol(rpot)]) ]
rna_3 <- rna_2 %>% mutate(avg_expr = rowMeans(log2(pick(everything())+1)),gene=rna$Hugo_Symbol) %>%
  .[!is.na(.$gene),] %>% select(gene,avg_expr)
rpot_2 <- rpot %>% select(gene,avg_expr)

total <- inner_join(rpot_2,rna_3, by="gene", suffix=c(".prot",".rna"))

ggplot(total,aes(x=avg_expr.prot,y=avg_expr.rna)) + geom_point()


rna_4 <- rna_2 %>% mutate(gene=rna$Hugo_Symbol)
total_all_samp <- inner_join(rpot,rna_4,by="gene",suffix=c(".prot",".rna")) %>% select(-gene)

protein <- read_tsv("/home/rstudio/staRfish/data_test/data_protein_quantification.txt") %>%
  separate(Composite.Element.REF,into=c("gene","gene_2")) %>% select(!gene_2) %>% na.omit()
rna <- read_tsv("/home/rstudio/staRfish/data_test/data_mrna_seq_fpkm.txt") %>% dplyr::rename("gene"=Hugo_Symbol)


cor_matrix <- function(rna,protein) {
  require(reshape2)
  require(ggplot2)
  require(viridis)
  gathered <- dplyr::inner_join(x=protein,y=rna,by="gene",suffix=c(".prot",".rna")) %>% select(-gene) %>% na.omit()
  cormat <- round(cor(gathered),2)
  mcormat <- melt(cormat)

  ggplot(data = mcormat, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    theme(axis.text.x=element_blank(), #remove x axis labels
          axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y=element_blank(),  #remove y axis labels
          axis.ticks.y=element_blank(),  #remove y axis ticks
          axis.title.x=element_blank(),
          axis.title.y=element_blank()
    ) + scale_fill_viridis(option="inferno")

}
