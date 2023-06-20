#require(tidyverse)
#require(scales)

#rna <- read_tsv("/home/rstudio/staRfish/data_test/data_mrna_seq_v2_rsem.txt")
#rpot <- read_tsv("/home/rstudio/staRfish/data_test/data_protein_quantification.txt") %>%
#  separate(Composite.Element.REF,into=c("gene","gene_2")) %>% select(!gene_2) %>% na.omit() %>%
#  mutate(avg_expr=rowMeans(pick(!contains("gene")),to=c((-1),1)))
#rna_2<- rna[3:ncol(rna)][colnames(rna[3:ncol(rna)]) %in% colnames(rpot[3:ncol(rpot)]) ]
#rna_3 <- rna_2 %>% mutate(avg_expr = rowMeans(log2(pick(everything())+1)),gene=rna$Hugo_Symbol) %>%
#  .[!is.na(.$gene),] %>% select(gene,avg_expr)
#rpot_2 <- rpot %>% select(gene,avg_expr)

#total <- inner_join(rpot_2,rna_3, by="gene", suffix=c(".prot",".rna"))

#ggplot(total,aes(x=avg_expr.prot,y=avg_expr.rna)) + geom_point()


#rna_4 <- rna_2 %>% mutate(gene=rna$Hugo_Symbol)
#total_all_samp <- inner_join(rpot,rna_4,by="gene",suffix=c(".prot",".rna")) %>% select(-gene)

#protein <- read_tsv("/home/rstudio/staRfish/data_test/data_protein_quantification.txt") %>%
#  separate(Composite.Element.REF,into=c("gene","gene_2")) %>% select(!gene_2) %>% na.omit()
#rna <- read_tsv("/home/rstudio/staRfish/data_test/data_mrna_seq_fpkm.txt") %>% dplyr::rename("gene"=Hugo_Symbol)
require(tidyverse)
### DATA PROCESSING ###
gather_rna_prot_data <- function(rna,protein) {
  gathered_data <- dplyr::inner_join(x=protein,y=rna,by="gene",suffix=c(".prot",".rna"))
  return(gathered_data)
}

create_rna_prot_correlation <- function(rna,protein) {
  r <- rna %>% na.omit()
  r$gene <- make.unique(r$gene)
  p <- protein %>% na.omit()
  p$gene <- make.unique(p$gene)

  common_genes <- dplyr::intersect(unique(r$gene),unique(p$gene)) %>% sort()
  rna_tib <- r %>% dplyr::filter(gene %in% common_genes) %>% arrange(gene)
  protein_tib <- p %>% dplyr::filter(gene %in% common_genes) %>% arrange(gene)

  trans_data <- as.matrix(rna_tib[,-1])
  prot_data <- as.matrix(protein_tib[,-1])

  gene_correlations <- diag(cor(t(trans_data), t(prot_data)))

  correlations_df <- data.frame(gene = common_genes, correlation = gene_correlations) %>% .[order(.$correlation),]
  return(correlations_df)
}

create_pw_df <- function(corr_df) {
  require(clusterProfiler)
  require(org.Hs.eg.db)

  genes <- corr_df$gene

  gene_ids <- mapIds(org.Hs.eg.db, keys=genes, column="ENTREZID", keytype="SYMBOL", multiVals="first")

  kegg_annotations <- enrichKEGG(gene         = gene_ids,
                                 organism     = 'hsa',
                                 pAdjustMethod = "none",
                                 minGSSize    = 1,
                                 maxGSSize    = 500)

  kegg_df <- as.data.frame(kegg_annotations) %>% select(ID,Description,geneID)

  expanded_kegg_df <- kegg_df %>%
    mutate(geneID = strsplit(as.character(geneID), "/")) %>%
    unnest(geneID) %>%
    group_by(ID) %>%
    mutate(gene_index = row_number()) %>%
    pivot_wider(names_from = gene_index, values_from = geneID, names_prefix = "gene_") %>%
    select(-ID) %>%
    pivot_longer(!Description) %>% ungroup() %>% select(-name) %>% na.omit() %>%
    mutate(value=getSYMBOL(as.character(value),data='org.Hs.eg')) %>%
    rename(gene=value)

  combined_pw_data <- corr_df %>%
    inner_join(expanded_kegg_df, by = "gene")
  return(combined_pw_data)

}

### START PLOTTING DATA ###
cor_matrix_samples <- function(gathered_data) {
  require(reshape2)
  require(ggplot2)
  require(viridis)
  gathered <- gathered_data %>% select(-gene) %>% na.omit()
  cormat <- round(cor(gathered),2)
  mcormat <- melt(cormat)

  p <- ggplot(data = mcormat, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    theme(axis.text.x=element_blank(), #remove x axis labels
          axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y=element_blank(),  #remove y axis labels
          axis.ticks.y=element_blank(),  #remove y axis ticks
          axis.title.x=element_blank(),
          axis.title.y=element_blank()
    ) + scale_fill_viridis(option="inferno")
  return(p)
}


correlation_plot_levels <- function(correlations_df) {
  correlations_df$gene <- factor(correlations_df$gene , levels = correlations_df$gene)
  Labels_top = rbind(head(arrange(correlations_df,desc(correlation)),10), head(arrange(correlations_df,correlation),10))
  correlations_df$label = if_else(correlations_df$gene %in% Labels_top$gene,
                                  correlations_df$gene, NA)
  correlations_df$color <- if_else(correlations_df$correlation<0,(correlations_df$correlation*2),correlations_df$correlation)

  p <- ggplot(correlations_df, aes(x = gene, y = correlation,color=color)) +
    geom_point(stat = "identity") +
    theme(axis.title.x= element_blank(),
          axis.text.x = element_blank(),
          panel.background = element_rect(fill="white"),
          legend.position = "none",
          plot.title= element_text(hjust = 0.5)) +
    geom_hline(yintercept=0)+
    ggrepel::geom_text_repel(aes(label = label),color="black",
                             size = 4, show.legend = FALSE, max.overlaps=10000,box.padding=0.75) +
    expand_limits(x= c(-50, length(levels(correlations_df$gene))*1.01 ))+
    colorspace::scale_color_continuous_divergingx(palette="PiYG") +
    labs(x = "Gene", y = "Correlation between Transcriptomics and Proteomics", title = "Gene Correlations")
  suppressWarnings(print(p))
  df <- rbind(correlations_df[1:10,],correlations_df[(length(correlations_df[[1]])-10):length(correlations_df[[1]]),])[1:2]
  return(df)
}

KEGG_correlation_plot <- function(combined_pw_data) {
  average_correlation_per_pathway <- combined_pw_data %>%
    group_by(Description) %>%
    summarise(average_correlation = mean(correlation, na.rm = TRUE))

  Labels_top = rbind(head(arrange(average_correlation_per_pathway,desc(average_correlation)),10), head(arrange(average_correlation_per_pathway,average_correlation),10))
  average_correlation_per_pathway$label = if_else(average_correlation_per_pathway$Description %in% Labels_top$Description,
                     average_correlation_per_pathway$Description, NA)
  average_correlation_per_pathway$color = if_else(average_correlation_per_pathway$Description %in% Labels_top$Description,
                                                  1, NA)

  p <- ggplot(average_correlation_per_pathway, aes(x=reorder(Description, average_correlation), y=average_correlation,color=color)) +
    geom_point(stat="identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggrepel::geom_text_repel(aes(label = label), color="black",
                             size = 4, show.legend = FALSE, max.overlaps=10000) +
    theme(axis.title.x= element_blank(),
          axis.text.x = element_blank(),
          panel.background = element_rect(fill="white"),
          legend.position = "none",
          plot.title= element_text(hjust = 0.5)) +
    labs(x = "KEGG Pathway", y = "Average Correlation", title = "Average Correlation of Protein/Transcript Levels in KEGG Pathways")
  suppressWarnings(print(p))
}

plot_gene <- function(gathered_data,gene_name) {
  require(ggpubr)
  stopifnot("Gene not available in dataset" = gene_name %in% gathered_data$gene)
  gene_tib <- gathered_data %>% filter(gene==gene_name) %>%
    pivot_longer(!gene,names_to=c("patient","type"),names_sep="\\.",values_to="expr") %>%
    pivot_wider(id_cols="patient",names_from="type",values_from="expr") %>% na.omit()

  ggplot(gene_tib,aes(x=prot,y=rna)) +
    geom_point() +
    geom_smooth(method='lm') +
    stat_regline_equation(label.y = max(gene_tib$rna), aes(label = ..rr.label..)) +
    ggtitle(gene_name) +
    theme(panel.background = element_rect(fill="white"),
          legend.position = "none",
          plot.title= element_text(hjust = 0.5)) +
    xlab("Protein expression") +
    ylab("RNA expression")

}
