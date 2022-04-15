#???Over-Representation Analysis with ClusterProfiler
#  https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/

# Install and load packages
  # BiocManager::install("clusterProfiler")
  # BiocManager::install("pathview")
  # install.packages("wordcloud")
  library(clusterProfiler)
  
# Annotations
  organism = "org.Dm.eg.db"
  #BiocManager::install(organism, character.only = TRUE)
  library(organism, character.only = TRUE)

##### Prepare Input #####
  # reading in input from deseq2
  df = read.csv("drosphila_example_de.txt", header=TRUE)
  
  # we want the log2 fold change 
  original_gene_list <- df$log2FoldChange
  
  # name the vector
  names(original_gene_list) <- df$X
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  # Exctract significant results (padj < 0.05)
  sig_genes_df = subset(df, padj < 0.05)
  
  # From significant results, we want to filter on log2fold change
  genes <- sig_genes_df$log2FoldChange
  
  # Name the vector
  names(genes) <- sig_genes_df$X
  
  # omit NA values
  genes <- na.omit(genes)
  
  # filter on min log2fold change (log2FoldChange > 2)
  genes <- names(genes)[abs(genes) > 2]

##### Create enrichGO object #####  
  # Create the object
  go_enrich <- enrichGO(gene = genes,
                        universe = names(gene_list),
                        OrgDb = organism, 
                        keyType = 'ENSEMBL',
                        readable = T,
                        ont = "BP",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)
  
##### Output #####
  ## Upset Plot
  # BiocManager::install("enrichplot")
  # BiocManager::install("ggupset")
  library(enrichplot)
  library(ggupset)
  upsetplot(go_enrich)

  ## Barplot
  barplot(go_enrich, 
          drop = TRUE, 
          showCategory = 10, 
          title = "GO Biological Pathways",
          font.size = 8)
  
  ## Dotplot
  dotplot(go_enrich)
  
  ## Encrichment map:
  emapplot(go_enrich)
  
  ## Enriched GO induced graph:
  goplot(go_enrich, showCategory = 10)
  
  ## Category Netplot
  # categorySize can be either 'pvalue' or 'geneNum'
  cnetplot(go_enrich, categorySize="pvalue", foldChange=gene_list)
  

# ##### Error part #####
#   
#   ##### Wordcloud #####
#   # install.versions(c('RColorBrewer'), c('1.1.2'))
#   library(wordcloud)
#   
#   ## Wordcloud
#   wcdf<-read.table(text=go_enrich$GeneRatio, sep = "/")[1]
#   wcdf$term<-go_enrich[,2]
#   wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)
#   
#       