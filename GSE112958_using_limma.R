library(limma)
library(GEOquery)
library(ggplot2)

gse <- getGEO("GSE112958")
exprs <- exprs(gse[[1]])

background_corrected_data <- backgroundCorrect(exprs, method="normexp")

norm_data <- normalizeBetweenArrays(background_corrected_data, method="quantile")

design <- model.matrix(~0+factor(c(rep("CTRL", 89), rep("EF", 89))))
colnames(design) <- levels(factor(c(rep("CTRL", 89), rep("EF", 89))))

fit <- lmFit(norm_data, design)
fit <- eBayes(fit)

all_genes <- topTable(fit, coef=1, n=Inf, adjust.method="holm", sort.by="P")

pvalues <- all_genes$adj.P.Val;

logFC <- all_genes$logFC

log_pvalues <- (-1)*(log10(pvalues))

# calculated the max and min values for choosing the thresholds

max(log_pvalues)  #309.6414
min(log_pvalues)  #4.168953

max(logFC)        #9.802698
min(logFC)        #0.06316902

all_genes <- cbind(all_genes,log_pvalues)

df_volcano_limma <- data.frame(log_pvalues,logFC)

# making a volcano plot to visualize the diff exp genes
df_volcano_limma$color <- ifelse(df_volcano_limma$log_pvalues > 201.301 & df_volcano_limma$logFC > 6, "red",
                           ifelse(df_volcano_limma$log_pvalues > 201.301 & df_volcano_limma$logFC < 3, "blue", "black"))


# Create the volcano plot
ggplot(df_volcano_limma, aes(x = logFC, y = log_pvalues, color = color)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("black", "red", "blue")) +
  theme_classic() +
  xlab("log2 fold change") +
  ylab("-log10 p-value") +
  ggtitle("Volcano plot") +
  geom_hline(yintercept = 201.301, linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = 3, linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = 6, linetype = "dashed", color = "gray30")

# Keep only the differential expressed genes
diff_exp_genes <- all_genes[all_genes$log_pvalues > 201.301 & (all_genes$logFC < 3 | all_genes$logFC >6),]


# storing the differential and all gene ids in vectors
diff_exp_gene_ids <- rownames(diff_exp_genes)
all_gene_ids <- rownames(all_genes)

BiocManager::install("illuminaHumanv4.db")

library(illuminaHumanv4.db)

# Convert ILMN probe IDs to Entrez gene IDs for the differential genes
ilmnToEntrez_limma <- select(illuminaHumanv4.db, keys=diff_exp_gene_ids, columns="ENTREZID", keytype="PROBEID")
entrezGenes_limma <- ilmnToEntrez_limma$ENTREZID
entrezGenes_limma

# Convert ILMN probe IDs for the gene universe , i.e. all 47231 genes
ilmnToEntrez_universe_limma <- select(illuminaHumanv4.db, keys=all_gene_ids, columns="ENTREZID", keytype="PROBEID")
entrezGenes_all_limma <- ilmnToEntrez_universe_limma$ENTREZID

library(org.Hs.eg.db)
library(clusterProfiler)

# omit the Null values from the vectors
entrezGenes_limma <- na.omit(entrezGenes_limma)
entrezGenes_all_limma <- na.omit(entrezGenes_all_limma)


# Performing enrichment analysis using enrichGO function
goEnrichment <- enrichGO(gene = entrezGenes_limma, 
                         universe = unique(entrezGenes_all_limma), 
                         OrgDb = org.Hs.eg.db, 
                         keyType = "ENTREZID", 
                         ont = c("BP","MF","CC","ALL"),
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.1,
                         readable = TRUE)

# converting the results in a data frame to analyse
enrichedTerms <- as.data.frame(goEnrichment)

# converting the dataframe to a csv file
write.csv(enrichedTerms,"enriched_genes.csv")


#########PLOTS########

# Subset the enriched genes dataframe to include only significant GO terms (p-value or q-value < 0.05)
enrichedTermsSig <- enrichedTerms[enrichedTerms$pvalue < 0.05, ]

# Create a bar plot of enriched GO terms
ggplot(enrichedTermsSig, aes(x = Description, y = -log10(pvalue), fill = p.adjust)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title = "Enriched GO terms", x = "GO term", y = "-log10(p-value)")


# Create a dot plot of enriched GO terms
ggplot(enrichedTermsSig, aes(x = -log10(pvalue), y = Description, color = p.adjust)) + 
  geom_point(size = 2) + 
  theme(axis.text.y = element_text(size = 8)) + 
  labs(title = "Enriched GO terms", x = "-log10(p-value)", y = "GO term")







