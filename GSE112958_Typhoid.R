library(GEOquery)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(limma)

# host gene signatures for distinguishing enteric fever(caused by salmonella typhi) 
# from general febrile diseases 
# diagnosis of enteric fever is not accurate in most cases

gse_id = 'GSE112958'
data <- getGEO(gse_id,GSEMatrix = TRUE)

#sample data : the characteristics_ch1 col contains the group name
#EF : enteric fever ,CTRL : control group, sEF : suspected enteric fever 

# data attributes
cat("Accession ID sample 1 : ",data$GSE112958_series_matrix.txt.gz$geo_accession[1])
cat("Platform : ",data$GSE112958_series_matrix.txt.gz$platform_id[1])
cat("type : ", data$GSE112958_series_matrix.txt.gz$type[1])
cat("source name : ", data$GSE112958_series_matrix.txt.gz$source_name_ch1[1])
cat(" organism name : ",data$GSE112958_series_matrix.txt.gz$organism_ch1[1])
cat(" growth protocol : ",data$GSE112958_series_matrix.txt.gz$growth_protocol_ch1[1])
cat("data processing : ",data$GSE112958_series_matrix.txt.gz$data_processing[1])
cat(" transcripts : ",data$GSE112958_series_matrix.txt.gz$data_row_count[1]) 

pdata <- pData(data[[1]])
pdata
#feature data : gene 
fdata <- fData(data[[1]])

# expression value of the data
data_matrix <- exprs(data[[1]])

#dimension of the matrix
dim(data_matrix)

max(data_matrix)  #15.8674
min(data_matrix)  #4.2624


#converted the matrix to a data frame
data_dataframe <- as.data.frame(data_matrix)

#making a box plot for x = samples/columns and y = intensity values(un-normalized data)

data_long <- pivot_longer(data_dataframe, everything(), names_to = "variable")
ggplot(data_long,aes(x = variable , y = value )) + geom_boxplot(fill = "#0099f8")

#performing background correction
data.bc <- backgroundCorrect(data_dataframe, method="normexp", offset=10)

#performing quantile normalization
data.norm <- normalizeQuantiles(data.bc)

#re-writing the corrected dataframe
data_dataframe <- as.data.frame(data.norm)

# making a box plot again
boxplot(data_dataframe, main="Normalized Intensity Distribution")

# as all the values 
log2(data_dataframe)
max(data_dataframe)   # 19.974
min(data_dataframe)   # 10.0054

# control samples in data and diseased are initialized to zero
control_samples <- 0
diseased_samples <- 0

# for all of the 178 samples we are checking if the sample is CTRL 
# or it is EF and then increasing the count
for(i in 1:178){
  if(identical(pdata$characteristics_ch1[i],"group: CTRL")){
      control_samples <- control_samples + 1
  }
  else if(identical(pdata$characteristics_ch1[i],"group: EF")){
      diseased_samples <- diseased_samples + 1
  }
}
# printing the samples count 
control_samples
diseased_samples

# Now we need to make a new dataframe which contains first 88 samples as normal 
# and next 50 samples as diseased and then perform t-test on each gene of new 
# data frame and calculate pvalues in a vector

#pvalues of all the genes where the index specify the gene number


final_dataframe <- data_dataframe[1]


for(i in 2:178){
  if(identical(pdata$characteristics_ch1[i],"group: CTRL")){
      final_dataframe[pdata$geo_accession[i]] <- data_dataframe[[i]]
  }
}

for(i in 2:178){
  if(identical(pdata$characteristics_ch1[i],"group: EF")){
      final_dataframe[pdata$geo_accession[i]] <- data_dataframe[[i]]
  }
}

#now the final dataframe must have first 88 samples (1-88) as CTRL and 
#next 50 samples as EF (89 - 138)

#now we can subset the final_dataframe and perform t-tests

#fold change array is also created

# These arrays will contain the fold chain values and pvalues 
# for each gene in order i.e 1 to 47231 genes that are in data frame
fold_changes <- c()
p_values <- c()


for(i in 1:47231){
   gene_i_normal <- final_dataframe[i,1:88]
   gene_i_disease <- final_dataframe[i,89:138]
   
   t_test <- t.test(gene_i_normal,gene_i_disease)
   
   mean1 <- t_test$estimate[[1]]
   mean2 <- t_test$estimate[[2]]
   
   f_c <- mean2/mean1
   fold_changes <- c(fold_changes,f_c)
   p_v <- t_test$p.value
   p_values <- c(p_values,p_v)
}


max(p_values)
min(p_values)
# now p_values vector must have all the p_values 
# now we can perform the holm correction

holm_adj <- p.adjust(p_values, method = "holm")
holm_adj

max(holm_adj)
min(holm_adj)


# now we can do the -log(pvalue) transform on our pvalues
#applied log2 negative transform to the corrected pvalues
log_holm_adj <- (-1)*(log10(holm_adj))
max(log_holm_adj)   # 43.08
min(log_holm_adj)   # 0

# and can form a log(FC) vector 
log_fold_change <- log2(fold_changes)
max(log_fold_change)    #0.54  
min(log_fold_change)    #-0.18

# then finally plot a volcano plot out of it
df_volcano <- data.frame(log_fold_change,log_holm_adj)

# FDR set to 0.1% false positives 

df_volcano$color <- ifelse(df_volcano$log_holm_adj > 4 & df_volcano$log_fold_change > 0.05, "red",
                   ifelse(df_volcano$log_holm_adj > 4 & df_volcano$log_fold_change < -0.05, "blue", "black"))

print(table(df_volcano$color))

# Create the volcano plot
ggplot(df_volcano, aes(x = log_fold_change, y = log_holm_adj, color = color)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("black", "blue", "red")) +
  theme_classic() +
  xlab("log2 fold change") +
  ylab("-log10 p-value") +
  ggtitle("Volcano plot") +
  geom_hline(yintercept = 4, linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = -0.05, linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "gray30")


transcript_ids <- rownames(final_dataframe)
# this matrix stores the transcript IDs from the final_dataframe in order
# and the respective -log10(pvalues) and log2(FC) values 

transcript_pv_fc <- cbind(transcript_ids,log_holm_adj,log_fold_change)

transcript_pv_fc_df <- as.data.frame(transcript_pv_fc)

# This vector will store the gene IDs which are differential
differential_genes_id <- c()


# selecting the transcript IDs with certain parameters
for(i in 1:nrow(transcript_pv_fc_df)){
    if(transcript_pv_fc_df[i,"log_holm_adj"] > 4 & (transcript_pv_fc_df[i,"log_fold_change"] > 0.05 | transcript_pv_fc_df[i,"log_fold_change"] < -0.05)){
          differential_genes_id <- c(differential_genes_id,transcript_pv_fc_df[i,"transcript_ids"])
    }
}

length(differential_genes_id)  #953?

differential_genes_id

# differential expressed gene IDs are stored in a text file

write.table(differential_genes_id, file = "diff_genes.txt", row.names = FALSE, col.names = FALSE)

all_genes <- fdata$ID

write.table(all_genes, file = "all_genes.txt", row.names = FALSE, col.names = FALSE)



########################################################################################
                                   #enrichment analysis

BiocManager::install("illuminaHumanv4.db")

library(illuminaHumanv4.db)

# Vector of ILMN probe IDs
differential_genes_id
all_genes

# Convert ILMN probe IDs to Entrez gene IDs for the differential genes
ilmnToEntrez <- select(illuminaHumanv4.db, keys=differential_genes_id, columns="ENTREZID", keytype="PROBEID")
entrezGenes <- ilmnToEntrez$ENTREZID
entrezGenes

# Convert ILMN probe IDs for the gene universe , i.e. all 47231 genes
ilmnToEntrez_universe <- select(illuminaHumanv4.db, keys=all_genes, columns="ENTREZID", keytype="PROBEID")
entrezGenes_all <- ilmnToEntrez_universe$ENTREZID

library(org.Hs.eg.db)
library(clusterProfiler)

# omit the Null values from the vectors
entrezGenes <- na.omit(entrezGenes)
entrezGenes_all <- na.omit(entrezGenes_all)


# Performing enrichment analysis using enrichGO function
goEnrichment <- enrichGO(gene = entrezGenes, 
                         universe = unique(entrezGenes_all), 
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


