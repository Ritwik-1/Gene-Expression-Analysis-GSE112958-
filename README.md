GSE112958 - CTRL group vs EF(Enteric fever) : Microarray data analysis

There are 2 code files , one of manual analysis and one with limma package:-
1) GSE112958_Typhoid.R
2) GSE112958_using_limma.R

In 1) MANUAL PROCESSING
 ->The data is downloaded using the getGEO function , then the data attributes are printed 
   like the accession ID , platform ID ,experiment type ,organism.
   Then the pheno data and feature data are extracted from the pData and fData functions

 ->The code then extracts the expression matrix from the data variable.
   I have converted the matrix into a dataframe to make a boxplot before preprocessing the data
   which is attached with this file.
 
 ->Then the dataframe is batch corrected and quantile normalized using appropriate R functions.
   I have made another boxplot showing the global variation.
   Then the dataframe is log2 transformed.
   
  ->The effect of log transformation is that all the intensity values have now come in the range of 10 to 19 instead of a skewed variation and is easy of visualize in     plots.

  ->The data has 3 groups , i.e. CTRL , EF and sEF , but for this experiment I have taken 88 samples of CTRL group and 50 samples of EF group 
    Hence the next step in the code is that I have calculated the sample counts of CTRL and 
    EF using a for loop and then modified the dataframe to contain the first 88 columns as CTRL group and the next 50 samples to be EF group.

  ->Now , I have simply traversed all the 47231 transcripts in the dataframe and grouped the 2 
    groups intensity values in two vectors to perform a simple t-test between the two.
    Alongside I have stored the p-values and the fold chain chain values for each transcript in seperate vectors.
  ->The pvalues are then Holm corrected using the function in R and then,
    the holm-pvalues and the fold chain vectors are then -log10 and log2 transformed respectively
    and cbinded with the gene transcript vector.

Then I have used these holm-pvalues and the logFC vectors to draw a volcano plot and chosen the significance values as:-

pvalue : 0.001 
(allowed only 0.1% false positives as we are testing our multiple hypothesis on 47231 transcripts)

logFC > 0.05 and logFC < -0.05 are the limits chosen for the x-axis
as the max(logFC) = 0.54  and min(logFC) = -0.18
this gives us the insight about the actual DEG's may lie in this range 


Keeping these significance thresholds I have selected the geneIDs from the 
"transcript_pv_fc_df" data frame and stored them in a vector.
There are 953 gene transcripts that are the DEGs coming out from these thresholds.


ENRICHMENT ANALYSIS:-

1)This final vector is used for the functional enrichment analysis using the 
  "clusterProfiler" package

2)The gene IDs are initially not in the entrez id format which is needed by the 
  enrichGO function hence the ILMN ids are converted to the entrez IDs for both the 
  gene of interests and the gene universe(using illuminaHumanv4.db package)

3)Then simply use the enrichGO function to get the enriched GO terms from the gene of 
  interests and the gene universe and convert the result into a dataframe
  The pvalue cutoff and qvalue cutoff is taken to be standard 0.05 and 0.1 here
  as now we have got significantly less number of genes as compared to the original set.

4) The dataframe contains alot of information like:-
   The GO ID
   The description of the pathway
   The gene ratio : i.e the proportion of gene present in the go terms that are in the gene of interests
   The Bgratio :  the proportion of genes in the go term that are present in the background set.
   the pvalue , adjusted p values and the qvalues
   
   
5) Then this dataframe is used to select the significant go terms based on the pvalues
   and a bar plot and a similar dot plot is used to visualize the final results
   
6) upon looking to the bar and the dot plot we can observe that some of the most significant pathways
   are:- 
   1) cytoplasmic translation (very high significance)
   2) ribosomal large subunit biogenesis
   3) ribosome biogenesis
   4) ribonucleoprotein complex biogenesis
   5) rRNA metabolic process
   
   All these processes are linked to the cellular processes
   
   eg :Ribosomal large subunit biogenesis refers to the process by which the large subunit of ribosomes is assembled in cells. The large subunit is made up of several rRNA molecules and a large number of ribosomal proteins, which are assembled in a highly coordinated manner.
   
   
OBSERVATION:-

The common focus is the process of translation in the cells and in details it will be linked with the
ribosomal biogenesis. As the process of mRNA translation happens within the ribosomes so it may be the
fact that the ribosomal subunits in the cell are not being able to assemble fully(because of the pathway 2,3,4,5 problems) which is affecting the overall cellular translation process of the cells.


RESULT:-
This may be a pathway that is being disrupted by in the cells when the salmonella typhi enters , which 
becomes the cause of differential expression in the cells of normal(CTRL group) and the EF(enteric fever group).


2) Using the limma package :-
   The data is downloaded in the same way and overall the code is doing the same thing with the
   limma package functions
   we create a design matrix here and use it in the lmfit function from limma package to asign the 
   normalized data a linear model fit
   Then the ebayes function is used to do the differential analysis directly and the parameter values
   are selected in the all_genes dataframe where the pvalues are corrected by holm correction.
   These genes in a similar way are used to make a box plot 
   where the -log(pvalue) cutoff is taken as 201.301 and the log(FC) value to be > 6, this is by
   looking at the max and min values of the respective values for the 47231 transcripts and which gives
   us the desired significant genes.
   The genes retreived are then in a similar way as manual used in the enrichment analysis and bar and 
   the scatter plots are plotted.


