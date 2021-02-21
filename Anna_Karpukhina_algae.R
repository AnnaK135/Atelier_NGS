###############################################################################
# Anna Karpukhina <anna.karpukhina@universite-paris-saclay.fr>
# Analyses statistiques de données génomiques
# AMI2B 2020
###############################################################################

#------------------------------------------------------------------------------
# This is a script performing differential expression analysis of
# the RNA-seq data for the marine green alga Ostreococcus tauri.
#
# The source publication of the data analysed:
# doi:10.1186/s12864-016-2666-6.
#
# The script was tested in R version 3.6.2

#------------------------------------------------------------------------------
# Loading the necessary libraries 
#------------------------------------------------------------------------------
library("DESeq2") # for differential expression analysis
library("mixOmics") # for distance map
library("FactoMineR") # for PCA
library("factoextra") # for PCA visualisation
library("RColorBrewer") # for pretty coloring
library("pheatmap") # for heatmaps 
library("VennDiagram") # for intersecting DEG lists and plotting diagramms

#------------------------------------------------------------------------------
# Loading the data 
#------------------------------------------------------------------------------

setwd("~/M2 materials/Analyse omique/TP_finale") #set you working directory

sampleInfo = read.table("sample_info.txt", header = T, row.names = 1)

#loading the table with raw counts
counts = read.table("mapping_rawdata_allGenes.txt", header = T, row.names = 1)
dim(counts)
#[1] 7699   48
counts <- counts[, 2:48] #removing the column 'length'
# there are actually 47 samples, not 48 (one triplicate of S7 missing)

#------------------------------------------------------------------------------
#######################   Exploring the data##################################
#------------------------------------------------------------------------------

# it is first useful to explore the data to detect possible problems

# -----looking at the distribution of the counts-------------------------------
hist(counts)
# we see that the distributions are strongly skewd to the left
# to make the graphs more visible we will log-transform the data
# the log-transformed counts are expected to be more or less normally distributed

#------------------------------------------------------------------------------
# !(this transformation is only for exploratry purposes, for the diffrential 
#  expression amalysis the rlog transformation embedded in the DESeq2 pipeline
#  will be used)
#------------------------------------------------------------------------------

logcounts=log(counts+1) # adding a pseudocount +1, since for some samples count values for a gene can be zero
hist(logcounts)
#the count distribution seems OK

# -----checking the library sizes ---------------------------------------------
colSums(counts) #shows the total number of counts for all genes for each sample
max(colSums(counts))
# 11010692  --- maximum library size
which.max(colSums(counts))
# HCA.29  --- the sample with max library size
min(colSums(counts))
# 2387232  --- minimum library size 
which.min(colSums(counts))
# HCA.11   --- the sample with min library size
boxplot(logcounts, col= c(rep("lightsalmon3", times=23), rep("aquamarine4", times=24))) #the distribution of pseoudocounts in each sample
legend(1,13, legend = c("ST", "LT" ), fill = c("lightsalmon3", "aquamarine4"))
# the library sizes are compatible
# though they are a little more variable for LT condition

# -----filtering for low counts ------------------------------------------------
vect_nb=rowSums(counts) #the vector with total number of counts for each gene in all samples 
table(vect_nb>0)
#FALSE  TRUE 
#4  7695        --- 4 genes with total 0 read count -> will remove them

# we might also want to remove the rows, where the number of counts is non zero, but extremely low
table(vect_nb>40) #we have 47 samples -> 40 is an adequate threshold
#FALSE  TRUE 
#72  7627       --- these 72 genes will be removed

counts_filtered=counts[which(vect_nb>40),]
dim(counts_filtered)
# 7627   47
dim(counts)
# 7699   47

# in general, the number of genes with low counts was not large, which is good
# we have not removed much -> the distributions should not change a lot
# we can proceed with the further analysis 


#-----------------------------------------------------------------------------------------------------
# Exploring the similarity between the samples
#-----------------------------------------------------------------------------------------------------

# this is useful to have and idea about how the samples cluster together
# and detect possible outliers/batch effects (together with the PCA analysis which will be done later)

#calculating sample-to-sample distance
mat.dist = logcounts  #using log-transformed counts
colnames(mat.dist) = paste(colnames(mat.dist), sampleInfo$Name, sep = " : ")
mat.dist = as.matrix(dist(t(mat.dist))) #calculating Euclidean distances between the samples
mat.dist = mat.dist/max(mat.dist) #normalizing
cim(mat.dist, symkey = FALSE, margins = c(9, 9)) #plotting the heatmap of distances
# blue colors (-> 0 values of normalized distance) indicate that the two samples are more similar in terms of gene expression
# red colors (-> 1 values of normalized distance) indicate that the two samples are more differernt

# in general, samples within one triplicate group are similar and cluster together
# this is good, no outliers are detected

#---------------------------------------------------------------------------------------------------------
# Building the DESeq2 dataset
#---------------------------------------------------------------------------------------------------------

colData <- data.frame(sampleInfo)
#creating the dataset
dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = colData,
                              design = ~ Light + Time + Condition + Iron)

#performing differerntial expression analysis
dds <- DESeq(dds) 
# we can now extract the results of the analysis, but first we will do the PCA (using the dds object)

#--------------------------------------------------------------------------------------------------------
# PCA analysis 
#--------------------------------------------------------------------------------------------------------

rld <- rlog(dds) #rlog transformation to compensate for the effect of extremely high or low count values
# (more comments on rlog transformation in the report)

# first we do the built-in DESeq2 PCA analysis
# (DESeq2 plotPCA takes into account only the most variable genes! -> results may differ from those obtained using other PCA tools)
plotPCA(rld, intgroup=c("Iron", "Light", "Time", "Condition")) #coloring the samples by bioiological replicate
# we see that triplicates cluster together well, there are no evident outliers -> good

plotPCA(rld, intgroup="Condition") #coloring by condition
plotPCA(rld, intgroup="Iron") #coloring by +/- iron
plotPCA(rld, intgroup="Light") # colorong by light

# we observe, that the PC1 and PC2 best distnguish the samples when contrasted by "Light"
# -> "Ligth" probably accounts for most of the variability in the data 

# to explore the influence of each of the three factors in more detail
# another PCA will be done with FactoMineR package

resPCA <- PCA(t(assay(rld)),graph=FALSE) #using the rlog transformed counts obtained before
barplot(resPCA$eig[,2],main="Variance explained by all PCs", ylab ="explained variance, %")
# we see that most of the variance is explained by the first ~10 components
# to obtain a more visually-friendly graph, we will show only the first 15
barplot(resPCA$eig[,2][1:15],main="Variance explained by the first 15 PCs", ylab ="explained variance, %", col="aquamarine4")

#plotting PC1 and PC2
#g1 <- fviz_pca_ind(resPCA, habillage = colData(dds)$Light, addEllipses=TRUE, geom=c("point", "text"), ellipse.level=0.8, palette = "Dark2", repel = "TRUE", title ="Light", pointsize = 2.5) #with labels
g1_light <- fviz_pca_ind(resPCA, habillage = colData(dds)$Light, addEllipses=TRUE, geom="point", ellipse.level=0.85, palette = "Dark2", repel = "TRUE", title ="Light", pointsize = 3) #without labels for a tidier graph
g1_light +theme(axis.text=element_text(size=12), axis.title=element_text(size=14),title=element_text(size=14))
# in PC1-PC2 coordinate space the samples cluster well by "Light", as observed before
# PC1 and PC2 explain most of the data variance (~52%) -> "Light' is likely the most influencial factor 

g1_iron <- fviz_pca_ind(resPCA, habillage = colData(dds)$Iron, addEllipses=TRUE, geom="point", ellipse.level=0.85, palette = "Dark2", repel = "TRUE", title ="Iron", pointsize = 3)
g1_iron +theme(axis.text=element_text(size=12), axis.title=element_text(size=14),title=element_text(size=14))
# no good clustering by "Iron"

g1_condition <- fviz_pca_ind(resPCA, habillage = colData(dds)$Condition, addEllipses=TRUE, geom="point", ellipse.level=0.85, palette = "Dark2", repel = "TRUE", title ="Condition", pointsize = 3) 
g1_condition +theme(axis.text=element_text(size=12), axis.title=element_text(size=14),title=element_text(size=14))
# no good clustering by "Condition"

g1_time <- fviz_pca_ind(resPCA, habillage = colData(dds)$Time, addEllipses=TRUE, geom="point", ellipse.level=0.85, palette = "Dark2", repel = "TRUE", title ="Time", pointsize = 3) 
g1_time +theme(axis.text=element_text(size=12), axis.title=element_text(size=14),title=element_text(size=14))
# good clustering by "Time"

#plotting PC2 and PC3 
g2_light <- fviz_pca_ind(resPCA, axes=c(2,3), habillage = colData(dds)$Light, addEllipses=TRUE, geom="point", ellipse.level=0.85, palette = "Dark2", repel = "TRUE", title ="Light", pointsize = 3)
g2_light+theme(axis.text=element_text(size=12), axis.title=element_text(size=14),title=element_text(size=14))
# no good clustering by "Light'

g2_iron <- fviz_pca_ind(resPCA, axes=c(2,3), habillage = colData(dds)$Iron, addEllipses=TRUE, geom="point", ellipse.level=0.85, palette = "Dark2", repel = "TRUE", title ="Iron", pointsize = 3)
g2_iron+theme(axis.text=element_text(size=12), axis.title=element_text(size=14),title=element_text(size=14))
# no good clustering by "Iron'

g2_condition <- fviz_pca_ind(resPCA, axes=c(2,3), habillage = colData(dds)$Condition, addEllipses=TRUE, geom="point", ellipse.level=0.85, palette = "Dark2", repel = "TRUE", title ="Condition", pointsize = 3)
g2_condition+theme(axis.text=element_text(size=12), axis.title=element_text(size=14),title=element_text(size=14))
# in the PC2-PC3 coordinates the samples cluster well by "Condition'
# -> "Condition" is likely the 2nd most influential factor (PC2 + PC3 ~32.1% of variance)

#plotting PC3 and PC4
g3_iron <- fviz_pca_ind(resPCA, axes=c(3,4), habillage = colData(dds)$Iron, addEllipses=TRUE, geom="point", ellipse.level=0.85, palette = "Dark2", repel = "TRUE", title ="Iron", pointsize = 3)
g3_iron+theme(axis.text=element_text(size=12), axis.title=element_text(size=14),title=element_text(size=14))
# better clustering by "Iron', but not good yet 

#plotting PC3 and PC3
g4_iron <- fviz_pca_ind(resPCA, axes=c(4,5), habillage = colData(dds)$Iron, addEllipses=TRUE, geom="point", ellipse.level=0.85, palette = "Dark2", repel = "TRUE", title ="Iron", pointsize = 3)
g4_iron+theme(axis.text=element_text(size=12), axis.title=element_text(size=14),title=element_text(size=14))
# Finally clustering by "Iron'!

# Based on the PCA analysis, we can suggest that the importance of the variables on gene expression
# decreases as follows: Light > Condition > Iron

#------------------------------------------------------------------------------------------------------
# Proceeding with the differential expression analysis 
#-------------------------------------------------------------------------------------------------------

# calling the results without specifying the contrast will produce the comparison 
# contrasting the last variable in our design formula, which is "Iron"
res <- results(dds, contrast=c("Iron", "NO", "YES"))
head(res)
#log2 fold change (MLE): Iron NO vs YES 
#Wald test p-value: Iron NO vs YES 
# ....

summary(res)
#out of 7627 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1305, 17%
#LFC < 0 (down)     : 1500, 20%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
#(mean count < 1)
res <- res[order(res$padj),]  #order by adjusted p-value
hist(res$pvalue) #no peaks near large p-values -> OK
write.csv(res, "FeNO_vs_FeYes.csv") #saving the results to Excel table

#these are the results contrasting Fe+ vs Fe- no matter all the other conditions

#if we want to extract the effect of Fe+ vs Fe- specifically for each sample group,
# we need to change the design formula including the group factor

# --- performing the analysis separately for the two Conditions

dds_cond <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = colData,
                              design = ~ Light + Time + Condition + Iron)

dds_cond$group <- factor(paste0(dds_cond$Condition, dds_cond$Iron))
design(dds_cond) <- ~ group
dds_cond <- DESeq(dds_cond)
resultsNames(dds_cond) 
#[1] "Intercept"  "group_LTYES_vs_LTNO" "group_STNO_vs_LTNO"  "group_STYES_vs_LTNO"

res_LT <- results(dds_cond, contrast=c("group", "LTNO", "LTYES"))
head(res_LT)
summary(res_LT)
#LFC > 0 (up)       : 534, 7%
#LFC < 0 (down)     : 430, 5.6%
res_LT <- res_LT[order(res_LT$padj),] #order by adjusted p-value
hist(res_LT$pvalue) #we expect ~uniform distribution between 0 and 1 with a peak near 0
write.csv(res_LT, "LT_FeNO_vs_FeYes.csv")

res_ST <- results(dds_cond, contrast=c("group", "STNO", "STYES"))
head(res_ST)
summary(res_ST)
#LFC > 0 (up)       : 103, 1.3%
#LFC < 0 (down)     : 22, 0.29%
#outliers [1]       : 0, 0%
#low counts [2]     : 296, 3.9%
res_ST <- res_ST[complete.cases(res_ST),]  #remove any rows with NA (row counts and outliers)
summary(res_ST)
#LFC > 0 (up)       : 103, 1.4%
#LFC < 0 (down)     : 22, 0.3%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
res_ST <- res_ST[order(res_ST$padj),] #order by adjusted p-value
hist(res_ST$pvalue) 
write.csv(res_ST, "ST_FeNO_vs_FeYes.csv")


# --- performing the analysis separately for each time point in each condition

dds_cond_time <- DESeqDataSetFromMatrix(countData = counts_filtered,
                                   colData = colData,
                                   design = ~ Light + Time + Condition + Iron)


dds_cond_time$group <- factor(paste0(dds_cond_time$Time, dds_cond_time$Condition, dds_cond_time$Iron))
design(dds_cond_time) <- ~ group
dds_cond_time <- DESeq(dds_cond_time)
resultsNames(dds_cond_time)

res_ST3H <- results(dds_cond_time, contrast=c("group", "3HSTNO", "3HSTYES"))
head(res_ST3H)
summary(res_ST3H)
#LFC > 0 (up)       : 66, 0.87%
#LFC < 0 (down)     : 30, 0.39%
res_ST3H <- res_ST3H[order(res_ST3H$padj),] #order by adjusted p-value
hist(res_ST3H$pvalue) # not so good
write.csv(res_ST3H, "ST3H_FeNO_vs_FeYes.csv")


res_ST6H <- results(dds_cond_time, contrast=c("group", "6HSTNO", "6HSTYES"))
head(res_ST6H)
summary(res_ST6H)
res_ST6H <- res_ST6H[complete.cases(res_ST6H),] #removing low counts
summary(res_ST6H)
#LFC > 0 (up)       : 160, 2.2%
#LFC < 0 (down)     : 184, 2.5%
res_ST6H <- res_ST6H[order(res_ST6H$padj),] #order by adjusted p-value
hist(res_ST6H$pvalue) 
write.csv(res_ST6H, "ST6H_FeNO_vs_FeYes.csv")


res_LT3H <- results(dds_cond_time, contrast=c("group", "3HLTNO", "3HLTYES"))
head(res_LT3H)
summary(res_LT3H)
#LFC > 0 (up)       : 372, 4.9%
#LFC < 0 (down)     : 132, 1.7%
res_LT3H <- res_LT3H[order(res_LT3H$padj),] #order by adjusted p-value
hist(res_LT3H$pvalue) 
write.csv(res_LT3H, "LT3H_FeNO_vs_FeYes.csv")

res_LT9H <- results(dds_cond_time, contrast=c("group", "9HLTNO", "9HLTYES"))
head(res_LT9H)
summary(res_LT9H)
res_LT9H <- res_LT9H[complete.cases(res_LT9H),] #removing low counts
summary(res_LT9H)
#LFC > 0 (up)       : 291, 3.9%
#LFC < 0 (down)     : 159, 2.1%
res_LT9H <- res_LT9H[order(res_LT9H$padj),] #order by adjusted p-value
hist(res_LT9H$pvalue) 
write.csv(res_LT9H, "LT9H_FeNO_vs_FeYes.csv")


res_LT15H <- results(dds_cond_time, contrast=c("group", "15HLTNO", "15HLTYES"))
head(res_LT15H)
summary(res_LT15H)
res_LT15H <- res_LT15H[complete.cases(res_LT15H),] #removing low counts
summary(res_LT15H)
#LFC > 0 (up)       : 482, 6.4%
#LFC < 0 (down)     : 124, 1.7%
res_LT15H <- res_LT15H[order(res_LT15H$padj),] #order by adjusted p-value
hist(res_LT15H$pvalue) 
write.csv(res_LT15H, "LT15H_FeNO_vs_FeYes.csv")


res_LT22H <- results(dds_cond_time, contrast=c("group", "22HLTNO", "22HLTYES"))
head(res_LT22H)
summary(res_LT22H)
#LFC > 0 (up)       : 864, 11%
#LFC < 0 (down)     : 312, 4.1%
res_LT22H <- res_LT22H[order(res_LT22H$padj),] #order by adjusted p-value
hist(res_LT22H$pvalue) 
write.csv(res_LT22H, "LT22H_FeNO_vs_FeYes.csv")

#-----------------------------------------------------------------------------------------------------------
# Result filtering and visualisation
#-----------------------------------------------------------------------------------------------------------

# we will apply the tresholds of padj < 0.01 and |log2FC| > 1 
# the genes meeting these criteria will be considered significantly differentially expressed 

# --- Global FeNo vs FeYes comparison ---------------------------------------------------------------------

# significant (padj < 0.01) up-regulated genes
up <- res[ res[,'log2FoldChange'] > 1, ]  
up <- up[ up[,'padj'] < 0.01, ] 
dim(up)
# 95 genes

#significant (padj < 0.01) down-regulated
down <- res[ res[,'log2FoldChange'] < -1, ]  
down <- down[ down[,'padj'] < 0.01, ] 
dim(down)
# 2 genes

updown <- c(row.names(up), row.names(down))
write.csv(res[updown,], "all_sign.csv") 


# ---volcano plot
# Summarizes the fold changes of the differentially expressed genes together with their statistical significance

alpha <- 0.01 # Threshold on the adjusted p-value (horizontal line)
cols <- densCols(res$log2FoldChange, -log10(res$padj))
plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),
     main="Fe -   vs   Fe +", xlab="log2 fold change", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.8)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

# The interesting genes are in upper right and left corners
# Much more genes are upregulated in Fe-poor condition

# ---heatmaps

# for significant upregulated genes
select_up =(assay(rld)[row.names(up), ])
pheatmap(select_up)
pheatmap(assay(rld)[select_up, ], scale = "row")
# for significant downregulated genes
select_down =(assay(rld)[row.names(down), ])
pheatmap(select_down)
pheatmap(assay(rld)[select_down, ], scale = "row")


# --- FeNo vs FeYes within Condition1 and Condition2---------------------------------------------------

# significant (padj < 0.01) up-regulated genes in Condition 1 (ST)
up_ST <- res_ST[ res_ST[,'log2FoldChange'] > 1, ]  
up_ST <- up_ST[ up_ST[,'padj'] < 0.01, ] 
dim(up_ST)
# 14 genes

# significant (padj < 0.01) down-regulated in Condition 1 (ST)
down_ST <- res_ST[ res_ST[,'log2FoldChange'] < -1, ]  
down_ST <- down_ST[ down_ST[,'padj'] < 0.01, ] 
dim(down_ST)
# 2 genes

updownST <- c(row.names(up_ST), row.names(down_ST))
write.csv(res_ST[updownST,], "ST_sign.csv") 

# significant (padj < 0.01) up-regulated genes in Condition 2 (LT)
up_LT <- res_LT[ res_LT[,'log2FoldChange'] > 1, ]  
up_LT <- up_LT[ up_LT[,'padj'] < 0.01, ] 
dim(up_LT)
# 155 genes

# significant (padj < 0.01) down-regulated in Condition 2 (LT)
down_LT <- res_LT[ res_LT[,'log2FoldChange'] < -1, ]  
down_LT <- down_LT[ down_LT[,'padj'] < 0.01, ] 
dim(down_LT)
# 13 genes

updownLT <- c(row.names(up_LT), row.names(down_LT))
write.csv(res_LT[updownLT,], "LT_sign.csv") 

# ---volcano plots 
par(mfrow = c(1, 2))

alpha <- 0.01 # Threshold on the adjusted p-value (horizontal line)
cols <- densCols(res_ST$log2FoldChange, -log10(res_ST$padj))
plot(res_ST$log2FoldChange, -log10(res_ST$padj), col=cols, panel.first=grid(),
     main="Short term response: Fe - vs Fe +", xlab="log2 fold change", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.8)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

alpha <- 0.01 # Threshold on the adjusted p-value (horizontal line)
cols <- densCols(res_LT$log2FoldChange, -log10(res_LT$padj))
plot(res_LT$log2FoldChange, -log10(res_LT$padj), col=cols, panel.first=grid(),
     main="Long term response: Fe - vs Fe +", xlab="log2 fold change", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.8)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

# we observe,that more genes are differentially expressed in response to Fe in Condition 2 (LT)
# more genes are upregulated in both conditions

# --- FeNo vs FeYes within Condition1 and Condition2 and different time points ------------------------

# significant (padj < 0.01) up-regulated genes in Condition 1 (ST), 3H 
up_ST3H <- res_ST3H[ res_ST3H[,'log2FoldChange'] > 1, ]  
up_ST3H <- up_ST3H[ up_ST3H[,'padj'] < 0.01, ] 
dim(up_ST3H)
# 13 genes

# significant (padj < 0.01) down-regulated in Condition 1 (ST), 3H 
down_ST3H <- res_ST3H[ res_ST3H[,'log2FoldChange'] < -1, ]  
down_ST3H <- down_ST3H[ down_ST3H[,'padj'] < 0.01, ] 
dim(down_ST3H)
# 1 gene

updownST3H <- c(row.names(up_ST3H), row.names(down_ST3H))
write.csv(res_ST3H[updownST3H,], "ST3H_sign.csv") 

# significant (padj < 0.01) up-regulated genes in Condition 1 (ST), 6H 
up_ST6H <- res_ST6H[ res_ST6H[,'log2FoldChange'] > 1, ]  
up_ST6H <- up_ST6H[ up_ST6H[,'padj'] < 0.01, ] 
dim(up_ST6H)
# 18 genes

# significant (padj < 0.01) down-regulated in Condition 1 (ST), 6H 
down_ST6H <- res_ST6H[ res_ST6H[,'log2FoldChange'] < -1, ]  
down_ST6H <- down_ST6H[ down_ST6H[,'padj'] < 0.01, ] 
dim(down_ST6H)

updownST6H <- c(row.names(up_ST6H), row.names(down_ST6H))
write.csv(res_ST6H[updownST6H,], "ST6H_sign.csv") 

# significant (padj < 0.01) up-regulated genes in Condition 2 (LT), 3H 
up_LT3H <- res_LT3H[ res_LT3H[,'log2FoldChange'] > 1, ]  
up_LT3H <- up_LT3H[ up_LT3H[,'padj'] < 0.01, ] 
dim(up_LT3H)
# 136 genes

# significant (padj < 0.01) down-regulated in Condition 2 (LT), 3H 
down_LT3H <- res_LT3H[ res_LT3H[,'log2FoldChange'] < -1, ]  
down_LT3H <- down_LT3H[ down_LT3H[,'padj'] < 0.01, ] 
dim(down_LT3H)
# 7 genes

updownLT3H <- c(row.names(up_LT3H), row.names(down_LT3H))
write.csv(res_LT3H[updownLT3H,], "LT3H_sign.csv") 

# significant (padj < 0.01) up-regulated genes in Condition 2 (LT), 9H 
up_LT9H <- res_LT9H[ res_LT9H[,'log2FoldChange'] > 1, ]  
up_LT9H <- up_LT9H[ up_LT9H[,'padj'] < 0.01, ] 
dim(up_LT9H)
# 109 genes

# significant (padj < 0.01) down-regulated in Condition 2 (LT), 9H 
down_LT9H <- res_LT9H[ res_LT9H[,'log2FoldChange'] < -1, ]  
down_LT9H <- down_LT9H[ down_LT9H[,'padj'] < 0.01, ] 
dim(down_LT9H)
# 10 genes

updownLT9H <- c(row.names(up_LT9H), row.names(down_LT9H))
write.csv(res_LT9H[updownLT9H,], "LT9H_sign.csv")

# significant (padj < 0.01) up-regulated genes in Condition 2 (LT), 15H 
up_LT15H <- res_LT15H[ res_LT15H[,'log2FoldChange'] > 1, ]  
up_LT15H <- up_LT15H[ up_LT15H[,'padj'] < 0.01, ] 
dim(up_LT15H)
# 142 genes

# significant (padj < 0.01) down-regulated in Condition 2 (LT), 15H 
down_LT15H <- res_LT15H[ res_LT15H[,'log2FoldChange'] < -1, ]  
down_LT15H <- down_LT15H[ down_LT15H[,'padj'] < 0.01, ] 
dim(down_LT15H)
# 2 genes

updownLT15H <- c(row.names(up_LT15H), row.names(down_LT15H))
write.csv(res_LT15H[updownLT15H,], "LT15H_sign.csv")

# significant (padj < 0.01) up-regulated genes in Condition 2 (LT), 22H 
up_LT22H <- res_LT22H[ res_LT22H[,'log2FoldChange'] > 1, ]  
up_LT22H <- up_LT22H[ up_LT22H[,'padj'] < 0.01, ] 
dim(up_LT22H)
# 326 genes

# significant (padj < 0.01) down-regulated in Condition 2 (LT), 22H 
down_LT22H <- res_LT22H[ res_LT22H[,'log2FoldChange'] < -1, ]  
down_LT22H <- down_LT22H[ down_LT22H[,'padj'] < 0.01, ] 
dim(down_LT22H)
# 43 genes

updownLT22H <- c(row.names(up_LT22H), row.names(down_LT22H))
write.csv(res_LT22H[updownLT22H,], "LT22H_sign.csv")

# --- plotting the number of significant DEGs at fifferent time points and different conditions -----------
par(mfrow = c(1, 1))
cond1 <- c(dim(up_ST3H)[1], dim(down_ST3H)[1], dim(up_ST6H)[1], dim(down_ST6H)[1])
cond2 <- c(dim(up_LT3H)[1], dim(down_LT3H)[1], dim(up_LT9H)[1], dim(down_LT9H)[1], dim(up_LT15H)[1], dim(down_LT15H)[1], dim(up_LT22H)[1], dim(down_LT22H)[1])
c <- c(cond1, cond2)
col = c("lightsalmon3", "lightsalmon", "lightsalmon3", "lightsalmon", "aquamarine4", "mediumaquamarine", "aquamarine4", "mediumaquamarine", "aquamarine4", "mediumaquamarine", "aquamarine4", "mediumaquamarine")
barplot(c, col=col, main = "Numbers of DEGs with padj < 0.01 and |log2Fold| > 1")
legend(1,300, legend = c("ST UP", "ST DOWN", "LT UP", "LT DOWN" ), fill = c("lightsalmon3", "lightsalmon", "aquamarine4", "mediumaquamarine"))

# --- intersecting the DEGs -----------------------------------------------------------------------------

par(mfrow = c(1, 1))
# upregulated DEGs
LT_degs_up <- c(row.names(up_LT3H), row.names(up_LT9H), row.names(up_LT15H), row.names(up_LT22H))
ST_degs_up <- c(row.names(up_ST3H), row.names(up_ST6H))
degs_up_intersect = intersect(LT_degs_up, ST_degs_up)
length(degs_up_intersect)
# 17 genes are upregulated in both Conditions (no matter the time point)
# extracting these genes
res_LT[degs_up_intersect,]
write.csv(res_LT[degs_up_intersect,], "up_degs_intersect.csv")

# upregulated DEGs
LT_degs_down <- c(row.names(down_LT3H), row.names(down_LT9H), row.names(down_LT15H), row.names(down_LT22H))
ST_degs_down <- c(row.names(down_ST3H), row.names(down_ST6H))
degs_down_intersect = intersect(LT_degs_down, ST_degs_down)
length(degs_down_intersect)
# 1 gene downregulated in both Conditions (no matter the time point)
res_LT[degs_down_intersect,]
write.csv(res_LT[degs_down_intersect,], "down_degs_intersect.csv")

# --- plotting Venn diagramms 

deg.venn_up <- list('intersect' = length(degs_up_intersect),
                 'Long-term response' = length(LT_degs_up),
                 'Short-term response' = length(ST_degs_up))

venn.plot_up <- draw.pairwise.venn(deg.venn_up$`Long-term response`, deg.venn_up$`Short-term response`, deg.venn_up$intersect,
                                category = c('Long-term response', 'Short-term response'), scaled = T,
                                fill = c("aquamarine4", "lightsalmon3"), alpha = rep(0.5, 2),
                                cat.pos = c(0, 0), cex = 2, cat.cex = 1.5,
                                fontfamily = "Arial",
                                cat.fontfamily = "Arial",
                                ext.text = FALSE)

deg.venn_down <- list('intersect' = length(degs_down_intersect),
                    'Long-term response' = length(LT_degs_down),
                    'Short-term response' = length(ST_degs_down))

venn.plot_down <- draw.pairwise.venn(deg.venn_down$`Long-term response`, deg.venn_down$`Short-term response`, deg.venn_down$intersect,
                                   category = c('Long-term response', 'Short-term response'), scaled = T,
                                   fill = c("aquamarine4", "lightsalmon3"), alpha = rep(0.5, 2),
                                   cat.pos = c(0, 0), cex = 2, cat.cex = 1.5,
                                   fontfamily = "Arial",
                                   cat.fontfamily = "Arial",
                                   ext.text = FALSE)


# These 17 + 1 genes which differerntially expressed in response to Fe starvation are probably the most
# interesting candidates to be involved in O.tauri's cellular adaptation to low Fe concentrations.
interesting <- c(degs_up_intersect, degs_down_intersect)

# These genes will be further analysed in the report.

# Of course, this criteria is rather stringent and other DEGs could also be explored 
# (eg. those, that are differentially expressed in only some of the groups or those with an interesting known function)
# But since there is not enough time to explore all of that till deadline, the analysis will stop here :)

# plotting the heatmaps for onteresting genes
select =(assay(rld)[interesting, ])
pheatmap(select)
pheatmap(select, scale = "row") #centering by row



# --------------------------------------------------------------------------------------------------------
# End of script. Thank you for reading! 
# --------------------------------------------------------------------------------------------------------










