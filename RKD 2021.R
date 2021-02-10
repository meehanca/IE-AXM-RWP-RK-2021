####################################################################
#Library
####################################################################
library(DESeq2)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(VennDiagram)
library(ggrepel)
library(rhdf5)
library(tximport)

####################################################################
# Set work environment and convert gene IDs to transcript IDs for txi import to DESeq2
####################################################################

setwd("~/Documents/Wheat/Transcriptomics_analysis/RKD_2021/data/3_aligned")
htseqDir<-getwd()
samples <- read.table(file='sampleTable.txt', header=T)

##  read in pre-constructed tx2gene table (transcript to gene table)
##  If putting new cDNA into library for kallisto remember to do this in the transcripts to genes file too

tx2gene <- read.table("../../../philippa_resources/transcripts_to_genes_RefSeqv1.0_annot_v1.1.txt", header=T)
head(tx2gene)

files <- c("YL104/abundance.tsv","YL107/abundance.tsv","YL117/abundance.tsv","YL118/abundance.tsv",
           "YL126/abundance.tsv","YL127/abundance.tsv","YL130/abundance.tsv","YL131/abundance.tsv","YL133/abundance.tsv")
names(files) <- c("IE_TM1","IE_TM2","IE_WE1","IE_WE2",
                  "Me_WN1","Me_WN2","Me_TN1","Me_TN2","Me_TN4")

##  If changing pseudoaligner i.e. Salmon, change type=""

txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
head(txi.kallisto.tsv$counts)

Quality_fastp <- read.table('Quality_fastp.txt', sep='\t',header = T)
####################################################################
# DESeq2 importing and normalising 
####################################################################

##  Read in the results from the analysis (the counts files)

ddsHTSeq<-DESeqDataSetFromTximport(txi=txi.kallisto.tsv,
                                   colData=samples,
                                   design =~condition)

##  design <- test everything in relations to condition
##  if you have more than one conditions you want to differentiate (for example different genotypes) you change design = ~  condition + genotype
##  And perform the analysis (details in the manual)

##  And perform the analysis (details in the manual)
##  Kallisto gives out TPM values but DESeq2 can only deal with integer count values so it will change them (or take the est. count column)

dds<-DESeq(ddsHTSeq)

## Dispersion plot shows you that the data (counts) have been fitted to a dispersion trend - ensures no outliers (normalisation)

plotDispEsts(dds, main="Dispersion plot")

####################################################################
# Do PCA
####################################################################

##  principal component analysis
##  this function does not allow you to plot PCs other than the first 2 but is much faster than other algorithm for calculating PCs that explain variance (dispersion)

##  All samples
vst = vst(ddsHTSeq)

v <- plotPCA(vst, intgroup=c("condition"))
v<- v+ geom_label_repel(aes(label = name),max.overlaps = 20)
v

##  Immature embryo samples
vst = vst(ddsHTSeq[,c(1,3,4,5)])

v <- plotPCA(vst, intgroup=c("condition"))
v<- v+ geom_label_repel(aes(label = name),max.overlaps = 20)
v

##  Meristem samples
vst = vst(ddsHTSeq[,c(7,8,10,11,12)])

v <- plotPCA(vst, intgroup=c("condition"))
v<- v+ geom_label_repel(aes(label = name),max.overlaps = 20)
v

####################################################################
#Plotting Reps
####################################################################

##  Plotting replicates function to get an idea of how well replicates within conditions match up with one another

plot_reps =  function(dds,x=1,y=2,cond_choice=1, cond='condition'){
  ##  Estimate the size factors for normalisation
  dds<-estimateSizeFactors(dds)
  
  ## Extract the normalised counts for the condition you want
  rep_values<- counts(dds, normalized=TRUE)[,dds[[cond]]==cond_choice]
  
  # Take logs of these values
  vals <- log2(rep_values[,c(x,y)] + 0.5)
  # And plot
  plot(vals,pch=16, cex=0.4,xlab=paste('rep',x),ylab=paste('rep',y))
  grid(col = "darkgray", lty = "solid",lwd = par("lwd"), equilogs = TRUE)
  title(paste("Comparison of",cond_choice,"replicates"))
}

##  This will take a while for bigger genomes (more rows of genes) as it involves plotting all the data points against each other

par(mfrow = c(3,3))

plot_reps(dds, x=1, y=2, cond_choice="TE")
plot_reps(dds, x=1, y=3, cond_choice="TE")
plot_reps(dds, x=2, y=3, cond_choice="TE")

plot_reps(dds, x=1, y=2, cond_choice="TM")
plot_reps(dds, x=1, y=3, cond_choice="TM")
plot_reps(dds, x=2, y=3, cond_choice="TM")

plot_reps(dds, x=1, y=2, cond_choice="NE")
plot_reps(dds, x=1, y=3, cond_choice="NE")
plot_reps(dds, x=2, y=3, cond_choice="NE")

####################################################################
#Immature Embryos DEGs
####################################################################

##  This filter_degs function must be changed whenever the alpha value in the results() function changes as it filters out DEGs from based on p-value

filter_degs <- function(res){
  summary(res)
  res2 = res[!(is.na(res$padj)),]
  res2 = res2[res2$padj < 0.01,]
  return(res2)
}

##  alpha=p.value, lfcThreshold=log fold change, pAdjust..=hypothesis correction
##  only really appropriate to change lfcThreshold in my opinion, maybe pAdjust but p-value is just changing the rate of false positives

resultsNames(dds)
IE_TE_vs_WE_DEGs = results(dds, contrast= c("condition", "IE_TE", "IE_WE"), alpha = 0.01, pAdjustMethod = "BH")
IE_TE_vs_WE_DEG = filter_degs(IE_TE_vs_WE_DEGs)

IE_TE_vs_TM_DEGs = results(dds, contrast= c("condition", "IE_TE", "IE_TM"), alpha = 0.01, pAdjustMethod = "BH")
IE_TE_vs_TM_DEG = filter_degs(IE_TE_vs_TM_DEGs)

IE_TM_vs_WE_DEGs = results(dds, contrast= c("condition", "IE_TM", "IE_WE"), alpha = 0.01, pAdjustMethod = "BH")
IE_TM_vs_WE_DEG = filter_degs(IE_TM_vs_WE_DEGs)

summary(IE_TE_vs_WE_DEG)
summary(IE_TE_vs_TM_DEG)
summary(IE_TM_vs_WE_DEG)



library(EnhancedVolcano)

EnhancedVolcano(IE_TM_vs_WE_DEG,
                lab = rownames(IE_TM_vs_WE_DEG),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5, 12),
                ylim = c(0,60))


####################################################################
#Meristem DEGs
####################################################################

##  This filter_degs function must be changed whenever the alpha value in the results() function changes as it filters out DEGs from based on p-value

filter_degs <- function(res){
  summary(res)
  res2 = res[!(is.na(res$padj)),]
  res2 = res2[res2$padj < 0.05,]
  return(res2)
}

##  alpha=p.value, lfcThreshold=log fold change, pAdjust..=hypothesis correction
##  only really appropriate to change lfcThreshold in my opinion, maybe pAdjust but p-value is just changing the rate of false positives

resultsNames(dds)
ME_TN_vs_WN_DEGs = results(dds, contrast= c("condition", "Me_TN", "Me_WN"), alpha = 0.05, pAdjustMethod = "BH")
ME_TN_vs_WN_DEG = filter_degs(ME_TN_vs_WN_DEGs)

ME_TE_vs_TM_DEGs = results(dds, contrast= c("condition", "Me_TE", "Me_TM"), alpha = 0.01, pAdjustMethod = "BH")
ME_TE_vs_TM_DEG = filter_degs(ME_TE_vs_TM_DEGs)

ME_TE_vs_WN_DEGs = results(dds, contrast= c("condition", "Me_TE", "Me_WN"), alpha = 0.01, pAdjustMethod = "BH")
ME_TE_vs_WN_DEG = filter_degs(ME_TE_vs_WN_DEGs)

ME_TM_vs_WN_DEGs = results(dds, contrast= c("condition", "Me_TE", "Me_WN"), alpha = 0.01, pAdjustMethod = "BH")
ME_TM_vs_WN_DEG = filter_degs(ME_TM_vs_WN_DEGs)

summary(ME_TN_vs_WN_DEG)
summary(ME_TE_vs_TM_DEG)
summary(ME_TE_vs_WN_DEG)
summary(ME_TM_vs_WN_DEG)


library(EnhancedVolcano)

EnhancedVolcano(ME_TN_vs_WN_DEG,
                lab = rownames(ME_TN_vs_WN_DEG),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5, 12),
                ylim = c(0,60))

####################################################################
#Heatmap - Immature embryos
####################################################################

counts = counts(dds , normalized = TRUE)

##  Removes genes that are not express and have no variance, will give error and heatmaps are meant to compare distinct expression profiles 

counts <-  counts[apply(counts, MARGIN = 1, FUN = function(x) sd(x) != 0 ),]
##  Select genes that we are interested in looking at across samples 
heatmap <- (counts[IE_TM_vs_WE_DEG_ID,1:5])

##  Heatmap has to be logarithmically transformed to shrink values that differ by orders of magnitude and because of this be +1 (log(0)=NA)
##  Scale so that rows with similar expression profiles are plotted beside each other

#  cluster_rows = F
#  cluster_cols = F
#  these arguments are by default set as T, if you change them to false interesting things can happen 

pheatmap((log2(heatmap+1)), scale = "row",border_color=NA,show_rownames = F, main = 'IE TM vs WE expression across samples')


####################################################################
#Heatmap - Meristem
####################################################################

counts = counts(dds , normalized = TRUE)

##  Removes genes that are not express and have no variance, will give error and heatmaps are meant to compare distinct expression profiles 

counts <-  counts[apply(counts, MARGIN = 1, FUN = function(x) sd(x) != 0 ),]
##  Select genes that we are interested in looking at across samples 
heatmap <- (counts[ME_TN_vs_WN_DEG_ID,6:10])

##  Heatmap has to be logarithmically transformed to shrink values that differ by orders of magnitude and because of this be +1 (log(0)=NA)
##  Scale so that rows with similar expression profiles are plotted beside each other

#  cluster_rows = F
#  cluster_cols = F
#  these arguments are by default set as T, if you change them to false interesting things can happen 

pheatmap((log2(heatmap+1)), scale = "row",border_color=NA,show_rownames = F, main = 'TN vs WN DEGs expression across samples')

####################################################################
#Counts of potentially interesting genes
####################################################################

library(ideal)
par(mfrow=c(2,2))
ggplotCounts(dds,'TraesCS7A02G135900',transform = F)
ggplotCounts(dds,"TaRKD2",transform = F)
ggplotCounts(dds,"TaRKD1-7B",transform = F)
ggplotCounts(dds,"TaRKD1-7D",transform = F)

## Meristem TN vs WN TF ##
library(gridExtra)
##AP2s
ggplotCounts(dds,"TraesCS1A02G058400",transform = F) #RSR1 A homeologue
ggplotCounts(dds,"TraesCS1A02G314300",transform = F) #OsERF#074 ortholog
ggplotCounts(dds,"TraesCS1B02G076300",transform = F) #RSR1 B homeologue
ggplotCounts(dds,"TraesCS7A02G376300",transform = F) #Os06ERF A homeologue
ggplotCounts(dds,"TraesCS7B02G277800",transform = F) #Os06ERF B homeologue
ggplotCounts(dds,"TraesCS7D02G372700",transform = F) #Os06ERF D homeologue

grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6)

#ARFs
ggplotCounts(dds,"TraesCS6B02G141300",transform = F) #OsARF5
ggplotCounts(dds,"TraesCS7B02G065800",transform = F) #ARF16 - B homeologue
ggplotCounts(dds,"TraesCS7D02G161900",transform = F) #ARF16 - D homeologue

grid.arrange(plot1,plot2,plot3)

ggplotCounts(dds,"TraesCS6A02G286800",transform = F) #SET Domain
ggplotCounts(dds,"TraesCS1D02G155200",transform = F) #OsHox9
ggplotCounts(dds,"TraesCS3B02G470000",transform = F) #MADS51

ggplotCounts(dds,"TraesCS3D02G531900",transform = F) #PCL Photoperiod control
ggplotCounts(dds,"TraesCS4B02G115100",transform = F) #PRR59 Photoperiod control
ggplotCounts(dds,"TraesCS5A02G320300",transform = F) #PRR95 Photoperiod control 

#bHLH
ggplotCounts(dds,"TraesCS3A02G252900",transform = F) #ICE1 - A homeologue
ggplotCounts(dds,"TraesCS3B02G284800",transform = F) #ICE1 - B homeologue
ggplotCounts(dds,"TraesCS3B02G288700",transform = F) #ICE1 - B homeologue paralogue?
ggplotCounts(dds,"TraesCS3D02G253700",transform = F) #ICE1 - D homeologue

## Immature Embryo NIL overlap

ggplotCounts(dds,"TraesCS5B02G349900",transform = F) #Histone H4 B
ggplotCounts(dds,"TraesCS5D02G354100",transform = F) #Histone H4 D
ggplotCounts(dds,"TraesCS4D02G042100",transform = F) #OsUCL8

## Immature Embryo Microspore overlap

ggplotCounts(dds,"TraesCS2A02G168200",transform = F) #SUS3 - sucrose synthatse
ggplotCounts(dds,"TraesCS3B02G068900",transform = F) #MYB methyltransferase regulator

## Immature Embryo TFs

ggplotCounts(dds,"TraesCS7B02G196900",transform = F) #NLT5 - NAC - suppressor of flowering *
ggplotCounts(dds,"TraesCS6B02G248400",transform = F) #OsCOL4 - Constans-like - suppressor of flowering
ggplotCounts(dds,"TraesCS4D02G159600",transform = F) #SUMO (small ubiquitin-related modifier) E3-ligase
ggplotCounts(dds,"TraesCS3A02G432900",transform = F) #MADS51
ggplotCounts(dds,"TraesCS4D02G231900",transform = F) #BBR
ggplotCounts(dds,"TraesCS4D02G083200",transform = F) #NLP5 *

ggplotCounts(dds,"TraesCS7B02G028200",transform = F) #Histone methyltransferase *
ggplotCounts(dds,"TraesCS4B02G054200",transform = F) #Histone H2B
ggplotCounts(dds,"TraesCS1A02G361600",transform = F) #Histone H2B
ggplotCounts(dds,"TraesCS5D02G478400",transform = F) #Histone H2A
ggplotCounts(dds,"TraesCS6A02G059000",transform = F) #Histone H4 x
ggplotCounts(dds,"TraesCS1D02G290000",transform = F) #Histone H4 x
ggplotCounts(dds,"TraesCS4D02G205000",transform = F) #Histone H2A x
ggplotCounts(dds,"TraesCS1B02G378400",transform = F) #Histone H3 x
ggplotCounts(dds,"TraesCS7B02G028200",transform = F) #Histone methyltransferase
ggplotCounts(dds,"TraesCS7B02G028200",transform = F) #Histone methyltransferase
ggplotCounts(dds,"TraesCS7B02G028200",transform = F) #Histone methyltransferase
ggplotCounts(dds,"TraesCS7B02G028200",transform = F) #Histone methyltransferase
ggplotCounts(dds,"TraesCS7B02G028200",transform = F) #Histone methyltransferase
ggplotCounts(dds,"TraesCS7B02G028200",transform = F) #Histone methyltransferase

###################################################################
#TF terms
###################################################################
TFs<- read.delim("~/Documents/R scripts and projects/Projects/RKD_grasses/philippa_resources/TFs/All.txt")

###################################################################
#Gage and list construction - Meristem
###################################################################

library(gage)

#Exclude lowly expressed genes for GSEA
DESeq2_negative_gene_IDs <- is.na(as.data.frame(ME_TN_vs_WN_DEGs$log2FoldChange))

###################################################################

list <- list()
for(i in 1:58){
  
  TF_class <- as.character(unique(TFs$superfamily))
  TF_class_name <- TF_class[i]
  
  list[[i]] <- TFs[grep(paste(TF_class_name),TFs$superfamily),2]
}
names(list)<-TF_class[1:58]

#Run GAGE command for all leaky and induced expressed transgenics

Enriched <- gage(counts(dds)[!DESeq2_negative_gene_IDs,],list,ref=c(6:7),samp=c(8:10),
                 rank.test = T,same.dir = T,
                 set.size=c(1,800), compare="unpaired")

Enriched_greater <- Enriched$greater[1:18,1:5]
Enriched_lesser <- Enriched$less[1:14,1:5]

Enriched_write <- rbind(Enriched_greater,Enriched_lesser)
write.table(Enriched$stats,file = '../../analysis/RKD_TF_gene_enrichment.txt', quote =F, sep= "\t")
write.table(Enriched$greater,file = '../../analysis/TF_greater_gene_enrichment.txt', quote =F, sep= "\t")

q.val <- -log10(Enriched_write[,4])
q.val[16:32] <- q.val[16:32]*-1

data<-data.frame(rownames(Enriched_write), q.val)
colnames(data) <- c("TF","q.val")

colours <- c(rep("indianred1",15), rep("royalblue",17))

library("RColorBrewer")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, 8))

library(wesanderson)

ggplot(data, aes(TF, q.val, fill=q.val)) + geom_bar(stat="identity") +
  scale_fill_continuous(low="yellow", high="red") +
  coord_flip()

ggplot(userData, aes(month, count, fill = count)) +
  geom_bar(stat = "identity") +
  scale_x_date() + 
  scale_fill_continuous(low="blue", high="red") +
  labs(x= "Time", y="Count")


###################################################################
#Gage and list construction - Immature embryo
###################################################################

library(gage)

#Exclude lowly expressed genes for GSEA
DESeq2_negative_gene_IDs <- is.na(as.data.frame(IE_TE_vs_WE_DEG_ID$log2FoldChange))

###################################################################

list <- list()
for(i in 1:58){
  
  TF_class <- as.character(unique(TFs$superfamily))
  TF_class_name <- TF_class[i]
  
  list[[i]] <- TFs[grep(paste(TF_class_name),TFs$superfamily),2]
}
names(list)<-TF_class[1:58]

#Run GAGE command for all leaky and induced expressed transgenics

Enriched <- gage(counts(dds)[!DESeq2_negative_gene_IDs,],list,ref=c(4:5),samp=c(1:3),
                 rank.test = T,same.dir = T,
                 set.size=c(1,800), compare="unpaired")

Enriched_greater <- Enriched$greater[1:18,1:5]
Enriched_lesser <- Enriched$less[1:14,1:5]

Enriched_write <- rbind(Enriched_greater,Enriched_lesser)
write.table(Enriched$stats,file = 'RKD_TF_gene_enrichment.txt', quote =F, sep= "\t")
write.table(Enriched$greater,file = 'TF_greater_gene_enrichment.txt', quote =F, sep= "\t")

q.val <- -log10(Enriched_write[,4])
q.val[16:32] <- q.val[16:32]*-1

data<-data.frame(rownames(Enriched_write), q.val)
colnames(data) <- c("TF","q.val")

colours <- c(rep("indianred1",15), rep("royalblue",17))

library("RColorBrewer")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, 8))

library(wesanderson)

ggplot(data, aes(TF, q.val, fill=q.val)) + geom_bar(stat="identity") +
  scale_fill_continuous(low="yellow", high="red") +
  coord_flip()

ggplot(userData, aes(month, count, fill = count)) +
  geom_bar(stat = "identity") +
  scale_x_date() + 
  scale_fill_continuous(low="blue", high="red") +
  labs(x= "Time", y="Count")

