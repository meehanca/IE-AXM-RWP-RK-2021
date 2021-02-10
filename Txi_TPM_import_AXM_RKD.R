
##### #1: Summarise counts per gene ########## 

#source("https://bioconductor.org/biocLite.R")
#biocLite("tximportData")
#install.packages("readr")
library(tximportData)
library(readr)
library(rhdf5)

# read in pre-constructed tx2gene table (transcript to gene table)
tx2gene <- read.table("../philippa_resources/transcripts_to_genes_RefSeqv1.0_annot_v1.1.txt", header=T)
head(tx2gene)

files <- c("YL100/abundance.tsv","YL101/abundance.tsv","YL103/abundance.tsv","YL104/abundance.tsv","YL105/abundance.tsv","YL107/abundance.tsv","YL117/abundance.tsv","YL118/abundance.tsv","YL120/abundance.tsv",
           "YL126/abundance.tsv","YL127/abundance.tsv","YL129/abundance.tsv","YL130/abundance.tsv","YL131/abundance.tsv","YL133/abundance.tsv","YL136/abundance.tsv","YL137/abundance.tsv","YL138/abundance.tsv",
           "YL141/abundance.tsv","YL144/abundance.tsv","YL145/abundance.tsv")
names(files) <- c("IE_TE1","IE_TE2","IE_TE4","IE_TM1","IE_TM2","IE_TM4","IE_WE1","IE_WE2","IE_WE3",
                  "Me_WN1","Me_WN2","Me_WN4","Me_TN1","Me_TN2","Me_TN4","Me_TM1","Me_TM2","Me_TM3",
                  "Me_TE2","Me_TE5","Me_TE6")

samples <- read.table(file.path(htseqDir,"sampleTable.txt"), header=T)

library(tximport)
# read in the files and sum per gene
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
names(txi)

samples <- read.table(file.path(htseqDir,"samples.txt"), header=T)


for(i in samples){
  # save counts summarised per gene
  write.table(txi$counts[,i], file=paste(i,"_count.tsv",sep=''),sep = "\t",quote=F,col.names = F)
}
# to see tpm summarised per gene
#head(txi$abundance)
#colnames(txi$abundance)

# save tpm summarised per gene
write.table(txi$abundance, file="_tpm.tsv",sep = "\t")

# see lengths summarised per gene
head(txi$length)

# calculate average gene length across all samples
gene_lengths <- as.data.frame(rowMeans(txi$length))
head(gene_lengths)
colnames(gene_lengths) <- c("length")
head(gene_lengths)
#save length per gene
write.csv(gene_lengths, file="_gene_lengths.csv")


