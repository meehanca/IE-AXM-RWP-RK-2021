###################################################################
#Preparing lists of diffrentially expressed genes for online tools
####################################################################

################## All DEGs ##################

# Extracting Gene IDs

ME_TN_vs_WN_DEG_ID <- rownames(ME_TN_vs_WN_DEG)
ME_TE_vs_TM_DEG_ID <- rownames(ME_TE_vs_TM_DEG)
ME_TE_vs_WN_DEG_ID <- rownames(ME_TE_vs_WN_DEG)
ME_TM_vs_WN_DEG_ID <- rownames(ME_TM_vs_WN_DEG)

# Writing tables - All DEG Counts

write.table(ME_TN_vs_WN_DEGs, file = "./DEGs/Meristem/DEGs/ME_TN_vs_WN_DEGs_counts.txt", quote = F, row.names = T, col.names = T, sep='\t')
write.table(ME_TE_vs_TM_DEGs, file = "./DEGs/Meristem/DEGs/ME_TE_vs_TM_DEGs_counts.txt", quote = F, row.names = T, col.names = T, sep='\t')
write.table(ME_TE_vs_WN_DEGs, file = "./DEGs/Meristem/DEGs/ME_TE_vs_WN_DEGs_counts.txt", quote = F, row.names = T, col.names = T, sep='\t')
write.table(ME_TM_vs_WN_DEGs, file = "./DEGs/Meristem/DEGs/ME_TM_vs_WN_DEGs_counts.txt", quote = F, row.names = T, col.names = T, sep='\t')

# Writing tables - All DEG gene IDs

write.table(ME_TN_vs_WN_DEG_ID, file = "./DEGS/Meristem/DEGs/ME_TN_vs_WN_DEG_ID.txt", quote = F, row.names = F, col.names = F)
write.table(ME_TE_vs_TM_DEG_ID, file = "./DEGs/Meristem/DEGs/ME_TE_vs_TM_DEG_ID.txt", quote = F, row.names = F, col.names = F)
write.table(ME_TE_vs_WN_DEG_ID, file = "./DEGS/Meristem/DEGs/ME_TE_vs_WN_DEG_ID.txt", quote = F, row.names = F, col.names = F)
write.table(ME_TM_vs_WN_DEG_ID, file = "./DEGs/Meristem/DEGs/ME_TM_vs_WN_DEG_ID.txt", quote = F, row.names = F, col.names = F)

################## Upregulated DEGs ##################

# Extracting upregulated genes

ME_TN_vs_WN_upregulated_DEG<-ME_TN_vs_WN_DEG[ME_TN_vs_WN_DEG$log2FoldChange>0,]
ME_TE_vs_TM_upregulated_DEG<-ME_TE_vs_TM_DEG[ME_TE_vs_TM_DEG$log2FoldChange>0,]
ME_TE_vs_WN_upregulated_DEG<-ME_TE_vs_WN_DEG[ME_TE_vs_WN_DEG$log2FoldChange>0,]
ME_TM_vs_WN_upregulated_DEG<-ME_TM_vs_WN_DEG[ME_TM_vs_WN_DEG$log2FoldChange>0,]

# Writing tables - Upregulated DEG counts

write.table(ME_TN_vs_WN_upregulated_DEG, file = "./DEGs/Meristem/Upregulated/TE_vs_NE_upregulated_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(ME_TE_vs_TM_upregulated_DEG, file = "./DEGs/Meristem/Upregulated/TE_vs_TM_upregulated_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(ME_TE_vs_WN_upregulated_DEG, file = "./DEGs/Meristem/Upregulated/TE_vs_WN_upregulated_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(ME_TM_vs_WN_upregulated_DEG, file = "./DEGs/Meristem/Upregulated/TM_vs_WN_upregulated_DEGs.txt", quote = F, row.names = T, col.names = T)

##  Gene IDs

ME_TN_vs_WN_upregulated_DEG_ID <- rownames(ME_TN_vs_WN_upregulated_DEG)
ME_TE_vs_TM_upregulated_DEG_ID <- rownames(ME_TE_vs_TM_upregulated_DEG)
ME_TE_vs_WN_upregulated_DEG_ID <- rownames(ME_TE_vs_WN_upregulated_DEG)
ME_TM_vs_WN_upregulated_DEG_ID <- rownames(ME_TM_vs_WN_upregulated_DEG)

# Writing tables - Upregulated DEG IDs

write.table(ME_TN_vs_WN_upregulated_DEG_ID , file = "./DEGs/Meristem/Upregulated/ME_TN_vs_WN_upregulated_DEG_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(ME_TE_vs_TM_upregulated_DEG_ID, file = "./DEGs/Meristem/Upregulated/ME_TE_vs_TM_upregulated_DEG_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(ME_TE_vs_WN_upregulated_DEG_ID , file = "./DEGs/Meristem/Upregulated/ME_TE_vs_WN_upregulated_DEG_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(ME_TM_vs_WN_upregulated_DEG_ID, file = "./DEGs/Meristem/Upregulated/ME_TM_vs_WN_upregulated_DEG_IDs.txt", quote = F, row.names = F, col.names = F)

################## Downregulated DEGs ##################

# Extracting downregulated genes

ME_TN_vs_WN_downregulated_DEG<-ME_TN_vs_WN_DEG[ME_TN_vs_WN_DEG$log2FoldChange<0,]
ME_TE_vs_TM_downregulated_DEG<-ME_TE_vs_TM_DEG[ME_TE_vs_TM_DEG$log2FoldChange<0,]
ME_TE_vs_WN_downregulated_DEG<-ME_TE_vs_WN_DEG[ME_TE_vs_WN_DEG$log2FoldChange<0,]
ME_TM_vs_WN_downregulated_DEG<-ME_TM_vs_WN_DEG[ME_TM_vs_WN_DEG$log2FoldChange<0,]

# Writing tables - Downregulated DEG counts

write.table(ME_TN_vs_WN_downregulated_DEG, file = "./DEGs/Meristem/Downregulated/ME_TN_vs_WN_downregulated_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(ME_TE_vs_TM_downregulated_DEG, file = "./DEGs/Meristem/Downregulated/ME_TE_vs_TM_downregulated_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(ME_TE_vs_WN_downregulated_DEG, file = "./DEGs/Meristem/Downregulated/ME_TN_vs_WN_downregulated_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(ME_TM_vs_WN_downregulated_DEG, file = "./DEGs/Meristem/Downregulated/ME_TE_vs_TM_downregulated_DEGs.txt", quote = F, row.names = T, col.names = T)

##  Gene IDs

ME_TN_vs_WN_downregulated_DEG_ID <- rownames(ME_TN_vs_WN_downregulated_DEG)
ME_TE_vs_TM_downregulated_DEG_ID <- rownames(ME_TE_vs_TM_downregulated_DEG)
ME_TE_vs_WN_downregulated_DEG_ID <- rownames(ME_TE_vs_WN_downregulated_DEG)
ME_TM_vs_WN_downregulated_DEG_ID <- rownames(ME_TM_vs_WN_downregulated_DEG)

# Writing tables - Downregulated DEG counts

write.table(ME_TN_vs_WN_downregulated_DEG_ID, file = "./DEGs/Meristem/Downregulated/ME_TN_vs_WN_downregulated_DEG_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(ME_TE_vs_TM_downregulated_DEG_ID, file = "./DEGs/Meristem/Downregulated/ME_TE_vs_TM_downregulated_DEG_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(ME_TN_vs_WN_downregulated_DEG_ID, file = "./DEGs/Meristem/Downregulated/ME_TN_vs_WN_downregulated_DEG_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(ME_TE_vs_TM_downregulated_DEG_ID, file = "./DEGs/Meristem/Downregulated/ME_TE_vs_TM_downregulated_DEG_IDs.txt", quote = F, row.names = F, col.names = F)

################## Shared DEGs betWNen (TE vs NE) and (TM vs NE) ##################

## Find shared DEGs

Shared_DEGs_IDs<- ME_TN_vs_WN_DEG_ID[ME_TN_vs_WN_DEG_ID %in% ME_TM_vs_WN_DEG_ID %in% ME_TE_vs_WN_DEG_ID]
Shared_upregulated_DEGs_IDs<- ME_TN_vs_WN_upregulated_DEG_ID[ME_TN_vs_WN_upregulated_DEG_ID %in% ME_TE_vs_WN_upregulated_DEG_ID %in% ME_TM_vs_WN_upregulated_DEG_ID]
Shared_downregulated_DEGs_IDs<- ME_TN_vs_WN_downregulated_DEG_ID[ME_TN_vs_WN_downregulated_DEG_ID %in% ME_TE_vs_TM_downregulated_DEG_ID %in% ME_TM_vs_WN_downregulated_DEG_ID]

# Writing tables - Shared DEG IDs (Can't do p-values as they will be different for both calculations so Shared ID is the best WN can do...)

write.table(Shared_DEGs_IDs, file = "./DEGs/Meristem/Shared/Shared_DEGs_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(Shared_upregulated_DEGs_IDs, file = "./DEGs/Meristem/Shared/Shared_upregulated_DEGs_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(Shared_downregulated_DEGs_IDs, file = "./DEGs/Meristem/Shared/Shared_downregulated_DEGs_IDs.txt", quote = F, row.names = F, col.names = F)

