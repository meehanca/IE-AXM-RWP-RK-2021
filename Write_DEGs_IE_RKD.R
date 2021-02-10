###################################################################
#Preparing lists of diffrentially expressed genes for online tools
####################################################################

################## All DEGs ##################

# Extracting Gene IDsA

IE_TE_vs_WE_DEG_ID <- rownames(IE_TE_vs_WE_DEG)
IE_TE_vs_TM_DEG_ID <- rownames(IE_TE_vs_TM_DEG)
IE_TM_vs_WE_DEG_ID <- rownames(IE_TM_vs_WE_DEG)

# Writing tables - All DEG Counts

write.table(IE_TE_vs_WE_DEGs, file = "./DEGs/Immature Embryo/DEGs/IE_TE_vs_WE_DEGs_counts.txt", quote = F, row.names = T, col.names = T, sep='\t')
write.table(IE_TE_vs_TM_DEGs, file = "./DEGs/Immature Embryo/DEGs/IE_TE_vs_TM_DEGs_counts.txt", quote = F, row.names = T, col.names = T, sep='\t')
write.table(IE_TM_vs_WE_DEGs, file = "./DEGs/Immature Embryo/DEGs/IE_TM_vs_WE_DEGs_counts.txt", quote = F, row.names = T, col.names = T, sep='\t')

# Writing tables - All DEG gene IDs

write.table(IE_TE_vs_WE_DEG_ID, file = "./DEGS/Immature Embryo/DEGs/IE_TE_vs_WE_DEG_ID.txt", quote = F, row.names = F, col.names = F)
write.table(IE_TE_vs_TM_DEG_ID, file = "./DEGs/Immature Embryo/DEGs/IE_TE_vs_TM_DEG_ID.txt", quote = F, row.names = F, col.names = F)
write.table(IE_TM_vs_WE_DEG_ID, file = "./DEGs/Immature Embryo/DEGs/IE_TM_vs_WE_DEG_ID.txt", quote = F, row.names = F, col.names = F)

################## Upregulated DEGs ##################

# Extracting upregulated genes

IE_TE_vs_WE_upregulated_DEG<-IE_TE_vs_WE_DEG[IE_TE_vs_WE_DEG$log2FoldChange>0,]
IE_TE_vs_TM_upregulated_DEG<-IE_TE_vs_TM_DEG[IE_TE_vs_TM_DEG$log2FoldChange>0,]
IE_TM_vs_WE_upregulated_DEG<-IE_TM_vs_WE_DEG[IE_TM_vs_WE_DEG$log2FoldChange>0,]

# Writing tables - Upregulated DEG counts

write.table(IE_TE_vs_WE_upregulated_DEG, file = "./DEGs/Immature Embryo/Upregulated/TE_vs_NE_upregulated_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(IE_TE_vs_TM_upregulated_DEG, file = "./DEGs/Immature Embryo/Upregulated/TE_vs_TM_upregulated_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(IE_TM_vs_WE_upregulated_DEG, file = "./DEGs/Immature Embryo/Upregulated/TM_vs_NE_upregulated_DEGs.txt", quote = F, row.names = T, col.names = T)

##  Gene IDs

IE_TE_vs_WE_upregulated_DEG_ID <- rownames(IE_TE_vs_WE_upregulated_DEG)
IE_TE_vs_TM_upregulated_DEG_ID <- rownames(IE_TE_vs_TM_upregulated_DEG)
IE_TM_vs_WE_upregulated_DEG_ID <- rownames(IE_TM_vs_WE_upregulated_DEG)

# Writing tables - Upregulated DEG IDs

write.table(IE_TE_vs_WE_upregulated_DEG_ID , file = "./DEGs/Immature Embryo/Upregulated/IE_TE_vs_WE_upregulated_DEG_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(IE_TE_vs_TM_upregulated_DEG_ID, file = "./DEGs/Immature Embryo/Upregulated/IE_TE_vs_TM_upregulated_DEG_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(IE_TM_vs_WE_upregulated_DEG_ID, file = "./DEGs/Immature Embryo/Upregulated/IE_TM_vs_WE_upregulated_DEG_IDs_txt", quote = F, row.names = F, col.names = F)

################## Downregulated DEGs ##################

# Extracting downregulated genes

IE_TE_vs_WE_downregulated_DEG<-IE_TE_vs_WE_DEG[IE_TE_vs_WE_DEG$log2FoldChange<0,]
IE_TE_vs_TM_downregulated_DEG<-IE_TE_vs_TM_DEG[IE_TE_vs_TM_DEG$log2FoldChange<0,]
IE_TM_vs_WE_downregulated_DEG<-IE_TM_vs_WE_DEG[IE_TM_vs_WE_DEG$log2FoldChange<0,]

# Writing tables - Downregulated DEG counts

write.table(IE_TE_vs_WE_downregulated_DEG, file = "./DEGs/Immature Embryo/Downregulated/IE_TE_vs_WE_downregulated_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(IE_TE_vs_TM_downregulated_DEG, file = "./DEGs/Immature Embryo/Downregulated/IE_TE_vs_TM_downregulated_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(IE_TM_vs_WE_downregulated_DEG, file ="./DEGs/Immature Embryo/Downregulated/IE_TM_vs_WE_downregulated_DEGs.txt", quote = F, row.names = T, col.names = T)

##  Gene IDs

IE_TE_vs_WE_downregulated_DEG_ID <- rownames(IE_TE_vs_WE_downregulated_DEG)
IE_TE_vs_TM_downregulated_DEG_ID <- rownames(IE_TE_vs_TM_downregulated_DEG)
IE_TM_vs_WE_downregulated_DEG_ID <- rownames(IE_TM_vs_WE_downregulated_DEG)

# Writing tables - Downregulated DEG counts

write.table(IE_TE_vs_WE_downregulated_DEG_ID, file = "./DEGs/Immature Embryo/Downregulated/IE_TE_vs_WE_downregulated_DEG_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(IE_TE_vs_TM_downregulated_DEG_ID, file = "./DEGs/Immature Embryo/Downregulated/IE_TE_vs_TM_downregulated_DEG_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(IE_TM_vs_WE_downregulated_DEG_ID, file = "./DEGs/Immature Embryo/Downregulated/IE_TM_vs_WE_downregulated_DEG_IDs_txt", quote = F, row.names = F, col.names = F)

################## Shared DEGs between (TE vs NE) and (TM vs NE) ##################

## Find shared DEGs

Shared_IE_DEGs_IDs<- IE_TE_vs_WE_DEG_ID[IE_TE_vs_WE_DEG_ID %in% IE_TM_vs_WE_DEG_ID]
Shared_IE_upregulated_DEGs_IDs<- IE_TE_vs_WE_upregulated_DEG_ID[IE_TE_vs_WE_upregulated_DEG_ID %in% IE_TM_vs_WE_upregulated_DEG_ID]
Shared_IE_downregulated_DEGs_IDs<- IE_TE_vs_WE_downregulated_DEG_ID[IE_TE_vs_WE_downregulated_DEG_ID %in% IE_TM_vs_WE_downregulated_DEG_ID]

# Writing tables - Shared DEG IDs (Can't do p-values as they will be different for both calculations so Shared ID is the best we can do...)

write.table(Shared_IE_DEGs_IDs, file = "./DEGs/Immature Embryo/Shared/Shared_DEGs_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(Shared_IE_upregulated_DEGs_IDs, file = "./DEGs/Immature Embryo/Shared/Shared_upregulated_DEGs_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(Shared_IE_downregulated_DEGs_IDs, file = "./DEGs/Immature Embryo/Shared/Shared_downregulated_DEGs_IDs.txt", quote = F, row.names = F, col.names = F)

