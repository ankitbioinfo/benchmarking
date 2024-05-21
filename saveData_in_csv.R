
Sys.setenv(LANG="en")


count <- readRDS("counts.Rds")
B1 = as.data.frame(as.matrix(count))
dim(B1)
E1=select(B1,contains("embryo1"))
dim(E1)
E2=select(B1,contains("embryo2"))
dim(E2)
E3=select(B1,contains("embryo3"))
dim(E3)

[1]   351 57536
[1]   351 19451
[1]   351 14891
[1]   351 23194

write.csv(E1, "embryo1/gene_by_cell.csv",quote=FALSE,row.names=TRUE)
write.csv(E2, "embryo2/gene_by_cell.csv",quote=FALSE,row.names=TRUE)
write.csv(E3, "embryo3/gene_by_cell.csv",quote=FALSE,row.names=TRUE)




#ct <- readRDS("metadata.Rds")


lapply(myList, function(x) { x["segmentation_vertices_x_global_affine"] <- NULL; x })
segmentation_vertices_y_global_affine


write.csv(d['uniqueID'], "metadata_1.csv",quote=FALSE,row.names=TRUE)
write.csv(d['embryo'], "metadata_2.csv",quote=FALSE,row.names=TRUE)
write.csv(d['pos'], "metadata_3.csv",quote=FALSE,row.names=TRUE)
write.csv(d['z'], "metadata_4.csv",quote=FALSE,row.names=TRUE)
write.csv(d['x_global'], "metadata_5.csv",quote=FALSE,row.names=TRUE)
write.csv(d['y_global'], "metadata_6.csv",quote=FALSE,row.names=TRUE)
write.csv(d['x_global_affine'], "metadata_7.csv",quote=FALSE,row.names=TRUE)
write.csv(d['y_global_affine'], "metadata_8.csv",quote=FALSE,row.names=TRUE)
write.csv(d['embryo_pos'], "metadata_9.csv",quote=FALSE,row.names=TRUE)
write.csv(d['Area'] , "metadata_10.csv",quote=FALSE,row.names=TRUE)
write.csv(d['UMAP1'] , "metadata_11.csv",quote=FALSE,row.names=TRUE)
write.csv(d['UMAP2'] , "metadata_12.csv",quote=FALSE,row.names=TRUE)
write.csv(d['celltype_mapped_refined'] , "metadata_13.csv",quote=FALSE,row.names=TRUE)










#cat("nrow = B1 ", nrow(ct), "ncol = ", ncol(ct), "\n");

#write.csv(ct,"metadata_info.csv",quote=FALSE,row.names=TRUE)

#ct <- readRDS("mRNA.Rds")
#write.csv(ct,"gene_by_cell.csv",quote=FALSE,row.names=TRUE)
