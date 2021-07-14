library(Seurat)
library(stringr)

preprocess_PBMC <- function(filepath) {
  
  print(paste('Preprocessing:', filepath, sep=' '))
  
  out.data <- Read10X(data.dir = filepath)
  out <- CreateSeuratObject(counts = out.data)
  rm(out.data)
  out[["percent.mt"]] <- PercentageFeatureSet(out, pattern = "^mt-|^MT-")
  expr1 <- FetchData(object = out, vars = "nCount_RNA")
  expr2 <- FetchData(object = out, vars = "percent.mt")
  out <- out[, which(x = expr1 >= 500 & expr2 <= 15)]
  
  out <- NormalizeData(object = out, normalization.method = "LogNormalize", scale.factor = 10000)
  out <- FindVariableFeatures(out, selection.method = "vst", nfeatures = 2000)
  out <- ScaleData(object = out)
  out <- RunPCA(out)
  out <- RunUMAP(out, dims = 1:30)
  
  out <- AddModuleScore(object = out, features = list(c('CD4', 'IL7R', 'CD3D', 'CTLA4')), name = 'T_CD4_markers', ctrl=5)
  out <- AddModuleScore(object = out, features = list(c('CD8A', 'CD8B', 'GZMB', 'CD3D')), name = 'T_CD8_markers', ctrl=5)
  out <- AddModuleScore(object = out, features = list(c('CD19', 'MS4A1', 'CD79A', 'CD79B', 'BLNK')), name = 'B_markers', ctrl=5)
  out <- AddModuleScore(object = out, features = list(c('FCGR3A', 'NCAM1', 'KLRB1', 'KLRC1', 'KLRD1', 'KLRF1', 'GNLY', 'NKG7')), name = 'NK_markers', ctrl=5)
  out <- AddModuleScore(object = out, features = list(c('CD14', 'LYZ')), name = 'Mono_CD14_markers', ctrl=5)
  out <- AddModuleScore(object = out, features = list(c('FCGR3A', 'MS4A7')), name = 'Mono_CD16_markers', ctrl=5)
  out <- AddModuleScore(object = out, features = list(c('IL3RA', 'CLEC4C', 'NRP1', 'FCER1A', 'CST3')), name = 'DC_markers', ctrl=5)

  saveRDS(out, file=paste(filepath, 'Processed.rds', sep='/'))

  write.table(t(out@assays$RNA@scale.data), 'NormalizedCounts.csv', col.names = F, sep=",", row.names = F)
  write.table(out@assays$RNA@var.features, 'PCA_Genes.csv', col.names = F, sep=",", row.names = F, quote=FALSE)
  
  write.table(out@reductions$pca@cell.embeddings, paste(filepath, 'PCA.csv', sep='/'), col.names = F, sep=",")
  write.table(out@reductions$umap@cell.embeddings, paste(filepath, 'UMap.csv', sep='/'), col.names = F, sep=",")
  write.table(out$nCount_RNA, paste(filepath, 'TotCounts.csv', sep='/'), col.names = F, sep=",")
  
  score_mat <- cbind(out$T_CD4_markers1, out$T_CD8_markers1, out$B_markers1, out$NK_markers1, out$Mono_CD14_markers1, out$Mono_CD16_markers1, out$DC_markers1)
  write.table(score_mat, 'CellTypeScore.csv', colnames<-c('CD4+ T', 'CD8+ T', 'B', 'NK', 'CD14+ Mono', 'CD16+ Mono', 'DC'), sep=",")
  
  return(out)
  
}
