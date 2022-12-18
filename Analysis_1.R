library(IBRAP)

objects <- readRDS('individual_samples.rds')

for(x in names(objects)[1:8]) { 
  
  objects[[x]] <- perform.sct(objects[[x]])
  objects[[x]] <- perfrom.scran(objects[[x]], vars.to.regress = 'RAW_total.counts')
  objects[[x]] <- perform.scanpy(objects[[x]], vars.to.regress = 'RAW_total.counts')
  
  objects[[x]] <- perform.pca(object = objects[[x]])
  
  objects[[x]] <- perform.umap(object = objects[[x]], assay = c('SCT','SCRAN','SCANPY'), reduction = 'PCA', dims.use = list(20))
  
  objects[[x]] <- perform.nn(object = objects[[x]], assay = c('SCT','SCRAN','SCANPY'), reduction = 'PCA', dims.use = list(20))
  
  # louvain clustering
  
  objects[[x]] <- perform.graph.cluster(object = objects[[x]], assay = c('SCT','SCRAN','SCANPY'), neighbours = 'PCA_NN', algorithm = 1)
  
  # louvain with multilevel refinement clustering
  
  objects[[x]] <- perform.graph.cluster(object = objects[[x]], assay = c('SCT','SCRAN','SCANPY'), neighbours = 'PCA_NN', algorithm = 2)
  
  # SLM clustering
  
  objects[[x]] <- perform.graph.cluster(object = objects[[x]], assay = c('SCT','SCRAN','SCANPY'), neighbours = 'PCA_NN', algorithm = 3)
  
  # leiden clustering
  
  objects[[x]] <- perform.graph.cluster(object = objects[[x]], assay = c('SCT','SCRAN','SCANPY'), neighbours = 'PCA_NN', algorithm = 4)
  
  # PAM clustering
  
  objects[[x]] <- perform.reduction.kmeans(object = objects[[x]], assay = c('SCT','SCRAN','SCANPY'), reduction = 'PCA_UMAP', 
                                           dims = list(2), k = 10:15, method = 'kmeans')
  
  # kmeans clustering 
  
  objects[[x]] <- perform.reduction.kmeans(object = objects[[x]], assay = c('SCT','SCRAN','SCANPY'), reduction = 'PCA_UMAP', 
                                           dims = list(2), k = 10:15, method = 'pam')
  
  objects[[x]] <- perform.sc3.slot.cluster(object = objects[[x]], assay = c('SCT','SCRAN','SCANPY'), slot = 'normalised', HVGs = T,  
                                           dims = list(2), ks = 10:15, n.core = 3)
  
  objects[[x]] <- benchmark.clustering(object = objects[[x]], assay = c('SCT','SCRAN','SCANPY'), 
                                           clustering = c('PCA_NN:LOUVAIN', 'PCA_NN:LOUVAINMLR', 'PCA_NN:SLM', 
                                           'PCA_NN:LEIDEN', 'PCA_UMAP:kmeans', 'PCA_UMAP:pam', 'normalised:SC3'), 
                                           reduction = c('PCA_UMAP', 'PCA_UMAP', 'PCA_UMAP', 'PCA_UMAP', 
                                                         'PCA_UMAP', 'PCA_UMAP', 'PCA_UMAP'), n.dims = 2, 
                                           ground.truth.column = 'celltype', verbose = T, threads = 6)
  
}


