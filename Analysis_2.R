library(IBRAP)
library(Matrix)

object <- readRDS('integrated_samples.rds')

for(x in names(object)) {
  
  object[[x]] <- perform.sct(object = object[[x]])
  object[[x]] <- perform.scran(object = object[[x]], vars.to.regress = 'RAW_total.counts')
  object[[x]] <- perform.scanpy(object = object[[x]], vars.to.regress = 'RAW_total.counts')
  
  object[[x]] <- perform.pca(object = object[[x]], assay = c('SCT','SCRAN','SCANPY'))
  
  object[[x]] <- perform.scanorama(object = object[[x]], assay = c('SCT','SCRAN','SCANPY'), batch = 'batch', n.dims = 20)
  
  object[[x]] <- perform.bbknn(object = object[[x]], assay = c('SCT','SCRAN','SCANPY'), reduction = 'PCA', batch = 'batch', n_pcs = list(20), verbose = T)
  
  object[[x]] <- perform.harmony(object = object[[x]], assay = c('SCT','SCRAN','SCANPY'), reduction = 'PCA', batch = 'batch', dims.use = list(20))
  
  object[[x]]_batches <- splitIBRAP(object = object[[x]], split.by = 'batch')
  
  object[[x]]_batches <- lapply(X = object[[x]]_batches, FUN = perform.sct)
  object[[x]]_batches <- lapply(X = object[[x]]_batches, FUN = perform.scran, vars.to.regress = 'RAW_total.counts')
  object[[x]]_batches <- lapply(X = object[[x]]_batches, FUN = perform.scanpy, vars.to.regress = 'RAW_total.counts')
  
  object[[x]] <- perform.seurat.integration(object = object[[x]], object.list = object[[x]]_batches, assay = c('SCT','SCRAN','SCANPY'), verbose = T)
  
  object[[x]] <- perform.nn(object = object[[x]], assay = c('SCT','SCRAN','SCANPY'), verbose = T, reduction = c('PCA_HARMONY', 'SCANORAMA', 'CCA', 'PCA'), dims.use = list(0,0,20,20))
  
  object[[x]] <- perform.graph.cluster(object = object[[x]], assay = c('SCT','SCRAN','SCANPY'), verbose = T, neighbours = c('PCA_HARMONY_NN', 'SCANORAMA_NN', 'CCA_NN', 'PCA_BBKNN_BBKNN', 'PCA_NN'))
  
  object[[x]] <- perform.umap(object = object[[x]], assay = c('SCT','SCRAN','SCANPY'), verbose = T, reduction = c('PCA_HARMONY', 'SCANORAMA', 'CCA', 'PCA'))
  
  object[[x]] <- perform.umap(object = object[[x]], assay = c('SCT','SCRAN','SCANPY'), verbose = T, graph = c('PCA_BBKNN_BBKNN'))
  
  object[[x]] <- perform.reduction.kmeans(object = object[[x]], assay = c('SCT','SCRAN','SCANPY'), reduction = c('PCA_HARMONY_UMAP', 'SCANORAMA_UMAP', 'CCA_UMAP', 'PCA_BBKNN_BBKNN:UMAP', 'PCA_UMAP'), dims = list(2,2,2,2,2), k = 2:6)
  
  object[[x]] <- perform.reduction.kmeans(object = object[[x]], method = 'pam', assay = c('SCT','SCRAN','SCANPY'), reduction = c('PCA_HARMONY_UMAP', 'SCANORAMA_UMAP', 'CCA_UMAP', 'PCA_BBKNN_BBKNN:UMAP', 'PCA_UMAP'), dims = list(2,2,2,2,2), k = 2:6)
  
  object[[x]] <- perform.sc3.reduction.cluster(object = object[[x]], assay = c('SCT','SCRAN','SCANPY'), reduction = c('PCA_HARMONY', 'SCANORAMA', 'CCA', 'PCA'), dims = list(20,20,20,20), ks = 2:6, n.core = 8)
  
  object[[x]] <- perform.sc3.slot.cluster(object = object[[x]], assay = c('SCT','SCRAN','SCANPY'), slot = 'normalised', HVGs = T, ks = 2:6, n.core = 8)
  
  object[[x]] <- benchmark.integration(object = object[[x]], batch = 'batch', assays = c('SCRAN','SCANPY','SCT'), 
                                       reduction = c('PCA_UMAP', 'PCA_HARMONY_UMAP', 'SCANORAMA_UMAP',
                                                     'CCA_UMAP', 'PCA_BBKNN_BBKNN:UMAP'), 
                                       result.names = c('uncorrected', 'harmony','scanorama','cca', 'bbknn'), threads = 2)
  
  object[[x]] <- benchmark.clustering(object = object[[x]], assay = c('SCT','SCRAN','SCANPY'), 
                                      clustering = c('PCA_HARMONY_NN:LOUVAIN', 'SCANORAMA_NN:LOUVAIN', 'CCA_NN:LOUVAIN', 
                                                     'PCA_BBKNN_BBKNN:LOUVAIN', 'PCA_HARMONY_UMAP:kmeans', 'SCANORAMA_UMAP:kmeans', 
                                                     'CCA_UMAP:kmeans', 'PCA_BBKNN_BBKNN:UMAP:kmeans', 'PCA_NN:LOUVAIN', 
                                                     'PCA_UMAP:kmeans', 'PCA_HARMONY:SC3', 'SCANORAMA:SC3', 'CCA:SC3', 
                                                     'PCA:SC3', 'normalised:SC3', 'PCA_HARMONY_UMAP:pam', 
                                                     'SCANORAMA_UMAP:pam', 'CCA_UMAP:pam', 'PCA_UMAP:pam', 
                                                     'PCA_BBKNN_BBKNN:UMAP:pam', 'PCA_HARMONY_NN:LOUVAINMLR', 'SCANORAMA_NN:LOUVAINMLR', 
                                                     'CCA_NN:LOUVAINMLR', 'PCA_NN:LOUVAINMLR', 'PCA_BBKNN_BBKNN:LOUVAINMLR', 
                                                     'PCA_HARMONY_NN:SLM', 'SCANORAMA_NN:SLM', 'CCA_NN:SLM', 
                                                     'PCA_NN:SLM', 'PCA_BBKNN_BBKNN:SLM', 'PCA_HARMONY_NN:LEIDEN', 
                                                     'SCANORAMA_NN:LEIDEN', 'CCA_NN:LEIDEN','PCA_NN:LEIDEN', 
                                                     'PCA_BBKNN_BBKNN:LEIDEN'), 
                                      reduction = c('PCA_HARMONY_UMAP', 'SCANORAMA_UMAP', 'CCA_UMAP',
                                                    'PCA_BBKNN_BBKNN:UMAP', 'PCA_HARMONY_UMAP', 'SCANORAMA_UMAP',
                                                    'CCA_UMAP', 'PCA_BBKNN_BBKNN:UMAP', 'PCA_UMAP',
                                                    'PCA_UMAP', 'PCA_HARMONY_UMAP', 'SCANORAMA_UMAP', 'CCA_UMAP',
                                                    'PCA_UMAP', 'PCA_UMAP', 'PCA_HARMONY_UMAP', 
                                                    'SCANORAMA_UMAP', 'CCA_UMAP', 'PCA_UMAP',
                                                    'PCA_BBKNN_BBKNN:UMAP', 'PCA_HARMONY_UMAP', 'SCANORAMA_UMAP',
                                                    'CCA_UMAP', 'PCA_UMAP', 'PCA_BBKNN_BBKNN:UMAP', 
                                                    'PCA_HARMONY_UMAP', 'SCANORAMA_UMAP', 'CCA_UMAP',
                                                    'PCA_UMAP', 'PCA_BBKNN_BBKNN:UMAP', 'PCA_HARMONY_UMAP',
                                                    'SCANORAMA_UMAP', 'CCA_UMAP', 'PCA_UMAP', 'PCA_BBKNN_BBKNN:UMAP'), 
                                      n.dims = 2, ground.truth.column = 'CellType', threads = 8, verbose = T)
  
}


