library(IBRAP)

objects <- readRDS('reference_samples.rds')

objects$celseq_query <- perform.singleR.annotation(object = objects$celseq_query, assay = c('SCT','SCRAN','SCANPY'), 
                                                   slot = 'counts', ref = objects$celseq_reference[['RAW']][['counts']], 
                                                   ref.labels = objects$celseq_reference[['celltype']])

objects$celseq2_query <- perform.singleR.annotation(object = objects$celseq2_query, assay = c('SCT','SCRAN','SCANPY'), 
                                                    slot = 'counts', ref = objects$celseq2_reference[['RAW']][['counts']], 
                                                    ref.labels = objects$celseq2_reference[['celltype']])

objects$smartseq2_query <- perform.singleR.annotation(object = objects$smartse2q_query, assay = c('SCT','SCRAN','SCANPY'), 
                                                      slot = 'counts', ref = objects$smartseq2_reference[['RAW']][['counts']], 
                                                      ref.labels = objects$smartseq2_reference[['celltype']])

objects$fluidigmc1_query <- perform.singleR.annotation(object = objects$fluidigmc1_query, assay = c('SCT','SCRAN','SCANPY'), 
                                                       slot = 'counts', ref = objects$fluidigmc1_reference[['RAW']][['counts']], 
                                                       ref.labels = objects$fluidigmc1_reference[['celltype']])

objects$indrop1_query <- perform.singleR.annotation(object = objects$indrop1_query, assay = c('SCT','SCRAN','SCANPY'), 
                                                    slot = 'counts', ref = objects$indrop1_reference[['RAW']][['counts']], 
                                                    ref.labels = objects$indrop1_reference[['celltype']])

objects$indrop2_query <- perform.singleR.annotation(object = objects$indrop2_query, assay = c('SCT','SCRAN','SCANPY'), 
                                                    slot = 'counts', ref = objects$indrop2_reference[['RAW']][['counts']], 
                                                    ref.labels = objects$indrop2_reference[['celltype']])

objects$indrop3_query <- perform.singleR.annotation(object = objects$indrop3_query, assay = c('SCT','SCRAN','SCANPY'), 
                                                    slot = 'counts', ref = objects$indrop3_reference[['RAW']][['counts']], 
                                                    ref.labels = objects$indrop3_reference[['celltype']])

objects$indrop4_query <- perform.singleR.annotation(object = objects$indrop4_query, assay = c('SCT','SCRAN','SCANPY'), 
                                                    slot = 'counts', ref = objects$indrop4_reference[['RAW']][['counts']], 
                                                    ref.labels = objects$indrop4_reference[['celltype']])

rownames(tmp@sample_metadata) <- colnames(individual_samples$smartseq2)

tmp@methods$SCT <- individual_samples$smartseq2@methods$SCT

tmp@sample_metadata$optimal_clusters <- tmp@methods$SCT@cluster_assignments$`PCA_NN:SLM`$neighbourhood_graph_res.0.1

tmp@sample_metadata$sctype_query<- tmp@methods$SCT@cluster_assignments$`PCA_NN:SLM`$neighbourhood_graph_res.0.1

perform.sctype(object = tmp, clust.method = 'metadata', column = 'sctype_query', tissue = 'Pancreas')

tmp@sample_metadata$optimal_clusters <- as.character(tmp@sample_metadata$optimal_clusters)
tmp@sample_metadata$optimal_clusters[tmp@sample_metadata$optimal_clusters=='0']='alpha_1'
tmp@sample_metadata$optimal_clusters[tmp@sample_metadata$optimal_clusters=='3']='alpha_2'
tmp@sample_metadata$optimal_clusters[tmp@sample_metadata$optimal_clusters=='2']='beta'
tmp@sample_metadata$optimal_clusters[tmp@sample_metadata$optimal_clusters=='6']='delta'
tmp@sample_metadata$optimal_clusters[tmp@sample_metadata$optimal_clusters=='1']='ductal'
tmp@sample_metadata$optimal_clusters[tmp@sample_metadata$optimal_clusters=='5']='acinar'
tmp@sample_metadata$optimal_clusters[tmp@sample_metadata$optimal_clusters=='7']='mixed_minor_celltypes'
tmp@sample_metadata$optimal_clusters[tmp@sample_metadata$optimal_clusters=='4']='gamma'

tmp <- perform.scran(reference_results$smartseq2_reference, verbose = T)
tmp <- perform.sct(reference_results$smartseq2_reference, verbose = T)
tmp <-perform.pca(tmp, assay='SCRAN', verbose = T)
tmp <- perform.bbknn(object = tmp, assay = c('SCT','SCRAN'), reduction = 'PCA', n_pcs = list(0), verbose = T, batch = 'dataset')
tmp <- perform.harmony(object = tmp, assay = 'SCT', batch = 'dataset', reduction = 'PCA', dims.use = list(10), verbose = T)
tmp <- perform.umap(tmp, assay = c('SCT','SCRAN'), reduction = 'PCA_HARMONY', dims.use = list(0))
tmp <- perform.umap(tmp, assay = c('SCT','SCRAN'), graph = 'PCA_BBKNN_BBKNN')
tmp<- perform.nn(object = tmp, assay = 'SCT', reduction = 'PCA_HARMONY', dims.use = list(10))
tmp <- perform.graph.cluster(object = tmp, assay = 'SCT', verbose = T, neighbours = 'PCA_HARMONY_NN')

tmp <- benchmark.clustering(object = tmp, assay = 'SCT', clustering = 'PCA_HARMONY_NN:LOUVAIN', 
                            n.dims = 2, reduction = 'PCA_HARMONY_UMAP', ground.truth.column = 'celltype')

tmp <- perform.sctype(object = tmp, clust.method = 'metadata', column = 'sctype_query', tissue = 'Pancreas')

