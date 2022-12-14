library(IBRAP)
library(SingleCellExperiment)
library(splatter)
library(SymSim)
library(phytools)
library(Seurat)
library(SeuratData)

### Section 1 ###

IBRAP_objects <- list()

#  Download Pancreas Samples

InstallData('panc8')
pancreas <- LoadData('panc8')

for(x in unique(pancreas$dataset)) {
  
  IBRAP_objects[[x]] <- createIBRAPobject(counts = pancreas@assays$RNA@counts[,pancreas$dataset == x], original.project = x)
  IBRAP_objects[[x]]@sample_metadata$batch <- x
  IBRAP_objects[[x]]@sample_metadata$batch <- pancreas[,pancreas$dataset == x]$celltype
  
}

#  Cell Lines

#  Download Cell Lines, the files can be found at: https://github.com/LuyiTian/sc_mixology/tree/master/data

load('sincell_with_class.RData')
load('sincell_with_class_5cl.RData')

IBRAP_objects[['sc_10x']] <- createIBRAPobject(counts = counts(sce_sc_10x_qc), original.project = 'sc_10x')
IBRAP_objects[['sc_10x']]@sample_metadata$batch <- 'sc_10x'
IBRAP_objects[['sc_10x']]@sample_metadata$celltype <- colData(sce_sc_10x_qc)$cell_line

IBRAP_objects[['sc_10x_5cl']] <- createIBRAPobject(counts = counts(sce_sc_10x_5cl_qc), original.project = 'sc_10x_5cl')
IBRAP_objects[['sc_10x_5cl']]@sample_metadata$batch <- 'sc_10x_5cl'
IBRAP_objects[['sc_10x_5cl']]@sample_metadata$celltype <- colData(sce_sc_10x_5cl_qc)$cell_line

IBRAP_objects[['sc_CELseq2']] <- createIBRAPobject(counts = counts(sce_sc_CELseq2_qc), original.project = 'sc_CELseq2')
IBRAP_objects[['sc_CELseq2']]@sample_metadata$batch <- 'sc_CELseq2'
IBRAP_objects[['sc_CELseq2']]@sample_metadata$celltype <- colData(sce_sc_CELseq2_qc)$cell_line

IBRAP_objects[['sc_Dropseq']] <- createIBRAPobject(counts = counts(sce_sc_Dropseq_qc), original.project = 'sc_Dropseq')
IBRAP_objects[['sc_Dropseq']]@sample_metadata$batch <- 'sc_Dropseq'
IBRAP_objects[['sc_Dropseq']]@sample_metadata$celltype <- colData(sce_sc_Dropseq_qc)$cell_line

IBRAP_objects[['Celseq2_5cl_p1']] <- createIBRAPobject(counts = counts(sc_Celseq2_5cl_p1), original.project = 'Celseq2_5cl_p1')
IBRAP_objects[['Celseq2_5cl_p1']]@sample_metadata$batch <- 'Celseq2_5cl_p1'
IBRAP_objects[['Celseq2_5cl_p1']]@sample_metadata$celltype <- colData(sc_Celseq2_5cl_p1)$cell_line

IBRAP_objects[['Celseq2_5cl_p2']] <- createIBRAPobject(counts = counts(sc_Celseq2_5cl_p2), original.project = 'Celseq2_5cl_p2')
IBRAP_objects[['Celseq2_5cl_p2']]@sample_metadata$batch <- 'Celseq2_5cl_p2'
IBRAP_objects[['Celseq2_5cl_p2']]@sample_metadata$celltype <- colData(sc_Celseq2_5cl_p2)$cell_line

IBRAP_objects[['Celseq2_5cl_p3']] <- createIBRAPobject(counts = counts(sc_Celseq2_5cl_p3), original.project = 'sc_CELseq2')
IBRAP_objects[['Celseq2_5cl_p3']]@sample_metadata$batch <- 'Celseq2_5cl_p3'
IBRAP_objects[['Celseq2_5cl_p3']]@sample_metadata$celltype <- colData(sc_Celseq2_5cl_p3)$cell_line

# Splatter simulations

params <- newSplatParams()

params <- setParam(params, "nGenes", 20000)

params <- setParams(params, mean.shape = 0.5, de.prob = 0.2)

params <- newSplatParams(lib.loc = 12, lib.scale = 0.7)

sim.8_unequal <- splatSimulate(params, 
                               batchCells = 1000, 
                               group.prob = c(0.20,0.15,0.1,0.05,0.30,0.02,0.07,0.11), 
                               method = "groups")

sim.4_unequal <- splatSimulate(params, 
                               batchCells = 1000, 
                               group.prob = c(0.3,0.18,0.42,0.1), 
                               method = "groups")

sim.8_equal <- splatSimulate(params, 
                             batchCells = 1000, 
                             group.prob = c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125), 
                             method = "groups")

sim.4_equal <- splatSimulate(params, 
                             batchCells = 1000, 
                             group.prob = c(0.25,0.25,0.25,0.25), 
                             method = "groups")

IBRAP_objects[['sim.8_unequal']] <- createIBRAPobject(counts = counts(sim.8_unequal), original.project = 'sim.8_unequal')
IBRAP_objects[['sim.8_equal']] <- createIBRAPobject(counts = counts(sim.8_equal), original.project = 'sim.8_unequal')
IBRAP_objects[['sim.4_unequal']] <- createIBRAPobject(counts = counts(sim.4_unequal), original.project = 'sim.8_unequal')
IBRAP_objects[['sim.4_equal']] <- createIBRAPobject(counts = counts(sim.4_equal), original.project = 'sim.8_unequal')

#  SymSim simulations

ngenes <- 10000

data(gene_len_pool)
gene_len <- sample(gene_len_pool, 10000, replace = FALSE)

counts <- SimulateTrueCounts(ncells_total=3000, min_popsize=50, ngenes=ngenes, nevf=30, 
                             evf_type="discrete", n_de_evf=12, vary="s", Sigma=0.4, 
                             phyla=Phyla5(), randseed=0)

UMI_highdepth_5 <- True2ObservedCounts(true_counts=counts[[1]], meta_cell=counts[[3]], 
                                       protocol="UMI", alpha_mean=0.04, alpha_sd=0.2, 
                                       gene_len=gene_len, depth_mean=5e7, depth_sd=3e4)

counts <- SimulateTrueCounts(ncells_total=3000, min_popsize=50, ngenes=ngenes, nevf=30, 
                             evf_type="discrete", n_de_evf=12, vary="s", Sigma=0.4, 
                             phyla=Phyla5(), randseed=0)

UMI_lowdepth_5 <- True2ObservedCounts(true_counts=counts[[1]], meta_cell=counts[[3]], 
                                      protocol="UMI", alpha_mean=0.04, alpha_sd=0.2, 
                                      gene_len=gene_len, depth_mean=5e5, depth_sd=3e3)

counts <- SimulateTrueCounts(ncells_total=3000, min_popsize=50, ngenes=ngenes, nevf=30, 
                             evf_type="discrete", n_de_evf=12, vary="s", Sigma=0.4, 
                             phyla=Phyla3(), randseed=0)

UMI_highdepth_3 <- True2ObservedCounts(true_counts=counts[[1]], meta_cell=counts[[3]], 
                                       protocol="UMI", alpha_mean=0.04, alpha_sd=0.2, 
                                       gene_len=gene_len, depth_mean=5e7, depth_sd=3e4)

counts <- SimulateTrueCounts(ncells_total=3000, min_popsize=50, ngenes=ngenes, nevf=30, 
                             evf_type="discrete", n_de_evf=12, vary="s", Sigma=0.4, 
                             phyla=Phyla3(), randseed=0)

UMI_lowdepth_3 <- True2ObservedCounts(true_counts=counts[[1]], meta_cell=counts[[3]], 
                                      protocol="UMI", alpha_mean=0.04, alpha_sd=0.2, 
                                      gene_len=gene_len, depth_mean=5e5, depth_sd=3e3)

counts <- UMI_highdepth_5$counts
rownames(counts) <- paste0('gene-', 1:10000)
colnames(counts) <- UMI_highdepth_5$cell_meta$cellid

IBRAP_objects[['UMI_highdepth_5']] <- createIBRAPobject(counts = counts, original.project = 'UMI_highdepth_5')
IBRAP_objects[['UMI_highdepth_5']]@sample_metadata$celltype <- UMI_highdepth_5$cell_meta$pop

counts <- UMI_lowdepth_5$counts
rownames(counts) <- paste0('gene-', 1:10000)
colnames(counts) <- UMI_lowdepth_5$cell_meta$cellid

IBRAP_objects[['UMI_lowdepth_5']] <- createIBRAPobject(counts = counts, original.project = 'UMI_lowdepth_5')
IBRAP_objects[['UMI_lowdepth_5']]@sample_metadata$celltype <- UMI_lowdepth_5$cell_meta$pop

counts <- UMI_highdepth_3$counts
rownames(counts) <- paste0('gene-', 1:10000)
colnames(counts) <- UMI_highdepth_3$cell_meta$cellid

IBRAP_objects[['UMI_highdepth_3']] <- createIBRAPobject(counts = counts, original.project = 'UMI_highdepth_3')
IBRAP_objects[['UMI_highdepth_3']]@sample_metadata$celltype <- UMI_highdepth_3$cell_meta$pop

counts <- UMI_highdepth_5$counts
rownames(counts) <- paste0('gene-', 1:10000)
colnames(counts) <- UMI_lowdepth_3$cell_meta$cellid

IBRAP_objects[['UMI_lowdepth_3']] <- createIBRAPobject(counts = counts, original.project = 'UMI_lowdepth_3')
IBRAP_objects[['UMI_lowdepth_3']]@sample_metadata$celltype <- UMI_lowdepth_3$cell_meta$pop

### SECTION 2 ###

integrated_samples <- list()

# Pancreas Samples

integrated_samples[['analysis_1']] <- merge(IBRAP_objects[[1]], c(IBRAP_objects[[2]], IBRAP_objects[[3]], IBRAP_objects[[4]],
                                                                  IBRAP_objects[[5]], IBRAP_objects[[6]], IBRAP_objects[[7]], IBRAP_objects[[8]]))

# Cell lines = sc_CELseq2, sc_Dropseq, and sc_10x

integrated_samples[['analysis_2']] <- merge(IBRAP_objects[[9]], c(IBRAP_objects[[10]], IBRAP_objects[[12]]))

# Cell lines = sc_CELseq2, sc_Dropseq, and sc_10x_5cl

integrated_samples[['analysis_3']] <- merge(IBRAP_objects[[9]], c(IBRAP_objects[[10]], IBRAP_objects[[11]]))

batch_dataset <- SimulateTrueCounts(ncells_total=3000, min_popsize=50, ngenes=ngenes, nevf=30, 
                                    evf_type="discrete", n_de_evf=12, vary="s", Sigma=0.4, 
                                    phyla=Phyla5(), randseed=0)

batch_dataset_highdepth_5 <- True2ObservedCounts(true_counts=batch_dataset[[1]], meta_cell=batch_dataset[[3]], 
                                                 protocol="UMI", alpha_mean=0.04, alpha_sd=0.2, 
                                                 gene_len=gene_len, depth_mean=5e7, depth_sd=3e4)

#  Symsim simulations

integrated_samples[['analysis_4']] <- DivideBatches(observed_counts_res = batch_dataset_highdepth_5, 
                                                    nbatch = 3, batch_effect_size = 1)

integrated_samples[['analysis_5']] <- DivideBatches(observed_counts_res = batch_dataset_highdepth_5, 
                                                    nbatch = 3, batch_effect_size = 5)

### SECTION 3 ###

supervised_clustering <- list()

reference <- panc[,panc$dataset!='celseq']
query <- panc[,panc$dataset=='celseq']

supervised_clustering[['celseq_reference']] <- createIBRAPobject(counts = round(reference@assays$RNA@counts), 
                                                                    original.project = 'celseq_reference', 
                                                                    meta.data = reference@meta.data[,c('celltype','dataset')])

supervised_clustering[['celseq_query']] <- createIBRAPobject(counts = round(query@assays$RNA@counts), 
                                                                original.project = 'celseq_query', 
                                                                meta.data = query@meta.data[,c('celltype','dataset')])

reference <- pancreas[,panc$dataset!='celseq2']
query <- pancreas[,panc$dataset=='celseq2']

supervised_clustering[['celseq2_reference']] <- createIBRAPobject(counts = round(reference@assays$RNA@counts), 
                                                                 original.project = 'celseq2_reference', 
                                                                 meta.data = reference@meta.data[,c('celltype','dataset')])

supervised_clustering[['celseq2_query']] <- createIBRAPobject(counts = round(query@assays$RNA@counts), 
                                                             original.project = 'celseq2_query', 
                                                             meta.data = query@meta.data[,c('celltype','dataset')])

reference <- pancreas[,pancreas$dataset!='smartseq2']
query <- pancreas[,pancreas$dataset=='smartseq2']

supervised_clustering[['smartseq2_reference']] <- createIBRAPobject(counts = round(reference@assays$RNA@counts), 
                                                                    original.project = 'smartseq2_reference', 
                                                                    meta.data = reference@meta.data[,c('celltype','dataset')])

supervised_clustering[['smartseq2_query']] <- createIBRAPobject(counts = round(query@assays$RNA@counts), 
                                                                original.project = 'smartseq2_query', 
                                                                meta.data = query@meta.data[,c('celltype','dataset')])

reference <- pancreas[,pancreas$dataset!='fluidigmc1']
query <- pancreas[,pancreas$dataset=='fluidigmc1']

supervised_clustering[['fluidigmc1_reference']] <- createIBRAPobject(counts = round(reference@assays$RNA@counts), 
                                                                    original.project = 'fluidigmc1_reference', 
                                                                    meta.data = reference@meta.data[,c('celltype','dataset')])

supervised_clustering[['fluidigmc1_query']] <- createIBRAPobject(counts = round(query@assays$RNA@counts), 
                                                                original.project = 'fluidigmc1_query', 
                                                                meta.data = query@meta.data[,c('celltype','dataset')])

reference <- pancreas[,pancreas$dataset!='indrop1']
query <- pancreas[,pancreas$dataset=='indrop1']

supervised_clustering[['indrop1_reference']] <- createIBRAPobject(counts = round(reference@assays$RNA@counts), 
                                                                     original.project = 'indrop1_reference', 
                                                                     meta.data = reference@meta.data[,c('celltype','dataset')])

supervised_clustering[['indrop1_query']] <- createIBRAPobject(counts = round(query@assays$RNA@counts), 
                                                                 original.project = 'indrop1_query', 
                                                                 meta.data = query@meta.data[,c('celltype','dataset')])

reference <- pancreas[,pancreas$dataset!='indrop1']
query <- pancreas[,pancreas$dataset=='indrop1']

supervised_clustering[['indrop1_reference']] <- createIBRAPobject(counts = round(reference@assays$RNA@counts), 
                                                                  original.project = 'indrop1_reference', 
                                                                  meta.data = reference@meta.data[,c('celltype','dataset')])

supervised_clustering[['indrop1_query']] <- createIBRAPobject(counts = round(query@assays$RNA@counts), 
                                                              original.project = 'indrop1_query', 
                                                              meta.data = query@meta.data[,c('celltype','dataset')])

reference <- pancreas[,pancreas$dataset!='indrop2']
query <- pancreas[,pancreas$dataset=='indrop2']

supervised_clustering[['indrop2_reference']] <- createIBRAPobject(counts = round(reference@assays$RNA@counts), 
                                                                  original.project = 'indrop2_reference', 
                                                                  meta.data = reference@meta.data[,c('celltype','dataset')])

supervised_clustering[['indrop2_query']] <- createIBRAPobject(counts = round(query@assays$RNA@counts), 
                                                              original.project = 'indrop2_query', 
                                                              meta.data = query@meta.data[,c('celltype','dataset')])


reference <- pancreas[,pancreas$dataset!='indrop3']
query <- pancreas[,pancreas$dataset=='indrop3']

supervised_clustering[['indrop3_reference']] <- createIBRAPobject(counts = round(reference@assays$RNA@counts), 
                                                                  original.project = 'indrop3_reference', 
                                                                  meta.data = reference@meta.data[,c('celltype','dataset')])

supervised_clustering[['indrop3_query']] <- createIBRAPobject(counts = round(query@assays$RNA@counts), 
                                                              original.project = 'indrop3_query', 
                                                              meta.data = query@meta.data[,c('celltype','dataset')])

reference <- pancreas[,pancreas$dataset!='indrop4']
query <- pancreas[,pancreas$dataset=='indrop4']

supervised_clustering[['indrop4_reference']] <- createIBRAPobject(counts = round(reference@assays$RNA@counts), 
                                                                  original.project = 'indrop4_reference', 
                                                                  meta.data = reference@meta.data[,c('celltype','dataset')])

supervised_clustering[['indrop4_query']] <- createIBRAPobject(counts = round(query@assays$RNA@counts), 
                                                              original.project = 'indrop4_query', 
                                                              meta.data = query@meta.data[,c('celltype','dataset')])

#  saving data

saveRDS(IBRAP_objects, 'individual_samples.rds')

saveRDS(integrated_samples, 'integrated_samples.rds')

saveRDS(supervised_clustering, 'reference_samples.rds')




