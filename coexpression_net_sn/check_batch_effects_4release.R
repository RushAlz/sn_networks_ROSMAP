library(tidyverse)
library(data.table)
library(WGCNA)
library(sva)
library(foreach)
library(doParallel)
library(lsa)

cosine_similarity_lower <- function(mat1, mat2) {
  idx <- lower.tri(mat1)
  v1 <- mat1[idx]
  v2 <- mat2[idx]
  sum(v1 * v2) / (sqrt(sum(v1^2)) * sqrt(sum(v2^2)))
}

cosine_similarity_lower_perm <- function(expr_A, expr_B, n_perm = 1000){
  # Observed correlation matrices
  cor_A <- cor(t(expr_A))
  cor_B <- cor(t(expr_B))
  
  # Observed cosine similarity
  obs_cosine <- cosine_similarity_lower(cor_A, cor_B)
  
  # Permutation loop to create a null distribution
  perm_cos = c()
  for(j in 1:1000){
    genes2use = which(rownames(tmm_voom)%in%rownames(tmm_voom)) # include genes in the module
    set.seed(j)
    genes_shuffled = sample(genes2use, length(modi_df$ensembl), replace = F)
    expr_perm <- tmm_voom[genes_shuffled,]
    expr_Batch_perm <- tmm_voom_rmBatch[genes_shuffled,]
    perm_cos = c(perm_cos,cosine_similarity_lower(cor(t(expr_perm)), cor(t(expr_Batch_perm))))
  }
  # One-sided p-value, testing if observed similarity is lower than random.
  pvalue <- sum(c(obs_cosine,perm_cos) <= obs_cosine) / length(c(obs_cosine,perm_cos))
  
  list(
    observed = obs_cosine,
    null_distribution = perm_cos,
    p_value = pvalue
  )
}

cell_names = c("Astrocytes", "Endothelial cells", "Excitatory Neurons", "Inhibitory Neurons", "Microglia", "Oligodendrocytes", "OPCs")
cell_names = setNames(cell_names, c("ast","end","ext","inh","mic","oli","opc"))

mypal = c("#1f77b4", # ast
          "#ff7f0e", # end
          "#2ca02c", # ext
          "#d62728", # inh
          "#9467bd", # mic
          "#8c564b", # opc
          "#e377c2"  # oli
)
mypal = setNames(mypal,c("ast","end","ext","inh","mic","opc","oli"))

se_results_dir = "/pastel/projects/speakeasy_dlpfc/SpeakEasy_singlenuclei/2nd_pass/snakemake-sn/results/"

# Load sn modules
sn_modules = data.frame()
for(cell_type in c("ext","inh","ast","oli","mic","end","opc")){
  # Filter modules with lest than 30 genes
  modules_file = read.table(paste0(se_results_dir,cell_type,"/geneBycluster.txt"), header = T)
  modules_file = modules_file[, c("ensembl", "gene_name", "cluster_lv3")]
  modules_file$module_clusters = paste0(cell_type,"_M", modules_file$cluster_lv3)
  modules_file$cluster_lv3 = NULL
  
  modules_file = modules_file %>% dplyr::group_by(module_clusters) %>%
    dplyr::mutate(module_size = n()) %>% dplyr::filter(module_size >= 30) %>%
    dplyr::mutate(cell_type = cell_type)
  
  sn_modules = rbind(sn_modules, modules_file)
}
# length(unique(sn_modules$module_clusters)) # 193
# table(unique(sn_modules[,c("module_clusters","cell_type")])$cell_type) # 7

batch_metadata = fread("/pastel/resources/20220203_snRNAseq_AMPAD/batch_metadata.txt", sep = "\t")
batch_metadata$projid = sprintf("%08d", batch_metadata$projid)
# length(unique(batch_metadata$batch)) # 60
# table(batch_metadata$batch)

load(paste0("/pastel/resources/20220203_snRNAseq_AMPAD/updated_annotations/celltype_exp_20220718.RData"))

results_by_cell_df = data.frame()

for(cell_i in cell_names){
  print(cell_i)
  
  gc()
  
  sn_modules_celli = sn_modules %>% filter(cell_type == names(cell_names)[cell_names == cell_i])
  donor_annotation_celli = batch_metadata %>% filter(cell.type == cell_i) %>% 
    dplyr::select(projid,batch) %>%
    distinct()
  
  # Correlation before adjustment
  tmm_voom = celltype_exp[[names(cell_names)[cell_names == cell_i]]]$tmm_voom[,donor_annotation_celli$projid]
  gene_cor_mat = WGCNA::cor(x = as.matrix(t(tmm_voom)))
  
  # Correlation after adjustment
  tmm_voom_rmBatch <- ComBat(dat = tmm_voom, batch = donor_annotation_celli$batch)
  gene_cor_mat_rmBatch = WGCNA::cor(x = as.matrix(t(tmm_voom_rmBatch)))

  # Loop over each module
  for(i in 1:length(unique(sn_modules_celli$module_clusters))){
    print(paste0(i,"/",length(unique(sn_modules_celli$module_clusters))))
    mod_i = unique(sn_modules_celli$module_clusters)[i]
    modi_df = sn_modules_celli %>% filter(module_clusters == mod_i)
    gene_cor_mat_modi = gene_cor_mat[modi_df$ensembl, modi_df$ensembl]

    perm_cosine_test_result = cosine_similarity_lower_perm(
      expr_A =  tmm_voom[modi_df$ensembl,],
      expr_B = tmm_voom_rmBatch[modi_df$ensembl,],
      n_perm = 1000
    )

    cat("Observed:", perm_cosine_test_result$observed, "\n")
    cat("Empirical p-value:", perm_cosine_test_result$p_value, "\n")

    results_by_cell_df = bind_rows(
      results_by_cell_df,
      data.frame(module = mod_i,
                 cosine_similarity = perm_cosine_test_result$observed,
                 cosine_similarity_p = perm_cosine_test_result$p_value,
                 module_size = nrow(modi_df)))
  }
}
results_by_cell_df$adj_p = p.adjust(results_by_cell_df$cosine_similarity_p, method = "bonferroni")
save(results_by_cell_df, file = "/pastel/projects/speakeasy_dlpfc/check_batch/check_batch_correction_4release.RData")

writexl::write_xlsx(results_by_cell_df[gtools::mixedorder(results_by_cell_df$module),], 
                    path = "/pastel/projects/speakeasy_dlpfc/check_batch/check_batch_correction_by_module_4release.xlsx")

  
