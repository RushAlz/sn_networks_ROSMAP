library(tidyverse)
library(WGCNA)
library(vegan)  
library(lsa)
library(tidyverse)
library(data.table)
library(sva)
library(writexl)
library(readxl)
 
GLOBAL_SEED = 12345
set.seed(GLOBAL_SEED)

# Cosine similarity using lower triangle only
cosine_similarity_lower <- function(mat1, mat2) {
  idx <- lower.tri(mat1)
  v1 <- mat1[idx]
  v2 <- mat2[idx]
  sum(v1 * v2) / (sqrt(sum(v1^2)) * sqrt(sum(v2^2)))
}

# Function to analyze a single module within a cell type
cosine_mantel_module <- function(expr_A, expr_B, module_genes, seed = 123) {
  set.seed(seed)
  
  # Subset expression matrices to module genes
  expr_A_sub <- expr_A[module_genes, , drop = FALSE]
  expr_B_sub <- expr_B[module_genes, , drop = FALSE]
  
  # Compute correlation matrices (samples x samples)
  cor_A <- cor(t(expr_A_sub))
  cor_B <- cor(t(expr_B_sub))
  
  # Cosine similarity (lower triangle only)
  cosine <- cosine_similarity_lower(cor_A, cor_B)
  
  # Mantel test (correlation of correlation matrices)
  mantel_result <- mantel(cor_A, cor_B, method = "pearson", permutations = 1000)
  
  list(
    observed_cosine = cosine,
    mantel_statistic = mantel_result$statistic,
    mantel_p = mantel_result$signif
  )
}

################################################################################

# Cell type names and palette
cell_names <- c("Astrocytes", "Endothelial cells", "Excitatory Neurons", 
                "Inhibitory Neurons", "Microglia", "Oligodendrocytes", "OPCs")
cell_ids = c("ast","end","ext","inh","mic","oli","opc")
cell_names <- setNames(cell_names, cell_ids)

mypal <- setNames(c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2"),
                  cell_ids)

se_results_dir <- "/pastel/projects/speakeasy_dlpfc/SpeakEasy_singlenuclei/2nd_pass/snakemake-sn/results/"

################################################################################

# Load modules
load_modules <- function() {
  cell_types <- cell_ids
  modules_list <- lapply(cell_types, function(cell_type) {
    modules_file <- fread(paste0(se_results_dir, cell_type, "/geneBycluster.txt"))
    modules_file[, module_clusters := paste0(cell_type, "_M", cluster_lv3)]
    modules_file[, cluster_lv3 := NULL]
    modules_file[, module_size := .N, by = module_clusters]
    modules_file <- modules_file[module_size >= 30]
    modules_file[, cell_type := cell_type]
    return(modules_file)
  })
  rbindlist(modules_list)
}

sn_modules <- load_modules()

load("/pastel/resources/20220203_snRNAseq_AMPAD/updated_annotations/celltype_exp_20220718.RData")

all_donors = colnames(celltype_exp[["ext"]]$tmm_voom)

syn_idmap = read_tsv("/pastel/resources/synapse/SynapseID_mapping.tsv", show_col_types = F)
batch_metadata = readxl::read_xlsx("/pastel/resources/20220203_snRNAseq_AMPAD/41586_2024_7871_MOESM5_ESM.xlsx", 
                                   sheet = "Participant inclusion QCs", col_names = TRUE)
batch_metadata = batch_metadata %>% left_join(syn_idmap, by = c("individualID" = "synapseid")) 
batch_metadata = batch_metadata %>% filter(Final_QC == "Pass") %>% select(projid, batch) %>% distinct() 
batch_metadata = batch_metadata %>% filter(projid %in% all_donors) %>% distinct() %>% as.data.table()

################################################################################

results_list = list()
for(cell_i in cell_names){
  cat("Processing:", cell_i, "\n")

  cell_key <- names(cell_names)[cell_names == cell_i]
  sn_modules_celli <- sn_modules[cell_type == cell_key]
  donor_annotation_celli <- batch_metadata[, .(projid, batch)]
  donor_annotation_celli <- unique(donor_annotation_celli)
  
  available_donors <- intersect(donor_annotation_celli$projid, 
                                colnames(celltype_exp[[cell_key]]$tmm_voom))
  donor_annotation_celli <- donor_annotation_celli[projid %in% available_donors]
  
  tmm_voom <- celltype_exp[[cell_key]]$tmm_voom[, donor_annotation_celli$projid, drop = FALSE]
  
  combat_seed <- GLOBAL_SEED + which(cell_ids == cell_key) * 1000
  set.seed(combat_seed)
  tmm_voom_rmBatch <- ComBat(dat = tmm_voom, batch = donor_annotation_celli$batch, mod = NULL)
  
  unique_modules <- unique(sn_modules_celli$module_clusters)
  n_modules <- length(unique_modules)
  
  cat("  Processing", n_modules, "modules\n")
  
  for(i in seq_along(unique_modules)){
    cat("    Module", i, "/", n_modules, "\n")
    
    mod_i <- unique_modules[i]
    modi_df <- sn_modules_celli[module_clusters == mod_i]
    
    available_genes <- intersect(modi_df$ensembl, rownames(tmm_voom))
    
    expr_mod <- tmm_voom[available_genes, , drop = FALSE]
    expr_mod_batch <- tmm_voom_rmBatch[available_genes, , drop = FALSE]
    
    module_seed <- GLOBAL_SEED + which(cell_ids == cell_key) * 1000 + i
    
    perm_result <- cosine_mantel_module(expr_A = tmm_voom, 
                                        expr_B = tmm_voom_rmBatch, 
                                        module_genes = available_genes, 
                                        seed = module_seed)
    
    results_list[[mod_i]] <- data.frame(
      module = mod_i,
      cell_type = cell_key,
      cosine_similarity = perm_result$observed_cosine,
      mantel_stat = perm_result$mantel_statistic,
      mantel_p = perm_result$mantel_p,
      module_size = length(available_genes),
      seed_used = module_seed,
      stringsAsFactors = FALSE
    )
  }
  
}

results_by_cell_df <- rbindlist(results_list[!sapply(results_list, is.null)])
results_by_cell_df[, adj_p := p.adjust(mantel_p, method = "fdr")]
results_ordered <- results_by_cell_df[gtools::mixedorder(results_by_cell_df$module)]

save(results_by_cell_df, file = "/pastel/projects/speakeasy_dlpfc/check_batch/check_batch_correction_4release.RData")
writexl::write_xlsx(results_ordered,
                    path = "/pastel/projects/speakeasy_dlpfc/check_batch/check_batch_correction_by_module_4release.xlsx")


# # Check if any of the modules assocaited with phenotypes are affected by batch
# load("/pastel/Github_scripts/SpeakEasy_dlpfc/sn_dlpfc/2nd_pass/eigen_reports/save_lr_adjcov/all_res_test_stats_SN.Rdata")
# all_stats = all_stats %>% mutate(module = paste0(network, "_M", gsub("AE","",module)))
# # Remove non-modules 
# all_stats_filt <- all_stats[all_stats$module %in% results_by_cell_df$module, ]
# all_stats_filt$network = as.factor(all_stats_filt$network)
# 
# pheno_variable = c("cogng_demog_slope","amyloid_sqrt","tangles_sqrt")
# pval_cutoff = 0.05
# all_stats_bonf = all_stats_filt %>% 
#   filter(phenotype %in% pheno_variable) %>%
#   group_by(phenotype) %>% # Adjust by phenotype separately
#   mutate(bonf_p = p.adjust(nom_p, method = "bonferroni")) %>% 
#   filter(bonf_p <= pval_cutoff) %>% # Bonferroni correction for all tests (modules X phenotypes; 1544 tests)
#   arrange(bonf_p)
# 
# results_ordered = results_ordered %>% mutate(is_assoc_with_pheno = module %in% all_stats_bonf$module)
# results_ordered[results_ordered$is_assoc_with_pheno,] %>% arrange(mantel_stat)

