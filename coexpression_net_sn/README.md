# sn_networks_ROSMAP

***************************************
Single-nucleus RNASeq networks (n = 424 ROSMAP participants) 
*************************************** 

:small_orange_diamond: Data preparation and normalization:
[Click here](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/sn_phenotype.html) for cell proportions by individuo, plot.

:small_orange_diamond: Coexpression modules by cell type:

- [Excitatory neurons](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/net_reports/Report_ext.html)

- [Inhibitory neurons](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/net_reports/Report_inh.html)

- [Oligodendrocytes](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/net_reports/Report_oli.html)

- [Endothelial cells](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/net_reports/Report_end.html)

- [Astrocytes](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/net_reports/Report_ast.html)

- [Microglia](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/net_reports/Report_mic.html)

- [OPCs](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/net_reports/Report_opc.html)

:small_orange_diamond: Functional Enrichment Analysis (FEA):

| Network | Level | Module size | Perc genes assigned | FEA | Module eigengene | Heatmap for gene lists |
| ------------- | ------ | ------ | ------ | ------ | ------ | ------ |
| Ext| 3 | 29/33 | 99% | [grofiler](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/FEA/gprofiler/Report_gprofiler_cl3_ext.html) | [ME and AE](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/eigen_reports/Report_eigen_2nd_ext_lv3.html) | [Link](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/FEA/fisher/Report_heatmap_fisher_cl3_ext.html) |
| Inh | 3 | 24/53 | 99% | [grofiler](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/FEA/gprofiler/Report_gprofiler_cl3_inh.html) | [ME and AE](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/eigen_reports/Report_eigen_2nd_inh_lv3.html) | [Link](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/FEA/fisher/Report_heatmap_fisher_cl3_inh.html) |
| Oli | 3 | 30/38 | 99% | [grofiler](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/FEA/gprofiler/Report_gprofiler_cl3_oli.html) | [ME and AE](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/eigen_reports/Report_eigen_2nd_oli_lv3.html) | [Link](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/FEA/fisher/Report_heatmap_fisher_cl3_oli.html) |
| End | 3 | 26/295 | 95% | [grofiler](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/FEA/gprofiler/Report_gprofiler_cl3_end.html) | [ME and AE](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/eigen_reports/Report_eigen_2nd_end_lv3.html) | [Link](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/FEA/fisher/Report_heatmap_fisher_cl3_end.html) |
| Ast | 3 | 26/29 | 99% | [grofiler](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/FEA/gprofiler/Report_gprofiler_cl3_ast.html) | [ME and AE](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/eigen_reports/Report_eigen_2nd_ast_lv3.html) | [Link](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/FEA/fisher/Report_heatmap_fisher_cl3_ast.html) |
| Mic | 3 | 30/274 | 95% | [grofiler](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/FEA/gprofiler/Report_gprofiler_cl3_mic.html) | [ME and AE](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/eigen_reports/Report_eigen_2nd_mic_lv3.html) | [Link](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/FEA/fisher/Report_heatmap_fisher_cl3_mic.html) |
| OPCs | 3 | 28/59 | 99% | [grofiler](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/FEA/gprofiler/Report_gprofiler_cl3_opc.html) | [ME and AE](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/eigen_reports/Report_eigen_2nd_opc_lv3.html) | [Link](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/FEA/fisher/Report_heatmap_fisher_cl3_opc.html) |

One term for each module. The spreadsheet with manually curated results is [here](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/FEA/FEA_byModule_SN.xlsx), for all cell types.

[Summarize](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/FEA/gprofiler/gprofiler_summarize_sn.html) gprofiler for all the coexpression modules. 

[Hub genes](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/hub_genes_sn.html) for all the coexpression modules. Includes histograms checking the networks **topology**.

[Regression](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/lr_mod_cov_sn_cl3_adj_cov.html) between the modules and AD-related traits. 

Do the modules recapitulate sub-cell types clusters/states? [Results here](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/emods_vs_subcells_jaccard.html).

:small_orange_diamond: Module preservation

- [All major cell types + bulk + Mostafavi](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/MP/02_mp_bulk_mostafavi_majorCelltypes.html) pairwise. Includes heatmap with the number of non-preserved modules. You can check the percentage of **preserved** [HERE](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn//MP/02_mp_proportion_Aug2023.html).

:small_orange_diamond: Replication with an external single-nucleus dataset 

 - [MIT filtered](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/replication/replication_vs_MIT_filt.html) single-nucleus. Includes association heatmaps by cell type for the samples NOT used to create the modules.

 - [Showing modules](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/replication/replication_vs_MIT_filt_summarized2.html) that replicates in MIT samples NOT included to build the modules.

:small_orange_diamond: Bayesian networks

- [BN for the SN networks](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/BN/app07_BN_snRNA_hprior.html). Criteria = with mic.12, mic.13 and ast.10, runt = 500, eMod selection = bonf 0.1, impute 10 PCs, phenotypes downstream of everything else.

- [BN mic_M46](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/BN/BN_singlenuc_mic_m46.html)

- [BN ast_M19](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/BN/BN_singlenuc_ast_m19.html) for the top 100 associated genes.

- [BN inh_M06](https://rushalz.github.io/sn_networks_ROSMAP/coexpression_net_sn/BN/BN_singlenuc_inh_m6.html) for the top 100 associated genes, runt = 5000.



