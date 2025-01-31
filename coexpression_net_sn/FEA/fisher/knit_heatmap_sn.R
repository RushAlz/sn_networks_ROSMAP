# Oct 24, 2022
# Katia Lopes
# Knit the gprofiler pages

git_dir = "/pastel/Github_scripts/SpeakEasy_dlpfc/sn_dlpfc/2nd_pass/FEA/heatmap/"

cell_type_list = c("ext", "inh", "oli", "ast", "mic", "opc", "end")
# cell_type_list = c("Glia", "Neuronal")

render_report = function(element_list){
  rmarkdown::render(
    paste0(git_dir, "base_heatmap_sn.Rmd"), params = list(
      cell_type=element_list # cell_type comes from the Rmd params. element_list is what changes. 
    ),
    output_file = paste0(git_dir, "Report_heatmap_fisher_cl3_", element_list, ".html")
  )
}

for (i in cell_type_list){
  print(i)
  render_report(i)
}
