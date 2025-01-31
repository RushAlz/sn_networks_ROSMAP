# July 7, 2022
# Katia Lopes
# Knit the Eigengene reports

git_dir = "/pastel/Github_scripts/SpeakEasy_dlpfc/sn_dlpfc/2nd_pass/eigen_reports/"
# cell_type_list = c("ext", "inh", "oli", "ast", "mic", "opc", "end")
cell_type_list = c("Glia", "Neuronal")

render_report = function(element_list){
  rmarkdown::render(
    paste0(git_dir, "eigen_base_glia_neuron_sn.Rmd"), params = list(
      cell_type=element_list # region_type comes from the Rmd params. element_list is what changes. 
    ),
    output_file = paste0(git_dir, "Report_eigen_2nd_", element_list, "_lv4.html")
  )
}

for (i in cell_type_list){
  print(i)
  render_report(i)
}


