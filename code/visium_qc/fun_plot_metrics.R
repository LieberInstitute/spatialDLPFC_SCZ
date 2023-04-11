plot_metric <- function(spe, var_name,
                        include_log = TRUE){
  
  col_df <- colData(spe) |> data.frame()
  
  p <- ggplot(col_df) +
    geom_violin(aes(x = sample_id, y = !!sym(var_name), color = in_tissue))
  
  p_log <- ggplot(col_df) +
    geom_boxplot(aes(x = sample_id, y = !!sym(var_name), color = in_tissue)) +
    scale_y_log10()
  
  p_ret <- ggpubr::ggarrange(p, p_log,
                             nrow = 2, common.legend =  TRUE)
  
  if(!include_log)
    p_ret <- p
  
  ggsave(
    filename = file.path(fldr_qc_plots, 
                         paste0(var_name, "_per_sample.pdf")
    ),
    p_ret
  )
}