plot_metric <- function(spe, var_name,
                        include_log = TRUE){
  
  thres_list <- metadata(spe)[grepl(paste0(var_name, "_thres"),
                                    names(metadata(spe)))]
  
  if(length(thres_list) > 1) {
    warning("More than 1 matched var_name in the metadata(spe).",
            "Default to first match ", names(thres_list)[1])
  }
  
  if(length(thres_list) > 0){
    
    thres <- thres_list[[1]]
    
    # browser()
    
    # Hard coding mt thres selection
    if(var_name == "expr_chrM_ratio"){
      thres <- thres["higher",] 
    } else {
      thres <- thres["lower",]
    }
    
    
    
    # Make thres matrix long form
    thres_df <- thres |> data.frame() |> 
      rownames_to_column("sample_id")
    
    
    # thres_point <- geom_point(
    #   aes(x = sample_id, y = thres, shape = "5"),
    #   data = thres_df
    # ) + scale_shape_manual(
    #   labels = "scuttle::isOutlier Thres"
    # )
  }
  
  col_df <- colData(spe) |> data.frame()
  
  # browser()
  
  p <- ggplot(col_df) +
    geom_violin(aes(x = sample_id, y = !!sym(var_name), color = in_tissue))
  
  
  
  p_log <- ggplot(col_df) +
    geom_boxplot(aes(x = sample_id, y = !!sym(var_name), color = in_tissue)) +
    scale_y_log10()
  
  if(length(thres_list) > 0){
    p <- p  +
      geom_point(
        aes(x = sample_id, y = thres, shape = "5"),
        data = thres_df
      ) + scale_shape_manual(
        name = "scuttle thres",
        values = c("5" = 5),
        labels = c("5" = "")
      )
    
    p_log <- p_log +
      geom_point(
        aes(x = sample_id, y = thres, shape = "5"),
        data = thres_df
      ) + scale_shape_manual(
        name = "scuttle thres",
        values = c("5" = 5),
        labels = c("5" = "")
      )
  }
  
  
  p_ret <- ggpubr::ggarrange(p, p_log,
                             nrow = 2, common.legend =  TRUE)
  
  if(!include_log){
    p_ret <- p
    # if(thre_exist){
    #   p_ret <- p +
    # }
    # 
  }
  
  ggsave(
    filename = file.path(fldr_qc_plots, 
                         paste0(var_name, "_per_sample.pdf")
    ),
    p_ret
  )
}