spe <- raw_spe


col_df <- colData(spe) |> data.frame()

thres_df <- col_df |> filter(in_tissue == FALSE) |> 
  group_by(sample_id) |> 
  summarize(
    sum_umi_thres = median(sum_umi)
  )

col_df <- col_df |> 
  left_join(thres_df, by = "sample_id") |> 
  mutate(outlier = factor(sum_umi <= sum_umi_thres))

colData(spe)$outlier <- col_df$outlier


library(escheR)
make_escheR(spe[, spe$sample_id == "Br2378_D1" & spe$in_tissue == TRUE]) |> 
  add_fill(var = "sum_umi") |> 
  add_ground(var = "outlier", stroke = 0.2) +
  scale_fill_viridis_c(trans="log10") +
  scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "transparent"))


make_escheR(spe[, spe$sample_id == "Br5182_C1"]) |> 
  add_fill(var = "sum_umi") |> 
  add_ground(var = "in_tissue", stroke = 0.2) +
  scale_fill_viridis_c(trans="log10") +
  scale_colour_manual(values = c("grey", "transparent"))


# * Out_tissue Comparison -------------------------------------------------

col_df <- colData(spe) |> data.frame() |> 
  filter(in_tissue == FALSE) |> 
  group_by(sample_id) |> 
  summarize(
    lib_median = median(sum_umi),
    lib_low_quart = quantile(sum_umi, probs = 0.25)
  )

tmp_df <- colData(spe_in_tissue) |> data.frame()
tmp_df <- tmp_df |> left_join(col_df, by = "sample_id") |> 
  mutate(dq_lib_median = factor(sum_umi <= lib_median),
         dp_lib_low_quart = factor(sum_umi <= lib_low_quart))

spe_in_tissue$dq_lib_median <- tmp_df$dq_lib_median
spe_in_tissue$dq_lib_low_quart <- tmp_df$dq_lib_low_quart

make_escheR(spe_in_tissue[, spe_in_tissue$sample_id == "Br1958_B1"]) |> 
  add_fill(var = "sum_umi") |> 
  add_ground(var = "dq_lib_median", stroke = 0.5) +
  scale_fill_viridis_c(trans="log10") +
  scale_colour_manual(values = c("transparent", "red"))