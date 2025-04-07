# Debug finalized SPE object


PRECAST_df <- readRDS(
  here(
    "processed-data/rds/03_visium_spatial_clustering",
    "PRECAST_label_df_semi_sup_k_2-16.rds"
  )
)

PRECAST_df <- PRECAST_df %>%
  separate(key, into = c("part1", "part2"), sep = "_", extra = "merge", remove = FALSE)


  table_result <- table(PRECAST_df$part2, PRECAST_df$PRECAST_07)
  print(table_result)


# This shows that two datas sets should have very similar rows.
identical(table(PRECAST_df$part2), table(spe$sample_id))


# This shows that the merge is not working as expected.
identical(
  table(PRECAST_df$part2, PRECAST_df$PRECAST_07),
  table(spe$sample_id , spe$PRECAST_07)
)

table(spe$sample_id , spe$PRECAST_07) |> rowSums()

setdiff(PRECAST_df$key, spe$key)


with(col_data_df, table(sample_id, PRECAST_07.y))

col_data_df |> filter(PRECAST_07.x != PRECAST_07.y) |>filter(sample_id == "V13M06-344_C1") |> 
select(key, sample_id, PRECAST_07.x, PRECAST_07.y)

spe |> colData() |> data.frame() |> filter(sample_id == "V13M06-344_C1") |> 
  select(key, sample_id, PRECAST_07) |> head()




colData(spe) |> data.frame() |> str()


spe <- readRDS(
  here(
    "processed-data/rds/01_build_spe",
    "fnl_spe_kept_spots_only.rds"
  )
)