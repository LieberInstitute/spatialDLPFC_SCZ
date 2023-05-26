############## NOTE ###################################################
#  The code is edited to only save SCZ space ranger results
####################################################################

library("here")
library("dplyr")
library("sessioninfo")

## Rscript -e 'sgejobs::job_single("01_read_spaceranger_metrics", create_shell = TRUE, queue = "bluejay")'

## Output dirs
dir_rdata <-
    here(
        "processed-data",
        "spaceranger",
        "02_read_spaceranger_summaries"
    )
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Pilot metrics have a different file structure
##
## Adapted from
## https://github.com/LieberInstitute/HumanPilot/blob/89e9002790a8b78c8c7ce06f5331809626386fd5/Analysis/check_metrics.R#L3-L8
##
# pilot_sample_dirs <-
#     dir(
#         "/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/10X",
#         pattern = "^151",
#         full.names = TRUE
#     )
# pilot_metric_csvs <- sapply(pilot_sample_dirs, function(x) {
#     list.files(x,
#         pattern = "metrics_summary_csv.csv",
#         full = TRUE
#     )
# })


# List all SpaceRanger output directories
# spaceranger_dirs <-
#     c(
#         "Visium_HB" = "/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Visium_Hb/processed-data/01_spaceranger",
#         "Visium_IF_AD" = "/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/spaceranger",
#         "spatialDLPFC" = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rerun_spaceranger",
#         "spatailDLPFC_IF" = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/01_spaceranger_IF",
#         "spatiall_DG_lifespan" = "/dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/processed-data/01_spaceranger",
#         "spatial_NAc" = "/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/01_spaceranger",
#         "spatial_hpc/spaceranger_novaseq" = "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/01_spaceranger/spaceranger_novaseq",
#         "spatial_hpc/spaceranger_2022-04-12_SPag033122" = "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/01_spaceranger/spaceranger_2022-04-12_SPag033122",
#         "locus-c/KMay_2021-07-09" = "/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/processed_data/spaceranger/KMay_2021-07-09",
#         "locus-c/Linda_2021-05-21" = "/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/processed_data/spaceranger/Linda_2021-05-21",
#         "locus-c/NextSeqMiSeq" = "/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/processed_data/spaceranger/NextSeqMiSeq"
#     )

SCZ_spaceranger_dirs <- "/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger"
spaceranger_metric_csvs <- lapply(SCZ_spaceranger_dirs, function(x) {
  # browser()
    message(Sys.time(), " processing directory ", x)
    sample_dirs <- dir(x, full.names = TRUE)
    potential_csvs <-
        file.path(sample_dirs, "outs", "metrics_summary.csv")
    potential_csvs[file.exists(potential_csvs)]
})

## Combine all paths to the metrics CSVs
# spaceranger_metric_csvs <-
#     c(
#         spaceranger_metric_csvs,
#         list("HumanPilot" = pilot_metric_csvs)
#     )
# lengths(spaceranger_metric_csvs)
#                                     Visium_HB
#                                             3
#                                  Visium_IF_AD
#                                            10
#                                  spatialDLPFC
#                                            30
#                               spatailDLPFC_IF
#                                             4
#                          spatiall_DG_lifespan
#                                            12
#                                   spatial_NAc
#                                            16
#               spatial_hpc/spaceranger_novaseq
#                                             8
# spatial_hpc/spaceranger_2022-04-12_SPag033122
#                                            24
#                       locus-c/KMay_2021-07-09
#                                             4
#                      locus-c/Linda_2021-05-21
#                                             3
#                          locus-c/NextSeqMiSeq
#                                             2
#                                    HumanPilot
#                                            12

## Read in all metrics files
SCZ_metrics <- lapply(spaceranger_metric_csvs, function(x) {
    do.call(rbind, lapply(x, function(y) {
        summary_info <- read.csv(y, header = TRUE)
        summary_info$summary_file <- y

        ## Find _invocation file
        invocation_file <-
            file.path(dirname(dirname(y)), "_invocation")
        if (file.exists(invocation_file)) {
            invocation <- readLines(invocation_file)
            capture_area <-
                invocation[grep("slide_serial_capture_area", invocation)]
            summary_info$slide_serial_capture_area <-
                gsub("\".*", "", gsub("^.+ \"", "", capture_area))
        } else {
            summary_info$slide_serial_capture_area <- NA
        }

        return(summary_info)
    }))
}) [[1]]

## Deal with the HumanPilot metrics separately
# pilot_metrics <- spaceranger_metrics$HumanPilot
# pilot_metrics$Sample.ID <- basename(rownames(pilot_metrics))
# rownames(pilot_metrics) <- NULL
## Adapted code from
## https://github.com/LieberInstitute/spatialDLPFC/blob/e44ee91370bdd8bc23928fa17aa641538510a322/analysis/04_sample_metrics.R

## Make it usable in R
# for (i in c(
#     "Estimated.Number.of.Spots",
#     "Mean.Reads.per.Spot",
#     "Median.Genes.per.Spot",
#     "Number.of.Reads",
#     "Total.Genes.Detected",
#     "Median.UMI.Counts.per.Spot"
# )) {
#     message("Converting ", i, " to a number")
#     pilot_metrics[[i]] <-
#         as.numeric(gsub(",", "", pilot_metrics[[i]]))
# }
# for (i in c(
#     "Valid.Barcodes",
#     "Sequencing.Saturation",
#     "Q30.Bases.in.Barcode",
#     "Q30.Bases.in.RNA.Read",
#     "Q30.Bases.in.Sample.Index",
#     "Q30.Bases.in.UMI",
#     "Reads.Mapped.to.Genome",
#     "Reads.Mapped.Confidently.to.Genome",
#     "Reads.Mapped.Confidently.to.Intergenic.Regions",
#     "Reads.Mapped.Confidently.to.Intronic.Regions",
#     "Reads.Mapped.Confidently.to.Exonic.Regions",
#     "Reads.Mapped.Confidently.to.Transcriptome",
#     "Reads.Mapped.Antisense.to.Gene",
#     "Fraction.Reads.in.Spots"
# )) {
#     message("Removing the % symbol and converting ", i, " to a number")
#     pilot_metrics[[i]] <-
#         as.numeric(gsub("%", "", pilot_metrics[[i]]))
# }
# for (i in which(
#     !colnames(pilot_metrics) %in% c(
#         "Sample.ID",
#         "Estimated.Number.of.Spots",
#         "Mean.Reads.per.Spot",
#         "Median.Genes.per.Spot",
#         "Number.of.Reads",
#         "Total.Genes.Detected",
#         "Median.UMI.Counts.per.Spot",
#         "summary_file",
#         "slide_serial_capture_area"
#     )
# )) {
#     message("Converting ", colnames(pilot_metrics)[i], " to a proportion")
#     pilot_metrics[[i]] <- pilot_metrics[[i]] / 100
# }


## All studies except the HumanPilot one have the same metric names
# colnames_nonpilot <-
#     lapply(spaceranger_metrics[-length(spaceranger_metrics)], colnames)
# # table(unlist(colnames_nonpilot))
# new_colnames <- unique(unlist(colnames_nonpilot))
# colnames(pilot_metrics)

## Re-name some old metric names
# colnames(pilot_metrics)[colnames(pilot_metrics) == "Estimated.Number.of.Spots"] <-
#     "Number.of.Spots.Under.Tissue"
# colnames(pilot_metrics)[colnames(pilot_metrics) == "Fraction.Reads.in.Spots"] <-
#     "Fraction.of.Spots.Under.Tissue"
## Drop and old metric we can't match and doens't seem to be "Valid.UMIs"
# pilot_metrics$Q30.Bases.in.Sample.Index <- NULL

## Compare old vs new metrics
# setdiff(new_colnames, colnames(pilot_metrics))
# setdiff(colnames(pilot_metrics), new_colnames)


## Combine new metrics
# all_spaceranger_metrics <-
#     do.call(rbind, spaceranger_metrics)
# all_spaceranger_metrics <-
#     dplyr::full_join(all_spaceranger_metrics, pilot_metrics)


## For how many do we have the capture area?
# all_spaceranger_metrics$study_name <-
#     rep(
#         names(spaceranger_metric_csvs),
#         lengths(spaceranger_metric_csvs)
#     )
# with(all_spaceranger_metrics, addmargins(table(
#     study_name,
#     "Has capture area info?" = !is.na(slide_serial_capture_area)
# )))
#                                                Has capture area info?
# study_name                                      FALSE TRUE Sum
#   HumanPilot                                       12    0  12
#   locus-c/KMay_2021-07-09                           0    4   4
#   locus-c/NextSeqMiSeq                              0    2   2
#   spatailDLPFC_IF                                   0    4   4
#   spatial_hpc/spaceranger_2022-04-12_SPag033122     0   24  24
#   spatial_hpc/spaceranger_novaseq                   0    8   8
#   spatial_NAc                                       0   16  16
#   spatialDLPFC                                      0   30  30
#   spatiall_DG_lifespan                              0   12  12
#   Visium_HB                                         0    3   3
#   Visium_IF_AD                                      0   10  10
#   Sum                                              12  113 125

## Export for later use
write.csv(
  SCZ_metrics,
    file.path(dir_rdata, "sequencing_metrics.csv"),
    row.names = FALSE
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
