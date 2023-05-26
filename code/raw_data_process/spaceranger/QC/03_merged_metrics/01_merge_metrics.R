library("here")
library("dplyr")
library("sessioninfo")

## Rscript -e 'sgejobs::job_single("01_merge_metrics_spaceranger", create_shell = TRUE, queue = "bluejay")'

## Output dirs
dir_rdata <-
    here("processed-data", "spaceranger", "03_merged_metrics")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)


## Read in previous metric tables.
presequencing_summary_metrics <- read.csv(
    here(
        "processed-data",
        "spaceranger",
        "01_read_experimental_summaries",
        "presequencing_summary.csv"
    ),
    check.names = FALSE
)
## Combine Slide serial number and capture area
presequencing_summary_metrics$slide_serial_capture_area <-
    toupper(paste0(
        presequencing_summary_metrics[["Slide #"]],
        "-",
        presequencing_summary_metrics[["Array #"]]
    ))

## Check for any duplications. This can happen when samples are included in
## more than one summary sheet
presequencing_summary_metrics[presequencing_summary_metrics$slide_serial_capture_area %in% presequencing_summary_metrics$slide_serial_capture_area[duplicated(presequencing_summary_metrics$slide_serial_capture_area) & !is.na(presequencing_summary_metrics$slide_serial_capture_area)], ]

## Manually drop the duplicated V10U24-093-D1 from LC and HB (it's an HB sample)
presequencing_summary_metrics <-
    presequencing_summary_metrics[!(
        presequencing_summary_metrics$slide_serial_capture_area == "V10U24-093-D1" &
            basename(presequencing_summary_metrics$sheet_file) == "Visium_LC_050321_Master_HD_Round2.xlsx"
    ), ]

## All should not be NAs
stopifnot(all(
    !is.na(presequencing_summary_metrics$slide_serial_capture_area)
))
## There should be no duplications at this stage
stopifnot(!any(
    duplicated(presequencing_summary_metrics$slide_serial_capture_area)
))
## And all should have two dashes
stopifnot(!any(
    grepl(
        "[[:alpha:]]+-[[:digit:]]+-[[:alpha:]]+",
        presequencing_summary_metrics$slide_serial_capture_area
    )
))

all_spaceranger_metrics <- read.csv(
    here(
        "processed-data",
        "spaceranger",
        "02_read_spaceranger_summaries",
        "sequencing_metrics.csv"
    ),
    check.names = FALSE
)

## Check for any duplications. This can happen when samples are processed
## on spaceranger for different FASTQ file inputs
all_spaceranger_metrics[all_spaceranger_metrics$slide_serial_capture_area %in% all_spaceranger_metrics$slide_serial_capture_area[duplicated(all_spaceranger_metrics$slide_serial_capture_area) & !is.na(all_spaceranger_metrics$slide_serial_capture_area)], ]


## There should be no duplications at this stage. Note that there are 12 NAs
## due to the HumanPilot study
stopifnot(!any(
    duplicated(all_spaceranger_metrics$slide_serial_capture_area) & !is.na(all_spaceranger_metrics$slide_serial_capture_area)
))

## Merge all the data into one large table
merged_metrics <-
    dplyr::full_join(presequencing_summary_metrics, all_spaceranger_metrics, by = "slide_serial_capture_area")

## Fix 0 "Ct" values
merged_metrics$Ct[merged_metrics$Ct == 0] <- NA

## Fix some "% Coverage Array" that are stored as proportions and not percents
which_prop <- which(merged_metrics[["% Coverage Array"]] <= 1)
merged_metrics[["% Coverage Array"]][which_prop] <- merged_metrics[["% Coverage Array"]][which_prop] * 100

## Shorten name (that is, remove the spaceranger sub-directory piece for some
## studies like locus-c and spatial_hpc)
merged_metrics$study_name_short <- gsub("/.*", "", merged_metrics$study_name)

## Add the study name based on the sheet files
merged_metrics$study_name_sheet <- "DLPFC_SCZ"
## Only the HumanPilot doesn't have summary sheet files
# stopifnot(sum(is.na(merged_metrics$study_name_sheet)) == 12)
# merged_metrics$study_name_sheet[is.na(merged_metrics$study_name_sheet)] <- "HumanPilot"


## Label samples into:
## * Seq: ok
## * Seq: fail
## * Pre: fail
## * Pre: ok?

## We don't know if these are ok or not
merged_metrics$sample_status <- "Pre: ok?"
merged_metrics$sample_status[!is.na(merged_metrics$Number.of.Spots.Under.Tissue)] <- "Seq: ok"
merged_metrics$sample_status[merged_metrics$slide_serial_capture_area %in% c(
    ## Anthony's spatial_DG_lifespan samples. Note that "D" is borderline
    paste0("V11D01-387-", LETTERS[1:4], "1"),
    ## 1 DLPFC and 1 LC: pilot blocks, with MiSeq sequencing that didn't look good
    ## TODO: create a small repo and run these 2 samples on spaceranger to get
    ## the post-sequencing metrics
    "V19B23-076-A1", "V19B23-076-D1"
)] <- "Seq: fail"
merged_metrics$sample_status[merged_metrics$slide_serial_capture_area %in% c(
    ## Sang Ho's Visium IF AD samples.
    ## Tissue had morphology problems for these 2 samples
    paste0("V10A27-004-", LETTERS[2:3], "1"),
    ## spatialDLPFC samples:
    ## * V10U24-091-A1: torn
    ## * V10U24-092-A1: wrinkled
    ## * V10B01-052-C1: very wrinkled
    ## * V10B01-002-A1: Abby decided not to sequence it: coverage is 20%
    "V10U24-091-A1", "V10U24-092-A1", "V10B01-052-C1", "V10B01-002-A1"
)] <- "Pre: fail"

## Check those that are "Pre: ok?" (aka, unlabeled by other rules)
# merged_metrics[merged_metrics$sample_status == "Pre: ok?", ]

## Tissue status: note that some samples could have decent/ok pre-sequencing
## metrics but were not sequenced due to issues with the tissue itself. Like
## major wrinkles, being torn, etc. We are thinking of dropping these samples
## from any prediction algorithm since the behavior of the pre-sequencing
## metrics likely doesn't predict these failed samples.
merged_metrics$tissue_status <- NA
merged_metrics$tissue_status[merged_metrics$slide_serial_capture_area %in% c(
    ## Sang Ho's Visium IF AD samples.
    ## Tissue had morphology problems for these 2 samples
    paste0("V10A27-004-", LETTERS[2:3], "1"),
    ## spatialDLPFC samples:
    ## * V10U24-091-A1: torn
    ## * V10U24-092-A1: wrinkled
    ## * V10B01-052-C1: very wrinkled
    ## * V10B01-002-A1: Abby decided not to sequence it: coverage is 20%
    "V10U24-091-A1", "V10U24-092-A1", "V10B01-052-C1", "V10B01-002-A1"
)] <- "Fail"

## Export for later use
write.csv(merged_metrics,
    file.path(dir_rdata, "merged_metrics.csv"),
    row.names = FALSE
)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
