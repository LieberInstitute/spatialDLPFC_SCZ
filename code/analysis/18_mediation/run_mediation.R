## command-line entry point; all project inputs are read-only.
args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
source(file.path(script_dir, "R", "core.R"))
opt <- list(stage = "audit", project_root = normalizePath(file.path(script_dir, "../../..")),
            config = file.path(script_dir, "config.R"), screens = "all", workers = "1", outdir = NULL, plots = NULL)
if ("--help" %in% args) {
  cat("Rscript run_mediation.R --stage audit|historical|screen|sensitivity|overlap|report|primary|all\n",
      "  --project-root PATH --config PATH --screens all|comma-separated-IDs\n",
      "  --outdir PATH --plots PATH --workers N --audit-only\n",
      "Stages 'all' and 'overlap' load the large raw spot object; use the 64 GB job.\n")
  quit(status = 0)
}
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--audit-only") { opt$stage <- "audit"; i <- i + 1L; next }
  key <- gsub("-", "_", sub("^--", "", args[i]))
  assert(key %in% names(opt) && i < length(args), paste("Unknown/incomplete argument", args[i]))
  opt[[key]] <- args[i + 1L]; i <- i + 2L
}
assert(opt$stage %in% c("audit", "historical", "screen", "sensitivity", "overlap", "report", "primary", "all"), "Unknown stage")
workers <- as.integer(opt$workers); assert(is.finite(workers) && workers >= 1L, "Invalid worker count")
packages <- c("SpatialExperiment", "edgeR", "limma", "data.table", "digest", "Matrix")
for (pkg in packages) assert(requireNamespace(pkg, quietly = TRUE), paste("Missing R package", pkg))
data.table::setDTthreads(1L)
source(opt$config)
source(file.path(script_dir, "R", "stages.R"))
source(file.path(script_dir, "R", "overlap.R"))
source(file.path(script_dir, "R", "report.R"))
set.seed(config$seed)
all_screens <- read_tsv(file.path(script_dir, "screens.tsv"))
unique_keys(all_screens$screen_id, "Screen IDs")
screens <- all_screens
if (opt$screens != "all") {
  ids <- strsplit(opt$screens, ",", fixed = TRUE)[[1]]
  assert(all(ids %in% screens$screen_id), "Unknown requested screen")
  screens <- screens[match(ids, screens$screen_id), ]
}
root <- normalizePath(path.expand(opt$project_root))
out <- if (is.null(opt$outdir)) file.path(root, "processed-data", "18_mediation") else path.expand(opt$outdir)
plots <- if (is.null(opt$plots)) file.path(root, "plots", "18_mediation") else path.expand(opt$plots)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
files <- unique(file.path(root, c(config$pb, config$historical, config$donor_meta)))
assert(all(file.exists(files)), paste("Unavailable input:", paste(files[!file.exists(files)], collapse = ", ")))
info <- file.info(files)
fp <- data.frame(path = files, bytes = info$size, modified = as.character(info$mtime), md5 = unname(tools::md5sum(files)))
code <- c(list.files(script_dir, pattern = "[.]R$", recursive = TRUE, full.names = TRUE), file.path(script_dir, "screens.tsv"), normalizePath(opt$config))
code <- unique(code)
code_fp <- data.frame(path = code, md5 = unname(tools::md5sum(code)))
versions <- vapply(packages, function(p) as.character(packageVersion(p)), "")
signature <- hash(list(config, all_screens, fp, code_fp, versions))
env <- list(root = root, out = out, plots = plots, signature = signature)
provenance_dir <- file.path(out, "provenance", signature)
dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(fp, file.path(provenance_dir, "input_fingerprints.tsv"))
write_tsv(code_fp, file.path(provenance_dir, "code_fingerprints.tsv"))
write_tsv(all_screens, file.path(provenance_dir, "screens.tsv"))
save_atomic(config, file.path(provenance_dir, "config.rds"))
writeLines(capture.output(sessionInfo()), file.path(provenance_dir, "sessionInfo.txt"))
commit <- tryCatch(system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE), error = function(e) "unavailable")
writeLines(commit, file.path(provenance_dir, "repository_commit.txt"))
write_tsv(data.frame(time = as.character(Sys.time()), stage = opt$stage, signature = signature,
                    screens = paste(screens$screen_id, collapse = ",")),
          file.path(provenance_dir, paste0("invocation_", opt$stage, ".tsv")))
inputs <- if (opt$stage != "report") load_inputs(config, root) else NULL
stages <- if (opt$stage == "all") c("audit", "historical", "screen", "sensitivity", "overlap", "report") else if (opt$stage == "primary") c("audit", "historical", "screen", "sensitivity", "report") else opt$stage
for (stage in stages) {
  message("Stage: ", stage, " signature=", substr(signature, 1, 12))
  switch(stage,
    audit = audit_stage(inputs, screens, config, env),
    historical = historical_stage(inputs, config, env, workers),
    screen = primary_stage(inputs, screens, all_screens, config, env, workers),
    sensitivity = sensitivity_stage(inputs, screens, config, env, workers),
    overlap = overlap_stage(inputs, screens, config, env),
    report = report_stage(screens, config, env))
}
message("Completed: ", paste(stages, collapse = ", "))
