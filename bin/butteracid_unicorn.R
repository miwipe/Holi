#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(scales)
library(data.table)
library(stringr)
library(rentrez)
library(taxonomizr)
library(tools)

source("postprocessing_functions.R")

#####################################################
# Load configuration from Bash-style config file
#####################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript postprocess_unicorn_metaDMG.R CONFIG_FILE", call. = FALSE)
}
config_file <- args[1]

read_config <- function(file) {
  lines <- readLines(file)
  lines <- lines[!grepl("^#", lines)]
  lines <- lines[nzchar(lines)]
  conf <- list()
  for (line in lines) {
    if (grepl("=", line)) {
      key <- sub("=.*", "", line)
      val <- sub("^[^=]*=", "", line)
      val <- gsub('^"|"$', "", val)
      conf[[key]] <- val
    }
  }
  return(conf)
}

config <- read_config(config_file)

#####################################################
# Define variables from config
#####################################################
LOG_FILE    <- config$LOG_FILE
TAX_PATH_NCBI <- config$TAX_PATH_NCBI
EUK_DIR     <- config$EUK_DIR
SAMPLE_LIST <- config$LIBRARY_LIST


#####################################################
# Logging helper
#####################################################
log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), file = LOG_FILE, append = TRUE)
  message(msg)
}

log_msg("Starting R postprocessing script...")
log_msg(paste("Using input directory:", EUK_DIR))

#####################################################
# Load taxonomy
#####################################################
sqlite_file <- "nameNode.sqlite"

if (!file.exists(sqlite_file)) {
  names_file <- file.path(TAX_PATH_NCBI, "taxdump", "names.dmp")
  nodes_file <- file.path(TAX_PATH_NCBI, "taxdump", "nodes.dmp")
  
  if (file.exists(names_file) && file.exists(nodes_file)) {
    log_msg("Importing taxonomy into SQLite...")
    read.names.sql(names_file, sqlFile = sqlite_file, overwrite = TRUE)
    read.nodes.sql(nodes_file, sqlFile = sqlite_file, overwrite = TRUE)
    log_msg("Taxonomy import complete.")
  } else {
    log_msg("Warning: taxonomy files missing!")
  }
} else {
  log_msg("SQLite taxonomy already exists. Skipping import.")
}

#####################################################
# Load sample list
#####################################################
sample_list_file <- SAMPLE_LIST
if (!file.exists(sample_list_file)) stop(paste("Sample list not found at", sample_list_file))
sample_list <- read.table(sample_list_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
lib_ids <- sample_list[[1]]
log_msg(paste("Loaded sample list with", length(lib_ids), "library IDs"))

#####################################################
# Helper functions
#####################################################
read_with_id <- function(file) {
  if (!file.exists(file) || file.info(file)$size == 0) {
    message(paste("Skipping missing or empty file:", file))
    return(NULL)
  }
  lib_id <- tools::file_path_sans_ext(basename(file)) %>% sub("\\..*$", "", .)
  fread(file, colClasses = "character") |> as.data.frame() |> mutate(library_id = lib_id)
}

# Read TSV if exists, else process raw files
read_or_process <- function(output_file, path = NULL, pattern = NULL) {
  if (file.exists(output_file)) {
    log_msg(paste("Reading existing file:", output_file))
    return(fread(output_file, sep = "\t", data.table = FALSE))
  } else {
    if (is.null(path) || is.null(pattern)) stop("Need path and pattern to process raw files")
    files <- list.files(path = path, pattern = pattern, full.names = TRUE, recursive = TRUE)
    if (length(files) == 0) stop(paste("No files found for pattern:", pattern))
    df <- map_dfr(files, read_with_id)
    fwrite(df, output_file, sep = "\t", quote = FALSE)
    log_msg(paste("Processed and saved:", output_file))
    return(df)
  }
}

OUTPUT_DIR <- "processed_data"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

#####################################################
# Load or process stats, LCA, acc2tax
#####################################################
stats_all_df <- read_or_process(file.path(OUTPUT_DIR, "stats_all.tsv"), EUK_DIR, ".comp\\.filtered\\.unicorn\\.refstats$")
lca_all_df   <- read_or_process(file.path(OUTPUT_DIR, "lca_all.tsv"), EUK_DIR, ".sort\\.comp\\.filtered\\.agg\\.stat\\.gz$")
raw_lca_all_df   <- read_or_process(file.path(OUTPUT_DIR, "raw_lca_all.tsv"), EUK_DIR, ".sort\\.comp\\.filtered\\.lca\\.gz$")
acc2tax_all_df <- read_or_process(file.path(OUTPUT_DIR, "acc2tax_all.tsv"), EUK_DIR, ".acc2tax$")



log_msg(paste("Stats rows:", nrow(stats_all_df)))
log_msg(paste("LCA rows:", nrow(lca_all_df)))
log_msg(paste("raw LCA rows:", nrow(raw_lca_all_df)))
log_msg(paste("Acc2Tax rows:", nrow(acc2tax_all_df)))

#####################################################
# Process metaDMG
#####################################################
metaDMG_file <- file.path(OUTPUT_DIR, "holi.tsv")
if (file.exists(metaDMG_file)) {
  metaDMG_data <- fread(metaDMG_file, sep = "\t", data.table = FALSE)
  log_msg(paste("Read existing metaDMG data:", metaDMG_file))
} else {
  metaDMG_data <- filter_metadmg(lca_all_df, unique(lca_all_df$library_id))
  fwrite(metaDMG_data, metaDMG_file, sep = "\t", quote = FALSE)
  log_msg(paste("Processed and saved metaDMG data:", metaDMG_file))
}

#####################################################
# Standardize column names
#####################################################
stats_all_df <- stats_all_df |> rename(accession = Id)
acc2tax_all_df <- acc2tax_all_df |> rename(accession = BAM_Reference, taxid = TaxID)

accs <- stats_all_df |> select(accession) |> distinct() |> pull(accession)
tax_ids <- acc2tax_all_df |> filter(accession %in% accs)
tax_ids_lst <- tax_ids |> select(taxid) |> distinct() |> pull(taxid)

stats_joined <- stats_all_df |>
  inner_join(acc2tax_all_df, by = c("accession", "library_id"))

tax_data <- getTaxonomy(
  tax_ids_lst,
  sqlite_file,
  desiredTaxa = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
) |> as_tibble() |> mutate(taxid = tax_ids_lst)

#####################################################
# Merge stats, taxonomy, and metaDMG
#####################################################
merged_file <- file.path(OUTPUT_DIR, "merged_stats_lca.tsv")
if (file.exists(merged_file)) {
  bf_md_data <- fread(merged_file, sep = "\t", data.table = FALSE)
  log_msg(paste("Read existing merged stats + metaDMG:", merged_file))
} else {
  merged_data <- stats_joined |> inner_join(tax_data, by = "taxid")
  bf_md_data <- merged_data |> inner_join(metaDMG_data, by = c("library_id" = "library_id", "taxid" = "tax_name")) |>
    rename(nreads = n_reads.x)
  fwrite(bf_md_data, merged_file, sep = "\t", quote = FALSE)
  log_msg(paste("Merged stats + metaDMG written:", merged_file))
}

#####################################################
# Aggregate data
#####################################################
agg_file <- file.path(OUTPUT_DIR, "genus_level.tsv")
if (file.exists(agg_file)) {
  bf_md_agg_data <- fread(agg_file, sep = "\t", data.table = FALSE)
  log_msg(paste("Read existing aggregated data:", agg_file))
} else {
  bf_md_agg_data <- aggregate_taxonomic_data_unicorn(bf_md_data)

}


log_msg("Finished R postprocessing for unicorn + metaDMG analysis.")
