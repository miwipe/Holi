#!/usr/bin/env Rscript

suppressPackageStartupMessages({
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
  library(future)
  library(furrr)
  library(tibble)
  library(DescTools)
})

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
  conf
}

config <- read_config(config_file)

#####################################################
# Define variables from config
#####################################################
LOG_FILE      <- config$LOG_FILE
TAX_PATH_NCBI <- config$TAX_PATH_NCBI
EUK_DIR       <- config$EUK_DIR
SAMPLE_LIST   <- config$LIBRARY_LIST

#####################################################
# Logging helper
#####################################################
log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), file = LOG_FILE, append = TRUE)
  message(msg)
}

#####################################################
# Cache helpers
#####################################################
cache_read <- function(path, read_fun, ...) {
  if (file.exists(path) && file.info(path)$size > 0) {
    log_msg(paste("Cache hit:", path))
    return(read_fun(path, ...))
  }
  NULL
}

cache_write <- function(df, path, write_fun, ...) {
  write_fun(df, path, ...)
  log_msg(paste("Wrote:", path))
  df
}

#####################################################
# Start
#####################################################
log_msg("Starting R postprocessing script...")
log_msg(paste("Using input directory:", EUK_DIR))

OUTPUT_DIR <- "processed_data_ss"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Final targets: if both exist, stop immediately (fastest rerun)
merged_species_file <- file.path(OUTPUT_DIR, "merged_stats_metaDMG_species_level.tsv")
merged_genus_file   <- file.path(OUTPUT_DIR, "merged_stats_metaDMG_genus_level.tsv")

if (file.exists(merged_species_file) && file.exists(merged_genus_file)) {
  log_msg("Merged outputs already exist; nothing to do. Exiting.")
  quit(save = "no", status = 0)
}

#####################################################
# Load taxonomy DB (only builds sqlite once)
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
# IO helpers (unchanged)
#####################################################
read_with_id <- function(file) {
  if (!file.exists(file) || file.info(file)$size == 0) {
    message(paste("Skipping missing or empty file:", file))
    return(NULL)
  }
  lib_id <- tools::file_path_sans_ext(basename(file)) %>% sub("\\..*$", "", .)
  fread(file, colClasses = "character") |>
    as.data.frame() |>
    mutate(library_id = lib_id)
}

read_or_process <- function(output_file, path = NULL, pattern = NULL) {
  if (file.exists(output_file) && file.info(output_file)$size > 0) {
    log_msg(paste("Reading existing file:", output_file))
    return(fread(output_file, sep = "\t", data.table = FALSE))
  } else {
    if (is.null(path) || is.null(pattern)) stop("Need path and pattern to process raw files")
    files <- list.files(path = path, pattern = pattern, full.names = TRUE, recursive = TRUE)
    if (length(files) == 0) stop(paste("No files found for pattern:", pattern))
    df <- map_dfr(files, read_with_id)
    fwrite(df, output_file, sep = "\t", quote = FALSE)
    log_msg(paste("Processed and saved:", output_file))
    df
  }
}

#####################################################
# metaDMG caches: avoid reading full lca_all unless needed
#####################################################
metaDMG_species_file <- file.path(OUTPUT_DIR, "metaDMG_species_level.tsv")
metaDMG_genus_file   <- file.path(OUTPUT_DIR, "metaDMG_genus_level.tsv")

metaDMG_species <- cache_read(metaDMG_species_file, fread, sep = "\t", data.table = FALSE)
metaDMG_genus   <- cache_read(metaDMG_genus_file,   fread, sep = "\t", data.table = FALSE)

if (is.null(metaDMG_species) || is.null(metaDMG_genus)) {
  log_msg("metaDMG rank caches missing; loading/processing lca_all.tsv ...")
  lca_all_df <- read_or_process(
    file.path(OUTPUT_DIR, "lca_all.tsv"),
    EUK_DIR,
    ".sort\\.comp\\.filtered\\.agg\\.stat\\.gz$"
  )
  log_msg(paste("LCA rows:", nrow(lca_all_df)))
  
  log_msg("Computing CCC damage-fit statistics for Eukaryota species rows (no fit classification)...")
  lca_all_df <- add_metadmg_ccc_stats(
    df      = lca_all_df,
    samples = unique(lca_all_df$library_id),
    ci      = "asymptotic",
    nperm   = 100,
    nproc   = 14
  )
  log_msg("CCC statistics joined onto lca_all_df.")
  

  if (is.null(metaDMG_species)) {
    metaDMG_species <- lca_all_df |>
      filter(rank == "species") |>
      mutate(
        library_id = as.character(library_id),
        taxid      = as.character(taxid)
      )|>
	  rename(total_n_reads_from_mD = nreads)
	  
    metaDMG_species <- cache_write(metaDMG_species, metaDMG_species_file, fwrite, sep = "\t", quote = FALSE)
  }

  if (is.null(metaDMG_genus)) {
	  metaDMG_genus <- lca_all_df |>
	    filter(rank == "genus") |>
	    mutate(
	      library_id = as.character(library_id),
	      taxid      = as.character(taxid)  # keep genus taxid!
	    ) |>
	    rename(total_n_reads_from_mD = nreads)
    metaDMG_genus <- cache_write(metaDMG_genus, metaDMG_genus_file, fwrite, sep = "\t", quote = FALSE)
  }
} else {
  log_msg(paste("Loaded metaDMG species rows:", nrow(metaDMG_species)))
  log_msg(paste("Loaded metaDMG genus rows:", nrow(metaDMG_genus)))
}

#####################################################
# Unicorn caches: avoid reading stats_all/acc2tax unless needed
#####################################################
stats_species_file <- file.path(OUTPUT_DIR, "unicorn_stats_species_level.tsv")
stats_genus_file   <- file.path(OUTPUT_DIR, "unicorn_stats_genus_level.tsv")

stats_species <- cache_read(stats_species_file, fread, sep = "\t", data.table = FALSE)
stats_genus   <- cache_read(stats_genus_file,   fread, sep = "\t", data.table = FALSE)

if (is.null(stats_species)) {
  log_msg("Species stats cache missing; loading/processing stats_all + acc2tax ...")

  stats_all_df <- read_or_process(
    file.path(OUTPUT_DIR, "stats_all.tsv"),
    EUK_DIR,
    ".comp\\.filtered\\.unicorn\\.refstats$"
  )
  acc2tax_all_df <- read_or_process(
    file.path(OUTPUT_DIR, "acc2tax_all.tsv"),
    EUK_DIR,
    ".acc2tax$"
  )

  log_msg(paste("Stats rows:", nrow(stats_all_df)))
  log_msg(paste("Acc2Tax rows:", nrow(acc2tax_all_df)))

  # Standardize column names
  stats_all_df   <- stats_all_df   |> rename(accession = Id, n_reads_from_unicorn = n_reads)
  acc2tax_all_df <- acc2tax_all_df |> rename(accession = BAM_Reference, taxid = TaxID)

  stats_joined <- stats_all_df |>
    inner_join(acc2tax_all_df, by = c("accession", "library_id"))

  stats_species <- summarise_stats_joined_all(stats_joined)
  stats_species <- cache_write(stats_species, stats_species_file, fwrite, sep = "\t", quote = FALSE)
} else {
  log_msg(paste("Loaded unicorn species stats rows:", nrow(stats_species)))
}


#####################################################
# Taxonomy cache: include BOTH unicorn + metaDMG taxids
#####################################################
tax_data_file <- file.path(OUTPUT_DIR, "tax_data_all_taxids.tsv")
tax_data <- cache_read(tax_data_file, fread, sep = "\t", data.table = FALSE)

if (is.null(tax_data)) {
  log_msg("Taxonomy cache missing; querying SQLite taxonomy for union of taxids...")

  # metaDMG may or may not exist yet depending on cache hits
  md_sp_taxids <- if (exists("metaDMG_species") && !is.null(metaDMG_species))
    metaDMG_species |> dplyr::distinct(taxid) |> dplyr::pull(taxid) else character()

  md_ge_taxids <- if (exists("metaDMG_genus") && !is.null(metaDMG_genus))
    metaDMG_genus |> dplyr::distinct(taxid) |> dplyr::pull(taxid) else character()

  uni_taxids <- stats_species |> dplyr::distinct(taxid) |> dplyr::pull(taxid)

  tax_ids_lst <- unique(c(as.character(uni_taxids), as.character(md_sp_taxids), as.character(md_ge_taxids)))
  tax_ids_lst <- tax_ids_lst[!is.na(tax_ids_lst) & nzchar(tax_ids_lst)]

  tax_data <- getTaxonomy(
    tax_ids_lst,
    sqlite_file,
    desiredTaxa = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  ) |>
    tibble::as_tibble() |>
    dplyr::mutate(taxid = tax_ids_lst)

  tax_data <- cache_write(tax_data, tax_data_file, fwrite, sep = "\t", quote = FALSE)
} else {
  tax_data <- tibble::as_tibble(tax_data)
  log_msg(paste("Loaded taxonomy rows:", nrow(tax_data)))
}

# helper: stable join key type
tax_data <- tax_data |> dplyr::mutate(taxid = readr::parse_integer(as.character(taxid)))


#####################################################
# Add taxonomy to unicorn species stats
#####################################################
stats_species_tax <- stats_species |>
  mutate(taxid = as.integer(taxid)) |>
  inner_join(tax_data |> mutate(taxid = as.integer(taxid)), by = "taxid")

log_msg(paste("Extracted species level rows from unicorn statistics (with taxonomy):", nrow(stats_species_tax)))

#####################################################
# Add taxonomy to metaDMG species + genus (always)
#####################################################
metaDMG_species <- metaDMG_species |>
  dplyr::mutate(taxid = readr::parse_integer(as.character(taxid))) |>
  dplyr::left_join(tax_data, by = "taxid")

metaDMG_genus <- metaDMG_genus |>
  dplyr::mutate(taxid = readr::parse_integer(as.character(taxid))) |>
  dplyr::left_join(tax_data, by = "taxid") |>
  dplyr::mutate(
    genus   = as.character(.data$genus),
    family  = as.character(.data$family),
    kingdom = as.character(.data$kingdom)
  )


#####################################################
# Genus stats: compute only if missing (from species_tax)
#####################################################
if (is.null(stats_genus)) {
  log_msg("Genus stats cache missing; aggregating from species-level stats ...")
  stats_genus <- summarise_stats_joined_all(
    stats_joined = stats_species_tax,
    group_cols   = c("library_id", "genus", "family", "kingdom"),
    weight_col   = "total_reads_from_unicorn",
    exclude_cols = c("library_id", "taxid", "total_reads_from_unicorn", "n_references", "genus"),
    sum_cols     = c("reference_length")
  )
  stats_genus <- cache_write(stats_genus, stats_genus_file, fwrite, sep = "\t", quote = FALSE)
} else {
  log_msg(paste("Loaded unicorn genus stats rows:", nrow(stats_genus)))
}

#####################################################
# Coerce join keys robustly for metaDMG species
#####################################################
metaDMG_species <- metaDMG_species |>
  mutate(
    library_id = as.character(library_id),
    taxid = readr::parse_integer(as.character(taxid))
  )

bad_taxid <- metaDMG_species |> filter(is.na(taxid))
if (nrow(bad_taxid) > 0) {
  log_msg(paste("Warning:", nrow(bad_taxid), "metaDMG species rows have non-integer taxid; they won't join."))
}

#####################################################
# Merge unicorn + metaDMG
#####################################################
bf_md_species <- stats_species_tax |>
  mutate(
    library_id = as.character(library_id),
    taxid = as.integer(taxid)
  ) |>
  inner_join(
    metaDMG_species |> mutate(taxid = as.integer(taxid)),
    by = c("library_id", "taxid")
  )

  bf_md_genus <- stats_genus |>
    dplyr::mutate(
      library_id = as.character(library_id),
      genus      = as.character(genus),
      family     = as.character(family),
      kingdom    = as.character(kingdom)
    ) |>
    dplyr::inner_join(
      metaDMG_genus |>
        dplyr::mutate(
          library_id = as.character(library_id),
          genus      = as.character(genus),
          family     = as.character(family),
          kingdom    = as.character(kingdom)
        ),
      by = c("library_id", "genus", "family", "kingdom")
    )


#####################################################
# Write merged outputs (cache)
#####################################################
if (!file.exists(merged_species_file) || file.info(merged_species_file)$size == 0) {
  fwrite(bf_md_species, merged_species_file, sep = "\t", quote = FALSE)
  log_msg(paste("Wrote species merged file:", merged_species_file))
} else {
  log_msg(paste("Species merged file already exists:", merged_species_file))
}

if (!file.exists(merged_genus_file) || file.info(merged_genus_file)$size == 0) {
  fwrite(bf_md_genus, merged_genus_file, sep = "\t", quote = FALSE)
  log_msg(paste("Wrote genus merged file:", merged_genus_file))
} else {
  log_msg(paste("Genus merged file already exists:", merged_genus_file))
}

log_msg("Finished R postprocessing for unicorn + metaDMG analysis.")

