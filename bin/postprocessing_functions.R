#!/usr/bin/env Rscript

### LIBRARIES LOAD
library(DescTools)

library(DescTools)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
library(remotes)

# Install perk from GitHub (repo contains perk/DESCRIPTION)
if (!requireNamespace("perk", quietly = TRUE)) {
  remotes::install_github("hyu-ub/perk", subdir = "perk", upgrade = "never")
}

library(perk)
library(tidyverse)
library(furrr)

####
dmg_fwd_CCC <- function(df, smp, ci = "z-transform", nperm = NULL, nproc = 1) {
  library(dplyr)
  library(tidyr)
  library(furrr)

  df <- df %>% filter(library_id %in% smp)

  # Make a list of per-library data frames (names are library_id)
  chunks <- split(df, df$library_id)

  plan(multisession, workers = nproc)

  results <- future_map_dfr(chunks, function(dX) {
    if (is.null(dX) || nrow(dX) == 0) return(NULL)

    X <- dX$library_id[1]  # library_id for labeling output

    dist_data <- dX %>%
      select(tax_name, library_id, starts_with("fwf")) %>%
      pivot_longer(-c(tax_name, library_id), names_to = "type", values_to = "f_fwd") %>%
      mutate(
        x = as.numeric(gsub("fwf", "", type)),
        f_fwd = suppressWarnings(as.numeric(f_fwd))
      ) %>%
      select(-type)

    dist_fit <- dX %>%
      select(tax_name, library_id, matches("^fwdx\\d+")) %>%
      pivot_longer(-c(tax_name, library_id), names_to = "type", values_to = "Dx_fwd") %>%
      mutate(
        x = as.numeric(gsub("fwdx", "", type)),
        Dx_fwd = suppressWarnings(as.numeric(Dx_fwd))
      ) %>%
      select(-type)

    data <- dist_data %>%
      inner_join(dist_fit, by = dplyr::join_by(tax_name, library_id, x)) %>%
      group_by(tax_name) %>%
      arrange(x, .by_group = TRUE)

    fits <- data %>%
      do({
        ccc <- DescTools::CCC(.$f_fwd, .$Dx_fwd, ci = ci)

        if (!is.null(nperm)) {
          ccc_perm <- perk_test(.$f_fwd, .$Dx_fwd, B = nperm, method = "ccc", alternative = "greater")
          tibble(
            rho_c = ccc$rho.c$est,
            rho_lwr_ci = ccc$rho.c$lwr.ci,
            rho_upr_ci = ccc$rho.c$upr.ci,
            C_b = ccc$C.b,
            l_shift = ccc$l.shift,
            s_shift = ccc$s.shift,
            library_id = X,
            rho_c_perm = ccc_perm$estimate,
            rho_c_perm_pval = ccc_perm$p.value
          )
        } else {
          tibble(
            rho_c = ccc$rho.c$est,
            rho_lwr_ci = ccc$rho.c$lwr.ci,
            rho_upr_ci = ccc$rho.c$upr.ci,
            C_b = ccc$C.b,
            l_shift = ccc$l.shift,
            s_shift = ccc$s.shift,
            library_id = X
          )
        }
      })

    fits
  }, .progress = TRUE)

  plan(sequential)
  results
}

### metaDMG filter function - variables can be changed here we filter at rank species, only Euks and split the euks into plants and animals
### we then estimate of a damage fit is good or bad estimated with the CCC model expected damage vs estimated damage fit 
filter_metadmg <- function(df, samples){
  holi_data_sp_euk <- df |>
    #filter(rank == "species") |> # 	NEEDS TO BE CONSISTNAT WITH METADMG AGG FILE 
    filter(grepl("Eukaryota", taxa_path)) |> # must be a Eukaryote 
    mutate(PlantAnimal = case_when( # define plant/animal based on taxpath 
      grepl("Viridiplantae", taxa_path) ~ "plant",
      grepl("Metazoa", taxa_path) ~ "animal",
    )) |>
    rename(tax_name = taxid, n_reads = nreads) # rename these for ease later
  
  
  # get the damage fits using CCC
  dat_all <- dmg_fwd_CCC(holi_data_sp_euk, samples, ci = "asymptotic", nperm = 100, nproc = 54)
  
  # define good and bad hits
  # bad if confidence interval is above 1 or below 0 (despite other params)
  # good if rho_c > 0.85 and C_b > 0.9 and rho_c p-value < 0.1
  dat_filt <- inner_join(dat_all, holi_data_sp_euk) |>
    mutate(fit = ifelse(rho_c >= 0.85 & C_b > 0.9 & round(rho_c_perm_pval, 3) < 0.1 & !is.na(rho_c), "good", "bad")) |>
    mutate(fit = ifelse(q_CI_h >= 1 | c_CI_l <= 0, "bad", fit))
}


penalized_weighted_median <- function(values, weights, ratios) {
  # Combine values, weights, and ratios into a single dataframe
  data <- data.frame(values, weights, ratios)
  
  # Remove rows with NA in any column
  data <- na.omit(data)
  
  # Create a new weighting scheme that penalizes the values
  # Here we multiply each value by its square and by the number of reads
  data$penalized_weights <- data$ratios^2 * data$weights
  
  # Order by values
  data <- data[order(data$values), ]
  
  # Calculate the cumulative penalized weights
  data$cum_weights <- cumsum(data$penalized_weights)
  
  # Total sum of penalized weights
  total_weights <- sum(data$penalized_weights)
  
  # Find the median based on the cumulative penalized weights
  median_index <- min(which(data$cum_weights >= total_weights / 2))
  return(data$values[median_index])
}


aggregate_taxonomic_data <- function(
    merged_data,
    species_tsv = "processed_data/species_level.tsv",
    genus_tsv   = "processed_data/genus_level.tsv"
) {
  # Filter strictly to species-level reads (if rank column exists)
  if ("rank" %in% colnames(merged_data)) {
    merged_data <- merged_data |> dplyr::filter(rank == "species")
    message("Species-level reads (head):")
    print(utils::head(merged_data))
  } else {
    warning("No 'rank' column detected; assuming species-level reads are correctly formatted.")
  }
  
  # ---- SPECIES-LEVEL AGGREGATION ----
  species_group_levels <- c("library_id", "kingdom", "phylum", "class",
                            "order", "family", "genus", "species", "fit") |>
    intersect(colnames(merged_data))
  
  species_aggregated <- merged_data |>
    dplyr::distinct(library_id, taxid, .keep_all = TRUE) |> 
    dplyr::group_by(dplyr::across(dplyr::all_of(species_group_levels))) |> 
    dplyr::summarise(
      n = dplyr::n_distinct(taxid),
      median_A = median(A, na.rm = TRUE),
      penalized_weighted_median_A = penalized_weighted_median(A, nreads, A),
	  median_c_b = median(c_b, na.rm = TRUE),
	  median_A_CI_l = median(A_CI_l, na.rm = TRUE),
	  median_A_CI_h = median(A_CI_h, na.rm = TRUE),
      mean_read_ani_median = mean(read_ani_median, na.rm = TRUE),
      mean_read_ani_std = mean(read_ani_std, na.rm = TRUE),
      median_reference_length = median(reference_length, na.rm = TRUE),
      sum_reference_length = sum(reference_length, na.rm = TRUE),
      mean_breadth_exp_ratio = mean(breadth_exp_ratio, na.rm = TRUE),
      median_breadth_exp_ratio = median(breadth_exp_ratio, na.rm = TRUE),
      penalized_weighted_median_entropy = penalized_weighted_median(norm_entropy, nreads, norm_entropy),
      penalized_weighted_median_gini = penalized_weighted_median(norm_gini, nreads, norm_gini),
      penalized_weighted_median_breadth_exp_ratio = penalized_weighted_median(breadth_exp_ratio, nreads, breadth_exp_ratio),
      median_entropy = median(norm_entropy, na.rm = TRUE),
      mean_entropy = mean(norm_entropy, na.rm = TRUE),
      median_gini = median(norm_gini, na.rm = TRUE),
      mean_gini = mean(norm_gini, na.rm = TRUE),
      median_n_reads = median(nreads, na.rm = TRUE),
      mean_n_reads = mean(nreads, na.rm = TRUE),
      total_n_reads = sum(nreads, na.rm = TRUE),
      n_reads_tad = sum(n_reads_tad, na.rm = TRUE),
      median_tax_abund_tad = median(tax_abund_tad, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Save species-level results
  readr::write_tsv(species_aggregated, species_tsv)
  message("Exported species-level aggregation to: ", species_tsv)
  
  # ---- GENUS-LEVEL AGGREGATION ----
  genus_group_levels <- c("library_id", "kingdom", "phylum", "class",
                          "order", "family", "genus", "fit") |>
    intersect(colnames(species_aggregated))
  
  genus_aggregated <- species_aggregated |>
    dplyr::group_by(dplyr::across(dplyr::all_of(genus_group_levels))) |> 
    dplyr::summarise(
      n_species = dplyr::n_distinct(species, na.rm = TRUE),
      median_A = median(median_A, na.rm = TRUE),
      penalized_weighted_median_A = penalized_weighted_median(median_A, total_n_reads, median_A),
      median_c_b = median(median_c_b, na.rm = TRUE),
	  median_A_CI_l = median(median_A_CI_l, na.rm = TRUE),
	  median_A_CI_h = median(median_A_CI_h, na.rm = TRUE),
      mean_read_ani_median = mean(mean_read_ani_median, na.rm = TRUE),
      mean_read_ani_std = mean(mean_read_ani_std, na.rm = TRUE),
      median_reference_length = median(median_reference_length, na.rm = TRUE),
      sum_reference_length = sum(sum_reference_length, na.rm = TRUE),
      mean_breadth_exp_ratio = mean(mean_breadth_exp_ratio, na.rm = TRUE),
      median_breadth_exp_ratio = median(median_breadth_exp_ratio, na.rm = TRUE),
      penalized_weighted_median_entropy = penalized_weighted_median(median_entropy, total_n_reads, median_entropy),
      penalized_weighted_median_gini = penalized_weighted_median(median_gini, total_n_reads, median_gini),
      penalized_weighted_median_breadth_exp_ratio = penalized_weighted_median(median_breadth_exp_ratio, total_n_reads, median_breadth_exp_ratio),
      median_entropy = median(median_entropy, na.rm = TRUE),
      mean_entropy = mean(mean_entropy, na.rm = TRUE),
      median_gini = median(median_gini, na.rm = TRUE),
      mean_gini = mean(mean_gini, na.rm = TRUE),
      median_n_reads = median(median_n_reads, na.rm = TRUE),
      mean_n_reads = mean(mean_n_reads, na.rm = TRUE),
      total_n_reads = sum(total_n_reads, na.rm = TRUE),
      n_reads_tad = sum(n_reads_tad, na.rm = TRUE),
      median_tax_abund_tad = median(median_tax_abund_tad, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Save genus-level results
  readr::write_tsv(genus_aggregated, genus_tsv)
  message("Exported genus-level aggregation to: ", genus_tsv)
  
  return(list(species_level = species_aggregated, genus_level = genus_aggregated))
}

aggregate_taxonomic_data_unicorn <- function(
    merged_data,
    species_tsv = "processed_data/species_level.tsv",
    genus_tsv   = "processed_data/genus_level.tsv"
) {
  # Filter strictly to species-level reads (if rank column exists)
  if ("rank" %in% colnames(merged_data)) {
    merged_data <- merged_data |> dplyr::filter(rank == "species")
    message("Species-level reads (head):")
    print(utils::head(merged_data))
  } else {
    warning("No 'rank' column detected; assuming species-level reads are correctly formatted.")
  }
  
  # ---- SPECIES-LEVEL AGGREGATION ----
  species_group_levels <- c("library_id", "kingdom", "phylum", "class",
                            "order", "family", "genus", "species", "fit") |>
    intersect(colnames(merged_data))
  
  species_aggregated <- merged_data |>
    dplyr::distinct(library_id, taxid, .keep_all = TRUE) |> 
    dplyr::group_by(dplyr::across(dplyr::all_of(species_group_levels))) |> 
    dplyr::summarise(
      n = dplyr::n_distinct(taxid),
      median_A = median(A, na.rm = TRUE),
      penalized_weighted_median_A = penalized_weighted_median(A, nreads, A),
      median_c_b = median(c_b, na.rm = TRUE),
	  median_A_CI_l = median(A_CI_l, na.rm = TRUE),
	  median_A_CI_h = median(A_CI_h, na.rm = TRUE),
      mean_read_ani_median = mean(md_alnani, na.rm = TRUE),
      mean_read_ani_std = mean(std_alnani, na.rm = TRUE),
      median_reference_length = median(Length, na.rm = TRUE),
      sum_reference_length = sum(Length, na.rm = TRUE),
      mean_breadth_exp_ratio = mean(breath_ratio, na.rm = TRUE),
      median_breadth_exp_ratio = median(breath_ratio, na.rm = TRUE),
      penalized_weighted_median_entropy = penalized_weighted_median(n_entropy, nreads, n_entropy),
      penalized_weighted_median_gini = penalized_weighted_median(n_gini, nreads, n_gini),
      penalized_weighted_median_breadth_exp_ratio = penalized_weighted_median(breath_ratio, nreads, breath_ratio),
      median_entropy = median(n_entropy, na.rm = TRUE),
      mean_entropy = mean(n_entropy, na.rm = TRUE),
      median_gini = median(n_gini, na.rm = TRUE),
      mean_gini = mean(n_gini, na.rm = TRUE),
      median_n_reads = median(nreads, na.rm = TRUE),
      mean_n_reads = mean(nreads, na.rm = TRUE),
      total_n_reads = sum(nreads, na.rm = TRUE),
      n_reads_tad = sum(tad80, na.rm = TRUE),
      median_tax_abund_tad = median(tad80, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Save species-level results
  readr::write_tsv(species_aggregated, species_tsv)
  message("Exported species-level aggregation to: ", species_tsv)
  
  # ---- GENUS-LEVEL AGGREGATION ----
  genus_group_levels <- c("library_id", "kingdom", "phylum", "class",
                          "order", "family", "genus", "fit") |>
    intersect(colnames(species_aggregated))
  
  genus_aggregated <- species_aggregated |>
    dplyr::group_by(dplyr::across(dplyr::all_of(genus_group_levels))) |> 
    dplyr::summarise(
      n_species = dplyr::n_distinct(species, na.rm = TRUE),
      median_A = median(median_A, na.rm = TRUE),
	  median_A_CI_l = median(median_A_CI_l, na.rm = TRUE),
	  median_A_CI_h = median(median_A_CI_h, na.rm = TRUE),
      penalized_weighted_median_A = penalized_weighted_median(median_A, total_n_reads, median_A),
      median_c_b = median(median_c_b, na.rm = TRUE),
      mean_read_ani_median = mean(mean_read_ani_median, na.rm = TRUE),
      mean_read_ani_std = mean(mean_read_ani_std, na.rm = TRUE),
      median_reference_length = median(median_reference_length, na.rm = TRUE),
      sum_reference_length = sum(sum_reference_length, na.rm = TRUE),
      mean_breadth_exp_ratio = mean(mean_breadth_exp_ratio, na.rm = TRUE),
      median_breadth_exp_ratio = median(median_breadth_exp_ratio, na.rm = TRUE),
      penalized_weighted_median_entropy = penalized_weighted_median(median_entropy, total_n_reads, median_entropy),
      penalized_weighted_median_gini = penalized_weighted_median(median_gini, total_n_reads, median_gini),
      penalized_weighted_median_breadth_exp_ratio = penalized_weighted_median(median_breadth_exp_ratio, total_n_reads, median_breadth_exp_ratio),
      median_entropy = median(median_entropy, na.rm = TRUE),
      mean_entropy = mean(mean_entropy, na.rm = TRUE),
      median_gini = median(median_gini, na.rm = TRUE),
      mean_gini = mean(mean_gini, na.rm = TRUE),
      median_n_reads = median(median_n_reads, na.rm = TRUE),
      mean_n_reads = mean(mean_n_reads, na.rm = TRUE),
      total_n_reads = sum(total_n_reads, na.rm = TRUE),
      n_reads_tad = sum(n_reads_tad, na.rm = TRUE),
      median_tax_abund_tad = median(median_tax_abund_tad, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Save genus-level results
  readr::write_tsv(genus_aggregated, genus_tsv)
  message("Exported genus-level aggregation to: ", genus_tsv)
  
  return(list(species_level = species_aggregated, genus_level = genus_aggregated))
}
