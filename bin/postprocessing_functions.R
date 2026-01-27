#!/bin/bash/Rscript 

### LIBRARIES LOAD
library(DescTools)

# Install remotes if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

library(remotes)

# Install perk from GitHub
#remotes::install_github("hyu-ub/perk/perk")

# Load perk and other libraries
library(perk)
library(tidyverse)
library(furrr)

####

dmg_fwd_CCC <- function(df, smp, ci = "z-transform", nperm = NULL, nproc = 1) {
  options(future.globals.maxSize = 2291289600)
  plan(multisession, workers = nproc)
  
  results <- future_map_dfr(smp, function(X) {
    dist_data <- df %>%
      filter(library_id == X) %>%
      select(tax_name, library_id, starts_with("fwf")) %>%
      pivot_longer(
        cols = -c(tax_name, library_id),
        names_to = "type",
        values_to = "f_fwd"
      ) %>%
      mutate(
        x = as.numeric(gsub("fwf", "", type)),
        f_fwd = suppressWarnings(as.numeric(f_fwd))  # <–– force numeric
      ) %>%
      select(-type)
    
    dist_fit <- df %>%
      filter(library_id == X) %>%
      select(tax_name, library_id, matches("^fwdx\\d+")) %>%
      pivot_longer(
        cols = -c(tax_name, library_id),
        names_to = "type",
        values_to = "Dx_fwd"
      ) %>%
      mutate(
        x = as.numeric(gsub("fwdx", "", type)),
        Dx_fwd = suppressWarnings(as.numeric(Dx_fwd))  # <–– force numeric
      ) %>%
      select(-type)
    
    data <- dist_data %>%
      inner_join(dist_fit, by = join_by(tax_name, library_id, x)) %>%
      group_by(tax_name) %>%
      arrange(x, .by_group = TRUE)
    
    if (nrow(dist_data) == 0) {
      return(NULL)
    }
    
    fits <- data %>%
      do({
        ccc <- DescTools::CCC(.$f_fwd, .$Dx_fwd, ci = ci)
        
        if (!is.null(nperm)) {
          ccc_perm <- perk_test(
            .$f_fwd, .$Dx_fwd,
            B = nperm, method = "ccc", alternative = "greater"
          )
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
  return(results)
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
  dat_all <- dmg_fwd_CCC(holi_data_sp_euk, samples, ci = "asymptotic", nperm = 100, nproc = 14)
  
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
summarise_stats_joined_all <- function(
    stats_joined,
    group_cols   = c("library_id", "taxid"),
    weight_col   = "n_reads_from_unicorn",
    exclude_cols = c("accession", "library_id", "taxid", "n_reads_from_unicorn", "n_alns"),
    sum_cols     = c("reference_length")
) {
  stopifnot(all(group_cols %in% names(stats_joined)))
  stopifnot(weight_col %in% names(stats_joined))
  
  # numeric columns
  num_cols <- names(stats_joined)[vapply(stats_joined, is.numeric, logical(1))]
  
  # columns we actually want to summarise
  target_cols <- setdiff(num_cols, c(group_cols, exclude_cols))
  
  # of those, which should be summed (and exist + numeric)
  sum_cols_present <- intersect(sum_cols, target_cols)
  
  # the rest: mean only
  mean_cols <- setdiff(target_cols, sum_cols_present)
  
  dplyr::as_tibble(stats_joined) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n_references = if ("accession" %in% names(stats_joined)) dplyr::n_distinct(accession) else dplyr::n(),
      total_reads_from_unicorn  = sum(.data[[weight_col]], na.rm = TRUE),
      
      # mean for all remaining numeric cols
      dplyr::across(
        dplyr::all_of(mean_cols),
        ~ mean(.x, na.rm = TRUE)
      ),
      
      # sum for reference_length (or any other sum_cols you pass)
      dplyr::across(
        dplyr::all_of(sum_cols_present),
        ~ sum(.x, na.rm = TRUE)
      ),
      
      .groups = "drop"
    )
}

summarise_bf_genus <- function(
    bf_md_data,
    group_cols = c("library_id", "genus"),
    reads_col = "n_reads",
    accession_col = "accession",
    # numeric columns to exclude from mean aggregation
    exclude_from_mean = c("taxid", "n_reads", "reference_length"),
    # numeric columns to sum
    sum_cols = c("reference_length"),
    # how to carry non-numeric columns: "first" or "collapse"
    nonnum_mode = c("first", "collapse")
) {
  nonnum_mode <- match.arg(nonnum_mode)
  
  stopifnot(all(group_cols %in% names(bf_md_data)))
  stopifnot(reads_col %in% names(bf_md_data))
  
  # numeric + non-numeric column sets
  num_cols <- names(bf_md_data)[vapply(bf_md_data, is.numeric, logical(1))]
  nonnum_cols <- setdiff(names(bf_md_data), num_cols)
  
  # don't summarise group cols or accession as carried columns
  nonnum_cols <- setdiff(nonnum_cols, c(group_cols, accession_col))
  
  # mean targets = all numeric except excluded + sums + group cols
  mean_cols <- setdiff(num_cols, c(group_cols, exclude_from_mean, sum_cols))
  
  # sum targets present
  sum_cols_present <- intersect(sum_cols, num_cols)
  
  bf_md_data %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      # requested totals
      total_n_reads = sum(.data[[reads_col]], na.rm = TRUE),
      n_references  = if (accession_col %in% names(bf_md_data)) {
        n_distinct(.data[[accession_col]])
      } else {
        n()
      },
      
      # sum columns (e.g. reference_length)
      across(all_of(sum_cols_present), ~ sum(.x, na.rm = TRUE)),
      
      # carry over ALL other numeric columns via mean
      across(all_of(mean_cols), ~ mean(.x, na.rm = TRUE)),
      
      # carry over ALL non-numeric columns
      across(
        all_of(nonnum_cols),
        ~ if (nonnum_mode == "first") {
          dplyr::first(stats::na.omit(.x))
        } else {
          paste(unique(stats::na.omit(.x)), collapse = "; ")
        }
      ),
      
      .groups = "drop"
    )
}



calculate_plot_grid <- function(n) {
  cols <- ceiling(sqrt(n))
  rows <- ceiling(n / cols)
  
  list(rows = rows, cols = cols)
}


get_dmg_decay_fit <- function(df, orient = "fwd", pos = 25, p_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
  df_dx_fwd <- df %>%
    select(tax_name, library_id, matches("^fwdx\\d+")) %>%
    pivot_longer(names_to = "type", values_to = "Dx_fwd", c(-tax_name, -library_id)) %>%
    mutate(x = gsub("fwdx", "", type)) %>%
    select(-type)
  
  df_dx_rev <- df %>%
    select(tax_name, library_id, matches("^bwdx\\d+")) %>%
    pivot_longer(names_to = "type", values_to = "Dx_rev", c(-tax_name, -library_id)) %>%
    mutate(x = gsub("bwdx", "", type)) %>%
    select(-type)
  
  df_dx_conf_fwd <- df %>%
    select(tax_name, library_id, starts_with("fwdxConf")) %>%
    pivot_longer(names_to = "type", values_to = "Dx_fwd_conf", c(-tax_name, -library_id)) %>%
    mutate(x = gsub("fwdxConf", "", type)) %>%
    select(-type)
  
  df_dx_conf_rev <- df %>%
    select(tax_name, library_id, starts_with("bwdxConf")) %>%
    pivot_longer(names_to = "type", values_to = "Dx_rev_conf", c(-tax_name, -library_id)) %>%
    mutate(x = gsub("bwdxConf", "", type)) %>%
    select(-type)
  
  
  df_fit_fwd <- df %>%
    select(tax_name, library_id, starts_with("fwf")) %>%
    pivot_longer(names_to = "type", values_to = "f_fwd", c(-tax_name, -library_id)) %>%
    mutate(x = gsub("fwf", "", type)) %>%
    select(-type)
  
  df_fit_rev <- df %>%
    select(tax_name, library_id, starts_with("bwf")) %>%
    pivot_longer(names_to = "type", values_to = "f_rev", c(-tax_name, -library_id)) %>%
    mutate(x = gsub("bwf", "", type)) %>%
    select(-type)
  
  dat <- df_dx_fwd %>%
    inner_join(df_dx_rev, by = join_by(tax_name, library_id, x)) %>%
    inner_join(df_dx_conf_fwd, by = join_by(tax_name, library_id, x)) %>%
    inner_join(df_dx_conf_rev, by = join_by(tax_name, library_id, x)) %>%
    inner_join(df_fit_fwd, by = join_by(tax_name, library_id, x)) %>%
    inner_join(df_fit_rev, by = join_by(tax_name, library_id, x)) %>%
    inner_join(df %>% select(library_id) %>% distinct(), by = join_by(library_id)) %>%
    mutate(x = as.numeric(x)) %>%
    filter(x <= pos) %>%
    # inner_join(cdata %>% select(library_id, member_unit)) %>%
    rowwise() %>%
    mutate(
      Dx_fwd_min = f_fwd - Dx_fwd_conf,
      Dx_fwd_max = f_fwd + Dx_fwd_conf,
      Dx_rev_min = f_rev - Dx_rev_conf,
      Dx_rev_max = f_rev + Dx_rev_conf
    )
  
  samples <- dat %>%
    select(library_id, tax_name) %>%
    distinct() %>%
    ungroup()
  
  smooth_it <- function(x, y, n = 1000, method = "natural") {
    t <- seq_along(x)
    new_t <- seq(min(t), max(t), length.out = n)
    new_x <- spline(t, x, xout = new_t, method = method)$y
    new_y <- spline(t, y, xout = new_t, method = method)$y
    data.frame(t = new_t, x = new_x, y = new_y)
  }
  
  
  if (orient == "fwd") {
    dat1 <- pmap_dfr(samples, function(...) {
      current <- tibble(...)
      d <- dat %>%
        filter(library_id == current$library_id, tax_name == current$tax_name)
      d <- smooth_it(d$x, d$Dx_fwd, n = 1000, method = "natural")
      d <- d %>% mutate(
        x = as.numeric(x),
        library_id = current$library_id,
        tax_name = current$tax_name,
      )
    })
    ggplot() +
      geom_ribbon(data = dat, aes(x, ymin = Dx_fwd_min, ymax = Dx_fwd_max, group = interaction(library_id, tax_name)), alpha = 0.1, fill = "#284B63") +
      geom_path(data = dat, aes(x, f_fwd, group = interaction(library_id, tax_name)), color = "#284B63", linewidth = 1.5) +
      geom_point(data = dat, aes(x, Dx_fwd), alpha = .3, size = 2, fill = "black") +
      xlab("Position") +
      ylab("Frequency") +
      theme_bw() +
      # scale_y_continuous(limits = c(y_min, y_max), breaks = p_breaks) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        text = element_text(size = 12),
      )
  } else {
    dat1 <- pmap_dfr(samples, function(...) {
      current <- tibble(...)
      d <- dat %>%
        filter(library_id == current$library_id, tax_name == current$tax_name)
      d <- smooth_it(d$x, d$Dx_rev, n = 1000, method = "natural")
      d <- d %>% mutate(
        x = as.numeric(x),
        library_id = current$library_id,
        tax_name = current$tax_name,
      )
    })
    ggplot() +
      geom_ribbon(data = dat, aes(x, ymin = Dx_rev_min, ymax = Dx_rev_max, group = interaction(library_id, tax_name)), alpha = 0.1, fill = "#B04035") +
      geom_path(data = dat, aes(x, f_rev, group = interaction(get_dmg_decay_fit <- function(df, orient = "fwd", pos = 25, p_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
  df_dx_fwd <- df %>%
    select(tax_name, library_id, matches("^fwdx\\d+")) %>%
    pivot_longer(names_to = "type", values_to = "Dx_fwd", c(-tax_name, -library_id)) %>%
    mutate(x = gsub("fwdx", "", type)) %>%
    select(-type)
  
  df_dx_rev <- df %>%
    select(tax_name, library_id, matches("^bwdx\\d+")) %>%
    pivot_longer(names_to = "type", values_to = "Dx_rev", c(-tax_name, -library_id)) %>%
    mutate(x = gsub("bwdx", "", type)) %>%
    select(-type)
  
  df_dx_conf_fwd <- df %>%
    select(tax_name, library_id, starts_with("fwdxConf")) %>%
    pivot_longer(names_to = "type", values_to = "Dx_fwd_conf", c(-tax_name, -library_id)) %>%
    mutate(x = gsub("fwdxConf", "", type)) %>%
    select(-type)
  
  df_dx_conf_rev <- df %>%
    select(tax_name, library_id, starts_with("bwdxConf")) %>%
    pivot_longer(names_to = "type", values_to = "Dx_rev_conf", c(-tax_name, -library_id)) %>%
    mutate(x = gsub("bwdxConf", "", type)) %>%
    select(-type)
  
  
  df_fit_fwd <- df %>%
    select(tax_name, library_id, starts_with("fwf")) %>%
    pivot_longer(names_to = "type", values_to = "f_fwd", c(-tax_name, -library_id)) %>%
    mutate(x = gsub("fwf", "", type)) %>%
    select(-type)
  
  df_fit_rev <- df %>%
    select(tax_name, library_id, starts_with("bwf")) %>%
    pivot_longer(names_to = "type", values_to = "f_rev", c(-tax_name, -library_id)) %>%
    mutate(x = gsub("bwf", "", type)) %>%
    select(-type)
  
  dat <- df_dx_fwd %>%
    inner_join(df_dx_rev, by = join_by(tax_name, library_id, x)) %>%
    inner_join(df_dx_conf_fwd, by = join_by(tax_name, library_id, x)) %>%
    inner_join(df_dx_conf_rev, by = join_by(tax_name, library_id, x)) %>%
    inner_join(df_fit_fwd, by = join_by(tax_name, library_id, x)) %>%
    inner_join(df_fit_rev, by = join_by(tax_name, library_id, x)) %>%
    inner_join(df %>% select(library_id) %>% distinct(), by = join_by(library_id)) %>%
    mutate(x = as.numeric(x)) %>%
    filter(x <= pos) %>%
    # inner_join(cdata %>% select(library_id, member_unit)) %>%
    rowwise() %>%
    mutate(
      Dx_fwd_min = f_fwd - Dx_fwd_conf,
      Dx_fwd_max = f_fwd + Dx_fwd_conf,
      Dx_rev_min = f_rev - Dx_rev_conf,
      Dx_rev_max = f_rev + Dx_rev_conf
    )
  
  samples <- dat %>%
    select(library_id, tax_name) %>%
    distinct() %>%
    ungroup()
  
  smooth_it <- function(x, y, n = 1000, method = "natural") {
    t <- seq_along(x)
    new_t <- seq(min(t), max(t), length.out = n)
    new_x <- spline(t, x, xout = new_t, method = method)$y
    new_y <- spline(t, y, xout = new_t, method = method)$y
    data.frame(t = new_t, x = new_x, y = new_y)
  }
  
  
  if (orient == "fwd") {
    dat1 <- pmap_dfr(samples, function(...) {
      current <- tibble(...)
      d <- dat %>%
        filter(library_id == current$library_id, tax_name == current$tax_name)
      d <- smooth_it(d$x, d$Dx_fwd, n = 1000, method = "natural")
      d <- d %>% mutate(
        x = as.numeric(x),
        library_id = current$library_id,
        tax_name = current$tax_name,
      )
    })
    ggplot() +
      geom_ribbon(data = dat, aes(x, ymin = Dx_fwd_min, ymax = Dx_fwd_max, group = interaction(library_id, tax_name)), alpha = 0.1, fill = "#284B63") +
      geom_path(data = dat, aes(x, f_fwd, group = interaction(library_id, tax_name)), color = "#284B63", linewidth = 1.5) +
      geom_point(data = dat, aes(x, Dx_fwd), alpha = .3, size = 2, fill = "black") +
      xlab("Position") +
      ylab("Frequency") +
      theme_bw() +
      # scale_y_continuous(limits = c(y_min, y_max), breaks = p_breaks) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        text = element_text(size = 12),
      )
  } else {
    dat1 <- pmap_dfr(samples, function(...) {
      current <- tibble(...)
      d <- dat %>%
        filter(library_id == current$library_id, tax_name == current$tax_name)
      d <- smooth_it(d$x, d$Dx_rev, n = 1000, method = "natural")
      d <- d %>% mutate(
        x = as.numeric(x),
        library_id = current$library_id,
        tax_name = current$tax_name,
      )
    })
    ggplot() +
      geom_ribbon(data = dat, aes(x, ymin = Dx_rev_min, ymax = Dx_rev_max, group = interaction(library_id, tax_name)), alpha = 0.1, fill = "#B04035") +
      geom_path(data = dat, aes(x, f_rev, group = interaction(librart_id, tax_name)), color = "#B04035", linewidth = 1.5, linejoin = "round") + # geom_point(data = dat, aes(x, f_rev), alpha = .3, size = 2, fill = "black") +
      # stat_summary(data = dat, aes(x, f_rev), fun.data = mean_sd, geom = "ribbon", alpha = .3, size = 1, fill = "black") +
      # stat_summary(data = dat, aes(x, f_rev), fun.data = mean_sd, geom = "errorbar", alpha = .3, size = 1, fill = "black") +
      # stat_summary(data = dat, aes(x, f_rev), fun = mean, geom = "line", size = 0.8, color = "black") +
      xlab("Position") +
      ylab("Frequency") +
      theme_bw() +
      scale_x_continuous(trans = "reverse") +
      scale_y_continuous(position = "right") +
      # scale_y_continuous(limits = c(y_min, y_max), position = "right", breaks = p_breaks) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        text = element_text(size = 12),
      )
  }
}
, tax_name)), color = "#B04035", linewidth = 1.5, linejoin = "round") + # geom_point(data = dat, aes(x, f_rev), alpha = .3, size = 2, fill = "black") +
      # stat_summary(data = dat, aes(x, f_rev), fun.data = mean_sd, geom = "ribbon", alpha = .3, size = 1, fill = "black") +
      # stat_summary(data = dat, aes(x, f_rev), fun.data = mean_sd, geom = "errorbar", alpha = .3, size = 1, fill = "black") +
      # stat_summary(data = dat, aes(x, f_rev), fun = mean, geom = "line", size = 0.8, color = "black") +
      xlab("Position") +
      ylab("Frequency") +
      theme_bw() +
      scale_x_continuous(trans = "reverse") +
      scale_y_continuous(position = "right") +
      # scale_y_continuous(limits = c(y_min, y_max), position = "right", breaks = p_breaks) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        text = element_text(size = 12),
      )
  }
}



