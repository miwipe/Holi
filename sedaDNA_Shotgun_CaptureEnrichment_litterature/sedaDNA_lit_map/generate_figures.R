# ===============================================================
# Ancient Metagenomic Publication Data Visualization Pipeline
# ===============================================================
# Description:
# This script processes and visualizes geospatial and temporal
# trends in ancient metagenomic studies using data from a Google Sheet.
#
# Main outputs:
#   - Map by molecular method
#   - Publication counts by year
#   - Cumulative publication counts by molecular method
#   - Yearly publication counts by molecular method
#   - Publications by reference database
#   - Publications by mapper
#   - Map colored by SampleType and shaped by TargetGroup
#
# Output figures are written to ../../figures/
# ===============================================================


# ---------------------------------------------------------------
# 1. Load required libraries
# ---------------------------------------------------------------

required_packages <- c(
  "sf",
  "ggplot2",
  "ggrepel",
  "rnaturalearth",
  "rnaturalearthdata",
  "readr",
  "readxl",
  "tidyverse",
  "googlesheets4",
  "scales",
  "viridis",
  "RColorBrewer"
)

installed_packages <- rownames(installed.packages())

for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

invisible(
  lapply(required_packages, library, character.only = TRUE)
)


# ---------------------------------------------------------------
# 2. Output directory
# ---------------------------------------------------------------

dir.create(
  "../../figures",
  recursive = TRUE,
  showWarnings = FALSE
)


# ---------------------------------------------------------------
# 3. Authenticate and read Google Sheet
# ---------------------------------------------------------------

gs4_auth()

coordinates <- read_sheet(
  "https://docs.google.com/spreadsheets/d/13cmBUi4cigUaTKtQeFLFvS0gXT8AeWxWKzHv2UcOBCI/edit?gid=0#gid=0",
  sheet = "Eukaryotes"
)


# ---------------------------------------------------------------
# 4. Basic data cleaning
# ---------------------------------------------------------------

coordinates <- coordinates %>%
  mutate(
    across(
      where(is.character),
      ~ trimws(.x)
    )
  ) %>%
  mutate(
    Latitude = na_if(as.character(Latitude), ""),
    Longitude = na_if(as.character(Longitude), ""),
    TargetGroup = na_if(as.character(TargetGroup), ""),
    TargetTaxa = na_if(as.character(TargetTaxa), ""),
    MolecularMethod = na_if(as.character(MolecularMethod), ""),
    SampleType = na_if(as.character(SampleType), ""),
    SiteName = na_if(as.character(SiteName), ""),
    DOI = na_if(as.character(DOI), ""),
    mapper = na_if(as.character(mapper), ""),
    reference_target = na_if(as.character(reference_target), "")
  )


# ---------------------------------------------------------------
# 5. Coordinate cleaning and automatic lat/lon swap detection
# ---------------------------------------------------------------

coordinates <- coordinates %>%
  mutate(
    Latitude_numeric = suppressWarnings(
      as.numeric(as.character(Latitude))
    ),
    Longitude_numeric = suppressWarnings(
      as.numeric(as.character(Longitude))
    ),
    
    # Obvious reversal: latitude impossible (>90) while longitude fits latitude range
    coords_swapped =
      !is.na(Latitude_numeric) &
      !is.na(Longitude_numeric) &
      abs(Latitude_numeric) > 90 &
      abs(Longitude_numeric) <= 90,
    
    Latitude_fixed = if_else(
      coords_swapped,
      Longitude_numeric,
      Latitude_numeric
    ),
    
    Longitude_fixed = if_else(
      coords_swapped,
      Latitude_numeric,
      Longitude_numeric
    )
  )

# Report automatically swapped coordinates
swapped_coordinates <- coordinates %>%
  filter(coords_swapped) %>%
  select(
    any_of(c(
      "Publication",
      "SiteName",
      "Latitude",
      "Longitude",
      "Latitude_fixed",
      "Longitude_fixed"
    ))
  )

if (nrow(swapped_coordinates) > 0) {
  message("\nCoordinates automatically swapped because latitude/longitude were reversed:")
  print(swapped_coordinates, n = Inf)
}

# Report rows that still cannot be mapped
bad_coordinates <- coordinates %>%
  filter(
    is.na(Latitude_fixed) |
      is.na(Longitude_fixed) |
      Latitude_fixed < -90 |
      Latitude_fixed > 90 |
      Longitude_fixed < -180 |
      Longitude_fixed > 180
  ) %>%
  select(
    any_of(c(
      "Publication",
      "SiteName",
      "Latitude",
      "Longitude"
    ))
  )

if (nrow(bad_coordinates) > 0) {
  message("\nRows excluded from maps because coordinates are missing or invalid:")
  print(bad_coordinates, n = Inf)
}


# ---------------------------------------------------------------
# 6. Keep rows with usable geographic coordinates
# ---------------------------------------------------------------

coordinates_noNA <- coordinates %>%
  filter(
    !is.na(Latitude_fixed),
    !is.na(Longitude_fixed),
    Latitude_fixed >= -90,
    Latitude_fixed <= 90,
    Longitude_fixed >= -180,
    Longitude_fixed <= 180,
    is.na(TargetGroup) | TargetGroup != "Microorganisms | NA"
  ) %>%
  mutate(
    Latitude = Latitude_fixed,
    Longitude = Longitude_fixed
  )


# ---------------------------------------------------------------
# 7. World map and projection
# ---------------------------------------------------------------

world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

robin_crs <- st_crs("+proj=robin")

world_proj <- st_transform(
  world,
  crs = robin_crs
)

coordinates_sf <- st_as_sf(
  coordinates_noNA,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)

coordinates_proj_sf <- st_transform(
  coordinates_sf,
  crs = robin_crs
)

xy <- st_coordinates(coordinates_proj_sf)

coordinates_proj <- coordinates_noNA %>%
  mutate(
    X = xy[, 1],
    Y = xy[, 2]
  )


# ---------------------------------------------------------------
# 8. Common filtered eukaryotic plotting dataset
# ---------------------------------------------------------------

coordinates_proj_SG_TE <- coordinates_proj %>%
  mutate(
    Lab = as.factor(Lab),
    year_published = suppressWarnings(
      as.integer(as.character(year_published))
    )
  ) %>%
  filter(
    !is.na(year_published),
    is.na(TargetGroup) | TargetGroup != "Microorganisms",
    is.na(TargetTaxa) |
      !TargetTaxa %in% c(
        "Microorganisms",
        "Prokaryotes"
      )
  )


# ---------------------------------------------------------------
# 9. Helper function for discrete shapes
# ---------------------------------------------------------------

make_shape_scale <- function(values) {
  
  levels_present <- sort(
    unique(
      values[
        !is.na(values) &
          values != ""
      ]
    )
  )
  
  available_shapes <- c(
    16, 17, 15, 18,
    3, 4, 8, 1, 2,
    0, 5, 6, 7, 9,
    10, 11, 12, 13, 14,
    19, 20, 21, 22, 23, 24, 25
  )
  
  shapes <- rep(
    available_shapes,
    length.out = length(levels_present)
  )
  
  names(shapes) <- levels_present
  
  shapes
}


# ===============================================================
# FIGURE 1
# Sites by Molecular Method
# ===============================================================

map_method_data <- coordinates_proj_SG_TE %>%
  filter(
    !is.na(MolecularMethod),
    MolecularMethod != ""
  )

method_shapes <- make_shape_scale(
  map_method_data$MolecularMethod
)

p_map_method <- ggplot() +
  
  geom_sf(
    data = world_proj,
    fill = "lightgrey",
    color = "black",
    linewidth = 0.25
  ) +
  
  geom_point(
    data = map_method_data,
    aes(
      x = X,
      y = Y,
      shape = MolecularMethod
    ),
    color = "black",
    size = 3,
    na.rm = TRUE
  ) +
  
  scale_shape_manual(
    values = method_shapes,
    drop = TRUE
  ) +
  
  
  coord_sf(
    crs = robin_crs
  ) +
  
  theme_minimal() +
  
  theme(
    panel.grid.major = element_line(
      color = "gray",
      linetype = "dashed"
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    legend.position = "right"
  ) +
  
  labs(
    title = "Ancient Metagenomic Study Sites",
    x = "Longitude",
    y = "Latitude",
    shape = "Molecular Method"
  )


ggsave(
  "../../figures/SG_TE_map_method.png",
  plot = p_map_method,
  width = 10,
  height = 7,
  dpi = 300
)


# ---------------------------------------------------------------
# 10. Publication-level plotting dataset
# ---------------------------------------------------------------

publication_data <- coordinates %>%
  mutate(
    year_published = suppressWarnings(
      as.integer(as.character(year_published))
    )
  ) %>%
  filter(
    !is.na(year_published),
    is.na(TargetGroup) | TargetGroup != "Microorganisms",
    is.na(TargetTaxa) |
      !TargetTaxa %in% c(
        "Microorganisms",
        "Prokaryotes"
      )
  )


# ===============================================================
# FIGURE 2
# Number of publications per year
# ===============================================================

publication_year_counts <- publication_data %>%
  filter(
    !is.na(DOI),
    DOI != "",
    !is.na(year_published)
  ) %>%
  distinct(
    year_published,
    DOI
  ) %>%
  count(
    year_published,
    name = "n"
  ) %>%
  arrange(year_published) %>%
  mutate(
    year_published = factor(
      year_published,
      levels = sort(unique(year_published))
    )
  )

p_publications_year <- ggplot(
  publication_year_counts,
  aes(
    x = year_published,
    y = n
  )
) +
  geom_col(
    width = 0.8,
    na.rm = TRUE
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_test() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  labs(
    title = "Ancient Metagenome Publication Counts by Year",
    x = "Year Published",
    y = "Number of Publications"
  )


ggsave(
  "../../figures/barplot_no_publications.png",
  plot = p_publications_year,
  width = 6,
  height = 4,
  dpi = 300
)


# ===============================================================
# FIGURE 3
# Cumulative publication trend by Molecular Method
# ===============================================================

cumulative_data <- publication_data %>%
  filter(
    !is.na(DOI),
    DOI != "",
    !is.na(MolecularMethod),
    MolecularMethod != "",
    !is.na(year_published)
  ) %>%
  distinct(
    year_published,
    MolecularMethod,
    DOI
  ) %>%
  count(
    year_published,
    MolecularMethod,
    name = "n"
  ) %>%
  tidyr::complete(
    year_published = seq(
      min(year_published, na.rm = TRUE),
      max(year_published, na.rm = TRUE),
      by = 1
    ),
    MolecularMethod,
    fill = list(n = 0)
  ) %>%
  arrange(
    MolecularMethod,
    year_published
  ) %>%
  group_by(MolecularMethod) %>%
  mutate(
    cumulative_n = cumsum(n)
  ) %>%
  ungroup() %>%
  mutate(
    year_published = factor(
      year_published,
      levels = sort(unique(year_published))
    )
  )

p_cumulative <- ggplot(
  cumulative_data,
  aes(
    x = year_published,
    y = cumulative_n,
    fill = MolecularMethod
  )
) +
  geom_col(
    width = 0.8,
    na.rm = TRUE
  ) +
  scale_fill_viridis_d(
    option = "D"
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_test() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  labs(
    title = "Cumulative Publications Over Time by Molecular Method",
    x = "Year Published",
    y = "Cumulative Number of Publications",
    fill = "Molecular Method"
  )


ggsave(
  "../../figures/barplot_cumsum_no_publications_methods.png",
  plot = p_cumulative,
  width = 7,
  height = 5,
  dpi = 300
)


# ===============================================================
# FIGURE 4
# Yearly publication counts by Molecular Method
# ===============================================================

year_method_counts <- publication_data %>%
  filter(
    !is.na(DOI),
    DOI != "",
    !is.na(MolecularMethod),
    MolecularMethod != "",
    !is.na(year_published)
  ) %>%
  distinct(
    year_published,
    MolecularMethod,
    DOI
  ) %>%
  count(
    year_published,
    MolecularMethod,
    name = "n"
  ) %>%
  arrange(year_published) %>%
  mutate(
    year_published = factor(
      year_published,
      levels = sort(unique(year_published))
    )
  )

p_year_method <- ggplot(
  year_method_counts,
  aes(
    x = year_published,
    y = n,
    fill = MolecularMethod
  )
) +
  geom_col(
    width = 0.8,
    na.rm = TRUE
  ) +
  scale_fill_viridis_d(
    option = "D"
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_test() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  labs(
    title = "Ancient Metagenome Publication Counts by Year and Molecular Method",
    x = "Year Published",
    y = "Number of Publications",
    fill = "Molecular Method"
  )


ggsave(
  "../../figures/barplot_no_publications_methods.png",
  plot = p_year_method,
  width = 8,
  height = 5,
  dpi = 300
)


# ===============================================================
# FIGURE 5
# Publications by Reference Database
# ===============================================================

reference_data <- publication_data %>%
  filter(
    !is.na(reference_target),
    reference_target != "",
    !is.na(DOI),
    DOI != ""
  ) %>%
  distinct(
    reference_target,
    DOI
  ) %>%
  count(
    reference_target,
    name = "n"
  ) %>%
  arrange(desc(n))

p_reference <- ggplot(
  reference_data,
  aes(
    x = reorder(reference_target, n),
    y = n
  )
) +
  geom_col(
    na.rm = TRUE
  ) +
  coord_flip() +
  scale_y_continuous(
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_test() +
  labs(
    title = "Publications by Reference Databases",
    x = "Reference Target",
    y = "Number of Publications"
  )


ggsave(
  "../../figures/barplot_N_databases.png",
  plot = p_reference,
  width = 8,
  height = 6,
  dpi = 300
)


# ===============================================================
# FIGURE 6
# Publications by Mapper
# ===============================================================

mapper_data <- publication_data %>%
  filter(
    !is.na(mapper),
    mapper != "",
    !is.na(DOI),
    DOI != ""
  ) %>%
  distinct(
    mapper,
    DOI
  ) %>%
  count(
    mapper,
    name = "n"
  ) %>%
  arrange(desc(n))

p_mapper <- ggplot(
  mapper_data,
  aes(
    x = reorder(mapper, n),
    y = n
  )
) +
  geom_col(
    na.rm = TRUE
  ) +
  coord_flip() +
  scale_y_continuous(
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_test() +
  labs(
    title = "Publications by Mapper",
    x = "Mapper",
    y = "Number of Publications"
  )


ggsave(
  "../../figures/barplot_N_mappers.png",
  plot = p_mapper,
  width = 8,
  height = 6,
  dpi = 300
)


# ===============================================================
# FIGURE 7
# Map colored by Sample Type and shaped by Target Group
# ===============================================================

map_target_data <- coordinates_proj_SG_TE %>%
  filter(
    !is.na(SampleType),
    SampleType != "",
    !is.na(TargetGroup),
    TargetGroup != ""
  )

message("\nSampleType categories used in final map:")
print(
  map_target_data %>%
    count(
      SampleType,
      sort = TRUE
    ),
  n = Inf
)

message("\nTargetGroup categories used in final map:")
print(
  map_target_data %>%
    count(
      TargetGroup,
      sort = TRUE
    ),
  n = Inf
)

target_shapes <- make_shape_scale(
  map_target_data$TargetGroup
)

p_map_target <- ggplot() +
  
  geom_sf(
    data = world_proj,
    fill = "lightgrey",
    color = "black",
    linewidth = 0.25
  ) +
  
  geom_point(
    data = map_target_data,
    aes(
      x = X,
      y = Y,
      shape = TargetGroup,
      color = SampleType
    ),
    size = 3,
    na.rm = TRUE
  ) +
  
  scale_shape_manual(
    values = target_shapes,
    drop = TRUE
  ) +
  
  # Dynamic discrete palette: works with any number of SampleType categories
  scale_color_viridis_d(
    option = "D"
  ) +
  
  coord_sf(
    crs = robin_crs
  ) +
  
  theme_minimal() +
  
  theme(
    panel.grid.major = element_line(
      color = "gray",
      linetype = "dashed"
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    legend.text = element_text(
      size = 9
    )
  ) +
  
  labs(
    title = "Ancient Metagenomic Study Sites",
    subtitle = "Global distribution of sampling locations by Sample Type and Studied Organism",
    x = "Longitude",
    y = "Latitude",
    shape = "Target Organism",
    color = "Sample Type"
  )


ggsave(
  "../../figures/SG_TE_map_SampleType_Target.png",
  plot = p_map_target,
  width = 10,
  height = 7,
  dpi = 300
)


# ===============================================================
# FINAL SUMMARY
# ===============================================================

message("\n---------------------------------------")
message("Finished generating figures.")
message("Rows read from Google Sheet: ", nrow(coordinates))
message("Rows with valid map coordinates: ", nrow(coordinates_noNA))
message(
  "Unique publications: ",
  n_distinct(
    coordinates$DOI[
      !is.na(coordinates$DOI) &
        coordinates$DOI != ""
    ]
  )
)
message(
  "Unique mapped sites: ",
  n_distinct(
    coordinates_noNA$SiteName[
      !is.na(coordinates_noNA$SiteName) &
        coordinates_noNA$SiteName != ""
    ]
  )
)
message("---------------------------------------\n")


#### END OF SCRIPT
