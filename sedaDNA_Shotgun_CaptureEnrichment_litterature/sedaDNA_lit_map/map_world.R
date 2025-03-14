
# Load required libraries
library(sf)
library(ggplot2)
library(ggrepel)
library(rnaturalearth)
library(rnaturalearthdata)
library(readr)
library(readxl)
library(tidyverse)
library(googlesheets4)
library(dplyr)
library(scales)
library(viridis) # For a better color palette



# Authenticate with Google (only required once)
gs4_auth()

# Read the Google Sheet
coordinates <- read_sheet("https://docs.google.com/spreadsheets/d/13cmBUi4cigUaTKtQeFLFvS0gXT8AeWxWKzHv2UcOBCI/edit?gid=0#gid=0", sheet = "Sheet1")


colnames(coordinates)
unique(coordinates$TargetGroup)
coordinates_noNA <- coordinates %>% filter(Latitude != "NA", TargetGroup != "Microorganisms | NA") 

# Load world dataset
world <- ne_countries(scale = "medium", returnclass = "sf")

# Reproject World to Robinson projection using +proj=robin
world_proj <- st_transform(world, crs = st_crs("+proj=robin"))


# Reproject coordinates to Robinson projection
coordinates_sf <- st_as_sf(coordinates_noNA, coords = c("Longitude", "Latitude"), crs = 4326)
coordinates_proj <- st_transform(coordinates_sf, crs = st_crs("+proj=robin"))
coordinates_proj <- cbind(coordinates_noNA, st_coordinates(coordinates_proj))


coordinates_proj_SG_TE <- coordinates_proj %>%
  filter(year_published != "NA", TargetGroup != "Microorganisms", TargetTaxa != "Microorganisms",  TargetTaxa != "Prokaryotes")

ggplot() +
  # Plot World (reprojected)
  geom_sf(data = world_proj, fill = "lightgrey", color = "black") +
  # Add points for the coordinates
  geom_point(data = coordinates_proj_SG_TE, 
             aes(x = X, y = Y, shape = MolecularMethod), # No color mapping here
             color = "black", size = 3) + # Set fixed color to black
  # Add text labels for sites
  geom_text_repel(data = coordinates_proj_SG_TE, 
                  aes(x = X, y = Y, label = SiteName),
                  color = "black", size = 3, fontface = "bold", box.padding = 0.3) +
  # Apply Robinson projection
  coord_sf(crs = st_crs("+proj=robin")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray", linetype = "dashed"),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Ancient Metagenomic Study Sites",
    x = "Longitude",
    y = "Latitude",
    shape = "Molecular Method" # Legend title for shape
  )

ggsave("../../figures/SG_TE_map_method.png", width = 10, height = 7, dpi = 300)


coordinates_proj_SG_TE <- coordinates_proj %>%
  mutate(Lab = as.factor(Lab)) +  # Ensure Lab is a factor
  filter(year_published != "NA", TargetGroup != "Microorganisms", TargetTaxa != "Microorganisms",  TargetTaxa != "Prokaryotes")


# count publications per year and make barplot
coordinates %>%
  filter(
    year_published != "NA", 
    TargetGroup != "Microorganisms", 
    TargetTaxa != "Microorganisms"
  ) %>%
  select(year_published, MolecularMethod, DOI) %>% # Include MolecularMethod
  unique() %>%
  group_by(year_published, MolecularMethod) %>% # Group by year and MolecularMethod
  count() %>%
  ungroup() %>%
  ggplot(aes(x = year_published, y = n)) + # Use fill for stacking
  geom_bar(stat = "identity") + # Default is position = "stack"
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5) # Ensure round number breaks
  ) +
  theme_test() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + # Rotate x-axis labels
  labs(
    title = "Ancient Metagenome Publication Counts by Year",
    x = "Year Published",
    y = "Number of Publications",
    fill = "Molecular Method" # Legend title
  )

ggsave("../../figures/barplot_no_publications.png", width = 6, height = 4, dpi = 300)



coordinates %>%
  filter(
    year_published != "NA", 
    TargetGroup != "Microorganisms", 
    TargetTaxa != "Microorganisms"
  ) %>%
  select(year_published, MolecularMethod, DOI) %>%
  unique() %>%
  group_by(year_published, MolecularMethod) %>%
  count() %>%
  ungroup() %>%
  arrange(year_published) %>% # Ensure chronological order
  group_by(MolecularMethod) %>%
  mutate(cumulative_n = cumsum(n)) %>% # Compute cumulative sum
  ggplot(aes(x = year_published, y = cumulative_n, fill = MolecularMethod)) + 
  geom_bar(stat = "identity") + 
  scale_fill_brewer(palette = "Set2") + # Use a vibrant color scheme
  scale_y_continuous(
    breaks = pretty_breaks(n = 5) 
  ) + 
  theme_minimal(base_size = 14) + # Use a cleaner theme with larger text
  theme_test() + 
  labs(
    title = "",
    x = "Year Published",
    y = "Cumulative Number of Publications",
    fill = "Molecular Method"
  )

ggsave("../../figures/barplot_cumsum_no_publications_methods.png", width = 6, height = 4, , dpi = 300)


coordinates %>%
  filter(
    year_published != "NA", 
    TargetGroup != "Microorganisms", 
    TargetTaxa != "Microorganisms"
  ) %>%
  select(year_published, MolecularMethod, DOI) %>%
  unique() %>%
  group_by(year_published, MolecularMethod) %>%
  count() %>%
  ungroup() %>%
  arrange(year_published) %>% # Ensure chronological order
  group_by(MolecularMethod) %>%
  mutate(cumulative_n = cumsum(n)) %>% # Compute cumulative sum
  ggplot(aes(x = year_published, y = cumulative_n)) + 
  geom_bar(stat = "identity") + 
  scale_fill_brewer(palette = "Set2") + # Use a vibrant color scheme
  scale_y_continuous(
    breaks = pretty_breaks(n = 5) 
  ) + 
  theme_minimal(base_size = 14) + # Use a cleaner theme with larger text
  theme_test() + 
  labs(
    title = "",
    x = "Year Published",
    y = "Cumulative Number of Publications",
    fill = "Molecular Method"
  )

ggsave("barplot_cumsum_no_publications.pdf", width = 6, height = 4)



coordinates %>%
  filter(
    year_published != "NA", 
    TargetGroup != "Microorganisms", 
    TargetTaxa != "Microorganisms"
  ) %>%
  select(year_published, MolecularMethod, DOI) %>% # Include MolecularMethod
  unique() %>%
  group_by(year_published, MolecularMethod) %>% # Group by year and MolecularMethod
  count() %>%
  ungroup() %>%
  ggplot(aes(x = year_published, y = n, fill = MolecularMethod)) + # Use fill for stacking
  geom_bar(stat = "identity") + # Default is position = "stack"
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5) # Ensure round number breaks
  ) +
  theme_test() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) + # Rotate x-axis labels
  labs(
    title = "Ancient Metagenome Publication Counts by Year and Molecular Method",
    x = "Year Published",
    y = "Number of Publications",
    fill = "Molecular Method" # Legend title
  )

ggsave("barplot_no_publications_methods.pdf", width = 8, height = 5)




# count publications per year with unique reference databases and make barplot
coordinates %>%
  filter(
    year_published != "NA", 
    TargetGroup != "Microorganisms", 
    TargetTaxa != "Microorganisms",  
    TargetTaxa != "Prokaryotes"
  ) %>%
  select(reference_target, DOI) %>%
  unique() %>%
  group_by(reference_target) %>%
  count() %>%
  ungroup() %>%
  ggplot(aes(reference_target, n)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + # Rotate text
  labs(
    title = "Publications by Reference Databases",
    x = "Reference Target",
    y = "Number of Publications"
  )
ggsave("barplot_N_databases.pdf", width = 8, height = 5)


# count publications per year with unique reference databases and make barplot
coordinates %>%
  filter(
    year_published != "NA", 
    TargetGroup != "Microorganisms", 
    TargetTaxa != "Microorganisms",  
    TargetTaxa != "Prokaryotes"
  ) %>%
  select(mapper, DOI) %>%
  unique() %>%
  group_by(mapper) %>%
  count() %>%
  ungroup() %>%
  ggplot(aes(mapper, n)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + # Rotate text
  labs(
    title = "Publications by mapper",
    x = "mapping strategy",
    y = "Number of Publications"
  )
ggsave("barplot_N_mappers.pdf", width = 8, height = 5)


coordinates %>%
  filter(
    year_published != "NA", 
    TargetGroup != "Microorganisms", 
    TargetTaxa != "Microorganisms",  
    TargetTaxa != "Prokaryotes"
  ) |> select(SiteName) |>
  unique()

# count publications per year with unique classifiers and make barplot
coordinates %>%
  filter(
    year_published != "NA", 
    TargetGroup != "Microorganisms", 
    TargetTaxa != "Microorganisms",  
    TargetTaxa != "Prokaryotes"
  ) %>%
  select(classifier, DOI) %>%
  unique() %>%
  group_by(classifier) %>%
  count() %>%
  ungroup() %>%
  ggplot(aes(classifier, n)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + # Rotate text
  labs(
    title = "Publications by Classifiers",
    x = "Classifier",
    y = "Number of Publications"
  )

ggsave("barplot_N_classifiers.pdf", width = 8, height = 5)

colnames(coordinates)
coordinates %>%
  filter(
    year_published != "NA", 
    TargetGroup != "Microorganisms", 
    TargetTaxa != "Microorganisms",  
    TargetTaxa != "Prokaryotes"
  ) %>%
  select(authenticator, authenticated_taxa, DOI) %>% # Include MolecularMethod
  unique() %>%
  group_by(authenticator, authenticated_taxa) %>% # Group by year and MolecularMethod
  count() %>%
  ungroup() %>%
  ggplot(aes(x = authenticator, y = n, fill = authenticated_taxa)) + # Use fill for stacking
  geom_bar(stat = "identity") + # Default is position = "stack"
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5) # Ensure round number breaks
  ) +
  theme_test() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + # Rotate x-axis labels
  labs(
    title = "Ancient Metagenome by Authentication Tool",
    x = "Authentication tool",
    y = "Number of Publications",
    fill = "Fraction of taxa" # Legend title
  )

ggsave("barplot_authenticaiton_methods.pdf", width = 8, height = 5)



coordinates %>%
  filter(
    year_published != "NA", 
    TargetGroup != "Microorganisms", 
    TargetTaxa != "Microorganisms",  
    TargetTaxa != "Prokaryotes"
  ) %>%
  select(year_published, reference_target, DOI) %>%
  unique() %>%
  group_by(year_published, reference_target) %>%
  count() %>%
  ungroup() %>%
  ggplot(aes(x = year_published, y = n, color = reference_target, group = reference_target)) +
  geom_line(size = 1) + # Add lines connecting points
  geom_point(size = 3) + # Add points for each observation
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1) # Rotate x-axis labels
  ) +
  labs(
    title = "Number of Publications by Year and Reference Target",
    x = "Year Published",
    y = "Number of Publications",
    color = "Reference Target" # Legend title
  )

ggsave("line_N_databases.pdf", width = 8, height = 5)


coordinates %>%
  filter(
    year_published != "NA", 
    TargetGroup != "Microorganisms", 
    TargetTaxa != "Microorganisms",  
    TargetTaxa != "Prokaryotes"
  ) %>%
  select(year_published, classifier, DOI) %>%
  unique() %>%
  group_by(year_published, classifier) %>%
  count() %>%
  ungroup() %>%
  ggplot(aes(x = year_published, y = n, color = classifier, group = classifier)) +
  geom_line(size = 1) + # Add lines connecting points
  geom_point(size = 3) + # Add points for each observation
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1) # Rotate x-axis labels
  ) +
  labs(
    title = "Number of Publications by Year by Classifier",
    x = "Year Published",
    y = "Number of Publications",
    color = "Classifier" # Legend title
  )


ggsave("lineplot_N_classifiers.pdf", width = 8, height = 5)


library(viridis)

coordinates %>%
  filter(
    year_published != "NA", 
    TargetGroup != "Microorganisms", 
    TargetTaxa != "Microorganisms",  
    TargetTaxa != "Prokaryotes"
  ) %>%
  select(authenticator, authenticated_taxa, DOI) %>% # Include MolecularMethod
  unique() %>%
  group_by(authenticator, authenticated_taxa) %>% # Group by year and MolecularMethod
  count() %>%
  ungroup() %>%
  ggplot(aes(x = authenticator, y = n, fill = authenticated_taxa)) + # Use fill for stacking
  geom_bar(stat = "identity") + # Default is position = "stack"
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5) # Ensure round number breaks
  ) +
  scale_fill_viridis_d(option = "D") + # Apply color-blind friendly palette
  theme_test() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + # Rotate x-axis labels
  labs(
    title = "Ancient Metagenome by Authentication Tool",
    x = "Authentication tool",
    y = "Number of Publications",
    fill = "Fraction of taxa" # Legend title
  )

ggsave("barplot_authenticaiton_methods.pdf", width = 8, height = 5)


