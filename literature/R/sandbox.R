
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
