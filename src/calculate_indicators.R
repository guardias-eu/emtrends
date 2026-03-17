# Load libraries
library(dplyr)
library(ggplot2)
library(here)
library(purrr)
library(readr)
library(rgbif)
library(stringr)
library(tidyr)

# Load the most recent species occurrence cube ####

## Load the list of cubes and filter for the most recent cube ####
list_cubes_path <- "https://raw.githubusercontent.com/guardias-eu/build-eu-cube/refs/heads/main/data/output/list_downloads.tsv"
list_cubes <- readr::read_tsv(list_cubes_path, na = "")
last_cube <- list_cubes %>%
  dplyr::filter(gbif_download_created == max(gbif_download_created, na.rm = TRUE))
last_cube_key <- last_cube$gbif_download_key

## Download the cube from GBIF ####
cube_path <- here::here("data", "input")
cube_zip <- rgbif::occ_download_get(key = last_cube_key, path = cube_path, overwrite = TRUE)

## Unzip and load the cube ####
unzip(zipfile = cube_zip, exdir = cube_path)
cube <- readr::read_tsv(
  here::here(cube_path, paste(last_cube_key, "csv", sep = ".")),
  col_types = readr::cols(
    specieskey = readr::col_integer(),
    year = readr::col_integer(),
    occurrences = readr::col_integer()
  ),
  na = ""
)

# Get the mapping of the grid cells to LMEs ####
grid_cells_with_lme_info_file <- "data/output/grid_cells_with_lme_info.csv"
grid_cells_with_lme_info <- readr::read_csv(grid_cells_with_lme_info_file, na = "")
lme_ids <- unique(grid_cells_with_lme_info$lme_id)
names(lme_ids) <- unique(grid_cells_with_lme_info$lme_name)

# Get the species list used to build the last cube ####

## Take the part of the string starting with `data/output/` ####
last_cube_species_list <- stringr::str_extract(
  last_cube$input_checklist,
  "data/output/.*"
)

## Define prefix and URL to load the species list, then load it ####
species_list_prefix <- "https://raw.githubusercontent.com/guardias-eu/build-eu-cube/refs/heads/main/"
species_list_url <- paste0(species_list_prefix, last_cube_species_list)
species_list <- readr::read_csv(species_list_url, na = "")

# Calculate indicators #### 

# Get number of occurrences and measured occupancy (n grid cells) for each species ####

# Define function
calc_em_indicator <- function(cube, key) {
  # If key is not present in the cube, return empty data.frame
  if (!key %in% cube$specieskey) {
    return(tidyr::tibble(
      specieskey = integer(),
      year = integer(),
      n_occurrences = integer(),
      n_grid_cells = integer()
    ))
  }
  species_cube <- cube %>%
    dplyr::filter(specieskey == key) %>%
    dplyr::group_by(specieskey, year) %>%
    dplyr::summarise(
      n_occurrences = sum(occurrences),
      n_grid_cells = dplyr::n(),
      .groups = "drop"
    )
  # min, max year
  min_year <- min(species_cube$year, na.rm = TRUE)
  max_year <- max(species_cube$year, na.rm = TRUE)
  # Add 0s for years with no occurrences
  species_cube <- species_cube %>%
    tidyr::complete(
      year = min_year:max_year,
      fill = list(
        specieskey = key,
        n_occurrences = 0,
        n_grid_cells = 0)
    )
  # Transform to long version, with indicator as variable and value as value
  species_cube <- species_cube %>%
    tidyr::pivot_longer(
      cols = c(n_occurrences, n_grid_cells),
      names_to = "indicator",
      values_to = "value"
    )
  return(species_cube)
}

# Apply function to all species in the species list and all LMEs ####

species_keys <- species_list$usageKey[1:50]
names(species_keys) <- species_keys
# Calculate indicators for all species
indicators_list <- purrr::imap(
  lme_ids,
  function(lme_id, lme_name) {
    message("Calcuating indicators for LME: ", lme_name, " (ID: ", lme_id, ")")
    # Get grid cells for this LME
    grid_cells <- grid_cells_with_lme_info %>%
      dplyr::filter(lme_id == !!lme_id) %>%
      dplyr::pull(cellCode)
    # Filter cube for these grid cells
    cube_lme <- cube %>%
      dplyr::filter(eeacellcode %in% grid_cells)
    # Calculate indicators for each species in this LME
    species_indicators <- purrr::map(
      species_keys,
      ~calc_em_indicator(cube = cube_lme, key = .),
      .progress = TRUE
    )
    names(species_indicators) <- names(species_keys)
    return(species_indicators)
  },
  .progress = TRUE
)

# Plot occurrences and occupancy as subplots next to each other for each species ####

# Function for plots
plot_indicators <- function(df, species_key) {
  if (nrow(df) == 0) return(NULL)
  else {
    df %>%
    ggplot(aes(x = year, y = value)) +
      geom_point() +
      geom_line() +
    facet_wrap(~indicator, scales = "free_y", ncol = 1) +
    labs(
      title = paste("Species key:", species_key, ". n_occurrences (top), n_grid_cells (bottom)"),
      x = "Year",
      y = "Value"
    ) +
    theme_minimal()
  }
}

# Generate and save plots for each species ####

# Delete existing plots in output folder to avoid confusion with old plots when saving new ones. It can be that we are creating indicators 
# for a new cube, where the species list has changed, so we want to make sure that old plots are not mixed with new ones.
output_folder <- here::here("data/output/indicators_plots/indicators_plots_png/")
existing_plots <- list.files(output_folder, pattern = "\\.png$", full.names = TRUE)
if (length(existing_plots) > 0) {
  message("Deleting existing plots in output folder: ", output_folder)
  file.remove(existing_plots)
}

# Generate plots for each species in each LME. The resulting list has the structure: list(LME_name = list(species_key = plot, ...), ...)
plot_list <- purrr::imap(
  indicators_list,
  function(lme_id, lme_name) {
    message("Plotting indicators for LME: ", lme_name, " (ID: ", lme_ids[names(lme_ids) == lme_name], ")")
    species_plots <- purrr::imap(
      lme_id,
      function(species_indicators, species_key) {
        if (is.null(species_indicators)) return(NULL)
        plot_indicators(species_indicators, species_key)
      },
      .progress = TRUE
    )
  },
  .progress = TRUE
)

# Save plots as .png in output folder
purrr::iwalk(
  plot_list,
  function(lme_plots, lme_name) {
    message("Saving png plots for LME: ", lme_name, " (ID: ", lme_ids[names(lme_ids) == lme_name], ")")
    purrr::iwalk(
      lme_plots,
      function(p, s) {
        if (!is.null(p)) {
          ggsave(
            filename = paste0("lme_", lme_name, "_species_", s, ".png"),
            plot = p,
            path = output_folder,
            width = 12,
            height = 6
          )
        }
      },
      .progress = TRUE
    )
  }
)
