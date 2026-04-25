# Load libraries
library(dplyr)
library(ggplot2)
library(here)
library(lubridate)
library(purrr)
library(readr)
library(rgbif)
library(stringr)
library(tidyr)
library(trias)
library(zip)

# Source utility functions
source(here::here("src/utils.R"))

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

# Define evaluation years ####

# For emerging trends, take two years before the current year as last year.
# E.g. if 2026, take 2024 as last year for emerging scores, to avoid issues with incomplete data 
# due to data publication delays.
last_eval_year <- lubridate::year(Sys.Date()) - 2
first_eval_year <- last_eval_year - 2
eval_years <- first_eval_year:last_eval_year

# Get appearing taxa ####

# Appearing species are species that have their first occurrences in one of the years from first
# evaluation year up to now. To be done for each LME, based on lme_ids and 
# grid_cells_with_lme_info, to get the appearing taxa for each LME.
appearing_species <- purrr::imap(
  lme_ids,
  function(lme_id, lme_name) {
    message("Getting appearing species for LME: ", lme_name, " (ID: ", lme_id, ")")
    # Get grid cells for this LME
    grid_cells <- grid_cells_with_lme_info %>%
      dplyr::filter(lme_id == !!lme_id) %>%
      dplyr::pull(cellCode)
    # Filter cube for these grid cells and evaluation years
    cube_lme <- cube %>%
      dplyr::filter(eeacellcode %in% grid_cells)
    # Get appearing species for this LME
    appearing_species_lme <- cube_lme %>%
      dplyr::group_by(specieskey, species) %>%
      dplyr::summarise(
        first_year = min(year),
        .groups = "drop"
      ) %>%
      dplyr::filter(first_year >= first_eval_year) %>%
      dplyr::mutate(
        lme_id = lme_id,
        lme_name = lme_name
      ) %>%
      dplyr::relocate(lme_id, lme_name, specieskey, species, first_year)
    return(appearing_species_lme)
  },
  .progress = TRUE) %>%
  purrr:::list_rbind()

## Save appearing species as csv in output folder ####
readr::write_csv(
  appearing_species,
  here::here("data", "output", "appearing_species.csv"),
  na = ""
)

# Get reappearing species ####

# Reappearing species are species that have their first occurrences before the 
# first evaluation year, but then have no occurrences for at least two years, 
# and then reappear in one of the evaluation years. To be done for each LME, 
# based on lme_ids and grid_cells_with_lme_info, to get the reappearing species 
# for each LME.

# Define lowest amount of years without occurrences to consider a species as reappearing.
latency_threshold <- 5 
reappearing_species <- purrr::imap(
  lme_ids,
  function(lme_id, lme_name) {
    message("Getting reappearing species for LME: ", lme_name, " (ID: ", lme_id, ")")
    # Get grid cells for this LME
    grid_cells <- grid_cells_with_lme_info %>%
      dplyr::filter(lme_id == !!lme_id) %>%
      dplyr::pull(cellCode)
    # Filter cube for these grid cells
    cube_lme <- cube %>%
      dplyr::filter(eeacellcode %in% grid_cells)
    if (nrow(cube_lme) == 0) return(NULL)
    # Get reappearing species for this LME
    reappearing_species_lme <- cube_lme %>%
      dplyr::group_by(specieskey, species) %>%
      dplyr::summarise(
        # set to max(year) if there are occurrences in the evaluation years, otherwise set to NA. This way we can filter for species that reappear in the evaluation years.
        last_year = ifelse(
          any(year >= first_eval_year),
          min(year[year >= first_eval_year]),
          NA_integer_
        ),
        gap_from_last_year = ifelse(
          any(year >= first_eval_year) & any(year < first_eval_year),
          min(year[year >= first_eval_year]) - max(year[year < first_eval_year]),
          NA_integer_
        ),
        .groups = "drop") %>%
      dplyr::filter(!is.na(last_year) & !is.na(gap_from_last_year) & gap_from_last_year >= latency_threshold) %>%
      dplyr::mutate(
        lme_id = lme_id,
        lme_name = lme_name
      ) %>%
      dplyr::relocate(lme_id, lme_name, specieskey, species, last_year, gap_from_last_year) %>%
      dplyr::rename(
        reappearance_year = last_year,
        years_without_occurrences = gap_from_last_year
      )
    return(reappearing_species_lme)
  },
  .progress = TRUE) %>%
  purrr:::list_rbind()

## Save reappearing species as csv in output folder ####
readr::write_csv(
  reappearing_species,
  here::here("data", "output", "reappearing_species.csv"),
  na = ""
)

# Calculate emergence trends indicators based on GAM and decision rules #### 

# Define function
calc_em_indicator <- function(cube, key) {
  # If `key` is not present in the `cube`, return `NULL`
  if (!key %in% cube$specieskey) {
    return(NULL)
  }
  species_cube <- cube %>%
    dplyr::filter(specieskey == key) %>%
    dplyr::group_by(specieskey, year) %>%
    dplyr::summarise(
      n_occurrences = sum(occurrences),
      n_grid_cells = dplyr::n(),
      .groups = "drop"
    )
  # Get minimum and maximum year
  min_year <- min(species_cube$year, na.rm = TRUE)
  max_year <- max(species_cube$year, na.rm = TRUE)
  # Add 0s for years between `min_year` and `max_year` with no occurrences
  species_cube <- species_cube %>%
    tidyr::complete(
      year = min_year:max_year,
      fill = list(
        specieskey = key,
        n_occurrences = 0,
        n_grid_cells = 0)
    )

  # Apply GAM modelling or decision rules to calculate emerging trends indicators.
  # Do it for both `n_occurrences` and `n_grid_cells`.
  if (!all(eval_years %in% species_cube$year)) {
    return(NULL)
  }
  variable <- list(
    "number of occurrences" = "n_occurrences",
    "number of grid cells (10x10km)" = "n_grid_cells"
  )
  purrr::imap(
    variable,
    calc_em_trend, # Defined in src/utils.R
    species_cube = species_cube,
    eval_years = eval_years,
    min_year = min_year,
    max_year = max_year,
    key = key
  )
}

# Apply function to all species in the species list and all LMEs ####
n_species <- nrow(species_list)
species_keys <- unique(species_list$usageKey[1:n_species])
names(species_keys) <- species_keys
species_names <- unique(species_list$canonicalName[match(species_keys, species_list$usageKey)])
names(species_names) <- species_keys
if (length(species_names) < length(species_keys)) {
  stop(paste(
    "Species names are not unique, please check the species list.",
    "Consider using the `scientificName` column instead of `canonicalName` if there are issues with non-unique names."
  ))
}
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

# Remove NULLs from each element in indicators_list, to avoid confusion 
# when saving the plots. NULL means that there were no occurrences for that 
# species in that LME, so no indicators (and so plots) were generated.
indicators_list <- purrr::map(
  indicators_list,
  function(ind) {
    purrr::compact(ind)
  }
)


# Save plots for each species ####
# Delete existing plots in output folder to avoid confusion with old plots when saving new ones. It can be that we are creating indicators 
# for a new cube, where the species list has changed, so we want to make sure that old plots are not mixed with new ones.

output_png_folder <- here::here("data/output/indicators_plots/png/")
existing_plots <- list.files(output_png_folder, pattern = "\\.png$", full.names = TRUE)
if (length(existing_plots) > 0) {
  message("Deleting ", length(existing_plots), " PNG plots in output folder: ", output_png_folder)
  file.remove(existing_plots)
}

# Save plots as .png in output folder
purrr::iwalk(
  indicators_list,
  function(lme_indicators, lme_name) {
    message("Saving png plots for LME: ", lme_name, " (ID: ", lme_ids[names(lme_ids) == lme_name], ")")
    purrr::iwalk(
      lme_indicators,
      function(ind, s) {
        purrr::imap(
          ind,
          function(trend_output, v) {
            p <- trend_output$plot
            if (v == "number of occurrences") {
              v <- "occs"
            } else if (v == "number of grid cells (10x10km)") {
              v <- "grid_cells"
            }
            # Save the plot as .png in output folder
             ggsave(
              filename = paste0("lme_", lme_name, "_species_", s, "_indicator_", v, ".png"),
              plot = p,
              path = output_png_folder,
              width = 12,
              height = 6,
              create.dir = TRUE
            )
          }
        )
      }
    )
  }, .progress = TRUE
)

# Save plots as ggplot2 obejcts in zip files in output folder. This allows us to test
# the reactivity of OJS to transform ggplot2 objects into plotly objects.
# One zip file per each LME or maximum of 100 species, containing the ggplot2 objects for all species in that LME.
# The structure of the zip file is: list(species_key = plot, ...)

# Delete existing plots in output folder to avoid confusion with old plots when saving new ones. It can be that we are creating indicators 
# for a new cube, where the species list has changed, so we want to make sure that old plots are not mixed with new ones.

output_ggplot_folder <- here::here("data/output/indicators_plots/ggplot")
existing_plots <- list.files(output_ggplot_folder, pattern = "\\.zip$", full.names = TRUE, recursive  = FALSE)
if (length(existing_plots) > 0) {
  message("Deleting ", length(existing_plots), " zip files plots in output folder: ", output_ggplot_folder)
  file.remove(existing_plots)
}

message("Save ggplot2 objects as zip files")
purrr::iwalk(
  indicators_list,
  function(lme_indicators, lme_name) {
    message("Saving ggplot2 plots for LME: ", lme_name, " (ID: ", lme_ids[names(lme_ids) == lme_name], ")")
    if (length(lme_indicators) == 0) return(NULL)
    # Split plots into chunks of 50 species and save each chunk as a separate zip file
    # Extract the ggplot2 objects for each species in this LME
    lme_plots <- purrr::imap(
      lme_indicators,
      function(ind, s) {
        purrr::imap(
          ind,
          function(trend_output, v) {
            p <- trend_output$plot
            if (v == "number of occurrences") {
              v <- "occs"
            } else if (v == "number of grid cells (10x10km)") {
              v <- "grid_cells"
            }
            return(p)
          }
        )
      }
    )
    plot_list_chunks <- split(lme_plots, ceiling(seq_along(lme_plots) / 50))
    purrr::imap(
      plot_list_chunks,
      function(chunk, i) {
        plot_list_zip_file <- here::here("data", "output", "indicators_plots", paste0("indicators_plots_ggplot2_", lme_name, "_chunk_", i, ".zip"))
        plot_list_rdata_file <- here::here("data", "output", "indicators_plots", paste0("indicators_plots_ggplot2_", lme_name, "_chunk_", i, ".RData"))
        save(chunk, file = plot_list_rdata_file)
        zip::zip(zipfile = plot_list_zip_file, files = plot_list_rdata_file, mode = "cherry-pick")
        file.remove(plot_list_rdata_file)
      }
    )
  },
.progress = TRUE
)


# Create a summary of all emerging indicators in a dataframe, to be able to 
# create a ranking list of species based on their emerging status in each LME.
# We have to extract the emerging status for each species, variable, yaer and LME, based on GAM or decision rules.
# Take into accout only the years in the evaluation years.
# This allow us to create a ranking list of species based on their emerging status in each LME.
em_df <- purrr::imap_dfr(
  indicators_list,
  function(lme_indicators, lme_name) {
    purrr::imap_dfr(
      lme_indicators,
      function(species_indicator, species_key) {
        if (!is.null(species_indicator)) {
          purrr::imap_dfr(
            species_indicator,
            function(trend_output, variable) {
              trend_output$em_summary %>%
                dplyr::mutate(
                  species_name = species_names[as.character(species_key)],
                  lme_id = as.integer(lme_ids[names(lme_ids) == lme_name]),
                  lme_name = lme_name,
                  variable = variable,
                  # if growth exists, hold it, otherwise set it to NA. Growth is only available for GAM, not for decision rules, but we want to have the column in the dataframe for both models to be able to compare them.
                  growth = ifelse("growth" %in% colnames(trend_output$em_summary), growth, NA_real_),
                  em_status = as.integer(em_status)
                ) %>%
                dplyr::select(lme_id, lme_name, specieskey, species_name, year, variable, model, em_status, growth)
            }
          )
        }
      }
    )
  }
) %>%
  dplyr::filter(year %in% eval_years) %>%
  dplyr::mutate(year = as.integer(year))

# Save emerging trends summary dataframe as csv in output folder
readr::write_csv(
  em_df,
  here::here("data/output/emerging_trends_summary.csv"),
  na = ""
)

# Assign  a weigth for each year/variable combination, to be able to calculate a weighted emerging status for each species in each LME.
# We give more weight to the most recent year and to the variable "number of grid cells" compared to "number of occurrences"
# because it is a more robust indicator of emerging trends, less affected by sampling effort.
weights <- dplyr::tibble(
  year = rep(eval_years, each = 2),
  variable = rep(c("number of occurrences", "number of grid cells (10x10km)"), times = length(eval_years)),
  weights = c(1, 1.5, 1.5, 2, 2, 2.5)
)

# Apply weigths to the emerging status and calculate a weighted emerging status for each species in each LME.
# Also, in case equal weighted emerging status, higher rank to species using GAM and use the growth column to give higher rank to species that have a  higher growth, as they are more likely to be emerging.
em_rank <- em_df %>%
  dplyr::left_join(weights, by = c("year", "variable")) %>%
  dplyr::mutate(weighted_em_status = em_status * weights) %>%
  dplyr::group_by(lme_id, lme_name, specieskey, species_name) %>%
  dplyr::summarise(
    weighted_em_status = sum(weighted_em_status, na.rm = TRUE),
    model = ifelse(any(model == "GAM"), "GAM", "decision rules"),
    growth = ifelse(any(model == "GAM"), growth[model == "GAM"][1], NA_real_),
    .groups = "drop"
  ) %>%
  dplyr::arrange(lme_id, desc(weighted_em_status), model, desc(growth))

# Save emerging trends ranking list as csv in output folder
readr::write_csv(
  em_rank,
  here::here("data/output/emerging_trends_ranking_list.csv"),
  na = ""
)

# For the dashboard, it's handy to have a csv with species keys, species names,
# LME IDs en LME names. Only existing combinations, e.g. not NULL plots.
species_lme_combinations <- purrr::imap_dfr(
  indicators_list,
  function(lme_indicators, lme_name) {
    purrr::imap_dfr(
      lme_indicators,
      function(species_indicator, species_key) {
        if (!is.null(species_indicator)) {
          tidyr::tibble(
            species_key = as.integer(species_key),
            species_name = species_names[[as.character(species_key)]],
            lme_id = as.integer(lme_ids[names(lme_ids) == lme_name]),
            lme_name = lme_name
          )
        }
      }
    )
  }
)
readr::write_csv(
  species_lme_combinations,
  here::here("data/output/species_lme_combinations.csv"),
  na = ""
)
