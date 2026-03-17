# Load packages
library(dplyr)
library(mregions2)
library(sf)


# Load the 10x10km grid of Europe from GBIF Custom Downloads
grid_file <- "https://download.gbif.org/grids/EEA/EEA-Reference-Grid-10km.gpkg"
grid <- sf::st_read(grid_file)

# Get the Large Marine Ecosystems (LME) layer
lme <- mregions2::mrp_get(layer = "lme")

# LME we are interested in
lme_ids_eu <- c(1,3,6,11,13,14,50,53,54)
lme_eu <- lme[lme$objectid %in% lme_ids_eu, ]

# Validate polygons
lme_eu <- sf::st_make_valid(lme_eu)

# Transform the LMEs to the same CRS as the grid (ETRS89 / LAEA Europe, EPSG:3035)
lme_eu <- sf::st_transform(lme_eu, crs = sf::st_crs(grid))

# For each lme region, get the grid cells that intersect or fall within it
grid_cells_in_lmes <- sf::st_intersects(lme_eu, grid)
names(grid_cells_in_lmes) <- lme_eu$lme_name

# Add the `lme_id` and `lme_name` columns to `grid`.
grid <- grid %>% 
  dplyr::mutate(cell_row = row_number()) %>%
  dplyr::left_join(
    data.frame(
      cell_row = unlist(grid_cells_in_lmes),
      lme_id = rep(lme_eu$objectid, sapply(grid_cells_in_lmes, length)),
      lme_name = rep(lme_eu$lme_name, sapply(grid_cells_in_lmes, length))
    ),
    by = "cell_row"
  ) %>%
   dplyr::select(-cell_row)

# Save grid cellcodes with LME information as csv file. Only keep cells that are assigned to an LME
# and remove geometry for easier handling in R later on.
output_file <- "data/output/grid_cells_with_lme_info.csv"
grid %>%
  dplyr::filter(!is.na(lme_id)) %>%
  st_set_geometry(NULL) %>%
  readr::write_csv(output_file)
