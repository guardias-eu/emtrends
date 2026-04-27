# Utility functions for the emerging status indicators

#' Function to summarize the number of occurrences and grid cells per year for a given species key
#' 
#' This function filters the input cube for the specified species key and 
#' summarizes the number of occurrences and grid cells per year. It then 
#' ensures that the resulting data frame includes all years between the 
#' minimum and maximum year, filling in 0s for any years with no occurrences.
#'
#' @param cube A data frame containing the species occurrence data with columns `specieskey`, `year`, and `occurrences`.
#' @param key The species key for which to get the timeseries data.
#' @return A data frame with columns `specieskey`, `year`, `n_occurrences`, 
#' and `n_grid_cells`, including all years between the minimum and maximum 
#' year for the specified species key.
summarise_species_timeseries <- function(cube, key) {
  # Filter the cube for the specified species key and summarize occurrences 
  # and grid cells per year
  species_timeseries <- cube %>%
    dplyr::filter(specieskey == key) %>%
    dplyr::group_by(specieskey, year) %>%
    dplyr::summarise(
      n_occurrences = sum(occurrences),
      n_grid_cells = dplyr::n(),
      .groups = "drop"
    )
  # Get minimum and maximum year
  min_year <- min(species_timeseries$year, na.rm = TRUE)
  max_year <- max(species_timeseries$year, na.rm = TRUE)
  # Add 0s for years between `min_year` and `max_year` with no occurrences
  species_timeseries %>%
    tidyr::complete(
      year = min_year:max_year,
      fill = list(
        specieskey = key,
        n_occurrences = 0,
        n_grid_cells = 0)
    )
}

#' Function to calculate the emerging status trend for a given variable and label
#' 
#' This function applies the GAM model to calculate the emerging status trend for a specified variable and label. If the GAM model cannot assess the emergence status, it falls back to using decision rules to determine the emergence status for each evaluation year. The resulting plot is adapted to reflect the emergence status based on decision rules when necessary.
#' @param v The variable for which to calculate the emerging status trend.
#' @param y_label The label for the y-axis in the resulting plot.
calc_em_trend <- function(v, y_label, species_ts, eval_years, min_year, max_year, key) {
    gam_output <- trias::apply_gam(
      df = species_ts,
      y_var = v,
      eval_years = eval_years,
      year = "year",
      taxonKey = "specieskey",
      type_indicator = "observations",
      baseline_var = NULL,
      x_label = "year",
      y_label = y_label,
      status_warning = FALSE
    )
    p <- gam_output$plot
    if (all(is.na(gam_output$em_summary$em_status))) {

      dr_output <- purrr::map_dfr(
        min_year:max_year,
        function(year) {
          trias::apply_decision_rules(
            df = species_ts,
            y_var = v,
            eval_year = year,
            year = "year",
            taxonKey = "specieskey"
          )
        }
      )
      dr_output <- dr_output %>% dplyr::mutate(model = "decision rules")
      # Adapt the plot with the colors based on `em_status` from decision rules output
      p <- adapt_plot(p, dr_output)
      trend_df <- dr_output
    } else {
      # Add model name to gam summary as well
      gam_output$em_summary <- gam_output$em_summary %>%
        dplyr::mutate(
          model = "GAM"
        )
      trend_df <- gam_output$em_summary
    }
    # Add model column to trend_occ and relocate it after `em_status`
    trend_df <- trend_df %>% dplyr::relocate(model, .after = em_status)
    # Adapt title of the plot to include species and model used
    p <- p + 
      ggplot2::ggtitle(paste0(
        "Species key: ", key, ". Model: ", unique(trend_df$model)
    ))
    # Bundle trend_df and the plot from GAM output in a list
    trend_output <- list(
      em_summary = trend_df,
      plot = p
    )
  }

#'  Colors mapping for emergence status
colors_mapping <- data.frame(
  em_status = c(c(3:0)),
  color = c("darkred", "orangered", "grey50", "darkgreen"),
  labels = c(
    "3" = "emerging (3)",
    "2" = "pot. emerging (2)",
    "1" = "unclear (1)",
    "0" = "not emerging (0)"
  )
)

#'  Function to adapt the plot with the colors based on `em_status`
#' 
#' This function is especially useful to adapt the plot created by 
#' `apply_gam()` where GAM modelling couldn't assess the emergence status.
#' Such plot gets the emergence status colors for evaluation years based on decision rules.
#' @param p A ggplot object, created by `apply_gam()`.
#' @param df A data frame that contains the `year` and `em_status` columns based on decision rules.
#' @param colors A data frame that contains the mapping of `em_status` to colors and labels.
#' @return A ggplot object with points colored according to `em_status`.
adapt_plot <- function(p, df, colors = colors_mapping) {
  # Join with plot data
  plot_data <- p$data %>%
    dplyr::left_join(df %>% select(year, em_status), by = "year")
  # Create the plot
  p + 
    ggplot2::geom_point(
      data = plot_data,
      aes(color = factor(em_status)),
      size = 3
    ) +
    ggplot2::scale_color_manual(
      values = colors$color, 
      breaks = colors$em_status,
      labels = colors$labels,
      na.value = "black"  # Set color for NA values
    ) +
    ggplot2::theme(legend.position = "right")
}
