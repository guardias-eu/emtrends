# emtrends

Thie repo contains the workflows to generate emerging trends for Europe's Large Marine Ecosystems. Such workflows are based on [species occurrence cubes](https://www.gbif.org/occurrence-cubes), a concept developed during the [TrIAS](https://osf.io/7dpgr/) project and made operational by [GBIF](https://www.gbif.org/) and the [B-Cubed](https://b-cubed.eu/) project.

The emerging trends indicators calculated in this repository have been conceived and applied at Belgian and regional level during the TrIAS project, see occurrence-based [indicators](https://trias-project.github.io/indicators/).

## Repo structure

```
├── src/
│   ├── assign_grid_cells_to_lme.R   # Assigns EEA 10×10 km grid cells to Large Marine Ecosystems (LMEs)
│   ├── calculate_indicators.R       # Downloads the latest species occurrence cube and calculates emerging-trend indicators per LME
│   └── utils.R                      # Shared utility functions used by the scripts above
├── data/
│   ├── input/                       # Species occurrence cubes downloaded from GBIF (CSV files)
│   └── output/
│       ├── grid_cells_with_lme_info.csv        # Mapping of EEA grid cell codes to LME identifiers
│       ├── species_lme_combinations.csv        # All species × LME combinations present in the cube
│       ├── emerging_trends_summary.csv         # Per-species/LME emerging-trend classification and statistics
│       ├── emerging_trends_ranking_list.csv    # Ranked list of emerging species per LME
│       ├── appearing_species.csv               # Species newly appearing in an LME
│       ├── reappearing_species.csv             # Species reappearing in an LME after absence
│       └── indicators_plots/                   # Indicator plots per species and LME (zip archives and PNG files)
└── .github/
    └── workflows/
        └── calculate-indicators-and-create-pr.yml  # GitHub Actions workflow that runs the full pipeline and opens a PR with the results
```

## Funding

This work is being developed in the framework of the [GuardIAS](https://guardias.eu/) prject. GuardIAS receives funding from the European Union’s Horizon Europe Research and Innovation Programme (ID No 101181413).
