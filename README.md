# Pesticide contamination in honey from organically managed colonies in the continental U.S. is comparable to certified organic imports

[![Zenodo DOI: 18321610](https://zenodo.org/badge/DOI/10.5281/zenodo.18321610.svg)](https://doi.org/10.5281/zenodo.18321610)

Beekeepers in the continental United States do not have the option of USDA organic certification due to unfeasible apiculture standards. Interim recommendations require more organic-certified land for foraging than is available surrounding most apiaries. To test whether typical organic farms in New York and Pennsylvania (USA) provide sufficient forage for organic honey production, we established 72 colonies using organic management practices on six organic farms in 2023. In 2024, we screened the resulting honey for 96 pesticide residues and tested whether landscape features were associated with pesticide contamination. Additionally, we compared the levels of pesticide contamination to 20 brands of store-bought organic and conventional honey.

## R Analysis

1.  We used R version 4.4.3.

2.  All data inputs required to run the analysis are located within the `data` directory.

3.  The R Markdown script, `Analysis.Rmd` generates all tables, figures, and statistics referenced in the manuscript. Results from this script can be viewed in either `Analysis.html` or `Analysis.md`.

4.  The `functions` directory contains some R scripts referenced by `Analysis.Rmd`. This is to keep code inside `Analysis.Rmd` concise.

5.  Running `Analysis.Rmd` creates a separate `cache` directory to store results of resource-intensive operations and save them for later use.

6.  The analysis script outputs figure images to the `figures` directory.

### Packages

| Package | Version | Purpose |
|------------------------|------------------------|------------------------|
| `tidyverse` | 2.0.0 | Easier data manipulation and visualization |
| `readxl` | 1.4.5 | Import Excel files directly into R |
| `sf` | 1.0-19 | Handle spatial data |
| `crawlUtils` | 0.1.62 | Provides function to expand `st_bbox` limits |
| `basemaps` | 0.0.8 | Downloads base map imagery for ggplot map |
| `ggspatial` | 1.1.9 | Provides legend scale and north arrow for map |
| `ggrepel` | 0.9.6 | Provides `geom_text_repel` function |
| `khroma` | 1.16.0 | Colors for heatmap |
| `patchwork` | 1.3.0 | Arrange plots in grids |
| `ggedit` | 0.4.2 | Remove a `geom` from a ggplot for visual abstract |
| `gt` | 1.0.0 | Formatted tables |
| `glmmTMB` | 1.1.10 | Generalized linear mixed effects models |
| `performance` | 1.13.0 | Provides `check_overdispersion` function |
| `emmeans` | 1.10.7 | Pairwise tests |
| `indicspecies` | 1.8.0 | Indicator Species Analysis |
| `vegan` | 2.6-10 | For ordinations |
| `ggordiplots` | 0.4.4 | Generates hulls for ordination plots |
| `colorspace` | 2.1-1 | Provides `darken` color function |
| `pairwiseAdonis` | 0.4.1 | Pairwise PERMANOVA |
| `units` | 0.8-6 | Support units in R vectors, for `set_units()` |
| `CropScapeR` | 1.1.5 | Download landscape data from the USDA Cropland Data Layer |
| `MuMIn` | 1.48.11 | AICC comparison for landscape models |
| `ggeffects` | 2.2.1 | GLMM predictions |
