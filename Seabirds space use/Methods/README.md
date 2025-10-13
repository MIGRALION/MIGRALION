# Migralion_ :bird: :hatched_chick:

## Instructions to run the different scripts :memo:
All scripts use functions implemented in `2.code/functions.R`.

### Steps to prepare the data :building_construction:

1. Run `2.code/prepare_environmental_data.R` to generate the prediction grid `1.data/covariates.rdata` using downloaded environmental covariates.
2. Run `2.code/prepare_colonies.R` to filter the colony data and generate `1.data/colonies.rdata` (hidden due to data confidentiality).
3. Run `2.code/prepare_count_data.R` to aggregate and filter the count data from multiple surveys, producing `1.data/all_seabirds_counts.rdata` (hidden due to data confidentiality).
4. Run `2.code/prepare_telemetry_data.R` to  filter the GPS data and generate `1.data/all_seabirds_telemetry.rdata` (hidden due to data confidentiality).

### Steps to run the models :computer:
1. Run `2.code/run_Nmix.R` to fit N-mixture models for all seabird species during the breeding or wintering season (with or without colony covariates). This will generate results saved as `3.results/species/species_season_nmix.rdata`. 
2. Run `2.code/run_RSF.R` to fit RSF models for four seabirds species during the breeding or wintering season (with or without colony covariates). This will generate results saved as `3.results/species/species_season_rsf.rdata`.
3. Run `2.code/run_integrated.R` to fit integrated models for four seabird species during the breeding or wintering season (with or without colony covariates). This will generate results saved as `3.results/species/species_season_int.rdata`.

### Create reports :chart_with_upwards_trend:
1. Run `4.reports/Nmixture.rmd` to generate `Nmixture.html` summarizing all the results from the N-mixture models.
