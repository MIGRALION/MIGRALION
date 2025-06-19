# Migralion :bird: :hatched_chick:

🌊 Welcome to the MIGRALION Work Package 6 GitHub! 

🪶 This repository hosts the tools, data, and analyses developed to understand how seabirds and migratory birds use the Gulf of Lion — a crucial step toward ecologically responsible planning of floating wind farms in the Mediterranean.

:mortar_board: Our work Package aimed at integrating multi-source of counts, radar and tracking data to uncover the spatiotemporal dynamics of avifauna in the Gulf of Lion.

## Project introduction 

The Mediterranean Sea is a vital biodiversity hotspot and migratory corridor, but the Gulf of Lion is under growing pressure from human activity and climate change. Plans for floating wind farms in this region raise concerns about their ecological impact, especially on seabirds and migratory birds. Understanding how these species use the Gulf in space and time is crucial for sustainable planning. MIGRALION work packages 3–5 deployed complementary technologies to collect bird movement data, but each has limitations on its own. Work package 6 aimed to integrate these datasets through advanced statistical methods to provide a comprehensive, four-dimensional picture of avifauna use in the Gulf of Lion.

## Three different sub-projects :memo:

1.	Seabirds space use : how do seabirds use the marine environment of the Gulf of Lion?
2.	Terrestrial birds migratory flux : which areas in the Gulf of Lion have the most intense migratory bird flows? 
3.	Terrestrial birds flight hight : at what altitude do land migrants fly when crossing the Gulf of Lion?
   
## Three different sub-projects :memo:

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
