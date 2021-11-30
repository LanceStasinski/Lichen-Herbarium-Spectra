# Lichen-Herbarium-Spectra

This repository contains R scripts, models, and figures constructed for my second thesis chapter called *Herbarium specimens: A source for building spectral libraries for lichens*.

## Main Scripts

* 01_clean_spectra.R - Cleans the raw spectra sorted by species in the `spectra` folder.

* 02_join_metadata.R - Adds taxonomic information to the metadata available in the `metadata` folder.

* 03_add_metadata.R - Adds the metadata to the spectra.

* 04_add_age.R - Adds specimen age to the metadata attached to the spectra.

* 05_plot_spectra.R - Plots the spectra in various configurations.

* 06_plsda_function.R - Function used to run PLS-DA on the spectra.

* 07_Classification.R - Classifies spectra into various taxonomic ranks using the function from `06_plsda_function`.

* 08_PLSDA_t-test.R - Compares accuracy metrics for the two PLS-DA model types for each taxonomic rank.

* 09_Compare_linear_mixed_models - Compares intercept-only, fixed slope - random intercept, and random slope - random intercept models using AIC and BIC for reflectance as a function of age.

* 10_LMM_taxon_as_random_effect - Computes linear mixed-effects model metrics for best model chosen by `09_Compare_linear_mixed_models`.

* 11_LMM_brightness_vs_age - Compares intercept-only, fixed slope - random intercept, and random slope - random intercept models using AIC for vector magnitude as a function of age and computes associated metrics.

## Dependencies

* caret
* corrplot
* dplyr
* lme4
* matrixStats
* mlbench
* naniar
* nlme
* optimx
* partR2
* rlist
* spectrolab
