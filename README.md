# TCR_diversity
Functions for preprocessing and downstream analysis of TCR sequencing data obtained with Adaptive Biotechnology.
The codes are specific for Dean et al dataset.

## PreProcessing
doing_PreProcessing is a function which allows for dataset filtration and curation based on rules provided in PreProcessing_human function. It also calculates Status Diversity and Sequences Diversity for Functional and NonFunctional set of sequences. 
PreProcessing_human is afunction which allows for dataset filtration, focusing on functional genes only.


## Subsampling for sample reconstitution regarding CD4/CD8 ratio
Samples obtained separately for CD4+ and CD8+ subsets of cells can be reconstituted to reflect the distribution of sequences obtained for blood sample with uncategorized CD4+ and CD8+ cells.
downsampling functional allows for the proprocessing and calclulation of Status/Sequence Diversity for the data after reconstitution.

## Linear regression modelling
regression_models_crosssectional function allows for constructing weighted regression models as well as segmented regression models on Dean et al dataset.
mixed-effectsModels are created to work with the longitudinal dataset results.

