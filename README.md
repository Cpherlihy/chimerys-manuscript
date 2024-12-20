# Unifying the analysis of bottom-up proteomics data with CHIMERYS

[Link to the bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2024.05.27.596040v2.full)

## Abstract

Proteomic workflows generate vastly complex peptide mixtures that are analyzed by liquid chromatography-tandem mass spectrometry (LC-MS/MS), creating thousands of spectra, most of which are chimeric and contain fragment ions from more than one peptide. Because of differences in data acquisition strategies such as data-dependent (DDA), data-independent (DIA) or parallel reaction monitoring (PRM), separate software packages employing different analysis concepts are used for peptide identification and quantification, even though the underlying information is principally the same. Here, we introduce CHIMERYS, a novel, spectrum-centric search algorithm designed for the deconvolution of chimeric spectra that unifies proteomic data analysis. Using accurate predictions of peptide retention time, fragment ion intensities and applying regularized linear regression, it explains as much fragment ion intensity as possible with as few peptides as possible. Together with rigorous false discovery rate control, CHIMERYS accurately identifies and quantifies multiple peptides per tandem mass spectrum in DDA, DIA and PRM experiments.

## Information

This repo describes how the data analysis and plots for all figures included in the manuscript were generated. To recreate the figures, make sure to download all input files (available on PRIDE via [PXD053241](https://www.ebi.ac.uk/pride/archive?keyword=PXD053241)), place them under `dataPath` (adjust in [load-dependencies.R](scripts/load-dependencies.R) to your own folder structure) and generate intermediate results in the linked `.R` scripts.
