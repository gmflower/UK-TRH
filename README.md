# Extreme heat and cause-specific risk of hospital admission in the adult population in England: a case time series analysis

Add DOI when available.

Partially reproducible code performing the analysis reported in the paper:

> Add full citation when available

## Scripts

The scripts are run in order. Script 01 should be run first to set up the
packages and parameters. Scripts 02 and 03 prepare the data and run the first stage
analysis, these require the HES data and cannot be reproduced here. Scripts 04 to 05 
can be run by loading results of the first stage saved beforehand. This reproduced
the main results from the study. Scripts 06 and 07 contain the sensitivity analyses
and can be reproduced from the second stage (line 166).

| Script | Descriptions |
| :--- | :--- |
| `01.pkg_param.R` | Load the necessary R libraries and defines all the analysis parameters. |
<<<<<<< HEAD
| `02.prepmean.R` | Load and prepare the environmental and hospital admissions data |
=======
| `02.prepmain.R` | Load and prepare the environmental and hospital admissions data |
>>>>>>> a0bb80198ef14fdce1f4d1c2a50cf18f3737edbc
| `03.firststage.R` | Centrepiece of the analysis. Loops through the causes and age groups to perform the case time series and produce the LAD level ERF.| 
|`04.secondstage.R` | Runs the meta analysis pooling the results to obtain a single national ERF.|
| `05.plots.R` | Produces plots featured in the main text of the article and supplementary plots related to the main analyses.|
| `06.sensitivity_knot.R` | Reproduces the stages in scripts 03-05 for the sensitivity analysis on knot placement.|
| `07.sensitivity_lag.R` | Reproduces the stages in scripts 03-05 for the sensitivity analysis on lag length.|
