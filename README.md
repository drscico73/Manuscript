# Allele Frequency and Selection Strength Analysis

## Overview
_AFnS_plots.R_ script processes pooled sequencing allele frequency data for the _Directional crosses_ to:

* Clean and restructure sample metadata
* Calculate logit-transformed allele frequencies
* Estimate selection strength between generations (F1 → F10)
* Perform genome-wide statistical comparisons between experimental populations
* Generate publication-style plots for:
  * Allele frequency trajectories
  * Selection strength across the genome
  * Statistical significance of selection differences

## Required Input
The script expects an `.rds` file containing allele frequency data:
```r
df <- readRDS("/Path/to/file/AlleleFreq.rds")
```
The dataset should include:
* Chromosome (`chr`)
* Genomic position (`pos`)
* Sample identifiers (`sample`)
* Allele frequency estimates (`Dmel_OregonR_non-inbred`)

## Main Analysis Steps

### 1. Data Cleaning
Sample metadata are extracted from the `sample` column using delimiter-based parsing. Chromosomes are ordered as:
```r
X, 2L, 2R, 3L, 3R
```
### 2. Logit Transformation
Allele frequencies are transformed using the logit function to stabilize variance and enable estimation of selection response.

### 3. Selection Strength Estimation
Selection strength (`s`) is calculated as:
```r
s = (logit(F10) - logit(F1)) / 10
```
where 10 represents the number of generations between F1 and F10.

### 4. Statistical Testing
A genome-wide t-test compares selection responses between populations (`Cross`) at each genomic position. The threshold for significance is determined based on the _threshold_ simulations _(See threshold_for_ttest.R)_ 

## Output Plots
### Allele Frequency Plot (`p1`)
Displays median allele frequencies across genomic positions for each population and generation.

### Selection Strength Plot (`p2`)
Shows estimated selection strength (`s`) across chromosomes with mean ± SD ribbons.

### t-test Significance Plot (`p3`)
Visualizes genome-wide statistical significance, with significant loci highlighted.


