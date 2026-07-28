# Selection Response Simulations and Allele Frequency Analysis

This repository contains R and SLiM scripts used to analyze allele frequency dynamics and simulate selection responses in experimental evolution populations.

The workflows include:

* Allele frequency estimation and visualization
* Selection response estimation
* Genome-wide statistical testing
* Forward-time SLiM simulations
* Empirical threshold estimation for genome-wide t-tests
* Testing for deviations from additive assumption

---

# Repository Structure

| File                        | Description                                                              |
| --------------------------- | ------------------------------------------------------------------------ |
| `AFnS_plots.R`              | Allele frequency and selection response analysis for directional crosses |
| `H0_Check.R`                | Single-chromosome SLiM simulations for directional crosses               |
| `H0_Check2Chrom.R`          | Two-chromosome SLiM simulations for directional crosses                  |
| `threshold_for_ttest.R`     | Empirical significance threshold estimation                              |
| `func_ee_analysis.R`       | Functions for analysis of experimental evolution data for 3x2 crosses              |
| `func_slim_analysis.R`      | Functions for analysis of SLiM simulations for 3x2 crosses               |
| `func_visual_ee_data.R`    | Functions for visualization of experimental evolution data for 3x2 crosses         |
| `func_visual_slim.R`        | Functions for visualization of SLiM simulation data for 3x2 crosses      |
| `main_script_ee_data.R`    | Main script for analysis of experimental evolution data for 3x2 crosses            |
| `main_script_slim_dom.R`    | Main script for dominance simulations in 3x2 crosses                     |
| `main_script_slim_one_chr.R`| Main script for single-chromosome SLiM simulations for 3x2 crosses       |
| `main_script_two_chr.R`     | Main script for two-chromosome SLiM simulations for 3x2 crosses          |
| `predef_win_haploFreq.R`    | Script for calculating allele frequencies using predefined windows       |
| `haploFreq.R`               | Script for calculating allele frequencies using step and window size     |
| `H0_Check.slim`             | SLiM model for single-chromosome simulations for directional crosses     |
| `H0_2Chrom.slim`            | SLiM model for two-chromosome simulations for directional crosses        |
| `Threshold.slim`            | SLiM model for null simulations                                          |
| `one_chr.slim`              | SLiM model for one-chromosome simulations for 3x2 crosses                |
| `two_chr.slim`              | SLiM model for two-chromosome simulations for 3x2 crosses                |
| `threshold_seeds.txt`       | Seeds for reproducible threshold simulations                             |
| `markers*.txt`              | Marker/window positions used in simulations                              |
| `ore_marker_file.txt`       | Marker SNPs for OregonR                                                  |
| `hel_marker_file.txt`       | Marker SNPs for Helsinki                                                 | 
| `sam_marker_file.txt`       | Marker SNPs for Samarkand                                                |
| `window_coordinates.txt`    | Predefined window coordinates for 3x2 crosses                            |

---

# Requirements

## R Packages

```r
data.table
dplyr
tidyr
ggplot2
ggpubr
scales
purrr
foreach
doParallel
poolSeq
boot
patchwork
cowplot
stringr
stats
emmeans
ggh4x
forcats
pROC
docopt
MASS
Matrix
matrixStats
quadprog
fuzzyjoin
```

## External Software

* SLiM (v5.01 or compatible)
* GNU command-line utilities:

  * `sed`
  * `cut`
  * `grep`

---

# 1. Allele Frequency and Selection Strength Analysis for Directional Crosses

## Script

`AFnS_plots.R`

## Overview

This script processes pooled sequencing allele frequency data from the directional crosses to:

* Clean and restructure sample metadata
* Calculate logit-transformed allele frequencies
* Estimate selection response between generations (F1 → F10)
* Perform genome-wide statistical comparisons between populations
* Generate publication-quality figures

## Required Input

```r
df <- readRDS("/Path/to/file/AlleleFreq.rds")
```

Required data:

* `chr` Chromosome
* `pos` Positions (Mean position of the window)
* `sample` The Cross direction - SO/OS
* `Dmel_OregonR_non-inbred` Oregon allele frequency

## Selection Response Calculation

```r
s = (logit(F10) - logit(F1)) / 10
```

## Statistical Testing

Genome-wide t-tests compare selection responses between populations (`Cross`) at each genomic position.

Significance thresholds are determined using simulations from:

```bash
threshold_for_ttest.R
```

## Output Plots
### `p1`
Allele frequency trajectories across genomic positions.
### `p2`
Selection response (`s`) across chromosomes.
### `p3`
Genome-wide significance plot showing statistically significant loci.

------

# 2. Single-Chromosome SLiM Simulations: Directional Crosses

## Scripts

* `H0_Check.R`
* `H0_Check.slim`

## Overview

These scripts perform forward-time SLiM simulations using SLiM to investigate how linkage distance, haplotype structure, and epistasis influence allele frequency change and selection response of two alleles associated with same chromsome (Specifically X) over time. The simulation configurations are controlled within the included SLiM script `H0_Check.slim`

## Simulation Parameters

| Parameter      | Description                     |
| -------------- | ------------------------------- |
| `popSize`      | Population size                 |
| `numReps`      | Number of replicate simulations |
| `target1`      | Reference selected locus        |
| `distances_mb` | Distance between selected loci  |
| `selCoefs`     | Selection coefficients          |
| `epis`         | Epistasis coefficient           |
| `domCoefs`     | Dominance coefficient           |

## Simulation Scenarios

The pipeline supports four evolutionary scenarios:

| Scenario | Epistasis                        | Haplotype Structure  | Required Modifications in SLiM script |
| -------- | -------------------------------- | -------------------- | --------------------------------------|
| I        | Epistatic         | Same haplotype       | Comment out the chunk from line `111 to 116` and use the chunk from line `120 to 125` of the SLiM script, set (`epis = 2.00`) | 
| II       | Neutral/Additive  | Same haplotype       | Comment out the chunk from line `111 to 116` and use the chunk from line `120 to 125` of the SLiM script, set (`epis = 1.00`) |
| III      | Epistatic     | Different haplotypes | Comment out the chunk from line `120 to 125` and use the chunk from line `111 to 116` of the SLiM script, set  (`epis = 2.00`)     |
| IV       | Neutral/Additive  | Different haplotypes | Comment out the chunk from line `120 to 125` and use the chunk from line `111 to 116` of the SLiM script, set (`epis = 1.00`) |

* The script needs to be run 4 times, once for each scenario by making necessary changes to the SLiM and/or R-script
  
## Important functions 
`initSlimSim()` creates an initial population state and writes a temporary binary population file.
`runSlimSim()` executes forward simulations to specified generations using SLiM.
`getSimFreqs2()` parses simulation output and calculates allele frequencies across replicates.

** These functions are common for all the SLiM simulations for directional cross scenarios, with minor changes in `getSimFreqs2()` based to the scenarios and chromosomes simulated.

## Output

The simulations generate:

* Selection response summaries
* Plots depicting the distance-dependent response

Simulation outputs are written as tab-delimited files.

---

# 3. Two-Chromosome SLiM Simulations: Directional Crosses

## Scripts

* `H0_Check2Chrom.R`
* `H0_2Chrom.slim`

## Overview

This workflow extends the previous simulations to selected loci located on separate chromosomes, enabling analysis of interchromosomal interactions, additive vs. epistatic fitness effects and chromosome-specific selection dynamics. The required SLiM script is `H0_2Chrom.slim` defining the evolutionary model used during simulations.

## Key Parameters

| Parameter            | Description            |
| -------------------- | ---------------------- |
| `popSize`            | Population size                 |
| `numReps`            | Number of replicate simulations |
| `target1`, `target2` | Selected loci          |
| `selCoefs`           | Selection coefficients          |
| `epis`               | Epistasis coefficient  |
| `domCoefs`           | Dominance coefficients |

## Simulation Scenarios
The script can be used to compare:

* Additive fitness effects (`epis = 1.00`)
  Alternatively, comment out the tick 2: late{} section in the SLiM script to run the additive scenario
  
* Epistatic fitness effects (`epis = 2.00`)


## Output
The script generates:
* Boxplots of selection response (`s`)
* Statistical comparisons between populations
* Additive vs. epistatic scenario comparisons

---

# 4. Empirical Significance Threshold Estimation

## Scripts

* `threshold_for_ttest.R`
* `Threshold.slim`

## Overview

This workflow estimates empirical significance thresholds for genome-wide t-tests using null simulations.

The simulations generate null distributions of allele frequency change across multiple chromosomes to evaluate the expected false-positive distribution of the t-tests. This proportion of the false positives is then used as as significance cutoff for the t-tests on empirical data.

The simulation uses the same window positions/marker positions as used in the allele frequency estimation to keep the number of windows consistent & test the multiple testing load due to the high number of windows.

The repository also includes the supporting files required:
* The marker `markers*.txt` files required for simulations and chromosome setup.
*  Seed files for reproducible simulations (`threshold_seeds.txt`)



## Population Setup

| Population | Effective Population Size | Initial Frequency |
| ---------- | ------------------------- | ----------------- |
| `OS`       | 243                       | 0.67              |
| `SO`       | 356                       | 0.33              |


## Threshold Estimation
The empirical significance cutoff is estimated as the 5th percentile of simulated p-values:
```r
quantile(all_results_df$p_value, 0.05)
```
Chromosome-specific cutoffs are also calculated separately.

## Output

The script generates:

* genome-wide p-value distributions,
* empirical p-value histograms,
* significance threshold estimates.


-----

# 5. Testing for Deviations from Additive Assumption: 3x2 Crosses

## Scripts

* `main_script_ee_data.R`: Main script for analysis of experimental evolution data
* `func_ee_analysis.R`: Functions for analysis of experimental evolution data
* `func_visual_ee_data.R`: Functions for visualization of experimental evolution data

## Overview

These scripts process pooled sequencing allele frequency data from 3x2 crosses to:

* Clean and restructure sample metadata
* Average allele-frequency estimates obtained from two independent marker sets
* Calculate logit-transformed allele frequencies as a measure of selection response
* Estimate deviations from additive expectations
* Perform genome-wide statistical test for deviations from additivity
* Correct for multiple testing using Benjamini-Hochberg procedure
* Classify genomic regions according to number of significant deviations
* Generate publication-quality figures

## Required Inputs

Six RDS files containing window-level allele-frequency estimates:

* hel_ho_out_500.rds
* ore_ho_out_500.rds
* hel_hs_out_500.rds
* sam_hs_out_500.rds
* ore_so_out_500.rds
* sam_so_out_500.rds

Required data fields include:

* `chr` Chromosome
* `win` Genomic window coordinates
* `sample` Sample/cross identity and experimental metadata
* `gen` Generation
* `rep` Replicate
*  Allele-frequency estimates for the corresponding parental/reference lines

## Data Processing

`process_data()`:

* Removes reference samples
* Parses sample, generation, and replicate information
* Standardizes sample and allele-frequency column names
* Calculates the midpoint `mpos` of each genomic window

For each pairwise comparison, allele frequencies estimated from two independent marker sets are averaged.

## Selection Response Calculation

Selection response is calculated with `logit_freq()` as the log ratio of allele frequencies: 

s= log(AF_A/AF_B)

## Additive Expectation and Deviation and Statistical Testing

`compute_p()`:

* Tests three additive expectations
* Correct for multiple testing
* Classify windows as significantly deviating (1) or not significantly deviating (0) from additive expectations
* Calculates deviation from additive expectations

Genome-wide statistical tests are performed independently for each genomic window using linear models: 

lm( s ~ group, data = .)

`emmeans` is used to estimate the expected contrast, which is compared with the observed selection response of the focal population.

FDR is estimated using the ***Benjamini-Hochberg*** method. Genomic windows with ***Adjusted P ≤ 0.05*** are classified as a significant deviation from additivity. 

The three additive expectations tested are: 

H(HO-OH) - S(SO-OS) = H(HS-SH) 

O(SO-OS) - H(HS-SH) = H(HO-OH) 

S(HS-SH) - O(HO-OH) = S(SO-OS)

The deviation from the additive expectation is calculated as: ***Deviation = Observed response - Expected response***
 

## Genomic Classification

Significance across the three additive tests is combined using `find_pattern()` to classify each window: 

|Classification|Description|
|------|-------
|0/3|No significant deviations comparisons|
|1/3|Significant deviations in one comparison|
|2/3|Significant deviations in two comparisons|
|3/3|Significant deviations in all three comparisons|

## Output Plots

### `af_comp`

Combined allele-frequency trajectories across genomic positions and chromosomes for the three pairwise comparisons.

### `plot_so_exp_obs`

Expected vs Observed selection response in SO-OS.

### `plot_class_res_10`

Genome-wide visualization of deviations from additive expectations, showing calculated versus expected logit allele-frequency differences and genomic regions classified according to the number of significant tests.

## Main Output Objects

### `merge_df_10`

Combined significance classifications across all three comparisons.

### `merge_dev_10`

Combined genome-wide deviation estimates across all three comparisons.

-----

# 6. Single-Chromosome Epistasis Simulations: 3x2 Crosses

## Scripts

* `main_script_slim_one_chr.R`
* `func_slim_analysis.R`
* `func_visual_slim.R`
* `one_chr.slim`

## Overview

These scripts perform forward-time SLiM simulations to investigate how genomic distance can affect the effect that epistasis can have on the response to selection when both the interacting targets are on the same chromosome but different haplotypes.. The simulation is described in one_chr.slim. 

## Simulation Parameters

|Parameter|Description|
|------|------|
|popSize|Population size|
|numReps|Number of replicate simulations|
|genomeSize|Size of the chromosome|
|domCoefs|Dominance coefficients|
|epis|Epistasis coefficient|
|targetPos0|Vector containing the marker positions of line A and line B|
|selCoefs|Vector containing the selection coefficients of markers of line A and line B|
|focalGen|Generation to generate the output|
|seed|Seed for the simulations|


## Simulation Scenarios

-----

# General Notes

* All the simulations are parallelized using `foreach` and `doParallel`.
* Temporary binary SLiM population files are automatically removed after execution.
*  Different frequency calculations are applied for the X chromosome loci and the autosomal loci to account for the sex-specific numbers of the respective chromosomes.
* The same simulation framework can be adapted for alternative chromosome architectures or fitness models.
* The simulations assume diploid genomes and converts mutation counts to frequencies accordingly.
* Using the provided seed files is recommended for reproducible simulations (Threshold.slim).
* Output paths in the scripts should be modified before execution. 

  
------------------------
