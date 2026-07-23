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
| `func_exp_analysis.R`       | Functions for analysis of experimental data for 3x2 crosses              |
| `func_slim_analysis.R`      | Functions for analysis of SLiM simulations for 3x2 crosses               |
| `func_visual_exp_data.R`    | Functions for visualization of experimental data for 3x2 crosses         |
| `func_visual_slim.R`        | Functions for visualization of SLiM simulation data for 3x2 crosses      |
| `main_script_exp_data.R`    | Main script for analysis of experimental data for 3x2 crosses            |
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
tidyverse
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
```

## External Software

* SLiM (v5.01 or compatible)
* GNU command-line utilities:

  * `sed`
  * `cut`

---

# 1. Allele Frequency and Selection Strength Analysis

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

Required columns:

* `chr`
* `pos`
* `sample`
* `Dmel_OregonR_non-inbred`

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

---

# 2. Single-Chromosome SLiM Simulations

## Scripts

* `H0_Check.R`
* `H0_Check.slim`

## Overview

These scripts perform forward-time SLiM simulations to evaluate how:

* linkage distance,
* haplotype structure,
* epistasis,
* and initial allele frequencies

influence selection response over time.

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

| Scenario | Epistasis                        | Haplotype Structure  |
| -------- | -------------------------------- | -------------------- |
| I        | Epistatic (`epis = 2.00`)        | Same haplotype       |
| II       | Neutral/Additive (`epis = 1.00`) | Same haplotype       |
| III      | Epistatic (`epis = 2.00`)        | Different haplotypes |
| IV       | Neutral/Additive (`epis = 1.00`) | Different haplotypes |

### Required SLiM Modifications

Scenario setup requires enabling/disabling specific sections in:

```bash
H0_Check.slim
```

* Lines `111–116` → different haplotypes
* Lines `120–125` → same haplotype
* Lines `148–163` → epistatic interaction block

## Output

The simulations generate:

* Selection response summaries
* Distance-dependent response curves
* Error-bar plots across replicate simulations

Simulation outputs are written as tab-delimited files.

---

# 3. Two-Chromosome SLiM Simulations

## Scripts

* `H0_Check2Chrom.R`
* `H0_2Chrom.slim`

## Overview

This workflow extends the simulations to selected loci located on separate chromosomes, enabling analysis of:

* interchromosomal interactions,
* additive vs. epistatic fitness effects,
* chromosome-specific selection dynamics.

## Key Parameters

| Parameter            | Description            |
| -------------------- | ---------------------- |
| `target1`, `target2` | Selected loci          |
| `epis`               | Epistasis coefficient  |
| `domCoefs`           | Dominance coefficients |

## Simulation Modes

### Additive

```r
epis = 1.00
```

Alternatively, disable the `tick 2: late{}` interaction block in the SLiM script.

### Epistatic

```r
epis = 2.00
```

## Output

The script generates:

* Selection response boxplots
* Population comparisons
* Statistical significance tests between crosses

---

# 4. Empirical Significance Threshold Estimation

## Scripts

* `threshold_for_ttest.R`
* `Threshold.slim`

## Overview

This workflow estimates empirical significance thresholds for genome-wide t-tests using null simulations.

The simulations:

* generate null allele frequency trajectories,
* calculate selection responses,
* perform genome-wide t-tests,
* estimate expected false-positive distributions.

The same marker/window positions used in allele frequency estimation are used here to preserve the multiple-testing structure of the dataset.

## Population Setup

| Population | Effective Population Size | Initial Frequency |
| ---------- | ------------------------- | ----------------- |
| `OS`       | 243                       | 0.67              |
| `SO`       | 356                       | 0.33              |

## Threshold Estimation

Selection response is calculated as:

```r
s = (logit(freq_10) - logit(freq_1)) / 10
```

Genome-wide t-tests are then performed:

```r
t.test(s ~ Cross)
```

The empirical significance cutoff is estimated as:

```r
quantile(all_results_df$p_value, 0.05)
```

Chromosome-specific cutoffs are also calculated separately.

## Output

The script generates:

* genome-wide p-value distributions,
* empirical p-value histograms,
* significance threshold estimates.

---

# General Notes

* Simulations are parallelized using `foreach` and `doParallel`.
* Temporary binary SLiM population files are automatically removed after execution.
* Chromosome-specific allele frequency calculations account for X-linked versus autosomal ploidy.
* Using the provided seed files is recommended for reproducible simulations.
* Output paths in the scripts should be modified before execution.

### t-test Significance Plot (`p3`)
Visualizes genome-wide statistical significance, with significant loci highlighted.

------------------------

# SLiM Simulation Pipeline for Selection Response Analysis
## Overview
_`H0_Check.R`_ performs forward-time population genetic simulations using SLiM to investigate how linkage distance, haplotype structure, and epistasis influence allele frequency change and selection response over time.

The workflow:
* Initializes SLiM simulations with user-defined selection parameters
* Simulates allele frequency trajectories across generations
* Calculates selection response (`s`) between generations
* Compares different population crosses and target configurations
* Summarizes and visualizes simulation outcomes across genomic distances

The repository also includes the required SLiM simulation script (`H0_Check.slim`), which controls the evolutionary model and haplotype setup.

### External Software

* SLiM (v5.01 or compatible)
* GNU utilities (`sed`, `cut`) for post-processing simulation output

## Simulation Design
The simulations evaluate the effect of:
* Distance between selected loci
* Epistatic vs. neutral interactions
* Same vs. separate haplotypes
* Initial allele frequency differences between populations

### Key Parameters

| Parameter      | Description                       |
| -------------- | --------------------------------- |
| `popSize`      | Population size                   |
| `numReps`      | Number of replicate simulations   |
| `target1`      | Reference selected locus          |
| `distances_mb` | Distance between selected targets |
| `selCoefs`     | Selection coefficients            |
| `epis`         | Epistasis coefficient             |
| `domCoefs`     | Dominance coefficient             |

---

## Simulation Scenarios & Required modifications

The script is designed to be run under four evolutionary scenarios:

1. ### Targets on the same haplotype — epistatic
   set `epis = 2.00`
   Comment out the chunk from line `111 to 116` and use the chunk from line `120 to 125` of the SLiM script
2. ### Targets on the same haplotype — neutral
   set `epis = 1.00` OR Comment out the `tick 2: late{}` chunk (lines `148 to 163`)
   Comment out the chunk from line `111 to 116` and use the chunk from line `120 to 125` of the SLiM script
3. ### Targets on different haplotypes — epistatic
   set `epis = 2.00`
   Comment out the chunk from line `120 to 125` and use the chunk from line `111 to 116` of the SLiM script
4. ### Targets on different haplotypes — neutral
   set `epis = 1.00` OR Comment out the `tick 2: late{}` chunk (lines `148 to 163`)
   Comment out the chunk from line `120 to 125` and use the chunk from line `111 to 116` of the SLiM script

These configurations are controlled within the included SLiM script:
```bash
H0_Check.slim
```
---

## Workflow Summary
`initSlimSim()` creates an initial population state and writes a temporary binary population file.
`runSlimSim()` executes forward simulations to specified generations using SLiM.
`getSimFreqs2()` parses simulation output and calculates allele frequencies across replicates.

### 4. Calculate Selection Response
Selection response is estimated as:
```r
s = (logit(freq_10) - logit(freq_1)) / 10
```
### 5. Summarize Across Distances
Mean selection response and standard errors are calculated for each:
* Target locus
* Population cross
* Genomic distance
---

## Output
### Simulation Results
Results are exported as tab-delimited text files:
```r
write.table(all_results, "filename.txt")
```
### Plots
The script generates:
* Selection response vs. genomic distance
* Error-bar summaries across replicates
* Multi-panel combined figures for all simulation scenarios
---

## Notes
* Temporary binary population files are automatically cleaned after execution.
* Distances are evaluated from near-complete linkage to 20 Mb separation.
* The script assumes diploid genomes and converts mutation counts to frequencies accordingly.
* The script needs to be run 4 times, once for each scenario by making necessary changes to the SLiM and/or R-script

----------------------------

# Two-Chromosome SLiM Simulation Pipeline

## Overview
_`H0_Check2Chrom.R`_ script performs forward-time population genetic simulations using SLiM to investigate selection response under a two-chromosome model.Unlike the single-chromosome simulations, this framework simulates selected loci located on separate chromosomes, allowing analysis of:

* Interchromosomal interactions
* Epistatic vs. additive fitness effects
* Cross-specific selection responses
* Allele frequency dynamics under differing starting frequencies

The required SLiM script (`H0_2Chrom.slim`) is included in the repository and defines the evolutionary model used during simulations.

---

### External Software

* SLiM (v5.01 or compatible)
* GNU command-line utilities (`sed`, `cut`) for simulation output parsing
---

## Simulation Design
The simulations model two selected targets positioned on separate chromosomes.

### Key Parameters

| Parameter            | Description                     |
| -------------------- | ------------------------------- |
| `popSize`            | Population size                 |
| `numReps`            | Number of replicate simulations |
| `target1`, `target2` | Selected target loci            |
| `selCoefs`           | Selection coefficients          |
| `epis`               | Epistasis coefficient           |
| `domCoefs`           | Dominance coefficients          |

---

## Simulation Scenarios

The script can be used to compare:

* Additive fitness effects (`epis = 1.00`)
  Alternatively, comment out the tick 2: late{} section in the SLiM script to run the additive scenario
  
* Epistatic fitness effects (`epis = 2.00`)


---

## Workflow Summary
`initSlimSim()` creates an initial population state using selected loci and starting allele frequencies.
`runSlimSim()` executes replicate simulations across generations using SLiM. 
`getSimFreqs2()` parses mutation counts and converts them to allele frequencies.

Different frequency calculations are applied for:
* X chromosome loci
* Autosomal loci
  
to account for chromosome-specific ploidy.

### 4. Calculate Selection Response
Selection response is estimated as:
```r
s = (logit(freq_10) - logit(freq_1)) / 10
```
where allele frequencies are compared between generations F1 and F10.

---

## Output
### Simulation Results
Simulation outputs are saved as tab-delimited files:
```r
write.table(all_results, "filename.txt")
```
### Plots
The script generates:
* Boxplots of selection response (`s`)
* Statistical comparisons between populations
* Additive vs. epistatic scenario comparisons

Significance testing is performed using t-tests via:
```r
stat_compare_means(method = "t.test")
```
---

## Notes
* Temporary SLiM population files are automatically removed after execution.
* The same simulation framework can be adapted for alternative chromosome architectures or fitness models.
* Parameters controlling additive versus epistatic interactions are defined both in the R script and the accompanying SLiM script.

-------------------------------

# Estimate empirical significance thresholds for genome-wide t-tests

## Overview

_`threshold_for_ttest.R`_ script uses forward-time SLiM simulations to estimate empirical significance thresholds for genome-wide t-tests comparing selection responses between experimental populations.

The workflow:
* Simulates allele frequency trajectories under a null evolutionary model
* Calculates selection response (`s`) between generations
* Performs genome-wide t-tests between populations
* Generates empirical distributions of p-values
* Estimates significance cutoffs for downstream statistical analyses

The repository also includes:
## Supporting Files
* The supporting SLiM script (`Threshold.slim`)
* The marker `.txt` files required for simulations and chromosome setup.
*  Seed files for reproducible simulations (`threshold_seeds.txt`)
---

### External Software

* SLiM (v5.01 or compatible)
* GNU command-line tools (`sed`, `cut`) for parsing simulation output
---

## Simulation Design
The simulations generate null distributions of allele frequency change across multiple chromosomes to evaluate the expected false-positive distribution of t-tests.
The simulation uses the same window positions/marker positions as used in the allele frequency estimation to keep the number of windows consistent & test the multiple testing load due to the high number of windows.

### Population Setup
Two experimental populations are simulated:
* `OS`
* `SO`

Each population uses:
* Different effective population sizes 
* Different initial allele frequencies

Example:
```r id="qzpryf"
OS: popSize = 243, initFreq = 0.67
SO: popSize = 356, initFreq = 0.33
```

---

## Workflow Summary
`initSlimSim()` creates initial SLiM populations and temporary binary population files.
`runSlimSim()` executes replicate simulations for specified generations.
`getSimFreqs2()` parses mutation counts and converts them to allele frequencies.

Chromosome-specific frequency calculations are applied for:
* X chromosome loci
* Autosomal loci
to account for differences in ploidy.

### 4. Calculate Selection Response
Selection response is estimated as:
```r 
s = (logit(freq_10) - logit(freq_1)) / 10
```
using allele frequencies from generations F1 and F10.

### 5. Perform Genome-wide t-tests
For each genomic position, a t-test compares selection responses between populations:
```r
t.test(s ~ Cross)
```
---

## Output

### Simulation Results
The script generates:
* Genome-wide p-values
* Empirical significance thresholds

### Plots
#### Genome-wide p-value Distribution
Scatterplots of:
```r
-log10(p_value)
```
across chromosomes.

#### Empirical p-value Histogram
Histogram of simulated p-values with estimated significance cutoff.

---

## Significance Threshold Estimation
The empirical significance threshold is calculated as the 5th percentile of simulated p-values:
```r
quantile(all_results_df$p_value, 0.05)
```
Chromosome-specific thresholds are also estimated separately.

---

## Notes
* Temporary SLiM population files are automatically removed after execution.
* The script is designed to estimate empirical null distributions for downstream genome-wide selection analyses.
* Using the provided seed file is recommended for reproducible threshold estimates.
