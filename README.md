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
