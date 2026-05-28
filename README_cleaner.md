# Selection Response Simulations and Allele Frequency Analysis

This repository contains R and SLiM scripts used to analyze allele frequency dynamics and simulate selection responses in experimental evolution populations.

The workflows include:

* Allele frequency estimation and visualization
* Selection response estimation
* Genome-wide statistical testing
* Forward-time SLiM simulations
* Empirical threshold estimation for genome-wide t-tests

---

# Repository Structure

| File                    | Description                                      |
| ----------------------- | ------------------------------------------------ |
| `AFnS_plots.R`          | Allele frequency and selection response analysis |
| `H0_Check.R`            | Single-chromosome SLiM simulations               |
| `H0_Check2Chrom.R`      | Two-chromosome SLiM simulations                  |
| `threshold_for_ttest.R` | Empirical significance threshold estimation      |
| `H0_Check.slim`         | SLiM model for single-chromosome simulations     |
| `H0_2Chrom.slim`        | SLiM model for two-chromosome simulations        |
| `Threshold.slim`        | SLiM model for null simulations                  |
| `threshold_seeds.txt`   | Seeds for reproducible threshold simulations     |
| `*.txt` marker files    | Marker/window positions used in simulations      |

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
