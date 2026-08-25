# Scripts for SLiM Simulations and Allele Frequency Analysis in R

This repository contains R and SLiM scripts used to analyze allele frequency dynamics and simulate selection responses in experimental evolution populations.

The workflows include:

* Allele frequency estimation and visualization.
* Selection response estimation.
* Genome-wide statistical testing.
* Forward-time SLiM simulations.
* Empirical threshold estimation for genome-wide t-tests.
* Testing for deviations from additive assumption.

---

# Repository Structure

| File                        | Description                                                              |
| --------------------------- | ------------------------------------------------------------------------ |
| `AFnS_plots.R`              | Allele frequency and selection response analysis for directional crosses |
| `H0_Check.R`                | Single-chromosome SLiM simulations for directional crosses               |
| `H0_Check2Chrom.R`          | Two-chromosome SLiM simulations for directional crosses                  |
| `threshold_for_ttest.R`     | Empirical significance threshold estimation                              |
| `process_mark.R` | Generate marker SNP catalogue for each inbred line |
| `cont_nt_pos.R`| Makes the windows continuous |
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
| `filtering_marker_snps.sh` | Filters the input data for marker positions |
| `merging_sample_markers.sh` |Generates merged sample and inbred line VCF files from cram files|
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
vcfR
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
  * `bcftools`
  * `awk`
    
* Softwares mentioned in [SNP calling pipline](https://zenodo.org/records/21506978).

---
# 1. SNP calling and window allele frequency estimation

## Scripts

* `merging_sample_markers.sh`
* `process_mark.R`
* `filtering_marker_snps.sh`
* `haploFreq.R`
* `cont_nt_pos.R`
* `predef_win_haploFreq.R`
* flag-short.awk
* info2fmt.awk
* post-merging.awk
* remapq.awk

*__NOTE: flag-short.awk, info2fmt.awk, post-merging.awk and remapq.awk are available at__* [SNP calling pipline](https://zenodo.org/records/21506978).

## Overview

These scripts process mapped pooled sequencing results to compute window allele frequency for marker SNPs. Marker SNPs were defined as SNPs that had a private allele for one parental line (either OregonR, Samarkand, or Helsinki) while the other lines shared the same alternate allele.

These scripts:

1. Call variants using bcftools pileup for all CRAM files.
2. Merge all the input BCFs.
3. Filter the merged BCF to keep only high-confidence SNPs in non-repeat regions and delete boundary SNPs around INDELs.
4. Make a merged inbred line file.
5. Make a catalogue of marker SNPs for each inbred line.
6. Filter the merged sample and inbred line files for marker SNPs.
7. Compute window coordinates to be used for estimating window allele frequencies.
8. Compute window allele frequencies.

## Required Inputs

* Aligned sequence data for sample and inbred lines as CRAM files
* _Drosophila melanogaster_ reference genome (dmel6.03)
* norepeats.bed containing non-repeat regions

## Order of the scripts/Workflow

CRAM files --> merging_sample_markers.sh --> process_mark.R --> filtering_marker_snps.sh --> haplofreq.R --> cont_nt_pos.R --> predef_win_haploFreq.R --> Downstream analysis

## Variant calling, Merging and Filtering

The first step in the pipeline is to call SNPs and merge all the SNP data using **merging_sample_markers.sh**. The script calls __*pileup.sh*__, which does variant calling with `bcftools pileup`. __*pileup.sh*__ also needs helper scripts __*info2fmt.awk*__, __*reampq.awk*__, and __*flag-short.awk*__. 

Variant calling is done for all the CRAM files present in the input folder. The resultant BCFs are stored in the BCF directory. All the BCF files are then merged using `bcftools merge`, and the script calls __*post-merging.awk*__ to perform post-merging processing. This adds additional information to the merged samples, including FORMAT/AF (observed allele frequencies), FORMAT/XF (expected allele frequencies, and FORMAT/SAD (sum of allelic depths), along with INFO/NUMALT. post-merging.awk also converts an alternate allele (one of all alternate alleles at multiallelic sites) to the reference if the reference allele is absent from all samples. 
For more details, go to: [SNP calling pipeline](https://zenodo.org/records/21506978).

The merged file is then filtered according to the following parameters using bcftools filter:
1. TYPE=\"snp\": Ensuring that we only have SNPs.
2. -T norepeats.bed: Removing repeat regions.
3. -g 5: Removing SNPs within 5 bp of INDELs.
4. INFO/DP > ((avg_depth / 2)) & INFO/DP < ((avg_depth * 2)): Removing positions with coverage depth more than 2 * average_depth or less than 1/2 * average_depth.

The filtered VCF file is then normalized.

The normalized file is then used to generate marker_flt.vcf.gz, which has only the inbred lines.

This script generates two outputs: 
1. mpileup_flt.vcf.gz: Merged sample and inbred lines VCF file, filtered and normalized. 
2. marker_flt.vcf.gz: Merged inbred lines VCF file, filtered and normalized.

## Identifying marker SNPs

Marker SNPs are identified using **process_mark.R**. It needs _marker.vcf.gz_ as the input and will produce a TXT file with the position and allelic configuration of the marker SNPs. One marker SNP catalogue for each line is generated. Marker SNPs are defined as SNPs that are fixed for the reference {0} or alternate {1} in one line and are fixed for the other allele in the other two lines. 

The three outputs are: 
1. ore_marker_file.txt
2. sam_marker_file.txt
3. hel_marker_file.txt

## Filtering for marker SNPs and samples

The outputs from merging_sample_markers.sh are filtered for marker SNPs, and sample-specific VCF files are generated using **filtering_marker_snps.sh**. The VCF files are annotated with an additional field that acts as a key to match the records from our marker catalogues. This is done using `bcftools annotate`. Filtering is done by matching the position and the SNP identity using `bcftools view`. The VCF files are filtered for markers of each inbred line, thus generating 6 VCF files.

The filtered VCF files are then filtered again for samples and corresponding inbred lines to make sample-specific merged marker and pileup files.  The resultant outputs are then used as inputs for calculating window allele frequencies.

The outputs are:
**_6 sample files_**
1. hel_ho_mpileup.vcf.gz: Merged sample and inbred line BCFs for sample HO-OH having marker SNPs of Helsinki.
2. hel_hs_mpileup.vcf.gz: Merged sample and inbred line BCFs for sample HS-SH having markers SNPs of Helsinki.
3. sam_hs_mpileup.vcf.gz: Merged sample and inbred line BCFs for sample HS-SH having markers SNPs of Samarkand.
4. sam_so_mpileup.vcf.gz: Merged sample and inbred line BCFs for sample SO-OS having markers SNPs of Samarkand.
5. ore_so_mpileup.vcf.gz: Merged sample and inbred line BCFs for sample SO-OS having markers SNPs of OregonR.
6. ore_ho_mpileup.vcf.gz: Merged sample and inbred line BCFs for sample HO-OH having markers SNPs of OregonR.
   
**_6 marker files_**
1. hel_ho_markers.vcf.gz: Merged inbred line BCFs for sample HO-OH having marker SNPs of Helsinki.
2. hel_hs_markers.vcf.gz: Merged inbred line BCFs for sample HS-SH having markers SNPs of Helsinki.
3. sam_hs_markers.vcf.gz: Merged inbred line BCFs for sample HS-SH having markers SNPs of Samarkand.
4. sam_so_markers.vcf.gz: Merged inbred line BCFs for sample SO-OS having markers SNPs of Samarkand.
5. ore_so_markers.vcf.gz: Merged inbred line BCFs for sample SO-OS having markers SNPs of OregonR.
6. ore_ho_markers.vcf.gz: Merged inbred line BCFs for sample HO-OH having markers SNPs of OregonR.

## Computing coordinates for windows

OregonR had the least number of markers; to ensure that it had enough SNPs in each window, we computed window ranges using its SNP data. To compute the windows in which approximately 500 OregonR marker SNPs were present, **haploFreq.R** was used: `haploFreq.R -i ore_so_mpileup.vcf.gz -m ore_so_markers.vcf.gz -o out_500.rds`.

The windows were made continuous using **cont_nt_pos.R**, and the output window_details.rds was used for subsequent analysis.

## Window allele frequency estimation

This is done using **predef_win_haploFreq.R**. The syntax is as follows: `predef_win_haploFreq.R -i <samples.vcf> -m <markers.vcf> -o <out.rds> -p <pre_win.rds>`

This generates 6 outputs. One for each sample:

1. hel_ho_out_500.rds
2. ore_ho_out_500.rds
3. hel_hs_out_500.rds
4. sam_hs_out_500.rds
5. ore_so_out_500.rds
6. sam_so_out_500.rds

## Directional Crosses

The variant calling and filtering for the SO & OS samples is carried out similarly to that of the _3x2 crosses (as mentioned above). Once the pileup VCF files for the samples are generated, a marker SNP VCF is created. For this, the marker SNP files `ore_marker_file.txt` and `sam_marker_file.txt` are combined and used as a regions file. The window-wise allele frequency estimates are then calculated for non-overlapping 500 SNP windows (`-w 500 -s 500`) using the R script `haplofreq.R`

`haploFreq.R -i <samples.vcf> -m <markers.vcf> -o <out.rds> -w 500 -s 500 -d 1,2`

-----

# 2. Testing for Deviations from Additive Assumption: 3x2 Crosses

## Scripts

* `main_script_ee_data.R`
* `func_ee_analysis.R`
* `func_visual_ee_data.R`

## Overview

These scripts process pooled sequencing allele frequency data from 3x2 crosses to:

* Clean and restructure sample metadata.
* Average allele-frequency estimates obtained from two independent marker sets.
* Calculate logit-transformed allele frequencies as a measure of selection response.
* Estimate deviations from additive expectations.
* Perform genome-wide statistical test for deviations from additivity.
* Correct for multiple testing using the Benjamini-Hochberg procedure.
* Classify genomic regions according to number of significant deviations.
* Generate publication-quality figures.

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

* Removes reference samples.
* Parses sample, generation, and replicate information.
* Standardizes sample and allele-frequency column names.
* Calculates the midpoint `mpos` of each genomic window.

For each pairwise comparison, allele frequencies estimated from two independent marker sets are averaged.

## Selection Response Calculation

Selection response is calculated with `logit_freq()` as the log ratio of allele frequencies: 

z = log[AF_A/(1-AF_A)]

where AF_A is the window allele frequency for line A

## Additive Expectation and Deviation and Statistical Testing

`compute_p()`:

* Tests three additive expectations.
* Correct for multiple testing.
* Classify windows as significantly deviating (1) or not significantly deviating (0) from additive expectations.
* Calculates deviation from additive expectations.

Genome-wide statistical tests are performed independently for each genomic window using linear models: 

lm( s ~ group, data = .)

`emmeans` is used to estimate the expected contrast, which is compared with the observed z of the focal population.

FDR is estimated using the ***Benjamini-Hochberg*** method. Genomic windows with ***Adjusted P ≤ 0.05*** are classified as a significant deviation from additivity. 

The three hypotheses tested are: 

H(HO-OH) - S(SO-OS) = H(HS-SH) 

O(SO-OS) - H(HS-SH) = H(HO-OH) 

S(HS-SH) - O(HO-OH) = S(SO-OS)

The deviation from the additive expectation is calculated as: ***Deviation = Observed response - Expected response***
 

## Genomic Classification

Significance across the three additive tests is combined using `find_pattern()` to classify each window: 

|Classification|Description|
|------|-------|
|0/3|No significant deviations comparisons|
|1/3|Significant deviations in one comparison|
|2/3|Significant deviations in two comparisons|
|3/3|Significant deviations in all three comparisons|

## Output Plots

### `af_comp`

Combined allele-frequency trajectories across genomic positions and chromosomes for the three pairwise comparisons.

### `plot_so_exp_obs`

Expected vs Observed `z` in SO-OS.

### `plot_class_res_10`

Genome-wide visualization of deviations from additive expectations, showing calculated versus expected logit allele-frequency differences and genomic regions classified according to the number of significant tests.

## Output Data

### `merge_df_10`

Combined significance classifications across all three comparisons.

### `merge_dev_10`

Combined genome-wide deviation estimates across all three comparisons.

-----

# 3. Epistasis Simulations: 3x2 Crosses

## Scripts

* `main_script_slim_one_chr.R`
* `main_script_slim_two_chr.R`
* `func_slim_analysis.R`
* `func_visual_slim.R`
* `one_chr.slim`
* `two_chr.slim`

## Overview

These scripts perform forward-time SLiM simulations to investigate how genomic distance affects the contribution of epistasis to the response to selection. It tests this in two cases:
1. When two interacting targets are located on the same chromosome but on different haplotypes. The simulation model is described in ***one_chr.slim***. The analysis is done using ***main_script_slim_one_chr.R***.
2. When two interacting targets are located on two different chromosomes. The simulation model is described in ***two_chr.slim***. The analysis is done using ***main_script_slim_two_chr.R***.

The pipeline tests the relationship:

**z_A(AB) = z_A(AC) − z_B(BC)**

where *z_A* and *z_B* represent the selection response of the marker allele from lines A and B, respectively. *AB*, *AC*, and *BC* represent crosses between lines A and B, lines A and C, and lines B and C, respectively. 

To investigate the effect of genomic distance, the genomic position of target A is kept constant, while the position of target B is systematically varied across 400 different genomic positions. For each position of target B, the simulation is run, and the relationship above is tested to assess how the genomic distance between the two interacting targets affects the ability to detect the epistatic contribution to the response to selection.


## Simulation Parameters

|Parameter|Description|
|------|------|
|`popSize` |Population size|
|`numReps` |Number of replicate simulations|
|`genomeSize` |Size of the chromosome|
|`domCoefs` |Dominance coefficients|
|`epis` |Epistasis coefficient|
|`targetPos0` |Vector containing the marker positions of line A and line B|
|`selCoefs` |Vector containing the selection coefficients of markers of line A and line B|
|`focalGen` |Generation to generate the output|
|`seed` |Seed for the simulations|


## Simulation Scenarios

The pipeline supports two scenarios: 
1. Scenario 1: Additive (epis = 1.00) 
2. Scenario 2: Epistasis (epis = 2.00)

## Functions for running simulations and analysis of simulation results

1. `initSlimSim()`: Function specifying the parameters for initializing the simulation and generating a temporary population file. The population file makes the subsequent simulations faster.
2. `runSlimSim()`: Function specifying the parameters for running the simulation.
3. `clean_data()` and `get_freq()`: Cleans the output from SLiM to get marker allele frequencies
4. `newGetSimFreqs2()`: Runs the simulation the specified number of times (`numReps`) and processes the output to generate marker allele frequencies for only the amrkers present in the focal population. Returns a dataframe containing marker allele frequencies for all the populations.
5. `logit_freq()` and `logit_sample()`: Calculates the selection response as logit-transformed allele frequencies.
6. `modify_result()`: Classifies deviations as significant(1) or insignificant (0).
 
## Output Plots

Both ***main_script_slim_one_chr.R*** and ***main_script_slim_two_chr.R*** generate the following plots:

### `p_dist_calc_add`

The effect of the distance between two selected targets (X-axis) and the statistical significance of deviations between predicted and calculated logit-transformed allele frequencies (Y-axis) for additive architecture.

### `p_dist_calc_int`

The effect of the distance between two selected targets (X-axis) and the statistical significance of deviations between predicted and calculated logit-transformed allele frequencies (Y-axis) for epistatic architecture.

### `s_dist_calc`

The effect of the distance between the selection targets (or markers; X-axis) on the deviation between the expected and observed response to selection (Y-axis).

### `roc_curve`

The receiver operating characteristic curve.

## Output Data

### `pred_vs_obs_calc`

The output containing the following information:
* `run`: Each run describes a unique position of target 2
* `arch`: The underlying simulated genetic architecture
* `pos`: Position of target 1
* `pred`: Predicted selection response
* `obs`: Observed selection response
* `diff_s`: pred-obs
* `diff_pos`: Distance between target 1 and target 2
* `adj_P`: FDR corrected p-values
* `p_value`: Raw p-values from hypothesis testing
* `sign`: Predicted architecture; 0- Additive and 1- Non-additive (epistatic)

-----
# 4. Dominance Simulations: 3x2 Crosses

## Scripts

* `main_script_slim_dom.R`
* `func_slim_analysis.R`
* `func_visual_slim.R`
* `one_chr.slim`

## Overview

These scripts perform forward-time SLiM simulations to investigate how genomic distance affects the contribution of dominance to the response to selection. It tests this when there are two targets located on the same chromosome but on different haplotypes. The simulation model is described in ***one_chr.slim***. The analysis is done using ***main_script_slim_dom.R***

The pipeline tests the relationship:

**z_A(AB) = z_A(AC) − z_B(BC)**

where *z_A* and *z_B* represent the selection response of the marker allele from lines A and B, respectively. *AB*, *AC*, and *BC* represent crosses between lines A and B, lines A and C, and lines B and C, respectively. To investigate the effect of genomic distance, the genomic position of target A is kept constant, while the position of target B is systematically varied across 400 different genomic positions. 


## Simulation Parameters

Described in **7. Epistasis Simulations: 3x2 crosses**

## Simulation Scenarios

The pipeline supports 9 scenarios: 

|Serial Number|Architecture|Dominance Coefficient Target 1|Dominance Coefficient Target 1|
|-----|-----|-----|----|
|1|rec-rec|0|0|
|2|rec-codom|0|0.5|
|3|codom-rec|0.5|0|
|4|codom-codom|0.5|0.5|
|5|codom-dom|0.5|1|
|6|dom-codom|1|0.5|
|7|dom-dom|1|1|
|8|dom-rec|1|0|
|9|rec-dom|0|1|

## Functions for running simulations and analysis of simulation results

Described in **6. Epistasis Simulations: 3x2 crosses**
 
## Output Data

### `pred_vs_obs_calc`

The output containing the following information:
* `run`: Each run describes a unique position of target 2
* `arch`: The underlying simulated genetic architecture
* `pos`: Position of target 1
* `pred`: Predicted selection response
* `obs`: Observed selection response
* `diff_s`: pred-obs
* `diff_pos`: Distance between target 1 and target 2
* `adj_P`: FDR corrected p-values
* `p_value`: Raw p-values from hypothesis testing
* `sign`: Predicted architecture; 0- Additive and 1- Non-additive (Dominant)
-----

# 5. Allele Frequency and Selection Strength Analysis for Directional Crosses

## Script

`AFnS_plots.R`

## Overview

This script processes pooled sequencing allele frequency data from the directional crosses to:

* Clean and restructure sample metadata.
* Calculate logit-transformed allele frequencies.
* Estimate selection response between generations (F1 → F10).
* Perform genome-wide statistical comparisons between populations.
* Generate publication-quality figures.

## Required Input

```r
df <- readRDS("/Path/to/file/AlleleFreq.rds")
```
Use the `DirCrossAF.rds` file in the `Supporting_Files` folder to regenerate the results in manuscript. You can also create your own _.rds_ file from the raw sequencing data following the SNP calling pipeline in Section 1, followed by the R script, `haploFreq.R` which generates allele frequency data from filtered _.bcf_ files.

Required data fields include:

* `chr` Chromosome
* `pos` Positions (Mean position of the window)
* `sample` The Cross direction - SO/OS
* `Dmel_OregonR_non-inbred` Oregon allele frequency

## Selection Strength Calculation

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

# 6. Single-Chromosome SLiM Simulations: Directional Crosses

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

### `summary_s`
The compilation of the selection response summaries for each scenario.

### `combined_plot`
The combined plot depicting the distance-dependent responses for each scenario.
plots need to be generated for each scenario separately before combining.

---

# 7. Two-Chromosome SLiM Simulations: Directional Crosses

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
  Alternatively, comment out the tick 2: late{} section in the SLiM script to run the additive scenario.
  
* Epistatic fitness effects (`epis = 2.00`)


## Output

### `all_results`
The tab-delimited file containing the simulation results

### `p2a` & `p2e`
Boxplots of selection strength (`s`) for additive vs. epistatic scenario comparisons

---

# 8. Empirical Significance Threshold Estimation

## Scripts

* `threshold_for_ttest.R`
* `Threshold.slim`

## Overview

This workflow estimates empirical significance thresholds for genome-wide t-tests using null simulations.

The simulations generate null distributions of allele frequency change across multiple chromosomes to evaluate the expected false-positive distribution of the t-tests. This proportion of false positives is then used as the significance cutoff for the t-tests on empirical data.

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

### `p1`
genome-wide p-value distributions

### `cutoff_all`
significance threshold estimate for combined p-values genome wide

### `cutoff_chrom`
Chromosome wise significance threshold estimates
  
### p2
Histogram of all the p-values with the cutoff value `cutoff_all`

-----


# General Notes

* All the simulations are parallelized using `foreach` and `doParallel`.
* Temporary binary SLiM population files are automatically removed after execution.
*  Different frequency calculations are applied for the X chromosome loci and the autosomal loci to account for the sex-specific numbers of the respective chromosomes.
* The same simulation framework can be adapted for alternative chromosome architectures or fitness models.
* The simulations assume diploid genomes and converts mutation counts to frequencies accordingly.
* Using the provided seed files is recommended for reproducible simulations (Threshold.slim).
* Input and output paths in the scripts should be modified before execution.

  
------------------------
