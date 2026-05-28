### Loading libraries ###

library("data.table")
library("foreach")
library("doParallel")
library("poolSeq")
library("tidyr")
library("dplyr")
library("ggplot2")
library("ggpubr")
library("scales")
library("purrr")
library("boot")

### Loading functions ###
asStr <- function(x) {
  if (is.numeric(x)) {
    if (length(x) > 1) {
      paste0("'c(", paste(x, collapse = ","), ")'")
    } else {
      as.character(x)
    }
  } else {
    paste0('"\'', x, '\'"')
  }
}

slimScript <- "/Path/to/file/Threshold.slim"

# initialize the simulations 
initSlimSim <- function(popFile, popSize, initFreq, seed, slimCmd, slimScript) {
  # fix small selCoefs that cannot be read by slim
  selCoefs[abs(selCoefs) < sqrt(.Machine$double.eps)] <- 0
  slimArgs <- paste0(
    "-seed ", seed,
    " -define initPop=T",
    " -define popFile=", asStr(popFile),
    " -define useBinaryPopfile=T",
    " -define popSize=", asStr(popSize),
    " -define initFreq=", asStr(initFreq)
  )
  slimCmdLine <- paste(slimCmd, slimArgs, slimScript)
  system(slimCmdLine)
}

# run simulation with initPop = FALSE to maxGen
runSlimSim <- function(popFile, seed, maxGen, writeOutput, popSize, initFreq,
                       slimCmd, slimScript,
                       slimPostProc = paste0("| sed -En -f /Users/vidyadheeshkelkar/Desktop/SLiM_Scripts/getcounts2.sed | cut -f4,9,10")) {
  # fix small selCoefs that cannot be read by slim
  selCoefs[abs(selCoefs) < sqrt(.Machine$double.eps)] <- 0
  if (is.null(writeOutput) || writeOutput <= 0) {
    writeOutput <- maxGen  # Default to maxGen if invalid
  }
  slimArgs <- paste0(
    "-seed ", seed,
    " -define initPop=F",
    " -define popFile=", asStr(popFile),
    " -define maxGen=", maxGen,
    " -define writeOutput=", asStr(writeOutput),
    " -define popSize=", asStr(popSize),
    " -define initFreq=", asStr(initFreq)
  )
  slimCmdLine <- paste(slimCmd, slimArgs, slimScript, slimPostProc)
  system(slimCmdLine, intern = TRUE)
}
# get allele frequencies from simulations
newGetSimFreqs2 <- function(slimCmd, slimScript) {
  function(atGen, popSize, initFreq, epis, domCoefs, seedBase = seed, numReps,
           saveSims = NULL, atGens = NULL) {
    popFile <- tempfile(pattern = "pop", fileext = ".bin")
    simFreqDT <- tryCatch({
      initSlimSim(
        popFile = popFile,
        popSize = popSize,
        initFreq = initFreq,
        seed = seedBase,
        slimCmd = slimCmd,
        slimScript = slimScript
      )
      message("SLiM initSim completed, starting runSim")
      foreach(
        repId = seq_len(numReps),
        .combine = rbind,
        .export = c("fread", "runSlimSim", "asStr", "seedBase", "maxGen", "atGens", "slimCmd", "popFile", "slimScript")
      ) %dopar% {
        fread(
          text = runSlimSim(
            popFile = popFile,
            seed = seedBase + repId,
            maxGen = atGen,
            writeOutput = (if (is.null(atGens)) atGen else atGens),
            popSize = popSize,
            initFreq = initFreq,
            slimCmd = slimCmd,
            slimScript = slimScript
          )
          , header = FALSE
          , col.names = c("pos", "count", "chrom"),
          key = c("pos")
        )[
          , rep := repId
        ]
      }
    },error = function(e) {
      message("Error in foreach iteration ", repId, ": ", e)
      return(data.table(pos = integer(), count = integer(), gen = integer(), state = character()))
    })
    if (file.exists(popFile)) {
      file.remove(popFile)
    } else {
      message("File ", popFile, " does not exist.")
    }
    if (is.null(simFreqDT) || ncol(simFreqDT) == 1) {
      return(NA)
    }
    # make 0-based positions 1-based
    simFreqDT[, pos := pos + 1L]
    simFreqDT[is.na(count), count := 0L]
    simFreqDT[
      # allele frequencies
      chrom == "X" , freq :=  (2/3) * count / popSize # 2 genomes / female & 1 genome / male!
    ]
    simFreqDT[
      # allele frequencies
      chrom != "X" , freq :=  (1/2) * count / popSize # 2 genomes / individual!
    ][
      , count := NULL
    ]
    setkey(simFreqDT, pos, rep)
    if (!is.null(saveSims)) {
      saveRDS(simFreqDT, file = saveSims)
    }
    return(simFreqDT)
  }
}
                          
getSimFreqs2 <- newGetSimFreqs2(slimCmd, slimScript)

# Run simulations remotely in SLiM
run_one_experiment <- function(seed_iter) {
  # popSize -> Effective population size calculated for each population
  OS_1 <- getSimFreqs2(atGen = 2L, popSize = 243, initFreq = 0.67,
                      seedBase = seed_iter, numReps = 5L)
  
  OS_10 <- getSimFreqs2(atGen = 11L, popSize = 243, initFreq = 0.67,
                     seedBase = seed_iter, numReps = 5L)
  
  SO_1 <- getSimFreqs2(atGen = 2L, popSize = 356, initFreq = 0.33,
                       seedBase = seed_iter+2000, numReps = 5L)
  
  SO_10 <- getSimFreqs2(atGen = 11L, popSize = 356, initFreq = 0.33,
                        seedBase = seed_iter+2000, numReps = 5L)
  
  bind_freqs <- function(F1, F10, cross) {
    left_join(
      F1 %>% mutate(gen = 1),
      F10 %>% mutate(gen = 10),
      by = c("pos", "rep", "chrom"),
      suffix = c("_1", "_10")
    ) %>%
      mutate(
        s = (logit(freq_10) - logit(freq_1)) / 10,
        Cross = cross
      )
  }
  
  results <- bind_rows(
    bind_freqs(OS_1, OS_10, "OS"),
    bind_freqs(SO_1, SO_10, "SO")
  )
  
  results %>%
    group_by(chrom, pos) %>%
    summarise(
      p_value = t.test(s ~ Cross)$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"))
}

### MAIN Script ###
set.seed(123)
seeds <- sample(1e8:9e8, 100) #(Use the file threshold_seeds.txt for consistent results)

all_results <- map(seeds, run_one_experiment)

# Plotting the results
all_results_df <- bind_rows(all_results, .id = "sim")
all_results_df$sim <- as.integer(all_results_df$sim)
all_results_df$chrom <- factor(all_results_df$chrom, levels = c("X", "2L", "2R", "3L", "3R"))

        
p1 <- ggplot(all_results_df) + 
  geom_point(aes(pos, -log10(p_value), group = sim, colour = chrom)) + 
  facet_grid(~chrom, scales = "free_x") +
  labs(x = "Position (Mb)", y = "-log10(p)", color = "Chromosome") +
  theme_bw()+
  theme(
    axis.text.x = element_text(vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1, "cm"),
    legend.position = "right",
    strip.text.x = element_text(size = 10)
  ) +
  scale_x_continuous(labels = label_number(scale = 1/1e6),
                     sec.axis = sec_axis(~ . , name = "Chromosome", breaks = NULL, labels = NULL))+
  scale_y_continuous(sec.axis = sec_axis(~ . , breaks = NULL, labels = NULL))


# Calculating the significance cutoff
cutoff_all <- quantile(all_results_df$p_value, 0.05, na.rm=TRUE)
cutoff_all

cutoffs_chr <- all_results_df %>%
  group_by(chrom) %>%
  summarise(cutoff = quantile(p_value, 0.05, na.rm=TRUE))

# Distribution of p-values & significance cutoff
p2 <- ggplot(all_results_df, aes(p_value)) +
  geom_histogram(bins=100,  fill = "grey70",color = "black", linewidth = 0.2) +
  geom_vline(xintercept=cutoff_all,color="red", linetype="dashed") +
  annotate("text", x = cutoff_all, y= Inf, label = paste0("Cutoff = ", signif(cutoff_all, 3)),
           color = "red", vjust = 2,hjust = -0.1) +
  theme_bw() + 
  labs(x="p-value", y="count")
                          
