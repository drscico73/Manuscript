###Loading libraries ###
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

registerDoParallel(cores = parallel::detectCores())
outFileRds <- "/Path/to/file/filename.rds"
numCores <- 12L
numReps <- 100L 
slimCmd <- "/usr/local/bin/slim"
slimScript <- "/Path/to/file/toy04.slim"
popSize <- 1250L
focalGen <- 10L
seed=sample(100000000:999999999,1)

### Load Functions ###
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

# initialize the simulation 
initSlimSim <- function(targetPos1, selCoefs, popFile, popSize, initFreq, epis, domCoefs, seed, slimCmd, slimScript) {
  # fix small selCoefs that cannot be read by slim
  selCoefs[abs(selCoefs) < sqrt(.Machine$double.eps)] <- 0
  slimArgs <- paste0(
    "-seed ", seed,
    " -define initPop=T",
    " -define targetPos0=", asStr(targetPos1 - 1L),
    " -define selCoefs=", asStr(selCoefs),
    " -define popFile=", asStr(popFile),
    " -define useBinaryPopfile=T",
    " -define popSize=", asStr(popSize),
    " -define initFreq=", asStr(initFreq),
    " -define epis=", asStr(epis),
    " -define domCoefs=", asStr(domCoefs)
  )
  slimCmdLine <- paste(slimCmd, slimArgs, slimScript)
  system(slimCmdLine)
}
# run simulation with initPop = FALSE to maxGen
runSlimSim <- function(targetPos1, selCoefs, popFile, seed, maxGen, writeOutput, popSize, initFreq, epis, domCoefs,
                       slimCmd, slimScript,
                       slimPostProc = paste0("| sed -En -f /Users/vidyadheeshkelkar/Desktop/SLiM_Scripts/getcounts.sed | cut -f4,9,11,13")) {
  # fix small selCoefs that cannot be read by slim
  selCoefs[abs(selCoefs) < sqrt(.Machine$double.eps)] <- 0
  if (is.null(writeOutput) || writeOutput <= 0) {
    writeOutput <- maxGen  # Default to maxGen if invalid
  }
  slimArgs <- paste0(
    "-seed ", seed,
    " -define initPop=F",
    " -define targetPos0=", asStr(targetPos1 - 1L),
    " -define selCoefs=", asStr(selCoefs),
    " -define popFile=", asStr(popFile),
    " -define maxGen=", maxGen,
    " -define writeOutput=", asStr(writeOutput),
    " -define popSize=", asStr(popSize),
    " -define initFreq=", asStr(initFreq),
    " -define epis=", asStr(epis),
    " -define domCoefs=", asStr(domCoefs)
  )
  slimCmdLine <- paste(slimCmd, slimArgs, slimScript, slimPostProc)
  system(slimCmdLine, intern = TRUE)
}

# get allele frequency data table from simulations
newGetSimFreqs2 <- function(slimCmd, slimScript) {
  function(atPos, atGen, targetPos, selCoefs, popSize, initFreq, epis, domCoefs, seedBase = seed, numReps,
           saveSims = NULL, atGens = NULL) {
    popFile <- tempfile(pattern = "pop", fileext = ".bin")
    simFreqDT <- tryCatch({
      initSlimSim(
        targetPos1 = targetPos,
        selCoefs = selCoefs,
        popFile = popFile,
        popSize = popSize,
        initFreq = initFreq,
        epis = epis,
        domCoefs = domCoefs,
        seed = seedBase,
        slimCmd = slimCmd,
        slimScript = slimScript
      )
      message("SLiM initSim completed, starting runSim")
      foreach(
        repId = seq_len(numReps),
        .combine = rbind,
        .export = c("fread", "runSlimSim", "asStr", "targetPos", "seedBase", "maxGen", "atGens", "slimCmd", "popFile", "slimScript")
      ) %dopar% {
        fread(
          text = runSlimSim(
            targetPos1 = targetPos,
            selCoefs = selCoefs,
            popFile = popFile,
            seed = seedBase + repId,
            maxGen = atGen,
            writeOutput = (if (is.null(atGens)) atGen else atGens),
            popSize = popSize,
            initFreq = initFreq,
            epis = epis,
            domCoefs = domCoefs,
            slimCmd = slimCmd,
            slimScript = slimScript
          )
          , header = FALSE
          , col.names = c("pos", "count", "gen","state"),
          key = c("gen", "pos")
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
    # simulations may fail if optim makes a step too wild
    if (is.null(simFreqDT) || ncol(simFreqDT) == 1) {
      return(NA)
    }
    # make 0-based positions 1-based
    simFreqDT[, pos := pos + 1L]
    # fixed mutations have a "wrong" count (is tick where the mut became fixed ...)
    simFreqDT[state == "F", count := (3/2) * popSize]
    simFreqDT[, state := NULL]
    # freqDT may be missing entries because of lost mutations
    simFreqDT <- simFreqDT[, .SD[data.table(pos = atPos), on = "pos"], keyby = .(gen, rep)]
    simFreqDT[is.na(count), count := 0L]
    simFreqDT[
      # allele frequencies
      , freq :=  (2/3) * count / popSize # 2 genomes / individual!
    ][
      , count := NULL
    ]
    setkey(simFreqDT, pos, gen, rep)
    if (!is.null(saveSims)) {
      saveRDS(simFreqDT, file = saveSims)
    }
    return(simFreqDT)
  }
}
getSimFreqs2 <- newGetSimFreqs2(slimCmd, slimScript)

# setup the target parameters for simulations 
distances_mb <- c(1e-6, seq(1,20, by=0.5))
distances_bp <- distances_mb * 1e6
target1 <- 1000000L

i <- 1L
# Setting up the simulations in SLiM remotely
run_distance <- function(dist_bp, seed) {
  
  allTargetPos <- c(target1, target1 + dist_bp)
  allSelCoefs <- c(0.05, 0.05)
  
  sortInfo <- sort(allTargetPos, index.return = TRUE)
  targetPos <- sortInfo$x
  selCoefs  <- allSelCoefs[sortInfo$ix]
  epis = 2.00 #1.00 for neutral scenario
  domCoefs <- c(1, 1)
  
  ## OS
  OS_1 <- getSimFreqs2(
    atPos = targetPos, atGen = 2L,
    targetPos = targetPos, selCoefs = selCoefs,
    popSize = popSize, initFreq = 0.67,
    epis = epis, domCoefs = domCoefs,
    seedBase = seed, numReps = 1000L
  )
  
  OS_10 <- getSimFreqs2(
    atPos = targetPos, atGen = 11L,
    targetPos = targetPos, selCoefs = selCoefs,
    popSize = popSize, initFreq = 0.67,
    epis = epis, domCoefs = domCoefs,
    seedBase = seed, numReps = 1000L
  )
  
  ## SO
  SO_1 <- getSimFreqs2(
    atPos = targetPos, atGen = 2L,
    targetPos = targetPos, selCoefs = selCoefs,
    popSize = popSize, initFreq = 0.33,
    epis = epis, domCoefs = domCoefs,
    seedBase = seed, numReps = 1000L
  )
  
  SO_10 <- getSimFreqs2(
    atPos = targetPos, atGen = 11L,
    targetPos = targetPos, selCoefs = selCoefs,
    popSize = popSize, initFreq = 0.33,
    epis = epis, domCoefs = domCoefs,
    seedBase = seed, numReps = 1000L
  )

  # combining the results from all simulations
  bind_freqs <- function(F1, F10, cross) {
    left_join(
      F1 %>% mutate(gen = 1),
      F10 %>% mutate(gen = 10),
      by = c("pos", "rep"),
      suffix = c("_1", "_10")
    ) %>%
      mutate(
        s = (logit(freq_10) - logit(freq_1)) / 10,
        Cross = cross,
        distance_mb = dist_bp / 1e6
      )
  }
  
  bind_rows(
    bind_freqs(OS_1, OS_10, "OS"),
    bind_freqs(SO_1, SO_10, "SO")
  )
}

### MAIN Script ###
set.seed(123)
seed <- sample(1e8:9e8, length(distances_bp))

# run simulations - 
# need to be run 4 times, changing the parameters to suit each case
# targets on same haplotype - epistatic & neutral
# targets on different haplotypes - epistatic & neutral
# changes need to be made in the SLiM Script

all_results <- map2_dfr(
  distances_bp,
  seed,
  run_distance
)
# save the results
write.table(all_results, "/Path/to/file/filename.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  
all_results <- all_results %>%
    mutate(
      Target = factor(
        if_else(pos == target1, 1L, 2L),
        levels = c(1, 2),
        labels = c("Target 1", "Target 2")
        )
      )

summary_s <- all_results %>%
  group_by(distance_mb, Target, Cross) %>%
  summarise(
    mean_s = mean(s, na.rm = TRUE),
    se_s   = sd(s, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
    )
                          
# plotting the results
ggplot(subset(summary_s, Target=="Target 1"), aes(x = distance_mb, y = mean_s, color = Cross, 
                      group = interaction(Target, Cross))) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = mean_s - se_s, ymax = mean_s + se_s), width = 0.03) +
  labs(x = "Distance between targets (Mb)", y = "Selection response (s)", color = "Cross", linetype = "Target") +
  scale_color_manual(values = c("darkgreen", "orange1")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text  = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12)
  )
  

combined_plot <- ((p1|p2) / (p3|p4)) +
  plot_layout(heights = c(2, 2), guides = "collect") &
  # labs(color = "New legend title") +
  plot_annotation(tag_levels = list(c("C", "D", "E", "F"))) &
  theme_bw() &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    axis.text.x = element_text(vjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "cm"),
    legend.position = "bottom",
    strip.text.x = element_text(size = 12)
  ) 
