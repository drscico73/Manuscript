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

registerDoParallel(cores = parallel::detectCores())
outFileRds <- "/Path/to/file/filename.rds"
numCores <- 12L
numReps <- 100L 
slimCmd <- "/usr/local/bin/slim"
slimScript <- "/Path/to/file/H0_2Chrom.slim"
popSize <- 1250L
focalGen <- 10L

### Loading functions ###
# initialize sumulations
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
                       slimPostProc = paste0("| sed -En -f /Users/vidyadheeshkelkar/Desktop/SLiM_Scripts/getcounts2.sed | cut -f4,9,10")) {
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
    # simulations may fail if optim makes a step too wild
    if (is.null(simFreqDT) || ncol(simFreqDT) == 1) {
      return(NA)
    }
    # make 0-based positions 1-based
    simFreqDT[, pos := pos + 1L]
    simFreqDT[is.na(count), count := 0L]
    simFreqDT[
      # allele frequencies
      chrom == "X" , freq :=  (2/3) * count / popSize # 2 genomes in females, only one in males!
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

distances_mb <- c(0)
distances_bp <- distances_mb * 1e6
target1 <- 1000000L
target2 <- 1000000L

run_once <- function(dist_bp, seed) {
  
  allTargetPos <- c(target1, target2 + dist_bp)
  allSelCoefs <- c(0.05, 0.05)
  
  sortInfo <- sort(allTargetPos, index.return = TRUE)
  targetPos <- sortInfo$x
  selCoefs  <- allSelCoefs[sortInfo$ix]
  epis <- 2.00 # 1.00 for additive scenario
  domCoefs <- c(0.5, 0.5)
  
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
  
  bind_freqs <- function(F1, F10, cross) {
    left_join(
      F1 %>% mutate(gen = 1),
      F10 %>% mutate(gen = 10),
      by = c("chrom", "pos", "rep"),
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
    bind_freqs(SO_1, SO_10, "SO") )
}

### MAIN Script ###
                        
set.seed(123)
seed_base <- sample(1e8:9e8,length(distances_bp))
                          
# run simulations
all_results <- map2_dfr(
  distances_bp,
  seed_base,
  run_once
)

# save results                          
write.table(all_results, "/Path/to/file/filename.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# load the result files                        
a <- read.table("/Path/to/file/filename.txt", sep = "\t", header = TRUE)  # additive scenario
e <- read.table("/Path/to/file/filename.txt", sep = "\t", header = TRUE)  # epistatic scenario

# plotting the results
# additive scenario                          
p2a <- ggplot(subset(a,chrom=="X" & pos == 1000000),
              aes(factor(Cross), s , group = Cross, fill = Cross)) + 
       geom_boxplot() +
       labs(x = "Distance between targets (Mb)", y = "Selection response (s)", color = "Cross", linetype = "Target") + 
       scale_fill_manual(values = c("darkgreen", "orange1")) +
       theme_bw() +
       theme(
         axis.title = element_text(size = 14, face = "bold"),
         axis.text  = element_text(size = 12),
         legend.title = element_text(size = 13),
         legend.text  = element_text(size = 12),
         legend.position = "bottom") +
      stat_compare_means(method = "t.test", label = "p.signif" , label.x = 1.4, label.y = 0.08, size = 6)

# epistatic scenario
p2e <- ggplot(subset(a,chrom=="X" & pos == 1000000),
              aes(factor(Cross), s , group = Cross, fill = Cross)) + 
       geom_boxplot() +
       labs(x = "Distance between targets (Mb)", y = "Selection response (s)", color = "Cross", linetype = "Target") + 
       scale_fill_manual(values = c("darkgreen", "orange1")) +
       theme_bw() +
       theme(
         axis.title = element_text(size = 14, face = "bold"),
         axis.text  = element_text(size = 12),
         legend.title = element_text(size = 13),
         legend.text  = element_text(size = 12),
         legend.position = "bottom") +
      stat_compare_means(method = "t.test", label = "p.signif" , label.x = 1.4, label.y = 0.20, size = 6)


  
