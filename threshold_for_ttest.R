############# Libraries #########

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

############## OLD threshold code ##########

registerDoParallel(cores = parallel::detectCores())
outFileRds <- "/Users/vidyadheeshkelkar/Desktop/s_h_estimation/QC/SO70_04.rds"
numCores <- 12L
numReps <- 100L 
slimCmd <- "/usr/local/bin/slim"
slimScript <- "/Users/vidyadheeshkelkar/Desktop/Warmup_Project/test/toy04.slim"
popSize <- 1250L
focalGen <- 10L
seed=sample(100000000:999999999,1)

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
newGetSimFreqs2 <- function(slimCmd, slimScript) {
  # lastAtPos <- lastTargetPos <- lastSelCoefs <- simFreqs <- NULL
  function(atPos, atGen, targetPos, selCoefs, popSize, initFreq, epis, domCoefs, seedBase = seed, numReps,
           saveSims = NULL, atGens = NULL) {
    # if (!is.null(lastAtPos) && length(atPos) == length(lastAtPos) && all(atPos == lastAtPos) && all(targetPos == lastTargetPos) && all(selCoefs == lastSelCoefs)) {
    #   return(simFreqs)
    # } else {
    #   lastAtPos <<- atPos
    #   lastTargetPos <<- targetPos
    #   lastSelCoefs <<- selCoefs
    # }
    # newGetSimFreqs <- function(atPos, atGen, targetPos, selCoefs, seedBase = seed, numReps = 100L,
    #                            slimCmd, slimScript, saveSims = NULL, atGens = NULL) {
    #   lastAtPos <- lastTargetPos <- lastSelCoefs <- simFreqs <- NULL
    #   if (!is.null(lastAtPos) && length(atPos) == length(lastAtPos) && all(atPos == lastAtPos) && all(targetPos == lastTargetPos) && all(selCoefs == lastSelCoefs)) {
    #     return(simFreqs)
    #   } else {
    #           lastAtPos <<- atPos
    #           lastTargetPos <<- targetPos
    #           lastSelCoefs <<- selCoefs
    #   }
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
    #simFreqs <<- matrixStats::colMedians(matrix(data = simFreqDT[.(atGen), freq], nrow = numReps))
    # numPos <- length(atPos)
    # numGens <- length(atGen)
    # array(data = simFreqDT[, freq], dim = c(numReps, numGens, numPos))
    return(simFreqDT)
  }
}
getSimFreqs2 <- newGetSimFreqs2(slimCmd, slimScript)

run_one_experiment <- function(seed_iter) {
  
  # ---- simulations ----
  OS_1  <- getSimFreqs2(atPos = allPos, atGen = 2L,
                        targetPos = allPos, selCoefs = selCoefs,
                        popSize = popSize, initFreq = 0.67,
                        epis = epis, domCoefs = domCoefs,
                        seedBase = seed_iter, numReps = 5L)
  
  OS_10 <- getSimFreqs2(atPos = allPos, atGen = 11L,
                        targetPos = allPos, selCoefs = selCoefs,
                        popSize = popSize, initFreq = 0.67,
                        epis = epis, domCoefs = domCoefs,
                        seedBase = seed_iter + 1000, numReps = 5L)
  
  SO_1  <- getSimFreqs2(atPos = allPos, atGen = 2L,
                        targetPos = allPos, selCoefs = selCoefs,
                        popSize = popSize, initFreq = 0.33,
                        epis = epis, domCoefs = domCoefs,
                        seedBase = seed_iter + 2000, numReps = 5L)
  
  SO_10 <- getSimFreqs2(atPos = allPos, atGen = 11L,
                        targetPos = allPos, selCoefs = selCoefs,
                        popSize = popSize, initFreq = 0.33,
                        epis = epis, domCoefs = domCoefs,
                        seedBase = seed_iter + 3000, numReps = 5L)
  
  bind_freqs <- function(F1, F10, cross) {
    left_join(
      F1 %>% mutate(gen = 1),
      F10 %>% mutate(gen = 10),
      by = c("pos", "rep"),
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
    group_by(pos) %>%
    summarise(
      p_value = t.test(s ~ Cross)$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      significant = p_adj < 0.05
    )
}



allPos <- as.integer(unique(obsDt[, pos]))
allSelCoefs <- c(runif(length(allPos), -0.05,0.05))
domCoefs = c(rep(0.5,length(allPos)))
epis=1.5

statsTable <- targetTable <- freqTable <- data.table()
popSize <- unique(t$popSize)

message("Starting iteration ", i, " at ", Sys.time())

# TOGETHER #
sortInfo <- sort(allPos, index.return = TRUE)
targetPos <- sortInfo$x
selCoefs <- allSelCoefs[sortInfo$ix]

focalPos <- c(targetPos)
atPos <- focalPos
atGen <- focalGen

set.seed(123)
seeds <- sample(1e8:9e8, 100)

all_results <- map(seeds, run_one_experiment)


library(tidyr)

combined <- bind_rows(all_results, .id = "experiment")

power_by_pos <- combined %>%
  group_by(pos) %>%
  summarise(
    times_significant = sum(significant),
    power = mean(significant),
    .groups = "drop"
  )

power_by_pos

highlight_pos <- c()
power_by_pos <- power_by_pos %>%
  mutate(highlight = ifelse(pos %in% highlight_pos,
                            "Epis", "Normal"))

ggplot(power_by_pos) +
  geom_point(aes(pos, power*100)) + 
  labs(x = "Position (Mb)", y = "% False Positives") +
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
  scale_x_continuous(labels = label_number(scale = 1/1e6), sec.axis = sec_axis(~ . , breaks = NULL, labels = NULL))+
  scale_y_continuous(sec.axis = sec_axis(~ . , breaks = NULL, labels = NULL))

sum(power_by_pos$power <= 0.05)

################## NEW Threshold Code ############

t <- read.table("/Users/vidyadheeshkelkar/Desktop/Warmup_Project/simsOS/targets5.txt", header = TRUE)

slimScript <- "/Users/vidyadheeshkelkar/Desktop/Warmup_Project/test/Threshold.slim"
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
newGetSimFreqs2 <- function(slimCmd, slimScript) {
  # lastAtPos <- lastTargetPos <- lastSelCoefs <- simFreqs <- NULL
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
    # simulations may fail if optim makes a step too wild
    if (is.null(simFreqDT) || ncol(simFreqDT) == 1) {
      return(NA)
    }
    # make 0-based positions 1-based
    simFreqDT[, pos := pos + 1L]
    # fixed mutations have a "wrong" count (is tick where the mut became fixed ...)
    # simFreqDT[state == "F", count := (3/2) * popSize]
    # simFreqDT[state == "F" & chrom != "X", count := (2) * popSize]
    # simFreqDT[, state := NULL]
    # freqDT may be missing entries because of lost mutations
    # simFreqDT <- simFreqDT[, .SD[data.table(pos = atPos), on = "pos"], keyby = .(gen, rep)]
    simFreqDT[is.na(count), count := 0L]
    simFreqDT[
      # allele frequencies
      chrom == "X" , freq :=  (2/3) * count / popSize # 2 genomes / individual!
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
    #simFreqs <<- matrixStats::colMedians(matrix(data = simFreqDT[.(atGen), freq], nrow = numReps))
    # numPos <- length(atPos)
    # numGens <- length(atGen)
    # array(data = simFreqDT[, freq], dim = c(numReps, numGens, numPos))
    return(simFreqDT)
  }
}
getSimFreqs2 <- newGetSimFreqs2(slimCmd, slimScript)


run_one_experiment <- function(seed_iter) {
  
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

set.seed(123)
seeds <- sample(1e8:9e8, 100)

all_results <- map(seeds, run_one_experiment)

all_results_df <- bind_rows(all_results, .id = "sim")
all_results_df$sim <- as.integer(all_results_df$sim)
all_results_df$chrom <- factor(all_results_df$chrom, levels = c("X", "2L", "2R", "3L", "3R"))

plot_chr <- function(chr) {
  ggplot(
    filter(all_results_df, chrom == chr),
    aes(x = pos, y = -log10(p_value))
  ) +
    geom_point(alpha = 0.5, size = 0.7) +
    labs(
      title = paste("Chromosome", chr),
      x = "Position",
      y = "-log10(p)"
    ) +
    theme_bw()
}

# plot_chr("X")
# plot_chr("2L")
# plot_chr("2R")
# plot_chr("3L")
# plot_chr("3R")

ggplot(all_results_df) + 
  geom_point(aes(pos, -log10(p_value), group = sim, colour = chrom)) + 
  # geom_hline(yintercept = -log10(0.05)) +
  facet_grid(~chrom, scales = "free_x") +
  labs(x = "Position (Mb)", y = "-log10(p)") +
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
  scale_x_continuous(labels = label_number(scale = 1/1e6), sec.axis = sec_axis(~ . , name = "Chromosome", breaks = NULL, labels = NULL))+
  scale_y_continuous(sec.axis = sec_axis(~ . , breaks = NULL, labels = NULL))


# chr_sizes <- all_results_df %>%
#   group_by(chrom) %>%
#   summarise(chr_len = max(pos))
# 
# chr_offsets <- chr_sizes %>%
#   mutate(offset = lag(cumsum(chr_len), default = 0))
# 
# manhattan_df <- all_results_df %>%
#   left_join(chr_offsets, by="chrom") %>%
#   mutate(global_pos = pos + offset)
# 
# ggplot(manhattan_df,
#        aes(global_pos, -log10(p_value), color=chrom)) +
#   geom_point(size=0.6, alpha=0.6) +
#   scale_colour_manual(values = c("red","blue","green","black","orange")) +
#   scale_x_continuous(labels = label_number(scale = 1/1e6), sec.axis = sec_axis(~ . , name = "Chromosome", breaks = NULL, labels = NULL))+
#   scale_y_continuous(sec.axis = sec_axis(~ . , breaks = NULL, labels = NULL))
#   labs(x="Genome position", y="-log10(p)")

cutoff_all <- quantile(all_results_df$p_value, 0.05, na.rm=TRUE)
cutoff_all

cutoffs_chr <- all_results_df %>%
  group_by(chrom) %>%
  summarise(cutoff = quantile(p_value, 0.05, na.rm=TRUE))

ggplot(all_results_df,
       aes(pos, -log10(p_value))) +
  geom_point(alpha=0.5) +
  facet_wrap(~chrom, scales="free_x") +
  geom_hline(
    data = cutoffs_chr,
    aes(yintercept = -log10(cutoff)),
    color="red",
    linetype="dashed"
  ) +
  theme_bw()


ggplot(all_results_df, aes(p_value)) +
  geom_histogram(bins=100,  fill = "grey70",color = "black", linewidth = 0.2) +
  geom_vline(xintercept=cutoff_all,color="red", linetype="dashed") +
  annotate("text", x = cutoff_all, y= Inf, label = paste0("Cutoff = ", signif(cutoff_all, 3)),
           color = "red", vjust = 2,hjust = -0.1) +
  theme_bw() + 
  labs(x="p-value", y="count")


########
cutoff_all_adj <- quantile(all_results_df$p_adj, 0.05, na.rm=TRUE)
cutoff_all_adj

cutoffs_chr_adj <- all_results_df %>%
  group_by(chrom) %>%
  summarise(cutoff = quantile(p_adj, 0.05, na.rm=TRUE))

ggplot(all_results_df,
       aes(pos, -log10(p_adj))) +
  geom_point(alpha=0.5) +
  facet_wrap(~chrom, scales="free_x") +
  geom_hline(
    data = cutoffs_chr_adj,
    aes(yintercept = -log10(cutoff)),
    color="red",
    linetype="dashed"
  ) +
  theme_bw()


ggplot(all_results_df, aes(p_adj)) +
  geom_histogram(bins=100,  fill = "grey70",color = "black", linewidth = 0.2) +
  geom_vline(xintercept=cutoff_all_adj,color="red", linetype="dashed") +
  annotate("text", x = cutoff_all_adj, y = Inf, label = paste0("Cutoff = ", signif(cutoff_all_adj, 3)),
           color = "red", vjust = 2, hjust = -0.1) +
  theme_bw() +
  labs(x="p-value", y="count")


