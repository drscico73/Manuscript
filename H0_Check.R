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
    #simFreqs <<- matrixStats::colMedians(matrix(data = simFreqDT[.(atGen), freq], nrow = numReps))
    # numPos <- length(atPos)
    # numGens <- length(atGen)
    # array(data = simFreqDT[, freq], dim = c(numReps, numGens, numPos))
    return(simFreqDT)
  }
}
getSimFreqs2 <- newGetSimFreqs2(slimCmd, slimScript)

distances_mb <- c(1e-6, seq(1,20, by=0.5))
distances_bp <- distances_mb * 1e6
target1 <- 1000000L

i <- 1L

run_distance <- function(dist_bp, seed) {
  
  allTargetPos <- c(target1, target1 + dist_bp)
  allSelCoefs <- c(0.05, 0.05)
  
  sortInfo <- sort(allTargetPos, index.return = TRUE)
  targetPos <- sortInfo$x
  selCoefs  <- allSelCoefs[sortInfo$ix]
  epis = 2.01
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

set.seed(123)
seed <- sample(1e8:9e8, length(distances_bp))

all_results <- map2_dfr(
  distances_bp,
  seed,
  run_distance
)

write.table(all_results, "/Users/vidyadheeshkelkar/Desktop/Warmup_Project/H0_diffHap_dom11.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#######################
  
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
  
# ggplot(summary_s,aes(x = distance_mb, y = mean_s, color = Cross, group = paste(Cross,pos))) +
#   geom_line(size = 1.2) +
#   geom_point(size = 2.5) +
#   geom_errorbar(aes(ymin = mean_s - se_s,ymax = mean_s + se_s), width = 0.03) +
#   facet_wrap(~ Cross, scales = "free_y") +
#     labs(x = "Distance between targets (Mb)", y = "Selection response (s)", color = "Cross") +
#     theme_bw() +
#     theme(
#       strip.text = element_text(size = 14, face = "bold"),
#       axis.title = element_text(size = 14, face = "bold"),
#       axis.text  = element_text(size = 12),
#       legend.title = element_text(size = 13),
#       legend.text  = element_text(size = 12))


ggplot(subset(summary_s, Target=="Target 1"), aes(x = distance_mb, y = mean_s, color = Cross, 
                      group = interaction(Target, Cross))) +
  geom_line(size = 1) + 
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = mean_s - se_s, ymax = mean_s + se_s), width = 0.03) +
  labs(x = "Distance between targets (Mb)", y = "Selection response (s)", color = "Cross", linetype = "Target") +
  scale_color_manual(values = c("darkgreen", "orange1")) +
  # scale_fill_manual(values = c("darkgreen", "orange1")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text  = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12)
  )
  
# p4

p5 <- p1 + theme(axis.title.x = element_blank())
p6 <- p2 + theme(axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.y  = element_blank())
p7 <- p3 + theme(axis.title.y = element_blank(),
                 axis.text.y  = element_blank())

combined_plot <- ((p5|p6) / (p7|p8)) +
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

common_y <- scale_y_continuous(
  breaks = scales::pretty_breaks(n = 6),
  labels = scales::label_number(accuracy = 0.001)
)

p5 <- p1 + common_y
p6 <- p2 + common_y
p7 <- p3 + common_y
p8 <- p4 + common_y

combined_plot <- ((p5 | p6) / (p7 | p8)) +
  plot_layout(guides = "collect") &
  theme_bw() &
  theme(
    axis.title = element_blank(),
    legend.position = "bottom"
  ) +
  plot_annotation(
    left = "Selection response (s)",
    bottom = "Distance between targets (Mb)"
  )

t.test(s ~ Cross, data = all_results, var.equal = FALSE) 

ttest_results <- all_results %>%
  group_by(distance_mb) %>%
  filter(n_distinct(Cross) == 2) %>%   # safety check
  do(tidy(t.test(s ~ Cross, data = ., var.equal = FALSE))) %>%
  ungroup()

ttest_results <- ttest_results %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"),
         estimate = -1*estimate,
         significance = ifelse(p_adj <= 0.05, "Significant", "Not significant"),
         log_pval = -log10(p_adj))

# locus1_results <- all_results %>%
#   filter(pos == target1)
# locus1_summary <- locus1_results %>%
#   group_by(distance_mb, Cross) %>%
#   summarise(
#     mean_s = mean(s, na.rm = TRUE),
#     se_s   = sd(s, na.rm = TRUE) / sqrt(n()),
#     .groups = "drop"
#   )
# 
# ggplot(locus1_summary,
#        aes(x = distance_mb, y = mean_s, color = Cross, group = Cross)) +
#   geom_line(size = 1.2) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(ymin = mean_s - se_s, ymax = mean_s + se_s), width = 0.03) +
#   labs(x = "Distance between targets (Mb)", y = "Selection response at locus 1 (s)",
#        color = "Population") +
#   theme_bw() +
#   facet_wrap(~Cross, scales="free_y", nrow = 2)
#   theme(
#     axis.title = element_text(size = 14, face = "bold"),
#     axis.text  = element_text(size = 12),
#     legend.title = element_text(size = 13),
#     legend.text  = element_text(size = 12)
#   )
#   

ggboxplot(all_results, x="Target", y="s", color="Cross") +
  labs(x = "Position", y = "Selection Strength", color = "Population") +
  scale_x_discrete(labels = c("Target1", "Target2")) +
  stat_compare_means(aes(pos,s,group=Cross,label = after_stat(p.signif)), method = "t.test") +
  theme_bw()+
  theme(
    axis.text.x = element_text(vjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "cm"),
    legend.position = "bottom",
    strip.text.x = element_text(size = 12)
  )

box_data <- all_results %>%
  mutate(
    Target = factor(
      if_else(pos == target1, 1L, 2L),
      labels = c("Target 1", "Target 2")
    )
  )

ggplot(box_data, aes(x = Target, y = s, color = Cross)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  stat_compare_means(aes(Target,s,group=Cross,label = after_stat(p.signif)), method = "t.test") +
  labs(x = "Target", y = "Selection response (s)",color = "Cross") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text  = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12)
  )

ggplot(box_data,aes(x = Target, y = s, fill = Cross)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(aes(Target,s,group=Cross,label = after_stat(p.signif)), method = "t.test") +
  facet_wrap(~ distance_mb, nrow = 2) +
  labs(x = "Target", y = "Selection response (s)") +
  theme_bw()

selected_distances <- all_results %>%
  distinct(distance_mb) %>%
  arrange(distance_mb) %>%
  slice(round(seq(1, n(), length.out = 10))) %>%
  pull(distance_mb)
box_data_10 <- all_results %>%
  filter(distance_mb %in% selected_distances) %>%
  mutate(
    Target = factor(
      if_else(pos == target1, 1L, 2L),
      labels = c("Target 1", "Target 2"))
  )

ggplot(box_data_10,aes(x = Target, y = s, color = Cross)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8)) +
  stat_compare_means(aes(Target,s,group=Cross,label = after_stat(p.signif)), method = "t.test") +
  facet_wrap(~ distance_mb, nrow = 2) +
  labs(
    x = "Target",
    y = "Selection response (s)",
   color = "Cross"
  ) +
  theme_bw()


#####################


t <- read.table("/Users/vidyadheeshkelkar/Desktop/Warmup_Project/simsSO/targets5.txt", header = TRUE)
t <- subset(t, numTargets == 25)
# obsDtOS <- as.data.table(obsDtm)
# obsDtSO <- as.data.table(obsDtm)
allPos <- as.integer(unique(obsDt[, pos]))

allTargetPos <- round((t$targetPos))
allSelCoefs <- c(rep(0,length(allPos)))
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
# obsFreqs <- obsDtOS[pos %in% focalPos,.(pos,freq)]
# obsFreqs <- dcast(obsFreqs, pos ~ rep, value.var = "freq")
# obsFreqs <- obsFreqs[,-c(1)]
# transSelCoefs <- transForth(selCoefs)

# simFreqX3 <- getSimFreqs2(atPos=allPos, atGen=10L, targetPos, selCoefs, as.integer(popSize), seedBase = seed, numReps = 5L,
                          # saveSims = NULL, atGens = NULL)
# simFreqX4 <- getSimFreqs2(atPos=allPos, atGen=10L, targetPos, selCoefs, as.integer(popSize), seedBase = seed, numReps = 100L,
                          # saveSims = NULL, atGens = NULL)


OS_1 <- getSimFreqs2(atPos = allPos, atGen = 2L,
  targetPos = allPos, selCoefs = selCoefs, popSize = popSize, 
  initFreq = 0.67, epis = epis, domCoefs = domCoefs,
  seedBase = seed, numReps = 5L)

OS_10 <- getSimFreqs2(atPos = allPos, atGen = 11L,
                     targetPos = allPos, selCoefs = selCoefs, popSize = popSize, 
                     initFreq = 0.67, epis = epis, domCoefs = domCoefs,
                     seedBase = seed, numReps = 5L)

SO_1 <- getSimFreqs2(atPos = allPos, atGen = 2L,
                     targetPos = allPos, selCoefs = selCoefs, popSize = popSize, 
                     initFreq = 0.33, epis = epis, domCoefs = domCoefs,
                     seedBase = seed, numReps = 5L)

SO_10 <- getSimFreqs2(atPos = allPos, atGen = 11L,
                      targetPos = allPos, selCoefs = selCoefs, popSize = popSize, 
                      initFreq = 0.33, epis = epis, domCoefs = domCoefs,
                      seedBase = seed, numReps = 5L)

############

run_dom_sim <- function(seed) {
  
  targetPos <- target1
  selCoefs  <- 0.1
  epis <- 2.01
  popSize <- 3000 
  
  initFreqs <- c(0.05, 0.33, 0.5, 0.67, 0.90)
  domVals   <- c(0, 0.5, 1)
  gens      <- c(2,11)
  
  results <- list()
  idx <- 1
  
  for (f in initFreqs) {
    for (d in domVals) {
      for (g in gens) {
        
        sim <- getSimFreqs2(
          atPos = targetPos,
          atGen = g,
          targetPos = targetPos,
          selCoefs = selCoefs,
          popSize = popSize,
          initFreq = f,
          epis = epis,
          domCoefs = d,
          seedBase = seed,
          numReps = 1000L
        )
        
        sim <- sim %>%
          mutate(
            gen = g,
            initFreq = f,
            domCoef = d
          )
        
        results[[idx]] <- sim
        idx <- idx + 1
      }
    }
  }
  
  bind_rows(results)
}

slimScript <- "/Users/vidyadheeshkelkar/Downloads/toy03.slim"

seed <- sample(1e8:9e8,1)
df <- run_dom_sim(seed)

# df <- df %>%
#   group_by(initFreq, domCoef) %>%
#   summarise(mean_s = mean(s), .groups="drop")

# ggplot(subset(df,gen == 10), aes(domCoef, freq, group = domCoef)) +
  # geom_boxplot() +
  # facet_wrap(~initFreq)

df <- df %>%
  group_by(gen, initFreq, domCoef) %>%
  mutate(gen = gen - 1)  # 2 → 1, 3 → 2, ..., 11 → 10

library(ggplot2)
library(dplyr)

df_summary <- df %>%
  group_by(initFreq, domCoef, gen) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

selCoeff <- subset(df, gen %in% c(1,10))
selCoeff <- selCoeff %>% 
  pivot_wider(names_from = gen, values_from = freq) 

selCoeff <- selCoeff %>%
  mutate(s = (logit(`10`) - logit(`1`))/10)

ggplot(df_summary, aes(x = gen, y = mean_freq, color = factor(domCoef))) +
  geom_line(size = 1.2) +
  facet_wrap(~initFreq, scales = "free_y") +
  labs(
    x = "Generation",
    y = "Allele frequency",
    color = "Dominance coefficient",
    title = "Allele frequency trajectories for different initial frequencies and dominance"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

library(ggplot2)
library(dplyr)

ggplot(subset(df, gen %in% c(1,10))) +
  # ribbon for mean ± sd
  stat_summary(aes(x = gen, y = freq, group = paste(initFreq, factor(domCoef)), fill = factor(domCoef)),
                fun.data = "mean_sd", geom = "ribbon",alpha = 0.2) +
  stat_summary(aes(x = gen, y= freq, group = paste(initFreq, factor(domCoef)), color = factor(domCoef)),
               fun = "median", geom = "line") +
  # facet_wrap(~initFreq, scales="free") +
  labs(
    x = "Generation",
    y = "Allele frequency",
    color = "Dominance coefficient",
    fill = "Dominance coefficient",
    title = "Allele frequency trajectories with mean ± SD"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

write.table(selCoeff, "/Users/vidyadheeshkelkar/Desktop/Warmup_Project/dom_auto6.txt", sep = "\t", row.names = FALSE, quote = FALSE)

ggplot(selCoeff, aes(x = factor(initFreq), y = s, fill = factor(domCoef))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(x = "initFreq", y = "s", fill = "domCoef")





X <- read.table("/Users/vidyadheeshkelkar/Desktop/Warmup_Project/dom_X5.txt", sep = "\t", header = TRUE)
X$chrom <- "X"
A <- read.table("/Users/vidyadheeshkelkar/Desktop/Warmup_Project/dom_auto5.txt", sep = "\t", header = TRUE)
A$chrom <- "A"

X <- X %>%
  mutate(AFC = X10 - X1)

A <- A %>%
  mutate(AFC = X10 - X1)

df2 <- bind_rows(X,A)
df2<- df2 %>% 
  pivot_longer(cols = c(X1,X10),
             names_to = "gen",
             values_to = "freq")

df2$gen<- as.numeric(gsub("X", "", df2$gen))
df2$domCoef <- as.numeric(df2$domCoef)

dfX<- df2 %>% 
  select(-c(5,9)) %>%
  pivot_wider(names_from = domCoef, values_from = AFC,
              names_prefix = "h_")

dfX <- dfX %>%
  group_by(rep,pos,initFreq,chrom,gen) %>%
  mutate(h_0 = h_0 / h_0.5,
         h_1 = h_1 / h_0.5,
         h_0.5 = h_0.5 / h_0.5) %>%
  ungroup()

dfX <- dfX %>% 
  pivot_longer(cols = c(h_0,h_0.5,h_1),
               names_to = "domCoef",
               values_to = "rel_AFC")

p1 <- ggplot(df2) +
  # ribbon for mean ± sd
  stat_summary(aes(x = gen, y = freq, group = paste(initFreq, factor(domCoef)), fill = factor(domCoef)),
               fun.data = "mean_sd", geom = "ribbon",alpha = 0.2) +
  stat_summary(aes(x = gen, y= freq, group = paste(initFreq, factor(domCoef)), color = factor(domCoef)),
               fun = "median", geom = "line") +
  facet_grid(chrom~initFreq, scales="free",
             labeller = labeller(initFreq = c("0.05" = "5%", "0.33" = "33%", "0.5" = "50%",
               "0.67" = "67%", "0.9" = "90%"), chrom = c("A" = "Autosome"))) +
  labs(x = "Generation", y = "Allele frequency", color = "Dominance Coefficient",
    fill = "Dominance Coefficient") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.position = "bottom"
  ) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Initial Frequency", breaks = NULL, labels = NULL), breaks = c(1,5,10))+
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Chromosome", breaks = NULL, labels = NULL), limits = c(0,1))


codom <- subset(df2, domCoef == 0.5)

dfX$domCoef <- as.numeric(gsub("h_", "", dfX$domCoef))

my_comparisons <- list(c("0", "0.5"), c("0.5", "1"))

p2 <- ggplot(subset(dfX, gen==10),
       aes(x = factor(domCoef), y = rel_AFC, group = factor(domCoef))) +
# ggplot(dfX) + 
  geom_boxplot() + 
  # stat_summary(aes(x = gen, y = rel_AFC, group = paste(initFreq, factor(domCoef)), fill = factor(domCoef)),
  #              fun.data = "mean_sd", geom = "ribbon",alpha = 0.2) +
  # stat_summary(aes(x = gen, y= rel_AFC, group = paste(initFreq, factor(domCoef)), color = factor(domCoef)),
  #              fun = "median", geom = "line") +
  facet_grid(chrom~initFreq, scales="free",
             labeller = labeller(initFreq = c("0.05" = "5%", "0.33" = "33%", "0.5" = "50%",
                                              "0.67" = "67%", "0.9" = "90%"), chrom = c("A" = "Autosome"))) +
  labs(x = "Dominance Coefficient", y = "Relative Allele Frequency Change", color = "Dominance Coefficient") +
  # scale_color_manual(values = c("red", "green", "blue")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Chromosome", breaks = NULL, labels = NULL))

# label.x = 1.4, label.y = 0.075, size = 6

(p1/p2) + 
  plot_layout(widths = c(2,2), guides = "collect") &
  # labs(color = "New legend title") +
  plot_annotation(tag_levels = "A") & #list(c("B", "C", "D", "E"))
  theme_bw() &
  theme(
    legend.key.size = unit(0.5, "cm"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 12),
    plot.title = element_text(size = 7, face = "bold"),
    legend.position = "bottom",
  )

