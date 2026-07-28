### Make a population file which will be used to initialise the populations in each iteration ###

initSlimSim <- function(targetPos1, selCoefs, popFile, popSize, initFreq, epis, domCoefs, seed, slimCmd, slimScript) {
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

### Run the simulation ###

runSlimSim <- function(targetPos1, selCoefs, popFile, seed, maxGen, writeOutput, popSize, initFreq, epis, domCoefs,
                       slimCmd, slimScript) {
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
  slimCmdLine <- paste(slimCmd, slimArgs, slimScript)
  system(slimCmdLine, intern = TRUE)
}

### Clean the raw data obtained from slim to retain only the relevant information: mutation frequencies in each population###

clean_data<- function(df, mut_start, ind_start, hap_start){
  cleaned_df<- df[(mut_start + 1):(hap_start - 1)] %>%
    as.data.frame() %>%
    dplyr::rename("line" = '.') %>%
    separate_wider_delim(cols="line",
                         delim = " ",
                         names = c("mut_id", "unk1", "mut_type", "pos", "selCoeff", "domCoeff",
                                   "pop", "extra", "count"),
                         too_many= "merge",
                         too_few = "align_start",
                         cols_remove = TRUE
    ) %>%
    filter(!is.na(mut_id) & mut_id != "Mutations:")%>%
    dplyr::mutate(across(count,as.integer), across(pos,as.integer), 
                  across(unk1, as.integer),across(selCoeff,as.numeric))%>%
    mutate(freq=count/2500)
  return(cleaned_df)
}

### Calculate marker AF ###

get_freq<- function(df_lines){
  mut_start <- grep("^Mutations:", df_lines)
  ind_start <- grep("^Individuals:", df_lines)
  hap_start <- grep("^Haplosomes:", df_lines)
  
  new_df <- clean_data(df_lines, mut_start, ind_start, hap_start)
}

### Get frequencies ###

newGetSimFreqs2 <- function(slimCmd, slimScript) {
  function(atPos, atGen, targetPos, selCoefs, popSize, initFreq, epis, domCoefs,
           seedBase = seed, numReps, run, saveSims = NULL, atGens = NULL) {
    
    # Temporary population file
    popFile <- tempfile(pattern = "pop", fileext = ".bin")
    
    # Initialize simulation
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
    
    # Run replicates in parallel
    results <- foreach(
      repId = seq_len(numReps),
      .combine = rbind,
      .export = c("runSlimSim", "get_freq", "clean_data")
    ) %dopar% {
      
      # 1️⃣ Run SLiM
      s_1 <- runSlimSim(
        targetPos1 = targetPos,
        selCoefs = selCoefs,
        popFile = popFile,
        seed = seedBase + repId,
        maxGen = atGen,
        writeOutput = if (is.null(atGens)) atGen else atGens,
        popSize = popSize,
        initFreq = initFreq,
        epis = epis,
        domCoefs = domCoefs,
        slimCmd = slimCmd,
        slimScript = slimScript
      )
      
      # 2️⃣ Parse mutation frequencies
      df_rep <- get_freq(s_1)
      
      # 3️⃣ Add replicate ID
      df_rep$rep <- repId
      
      # 4️⃣ Return dataframe
      df_rep
    }
    
    marker1<- targetPos[1]-1
    marker2<- targetPos[2]-1

    #Identify the source line for mutation in the population

    sample1<- results%>%filter(pop=="p1")
    sample1<- sample1%>%mutate(line=ifelse(pos %in% marker1 & mut_type!="m1", "1", 
                                           ifelse(pos %in% marker2 & mut_type!="m1", "2", 
                                                  ifelse(pos %in% marker2 & mut_type=="m1","1",
                                                         ifelse(pos %in% marker1 & mut_type=="m1", "2",NA)))))
    sample2<- results%>%filter(pop=="p2")
    sample2<- sample2%>%mutate(line=ifelse(pos %in% marker2 & mut_type=="m1", "1", 
                                           ifelse(pos %in% marker1 & mut_type=="m2", "1", 
                                                  ifelse(pos %in% marker1 & mut_type=="m4","3",
                                                         ifelse(pos %in% marker2 & mut_type=="m4", "3",NA)))))
    sample3<-results%>%filter(pop=="p3")
    sample3<- sample3%>%mutate(line=ifelse(pos %in% marker2 & mut_type=="m3", "2", 
                                           ifelse(pos %in% marker1 & mut_type=="m1", "2", 
                                                  ifelse(pos %in% marker1 & mut_type=="m4","3",
                                                         ifelse(pos %in% marker2 & mut_type=="m4", "3",NA)))))
    
    sample1_avg <- sample1 %>%
      group_by(pos,pop, rep, line) %>%
      summarise(freq = mean(freq), .groups = "drop")
    sample2_avg <- sample2 %>%
      group_by(pos,pop,  rep, line) %>%
      summarise(freq = mean(freq), .groups = "drop")
    sample3_avg <- sample3 %>%
      group_by(pos,pop, rep, line) %>%
      summarise(freq = mean(freq), .groups = "drop")
    
    df_sample<- rbind(sample1_avg, sample2_avg, sample3_avg)
    df_sample$run<- run
    
    # Cleanup
    if (file.exists(popFile)) file.remove(popFile)
    
    # Return combined results
    return(df_sample)
  }
}

### Calculating logit-transformed AF ###

logit_freq<- function(samp_df){
  line_a<- colnames(samp_df)[ncol(samp_df)-1]
  line_b<- colnames(samp_df)[ncol(samp_df)]
  col_1<- paste0("logit_", line_a)
  col_2<- paste0("logit_", line_b)
  temp_df<-samp_df%>%mutate(
    !!col_1 := log(.data[[line_a]] / .data[[line_b]])
  )%>%
    mutate(
      !!col_2 := log(.data[[line_b]] / .data[[line_a]])
    )
  temp_df<- temp_df%>%pivot_longer(cols =c(.data[[col_1]], .data[[col_2]]), 
                                   names_to = "line", values_to = "s")%>%
    separate_wider_delim(line, delim ="_", names = c("drop", "line"))%>%dplyr::select(!drop)%>%
    dplyr::select(!c(.data[[line_a]], .data[[line_b]]))%>%
    dplyr::rename("group"=pop)
  return(temp_df)
}

### Calculate response to selection ###

logit_sel<- function(df){
  sample1<- df%>%filter(pop=="p1")
  sample2<- df%>%filter(pop=="p2")
  sample3<- df%>%filter(pop=="p3")
  sample1_avg <- sample1 %>%
    group_by(pos,pop, run, rep, line) %>%
    summarise(freq = mean(freq), .groups = "drop")
  df_sample1<- sample1_avg%>%group_by(pos, pop,run, rep)%>%
    pivot_wider(names_from = line, values_from = freq)%>%ungroup()
  sample2_avg <- sample2 %>%
    group_by(pos,pop, run, rep, line) %>%
    summarise(freq = mean(freq), .groups = "drop")
  df_sample2<- sample2_avg%>%group_by(pos, pop,run, rep)%>%
    pivot_wider(names_from = line, values_from = freq)%>%ungroup()
  sample3_avg <- sample3 %>%
    group_by(pos,pop, run, rep, line) %>%
    summarise(freq = mean(freq), .groups = "drop")
  df_sample3<- sample3_avg%>%group_by(pos, pop,run, rep)%>%
    pivot_wider(names_from = line, values_from = freq)%>%ungroup()
  logit_sample1<- logit_freq(df_sample1)
  logit_sample2<- logit_freq(df_sample2)
  logit_sample3<- logit_freq(df_sample3)
  logit_sample<- rbind(logit_sample1, logit_sample2, logit_sample3)
  return(logit_sample)
}

### Classifies the result from the hypothesis testing ###

modify_result<- function(df, target){
  pred_df<- df%>%dplyr::select(pos, run, contrast, s_val, arch)%>%
    dplyr::rename("pred"=contrast, "obs"=s_val)%>%mutate(diff_s=pred-obs)
  pred_df<- pred_df%>%merge(target%>%dplyr::select(run, pos2), by="run")%>%
    mutate(across(pos, as.integer))%>%mutate(diff_pos=abs(pos-pos2))%>%dplyr::select(!pos2)
  pred_df<- pred_df%>%
    merge(df%>%dplyr::select(arch,run, adj_P, p_value), by=c("run",  "arch"))
  
  mod_df<- pred_df%>%mutate(sign=ifelse(adj_P<=0.05, "1", "0"))%>%
    mutate(across(sign, as.factor))
  return(mod_df)
}
