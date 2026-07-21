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
library("forcats")
library("emmeans")
library("pROC")

### Initialization of parameters
registerDoParallel(cores = parallel::detectCores())
outFileRds <- "path_to_output_directory/out_slim.rds"
numCores <- 12L
numReps <- 100L 
slimCmd <- "/usr/local/bin/slim"
slimScript <- "path_to_slim_script/single_chr.slim"
popSize <- 1250L
focalGen <- 10L

af_df_add<- list()
sel_df_add<-list()
af_df_int<- list()
sel_df_int<-list()
target_info<- data.frame(matrix(ncol=5, nrow=400))

getSimFreqs2 <- newGetSimFreqs2(slimCmd, slimScript)
genomeSize <- 32079331
base_seed<- 1234
set.seed(base_seed)
target_1 <- sample(1:genomeSize, 1, replace = F) 
allowed_2 <- setdiff(1:genomeSize, unique(target_1))
selCoefs_1 <- runif(1, min=0, max=0.1)
selCoefs_2 <- runif(1, min=0, max=0.1)
selCoefs<- c(selCoefs_1, selCoefs_2)
# sortInfo <- sort(allTargetPos, index.return = TRUE)
# targetPos <- sortInfo$x
# selCoefs  <- allSelCoefs[sortInfo$ix]
epis = 2.00
domCoefs <- c(0.5, 0.5)

### Run the simulation for 400 different positions of marker 2 ###

for(i in 1:400){
  
  set.seed(base_seed+i)
  target_2 <- sample(allowed_2, 1, replace = F)
  targetPos<- c(target_1, target_2)
  af_df_add[[i]]<-getSimFreqs2(
    atPos = targetPos, atGen = 11L,
    targetPos = targetPos, selCoefs = selCoefs,
    popSize = popSize, initFreq = 0.5,
    epis = 1.0, domCoefs = domCoefs,
    seedBase = (base_seed+i), numReps = 5L, run=i
  )
  af_df_int[[i]]<-getSimFreqs2(
    atPos = targetPos, atGen = 11L,
    targetPos = targetPos, selCoefs = selCoefs,
    popSize = popSize, initFreq = 0.5,
    epis = epis, domCoefs = domCoefs,
    seedBase = (base_seed+i), numReps = 5L, run=i
  )
  sel_df_add[[i]]<- logit_sel(af_df_add[[i]])
  sel_df_int[[i]]<- logit_sel(af_df_int[[i]])
  
  target_info$run[i]<- paste0(i)
  target_info$pos1[i] <- target_1-1
  target_info$pos2[i] <- target_2-1
  target_info$coef1[i] <- selCoefs_1
  target_info$coef2[i] <- selCoefs_2
  print(i)
}

### Combine and filter the results ###
all_runs_af_add<- bind_rows(af_df_add)
all_runs_af_add$arch<- paste("add")
all_runs_sel_add<- bind_rows(sel_df_add)
all_runs_sel_add$arch<- paste("add")
all_runs_af_int<- bind_rows(af_df_int)
all_runs_af_int$arch<- paste("non-add")
all_runs_sel_int<- bind_rows(sel_df_int)
all_runs_sel_int$arch<- paste("non-add")
all_runs_af<- rbind(all_runs_af_add, all_runs_af_int)
all_runs_sel<- rbind(all_runs_sel_add, all_runs_sel_int)


df_target_1_add<- target_info%>%dplyr::select(run, pos1, coef1)%>%dplyr::rename("pos"=pos1, "s"=coef1)%>%
  mutate("arch"="add")
df_target_1_int<- df_target_1_add%>%mutate("arch"="non-add")
df_target_1<- rbind(df_target_1_add, df_target_1_int)
df_target_2_add<- target_info%>%dplyr::select(run, pos2, coef2)%>%dplyr::rename("pos"=pos2, "s"=coef2)%>%
  mutate("arch"="add")
df_target_2_int<- df_target_2_add%>%mutate("arch"="non-add")
df_target_2<- rbind(df_target_2_add, df_target_2_int)
sel_sample1<- all_runs_sel%>%filter(group=="p1")
sel_sample2<- all_runs_sel%>%filter(group=="p2")
sel_sample3<- all_runs_sel%>%filter(group=="p3")

### Hypothesis testing: s_1(AB)=s_1(AC)-s_1(BC) ###
test1_lhs<- sel_sample2%>%dplyr::filter(line=="A")
test1_rhs<- sel_sample3%>%dplyr::filter(line=="B")
data_pred_1<- rbind(test1_lhs, test1_rhs)
data_pred_1<- data_pred_1%>%dplyr::filter(pos %in% unique(df_target_1$pos))
data_obs_1<- sel_sample1%>%dplyr::filter(line=="A")
data_obs_1<- data_obs_1%>%dplyr::filter(pos %in% unique(df_target_1$pos))
result_test_1<-compute_p(data_pred_1, data_obs_1, "p2")
data_obs_1<- data_obs_1 %>%
  group_by(pos, run, arch) %>%
  summarise(val = mean(s, na.rm = TRUE))

df_target_info<- rbind(df_target_1, df_target_2)
df_result_calc<- result_test_1
df_af_sample<- all_runs_af
df_sel_sample<- all_runs_sel

### Visualization ###
pred_vs_obs_calc<- modify_result(df_result_calc, target_info) 
pred_vs_obs_calc$sign<- as.factor(pred_vs_obs_calc$sign)
p_dist_calc_add<- plot_sig(pred_vs_obs_calc%>%filter(arch=="add"), "Additive")
p_dist_calc_int<- plot_sig(pred_vs_obs_calc%>%filter(arch=="non-add"), "Epistasis")
print(p_dist_calc_add)
print(p_dist_calc_int)

s_dist_calc<- plot_diff_dist(pred_vs_obs_calc,"Calculated")
print(s_dist_calc)

roc_curve<- plot_roc(roc_df, 0.05)
print(roc_curve)
