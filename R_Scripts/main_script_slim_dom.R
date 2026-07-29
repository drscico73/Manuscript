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

### Include the functions ###
source("path_to_directory/func_slim_analysis.R")
source("path_to_directory/func_visual_slim.R")

### Initialization of parameters ###
registerDoParallel(cores = parallel::detectCores())
outFileRds <- "path_to_output_directory/out_slim.rds"
numCores <- 12L
slimCmd <- "/usr/local/bin/slim"
slimScript <- "path_to_slim_script/one_chr.slim"
popSize <- 1250L
focalGen <- 10L

af_df_rec_rec<- list()
af_df_dom_dom<- list()
af_df_dom_rec<- list()
af_df_rec_dom<- list()
af_df_codom_dom<- list()
af_df_codom_codom<- list()
af_df_dom_codom<- list()
af_df_rec_codom<- list()
af_df_codom_rec<- list()
sel_df_rec_rec<- list()
sel_df_dom_dom<- list()
sel_df_dom_rec<- list()
sel_df_rec_dom<- list()
sel_df_codom_dom<- list()
sel_df_dom_codom<- list()
sel_df_codom_codom<- list()
sel_df_rec_codom<- list()
sel_df_codom_rec<- list()
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
epis = 1.00

### Run the simulation for 400 different positions of marker 2 ###
for(i in 1:400){
  
  set.seed(base_seed+i)
  target_2 <- sample(allowed_2, 1, replace = F)
  targetPos<- c(target_1, target_2)
  af_df_rec_rec[[i]]<-getSimFreqs2(
    atPos = targetPos, atGen = 11L,
    targetPos = targetPos, selCoefs = selCoefs,
    popSize = popSize, initFreq = 0.5,
    epis = 1.0, domCoefs =c(0,0),
    seedBase = (base_seed+i), numReps = 5L, run=i
  )
  af_df_rec_dom[[i]]<-getSimFreqs2(
    atPos = targetPos, atGen = 11L,
    targetPos = targetPos, selCoefs = selCoefs,
    popSize = popSize, initFreq = 0.5,
    epis = 1.0, domCoefs =c(0,1),
    seedBase = (base_seed+i), numReps = 5L, run=i
  )
  af_df_dom_rec[[i]]<-getSimFreqs2(
    atPos = targetPos, atGen = 11L,
    targetPos = targetPos, selCoefs = selCoefs,
    popSize = popSize, initFreq = 0.5,
    epis = 1.0, domCoefs =c(1,0),
    seedBase = (base_seed+i), numReps = 5L, run=i
  )
  af_df_dom_dom[[i]]<-getSimFreqs2(
    atPos = targetPos, atGen = 11L,
    targetPos = targetPos, selCoefs = selCoefs,
    popSize = popSize, initFreq = 0.5,
    epis = 1.0, domCoefs =c(1,1),
    seedBase = (base_seed+i), numReps = 5L, run=i
  )
  af_df_codom_dom[[i]]<-getSimFreqs2(
    atPos = targetPos, atGen = 11L,
    targetPos = targetPos, selCoefs = selCoefs,
    popSize = popSize, initFreq = 0.5,
    epis = 1.0, domCoefs =c(0.5,1),
    seedBase = (base_seed+i), numReps = 5L, run=i
  )
  af_df_dom_codom[[i]]<-getSimFreqs2(
    atPos = targetPos, atGen = 11L,
    targetPos = targetPos, selCoefs = selCoefs,
    popSize = popSize, initFreq = 0.5,
    epis = 1.0, domCoefs =c(1,0.5),
    seedBase = (base_seed+i), numReps = 5L, run=i
  )
  af_df_codom_rec[[i]]<-getSimFreqs2(
    atPos = targetPos, atGen = 11L,
    targetPos = targetPos, selCoefs = selCoefs,
    popSize = popSize, initFreq = 0.5,
    epis = 1.0, domCoefs =c(0.5,0),
    seedBase = (base_seed+i), numReps = 5L, run=i
  )
  af_df_rec_codom[[i]]<-getSimFreqs2(
    atPos = targetPos, atGen = 11L,
    targetPos = targetPos, selCoefs = selCoefs,
    popSize = popSize, initFreq = 0.5,
    epis = 1.0, domCoefs =c(0,0.5),
    seedBase = (base_seed+i), numReps = 5L, run=i
  )
  af_df_codom_codom[[i]]<-getSimFreqs2(
    atPos = targetPos, atGen = 11L,
    targetPos = targetPos, selCoefs = selCoefs,
    popSize = popSize, initFreq = 0.5,
    epis = 1.0, domCoefs =c(0.5,0.5),
    seedBase = (base_seed+i), numReps = 5L, run=i
  )
  sel_df_rec_rec[[i]]<- logit_sel(af_df_rec_rec[[i]])
  sel_df_rec_codom[[i]]<- logit_sel(af_df_rec_codom[[i]])
  sel_df_codom_rec[[i]]<- logit_sel(af_df_codom_rec[[i]])
  sel_df_rec_dom[[i]]<- logit_sel(af_df_rec_dom[[i]])
  sel_df_dom_rec[[i]]<- logit_sel(af_df_dom_rec[[i]])
  sel_df_dom_codom[[i]]<- logit_sel(af_df_dom_codom[[i]])
  sel_df_dom_dom[[i]]<- logit_sel(af_df_dom_dom[[i]])
  sel_df_codom_dom[[i]]<- logit_sel(af_df_codom_dom[[i]])
  sel_df_codom_codom[[i]]<- logit_sel(af_df_codom_codom[[i]])
  
  target_info$run[i]<- paste0(i)
  target_info$pos1[i] <- target_1-1
  target_info$pos2[i] <- target_2-1
  target_info$coef1[i] <- selCoefs_1
  target_info$coef2[i] <- selCoefs_2
  print(i)
}

### Combine and filter the results ###

all_runs_af_rec_rec<- bind_rows(af_df_rec_rec)
all_runs_af_rec_rec$arch<-paste("rec-rec")
all_runs_af_rec_dom<- bind_rows(af_df_rec_dom)
all_runs_af_rec_dom$arch<-paste("rec-dom")
all_runs_af_dom_rec<- bind_rows(af_df_dom_rec)
all_runs_af_dom_rec$arch<-paste("dom-rec")
all_runs_af_rec_codom<- bind_rows(af_df_rec_codom)
all_runs_af_rec_codom$arch<-paste("rec-codom")
all_runs_af_codom_rec<- bind_rows(af_df_codom_rec)
all_runs_af_codom_rec$arch<-paste("codom-rec")
all_runs_af_codom_dom<- bind_rows(af_df_codom_dom)
all_runs_af_codom_dom$arch<-paste("codom-dom")
all_runs_af_dom_codom<- bind_rows(af_df_dom_codom)
all_runs_af_dom_codom$arch<-paste("dom-codom")
all_runs_af_codom_codom<- bind_rows(af_df_codom_codom)
all_runs_af_codom_codom$arch<-paste("codom-codom")
all_runs_af_dom_dom<- bind_rows(af_df_dom_dom)
all_runs_af_dom_dom$arch<-paste("dom-dom")

all_runs_sel_rec_rec<- bind_rows(sel_df_rec_rec)
all_runs_sel_rec_rec$arch<-paste("rec-rec")
all_runs_sel_rec_dom<- bind_rows(sel_df_rec_dom)
all_runs_sel_rec_dom$arch<-paste("rec-dom")
all_runs_sel_dom_rec<- bind_rows(sel_df_dom_rec)
all_runs_sel_dom_rec$arch<-paste("dom-rec")
all_runs_sel_rec_codom<- bind_rows(sel_df_rec_codom)
all_runs_sel_rec_codom$arch<-paste("rec-codom")
all_runs_sel_codom_rec<- bind_rows(sel_df_codom_rec)
all_runs_sel_codom_rec$arch<-paste("codom-rec")
all_runs_sel_codom_dom<- bind_rows(sel_df_codom_dom)
all_runs_sel_codom_dom$arch<-paste("codom-dom")
all_runs_sel_dom_codom<- bind_rows(sel_df_dom_codom)
all_runs_sel_dom_codom$arch<-paste("dom-codom")
all_runs_sel_codom_codom<- bind_rows(sel_df_codom_codom)
all_runs_sel_codom_codom$arch<-paste("codom-codom")
all_runs_sel_dom_dom<- bind_rows(sel_df_dom_dom)
all_runs_sel_dom_dom$arch<-paste("dom-dom")

all_runs_af<- rbind(all_runs_af_rec_rec,all_runs_af_rec_dom,all_runs_af_dom_rec,
                    all_runs_af_rec_codom,all_runs_af_codom_rec,all_runs_af_codom_dom,
                    all_runs_af_dom_codom,all_runs_af_dom_dom,all_runs_af_codom_codom)
all_runs_sel<- rbind(all_runs_sel_rec_rec,all_runs_sel_rec_dom,all_runs_sel_dom_rec,
                    all_runs_sel_rec_codom,all_runs_sel_codom_rec,all_runs_sel_codom_dom,
                    all_runs_sel_dom_codom,all_runs_sel_dom_dom, all_runs_sel_codom_codom)

df_target_1_rec_rec<- target_info%>%dplyr::select(run, pos1, coef1)%>%dplyr::rename("pos"=pos1, "s"=coef1)%>%
  mutate("arch"="rec-rec")
df_target_1_rec_dom<- target_info%>%dplyr::select(run, pos1, coef1)%>%dplyr::rename("pos"=pos1, "s"=coef1)%>%
  mutate("arch"="rec-dom")
df_target_1_dom_rec<- target_info%>%dplyr::select(run, pos1, coef1)%>%dplyr::rename("pos"=pos1, "s"=coef1)%>%
  mutate("arch"="dom-rec")
df_target_1_dom_dom<- target_info%>%dplyr::select(run, pos1, coef1)%>%dplyr::rename("pos"=pos1, "s"=coef1)%>%
  mutate("arch"="dom-dom")
df_target_1_rec_codom<- target_info%>%dplyr::select(run, pos1, coef1)%>%dplyr::rename("pos"=pos1, "s"=coef1)%>%
  mutate("arch"="rec-codom")
df_target_1_codom_rec<- target_info%>%dplyr::select(run, pos1, coef1)%>%dplyr::rename("pos"=pos1, "s"=coef1)%>%
  mutate("arch"="codom-rec")
df_target_1_codom_dom<- target_info%>%dplyr::select(run, pos1, coef1)%>%dplyr::rename("pos"=pos1, "s"=coef1)%>%
  mutate("arch"="codom-dom")
df_target_1_codom_codom<- target_info%>%dplyr::select(run, pos1, coef1)%>%dplyr::rename("pos"=pos1, "s"=coef1)%>%
  mutate("arch"="codom-codom")
df_target_1_dom_codom<- target_info%>%dplyr::select(run, pos1, coef1)%>%dplyr::rename("pos"=pos1, "s"=coef1)%>%
  mutate("arch"="dom-codom")

df_target_1<- rbind(df_target_1_rec_rec, df_target_1_dom_codom, df_target_1_codom_codom,
                    df_target_1_codom_dom, df_target_1_codom_rec, df_target_1_rec_codom, 
                    df_target_1_dom_dom, df_target_1_dom_rec, df_target_1_rec_dom)

df_target_2_rec_rec<- target_info%>%dplyr::select(run, pos2, coef2)%>%dplyr::rename("pos"=pos2, "s"=coef2)%>%
  mutate("arch"="rec-rec")
df_target_2_rec_dom<- target_info%>%dplyr::select(run, pos2, coef2)%>%dplyr::rename("pos"=pos2, "s"=coef2)%>%
  mutate("arch"="rec-dom")
df_target_2_dom_rec<- target_info%>%dplyr::select(run, pos2, coef2)%>%dplyr::rename("pos"=pos2, "s"=coef2)%>%
  mutate("arch"="dom-rec")
df_target_2_dom_dom<- target_info%>%dplyr::select(run, pos2, coef2)%>%dplyr::rename("pos"=pos2, "s"=coef2)%>%
  mutate("arch"="dom-dom")
df_target_2_rec_codom<- target_info%>%dplyr::select(run, pos2, coef2)%>%dplyr::rename("pos"=pos2, "s"=coef2)%>%
  mutate("arch"="rec-codom")
df_target_2_codom_rec<- target_info%>%dplyr::select(run, pos2, coef2)%>%dplyr::rename("pos"=pos2, "s"=coef2)%>%
  mutate("arch"="codom-rec")
df_target_2_codom_dom<- target_info%>%dplyr::select(run, pos2, coef2)%>%dplyr::rename("pos"=pos2, "s"=coef2)%>%
  mutate("arch"="codom-dom")
df_target_2_codom_codom<- target_info%>%dplyr::select(run, pos2, coef2)%>%dplyr::rename("pos"=pos2, "s"=coef2)%>%
  mutate("arch"="codom-codom")
df_target_2_dom_codom<- target_info%>%dplyr::select(run, pos2, coef2)%>%dplyr::rename("pos"=pos2, "s"=coef2)%>%
  mutate("arch"="dom-codom")

df_target_2<- rbind(df_target_2_rec_rec, df_target_2_dom_codom, df_target_2_codom_codom,
                    df_target_2_codom_dom, df_target_2_codom_rec, df_target_2_rec_codom, 
                    df_target_2_dom_dom, df_target_2_dom_rec, df_target_2_rec_dom)

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
pred_vs_obs_calc<- modify_result(result_test_1, target_info)
pred_vs_obs_calc$sign<- as.factor(pred_vs_obs_calc$sign)

### Get the number of predicted non-additive interactions (dominance) ###
table(pred_vs_obs_calc$sign)
