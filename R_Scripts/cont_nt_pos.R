#IN USE
#MAKE SOEM CHANGES
####making continuous windows of nucelotide positions using ore SNP information.
#To be Edited: 02.09.25
library(dplyr)
library(tidyr)
a<- readRDS("~/Desktop/manuscript_files/out_600.rds")
# Make sure start and end are numeric
win_ore <- a[,"window"] %>%
  distinct()
win_ore$old_win<- win_ore$window
win_ore<- win_ore%>%
  separate_wider_delim(window, ":", names = c("chr", "win")) %>%
  separate_wider_delim(win, "-", names = c("start", "end")) %>%
  mutate(across(c(start, end), as.numeric))

# Initialize new columns
win_ore$new_start <- NA_integer_
win_ore$new_end <- NA_integer_

adjust_windows <- function(df_chr) {
  residual <- 0  # Reset for each chromosome
  for (i in 1:nrow(df_chr)) {
    if (i == 1) {
      df_chr[i, "new_start"] <- df_chr[i, "start"]
      mid <- (df_chr[i, "end"] + df_chr[i + 1, "start"]) / 2 + residual
      rounded_mid <- round(mid)
      df_chr[i, "new_end"] <- rounded_mid
      residual <- mid - rounded_mid
    } else if (i == nrow(df_chr)) {
      mid <- (df_chr[i - 1, "end"] + df_chr[i, "start"]) / 2 + residual
      rounded_mid <- round(mid)
      df_chr[i, "new_start"] <- rounded_mid
      df_chr[i, "new_end"] <- df_chr[i, "end"]
    } else {
      mid_start <- (df_chr[i - 1, "end"] + df_chr[i, "start"]) / 2 + residual
      new_start <- round(mid_start)
      residual <- mid_start - new_start
      
      mid_end <- (df_chr[i, "end"] + df_chr[i + 1, "start"]) / 2 + residual
      new_end <- round(mid_end)
      residual <- mid_end - new_end
      
      df_chr[i, "new_start"] <- new_start
      df_chr[i, "new_end"] <- new_end
    }
  }
  return(df_chr)
}

# Apply to each chromosome
win_ore_updated <- win_ore %>%
  group_split(chr) %>%
  map_df(adjust_windows)
win_ore<- win_ore_updated
win_ore$new_win<- paste(win_ore$chr, win_ore$new_start, sep=":")
win_ore$new_win<- paste(win_ore$new_win, win_ore$new_end, sep="-")
saveRDS(win_ore, file="~/Desktop/manuscript_files/window_details.rds")
a_updated <- a %>%
  left_join(win_ore %>% dplyr::select(old_win, new_win), by = c("window" = "old_win"))%>%as.data.frame()
a_updated <- a_updated %>%
  dplyr::select(-window) %>%
  rename(window = new_win)
a_new<- a_updated[,c(1, ncol(a_updated), 2: (ncol(a_updated)-1))]
saveRDS(a_new, "~/Desktop/Projects/new_SOH/data_window_analysis/win_ore.rds")
