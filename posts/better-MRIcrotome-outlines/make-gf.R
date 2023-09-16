
library(readxl)

# read the excel file describing the data
running <- read_excel("/well/lerch/users/yrf023/plasticity/analyses/Running_Data_Output_10s_Claire.xlsx")

# there's a bunch of weirdness in there to merge filenames and munge filenames
library(tidyverse)
gf <- data.frame(dirs=list.files("/well/lerch/users/yrf023/plasticity/plasticity-2023-03-10_processed/")) %>%
  mutate(id = str_replace(dirs, '(M\\d+)_.+', '\\1'),
         timepoint = str_replace(dirs, '.+w(\\d+)[-_].+', '\\1'),
         timepoint = case_when(
           str_detect(timepoint, "baseline") ~ -2,
           str_detect(timepoint, "1w") ~ 1,
           str_starts(timepoint, "1-") ~ 1,
           TRUE ~ as.numeric(timepoint)
         ),
         timepoint = ifelse(is.na(timepoint), 0, timepoint))

# merge the excel file and the dataframe based on the directories
gf <- left_join(gf %>% mutate(id = str_replace(id, '0119', '019')), 
                running %>% select(`animal id`, CuDi_Day13, sex), 
                by=c("id" = "animal id")) %>%
  mutate(group = ifelse(is.na(CuDi_Day13), "control", "running"),
         labels = paste0("/well/lerch/users/yrf023/plasticity/plasticity-2023-03-10_processed/", dirs, "/", dirs, "_N_I_lsq6_voted.mnc"),
         reljacs = paste0("/well/lerch/users/yrf023/plasticity/plasticity-2023-03-10_processed/", dirs, "/stats-volumes/", dirs, "_N_I_lsq6_lsq12_and_nlin__concat_inverted_linear_part_displ_log_det_rel_fwhm0.2.mnc"))

# just the week two timepoint for sake of simplicity
gfw2 <- gf %>% filter(timepoint==2)

write_csv(gfw2, file = "exercise-gfw2.csv")
