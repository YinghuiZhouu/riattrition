# code to prepare `bursztyn_et_al_2020` dataset

# Step 1: run bursztyn_et_al_2020/bursztyn_et_al_2020aer_clean.do
# Output: bursztyn_et_al_2020/data/clean_data/02_follow_up_clean.dta

# Step 2: run the following R code
# Output: bursztyn_et_al_2020.rda

library(haven)
library(tidyverse)

bursztyn_et_al_2020 <- read_dta("data-raw/bursztyn_et_al_2020/data/clean_data/02_follow_up_clean.dta")
bursztyn_et_al_2020 <- bursztyn_et_al_2020 %>%
  zap_label() %>%
  as.data.frame()

# Generate Kling, Liebman, Katz (2007) index
index <- c(
  "employed_now_out_fl",
  "applied_out_fl",
  "interviewed_out_fl",
  "interview_sched_out_fl",
  "driving_fl",
  "outside_others_fl2"
)

bursztyn_et_al_2020 <- bursztyn_et_al_2020 %>%
  mutate(
    across(
      all_of(index),
      ~ as.numeric(scale(.x)),
      .names = "std_{.col}"
    ),
    klk_index =
      (
        std_employed_now_out_fl +
          std_applied_out_fl +
          std_interviewed_out_fl +
          std_interview_sched_out_fl +
          std_driving_fl +
          std_outside_others_fl2
      ) / 6
  )

usethis::use_data(bursztyn_et_al_2020, overwrite = TRUE, compress = "xz")
