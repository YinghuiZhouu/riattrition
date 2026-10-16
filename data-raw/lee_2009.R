# code to prepare `lee_2009` dataset

# Source of raw data: <https://www.openicpsr.org/openicpsr/project/113269/version/V1/view>
# Source of Lee (2009) replication code: <https://academic.oup.com/restud/article/76/3/1071/1590707?guestAccessKey=#supplementary-data>

# Step 1: run lee_2009/trimdata1.do, lee_2009/trimdata2.do, lee_2009/trimdata3.do.
# Output: trimdata3.dta

# Step 2: run the following R code
# Output: lee_2009.rda

library(haven)
library(tidyverse)
jc <- read_dta("data-raw/lee_2009/trimdata3.dta")
jc <- jc %>%
  zap_label() %>%
  as.data.frame()

jc_treat <- jc %>%
  group_by(mprid) %>%
  summarize(treatmnt = first(treatmnt), .groups = "drop")

jc_wide <- jc %>%
  mutate(lhrwage = log(earnh / hwh)) %>%
  select(mprid, week, lhrwage) %>%
  pivot_wider(
    id_cols = mprid,
    names_from = week,
    values_from = c(lhrwage),
    names_glue = "{.value}{week}"
  )

jc_wide <- left_join(jc_treat, jc_wide, by = "mprid")

jc_wide <- jc_wide %>%
  mutate(across(
    where(is.numeric),
    ~ ifelse(is.nan(.) | is.infinite(.), NA, .)
  ))

lee_2009 <- jc_wide

usethis::use_data(lee_2009, overwrite = TRUE, compress = "xz")
