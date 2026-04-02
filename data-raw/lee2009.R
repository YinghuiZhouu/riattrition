## code to prepare `lee2009` dataset goes here

# Source of raw data: <https://www.openicpsr.org/openicpsr/project/113269/version/V1/view?flag=follow&pageSize=100&sortOrder=(?title)&sortAsc=true>
# Source of Lee (2009) replication code: <https://academic.oup.com/restud/article/76/3/1071/1590707?guestAccessKey=#supplementary-data>

# Step 1: run trimdata1.do, trimdata2.do, trimdata3.do.
# Output: trimdata3.dta

# Step 2: run the following R code
# Output: lee2009.rda
library(haven)
library(tidyverse)
jc <- read_dta("data-raw/trimdata3.dta")

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

lee2009 <- jc_wide %>%
  mutate(across(
    where(is.numeric),
    ~ ifelse(is.nan(.) | is.infinite(.), NA, .)
  ))

usethis::use_data(lee2009, overwrite = TRUE)

