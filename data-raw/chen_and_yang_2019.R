# code to prepare `chen_and_yang_2019` dataset

# Step 1: run chen_and_yang/chen_and_yang_2019aer_clean.do
# Output: chen_yang_2019aer_clean.dta

# Step 2: run the following R code
# Output: chen_and_yang_2019.rda

library(haven)
library(tidyverse)

chen_and_yang_2019 <- read_dta("data-raw/chen_and_yang_2019/chen_yang_2019aer_clean.dta")
chen_and_yang_2019 <- chen_and_yang_2019 %>%
  zap_label() %>%
  as.data.frame()

usethis::use_data(chen_and_yang_2019, overwrite = TRUE, compress = "xz")

