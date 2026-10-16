# code to prepare `noy_and_zhang_2023` dataset

noy_and_zhang_2023 <- read.csv("data-raw/noy_and_zhang_2023/noy_zhang_linconv_clean.csv")

usethis::use_data(noy_and_zhang_2023, overwrite = TRUE, compress = "xz")
