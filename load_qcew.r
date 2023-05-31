library(tidyverse)

qcewGetIndustryData <- function(year, qtr, industry) {
  url <- "http://data.bls.gov/cew/data/api/YEAR/QTR/industry/INDUSTRY.csv"
  url <- sub("YEAR", year, url, ignore.case = FALSE)
  url <- sub("QTR", tolower(qtr), url, ignore.case = FALSE)
  url <- sub("INDUSTRY", industry, url, ignore.case = FALSE)
  read.csv(url,
    header = TRUE,
    sep = ",", quote = "\"", dec = ".", na.strings = " ", skip = 0
  )
}

qcewGetAreaData <- function(year, qtr, area) {
  url <- "http://data.bls.gov/cew/data/api/YEAR/QTR/area/AREA.csv"
  url <- sub("YEAR", year, url, ignore.case = FALSE)
  url <- sub("QTR", tolower(qtr), url, ignore.case = FALSE)
  url <- sub("AREA", toupper(area), url, ignore.case = FALSE)
  read.csv(url, header = TRUE,
    sep = ",", quote = "\"", dec = ".", na.strings = " ", skip = 0)
}

## load data
year_list <- 2014:2021
quarter_list <- 1:4
## area_list <- c(24037, 24015, 24013, 24999)
a <- 24037
industry_list <- c(448, 4481, 4482, 4483)
## initialize data frame
qcew_data <- as.data.frame(matrix(nrow = length(industry_list),
  ncol = 2 + length(year_list) * 5))
colnames(qcew_data) <- c("naics", "area",
  paste0(rep("wage", length(year_list) * 5),
    rep(year_list - year_list[1] + 1, each = 5),
    rep("-", length(year_list) * 5),
    rep(c(1:4, "a"), length(year_list))))
## assign data
qcew_data$naics <- industry_list
qcew_data$area <- a
for (y in year_list) {
  for (q in quarter_list) {
    quarter_data <- qcewGetAreaData(
      as.character(y), as.character(q), as.character(a)) %>%
      filter(industry_code %in% industry_list) %>%
      select(area_fips, qtr, industry_code, total_qtrly_wages)
    qcew_data[, sprintf("wage%s-%s", y - year_list[1] + 1, q)] <-
      quarter_data$total_qtrly_wages
  }
  annual_data <- qcewGetAreaData(
    as.character(y), "A", as.character(a)) %>%
    filter(industry_code %in% industry_list) %>%
    select(area_fips, qtr, industry_code, total_annual_wages)
  qcew_data[, sprintf("wage%s-%s", y - year_list[1] + 1, "a")] <-
    annual_data$total_annual_wages
}

save(qcew_data, file = "./data/qcew_data")
