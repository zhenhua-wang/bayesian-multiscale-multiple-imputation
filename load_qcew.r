library(tidyverse)

qcewGetAreaData <- function(year, area, annual = FALSE) {
  if (annual) {
    url <- "https://data.bls.gov/cew/data/files/%s/csv/%s_annual_by_area.zip"
  } else {
    url <- "https://data.bls.gov/cew/data/files/%s/csv/%s_qtrly_by_area.zip"
  }
  url <- sprintf(url, year, year)
  save_dir <- file.path(".", "tmp")
  save_path <- file.path(save_dir, sprintf("%s.zip", year))
  dir.create(save_dir)
  download.file(url, destfile = save_path, method = "curl")
  csv_name <- grep(sprintf(".+%s.+", area), unzip(save_path, list = TRUE)$Name,
    ignore.case = TRUE, value = TRUE)
  unzip(save_path, files = csv_name, junkpaths = TRUE, exdir = save_dir)
  file.remove(save_path)
  csv_path <- file.path(save_dir, last(strsplit(csv_name, "/")[[1]]))
  csv <- read.csv(csv_path)
  file.remove(csv_path)
  unlink(save_dir, recursive = TRUE)
  return(csv)
}

## load data
year_list <- 1990:2022
quarter_list <- 1:4
save_dir <- "./data"
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
  quarter_data <- qcewGetAreaData(y, a, annual = FALSE) %>%
    filter(industry_code %in% industry_list) %>%
    select(area_fips, qtr, industry_code, total_qtrly_wages)
  for (q in quarter_list) {
    qcew_data[, sprintf("wage%s-%s", y - year_list[1] + 1, q)] <-
      quarter_data[quarter_data$qtr == q, "total_qtrly_wages"]
  }
  annual_data <- qcewGetAreaData(y, a, annual = TRUE) %>%
    filter(industry_code %in% industry_list) %>%
    select(area_fips, qtr, industry_code, total_annual_wages)
  qcew_data[, sprintf("wage%s-%s", y - year_list[1] + 1, "a")] <-
    annual_data$total_annual_wages
}
dir.create(save_dir)
save(qcew_data, file = file.path(save_dir, "qcew_data"))
