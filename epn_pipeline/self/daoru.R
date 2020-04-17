library(xlsx)
test_xlsx <- "/home/denglf/r_neo/colon_project_sample_info.xlsx"
input_file <- read.xlsx(test_xlsx,1,header = TRUE,encoding = "UTF-8")

