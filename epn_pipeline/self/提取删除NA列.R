script_dir = "D:/R_Neo/"
source(paste0( script_dir ,"daoru.R"),encoding="utf-8")
input_file <- read.xlsx(test_xlsx,1,header = TRUE,encoding = "UTF-8")
del1 <- which(is.na(input_file$测序编号)==TRUE)
input_file <- input_file[-del1,]
input_file$seq_id <- paste0("SEQ", input_file$测序编号)

