del1 <- which(is.na(input_file$测序编号)==TRUE)
input_file <- input_file[-del1,]

write.table (input_HEA_df[1], file ="input_HEA.txt", row.names =FALSE, col.names =FALSE, quote =FALSE)
