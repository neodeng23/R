setwd("/home/denglf/r_neo/test_res")
script_dir = "/home/denglf/r_neo/SciencePipeline/"
source(paste0( "/home/denglf/r_neo/self/daoru.R"),encoding="utf-8")
source(paste0( script_dir ,"plot.R"),encoding="utf-8")
source(paste0( script_dir ,"SciencePipeline.R"),encoding="utf-8")
source(paste0( script_dir ,"epnInputs.R"),encoding="utf-8")

get_init_C_input <- function(input_sql, disease_type = NULL ) {
  if ( is.null(disease_type) ){
    if ( is.null(input_sql$detail) ){
      print( "请输入disease_type" )
    } else{
      input_sql$type <- input_sql$detail
    }
  } else {
    input_sql$type <- disease_type
    input_sql$detail <- input_sql$type
  }
  del1 <- which(is.na(input_sql$测序编号)==TRUE)
  input_sql <- input_sql[-del1,]
  input_sql$seq_id <- paste0("SEQ", input_sql$测序编号)
  rownames(input_sql) <- input_sql$seq_id
  input_sql$feature_count_location <- paste0("/data_bak/rawdata/featurecounts/SEQ",input_sql$测序编号,"/SEQ",input_sql$测序编号,".genebody")
  output_col <- c( "seq_id" , "年龄" ,"type" , "detail" ,"feature_count_location" , "genebody_chr", "duplicate_ratio" , "文库浓度" )
  input <- input_sql[,output_col]
  names(input)[2] <- "age"
  names(input)[8] <- "lib_conc_ngul"
  return( input )
}   #中文翻译版
input_CAC <- get_init_C_input(input_file,disease_type ="CAC")

sql_statement2 <- "SELECT * FROM All_Sample_Info where diagnosis LIKE '健康';"

input_HEA_df <- get_sqloutput( user = "denglf" , password = "550908392" , dbname = "BioInfo", sql_statement2 )
input_HEA <- get_init_input(input_HEA_df, "HEA")
input_tot <- rbind( input_HEA, input_CAC )
input_tot <- input_tot[ input_tot$genebody_chr > 0.56 & input_tot$genebody_chr < 0.586 , ]  #筛选genebody_chr
input_CAC <- input_tot[ input_tot$type == "CAC" , ]
input_HEA <- input_tot[ input_tot$type == "HEA" , ]
countData <- get_featurecount(input_tot)              #获取countData

gene_profile_deseq( countData , input_tot , c( "HEA" , "CAC" ) , case_num = "first" , volcano_cutoff = 0.3 , sort_method = "log2FoldChange" )
