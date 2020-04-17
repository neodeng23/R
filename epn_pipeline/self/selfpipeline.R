setwd("/home/denglf/r_neo/test_res")
script_dir = "/home/denglf/r_neo/SciencePipeline/"
source(paste0( "/home/denglf/r_neo/self/daoru.R"),encoding="utf-8")
source(paste0( script_dir ,"plot.R"),encoding="utf-8")
source(paste0( script_dir ,"SciencePipeline.R"),encoding="utf-8")
source(paste0( script_dir ,"epnInputs.R"),encoding="utf-8")
input_CAC <- get_init_C_input(input_file,disease_type ="CAC")

#sql_statement1 <- "SELECT * FROM Colon_Sample_Info;"
#input_CAC_df <- get_sqloutput( user = "denglf" , password = "550908392" , dbname = "BioInfo", sql_statement1 )
#input_CAC <- get_init_input(input_CAC_df,disease_type ="CAC")

sql_statement2 <- "SELECT * FROM All_Sample_Info where diagnosis LIKE '健康';"
input_HEA_df <- get_sqloutput( user = "denglf" , password = "550908392" , dbname = "BioInfo", sql_statement2 )
input_HEA <- get_init_input(input_HEA_df, "HEA")

save( input_CAC , input_HEA , file = "conlon_input.rda")

input_tot <- rbind( input_HEA, input_CAC )
input_tot <- input_tot[ input_tot$genebody_chr > 0.56 & input_tot$genebody_chr < 0.586 , ]
input_CAC <- input_tot[ input_tot$type == "CAC" , ]
input_HEA <- input_tot[ input_tot$type == "HEA" , ]

countData <- get_featurecount(input_tot)
colData <-  data.frame("type" = input_tot[,c("type")] , row.names = rownames(input_tot) , check.names = FALSE )
countData <- countData[ rowSums(countData) > 10 , ]     # filter low counts
geoMeans <- apply( countData ,1,function(x) exp(mean(log(x))))
countData <- countData[ names(geoMeans),  ]
get_result_list(input_tot,"res_rda",model_type_set=c("HEA","CAC")) #耗时比较长

load("res_rda")
#################################################################
#质控图
#################################################################
plot_sample_cor(input_tot,countData,"tot.png",cutoff=0.97)        #correlation图
pca_plot( vst_all , input_tot, unique_id=NULL, "11.png", file_name_pr = NULL )  #PCA图
