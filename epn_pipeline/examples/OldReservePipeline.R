# 保留pipeline示例1

# 获取script路径，加载相应文件
script_dir = "~/work_data/SciencePipeline/SciencePipeline/"
source(paste0( script_dir ,"SciencePipeline.R"))
source(paste0( script_dir ,"epnInputs.R"))
source(paste0( script_dir ,"ScoreUtils.R"))
source(paste0( script_dir ,"plot.R"))

# setwd("/data/group/huyl/DATA/LUNG_FS_PROJECT3")


set.seed(1)
# get inputs
sql_statement1 <- "SELECT * FROM Lung_Sample_FS_Project;"

input_CAC_df <- get_sqloutput( user = "huyl" , password = "fpHSibuaDo955h" , dbname = "BioInfo",
                            sql_statement1 )
input_CAC <- get_init_input(input_CAC_df)

sql_statement2 <- "SELECT * FROM All_Sample_Info where epn_id like 'SZ%' AND diagnosis LIKE '健康';"
input_HEA_df <- get_sqloutput( user = "huyl" , password = "fpHSibuaDo955h" , dbname = "BioInfo",
                            sql_statement2 )
input_HEA <- get_init_input(input_HEA_df, "HEA")

input_pathology <- get_pathology( input_CAC_df ,input_HEA_df  )

save( input_CAC , input_HEA , file = "original_input.rda")



case_num <-  1

# QC
model_type_set <- c( "HEA","CAC" )
input_NOR <- input_CAC[input_CAC$detail == "NOR",]
input_OTH <- input_CAC[input_CAC$detail == "OTHERs",]
input_CAC_OTHER <- input_CAC[ input_CAC$detail %in% c( "CAC_t0" , "CAC_t1" ) , ]
input_CAC <- input_CAC[input_CAC$detail %in% c( "CAC_t2" , "CAC_t25" ) == TRUE,]
input_CAC$detail <- "CAC"
input_CAC$type <- "CAC"
input_HEA <- input_HEA[ (dim(input_HEA)[1]-200):dim(input_HEA)[1] , ]



input_tot <- rbind( input_HEA, input_CAC )
input_tot <- input_tot[ input_tot$genebody_chr > 0.56 & input_tot$genebody_chr < 0.586 , ]
input_CAC <- input_tot[ input_tot$type == "CAC" , ]
input_HEA <- input_tot[ input_tot$type == "HEA" , ]

input <- rbind( input_CAC[ sample(dim(input_CAC)[1], size = 100 ) , ] , input_HEA[ sample(dim(input_HEA)[1] , size = 100 ) , ] )
input$detail <- factor( input$detail , levels = model_type_set )


input_validation <- input_tot[ rownames(input_tot) %in% rownames( input ) == FALSE , ]
input_validation <- input_validation[ input_validation$seq_id %in% c("SEQ10793", "SEQ11521", "SEQ9498", 
                                                                     "SEQ10187", "SEQ10065", "SEQ11676", 
                                                                     "SEQ8858", "SEQ8850", "SEQ11193", "SEQ9773", 
                                                                     "SEQ11868", "SEQ8807", "SEQ8759", "SEQ10028", 
                                                                     "SEQ11184", "SEQ10691", "SEQ11179", "SEQ10749", 
                                                                     "SEQ11675", "SEQ8941") == FALSE , ]
# get correlation
plotCor( input , "LUNG_FS_project" , 0.98 )


# cutoff/lof
input <- input[ input$seq_id %in% c( "SEQ7338" , "SEQ7363" , "SEQ11748")  == FALSE , ]

# fc / rpm
countData <- get_featurecount(input)
colData <- input$detail
gene_profile_deseq(countData = countData , input = input ,model_type_set = model_type_set ,case_num = case_num , sort_method = "padj" , volcano_cutoff = 0.15 )

load("/data/group/huyl/DATA/LUNG_FS_PROJECT3/cvfit_model.rda")
load("/data/group/huyl/DATA/LUNG_FS_PROJECT3/enrichkegg.rda")
load("/data/group/huyl/DATA/LUNG_FS_PROJECT3/gene_result.rda")
load("/data/group/huyl/DATA/LUNG_FS_PROJECT3/deseq_res1.rda")
# cvfit
score <- get_score( cvfit , dds , geoMeans , gene_list , input_validation)
# 带上病理信息的df
write.table( merge(score ,input_pathology , by = 0) , sep = ",", file = "validation_score.csv" )
score$score <- score$X1
auc_plot( score , data_name = "LUNG_FS" , plot_name = "LUNG_FS", file_name = "LUNG_FS.pdf")

score_CAC_OTHER <- get_score( cvfit , dds , geoMeans , gene_list , input_CAC_OTHER  )
write.table( merge( score_CAC_OTHER ,input_pathology , by = 0) , sep = ",", file = "validation_CAC_OTHER.csv" ) 

score_NOR <- get_score( cvfit , dds , geoMeans , gene_list , input_NOR , filename = "validation_NOR.csv" )
write.table( merge( score_NOR ,input_pathology , by = 0) , sep = ",", file = "validation_NOR.csv" ) 

score_OTH <- get_score( cvfit , dds , geoMeans , gene_list , input_OTH , filename = "validation_OTH.csv" )
write.table( merge( score_OTH ,input_pathology , by = 0) , sep = ",", file = "validation_OTH.csv" ) 


# 保留pipeline示例2

# 获取script路径，加载相应文件
script_dir = "~/work_data/SciencePipeline/SciencePipeline/"
source(paste0( script_dir ,"SciencePipeline.R"))
source(paste0( script_dir ,"epnInputs.R"))
source(paste0( script_dir ,"ScoreUtils.R"))
source(paste0( script_dir ,"plot.R"))

# 设置工作路径
setwd("/data/group/huyl/DATA/LUNG_FS_PROJECT2")

set.seed(4096)
# get inputs

sql_statement1 <- "SELECT * FROM Lung_Sample_FS_Project;"

input_CAC <- get_sqloutput( user = "huyl" , password = "fpHSibuaDo955h" , dbname = "BioInfo",
                            sql_statement1 )
input_CAC <- get_init_input(input_CAC)

sql_statement2 <- "SELECT * FROM All_Sample_Info where epn_id like 'SZ%' AND diagnosis LIKE '健康';"
input_HEA <- get_sqloutput( user = "huyl" , password = "fpHSibuaDo955h" , dbname = "BioInfo",
                            sql_statement2 )
input_HEA <- get_init_input(input_HEA, "HEA")

save( input_CAC , input_HEA , file = "original_input.rda")

case_num <-  1
model_type_set <- c( "HEA","CAC_t2" )
input_NOR <- input_CAC[input_CAC$detail == "NOR",]
input_OTH <- input_CAC[input_CAC$detail == "OTHERs",]
input_CAC <- input_CAC[input_CAC$detail == "CAC_t2",]
input_HEA <- input_HEA[ (dim(input_HEA)[1]-150):dim(input_HEA)[1] , ]
input <- rbind( input_CAC[ sample(dim(input_CAC)[1], size = 100 ) , ] , input_HEA[ sample(dim(input_HEA)[1] , size = 100 ) , ] )
input$detail <- factor( input$detail , levels = model_type_set )

input_tot <- rbind( input_HEA, input_CAC )
input_validation <- input_tot[ rownames(input_tot) %in% rownames( input ) == FALSE , ]
# get correlation
plotCor( input , "LUNG_FS_project" , 0.98 )


# cutoff/lof
input <- input[ input$seq_id %in% c( "SEQ7338" , "SEQ7363" , "SEQ11748")  == FALSE , ]

# fc / rpm
countData <- get_featurecount(input)
colData <- input$detail
gene_profile_deseq(countData = countData , input = input ,model_type_set = model_type_set , case_num = case_num , sort_method = "padj" )

# cvfit
score <- get_score( cvfit , dds , geoMeans , gene_list , input_validation, filename = "validation_score.csv" )
score$score <- score$X1
auc_plot( score , data_name = "LUNG_FS" , plot_name = "LUNG_FS", file_name = "LUNG_FS.pdf")

score_NOR <- get_score( cvfit , dds , geoMeans , gene_list , input_NOR , filename = "validation_NOR.csv" )
score_OTH <- get_score( cvfit , dds , geoMeans , gene_list , input_OTH , filename = "validation_OTH.csv" )

# 保留pipeline示例3

# 获取script路径，加载相应文件
script_dir = "~/work_data/SciencePipeline/SciencePipeline/"
source(paste0( script_dir ,"SciencePipeline.R"))
source(paste0( script_dir ,"epnInputs.R"))
source(paste0( script_dir ,"ScoreUtils.R"))
source(paste0( script_dir ,"plot.R"))

setwd("~/DATA/LUNG_tissue")


set.seed(4096)
# get inputs

sql_statement1 <- "SELECT * FROM Analy_Tissue_Data where epn_id like 'FSATU%';"
input_CAC <- get_sqloutput( user = "huyl" , password = "password" , dbname = "BioInfo",
                            sql_statement1 )
input_CAC <- get_init_input(input_CAC , disease_type = "CAC")

sql_statement2 <- "SELECT * FROM Analy_Tissue_Data where epn_id like 'FSATI%';"
input_NOR <- get_sqloutput( user = "huyl" , password = "password" , dbname = "BioInfo",
                            sql_statement2 )
input_NOR <- get_init_input(input_NOR , disease_type = "NOR")
input <- rbind( input_NOR , input_CAC )
model_type_set <- c( "NOR" , "CAC" )
countData <- get_featurecount(input)

# colData: columns correspond to type condition, and rows correspond to the samples
colData <-  data.frame("type" = input[,c("type")] , row.names = rownames(input) , check.names = FALSE )
colData$type <- factor(colData$type , levels = model_type_set)
# filter low counts
countData <- countData[ rowSums(countData) > 10 , ]
geoMeans <- apply( countData ,1,function(x) exp(mean(log(x+0.1))))
countData <- countData[ names(geoMeans),  ]

# DESeq分析
# 返回genes pvalue, padj, log2FoldChange
dds <- DESeqDataSetFromMatrix( countData = countData , colData = colData, design = ~type )
dds <- DESeq(dds,parallel=TRUE)
res <- results(dds)
# 选择显著基因，过滤na log2FoldChange
resOrdered <- res[which(res$pvalue<0.05),]
resOrdered <- resOrdered[ !is.na(resOrdered$log2FoldChange),]
resOrdered <- resOrdered[order(resOrdered$pvalue), ]
save( countData , dds , res , resOrdered, input ,geoMeans , file = "deseq_res.rda" )
