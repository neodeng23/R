######## 建模差异分析示例程序1（分析中山医院健康与胃癌之间的差异）
######## 通过seq_file控制患病样本的取样，需将GastricExample_seq.txt文件放至工作目录
rm(list = ls())
######## 修改为自己的工作目录
setwd("~/science_service/develope/all_model/")
getwd()
source("~/science_service/sciencepipeline/SciencePipeline/epnInputs.R")
source("~/science_service/sciencepipeline/SciencePipeline/plot.R")
source("~/science_service/sciencepipeline/SciencePipeline/SciencePipeline.R")
source("~/science_service/sciencepipeline/SciencePipeline/ScoreUtils.R")
# sql_statement1为从数据库取出中山医院健康样本
sql_statement1 <- "SELECT seq_num,featurecount_folder,age,gender,stage FROM BioInfo.All_Sample_Info WHERE epn_id LIKE 'SZ%HEA'"
# sql_statement2为从数据库取出所有患病样本
sql_statement2 <- "SELECT seq_num, featurecount_folder,age,gender,stage FROM BioInfo.All_Sample_Info"
# sql_stage为从数据库取出患病的分期
#sql_stage <- "SELECT seq_num, stage FROM BioInfo.All_Sample_Info"
# 健康CAC（label为0）,患病CAC（label为1）
# 修改user和password参数
model_analysis(user = "wangbw",password = "password",sql_statement1,sql_statement2,sql_valid1=NULL,sql_valid2=NULL,seq_file="GastricExample_seq.txt",HEA_num=115,CAC_num = 115,
               seed_num = 1,type_list=c("HEA","CAC"),case_num = "gastric",cutoff_list=c(0.05,0.05,0.26),
               flag_list=c(TRUE,TRUE),up_down=c(100,0),volcano_cutoff = 0.3,GO_flag="all",pathway_num=20,pvalueCut=0.9,qvalueCut=0.9,cor_cutoff=0.95)

######## 建模差异分析示例程序2（分析中山医院健康与肾癌之间的差异）
rm(list = ls())
######## 修改为自己的工作目录，需将Kidney_seq.txt文件放至工作目录
setwd("~/science_service/develope/kidney_model/")
getwd()
source("~/science_service/sciencepipeline/SciencePipeline/epnInputs.R")
source("~/science_service/sciencepipeline/SciencePipeline/plot.R")
source("~/science_service/sciencepipeline/SciencePipeline/SciencePipeline.R")
source("~/science_service/sciencepipeline/SciencePipeline/ScoreUtils.R")
# 从数据库中山健康样并抽样取100例
SZ_HEA <- get_sqloutput(user = "wangbw",password = "password","SELECT seq_num,featurecount_folder,age,gender,stage FROM BioInfo.All_Sample_Info where epn_id like'SZ%HEA'")
set.seed(1)
SZ_HEA_sample <- SZ_HEA[sample(1:nrow(SZ_HEA),100),]
# 从100例健康中抽50个作为训练集，剩余50例作为验证集
SZ_HEA_model <- SZ_HEA_sample[sample(1:nrow(SZ_HEA_sample),50),]
SZ_HEA_valid <- SZ_HEA_sample[!SZ_HEA_sample$seq_num %in% SZ_HEA_model$seq_num,]
# 从所提供的肾癌中抽样取31例作为训练集，剩余55例作为验证集
SZ_kidney <- data.frame(seq_num=read.table("Kidney_seq.txt",header = TRUE,sep = '\n')[,1])
kidney_model_index <- sample(1:nrow(SZ_kidney),31)
SZ_kidney_model <- data.frame(seq_num=SZ_kidney[kidney_model_index,])
SZ_kidney_valid <- data.frame(seq_num=SZ_kidney[-kidney_model_index,])
# sql_statement1/sql_valid1为从数据库取出建模健康样本和验证集健康样本
sql_statement1 <- paste0("SELECT seq_num,featurecount_folder,age,gender,stage FROM BioInfo.All_Sample_Info WHERE seq_num IN ('",paste(SZ_HEA_model$seq_num,collapse = "','"),"')")
sql_valid1 <- paste0("SELECT seq_num,featurecount_folder,age,gender,stage FROM BioInfo.All_Sample_Info WHERE seq_num IN ('",paste(SZ_HEA_valid$seq_num,collapse = "','"),"')")
# sql_statement2/sql_valid2为从数据库取出建模肾癌样本和验证集肾癌样本
sql_statement2 <- paste0("SELECT seq_num,featurecount_folder,age,gender,stage FROM BioInfo.All_Sample_Info WHERE seq_num IN ('",paste(SZ_kidney_model$seq_num,collapse = "','"),"')")
sql_valid2 <- paste0("SELECT seq_num,featurecount_folder,age,gender,stage FROM BioInfo.All_Sample_Info WHERE seq_num IN ('",paste(SZ_kidney_valid$seq_num,collapse = "','"),"')")
# 健康HEA（label为0），患病CAC（label为1）
# 修改user和password参数
model_analysis(user = "wangbw",password = "password",sql_statement1,sql_statement2,sql_valid1,sql_valid2,seq_file=NULL,HEA_num=50,CAC_num = 31,
               seed_num = 1,type_list=c("HEA","CAC"),case_num = "kidney",cutoff_list=c(0.05,0.05,0.26),
               flag_list=c(FALSE,TRUE),up_down=c(41,47),volcano_cutoff = 0.26,GO_flag="all",pathway_num=20,pvalueCut=0.9,qvalueCut=0.9,cor_cutoff=0.95)