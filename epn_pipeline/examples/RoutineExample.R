######## 常规差异分析示例程序1（分析结直肠癌中右半比左半之间有无差异）
######## 通过sql_file1、sql_file2控制两类样本的取样，需将ColonPositon_right.txt、ColonPositon_left.txt文件放至工作目录
rm(list = ls())
######## 修改为自己的工作目录
setwd("/home/denglf/r_neo/test_res")
getwd()
source("/home/denglf/backup/science_service/SciencePipeline/epnInputs.R")
source("/home/denglf/backup/science_service/SciencePipeline/plot.R")
source("/home/denglf/backup/science_service/SciencePipeline/SciencePipeline.R")
source("/home/denglf/backup/science_service/SciencePipeline/ScoreUtils.R")
# sql_statement1为从数据库取出所有患病样本（第一类将定义为右半）
sql_statement1 <- "SELECT seq_num, featurecount_folder,age,gender FROM BioInfo.All_Sample_Info"
# sql_statement2同sql_statement1，即从数据库取出所有患病样本，（第二类将定义为左半）
sql_statement2 <- sql_statement1
# 右半RCAC（label为0）,左半LCAC（label为1）
# 修改user和password参数
routine_analysis(user = "denglf",password = "550908392",sql_statement1,sql_statement2,type_list=c("RCAC","LCAC"),
                 sql_file1="ColonPositon_right.txt",sql_file2="ColonPositon_left.txt",case_num = "colon_position",cutoff_list=c(0.05,0.05,0.26),
                 flag_list=c(FALSE,TRUE),up_down=c(16,29),volcano_cutoff = 0.1,GO_flag="all",pathway_num=20,pvalueCut=0.9,qvalueCut=0.9,cor_cutoff=0.95)


######## 常规差异分析示例程序2（分析肺的癌旁组织比癌组织之间有无差异）
rm(list = ls())
######## 修改为自己的工作目录
setwd("/home/denglf/r_neo/test_res")
getwd()
source("/home/denglf/backup/science_service/SciencePipeline/epnInputs.R")
source("/home/denglf/backup/science_service/SciencePipeline/plot.R")
source("/home/denglf/backup/science_service/SciencePipeline/SciencePipeline.R")
source("/home/denglf/backup/science_service/SciencePipeline/ScoreUtils.R")
# sql_statement1为从数据库取出癌旁组织TI（轻症）
sql_statement1 <- "SELECT seq_num,featurecount_folder FROM BioInfo.gDNA_Sample_Info
                          WHERE seq_num IN(11319,11320,11321,11322,11416)"
# sql_statement2为从数据库取出癌组织TU（重症）
sql_statement2 <- "SELECT seq_num, featurecount_folder FROM BioInfo.gDNA_Sample_Info
                          WHERE seq_num IN(11298,11324,11325,11326,11327)"
# 癌旁组织TI（轻症label为0）,癌组织TU（重症label为1）
# 修改user和password参数
routine_analysis(user = "denglf",password = "550908392",sql_statement1,sql_statement2,type_list=c("TI","TU"),
                 sql_file1=NULL,sql_file2=NULL,case_num = "lung",cutoff_list=c(0.05,0.05,0.26),
                 flag_list=c(FALSE,TRUE),up_down=c(100,0),volcano_cutoff = 0.26,GO_flag="all",pathway_num=20,pvalueCut=0.9,qvalueCut=0.9,cor_cutoff=0.95)
case_num = "lung"
load(paste0("input_list_",case_num,".rda"))
load(paste0("deseq_res_",case_num,".rda"))
load(paste0("gene_list_",case_num,".rda"))

