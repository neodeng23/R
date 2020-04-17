library(doParallel)
library(foreach)
library(RMySQL)

#####################################################################################
# 从mysql获取数据
#####################################################################################
get_sqloutput <- function(user, password, sql_statement,dbname = "BioInfo",host= "192.168.2.209",port = 6666){
  # Retrieve data from the database
  # Args:
  #   user: username of your sql
  #   password: passwd of your sql
  #   sql_statment: sql statment 
  #   dbname: default "bioinfo" ,OnePiece" or 'epican_v1'
  #   host: default "192.168.1.142"
  #   port: 6666
  # Returns:
  #   results of the query set as a data frame object
  
  dbcon<-dbConnect( MySQL(), host= host ,dbname = dbname , user = user, password = password , port = port ) # create a database connection object
  dbSendQuery( dbcon,"SET NAMES utf8" ) # encoding character utf8 
  input_sql <- dbGetQuery( dbcon, sql_statement ) # read data into dataframe
  dbDisconnect( dbcon ) # close connection
  
  return( input_sql ) 
}

#####################################################################################
# 数据整形函数，可将样本的质量指标分组用于绘制质控图（不在分析流程中）
#####################################################################################
shaping_data <- function(origin_data,file_type,isresearch=FALSE){
  # file_type:取值可为genebody/bin/chr/promoter
  plot_data <- origin_data[!is.na(origin_data$featurecount_folder),]
  plot_data$seq_id <- paste0("SEQ",plot_data$seq_num)
  if (isresearch){
    plot_data$feature_count_location <- paste0(plot_data$featurecount_folder,"_q","/",plot_data$seq_id,".",file_type)
  } else {
    plot_data$feature_count_location <- paste0(plot_data$featurecount_folder,"/",plot_data$seq_id,".",file_type)
    #plot_data$feature_count_location <- paste0("/data_bak/rawdata/featurecounts/",plot_data$seq_id,"/",plot_data$seq_id,".",file_type)
  }
  
  rownames(plot_data) <- plot_data$seq_id
  plot_data$binding_rate_group <- cut(plot_data$binding_rate,
                                      breaks = c(-Inf,4,5,6,7,8,9,10,Inf),
                                      labels = c("4-","4-5","5-6","6-7","7-8","8-9","9-10","10+"),
                                      right = FALSE)
  plot_data$elution_rate_group <- cut(plot_data$elution_rate,
                                      breaks =c(-Inf,10,11,12,13,Inf),
                                      lables =c("10-","10-11","11-12","12-13","13+"),
                                      right = FALSE)
  plot_data$enrichment_rate_group <- cut(plot_data$enrichment_rate,
                                         breaks = c(-Inf,4,5,6,7,8,9,10,Inf),
                                         lables = c("4-","4-5","5-6","6-7","7-8","8-9","9-10","10+"),
                                         right = FALSE)
  
  plot_data$lib_conc_ngul_group <- cut(plot_data$lib_conc_ngul,
                                         breaks = c(-Inf,1.632,5,9,13,17,21,Inf),
                                         lables = c("1.632-","1.632-5","5-9","9-13","13-17","17-21","21+"),
                                         right = FALSE)
  plot_data$trimmomatic_group <- cut(plot_data$trimmomatic,
                                         breaks = c(0,0.8909,0.93,0.96,0.99,1),
                                         lables = c("0-0.8909","0.8909-0.93","0.93-0.96","0.96-0.99","0.99-1"),
                                         right = FALSE)
  plot_data$duplicate_ratio_group <- cut(plot_data$duplicate_ratio,
                                         breaks = c(0,0.2,0.4,0.6,0.8,1),
                                         lables = c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1"),
                                         right = FALSE)
  plot_data$genebody_chr_group <- cut(plot_data$genebody_chr,
                                         breaks = c(-Inf,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,Inf),
                                         lables = c("0.52-","0.52-0.53","0.53-0.54","0.54-0.55","0.55-0.56","0.56-0.59","0.59+"),
                                         right = FALSE)
  plot_data$DNAconc_ngul_group <- cut(plot_data$DNAconc_ngul,
                                      breaks = c(0,1,2,3,4,5,6,7,Inf),
                                      lables = c("0-1","1-2","2-3","3-4","4-5","5-6","6-7","7+"),
                                      right = FALSE)
  plot_data$fastq_size_group <- cut(plot_data$fastq_size,
                                      breaks = c(0,300,600,900,1200,1500,1800,2100,Inf),
                                      lables = c("0-300","300-600","600-900","900-1200","1200-1500","1500-1800","1800-2100","2100+"),
                                      right = FALSE)
  
  plot_data$lib_build_kits_group <- plot_data$lib_build_kits
  return(plot_data)
}


#####################################################################################
# 将input_sql初始化，赋值type与detail（未使用）
#####################################################################################
get_init_input <- function(input_sql, disease_type = NULL ) {
  # Attach columns of seq_id, type, feature_count_location and set rownames as seq_id
  # Arg:  
  #   input_sql: sql data pulled from sql
  #   disease_type： Health or Cancer or any symbol to indicate the type of the disease
  # Returns: 
  #   dataframe with new colums added
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
  input_sql$seq_id <- paste0("SEQ", input_sql$seq_num)
  rownames(input_sql) <- input_sql$seq_id
  input_sql$feature_count_location <- paste0(input_sql$featurecount_folder,"/", "SEQ",input_sql$seq_num,".genebody")
  output_col <- c( "seq_id" , "age" , "type" ,"detail" , "feature_count_location" , "genebody_chr", "duplicate_ratio" , "lib_conc_ngul" )
  input <- input_sql[ , output_col ]
  return( input )
}

#####################################################################################
# 通过inputs_df获取featurecount
#####################################################################################
get_featurecount <- function( input_list, unique_id=NULL, workers_number = NULL ){
  # Pharse the featurecount info from file , usually used to get correlation of samples
  # Args:
  #   input_list contains "seq_id" "feature_count_location"
  # Returns:
  #   sample rows & feature count column
  get_current_featurecount <- function( input ){
    return( read.table( as.character(input$feature_count_location) , head=TRUE , sep = "\t" , skip = 1)[,7] )
  }
  
  if ( is.null(workers_number)){
    cl <- makeCluster(25)
  } else{
    cl <- makeCluster(workers_number)
  }
  
  registerDoParallel(cl)
  fc<-foreach(x=1:dim(input_list)[1], .combine="cbind") %dopar% get_current_featurecount(input_list[x,])
  stopCluster(cl)
  
  gene_id <- as.character(read.table(as.character(input_list$feature_count_location[1]) , head=TRUE , sep = "\t" , skip = 1 )[,1])
  rownames(fc) <- gene_id
  if(!is.null(unique_id)){
    colnames(fc) <- input_list[,unique_id]
  } else {
    colnames(fc) <- input_list$seq_id
  }
  return( data.frame(fc , check.names = FALSE) )
}

#####################################################################################
# 获取病理信息（未使用）
#####################################################################################
get_pathology <- function( ... ){
  input_list <- list(...)
  columns <- c("seq_num" ,"diagnosis" ,"spec_diag1" ,"spec_diag2", "pathology" ,  "metastasis_site", "tumor_num", "max_tumor_size")
  for ( i in length( input_list) ){
    if (!exists("input")){
      input <- input_list[[i]][ , columns]
    } else {
      input <- rbind( input , input_list[[i]][ , columns] )
    }
  }
  return( input )
}

#####################################################################################
# 常规差异分析时用的input function
#####################################################################################
get_easy_input <- function(slight_input,serious_input,slight_file=NULL,serious_file=NULL,type_list=c("slight","serious"),num,file_name="input_list.rda"){
  # 获取加type后的两组样本总列表，注意将样本与type正确对应
  # Arg:
  #   slight_input/serious_input:两组样本，需列数、列名一致
  #   slight_file/serious_file:两组样本的seq_num文件，txt文件即可，注意要一一对应，默认为NULL
  #   type_list：需与get_result_list()的model_type_set一致，需视情况修改默认值
  #   file_name：保存文件名
  # Return:
  #   input_list:总样本列表
  if(!is.null(slight_file)){
    slight_data <- read.table(slight_file,header = FALSE,sep = '\n')[,1]
    slight_input <- slight_input[slight_input$seq_num %in% slight_data,]
  }
  if(!is.null(serious_file)){
    serious_data <- read.table(serious_file,header = FALSE,sep = '\n')[,1]
    serious_input <- serious_input[serious_input$seq_num %in% serious_data,]
  }
  slight_input$type <- type_list[1]
  serious_input$type <- type_list[2]
  input_list <- rbind(slight_input,serious_input)
  input_list$type <- factor(input_list$type,levels = type_list)
  input_list$seq_id <- paste0("SEQ",input_list$seq_num)
  #input_list$feature_count_location <- paste0(input_list$featurecount_folder,"_q/",input_list$seq_id,".genebody")
  #input_list$feature_count_location <- paste0(input_list$featurecount_folder,"/",input_list$seq_id,".genebody")
  input_list$feature_count_location <- paste0("/data_bak/rawdata/featurecounts/",input_list$seq_id,"/",input_list$seq_id,".genebody")
  rownames(input_list) <-input_list$seq_id
  save(input_list,file = file_name)
  return(input_list)
}

#####################################################################################
# 建模差异分析时用的input function
#####################################################################################
get_input_list <- function(input_type1,input_type2,seq_valid1=NULL,seq_valid2=NULL,seq_file=NULL,HEA_num,CAC_num,seed_num=1,
                           type_list=c("HEA","CAC"),file_name="input_list.rda"){
  # 获取分好T、V组后的原输入信息
  # Arg:  
  #   HEA/input：从数据库获取的健康/患病样本，需列数、列名一致
  #   seq_file：本次科研服务提供的患病样本的seq_num文件，一个txt文件即可
  #   HEA_num/CAC_num：训练集中健康/患病样本数设置
  #   seed_num：种子，修改可进行不同的T组采样
  #   stage：默认为NULL，返回分好组后的原输入信息；传入患病样本分期信息时，函数返回分好组后的分期信息
  # Returns: 
  #   分好组的input_list
  if (!is.null(seq_file)){
    seq_data <- read.table(seq_file,header = TRUE,sep = '\n')[,1]
    CAC <- input_type2[input_type2$seq_num %in% seq_data,]
  } else {
    CAC <- input_type2
  }
  HEA <- input_type1
  HEA$type <- type_list[1]
  CAC$type <- type_list[2]
  set.seed(seed_num)
  train_HEA <- sample(1:nrow(HEA),HEA_num)
  HEA_T <- HEA[train_HEA,]
  train_CAC <- sample(1:nrow(CAC),CAC_num)
  CAC_T <- CAC[train_CAC,]
  origin_T <- rbind(HEA_T,CAC_T)
  origin_T$grouping <- "T1V1"
  
  if (!is.null(seq_valid1) & !is.null(seq_valid2)){
    seq_valid1$type <- type_list[1]
    seq_valid2$type <- type_list[2]
    origin_V <- rbind(seq_valid1,seq_valid2)
  }
  
  if (!is.null(seq_valid1) & is.null(seq_valid2)){
    seq_valid1$type <- type_list[1]
    origin_V <- seq_valid1
  }
  
  if (is.null(seq_valid1) & !is.null(seq_valid2)){
    seq_valid2$type <- type_list[2]
    origin_V <- seq_valid2
  }
  
  if (is.null(seq_valid1) & is.null(seq_valid2)){
    origin_V <- rbind(HEA[-train_HEA,],CAC[-train_CAC,])
  }

  origin_V$grouping <- "V2"
  input_list <- rbind(origin_T,origin_V)
  input_list$type <- factor(input_list$type,levels = type_list)
  input_list$seq_id <- paste0("SEQ",input_list$seq_num)
  #input_list$feature_count_location <- paste0(input_list$featurecount_folder,"/",input_list$seq_id,".genebody")
  input_list$feature_count_location <- paste0("/data_bak/rawdata/featurecounts/",input_list$seq_id,"/",input_list$seq_id,".genebody")
  rownames(input_list) <-input_list$seq_id
  save(input_list,file = file_name)
  return( input_list )
}

#####################################################################################
# 通用函数，获取result和gene_list的function。执行较花时间，请耐心等待。
#####################################################################################
get_result_list <- function(T_list,file_name,model_type_set=c("HEA","CAC"),unique_id=NULL){
  # Arg:
  #   T_list:常规分析中的总样本/建模分析中的训练集样本！！！
  #   file_name：保存文件名
  #   model_type_set：与type_list对应，轻症/健康在前标记为0，重症/患病在后标记为1，切勿混淆
  # Output:
  #   file_name:包含各项结果的rda文件，执行完成后需在主程序执行load(file_name)将结果加载到工作空间 
  if(!is.null(unique_id)){
    colData <- data.frame("type"=T_list$type,
                          row.names = T_list[,unique_id],check.names = FALSE)
  } else {
    colData <- data.frame("type"=T_list$type,
                          row.names = rownames(T_list),check.names = FALSE)
  }
  
  colData$type <- factor(colData$type,levels = model_type_set)
  countData <- get_featurecount(T_list,unique_id)
  # countData去除ASB4基因
  #countData <- countData[!(rownames(countData) %in% "ASB4"),]
  countData <- countData[rowSums(countData)>10,]
  geoMeans <- apply(countData,1,function(x) exp(mean(log(x))))
  dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~type)
  dds <- DESeq(dds,parallel = TRUE)
  res <- results(dds)
  
  vst <- varianceStabilizingTransformation(dds,blind = FALSE)
  vst_all <- assay(vst)
  save(countData,colData,dds,res,vst_all,vst,geoMeans,file = file_name)
}

#####################################################################################
# 根据条件筛选gene_list
#####################################################################################
get_gene_list <- function(res,cutoff_list=c(0.05,0.05,0.26),flag_list=c(TRUE,TRUE),
                          up_down=c(100,0),res_file="resOrdered.csv",gene_file="gene_list.rda"){
  # 获取最终选取的top差异基因。函数可打印出可选的上升/下降基因数，可据此调整参数值
  # 若条件筛选后没有可供选择的基因会抛出错误提示，此时需合理调整cutoff_list或flag_list
  # Arg:
  #   cutoff_list:三个值分别代表pvalue阈值，padj阈值，log2FoldChange阈值
  #   flag_list：两个布尔型值分别代表是否启用padj阈值筛选，是否启用log2FoldChange阈值
  #   up_down：指定选取的上升、下降基因数。如c(60,40)代表从上升基因取60，下降基因取40
  # Output:
  #   resOrdered.csv，gene_list.rda  
  # Retuen:
  #   gene_list:选取的差异基因
  resOrdered <- res[which(res$pvalue<cutoff_list[1]),]
  #resOrdered <- res[res$pvalue < cutoff_list[1],]
  resOrdered <- resOrdered[!is.na(resOrdered$log2FoldChange),]
  if (flag_list[1] == TRUE){
    resOrdered <- resOrdered[!is.na(resOrdered$padj),]
    resOrdered <- resOrdered[resOrdered$padj<cutoff_list[2],]
    resOrdered <- resOrdered[order(resOrdered$padj),]
  } else {
    resOrdered <- resOrdered[order(resOrdered$pvalue),]
  }
  
  if (flag_list[2] == TRUE){
    resOrdered <- resOrdered[abs(resOrdered$log2FoldChange)>cutoff_list[3],]
  }
  resOrdered_high <- resOrdered[resOrdered$log2FoldChange>0,]
  resOrdered_low <- resOrdered[resOrdered$log2FoldChange<0,]
  print(paste("可选上升基因数：",length(rownames(resOrdered_high))))
  print(paste("可选下降基因数：",length(rownames(resOrdered_low))))
  high_num <- length(rownames(resOrdered_high))
  low_num <- length(rownames(resOrdered_low))
  #save(high_num,low_num,file = "num.rda")
  if(length(rownames(resOrdered_high))<up_down[1] | length(rownames(resOrdered_low))<up_down[2]){
    print("可选的上升或下降基因数不能满足up_down参数的指定，请知悉")
  }
  gene_list_high <- matrix(rownames(resOrdered_high)[0:up_down[1]])
  gene_list_low <- matrix(rownames(resOrdered_low)[0:up_down[2]])
  gene_list <- rbind(gene_list_high,gene_list_low)
  if (length(gene_list) == 0){
    stop("No genes to choose from.")
  }else {
    write.csv(resOrdered,file = res_file)
    save(gene_list,high_num,low_num,file = gene_file)
    return( gene_list )}
}

#####################################################################################
# 获取与gene_list相匹配的vst
#####################################################################################
get_vst <- function(vst,gene_list,geoMeans){
  vst <- assay(vst)[match(gene_list,names(geoMeans)),]
  if(length(data.frame(vst))==1){
    vst <- data.frame(t(vst))
    rownames(vst) <- gene_list
  }
  return(vst)
}

#####################################################################################
# 建模差异分析时用的建模function
#####################################################################################
get_cvfit <- function(vst,colData,file_name="cvfit_model.rda"){
  # 获取模型，打印模型的lambda.1se参数，该参数可供参考选较佳模型
  # Output:
  #   file_name：保存模型，需要时可load(file_name)加载至工作空间
  cvfit <- cv.glmnet(x=t(vst),y=colData$type,family="binomial",type.measure="class",alpha=1)
  print(paste("lambda.1se:",round(cvfit$lambda.1se,7)))
  save(cvfit,file = file_name)
  return(cvfit)
}

#####################################################################################
# 其他功能函数，取rpm、rpkm的函数（不在分析流程中，供参考）
#####################################################################################
get_rpm <- function( input_list, target_gene_list = NULL, merge_flag = FALSE, workers_number = NULL){
  # Computes the RPM value by feature count file specified in the input_list.
  # Args:
  #   input_list: must contains "seq_id", "feature_count_location".
  #   merge_flag: merge rpm matrix and input_list or not(by seq_id)
  # Returns:
  #    rpm df
  print(" please make sure your are using rpm. ensure you need genecounts/allcount ")
  get_current_rpm <- function(input){
    current_count_data <- read.table( as.character(input$feature_count_location) , head=TRUE , sep = "\t" , skip = 1)[,7]
    countData_sum <- sum(current_count_data)
    rpm <- current_count_data / countData_sum
    return(rpm)
  }
  
  if ( is.null(workers_number) ){
    cl <- makeCluster(10)
  } else{
    cl <- makeCluster(workers_number)
  }
  
  registerDoParallel(cl)
  rpm<-foreach(x=1:dim(input_list)[1], .combine="cbind") %dopar% get_current_rpm(input_list[x,])
  stopCluster(cl)
  gene_id <- as.character(read.table(as.character(input_list$feature_count_location[1]) , head=TRUE , sep = "\t" , skip = 1 )[,1])
  rownames(rpm) <- gene_id
  colnames(rpm) <- input_list$seq_id
  # get predictor and corresponding rpm value
  if ( is.null(target_gene_list) ){
    genes_rpm <- rpm
  } else if ( is.character(target_gene_list) ){
    genes_rpm <- rpm[ target_gene_list ,]
  } else {
    quit()
  }
  genes_rpm <- data.frame(t( genes_rpm ))
  if(merge_flag){
    genes_rpm$seq_id <- rownames(genes_rpm)
    genes_rpm <- inner_join(input_list,genes_rpm)
    rownames(genes_rpm) <- genes_rpm$seq_id
  }
  return(genes_rpm)
}

get_rpkm <- function( input_list, target_gene_list = NULL, merge_flag = FALSE, workers_number = NULL){
  print(" please make sure your are using rpkm.ensure you need genecounts/(allcount * gene_length) ")
  load( '/work_data/genes_length.rda' )
  get_current_rpm <- function(input){
    current_count_data <- read.table( as.character(input$feature_count_location) , head=TRUE , sep = "\t" , skip = 1)[,7]
    countData_sum <- sum(current_count_data)
    rpm <- current_count_data / countData_sum
    return(rpm)
  }
  
  if ( is.null(workers_number) ){
    cl <- makeCluster(10)
  } else{
    cl <- makeCluster(workers_number)
  }
  
  registerDoParallel(cl)
  rpm<-foreach(x=1:dim(input_list)[1], .combine="cbind") %dopar% get_current_rpm(input_list[x,])
  stopCluster(cl)
  gene_id <- as.character(read.table(as.character(input_list$feature_count_location[1]) , head=TRUE , sep = "\t" , skip = 1 )[,1])
  rownames(rpm) <- gene_id
  colnames(rpm) <- input_list$seq_id
  # get predictor and corresponding rpm value
  input_rpm <- t(rpm)
  input_rpkm <- data.frame(sapply( 1:19100 , function(x)input_rpm[,x]*10^6/genes_length[x]))
  rownames(input_rpkm) <- input_list$seq_id
  colnames(input_rpkm) <- gene_id
  if ( is.null(target_gene_list) ){
    genes_rpkm <- input_rpkm
  } else if ( is.character(target_gene_list) ){
    genes_rpkm <- input_rpkm[ ,target_gene_list]
  } else {
    quit()
  }
  if(merge_flag){
    genes_rpkm$seq_id <- rownames(genes_rpkm)
    genes_rpkm <- inner_join(input_list,genes_rpkm)
    rownames(genes_rpkm) <- genes_rpkm$seq_id
  }
  return(genes_rpkm)
}

