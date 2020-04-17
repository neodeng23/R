library(DESeq2)
library(pROC)
library(ggplot2)
library(pheatmap)
library(glmnet)
library(annotate)
library(clusterProfiler)
library(org.Hs.eg.db)

#####################################################################################
# 流程（未使用）
#####################################################################################
gene_profile_deseq <- function( countData , input , model_type_set , case_num = "default" , volcano_cutoff = 0.3 , sort_method = "log2FoldChange" ){
  # get all results
  # Arg:
  #   countData: raw counts table (no normalization)
  #   input: 输入的原数据，从sql中提取
  #   model_type_set: 标签对应的顺序。c( "HEA" , "CAC" )
  #   case_num: 不同case可以输出不同名字
  #   volcano_cutoff = 0.3 火山图cutoff
  #   sort_method = "log2FoldChange" 通过log2FoldChange排序resOrdered
  # output: 
  #   Top 100 most differentially expressed genes list 
  #   PCA plot
  #   Volcano plot
  #   Heatmap for DE genes with hierarchical clustering
  #   AUC plot
  #   GO analysis plot regarding BP ontology
  #   GO analysis plot with "MF", "BP", and "CC" subontologies
  #   KEGG analysis result
  
  # colData: columns correspond to type condition, and rows correspond to the samples
  colData <-  data.frame("type" = input[,c("type")] , row.names = rownames(input) , check.names = FALSE )
  colData$type <- factor(colData$type , levels = model_type_set)
  # filter low counts
  countData <- countData[ rowSums(countData) > 10 , ]
  geoMeans <- apply( countData ,1,function(x) exp(mean(log(x))))
  countData <- countData[ names(geoMeans),  ]
  # DESeq分析
  # 返回genes pvalue, padj, log2FoldChange
  dds <- DESeqDataSetFromMatrix( countData = countData , colData = colData, design = ~type )
  dds <- DESeq(dds,parallel=TRUE)
  res <- results(dds)
  # 选择显著基因，过滤na log2FoldChange
  resOrdered <- res[which(res$pvalue<0.05),]
  resOrdered <- resOrdered[ !is.na(resOrdered$log2FoldChange),]
  # 通过指定的方法排序
  if ( sort_method == "log2FoldChange"){
    resOrdered <- resOrdered[order(resOrdered$log2FoldChange , decreasing = TRUE), ]
    # 选取+50 和 -50个基因
    gene_list <- rownames( resOrdered[c(1:50,(dim(resOrdered)[1]-49):dim(resOrdered)[1]) , ] )
    write( gene_list , file = paste0( "gene_list", case_num , ".csv" ) )
  } else if (sort_method == "padj" ){
    resOrdered <- resOrdered[order(resOrdered$padj), ]
    gene_list <- c( rownames( resOrdered[ resOrdered$log2FoldChange > 0 , ] )[1:50] , rownames( resOrdered[ resOrdered$log2FoldChange < 0 , ] )[1:50])
    write( gene_list , file = paste0( "gene_list", case_num , ".csv" ) )
  }
  # 保存数据
  save( countData , dds , res , resOrdered, input ,geoMeans , file = paste0("deseq_res",case_num,".rda") )
  write.table(resOrdered, file = paste0( "resOrdered" , case_num , ".csv" ) , sep = ",")
  # 做 vst
  vst <-  varianceStabilizingTransformation(dds , blind = FALSE)
  vst_all <- assay(vst)
  vst <- assay(vst)[match(gene_list, names(geoMeans)),]
  save( gene_list ,vst_all , vst , file = "gene_result.rda" )
  # 火山图
  volcano_plot(res , paste0("volcano_plot",case_num,".png") , volcano_cutoff )
  # vst <- vst[,order(fc$type)]
  # 热图
  plot_heatmap( vst  , gene_list , input , paste0( "topheatmap",case_num , ".pdf" ))
  # 全部基因pca和差异基因pca
  pca_plot( vst_all , input , paste0("PCA",case_num,".pdf") , file_name_pr = paste0("PCApr",case_num,".csv") )
  pca_plot_deg( vst , gene_list ,  input , paste0("PCA_deg_",case_num,".pdf") , file_name_pr = paste0("PCApr_deg_",case_num,".csv") )
  # cvfit建模
  cvfit <- cv.glmnet(x = t(vst), y = colData$type ,family = "binomial", type.measure = "class")
  save( cvfit , file = "cvfit_model.rda" )
  # GO分析
  GOanaly_plot(res , filename = paste0("GO",case_num,".pdf"), data_filename = paste0("GO",case_num,".txt"))
  # KEGG分析
  keggplot( res , gene_list )
  return( list( "resOrdered"=resOrdered , "cvfit"= cvfit , "res" = res , "gene_list" = gene_list , "geoMeans" = geoMeans , "dds" = dds , "vst" = vst) )
}

#####################################################################################
# 常规差异分析pipeline
#####################################################################################
routine_analysis <- function(user,password,sql_statement1,sql_statement2,type_list,sql_file1=NULL,sql_file2=NULL,unique_id=NULL,
                            case_num = "default",cutoff_list=c(0.05,0.05,0.26),flag_list=c(TRUE,TRUE),up_down=c(100,0),
                            volcano_cutoff = 0.26,GO_flag="all",pathway_num=20,pvalueCut=0.5,qvalueCut=0.5,cor_cutoff=0.95){
  # get all results
  # Arg:
  #   user/password: 数据库用户名/密码
  #   sql_statement1/sql_statement2: 数据库提取两类数据，sql_statement1为轻症/健康，sql_statement2为重症/患病
  #   sql_file1/sql_file2：两类数据的seq_num文件，txt文件即可，注意要一一对应。如不需此项数据，置为NULL即可。
  #   type_list：两类数据的类型标记，轻症/健康在前标记为0，重症/患病在后标记为1
  #   case_num: 不同case可以起不同名字作为标识
  #   cutoff_list:三个值分别代表pvalue阈值，padj阈值，log2FoldChange阈值
  #   flag_list：两个布尔型值分别代表是否启用padj过滤，是否启用log2FoldChange过滤
  #   up_down：指定选取的上升、下降基因数。如c(60,40)代表从上升基因取60，下降基因取40
  #   volcano_cutoff：火山图的cutoff
  #   GO_flag:只可以为"all"、"bp"、"cc"、"mf"
  #   pathway_num/pvalueCut/qvalueCut：控制kegg气泡图里数据的参数
  #   cor_cutoff：根据featurecount做样本间correlation图的cutoff
  # output: 
  #   Top most differentially expressed genes list
  #   Volcano plot、Heatmap plot、PCA plot
  #   GO analysis plot、KEGG plot
  #   Correlation plot、Cluster plot
  input_type1 <- get_sqloutput(user,password,sql_statement1)
  input_type2 <- get_sqloutput(user,password,sql_statement2)
  input_list <- get_easy_input(input_type1,input_type2,slight_file=sql_file1,serious_file=sql_file2,type_list,case_num,file_name=paste0("input_list_",case_num,".rda"))
  
  ######  其他分布图形
  if("age" %in% colnames(input_list)){
    if (length(na.omit(input_list[,c("age")]))== dim(input_list)[1]){
      plot_age_distr(input_list,type_list,paste0("age_",case_num,".pdf"))
    }
    if (length(na.omit(input_list[,c("age")]))!= dim(input_list)[1]){
      omitna_data <- input_list[complete.cases(input_list[,"age"]),]
      plot_age_distr(omitna_data,type_list,paste0("age_subset_",case_num,".pdf"))
    }
  }
  
  if("gender" %in% colnames(input_list)){
    if (length(na.omit(input_list[,c("gender")]))==dim(input_list)[1]){
      plot_gender_distr(input_list,type_list,paste0("gender_",case_num,".pdf"))
    }
    if (length(na.omit(input_list[,c("gender")]))!= dim(input_list)[1]){
      omitna_data <- input_list[complete.cases(input_list[,"gender"]),]
      plot_gender_distr(omitna_data,type_list,paste0("gender_subset_",case_num,".pdf"))
    }
  }
  
  if("stage" %in% colnames(input_list)){
    stage_type2 <- input_list[input_list$type %in% type_list[2],]
    plot_stage_distr(stage_type2,paste0("stage_",case_num,".png"))
  }

  ######  差异基因分析
  get_result_list(input_list,paste0("deseq_res_",case_num,".rda"),type_list,unique_id)
  load(paste0("deseq_res_",case_num,".rda"))
  gene_list <- get_gene_list(res,cutoff_list,flag_list,up_down,res_file=paste0("resOrdered_",case_num,".csv"),gene_file=paste0("gene_list_",case_num,".rda"))
  vst <- get_vst(vst,gene_list,geoMeans)  
  ######  重要图形输出
  volcano_plot(res,paste0("volcano_",case_num,".png"),volcano_cutoff)
  plot_heatmap(vst,gene_list,input_list,unique_id,paste0("topheatmap_",case_num,".pdf"))
  pca_plot(vst_all,input_list,unique_id,paste0("PCA_",case_num,".png"),file_name_pr = paste0("PCApr_",case_num,".csv"))
  pca_plot_deg(vst,gene_list,input_list,unique_id,paste0("PCA_deg_",case_num,".png"),file_name_pr = paste0("PCApr_deg_",case_num,".csv"))
  ######  其他图形输出
  plot_sample_cor(input_list,countData,paste0("correlation_",case_num,".pdf"),cor_cutoff)
  if(nrow(input_list) <=50){
    plot_clust(vst,colData,paste0("cluster_",case_num,".pdf"))
    }
  ######  生信分析图形
  GOanalysis_plot(res,paste0("GO_",case_num,".pdf"),flag=GO_flag,data_file=paste0("GO_detail_",case_num,".csv"))
  kegg_bubble_plot(res, pathway_num, paste0("kegg_plot_",case_num,".pdf"),pvalueCut,qvalueCut)
  keggplot(res,gene_list)
}

#####################################################################################
# 建模差异分析pipeline
#####################################################################################
model_analysis <- function(user,password,sql_statement1,sql_statement2,sql_valid1=NULL,sql_valid2=NULL,seq_file,HEA_num,CAC_num,seed_num,type_list,
                           unique_id=NULL,case_num = "default",cutoff_list=c(0.05,0.05,0.26),flag_list=c(TRUE,TRUE),
                           up_down=c(100,0),volcano_cutoff = 0.26,GO_flag="all",pathway_num=20,pvalueCut=0.5,qvalueCut=0.5,cor_cutoff=0.97){
  # get all results
  # Arg:
  #   user/password: 数据库用户名/密码
  #   sql_statement1/sql_statement2: 数据库提取两类数据，sql_statement1为轻症/健康，sql_statement2为重症/患病
  #   sql_valid1/sql_valid2: 指定验证组中的轻症/健康和重症/患病。有这两个参数中的任一参数则采样剩余的组不能作为验证集中的数据。
  #   seq_file：本次科研服务提供的患病样本的seq_num文件，一个txt文件即可
  #   HEA_num/CAC_num/seed_num：训练集中健康/患病样本数设置/种子，修改可进行不同的T组采样
  #   type_list：两类数据的类型标记，轻症/健康在前标记为0，重症/患病在后标记为1
  #   case_num: 不同case可以起不同名字作为标识
  #   cutoff_list:三个值分别代表pvalue阈值，padj阈值，log2FoldChange阈值
  #   flag_list：两个布尔型值分别代表是否启用padj过滤，是否启用log2FoldChange过滤
  #   up_down：指定选取的上升、下降基因数。如c(60,40)代表从上升基因取60，下降基因取40
  #   volcano_cutoff：火山图的cutoff
  #   GO_flag:只可以为"all"、"bp"、"cc"、"mf"
  #   pathway_num/pvalueCut/qvalueCut：控制kegg气泡图里数据的参数
  #   cor_cutoff：根据featurecount做样本间correlation图的cutoff
  # output: 
  #   Top most differentially expressed genes list
  #   Volcano plot、Heatmap plot、PCA plot
  #   AUC plot、Logit plot
  #   GO analysis plot、KEGG plot
  #   Correlation plot、Cluster plot
  #   Age/Gender/Stage distribution
  # 数据库取数
  input_type1 <- get_sqloutput(user,password,sql_statement1)
  input_type2 <- get_sqloutput(user,password,sql_statement2)
  
  if (!is.null(sql_valid1)){
    seq_valid1 <- get_sqloutput(user,password,sql_valid1)
  } else {
    seq_valid1 <- NULL
  }
  
  if (!is.null(sql_valid2)){
    seq_valid2 <- get_sqloutput(user,password,sql_valid2)
  } else {
    seq_valid2 <- NULL
  }

  # 采样并整形获取input_list
  input_list <- get_input_list(input_type1,input_type2,seq_valid1,seq_valid2,seq_file,HEA_num,CAC_num,seed_num,
                                 type_list,file_name=paste0("input_list_",case_num,".rda"))
  
  ######  其他分布图形
  origin_copy <- input_list
  origin_copy$mold <- paste0(origin_copy$type,origin_copy$grouping)
  
  if("age" %in% colnames(input_list)){
    if (length(na.omit(input_list[,c("age")]))== dim(input_list)[1] & length(levels(as.factor(origin_copy$mold)))==4){
      plot_age_distr(origin_copy,type_list,paste0("age_",case_num,".pdf"))
    }
    if (length(na.omit(input_list[,c("age")]))!= dim(input_list)[1] & length(levels(as.factor(origin_copy$mold)))==4){
      omitna_data <- origin_copy[complete.cases(origin_copy[,"age"]),]
      plot_age_distr(omitna_data,type_list,paste0("age_subset_",case_num,".pdf"))
    }
  }

  if("gender" %in% colnames(input_list)){
    if (length(na.omit(input_list[,c("gender")]))==dim(input_list)[1] & length(levels(as.factor(origin_copy$mold)))==4){
      plot_gender_distr(origin_copy,type_list,paste0("gender_",case_num,".pdf"))
    }
    if (length(na.omit(input_list[,c("gender")]))!= dim(input_list)[1] & length(levels(as.factor(origin_copy$mold)))==4){
      omitna_data <- origin_copy[complete.cases(origin_copy[,"gender"]),]
      plot_gender_distr(omitna_data,type_list,paste0("gender_subset_",case_num,".pdf"))
    }
  }
  
  if ("stage" %in% colnames(input_list)){
    CAC_stage_T <- input_list[input_list$type == type_list[2] & input_list$grouping == "T1V1",]
    CAC_stage_T$type <- paste0(type_list[2],"_training")
    CAC_stage_V <- input_list[input_list$type == type_list[2] & input_list$grouping == "V2",]
    CAC_stage_V$type <- paste0(type_list[2],"_validation")
    origin_stage <- rbind(CAC_stage_T,CAC_stage_V)
    plot_stage_distr(origin_stage,paste0("stage_",case_num,".png"))
  }

  ######  对训练集进行差异基因分析
  T_group <- input_list[input_list$grouping=="T1V1",]
  V_group <- input_list[input_list$grouping=="V2",]
  get_result_list(T_group,paste0("deseq_res_",case_num,".rda"),type_list,unique_id)
  load(paste0("deseq_res_",case_num,".rda"))
  gene_list <- get_gene_list(res,cutoff_list,flag_list,up_down,res_file=paste0("resOrdered_",case_num,".csv"),gene_file=paste0("gene_list_",case_num,".rda"))
  vst <- get_vst(vst,gene_list,geoMeans)
  ######  重要图形输出
  volcano_plot(res,paste0("volcano_",case_num,".png"),volcano_cutoff)
  plot_heatmap(vst,gene_list,input_list,unique_id,paste0("topheatmap_",case_num,".pdf"))
  pca_plot(vst_all,T_group,unique_id,paste0("PCA_",case_num,".png"),file_name_pr = paste0("PCApr_",case_num,".csv"))
  if(dim(vst)[1]!=1){
    pca_plot_deg(vst,gene_list,T_group,unique_id,paste0("PCA_deg_",case_num,".png"),file_name_pr = paste0("PCApr_deg_",case_num,".csv"))
    }
  ######  其他图形输出
  plot_sample_cor(T_group,countData,paste0("correlation_",case_num,".pdf"),cor_cutoff)
  if(nrow(input_list) <=50){
    plot_clust(vst,colData,paste0("cluster_",case_num,".pdf"))
    }
  ########  建模并绘制判分图形
  cvfit <- get_cvfit(vst,colData,paste0("cvfit_",case_num,".rda"))
  load(paste0("cvfit_",case_num,".rda"))
  T_score <- get_pred_score(cvfit, dds, geoMeans, gene_list, input = T_group)
  V_score <- get_pred_score(cvfit, dds, geoMeans, gene_list, input = V_group)
  TV_score <- rbind(T_score,V_score)
  TV_score$seq_id <- rownames(TV_score)
  score <- inner_join(input_list,TV_score[,c("seq_id","score")])
  write.csv(score,file = paste0("score_",case_num,".csv"),row.names = FALSE)
  cutoff <- coords(roc(V_score$type,V_score$score,levels=type_list), "b",transpose = FALSE)$threshold[1]
  print(paste("Youden index threshold(V group)：",cutoff))
  print(paste(c("T sensitivity：","T specificity：","T accuracy："),coords(roc(T_score$type,T_score$score,levels=type_list),cutoff, ret=c("sen","spec","accu"),transpose = TRUE)))
  print(paste(c("V sensitivity：","V specificity：","V accuracy："),coords(roc(V_score$type,V_score$score,levels=type_list),cutoff, ret=c("sen","spec","accu"),transpose = TRUE)))
  if(length(levels(V_score$type))==2){
    auc_plot(T_score,V_score,data_name = c("Training","Validation"),type_list=type_list,file_name=paste0("ROC_",case_num,".pdf"))
  } else {
    auc_plot(T_score,data_name = c("Training"),type_list=type_list,file_name=paste0("ROC_",case_num,".pdf"))
  }
  logit_plot(T_score,"T Group Logit Plot",paste0("T_logit_",case_num,".png"))
  logit_plot(V_score,"V Group Logit Plot",paste0("V_logit_",case_num,".png"))
  ######  生信分析图形
  GOanalysis_plot(res,paste0("GO_",case_num,".pdf"),flag=GO_flag,data_file=paste0("GO_detail_",case_num,".csv"))
  kegg_bubble_plot(res, pathway_num, paste0("kegg_plot_",case_num,".pdf"),pvalueCut,qvalueCut)
  keggplot(res,gene_list)
}
