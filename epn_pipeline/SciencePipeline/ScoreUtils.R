library(doParallel)
library(parallel)
library(DESeq2)

#####################################################################################
# 留一法建模及判分
#####################################################################################
loo <- function(input, countData, model_type_set){
  fc <- data.frame(t(countData))
  fc$type <- input$type
  cl <- makeCluster( 30L )
  registerDoParallel(cl)
  score_df <- foreach( i=1:dim(fc)[1] , .combine = "rbind" , .packages = c("glmnet", "DESeq2" )) %dopar% {
    new_fc <- fc[-i,]
    countData <- t(new_fc[,1:(dim(new_fc)[2]-1)])
    countData_test <- t( fc[ i,1:(dim(fc)[2]-1)] )
    colData <- data.frame("type" = new_fc[ ,c("type")], row.names = rownames(new_fc))
    
    dds <- DESeqDataSetFromMatrix(countData=countData , colData = colData , design = ~type )
    dds <- DESeq(dds)
    res <- results(dds)
    resOrdered <- res[order(res$log2FoldChange),]
    resOrdered <- resOrdered[ !is.na(resOrdered$log2FoldChange),]
    gene_list <- rownames( resOrdered[c(1:50,(dim(resOrdered)[1]-49):dim(resOrdered)[1]) , ] )
    geoMeans <- apply(countData,1,function(x) exp(mean(log(x))))
    
    vst <-  varianceStabilizingTransformation(dds , blind = FALSE)
    vst <- assay(vst)[match(gene_list, names(geoMeans)),]
    cvfit <- cv.glmnet(x = t(vst), y = colData$type ,family = "binomial", type.measure = "class")
    colnames(countData_test) <- "counts"
    test_df <- data.frame(countData_test,"geo"=geoMeans)
    test_df <- test_df[test_df$geo != 0,]
    test_df$div <- test_df$counts/test_df$geo
    test_dds <- DESeqDataSetFromMatrix(countData=countData_test, colData=data.frame(condition="unknown"), design=~1)
    #test_dds <- estimateSizeFactors(test_dds, geoMeans=geoMeans)
    sizeFactors(test_dds) <- median(test_df$div)
    dispersionFunction(test_dds) <- dispersionFunction(dds)
    test_vst <- varianceStabilizingTransformation(test_dds, blind=FALSE)
    test_vst <- assay(test_vst)[match(gene_list, names(geoMeans)),]
    score <- data.frame( predict( cvfit , newx = t(test_vst) , s="lambda.1se", type="response" ) , as.character( fc[i,"type"] ))
    rownames(score) <- input[,"seq_id"][i]
    colnames(score) <- c("score","type")
    return( score )
  }
  stopCluster(cl)
  score_df$type <- factor(score_df$type,levels = model_type_set)
  return( score_df )
}

#####################################################################################
# 判分函数
#####################################################################################
get_pred_score <- function( cvfit, dds, geoMeans, gene_list, input_list, workers_number = NULL, filename = NULL){
  # 通过cvfit和数据产生的dds,geoMeans,gene_list 对input组判分
  # Args:
  #   cvfit:  cv.glmnet产生的class
  #   dds:  Deseq输出
  #   geoMeans:   Pipeline输出
  #   gene_list:  Pipeline输出
  #   input:  需要进行判分的input
  # Output:
  #   filename:保存判分文件，一个csv文件即可
  # return:
  #   score_df
  pred <- function(input){
    counts <- read.delim(as.character(input$feature_count_location) , skip=1, row.names=1)
    #去除ASB4基因
    #counts <- counts[!(rownames(counts) %in% "ASB4"),]
    counts <- counts[names(geoMeans),6]
    test_df <- data.frame("counts"=counts,"geo"=geoMeans)
    test_df <- test_df[test_df$geo != 0,]
    test_df$div <- test_df$counts/test_df$geo
    #dds1 <- estimateSizeFactors(dds1, geoMeans=geoMeans) 
    dds1 <- DESeqDataSetFromMatrix(counts, colData=data.frame(condition="unknown"), design=~1) 
    sizeFactors(dds1) <- median(test_df$div)
    #mcols(dds1)$dispFit <- dispFit
    dispersionFunction(dds1)<-dispersionFunction(dds)
    rld2 <- varianceStabilizingTransformation(dds1, blind=F) 
    rld3 <- assay(rld2) [match(gene_list, names(geoMeans)),]
    sub_score <- data.frame("score"=predict( cvfit, newx = t(rld3), s="lambda.1se", type="response"), "type" = input$type)
    #sub_score <- data.frame(predict(cvfit, newx=t(rld3), s="lambda.1se", type="response"))
    #colnames(sub_score) <- input$seq_id
    rownames(sub_score) <- input$seq_id
    colnames(sub_score) <- c("score","type")
    return(sub_score)
  }
  
  if ( is.null(workers_number)){
    cl <- makeCluster(40)
  } else{
    cl <- makeCluster(workers_number)
  }

  registerDoParallel(cl)
  score_df<-foreach(x=1:dim(input_list)[1], .combine="rbind",.packages = c("DESeq2", "glmnet")) %dopar% pred(input_list[x,])
  stopCluster(cl)
  
  if (!is.null( filename )){
    write.csv(score_df, file = filename )
    #write.table( score_df , file = filename , sep = "," , col.names = TRUE)
  }
  return(score_df)
}

#####################################################################################
# estimateSizeFactors()函数对sizeFactor矫正了两次（弃用判分函数，未使用）
#####################################################################################
get_score <- function( cvfit , dds , geoMeans , gene_list , input , filename = NULL ){
  fc_test <- data.frame( get_featurecount( input ) , check.names = FALSE)
  countData_test <- fc_test[ names( geoMeans ) ,]
  
  test_dds <- DESeqDataSetFromMatrix(countData=countData_test, colData=data.frame( rep( "unknown" , dim(input)[1] )), design=~1) 
  test_dds <- estimateSizeFactors(test_dds, geoMeans=geoMeans)
  dispersionFunction(test_dds) <- dispersionFunction(dds)
  test_vst <- varianceStabilizingTransformation(test_dds, blind=FALSE)
  test_vst <- assay(test_vst)[match(gene_list, names(geoMeans)),]
  
  score_df <- data.frame("score"=predict( cvfit, newx = t(test_vst), s="lambda.1se", type="response"), "type" = input$type)
  colnames(score_df) <- c("score","type")
  if (!is.null( filename )){
    write.csv(score_df, file = filename )
    #write.table( score_df , file = filename , sep = "," , col.names = TRUE)
  }
  return( score_df )
}