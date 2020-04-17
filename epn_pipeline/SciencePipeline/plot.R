library(pROC)
library(ggplot2)
library(pheatmap)
library(annotate)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(plyr)

##########################################################################################
# 分别根据gc值、duplication、文库浓度绘制样本间的correlation图（不在分析流程中，可供参考）
##########################################################################################
plotCor <- function( input, filename_diag, qc_thres, detail_sort = FALSE ){
  input$genebody_chr <- cut(input$genebody_chr,
                            breaks=c(-Inf,0.56,0.586,0.59,0.6,Inf),
                            labels=c("0.56-","0.56-0.586","0.586-0.59","0.59-0.6","0.6+"),
                            right=FALSE)
  input$duplicate_ratio <- cut(input$duplicate_ratio,
                               breaks=c(-Inf,0.2,0.8,Inf),
                               labels=c("0.2-","0.2-0.8","0.8+"),
                               right=FALSE)
  input$lib_conc_ngul <- cut(input$lib_conc_ngul,
                             breaks=c(-Inf,1.632,13.56,15,Inf),
                             labels=c("1.632-","1.632-13.56","13.56-15","15+"),
                             right=FALSE)

  for ( annotate_type in c("genebody_chr", "duplicate_ratio" , "lib_conc_ngul") ){
    input <- input[ order( input[annotate_type] ) , ]
    if ( detail_sort ){
      input <- input[ order( input["detail"] ), ]
    }
    fc  <- get_featurecount(input)
    fc_cor <- cor(fc)
    for (i in 1:dim(fc_cor)[1]){
      fc_cor[i,i] <- 0
    }
    fc_cor[ fc_cor < qc_thres] <- qc_thres
    pheatmap( fc_cor ,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F, 
             annotation_names_row = FALSE , annotation_names_col = FALSE , 
             annotation_col = data.frame(row.names = input$seq_id ,input[,c("detail" , annotate_type)]),
             annotation_row = data.frame(row.names = input$seq_id ,input[,c("detail" , annotate_type)]) ,
             filename = paste0( filename_diag , "_" , annotate_type ,".png"),
             width = 10 , height= 7.69)
  }
}

##########################################################################################
# 绘制质控图（不在分析流程中，可根据实际需求参考并修改使用）
##########################################################################################
plot_picture <- function(plot_data,countData,file_path,cutoff_down=NULL,cutoff_up=NULL){
  var_first <- c("binding_rate","enrichment_rate","elution_rate","lib_conc_ngul","trimmomatic","duplicate_ratio","genebody_chr","DNAconc_ngul","fastq_size")
  var_second <- c(var_first,"lib_build_kits")
  cor_matrix <- cor(countData)
  if(!is.null(cutoff_down)){
    cor_matrix[cor_matrix<cutoff_down] <- cutoff_down
  }
  if(!is.null(cutoff_up)){
    cor_matrix[cor_matrix>=cutoff_up] <- cutoff_up
  }
  cor_matrix_median <- apply(cor_matrix,1,median)
  cor_matrix_median_info <- cbind(plot_data[,var_second],cor_matrix_median)
  write.table(cor_matrix_median_info,file=paste0(file_path,"/cor_median_",Sys.Date(),".txt"),quote = F,sep = "\t")
  
  for ( heatmap_var in var_second){
    plot_data <- plot_data[order(plot_data[,heatmap_var]),]
    this_var_group <- paste0(heatmap_var,"_group")
    pheatmap( cor_matrix[rownames(plot_data),rownames(plot_data)],cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F, 
              annotation_names_row = FALSE , annotation_names_col = FALSE ,
              annotation_col = data.frame(row.names = rownames(plot_data),var_group=plot_data[,this_var_group]),
              annotation_row = data.frame(row.names = rownames(plot_data),var_group=plot_data[,this_var_group]),
              main = paste0("Heatmap group by ",heatmap_var),
              filename = paste0(file_path,"/QC_orderby_",heatmap_var,".png"),fontsize = 4,
              #width = 10 , height= 7.69
              cellwidth = 0.5, cellheight= 0.5)
  }
  
  for (point_var in var_first){
    png(file = paste0(file_path,"/cor_median_",point_var,".png"), width = 350,height=200)
    print(ggplot(data=cor_matrix_median_info,aes_string(x=point_var,y=cor_matrix_median))+geom_point()+ylab("cor_matrix_median")+
            geom_smooth(span=0.2,se=FALSE)+theme(axis.text = element_text(size=12))+
            labs(title=paste0("Cor Median order by ",point_var))+
            theme(legend.position = "none",plot.title = element_text(hjust = 0.5),axis.title = element_text(size=15)))
    dev.off()
  } 
}

#####################################################################################
# 重要输出图形：火山图、差异基因热图、PCA图、PCA_deg图
#####################################################################################
volcano_plot<- function( res, file_name, plot_cutoff = 0.1 ){
  # volcano plot
  # Arg:
  #   res: result from DESeq analysis
  #   plot_cutoff: cutoff for log2FoldChange
  # output: 
  #   volcano plot
  resOrdered_forvolcano <- res[order(res$pvalue),]
  P.value=resOrdered_forvolcano$pvalue
  FC=resOrdered_forvolcano$log2FoldChange
  df_forvolcano=data.frame(P.value,FC)
  df_forvolcano$change=as.factor(ifelse(df_forvolcano$P.value < 0.05 & abs(df_forvolcano$FC) >=plot_cutoff,ifelse(df_forvolcano$FC > plot_cutoff,'Up','Down'),'Not'))
  df_forvolcano=droplevels(df_forvolcano)
  ggplot(data=na.omit(df_forvolcano),aes(x=FC, y =-log10(P.value),colour=change)) +
    geom_point(alpha=0.4, size=1.75) + 
    #xlim(c(-2, 2)) +
    xlab("log2 fold change") + ylab("-log10 p-value") +
    labs(title="Volcano picture")+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_vline(xintercept = c(-plot_cutoff,plot_cutoff),lty=4,col="grey",lwd=0.5)+
    geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.5)+
    scale_color_manual(values=c("blue", "grey","red"))
  ggsave(file_name,width = 5,height = 4)
}

plot_heatmap <- function( vst, gene_list, input, unique_id=NULL, file_heatmap_name ){
  # 做heatmap必须把gene_list制作完成。行列的聚类默认使用的都是欧氏距离。
  plot_mat <- vst[ gene_list, ]
  diagnosis <- data.frame( "type" = input[colnames(plot_mat), "type"] , row.names = colnames(plot_mat) )
  c1=colorRampPalette(c("blue","white"))(100)
  c2=colorRampPalette(c("white","red"))(100)
  color_now=cbind(c1,c2)
  if(dim(plot_mat)[1]!=1){
    hc <- hclust(dist(t(plot_mat)))
    callback=function(hc,plot_mat){
      sv=svd(t(plot_mat))$v[,1]
      dend=reorder(as.dendrogram(hc),wts = sv)
      as.hclust(dend)
    }
    if(!is.null(unique_id)){
      rownames(input) <- input[,unique_id]
    }
    pdf( file_heatmap_name)
    par( mar=c( 5,4,4,4 ) )
    pheatmap(plot_mat, 
             scale="row",
             # clustering_distance_row="correlation",
             color = color_now,
             annotation=diagnosis,
             show_rownames = T,
             fontsize=7,
             fontsize_row=5,
             width = 100,
             height = 100,
             clustering_callback = callback,
             cluster_cols = TRUE,
             cluster_rows = TRUE,
             show_colnames = F
    )
    dev.off()
  }
}

pca_plot <- function( rld, input, unique_id=NULL, file_name_pca, file_name_pr = NULL ){
  # 全部基因进行PCA并绘图
  # PCA plot
  # Arg:
  #   rld: DESeqTransform data, rld or vst for all genes 
  #   input : A dataframe including factor 
  #   file_name_pr: name of this pc value file for save
  # output: 
  #   PCA plot
  input$type <- factor(input$type,levels = rev(levels(input$type)))
  # filter low values
  rld_forpca=rld[rowSums(rld)>1,]
  # PCA analysis result
  pr=prcomp(t(rld_forpca))
  # summary info for each component
  percent=summary(pr)
  # PC1 and PC2
  percent=percent$importance[2,1:2]
  if(!is.null(unique_id)){
    label_list <- input[,unique_id]
  } else {
    label_list <- input$seq_id
  }
  ggplot(as.data.frame(pr$x),aes(x=PC1, y=PC2 )) +
    geom_point(aes(color=input$type)) + theme(plot.title = element_text(hjust = 0.5)) +
    #geom_text(label=label_list,size=3) +
    labs(title = "PCA plot of all gene",color='type',x=paste("PC1","(",percent[1]*100,"%",")",sep = ""),y=paste("PC2","(",percent[2]*100,"%",")",sep = ""))
  ggsave(file_name_pca,width = 5,height = 4)
  
  if (!is.null(file_name_pr)){
    pr_forsave=as.data.frame(pr$x)
    write.csv(pr_forsave, file = file_name_pr)  
  }
}

pca_plot_deg <- function( rld, gene_list, input, unique_id=NULL, file_name_pca, file_name_pr = NULL ){
  # 选出的差异基因进行PCA并绘图
  # PCA plot for DE genes
  # Arg:
  #   rld: DESeqTransform data, rld or vst for DE genes 
  #   gene_list: top 100 DE gene list 
  #   input : A dataframe including factor 
  #   file_name_pr: name of file for save
  input$type <- factor(input$type,levels = rev(levels(input$type)))
  rld_forpca=rld[ gene_list ,]
  pr=prcomp(t(rld_forpca))
  percent=summary(pr)
  percent=percent$importance[2,1:2]
  if(!is.null(unique_id)){
    label_list <- input[,unique_id]
  } else {
    label_list <- input$seq_id
  }
  ggplot(as.data.frame(pr$x),aes(x=PC1, y=PC2))+
    geom_point(aes(color=input$type)) + theme(plot.title = element_text(hjust = 0.5)) +
    #geom_text(label=label_list,size=3) +
    labs(title = "PCA plot of differential gene",color='type',x=paste("PC1","(",percent[1]*100,"%",")",sep = ""),y=paste("PC2","(",percent[2]*100,"%",")",sep = ""))
  ggsave( file_name_pca,width = 5,height = 4)
  
  if (!is.null(file_name_pr)){
    pr_forsave=as.data.frame(pr$x)
    write.csv(pr_forsave, file = file_name_pr) 
  }
}

#####################################################################################
# 生信分析图形：GO图、KEGG图
#####################################################################################
GOanalysis_plot <- function(res,file_name,flag="all",data_file="GO_detail.csv"){
  # 采用pvalue<0.05的基因列表。
  # GO图绘制，可以单独选择bp、cc、mf，默认全部绘制到一张直方图上
  # Arg:
  #   flag:只可以为"all"、"bp"、"cc"、"mf"。
  # output: GO直方图、绘图数据文件
  local_flag <- c(FALSE,FALSE,FALSE)
  mode <- switch(flag, all = 1, bp = 2, cc=3, mf=4)
  res_new <- res[ !is.na(res$pvalue), ]
  resOrdered <- res_new[ res_new$pvalue <0.05 , ]
  gene_list <- rownames(resOrdered)
  # Get gene entrez id
  keytypes(org.Hs.eg.db)
  keys(org.Hs.eg.db,keytypes="ENTREZID") %>% head
  gene_list_entrezid=AnnotationDbi::select(org.Hs.eg.db,keys = gene_list ,columns = c("ENTREZID") , keytype = "SYMBOL")
  gene_list_entrezid=gene_list_entrezid$ENTREZID
  # extract unique genes
  target_gene_id <- unique(gene_list_entrezid)
  display_number = c(10, 10, 10)
  if (as.numeric(mode) == 1 | as.numeric(mode) == 2){
    ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                       gene = target_gene_id,
                       pvalueCutoff = 0.05,
                       ont = "BP",
                       readable=TRUE)
    ego_result_BP <- na.omit(as.data.frame(ego_BP)[1:display_number[3], ])
    ego_result_BP <- ego_result_BP[order(ego_result_BP$pvalue),]
    local_flag[1]=TRUE
    if(dim(ego_result_BP)[1] < display_number[1]){
      display_number[1] <- dim(ego_result_BP)[1]
    }
    go_enrich_df <- data.frame(ID=ego_result_BP$ID,
                               Description=ego_result_BP$Description,
                               GeneCount=ego_result_BP$Count,
                               pvalue=ego_result_BP$pvalue,
                               type=factor(rep("Biological Process", display_number[1])))
    CPCOLS <- "#66C3A5"
  }
  if (as.numeric(mode) == 1 | as.numeric(mode) == 3){
    ego_CC <- enrichGO(OrgDb="org.Hs.eg.db",
                       gene = target_gene_id,
                       pvalueCutoff = 0.05,
                       ont = "CC",
                       readable=TRUE)
    ego_result_CC <- na.omit(as.data.frame(ego_CC)[1:display_number[2], ])
    ego_result_CC <- ego_result_CC[order(ego_result_CC$pvalue),]
    local_flag[2]=TRUE
    if(dim(ego_result_CC)[1] < display_number[2]){
      display_number[2] <- dim(ego_result_CC)[1]
    }
    go_enrich_df <- data.frame(ID=ego_result_CC$ID,
                               Description=ego_result_CC$Description,
                               GeneCount=ego_result_CC$Count,
                               pvalue=ego_result_CC$pvalue,
                               type=factor(rep("Cellular Component", display_number[2])))
    CPCOLS <- "#FD8D62"
    
  }
  if (as.numeric(mode) == 1 | as.numeric(mode) == 4){
    ego_MF <- enrichGO(OrgDb="org.Hs.eg.db",
                       gene = target_gene_id,
                       pvalueCutoff = 0.05,
                       ont = "MF",
                       readable=TRUE)
    ego_result_MF <- na.omit(as.data.frame(ego_MF)[1:display_number[1], ])
    ego_result_MF <- ego_result_MF[order(ego_result_MF$pvalue),]
    local_flag[3]=TRUE
    if(dim(ego_result_MF)[1] < display_number[3]){
      display_number[3] <- dim(ego_result_MF)[1]
    }
    go_enrich_df <- data.frame(ID=ego_result_MF$ID,
                               Description=ego_result_MF$Description,
                               GeneCount=ego_result_MF$Count,
                               pvalue=ego_result_MF$pvalue,
                               type=factor(rep("Molecular Function", display_number[3])))
    CPCOLS <- "#8DA1CB"
  }
  if (local_flag[1] & local_flag[2] & local_flag[3]) {
    go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
                               Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
                               GeneCount=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
                               pvalue=c(ego_result_BP$pvalue,ego_result_CC$pvalue, ego_result_MF$pvalue),
                               type=factor(c(rep("Biological Process", display_number[1]), rep("Cellular Component", display_number[2]),
                                             rep("Molecular Function", display_number[3])), levels=c("Biological Process", "Cellular Component", "Molecular Function")))
    # colors for bar // green,orange,blue
    CPCOLS <- c("#66C3A5", "#FD8D62", "#8DA1CB")
  }
  go_enrich_df <- go_enrich_df[!is.na(go_enrich_df[,1]),]
  if ( nrow(go_enrich_df)==0 ) {
    print (paste("no",toupper(flag),"enriched GO terms"))
  }else{
    write.csv(go_enrich_df,file = data_file)
    go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
    go_enrich_df$pvalue <- format(go_enrich_df$pvalue,digit=3)
    shorten_names <- function(x, n_word=6, n_char=60){
      if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 60))
      {
        if (nchar(x) > 60) x <- substr(x, 1, 60)
        x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                         collapse=" "), "...", sep="")
        return(x)
      }
      else
      {
        return(x)
      }
    }
    labels=(sapply(
      levels(go_enrich_df$Description)[as.numeric(go_enrich_df$Description)],
      shorten_names))
    #names(labels) = rev(1:nrow(go_enrich_df))
    ggplot(data=go_enrich_df, aes(x=number, y=GeneCount, fill=type)) +
      geom_bar(stat="identity", width=0.8) + coord_flip() +
      scale_fill_manual(values = CPCOLS) + theme_bw() +
      scale_x_discrete(labels=rev(go_enrich_df$pvalue)) +
      geom_text(aes(label = labels,y=0),hjust=0,cex=3)+
      labs(x="p-value",title = "The Most Enriched GO Terms", subtitle="differential gene:pvalue<0.05") +
      theme(panel.grid =element_blank(),plot.title=element_text(hjust=0.5),plot.subtitle=element_text(size=10, hjust=1))
      #theme(axis.text=element_text(face = "bold", color="gray50"))
    ggsave(file_name,width = 7,height = 6) 
  }
}

kegg_bubble_plot <- function(res, pathway_num, file_name, pvalueCut=0.9, qvalueCut=0.9){
  # 采用pvalue<0.05的基因列表。
  # ID convertion
  library(org.Hs.eg.db)
  library(annotate)
  library(clusterProfiler)
  require(pathview)
  res_new <- res[ !is.na(res$pvalue), ]
  resOrdered <- res_new[ res_new$pvalue <0.05 , ]
  gene_list <- rownames(resOrdered)
  gene_list_entrezid_df = AnnotationDbi::select( org.Hs.eg.db , keys = gene_list ,columns = c( "ENTREZID" ) , keytype = "SYMBOL" )
  gene_list_entrezid = gene_list_entrezid_df$ENTREZID
  enrichkegg=enrichKEGG(gene=gene_list_entrezid,organism = "hsa",pvalueCutoff = pvalueCut, pAdjustMethod = "BH",qvalueCutoff = qvalueCut ) #
  enrichkegg_summary <- as.data.frame( enrichkegg )
  enrichkegg_summary <- enrichkegg_summary[order(enrichkegg_summary$pvalue),]
  plot_data <- data.frame(enrichkegg_summary)
  if(nrow(plot_data)==0){
    print (paste("no enriched pathway"))
    } else {
      for(i in 1:nrow(plot_data)){
        plot_data$All_Unigene[i] <- unlist(strsplit(plot_data$BgRatio[i],split="/"))[1]
        plot_data$richFactor[i] <- plot_data$Count[i]/as.numeric(plot_data$All_Unigene[i])
        }
      if(length(plot_data$ID)>=pathway_num){
        plot_data <- plot_data[1:pathway_num,]
        } else {
          plot_data <- plot_data[1:nrow(plot_data),]
          }
      plot_data$number <- factor(rev(1:nrow(plot_data)))
      #plot_data$de <- paste0(plot_data$Description,plot_data$pvalue)
      ggplot(plot_data,aes(richFactor,number))+ 
        geom_point(aes(size=Count,color=-1*log10(qvalue))) +
        scale_y_discrete(labels=rev(plot_data$Description)) +
        #scale_colour_gradient(low="green",high="red") +
        labs(color=expression(-log[10](qvalue)),size="Gene number",x="Rich factor",y=" ",
             title=paste0("Top",dim(plot_data)[1]," of pathway enrichment"),subtitle="differential gene:pvalue<0.05") + 
        theme_bw() + theme(panel.grid =element_blank(),plot.title=element_text(hjust=0.5),plot.subtitle=element_text(size=10, hjust=1))
      ggsave(file_name,width = 6.5,height = 6)
      save( enrichkegg_summary , file = "enrichkegg.rda")
      write.csv(enrichkegg_summary[,"ID"!=colnames(enrichkegg_summary)],file = "KEGG_data.csv")
    }
  }

keggplot <- function( res, gene_list ){
  # 采用提供的gene_list列表。
  # ID convertion
  library(org.Hs.eg.db)
  library(annotate)
  library(clusterProfiler)
  require(pathview)
  gene_list_entrezid_df = AnnotationDbi::select( org.Hs.eg.db , keys = gene_list ,columns = c( "ENTREZID" ) , keytype = "SYMBOL" )
  # KEGG figure
  enrichkegg=enrichKEGG(gene=gene_list_entrezid_df$ENTREZID,organism = "hsa",pvalueCutoff = 0.9, pAdjustMethod = "BH",qvalueCutoff = 0.9 )
  enrichkegg_summary <- as.data.frame( enrichkegg )
  gene_fc<-data.frame(log2FoldChange=res[match(gene_list,rownames(res)),]$log2FoldChange)
  gene_fc$SYMBOL <- rownames(res[match(gene_list,rownames(res)),])
  gene_fc <- inner_join(gene_fc,gene_list_entrezid_df)
  gene_fc <- na.omit(gene_fc)
  genedata <- as.matrix(gene_fc[,1])
  rownames(genedata) <- gene_fc$ENTREZID
  if(nrow(enrichkegg_summary)==0){
    print (paste("no enriched pathway"))
  } else {
    if(length(enrichkegg_summary$ID)>=3){
      for (i in 1:3){
        trynext=try(pathview(gene.data = genedata,pathway.id = enrichkegg_summary$ID[i],species="hsa"),silent=T)
        if ('try-error' %in% class(trynext)) next
        }
      } else {
        for (i in 1:length(enrichkegg_summary$ID)){
          trynext=try(pathview(gene.data = genedata,pathway.id = enrichkegg_summary$ID[i],species="hsa"),silent=T)
          if ('try-error' %in% class(trynext)) next
          }
        }
    save( enrichkegg_summary , file = "enrichkegg.rda")
    write.csv(enrichkegg_summary[,"ID"!=colnames(enrichkegg_summary)],file = "KEGG_data.csv")
    }
  }

#####################################################################################
# 建模判分图形：ROC曲线AUC图、logit图
#####################################################################################
auc_plot <- function(... , data_name, type_list, file_name){
  # ROC曲线绘制到一张图里进行对比
  # input：
  #   ...：判分结果数据，如score1,score2...
  #   data_name：与...的一一对应，如c("score1_AUC","score2_AUC")
  #   file_name：保存的文件名
  input <- list(...)
  names(input) <- data_name
  pdf(file=file_name,width = 4, height = 4)
  if (length(input) == 1){
    input[[1]]$type <- factor(input[[1]]$type,levels=type_list)
    plot.roc(input[[1]]$type,input[[1]]$score ,main = "ROC curve and AUC value" )
    theauc <- pROC::auc( input[[1]]$type,input[[1]]$score )
    youden <-round( unname(coords(roc(input[[1]]$type,input[[1]]$score),"b", ret=c("sen","spec","t","accu"),transpose = TRUE)),3)
    spec_sens <- data.frame( "spec" = attr(theauc , "roc")$specificities , "sens" = attr(theauc , "roc")$sensitivities , "thres" = attr(theauc , "roc")$thresholds)
    legend("bottomright",legend=paste0(names(input)[1]," AUC:",round(theauc,3)),cex=0.7)
    #legend("bottomright",legend=paste0(names(input)[1]," AUC:",round(theauc,3),",Ysens:",youden[1],",Yspec:",youden[2],",Ythres:",youden[3]))
    } else {
      input[[1]]$type <- factor(input[[1]]$type,levels=type_list)
      p1 <- plot.roc(input[[1]]$type,input[[1]]$score,col="1",main = "ROC curve and AUC value")
      theauc <- pROC::auc( input[[1]]$type,input[[1]]$score )
      youden <-round( unname(coords(roc(input[[1]]$type,input[[1]]$score),"b", ret=c("sen","spec","t","accu"),transpose = TRUE)),3)
      spec_sens <- data.frame( "spec" = attr(theauc , "roc")$specificities , "sens" = attr(theauc , "roc")$sensitivities , "thres" = attr(theauc , "roc")$thresholds)
      legend_content <- paste0(names(input)[1]," AUC:",round(theauc,3))
      #legend("bottomright",legend=paste0(names(input)[1]," AUC:",round(theauc,3),",Ysens:",youden[1],",Yspec:",youden[2],",Ythres:",youden[3]))
      for ( i in 2:length(input) ){
        input[[i]]$type <- factor(input[[i]]$type,levels=type_list)
        lines.roc(input[[i]]$type,input[[i]]$score,col=as.character(i))
        theauc <- pROC::auc( input[[i]]$type,input[[i]]$score)
        youden <-round( unname(coords(roc(input[[i]]$type,input[[i]]$score),"b", ret=c("sen","spec","t","accu"),transpose = TRUE)),3)
        spec_sens <- data.frame( "spec" = attr(theauc , "roc")$specificities , "sens" = attr(theauc , "roc")$sensitivities , "thres" = attr(theauc , "roc")$thresholds )
        legend_content <- c(legend_content,paste0(names(input)[i]," AUC:",round(theauc,3)))
        #legend("bottomright",legend=paste0(names(input)[i]," AUC:",round(theauc,3),",Ysens:",youden[1],",Yspec:",youden[2],",Ythres:",youden[3]))
        }
      legend("bottomright",legend = legend_content,col=as.character(1:length(input)),lwd=2,cex=0.7)
    }
  dev.off()
  }

logit_plot <- function(score_list, plot_name, file_name){
  # logit plot
  # input:
  #   score_list：判分函数获取到的判分结果
  #   plot_name：图名
  score_list <- score_list[order(score_list$score),]
  score_list$type <- as.character(score_list$type)
  score_list <- score_list[order(score_list$type),]
  score_list$number <- (1:nrow(score_list))
  ggplot(score_list,aes(x=number,y=score,colour=type)) + 
    facet_grid(~type) +
    labs(title=plot_name) +
    geom_point() +
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
  ggsave(file_name,width = 4,height = 4)
  }

#####################################################################################
# 其他输出图形：根据样本类型分类的correlation图、cluster图
#####################################################################################
plot_sample_cor <- function(input_data,countData,file_name,cutoff=0.97){
  # 根据所有基因上featurecount做样本间的correlation图。采用皮尔逊距离计算相关性。
  input_data <- input_data[order(input_data$type),]
  cor_matrix <- cor(countData)
  cor_matrix[cor_matrix<cutoff] <- cutoff
  pheatmap(cor_matrix[rownames(input_data),rownames(input_data)],main = "Correlation Heatmap",
           annotation_col = data.frame(row.names = rownames(input_data),"type"=input_data$type),
           annotation_row = data.frame(row.names = rownames(input_data),"type"=input_data$type),
           annotation_names_row = FALSE , annotation_names_col = FALSE, 
           cluster_cols = FALSE,cluster_rows = FALSE,fontsize = 7,
           show_rownames = F,show_colnames = F,filename = file_name)
  }

plot_clust <- function(vst,colData,file_name,dist_num=1,clust_num=3){
  # 采用经DESeq2处理的差异基因结果（vst）做聚类图，此处用的是欧氏距离。
  dist_method <- list("euclidean","maximum","manhattan","canberra","minkowski","binary ")
  clust_method <- list("ward.D","single","complete","median","mcquitty","average","centroid")
  clust_input <- t(vst)
  Samples <- dist(clust_input,method = dist_method[dist_num])
  out_hclust <- hclust(Samples,method = clust_method[clust_num])
  
  colLab <- function(n) {
    if (is.leaf(n)) {
      a <- attributes(n)
      labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
      attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
      }
    n
    }
  hcd <-  as.dendrogram(out_hclust)
  labelColors = c("#33CCFF", "#FF3399")
  clusMember = as.integer(factor(colData$type,labels=c(1,2)))
  names(clusMember) <- rownames(colData)
  clusDendro = dendrapply(hcd, colLab)
  # you can revise the size of the labels:
  #library(dendextend)
  #clusDendro <- set(clusDendro, "labels_cex", 0.5)
  # make plot
  pdf(file_name)
  plot(clusDendro, main = "clustering")
  legend("topright",legend = unique(colData$type),fill=labelColors)
  dev.off()
  }

#####################################################################################
# 其他分析图形：按不同因素分组的判分的95%CI图（供参考）
#####################################################################################
plot_score_ci <- function(... , data_name,ylab_name='',file_name){
  # input：
  #   ...：分组的判分结果数据，如data1,data2...
  #   data_name：与...的一一对应，如c("data1","data2")
  #   ylab_name：纵轴名称
  #   file_name：保存的文件名，推荐使用png格式
  input <- list(...)
  names(input) <- data_name
  png(file=file_name)
  if (length(input) == 1){
    t_input <- t.test(input[[1]]$score)
    plot(0,xlim=c(t_input$conf.int[1]-0.1,t_input$conf.int[2]+0.1),ylim=c(0,length(data_name)),type="n",xlab="score(95%CI)",ylab=ylab_name,yaxt="n",bty="l")
    lines(x=c(t_input$conf.int[1],t_input$conf.int[2]),y=rep(1,2))
    lines(x=rep(t_input$conf.int[1],2),y=c(0.97,1.03))
    lines(x=rep(t_input$conf.int[2],2),y=c(0.97,1.03))
    points(x=mean(c(t_input$conf.int[1],t_input$conf.int[2])),y=rep(1),pch=16)
    axis(side=2,las=2,at=seq(1:length(data_name)),labels = data_name)
    abline(v=mean(c(t_input$conf.int[1],t_input$conf.int[2])),col="red",lty=2)
  } else {
    t_df <- list()
    x_loc1 <- c()
    x_loc2 <- c()
    for (i in 1:length(input)){
      t_sub <- t.test(input[[i]]$score)
      t_df[[i]] <- t_sub
      x_loc1 <- c(x_loc1,t_sub$conf.int[1])
      x_loc2 <- c(x_loc1,t_sub$conf.int[2])
    }
    plot(0,xlim=c(min(x_loc1)-0.1,max(x_loc1)+0.1),ylim=c(0,length(data_name)),type="n",xlab="score(95%CI)",ylab=ylab_name,yaxt="n",bty="l")
    for (i in rev(1:length(input))){
      j <- length(input)+1-i
      lines(x=c(t_df[[i]]$conf.int[1],t_df[[i]]$conf.int[2]),y=rep(j,2))
      lines(x=rep(t_df[[i]]$conf.int[1],2),y=c(j-0.03,j+0.03))
      lines(x=rep(t_df[[i]]$conf.int[2],2),y=c(j-0.03,j+0.03))
      points(x=mean(c(t_df[[i]]$conf.int[1],t_df[[i]]$conf.int[2])),y=rep(j),pch=16)
      abline(v=mean(c(t_df[[i]]$conf.int[1],t_df[[i]]$conf.int[2])),col="red",lty=2)
    }
    axis(side=2,las=2,at=seq(1:length(data_name)),labels = rev(data_name))
  }
  dev.off()
}

#####################################################################################
# 其他分布图形：年龄、性别、分期分布图
#####################################################################################
plot_age_distr <- function(origin_input,type_list,file_name){
  if("mold" %in% colnames(origin_input)){
    origin_data_T <- origin_input[origin_input$grouping == "T1V1",]
    origin_data_V <- origin_input[origin_input$grouping == "V2",]
    origin_data_T$mold <- factor(origin_data_T$type,levels = type_list,labels = paste0(type_list,"_training"))
    origin_data_V$mold <- factor(origin_data_V$type,levels = type_list,labels = paste0(type_list,"_validation"))
    revised_data <- rbind( origin_data_T,origin_data_V)
    revised_data$type <- factor(revised_data$type,levels = type_list)
    revised_data$grouping <- factor(revised_data$grouping,levels = c("T1V1","V2"),labels = c("Training","Validation"))
    calculate_data <- revised_data[,c("mold","age")]
    result_mean <- aggregate(.~mold,calculate_data, mean)
    result_sd <- aggregate(.~mold,calculate_data, sd)
    ggplot(revised_data, aes(x = age, fill = mold)) +
      #facet_grid(grouping~type)
      facet_wrap(~mold,ncol=2) + labs(title="Age Distribution") + 
      theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) +
      geom_histogram(position = "identity",alpha = 0.2, bin = 0.5,aes(y = ..density..)) +
      geom_line(stat = "density",lty=2,lwd=0.3,col = "grey18") +
      geom_rug(col="grey", sides = "b") + scale_fill_manual(values = c("red","blue","orange","green")) +
      geom_text(x=-Inf, y=Inf,vjust=1,hjust=-0.1,label=paste0("mean:",round(result_mean[,2],1)),data=result_mean,cex=3)+
      geom_text(x=-Inf, y=Inf,vjust=3,hjust=-0.2,aes(label=paste0("sd:",round(result_sd[,2],1))),data=result_sd,cex=3)
      #annotate("text", x=-Inf, y=Inf,vjust=1,hjust=-0.1, label=paste0("mean:",round(result_mean[,4],1))) +
      #annotate("text", x=-Inf, y=Inf,vjust=3,hjust=-0.2, label=paste0("sd:",round(result_sd[,4],1)))
    ggsave(file_name)
    } else {
      origin_input$type <- factor(origin_input$type,levels = type_list)
      calculate_data <- origin_input[,c("type","age")]
      result_mean <- aggregate(.~type,calculate_data, mean)
      result_sd <- aggregate(.~type,calculate_data, sd)
      ggplot(origin_input, aes(x = age, fill = type)) + 
        facet_grid(~type) + labs(title="Age Distribution") + 
        theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) +
        geom_histogram(position = "identity",alpha = 0.2,bin = 0.5, aes(y = ..density..)) +
        geom_line(stat = "density",lty=2,lwd=0.3,col = "grey18") +
        geom_rug(col="grey", sides = "b") + scale_fill_manual(values = c("red","blue","orange","green")) +
        geom_text(x=-Inf, y=Inf,vjust=1,hjust=-0.1,aes(label=paste0("mean:",round(result_mean[,2],1))),data=result_mean,cex=3)+
        geom_text(x=-Inf, y=Inf,vjust=3,hjust=-0.2,aes(label=paste0("sd:",round(result_sd[,2],1))),data=result_sd,cex=3)
        #annotate("text", x=-Inf, y=Inf,vjust=1,hjust=-0.1, label=paste0("mean:",round(result_mean[,2],1))) +
        #annotate("text", x=-Inf, y=Inf,vjust=3,hjust=-0.2, label=paste0("sd:",round(result_sd[,2],1)))
      ggsave(file_name,width = 7,height = 4)
    }
  }

plot_gender_distr <- function(origin_input,type_list,file_name){
  if("mold" %in% colnames(origin_input)){
    origin_data_T <- origin_input[origin_input$grouping == "T1V1",]
    origin_data_V <- origin_input[origin_input$grouping == "V2",]
    t_result <- t.test(as.numeric(origin_data_T$gender),as.numeric(origin_data_V$gender))
    TV_pvalue <- t_result$p.value
    origin_data_T$mold <- factor(origin_data_T$type,levels = type_list,labels = paste0(type_list,"_training"))
    origin_data_V$mold <- factor(origin_data_V$type,levels = type_list,labels = paste0(type_list,"_validation"))
    revised_data <- rbind( origin_data_T,origin_data_V)
    revised_data$type <- factor(revised_data$type,levels = type_list)
    revised_data$grouping <- factor(revised_data$grouping,levels = c("T1V1","V2"),labels = c("Training","Validation"))
    revised_data$gender <- factor(revised_data$gender,levels = c("1","0"),labels = c("male","female"))
    calculate_data <- revised_data[,c("gender","mold","type","grouping")]
    result_data <- as.data.frame(table(calculate_data))
    colnames(result_data) <- c("gender","mold","type","grouping","num")
    result_data <- result_data[result_data$num != 0,]
    plot_data <- ddply(result_data,.(mold),transform,percent=num/sum(num))
    plot_data <- ddply(plot_data,.(mold),transform,position=sum(percent)-(percent/2 + c(0, cumsum(percent)[-length(percent)])))
    ggplot(plot_data,aes(x="",y=percent,fill=gender))+
      facet_wrap(~mold,ncol=2)+
      #facet_grid(grouping~type)+
      geom_bar(stat = "identity",position="stack",width = 1)+coord_polar(theta = "y")+
      labs(x = "", y = "",title=paste0("Gender Distribution","\n","TV_pvalue=",round(TV_pvalue,3))) +
      theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.grid=element_blank()) +
      scale_fill_manual(values = c("lightblue","pink")) +
      theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+
      geom_text(aes(x=1,y=plot_data$position,label=paste(gender,"\n",round(percent*100,1),"%")),cex=3)
    ggsave(file_name)
    } else {
      origin_input$type <- factor(origin_input$type,levels = type_list)
      origin_data_1 <- origin_input[origin_input$type == levels(origin_input$type)[1],]
      origin_data_2 <- origin_input[origin_input$type == levels(origin_input$type)[2],]
      t_result <- t.test(as.numeric(origin_data_1$gender),as.numeric(origin_data_2$gender))
      type_pvalue <- t_result$p.value
      origin_input$gender <- factor(origin_input$gender,levels = c("1","0"),labels = c("male","female"))
      calculate_data <- origin_input[,c("gender","type")]
      result_data <- as.data.frame(table(calculate_data))
      colnames(result_data) <- c("gender","type","num")
      result_data <- result_data[result_data$num != 0,]
      plot_data <- ddply(result_data,.(type),transform,percent=num/sum(num))
      plot_data <- ddply(plot_data,.(type),transform,position=sum(percent)-(percent/2 + c(0, cumsum(percent)[-length(percent)])))
      ggplot(plot_data,aes(x="",y=percent,fill=gender))+
        facet_grid(~type)+
        geom_bar(stat = "identity",position="stack",width = 1)+coord_polar(theta = "y")+
        labs(x = "", y = "",title=paste0("Gender Distribution","\n","type_pvalue=",round(type_pvalue,3))) +
        theme(legend.position = "none",plot.title = element_text(hjust = 0.5),panel.grid=element_blank()) +
        scale_fill_manual(values = c("lightblue","pink")) +
        theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+
        geom_text(aes(x=1,y=plot_data$position,label=paste(gender,"\n",round(percent*100,1),"%")),cex=3)
      ggsave(file_name)
    }
  }

plot_stage_distr <- function(origin_stage,file_name){
  # origin_stage里只包括seq_num、stage和type即可，此处type为“CAC_T”和“CAC_V”
  # 该图片需保存为png格式！！！用pdf格式将不能正确显示分期或会文件错误。
  origin_stage[is.na(origin_stage)] <- "Unknown"
  origin_stage[origin_stage==""] <- "Unknown"
  if("type" %in% colnames(origin_stage)){
    #origin_stage$mold <- factor(origin_stage$type,levels = c("CAC_T","CAC_V"),labels = c("Cancer_training","Cancer_validation"))
    stage_data<- origin_stage[,c("stage","type")]
    stage_data <- as.data.frame(table(stage_data))
    colnames(stage_data) <- c("stage","type","num")
    plot_data <- ddply(stage_data,.(type),transform,percent=round(num/sum(num) *100,1))
    plot_data$label <- paste0(plot_data$percent,"%")
    ggplot(plot_data,aes(x=stage,y=percent,fill=type)) + 
      facet_grid(~type) + labs(title="Stage Distribution") +
      theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) +
      geom_bar(stat="identity",position="dodge",width=0.7) + 
      geom_text(aes(label = label),vjust = -0.5,position = position_dodge(.9),cex=3) +
      scale_fill_manual(values = c("pink","lightblue"))
    ggsave(file_name,width = 7, height = 4)
    } else {
      origin_stage[is.na(origin_stage)] <- "Unknown"
      origin_stage[origin_stage==""] <- "Unknown"
      plot_data <- as.data.frame(table(origin_stage$stage))
      colnames(plot_data) <- c("stage","num")
      plot_data$percent <- round((plot_data$num / sum(plot_data$num))*100,1)
      plot_data$label <- paste0(plot_data$percent,"%")
      ggplot(plot_data,aes(x=stage,y=percent)) + 
        labs(title="Stage Distribution") +
        theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) +
        geom_bar(stat="identity",position="dodge",width=0.7) + 
        geom_text(aes(label = label),vjust = -0.5,position = position_dodge(.9),cex=3)
      ggsave(file_name,width = 4, height = 4)
    }
  }
