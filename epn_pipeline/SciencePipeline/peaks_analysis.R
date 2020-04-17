#####################################################################################
# Peaks相关分析
#####################################################################################
#此函数用于2~5个样本的peaks相关分析，使用了findOverlapsOfPeaks和ChIPseeker两个package中函数，分析完成可得8张图表
#最后四行为2~5个样本的分析示例
#图表说明：
#OverlappingPeaksVenn：以韦恩图可视化样本之间peaks的overlap
#AvgProf-samples：以谱图的形式展示所有样本的peaks在TSS区域的联合结合强度
#AvgProf-combined:以谱图的形式展示各个样本的peaks在TSS区域的结合强度
#AvgProf-splited：在上一张图的基础上，以色带可视化各个样本的peaks在TSS区域结合强度的置信区间
#tagHeatmap：看peaks在某个窗口上的结合谱图。在一个固定的窗口里（ex:启动子区域，使用转录起始位点，然后指定上下游，再把peaks比对到这个窗口）把peaks全部排列起来，生成矩阵供后续分析和可视化
#FeatureDistBar：可视化peaks在features上比对情况的占比
#DistToTSS：可视化离peaks最近基因的距离分布
#AnnoGenesVenn：看注释所得最近基因在不同样本中的overlap

library(ChIPpeakAnno)
library(GenomicRanges)
library(GenomeInfoDb)
library(VennDiagram)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(DOSE)


Peaks_Analysis <- function(input_file , data_name, file_name,upstream=3000, downstream=3000){
  # input_file： peak结果数据，文件第一列为seq_id(如SEQ01、SEQ02),第二列为location(.broadPeak文件路径)
  # data_name：与...的一一对应，如c("SEQ01","SEQ02")
  # file_name：保存的文件名
  # upstream:TSS区域上游距离
  # downstream:TSS区域下游距离

  input<-as.character(read.table( input_file , head=FALSE , sep = "\t" , skip = 1)[,2])
  peaklist <- list()
  name_list<-c()
  color=c("orange", "blue", "green", "red", "yellow")
  
  for (i in 1:length(input)){
    temp<-readPeakFile(input[i], head=FALSE)
    peaklist<-c(peaklist, temp)
    name_list<-c(name_list,data_name[i])
  }
  names(peaklist)[1:length(input)]<-name_list
  
  #///绘制样本间overlapping peaks韦恩图"
  result<-findOverlapsOfPeaks(unlist(peaklist))
  pdf(file=paste0(file_name,"-OverlappingPeaksVenn",".pdf"))
  trynext=try(makeVennDiagram(result,fill=color[1:length(peaklist)],NameOfPeaks=data_name)+title(file_name),silent=T)
  if ('try-error' %in% class(trynext)) {dev.off()}
  
  #///样本间peaks的其他比较
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  promoter <- getPromoters(TxDb=txdb,upstream=upstream,downstream=downstream)
  
  #plot-samples Average Profile of ChIP peaks binding to TSS region
  merge.data = read.table(file =input[1],header=FALSE)
  for (i in 2:length(input)){
    new.data = read.table(file = input[i], header=FALSE)
    merge.data = rbind(merge.data,new.data)
  }
  write.table(merge.data,file = "merge_peaks.txt",row.names=F,col.names = FALSE,quote=FALSE,sep = "\t")  
  samples_Matrix <- getTagMatrix("merge_peaks.txt", windows=promoter)
  pdf(file=paste0(file_name,"-AvgProf-samples",".pdf"))
  print(plotAvgProf(samples_Matrix, xlim=c(-downstream, upstream)))
  dev.off()
  
  #plot-combined Average Profile of ChIP peaks binding to TSS region
  tagMatrixList <- lapply(peaklist, getTagMatrix, windows=promoter)
  pdf(file=paste0(file_name,"-AvgProf-combined",".pdf"))
  print(plotAvgProf(tagMatrixList, xlim=c(-downstream, upstream)))
  dev.off()
  
  #plot-splited Average Profile of ChIP peaks binding to TSS region
  pdf(file=paste0(file_name,"-AvgProf-splited",".pdf"))
  print(plotAvgProf(tagMatrixList, xlim=c(-downstream, upstream), conf=0.95,resample=500, facet="row"))
  dev.off()
  
  #plot-Heatmap
  pdf(file=paste0(file_name,"-tagHeatmap",".pdf"))
  print(tagHeatmap(tagMatrixList, xlim=c(-downstream, upstream), color=NULL))
  dev.off()
  
  #///ChIP peak annotation comparision
  peakAnnoList <- lapply(peaklist, annotatePeak, TxDb=txdb,
                         tssRegion=c(-downstream, upstream), verbose=FALSE)
  
  #plot-feature distribution
  pdf(file=paste0(file_name,"-FeatureDistBar",".pdf"))
  print(plotAnnoBar(peakAnnoList))
  dev.off()
  
  #plot-Distribution of transcription factor-binding loci relative to TSS
  pdf(file=paste0(file_name,"-DistToTSS",".pdf"))
  print(plotDistToTSS(peakAnnoList))
  dev.off()
  
  #///Overlap of peaks and annotated genes
  genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
  #plot-venn
  pdf(file=paste0(file_name,"-AnnoGenesVenn",".pdf"))
  print(vennplot(genes))
  dev.off()
  
}  
