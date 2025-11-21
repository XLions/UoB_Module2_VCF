#清理内存，设置工作目录
#Clear memory and set working directory
setwd('~/UoB_Work/M2_Group_vcf/')
rm(list = ls());gc()

#加载R包
#Load R Packages
library(vcfR)
library(stringr)
library(progress)
library(tidyverse)
library(biomaRt)
library(ggrepel)
library(RCircos)
library("rtracklayer")
library(dplyr)
library(ggplot2)
data("UCSC.HG38.Human.CytoBandIdeogram")

#定义必要函数
#Difine Necessary Funnctions
##判断字符串是否能为数字
##Check if a string can be converted to a number
is_numeric<-function(values){
  assign('result',c())
  for (i in 1:length(values)) {
    if(grepl("^-?(\\d\\.)?\\d+$",values[i])){
      result<-c(result,T)
    }else{
      result<-c(result,F)
    } 
  }
  return(result)
}

#A.数据处理与初步探索
#A.Data Processing & Preliminary Exploration
#1.读取数据
#1.Read data
{
  ##突变数据
  ##Mutation data
  vcf <- read.table('./Data/DellyVariation.vcf', comment.char = '#' )
  colnames(vcf)<-
    c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","HG00512",
      "HG00513","HG00514","HG00731","HG00732","HG00733","NA19238","NA19239","NA19240")
  vcf <- vcf %>%
    dplyr::select(-c("HG00512","HG00513","HG00514","HG00731","HG00732","HG00733"))
  
  #注释
  #Annotation
  #读取来源于GENCODE的基因组注释数据
  #Read genome annotation data from GENCODE
  ann_data_gff3<-read.table('./Data/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz')
  #按照GFF3标准格式重命名列
  #Rename the column according to GFF3 standard format
  colnames(ann_data_gff3)<-c(
    'seqid','source','type','start','end','score','strand','phase','attributes'
  )
  ##定义函数：拆分attributes
  ##Define function: Split attributes
  ##输出：包含向量Item和Value及对应数量的数据框
  #Output: A data box containing vectors Item, Value, and their corresponding quantities
  separateAttributesCharToListIncludingItemValue<-function(attribute_char){
    singleChars<-stringr::str_split_1(attribute_char,';')
    assign('Items',c())
    assign('Values',c())
    for(i in 1:length(singleChars)){
      Items<-c(Items,
               stringr::str_split_1(singleChars[i],'=')[1])
      Values<-c(Values,
                stringr::str_split_1(singleChars[i],'=')[2])
    }
    return(data.frame(Item=Items,Value=Values,Length=length(singleChars)))
  }
  #注释vcf【【【【十分费时间】】】】
  #Annotate VCF (computationally intensive process)
  vcf$SYMBOL<-NA
  pb <- progress_bar$new(total = nrow(vcf))#创建一个进度条对象 #Create a progress bar object
  for(i in 1:nrow(vcf)){
    selectedAnnoData<-
      ann_data_gff3 %>%
      dplyr::filter(seqid==vcf$CHROM[i]) %>%
      dplyr::filter(start<=vcf$POS[i]) %>%
      dplyr::filter(end>=vcf$POS[i]) %>%
      dplyr::filter(type %in% c('gene'))
    if(nrow(selectedAnnoData)==0){
      selectedAttributes<-NA
    }else{
      selectedAttributes<-
        (separateAttributesCharToListIncludingItemValue(selectedAnnoData$attributes[1]) %>%
        dplyr::filter(Item=='gene_name'))$Value
    }
    vcf$SYMBOL[i]<-selectedAttributes
    if((i %% 100)==0){
      print(i)
    }
    pb$tick()  #更新进度条#Update progress bar
  }
  #重新排序列 Rearrange columns
  vcf<-vcf %>% dplyr::select(c("CHROM","POS","ID","SYMBOL","REF","ALT","QUAL",
                               "FILTER","INFO","FORMAT","NA19238","NA19239","NA19240"))
  vcf<-vcf %>% dplyr::filter(FILTER=='PASS')
  saveRDS(vcf,'./Data/vcf.RDS')
  #meta数据 MetaData
  vcf_format_meta<-data.frame(
    abbrs=c('GT','GL','FT','RC','RCL','RCR','CN','DR','DV','RR','RV'),
    descriptions=c(
      "Genotype","Log10-scaled genotype likelihoods for RR,RA,AA genotypes",
      "Per-sample genotype filter","Raw high-quality read counts for the SV",
      "Raw high-quality read counts for the left control region",
      "Raw high-quality read counts for the right control region",
      "Read-depth based copy-number estimate for autosomal sites",
      "high-quality reference pairs","high-quality variant pairs",
      "high-quality reference junction reads","high-quality variant junction reads" 
    )
  )
  vcf_info_meta<-data.frame(
    abbrs=c("CIEND","CIPOS","CHR2","END","PE","MAPQ","SR","SRQ",
            "CONSENSUS","CE","CT","IMPRECISE","PRECISE",
            "SVTYPE","SVMETHOD","INSLEN","HOMLEN"),
    descriptions=
      c("PE confidence interval around END","PE confidence interval around POS",
        "Chromosome for END coordinate in case of a translocation","End position of the structural variant",
        "Paired-end support of the structural variant","Median mapping quality of paired-ends",
        "Split-read support","Split-read consensus alignment quality",
        "Split-read consensus sequence","Consensus sequence entropy",
        "Paired-end signature induced connection type","Imprecise structural variation",
        "Precise structural variation","Type of structural variant",
        "Type of approach used to detect SV","Predicted length of the insertion",
        "Predicted microhomology length using a max. edit distance of ")
  )
  vcf_meta<-rbind(vcf_info_meta,vcf_format_meta) %>%
    dplyr::arrange(abbrs)
}

#2.整理数据：编号、分割 Organize data: Indexing and splitting
##定义分割的函数 Define splitting functions
{
  SplitINFOFunction<-function(info_char,sample){
    library(stringr)
    splitedchars<-str_split_1(info_char,';')
    assign('Items',c());assign('Values',c())
    for(i in 1:length(splitedchars)){
      if(str_detect(splitedchars[i],'=')){
        Items<-c(Items,
                 str_sub(splitedchars[i],1,
                         (str_locate_all(splitedchars[i],'=')[[1]][1,1]-1)))
        Values<-c(Values,
                  str_sub(splitedchars[i],
                          (str_locate_all(splitedchars[i],'=')[[1]][1,1]+1),
                          (nchar(splitedchars[i]))))
      }else if(!str_detect(splitedchars[i],'=')){
        Items<-c(Items,splitedchars[i])
        Values<-c(Values,splitedchars[i])
      }
    }
    result<-data.frame(ITEM=Items,VALUE=Values)
    result$SAMPLE<-sample
    return(
      result
    )
  }
  SplitFORMATFunction<-function(format_char,values,sample){
    library(stringr)
    result<-data.frame(
      ITEM=str_split_1(format_char,':'),
      VALUE=str_split_1(values,':')
    )
    result$SAMPLE<-sample
    return(result)
  }
}
##建立新变量 Create new variables
VCF_SplitedDetailedLongerDataFrame<-data.frame(
  matrix(ncol=13,nrow=0)
)
colnames(VCF_SplitedDetailedLongerDataFrame)<-c(
  colnames(vcf)[1:8],'INFO_ITEM','INFO_VALUES','FORMAT_ITEM','SAMPLE_VALUE','SAMPLE'
)
pb2 <- progress_bar$new(total = nrow(vcf))#创建一个进度条对象 #Create a progress bar object
for(i in 1:nrow(vcf)){
  ITEM_Frame_i<-rbind(
    SplitINFOFunction(vcf$INFO[i],'NA19238'),
    SplitINFOFunction(vcf$INFO[i],'NA19239'),
    SplitINFOFunction(vcf$INFO[i],'NA19240')
  )
  FORMAT_Frame_i<-rbind(SplitFORMATFunction(vcf$FORMAT[i],
                                            vcf$NA19238[i],'NA19238'),
                        SplitFORMATFunction(vcf$FORMAT[i],
                                            vcf$NA19239[i],'NA19239'),
                        SplitFORMATFunction(vcf$FORMAT[i],
                                            vcf$NA19240[i],'NA19240'))
  Frame_Value_i<-rbind(ITEM_Frame_i,FORMAT_Frame_i)
  Frame_i<-vcf[i,(1:8)]
  while (nrow(Frame_i)<nrow(Frame_Value_i)) {
    Frame_i<-rbind(Frame_i,vcf[i,(1:8)])
  }
  VCF_SplitedDetailedLongerDataFrame<-rbind(
    VCF_SplitedDetailedLongerDataFrame,
    cbind(Frame_i,Frame_Value_i)
  )
  if((i %% 100)==0){
    print(i)
  }
  pb2$tick()  #更新进度条#Update progress bar
}

saveRDS(VCF_SplitedDetailedLongerDataFrame,
        './Data/VCF_SplitedDetailedLongerDataFrame.RDS')

#3.绘图 Plotting
colors_values<-c('#14517C','#2F7FC1','#E7EFFA','#96C37D','#F3D266','#D8383A',
                 '#F7E1ED','#F8F3F9','#C497B2','#A9B8C6','#8E8BFE','#63E398')
vcf_meta_plot<-vcf_meta %>%
  dplyr::filter(abbrs %in% c('CHR2','CN','MAPQ','DR','DV','GL','GT','INSLEN',
                             'RC','RCL','RCR','RR','RV','SVTYPE'))
for(i in 1:nrow(vcf_meta_plot)){
  selectedframe<-(VCF_SplitedDetailedLongerDataFrame %>%
                    dplyr::filter(ITEM==vcf_meta_plot$abbrs[i]))
  #染色体因子化确保排序 #Chromosome Character factorization ensures sorting
  if(vcf_meta_plot$abbrs[i]=='CHR2'){
    selectedframe$VALUE<-factor(
      selectedframe$VALUE,levels = unique(selectedframe$VALUE)
    )
  }
  #绘图 Plot
  if(length(which(is_numeric(selectedframe$VALUE)==F))>0){#数值型变量 Numerical variables
    plot<-ggplot(data=selectedframe,
                 mapping = aes(x=VALUE))+
      geom_histogram(stat = "count",alpha=0.7,color='black',
                     fill=sample(colors_values,1))+
      labs(x=vcf_meta_plot$abbrs[i],
           y=paste0(vcf_meta_plot$abbrs[i],'\n(',vcf_meta_plot$descriptions[i],')'))+
      theme_classic()+
      theme(axis.title = element_text(face='bold'),
            plot.title = element_text(face='bold',hjust=0.5),
            axis.text.x = element_text(angle = 90, vjust=0.5))
    ggsave(paste0('./Figs/A_PrimaryExploration/',i,'. ',vcf_meta_plot$abbrs[i],'.pdf'),width = 8,height = 6) 
  }else if(length(which(is_numeric(selectedframe$VALUE)==F))==0){#可全转为数字类型 Can be fully converted to numeric type
    plot<-ggplot(data=selectedframe,
                   mapping = aes(x=as.numeric(VALUE)))+
      geom_histogram(stat = "bin",alpha=0.7,color='black',
                     fill=sample(colors_values,1))+
      labs(x=vcf_meta_plot$abbrs[i],
           y=paste0(vcf_meta_plot$abbrs[i],'\n(',vcf_meta_plot$descriptions[i],')'))+
      theme_classic()+
      theme(axis.title = element_text(face='bold'),
            plot.title = element_text(face='bold',hjust=0.5),
            axis.text.x = element_text(angle = 90, vjust=0.5))
    ggsave(paste0('./Figs/A_PrimaryExploration/',i,'. ',vcf_meta_plot$abbrs[i],'.pdf'),width = 8,height = 6)
  }
}

#B.确定父母样本 Identifying Parental Samples
#目标：确认孩子对应的样本
#Goal: Confirm samples corresponding to child
#通过每一个变异来计算可能的情况
#Calculate possible scenarios through each mutation
{
  #1.建立函数基于表型寻找可能的小孩
  #1.Define function to identify potential child based on genotype
  FindWhoIsChild<-function(GTs){
    #列举可能的表型匹配 List possible genotype matches
    {
      possibleGTs<-data.frame(
        Parent1=c('0/0','0/0','0/0','0/0','0/1','0/1','0/1','0/1','0/1',
                  '0/1','0/1','1/1','1/1','1/1','1/1'),
        Parent2=c('0/0','0/1','0/1','1/1','0/0','0/0','0/1','0/1','0/1',
                  '1/1','1/1','0/0','0/1','0/1','1/1'),
        Child =c('0/0','0/0','0/1','0/1','0/0','0/1','0/0','0/1','1/1',
                 '0/1','1/1','0/1','1/1','0/1','1/1')
      )
    }
    #建立假设检验结果变量 Initialize hypothesis testing result variable
    assign('AsumeCheck',c())
    #建立循环逐个查询 Loop through each variant for validation
    for(i in 1:length(GTs)){
      #假设选中的为小孩的表型 #Assuming the selected phenotype is a child
      AsumedChildGT<-GTs[i]
      #对应的假设中父母的表型 The phenotype of the parents in the corresponding hypothesis
      AsumedParentsGT<-GTs[-i]
      #筛选出对应可能的父母表型 Screen out corresponding possible parental phenotypes
      AsumedPossibleParentsGTs<-possibleGTs %>%
        dplyr::filter(Child==AsumedChildGT) %>%
        dplyr::select(c('Parent1','Parent2'))
      #假设父母的第一个表型假设为家长1的可能情况数
      ##Assuming the first phenotype of the parents is assumed to be the number of possible scenarios for parent 1
      PossibleSituationCount1<-
        length(which(AsumedPossibleParentsGTs$Parent1==AsumedParentsGT[1] &
                       AsumedPossibleParentsGTs$Parent2==AsumedParentsGT[2]))
      #假设父母的第二个表型假设为家长1的可能情况数
      ##Assuming the second phenotype of parents is the number of possible scenarios for parent 1
      PossibleSituationCount2<-
        length(which(AsumedPossibleParentsGTs$Parent1==AsumedParentsGT[2] &
                       AsumedPossibleParentsGTs$Parent2==AsumedParentsGT[1]))
      #如果两个假设中有一个可能即为假设成立
      ##If one of the two hypotheses is possible, then the hypothesis is valid
      if((PossibleSituationCount1+PossibleSituationCount2)>0){
        AsumeCheck<-c(AsumeCheck,TRUE)
      }else{
        AsumeCheck<-c(AsumeCheck,FALSE)
      }
    }
    return(AsumeCheck)
  }
  #2.摘取VCF表格为新变量 Extract VCF table into a new variable
  FindChildFrame<-data.frame(matrix(ncol=6,nrow=0))
  colnames(FindChildFrame)<-c("CHROM","POS","ID","SYMBOL","SAMPLE","CHILD_CHECK_RESULT")
  #3.建立循环检验【【【【十分费时间】】】】
  #3.Perform loop validation (computationally intensive process)
  pb3 <- progress_bar$new(total = nrow(vcf))#创建一个进度条对象 #Create a progress bar object
  for(i in 1:nrow(vcf)){
    FindChildData_i<-VCF_SplitedDetailedLongerDataFrame %>%
      dplyr::filter(ID==vcf$ID[i]) %>%
      dplyr::filter(ITEM=='GT')
    if(length(which(FindChildData_i$VALUE=='./.'))>0){
      ChildCheckResults_i<-rep(NA,3)
    }else{
      ChildCheckResults_i<-
        FindWhoIsChild(FindChildData_i$VALUE)
    }
    FindChildFrame_i<-vcf[i,1:4]
    while (nrow(FindChildFrame_i)<nrow(FindChildData_i)) {
      FindChildFrame_i<-rbind(FindChildFrame_i,vcf[i,(1:4)])
    }
    FindChildFrame_i$SAMPLE<-FindChildData_i$SAMPLE
    FindChildFrame_i$CHILD_CHECK_RESULT<-ChildCheckResults_i
    FindChildFrame<-rbind(FindChildFrame,FindChildFrame_i)
    pb3$tick()  #更新进度条#Update progress bar
    if((i %% 1000)==0){
      print(i)
    }
  }
  #统计结果Statistical results
  FindChildFrame_Seq<-data.frame(
    SAMPLE=c('NA19238','NA19239','NA19240'),
    Freq=c(
      length(which(FindChildFrame$SAMPLE=='NA19238' &
                     FindChildFrame$CHILD_CHECK_RESULT==T)),
      # [1] 14089
      length(which(FindChildFrame$SAMPLE=='NA19239' &
                     FindChildFrame$CHILD_CHECK_RESULT==T)),
      # [1] 14291
      length(which(FindChildFrame$SAMPLE=='NA19240' &
                     FindChildFrame$CHILD_CHECK_RESULT==T))
      # [1] 17107
    )
  )
  # 确定NA19240为小孩 Confirm NA19240 is the child
  saveRDS(FindChildFrame,'./Figs/B_FindChild/FindChildFrame.RDS')
  
  #绘图 Plot
  library(ggbreak)
  ggplot()+
    geom_col(data=FindChildFrame_Seq,
             mapping=aes(x=SAMPLE,y=Freq),
             fill='skyblue',color='black',
             width = 0.8)+
    geom_text(data=FindChildFrame_Seq,
              mapping=aes(x=SAMPLE,y=Freq+300,
                          label=Freq))+
    labs(x='Sample',
         y='Frequence')+
    theme_classic()+
    scale_y_break(breaks = c(1000,13000),space = 1,scales = "fixed",
                  expand = T)+
    theme(axis.title = element_text(face = 'bold'))+
    coord_flip()
  ggsave('./Figs/B_FindChild/1.Sample_Freq.pdf',
         width = 8,height = 6)
}

#C.确定感兴趣的基因 Identify genes of interest
##1.基因内出现突变的频数 Frequency of variants within genes
{
  symbol_seq<-table(vcf$SYMBOL) %>%
    t() %>%
    data.frame() %>%
    dplyr::arrange(desc(as.numeric(Freq))) %>%
    dplyr::select(c('Var2','Freq')) %>%
    rename('Var2'='SYMBOL',
           'Freq'='Freq')
  ##因子化确保顺序 Factorize to ensure ordering
  symbol_seq$SYMBOL<-factor(symbol_seq$SYMBOL,levels=symbol_seq$SYMBOL)
  ##绘制频数的频数图 Plot frequency distribution of variant counts
  symbol_seq_seq<-table(symbol_seq$Freq) %>%
    t() %>%
    data.frame() %>%
    dplyr::arrange(desc(as.numeric(Freq))) %>%
    dplyr::select(c('Var2','Freq')) %>%
    clusterProfiler::rename('Freq1'='Var2',
           'Freq2'='Freq')
  symbol_seq_se_plot<-
    ggplot()+
    geom_col(data=symbol_seq_seq,
             mapping=aes(x=Freq1,y=Freq2),
             fill='pink',color='black',
             width = 0.8)+
    geom_text(data=symbol_seq_seq,
              mapping=aes(x=Freq1,y=Freq2+200,
                          label=Freq2))+
    labs(x='Variant Number of Genes',
         y='Frequence')+
    theme_classic()+
    theme(axis.title = element_text(face = 'bold'))+
    coord_flip()
  ggsave('./Figs/C_InterestedGenes/1.Symbol_seq_seq.pdf',
         symbol_seq_se_plot,
         width = 8,height = 6)
}
##2.频数Top18的频数图 Frequency plot for top 18 genes
{
  symbol_seq_top18<-symbol_seq[1:18,]
  symbol_seq_top18_plot<-
    ggplot()+
    geom_col(data=symbol_seq_top18,
             mapping=aes(x=SYMBOL,y=Freq),
             fill='#40E0D0',color='black',
             width = 0.8)+
    geom_text(data=symbol_seq_top18,
              mapping=aes(x=SYMBOL,y=Freq+1,
                          label=Freq))+
    labs(x='Genes',
         y='Frequence')+
    theme_classic()+
    theme(axis.title = element_text(face = 'bold'),
          axis.text.x = element_text(angle=90,vjust=0.5))
  ggsave('./Figs/C_InterestedGenes/2.Symbol_seq_top18_freq.pdf',
         symbol_seq_top18_plot,
         width = 6,height = 6)
}
##3.频数Top18的染色体定位 Chromosomal localization of top 18 genes
{
  symbol_seq_top18_genes<-symbol_seq$SYMBOL[1:18]
  #参数设置 Parameter settings
  RCircos.Set.Core.Components(UCSC.HG38.Human.CytoBandIdeogram,#这是上面load的基因组文件 #This is the genome file that was previously loaded
                              chr.exclude<- NULL, #无排除的染色体 No Chromosomes exclused
                              tracks.inside=10, #染色体圆圈内部一共要画10个圆圈 #There are a total of 10 circles to be drawn inside the chromosome circle
                              tracks.outside=0)#在外部画0个圆圈 #Draw 0 circles on the outside
  settings<-RCircos.Get.Plot.Parameters()
  settings$text.size<-1
  RCircos.Reset.Plot.Parameters(settings)
  
  #骨架 Frame
  RCircos.Set.Plot.Area()#建立一个画板 #Establish a drawing board
  RCircos.Chromosome.Ideogram.Plot()#在当前的画板上画基因组的圆圈骨架 #Draw the circular skeleton of the genome on the current canvas
  
  #候选基因位置信息 Location information for candidate genes
  #读取基因组注释文件，时间略长 Read genome annotation file (time-consuming)
  Hs_GRCh38.p14<-import('./Data/Homo_sapiens.GRCh38.114.chr.gtf.gz')
  #转换成数据框 Transfer to DataFrame
  Hs_GRCh38.p14_frame<-data.frame(Hs_GRCh38.p14)
  #提取候选基因位置信息 #Extract candidate gene location information
  candg_location<-Hs_GRCh38.p14_frame %>% 
    filter(type=='gene') %>% #保留的是完整基因的位置 #Preserving the location of complete genes
    filter(gene_name %in% symbol_seq_top18_genes) %>% #筛选出候选基因 #Screen out candidate genes
    dplyr::select('seqnames','start','end','gene_name') %>% #保留有用的列 #Keep useful columns
    rename('seqnames'='Chromosome',
           'start'='chromStart',
           'end'='chromEnd',
           'gene_name'='Gene')#重命名列 Rename Columns
  candg_location$Chromosome<-paste0('chr',candg_location$Chromosome)
  name.col <- 4 
  side <- "in" #画在基因组骨架的内侧 Inner Frame
  track.num <- 1 
  RCircos.Gene.Connector.Plot(candg_location,
                              + track.num, side)#画connector（连接基因名称和基因组位置）
  track.num <- 2
  RCircos.Gene.Name.Plot(candg_location,
                         + name.col,track.num, side)#加基因名称 Add symbol
  geneChrLocation<-recordPlot()
  pdf('./Figs/C_InterestedGenes/3.Top18GeneChrLocation.pdf',width=8,height=8)
  geneChrLocation
  dev.off()
}
##4.频数Top18的富集 Enrichment for top 18 genes
{
  setwd('./Figs/C_InterestedGenes/')
  folder_name<-"4_Enrichment";if (!dir.exists(folder_name)) {dir.create(folder_name)}
  setwd(folder_name)
  # 加载R包
  suppressPackageStartupMessages(library(clusterProfiler))
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  suppressPackageStartupMessages(library(enrichplot))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(GOplot))
  suppressPackageStartupMessages(library(cowplot))
  suppressPackageStartupMessages(library(patchwork))
  suppressPackageStartupMessages(library(tidyverse))
  library(export)
  library(circlize)
  library(grid)
  library(graphics)
  library(ComplexHeatmap)
  select=dplyr::select
  #候选基因 Candidate Genes
  candg <- symbol_seq_top18_genes
  ## id转换 Convert Gene ID
  symbol2entrezid <- bitr(geneID = candg,
                          fromType = 'SYMBOL',
                          toType = 'ENTREZID',
                          OrgDb = 'org.Hs.eg.db')
  # GO富集分析 GO Enrichment####---------------------
  pvalueCutoff <- 1
  qvalueCutoff <- 1
  
  ego <- enrichGO(gene = as.numeric(symbol2entrezid$ENTREZID),
                  keyType = "ENTREZID",
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff = pvalueCutoff, 
                  qvalueCutoff = qvalueCutoff,
                  ont="ALL",
                  readable =T)
  saveRDS(ego,'ego.rds')
  as.data.frame(ego) %>% filter(pvalue <= 0.05) %>% 
    group_by(ONTOLOGY) %>% dplyr::count()
  # # A tibble: 3 × 2
  # # Groups:   ONTOLOGY [3]
  # ONTOLOGY     n
  # <chr>    <int>
  #   1 BP         368
  # 2 CC          32
  # 3 MF          56
  write_csv(as.data.frame(ego) %>% filter(pvalue <= 0.05),"1.Rich_GO_enrich_sig.csv")
  go.df<-as.data.frame(ego) %>% filter(pvalue <= 0.05)
  go.df <- go.df %>% group_by(ONTOLOGY) %>% slice_head(n=5)
  # 使画出的GO term的顺序与输入一致
  # Make sure the order in plot is same with frame
  go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))
  # 绘图 plot
  GO_Plot<-
    ggplot(data = go.df)+ # 绘图使用的数据
    geom_point(aes(x = Description, y=(-log10(pvalue)), 
                   size = Count,color = ONTOLOGY))+
    scale_color_manual(values = c("#0000CD","orange","#43CD80"))+
    coord_flip()+# 横纵坐标反转
    theme_bw()+ #去除背景色
    scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
    labs(x = "GO terms",y = paste0('-log10(P)'))+ # 设置坐标轴标题及标题
    ggtitle('GO Enrichment')+
    ggplot2::theme(axis.title = element_text(size = 13), # 坐标轴标题大小
                   axis.text = element_text(size = 11), # 坐标轴标签大小
                   plot.title = element_text(size = 10,hjust = 0.5,face = "bold"), # 标题设置
                   legend.title = element_text(size = 10), # 图例标题大小
                   legend.text = element_text(size = 10) # 图例标签大小
    )
  ggsave('2.GO_bubble.pdf',GO_Plot,width = 10,height=6)
  
  # KEGG富集分析 KEGG enrichment####---------------------
  ekegg <- enrichKEGG(gene = symbol2entrezid$ENTREZID ,
                      keyType = "kegg",
                      organism = "hsa",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)
  saveRDS(ekegg,'ekegg.rds')
  
  ekegg2 <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  write_csv(as.data.frame(ekegg2),"3.KEGG_enrich.csv")
  
  kegg.df <- ekegg2[order(ekegg2$pvalue),]
  kegg.df$Y_Axis_Value<-kegg.df[,which(colnames(kegg.df)=='pvalue')]
  #计算数值型GeneRatio Colculate GeneRatio
  kegg.df$GeneRatio_Number<-NA
  for (i in 1:nrow(kegg.df)) {
    kegg.df$GeneRatio_Number[i]<-
      as.numeric(
        str_sub(kegg.df$GeneRatio[i],1,str_locate_all(kegg.df$GeneRatio[i],'/')[[1]][1,1]-1)
      )/as.numeric(
        str_sub(kegg.df$GeneRatio[i],
                str_locate_all(kegg.df$GeneRatio[i],'/')[[1]][1,1]+1,
                nchar(kegg.df$GeneRatio[i]))
      )
  }
  
  kegg_df_top<-kegg.df %>% 
    slice_head(n=10)
  # 使画出的kegg term的顺序与输入一致
  kegg_df_top$Description <- factor(kegg_df_top$Description,levels = rev(kegg_df_top$Description))
  # 绘图
  KEGG_Plot<-
    ggplot(data = kegg_df_top)+ # 绘图使用的数据
    geom_point(aes(x = Description, y=GeneRatio_Number, size = Count,color = (-log10(pvalue))))+
    scale_color_gradient(low="blue",high="red")+
    coord_flip()+# 横纵坐标反转
    theme_bw()+ #去除背景色
    scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
    labs(x = "KEGG terms",y = "GeneRatio", 
         color = paste0('-log10(P)'))+ # 设置坐标轴标题及标题
    ggtitle('KEGG Enrichment')+
    ggplot2::theme(axis.title = element_text(size = 13), # 坐标轴标题大小
                   axis.text = element_text(size = 11), # 坐标轴标签大小
                   plot.title = element_text(size = 10,hjust = 0.5,face = "bold"), # 标题设置
                   legend.title = element_text(size = 10), # 图例标题大小
                   legend.text = element_text(size = 10) # 图例标签大小
    )
  ggsave('4.KEGG_bubble.pdf',KEGG_Plot,width = 10,height=6)
}
















