##############################函数描述##############################
# “Facets.SNP.Pileup”通过R函数传参调用SNP Pileup统计来自同一个体的正常样本和肿瘤样本在SNP位点处的[参考、替代、错误、缺失]的read数
# "Facets.CNV.Calling"通过Facets包通过SNP Pileup结果推断等位特异拷贝数
####################################################################


##' @description 通过R函数传参调用SNP Pileup统计来自同一个体的正常样本和肿瘤样本在SNP位点处的[参考、替代、错误、缺失]的read数
##' @param Tumor.Bam character 来自同一个体的肿瘤样本对应的bam文件
##' @param Normal.Bam character 来自同一个体的正常样本对应的bam文件
##' @param Common.Vcf character 常见的多态SNP的VCF文件(一个很好的来源是dbSNP的common_all.vcf.gz)
##' @param Output.Prefix character 结果文件(csv格式)前缀[可携带路径]; 默认为当前工作目录下的"Facets.SNP-Pileup"
##' @param Sort.Operation character 设置对Tumor.Bam、Normal.Bam、Common.Vcf进行排序操作的方式, 可选("Auto", "Sort", "None"); 默认"Auto"即交由程序自行判断并排序
##' @param Show.Progress logical 设置是否显示SNP Pileup的进度条(需要先统计SNP数量，增加程序的运行时间); 默认为FALSE
##' @param Compress.Output logical 设置是否对对结果文件进行压缩(csv.gz格式); 默认TRUE
##' @param Skip.Anomalous logical 设置是否跳过异常的异常读取对; 默认TRUE
##' @param Check.Overlaps logical 设置是否启用读取对重叠检测; 默认TRUE
##' @param Max.Depth numeric 设置最大深度, 即每个位置允许的最大read数; 默认0
##' @param Pseudo.Snps numeric 设置虚拟SNP间隔, 即若每Pseudo.Snps个碱基位置上没有SNP则插入一个标记该区间总read数的空记录; 默认0
##' @param Min.Map.Quality numeric 设置比对质量的最小阈值(针对MAPQ信息); 默认0
##' @param Min.Base.Quality numeric 设置序列质量的最小阈值(针对QUAL信息); 默认0
##' @param Min.Read.Counts numeric[] 设置输出位点在正常样本和肿瘤样本中最小read数; 默认c(0, 0)
Facets.SNP.Pileup <- function(Tumor.Bam, Normal.Bam, Common.Vcf, 
                              Output.Prefix = NULL, Sort.Operation  = c("Auto", "Sort", "None"), 
                              Show.Progress = FALSE, Compress.Output = TRUE, Skip.Anomalous = TRUE, Check.Overlaps = TRUE, 
                              Max.Depth = 4000, Pseudo.Snps = 0, Min.Map.Quality = 0, Min.Base.Quality = 0, Min.Read.Counts = c(0, 0), 
                              System.Pileup.Alias = "snp-pileup", System.Samtools.Alias = "samtools", System.Bedtools.Alias = "bedtools"){
  
  System.Pileup.Alias <- as.character(System.Pileup.Alias)
  if(length(System.Pileup.Alias) == 1){
    # 判断System.Pileup.Alias在系统中是否存在
    if(nchar(Sys.which(System.Pileup.Alias)) > 0){
      
      # 文件参数的判断
      Tumor.Bam <- as.character(Tumor.Bam)
      Normal.Bam <- as.character(Normal.Bam)
      Common.Vcf <- as.character(Common.Vcf)
      if((length(Tumor.Bam) == 1 && file.exists(Tumor.Bam)) && 
         (length(Normal.Bam) == 1 && file.exists(Normal.Bam)) && 
         (length(Common.Vcf) == 1 && file.exists(Common.Vcf))){
        Tumor.Bam <- normalizePath(Tumor.Bam, winslash = "/", mustWork = TRUE)
        Normal.Bam <- normalizePath(Normal.Bam, winslash = "/", mustWork = TRUE)
        Common.Vcf <- normalizePath(Common.Vcf, winslash = "/", mustWork = TRUE)
      }else{
        stop("'Tumor.Bam'、'Normal.Bam'、'Common.Vcf'应为单一的文件路径 ...")
      }
      
      # 初始化SNP Pileup指令
      SNP.Pileup.Command <- sprintf("echo \"=>=>=>正在进行SNP Pileup ...\"\n\"%s\"", System.Pileup.Alias)
      
      # 配置Show.Progress[--progress / -p]
      Show.Progress <- as.logical(Show.Progress)
      if(length(Show.Progress) == 1){
        if(Show.Progress){ SNP.Pileup.Command <- sprintf("%s --progress", SNP.Pileup.Command) } 
      }else{
        stop("'Show.Progress'应为单一的logical值 ...")
      }
      
      # 配置Compress.Output[--gzip / -g]
      Compress.Output <- as.logical(Compress.Output)
      if(length(Compress.Output) == 1){
        if(Compress.Output){ SNP.Pileup.Command <- sprintf("%s --gzip", SNP.Pileup.Command) } 
      }else{
        stop("'Compress.Output'应为单一的logical值 ...")
      }
      
      # 配置Skip.Anomalous[--count-orphans / -A]
      Skip.Anomalous <- as.logical(Skip.Anomalous)
      if(length(Skip.Anomalous) == 1){
        if(! Skip.Anomalous){ SNP.Pileup.Command <- sprintf("%s --count-orphans", SNP.Pileup.Command) } 
      }else{
        stop("'Skip.Anomalous'应为单一的logical值 ...")
      }
      
      # 配置Check.Overlaps[--ignore-overlap / -x]
      Check.Overlaps <- as.logical(Check.Overlaps)
      if(length(Check.Overlaps) == 1){
        if(! Check.Overlaps){ SNP.Pileup.Command <- sprintf("%s --ignore-overlap", SNP.Pileup.Command) } 
      }else{
        stop("'Check.Overlaps'应为单一的logical值 ...")
      }
      
      # 配置Max.Depth[--max-depth / -d]
      Max.Depth <- as.numeric(Max.Depth)
      if(length(Max.Depth) == 1 && Max.Depth > 0 && Max.Depth %% 1 == 0){
        SNP.Pileup.Command <- sprintf("%s --max-depth %s", SNP.Pileup.Command, Max.Depth)
      }else{
        stop("'Max.Depth'应为单一的大于0的整型numeric值 ...")
      }
      
      # 配置Pseudo.Snps[--pseudo-snps / -P]
      Pseudo.Snps <- as.numeric(Pseudo.Snps)
      if(length(Pseudo.Snps) == 1 && Pseudo.Snps >= 0 && Pseudo.Snps %% 1 == 0){
        SNP.Pileup.Command <- sprintf("%s --pseudo-snps %s", SNP.Pileup.Command, Pseudo.Snps)
      }else{
        stop("'Pseudo.Snps'应为单一的大于等于0的整型numeric值 ...")
      }
      
      # 配置Min.Map.Quality[--min-map-quality / -q]
      Min.Map.Quality <- as.numeric(Min.Map.Quality)
      if(length(Min.Map.Quality) == 1 && Min.Map.Quality >= 0 && Min.Map.Quality %% 1 == 0){
        SNP.Pileup.Command <- sprintf("%s --min-map-quality %s", SNP.Pileup.Command, Min.Map.Quality)
      }else{
        stop("'Min.Map.Quality'应为单一的大于等于0的整型numeric值 ...")
      }
      
      # 配置Min.Base.Quality[--min-base-quality / -Q]
      Min.Base.Quality <- as.numeric(Min.Base.Quality)
      if(length(Min.Base.Quality) == 1 && Min.Base.Quality >= 0 && Min.Base.Quality %% 1 == 0){
        SNP.Pileup.Command <- sprintf("%s --min-base-quality %s", SNP.Pileup.Command, Min.Base.Quality)
      }else{
        stop("'Min.Base.Quality'应为单一的大于等于0的整型numeric值 ...")
      }
      
      # 配置Min.Read.Counts[--min-read-counts / -r]
      Min.Read.Counts <- as.numeric(Min.Read.Counts)
      if(length(Min.Read.Counts) == 2 && all(Min.Read.Counts >= 0) && all(Min.Read.Counts %% 1 == 0)){
        SNP.Pileup.Command <- sprintf("%s --min-read-counts %s,%s", SNP.Pileup.Command, Min.Read.Counts[1], Min.Read.Counts[2])
      }else{
        stop("'Min.Read.Counts'应为大于等于0的包含两个元素的整型numeric向量 ...")
      }
      
      # 配置Output.Prefix
      Output.Prefix <- as.character(Output.Prefix)
      if(length(Output.Prefix) == 0){
        Output.Prefix <- sprintf("%s/Facets.SNP-Pileup", getwd())
      }else if(length(Output.Prefix) == 1){
        Output.Prefix <- normalizePath(sub(pattern = "(\\\\*/*|/*\\\\*)$", replacement = "", x = trimws(Output.Prefix)), winslash = "/", mustWork = FALSE)
        dir.create(dirname(Output.Prefix), recursive = TRUE, showWarnings = FALSE)
      }else{
        stop("'Output.Prefix'应为Null或单一的character值 ...")
      }
      unlink((SNP.Pileup.Output <- sprintf("%s.%s", Output.Prefix, ifelse(Compress.Output, "csv.gz", "csv"))), force = TRUE)
      
      # 对bam文件与vcf文件进行排序处理
      Sort.Operation <- match.arg(Sort.Operation)
      if(Sort.Operation == "None"){
        Tumor.Sorted.Bam <- Tumor.Bam; Normal.Sorted.Bam <- Normal.Bam; Common.Sorted.Vcf <- Common.Vcf
      }else{
        System.Samtools.Alias <- as.character(System.Samtools.Alias)
        System.Bedtools.Alias <- as.character(System.Bedtools.Alias)
        if(length(System.Samtools.Alias) == 1 && length(System.Bedtools.Alias) == 1 ){
          # 判断System.Samtools.Alias与System.Bedtools.Alias在系统中是否存在
          if(nchar(Sys.which(System.Samtools.Alias)) > 0){
            if(nchar(Sys.which(System.Bedtools.Alias)) > 0){
              Tumor.Bam.SeqName.Order <- system(sprintf("\"%s\" view -H \"%s\" | grep ^@SQ | cut -f 2 | tr -d \"SN:\" | uniq", System.Samtools.Alias, Tumor.Bam), intern = TRUE)
              Normal.Bam.SeqName.Order <- system(sprintf("\"%s\" view -H \"%s\" | grep ^@SQ | cut -f 2 | tr -d \"SN:\" | uniq", System.Samtools.Alias, Normal.Bam), intern = TRUE)
              if(all(Tumor.Bam.SeqName.Order == Normal.Bam.SeqName.Order)){
                Common.Vcf.SeqName.Order <- system(sprintf("awk '{print $1}' \"%s\" | grep -v \"^#\" | uniq", Common.Vcf), intern = TRUE)
                if(any(grepl("^chr([1-9]|1[0-9]|2[0-2]|[MXY]|MT)$", Common.Vcf.SeqName.Order))){
                  Bam.SeqName.Order <- intersect(Tumor.Bam.SeqName.Order, Normal.Bam.SeqName.Order)
                }else{
                  Bam.SeqName.Order <- gsub(pattern = "^chr", "", intersect(Tumor.Bam.SeqName.Order, Normal.Bam.SeqName.Order))
                }
                switch(Sort.Operation, 
                       Auto = {
                         if(system(sprintf("\"%s\" view -H \"%s\" | grep ^@HD.*SO:coordinate.*$", System.Samtools.Alias, Tumor.Bam), ignore.stdout = TRUE) == 0){
                           Tumor.Sorted.Bam <- Tumor.Bam
                         }else{
                           Tumor.Sorted.Bam <- sprintf("%s/(Coordinate.Sorted)%s", dirname(Tumor.Bam), basename(Tumor.Bam))
                           File.Order.Command <- sprintf("echo \"=>=>=>正在对'%s'进行排序操作 ...\"\n\"%s\" sort -o \"%s\" \"%s\"\n", Tumor.Bam, System.Samtools.Alias, Tumor.Sorted.Bam, Tumor.Bam)
                         }
                         if(system(sprintf("\"%s\" view -H \"%s\" | grep ^@HD.*SO:coordinate.*$", System.Samtools.Alias, Normal.Bam), ignore.stdout = TRUE) == 0){
                           Normal.Sorted.Bam <- Normal.Bam
                         }else{
                           Normal.Sorted.Bam <- sprintf("%s/(Coordinate.Sorted)%s", dirname(Normal.Bam), basename(Normal.Bam))
                           File.Order.Command <- sprintf("%secho \"=>=>=>正在对'%s'进行排序操作 ...\"\n\"%s\" sort -o \"%s\" \"%s\"\n", ifelse(exists("File.Order.Command"), File.Order.Command, ""), Normal.Bam, System.Samtools.Alias, Normal.Sorted.Bam, Normal.Bam)
                         }
                         if(all(intersect(Bam.SeqName.Order, Common.Vcf.SeqName.Order) == intersect(Common.Vcf.SeqName.Order, Bam.SeqName.Order))){
                           Common.Sorted.Vcf <- Common.Vcf
                         }else{
                           Ref.Order <- c(Bam.SeqName.Order, setdiff(Common.Vcf.SeqName.Order, Bam.SeqName.Order))
                           Ref.Order.File <- sprintf("%s/Ref.Order.txt", getwd())
                           write(Ref.Order, file = Ref.Order.File, sep = "\n")
                           on.exit({unlink(Ref.Order.File, force = TRUE)}, add = TRUE)
                           Common.Sorted.Vcf <- sprintf("%s/(Sorted.Ref.Bam)%s", dirname(Common.Vcf), basename(Common.Vcf))
                           File.Order.Command <- sprintf("%secho \"=>=>=>正在对'%s'进行排序操作 ...\"\n\"%s\" sort -header -faidx \"%s\" -i \"%s\" > \"%s\"\n", ifelse(exists("File.Order.Command"), File.Order.Command, ""), Common.Vcf, System.Bedtools.Alias, Ref.Order.File, Common.Vcf, Common.Sorted.Vcf)
                         }
                       }, 
                       Sort = {
                         message(sprintf("将对'%s'、'%s'、'%s'进行排序操作 ...", Tumor.Bam, Normal.Bam, Common.Vcf))
                         Tumor.Sorted.Bam <- sprintf("%s/(Coordinate.Sorted)%s", dirname(Tumor.Bam), basename(Tumor.Bam))
                         File.Order.Command <- sprintf("echo \"=>=>=>正在对'%s'进行排序操作 ...\"\n\"%s\" sort -o \"%s\" \"%s\"\n", Tumor.Bam, System.Samtools.Alias, Tumor.Sorted.Bam, Tumor.Bam)
                         Normal.Sorted.Bam <- sprintf("%s/(Coordinate.Sorted)%s", dirname(Normal.Bam), basename(Normal.Bam))
                         File.Order.Command <- sprintf("%secho \"=>=>=>正在对'%s'进行排序操作 ...\"\n\"%s\" sort -o \"%s\" \"%s\"\n", File.Order.Command, Normal.Bam, System.Samtools.Alias, Normal.Sorted.Bam, Normal.Bam)
                         Ref.Order <- c(Bam.SeqName.Order, setdiff(Common.Vcf.SeqName.Order, Bam.SeqName.Order))
                         Ref.Order.File <- sprintf("%s/Ref.Order.txt", getwd())
                         write(Ref.Order, file = Ref.Order.File, sep = "\n")
                         on.exit({unlink(Ref.Order.File, force = TRUE)}, add = TRUE)
                         Common.Sorted.Vcf <- sprintf("%s/(Sorted.Ref.Bam)%s", dirname(Common.Vcf), basename(Common.Vcf))
                         File.Order.Command <- sprintf("%secho \"=>=>=>正在对'%s'进行排序操作 ...\"\n\"%s\" sort -header -faidx \"%s\" -i \"%s\" > \"%s\"\n", File.Order.Command, Common.Vcf, System.Bedtools.Alias, Ref.Order.File, Common.Vcf, Common.Sorted.Vcf)
                       })
              }else{
                stop(sprintf("'%s'与'%s'的参考基因组信息(@SQ)不一致 ...", Tumor.Bam, Normal.Bam))
              }
            }else{
              stop(sprintf("非系统的可执行命令'%s' ...", System.Samtools.Alias))
            }
          }else{
            stop(sprintf("非系统的可执行命令'%s ...", System.Bedtools.Alias))
          }
        }else{
          stop("'System.Samtools.Alias'与'System.Bedtools.Alias'应为单一的character值 ...")
        }
      }
      
      # 完成最终指令的拼接
      SNP.Pileup.Command <- sprintf("%s%s \"%s\" \"%s.csv\" \"%s\" \"%s\"", ifelse(exists("File.Order.Command"), File.Order.Command, ""), SNP.Pileup.Command, Common.Sorted.Vcf, Output.Prefix, Normal.Sorted.Bam, Tumor.Sorted.Bam)
      
      # 根据操作系统环境设置脚本内容
      SNP.Pileup.Command <- sprintf(ifelse(Sys.info()["sysname"] == "Windows", "@echo off\n%s", "#!/bin/sh\n%s"),  SNP.Pileup.Command)
      # 根据操作系统环境设置脚本文件
      SNP.Pileup.Command.File <- sprintf(ifelse(Sys.info()["sysname"] == "Windows", "%s/SNP.Pileup.Command.bat", "%s/SNP.Pileup.Command.sh"), getwd())
      # 将指令写入对应系统的脚本文件
      write(SNP.Pileup.Command, SNP.Pileup.Command.File)
      # 赋予脚本文件读写以及可执行权限
      Sys.chmod(SNP.Pileup.Command.File)
      
      # 运行SNP Pileup
      tryCatch(
        {
          message("<<====== RUNNING MESSAGE ======>>")
          # 执行脚本文件SNP.Pileup.Command.File
          SNP.Pileup.Command.Run <- system(sprintf("\"%s\"", SNP.Pileup.Command.File))
        },
        error = function(e){ # 抛出错误信息
          message(sprintf("<<====== ERROR MESSAGE ======>>\n%s", e))
        },
        warning = function(w){ # 抛出警告信息
          message(sprintf("Warning: %s ...", trimws(gsub(".*\\)\\:", "", w))))
        },
        finally = {
          # 如果脚本文件顺利执行，则SNP.Pileup.Command.Run返回的状态信息为0
          if(exists("SNP.Pileup.Command.Run") && SNP.Pileup.Command.Run == 0){
            message(sprintf("<<===== SUCCESS MESSAGE =====>>\nSNP.Pileup.Command执行成功, 结果已输出至文件'%s' ...", normalizePath(SNP.Pileup.Output, winslash = "/", mustWork = TRUE)))
            unlink(SNP.Pileup.Command.File, force = TRUE)
          }else{
            message(sprintf("<<====== ERROR MESSAGE ======>>\nSNP.Pileup.Command执行过程中发生了错误, 请通过查看'RUNNING MESSAGE'中的信息或通过控制台运行脚本文件'%s'来查看具体错误 ...", SNP.Pileup.Command.File))
          }
        }
      )
      
    }else{
      stop(sprintf("非系统的可执行命令'%s' ...", System.Pileup.Alias))
    }
  }else{
    stop("'System.Pileup.Alias'应为单一的character值 ...")
  }
}


##' @description 通过Facets包通过SNP Pileup结果推断等位特异拷贝数
##' @param SNP.Pileup.Input character SNP Pileup的结果文件, 可通过函数Facets.SNP.Pileup获得
##' @param Err.Thresh numeric 设置SNP位点发生错误替代的read计数阈值, SNP Pileup结果中超过该值的记录将被去除; 默认Inf
##' @param Del.Thresh numeric 设置SNP位点发生缺失的read计数阈值, SNP Pileup结果中超过该值的记录将被去除; 默认Inf
##' @param Min.Depth numeric 设置最小深度, 数据预处理时总read计数小于该值的记录将被去除; 默认35
##' @param Max.Depth numeric 设置最大深度, 数据预处理时总read计数大于该值的记录将被去除; 默认1000
##' @param SNP.Het.VAF numeric 设置SNP位点被判定为杂合位点的阈值, 正常样本VAF介于(SNP.Het.VAF, 1 - SNP.Het.VAF)之间的SNP被判定为杂合; 默认0.25
##' @param Bin.Size numeric 设置窗口大小, 每个窗口将随机选择一个SNP位点(优先选择杂合位点)进行后续分析(基因组中的SNP不是均匀分布的, 使用所有位点将导致数据中的序列相关)
##' @param Segment.Min.Het numeric 设置分析次拷贝数时基因组片段至少应该具有的杂合性SNP位点的数量(当一个片段的杂 SNP少于Segment.Min.Het时，可能无法可靠地估计次要拷贝数，因此将返回); 默认15
##' @param EM.Max.Iter numeric 设置期望最大化算法最大的迭代次数; 默认10
##' @param EM.Con.Thresh numeric 设置在EM.Max.Iter内达到终止条件的收敛阈值; 默认0.001
##' @param Genome.Assemblies character 设置基因组版本号, 可选c("hg18", "hg19", "hg38", "mm9", "mm10", "udef")
##' @param Show.Plot logica 设置是否绘制基因组图示(包括LogR、LogOR、CNV)
##' @return list 包含肿瘤纯度、倍性、位点的估计信息[SeqName、Position、LogR、LogOR]以及片段的估计信息[SeqName、Position.Start、Position.End、LogR.Mean、LogOR.Mean.Abs、CN.Total、CN.Minor]
Facets.CNV.Calling <- function(SNP.Pileup.Input, 
                               Genome.Assemblies = c("hg18", "hg19", "hg38", "mm9", "mm10", "udef"), Err.Thresh = Inf, Del.Thresh = Inf, 
                               Min.Depth = 35, Max.Depth = 1000, SNP.Het.VAF = 0.25, Bin.Size = 250, Segment.Min.Het = 15, EM.Max.Iter = 10, EM.Con.Thresh = 0.001, Show.Plot = TRUE){
  library(facets)
  SNP.Pileup.Input <- as.character(SNP.Pileup.Input)
  if(length(SNP.Pileup.Input) == 1 && file.exists(SNP.Pileup.Input)){
    
    ############
    ## 1.参数判断
    ############
    Err.Thresh <- as.numeric(Err.Thresh)
    if(length(Err.Thresh) != 1 || ! (is.infinite(Err.Thresh) || (Err.Thresh >= 0 && Err.Thresh %% 1 == 0))){
      stop("'Err.Thresh'应为Inf或单一的大于等于0的整型numeric值 ...")
    }
    
    Del.Thresh <- as.numeric(Del.Thresh)
    if(length(Del.Thresh) != 1 || ! (is.infinite(Del.Thresh) || (Del.Thresh >= 0 && Del.Thresh %% 1 == 0))){
      stop("'Del.Thresh'应为Inf或单一的大于等于0的整型numeric值 ...")
    }
    
    Min.Depth <- as.numeric(Min.Depth)
    if(length(Min.Depth) != 1 || Min.Depth < 0 || Min.Depth %% 1 != 0){
      stop("'Min.Depth'应为单一的大于等于0的整型numeric值 ...")
    }
    
    Max.Depth <- as.numeric(Max.Depth)
    if(length(Max.Depth) != 1 || ! (is.infinite(Max.Depth) || (Max.Depth >= 0 && Max.Depth %% 1 == 0))){
      stop("'Max.Depth'应为Inf或单一的大于等于0的整型numeric值 ...")
    }
    
    SNP.Het.VAF <- as.numeric(SNP.Het.VAF)
    if(length(SNP.Het.VAF) != 1 || SNP.Het.VAF < 0 || SNP.Het.VAF >= 0.5){
      stop("'SNP.Het.VAF'应为单一的介于[0,0.5)之间的numeric值 ...")
    }
    
    Bin.Size <- as.numeric(Bin.Size)
    if(length(Bin.Size) != 1 || Bin.Size < 0 || Bin.Size %% 1 != 0){
      stop("'Bin.Size'应为单一的大于等于0的整型numeric值 ...")
    }
    
    Segment.Min.Het <- as.numeric(Segment.Min.Het)
    if(length(Segment.Min.Het) != 1 || Segment.Min.Het <= 0 || Segment.Min.Het %% 1 != 0){
      stop("'EM.Con.Thresh'应为单一的大于0的整型numeric值 ...")
    }
    
    EM.Max.Iter <- as.numeric(EM.Max.Iter)
    if(length(EM.Max.Iter) != 1 || EM.Max.Iter < 0 || EM.Max.Iter %% 1 != 0){
      stop("'EM.Max.Iter'应为单一的大于等于0的整型numeric值 ...")
    }
    
    EM.Con.Thresh <- as.numeric(EM.Con.Thresh)
    if(length(EM.Con.Thresh) != 1 || EM.Con.Thresh <= 0){
      stop("'EM.Con.Thresh'应为单一的大于0的numeric值 ...")
    }
    
    Show.Plot <- as.logical(Show.Plot)
    if(length(Show.Plot) != 1){
      stop("'Show.Plot'应为单一的logical值 ...")
    }
    
    ############
    ## 2.数据预处理以及拷贝数的估计
    ############
    # 读取snp-pileup生成SNP位点read计数矩阵, 统计SNP位点在正常样本与肿瘤样本中的总read数以及匹配到参考等位的read数
    SNP.Read.Mtr <- readSnpMatrix(SNP.Pileup.Input, err.thresh = Err.Thresh, del.thresh = Del.Thresh)
    Chrom.Is.Numeric <- grepl("^\\d*$", SNP.Read.Mtr$Chromosome)
    Mtr.Chrom.Numeric <- SNP.Read.Mtr[Chrom.Is.Numeric, ]
    Mtr.Chrom.Character <- SNP.Read.Mtr[! Chrom.Is.Numeric, ]
    SNP.Read.Mtr <- rbind(Mtr.Chrom.Numeric[order(as.numeric(Mtr.Chrom.Numeric$Chromosome)), ], Mtr.Chrom.Character[order(Mtr.Chrom.Character$Chromosome), ])
    # 数据预处理(计算LogR与LogOR值, 将SNP位点片段化) 
    Mtr.Processes <- preProcSample(SNP.Read.Mtr, gbuild = match.arg(Genome.Assemblies), ndepth = Min.Depth, ndepthmax = Max.Depth, het.thresh = SNP.Het.VAF, snp.nbhd = Bin.Size)
    # 对片段聚类并估计初始的等位特异拷贝数 
    NV.Fit <- procSample(Mtr.Processes, cval = 150, min.nhet = Segment.Min.Het)
    # 基于期望最大化算法估计等位特异拷贝数和细胞分数以及肿瘤纯度和倍性
    EM.Fit <- emcncf(NV.Fit, min.nhet = Segment.Min.Het, maxiter = EM.Max.Iter, eps = EM.Con.Thresh)
    
    ############
    ## 3.绘图
    ############
    if(Show.Plot){
      plotSample(x = NV.Fit, emfit = EM.Fit)
    }
    
    ############
    ## 3.结果封装
    ############
    return(list(
      Purity = EM.Fit[["purity"]],
      Ploidy = EM.Fit[["ploidy"]],
      Point.Data = data.frame(SeqName = NV.Fit$jointseg$chrom, Position = NV.Fit$jointseg$maploc, LogR = NV.Fit$jointseg$cnlr, LogOR = NV.Fit$jointseg$valor),
      Segment.Data = data.frame(SeqName = EM.Fit$cncf$chrom,  Position.Start = EM.Fit$cncf$start, Position.End = EM.Fit$cncf$end, LogR.Mean = EM.Fit$cncf$cnlr.median, LogOR.Mean.Abs = sqrt(EM.Fit$cncf$mafR), CN.Total = EM.Fit$cncf$tcn.em, CN.Minor = EM.Fit$cncf$lcn.em)
    ))
    
  }else{
    stop("'SNP.Pileup.Input'应为单一的character值 ...")
  }
}
