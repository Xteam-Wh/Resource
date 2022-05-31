##############################函数描述##############################
# "ASCAT.CNV.Calling"结合SNP位点的BAF值与LogR值推断等位特异拷贝数
# “ASCAT.Extract.LogR.BAF”从高通量测序数据中计算SNP位点处的LogR值与BAF值
####################################################################



##' @description 结合SNP位点的BAF值与LogR值推断等位特异拷贝数
##' @param Tumor.BAF.File character 肿瘤样本的BAF文件
############' 要求制表符分隔且包含行名列名, 行名为SNP位点的标识符, 各列所包含的信息依次为[SeqName、Position、Sample1.BAF、Sample2.BAF、Sample3.BAF、...], 其中[SeqName]信息不能包含"chr"前缀, 列名整体前端要保留一个制表符
##' @param Tumor.LogR.File character 肿瘤样本的LogR文件
############' 要求制表符分隔且包含行名列名, 行名为SNP位点的标识符, 各列所包含的信息依次为[SeqName、Position、Sample1.LogR、Sample2.LogR、Sample3.LogR、...], 其中[SeqName]信息不能包含"chr"前缀, 列名整体前端要保留一个制表符
##' @param Normal.BAF.File character 正常样本的BAF文件; 默认NULL
############' 要求制表符分隔且包含行名列名, 行名为SNP位点的标识符, 各列所包含的信息依次为[SeqName、Position、Sample1.BAF、Sample2.BAF、Sample3.BAF、...], 其中[SeqName]信息不能包含"chr"前缀, 列名整体前端要保留一个制表符
##' @param Normal.LogR.File character 正常样本的LogR文件; 默认NULL
############' 要求制表符分隔且包含行名列名, 行名为SNP位点的标识符, 各列所包含的信息依次为[SeqName、Position、Sample1.LogR、Sample2.LogR、Sample3.LogR、...], 其中[SeqName]信息不能包含"chr"前缀, 列名整体前端要保留一个制表符
##' @param GC.Content.File character 用于矫正肿瘤样本LogR值的GC含量文件, GC含量文件中不存在的位点将被剔除; 默认NULL
############' 要求制表符分隔且包含行名列名, 行名为SNP位点的标识符, 各列所包含的信息依次为[SeqName、Position、Windows1、Windows2、Windows3、...], 其中[SeqName]信息不能包含"chr"前缀, 且列名必须为"Chr", 各窗口大小在1Mb以内且依次递增, 必须包含1kb大小和1Mb大小的窗口, 窗口命名规则为整数开头并以"bp/kb/Mb"为单位, 列名整体前端要保留一个制表符
##' @param Replication.Timing.File character 用于矫正肿瘤样本LogR值各细胞系的DNA Replication Timing文件, 需要包含GC含量文件中的所有位点; 默认NULL
############' 要求制表符分隔且包含行名列名, 行名为SNP位点的标识符, 各列所包含的信息依次为[SeqName、Position、Cell.Line1、Cell.Line2、Cell.Line3、...], 其中[SeqName]信息不能包含"chr"前缀, 列名整体前端要保留一个制表符
##' @param Sample.Genders character[] 声明各样本的性别的字符向量, "XX"表示女性, "XY"表示男性; 默认NULL, 即所有样本均视为女性
##' @param Genome.Version character 设置基因组版本号, 可选值为("hg19", "hg38"); 默认"hg19"
##' @param Min.Ploidy numeric 肿瘤样本倍性解搜索空间的最小边界, 倍性小于该值的样本视为分析失败; 默认1.5
##' @param Max.Ploidy numeric 肿瘤样本倍性解搜索空间的最大边界, 倍性大于该值的样本视为分析失败; 默认5.5
##' @param Technology.Parameter.Gamma numeric 平台规范参数, 对于大多数SNP Arrays(例如: Illumina 109k和Affymetrix SNP 6.0)该值为0.55, 对于高通量测序数据, 该值必须为1; 默认0.55
##' @param Array.Platform character SNP Arrays对应的平台, 这将用于在缺少配对的正常样本数据时, 估计正常样本的基因型数据
##' @return list 包含每个样本的肿瘤纯度[Purity]、倍性[Ploidy]、位点的估计信息[SeqName、Position、BAF、LogR]以及片段的估计信息[SeqName、Position.Start、Position.End、CN.Major、CN.Minor]
ASCAT.CNV.Calling <- function(Tumor.BAF.File, Tumor.LogR.File, 
                              Normal.BAF.File = NULL, Normal.LogR.File = NULL, 
                              GC.Content.File = NULL, Replication.Timing.File = NULL, 
                              Sample.Genders = NULL, Genome.Version = c("hg19","hg38"), Min.Ploidy = 1.5, Max.Ploidy = 5.5, Technology.Parameter.Gamma = 0.55, 
                              Array.Platform = c("Affy10k", "AffySNP6", "Affy100k", "Custom10k", "Illumina1M", "IlluminaASA", "HumanCore12", "Illumina109k", "Illumina610k", "Illumina660k", "Illumina700k", "Illumina2.5M", "Affy250k_sty", "Affy250k_nsp", "AffyOncoScan", "IlluminaGSAv3", "IlluminaOmni5", "AffyCytoScanHD", "IlluminaCytoSNP", "HumanCNV370quad", "HumanCoreExome24", "HumanOmniExpress12", "IlluminaCytoSNP850k", "IlluminaOmniExpressExome")){
  library(ASCAT)
  ############
  ## 1.参数判断与预处理
  ############
  Tumor.BAF.File <- as.character(Tumor.BAF.File)
  if(length(Tumor.BAF.File) == 1 && file.exists(Tumor.BAF.File)){
    Tumor.BAF.File <- normalizePath(Tumor.BAF.File, winslash = "/", mustWork = TRUE)
    Tumor.BAF.Sample.Num <- length(unlist(strsplit(trimws(readLines(con = Tumor.BAF.File, n = 1)), "\t"))) - 2
  }else{
    stop("'Tumor.BAF.File'应为单一且存在的文件路径 ...")
  }
  
  Tumor.LogR.File<- as.character(Tumor.LogR.File)
  if(length(Tumor.LogR.File) == 1 && file.exists(Tumor.LogR.File)){
    Tumor.LogR.File <- normalizePath(Tumor.LogR.File, winslash = "/", mustWork = TRUE)
    Tumor.LogR.Sample.Num <- length(unlist(strsplit(trimws(readLines(con = Tumor.LogR.File, n = 1)), "\t"))) - 2
  }else{
    stop("'Tumor.LogR.File'应为单一且存在的文件路径 ...")
  }
  
  if(Tumor.BAF.Sample.Num == Tumor.LogR.Sample.Num){
    Sample.Num <- Tumor.BAF.Sample.Num <- Tumor.LogR.Sample.Num
  }else{
    stop("'Tumor.BAF.File'与'Tumor.LogR.File'应包含的同等的样本数量 ...")
  }
  
  Normal.BAF.File <- as.character(Normal.BAF.File)
  if(length(Normal.BAF.File) == 0){
    Normal.BAF.File <- NULL
  }else if(length(Normal.BAF.File) == 1 && file.exists(Normal.BAF.File)){
    Normal.BAF.File <- normalizePath(Normal.BAF.File, winslash = "/", mustWork = TRUE)
    Normal.BAF.Sample.Num <- length(unlist(strsplit(trimws(readLines(con = Normal.BAF.File, n = 1)), "\t"))) - 2
    if(Normal.BAF.Sample.Num != Sample.Num){
      stop("'Normal.BAF.File'应与'Tumor.BAF.File'以及'Tumor.LogR.File'包含的同等的样本数量 ...")
    }
  }else{
    stop("'Normal.BAF.File'应为NULL或单一且存在的文件路径 ...")
  }
  
  Normal.LogR.File <- as.character(Normal.LogR.File)
  if(length(Normal.LogR.File) == 0){
    Normal.LogR.File <- NULL
  }else if(length(Normal.LogR.File) == 1 && file.exists(Normal.LogR.File)){
    Normal.LogR.File <- normalizePath(Normal.LogR.File, winslash = "/", mustWork = TRUE)
    Normal.LogR.Sample.Num <- length(unlist(strsplit(trimws(readLines(con = Normal.LogR.File, n = 1)), "\t"))) - 2
    if(Normal.LogR.Sample.Num != Sample.Num){
      stop("'Normal.LogR.File'应与'Tumor.BAF.File'以及'Tumor.LogR.File'包含的同等的样本数量 ...")
    }
  }else{
    stop("'Normal.LogR.File'应为NULL或单一且存在的文件路径 ...")
  }
  
  GC.Content.File <- as.character(GC.Content.File)
  if(length(GC.Content.File) == 0){
    GC.Content.File <- NULL
  }else if(length(GC.Content.File) == 1 && file.exists(GC.Content.File)){
    GC.Content.File <- normalizePath(GC.Content.File, winslash = "/", mustWork = TRUE)
  }else{
    stop("'GC.Content.File'应为NULL或单一且存在的文件路径 ...")
  }
  
  Replication.Timing.File <- as.character(Replication.Timing.File)
  if(length(Replication.Timing.File) == 0){
    Replication.Timing.File <- NULL
  }else if(length(Replication.Timing.File) == 1 && file.exists(Replication.Timing.File)){
    Replication.Timing.File <- normalizePath(Replication.Timing.File, winslash = "/", mustWork = TRUE)
  }else{
    stop("'Replication.Timing.File'应为NULL或单一且存在的文件路径 ...")
  }
  
  Sample.Genders <- as.character(Sample.Genders)
  if(length(Sample.Genders) == 0){
    Sample.Genders <- NULL
  }else if(length(Sample.Genders) != Sample.Num || ! all(Sample.Genders %in% c("XX", "XY"))){
    stop(sprintf("'Sample.Genders'应为null或与样本个数等长的character向量, 且全部元素应属于(%s) ...", paste0(c("XX", "XY"), collapse = ",")))
  }
  
  Min.Ploidy <- as.numeric(Min.Ploidy)
  if(length(Min.Ploidy) != 1){
    stop("'Min.Ploidy'应为单一numeric值 ...")
  }
  
  Max.Ploidy <- as.numeric(Max.Ploidy)
  if(length(Max.Ploidy) != 1){
    stop("'Max.Ploidy'应为单一numeric值 ...")
  }
  
  Technology.Parameter.Gamma <- as.numeric(Technology.Parameter.Gamma)
  if(length(Technology.Parameter.Gamma) != 1 || Technology.Parameter.Gamma <= 0 || Technology.Parameter.Gamma > 1){
    stop("'Technology.Parameter.Gamma'应为介于(0,1]之间的单一numeric值 ...")
  }
  
  ############
  ## 2.ASCAT分析
  ############
  # 创建一个临时目录用于存放整个过程产生的临时文件, 当程序结束或意外终止时自动删除该目录
  ASCAT.Temp.Dir <- normalizePath(tempdir(check = TRUE), winslash = "/", mustWork = TRUE)
  on.exit(unlink(ASCAT.Temp.Dir, recursive = TRUE, force = TRUE))
  # 载入数据
  ASCAT.Object <- ascat.loadData(gender = Sample.Genders, 
                                genomeVersion = match.arg(Genome.Version), 
                                Tumor_BAF_file = Tumor.BAF.File, Tumor_LogR_file = Tumor.LogR.File, 
                                Germline_BAF_file = Normal.BAF.File, Germline_LogR_file = Normal.LogR.File)
  
  # 矫正校正肿瘤样本的LogR值
  if(! is.null(GC.Content.File)){
    ASCAT.Object <- ascat.correctLogR(ASCAT.Object, GCcontentfile = GC.Content.File, replictimingfile = Replication.Timing.File)
  }
  # 若缺少配对的正常样本数据, 需要依据SNP阵列平台预测正常样本的基因型数据
  if(is.null(Normal.BAF.File) || is.null(Normal.LogR.File)){
    # 预测正常样本的因型
    ASCAT.Germline.Genotype <- ascat.predictGermlineGenotypes(ASCAT.Object, platform = match.arg(Array.Platform))
  }else{
    ASCAT.Germline.Genotype <- NULL
  }
  # 片段分割
  ASCAT.Object <- ascat.aspcf(ASCAT.Object, ascat.gg = ASCAT.Germline.Genotype, out.dir = ASCAT.Temp.Dir)
  # 估计肿瘤纯度、倍性、等位特异拷贝数等信息
  ASCAT.Output <- ascat.runAscat(ASCAT.Object, gamma = Technology.Parameter.Gamma, pdfPlot = TRUE, min_ploidy = Min.Ploidy, max_ploidy = Max.Ploidy, img.dir = ASCAT.Temp.Dir)
  
  ############
  ## 3.结果封装
  ############
  Analysis.Samples <- ASCAT.Object$samples
  ASCAT.Result <- lapply(Analysis.Samples, function(Analysis.Sample){
    if(Analysis.Sample %in% ASCAT.Output$failedarrays){
      return(
        list(Purity = NULL, 
             Ploidy = NULL, 
             Point.Data = data.frame(SeqName = NULL, Position = NULL, BAF = NULL, LogR = NULL), 
             Segment.Data = data.frame(SeqName = NULL, Position.Start = NULL, Position.End = NULL, CN.Major = NULL, CN.Minor = NULL)
        )
      )
    }else{
      Analysis.Sample.Segment.Data <- ASCAT.Output$segments[ASCAT.Output$segments$sample %in% Analysis.Sample, ]
      return(
        list(
          Purity = as.numeric(ASCAT.Output$purity[Analysis.Sample]), 
          Ploidy = as.numeric(ASCAT.Output$ploidy[Analysis.Sample]), 
          Point.Data = data.frame(SeqName = as.character(ASCAT.Object$SNPpos[, 1]), Position = as.numeric(ASCAT.Object$SNPpos[, 2]), BAF = as.numeric(ASCAT.Object$Tumor_BAF[[Analysis.Sample]]), LogR = as.numeric(ASCAT.Object$Tumor_LogR[[Analysis.Sample]])), 
          Segment.Data = data.frame(SeqName = Analysis.Sample.Segment.Data$chr, Position.Start = Analysis.Sample.Segment.Data$startpos, Position.End = Analysis.Sample.Segment.Data$endpos, CN.Major = Analysis.Sample.Segment.Data$nMajor, CN.Minor = Analysis.Sample.Segment.Data$nMinor)
        )
      )
    }
  })
  names(ASCAT.Result) <- Analysis.Samples
  return(ASCAT.Result)
}


##' @description 从高通量测序数据中计算LogR值与BAF值
##' @param Loci.Prefix character SNP位点坐标数据文件的前缀[可携带路径], 每条染色体对应一个坐标文件, 文件的全路径格式为“[Loci.Prefix][Chromosome][.txt]”
############' 每一个SNP位点坐标数据文件要求是没有列名并以制表符分隔的文件, 包含[Chromosome、Position]两列信息, 其中"Chromosome"信息应与BAM文件或CRAM文件中使用得格式统一, 即是否包含前缀"chr"
##' @param Alleles.Prefix character SNP位点的等位基因数据文件的前缀[可携带路径], 每条染色体对应一个坐标文件, 文件的全路径格式为“[Alleles.Prefix][Chromosome][.txt]”
############'  每一个SNP位点的等位基因数据文件要求是有列名并以制表符分隔的文件, 包含[Position、Ref.Allele、Alt.Allele]三列信息, 其中[Position]信息对应的列名必须为"position", [Ref.Allele、Alt.Allele]分别用数值1、2、3、4表示A、C、G、T四种核苷酸
##' @param HTS.Tumor.File character 肿瘤样本的高通量测序数据, 要求是BAM文件或CRAM文件格式
##' @param HTS.Normal.File character 正常样本的高通量测序数据, 要求是BAM文件或CRAM文件格式
##' @param Tumor.Name character 设置肿瘤样本的名称; 默认NULL, 即HTS.Tumor.File的文件名
##' @param Normal.Name character 设置正常样本的名称; 默认NULL, 即HTS.Normal.File的文件名
##' @param Intervals.Bed character 限制对基因组区域的子集进行操作的基因组区间Bed文件, 仅考虑落在其中的SNP位点; 默认NULL
##' @param Genome.Refence character 参考基因组文件(FASTA格式), 当BCRAM文件的Header中找不到对应的参考基因组文件信息时需要进行指定; 默认NULL
##' @param OutPut.Dir character 结果文件的输出目录; 默认NULL, 即当前工作目录
##' @param System.Allelecounter.Alias character AlleleCounter软件在系统中的可执行命令名; 默认"alleleCounter"
##' @param N.Threads numeric 并行运算所使用的线程数量; 默认1
##' @param Min.Depth numeric 保留SNP位点的最小限制深度, 正常样本中reads数小于该值的SNP位点将被剔除; 默认10
##' @param Min.Map.Quality numeric 设置比对质量的最小阈值(针对MAPQ信息); 默认35
##' @param Min.Base.Quality numeric 设置序列质量的最小阈值(针对QUAL信息); 默认20
##' @param Chromosomes character[] 设置要使用的染色体; 默认c(1:22, "X")
##' @param Sample.Gender character 设置样本性别, 用两条性染色体表示, 可选值为("XX", "XY"); 默认"XX"
##' @param Genome.Version character 设置基因组版本号, 可选值为("hg19", "hg38"); 默认"hg19"
ASCAT.Extract.LogR.BAF <- function(Loci.Prefix, Alleles.Prefix, 
                                   HTS.Tumor.File, HTS.Normal.File, 
                                   Tumor.Name = NULL, Normal.Name = NULL, 
                                   System.Allelecounter.Alias = "alleleCounter", 
                                   Intervals.Bed = NULL, Genome.Refence = NULL, OutPut.Dir = NULL, 
                                   N.Threads = 1, Min.Depth = 10, Min.Map.Quality = 35, Min.Base.Quality = 20, 
                                   Chromosomes = c(1:22, "X"), Sample.Gender = c("XX", "XY"), Genome.Version = c("hg19","hg38")){
  library(ASCAT)
  System.Allelecounter.Alias <- as.character(System.Allelecounter.Alias)
  if(length(System.Allelecounter.Alias) == 1){
    # 判断System.Allelecounter.Alias在系统中是否存在
    if(nchar(Sys.which(System.Allelecounter.Alias)) > 0){
      
      ############
      ## 1.参数判断与预处理
      ############
      Loci.Prefix <- as.character(Loci.Prefix)
      if(length(Loci.Prefix) != 1){
        stop("'Loci.Prefix'应为单一的character值 ...")
      }
      
      Alleles.Prefix <- as.character(Alleles.Prefix)
      if(length(Alleles.Prefix) != 1){
        stop("'Alleles.Prefix'应为单一的character值 ...")
      }
      
      HTS.Tumor.File <- as.character(HTS.Tumor.File)
      if(length(HTS.Tumor.File) == 1 && file.exists(HTS.Tumor.File)){
        HTS.Tumor.File <- normalizePath(HTS.Tumor.File, winslash = "/", mustWork = TRUE)
      }else{
        stop("'HTS.Tumor.File'应为单一且存在的文件路径 ...")
      }
      
      HTS.Normal.File <- as.character(HTS.Normal.File)
      if(length(HTS.Normal.File) == 1 && file.exists(HTS.Normal.File)){
        HTS.Normal.File <- normalizePath(HTS.Normal.File, winslash = "/", mustWork = TRUE)
      }else{
        stop("'HTS.Normal.File'应为单一且存在的文件路径 ...")
      }
      
      Tumor.Name <- as.character(Tumor.Name)
      if(length(Tumor.Name) == 0){
        Tumor.Name <- basename(HTS.Tumor.File)
      }else if(length(Tumor.Name) != 1){
        stop("'Tumor.Name'应为单一的character值 ...")
      }
      
      Normal.Name <- as.character(Normal.Name)
      if(length(Normal.Name) == 0){
        Normal.Name <- basename(HTS.Normal.File)
      }else if(length(Normal.Name) != 1){
        stop("'Normal.Name'应为单一的character值 ...")
      }
      
      Intervals.Bed <- as.character(Intervals.Bed)
      if(length(Intervals.Bed) == 0){
        Intervals.Bed <- NA
      }else if(length(Intervals.Bed) == 1 && file.exists(Intervals.Bed)){
        Intervals.Bed <- normalizePath(Intervals.Bed, winslash = "/", mustWork = TRUE)
      }else{
        stop("'Intervals.Bed'应为NULL或单一且存在的文件路径 ...")
      }
      
      Genome.Refence <- as.character(Genome.Refence)
      if(length(Genome.Refence) == 0){
        Genome.Refence <- NA
      }else if(length(Genome.Refence) == 1 && file.exists(Genome.Refence)){
        Genome.Refence <- normalizePath(Genome.Refence, winslash = "/", mustWork = TRUE)
      }else{
        stop("'Genome.Refence'应为NULL或单一且存在的文件路径 ...")
      }
      
      OutPut.Dir <- as.character(OutPut.Dir)
      if(length(OutPut.Dir) == 0){
        OutPut.Dir <- getwd()
      }else if(length(OutPut.Dir) == 1 && dir.exists(OutPut.Dir)){
        OutPut.Dir <- normalizePath(OutPut.Dir, winslash = "/", mustWork = TRUE)
      }else{
        stop("'OutPut.Dir'应为NULL或单一且存在的目录路径 ...")
      }
      
      N.Threads <- as.numeric(N.Threads)
      if(length(N.Threads) != 1 || N.Threads < 1 || N.Threads  %% 1 != 0){
        stop("'N.Threads'应为单一的大于等于1的整型numeric值 ...")
      }
      
      Min.Depth <- as.numeric(Min.Depth)
      if(length(Min.Depth) != 1 || Min.Depth < 0 || Min.Depth  %% 1 != 0){
        stop("'Min.Depth'应为单一的大于等于0的整型numeric值 ...")
      }
      
      Min.Map.Quality <- as.numeric(Min.Map.Quality)
      if(length(Min.Map.Quality) != 1 || Min.Map.Quality < 0 || Min.Map.Quality  %% 1 != 0){
        stop("'Min.Map.Quality'应为单一的大于等于0的整型numeric值 ...")
      }
      
      Min.Base.Quality <- as.numeric(Min.Base.Quality)
      if(length(Min.Base.Quality) != 1 || Min.Base.Quality < 0 || Min.Base.Quality  %% 1 != 0){
        stop("'Min.Base.Quality'应为单一的大于等于0的整型numeric值 ...")
      }
      
      Chromosomes <- unique(as.character(Chromosomes))
      if(length(Chromosomes) == 0 || ! all(Chromosomes %in% c(1:22, "X", "Y"))){
        stop(sprintf("'Chromosomes'应为至少包含一个元素的character向量, 且全部元素应属于(%s) ...", paste0(c(1:22, "X", "Y"), collapse = ",")))
      }
      
      ############
      ## 2.计算LogR值与BAF值
      ############
      ascat.prepareHTS(allelecounter_exe = System.Allelecounter.Alias, 
                       BED_file = Intervals.Bed, ref.fasta = Genome.Refence, 
                       loci.prefix = Loci.Prefix, alleles.prefix = Alleles.Prefix,
                       tumourseqfile = HTS.Tumor.File, normalseqfile = HTS.Normal.File, tumourname = Tumor.Name,  normalname = Normal.Name, 
                       tumourBAF_file = sprintf("%s/%s.BAF", OutPut.Dir, Tumor.Name), normalBAF_file = sprintf("%s/%s.BAF", OutPut.Dir, Normal.Name), 
                       tumourLogR_file = sprintf("%s/%s.LogR", OutPut.Dir, Tumor.Name), normalLogR_file = sprintf("%s/%s.LogR", OutPut.Dir, Normal.Name),
                       chrom_names = Chromosomes, gender = match.arg(Sample.Gender),genomeVersion = match.arg(Genome.Version), nthreads = N.Threads, minCounts = Min.Depth, min_map_qual = Min.Map.Quality, min_base_qual = Min.Base.Quality)

      ############
      ## 3.删除过程中产生的临时文件(等位基因计数文件)
      ############
      Temp.File.Regular <- sprintf("(%s|%s)_alleleFrequencies_chr(%s)[.]txt", Tumor.Name, Normal.Name, paste0(c(1:22, "X"), collapse = "|"))
      unlink(list.files(pattern = Temp.File.Regular, full.names = T), force = TRUE)
      
    }else{
      stop(sprintf("非系统的可执行命令'%s' ...", System.Allelecounter.Alias))
    }
  }else{
    stop("'System.Allelecounter.Alias'应为单一的character值 ...")
  }
}