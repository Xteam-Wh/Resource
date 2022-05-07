##############################函数描述##############################
# "ASCAT.CNV.Calling"结合LogR值与BAF值进行拷贝数推断
# “ASCAT.Extract.LogR.BAF”从高通量测序数据中计算SNP位点处的LogR值与BAF值
####################################################################

##' @description 从高通量测序数据中计算LogR值与BAF值
##' @param Loci.Prefix character SNP位点坐标数据文件的前缀[可携带路径], 每条染色体对应一个坐标文件, 文件的全路径格式为“[Loci.Prefix][chr][Chromosome][.txt]”
############' 每一个SNP位点坐标数据文件要求是没有列名并以制表符分隔的文件, 包含[Chromosome、Position]两列信息, 其中"Chromosome"信息应与BAM文件或CRAM文件中使用得格式统一, 即是否包含前缀"chr"
##' @param Alleles.Prefix character SNP位点的等位基因数据文件的前缀[可携带路径], 每条染色体对应一个坐标文件, 文件的全路径格式为“[Alleles.Prefix][Chromosome][.txt]”
############'  每一个SNP位点的等位基因数据文件要求是有列名并以制表符分隔的文件, 包含[Position、Ref.Allele、Alt.Allele]三列信息, 其中[Position]信息对应的列名必须为"position", [Ref.Allele、Alt.Allele]分别用数值1、2、3、4表示A、C、G、T四种核苷酸
##' @param HTS.Tumor.File character 肿瘤样本的高通量测序数据, 要求是BAM文件或CRAM文件格式
##' @param HTS.Normal.File character 正常样本的高通量测序数据, 要求是BAM文件或CRAM文件格式
##' @param Genome.Refence character 参考基因组文件(FASTA格式), 当输入的通量测序数据类型为CRAM格式时, 必须设置该参数; 默认NULL
##' @param OutPut.Dir character 结果文件的输出目录; 默认NULL, 即当前工作目录
##' @param System.Allelecounter.Alias character AlleleCounter软件在系统中的可执行命令名; 默认"alleleCounter"
##' @param N.Threads numeric 并行运算所使用的线程数量; 默认1
##' @param Min.Depth numeric 保留SNP位点的最小限制深度, 正常样本中reads数小于该值的SNP位点将被剔除; 默认10
##' @param Min.Map.Quality numeric 设置比对质量的最小阈值(针对MAPQ信息); 默认35
##' @param Min.Base.Quality numeric 设置序列质量的最小阈值(针对QUAL信息); 默认20
##' @param Chromosomes character[] 设置要使用的染色体; 默认c(1:22, "X")
##' @param Gender character 设置样本性别, 用两条性染色体表示, 可选值为("XX", "XY"); 默认"XX"
##' @param Genome.Version character 设置基因组版本号, 可选值为("hg19", "hg38"); 默认"hg19"
ASCAT.Extract.LogR.BAF <- function(Loci.Prefix, Alleles.Prefix, 
                                   HTS.Tumor.File, HTS.Normal.File, 
                                   Genome.Refence = NULL, OutPut.Dir = NULL, 
                                   System.Allelecounter.Alias = "alleleCounter",
                                   N.Threads = 1, Min.Depth = 10, Min.Map.Quality = 35, Min.Base.Quality = 20, 
                                   Chromosomes = c(1:22, "X"), Gender = c("XX", "XY"), Genome.Version = c("hg19","hg38")){
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
      if(length(HTS.Tumor.File) == 1){
        HTS.Tumor.File <- normalizePath(HTS.Tumor.File, winslash = "/", mustWork = TRUE)
      }else{
        stop("'HTS.Tumor.File'应为单一且存在的文件路径 ...")
      }
      
      HTS.Normal.File <- as.character(HTS.Normal.File)
      if(length(HTS.Normal.File) == 1){
        HTS.Normal.File <- normalizePath(HTS.Normal.File, winslash = "/", mustWork = TRUE)
      }else{
        stop("'HTS.Tumor.File'应为单一且存在的文件路径 ...")
      }
      
      Genome.Refence <- as.character(Genome.Refence)
      if(length(HTS.Tumor.File) == 0){
        Genome.Refence <- NA
      }else if(length(HTS.Tumor.File) == 1){
        HTS.Tumor.File <- normalizePath(HTS.Tumor.File, winslash = "/", mustWork = TRUE)
      }else{
        stop("'HTS.Tumor.File'应为NULL或单一且存在的文件路径 ...")
      }
      
      OutPut.Dir <- as.character(OutPut.Dir)
      if(length(OutPut.Dir) == 0){
        OutPut.Dir <- getwd()
      }else if(length(OutPut.Dir) == 1){
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
      ascat.prepareHTS(loci.prefix = Loci.Prefix, 
                       alleles.prefix = Alleles.Prefix,
                       tumourseqfile = HTS.Tumor.File,
                       normalseqfile = HTS.Normal.File,
                       tumourname = basename(HTS.Tumor.File), 
                       normalname = basename(HTS.Normal.File),
                       allelecounter_exe = System.Allelecounter.Alias,
                       gender = match.arg(Gender), genomeVersion = match.arg(Genome.Version), 
                       tumourBAF_file = sprintf("%s/%s.BAF", OutPut.Dir, basename(HTS.Tumor.File)),
                       normalBAF_file = sprintf("%s/%s.BAF", OutPut.Dir, basename(HTS.Normal.File)),
                       tumourLogR_file = sprintf("%s/%s.LogR", OutPut.Dir, basename(HTS.Tumor.File)),
                       normalLogR_file = sprintf("%s/%s.LogR", OutPut.Dir, basename(HTS.Normal.File)),
                       chrom_names = Chromosomes, nthreads = N.Threads, minCounts = Min.Depth, min_map_qual = Min.Map.Quality, min_base_qual = Min.Base.Quality)
      
      
      ############
      ## 3.删除过程中产生的临时文件(等位基因计数文件)
      ############
      Temp.File.Regular <- sprintf("(%s|%s)_alleleFrequencies_chr(%s)[.]txt", basename(HTS.Tumor.File), basename(HTS.Normal.File), paste0(c(1:22, "X"), collapse = "|"))
      unlink(list.files(pattern = Temp.File.Regular, full.names = T), force = TRUE)
      
    }else{
      stop(sprintf("非系统的可执行命令'%s' ...", System.Allelecounter.Alias))
    }
  }else{
    stop("'System.Allelecounter.Alias'应为单一的character值 ...")
  }
}


ASCAT.Extract.LogR.BAF <- function(Tumor.BAF.File, Tumor.LogR.File, 
                                   Normal.BAF.File = NULL, Normal.LogR.File = NULL, 
                                   GC.Model.File = NULL, Replication.Timing.File = NULL, 
                                   Gamma = 1, Gender = c("XX", "XY"), Genome.Version = c("hg19","hg38")){
  
}