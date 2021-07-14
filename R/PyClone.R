##############################函数描述##############################
# "PyClone.Run"通过R函数传参调用PyClone进行分析
####################################################################


##' @description 通过R函数传参调用PyClone进行分析
##' @author Xteam.Wh
##' @param PyClone.Input characcter[] tsv格式文件集合, 每个tsv文件都应该至少包含以下六列信息:
###############' mutation_id ——> 突变位点的唯一标识, 跨数据集应该是相同的(个人建议格式为: "gene:Chromosomes:position")
###############' ref_counts  ——> 突变位点与reference allele相匹配的locus的reads数
###############' var_counts  ——> 突变位点与variant allele相匹配的locus的reads数
###############' normal_cn   ——> 正常细胞该位点对应locus的拷贝数, 对于人类常染色体来说是2(性染色体除外)
###############' minor_cn    ——> 肿瘤细胞该位点对应minor allele的拷贝数
###############' major_cn    ——> 肿瘤细胞该位点对应major allele的拷贝数
##' @param PyClone.Output.Dir characcter PyClone分析结果的存放目录; 默认为当前工作路径下的"PyClone.Output.Dir"目录
##' @param Config.Extras.File characcter 具有用于分析的额外参数的yaml配置文件的路径
##' @param Samples characcter[] 设置样本的名称, 应与"PyClone.Input"保持一一对应的关系
##' @param Seed integer 随机种子, 便于结果的重现
##' @param Max.Clusters integer 最大聚类簇数量, 设置该参数后, PyClone聚类结果的聚类簇数目不会超过该值
##' @param Min.Cluster.Size integer 最小聚类簇容量, 设置该参数后, PyClone聚类结果中样本含量低于该值的簇将不参与可视化
##' @param Tumour.Contents nurmic 肿瘤纯度, 应与"PyClone.Input"保持一一对应的关系; 若只给定一个值, 则认为所有样本的肿瘤纯度一致
##' @param Num.Iters integer MCMC总迭代次数; 默认10000
##' @param Burnin integer MCMC不稳定迭代的数量(最初的一些MCMC迭代是不稳定的, 在后续分析时需要去除), 一般设置为MCMC总迭代次数的10%; 默认0
##' @param Thin integer MCMC迭代稀薄数, 在去除前"Burnin"个不稳定迭代后, 从剩余的MCMC迭代中, 每"Thin"个迭代中取第一个用于后续分析; 默认1
##' @param Mesh.Size integer 用于近似聚类后验的点数; 默认101
##' @param System.PyClone.Alias character PyClone软件在系统中的可执行命令名, 默认为"PyClone"
##' @param Plot.File.Format character 图片的输出格式, 支持"pdf"与"svg"; 默认"pdf"
##' @param Init.Method character 初始化Dirichlet Process(DP)聚类算法, 可选("disconnected", "connected"), "disconnected"将每个突变位点单独视为一个集群, "connected"将所有位点视为一个集群; 默认"disconnected"
##' @param Density character 计算后验概率采用的密度分布函数, 可选("pyclone_beta_binomial", "pyclone_binomial"), "pyclone_beta_binomial"为beta二项分布, "pyclone_binomial"为二项分布;  默认"pyclone_beta_binomial"
##' @param Prior character 设置拷贝数的先验方式, 可选("major_copy_number", "parental_copy_number", "total_copy_number"), 当minor_cn与major_cn已知时使用"major_copy_number" 或 "parental_copy_number", 当minor_cn与major_cn未知且采用作者的建议将minor_cn设为0、major_cn设置为total_cn时, 使用"total_copy_number"; 
###############' 在模糊程度(突变拷贝的不确定性)上: "parental" < "major" < "total" 即"parental"假设突变位点数量为[1 | mimor_cn(>0) | major_cn], "major"假设突变位点数量为[1 ~ major_cn]之间, "total"假设突变位点数量为[1 ~ total_cn]之间.
PyClone.Run <- function(PyClone.Input, 
                        PyClone.Output.Dir = NULL, Config.Extras.File = NULL, 
                        Samples = NULL, Seed = NULL, Max.Clusters = NULL,  Min.Cluster.Size = NULL, 
                        Tumour.Contents = 1, Num.Iters = 10000, Burnin = 0, Thin = 1, Mesh.Size = 101, 
                        System.PyClone.Alias = "PyClone", 
                        Plot.File.Format = c("pdf", "svg"), 
                        Init.Method = c("disconnected", "connected"), 
                        Density = c("pyclone_beta_binomial", "pyclone_binomial"), 
                        Prior = c("major_copy_number", "parental_copy_number", "total_copy_number")){
  System.PyClone.Alias <- as.character(System.PyClone.Alias)
  if(length(System.PyClone.Alias) == 1){
    # 判断System.PyClone.Alias在系统中是否存在
    if(nchar(Sys.which(System.PyClone.Alias)) > 0){
      
      # 初始化PyClone指令
      PyClone.Command <- sprintf("\"%s\" run_analysis_pipeline", System.PyClone.Alias)
      
      # 配置PyClone.Input[--in_files]
      PyClone.Input <- as.character(PyClone.Input)
      if(length(PyClone.Input) > 0){
        PyClone.Command <- sprintf("%s --in_files %s", PyClone.Command, paste0(sprintf("\"%s\"", normalizePath(PyClone.Input, winslash = "/", mustWork = TRUE)), collapse = " "))
      }else{
        stop("'PyClone.Input'应为至少包含一个元素的文件集合, 且各文件应已经存在 ...")
      }
      
      
      # 配置Samples[--samples]
      Samples <- as.character(Samples)
      if(length(Samples) > 0){
        if(length(Samples) == length(PyClone.Input)){
          PyClone.Command <- sprintf("%s --samples %s", PyClone.Command, paste0(sprintf("\"%s\"", Samples), collapse = " "))
        }else{
          stop("'Samples'应为与'PyClone.Input'等长的character向量 ...")
        }
      }
      
      # 配置Tumour.Contents[--tumour_contents]
      Tumour.Contents <- as.numeric(Tumour.Contents)
      if(length(Tumour.Contents) == 1 && Tumour.Contents> 0 && Tumour.Contents <= 1){
        PyClone.Command <- sprintf("%s --tumour_contents %s", PyClone.Command, paste0(rep(format(Tumour.Contents, scientific = FALSE), length(PyClone.Input)), collapse = " ")) 
      }else if(length(Tumour.Contents) == length(PyClone.Input) && all(Tumour.Contents> 0 & Tumour.Contents <= 1)){
        PyClone.Command <- sprintf("%s --tumour_contents %s", PyClone.Command, paste0(format(Tumour.Contents, scientific = FALSE), collapse = " ")) 
      }else{
        stop("'Tumour_contents'应为介于(0,1]之间的单一数值或与'PyClone.Input'等长的numeric向量 ...")
      }
      
      # 配置Seed[--seed]
      Seed <- as.numeric(Seed)
      if(length(Seed) > 0){
        if(length(Seed) == 1 && Seed %% 1 == 0){
          PyClone.Command <- sprintf("%s --seed %s", PyClone.Command, format(Seed, scientific = FALSE))
        }else{
          stop("'Seed'应为NULL或单一的整型numeric值 ...")
        }
      }
      
      # 配置Max.Clusters[--max_clusters]
      Max.Clusters <- as.numeric(Max.Clusters)
      if(length(Max.Clusters) > 0){
        if(length(Max.Clusters) == 1 && Max.Clusters > 0 && Max.Clusters %% 1 == 0){
          PyClone.Command <- sprintf("%s --max_clusters %s", PyClone.Command, format(Max.Clusters, scientific = FALSE))
        }else{
          stop("'Max.Clusters'应为NULL或单一的大于0的整型numeric值 ...")
        }
      }
      
      # 配置Min.Cluster.Size[--min_cluster_size]
      Min.Cluster.Size <- as.numeric(Min.Cluster.Size)
      if(length(Min.Cluster.Size) > 0){
        if(length(Min.Cluster.Size) == 1 && Min.Cluster.Size > 0 && Min.Cluster.Size %% 1 == 0){
          PyClone.Command <- sprintf("%s --min_cluster_size %s", PyClone.Command, format(Min.Cluster.Size, scientific = FALSE))
        }else{
          stop("'Min.Cluster.Size'应为NULL或单一的大于0的整型numeric值 ...")
        }
      }
      
      # 配置Num.Iters[--num_iters]
      Num.Iters <- as.numeric(Num.Iters)
      if(length(Num.Iters) == 1 && Num.Iters > 0 && Num.Iters %% 1 == 0){
        PyClone.Command <- sprintf("%s --num_iters %s", PyClone.Command, format(Num.Iters, scientific = FALSE))
      }else{
        stop("'Num.Iters'应为单一的大于0的整型numeric值 ...")
      }
      
      # 配置Burnin[--burnin]
      Burnin <- as.numeric(Burnin)
      if(length(Burnin) == 1 && Burnin >= 0 && Burnin < Num.Iters && Burnin %% 1 == 0){
        PyClone.Command <- sprintf("%s --burnin %s", PyClone.Command, format(Burnin, scientific = FALSE))
      }else{
        stop("'Burnin'应为单一的大于0且小于'Num.Iters'的整型numeric值 ...")
      }
      
      # 配置Thin[--thin]
      Thin <- as.numeric(Thin)
      if(length(Thin) == 1 && Thin <= Num.Iters - Burnin && Thin %% 1 == 0){
        PyClone.Command <- sprintf("%s --thin %s", PyClone.Command, format(Thin, scientific = FALSE))
      }else{
        stop("'Thin'应为小于等于'Num.Iters - Burnin'的单一整型numeric值  ...")
      }
      
      # 配置参数Mesh.Size[--mesh_size]
      Mesh.Size <- as.numeric(Mesh.Size)
      if(length(Mesh.Size) == 1 && Mesh.Size > 0 && Mesh.Size %% 1 == 0){
        PyClone.Command <- sprintf("%s --mesh_size %s", PyClone.Command, format(Mesh.Size, scientific = FALSE))
      }else{
        stop("'Mesh.Size'应为单一的大于0的整型numeric值 ...")
      }
      
      # 配置Prior[--prior]、Density[--density]、Init.Method[--init_method]、Plot.File.Format[--plot_file_format]
      PyClone.Command <- sprintf("%s --prior \"%s\" --density \"%s\" --init_method \"%s\" --plot_file_format \"%s\"", PyClone.Command, match.arg(Prior), match.arg(Density), match.arg(Init.Method), match.arg(Plot.File.Format))
      
      # 配置Config.Extras.File[--config_extras_file]
      Config.Extras.File <- as.character(Config.Extras.File)
      if(length(Config.Extras.File) > 0){
        if(length(Config.Extras.File) == 1){
          PyClone.Command <- sprintf("%s --config_extras_file %s", PyClone.Command, normalizePath(Config.Extras.File, winslash = "/", mustWork = TRUE))
        }else{
          stop("'Config.Extras.File'应为NULL或单一且存在的文件路径 ...")
        }
      }
      
      # 配置PyClone.Output.Dir[--working_dir]
      PyClone.Output.Dir <- as.character(PyClone.Output.Dir)
      if(length(PyClone.Output.Dir) == 0){
        PyClone.Output.Dir <- sprintf("%s/PyClone.Output.Dir", getwd())
      }else if(length(PyClone.Output.Dir) == 1){
        PyClone.Output.Dir <- normalizePath(sub(pattern = "(\\\\*/*|/*\\\\*)$", replacement = "", x = trimws(PyClone.Output.Dir)), winslash = "/", mustWork = FALSE)
      }else{
        stop("'PyClone.Output.Dir'应为单一的目录路径 ...")
      }
      dir.create(PyClone.Output.Dir, recursive = TRUE, showWarnings = FALSE)
      PyClone.Output.Dir <- normalizePath(PyClone.Output.Dir, winslash = "/", mustWork = TRUE)
      PyClone.Command <- sprintf("%s --working_dir \"%s\"", PyClone.Command, PyClone.Output.Dir)
      
      # 运行指令
      System.Command.Run(System.Command = PyClone.Command, Success.Message = sprintf("结果文件已输出至目录'%s'  ...", PyClone.Output.Dir))
      
    }else{
      stop(sprintf("非系统的可执行命令'%s' ...", System.PyClone.Alias))
    }
  }else{
    stop("'System.PyClone.Alias'应为单一的character值 ...")
  }
}