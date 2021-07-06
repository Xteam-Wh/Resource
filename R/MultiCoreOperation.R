##############################函数描述##############################
# "Parallel.Execut"基于parallel包进行并行运算
# "Foreach.Execut"基于foreach包与doParallel包进行并行运算
####################################################################


##' @description 基于parallel包进行并行运算
##' @author Xteam.Wh
##' @param Data 用于进行并行运算的数据(数据的适用格式与迭代方式式均取决于参数"Parallel.Type")
##' @param Function funcction 要使用的函数, 函数名或直接进行自定义函数的编写
##' @param Cores.Need numeric 预计要使用的线程(CPU内核)数, 默认2个
##' @param Cores.Type Cores.Type 使用的线程(CPU内核)的类型(逻辑型, 物理型), 可选("Logical", "Physical")
##' @param Common.Packages character[] 并行环境需要共享使用的R包集合, 当且仅当操作系统为Windows时生效, 默认为当前环境已载入的所有R包
##' @param Common.Objects character[] 并行环境需要共享使用的变量集合, 当且仅当操作系统为Windows时生效, 默认为.GlobalEnv环境中存在的所有变量
##' @param Parallel.Type character 进行并行运算的方式, 可选("ParApply", "ParLapply", "ParSapply")效果等同于基础包的("apply", "lapply", "sapply")
##' @param ParApply.MARGIN numeric[] 等同于函数“apply”的参数MARGIN, 当且仅当Parallel.Type = "ParApply"时, 必须设置该参数
##' @param ParSapply.USE.NAMES logical 等同于函数“sapply”的参数USE.NAMES, 当且仅当Parallel.Type = "ParLapply"时有效
##' @return 根据参数Parallel.Type的设定返回与("apply", "lapply", "sapply")相对应的处理结果
Parallel.Execut <- function(Data, Function, Cores.Need = 2, Cores.Type = c("Logical", "Physical"), 
                            Common.Packages = .packages(), Common.Objects = ls(all.names = TRUE, envir = .GlobalEnv), 
                            Parallel.Type = c("ParApply", "ParLapply", "ParSapply"), ParApply.MARGIN, ParSapply.USE.NAMES = TRUE){
  if(! (missing(Data) || is.null(Data))){
    Cores.Need <- as.numeric(Cores.Need)
    if(length(Cores.Need) == 1 && Cores.Need > 1 && Cores.Need %% 1 == 0){
      ParSapply.USE.NAMES <- as.logical(ParSapply.USE.NAMES)
      if(length(ParSapply.USE.NAMES) == 1){
        library(parallel)
        Function <- match.fun(Function)
        Sys.Type <- Sys.info()["sysname"] # 获取当前操作系统的类型
        Cores.Type <- match.arg(Cores.Type)
        Parallel.Type <- match.arg(Parallel.Type)
        Sys.Cores <- detectCores(logical = switch(Cores.Type, Logical = TRUE, Physical = FALSE)) # 获取当前系统可用的线程总数
        if(Cores.Need < Sys.Cores){
          Cores.Use <- Cores.Need
        }else{
          Cores.Use <- Sys.Cores - 1
          message(sprintf("需要的线程数'Cores.Need'不应超过或等于系统的可用线程数[%s], 将自动分配[%s]个线程用于集群的创建 ...", Sys.Cores, Cores.Use))
        }
        # 创建集群(并行运算的运行环境)
        Cluster <- makeCluster( spec = Cores.Use, type = ifelse(Sys.Type == "Windows", "PSOCK", "FORK"))
        on.exit(stopCluster(Cluster), add = TRUE) # 当函数运行结束或发生错误被终止时，自动关闭集群
        # 由于Windows操作系统各线程之间无法共享内存，因此需要将用到的R包和变量在每个线程中均分配一份
        if(Sys.Type == "Windows"){
          Common.Packages <- as.character(Common.Packages)
          Common.Objects <- as.character(Common.Objects)
          assign(".Common.Packages", Common.Packages, envir = globalenv())
          on.exit(remove(.Common.Packages, envir = .GlobalEnv), add = TRUE)
          clusterExport(Cluster, varlist = c(Common.Objects, ".Common.Packages"))
          clusterEvalQ(Cluster, eval(parse(text = paste0(sprintf("library(%s)", .Common.Packages), collapse = ";"))))
        }
        # 依据Type类型，进行并行运算的操作并获取返回结果
        Parallel.Result <- switch(Parallel.Type,
                                  ParApply = parApply(Cluster, Data, MARGIN = ParApply.MARGIN, Function),
                                  ParLapply = parLapply(Cluster, Data, Function),
                                  ParSapply = parSapply(Cluster, Data, Function, USE.NAMES = ParSapply.USE.NAMES))
        return(Parallel.Result)
      }else{
        stop("'ParSapply.USE.NAMES'应为单一的logical值 ...")
      }
    }else{
      stop("'Cores.Need'应为单一的大于1的整型numeric值 ...")
    }
  }else{
    stop("必须指定参数'Data'且不能为空 ...")
  }
  
}


##' @description 基于foreach包与doParallel包进行并行运算
##' @author Xteam.Wh
##' @param Data 用于进行并行运算的数据(对于vector和list对每一个成员进行迭代, 对于matrix和data.frame对每列进行迭代, 对于array则对每个元素进行迭代)
##' @param Function funcction 要使用的函数, 函数名或直接进行自定义函数的编写
##' @param Cores.Need numeric 预计要使用的线程(CPU内核)数, 默认2个
##' @param Cores.Type Cores.Type 使用的线程(CPU内核)的类型(逻辑型, 物理型), 可选("Logical", "Physical")
##' @param Common.Packages character[] 并行环境需要共享使用的R包集合, 当且仅当操作系统为Windows时生效, 默认为当前环境已载入的所有R包
##' @param Common.Objects character[] 并行环境需要共享使用的变量集合, 当且仅当操作系统为Windows时生效, 默认为.GlobalEnv环境中存在的所有变量
##' @return list 每条迭代数据通过函数"Function"操作后返回的结果
Foreach.Execut <- function(Data, Function, Cores.Need = 2, Cores.Type = c("Logical", "Physical"), 
                           Common.Packages = .packages(), Common.Objects = ls(all.names = TRUE, envir = .GlobalEnv)){
  if(! (missing(Data) || is.null(Data))){
    Cores.Need <- as.numeric(Cores.Need)
    if(length(Cores.Need) == 1 && Cores.Need > 1 && Cores.Need %% 1 == 0){
      library(foreach)
      library(doParallel)
      Function <- match.fun(Function)
      Sys.Type <- Sys.info()["sysname"] # 获取当前操作系统的类型
      Cores.Type <- match.arg(Cores.Type)
      Sys.Cores <- detectCores(logical = switch(Cores.Type, Logical = TRUE, Physical = FALSE)) # 获取当前系统可用的线程总数
      if(Cores.Need < Sys.Cores){
        Cores.Use <- Cores.Need
      }else{
        Cores.Use <- Sys.Cores - 1
        message(sprintf("需要的线程数'Cores.Need'不应超过或等于系统的可用线程数[%s], 将自动分配[%s]个线程用于集群的创建 ...", Sys.Cores, Cores.Use))
      }
      # 创建集群(并行运算的运行环境)
      Cluster <- makeCluster( spec = Cores.Use, type = ifelse(Sys.Type == "Windows", "PSOCK", "FORK"))
      on.exit(stopCluster(Cluster), add = TRUE) # 当函数运行结束或发生错误被终止时，自动关闭集群
      registerDoParallel(Cluster) # 进行集群注册
      # 进行并行运算并返回结果
      Common.Packages <- as.character(Common.Packages)
      Common.Objects <- as.character(Common.Objects)
      Foreach.Result <- foreach(Iterative.Data = Data, .packages = Common.Packages, .export = Common.Objects) %dopar% Function(Iterative.Data)
      return(Foreach.Result)
    }else{
      stop("'Cores.Need'应为单一的大于1的整型numeric值 ...")
    }
  }else{
    stop("必须指定参数'Data'且不能为空 ...")
  }
}