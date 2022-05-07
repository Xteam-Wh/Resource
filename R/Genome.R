##############################函数描述##############################
# “Calculate.GC.Content”计算位点在参考序列中的GC含量
####################################################################


##' @description 计算位点在参考序列中的GC含量
##' @author Xteam.Wh
##' @param Data matrix | data.frame 包含位点信息的矩阵, 要求至少包含[SeqName、Position], 其中"SeqName"表示位点所在参考序列的片段名称, "Position"表示位点所在参考序列中的位置
##' @param Genome.Refence character 参考基因组文件(FASTA格式)
##' @param Seq.Type character 序列类型, 可选("DNA", "RNA"), 默认"DNA"
##' @param Window.Size numeric 计算GC含量的窗口大小, 以碱基为单位, 默认取位点两侧各"Window.Size/2"的序列进行计算; 默认25
##' @return data.frame 在Data的基础上添加列"GC.Content"
Calculate.GC.Content <- function(Data, Genome.Refence, 
                                 Seq.Type = c("DNA", "RNA"), Window.Size = 25){
  Seq.Type <- match.arg(Seq.Type)
  Window.Size <- as.numeric(Window.Size)
  Genome.Refence <- as.character(Genome.Refence)
  Data <- as.data.frame(Data, row.names = NULL)
  if(all(c("SeqName", "Position") %in% colnames(Data))){
    if(length(Window.Size) == 1 && Window.Size > 1 && Window.Size %% 1 == 0 ){
      if(length(Genome.Refence) == 1 && file.exists(Genome.Refence)){
        library(Biostrings)
        # 将位点信息按序列类型进行分割
        Data$SeqName <- as.character(Data$SeqName)
        Data$Position <- as.numeric(Data$Position)
        Data <- split(Data, Data$SeqName)
        # 载入参考基因组的序列信息
        Genome.Refence.Info <- switch(Seq.Type, DNA = readDNAStringSet(Genome.Refence), RNA = readRNAStringSet(Genome.Refence))
        # 格式化参考基因组的序列名称
        names(Genome.Refence.Info) <- gsub(pattern = " .*$", "", trimws(names(Genome.Refence.Info)))
        # 遍历各染色体, 分别计算位点在大小为Window.Size的窗口内的GC含量
        return(Reduce(rbind, lapply(names(Data), function(SeqName){
          SeqName.Data <- Data[[SeqName]]
          Genome.Refence.SeqName.Info <- Genome.Refence.Info[[SeqName]]
          # 在参考序列两端分别拼接"Window.Size/2"个未定义碱基"N", 以便计算参考序列首末部位碱基的GC含量
          Genome.Refence.SeqName.Info <- xscat(paste0(rep("N", trunc(Window.Size/2)), collapse = ""), Genome.Refence.SeqName.Info, paste0(rep("N", trunc(Window.Size/2)), collapse = ""))
          # 计算参考序大小为Window.Size的窗口内各A、C、G、T(U)四种碱基的数量(以固定窗口每次向后移动一个碱基)
          Base.Type <- switch(Seq.Type, DNA = c("A", "C", "G", "T"), RNA = c("A", "C", "G", "U"))
          Window.Base.Content <- letterFrequencyInSlidingView(Genome.Refence.SeqName.Info, view.width = ifelse(Window.Size%%2 != 0, Window.Size, Window.Size + 1), letters = Base.Type)
          # 计算参考序列每个碱基在大小为Window.Size的窗口内的GC含量
          Window.GC.Content <- rowSums(Window.Base.Content[, c("C", "G")])/rowSums(Window.Base.Content)
          # 提取各位点的在大小为Window.Size的窗口内的GC含量
          SeqName.Data$GC.Content <- Window.GC.Content[SeqName.Data$Position]
          return(SeqName.Data)
        })))
      }else{
        stop("'Genome.Refence'应为单一且存在的文件路径 ...")
      }
    }else{
      stop("'Window.Size'应为单一的大于等于2的整型numeric值 ...")
    }
  }else{
    stop("'Data'应至少包含位点'SeqName'和'Position'两列信息 ...")
  }
}

# https://hgdownload.soe.ucsc.edu/downloads.html
Coordinate.Convert <- function(Data, Chain.File){
  Data <- as.data.frame(Data)
  if(all(c("SeqName", "Position.Start", "Position.End") %in% colnames(Data))){
    Chain.File <- as.character(Chain.File)
    if(length(Chain.File) == 1 && file.exists(Chain.File)){
      Chain <- import.chain(con = Chain.File)
      Chain.Allow.SeqNames <- names(Chain)
      if(all(Data$SeqName %in% Chain.Allow.SeqNames)){
        
        
      }else{
        stop(sprintf("'Data$SeqName'中的元素均应属于(%s) ...", paste0(Chain.Allow.SeqNames, collapse = ",")))
      }
    }else{
      stop("'Chain.File'应为单一且存在的文件路径 ...")
    }
  }else{
    stop("'Data'应至少包含位点'SeqName'、'Position.Start'以及'Position.End'三列信息 ...")
  }
}