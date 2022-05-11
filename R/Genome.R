##############################函数描述##############################
# "Genome.Coord.Convert"基因组坐标转换
# “Calculate.GC.Content”计算特定位点在参考序列中的GC含量
####################################################################


##' @description 特定位点在不同基因组之间的坐标转换
##' @author Xteam.Wh
##' @param Data matrix | data.frame 包含位点信息的矩阵, 要求至少包含[SeqName(序列名), Position.Start(所在序列的起始位点), Position.End(所在序列的结束位点)]三列信息
##' @param Chain.File character 用于基因组坐标转换的Chain文件, 可从UCSC进行下载(https://hgdownload.soe.ucsc.edu/downloads.html)
############' 进入下载界面后点击当前基因组对应的"LiftOver files"连接, 之后将与目标基因组进行坐标转换的文件下载到本地并解压
##' @return list 包含Data每行数据的转换结果, 每一个元素为一个data.frame, 包含转换后的[SeqName(序列名), Position.Start(所在序列的起始位点), Position.End(所在序列的结束位点)]三列信息
Genome.Coord.Convert <- function(Data, Chain.File){
  library(liftOver)
  Data <- as.data.frame(Data)
  if(all(c("SeqName", "Position.Start", "Position.End") %in% colnames(Data))){
    Chain.File <- as.character(Chain.File)
    if(length(Chain.File) == 1 && file.exists(Chain.File)){
      # 导入用于基因组坐标转换的Chain文件
      Chain <- import.chain(con = Chain.File)
      Chain.Allow.SeqNames <- names(Chain)
      if(all(Data$SeqName %in% Chain.Allow.SeqNames)){
        # 格式化代转换的基因组位点坐标信息
        Genome.Ranges <- GRanges(seqnames = Data$SeqName, ranges=IRanges(start = Data$Position.Start, end = Data$Position.End, names = rownames(Data)))
        # 进行基因组位点坐标转换
        Genome.Converted.Ranges <- liftOver(Genome.Ranges, chain = Chain)
        # 格式化输出结果
        return(
          lapply(Genome.Converted.Ranges, function(Genome.Converted.Range){
            Temp <- as.data.frame(Genome.Converted.Range)
            if(nrow(Temp) > 0){
              return(data.frame(SeqName = Temp$seqnames, Position.Start = Temp$start, Position.End = Temp$end))
            }else{
              return(data.frame(SeqName = NA, Position.Start = NA, Position.End = NA))
            }
          })
        )
      }else{
        stop(sprintf("'Data$SeqName'中的元素均应属于(%s) ...", paste0(Chain.Allow.SeqNames, collapse = ",")))
      }
    }else{
      stop("'Chain.File'应为单一且存在的文件路径 ...")
    }
  }else{
    stop("'Data'应至少包含[SeqName(序列名), Position.Start(所在序列的起始位点), Position.End(所在序列的结束位点)]三列信息 ...")
  }
}


##' @description 计算特定位点在参考序列中的GC含量
##' @author Xteam.Wh
##' @param Data matrix | data.frame 包含位点信息的矩阵, 要求至少包含[SeqName(序列名), Position(所在序列的位点)]两列信息
##' @param Genome.Refence character 参考基因组文件(FASTA格式)
##' @param Seq.Type character 序列类型, 可选("DNA", "RNA"), 默认"DNA"
##' @param Window.Size numeric 计算GC含量的窗口大小, 以碱基为单位, 默认取位点两侧各"Window.Size/2"的序列进行计算; 默认25
##' @return data.frame 在Data的基础上添加列"GC.Content"
Genome.GC.Content <- function(Data, Genome.Refence, 
                              Seq.Type = c("DNA", "RNA"), Window.Size = 25){
  library(Biostrings)
  Seq.Type <- match.arg(Seq.Type)
  Window.Size <- as.numeric(Window.Size)
  Genome.Refence <- as.character(Genome.Refence)
  Data <- as.data.frame(Data, row.names = NULL)
  if(all(c("SeqName", "Position") %in% colnames(Data))){
    if(length(Window.Size) == 1 && Window.Size > 1 && Window.Size %% 1 == 0 ){
      if(length(Genome.Refence) == 1 && file.exists(Genome.Refence)){
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
    stop("'Data'应至少包含[SeqName(序列名), Position(所在序列的位点)]两列信息 ...")
  }
}



