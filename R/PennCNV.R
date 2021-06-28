##############################函数描述##############################
# "PennCNV.MakePFB"通过PennCNV的compile_pfb.pl从包含BAF值的多个信号强度文件编译PFB文件
# "PennCNV.Calling"通过PennCNV的detect_cnv.pl结合HMM文件与PFB文件对给定的信号强度文件进行CNV识别
####################################################################

##' @description 通过PennCNV的compile_pfb.pl从包含BAF值的多个信号强度文件编译PFB文件
##' @author Xteam.Wh
##' @param PennCNV.Input character[] 制表符分割的信号强度文件集合, 每个文件必须包含三列信息[Name(标记名称, 必须为首列), *.B Allele Freq(BAF)]; 若未提供Position.File, 则还需包含信息[Chr(染色体), Position(位置)]; 列名需严格按照要求设置
##' @param Position.File character 指定标记坐标信息文件, 文件必须包含三列信息[Name(标记名称, 必须为首列), Chr(染色体), Position(位置)]; 列名需严格按照要求设置, 不存在于该文件中的标记信息不会计算PFB值
############' [注: 若设置了Position.File, 则PennCNV默认仅处理信号强度文件与Position.File共有的标记, 即共有的Name信息所在行]
##' @param Output.Prefix character 结果文件的前缀名[可携带路径]; 默认为当前工作目录下的"PennCNV.MakePFB.Result"
##' @param PennCNV.Dir character PennCNV的安装目录(必要pl文件所存在的目录), 设为NULL则默认安装目录为当前工作路径
##' @param System.Perl.Alias character Perl软件在系统中的可执行命令名; 默认为"perl"
PennCNV.MakePFB <- function(PennCNV.Input, 
                            Position.File = NULL, Output.Prefix = NULL, 
                            PennCNV.Dir = NULL, System.Perl.Alias = "perl"){
  System.Perl.Alias <- as.character(System.Perl.Alias)
  if(length(System.Perl.Alias) == 1){
    # 判断System.Perl.Alias在系统中是否存在
    if(nchar(Sys.which(System.Perl.Alias)) > 0){
      
      # 配置PennCNV.Dir
      PennCNV.Dir <- as.character(PennCNV.Dir)
      if(length(PennCNV.Dir) == 0){
        # 若PennCNV.Dir未配置, 则将其设置为当前工作目录
        PennCNV.Dir <- getwd()
      }else if(length(PennCNV.Dir) == 1){
        PennCNV.Dir <- normalizePath(PennCNV.Dir, winslash = "/", mustWork = TRUE)
      }else{
        stop("'PennCNV.Dir'应为单一且存在的目录路径 ...")
      }
      PennCNV.Command <- sprintf("\"%s\" \"%s/compile_pfb.pl\"", System.Perl.Alias, PennCNV.Dir)
      
      # 配置PennCNV.Input[--listfile]
      PennCNV.Input <- as.character(PennCNV.Input)
      if(length(PennCNV.Input) > 0){
        Input.List.File <- sprintf("%s/Input.List.File", getwd())
        write(normalizePath(PennCNV.Input, winslash = "/", mustWork = TRUE), file = Input.List.File, sep = "\n")
        PennCNV.Command <- sprintf("%s --listfile \"%s\"", PennCNV.Command, normalizePath(Input.List.File, winslash = "/", mustWork = TRUE))
      }else{
        stop("'PennCNV.Input'应为至少包含一个元素的文件集合, 且各文件应已经存在 ...")
      }
      
      # 配置Position.File[--snpposfile]
      Position.File <- as.character(Position.File)
      if(length(Position.File) > 0){
        if(length(Position.File) == 1){
          PennCNV.Command <- sprintf("%s --snpposfile \"%s\"", PennCNV.Command, normalizePath(Position.File, winslash = "/", mustWork = TRUE))
        }else{
          stop("'Position.File'应为NULL或单一且存在的文件路径 ...")
        }
      }
      
      # 配置Output.Prefix[--output]
      Output.Prefix <- as.character(Output.Prefix)
      if(length(Output.Prefix) == 0){
        Output.Prefix <- sprintf("%s/PennCNV.MakePFB.Result", getwd())
      }else if(length(Output.Prefix) == 1){
        Output.Prefix <- normalizePath(sub(pattern = "(\\\\*/*|/*\\\\*)$", replacement = "", x = trimws(Output.Prefix)), winslash = "/", mustWork = FALSE)
        dir.create(dirname(Output.Prefix), recursive = TRUE, showWarnings = FALSE)
      }else{
        stop("'Output.Prefix'应为Null或单一且存在的文件路径 ...")
      }
      PennCNV.Output <- sprintf("%s.pfb", Output.Prefix)
      file.create(PennCNV.Output, showWarnings = FALSE)
      PennCNV.Output <- normalizePath(PennCNV.Output, winslash = "/", mustWork = TRUE)
      PennCNV.Command <- sprintf("%s --output \"%s\"", PennCNV.Command, PennCNV.Output)
      
      # 根据操作系统环境设置脚本内容
      PennCNV.Command <- sprintf(ifelse(Sys.info()["sysname"] == "Windows", "@echo off\n%s", "#!/bin/sh\n%s"),  PennCNV.Command)
      # 根据操作系统环境设置脚本文件
      PennCNV.Command.File <- sprintf(ifelse(Sys.info()["sysname"] == "Windows", "%s/PennCNV.Command.bat", "%s/PennCNV.Command.sh"), getwd())
      # 将指令写入对应系统的脚本文件
      write(PennCNV.Command, PennCNV.Command.File)
      # 赋予脚本文件读写以及可执行权限
      Sys.chmod(PennCNV.Command.File)
      
      # 运行PennCNV
      tryCatch(
        {
          message("<<====== RUNNING MESSAGE ======>>")
          # 执行脚本文件PennCNV.Command.File
          PennCNV.Command.Run <- system(PennCNV.Command.File)
        },
        error = function(e){ # 抛出错误信息
          message(sprintf("<<====== ERROR MESSAGE ======>>\n%s", e))
        },
        warning = function(w){ # 抛出警告信息
          message(sprintf("Warning: %s ...", trimws(gsub(".*\\)\\:", "", w))))
        },
        finally = {
          # 如果脚本文件顺利执行，则PennCNV.Command.Run返回的状态信息为0
          if(exists("PennCNV.Command.Run") && PennCNV.Command.Run == 0){
            message(sprintf("<<===== SUCCESS MESSAGE =====>>\nPennCNV.Command执行成功, 结果已输出至文件'%s' ...", PennCNV.Output))
            unlink(c(Input.List.File, PennCNV.Command.File), force = TRUE)
          }else{
            message(sprintf("<<====== ERROR MESSAGE ======>>\nPennCNV.Command执行过程中发生了错误, 请通过查看'RUNNING MESSAGE'中的信息或通过控制台运行脚本文件'%s'来查看具体错误 ...", PennCNV.Command.File))
          }
        }
      )
      
    }else{
      stop(sprintf("非系统的可执行命令'%s' ...", System.Perl.Alias))
    }
  }else{
    stop("'System.Perl.Alias'应为单一的character值 ...")
  }
}


##' @description 通过PennCNV的detect_cnv.pl结合hmm文件与pfb文件对给定的信号强度文件进行CNV识别
##' @author Xteam.Wh
##' @param PennCNV.Input character[] 制表符分割的信号强度文件集合, 每个文件必须包含三列信息[Name(标记名称, 必须为首列), *.Log R Ratio(LRR), *.B Allele Freq(BAF)]; 若Coordinate.From.Input = FALSE, 则还需包含信息[Chr(染色体), Position(位置)]; 列名需严格按照要求设置
##' @param HMM.File character 指定HMM(Hidden Markov Model)文件, 具体形式参考(http://penncnv.openbioinformatics.org/en/latest/user-guide/input/#hmm-file)
##' @param PFB.File character 指定PFB(Population frequency of B allele)文件, 提供每个标记的PFB信息, 必须包含四列信息[Name(标记名称), Chr(染色体), Position(位置), PFB(B等位基因的群体频率)], 列名可有可无, 但四列信息顺序需严格按照要求设置
############' [注: PennCNV默认仅处理信号强度文件与PFB文件共有的标记, 即共有的Name信息所在行]
##' @param Target.Chr character 处理的染色体类型, 可选("Autosome", "Chr.X", "Chr.Y"); 默认"Autosome"
##' @param PennCNV.Dir character PennCNV的安装目录(必要pl文件所存在的目录), 设为NULL则默认安装目录为当前工作路径
##' @param Output.Prefix character 结果文件的前缀名[可携带路径]; 默认为当前工作目录下的"PennCNV.Calling.Result"
##' @param Output.Tab logical 是否将结果格式化为以制表符分隔的文件; 默认FALSE
##' @param GC.Model.File character 指定GCModel文件, 提供每个标记记周围1Mb基因组区的GC含量(每边500kb), 用于矫正基因组波(genomic waves)的影响, 必须包含四列信息[Name(标记名称), Chr(染色体), Position(位置), GC(GC含量百分比, 0-100)], 列名可有可无, 但四列信息顺序需严格按照要求设置
##' @param Sample.Sex character[] 指定信号强度文件对应样本的性别(male/female), 用于X染色体CNV的推断, 与'PennCNV.Input'保持一一对应的关系, 则认为所有样本性别一致; 也可由程序自动识别性别, 但推荐通过该参数指定
##' @param With.Confidence logical 是否为每个识别到的CNV计算一个置信度分数; 默认FALSE
##' @param PennCNV.Analysis.Type character 指定分析的信号强度文件的类型, 可选("CHIP","WGS"), "CHIP"表示数据来自芯片平台, "WGS"表示数据由PennCNV-Seq从全基因组测序数据转换而来的; 默认"CHIP"
##' @param Coordinate.From.Input character 是否从输入的信号强度文件中来获取标记位点的位置信息; 默认FALSE, 即从PFB文件中获取标记位点的位置信息
##' @param System.Perl.Alias character Perl软件在系统中的可执行命令名; 默认为"perl"
##' @param Other.Options.Settings 其他参数的设置, 将被拼接到指令中进行调用
PennCNV.Calling <- function(PennCNV.Input, 
                        HMM.File, PFB.File, 
                        Target.Chr = c("Autosome", "Chr.X", "Chr.Y"),
                        PennCNV.Dir = NULL, Output.Prefix = NULL, Output.Tab = FALSE,
                        GC.Model.File = NULL, Sample.Sex = NULL, With.Confidence = FALSE, 
                        Coordinate.From.Input = FALSE, System.Perl.Alias = "perl", Other.Options.Settings = NULL){
  System.Perl.Alias <- as.character(System.Perl.Alias)
  if(length(System.Perl.Alias) == 1){
    # 判断System.Perl.Alias在系统中是否存在
    if(nchar(Sys.which(System.Perl.Alias)) > 0){
      
      # 配置PennCNV.Dir
      PennCNV.Dir <- as.character(PennCNV.Dir)
      if(length(PennCNV.Dir) == 0){
        # 若PennCNV.Dir未配置, 则将其设置为当前工作目录
        PennCNV.Dir <- getwd()
      }else if(length(PennCNV.Dir) == 1){
        PennCNV.Dir <- normalizePath(PennCNV.Dir, winslash = "/", mustWork = TRUE)
      }else{
        stop("'PennCNV.Dir'应为单一且存在的目录路径 ...")
      }
      PennCNV.Command <- sprintf("\"%s\" \"%s/detect_cnv.pl\" --test", System.Perl.Alias, PennCNV.Dir)
      
      # 配置Target.Chr[--chrx / --chry]
      PennCNV.Command <- trimws(sprintf("%s %s", PennCNV.Command, switch(match.arg(Target.Chr), Autosome = "", Chr.X = "--chrx", Chr.Y = "--chry")))
      
      # 配置HMM.File[--hmmfile]
      HMM.File <- as.character(HMM.File)
      if(length(HMM.File) == 1){
        PennCNV.Command <- sprintf("%s --hmmfile \"%s\"", PennCNV.Command, normalizePath(HMM.File, winslash = "/", mustWork = TRUE)) 
      }else{
        stop("'HMM.File'应为单一且存在的文件路径 ...")
      }
      
      # 配置PFB.File[--pfbfile]
      PFB.File <- as.character(PFB.File)
      if(length(PFB.File) == 1){
        PennCNV.Command <- sprintf("%s --pfbfile \"%s\"", PennCNV.Command, normalizePath(PFB.File, winslash = "/", mustWork = TRUE)) 
      }else{
        stop("'PFB.File'应为单一且存在的文件路径 ...")
      }
      
      # 配置GC.Model.File[--gcmodelfile]
      GC.Model.File <- as.character(GC.Model.File)
      if(length(GC.Model.File) > 0){
        if(length(GC.Model.File) == 1){
          PennCNV.Command <- sprintf("%s --gcmodelfile \"%s\"", PennCNV.Command, normalizePath(GC.Model.File, winslash = "/", mustWork = TRUE)) 
        }else{
          stop("'GC.Model.File'应为NULL或单一且存在的文件路径 ...")
        }
      }
      
      # 配置PennCNV.Input[--listfile]
      PennCNV.Input <- as.character(PennCNV.Input)
      if(length(PennCNV.Input) > 0){
        Input.List.File <- sprintf("%s/Input.List.File", getwd())
        write(normalizePath(PennCNV.Input, winslash = "/", mustWork = TRUE), file = Input.List.File, sep = "\n")
        PennCNV.Command <- sprintf("%s --listfile \"%s\"", PennCNV.Command, normalizePath(Input.List.File, winslash = "/", mustWork = TRUE))
      }else{
        stop("'PennCNV.Input'应为至少包含一个元素的文件集合, 且各文件应已经存在 ...")
      }
      
      # 配置Sample.Sex[--sexfile]
      Sample.Sex <- as.character(Sample.Sex)
      if(length(Sample.Sex) > 0){
        if(all(Sample.Sex %in% c("male", "female"))){
          if(length(Sample.Sex) == 1 || length(Sample.Sex) == length(PennCNV.Input)){
            Sample.Sex.File <- sprintf("%s/Sample.Sex.File", getwd())
            write(sprintf("%s\t%s", normalizePath(PennCNV.Input, winslash = "/", mustWork = TRUE), Sample.Sex), file = Sample.Sex.File, sep = "\n")
            PennCNV.Command <- sprintf("%s --sexfile \"%s\"", PennCNV.Command, normalizePath(Sample.Sex.File, winslash = "/", mustWork = TRUE))
          }else{
            stop("'Sample.Sex'应为NULL或单一的character值或与'PennCNV.Input'等长的character向量 ...")
          }
        }else{
          stop("'Sample.Sex'的所有元素均应选自('male', 'female') ...")
        }
      }
      
      # 配置Coordinate.From.Input[--coordinate_from_input]
      Coordinate.From.Input <- as.logical(Coordinate.From.Input)
      if(length(Coordinate.From.Input) == 1){
        if(Coordinate.From.Input){
          PennCNV.Command <- sprintf("%s --coordinate_from_input", PennCNV.Command)
        }
      }else{
        stop("'Coordinate.From.Input'应为单一的logical值 ...")
      }
      
      # 配置With.Confidence[--confidence]
      With.Confidence <- as.logical(With.Confidence)
      if(length(With.Confidence) == 1){
        if(With.Confidence){
          PennCNV.Command <- sprintf("%s --confidence", PennCNV.Command)
        }
      }else{
        stop("'With.Confidence'应为单一的logical值 ...")
      }
      
      # 配置Other.Options.Settings[其他参数]
      Other.Options.Settings <- as.character(Other.Options.Settings)
      if(length(Other.Options.Settings) > 0){
        if(length(Other.Options.Settings) == 1){
          PennCNV.Command <- sprintf("%s %s", PennCNV.Command, trimws(Other.Options.Settings))
        }else{
          stop("'Other.Options.Settings'应为Null或单一的character值 ...")
        }
      }
      
      # 配置Output.Tab[--tabout]
      Output.Tab <- as.logical(Output.Tab)
      if(length(Output.Tab) == 1){
        if(Output.Tab){
          PennCNV.Command <- sprintf("%s --tabout", PennCNV.Command)
        }
      }else{
        stop("'Output.Tab'应为单一的logical值 ...")
      }
      
      # 配置Output.Prefix[--output]
      Output.Prefix <- as.character(Output.Prefix)
      if(length(Output.Prefix) == 0){
        Output.Prefix <- sprintf("%s/PennCNV.Calling.Result", getwd())
      }else if(length(Output.Prefix) == 1){
        Output.Prefix <- normalizePath(sub(pattern = "(\\\\*/*|/*\\\\*)$", replacement = "", x = trimws(Output.Prefix)), winslash = "/", mustWork = FALSE)
        dir.create(dirname(Output.Prefix), recursive = TRUE, showWarnings = FALSE)
      }else{
        stop("'Output.Prefix'应为Null或单一且存在的文件路径 ...")
      }
      PennCNV.Output <- sprintf("%s.txt", Output.Prefix)
      file.create(PennCNV.Output, showWarnings = FALSE)
      PennCNV.Output <- normalizePath(PennCNV.Output, winslash = "/", mustWork = TRUE)
      PennCNV.Command <- sprintf("%s --output \"%s\"", PennCNV.Command, PennCNV.Output)
      
      # 根据操作系统环境设置脚本内容
      PennCNV.Command <- sprintf(ifelse(Sys.info()["sysname"] == "Windows", "@echo off\n%s", "#!/bin/sh\n%s"),  PennCNV.Command)
      # 根据操作系统环境设置脚本文件
      PennCNV.Command.File <- sprintf(ifelse(Sys.info()["sysname"] == "Windows", "%s/PennCNV.Command.bat", "%s/PennCNV.Command.sh"), getwd())
      # 将指令写入对应系统的脚本文件
      write(PennCNV.Command, PennCNV.Command.File)
      # 赋予脚本文件读写以及可执行权限
      Sys.chmod(PennCNV.Command.File)
      
      # 运行PennCNV
      tryCatch(
        {
          message("<<====== RUNNING MESSAGE ======>>")
          # 执行脚本文件PennCNV.Command.File
          PennCNV.Command.Run <- system(PennCNV.Command.File)
        },
        error = function(e){ # 抛出错误信息
          message(sprintf("<<====== ERROR MESSAGE ======>>\n%s", e))
        },
        warning = function(w){ # 抛出警告信息
          message(sprintf("Warning: %s ...", trimws(gsub(".*\\)\\:", "", w))))
        },
        finally = {
          # 如果脚本文件顺利执行，则PennCNV.Command.Run返回的状态信息为0
          if(exists("PennCNV.Command.Run") && PennCNV.Command.Run == 0){
            message(sprintf("<<===== SUCCESS MESSAGE =====>>\nPennCNV.Command执行成功, 结果已输出至文件'%s' ...", PennCNV.Output))
            unlink(na.omit(c(Input.List.File, PennCNV.Command.File, ifelse(length(Sample.Sex) > 0, Sample.Sex.File, ""))), force = TRUE)
          }else{
            message(sprintf("<<====== ERROR MESSAGE ======>>\nPennCNV.Command执行过程中发生了错误, 请通过查看'RUNNING MESSAGE'中的信息或通过控制台运行脚本文件'%s'来查看具体错误 ...", PennCNV.Command.File))
          }
        }
      )
      
    }else{
      stop(sprintf("非系统的可执行命令'%s' ...", System.Perl.Alias))
    }
  }else{
    stop("'System.Perl.Alias'应为单一的character值 ...")
  }
}