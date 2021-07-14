##############################函数描述##############################
# "Annovar.Download"用于下载Annovar注释所需数据库的相关信息
# "Annovar.Run"通过Annovar的table_Annovar.pl对给定的avinput或vcf文件进行注释
####################################################################


##' @description 用于下载Annovar注释所需数据库的相关信息
##' @author Xteam.Wh
##' @param Database character 要下载的注释数据库名称
##' @param Annovar.Dir character Annovar的安装目录(必要pl文件所存在的目录); 默认安装目录为当前工作路径
##' @param Database.Dir character 注释数据库信息的下载目录; 默认安装目录为当前"Annovar.Dir"下的"humandb"目录
##' @param System.Perl.Alias character Perl软件在系统中的可执行命令名; 默认为"perl"
##' @param Database.Buildver character 注释数据库对应的基因组版本号, 目前可选则的有("hg18", "hg19", "hg38"); 默认为"hg18"
##' @param Webfrom character 注释数据库的下载源, 可选("annovar", "ucsc"); 默认"annovar"
Annovar.Download <- function(Database, 
                             Annovar.Dir = NULL, Database.Dir = NULL, System.Perl.Alias = "perl", 
                             Database.Buildver = c("hg18", "hg19", "hg38"), Webfrom = c("annovar", "ucsc")){
  System.Perl.Alias <- as.character(System.Perl.Alias)
  if(length(System.Perl.Alias) == 1){
    # 判断System.Perl.Alias在系统中是否存在
    if(nchar(Sys.which(System.Perl.Alias)) > 0){
      
      # 配置Annovar.Dir
      Annovar.Dir <- as.character(Annovar.Dir)
      if(length(Annovar.Dir) == 0){
        # 若Annovar.Dir未配置, 则将其设置为当前工作目录
        Annovar.Dir <- getwd()
      }else{
        if(length(Annovar.Dir) == 1){
          Annovar.Dir <- normalizePath(Annovar.Dir, winslash = "/", mustWork = TRUE)
        }else{
          stop("'Annovar.Dir'应为NULL或单一且存在的目录路径 ...")
        }
      }
      Annovar.Command <- sprintf("\"%s\" \"%s/annotate_variation.pl\"", System.Perl.Alias, Annovar.Dir)
      
      # 配置Database[<table-name>]
      Database <- as.character(Database)
      if(length(Database) == 1){
        Annovar.Command <- sprintf("%s \"%s\"", Annovar.Command, Database)
      }else{
        stop("'Database'应为单一的character值 ...")
      }
      
      # 配置Database.Dir[<database-location>]
      Database.Dir <- as.character(Database.Dir)
      if(length(Database.Dir) == 0){
        Database.Dir <- sprintf("%s/humandb", Annovar.Dir)
      }else{
        if(length(Database.Dir) == 1){
          Database.Dir <- sub(pattern = "(\\\\*/*|/*\\\\*)$", replacement = "", x = Database.Dir)
        }else{
          stop("'Database.Dir'应为NULL或单一且存在的目录路径 ...")
        }
      }
      dir.create(Database.Dir, recursive = TRUE, showWarnings = FALSE)
      Database.Dir <- normalizePath(Database.Dir, winslash = "/", mustWork = TRUE)
      Annovar.Command <- sprintf("%s \"%s\"", Annovar.Command, Database.Dir)
      
      # 运行指令
      System.Command.Run(System.Command = Annovar.Command, Success.Message = sprintf("注释数据库'%s'相关信息已下载至目录'%s' ...", Database, Database.Dir))
      
    }else{
      stop(sprintf("非系统的可执行命令'%s'", System.Perl.Alias))
    }
  }else{
    stop("'System.Perl.Alias'应为单一的character值 ...")
  }
}


##' @description 通过Annovar的table_Annovar.pl对给定的avinput或vcf文件进行注释
##' @author Xteam.Wh
##' @param Annovar.Input character 需要进行注释的文件路径, 必须为avinput或vcf格式文件
############' 其中avinput文件要求是没有列名并以制表符分隔的文件, 且前5列分别为突变位点的[染色体(前缀"chr"可省)、起始位置、终止位置、参考核苷酸、观测核苷酸]信息
##' @param Database character[] 注释所要使用的数据库集合
##' @param Input.Type 需要进行注释的文件类型, 可选("avinput", "vcf"); 默认"avinput"
##' @param Annovar.Dir character Annovar的安装目录(必要pl文件所存在的目录), 设为NULL则默认安装目录为当前工作路径
##' @param Database.Dir character 注释数据库信息的存放目录, 设为NULL则默认为"Annovar.Dir"下的"humandb"目录
##' @param Database.Buildver character 基因组版本号, 可选值为("hg18", "hg19", "hg38"); 默认为"hg18"
##' @param Operation character[] 要进行的注释操作集合("g" for gene-based, "gx" for gene-based with cross-reference annotation, "r" for region-based and "f" for filter-based), 与"Database"一一对应
##' @param Xref.File character 用于基因的交叉引用注释的文件(由制表符分隔, 第一列必须为gene名, 其余列为对应基因的其他注释信息, 文件如有列名, 首行前需添加"#'作为列名行的识别符), 若"Operation"不包含"gx", 则该参数将失去作用
##' @param Do.Polish logical 是否要通过"encoding_change.pl"重新计算蛋白质序列来修饰indel蛋白质注释, 默认TRUE
##' @param Nastring character 设置对于缺失信息使用的填充符号, 默认为"."(对于vcf文件必须使用".", 因此当"Input.Type"为"vcf"时, 该参数将失去作用)
##' @param Output.Prefix character 注释结果文件的前缀名[可携带路径]; 默认为当前工作目录下的"Annovar.Result"
##' @param With.Otherinfo logical 对于传入的avinput文件, 注释结果是否保留前五列以外的信息, 默认FALS(当"Input.Type"为"vcf"时, 该参数将失去作用）
##' @param Temp.Remove logical 是否删除注释过程中产生的临时文件, 默认TRUE
##' @param Csvout logical 注释结果是否以csv格式输出; 默认FALSE, 即默认结果是以制表符分隔的txt文件输出(对于vcf文件将输出转化后的avinput文件以及txt与vcf两种格式的结果文件, 因此当"Input.Type"为"vcf"时, 该参数将失去作用)
##' @param System.Perl.Alias character Perl软件在系统中的可执行命令名; 默认为"perl"
##' @param Other.Options.Settings 其他参数的设置, 将被拼接到指令中进行调用
Annovar.Run <- function(Annovar.Input, Database, 
                        Input.Type = c("avinput", "vcf"),
                        Annovar.Dir = NULL, Database.Dir = NULL, 
                        Database.Buildver = c("hg18", "hg19", "hg38"), 
                        Operation = NULL, Xref.File = NULL, Do.Polish = TRUE, Nastring = ".", Output.Prefix = NULL, 
                        With.Otherinfo = FALSE, Temp.Remove = TRUE, Csvout = FALSE, System.Perl.Alias = "perl", Other.Options.Settings = NULL){
  System.Perl.Alias <- as.character(System.Perl.Alias)
  if(length(System.Perl.Alias) == 1){
    # 判断System.Perl.Alias在系统中是否存在
    if(nchar(Sys.which(System.Perl.Alias)) > 0){
      
      # 配置Annovar.Dir
      Annovar.Dir <- as.character(Annovar.Dir)
      if(length(Annovar.Dir) == 0){
        # 若Annovar.Dir未配置, 则将其设置为当前工作目录
        Annovar.Dir <- getwd()
      }else if(length(Annovar.Dir) == 1){
        Annovar.Dir <- normalizePath(Annovar.Dir, winslash = "/", mustWork = TRUE)
      }else{
        stop("'Annovar.Dir'应为NULL或NULL或单一且存在的目录路径 ...")
      }
      Annovar.Command <- sprintf("\"%s\" \"%s/table_annovar.pl\"", System.Perl.Alias, Annovar.Dir)
      
      # 配置Annovar.Input[<query-file>]
      Annovar.Input <- as.character(Annovar.Input)
      if(length(Annovar.Input) == 1){
        Annovar.Command <- sprintf("%s \"%s\"", Annovar.Command, normalizePath(Annovar.Input, winslash = "/", mustWork = TRUE))
      }else{
        stop("'Annovar.Input'应为单一且存在的文件路径 ...")
      }
      
      # 配置Database.Dir[<database-location>]
      Database.Dir <- as.character(Database.Dir)
      if(length(Database.Dir) == 0){
        Annovar.Command <- sprintf("%s \"%s\"", Annovar.Command, normalizePath(sprintf("%s/humandb", Annovar.Dir), winslash = "/", mustWork = TRUE))
      }else if(length(Database.Dir) == 1){
        Annovar.Command <- sprintf("%s \"%s\"", Annovar.Command, normalizePath(Database.Dir, winslash = "/", mustWork = TRUE))
      }else{
        stop("'Database.Dir'应为NULL或单一且存在的目录路径 ...")
      }
      
      # 配置Database.Buildver[--buildver]
      Annovar.Command <- sprintf("%s --buildver \"%s\"", Annovar.Command, match.arg(Database.Buildver))
      
      # 配置Database[--protocol]
      Database <- as.character(Database)
      Database.Num <- length(Database)
      if(Database.Num > 0){
        Annovar.Command <- sprintf("%s --protocol \"%s\"", Annovar.Command, paste0(Database, collapse = ","))
      }else{
        stop("'Database'应为至少包含一个元素的character向量 ...")
      }
      
      # 配置Operation[--operation]
      Operation <- as.character(Operation)
      if(length(Operation) > 0 && all(Operation %in% c("g", "gx", "r", "f"))){
        if(length(Operation) == 1){
          Annovar.Command <- sprintf("%s --operation \"%s\"", Annovar.Command, paste0(rep(Operation, Database.Num), collapse = ","))
        }else if(length(Operation) == Database.Num){
          Annovar.Command <- sprintf("%s --operation \"%s\"", Annovar.Command, paste0(Operation, collapse = ","))
        }else{
          stop("'Operation'应为单一的character值或与'Database'等长的character向量 ...")
        }
      }else{
        stop("'Operation'应为至少包含一个元素的character向量, 且所有元素均应选自('g' for gene-based, 'gx' for gene-based with cross-reference annotation, 'r' for region-based and 'f' for filter-based) ...")
      }
      
      # 配置Xref.File[--xreffile]
      Xref.File <- as.character(Xref.File)
      if(grepl("gx", Operation)){
        if(length(Xref.File) > 0){
          if(length(Xref.File) == 1){
            Annovar.Command <- sprintf("%s --xreffile \"%s\"", Annovar.Command, normalizePath(Xref.File, winslash = "/", mustWork = TRUE))
          }else{
            "'Xref.File'应为NULL或单一的文件路径 ..."
          }
        }else{
          stop("'Operation'中含有'gx', 应通过'Xref.File'为其设置要使用的交叉引用注释文件 ...")
        }
      }
      
      # 配置Do.Polish[--polish]
      Do.Polish <- as.logical(Do.Polish)
      if(length(Do.Polish) == 1){
        if(Do.Polish){
          Annovar.Command <- sprintf("%s --polish", Annovar.Command)
        }else{
          Annovar.Command <- sprintf("%s --nopolish", Annovar.Command)
        }
      }else{
        stop("'Do.Polish'应为单一的logical值 ...")
      }
      
      # 配置Nastring[--nastring]
      Nastring <- as.character(Nastring)
      if(match.arg(Input.Type) == "avinput"){
        if(length(Nastring) == 1){
          Annovar.Command <- sprintf("%s --nastring \"%s\"", Annovar.Command, Nastring)
        }else{
          stop("'Nastring'应为单一的character值 ...")
        }
      }
      
      # 配置With.Otherinfo[--otherinfo]
      With.Otherinfo <- as.logical(With.Otherinfo)
      if(match.arg(Input.Type) == "avinput"){
        if(length(With.Otherinfo) == 1){
          if(With.Otherinfo){
            Annovar.Command <- sprintf("%s --otherinfo", Annovar.Command)
          }
        }else{
          stop("'With.Otherinfo'应为单一的logical值 ...")
        } 
      }
      
      # 配置Temp.Remove[--remove]
      Temp.Remove <- as.logical(Temp.Remove)
      if(length(Temp.Remove) == 1){
        if(Temp.Remove){
          Annovar.Command <- sprintf("%s --remove", Annovar.Command)
        }
      }else{
        stop("'ATemp.Remove'应为单一的logical值 ...")
      }
      
      # 配置[--vcfinput]、Csvout[--csvout]
      Csvout <- as.logical(Csvout)
      if(length(Csvout) == 1){
        if(match.arg(Input.Type) == "vcf"){
          Annovar.Command <- sprintf("%s --vcfinput", Annovar.Command)
        }else{
          if(Csvout){
            Annovar.Command <- sprintf("%s --csvout", Annovar.Command)
          }
        }
      }else{
        stop("'Csvout'应为单一的logical值 ...")
      }
      
      # 配置Other.Options.Settings[其他参数]
      Other.Options.Settings <- as.character(Other.Options.Settings)
      if(length(Other.Options.Settings) > 0){
        if(length(Other.Options.Settings) == 1){
          Annovar.Command <- sprintf("%s %s", Annovar.Command, trimws(Other.Options.Settings))
        }else{
          stop("'Other.Options.Settings'应为Null或单一的character值 ...")
        }
      }
      
      # 配置Output.Prefix[--outfile]
      Output.Prefix <- as.character(Output.Prefix)
      if(length(Output.Prefix) == 0){
        Output.Prefix <- sprintf("%s/Annovar.Result", getwd())
      }else if(length(Output.Prefix) == 1){
        Output.Prefix <- normalizePath(sub(pattern = "(\\\\*/*|/*\\\\*)$", replacement = "", x = trimws(Output.Prefix)), winslash = "/", mustWork = FALSE)
        dir.create(dirname(Output.Prefix), recursive = TRUE, showWarnings = FALSE)
      }else{
        stop("'Output.Prefix'应为Null或单一的character值 ...")
      }
      Annovar.Command <- sprintf("%s --outfile \"%s\"", Annovar.Command, Output.Prefix)
      
      # 运行指令
      System.Command.Run(System.Command = Annovar.Command, Success.Message = sprintf("结果文件已输出至'%s.Xxx' ...", Output.Prefix))
      
    }else{
      stop(sprintf("非系统的可执行命令'%s' ...", System.Perl.Alias))
    }
  }else{
    stop("'System.Perl.Alias'应为单一的character值 ...")
  }
}