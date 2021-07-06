##############################函数描述##############################
# "Humam.Gene.Conversion.Local"基于org.Hs.eg.db包进行人类基因组的基因转换
# "Humam.Gene.Conversion.Internet"基于biomaRt包进行人类基因组的基因转换(并非支持全部待转类型), 联网操作
####################################################################


##' @description 基于org.Hs.eg.db包进行人类基因组的基因转换
##' @author Xteam.Wh
##' @param Values character[] 待转换的基因集合
##' @param Values.Type character 待转换的基因集合的基因类型
##' @param Conversion.Type character[] 期望转换的基因类型(支持多类别转化)
##' @param Run.Model character 运行模式(自动模式或引导模式), 可选("Auto", "Lead")
##' @return data.frame 对应'Values'集合以及对应各'Conversion.Type'类型的转换后的结果, 未能匹配的基因转换结果为NA
Humam.Gene.Conversion.Local <- function(Values, 
                                        Values.Type, 
                                        Conversion.Type, 
                                        Run.Model = c("Auto", "Lead")){
  library(org.Hs.eg.db)
  Values <- unique(as.character(Values))
  Run.Model <- match.arg(Run.Model)
  Values.Allow.Type <- keytypes(org.Hs.eg.db)
  Conversion.Allow.Type <- columns(org.Hs.eg.db)
  switch(Run.Model,
         Auto = {
           # 判断用户输入的Values.Type是否符合条件
           Values.Type <- as.character(Values.Type)
           if(length(Values.Type) == 1){
             if(! Values.Type %in% Values.Allow.Type){
               stop("'Values.Type'为不适用于本函数的类别名 ...")
             }
           }else{
             stop("'Values.Type'应为单一的character值 ...")
           }
           # 判断用户输入的Conversion.Type是否符合条件
           Conversion.Type <- as.character(Conversion.Type)
           if(length(Conversion.Type) > 0){
             if(! all(Conversion.Type %in% Conversion.Allow.Type)){
               stop("'Conversion.Type'中存在不适用于本函数的类别名 ...")
             }
           }else{
             stop("'Conversion.Type'应为至少存在一个元素的character向量 ...")
           }
         },
         Lead = {
           # 引导用户选择Values.Type
           message(paste0(sprintf("[%s]%s", 1:length(Values.Allow.Type), Values.Allow.Type), collapse = "\n"))
           Values.Type.Index <- readline(prompt = "Tip >> 请从上述列表中选择待转换元素所属类型的编号: ")
           while (! Values.Type.Index %in% 1:length(Values.Allow.Type)) {
             Values.Type.Index  <- readline(prompt = sprintf("Tip >> [%s]不在可选编号中, 请重新选择择待转换元素所属类型的编号: ", Values.Type.Index))
           }
           Values.Type <- Values.Allow.Type[as.numeric(Values.Type.Index)]
           # 引导用户选择Conversion.Type
           message(paste0(sprintf("[%s]%s", 1:length(Conversion.Allow.Type), Conversion.Allow.Type), collapse = "\n"))
           Conversion.Type.Indexs <- unlist(strsplit(readline(prompt = "Tip >> 请从上述信息中选择期望转换类型的编号(多类型各编号之间用空格分隔): "), split = "\\s+", perl = TRUE))
           while (! all(Is.Match <- Conversion.Type.Indexs %in% 1:length(Conversion.Allow.Type))) {
             Conversion.Type.Indexs  <- unlist(strsplit(readline(prompt = sprintf("Tip >> [%s]不在可选编号中, 请重新选择期望转换类型的编号(多类型各编号之间用空格分隔): ", paste0(Conversion.Type.Indexs[! Is.Match], collapse = ",") )), split = "\\s+", perl = TRUE))
           }
           Conversion.Type <- Conversion.Allow.Type[as.numeric(Conversion.Type.Indexs)]
         })
  Conversion.Result <- AnnotationDbi::select(org.Hs.eg.db, keys = Values, columns = Conversion.Type, keytype = Values.Type)
  return(Conversion.Result)
}


##' @description 基于biomaRt包进行人类基因组的基因转换(并非支持全部待转类型), 联网操作
##' @author Xteam.Wh
##' @param Values character[] 待转换的基因集合
##' @param Values.Type character 待转换的基因集合的基因类型
##' @param Conversion.Type character[] 期望转换的基因类型
##' @param Run.Model character 运行模式(自动模式或引导模式), 可选("Auto", "Lead")
##' @return data.frame 对应'Values'集合以及对应各'Conversion.Type'类型的转换后的结果, 未能匹配的基因转换结果为NA
Humam.Gene.Conversion.Internet <- function(Values, 
                                           Values.Type, 
                                           Conversion.Type, 
                                           Run.Model = c("Auto", "Lead")){
  library(biomaRt)
  Values <- unique(as.character(Values))
  Run.Model <- match.arg(Run.Model)
  # 创建Mart对象(选择使用ENSEMBL数据库以及人类基因组数据集)
  Mart.Use <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  Values.Allow.Type <- listFilters(Mart.Use)
  Conversion.Allow.Type <- listAttributes(Mart.Use, what = c("name", "description"))
  Values.Allow.Type <- Values.Allow.Type[Values.Allow.Type$name %in% Conversion.Allow.Type$name, ]
  switch(Run.Model,
         Auto = {
           # 判断用户输入的Values.Type是否符合条件
           Values.Type <- as.character(Values.Type)
           if(length(Values.Type) == 1){
             if(! Values.Type %in% Values.Allow.Type$name){
               stop("'Values.Type'为不适用于本函数的类别名 ...")
             }
           }else{
             stop("'Values.Type'应为单一的character值 ...")
           }
           # 判断用户输入的Conversion.Type是否符合条件
           Conversion.Type <- as.character(Conversion.Type)
           if(length(Conversion.Type) > 0){
             if(! all(Conversion.Type %in%  Conversion.Allow.Type$name)){
               stop("'Conversion.Type'中存在不适用于本函数的类别名 ...")
             }
           }else{
             stop("'Conversion.Type'应为至少存在一个元素的character向量 ...")
           }
         },
         Lead = {
           # 引导用户选择Values.Type
           Values.Allow.Type.File <- normalizePath(tempfile(pattern = "Values_Allow_Type[", tmpdir = getwd(), fileext = "].txt"), winslash = "/", mustWork = FALSE)
           on.exit(unlink(Values.Allow.Type.File, force = TRUE), add = TRUE)
           writeLines(sprintf("[%s]%s ======>> %s", 1:nrow(Values.Allow.Type), Values.Allow.Type$name, Values.Allow.Type$description), con = Values.Allow.Type.File)
           Values.Type.Index <- readline(prompt = sprintf("Tip >> 请从文件'%s'中选择待转换元素所属类型的编号: ", Values.Allow.Type.File))
           while (! Values.Type.Index %in% 1:nrow(Values.Allow.Type)) {
             Values.Type.Index  <- readline(prompt = sprintf("Tip >> [%s]不在可选编号中, 请重新选择择待转换元素所属类型的编号: ", Values.Type.Index))
           }
           Values.Type <- Values.Allow.Type$name[as.numeric(Values.Type.Index)]
           # 引导用户选择Conversion.Type
           Conversion.Allow.Type.File <- normalizePath(tempfile(pattern = "Conversion_Allow_Type[", tmpdir = getwd(), fileext = "].txt"), winslash = "/", mustWork = FALSE)
           on.exit(unlink(Conversion.Allow.Type.File, force = TRUE), add = TRUE)
           writeLines(sprintf("[%s]%s ======>> %s", 1:nrow(Conversion.Allow.Type), Conversion.Allow.Type$name, Conversion.Allow.Type$description), con = Conversion.Allow.Type.File)
           Conversion.Type.Indexs <- unlist(strsplit(readline(prompt = sprintf("Tip >> 请从文件'%s'中选择期望转换类型的编号(多类型各编号之间用空格分隔): ", Conversion.Allow.Type.File)), split = "\\s+", perl = TRUE))
           while (! all(Is.Match <- Conversion.Type.Indexs %in% 1:nrow(Conversion.Allow.Type))) {
             Conversion.Type.Indexs  <- unlist(strsplit(readline(prompt = sprintf("Tip >> [%s]不在可选编号中, 请重新选择期望转换类型的编号(多类型各编号之间用空格分隔): ", paste0(Conversion.Type.Indexs[! Is.Match], collapse = ",") )), split = "\\s+", perl = TRUE))
           }
           Conversion.Type <- Conversion.Allow.Type$name[as.numeric(Conversion.Type.Indexs)]
         })
  # 获取转换结果
  Conversion.Result <- getBM(values = Values, filters = Values.Type, attributes = c(Values.Type, Conversion.Type), mart = Mart.Use)
  return(Conversion.Result)
}
