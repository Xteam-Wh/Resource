##############################函数描述##############################
# "Genome.View"对基因组上的信号特征(点和线段)进行可视化(至少可视化一条基因组版本包含的序列)
####################################################################


##' @description 对基因组上的信号特征(点和线段)进行可视化(至少可视化一条基因组版本包含的序列)
##' @author Xteam.Wh
##' @param ... list 每个list包含以下元素(其中至少要包含Point.Data与Segment.Data其中的一项, 否则将被从队列中移除)：
############' $Feature.Name character 特征名，将作为对应的纵坐标title属性
############' $Point.Data data.frame 包含必要列[SeqName(序列名), Position(所在序列的位点), Feature.Value(特征信号值)], 可选列[Feature.Type(特征信号所属类别)]
############' $Segment.Data data.frame 包含必要列[SeqName(序列名), Position.Start(所在序列的起始位点), Position.End(所在序列的结束位点), Feature.Value(特征信号值)], 可选列[Feature.Type(特征信号所属类别)]
############' $Point.Size numeric 设置点的尺寸
############' $Point.Shape numeric | character 设置点型
############' $Point.Alpha numeric 设置点的透明度
############' $Point.Color character 设置点的颜色
############' $Point.Fill character 设置点的填充色
############' $Point.Stroke numeric 设置点的边缘线尺寸
############' $Segment.Size numeric 设置线段的尺寸
############' $Segment.Alpha numeric 设置线段的透明度
############' $Segment.Color character 设置线段的颜色
############' $Segment.LineType numeric | character 设置线段的线型
############' $Color.Map characte[] 颜色映射集合, 与"Feature.Type"包含的元素种类相对应
##' @param Feature.List.Data list 特征list数据集合，每个特征信号对应一个list, 每个list包含的元素与可变参数(...)传入的每个list一致
##' @param Auto.Marge logical 是否自动对各图表进行合并
##' @param SeqName.Ratio numeric 指定组合图标中基因组条带图所占的比例, 当且仅当Auto.Marge = TRUE是生效
##' @param Genome.Assemblies character 基因组版本号, 可选unique(c(GenomeInfoDb::registered_UCSC_genomes()$genome, GenomeInfoDb::registered_NCBI_assemblies()$assembly))
##' @return 若Auto.Marge = TRUE, 则返回组合后的绘图信息; 若Auto.Marge = FALSE, 则返每个特征信号的绘图信息以及上下基因组条段绘图信息
Genome.View <- function(..., Feature.List.Data = NULL, Auto.Marge = TRUE, SeqName.Ratio = 0.125,  
                        Genome.Assemblies = unique(c(GenomeInfoDb::registered_UCSC_genomes()$genome, GenomeInfoDb::registered_NCBI_assemblies()$assembly))){
  
  
  ############
  ## 0.特征信号数据的整合
  ############
  library(ggplot2)
  Feature.List.Data <- c(list(...), as.list(Feature.List.Data))
  
  ############
  ## 1.导入相应的基因组版本相关信息
  ############
  # 匹配基因组版本
  Genome.Assemblies <- match.arg(Genome.Assemblies)
  # 查询基因组版本对应的基因组序列信息(名称、长度、...)
  Genome.Seqinfo <- GenomeInfoDb::Seqinfo(genome = Genome.Assemblies)
  Genome.Seqinfo.SeqName <- Genome.Seqinfo@seqnames
  Genome.Seqinfo.SeqLength <- Genome.Seqinfo@seqlengths
  
  ############
  ## 3.提取绘图数据所包含的基因组序列， 并判断这些序列是否存在与所选取的基因组版本中
  ############
  # 提取特征信号数据中包含的基因组序列名称
  Common.SeqName <- unique(unlist(lapply(Feature.List.Data, function(Feature.Data){return(unique(c(as.list(Feature.Data)$Point.Data$SeqName, as.list(Feature.Data)$Segment.Data$SeqName)))})))
  if(all(Common.SeqName %in% Genome.Seqinfo.SeqName)){
    
    ############
    ## 4.格式化需要绘制的基因组序列相关的位置信息
    ############
    Common.SeqName.Info <- data.frame(SeqName = Genome.Seqinfo.SeqName[Genome.Seqinfo.SeqName %in% Common.SeqName], SeqLength = Genome.Seqinfo.SeqLength[Genome.Seqinfo.SeqName %in% Common.SeqName])
    Common.SeqName.Info$Accumulate.SeqLength  <- cumsum(as.numeric(Common.SeqName.Info$SeqLength))
    Common.SeqName.Info$Label.Position <- (2*Common.SeqName.Info$Accumulate.SeqLength - Common.SeqName.Info$SeqLength)/2
    Common.SeqName.Accumulate.Before.Map <- setNames(object = Common.SeqName.Info$Accumulate.SeqLength - Common.SeqName.Info$SeqLength, nm = Common.SeqName.Info$SeqName)
    
    ############
    ## 5.格式化绘图数据基因组位置信息, 以及美学映射信息
    ############
    Feature.List.Data  <- lapply(Feature.List.Data, function(Feature.Data){
      Feature.Data <- as.list(Feature.Data)
      if(is.null(Feature.Data$Point.Data) && is.null(Feature.Data$Segment.Data)){
        warning("当前数据不存在'Point.Data'或'Segment.Data', 不符合要求, 已将其从队列中移 ...", call. = FALSE)
        return(NA)
      }else{
        if(is.null(Feature.Data$Feature.Name)){
          warning("当前数据未设置'Feature.Name'属性, 该属性将作为纵坐标title属性 ...", call. = FALSE)
        }
        if(! is.null(Feature.Data$Point.Data)){
          Feature.Data$Point.Data <- as.data.frame(Feature.Data$Point.Data)
          if(is.null(Feature.Data$Point.Size)){ Feature.Data$Point.Size = 1 }
          if(is.null(Feature.Data$Point.Shape)){ Feature.Data$Point.Shape = 20 }
          if(is.null(Feature.Data$Point.Stroke)){ Feature.Data$Point.Stroke = 1 }
          if(is.null(Feature.Data$Point.Alpha)){ Feature.Data$Point.Alpha = 0.66 }
          if(is.null(Feature.Data$Point.Color)){ Feature.Data$Point.Color = "DarkSlateGray" }
          if(is.null(Feature.Data$Point.Fill)){ Feature.Data$Point.Fill = Feature.Data$Point.Color }
          Feature.Data$Point.Data$Position  <- Feature.Data$Point.Data$Position + Common.SeqName.Accumulate.Before.Map[Feature.Data$Point.Data$SeqName]
        }
        if(! is.null(Feature.Data$Segment.Data)){
          Feature.Data$Segment.Data <- as.data.frame(Feature.Data$Segment.Data)
          if(is.null(Feature.Data$Segment.Size)){ Feature.Data$Segment.Size = 1 }
          if(is.null(Feature.Data$Segment.Alpha)){ Feature.Data$Segment.Alpha = 0.66 }
          if(is.null(Feature.Data$Segment.LineType)){ Feature.Data$Segment.LineType = "solid" }
          if(is.null(Feature.Data$Segment.Color)){ Feature.Data$Segment.Color = ifelse(is.null(Feature.Data$Point.Data), "DarkSlateGray", "OrangeRed") }
          Feature.Data$Segment.Data$Position.End  <- Feature.Data$Segment.Data$Position.End + Common.SeqName.Accumulate.Before.Map[Feature.Data$Segment.Data$SeqName]
          Feature.Data$Segment.Data$Position.Start  <- Feature.Data$Segment.Data$Position.Start + Common.SeqName.Accumulate.Before.Map[Feature.Data$Segment.Data$SeqName]
        }
        return(Feature.Data)
      }
    })
    
    ############
    ## 6.绘制基因组序列条带(上下两条)
    ############
    SeqName.Plots <- lapply(c(top = "top", bottom = "bottom"), function(Postion){
      Common.SeqName.Info.Num <- nrow(Common.SeqName.Info)
      if(Common.SeqName.Info.Num == 1){
        Label.Position.Index <- 1
      }else{
        Label.Position.Index <- switch(Postion, top = seq(1, Common.SeqName.Info.Num, 2), bottom = seq(2, Common.SeqName.Info.Num, 2))
      }
      SeqName.Plot <- ggplot() +
        geom_rect(data = Common.SeqName.Info, aes(xmin = Accumulate.SeqLength - SeqLength + 1, xmax = Accumulate.SeqLength, ymin = -Inf, ymax = Inf, fill = SeqName), alpha = 0.5, show.legend = FALSE) + 
        scale_y_continuous(limits = c(-1, 1), expand = expansion(0)) +
        scale_x_continuous(limits = c(1, max(Common.SeqName.Info$Accumulate.SeqLength)), breaks = Common.SeqName.Info$Label.Position[Label.Position.Index], labels = Common.SeqName.Info$SeqName[Label.Position.Index], expand = expansion(0), position = Postion) + 
        scale_fill_manual(values = setNames(sapply(seq_along(Common.SeqName.Info$SeqName), function(i){ifelse(i %% 2 == 0, 'gray', 'black')}), nm = Common.SeqName.Info$SeqName)) +
        theme_test() +
        theme(
          axis.title = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(face = "bold"), 
          plot.margin = margin(t = switch(Postion, top = 10, bottom = 2.5), r = 10, b = switch(Postion, top = 2.5, bottom = 10), l = 10)
        )
      return(SeqName.Plot)
    })
    
    ############
    ## 7.绘制各特征信号(点、线段)
    ############
    Feature.Plots <- lapply(Feature.List.Data[!is.na(Feature.List.Data)], function(Feature.Data){
      Point.Data <- Feature.Data$Point.Data
      Segment.Data <- Feature.Data$Segment.Data
      Common.SeqName.Info.Num <- nrow(Common.SeqName.Info)
      if(Common.SeqName.Info.Num == 1){
        Top.Position.Index <- 1
        Bottom.Position.Index <- 1
      }else{
        Top.Position.Index <- seq(1, Common.SeqName.Info.Num, 2)
        Bottom.Position.Index <- seq(2, Common.SeqName.Info.Num, 2)
      }
      Feature.Plot <- ggplot() +
        geom_rect(data = Common.SeqName.Info, aes(xmin = Accumulate.SeqLength - SeqLength + 1, xmax = Accumulate.SeqLength, ymin = -Inf, ymax = Inf, fill = SeqName), alpha = 0.05, show.legend = FALSE) +
        scale_x_continuous(limits = c(1, max(Common.SeqName.Info$Accumulate.SeqLength)), breaks = Common.SeqName.Info$Label.Position[Bottom.Position.Index], labels = Common.SeqName.Info$SeqName[Bottom.Position.Index], expand = expansion(0), sec.axis = dup_axis(breaks = Common.SeqName.Info$Label.Position[Top.Position.Index], labels = Common.SeqName.Info$SeqName[Top.Position.Index])) + 
        scale_fill_manual(values = setNames(sapply(seq_along(Common.SeqName.Info$SeqName), function(i){ifelse(i %% 2 == 0, 'gray', 'black')}), nm = Common.SeqName.Info$SeqName))
      if(! is.null(Point.Data)){
        if(is.null(Point.Data$Feature.Type)){
          Feature.Plot <- Feature.Plot + geom_point(data = Point.Data, aes(x = Position, y = Feature.Value), color = Feature.Data$Point.Color, shape = Feature.Data$Point.Shape, fill = Feature.Data$Point.Fill, stroke = Feature.Data$Point.Stroke, size = Feature.Data$Point.Size, alpha = Feature.Data$Point.Alpha)
        }else{
          Feature.Plot <- Feature.Plot + geom_point(data = Point.Data, aes(x = Position, y = Feature.Value, color = Feature.Type), shape = Feature.Data$Point.Shape, fill = Feature.Data$Point.Fill, stroke = Feature.Data$Point.Stroke, size = Feature.Data$Point.Size, alpha = Feature.Data$Point.Alpha)
        }
      }
      if(! is.null(Segment.Data)){
        if(is.null(Segment.Data$Feature.Type)){
          Feature.Plot <- Feature.Plot + geom_segment(data = Segment.Data, aes(x = Position.Start, y = Feature.Value, xend = Position.End, yend = Feature.Value), color = Feature.Data$Segment.Color, size = Feature.Data$Segment.Size, alpha = Feature.Data$Segment.Alpha, linetype = Feature.Data$Segment.LineType)
        }else{
          Feature.Plot <- Feature.Plot + geom_segment(data = Segment.Data, aes(x = Position.Start, y = Feature.Value, xend = Position.End, yend = Feature.Value, color = Feature.Type), size = Feature.Data$Segment.Size, alpha = Feature.Data$Segment.Alpha, linetype = Feature.Data$Segment.LineType)
        }
      }
      if(!is.null(Feature.Data$Color.Map) && !(is.null(Point.Data$Feature.Type) && is.null(Segment.Data$Feature.Type))){
        if("character" %in% class(Feature.Data$Color.Map) && length(Feature.Data$Color.Map) == length(unique(c(Point.Data$Feature.Type, Segment.Data$Feature.Type)))){
          Feature.Plot <- Feature.Plot + scale_color_manual(values = Feature.Data$Color.Map)
        }else{
          stop("当前数据'Color.Map'应为NULL或与‘Feature.Type’包含的元素种类相对应的character向量 ...")
        }
      }
      Feature.Plot <- Feature.Plot + labs(x = NULL, y = Feature.Data$Feature.Name) +  coord_cartesian(clip = "off") + theme_test()
      return(Feature.Plot)
    })
    
    ############
    ## 7.格式话返回的可视化结果
    ############
    Auto.Marge <- as.logical(Auto.Marge)
    if(length(Auto.Marge) == 1){
      if(Auto.Marge){
        Feature.Plots <- lapply(Feature.Plots, function(Feature.Plot){
          Feature.Plot <- Feature.Plot + theme(
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.length.x = unit(0, "pt"),
            legend.title = element_blank(),
            plot.margin = margin(t = 2.5, r = 10, b = 2.5, l = 10)
          )
          return(Feature.Plot)
        })
        SeqName.Ratio <- as.numeric(SeqName.Ratio)
        if(length(SeqName.Ratio) == 1 && SeqName.Ratio > 0 && SeqName.Ratio < 1){
          return(cowplot::plot_grid(plotlist = c(SeqName.Plots[1], Feature.Plots, SeqName.Plots[2]), ncol = 1, rel_heights = c(SeqName.Ratio/2, rep((1 - SeqName.Ratio)/length(Feature.Plots), length(Feature.Plots)), SeqName.Ratio/2), align = "v", axis = "lr"))
        }else{
          stop("'SeqName.Ratio'应为单一的介于(0, 1)之间的numeric值 ...")
        }
        
      }else{
        return(list(SeqName.Plots = SeqName.Plots, Feature.Plots = Feature.Plots))
      }
    }else{
      stop("'Auto.Marge'应为单一的logical值 ...")
    }
    
  }else{
    stop(sprintf("'SeqName'信息应均存在于(%s)", paste0(Genome.Seqinfo.SeqName, collapse = ", ")))
  }
}
