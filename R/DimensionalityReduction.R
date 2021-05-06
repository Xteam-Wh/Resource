##############################函数描述##############################
# "DR.PCA"通过FactoMineR包进行主成分分析
####################################################################


##' @description 通过FactoMineR包进行主成分分析
##' @author Xteam.Wh
##' @param Data matrix | data.frame 特征矩阵, 要求行为样本, 列为特征
##' @param Is.Scale logical 是否对特征矩阵进行标准化(默认进行中心化); 默认TRUE
##' @param Min.Contribution numeric 最小累积贡献度限制, [0,100]; 默认80
##' @return list 包含达到Min.Contribution的主成分的样本得分矩阵、特征载荷矩阵以及特征贡献得分矩阵
DR.PCA <- function(Data, Is.Scale = T, Min.Contribution = 80){
  library(FactoMineR)
  Data <- as.data.frame(Data)
  Is.Scale <- as.logical(Is.Scale)
  Min.Contribution <- as.numeric(Min.Contribution)
  if(length(Is.Scale) == 1){
    if(length(Min.Contribution) == 1 && Min.Contribution >= 0 && Min.Contribution <= 100){
      # 进行主成分分析
      PCA.Result <- PCA(Data, ncp = NULL, scale.unit = Is.Scale, graph = FALSE)
      # 提取累积方差贡献度
      Contribution <- PCA.Result[["eig"]][, "percentage of variance"]
      Cumulative.Contribution <- PCA.Result$eig[, "cumulative percentage of variance"]
      # 计算达到最小累积方差贡献度时的主成分个数
      Pc.Number <- min(which(Cumulative.Contribution >= Min.Contribution))
      # 保留Pc.Number个主成分的样本得分矩阵、特征载荷矩阵以及特征贡献得分矩阵
      Pc.Scores <- as.data.frame(PCA.Result[["ind"]][["coord"]][, 1:Pc.Number])
      Pc.Loadings <- as.data.frame(PCA.Result[["svd"]][["V"]][, 1:Pc.Number])
      Pc.Contribution <- as.data.frame(PCA.Result[["var"]][["contrib"]][, 1:Pc.Number])
      # 格式化行名列名
      Pc.Names <- paste(sprintf("Dim.%s ", 1:Pc.Number), ' (', signif(Contribution)[1:Pc.Number], "% variance)", sep='')
      colnames(Pc.Scores) <-  Pc.Names
      colnames(Pc.Loadings) <-  Pc.Names
      colnames(Pc.Contribution) <- Pc.Names
      rownames(Pc.Loadings) <- rownames(Pc.Contribution)
      # 返回相应的主成分的样本得分矩阵、特征载荷矩阵以及特征贡献得分矩阵
      return(
        list(Pc.Scores = Pc.Scores, Pc.Loadings = Pc.Loadings, Pc.Contribution = Pc.Contribution)
      )
    }else{
      stop("'Min.Contribution'应为介于[0,1]之间的单一的numeric值 ...")
    }
  }else{
    stop("'Is.Scale'应为单一的logica值 ...")
  }
}