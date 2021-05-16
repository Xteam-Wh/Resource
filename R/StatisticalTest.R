##############################函数描述##############################
# "Statistical.Test.TwoDiff"计算两组数据之间的差异显著性
# "Statistical.Test.TableDiff"对给定四格表计算差异显著性
# "Statistical.Test.MultiDiff"给定分组信息, 计算各组分数据之间差异显著性
# "Statistical.Test.Survdiff"给定分组信息, 计算各组分之间样本生存差异显著性
####################################################################


##' @description 计算两组数据之间的差异显著性
##' @author  Xteam.Wh
##' @param Data1 numeric[] 样本数据值向量
##' @param Data2 numeric[] 样本数据值向量
##' @param Alternative character 指定差异检验的检验方式(双边|单侧), 可选("two.sided", "less", "greater"); 默认"two.sided"
##' @param Paired logical 是否进行配对差异检验; 默认FALSE
##' @param Normality logical 各组分数据是否符合正态分布, 若符合正态则进行参数检验(T检验)，否则使用非参数检验(Wilcoxon检验); 默认TRUE
##' @param Var.Equal logical 各组分数据是否满足方差齐性, 仅适用于参数检验, 若满足则进行Student T 检验, 否则进行Welch T 检验; 默认TRUE
##' @return numeric 差异显著性P值
Statistical.Test.TwoDiff <- function(Data1, Data2, 
                                     Alternative = c("two.sided", "less", "greater"), 
                                     Paired = FALSE, Normality = TRUE, Var.Equal = TRUE){
  Data1 <- as.numeric(Data1)
  Data2 <- as.numeric(Data2)
  Paired <- as.logical(Paired)
  Normality <- as.logical(Normality)
  Var.Equal <- as.logical(Var.Equal)
  Alternative <- match.arg(Alternative)
  if(length(Paired) == 1 && length(Normality) == 1 && length(Var.Equal) == 1 ){
    if(Normality){
      return(t.test(x = Data1, y = Data2, alternative = Alternative, paired = Paired, var.equal = Var.Equal)$p.value)
    }else{
      return(wilcox.test(x = Data1, y = Data2, alternative = Alternative, paired = Paired)$p.value)
    }
  }else{
    stop("'Paired'、'Normality'与'Var.Equal'均应为单一的logical值 ...")
  }
}


##' @description 对给定四格表计算差异显著性
##' @author  Xteam.Wh
##' @param Data matrix | data.frame 待检验的四格表
##' @param Method character 指定对四格表进行差异检验的方法, 可选("Chisq", "Fisher"); 默认"Chisq"
##' @param Chisq.Correct logical 是否使用连续性矫正卡方检验, 当且仅当Method = "Chisq"时有效; 默认TRUE
##' @param Fisher.Alternative character 指定Fisher精确检验的检验方式(双边|单侧), 当且仅当Method = "Fisher"时有效, 可选("two.sided", "less", "greater"); 默认"two.sided"
##' @return numeric 差异显著性P值
Statistical.Test.TableDiff <- function(Data, 
                                       Method = c("Chisq", "Fisher"), Chisq.Correct = TRUE, 
                                       Fisher.Alternative = c("two.sided", "less", "greater")){
  Data <- as.matrix(Data)
  Method <- match.arg(Method)
  Chisq.Correct <- as.logical(Chisq.Correct)
  Fisher.Alternative <- match.arg(Fisher.Alternative)
  if(all(dim(Data) == 2)){
    if(length(Chisq.Correct) == 1){
      N <- sum(Data)
      T.Min <- min(rowSums(Data))*min(colSums(Data))/N
      if(N >= 40 && T.Min > 5){
        message(sprintf("由于样本量N = %s >= 40, 且T.Min = %s > 5; 建议使用Pearson卡方检验", N, T.Min))
      }else if(N >= 40 && T.Min >= 1 && T.Min < 5){
        message(sprintf("由于样本量N = %s >= 40, 且1 <= T.Min = %s < 5; 建议使用连续性矫正卡方检验", N, T.Min))
      }else{
        message(sprintf("由于样本量N = %s < 40, 且T.Min = %s < 1; 建议用Fisher精确检验", N, T.Min))
      }
      return(
        switch (Method,
                Chisq = chisq.test(Data, correct = Chisq.Correct)$p.value,
                Fisher = fisher.test(Data, alternative = Fisher.Alternative)$p.value
        )
      )
    }else{
      stop("'Chisq.Correct'应为单一的logical值 ...")
    }
  }else{
    stop("'Data'应为一个2*2的matrix或data.frame ...")
  }
}


##' @description 给定分组信息, 计算各组分数据之间差异显著性
##' @author  Xteam.Wh
##' @param Data numeric[] 样本数据值向量
##' @param Group character[] 样本分组信息, 即每个样本对应的组分
##' @param Normality logical 各组分数据是否符合正态分布, 若符合正态则进行参数检验(Anova|Welch)，否则使用非参数检验(Kruskal-Wallis); 默认TRUE
##' @param Var.Equal logical 各组分数据是否满足方差齐性, 仅适用于参数检验, 若满足则进行单因素方差分析F检验, 否则进行Welch检验; 默认TRUE
##' @return numeric 差异显著性P值
Statistical.Test.MultiDiff <- function(Data, Group, Normality = TRUE, Var.Equal = TRUE){
  Data <- as.numeric(Data)
  Group <- as.character(Group)
  Normality <- as.logical(Normality)
  Var.Equal <- as.logical(Var.Equal)
  if(length(Normality) == 1 && length(Var.Equal) == 1 ){
    if(length(Data) == length(Group)){
      if(Normality){
        return(oneway.test(formula = Data ~ Group, var.equal = Var.Equal)$p.value)
      }else{
        return(kruskal.test(formula = Data ~ Group)$p.value)
      }
    }else{
      stop("'Data'与‘Group'应分别为等长的的numeric与character向量 ...")
    }
  }else{
    stop("'Normality'与'Var.Equal'均应为单一的logical值 ...")
  }
}


##' @description 给定分组信息, 计算各组分之间样本生存差异显著性
##' @author  Xteam.Wh
##' @param Time numeric[] 样本(患者)无复发或死亡的生存时间
##' @param Group character[] 样本分组信息, 即每个样本对应的组分
##' @param Is.Event logical[] 样本(患者)是否发生复发或死亡等截断事件
##' @param With.Weight logical 是否使用加权检验(Gehan Wilcoxon); 默认FALSE, 即使用标准的Log Rank检验
##' @return numeric 差异显著性P值
Statistical.Test.Survdiff <- function(Time, Group, Is.Event, With.Weight = FALSE){
  Time <- as.numeric(Time)
  Group <- as.character(Group)
  Is.Event <- as.logical(Is.Event)
  With.Weight <- as.logical(With.Weight)
  if(length(With.Weight) == 1){
    if(length(Time) == length(Group) && length(Time) == length(Is.Event)){
      library(survival)
      Sample.Status <- Surv(Time, Is.Event) # 构建生存对象——判断样本的的生存状态
      Survdiff.Result <- survdiff(formula = Sample.Status ~ Group, rho = ifelse(With.Weight, 1, 0))
      P.Value <- 1 - pchisq(Survdiff.Result$chisq, length(Survdiff.Result$n) - 1)
      return(P.Value)
    }else{
      stop("'Time'、'Group'以及'Is.Event'应分别为等长的的numeric、character以及logical向量 ...")
    }
  }else{
    stop("'With.Weight'应为单一的logical值 ...")
  }
}