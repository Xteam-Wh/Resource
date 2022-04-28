##############################函数描述##############################
# "Feature.Selection.Stepwise"使用逐步回归的三种方法通过AIC筛选特征
# "Feature.Selection.Cox"通过单/多因素COX回归分析, 筛选与预后相关的特征
####################################################################


##' @description 使用逐步回归的三种方法通过AIC筛选特征
##' @author  Xteam.Wh
##' @param Data matrix | data.frame 特征矩阵(行为样本列为特征,最后一列为因变量)
##' @param Is.Reg logical 是否为回归模型筛选特征, 否则所选特征应用于分类模型, 默认True
##' @return 逐步回归三种方法筛选的特征集合
Feature.Selection.Stepwise <- function(Data, Is.Reg = TRUE){
  library(MASS)
  Data <- as.data.frame(Data)
  Category <- colnames(Data)[ncol(Data)] # 获取因变量在特征矩阵中的名称
  Is.Reg <- as.logical(Is.Reg)
  if(length(Is.Reg) == 1){
    if(Is.Reg){
      Null.Model <- lm(as.formula(paste0(Category," ~ NULL")), data = Data) # 构建空特征线性回归模型
      Full.Model <- lm(as.formula(paste0(Category," ~ .")), data = Data)# 构建全特征线性回归模型
    }else{
      Null.Model <- glm(as.formula(paste0(Category," ~ NULL")), family = binomial(link='logit'), data = Data) # 构建空特征罗杰斯特回归模型
      Full.Model <- glm(as.formula(paste0(Category," ~ .")), family = binomial(link='logit'), data = Data)# 构建全特征罗杰斯特回归模型
    }
  }else{
    stop("'Is.Reg'应为单一的logical值")
  }
  # 逐步回归进行特征选择
  Forward.Model <- stepAIC(Null.Model, trace = FALSE, direction = "forward", scope = list(lower=Null.Model, upper=Full.Model))
  Backward.Model <- stepAIC(Full.Model, trace = FALSE, direction = "backward")
  Both.Model <- stepAIC(Full.Model, trace = FALSE, direction = "both")
  # 返回三种方式筛选出的特征
  return(
    list(
      Forward.Features = names(Forward.Model$coefficients)[-1],
      Backward.Features = names(Backward.Model$coefficients)[-1],
      Both.Features = names(Both.Model$coefficients)[-1]
    )
  )
}


##' @description 通过单/多因素COX回归分析, 筛选与预后相关的特征
##' @author  Xteam.Wh
##' @param Data matrix | data.frame 特征矩阵(行为样本列为特征, 无因变量)
##' @param Time numeric[] 样本(患者)无复发或死亡的生存时间
##' @param Is.Event logical[] 样本(患者)是否发生复发或死亡等截断事件
##' @param Ties.Method character 指定估计生存率所用的方法, 可选("efron","breslow","exact"); 默认"efron"
##' @param Analysis.Mode character 指定分析模式(单因素/多因素), 可选("Single", "Multiple"); 默认"Single"
##' @return data.frame 包含所有特征的信息, 包含"Characteristics、HR、%95 CI、P.Value"
Feature.Selection.Cox <- function(Data, Time, Is.Event, 
                                  Ties.Method = c("efron","breslow","exact"), Analysis.Mode = c("Single", "Multiple")){
  Data <- as.data.frame(Data)
  Time <- as.numeric(Time)
  Is.Event <- as.logical(Is.Event)
  if(nrow(Data) == length(Time) && nrow(Data) == length(Is.Event)){
    library(survival)
    Features <- colnames(Data) # 要进行筛选的特征
    Sample.Status <- Surv(Time, Is.Event) # 构建生存对象——判断样本的的生存状态
    Ties.Method <- match.arg(Ties.Method)
    Analysis.Mode <- match.arg(Analysis.Mode)
    # 对各特征进行单因素分析
    Features.Cox.Result <- switch(Analysis.Mode,
                                  Single = as.data.frame(
                                    t(
                                      sapply(Features, function(Feature){
                                        # 进行单因素cox回归分析, 并获取分析结果的统计信息
                                        Formula <- as.formula(sprintf("Sample.Status ~ `%s`", Feature))
                                        Model <- coxph(Formula, data = Data, ties = Ties.Method)
                                        Statistics <- summary(Model)
                                        return(
                                          c(
                                            HR = Statistics$coefficients[1, "exp(coef)"],
                                            `Lower.%95CI` = Statistics$conf.int[1, "lower .95"],
                                            `Upper.%95CI` = Statistics$conf.int[1, "upper .95"],
                                            P.Value = Statistics$coefficients[1, "Pr(>|z|)"]
                                          )
                                        )
                                      })
                                    )
                                  ),
                                  Multiple = {
                                    # 进行多因素cox回归分析, 并获取分析结果的统计信息
                                    Formula <- as.formula(sprintf("Sample.Status ~ `%s`", paste0(Features, collapse = "` + `")))
                                    Model <- coxph(Formula, data = Data, ties = Ties.Method)
                                    Statistics <- summary(Model)
                                    # 获取风险比(HR-hazard ratio)以、%95置信区间以及显著性检验的P值
                                    return(
                                      data.frame(
                                        HR = Statistics$coefficients[, "exp(coef)"],
                                        `Lower.%95CI` = Statistics$conf.int[, "lower .95"],
                                        `Upper.%95CI` = Statistics$conf.int[, "upper .95"],
                                        P.Value = Statistics$coefficients[, "Pr(>|z|)"],
                                        check.names = FALSE
                                      )
                                    )
                                  }
    )
    # 返回结果
    return(Features.Cox.Result)
  }else{
    stop("'Time'与'Is.Event'应分别为与'Data'行数相等的numeric向量以及logical向量 ...")
  }
}