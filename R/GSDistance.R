##############################函数描述##############################
# "Distance.MM"计算给定的两个集合之间的Meet/Min distance(受两集合重叠基因数目的影响)
# "Distance.KAPPA"计算给定的两个集合之间的Cohen's kappa distance(除了考虑了重叠基因数目外还受背景基因数目的影响)
# "Distance.PMM"计算给定的两个集合之间的PPI-weighted Meet/Min distance(除了考虑了重叠基因数目外还考虑了基因集中的基因在PPI网络中存在互作关系的影响)
# "Distance.Free"计算给定的两个集合之间的距离(可选: Meet/Min distance、Cohen's kappa distance、PPI-weighted Meet/Min distance), 并支持基于随机抽样法对计算结果进行显著性检验
####################################################################


##' @description 计算给定的两个集合之间的Meet/Min distance(受两集合重叠基因数目的影响)
##' @author Xteam.Wh
##' @param Set.A character[] 集合A
##' @param Set.B character[] 集合B
Distance.MM <- function(Set.A, Set.B){
  Set.A <- as.character(Set.A)
  Set.B <- as.character(Set.B)
  if(length(Set.A)*length(Set.B) > 0){
    Num.A <- length(unique(Set.A))
    Num.B <- length(unique(Set.B))
    Num.AB.Overlap <- length(intersect(Set.A, Set.B))
    MM <- 1 - Num.AB.Overlap/min(Num.A, Num.B)
    return(MM)
  }else{
    stop("'Set.A'与'Set.B'均需至少包含一个元素 ...")
  }
}


##' @description 计算给定的两个集合之间的Cohen's kappa distance(除了考虑了重叠基因数目外还受背景基因数目的影响)
##' @author Xteam.Wh
##' @param Set.A character[] 集合A
##' @param Set.B character[] 集合B
##' @param Set.Bg character[] 背景集合(受至于研究背景，且需要包含Set.A与Set.B中的所有元素)
Distance.KAPPA <- function(Set.A, Set.B, Set.Bg){
  Set.A <- as.character(Set.A)
  Set.B <- as.character(Set.B)
  Set.Bg <- as.character(Set.Bg)
  if(length(Set.A)*length(Set.B)*length(Set.Bg) > 0){
    if(all(c(Set.A, Set.B) %in% Set.Bg)){
      Num.A <- length(unique(Set.A))
      Num.B <- length(unique(Set.B))
      Num.Bg <- length(unique(Set.Bg))
      Num.AB.Union <- length(union(Set.A, Set.B))
      Num.AB.Overlap <- length(intersect(Set.A, Set.B))
      Num.A.Complement <- length(setdiff(Set.Bg, Set.A))
      Num.B.Complement <- length(setdiff(Set.Bg, Set.B))
      Num.AB.Complement <- length(setdiff(Set.Bg, union(Set.A, Set.B)))
      O <- (Num.AB.Overlap + Num.AB.Complement)/Num.Bg
      E <- (Num.A*Num.B + Num.A.Complement*Num.B.Complement)/(Num.Bg^2)
      KAPPA <- 1 - (O - E)/(1 - E)
      return(KAPPA)
    }else{
      stop("'Set.A'与'Set.B'中的元素需均存在于'Set.Bg'中 ...")
    }
  }else{
    stop("'Set.A'、'Set.B'以及'Set.Bg'均需至少包含一个元素 ...")
  }
}


##' @description 计算给定的两个集合之间的PPI-weighted Meet/Min distance(除了考虑了重叠基因数目外还考虑了基因集中的基因在PPI网络中存在互作关系的影响)
##' @author Xteam.Wh
##' @param Set.A character[] 集合A
##' @param Set.B character[] 集合B
##' @param Alpha numeric 网络的重要性/质量评分, 范围(0,1], 默认1
##' @param PPI.Score.Mtr matrix|data.frame 网络的邻接矩阵(网络中要包含Set.A与Set.B中的所有元素)
Distance.PMM <- function(Set.A, Set.B, PPI.Score.Mtr, Alpha = 1){
  Set.A <- as.character(Set.A)
  Set.B <- as.character(Set.B)
  if(length(Set.A)*length(Set.B) > 0){
    PMM <- Distance.MM(Set.A, Set.B)
    Alpha <- as.numeric(Alpha)
    if(length(Alpha) == 1 && Alpha > 0 && Alpha <= 1){
      PPI.Score.Mtr <- as.data.frame(PPI.Score.Mtr)
      if(all(c(Set.A, Set.B) %in% colnames(PPI.Score.Mtr)) && all(c(Set.A, Set.B) %in% rownames(PPI.Score.Mtr))){
        if(all(PPI.Score.Mtr >= 0) && all(PPI.Score.Mtr <= 1)){
          Num.A <- length(unique(Set.A))
          Num.B <- length(unique(Set.B))
          Coefficient <- Alpha/min(Num.A, Num.B)
          w <- min(Num.A, Num.B)/(Num.A + Num.B)
          A.Not.B <- setdiff(Set.A, Set.B)
          B.Not.A <- setdiff(Set.B, Set.A)
          A.And.B <- intersect(Set.A, Set.B)
          if(length(A.Not.B) > 0){
            A.TO.B <- sum(
              sapply(A.Not.B, function(x){
                return(w*sum(PPI.Score.Mtr[x, A.And.B]) + sum(PPI.Score.Mtr[x, B.Not.A]))
              })
            )/(max(PPI.Score.Mtr)*(w*length(A.And.B) + length(B.Not.A)))
          }else{
            A.TO.B <- 0
          }
          if(length(B.Not.A) > 0){
            B.TO.A <- sum(
              sapply(B.Not.A, function(x){
                return(w*sum(PPI.Score.Mtr[x, A.And.B]) + sum(PPI.Score.Mtr[x, A.Not.B]))
              })
            )/(max(PPI.Score.Mtr)*(w*length(A.And.B) + length(A.Not.B)))
          }else{
            B.TO.A <- 0
          }
          PMM <- PMM - Coefficient*max(A.TO.B, B.TO.A)
        }else{
          stop("'PPI.Score.Mtr'中的元素均需为介于[0,1]之间的numeric值 ...")
        }
      }else{
        stop("'PPI.Score.Mtr'的行名与列名需包含'Set.A'与'Set.B'中的所有元素 ...")
      }
    }else{
      stop("'Alpha'应该是单一的介于(0,1]之间的numeric值 ...")
    }
    return(PMM)
  }else{
    stop("'Set.A'与'Set.B'均需至少包含一个元素 ...")
  }
}


##' @description 计算给定的两个集合之间的距离(可选: Meet/Min distance、Cohen's kappa distance、PPI-weighted Meet/Min distance), 并支持基于随机抽样法对计算结果进行显著性检验
##' @author Xteam.Wh
##' @param Set.A character[] 集合A
##' @param Set.B character[] 集合B
##' @param Type character 指定要计算的距离类型, 可选("MM", "KAPPA", "PMM")
##' @param Sign.Test logical 是否对计算的距离进行显著性检验(随机抽样法), 默认TRUE
##' @param Test.Num numeric 进行显著性检验的随机检验次数, 默认1000次
##' @param Set.Bg character[] 背景集合，Type = "KAPPA"时或Type = "MM"/"KAPPA"并且Sign.Test = TRUE时必须设置该项(需要包含Set.A与Set.B中的所有元素)
##' @param Alpha numeric 网络的重要性/质量评分, 范围(0,1], 默认1
##' @param PPI.Score.Mtr character[] 网络的邻接矩阵(网络中要包含Set.A与Set.B中的所有元素)
Distance.Free <- function(Set.A, Set.B, 
                          Type = c("MM", "KAPPA", "PMM"), 
                          Sign.Test = TRUE, Test.Num = 1000, 
                          Set.Bg = NULL, Alpha = 1, PPI.Score.Mtr = NULL){
  Type <- match.arg(Type)
  Distance <- switch (Type,
                      MM = Distance.MM(Set.A, Set.B),
                      KAPPA = Distance.KAPPA(Set.A, Set.B, Set.Bg),
                      PMM = Distance.PMM(Set.A, Set.B, Alpha = Alpha, PPI.Score.Mtr = PPI.Score.Mtr)
  )
  if(length(Sign.Test) == 1 && is.logical(Sign.Test)){
    if(length(Test.Num) == 1 && is.numeric(Test.Num) && Test.Num > 0){
      if(Sign.Test){
        Pool <- switch (Type,
                        MM = Set.Bg,
                        KAPPA = Set.Bg,
                        PMM = colnames(PPI.Score.Mtr)
        )
        if(all(c(Set.A, Set.B) %in% Pool)){
          Distance.Randoms <- sapply(1:Test.Num, function(seed){
            print(seed)
            set.seed(seed = seed)
            Temp.A <- sample(Pool, length(unique(Set.A)))
            Temp.B <- sample(Pool, length(unique(Set.B)))
            Distance.Random <- switch (Type,
                                       MM = Distance.MM(Temp.A, Temp.B),
                                       KAPPA = Distance.KAPPA(Temp.A, Temp.B, Set.Bg),
                                       PMM = Distance.PMM(Temp.A, Temp.B, Alpha = Alpha, PPI.Score.Mtr = PPI.Score.Mtr)
            )
            return(Distance.Random)
          })
          P.Value <- sum(Distance.Randoms <= Distance)/Test.Num
          return(
            list(
              Distance = Distance,
              P.Value = P.Value
            )
          )
        }else{
          stop(sprintf("'Set.A'与'Set.B'中的元素需均存在于'%s' ...", ifelse(Type != "PMM", "Set.Bg", "names(PPI.Score.Mtr)")))
        }
      }else{
        return(Distance)
      }
    }else{
      stop("'Test.Num'应为单一的大于0的numeric值 ...")
    }
  }else{
    stop("'Sign.Test'应为单一的logical值 ...")
  }
}