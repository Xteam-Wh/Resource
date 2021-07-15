##############################函数描述##############################
# "System.Command.Run"通过R执行给定的系统终端命令
####################################################################


##' @description 通过R执行给定的系统终端命令
##' @param System.Command character 设置要执行的终端命令
##' @param Success.Message character 设置终端命令成功运行后给出的提示信息
System.Command.Run <- function(System.Command, Success.Message = NULL){
  
  # 参数判断
  System.Command <- as.character(System.Command)
  Success.Message <- as.character(Success.Message)
  if(length(System.Command) != 1){
    stop("'Command'应为单一的character值 ...")
  }
  if(length(Success.Message) > 0){
    if(length(Success.Message) != 1){
      stop("'Success.Message'应为NULL或单一的character值 ...")
    }
  }else{
    Success.Message <- ""
  }
  
  # 根据操作系统环境设置脚本文件
  Command.Script.File <- sprintf(ifelse(Sys.info()["sysname"] == "Windows", "%s/Command.Script.bat", "%s/Command.Script.sh"), getwd())
  # 将指令写入对应系统的脚本文件
  write(sprintf(ifelse(Sys.info()["sysname"] == "Windows", "@echo off\n%s", "#!/bin/sh\n%s"),  System.Command), Command.Script.File)
  # 赋予脚本文件读写以及可执行权限
  Sys.chmod(Command.Script.File)
  
  # 执行脚本
  tryCatch(
    expr = {
      message(sprintf("<<<<<<······ %s ······>>>>>>", Sys.time()))
      message(sprintf("正在执行终端命令 ...\n*-*-*-*-*-*\n%s\n*-*-*-*-*-*", System.Command))
      Command.Run.Status <- system(sprintf("\"%s\"", Command.Script.File))
    },
    error = function(Error.Message){
      stop(Error.Message, call. = FALSE)
    },
    warning = function(Warning.Message){
      warning(Warning.Message, call. = FALSE)
    },
    finally = {
      unlink(Command.Script.File, force = TRUE)
      if(exists("Command.Run.Status") && Command.Run.Status == 0){
        message("当前终端命令执行成功 ...")
        if(nchar(Success.Message) > 0){message(sprintf("Tip: %s", Success.Message))}
      }else{
        stop("当前终端命令执行失败 ...", call. = FALSE)
      }
    }
  )
}
