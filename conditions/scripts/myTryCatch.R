## Copy and pasted from stack over flow
## SAve the warning message
## This code is adapted from the code that was on this following website:
##https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function/24569739#24569739

#
# myTryCatch2 <- function(expr) {
#   warn <- err <- NULL
#   value <- withCallingHandlers(
#     tryCatch(expr, error=function(e) {
#       err <<- e                        # <<- is not typo
#       print("error")
#     }), warning=function(w) {
#       warn <<- w                      # <<- is not typo
#       invokeRestart("muffleWarning")
#     })
#   if(is.null(warn)){
#     warn<-'Nothing'
#   }
#   if(is.null(err)){
#     err<-'Nothing'
#   }
#
#   list(value=value, warning=warn, error=err)
# }


library(NLP)

#
# myTryCatch <- function(expr) {
#   require(NLP)
#   warn <- err <- NULL
#   value <- withCallingHandlers(
#     tryCatch(expr, error=function(e) {
#       err <<- e                        # <<- is not typo
#       print("error")
#     }), warning=function(w) {
#       warn <<- w                      # <<- is not typo
#       invokeRestart("muffleWarning")
#     })
#   if(is.null(warn)){
#     warn<-'Nothing'
#   }
#   if(is.null(err)){
#     err<-'Nothing'
#   }
#   warning=warn
#   if(!is.na(warning[["message"]][2])){
#    out<- paste0(warning[["message"]][1],"_",warning[["message"]][2])
#   }
#   out<-warning[["message"]][1]
#   NLP::as.String(out)
#
#   return(out)
# }
#

#
# myTryCatch <- function(expr) {
#   warn <- err <- NULL
#   value <- withCallingHandlers(
#     tryCatch(expr, error=function(e) {
#       err <<- e                        # <<- is not typo
#       print("error")
#     }), warning=function(w) {
#       warn <<- w                      # <<- is not typo
#       invokeRestart("muffleWarning")
#     })
#   if(is.null(warn)){
#     warn<-'Nothing'
#     out<-warn
#   }
#   if(!is.null(warn)){
#     out<-warn[["message"]]
#   }
#   if(is.null(err)){
#     err<-'Nothing'
#   }
#
#   value=value
#   warning=warn
#   error=err
#
#   require(NLP)
#   NLP::as.String(out)
#   return(out)
# }





myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e                        # <<- is not typo
      print("error")
    }), warning=function(w) {
      warn <<- w                      # <<- is not typo
      invokeRestart("muffleWarning")
    })
  if(is.null(warn)){
    warn<-'Nothing'
  }
  if(is.null(err)){
    err<-'Nothing'
  }
  warning=warn
  paste(unlist(warning),collapse="")


}

