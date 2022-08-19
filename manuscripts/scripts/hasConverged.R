# this funtion determines whether the model has been converged or not

hasConverged <- function (mm) {

  if ( !inherits(mm, "merMod")) stop("Error: must pass a lmerMod object")

  retval <- NULL

  if(is.null(unlist(mm@optinfo$conv$lme4))) {
    retval = 1
  }
  else {
    if (isSingular(mm)) {
      retval = 1
    } else {
      retval = 0
    }
  }
  return(retval)
}
# success = 1
# failed to converge = 0







#
#
# # helper function
# # Has the model converged ?
#
# hasConverged <- function (mm) {
#
#   if ( !inherits(mm, "merMod")) stop("Error: must pass a lmerMod object")
#
#   retval <- NULL
#
#   if(is.null(unlist(mm@optinfo$conv$lme4))) {
#     retval = 1
#   }
#   else {
#     if (isSingular(mm)) {
#       retval = 0
#     } else {
#       retval = -1
#     }
#   }
#   return(retval)
# }
## the source of this function is below, and I adapted it
# which returns 1 if the model converged normally ie not to a singular fit, 0 if it converges
# to a singular fit and -1 if it fails to converge.

#https://stackoverflow.com/questions/72090177/how-can-i-know-whether-the-model-is-converged-or-failed-to-converge-in-lme4-with
