##############################################################################
#RmRMR package auxiliary functions
##############################################################################

# An auxiliary function for computing inter-attribute dependencies
scoringFunction = function(testedAttr, score, DT, dependencyF, nCores = 1) {

  varDependencyScores = as.numeric(DT[,
                                      unlist(
                                      mclapply(.SD, dependencyF,
                                               testedAttr, mc.cores = nCores)),
                                      .SDcols=1:ncol(DT)])
  return(score - max(varDependencyScores))
}

#' This is an exemplary function for measuring inter-attribute dependency
#' and dependency between attributes and the dacisions. It is used by default
#' by the provided implementation of the mRMR framework.
#'
#' @title Computation of indiscernibility classes based on the rough set theory
#' @author Andrzej Janusz
#'
#' @param x,y two numeric attributes or an attribute and a numeric decision
#'            (binary and integer valued decisions are also accepted).
#' @param ... optional arguments (currently omitted)
#'
#' @return a numeric value expressing linear dependency between \code{x} and \code{y}
#'
#' @references
#' Mark Hall. Correlation-based Feature Selection for Machine Learning. PhD thesis,
#' University of Waikato, 1999.
#'
#' @examples
#' #############################################
#' data(methaneSampleData)
#'
#' ## an experiment on a sample from methane competition data (IJCRS'16)
#' mrmrAttrs = mRMRfs(dataT = methaneData$methaneTraining,
#'                    target = methaneData$methaneTrainingLabels[, as.integer(V2 == 'warning')],
#'                    dependencyF = corrDependency)
#'
#' mrmrAttrs
#'
#' @export
corrDependency = function(x,y, ...) abs(cor(x,y))
