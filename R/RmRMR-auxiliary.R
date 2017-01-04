##############################################################################
#RmRMR package auxiliary functions
##############################################################################

# An auxiliary function for computing inter-attribute dependencies
scoringFunction = function(testedAttr, score, DT, dependencyF, nCores = 1) {

  varDependencyScores = as.numeric(DT[,
                                      unlist(
                                      mclapply(.SD, dependencyF,
                                               unlist(testedAttr), 
                                               mc.cores = nCores)),
                                      .SDcols=1:ncol(DT)])
  return(score - max(varDependencyScores))
}

#' An auxiliary function for performing a permutation test of randomness
#' 
#' @title Permutation test of randomness
#' @author Andrzej Janusz
#' 
#' @param testedAttr tested attribute
#' @param target the target attribute
#' @param attrScore tested attribute's base score value
#' @param dependencyF a dependency function used by the test
#' @param Nprobes a number of generated random probes
#' @param nCores an integer specifying the number of available processor cores for parallel
#'               computations using forking. The default is \code{1} for compatibility with
#'               Windows systems.
#' @param ... additional parameters (currently omitted)
#' 
#' @return a numeric value - an estimation of a probability that a random probe could obtain
#'         a score equal or greater than \code{attrScore}.
#'         
#' @references
#' XY
#' 
#' @examples
#' #############################################
#' 
#' @export
permutationTest = function(testedAttr, target, attrScore, dependencyF, Nprobes = 1000, nCores = 1, ...) {
  
  probeScores = mclapply(1:Nprobes,
                         function(x) {y = sample(testedAttr); dependencyF(y, unlist(target))},
                         mc.cores = nCores)
  
  return(mean(unlist(probeScores) >= attrScore))
}

#' This is an exemplary function for measuring inter-attribute dependency
#' and dependency between attributes and the dacisions. It is used by default
#' by the provided implementation of the mRMR framework.
#'
#' @title Linear correlation-based dependency measure
#' @author Andrzej Janusz
#'
#' @param x,y two numeric attributes or an attribute and a numeric decision
#'            (binary and integer valued decisions are also accepted).
#' @param ... optional arguments to the \code{cor} function.
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
#' ## an experiment on a sample from the data used in a data mining competition - 
#' ## IJCRS'15 Data Challenge: Mining Data from Coal Mines 
#' ## (https://knowledgepit.fedcsis.org/contest/view.php?id=109).
#' ## The whole data set can be downloaded from the competition web page.
#' 
#' mrmrAttrs = mRMRfs(dataT = methaneData$methaneTraining,
#'                    target = methaneData$methaneTrainingLabels[, as.integer(V2 == 'warning')],
#'                    dependencyF = corrDependency)
#'
#' mrmrAttrs
#'
#' @export
corrDependency = function(x,y, ...) abs(cor(x, y, ...))


#' A function for measuring ordering dependency between two vectors (the ordering measure).
#'
#' @title Ordering-based dependency measure
#' @author Andrzej Janusz
#'
#' @param attrVec,target two numeric attributes or an attribute and a decision vector
#'            (not necessarily ordered). At least the first argument has to be numeric.
#' @param numericTarget logical indicating whether the \code{target} is ordered and numeric. The 
#'            default is \code{NULL} in which case the function makes a guess using a simple heuristic.            
#' @param ... optional arguments (currently omitted).
#'
#' @return a numeric value expressing ordering dependency between \code{attrVec} and \code{target}
#'
#' @references
#' Andrzej Janusz and Marek Grzegorowski. Efficient Attribute Quality Assessment Using 
#' a Decision Ordering Measure. 2017.
#'
#' @examples
#' #############################################
#' data(methaneSampleData)
#'
#' ## an experiment on a sample from the data used in a data mining competition - 
#' ## IJCRS'15 Data Challenge: Mining Data from Coal Mines 
#' ## (https://knowledgepit.fedcsis.org/contest/view.php?id=109).
#' ## The whole data set can be downloaded from the competition web page.
#' 
#' mrmrAttrs = mRMRfs(dataT = methaneData$methaneTraining,
#'                    target = methaneData$methaneTrainingLabels[, V2],
#'                    dependencyF = orderingDependency)
#'
#' mrmrAttrs
#'
#' @export
orderingDependency <- function(attrVec, target, numericTarget = NULL, ...)  {
  
  tmpIdx = order(unlist(attrVec))
  if(is.null(numericTarget)) {
    numericTarget = (is.numeric(target) & (length(unique(target)) > 3)) | is.ordered(target)
  }
  
  if(!numericTarget) {
    orderingMeasureValue = .C("computeSwitchesC",
                              target = as.double(factor(target[tmpIdx])),
                              N = as.integer(length(target)),
                              value = as.integer(0))
  } else {
    tmp = target[tmpIdx]
    target = as.integer((tmp[2:length(tmp)] - tmp[1:(length(tmp) - 1)]) >= 0)
    orderingMeasureValue = .C("computeSwitchesC",
                              target = as.double(target),
                              N = as.integer(length(target)),
                              value = as.integer(0))
  }
  
  return((length(target) - orderingMeasureValue$value)/length(target))
}

#' A balanced version of the function for measuring ordering dependency between two vectors 
#' (the ordering measure).
#'
#' @title Balanced ordering-based dependency measure
#' @author Andrzej Janusz
#'
#' @param attrVec,target two numeric attributes or an attribute and a decision vector
#'            (not necessarily ordered). At least the first argument has to be numeric.
#' @param numericTarget logical indicating whether the \code{target} is ordered and numeric. The 
#'            default is \code{NULL} in which case the function makes a guess using a simple heuristic.            
#' @param ... optional arguments (currently omitted).
#'
#' @return a numeric value expressing ordering dependency between \code{attrVec} and \code{target}
#'
#' @references
#' Andrzej Janusz and Marek Grzegorowski. Efficient Attribute Quality Assessment Using 
#' a Decision Ordering Measure. 2017.
#'
#' @examples
#' #############################################
#' data(methaneSampleData)
#'
#' ## an experiment on a sample from the data used in a data mining competition - 
#' ## IJCRS'15 Data Challenge: Mining Data from Coal Mines 
#' ## (https://knowledgepit.fedcsis.org/contest/view.php?id=109).
#' ## The whole data set can be downloaded from the competition web page.
#' 
#' mrmrAttrs = mRMRfs(dataT = methaneData$methaneTraining,
#'                    target = methaneData$methaneTrainingLabels[, V2],
#'                    dependencyF = balancedOrderingDependency, sampleSize = 10)
#'
#' mrmrAttrs
#'
#' @export
balancedOrderingDependency <- function(attrVec, target, numericTarget = NULL, nSamples = 100, sampleSize = 10, ...)  {
  
  #tmpIdx = order(unlist(attrVec))
  if(is.null(numericTarget)) {
    numericTarget = (is.numeric(target) & (length(unique(target)) > 3)) | is.ordered(target)
  }
  
  if(!numericTarget) {
    decisionClasses = split(1:length(target), target)
    samples = replicate(nSamples, 
                        {sample(unlist(lapply(decisionClasses, function(x) sample(x, min(length(x), sampleSize)))))},
                        simplify = FALSE)
    score = mean(sapply(samples, 
                        function(x) orderingDependency(attrVec[x], target[x], 
                                                       numericTarget = FALSE)))
    orderingMeasureValue = list(value = length(target) - score*length(target))
  } else {
    tmpIdx = order(unlist(attrVec))
    tmp = target[tmpIdx]
    target = as.integer((tmp[2:length(tmp)] - tmp[1:(length(tmp) - 1)]) >= 0)
    orderingMeasureValue = .C("computeSwitchesC",
                              target = as.double(target),
                              N = as.integer(length(target)),
                              value = as.integer(0))
  }
  
  return((length(target) - orderingMeasureValue$value)/length(target))
}
