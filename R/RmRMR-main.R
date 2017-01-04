##############################################################################
#RmRMR package main functions
##############################################################################

#' An implementation of the mRMR feature selection freamework for the purpose
#' of DISESOR project. In this version the random probe test is used as the main
#' stopping criteria.
#'
#' @title Feature Subset Selection using mRMR
#' @author Andrzej Janusz
#'
#' @param dataT a data table in \code{data.table} format. Columns of the table shoud
#'              correspond to features (attributes) and rows shoul represent cases (objects).
#' @param target a vector of target values for data in \code{dataT}. For the default
#'               dependency function it is assumed that values of \code{target} are numeric,
#'               integer or binary.
#' @param dependencyF a function for computing dependencies between attributes and between
#'                    attributes and the decisions. The default is an absolute value of
#'                    Pearson's correlation (function \code{corrDependency}).
#' @param randomnessTest a function implementing a randomness test used as a stopping criteria.
#'                       The default (\code{permutationTest}) is a permutation test based on 
#'                       random probes.
#' @param Nprobes an integer specifying the number of probes to use in estimation of
#'                attribute irrelevance (for stopping criteria). The default is \code{1000}.
#' @param allowedRandomness a numeric value specifying allowed attribute irrelevance
#'                          probability. The default is \code{0.01}.
#' @param nMax an integer specifying maximal number of features that can be returned by
#'             the function. The default is \code{20}.
#' @param nCores an integer specifying the number of available processor cores for parallel
#'               computations using forking. The default is \code{1} for compatibility with
#'               Windows systems.
#' @param ... optional arguments (currently omitted).
#'
#' @return an integer vector representing indexes of attributes from the selected subset.
#'
#' @references
#' Hanchuan Peng, Fuhui Long, and Chris Ding. Feature selection based
#' on mutual information: Criteria of max-dependency, max-relevance, and
#' min-redundancy. IEEE Trans. Pattern Anal. Mach. Intell., 27(8):1226–1238
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
#' regModel = glm(targets ~.,
#'                cbind(methaneData$methaneTraining[, mrmrAttrs, with = FALSE],
#'                      targets = methaneData$methaneTrainingLabels[, as.integer(V2 == 'warning')]),
#'                family = gaussian(link = "identity"))
#'
#' preds = predict(regModel, methaneData$methaneTest, type = "response")
#' caTools::colAUC(preds, methaneData$methaneTestLabels[, V2])
#'
#' @export
mRMRfs = function(dataT, target, 
                  dependencyF = corrDependency,
                  randomnessTest = permutationTest,
                  Nprobes = 1000, allowedRandomness = 0.01, nMax = 20,
                  nCores = 1, ...) {

  if(any(is.na(dataT)) | any(sapply(as.list(dataT), function(x) any(is.infinite(x)))))
    stop("Data table contains missing or infinite values - computations aborted.")
  
  selectedAttrsIdx = integer()
  attrScores = numeric()
  attrsIdxs = 1:ncol(dataT)
  varScores = as.numeric(dataT[, unlist(mclapply(.SD, dependencyF,
                                                 unlist(target), 
                                                 mc.cores = nCores)),
                               .SDcols=attrsIdxs])

  #print(colnames(dataT)[order(varScores, decreasing = TRUE)][1:nMax])
  #select the first attribute
  selectedAttrsIdx[1] = attrsIdxs[which.max(varScores)]
  attrsIdxs = attrsIdxs[-selectedAttrsIdx[length(selectedAttrsIdx)]]
  varScores = varScores[-selectedAttrsIdx[length(selectedAttrsIdx)]]

  #test other attributes
  endFlag = FALSE
  while(!endFlag) {
    tmpDT = dataT[, selectedAttrsIdx, with = FALSE]
    dependencyScores = dataT[,
                             mcmapply(scoringFunction,
                                      .SD, varScores,
                                      MoreArgs = list(DT = tmpDT,
                                                      dependencyF = dependencyF,
                                                      nCores = nCores),
                                      mc.cores = nCores),
                             .SDcols = attrsIdxs]
    maxScoreIdx = which.max(dependencyScores)
    
    #sprawdzić kryterium stopu (probe-y)
    probeScore = randomnessTest(unlist(dataT[sample(nrow(dataT)),
                                             attrsIdxs[maxScoreIdx],
                                             with = FALSE]),
                                unlist(target),
                                varScores[maxScoreIdx],
                                dependencyF, ...)

    if(probeScore > allowedRandomness ||
       dependencyScores[maxScoreIdx] < 0) {
      endFlag = TRUE
    } else {
      selectedAttrsIdx[length(selectedAttrsIdx) + 1] = attrsIdxs[maxScoreIdx]

      attrsIdxs = attrsIdxs[-maxScoreIdx]
      varScores = varScores[-maxScoreIdx]
      
      if(length(selectedAttrsIdx) >= nMax) endFlag = TRUE
    }
  }

  names(selectedAttrsIdx) <- colnames(dataT)[selectedAttrsIdx]

  class(selectedAttrsIdx) <- unique(c("FeatureSubset", class(selectedAttrsIdx)))
  return(selectedAttrsIdx)
}
