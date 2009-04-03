dyebias.umcu.proper.estimators <- function
(reporter.info,                          
 verbose=FALSE
 ) {
  here <- "dyebias.umcu.proper.estimators"
  
### any columns not available must be ignored, i.e. have to be TRUE:
  some.missing <- FALSE
  subset <- TRUE

  if( is.null(reporter.info$reporterGroup )) {
    some.missing <- TRUE
  } else {
    subset <-  subset & (reporter.info$reporterGroup == "gene")
  }

  if( is.null(reporter.info$chromosomeName)) {
    some.missing <- TRUE
  } else {
    subset <-  subset & (reporter.info$chromosomeName != "MT")
  }
  
  if( is.null(reporter.info$bioSeqType)) {
    some.missing <- TRUE
  } else {
    subset <-  subset & (reporter.info$bioSeqType != "transposable_element_gene")
  }

  if( is.null(reporter.info$crosshybRank)) {
    some.missing <- TRUE
  } else {
    subset <-  subset & (reporter.info$crosshybRank == 1)
  }

  if( is.null(reporter.info$reporterSequence) ) {
    some.missing <- TRUE
  } else {
    subset <-  subset & (!duplicated(reporter.info$reporterSequence))
  }
  
  if(verbose  & some.missing) {
    warning("reporter.info is missing >=1 of the columns\n\
reporterGroup, chromosomeName, bioSeqType, crosshybRank, reporterSequence\n\
as a result, unsuitable reporters maybe included in the estimate, in", here, call.=TRUE)
  }

  if (length(subset)==0 || sum(subset)==0) {
    if(verbose) {
      warning("Could not find suitable genes to base dyebias extremes, so using all of them, in", here, call.=TRUE)
    }
    return(TRUE)                        #i.e. all
  }
  return(subset)
}                                       # dyebias.umcu.proper.estimators
