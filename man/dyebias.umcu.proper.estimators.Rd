\name{dyebias.umcu.proper.estimators}
\alias{dyebias.umcu.proper.estimators}

\title{ Determine which spots should not be ruled out as slide bias estimators }
\description{
  Some spots (reporters/probes) should not be used when estimating the slide
  bias. Typical examples are mitochondrial genes and spots known to
  cross-hybridize. This function finds the ones that are OK to use. 
}

\synopsis{
  dyebias.umcu.proper.estimators(reporter.info, verbose = FALSE)
}

\arguments{

  \item{reporter.info}{
    A data.frame, one row per spot, with (at least) columns \code{reporterId}
    (e.g. gene id or oligo id) and any of the following
    characteristics: \code{reporterGroup, chromosomeName, bioSeqType,
    crosshybRank and reporterSequence}.  They are used to get
    rid of reporters that are not suitable when estimating the slide
    bias.
  }

  \item{verbose}{ Logical speficying whether to be verbose or not }
}

\details{

  This function is particular to the slides and database set-up at the
  Holstege lab, but might serve as inspiration.
}

\value{
  Returns and index vector that can be
  used as the \code{estimator.subset}-argument to
  \code{dyebias.application.subset}.
}

\references{
  Margaritis, T., Lijnzaad, P., van~Leenen, D., Bouwmeester, D.,
  Kemmeren, P., van~Hooff, S.R and Holstege, F.C.P. (2009)
  Adaptable gene-specific dye bias correction for two-channel DNA microarrays.
  \emph{Molecular Systems Biology, submitted}
}

\author{Philip Lijnzaad \email{p.lijnzaad@umcutrecht.nl} }

\seealso{
  \code{\link[dyebias]{dyebias.apply.correction}}
}

\examples{

  \dontshow{
     library(dyebias)
     library(dyebiasexamples)
     data(data.raw)
     data(data.norm)

     ### obtain estimate for the iGSDBs:
     iGSDBs.estimated <- dyebias.estimate.iGSDBs(data.norm,
                                                 is.balanced=TRUE,
                                                 verbose=TRUE)

  }

  ### choose the estimators and which spots to correct:
  estimator.subset <- dyebias.umcu.proper.estimators(maInfo(maGnames(data.norm)))

  summary(estimator.subset)

  ### do the correction
  \dontrun{
     correction <- dyebias.apply.correction(data.norm=data.norm,
                                         iGSDBs = iGSDBs.estimated,
                                         estimator.subset=estimator.subset,
                                         application.subset = TRUE,
                                         verbose=TRUE)
   }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}                          % silly but compulsary
