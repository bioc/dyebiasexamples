\name{dyebias.geo2marray}
\alias{dyebias.geo2marray}

\title{ convenience function to convert GEO objects to marray objects }

\description{
  convenience function to convert GEO objects to marray objects
}

\synopsis{

dyebias.geo2marray(gse, slide.ids = NULL, type = "norm",
                   gene.selector = function(table) {TRUE},
                   reporterid.name, cy3.name, cy5.name,
                   R.name = NULL, G.name = NULL, M.name = NULL, A.name =NULL, 
                   Rf.name = NULL, Gf.name = NULL, Rb.name = NULL, Gb.name = NULL)
}

\arguments{

  \item{gse}{ GSE data set }

  \item{slide.ids}{ Return only the slides with these ids. If NULL,
    return all.}

\item{type}{ what to extract; must be either "norm" or "raw". }

\item{gene.selector}{ function(table) acting on Table(GPL) giving back 
                      an index with the rows considered to be genes. }

\item{reporterid.name}{column containing the reporter.id, in Table(gpl).}

\item{cy3.name}{ The column name containing the factor value for the Cy3
  (green) channel }

\item{cy5.name}{ The column name containing the factor value for the Cy5
  (red) channel }

 \item{R.name}{ column name for extracting the R data from Table(gsm) }

  \item{G.name}{ column name for extracting the G data from Table(gsm) }

  \item{M.name}{ column name for extracting the M data from Table(gsm) }

  \item{A.name}{ column name for extracting the A data from Table(gsm) }

  \item{Rf.name}{ column name for extracting the Rf data from Table(gsm) }

  \item{Gf.name}{ column name for extracting the Gf data from Table(gsm) }

  \item{Rb.name}{ column name for extracting the Rb data from Table(gsm) }

  \item{Gb.name}{ column name for extracting the Gb data from Table(gsm) }
}

\details{
  The \code{XYZ.name} mechanism is the same as that used in 
\code{\link[marray]{read.marrayRaw}}; i.e. you specify the name of the
  column that contains the desired data.
}

\value{
  A full-fledged  marrayRaw (if
  \code{type} was "raw") or marrayNorm (if \code{type} was "norm") is returned.
  
}


\author{Philip Lijnzaad}

\note{At some point, this functionality should be merged into the
  convert package.
}

\examples{

\dontshow{
  cat("Running this example takes too much time; if you want that, see
the second example in the vignette\n")
}                                       % dontshow

\dontrun{
  ## Running this example takes too much time; if you want that, see the
  ## second example in the vignette
}                                       % dontrun

}

\references{

  Davis, S. and Meltzer, P.S (2007). GEOquery: a bridge
  between the Gene Expression Omnibus (GEO) and BioConductor.
  \emph{Bioinformatics} 23, 1846--1847 (doi:10.1093/bioinformatics/btm254).

  Dudoit, S. and Yang, Y.H. (2002)  Bioconductor R packages for
  exploratory analysis and normalization of cDNA microarray data. In:
  Parmigiani, G., Garrett, E.S. , Irizarry, R.A., and Zeger, S.L. (eds.)
  \emph{The Analysis of Gene Expression Data: Methods and Software},
  New~York: Springer

  Chen,S., de~Vries, M.A. and Bell, S.P. (2007) \emph{Genes Dev.} 21, 2897--2907
  "Orc6 is required for dynamic recruitment of Cdt1 during repeated Mcm2-7 loading"
  (doi:10.1101/gad.1596807)

}

\keyword{misc}                          % silly but compulsary


















