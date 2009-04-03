\name{data.raw}
\alias{data.raw}
\alias{data.norm}

\docType{data}

\title{Example data for the dyebias package}

\description{

  The \code{dyebias}-package, described in Margaritis et al. (2009) can
  be used to get rid of dye bias in two-colour microarrays. The
  \code{data.raw} and \code{data.norm} objects are used in its examples.

  The objects represent four hybridizations of identical mRNA, with
  increasing Cy3 and Cy5 labeling percentages (identical per slide) and
  differently spiked-in external controls to judge the process of
  dyebias correction.

}

\usage{
       data(data.raw)
       data(data.norm)
}

\format{
       
  The data uses the \code{marray}-package by Dudoit and Yang (2002).
  \code{data.raw} is a \code{marrayRaw} object, \code{data.norm} is a
  \code{marrayNorm} object derived from it by print-tip LOESS
  normalization. Neither is dyebias-corrected yet. 
  
}

\source{

  All accession numbers below refer to ArrayExpress (http://www.ebi.ac.uk/microarray).

  This two-colour microarrray data was obtained from identical mRNA
  extracts (protocol P-UMCU-37), spiked with external controls,
  dUTP-labeled with Cy3 and Cy5 (protocol P-UMCU-38). This was
  hybridized (protocol P-UMCU-39) onto self-spotted slides containing
  70-mer oligonucleotides (2 replicates per oligo, Operon "Array-Ready",
  and including 2838 control features; protocol P-UMCU-34).  Scanning
  was done with an Agilent G2565AA scanner (protocol P-UMCU-40) and
  images were quantified with BioDiscovery's ImaGene 7.x (protocol
  P-UMCU-42)

}
  
\details{
  The column \code{R.group} of \code{maInfo(maTargets(data.norm))}
  shows the details. Eg.,  \code{4\%_2EC} indicates that the labeling
  (of both channels) was at 4\%, and the external controls were spiked
  in at a concentration twice that of the green channel. See Margaritis
  et~al. (2009) for details.
}

\examples{
  data(data.raw)
  data(data.norm)
}

\note{
  The Tuteja data is also included in this package under the
  (inst)/doc directory, as this data is not proper rda, tab or csv
  data. For details, refer to the original publication and/or the
  \code{dyebias} vignette.
}


\author{Philip Lijnzaad}

\references{

  Margaritis, T., Lijnzaad, P., van~Leenen, D., Bouwmeester, D.,
  Kemmeren, P., van~Hooff, S.R and Holstege, F.C.P. (2009).
  Adaptable gene-specific dye bias correction for two-channel DNA microarrays.
  \emph{Molecular Systems Biology, submitted}

  Dudoit, S. and Yang, Y.H. (2002)  Bioconductor R packages for
  exploratory analysis and normalization of cDNA microarray data. In:
  Parmigiani, G., Garrett, E.S. , Irizarry, R.A., and Zeger, S.L. (eds.)
  \emph{The Analysis of Gene Expression Data: Methods and Software},
  New~York: Springer


}

\keyword{datasets}
