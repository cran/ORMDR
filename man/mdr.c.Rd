\name{mdr.c}
\alias{mdr.c}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MDR : Multifactor Dimensionality Reduction}
\description{ MDR method classifies the combination of multilocus tgenotypes into high-risk and low-risk groups based on a simple comparison of the ratios of the number of cases and controls.
}
\usage{
mdr.c(dataset,colresp,cs,combi,cv.fold=10,randomize=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{ SNP data set with class variable}
  \item{colresp}{ location of class variable in dataset}
  \item{cs}{ how to code "case"  class}
  \item{combi}{ number of genes for combination}
  \item{cv.fold}{ number of fold in cross validation}
  \item{randomize}{ logical. If 'TRUE' (the default) the cross validation sets are randomized}
}

\value{
  \item{min.comb}{combination with minimum error rate in each cross validation}  
  \item{train.erate}{training error}  
  \item{test.erate}{test error}  
  \item{data}{dataset}  
  \item{best.combi}{best combination}  
}
%\references{ ~put references to the literature/web site here ~ }
\examples{
data(sample.dat)
mdr.c(sample.dat, colresp=21, cs=1, combi=2, cv.fold = 10)
}
\keyword{distribution}
