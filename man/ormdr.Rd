\name{ormdr}
\alias{ormdr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ORMDR : Odds ratio based multifactor-dimensioinality reduction method}
\description{ORMDR usesthe odds ratio as a new quantitative measure of disease risk.
}
\usage{
ormdr(dataset,bestcombi,cs,colresp,CI.Asy,CI.Boot,B)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{ SNP data set with class variable}
  \item{bestcombi}{best combination from mdr.c function}
  \item{cs}{ how to code "case"  class}
  \item{colresp}{location of class variable in dataset}
  \item{CI.Asy}{indicator whether asymptotic CI is calculated or not}
  \item{CI.Boot}{indicator whether bootstrap CI is calculated or not}
  \item{B}{number of bootstrap samples}
}

\value{
  ...
}
%\references{ ~put references to the literature/web site here ~ }
\examples{
data(sample.dat)
z<-mdr.c(sample.dat, colresp=21, cs=1, combi=2, cv.fold = 10)
ormdr(sample.dat,bestcombi=as.numeric(z$best.combi),cs=1,colresp=21,CI.Asy=TRUE,CI.Boot=TRUE,B=1000)

}
\keyword{distribution}
