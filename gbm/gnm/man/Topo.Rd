\name{Topo}
\alias{Topo}
\title{ Topological Interaction of Factors }
\description{
  Given two or more factors \code{Topo} creates an interaction factor
  as specified by an array of levels, which may be arbitrarily
  structured.
}
\usage{
Topo(..., spec = NULL)
}
\arguments{
  \item{\dots}{ two or more factors }
  \item{spec}{ an array of levels, with dimensions corresponding to the
    number of levels of each factor in the interaction }
}
\value{
  A factor of levels extracted from the levels array given in
  \code{spec}, using the given factors as index variables.
}
\references{
  Erikson, R., GoldThorpe, J. H. and Portocarero, L. (1982) Social
  Fluidity in Industrial Nations: England, France and
  Sweden. \emph{Brit. J. Sociol.} \bold{33(1)}, 1-34.

  Xie, Y. (1992) The Log-multiplicative Layer Effect Model for Comparing
  Mobility Tables. \emph{Am. Sociol. Rev.} \bold{57(3)}, 380-395.
}
\author{ David Firth }
\seealso{ \code{\link{Symm}} and \code{\link{Diag}} for special cases }
\examples{
set.seed(1)
data(erikson)

### Collapse to 7 by 7 table as in Erikson (1982)

erikson <- as.data.frame(erikson)
lvl <- levels(erikson$origin)
levels(erikson$origin) <- levels(erikson$destination) <-
    c(rep(paste(lvl[1:2], collapse = " + "), 2), lvl[3],
      rep(paste(lvl[4:5], collapse = " + "), 2), lvl[6:9])
erikson <- xtabs(Freq ~ origin + destination + country, data = erikson)

### Create array of interaction levels as in Table 2 of Xie (1992)

levelMatrix <- matrix(c(2, 3, 4, 6, 5, 6, 6,
                        3, 3, 4, 6, 4, 5, 6,
                        4, 4, 2, 5, 5, 5, 5,
                        6, 6, 5, 1, 6, 5, 2,
                        4, 4, 5, 6, 3, 4, 5,
                        5, 4, 5, 5, 3, 3, 5,
                        6, 6, 5, 3, 5, 4, 1), 7, 7, byrow = TRUE)

### Fit the levels models given in Table 3 of Xie (1992)

## Null association between origin and destination
nullModel <- gnm(Freq ~ country:origin + country:destination,
                 family = poisson, data = erikson)

## Interaction specified by levelMatrix, common to all countries
commonTopo <- update(nullModel, ~ . +
                     Topo(origin, destination, spec = levelMatrix))

## Interaction specified by levelMatrix, different multiplier for
## each country
multTopo <- update(nullModel, ~ . +
                   Mult(Exp(-country),
                        Topo(origin, destination, spec = levelMatrix)))

## Interaction specified by levelMatrix, different effects for
## each country
separateTopo <- update(nullModel, ~ . +
                       country:Topo(origin, destination,
                                    spec = levelMatrix))
}
\keyword{ models }
