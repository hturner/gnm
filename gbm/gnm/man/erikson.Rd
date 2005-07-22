\name{erikson}
\alias{erikson}
\docType{data}
\title{Intergenerational Class Mobility in England/Wales, France and Sweden}
\description{
  Intergenerational class mobility among the male populations of England
  and Wales; France, and Sweden.
}
\usage{data(erikson)}
\format{
    A table of counts, with classifying factors \code{origin} (father's
    class; levels \code{I}, \code{II}, \code{III}, \code{IVa},
      \code{IVb}, \code{IVc}, \code{V/VI}, 
      \code{VIIa}, \code{VIIb}) \code{destination} (son's
    class; levels as before), and \code{country} (son's country of
    residence; levels \code{EW}, \code{F}, \code{S}).
}
\source{
  Hauser, R. M. (1984) Vertical Class Mobility in England, France and
  Sweden. \emph{Acta Sociol.}, \bold{27(2)}, 87-110.
}
\references{
  Erikson, R., GoldThorpe, J. H. and Portocarero, L. (1982) Social
  Fluidity in Industrial Nations: England, France and
  Sweden. \emph{Brit. J. Sociol.} \bold{33(1)}, 1-34.

  Xie, Y. (1992) The Log-multiplicative Layer Effect Model for Comparing
  Mobility Tables. \emph{Am. Sociol. Rev.} \bold{57(3)}, 380-395.
}
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

### Fit the models given in first half of Table 3 of Xie (1992)

## Null association between origin and destination
nullModel <- gnm(Freq ~ country:origin + country:destination,
                 family = poisson, data = erikson)

## Full interaction, common to all countries
commonInteraction <- update(nullModel, ~ . + origin:destination)

## Full Interaction, different multiplier for each country
multInteraction <- update(nullModel, ~ . + Mult(country, origin:destination))

### Create array of interaction levels as in Table 2 of Xie (1992)

levelMatrix <- matrix(c(2, 3, 4, 6, 5, 6, 6,
                        3, 3, 4, 6, 4, 5, 6,
                        4, 4, 2, 5, 5, 5, 5,
                        6, 6, 5, 1, 6, 5, 2,
                        4, 4, 5, 6, 3, 4, 5,
                        5, 4, 5, 5, 3, 3, 5,
                        6, 6, 5, 3, 5, 4, 1), 7, 7, byrow = TRUE)

### Fit models in second half of Table 3 in Xie (1992)

## Interaction specified by levelMatrix, common to all countries
commonTopo <- update(nullModel, ~ . +
                     Topo(origin, destination, spec = levelMatrix))

## Interaction specified by levelMatrix, different multiplier for each country
multTopo <- update(nullModel, ~ . +
                   Mult(country, Topo(origin, destination, spec = levelMatrix)))

## Interaction specified by levelMatrix, different effects for each country
separateTopo <- update(nullModel, ~ . +
                       country:Topo(origin, destination, spec = levelMatrix))
}
\keyword{datasets}
