\name{Nonlin}
\alias{Nonlin}
\title{Specify a Special Nonlinear Term in a Model Formula}
\description{
  A symbolic wrapper, for use in the formula argument to
  \code{\link{gnm}}, to indicate a nonlinear term that requires special
  handling.
}
\usage{
Nonlin(functionCall, data = NULL)
}
\arguments{
  \item{functionCall}{ a call to a \dQuote{plug-in} function that provides the
    objects required to fit the special nonlinear term. }
  \item{data}{ (optional) a data.frame containing any variables in
    \bold{specified} arguments of \code{functionCall}. }
}
\details{
    Nonlinear terms that cannot be specified by \code{Mult} require special
    handling by \code{gnm}.  \code{Nonlin} provides a way of specifying
    such terms in a \code{gnm} formula.  

    To fit the special terms, \code{gnm} requires certain objects and
    functions, which must be produced by a plug-in function.  The
    \code{functionCall}
    argument of \code{Nonlin} specifies the call to this function that is
    necessary to fit the desired term.  

    The \pkg{gnm} package includes plug-in functions for multiplicative
    interactions with homogeneous effects (\code{MultHomog}) and diagonal
    reference terms (\code{Dref}).  Users may also define their own plug-in
    functions which should return a list with at least the following
    three components:
    \describe{
      \item{labels}{a character vector of labels for the
        parameters (to which \code{gnm} will prefix the call to the
	plug-in function).}
      \item{predictor}{a function which takes a vector of parameter
        estimates and returns either a vector of
        fitted values or a matrix whose columns are
        additive components of the fitted values.}
      \item{localDesignFunction}{a function which takes the
        arguments \code{coef} (a vector of parameter estimates) and 
        \code{predictor} (the result of the predictor function), and
        returns the local design matrix.}
    }
    and optionally one further component
    \describe{
        \item{start}{either a vector of default starting values
        for the parameters or a function which takes
	the number of parameters and returns a vector of default
        starting values.}
    }
    
    Any variables passed to unspecified arguments of the plug-in function will
    be automatically added to the model frame used in \code{gnm}. It is
    assumed that there will be at least one such variable, to form the
    basis of the term. The model frame can be accessed within a plug-in
    function using \code{\link{getModelFrame}}.   
    
    The \code{data} argument of \code{Nonlin} can be used to provide any
    variables required by specified arguments of the plug-in function.
    ALL the variables given in \code{data} will be added to the model frame.
}
\value{
    An object of class \code{"Nonlin"} which is a list of deparsed expressions
    representing the variables passed to unspecified arguments of the plug-in
    function, with the attributes
  \item{call}{ the call to make to the plug-in function. }
  \item{extraData}{ if supplied, the data frame passed to the \code{data}
                   argument of \code{Nonlin}. }
	       }
\references{
  Goodman, L. A. (1979) Simple Models for the Analysis of Association in
  Cross-Classifications having Ordered
  Categories. \emph{J. Am. Stat. Assoc.}, \bold{74(367)}, 537-552.
}
\author{Heather Turner}
\seealso{\code{\link{gnm}}, \code{\link{getModelFrame}},
  \code{\link{MultHomog}}, \code{\link{Dref}}, \code{\link{Mult}}}
\examples{
set.seed(1)
data(occupationalStatus)

##  Fit an association model with homogeneous row-column effects
RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               Nonlin(MultHomog(origin, destination)), family = poisson,
               data = occupationalStatus)
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
