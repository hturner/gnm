\name{Nonlin}
\alias{Nonlin}
\title{Specify a Special Nonlinear Term in a Model Formula}
\description{
  A symbolic wrapper, for use in the formula argument to \code{gnm}, to
  indicate anonlinear term that requires special handling.  
}
\usage{
Nonlin(functionCall, data = NULL)
}
\arguments{
  \item{functionCall}{a call to a ``plug-in'' function that provides the
    objects required to fit the special nonlinear term.}
  \item{data}{a data.frame containing any variables required by the plug-in
    function that are ONLY passed to specified arguments by
    \code{functionCall}.}
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
    functions which should return a list with the following components:
    \describe{
      \item{start}{(optional) either a vector of starting values
        for the parameters or a function which takes
	the number of parameters and returns a vector of
        starting values.}
      \item{labels}{a character vector of labels for the
        parameters.}
      \item{predictor}{a function which takes a vector of parameter
        estimates and returns either a vector of
        fitted values or a matrix whose columns are
        additive components of the fitted values.}
      \item{localDesignFunction}{a function which takes the
        arguments \code{coef} (a vector of
        parameter estimates) and/or
        \code{predictor} (the result
        of the predictor function), and
        returns the local design matrix.}
    }
    The plug-in function is evaluated in the environment of the model frame, so
    variables from the model frame can be accessed within the plug-in function
    using parent.frame() (see \code{MultHomog} or \code{Dref}).  

    Any variables passed to unspecified arguments of the plug-in function will
    be automatically added to the model frame used in \code{gnm}.
    Any variables passed to specified arguments must be given in the
    \code{data} argument of \code{Nonlin}.  ALL the variables given in
    \code{data} will be added to the model frame.  
}
\value{
    An object of class ``Nonlin'' which is a list of deparsed expressions
    representing the variables passed to unspecified arguments of the plug-in
    function, with the attributes
  \item{call}{the call to make to the plug-in function}
  \item{extraData}{if supplied, the data frame passed to the \code{data}
                   argument of \code{Nonlin}}
}
\author{Heather Turner}
\seealso{\code{\link{gnm}}, \code{\link{MultHomog}}, \code{\link{Dref}},
         \code{\link{Mult}}}
\examples{
data(occupational.status)
test <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
            Nonlin(MultHomog(origin, destination)), family = poisson,
            data = occupational.status)
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
