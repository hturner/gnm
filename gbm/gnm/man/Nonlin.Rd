\name{Nonlin}
\alias{Nonlin}
\title{Specify a Special Nonlinear Term in a Model Formula}
\description{
  A symbolic wrapper, for use in the formula argument to
  \code{\link{gnm}}, to indicate a nonlinear term that requires special
  handling.
}
\usage{
Nonlin(functionCall)
}
\arguments{
  \item{functionCall}{ a call to a plug-in function that provides the
    objects required to fit the special nonlinear term. }
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
    reference terms (\code{Dref}).

    Users may also define their own plug-in functions which should
    return a list with at least the following three components:
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
        returns the local design matrix. If the plug-in function does
	not return a 'start' component, this function must also take the
	argument \code{ind}, the index of a column to be returned
	rather than the full matrix.}  
    }
    and optionally one further component
    \describe{
        \item{start}{a vector of default starting values for the
	  parameters. \code{NA} may be used to indicate parameters which
	  may be treated as linear for the purpose of finding starting
	  values, given the non-\code{NA} values. } 
      }

    \code{Nonlin} identifies variables in the term from the call to the
    plug-in function and returns deparsed expressions representing these
    variables, so that they can be added to the model frame. By default,
    expressions passed to unspecified arguments of the plug-in function
    are deparsed.

    If the default action of \code{Nonlin} will not capture the required
    variables, a companion function must exist in the environment of the
    plug-in function, which takes a call to the plug-in function and
    returns the necessary deparsed expressions. The name of this
    function must be the name of the plug-in function suffixed with
    \code{Variables}, e.g. \code{PlugInFunctionVariables} would be the
    companion function for \code{PlugInFunction}. 

    The call to the plug-in function is evaluated in the environment of
    the model frame and in the enclosing environment of the parent frame
    of the call to \code{gnm}. To access the model frame directly from
    within a plug-in function (e.g. for use with \code{model.matrix}),
    use \code{\link{getModelFrame}}.
}
\value{
    An object of class \code{"Nonlin"} which is a list of deparsed expressions
    representing the variables in the nonlinear term, with the call to
    the plug-in function as an attribute.  
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
