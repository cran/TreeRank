\name{varImportance}
\alias{varImportance}
\alias{varImportance.TR_LRCart}
\alias{varImportance.TR_LRforest}
\alias{varImportance.TR_forest}
\alias{varImportance.TR_TreeRank}
\alias{varDep}

\title{Variable Importance and Variable Dependance analysis}
\description{Variable Importance and Dependance functions for Cart- and Random Forest-based LeafRank models}
\usage{
varImportance(obj,norm=TRUE)
varDep(obj,data,varx,vary,vminx=min(data[varx]),vmaxx=max(data[varx]),vminy=min(data[vary]),vmaxy=max(data[vary]),subdivx=100,subdivy=subdivx)
}

\arguments{
\item{obj}{an object of \code{\link{TR_LRCart}}, \code{\link{TR_LRforest}},  \code{\link{TR_forest}} or  \code{\link{TR_TreeRank}} class.}
\item{norm}{a boolean indicating when the result is normalized.}
\item{data}{a dataframe.}
\item{varx}{the variable name of the first variable.}
\item{vary}{the variable name of the second variable.}
\item{vminx}{the minimal value of the first variable.}
\item{vmaxx}{the maximal value of the first variable.}
\item{vminy}{the minimal value of the second variable.}
\item{vmaxy}{the maximal value of the second variable.}
\item{subdivx}{the number of subdivision for the range of the first variable.}
\item{subdivy}{the number of subdivision for the range of the second variable.}
}

\value{

\code{varImportance} returns an array of length the number of model variables with the importance measure of each variable.
\code{varDep} returns a matrix with the variable dependance measures.
}


\examples{
data("TRdata");

## Building a LeafRank Random Forest model
cmod <- LRforest(class~., Gauss2D.learn, bestresponse = 1)

## Importance measures computation
vi <- varImportance(cmod);

}

\author{N. Baskiotis}

\references{
S. Clemencon and N. Vayatis. \emph{Tree-based ranking methods.}
 IEEE Transactions on Information Theory, 55(9):4316-4336, 200
}


