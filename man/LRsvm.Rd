\name{LRsvm}
\alias{LRsvm}
\alias{TR_LRsvm}
\alias{predict.TR_LRsvm}

\title{SVM based implementation of  procedure}
\description{SVM based implementation of  \code{\link{LeafRank}} procedure}

\usage{
LRsvm(formula, data, bestresponse,wpos=0.5, weights = NULL, \dots)
predict.TR_LRsvm(object,newdata,\dots)
}

\arguments{
\item{formula}{an object of class \code{formula}: a symbolic description of the model to be fit.}
\item{data}{a data frame in which to interpret the variables named in the formula.}
\item{bestresponse}{The response/class/label value denoting the positive examples.}
\item{wpos}{weight of positive (bestresponse) examples in the data.}
\item{weights}{Unused in the current version. An optional vector of weights.}
\item{object}{an object of class TR_LRCart.}
\item{newdata}{a dataframe.}
\item{\dots}{optional parameters for ksvm function.}
}
\details{

LRsvm  is the svm-based version of the \code{\link{LeafRank}} procedure. This version uses the \code{kernlab} package to compute a svm classifier from the data. The \code{\link{predict}} function can be used to predict label from new data. The svm parameters are identical to the \code{\link{ksvm}} function of the \code{kernlab} package.

}

\value{
\code{LRsvm} returns an object of class \code{TR_LRsvm}.
}

\examples{
data("TRdata");
## Building a LeafRank SVM model with a polynomial kernel of degree 2 
csvm <- LRsvm(class~., Gauss2D.learn, bestresponse = 1,
kernel="polydot", kpar = list(degree = 2))

}

\author{N. Baskiotis}

\references{
S. Clemencon and N. Vayatis. \emph{Tree-based ranking methods.}
 IEEE Transactions on Information Theory, 55(9):4316-4336, 200
}



