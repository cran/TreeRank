\name{LeafRank}
\alias{LeafRank}
\alias{LRsvm}
\alias{LRCart}
\alias{predict.TR_LRCart}
\alias{predict.LRsvm}
\alias{TR_LRCart}
\alias{varImportance.TR_LRCart}


\title{LeafRank functions for internal split of TreeRank master trees}
\description{Documentation on the two native LeafRank procedures for \code{\link{TreeRank}} and
on how to implement new LeafRank procedures}

\usage{
LRCart(formula,data,bestresponse, weights = NULL, maxdepth = 10,
minsplit = 50, nfcv = 0, prcSplitVar = 1,
prcSplitData = 1, mincrit = 0)

LRsvm(formula, data, bestresponse, weights = NULL, prcSplitVar = 1,
prcSplitData = 1,\dots)

\method{varImportance}{TR_LRCart}(obj,norm =TRUE);

\method{predict}{TR_LRCart}(object, newdata,\dots)

\method{predict}{LRsvm}(object, newdata,\dots)
}

\arguments{
\item{formula}{an object of class \code{formula}: a symbolic description of the model to be fit.}
\item{data}{a data frame containing the training data.}
\item{bestresponse}{The response/class value denoting the positive examples.}
\item{weights}{an optional vector of weights; should be \code{NULL} or a numeric vector.}
\item{nfcv}{the number of sets for the n-fold cross validation for the pruning procedure. If it is set to 0 or 1, no pruning is done.}\item{prcSplitVar}{the percentage of variables to be used to compute the model, drawn randomly. Must be a real value between 0 and 1.}
\item{prcSplitData}{the percentage of examples to be used to compute the model, drawn randomly. Must be a real value between 0 and 1.}
\item{maxdepth}{the maximal depth for the output tree.}
\item{minsplit}{number of observations/examples in a node at which no further split is desired in the growing procedure.}
\item{mincrit}{the minimum value of the criterion to split a node.}
\item{obj}{a \code{TR_LRCart} object.}
\item{object}{a \code{TR_LRCart} or \code{LRsvm} object.}
\item{newdata}{an dataframe containing data for prediction.}
\item{norm}{a boolean value; if it is set to true, the returned vector is  normalized.}
\item{\dots}{additional arguments.}
}

\details{

LeafRank is an abstract procedure in charge of the internal split of the \code{\link{TreeRank}} procedure. Basically, the procedure takes as input a \code{formula}, a dataframe \code{data} and
a value \code{bestresponse} denoting the class/response value of the positive examples.
The output is an object with at least a  \code{\link{predict}} method with real output values. The convention for internal split of \code{\link{TreeRank}} is that examples with negative predicted values have a best score than the ones with positive value, and so classified on the left branch of the master \code{\link{TreeRank}} tree.

Two version are provided with the \code{\link{TreeRank}} package: \code{LRCart}, based on a cost-sensitive version of the CART algorithm, using the \code{\link{rpart}} package, and  \code{LRsvm},
based on cost-sensitive SVMs, using the \code{kernlab} package and \code{\link{ksvm}} procedure.
For the \code{LRsvm} procedure, additional parameters for the \code{\link{ksvm}} procedure can be passed through the \code{\dots} argument.
A \code{TR_LRCart} object can be plotted with the \code{\link{TRplot}} procedure.

The \code{varImportance} method is used to compute the variable importance of the \code{LeafRank} model when a such computation is available.

In order to implement new \code{LeafRank} procedures, two procedure have to be implemented: the core one taking as input at least the \code{formula}, \code{data} and \code{bestresponse} arguments and returning a model; a \code{\link{predict}} method taking as input a model and a dataframe \code{newdata} and returning a real value. 
The custom \code{LeafRank} procedure can be used by passing the name of the procedure to the \code{LeafRank} argument of the \code{\link{TreeRank}} procedure.
See the source code of \code{LRsvm} and \code{LRCart} for more details.

Optionaly, a \code{varImportance} method can be implemented, which returns a named list where each name is a variable name and the value is the score of the variable.

See also \code{\link{TRGui}} and \code{\link{TRplot}} in order to integrate a custom \code{LeafRank} procedure to the graphical interface.
}

\seealso{\code{\link{TreeRank}}, \code{\link{TRplot}}, \code{\link{TRGui}}} 


\value{
  \code{LRCart} returns an object of class \code{TR_LRCart}.
  
  \code{LRsvm} returns an object of class \code{LRsvm}.
  
\code{predict} returns a real value vector.

\code{varImportance.TR_LRCart} returns a named list where each name is a variable name and the value is the score of the variable.
}


\examples{

data("TRdata");

## Building a LeafRank Cart model
cmod <- LRCart(class~., Gauss2D.learn, bestresponse = 1)

## variable importance of the computed model
varImportance(cmod)

## Building a LeafRank SVM model with a polynomial kernel of degree 2 
csvm <- LRsvm(class~., Gauss2D.learn, bestresponse = 1,
kernel="polydot", kpar = list(degree = 2))


}

\author{N. Baskiotis}

\references{
S. Clemencon and N. Vayatis. \emph{Tree-based ranking methods.}
 IEEE Transactions on Information Theory, 55(9):4316-4336, 200
}

