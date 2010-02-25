\name{LRCart}
\alias{LRCart}
\alias{TR_LRCart}
\alias{print.TR_LRCart}
\alias{predict.TR_LRCart}
\title{Cart implementation of  LeafRank procedure}
\description{Cart implementation of  \code{\link{LeafRank}} procedure}

\usage{
LRCart(formula,data,bestresponse, wpos=0.5,weights = NULL, maxdepth = 10,minsplit = 50, mincrit = 0, nfcv = 0)
predict.TR_LRCart(object,newdata,\dots)
}

\arguments{
\item{formula}{an object of class \code{formula}: a symbolic description of the model to be fit.}
\item{data}{a data frame in which to interpret the variables named in the formula.}
\item{bestresponse}{The response/class/label value denoting the positive examples.}
\item{wpos}{weight of positive (bestresponse) examples in the data.}
\item{weights}{Unused in the current version. An optional vector of weights.}
\item{maxdepth}{The maximal depth of any node of the final tree.}
\item{minsplit}{The minimum number of observations that must exist in a node in order for a split to be attempted.}
\item{mincrit}{The minimum value of the split criterion in order for a split to be attempted.}
\item{nfcv}{The number n of sets for the n-fold cross validation used for the pruning procedure. If n<=1, no pruning is done.} 
\item{object}{an object of class TR_LRCart.}
\item{newdata}{a dataframe.}
\item{\dots}{optional parameters for rpart function.} 
}

\details{

LRCart is the Cart-based version of the \code{\link{LeafRank}} procedure. This version uses the \code{\link{rpart}} package to compute a decision tree from the data. The computed decision tree is next transformed in a TR_LRCart object. The predict function can be used to predict label from new data. A TR_LRCart object can be drawn by the \code{\link{TRplot}} function. \code{\link{varImportance}} can be used to compute the variable importance of the model.

}

\value{
\code{LRCart} returns an object of class \code{TR_LRCart}.
}
\seealso{LeafRank,varImportance,predict,TRplot}

\examples{
data("TRdata");

## Building a LeafRank Cart model
cmod <- LRCart(class~., Gauss2D.learn, bestresponse = 1)

}



