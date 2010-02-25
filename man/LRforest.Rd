\name{LRforest}
\alias{LRforest}
\alias{TR_LRforest}
\alias{predict.TR_LRforest}
\title{Random Forest implementation of  LeafRank procedure}
\description{Random Forest implementation of  \code{\link{LeafRank}} procedure}

\usage{
LRforest(formula,data,bestresponse, wpos=0.5,mtry=(ncol(data)-1),prcsize=1,\dots)

predict.TR_LRforest(object,newdata,\dots)
}

\arguments{
\item{formula}{an object of class \code{formula}: a symbolic description of the model to be fit.}
\item{data}{a data frame in which to interpret the variables named in the formula.}
\item{bestresponse}{The response/class/label value denoting the positive examples.}
\item{wpos}{weight of positive (bestresponse) examples in the data.}
\item{mtry}{Number of variables randomly sampled as candidates at each split.}
\item{prcsize}{Percentage of data to draw to compute a tree.} 
\item{object}{an object of class TR_LRCart.}
\item{newdata}{a dataframe.}
\item{\dots}{optional parameters for randomForest function.}
}

\details{

LRforest is the Random Forest-based version of the \code{\link{LeafRank}} procedure. This version uses the \code{\link{randomForest}} package to compute a random forest classifier  from the data. The predict function can be used to predict label from new data. \code{\link{varImportance}} can be used to compute the variable importance of the model. The random forest parameters are identical to the \code{\link{randomForest}} function.

}

\value{
\code{LRforest} returns an object of class \code{TR_LRforest}.
}
\seealso{LeafRank,varImportance,predict}

\examples{
data("TRdata");

## Building a LeafRank Random Forest model
cmod <- LRforest(class~., Gauss2D.learn, bestresponse = 1)

}


\author{N. Baskiotis}

\references{
S. Clemencon and N. Vayatis. \emph{Tree-based ranking methods.}
 IEEE Transactions on Information Theory, 55(9):4316-4336, 200
}


