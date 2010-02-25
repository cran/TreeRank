\name{TreeRank}
\alias{TreeRank}
\title{TreeRank procedures}
\alias{TreeRankForest}
\alias{TR_TreeRank}
\alias{growing_ctrl}
\alias{predict.TR_TreeRank}
\alias{predict.TR_forest}
\alias{print.TR_TreeRank}
\alias{TR_forest}
\description{TreeRank and TreeRankForest procedures for bipartite ranking.}
\usage{
TreeRank(formula,data,bestresponse, weights=NULL,growing=growing_ctrl(),LeafRank=LRCart,nfcv=0,varsplit=1)
TreeRankForest(formula,data,bestresponse,ntree=5,replace=TRUE,sampsize=0.5,varsplit=1,\dots)

predict.TR_TreeRank(object, newdata = NULL, type = "score",\dots)
predict.TR_forest(object,newdata=NULL,\dots)
growing_ctrl(minsplit=50,maxdepth=10,mincrit=0)
}

\arguments{
\item{formula}{an object of class \code{formula}: a symbolic description of the model to be fit.}
\item{data}{a dataframe containing the data.}
\item{newdata}{a dataframe containing the data.}
\item{bestresponse}{The response/class value denoting the positive examples.}
\item{weights}{Unused.}
\item{growing}{parameters for the tree growing procedure as returned by \code{growing_ctrl}}.
\item{LeafRank}{a \code{\link{LeafRank}} procedure.}
\item{nfcv}{the number of set for the n-fold cross validation pruning procedure. If nfcv<2, no pruning is done.}
\item{varsplit}{the percentage of variables to be used for the computation of each internal node in the TreeRank tree, drawn randomly. Must be a real value between 0 and 1.}
\item{ntree}{the number of trees to be computed for a forest.}
\item{sampsize}{For bagging and forest purposes, the percentage of examples to be used for the computation of a tree. Must be a real value between 0 and 1.}
\item{replace}{If sampsize < 1, it denotes when the examples are drawn with or without replacement.}
\item{maxdepth}{The maximal depth of any node of the final tree.}
\item{minsplit}{The minimum number of observations that must exist in a node in order for a split to be attempted.}
\item{mincrit}{The minimum value of the split criterion in order for a split to be attempted.}
\item{object}{a TR_TreeRank object.}
\item{\dots}{additinal arguments}
\item{type}{if this argument is fixed to the string "node", the
  prediction function returns the id of the leaf where the example falls.}
}
\details{

TreeRank is a learning algorithm tailored for ROC optimization in the context of bipartite ranking. 

The main procedure \code{TreeRank} uses the argument \code{\link{LeafRank}} to compute the internal splits of the master ranking Tree.

\code{TreeRankForest} is the forest/bagging version of the TreeRank procedure, computing and aggregationg \code{ntree} TreeRank trees. \code{sampsize} denotes the percentage of the examples to be drawn randomly (with replacement if \code{replace} is set to true) to compute each tree. \code{varsplit} denotes the percentage of variables to be used at each internal node. The other parameters are the same as \code{TreeRank}.

\code{growing_ctrl} returns a list of parameters for the growing procedure.

\code{predict} returns the predicted score for each example in the submitted dataframe \code{newdata}.
The \code{predict.TR_TreeRank} method has an optional argument \code{type} for returning the leaf id where the example falls when it is fixed to "node".


}

\value{
\code{TreeRank} returns an object of class \code{TR_TreeRank}.

\code{TreeRankForest} returns an object of class \code{TR_forest}.

\code{growing_ctrl} returns a list.
}

\examples{

data("TRdata");

## TreeRank learning with default parameters and Cart-like LeafRank procedure.
tree <- TreeRank(class~., Gauss2D.learn, bestresponse = 1)


## Predict scores on training and test sets
scoretrain <- predict(tree,Gauss2D.learn);
scoretest <- predict(tree,Gauss2D.test);



## TreeRank learning with svm-like LeafRank procedure, using  70%
##of the variables drawn randomly and with a maximal depth
##of 5 for the master tree.  
treesvm <- TreeRank(class~., Gauss20DFar.learn, bestresponse = 1,
growing=growing_ctrl(maxdepth=5), nfcv = 5, LeafRank = LRsvm);

}

\author{N. Baskiotis}

\references{
S. Clemencon and N. Vayatis. \emph{Tree-based ranking methods.}
 IEEE Transactions on Information Theory, 55(9):4316-4336, 200
}





