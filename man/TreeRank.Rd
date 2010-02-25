\name{TreeRank}
\title{TreeRank}
\alias{TreeRank}
\alias{TreeRankForest}
\alias{TR_TreeRank}
\alias{growing_ctrl}
\alias{predict.TR_TreeRank}
\alias{predict.TR_forest}
\alias{varImportance.TR_TreeRank}
\alias{varImportance.TR_forest}
\alias{auc}
\alias{getROC}
\alias{getROC.TR_TreeRank}
\alias{getROC.TR_forest}
\alias{getClassifier}
\alias{subTreeRank}
\alias{varImportance}
\alias{print.TR_LRCart}
\alias{print.TR_TreeRank}
\alias{plotROC}
\description{TreeRank procedures for bipartite ranking.}

\usage{
TreeRank(formula,data,bestresponse, weights = NULL,
growing = growing_ctrl(),LeafRank=LRCart, nfcv = 0,
cbFct = NULL, prcSplitVar = 1, prcSplitData = 1)

TreeRankForest(nbforest = 5,  \dots)

growing_ctrl(minsplit = 50, maxdepth = 10, mincrit = 0)

\method{varImportance}{TR_TreeRank}(obj,norm =TRUE);
\method{varImportance}{TR_forest}(obj,norm =TRUE);

\method{predict}{TR_forest}(object, newdata = NULL,\dots)
\method{predict}{TR_TreeRank}(object, newdata = NULL,type = "score",\dots)
auc(roc)

\method{getROC}{TR_TreeRank}(obj,data = NULL);
\method{getROC}{TR_forest}(obj,data = NULL);

getClassifier(tree,id)

subTreeRank(tree,listnode)

plotROC(rocs,colorlist=NULL,points=NULL)
}

\arguments{
\item{formula}{an object of class \code{formula}: a symbolic description of the model to be fit.}
\item{data}{a dataframe containing the data.}
\item{bestresponse}{The response/class value denoting the positive examples.}
\item{weights}{an optional vector of weights; should be \code{NULL} or a numeric vector.}
\item{growing}{parameters for the tree growing procedure as returned by \code{growing_ctrl}.}
\item{LeafRank}{a \code{\link{LeafRank}} procedure.}
\item{nfcv}{the number of set for the n-fold cross validation for the pruning procedure. If it is set to 0 or 1, no pruning is done.}
\item{cbFct}{an optional CallBack function, called after each TreeRank node computation. The input arguments are: the current node id,
the depth of the current node, the number of positive and negative examples matching the current node.}
\item{prcSplitVar}{the percentage of variables to be used to compute the model, drawn randomly. Must be a real value between 0 and 1.}
\item{prcSplitData}{the percentage of examples to be used to compute the model, drawn randomly. Must be a real value between 0 and 1.} 
\item{nbforest}{the number of trees to be computed for a forest.}
\item{maxdepth}{the maximal depth of the computed tree.}
\item{minsplit}{number of observations/examples in a node at which no further split is desired in the growing procedure.}
\item{mincrit}{the minimum value of the criterion to split a node.}
\item{obj}{a \code{TR_TreeRank} or \code{TR_forest} object.}
\item{object}{a \code{TR_TreeRank} or \code{TR_forest} object.}
\item{newdata}{an optional data frame containing data for prediction.}
\item{type}{set to \code{score} or \code{node}.}
\item{norm}{a boolean value; if it is set to true, the returned vector is  normalized.}
\item{tree}{a \code{TR_TreeRank} or \code{TR_forest} object.}
\item{id}{ a node id of the \code{TR_TreeRank} object.}
\item{listnode}{a list of node id of the \code{TR_TreeRank} object.}
\item{roc}{a ROC as outputed by getROC (2D-matrix with the coordinates of the points).}
\item{rocs}{ a ROC list, as outputed by getROC.}
\item{colorlist}{a color list.}
\item{points}{a list of points coordinates to emphasize}.
\item{\dots}{additional arguments.}
}

\details{

TreeRank is a learning algorithm tailored for ROC optimization in the context of bipartite ranking.

The main procedure \code{TreeRank} uses the argument \code{\link{LeafRank}} to compute the internal splits
of the master ranking tree. 

\code{TreeRankForest} is the forest version of the TreeRank procedure, computing \code{nbforest} TreeRank trees. The other parameters are the same as \code{TreeRank}.

\code{growing_ctrl} returns a list of parameters for the growing procedure.

\code{varImportance} is used to compute the variable importance of the \code{\link{LeafRank}} model when a such computation is available.

\code{predict} returns the predicted score for each example in the submitted dataframe \code{newdata}.
The \code{predict.TR_TreeRank} method has an optional argument \code{type} for returning the leaf id where the example falls when it is fixed to "node".

\code{getClassifier} returns the classifier of the \code{TreeRank} node
id.

\code{subTreeRank} computes the subtree wioth the same root as the \code{TreeRank} \code{tree} and with
leaves in \code{listnode}.

\code{plotROC} plots a list of ROCs, with the specified colors and
emphazing the specified points. The argument has to be a list even if it
is a list of one element.

\code{getROC} computes the ROC corresponding to the predicted scores for \code{data} by the object \code{obj}.
\code{auc} computes the AUC of the \code{roc}.
}

\value{

  \code{TreeRank}  returns an object of class \code{TR_TreeRank}.
  
  \code{TreeRankForest} returns an object of class \code{TR_forest}.
  
  \code{growing_ctrl} returns a list.
  
\code{predict} returns a real value vector corresponding to the score of
each example.

\code{varImportance} returns a named list where each name is a variable
name and the value is the score of the variable.

\code{getROC} returns a 2D-matrix with the points coordinates of the
ROC.

\code{auc} returns a real value.

\code{getClassifier} returns a \code{\link{LeafRank}} object.

\code{subTreeRank} returns a \code{TreeRank} object.
}


\seealso{\code{\link{LeafRank}},\code{\link{TRplot}},\code{\link{TRGui}}}

\examples{


data("TRdata");

## TreeRank learning with default parameters and Cart-like LeafRank procedure.
tree <- TreeRank(class~., Gauss2D.learn, bestresponse = 1)

## Compute the variable importance of the computed model
varImportance(tree)

## Compute the ROC curve on training and test sets

roctrain <-getROC(tree);
roctest <- getROC(tree,Gauss2D.test);

## Predict scores on training and test sets
scoretrain <- predict(tree,Gauss2D.learn);
scoretest <- predict(tree,Gauss2D.test);

## Compute the AUC for the training and test sets
auctrain <- auc(getROC(tree));
auctest <- auc(getROC(tree,Gauss2D.test));

## Plot the ROC curve for the training and test sets
plotROC(list(getROC(tree)));
plotROC(list(getROC(tree),getROC(tree,Gauss2D.test)));

## Return the LeafRank classifier of the root node
c1 <- getClassifier(tree,1)


## TreeRank learning with svm-like LeafRank procedure, using  70%
##of the variables drawn randomly and with a maximal depth
##of 5 for the master tree.  
treesvm <- TreeRank(class~., Gauss20DFar.learn, bestresponse = 1,
growing=growing_ctrl(maxdepth=5), nfcv = 5, LeafRank = LRsvm,
prcSplitVar = 0.70);

## Plot the tree
TRplot(tree);

}

\author{N. Baskiotis}

\references{
S. Clemencon and N. Vayatis. \emph{Tree-based ranking methods.}
 IEEE Transactions on Information Theory, 55(9):4316-4336, 200
}





