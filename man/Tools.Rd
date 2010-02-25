\name{Tools}
\title{TreeRank tools}
\alias{getROC}
\alias{getClassifier}
\alias{subTreeRank}
\alias{getCurves}
\alias{getCurves.TR_forest}
\alias{getCurves.TR_TreeRank}
\alias{getPREC}
\alias{getPREC.TR_forest}
\alias{getPREC_TR_TreeRank}
\alias{getROC.TR_forest}
\alias{getROC.TR_TreeRank}
\alias{auc}
\alias{plotROC}

\description{Tools for the manipulation and the analysis of TreeRank models}

\usage{

auc(roc)
getROC(obj,data = NULL);
getPREC(obj,data=NULL);
getCurves(obj,data=NULL);
getClassifier(tree,id)
subTreeRank(tree,listnode)
plotROC(rocs,colorlist=NULL,points=NULL)
}

\arguments{
\item{data}{a dataframe containing the data.}
\item{obj}{a \code{TR_TreeRank} or \code{TR_forest} object.}
\item{tree}{a \code{TR_TreeRank} or \code{TR_forest} object.}
\item{id}{ a node id of the \code{TR_TreeRank} object.}
\item{listnode}{a list of node id of the \code{TR_TreeRank} object.}
\item{roc}{a ROC as outputed by getROC (2D-matrix with the coordinates of the points).}
\item{rocs}{ a ROC list, as outputed by getROC.}
\item{colorlist}{a color list.}
\item{points}{a list of points coordinates to emphasize}.
}

\details{
\code{getClassifier} returns the classifier of the \code{TreeRank} node id.

\code{subTreeRank} computes the subtree with the same root as the \code{TreeRank} \code{tree} and with
leaves in \code{listnode}.

\code{getROC} computes the ROC corresponding to the predicted scores for \code{data} by the object \code{obj}.
\code{getPREC} computes the precision/recall curve.
\code{getCurves} computes ROC and precision/recall curves.

\code{plotROC} plots a list of ROCs, with the specified colors and
emphazing the specified points. The argument has to be a list even if it
is a list of one element.

\code{auc} computes the AUC of the \code{roc}.
}


\value{

\code{getROC} and \code{getPREC} returns a 2D-matrix with the points coordinates of the
ROC.

\code{auc} returns a real value.

\code{getClassifier} returns a \code{\link{LeafRank}} object.

\code{subTreeRank} returns a \code{TreeRank} object.
}

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


}


\author{N. Baskiotis}

\references{
S. Clemencon and N. Vayatis. \emph{Tree-based ranking methods.}
 IEEE Transactions on Information Theory, 55(9):4316-4336, 200
}

