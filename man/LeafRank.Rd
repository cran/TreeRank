\name{LeafRank}
\alias{LeafRank}

\title{LeafRank abstract function for internal splits of TreeRank master tree}
\description{Structure and implementation of a generic LeafRank procedure}


\details{
LeafRank is the abstract procedure used for the computation of internal splits of TreeRank master trees. Basically, a LeafRank procedure computes a binary classifier from a dataset and a formula in order to split the input space in two regions, one with the "best" instances and the other one with the remaining instances. A custom LeafRank procedure has to take as input at least the following arguments : formula (of type R formula), data (a R dataframe), bestresponse (the value of the wanted positive label), wpos (the weight of positive instances  in the dataframe).
The output of the LeafRank procedure is a predictable model with the following convention : the instances with a predicted value <0 are considered most positive than the instances with a predicted value >0 (and  conventionaly are classified on the left branch of the TreeRank master tree).  
If a custom LeafRank procedure is used, a custom \code{\link{predict}} function for the considered custom model has also to be provide. 

One can see the implementation of \code{\link{LRsvm}} for a simple example of LeafRank implementation.
}

\seealso{LRCart,LRforest,LRsvm}
\author{N. Baskiotis}

\references{
S. Clemencon and N. Vayatis. \emph{Tree-based ranking methods.}
 IEEE Transactions on Information Theory, 55(9):4316-4336, 200
}


