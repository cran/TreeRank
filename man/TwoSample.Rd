\name{TwoSample}

\title{TreeRank procedure for the two sample problem.}

\alias{TwoSample}
\alias{TwoSampleGui}

\description{TreeRank procedure for the two sample problem.}

\usage{
TwoSample(x,y, split = 60, TRalgo = TreeRank, alpha= 5);
TwoSampleGui();
}

\arguments{
\item{x}{ a data frame representing the first sample.}
\item{y}{ a second data frame with the same variables as x, representing the second sample.}
\item{split}{the split percentage for the training and test set.}
\item{TRalgo}{the TreeRank procedure to be used.}
\item{alpha}{the confidence level.}
}

\details{

\code{TwoSample} is an adaptation of the TreeRank procedure for the two sample problem.
\code{TwoSampleGui} launches the interface for the TwoSample procedure.
}

\value{

\code{TwoSample} returns a list composed by a \code{wilcox_test} object (from package coin), the computed TreeRank object, the train and test drawn sets.
}

\examples{


## Launching TwoSample GUI
TwoSampleGui()

## A two sample test.
data("TRdata");
TwoSample(Gauss2D.learn[1:500,],Gauss2D.learn[501:1000,]);
}


\author{N. Baskiotis}

\references{
S. Clemencon and N. Vayatis. \emph{Tree-based ranking methods.}
 IEEE Transactions on Information Theory, 55(9):4316-4336, 200

S. Clemencon and M. Depecker and N. Vayatis
\emph{AUC optimization and the two-sample problem},
Advances in Neural Information Processing Systems, 2009
}

