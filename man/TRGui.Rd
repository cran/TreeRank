\name{TRGui}
\alias{TRGui}
\alias{TRplot}

\title{Graphical interfaces for the TreeRank package.}

\description{Graphical interfaces for the TreeRank package.}

\usage{
TRGui()
TRplot(tree, top = NULL)

}

\arguments{
\item{tree}{a \code{\link{TR_TreeRank}} or  \code{\link{TR_LRCart}} object.}
\item{top}{a tkframe to draw the object.}
}

\details{

\code{TRGui} starts the graphical user interface for launching \code{\link{TreeRank}} or \code{\link{TreeRank}}.


\code{TRplot} starts the graphical user interface to interact with the resulting model of \code{\link{TreeRank}} and \code{\link{TreeRankForest}} computation.

All the functionalities are detailled in the \emph{Getting Started} documentation.

Custom \code{\link{LeafRank}} procedures can be automatically integrated in \code{TRGui} by implementing 
an optional object in the global environment. To be recognize, the
object name has to start  with \emph{LR}
and end by \emph{def} (for example \emph{LRCart.def} is the object corresponding to the \code{\link{LRCart}} LeafRank procedure). 

This object has to be a list with the following items:
\itemize{
  \item{\code{name}: the name of the associated custom \code{\link{LeafRank}} procedure.}
\item{\code{fun}: the name of the implemented \code{\link{LeafRank}}
procedure in R.}
\item{\code{opt}: a list defining the options for the custom
\code{\link{LeafRank}}. Each item of the list has to contain the fields:
\itemize{
  \item{\code{name}: denoting the name of the option}
  \item{\code{default}:denoting the default value of the
    option}
    \item{\code{type}: which can be \emph{listbox},
      \emph{entry} or \emph{check} denoting the tcltk widget type of the
     option}
   \item{\code{optName}: denoting the variable name of the
     options.}
 }
 }
}


 When the type is set to \emph{listbox}, two another fields are required, the field \emph{choiceName} corresponding to a list with the name of each listbox choice, and the field \code{choice}
corresponding to a list with the variable value of each listbox choice.

\code{opt2cmd} a function which takes as input a list composed by the  set of \code{optName} specified in the \code{opt} list and the corresponding  value entered by the user. The result of this function has to be the string  representing the arguments passed to the \code{\link{LeafRank}} procedure. 

Refer to the source code for two examples. 


}


\examples{
## Loading example datasets

data("TRdata");

## Launching TreeRank GUI
TRGui()

## Plot a TreeRank tree

TRplot(TreeRank(class~.,Gauss2D.learn,bestresponse = 1))
}



\author{N. Baskiotis}

\references{
S. Clemencon and N. Vayatis. \emph{Tree-based ranking methods.}
 IEEE Transactions on Information Theory, 55(9):4316-4336, 200
}

