#################################################################################
#
# LeafRankCart functions : implementation of a LeafRank classifier,
#                          a modified Cart classifier.
#
#################################################################################






#################################################################################
#
# predict.TR_LRCart(object, newdata, type) :
#                           object : LRCart object
#                           newdata :
#                           type : "node" to predict terminal nodes
#                                  otherwise -1 for best label, +1 for worst
#
#################################################################################


predict.TR_LRCart <- function(object, newdata,...)
{
  if (!inherits(object,"TR_LRCart"))
    stop("object not of class TR_LRCart");
  #tmp <- predict.party(object,newdata,type = "node")
  #if(!(is.null(type)))
   # if (type == "node")
   #   return(tmp);
  #pred <- array(sapply(tmp, function(x) {if (x %in% object$info$Lnode) -1 else 1}))
  #pred
  return(predict(object$fustree,newdata));
}


#################################################################################
#
# LRCartFusion(tree, pcInit,ncInit)
#    Reorder leafs for concavification purpose and merge leafs in two subsets
#
#    return : 2 lists, LeftNode set and RightNode set
#
#################################################################################


LRCartFusion <- function(tree,pcInit,ncInit){


  # Retrieve  informations on true postive and true negative
 
  listLeafs <- (1:tree$nbNode)[tree$isleaf];
  if (length(listLeafs) <2)
    {return(list(Lnode = tree$root, Rnode = 0))}
  betaList <- tree$pcount[listLeafs]
  alphaList <- tree$ncount[listLeafs]
   
  betaList <- betaList/pcInit;
  alphaList <- alphaList/ncInit;

  #Ratio list of  true positive and  false negative for each leaf
  crVec <- betaList/alphaList;

  #Reorder leafs
  listIndex <- order(-crVec)

  #Compute the cumulate vector of ratio TP and FN for each leaf
  alphaListOrd <- cumsum(alphaList[listIndex])
  betaListOrd <- cumsum(betaList[listIndex])
  
  #Compute entropy, find the best and divide leafs between right and left subset
  entropy <- betaListOrd-alphaListOrd
  Lnode <- listLeafs[listIndex[1:which.max(entropy)]]
  Rnode <- listLeafs[listIndex[as.integer(which.max(entropy)+1):length(listIndex)]]
  list(Lnode = Lnode, Rnode = Rnode)
}



#################################################################################
#
# LRCart(formula, data, bestresponse, weights, growing, pruning) : Main Cart function
#
#        formula, data, weights :
#        bestresponse: value of the label considered as best
#       growing, pruning :
#
#       Return : a LeafRank Cart classifier, type TR_LRCart
#################################################################################
LRCart <- function(formula, data,bestresponse, weights = NULL,maxdepth=10,minsplit=50,nfcv=0,prcSplitVar = 1, prcSplitData = 1,mincrit=0){
  
  evaluation <- function(y, wt, parms){
    idx <- y== bestresponse;
    pc <- sum(wt[idx])*(1-parms$pInit);
    nc <- sum(wt[!idx])*parms$pInit;
    label <- bestresponse;
    miss <- nc;
    if (nc>pc){
      label <-parms$neglab;
      miss <- pc;
    }
    return(list(label=label,deviance=miss))
  }

  split <- function(y,wt,x,parms,continuous,...){
    n <- length(y)
    
    pInit<-parms$pInit;
    if (!continuous){
      stop("Only continuous variables are handled");
    }
    pvec <- y==bestresponse;
    nvec <- !pvec;
    pcount <- sum(pvec*wt)
    ncount <- sum(nvec*wt)
    leftpos <- cumsum(pvec*wt)[-n];
    leftneg <- cumsum(nvec*wt)[-n];
    WERMLess <- 2-(2*pInit*leftneg/(pcount+ncount)+2*(1-pInit)*(pcount-leftpos)/(pcount+ncount));
    WERMGreat <- 2-(2*pInit*(ncount-leftneg)/(pcount+ncount)+2*(1-pInit)*(leftpos)/(pcount+ncount));
    ret <- list(goodness= WERMLess,direction=rep(-1,(n-1)))
    if (max(WERMLess)<max(WERMGreat)){
      ret <-list(goodness=WERMGreat,direction=rep(-1,(n-1)))
    }
    return(ret);
  }

  init <- function(y,offset,parms,wt){
    pcInit=sum((y==bestresponse)*wt);
    ncInit = sum((y!=bestresponse)*wt);
    ntot = pcInit+ncInit;
    pInit = pcInit/ntot;
    neglab<-y[which(y!=bestresponse)][[1]]
    list(y=y, parms=list(pInit=pInit,ntot=ntot,poslab=bestresponse,
                neglab=neglab), numresp=1, numy=1,
	      summary= function(yval, dev, wt, ylevel, digits ) {
		  paste("  mean=", format(signif(yval, digits)),
			", MSE=" , format(signif(dev/wt, digits)),
			sep='')})
  }

  #initiallisation
  if (missing(data))
    data <- environment(formula)
 
  mf <- call <-  match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"),
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- FALSE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  y <-model.response(mf)
                                        #  if (is.null(weights)) weights <- rep.int(1, length(y))
  rn <- names(mf)[1] ### model.response(mf)
  x <- mf[,colnames(mf)!="(weights)"]
  inputs <- which(!(colnames(x) %in% rn));
  
  idxVar <- colnames(x);
  idxData <- 1:nrow(data);
  if (prcSplitVar<1){
    idxVar <- c(which(colnames(x)%in%rn),sample(inputs,ceiling(length(inputs)*prcSplitVar)));
  }
   if (prcSplitData<1){
     idxData <- sample(nrow(data),ceiling(nrow(data)*prcSplitData));
  }
  pcInit <- y[idxData] == bestresponse;
  ncInit <- sum(!pcInit);
  pcInit <- sum(pcInit);
  pInit <- pcInit/(pcInit+ncInit);
  
  if ((pInit == 1) || (pInit == 0))
    {
      return(NULL);
   }
  alist <- list(eval=evaluation,split=split,init=init);
  tmpCRdata <- mf[idxData,idxVar];
  rtree <- rpart(formula,tmpCRdata,method=alist,control=rpart.control(cp=0,maxdepth=maxdepth,minsplit=minsplit,maxsurrogate=0,maxcompete=0,minbucket=1),model=TRUE,y=TRUE,xval=0);

  tree <- rpart2TR(rtree,bestresponse);
  class(tree) <- "TR_LRCart";
  tree$terms <- rtree$terms;
  tree$rpart <- rtree;
  tree$formula <- formula;
  tree$call <- call;
  tree$bestresponse <- bestresponse;
  if (nfcv>0){
    oldtree <- tree;
    rtree <- .pruneRpart(rtree,bestresponse,pInit,nfcv)
    tree <- rpart2TR(rtree,bestresponse);
    class(tree) <- "TR_LRCart"
    tree$unpruned  <- oldtree;
    tree$rpart <-rtree;
    tree$formula <- formula;
    tree$call <- call;
    tree$terms <- rtree$terms;
    tree$bestresponse <- bestresponse;
  }
  splitNode <- LRCartFusion(tree,pcInit,ncInit);
  tree$Lnode <-splitNode$Lnode;
  tree$Rnode <- splitNode$Rnode;
  fustree <- tree$rpart;
  for (i in splitNode$Lnode){
    fustree$frame[i,"yval"]<- -1;
  }
  for (i in splitNode$Rnode){
    fustree$frame[i,"yval"]<-1;
  }
  tree$fustree <- fustree;
  return(tree);
}

.pruneRpart <- function(tree,bestresponse,pInit,nfcv){
  if (nrow(tree$frame)<4)
    return(tree)
  
    cvinfo <- xpred.rpart(tree,xval = nfcv)
    cverr <- apply(cvinfo,2,function(x) {
      idx <- which(tree$y != x);
      err <- (pInit)*sum(tree$y[idx]!=bestresponse);
      err <- err+ (1-pInit)*sum(tree$y[idx]==bestresponse);
      return(err)
    })
    copt <- which.min(cverr[-1])+1
  
    rtree <- prune.rpart(tree,cp = tree$cptable[[copt]]);
   return(rtree)
  }
  



print.TR_LRCart<- function(x,...){
object <- x;
  if(!inherits(object,"TR_LRCart"))
    stop("object not of class TR_LRCart");

  id <- object$root;
  nodestack<-id;
  cat("TreeRank tree\n   id) var <|>= threshold #pos:#neg dauc \n\n");

  while(length(nodestack)>0){
    id <- nodestack[[1]];
    nodestack <-nodestack[-1];
    s <- "";
    sp <-"root";
    if (id != object$root){
       parent <- object$parentslist[id];
       if (object$kidslist[[parent]][[1]] == id)
          sp <- paste(object$split[[parent]]["name"],"<",format(object$split[[parent]]["breaks"],digits=3))
        else sp <- paste(object$split[[parent]]["name"],">=",format(object$split[[parent]]["breaks"],digits=3));
     }
    s <- paste(cat(rep(' ',2*object$depth[id])),id,"| ", sp,"  ",
               object$pcount[id],":",
               object$ncount[id]," ",format(object$ldauc[id],digits=3),sep="")
    
    if (!(object$isleaf[id])){
      nodestack <- c(object$kidslist[[id]][[1]],object$kidslist[[id]][[2]],nodestack);
    }else{s<- paste(s,"*");}
    cat(paste(s,"\n"));
  }
}
    
    

varImportance.TR_LRCart <- function(obj,norm=TRUE){
  nr <- attributes(obj$terms)$term.labels
  res <- array(0,length(nr))
  names(res)<- nr;
  listvar <- array(unlist(lapply(obj$inner,function(x){obj$split[[x]]["name"]})));
  for (i in 1:length(obj$inner))
     res[listvar[i]] <- res[listvar[i]]+obj$ldauc[obj$inner[i]]^2;
  if (norm) ret <- res/max(res)
  else ret <- res;
  ret
}



LRsvm <- function(formula,data,bestresponse,weights=NULL,prcSplitVar=1,prcSplitData=1,...){
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"),
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- FALSE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  y <-model.response(mf)
  if (is.null(weights)) weights <- rep.int(1, length(y))
  rn <- names(mf)[1] ### model.response(mf)

  x <- mf[,colnames(mf)!="(weights)"]
  inputs <- which(!(colnames(x) %in% rn));
  idxVar <- colnames(x);
  idxData <- 1:nrow(data);
  if (prcSplitVar<1){
    idxVar <- c(which(colnames(x)%in%rn),sample(inputs,ceiling(length(inputs)*prcSplitVar)));
  }
   if (prcSplitData<1){
     idxData <- sample(nrow(data),ceiling(nrow(data)*prcSplitData));
  }

  pc <- sum(weights[idxData][y[idxData] == bestresponse]);
  nc <- sum(weights[idxData][y[idxData]!=bestresponse]);
  classw <- c(nc/(pc+nc),pc/(pc+nc));
  nclass <- unique(data[,rn]);

  data[,rn] <- as.factor(data[,rn]);

  if (bestresponse !=nclass[[1]])
    {names(classw) <- nclass[c(2,1)];}
  else
    {names(classw) <- nclass;}
  
  if (classw[[1]]*classw[[2]] == 0)
    return(NULL);
  
  tmpSRdata <- data[idxData,idxVar];
#  m <- ksvm(formula,data,type="C-svc",class.weights=classw,prob.model=TRUE);
  m <- ksvm(formula,tmpSRdata,class.weights=classw,type="C-svc",...);
  
  ret <- list(svm=m,bestresponse= bestresponse);
  class(ret)<-"LRsvm";
  ret;
}



predict.LRsvm <- function(object,newdata=NULL,...){
  p <- predict(object$svm,newdata);
  ret <- rep(1,nrow(newdata));
  ret[p==object$bestresponse] <- -1;
  ret;
}
  


r2weka <- function(data,dest,lab){
  dest<-file(dest,open="w");
  writeLines("@relation 'cpu'",dest);
  for (i in names(data)){
    if (i !=lab)
    {writeLines(paste("@attribute",i,"real"),dest)}
    else{
      cat("@attribute",i,"{",file=dest);
      cat(paste(unlist(unique(data[,lab])),",",sep=""),file=dest);
      writeLines("}",dest);
    }
    
  }
  writeLines("@DATA",dest);
  for (i in 1:nrow(data)){
    cat(paste(unlist(data[i,]),",",sep=""),sep="",file=dest);
    cat("\n",file=dest);
  }
  close(dest)
}
    


