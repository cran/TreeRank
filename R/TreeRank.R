
#################################################################################
#
# TreeRank functions : implementation of the TreeRank algorithm
#
#################################################################################


#################################################################################
#
# predict.TR_TreeRank(object, newdata, type) : predict node or score for newdata
#             object : TreeRank object
# 	      newdata : data for the prediction
#      	      type: node : node id of the prediction
#                   score : normalized score (default)
#
#################################################################################



predict.TR_TreeRank<- function(object, newdata = NULL, type = "score",...)
{
  if (!inherits(object,"TR_TreeRank"))
    stop("object not of class TR_TreeRank");

 if (is.null(newdata))
    return(object$fitted[["(fitted)"]])
  
  retid <- rep(object$root,nrow(newdata));
  indextab <- list(rep(1:nrow(newdata)));
  nodestack <- list(object$root);
  
  while(length(nodestack)>0){
	id <- nodestack[[1]];

	nodestack <- nodestack[-1];
        if (object$isleaf[id]) {
		retid[indextab[[id]]]<- id;
 		next;
	}
        if (length(indextab[[id]]) >0){
	tmp <- predict(getClassifier(object,id),newdata[indextab[[id]],]);
        kids <- object$kidslist[[id]];
        indextab[kids[1]] <- list(indextab[[id]][tmp<0]);
        indextab[kids[2]] <- list(indextab[[id]][tmp>0]);          
        nodestack <- c(nodestack,kids[1],kids[2]);
}
  }
  if (type == "node")
      return(retid)
  return(object$score[retid]);
}




predict.TR_forest<- function(object, newdata = NULL,...){
  if (!inherits(object,"TR_forest"))
    stop("Object not of class TR_forest");
  
  if (is.null(newdata))
    newdata <- object$forest[[1]]$data;
  res <- array(0,dim=nrow(newdata))
  for (i in 1:object$ntree	)
  {
    tmp <- predict(object$forest[[i]],newdata);
    res <- res + tmp;
  }
  res <- res / object$ntree;
  res
}


#################################################################################
#
# getClassifier(tree, id)
#    tree: a TreeRank object
#    id : id of the node
#
# Return : the LeafRank Classifier associated to the node id
#
#################################################################################


getClassifier <- function(tree,id){
  if (!inherits(tree,"TR_TreeRank")){
    stop("object not of class TR_TreeRank");
  }
  tree$LRList[[id]];
}

#################################################################################
#
# TreeRankRec(id,formula,data,response, bestresponse, weights,growing,depth,LeafRank)
#   Main TreeRank recursive function
#    id : id of current underconstruction node
#    formula,data,response,weights:
#    bestresponse : best label value
#    growing : growing control
#    depth : current depth
#    LeafRank : classifier to use to do the split
#
# Return : a TreeRank tree
#
#################################################################################


TreeRankRec <- function(formula, data,  bestresponse, growing = growing_ctrl(), LeafRank = LRCart,varsplit=1){
   
  if (missing(data))
        data <- environment(formula)
  call <- mf  <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights"),names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- FALSE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Terms <- attr(mf,"terms");
  response <- model.response(mf);
  w <- model.extract(mf,"weights");
  if (length(w)==0L) w <- rep.int(1,length(response));
  rn <- names(mf)[1] 
 x <- mf[,colnames(mf) != "(weights)"]

  inputs <- which(!(colnames(x) %in% rn ))
  
   print(paste("Computing Master Tree - :"));
   pcInit <- sum(w[response == bestresponse]);
   ncInit <- sum(w[response != bestresponse]);
   #initialise root node
   id <- 1;
   nextnode <- 1;
   wtmp <- list(seq(1,length(response)));
   kidslist <- list();
   parentslist <- array();
   LRList <- list();
   isleaf <-array();
   ldauc <- array();
   score <- array();
   lalpha <- 0;
   lbeta <- 0;
   ncount <- ncInit;
   pcount <- pcInit;
   depth <- c(1);
   nodestack <- list(id);
   nodeorder <- array(); 
   curscore <- 1;
   #recursive while
   while(length(nodestack)>0){
     cat(".");
     id <- nodestack[[1]];
     isleaf[id] <- TRUE;
     nodestack <- nodestack[-1];
     tmpdata <- x[wtmp[[id]],];
     tmpweights <- w[wtmp[[id]]];
     tmpresponse <- response[wtmp[[id]]];
     pcount[id] <- sum(tmpweights[tmpresponse == bestresponse]);
     ncount[id] <- sum(tmpweights[tmpresponse != bestresponse]);
     ldauc[id] <- 0;
     nodeorder[id] <- curscore;
     curscore <- curscore+1;
     if ((depth[id] >=growing$maxdepth)|| ((pcount[id]+ncount[id])<growing$minsplit)|| (!checkImpure(tmpweights,tmpresponse))){

 	kidslist[id] <- NA;
 	LRList[id] <- NA;
        next;
    }
	if (varsplit <1)
	xid <- c(1,sample(2:ncol(tmpdata),ceiling(varsplit*(ncol(tmpdata)-1))))
	else xid <- 1:ncol(tmpdata);
	tmpdata <- tmpdata[,xid];
    #Build the LeafRank  classifier for the current node
    if (is.function(LeafRank)){
	lrTree <- LeafRank(formula = formula,data = tmpdata, bestresponse = bestresponse);
	}        
	else{
	LRlist <- LeafRank$LRlist;
	splitD <- LeafRank$split;
	idLearn <- sample(nrow(tmpdata),ceiling(nrow(tmpdata)*splitD));
	listTree <- list();
	listAUC <- list();
	maxauc <- 0;
	maxtree<- NULL;
	print("begin")
	for (alg in LRlist){
#browser()
		tmptree <- alg(formula=formula,data=tmpdata[idLearn,],bestresponse=bestresponse);
		lrResponse <- predict(tmptree,tmpdata[-idLearn,]); 
	        left <- lrResponse <=0;
	        right <- lrResponse >0;
		npc <- sum(tmpweights[(tmpresponse[-idLearn] == bestresponse) & left]);
		nnc <- sum(tmpweights[(tmpresponse[-idLearn] != bestresponse) & left]);
		auctmp = -((pcount[id]/pcInit)*(nnc/ncInit)-(ncount[id]/ncInit)*(npc/pcInit))/2
		if (auctmp>maxauc)
			{
			maxtree <- tmptree;
			maxauc <- auctmp;
			maxalg <- alg;
			}
		
	}
	lrTree <-  alg(formula=formula,data=tmpdata,bestresponse=bestresponse);
	print(class(lrTree))
     }

     if (is.null(lrTree))
       {
         kidslist[id]<-NA;
         LRList[id]<-NA;
         next;
       }
    lrResponse <- predict(lrTree,tmpdata);
   
  
    #Compute left and right node
    left <- lrResponse <=0;
    right <- lrResponse >0;
    npc <- sum(tmpweights[(tmpresponse == bestresponse) & left]);
    nnc <- sum(tmpweights[(tmpresponse != bestresponse) & left]);
    #If a kid is empty, the current node is a leaf
    if (((npc+nnc) == 0) || ((pcount[id]+ncount[id]-(npc+nnc)) == 0)){
       kidslist[id]<-NA;
       LRList[id] <- NA;
       next;
    }
    LRList[id] <- list(lrTree);
    ldauc[id] <- -((pcount[id]/pcInit)*(nnc/ncInit)-(ncount[id]/ncInit)*(npc/pcInit))/2
    nalpha <-  lalpha[id] + nnc/ncInit;
    nbeta <- lbeta[id] + npc/pcInit;
    #build kids
    kidslist[id] <- list(c(nextnode+1,nextnode+2));
    parentslist[nextnode+1] <- parentslist[nextnode+2]<-id;
    depth[nextnode+1]<- depth[nextnode+2] <- depth[id]+1;
    #left node
    lalpha[nextnode+1] <- lalpha[id]
    lbeta[nextnode+1] <- lbeta[id];
    wtmp[nextnode+1] <- list(wtmp[[id]][left]);
    #right node
    lalpha[nextnode+2] <- nalpha;
    lbeta[nextnode+2] <- nbeta;
    wtmp[nextnode+2] <- list(wtmp[[id]][right]);
    isleaf[id] <- FALSE;
    nodestack <- c(nextnode+1,nextnode+2,nodestack);
    nextnode <- nextnode+2;
   


   }
   cat("\n");

   nbNode <- nextnode;
   score <- array(0,nbNode);
   leaf <- which(isleaf);
   leafOrdered <- order(nodeorder);
   leafOrdered <- leafOrdered[leafOrdered %in% leaf];
   nbLeaf <- length(leafOrdered); 
   for (i in   0:(nbLeaf-1))
      score[leafOrdered[i+1]] <- ((nbLeaf -i)/nbLeaf)

   ret <- list();
   ret$leafOrdered <- leafOrdered;
   ret$nodeorder <- nodeorder;
   ret$nodes <- (1:nbNode)
   ret$inner <- ret$nodes[!(isleaf)]
   ret$root <- 1L;
   ret$parentslist <- parentslist;
   ret$kidslist <- kidslist;
   ret$lalpha <- lalpha;
   ret$lbeta <- lbeta;
   ret$ldauc <- ldauc;
   ret$isleaf <- isleaf;
   ret$bestresponse <- bestresponse;
   ret$pcount <- pcount;
   ret$ncount <- ncount;
   ret$depth <- depth;
   ret$LRList <- LRList;
   ret$nbNode <- nextnode;
   ret$score <- score;
   return(ret); 
}



#################################################################################
#
# TreeRankForest(formula,data,bestresponse,ntree,replace,sampsize,varsplit,...)
# TreeRank Forest version
#
#    formula,data 
#    bestresponse : value of the best label 
#    ntree : number of trees to be computed
#    sampsize : percent of the data to use for each tree
#    replace : drawing examples from data with replacement or not
#    varsplit : percent of variables to be used for each internal node of trees.
#
#################################################################################

TreeRankForest <- function(formula,data,bestresponse,ntree=5,replace=TRUE,sampsize=0.5,varsplit=1,...){
  forest <- list();
 if (missing(data))
        data <- environment(formula)
  call <- mf  <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights"),names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- FALSE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Terms <- attr(mf,"terms");
  for (i in 1:ntree){
    print(paste("Compute tree # ",i,"/",ntree,sep=""));
    tmpdata <- mf[sample(nrow(mf),ceiling(sampsize*nrow(mf)),replace=replace),]
    forest <- c(forest,list(TreeRank(formula,tmpdata,bestresponse,varsplit=varsplit,...)));
  }
  ret <- list(forest= forest,ntree = ntree);
  class(ret) <- "TR_forest";
  ret
}


#################################################################################
#
# TreeRank(formula,data,bestresponse,weights,growing,LeafRank,nfcv)
# Main TreeRank function
#
#    bestresponse : value of the best label 
#    weights : weights on the data
#    growing : growing control of the TreeMaster
#    LeafRank : classifier to use at each node
#               classifier has to take at least args : formula, data, bestresponse, weights
#                   and returns -1/+1 : best/worst ones
#    pruning : pruning control
#    nfcv :  number of cases for  n-fold cross pruning procedure , 0 or 1 = no pruning.
#    varsplit : percent of variables to be used for each node.
#
#################################################################################



TreeRank <- function(formula, data,bestresponse, weights=NULL,growing = growing_ctrl(), LeafRank = LRCart,nfcv=0,varsplit=1)
{


  #Prepare the call to the recursive TreeRank function
  if (missing(data))
        data <- environment(formula)
  call <- mf  <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights"),names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- FALSE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Terms <- attr(mf,"terms");
  y <- model.response(mf);
  w <- model.extract(mf,"weights");
  if (length(w)==0L) w <- rep.int(1,length(y));
  treeMaster <- TreeRankRec(formula = formula, data = mf, bestresponse  = bestresponse,growing = growing, LeafRank =LeafRank,varsplit=varsplit)
   
#  if (prcSplitVar<1){
    #idxVar <- c(which(colnames(x)%in%rn),sample(inputs,ceiling(length(inputs)*prcSplitVar)));
  #}

  #if ((prcSplitData<1) || dataRepl){
    #idxData <-sample(nrow(data),ceiling(nrow(data)*prcSplitData),replace=dataRepl)
  #}
  
#  colname.list <-c(names(data),".TRvar");
#  cols <-lapply(colname.list,function(x)numeric(0));
#  names(cols)<-colname.list;
#  dataN <- do.call("data.frame",cols);

    #Build the Master TreeRank
 #tmpTRdata <- mf[idxData,idxVar];
  #tmpTRy <- y[idxData];
 
  #Compute the TreeRank object
  treeMaster$terms <- Terms;
  treeMaster$formula <- formula;
  treeMaster$call <-call;
  treeMaster$data <-data;
  ret <- treeMaster;
  class(ret) <- c("TR_TreeRank")
  #Do the pruning
  if (nfcv>1){
      print("pruning master tree...");
      tmpTR <- function(...){TreeRank(LeafRank = LeafRank,growing = growing,...)}
      ret <- pruneCV(ret,formula,mf,bestresponse,DTalgo = tmpTR,nfcv);
    }
  ret
}



print.TR_TreeRank <- function(x,...){
object <- x;
  if(!inherits(object,"TR_TreeRank"))
    stop("objectnot of class TR_TreeRank");

  id <- object$root;
  nodestack<-id;
  cat("TreeRank tree\n   id) #pos:#neg dauc score\n\n");

  while(length(nodestack)>0){
    id <- nodestack[[1]];
    nodestack <-nodestack[-1];
    s <- "";
    s <- paste(cat(rep(' ',2*object$depth[id])),id,"| ",object$pcount[id],":",
               object$ncount[id]," ",format(object$ldauc[id],digits=3)," ",format(object$score[id],digits=3),sep="")
    
    if (!(object$isleaf[id])){
      nodestack <- c(object$kidslist[[id]][[1]],object$kidslist[[id]][[2]],nodestack);
    }else{s<- paste(s,"*");}
    cat(paste(s,"\n"));
  }
}
    
    
