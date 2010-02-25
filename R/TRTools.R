#################################################################################
#
# Many common auxiliary functions for TreeRank and LeafRank:
#   - Pruning
#   - ROC computing
#   - Tree information computing
#
#################################################################################

#################################################################################
#
# checkImpure(weights,y) : check if at least two labels with weights != 0
#                           exists in data y
#
#################################################################################

checkImpure <- function(w,y){
  if (length(y)<2)
    return(FALSE);
  ret <- FALSE
	ytmp <- y[w >0]	
  def <- ytmp[1]
	for (i in 2:length(ytmp)){
		if (ytmp[i] != def) 
		{
			ret <- TRUE
			break;
		}
	}
	ret
}



#################################################################################
#
# growing_ctrl(minsplit,maxdepth,mincrit): growing tree control function
#
#
#################################################################################

growing_ctrl <- function(minsplit = 50,maxdepth =10, mincrit = 0) {

    ret <- list(minsplit = minsplit, maxdepth = maxdepth, mincrit= mincrit)
    return(ret)
}



#################################################################################
#################################################################################
##
## ROC manipulation
##
#################################################################################
#################################################################################

#################################################################################
#
#  auc(roc) :
#  Compute the auc from points list describing the roc
#
#################################################################################


auc <- function(roc){
  l = nrow(roc);
  auc <- sum((roc[2:l,1] - roc[1:(l-1),1])*(roc[2:l,2] +roc[1:(l-1),2]))/2
  auc
}

#################################################################################V1.5
#
#  getROC(tree,data) :
#  Compute the roc. If data is null, return the learning roc.
#
#################################################################################


getROC <- function(obj,data =NULL){
  UseMethod("getROC");
}


getPREC <- function(obj,data =NULL){
  UseMethod("getPREC");
}


getCurves <- function(obj,data=NULL){
UseMethod("getCurves");
}

getCurves.TR_TreeRank <- function(obj,data=NULL){
tree <- obj;
  if (!inherits(tree,"TR_TreeRank"))
    stop("object not of class TR_TreeRank");
  
  if(is.null(data))
data <- obj$data;
  score <- predict(tree,data);
  #Extract the name of the response variable
  resname <-all.vars(terms(tree))[[attr(terms(tree),"response")]];
  response <- data[[resname]];
  #compute the predicted score for data
   prec <- getPRECfromScore(score,response,tree$bestresponse);
  roc<- getROCfromScore(score,response,tree$bestresponse);
  c(list(roc),list(prec))
}          



getPREC.TR_TreeRank <- function(obj, data=NULL){
tree <- obj;
  if (!inherits(tree,"TR_TreeRank"))
    stop("object not of class TR_TreeRank");
  
  if(is.null(data)){
  data <- tree$data
  }
  score <- predict(tree,data);
  #Extract the name of the response variable
  resname <-all.vars(terms(tree))[[attr(terms(tree),"response")]];
  response <- data[[resname]];
  #compute the predicted score for data
   res <- getPRECfromScore(score,response,tree$bestresponse);
  res
}          

getROC.TR_TreeRank <- function(obj, data=NULL){
tree <- obj;
  if (!inherits(tree,"TR_TreeRank"))
    stop("object not of class TR_TreeRank");
  
  if(is.null(data)){
    return(matrix(c(tree$lalpha[tree$leafOrdered],1,tree$lbeta[tree$leafOrdered],1),length(tree$leafOrdered)+1))
  }
  score <- predict(tree,data);
  #Extract the name of the response variable
  resname <-all.vars(terms(tree))[[attr(terms(tree),"response")]];
  response <- data[[resname]];
  #compute the predicted score for data

  res <- getROCfromScore(score,response,tree$bestresponse);
  res
}

getPREC.TR_forest <- function(obj,data=NULL){
tree <- obj;
  if (!(inherits(tree,"TR_forest")))
    stop("objectnot of class TR_forest");
  if (is.null(data))
    data <- tree$forest[[1]]$data
    score <- predict(tree,data);
    resname<- all.vars(terms(tree$forest[[1]]))[[attr(terms(tree$forest[[1]]),"response")]];
    response <- data[[resname]];
      roc <- getPRECfromScore(score,response,tree$forest[[1]]$bestresponse);
     roc
}

getROC.TR_forest <- function(obj,data=NULL){
tree <- obj;
  if (!(inherits(tree,"TR_forest")))
    stop("objectnot of class TR_forest");
  if (is.null(data))
    data <- tree$forest[[1]]$data
    score <- predict(tree,data);
    resname<- all.vars(terms(tree$forest[[1]]))[[attr(terms(tree$forest[[1]]),"response")]];
    response <- data[[resname]];

    roc <- getROCfromScore(score,response,tree$forest[[1]]$bestresponse);
  roc
}




getCurves.TR_forest <- function(obj,data=NULL){
tree <- obj;
  if (!(inherits(tree,"TR_forest")))
    stop("objectnot of class TR_forest");
  if (is.null(data))
    data <- tree$forest[[1]]$data
    score <- predict(tree,data);
    resname<- all.vars(terms(tree$forest[[1]]))[[attr(terms(tree$forest[[1]]),"response")]];
    response <- data[[resname]];

    roc <- getROCfromScore(score,response,tree$forest[[1]]$bestresponse);
    prec <- getPRECfromScore(score,response,tree$forest[[1]]$bestresponse);
  c(list(roc),list(prec))
}



#################################################################################
#
#  getROCfromScore(score,y,bestresponse) :
#  Compute the roc from the score list.
#
#################################################################################




getROCfromScore<- function (score, y, bestresponse){
	orderedIndex <- order(-score)
        alphaList <-  0
	betaList <- 0
	curScore <- score[[orderedIndex[1]]]
	curAlpha <- 0
	curBeta <- 0
	pcount <- 0
	ncount <- 0
	for (i in 1:length(orderedIndex)){
		if(curScore != score[[orderedIndex[i]]])
		{
		   alphaList <- c(alphaList,curAlpha)
		   betaList <- c(betaList,curBeta)
                   		   curScore <- score[[orderedIndex[i]]]	
		}
	if (y[orderedIndex[i]] == bestresponse)
	{
		curBeta <- curBeta +1
		pcount <- pcount +1	
	}
	else {
		curAlpha <- curAlpha +1
		ncount <- ncount +1 	
	}
	}
	alphaList <- c(alphaList,curAlpha)
	betaList <- c(betaList,curBeta)
	alphaList <- alphaList/ncount
	betaList <- betaList/pcount
	matrix(c(alphaList,betaList),length(alphaList))
}

getPRECfromScore <- function(score,y,bestresponse){
  orderedIndex <- order(-score);
  alphaList <-  0
  betaList <- 0
  curScore <- score[[orderedIndex[1]]]
  curAlpha <- 0
  curBeta <- 0
  pcount <- 0
  ncount <- 0
  for (i in 1:length(orderedIndex)){
    if(curScore != score[[orderedIndex[i]]])
      {
        alphaList <- c(alphaList,curAlpha)
        betaList <- c(betaList,curBeta)
        curScore <- score[[orderedIndex[i]]]	
      }
    if (y[orderedIndex[i]] == bestresponse)
      {
        curBeta <- curBeta +1
        pcount <- pcount +1	
      }
    else {
      curAlpha <- curAlpha +1
      ncount <- ncount +1 	
    }
  }

  	
  alphaList <- c(alphaList[-1],curAlpha)
  betaList <- c(betaList[-1],curBeta)
  prec <- c(0,betaList/pcount,1);
  rec <-  c(1,betaList/(alphaList+betaList),0);
  matrix(c(prec,rec),length(prec))
}
#prec TP/TP+FP
# recall TP/TP+FN







#################################################################################
#################################################################################v1
##
##  Tree manipulation
##
#################################################################################
#################################################################################

rpart2TR <- function(tree,bestresponse){
  nbnode <- nrow(tree$frame);

  frame <- tree$frame;
  splits <- tree$splits
  sortIdx <- order(as.numeric(rownames(frame)),decreasing=TRUE);
  sortSplit <- as.numeric(rownames(frame))[which(frame$var != "<leaf>")];
  sortSplit <- order(sortSplit);
  newnodestack <- list();
  pcInit <- sum(tree$y == bestresponse);
  ncInit <- length(tree$y) -pcInit;
  pInit <- pcInit/(pcInit+ncInit);
  idsplit <- length(sortSplit);

  kidslist <- list();
  parentslist <- array();
  split <- list();
  isleaf <- array();
  ldauc <- array();
  lalpha <- array();
  lbeta <- array();
  ncount <- ncInit;
  pcount <- pcInit;
  d <-c(1);
  depth <-array();
  sign <- array();
  parentslist[1] <- NA;
  if (nbnode <2){
    kidslist <- list();
    nbNode <- 1;
    # return(party(ret,data = tree$model,terms=tree$terms));
  }
else{

  for (i in sortIdx){
    d <- floor(log2(as.numeric(rownames(frame)[[i]])))+1;
    if (frame$var[[i]] == "<leaf>"){
      pcount[i] <- sum(tree$where==i & tree$y == bestresponse);
      ncount[i] <- sum(tree$where==i & tree$y !=bestresponse);
      ldauc[i]<-0;
      depth[i]<-d;
      split[i]<-NA;
      kidslist[i]<-NA;
      isleaf[i]<-TRUE;
    }
    else{
      if ((length(newnodestack)<2) || (idsplit<1)){
        stop("Error converting rpart tree to TR tree");
      }

      lnode <- newnodestack[[2]];
      rnode <- newnodestack[[1]];
      pcount[i]<-pc <- pcount[lnode]+pcount[rnode];
      ncount[i]<-nc <- ncount[lnode] + ncount[rnode];
      npc <- pcount[lnode];
      nnc <- ncount[lnode];
      ldauc[i]<- dauc <- -((pc/pcInit)*(nnc/ncInit)-(nc/ncInit)*(npc/pcInit))/2
      sign[i] <- -1;
      if (dauc<0){
        ldauc[i]<- dauc <- -dauc;
        sign[i] <- 1;
      }
      idVar <- which(names(tree$model)==frame$var[[i]])
      tmpname <- as.character(frame[i,"var"]);

      if (attributes(tree$terms)$dataClasses[tmpname] == "numeric"){
        br <- as.numeric(splits[,"index"][[sortSplit[[idsplit]]]])
      split[i]<-list(list(idVar=idVar,name=tmpname,breaks = br,type=0));
	}
      else
	{
 	br <- splits[,"index"][[sortSplit[[idsplit]]]];
	lev <- which(tree$csplit[br,] == 3);
	levVal <- attr(tree,"xlevels")[tmpname][[1]];
	# print(lev)
	lev <- strsplit(levVal, " ")[lev];
      	split[i]<-list(list(idVar=idVar,name=tmpname,breaks = lev,type=1));
	}
      #      sp <- partysplit(idVar, breaks = splits[,"index"][[sortSplit[[idsplit]]]],right=FALSE);
      kidslist[i] = list(c(lnode,rnode));
      parentslist[lnode]<-parentslist[rnode]<-i;
      depth[i]<-d;
      isleaf[i]<-FALSE;
      newnodestack <- newnodestack[c(-1,-2)]
      idsplit <- idsplit-1;
    }
    newnodestack <- c(newnodestack,i);
  }
  if ((length(newnodestack)>1) || (idsplit>0)){
    stop("Error converting rpart tree to TR");
  }
}
   ret <- list();
   ret$nbNode <- max(sortIdx);
   ret$root <- 1L;
   ret$nodes <- 1:ret$nbNode;
   ret$inner <- ret$nodes[!(isleaf)]
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
  ret$split <-  split;

#  ret <- data = tree$model,terms=tree$terms);
  ret;
}



#################################################################################
#
# ComputeTreeInfo(tree,pcount,ncount) : Compute auc, delta-auc, TP and FP rate
#         if pcount and ncount are both not null, the computation is done wrt these  lists
#           and not the internal pcount and ncount tree lists.
#
#################################################################################



ComputeTreeInfo <- function(tree,pcount = NULL, ncount = NULL){
  nbnode <- tree$nbNode;
  isleaf <- tree$isleaf;
  lalpha <- rep(0,nbnode);
  lbeta <- rep(0,nbnode);
  lauc <- rep(0,nbnode);
  ldauc <- rep(0,nbnode);
  if (is.null(pcount)){
    pcount <- tree$pcount;
    ncount <- tree$ncount;
  }
  pcinit <- pcount[[1]];
  ncinit <- ncount[[1]];

  nodestack <- tree$root;
  ldauc[[tree$root]] <- 0.;
  lauc[[tree$root]] <- 0.5;
  auc <- 0.5;
  dauc <- 0;

  # Breadth-first exploration of the tree
  while(length(nodestack)>0){
    dauc <- 0
    curnode <- nodestack[[1]];
    nodestack <- nodestack[-1];
     if (!(isleaf[[curnode]])){
        lkid<- tree$kidslist[[curnode]][[1]];
        rkid<- tree$kidslist[[curnode]][[2]];
      alphaR <- 1;
      betaR <- 1;

    #Computing new alpha and beta, cumulative TP and FP rate, auc and dauc
        nalpha <- lalpha[[curnode]] + ncount[[lkid]]/ncinit;
      nbeta <-  lbeta[[curnode]] + pcount[[lkid]]/pcinit;
      lalpha[[lkid]] <-lalpha[[curnode]];
      lbeta[[lkid]] <- lbeta[[curnode]];
      lalpha[[rkid]] <- nalpha;
      lbeta[[rkid]] <- nbeta;
      dauc <-  -((pcount[[curnode]]/pcinit)*(ncount[[lkid]]/ncinit) -(pcount[[lkid]]/pcinit)*(ncount[[curnode]]/ncinit))/2
      auc <-auc + dauc;
      lauc[[curnode]] <- auc;
      ldauc[[curnode]] <- dauc;
      nodestack <- c(nodestack, lkid,rkid);
  }

  }
  if (inherits(tree,"TR_LRCart")){
    sign <- unlist(lapply(tree$sign,function(x){if (is.null(x$sign))
                                                            return(0);
                                                            return(-1*x$sign);
                                                          }))
    ldauc <- sign*ldauc;
    return(list(lalpha = lalpha, lbeta = lbeta, lauc = lauc, ldauc = ldauc,sign=sign));

  }                       
  return(list(lalpha = lalpha, lbeta = lbeta, lauc = lauc, ldauc = ldauc));
}


#################################################################################
#
# subTreeRank : extract a subTree from the list of the new leaves and a tree
#
#################################################################################
 

subTreeRank <- function(tree,listnode)
{
   id <- tree$root;
   nextnode <- 1;
   kidslist <- list();
   parentslist <- array();
   LRList <- list();
   isleaf <-array();
   ldauc <- array();
   lalpha <- array();
   lbeta <- array();
   ncount <- array();
   pcount <- array();
   depth <- c(1);
   nodestack <- list(id);
   nodestackO <- list(id);
   curscore <- 1;
   nodeorder <- array();
   parentslist[1] <- NA;
   while(length(nodestackO)>0){
     idO <- nodestackO[[1]];
     nodestackO <- nodestackO[-1];
     id <- nodestack[[1]];
     nodestack <- nodestack[-1];
     isleaf[id]<-TRUE;
     lalpha[id]<-tree$lalpha[idO];
     lbeta[id] <-tree$lbeta[idO];
     ncount[id] <-tree$ncount[idO];
     pcount[id]<-tree$pcount[idO];
     depth[id] <- tree$depth[idO];
     nodeorder[id] <- curscore;
     curscore <- curscore +1;
     if ((tree$isleaf[idO]) || (idO %in% listnode)){
       kidslist[id]<-NA;
       ldauc[id]<-0;
     	LRList[id]<-NA;  
       next;
     }
     nodestack <- c(nextnode+1,nextnode+2,nodestack);
     nodestackO <- c(tree$kids[[idO]][[1]],
                     tree$kids[[idO]][[2]],nodestackO);
     parentslist[nextnode+1] <- id;
     parentslist[nextnode+2] <- id;
     kidslist[id] <- list(c(nextnode+1,nextnode+2));
     LRList[id] <-tree$LRList[idO];
     nextnode <- nextnode +2;
     ldauc[id]<-tree$ldauc[idO];
     isleaf[id]<- FALSE;
   }

   ret <- tree;
   ret$nbNode <- nextnode;
   score <- array(0,ret$nbNode);
   leaf <- which(isleaf);
   leafOrdered <- order(nodeorder);
   leafOrdered <- leafOrdered[leafOrdered %in% leaf];
   nbLeaf <- length(leafOrdered); 
   for (i in   0:(nbLeaf-1))
      score[leafOrdered[i+1]] <- ((nbLeaf -i)/nbLeaf)
   ret$score <- score;
   ret$nodes <- (1:ret$nbNode);
   ret$leafOrdered <- leafOrdered;
   ret$inner <- ret$nodes[!(isleaf)];
   ret$root <- 1L;
   ret$parentslist <- parentslist;
   ret$kidslist <- kidslist;
   ret$lalpha <- lalpha;
   ret$lbeta <- lbeta;
   ret$ldauc <- ldauc;
   ret$isleaf <- isleaf;
   ret$pcount <- pcount;
   ret$ncount <- ncount;
   ret$depth <- depth;
   ret$LRList <- LRList;
   ret$data <- tree$data;
   return(ret); 
          
}



#################################################################################
#################################################################################
##
##  pruning stuff
##
#################################################################################
#################################################################################


###############################################################
#
# pruneInfo(tree) : compute the informations regarding the pruning process :
#  return : sorted list by pruning level of : node ids of pruned nodes,
#                                             alpha complexity list
#                                             # of terminal nodes
#
###############################################################

getMiss <- function(pcount,ncount){
  if (pcount>ncount){
    return(ncount)
  }
  return(pcount)
}

pruneInfo <- function(tree){
  kidslist <- tree$kidslist
  listnbleaves <- as.double(tree$isleaf)
  cumdauc <- tree$ldauc;
  miss <- array(0,length(tree$pcount))
  pc <- tree$pcount
  nc <- tree$ncount
  isleaftmp <- tree$isleaf;
  listnode <- tree$nodes[tree$isleaf];

  #Compute for each node the cumulative delta-auc of the node + his kids and the number of leaves
  # by a reverse breadth-first exploration
  while (length(listnode)>0){
    id <- listnode[[1]];
    listnode <- listnode[-1];
    if (tree$isleaf[[id]]){
      miss[[id]] <- getMiss(pc[[id]],nc[[id]]);}
    pred <- tree$parentslist[[id]];
    isleaftmp <- c(isleaftmp,id);
    if((!(is.na(pred))) && (pred>0)){
      listnbleaves[[pred]] <- listnbleaves[[pred]]+listnbleaves[[id]];
      cumdauc[[pred]] <- cumdauc[[pred]]+ cumdauc[[id]];
      miss[[pred]] <- miss[[pred]]+miss[[id]];
      if (min(kidslist[[pred]] %in% isleaftmp) == 1)
        listnode <- c(listnode,pred);
    }
  }
  isleaftmp <- tree$isleaf;
  idpruned <- vector(mode="double",0);
  alphalist <- vector(mode="double",0);
  ntermnodes <- vector(mode="double",0);
  listnode <-tree$nodes;
  activenode <- listnode;
  activenode <-activenode[ -which(tree$isleaf)]
  
  while(length(activenode > 0)){
   #Find  next nodes to be pruned wrt (miss-classified / complexity)

#    alpha <- pmax(0,1-cumdauc[activenode])/ pmax(1,listnbleaves[activenode]-1);
#    missloc <- unlist(sapply(activenode,function(x){(getMiss(pc[[x]],nc[[x]])-miss[[x]])/(pc[[x]]+nc[[x]])}));
    
 #   alpha <- missloc/pmax(1,listnbleaves[activenode]-1);
    alpha <- pmax(0,cumdauc[activenode])/pmax(1,listnbleaves[activenode]-1);
    idprune <- which.min(alpha);
    bestalpha <- alpha[[idprune]];
    idprune <- which(alpha == bestalpha)
    idprune <- activenode[idprune];
    kidstmp <- unlist(tree$kidslist[idprune]);

    #desactivate all kids of  next nodes to be pruned
    while (length(kidstmp)> 0){
      idtmp <- kidstmp[[1]];
      kidstmp <- kidstmp[-1];
      if (idtmp %in% activenode){
        isleaftmp[[idtmp]] <- TRUE;
        ltmp <- kidslist[[idtmp]];
        if (length(ltmp)>0){
          kidstmp <- c(kidstmp,ltmp[which(!isleaftmp[ltmp])]);
          activenode <- activenode[-which(activenode ==idtmp)];
        }
      }
    }

    #Update the complexity and the delta-auc of remaining nodes
    for (idtmp in idprune){
      isleaftmp[idtmp] <- TRUE;
      if (idtmp %in% activenode){
        activenode <- activenode[-which(idtmp== activenode)];
        daucdiff <- cumdauc[[idtmp]];
        missdiff <- miss[[idtmp]]-getMiss(pc[[idtmp]],nc[[idtmp]])
        nbleavesdiff <- listnbleaves[[idtmp]] -1;
        idupdate <- idtmp;
        while ((!(is.na(idupdate))) && (idupdate>0)){
          listnbleaves[[idupdate]] <- listnbleaves[[idupdate]]-nbleavesdiff;
          cumdauc[[idupdate]] <- cumdauc[[idupdate]] - daucdiff;
          miss[[idupdate]] <- miss[[idupdate]]-missdiff
          idupdate <- tree$parentslist[[idupdate]];
        }
      }
    }
    #Update return lists information
    idpruned <- c(idpruned,list(idprune));
    alphalist <- c(alphalist,bestalpha);
    ntermnodes <- c(ntermnodes,listnbleaves[[1]]);
  }
  return(list(idpruned = idpruned, alphalist = alphalist, ntermnodes = ntermnodes));
}


###############################################################
#v1
# pruneCV(tree,formula,data, bestresponse, DTalgo = LeafRankCart,pruning)
#
# Main pruning function :
#      DTalgo : function used to grow the tree
#
###############################################################


pruneCV <- function(tree,formula, data, bestresponse, DTalgo = LRCart, nfcv){

  if (nfcv <=1)
    return(tree);
  # if not enough data or tree too small, don't prune
  if (nrow(data)/nfcv <5)return(tree)
  if (tree$nbNode<4)
    return(tree)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data"),names(mf),0)
  mf <- mf[c(1,m)]
  mf$drop.unused.levels <-FALSE;
  mf[[1]] <- as.name("model.frame")
  mf <-eval(mf,parent.frame());
  response <- model.response(mf);
  rn <- names(mf)[1]
  response <- model.response(mf);
  x <- mf[,colnames(mf) != "(weights)"]
  inputs <- which(!(colnames(x) %in% rn));
  if (is.null(weights)) weights <-rep.int(1,length(response));

  #Get pruning information on  current tree
  treePI <- pruneInfo(tree);

  #Build folds for cross validation
  index <- array(1:length(data[,1]));
  nc = array(dim=nfcv);
  
#  if (pruning$strat){
 #   indexpos <- sample(array(which(response == bestresponse)));
  # indexneg <- sample(array(index[-indexpos]));
   # nfcpossize <- floor(length(indexpos)/nfcv);
   # nfcnegsize <- floor(length(indexneg)/nfcv);
    #for (i in 1:(nfcv-1)){
     # nc[i] <- list(c(indexpos[(nfcpossize*(i-1)+1):(nfcpossize*i)],indexneg[(nfcnegsize*(i-1)+1):(nfcnegsize*i)]));
    #}
    #nc[nfcv] <- list(c(indexpos[(nfcpossize*(nfcv-1)+1):length(indexpos)],indexneg[(nfcnegsize*(nfcv-1)+1):length(indexneg)]));
  #}
  #else{
    nfcsize <- length(data[,1]) / nfcv;
    index <- sample(length(data[,1]));
    for (i in 1:(nfcv-1))
      nc[i] <- list(index[(nfcsize*(i-1)+1):(nfcsize*i)]);
    nc[nfcv] <- list(index[(nfcsize*(nfcv-1)+1):length(index)]);
  #}

#Do the cross validation
  if (length(treePI$alphalist)<2){
    avgalpha <- c(treePI$alphalist[[1]],1)
  }else
  {
    avgalpha <- c(0,sqrt(treePI$alphalist[2:length(treePI$alphalist)]*treePI$alphalist[1:(length(treePI$alphalist)-1)]),1);
  }
  aucnfcv = array(dim=c(nfcv,length(avgalpha)));
  for (i in 1:nfcv){
    print(c("--- cv learn pass: ",i));
    dtmp<- data[-nc[[i]],]
    rtmp <- response[-nc[[i]]];
    treetmp <- DTalgo(formula = formula,data=dtmp,
                      bestresponse= bestresponse);
    aucnfcv[i,] <- aucTestInfo(treetmp,dtmp,rtmp,bestresponse,avgalpha)
    
  }
  auc <- array();
  for (i in 1:length(avgalpha)){
    auc[i] <-mean(aucnfcv[,i]);
  }
  #aucmax <- max(auc);
  aucmax <- max(auc)
  if (is.na(aucmax))
    return(ret);
  prlevel <- max(which(auc==aucmax));
  ret <- tree;
  if ((length(prlevel)>0) &&(prlevel>1)){
    if(length(unlist(treePI$idpruned[1:(prlevel-1)]))>0){
      ret <- subTreeRank(tree,unlist(treePI$idpruned[1:(prlevel-1)]))
    }
  }
  ret$unpruned <- tree;
  ret;
}

aucTestInfo <- function(tree, data,y,bestresponse,avgalpha){
  if (!(checkImpure(rep(1,length(y)),y))){
    return(array(dim=length(avgalpha)));}
    
  res <- predict(tree,data,type="node");
  treePI <- pruneInfo(tree);
    

  idlist <- tree$nodes;
  ncount <- array(0,dim=length(idlist));
  pcount <- array(0,dim=length(idlist));
  
  #Compute pos and neg count for tree nodes
  for (i in 1:length(res))
    ifelse(y[[i]] == bestresponse,pcount[res[i]] <- pcount[res[i]] +1, ncount[res[i]] <- ncount[res[i]] +1);
  nodestack <- (1:tree$nbNode)[tree$isleaf];
  treated <- vector(mode="integer");
  while (length(nodestack) >0){
    idcur <- nodestack[[1]];
    nodestack <- nodestack[-1];
    idpar <- tree$parentslist[idcur];
    pcount[idpar] <- pcount[idpar] + pcount[idcur];
    ncount[idpar] <- ncount[idpar] + ncount[idcur];
    treated <- c(treated,idcur);
    if ((!(is.na(idpar)))&&(idpar>0))
      if (min(tree$kidslist[[idpar]] %in% treated) == 1)
        nodestack <- c(nodestack,idpar);
  }
  aucinfo <- ComputeTreeInfo(tree,pcount,ncount);
  misscount <- unlist(sapply(1:tree$nbNode,function(x){getMiss(tree$pcount[[x]],tree$ncount[[x]])}));
  auctestlist <- array(dim=length(avgalpha));
  for (i in 1:length(avgalpha)){
    itmp <- which(treePI$alphalist <=avgalpha[i]);
    if (length(itmp)>0){
      idtmp <- unlist(treePI$idpruned[1:max(itmp)]);
      auctestlist[i] <- getAUCsubTree(tree,idtmp,aucinfo$ldauc);
 #     auctestlist[i]<-getMissSubTree(tree,idtmp,misscount);
    }else{
      auctestlist[i] <- 0.5+sum(aucinfo$ldauc);
  #    auctestlist[i] <- getMissSubTree(tree,nodeids(tree,terminal=TRUE),misscount);
        }
  }
  
  return(auctestlist);
}


getMissSubTree <- function(tree,listnode,misscount){
  nodestack<- 1;
  miss <-0;
  while(length(nodestack)>0){
    cur <- nodestack[[1]]
    nodestack <- nodestack[-1];
    if (!(cur %in% listnode)){
      kids <- tree$kidslist[[cur]];
      nodestack <- c(nodestack,unlist(kids));
    }
    else{
      miss <- miss + misscount[[cur]]
    }
  }
  return(miss/(tree$pcount[[1]]+tree$ncount[[1]]));
}


getAUCsubTree <- function(tree,listnode,ldauc){
  nodestack <- 1;
  dauc <-0.5;
  while(length(nodestack)>0){
    cur <- nodestack[[1]]
    nodestack <- nodestack[-1]
    if (!(cur %in% listnode)){
      dauc <- dauc+ldauc[[cur]];
      kids <-tree$kidslist[[cur]];
      if (!(sum(is.na(kids)))){
        nodestack <- c(nodestack,unlist(kids));
      }
    }
  }
  return(dauc)
}
########################################################
#
#  importance Variable
#
########################################################

varImportance <- function(obj,norm=TRUE){
  UseMethod("varImportance");
}

varImportance.default <- function(obj,norm=TRUE){
 NULL;
}

varImportance.TR_forest <- function(obj,norm=TRUE){
  tree <- obj$forest[[1]]
  nr <- attributes(tree$terms)$term.labels;
  res <- array(0,length(nr))
  names(res) <- nr
  for (i in obj$forest){

    restmp <- varImportance(i,FALSE);
    if (is.null(restmp)){return(NULL);}
    for (j in 1:length(restmp)){
      id <- which(names(restmp[j])==names(res));
      res[id] <- res[id]+restmp[j];
    }
  }
  if (norm) ret <- res/max(res)
  else ret <- res;
  ret;
}
    

varImportance.TR_TreeRank<- function(obj,norm=TRUE){
  nr <- attributes(obj$terms)$term.labels;
  res <- array(0,length(nr))
  names(res) <- nr
  
  for (i in obj$inner){

    restmp <- varImportance(getClassifier(obj,i),FALSE)
    if (is.null(restmp)){return(NULL);}
    for (j in 1:length(restmp)){
      id <- which(names(restmp[j])==names(res));
      res[id] <- res[id]+restmp[j]*((obj$pcount[i]*obj$ncount[i])/(obj$pcount[obj$root]*obj$ncount[obj$root]))^2;
    }
  }
  if (norm) ret <- res/max(res)
  else ret <- res;
  ret;
}


plotROC <- function( rocs,colorlist = NULL,points=NULL){
  if (is.null(colorlist)){  colorlist <- list("black","blue","red","green","yellow")}
  plot(function(x){x},0:1, col ="black", xlab = "FP rate", ylab = "TP rate")
  if (length(rocs)>0){
    for (i in 1:length(rocs)){
      curr <- rocs[[i]];
      if (is.list(curr)){
        for(j in 1:length(curr)){
          par(new = T);
          if (nrow(curr[[j]])<100){type <- "b"}
          else{type <-"l"}
          plot(curr[[j]], type=type, col= colorlist[[(i-1) %% length(colorlist) + 1]], xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
        }
      }else{
        par(new=T);
        if (nrow(rocs[[i]])<100){type <- "b"}
        else{type <-"l"}
        plot(rocs[[i]], type=type, col= colorlist[[(i-1) %% length(colorlist) + 1]], xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
      }
    }
    
    if (!(is.null(points))){
      points<-c(list(c(0,0)),points,list(c(1,1)))
      tmpm<- matrix(unlist(points),nrow=2);
      par(new=T);
      
      plot(x=tmpm[1,],y=tmpm[2,],xlab="",ylab="",xlim=c(0,1),ylim=c(0,1));
      
    }
    
  }
}


##############################################################
#
# Bagging
#
##############################################################


TreeRankBagging <- function(forest){
  
  res <- list(forest= forest, ntree = length(forest));
  class(res) <- "TR_forest"
  res
}



varDep <- function(obj,data,varx,vary,vminx=min(data[varx]),vmaxx=max(data[varx]),vminy=min(data[vary]),vmaxy=max(data[vary]),subdivx=100,subdivy=subdivx){

   nbex <- nrow(data);
nc <- which(colnames(data) %in%c(varx,vary));
seqx <-  seq(vminx,vmaxx,length.out=subdivx);
seqy <- seq(vminy,vmaxy,length.out=subdivy);
ret <- matrix(nrow=subdivx,ncol=subdivy,dimnames=list(seqx,seqy));
  for (i in 1:length(seqx)){
    for(j in 1:length(seqy)){
   data[,nc[1]] <- seqx[i];
   data[,nc[2]] <- seqy[j];
   pre <-  sum(predict(obj,data)/nbex);
   ret[i,j] <-pre;
  }
  }
ret;
}


TwoSample<- function(x,y,split=60,TRalgo=TreeRank,alpha= 5){
  xn <-nrow(x);
  yn <-nrow(y);
  split <- split/100;
  alpha <- 1-(alpha/100);
  xset <- data.frame(x,"twTRclass"=rep(1,xn));
  yset <- data.frame(y,"twTRclass"=rep(-1,yn));
  xsplit <- ceiling(xn*split);
  ysplit <- ceiling(yn*split);
  idx <- sample(xn,xsplit)
  idy <- sample(yn,ysplit);
  train <- rbind(xset[idx,],yset[idy,]);
  test <- rbind(xset[-idx,],yset[-idy,]);
  tree <- TRalgo(twTRclass ~.,train,bestresponse=1);

  
  score <- predict(tree,test);
  dft <- data.frame(score = score,class=as.factor(test[,"twTRclass"]))
  wtest <- wilcox_test(score~class,dft,conf.level=alpha,conf.int=TRUE);
  return(list(wtest = wtest, tree= tree,train=train,test=test,dft=dft))

}
