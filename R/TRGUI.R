#########################################################
#
# GUI functions, largely inspired by igraph tkiplot function
#
#
#########################################################




#########################################################
# trg : Environment var, fields:
# tree, coords, treeCanvas,minNodeSize,maxNodeSize,minNodePad

if (!exists(".TRgui.env")){
  .TRgui.env <- new.env();
  assign(".next",1,.TRgui.env);
}


#########################################################
# LeafRank functions definition
#########################################################

cartOpt <- list(
                list(name="Minimum Split",default=50,type="entry",optName="minsplit"),

                list(name="Min. Crit.",default=0,type="entry",optName="mincrit"),
                list(name="Maximum Depth",default=10,type="entry",optName="maxdepth"),
                list(name="n-fold",default=0,type="entry",optName="nfcv"),
                list(name="%Var. Split",default=100,type="entry",optName="varsplit"),
                list(name="%Data Split",default=100,type="entry",optName="datasplit")
                );
LRCartopt2cmd <-function(obj){
  return(paste("maxdepth=",obj$maxdepth,
               ",minsplit=",obj$minsplit,",mincrit=",obj$mincrit,
               ",prcSplitVar=",(obj$varsplit/100),
               ",prcSplitData=",(obj$datasplit/100),
               ",nfcv=",obj$nfcv,sep=""));
}
LRsvmopt2cmd <- function(obj){
  tmp <- "";
  if (!(obj$auto)){
    if (obj$kernel == "polydot"){
      tmp <-  paste(",kpar=list(","scale=",obj$scale,",offset=",obj$offset,
    ",degre=",obj$sigma,")",sep="");
    }
    if(obj$kernel == "tanhdot"){
      tmp <-  paste(",kpar=list(","scale=",obj$scale,",offset=",obj$offset,
                    ")",sep="");
    }
    if (obj$kernel == "rbfdot"){
      tmp <-  paste(",kpar=list(sigma=",obj$sigma,")",sep="");
    }
  }
  
  return(paste("C=",obj$C,
               ",kernel=\"",obj$kernel,
               "\",prcSplitVar=",(obj$varsplit/100),
               ",prcSplitData=",(obj$datasplit/100),tmp,sep=""));
}

LRCart.def <- list(name="Cart",fun="LRCart",
                   opt=cartOpt,opt2cmd=LRCartopt2cmd);


LRsvm.def <- list(name="svm",fun="LRsvm",
                     opt=list(
                       list(name="C",default=1,type="entry",optName="C"),
                       list(name="Degre/Sigma",default=1,type="entry",optName="sigma"),
                       list(name="scale",default=1,type="entry",optName="scale"),
                       list(name="offset",default=1,type="entry",optName="offset"),
                       list(name="%Var. Split",default=100,type="entry",optName="varsplit"),
                       list(name="%Data Split",default=100,type="entry",optName="datasplit"),
                       list(name="Kernel",default=0,type="listbox",optName="kernel",choiceName=c("RBF","Poly.","TanHyp"),choice=c("rbfdot","polydot","tanhdot")),
                       list(name="Auto. Parameters",default=TRUE,type="check",optName="auto")),
                  opt2cmd=LRsvmopt2cmd);


#########################################################
#
# Main launching Frame:
#  composed by 3 frames  : infoFrame (dataset definition)
#                          LR frame (for LeafRank)
#                          TR frame (for TreeRank)
#
#########################################################


TRGui <- function(){

  top <- tktoplevel(height="200",width="100")
  mainGuiEnv <- new.env();
  tktitle(top) <- "TreeRank GUI"
  mainFrame <- tkframe(top,height="200",width="100",borderwidth=2);

  infoFrame <- tkframe(mainFrame,borderwidth=5,width=100);
  LRFrame <- tkframe(mainFrame,borderwidth=5,width=100,relief="groove");
  TRFrame <-tkframe(mainFrame,borderwidth=5,width=100,relief="groove");
  
  optTRFrame <- tkframe(TRFrame,borderwidth=2,width=100);
  optLRFrame <- tkframe(LRFrame,borderwidth=2,width=100);
  

  ## info Frame definition
  dataset <- tclVar("");
  resname<-tclVar("class");
  bestresponse <-tclVar("1");

  entry.dataset <- tkentry(infoFrame,width="20",textvariable=dataset);
  entry.resname <- tkentry(infoFrame,width="20",textvariable=resname);
  entry.bestresponse <- tkentry(infoFrame,width="5",textvariable=bestresponse);
  tkgrid(tklabel(infoFrame,text="Data set: "),entry.dataset);
  tkgrid(tklabel(infoFrame,text="Label name: "),entry.resname);
  tkgrid(tklabel(infoFrame,text="Best Label: "),entry.bestresponse);
  
  
  #TreeRank options frame
  dftOpt <- growing_ctrl();
  minOptTR <- tclVar(dftOpt$minsplit);
  maxOptTR <- tclVar(dftOpt$maxdepth);
  mincritOptTR <- tclVar(dftOpt$mincrit);
  forestOptTR<-tclVar(0);
  varSplitTR <- tclVar(100);
  dataSplitTR <- tclVar(100);
  nfcvOptTR <-tclVar("0");
  entry.minOptTR <-tkentry(optTRFrame,width="5",textvariable=minOptTR);
  entry.mincritOptTR <-tkentry(optTRFrame,width="5",textvariable=mincritOptTR);
  entry.maxOptTR <-tkentry(optTRFrame,width="5",textvariable=maxOptTR);
  entry.forestOptTR <- tkentry(optTRFrame,width=5,textvariable=forestOptTR);
  entry.varSplitTR <- tkentry(optTRFrame,width="5",textvariable=varSplitTR);
  entry.dataSplitTR <- tkentry(optTRFrame,width="5",textvariable= dataSplitTR);
  entry.nfcvOptTR <- tkentry(optTRFrame,width=5,textvariable=nfcvOptTR);

  tkgrid(tklabel(optTRFrame,text="Minimum Split"),entry.minOptTR,
         tklabel(optTRFrame,text="Min. Criteria"),entry.mincritOptTR);
  tkgrid(tklabel(optTRFrame,text="Maximum Depth"),entry.maxOptTR,
         tklabel(optTRFrame,text="Forest"),entry.forestOptTR);
  tkgrid(tklabel(optTRFrame,text="%Var. split"),entry.varSplitTR,
         tklabel(optTRFrame,text="%Data. split"),entry.dataSplitTR);
  tkgrid(tklabel(optTRFrame,text="n-fold Cross Validation "),entry.nfcvOptTR);
  
  tkgrid(tklabel(TRFrame,text="TreeRank Options"));
  tkgrid(optTRFrame);

  
  
  
  #LeafRank
  ## Scan the environment for new LeafRank functions definition
  LRlist <- ls(pattern="LR.*.def",name=globalenv());
  nameLR <- list(LRCart.def$name,LRsvm.def$name);
  funNameLR <- list(LRCart.def$fun,LRsvm.def$fun);
  optionLR <- list(LRCart.def$opt,LRsvm.def$opt);
  opt2cmdLR <- list(LRCart.def$opt2cmd,LRsvm.def$opt2cmd);
  lbLR <- tklistbox(LRFrame,height=3,width=10,selectmode="single");
  tkinsert(lbLR,"end",LRCart.def$name);
  tkinsert(lbLR,"end",LRsvm.def$name);
  for (x in LRlist){
    obj <- eval(parse(text=x),globalenv())
    nameLR <- c(nameLR,obj$name);
    funNameLR <- c(funNameLR,obj$fun);
    optionLR <- c(optionLR,list(obj$opt));
    opt2cmdLR <- c(opt2cmdLR,list(obj$opt2cmd));
    tkinsert(lbLR,"end",obj$name)
  }

  tkselection.set(lbLR,0);
  assign("curLeafRank",1,mainGuiEnv);
  tkgrid(tklabel(LRFrame,text="LeafRank: "),lbLR);

  assign("listVarLR",list(),mainGuiEnv);
  assign("listWidgetLR",list(),mainGuiEnv);

  ## LR Frame building  
  buildLRopt <- function(){
    listVarLR <- list()
    id <- as.numeric(tkcurselection(lbLR))+1;
    if (length(id)<1)
      id <- 1;
    listOptLR <- optionLR[[id]];
    listSlaves <- tclvalue(tkgrid.slaves(optLRFrame));
    if (listSlaves[[1]] !=""){
      sapply(unlist(strsplit(listSlaves,split=" ")),tkgrid.remove)
    }
    listWidgetLR <- list()
    listGrid <- list();
    for (i in 1:length(listOptLR)){
      if (listOptLR[[i]]$type=="entry"){
        tmp <- tclVar(listOptLR[[i]]$default)
        listVarLR <- c(listVarLR,list(tmp))
        listWidgetLR<-c(listWidgetLR,list(tkentry(optLRFrame,width="5",textvariable=tmp)))
      }
      if (listOptLR[[i]]$type=="check"){
        tmp <- tclVar(listOptLR[[i]]$default);
        listVarLR <- c(listVarLR,list(tmp))
        listWidgetLR<-c(listWidgetLR,list(tkcheckbutton(optLRFrame,variable=tmp)))
      }
      if (listOptLR[[i]]$type=="listbox"){
        tmplb <- tklistbox(optLRFrame,height=3,width=6,selectmode="single");
        for (na in listOptLR[[i]]$choiceName){
          tkinsert(tmplb,"end",na);
        }
        tkselection.set(tmplb,listOptLR[[i]]$default);
        listVarLR <- c(listVarLR,list(tmplb));
        listWidgetLR<-c(listWidgetLR,list(tmplb))

      }
      tmpW <- listWidgetLR[[length(listWidgetLR)]]
    }
    for (i in (1:(length(listWidgetLR)/2))){
      tkgrid(tklabel(optLRFrame,text=listOptLR[[(i*2-1)]]$name),
             listWidgetLR[[(i*2-1)]],
             tklabel(optLRFrame,text=listOptLR[[(i*2)]]$name),
             listWidgetLR[[(i*2)]]);
    }
    if (length(listWidgetLR) %% 2 >0)
      tkgrid(tklabel(optLRFrame,text=listOptLR[[length(listOptLR)]]$name),
             listWidgetLR[[length(listOptLR)]]);
    assign("listWidgetLR",listWidgetLR,mainGuiEnv)
    assign("listVarLR",listVarLR,mainGuiEnv)
    assign("curLeafRank",id,mainGuiEnv)
    
  }
  
  buildLRopt();
  tkgrid(optLRFrame)
  tkbind(lbLR,"<<ListboxSelect>>",buildLRopt);

  ##Command Line building
  
  cmdLine <- function(){
    datasetV <-  tclvalue(dataset);
    resnameV <- tclvalue(resname);
    bestresponseV <- tclvalue(bestresponse);
    idLR <- get("curLeafRank",mainGuiEnv);
    #LeafRank options
    listOptLR <- optionLR[[idLR]]
    optPar <- list();
    listVarLR <- get("listVarLR",mainGuiEnv)
    for(i in 1:length(listVarLR)){
      
      if (listOptLR[[i]]$type=="listbox"){
        id <- as.numeric(tkcurselection(listVarLR[[i]]))+1;
        tmp<- paste("c(optPar,",listOptLR[[i]]["optName"],"=\"",
                    listOptLR[[i]]$choice[id],"\")",sep="");
      }else{
        tmp <- paste("c(optPar,",listOptLR[[i]]["optName"],"=",tclvalue(listVarLR[[i]]),")",sep="");
      }
      optPar <- eval(parse(text=tmp));
      
    }
    optLR <- opt2cmdLR[[idLR]](optPar);
                                        #TreeRank options
    minOptTRV <- as.double(tclvalue(minOptTR));
    maxOptTRV <- as.double(tclvalue(maxOptTR));
    mincritOptTRV <-as.double(tclvalue(mincritOptTR));
    nfcvOptTRV <- as.integer(tclvalue(nfcvOptTR));
    varSplitTRV <- as.numeric(tclvalue(varSplitTR));
    dataSplitTRV <- as.numeric(tclvalue(dataSplitTR));
    forestTRV <- as.numeric(tclvalue(forestOptTR));
    if (is.null(datasetV)||is.null(resnameV)||is.null(bestresponseV)
        ||is.null(minOptTRV)||is.null(maxOptTRV)||is.null(mincritOptTRV)
        ||is.null(nfcvOptTRV)
        ||is.null(varSplitTRV))
      {
        errorGui("undefined value");
        return
      }

    grTRctrl <- paste("growing_ctrl(maxdepth=",maxOptTRV,",mincrit=",mincritOptTRV,",minsplit=",minOptTRV,")",sep="");
                                        #    prTRctrl <- paste("pruning_ctrl(prune=",pruneOptTRV,",strat=",stratOptTRV,",nfcv=",nfcvOptTRV,")",sep="");
 #   prTRctrl <- paste("pruning_ctrl(prune=TRUE,strat=FALSE,nfcv=",nfcvOptTRV,")",sep="");
    listx <- 0;
    listy <- 0;
    pcinit <-0;
    ncinit <-0;
                                        #cb <-function(id,depth,pc,nc,pck,nck){
                                        # if (length(listx)==0){
                                        #  pcinit <- pc;
                                        #  ncinit <- nc;
                                        #}
    
                                        #}
    
                                        # print(paste(id,depth,pc,nc,pck,nck,sep=" "));
                                        #}
    if  (forestTRV>1){
    cmd <-  paste("TreeRankForest(formula=",resnameV,"~., data=",datasetV,",bestresponse=\"",bestresponseV,"\",growing=",grTRctrl,",nfcv=",nfcvOptTRV,",LeafRank= function(...){",funNameLR[[idLR]],"(",optLR,",...)},prcSplitVar=",(varSplitTRV/100),",prcSplitData=",(dataSplitTRV/100),",nbforest=",forestTRV,")",sep="");
  }
    else     cmd <-  paste("TreeRank(formula=",resnameV,"~., data=",datasetV,",bestresponse=\"",bestresponseV,"\",growing=",grTRctrl,",nfcv=",nfcvOptTRV,",LeafRank= function(...){",funNameLR[[idLR]],"(",optLR,",...)},prcSplitVar=",(varSplitTRV/100),",prcSplitData=",(dataSplitTRV/100),")",sep="");
    return(cmd);
    
  }
  
  runGui <- function(){
    tree <- eval(parse(text=cmdLine()),.GlobalEnv);
    assign("tree",tree,.GlobalEnv);
    TRplot(tree);
    print("Computation done");
  }
  
  exportCmd <- function(){
    cmd <- cmdLine();
    tt <- tktoplevel()
    name <- tclVar("funName");
    entry.Name <- tkentry(tt,width="20",textvariable=name)
    tkgrid(tklabel(tt,text="Fun name"));
    OnOk <- function()
      {
        NameVal <- tclvalue(name);
        tkdestroy(tt)
        assign(NameVal,cmd,.GlobalEnv);
      }
    OK.but <- tkbutton(tt,text=" set ",command=OnOk);
    tkbind(entry.Name,"<Return>",OnOk)
    tkgrid(entry.Name)
    tkgrid(OK.but)
    tkfocus(tt);
  }

  but.run <-tkbutton(mainFrame,text="run",command =runGui);
  #but.exp <- tkbutton(mainFrame,text="export cmd",command=exportCmd);
  tkgrid(infoFrame);
  
  tkgrid(LRFrame);
  tkgrid(TRFrame);
  tkgrid(but.run);
  tkbind(mainFrame,"<<Return>>",runGui);
                                        #tkgrid(but.exp);
  tkgrid(mainFrame);
  tkfocus(top);
  invisible(NULL);
}












#############################################################

TRplot <- function(tree,top= NULL){
  if (is.null(top))
    top <- tktoplevel(width="600",height="600");
  treeType <- -1;
  tktitle(top)<- "TreeRank GUI"
  if (inherits(tree,"TR_LRCart")){
    treeType <- 1;
    
    popup.menu<-
      list(list(label="View Unpruned",func=.TRgui.viewUnpruned),
           list(label="Save Tree",func = .TRgui.saveTree),
           list(label="Export Tree",func = .TRgui.exportTree));
    node.popup.menu <-
      list(list(label="View Unpruned",func=.TRgui.viewUnpruned),
           list(label="Save Tree",func = .TRgui.saveTree),
           list(label="Export Tree",func = .TRgui.exportTree));
  }
  if ((inherits(tree,"TR_TreeRank"))||(inherits(tree,"TR_forest"))){
    treeType <- 0;
    if (inherits(tree,"TR_forest"))
      treeType<-2;
    popup.menu<-
      list(list(label="View Subtree", func = .TRgui.viewSubTree),
           list(label="View Unpruned",func=.TRgui.viewUnpruned),
           list(label="Plot Subtree ROC", func = .TRgui.plotROCSubTree),
           list(label="Plot Unpruned ROC",func=.TRgui.plotROCUnpruned),
           list(label="Add test set", func = .TRgui.submitTestSet),
          list(label="Add ROC",func=.TRgui.addExtROC),
           list(label="Save Tree",func = .TRgui.saveTree),
           list(label="Export ROC",func = .TRgui.exportROC),
           list(label="Export Tree",func = .TRgui.exportTree));

    node.popup.menu <-
      list(list(label = "View LeafRank", func = .TRgui.viewLeafRank),
           list(label="View Subtree", func = .TRgui.viewSubTree),
           list(label="View Unpruned",func=.TRgui.viewUnpruned),
           list(label="Plot Subtree ROC", func = .TRgui.plotROCSubTree),
           list(label="Plot Unpruned ROC",func=.TRgui.plotROCUnpruned),
           list(label="Add test set", func = .TRgui.submitTestSet),
           list(label="Add ROC",func=.TRgui.addExtROC),
           list(label="Save Tree",func = .TRgui.saveTree),
            list(label="Export ROC",func = .TRgui.exportROC),
           list(label="Export Tree",func = .TRgui.exportTree));
  }
  if (treeType == -1){
    print("Bad class for tree");
    stop();
  }
  trg.id <- .TRgui.new(list(top=top,nodeSize = 10,panX=50,panY=50,maxNodeSize=30,minNodeSize=5,tree=tree,treeType = treeType));
    if (treeType==2){
      .TRgui.set(trg.id,"forest",tree);
      .TRgui.set(trg.id,"nbforest",tree$nbforest);
      .TRgui.set(trg.id,"tree",tree$forest[[1]]);
      .TRgui.set(trg.id,"forestCur",1);
    }

  
  treeFrame <- .TRgui.treeFrame(trg.id,popup.menu,node.popup.menu);
  infoFrame <- .TRgui.infoFrame(trg.id);
  if(treeType != 1){
    rightFrame <- .TRgui.rightFrame(trg.id);
  }
  tkgrid(infoFrame,row=0,column=0,sticky="nesw");
  tkgrid(treeFrame,row=0,column=1,sticky="nesw");
  if(treeType!=1){
    tkgrid(rightFrame,row=0,column=2,sticky="nesw");
  }
  tkgrid.columnconfigure(top,1);

  if (as.numeric(tclvalue(tkwinfo("viewable",top))!=1))
    tkwait.visibility(top);
  .TRgui.computeCoords(trg.id);
  .TRgui.scale(trg.id);
  .TRgui.create.nodes(trg.id);
  .TRgui.create.edges(trg.id);
  .TRgui.update.nodes(trg.id);
  .TRgui.refresh(trg.id);
  tkbind(top, "<Destroy>", function() TRgui.close(trg.id, FALSE))

  main.menu <- tkmenu(top)
  tkadd(main.menu, "command", label="Close", command=function() {
    TRgui.close(trg.id, TRUE)})
  sapply(node.popup.menu,function(x){tkadd(main.menu,"command",label=x$label,command=function(){x$f(trg.id)})})
  tkconfigure(top, "-menu", main.menu)
  
  trg.id;
 }


.TRgui.info2str <- function(trg.id,i = NULL){
  trg <- .TRgui.get(trg.id);
  tree <- .TRgui.get(trg.id,"tree")
  if (is.null(i)){
    ret <- "                                \nNode: \nDauc:\n#Pos:\n#Neg.: \n\n";
    return(ret);
  }
  else{
    ret <- paste("                                \nNode: ",i,
                 "\nDauc: ",format(tree$ldauc[[i]],digits=3),
                 "\n#Pos: ",tree$pcount[[i]],
                 "\n#Neg: ",tree$ncount[[i]],"\n",sep="");
    
  }
  if (trg$treeType != 1){
   if (tree$isleaf[i]){
  ret <- paste(ret,"Score: ",format(tree$score[[i]],digits=3),"\n",sep="");
   }else ret <- paste(ret,"\n");
   }
  ret;
}

.TRgui.infoFrame <- function(trg.id){
  top <- .TRgui.get(trg.id,"top");
  trg <- .TRgui.get(trg.id);
  tree <- .TRgui.get(trg.id,"tree");
  frame <- tkframe(top,borderwidth=2,width=120,relief="groove");
  .TRgui.set(trg.id,"infoFrame",frame);
  txt <- .TRgui.info2str(trg.id);
  infoLabel <-tklabel(frame,text= txt );
  tkgrid(infoLabel,sticky="w");
  .TRgui.set(trg.id,"infoLabel",infoLabel);
  lbVI <-  tklistbox(frame,height=10,yscrollcommand=function(...)tkset(scr,...));
  scr <-tkscrollbar(frame,command=function(...)tkyview(lbVI,...))
  if (.TRgui.get(trg.id,"treeType")==2){
    vi <- varImportance(.TRgui.get(trg.id,"forest"))
  }else{
  vi <- varImportance(tree);}
  if (!(is.null(vi))){
    tkgrid(tklabel(frame,text="Var. importance:"));
    for (i in 1:length(vi))
      tkinsert(lbVI,"end",paste(names(vi)[[i]],": ",format(vi[[i]],digits=4),sep=""));
    tkgrid(lbVI,scr);
    tkgrid.configure(scr,sticky="nsw")
  }
    .TRgui.set(trg.id,"lbVI",lbVI);
  if (trg$treeType !=1){
    frameListROC <- tkframe(frame,borderwidth=2,width=200,height=200,relief="groove")
        if (trg$treeType==2){
          
          foresttl<-tklistbox(frameListROC,height=5,selectmode="single",yscrollcommand=function(...)tkset(forestscr,...))
          forestscr <- tkscrollbar(frameListROC,
                                   command=function(...)tkyview(foresttl,...))
          tkgrid(tklabel(frameListROC,text="Forest"));
          for (i in (1:.TRgui.get(trg.id,"nbforest")))
            {
              tkinsert(foresttl,"end",i)
            }
          tkselection.set(foresttl,0);
          .TRgui.set(trg.id,"foresttl",foresttl);
        changeTree <- function(){
          tl <- .TRgui.get(trg.id,"foresttl");
          id <-as.numeric(tkcurselection(tl))+1;
          forest <- .TRgui.get(trg.id,"forest");
          .TRgui.set(trg.id,"tree",forest$forest[[id]]);
          .TRgui.set(trg.id,"forestCur",id);
          .TRgui.clean(trg.id);
          .TRgui.computeCoords(trg.id);
          .TRgui.scale(trg.id);
          .TRgui.create.nodes(trg.id);
          .TRgui.create.edges(trg.id);
          .TRgui.update.nodes(trg.id);
          .TRgui.refresh(trg.id);
          .TRgui.addROC(trg.id,getROC(forest$forest[[id]]),paste("Learn. ",id));
          .TRgui.updateROC(trg.id);
          .TRgui.viFrameReplot(trg.id);
          
        }
          tkgrid(foresttl,forestscr)
          tkgrid.configure(forestscr,sticky="nsw")
          tkbind(foresttl,"<<ListboxSelect>>",changeTree);
        }
    tkgrid(tklabel(frameListROC,text="ROC List"),tklabel(frameListROC,text="auc"));
    tkgrid(frameListROC,sticky="nsew")
    .TRgui.set(trg.id,"frameListROC",frameListROC);

    
    }
             #txt <- paste(txt,names(vi)[[i]],": ",vi[[i]],"\n",sep="");
  #tkgrid(tklabel(frame,text=txt),sticky="w");
  return(frame);
}

.TRgui.rightFrame <- function(trg.id){
  top <- .TRgui.get(trg.id,"top");
  frame <-tkframe(top,borderwidth=2,relief="groove");
  .TRgui.set(trg.id,"rightFrame",frame);
  frameROC <- .TRgui.frameROC(trg.id);
  viFrame <-.TRgui.viFrame(trg.id);
  frame;
}



.TRgui.viFrame <- function(trg.id){
  top<- .TRgui.get(trg.id,"rightFrame");
  vi <- varImportance(.TRgui.get(trg.id,"tree"));
  if (is.null(vi))return(tkframe(top));
  viFrame <- tkrplot(top,function() barplot(vi,col=heat_hcl(length(vi))[rank(-vi)]),vscale=0.7);
  tkgrid(viFrame,sticky="nesw");
  .TRgui.set(trg.id,"viFrame",viFrame);
  viFrame;
}

.TRgui.viFrameReplot <- function(trg.id){
  vi <- varImportance(.TRgui.get(trg.id,"tree"));
  viFrame <- .TRgui.get(trg.id,"viFrame");
  if (is.null(vi))return(tkframe(viFrame));
  tkrreplot(viFrame,function() barplot(vi,col=heat_hcl(length(vi))[rank(-vi)]),vscale=0.7);
  invisible(NULL);
}


.TRgui.treeFrame <- function(trg.id,popup.menu,node.popup.menu){
  top <- .TRgui.get(trg.id,"top");
  frame <- tkframe(top,relief="groove",borderwidth=2);
  canvas <- tkcanvas(frame,width=500,height=700);
  .TRgui.set(trg.id,"treeCanvas",canvas);
  .TRgui.set(trg.id,"treeFrame",frame);
  tmp<-tkframe(frame)
  yscr <- tkscrollbar(frame,repeatinterval=5, command = function(...) tkyview(canvas,...));
  tkconfigure(canvas, yscrollcommand= function(...) tkset(yscr,...));
  xscr <- tkscrollbar(frame, repeatinterval=5,orient = "horizontal",command = function(...) tkxview(canvas,...));
  tkconfigure(canvas, xscrollcommand= function(...) tkset(xscr,...));
  but.zoom <- tkbutton(tmp,text="+",command=function(...).TRgui.zoom(trg.id,1.2,...));
  but.dezoom <- tkbutton(tmp,text="-",command=function(...) .TRgui.zoom(trg.id,1/1.2,...));
  tkgrid(canvas,yscr,stick="ns");
  tkgrid(xscr,sticky="ew");
  tkgrid(but.dezoom,but.zoom)
  tkgrid(tmp);

  #zoom
  tkbind(canvas,"<4>",function().TRgui.zoom(trg.id,1.2,TRUE));
  tkbind(canvas,"<5>",function().TRgui.zoom(trg.id,1/1.2,TRUE));
  tkbind(canvas,"<Configure>",function() .TRgui.reinitGUI(trg.id));
  #node (de/)selection
  tkitembind(canvas,"node||nodelabel","<1>",function(x,y){
    trg <- .TRgui.get(trg.id);
    .TRgui.deselect.all(trg.id)
    .TRgui.select.current(trg.id)
    
  })
  
  tkitembind(canvas,"node||nodelabel","<Control-1>",function(x,y){
    treeCanvas <-.TRgui.get(trg.id,"treeCanvas");
    curtags <- as.character(tkgettags(canvas,"current"));
    seltags <- as.character(tkgettags(canvas,"selected"));
    .TRgui.select.current(trg.id);
  });


  pm <- tkmenu(canvas);
  node.pm <- tkmenu(canvas);
  sapply(popup.menu,function(x){tkadd(pm,"command",label = x$label, command = function(){x$f(trg.id)})})
  sapply(node.popup.menu,function(x){tkadd(node.pm,"command",label = x$label, command = function(){x$f(trg.id)})})
  
  tkbind(canvas,"<3>",function(x,y){
    treeCanvas <- .TRgui.get(trg.id,"treeCanvas");
    tags <- as.character(tkgettags(treeCanvas,"current"));
    if (!("selected" %in% tags)){
      .TRgui.deselect.all(trg.id)
      .TRgui.select.current(trg.id);
    }
    tags <- as.character(tkgettags(treeCanvas,"selected"));
    if ("node" %in% tags){
      menu <- node.pm;
    }else{menu <- pm};
        
    
    x <- as.integer(x) + as.integer(tkwinfo("rootx", canvas))
    y <- as.integer(y) + as.integer(tkwinfo("rooty", canvas))
    .Tcl(paste("tk_popup", .Tcl.args(menu, x, y)))

  })
  return(frame);
}


.TRgui.zoom<- function(trg.id,zoom,m=FALSE){
  trg <- .TRgui.get(trg.id);
  x <-0;
  y<-0;
  if (m){
    x <- as.numeric(tclvalue(tkwinfo("pointerx",trg$treeCanvas)))-as.numeric(tclvalue(tkwinfo("rootx",trg$treeCanvas)));
    y <- as.numeric(tclvalue(tkwinfo("pointery",trg$treeCanvas)))-as.numeric(tclvalue(tkwinfo("rooty",trg$treeCanvas)));
  }
   
  tcl(trg$treeCanvas,"scale","all",x,y,zoom,zoom);
.TRgui.refresh(trg.id);
  invisible(NULL);
  
}
    
.TRgui.frameROC<- function(trg.id){
  top <- .TRgui.get(trg.id,"rightFrame");
  tree <- .TRgui.get(trg.id,"tree");
  frameROC <-tkrplot(top,function() 0);
    tabc = rainbow(10);
  .TRgui.set(trg.id,"frameROC",frameROC);
  .TRgui.set(trg.id,"rocList",list());
  .TRgui.set(trg.id,"aucList",list());
  .TRgui.set(trg.id,"rocListColor",tabc);
  .TRgui.set(trg.id,"rocColor",list());
  .TRgui.set(trg.id,"indexColorROC",1);
  .TRgui.set(trg.id,"listVB",list());
  .TRgui.set(trg.id,"listCB",list());
  if (.TRgui.get(trg.id,"treeType")==2){
    .TRgui.addROC(trg.id,getROC(.TRgui.get(trg.id,"forest")),"Learn. forest");
   # .TRgui.addROC(trg.id,getROC(tree),"Learn. 1")
  }
  .TRgui.addROC(trg.id,getROC(tree),"Learn. 1")
  tkgrid(frameROC,sticky="nsew");
  .TRgui.updateROC(trg.id);
  frameROC;
}



.TRgui.addROC <- function(trg.id,roc,txt)
{
  trg <- .TRgui.get(trg.id)
  .TRgui.set(trg.id,"rocList",c(trg$rocList,list(roc)));
  .TRgui.set(trg.id,"aucList",c(trg$aucList,auc(roc)));
  .TRgui.set(trg.id,"rocColor",c(trg$rocColor,trg$rocListColor[[trg$indexColorROC]]));
  .TRgui.set(trg.id,"indexColorROC",(((trg$indexColorROC) %% length(trg$rocListColor))+1));
  frameListROC <- .TRgui.get(trg.id,"frameListROC");
  vb <- tclVar("1");
  cb <- tkcheckbutton(frameListROC);
  tkconfigure(cb,variable=vb,command=function(...){.TRgui.updateROC(trg.id)});
  .TRgui.set(trg.id,"listVB",c(trg$listVB,list(vb)));
  .TRgui.set(trg.id,"listCB",c(trg$listCB,list(cb)));
  tkgrid(tklabel(frameListROC,text=txt),tklabel(frameListROC,text=format(auc(roc),digits=3)),cb);
  
  .TRgui.updateROC(trg.id);
  invisible(NULL);
}


.TRgui.updateROC <- function(trg.id){
  frame <- .TRgui.get(trg.id,"frameROC");
  rocList <- .TRgui.get(trg.id,"rocList")
  rocColor <- .TRgui.get(trg.id,"rocColor")
  #rocStyle <- .TRgui.get(trg.id,"rocStyle")
  listVB <- as.numeric(sapply(.TRgui.get(trg.id,"listVB"),function(x) tclvalue(x)));
  tree <- .TRgui.get(trg.id,"tree")
  vids <- .TRgui.get.selected.nodes(trg.id)
  points <- list();
  for (i in vids)
    {
      if (tree$isleaf[i]){
        points <- c(points,list(c(tree$lalpha[i],tree$lbeta[i])))
      }
      else{points <- c(points,list(c(tree$lalpha[tree$kidslist[[i]][2]],
                                     tree$lbeta[tree$kidslist[[i]][[2]]])));
         }
    }
#  plotROC(rocList[listVB==1],rocColor[listVB==1],points)
  tkrreplot(frame,function()plotROC(rocList[listVB==1],rocColor[listVB==1],points))
  invisible(NULL);
}

  


.TRgui.viewLeafRank <- function(id){
  tree <- .TRgui.get(id,"tree")
  vids <- .TRgui.get.selected.nodes(id)
  for (i in vids)
    TRplot(getClassifier(tree,i));
  invisible(NULL)
}


.TRgui.viewSubTree <- function(id){
  tree <- .TRgui.get(id,"tree")
  vids <- .TRgui.get.selected.nodes(id)
  if (length(vids)==0) return(FALSE)
  TRplot(subTreeRank(tree,vids))
  invisible(NULL)
}


.TRgui.viewUnpruned <- function(id){
  tree <- .TRgui.get(id,"tree")
  if (is.null(tree$unpruned)){
    tkmessageBox(title= "TreeRank",message="Not pruned tree",type="ok")
  
  }else
  {TRplot(tree$unpruned)}
  invisible(NULL)
}

.TRgui.plotROCUnpruned <- function(id){
  tree <- .TRgui.get(id,"tree");
 
  if (is.null(tree$unpruned)){
    tkmessageBox(title="TreeRank",message="Not pruned tree",type="ok")
  }else{
  .TRgui.addROC(id,getROC(tree$unpruned),"Unpruned Tree");
  .TRgui.updateROC(id);}
  invisible(NULL)
}


.TRgui.plotROCSubTree <- function(id){
  trg <- .TRgui.get(id);
  vids <- .TRgui.get.selected.nodes(id)
  if (length(vids)==0) return(FALSE)
  newtree <-subTreeRank(trg$tree,vids)
  .TRgui.addROC(id,getROC(newtree),"Subtree");
  .TRgui.updateROC(id);
  invisible(NULL)
} 

.TRgui.submitTestSet <- function(id){
  tree <- .TRgui.get(id,"tree")	
  tt<-tktoplevel()
  Name <- tclVar("testset")
  entry.Name <-tkentry(tt,width="20",textvariable=Name)
  tkgrid(tklabel(tt,text="Test Set"))
  tkgrid(entry.Name)
  OnOK <- function()
    {
      NameVal <- tclvalue(Name)
      tkdestroy(tt)
      if ((.TRgui.get(id,"treeType")==2)){
  #      if (as.numeric(tclvalue(vb))!=0){
        forest <- .TRgui.get(id,"forest");
        roc <- eval(parse(text=paste("getROC(forest,",NameVal,")",sep="")));
        .TRgui.addROC(id,roc,paste(NameVal,"forest"));
      }
      roc <- eval(parse(text=paste("getROC(tree,",NameVal,")",sep="")))
      if ((.TRgui.get(id,"treeType")==2))
        .TRgui.addROC(id,roc,paste(NameVal,.TRgui.get(id,"forestCur")))
      else      .TRgui.addROC(id,roc,NameVal);
      .TRgui.updateROC(id);
    }
  OK.but <-tkbutton(tt,text="  set  ",command=OnOK)
  tkbind(entry.Name, "<Return>",OnOK)
  tkgrid(OK.but)
  tkfocus(tt)
  invisible(NULL)
}

.TRgui.addExtROC <- function(id){
  tt <- tktoplevel()
  Name <-tclVar("roc")
  entry.Name <- tkentry(tt,width="20",textvariable=Name)
  tkgrid(tklabel(tt,text="ROC var. name"))
  tkgrid(entry.Name)
  OnOK <- function()
    {
      NameVal <- tclvalue(Name)
      tkdestroy(tt)
      roc <- eval(parse(text=NameVal));
      trg <- .TRgui.get(id);
      .TRgui.addROC(id,as.matrix(roc),NameVal)
      .TRgui.updateROC(id);
    }
  OK.but <- tkbutton(tt,text="  set  ",command=OnOK)
  tkbind(entry.Name,"<Return>",OnOK)
  tkgrid(OK.but)
  tkfocus(tt)
  invisible(NULL)
}


.TRgui.exportTree <- function(trg.id){
  tt <- tktoplevel()
  frame <- .TRgui.get(trg.id,"treeCanvas");
  name <- tclVar("file.eps")
  entry.Name <- tkentry(tt,width="20",textvariable=name)
  tkgrid(tklabel(tt,text="File Name"))
  tkgrid(entry.Name);
  OnOk <- function(){
    NameVal <- tclvalue(name);
    tkdestroy(tt);
    tkpostscript(frame,file=NameVal);
    tkmessageBox(title= "Info",message="Tree saved",type="ok")
  }
  OK.but <- tkbutton(tt,text=" save ",command=OnOk);
  tkbind(entry.Name,"<Return>",OnOk);
  tkgrid(OK.but);
  tkfocus(tt)
  invisible(NULL);
}


 .TRgui.exportROC <- function(trg.id){
  tt <- tktoplevel()
  rocList <- .TRgui.get(trg.id,"rocList")
  rocColor <- .TRgui.get(trg.id,"rocColor");
  listVB <-  as.numeric(sapply(.TRgui.get(trg.id,"listVB"),function(x) tclvalue(x)));
  name <- tclVar("file.eps")
  entry.Name <- tkentry(tt,width="20",textvariable=name)
  tkgrid(tklabel(tt,text="File Name"))
  tkgrid(entry.Name);
  OnOk <- function(){
    NameVal <- tclvalue(name);
    tkdestroy(tt);
    postscript(NameVal)
    plotROC(rocList[listVB==1],rocColor[listVB==1]);
    dev.off()
    tkmessageBox(title= "Info",message="ROC saved",type="ok")

  }
  OK.but <- tkbutton(tt,text=" save ",command=OnOk);
  tkbind(entry.Name,"<Return>",OnOk);
  tkgrid(OK.but);
  tkfocus(tt)
  invisible(NULL);
}
.TRgui.saveTree <- function(id){
  tr <- .TRgui.get(id,"tree")
  tt <- tktoplevel()
  name <- tclVar("varname")
  entry.Name <- tkentry(tt,width="20",textvariable=name)
  tkgrid(tklabel(tt,text="Var name"))
  tkgrid(entry.Name);
  OnOk <- function()
    {
      NameVal <- tclvalue(name)
      tkdestroy(tt)
      assign(NameVal,tr,.GlobalEnv)
    }
  OK.but <- tkbutton(tt,text=" set ",command=OnOk)
  tkbind(entry.Name,"<Return>",OnOk)
  tkgrid(OK.but)
  tkfocus(tt)
  invisible(NULL)
}
  
.TRgui.refresh <-function(trg.id){
  trg <- .TRgui.get(trg.id);
  bbox <- tclvalue(tcl(trg$treeCanvas,"bbox","all"));
  bbox <- as.numeric(unlist(strsplit(bbox,split=" ")));
  if(length(bbox>1)){
  tkconfigure(trg$treeCanvas, scrollregion =c(bbox[[1]],bbox[[2]],bbox[[3]],bbox[[4]]));
  .TRgui.set(trg.id,"treeScrollRegion",bbox);}
  invisible(NULL)
}


.TRgui.create.node <- function(trg.id, id, label=NULL, x=0, y=0) {
  trg <- .TRgui.get(trg.id)
  node.size <-trg$sizeNodeList[[id]];
  node.color <- trg$colorNode[[id]];
  node.out.color <- trg$outColorNode[[id]];
  item <- tkcreate(trg$treeCanvas, "oval", x-node.size, y-node.size,
                   x+node.size, y+node.size, width=1,
                   outline=node.out.color,  fill=node.color)
  tkaddtag(trg$treeCanvas, "node", "withtag", item)
  tkaddtag(trg$treeCanvas, paste("n-", id, sep=""), "withtag", item)
  if (!(is.null(label))){
    litem <- tkcreate(trg$treeCanvas, "text", x,y,text=label)
    tkaddtag(trg$treeCanvas, "nodelabel", "withtag", litem)
    tkaddtag(trg$treeCanvas, paste("nlabel-", id, sep=""), "withtag", litem)
  }
    item;
}


.TRgui.create.edge <- function(trg.id, from, to, id,label=NULL) {
  trg <- .TRgui.get(trg.id);
  from.c <- trg$treeCoords[from,]
  to.c   <- trg$treeCoords[to,]
  edge.color <- "black";
  edge.width <-"2";
  phi <- atan2(to.c[2]-from.c[2], to.c[1]-from.c[1])
  r <- sqrt( (to.c[1]-from.c[1])^2 + (to.c[2]-from.c[2])^2 )
  to.c[1] <- from.c[1] + (r-trg$sizeNodeList[[to]])*cos(phi)
    to.c[2] <- from.c[2] + (r-trg$sizeNodeList[[to]])*sin(phi)
  from.c[1] <- from.c[1] +trg$sizeNodeList[[from]]*cos(phi)
  from.c[2] <- from.c[2] + trg$sizeNodeList[[from]]*sin(phi)

  coords <- c(from.c[1], from.c[2], to.c[1], to.c[2])
  args <- c(list(trg$treeCanvas, "line"),
            coords, 
            list(fill=edge.color, activefill="red", 
                   tags=c("edge", paste(sep="", "edge-", id),
                     paste(sep="", "from-", from),
                     paste(sep="", "to-", to))))
  do.call(tkcreate, args)
  if (!(is.null(label))){
    label.x <- (to.c[1]+from.c[1])/2
    label.y <- (to.c[2]+from.c[2])/2
    litem <- tkcreate(trg$treeCanvas,"text",label.x,label.y,
                      text=as.character(label));
    tkaddtag(trg$treeCanvas,"label","withtag",litem);
    tkaddtag(trg$treeCanvas,paste(sep="","elabel-",id),"withtag",litem)
  }
  invisible(NULL)
  
}

.TRgui.create.nodes <- function(trg.id) {
  trg <- .TRgui.get(trg.id);
  tree <- trg$tree;
  treeType <- .TRgui.get(trg.id,"treeType");
  nodeLab <- NULL

  if (treeType == 1){
    nodeLab <-  unlist(lapply(tree$nodes,function(x) {
    
    if (!(sum(is.na(tree$split[[x]])))){
        return(tree$split[[x]]$name)
    }
    else return("");
  }))
  }
  .TRgui.set(trg.id,"nodeLab",nodeLab);
  if (!(is.null(nodeLab))){
    mapply(function(n, l, x, y) .TRgui.create.node(trg.id, n, l, x, y),
           tree$nodes, nodeLab, trg$treeCoords[,1], trg$treeCoords[,2])
  }
  else{
    mapply(function(n,  x, y) .TRgui.create.node(trg.id, n, NULL, x, y),
           tree$nodes,  trg$treeCoords[,1], trg$treeCoords[,2])
 }
  invisible(NULL)
}


.TRgui.create.edges <- function(trg.id) {
  trg <- .TRgui.get(trg.id)
  tree <- trg$tree;
  edgematrix <- .TRgui.get(trg.id,"edgesList");

  edgeLab <- NULL;
  if (.TRgui.get(trg.id,"treeType") == 1){
    edgeLab <- unlist(lapply(tree$inner,function(x)
                                {sp <- tree$split[[x]];
                                 if (!(is.null(sp))){
                                   lab <-c(paste("<",format(sp$breaks,digits=3)),paste(">=",format(sp$breaks,digits=3)));
                                   return(strtrim(lab,8));
                                 }
                                 else{return(c("",""))}}))
  }
                               
  .TRgui.set(trg.id,"edgeLab",edgeLab);
  if (!(is.null(edgeLab))){
    mapply(function(from, to, id,l) .TRgui.create.edge(trg.id, from, to, id,l),
           edgematrix[,1],
           edgematrix[,2], 1:nrow(edgematrix),edgeLab)
  }
  else{
   mapply(function(from, to, id) .TRgui.create.edge(trg.id, from, to, id,NULL),
           edgematrix[,1],
           edgematrix[,2], 1:nrow(edgematrix))
 }
  invisible(NULL)
}





.TRgui.update.node <- function(trg.id, id, x, y) {
  trg <- .TRgui.get(trg.id)
  node.size <- trg$sizeNodeList[[id]];
  node.color <- trg$colorNode[[id]];
  node.out.color <- trg$outColorNode[[id]];
 
  # Vertex
  tkcoords(trg$treeCanvas, paste("node&&n-", id, sep=""),
           x-node.size, y-node.size,
           x+node.size, y+node.size)
  # Label
  #  .TRgui.update.label(trg.id, id, x, y)
  
  # Edges
  edge.from.ids <- as.numeric(tkfind(trg$treeCanvas, "withtag",
                                     paste("from-", id, sep="")))
  edge.to.ids <- as.numeric(tkfind(trg$treeCanvas, "withtag",
                                   paste("to-", id, sep="")))
  for (i in seq(along=edge.from.ids)) {
    .TRgui.update.edge(trg.id, edge.from.ids[i])
  }
  invisible(NULL)
}

.TRgui.update.nodes <- function(trg.id) {
  trg <- .TRgui.get(trg.id)
  n <- trg$tree$nbNode;
  mapply(function(v, x, y) .TRgui.update.node(trg.id, v, x, y), 1:n,
         trg$treeCoords[,1], trg$treeCoords[,2])
  invisible(NULL);
}
.TRgui.update.edges <- function(trg.id) {
  trg <- .TRgui.get(trg.id)
  n <- length(length(trg$tree$inner)*2);
  mapply(function(v) .TRgui.update.edgeById(trg.id, v), 1:n);
  invisible(NULL);
}

# Update an edge with given itemid (not edge id!)
.TRgui.update.edge <- function(trg.id, itemid) {
  trg <- .TRgui.get(trg.id)
  tags <- as.character(tkgettags(trg$treeCanvas, itemid))
  from <- as.numeric(substring(grep("from-", tags, value=TRUE, fixed=TRUE),6))
  to <- as.numeric(substring(grep("to-", tags, value=TRUE, fixed=TRUE),4))
  from.c <- trg$treeCoords[from,]
  to.c <- trg$treeCoords[to,]
    phi <- atan2(to.c[2]-from.c[2], to.c[1]-from.c[1])
  r <- sqrt( (to.c[1]-from.c[1])^2 + (to.c[2]-from.c[2])^2 )
  to.c[1] <- from.c[1] + (r-trg$sizeNodeList[[to]])*cos(phi)
    to.c[2] <- from.c[2] + (r-trg$sizeNodeList[[to]])*sin(phi)
  from.c[1] <- from.c[1] +trg$sizeNodeList[[from]]*cos(phi)
  from.c[2] <- from.c[2] + trg$sizeNodeList[[from]]*sin(phi)

  tkcoords(trg$treeCanvas, itemid, from.c[1], from.c[2], to.c[1], to.c[2]);
  invisible(NULL);
}  


.TRgui.update.edgeById <- function(trg.id, id) {
  
  trg <- .TRgui.get(trg.id)
  itemid <- as.numeric(tkfind(trg$treeCanvas,"withtag",paste("edge-",id,sep="")));
  .TRgui.update.edge(trg.id,itemid);
  invisible(NULL)
}  





###################################################################
# Internal functions, handling data about layout
###################################################################

.TRgui.scale <- function(trg.id){
  trg <- .TRgui.get(trg.id);
  width <- as.numeric(tkwinfo("width",trg$treeCanvas));
  height <- as.numeric(tkwinfo("height",trg$treeCanvas));
  treeCoords <- trg$treeCoords;
  treeCoords[,1] <- trg$treeCoordsNorm[,1]*(width-2*trg$panX)+trg$panX;
  treeCoords[,2] <- trg$treeCoordsNorm[,2]*(height-2*trg$panY)+trg$panY;
  .TRgui.set(trg.id,"width",width);
  .TRgui.set(trg.id,"height",height);
  .TRgui.set(trg.id,"treeCoords",treeCoords);
  invisible(NULL)
  
}
  
.TRgui.computeCoords <- function(trg.id){
  trg <- .TRgui.get(trg.id);
  tree <- trg$tree;
  nodeSize <-trg$nodeSize;
  listedge <- lapply(tree$inner,function(x){
    return(c(x,tree$kidslist[[x]][[1]],x,tree$kidslist[[x]][[2]]))
  });
  edgematrix <- matrix(unlist(listedge),ncol=2,byrow=TRUE);
  treeCoords <- layout.norm(layout.reingold.tilford(graph.edgelist(edgematrix-1),root=0,mode="all"),xmin=0,xmax=1,ymin=0,ymax=1);

  propPos <- tree$pcount/(tree$pcount+tree$ncount);
  propNb <- (tree$pcount+tree$ncount)/(tree$pcount[[1]]+tree$ncount[[1]]);
  scaleM <- log(2)/log(1+1/min(propNb));
  propNb <- (trg$maxNodeSize-trg$minNodeSize)/(1-scaleM)*(log(2)/log(1+1/propNb)-1)+trg$maxNodeSize;
  pal <- diverge_hcl(100,h=c(0,120),c=100,l=c(70,100),power=1,gamma=1);

  propPos <- pal[ceiling(propPos*99)+1];
  outColor <- rep("black",length(propPos))
  if (trg$treeType == 1){
    outColor[trg$tree$Lnode] <- "green"
    outColor[trg$tree$Rnode] <- "red"
  }
  
  .TRgui.set(trg.id,"sizeNodeList",propNb);
  .TRgui.set(trg.id,"colorNode",propPos);
  .TRgui.set(trg.id,"outColorNode",outColor);
  .TRgui.set(trg.id,"treeCoords",treeCoords);
  .TRgui.set(trg.id,"treeCoordsNorm",treeCoords);
  .TRgui.set(trg.id,"edgesList",edgematrix)
  invisible(NULL)
}  


.TRgui.reinitGUI <- function(trg.id){
  .TRgui.computeCoords(trg.id);
  .TRgui.scale(trg.id);
  .TRgui.update.nodes(trg.id);
  .TRgui.refresh(trg.id);
  invisible(NULL)
}


###################################################################
# Internal functions, handling the internal environment
###################################################################

.TRgui.new <- function(trg) {
  id <- get(".next", .TRgui.env)
  assign(".next", id+1, .TRgui.env)
  assign("tmp", trg, .TRgui.env)
  cmd <- paste("trg.", id, "<- tmp", sep="")
  eval(parse(text=cmd), .TRgui.env)
  rm("tmp", envir=.TRgui.env)
  id
}

.TRgui.get <- function(trg.id, what=NULL) {
  if (is.null(what)) {
    get(paste("trg.", trg.id, sep=""), .TRgui.env)
  } else {
    cmd <- paste("trg.", trg.id, "$", what, sep="")
    eval(parse(text=cmd), .TRgui.env)
  }
}

.TRgui.set <- function(trg.id, what, value) {
  assign("tmp", value, .TRgui.env)
  cmd <- paste(sep="", "trg.", trg.id, "$", what, "<-tmp")
  eval(parse(text=cmd), .TRgui.env)
  rm("tmp", envir=.TRgui.env)
  invisible(NULL)
}


.TRgui.deselect.all <- function(trg.id) {
  treeCanvas <- .TRgui.get(trg.id, "treeCanvas")
  ids <- as.numeric(tkfind(treeCanvas, "withtag", "selected"))
  for (i in ids) {
    .TRgui.deselect.this(trg.id, i)
  }
  invisible(NULL)
}

.TRgui.select.all.nodes <- function(trg.id) {
  treeCanvas <- .TRgui.get(trg.id, "treeCanvas")
  nodes <- as.numeric(tkfind(treeCanvas, "withtag", "node"))
  for (i in nodes) {
    .TRgui.select.node(trg.id, i)
  }
  invisible(NULL)
}

.TRgui.select.some.nodes <- function(trg.id, vids) {
  treeCanvas <- .TRgui.get(trg.id, "treeCanvas")
  vids <- unique(vids)
  for (i in vids) {
    tkid <- as.numeric(tkfind(treeCanvas, "withtag",
                              paste(sep="", "node&&n-", i)))
    .TRgui.select.node(trg.id, tkid)
  }
  invisible(NULL)
}


.TRgui.select.node <- function(trg.id, tkid) {
  treeCanvas <- .TRgui.get(trg.id, "treeCanvas")
  tkaddtag(treeCanvas, "selected", "withtag", tkid)
  tkitemconfigure(treeCanvas, tkid, "-outline", "blue",
                  "-width", 2);
  
  .TRgui.updateROC(trg.id);
  invisible(NULL)
}

.TRgui.deselect.node <- function(trg.id, tkid) {
  treeCanvas <- .TRgui.get(trg.id, "treeCanvas")
  outColorNode <- .TRgui.get(trg.id,"outColorNode");
  tkdtag(treeCanvas, tkid, "selected")
  trg <- .TRgui.get(trg.id)
  tags <- as.character(tkgettags(treeCanvas, tkid))
  id <- as.numeric(substring(tags[pmatch("n-", tags)], 3))
  tkitemconfigure(treeCanvas, tkid, "-outline",outColorNode[[id]],
                  "-width", 1)
  .TRgui.updateROC(trg.id);
  invisible(NULL)
}

.TRgui.select.current <- function(trg.id) {
  treeCanvas <- .TRgui.get(trg.id, "treeCanvas")
  tkid <- as.numeric(tkfind(treeCanvas, "withtag", "current"))
  .TRgui.select.this(trg.id, tkid)
  invisible(NULL)
}

.TRgui.deselect.current <- function(trg.id) {
  treeCanvas <- .TRgui.get(trg.id, "treeCanvas")
  tkid <- as.numeric(tkfind(treeCanvas, "withtag", "current"))
  .TRgui.deselect.this(trg.id, tkid)
  invisible(NULL)
}

.TRgui.select.this <- function(trg.id, tkid) {
  treeCanvas <- .TRgui.get(trg.id, "treeCanvas")
  trg <- .TRgui.get(trg.id);
  tags <- as.character(tkgettags(treeCanvas, tkid))
  id <- 1
  if ("node" %in% tags) {
    id<- as.numeric(substring(tags[pmatch("n-",tags)],3));
    .TRgui.select.node(trg.id, tkid)
  }else if ("nodelabel" %in% tags){
    id <- as.numeric(substring(tags[pmatch("nlabel-",tags)],8));
    tkid <- as.character(tkfind(treeCanvas,"withtag",paste(sep="","node&&n-",id)))
    .TRgui.select.node(trg.id,tkid)
  }
  txt <- .TRgui.info2str(trg.id,id);
  tkconfigure(trg$infoLabel,text=txt);
  invisible(NULL)
}

.TRgui.deselect.this <- function(trg.id, tkid) {
  treeCanvas <- .TRgui.get(trg.id, "treeCanvas")
  tags <- as.character(tkgettags(treeCanvas, tkid))
  if ("node" %in% tags) {
    .TRgui.deselect.node(trg.id, tkid)
  }
  invisible(NULL)
}
.TRgui.get.selected.nodes <- function(trg.id) {
  treeCanvas <- .TRgui.get(trg.id, "treeCanvas")
  tkids <- as.numeric(tkfind(treeCanvas, "withtag", "node&&selected"))

  ids <- sapply(tkids, function(tkid) {
    tags <- as.character(tkgettags(treeCanvas, tkid))
    id <- as.numeric(substring(tags [pmatch("n-", tags)], 3))
    id})

  ids
}

.TRgui.clean <- function(trg.id){
  treeCanvas <- .TRgui.get(trg.id,"treeCanvas");
  tkdelete(treeCanvas,"all");
}


TRgui.close <- function(trg.id, window.close=TRUE) {
  if (window.close) {
    cmd <- paste(sep="", "trg.", trg.id, "$top")
    top <- eval(parse(text=cmd), .TRgui.env)
    tkbind(top, "<Destroy>", "")
    tkdestroy(top)
  }
  cmd <- paste(sep="", "trg.", trg.id)
#  rm(list=cmd, envir=.TRgui.env)
  invisible(NULL)
}




################################################################################
# Two Sample interface






TwoSampleGui <- function(){

  top <- tktoplevel(height="200",width="100")
  mainGuiEnv <- new.env();
  tktitle(top) <- "TwoSample GUI"
  mainFrame <- tkframe(top,height="200",width="100",borderwidth=2);

  infoFrame <- tkframe(mainFrame,borderwidth=5,width=100);
  LRFrame <- tkframe(mainFrame,borderwidth=5,width=100,relief="groove");
  TRFrame <-tkframe(mainFrame,borderwidth=5,width=100,relief="groove");
  
  optTRFrame <- tkframe(TRFrame,borderwidth=2,width=100);
  optLRFrame <- tkframe(LRFrame,borderwidth=2,width=100);

  #info frame
  dataset1 <- tclVar("");
  dataset2 <- tclVar("");
  alpha <- tclVar("5");
  split <- tclVar("60");
  
  entry.dataset1 <- tkentry(infoFrame,width="10",textvariable=dataset1);
  entry.dataset2 <- tkentry(infoFrame,width="10",textvariable=dataset2);
  entry.alpha <- tkentry(infoFrame,width="5",textvariable=alpha);
  entry.split <- tkentry(infoFrame,width="5", textvariable = split);
  tkgrid(tklabel(infoFrame,text="1st Data set: "),entry.dataset1);
  tkgrid(tklabel(infoFrame,text="2nd Data set: "),entry.dataset2);
  tkgrid(tklabel(infoFrame,text="Learning % size: "),entry.split);
  tkgrid(tklabel(infoFrame,text="Confidence level: "),entry.alpha);


  #TreeRank options frame
  dftOpt <- growing_ctrl();
  minOptTR <- tclVar(dftOpt$minsplit);
  maxOptTR <- tclVar(dftOpt$maxdepth);
  mincritOptTR <- tclVar(dftOpt$mincrit);
  varSplitTR <- tclVar(100);
  dataSplitTR <- tclVar(100);
  nfcvOptTR <-tclVar("0");
  entry.minOptTR <-tkentry(optTRFrame,width="5",textvariable=minOptTR);
  entry.maxOptTR <-tkentry(optTRFrame,width="5",textvariable=maxOptTR);
  entry.mincritOptTR <-tkentry(optTRFrame,width="5",textvariable=mincritOptTR);
  entry.varSplitTR <- tkentry(optTRFrame,width="5",textvariable=varSplitTR);
  entry.dataSplitTR <- tkentry(optTRFrame,width="5",textvariable=dataSplitTR);
  
  entry.nfcvOptTR <- tkentry(optTRFrame,width=5,textvariable=nfcvOptTR);
  

  tkgrid(tklabel(optTRFrame,text="Minimum Split"),entry.minOptTR,
         tklabel(optTRFrame,text="Min. Criteria"),entry.mincritOptTR);
  tkgrid(tklabel(optTRFrame,text="Maximum Depth"),entry.maxOptTR,
         tklabel(optTRFrame,text="%Var. split"),entry.varSplitTR);
  tkgrid(tklabel(optTRFrame,text="%Data split"),entry.dataSplitTR,
         tklabel(optTRFrame,text="n-fold Cross Validation "),entry.nfcvOptTR);
  
  tkgrid(tklabel(TRFrame,text="TreeRank Options"));
  tkgrid(optTRFrame);
  
  
  
  #LeafRank
  ## Scan the environment for new LeafRank functions definition
  LRlist <- ls(pattern="LR.*.def",name=globalenv());
  nameLR <- list(LRCart.def$name,LRsvm.def$name);
  funNameLR <- list(LRCart.def$fun,LRsvm.def$fun);
  optionLR <- list(LRCart.def$opt,LRsvm.def$opt);
  opt2cmdLR <- list(LRCart.def$opt2cmd,LRsvm.def$opt2cmd);
  lbLR <- tklistbox(LRFrame,height=3,width=10,selectmode="single");
  tkinsert(lbLR,"end",LRCart.def$name);
  tkinsert(lbLR,"end",LRsvm.def$name);
  for (x in LRlist){
    obj <- eval(parse(text=x),globalenv())
    nameLR <- c(nameLR,obj$name);
    funNameLR <- c(funNameLR,obj$fun);
    optionLR <- c(optionLR,list(obj$opt));
    opt2cmdLR <- c(opt2cmdLR,list(obj$opt2cmd));
    tkinsert(lbLR,"end",obj$name)
  }

  tkselection.set(lbLR,0);
  assign("curLeafRank",1,mainGuiEnv);
  tkgrid(tklabel(LRFrame,text="LeafRank: "),lbLR);

  assign("listVarLR",list(),mainGuiEnv);
  assign("listWidgetLR",list(),mainGuiEnv);

  ## LR Frame building  
  buildLRopt <- function(){
    listVarLR <- list()
    id <- as.numeric(tkcurselection(lbLR))+1;
    if (length(id)<1)
      id <- 1;
    listOptLR <- optionLR[[id]];
    listSlaves <- tclvalue(tkgrid.slaves(optLRFrame));
    if (listSlaves[[1]] !=""){
      sapply(unlist(strsplit(listSlaves,split=" ")),tkgrid.remove)
    }
    listWidgetLR <- list()
    listGrid <- list();
    for (i in 1:length(listOptLR)){
      if (listOptLR[[i]]$type=="entry"){
        tmp <- tclVar(listOptLR[[i]]$default)
        listVarLR <- c(listVarLR,list(tmp))
        listWidgetLR<-c(listWidgetLR,list(tkentry(optLRFrame,width="5",textvariable=tmp)))
      }
      if (listOptLR[[i]]$type=="check"){
        tmp <- tclVar(listOptLR[[i]]$default);
        listVarLR <- c(listVarLR,list(tmp))
        listWidgetLR<-c(listWidgetLR,list(tkcheckbutton(optLRFrame,variable=tmp)))
      }
      if (listOptLR[[i]]$type=="listbox"){
        tmplb <- tklistbox(optLRFrame,height=3,width=6,selectmode="single");
        for (na in listOptLR[[i]]$choiceName){
          tkinsert(tmplb,"end",na);
        }
        tkselection.set(tmplb,listOptLR[[i]]$default);
        listVarLR <- c(listVarLR,list(tmplb));
        listWidgetLR<-c(listWidgetLR,list(tmplb))

      }
      tmpW <- listWidgetLR[[length(listWidgetLR)]]
    }
    for (i in (1:(length(listWidgetLR)/2))){
      tkgrid(tklabel(optLRFrame,text=listOptLR[[(i*2-1)]]$name),
             listWidgetLR[[(i*2-1)]],
             tklabel(optLRFrame,text=listOptLR[[(i*2)]]$name),
             listWidgetLR[[(i*2)]]);
    }
    if (length(listWidgetLR) %% 2 >0)
      tkgrid(tklabel(optLRFrame,text=listOptLR[[length(listOptLR)]]$name),
             listWidgetLR[[length(listOptLR)]]);
    assign("listWidgetLR",listWidgetLR,mainGuiEnv)
    assign("listVarLR",listVarLR,mainGuiEnv)
    assign("curLeafRank",id,mainGuiEnv)
    
  }
  
  buildLRopt();
  tkgrid(optLRFrame)
  tkbind(lbLR,"<<ListboxSelect>>",buildLRopt);


  ##Command Line building
  
  cmdLine <- function(){
    dataset1V <-  tclvalue(dataset1);
    dataset2V <- tclvalue(dataset2);
    alphaV <- tclvalue(alpha);
    splitV <- tclvalue(split);
    
    idLR <- get("curLeafRank",mainGuiEnv);
                                        #LeafRank options
    listOptLR <- optionLR[[idLR]]
    optPar <- list();
    listVarLR <- get("listVarLR",mainGuiEnv)
    for(i in 1:length(listVarLR)){
      
      if (listOptLR[[i]]$type=="listbox"){
        id <- as.numeric(tkcurselection(listVarLR[[i]]))+1;
        tmp<- paste("c(optPar,",listOptLR[[i]]["optName"],"=\"",
                    listOptLR[[i]]$choice[id],"\")",sep="");
      }else{
        tmp <- paste("c(optPar,",listOptLR[[i]]["optName"],"=",tclvalue(listVarLR[[i]]),")",sep="");
      }
      optPar <- eval(parse(text=tmp));
      
    }
    optLR <- opt2cmdLR[[idLR]](optPar);
                                        #TreeRank options
    minOptTRV <- as.double(tclvalue(minOptTR));
    maxOptTRV <- as.double(tclvalue(maxOptTR));
    mincritOptTRV <-as.double(tclvalue(mincritOptTR));
    nfcvOptTRV <- as.integer(tclvalue(nfcvOptTR));
    varSplitTRV <- as.numeric(tclvalue(varSplitTR));
    dataSplitTRV <- as.numeric(tclvalue(dataSplitTR));
    grTRctrl <- paste("growing_ctrl(maxdepth=",maxOptTRV,",mincrit=",mincritOptTRV,",minsplit=",minOptTRV,")",sep="");
                                        #    prTRctrl <- paste("pruning_ctrl(prune=",pruneOptTRV,",strat=",stratOptTRV,",nfcv=",nfcvOptTRV,")",sep="");


    listx <- 0;
    listy <- 0;
    pcinit <-0;
    ncinit <-0;
    TRcmd <-  paste("TreeRank(growing=",grTRctrl,",nfcv=",nfcvOptTRV,",LeafRank= function(...){",funNameLR[[idLR]],"(",optLR,",...)},prcSplitVar=",(varSplitTRV/100),",prcSplitData=",(dataSplitTRV/100),",...)",sep="");
    cmd <- paste("TwoSample(x=",dataset1V,",y=",dataset2V,",split=",splitV,
                 ",alpha=",alphaV,",TRalgo=function(...) ",TRcmd,")",sep="");
  return(cmd);
  
  }
  
  runGui <- function(){
    wobj <- eval(parse(text=cmdLine()));
   assign("wobj",wobj,.GlobalEnv);
   id <- TRplot(wobj$tree);
    .TRgui.addROC(id,getROC(wobj$tree,wobj$test),"Test sample");
                                        #  return(eval(parse(text=cmdLine())));
    print(wobj$wtest);
    
  }

  exportCmd <- function(){
    cmd <- cmdLine();
    tt <- tktoplevel()
    name <- tclVar("funName");
    entry.Name <- tkentry(tt,width="20",textvariable=name)
    tkgrid(tklabel(tt,text="Fun name"));
    OnOk <- function()
      {
        NameVal <- tclvalue(name);
        tkdestroy(tt)
        assign(NameVal,cmd,.GlobalEnv);
      }
    OK.but <- tkbutton(tt,text=" set ",command=OnOk);
    tkbind(entry.Name,"<Return>",OnOk)
    tkgrid(entry.Name)
    tkgrid(OK.but)
    tkfocus(tt);
  }

  but.run <-tkbutton(mainFrame,text="run",command =runGui);
  #but.exp <- tkbutton(mainFrame,text="export cmd",command=exportCmd);
  tkgrid(infoFrame);
  
  tkgrid(LRFrame);
  tkgrid(TRFrame);
  tkgrid(but.run);
  #tkgrid(but.exp);
  tkgrid(mainFrame);
  tkfocus(top);
}

errorGui <- function(err){
print(err);
}
