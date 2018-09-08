##########################################
###    msmsTests - GUI - Version 1.1   ###
##########################################

library(msmsTests)
library(MASS)
options(width=120)
require(gWidgets)
options(guiToolkit="RGtk2")

help.text <- c(
    "\nPlease proceed as follows:\n",
    "\n 1- Set experiment name",    
    "\n 2- Select a samples description file name",
    "\n    (defaults will be set and displayed)",
    "\n 3- Select a spectral counts file name",
    "\n 4- Select a normalization method",
    "\n 5- Select the statistical test for differential expression",
    "\n 6- Select significance and a p-value multitest correction method",
    "\n 7- Give the null model formula as factor names separed by '+'",
    "\n 8- Give the alternative model formula as before",
    "\n    (note that null and alternative should differ in a single factor)",
    "\n 9- Specify the treatment condition",
    "\n19- Specify the control condition",
    "\n11- Specify the post-test filter thresholds",
    "\n12- Specify the volcano plot limits",
    "\n13- Check everything again and click on 'Execute' to start\n")

##  Data global to controls
samples.flnm <- data.flnm <- NULL
tit="MSMS"

###  Construct a MSnSet object from SpC matrix and factors
SpCm2MSnSet <- function(counts,facts)
{ e <- new("MSnSet", exprs=counts)
  pData(e) <- facts
  e <- pp.msms.data(e)
  return(e)
}  

###  EXECUTION
msms.tests <- function()
{
  ###  Unique (experiment) name 
  tit <- get.expnm()
  if( nchar(tit)==0 ) 
  { tit <- paste( paste(sample(LETTERS,5),collapse=""), 
                          format(Sys.time(),"%d%b%Y"),sep="." )
    svalue(g.exp.nm) <- tit
  }
  etit <- tit
  pdf.flnm <- paste(tit,"_vplot.pdf",sep="")
  hmpdf.flnm <- paste(tit,"_heatmap.pdf",sep="")
  fhmpdf.flnm <- paste(tit,"_big_heatmap.pdf",sep="")
  txt.flnm <- paste(tit,"_tests.txt",sep="")

  ###  Read samples description
  samples.flnm <- svalue(g.smpls.flnm)
  if( !file.exists(samples.flnm) )
  { galert("Samples file not found!",cont=window)
    reset.smpls()
    return()
  }
  samples <- read.table(file=samples.flnm,sep="\t",header=TRUE,dec=",")
  cnms <- colnames(samples)
  if(cnms[1]!="Sample")
  { galert("Wrong samples file: No 'Sample' column!",cont=window)
    reset.smpls()
    return()
  }
  samples$Sample <- as.character(samples$Sample)
  if(!is.null(samples$offsets))
    samples$offsets <- as.numeric(as.character(samples$offsets))
    
  ###  Check if cell counts / secreted protein given
  jdx <- 1
  cell.j <- which(cnms=="offsets")
  if(length(cell.j)>0) jdx <- c(jdx,cell.j)

  ###  Build factors matrix
  facs <- data.frame(samples[,-jdx,drop=FALSE])
  rownames(facs) <- as.character(samples$Sample)

  ###  Read dataset
  data.flnm <- svalue(g.msms.flnm)  
  if( !file.exists(data.flnm) )
  { galert("Counts file not found!",cont=window)
    reset.msms()
    return()
  }
  
  msms <- read.table(file=data.flnm,header=TRUE,sep="\t",
                     stringsAsFactors=FALSE)
 
  ###  Check for data consistency
  if( length(intersect(colnames(msms),samples$Sample)) <
      nrow(samples) )
  { galert("\nSamples description file does not match data file.\n",
           cont=window)
    reset.msms()
    return()
  }
  if( ! "Accession" %in% colnames(msms) )
  { galert("\nWrong counts file: no 'Accession' column in data file.\n",
           cont=window)
    reset.msms()
    return()
  }
  if( ! "Proteins" %in% colnames(msms) )
  { if( !gconfirm("\nNo 'Proteins' column in data file. Proceed ?\n",
             title="WARNING",icon=c("question"),cont=window) )
    { reset.msms()
      return()
    }
  }

  if( (svalue(g.norm) == "Bio" | svalue(g.norm) == "Size & Bio" ) & 
                   is.null(samples$offsets) )
  { galert("\nNo offsets given for normalization.\n",
           cont=window)
    return()
  }

  ###  Read / check parameters
  test.nm <- svalue(g.test)
  alpha <- as.numeric(svalue(g.alpha))
  svalue(g.alpha) <- alpha
  mt.corr <- svalue(g.mtc)
  form.Ha <- paste("y~",svalue(g.Ha),sep="")
  form.Ho <- paste("y~",svalue(g.Ho),sep="")
  minSpC <- as.numeric(svalue(g.mnspc))
  if(is.na(minSpC) | minSpC<0) minSpC <- 0
  svalue(g.mnspc) <- minSpC
  minLFC <- as.numeric(svalue(g.mnlfc))
  if(is.na(minLFC) | minLFC<0) minLFC<- 0
  svalue(g.mnlfc) <- minLFC
  condA <- svalue(g.condA)
  condB <- svalue(g.condB)

  ###  Check factors and identify treatment factor (fnm)
  fnm <- svalue(g.Ho)
  n.fnms <- c("1",strsplit(fnm,"\\+")[[1]])
  if(length(setdiff(n.fnms,c("1",colnames(facs)))))
  { galert("\nUnknown factor in null model.\nPlease check.",
           cont=window)
    return()
  }
  fnm <- svalue(g.Ha)
  a.fnms <- c("1",strsplit(fnm,"\\+")[[1]])
  if(length(setdiff(a.fnms,c("1",colnames(facs)))))
  { galert("\nUnknown factor in alternative model.\nPlease check.",
           cont=window)
    return()
  }
  fnm <- setdiff(a.fnms,n.fnms)
  if(length(fnm)>1)
  { galert("\nNull and alternative models should differ in a single factor.\nPlease check.",
           cont=window)
    return()
  }
  ###  Make sure the levels to contrast exist.
  if( sum(c(condA,condB) %in% levels(facs[,fnm])) < 2 )
  { galert("\nWrong conditions to compare.\nPlease check.",
           cont=window)
    return()
  }
  ###  Set control condition as ref.
  facs[,fnm] <- relevel(facs[,fnm],ref=condB)

  ### Subset to needed samples
  fl <- (facs[,fnm]==condA | facs[,fnm]==condB)
  msms.counts <- data.matrix(msms[ ,samples$Sample[fl]])
  rownames(msms.counts) <- msms$Accession
  facs <- facs[fl,,drop=FALSE]
  nsmpl <- sum(fl)
  smpls <- samples[fl,,drop=FALSE]
  
  ###  Short sample names to be used in graphics
  short.nms <- paste("X",1:nsmpl,sep="")

  ###  Create MSnSet object, and preprocess it by setting NAs to 0,
  ###    and removing all 0 rows and -R proteins.
  e <- SpCm2MSnSet(msms.counts,facs)
  if(!validObject(e))
   { galert("\nInvalid MSnSet object.\nPlease check.",
           cont=window)
    return()
  }
  msms.counts <- exprs(e)
  
  ###  Check for a minimum dimensionality (just warn)
  if( nrow(msms.counts) < 100 )
    if ( ! gconfirm("\nExpression matrix with less than 100 rows.\nProceed ?",
              title="WARNING",icon=c("question"),cont=window) ) 
                     return()

  if( !gconfirm("Data checked.\nConfirm execution, please.  ",
       title="DATA CHECKED", icon=c("question")) )
          return()

  insert(output,"\n\nStarting execution")
  insert(output,"------------------\n")
  insert(output,paste("Normalization:",svalue(g.norm)))
  insert(output,paste("Statistic test:",test.nm))
  insert(output,paste("Null hypothesis model:",form.Ho))
  insert(output,paste("Alternative hypothesis model:",form.Ha))
  insert(output,paste("Control:",condB))
  insert(output,paste("Treatment:",condA))
  insert(output,paste("Adjusted significance:",alpha))
  insert(output,paste("Multest correction method:",mt.corr))
  insert(output,paste("Post-test filter:  Spc",minSpC,"  LogFC",minLFC))

  ###  Gene names and accession table  
  gn.tbl <- NULL
  if(!is.null(msms$Protein))
    gn.tbl <- gene.table(msms$Accession,msms$Protein,
                         pat="GN=[A-Z0-9_]*",off=3)

  ###  Compute normalizing factors
  norm <- svalue(g.norm)
  # None
  div <- rep(1,nsmpl)
  # Size
  if( norm=="Size" | norm=="Size & Bio")
    div <- apply(msms.counts,2,sum)
  # Bio
  if( norm=="Bio" )
    div <- smpls$offsets
  # producció de proteïna per cel•lula
  if( norm=="Size & Bio" )
    div <- div*as.numeric(smpls$offsets)

  ###  Estadístics de normalitzacuó
  if( norm!="None" )
  { insert(output,"\nNormalizing divisors statistics:\n")
    div.m <- tapply(div,facs[,fnm],mean)
    div.max <- tapply(div,facs[,fnm],max)
    div.min <- tapply(div,facs[,fnm],min)
    div.sd <- tapply(div,facs[,fnm],sd)
    div.rg <- round((div.max-div.min)/div.min*100,2)
    df <- data.frame(min=div.min,max=div.max,mean=div.m,sd=div.sd,
                     pct.range=div.rg)
    sink("tmp.txt")
    print(df)
    sink()
    show.file("tmp.txt")
  }
  
  ###  Statistic tests
  insert(output,"\nRunning tests")
  if(test.nm=="QL")
    tres <- msms.glm.qlll(e,form.Ha,form.Ho,facs,div)
  if(test.nm=="Poisson")
    tres <- msms.glm.pois(e,form.Ha,form.Ho,facs,div)
  if(test.nm=="NB-EdgeR")
    tres <- msms.edgeR(e,form.Ha,form.Ho,facs,div,fnm)

  lres <- test.results(tres,e,facs[,fnm],condA,condB,div,
                    alpha=alpha,minSpC=minSpC,minLFC=minLFC,method=mt.corr)
  nms <- rownames(lres$tres[lres$tres$adjp<=alpha,])
  DEP.tbl <- lres$tres[nms,]
  if(!is.null(gn.tbl))
    DEP.tbl <- data.frame(Prot.Nm=gn.tbl[nms],lres$tres[nms,],
                     stringsAsFactors=FALSE)
  ###  Salvar en format CSV
  rflnm <- paste(tit,"-results.csv",sep="")
  write.table(DEP.tbl,file=rflnm,row.names=TRUE,col.names=NA,
              dec=",",sep="\t")

  ###  Taula de LogFC per p-valor ajustat
  sink("tmp.txt")
  print( sig.tbl <- pval.by.fc(lres$tres$adjp,lres$tres$LogFC) )
  sink()
  insert(output,"\nTable of LogFC by p-values:\n")
  show.file("tmp.txt")

  deps <- sum(lres$tres$adjp<alpha)
  insert(output,paste("\n",deps," DEPs at ",alpha,sep=""))
  deps <- sum(lres$tres$DEP)
  insert(output,paste(deps," DEPs passing the post test filter",
         sep=""))

  ###  Show top 100
  sink("tmp.txt")
  print(lres$tres[1:100,])
  sink()
  insert(output,"\nTop 100 features:\n")
  show.file("tmp.txt")

  ###  Save data for later use
  rflnm <- paste(tit,"-test.data.RData",sep="")
  save(samples,norm,test.nm,form.Ho,form.Ha,condA,condB,alpha,mt.corr,minSpC,
       minLFC,msms.counts,gn.tbl,facs,div,lres,sig.tbl,tit,file=rflnm)
  file.remove("tmp.txt")

  ###  Save output
  sv <- unlist(strsplit(svalue(output),"\\n"))
  sink(txt.flnm)
  cat("\nLC-MS/MS EXPERIMENT:",tit,"\n")
  i0 <- length(help.text)+4
  for(i in i0:length(sv)) cat(sv[i],"\n")
  sink()

  ###  Draw heatmap for significant proteins using the normalized signal.
  msms.sig <- msms.counts[rownames(lres$tres[lres$tres$DEP,]),]
  msms.sig <- sweep(msms.sig,MARGIN=2,STATS=div,FUN="/")
  pdf(file=hmpdf.flnm,paper="a4",width=7.5,height=11)
  library(gplots)
  gtit <- paste("Heat map - ",tit,sep="")
  heatmap(t(scale(t(msms.sig))),col=greenred(255),labRow=NA,
              cexCol=0.7)
  title(main=gtit,line=1,cex=1)
  dev.off()

  h <- nrow(msms.sig)/(2.54/0.3)
  pdf(file=fhmpdf.flnm,width=7,height=h)
  gtit <- paste("Heat map - ",tit,sep="")
  heatmap.2(t(scale(t(msms.sig))),col=greenred(255),trace="none",key=FALSE,
            cexRow=0.6,cexCol=0.7,margins=c(5,6),dendrogram="row",
            lhei=c(1,h-1))
  title(main=gtit,line=1,cex=1)
  dev.off()

  ###  Draw volcano plot
  pdf(file=pdf.flnm,paper="a4",height=5.5,width=5.5)
  maxx <- as.numeric(svalue(g.vpxlm))
  if(is.na(maxx)) maxx <- 2
  if(maxx<2) maxx <- 2
  svalue(g.vpxlm) <- as.character(maxx)
  maxy <- as.numeric(svalue(g.vpylm))
  if(is.na(maxy)) maxy <- 2
  if(maxy<1) maxy <- 2
  svalue(g.vpylm) <- as.character(maxy)
  res.volcanoplot(lres$tres,max.pval=alpha,min.LFC=minLFC,
                  maxx=maxx,maxy=maxy,ylbls=100)
  dev.off()

  save(lres,file="results.tmp")
  enabled(a.replot) <- TRUE

  insert(output,"\n = = =  D O N E  = = =")
}

##  Functions used by handlers
show.file <- function(flnm,n=-1L)
{ insert(output,readLines(flnm,n),
     font.attr=c(family="monospace",style="tty"))
}

get.models <- function(samples.flnm)
{
  samples <- read.table(file=samples.flnm,sep="\t",header=TRUE,dec=",")
  jdx <- 1
  cnms <- colnames(samples)
  cell.j <- which(cnms=="offsets")
  if(length(cell.j)>0) jdx <- c(jdx,cell.j)
  cnms <- cnms[-jdx]
  if(length(cnms) < 2)
  { svalue(g.Ho) <- "1"
  } else svalue(g.Ho) <- paste(cnms[-1],collapse="+")
  svalue(g.Ha) <- paste(cnms,collapse="+")
  lev <- levels(factor(samples[,cnms[1]]))
  if(length(lev)>1)
  { svalue(g.condA) <- lev[2]
    svalue(g.condB) <- lev[1]
  }
}

reset.smpls <- function()
{
  enabled(a.exec) <- FALSE
  enabled(a.replot) <- FALSE
  svalue(g.smpls.flnm) <- ""
  svalue(g.Ho) <- ""
  svalue(g.Ha) <- ""
}

reset.msms <- function()
{
  enabled(a.exec) <- FALSE
  enabled(a.replot) <- FALSE
  svalue(g.msms.flnm) <- ""
}


##  Get current contents of g.exp.nm gedit control
get.expnm <- function()
{ svalue(g.exp.nm) }

##  Handlers
h.exec <- function(h,...)
{ msms.tests() }

h.about <- function(h,...)
{ cms <- c("Lable-free LC-MS/MS SpC-based Biomarker Discovery     ",
           "",
           "                     StatsCountsMS_Tests v 1.30",
           "",
           "See:",
		   "   Gregori J. et al. (2012) J. Proteomics 75, 3938-3951",
		   "   Gregori J. et al. (2013) J. Proteomics 95, 55-65",
		   "   Gregori J. et al. (2014) J. Proteome Res. Jun 4",
		   "",
		   "(C) J. Gregori, J. Villanueva, A. Sanchez, UB - VHIO, 2014-2018")
  amsg <- paste(cms,collapse="\r\n")		   
  gmessage(amsg,title="About",icon="info") 
}

h.smpls <- function(h,...)
{ samples.flnm <<- svalue(h$obj)
  if( file.exists(samples.flnm) )
  { insert(output,"\n\nSamples description:\n")
    show.file(samples.flnm)
    get.models(samples.flnm)
  }
  if(!is.null(samples.flnm) && file.exists(samples.flnm) &&
     !is.null(data.flnm) && file.exists(data.flnm) )
    enabled(a.exec) <- TRUE
}

h.msms <-function(h,...)
{ data.flnm <<- svalue(h$obj)

  if( file.exists(data.flnm) )
  { insert(output,"\n\nExpression matrix column names:\n")
    show.file(data.flnm,1)
  }

  if(!is.null(samples.flnm) && file.exists(samples.flnm) &&
     !is.null(data.flnm) && file.exists(data.flnm) )
    enabled(a.exec) <- TRUE
}    

h.clear <- function(h,...)
{ svalue(output) <- ""
  svalue(g.smpls.flnm) <- "Select a file ..."
  svalue(g.msms.flnm) <- "Select a file ..."
  samples.flnm <<- NULL
  data.flnm <<- NULL
  enabled(a.exec) <- FALSE
  enabled(a.replot) <- FALSE
  txt <- paste(help.text,collapse="") 
  svalue(output) <- txt
}

h.replot <- function(h,...)
{ if(!file.exists("results.tmp")) return
  load("results.tmp")
  if(is.null(lres)) return
  tit <- get.expnm()
  if( nchar(tit)==0 ) tit <- "MSMS"
  pdf.flnm <- paste(tit,"_vplot.pdf",sep="")
  pdf(file=pdf.flnm,paper="a4",height=5.5,width=5.5)
  maxx <- as.numeric(svalue(g.vpxlm))
  if(is.na(maxx)) maxx <- 2
  if(maxx<2) maxx <- 2
  svalue(g.vpxlm) <- as.character(maxx)
  maxy <- as.numeric(svalue(g.vpylm))
  if(is.na(maxy)) maxy <- 2
  if(maxy<1) maxy <- 2
  svalue(g.vpylm) <- as.character(maxy)
  res.volcanoplot(lres$tres,maxx=maxx,maxy=maxy,ylbls=100)
  dev.off()
}

##  Main window
window <- gwindow("StatsCountsMS - SpC Differential Expression", 
                   visible=FALSE)

### Per encabir missatge a peu de la GUI
topg <- ggroup(cont = window, horizontal = FALSE)

###  La part operativa de la GUI
outg <- ggroup(cont = topg, horizontal = TRUE)

maing <- ggroup(cont = outg, horizontal = FALSE,expand=TRUE)

## top group organized horizontally
group <- ggroup(cont = maing, horizontal = TRUE)

## Group for file selection organized vertically
grfl <- ggroup(cont = group, horizontal = FALSE, spacing=5)


glabel("Samples file name:",cont = grfl, anchor = c(-1,0))
g.smpls.flnm <- gfilebrowse(text = "Select a file ...",width=55,
                        quote = FALSE, handler = h.smpls,
                        type = "open", cont = grfl)
addSpace(grfl,12)
glabel("Counts file name:",cont = grfl, anchor = c(-1,0))
g.msms.flnm <- gfilebrowse(text = "Select a file ...",width=55,
                        quote = FALSE, handler = h.msms,
                        type = "open", cont = grfl)
addSpace(grfl,18)

##  Separator
addSpace(group,12)
gseparator(horizontal=FALSE,cont=group)
addSpace(group,12)
#addSpring(group)

## Group for buttons and experiment name
grrg <- ggroup(cont = group, horizontal = FALSE)

## Group for the experiment name
#addSpace(grrg,5,horizontal = FALSE)
grenm <- ggroup(cont = grrg, horizontal = TRUE)
glabel("Experiment name:",cont = grenm, anchor = c(-1,0))
g.exp.nm <- gedit(text=tit,width=15,cont=grenm)
##  el valor es recuperarà dins l'execució com svalue(g.exp.nm)

## Group for the normalization method
addSpace(grrg,5,horizontal = FALSE)
grnrm <- ggroup(cont = grrg, horizontal = TRUE)
addSpace(grnrm,24,horizontal = TRUE)
glabel("Normalization:",cont = grnrm, anchor = c(-1,0))
g.norm <- gcombobox(c("None","Size","Bio","Size & Bio"),cont=grnrm,selected=1)

## Group for clear, exec & about button
addSpring(grrg,horizontal=FALSE)
grbt <- ggroup(cont = grrg, horizontal = TRUE)
addSpace(grbt,10)
## A button to clear dialog
clr_button <- gbutton("Clear",handler=h.clear, cont = grbt)

## A button to initiate the execution
addSpace(grbt,20)
a.exec <- gaction("Execute",handler=h.exec)
g.exec <- gbutton(text="Execute",action=a.exec, cont = grbt)
enabled(a.exec) <- FALSE

## A button to open the about dialog
addSpace(grbt,20)
about <- gaction("About",handler=h.about)
about <- gbutton(action=about, cont = grbt)
addSpace(grbt,10)

##  Separator
addSpace(group,18)

## Area for output
frame <- gframe(" Output & Status ", cont = maing, horizontal = FALSE)
txt <- paste(help.text,collapse="") 
output <- gtext(txt, cont = frame, expand = TRUE,
                 font.attr=c(family="monospace"))
size(output) <- c(500, 700)


## Group for execution parametets
prfr.dec <- gframe("", cont = outg, horizontal = TRUE, expand=FALSE)
addSpace(prfr.dec,8)
prfr  <- ggroup(cont=prfr.dec, horizontal = FALSE, expand=FALSE)

addSpace(prfr,12)
g1 <- ggroup(cont = prfr, horizontal = TRUE)
glabel("Statistic test",cont=g1,anchor=c(-1,0))
g.test <- gcombobox(c("Poisson","QL","NB-EdgeR"),cont=g1,selected=1)

addSpace(prfr,12)
g2 <- ggroup(cont = prfr, horizontal = TRUE)
glabel("Significance",cont = g2, anchor = c(-1,0))
g.alpha <- gedit(text="0.05",width=5,cont=g2)
addSpring(g2)

addSpace(prfr,12)
g3 <- ggroup(cont = prfr, horizontal = TRUE)
glabel("MT correction",cont=g3,anchor=c(-1,0))
g.mtc <- gcombobox(c("BH"),cont=g3)

addSpace(prfr,12)
g11 <- ggroup(cont = prfr, horizontal = FALSE)
glabel("Ho model",cont = g11, anchor = c(-1,0))
g.Ho <- gedit(text="",width=25,cont=g11)

addSpace(prfr,12)
g12 <- ggroup(cont = prfr, horizontal = FALSE)
glabel("Ha model",cont = g12, anchor = c(-1,0))
g.Ha <- gedit(text="",width=25,cont=g12)

addSpace(prfr,12)
g13 <- ggroup(cont = prfr, horizontal = FALSE)
glabel("Treatment condition",cont = g13, anchor = c(-1,0))
g.condA <- gedit(text="",width=25,cont=g13)

addSpace(prfr,12)
g14 <- ggroup(cont = prfr, horizontal = FALSE)
glabel("Control condition",cont = g14, anchor = c(-1,0))
g.condB <- gedit(text="",width=25,cont=g14)


##  Separator
addSpace(prfr,18)
gseparator(horizontal=TRUE,cont=prfr)
addSpace(prfr,18)

glabel("Post test filter",cont=prfr,anchor=c(-1,0))

addSpace(prfr,6, horizontal = FALSE)
glabel("Most abundant min mean SpC",cont=prfr,anchor=c(-1,0))
g4 <- ggroup(cont = prfr, horizontal = TRUE)
g.mnspc <- gspinbutton(from=0,to=5,by=0.01,cont=g4,value=2)
addSpring(g4)

addSpace(prfr,6, horizontal = FALSE)
glabel("Min absolute log fold change",cont=prfr,anchor=c(-1,0))
g5 <- ggroup(cont = prfr, horizontal = TRUE)
g.mnlfc <- gspinbutton(from=0,to=1.5,by=0.05,cont=g5,value=0.8)
addSpring(g5)

##  Separator
addSpace(prfr,18)
gseparator(horizontal=TRUE,cont=prfr)
addSpace(prfr,18)

glabel("Volcanoplot limits",cont=prfr,anchor=c(-1,0))

addSpace(prfr,6, horizontal = FALSE)

glabel("Log FC",cont=prfr,anchor=c(-1,0))
g6 <- ggroup(cont = prfr, horizontal = TRUE)
g.vpxlm <- gedit(text="5",width=4,cont=g6)
addSpring(g6)

addSpace(prfr,6, horizontal = FALSE)
glabel("Log p-value",cont=prfr,anchor=c(-1,0))
g7 <- ggroup(cont = prfr, horizontal = TRUE)
g.vpylm <- gedit(text="10",width=4,cont=g7)
addSpring(g7)
a.replot <- gaction("Replot",handler=h.replot)
g.replot <- gbutton(text="Replot",action=a.replot, cont = g7)
enabled(a.replot) <- FALSE
#addSpring(g7)

addSpace(prfr,18)
gseparator(horizontal=TRUE,cont=prfr)
addSpace(prfr,18)

g8 <- ggroup(cont = prfr, horizontal = TRUE)
addSpring(g8)
gbutton("Quit",cont=g8, handler = function(h,...) dispose(window))
addSpring(g8)

addSpace(prfr.dec,8)

###  Copyright message
glabel("(C) Josep Gregori, UB - VHIO, 2014-2018",cont = topg, anchor = c(-1,0))

visible(window) <- TRUE

