library(marray)
library(GEOquery)
library(Biobase)

dyebias.geo2marray <- function
(gse,                                   #data set
 slide.ids=NULL,                        #return only slides with these ids
 type="norm",                           #what to extract
 gene.selector=function(table){T},      #function acting on Table(GPL) giving back 
                                        # an index with the rows considered to be genes
 reporterid.name,                       # column containing the reporter.id, in Table(gpl)
 cy3.name,                              # the column names containing the factor values
 cy5.name,

 
 ## column names for extracting data from Table(gsm):
 R.name=NULL,
 G.name=NULL,
 M.name=NULL,
 A.name=NULL,
 Rf.name=NULL,
 Gf.name=NULL,
 Rb.name=NULL,
 Gb.name=NULL
 ) {                                    #to corresponding columns

  here <- ", in geoutils.R:geo2marray"

  my.is <- function(thing, expected.class){      # is() from methods gives load errors ...
    if(is.null(thing)){return(FALSE);}
    class <- attributes(thing)$class
    if(is.null(class)){return(FALSE)}
    return (class==expected.class)
  }

  if( ! my.is(gse,"GSE") ) {
    stop('Need a GSE object here (as returned from getGEO("GSExyz", GSEMatrix=FALSE)',
         here, call. =TRUE)
  }

  ## layout and gene info
  if (is.null(slide.ids)) {
    slide.subset <- T
    n.slides=length(GSMList(gse))
  } else {
    if(any(duplicated(slide.ids))) { stop("slide.subset contains duplicates", here, call.=TRUE) }
    slide.subset=(sapply( GSMList(gse),
                         function(x){Meta(x)$geo_accession}) %in% slide.ids)
  }

  # use GPL of first of the slides in this set (lets hope the rest uses the same ... not checked)
  gpl.id <- Meta((GSMList(gse)[slide.subset])[[1]])$platform_id
  if (is.null(gpl.id)) { stop("gpl.id not given", here, call.  =TRUE) }
  gpl <- GPLList(gse)[[gpl.id]]

  ## n.spots <- as.numeric(Meta(gpl)[["data_row_count"]]) # or flength(Table(gpl))??
  ## may not match; GPL884 and GPL2883 both claim 101*107=10807 rows, but have 10806!
  table <- Table(gpl)
  table$ID <- as.character(table$ID)    #they meant ID_REF ...
  table <- table[order(table$ID),]
  
  ## HOWEVER: not all spots need be present in Table(gsm), so filter them out:
  found.ids <- sort(as.character((Table(GSMList(gse)[slide.subset][[1]]))$ID_REF))
  table <- table[ table$ID %in% found.ids, ]
  
  ## filter out the non-Genes
  genes <- gene.selector(table)
  table <- table[genes,]
  spot.ids <- table$ID

  found.ids <- found.ids[ found.ids %in% spot.ids] # 

  reporter.info <- merge(data.frame(ID=found.ids), table, by="ID", sort=TRUE) # now one row per valid spot
  
  spot.ids <- found.ids                  #now has the correct ordering

  n.spots<- length(reporter.info[[1]])
  
  n.cols <- 1                              # was length(unique(table$COL))
  n.rows <- n.spots                        # was length(unique(table$ROW))

  layout <- read.marrayLayout(nsr=1, nsc=1, ngr=n.rows, ngc=n.cols,
                              notes=paste("converted from GSE", Accession(gpl)))
  maControls(layout) <- rep("Gene", n.spots) #only thing left

  gnames <- new("marrayInfo")
  reporter.ids <- as.character(reporter.info[[reporterid.name]])
  reporter.info$reporterId <- reporter.ids
  maLabels(gnames) <- reporter.ids
  maInfo(gnames) <- reporter.info
  
  # targets:
  pdata <- .gsmlist2targetinfo(gse, slide.subset)
  
  ## cols=c("id", "type", "molecule_ch1", "characteristics_ch1", "label_ch1", "molecule_ch2", "characteristics_ch2", "label_ch2", "description")
  cols <- c("id", cy3.name, cy5.name, "type", "description")
  if ( length( setdiff(cols, names(pdata))) > 0 ) {
    missing <- setdiff(cols, names(pdata))
    stop("Could not find these column(s):", paste(missing, collapse=","), here, call. =TRUE)
  }
  
  target.info = pdata[ , cols]
  names(target.info)=c("filename", "Cy3", "Cy5", "type", "description")
  
  target.info$slide <- 1:length(target.info[[1]])
  targets <- new("marrayInfo")
  maInfo(targets) <- target.info
  maLabels(targets) <- target.info$filename

  ## measurements:
  if (type=="norm") { 
    if (  !is.null(R.name) && !is.null(G.name)
        && is.null(M.name) && is.null(A.name) ) {
      ## got R and G, calculated M and A
      R <- .extract.measurements(gse=gse,
                                slide.subset=slide.subset,
                                column.name=R.name,
                                spot.ids=spot.ids)
      
      G <- .extract.measurements(gse=gse,
                                slide.subset=slide.subset,
                                column.name=G.name,
                                spot.ids=spot.ids)
      
      M <- log2(R)-log2(G)
      A <- (log2(R)+log2(G))/2
      W <- matrix(1, nrow=nrow(A), ncol=ncol(A))
      data<-new("marrayNorm", maR=R, maG=G, maM=M, maA=A, maW=W )
    } else if (   !is.null(M.name) && !is.null(A.name)
             && is.null(R.name) && is.null(G.name) ) {
      ## got M and A, calculate R and G
      M <- .extract.measurements(gse=gse,
                                slide.subset=slide.subset,
                                column.name=M.name, #e.g. "VALUE"
                                spot.ids=spot.ids)
      A <- .extract.measurements(gse=gse,
                                slide.subset=slide.subset,
                                column.name=A.name,
                                spot.ids=spot.ids)
      W <- matrix(1, nrow=nrow(A), ncol=ncol(A))
      expA=2^A
      expM2=2^(M/2)
      R= expA*expM2
      G= expA/expM2
      data<-new("marrayNorm", maA=A, maM=M, maR=R, maG=G, maW=W)
    } else if ( !is.null(M.name) && is.null(A.name) ) {    #only M given: calculate A later from raw object

      M <- .extract.measurements(gse=gse,
                                slide.subset=slide.subset,
                                column.name=M.name, #e.g. "VALUE"
                                spot.ids=spot.ids)

      W <- matrix(1, nrow=nrow(M), ncol=ncol(M))
      data<-new("marrayNorm", maM=M, maW=W)
      warning("Norm object: only M given, add A later from a raw object",here, call. = TRUE)
    } else {
      stop("Need either both R.name and G.name, or both A and M.name", here, call.=TRUE)
    }
  } else if (type == "raw") {

    if(is.null(Rf.name) || is.null(Rf.name)) {
      stop("For raw data, need an Rf.name and a Gf.name arguments", here, call.=TRUE)
    }

    Rf <- .extract.measurements(gse=gse,
                               slide.subset=slide.subset,
                               column.name=Rf.name,
                               spot.ids=spot.ids)
    
    Gf <- .extract.measurements(gse=gse,
                               slide.subset=slide.subset,
                               column.name=Gf.name,
                               spot.ids=spot.ids)
    W <- matrix(1, nrow=nrow(Rf), ncol=ncol(Rf))
    

    if(! is.null(Rb.name) && !is.null(Gb.name)  ) {
      Rb <- .extract.measurements(gse=gse,
                                 slide.subset=slide.subset,
                                 column.name=Rb.name,
                                 spot.ids=spot.ids)
      
      Gb <- .extract.measurements(gse=gse,
                                 slide.subset=slide.subset,
                                 column.name=Gb.name,
                                 spot.ids=spot.ids)

      if(.have.umcu.version())
        data<-new("marrayRaw", maRf=Rf, maGf=Gf, maRb=Rb, maGb=Gb, maW=W, maBgSubtract=FALSE)
      else
        data<-new("marrayRaw", maRf=Rf, maGf=Gf, maRb=Rb, maGb=Gb, maW=W)
      
    } else { 

      if(.have.umcu.version())
        data<-new("marrayRaw", maRf=Rf, maGf=Gf,                   maW=W, maBgSubtract=FALSE)
      else
        data<-new("marrayRaw", maRf=Rf, maGf=Gf,                   maW=W)

      n.cols=ncol(maRf(data))
      n.rows=nrow(maRf(data))
      maRb(data) = matrix(0, ncol=n.cols, nrow=n.rows)
      maGb(data) = matrix(0, ncol=n.cols, nrow=n.rows)
      warning(here,": setting missing backgrounds of raw data to 0", call. =TRUE)
    }
  } else {
    stop("Unknown type, must be 'raw' or 'norm', got: ", type, here, call. =TRUE)
  }


  ## put everything together:
  maTargets(data) <- targets
  maLayout(data) <- layout
  maGnames(data) <- gnames
  
  data
}                                       # geo2marray

## ===========  support functions ======================================

.extract.meta <- function
## to extract one Meta 'column' from the gse's list of gsm's and return as a vector
(gse, slide.subset=TRUE, column.name) {
  here="geoutils.R:extract.meta"
  sapply(GSMList(gse)[slide.subset],
         function(x) {
           thing=Meta(x)[[column.name]]
           # if(is.null(thing)) { stop("No such column:",column.name,", in:",here, call.=TRUE) }
           if(is.null(thing)) { return(NA) }
           if(length(thing) != 1) { return(NA) }
           if(is.factor(thing)) { return(as.character(thing))
                                } else {return(thing)}})
}                                       # .extract.meta

.extract.measurements <- function
## to extract one measurement 'column' from the gse's list of gsm's and return as a matrix
## of the type needed by marray. Note that the spots are orderd by ID_REF, since
## it turns out that there are data sets with different spot order per slide (even if the
## platform is identical !@!@!)
(gse, slide.subset=TRUE, spot.ids=NULL, column.name) {
  here <- ", in:  geoutils.R:.extract.measurements"
  if(is.null(column.name)) { stop("column name is null",here, call. =TRUE) }
  list <- sapply( GSMList(gse)[slide.subset],
         function(x, spot.ids) {
           id.name <- "ID_REF"          # fixed name , it seems
           table <- (Table(x))[order((Table(x))[, id.name]),] # ordering may not be homogeneous ?!?! 
           if(is.null(spot.ids)) {
             spot.subset <- T
           } else {
             spot.subset <- table[, id.name] %in% spot.ids
           }
           result <- as.numeric(as.character(table[spot.subset, column.name]))
           if(is.null(result)) { stop("No such column (or wrong type):",column.name,here, call.=TRUE) }
           result #  NA's will give errors in the warnings list
         },
         spot.ids)                      #arg to the func.
  try( frame <- as.data.frame(list), silent=FALSE) 
  if ( class(frame)=="try-error") {      # unequal number of columns
    stop("Failed to cast ", column.name, " of all slides to decent data.frame: ", frame,
         "\nyou prolly have different GPL ids with differnt numbers of spots\n",
         "maybe use spot.ids?\n",here,
         call. =TRUE)
  }
  as.matrix(frame)
}                                       # .extract.measurements

.gsmlist2targetinfo <- function(gse, slide.subset=TRUE) {
  ## take gse, and produce data.frame with all information that is singular
  ## (i.e. not vectors of strings describing protocols). 
  ## for now, it is assumed that the columns found for the first slide are
  ## valid and useful for all slides in the list; this need not be the case.
  ## (in which case they should end up being NA)

  gsmlist=GSMList(gse)[slide.subset]

  l=list(id=names(gsmlist))             #id duplicates geo_accession
  ## id=sapply(GSMList(gse), Accession)

  columns=names(Meta(gsmlist[[1]]))
  for (col in columns) {
    l[[col]] <- .extract.meta(gse=gse, slide.subset=slide.subset, column.name=col)
  }
  
  useless <- names(which(sapply(l,
                                function(x) {
                                  all(is.na(x)) | all(is.null(x)) | all(x=="")
                                })))
  for(col in useless) {
    l[[col]]=NULL
  }
  
  return(as.data.frame(l))
}                                       # .gsmlist2targetinfo

.have.umcu.version <- function() {
## tell if we're running our own modified version of marray or not
  return(package.version(pkg="marray")=="1.5.8")
}
