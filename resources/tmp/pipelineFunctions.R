#library(seqMeta) removing dependence on seqMeta
library(plyr)
library(gap)


# This function gets the genotype matrices from a specified location
#
# path: a character vector of full path names; the default corresponds to the 
#       working directory, getwd(). Tilde expansion (see path.expand) is 
#       performed. Missing values will be ignored.
#
# pattern: an optional regular expression. Only file names which match the 
#          regular expression will be returned.
#
# returns: file names with their corresponding paths
#
# GetGenotypeFiles <- function(path=".", pattern=NULL) {
#   list.files(path=path, pattern=pattern, full.names=TRUE)
# }

# This function provides a basic check that the genotype matrix is a numeric
# matrix suitable for seqMeta.
#
checkGenotypeMatrix <- function(Z) {
  if (!is.matrix(Z)) {
    stop("Genotype must be in a numeric matrix.")
  }
  
  # make sure the data is numeric or integer.
  if ((typeof(Z) != "integer") & (typeof(Z) != "double")) {
    stop("Genotype matrix must be numeric (integer or double).")    
  }
  
  return(invisible(NULL)) 
  
}

# This function provides a basic check that the SNPInfo file is a suitable for
# seqMeta.  It check that the required column names are in the SNPInfo data
# frame and that they are of the correct type. 
#
checkSNPInfo <- function(snpinfo, snpNames=NULL, aggregateBy=NULL, filterBy=NULL, chrName=NULL) {
  if (is.null(snpNames) & is.null(aggregateBy)) {
    #check that there are no factors 
    if (any(sapply(snpinfo, is.factor))) {
      stop("Factor columns not allowed in snpinfo.  Convert to characters.")
    }
  }
  
  # check snpNames
  if (!is.null(snpNames)) {
    if (!(snpNames %in% colnames(snpinfo))) {
      stop(paste0("snpNames: ", snpNames, " is not a column name of the provided SNPInfo file."))
    }
    if (typeof(snpinfo[ , snpNames]) != "character") {
      stop(paste0("snpNames: ", snpNames, "is of type ", typeof(snpinfo[ , snpNames]), ".  Must be of type character"))
    }
  }
  
  # check aggregateBy
  if (!is.null(aggregateBy)) {
    if (!(aggregateBy %in% colnames(snpinfo))) {
      stop(paste0("aggregateBy: ", aggregateBy, " is not a column name of the provided SNPInfo file."))
    }
    if (typeof(snpinfo[ , aggregateBy]) != "character") {
      stop(paste0("aggregateBy: ", aggregateBy, "is of type ", typeof(snpinfo[ , aggregateBy]), ".  Must be of type character"))
    }
  }
  
  # check filterBy
  if (!is.null(filterBy)) {
    if (any(!(filterBy %in% colnames(snpinfo)))) {
      stop(paste0("filterBy: ", filterBy[!(filterBy %in% colnames(snpinfo))], " is not a column name of the provided SNPInfo file."))
    }
    if (length(filterBy) == 1L) {
      if (typeof(snpinfo[ , filterBy]) != "logical") {
        stop(paste0("filterBy: ", filterBy, "is of type ", typeof(snpinfo[ , chrName]), ".  Must be of type logical"))
      }      
    } else if (length(filterBy) >= 1L) {
      test <- sapply(snpinfo[, snpinfo.filterBy], typeof)
      if (any(test != "logical")) {
        stop(paste0("filterBy: ", names(test)[test != "logical"], "is of type ", test[test != "logical"], ".  Must be of type logical"))
      }      
    }
  }  
  
  #check chrName
  if (!is.null(chrName)) {
    if (!(chrName %in% colnames(snpinfo))) {
      stop(paste0("chrName: ", chrName, " is not a column name of the provided SNPInfo file."))
    }
    if (typeof(snpinfo[ , chrName]) != "character") {
      stop(paste0("chrName: ", chrName, "is of type ", typeof(snpinfo[ , chrName]), ".  Must be of type character"))
    }
  }  
  
  return(invisible(NULL))   
}


# This function provides a basic check that the phenotype file is a suitable
# for seqMeta.  It checks that the required column names are in the phenotype
# data frame. 
#
checkPhenotype <- function(p, pformula, idCol=NULL, genderCol=NULL) {
  if (!is.null(idCol)) {
    if (any(duplicated(p[ , idCol]))) {
      stop("Duplicated phenotype ids.")
    }
  }
  
  if(!is.null(genderCol)) {
    if (!(genderCol %in% colnames(p))) {
      msg <- paste(genderCol, "not found in phenotype file.", sep=" ")
      stop(msg)
    }
    gtype <- typeof(p[ , genderCol])
    g <- unique(p[ , genderCol])
    if(gtype == "integer") {
      if (!all(g %in% c(0, 1))) {
        stop("Gender must be (0/1 or F/T) indicating female/male.")
      }      
    } else if(gtype == "character") {
      if (!all(g %in% c("F", "M"))) {
        stop("Gender must be (0/1 or F/T) indicating female/male.")
      }        
    }     
  } else {
    wmsg <- "No column given to identify Males in phenotype file.  No special 
         handling of the X chromosme."
    warning(wmsg)
  }
  
  tf <- terms(as.formula(pformula))
  cn <- rownames(attr(tf, "factors"))
  
  if (any(!(cn %in% colnames(p)))) {
    msg <- paste("Formula varaibles:", cn[!(cn %in% colnames(p))], "not found in phenotype file.", sep=" ")
    stop(msg)
  } 
  return(invisible(NULL)) 
}



# This function reduces a data set to only the variables used in a model
# removing subjects with missing data.  Also, it makes the row names of
# the resulting data fram the subject identifier
#
# JAB addition: subsets to complete cases (i.e. no NAs in outcome or covariates)
#
# p: a data frame containing the variables in the model
#
# formula: a character vector which can be coered to an object of class 
#          "formula" (with as.formula): a symbolic description of the model to
#          be fitted. The details of model specification are given under 
#          'Details' in the "lm" help file.
#
# id: (optional) colunm name identifier of the subjects
#
# gender: (optional) colunm name identifier for the gender classification of
#         the subjects.
#
# returns: data frame with only the columns specified in the formula and with
#          the (optional) row names as the subject identifier.
#

reducePheno <- function(p, pformula, id=NULL, gender=NULL) {
  checkPhenotype(p, pformula, idCol=id, genderCol=gender)   
  
  if (!is.null(id)) {
    rownames(p) <- p[ , id]
  }
  
  if (!is.null(gender)) {
    gtype <- typeof(p[ , gender])
    g <- unique(p[ , gender])
    if(gtype == "character") {
      #if (all.equal(g, c("F", "M"))) {
      if (all.equal(sort(g), c("F", "M"))) {
        MF <- p[, gender] == "M"
        p[, gender] <- MF
      } else {
        stop("Unable to convert gender.  Gender must be (0/1 or F/T) indicating female/male.")
      }       
    }                                   
  }

  
  tf <- terms(as.formula(pformula))
  cn <- rownames(attr(tf, "factors"))
  
  if (is.null(cn)) {
    cn <- all.vars(pformula)[1]
  }

  if (!(gender %in% cn)) {
    cn <- c(cn, gender)
  }
  
  if (any(!(cn %in% colnames(p)))) {
    msg <- paste("Formula variables:", cn[!(cn %in% colnames(p))], "not found in phenotype file.", sep=" ")
    stop(msg) 
  } else {
    if (!is.null(cn)) {
      pheno <- na.omit(p[, cn,drop=F])
    } else {
      pheno <- na.omit(p)
    }    
  }

  return(pheno)
}


# This Calculates the MAC without imputing missing values from the cohort
# object and the phenotype data used.  This is done by first calling
# AddMACtoSNPInfo and adding the MAC from the current model to the snpinfo file
#
CalculateMAC <- function(cohort, snpinfo, snpNames, aggregateBy, mafRange) {
  #  check we have deqMeta or skatMeta cohort objects
  if (class(cohort) != "seqMeta" & class(cohort) != "skatCohort") {
    stop("cohort is not a seqMeta object!")
  }
  # check we have the correct columns in our snpinfo file
  if (!(snpNames %in% colnames(snpinfo))) {
    msg <- paste("snapNames: ", snpNames, " does not match any column name in SNPInfo", sep='')
    stop(msg)
  }
  if (!(aggregateBy %in% colnames(snpinfo))) {
    msg <- paste("aggregateBy: ", aggregateBy, " does not match any column name in SNPInfo", sep='')
    stop(msg)
  }
  
  genes <- names(cohort)
  MAC <- integer(length(genes))
  names(MAC) <- genes
  MAC[]<-NA
  
  if ("MAC" %in% colnames(snpinfo)){
    for (gene in genes) {
      x<-cohort[[gene]]
      si <- snpinfo[(snpinfo[, aggregateBy] %in% gene), c(snpNames, aggregateBy, "MAC")]
      
      snps <- unique(names(x$maf)[(x$maf > min(mafRange)) & (x$maf <= max(mafRange))])
      if (length(snps) > 0L) {
        idx <- match(snps, si[, snpNames])
        MAC[gene] <- sum(si[idx, "MAC"], na.rm=TRUE)
      }  
    }
  }  
 
  return(MAC)
  
}



