
args<-commandArgs(TRUE)
print(args)
#===mandatory parameters
phenotype.file <- args[1]
phenotype.formula <- args[2]
snpinfo.file <- args[3]
genotype.files <- args[4]
outputfilename <- args[5]

#==optional parameters
kinship.matrix <- args[6]
genderCol <- args[7]
pheno.id <- args[8]
snpNames <- args[9]
aggregateBy <- args[10]
nsmatch <- args[11]


# added these to JSON
BUFFER <- as.numeric(args[12]) #jb
gene.file <- args[13] #jb
snp.filter <- args[14] #jb
gene.filter <- args[15] #jb
top.maf <- as.numeric(args[16]) #jb
mmax <-  as.numeric(args[17])
TESTTYPE <-  args[18]
TESTFAMILY <-  args[19]



cat('outputfilename',outputfilename,'\n')
cat('kinship.matrix',kinship.matrix,'\n')
cat('genderCol',genderCol,'\n')
cat('aggregateBy',aggregateBy,'\n')
cat('nsmatch',nsmatch,'\n')
cat('buffer',BUFFER,'\n')
cat('gene.file',gene.file,'\n')
cat('snp.filter',snp.filter,'\n')
cat('top.maf',top.maf,'\n')
cat('mmax',mmax,'\n')
cat('TESTTYPE',TESTTYPE,'\n')
cat('TESTFAMILY',TESTFAMILY,'\n')


#library(svd)
library(survey)

#library(seqMeta)
library(SeqArray)
library(SeqVarTools)
library(gap)
library(Matrix)
library(plyr)
library(bdsmatrix)
library(coxme)
library(CompQuadForm)
library(parallel)
library(GENESIS)
library(GWASTools)
#setMKLthreads(1)

library(data.table)
library(doMC)
num.cores <- detectCores(logical=TRUE)
registerDoMC(cores=num.cores-1)
cat('Number of cores', num.cores,'\n')

## Setup
source("/tmp/seqMetaPipelineFunctions.R")
source("/tmp//pchisqsum2.R")
source("/tmp/utils.R")

##  Gene list
cat('reading gene file...')
kg = read.table(gene.file,as.is=T,sep=',',header=T) #jb
cat('GENE Filter',gene.filter,'\n')
kg = eval(parse(text= paste0('subset(kg,',gene.filter,')')))

cat(NROW(kg),'done\n')
genes <- kg$name


## snp info
cat('Reading snpinfo....')
snpinfo <- fread(snpinfo.file)
cat('done\n')

## Filtering snp info
## Ive heard this is a horrible was to program things... but I do this sometimes -- feel free to edit
cat('Input SNPINFO N=',nrow(snpinfo),'\n')
snpinfo = eval(parse(text= paste0('subset(snpinfo,',snp.filter,')')))
cat('Output SNPINFO N=',nrow(snpinfo),'\n')


## phenotype 
phenotype.formula <- as.formula(phenotype.formula)
p <- read.csv(phenotype.file, header=TRUE, as.is=TRUE)
cat('Input pheno N=',nrow(p),'\n')
pheno <- reducePheno(p, phenotype.formula, pheno.id, genderCol)
cat('Output pheno N=',nrow(pheno),'\n')

## TBD: Report who was dropped between p and pheno


# For GDS files
f <- seqOpen(genotype.files)
sample.ids <- seqGetData(f, "sample.id")
pheno = pheno[row.names(pheno) %in% sample.ids,,drop=F]
#subset to phenotyped samples
seqSetFilter(f,sample.id = row.names(pheno))

# order pheno to the GDS subject order
sample.ids <- seqGetData(f, "sample.id")
pheno <- pheno[match(sample.ids,row.names(pheno)),,drop=F]


# get position list before any variant filters are applied
pos = seqGetData(f, "position")

## Load KINSHIP matrix
## Kinship doesn't contain all samples
kmatr = get(load(kinship.matrix))
pheno = pheno[row.names(pheno) %in% row.names(kmatr),,drop=F]
kmatr = kmatr[row.names(kmatr) %in% row.names(pheno),colnames(kmatr) %in% row.names(pheno)]
cat('Output pheno in Kinship N=',nrow(pheno),'\n')
kmatr = kmatr[match(row.names(pheno),row.names(kmatr)),match(row.names(pheno),colnames(kmatr))]
table(row.names(kmatr) %in% sample.ids ) # check that these match the genos
identical(sample.ids,row.names(pheno))
identical(row.names(kmatr),row.names(pheno))
cohort.file <- outputfilename

seqClose(f)
mydat <- data.frame(scanID = row.names(pheno), pheno = pheno)
mydat2 <- data.frame(sample.id = row.names(pheno), scanID = row.names(pheno), pheno = pheno)
scanAnnot <- ScanAnnotationDataFrame(mydat)
###################
## NULL MODEL
##################
tf <- terms(as.formula(phenotype.formula))
cn <- rownames(attr(tf, "factors"))
cat('start fit....\n')
kmatr_ns = as.matrix(kmatr)
nullmod <- fitNullMM(scanData = scanAnnot, outcome = paste0('pheno.',cn[1]), covars = paste0('pheno.',cn[2:length(cn)]), family = gaussian, covMatList = kmatr_ns)


if(!exists("null.fit")){
  cat('No Fit option for testtype:',TESTTYPE,' of family:',TESTFAMILY,'\n')
}
cat('done\n')


## For each aggregation unit - nemed 'gene' in code, but can be any start-stop region
## Work for each gene is parallelized over the cores
mcoptions <- list(preschedule=FALSE, set.seed=FALSE)
sm_obj <- 
foreach (cgene=genes, .combine='c',.inorder=FALSE, .options.multicore=mcoptions) %dopar% {

  ##############
  ## Apply variant filters to get a list of variants in each gene
  ###############
  
  gidx = which(kg$name == cgene)
  geneSNPinfo = subset(snpinfo, (POS > (kg[gidx,]$start - BUFFER) & POS < (kg[gidx,]$stop + BUFFER)))


  ## These nsmatch variables are too analysis specific and we may need to move to a array of variants per aggregation unit method
  ## NS=1 Matches nonsyn snps *only* if they match the named gene
  if(nsmatch == '1'){
    geneSNPinfo = subset(geneSNPinfo,  ! ( gene != kg[gidx,]$gene & sc_nonsynSplice == 1))
  ## NS=2 Matches nonsyn snps *only* if they match the named gene AND other variants *only* if they match 'gene2' column in annotation
  }else if(nsmatch == '2'){
    geneSNPinfo = subset(geneSNPinfo,  !( gene != kg[gidx,]$gene & sc_nonsynSplice == 1) || !( gene2 != kg[gidx,]$gene & isCorrDHS == 1) )
  }

  
  snp_idx = which(pos %in%  geneSNPinfo$POS)
  cat(cgene,'NEW \t',length(snp_idx),'\n')

  res = list()
  res[[cgene]]=list(generes=data.frame())
  if(length(snp_idx) > 0){
    
    ## extract genotypes
    f <- seqOpen(genotype.files)
    seqSetFilter(f,sample.id = row.names(pheno),verbose=FALSE)
    seqSetFilter(f,variant.id = snp_idx,verbose=FALSE)
    dose = altDosage(f, use.names=TRUE)
    colnames(dose) = paste0(seqGetData(f, "chromosome"),':',seqGetData(f, "position"))

    cat('+')
    cat(TESTTYPE)
    
    
    ## filter to maf
    dose = dose[,which(colMeans(dose,na.rm=T)/2 < top.maf),drop=F]

    ## flip to minor for burden tests
    if(TESTTYPE == 'B'){
      dose[,which(colMeans(dose,na.rm=T) > 1)] = 2-dose[,which(colMeans(dose,na.rm=T) > 1),drop=F]
    }

    NSNP = NCOL(dose)

    # Remove monomorphic varaints
    dose = dose[,which(colMeans(dose,na.rm=T)/2 >0),drop=F]
    NpolySNP = NCOL(dose)
    
   

    generes = data.frame('gene'=cgene,'cmaf'=sum(colSums(dose,na.rm=T)/2,na.rm=T)/NROW(dose),'npsnps'=NpolySNP,'nsnps'=NSNP,'p'=NA)
    cat(NpolySNP)
    if(NpolySNP > 0){
      cat('.')
      if(TESTTYPE == 'K'  & TESTFAMILY %in% c('C')){
        
        cat('-')
        class(f) = "SeqVarData"
        f@sampleData <- AnnotatedDataFrame(mydat2)
         # filter to MAC >=3
        seqSetFilter(f,variant.id = snp_idx[which(colSums(dose,na.rm=T) >= 3)],verbose=TRUE)
        
        system.time({generes <- assocTestMM(genoData = f, nullMMobj = nullmod, test = "Score")})
        generes$SNP = paste0(seqGetData(f, "chromosome"),':',seqGetData(f, "position"))
       
      }

      if(TESTTYPE == 'GENESIS_SKAT'  & TESTFAMILY %in% c('C')){
        
        cat('-')
        class(f) = "SeqVarData"
        f@sampleData <- AnnotatedDataFrame(mydat2)
         # filter to MAC >=3
        seqSetFilter(f,variant.id = snp_idx[which(colSums(dose,na.rm=T) >= 3)],verbose=TRUE)
        xlist = list()
        xlist[[1]] = data.frame('variant.id'=seqGetData(f,"variant.id"),
                                'allele.index'=rep(1,length(seqGetData(f,"variant.id")))
                                )
        
        system.time({
               test.out <- assocTestSeq(genoData = f, nullMMobj = nullmod,aggVarList=xlist)
                   })
        save(test.out,file=paste0(cgene,".RData"))
        #generes$SNP = paste0(seqGetData(f, "chromosome"),':',seqGetData(f, "position"))
       
      }


      
      ## MatrixFree SKAT - UNRELATED - TODO: update code for pedigree
      ## Fails back to regular SKAT if there is a warnin
      if(TESTTYPE == 'MF'  & TESTFAMILY %in% c('D','C')){
        skat.mf <- SKAT.matrixfree(dose, model=model)
        Qmf<-sum(skat.mf$tmult(resid(model))^2)/summary(model)$sigma^2
        generes$p = tryCatch({pQF(Qmf,skat.mf, neig=100, convolution.method="integration")},error = function(e) NA, warning=function(w) -9)
        if(p == -9 & NSNP < mmax){
          generes$p = tryCatch({SKAT(dose, null.fit)$p.value},error = function(e) NA)
        }
      }
      
      ## SKAT
      if(TESTTYPE == 'SK'  & TESTFAMILY %in% c('D','C') & NSNP < mmax){
        generes$p = tryCatch({SKAT(dose, null.fit)$p.value},error = function(e) NA)
      }
   
      
      ## BURDEN
      if(TESTTYPE == 'B' & TESTFAMILY %in% c('D','C')){

        pheno$burden = colSums(dose,na.rm=T)
        fit1<-update(null.fit, ~ . + burden,data=pheno)
        f = tryCatch(summary(fit1)$coef['burden',c('Estimate','Std. Error')], error = function(e) c(NA,NA))

        p = 2*pnorm(-abs(beta)/se)
        generes$beta = f[1]
        generes$se = f[2]
        generes$p = 2*pnorm(-abs(generes$beta)/generes$se)
        
      }
    }
    seqClose(f)
  }

  
  if(!exists("generes")){
    generes= data.frame('gene'=cgene,'cmaf'=NA,'npsnps'=NA,'nsnps'=NA,'p'=NA)
  }
  
  res[[cgene]]$generes = generes
  message(paste(cgene,"SNPS:",length(snp_idx),"pct:",round(which(genes == cgene)/length(genes),2),"\t",generes$gene[1],"finished."))
  res
}
  
cat('loop finished\n')

mgenes = names(sm_obj)
print(head(mgenes))
outres = list()
  
for(i in 1:length(mgenes)){
  outres[[i]] = sm_obj[[mgenes[i]]]$generes
}


#change to data.table - way faster
result = rbindlist(outres,fill=TRUE)


write.csv(result,file=outputfilename)
