library("devtools")
install_github("smgogarten/GENESIS")
install_github("smgogarten/SeqVarTools")

## Paths to files for CWS
phenotype.file <- "Commons_WES_FakePheno1_allCohort.csv"
phenotype.formula <- "pheno1~cohort+pheno2"
snpinfo.file <- "WGS_Freeze3_UCSCAnnovar_SNPlist_lipidAnnot_chr19.txt"
gene.file <- "AggUnit_chr19_50KBwindow.csv" # uploaded to annotation on DNAnexus	
genotype.files <- "CHGS_WGS_frz3_chr19.postqc.recode.gds"
snp.filter = 'isLipidAnnot_v1==1|sc_nonsynSplice==TRUE'
kinship.matrix <- "WGS_kinship_gwas_id.Rdata"
pheno.id <- "idno"
snpNames <- "SNP"
top.maf = 1
gene.filter = "TRUE"


library(GENESIS)
library(GWASTools)
library(svd)
library(survey)
library(SeqArray)
library(SeqVarTools)
library(gap)
library(Matrix)
library(plyr)
library(bdsmatrix)
library(coxme)
library(CompQuadForm)
library(parallel)
library(SKAT)
library(data.table)
library(doMC)

(num_cores <- detectCores(logical=TRUE))
registerDoMC(cores=num_cores-1)



## SNPINFO
cat('Reading snpinfo....')
snpinfo <- fread(snpinfo.file)
cat('done\n')
## Ive heard this is a horrible was to program things... but I do this sometimes -- feel free to edit
cat('Input SNPINFO N=',nrow(snpinfo),'\n')
snpinfo = eval(parse(text= paste0('subset(snpinfo,',snp.filter,')')))
cat('Output SNPINFO N=',nrow(snpinfo),'\n')



## get list of GENES
cat('reading gene file...')
kg = read.table(gene.file,as.is=T,sep=',',header=T) #jb
cat('GENE Filter',gene.filter,'\n')
kg = eval(parse(text= paste0('subset(kg,',gene.filter,')')))

cat(NROW(kg),'done\n')
genes <- kg$name

## PHENO
phenotype.formula <- as.formula(phenotype.formula)
pheno <- read.csv(phenotype.file, header=TRUE, as.is=TRUE) 


# For GDS files
f <- seqOpen(genotype.files)
sample.ids <- seqGetData(f, "sample.id")
pheno = pheno[pheno$idno %in% sample.ids,,drop=F]

# get position list before any variant filters are applied
pos = seqGetData(f, "position")

mmat = lm(phenotype.formula,data= pheno)
if (!is.null(mmat$na.action)) pheno <- pheno[-mmat$na.action,]


#subset to phenotyped samples
seqSetFilter(f,sample.id = pheno$idno)
sample.ids <- seqGetData(f, "sample.id")
pheno <- pheno[match(sample.ids,pheno$idno),,drop=F]


head(pheno)

# just doing it both ways for now, not sure why it isn't always the same
mydat <- data.frame(scanID = pheno$idno, pheno = pheno)
mydat2 <- data.frame(sample.id = pheno$idno, scanID = pheno$idno, pheno = pheno)


scanAnnot <- ScanAnnotationDataFrame(mydat)


# just grab a few SNPs for testing - will eventually loop over the whole set
gidx = which(kg$name == genes[100])	
geneSNPinfo = subset(snpinfo, (POS > (kg[gidx,]$start) & POS < (kg[gidx,]$stop)))
snp_idx = which(pos %in%  geneSNPinfo$POS)


## Load KINSHIP matrix
## Kinship doesn't contain all samples
kmatr = get(load(kinship.matrix))
pheno = pheno[pheno$idno %in% row.names(kmatr),,drop=F]
kmatr = kmatr[row.names(kmatr) %in% pheno$idno,colnames(kmatr) %in% pheno$idno]
cat('Output pheno in Kinship N=',nrow(pheno),'\n')

kmatr = kmatr[match(pheno$idno,row.names(kmatr)),match(pheno$idno,colnames(kmatr))]
table(row.names(kmatr) %in% sample.ids ) # check that these match the genos
identical(sample.ids,pheno$idno)
identical(row.names(kmatr),pheno$idno)




# need to be converted from sparse matrix
kmatr_ns = as.matrix(kmatr)

 
seqSetFilter(f,variant.id = snp_idx) 
seqSetFilter(f,sample.id = pheno$idno)

#phenotype.formula <- "pheno1~cohort+pheno2"
nullmod <- fitNullMM(scanData = scanAnnot, outcome = "pheno.pheno1", covars = c('pheno.cohort','pheno.pheno2'), family = gaussian, covMatList = kmatr_ns)

#genoData <- SeqVarData(f, sampleData=AnnotatedDataFrame(mydat2)) # didn't work for me - just hacking it together below
class(f) = "SeqVarData"
f@sampleData <- AnnotatedDataFrame(mydat2)

## This should be faster when new R build is ready
system.time({assoc <- assocTestMM(genoData = f, nullMMobj = nullmod, test = "Score")})
