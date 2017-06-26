# Read a gds file, filter, and output a vcf

LOCAL <- FALSE
if(LOCAL) {
  gds.file <- "/home/mbrown/DNAnexus/data/CHGS_WGS_frz3_chr22.postqc.recode.gds"
  start <- 16050150
  stop <- 16050950
} else {
  args <- commandArgs(trailing=TRUE)
  gds.file <- args[1]
  start <- args[2]
  stop <- args[3]
}

vcf.file <- "tmp.vcf"

library(SeqArray)
print(sessionInfo())

f <- seqOpen(gds.file)
chr <- seqGetData(f, "chromosome")
pos <- seqGetData(f, "position")
id <- paste(chr,pos,sep=":")
vars <- seqGetData(f, "variant.id")
wh <- vars[pos>=start & pos<=stop]

id.wh <- data.frame(ID=id[wh],stringsAsFactors=F)
print(head(id.wh))
write.table(id.wh,"list.txt",row.names=F,quote=F)

seqSetFilter(f, variant.id=wh)
seqGDS2VCF(f, vcf.file, info.var=character(), fmt.var=character())
seqClose(f)
