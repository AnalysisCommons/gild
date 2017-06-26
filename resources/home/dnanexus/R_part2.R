# Haven't got the gene based test to work yet
xlist = list()
xlist[[1]] = data.frame('variant.id='=seqGetData(f, "variant.id"),'allele.index'=rep(1,length(seqGetData(f, "variant.id")))) 
x =assocTestSeq(f, nullmod, aggVarList = xlist)
