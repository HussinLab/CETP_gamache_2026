###############################################################################
# helpers_GTEx.R
#
# Functions shared by the GTEx screening scripts
#   01_GTEx.Screening.Gene_Isoforms.R
#   02a_GTEx.Majiq.Screening.MR.sQTL_Fstat.R
# Load with: source('helpers_GTEx.R')   (the scripts are run from their own folder)
###############################################################################

library(reshape2) # colsplit


# Short GTEx subject ID: "GTEX-1117F" or "GTEX-1117F-0226-SM-5GZZ7" -> "1117F"
shortID=function(x){
	return(vapply(strsplit(as.character(x),'-',fixed=TRUE),function(p) p[2],character(1)))
}

# "0|0" -> 0, "0|1" or "1|0" -> 1, "1|1" -> 2, anything else (e.g. ".|.") -> NA
recodeGeno=function(x){
	y=rep(NA_real_,length(x))
	y[x %in% '0|0']=0
	y[x %in% c('0|1','1|0')]=1
	y[x %in% '1|1']=2
	return(y)
}

# Remove the columns of `cols` that have a single value
dropConstantCols=function(d,cols){
	for(v in cols){
		if(length(unique(d[[v]]))==1){
			d[[v]]=NULL
		}
	}
	return(d)
}

# Frequencies of the "other" (genotype 0) and "effect" (genotype 2) alleles,
# and number of different genotypes, for a vector of genotypes coded 0/1/2
alleleFreq=function(snp){
	g=table(factor(snp,levels=0:2))
	n=sum(g)
	return(list(other=unname((g[1]+g[2]/2)/n),effect=unname((g[3]+g[2]/2)/n),nGenotypes=sum(g>0)))
}

# Read the VCF. Its header line is read as a data row (V2 == 'POS'), so after the
# transposition the column 'POS' holds the sample IDs.
# Returns: geno      numeric matrix (samples x positions) coded 0/1/2
#          ids       short subject ID of each row of geno
#          positions SNP positions (column names of geno)
#          other, effect : named vectors (by position) of the REF and ALT alleles
loadVCF=function(file){
	raw=read.table(file)
	rownames(raw)=raw[,2]
	alleles=list(other=setNames(as.character(raw$V4),rownames(raw)),effect=setNames(as.character(raw$V5),rownames(raw)))

	geno=t(raw[,c(10:ncol(raw))])
	ids=shortID(geno[,'POS'])
	positions=setdiff(colnames(geno),'POS')
	geno=apply(geno[,positions,drop=FALSE],2,recodeGeno)
	rownames(geno)=ids

	return(list(geno=geno,ids=ids,positions=positions,other=alleles$other,effect=alleles$effect))
}

# Covariates (all except PEER factors)
#   nCols : number of columns of the principal components file to keep (IID + PCs); NULL = all
# Returns: allCov    one row per subject (age, sex, race, PCs, time from death)
#          allCovSam one row per sample (sequencing platform, collection center)
readCovariates=function(dirData,nCols=NULL){
	readCov=function(file,...){
		return(read.table(paste0(dirData,'Covariables/',file),header=TRUE,fill=TRUE,...))
	}

	pc=read.table(paste0(dirData,'Covariables/pcs.txt'),header=TRUE)
	pc=pc[,-1]
	if(!is.null(nCols)){
		pc=pc[,1:nCols]
	}
	colnames(pc)[1]='SUBJID'

	allCov=merge(readCov('AGE.txt'),readCov('SEX.txt'),by='SUBJID')
	allCov=merge(allCov,readCov('RACE.txt'),by='SUBJID')
	allCov=merge(allCov,pc,by='SUBJID')
	allCov=merge(allCov,readCov('TRISCHD.txt'),by='SUBJID') # time from death
	allCov$SUBJID=shortID(allCov$SUBJID)

	allCovSam=merge(readCov('SMGEBTCHT.txt',sep='\t'),readCov('SMCENTER.txt',sep='\t'),by='SAMPID') # sequencing platform, collection center

	return(list(allCov=allCov,allCovSam=allCovSam))
}

# Expression of `gene` adjusted for all covariates and PEER factors (residuals).
# Returns a data.frame (SUBJID, <gene>) that can be merged with the table of a tissue.
residualExpression=function(tissue,gene,dirData,covariates){
	expr=read.table(paste0(dirData,'Expression/OnlyWhite/Normalized/',gene,'/',tissue,'.Gene.Filter6_20.OnlyWhite.normalizedLimmaVoom.',gene,'.txt'))
	rownames(expr)=expr[,1]
	expr=expr[,-1]
	expr=t(expr)
	colnames(expr)[2]=gene # the second column holds the expression values

	peer=read.table(paste0(dirData,'PEER/OnlyWhite/',tissue,'.OnlyWhite.PEER.txt'),header=TRUE)
	peer=merge(expr,peer,by=1)
	peer=merge(peer,covariates$allCovSam,by=1)
	peer[,1]=shortID(peer[,1])

	cov=merge(covariates$allCov,peer,by=1)
	cov[[gene]]=as.numeric(as.character(cov[[gene]]))

	constant=sapply(cov,function(x) length(unique(x))==1)
	cov=cov[,!constant]

	cov[[gene]]=resid(lm(as.formula(paste(gene,'~ . - SUBJID')),data=cov))
	return(cov[,c('SUBJID',gene)])
}
