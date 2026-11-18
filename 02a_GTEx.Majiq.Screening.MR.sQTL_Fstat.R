###############################################################################
# 02a_GTEx.Majiq.Screening.MR.sQTL_Fstat.R
#
# sQTL screening of the CETP splicing junctions (MAJIQ PSI values).
# For each tissue, junction (Jct1up, Jct2up, Jct3up) and SNP of the VCF, test the
# association between the genotype (0/1/2) and the PSI, adjusted for covariates and
# PEER factors.
#
# Step 1 (only if it does not exist yet): build Results/Splicing/TableForLM.<tissue>.Center.txt
#        from the individual MAJIQ files (one row per sample and junction).
# Step 2: screening. Outputs, in Results/Splicing/Screening/ :
#   GTEx.<tissue>.Splicing.MaxPEER.LM.ForMR<tag>.txt            beta / se / p-value of the SNP (input of the MR)
#   GTEx.<tissue>.Splicing.Fstat.CompareFullModel<tag>.txt      F-statistic of the SNP (full vs. reduced model)
###############################################################################

source('helpers_GTEx.R')


#### 1. PARAMETERS ############################################################

addExpressionAsCovariate=FALSE # add the (residual) expression of CETP as a covariate
interactionSex=FALSE           # model psi ~ SNP * SEX
outlierPSILung=c('1LG7Y')      # individuals removed from the analysis; if not NULL, only the Lung is analysed

minSamples=30     # minimal number of individual MAJIQ files in a tissue
minSampleSNP=30   # minimal number of individuals for a regression
minMAF=0.01       # minimal allele frequency of the SNP
nPeer=16          # number of columns kept in the PEER file (sample ID + 15 factors)

junctions=c('Jct1up','Jct2up','Jct3up')

# MAJIQ junctions of CETP that are analysed (exon 9 region) and the splice sites that are kept
allJctPossible=c('gene:ENSG00000087237:t:56963010-56963124','gene:ENSG00000087237:s:56971992-56972083','gene:ENSG00000087237:t:56975101-56975151')
wantedPositions=c('56972083-56973331','56972083-56975101','56973510-56975101','56962097-56963010','56962296-56963010')

dirData='../../Data/'
dirSplicing='../../Results/Splicing/'
dirScreening=paste0(dirSplicing,'Screening/')
dir.create(dirScreening,recursive=TRUE,showWarnings=FALSE)

# Suffix of the output files
if(!is.null(outlierPSILung)){
	tag='.OutlierPSILung'
}else if(addExpressionAsCovariate && interactionSex){
	tag='.CETPasCov.IntSex'
}else if(interactionSex){
	tag='.IntSex'
}else if(addExpressionAsCovariate){
	tag='.CETPasCov'
}else{
	tag=''
}


#### 2. COMMON DATA ###########################################################

vcf=loadVCF(paste0(dirData,'vcf/GTEx.Hg38.ADCY9_CETP.MAF01.NoIndel.White.005.Modified.txt'))
covariates=readCovariates(dirData,nCols=6) # IID + 5 principal components
allCov=covariates$allCov
allCovSam=covariates$allCovSam


#### 3. FUNCTIONS #############################################################

# Keep only the wanted sub-junctions in the ';' separated strings (PSI, SD, positions) of a row
keepWantedJunctions=function(psi,sd,positions){
	keep=which(as.character(unlist(colsplit(positions,';',paste0('V',1:5))))%in%wantedPositions)
	subset=function(string){
		return(paste(as.character(unlist(colsplit(string,';',paste0('V',1:5))[keep])),collapse=';'))
	}
	return(c(subset(psi),subset(sd),subset(positions)))
}

# Read the ASpli file of one individual: one row per junction of allJctPossible
readMajiqSample=function(file,id){
	oneSample=read.table(file,sep='\t',header=TRUE)
	oneSample=oneSample[oneSample[,2]%in%allJctPossible,]
	if(nrow(oneSample)==0){
		return(NULL)
	}
	jctName=data.frame(JctName=oneSample[,2])
	oneSample=oneSample[,c(4,5,11)] # PSI, PSI standard deviation, positions of the sub-junctions
	oneSample[]=lapply(oneSample,as.character)

	for(row in seq_len(nrow(oneSample))){
		oneSample[row,]=keepWantedJunctions(oneSample[row,1],oneSample[row,2],oneSample[row,3])
	}

	psiValue=colsplit(oneSample[,1],';',c('Jct1PSI','Jct2PSI','Jct3PSI'))
	psiSD=colsplit(oneSample[,2],';',c('Jct1SD','Jct2SD','Jct3SD'))
	pos=colsplit(oneSample[,3],';',c('Jct1Pos','Jct2Pos','Jct3Pos'))

	oneSample=cbind(jctName,psiValue,psiSD,pos)
	oneSample$ID=id
	return(oneSample)
}

# Length of a junction given as "start-end"
junctionLength=function(position){
	ab=as.numeric(strsplit(position,'-')[[1]])
	return(ab[2]-ab[1])
}

# Table of one tissue for the regressions (one row per sample and junction: psi, covariates, PEER)
buildTableForLM=function(tissue,dirIndividual,files){

	allSamples=NULL
	for(file in files){
		id=sub('-.*$','',file) # ID = file name before the first '-'
		allSamples=rbind(allSamples,readMajiqSample(paste0(dirIndividual,file),id))
	}
	if(is.null(allSamples)){
		return(NULL)
	}

	allJct=NULL
	for(o in seq_along(allJctPossible)){ # o = junction number
		onlyThisJct=allSamples[allSamples$JctName==allJctPossible[o],]
		if(nrow(onlyThisJct)==0){
			next
		}

		# Order of the "up" and "down" junctions depends on their length
		if(junctionLength(onlyThisJct[1,'Jct1Pos'])>junctionLength(onlyThisJct[1,'Jct2Pos'])){
			p=c('up','down')
		}else{
			p=c('down','up')
		}

		# One table per sub-junction
		subJunction=function(k,name,type){
			d=onlyThisJct[,grep(paste0('Jct',k),colnames(onlyThisJct))]
			colnames(d)=c('psi','sd.psi','Position')
			d$Jct=name
			d$JctT=type
			d$ID=onlyThisJct$ID
			return(d)
		}
		jct=rbind(subJunction(1,paste0('Jct',o,p[1]),p[1]),
		          subJunction(2,paste0('Jct',o,p[2]),p[2]),
		          subJunction(3,paste0('Jct',o,'More',p[2]),paste0('More',p[2])))
		jct=jct[!is.na(jct$psi),]
		allJct=rbind(allJct,jct)
	}

	# Add PEER factors, sample covariates and subject covariates
	peer=read.table(paste0(dirData,'PEER/OnlyWhite/',tissue,'.OnlyWhite.PEER.txt'),header=TRUE)
	peer=peer[,1:nPeer]
	peer=merge(peer,allCovSam,by=1,sort=FALSE)
	peer[,1]=shortID(peer[,1])

	forLM=merge(allJct,allCov,by.x='ID',by.y='SUBJID')
	forLM=merge(forLM,peer,by=1)
	forLM=forLM[,!(colnames(forLM)%in%c('sd.psi','IID','Position','Diff','PC6','PC7','PC8','PC9','PC10'))]
	forLM=dropConstantCols(forLM,c('SMGEBTCHT','SMCENTER'))
	return(forLM)
}


#### 4. SCREENING #############################################################

listOrgan=list.files(dirSplicing)
listOrgan=listOrgan[!grepl('txt|pdf',listOrgan)]
if(!is.null(outlierPSILung)){
	listOrgan='Lung'
}

for(tissue in listOrgan){
	print(tissue)

	fileQTL=paste0(dirScreening,'GTEx.',tissue,'.Splicing.MaxPEER.LM.ForMR',tag,'.txt')
	fileFstat=paste0(dirScreening,'GTEx.',tissue,'.Splicing.Fstat.CompareFullModel',tag,'.txt')

	# ASpli files of the individuals ('Individuel' is the name of the folder)
	dirIndividual=paste0(dirSplicing,tissue,'/psi/Individuel/')
	files=list.files(dirIndividual)
	files=files[grepl('tsv',files)]
	if(length(files)<minSamples){
		print('Not enough ASpli files')
		next
	}

	# Table for the regressions
	fileTable=paste0(dirSplicing,'TableForLM.',tissue,'.Center.txt')
	if(file.exists(fileTable)){
		forLM=read.table(fileTable,header=TRUE,sep='\t')
	}else{
		forLM=buildTableForLM(tissue,dirIndividual,files)
		if(is.null(forLM)){
			print('No junction found')
			next
		}
		write.table(forLM,fileTable,col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')
	}

	if(!is.null(outlierPSILung)){
		forLM=forLM[!(forLM$ID%in%outlierPSILung),]
	}

	if(addExpressionAsCovariate){
		forLM=merge(forLM,residualExpression(tissue,'CETP',dirData,covariates),by=1,sort=FALSE)
		forLM$CETP=as.numeric(as.character(forLM$CETP))
	}

	rows=list()
	rowsF=list()

	for(pos in vcf$positions){ # for each SNP

		expression=forLM
		expression$SNP=vcf$geno[match(expression$ID,vcf$ids),pos]
		expression=expression[!is.na(expression$SNP),]
		expression$psi=as.numeric(as.character(expression$psi))

		for(jct in junctions){

			d=expression[expression$Jct==jct,]
			d=d[,!(colnames(d)%in%c('ID','Jct','JctT'))]
			d=dropConstantCols(d,c('SMGEBTCHT','SMCENTER','SEX','RACE'))

			freq=alleleFreq(d$SNP)
			if(freq$nGenotypes==1 || min(freq$other,freq$effect)<minMAF || nrow(d)<minSampleSNP){
				next
			}

			alleles=data.frame(samplesize.exposure=nrow(d),other_allele.exposure=vcf$other[[pos]],effect_allele.exposure=vcf$effect[[pos]],oaf.exposure=freq$other,eaf.exposure=freq$effect)

			if(interactionSex){
				co=summary(lm(psi~SNP*SEX+.,data=d))$coefficients
				if(sum(rownames(co)=='SNP:SEX')==0){
					next
				}
				rows[[length(rows)+1]]=cbind(data.frame(Tissus=tissue,Junction=jct,Position=as.numeric(pos),
					beta.exposure=co['SNP',1],se.exposure=co['SNP',2],pval.exposure=co['SNP',4],
					beta.sex.exposure=co['SNP:SEX',1],se.sex.exposure=co['SNP:SEX',2],pval.sex.exposure=co['SNP:SEX',4]),alleles)
			}else{
				fit=lm(psi~.,data=d)
				fitNoSNP=lm(psi~.-SNP,data=d)
				co=summary(fit)$coefficients
				if(sum(rownames(co)=='SNP')==0){
					next
				}
				rows[[length(rows)+1]]=cbind(data.frame(Tissus=tissue,Junction=jct,Position=as.numeric(pos),
					beta.exposure=co['SNP',1],se.exposure=co['SNP',2],pval.exposure=co['SNP',4]),alleles)
				rowsF[[length(rowsF)+1]]=data.frame(Tissus=tissue,Transcript=jct,Position=as.numeric(pos),
					'F stat'=anova(fit,fitNoSNP)[2,'F'],samplesize.exposure=nrow(d),check.names=FALSE)
			}
		}
	}

	if(length(rows)==0){
		next
	}
	write.table(do.call(rbind,rows),fileQTL,col.names=TRUE,quote=FALSE,row.names=FALSE,sep='\t')
	if(!interactionSex){
		write.table(do.call(rbind,rowsF),fileFstat,col.names=TRUE,quote=FALSE,row.names=FALSE,sep='\t')
	}
}
