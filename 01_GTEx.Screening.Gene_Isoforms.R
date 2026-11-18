###############################################################################
# 01_GTEx.Screening.Gene_Isoforms.R
#
# eQTL screening. For each tissue and each SNP of the VCF, test the association
# between the genotype (0/1/2) and the expression of CETP (gene or isoforms),
# adjusted for covariates and PEER factors.
#
# Outputs (per tissue), in Results/Expression/GTEx_ScreeningeQTL/<gene>/ :
#   GTEx.<tissue>.<gene>.MaxPEER.LM.ForMR[.SD].txt      beta / se / p-value of the SNP (input of the MR)
#   GTEx.<tissue>.<gene>.Fstat[.SD].CompareFullModel.txt F-statistic of the SNP (full vs. reduced model)
# and one file with all the tissues: GTEx.<gene>.AllTissues[.SD].txt
# With interactionSex=TRUE: model with a SNP x SEX interaction, outputs in <gene>/InteractionSex/
###############################################################################

source('helpers_GTEx.R')


#### 1. PARAMETERS ############################################################

allGenes=c('CETP_Isoforms','CETP') # 'CETP' (gene level) or 'CETP_Isoforms' (transcripts)
interactionSex=FALSE
standardDeviation=TRUE # standardize (scale) the expression before the regression

minSampleSNP=10    # minimal number of individuals with a genotype
minSample=30       # minimal number of individuals for a regression
minSampleTissue=50 # minimal number of samples in the expression file
minMAF=0.01        # minimal allele frequency of the SNP

# Ensembl IDs -> names used in the analysis
featureNames=c('ENST00000200676.7'='CETP201','ENST00000379780.6'='CETP202','ENST00000566128.1'='CETP203','ENST00000569082.1'='CETP204','ENSG00000087237.10'='CETP')
featuresToTest=c('CETP','CETP201','CETP202','CETP203')

dirData='../../Data/'
dirOut='../../Results/Expression/GTEx_ScreeningeQTL/'


#### 2. COMMON DATA ###########################################################

vcf=loadVCF(paste0(dirData,'vcf/GTEx.Hg38.ADCY9_CETP.MAF01.NoIndel.White.005.Modified.txt'))
covariates=readCovariates(dirData)
allCov=covariates$allCov
allCovSam=covariates$allCovSam


#### 3. SCREENING #############################################################

for(gene in allGenes){

	if(gene=='CETP_Isoforms'){
		dirExpression=paste0(dirData,'Expression/OnlyWhite/Normalized/CETP/Isoform/')
		listOrgan=list.files(dirExpression)
		suffixExpression='.Isoforms.Filter6_20.normalizedVoomWhite.CETP.txt'
	}else{
		dirExpression=paste0(dirData,'Expression/OnlyWhite/Normalized/',gene,'/')
		listOrgan=list.files(dirExpression,pattern='.txt')
		suffixExpression=paste0('.Gene.Filter6_20.OnlyWhite.normalizedLimmaVoom.',gene,'.txt')
	}

	# Output names (the SD tag is not used for the sex-interaction analysis)
	dirGene=paste0(dirOut,gene,'/')
	sdTag=if(standardDeviation) '.SD' else ''
	if(interactionSex){
		dirGene=paste0(dirGene,'InteractionSex/')
		sdTag=''
	}
	dir.create(dirGene,recursive=TRUE,showWarnings=FALSE)

	allResult=NULL

	for(organ in listOrgan){
		print(organ)

		# The tissue name is the file name before the first '.'
		tissue=sub('\\..*$','',organ)
		fileExpression=paste0(dirExpression,tissue,suffixExpression)
		tissue=gsub('\\)','',tissue)
		tissue=gsub('\\(','_',tissue)

		fileQTL=paste0(dirGene,'GTEx.',tissue,'.',gene,'.MaxPEER.LM.ForMR',sdTag,'.txt')
		fileFstat=paste0(dirGene,'GTEx.',tissue,'.',gene,'.Fstat',sdTag,'.CompareFullModel.txt')

		if(!file.exists(fileExpression)){
			print(paste('Expression file not found:',tissue))
			next
		}
		if(file.size(fileExpression)==0L){
			next
		}

		# Expression: first row of the transposed matrix = sample ID column + feature names
		expression=t(read.table(fileExpression))
		if(gene=='CETP_Isoforms' && sum(expression[1,]%in%names(featureNames)[1:4])==0){
			next
		}
		isKnown=expression[1,]%in%names(featureNames)
		expression[1,isKnown]=featureNames[expression[1,isKnown]]

		if(nrow(expression)<minSampleTissue){
			print(paste('Not enough samples:',tissue))
			next
		}

		# Add covariates (sample covariates, PEER factors, subject covariates)
		peer=read.table(paste0(dirData,'PEER/OnlyWhite/',tissue,'.OnlyWhite.PEER.txt'),header=TRUE)
		colnames(expression)=expression[1,]
		expression=expression[-1,]
		featureCols=colnames(expression)[-1]

		expression=merge(expression,allCovSam,by=1)
		expression=merge(expression,peer,by=1)
		expression[,1]=shortID(expression[,1])
		colnames(expression)[1]='ID'
		expression=merge(expression,allCov,by=1)

		if(interactionSex && length(unique(expression$SEX))==1){
			next
		}
		expression=dropConstantCols(expression,c('SMGEBTCHT','SMCENTER','SEX','RACE'))

		allPosition=NULL
		allFvalue=NULL

		if(file.exists(fileQTL)){
			# Already computed
			allPosition=read.table(fileQTL,header=TRUE,sep='\t')
			allPosition$Organ=tissue
		}else{
			rows=list()
			rowsF=list()

			for(pos in vcf$positions){ # for each SNP

				expressionSNP=expression
				expressionSNP$SNP=vcf$geno[match(expressionSNP$ID,vcf$ids),pos]
				expressionSNP=expressionSNP[!is.na(expressionSNP$SNP),]

				for(feature in featuresToTest[featuresToTest%in%colnames(expressionSNP)]){

					# Keep the tested feature (renamed GeneData) and drop the other ones
					d=expressionSNP[,!(colnames(expressionSNP)%in%c('ID',featureCols[featureCols!=feature]))]
					colnames(d)[colnames(d)==feature]='GeneData'
					d$GeneData=as.numeric(as.character(d$GeneData))

					if(nrow(d)<minSampleSNP){
						print(paste(organ,': low number of individuals'))
						next
					}

					freq=alleleFreq(d$SNP)
					if(freq$nGenotypes==1 || min(freq$other,freq$effect)<minMAF || nrow(d)<minSample){
						next
					}

					if(standardDeviation){
						d$GeneData=scale(d$GeneData)
					}

					info=data.frame(Tissus=tissue,Transcript=feature)
					alleles=data.frame(samplesize.exposure=nrow(d),other_allele.exposure=vcf$other[[pos]],effect_allele.exposure=vcf$effect[[pos]],oaf.exposure=freq$other,eaf.exposure=freq$effect,row.names=NULL)

					if(interactionSex){
						co=summary(lm(GeneData~SNP*SEX+.,data=d))$coefficients
						if(sum(rownames(co)%in%c('SNP','SNP:SEX'))<2){
							next
						}
						for(term in c('SNP','SNP:SEX')){
							rows[[length(rows)+1]]=cbind(info,data.frame(Terme=term,Position=as.numeric(pos),beta.exposure=co[term,1],se.exposure=co[term,2],pval.exposure=co[term,4]),alleles)
						}
					}else{
						fit=lm(GeneData~.,data=d)
						fitNoSNP=lm(GeneData~.-SNP,data=d)
						co=summary(fit)$coefficients
						if(!('SNP'%in%rownames(co))){
							next
						}
						rows[[length(rows)+1]]=cbind(info,data.frame(Position=as.numeric(pos),beta.exposure=co['SNP',1],se.exposure=co['SNP',2],pval.exposure=co['SNP',4]),alleles)
						rowsF[[length(rowsF)+1]]=data.frame(Tissus=tissue,Transcript=feature,Position=as.numeric(pos),'F stat'=anova(fit,fitNoSNP)[2,'F'],samplesize.exposure=nrow(d),check.names=FALSE)
					}
				}
			}

			if(length(rows)==0){
				next
			}
			allPosition=do.call(rbind,rows)
			allPosition$Organ=tissue
			if(length(rowsF)>0){
				allFvalue=do.call(rbind,rowsF)
			}
			write.table(allPosition,fileQTL,col.names=TRUE,quote=FALSE,row.names=FALSE,sep='\t')
			if(!is.null(allFvalue)){
				write.table(allFvalue,fileFstat,col.names=TRUE,quote=FALSE,row.names=FALSE,sep='\t')
			}
		}

		allResult=rbind(allResult,allPosition)
	}

	if(!is.null(allResult)){
		write.table(allResult,paste0(dirGene,'GTEx.',gene,'.AllTissues',sdTag,'.txt'),quote=FALSE,row.names=FALSE,col.names=TRUE,sep='\t')
	}
}
