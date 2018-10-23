# code to make drug response prediction in TCGA

library(data.table)
library(parallel)

### set working directory 
#setwd("/cbcb/project2-scratch/jooslee/SL/revision/github8/")

### load the discretized mRNA expression data 
load("data/TCGA.drug.response/prob.TCGA.drug.response.RData")
### prob$genes  : gene symbols
### prob$samples: sample IDs
### prob$types  : cancer types
### prob$drug.patient: drug-patient table
### prob$drug.target : drug-target table
### prob$mRNA: gene expression level
### prob$mRNAq2: discretized gene expression level
### 0: lowly expressed (bottom 1/3 -quantile in each cancer type)
### 1: mediocre expressed (middle 1/3 -quantile in each cancer type)
### 2: highly expressed (top 1/3 -quantile in each cancer type)
### (see isle.r for further details)
### here we consider 9041 samples that is larger than 8749 samples used for the ISLE inference
### this is because for ISLE inferene we used the samples that have both mRNA and SCNA data available

### the outcome of ISLE pipeline with the [drug target] X [all genes] candidate SL pairs as input
load("data/TCGA.drug.response/sl.pairs.patient.drug.response.RData")
##########################################################################################
sl=sl.tot[,1:2]

n.targets=table(prob$drug.target$drugs)					#identify number of drug targets
drugs=names(n.targets)
drugs=drugs[n.targets<=3] 								#choose drugs with focused targets

n=m=rep(NA,length(drugs))
for (i in seq(length(drugs))) {
	idrug=which(prob$drug.patient$drug_name.adj %in% drugs[i])						#identify the samples treated by ith drug
	dpm=prob$drug.patient[idrug,]
	isample=which(prob$samples %in% prob$drug.patient$bcr_patient_barcode[idrug])
	n[i]=length(isample)
	i.dpm=match(prob$samples[isample],dpm$bcr_patient_barcode)						#pull up the response info of ith drug
	dpm=dpm[i.dpm,]
	ires=isample[which(dpm$response %in% c("Complete Response","Partial Response","Stable Disease"))]    # responder
	iirs=isample[which(dpm$response %in% c("Clinical Progressive Disease"))]							 # non-responder

	tb=table((prob$types[c(ires,iirs)]))
	sss2=names(tb)[tb>12]								# consider cancer types >12 samples 
	ires2=ires[which(prob$types[ires] %in% sss2)]		# responder with cancer types of sufficient sample size
	iirs2=iirs[which(prob$types[iirs] %in% sss2)]		# non-responder with cancer types of sufficient sample size
	m[i]=length(ires2)+length(iirs2)  					# number of samples in drug response annotation with cancer types of sufficient sample size
}
ix=order(m,decreasing=T)
drugs=drugs[ix]
m=m[ix]
n=n[ix]
drugs=drugs[m>12]										# choose drugs with sufficient sample size
FDR=0.2

res=sl.network=score2T=drugsT=groupT=NULL
drugs2=score2=group2=NULL
for (i in seq(length(drugs))) {
	idrug=which(prob$drug.patient$drug_name.adj %in% drugs[i])						#identify the samples treated by ith drug
	dpm=prob$drug.patient[idrug,]
	isample=which(prob$samples %in% prob$drug.patient$bcr_patient_barcode[idrug])	
	i.dpm=match(prob$samples[isample],dpm$bcr_patient_barcode)						#pull up the response info of ith drug
	dpm=dpm[i.dpm,]

	ires=isample[which(dpm$response %in% c("Complete Response","Partial Response","Stable Disease"))] #responders
	iirs=isample[which(dpm$response %in% c("Clinical Progressive Disease"))]						  #non-responders	
	tb=table((prob$types[c(ires,iirs)]))
	sss2=names(tb)[tb>12]							#filter out cancer types with marginal sample size (to remove noise)

	ires2=ires[which(prob$types[ires] %in% sss2)]	#responders for cancer types of sufficient sample size
	iirs2=iirs[which(prob$types[iirs] %in% sss2)]	#non-responders for cancer types of sufficient sample size

	idrug=which(prob$drug.target$drugs==drugs[i])	
	targetIDs=as.numeric(prob$drug.target$targetIDs[idrug])		#identify targets of ith drug
	sl.ess1=sl.tot[sl.tot[,1] %in% targetIDs,]		#identify sl interactions of the targets of ith drug
	flag=0

	### determine SL partners
	ix=which(p.adjust(sl.ess1[,3],"BH")<FDR)										# initial pool based on in vitro essentiality/SL screens
	if (length(ix)>1) {
		sl.ess1=sl.ess1[ix,]												
		iy=which(p.adjust(sl.ess1[,5],"BH")<FDR & p.adjust(sl.ess1[,6],"BH")<FDR)	# step I: underrepresented
		if (length(iy)>1) {
			sl.ess1=sl.ess1[iy,]													
			iz1=which(sl.ess1[,7]<0  & p.adjust(sl.ess1[,8],"BH")<FDR)			# step II: clinical relevance
			iz2=which(sl.ess1[,15]<0  & p.adjust(sl.ess1[,16],"BH")<FDR)
			iz1=iz1[which(sl.ess1[iz1,13]<0  & p.adjust(sl.ess1[iz1,14],"BH")<FDR)]
			iz2=iz2[which(sl.ess1[iz2,21]<0  & p.adjust(sl.ess1[iz2,22],"BH")<FDR)]	
			iz=union(iz1,iz2)
			if (length(iz)>1) {
				sl.ess1=sl.ess1[iz,]		
				ix=which(sl.ess1[,24]< .5)											# step III: phylogenetic
				if (length(ix)>0) flag=1											# flag: marker of SL partners
			}
		}
	}
	if (flag==1 & length(isample)>0){
		if (length(ix) >1 ) sl.network=rbind(sl.network,cbind(i,sl.ess1[ix,]))		# store SL pairs
		if (length(ix)==1 ) sl.network=rbind(sl.network,c(i,sl.ess1[ix,]))
		st=1

		if (length(ix)>1 & length(ires2)>1 & length(iirs2)>1){						# in case of multiple SL partners
			score.res=colSums(1*(prob$mRNAq2[sl.ess1[ix,2],ires2]==0),na.rm=T)		# cSL-score of responders
			score.irs=colSums(1*(prob$mRNAq2[sl.ess1[ix,2],iirs2]==0),na.rm=T)		# cSL-score of non-responders
			st=wilcox.test(score.res,score.irs,alternative="greater",exact=T)$p.value		# Wilcoxon test
		}
		if (length(ix)==1 & length(ires2)>1 & length(iirs2)>1){						# in case of 1 SL partner
			score.res=1*(prob$mRNAq2[sl.ess1[ix,2],ires2]==0)
			score.irs=1*(prob$mRNAq2[sl.ess1[ix,2],iirs2]==0)
			st=wilcox.test(score.res,score.irs,alternative="greater",exact=T)$p.value
		}
		
		drugs2=c(drugs2,rep(drugs[i],length(c(score.res,score.irs))))
		score2=c(score2,c(score.res,score.irs)/max(c(score.res,score.irs)))
		group2=c(group2,c(rep("responders",length(score.res)),rep("non-responders",length(score.irs))))
		res=rbind(res,c(i,length(ix),length(ires2),length(iirs2),st))	# collect outcome
	}else{
		res=rbind(res,c(i,0,0,0,1))
	}
}
	
colnames(res)=c("drug","n.partner","n.res","n.irs","pv")	
drugs=drugs[res[,5]!=1]																# consider drugs that have significant SL partners
res=res[res[,5]!=1,]
i.sig.drugs=which(p.adjust(res[,5],"BH")<0.2)						# FDR <=0.2
rownames(res)=drugs
res=cbind(res,p.adjust(res[,5],"BH"))
colnames(res)[6]="fdr"
