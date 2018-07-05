##### code to predict clinical drug response in
##### 1. Hatzis et al (2012) - Taxane+anthracycline treatment in 508 breast cancer patients
##### 2. Bryer et al (2013) - Erlotinib treatment in 25 lung cancer patients
##### 3. Patch et al (2015) - Taxane (+cisplatin) treatment in 80 ovarian cancer patients

library(data.table)
library(Rcpp)
library(parallel)
library(survival)

##### set up the working directory
#setwd("/cbcb/project2-scratch/jooslee/SL/revision/github3")

##### define basic functions
qnorm.array <- function(mat)     
{
	mat.back = mat 
	mat = mat[!is.na(mat)]
    mat = rank(mat, ties.method = "average");
    mat = qnorm(mat / (length(mat)+1));
    mat.back[!is.na(mat.back)] = mat 
    mat.back
}
rank.array <- function(mat)     
{
	mat.back = mat 
	mat = mat[!is.na(mat)]
    mat = rank(mat, ties.method = "average")/length(mat);
    mat.back[!is.na(mat.back)] = mat 
    mat.back
}

##### load source file
sourceCpp("src/logRankGene.cpp")

##########################################################################################
##### 1. Hatzis et al (2012) -Taxane+anthracycline treatment in 508 breast cancer patients
##########################################################################################
cox.sl.hatzis = function(surv.dt1,age1,score1,types1)
{
	dt1 = data.frame(cbind(surv.dt1, score1, types1))
	cox.out = coxph(Surv(time,status) ~ score1 + age1+strata(types1), data=dt1)
	aa  = summary(cox.out)
	return(aa$coefficients["score1",])
} 

## load ISLE results for drug targets
load("data/TCGA.drug.response/sl.pairs.patient.drug.response.RData")

## load drug prob
load("data/GSE25055.RData")
targetIDs=prob$targetIDs
ix=which(sl.tot[,1] %in% targetIDs)
sl.tot=sl.tot[ix,]
sl=sl.tot[,1:2]
sl.ess1=sl.tot
FDR=0.2
library(nestedRanksTest)
flag=0
ix=which(p.adjust(sl.ess1[,3],"BH")<FDR)
if (length(ix)>1) {
	sl.ess1=sl.ess1[ix,]
	iy=which(p.adjust(sl.ess1[,5],"BH")<FDR & p.adjust(sl.ess1[,6],"BH")<FDR)
	if (length(iy)>1) {
		sl.ess1=sl.ess1[iy,]
		iz1=which(sl.ess1[,7]<0  & p.adjust(sl.ess1[,8],"BH")<FDR)
		iz2=which(sl.ess1[,15]<0  & p.adjust(sl.ess1[,16],"BH")<FDR)
		iz1=iz1[which(sl.ess1[iz1,13]<0  & p.adjust(sl.ess1[iz1,14],"BH")<FDR)]
		iz2=iz2[which(sl.ess1[iz2,21]<0  & p.adjust(sl.ess1[iz2,22],"BH")<FDR)]	
		iz=union(iz1,iz2)
		if (length(iz)>1) {
			sl.ess1=sl.ess1[iz,]		
			ix=which(sl.ess1[,24]< .5)						
			if (length(ix)>0) 
			{
				flag=1
				sl.ess1=sl.ess1[ix,]
			}
		}
	}
}
	
g1=sl.ess1[,2]      # cSL-partners
types1=prob$subtypes;surv.dt1=prob$surv.dt
age1=qnorm.array(as.numeric(prob$age))
score1=qnorm.array(as.numeric(colSums((prob$mRNAq[g1,]==0)*1,na.rm=T))) # cSL-score

## cox regression based on cSL-score controlling for age and breast cancer subtypes
cox.out=cox.sl.hatzis(surv.dt1,age1,score1,types1)[c(1,5)]
cox.beta=exp(-(cox.out[1]))
cox.p=cox.out[2]
	
## logrank test comparing the survival of the patients with high cSL-score vs low cSL-score	
ix1=(score1>median(score1) & !is.na(prob$survival[,1]))
ix2=(score1<=median(score1) & !is.na(prob$survival[,1]))
times1=prob$survival[ix1,]
times2=prob$survival[ix2,]
out=logRank(times1,times2)											
logrank.p=out[1]
logrank.deltaAUC=out[7]-out[8]

## compare the cSL-score of responders vs nonresponders controlling for breast cancer subtypes
dt=data.table(score=score1,subtypes=prob$subtypes,response=prob$response)
dt=dt[!is.na(dt$response),]
dt1=dt
dt1$response1=ifelse(dt1$response=="pCR",1,0)		
ranksum.p=nestedRanksTest(score ~ response1 |subtypes, data=dt1, alternative="greater")$p.value	

## control:ncSL partners
ig=sl.tot[order(sl.tot[,3])[seq(nrow(sl.ess1))],2]
score.ncSL=colSums(prob$mRNAq[ig,]==0,na.rm=T)
ix1=(score.ncSL>median(score.ncSL) & !is.na(prob$survival[,1]))
ix2=(score.ncSL<=median(score.ncSL) & !is.na(prob$survival[,1]))
times1=prob$survival[ix1,]
times2=prob$survival[ix2,]
out=logRank(times1,times2)											
logrank.p.ncSL=out[1] 
logrank.deltaAUC.ncSL=out[7]-out[8]
#logrank.p.cnSL=0.37

##### control:DAISY SL partners
load("data/prob.TCGA.RData")
sl=sl.tot[,1:2]

mRNA=prob$mRNA;mRNA.rank2=prob$mRNA.rank2;scna=prob$scna
dmol1=do.call(rbind,mclapply(seq(nrow(sl)),function(X){
	gene1=sl[X,1];gene2=sl[X,2]
	q1=quantile(mRNA[gene1,],0.1,na.rm=T)
	dn=which(scna[gene1,]< -0.2 & mRNA[gene1,]< q1)
	up=as.numeric(which(scna[gene1,]>0))
	st.mRNA1=st.scna1=st.mRNA2=st.scna2=NA
	if (sum(!is.na(dn))>100 & sum(!is.na(up))>100 & sum(!is.na(mRNA[gene2,]))>100 & sum(!is.na(scna[gene2,]))>100 ){
		st.mRNA1=wilcox.test(mRNA[gene2,dn],mRNA[gene2,up],alternative="greater")$p.value
		st.scna1=wilcox.test(scna[gene2,dn],scna[gene2,up],alternative="greater")$p.value
	}
	q1=quantile(mRNA[gene2,],0.1,na.rm=T)
	dn=which(scna[gene2,]< -0.2 & mRNA[gene1,]< q1)
	up=as.numeric(which(scna[gene2,]>0))
	if (sum(!is.na(dn))>100 & sum(!is.na(up))>100 & sum(!is.na(mRNA[gene1,]))>100 & sum(!is.na(scna[gene1,]))>100 ){
		st.mRNA2=wilcox.test(mRNA[gene1,dn],mRNA[gene1,up],alternative="greater")$p.value
		st.scna2=wilcox.test(scna[gene1,dn],scna[gene1,up],alternative="greater")$p.value
	}
	print(X)
	return(c(st.mRNA1,st.scna1,st.mRNA2,st.scna2))
},mc.cores=20))

cr1=do.call(rbind,mclapply(seq(nrow(sl)),function(X){
	gene1=sl[X,1];gene2=sl[X,2]
	cr=cor.test(prob$mRNA[gene1,],prob$mRNA[gene2,],method="spearman")$estimate
},mc.cores=20))

p.mRNA=apply(dmol1[,c(1,3)],1,max,na.rm=T)
p.scna=apply(dmol1[,c(2,4)],1,max,na.rm=T)

ig=which(p.adjust(sl.tot[,3],"BH")<0.05 & p.adjust(p.mRNA,"BH")<0.05 & p.adjust(p.scna,"BH")<0.05 & cr1>0.3) 
ig=unique(match(unique(prob$genes[sl.tot[ig,2]]),prob$genes))

score.daisy=colSums(prob$mRNAq[ig,]==0,na.rm=T)
ix1=(score.daisy>median(score.daisy) & !is.na(prob$survival[,1]))
ix2=(score.daisy<=median(score.daisy) & !is.na(prob$survival[,1]))
times1=prob$survival[ix1,]
times2=prob$survival[ix2,]
out=logRank(times1,times2)											
logrank.p.daisy=out[1] 
logrank.deltaAUC.daisy=out[7]-out[8]
#logrank.p.daisy=0.046 (deltaAUC=-0.03) -> significance in the opposite direction

##########################################################################################
##### 2. Bryer et al (2013) - Erlotinib treatment in 25 lung cancer patients #############
##########################################################################################
cox.sl.erlotinib = function(surv.dt1,score1,types1)
{
	dt1 = data.frame(cbind(surv.dt1, score1, types1))
	cox.out = coxph(Surv(time,status) ~ score1 + strata(types1), data=dt1)
	aa  = summary(cox.out)
	aa$coefficients["score1",]
}

##### upload drug-cSL network
load("data/TCGA.drug.response/sl.pairs.patient.drug.response.RData")
sl=sl.tot[,1:2]
ii=which(sl.tot[,1] %in% which(prob$genes=="EGFR"))	# select candidate SL-partners of EGFR
sl.ess1=sl.tot=sl.tot[ii,]

##### calculate cSL-pair-score
ess.r=1-rank.array(sl.ess1[,3])									#r_initial
mol.r=1-rank.array(apply(sl.ess1[,5:6],1,max,na.rm=T))			#r_step I
cln.r=1-rank.array(												
	apply(cbind(-apply(cbind(sl.tot[, 7],sl.tot[,13]),1,min,na.rm=T),
	-apply(cbind(sl.tot[,15],sl.tot[,21]),1,min,na.rm=T)),1,max,na.rm=T))
phy.r=1-sl.ess1[,24]											#r_step III

ps=apply(cbind(ess.r,mol.r,cln.r,phy.r),1,mean,na.rm=T)			#cSL-pair-score
ii=order(ps,decreasing=T)
sl.cox=sl.ess1=sl.ess1[ii[1:76],]										#select top 76-genes
								
load("data/GSE33072.RData")	# erlotinib treated cohort
types1=prob$kras;	surv.dt1=prob$surv.dt
mRNA1=prob$mRNA;	mRNAq1=prob$mRNAq;	samples1=prob$sample
survival1=prob$survival;	targetID=prob$targetIDs

load("data/GSE33072.sorafenib.RData") # sorafenib treated cohort
mRNAq2=prob$mRNAq;survival2=prob$survival

g1=sl.cox[,2]														# identify SL-partners
score1=qnorm.array(as.numeric(colSums((mRNAq1[g1,]==0)*1,na.rm=T))) # calculate cSL-score for erlotinib treated cohort
score2=qnorm.array(as.numeric(colSums((mRNAq2[g1,]==0)*1,na.rm=T))) # calculate cSL-score for sorafenib treated cohort

cox.out=cox.sl.erlotinib(surv.dt1,score1,types1)[c(1,5)]	# perform Cox regression in erlotinib cohort
cox.beta=exp(abs(cox.out[1]))
cox.p=cox.out[2]

ix1=(score1>median(score1) & !is.na(survival1[,1]))
ix2=(score1<=median(score1) & !is.na(survival1[,1]))
times1=survival1[ix1,]
times2=survival1[ix2,]
out1=logRank(times1,times2)											# perform logrank test based on cSL-score in erlotinib cohort
logrank.p=out1[1]
logrank.deltaAUC=out1[7]-out1[8]

ix1=(score2>median(score2) & !is.na(survival2[,1]))
ix2=(score2<=median(score2) & !is.na(survival2[,1]))
times1=survival2[ix1,]
times2=survival2[ix2,]
out2=logRank(times1,times2)											# perform logrank test based on cSL-score in sorafenib cohort
logrank.p.sorafenib=out2[1]
logrank.deltaAUC.sorafenib=out2[7]-out2[8]

cr=cor.test(score1,survival1[,1],method="spearman")					# correlation between cSL-score and DFS
cr.kras=cor.test(score1[types1!="Mutant"],
				 survival1[types1!="Mutant",1],method="spearman")	# correlation between cSL-score and DFS in KRAS-WT population

### 8-week treatment
load("data/GSE33072.RData")
eight=8*7/30
surv=prob$survival[,1]
t.test(score1[surv>=eight],score1[surv<eight],alternative="greater")

## control:ncSL partners
ig=sl.tot[order(sl.tot[,3])[seq(nrow(sl.ess1))],2]
score.ncSL=colSums(mRNAq1[ig,]==0,na.rm=T)
ix1=(score.ncSL>median(score.ncSL))
ix2=(score.ncSL<=median(score.ncSL))
times1=survival1[ix1,]
times2=survival1[ix2,]
out=logRank(times1,times2)											
logrank.p.ncSL=out[1] 
logrank.deltaAUC.ncSL=out[7]-out[8]
cr.ncSL=cor.test(score.ncSL,survival1[,1],method="spearman")
#logrank.p.cnSL=0.56

##### control:DAISY SL partners
load("data/prob.TCGA.RData")
sl=sl.tot[,1:2]

mRNA=prob$mRNA;mRNA.rank2=prob$mRNA.rank2;scna=prob$scna
dmol1=do.call(rbind,mclapply(seq(nrow(sl)),function(X){
	gene1=sl[X,1];gene2=sl[X,2]
	q1=quantile(mRNA[gene1,],0.1,na.rm=T)
	dn=which(scna[gene1,]< -0.2 & mRNA[gene1,]< q1)
	up=as.numeric(which(scna[gene1,]>0))
	st.mRNA1=st.scna1=st.mRNA2=st.scna2=NA
	if (sum(!is.na(dn))>100 & sum(!is.na(up))>100 & sum(!is.na(mRNA[gene2,]))>100 & sum(!is.na(scna[gene2,]))>100 ){
		st.mRNA1=wilcox.test(mRNA[gene2,dn],mRNA[gene2,up],alternative="greater")$p.value
		st.scna1=wilcox.test(scna[gene2,dn],scna[gene2,up],alternative="greater")$p.value
	}
	q1=quantile(mRNA[gene2,],0.1,na.rm=T)
	dn=which(scna[gene2,]< -0.2 & mRNA[gene1,]< q1)
	up=as.numeric(which(scna[gene2,]>0))
	if (sum(!is.na(dn))>100 & sum(!is.na(up))>100 & sum(!is.na(mRNA[gene1,]))>100 & sum(!is.na(scna[gene1,]))>100 ){
		st.mRNA2=wilcox.test(mRNA[gene1,dn],mRNA[gene1,up],alternative="greater")$p.value
		st.scna2=wilcox.test(scna[gene1,dn],scna[gene1,up],alternative="greater")$p.value
	}
	print(X)
	return(c(st.mRNA1,st.scna1,st.mRNA2,st.scna2))
},mc.cores=20))

cr1=do.call(rbind,mclapply(seq(nrow(sl)),function(X){
	gene1=sl[X,1];gene2=sl[X,2]
	cr=cor.test(prob$mRNA[gene1,],prob$mRNA[gene2,],method="spearman")$estimate
},mc.cores=20))

p.mRNA=apply(dmol1[,c(1,3)],1,max,na.rm=T)
p.scna=apply(dmol1[,c(2,4)],1,max,na.rm=T)

ig=which(p.adjust(sl.tot[,3],"BH")<0.05 & p.adjust(p.mRNA,"BH")<0.05 & p.adjust(p.scna,"BH")<0.05 & cr1>0.3) 
ig=unique(match(unique(prob$genes[sl.tot[ig,2]]),prob$genes))

score.daisy=colSums(mRNAq1[ig,]==0,na.rm=T)
ix1=(score.daisy>median(score.daisy) & !is.na(prob$survival[,1]))
ix2=(score.daisy<=median(score.daisy) & !is.na(prob$survival[,1]))
times1=prob$survival[ix1,]
times2=prob$survival[ix2,]
out=logRank(times1,times2)											
logrank.p.daisy=out[1] 
logrank.deltaAUC.daisy=out[7]-out[8]
#logrank.p.daisy=0.67
cr.daisy=cor.test(score.daisy,survival1[,1],method="spearman")


##########################################################################################
##### 3. Patch et al (2015) - Taxane (+cisplatin) treatment in 80 ovarian cancer patients
##########################################################################################
##### load ovarian cancer data
load("data/prob.ICGC.ovarian.RData")
prob.ov=prob

##### load pancancer data just to match gene indices
load("data/TCGA.drug.response/prob.TCGA.drug.response.RData")

##### load drug-cSL network
load("data/TCGA.drug.response/sl.pairs.patient.drug.response.RData")
targets=c("TUBB1","BCL2","MAPT","MAP2","MAP4") #Taxane targets
ix=which(sl.tot[,1] %in% match(targets,prob$genes))
sl.tot=sl.tot[ix,]

##### gene ID mapping to ovarian cancer dataset
sl=cbind(prob$genes[sl.tot[,1]],prob$genes[sl.tot[,2]])
sl=cbind(match(sl[,1],prob.ov$genes),match(sl[,2],prob.ov$genes))
sl.tot=sl.tot[!is.na(sl[,2]),]
sl=sl[!is.na(sl[,2]),]
sl.tot[,1:2]=sl						# change gene IDs based on the ovarian cancer dataset

##### define Taxane's SL partners
FDR=0.2
sl.ess1=sl.tot
ix=which(p.adjust(sl.ess1[,3],"BH")<FDR)
if (length(ix)>1) {
	sl.ess1=sl.ess1[ix,]
	iy=which(p.adjust(sl.ess1[,5],"BH")<FDR & p.adjust(sl.ess1[,6],"BH")<FDR)
	if (length(iy)>1) {
		sl.ess1=sl.ess1[iy,]
		iz1=which(sl.ess1[,7]<0  & p.adjust(sl.ess1[,8],"BH")<FDR)
		iz2=which(sl.ess1[,15]<0  & p.adjust(sl.ess1[,16],"BH")<FDR)
		iz1=iz1[which(sl.ess1[iz1,13]<0  & p.adjust(sl.ess1[iz1,14],"BH")<FDR)]
		iz2=iz2[which(sl.ess1[iz2,21]<0  & p.adjust(sl.ess1[iz2,22],"BH")<FDR)]	
		iz=union(iz1,iz2)
		if (length(iz)>1) {
			sl.ess1=sl.ess1[iz,]		
			ix=which(sl.ess1[,24]< .5)
			if (length(ix)>0) flag=1
		}
	}
}
n.partners=nrow(sl.ess1)

##### obtain cSL-score and responder/nonresponder information
score1=colSums(prob.ov$mRNAq[sl.ess1[,2],]==0,na.rm=T)								#cSL-score
ires= which(prob.ov$tp=="primary" & prob.ov$rp=="sensitive" )						#responders
iirs= which(prob.ov$tp=="primary" & prob.ov$rp %in% c("resistant","refractory"))	#nonresponders
ires=ires[-c(9,33,34)]  # remove overlaping samples from responders
iirs=iirs[-4]			# remove overlaping samples from nonresponders

##### cSL-score distinguishes responders vs nonresponders
wilcox.test(score1[ires],score1[iirs],alternative="greater")						

##### control:randomly selected SL partners
st.random=1
for (i in seq(1000)){
score.random=colSums(prob.ov$mRNAq[sample(which(!is.na(prob.ov$mRNAq[,1])),nrow(sl.ess1)),]==0,na.rm=T)
st.random[i]=wilcox.test(score.random[ires],score.random[iirs],alternative="greater")$p.value}
sum(st.random<0.009062)/1000
# empirical P < 0.043

##### control:SL partners of other drugs
load("data/TCGA.drug.response/sl.pairs.patient.drug.response.RData")
targets=c("TUBB1","BCL2","MAPT","MAP2","MAP4")
ix=which(sl.tot[,1] %in% match(targets,prob$genes))
sl.tot=sl.tot[-ix,]
sl.ess1=sl.tot
ix=which(p.adjust(sl.ess1[,3],"BH")<FDR)
if (length(ix)>1) {
	sl.ess1=sl.ess1[ix,]
	iy=which(p.adjust(sl.ess1[,5],"BH")<FDR & p.adjust(sl.ess1[,6],"BH")<FDR)
	if (length(iy)>1) {
		sl.ess1=sl.ess1[iy,]
		iz1=which(sl.ess1[,7]<0  & p.adjust(sl.ess1[,8],"BH")<FDR)
		iz2=which(sl.ess1[,15]<0  & p.adjust(sl.ess1[,16],"BH")<FDR)
		iz1=iz1[which(sl.ess1[iz1,11]<0  & p.adjust(sl.ess1[iz1,12],"BH")<FDR)]
		iz2=iz2[which(sl.ess1[iz2,19]<0  & p.adjust(sl.ess1[iz2,20],"BH")<FDR)]	
		iz=union(iz1,iz2)
		if (length(iz)>1) {
			sl.ess1=sl.ess1[iz,]		
			ix=which(sl.ess1[,24]< .5)
			if (length(ix)>0){
				 flag=1;sl.ess1=sl.ess1[ix,]
			}
		}
	}
}
sl=cbind(prob$genes[sl.ess1[,1]],prob$genes[sl.ess1[,2]])
sl1=cbind(match(sl[,1],prob.ov$genes),match(sl[,2],prob.ov$genes))

st.random=1
for (i in seq(1000)){
score.random=colSums(prob.ov$mRNAq[sample(sl1[,2],n.partners),]==0,na.rm=T)
st.random[i]=wilcox.test(score.random[ires],score.random[iirs],alternative="greater")$p.value}
sum(st.random<0.009062)/1000
# empirical P < 0.023

##### control:ncSL partners
load("data/TCGA.drug.response/sl.pairs.patient.drug.response.RData")
targets=c("TUBB1","BCL2","MAPT","MAP2","MAP4") #Taxane targets
ix=which(sl.tot[,1] %in% match(targets,prob$genes))
sl.tot=sl.tot[ix,]

##### gene ID mapping to ovarian cancer dataset
sl=cbind(prob$genes[sl.tot[,1]],prob$genes[sl.tot[,2]])
sl=cbind(match(sl[,1],prob.ov$genes),match(sl[,2],prob.ov$genes))
sl.tot=sl.tot[!is.na(sl[,2]),]
sl=sl[!is.na(sl[,2]),]
sl.tot[,1:2]=sl	

ig=sl.tot[order(sl.tot[,3])[seq(n.partners)],2]
score.ncSL=colSums(prob.ov$mRNAq[ig,]==0,na.rm=T)
st.ncSL=wilcox.test(score.ncSL[ires],score.ncSL[iirs],alternative="greater")$p.value
# ncSL P < 0.7854895

##### control:DAISY SL partners
load("data/prob.TCGA.RData")
load("data/TCGA.drug.response/sl.pairs.patient.drug.response.RData")
targets=c("TUBB1","BCL2","MAPT","MAP2","MAP4")
ix=which(sl.tot[,1] %in% match(targets,prob$genes))
sl.tot=sl.tot[ix,]
sl=sl.tot[,1:2]

mRNA=prob$mRNA;mRNA.rank2=prob$mRNA.rank2;scna=prob$scna

dmol1=do.call(rbind,mclapply(seq(nrow(sl)),function(X){
	gene1=sl[X,1];gene2=sl[X,2]
	q1=quantile(mRNA[gene1,],0.1,na.rm=T)
	dn=which(scna[gene1,]< -0.2 & mRNA[gene1,]< q1)
	up=as.numeric(which(scna[gene1,]>0))
	st.mRNA1=st.scna1=st.mRNA2=st.scna2=NA
	if (sum(!is.na(dn))>100 & sum(!is.na(up))>100 & sum(!is.na(mRNA[gene2,]))>100 & sum(!is.na(scna[gene2,]))>100 ){
		st.mRNA1=wilcox.test(mRNA[gene2,dn],mRNA[gene2,up],alternative="greater")$p.value
		st.scna1=wilcox.test(scna[gene2,dn],scna[gene2,up],alternative="greater")$p.value
	}
	q1=quantile(mRNA[gene2,],0.1,na.rm=T)
	dn=which(scna[gene2,]< -0.2 & mRNA[gene1,]< q1)
	up=as.numeric(which(scna[gene2,]>0))
	if (sum(!is.na(dn))>100 & sum(!is.na(up))>100 & sum(!is.na(mRNA[gene1,]))>100 & sum(!is.na(scna[gene1,]))>100 ){
		st.mRNA2=wilcox.test(mRNA[gene1,dn],mRNA[gene1,up],alternative="greater")$p.value
		st.scna2=wilcox.test(scna[gene1,dn],scna[gene1,up],alternative="greater")$p.value
	}
	print(X)
	return(c(st.mRNA1,st.scna1,st.mRNA2,st.scna2))
},mc.cores=20))

cr1=do.call(rbind,mclapply(seq(nrow(sl)),function(X){
	gene1=sl[X,1];gene2=sl[X,2]
	cr=cor.test(prob$mRNA[gene1,],prob$mRNA[gene2,],method="spearman")$estimate
},mc.cores=20))

p.mRNA=apply(dmol1[,c(1,3)],1,max,na.rm=T)
p.scna=apply(dmol1[,c(2,4)],1,max,na.rm=T)

ig=which(p.adjust(sl.tot[,3],"BH")<0.05 & p.adjust(p.mRNA,"BH")<0.05 & p.adjust(p.scna,"BH")<0.05 & cr1>0.3) 
ix=match(unique(prob$genes[sl.tot[ig,2]]),prob.ov$genes)
ig=ix[!is.na(ix)]

score.daisy=colSums(prob.ov$mRNAq[ig,]==0,na.rm=T)
st.daisy=wilcox.test(score.daisy[ires],score.daisy[iirs],alternative="greater")$p.value
# DAISY P < 0.1851724

