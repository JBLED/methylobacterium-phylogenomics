###Methylobacterium core genome determination and production of alignments - Starting for 213 Methylobacteriaceae genomes used in Leducq et al. 2022b (doi.org/10.1093/gbe/evac123) and enriched with 21 newly released genomes

####STEP 1 : Grep all genomes information from annotation files

#Package seqinr to handle nucleotide sequences
library(seqinr)

#Work directory
pR="/Users/jean-baptisteleducq/Desktop/METHYLOBACTERIUM/Github-R-Scripts-Methylo-Genomes/2-Determination-Methylobacterium-core-genome"

setwd(pR)

#File: List of 213 genomes analyzed in Leducq et al. 2022b (doi.org/10.1093/gbe/evac123) 
#Column 1: strain name
#Column 2: short (simplified) name
#Column 3: species as determined by ANI (threshold 97%) in Leducq et al. 2022b 
#Column 4: ??? 
#Column 5: group as determined in Leducq et al. 2022a (doi.org/10.1128/mbio.03175) 
#Column 6: group as determined in Aleesa et al 2022 (10.3389/fmicb.2021.740610 )
#Column 7: group as determined in Leducq et al. 2022b
#Column 8: source of the strain (biome)
#Column 9: myRAST annotation file
#Column 10: strain order in the SVQquartet lineage tree (Leducq et al. 2022b; Figures 1c and S3b)
#Column 11: strain order in the Astral lineage tree (Leducq et al. 2022b; Figures 1b and S3a)
#Column 12: strain order in the lineage tree based on gene content (Leducq et al. 2022b; Figures 2d and S9a)
#Column 13: strain order in the tree based on core genome synteny (Leducq et al. 2022b; Figures 3b and S9b)
#Column 14: composed name for labels (columns 3 and 1)
LIN<-read.table("Lineages.txt",header=T)

#File: List of 124 Methylobacteriaceae species determined in Leducq et al. 2022b by ANI (with species tree information)
#Column 1: species order in the Astral species tree (Leducq et al. 2022b; Figure S5a)
#Column 2: species order in the SVDquarted species tree (Leducq et al. 2022b; Figure S5b)
#Column 3: group as determined in Leducq et al. 2022b 
#Column 4: species as determined in Leducq et al. 2022b 
#Column 5: strain short name(s) assigned to each species
SPE<-read.table("Species.txt",header=T)

#List of 21 newly released Methylobacterium genomes, annotated with myRast on the 06-12-2022
NEW<-read.table("New-Genomes.txt",header=T)

#Grap only relevant annotation files and calculate genome statistics
ANN<-as.vector(LIN[,9])
NEW <-subset(NEW,NEW[,7]=="no")
ANNnew<-NEW[,5]

#Annotation statistics for the 213 original genomes (!!!!!!! STORED IN 0-Original-codes-and-ressources !!!!!!)
SUMann<-t(sapply(ANN,function(x){
	setwd(paste(pR,"Annotated-Genomes",sep="/"))
	ANNt<-read.delim(x,header=T)
	#print(dim(ANNt))
	
	#scaffold size based on annotation
	Sx<-sort(sapply(unique(ANNt[,1]),function(x){
		return(max(subset(ANNt[,c(5:6)],ANNt[,1]==x)))
	}),decreasing=T)
	#N50
	N50<-subset(Sx,ifelse(cumsum(Sx)>=(sum(Sx)/2),1,0)==1)[1]
	
	GC<-GC(unlist(strsplit(ANNt[,12],split="")))
	
	#return number of scaffolds, number of annotations, number of unique annotations, total genome size, N50 and GC content (coding)
	setwd(paste(pR,"Annotated-Genomes-to-use",sep="/")) ####YOU HAVE TO CREATE THIS FOLDER AND UPLOAD ALL ANNOTATION FILES IN THERE
	write.table(ANNt,x,col.names=T)
	
	return(c(
	length(unique(ANNt[,1])),
	length(ANNt[,8]),
	length(unique(ANNt[,8])),
	sum(Sx),N50,GC
	))	
	
}))

#Annotation statistics for the 21 newly released and annotated genomes
SUMannNew<-t(sapply(ANNnew,function(x){
	setwd(paste(pR,"New-Genomes-06-12-2023",sep="/"))
	ANNt<-read.delim(x,header=T)
	#print(dim(ANNt))
	
	#scaffold size based on annotation
	Sx<-sort(sapply(unique(ANNt[,1]),function(x){
		return(max(subset(ANNt[,c(5:6)],ANNt[,1]==x)))
	}),decreasing=T)
	#N50
	N50<-subset(Sx,ifelse(cumsum(Sx)>=(sum(Sx)/2),1,0)==1)[1]
	
	GC<-GC(unlist(strsplit(ANNt[,12],split="")))
	
	#return number of scaffolds, number of annotations, number of unique annotations, total genome size, N50 and GC content (coding)
	setwd(paste(pR,"Annotated-Genomes-to-use",sep="/"))
	write.table(ANNt,x,col.names=T)
	
	return(c(
	length(unique(ANNt[,1])),
	length(ANNt[,8]),
	length(unique(ANNt[,8])),
	sum(Sx),N50,GC
	))	
	
}))

colnames(SUMann)<-c("Contigs","Genes","Annotations",
	"Size","N50","GC")
colnames(SUMannNew)<-c("Contigs","Genes","Annotations",
	"Size","N50","GC")

MaxSt<-max(as.numeric(setdiff(unlist(strsplit(LIN[,2],split="st")),"")))
NewSt<-paste("st",c(1:length(NEW[,1]))+ MaxSt,sep="")
	
SUMannNew<-cbind(NEW[,c(3,4,5)], NewSt,"New", SUMannNew)	
SUMann<-cbind(LIN[,c(1,8,9,2,7)], SUMann)

colnames(SUMannNew)<-colnames(SUMann)

SUMann<-rbind(SUMann, SUMannNew)

setwd(pR)

write.table(SUMann,"Annotation-Summary.txt",row.names=F)

####STEP 2 - Calculate raw gene abundance per group

setwd(pR)

SUMann <-read.table("Annotation-Summary.txt",header=T)
#SUMann <-subset(SUMann, SUMann[,5]!="M")
#SUMann <-subset(SUMann, SUMann[,5]!="E")

pdf("Genome-summary.pdf",width=5,height=5)

par(mar=c(4,4,1,1),bty="n")

plot(SUMann[,9], SUMann[,11],
	col=ifelse(SUMann[,5]=="A","red",
		ifelse(SUMann[,5]=="B","purple",
		ifelse(SUMann[,5]=="C","green",
		ifelse(SUMann[,5]=="D","blue",
		ifelse(SUMann[,5]=="New","yellow2","grey"))))),
	pch=19,cex=0.5,las=1,cex.axis=0.8,xlab="Genome Size",
	ylab="GC content")

legend("topright",cex=0.6, horiz=F,
	legend= c("A","B","C","D",
	"New genomes","Outgroups"),
	text.col= "black",box.col=NA,
	border=NA,col= c("red","purple","green2","blue",
	"yellow2","grey"),
	title="Groups",title.col="black",ncol=1,pch=19)

dev.off()

##List of unique gene names

ANNx<-as.vector(SUMann[,3])

setwd(paste(pR,"Annotated-Genomes-to-use",sep="/"))

Fu<-sort(unique(unlist(sapply(ANNx,function(x){
	print(x)
	ANNt<-read.table(x,header=T)
	Fx<-as.vector(ANNt[,8])
	return(unique(Fx))
}))))

#Remove hypothetical proteins, repeat regions and Mobile element proteins
Fu<-setdiff(Fu,c("hypothetical protein","repeat region","Mobile element protein"))

#Retrieve the raw abundance of each gene per genome
Fall<-sapply(ANNx,function(x){
	print(x)
	ANNt<-read.table(x,header=T)
	Fx<-as.vector(ANNt[,8])	
	return(sapply(Fu,function(x){
		return(length(subset(Fx,Fx==x)))
	}))	
})

setwd(pR)

write.table(Fall,"Gene-copy-number-raw.txt")

#STEP 3: Proportion of genes with 1,2,3,4 or more copies in fonction of the number of scaffolds 
setwd(paste(pR,"Annotated-Genomes-to-use",sep="/"))	

GeneCopy<-t(sapply(c(1:length(ANNx)),function(x){
	print(x)
	Gx=x
	ANNt<-read.table(ANNx[Gx],header=T)
	#number of scaffolds
	NS=length(unique(ANNt[,1]))
	#Number of genes features
	Gall<-length(subset(Fu,Fall[,Gx]!=0))
	#Number of genes that are in 1,2,3,4 or more copies
	G1<-length(subset(Fu,Fall[,Gx]==1))/Gall
	G2<-length(subset(Fu,Fall[,Gx]==2))/Gall
	G3<-length(subset(Fu,Fall[,Gx]==3))/Gall
	G4<-length(subset(Fu,Fall[,Gx]==4))/Gall
	G5<-length(subset(Fu,Fall[,Gx]==5))/Gall
	Gn<-length(subset(Fu,Fall[,Gx]>5))/Gall
	G0<-length(subset(Fu,Fall[,Gx]==0))/Gall	
	return(c(NS,Gall,c(G1,G2,G3,G4,G5,Gn,G0)))
}))

rownames(GeneCopy)<-ANNx

GeneCopy<-GeneCopy[order(GeneCopy[,1]),]


colnames(GeneCopy)<-c("Scaffolds","Genes","n=1_copy",paste("n=",c(2:5),"_copies",sep=""),"n>5_copies","n=0_copy")

setwd(pR)

write.table(GeneCopy,"Gene-copy.txt")

#Figure showing per genome the proportion of genes in 1 to 5 copy in function of genome assembly quality as measured by the number of scaffolds (Fig S1; Leducq et al. 2022b)
pdf("GeneCopy-function-scaffold.pdf",height=5,width=4)

par(mar=c(4,4,1,1),bty="n")

plot(GeneCopy[,1], GeneCopy[,3],xlab="number of scaffolds",ylab="proportion of genes",las=1,cex.axis=0.8,log="x",ylim=c(0,1),pch=19,col=rgb(0,0,1,0.4),cex=0.5)

points(GeneCopy[,1], GeneCopy[,4],pch=19,col=rgb(0,1,0,0.4),cex=0.5)

points(GeneCopy[,1], GeneCopy[,5],pch=19,col=rgb(1,0.7,0,0.4),cex=0.5)

points(GeneCopy[,1], GeneCopy[,6],pch=19,col=rgb(1,0,0,0.4),cex=0.5)

points(GeneCopy[,1], GeneCopy[,7],pch=19,col=rgb(0.7,0,1,0.4),cex=0.5)

text(seq(5,5,length.out=5),seq(0.6,0.4,length.out=5),c("Single copy","2 copies","3 copies","4 copies","5 copies"),col=c("blue","green","orange","red","purple"),cex=0.8)

dev.off()

 
#STEP 4: estimate the size of Methylobacterium core genome by random ressampling of genomes within groups A, B, C and D and counting New genomes as an independant category. Number of core genes in function of the number of genomes (per group, in total in Methylobacterium, without accounting for New genomes)

setwd(pR)

Fall<-read.table("Gene-copy-number-raw.txt",header=T)
row.names(Fall)<-NULL
SUMann <-read.table("Annotation-Summary.txt",header=T)
#SUMann <-subset(SUMann, SUMann[,5]!="M")
#SUMann <-subset(SUMann, SUMann[,5]!="E")

#Number of genomes in each group
Ng<-summary(as.factor(SUMann[,5]))

#Remove new genomes (not classified) and Enterovirga (only 2 genomes)
Ng<-subset(Ng,names(Ng)!="E")
Ng<-subset(Ng,names(Ng)!="New")

#Maximum number of genomes per group to consider
Gmax<-min(Ng)

#Groups to consider
MG=names(Ng)

#Number of genomes to ressample per group
ITg<-c(2:Gmax)

#Number of iterations
K=100

CORE<-sapply(ITg,function(x){ #iteration level
	print(x)
	X=x
	
	NX<-sapply(c(1:K),function(x){ #replicate level
		#print(x)
		#Sample X genomes per group
		SG<-sapply(MG,function(x){
			return(sample(grep(x, SUMann[,5]),X))
		})
		
		#Number of core genes per group
		NCx<-sapply(c(1:length(MG)),function(x){ #Group level
			Fx<-Fall[,SG[,x]]
			return(dim(subset(Fx,
				paste(apply(Fx,1,max), 
					apply(Fx,1,min),
				sep="")== "11"))[1])
		})
		#Number of core genes in Methylobacterium 
		Fx<-Fall[,as.vector(SG[,c(1:4)])]
		NMx<-dim(subset(Fx,
			paste(apply(Fx,1,max), 
			apply(Fx,1,min),
			sep="")== "11"))[1]
			
		#Number of core genes in ABD
		Fx<-Fall[,as.vector(SG[,c(1,2,4)])]
		ABDx<-dim(subset(Fx,
			paste(apply(Fx,1,max), 
			apply(Fx,1,min),
			sep="")== "11"))[1]
			
		#Number of core genes in AB
		Fx<-Fall[,as.vector(SG[,c(1,2)])]
		ABx<-dim(subset(Fx,
			paste(apply(Fx,1,max), 
			apply(Fx,1,min),
			sep="")== "11"))[1]	

		#Number of core genes in AD
		Fx<-Fall[,as.vector(SG[,c(1,4)])]
		ADx<-dim(subset(Fx,
			paste(apply(Fx,1,max), 
			apply(Fx,1,min),
			sep="")== "11"))[1]	

		#Number of core genes in BD
		Fx<-Fall[,as.vector(SG[,c(2,4)])]
		BDx<-dim(subset(Fx,
			paste(apply(Fx,1,max), 
			apply(Fx,1,min),
			sep="")== "11"))[1]	
			
		#Number of core genes in Microvirga
		Fx<-Fall[,as.vector(SG[,5])]
		Nmx<-dim(subset(Fx,
			paste(apply(Fx,1,max), 
			apply(Fx,1,min),
			sep="")== "11"))[1]	
			
		#Number of core genes in Methylobacteriaceae
		Fx<-Fall[,as.vector(SG)]
		NAx<-dim(subset(Fx,
			paste(apply(Fx,1,max), 
			apply(Fx,1,min),
			sep="")== "11"))[1]															
						
		return(c(NCx,NMx,ABDx,ABx,ADx,BDx,Nmx,NAx))
		
	})	
		
	return(round(c(apply(NX,1,mean),apply(NX,1,sd)),3))
})

setwd(pR)

pdf(paste("Annot-Raw-Rarefaction_K=",K,".pdf",sep=""),height=7,width=5)

plot(100000000,100000000,
	xlim=c(0.4,max(ITg))*5,
	ylim=c(200,max(CORE)),log="x",
	xlab="Genomes (log scale)",ylab="Core genome size (genes)",cex.axis=0.8,las=1,
	main="")

#Graphique parameters
COLg<-cbind(setdiff(c(1:12),11),
	c("red","purple","green2","blue",
	"grey","brown","yellow2","black",
	"black","black","black"),
	c("A","B","C","D",
	"Micr","Meth","ABD","AB",
	"AD","BD","All"),
	c(1,1,1,1,1,1,1,1,2,3,1), #line type (lty)
	c(1,1,1,1,1,1,1,1,1,1,1.5), #line width
	c(1,1,1,1,1,4,3,2,2,2,5)) #coefficient multiplication
	
sapply(COLg[,1],function(x){
	COLx<-as.vector(subset(COLg, COLg[,1]==x))
	points(ITg*as.numeric(COLx[6]),
		CORE[as.numeric(x),],col= COLx[2],
		type="l",lty=as.numeric(COLx[4]),
		lwd=as.numeric(COLx[5]))
})

sapply(COLg[,1],function(x){
	COLx<-as.vector(subset(COLg, COLg[,1]==x))
	points(max(ITg)*as.numeric(COLx[6]),
		CORE[as.numeric(x),(length(CORE[1,]))],
		col= COLx[2],
		cex=2.5,pch=21,bg="white",
		lwd=as.numeric(COLx[5]))
	text(max(ITg)*as.numeric(COLx[6]),
		CORE[as.numeric(x),(length(CORE[1,]))],
		cex=0.5,labels=COLx[3])				
})

dev.off()

colnames(CORE)<-paste("n=", ITg,sep="")
rownames(CORE)<-c(paste("mean",c("A","B","C","D",
	"Micr","Meth","ABD","AB",
	"AD","BD","Micr","All"),sep=""),
	paste("sd",c("A","B","C","D",
	"Micr","Meth","ABD","AB",
	"AD","BD","Micr","All"),sep=""))
	
write.table(CORE,paste("Annot-Raw-Rarefaction_K=",K,".txt",sep=""))	
	
#legend("topright",cex=0.6, horiz=F,
#	legend= c("A","B","C","D",
#	"Methylobacterium","ABD","AB",
#	"AD","BD"),
#	text.col= "black",box.col=NA,
#	border=NA,col= c("red","purple","green2","blue",
#	"grey","yellow3","black","black","black"),
#	title="",title.col="black",ncol=1,
#	lty=c(1,1,1,1,1,1,1,2,3))

# STEP 5: look at the number of core genes per group in function of different threshold (>=0.8)

TH<-seq(0.8,1,length.out=21)

#In each group separately
Pcore<-sapply(c(1:length(MG)),function(x){ #Group level
	print(MG[x])
	SG<-grep(MG[x], SUMann[,5]) #column numbers for this group
	Fx<-Fall[,SG] #gene occurence in this group
	Nc<-sapply(c(1:length(Fx[,1])),function(x){
		return(length(subset(as.numeric(Fx[x,]),
		as.numeric(Fx[x,])==1)))
	})/length(SG) #proportion of genes with 1 copy in this group
	
	Nc<-sapply(TH,function(x){
		return(length(subset(Nc,Nc>=x)))
	})
	return(Nc)
})

####In Methylobacterium
SG<-c(grep("A", SUMann[,5]),
	grep("B", SUMann[,5]),
	grep("C", SUMann[,5]),
	grep("D", SUMann[,5]))
Fx<-Fall[,SG] #gene occurence in this group
Nc<-sapply(c(1:length(Fx[,1])),function(x){
	return(length(subset(as.numeric(Fx[x,]),
	as.numeric(Fx[x,])==1)))
})/length(SG) #proportion of genes with 1 copy in this group
	
NcMeth<-sapply(TH,function(x){
	return(length(subset(Nc,Nc>=x)))
})	

####In Methylobacteriaceae (All genomes)

Nc<-sapply(c(1:length(Fall[,1])),function(x){
	return(length(subset(as.numeric(Fall[x,]),
	as.numeric(Fall[x,])==1)))
})/length(Fall[1,]) #proportion of genes with 1 copy in this group
	
NcAll<-sapply(TH,function(x){
	return(length(subset(Nc,Nc>=x)))
})
	
Pcore<-cbind(Pcore, NcMeth, NcAll)

colnames(Pcore)<-c(MG[1:4],"Micr","Meth","All")
rownames(Pcore)<-paste("P>=",TH,sep="")

###There are 995 candidate core genes in Methylobacterium, assuming genes having 1:1 copy in 95% of genomes 
	
write.table(Pcore ,"Annot-Numb-Core-Gene-per-group.txt")
	
#STEP 6: refine the copy number estimation for these genes in all genomes using the sequence size normalized by size in complete genomes - look at distributions and define threshold to keep only 1:1 genes within groups AND within Methylobacterium

library(seqinr)

setwd(pR)

SUMann <-read.table("Annotation-Summary.txt",header=T)
Fall<-read.table("Gene-copy-number-raw.txt",header=T)

#Simflify gene names
Fu<-rownames(Fall)
ng=nchar(length(Fu))

FuMOD<-sapply(c(1:length(Fu)),function(x){
	paste(c("Gene",seq(0,0,length.out=(ng-nchar(x))),x),collapse="")
})

write.table(cbind(FuMOD,Fu),"Gene-modified-names.txt",row.names=F,col.names=c("FuMOD","Fu"))

#Complete genome indexes ( with N50 > 3.10^6)
COMPi<-subset(c(1:length(SUMann[,1])),SUMann[,10]>=(3*10^6))

#candidate core genes in Methylobacterium, assuming genes having 1:1 copy in NCx% of Methylobacterium genomes

NCx=0.90

####In Methylobacterium
SG<-c(grep("A", SUMann[,5]),
	grep("B", SUMann[,5]),
	grep("C", SUMann[,5]),
	grep("D", SUMann[,5]))
Fx<-Fall[,SG] #gene occurence in this group
Nc<-sapply(c(1:length(Fx[,1])),function(x){
	return(length(subset(as.numeric(Fx[x,]),
	as.numeric(Fx[x,])==1)))
})/length(SG) #proportion of genes with 1 copy in this group

Fcore<-subset(c(1:length(Fu)), Nc>= NCx)

# Average size of each gene among complete genomes only
setwd(paste(pR,"Annotated-Genomes-to-use",sep="/"))

FcoreSizeComp<-sapply(COMPi,function(x){
	print(x)
	Gx=x
	ANNt<-read.table(SUMann[Gx,3],header=T)
	return(as.numeric(unlist(
		sapply(Fu[Fcore],function(x){
			return(c(nchar(subset(as.vector(ANNt[,12]),
			as.vector(ANNt[,8])==x)),NA)[1])
	}))))
})

AVEcomp<-apply(FcoreSizeComp,1,mean,na.rm=T)
SDcomp<-apply(FcoreSizeComp,1,sd,na.rm=T)

#IN ALL GENOMES (complete and draft): Average size of sequences from a same gene, normalized by the average size of single copies found in complete genomes)
FcoreSize<-sapply(c(1:length(SUMann[,1])),function(x){
	print(x)
	Gx=x
	ANNt<-read.table(SUMann[Gx,3],header=T)
	z<-as.vector(sapply(Fu[Fcore],function(x){		
		mean(nchar(subset(as.vector(ANNt[,12]),
			as.vector(ANNt[,8])==x)),na.rm=T)	
	}))/AVEcomp
	return(z)
})

#Raw gene copy number
FcoreCount<-sapply(c(1:length(SUMann[,1])),function(x){
	print(x)
	Gx=x
	ANNt<-read.table(SUMann[Gx,3],header=T)
	nz<-as.vector(sapply(Fu[Fcore],function(x){	#	
		length(subset(as.vector(ANNt[,12]),
			as.vector(ANNt[,8])==x))	
	}))
	return(nz)
})

rownames(FcoreCount)<-FuMOD[Fcore]
rownames(FcoreSize)<-FuMOD[Fcore]

setwd(pR)
write.table(FcoreCount,paste("RawFcoreCount_TH>=",NCx,".txt",sep=""),row.names=T,col.names=F)
write.table(FcoreSize,paste("RawFcoreSize_TH>=",NCx,".txt",sep=""),row.names=T,col.names=F)


vx<-seq(-0.25,0.25,length.out=101)

###Figure showing gene copy size (normalized by average size in complete genomes) in function of the observed number of copy in each genome before filtering true duplicates (Figure S2 in Leducq et al 2022b)

pdf(paste("NormalizedSize-vs-CopyNumber_TH>=",NCx,".pdf",sep=""),width=4,height=5)

plot(-10,-10,xlim=c(0,max(as.numeric(FcoreCount))),
	ylim=c(0,max(as.numeric(FcoreSize),na.rm=T)),
	las=1,cex.axis=0.8,
	xlab="observed copy number",
	ylab="normalized gene size",
	#main=paste("Before filtering (",length(Fcore)," genes)",sep="")
	)


ne<-seq(1,max(FcoreCount),length.out=101)

points(ne,1/ne,type="l",
	col="black",lwd=2,lty=1)
	
points(ne,2/ne,type="l",
	col="black",lwd=2,lty=2)
	
points(ne,3/ne,type="l",
	col="black",lwd=1.5,lty=3)	

points(ne,4/ne,type="l",
	col="black",lwd=1,lty=3)	

sapply(c(1:length(Fcore)),function(x){
	
	points(FcoreCount[x,]+sample(vx,1),FcoreSize[x,],pch=19,cex=0.1,col=rgb(0.7,0.7,0.7,0.05))
})


legend(x = "topright",cex=0.8, horiz=F,legend= paste("n=",c(1:4),sep=""),text.col= "black",box.col=NA,col= c("black"),border=NA,text.font=1,lty=c(1,2,3,3),title="Real copy number",pt.lwd=c(2,2,1.5,1))

dev.off()

#The real gene copy number can be approximated by FcoreCount*FcoreSize (assuming very few or no overlap between contigs in draft genomes resulting in false duplicates)


#For each gene, count the number of genomes with at least 2 copies with an average relative size of at least 0.75  (potentially true duplicated gene if counted in at least one genome) === Remove outgroup genomes

METHi<-c(grep("A", SUMann[,5]),
	grep("B", SUMann[,5]),
	grep("C", SUMann[,5]),
	grep("D", SUMann[,5]),
	grep("New", SUMann[,5]))

DUP<-apply(ifelse(FcoreCount[, METHi]<2,0,ifelse(is.na(FcoreSize[, METHi])==T,0,ifelse(FcoreSize[, METHi]>0.75,1,0))),1,sum) 

CORE<-apply(ifelse(round(FcoreCount[, METHi]*FcoreSize[, METHi],0)==1,1,0),1,sum,na.rm=T) 


#Remove potentially duplicated genes in at least one genome (check that "DNA-directed RNA polymerase beta subunit (EC 2.7.7.6)"==rpoB is conserved after this step)
FcoreNEW<-subset(Fcore, DUP<1)
FcoreCountNEW<-subset(FcoreCount, DUP<1)
FcoreSizeNEW<-subset(FcoreSize, DUP<1)

#Figure showing gene copy size (normalized by average size in complete genomes) in function of the observed number of copy in each genome after filtering true duplicates 

setwd(pR)

pdf(paste("Filtered-NormalizedSize-vs-CopyNumber_TH>=",NCx,".pdf",sep=""),width=4,height=5)

plot(-10,-10,xlim=c(0,max(FcoreCountNEW)),ylim=c(0,max(FcoreSizeNEW,na.rm=T)),las=1,cex.axis=0.8,xlab="copy number",ylab="normalized size",main=paste("After filtering (",length(FcoreNEW)," genes)",sep=""))

points(ne,1/ne,type="l",
	col="red",lwd=2,lty=1)

sapply(c(1:length(FcoreNEW)),function(x){
	
	points(FcoreCountNEW[x, METHi]+sample(vx,1), FcoreSizeNEW[x, METHi],pch=19,cex=0.1,col=rgb(0,0,0,0.1))
})


dev.off()

#For each remaining gene and each genome, only keep single copies that have a relative size of more than 0.7 and less than 1.3 

SINGcountNEW<-ifelse(is.na(FcoreSizeNEW)==T,0,ifelse(FcoreSizeNEW<0.7,0,ifelse(FcoreSizeNEW>1.3,0,1)))

SINGcountNEW<-ifelse(FcoreCountNEW==1, SINGcountNEW,0)

#Graphique parameters
COLg<-cbind(setdiff(c(1:12),11),
	c("red","purple","green2","blue",
	"grey","brown","yellow2","black",
	"black","black","black"),
	c("A","B","C","D",
	"Micr","Meth","ABD","AB",
	"AD","BD","All"),
	c(1,1,1,1,1,1,1,1,2,3,1), #line type (lty)
	c(1,1,1,1,1,1,1,1,1,1,1.5), #line width
	c(1,1,1,1,1,4,3,2,2,2,5)) #coefficient multiplication


COLcla<-sapply(SUMann[,5],function(x){
	return(c(subset(COLg[,2], COLg[,3]==x),ifelse(x=="New","yellow2","grey"))[1])
})

pdf("Gene-occurence-after-filtering.pdf",width=5,height=5)

Ng=200

plot(as.numeric(SUMann[,10])/10^6,apply(SINGcountNEW,2,sum),log="x",las=1,cex.axis=0.8,xlab="N50 (Mb)",ylab="Core genes",bg=COLcla,pch=21,main=paste("Before removing scattered genes (",length(FcoreNEW),")",sep=""))

hist(apply(SINGcountNEW[,METHi],1,sum),breaks=50,las=1,cex.axis=0.8,xlab="Occurence in genomes (Methylobacterium only)",ylab="Gene count",main="Distribution of single copies after filtering",cex.main=0.8)

segments(Ng,0, Ng,1000,lwd=2,lty=2,col="red")

#heatmap.2(subset(SINGcountNEW,apply(SINGcountNEW,1,sum)>=175),scale="none",labCol=NAMES[,2],labRow="",cexCol=0.8,margins=c(10,2),trace="none",key=F,col=c("white",rgb(0,1,0,1)),ColSideColors= COLcla)

#Remove genes for which we have a copy in less than Ng genomes


FcoreNEW<-subset(FcoreNEW, apply(SINGcountNEW[,METHi],1,sum)>= Ng)
FcoreCountNEW<-subset(FcoreCountNEW, apply(SINGcountNEW[,METHi],1,sum)>= Ng)
FcoreSizeNEW<-subset(FcoreSizeNEW, apply(SINGcountNEW[,METHi],1,sum)>= Ng)
SINGcountNEW<-subset(SINGcountNEW, apply(SINGcountNEW[,METHi],1,sum)>= Ng)

plot(as.numeric(SUMann[,10])/10^6,apply(SINGcountNEW,2,sum),log="x",las=1,cex.axis=0.8,xlab="N50 (Mb)",ylab="Core genes",bg=COLcla,pch=21,main=paste("After removing scattered genes (",length(FcoreNEW),")",sep=""))

dev.off()

#########For each genome, extract the nucleotide sequence of single copy genes present in at least Ng genomes and that have a relative sequence size in the range 0.7-1.3 from GFF annotation file


NcorePerGenome<-
sapply(c(1:length(SUMann[,1])),function(x){
	#print(x)
	Gx=x
	Sx=SUMann[x,4]
	setwd(paste(pR,"Annotated-Genomes-to-use",sep="/"))
	ANNt<-read.table(SUMann[Gx,3],header=T)
	#Only keep Fcore genes that are present in this genome
	Fx=Fu [subset(FcoreNEW ,SINGcountNEW[,Gx]==1)]
	#Retrieve modified gene names (without space)
	FxMOD<-as.vector(sapply(Fx,function(x){
		return(subset(FuMOD,Fu==x))
	}))
	SEQx<-sapply(Fx,function(x){
		GN=x
		return(list(unlist(strsplit(
			subset(as.vector(ANNt[,12]),
			as.vector(ANNt[,8])==GN),split=""))))
	})		
	setwd(paste(pR,"/OUT-FASTA",sep="")) #YOU HAVE TO CREATE THIS FOLDER
	write.fasta(SEQx,names= paste(Sx, FxMOD,sep="_"),paste("Core-Genes-", Sx,".fas",sep=""))
	print(paste(c("Strain",Sx,":",length(Fx),"core genes",";",length(SEQx),"sequences"),collapse=" "))
	return(length(Fx))
})

setwd(pR)
pdf("CoreGene-per-genome-N50.pdf",width=5,height=5)

plot(as.numeric(SUMann[,10]), NcorePerGenome,col= COLcla,log="x",las=1,cex.axis=0.8,ylab="Core gene sequences retrieved per genome",xlab="N50 (pb, log scale)",cex=1,pch=ifelse(SUMann[,5]=="M",4,ifelse(SUMann[,5]=="E",4,1)))

dev.off()

#Step 7: For each core gene, produce an alignement with genomes for which sequences are available

setwd(paste(pR,"/OUT-FASTA",sep="")) 

FcoreCountGen<-sapply(FcoreNEW,function(x){
	Fx= FuMOD[x]	
	
	#count the number of the gene occurence per genome
	Nx=as.vector(subset(SINGcountNEW, FcoreNEW== x))
	
	#remove genomes where the gene is absent
	ANx<-subset(c(1:length(SUMann[,1])),Nx==1)
	SEQ<-sapply(ANx,function(x){
		#print(x)
		Gx=x
		Sx=SUMann[x,4]
		SEQx<-read.fasta(paste("Core-Genes-", Sx,".fas",sep=""))
		NS=paste(Sx, Fx,sep="_")
		return(subset(SEQx,names(SEQx)==NS))
	})
	write.fasta(SEQ,paste(Fx,".fas",sep=""),names=as.vector(SUMann[ANx,4]))
	print(paste(c("Gene",Fx,": in",sum(Nx),"genomes"),collapse=" "))
	return(sum(Nx))

})

AVEcompNew<-sapply(FcoreNEW,function(x){
	return(subset(AVEcomp,Fcore==x))
})

FcoreCountGen <-as.data.frame(cbind(Fu[FcoreNEW],FuMOD[FcoreNEW], FcoreCountGen, AVEcompNew))

colnames(FcoreCountGen)<-c("Annotation","Fasta-name","Occurence","AverageSize(bp)")

setwd(pR)

write.table(FcoreCountGen,"Summary-per-core-gene.txt",row.names=F)

###Summarize the occurence of core gene per genome in the final dataset

setwd(pR)

pdf("summary-core-genes-final.pdf")

heatmap(SINGcountNEW,scale="none",ColSideColors=COLcla,
	col=c("black","green2"))

dev.off()

write.table(SINGcountNEW,
	paste("FinalFcoreCount_TH>=",NCx,".txt",sep=""),
	col.names=F,row.names=T)

###############
###STEP 8: align sequences for each core gene

library(seqinr)
library(msa)

pR="/Users/jean-baptisteleducq/Desktop/Removed-DB-08-08-2022/PostDoc-Moscow-ID/Project/Article-Methylo-genome-evolution/for-Github/2-Determination-Methylobacterium-core-genome"

setwd(pR)

NCx=0.9

#List of genomes
SUMann <-read.table("Annotation-Summary.txt",header=T)

#List of core genes
FcoreCountGen <-read.table("Summary-per-core-gene.txt",header=T)

#Occurence of core genes per genome
SINGcountNEW<-read.table(paste("FinalFcoreCount_TH>=",
	NCx,".txt",sep=""),row.names=1)

#Type of gene (coding(peg)/rna); use a genome with all core genes, and retrieve the information in annotation file

REFann=subset(SUMann[,3],
	apply(SINGcountNEW,2,sum)==
	length(FcoreCountGen[,1]))[1]

setwd(paste(pR,"Annotated-Genomes-to-use",sep="/"))
	
ANNOx<-read.table(REFann,header=T)	

TYPE<-sapply(c(1:length(FcoreCountGen[,1])),
	function(x){
	Ax<-FcoreCountGen[x,1]
	return(unique(subset(ANNOx[,3], ANNOx[,8]==Ax)))
})

#Check that each gene alignement size is a multiple of 3
setwd(paste(pR,"OUT-FASTA",sep="/"))

SUM3<-sapply(c(1:length(FcoreCountGen[,1])),
	function(x){
	Ix<-FcoreCountGen[x,2] #short gene name
	print(Ix)
	#Get the unaligned nucleotide sequence
	SEQ<-read.fasta(paste(Ix,
		".fas",sep=""))
	z<-unique(as.numeric(summary(SEQ)[,1]))
	return(unique((3*round(z/3))==z))
})

subset(TYPE,SUM3==F) #only rna shouldn't be multiple of 3

####Align non protein coding genes
Irna<-subset(c(1:length(FcoreCountGen[,1])), TYPE=="rna")

setwd(paste(pR,"OUT-FASTA",sep="/"))

sapply(Irna,function(x){	
	Ix<-FcoreCountGen[x,2] #short gene name
	print(Ix)
	SEQw<-msa(paste(Ix,
		".fas",sep=""),method="ClustalW",type="dna",order="input")
	#Retrieve aligned sequence in seqinr format
	SEQn<-msaConvert(SEQw,type="seqinr::alignment")$seq
	#Get sequence names from the unaligned nucleotide sequence
	NX=names(read.fasta(paste(Ix,
		".fas",sep="")))
	#Convert
	SEQn<-sapply(c(1:length(SEQn)),function(x){
		return(strsplit(SEQn[x],split=""))
	})
	write.fasta(SEQn, names=NX,
		file.out =paste(Ix,
		"_aligned.fas",sep=""))
})


####Align protein encoding genes
Ipeg<-subset(c(1:length(FcoreCountGen[,1])), TYPE=="peg")

setwd(paste(pR,"OUT-FASTA",sep="/"))

sapply(Ipeg,function(x){
	
	Ix<-FcoreCountGen[x,2] #short gene name
	print(Ix)

	#Get the unaligned nucleotide sequence
	SEQ<-read.fasta(paste(Ix,
		".fas",sep=""))

	#Sequence names
	NX=names(SEQ)	
	
	#Remove gaps
	SEQ<-sapply(c(1:length(SEQ)),function(x){
		return(list(subset(unlist(SEQ[[x]]),
			unlist(SEQ[[x]])!="-")))
	})
	
	#Sequence sizes
	Sx<-as.numeric(summary(SEQ)[,1])
	Sx<-c(Sx,min(Sx)-3)
	
	#Add a vector of "Ns" for potential missing data at the end and begining of each sequence (length: maximum-minimum sequence size)
	nx<-sapply(c(1:(max(Sx)-min(Sx))),function(x){
		return("N")
	})
	
	SEQ<-sapply(c(1:length(SEQ)),function(x){
		return(list(c(nx,as.vector(SEQ[[x]]),nx)))
	})
	
	#Convert into amino-acid sequences
	SEQaa<-getTrans(SEQ)

	#Write a temporary file with amino acid sequences unaligned
	write.fasta(SEQaa,names=NX,"tempAA.fas")


	#Align AA sequences together using ClustalW - IMPORTANT: keep same order as input, since sequence names are removed during analysis, and retrieve sequence names in the original unaligned fasta file


	SEQw<-msa("tempAA.fas",method="ClustalW",type="protein",order="input")

	#Retrieve aligned sequence in seqinr format
	SEQn<-msaConvert(SEQw,type="seqinr::alignment")$seq
	
	#Convert back in nucleotide sequence

	SEQn<-sapply(c(1:length(NX)),function(x){
		#Unaligned nucleotide sequence
		seqx<-as.vector(SEQ[[x]])
		#Unaligned protein sequence
		sequ<-unlist(strsplit(SEQaa[[x]],split=""))
		#Aligned protein sequence
		seqa<-unlist(strsplit(SEQn[[x]],split=""))
		#Unaligned nucleotide sequence coordinates
		coordx<-c(1:length(seqx))
		#Unaligned proteine sequence coordinates
		coordu<-matrix(coordx,ncol=3,byrow=T)
		#Removed stop codon (because it's automaticaly removed by msa alignement)
		coordu<-subset(coordu, sequ!="*")
		#Aligned protein sequence coordinates
		coorda<-matrix(c(1:(length(seqa)*3)),
			ncol=3,byrow=T)
		#Ungapped positions in Aligned protein sequence coordinates
		coordau<-sort(as.vector(subset(coorda, 
			seqa!="-")))
		#Gapped position in Aligned protein sequence coordinates
	coordag<-sort(as.vector(subset(coorda, 
		seqa=="-")))		
		#Nucleotide sequence (ungapped position)	
		seqx<-seqx[sort(as.vector(coordu))]
		#Nucleotide sequence (Gapped position)
		seqg<-ifelse(seq(0,0,length.out=
			length(coordag))==0,"-","")
		return(list(c(seqg, seqx)
			[order(c(coordag, coordau))]))
	})
	
	#Remove position that only contain missing data
	SEQn<-sapply(c(1:length(NX)),function(x){
		return(as.vector(unlist(	SEQn[[x]])))
	})

	SEQn<-subset(SEQn,sapply(c(1:length(SEQn[,1])),function(x){
		seqx<-SEQn[x,]
		return(length(subset(seqx, seqx=="-"))+
		length(subset(seqx, seqx=="N")))
	})!=length(NX))
	
	SEQn<-sapply(c(1:length(NX)),function(x){
		return(list(SEQn[,x]))
	})
		
	write.fasta(SEQn, names=NX,
		file.out =paste(Ix,
		"_aligned.fas",sep=""))
})


#Review each protein-coding alignment and replace "NNN(---)NNN" motifs by "NNN(NNN)NNN"

setwd(paste(pR,"OUT-FASTA",sep="/"))

sapply(c(Ipeg,Irna),function(x){
	
	Ix<-FcoreCountGen[x,2] #short gene name
	print(Ix)

	#Get the aligned nucleotide sequence
	SEQ<-read.fasta(paste(Ix,
		"_aligned.fas",sep=""))
	
	SEQn<-sapply(c(1:length(SEQ)),function(x){
		#print(x)
		seqx<-c(as.vector(SEQ[[x]]),"x","n","-","n","x","n","-","n")
		lx=length(seqx)
		z<-c(0,sapply(c(2: lx),function(x){
			ifelse(seqx[x-1]==seqx[x],1,0)
		}))
		seqxx<-subset(seqx,z==0)
		z<-subset(c(1:length(seqx)),z==0)
		
		y<-t(sapply(c(3:length(seqxx)),function(x){
			return(c(
			paste(z[(x-2):x],collapse="_"),
			paste(seqxx[(x-2):x],collapse="_")
			))
		}))
		
		y<-subset(y,y[,2]=="n_-_n")
		
		#positions to replace by 'Ns'
		pn<-unique(unlist(sapply(c(1:length(y[,1])),function(x){
			yx<-as.numeric(unlist(strsplit(y[x,1],split="_")))
			return(c(min(yx):max(yx)))
		})))
		sn<-sapply(pn,function(x){return("n")})
		#positions to keep
		pk<-setdiff(c(1:lx),pn)
		sk<-seqx[pk]
		
		seqn<-c(sn, sk)[order(c(pn,pk))]
		return(list(toupper(seqn[1:(lx-8)])))
	})
	
	write.fasta(SEQn,paste(Ix,"_alignedN.fas",sep=""),names=names(SEQ))
	
})	


#Check that gap replacement by Ns did not affected alignments

setwd(paste(pR,"OUT-FASTA",sep="/"))

#only "n" should be returned

unique(unlist(
	sapply(c(1:length(FcoreCountGen[,1])),function(x){
	
	Ix<-FcoreCountGen[x,2] #short gene name
	print(Ix)
	#Get the aligned nucleotide sequence (before and after N correcting)
	SEQ<-read.fasta(paste(Ix,
		"_aligned.fas",sep=""))	
	SEQn<-read.fasta(paste(Ix,
		"_alignedN.fas",sep=""))
	sumx<-unlist(sapply(c(1:length(SEQ)),function(x){
		seq<-as.vector(SEQ[[x]])
		seqn<-as.vector(SEQn[[x]])
		return(subset(seqn,seq!=seqn))
	}))
	print(c(x,unique(sumx))	)
	return(sumx)
})))


#Proportion of missing data (only for protein coding)
setwd(paste(pR,"OUT-FASTA",sep="/"))

PPX<-unlist(sapply(Ipeg,function(x){
	
	Ix<-FcoreCountGen[x,2] #short gene name
	print(Ix)
	#Get the aligned nucleotide sequence
	SEQn<-read.fasta(paste(Ix,
		"_alignedN.fas",sep=""))
	NX<-names(SEQn)
	
	#Convert into amino-acid sequences
	SEQaa<-getTrans(SEQn)
	
	#Convert in matrix
	SEQn<-sapply(c(1:length(SEQn)),function(x){
		return(as.vector(SEQaa[[x]]))
	})
	
	#Count the proportion of X per position
	pn<-sapply(c(1:length(SEQn[,1])),function(x){
		length(subset(SEQn[x,],SEQn[x,]=="X"))
	})/length(NX)
	return(pn)
}))

#Remove aa positions with more than ppn% of missing data (only for protein-coding)
ppn=0.9

setwd(paste(pR,"OUT-FASTA",sep="/"))

#Ony "n" should be returned


sumNEW<-t(sapply(Ipeg,function(x){
	
	Ix<-FcoreCountGen[x,2] #short gene name
	print(Ix)
	#Get the aligned nucleotide sequence
	SEQn<-read.fasta(paste(Ix,
		"_alignedN.fas",sep=""))
	NX<-names(SEQn)
	
	#Convert into amino-acid sequences
	SEQaa<-getTrans(SEQn)
	
	
	#Convert both in matrix
	SEQn<-sapply(c(1:length(SEQn)),function(x){
		return(as.vector(SEQn[[x]]))
	})
	SEQaa <-sapply(c(1:length(SEQaa)),
		function(x){
		return(as.vector(SEQaa[[x]]))
	})
	#Count the proportion of X per position
	pn<-sapply(c(1:length(SEQaa[,1])),
		function(x){
		length(subset(SEQaa[x,], SEQaa[x,]=="X"))
	})/length(NX)
	
	pn <-as.vector(matrix(t(cbind(pn,pn,pn)),
		ncol=1,byrow=T))
	
	SEQn<-subset(SEQn, pn<(1-ppn))
	
	SEQn<-sapply(c(1:length(NX)),function(x){
		return(list(SEQn[,x]))
	})
	
	write.fasta(SEQn,paste(Ix,
		"_alignedNF.fas",sep=""),names=NX)
	return(c(x,Ix,length(SEQn),length(pn),
		length(SEQn[[1]])))	
		
}))	

#Same with rna

sumNEWrna<-t(sapply(Irna,function(x){
	
	Ix<-FcoreCountGen[x,2] #short gene name
	print(Ix)
	#Get the aligned nucleotide sequence
	SEQn<-read.fasta(paste(Ix,
		"_alignedN.fas",sep=""))
	NX<-names(SEQn)
	
	write.fasta(SEQn,paste(Ix,
		"_alignedNF.fas",sep=""),names=NX)
	return(c(x,Ix,length(SEQn),length(SEQn[[1]]),
		length(SEQn[[1]])))	
		
}))	


sumNEWrna<-cbind(sumNEWrna,"rna")
sumNEW<-cbind(sumNEW,"peg")
sumNEW<-rbind(sumNEWrna, sumNEW)

setwd(pR)
colnames(sumNEW)<-c("Name1","Name2","Sequences","Size1(bp)","Size2(bp)","Type")
write.table(sumNEW,"Summary-per-core-gene-after-alignment.txt",row.names=F)	

