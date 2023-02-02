####Grep all genomes information from .fna files

library(seqinr)

pR="/Users/jean-baptisteleducq/Dropbox/PostDoc-Moscow-ID/Project/New-Genome-Analysis-Methylobacterium-ss"

setwd(pR)

#list of genomes downloaded on NCBI on sept 21, 2021 
GEN<-as.vector(read.table("Genome-file.txt",header=T)[,1])

pG<-"/Users/jean-baptisteleducq/Desktop/New-Methylo=Genomes-21-09-2021/Unzip-for-new-annotation"

setwd(pG)

ATo<-c("genome","isolate","MAG","sequence","shotgun","whole","complete","contig","scaffold","genomic","assembly","strain","DNA","WGS","project","data","chromosome","uncultured","new")

GENn<-t(sapply(GEN,function(x){
	print(x)
	Gx<-unlist(strsplit(x,split=".gz"))
	SEQx<-read.fasta(Gx)
	ATx<-sapply(c(1:length(SEQx)),function(x){
		return(attributes(SEQx[[x]])$Annot)
	})
	ATu<-unlist(strsplit(ATx,split=" "))
	ATu<-unlist(strsplit(ATu,split=","))
	ATu<-unlist(strsplit(ATu,split=":"))
	ATu<-summary(as.factor(ATu),maxsum=1000000)
	ATu<-subset(names(ATu),ATu>=0.9*length(SEQx))
	ATu<-setdiff(ATu,ATu[grep(">",ATu)])
	ATu<-setdiff(ATu,ATu[grep("scaffold",ATu)])	
	ATu<-setdiff(ATu,ATo)
	AT1<-unlist(strsplit(ATx[1],split=" "))
	AT1<-unlist(strsplit(AT1,split=","))
	AT1<-unlist(strsplit(AT1,split=":"))
	
	AT1<-AT1[sort(sapply(ATu,function(x){
		return(min(grep(x,AT1)))
	}))]
	AT1<-c(AT1[1],AT1[2],
		paste(AT1[3:length(AT1)],collapse=" "))
		
	nSCAx<-length(SEQx)
	GCx<-round(GC(unlist(SEQx)),4)
	N50x<-sort(as.numeric(summary(SEQx)[,1])	,
		decreasing=T)
	Sx=sum(N50x)
	N50x<-subset(N50x, cumsum(N50x)>=(0.5* Sx))[1]
		
	return(c(x,AT1, nSCAx, GCx, Sx, N50x))
}))

rownames(GENn)<-NULL

colnames(GENn)<-c("Assembly-file","Genus","Species","Strain","Scaffolds","GCcontent","TotalSize","N50")

setwd(pR)

write.table(GENn,"New-genome-collection.txt",row.names=F)

#install.packages("http://www.bioconductor.org/packages/release/bioc/bin/macosx/contrib/4.0/msa_1.22.0.tgz")


###GENOME ANNOTATION OF METHYLOBACTERIACEAE GENOMES WITH myRAST - Available genomes in october 2020

library(seqinr)
library(gplots)
library(ade4)

#work directory for inputs
pS="/Users/jean-baptisteleducq/Dropbox/PostDoc-Moscow-ID/Project/myRAST-gene_annotation"
#work directory for outputs
pR="/Users/jean-baptisteleducq/Dropbox/PostDoc-Moscow-ID/Project/New-Genome-Analysis-Methylobacterium-ss"

#work directory for annotation files (GFF3 myRAST outputs)
pO="/Users/jean-baptisteleducq/Desktop/myRAST-out"


setwd(pS)

##126 Genomes anotated by Malia and David
#List of MyRAST output files
desc<-read.table("description.txt",header=T)
#Genome available information (phylogenetic order and clade according to rpoB gene phylogeny, ecology, genus and species)
ord<-read.table("order_phylo-update.txt",header=T)
ordo<-read.table("order_phylo-old.txt",header=T)

##Step 1: Link GFF files with genomes (DAVID AND MALIA anotations)

setwd(pO)
ANNx<-as.vector(desc[,1])

#Retrieve genome NCBI number for each annotation file
GENc<-as.vector(sapply(ANNx,function(x){
	print(x)
	ANNt<-read.delim(x)
	z<-paste(unique(c("none",
		intersect(ordo[,2],as.vector(unique(ANNt[,1]))),
		intersect(ord[,2],as.vector(unique(ANNt[,1]))))),
		collapse="-")
	return(z)
}))

#Remove annotation files without NCBI number match (correspond to duplicates)
ANNx<-subset(ANNx,GENc!="none")
GENc<-subset(GENc,GENc!="none")

##Genome characteristics

PHYLOo<-c(as.vector(ord[,1]),seq(0,0,length.out=length(ordo[,1]))) #Phylogenetic order in rpoB tree
NCBIo<-c(as.vector(ord[,2]),as.vector(ordo[,2])) #NCBI number
STRAINo<-c(as.vector(ord[,3]),seq(0,0,length.out=length(ordo[,1]))) #Strain
GENo<-c(as.vector(ord[,5]),as.vector(ordo[,4])) #Genus
SPEo<-c(as.vector(ord[,6]),as.vector(ordo[,5])) #Species
ENVo<-c(as.vector(ord[,7]),as.vector(ordo[,6])) #Origin of strain isolation
CLAo<-c(as.vector(ord[,8]),as.vector(ordo[,7])) #Clade according to the rpoB phylogeny

#Metadata summary for each annotation table
NAMES<-t(sapply(GENc,function(x){
	nx<-setdiff(unlist(strsplit(x,split="-")),"none")
	phylox<-max(unique(unlist(sapply(nx,function(x){
		return(subset(PHYLOo, NCBIo ==x))
	}))))
	strainx<-c(setdiff(unique(unlist(sapply(nx,function(x){
		return(subset(STRAINo, NCBIo ==x))
	}))),"0"),"unknown")[1]
	genx<-c(setdiff(unique(unlist(sapply(nx,function(x){
		return(subset(GENo, NCBIo ==x))
	}))),"0"),"unknown")[1]	
	spex<-c(setdiff(unique(unlist(sapply(nx,function(x){
		return(subset(SPEo, NCBIo ==x))
	}))),"0"),"unknown")[1]	
	envx<-c(setdiff(unique(unlist(sapply(nx,function(x){
		return(subset(ENVo, NCBIo ==x))
	}))),"0"),"unknown")[1]	
	clax<-c(setdiff(unique(unlist(sapply(nx,function(x){
		return(subset(CLAo, NCBIo ==x))
	}))),"0"),"unknown")[1]					
	return(c(phylox, strainx, genx, spex, envx, clax))
}))

colnames(NAMES)<-c("Phylo-order","Strain","Genus","species","Environment","Clade")

#Sort by phylogenetic order for complete rpoB sequence phylogeny
ANNx<-ANNx[order(as.numeric(NAMES[,1]))]
NAMES<-NAMES[order(as.numeric(NAMES[,1])),]

#Strains for which genome annotations are missing
MISst<-sort(setdiff(ord[,3],NAMES[,2]))

MISord<-ord[as.vector(sapply(MISst,function(x){
	return(as.vector(subset(c(1:length(ord[,1])),as.vector(ord[,3])==x))[1])
})),]


setwd(pR)
write.table(MISord,"Additional-genomes.txt")

##Step 2: Link GFF files with genomes (JB annotations): 63 genomes missing in M&D annotations (mostly outgroup Microvirga) and 24 annitional genomes from Kembel lab assembled, filtered and annotated by JB (see script assemblies.R)

setwd(pS)
ordn<-read.table("Additional-genomes.txt",header=T)
descn <-read.table("description-additionnal.txt",header=T)
#work directories for genome fasta files
pF="/Users/jean-baptisteleducq/Desktop/myRAST-additional/fasta"


#Scaffold names in GFF files
setwd(pO)

ANNn<-as.vector(descn[,1])

Ngff<-sapply(ANNn,function(x){
	print(x)
	ANNt<-read.delim(x)
	return(paste(as.vector(unique(ANNt[,1])),collapse="=="))
	
})

#Find connection with scaffold names in fasta file
setwd(pF)

MATCH<-sapply(c(1:length(ordn[,1])),function(x){
	print(x)
	#Genome scaffold(s)
	Gx<-unlist(strsplit(as.vector(ordn[x,4]),split=".gz"))
	SEQx<-read.fasta(Gx)
	Nx<-names(SEQx)

	MATCHx<-sapply(Ngff,function(x){
		return(length(intersect(Nx,unlist(strsplit(x,split="==")))))
	})
	return(subset(ANNn, MATCHx==max(MATCHx)))
	
})

##Genome characteristics

PHYLOn<-as.vector(ordn[,1]) #Phylogenetic order in rpoB tree
NCBIn<-as.vector(ordn[,2]) #NCBI number
STRAINn<-as.vector(ordn[,3]) #Strain
GENn<-as.vector(ordn[,5]) #Genus
SPEn<-as.vector(ordn[,6]) #Species
ENVn<-as.vector(ordn[,7]) #Origin of strain isolation
CLAn<-as.vector(ordn[,8]) #Clade according to the rpoB phylogeny

#Metadata summary for each annotation table
NAMESn<-cbind(PHYLOn, STRAINn, GENn, SPEn, ENVn, CLAn)

colnames(NAMESn)<-c("Phylo-order","Strain","Genus","species","Environment","Clade")

##Merge data
NAMES<-rbind(NAMES,NAMESn)
ANNx<-c(ANNx, MATCH)
rownames(NAMES)<-ANNx

#remove 6666666.688705.txt (two microvirga genomes in the same file)

NAMES<-subset(NAMES,NAMES[,2]!="R24845")

#Remove "unknown" (duplicated CBMB27)
NAMES<-subset(NAMES,NAMES[,2]!="unknown")

#Remove duplicate genome names
NAMEi<-as.vector(sapply(unique(NAMES[,2]),function(x){
	return(subset(c(1:length(NAMES[,1])),NAMES[,2]==x)[1])
}))

NAMES<-NAMES[NAMEi,]

setwd(pR)

write.table(NAMES,"Annotated-genomes.txt")

##ARRIVÉ ICI 22092021

##Color for Methylobacterium clades (Outgroups in white)
COLcla<-ifelse(NAMES[,6]=="A1","blue",
	ifelse(NAMES[,6]=="A2","cyan",
	ifelse(NAMES[,6]=="A3","blue3",
	ifelse(NAMES[,6]=="A4","brown",
	ifelse(NAMES[,6]=="A5","yellow2",
	ifelse(NAMES[,6]=="A6","orange",
	ifelse(NAMES[,6]=="A7","pink",
	ifelse(NAMES[,6]=="A8","purple",
	ifelse(NAMES[,6]=="A9","red",
	ifelse(NAMES[,6]=="A10","green2",
	ifelse(NAMES[,6]=="B","black",
	ifelse(NAMES[,6]=="C","grey",
	ifelse(NAMES[,6]=="Microvirga","white",
	ifelse(NAMES[,6]=="Enterovirga",rgb(0.9,0.9,0.9),
	"white"))))))))))))))

##List of unique gene names

setwd(pO)
Fu<-sort(unique(unlist(sapply(ANNx,function(x){
	print(x)
	ANNt<-read.delim(x)
	Fx<-as.vector(ANNt[,8])
	return(unique(Fx))
}))))

#Remove genomes with low assembly quality (based on the number of scaffolds)
ANNx<-rownames(NAMES)

setwd(pO)
NSCA<-t(sapply(ANNx,function(x){
	print(x)
	ANNt<-read.delim(x)
	#scaffold size based on annotation
	Sx<-sort(sapply(unique(ANNt[,1]),function(x){
		return(max(subset(ANNt[,c(5:6)],ANNt[,1]==x)))
	}),decreasing=T)
	#N50
	N50<-subset(Sx,ifelse(cumsum(Sx)>=(sum(Sx)/2),1,0)==1)[1]
	N80<-subset(Sx,ifelse(cumsum(Sx)>=(0.8*sum(Sx)),1,0)==1)[1]
	N20<-subset(Sx,ifelse(cumsum(Sx)>=(0.2*sum(Sx)),1,0)==1)[1]
	#return number of scaffolds, number of annotations, number of unique annotations, total genome size, and N50
	return(c(
	length(unique(ANNt[,1])),
	length(ANNt[,8]),
	length(unique(ANNt[,8])),
	sum(Sx),N50,N80,N20
	))
}))

colnames(NSCA)<-c("Scaffolds","Annotations","Annot(unique)","GenomeSize","N50","N80","N20")

setwd(pO)

pdf("Methylo-Summary-genome-assembly.pdf",height=5,width=5)

plot(as.factor(NAMES[,6]), NSCA[,4]/10^6,las=2,cex.axis=0.8,ylab="Genome size (Mb)",xlab="Clades")

plot(as.factor(NAMES[,6]), NSCA[,3],las=2,cex.axis=0.8,ylab="Unique annotations",xlab="Clades")



plot(NSCA[,4]/10^6, NSCA[,3],log="xy",las=1,cex.axis=0.8,xlab="Genome size (Mb)",ylab="Unique annotations",bg= COLcla,pch=21,cex=1)


plot(NSCA[,4]/10^6, NSCA[,5]/10^6,log="y",las=1,cex.axis=0.8,xlab="Genome size (Mb)",ylab="N50 (Mb)",bg= COLcla,pch=21,cex=0.8)

plot(NSCA[,2]/NSCA[,1],NSCA[,5]/10^6,log="xy",bg= COLcla,pch=21,cex=0.8,las=1,cex.axis=0.8,xlab="Gene per scaffold",ylab="N50 (Mb)")

dev.off()

#Remove genomes with low N50 and low number of core genes (determined later in SINGcountNEW)

#NAMES<-subset(NAMES,NSCA[,5]>=10000)
LOW<-c("WL116","WL103","WL18","B1","B34","CP3")

#Remove genomes that have a high proportion of Gapped alignments (Determined later in step 9), probably low quality assembly

OUTg<-c("WSM2598", "MAMP4754","DSM16961","DSM13060","TK0001","NCCP-1258","JC119","MGYG-HGUT-02310","DSM16371","WL7","WL1","c27j1","CDVBN77","SE3.6","DSM25844","JCM14648","SR1.6/6")

NAMES<-subset(NAMES,sapply(NAMES[,2],function(x){
	return(length(intersect(x,c(LOW, OUTg))))
})==0)

setwd(pS)
write.table(NAMES,"Genome-description.txt")

###############
###############
###############
##Step 3: Calculate gene occurence per genome

#Number of scafolds with annotation per genomes and maximum scaffold size

library(seqinr)
library(gplots)
library(ade4)

#work directory for analyzes
pS="/Users/jean-baptisteleducq/Dropbox/PostDoc-Moscow-ID/Project/myRAST-gene_annotation"
#work directory for annotation files (GFF3 myRAST outputs)
pO="/Users/jean-baptisteleducq/Desktop/myRAST-out"

setwd(pS)
NAMES<-read.table("Genome-description.txt",header=T)
ANNx<-rownames(NAMES)

setwd(pO)
SCAc<-t(sapply(ANNx,function(x){
	print(x)
	ANNt<-read.delim(x)
	#Annotations
	Fx<-summary(as.factor(ANNt[,8]),maxsum=length(unique(ANNt[,8])))
	#Scaffolds
	Sx<-as.vector(unique(ANNt[,1]))
	
	#Number of scaffolds
		Ns=length(unique(Sx))
	#Number of annotations
		Na=sum(Fx)
	#Number of hypothetical proteins, repeat regions and Mobile element proteins
		Nh=subset(Fx, names(Fx)=="hypothetical protein")
		Nr=subset(Fx, names(Fx)=="repeat region")
		Nm=subset(Fx, names(Fx)=="Mobile element protein")
	#Number of unique annotations, excluding hypothetical proteins, repeat regions and Mobile element proteins
		Nu=length(setdiff(names(Fx),c("hypothetical protein","repeat region","Mobile element protein")))
	#Average number of copies per annotation, excluding hypothetical proteins, repeat regions and Mobile element proteins
		Fx<-subset(Fx,names(Fx)!="hypothetical protein")
		Fx<-subset(Fx,names(Fx)!="repeat region")
		Fx<-subset(Fx,names(Fx)!="Mobile element protein")
		Mx=round(mean(Fx),3)
		SDx=round(sd(Fx),3)
	#Scaffold size
	SIZE<-sort(as.vector(sapply(Sx,function(x){
		return(max(subset(ANNt[,c(5:6)],ANNt[,1]==x)))
	})),decreasing=T)
	#Total size
	Tx=sum(SIZE)
	#N50
	N50<-subset(SIZE,cumsum(SIZE)>=0.5* Tx)[1]
	
	return(c(x,Ns,Tx,N50, Na,Nh,Nr,Nm,Nu,Mx,SDx))
		
}))

colnames(SCAc)<-c("GFF_file","Scaffolds","Total_Size","N50","Annotations","hypoth.","repeat","Mobile","other_unique_annotation","average_copy(other_unique)","sd_copy(other_unique)")

setwd(pS)

write.table(cbind(NAMES,SCAc),"Annotation-Assembly-summary.txt",row.names=F)

COUNTdraft<-sapply(sort(unique(NAMES[,3])),function(x){
	SCAx<-subset(ifelse(as.numeric(SCAc[,4])>3*10^6,"complete","draft"),NAMES[,3]==x)
	return(sapply(c("complete","draft"),function(x){
		return(length(subset(SCAx, SCAx==x)))
	}))
})
colnames(COUNTdraft)<-sort(unique(NAMES[,3]))

barplot(COUNTdraft,las=1,cex.axis=0.8,ylab="genomes",legend=T)

##List of unique gene names

setwd(pO)
Fu<-sort(unique(unlist(sapply(ANNx,function(x){
	print(x)
	ANNt<-read.delim(x)
	Fx<-as.vector(ANNt[,8])
	return(unique(Fx))
}))))

#Remove hypothetical proteins, repeat regions and Mobile element proteins
Fu<-setdiff(Fu,c("hypothetical protein","repeat region","Mobile element protein"))

#Retrieve the abundance of each gene per genome
Fall<-sapply(ANNx,function(x){
	print(x)
	ANNt<-read.delim(x)
	Fx<-as.vector(ANNt[,8])	
	return(sapply(Fu,function(x){
		return(length(subset(Fx,Fx==x)))
	}))	
})

#Color clade for Methylobacterium clades (Outgroups in white)
COLcla<-ifelse(NAMES[,6]=="A1","blue",
	ifelse(NAMES[,6]=="A2","cyan",
	ifelse(NAMES[,6]=="A3","blue3",
	ifelse(NAMES[,6]=="A4","brown",
	ifelse(NAMES[,6]=="A5","yellow2",
	ifelse(NAMES[,6]=="A6","orange",
	ifelse(NAMES[,6]=="A7","pink",
	ifelse(NAMES[,6]=="A8","purple",
	ifelse(NAMES[,6]=="A9","red",
	ifelse(NAMES[,6]=="A10","green2",
	ifelse(NAMES[,6]=="B","black",
	ifelse(NAMES[,6]=="C","grey",
	ifelse(NAMES[,6]=="Microvirga","white",
	ifelse(NAMES[,6]=="Enterovirga",rgb(0.9,0.9,0.9),
	"white"))))))))))))))


setwd(pS)

#Matrix of gene occurence (row) per genome (column). Abundance was set as 2 when abundance >=2

pdf("Matrix-gene.pdf",height=20,width=25)

heatmap.2(
	ifelse(Fall>=2,2,Fall),scale="none",labCol=NAMES[,2],labRow="",cexCol=0.8,margins=c(10,2),trace="none",key=F,breaks=4,col=c("white",rgb(0,1,0,1),rgb(0,0.3,0,1)),ColSideColors= COLcla)

dev.off()



Fsum<-t(sapply(c(1:length(Fall[,1])),function(x){
	fx<-as.vector(Fall[x,])
	return(c(length(subset(fx,fx!=0)),
		mean(fx),
		min(fx),
		max(fx)))
}))
rownames(Fsum)<-Fu
rownames(Fall)<-Fu

#Sort genes according to occurence in genomes: 1:1 copies first

ORDO2<-order((abs(as.numeric(Fsum[,2])-1)),decreasing=F)

Fall<-Fall[ORDO2,]
Fsum <-Fsum[ORDO2,]
Fu<-Fu[ORDO2]

ORDO1<-order(as.numeric(Fsum[,1]),decreasing=T)

Fall<-Fall[ORDO1,]
Fsum <-Fsum[ORDO1,]
Fu<-Fu[ORDO1]

Fsum<-cbind(c(1:length(Fsum[,1])),Fsum)

colnames(Fsum)<-c("ann","occurence","mean","min","max")

write.table(Fall,"Matrix-gene.txt")
write.table(Fsum,"Summary-per-gene.txt")

#Display PCA
ACP=dudi.pca(t(Fall), scannf=F)
part<-round(100* ACP $eig/sum(ACP $eig),2)


x1=min(ACP $li[,1])
x2=max(ACP $li[,1])
y1=min(ACP $li[,2])
y2=max(ACP $li[,2])

CLAu<-sort(unique(NAMES[,6]))
COLclau <-sapply(CLAu,function(x){
	return(unique(subset(COLcla, NAMES[,6]==x)))
})

pdf("PCA-Gene-Copy-number-per-clade.pdf",height=5,width=5)
par(mar=c(4,4,1,1),bty="n")
plot(0,0,col=NA,ylim=c(y1,y2),xlim=c(x1,x2),
		las=1,cex.lab=1,cex.axis=0.8,font.lab=2,
		xlab=paste("Axis 1 (",part[1],"%)",sep=""),
		ylab=paste("Axis 2 (",part[2],"%)",sep=""),
		main="",cex.main=1)
		
points(ACP $li[,c(1,2)],pch=21,bg= COLcla,col="black")	

legend(x = "bottomright",cex=0.8, 
	horiz=F,legend= sort(unique(NAMES[,6])),
	text.col= "black",pt.bg=COLclau,
	box.col=NA,pch=21,border=NA,
	col="black",text.font=3,title="Clades",
	ncol=2)	

dev.off() 


#Step 4: Proportion of genes with 1,2,3,4 or more copies in fonction of the number of scaffolds 
setwd(pO)	

GeneCopy<-t(sapply(c(1:length(ANNx)),function(x){
	print(x)
	Gx=x
	ANNt<-read.delim(ANNx[Gx])
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

setwd(pS)

write.table(GeneCopy,"Gene-copy.txt")

setwd(pS)

pdf("GeneCopy-function-scaffold.pdf",height=5,width=4)

par(mar=c(4,4,1,1),bty="n")

plot(GeneCopy[,1], GeneCopy[,3],xlab="number of scaffolds",ylab="proportion of genes",las=1,cex.axis=0.8,log="x",ylim=c(0,1),pch=19,col=rgb(0,0,1,0.4),cex=0.5)

points(GeneCopy[,1], GeneCopy[,4],pch=19,col=rgb(0,1,0,0.4),cex=0.5)

points(GeneCopy[,1], GeneCopy[,5],pch=19,col=rgb(1,0.7,0,0.4),cex=0.5)

points(GeneCopy[,1], GeneCopy[,6],pch=19,col=rgb(1,0,0,0.4),cex=0.5)

points(GeneCopy[,1], GeneCopy[,7],pch=19,col=rgb(0.7,0,1,0.4),cex=0.5)

text(seq(5,5,length.out=5),seq(0.6,0.4,length.out=5),c("Single copy","2 copies","3 copies","4 copies","5 copies"),col=c("blue","green","orange","red","purple"),cex=0.8)

dev.off()



#Step 5: For each annotation, return the number of genomes with 0, 1 or more copies


triangle.plot <- function(x, ...) {
	cat<-colnames(x)
	names=rownames(x)
	col.cat <-c("orange2","blue","green2")
	n.axis=10
	col.names <-rgb(0,0,0,0.05)
	nd=n.axis+1	
	#graphique
	par(mar=c(1,1,1,1),bty="n")
	#Grid
	zv<-round(seq(0,1,length.out=nd),nchar(nd))
	plot(-100,-100,xlim=c(-0.2, 1.2),ylim=c(-0.2,1.2),xlab="",ylab="",cex.axis=0.8,xaxt="n",yaxt="n",col.lab="black")
	sapply(zv,function(x){
		X1=0.5*x
		X2=1-0.5*x
		segments(X1, x,X2, x,col= col.cat[1],lwd=0.5)
		segments(X1, x,x,0,col= col.cat[2],lwd=0.5)
		segments(X2, x,1-x,0,col= col.cat[3],lwd=0.5)
	})
	text(zv-0.01,-0.02,zv,cex=0.8,col= col.cat[3],srt=60, pos=2, offset=0)
	text(0.5*zv-0.01,zv+0.02,1-zv,cex=0.8,col= col.cat[2],srt=-60, pos=2, offset=0)
	text(1-0.5*zv+0.05,zv,zv,cex=0.8,col= col.cat[1])
	text(0.45,-0.15, cat[3],cex=0.8,col= col.cat[3],font=2)
	text(0.1,0.6, cat[2],cex=0.8,col= col.cat[2],font=2)
	text(0.9,0.55, cat[1],cex=0.8,col= col.cat[1],font=2)
	#Observed values
	Y0=x[,1]/apply(x,1,sum)
	X0=x[,3]/apply(x,1,sum)
	X0=X0+0.5*Y0
	par(new=T)
	plot(X0,Y0,xlim=c(-0.2, 1.2),ylim=c(-0.2,1.2),xlab="",ylab="",xaxt="n",yaxt="n",pch=19,col= col.names,lwd=1.5,cex=1)	
	text(X0,Y0+0.05, names,cex=0.6,col= col.names,font=1)
}


GeneCount<-t(sapply(c(1:length(Fall[,1])),function(x){
	Fx<-as.vector(Fall[x,])
	return(c(length(subset(Fx,Fx==0)),
	length(subset(Fx,Fx==1)),
	length(subset(Fx,Fx>1))))
}))

colnames(GeneCount)<-c("Absent","1 copy","> 1 copies")

write.table(GeneCount,"GeneCount.txt",row.names=F)

pdf("GeneCount-TrianglePlot.pdf")
triangle.plot(GeneCount)
dev.off()


##First option: Define candidate core genes as genes that are present in one copy in at least pG% of genomes

pG=0.9

setwd(pS)

Fcore<-subset(c(1:length(Fu)), GeneCount[,2]>=(pG*length(ANNx)))

pdf("Matrix-gene-core.pdf",height=20,width=25)

heatmap.2(
	ifelse(Fall[Fcore,]>=2,2, Fall[Fcore,]),scale="none",labCol=NAMES[,2],labRow="",cexCol=0.8,margins=c(10,2),trace="none",key=F,breaks=4,col=c("white",rgb(0,1,0,1),rgb(0,0.3,0,1)),ColSideColors= COLcla)

dev.off()


#Second option: Define candidate core genes that are present in one copy in all complete genomes (defined as genomes with N50 > 3.10^6)

COMP<-ifelse(as.numeric(SCAc[,4])>3*10^6,"complete","draft")
COMPi<-subset(c(1:length(COMP)),COMP=="complete")

GeneCount2<-t(sapply(c(1:length(Fall[,1])),function(x){
	Fx<-as.vector(Fall[x, COMPi])
	return(c(length(subset(Fx,Fx==0)),
	length(subset(Fx,Fx==1)),
	length(subset(Fx,Fx>1))))
}))

Fcore2<-subset(c(1:length(Fu)), GeneCount2[,2]==length(COMPi))

pdf("GeneCount-Complete-genomes-TrianglePlot.pdf")
triangle.plot(GeneCount2)
dev.off()


###OPTIONAL
#Fcore<-intersect(Fcore,Fcore2) #The intersect between Fcore and Fcore2 is almost the same than Fcore2 alone, suggesting that the analysis of only complete genomes can already give a robust idea of the core genome and how the pangenome is distributed withing Methylobacteriaceae. We will use this conservative list of gene to study the core genome
###

#Step 6: For each genome produce a fasta file only with core genes (Fcore = intersect among all complete genomes)

#At this point, we selected genes that are present in one copy in all complete genomes. But some of these genes were absent, fragmented and/or duplicated in some draft genomes. We thus need to define criteria to keep or reject sequences before generating alignments. 

#Simflify gene names

ng=nchar(length(Fu))

FuMOD<-sapply(c(1:length(Fu)),function(x){
	paste(c("Gene",seq(0,0,length.out=(ng-nchar(x))),x),collapse="")
})

write.table(cbind(FuMOD,Fu),"Gene-modified-names.txt",row.names=F,col.names=c("FuMOD","Fu"))

Ax<-c(1:length(ANNx))

#Number of copies per gene per genome
FcoreCount<-sapply(Ax,function(x){
	Gx=x
	setwd(pO)
	ANNt<-read.delim(ANNx[Gx])
	ANNtx<-unlist(strsplit(ANNx[Gx],split=".txt"))[1]
	NAMEx<-as.vector(NAMES[Gx,2])
	print(NAMEx)
	return(
	as.vector(sapply(Fu[Fcore],function(x){
		return(length(subset(c(1:length(ANNt[,1])),as.vector(ANNt[,8])==x)))
	}))
	)
})

heatmap.2(
	ifelse(FcoreCount>=2,2, FcoreCount),scale="none",labCol=NAMES[,2],labRow="",cexCol=0.8,margins=c(10,2),trace="none",key=F,breaks=4,col=c("white",rgb(0,1,0,1),rgb(0,0.3,0,1)),ColSideColors= COLcla)

# Average size of each gene among complete genomes only

FcoreSizeComp<-sapply(COMPi,function(x){
	Gx=x
	setwd(pO)
	ANNt<-read.delim(ANNx[Gx])
	ANNtx<-unlist(strsplit(ANNx[Gx],split=".txt"))[1]
	NAMEx<-as.vector(NAMES[Gx,2])
	print(NAMEx)
	return(as.numeric(unlist(sapply(Fu[Fcore],function(x){
		return(c(nchar(subset(as.vector(ANNt[,12]),as.vector(ANNt[,8])==x)),NA)[1])
	}))))
})

AVEcomp<-apply(FcoreSizeComp,1,mean,na.rm=T)

#Average size of sequences from a same gene, normalized by the average size of single copies found in complete genomes)
FcoreSize<-sapply(Ax,function(x){
	Gx=x
	setwd(pO)
	ANNt<-read.delim(ANNx[Gx])
	ANNtx<-unlist(strsplit(ANNx[Gx],split=".txt"))[1]
	NAMEx<-as.vector(NAMES[Gx,2])
	print(NAMEx)
	
	z<-as.vector(sapply(Fu[Fcore],function(x){		
		mean(nchar(subset(as.vector(ANNt[,12]),
			as.vector(ANNt[,8])==x)),na.rm=T)	
	}))/AVEcomp
	
	#nz<-as.vector(sapply(Fu[Fcore],function(x){	#	
	#	length(subset(as.vector(ANNt[,12]),
	#		as.vector(ANNt[,8])==x))	
	#}))
	
	return(z)
})


write.table(FcoreCount,"FcoreCount.txt",row.names=F,col.names=F)
write.table(FcoreSize,"FcoreSize.txt",row.names=F,col.names=F)


vx<-seq(-0.25,0.25,length.out=101)

setwd(pS)

pdf("NormalizedSize-vs-CopyNumber.pdf",width=4,height=5)

plot(-10,-10,xlim=c(0,max(FcoreCount)),
	ylim=c(0,max(FcoreSize,na.rm=T)),
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
	
	points(FcoreCount[x,]+sample(vx,1),FcoreSize[x,],pch=19,cex=0.1,col=rgb(0.7,0.7,0.7,0.1))
})


legend(x = "topright",cex=0.8, horiz=F,legend= paste("n=",c(1:4),sep=""),text.col= "black",box.col=NA,col= c("black"),border=NA,text.font=1,lty=c(1,2,3,3),title="Real copy number",pt.lwd=c(2,2,1.5,1))

dev.off()

#For each gene, count the number of genomes with at least 2 copies with an average relative size of at least 0.75  (potentially true duplicated gene if counted in at least one genome)

DUP<-apply(ifelse(FcoreCount<2,0,ifelse(is.na(FcoreSize)==T,0,ifelse(FcoreSize>0.75,1,0))),1,sum) 

#Remove potentially duplicated genes in at least one genome (check that "DNA-directed RNA polymerase beta subunit (EC 2.7.7.6)"==rpoB is conserved after this step)
FcoreNEW<-subset(Fcore, DUP<1)
FcoreCountNEW<-subset(FcoreCount, DUP<1)
FcoreSizeNEW<-subset(FcoreSize, DUP<1)

#plot(-10,-10,xlim=c(0,max(FcoreCountNEW)),ylim=c(0,max(FcoreSizeNEW,na.rm=T)),las=1,cex.axis=0.8,xlab="copy number",ylab="normalized size",main=paste("After filtering (",length(FcoreNEW)," genes)",sep=""))
sapply(c(1:length(FcoreNEW)),function(x){
	
	points(FcoreCountNEW[x,]+sample(vx,1), FcoreSizeNEW[x,],pch=19,cex=0.1,col=rgb(0,0,0,0.1))
})


dev.off()

#For each remaining gene and each genome, only keep single copies that have a relative size of more than 0.7 and less than 1.3 

SINGcountNEW<-ifelse(is.na(FcoreSizeNEW)==T,0,ifelse(FcoreSizeNEW<0.7,0,ifelse(FcoreSizeNEW>1.3,0,1)))

SINGcountNEW<-ifelse(FcoreCountNEW==1, SINGcountNEW,0)

pdf("Gene-occurence-after-filtering.pdf",width=5,height=5)

Ng=181

plot(as.numeric(SCAc[,4])/10^6,apply(SINGcountNEW,2,sum),log="x",las=1,cex.axis=0.8,xlab="N50 (Mb)",ylab="Core genes",bg=COLcla,pch=21,main=paste("Before removing scattered genes (",length(FcoreNEW),")",sep=""))

hist(apply(SINGcountNEW,1,sum),breaks=30,las=1,cex.axis=0.8,xlab="Occurence in genomes",ylab="Gene count",main="Distribution of single copies after filtering",cex.main=0.8)

segments(Ng,0, Ng,1000,lwd=2,lty=2,col="red")

#heatmap.2(subset(SINGcountNEW,apply(SINGcountNEW,1,sum)>=175),scale="none",labCol=NAMES[,2],labRow="",cexCol=0.8,margins=c(10,2),trace="none",key=F,col=c("white",rgb(0,1,0,1)),ColSideColors= COLcla)

#Remove genes for which we have a copy in less than Ng genomes


FcoreNEW<-subset(FcoreNEW, apply(SINGcountNEW,1,sum)>= Ng)
FcoreCountNEW<-subset(FcoreCountNEW, apply(SINGcountNEW,1,sum)>= Ng)
FcoreSizeNEW<-subset(FcoreSizeNEW, apply(SINGcountNEW,1,sum)>= Ng)
SINGcountNEW<-subset(SINGcountNEW, apply(SINGcountNEW,1,sum)>= Ng)

plot(as.numeric(SCAc[,4])/10^6,apply(SINGcountNEW,2,sum),log="x",las=1,cex.axis=0.8,xlab="N50 (Mb)",ylab="Core genes",bg=COLcla,pch=21,main=paste("After removing scattered genes (",length(FcoreNEW),")",sep=""))

dev.off()

#########For each genome, extract the nucleotide sequence of single copy genes present in at least Ng genomes and that have a relative sequence size in the range 0.7-1.3 from GFF annotation file


NcorePerGenome<-sapply(Ax,function(x){
	Gx=x
	setwd(pO)
	ANNt<-read.delim(ANNx[Gx])
	ANNtx<-unlist(strsplit(ANNx[Gx],split=".txt"))[1]
	NAMEx<-as.vector(NAMES[Gx,2])
	#Only keep Fcore genes that are present in this genome
	Fx=Fu [subset(FcoreNEW ,SINGcountNEW[,Gx]==1)]
	#Retrieve modified gene names (without space)
	FxMOD<-as.vector(sapply(Fx,function(x){
		return(subset(FuMOD,Fu==x))
	}))
	SEQx<-sapply(Fx,function(x){
		GN=x
		return(list(unlist(strsplit(subset(as.vector(ANNt[,12]),
			as.vector(ANNt[,8])==GN),split=""))))
	})		
	setwd(paste(pO,"/OUT-FASTA",sep=""))
	write.fasta(SEQx,names= paste("Strain",NAMEx,"GFF" ,ANNtx, "Gene", FxMOD,sep="=="),paste("Core-Genes-", ANNtx,".fas",sep=""))
	print(paste(c("Strain",NAMEx,":",length(Fx),"core genes",";",length(SEQx),"sequences"),collapse=" "))
	return(length(Fx))
})

setwd(pS)
pdf("CoreGene-per-genome-N50.pdf",width=5,height=5)

plot(as.numeric(SCAc[,4]), NcorePerGenome,pch=1,col= COLcla,log="x",las=1,cex.axis=0.8,ylab="Core gene sequences retrieved per genome",xlab="N50 (pb, log scale)",cex=1)

dev.off()


#Step 7: For each core gene, produce an alignement with genomes for which sequences are available

setwd(paste(pO,"/OUT-FASTA",sep=""))

FcoreCountGen<-sapply(FcoreNEW,function(x){
	Fx= FuMOD[x]	
	
	#count the number of the gene occurence per genome
	Nx=as.vector(subset(SINGcountNEW, FcoreNEW== x))
	
	#remove genomes where the gene is absent
	ANx<-subset(Ax,Nx==1)
	SEQ<-sapply(ANx,function(x){
		#print(x)
		Gx=x
		ANNtx<-unlist(strsplit(ANNx[Gx],split=".txt"))[1]
		NAMEx<-as.vector(NAMES[Gx,2])
		SEQx<-read.fasta(paste("Core-Genes-", ANNtx,".fas",sep=""))
		NS=paste("Strain",NAMEx,"GFF" ,ANNtx, "Gene",Fx,sep="==")
		return(subset(SEQx,names(SEQx)==NS))
	})
	write.fasta(SEQ,paste(Fx,".fas",sep=""),names=as.vector(NAMES[ANx,2]))
	print(paste(c("Gen",Fx,": in",sum(Nx),"genomes"),collapse=" "))
	return(sum(Nx))

})

#Genes manually modified: 
# Gene0006 (17SD2-17 manually added from genome)
# Gene0032 (NI91,  CLZ, DM1 manually retrieved from annotation file)
# Gene0043 (WL9 manually retrieved from genome)
# Gene0050 (L1A1 manually retrieved from annotation file)
# Gene0057 (DSM25903 and DB1703  manually retrieved from annotation file)

AVEcompNew<-sapply(FcoreNEW,function(x){
	return(subset(AVEcomp,Fcore==x))
})

FcoreCountGen <-as.data.frame(cbind(Fu[FcoreNEW],FuMOD[FcoreNEW], FcoreCountGen, AVEcompNew))

colnames(FcoreCountGen)<-c("Annotation","Fasta-name","Occurence","AverageSize(bp)")

setwd(pS)

write.table(FcoreCountGen,"Summary-per-core-gene.txt",row.names=F)

###############
###Step 8: align sequence for each core gene


#work directory for analyzes
pS="/Users/jean-baptisteleducq/Dropbox/PostDoc-Moscow-ID/Project/myRAST-gene_annotation"
#work directory for annotation files (GFF3 myRAST outputs)
pO="/Users/jean-baptisteleducq/Desktop/myRAST-out"


setwd(pS)


FcoreNEW<-as.vector(read.table("Summary-per-core-gene.txt",header=T))[,2]


library(msa)
library(seqinr)


setwd(paste(pO,"/OUT-FASTA",sep=""))


sapply(FcoreNEW,function(x){
	print(x)
	Fx= x

	#Get the unaligned nucleotide sequence
	SEQ<-read.fasta(paste(Fx,".fas",sep=""))

	#Sequence names
	NX=names(SEQ)	
	
	#Sequence sizes
	Sx<-as.numeric(summary(SEQ)[,1])
	
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

	#SEQmuscle<-msa(paste(Fx,".fas",sep=""),method="Muscle",type="dna")

	#Align AA sequences together using ClustalW - IMPORTANT: keep same order as input, since sequence names are removed during analysis, and retrieve sequence names in the original unaligned fasta file
	#SEQw<-msa(paste(Fx,".fas",sep=""),method="ClustalW",type="dna",order="input")

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
	
	Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
	
	write.fasta(SEQn, names=NX,
		file.out =paste(Fx,"_aligned.fas",sep=""))
})


#Review each alignment and replace "NNN(---)NNN" motifs by "NNN(NNN)NNN"

setwd(paste(pO,"/OUT-FASTA",sep=""))

sapply(FcoreNEW,function(x){
	print(x)
	Fx= x
	Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
	#Get the aligned nucleotide sequence
	SEQ<-read.fasta(paste(Fx,"_aligned.fas",sep=""))
	
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
		
		#positions tu replace by 'Ns'
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
	
	write.fasta(SEQn,paste(Fx,"_alignedN.fas",sep=""),names=names(SEQ))
	
})	

#Check that gap replacement by Ns did not affected alignments

setwd(paste(pO,"/OUT-FASTA",sep=""))

#Ony "n" should be returned
sum<-sapply(FcoreNEW,function(x){
	Fx= x
	Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
	#Get the aligned nucleotide sequence
	SEQ<-read.fasta(paste(Fx,"_aligned.fas",sep=""))
	SEQn<-read.fasta(paste(Fx,"_alignedN.fas",sep=""))
	sumx<-unlist(sapply(c(1:length(SEQ)),function(x){
		seq<-as.vector(SEQ[[x]])
		seqn<-as.vector(SEQn[[x]])
		return(subset(seqn,seq!=seqn))
	}))
	print(c(x,unique(sumx))	)
	return(sumx)
})

unique(unlist(sum))

#Remove positions with more than ppn% of missing data
ppn=0.9

setwd(paste(pO,"/OUT-FASTA",sep=""))

#Ony "n" should be returned
sumNEW<-t(sapply(FcoreNEW,function(x){
	print(x)
	Fx= x
	Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
	#Get the aligned nucleotide sequence
	SEQn<-read.fasta(paste(Fx,"_alignedN.fas",sep=""))
	NX<-names(SEQn)
	#Convert in matrix
	SEQn<-sapply(c(1:length(SEQn)),function(x){
		return(as.vector(SEQn[[x]]))
	})
	#Count the proportion of N per position
	pn<-sapply(c(1:length(SEQn[,1])),function(x){
		length(subset(SEQn[x,],SEQn[x,]=="n"))
	})/length(NX)
	
	SEQn<-subset(SEQn, pn<(1-ppn))
	
	SEQn<-sapply(c(1:length(NX)),function(x){
		return(list(SEQn[,x]))
	})
	
	write.fasta(SEQn,paste(Fx,
		"_alignedNF.fas",sep=""),names=NX)
	return(c(x,Fx,length(SEQn),length(SEQn[[1]])))	
		
}))	

colnames(sumNEW)<-c("Name1","Name2","Sequences","Size(bp)")
write.table(sumNEW,"Summary-per-gene-after-alignment.txt",row.names=F)


##Step 9: #Review each alignment and identify strains that decrease significantly the number of comparable sites (incomplete sequences, wrong alignements): RAxML and most phylogenetic approaches using ML do not take into account gaps and missing data (ns), hence decreasing dramatically the information available in an alignment. For each alignment, count the number of Ns and gaps overall and induced by each strain

library(ade4)
library(gplots)
library(seqinr)

#work directory for analyzes
pS="/Users/jean-baptisteleducq/Dropbox/PostDoc-Moscow-ID/Project/myRAST-gene_annotation"
#work directory for annotation files (GFF3 myRAST outputs)
pO="/Users/jean-baptisteleducq/Desktop/myRAST-out"

setwd(pS)

SUMgen<-read.table("Summary-per-core-gene.txt",header=T)
FcoreNEW<-as.vector(SUMgen[,2])
ANNx<-as.vector(SUMgen[,1])
NAMES<-read.table("Annotation-Assembly-summary.txt",header=T)


#Proportion of Ns per strain and per gene
setwd(paste(pO,"/OUT-FASTA",sep=""))

OUTs<-c("")

OUTg<-c("")

#Strain that should not be removed because representing rare clades
noOUTs<-c(
	subset(NAMES[,2],NAMES[,6]=="A7"),
	subset(NAMES[,2],NAMES[,6]=="A8"),
	subset(NAMES[,2],NAMES[,6]=="A10"),
	subset(NAMES[,2],NAMES[,6]=="Enterovirga")
)

OUTs<-setdiff(OUTs,noOUTs)

#clades
Cx=sapply(setdiff(NAMES[,2], OUTs),function(x){
	subset(NAMES[,6],NAMES[,2]==x)
})

PNs<-sapply(setdiff(FcoreNEW, OUTg),function(x){
	print(x)
	Ax<-subset(ANNx, FcoreNEW==x)
	Fx= x
	Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
	#Get the aligned nucleotide sequence
	SEQ<-read.fasta(paste(Fx,"_alignedNF.fas",sep=""))
	#Strain names
	Nx<-names(SEQ)
		
	#Convert in matrix
	SEQ<-sapply(c(1:length(SEQ)),function(x){
		return(as.vector(SEQ[[x]]))
	})
		
	#Proportion of Ns positions per strain
	PNx<-sapply(c(1:length(SEQ[1,])),
		function(x){
			return(length(subset(SEQ[,x],
			SEQ[,x]=='n')))
		})/length(SEQ[,1])
	
	#Complete missing gene (replace by 1)
	PNx<-sapply(setdiff(NAMES[,2], OUTs),function(x){
		return(c(subset(PNx,Nx==x),1)[1])
	})
	
	#Normalize values per clade
	PNx<-sapply(c(1:length(PNx)),function(x){
		return(PNx[x]/mean(subset(PNx, Cx==Cx[x])))
	})
	
	return(ifelse(is.nan(PNx)==T,0, PNx))
	
})	


#Display PCA

par(mar=c(4,4,1,1),bty="n",mfrow=c(1,2))

ACP=dudi.pca(PNs, scannf=F)
part<-round(100* ACP $eig/sum(ACP $eig),2)

x1=min(ACP $li[,1])*1.1
x2=max(ACP $li[,1])*1.1

y1=min(ACP $li[,2])*1.1
y2=max(ACP $li[,2])*1.1

par(mar=c(4,4,1,1),bty="n")
plot(0,0,col=NA,ylim=c(y1,y2),xlim=c(x1,x2),
		las=1,cex.lab=1,cex.axis=0.8,font.lab=2,
		xlab=paste("Axis 1 (",part[1],"%)",sep=""),
		ylab=paste("Axis 2 (",part[2],"%)",sep=""),
	cex.main=1,main="Per strain")		
text(ACP $li[,c(1,2)],col="black",cex=0.5,labels= rownames(PNs),font=1)

ACP=dudi.pca(t(PNs), scannf=F)
part<-round(100* ACP $eig/sum(ACP $eig),2)

x1=min(ACP $li[,1])*1.1
x2=max(ACP $li[,1])*1.1

y1=min(ACP $li[,2])*1.1
y2=max(ACP $li[,2])*1.1

par(mar=c(4,4,1,1),bty="n")
plot(0,0,col=NA,ylim=c(y1,y2),xlim=c(x1,x2),
		las=1,cex.lab=1,cex.axis=0.8,font.lab=2,
		xlab=paste("Axis 1 (",part[1],"%)",sep=""),
		ylab=paste("Axis 2 (",part[2],"%)",sep=""),
		cex.main=1,main="Per gene")		

text(ACP $li[,c(1,2)],col="black",cex=0.5,labels= colnames(PNs))		

####Display alignment and produce fasta files for each gene after removing outlier genes and genomes

COLcla<-ifelse(NAMES[,6]=="A1","blue",
	ifelse(NAMES[,6]=="A2","cyan",
	ifelse(NAMES[,6]=="A3","blue3",
	ifelse(NAMES[,6]=="A4","brown",
	ifelse(NAMES[,6]=="A5","yellow2",
	ifelse(NAMES[,6]=="A6","orange",
	ifelse(NAMES[,6]=="A7","pink",
	ifelse(NAMES[,6]=="A8","purple",
	ifelse(NAMES[,6]=="A9","red",
	ifelse(NAMES[,6]=="A10","green2",
	ifelse(NAMES[,6]=="B","black",
	ifelse(NAMES[,6]=="C","grey",
	ifelse(NAMES[,6]=="Microvirga","white",
	ifelse(NAMES[,6]=="Enterovirga",rgb(0.9,0.9,0.9),
	"white"))))))))))))))

setwd(paste(pO,"/OUT-FASTA",sep=""))

GapNs<-t(sapply(
	setdiff(FcoreNEW,OUTg),function(x){
	print(x)
	Ax<-subset(ANNx, FcoreNEW==x)
	Fx= x
	Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
	#Get the aligned nucleotide sequence
	SEQ<-read.fasta(paste(Fx,"_alignedNF.fas",sep=""))
	#Strain names
	Nx<-names(SEQ)
	
	#Remove outliers
	#keepx<-sapply(Nx,function(x){
 	#	return(length(
 	#		intersect(x,
 	#		setdiff(NAMES[,2], OUTs))))
 	#})
 	#SEQ<-subset(SEQ, keepx==1)
 	#Nx<-subset(Nx,keepx==1)		
	
	#write.fasta(SEQ,names=Nx,
	#	paste(Fx,"_alignedNF.fas",sep=""))
		
	#Convert in matrix and display polymorphism
	SEQ<-sapply(c(1:length(SEQ)),function(x){
		return(as.vector(SEQ[[x]]))
	})
	
	cx<-as.vector(sapply(Nx,function(x){
			return(subset(NAMES[,6],NAMES[,2]==x))
		}))
	
	colc<-as.vector(sapply(Nx,function(x){
			return(subset(COLcla,NAMES[,2]==x))
		}))
	
	SEQi<-t(as.matrix(ifelse(SEQ=="n",NA,
		ifelse(SEQ=="t",1,
		ifelse(SEQ=="a",2,
		ifelse(SEQ=="g",3,
		ifelse(SEQ=="c",4,0
		)))))))
	
	#jpeg(paste(Fx,"_alignedNF.jpeg",sep=""),
	#	height=480,width= 1440,quality=200)
	#heatmap.2(SEQi,Colv=NA,scale="none",
	#	labCol='',labRow=Nx,key=F,
	#	trace="none",
	#	col=c("yellow","red","green2",
	#	"violet","blue"),
	#	RowSideColors= colc,margin=c(2,8),
	#	main=paste(x, Ax,sep=":"),cex.main=0.8)
	
	#dev.off()
	
	#Summary per position	
	sump<-t(sapply(c(1:length(SEQ[,1])),function(x){
		seqx<-SEQ[x,]
		return(sapply(c("-","n","a","t","g","c"),function(x){
			return(length(subset(seqx, seqx==x)))
		}))	
	}))
	
	#Number of positions
	NP=length(sump[,1])
	#Number of gapped position
	NG<-length(subset(sump[,1], sump[,1]!=0))
	#Number of Ns position
	NS<-length(subset(sump[,1], sump[,2]!=0))	
	#Number of Ns + gap position
	NSG<-length(subset(sump[,1], sump[,1]+sump[,2]!=0))	
	#Number and identity of strains with more than 5% of Ns positions
	STg<-subset(Nx,sapply(
		c(1:length(SEQ[1,])),function(x){
			return(length(subset(SEQ[,x],SEQ[,x]=='n')))
		})>=(0.05* NP))
		
	return(c(length(Nx),Fx,NP,NG,NS,NSG,
		length(STg),paste(STg,collapse="__")))	
}))	


colnames(GapNs)<-c("Strains","Name","Size","Gaps","Missing","Gap+Missing","NStWith>5%Missing","StWith>5%Missing")

setwd(pS)
write.table(GapNs,"Summary-Genes-After-GapNs-filtering.txt",row.names=F)


######################	
###calculate average pairwise similarity (complete deletion) per gene 

#Genes to remove (see end of step 9)
#OUT<-c("Gene8266","Gene0219","Gene4052","Gene1216")


setwd(pS)


FcoreNEW<-as.vector(read.table("Summary-per-core-gene.txt",header=T)[,2])
#FcoreNEW<-setdiff(FcoreNEW,OUT)

library(seqinr)

setwd(paste(pO,"/OUT-PS",sep=""))

pdf("PS-dist-per-gene.pdf")

par(mar=c(2,2,1,1),bty="n",mfrow=c(6,6))

PS<-t(sapply(FcoreNEW,function(x){
	Fx= x
	Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
	print(Fx)
	setwd(paste(pO,"/OUT-FASTA",sep=""))
	#Get the aligned nucleotide sequence
	SEQ<-read.fasta(paste(Fx,"_alignedNF.fas",sep=""))
	Nx=names(SEQ)
	nx=length(SEQ)
	lx=length(SEQ[[1]])

	#Pairwise similarity matrix
	PSx<-matrix(unlist(sapply(c(1:(nx-1)),function(x){
		X=x
		seqx<-as.vector(SEQ[[X]])
		#positions to keep (no gap nor missing data)
		p1<-intersect(subset(c(1: lx), seqx!="-"),
			subset(c(1: lx), seqx!="n"))
		
		return(sapply(c((X+1):nx),function(x){
			seqy<-as.vector(SEQ[[x]])
			#positions to keep (no gap nor missing data)
			p12<-intersect(p1,
				intersect(subset(c(1: lx), seqy!="-"),
				subset(c(1: lx), seqy!="n")))
			
			return(c(X,x,length(subset(p12, 
				seqx[p12]==seqy[p12]))/length(p12)))
		}))
	})),ncol=3,byrow=T)
	
	z<-c(nx,lx,median(PSx[,3]),
		mean(PSx[,3]),sd(PSx[,3]),
		max(PSx[,3]),min(PSx[,3]))
	
	PSx<-cbind(Nx[PSx[,1]],Nx[PSx[,2]],PSx[,3])
	colnames(PSx)<-c("St1","St2","PS")
	
	setwd(paste(pO,"/OUT-PS",sep=""))
	
	write.table(PSx,
		paste("PairwiseSimilarity_",Fx,".txt",seq=""),
		row.names=F)
	
	hist(as.numeric(PSx[,3]),breaks=100,main=Fx,las=1,cex.axis=0.6,xlab="",ylab="")
	
	return(z)
}))

dev.off()

setwd(pS)
colnames(PS)<-c("sequences","size","med.PS","av.PS","sd.PS","max.PS","min.PS")

write.table(PS,"PS-core-genes.txt")	


###############
##Step 10
##For each possible PS between two strains, retrieve all possible PS values calculated from individual genes and calculate an average

library(ade4)
library(gplots)

#work directory for analyzes
pS="/Users/jean-baptisteleducq/Dropbox/PostDoc-Moscow-ID/Project/myRAST-gene_annotation"
#work directory for annotation files (GFF3 myRAST outputs)
pO="/Users/jean-baptisteleducq/Desktop/myRAST-out"

setwd(pS)

FcoreNEW<-as.vector(read.table("Summary-per-core-gene.txt",header=T)[,2])

#Genes to remove (see end of step 9)
OUT<-c("")
FcoreNEW<-setdiff(FcoreNEW,OUT)

library(seqinr)



NAMES<-read.table("Annotation-Assembly-summary.txt",header=T)

Nx<-as.vector(NAMES[,2])
nx<-length(Nx)


PSo<-cbind('na','na',NA)
colnames(PSo)<-c("St1","St2","PS")

PSall<-sapply(FcoreNEW,function(x){
	print(x)
	Fx= x
	Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
	setwd(paste(pO,"/OUT-PS",sep=""))
	PS<-as.matrix(read.table(
		paste("PairwiseSimilarity_",Fx,".txt",seq=""),
		header=T))
	PS<-unlist(sapply(c(1:(nx-1)),function(x){
		X=x
		sx= Nx[X]
		PSx<-rbind(subset(PS,PS[,1]==sx),
			subset(PS,PS[,2]==sx),PSo)
		sapply(c((X+1):nx),function(x){
			sy=Nx[x]
			PSxx<-rbind(subset(PSx,PSx[,1]== sy),
				subset(PSx,PSx[,2]== sy),PSo)
			return(as.numeric(PSxx[1,3]))
		})
	}))
	return(PS)
})

#Return median value accross genes for each PS

PSmed<-sapply(c(1:length(PSall[,1])),function(x){
	return(median(PSall[x,],na.rm=T))
})

#Format in 2D matrix

COMP<-matrix(unlist(sapply(c(1:(nx-1)),function(x){
	X=x
	sapply(c((X+1):nx),function(x){
		return(c(X,x))
	})
})),ncol=2,byrow=T)

PSmed<-as.matrix(cbind(COMP, PSmed))


PSmed2<-sapply(c(1:nx),function(x){
	X=x
	PSx<-rbind(subset(PSmed[,c(2,3)], PSmed[,1]==X),
		subset(PSmed[,c(1,3)], PSmed[,2]==X),
		cbind(X,1))
	return(PSx[order(PSx[,1]),2])
})

write.table(PSmed2,"Median-Pairwise-Nuc-Matrix.txt",col.names=Nx)

COLcla<-ifelse(NAMES[,6]=="A1","blue",
	ifelse(NAMES[,6]=="A2","cyan",
	ifelse(NAMES[,6]=="A3","blue3",
	ifelse(NAMES[,6]=="A4","brown",
	ifelse(NAMES[,6]=="A5","yellow2",
	ifelse(NAMES[,6]=="A6","orange",
	ifelse(NAMES[,6]=="A7","pink",
	ifelse(NAMES[,6]=="A8","purple",
	ifelse(NAMES[,6]=="A9","red",
	ifelse(NAMES[,6]=="A10","green2",
	ifelse(NAMES[,6]=="B","black",
	ifelse(NAMES[,6]=="C","grey",
	ifelse(NAMES[,6]=="Microvirga","white",
	ifelse(NAMES[,6]=="Enterovirga",rgb(0.9,0.9,0.9),
	"white"))))))))))))))

setwd(pS)
pdf("MedianPS-core-genes.pdf",width=12,height=12)

zall<-heatmap.2(PSmed2,scale="none",trace="none",
	ColSideColors= COLcla,
	RowSideColors= COLcla,labCol=Nx,breaks=200,
	cexCol=0.4,cexRow=0.4,
	labRow=Nx,key.xlab="Median PS",key.ylab="",
	key.title = "",
	key.par=list(las=1,cex.axis=0.8,cex.lab=1.2))

dev.off()

###Check correlation between PSmed and PS for each gene taken individually (indicator of potential HGT)

COR<-as.vector(sapply(c(1:length(FcoreNEW)),function(x){
	PSx<-PSall[,x]
	return(cor.test(PSmed[,3], PSx)$estimate)
}))

#Redo the heatmap only xwwith genes having at least 95% of correlation

i95<-subset(c(1:length(COR)),COR>=0.95)

PSmed95<-sapply(c(1:length(PSall[,1])),function(x){
	return(median(PSall[x, i95],na.rm=T))
})
PSmed95 <-as.matrix(cbind(COMP, PSmed95))

PSmed95 <-sapply(c(1:nx),function(x){
	X=x
	PSx<-rbind(subset(PSmed95[,c(2,3)], PSmed95[,1]==X),
		subset(PSmed95[,c(1,3)], PSmed95[,2]==X),
		cbind(X,1))
	return(PSx[order(PSx[,1]),2])
})

setwd(pS)
pdf("MedianPS-core-genes-COR>0.95.pdf",width=12,height=12)

heatmap.2(PSmed95,scale="none",trace="none",
	ColSideColors= COLcla,
	RowSideColors= COLcla,labCol=Nx,breaks=200,
	cexCol=0.4,cexRow=0.4,
	labRow=Nx,key.xlab="Median PS",key.ylab="",
	key.title = "",
	key.par=list(las=1,cex.axis=0.8,cex.lab=1.2))

dev.off()


#Display PCA
ACP=dudi.pca(PSmed2, scannf=T)
3
part<-round(100* ACP $eig/sum(ACP $eig),2)


x1=min(ACP $li[,1])
x2=max(ACP $li[,1])

y1=min(ACP $li[,2])
y2=max(ACP $li[,2])

z1=min(ACP $li[,3])
z2=max(ACP $li[,3])


pdf("PCA-PS-core-genes.pdf",height=8,width=4)
par(mar=c(4,4,1,1),bty="n",mfrow=c(2,1))
plot(0,0,col=NA,ylim=c(y1,y2),xlim=c(x1,x2),
		las=1,cex.lab=1,cex.axis=0.8,font.lab=2,
		xlab=paste("Axis 1 (",part[1],"%)",sep=""),
		ylab=paste("Axis 2 (",part[2],"%)",sep=""),
		main="",cex.main=1)
		
points(ACP $li[,c(1,2)],pch=21,bg= COLcla,col="black")		
plot(0,0,col=NA,xlim=c(y1,y2),ylim=c(z1,z2),
		las=1,cex.lab=1,cex.axis=0.8,font.lab=2,
		xlab=paste("Axis 2 (",part[2],"%)",sep=""),
		ylab=paste("Axis 3 (",part[3],"%)",sep=""),
		main="",cex.main=1)
		
points(ACP $li[,c(2,3)],pch=21,bg= COLcla,col="black")		
dev.off() 

#####################
###Step 11: produce a fasta file with concatenated aligned sequences for Core gene present in all genomes for partitionning analysis

library(seqinr)
library(gplots)
library(ade4)

#work directory for analyzes
pS="/Users/jean-baptisteleducq/Dropbox/PostDoc-Moscow-ID/Project/myRAST-gene_annotation"
#work directory for annotation files (GFF3 myRAST outputs)
pO="/Users/jean-baptisteleducq/Desktop/myRAST-out"

setwd(pS)

FcoreNEW<-read.table("Summary-per-core-gene.txt",header=T)
OUT<-c("")

NAMES<-read.table("Annotation-Assembly-summary.txt",header=T)

Nx<-as.vector(NAMES[,2])
nx<-length(Nx)

#Recalculate occurence
FcoreNEW<-setdiff(FcoreNEW[,2],OUT)

OCC<-as.vector(sapply(FcoreNEW,function(x){
	Fx= x
	Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
	setwd(paste(pO,"/OUT-FASTA",sep=""))
	return(length(read.fasta(paste(Fx,"_alignedNF.fas",sep=""))))
}))			

Fcomp<-subset(FcoreNEW,OCC==nx)

SEQ<-sapply(Nx,function(x){
	print(x)
	X=subset(c(1:length(Nx)),Nx==x)
	SEQs<-as.vector(unlist(sapply(Fcomp,function(x){
		Fx= x
		Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
		setwd(paste(pO,"/OUT-FASTA",sep=""))
		#Get the aligned nucleotide sequence
		SEQx<-read.fasta(paste(Fx,"_alignedNF.fas",sep=""))
		return(SEQx[[X]])
	})))	
	return(list(SEQs))
})

setwd(paste(pO,"/OUT-FASTA",sep=""))

write.fasta(SEQ,names=Nx,paste("Concatenated-core-genes_n=",length(Fcomp),".fas",sep=""))

#Generate a partition file for analysis in RAxML

#DNA, p1=1-30
#DNA, p2=31-60
	
PART<-matrix(unlist(sapply(c(1:length(Fcomp)),function(x){
	Fx= Fcomp [x]
	Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
	setwd(paste(pO,"/OUT-FASTA",sep=""))
	#Get the aligned nucleotide sequence
	SEQx<-read.fasta(paste(Fx,"_alignedNF.fas",sep=""))
	return(t(cbind(x,seq(1,1,length.out=length(SEQx[[1]])))))
})),ncol=2,byrow=T)

PART2<-cumsum(PART[,2])
PART<-sapply(unique(PART[,1]),function(x){
	return(paste("DNA,",Fcomp[x],"=",
	min(subset(PART2,PART[,1]==x)),"-",
	max(subset(PART2,PART[,1]==x)),sep=""))	
})

write.table(PART,paste("Partition-concatenated-core-genes_n=",length(Fcomp),".txt",sep=""),col.names=F,row.names=F,quote=F)

####Produce a concatenated alignment only with the 64 or 128 less or more polymorphic CORE genes present in all genomes

##Calculate average nucleotide diversity per gene (complete deletion)
setwd(paste(pO,"/OUT-FASTA",sep=""))

DIV<-(1-sapply(Fcomp,function(x){
	print(x)
	###aligned fasta file for this gene
	namex<-paste("Gene",as.numeric(unlist(
		strsplit(x,split="Gene")))[2],sep="")
	
	SEQ<-read.fasta(paste(namex,
		"_alignedNF.fas",sep=""))
	#Use the first sequence as reference
	SEQr<-as.vector(SEQ[[1]])
	#Identity with other sequence
	return(mean(apply(
	sapply(c(2:length(SEQ)),function(x){
		SEQx<-as.vector(SEQ[[x]])
		return(ifelse(SEQr=="-",NA,
			ifelse(SEQx=="-",NA,
			ifelse(SEQx==SEQr,1,0))))
	}),1,sum,na.rm=T)/(length(SEQ)-1),na.rm=T))
}))


sapply(c(64,128),function(x){

X=x

#the X Genes with the highest polymoprhism
Fdiv<-Fcomp[order(DIV,decreasing=T)][1:X]

SEQ<-sapply(Nx,function(x){
	print(x)
	X=subset(c(1:length(Nx)),Nx==x)
	SEQs<-as.vector(unlist(sapply(Fdiv,function(x){
		Fx= x
		Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
		setwd(paste(pO,"/OUT-FASTA",sep=""))
		#Get the aligned nucleotide sequence
		SEQx<-read.fasta(paste(Fx,"_alignedNF.fas",sep=""))
		return(SEQx[[X]])
	})))	
	return(list(SEQs))
})

setwd(paste(pO,"/OUT-FASTA",sep=""))

write.fasta(SEQ,names=Nx,paste("Concatenated-core-genes_n=",length(Fdiv),"-MostPolymorph.fas",sep=""))

#Generate a partition file for analysis in RAxML

#DNA, p1=1-30
#DNA, p2=31-60
	
PART<-matrix(unlist(sapply(c(1:length(Fdiv)),function(x){
	Fx= Fdiv [x]
	Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
	setwd(paste(pO,"/OUT-FASTA",sep=""))
	#Get the aligned nucleotide sequence
	SEQx<-read.fasta(paste(Fx,"_alignedNF.fas",sep=""))
	return(t(cbind(x,seq(1,1,length.out=length(SEQx[[1]])))))
})),ncol=2,byrow=T)

PART2<-cumsum(PART[,2])
PART<-sapply(unique(PART[,1]),function(x){
	return(paste("DNA,",Fdiv[x],"=",
	min(subset(PART2,PART[,1]==x)),"-",
	max(subset(PART2,PART[,1]==x)),sep=""))	
})

write.table(PART,paste("Partition-concatenated-core-genes_n=",length(Fdiv),"-MostPolymorph.txt",sep=""),col.names=F,row.names=F,quote=F)

})

sapply(c(64,128),function(x){

X=x

#the X Genes with the highest polymoprhism
Fdiv<-Fcomp[order(DIV,decreasing=F)][1:X]

SEQ<-sapply(Nx,function(x){
	print(x)
	X=subset(c(1:length(Nx)),Nx==x)
	SEQs<-as.vector(unlist(sapply(Fdiv,function(x){
		Fx= x
		Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
		setwd(paste(pO,"/OUT-FASTA",sep=""))
		#Get the aligned nucleotide sequence
		SEQx<-read.fasta(paste(Fx,"_alignedNF.fas",sep=""))
		return(SEQx[[X]])
	})))	
	return(list(SEQs))
})

setwd(paste(pO,"/OUT-FASTA",sep=""))

write.fasta(SEQ,names=Nx,paste("Concatenated-core-genes_n=",length(Fdiv),"-LessPolymorph.fas",sep=""))

#Generate a partition file for analysis in RAxML

#DNA, p1=1-30
#DNA, p2=31-60
	
PART<-matrix(unlist(sapply(c(1:length(Fdiv)),function(x){
	Fx= Fdiv [x]
	Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
	setwd(paste(pO,"/OUT-FASTA",sep=""))
	#Get the aligned nucleotide sequence
	SEQx<-read.fasta(paste(Fx,"_alignedNF.fas",sep=""))
	return(t(cbind(x,seq(1,1,length.out=length(SEQx[[1]])))))
})),ncol=2,byrow=T)

PART2<-cumsum(PART[,2])
PART<-sapply(unique(PART[,1]),function(x){
	return(paste("DNA,",Fdiv[x],"=",
	min(subset(PART2,PART[,1]==x)),"-",
	max(subset(PART2,PART[,1]==x)),sep=""))	
})

write.table(PART,paste("Partition-concatenated-core-genes_n=",length(Fdiv),"-LessPolymorph.txt",sep=""),col.names=F,row.names=F,quote=F)

})

####Produce a concatenated alignment only with  CORE genes with the most/less conserved synteny: 64 or 128 more or less conserved synteny

##Average nucleotide diversity per gene (complete deletion) - for all core genes
setwd(paste(pO,"/OUT-FASTA",sep=""))

DIVall<-(1-sapply(FcoreNEW,function(x){
	print(x)
	###aligned fasta file for this gene
	namex<-paste("Gene",as.numeric(unlist(
		strsplit(x,split="Gene")))[2],sep="")
	
	SEQ<-read.fasta(paste(namex,
		"_alignedNF.fas",sep=""))
	#Use the first sequence as reference
	SEQr<-as.vector(SEQ[[1]])
	#Identity with other sequence
	return(mean(apply(
	sapply(c(2:length(SEQ)),function(x){
		SEQx<-as.vector(SEQ[[x]])
		return(ifelse(SEQr=="-",NA,
			ifelse(SEQx=="-",NA,
			ifelse(SEQx==SEQr,1,0))))
	}),1,sum,na.rm=T)/(length(SEQ)-1),na.rm=T))
}))

write.table(DIVall,"diversity-per-core-gene.txt")

#synteny directory
pP="/Users/jean-baptisteleducq/Desktop/myRAST-out/OUT-synteny"

#Average synteny per gene
setwd(pP)
SUMscoreG <-read.table("Summary-SI-per-gene.txt",header=T)

pdf("Synteny-vs-polymoph.core-genes.pdf",width=4,height=4)
par(mar=c(4,4,1,1),bty="n")
plot(SUMscoreG[,3],DIVall,las=1,xlab="Synteny per core gene",ylab="Nuc. polymorphism per core gene",cex=0.8,pch=19,cex.axis=0.8,col=rgb(0,0,0,0.3))
dev.off()

#Format names from core gene present in all genomes
Fcomp2<- as.vector(sapply(Fcomp,function(x){
	return(paste("Gene",as.numeric(unlist(
		strsplit(x,split="Gene")))[2],sep=""))
}))

#Retrieve their synteny
SI<-sapply(Fcomp2,function(x){
	return(subset(SUMscoreG[,3],SUMscoreG[,1]==x))
})

sapply(c(64,128),function(x){

X=x

#the X Genes with the highest synteny
Fsi<-Fcomp[order(SI,decreasing=T)][1:X]

SEQ<-sapply(Nx,function(x){
	print(x)
	X=subset(c(1:length(Nx)),Nx==x)
	SEQs<-as.vector(unlist(sapply(Fsi,function(x){
		Fx= x
		Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
		setwd(paste(pO,"/OUT-FASTA",sep=""))
		#Get the aligned nucleotide sequence
		SEQx<-read.fasta(paste(Fx,"_alignedNF.fas",sep=""))
		return(SEQx[[X]])
	})))	
	return(list(SEQs))
})

setwd(paste(pO,"/OUT-FASTA",sep=""))

write.fasta(SEQ,names=Nx,paste("Concatenated-core-genes_n=",length(Fsi),"-HighestSynteny.fas",sep=""))

#Generate a partition file for analysis in RAxML

#DNA, p1=1-30
#DNA, p2=31-60
	
PART<-matrix(unlist(sapply(c(1:length(Fsi)),function(x){
	Fx= Fsi [x]
	Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
	setwd(paste(pO,"/OUT-FASTA",sep=""))
	#Get the aligned nucleotide sequence
	SEQx<-read.fasta(paste(Fx,"_alignedNF.fas",sep=""))
	return(t(cbind(x,seq(1,1,length.out=length(SEQx[[1]])))))
})),ncol=2,byrow=T)

PART2<-cumsum(PART[,2])
PART<-sapply(unique(PART[,1]),function(x){
	return(paste("DNA,",Fsi[x],"=",
	min(subset(PART2,PART[,1]==x)),"-",
	max(subset(PART2,PART[,1]==x)),sep=""))	
})

write.table(PART,paste("Partition-concatenated-core-genes_n=",length(Fsi),"-HighestSynteny.txt",sep=""),col.names=F,row.names=F,quote=F)


#the X Genes with the lowest synteny
Fsi<-Fcomp[order(SI,decreasing=F)][1:X]

SEQ<-sapply(Nx,function(x){
	print(x)
	X=subset(c(1:length(Nx)),Nx==x)
	SEQs<-as.vector(unlist(sapply(Fsi,function(x){
		Fx= x
		Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
		setwd(paste(pO,"/OUT-FASTA",sep=""))
		#Get the aligned nucleotide sequence
		SEQx<-read.fasta(paste(Fx,"_alignedNF.fas",sep=""))
		return(SEQx[[X]])
	})))	
	return(list(SEQs))
})

setwd(paste(pO,"/OUT-FASTA",sep=""))

write.fasta(SEQ,names=Nx,paste("Concatenated-core-genes_n=",length(Fsi),"-LowestSynteny.fas",sep=""))

#Generate a partition file for analysis in RAxML

#DNA, p1=1-30
#DNA, p2=31-60
	
PART<-matrix(unlist(sapply(c(1:length(Fsi)),function(x){
	Fx= Fsi [x]
	Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
	setwd(paste(pO,"/OUT-FASTA",sep=""))
	#Get the aligned nucleotide sequence
	SEQx<-read.fasta(paste(Fx,"_alignedNF.fas",sep=""))
	return(t(cbind(x,seq(1,1,length.out=length(SEQx[[1]])))))
})),ncol=2,byrow=T)

PART2<-cumsum(PART[,2])
PART<-sapply(unique(PART[,1]),function(x){
	return(paste("DNA,",Fsi[x],"=",
	min(subset(PART2,PART[,1]==x)),"-",
	max(subset(PART2,PART[,1]==x)),sep=""))	
})

write.table(PART,paste("Partition-concatenated-core-genes_n=",length(Fsi),"-LowestSynteny.txt",sep=""),col.names=F,row.names=F,quote=F)

})

})

#####################
###Step 12: produce a fasta file with concatenated aligned sequences for all Core gene for partitionning analysis - replace missing sequences by missing data "---"

setwd(pS)

FcoreNEW<-read.table("Summary-per-core-gene.txt",header=T)[,2]


OCC2<-sapply(FcoreNEW,function(x){
	Fx= x
	Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
	setwd(paste(pO,"/OUT-FASTA",sep=""))
	seqx<-read.fasta(paste(Fx,"_alignedNF.fas",sep=""))
	lx=length(seqx[[1]])
	
	return(lx *sapply(Nx,function(x){
		return(length(intersect(names(seqx),x)))
	}))
})		

#Missing data (formated for RAXML)
SEQo<-sapply(as.vector(apply(OCC2,2,max)),
	function(x){
		X=x
		return(list(sapply(c(1:X),function(x){
			return("n")
		})))
})


SEQ<-sapply(Nx,function(x){
	print(x)
	X=x
	SEQs<-sapply(FcoreNEW,function(x){
		Fx= x
		Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
		setwd(paste(pO,"/OUT-FASTA",sep=""))
		#Get the aligned nucleotide sequence
		SEQx<-read.fasta(paste(Fx,
			"_alignedNF.fas",sep=""))
		return(as.vector(unlist(
			subset(SEQx,names(SEQx)==X))))
	})
	
	sox<-as.numeric(summary(SEQs)[,1])
	SEQox<-subset(SEQo ,sox==0)
	Iox<-subset(c(1:length(FcoreNEW)),sox==0)
	
	SEQs<-subset(SEQs, sox>0)
	Isx<-subset(c(1:length(FcoreNEW)),sox>0)
	
	SEQs<-c(SEQs, SEQox)[order(c(Isx, Iox))]
	
	return(list(as.vector(unlist(SEQs))))
})

setwd(paste(pO,"/OUT-FASTA",sep=""))

write.fasta(SEQ,names=Nx,paste("Concatenated-core-genes_n=",length(FcoreNEW),".fas",sep=""))

#Generate a partition file for analysis in RAxML

#DNA, p1=1-30
#DNA, p2=31-60
	
PART<-matrix(unlist(sapply(c(1:length(FcoreNEW)),function(x){
	Fx= FcoreNEW [x]
	Fx<-paste("Gene",as.numeric(unlist(strsplit(Fx,split="Gene"))[2]),sep="")
	setwd(paste(pO,"/OUT-FASTA",sep=""))
	#Get the aligned nucleotide sequence
	SEQx<-read.fasta(paste(Fx,"_alignedNF.fas",sep=""))
	return(t(cbind(x,seq(1,1,length.out=length(SEQx[[1]])))))
})),ncol=2,byrow=T)

PART2<-cumsum(PART[,2])
PART<-sapply(unique(PART[,1]),function(x){
	return(paste("DNA,", FcoreNEW[x],"=",
	min(subset(PART2,PART[,1]==x)),"-",
	max(subset(PART2,PART[,1]==x)),sep=""))	
})

write.table(PART,paste("Partition-concatenated-core-genes_n=",length(FcoreNEW),".txt",sep=""),col.names=F,row.names=F,quote=F)


####Produce a concatenated alignment only with  CORE genes with the most conserved synteny

#synteny directory
pP="/Users/jean-baptisteleducq/Desktop/myRAST-out/OUT-synteny"

#Average synteny per gene
setwd(pP)
SUMscoreG <-read.table("Summary-SI-per-gene.txt",header=T)



