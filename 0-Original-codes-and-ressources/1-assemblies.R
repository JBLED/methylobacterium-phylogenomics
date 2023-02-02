library(seqinr)
library(gplots)

 p0<-"/Users/jean-baptisteleducq/Desktop/Divers-related-to-PosdocS-Methylo/Megahit_assembler"

setwd(p0)

#Genome names
STRAIN<-c(
"J-092","J-090","J-088","J-078",
"J-077","J-076","J-072","J-070",
"J-068","J-067","J-059","J-048",
"J-043","J-030","J-026","J-001",
"E-066","E-065","E-046","E-045",
"E-041","E-025","E-016","E-005")

pdf("summary-megahit.pdf",width=9,height=3)

SUM<-t(sapply(STRAIN,function(x){
	
	print(x)	
	SEQ<-read.fasta(paste(x,"_final.contigs.fa",sep=""))

	ATT<-t(sapply(c(1:length(SEQ)),function(x){
		attx<-unlist(strsplit(
			attributes(SEQ[[x]])$Annot,split=" "))[2:4]
		return(c(as.numeric(matrix(unlist(
			strsplit(attx,split="=")),
			ncol=2,byrow=T)[,2]),GC(SEQ[[x]])))
	}))

	ATT<-cbind(names(SEQ),ATT)

	colnames(ATT)<-c("Contig","flag","cov","len","GC")

	N1x=length(ATT[,1])

	S1x=sum(as.numeric(ATT[,4]))

	C1x=sum(as.numeric(ATT[,3])*as.numeric(ATT[,4]))/S1x

	SEQ<-subset(SEQ,as.numeric(ATT[,3])>10)
	ATT<-subset(ATT,as.numeric(ATT[,3])>10)

	write.table(ATT,paste(x,"_summary.txt",sep=""))

	write.fasta(SEQ,names=names(SEQ),paste(x,"_filtered.fa",sep=""))

	#plot(as.numeric(ATT[,4]),as.numeric(ATT[,3]),log="",xlab="length",ylab="coverage")

	N2x=length(ATT[,1])

	S2x=sum(as.numeric(ATT[,4]))

	C2x=sum(as.numeric(ATT[,3])*as.numeric(ATT[,4]))/S1x

	Cmdx=median(as.numeric(ATT[,3]))

	M2x=max(as.numeric(ATT[,4]))

	par(mar=c(4,4,1,1),bty="n",mfrow=c(1,3))

	COVm=max(as.numeric(ATT[,3]))

	plot(as.numeric(ATT[,4])/1000,
		as.numeric(ATT[,3]),log="x",
		xlab="length (Kb, log scale)",
		ylab="coverage",cex.axis=0.8,
		las=1,main=x,pch=19,cex=0.6,col=rgb(0,0,0,0.5))

	segments(0.001, Cmdx,M2x, Cmdx,lwd=2,col="red")

	text(0.0001*M2x,0.9* COVm,
		paste("cov. = ",round(Cmdx,0),
		" ; max = ",round(M2x/1000)," Kb",sep=""),
		cex=0.7,font=2)

	text(0.0001*M2x,0.8* COVm,
		paste(N2x," contigs ; total size = ",
		round(S2x/1000000,3)," Mb",sep=""),
		cex=0.7,font=2)

	par(mar=c(4,4,1,1),bty="n")

	plot(cumsum(as.numeric(
		ATT[order(as.numeric(ATT[,4]),
		decreasing=T),4]))/1000000,type="l",
		xlab="contigs (log scale)",
		ylab="cumulated size (Kb)",las=1,cex.axis=0.8,log="x",
		col="black",lwd=2)

	Splot<-as.numeric(ATT[,4])/max(as.numeric(ATT[,4]))

	Splot =4*sqrt(Splot/pi)

	par(mar=c(4,4,1,1),bty="n")
	plot(as.numeric(ATT[,3]),
		as.numeric(ATT[,5]),pch=19,
		col=rgb(0,0,0,0.2),
		cex.axis=0.8,las=1,xlab="coverage",
		ylab="GC content",cex= Splot)


	GCd<-unlist(sapply(c(1:length(ATT[,1])),
		function(x){
		lx=as.numeric(ATT[x,4])
		gcx=as.numeric(ATT[x,5])
		return(seq(gcx,gcx,length.out=lx))	
	}))

	return(c(x,N1x,S1x,C1x,N2x,S2x,
		C2x, Cmdx ,M2x,mean(GCd),
		median(GCd),sd(GCd)))
	
}))

dev.off()


colnames(SUM)<-c("Strain","Contigs","Size","Mean.cov",
	"Contigs(cov>10)","Size(cov>10)","Mean.cov(cov>10)",
	"Med.cov(cov>10)","Max.size","Mean.GC","Median.GC","sd.GC")

	write.table(SUM,"summary-megahit.txt",row.names=F)

	SUM<-read.table("summary-megahit.txt",header=T)

##For contaminated genomes (J-059, E-025 and E-016), generate additionnal filtered fasta file based on coverage

pdf("Cleaning-contamined-genomes.pdf",height=5,width=5)

	STRAINc<-c("J-059", "E-025", "E-016")
	COVc<-c(100,100,150)
	setwd(p0)
SUMclean<-t(sapply(STRAINc,function(x){
	
	print(x)	
	SEQ<-read.fasta(paste(x,"_filtered.fa",sep=""))
	ATT<-read.table(paste(x,"_summary.txt",sep=""),header=T)
	Splot<-as.numeric(ATT[,4])/max(as.numeric(ATT[,4]))

	Splot =4*sqrt(Splot/pi)
	COVx<-subset(COVc,STRAINc==x)
	par(mar=c(4,4,1,1),bty="n")
	plot(ATT[,4]/1000,ATT[,3],pch=19,cex= Splot,col= ifelse(ATT[,3]<=COVx,rgb(1,0,0,0.3),rgb(0,0,0,0.3)),log="x",cex.axis=0.8,las=1,xlab="length (Kb, log scale)",ylab="coverage",main=x)
	
	SEQf<-subset(SEQ,ATT[,3]>COVx)
	namesf<-subset(names(SEQ),ATT[,3]>COVx)
	ATTf<-subset(ATT,ATT[,3]>COVx)
	
	SEQc<-subset(SEQ,ATT[,3]<=COVx)
	namesc<-subset(names(SEQ),ATT[,3]<=COVx)
	ATTc<-subset(ATT,ATT[,3]<=COVx)
	
	GCf<-unlist(sapply(c(1:length(ATTf[,1])),
		function(x){
		lx=as.numeric(ATTf[x,4])
		gcx=as.numeric(ATTf[x,5])
		return(seq(gcx,gcx,length.out=lx))	
	}))
	
	GCf<-c(mean(GCf),sd(GCf))
	
	GCc<-unlist(sapply(c(1:length(ATTc[,1])),
		function(x){
		lx=as.numeric(ATTc[x,4])
		gcx=as.numeric(ATTc[x,5])
		return(seq(gcx,gcx,length.out=lx))	
	}))
	
	GCc<-c(mean(GCc),sd(GCc))	
	
	write.fasta(SEQf,names= namesf,
		paste(x,"_filtered-clean.fa",sep=""))
	write.fasta(SEQc,names= namesc,
		paste(x,"_filtered-contaminant.fa",sep=""))

	Nfx=length(ATTf[,1])
	Sfx=sum(as.numeric(ATTf[,4]))
	Cfx=sum(as.numeric(ATTf[,3])*as.numeric(ATTf[,4]))/Sfx
	Cmfx=median(as.numeric(ATTf[,3]))

	Ncx=length(ATTc[,1])
	Scx=sum(as.numeric(ATTc[,4]))
	Ccx=sum(as.numeric(ATTc[,3])*as.numeric(ATTc[,4]))/Scx
	Cmcx=median(as.numeric(ATTc[,3]))
		
	return(c(Nfx,Sfx,Cfx,GCf,Ncx,Scx,Ccx,GCc))
		
		
}))

dev.off()

############## CHECK FOR ERROR FILES WHEN GENOMES WERE SUBMITTED TO NCBI (FEBRUARY 2022): CONTAMINATIONS WITH ADAPTORS IN SCAFFOLDS

###DATA FILE WITH LIST OF ERRORS

 p0<-"/Users/jean-baptisteleducq/Desktop/Divers-related-to-PosdocS-Methylo/Megahit_assembler/Error-Submission-NCBI"

setwd(p0)

ERRORS<-read.table("data.txt",header=T)

###Source folder for MEGAHIT assemblies before and after filtering

pS<-"/Users/jean-baptisteleducq/Desktop/Divers-related-to-PosdocS-Methylo/Megahit_assembler/" 

###Source folder for Annotation files

pA<-"/Users/jean-baptisteleducq/Desktop/myRAST-out/"

library(seqinr)

#genome names
ST<-as.vector(ERRORS[,1])


SUMall<-sapply(ST,function(x){
	STx=x
	
	setwd(pS)
	#original assembly (MEGAHIT OUTPUT)
	Nm<-subset(ERRORS[,2],ERRORS[,1]==STx)
	SEQm<-read.fasta(Nm)
	#filtered assembly (FROM WHICH ANALYSES WERE CONDUCTED)
	NF<-subset(ERRORS[,3],ERRORS[,1]==STx)
	SEQf<-read.fasta(NF)	
	
	#Annotation file
	setwd(pA)
	AN<-subset(ERRORS[,4],ERRORS[,1]==STx)
	ANNO<-read.delim(AN,header=T)
	
	setwd(p0)
	#ERROR FILE: RemainingContamination
	ER<-subset(ERRORS[,5],ERRORS[,1]==STx)
	CONTer<-read.delim(ER)	
	#ERROR FILE: FixedForeignContaminations
	EF<-subset(ERRORS[,6],ERRORS[,1]==STx)
	CONTef<-read.delim(EF)
	#ERROR FILE: Contamination
	EC<-subset(ERRORS[,7],ERRORS[,1]==STx)
	CONTec<-read.delim(EC)	
	
	#Check for the occurence of scaffold from the clean (analysed) genome file in NCBI error files
	
	Check<-t(sapply(c(1:length(SEQf)),
		function(x){
		nx=names(SEQf)[x]
		ERx<-length(intersect(CONTer[,1],nx))
		EFx<-length(intersect(CONTef[,1],nx))
		ECx<-length(intersect(CONTec[,1],nx))
		return(c(x, ERx, EFx, ECx))	
	}))
	rownames(Check)<-names(SEQf)
	Check<-subset(Check,
		apply(Check[,c(2:4)],1,max)!=0)
		
	#Scaffold characteristics (size and occurence in annotation file)
	
	SIZE<-as.numeric(summary(SEQf)[Check[,1],1])
	
	ANNOT<-sapply(rownames(Check),function(x){
		return(length(subset(ANNO[,1],
			ANNO[,1]==x)))
	})
	
	SPAN<-sapply(rownames(Check),function(x){
		return(unique(subset(CONTec[,3], CONTec[,1]==x)))
		
	})
		
	SUM<-cbind(Check, SIZE, ANNOT, SPAN)
	
	colnames(SUM)<-c("Line","Remaining","Fixed","All","ScafSize","Annotations","Span/Foreign")
	print(STx)
	print(SUM)
	return(list(SUM))
	
})

##PolyA and T check

POLYa<-sapply(ST,function(x){
	STx=x
	print(x)
	setwd(pS)
	#filtered assembly (FROM WHICH ANALYSES WERE CONDUCTED)
	NF<-subset(ERRORS[,3],ERRORS[,1]==STx)
	SEQf<-read.fasta(NF)	
	nx=names(SEQf)
	SEQf<-sapply(c(1:length(SEQf)),function(x){
		return(paste(SEQf[[x]],collapse=""))
	})
	
	
	return(nx [c(grep("aaaaaaaaa", SEQf),
		grep("ttttttttt", SEQf))])
})

pF<-"/Users/jean-baptisteleducq/Dropbox/PostDoc-Moscow-ID/Project/Genome-Submission-NCBI/Manually-Fixed-Assemblies-For-NCBI"

setwd(pF)

FIX<-read.table("Report-manual-fixation.txt",header=T)

SS=50

SEQer<-sapply(c(1:length(FIX[,1])),function(x){
	print(x)
	setwd(pS)
	SEQx<-read.fasta(paste(FIX[x,1],".fa",sep=""))
	SEQx<-as.vector(SEQx[[subset(
		c(1:length(SEQx)),
		names(SEQx)==FIX[x,2])]])
	cx<-as.numeric(unlist(strsplit(FIX[x,4],
		split="[..]")))
	cx<-subset(cx,is.na(cx)==F)
	sx=as.numeric(FIX[x,3])
	
	x1<-cx[1]-SS
	x2<-cx[2]+SS
	
	x1<-ifelse(is.na(x1)==T,1,
		ifelse(x1<1,1,x1))
	x2<-ifelse(is.na(x2)==T, sx,
		ifelse(x2>sx,sx,x2))
	
	return(SEQx[x1:x2])
	
})

setwd(pF)

write.fasta(SEQer,names=paste(FIX[,1],"-",FIX[,2],sep=""),"Errors.fas")

ix<-subset(c(1:length(FIX[,1])),
	is.na(FIX[,4])==F)

SEQfix<-sapply(ix,function(x){
	print(x)
	setwd(pS)
	SEQx<-read.fasta(paste(FIX[x,1],".fa",sep=""))
	SEQx<-as.vector(SEQx[[subset(
		c(1:length(SEQx)),
		names(SEQx)==FIX[x,2])]])
	cx<-as.numeric(unlist(strsplit(FIX[x,4],
		split="[..]")))
	cx<-subset(cx,is.na(cx)==F)
	sx=as.numeric(FIX[x,3])
	
	Ns<-paste(SEQx[c(cx[1]:cx[2])],collapse="")
	Rs<-paste(ifelse(seq(1,1,
		length.out=nchar(Ns))==1,
		"n",""),collapse="")
	SEQx<-paste(SEQx,collapse="")
	SEQx<-gsub(Ns,Rs, SEQx)	

	x1<-cx[1]-SS
	x2<-cx[2]+SS
	
	x1<-ifelse(is.na(x1)==T,1,
		ifelse(x1<1,1,x1))
	x2<-ifelse(is.na(x2)==T, sx,
		ifelse(x2>sx,sx,x2))
		
	return(as.vector(unlist(
		strsplit(SEQx,split="")))[x1:x2])		
})

setwd(pF)

write.fasta(SEQfix,names=paste(FIX[ix,1],"-",FIX[ix,2],sep=""),"Errors-fixed.fas")

###NOW FIXING THESE DAMN GENOMES TO PLEASE NCBI


sapply(unique(FIX[,1]),function(x){
	setwd(pS)
	print(x)
	SEQx<-read.fasta(paste(x,".fa",sep=""))
	
	scax<-unique(subset(FIX[,2],FIX[,1]==x))
	spanx<-sapply(scax,function(x){
		return(unique(subset(FIX[,4],
			FIX[,2]==x)))
	})
	
	scaOK<-setdiff(names(SEQx), scax)
	
	#remove scaffolds marked as contaminated
	scax<-subset(scax,is.na(spanx)==F)
	spanx <-subset(spanx,is.na(spanx)==F)
	
	SEQfix<-sapply(scax,function(x){
		seqx<-SEQx[[subset(c(1:length(SEQx)),
			names(SEQx)==x)]]
		cx<-subset(spanx,scax==x)
		cx<-as.numeric(unlist(strsplit(cx,
			split="[..]")))
		cx<-subset(cx,is.na(cx)==F)
		sx=length(seqx)
	
		Ns<-paste(seqx[c(cx[1]:cx[2])],
			collapse="")
		Rs<-paste(ifelse(seq(1,1,
			length.out=nchar(Ns))==1,
			"n",""),collapse="")
		seqx <-paste(seqx,collapse="")
		seqx <-gsub(Ns,Rs, seqx)		
		return(list(as.vector(unlist(
			strsplit(seqx,split="")))))
	})
	
	SEQok<-sapply(scaOK,function(x){
		seqx<-SEQx[[subset(c(1:length(SEQx)),
			names(SEQx)==x)]]
		return(list(seqx))	
	})
	
	setwd(pF)
	
	write.fasta(c(SEQok, SEQfix),
		names=c(scaOK, scax),
		paste(x,".fa",sep=""))
		
	
})

####Do manually: remove small scaffolds with strenght of N and trim large scaffolds with ends with Ns

##Then compare genome statistics before and after correction

STAT<-t(sapply(unique(FIX[,1]),function(x){
	print(x)
	setwd(pS)
	SEQx<-read.fasta(paste(x,".fa",sep=""))
	setwd(pF)
	SEQf<-read.fasta(paste(x,".fa",sep=""))
	
	return(c(
		length(SEQx),
		length(SEQf),
		sum(as.numeric(summary(SEQx)[,1])),
		sum(as.numeric(summary(SEQf)[,1]))	
	))
}))	

###### sample K mers in contaminated genomes and look for them in clean genomes

Nk=100
Ns=16

setwd(pS)

CONTx<-read.fasta("J-059_filtered-contaminant.fa")

SIZE<-as.numeric(summary(CONTx)[,1])
SIZEK<-ceiling(Nk*SIZE/sum(SIZE))

Km<-matrix(unlist(sapply(
	c(1:length(CONTx)),function(x){
	#print(x)
	seqx<-as.vector(CONTx[[x]])
	Kx=sort(sample(c(1:(length(seqx)-Ns+1)),
		size= SIZEK[x],replace=F))
	return(t(cbind(x,Kx,t(sapply(Kx,function(x){
		return(c(
		paste(seqx[x:(x+ Ns-1)],collapse=""),
		paste(rev(comp(seqx[x:(x+ Ns-1)])),
		collapse="")
		))
	})))))	
})),ncol=4,byrow=T)


CLEANx<-read.fasta("J-059_filtered-clean.fa")

SIZEc<-as.numeric(summary(CLEANx)[,1])
SIZEKc<-ceiling(Nk* SIZEc/sum(SIZEc))

Kmc<-matrix(unlist(sapply(
	c(1:length(CLEANx)),function(x){
	#print(x)
	seqx<-as.vector(CLEANx[[x]])
	Kx=sort(sample(c(1:(length(seqx)-Ns+1)),
		size= SIZEKc[x],replace=F))
	return(t(cbind(x,Kx,t(sapply(Kx,function(x){
		return(c(
		paste(seqx[x:(x+ Ns-1)],collapse=""),
		paste(rev(comp(seqx[x:(x+ Ns-1)])),
		collapse="")
		))
	})))))	
})),ncol=4,byrow=T)


sapply(ST,function(x){
	STx=x
	#filtered assembly (FROM WHICH ANALYSES WERE CONDUCTED)
	NF<-subset(ERRORS[,3],ERRORS[,1]==STx)
	SEQf<-read.fasta(NF)
	
	SEQf<-sapply(c(1:length(SEQf)),function(x){
		return(paste(SEQf[[x]],collapse=""))
	})
	
	Kscore<-sapply(c(1:length(Km[,1])),
	function(x){
		print(x)
		return(c(
			length(grep(Km[x,3], SEQf)),
			length(grep(Km[x,4], SEQf))
		))
	})
	Kscorec<-sapply(c(1:length(Kmc[,1])),
	function(x){
		print(x)
		return(c(
			length(grep(Kmc[x,3], SEQf)),
			length(grep(Kmc[x,4], SEQf))
		))
	})	
	
})

#######################################################CHECK quickly contaminated genomes

setwd("/Users/jean-baptisteleducq/Dropbox/PostDoc-Moscow-ID/Project/Genome-Submission-NCBI/Contaminant_Genomes")

library(seqinr)

SeqJ59<-read.fasta("J-059_filtered-contaminant.fa")

GCc<-sapply(c(1:length(SeqJ59)),
	function(x){
	return(GC(SeqJ59[[x]]))
})

SIZE<-as.numeric(summary(SeqJ59)[,1])

MAXseq<-subset(c(1:length(SeqJ59)), SIZE==max(SIZE))

MAXseq<-unlist(sapply(c(1:length(SeqJ59)),
	function(x){
	return(unlist(SeqJ59[[x]]))
}))

nSEQ=1000

SEQi<-round(seq(1,length(MAXseq)+1,length.out=(nSEQ+1)))
SEQi<-cbind(SEQi[1: nSEQ],SEQi[2:(nSEQ+1)]-1)

GCmax<-sapply(c(1: nSEQ),function(x){
	x1= SEQi[x,1]
	x2= SEQi[x,2]
	seqx<-MAXseq[x1:x2]
	return(GC(seqx))
})

plot(GCmax)

##############################################################################################################################################################################################################################################################################################OLD
####Tetranucleotide diversity to distinguish between scafffolds from different data

#Size of Kmers
f=4
#Size of sub-regions
s=100000

x="E-016"


SEQ<-read.fasta(paste(x,"_filtered.fa",sep=""))
ATT<-read.table(paste(x,"_summary.txt",sep=""),header=T)
sx=length(SEQ)

dir.create("out")
setwd(paste(p0,"out",sep="/"))
	print("scaffold and subregions coordinates")
	COORD<-matrix(unlist(sapply(c(1:sx),function(x){ #scaffold level
		X=x
		SCAx<-SEQ[[x]]
		SCAx <-subset(SCAx, SCAx!="-")
		M<-length(SCAx)
		N<-ceiling(M/s)
		W<-round(seq(0,M,length.out=N+1))
		Ws<-W[1:N]+1
		We<-W[2:(N+1)]
		return(t(cbind(X,c(1:N),Ws,We)))		
	})),ncol=4,byrow=T)
	
	write.table(COORD,"coordinates",row.names=F,col.names=F)
	
	#Number of subregions per chromosome
	NR<-sapply(unique(COORD[,1]),function(x){
		return(length(subset(COORD[,1],COORD[,1]==x)))	
	})
	
	print("Generate windows, unique and redundant per region")
	sapply(c(1:length(COORD[,1])),function(x){ #chromosome level
		X=COORD[x,1]
		R=COORD[x,2]
		print(paste("Sca",X,"reg",R))
		st<-COORD[x,3]
		en<-COORD[x,4]
		SCAxx<-c(SEQ[[X]], SEQ[[X]][1:f])[st:(en+f-1)]
		Mw<-length(SCAxx)	
		S<-c(1:(Mw-f+1))
		E<-c((f): Mw)
		#progx<-round(seq(1, Mw,length.out=101),0)
		KMER<-sapply(c(1:length(S)),function(x){		
			#ifelse(subset(progx,progx==x)=="","",print(subset(progp,progx==x)))
			return(c(paste(SCAxx[S[x]:E[x]],collapse=""),
			paste(rev(comp(SCAxx[S[x]:E[x]])),collapse="")))
		})	
		KMER <-summary(as.factor(c(KMER[1,], KMER[2,])), 	
			maxsum=length(unique(c(KMER[1,], KMER[2,]))))
		write.table(t(t(KMER)),
			paste("kmers",X,R,sep="_"),col.names=F)
		
	})	
	
	SEQr<-t(0)
	rownames(SEQr)<-as.character(c("n"))
	
	print("merge kmer abundance per scaffold")	
	sapply(unique(COORD[,1]),function(x){ #scaffold level
		X=x
		COORDx<-subset(COORD,COORD[,1]==X)
		write.table(SEQr,paste("kmers",X,sep="_"),
			col.names=F)	
		
		sapply(unique(COORDx[,2]),function(x){ #subregion level
			R=x
			print(x)	
			KMER <-read.table(paste("kmers",X,R,sep="_"))
			z<-read.table(paste("kmers",X,sep="_"))
			
			#All Kmers
			kall<-unique(c(as.vector(KMER[,1]),as.vector(z[,1])))		
			ksum<-t(t(sapply(kall,function(x){
				return(sum(c(subset(KMER[,2], KMER[,1]==x),
					subset(z[,2],z[,1]==x))))
			})))
			write.table(ksum,paste("kmers",X,sep="_"),
			col.names=F)			
			file.remove(paste("kmers",X,R,sep="_"))
		})
	})	
		
	####unique pairs of kmers and their Reverse complement 
	X=subset(c(1:sx),as.numeric(summary(SEQ)[,1])==
		max(as.numeric(summary(SEQ)[,1])))
	z<-read.table(paste("kmers",X,sep="_"))
	z<-subset(z,z[,1]!="n")
	
	zu<-matrix(unlist(strsplit(unique(sapply(as.vector(z[,1]),function(x){
		zrc<-paste(rev(comp(unlist(strsplit(x,split="")))),
			collapse="")
		return(paste(sort(c(x, zrc)),collapse="_"))
	})),split="_")),ncol=2,byrow=T)
	
	###Build a matrix of kmer abundance for each scaffold
	
	MATk<-sapply(c(1:sx),function(x){ #scaffold level
		X=x	
		z<-read.table(paste("kmers",X,sep="_"))
		z<-sapply(zu[,1],function(x){
			return(c(subset(z[,2],z[,1]==x),0)[1])
		})
		return(z/sum(z))
	})
		
	colr<-c(
	seq(0,0,length.out= 26),
	seq(0,1,length.out= 26),
	seq(1,1,length.out= 26),
	seq(1,1,length.out= 26)
	)

	colg<-c(
	seq(0,1,length.out= 26),
	seq(1,1,length.out= 26),
	seq(1,0.5,length.out= 26),
	seq(0.5,0,length.out= 26)
	)
	
	colb<-c(
	seq(1,1,length.out= 26),
	seq(1,0,length.out= 26),
	seq(0,0,length.out= 26),
	seq(0,0,length.out= 26)
	)


Splot<-as.numeric(ATT[,4])/max(as.numeric(ATT[,4]))

Splot =4*sqrt(Splot/pi)



	#Contigs >Nx bp
	Nx=20000
	C1kb<-subset(c(1:sx),ATT[,4]>= Nx)

	COVx<-ATT[C1kb,3]-min(ATT[C1kb,3])
	COVx<-round(100*COVx/max(COVx))+1
	
	COLx<-unique(rgb(colr,colg,colb))[COVx]
	
	plot(ATT[C1kb,4],ATT[C1kb,3],pch=1,cex= Splot[C1kb],col= COLx,log="x")		

	
	heatmap.2(MATk[, C1kb],scale="row",trace="none",ColSideColors= COLx)
	

library(ade4)

#Display PCA
ACPk=dudi.pca(t(MATk[, C1kb]), scannf=T)
2
part<-round(100* ACPk $eig/sum(ACPk $eig),2)


x1=min(ACPk $li[,1])
x2=max(ACPk $li[,1])
y1=min(ACPk $li[,2])
y2=max(ACPk $li[,2])

par(mar=c(4,4,1,1),bty="n")
plot(ACPk $li[,1],ACPk $li[,2],col= COLx,ylim=c(y1,y2),xlim=c(x1,x2),
		las=1,cex.lab=1,cex.axis=0.8,font.lab=2,
		xlab=paste("Axis 1 (",part[1],"%)",sep=""),
		ylab=paste("Axis 2 (",part[2],"%)",sep=""),
		main="",cex.main=1,pch=1,cex= Splot[C1kb])
