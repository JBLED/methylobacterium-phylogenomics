library(seqinr)
library(gplots)

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
