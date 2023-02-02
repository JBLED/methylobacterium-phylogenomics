p0<-"/Users/jean-baptisteleducq/Desktop/New-Methylo-Genomes/Alessa-et-al-2021-genomes/"

setwd(paste(p0,"Synteny","Synteny-Maps",sep="/"))

#Nodes in the astral species tree corresponding to main groups
Nodes<-c(20,59,81,99,116,117)
Group<-c("B","A","D","C","M","E")

#Reference species for mapping of links on core genome
REF="Meth_planium"

#Links in this species genome
MapRef<-read.table(paste("Consensus-Map-",REF,".txt",sep=""),header=T)

LINKref<-paste(MapRef[,1], MapRef[,2],sep="_")

#Links in all groups
LINKgroup<-unique(unlist(sapply(Nodes,
	function(x){
	MAPx<-read.table(paste(
		"Consensus-Map-",x,".txt",sep=""),
		header=T)
	return(list(paste(MAPx[,1], 
		MAPx[,2],sep="_")))	
})))

#Occurence of these links in the reference genome
OCCref<-as.vector(sapply(LINKgroup,function(x){
	return(length(intersect(x, LINKref)))
}))

#Occurence of these links in each group

OCCmeth<-sapply(Nodes[1:4],
	function(x){
	MAPx<-read.table(paste(
		"Consensus-Map-",x,".txt",sep=""),
		header=T)
	Gx=subset(Group,Nodes==x)
	LINKx<-paste(MAPx[,1], 
		MAPx[,2],sep="_")
	OCCx<-as.vector(sapply(LINKgroup,
		function(x){
		return(c(subset(MAPx[,3], 
		LINKx==x),0)[1]
		)
	}))
	#return(ifelse(OCCx==0,"",
	#	ifelse(OCCx>=0.5,Gx,tolower(Gx))))
	return(ifelse(OCCx>0.5,Gx,""))		
})

OCCmeth <-apply(OCCmeth,1,paste,collapse="")

OCCout<-sapply(Nodes[5:6],
	function(x){
	MAPx<-read.table(paste(
		"Consensus-Map-",x,".txt",sep=""),
		header=T)
	Gx=subset(Group,Nodes==x)
	LINKx<-paste(MAPx[,1], 
		MAPx[,2],sep="_")
	OCCx<-as.vector(sapply(LINKgroup,
		function(x){
		return(c(subset(MAPx[,3], 
		LINKx==x),0)[1]
		)
	}))
	#return(ifelse(OCCx==0,"",
	#	ifelse(OCCx>=0.5,Gx,tolower(Gx))))
	return(ifelse(OCCx>0.5,Gx,""))		
})

OCCout <-apply(OCCout,1,paste,collapse="")

OCCall<-paste(OCCmeth, OCCout,sep="")

OCCall<-ifelse(OCCall=="BADCME","All",
	ifelse(OCCall=="BADC","Meth",
	ifelse(OCCall=="BADCE","MethEnte",
	ifelse(OCCall=="BADCM","MethMicro",
	ifelse(OCCall=="M","Micro",
	ifelse(OCCall=="E","Ente",
	ifelse(OCCall=="ME","MicroEnte","")))))))

LINKgroup<-matrix(unlist(strsplit(LINKgroup,
	split="_")),ncol=2,byrow=T)
	
SUM<-cbind(LINKgroup, OCCref, OCCall, OCCmeth)	
	
colnames(SUM)<-c("Gene1","Gene2",REF,"Genus","Methylo")

SUM <-subset(SUM,apply(SUM[,3:5],1,
	paste,collapse="")!="0")

SUM=ifelse(SUM=="","o", SUM)
	
write.table(SUM,"Consensus-Synteny.txt",quote=F,
	row.names=F)	
	
