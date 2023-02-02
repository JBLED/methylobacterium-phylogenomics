

library(seqinr)

#Alessa genomes
p0<-"/Users/jean-baptisteleducq/Desktop/New-Methylo-Genomes/Alessa-et-al-2021-genomes/"

#Aligned core genes
pA<-"/Users/jean-baptisteleducq/Desktop/myRAST-out/OUT-FASTA/"

setwd(p0)

#List of genomes
LIST<-as.vector(read.table("list-genomes.txt")[,1])

#List of core genes
GEN<-read.table("Core-Gene-REF=YIM132548-Cla=A5b.txt",header=T)

##convert gff file in fasta (do M. durans, Mb_radiodurans, Mic_thermotolerans.gff, M terrae manually)

sapply(setdiff(LIST,c("Mb_durans.gff","Mb_radiodurans.gff","Mb_terrae.gff","Mic_thermotolerans.gff","Mr_extorquens.gff","Mr_populi.gff")),function(x){
	print(x)
	setwd(paste(p0,"Alessa-genomes",
		sep="/"))
	genx<-as.vector(read.delim(x)[,1])
	genx<-genx[(grep("##FASTA", genx)+1):length(genx)]
	nx=paste(strsplit(x,
		split="[.]")[[1]][1],
		".fas",sep="")
	write.table(genx,row.names=F,
		col.names=F,quote=F,nx)	
})


#Check core gene abundance in each new genome

setwd(p0)

ANNO<-read.table("Genomes-to-annotate.txt",header=T)


setwd(paste(p0,"Annotations-16-11-2021",
		sep="/"))

COREn<-sapply(c(1:length(ANNO[,1])),
	function(x){
	print(ANNO[x,1])
	ANNOx<-read.delim(ANNO[x,5],	header=T)	
	ID<-ANNOx[,8]	
	return(sapply(c(1:length(GEN[,1])),
		function(x){
		IDx<-GEN[x,3]
		nx<-GEN[x,4]
		return(length(subset(ID,ID== nx)))
	}))
})

#Ribosomal genes missing in M. haplocladii and M. jeotgali, otherwise everything looks OK: only keep core gene in 1:1 copy and report size of the nucleotide sequence

COREs<-sapply(c(1:length(ANNO[,1])),
	function(x){
	print(ANNO[x,1])
	ANNOx<-read.delim(ANNO[x,5],	header=T)	
	ID<-ANNOx[,8]
	Sx<-	ANNOx[,12]
	return(sapply(c(1:length(GEN[,1])),
		function(x){
		IDx<-GEN[x,3]
		nx<-GEN[x,4]
		return(sum(nchar(subset(Sx,ID== nx))))
	}))
})
	
COREs<-ifelse(COREn==1, COREs,NA)

#Normalize size for each core gene
COREsn<-(COREs)/apply(COREs,1,mean,na.rm=T)

#Remove gene copies with normalized size lower than 0.7 and higher than 1.3
COREnn<-ifelse(COREsn<0.7,0,
	ifelse(COREsn>1.3,0,COREn))
COREnn<-ifelse(is.na(COREnn)==T,0, COREnn)
	
#For each core gene, retrieve nucleotide sequences from each genome and align it with 184 genomes that where already analyzed

#Step 1: fasta file from Alessa et al 2021 30 new genomes + 184 genomes already analyzed
sapply(c(1:length(GEN[,1])),function(x){
	
	COREnnx<-COREnn[x,]
	vx<-subset(c(1:length(ANNO[,1])), COREnnx==1)
	Nx<-GEN[x,4]
	Ix<-GEN[x,3]
	print(Ix)
	Sx<-subset(paste(ANNO[,1], 
		ANNO[,2],sep="_"),COREnnx==1)
	setwd(paste(p0,"Annotations-16-11-2021",
		sep="/"))
	seqx<-sapply(vx,function(x){
		#print(ANNO[x,1])
		ANNOx<-read.delim(ANNO[x,5],	header=T)	
		return(list(unlist(strsplit(
			subset(ANNOx[,12],
			ANNOx[,8]== Nx),split=""))))
	})
	names(seqx)<-Sx	
	
	#retrieve unaligned nucleotide sequences from genomes already annotated
	ix<-unlist(strsplit(Ix,split="Gene"))[2]
	
	ix<-paste(c("Gene",
		seq(0,0,length.out=(4-nchar(ix))),ix),
		collapse="")
	
	setwd(pA)
	GENx<-read.fasta(paste(ix,
		".fas",sep=""))
	GENn<-names(GENx)
	
	seqn<-c(GENx, seqx)
	
	setwd(paste(p0,"New-Align",sep="/"))
	write.fasta(seqn,names=c(names(GENx), Sx),
		paste(Ix,"-New.fas",sep=""))
})

###Reformat strain names to avoid fordidden characters in different softwares

ST<-sort(unique(unlist(sapply(
	c(1:length(GEN[,1])),function(x){
	Ix<-GEN[x,3]
	setwd(paste(p0,"New-Align",sep="/"))
	seqx<-read.fasta(paste(Ix,"-New.fas",sep=""))
	return(names(seqx))
}))))

STn<-paste("st",c(1:length(ST)),sep="")
setwd(p0)
write.table(cbind(STn,ST),row.names=F,col.names=c("Simple-name","Strain-name"),"Strain-names.txt")

sapply(c(1:length(GEN[,1])),function(x){
	print(x)
	Ix<-GEN[x,3]
	setwd(paste(p0,"New-Align",sep="/"))
	seqx<-read.fasta(paste(Ix,"-New.fas",sep=""))
	sx<-names(seqx)
	snx<-as.vector(sapply(sx,function(x){
		return(subset(STn,ST==x))
	}))
	write.fasta(seqx,names= snx,
		paste(Ix,"-renammed.fas",sep=""))
})


#Step 2: align fasta files
library(msa)
library(seqinr)

#Alessa genomes
p0<-"/Users/jean-baptisteleducq/Desktop/New-Methylo-Genomes/Alessa-et-al-2021-genomes/"

setwd(p0)

#List of genomes
LIST<-as.vector(read.table("list-genomes.txt")[,1])

#List of core genes
GEN<-read.table("Core-Gene-REF=YIM132548-Cla=A5b.txt",header=T)

#Type of gene (coding(peg)/rna)
ANNO<-read.table("Genomes-to-annotate.txt",header=T)

setwd(paste(p0,"Annotations-16-11-2021",
	sep="/"))
ANNOx<-read.delim(ANNO[1,5],	header=T)

TYPE<-sapply(c(1:length(GEN[,1])),function(x){
	Ax<-GEN[x,4]
	return(unique(subset(ANNOx[,3], ANNOx[,8]==Ax)))
})	

#Check alignement size is a multiple of 3
setwd(paste(p0,"New-Align",sep="/"))

SUM3<-sapply(c(1:length(GEN[,1])),function(x){
	Ix<-GEN[x,3]
	print(Ix)
	#Get the unaligned nucleotide sequence
	SEQ<-read.fasta(paste(Ix,
		"-renammed.fas",sep=""))
	z<-unique(as.numeric(summary(SEQ)[,1]))
	return(unique((3*round(z/3))==z))
})

####Align non protein coding genes
Irna<-subset(c(1:length(GEN[,1])), TYPE=="rna")

setwd(paste(p0,"New-Align",sep="/"))

sapply(Irna,function(x){	
	Ix<-GEN[x,3]
	print(Ix)
	SEQw<-msa(paste(Ix,"-renammed.fas",sep=""),method="ClustalW",type="dna",order="input")
	#Retrieve aligned sequence in seqinr format
	SEQn<-msaConvert(SEQw,type="seqinr::alignment")$seq
	#Get sequence names from the unaligned nucleotide sequence
	NX=names(read.fasta(paste(Ix,
		"-renammed.fas",sep="")))
	#Convert
	SEQn<-sapply(c(1:length(SEQn)),function(x){
		return(strsplit(SEQn[x],split=""))
	})
	write.fasta(SEQn, names=NX,
		file.out =paste("New-",Ix,
		"_aligned.fas",sep=""))
})

####Align protein coding genes
Ipeg<-subset(c(1:length(GEN[,1])), TYPE=="peg")

setwd(paste(p0,"New-Align",sep="/"))

sapply(Ipeg,function(x){
	
	Ix<-GEN[x,3]
	print(Ix)

	#Get the unaligned nucleotide sequence
	SEQ<-read.fasta(paste(Ix,"-renammed.fas",sep=""))

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
		file.out =paste("New-",Ix,
		"_aligned.fas",sep=""))
})


#Review each protein-coding alignment and replace "NNN(---)NNN" motifs by "NNN(NNN)NNN"

setwd(paste(p0,"New-Align",sep="/"))

sapply(c(Ipeg,Irna),function(x){
	
	Ix<-GEN[x,3]
	print(Ix)

	#Get the aligned nucleotide sequence
	SEQ<-read.fasta(paste("New-",Ix,
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

setwd(paste(p0,"New-Align",sep="/"))

#only "n" should be returned

unique(unlist(sapply(c(1:length(GEN[,1])),function(x){
	
	Ix<-GEN[x,3]
	print(Ix)
	#Get the aligned nucleotide sequence (before and after N correcting)
	SEQ<-read.fasta(paste("New-",Ix,
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
setwd(paste(p0,"New-Align",sep="/"))

PPX<-unlist(sapply(Ipeg,function(x){
	
	Ix<-GEN[x,3]
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

setwd(paste(p0,"New-Align",sep="/"))

#Ony "n" should be returned
sumNEW<-t(sapply(Ipeg,function(x){
	
	Ix<-GEN[x,3]
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

sumNEWrna<-t(sapply(Irna,function(x){
	
	Ix<-GEN[x,3]
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

setwd(p0)
colnames(sumNEW)<-c("Name1","Name2","Sequences","Size1(bp)","Size2(bp)","Type")
write.table(sumNEW,"Summary-per-gene-after-alignment.txt",row.names=F)	

####STEP 3: format data for phylogenetic analyses

####A) RAxML concatenates tree
#A1) Need a concatenated alignement of the 384 core genes

#List of strains
setwd(p0)
ST<-read.table("Strain-names.txt",header=T)

#Retrieve missing genome for each gene (to replace missing sequences by NNNs in the concatenated alignement)
setwd(paste(p0,"New-Align",sep="/"))

OCC<-as.vector(sapply(c(1:length(GEN[,1])),function(x){	
	Ix<-GEN[x,3]
	print(Ix)
	#Get the aligned nucleotide sequence
	SEQn<-read.fasta(paste(Ix,
		"_alignedNF.fas",sep=""))
	return(setdiff(ST[,1],names(SEQn)))
}))			

#Create temporary alignemnts with "nn" sequences for missing data
sapply(c(1:length(GEN[,1])),function(x){	
	Ix<-GEN[x,3]
	print(Ix)
	#Get the aligned nucleotide sequence
	SEQn<-read.fasta(paste(Ix,
		"_alignedNF.fas",sep=""))
	sx=length(SEQn[[1]])
	occx<-c(OCC[[x]],0)	
	ns<-names(SEQn)
	Nx<-sapply(c(1:sx),function(x){
		return("n")
	})
	Nx<-sapply(c(1:length(occx)),function(x){
		return(list(Nx))
	})
	names(Nx)<-occx
	SEQn<-c(SEQn,Nx)
	SEQn<-SEQn[order(c(ns ,occx))]
	write.fasta(SEQn,paste(Ix,
		"_aligned-temp.fas",sep=""),
		names=sort(c(ns ,occx)))
})

#Construct the concatenated file

CONC<-sapply(ST[,1],function(x){
	STx=x
	print(STx)
	
	STn<-names(read.fasta(paste(Ix[1],
		"_aligned-temp.fas",sep="")))
	vx<-subset(c(1:length(STn)), STn ==STx)
	
	CONCx<-unlist(sapply(c(1:length(GEN[,1])),
	function(x){	
		Ix<-GEN[x,3]
		#print(Ix)
		#Get the aligned nucleotide sequence
		SEQn<-read.fasta(paste(Ix,
			"_aligned-temp.fas",sep=""))
		return(as.vector(SEQn[[vx]]))
	}))
	return(list(CONCx))
})

setwd(p0)

write.fasta(CONC,names= ST[,1],"Concatenated-genes-nov21.fas")

#A2) Need a partition file formated for IQ-tree; partitions will be reduced to group together genes with similar substitution models

#Generate a partition file for analysis in RAxML

#DNA, p1=1-30
#DNA, p2=31-60
	
PART<-matrix(unlist(sapply(c(1:length(GEN[,1])),
	function(x){
	Ix<-GEN[x,3]
	print(Ix)
	setwd(paste(p0,"New-Align",sep="/"))
	#Get the aligned nucleotide sequence
	SEQx<-read.fasta(paste(Ix,
		"_alignedNF.fas",sep=""))
	return(t(cbind(x,
		seq(1,1,length.out=
		length(SEQx[[1]])))))
})),ncol=2,byrow=T)

PART2<-cumsum(PART[,2])
PART<-sapply(unique(PART[,1]),function(x){
	return(paste("DNA,", GEN[x,3],"=",
	min(subset(PART2,PART[,1]==x)),"-",
	max(subset(PART2,PART[,1]==x)),sep=""))	
})

setwd(p0)

write.table(PART,
	"Partition-Concatenated-genes-nov21.txt",
	col.names=F,row.names=F,quote=F)


#A3) Partition finder in IQ-tree (on the server: IQ-tree-384-MERGER-rcluster.slurm)

#iqtree2 -s Concatenated-genes-nov21.fas -p Partition-Concatenated-genes-nov21.txt -m TESTMERGEONLY -rcluster 10

#Replace DNAF by DNA in the partition file (Partition-Concatenated-genes-nov21.best_scheme.txt)

#A4) ML tree in RAxML (GTRCAT model, 1,000 replicates; on the server): RAXML-GTRCAT-partitioning-n=384-1000.slurm

#mpirun raxmlHPC-MPI -f a -x 12345 -m GTRCAT -p 12345 -# 1000 -q Partition-Concatenated-genes-nov21.best_scheme.txt -s Concatenated-genes-nov21.fas -n Partition-Concatenated-genes-nov21-b=1000

#A5) Calculate average nucleotide identity from the contatenated alignement, in order to classify genomes in groups (similar to ANI), and construct a maping file

library(seqinr)

#Alessa genomes
p0<-"/Users/jean-baptisteleducq/Desktop/New-Methylo-Genomes/Alessa-et-al-2021-genomes/"

setwd(p0)

#List of genomes
LIST<-as.vector(read.table("list-genomes.txt")[,1])

#List of core genes
GEN<-read.table("Core-Gene-REF=YIM132548-Cla=A5b.txt",header=T)

#Type of gene (coding(peg)/rna)
ANNO<-read.table("Genomes-to-annotate.txt",header=T)

setwd(paste(p0,"Annotations-16-11-2021",
	sep="/"))
ANNOx<-read.delim(ANNO[1,5],	header=T)

TYPE<-sapply(c(1:length(GEN[,1])),function(x){
	Ax<-GEN[x,4]
	return(unique(subset(ANNOx[,3], ANNOx[,8]==Ax)))
})	

setwd(p0)

sumNEW<-read.table("Summary-per-gene-after-alignment.txt",header=T)
ST<-read.table("Strain-names.txt",header=T)
CONC<-read.fasta("Concatenated-genes-nov21.fas")

####Summary of core genes retrieved per genome

setwd(paste(p0,"New-Align",sep="/"))

SUMcore<-sapply(sumNEW[,2],function(x){
	print(x)
	seqx<-read.fasta(paste(x,
		"_aligned-temp.fas",sep=""))
	STn<-names(seqx)
	seqx<-sapply(c(1:length(seqx)),
	function(x){
		return(length(subset(seqx[[x]],
			seqx[[x]]=="n"))/length(seqx[[x]]))
	})
	return(sapply(ST[,1],function(x){
		return(round(subset(seqx, STn==x)))
	}))
})

setwd(p0)

write.table(cbind(ST, apply(SUMcore,1,sum)),"missing-core-gene-seq-per-genomes.txt")

nST<-length(ST[,1])

PX=c(1:length(CONC[[1]]))

ANI<-sapply(c(1:(nST-1)),function(x){
	X=x
	print(x)
	seq1<-as.vector(CONC[[X]])
	SUP1<-c(grep("n",seq1),grep("-",seq1))
	return(sapply(c((X+1): nST),function(x){
		seq2<-as.vector(CONC[[x]])
		SUP2<-c(grep("n",seq2),grep("-",seq2))
		seq12<-cbind(seq1,seq2)[c(setdiff(PX,
			c(SUP1, SUP2))),]
		return(length(subset(seq12[,1],
			seq12[,1]==seq12[,2]))/length(
			seq12[,1]))	
	}))
})

Pairs<-sapply(c(1:(nST-1)),function(x){
	X=x
	return(sapply(c((X+1): nST),function(x){
		return(paste(ST[X,1],ST[x,1],sep="_"))
	}))
})

ANI<-cbind(unlist(Pairs),unlist(ANI))

setwd(p0)

write.table(ANI,"ANI.txt",row.names=F,col.names=F,quote=F)

ANI<-read.table("ANI.txt")

hist(as.numeric(ANI[,2]),breaks=1000,
	main="",cex.axis=0.8,las=1,xlab="ANI")

TH0=0.83	
TH1=0.87
TH2=0.91
TH3=0.97
TH4=0.965 #to compare with TH3; threshold used in (Mende et al., 2013; Chun and Rainey)
	
segments(TH1,0, TH1,100000,col="red")
segments(TH2,0, TH2,100000,col="red")	
segments(TH3,0, TH3,100000,col="red")
segments(TH4,0, TH3,100000,col="red")
segments(TH0,0, TH0,100000,col="red")

###FUNCTIONS
library(parallel)

num.strsplit.u<-function(x){return(as.numeric(unlist(strsplit(x,split= "_"))))}

reduce.group<-function(groups,ncor) {
	reduce.gx<-function(i){
		gx<-num.strsplit.u(groups[i])
		z<-sapply(groups,function(x){
			return(length(intersect(
				gx,num.strsplit.u(x))))
		})
		return(paste(sort(unique(
			num.strsplit.u(names(subset(z,z!=0))))),
			collapse="_"))
	}
	return(unlist(mclapply(c(1:length(groups)), 
		reduce.gx,mc.cores=ncor)))
}

#Group strains with a different thresholds

TH<-c(TH0,TH1,TH2,TH3,TH4)

GROUPS<-sapply(TH,function(x){
	print(x)

	GROUPx<-subset(ANI[,1],
		as.numeric(ANI[,2])>=x)
	GROUPx <-gsub(pattern="st",
		replacement="", GROUPx)	

	GROUPx <-unique(reduce.group(
		unique(reduce.group(
		unique(reduce.group(
		unique(reduce.group(
		unique(reduce.group(
		unique(reduce.group(
		unique(reduce.group(
		unique(reduce.group(
	GROUPx,ncor=2)),ncor=2)),
		ncor=2)),ncor=2)),
		ncor=2)),ncor=2)),
		ncor=2)),ncor=2))

	GROUPx <-c(GROUPx,
		setdiff(c(1:length(ST[,1])),
		num.strsplit.u(GROUPx)))
	print(length(GROUPx))
	return(GROUPx)	
})

#For each threshold and each strain, report a group 

CLASS<-sapply(c(1:length(GROUPS)),function(x){
	Gx=x
	Nx=length(GROUPS[[x]])
	GPx<-matrix(unlist(sapply(c(1:Nx),function(x){
		return(t(cbind(x,num.strsplit.u(
			GROUPS[[Gx]][x]))))
	})),ncol=2,byrow=T)
	return(GPx[order(GPx[,2]),1])	
})

rownames(CLASS)<-ST[,1]
colnames(CLASS)<-paste("TH=",TH,sep="")

setwd(p0)

write.table(CLASS,"CLASSIFICATION-ANI.txt")

CLASS<-read.table("CLASSIFICATION-ANI.txt")

#Retrieve strain classification from previous studies
TAXA<-read.table("Table-Correspondance-ID.txt",
	header=T)

#Summary of classification in "species" with a 97% threshold on pairwise nucleotide similarity 
CLASSn<-t(sapply(unique(CLASS[,4]),function(x){
	TAXx<-subset(TAXA,CLASS[,4]==x)
	return(c(	
		length(TAXx[,1]), #number of genomes in this species
		paste(setdiff(unique(
			TAXx[,4]),"-"),collapse="_"), #clades in Leducq 2021
		paste(setdiff(unique(
			TAXx[,5]),"-"),collapse="_"), #clades in Alessa 2021
		paste(unique(
			TAXx[,6]),collapse="_"), #species names
		paste(unique(
			TAXx[,2]),collapse="_"), #strains in this species (phylogeny format)
		paste(unique(
			TAXx[,3]),collapse=";") #strains in this species (actual names)	
	))	
}))

#Calculate average, sd and range of pairwise nucleotide diversity per species (only in species with more than one strain)

gx<-subset(c(1:length(CLASSn[,1])), as.numeric(CLASSn[,1])>1)

ANIs<-t(sapply(gx,function(x){
	stx<-unlist(strsplit(CLASSn[x,5],split="_"))
	nx<-length(stx)
	px<-unlist(sapply(c(1:(nx-1)),function(x){
		X=x
		sapply(c((X+1):(nx)),function(x){
			return(paste(stx[X],stx[x],sep="_"))
		})	
	}))
	anix<-sapply(px,function(x){
		return(subset(ANI[,2],ANI[,1]==x))
	})
	return(c(mean(anix),
		sd(anix),
		min(anix),
		max(anix)))
}))

CLASSs<-subset(CLASSn,CLASSn[,1]!="1")
CLASSn<-subset(CLASSn,CLASSn[,1]=="1")
CLASSs<-cbind(CLASSs,ANIs)

ANIn<-seq(1,1,length.out=length(CLASSn[,1]))
ANIn<-cbind(ANIn, ANIn*NA, ANIn, ANIn)
CLASSn<-cbind(CLASSn,ANIn)

CLASSn<-rbind(CLASSn,CLASSs)

CLASSn<-CLASSn[order(CLASSn[,4],decreasing=T),]
CLASSn<-CLASSn[order(CLASSn[,3]),]
CLASSn<-CLASSn[order(CLASSn[,2]),]

colnames(CLASSn)<-c("Genomes","Leducq2021","Alessa2021","Species","NamePhylo","NameStrains","AvPS","SdPS","MinPS","MaxPS")

#Report biome of origine for each species
BIOME<-sort(unique(unlist(
	strsplit(TAXA[,7],split="/"))))
	
BIOMEf<- t(sapply(c(1:length(CLASSn[,1])),function(x){
	stx<-as.vector(unlist(strsplit(
		CLASSn[x,5],split="_")))
	bx<-matrix(unlist(sapply(stx,function(x){
		bxx<-unlist(strsplit(subset(TAXA[,7],
			TAXA[,2]==x),split="/"))
		fxx<-seq(1,1,length.out=
			length(bxx))
		return(t(cbind(bxx, fxx)))			
	})),ncol=2,byrow=T)
	
	return(sapply(BIOME,function(x){
		return(sum(subset(as.numeric(bx[,2]),
			bx[,1]==x))/length(stx))
	}))
}))	

setwd(p0)

write.table(cbind(CLASSn, BIOMEf),"CLASSIFICATION-ANI2.txt",row.names=F)

#Convert CLASSn in Excell to Mapping-ASTRAL.txt; see https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md section "Running on a multi-individual datasets" for formating


####B) ASTRAL tree
#B1) Need a ML tree for each of the 384 core genes: ML tree in RAxML (GTRGAMMA or GRTCAT? model, 1,000 replicates per tree) ; use GeneXXX_alignedNF.fas files (on the server)
#B2) combine all tree in the same tree file (see how to do in script format-tree.R)
#B2) ASTRAL

#package ape for trees
#install.packages("https://cran.r-project.org/bin/macosx/contrib/4.0/ape_5.5.tgz")

#Package ips for collapseUnsupportedEdges
#install.packages("https://cran.r-project.org/bin/macosx/contrib/4.1/ips_0.0.11.tgz")

#ips dependancies
#install.packages("https://cran.r-project.org/bin/macosx/contrib/4.1/XML_3.99-0.8.tgz")
#install.packages("https://cran.r-project.org/bin/macosx/contrib/4.1/phangorn_2.8.0.tgz")
#install.packages("https://cran.r-project.org/bin/macosx/contrib/4.1/fastmatch_1.1-3.tgz")
#install.packages("https://cran.r-project.org/bin/macosx/contrib/4.1/igraph_1.2.8.tgz")
#install.packages("https://cran.r-project.org/bin/macosx/contrib/4.1/pkgconfig_2.0.3.tgz")
#install.packages("https://cran.r-project.org/bin/macosx/contrib/4.1/quadprog_1.5-8.tgz")
#install.packages("https://cran.r-project.org/bin/macosx/contrib/4.1/plyr_1.8.6.tgz")

library(ape)
library(ips)

p0<-"/Users/jean-baptisteleducq/Desktop/New-Methylo-Genomes/Alessa-et-al-2021-genomes/"

setwd(paste(p0,"For-RAxML-gene-trees",sep="/"))

LISTt<-as.vector(read.table("list-tree.txt")[,1])

#Bootstrapp support are in brackets and must be extracted and modified to be readable by ape

mBS=10 #minimum node support to keep node

#Raxml best trees

TREE<-sapply(LISTt,function(x){
	print(x)
	#read tree as character
	treex<-as.vector(read.table(paste(
		"RAxML_bipartitionsBranchLabels.",
		x,"_fa-#1000_mGTRGAMMA",sep="")))
	#split nodes
	treex<-unlist(strsplit(
		as.character(treex),split=":"))
	
	#identify nodes with brackets (node supports)
	bi<-grep("[[]",treex)
	bo<-setdiff(c(1:length(treex)), bi)
	#extract node support and convert in readable tree format: #BranchLength[NodeSupport])# converted in #NodeSupport:BranchLength)# 
	ns<-sapply(bi,function(x){
		nx<-unlist(strsplit(treex[x],split="[[]"))
		nx<-unlist(strsplit(nx,split="[]]"))
		return(paste(nx[2],":",nx[1],nx[3],sep=""))
	})
	treex<-c(ns, treex[bo])[order(c(bi,bo))]
	
	#if last sign is not ")" or ";", add ":"
	bi<-grep("[)]", treex)
	bo<-setdiff(c(1:length(treex)), bi)
	
	ns<-paste(treex[bo],":",sep="")
	
	treex<-c(ns, treex[bi])[order(c(bo,bi))]
	
	treex<-paste(treex,collapse="")
	write.table(treex,
		paste("RAxML_tree.",x,"_reformated",sep=""),
		row.names=F,col.names=F,quote=F)
	
	treex<-read.tree(paste("RAxML_tree.",
		x,"_reformated",sep=""))	
	
	treex <-collapseUnsupportedEdges(treex, 
		value = "node.label", cutoff= mBS)	
	write.tree(treex,paste("RAxML_tree.",
		x,"_reformated",sep=""))
	return(list(unroot(treex)))
})

class(TREE) <- "multiPhylo"

write.tree(TREE,"Multi-coregene.tre")

#Then in shell, compute lineage tree with Astral (-Xmx3000M to allow the use of more memory by java)

cd Desktop/
cd myRAST-out
cd Astral.5.7.8

java -Xmx3000M -jar astral.5.7.8.jar -i With-genomes-Alessa-nov-2021/Multi-coregene.tre -o With-genomes-Alessa-nov-2021/Astral-lineage-tree.tre  

#Also in shell, compute species tree with Astral (-Xmx3000M to allow the use of more memory by java)

cd Desktop/
cd myRAST-out
cd Astral.5.7.8

java -Xmx3000M -jar astral.5.7.8.jar -i With-genomes-Alessa-nov-2021/Multi-coregene.tre -a With-genomes-Alessa-nov-2021/Mapping-ASTRAL.txt -o With-genomes-Alessa-nov-2021/Astral-species-tree.tre  

####C) SVDquartet lineage and species tree tree
#Needs the concatenated alignement (A1) ===> DAVID


####STEP 4: Calculate RF distances (OR DO THIS IN PAUP)

#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/TreeDist_2.2.0.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/Rdpack_2.1.2.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/rbibutils_2.2.4.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/TreeTools_1.5.1.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/R.cache_0.15.0.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/R.methodsS3_1.8.1.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/R.oo_1.24.0.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/R.utils_2.11.0.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/digest_0.6.28.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/bit64_4.0.5.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/bit_4.0.4.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/colorspace_2.0-2.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/shiny_1.7.1.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/htmltools_0.5.2.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/httpuv_1.6.3.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/later_1.3.0.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/promises_1.2.0.1.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/mime_0.12.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/xtable_1.8-4.tgz")
#install.packages("https://cran.rstudio.com/bin/macosx/contrib/4.1/shinyjs_2.0.0.tgz")

p0<-"/Users/jean-baptisteleducq/Desktop/New-Methylo-Genomes/Alessa-et-al-2021-genomes/"

library(ape)
library(TreeDist)

setwd(paste(p0,"For-RAxML-concatenated",sep="/"))

RAxML<-read.tree("RAxML_bipartitionsBranchLabels.Partition-Concatenated-genes-nov21-b=500") #best RAxML tree from concatenated alignement (500 replicates)

MULTIrax<-read.tree("RAxML_bootstrap.Partition-Concatenated-genes-nov21-b=500") #(replicates)

#MULTIrax<-read.tree("RAxML_bootstrap.Partition-Concatenated-genes-nov21-b=1000") #bootstrapp replicates from the RAxML search on contatenated alignments

#MJraxml<-read.tree("RAxML_bootstrap.Partition-Concatenated-genes-nov21-b=1000-consensus-MJ50-PAUP2") #Majority rule consensus tree estimated in PAUP from the 1,008 RAxML replicated trees 

#IQtree<-read.tree("Partition-Concatenated-genes-nov21.txt.treefile") #tree from IQtree (partition search)

#DISTiqrep<-TreeDistance(IQtree, MULTIrax)

#DISTraxrep<-TreeDistance(MJraxml, MULTIrax)

#Astral trees: carreful: nodal support are in local porterior probability 

setwd(paste(p0,"For-RAxML-gene-trees",sep="/"))

#MULTI<-read.tree("Multi-coregene.tre") #individual RAxML gene trees combined in a single file
ASTRAL<-read.tree("Astral-lineage-tree.tre") #lineage tree
ASTRALspe<-read.tree("Astral-species-tree.tre") #species tree

#SVDquartet trees (carefull: nodal support were included in branch length, need to reformat the trees)
setwd(paste(p0,"For-SVDquartet",sep="/"))
SVD<-read.tree("Concat_n=213_no_map.nwk") #lineage tree
SVDspe<-read.tree("Concat_n=213_map.nwk") #species tree

SVD$node.label<-SVD$edge.length
SVD$edge.length<-NULL
SVDspe $node.label<-SVDspe $edge.length
SVDspe $edge.length<-NULL

write.tree(SVDspe,"Concat_n=213_map-reformated.nwk")
write.tree(SVD,"Concat_n=213_no_map-reformated.nwk")


#RF distances between SVDq and ASTRAL trees
#TreeDistance(SVDspe, ASTRALspe)
#RFastsvd <- TreeDistance(ASTRAL, SVD)
#RFastrax <-TreeDistance(ASTRAL, RAxML)
#RFsvdrax <-TreeDistance(RAxML, SVD)

#setwd(p0)

#pdf("RF-dist-213-genomes.pdf",height=4,width=4)

#hist(DISTraxrep,breaks=50,border=NA,
#	col="grey",xlim=c(0,
#	max(c(RFsvdrax, RFastrax, RFastsvd))),
#	main="213 genomes",
#	cex.axis=0.8,las=1,xlab="RFnorm")
	
#points(RFsvdrax,0,pch=19,cex=1,col="blue")
#points(RFastrax,0,pch=19,cex=1,col="green2")
#points(RFastsvd,0,pch=19,cex=1,col="red")

#legend("topright",cex=0.5, horiz=F,
#	legend= c("RAxML conc. vs. rep.",
#	"RAxML conc. vs. SVDq.",
#	"RAxML conc. vs. Astral",
#	"SVDq. vs. Astral"),
#	text.col= "black",box.col=NA,pch=c(15,19,19,19),
#	border=NA,col= c("grey","blue","green2","red"),
#	title="",title.col="black",ncol=1)

#dev.off()

###Generate multitree files for RF distance calculation in PAUP

#1 For RF calculation between RAxML replicate trees and the RAxML majority rule consensus tree

RAxMLtrees<-c(list(RAxML),MULTIrax)

class(RAxMLtrees) <- "multiPhylo"

#2 For RF calculation between the tree lineage trees (RAxML majority rule consensus tree; ASTRAL, SVDquartet) 

LINtrees<-c(list(RAxML),list(ASTRAL),list(SVD))

class(LINtrees) <- "multiPhylo"

#3 Ror RF calculation between the two species trees

SPEtrees<-c(list(ASTRALspe),list(SVDspe))

class(SPEtrees) <- "multiPhylo"

setwd(paste(p0,"Trees-RF",sep="/"))

write.tree(LINtrees,"Lineage-trees-RAS")
write.tree(RAxMLtrees,"Lineage-RAxML-replicates")
write.tree(SPEtrees,"Species-trees-AS")

#Open these trees in PAUP 
#Compute tree-to-tree distances
### for "Lineage-RAxML-replicates"
#Choose tree 1 to all other trees (RAxMLtrees)
#Symmetric differences (RFx2)
#Save distances to filesLineage-RAxML-replicates-treedist.txt

### for "Lineage-trees-RAS"
#Choose all pairwise comparisons (LINtrees)
#Symmetric differences (RFx2)
#Save distances to Lineage-trees-RAS-treedist.txt

#Open distance files (after manual space removing and reformating)

LINtreesRF<-read.table("Lineage-trees-RAS-treedist.txt",header=T)
REPtreesRF<-read.table("Lineage-RAxML-replicates-treedist.txt",header=T)

NORM=2*length(SVD$tip.label)

pdf("RF-dist-213-genomes-PAUP.pdf",height=4,width=4)

hist(REPtreesRF[,2]/NORM,breaks=20,border=NA,
	col="grey",xlim=c(0,
	max(LINtreesRF[,3])/NORM),
	main="",
	cex.axis=0.8,las=1,xlab="Norm. RF distance")
	
points(LINtreesRF[1,3]/NORM,0,pch=19,cex=1,col="blue")
points(LINtreesRF[2,3]/NORM,0,pch=19,cex=1,col="green2")
points(LINtreesRF[3,3]/NORM,0,pch=19,cex=1,col="red")

legend("topright",cex=0.7, horiz=F,
	legend= c("RAxML conc. vs. replicates",
	"RAxML conc. vs. Astral",
	"RAxML conc. vs. SVDquartet",
	"SVDquartet vs. Astral"),
	text.col= "black",box.col=NA,pch=c(15,19,19,19),
	border=NA,col= c("grey","blue","green2","red"),
	title="",title.col="black",ncol=1)

dev.off()

LINtreesRF[,3]/NORM
# 0.1807512 0.2253521 0.2887324
mean(REPtreesRF[,2]/NORM)
# 0.06924883
 sd(REPtreesRF[,2]/NORM)
 # 0.01524204
 min(REPtreesRF[,2]/NORM)
# 0.02816901
 max(REPtreesRF[,2]/NORM)
# 0.1126761

###Reformat and rename lineage trees tip labels with actual strain names and species name when available in the RAxML tree

setwd(paste(p0,"For-RAxML-concatenated",sep="/"))

	#read tree as character
	treex<-as.vector(read.table("RAxML_bipartitionsBranchLabels.Partition-Concatenated-genes-nov21-b=500"))
	#split nodes
	treex<-unlist(strsplit(
		as.character(treex),split=":"))
	
	#identify nodes with brackets (node supports)
	bi<-grep("[[]",treex)
	bo<-setdiff(c(1:length(treex)), bi)
	#extract node support and convert in readable tree format: #BranchLength[NodeSupport])# converted in #NodeSupport:BranchLength)# 
	ns<-sapply(bi,function(x){
		nx<-unlist(strsplit(treex[x],split="[[]"))
		nx<-unlist(strsplit(nx,split="[]]"))
		return(paste(nx[2],":",nx[1],nx[3],sep=""))
	})
	treex<-c(ns, treex[bo])[order(c(bi,bo))]
	
	#if last sign is not ")" or ";", add ":"
	bi<-grep("[)]", treex)
	bo<-setdiff(c(1:length(treex)), bi)
	
	ns<-paste(treex[bo],":",sep="")
	
	treex<-c(ns, treex[bi])[order(c(bo,bi))]
	
	treex<-paste(treex,collapse="")
	write.table(treex,"RAxML_bipartitionsBranchLabels.Partition-Concatenated-genes-nov21-b=500",
		row.names=F,col.names=F,quote=F)
	
	treex<-read.tree("RAxML_bipartitionsBranchLabels.Partition-Concatenated-genes-nov21-b=500")	


setwd(p0)

ID<-read.table("Table-Correspondance-ID.txt",header=T)

IDr<-treex$tip.label

IDr <-as.vector(sapply(IDr,function(x){
	nx<-subset(ID[,c(6,3)],ID[,2]==x)
	return(paste(nx,collapse="_"))
}))

RAxMLr<-treex
RAxMLr $tip.label<-IDr

setwd(paste(p0,"Trees-RF",sep="/"))

write.tree(RAxMLr,"RAxML.renammed-tree")

#RF distances between the ASTRAL lineage tree and gene trees
#RFastral<-sapply(c(1:length(MULTI)),function(x){
#	z<-MULTI[[x]]$tip.label
#	treex<-keep.tip(ASTRAL,z)
#	return(TreeDistance(MULTI[[x]], treex))
#})

#RF distances between the RAXML lineage tree and gene trees
#RFiqtree<-sapply(c(1:length(MULTI)),function(x){
#	z<-MULTI[[x]]$tip.label
#	treex<-keep.tip(IQtree,z)
#	return(TreeDistance(MULTI[[x]], treex))
#})

#RF distances between the SVDq lineage tree and gene trees
#RFsvd<-sapply(c(1:length(MULTI)),function(x){
#	z<-MULTI[[x]]$tip.label
#	treex<-keep.tip(SVD,z)
#	return(TreeDistance(MULTI[[x]], treex))
#})

	
#hist(RFastral,breaks=50,
#	col=rgb(1,0,0,0.2),border=NA,lwd=0.5,
#	main="",las=1,cex.axis=0.8,xlim=c(0,1),
#	ylim=c(0,70))
#par(new=T)	
#hist(RFiqtree,breaks=50,
#	col=rgb(0,1,0,0.2),border=NA,lwd=0.5,
#	main="",las=1,cex.axis=0.8,xlim=c(0,1),
#	ylim=c(0,70))
#par(new=T)	
#hist(RFsvd,breaks=50,
#	col=rgb(0,0,1,0.2),border=NA,lwd=0.5,
#	main="",las=1,cex.axis=0.8,xlim=c(0,1),
#	ylim=c(0,70))


#######STEP 5: ANNOTATIONS

### Step 5a: annotation summary
#Retrieve annotations for all genomes

library(seqinr)

#Alessa genomes
p0<-"/Users/jean-baptisteleducq/Desktop/New-Methylo-Genomes/Alessa-et-al-2021-genomes/"

setwd(p0)

#List of genomes
LIST<-as.vector(read.table("list-genomes.txt")[,1])
#Description of genomes
ANNO<-read.table("Genomes-to-annotate.txt",header=T)

#Previously annotated genomes
#MyRAST annotations (see script Analysis-Gene-annotation.R)
pS="/Users/jean-baptisteleducq/Dropbox/PostDoc-Moscow-ID/Project/myRAST-gene_annotation"

setwd(pS)

#Core gene names
ANNOf<-read.table("Summary-Genes-After-GapNs-filtering.txt",header=T) # after filtering
#gene names (retrieve short names from original set)

GENES<-read.table("Gene-modified-names.txt",header=T)


#Genome description
NAMES<-read.table("Genome-description.txt",header=T)

#Annotations in Alessa genomes

setwd(paste(p0,"Annotations-16-11-2021",
		sep="/"))
		
#Number of scaffold with annotation and scaffold size per genome		

Sa<-t(sapply(c(1:length(ANNO[,1])),
	function(x){
	print(ANNO[x,1])
	ANNOx<-read.delim(ANNO[x,5],	header=T)
	GCu<-summary(as.factor(unlist(
		strsplit(ANNOx[,12],split=""))))
	GCu<-(subset(GCu,names(GCu)=="c")+subset(GCu,names(GCu)=="g"))/sum(GCu)	
	SCAu<-sort(unique(ANNOx[,1]))	
	SCAs<-sapply(SCAu,function(x){
		sx<-subset(ANNOx[,c(5:6)],ANNOx[,1]==x)
		return(max(sx))
	})			
	return(c(length(SCAu),
		sum(SCAs),max(SCAs), GCu))
}))

#Unique annotation
Fua<-sort(unique(unlist(
	sapply(c(1:length(ANNO[,1])),
	function(x){
	print(ANNO[x,1])
	ANNOx<-read.delim(ANNO[x,5],	header=T)	
	ID<-ANNOx[,8]	
	return(unique(ID))
}))))

#Absolute gene abundances
Fra<-sapply(c(1:length(ANNO[,1])),
	function(x){
	print(ANNO[x,1])
	ANNOx<-read.delim(ANNO[x,5],	header=T)	
	ID<-ANNOx[,8]	
	return(summary(as.factor(ID),
		maxsum=length(Fua)))
})

#Annotations in previously annotated genomes



pR<-"/Users/jean-baptisteleducq/Desktop/myRAST-out"

setwd(pR)

Su<-t(sapply(c(1:length(NAMES[,1])),
	function(x){
	print(NAMES[x,2])
	ANNOx<-read.delim(rownames(NAMES)[x],header=T)
	GCu<-summary(as.factor(unlist(
		strsplit(ANNOx[,12],split=""))))
	GCu<-(subset(GCu,names(GCu)=="c")+subset(GCu,names(GCu)=="g"))/sum(GCu)	
		
	SCAu<-sort(unique(ANNOx[,1]))	
	SCAs<-sapply(SCAu,function(x){
		sx<-subset(ANNOx[,c(5:6)],ANNOx[,1]==x)
		return(max(sx))
	})			
	return(c(length(SCAu),
		sum(SCAs),max(SCAs), GCu))
}))


#Unique annotation
Fu<-sort(unique(unlist(
	sapply(c(1:length(NAMES[,1])),
	function(x){
	print(NAMES[x,2])
	ANNOx<-read.delim(rownames(NAMES)[x],header=T)	
	ID<-ANNOx[,8]	
	return(unique(ID))
}))))

#Absolute gene abundances
Fr<-sapply(c(1:length(NAMES[,1])),
	function(x){
	print(NAMES[x,2])
	ANNOx<-read.delim(rownames(NAMES)[x],header=T)	
	ID<-ANNOx[,8]	
	return(summary(as.factor(ID),
		maxsum=length(Fu)))
})

#Merge abundance from both datasets
Fu<-unique(c(Fu,Fua)) #annotations
Fr<-c(Fr,Fra) #annotations abundances
ST<-c(NAMES[,2],paste(ANNO[,1],ANNO[,2],sep="_")) #strain names
ANx<-c(rownames(NAMES),ANNO[,5]) #annotation files
Su<-rbind(Su, Sa) #number of scaffolds, scaffold size and genome

setwd(p0)
STn<-read.table("Strain-names.txt",header=T)
STn<-as.vector(sapply(ST,function(x){
	return(subset(STn[,1],STn[,2]==x))
})) #reformated strain names for phylogenies

#Reformate gene abundance in a matrix
Fr<-sapply(c(1:length(ST)),
	function(x){
	print(ST[x])
	Fx<-Fr[[x]]	
	F0<-setdiff(Fu,names(Fx))
	Fn<-seq(0,0,length.out=length(F0))
	names(Fn)<-F0
	Fx<-c(Fx,Fn)
	return(Fx[order(names(Fx))])
})

colnames(Fr)<-STn
Fu<-rownames(Fr)

#Summary per gene 

SUMgene<-t(sapply(c(1:length(Fu)),function(x){
	#Long name of this gene (as in annotation raw files from myRAST)
	Fux=Fu[x]
	print(x)
	#Is this gene was in the original set? if yes, report short name, otherwise report "NA"
	GENESx<-c(subset(GENES[,1],GENES[,2]== Fux),NA)[1]
	#Convert this name as it should appear in alignments and RAxML gene tree for core genes
	Nx<-paste("Gene",
		setdiff(as.numeric(unlist(strsplit(
		GENESx,split="Gene"))),NA),sep="")
	Nx<-ifelse(is.na(GENESx)==T,NA, Nx)	
	#Is this gene is one of the 384 core genes? if yes, report short name, otherwise report "NA"
	COREx<-c(subset(ANNOf[,2],ANNOf[,2]==Nx),NA)[1]
	#Abundance of this gene accross genomes
	Fx<-as.vector(subset(Fr,Fu==Fux))
	#Number of genomes were this gene is present
	Px<-length(subset(Fx,Fx!=0))
	#Average abundance of this gene accross genomes
	Mx=mean(Fx)
	#SD abundance of this gene accross genomes
	SDx<-sd(Fx)
	#Average abundance of this gene accross genomes where it is present
	MPx=mean(subset(Fx,Fx!=0))
	#SD abundance of this gene accross genomeswhere it is present
	SDPx<-sd(subset(Fx,Fx!=0))
	
	return(c(Fux, GENESx, COREx, Px, Mx, SDx, MPx, SDPx))
				
}))

colnames(SUMgene)<-c("Annotation","Short-Name","Core-gene","Occurence","Mean-abundance","SD-abundance","Mean-abundance-if-present","SD-abundance-if-present")

#Summary per genome

SUMgenome<-t(sapply(c(1:length(ST)),function(x){
	print(x)
	STx<-ST[x]
	STnx<-STn[x]
	Fx<-Fr[,x]
	anx<-ANx[x]
	Sux<-Su[x,]
	#Number of annotation, including duplications
	Sx=sum(Fx)
	#Number of "hypothetical protein"
	HPx<-subset(Fx,names(Fx)=="hypothetical protein")
	#Number of "repeat region"
	RRx<-subset(Fx,names(Fx)=="repeat region")	
	#Number of "Mobile element protein"
	MEx<-subset(Fx,names(Fx)=="Mobile element protein")
	#Remove these abundant elements from Fx
	Fx<-subset(Fx,names(Fx)!="hypothetical protein")
	Fx<-subset(Fx,names(Fx)!="repeat region")
	Fx<-subset(Fx,names(Fx)!="Mobile element protein")
	#Remove absent annotations
	Fx<-subset(Fx,Fx!=0)
	#Number of unique annotations after removing HP, RR and ME
	Lx=length(Fx)
	#Sum of unique annotations after removing HP, RR and ME
	SLx=sum(Fx)
	#Mean and SD copy number
	CNx<-c(mean(Fx),sd(Fx))
	
	return(c(STx,STnx, anx, Sux[1:3] ,Sx, HPx, RRx, MEx, Lx, SLx, CNx,Sux[4]))
	
}))

colnames(SUMgenome)<-c("Strain","Strain-simp","annotation-file","scaffolds","genome-size","max-scaffold-size","annotations","hypothetical-protein","repeat-region","Mobile-element-protein","other-annotations-unique","other-annotations-sum","other-annotations-mean-copy-number","other-annotations-sd-copy-number","GC-content")

setwd(paste(p0,"Gene-content",sep="/"))

write.table(SUMgenome,"Annotation-Summary-per-genome.txt",row.names=F)
write.table(SUMgene,"Annotation-Summary-per-gene.txt",row.names=F)
write.table(Fr,"Annotation-Abundance.txt")

#Convert annotation in 01 matrix and then in fasta file for phylogenetic analyses (1: gene present; 0: gene absent)

library(seqinr)

SEQan<-sapply(c(1:length(ST)),function(x){
	Fx<-as.vector(Fr[,x])
	Fx<-ifelse(Fx==0,0,1)
	return(Fx)
})

#Remove annotations present in all genomes (unimformative)
SEQan<-subset(SEQan,apply(SEQan,1,sum)!=length(ST))

SEQan<-sapply(c(1:length(ST)),function(x){
	Fx<-as.vector(SEQan[,x])
	Fx<-ifelse(Fx==0,0,1)
	return(list(Fx))
})

write.fasta(SEQan,names=STn,"Annotation-occurence.fas")


#####Step 5b: reformat ANNOTATION tree en calculate RF distances


p0<-"/Users/jean-baptisteleducq/Desktop/New-Methylo-Genomes/Alessa-et-al-2021-genomes/"

setwd(paste(p0,"Gene-content",sep="/"))

#Retrieve phylogenetic tree based on annotation occurence
	#read tree as character
	treex<-as.vector(read.table("RAxML_bipartitionsBranchLabels.Annotation-occurence-1000_mBINCAT"))
	#split nodes
	treex<-unlist(strsplit(
		as.character(treex),split=":"))
	
	#identify nodes with brackets (node supports)
	bi<-grep("[[]",treex)
	bo<-setdiff(c(1:length(treex)), bi)
	#extract node support and convert in readable tree format: #BranchLength[NodeSupport])# converted in #NodeSupport:BranchLength)# 
	ns<-sapply(bi,function(x){
		nx<-unlist(strsplit(treex[x],split="[[]"))
		nx<-unlist(strsplit(nx,split="[]]"))
		return(paste(nx[2],":",nx[1],nx[3],sep=""))
	})
	treex<-c(ns, treex[bo])[order(c(bi,bo))]
	
	#if last sign is not ")" or ";", add ":"
	bi<-grep("[)]", treex)
	bo<-setdiff(c(1:length(treex)), bi)
	
	ns<-paste(treex[bo],":",sep="")
	
	treex<-c(ns, treex[bi])[order(c(bo,bi))]
	
	treex<-paste(treex,collapse="")
	write.table(treex,"RAxML_bipartitionsBranchLabels.Annotation-occurence-1000_mBINCAT_reformated",
		row.names=F,col.names=F,quote=F)
	
	library(ape)
	treex<-read.tree("RAxML_bipartitionsBranchLabels.Annotation-occurence-1000_mBINCAT_reformated")	

#Retrieve species names

setwd(p0)

SPE <-read.table("Species-tree-order.txt",header=T)
LIN <-read.table("Lineage-tree-order.txt",header=T)

#check that each species is monophyletic in the annotation tree
MONOspe<-sapply(c(1:length(SPE[,1])),function(x){
	stx<-unlist(strsplit(
		SPE[x,3],split='_'))
	is.monophyletic(treex, stx)
})

#check that each group is monophyletic in the annotation tree
MONOgroup<-sapply(unique(SPE[,4]),function(x){
	stx<-unlist(strsplit(
		subset(SPE[,3],
		SPE[,4]==x),split='_'))
	is.monophyletic(treex, stx)
})


#Compare RF distance
# 1 between the annotation BINCAT best tree and replicates
# 2 between the annotation BINCAT best tree and lineage trees

setwd(paste(p0,"Gene-content",sep="/"))

REP<-read.tree("RAxML_bootstrap.Annotation-occurence-1000_mBINCAT")

REP<-c(list(treex), REP)

class(REP) <- "multiPhylo"

setwd(paste(p0,"Trees-RF",sep=""))

write.tree(REP,"Gene-content-RAxML-replicates.tre")

LINtree<-read.tree("Lineage-trees-RAS")

LINtree<-c(list(treex),LINtree)

class(LINtree) <- "multiPhylo"

write.tree(LINtree,"Lineage-trees-RAS-Gene-content")

#Get tree distances (replace scales by underscores in files)

RFrep<-read.table("Gene-content-RAxML-replicates-treedist.txt",header=T)
RFlin<-read.table("Lineage-trees-RAS-Gene-content-treedist.txt",header=T)

mean(RFrep[,2]/NORM)
sd(RFrep[,2]/NORM)
min(RFrep[,2]/NORM)
max(RFrep[,2]/NORM)
RFlin[,2]/NORM

NORM<-length(treex$tip.label)*2

pdf("RF-dist-gene-content-213-genomes-PAUP.pdf",height=3,width=4)

par(mar=c(4,4,1,1),bty="n")

hist(RFrep[,2]/NORM,breaks=20,border=NA,
	col="grey",xlim=c(min(RFrep[,2]/NORM),
	max(RFlin[,2])/NORM),
	main="",
	cex.axis=0.8,las=1,xlab="Norm. RF distance")
	
points(RFlin[,2]/NORM,c(0,0,0),pch=19,
	cex=1,col=c("blue","green2","red"))

legend("topright",cex=0.5, horiz=F,
	legend= c("Gene content vs. replicates",
	"RAxML conc. vs. Gene content",
	"Astral vs. Gene content",
	"SVDquartet vs. Gene content"),
	text.col= "black",box.col=NA,pch=c(15,19,19,19),
	border=NA,col= c("grey","blue","green2","red"),
	title="",title.col="black",ncol=1)

dev.off()

#####Step 5c1: Stat on Genome characteristics per group 

p0<-"/Users/jean-baptisteleducq/Desktop/New-Methylo-Genomes/Alessa-et-al-2021-genomes/"

setwd(paste(p0,"Gene-content",sep="/"))

#Retrieve genome summaries
SUMgenome<-read.table("Annotation-Summary-per-genome.txt",header=T)

#Color code for groups
COLgroup<-cbind(
	c("A","B","C","D",
	"Enterovirga","Microvirga"),
	c(rgb(1,0,0,0.1),
	rgb(0.8,0,1,0.1),
	rgb(0,0.8,0,0.1),
	rgb(0,0.5,1,0.1),
	rgb(0.7,0.7,0.7,0.1),
	rgb(0.8,0.8,0.8,0.1)),
	c(rgb(1,0,0,0.4),
	rgb(0.8,0,1,0.4),
	rgb(0,0.8,0,0.4),
	rgb(0,0.5,1,0.4),
	rgb(0.7,0.7,0.7,0.4),
	rgb(0.8,0.8,0.8,0.4)),
	c(rgb(1,0,0,1),
	rgb(0.8,0,1,1),
	rgb(0,0.8,0,1),
	rgb(0,0.5,1,1),
	rgb(0.7,0.7,0.7,1),
	rgb(0.8,0.8,0.8,1)))
	
COLgroup<-	COLgroup[c(1,4,2,3,5,6),]


##Genome summary per group

setwd(p0)

SPE <-read.table("Species-tree-order.txt",header=T)
SPE<-SPE[order(SPE[,1]),]

ID<-read.table(
	"Table-Correspondance-ID.txt",header=T)

STx<-SUMgenome[,2]

GroupLIN<-sapply(STx,function(x){
	stx<-x
	return(SPE[grep(paste(stx,"_",sep=""),
		paste(SPE[,3],"_",sep="")),4])
})

SpeLIN<-sapply(STx,function(x){
	stx<-x
	return(SPE[grep(paste(stx,"_",sep=""),
		paste(SPE[,3],"_",sep="")),2])
})

BIOME<-sapply(STx,function(x){
	return(subset(ID[,7],ID[,2]==x))
})

SUMgroup<-t(sapply(sort(unique(GroupLIN)),function(x){
	sumx<-matrix(as.numeric(as.matrix(
		as.data.frame(subset(SUMgenome, 
		GroupLIN==x))[c(5,7,11,12,13,15)])),
		ncol=6)
	nSPE<-length(unique(subset(SpeLIN, GroupLIN==x)))
		
	BIOMEx<-subset(BIOME,GroupLIN==x)	
	return(c(
		length(sumx[,1]), nSPE,
		paste(round(mean(sumx[,1]/10^6),2),
			round(sd(sumx[,1]/10^6),2),
			sep=" ± "),
		paste(round(mean(sumx[,2])),
			round(sd(sumx[,2])),
			sep=" ± "),
		paste(round(mean(sumx[,3])),
			round(sd(sumx[,3])),
			sep=" ± "),	
		paste(round(mean(sumx[,5]),3),
			round(sd(sumx[,5]),3),
			sep=" ± "),
		paste(round(mean(sumx[,6]),3),
			round(sd(sumx[,6]),3),
			sep=" ± ")	
		))			
}))

SUMgroup2<-t(sapply(sort(unique(GroupLIN)),function(x){
	sumx<-matrix(as.numeric(as.matrix(
		as.data.frame(subset(SUMgenome, 
		GroupLIN==x))[c(9,10)])),
		ncol=2)
	nSPE<-length(unique(subset(SpeLIN, GroupLIN==x)))
		
	BIOMEx<-subset(BIOME,GroupLIN==x)	
	return(c(
		length(sumx[,1]), nSPE,	
		paste(round(mean(sumx[,1]),0),
			round(sd(sumx[,1]),0),
			sep=" ± "),
		paste(round(mean(sumx[,2]),0),
			round(sd(sumx[,2]),0),
			sep=" ± ")	
		))			
}))

SUMgroup<-cbind(SUMgroup[,c(1:6)], SUMgroup2[,c(3:4)],SUMgroup[,c(7)])

colnames(SUMgroup)<-c("Genomes","Species","Size(Mb)","Annotations","Unique-Annotations","Est-Copy-Number","Repeat-Elements","Mobile-Elements","GC-content")

setwd(paste(p0,"Gene-content",sep="/"))

write.table(SUMgroup,"Summary-genome-per-group.txt")

#Tukey test among groups for different genome statistics

pdf("Stat-genome-per-group.pdf")

par(bty="n",mfrow=c(2,3),mar=c(5,4,1,1))

TukeyGroup<-sapply(c(5,7,9,10,11,12,13,15),function(x){
	varx<-as.numeric(unlist(SUMgenome[,x]))
	Nx<-colnames(SUMgenome)[x]
	TUKEYx<-TukeyHSD(aov(
		glm(varx ~ GroupLIN)))$GroupLIN
	
	plot(as.factor(GroupLIN), varx,
		cex.axis=0.8,las=2,xlab="",ylab=Nx,
		col= COLgroup[order(COLgroup[,1]),4])
	return(TUKEYx[,4])
})
colnames(TukeyGroup)<-colnames(SUMgenome)[c(5,7,9,10,11,12,13,15)]

dev.off()

write.table(TukeyGroup,"Padj-Tukey-genome-per-group.txt")


library(vegan)

COLg<-sapply(GroupLIN,function(x){
	return(subset(COLgroup[,4], COLgroup[,1]==x))
})

DATA<-cbind(as.numeric(unlist(
	SUMgenome[,5]))/10^6,
	as.numeric(unlist(SUMgenome[,12])),
	as.numeric(unlist(SUMgenome[,13])),
	as.numeric(unlist(SUMgenome[,15]))*100)
DATA<-subset(DATA, GroupLIN!="Enterovirga")	
COLdata<-subset(COLg, GroupLIN!="Enterovirga")
BIOMEdata<-subset(BIOME,GroupLIN!="Enterovirga")
GROUPdata<-subset(GroupLIN, GroupLIN!="Enterovirga")
DATA<-subset(DATA, GROUPdata!="Microvirga")	
COLdata<-subset(COLdata, GROUPdata!="Microvirga")
BIOMEdata<-subset(BIOMEdata, GROUPdata!="Microvirga")
GROUPdata <-subset(GROUPdata, GROUPdata!="Microvirga")


labx=c("Genome size (Mb)","Unique annotations","Copy number per annotation","Annotation GC %")

pdf("Genome-Characteristics.pdf",width=8,height=6)

par(mar=c(4,4,1,1),bty="n",mfrow=c(2,3))

sapply(c(1:3),function(x){

X=x

sapply(c((X+1):4),function(x){

Y=x
plot(1,1,log="",las=1,cex.axis=0.8,
	xlab= labx[X],
	ylab= labx[Y],
	ylim=c(min(DATA[,Y]),max(DATA[,Y])),
	xlim=c(min(DATA[,X]),max(DATA[,X])))

ordispider(
	DATA[,c(X,Y)],
	groups= GROUPdata,
	#border=NA,
	#border= COLgroup[order(COLgroup[,1]),4][1:4],
	col= COLgroup[order(COLgroup[,1]),3][1:4],
	label=F,lwd=0.5)


ordiellipse(
	DATA[,c(X,Y)],
	groups= GROUPdata,
	border=NA,
	#border= COLgroup[order(COLgroup[,1]),4][1:4],
	col= COLgroup[order(COLgroup[,1]),2][1:4],
	label=F,alpha=30,draw="polygon",
	lwd=1,kind="sd",conf=0.5)
	
points(DATA[,c(X,Y)],
	col= COLdata,bg="white",
	pch=21,cex=0.6)
	
sapply(unique(GROUPdata),function(x){
	text(mean(subset(DATA[,X], GROUPdata==x)),
		mean(subset(DATA[,Y], GROUPdata==x)),
		labels=x,cex=1.2,font=1,
		#col=subset(COLgroup[,4],COLgroup[,1]==x)
		col="black")
})	
	
	})
	})

dev.off()

####Step 5c2: Matrix of gene occurence per species aligned on the Astral species tree


p0<-"/Users/jean-baptisteleducq/Desktop/New-Methylo-Genomes/Alessa-et-al-2021-genomes/"

setwd(p0)

SPE <-read.table("Species-tree-order.txt",header=T)

SPE<-SPE[order(SPE[,1]),]

setwd(paste(p0,"Gene-content",sep="/"))

#Annotation abundance per genome
Fall<-as.matrix(read.table("Annotation-Abundance.txt",header=T))

#Remove Annotations present in average in more than 50 copies in genomes
Fall<-subset(Fall,apply(Fall,1,mean)<50)

#Annotations occurence
OCC<-length(Fall[1,])-sapply(c(1:length(Fall[,1])),function(x){
	return(length(subset(Fall[x,],Fall[x,]==0)))
})

#Average annotation occurence per species 

Fsep<-t(sapply(SPE[,2],function(x){
	SPEx=x
	stx<-unlist(strsplit(
		subset(SPE[,3],SPE[,2]==x),
		split='_'))	
	Fx<-apply(sapply(c(stx,stx),function(x){
		return(as.vector(subset(
			t(Fall),colnames(Fall)==x)))
	}),1,mean)
	return(Fx)
}))


#Color code for groups
COLgroup<-cbind(
	c("A","B","C","D",
	"Enterovirga","Microvirga"),
	c(rgb(1,0,0,0.1),
	rgb(0.8,0,1,0.1),
	rgb(0,0.8,0,0.1),
	rgb(0,0.5,1,0.1),
	rgb(0.7,0.7,0.7,0.1),
	rgb(0.8,0.8,0.8,0.1)),
	c(rgb(1,0,0,0.4),
	rgb(0.8,0,1,0.4),
	rgb(0,0.8,0,0.4),
	rgb(0,0.5,1,0.4),
	rgb(0.7,0.7,0.7,0.4),
	rgb(0.8,0.8,0.8,0.4)),
	c(rgb(1,0,0,1),
	rgb(0.8,0,1,1),
	rgb(0,0.8,0,1),
	rgb(0,0.5,1,1),
	rgb(0.7,0.7,0.7,1),
	rgb(0.8,0.8,0.8,1)))
	
COLgroup<-	COLgroup[c(1,4,2,3,5,6),]

###To match with the annotation tree: strain to species
#For heatmap: color codes for average abundances

COLback<-rgb(0,0,0,1)

colSPE<-cbind(COLback,
t(sapply(SPE[,2][length(SPE[,1]):1],function(x){
	GROUPx<-subset(SPE[,4], SPE[,2]==x)
	return(subset(COLgroup[,c(2:4)], 
		COLgroup[,1]== GROUPx))
})))


#Gene abundance per GROUP weighted by species

Fgroup<-sapply(COLgroup[,1],function(x){
	SPEx<-subset(SPE[,2],SPE[,4]==x)
	return(apply(sapply(SPEx,function(x){
		return(as.numeric(subset(
			as.matrix(Fsep),
			rownames(Fsep)==x)))
	}),1,mean))
})

FgroupO<-ifelse(round(Fgroup,1)<0.5,0,
	ifelse(round(Fgroup)==0,1,
	ifelse(round(Fgroup)==1,2,3)))

Fgroup<-apply(FgroupO,1,paste,collapse="")

Gs<-sapply(COLgroup[,1],function(x){
	return(unlist(strsplit(x,split=""))[1])
})

Fgroup2<-as.vector(sapply(Fgroup,function(x){
	Fx<-as.numeric(unlist(
		strsplit(x,split="")))
	Fx<-ifelse(Fx>=2,2,0)	
	return(paste(subset(Gs,Fx==2),
		collapse=""))	
}))

Fgroup3<-as.vector(sapply(Fgroup,function(x){
	Fx<-as.numeric(unlist(
		strsplit(x,split="")))[1:4]
	Fx<-ifelse(Fx>=2,2,0)	
	return(paste(subset(Gs[1:4],Fx==2),
		collapse=""))	
}))

z<-order(Fgroup,decreasing=T)

#Group average annotation abundance Y axis
Ygroup<-seq(4,24,length.out=length(FgroupO[1,]))

#Width of the plot
W=100

#Space for the tree
SpT=10

setwd(paste(p0,"Gene-content",sep="/"))

pdf("Gene-content-matrix.pdf",height=4,width=10)

par(mar=c(0.5,0,0.5,0.5),bty="n")

plot(-100,-100,
	xlim=c(-1, W),ylim=c(0,length(SPE[,1])+max(Ygroup)),
	xaxt="n",yaxt='n',xlab="",ylab="")

Lz<-length(z)

FsepO<-ifelse(round(Fsep,1)<0.5,0,
	ifelse(round(Fsep)==0,1,
	ifelse(round(Fsep)==1,2,3)))


sapply(c(Lz:1),function(x){
	print(x)
	X=x
	Fx<-as.vector(FsepO[,z[X]])
	#Fx<-ifelse(Fx==0,0,
	#	ifelse(Fx==1,1,2))
	ax= SpT +(W-SpT)*(X/Lz)
	bx= SpT +(W-SpT)*((X+1)/Lz)
	Fx<-Fx[length(Fx):1]
	#colx<-ifelse(Fx==0,"black",
	#	ifelse(Fx==1,"orange",
	#	"orange4"))
		
	sapply(c(1:length(Fx)),function(x){
		polygon(c(ax,ax,bx,bx),
			c(x,x-1,x-1,x),
			border=NA,
			col= colSPE [x,(Fx[x]+1)])
			#col= colx[x])
	})
	return("")
	
})

#Average annotation abundance per group

sapply(c(1:Lz),function(x){
	print(x)
	X=x
	Fx<-round(as.vector(FgroupO[z[X],]))
	#Fx<-ifelse(Fx==0,0,
	#	ifelse(Fx==1,1,2))
	ax= SpT +(W-SpT)*(X/Lz)
	bx= SpT +(W-SpT)*((X+1)/Lz)
	#Fx<-Fx[length(Fx):1]
	#colx<-ifelse(Fx==0,"black",
	#	ifelse(Fx==1,"orange",
	#	"orange4"))
	COLx<-cbind(COLgroup[,1],
		COLback, COLgroup[,c(2:4)])
		
	sapply(c(1:length(Fx)),function(x){
		Y<-Ygroup[x]+length(SPE[,1])
		Gx= COLx[x,1]
		colx<-COLx[x,(Fx[x]+2)]
		polygon(c(ax,ax,bx,bx),
			c(Y,Y-3,Y-3,Y),
			border=NA,
			col= colx)
			#col= colx[x])
	})
	return("")
	
})

sapply(c(1:length(COLgroup[,1])),function(x){
	Y<-Ygroup[x]+length(SPE[,1])
	Gx= paste(ifelse(
		COLgroup[,1]=="Enterovirga",
		"",ifelse(
		COLgroup[,1]=="Microvirga","",
		"Methylobacterium ")),
		COLgroup[,1],sep="")[x]
	colx<-COLgroup[x,4]
	text(SpT,Y-1.5, Gx,pos=2,
		adj=0,cex=0.6,font=3,col= colx,
		offset=1)
})

legend("topright",cex=0.45, horiz=F,
	legend= c("n=0.5","",
		"","","","",
		"n=1","","","",
		"","","n>2","",
		"","","",""),
	text.col= "black",box.col="black",pch=15,
	border=NA,col= as.vector(COLgroup[,c(2:4)]),
	title="Average annotation copy number",title.col="black",ncol=3,
	text.font=3,pt.cex=1,bg="white")

dev.off()


##Calculate pan genome size per group (genes presents regardless their abundance)

K =100  ##Number of replicates 


pdf(paste("Pan-Core-Genome-Size-Estimation-K=",K,".pdf",sep=""),height=4,width=7)

par(mar=c(4,4,1,1),bty="n",mfrow=c(1,2))

plot(100000000,100000000,
	xlim=c(1,max(summary(as.factor(SPE[,4])))),
	ylim=c(2400,6400),log="",
	xlab="Species per group",ylab="Pan genes per group",cex.axis=0.8,las=1,
	main="")
	
points(c(15,15),c(0,6500),lty=2,type="l")	
	
legend("bottomright",cex=0.8, horiz=F,
	legend= COLgroup[c(1:4,6),1],
	text.col= "black",box.col=NA,
	border=NA,fill= as.vector(COLgroup[c(1:4,6),4]),
	title="Groups",title.col="black",ncol=1,
	text.font=3,pt.cex=1,bg="white")	

PAN<-sapply(c(6,1:4),function(x){
	Gx=COLgroup[x,1]
	SPEx<-subset(SPE[,2],SPE[,4]== Gx)
	Fx<-sapply(SPEx,function(x){
		return(as.numeric(subset(
			as.matrix(Fsep),
			rownames(Fsep)==x)))
	})
	Fx<-ifelse(Fx==0,0,1)
	Fx<-subset(Fx,apply(Fx,2,sum)!=0)
	Pk<-sapply(c(1:length(SPEx)),function(x){
		print(x)
		nx=x
		Px<-sapply(c(1:K),function(x){
			z<-sample(c(1:length(SPEx)),size= nx)
			fz<-Fx[,c(z,z)]
			return(length(subset(fz[,1],
				apply(fz,1,sum)!=0)))
		})
		return(c(mean(Px),sd(Px)))
	})
	polygon(c(c(1:length(SPEx)),c(length(SPEx):1)),
		c(Pk[1,]-Pk[2,],
		Pk[1,c(length(SPEx):1)]+
		Pk[2,c(length(SPEx):1)]),
		border=NA,col=COLgroup[x,2])
	points(c(1:length(SPEx)),
		Pk[1,],type="l",lwd=2,
		col=COLgroup[x,4])
	return(t(cbind(x,
		c(1:length(SPEx)),t(Pk))))
})

plot(100000000,100000000,
	xlim=c(1,max(summary(as.factor(SPE[,4])))),
	ylim=c(400,2500),log="",
	xlab="Species per group",ylab="Core genes per group",cex.axis=0.8,las=1,
	main="")

points(c(15,15),c(0,6500),lty=2,type="l")	

CORE<-sapply(c(6,1:4),function(x){
	Gx=COLgroup[x,1]
	SPEx<-subset(SPE[,2],SPE[,4]== Gx)
	Fx<-sapply(SPEx,function(x){
		return(as.numeric(subset(
			as.matrix(Fsep),
			rownames(Fsep)==x)))
	})
	Fx<-ifelse(round(Fx,1)==1,1,0)
	Pk<-sapply(c(1:length(SPEx)),function(x){
		print(x)
		nx=x
		Px<-sapply(c(1:K),function(x){
			z<-sample(c(1:length(SPEx)),size= nx)
			fz<-Fx[,c(z,z)]
			return(length(subset(fz[,1],
				apply(fz,1,sum)==(2*nx))))
		})
		return(c(mean(Px),sd(Px)))
	})
	polygon(c(c(1:length(SPEx)),c(length(SPEx):1)),
		c(Pk[1,]-Pk[2,],
		Pk[1,c(length(SPEx):1)]+
		Pk[2,c(length(SPEx):1)]),
		border=NA,col=COLgroup[x,2])
	points(c(1:length(SPEx)),
		Pk[1,],type="l",lwd=2,
		col=COLgroup[x,4])
	return(t(cbind(x,
		c(1:length(SPEx)),t(Pk))))
})

dev.off()

CORE<-matrix(unlist(CORE),ncol=4,byrow=T)
CORE<-cbind(COLgroup[CORE[,1],1], CORE)
colnames(CORE)<-c("Group","Index","Species","Mean","SD")

PAN<-matrix(unlist(PAN),ncol=4,byrow=T)
PAN <-cbind(COLgroup[PAN[,1],1], PAN)
colnames(PAN)<-c("Group","Index","Species","Mean","SD")

write.table(CORE,paste("Core-Genome-Size-Estimation-K=",K,".txt",sep=""),row.names=F)
write.table(PAN,paste("Pan-Genome-Size-Estimation-K=",K,".txt",sep=""),row.names=F)

#Venn diagram on pan genome assuming 15 species sampled

K=100
SS=15

#Categories
CAT<-sapply(c(1:4),function(x){
	Gx=COLgroup[x,1]
	SPEx<-subset(SPE[,2],
		SPE[,4]== Gx)
	Fx<-sapply(SPEx,function(x){
		return(as.numeric(subset(
			as.matrix(Fsep),
		rownames(Fsep)==x)))
	})
	Fx<-apply(ifelse(round(Fx,1)==1,1,0),1,sum)
	return(ifelse(Fx==length(SPEx),Gx,""))
})
CAT <-summary(as.factor(apply(CAT,
	1,paste,collapse="")))


COREss<-sapply(c(1:K),function(x){
	print(x)
	Kx=x
	Tx<-sapply(c(1:4),function(x){
		Gx=COLgroup[x,1]
		SPEx<-sample(subset(SPE[,2],
			SPE[,4]== Gx),size=SS)
		Fx<-sapply(SPEx,function(x){
			return(as.numeric(subset(
				as.matrix(Fsep),
			rownames(Fsep)==x)))
		})
		Fx<-apply(ifelse(round(Fx,1)==1,1,0),1,sum)
		return(ifelse(Fx==SS,Gx,""))
	})
	Tx<-summary(as.factor(apply(Tx,
		1,paste,collapse="")))	
	return(sapply(names(CAT),function(x){
		return(subset(Tx,names(Tx)==x))
	}))	
	
})

rownames(COREss)<-names(CAT)
COREss<-subset(COREss,rownames(COREss)!="")

mean(apply(COREss,2,sum))
sd(apply(COREss,2,sum))
round(apply(COREss,1,mean))

t(sapply(c(1:4),function(x){
	gx=COLgroup[x,1]
	cx<-apply(COREss [grep(gx,
		rownames(COREss)),],2,sum)
	return(c(gx,round(c(mean(cx),sd(cx)),0)))
}))


PANss<-sapply(c(1:K),function(x){
	print(x)
	Kx=x
	Tx<-sapply(c(1:4),function(x){
		Gx=COLgroup[x,1]
		SPEx<-sample(subset(SPE[,2],
			SPE[,4]== Gx),size=SS)
		Fx<-sapply(SPEx,function(x){
			return(as.numeric(subset(
				as.matrix(Fsep),
			rownames(Fsep)==x)))
		})
		Fx<-apply(ifelse(Fx==0,0,1),1,sum)
		return(ifelse(Fx==0,"",Gx))
	})
	Tx<-summary(as.factor(apply(Tx,
		1,paste,collapse="")))	
	return(sapply(names(CAT),function(x){
		return(subset(Tx,names(Tx)==x))
	}))	
	
})

rownames(PANss)<-names(CAT)
PANss <-subset(PANss,rownames(PANss)!="")

mean(apply(PANss,2,sum))
sd(apply(PANss,2,sum))
round(apply(PANss,1,mean))
round(apply(PANss,1,sd))

t(sapply(c(1:4),function(x){
	gx=COLgroup[x,1]
	cx<-apply(PANss [grep(gx,
		rownames(PANss)),],2,sum)
	return(c(gx,round(c(mean(cx),sd(cx)),0)))
}))

write.table(t(COREss),paste("Venn-Core-Genome-Size-Estimation-K=",K,"-ss=",SS,".txt",sep=""),row.names=F)
write.table(t(PANss),paste("Venn-Pan-Genome-Size-Estimation-K=",K,"-ss=",SS,".txt",sep=""),row.names=F)

#####Step 5c3: 
#Matrix distance in annotation occurence using bray distance

#install.packages("https://cran.r-project.org/bin/macosx/contrib/4.1/vegan_2.5-7.tgz")
#install.packages("https://cran.r-project.org/bin/macosx/contrib/4.1/permute_0.9-5.tgz")

library(vegan)

p0<-"/Users/jean-baptisteleducq/Desktop/New-Methylo-Genomes/Alessa-et-al-2021-genomes/"

setwd(paste(p0,"Gene-content",sep="/"))

Fall<-as.matrix(read.table("Annotation-Abundance.txt",header=T))

#Remove Annotations present in average in more than 50 copies in genomes
Fall<-subset(Fall,apply(Fall,1,mean)<50)

#Hellinger transformation to minimize the effect of false gene duplication due to genome incompleteness

FallH <-as.data.frame(decostand(Fall ,method="hellinger", MARGIN=2))

ANNdist<-sapply(c(1:length(Fall[1,])),function(x){
	print(x)
	F1<-as.vector(FallH[,x])
	sapply(c(1:length(Fall[1,])),function(x){
		F2<-as.vector(FallH[,x])
		return(vegdist(t(cbind(F1, F2)),method="bray"))
		#return(vegdist(t(ifelse(cbind(F1, F2)==0,
		#	0,1/cbind(F1, F2))),method="bray"))
	})
})

colnames(ANNdist)<-colnames(Fall)
rownames(ANNdist)<-colnames(Fall)

write.table(ANNdist,"Anno-Bray.txt",quote=F)

pdf("Anno-Bray.pdf",width=7,height=7)

heatmap(1-ANNdist,scale="none")

dev.off()

###Calculate average Bray index per species (display as in the same order as in ASTRAL species tree)

setwd(p0)

SpeTreeOrder<-read.table("Species-tree-order.txt",header=T)
SpeTreeOrder<-SpeTreeOrder[order(as.numeric(SpeTreeOrder[,1]),decreasing=T),]

ANNdistSPE<-sapply(SpeTreeOrder[,2],function(x){
	SPEx<-x
	STx<-unlist(strsplit(subset(SpeTreeOrder[,3], 
		SpeTreeOrder[,2]== SPEx),split="_"))
	ix<-sapply(STx,function(x){
		return(grep(paste(x,"-"),
			paste(colnames(ANNdist),"-")))
	})
	sapply(SpeTreeOrder[,2],function(x){
		SPEy<-x
		STy<-unlist(strsplit(subset(SpeTreeOrder[,3], 
			SpeTreeOrder[,2]== SPEy),split="_"))
		iy<-sapply(STy,function(x){
			return(grep(paste(x,"-"),
				paste(colnames(ANNdist),"-")))
		})
	
		ANNdistxy<-ANNdist[ix,iy]
		return(mean(ANNdistxy))
	})
})

#Color legend

COL<-cbind(
c(seq(1,0,length.out=31),seq(0,0.25,length.out=40),seq(0.25,0,length.out=30)),
c(seq(1,1,length.out=31),seq(1,0.25,length.out=40),seq(0.25,0,length.out=30)),
c(seq(1,1,length.out=70),seq(1,0,length.out=31))
)

COLa<-unique(sapply(
	c(1:length(COL[,1])),function(x){
	COLx<-COL[x,]
	return(rgb(COLx[1], COLx[2], COLx[3]))
}))

COLa<-COLa[length(COLa):1]


setwd(paste(p0,"Gene-content",sep="/"))

pdf("legend-AN.pdf",height=2.5,width=6)
par(bty="n")

plot(-100,-100,yaxt='n',xaxt="n",
	xlab="Pairwise gene abund. dissimilarity",ylab="",
	xlim=c(1,101),
	ylim=c(0,1),cex.lab=2)
axis(1,c(0,25,50,75,100),
	round(seq(max(ANNdistSPE),min(ANNdistSPE),
	length.out=5),2),cex.axis=1.5)	

sapply(c(1:length(COLa)),function(x){	
	polygon(c(x-1,x-1,x,x),
		c(0,1,1,0),
		border=NA,
		col= COLa[x])
})

dev.off()

pdf("Anno-Bray-Species-ASTRAL.pdf")

heatmap(1-ANNdistSPE,scale="none",Rowv=NA,Colv=NA,col=COLa,margin=c(5,5),cexRow=0.3,cexCol=0.3)

dev.off()

#Step 5d: Calculate average bray dissimilarity index among groups and within group among species

#Alessa genomes annotation file location
p0<-"/Users/jean-baptisteleducq/Desktop/New-Methylo-Genomes/Alessa-et-al-2021-genomes/"

setwd(p0)

#Strain characteristics
ID<-read.table("Table-Correspondance-ID.txt",header=T)
SPE<-read.table("Species-tree-order.txt",header=T)
SPE<-SPE[order(SPE[,1]),]

#Synteny index
setwd(paste(p0,"Gene-content",sep="/"))
ANNO<-read.table("Anno-Bray.txt",header=T)

#Type of comparisons

ANNOcat<-unlist(sapply(c(1:(length(SI[,1])-1)),function(x){
	X=x
	stx<-row.names(ANNO)[x]
	spx<-SPE[grep(paste(stx,"_",sep=""),
			paste(SPE[,3],"_",sep="")),2]
	gx<-SPE[grep(paste(stx,"_",sep=""),
			paste(SPE[,3],"_",sep="")),4]
	sapply(c((X+1):(length(ANNO[,1]))),function(x){
		Y=x
		sty<-row.names(ANNO)[x]
		spy<-SPE[grep(paste(sty,"_",sep=""),
			paste(SPE[,3],"_",sep="")),2]
		gy<-SPE[grep(paste(sty,"_",sep=""),
			paste(SPE[,3],"_",sep="")),4]
		return(paste(
			ifelse(spx==spy,paste("spe",gx,sep="-"),
			ifelse(gx==gy,gx,
			paste(sort(c(gx,gy)),collapse="-"))),
			unlist(ANNO[X,Y]),sep="="))
	})		
}))

ANNOcat <-matrix(unlist(strsplit(ANNOcat,
	split="=")),ncol=2,byrow=T)
	
ANNOgs<-t(sapply(sort(unique(SPE[,4])),function(x){
	spex<-as.numeric(subset(ANNOcat[,2],
		ANNOcat[,1]==paste("spe-",x,sep="")))
	gx<-as.numeric(subset(ANNOcat[,2],
		ANNOcat[,1]==x))
	return(c(
		paste(round(c(mean(spex),
			sd(spex)),2),collapse=" ± "),
		paste(round(c(mean(gx),
			sd(gx)),2),collapse=" ± ")))
}))

colnames(ANNOgs)<-c("within_sp","among_sp")

CATg<-sapply(sort(unique(SPE[,4])),function(x){
	catx<-paste(x,sort(unique(SPE[,4])),sep="-")
	return(sapply(catx,function(x){
		gx<-as.numeric(subset(
			ANNOcat[,2], ANNOcat[,1]==x))
		return(paste(round(c(mean(gx),
			sd(gx)),2),collapse=" ± "))
	}))
})

rownames(CATg)<-rownames(ANNOgs)

ANNOgs <-cbind(ANNOgs,CATg)

ANNOgs <-ifelse(ANNOgs =="NaN ± NA","", ANNOgs)

setwd(paste(p0,"Gene-content",sep="/"))

write.table(ANNOgs,"Summary-ANNO-per-Group.txt")


#######STEP 6a: SYNTENY


#Alessa genomes annotation file location
p0<-"/Users/jean-baptisteleducq/Desktop/New-Methylo-Genomes/Alessa-et-al-2021-genomes/"

pA<-paste(p0,"Annotations-16-11-2021",sep="/")

#Previously annotated genomes annotation file location
#MyRAST annotations (see script Analysis-Gene-annotation.R)
pS="/Users/jean-baptisteleducq/Desktop/myRAST-out"

#Genome and gene description
setwd(paste(p0,"Gene-content",sep="/"))

SUMgenome<-read.table("Annotation-Summary-per-genome.txt",header=T)

ST<-SUMgenome[,2]

SUMgene <-read.table("Annotation-Summary-per-gene.txt",header=T)

#Only keep core genes
SUMgene<-subset(SUMgene, is.na(SUMgene[,3])==F)

#Core gene simplified names (used in RAxML trees)
COREn<-SUMgene[,3]
#Core gene complete name
Fn<-SUMgene[,1]

####Retrieve core gene position in each genome

dataset<-as.numeric(matrix(unlist(
	strsplit(SUMgenome[,3],split="[.]")),
	ncol=3,byrow=T)[,2])

#Coordinates of core genes in each genome
COORD<-sapply(ST,function(x){
	print(x)
	GFFx<-subset(SUMgenome[,3], SUMgenome[,2]==x)
	datasetx<-subset(dataset, SUMgenome[,2]==x)
	setwd(ifelse(datasetx > 800000 , pA, pS))
	ANNx<-read.delim(GFFx)	
	return(sapply(c(1:length(COREn)),
		function(x){
		ANNOx<-Fn[x]
		return(paste(c(x,unlist(
			subset(ANNx,ANNx[,8]==
			ANNOx)[1,c(1,5,6)]),
			NA,NA,NA)[1:4],collapse="="))	
	}))
})

###Number of scaffolds with core gene per genome
NSCA<-sapply(c(1:length(ST)),function(x){
	STx<-ST[x]
	print(STx)
	#core gene coordinates
	coordx<-matrix(unlist(strsplit(
		as.vector(COORD[,x]),split="=")),
		ncol=4,byrow=T)
	return(length(unique(coordx[,2])))	
})

#For each genome reconstruct scaffolds with core gene order, writen in a separate file. Put orphan gene and missing gene in separate scaffolds

pY<-paste(p0,"Synteny","coord",sep="/")

setwd(pY)
sapply(c(1:length(ST)),function(x){
	STx<-ST[x]
	print(STx)
	#core gene coordinates
	coordx<-matrix(unlist(strsplit(
		as.vector(COORD[,x]),split="=")),
		ncol=4,byrow=T)
	coordx<-subset(coordx,coordx[,2]!='NA')	
	#Linkage groups (genes in the same order as in the scaffold)
	Gx<-as.vector(sapply(unique(coordx[,2]),
		function(x){
		coord<-subset(coordx, coordx[,2]==x)
		return(list(
			as.numeric(coord[order(apply(cbind(
			as.numeric(coord[,3]),
			as.numeric(coord[,4])),
			1,min,na.rm=T)),1])))
	}))
	#Orphan and missing genes (considered in separate scaffolds)
	Gna<-setdiff(c(1:length(COREn)),
		unlist(Gx))
	Gna<-sapply(Gna,function(x){
		return(list(x))
	})
	Gx<-c(Gx,Gna)	
	Gx<-sapply(c(1:length(Gx)),function(x){
		return(paste(Gx[[x]],collapse="_"))
	})
	write.table(Gx,paste("strain=",
		x,".txt",sep=""),
		row.names=F,col.names=F)
})


####Reconstruct core gene organisation for incomplete genomes

#Iteratively align draft genomes on scaffold genomes
#Draft genome indexes
Di<-subset(c(1:length(ST)),NSCA>1)

#Complete genome indexes
Ci<-subset(c(1:length(ST)),NSCA==1)

#Sort draft genome indexs from the lowest to the highest NSCA
Di<-Di[order(NSCA[Di],decreasing =F)]

###FUNCTIONS
num.strsplit.u<-function(x){return(as.numeric(unlist(strsplit(x,split= "_"))))}

#align each draft genome scaffold on each complete genome and rearanged genomes with highest completeness

sapply(c(1:length(Di)),function(x){
	Dix=Di[x]
	print(x)
	STx<-ST[Dix]
	
	#Index of genomes with lowest NSCA
 	Cix<-c(Ci,Di[subset(c(1:length(Di)),
 		c(1:length(Di))<x)])

	#Links in genomes with lowest NSCA
	LINKc<-sapply(Cix,function(x){
 		STy<-ST[x]
 		z<-ifelse(length(intersect(x,Di))==1,
 			"REA-","")	
 		GKy<-as.vector(read.table(
 			paste(z,"strain=",x,".txt",sep=""))[,1])
 		#gene initial order
 		ORDi<-num.strsplit.u(GKy)
 		#Links in this genome
 		ORDy<-c(sapply(c(2:length(ORDi)),
 			function(x){
 			return(paste(sort(
 			c(ORDi[x-1], ORDi[x])),
 			collapse="_"))
 		}),
 		paste(sort(c(ORDi[1], 
 			ORDi[length(ORDi)])),	
 			collapse="_"))		
 		return(ORDy)	
	})

 	#Scaffolds in the draft genome
 	Gx<-as.vector(read.table(paste(
 		"strain=", Dix,".txt",sep=""))[,1])
 	
 	#Draft genome scafffold index (first column) and gene name according to order in the scaffold (second column)	
 	Sx<-c(1:length(Gx))
 	Sx<-matrix(unlist(sapply(Sx,function(x){
 		return(t(cbind(x,num.strsplit.u(Gx[x]))))
 	})),ncol=2,byrow=T)
 	
 	#Start draft genome alignment on reference genomes and SI calculation after rearangement
	SIt<-sapply(Cix,function(x){
		STy<-ST[x]
		STY=x
 		z<-ifelse(length(intersect(STY,Di))==1,
 			"REA-","")	
 		#Reference genome: if already rearange, open the REA file, otherwise open the original file (for complete genomes only)
 		Gy<-as.vector(read.table(paste(
 			z,"strain=", STY,".txt",sep=""))[,1]) 	
 		Sy<-c(1:length(Gy))
 		Gy<-num.strsplit.u(Gy)
 		
 		#links in the reference genome
 		Lx<-as.vector(t(subset(t(LINKc), Cix== STY)))
 		
 		#For each gene in this reference genome,report the scaffold index of the draft genome
		Ixy<-unique(sapply(Gy,function(x){
			return(subset(Sx[,1],Sx[,2]==x))
		}))
		
		#Sort draft genome scaffolds according to this index
		Gxy<-Gx[Ixy]
		
		#For each draft scaffold, calculate a score with neighor scaffolds by looking at edges and adjust the scaffold orientation accordingly
		
		Gxy<-sapply(c(1:length(Gxy)),function(x){
			gx<-num.strsplit.u(Gxy[x])
			w=ifelse((x-1)<1,length(Gxy),x-1)
			y=ifelse((x+1)>length(Gxy),1,x+1)
			gw<-num.strsplit.u(Gxy[w])
			gy<-num.strsplit.u(Gxy[y])
			#Edges of this scaffold 
			x1<-gx[1]
			x2<-gx[length(gx)]
			#Edges of neigbor scaffolds
			y1<-gy[1]
			y2<-gy[length(gy)]			
			w1<-gw[1]
			w2<-gw[length(gw)]
			#links if the scaffold do not need to be inverted and occurence in the complete genome
			l1<-paste(c(x1,x1,x2,x2),
				c(w1,w2,y1,y2),sep="_")
			l1<-sum(sapply(l1,function(x){
				lx<-paste(sort(num.strsplit.u(x)),
					collapse="_")
				return(length(intersect(lx,Lx)))
			}))
			#links if the scaffold needs to be inverted and occurence in the complete genome
			l2<-paste(c(x1,x1,x2,x2),
				c(y1,y2,w1,w2),sep="_")
			l2<-sum(sapply(l2,function(x){
				lx<-paste(sort(num.strsplit.u(x)),
					collapse="_")
				return(length(intersect(lx,Lx)))
			}))
			#If sum(l2)>sum(l1), invert the scaffold
			return(paste(gx[order(c(1:length(gx)),
				decreasing=ifelse(l2>l1,T,F))],
				collapse="_"))			
		})
		#Convert scaffold in links
		#gene order
 		ORDi<-num.strsplit.u(Gxy)
 		#Links in this genome
 		ORDy<-c(sapply(c(2:length(ORDi)),
 			function(x){
 			return(paste(sort(
 			c(ORDi[x-1], ORDi[x])),
 			collapse="_"))
 		}),
 		paste(sort(c(ORDi[1], 
 			ORDi[length(ORDi)])),	
 			collapse="_"))		
 		#Calculate synteny with complete genomes
 		SIxy<-sapply(Cix,function(x){	
 			#links in the complete genome
 			Lx<-as.vector(t(subset(t(LINKc), Cix==x)))
 			return(length(intersect(ORDy,Lx)))
 		})/length(COREn)
 		
 		write.table(Gxy,paste("TEMP-", STY,sep=""),
 			row.names=F,col.names=F,quote=F)
 		
 		return(c(SIxy))
		
	})
	
	#Retrieve the rearangement with the highest average synteny index
	Cib<-subset(Cix,apply(SIt,2,max)==
		max(apply(SIt,2,max)))[1]
		
	Gnx<-as.vector(read.table(
		paste("TEMP-", Cib,sep=""))[,1])	
	write.table(Gnx,
		paste("REA-strain=", Dix,".txt",sep=""),
		row.names=F,col.names=F)
})


#Links in complete and rearanged draft genomes
LINKr<-sapply(c(1:length(ST)),function(x){
 	STy<-ST[x]
 	z<-ifelse(length(intersect(x,Di))==1,"REA-","")
 	
 	GKy<-as.vector(read.table(
 		paste(z,"strain=",x,".txt",sep=""))[,1])
 	#gene initial order
 	ORDi<-num.strsplit.u(GKy)
 	#Links in this genome
 	ORDy<-c(sapply(c(2:length(ORDi)),
 		function(x){
 		return(paste(sort(
 		c(ORDi[x-1], ORDi[x])),
 		collapse="_"))
 	}),
 	paste(sort(c(ORDi[1], 
 		ORDi[length(ORDi)])),	
 		collapse="_"))		
 	return(ORDy)	
})

colnames(LINKr)<-ST

setwd(paste(p0,"Synteny",sep="/"))

write.table(LINKr,"Rearanged-core-genomes.txt",row.names=F)

GENdesc<-cbind(c(1:length(COREn)), COREn,Fn)
colnames(GENdesc)<-c("Code-link","Short-name","Annotation")

write.table(GENdesc,"Corresp-core-gene-names.txt",row.names=F)

###Calculate a synteny index among all genomes after rearangement of draft genomes

SI<-sapply(c(1:length(ST)),function(x){
	X=x
	return(sapply(c(1:length(ST)),function(x){
		return(length(intersect(LINKr[,X],LINKr[,x])))
	}))
})/length(COREn)

colnames(SI)<-ST
rownames(SI)<-ST


write.table(SI,"Synteny-Index.txt",quote=F)

pdf("Synteny.pdf",width=7,height=7)

heatmap(1-SI,scale="none")

dev.off()

###Calculate average Synteny index per species (display as in the same order as in ASTRAL species tree)

setwd(p0)

SpeTreeOrder<-read.table("Species-tree-order.txt",header=T)
SpeTreeOrder<-SpeTreeOrder[order(as.numeric(SpeTreeOrder[,1]),decreasing=T),]

SyntenySPE<-sapply(SpeTreeOrder[,2],function(x){
	SPEx<-x
	STx<-unlist(strsplit(subset(SpeTreeOrder[,3], 
		SpeTreeOrder[,2]== SPEx),split="_"))
	ix<-sapply(STx,function(x){
		return(grep(paste(x,"-"),
			paste(colnames(SI),"-")))
	})
	sapply(SpeTreeOrder[,2],function(x){
		SPEy<-x
		STy<-unlist(strsplit(subset(SpeTreeOrder[,3], 
			SpeTreeOrder[,2]== SPEy),split="_"))
		iy<-sapply(STy,function(x){
			return(grep(paste(x,"-"),
				paste(colnames(SI),"-")))
		})
	
		SIxy<-SI[ix,iy]
		return(mean(SIxy))
	})
})

#Color legend


COL<-cbind(
seq(1,1,length.out=102),
c(seq(0,1,length.out=51),seq(1,1,length.out=51)),
c(seq(0,0,length.out=51),seq(0,1,length.out=51))
)

COLg<-unique(sapply(
	c(1:length(COL[,1])),function(x){
	COLx<-COL[x,]
	return(rgb(COLx[1], COLx[2], COLx[3]))
}))

setwd(paste(p0,"Synteny",sep="/"))

pdf("legend-SI.pdf",height=2.5,width=6)
par(bty="n")

plot(-100,-100,yaxt='n',xaxt="n",
	xlab="Synteny index",ylab="",
	xlim=c(1,101),
	ylim=c(0,1),cex.lab=2)
axis(1,c(0,25,50,75,100),
	round(seq(min(SyntenySPE),max(SyntenySPE),
	length.out=5),2),cex.axis=1.5)	

sapply(c(1:length(COLg)),function(x){	
	polygon(c(x-1,x-1,x,x),
		c(0,1,1,0),
		border=NA,
		col= COLg[x])
})

dev.off()


pdf("Synteny-Species-ASTRAL.pdf")

heatmap(SyntenySPE,scale="none",Rowv=NA,Colv=NA,col= COLg,margin=c(5,5),cexRow=0.3,cexCol=0.3)

dev.off()

####


#Convert Synteny in 01 matrix and then in fasta file for phylogenetic analyses (1: link present; 0: link absent)

library(seqinr)

#Unique LINKs
LINKu<-unique(as.vector(LINKr))

SEQsi<-sapply(c(1:length(ST)),function(x){
	print(x)
	Lx<-as.vector(LINKr[,x])
	return(
	as.vector(sapply(LINKu,function(x){
		return(length(subset(Lx, Lx==x)))
	})))
})

#Remove annotations present in all genomes (unimformative)
SEQsi <-subset(SEQsi,apply(SEQsi,1,sum)!=length(ST))

SEQsi <-sapply(c(1:length(ST)),function(x){
	Fx<-as.vector(SEQsi[,x])
	return(list(Fx))
})

setwd(paste(p0,"Synteny",sep="/"))

write.fasta(SEQsi,names=ST,"Synteny.fas")

#### STEP 6b: Built consensus core genome architecture for species, main clades and groups using the Astral species tree as a guide, produce synteny maps for use in cytoscape (gene 1, gene2, weigth of the link)

p0<-"/Users/jean-baptisteleducq/Desktop/New-Methylo-Genomes/Alessa-et-al-2021-genomes/"

setwd(p0)

SPE <-read.table("Species-tree-order.txt",header=T)
SPE<-SPE[order(SPE[,1]),]
LIN <-read.table("Lineage-tree-order.txt",header=T)
ID<-read.table("Table-Correspondance-ID.txt",header=T)

setwd(paste(p0,"Gene-content",sep="/"))

SUMgenome<-read.table("Annotation-Summary-per-genome.txt",header=T)

setwd(paste(p0,"Synteny",sep="/"))

LINKr <-read.table("Rearanged-core-genomes.txt",header=T)
GENdesc<-read.table("Corresp-core-gene-names.txt",header=T)

### Unique links
LINKu<-sort(unique(as.vector(as.matrix(LINKr))))

setwd(paste(p0,"Synteny","Synteny-Maps",sep="/"))

###Average link frequency per species
LINKspe<-sapply(c(1:length(SPE[,1])),function(x){
	SPEx<-SPE[x,2]
	print(SPEx)
	STx<-unlist(strsplit(SPE[x,3],split="_"))
	ix<-sapply(STx,function(x){
		return(grep(paste(x,"_",sep=""),
			paste(colnames(LINKr),"_",sep="")))
	})
	#links relative frequency in this species
	LINKx<-summary(as.factor(
		as.vector(as.matrix(LINKr[,ix]))),
		maxsum=length(LINKu))/length(ix)
	#Produce a synteny consensus map for this species for network display in cytoscape
	LINKm<-cbind(matrix(unlist(strsplit(names(LINKx),split="_")),ncol=2,byrow=T), LINKx)
	colnames(LINKm)<-c("Gene1","Gene2","F")
	write.table(LINKm,paste("Consensus-Map-",SPEx,".txt",sep=""),quote=F,row.names=F)
	#Return link abundance
	return(sapply(LINKu,function(x){
		return(c(subset(LINKx,
			names(LINKx)==x),0)[1])
	}))	
})

colnames(LINKspe)<-SPE[,2]

write.table(LINKspe,"Consensus-synteny-per-species.txt")

#Using the Astral species tree as a guide, reconstructe the consensus synteny map for each node of the tree 

setwd(paste(p0,"For-RAxML-gene-trees",sep="/"))

ASTRALspe<-read.tree("Astral-species-tree.tre") 

#Convert the ASTRAL species tree in a bifurcating matrix (manually built in excell). In this file, nodes are numbered from the tip to the root, assuming a tree rooted on Microvirga + Enterovirga

setwd(p0)

MAPtree<-read.table("Species-tree-MAP-ASTRAL.txt",header=T)

###For each node (increasing order), build a consensus synteny map, using species synteny maps as origins

MAP<-MAPtree[,5:length(MAPtree[1,])]
MAP<-cbind(0,MAP)

#Nodes
ND<-sort(as.numeric(setdiff(
	unique(as.matrix(as.vector(MAP))),0)))

setwd(paste(p0,"Synteny","Synteny-Maps",sep="/"))

t(sapply(ND,function(x){
	NDx=x
	#Identify species in this node
	SPEx<-subset(MAPtree[,2],
		apply(ifelse(MAP==NDx,1,0),1,sum)!=0)
	#identify daugther nodes
	MAPx<-(t(subset(
		MAPtree[,c(5:length(MAPtree[1,]))],
		apply(ifelse(MAP==NDx,1,0),1,sum)!=0)))
	MAPx<-ifelse(MAPx>= NDx,0, MAPx)
	MAPx<-apply(MAPx,2,max)
	MAPx<-unique(ifelse(MAPx==0,SPEx,MAPx))
	#return(cbind(x,length(MAPx)))#control line to check that only 2 dauchter node are found
	#retrieve the two synteny maps
	MAP2<-as.vector(unlist(sapply(MAPx,function(x){
		mapx<-read.table(paste(
			"Consensus-Map-",x,".txt",sep=""),
			header=T)
		return(list(paste(paste(mapx[,1], mapx[,2],sep="_"),mapx[,3],sep="=")))
	})))
	MAP2<-matrix(unlist(strsplit(MAP2,split="=")),
		ncol=2,byrow=T)
		
	MAP2 <-sapply(unique(MAP2[,1]),function(x){
		return(sum(subset(as.numeric(MAP2[,2]), MAP2[,1]==x)))
	})/2
	print(c(x,sum(MAP2)))	
	MAP2 <-cbind(matrix(unlist(strsplit(names(MAP2),split="_")),ncol=2,byrow=T), round(MAP2,4))
	colnames(MAP2)<-c("Gene1","Gene2","F")
	write.table(MAP2,paste("Consensus-Map-", NDx,".txt",sep=""),quote=F,row.names=F)
	
}))

#Identify the genome with core genome architecture the closer from the consensus, use it as a reference and map synteny conservation of other genomes on it, ordered according to the synteny BINCAT tree (sup figure) 

#Methylobacterium consensus synteny map: 

METHmap<-read.table("Consensus-Map-120.txt",header=T)

#Links  
Lmeth<-paste(METHmap[,1], METHmap[,2],sep="_")

#Retrieve abundance per genome
AB<-sapply(c(1:length(LINKr[1,])),function(x){
	print(x)
	Lx<-LINKr[,x]
	return(sum(unlist(sapply(Lx,function(x){
		subset(METHmap[,3], Lmeth==x)
	}))))
})

#Strain with the core genome architecture the closestfrom methylobacterium

iref<-subset(c(1:length(AB)),AB==max(AB))[1]
STref<-colnames(LINKr)[iref]
LINKref<-LINKr[, iref]

#Colors codes for genes, ordered according to the reference genome	
	
setwd(p0)

COLg<-as.matrix(read.table("COL.txt"))
	
COLg<-unique(as.vector(sapply(c(2:length(COLg[,1])),function(x){
	col1<-as.vector(COLg[(x-1),])
	col2<-as.vector(COLg[(x),])
	rx<-seq(col1[1], col2[1],length.out=101)
	gx<-seq(col1[2], col2[2],length.out=101)
	bx<-seq(col1[3], col2[3],length.out=101)
	return(rgb(rx,gx,bx))
})))

COLg<-COLg[round(seq(1,length(COLg),length.out=length(LINKref)))]

#Genes as ordered in the reference genome
GENref<-
	unlist(sapply(c(2:(length(LINKref))),function(x){
		Lx<-as.numeric(unlist(strsplit(
			LINKref[x],split="_")))
		Ix<-as.numeric(unlist(strsplit(
			LINKref[(x-1)],split="_")))
		return(setdiff(Lx,Ix))	
	}))
GENref<-c(GENref,setdiff(
	as.numeric(unlist(strsplit(LINKref,split="_"))),
	GENref))
	
GENref<-cbind(GENref, COLg)
colnames(GENref)<-c("Gene","Color")

setwd(paste(p0,"Synteny","Synteny-Maps",sep="/"))	
	
write.table(GENref,paste("Color-Code-Synteny-Reference=", STref,".txt",sep=""),quote=F,row.names=F)		

##Plot Synteny conservation along the reference genome (species) inthe same order as in the ASTRAL species 

Lspe<-rownames(LINKspe)
Lspe<-as.vector(sapply(Lspe,function(x){
	return(unique(unlist(strsplit(x,split="[.]"))))
}))

#Average link conservation per species (only for links present in the reference genome)
CONspe<-sapply(c(1:length(SPE[,1])),function(x){
	#Average link frequency in this species
	Lx<-LINKspe[,x]
	#Retrieve only frequency for links shared with the reference genome
	Ls<-as.vector(sapply(LINKref,function(x){
		return(c(subset(Lx, Lspe ==x),0)[1])
	}))
	return(as.vector(Ls))
})

#Limits of the average SI plot / GROUP / LINK
S1=(-20)
S2=(-5)

#Limits of the average SI plot / SPECIES
L1= length(LINKref)+5
L2= length(LINKref)+20

#Color codes for groups 

COLgroup<-cbind(
	c("A","B","C","D",
	"Enterovirga","Microvirga"),
	c(rgb(1,0,0,1),
	rgb(0.8,0,1,1),
	rgb(0,0.8,0,1),
	rgb(0,0.5,1,1),
	rgb(0.7,0.7,0.7,1),
	rgb(0.8,0.8,0.8,1)),
	c(rgb(1,0,0,0.2),
	rgb(0.8,0,1,0.2),
	rgb(0,0.8,0,0.2),
	rgb(0,0.5,1,0.2),
	rgb(0.7,0.7,0.7,0.2),
	rgb(0.8,0.8,0.8,0.2)))
	
COLgroup2<-sapply(SPE[,4],function(x){
	return(subset(COLgroup[,2], COLgroup[,1]==x))
})	


setwd(paste(p0,"Synteny",sep="/"))

pdf("Synteny-HM.pdf",height=3,width=8)

par(mar=c(1,1,1,1),bty="n")

plot(-100,-100,xlim=c(-5, L2),
	ylim=c(length(SPE[,1]),S1),
	xlab="",ylab="",
	xaxt="n",yaxt="n")
	
polygon(c(0,length(LINKref),length(LINKref),0),
	c(0,0,length(SPE[,1]),length(SPE[,1])),
	border=NA,col=rgb(0,0,0,1))
	
#polygon(c(0,length(LINKref),length(LINKref),0),
#	c(S2,S2,S1,S1),
#	border=NA,col=rgb(0,0,0,1))	

#Plot average conservation of each link using average color code for the two linked genes, shaded by frequency in the species

segments(-1,S2,-1,S1,lwd=0.5)
segments(-2,S2,-1,S2,lwd=0.5)
segments(-2,S1,-1,S1,lwd=0.5)
text(c(-2,-2),c(S1,S2),c(1,0),cex=0.5,pos=2,offset=0.1)

sapply(c(1:length(LINKref)),function(x){
	Lx=x
	#Genes in this link
	#Gx<-as.numeric(unlist(
	#	strsplit(LINKref[Lx],split='_')))
	#Color codes
	#Cx<-sapply(Gx,function(x){
	#	return(subset(GENref[,2], GENref[,1]==x))
	#})
	#COlor codes combined
	#Cx<-apply(col2rgb(Cx),1,mean)/255
	#Cx<-1-(Cx * mean(Cx)/sum(Cx))
	#Grey code (for the background)
	#Gx=(1-mean(Cx))/2
	#Gx<-rgb(Gx, Gx, Gx)
	#polygon(c(Lx,Lx,Lx-1,Lx-1),
	#	c(length(SPE[,1]),0,0,length(SPE[,1])),
	#	col=Gx,
	#	border=NA)

	#For each species
	sapply(c(1:length(SPE[,1])),function(x){
		#Group
		GRx<-SPE[x,4]
		#Groups color
		colx<-col2rgb(subset(COLgroup[,2], 
			COLgroup[,1]==GRx))[,1]/255
		#Frequency of this link in this species
		Fx= CONspe[Lx,x]
		#Color code
		#cx<-rgb(Cx[1],Cx[2],Cx[3],Fx)
		cx<-rgb(colx[1],colx[2],colx[3],Fx)
		polygon(c(Lx,Lx,Lx-1,Lx-1),
			c(x,x-1,x-1,x),
			col=cx,
			border=NA)
	})
	
	return("")
})

#Average SI compared to the reference per species (colored according to groups)

AVsi<-apply(CONspe,2,sum)/length(LINKref)
AVsi<-L1+(AVsi*(L2-L1))

segments(L1,-1,L2,-1,lwd=0.5)
segments(L1,-1,L1,-2,lwd=0.5)
segments(L2,-1,L2,-2,lwd=0.5)
text(c(L1,L2),c(-2,-2),c(0,1),cex=0.5,pos=3,offset=0.1)


sapply(c(1:length(AVsi)),function(x){
	polygon(c(L1,AVsi[x],AVsi[x],L1),
		c(x-1,x-1,x,x),
		col=COLgroup2[x],border=NA)
})



#Average SI compared to the reference per link, per group (colored according to groups)

Fcons<-sapply(c("A","D","B","C")[4:1],function(x){
	#Node
	MAPx<-subset(MAP,MAPtree[,4]==x)
	SDx<-apply(MAPx,2,sd)
	MAPx<-apply(MAPx,2,min)
	SDx <-subset(SDx, MAPx!=0)
	MAPx<-subset(MAPx, MAPx!=0)
	NDx<-min(subset(MAPx,SDx==0))
	setwd(paste(p0,"Synteny",
		"Synteny-Maps",sep="/"))
	GROUPmap<-read.table(paste("Consensus-Map-", 
		NDx,".txt",sep=""),header=T)
	#Frequency of links from the reference
	Fx<-as.vector(sapply(LINKref,function(x){
		return(c(subset(GROUPmap[,3],
			paste(GROUPmap[,1], 
			GROUPmap[,2],sep="_")==x),0)[1])
	}))
	COLx<-subset(COLgroup[,2],COLgroup[,1]==x)
	COL2x<-subset(COLgroup[,3],COLgroup[,1]==x)	
	polygon(c(length(LINKref),0:length(LINKref)),
		c(0,0,Fx)*(S1-S2)+S2,
		col=COL2x,border=COLx,lwd=0.25)
	return(Fx)
})

z<-setdiff(subset(LINKref,Fcons[,3]==apply(Fcons,1,max)),
c(subset(LINKref,Fcons[,1]==apply(Fcons,1,max)),
	subset(LINKref,Fcons[,2]==apply(Fcons,1,max)),
	subset(LINKref,Fcons[,4]==apply(Fcons,1,max))))

zl<-sapply(LINKref,function(x){
	return(length(subset(z,z==x)))
})

points(c(1:length(LINKref)),ifelse(zl==1,(S1-S2)+S2,-100000)-2,pch=19,cex=0.2,col="blue")



dev.off()



z<-summary(as.factor(unlist(strsplit(z,split="_"))),maxsum=384)

z<-subset(names(z),z==2)
#z<-names(z)

setwd(paste(p0,"Synteny",sep="/"))

SUMgen<-read.table("Corresp-core-gene-names.txt",header=T)

SUMgenz<-t(sapply(z,function(x){
	return(subset(SUMgen, SUMgen[,1]==x))
}))

Zx<-SPE[c(23:27),3]
Ax<-setdiff(unlist(strsplit(subset(SPE[,3],SPE[,4]=="A"),split="_")), Zx)
Dx<-unlist(strsplit(subset(SPE[,3],SPE[,4]=="D"),split="_"))


setwd(paste(p0,"For-RAxML-gene-trees",sep="/"))

ADz<-t(sapply(SUMgenz[,2],function(x){
	treex<-unroot(read.tree(paste(
		"RAxML_bipartitionsBranchLabels.",
		x,"_fa-#1000_mGTRGAMMA",sep="")))
	return(list(treex))
}))

class(ADz)<-"multiPhylo"

write.tree(ADz,"AD-sim-synteny.tree")

#Then in shell, compute lineage tree with Astral (-Xmx3000M to allow the use of more memory by java)

cd Desktop/
cd myRAST-out
cd Astral.5.7.8

java -Xmx3000M -jar astral.5.7.8.jar -i With-genomes-Alessa-nov-2021/AD-sim-synteny.tree -o With-genomes-Alessa-nov-2021/Astral-AD-tree.tre  

#Also in shell, compute species tree with Astral (-Xmx3000M to allow the use of more memory by java)

cd Desktop/
cd myRAST-out
cd Astral.5.7.8

java -Xmx3000M -jar astral.5.7.8.jar -i With-genomes-Alessa-nov-2021/AD-sim-synteny.tree -a With-genomes-Alessa-nov-2021/Mapping-ASTRAL.txt -o With-genomes-Alessa-nov-2021/Astral-AD-species-tree.tre 

#Synteny plot for 25 phylogenetically representative species 

setwd(p0)

LISTsynt<-as.vector(read.table("List-Species-Synteny-plot.txt")[c(6:18),1])

LISTsynt<-c("Meth_trifolii", LISTsynt)
LISTsynt<-setdiff(LISTsynt,"Meth_hispanicum")

#For each species, report group and LINKr for the strain with the most complete genome

LISTsynt<-t(sapply(LISTsynt,function(x){
	SPEx=x
	Ox<-subset(SPE[,1],SPE[,2]==x)
	Gx<-subset(SPE[,4],SPE[,2]==x)
	COLx<-as.vector(subset(COLgroup, COLgroup[,1]==Gx)[1,])
	STx<-unlist(strsplit(subset(SPE[,3],SPE[,2]==x),split="_"))
	NSCA<-sapply(STx,function(x){
		return(subset(SUMgenome[,4], SUMgenome[,2]==x))
	})
	STx<-subset(STx, NSCA ==min(NSCA))[1]
	STn<-subset(SUMgenome[,1], SUMgenome[,2]==STx)
	STn <-gsub("Mb.","Meth_", STn)
	STn <-gsub("Mr.","Meth_", STn)
	STn <-paste(
		setdiff(unlist(strsplit(STn,split="_")),
		unlist(strsplit(SPEx,split="_"))),
		collapse="_")
	
	return(c(SPEx,Ox, COLx,STx, STn,min(NSCA)))
}))

vx<-seq(0,1,length.out=length(LINKref))

#Color legend

#SI range
SIr<-
unlist(sapply(c(1:(length(LISTsynt[,1])-1)),function(x){
	X=x
	STx<-LISTsynt[x,6]
	LX<-LINKr[,c(subset(c(1:length(LINKr[1,])),
		colnames(LINKr)==STx))]
	sapply(c((X+1):length(LISTsynt[,1])),function(x){
		Y=x
		STy <-LISTsynt[x,6]
		LY<-LINKr[,c(subset(c(1:length(LINKr[1,])),
			colnames(LINKr)==STy))]	
		return(
			length(intersect(LX,LY))/length(LX))
	})
}))


COL<-cbind(
seq(1,1,length.out=102),
c(seq(0,1,length.out=41),seq(1,1,length.out=61)),
c(seq(0,0,length.out=41),seq(0,1,length.out=61))
)

COLa<-unique(sapply(
	c(1:length(COL[,1])),function(x){
	COLx<-COL[x,]
	return(rgb(COLx[1], COLx[2], COLx[3]))
}))

SIr<-seq(min(SIr),max(SIr),length.out=length(COLa))

setwd(paste(p0,"Synteny",sep="/"))

pdf("Synteny-plot-AD.pdf",width=5,height=5)

par(mfcol=c(length(LISTsynt[,1]),length(LISTsynt[,1])),mar=c(0,0,0,0),bty="n")

sapply(c(1:length(LISTsynt[,1])),function(x){
	X=x
	STx<-LISTsynt[x,6]
	LX<-LINKr[,c(subset(c(1:length(LINKr[1,])),
		colnames(LINKr)==STx))]
	#Genes as ordered in this genome
	GX<-unlist(sapply(c(2:(length(LX))),
			function(x){
		Lx<-as.numeric(unlist(strsplit(
			LX[x],split="_")))
		Ix<-as.numeric(unlist(strsplit(
			LX[(x-1)],split="_")))
		return(setdiff(Lx,Ix))	
	}))
	GX <-c(GX,
		setdiff(
		as.numeric(unlist(strsplit(
		LX,split="_"))), GX))	
		
	
	sapply(c(1:length(LISTsynt[,1])),function(x){
		Y=x
		
		STN<-unlist(strsplit(LISTsynt[X,1],
			split="_"))
		STN<-paste(
		paste(unlist(strsplit(STN[1],
			split=""))[1],collapse=""),	
		paste(unlist(strsplit(STN[length(STN)],
			split=""))[1:3],collapse=""),
			sep=". ")
		
		
		plot(-10,-10,xlim=c(0,1),ylim=c(0,1),
			xaxt="n",yaxt="n",xlab="",ylab="")	
		text(0.5,0.5,ifelse(X==Y, STN,""),
			cex=1,font=4,
			col=LISTsynt[X,4],srt=0)
		STy <-LISTsynt[x,6]
		LY<-LINKr[,c(subset(c(1:length(LINKr[1,])),
			colnames(LINKr)==STy))]	
		
		SIx<-
			length(intersect(LX,LY))/length(LX)
			
		polygon(c(0,0,1,1),c(0,1,1,0),
			border=NA,lwd=0.01,
			col=ifelse(X>=Y,NA,
			subset(COLa,
			abs(SIx-SIr)==min(abs(SIx-SIr)))[1]))	
			
		#Genes as ordered in this genome
		GY<-unlist(sapply(c(2:(length(LY))),
				function(x){
			Lx<-as.numeric(unlist(strsplit(
				LY[x],split="_")))
			Ix<-as.numeric(unlist(strsplit(
				LY[(x-1)],split="_")))
			return(setdiff(Lx,Ix))	
		}))
		GY <-c(GY,
			setdiff(
			as.numeric(unlist(strsplit(
			LY,split="_"))), GY))
			
		ORDx<-sapply(GX,function(x){
			return(subset(c(1:length(GY)),GY==x))
		})/length(GX)
		points(vx,ORDx,pch=19,cex=0.02,
			col=ifelse(X>=Y,NA,"black"))			
				
	})		
})

dev.off()



#### STEP 6c: Format Synteny tree and calculate RF distances
# 1 # between the synteny best tree and replicates
# 2 # between the synteny best tree, the gene content tree and lineage trees

library(ape)

p0<-"/Users/jean-baptisteleducq/Desktop/New-Methylo-Genomes/Alessa-et-al-2021-genomes/"

setwd(paste(p0,"Synteny",sep="/"))

#Retrieve phylogenetic tree based on annotation occurence
	#read tree as character
	treex<-as.vector(read.table("RAxML_bipartitionsBranchLabels.Synteny-1000_mBINCAT"))
	#split nodes
	treex<-unlist(strsplit(
		as.character(treex),split=":"))
	
	#identify nodes with brackets (node supports)
	bi<-grep("[[]",treex)
	bo<-setdiff(c(1:length(treex)), bi)
	#extract node support and convert in readable tree format: #BranchLength[NodeSupport])# converted in #NodeSupport:BranchLength)# 
	ns<-sapply(bi,function(x){
		nx<-unlist(strsplit(treex[x],split="[[]"))
		nx<-unlist(strsplit(nx,split="[]]"))
		return(paste(nx[2],":",nx[1],nx[3],sep=""))
	})
	treex<-c(ns, treex[bo])[order(c(bi,bo))]
	
	#if last sign is not ")" or ";", add ":"
	bi<-grep("[)]", treex)
	bo<-setdiff(c(1:length(treex)), bi)
	
	ns<-paste(treex[bo],":",sep="")
	
	treex<-c(ns, treex[bi])[order(c(bo,bi))]
	
	treex<-paste(treex,collapse="")
	write.table(treex,"RAxML_bipartitionsBranchLabels.Synteny-1000_mBINCAT_reformated",
		row.names=F,col.names=F,quote=F)
	
	library(ape)
	treex<-read.tree("RAxML_bipartitionsBranchLabels.Synteny-1000_mBINCAT_reformated")	



#Compare RF distance
# 1 between the annotation BINCAT best tree and replicates
# 2 between the annotation BINCAT best tree and lineage trees

#Replicate trees

REP<-read.tree("RAxML_bootstrap.Synteny-1000_mBINCAT")

REP<-c(list(treex), REP)

class(REP) <- "multiPhylo"

setwd(paste(p0,"Trees-RF",sep=""))

write.tree(REP,"Synteny-RAxML-replicates.tre")

LINtree<-read.tree("Lineage-trees-RAS")

LINtree<-c(list(treex),LINtree)

class(LINtree) <- "multiPhylo"

write.tree(LINtree,"Lineage-trees-RAS-Synteny.tre")

#Get tree distances (replace scales by underscores in files)

RFrep<-read.table("Synteny-RAxML-replicates-treedist.txt",header=T)
RFlin<-read.table("Lineage-trees-RAS-Synteny-treedist.txt",header=T)

NORM<-length(treex$tip.label)*2

pdf("RF-dist-Synteny-213-genomes-PAUP.pdf",height=2,width=2)

par(mar=c(4,4,1,1),bty="n")

hist(RFrep[,2]/NORM,breaks=20,border=NA,
	col="grey",xlim=c(min(RFrep[,2]/NORM),
	max(RFlin[,2])/NORM),
	main="",
	cex.axis=0.6,las=1,xlab="Norm. RF distance",cex.lab=0.8,cex.axis=0.8)
	
points(RFlin[,2]/NORM,c(0,0,0),pch=19,
	cex=0.5,col=c("blue","green2","red"))

legend("topright",cex=0.3, horiz=F,
	legend= c("Synteny vs. replicates",
	"RAxML conc. vs. Synteny",
	"Astral vs. Synteny",
	"SVDquartet vs. Synteny"),
	text.col= "black",box.col=NA,pch=c(15,19,19,19),
	border=NA,col= c("grey","blue","green2","red"),
	title="",title.col="black",ncol=1)

dev.off()

#Step 6d: Calculate average synteny index among groups and within group among species

#Alessa genomes annotation file location
p0<-"/Users/jean-baptisteleducq/Desktop/New-Methylo-Genomes/Alessa-et-al-2021-genomes/"

setwd(p0)

#Strain characteristics
ID<-read.table("Table-Correspondance-ID.txt",header=T)
SPE<-read.table("Species-tree-order.txt",header=T)
SPE<-SPE[order(SPE[,1]),]

#Synteny index
setwd(paste(p0,"Synteny",sep="/"))
SI<-read.table("Synteny-Index.txt",header=T)





#Type of comparisons


SIcat<-unlist(sapply(c(1:(length(SI[,1])-1)),function(x){
	X=x
	stx<-row.names(SI)[x]
	spx<-SPE[grep(paste(stx,"_",sep=""),
			paste(SPE[,3],"_",sep="")),2]
	gx<-SPE[grep(paste(stx,"_",sep=""),
			paste(SPE[,3],"_",sep="")),4]
	sapply(c((X+1):(length(SI[,1]))),function(x){
		Y=x
		sty<-row.names(SI)[x]
		spy<-SPE[grep(paste(sty,"_",sep=""),
			paste(SPE[,3],"_",sep="")),2]
		gy<-SPE[grep(paste(sty,"_",sep=""),
			paste(SPE[,3],"_",sep="")),4]
		return(paste(
			ifelse(spx==spy,paste("spe",gx,sep="-"),
			ifelse(gx==gy,gx,
			paste(sort(c(gx,gy)),collapse="-"))),
			unlist(SI[X,Y]),sep="="))
	})		
}))

SIcat<-matrix(unlist(strsplit(SIcat,
	split="=")),ncol=2,byrow=T)
	
SIgs<-t(sapply(sort(unique(SPE[,4])),function(x){
	spex<-as.numeric(subset(SIcat[,2],
		SIcat[,1]==paste("spe-",x,sep="")))
	gx<-as.numeric(subset(SIcat[,2],
		SIcat[,1]==x))
	return(c(
		paste(round(c(mean(spex),
			sd(spex)),2),collapse=" ± "),
		paste(round(c(mean(gx),
			sd(gx)),2),collapse=" ± ")))
}))

colnames(SIgs)<-c("within_sp","among_sp")

CATg<-sapply(sort(unique(SPE[,4])),function(x){
	catx<-paste(x,sort(unique(SPE[,4])),sep="-")
	return(sapply(catx,function(x){
		gx<-as.numeric(subset(
			SIcat[,2],SIcat[,1]==x))
		return(paste(round(c(mean(gx),
			sd(gx)),2),collapse=" ± "))
	}))
})

rownames(CATg)<-rownames(SIgs)

SIgs<-cbind(SIgs,CATg)

SIgs<-ifelse(SIgs=="NaN ± NA","",SIgs)

setwd(paste(p0,"Synteny",sep="/"))

write.table(SIgs,"Summary-SI-per-Group.txt")

#######STEP 7: SUMMARY FIGURE (species tree + annotation + synteny + eco distribution)

p0<-"/Users/jean-baptisteleducq/Desktop/New-Methylo-Genomes/Alessa-et-al-2021-genomes/"

#Synteny index
setwd(paste(p0,"Synteny",sep="/"))
SI<-read.table("Synteny-Index.txt",header=T)

#Annotation dissimilarity
setwd(paste(p0,"Gene-content",sep="/"))
ANNO<-read.table("Anno-Bray.txt",header=T)

#Species ordered according to the ASTRAL tree
setwd(p0)
SPE<-read.table("Species-tree-order.txt",header=T)
SPE<-SPE[order(SPE[,1]),]


#Strain characteristics
ID<-read.table("Table-Correspondance-ID.txt",header=T)

#Calculate average synteny index per pairwise species comparison
SIspe<-sapply(c(1:length(SPE[,1])),function(x){
	SPEx<-SPE[x,2]
	STx<-unlist(strsplit(SPE[x,3],split="_"))
	ix<-sapply(STx,function(x){
		return(subset(c(1:length(SI[,1])),
			rownames(SI)==x))
	})
	return(sapply(c(1:length(SPE[,1])),function(x){
		SPEy<-SPE[x,2]
		STy<-unlist(strsplit(SPE[x,3],split="_"))
		iy<-sapply(STy,function(x){
			return(subset(c(1:length(SI[,1])),
				rownames(SI)==x))
		})
		return(mean(unlist(SI[ix,iy])))
	}))
})

Nspe<-paste(SPE[,4],SPE[,2],sep="=")
colnames(SIspe)<-Nspe
rownames(SIspe)<-Nspe

setwd(paste(p0,"Synteny",sep="/"))

write.table(SIspe,"SI-per-species.txt")

OUT<-c("Meth_jeotgali","Meth_planium","Meth_soli","Meth_oxalidis","Meth_durans","Meth_segetis","Meth_trifolii")

Ax<-grep("A=", Nspe)
Bx<-grep("B=", Nspe)
Cx<-grep("C=", Nspe)
Dx<-grep("D=", Nspe)
Ex<-grep("Enterovirga=", Nspe)
Mx<-grep("Microvirga=", Nspe)

Ox<-sapply(OUT,function(x){
	return(grep(x, Nspe))
})	
	
apply(SIspe[Ox,setdiff(Ax, Ox)],1,mean)
apply(SIspe[Ox,setdiff(Ax, Ox)],1,sd)

mean(SIspe[Ox,setdiff(Dx, Ox)])
sd(SIspe[Ox,setdiff(Dx, Ox)])

#Calculate average BC index of dissimilarity in gene content per pairwise species comparison
ANNOspe<-sapply(c(1:length(SPE[,1])),function(x){
	SPEx<-SPE[x,2]
	STx<-unlist(strsplit(SPE[x,3],split="_"))
	ix<-sapply(STx,function(x){
		return(subset(c(1:length(ANNO[,1])),
			rownames(SI)==x))
	})
	return(sapply(c(1:length(SPE[,1])),function(x){
		SPEy<-SPE[x,2]
		STy<-unlist(strsplit(SPE[x,3],split="_"))
		iy<-sapply(STy,function(x){
			return(subset(c(1:length(ANNO[,1])),
				rownames(SI)==x))
		})
		return(mean(unlist(ANNO[ix,iy])))
	}))
})

colnames(ANNOspe)<-Nspe
rownames(ANNOspe)<-Nspe


setwd(paste(p0,"Gene-content",sep="/"))

write.table(ANNOspe,"BC-per-species.txt")

#Color scales

#SI

COLs<-cbind(c(seq(0,1,length.out=26),
		seq(1,1,length.out=26),
		seq(1,1,length.out=51)),
	c(seq(0,0,length.out=26),
		seq(0,1,length.out=26),
		seq(1,1,length.out=51)),
	c(seq(0,0,length.out=52),
		seq(0,1,length.out=51))
)

COLs <-unique(sapply(
	c(1:length(COLs[,1])),function(x){
	COLx<-COLs[x,]
	return(rgb(COLx[1], COLx[2], COLx[3]))
}))

#SI range
SIr<-seq(min(SIspe),max(SIspe),length.out=length(COLs))

#Annotations

COLa<-cbind(
c(seq(1,0,length.out=31),seq(0,0.25,length.out=31),seq(0.25,0,length.out=41)),
c(seq(1,1,length.out=31),seq(1,0.25,length.out=31),seq(0.25,0,length.out=41)),
c(seq(1,1,length.out=62),seq(1,0,length.out=41))
)

COLa<-unique(sapply(
	c(1:length(COLa[,1])),function(x){
	COLx<-COLa[x,]
	return(rgb(COLx[1], COLx[2], COLx[3]))
}))

#ANNOTATION range
ANNOr<-seq(min(ANNOspe),max(ANNOspe),length.out=length(COLa))

#Size of rigth and bottom margin

MAR=15

##Define categories of biomes

BIOME<-cbind(
c("Air","Anthropo","Endophyte","Food","Fungus","Lake","Lichen","Microbiome","Ocean","Phyllosphere","Plant","Rhizhosphere","Sediments","Seed","Skin","Soil","Space","Unknown","Water"),
c("Air","Human","Endoph.","Human","Fung./Lich.","Sed./Water","Fung./Lich.","Human","Ocean","Phyllo.","Plant","Rhizho.","Sed./Water","Seed","Human","Soil","Human","Unknown","Sed./Water")
)

COLbiome<-cbind(
c("Air","Human","Endoph.","Fung./Lich.","Ocean","Phyllo.","Plant","Rhizho.","Sed./Water","Seed","Soil","Unknown"),
c("cyan","black","yellow4","orange","blue4","green","green3","red","blue","yellow2","brown","grey")
)

#Relative Abundance of each species in each biome

Fbiomes<-t(sapply(c(1:length(SPE[,1])),function(x){
	SPEx<-SPE[x,2]
	#strains in this species
	STx<-unlist(strsplit(SPE[x,3],split="_"))
	Bx<-apply(sapply(c(STx,STx),function(x){
		#biome(s) where this strain come from
		bx<-unlist(strsplit(subset(ID[,7],ID[,2]==x),split="/"))
		#relative abundance
		fx<-unique(ifelse(bx=='-',0,1)/length(bx))
		bx<-sapply(BIOME[,1],function(x){
			return(length(intersect(x,bx))* fx)
		})
		bx<-sapply(COLbiome[,1],function(x){
			return(sum(subset(bx,BIOME[,2]==x)))
		})
		
	}),1,mean)
	return(Bx)
}))

write.table(cbind(SPE,Fbiomes),"Relative-Biome-Fr-per-species.txt",row.names=F)

Gst<-unlist(sapply(c(1:length(SPE[,1])),function(x){
	return(paste(SPE[x,4],SPE[x,2],unlist(strsplit(SPE[x,3],split="_")),sep="="))
}))

Gst<-matrix(unlist(strsplit(Gst,
	split="=")),ncol=3,byrow=T)


PLANT<-unique(c(
	grep("Endophyte",ID[,7]),
	grep("Phyllosphere",ID[,7]),
	grep("Plant",ID[,7]),
	grep("Rhizhosphere",ID[,7]),
	grep("Seed",ID[,7])
))
PLANTst<-ID[PLANT,2]
PLANTg<-sapply(PLANTst,function(x){
	return(subset(Gst[,1], Gst[,3]==x))	
})
PLANTg<-sapply(c("A","B","C","D"),function(x){
	return(length(subset(PLANTg, PLANTg ==x)))
})/length(PLANTg)


BIOMEg<-sapply(BIOME[,1],function(x){
	bx<-grep(x,ID[,7])
	stx<-ID[bx,2]
	Gx<-sapply(stx,function(x){
		return(subset(Gst[,1], Gst[,3]==x))	
	})
	return(sapply(c("A","B","C","D"),function(x){
		return(length(subset(Gx, Gx==x)))
	})/length(Gx))
})

#Relative biome abundance per group (Proportion of species that were isolated at least once in each biome)

FbiomesG<-t(sapply(sort(unique(SPE[,4])),function(x){
	Fx<-subset(Fbiomes,SPE[,4]==x)
	Fx<-ifelse(Fx==0,0,1)
	return(apply(Fx,2,mean))
}))

FbiomesG2<-t(sapply(sort(unique(SPE[,4])),function(x){
	Fx<-subset(Fbiomes,SPE[,4]==x)
	return(apply(Fx,2,mean))
}))

write.table(round(FbiomesG,4),"Relative-Biome-Fr-per-Group.txt")

#Limits for biome plot
B1=(-MAR-2)
B2=(-6)

#Summary figure
setwd(p0)

pdf("Summary-Synteny-Annotations.pdf",width=10,height=10)

par(mar=c(0,0,0,0),bty="n")

plot(-10,-10,xlim=c(length(SPE[,1]),-MAR),
	ylim=c(length(SPE[,1])+MAR,0),
	xaxt="n",yaxt="n",xlab="",ylab="")

#BIOMES

#sapply(c(1:length(SPE[,1])),function(x){
#	X=x
#	SPEx<-SPE[x,2]
#	z=length(grep("Meth_sp_",SPEx))
#	z=ifelse(z==1,"*","")
#	text(B2+1,X-0.5,z,cex=1.1,font=2,srt=90)
#	Fx<-cumsum(Fbiomes[x,])
#	Gx<-c(0,Fx[1:(length(Fx)-1)])
#	Fx=B1+(B2-B1)*Fx
#	Gx=B1+(B2-B1)*Gx
#	sapply(c(1:length(Fx)),function(x){
#		colx<-COLbiome[x,2]
#		polygon(c(Fx[x],Fx[x],Gx[x],Gx[x]),
#			c(X-1,X,X,X-1),border=NA,col=colx)
#	})
#})
#Legends



Y1=6+length(SPE[,1])
Y2=10+length(SPE[,1])
s2=25
s1=50
a2=75
a1=100
xs<-seq(s1,s2,length.out=(length(COLs)+1))
ys<-c(Y2,Y1)
sapply(c(length(COLs):1),function(x){
	polygon(c(xs[x],xs[x],xs[x+1],xs[x+1]),
		c(ys[1],ys[2],ys[2],ys[1]),
		col=COLs[x],border=NA)
})
polygon(c(s1,s1,s2,s2),
		c(Y1,Y2,Y2,Y1),
		col=NA,border="black",lwd=0.5)
text(c(s1,s2),c(mean(ys),mean(ys)),
	round(c(min(SIr),max(SIr)),2),
	cex=1,pos=c(2,4),offset=0.5)		

xa<-seq(a1,a2,length.out=(length(COLa)+1))
ya<-c(Y2,Y1)
sapply(c(1:length(COLa)),function(x){
	polygon(c(xa[x],xa[x],xa[x+1],xa[x+1]),
		c(ya[1],ya[2],ya[2],ya[1]),
		col=COLa[x],border=NA)
})
polygon(c(a1,a1,a2,a2),
		c(Y1,Y2,Y2,Y1),
		col=NA,border="black",lwd=0.5)
text(c(a1,a2),c(mean(ya),mean(ya)),
	round(c(min(ANNOr),max(ANNOr)),2),
	cex=1,pos=c(2,4),offset=0.5)

sapply(c(1:length(SPE[,1])),function(x){
	X=x
	SIx<-SIspe[,x]
	ANNOx<-ANNOspe[,x]
		
	sapply(c(1:length(SPE[,1])),function(x){
		Y=x
		SIy<-SIx[x]
		ANNOy<-ANNOx[x]
		
		cols<-subset(COLs,
			abs(SIy-SIr)==min(
			abs(SIy-SIr)))[1]
		cola<-subset(COLa,
			abs(ANNOy-ANNOr)==min(
			abs(ANNOy-ANNOr)))[1]		
			
		polygon(c(X-1,X-1,X,X),c(Y-1,Y,Y,Y-1),
			border=NA,lwd=0.01,
			col=ifelse(X<Y, cols,
				ifelse(Y<X, cola,NA)))			
				
	})	
	return("")	
})


#par(new=T,mar=c(3,0,0,0),bty="n")

#plot(-10,-10,xlim=c(0,1),
	ylim=c(0,1),
	xaxt="n",yaxt="n",xlab="",ylab="")

#legend("bottomright",cex=0.95, horiz=F,
#	legend= COLbiome[,1],
#	text.col= "black",box.col=NA,pch=15,
#	border=NA,col= COLbiome[,2],
#	title="Biomes",title.col="black",ncol=1,
#	y.intersp=0.75)

dev.off()

############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
