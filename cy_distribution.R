ngs=read.table("AML_001_WGS_novobreak_SV_sort.bed",sep="\t",stringsAsFactors=F)

#
ngs = inter[,4]
ngs[ngs=="<DEL>"]="Deletion"
ngs[ngs=="<DUP>"]="Duplication"
ngs[ngs=="<INV>"]="Inversion"
ngs[ngs=="<TRA>"]="Translocation"
###########


ngs_table=c(ngs_table[1:3],ngs_table[4],0)
names(ngs_table)=c("Deletion","Duplication","Inversion","Trans","Nov_Ins")

pdf("NGS_SV")
barplot(ngs_table,border=NA,col="salmon",main="NGS",cex.main=2,ylab="Number of events",cex.lab=1.3,ylim=c(0,130000))
dev.off()

#####################################################################################################



nano=read.table("tmpfnaif89jks0c_sort_4-cov-70_sort.bed",sep="\t",stringsAsFactors=F)

nano = nano[,6]

nano[nano=="E-Nov_Ins_bp"]="Novel Insertion"
nano[nano=="S-Nov_Ins_bp"]="Novel Insertion"
nano[nano=="Nov_Ins"]="Novel Insertion"

nano[nano=="Inter-Ins(1)"]="Translocation"
nano[nano=="Inter-Ins(2)"]="Translocation"
nano[nano=="InterTx"]="Translocation"
nano[nano=="Intra-Ins"]="Translocation"
nano[nano=="Intra-Ins(1)"]="Translocation"
nano[nano=="Intra-Ins(2)"]="Translocation"

#nano[nano=="Inter-Ins(1)"]="Insertion"
#nano[nano=="Inter-Ins(2)"]="Insertion"
#nano[nano=="InterTx"]="Tra/Ins"
#nano[nano=="Intra-Ins"]="Tra/Ins"
#nano[nano=="Intra-Ins(1)"]="Insertion"
#nano[nano=="Intra-Ins(2)"]="Insertion"

#nano[nano=="Inter-Ins(2)"]="TRASH"
#nano[nano=="Intra-Ins(2)"]="TRASH"
#nano[nano=="S-Nov_Ins"]="TRASH"

nano[nano=="Del"]="Deletion"
nano[nano=="Inv"]="Inversion"
nano[nano=="Inv(1)"]="Inversion"
nano[nano=="Inv(2)"]="Inversion"
nano[nano=="TDupl"]="Duplication"
#

###########

nano_table=c(nano_table[1:3],nano_table[5],nano_table[4])
names(nano_table)=c("Deletion","Duplication","Inversion","Trans","Nov_Ins")

pdf("Nanopore_SV")
barplot(nano_table,border=NA,col="darkturquoise",main="Hybrid-Seq",cex.main=2,ylab="Number of events",cex.lab=1.3,ylim=c(0,700))
dev.off()

###################################################################################################################
inter=read.table("AML_001_run1_nanoSV-WGS-intersect.tsv",sep="\t",stringsAsFactors=F)

inter_nano = inter[,6]

inter_nano[inter_nano=="E-Nov_Ins_bp"]="Novel Insertion"
inter_nano[inter_nano=="S-Nov_Ins_bp"]="Novel Insertion"
inter_nano[inter_nano=="Nov_Ins"]="Novel Insertion"

inter_nano[inter_nano=="Inter-Ins(1)"]="Translocation"

inter_nano[inter_nano=="Inter-Ins(1)"]="Translocation"
inter_nano[inter_nano=="Inter-Ins(2)"]="Translocation"
inter_nano[inter_nano=="InterTx"]="Translocation"
inter_nano[inter_nano=="Intra-Ins"]="Translocation"
inter_nano[inter_nano=="Intra-Ins(1)"]="Translocation"
inter_nano[inter_nano=="Intra-Ins(2)"]="Translocation"

inter_nano[inter_nano=="Del"]="Deletion"
inter_nano[inter_nano=="Inv"]="Inversion"
inter_nano[inter_nano=="Inv(1)"]="Inversion"
inter_nano[inter_nano=="Inv(2)"]="Inversion"
inter_nano[inter_nano=="TDupl"]="Duplication"

#
inter_ngs = inter[,11]
inter_ngs[inter_ngs=="<DEL>"]="Deletion"
inter_ngs[inter_ngs=="<DUP>"]="Duplication"
inter_ngs[inter_ngs=="<INV>"]="Inversion"
inter_ngs[inter_ngs=="<TRA>"]="Translocation"
###########
colu=cbind(inter_ngs,inter_nano)

#colu_noTrash=colu[!(colu[,2]=="TRASH"),]

inter_nano_table=table(colu[,2])
inter_ngs_table=table(colu[,1])

inter_ngs_table=c(inter_ngs_table[1:4],0)
names(inter_ngs_table)=c("Deletion","Duplication","Inversion","Trans","Nov_Ins")

inter_nano_table=c(inter_nano_table[1:3],inter_nano_table[5],inter_nano_table[4])
names(inter_nano_table)=c("Deletion","Duplication","Inversion","Trans","Nov_Ins")

pdf("ngs_intersect")
#par(mfrow=c(1,2))
barplot(inter_ngs_table,border=NA,col="salmon",main="NGS Intersection",cex.main=2,ylab="Number of events",cex.lab=1.3,ylim=c(0,120))
dev.off()

pdf("nano_intersect")
barplot(inter_nano_table,border=NA,col="darkturquoise",main="Hybrid-Seq Intersection",cex.main=2,ylab="Number of events",cex.lab=1.3,ylim=c(0,120))
dev.off()
