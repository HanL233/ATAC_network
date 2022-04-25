library("ggplot2")
library("ggsignif")
library("FactoMineR")
library("factoextra")
library("ggfortify")
library("survival")
library("pheatmap")
library('ggthemes')
require("RColorBrewer")
require("graphics")

markOT <- function(arr){
	cat(paste('length=',length(arr),'\n'))
        ql <- quantile(arr, probs = 0.25)
        qu <- quantile(arr, probs = 0.75)
        ul <- qu-ql
        m <- rep(T,length(arr))
        m[arr>qu+1.5*ul] <- F
        m[arr<ql-1.5*ul] <- F
        return(m)
}

data <- read.table('./MEAN.txt',header=T)
row.names(data) <- data[,1]
data <- data[,-1]
data <- t(data)

clin <- read.table('./clinical.all.MEAN.txt',header=T,sep="\t")
row.names(clin) <- clin[,1]
clin <- clin[,-1]

sam <- intersect(row.names(data),row.names(clin))
clin <- clin[match(sam,row.names(clin)),]
data <- data[match(sam,row.names(data)),]
clin$PFS <- unlist(lapply(seq(nrow(clin)),function(x){if(clin$OS[x]==1 | clin$PFI[x]==1){1}else{0}}))
clin$PFS.time <- unlist(lapply(seq(nrow(clin)),function(x){if(clin$OS[x]==0 & clin$PFI[x]==0){max(clin$OS.time[x],clin$PFI.time[x])}else{if(clin$OS[x]==1 & clin$PFI[x]==0){clin$OS.time[x]}else{if(clin$OS[x]==0 & clin$PFI[x]==1){clin$PFI.time[x]}else{min(clin$OS.time[x],clin$PFI.time[x])}}}}))

type <- unlist(lapply(row.names(data),function(x){unlist(strsplit(x,"_"))[1]}))

res.pca <- PCA(data,graph=F)
pic <- fviz_pca_ind(res.pca,geom.ind="point",col.ind=type,pointshape=20,pointsize=2) + theme(legend.position = "alle") + scale_x_discrete() + scale_y_discrete()+theme_minimal()
mycolor <- ggplot_build(pic)$data[[2]]
utype <- unique(type)
for(i in seq(length(mycolor$x))){
	pic <- pic + annotate("text",x=mycolor$x[i],y=mycolor$y[i],label=utype[i],hjust=0.5,vjust=0.5,cex=3.5,colour=mycolor$colour[i])
}
ggsave(file=paste("pca.panCancer.pdf",sep="."),width=4,height=3.25, plot = pic)
write.csv(file="./contribution.panCancer.csv",res.pca$var$contrib[order(res.pca$var$contrib[,1],decreasing=T),])

clin <- read.csv('./clinical.LUAD.csv',header=T)
row.names(clin) <- clin[,1]
clin <- clin[,-1]
sam <- intersect(row.names(data),row.names(clin))
clin <- clin[match(sam,row.names(clin)),]
data <- data[match(sam,row.names(data)),]
clin$PFS <- unlist(lapply(seq(nrow(clin)),function(x){if(clin$OS[x]==1 | clin$PFI[x]==1){1}else{0}}))
clin$PFS.time <- unlist(lapply(seq(nrow(clin)),function(x){if(clin$OS[x]==0 & clin$PFI[x]==0){max(clin$OS.time[x],clin$PFI.time[x])}else{if(clin$OS[x]==1 & clin$PFI[x]==0){clin$OS.time[x]}else{if(clin$OS[x]==0 & clin$PFI[x]==1){clin$PFI.time[x]}else{min(clin$OS.time[x],clin$PFI.time[x])}}}}))

type <- unlist(lapply(row.names(data),function(x){unlist(strsplit(x,"_"))[1]}))



data2 <- data[which(type=="LUAD"),]
A <- as.numeric(as.character(clin$tobacco_smoking_pack_years_smoked))[type=="LUAD"]
A1 <- ((!is.na(A) & A<20) | clin$tobacco_smoking_history==1)
A2 <- ((!is.na(A) & A>=20) & clin$tobacco_smoking_history!=1)
res.pca <- PCA(data2,graph=F)
type2 <- rep("NA",sum(type=="LUAD"))
type2[A1] <- "Light smoker"
type2[A2] <- "Heavy smoker"
x <- res.pca$ind$coord[,1]
y <- res.pca$ind$coord[,2]
B <- abs(y)
print(mean(B))
pic <- fviz_pca_ind(res.pca,geom.ind="point",col.ind=type2,pointshape=20,pointsize=5)+geom_hline(yintercept=c(-1*mean(B),mean(B)),color="orange")+theme_minimal()
ggsave(file=paste("pca.PC2.LUAD.pdf",sep="."),width=4,height=3.25, plot = pic)

tp <- rep(NA,length(y))
tp[A1] <- "Light smoker"
tp[A2] <- "Heavy smoker"
A <- data.frame(distance=B[!is.na(tp)],type=tp[!is.na(tp)])
A <- A[-1*which((A$distance==A$distance[A$type=="Light smoker"][!markOT(A$distance[A$type=="Light smoker"])]) & A$type=="Light smoker"),]
print(table(tp))
pic <- ggplot(A,aes(y=distance,x=type,fill=type))+geom_boxplot(outlier.shape = NA) + geom_signif(aes(y=distance,x=type),comparisons = list(c("Light smoker", "Heavy smoker")),test=t.test)+theme_minimal()
ggsave(file=paste("./boxplot.PC2.LUAD.pdf",sep="."),width=3,height=3, plot = pic)
t.test(A$distance[A$type=="Light smoker"],A$distance[A$type=="Heavy smoker"])

type2 <- clin$PFS.time[which(type=="LUAD")]
ind <- rep(T,length(y))#markOT(y) & markOT(type2)
OS <- data.frame(OS=clin$PFS[type=='LUAD'][ind],OS.time=clin$PFS.time[type=='LUAD'][ind],OS2=clin$OS[type=='LUAD'][ind],OS2.time=clin$OS.time[type=='LUAD'][ind],type=B[ind]>=mean(B))
OS$type[OS$type==T] <- "Far"
OS$type[OS$type==F] <- "Near"
A <- survdiff(formula = Surv(OS.time, OS) ~ type, data = OS)
p <- pchisq(A$chisq,df=1,lower.tail=F)
pic <- autoplot(survfit(Surv(OS.time, OS) ~ type, data = OS))+annotate("text",x=1000,y=0.2,label=paste("p=",format(p,digits=2),sep=""),hjust=0,vjust=0,cex=4)+annotate("text",x=1000,y=0.12,label=paste("n(Near)=",sum(OS$type=="Near"),sep=""),hjust=0,vjust=0,cex=4)+annotate("text",x=1000,y=0.04,label=paste("n(Far)=",sum(OS$type=="Far"),sep=""),hjust=0,vjust=0,cex=4)+theme_minimal()
ggsave(file=paste("./surv.PC2.OS.LUAD.pdf",sep="."),width=3,height=3, plot = pic)

A <- survdiff(formula = Surv(OS2.time, OS2) ~ type, data = OS)
p <- pchisq(A$chisq,df=1,lower.tail=F)
pic <- autoplot(survfit(Surv(OS2.time, OS2) ~ type, data = OS))+annotate("text",x=1000,y=0.2,label=paste("p=",format(p,digits=2),sep=""),hjust=0,vjust=0,cex=4)+annotate("text",x=1000,y=0.12,label=paste("n(Near)=",sum(OS$type=="Near"),sep=""),hjust=0,vjust=0,cex=4)+annotate("text",x=1000,y=0.04,label=paste("n(Far)=",sum(OS$type=="Far"),sep=""),hjust=0,vjust=0,cex=4)+theme_minimal()
ggsave(file=paste("./surv.PC2.PFS.LUAD.pdf",sep="."),width=3,height=3, plot = pic)


B <- abs(x)
res.pca <- PCA(data2,graph=F)
type2 <- rep("NA",sum(type=="LUAD"))
type2[A1] <- "Light smoker"
type2[A2] <- "Heavy smoker"

pic <- fviz_pca_ind(res.pca,geom.ind="point",col.ind=type2,pointshape=20,pointsize=5)+geom_vline(xintercept=c(-1*mean(B),mean(B)),color="orange")+theme_minimal()
ggsave(file=paste("pca.PC1.LUAD.pdf",sep="."),width=4,height=3.25, plot = pic)
A <- data.frame(distance=B[!is.na(tp)],type=tp[!is.na(tp)])
#A <- A[-1*which((A$distance==A$distance[A$type=="Light smoker"][!markOT(A$distance[A$type=="Light smoker"])]) & A$type=="Light smoker"),]
A <- A[-1*which((A$distance==A$distance[A$type=="Heavy smoker"][!markOT(A$distance[A$type=="Heavy smoker"])]) & A$type=="Heavy smoker"),]
print(table(tp))
#pic <- ggplot(A)+geom_boxplot(aes(y=distance,x=type,fill=type)) + geom_signif(aes(y=distance,x=type),comparisons = list(c("Light smoker", "Heavy smoker")),test='t.test')
#pic <- ggplot(A,aes(y=distance,x=type,fill=type))+geom_boxplot()+geom_jitter()+geom_point() + geom_signif(aes(y=distance,x=type),comparisons = list(c("Light smoker", "Heavy smoker")),test=t.test)
pic <- ggplot(A,aes(y=distance,x=type,fill=type))+geom_boxplot() + geom_signif(aes(y=distance,x=type),comparisons = list(c("Light smoker", "Heavy smoker")),test=t.test)+theme_minimal()

ggsave(file=paste("./boxplot.PC1.LUAD.pdf",sep="."),width=3,height=3, plot = pic)
t.test(A$distance[A$type=="Light smoker"],A$distance[A$type=="Heavy smoker"])

type2 <- clin$PFS.time[which(type=="LUAD")]
#ind <- !is.na(tp)#rep(T,length(y))#markOT(y) & markOT(type2)
ind <- rep(T,length(y))#markOT(y) & markOT(type2)
OS <- data.frame(OS=clin$PFS[type=='LUAD'][ind],OS.time=clin$PFS.time[type=='LUAD'][ind],OS2=clin$OS[type=='LUAD'][ind],OS2.time=clin$OS.time[type=='LUAD'][ind],type=B[ind]>=mean(B))
OS$type[OS$type==T] <- "Far"
OS$type[OS$type==F] <- "Near"
A <- survdiff(formula = Surv(OS.time, OS) ~ type, data = OS)
p <- pchisq(A$chisq,df=1,lower.tail=F)
pic <- autoplot(survfit(Surv(OS.time, OS) ~ type, data = OS))+annotate("text",x=1500,y=0.7,label=paste("p=",format(p,digits=2),sep=""),hjust=0,vjust=0,cex=4)+annotate("text",x=1500,y=0.52,label=paste("n(Near)=",sum(OS$type=="Near"),sep=""),hjust=0,vjust=0,cex=4)+annotate("text",x=1500,y=0.405,label=paste("n(Far)=",sum(OS$type=="Far"),sep=""),hjust=0,vjust=0,cex=4)+theme_minimal()
ggsave(file=paste("./surv.PC1.OS.LUAD.pdf",sep="."),width=3,height=3, plot = pic)

A <- survdiff(formula = Surv(OS2.time, OS2) ~ type, data = OS)
p <- pchisq(A$chisq,df=1,lower.tail=F)
pic <- autoplot(survfit(Surv(OS2.time, OS2) ~ type, data = OS))+annotate("text",x=1500,y=0.7,label=paste("p=",format(p,digits=2),sep=""),hjust=0,vjust=0,cex=4)+annotate("text",x=1500,y=0.52,label=paste("n(Near)=",sum(OS$type=="Near"),sep=""),hjust=0,vjust=0,cex=4)+annotate("text",x=1500,y=0.405,label=paste("n(Far)=",sum(OS$type=="Far"),sep=""),hjust=0,vjust=0,cex=4)+theme_minimal()
ggsave(file=paste("./surv.PC1.PFS.LUAD.pdf",sep="."),width=3,height=3, plot = pic)





write(file="contribution.LUAD.txt",row.names(res.pca$var$contrib)[order(res.pca$var$contrib[,1],decreasing=T)])
