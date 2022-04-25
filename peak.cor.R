args <- commandArgs();
if(length(args) < 6 || length(args) > 8 || !file.exists(args[6])){
	write(paste("Rscript",unlist(strsplit(args[4],"="))[2],"<FI:data.txt>","[NUM:from=1]","[NUM:to=nrow]"),stderr());
        q();
}
data_file <- args[6]
FROM <- 1
TO <- -1
if(length(args) > 6){
        FROM <- as.numeric(args[7])
}
if(length(args) > 7){
        TO <- as.numeric(args[8])
}
N <- 100 #genes with top N highest correlation (both positive & negative) will be listed

markOT <- function(arr){
        ql <- quantile(arr, probs = 0.25)
        qu <- quantile(arr, probs = 0.75)
        ul <- qu-ql
        m <- rep(T,length(arr))
        m[arr>qu+1.5*ul] <- F
        m[arr<ql-1.5*ul] <- F
        return(m)
}

data <- read.table(data_file,header=T)
row.names(data) <- data[,1]
data <- data[,c(-1)]

FROM <- max(FROM,1)
FROM <- min(FROM,nrow(data))
TO <- min(TO,nrow(data))
if(TO<0){
        TO <- nrow(data)
}


idx <- array()
for(i in seq(nrow(data))){
	idx <- rbind(idx,markOT(as.numeric(data[i,])))
}
idx <- idx[-1,]

cat(paste(date(),"\n",sep=""),file=stderr())
rng <- seq(nrow(idx))
re <- lapply(seq(from=FROM,to=TO),function(x){
	co <- unlist(lapply(rng,function(y){
		ix <- which(idx[x,] & idx[y,])
		A <- as.numeric(data[x,])[ix]
		B <- as.numeric(data[y,])[ix]
		cor.test(A,B)$estimate
	}))
	od <- order(co,decreasing=T)
	rd <- c(head(od,n=N),tail(od,n=N))
	write(paste("#",paste(row.names(data)[rd],collapse=","),sep=","),stdout())
	write(paste(co[rd],collapse=","),stdout())
	cat(paste(date(),"\n",sep=""),file=stderr())
})

