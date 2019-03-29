RC_DT<- read.table('niptest',sep='\t',head=TRUE)
RC_DT<- RC_DT[order(RC_DT$GC),]
gcCount.loess <- loess(RC~GC,data=RC_DT,control = loess.control(surface = 'direct'),degree=2)
predictions1<- predict(gcCount.loess,RC_DT$GC)
m <- median(RC_DT$RC)
f<- m/predictions1
RR <- RC_DT$RC*f
RC_DT$RC<- RR
RC_DT<- RC_DT[order(RC_DT$MP),]
mpCount.loess <- loess(RC~MP,data=RC_DT,control = loess.control(surface = 'direct'),degree=2)
predictions2<- predict(mpCount.loess,RC_DT$MP)
m <- median(RC_DT$RC)
f<- m/predictions2
RR <- RC_DT$RC*f
RC_DT$RC<- RR
write.table(RC_DT,'nipt.txt',sep='\t',row.names = FALSE,quote =FALSE)
