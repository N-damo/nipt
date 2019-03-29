library(stats)
result<- data.frame(chr=c('chr13','chr18','chr21'),pvalue=c(0,0,0))
chr13<- read.table('chr13',sep=',',head=TRUE)
artificialchr13<- read.table('artificialchr13',sep=',',head=TRUE)
chr13<- wilcox.test(chr13$RC, artificialchr13$RC)
result[result$chr=='chr13',]$pvalue=chr13$p.value
chr18<- read.table('chr18',sep=',',head=TRUE)
artificialchr18<- read.table('artificialchr18',sep=',',head=TRUE)
chr18<- wilcox.test(chr18$RC, artificialchr18$RC)
result[result$chr=='chr18',]$pvalue=chr18$p.value
chr21<- read.table('chr21',sep=',',head=TRUE)
artificialchr21<- read.table('artificialchr21',sep=',',head=TRUE)
chr21<- wilcox.test(chr21$RC, artificialchr21$RC)
result[result$chr=='chr21',]$pvalue=chr21$p.value
write.table(result,'trisomy.txt',sep='\t',row.names = FALSE,quote =FALSE)
