mydata2 = read.csv(file= "file2.txt", sep="\t", header=F)
p = mydata2$V2
a = p.adjust(p, method = "fdr", n = length(p))
mydata2$V3 = a
write.table(file="file2_adj.txt", mydata2, sep="\t", row.names=F, col.names=F, quote=F)
