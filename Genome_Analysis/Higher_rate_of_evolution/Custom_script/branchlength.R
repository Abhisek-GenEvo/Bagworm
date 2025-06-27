args = commandArgs(trailingOnly=TRUE)
require(ape)
t = read.tree(args[1])
name = as.character()
bl = as.numeric()
require("adephylo")
for(i in c(1:length(t$tip.label))) {
    name = c(name, t$tip.label[i])
    val = distRoot(t, t$tip.label[i])
    bl = c(bl, val)
}
data = data.frame(sp = name, brl = bl)
write.table(data, file=paste(args[1], "tmp", sep="_"), quote=F, col.names=F, row.names=F, sep="\t")
