#!/usr/bin/Rscript

library(ape)

con = file("~/yytao/16S/16S.names.txt", "r")
while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    tr <- read.tree(file = paste("~/yytao/16S/remove1/epang_insert1_bi.nwk.",line,sep=""))
	newtr <- di2multi(tr)
    write.tree(newtr, file = paste("~/yytao/16S/remove1/epang_insert1_multi.nwk.",line,sep=""))
}
close(con)
