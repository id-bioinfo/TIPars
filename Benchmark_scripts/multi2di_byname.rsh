#!/usr/bin/Rscript

library(ape)

con = file("~/yytao/16S/16S.names.txt", "r")
while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    tr <- read.tree(file = paste("~/yytao/16S/remove1/remove1.nwk.",line,sep=""))
	newtr <- multi2di(tr)
    write.tree(newtr, file = paste("~/yytao/16S/remove1/remove1_binary.nwk.",line,sep=""))
}
close(con)
