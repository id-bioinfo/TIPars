#!/usr/bin/env Rscript


if (is.null(tryCatch(packageVersion("optparse"), error=function(e) NULL))) {
    stop("'optparse' should be installed ...")
}

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("-t", "--tree"), help="tree file, in Newick format"),
    make_option(c("-s", "--sequence"), help="fasta file contains aligned taxa sequences"),
    make_option(c("-a", "--ancseq"), help="fasta file contains aligned ancestral sequences"),
    make_option(c("-q", "--query"), help="fasta file contains 1 query seqence"),
    make_option(c("-m", "--model"), default="LE", help="one of 'LE' (by default), 'JC69' or 'K2P'"),
    make_option(c("-g", "--gap"), default="ignore", help="one of 'ignore', 'inner' or 'all'"),
    make_option(c("-o", "--output"), default="TIPars_output.tree", help="output tree file name"),
    make_option(c("-p", "--type"), default="insertion", help="one of 'insertion' or 'placement'")
    )

opt <- parse_args(OptionParser(option_list=option_list))

## check parameters

if (is.null(opt$tree))     stop("tree file must be provided...")
if (is.null(opt$sequence)) stop("aligned taxa sequence (fasta) file must be provided...")
if (is.null(opt$ancseq))   stop("aligned ancestral sequence (fasta) file must be provided...")
if (is.null(opt$query))    stop("query seqence (fasta) file must be provided...")


cmd_args <- commandArgs(trailingOnly = FALSE)
exec_script <- cmd_args[which(grepl("--file", cmd_args))[1]]
exec_script <- sub("--file=", "", exec_script)
exec_dir <- dirname(exec_script)

TIPars <- paste0("java -jar ", exec_dir, "/", "TIPars.jar")
cmd <- paste(TIPars, opt$tree, opt$sequence, opt$ancseq, opt$query, opt$model, opt$gap, opt$output, opt$type, "0")

cat("\n\n")
cat(format(Sys.time(), "%Y-%m-%d %X"), "\n")
cat(cmd, "\n")
system(cmd)



