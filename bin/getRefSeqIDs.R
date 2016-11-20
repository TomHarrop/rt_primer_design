#!/usr/bin/env Rscript

PrintDate <- function(){
  format(Sys.time(), "[%a %b %d %H:%M:%S %Y] ")
}

# sandbox

#genes <- data.frame('MSU.ID' = c('LOC_Os01g09590', 'LOC_Os01g50700', 'LOC_Os03g07570', 'LOC_Os03g29190', 'LOC_Os03g42760', 'LOC_Os04g50920', 'LOC_Os04g53210', 'LOC_Os05g06560', 'LOC_Os05g11414', 'LOC_Os07g32570', 'LOC_Os08g07480', 'LOC_Os08g39694', 'LOC_Os08g41950', 'LOC_Os09g26900', 'LOC_Os10g21784', 'LOC_Os10g35660'))
#inputFile <- 'sandbox'

## get file from command line input (this will be sent by python)

args <- commandArgs(TRUE)
input.file <- args[1]
out.dir <- args[2]
#input.file <- "data/primers.txt"
genes <- read.table(input.file, col.names = 'MSU.ID', stringsAsFactors = FALSE)
#genes <- genes[1:2001,,drop = FALSE]
number.of.records <- dim(genes)[1]
genes <- unique(genes)

cat(PrintDate(), number.of.records, ' records in ', input.file,' including ',
    number.of.records - dim(genes)[1], ' duplicates.\nFinding Rap-DB IDs for ',
    dim(genes)[1], ' records.\nThis could take a while.\n', sep = '')

RAP.MSU.refseq <- oryzr::LocToRefSeq(genes$MSU.ID, useBiomart = TRUE)

if (dim(RAP.MSU.refseq$notMatched)[1] > 0){
  write.csv(RAP.MSU.refseq$notMatched,
            file = paste0(out.dir, "/excludedRecords.csv"), na = '',
            row.names = FALSE, quote = FALSE)
}

write.csv(RAP.MSU.refseq$matched,
          file = paste0(out.dir, "/RAP.MSU.refseq.csv"), na = '',
          row.names = FALSE, quote = FALSE)
