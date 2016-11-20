#!/usr/bin/env Rscript

RAP.MSU.refseq <- read.csv(file = "RAP.MSU.refseq.csv")
primerSummary <- read.csv(file = 'primerSummary.csv')

mergedOutput <- merge(primerSummary, RAP.MSU.refseq[,c('MSU.ID', 'symbol')])
mergedOutput$PrimerF.length = nchar(as.character(mergedOutput$PrimerF))
mergedOutput$PrimerR.length = nchar(as.character(mergedOutput$PrimerR))
mergedOutput <- mergedOutput[,c('MSU.ID', 'RefSeqID', 'symbol', 'Status', 'PrimerF', 'PrimerF.length', 'PrimerF.TM', 'PrimerR', 'PrimerR.length', 'PrimerR.TM', 'ProductSize', 'IntronSize')]
mergedOutput[mergedOutput == 0]  <- NA
mergedOutput[mergedOutput == '']  <- NA
write.csv(mergedOutput, file = "primerSummary.csv", na = '', row.names = FALSE, quote = FALSE)
