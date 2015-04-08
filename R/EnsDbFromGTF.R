####
## function to create a EnsDb object (or rather the SQLite database) from
## a Ensembl GTF file.
## Limitation:
## + There is no way to get the Entrezgene ID from this file.
## + Assuming that the element 2 in a row for a transcript represents its biotype, since
##   there is no explicit key transcript_biotype in element 9.
## + The CDS features in the GTF are somewhat problematic, while we're used to get just the
##   coding start and end for a transcript from the Ensembl perl API, here we get the coding
##   start and end for each exon.
ensDbFromGtf <- function(gtf, verbose=FALSE){
    options(useFancyQuotes=FALSE)
    if(verbose)
        cat("importing gtf file...")
    wanted.features <- c("gene", "transcript", "exon", "CDS")
    GTF <- import(con=gtf, format="gtf", feature.type=wanted.features, asRangedData=FALSE)
    if(verbose)
        cat("done\n")
    ## check what we've got...
    ## all wanted features?
    if(!(wanted.features %in% levels(GTF$type)))
        stop(paste0("One or more required types are not in the gtf file. Need ",
                    paste(wanted.features, collapse=","), " but got only ",
                    paste(wanted.features[wanted.features %in% level(GTF$type)], collapse=","),
                    "."))
    ## transcript biotype?
    if(any(colnames(mcols(GTF)))=="transcript_biotype"){
        txBiotypeCol <- "transcript_biotype"
    }else{
        ## that's a little weird, but it seems that certain gtf files from Ensembl
        ## provide the transcript biotype in the element "source"
        txBiotypeCol <- "source"
    }
    ## create database etc
    dbname <- dummy.sqlite
    con <- dbConnect(dbDriver("SQLite"), dbname=dbname)
    ## ----------------------------
    ##
    ## process genes
    ## we're lacking NCBI Entrezids and also the coord system, but these are not
    ## required columns anyway...
    if(verbose){
        cat("processing genes...")
    }
    ## want to have: gene_id, gene_name, entrezid, gene_biotype, gene_seq_start,
    ##               gene_seq_end, seq_name, seq_strand, seq_coord_system.
    reqCols <- c("gene_id", "gene_name", "gene_biotype")
    if(!any(regCols %in% colnames(mcols(GTF))))
        stop(paste0("One or more required fields are not defined in the GTF! Need ",
                    paste(reqCols, collapse=","), " but got only ",
                    paste(reqCols[reqCols %in% colnames(mcols(GTF))], collapse=","),
                    "."))
    genes <- as.data.frame(GTF[GTF$type == "gene", reqCols])
    genes <- cbind(genes, entrezid=rep(NA, nrow(genes)),
                   seq_coord_system=rep(NA, nrow(genes)))
    colnames(genes)[1:5] <- c("seq_name", "gene_seq_start", "gene_seq_end", "width",
                              "seq_strand")
    ## transforming seq_strand from +/- to +1, -1.
    strand <- rep(0, nrow(genes))
    strand[as.character(genes$seq_strand) == "+"] <- 1
    strand[as.character(genes$seq_strand) == "-"] <- -1
    genes[ , "seq_strand"] <- strand
    ## rearranging data.frame...
    genes <- genes[ , c("gene_id", "gene_name", "entrezid", "gene_biotype",
                        "gene_seq_start", "gene_seq_end", "seq_name",
                        "seq_strand", "seq_coord_system")]
    dbWriteTable(con, name="gene", genes, overwrite=TRUE, row.names=FALSE)
    rm(genes)   ## don't need that anymore...
    if(verbose){
        cat("OK\n")
    }
    ## ----------------------------
    ##
    ## process transcripts
    if(verbose)
        cat("processing transcripts...")
    ## want to have: tx_id, tx_biotype, tx_seq_start, tx_seq_end, tx_cds_seq_start,
    ##               tx_cds_seq_end, gene_id
    reqCols <- c("transcript_id", "gene_id", txBiotypeCol)
    if(!any(regCols %in% colnames(mcols(GTF))))
        stop(paste0("One or more required fields are not defined in the GTF! Need ",
                    paste(reqCols, collapse=","), " but got only ",
                    paste(reqCols[reqCols %in% colnames(mcols(GTF))], collapse=","),
                    "."))
    if(verbose){
        cat("OK\n")
    }
}

.prepareGtf <-

### timing: 178, 0.6, 179
ensDbFromGtf.old <- function(gtf, verbose=FALSE, fast=FALSE){
    options(useFancyQuotes=FALSE)
    if(length(grep(gtf, pattern="gz$"))){
        gtfin <- gzfile(gtf, open="r")
    }else{
        gtfin <- file(gtf, open="r")
    }
    i <- 0
    on.exit(close(gtfin))
    Header <- matrix(ncol=2, nrow=0)
    colnames(Header) <- c("name", "value")
    ## what attributes do we need for our database tables?
    ## GENE: gene_id, gene_name, entrezid (empty), gene_biotype, gene_seq_start,
    ##       gene_seq_end, seq_name, seq_strand, seq_coord_system.
    ## TRANSCRIPT: tx_id, tx_biotype, tx_seq_start, tx_seq_end, tx_cds_seq_start,
    ##             tx_cds_seq_end, gene_id
    ## EXON: exon_id, exon_seq_start, exon_seq_end
    ## T2E: tx_id, exon_id, exon_idx
    ## CHR: seq_name, seq_length, is_circular
    NR <- 1000
    geneBuffer <- data.frame(matrix(ncol=9, nrow=NR, NA))
    colnames(geneBuffer) <- c("gene_id", "gene_name", "entrezid", "gene_biotype",
                              "gene_seq_start", "gene_seq_end", "seq_name",
                              "seq_strand", "seq_coord_system")
    geneIndex <- 0
    txBuffer <- data.frame(matrix(ncol=7, nrow=NR, NA))
    colnames(txBuffer) <- c("tx_id", "tx_biotype", "tx_seq_start", "tx_seq_end",
                            "tx_cds_seq_start", "tx_cds_seq_end", "gene_id")
    txIndex <- 0
    exonBuffer <- data.frame(matrix(ncol=3, nrow=NR, NA))
    colnames(exonBuffer) <- c("exon_id", "exon_seq_start", "exon_seq_end")
    exonIndex <- 0
    t2eBuffer <- data.frame(matrix(ncol=3, nrow=NR, NA))
    colnames(t2eBuffer) <- c("tx_id", "exon_id", "exon_idx")
    t2eIndex <- 0
    exonEnv <- new.env()    ## will use this one to store the names of the exons...
    dbname <- "test.sqlite"
    con <- initializeEnsDb(dbname)
    on.exit(dbDisconnect(con))
    while(length(Line <- readLines(gtfin, n=1)) > 0){
        i <- i+1
        if(verbose){
            if(floor(i/10000) == ceiling(i/10000))
                message(paste0(" ...processed ", i, " lines..."))
        }
        if( i == 30000)
            break
        ## check if we have to save one of the tables...
        if(geneIndex == NR){
            geneIndex <- 0
            dbWriteTable(con, name="gene", value=geneBuffer, append=TRUE)
        }
        if(txIndex == NR){
            txIndex <- 0
            dbWriteTable(con, name="tx", value=txBuffer, append=TRUE)
        }
        if(exonIndex == NR){
            exonIndex <- 0
            dbWriteTable(con, name="exon", value=exonBuffer, append=TRUE)
        }
        if(t2eIndex == NR){
            t2eIndex <- 0
            dbWriteTable(con, name="tx2exon", value=t2eBuffer, append=TRUE)
        }
        ## header.
        if(length(grep(Line, pattern="^#")) > 0 ){
            Header <- rbind(Header, .getHeaderInfo(Line))
        }else{
            ## OK, now processing the tab separated line
            ## 1:seqname, 2:source, 3:feature, 4:start, 5:end, 6:score
            ## 7:strand, 8:frame, 9:attribute (; separated )
            elms <- unlist(strsplit(Line, split="\t"))
            if(elms[3] == "gene"){
                geneIndex <- geneIndex+1
                geneBuffer[geneIndex, ] <- buildGeneList(elms)
            }else if(elms[3] == "transcript"){
                ## if we have the buffer almost full save and reset it, as we're still waiting for
                ## eventual CDS information...
                if(txIndex == (NR-1)){
                    txIndex <- 0
                    dbWriteTable(con, name="tx", value=txBuffer, append=TRUE)
                }
                txIndex <- txIndex + 1
                txBuffer[txIndex, ] <- buildTxList(elms)
            }else if(elms[3] == "exon"){
                ## Note: exons are stored for each transcript in the gtf, while we
                ## only want to have a unique list of exons in the exon table.
                exonAttrs <- processAttrs(elms[9])
                t2eIndex <- t2eIndex + 1
                t2eBuffer[t2eIndex, ] <- list(exonAttrs["transcript_id"], exonAttrs["exon_id"],
                                              as.numeric(exonAttrs["exon_number"]))
                ## check if we do already have the exon, if so, don't add!
                if(!(any(ls(envir=exonEnv) == exonAttrs["exon_id"]))){
                    assign(exonAttrs["exon_id"], 1)
                    exonIndex <- exonIndex + 1
                    exonBuffer[exonIndex, ] <- list(exonAttrs["exon_id"], as.numeric(elms[4]),
                                                    as.numeric(elms[5]))
                }
            }else if(elms[3] == "CDS"){
                ## we have to update the tx_cds_seq_start and tx_cds_seq_end attribute of the
                ## transcript (which we do already have!!)
                ## record the min start and max end for a transcript.
                cdsAttrs <- processAttrs(elms[9])
                idx <- which(txBuffer$tx_id==cdsAttrs["transcript_id"])
                if(length(idx) == 0)
                    stop(paste0("Something went wrong! Got a CDS without a transcript! Transcript ID: ",
                                cdsAttrs["transcript_id"]))
                txBuffer[idx, "tx_cds_seq_start"] <- min(c(as.numeric(txBuffer[idx, "tx_cds_seq_start"]),
                                                           as.numeric(elms[4])), na.rm=TRUE)
                txBuffer[idx, "tx_cds_seq_end"] <- max(c(as.numeric(txBuffer[idx, "tx_cds_seq_end"]),
                                                           as.numeric(elms[5])), na.rm=TRUE)
            }else{
                ##cat("ELSE: ", Line, "\n")
            }
        }
    }
    ## save the buffers...
    dbWriteTable(con, name="gene", value=geneBuffer[1:geneIndex, ], append=TRUE)
    dbWriteTable(con, name="tx", value=txBuffer[1:txIndex, ], append=TRUE)
    dbWriteTable(con, name="exon", value=exonBuffer[1:exonIndex, ], append=TRUE)
    dbWriteTable(con, name="tx2exon", value=t2eBuffer[1:t2eIndex, ], append=TRUE)
    if(verbose){
        cat("...done.")
    }
    return(Header)
}


### timing: 181, 0.6, 182.
buffer.ensDbFromGtf <- function(gtf, verbose=FALSE, fast=FALSE){
    options(useFancyQuotes=FALSE)
    if(length(grep(gtf, pattern="gz$"))){
        gtfin <- gzfile(gtf, open="r")
    }else{
        gtfin <- file(gtf, open="r")
    }
    i <- 0
    on.exit(close(gtfin))
    Header <- matrix(ncol=2, nrow=0)
    colnames(Header) <- c("name", "value")
    ## what attributes do we need for our database tables?
    ## GENE: gene_id, gene_name, entrezid (empty), gene_biotype, gene_seq_start,
    ##       gene_seq_end, seq_name, seq_strand, seq_coord_system.
    ## TRANSCRIPT: tx_id, tx_biotype, tx_seq_start, tx_seq_end, tx_cds_seq_start,
    ##             tx_cds_seq_end, gene_id
    ## EXON: exon_id, exon_seq_start, exon_seq_end
    ## T2E: tx_id, exon_id, exon_idx
    ## CHR: seq_name, seq_length, is_circular
    ## create one temp file per table, open the connection and use writeLines.
    geneFile <- file("ens_gene.txt", open="a")
    txFile <- file("ens_tx.txt", open="a")
    exonFile <- file("ens_exon.txt", open="a")
    t2eFile <- file("ens_tx2exon.txt", open="a")
    closeAll <- function(x){
        close(geneFile)
        close(txFile)
        close(exonFile)
        close(t2eFile)
    }
    on.exit(
        closeAll()
        )
    NR <- 1000
    geneBuffer <- character(NR)
    ## geneBuffer <- matrix(ncol=9, nrow=NR, NA)
    ## colnames(geneBuffer) <- c("gene_id", "gene_name", "entrezid", "gene_biotype",
    ##                           "gene_seq_start", "gene_seq_end", "seq_name",
    ##                           "seq_strand", "seq_coord_system")
    geneIndex <- 0
    txBuffer <- matrix(ncol=7, nrow=NR, NA)
    colnames(txBuffer) <- c("tx_id", "tx_biotype", "tx_seq_start", "tx_seq_end",
                            "tx_cds_seq_start", "tx_cds_seq_end", "gene_id")
    txIndex <- 0
    ## exonBuffer <- matrix(ncol=3, nrow=NR, NA)
    ## colnames(exonBuffer) <- c("exon_id", "exon_seq_start", "exon_seq_end")
    exonBuffer <- character(NR)
    exonIndex <- 0
    ## t2eBuffer <- matrix(ncol=3, nrow=NR, NA)
    ## colnames(t2eBuffer) <- c("tx_id", "exon_id", "exon_idx")
    t2eBuffer <- character(NR)
    t2eIndex <- 0
    exonEnv <- new.env()    ## will use this one to store the names of the exons...
    while(length(Line <- readLines(gtfin, n=1)) > 0){
        i <- i+1
        if(verbose){
            if(floor(i/10000) == ceiling(i/10000))
                message(paste0(" ...processed ", i, " lines..."))
        }
        if( i == 30000)
            break
        ## check if we have to save one of the tables...
        if(geneIndex == NR){
            geneIndex <- 0
            writeLines(geneBuffer, con=geneFile)
        }
        if(txIndex == NR){
            txIndex <- 0
            writeLines(apply(txBuffer, MARGIN=1, paste, collapse="\t"), con=txFile)
        }
        if(exonIndex == NR){
            exonIndex <- 0
            writeLines(exonBuffer, con=exonFile)
        }
        if(t2eIndex == NR){
            t2eIndex <- 0
            writeLines(t2eBuffer, con=t2eFile)
        }
        ## header.
        if(length(grep(Line, pattern="^#")) > 0 ){
            Header <- rbind(Header, .getHeaderInfo(Line))
        }else{
            ## OK, now processing the tab separated line
            ## 1:seqname, 2:source, 3:feature, 4:start, 5:end, 6:score
            ## 7:strand, 8:frame, 9:attribute (; separated )
            elms <- unlist(strsplit(Line, split="\t"))
            if(elms[3] == "gene"){
                geneIndex <- geneIndex+1
                geneBuffer[geneIndex] <- paste(buildGeneRow(elms), collapse="\t")
            }else if(elms[3] == "transcript"){
                ## if we have the buffer almost full save and reset it, as we're still waiting for
                ## eventual CDS information...
                if(txIndex == (NR-1)){
                    writeLines(apply(txBuffer[1:(NR-1), ], MARGIN=1, paste, collapse="\t"), con=txFile)
                    txIndex <- 0
                }
                txIndex <- txIndex + 1
                txBuffer[txIndex, ] <- buildTxRow(elms)
            }else if(elms[3] == "exon"){
                ## Note: exons are stored for each transcript in the gtf, while we
                ## only want to have a unique list of exons in the exon table.
                exonAttrs <- processAttrs(elms[9])
                t2eIndex <- t2eIndex + 1
                t2eBuffer[t2eIndex] <- paste(c(exonAttrs["transcript_id"], exonAttrs["exon_id"],
                                              exonAttrs["exon_number"]), collapse="\t")
                ## check if we do already have the exon, if so, don't add!
                if(!(any(ls(envir=exonEnv) == exonAttrs["exon_id"]))){
                    assign(exonAttrs["exon_id"], 1)
                    exonIndex <- exonIndex + 1
                    exonBuffer[exonIndex] <- paste(c(exonAttrs["exon_id"], as.numeric(elms[4]),
                                                    elms[5]), collapse="\t")
                }
            }else if(elms[3] == "CDS"){
                ## we have to update the tx_cds_seq_start and tx_cds_seq_end attribute of the
                ## transcript (which we do already have!!)
                ## record the min start and max end for a transcript.
                cdsAttrs <- processAttrs(elms[9])
                idx <- which(txBuffer[, "tx_id"]==cdsAttrs["transcript_id"])
                if(length(idx) == 0)
                    stop(paste0("Something went wrong! Got a CDS without a transcript! Transcript ID: ",
                                cdsAttrs["transcript_id"]))
                txBuffer[idx, "tx_cds_seq_start"] <- min(c(as.numeric(txBuffer[idx, "tx_cds_seq_start"]),
                                                           as.numeric(elms[4])), na.rm=TRUE)
                txBuffer[idx, "tx_cds_seq_end"] <- max(c(as.numeric(txBuffer[idx, "tx_cds_seq_end"]),
                                                           as.numeric(elms[5])), na.rm=TRUE)
            }else{
                ##cat("ELSE: ", Line, "\n")
            }
        }
    }
    writeLines(geneBuffer[1:geneIndex ], con=geneFile)
    writeLines(apply(txBuffer[1:txIndex, ], MARGIN=1, paste, collapse="\t"), con=txFile)
    writeLines(exonBuffer[1:exonIndex ], con=exonFile)
    writeLines(t2eBuffer[1:t2eIndex ], con=t2eFile)

    if(verbose){
        cat("...done.")
    }
    ## Still have to a) transform to data.frame, save as table...
    return(Header)
}

### timing: 146, 1, 150.
### full data:
files.ensDbFromGtf <- function(gtf, verbose=FALSE){
    options(useFancyQuotes=FALSE)
    if(length(grep(gtf, pattern="gz$"))){
        gtfin <- gzfile(gtf, open="r")
    }else{
        gtfin <- file(gtf, open="r")
    }
    i <- 0
    Header <- matrix(ncol=2, nrow=0)
    colnames(Header) <- c("name", "value")
    ## what attributes do we need for our database tables?
    ## GENE: gene_id, gene_name, entrezid (empty), gene_biotype, gene_seq_start,
    ##       gene_seq_end, seq_name, seq_strand, seq_coord_system.
    ## TRANSCRIPT: tx_id, tx_biotype, tx_seq_start, tx_seq_end, tx_cds_seq_start,
    ##             tx_cds_seq_end, gene_id
    ## EXON: exon_id, exon_seq_start, exon_seq_end
    ## T2E: tx_id, exon_id, exon_idx
    ## CHR: seq_name, seq_length, is_circular
    ## create one temp file per table, open the connection and use writeLines.
    geneFile <- file("ens_gene.txt", open="a")
    txFile <- file("ens_tx.txt", open="a")
    exonFile <- file("ens_exon.txt", open="a")
    t2eFile <- file("ens_tx2exon.txt", open="a")
    closeAll <- function(x){
        close(geneFile)
        close(txFile)
        close(exonFile)
        close(t2eFile)
        close(gtfin)
    }
    on.exit(
        closeAll()
        )
    NR <- 10
    txBuffer <- matrix(ncol=7, nrow=NR, NA)
    txIdx <- 0
    colnames(txBuffer) <- c("tx_id", "tx_biotype", "tx_seq_start", "tx_seq_end",
                            "tx_cds_seq_start", "tx_cds_seq_end", "gene_id")
    exonEnv <- new.env()    ## will use this one to store the names of the exons...
    while(length(Line <- readLines(gtfin, n=1)) > 0){
        i <- i+1
        if(verbose){
            if(floor(i/10000) == ceiling(i/10000))
                message(paste0(" ...processed ", i, " lines..."))
        }
        ## if( i == 30000)
        ##     break
        ## header.
        if(length(grep(Line, pattern="^#")) > 0 ){
            Header <- rbind(Header, .getHeaderInfo(Line))
        }else{
            ## OK, now processing the tab separated line
            ## 1:seqname, 2:source, 3:feature, 4:start, 5:end, 6:score
            ## 7:strand, 8:frame, 9:attribute (; separated )
            elms <- unlist(strsplit(Line, split="\t"))
            if(elms[3] == "gene"){
                writeLines(paste(buildGeneRow(elms), collapse="\t"), con=geneFile)
            }else if(elms[3] == "transcript"){
                ## if we have the buffer almost full save and reset it, as we're still waiting for
                ## eventual CDS information...
                if(txIdx > 0){
                    writeLines(paste(txBuffer[txIdx, ], collapse="\t"), con=txFile)
                    txIdx <- 0
                }
                txIdx <- txIdx + 1
                txBuffer[txIdx, ] <- buildTxRow(elms)
            }else if(elms[3] == "exon"){
                ## Note: exons are stored for each transcript in the gtf, while we
                ## only want to have a unique list of exons in the exon table.
                exonAttrs <- processAttrs(elms[9])
                writeLines(paste(c(exonAttrs["transcript_id"], exonAttrs["exon_id"],
                                   exonAttrs["exon_number"]), collapse="\t"), con=t2eFile)
                ## check if we do already have the exon, if so, don't add!
                if(!(any(ls(envir=exonEnv) == exonAttrs["exon_id"]))){
                    assign(exonAttrs["exon_id"], 1)
                    writeLines(paste(c(exonAttrs["exon_id"], as.numeric(elms[4]),
                                       elms[5]), collapse="\t"), con=exonFile)
                }
            }else if(elms[3] == "CDS"){
                ## we have to update the tx_cds_seq_start and tx_cds_seq_end attribute of the
                ## transcript (which we do already have!!)
                ## record the min start and max end for a transcript.
                cdsAttrs <- processAttrs(elms[9])
                idx <- which(txBuffer[, "tx_id"]==cdsAttrs["transcript_id"])
                if(length(idx) == 0)
                    stop(paste0("Something went wrong! Got a CDS without a transcript! Transcript ID: ",
                                cdsAttrs["transcript_id"]))
                txBuffer[idx, "tx_cds_seq_start"] <- min(c(as.numeric(txBuffer[idx, "tx_cds_seq_start"]),
                                                           as.numeric(elms[4])), na.rm=TRUE)
                txBuffer[idx, "tx_cds_seq_end"] <- max(c(as.numeric(txBuffer[idx, "tx_cds_seq_end"]),
                                                           as.numeric(elms[5])), na.rm=TRUE)
            }else{
                ##cat("ELSE: ", Line, "\n")
            }
        }
    }
    if(txIdx > 0){
        writeLines(paste(txBuffer[txIdx, ], collapse="\t"), con=txFile)
        txIdx <- 0
    }
    if(verbose){
        cat("...done.")
    }
    ## building metadata table:
    MetaData <- data.frame(matrix(ncol=2, nrow=11))
    colnames(MetaData) <- c("name", "value")
    MetaData[1, ] <- c("Db type", "EnsDb")
    MetaData[2, ] <- c("Type of Gene ID", "Ensembl Gene ID")
    MetaData[3, ] <- c("Supporting package", "ensembldb")
    MetaData[4, ] <- c("Db created by", "ensembldb package from Bioconductor")
    MetaData[5, ] <- c("script_version", "0.0.1")
    MetaData[6, ] <- c("Creation time", date())
    tmp <- unlist(strsplit(gtf, split=.Platform$file.sep, fixed=TRUE))
    splitty <- unlist(strsplit(tmp[length(tmp)], split=".", fixed=TRUE))
    MetaData[7, ] <- c("ensembl_version", splitty[(grep(splitty, pattern="gtf")-1)])
    MetaData[8, ] <- c("ensembl_host", tmp[length(tmp)])
    MetaData[9, ] <- c("Organism", organismFromFilename(gtf) )
    MetaData[10, ] <- c("genome_build", Header[ Header[, "name"]=="genome-build", "value"])
    MetaData[11, ] <- c("DBSCHEMAVERSION", "1.0")
    write.table(MetaData, "ens_metadata.txt", sep="\t", row.names=FALSE)
    ## Still have to a) transform to data.frame, save as table...
    return(Header)
}



initializeEnsDb <- function(dbname){
    con <- dbConnect(dbDriver("SQLite"), dbname=dbname)
    dbSendQuery(con , "CREATE TABLE chromosome ('seq_name' TEXT, 'seq_length' INTEGER, 'is_circular' INTEGER);")
    dbSendQuery(con, "CREATE TABLE gene ('gene_id' TEXT, 'gene_name' TEXT, 'entrezid' TEXT, 'gene_biotype' TEXT, 'gene_seq_start' INTEGER, 'gene_seq_end' INTEGER, 'seq_name' TEXT, 'seq_strand' INTEGER, 'seq_coord_system' TEXT);")
    dbSendQuery(con, "CREATE TABLE tx ('tx_id' TEXT, 'tx_biotype' TEXT, 'tx_seq_start' INTEGER, 'tx_seq_end' INTEGER, 'tx_cds_seq_start' TEXT, 'tx_cds_seq_end' TEXT, 'gene_id' TEXT);")
    dbSendQuery(con, "CREATE TABLE exon ('exon_id' TEXT, 'exon_seq_start' INTEGER, 'exon_seq_end' INTEGER);")
    dbSendQuery(con, "CREATE TABLE tx2exon ('tx_id' TEXT, 'exon_id' TEXT, 'exon_idx' INTEGER);")
    return(con)
}

insertExonIfNotPresent <- function(con, exonId, start, end){
    exCount <- dbGetQuery(con, paste0("select count(*) from tx2exon where exon_id='",
                                      exonId, "'"))
    if(exCount[1,1] == 0){
        dbSendQuery(con, paste0("insert into exon values ('",
                                exonId,
                                "', ", start,
                                ", ", end, ");"))
    }
}

## the elements of a transcript:
## gene_id, transcript_id, gene_name, gene_source, gene_biotype, transcript_name,
## transcript_source, tag.
buildTxQuery <- function(txId, txBiotype, start, end, cdsStart, cdsEnd, geneId){
    query <- paste0("insert into tx values ('", txId,
                    "', '", txBiotype,
                    "', ", start,
                    ", ", end,
                    ", ", cdsStart,
                    ", ", cdsEnd,
                    ", '", geneId,
                    "');")
    ##cat(query, "\n")
    return(query)
}

## process the elements and save the gene
buildGeneQuery <- function(elms){
    attrs <- processAttrs(elms[9])
    query <-  paste0("insert into gene values ('",
                     attrs["gene_id"],
                     "', '", attrs["gene_name"],
                     "', ", "NULL",
                     ", '", attrs["gene_biotype"],
                     "', ", elms[4],
                     ", ", elms[5],
                     ", '", elms[1],
                     "', ", strand2Int(elms[7]),
                     ", ", "NULL",
                     ");")
    return(query)
}

## transform a line from the GTF to a character vector.
buildGeneRow <- function(elms){
    attrs <- processAttrs(elms[9])
    return(c(attrs["gene_id"], attrs["gene_name"], entrezid=NA,
             attrs["gene_biotype"], gene_seq_start=elms[4],
             gene_seq_end=elms[5], seq_name=elms[1], seq_strand=strand2Int(elms[7]),
             seq_coord_system=NA))
}
buildGeneList <- function(elms){
    attrs <- processAttrs(elms[9])
    return(list(gene_id=attrs["gene_id"], gene_name=attrs["gene_name"], entrezid=NA,
                gene_biotype=attrs["gene_biotype"], gene_seq_start=as.numeric(elms[4]),
                gene_seq_end=as.numeric(elms[5]), seq_name=elms[1], seq_strand=strand2Int(elms[7]),
                seq_coord_system=NA))
}
buildTxRow <- function(elms){
    ## Note: there is no transcript biotype, but it seems that elms[2]
    ## is the transcript biotype...
    attrs <- processAttrs(elms[9])
    return(c(tx_id=unname(attrs["transcript_id"]), tx_biotype=elms[2], tx_seq_start=elms[4],
             tx_seq_end=elms[5], tx_cds_seq_start=NA, tx_cds_seq_end=NA,
             attrs["gene_id"]))
}
buildTxList <- function(elms){
    ## Note: there is no transcript biotype, but it seems that elms[2]
    ## is the transcript biotype...
    attrs <- processAttrs(elms[9])
    return(list(tx_id=attrs["transcript_id"], tx_biotype=elms[2],
                tx_seq_start=as.numeric(elms[4]), tx_seq_end=as.numeric(elms[5]),
                tx_cds_seq_start=NA, tx_cds_seq_end=NA,
                gene_id=attrs["gene_id"]))
}

strand2Int <- function(x){
    if(x=="+"){
        return(1)
    }else if(x=="-"){
        return(-1)
    }
}

processAttrs <- function(x){
    ##if(remQuotes)
    x <- gsub(x, pattern="\"", replacement="")
    pairs <- unlist(strsplit(x, split="; "))
    ## last one can still contain a ;
    pairs[length(pairs)] <- gsub(pairs[length(pairs)], pattern=";", replacement="")
    procPairs <- strsplit(pairs, split=" ", fixed=TRUE)
    attrs <- unlist(lapply(procPairs, function(x)x[2]))
    names(attrs) <- unlist(lapply(procPairs, function(x)x[1]))
    return(attrs)
}

.getHeaderInfo <- function(x){
    x <- sub(x, pattern="#", replacement="")
    x <- sub(x, pattern="!", replacement="")
    return(unlist(strsplit(x, split=" ", fixed=TRUE)))
}


## compare the contents of the EnsDb sqlite database generated from a GTF (file name submitted
## with x ) with the one provided by package "lib".
evaluateGtfEnsDb <- function(x, lib){
}

organismFromFilename <- function(x){
    tmp <- unlist(strsplit(x, split=.Platform$file.sep, fixed=TRUE))
    splitty <- unlist(strsplit(tmp[length(tmp)], split=".", fixed=TRUE))
    return(splitty[1])
}



