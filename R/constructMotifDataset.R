#----------------------------------------------------------------------------------------------------
#' Construct the Lymphoblast Datast
#'
#' Using the included ChIPSeq data, construct the lymphoblast dataset 
#'
#' @param distinctFlag A Boolean flag indicating whether or not to take just the unique rows
#' of the data frame (default = TRUE)
#' @param sampleSize An integer value indicating how many rows of the data frame to take. If this
#' argument is not specified or make to be NULL, the entire data frame will be returned
#' (default = NULL)
#' @param isTest A Boolean flag indicating whether this script is being run as a test. If TRUE,
#' it will only use 2 TFs instead of the full set to save on time (default = FALSE)
#' @return The complete motif/ChIPSeq dataset for lymphoblast
#'
#' @export

constructLymphoblastDataset <- function(distinctFlag = TRUE, sampleSize = NULL,
                                        isTest = FALSE){
    
    # Read the ChIPseq data from the local database
    db.chipseq <- DBI::dbConnect(drv=RPostgreSQL::PostgreSQL(),
                                 user = "trena",
                                 password = "trena",
                                 dbname = "chipseq",
                                 host= "localhost")                           

    # Grab the hits table as is
    chipseq.hits <- DBI::dbGetQuery(db.chipseq, "select * from hits")
    chipseq.hits <- tibble::as_tibble(chipseq.hits)

    # Grab the regions table and modify the chrom column
    chipseq.regions <- DBI::dbGetQuery(db.chipseq, "select * from regions")
    chipseq.regions <- tibble::as_tibble(chipseq.regions)
    chr.list <- chipseq.regions$chrom
    cutoff <- nchar("chr")+1
    no.chr.list <- substring(chr.list,cutoff)
    chipseq.regions$chrom <- no.chr.list

    # Create the TF-Motif mapping using the correct function
    TFs.to.motifs <- mapTFsToMotifs(chipseq.hits)
        
    # Run it in parallel
    sorted.TF.names <- sort(names(TFs.to.motifs))

    # Cut it down to 2 TFs with 1 Motif each if it's just a test
    if (isTest){
        sorted.TF.names <- c("BCL3", "WRNIP1")
    }    
    
    BiocParallel::register(BiocParallel::MulticoreParam(workers = 62,
                                                        stop.on.error = FALSE,
                                                        log = TRUE),
                           default = TRUE)
    all.TF.df <- BiocParallel::bplapply(sorted.TF.names, createTfDf,
                                        TFs.to.motifs = TFs.to.motifs,
                                        chipseq.regions = chipseq.regions,
                                        chipseq.hits = chipseq.hits,
                                        verbose = TRUE)
    all.TF.df <- dplyr::bind_rows(all.TF.df)
    
    # If you want to take just distinct rows, do it here
    if(distinctFlag){
        all.TF.df <- all.TF.df %>% dplyr::distinct()
    }

    if(!is.null(sampleSize)){

        # Take a sample of the correct size
        all.TF.df <- sampleTfDataset(all.TF.df, sampleSize)
    }    
    
    return(all.TF.df)
    
} # createLymphoblastDataset
#----------------------------------------------------------------------------------------------------
#' Create a TF-Motif Map in List Form
#'
#' Given a set of validation data (generally ChIPseq data), return a mapping of all TFs in the data
#' to their JASPAR motifs. This will only return TFs for which there is a mapping to a JASPAR motif
#' via our TOMTOM motif matching. If there is no corresponding JASPAR motif, the TF will be dropped.
#'
#' @param chipseq.hits A table of chipseq hits
#'
#' @return A list where each name is a TF from the supplied validation set and the each entry is
#' a character vector of motifs that map to that TF
#'
#' @export

mapTFsToMotifs <- function(chipseq.hits){

    # Load our included TF mapping
    TF.motif.pairs <- readRDS(system.file(package="FPML", "extdata/2017_08_23_Motif_TF_Map.RDS"))
    names(TF.motif.pairs) <- c("motif", "tfs")

    # Grab the unique set of chipseq TFs
    unique.cs.tfs <- unique(chipseq.hits$name)
    cs.tf.matches <- unique.cs.tfs[unique.cs.tfs %in% TF.motif.pairs$tfs]

    # Grab the correct motifs for each match
    getMotifsForTF <- function(TF, tf.pairs){
        TF.df <- dplyr::filter(tf.pairs, tfs %in% TF)
        return(TF.df$motif)
    }
    
    TF.motif.map <- lapply(cs.tf.matches, getMotifsForTF, tf.pairs = TF.motif.pairs)
    names(TF.motif.map) <- cs.tf.matches

    return(TF.motif.map)
    
} # mapTFsToMotifs
#----------------------------------------------------------------------------------------------------
#' Sample the FIMO/ChIPSeq Dataset
#'
#' Take a random sample of some size from the FIMO/ChIPSeq dataset
#'
#' @param all.TF.df The FIMO/ChIPSeq dataset, a data frame
#' @param sampleSize An integer indicating the sample size to be taken
#'
#' @return A subset of the supplied dataframe, with row number equal to sampleSize
#'
#' @export

sampleTfDataset <- function(all.TF.df, sampleSize){
    
    # Create a sample of the correct size
    sample.df <- sample(1:nrow(all.TF.df), size = sampleSize)
    
    # Return the sampled dataframe
    return(all.TF.df[sample.df,])
    
} #sampleTfDataset
#----------------------------------------------------------------------------------------------------
#' Pull all FIMO Motifs for a Given TF
#'
#' Given a particular TF, pull all FIMO data for the motifs that map to that TF
#'
#' @param TF A transcription factor of interest
#' @param TFs.to.motifs A list where each entry is a character vector of motifs and the name of each entry
#' is a transcription factor
#' @param chipseq.regions A table of ChIPseq regions
#' @param chipseq.hits A table of ChIPseq hits
#' @param verbose A Boolean flag indicating whether the function should print output (default = FALSE)
#'
#' @return The data frame containing all FIMO motifs that correspond to the supplied TF
#'
#' @export

createTfDf <- function(TF, TFs.to.motifs,
                       chipseq.regions,
                       chipseq.hits,
                       verbose = FALSE){

    # Make the database connection
    db.fimo.dplyr <- DBI::dbConnect(drv = RPostgreSQL::PostgreSQL(),
                                    user = "trena",
                                    password = "trena",
                                    dbname = "fimo",
                                    host = "localhost")
    tbl.fimo.dplyr <- dplyr::tbl(db.fimo.dplyr, "fimo_hg38")
    
    # Grab all hits for a TF, then grab the regions of those hits
    chipseq.hits.TF <- subset(chipseq.hits, name == TF)
    locs.TF <- chipseq.hits.TF$loc
    chipseq.regions.TF <- subset(chipseq.regions, loc %in% locs.TF)
    
    # this is the slow step -- doing SQL queries on tbl.fimo.dplyr = call to whole fimo database
    # need branch since %in% conversion to SQL doesn't work on length == 1
    ## Basically: we find all instances of the TF's motif(s) in FIMO
    if (length(TFs.to.motifs[[TF]]) > 1 ) {
        fimo.motifs.for.TF <- tbl.fimo.dplyr %>%
            dplyr::filter(motifname %in% TFs.to.motifs[[TF]]) %>%
            tibble::as_tibble()
    } else {
        fimo.motifs.for.TF <- tbl.fimo.dplyr %>%
            dplyr::filter(motifname == TFs.to.motifs[[TF]]) %>%
            tibble::as_tibble()
    }
    # Print if requested
    if (verbose == TRUE) {
        if (length(TFs.to.motifs[[TF]])==1) {
            message(paste(TF, "- querying fimo database for", length(TFs.to.motifs[[TF]]), "motif"))
        } else {
            message(paste(TF, "- querying fimo database for", length(TFs.to.motifs[[TF]]), "motifs"))
        }
    }
    
    # find intersect using fast genomic ranges data structure
    # We make GR objects for the fimo motifs we just found and for the chipseq regions,
    # then find their overlaps
    browser()
    gr.fimo.TF <- with(fimo.motifs.for.TF,
                       GenomicRanges::GRanges(chrom,
                                              IRanges::IRanges(start=start,
                                                               end=endpos)))
    gr.chipseq.TF <- with(chipseq.regions.TF,
                          GenomicRanges::GRanges(chrom,
                                                 IRanges::IRanges(start=start,
                                                                  end=endpos)))
    overlaps.gr.TF <- GenomicRanges::findOverlaps(gr.chipseq.TF, gr.fimo.TF, type="any")
    overlaps.TF <- tibble::as_tibble(overlaps.gr.TF)
    
    # row numbers in fimo.motifs.for.TF where motifs overlap with chipseq peaks
    positive.fimo.examples.rows.TF <- unique(overlaps.TF$subjectHits)
    positive.fimo.examples.TF.df <- fimo.motifs.for.TF[positive.fimo.examples.rows.TF,]
    
    # Simply take the other rows as the negative fimo examples
    negative.fimo.examples.TF.df <- dplyr::setdiff(fimo.motifs.for.TF, positive.fimo.examples.TF.df)    
        
    # annotate and collect all samples
    positive.fimo.examples.TF.df <- tibble::as_tibble(dplyr::bind_colscbind(positive.fimo.examples.TF.df,
                                                                            "cs_hit"=1))
    negative.fimo.examples.TF.df <- tibble::as_tibble(dplyr::bind_cols(negative.fimo.examples.TF.df,
                                                                       "cs_hit"=0))
    all.fimo.examples.TF.df <- tibble::as_tibble(dplyr::bind_rows(positive.fimo.examples.TF.df,
                                                                  negative.fimo.examples.TF.df))

    return(all.fimo.examples.TF.df)

} # createTfDf
#----------------------------------------------------------------------------------------------------
