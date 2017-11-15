library(RPostgreSQL)
library(dplyr)
library(GenomicRanges)
library(doParallel)
library(dbplyr)
library(BiocParallel)
library(fst)

source("/scratch/github/BDDS/footprints/testdb/src/dbFunctions.R")


constructLymphoblastDataset <- function(){
    # Load the chipseq data locally
    load(system.file(package = "FPML", "extdata/Rdata_files/chipSeqData.Rdata"))
    
    # Load the TF-Motif mapping data
    load(system.file(package = "FPML", "extdata/Tfmotifmap.Rdata"))
    
    # Grab all motifs for the TFs we have
    allmots <- c()
    
    for (TFname in names(TFs.to.motifs)) {
        allmots  <-  c(allmots, TFs.to.motifs[[TFname]])
    }
    # length(unique(allmots))
    
    # Run it in parallel
    sorted.TF.names <- sort(names(TFs.to.motifs))
    
    BiocParallel::register(BiocParallel::MulticoreParam(workers = 62,
                                                        stop.on.error = FALSE,
                                                        log = TRUE),
                           default = TRUE)
    all.TF.df <- BiocParallel::bplapply(sorted.TF.names, createTfDf, verbose = TRUE)
    all.TF.df <- dplyr::bind_rows(all.TF.df)
    
    # If you want to take just distinct rows, do it here
    if(distinctFlag){
        all.TF.df <- all.TF.df %>% dplyr::distinct()
    }
    
    return(all.TF.df)
    
} # createLymphoblastDataset
#----------------------------------------------------------------------------------------------------
sampleTfDataset <- function(all.TF.df, sampleSize){
    
    # Create a sample of the correct size
    sample.df <- sample(1:nrow(all.TF.df), size = sampleSize)
    
    # Return the sampled dataframe
    return(all.TF.df[sample.df,])
    
} #sampleTfDataset
#----------------------------------------------------------------------------------------------------
# Define the function to pull motifs for a TF
createTfDf <- function(TF, verbose = FALSE){
    
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
