#----------------------------------------------------------------------------------------------------
#' Merge Footprints for Lymphoblast Data
#'
#' Merge the HINT and Wellington footprints into the dataframe of FIMO motifs and ChIPSeq Data for
#' all chromosomes
#'
#' @param fimo.df A data frame containing FIMO motifs and ChIPSeq data, generally constructed
#' using the \code{constructLymphoblastDataset.R} function
#' @param seedNum A number (either 16 or 20) denoting which seed size to use for the merge
#' @param host A string denoting the database host (default = "localhost")
#' @param verbose A Boolean value that determines whether the function should print output during
#' the merge (default = TRUE)
#'
#' @return The supplied FIMO/ChIPSeq data frame with added columns for HINT count, HINT maximum
#' score, Wellington count, and Wellington minimum score
#'
#' @export

mergeLymphoblastFootprints <- function(fimo.df, seedNum, host = "localhost", verbose = TRUE){

    # Create database name strings
    db.hint <- paste0("lymphoblast_hint_", seedNum)
    db.well <- paste0("lymphoblast_wellington_", seedNum)
    
    # Read data from lymphoblast
    db_lymph_hint <- DBI::dbConnect(drv = RPostgreSQL::PostgreSQL(),
                                    user="trena",
                                    password="trena",
                                    dbname = db.hint,
                                    host = host)
    db_lymph_well <- DBI::dbConnect(drv = RPostgreSQL::PostgreSQL(),
                                    user="trena",
                                    password="trena",
                                    dbname = db.well,
                                    host = host)
    hint_regions <- dplyr::tbl(db_lymph_hint, "regions")
    hint_hits    <- dplyr::tbl(db_lymph_hint, "hits")
    well_regions <- dplyr::tbl(db_lymph_well, "regions")
    well_hits    <- dplyr::tbl(db_lymph_well, "hits")

    # Run the big function using a for loop
    chromosomes <- as.character(c(1:22, "X","Y"))

    big_list <- list()
    counter <- 1
    now.time <- Sys.time()

    for (chr_str in chromosomes) {
        # Print if it's verbose
        if(verbose){
            message(paste("working on chromosome", chr_str))
        }
        
        big_list[[counter]] <- mergeFootprintsOneChrom(chr_str,
                                                       fimo.df,
                                                       hint_regions,
                                                       hint_hits,
                                                       well_regions,
                                                       well_hits)
        counter <- counter + 1

        if(verbose){
            message(paste("Time elapsed:", Sys.time() - now.time))
        }
    }
    
    # Combine them
    merged.df <- dplyr::bind_rows(big_list)
    rm(big_list)
    
    # Return them
    return(merged.df)
    
} # mergeLymphoblastFootprints
#----------------------------------------------------------------------------------------------------
#' Merge Footprints for 1 Chromosome
#'
#' Merge HINT and Wellington footprints into the dataframe of FIMO motifs and ChipSeq for
#' 1 chromosome
#'
#' @param chrom_str A string representing one chromosome; must be one of c(1:22, "X", "Y")
#' @param fimo_tbl A dataframe of FIMO motifs and ChIPSeq data
#' @param hint_regions_tbl A database table of HINT regions
#' @param hint_hits_tbl A database table of HINT hits
#' @param well_regions_tbl A database table of Wellington regions
#' @param well_hits_tbl A database table of Wellington hits
#'
#' @return The original FIMO dataframe with added columns for footprint counts and the maximum or
#' minimum score
#' 
#' @export

mergeFootprintsOneChrom <- function(chrom_str,
                                    fimo_tbl,
                                    hint_regions_tbl,
                                    hint_hits_tbl,
                                    well_regions_tbl,
                                    well_hits_tbl
                                    ) {

    # some tables use chr22 and some just use 22
    chrom_long_str = paste("chr",chrom_str, sep="")

    # select one chromosome from my data table
    fimo_tbl %>%
        dplyr::filter(chrom==chrom_str) %>%
        dplyr::select(-empty) ->
        chrom_all_tf_df

        # select one chromosome from hint
    hint_regions_tbl %>%
        dplyr::filter(chrom==chrom_long_str) %>%
        dplyr::left_join(hint_hits_tbl, by="loc") %>%
        dplyr::select(start, endpos, strand, name, h_score = score1) %>%
        dplyr::group_by(start, endpos, name, strand) %>%
        dplyr::mutate(h_count = n(), h_max_score = max(h_score)) %>%
        dplyr::select(-h_score) %>%
        tibble::as_tibble() ->
        chrom_hint_all_tbl

    # Grab only distinct ones
    chrom_hint_all_tbl %>%
        dplyr::distinct(start, endpos, name, strand, .keep_all = TRUE) ->
        chrom_hint_unique_tbl
    rm(chrom_hint_all_tbl)

    # select one chromosome from wellington
    well_regions_tbl %>%
        dplyr::filter(chrom==chrom_long_str) %>%
        dplyr::left_join(well_hits_tbl, by="loc") %>%
        dplyr::select(start, endpos, strand, name, w_score = score1) %>%
        dplyr::group_by(start, endpos, name, strand) %>%
        dplyr::mutate(w_count = n(), w_min_score = min(w_score)) %>%
        dplyr::select(-w_score) %>%
        tibble::as_tibble() ->
        chrom_well_all_tbl

    # keep only min wellington score but count total nontrivial scores
    chrom_well_all_tbl %>%
        dplyr::distinct(start, endpos, name, strand, .keep_all = TRUE) ->
        chrom_well_unique_tbl
    rm(chrom_well_all_tbl)
    
    # merge hint and wellington into my table
    chrom_all_tf_df %>%
        dplyr::left_join(chrom_hint_unique_tbl,
                         by=c("start", "endpos", "strand", "motifname"="name")) %>%
        dplyr::left_join(chrom_well_unique_tbl,
                         by=c("start", "endpos", "strand", "motifname"="name")) %>%
        tidyr::replace_na(list(h_count=0, w_count=0, h_max_score=0, w_min_score=0)) ->
        chrom_all_tf_df_merged
    
    return(chrom_all_tf_df_merged)

}# mergeFootprintsOneChrom
#----------------------------------------------------------------------------------------------------
