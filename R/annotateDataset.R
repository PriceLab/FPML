#----------------------------------------------------------------------------------------------------
#' Annotate the Footprint Data
#'
#' Annotate the footprint and FIMO dataframe by adding TF class, GC content, and TSS distance
#'
#' @param fp.df The data frame of FIMO motifs and footprints
#'
#' @return The supplied data frame with added columns for the one-hot coded TF class data, GC
#' content, and asinh-transformed distance from the nearest transcription start site
#'
#' @export

annotateFootprintData <- function(fp.df){

    # Get the motif data
    motif_class_hot <- createMotifClassMap(fp.df)
    
    # Add the motif class to the dataframe
    annotated.df <- dplyr::left_join(fp.df,
                              motif_class_hot)

    # Add GC Content
    annotated.df <- annotated.df %>% 
        dplyr::mutate("gc_content" = getGCContent(start,endpos,chrom))

    # Add TSS distance
    annotated.df <- addTSSDistance(annotated.df)
    
    # Change counts to fracs
    annotated.df %>%
        dplyr::mutate(h_frac = h_count/max(h_count)) %>%
        dplyr::mutate(w_frac = w_count/max(w_count)) %>%
        dplyr::select(-one_of("h_count","w_count")) %>%
        dplyr::select(motifname:w_min_score,
                      h_frac, w_frac, gc_content,
                      asinh_tss_dist,
                      dplyr::everything()) ->
        annotated.df

    return(annotated.df)

} # annotateFootprintData
#----------------------------------------------------------------------------------------------------
#' Create Motif Class Map
#'
#' Create the one-hot form of the JASPAR motifs and their classes
#'
#' @param fp.df The data frame of FIMO motifs and footprints. It is used for filtering out only those
#' motifs that appear in the data frame.
#'
#' @return The one-hot data frame of motifs and their membership in different classes
#'
#' @export

createMotifClassMap <- function(fp.df){

    # Load data on motif mapping
    load(system.file(package="FPML", "extdata/Tfmotifmap.Rdata"))
    load(system.file(package="FPML", "extdata/motif_class_pairs.Rdata"))

    # Make Jaspar translation table
    jaspar.motifs <- MotifDb::subset(MotifDb::MotifDb, dataSource == "jaspar2016")
    jaspar.df <- dplyr::data_frame(Long.Name = names(jaspar.motifs),
                                   Short.Name = trimws(S4Vectors::values(jaspar.motifs)$providerName)
                                   )

    # Add long motif names to the Jaspar data
    fixed.motif.class <- motif.class %>%
        dplyr::left_join(jaspar.df, by = c("motif" = "Short.Name")) %>%
        dplyr::select("motifname" = "Long.Name", class)

    # Filter the motif class using the dataframe
    filtered.motif.class <- dplyr::semi_join(fixed.motif.class,
                                             fp.df,
                                             by = "motifname")    

    # Make it into a one-hot form
    filtered.motif.class %>%
        # clean up and subset to only relevant motifs
        dplyr::mutate_all(stringr::str_trim) %>%
        # fix double classes
        dplyr::mutate(class = stringr::str_split(class, "::")) %>%
        tidyr::unnest(class) %>%
        # create one-hot(ish, some double matches) version
        dplyr::mutate(dummy_yesno = 1) %>%
        dplyr::distinct() %>%
        tidyr::spread(class, dummy_yesno, fill = 0) ->
        motif_class_hot

    return(motif_class_hot)

} # createMotifClassMap
#----------------------------------------------------------------------------------------------------
#' Add GC content
#'
#' Add GC content to a row of the data frame containing FIMO motifs and footprint data
#'
#' @param start_col The name of the column containing the starting location of each motif
#' @param end_col The name of the column containing the ending location of each motif
#' @param chrom_col The name of the column containing the chromosome for each motif
#' @param shoulder An integer indicating the distance on either side of the center of the sequence
#' to use as a window (default = 100)
#'
#' @return The GC content of the supplied region
#'
#' @export

getGCContent <- function(start_col, end_col, chrom_col, shoulder=100) {

    window_center <- round((start_col + end_col)/2)
    windows <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                                  paste0("chr",chrom_col),
                                  window_center-shoulder,
                                  window_center+shoulder)

    alph_freq <- Biostrings::alphabetFrequency(windows)
    gc_content <- rowSums(alph_freq[,c("C","G")])/(2*shoulder+1)

    return(gc_content)

} # getGCContent
#----------------------------------------------------------------------------------------------------
#' Add Transcription Start Site Distance
#'
#' Add transcription start site distance to a dataframe of FIMO motifs and footprints data
#'
#' @param annotated.df A dataframe containing FIMO motifs, ChIPSeq data, and footprints data
#' @param host A string indicating the location of the hg38 database (default = "localhost")
#'
#' @return The dataframe with an added column containing the asinh-transformed distance from each
#' motif to the nearest transcription start site
#'
#' @export

addTSSDistance <- function(annotated.df, host = "localhost"){

    db_gtf <- DBI::dbConnect(RPostgreSQL::PostgreSQL(),
                             user = "trena",
                             password = "trena",
                             dbname = "hg38",
                             host = host)
    query <- "select * from gtf where moleculetype='gene' and gene_biotype='protein_coding'"
    tss_raw_table <- DBI::dbGetQuery(db_gtf,
                                     query)[, c("chr", "gene_name", "start", "endpos","strand")]
    
    tss_raw_table %>%
        dplyr::mutate(ref = ifelse(strand == '+', start, endpos)) %>%
        dplyr::select("chrom" = "chr", "ts_start" = "ref") %>%
        dplyr::filter(chrom != 'chrMT') %>%
        dplyr::mutate(chrom=stringr::str_sub(chrom,  start = 4)) ->
        tss_tbl

    motif_gr <- GenomicRanges::makeGRangesFromDataFrame(annotated.df,
                                                        start.field="start",
                                                        end.field="endpos")
    tss_gr <- GenomicRanges::makeGRangesFromDataFrame(annotated.df,
                                                      start.field="ts_start",
                                                      end.field="ts_start")
    dist_to_nearest_tss <- GenomicRanges::distanceToNearest(motif_gr,
                                                            tss_gr,
                                                            select="arbitrary")
    tss_dists <- mcols(dist_to_nearest_tss)[,1]

    annotated.df %>%
        dplyr::mutate(asinh_tss_dist = asinh(tss_dists)) ->
        df.with.tss.dist

    return(df.with.tss.dist)

} # addTSSDistance


   




