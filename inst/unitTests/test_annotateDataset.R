library(FPML)
library(RUnit)
#----------------------------------------------------------------------------------------------------
# Load the test data (merged.df)
load(system.file(package="FPML", "extdata/annotateTestInput.Rdata"))

# Get the chipseq hit data locally
db.chipseq <- DBI::dbConnect(drv=RPostgreSQL::PostgreSQL(),
                             user = "trena",
                             password = "trena",
                             dbname = "chipseq",
                             host= "localhost")                           

chipseq.hits <- DBI::dbGetQuery(db.chipseq, "select * from hits")
chipseq.hits <- tibble::as_tibble(chipseq.hits)

runTests <- function(){

    test_annotateFootprintData()
    test_createMotifClassMap()
    test_getGCContent()
    test_addTSSDistance()    

    } # runTests

#----------------------------------------------------------------------------------------------------
# Test single TF function

test_annotateFootprintData <- function(){

    message("---test_annotateFootprintData")
    
    # Get the results of annotating
    annotated.df <- annotateFootprintData(merged.df, chipseq.hits)

    # Check that the dimensions are correct
    checkEquals(ncol(annotated.df), 21)
    checkEquals(nrow(annotated.df), 7)

    # Check that we added 7 columns
    checkEquals(ncol(merged.df) + 7, ncol(annotated.df))

    # Check that column names are correct
    expected.names = sort(c("motifname",
                            "chrom",
                            "start",                                
                            "endpos",                               
                            "strand",                               
                            "motifscore",                           
                            "pval",
                            "sequence",
                            "loc",
                            "cs_hit",
                            "h_max_score",
                            "w_min_score",
                            "h_frac",                               
                            "w_frac",                               
                            "gc_content",
                            "asinh_tss_dist",
                            "Basic helix-loop-helix factors (bHLH)",
                            "C2H2 zinc finger factors",    
                            "Fork head / winged helix factors",
                            "Homeo domain factors",                 
                            "Zinc-coordinating"))
    checkEquals(sort(names(annotated.df)), expected.names)
    
} # test_annotateFootprintData
#----------------------------------------------------------------------------------------------------
# Test Motif Class Mapping

test_createMotifClassMap <- function(){

    message("---test_createMotifClassMap")

    # Grab the proper mapping
    TFs.to.motifs <- mapTFsToMotifs(chipseq.hits)
    
    # Create the map for my merged data; it's in one-hot form
    hot.map <- createMotifClassMap(merged.df, TFs.to.motifs)

    # The mapping should have the same number of rows as the input
    checkEquals(nrow(hot.map), nrow(merged.df), 7)
    checkEquals(ncol(hot.map), 6)

    # Check that the 6 columns are the ones we expect
    expected.names <- sort(c("motifname",
                             "Basic helix-loop-helix factors (bHLH)",
                             "C2H2 zinc finger factors",
                             "Fork head / winged helix factors",
                             "Homeo domain factors",
                             "Zinc-coordinating"))
    checkEquals(sort(names(hot.map)), expected.names)

} # test_createMotifClassMap
#----------------------------------------------------------------------------------------------------
test_getGCContent <- function(){

    message("---test_getGCContent")
        
    # Add a column of GC content
    gc.added <- merged.df %>%
        dplyr::mutate("gc_content" = getGCContent(start, endpos, chrom))

    # Check that there's 1 more column and the same number of rows
    checkTrue(ncol(gc.added) == ncol(merged.df) + 1)
    checkTrue(nrow(gc.added) == nrow(merged.df))

    # Check that it's numeric
    checkEquals(class(gc.added$gc_content), "numeric")
    
} # test_getGCContent
#----------------------------------------------------------------------------------------------------
test_addTSSDistance <- function(){

    message("---test_addTSSDistance")

    # Add a column of TSS distance, transformed by asinh
    tss.added <- addTSSDistance(merged.df)

    # Make sure we added a column, but not a row
    checkTrue(ncol(tss.added) == ncol(merged.df) + 1)
    checkTrue(nrow(tss.added) == nrow(merged.df))

    # Check that it's a numeric
    checkEquals(class(tss.added$asinh_tss_dist), "numeric")
    
} # test_addTSSDistance
#----------------------------------------------------------------------------------------------------
