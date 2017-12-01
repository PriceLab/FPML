library(FPML)
library(RUnit)
#----------------------------------------------------------------------------------------------------

runTests <- function(){

    test_mergeFootprintsOneChrom()

    } # runTests

#----------------------------------------------------------------------------------------------------
# Test lymphoblast construction function

## Too intensive; it just calls the below function, so do that for now


#----------------------------------------------------------------------------------------------------
# Test single TF function

test_mergeFootprintsOneChrom <- function(){

    message("---test_mergeFootprintsOneChrom")    


    # Load the 5 data subsets
    load(system.file(package="FPML", "extdata/mergeTestData.1000.Rdata"))
    load(system.file(package="FPML", "extdata/fpSubTables.Rdata"))

    # Run the one-chromosome function
    results <- mergeFootprintsOneChrom(chrom_str = "Y",
                                       fimo_tbl = sample.data,
                                       hint_regions_tbl = hint.regions.Y,
                                       hint_hits_tbl = hint.hits.Y,
                                       well_regions_tbl= well.regions.Y,
                                       well_hits_tbl= well.hits.Y)

    # Should be 7 x 14
    checkTrue(ncol(results) == 16)
    checkTrue(nrow(results) == 7)

    # Check that the columns are all right
    expected.names <- sort(c("motifname", "chrom", "start", "endpos", "strand", "motifscore",
                             "pval", "sequence", "loc", "cs_hit", "h_count", "h_max_score",
                             "h_percent_overlap", "w_count", "w_min_score", "w_percent_overlap"))
    
    checkEquals(sort(names(results)), expected.names)
    
} # test_mergeFootprintsOneChrom
#----------------------------------------------------------------------------------------------------
