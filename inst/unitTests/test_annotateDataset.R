library(FPML)
library(RUnit)
#----------------------------------------------------------------------------------------------------

runTests <- function(){

    test_annotateFootprintData()

    } # runTests

#----------------------------------------------------------------------------------------------------
# Test single TF function

test_annotateFootprintData <- function(){

    message("---test_annotateFootprintData")

    # Load the play merged data (merged.df)
    load(system.file(package="FPML", "extdata/annotateTestInput.Rdata"))

    # Get the results of annotating
    annotated.df <- annotateFootprintData(merged.df)

    # Check that the dimensions are correct


    # Check that column names are correct

    

} # test_createTfDf
#----------------------------------------------------------------------------------------------------
# Test TF-mapping function

test_mapTFsToMotifs <- function(){

    message("---test_mapTFsToMotifs")

    # Create a small set of TFs to map
    my.tfs <- c("E2F4", "BRCA1", "SRF", "NRF1", "YY1",
                "TBL1XR1", "IRF3", "GABPB1", "TAF1","BATF"   )
    fake.df <- data_frame(name = my.tfs)
    fake.map <- mapTFsToMotifs(fake.df)

    # The mapping doesn't find one of them
    checkEquals(length(fake.map), 9)
    checkEquals(setdiff(my.tfs, names(fake.map)), "GABPB1")

    # Check the data type
    all.types <- sapply(fake.map, class)
    checkTrue(all(all.types == "character"))

    # Check a couple elements
    checkEquals(length(fake.map[[1]]), 12)
    checkEquals(length(fake.map[[8]]), 5)
    

} # test_mapTFsToMotifs
#----------------------------------------------------------------------------------------------------
# Test Sampling Function
test_sampleTfDataset <- function(){

    message("---test_sampleTfDataset")

    # Make a fake dataframe and sample it
    fake.df <- data_frame(x = sample(1:1e6), y = sample(1:1e6))

    fake.subset <- sampleTfDataset(fake.df, 1000)

    checkEquals(nrow(fake.subset), 1000)
    checkEquals(ncol(fake.subset), 2)
    checkEquals(dplyr::intersect(fake.df, fake.subset), fake.subset)
    
} # test_sampleTfDataset
#----------------------------------------------------------------------------------------------------
