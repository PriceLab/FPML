library(FPML)
library(RUnit)
#----------------------------------------------------------------------------------------------------
# Load the data
load(system.file(package = "FPML", "extdata/buildTestInput.Rdata")


runTests <- function(){

    test_mapTFsToMotifs()
    test_sampleTfDataset()

    } # runTests

#----------------------------------------------------------------------------------------------------
# Test data prep

test_prepModelData <- function(){

    message("---test_prepModelData")

    # Simply use the random seed 20 data and prep with defaults
    prepped.data <- prepModelData(random.20, "20")
    

    # Check for 6 list elements with correct names
    checkEquals(length(prepped.data), 6)
    expected.names <- sort(c("X_train", "y_train", "X_test", "y_test", "X_val", "y_val"))
    checkEquals(sort(names(prepped.data)), expected.names)

    # Check dimensions of training predictors
    checkEquals(nrow(prepped.data$X_train), 637)
    checkEquals(ncol(prepped.data$X_train), 30)

    # Now do it with bias


    # Now do it with cutoff for HINT


    # Now do it with cutoff for Wellington


    # Now do it with only motifs
}




#----------------------------------------------------------------------------------------------------
# Test single TF function

test_joinData <- function(){

    message("---test_createTfDf")

    # Join the data frames together
    joined.df <- joinModelData(sub.16, sub.20)

    # Check the dimensions
    checkEquals(nrow(joined.df), 1000)
    checkEquals(ncol(joined.df), 39)    

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
