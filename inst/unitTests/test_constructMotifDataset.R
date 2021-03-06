library(FPML)
library(RUnit)
#----------------------------------------------------------------------------------------------------
runTests <- function(){

    test_constructLymphoblastDataset # Takes a while, ~5 min
    test_createTfDf() # Takes a while, ~5 min
    test_mapTFsToMotifs()
    test_sampleTfDataset()

    } # runTests

#----------------------------------------------------------------------------------------------------
# Test lymphoblast construction function
test_constructLymphoblastDataset <- function(){

    # Run the function using the "isTest" flag
    result <- constructLymphoblastDataset(isTest = TRUE, distinctFlag = FALSE)

    # Check rows in the result
    checkTrue(is.data.frame(result))
    checkEquals(nrow(result), 12604861)
    
    # Check unique motifs (there should be only 2)
    checkEquals(length(unique(result$motifname)), 2)

    # Check columns names
    column.names <- c("motifname", "chrom", "start", "endpos", "strand",
                      "motifscore", "pval", "empty", "sequence", "loc", "cs_hit")
    checkTrue(all(column.names %in% names(result)))

} # test_constructLymphoblstDataset
#----------------------------------------------------------------------------------------------------
# Test single TF function
test_createTfDf <- function(){

    message("---test_createTfDf")

    # Bring in the chipseq data
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

    # Run on a TF with only 1 motif
    my.TF <- "BCL3"
    result <- createTfDf(my.TF, TFs.to.motifs, chipseq.regions, chipseq.hits, TRUE)

    # Run tests
    checkTrue(is.data.frame(result))
    checkEquals(nrow(result), 2843902)

    # Split Positive/negative
    pos <- result %>% filter(cs_hit == 1)
    neg <- result %>% filter(cs_hit == 0)

    checkEquals(nrow(pos), 33752)
    checkEquals(nrow(neg), 2810150)

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
