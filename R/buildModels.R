#----------------------------------------------------------------------------------------------------
#' Build a Linear Models
#'
#' Build a linear model with all features and linear models with each individual feature for the
#' given list of data and return the stats data framed for inspection and plotting
#'
#' @param dataList The list of prepared data, generally obtained by running the \code{prepModelData}
#' function. It must contain at least these 4  matrices: X_train, y_train, X_test, y_test
#'
#' @return A named list containing two items: 1) a  stats dataframe with statistics for the both the
#' full linear model and for the individual linear models; 2) the full linear model
#'
#' @export

buildLinearModels <- function(dataList){

    # Load the motif data
    load(system.file(package="FPML", "extdata/motif_class_pairs.Rdata"))
    
    # Make data ready for linear models
    X_train_lin <- dataList$X_train
    y_train_lin <- dataList$y_train
    X_test_lin  <- dataList$X_test
    y_test_lin  <- dataList$y_test
    
    # Make the Full linear model
    tf.regressors <- colnames(X_train_lin)[colnames(X_train_lin) %in% unique(motif.class$class)]
    non.tf.regressors <-  colnames(X_train_lin)[!colnames(X_train_lin) %in% unique(motif.class$class)]

    tf.regressors.formula <- paste("as.factor(",
                                   paste(tf.regressors, collapse=") + as.factor("), ")")
    non.tf.regressors.formula <- paste(non.tf.regressors, collapse=" + ")
    all.regressors.formula <- paste(non.tf.regressors.formula, tf.regressors.formula, sep=" + ")

    glm.formula <- paste("ChIPseq.bound ~ ", all.regressors.formula, sep='')

    glm.df.train <- as.data.frame(cbind(y_train_lin, X_train_lin)) %>%
        dplyr::rename("ChIPseq.bound" = "cs_hit")
    glm.df.test <-  as.data.frame(cbind(y_test_lin, X_test_lin)) %>%
        dplyr::rename("ChIPseq.bound" = "cs_hit")

    glm.all <- stats::glm(as.formula(glm.formula), data=glm.df.train, family=binomial)
    glm.all$Model.Name <- "linear model (all regressors)"

    # Gather full linear stats
    glm.pred.df <- make.pred.df.from.glm(glm.all, glm.df.test)
    glm.stat.df <- make.stats.df.from.preds(glm.pred.df)

    # Make individual linear models
    stats.regressors.df <- data.frame()
    
    for (this.regressor in colnames(X_train_lin)) {
        
        if (this.regressor %in% unique(motif.class$class)) {
            glm.formula <- paste("ChIPseq.bound ~ ", "as.factor(", this.regressor, ")",sep='')
        } else {
            glm.formula <- paste("ChIPseq.bound ~ ", this.regressor,sep='')
        }
        
        glm.df.train <- as.data.frame(cbind(y_train_lin, X_train_lin)) %>%
            dplyr::rename("ChIPseq.bound" = "cs_hit")
        glm.df.test <-  as.data.frame(cbind(y_test_lin, X_test_lin)) %>%
            dplyr::rename("ChIPseq.bound" = "cs_hit")

        glm.single <- stats::glm(as.formula(glm.formula),
                                  data=glm.df.train, family=binomial)
        glm.single$Model.Name <- paste("glm ", this.regressor, sep='')

        glm.pred.single.df <- make.pred.df.from.glm(glm.single, glm.df.test)
        glm.stat.single.df <- make.stats.df.from.preds(glm.pred.single.df)
    
        stats.regressors.df <- bind_rows(stats.regressors.df, glm.stat.single.df)        
    }

    linear.stat.df <- dplyr::bind_rows(glm.stat.df, stats.regressors.df)
    
    return(list(LinearStats = linear.stat.df,
                LinearModel = glm.all))
    
} # buildLinearModels
#----------------------------------------------------------------------------------------------------
#' Plot the Importance Matrix for a Gradient Boosted Model
#'
#' Plot the importance matrix, showing the features with the largest gain, for a given gradient
#' boosted model and the full path to a desired .PNG file containing the figure
#'
#' @param gbModel A gradient boosted model, generally obtained using the \code{buildBoostedModel}
#' function
#' @param X_train The matrix of training data used to created the supplied model. This is generally
#' created as an element of the list returned by the \code{prepModelData} function
#' @param filePath The path to the .PNG file that will be created by the function. This file will
#' house the figure of the importance matrix
#'
#' @return A .PNG file containing the importance matrix figure and located at the supplied file path
#'
#' @export

plotImportanceMatrix <- function(gbModel, X_train, filePath){

    # Load motifs as "motif.class"
    load(system.file(package="FPML", "extdata/motif_class_pairs.Rdata"))

    # Make the importance matrix figure
    motif.class$class <- lapply(motif.class$class, make.names, unique=TRUE)
    importance_matrix <- xgboost::xgb.importance(colnames(X_train),model=gbModel)
    
    df <- tibble::as_data_frame(importance_matrix)
    df.tf <- subset(df, Feature %in% unique(motif.class$class))
    df.notf <- subset(df, !(Feature %in% unique(motif.class$class)))
    tfclass.row <- c("TF_class", unname(as.list(colSums(df.tf[!(colnames(df.tf) %in% c("Feature"))]))) )
    names(tfclass.row) <- colnames(df)
    df.sum <- dplyr::bind_rows(df.notf,tfclass.row)


    g <- ggplot2::ggplot(data=df.sum, ggplot2::aes(x=reorder(Feature, Gain), y=Gain)) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::coord_flip() +
        ggplot2::theme_minimal(base_size = 30) +
        ggplot2::labs(x = "Feature", y="Gain")

    ggplot2::ggsave(filename = filePath,
                    plot = g)
    
} # plotImportanceMatrix
#----------------------------------------------------------------------------------------------------
#' Build a Gradient Boosted Model
#'
#' Build the gradient boosted model for the given list of data and return the stats data frame for
#' inspection and plotting
#'
#' @param dataList The list of prepared data, generally obtained by running the \code{prepModelData}
#' function. It must contain at least these 4  matrices: X_train, y_train, X_test, y_test
#'
#' @return The stats dataframe from the model and the model itself
#'
#' @export

buildBoostedModel <- function(dataList){

    # Pull out the relevant data
    X_train <- dataList$X_train
    y_train <- dataList$y_train
    X_test <- dataList$X_test
    y_test <- dataList$y_test
    
    # Train the boosted model
    param <- list("objective" = "binary:logistic",
                  "max.depth" = 7,
                  "eta" = 0.005,
                  "eval.metric" = "auc"
                  )

    gbdt_medium <- xgboost::xgboost(params = param,
                                    data = X_train,
                                    label = y_train,
                                    nround = 200,
                                    verbose = FALSE,
                                    missing = NA
                                    )
    gbdt_medium$Model.Name <- "gradient boosted model"

    # Gather stats from boosted model
    medium_pred_df <- make.pred.df.from.model(gbdt_medium, X_test, y_test)
    colnames(medium_pred_df)[1] <- "ChIPseq.bound"
    medium_stat_df <- make.stats.df.from.preds(medium_pred_df)


    return(list(BoostedStats = medium_stat_df,
                BoostedModel = gbdt_medium))
}
#----------------------------------------------------------------------------------------------------
#' Prepare Dataset for Modeling
#'
#' Prepare the fully annotated FIMO motif dataset for running models by filtering it accordingly
#'
#' @param annotated.df The fully annotated FIMO motif dataset. This must be a data frame, or else
#' the function will stop and give an error.
#' @param seed A character string that should be one of c("16","20","both"). This seed should
#' match the data in the data frame
#' @param hintCutoff A cutoff value used as a threshold on the HINT footprint score. Only scores
#' above this threshold will be kept (default = -Inf)
#' @param wellCutoff A cutoff value used as a threshold on the Wellington footprint score. Only
#' scores below this threshold will be kept (default = Inf)
#' @param motifsOnly A Boolean value indicating whether to prepare data such that it uses only motifs,
#' neglecting the footprint data (default = FALSE)
#'
#' @return A list in which the filtered data frame has been split into 6 pieces for X/Y and
#' training/testing/validation sets. These sets are completely ready to run gradient boosted
#' and linear models and should be passed as a list to those commands accordingly
#'
#' @export

prepModelData <- function(annotated.df, seed,
                          hintCutoff = -Inf,
                          wellCutoff = Inf,
                          motifsOnly = FALSE,
                          biasTraining = FALSE,
                          biasRatio = 9){

    # Catch non-data frames
    stopifnot("data.frame" %in% class(annotated.df))

    # Fix columns and filter
    colnames(annotated.df) <- make.names(colnames(annotated.df), unique=TRUE)

    if(seed == "both"){
        annotated.df %>%
            dplyr::filter(h_max_score_16 > hintCutoff |
                          h_max_score_20 > hintCutoff |
                          w_min_score_16 < wellCutoff |
                          w_min_score_20 < wellCutoff ) ->
            filtered.df
    } else {
        annotated.df %>%
            dplyr::filter(h_max_score > hintCutoff |
                          w_min_score < wellCutoff ) ->
            filtered.df
    }

    rm(annotated.df)

    # Denote columns to drop based on the motifsOnly flag
    if(motifsOnly){

        cols_to_drop <- c('motifname', 'chrom', 'strand', 'loc',
                          'h_frac_16', 'h_frac_20', 'h_max_score_16', 'h_max_score_20',
                          'w_frac_16', 'w_frac_20', 'w_min_score_16', 'w_min_score_20',
                          'h_frac', 'w_frac', 'h_max_score', 'w_min_score',
                          'h_percent_overlap', 'w_percent_overlap')
    } else{
        cols_to_drop <- c('motifname', 'chrom', 'start', 'endpos',
                          'strand', 'pval', 'sequence','loc')
    }

    # Intersect the columns to drop with the actual columns to avoid warnings
    cols_to_drop <- intersect(names(filtered.df), cols_to_drop)
        
    # Split up the data into training/testing/validation    
    filtered.df %>%
        dplyr::filter(chrom %in% c("2","4")) %>%
        dplyr::select(-dplyr::one_of(cols_to_drop)) ->
        val_df

    filtered.df %>%
        dplyr::filter(chrom %in% c("1","3","5")) %>%
        dplyr::select(-dplyr::one_of(cols_to_drop)) ->
        test_df

    # Split for training, allowing for training bias
    if(biasTraining){
        filtered.df %>%
            dplyr::filter(!(chrom %in% c("1","2","3","4","5"))) ->
            train.init

        # Bias all the motifs
        train.list <- lapply(unique(train.init$motifname),
                                    createBiasOneMotif,
                                    train.df = train.init,
                                    negPosRatio = biasRatio)
        train.list %>%
            dplyr::bind_rows() %>%
            dplyr::select(-dplyr::one_of(cols_to_drop)) ->
            train_df
        # Remove the intermediate bits
        rm(train.list); rm(train.init)
        
    } else {
        filtered.df %>%
            dplyr::filter(!(chrom %in% c("1","2","3","4","5"))) %>%
            dplyr::select(-dplyr::one_of(cols_to_drop)) ->
            train_df
        }

    remove(filtered.df)    

    # Split predictors and responses
    val_df %>%
        dplyr::select(-cs_hit) %>%
        as.matrix ->
        X_val

    val_df %>%
        dplyr::select(cs_hit) %>%
        as.matrix ->
        y_val

    test_df %>%
        dplyr::select(-cs_hit) %>%
        as.matrix ->
        X_test

    test_df %>%
        dplyr::select(cs_hit) %>%
        as.matrix ->
        y_test

    train_df %>%
        dplyr::select(-cs_hit) %>%
        as.matrix ->
        X_train

    train_df %>%
        dplyr::select(cs_hit) %>%
        as.matrix ->
        y_train

    remove(val_df, test_df, train_df)


    # Return the prepared data as a list
    return(list(X_train = X_train,
                y_train = y_train,
                X_test = X_test,
                y_test = y_test,
                X_val = X_val,
                y_val = y_val))

} # prepModelData
#----------------------------------------------------------------------------------------------------
#' Create a biased training set for one motif
#'
#' Create a biased training set for one motif by reducing the ratio of ChIPseq non-hits to hits to a
#' given ratio. This should only be run on training data, as running it on testing or validation data
#' is not good practice. 
#'
#' @param train.df The data frame of training data
#' @param motif The motif of interest; this MUST be part of the dataset, else you'll get an error.
#' @param negPosRatio The desired ratio of negatives to positives, must be expressed as an integer
#' (default = 9)
#'
#' @return A reduced training set for the supplied motif, where the ratio of negatives to positives
#' does not exceed that supplied to the function. It may fall short if there are motifs where there
#' are insufficient negative points to fulfill the desired ratio, as it does not discard any
#' positives
#'
#' @export

createBiasOneMotif <- function(train.df, motif, negPosRatio = 9){

    # Separate out the negatives from the positives
    all.hits <- train.df %>% dplyr::filter(motifname == motif, cs_hit == 1)
    all.non.hits <- train.df %>% dplyr::filter(motifname == motif, cs_hit == 0)

    # Check to see how many rows and if ratio is already sufficient, return the df as is
    if(nrow(all.non.hits) <= negPosRatio * nrow(all.hits)) return(train.df)

    # Otherwise, use sampling to get the right number of negatives
    sub.non.hits <- sampleTfDataset(all.non.hits, negPosRatio * nrow(all.hits))

    # Combine and return them
    return(dplyr::bind_rows(all.hits, sub.non.hits))
    
} # createBiasOneMotif
#----------------------------------------------------------------------------------------------------
#' Join the Seed 16 and Seed 20 Data
#'
#' Join the seed 16 and seed 20 data in preparation for building models on the full dataset. This
#' joined data can be used to run models that either include or exclude footprint data.
#'
#' @param full.16 The annotated data frame containing the seed 16 data
#' @param full.20 The annotated data frame containing the seed 20 data
#'
#' @return A joined data frame, containing footprints from both seeds. The footprint columns of this
#' data frame have added suffixes to denote which seed they come from, as there are now 2 of each
#' column, one originating from each seed size. For example, "h_frac" becomes "h_frac_20" and
#' "h_frac_16"
#'
#' @export

joinModelData <- function(full.16, full.20){

    # Remove designated columns
    sub.20 <- dplyr::select(full.20, -pval, -sequence, -start, -endpos)
    sub.16 <- dplyr::select(full.16, -pval, -sequence, -start, -endpos)

    # Join them together
    my.keys <- setdiff(names(sub.20),
                       c("h_max_score","w_min_score","h_frac","w_frac",
                         "h_percent_overlap", "w_percent_overlap"))
    joined.df <- dplyr::full_join(sub.20,
                                  sub.16,
                                  by = my.keys,
                                  suffix = c("_20","_16"))
    return(joined.df)
    
} # joinModelData
#----------------------------------------------------------------------------------------------------    
#' Make plots for MCC, ROC, and Prec-Rec
#'
#' Given the stats dataframes for boosted and linear models, plus a flag denoting the data type,
#' plot the 3 relevant curves: the Matthews Correlation Coefficient curve, the ROC curve, and the
#' Precision-Recall curve. Also requires specification of file paths to store the .PNG files for
#' these curves
#'
#' @param boostedStats The dataframe containing statistics from the boosted model. Generally, these
#' are generated using the \code{buildBoostedModel} function
#' @param linearStatList The list containing 2 dataframes, one with the statistics from the full
#' linear model, the other with statistics from the individual linear models. Generally, these are
#' generated using the \code{buildLinearModels} function
#' @param dataType A string that must be one of c("seed16","seed20","both","motifsOnly"). This
#' string denotes what labels will be given for the plots within the figures themselves
#' @param filePaths A character string of length = 3, giving the 3 filepaths to which the function
#' will save the curves. They should be in the following format: c(MCC_path, ROC_path, PrecRec_path)
#'
#' @return 3 .PNG files, each containing one of the 3 figures for the given data.
#'
#' @export

plotStatCurves <- function(boostedStats, linearStats, dataType, filePaths){

    # Catch bad datatypes
    stopifnot(dataType %in% c("seed16", "seed20", "both", "motifsOnly"))

    # Bind together the input data
    all.stats.df <- bind_rows(boostedStats, linearStats)

    # Define the general color palette
    myPalette <- c("black","red","blue","green2","darkorange",
                   "purple","magenta","brown","cyan")

    # Based on dataType, change the following:
    # 1) List of models to filter out
    # 2) Color palette
    # 3) Individual linear model names
    # 4) Order of models in the legend
    if(dataType == "both"){       
        modelList <- c("glm h_max_score_16",
                       "glm w_min_score_16",
                       "glm h_max_score_20",
                       "glm w_min_score_20",
                       "gradient boosted model",
                       "linear model (all regressors)")

        # Filter the data frame accordingly
        all.stats.df <- all.stats.df %>%
            dplyr::filter(Model.Name %in% modelList)

        myColors <- myPalette[1:6]
        names(myColors) <- c('gradient boosted model',
                             'linear model (all regressors)',
                             'HINT best score seed 16',
                             'HINT best score seed 20',
                             'Wellington best score seed 16',
                             'Wellington best score seed 20')
        all.stats.df$Model.Name[all.stats.df$Model.Name== 'glm h_max_score_16'] <-
            'HINT best score seed 16'
        all.stats.df$Model.Name[all.stats.df$Model.Name== 'glm h_max_score_20'] <-
            'HINT best score seed 20'
        all.stats.df$Model.Name[all.stats.df$Model.Name== 'glm w_min_score_16'] <-
            'Wellington best score seed 16'
        all.stats.df$Model.Name[all.stats.df$Model.Name== 'glm w_min_score_20'] <-
            'Wellington best score seed 20'

        
        all.stats.df$Model.Name <- factor(all.stats.df$Model.Name,
                                          levels = c('gradient boosted model',
                                                     'linear model (all regressors)',
                                                     'HINT best score seed 16',
                                                     'HINT best score seed 20',
                                                     'Wellington best score seed 16',
                                                     'Wellington best score seed 20'
                                                     )
                                          )
    } else if(dataType == "motifsOnly"){
        modelList <- c("glm motifscore",
                       "glm gc_content",
                       "glm asinh_ts_dist",
                       "gradient boosted model",
                       "linear model (all regressors)")

        # Filter the data frame accordingly
        all.stats.df <- all.stats.df %>%
            dplyr::filter(Model.Name %in% modelList)
        
        myColors <- myPalette[c(1,2,7,8,9)]
        names(myColors) <- c('gradient boosted model',
                             'linear model (all regressors)',
                             'GC content',
                             'motif score',
                             'arcsinh(TSS distance)'
                             )
        all.stats.df$Model.Name[all.stats.df$Model.Name== 'glm gc_content'] <-
            'GC content'
        all.stats.df$Model.Name[all.stats.df$Model.Name== 'glm motifscore'] <-
            'motif score'
        all.stats.df$Model.Name[all.stats.df$Model.Name== 'glm asinh_tss_dist'] <-
            'arcsinh(TSS distance)'
        all.stats.df$Model.Name <- factor(all.stats.df$Model.Name,
                                          levels = c('gradient boosted model',
                                                     'linear model (all regressors)',
                                                     'GC content',
                                                     'motif score',
                                                     'arcsinh(TSS distance)'
                                                     )
                                          )
        
    } else if(dataType == "seed16"){
        modelList <- c("gradient boosted model",
                       "linear model (all regressors)",
                       "glm h_max_score",
                       "glm w_min_score")

        # Filter the data frame accordingly
        all.stats.df <- all.stats.df %>%
            dplyr::filter(Model.Name %in% modelList)

        myColors <- myPalette[1:4]
        names(myColors) <- c('gradient boosted model',
                             'linear model (all regressors)',
                             'HINT best score seed 16',
                             'Wellington best score seed 16'
                             )
        all.stats.df$Model.Name[all.stats.df$Model.Name== 'glm h_max_score'] <-
            'HINT best score seed 16'
        all.stats.df$Model.Name[all.stats.df$Model.Name== 'glm w_min_score'] <-
            'Wellington best score seed 16'
        all.stats.df$Model.Name <- factor(all.stats.df$Model.Name,
                                          levels = c("gradient boosted model",
                                                     "linear model (all regressors)",
                                                     "HINT best score seed 16",
                                                     "Wellington best score seed 16"
                                                     )
                                          )
    } else {
        modelList <- c("gradient boosted model",
                       "linear model (all regressors)",
                       "glm h_max_score",
                       "glm w_min_score")

        # Filter the data frame accordingly
        all.stats.df <- all.stats.df %>%
            dplyr::filter(Model.Name %in% modelList)

        myColors <- myPalette[1:4]
        names(myColors) <- c('gradient boosted model',
                             'linear model (all regressors)',
                             'HINT best score seed 20',
                             'Wellington best score seed 20'
                             )
        all.stats.df$Model.Name[all.stats.df$Model.Name== 'glm h_max_score'] <-
            'HINT best score seed 20'
        all.stats.df$Model.Name[all.stats.df$Model.Name== 'glm w_min_score'] <-
            'Wellington best score seed 20'
        all.stats.df$Model.Name <- factor(all.stats.df$Model.Name,
                                          levels = c("gradient boosted model",
                                                     "linear model (all regressors)",
                                                     "HINT best score seed 20",
                                                     "Wellington best score seed 20"
                                                     )
                                          )
    } 
    

    # Make the color palette
    colScale <- ggplot2::scale_colour_manual(name = "Model.Name",values = myColors)

#    browser()
    # Plot the Matt CC Curve
    p1 <- plot.mattcc.curve(all.stats.df) +
        ggplot2::theme_minimal(base_size = 15) + colScale
    ggplot2::ggsave(filename = filePaths[1],
                    plot = p1)

    # Plot the ROC curve
    p2 <- plot.roc.curve(all.stats.df, TRUE) +
        ggplot2::theme_minimal(base_size = 15) + colScale
    ggplot2::ggsave(filename = filePaths[2],
                    plot = p2)
    
    # Plot the Prec-Rec curve
    p3 <- plot.precrecall.curve(all.stats.df) +
        ggplot2::theme_minimal(base_size = 15) + colScale
    ggplot2::ggsave(filename = filePaths[3],
                    plot = p3)
}

#----------------------------------------------------------------------------------------------------
#' Extract the Maximum MCC Values
#'
#' Using the statistics data frames from building models, extract the maximum MCC values for
#' each model
#'
#' @param statsDF A data frame of statistics created by one of the model building functions,
#' \code{buildBoostedModel} or \code{buildLinearModels}. This can only be a data frame, not a
#' list of data frames. In order to use multiple data frames, bind them together into one data
#' frame first
#'
#' @return A table summarizing the maximum MCC for each model in the supplied data frame
#'
#' @export

extractMaxMCC <- function(boostedStats, linearStats){

    # Check to make sure both are dataframes
    stopifnot("data.frame" %in% class(boostedStats))
    stopifnot("data.frame" %in% class(linearStats))

    # Bind them together
    statsDF <- bind_rows(boostedStats, linearStats)


    # Use a dplyr summarization to get the correct table
    statsDF %>% dplyr::group_by(Model.Name) %>%
        summarise(Max.MCC = max(MattCC)) %>%
        dplyr::arrange(desc(Max.MCC)) ->
            mccDF

    #Turn it into a vector for each calling
    mccs <- mccDF$Max.MCC
    names(mccs) <- mccDF$Model.Name

    return(mccs)
}
#----------------------------------------------------------------------------------------------------
#' Create the "truth plot"
#'
#' Create a "truth" plot, which plots HINT v. Wellington scores given a model and a threshold
#' to use for determining a "positive"
#'
#' @param model A model trained on the annotated dataset, either a linear or gradient-boosted model
#' @param modelType Either "boosted" or "linear" to indicate the type of model provided
#' @param seed One of "16", "20", or "both", denoting which seed is being used
#' @param threshold A numeric from 0-1 indicating the cutoff for what constitutes a "positive"
#' probability. This is an INCLUSIVE value; positives are greater than or equal to the threshold
#' @param X_test The matrix of test data for predictors, generally created using "prepModelData"
#' @param y_test The matrix of test data for response, generally created using "prepModelData"
#' @param plotFile A location to use for creating the png plot file
#'
#' @return The "truth data frame", containing the Wellington/HINT scores for each point in the
#' test set, the predicted probability of a ChIPseq hit, the binary prediction based on that
#' probability and the supplied threshold, the true ChIPseq hit value, and the corresponding
#' label of TP/TN/FP/FN.
#'
#' It also constructs and saves a plot of HINT vs Wellington scores, where both have been
#' transformed via asinh() and the Wellington scores have been converted via absolute value.
#' Points are coded by confusion matrix outcome. 
#'
#' @export

createTruthPlot <- function(model, modelType, seed, threshold,
                            X_test, y_test, plotFile){

    # Combine the test data and make the predictions with the correct function
    if(modelType == "boosted"){

        pred.df <- make.pred.df.from.model(model, X_test, y_test)
    } else {

        glm.df.test <-  as.data.frame(cbind(y_test, X_test)) %>%
            dplyr::rename("ChIPseq.bound" = "cs_hit")

        pred.df <- make.pred.df.from.glm(model, glm.df.test)
    }

    # Use seed to determine columns to look for
    if(seed == "both"){
        hint.col <- "h_max_score_20"
        well.col <- "w_min_score_20"
    } else {
        hint.col <- "h_max_score"
        well.col <- "w_min_score"
    }

    # Fix the column name for the prediction dataframe
    names(pred.df)[[1]] <- "ChIPseq.bound"

    # Add the predictions to X footprint data to make a dataframe
    truth.df <- dplyr::data_frame(abs_w_min_score = asinh(abs(X_test[,well.col])),
                                  h_max_score = asinh(X_test[,hint.col]),
                                  true_cs_hit = pred.df$ChIPseq.bound,
                                  prediction = pred.df$Prediction)

    # Add the predictions for the threshold
    # If its greater than or equal to the threshold, make it positive
    truth.df <- truth.df %>%
        dplyr::mutate(pred_cs_hit = ifelse(prediction >= threshold, 1, 0))

    # Create a labeling function
    addLabel <- function(prediction, truth){

        if(prediction){
            if(truth){ label <- "TP"
            } else { label <- "FP"}
        } else {
            if(truth){ label <- "FN"
            } else { label <- "TN"}
        }
        return(label)
    }

    # Add label to df
    truth.df <- truth.df %>%
        dplyr::rowwise() %>%
        dplyr::mutate( Label = addLabel(pred_cs_hit, true_cs_hit))

    # Create a plot and save it to the file
    truth.plot <- ggplot2::ggplot(truth.df) +
        ggplot2::geom_point(ggplot2::aes(x = abs_w_min_score,
                                         y = h_max_score,
                                         color = Label)) +
        ggplot2::theme_minimal(base_size = 15)

    # Save it
    ggplot2::ggsave(filename = plotFile,
                    plot = truth.plot)

    return(truth.df)
    
} # createTruthPlot
#----------------------------------------------------------------------------------------------------
#' Create a threshold plot
#'
#' Create a "threshold" plot, which plots the 3 relevant prediction metrics against threshold to
#' enable identification of the ideal "threshold" for ChIPseq hit probability.
#'
#' @param stats.df A dataframe of statistics generated by one of the model runs
#' @param modelName A string containing a model name from the stats data frame. To see the options
#' available, use unique(stats.df$Model.Name)
#' @param plotFile A path to a file for saving the generated plot as a .png file
#'
#' @return A .png file containing a plot of 3 metrics (sensitivity, specificity, ppv) as a function
#' of threshold
#'
#' @export

createThresholdPlot <- function(stats.df, modelName, plotFile){

    # These parts all failed because of weirdness with ggplot2
    # Melt the dataframe
    #thresh.df <- stats.df %>% 
    #  dplyr::filter(Model.Name == modelName) %>% 
    #  reshape2::melt(id.vars = "threshold",
    #                 variable.name = "Metric",
    #                 value.name = "Value") %>%
    #  dplyr::filter(Metric %in% c("sensitivity", "specificity", "ppv", "npv"))
    
    # Format the Value column to be non-scientific
    #thresh.df$Value <- format(thresh.df$Value, scientific = FALSE)
    # Create the plot
    #thresh.plot <- ggplot2::ggplot(thresh.df) +
    #    ggplot2::geom_point(ggplot2::aes(x = threshold,
    #                                    y = Value,
    #                                    color = Metric)) +
    #    ggplot2::theme_minimal(base_size = 15)

    # Save it
    #ggplot2::ggsave(filename = plotFile,
    #                plot = thresh.plot)

    # Do it with Base plot
    filtered.df <- stats.df %>% dplyr::filter(Model.Name == modelName)
    
    # Plot the 4 metrics (not just 3)
    png(filename = plotFile)
    with(filtered.df, plot(threshold, sensitivity, col = "blue", pch = 20,
                           xlab = "Threshold", ylab = "Metric"))
    with(filtered.df, points(threshold, specificity, col = "red", pch = 20))
    with(filtered.df, points(threshold, ppv, col = "green", pch = 20))
    with(filtered.df, points(threshold, npv, col = "orange", pch = 20))
    
    legend("topright", legend = c("Sens. (Recall)", "Specificity",
                                  "PPV (Precision)", "NPV"),
           pch = c(20, 20, 20, 20),
           col = c("blue", "red", "green", "orange"))
    dev.off()
  
} # createThresholdPlot
#----------------------------------------------------------------------------------------------------


